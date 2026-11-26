// Integrator study: explicit Euler, semi-implicit (symplectic) Euler,
// velocity Verlet and classical RK4 on three systems with known periods:
// a harmonic oscillator, a large-amplitude pendulum and an eccentric Kepler
// orbit. Each run covers exactly 100 periods, so the exact final state equals
// the initial state. Prints JSON.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "rbsim/math.hpp"

using rbsim::Vec2;

namespace {

constexpr double kPi = 3.14159265358979323846;

struct System {
    std::string name;
    Vec2 q0, v0;
    double period;
    Vec2 (*accel)(Vec2);
    double (*energy)(Vec2, Vec2);
};

Vec2 sho_accel(Vec2 q) { return -1.0 * q; }
double sho_energy(Vec2 q, Vec2 v) { return 0.5 * (dot(v, v) + dot(q, q)); }

Vec2 pendulum_accel(Vec2 q) { return {-std::sin(q.x), 0.0}; }
double pendulum_energy(Vec2 q, Vec2 v) { return 0.5 * v.x * v.x - std::cos(q.x); }

Vec2 kepler_accel(Vec2 q) {
    const double r = rbsim::length(q);
    return (-1.0 / (r * r * r)) * q;
}
double kepler_energy(Vec2 q, Vec2 v) { return 0.5 * dot(v, v) - 1.0 / rbsim::length(q); }

double pendulum_period(double amplitude) {
    // T = 4 K(sin(a/2)) for w0 = 1, with K from the arithmetic-geometric mean.
    const double k = std::sin(0.5 * amplitude);
    double a = 1.0, b = std::sqrt(1.0 - k * k);
    for (int i = 0; i < 60 && std::fabs(a - b) > 1e-16; ++i) {
        const double an = 0.5 * (a + b);
        b = std::sqrt(a * b);
        a = an;
    }
    return 4.0 * (kPi / (2.0 * a));
}

struct State {
    Vec2 q, v;
};

using Stepper = State (*)(const System&, State, double);

State explicit_euler(const System& s, State x, double h) {
    const Vec2 a = s.accel(x.q);
    return {x.q + h * x.v, x.v + h * a};
}

State semi_implicit_euler(const System& s, State x, double h) {
    const Vec2 v = x.v + h * s.accel(x.q);
    return {x.q + h * v, v};
}

// Kick-drift-kick form. The trailing acceleration is the next step's leading
// one, so a real implementation evaluates the force once per step; counting
// and timing below treat it that way by caching.
struct VerletCache {
    bool valid = false;
    Vec2 q, a;
};
VerletCache g_verlet;

State velocity_verlet(const System& s, State x, double h) {
    Vec2 a0;
    if (g_verlet.valid && g_verlet.q.x == x.q.x && g_verlet.q.y == x.q.y) {
        a0 = g_verlet.a;
    } else {
        a0 = s.accel(x.q);
    }
    const Vec2 vh = x.v + (0.5 * h) * a0;
    const Vec2 q1 = x.q + h * vh;
    const Vec2 a1 = s.accel(q1);
    g_verlet = {true, q1, a1};
    return {q1, vh + (0.5 * h) * a1};
}

State rk4(const System& s, State x, double h) {
    const Vec2 k1q = x.v, k1v = s.accel(x.q);
    const Vec2 k2q = x.v + (0.5 * h) * k1v, k2v = s.accel(x.q + (0.5 * h) * k1q);
    const Vec2 k3q = x.v + (0.5 * h) * k2v, k3v = s.accel(x.q + (0.5 * h) * k2q);
    const Vec2 k4q = x.v + h * k3v, k4v = s.accel(x.q + h * k3q);
    return {x.q + (h / 6.0) * (k1q + 2.0 * k2q + 2.0 * k3q + k4q),
            x.v + (h / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v)};
}

struct Method {
    std::string name;
    Stepper step;
    int evals_per_step;
};

double g_sink = 0.0;

}  // namespace

int main() {
    const double pend_amp = 2.0;  // radians, about 115 degrees
    const double ecc = 0.5;
    const std::vector<System> systems = {
        {"harmonic_oscillator", {1.0, 0.0}, {0.0, 0.0}, 2.0 * kPi, sho_accel, sho_energy},
        {"pendulum_115deg", {pend_amp, 0.0}, {0.0, 0.0}, pendulum_period(pend_amp),
         pendulum_accel, pendulum_energy},
        {"kepler_e0.5", {1.0 - ecc, 0.0}, {0.0, std::sqrt((1.0 + ecc) / (1.0 - ecc))}, 2.0 * kPi,
         kepler_accel, kepler_energy},
    };
    const std::vector<Method> methods = {
        {"explicit_euler", explicit_euler, 1},
        {"semi_implicit_euler", semi_implicit_euler, 1},
        {"velocity_verlet", velocity_verlet, 1},
        {"rk4", rk4, 4},
    };
    // Force evaluations per period; RK4 takes a quarter as many steps, so
    // every method spends the same work at a given budget.
    const std::vector<int> evals_per_period = {200, 800, 3200};
    const int periods = 100;

    std::printf("{\n  \"periods\": %d,\n  \"runs\": [\n", periods);
    bool first = true;
    for (const System& sys : systems) {
        const double e0 = sys.energy(sys.q0, sys.v0);
        for (const Method& m : methods) {
            for (int budget : evals_per_period) {
                const int steps_per_period = budget / m.evals_per_step;
                const double h = sys.period / steps_per_period;
                const long total = static_cast<long>(steps_per_period) * periods;
                g_verlet.valid = false;
                State x{sys.q0, sys.v0};
                double max_de = 0.0;
                bool blew_up = false;
                for (long i = 0; i < total; ++i) {
                    x = m.step(sys, x, h);
                    const double de = std::fabs(sys.energy(x.q, x.v) - e0) / std::fabs(e0);
                    if (!(de < 1e6)) {
                        blew_up = true;
                        break;
                    }
                    if (de > max_de) max_de = de;
                }
                g_sink += x.q.x;
                const double pos_err = blew_up ? -1.0 : rbsim::length(x.q - sys.q0);
                // Cost: time the bare stepping loop (no energy check), best of 5.
                double ns = 1e300;
                for (int rep = 0; rep < 5; ++rep) {
                    g_verlet.valid = false;
                    State y{sys.q0, sys.v0};
                    const long n = total < 200000 ? total : 200000;
                    const auto t0 = std::chrono::steady_clock::now();
                    for (long i = 0; i < n; ++i) y = m.step(sys, y, h);
                    const auto t1 = std::chrono::steady_clock::now();
                    g_sink += y.q.x;
                    const double t = std::chrono::duration<double, std::nano>(t1 - t0).count() / n;
                    if (t < ns) ns = t;
                }
                std::printf("%s    {\"system\": \"%s\", \"method\": \"%s\", \"evals_per_period\": %d, "
                            "\"steps_per_period\": %d, \"dt\": %.6g, \"max_rel_energy_error\": %.3e, "
                            "\"final_position_error\": %.3e, \"blew_up\": %s, "
                            "\"ns_per_step\": %.2f}",
                            first ? "" : ",\n", sys.name.c_str(), m.name.c_str(), budget,
                            steps_per_period, h, blew_up ? -1.0 : max_de, pos_err,
                            blew_up ? "true" : "false", ns);
                first = false;
            }
        }
    }
    std::printf("\n  ],\n");

    // Long horizon: harmonic oscillator at 200 force evaluations per period
    // for 10,000 periods. Symplectic methods keep a bounded energy error;
    // RK4's error grows steadily.
    std::printf("  \"long_horizon\": {\"system\": \"harmonic_oscillator\", \"evals_per_period\": 200, "
                "\"max_rel_energy_error_by_periods\": [\n");
    first = true;
    for (const Method& m : methods) {
        if (m.name == "explicit_euler") continue;  // diverges within 100 periods
        const System& sys = systems[0];
        const double e0 = sys.energy(sys.q0, sys.v0);
        const int spp = 200 / m.evals_per_step;
        const double h = sys.period / spp;
        g_verlet.valid = false;
        State x{sys.q0, sys.v0};
        double max_de = 0.0;
        double at[3] = {0.0, 0.0, 0.0};
        const long marks[3] = {100, 1000, 10000};
        int mi = 0;
        for (long i = 1; i <= static_cast<long>(spp) * 10000; ++i) {
            x = m.step(sys, x, h);
            const double de = std::fabs(sys.energy(x.q, x.v) - e0) / e0;
            if (de > max_de) max_de = de;
            if (i == marks[mi] * spp) at[mi++] = max_de;
        }
        std::printf("%s    {\"method\": \"%s\", \"p100\": %.3e, \"p1000\": %.3e, \"p10000\": %.3e}",
                    first ? "" : ",\n", m.name.c_str(), at[0], at[1], at[2]);
        first = false;
    }
    std::printf("\n  ]},\n  \"sink\": %.3f\n}\n", g_sink);
    return 0;
}
