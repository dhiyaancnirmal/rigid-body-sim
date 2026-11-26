// C++ tests for contracts that live below the Python layer: mass properties,
// bit-exact replay, broadphase equivalence and the cross-machine golden hash.
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "rbsim/scenes.hpp"

using namespace rbsim;

namespace {

int failures = 0;

#define CHECK(cond)                                                              \
    do {                                                                         \
        if (!(cond)) {                                                           \
            std::printf("  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);        \
            ++failures;                                                          \
        }                                                                        \
    } while (0)

bool near(double a, double b, double tol) { return std::fabs(a - b) <= tol; }

void polygon_mass_matches_closed_form() {
    // Solid rectangle: I = m (w^2 + h^2) / 12.
    const Shape box = make_box(1.5, 0.25);
    const MassData md = compute_mass(box, 2.0);
    const double m = 2.0 * 3.0 * 0.5;
    CHECK(near(md.mass, m, 1e-12));
    CHECK(near(md.inertia, m * (9.0 + 0.25) / 12.0, 1e-12));

    // Right triangle with legs 3 and 3, given off-centre: centroid at (1, 1)
    // relative to the right angle; I about centroid = m (a^2 + b^2) / 18.
    Vec2 centroid;
    const Shape tri = make_polygon({{10.0, 0.0}, {13.0, 0.0}, {10.0, 3.0}}, &centroid);
    CHECK(near(centroid.x, 11.0, 1e-12) && near(centroid.y, 1.0, 1e-12));
    const MassData mt = compute_mass(tri, 1.0);
    CHECK(near(mt.mass, 4.5, 1e-12));
    CHECK(near(mt.inertia, 4.5 * 18.0 / 18.0, 1e-12));
}

void replay_is_bit_exact() {
    World a = scenes::pile(150);
    World b = scenes::pile(150);
    bool same = true;
    for (int i = 0; i < 600; ++i) {
        a.step(1.0 / 60.0);
        b.step(1.0 / 60.0);
        same = same && a.state_hash() == b.state_hash();
    }
    CHECK(same);
}

void broadphases_give_identical_trajectories() {
    World brute = scenes::pile(150, Broadphase::BruteForce);
    World sap = scenes::pile(150, Broadphase::SweepAndPrune);
    bool same = true;
    for (int i = 0; i < 600; ++i) {
        brute.step(1.0 / 60.0);
        sap.step(1.0 / 60.0);
        same = same && brute.state_hash() == sap.state_hash() &&
               brute.pair_count() == sap.pair_count();
    }
    CHECK(same);
}

void golden_hash_matches() {
    World w = scenes::pile(100);
    for (int i = 0; i < 600; ++i) w.step(1.0 / 60.0);
    const std::uint64_t h = w.state_hash();
    std::printf("  pile(100) hash after 600 steps: 0x%016" PRIx64 "\n", h);
    CHECK(h == scenes::kPile100Hash600);
}

void stack_of_ten_settles() {
    World w = scenes::stack(10);
    for (int i = 0; i < 600; ++i) w.step(1.0 / 60.0);
    double max_dx = 0.0;
    for (std::size_t i = 1; i < w.bodies.size(); ++i) {
        max_dx = std::fmax(max_dx, std::fabs(w.bodies[i].x.x));
    }
    const Body& top = w.bodies.back();
    std::printf("  stack(10): max sideways drift %.2e m, top y %.4f m, KE %.2e J\n", max_dx,
                top.x.y, w.kinetic_energy());
    CHECK(max_dx < 1e-3);
    CHECK(near(top.x.y, 9.5, 0.05));  // sinks at most linear_slop per contact
    CHECK(w.kinetic_energy() < 1e-3);
}

}  // namespace

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"polygon_mass_matches_closed_form", polygon_mass_matches_closed_form},
        {"replay_is_bit_exact", replay_is_bit_exact},
        {"broadphases_give_identical_trajectories", broadphases_give_identical_trajectories},
        {"golden_hash_matches", golden_hash_matches},
        {"stack_of_ten_settles", stack_of_ten_settles},
    };
    for (const auto& [name, fn] : tests) {
        const int before = failures;
        fn();
        std::printf("%s %s\n", failures == before ? "ok  " : "FAIL", name.c_str());
    }
    std::printf("%d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
