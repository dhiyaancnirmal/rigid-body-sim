#include "rbsim/world.hpp"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace rbsim {
namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point t0, Clock::time_point t1) {
    return std::chrono::duration<double>(t1 - t0).count();
}

std::uint64_t pair_key(int a, int b) {
    const auto lo = static_cast<std::uint64_t>(a < b ? a : b);
    const auto hi = static_cast<std::uint64_t>(a < b ? b : a);
    return (lo << 32) | hi;
}

bool overlaps(const AABB& a, const AABB& b) {
    return a.lower.x <= b.upper.x && b.lower.x <= a.upper.x && a.lower.y <= b.upper.y &&
           b.lower.y <= a.upper.y;
}

template <class B>
Vec2 velocity_at(const B& body, Vec2 r) { return body.v + cross(body.w, r); }

template <class B>
void apply_impulse(B& a, B& b, Vec2 ra, Vec2 rb, Vec2 p) {
    a.v -= a.inv_mass * p;
    a.w -= a.inv_inertia * cross(ra, p);
    b.v += b.inv_mass * p;
    b.w += b.inv_inertia * cross(rb, p);
}

// Solve both normal impulses of a two-point manifold at once, as in Box2D's
// b2ContactSolver: find x >= 0 with w = K x + b >= 0 and x_i w_i = 0 by
// trying the four complementarity cases in turn.
void solve_block(Manifold& m, SolverBody& a, SolverBody& b) {
    ContactPoint& c1 = m.points[0];
    ContactPoint& c2 = m.points[1];
    const Vec2 n = m.normal;
    const double a1 = c1.normal_impulse, a2 = c2.normal_impulse;
    const double vn1 = dot(velocity_at(b, c1.rb) - velocity_at(a, c1.ra), n);
    const double vn2 = dot(velocity_at(b, c2.rb) - velocity_at(a, c2.ra), n);
    // b' = vn - target - K * a  (velocities as if the accumulated impulse were removed)
    const double b1 = vn1 - c1.target_velocity - (m.k11 * a1 + m.k12 * a2);
    const double b2 = vn2 - c2.target_velocity - (m.k12 * a1 + m.k22 * a2);

    double x1 = 0.0, x2 = 0.0;
    bool found = false;
    // Case 1: both points push.
    const double det = m.k11 * m.k22 - m.k12 * m.k12;
    const double inv_det = 1.0 / det;
    x1 = -inv_det * (m.k22 * b1 - m.k12 * b2);
    x2 = -inv_det * (m.k11 * b2 - m.k12 * b1);
    if (x1 >= 0.0 && x2 >= 0.0) found = true;
    // Case 2: only point 1 pushes.
    if (!found) {
        x1 = -b1 / m.k11;
        x2 = 0.0;
        if (x1 >= 0.0 && m.k12 * x1 + b2 >= 0.0) found = true;
    }
    // Case 3: only point 2 pushes.
    if (!found) {
        x1 = 0.0;
        x2 = -b2 / m.k22;
        if (x2 >= 0.0 && m.k12 * x2 + b1 >= 0.0) found = true;
    }
    // Case 4: neither pushes.
    if (!found) {
        x1 = 0.0;
        x2 = 0.0;
        if (b1 >= 0.0 && b2 >= 0.0) found = true;
    }
    if (!found) return;  // degenerate; keep last iterate

    const double d1 = x1 - a1, d2 = x2 - a2;
    apply_impulse(a, b, c1.ra, c1.rb, d1 * n);
    apply_impulse(a, b, c2.ra, c2.rb, d2 * n);
    c1.normal_impulse = x1;
    c2.normal_impulse = x2;
}

}  // namespace

World::World(const Settings& s) : settings(s) {}

int World::add_body(const BodyDef& def, const Shape& shape) {
    Body b;
    b.shape = shape;
    b.x = def.position;
    b.q = def.rotation;
    b.friction = def.friction;
    b.restitution = def.restitution;
    if (!def.is_static) {
        if (!(def.density > 0.0)) throw std::invalid_argument("dynamic body needs positive density");
        const MassData md = compute_mass(shape, def.density);
        b.mass = md.mass;
        b.inv_mass = 1.0 / md.mass;
        b.inertia = md.inertia;
        b.inv_inertia = 1.0 / md.inertia;
        b.v = def.velocity;
        b.w = def.angular_velocity;
    }
    bodies.push_back(b);
    return static_cast<int>(bodies.size()) - 1;
}

int World::add_distance_joint(int a, int b, Vec2 anchor_a, Vec2 anchor_b) {
    if (a == b || a < 0 || b < 0 || a >= static_cast<int>(bodies.size()) ||
        b >= static_cast<int>(bodies.size())) {
        throw std::invalid_argument("joint needs two distinct existing bodies");
    }
    DistanceJoint j;
    j.a = a;
    j.b = b;
    j.local_a = apply_inv(bodies[a].transform(), anchor_a);
    j.local_b = apply_inv(bodies[b].transform(), anchor_b);
    j.length = length(anchor_b - anchor_a);
    distance_joints.push_back(j);
    jointed_.insert(std::upper_bound(jointed_.begin(), jointed_.end(), pair_key(a, b)),
                    pair_key(a, b));
    return static_cast<int>(distance_joints.size()) - 1;
}

int World::add_revolute_joint(int a, int b, Vec2 anchor) {
    if (a == b || a < 0 || b < 0 || a >= static_cast<int>(bodies.size()) ||
        b >= static_cast<int>(bodies.size())) {
        throw std::invalid_argument("joint needs two distinct existing bodies");
    }
    RevoluteJoint j;
    j.a = a;
    j.b = b;
    j.local_a = apply_inv(bodies[a].transform(), anchor);
    j.local_b = apply_inv(bodies[b].transform(), anchor);
    revolute_joints.push_back(j);
    jointed_.insert(std::upper_bound(jointed_.begin(), jointed_.end(), pair_key(a, b)),
                    pair_key(a, b));
    return static_cast<int>(revolute_joints.size()) - 1;
}

void World::step(double dt) {
    if (!(dt > 0.0)) throw std::invalid_argument("dt must be positive");
    const auto t0 = Clock::now();

    // Semi-implicit (symplectic) Euler: velocities first, positions later
    // from the solved velocities.
    for (Body& b : bodies) {
        if (b.is_static()) continue;
        b.v += dt * (settings.gravity + b.inv_mass * b.force);
        b.w += dt * b.inv_inertia * b.torque;
        b.force = {0.0, 0.0};
        b.torque = 0.0;
    }

    find_pairs(dt);
    const auto t1 = Clock::now();
    narrowphase();
    const auto t2 = Clock::now();
    prepare(dt);
    solver_.resize(bodies.size());
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        const Body& b = bodies[i];
        solver_[i] = {b.v, b.w, b.inv_mass, b.inv_inertia, 0.0};
    }
    warm_start();
    for (int i = 0; i < settings.velocity_iterations; ++i) solve_velocities();
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        bodies[i].v = solver_[i].v;
        bodies[i].w = solver_[i].w;
    }
    const auto t3 = Clock::now();
    integrate_positions(dt);
    const auto t4 = Clock::now();

    timings.broadphase_s += seconds_since(t0, t1);
    timings.narrowphase_s += seconds_since(t1, t2);
    timings.solve_s += seconds_since(t2, t3);
    timings.integrate_s += seconds_since(t3, t4);
    timings.steps += 1;
}

void World::find_pairs(double dt) {
    const int n = static_cast<int>(bodies.size());
    aabbs_.resize(n);
    margins_.resize(n);
    for (int i = 0; i < n; ++i) {
        const Body& b = bodies[i];
        // Distance this body's surface can travel in one step.
        double margin = 0.0;
        if (settings.speculative && !b.is_static()) {
            const double w = b.w < 0.0 ? -b.w : b.w;
            margin = dt * (length(b.v) + w * b.shape.bounding_radius);
        }
        margins_[i] = margin;
        AABB box = compute_aabb(b.shape, b.transform());
        const Vec2 grow{margin + settings.linear_slop, margin + settings.linear_slop};
        box.lower -= grow;
        box.upper += grow;
        aabbs_[i] = box;
    }

    pairs_.clear();
    auto consider = [&](int a, int b) {
        if (a > b) std::swap(a, b);
        if (bodies[a].is_static() && bodies[b].is_static()) return;
        if (!overlaps(aabbs_[a], aabbs_[b])) return;
        if (std::binary_search(jointed_.begin(), jointed_.end(), pair_key(a, b))) return;
        pairs_.push_back({a, b, margins_[a] + margins_[b] + settings.linear_slop});
    };

    if (settings.broadphase == Broadphase::BruteForce) {
        for (int a = 0; a < n; ++a) {
            for (int b = a + 1; b < n; ++b) consider(a, b);
        }
    } else {
        // Sweep and prune on x. Ties break on body index, so the order is a
        // total order and any correct sort gives the same result.
        auto before = [&](int i, int j) {
            if (aabbs_[i].lower.x != aabbs_[j].lower.x) return aabbs_[i].lower.x < aabbs_[j].lower.x;
            return i < j;
        };
        if (static_cast<int>(order_.size()) != n) {
            order_.resize(n);
            for (int i = 0; i < n; ++i) order_[i] = i;
            std::sort(order_.begin(), order_.end(), before);
        } else {
            // Bodies move little per step, so last step's order is nearly
            // sorted and insertion sort runs in close to linear time.
            for (int s = 1; s < n; ++s) {
                const int id = order_[s];
                int t = s - 1;
                while (t >= 0 && before(id, order_[t])) {
                    order_[t + 1] = order_[t];
                    --t;
                }
                order_[t + 1] = id;
            }
        }
        // Copy the boxes into sweep order so the inner loop reads contiguous
        // memory, and test the cheap overlap before touching anything else.
        sorted_.resize(n);
        for (int s = 0; s < n; ++s) {
            const int i = order_[s];
            sorted_[s] = {aabbs_[i], i, bodies[i].is_static()};
        }
        for (int s = 0; s < n; ++s) {
            const SweepBox& bi = sorted_[s];
            const double hi = bi.box.upper.x;
            for (int t = s + 1; t < n; ++t) {
                const SweepBox& bj = sorted_[t];
                if (bj.box.lower.x > hi) break;
                if (bj.box.lower.y > bi.box.upper.y || bi.box.lower.y > bj.box.upper.y) continue;
                if (bi.is_static && bj.is_static) continue;
                const int a = bi.id < bj.id ? bi.id : bj.id;
                const int b = bi.id < bj.id ? bj.id : bi.id;
                if (!jointed_.empty() &&
                    std::binary_search(jointed_.begin(), jointed_.end(), pair_key(a, b))) {
                    continue;
                }
                pairs_.push_back({a, b, margins_[a] + margins_[b] + settings.linear_slop});
            }
        }
        // Same pair order as brute force, so the solver sees identical input.
        std::sort(pairs_.begin(), pairs_.end(), [](const Pair& p, const Pair& q) {
            return p.a != q.a ? p.a < q.a : p.b < q.b;
        });
    }
}

void World::narrowphase() {
    std::swap(previous_, manifolds_);
    manifolds_.clear();
    std::size_t cursor = 0;  // previous_ and pairs_ are both sorted by (a, b)
    for (const Pair& pair : pairs_) {
        const Body& ba = bodies[pair.a];
        const Body& bb = bodies[pair.b];
        const LocalManifold lm =
            collide(ba.shape, ba.transform(), bb.shape, bb.transform(), pair.max_separation);
        if (lm.count == 0) continue;

        Manifold m;
        m.a = pair.a;
        m.b = pair.b;
        m.normal = lm.normal;
        m.count = lm.count;
        m.friction = std::sqrt(ba.friction * bb.friction);
        m.restitution = ba.restitution > bb.restitution ? ba.restitution : bb.restitution;

        while (cursor < previous_.size() &&
               pair_key(previous_[cursor].a, previous_[cursor].b) < pair_key(m.a, m.b)) {
            ++cursor;
        }
        const Manifold* old = nullptr;
        if (cursor < previous_.size() && previous_[cursor].a == m.a && previous_[cursor].b == m.b) {
            old = &previous_[cursor];
        }

        for (int i = 0; i < lm.count; ++i) {
            ContactPoint& cp = m.points[i];
            cp.point = lm.points[i].point;
            cp.separation = lm.points[i].separation;
            cp.id = lm.points[i].id;
            if (old && settings.warm_starting) {
                for (int k = 0; k < old->count; ++k) {
                    if (old->points[k].id == cp.id) {
                        cp.normal_impulse = old->points[k].normal_impulse;
                        cp.tangent_impulse = old->points[k].tangent_impulse;
                        break;
                    }
                }
            }
        }
        manifolds_.push_back(m);
    }
}

void World::prepare(double dt) {
    const double inv_dt = 1.0 / dt;
    const Settings& s = settings;

    for (Manifold& m : manifolds_) {
        const Body& a = bodies[m.a];
        const Body& b = bodies[m.b];
        const Vec2 n = m.normal;
        const Vec2 t{n.y, -n.x};
        for (int i = 0; i < m.count; ++i) {
            ContactPoint& cp = m.points[i];
            cp.ra = cp.point - a.x;
            cp.rb = cp.point - b.x;
            const double rna = cross(cp.ra, n);
            const double rnb = cross(cp.rb, n);
            const double kn = a.inv_mass + b.inv_mass + a.inv_inertia * rna * rna +
                              b.inv_inertia * rnb * rnb;
            cp.normal_mass = 1.0 / kn;
            const double rta = cross(cp.ra, t);
            const double rtb = cross(cp.rb, t);
            const double kt = a.inv_mass + b.inv_mass + a.inv_inertia * rta * rta +
                              b.inv_inertia * rtb * rtb;
            cp.tangent_mass = 1.0 / kt;

            // Normal velocity before any impulse this step; restitution uses
            // this so the bounce is not lost when a speculative contact has
            // already slowed the bodies down.
            const double vn0 = dot(velocity_at(b, cp.rb) - velocity_at(a, cp.ra), n);
            const bool bounces = m.restitution > 0.0 && vn0 < -s.restitution_threshold;
            const double sep = cp.separation;
            if (sep > 0.0) {
                if (bounces && sep + vn0 * dt < 0.0) {
                    cp.target_velocity = -m.restitution * vn0;  // impact lands inside this step
                } else if (m.restitution > 0.0) {
                    // Elastic contact not expected to close this step: leave it
                    // free. If another impulse knocks the bodies together anyway,
                    // capping the approach here would stop them dead and lose the
                    // bounce; next step it is a touching contact and bounces.
                    cp.target_velocity = -std::numeric_limits<double>::infinity();
                } else {
                    cp.target_velocity = -sep * inv_dt;  // may close the gap, never more
                }
            } else {
                const double push = sep + s.linear_slop < 0.0 ? -(sep + s.linear_slop) : 0.0;
                cp.target_velocity = s.baumgarte * inv_dt * push;
                if (bounces && -m.restitution * vn0 > cp.target_velocity) {
                    cp.target_velocity = -m.restitution * vn0;
                }
            }
        }

        m.block = false;
        if (m.count == 2) {
            const ContactPoint& c1 = m.points[0];
            const ContactPoint& c2 = m.points[1];
            const double rn1a = cross(c1.ra, n), rn1b = cross(c1.rb, n);
            const double rn2a = cross(c2.ra, n), rn2b = cross(c2.rb, n);
            const double ms = a.inv_mass + b.inv_mass;
            m.k11 = ms + a.inv_inertia * rn1a * rn1a + b.inv_inertia * rn1b * rn1b;
            m.k22 = ms + a.inv_inertia * rn2a * rn2a + b.inv_inertia * rn2b * rn2b;
            m.k12 = ms + a.inv_inertia * rn1a * rn2a + b.inv_inertia * rn1b * rn2b;
            // Use the block solver only when K is well conditioned.
            // Use the block solver only when K is well conditioned and both
            // points are constrained (a free point has an infinite target).
            const double lowest = -std::numeric_limits<double>::infinity();
            m.block = m.k11 * m.k11 < 1000.0 * (m.k11 * m.k22 - m.k12 * m.k12) &&
                      c1.target_velocity != lowest && c2.target_velocity != lowest;
        }
    }

    for (DistanceJoint& j : distance_joints) {
        const Body& a = bodies[j.a];
        const Body& b = bodies[j.b];
        j.ra = rotate(a.q, j.local_a);
        j.rb = rotate(b.q, j.local_b);
        const Vec2 d = (b.x + j.rb) - (a.x + j.ra);
        const double len = length(d);
        j.u = len > 1e-12 ? (1.0 / len) * d : Vec2{1.0, 0.0};
        const double cra = cross(j.ra, j.u);
        const double crb = cross(j.rb, j.u);
        const double k = a.inv_mass + b.inv_mass + a.inv_inertia * cra * cra +
                         b.inv_inertia * crb * crb;
        j.mass = k > 0.0 ? 1.0 / k : 0.0;
        j.bias = s.baumgarte * inv_dt * (len - j.length);
        if (!s.warm_starting) j.impulse = 0.0;
    }

    for (RevoluteJoint& j : revolute_joints) {
        const Body& a = bodies[j.a];
        const Body& b = bodies[j.b];
        j.ra = rotate(a.q, j.local_a);
        j.rb = rotate(b.q, j.local_b);
        const double ma = a.inv_mass, mb = b.inv_mass;
        const double ia = a.inv_inertia, ib = b.inv_inertia;
        j.k11 = ma + mb + ia * j.ra.y * j.ra.y + ib * j.rb.y * j.rb.y;
        j.k12 = -ia * j.ra.x * j.ra.y - ib * j.rb.x * j.rb.y;
        j.k22 = ma + mb + ia * j.ra.x * j.ra.x + ib * j.rb.x * j.rb.x;
        j.bias = (s.baumgarte * inv_dt) * ((b.x + j.rb) - (a.x + j.ra));
        if (!s.warm_starting) j.impulse = {0.0, 0.0};
    }
}

void World::warm_start() {
    for (Manifold& m : manifolds_) {
        SolverBody& a = solver_[m.a];
        SolverBody& b = solver_[m.b];
        const Vec2 t{m.normal.y, -m.normal.x};
        for (int i = 0; i < m.count; ++i) {
            const ContactPoint& cp = m.points[i];
            apply_impulse(a, b, cp.ra, cp.rb, cp.normal_impulse * m.normal + cp.tangent_impulse * t);
        }
    }
    for (DistanceJoint& j : distance_joints) {
        apply_impulse(solver_[j.a], solver_[j.b], j.ra, j.rb, j.impulse * j.u);
    }
    for (RevoluteJoint& j : revolute_joints) {
        apply_impulse(solver_[j.a], solver_[j.b], j.ra, j.rb, j.impulse);
    }
}

void World::solve_velocities() {
    for (DistanceJoint& j : distance_joints) {
        SolverBody& a = solver_[j.a];
        SolverBody& b = solver_[j.b];
        const double cdot = dot(velocity_at(b, j.rb) - velocity_at(a, j.ra), j.u);
        const double lambda = -j.mass * (cdot + j.bias);
        j.impulse += lambda;
        apply_impulse(a, b, j.ra, j.rb, lambda * j.u);
    }

    for (RevoluteJoint& j : revolute_joints) {
        SolverBody& a = solver_[j.a];
        SolverBody& b = solver_[j.b];
        const Vec2 rhs = velocity_at(b, j.rb) - velocity_at(a, j.ra) + j.bias;
        // Solve K * impulse = -rhs for the symmetric 2x2 K.
        const double det = j.k11 * j.k22 - j.k12 * j.k12;
        const double inv_det = det != 0.0 ? 1.0 / det : 0.0;
        const Vec2 p{-inv_det * (j.k22 * rhs.x - j.k12 * rhs.y),
                     -inv_det * (j.k11 * rhs.y - j.k12 * rhs.x)};
        j.impulse += p;
        apply_impulse(a, b, j.ra, j.rb, p);
    }

    for (Manifold& m : manifolds_) {
        SolverBody& a = solver_[m.a];
        SolverBody& b = solver_[m.b];
        const Vec2 n = m.normal;
        const Vec2 t{n.y, -n.x};

        // Friction first, bounded by the current normal impulse (Coulomb cone).
        for (int i = 0; i < m.count; ++i) {
            ContactPoint& cp = m.points[i];
            const double vt = dot(velocity_at(b, cp.rb) - velocity_at(a, cp.ra), t);
            double lambda = -cp.tangent_mass * vt;
            const double max_f = m.friction * cp.normal_impulse;
            const double old = cp.tangent_impulse;
            double next = old + lambda;
            if (next > max_f) next = max_f;
            if (next < -max_f) next = -max_f;
            cp.tangent_impulse = next;
            lambda = next - old;
            apply_impulse(a, b, cp.ra, cp.rb, lambda * t);
        }

        if (m.block) {
            solve_block(m, a, b);
            continue;
        }

        // Non-penetration: accumulated impulse stays non-negative (push only).
        for (int i = 0; i < m.count; ++i) {
            ContactPoint& cp = m.points[i];
            const double vn = dot(velocity_at(b, cp.rb) - velocity_at(a, cp.ra), n);
            double lambda = -cp.normal_mass * (vn - cp.target_velocity);
            const double old = cp.normal_impulse;
            const double next = old + lambda > 0.0 ? old + lambda : 0.0;
            cp.normal_impulse = next;
            lambda = next - old;
            apply_impulse(a, b, cp.ra, cp.rb, lambda * n);
        }
    }
}

void World::integrate_positions(double dt) {
    for (Body& b : bodies) {
        if (b.is_static()) continue;
        b.x += dt * b.v;
        b.q = integrate_rotation(b.q, dt * b.w);
    }
}

std::uint64_t World::state_hash() const {
    std::uint64_t h = 1469598103934665603ull;
    auto mix = [&h](double value) {
        std::uint64_t bits;
        std::memcpy(&bits, &value, sizeof bits);
        for (int i = 0; i < 8; ++i) {
            h ^= (bits >> (8 * i)) & 0xffu;
            h *= 1099511628211ull;
        }
    };
    for (const Body& b : bodies) {
        mix(b.x.x);
        mix(b.x.y);
        mix(b.q.c);
        mix(b.q.s);
        mix(b.v.x);
        mix(b.v.y);
        mix(b.w);
    }
    return h;
}

double World::kinetic_energy() const {
    double e = 0.0;
    for (const Body& b : bodies) e += 0.5 * (b.mass * dot(b.v, b.v) + b.inertia * b.w * b.w);
    return e;
}

double World::potential_energy() const {
    double e = 0.0;
    for (const Body& b : bodies) e -= b.mass * dot(settings.gravity, b.x);
    return e;
}

Vec2 World::linear_momentum() const {
    Vec2 p{0.0, 0.0};
    for (const Body& b : bodies) p += b.mass * b.v;
    return p;
}

double World::angular_momentum() const {
    double l = 0.0;
    for (const Body& b : bodies) l += b.inertia * b.w + b.mass * cross(b.x, b.v);
    return l;
}

}  // namespace rbsim
