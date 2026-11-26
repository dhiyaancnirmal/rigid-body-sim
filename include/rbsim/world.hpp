#pragma once

#include <cstdint>
#include <vector>

#include "rbsim/collide.hpp"
#include "rbsim/shape.hpp"

namespace rbsim {

enum class Broadphase : std::uint8_t { BruteForce, SweepAndPrune };

struct Settings {
    Vec2 gravity{0.0, -9.81};
    int velocity_iterations = 8;
    Broadphase broadphase = Broadphase::SweepAndPrune;
    bool speculative = true;   // keep contacts that may close within the step
    bool warm_starting = true;
    double baumgarte = 0.2;    // fraction of penetration removed per step
    double linear_slop = 0.005;  // allowed penetration, metres
    double restitution_threshold = 0.5;  // m/s; slower impacts do not bounce
};

struct BodyDef {
    Vec2 position;  // of the shape's centroid
    Rot rotation;
    Vec2 velocity;
    double angular_velocity = 0.0;
    double density = 1.0;
    double friction = 0.5;
    double restitution = 0.0;
    bool is_static = false;
};

struct Body {
    Shape shape;
    Vec2 x;
    Rot q;
    Vec2 v;
    double w = 0.0;
    Vec2 force;
    double torque = 0.0;
    double mass = 0.0, inv_mass = 0.0;
    double inertia = 0.0, inv_inertia = 0.0;
    double friction = 0.5;
    double restitution = 0.0;

    Transform transform() const { return {x, q}; }
    bool is_static() const { return inv_mass == 0.0; }
};

// The four numbers the velocity solver touches, packed into 48 bytes. The
// solver iterates over these instead of the full Body (about 400 bytes with
// its shape), so a cache line holds a body's whole solver state.
struct SolverBody {
    Vec2 v;
    double w = 0.0;
    double inv_mass = 0.0;
    double inv_inertia = 0.0;
    double pad_ = 0.0;
};

struct ContactPoint {
    Vec2 point;
    double separation = 0.0;
    std::uint32_t id = 0;
    Vec2 ra, rb;
    double normal_mass = 0.0, tangent_mass = 0.0;
    double normal_impulse = 0.0, tangent_impulse = 0.0;
    double target_velocity = 0.0;  // the solver keeps normal velocity >= this
};

struct Manifold {
    int a = 0, b = 0;  // a < b
    Vec2 normal;       // from a to b
    int count = 0;
    ContactPoint points[2];
    double friction = 0.0;
    double restitution = 0.0;
    // Two-point manifolds solve both normal impulses together (a 2x2 LCP),
    // which keeps stacks symmetric. K is the effective mass matrix.
    bool block = false;
    double k11 = 0.0, k12 = 0.0, k22 = 0.0;
};

// Keeps the distance between two anchor points fixed (a rigid rod).
struct DistanceJoint {
    int a = 0, b = 0;
    Vec2 local_a, local_b;
    double length = 0.0;
    Vec2 u, ra, rb;
    double mass = 0.0, bias = 0.0, impulse = 0.0;
};

// Pins two bodies together at a shared point (a hinge).
struct RevoluteJoint {
    int a = 0, b = 0;
    Vec2 local_a, local_b;
    Vec2 ra, rb;
    double k11 = 0.0, k12 = 0.0, k22 = 0.0;  // effective mass matrix (symmetric)
    Vec2 bias;
    Vec2 impulse;
};

// Wall-clock time spent in each phase of step(), summed over all steps.
struct StepTimings {
    double broadphase_s = 0.0;
    double narrowphase_s = 0.0;
    double solve_s = 0.0;
    double integrate_s = 0.0;
    std::int64_t steps = 0;
};

class World {
public:
    explicit World(const Settings& settings = Settings{});

    int add_body(const BodyDef& def, const Shape& shape);
    // Anchors are given in world coordinates at the current pose.
    int add_distance_joint(int a, int b, Vec2 anchor_a, Vec2 anchor_b);
    int add_revolute_joint(int a, int b, Vec2 anchor);

    void step(double dt);

    // FNV-1a over the exact bit patterns of every body's pose and velocity.
    std::uint64_t state_hash() const;

    double kinetic_energy() const;
    double potential_energy() const;  // gravity, relative to y = 0
    Vec2 linear_momentum() const;
    double angular_momentum() const;  // about the world origin

    Settings settings;
    std::vector<Body> bodies;
    std::vector<DistanceJoint> distance_joints;
    std::vector<RevoluteJoint> revolute_joints;
    const std::vector<Manifold>& manifolds() const { return manifolds_; }
    std::size_t pair_count() const { return pairs_.size(); }
    StepTimings timings;

private:
    void find_pairs(double dt);
    void narrowphase();
    void prepare(double dt);
    void warm_start();
    void solve_velocities();
    void integrate_positions(double dt);

    struct Pair {
        int a, b;
        double max_separation;
    };
    std::vector<AABB> aabbs_;
    std::vector<double> margins_;
    std::vector<int> order_;  // body indices in sweep order, kept between steps
    struct SweepBox {
        AABB box;
        int id;
        bool is_static;
    };
    std::vector<SweepBox> sorted_;
    std::vector<Pair> pairs_;
    std::vector<std::uint64_t> jointed_;  // sorted body-pair keys that never collide
    std::vector<Manifold> manifolds_;
    std::vector<Manifold> previous_;
    std::vector<SolverBody> solver_;
};

}  // namespace rbsim
