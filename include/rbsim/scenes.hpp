// Reference scenes shared by the tests and benchmarks. They use no libm calls
// during setup, so their state hashes do not depend on the platform's libm.
#pragma once

#include "rbsim/world.hpp"

namespace rbsim::scenes {

// State hash of pile(100) after 600 steps of 1/60 s. Every build with this
// repo's floating-point flags must reproduce it, from C++ and from Python;
// scripts/check_determinism.sh checks several optimization levels and
// instruction sets. If you change the physics on purpose, update it and say why in
// the commit.
constexpr std::uint64_t kPile100Hash600 = 0x778f89916d94c2deull;

// Static floor and two walls forming an open box `width` wide.
inline void add_container(World& w, double width, double height) {
    BodyDef d;
    d.is_static = true;
    d.friction = 0.6;
    d.position = {0.0, -0.5};
    w.add_body(d, make_box(0.5 * width + 1.0, 0.5));
    d.position = {-0.5 * width - 0.5, 0.5 * height};
    w.add_body(d, make_box(0.5, 0.5 * height));
    d.position = {0.5 * width + 0.5, 0.5 * height};
    w.add_body(d, make_box(0.5, 0.5 * height));
}

// `n` dynamic bodies (alternating boxes, circles and triangles) dropped on a
// grid into a container. Column count grows with sqrt(n) so the pile stays
// roughly square.
inline World pile(int n, Broadphase broadphase = Broadphase::SweepAndPrune) {
    Settings s;
    s.broadphase = broadphase;
    World w(s);
    int cols = 1;
    while (cols * cols < n) ++cols;
    const double spacing = 1.1;
    const double width = cols * spacing + 1.0;
    add_container(w, width, 4.0 + n / cols * spacing);
    const Shape box = make_box(0.4, 0.4);
    const Shape circle = make_circle(0.45);
    const Shape tri = make_polygon({{-0.45, -0.35}, {0.45, -0.35}, {0.0, 0.45}});
    for (int i = 0; i < n; ++i) {
        const int r = i / cols;
        const int c = i % cols;
        BodyDef d;
        // Odd rows are offset by a quarter spacing so the pile does not stack
        // in perfect columns.
        d.position = {-0.5 * width + 1.0 + c * spacing + (r % 2) * 0.25, 1.0 + r * spacing};
        d.friction = 0.5;
        d.angular_velocity = (i % 5 - 2) * 0.5;
        const Shape& sh = i % 3 == 0 ? box : (i % 3 == 1 ? circle : tri);
        w.add_body(d, sh);
    }
    return w;
}

// A vertical stack of `n` unit boxes on a static floor.
inline World stack(int n) {
    World w;
    BodyDef ground;
    ground.is_static = true;
    ground.position = {0.0, -0.5};
    w.add_body(ground, make_box(20.0, 0.5));
    for (int i = 0; i < n; ++i) {
        BodyDef d;
        d.position = {0.0, 0.5 + i * 1.0};
        d.friction = 0.6;
        w.add_body(d, make_box(0.5, 0.5));
    }
    return w;
}

}  // namespace rbsim::scenes
