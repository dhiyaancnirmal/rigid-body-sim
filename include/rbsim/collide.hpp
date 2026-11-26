#pragma once

#include <cstdint>

#include "rbsim/shape.hpp"

namespace rbsim {

struct ManifoldPoint {
    Vec2 point;         // world, midway between the two surfaces
    double separation;  // negative when penetrating
    std::uint32_t id;   // feature key, stable across steps for warm starting
};

struct LocalManifold {
    Vec2 normal;  // world, points from shape A to shape B
    int count = 0;
    ManifoldPoint points[2];
};

// Contact points between two shapes. Points whose separation exceeds
// `max_separation` are dropped; a positive value keeps "speculative" points
// that are not touching yet but could be by the end of the step.
LocalManifold collide(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb,
                      double max_separation);

}  // namespace rbsim
