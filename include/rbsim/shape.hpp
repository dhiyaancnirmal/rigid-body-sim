#pragma once

#include <cstdint>
#include <vector>

#include "rbsim/math.hpp"

namespace rbsim {

constexpr int kMaxPolygonVertices = 8;

enum class ShapeType : std::uint8_t { Circle, Polygon };

// A shape in its body's local frame, with the centroid at the origin.
struct Shape {
    ShapeType type = ShapeType::Circle;
    double radius = 0.0;  // circles only
    int count = 0;        // polygons only
    Vec2 vertices[kMaxPolygonVertices];
    Vec2 normals[kMaxPolygonVertices];
    double area = 0.0;
    double unit_inertia = 0.0;  // polar moment about the centroid per unit density
    double bounding_radius = 0.0;
};

struct MassData {
    double mass = 0.0;
    double inertia = 0.0;
};

Shape make_circle(double radius);
Shape make_box(double half_width, double half_height);
// Vertices must be convex and counter-clockwise; 3..8 of them. The polygon is
// shifted so its centroid sits at the origin. Throws std::invalid_argument.
// If centroid_out is given, it receives the centroid of the input vertices.
Shape make_polygon(const std::vector<Vec2>& vertices, Vec2* centroid_out = nullptr);

MassData compute_mass(const Shape& shape, double density);

struct AABB {
    Vec2 lower;
    Vec2 upper;
};

AABB compute_aabb(const Shape& shape, const Transform& xf);

}  // namespace rbsim
