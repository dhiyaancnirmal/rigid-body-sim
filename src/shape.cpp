#include "rbsim/shape.hpp"

#include <stdexcept>

namespace rbsim {

Shape make_circle(double radius) {
    if (!(radius > 0.0)) throw std::invalid_argument("circle radius must be positive");
    Shape s;
    s.type = ShapeType::Circle;
    s.radius = radius;
    // Area and polar moment of a disk; pi is a constant, so no libm call.
    constexpr double kPi = 3.14159265358979323846;
    s.area = kPi * radius * radius;
    s.unit_inertia = 0.5 * s.area * radius * radius;
    s.bounding_radius = radius;
    return s;
}

Shape make_box(double hw, double hh) {
    if (!(hw > 0.0 && hh > 0.0)) throw std::invalid_argument("box half extents must be positive");
    return make_polygon({{-hw, -hh}, {hw, -hh}, {hw, hh}, {-hw, hh}});
}

Shape make_polygon(const std::vector<Vec2>& in, Vec2* centroid_out) {
    const int n = static_cast<int>(in.size());
    if (n < 3 || n > kMaxPolygonVertices) {
        throw std::invalid_argument("polygon needs 3 to 8 vertices");
    }
    for (int i = 0; i < n; ++i) {
        const Vec2 e1 = in[(i + 1) % n] - in[i];
        const Vec2 e2 = in[(i + 2) % n] - in[(i + 1) % n];
        if (!(cross(e1, e2) > 0.0)) {
            throw std::invalid_argument("polygon must be strictly convex and counter-clockwise");
        }
    }

    // Area, centroid and inertia by fan triangulation from the first vertex.
    const Vec2 origin = in[0];
    double area = 0.0;
    Vec2 c{0.0, 0.0};
    double inertia_origin = 0.0;  // about `origin`
    for (int i = 1; i + 1 < n; ++i) {
        const Vec2 e1 = in[i] - origin;
        const Vec2 e2 = in[i + 1] - origin;
        const double d = cross(e1, e2);
        const double tri_area = 0.5 * d;
        area += tri_area;
        c += (tri_area / 3.0) * (e1 + e2);
        const double intx2 = e1.x * e1.x + e2.x * e1.x + e2.x * e2.x;
        const double inty2 = e1.y * e1.y + e2.y * e1.y + e2.y * e2.y;
        inertia_origin += (0.25 / 3.0) * d * (intx2 + inty2);
    }
    c = (1.0 / area) * c;  // relative to origin

    Shape s;
    s.type = ShapeType::Polygon;
    s.count = n;
    s.area = area;
    // Parallel axis theorem: move from `origin` to the centroid.
    s.unit_inertia = inertia_origin - area * dot(c, c);
    const Vec2 centroid = origin + c;
    if (centroid_out) *centroid_out = centroid;
    double r2 = 0.0;
    for (int i = 0; i < n; ++i) {
        s.vertices[i] = in[i] - centroid;
        const double d2 = length_sq(s.vertices[i]);
        if (d2 > r2) r2 = d2;
    }
    for (int i = 0; i < n; ++i) {
        const Vec2 e = s.vertices[(i + 1) % n] - s.vertices[i];
        s.normals[i] = normalize(Vec2{e.y, -e.x});
    }
    s.bounding_radius = std::sqrt(r2);
    return s;
}

MassData compute_mass(const Shape& shape, double density) {
    return {density * shape.area, density * shape.unit_inertia};
}

AABB compute_aabb(const Shape& shape, const Transform& xf) {
    if (shape.type == ShapeType::Circle) {
        const Vec2 r{shape.radius, shape.radius};
        return {xf.p - r, xf.p + r};
    }
    Vec2 lo = apply(xf, shape.vertices[0]);
    Vec2 hi = lo;
    for (int i = 1; i < shape.count; ++i) {
        const Vec2 v = apply(xf, shape.vertices[i]);
        lo = vmin(lo, v);
        hi = vmax(hi, v);
    }
    return {lo, hi};
}

}  // namespace rbsim
