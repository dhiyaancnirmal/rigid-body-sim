#include "rbsim/collide.hpp"

namespace rbsim {
namespace {

LocalManifold circle_circle(const Shape& a, const Transform& xa, const Shape& b,
                            const Transform& xb, double max_sep) {
    LocalManifold m;
    const Vec2 d = xb.p - xa.p;
    const double dist = length(d);
    const double sep = dist - a.radius - b.radius;
    if (sep > max_sep) return m;
    const Vec2 n = dist > 1e-12 ? (1.0 / dist) * d : Vec2{0.0, 1.0};
    m.normal = n;
    m.count = 1;
    m.points[0] = {xa.p + (a.radius + 0.5 * sep) * n, sep, 0u};
    return m;
}

// Polygon A against circle B. Normal points from A to B.
LocalManifold polygon_circle(const Shape& poly, const Transform& xa, const Shape& circle,
                             const Transform& xb, double max_sep) {
    LocalManifold m;
    const Vec2 c = apply_inv(xa, xb.p);  // circle centre in polygon frame
    const double r = circle.radius;

    int face = 0;
    double s = -1e300;
    for (int i = 0; i < poly.count; ++i) {
        const double si = dot(poly.normals[i], c - poly.vertices[i]);
        if (si > s) {
            s = si;
            face = i;
        }
    }
    if (s - r > max_sep) return m;

    const Vec2 v1 = poly.vertices[face];
    const Vec2 v2 = poly.vertices[(face + 1) % poly.count];
    const double u1 = dot(c - v1, v2 - v1);
    const double u2 = dot(c - v2, v1 - v2);

    Vec2 n_local;
    double sep;
    std::uint32_t id;
    if (u1 < 0.0 && s > 1e-12) {
        // Vertex region of v1.
        const Vec2 d = c - v1;
        const double dist = length(d);
        sep = dist - r;
        if (sep > max_sep) return m;
        n_local = (1.0 / dist) * d;
        id = 0x100u | static_cast<std::uint32_t>(face);
    } else if (u2 < 0.0 && s > 1e-12) {
        const Vec2 d = c - v2;
        const double dist = length(d);
        sep = dist - r;
        if (sep > max_sep) return m;
        n_local = (1.0 / dist) * d;
        id = 0x100u | static_cast<std::uint32_t>((face + 1) % poly.count);
    } else {
        n_local = poly.normals[face];
        sep = s - r;
        id = static_cast<std::uint32_t>(face);
    }

    const Vec2 on_circle = c - r * n_local;
    const Vec2 mid = on_circle - (0.5 * sep) * n_local;
    m.normal = rotate(xa.q, n_local);
    m.count = 1;
    m.points[0] = {apply(xa, mid), sep, id};
    return m;
}

struct WorldPoly {
    int count;
    Vec2 v[kMaxPolygonVertices];
    Vec2 n[kMaxPolygonVertices];
};

WorldPoly to_world(const Shape& s, const Transform& xf) {
    WorldPoly w;
    w.count = s.count;
    for (int i = 0; i < s.count; ++i) {
        w.v[i] = apply(xf, s.vertices[i]);
        w.n[i] = rotate(xf.q, s.normals[i]);
    }
    return w;
}

// Largest separation of B from any face of A (separating axis test).
double max_separation(const WorldPoly& a, const WorldPoly& b, int* edge) {
    double best = -1e300;
    int best_i = 0;
    for (int i = 0; i < a.count; ++i) {
        double smin = 1e300;
        for (int j = 0; j < b.count; ++j) {
            const double sj = dot(a.n[i], b.v[j] - a.v[i]);
            if (sj < smin) smin = sj;
        }
        if (smin > best) {
            best = smin;
            best_i = i;
        }
    }
    *edge = best_i;
    return best;
}

struct ClipVertex {
    Vec2 v;
    std::uint32_t id;
};

// Keep the part of segment `in` where dot(normal, p) <= offset.
int clip_segment(ClipVertex out[2], const ClipVertex in[2], Vec2 normal, double offset,
                 std::uint32_t clip_id) {
    int n = 0;
    const double d0 = dot(normal, in[0].v) - offset;
    const double d1 = dot(normal, in[1].v) - offset;
    if (d0 <= 0.0) out[n++] = in[0];
    if (d1 <= 0.0) out[n++] = in[1];
    if (d0 * d1 < 0.0) {
        const double t = d0 / (d0 - d1);
        out[n].v = in[0].v + t * (in[1].v - in[0].v);
        out[n].id = clip_id;
        ++n;
    }
    return n;
}

LocalManifold polygon_polygon(const Shape& sa, const Transform& xa, const Shape& sb,
                              const Transform& xb, double max_sep) {
    LocalManifold m;
    const WorldPoly a = to_world(sa, xa);
    const WorldPoly b = to_world(sb, xb);

    int edge_a = 0;
    const double sep_a = max_separation(a, b, &edge_a);
    if (sep_a > max_sep) return m;
    int edge_b = 0;
    const double sep_b = max_separation(b, a, &edge_b);
    if (sep_b > max_sep) return m;

    // Prefer A as the reference face unless B is clearly better; the bias keeps
    // the choice from flickering between steps.
    const bool flip = sep_b > sep_a + 0.0005;
    const WorldPoly& ref = flip ? b : a;
    const WorldPoly& inc = flip ? a : b;
    const int ref_edge = flip ? edge_b : edge_a;
    const Vec2 n = ref.n[ref_edge];

    // Incident edge: the face of `inc` most anti-parallel to the reference normal.
    int inc_edge = 0;
    double min_dot = 1e300;
    for (int i = 0; i < inc.count; ++i) {
        const double d = dot(n, inc.n[i]);
        if (d < min_dot) {
            min_dot = d;
            inc_edge = i;
        }
    }
    const int i1 = inc_edge;
    const int i2 = (inc_edge + 1) % inc.count;
    const ClipVertex incident[2] = {{inc.v[i1], static_cast<std::uint32_t>(i1)},
                                    {inc.v[i2], static_cast<std::uint32_t>(i2)}};

    const Vec2 v1 = ref.v[ref_edge];
    const Vec2 v2 = ref.v[(ref_edge + 1) % ref.count];
    const Vec2 t = normalize(v2 - v1);

    ClipVertex clip1[2];
    ClipVertex clip2[2];
    if (clip_segment(clip1, incident, -t, -dot(t, v1), 0x40u) < 2) return m;
    if (clip_segment(clip2, clip1, t, dot(t, v2), 0x80u) < 2) return m;

    // Feature ids name the reference edge, the incident edge and the slot
    // (which end of the reference edge the point is nearer). They do not
    // depend on whether a point was clipped, so a box hanging a hair over the
    // edge below it keeps its warm-start impulse from step to step.
    const std::uint32_t base = (flip ? 0x10000u : 0u) |
                               (static_cast<std::uint32_t>(ref_edge) << 8) |
                               (static_cast<std::uint32_t>(inc_edge) << 4);
    const bool swap = dot(t, clip2[0].v) > dot(t, clip2[1].v);
    m.normal = flip ? -n : n;
    for (int k = 0; k < 2; ++k) {
        const ClipVertex& cv = clip2[swap ? 1 - k : k];
        const double sep = dot(n, cv.v - v1);
        if (sep > max_sep) continue;
        m.points[m.count++] = {cv.v - (0.5 * sep) * n, sep, base | static_cast<std::uint32_t>(k)};
    }
    return m;
}

}  // namespace

LocalManifold collide(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb,
                      double max_sep) {
    if (a.type == ShapeType::Circle && b.type == ShapeType::Circle) {
        return circle_circle(a, xa, b, xb, max_sep);
    }
    if (a.type == ShapeType::Polygon && b.type == ShapeType::Circle) {
        return polygon_circle(a, xa, b, xb, max_sep);
    }
    if (a.type == ShapeType::Circle && b.type == ShapeType::Polygon) {
        LocalManifold m = polygon_circle(b, xb, a, xa, max_sep);
        m.normal = -m.normal;
        return m;
    }
    return polygon_polygon(a, xa, b, xb, max_sep);
}

}  // namespace rbsim
