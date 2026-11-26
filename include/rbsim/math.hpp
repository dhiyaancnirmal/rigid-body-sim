// Small 2D math kit. Everything here uses only IEEE-754 basic operations and
// sqrt, which are correctly rounded on every conforming platform. That is what
// lets a simulation replay bit-for-bit across machines: no libm sin/cos/atan
// runs inside World::step.
#pragma once

#include <cmath>

namespace rbsim {

struct Vec2 {
    double x = 0.0;
    double y = 0.0;
};

inline Vec2 operator+(Vec2 a, Vec2 b) { return {a.x + b.x, a.y + b.y}; }
inline Vec2 operator-(Vec2 a, Vec2 b) { return {a.x - b.x, a.y - b.y}; }
inline Vec2 operator-(Vec2 a) { return {-a.x, -a.y}; }
inline Vec2 operator*(double s, Vec2 a) { return {s * a.x, s * a.y}; }
inline Vec2& operator+=(Vec2& a, Vec2 b) { a.x += b.x; a.y += b.y; return a; }
inline Vec2& operator-=(Vec2& a, Vec2 b) { a.x -= b.x; a.y -= b.y; return a; }

inline double dot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }
inline double cross(Vec2 a, Vec2 b) { return a.x * b.y - a.y * b.x; }
// cross(w, v) for scalar angular velocity w: w * perp(v)
inline Vec2 cross(double w, Vec2 v) { return {-w * v.y, w * v.x}; }
inline Vec2 cross(Vec2 v, double s) { return {s * v.y, -s * v.x}; }
inline Vec2 perp(Vec2 v) { return {-v.y, v.x}; }
inline double length_sq(Vec2 v) { return dot(v, v); }
inline double length(Vec2 v) { return std::sqrt(dot(v, v)); }
inline Vec2 normalize(Vec2 v) {
    const double len = length(v);
    if (len < 1e-300) return {0.0, 0.0};
    return (1.0 / len) * v;
}
inline Vec2 vmin(Vec2 a, Vec2 b) { return {a.x < b.x ? a.x : b.x, a.y < b.y ? a.y : b.y}; }
inline Vec2 vmax(Vec2 a, Vec2 b) { return {a.x > b.x ? a.x : b.x, a.y > b.y ? a.y : b.y}; }

// Rotation stored as a unit complex number (cos, sin). Rotating, inverse
// rotating and integrating it need no trig functions.
struct Rot {
    double c = 1.0;
    double s = 0.0;
};

inline Vec2 rotate(Rot q, Vec2 v) { return {q.c * v.x - q.s * v.y, q.s * v.x + q.c * v.y}; }
inline Vec2 inv_rotate(Rot q, Vec2 v) { return {q.c * v.x + q.s * v.y, -q.s * v.x + q.c * v.y}; }

// First-order rotor update, as in Box2D v3 (b2IntegrateRotation): multiply by
// (1, w*dt) and renormalize. The body turns by atan(w*dt) per step instead of
// w*dt; the relative error is (w*dt)^2/3, about 1e-4 at 1 rad per 60 Hz step.
inline Rot integrate_rotation(Rot q, double delta_angle) {
    Rot r{q.c - delta_angle * q.s, q.s + delta_angle * q.c};
    const double inv = 1.0 / std::sqrt(r.c * r.c + r.s * r.s);
    return {r.c * inv, r.s * inv};
}

struct Transform {
    Vec2 p;
    Rot q;
};

inline Vec2 apply(const Transform& xf, Vec2 v) { return rotate(xf.q, v) + xf.p; }
inline Vec2 apply_inv(const Transform& xf, Vec2 v) { return inv_rotate(xf.q, v - xf.p); }

}  // namespace rbsim
