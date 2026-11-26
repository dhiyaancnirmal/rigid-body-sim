#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <tuple>
#include <vector>

#include "rbsim/scenes.hpp"

namespace py = pybind11;
using namespace rbsim;

namespace {

// Angle 0 maps to the identity without calling libm, so scenes built without
// rotations stay bit-identical across platforms.
Rot rot_from_angle(double angle) {
    if (angle == 0.0) return {1.0, 0.0};
    return {std::cos(angle), std::sin(angle)};
}

BodyDef make_def(double x, double y, double angle, double vx, double vy, double w, double density,
                 double friction, double restitution, bool is_static) {
    BodyDef d;
    d.position = {x, y};
    d.rotation = rot_from_angle(angle);
    d.velocity = {vx, vy};
    d.angular_velocity = w;
    d.density = density;
    d.friction = friction;
    d.restitution = restitution;
    d.is_static = is_static;
    return d;
}

py::array_t<double> columns(const World& w, int cols, void (*fill)(const Body&, double*)) {
    const auto n = static_cast<py::ssize_t>(w.bodies.size());
    py::array_t<double> out({n, static_cast<py::ssize_t>(cols)});
    auto view = out.mutable_unchecked<2>();
    for (py::ssize_t i = 0; i < n; ++i) fill(w.bodies[i], view.mutable_data(i, 0));
    return out;
}

}  // namespace

PYBIND11_MODULE(_rbsim, m) {
    m.doc() = "Deterministic 2D rigid-body simulator (C++ core)";

    py::enum_<Broadphase>(m, "Broadphase")
        .value("BRUTE_FORCE", Broadphase::BruteForce)
        .value("SWEEP_AND_PRUNE", Broadphase::SweepAndPrune);

    py::class_<Settings>(m, "Settings")
        .def(py::init<>())
        .def_property(
            "gravity", [](const Settings& s) { return std::make_tuple(s.gravity.x, s.gravity.y); },
            [](Settings& s, std::tuple<double, double> g) {
                s.gravity = {std::get<0>(g), std::get<1>(g)};
            })
        .def_readwrite("velocity_iterations", &Settings::velocity_iterations)
        .def_readwrite("broadphase", &Settings::broadphase)
        .def_readwrite("speculative", &Settings::speculative)
        .def_readwrite("warm_starting", &Settings::warm_starting)
        .def_readwrite("baumgarte", &Settings::baumgarte)
        .def_readwrite("linear_slop", &Settings::linear_slop)
        .def_readwrite("restitution_threshold", &Settings::restitution_threshold);

    py::class_<World>(m, "World")
        .def(py::init<const Settings&>(), py::arg("settings") = Settings{})
        .def_readwrite("settings", &World::settings)
        .def(
            "add_circle",
            [](World& w, double x, double y, double radius, double angle, double vx, double vy,
               double av, double density, double friction, double restitution, bool is_static) {
                return w.add_body(make_def(x, y, angle, vx, vy, av, density, friction, restitution,
                                           is_static),
                                  make_circle(radius));
            },
            py::arg("x"), py::arg("y"), py::arg("radius"), py::kw_only(), py::arg("angle") = 0.0,
            py::arg("vx") = 0.0, py::arg("vy") = 0.0, py::arg("w") = 0.0,
            py::arg("density") = 1.0, py::arg("friction") = 0.5, py::arg("restitution") = 0.0,
            py::arg("static") = false)
        .def(
            "add_box",
            [](World& w, double x, double y, double hw, double hh, double angle, double vx,
               double vy, double av, double density, double friction, double restitution,
               bool is_static) {
                return w.add_body(make_def(x, y, angle, vx, vy, av, density, friction, restitution,
                                           is_static),
                                  make_box(hw, hh));
            },
            py::arg("x"), py::arg("y"), py::arg("half_width"), py::arg("half_height"),
            py::kw_only(), py::arg("angle") = 0.0, py::arg("vx") = 0.0, py::arg("vy") = 0.0,
            py::arg("w") = 0.0, py::arg("density") = 1.0, py::arg("friction") = 0.5,
            py::arg("restitution") = 0.0, py::arg("static") = false)
        .def(
            "add_polygon",
            [](World& w, const std::vector<std::tuple<double, double>>& verts, double vx,
               double vy, double av, double density, double friction, double restitution,
               bool is_static) {
                std::vector<Vec2> v;
                for (const auto& p : verts) v.push_back({std::get<0>(p), std::get<1>(p)});
                Vec2 centroid;
                const Shape shape = make_polygon(v, &centroid);
                return w.add_body(make_def(centroid.x, centroid.y, 0.0, vx, vy, av, density,
                                           friction, restitution, is_static),
                                  shape);
            },
            "Add a convex polygon from counter-clockwise world-space vertices.",
            py::arg("vertices"), py::kw_only(), py::arg("vx") = 0.0, py::arg("vy") = 0.0,
            py::arg("w") = 0.0, py::arg("density") = 1.0, py::arg("friction") = 0.5,
            py::arg("restitution") = 0.0, py::arg("static") = false)
        .def(
            "add_distance_joint",
            [](World& w, int a, int b, std::tuple<double, double> pa,
               std::tuple<double, double> pb) {
                return w.add_distance_joint(a, b, {std::get<0>(pa), std::get<1>(pa)},
                                            {std::get<0>(pb), std::get<1>(pb)});
            },
            py::arg("a"), py::arg("b"), py::arg("anchor_a"), py::arg("anchor_b"))
        .def(
            "add_revolute_joint",
            [](World& w, int a, int b, std::tuple<double, double> p) {
                return w.add_revolute_joint(a, b, {std::get<0>(p), std::get<1>(p)});
            },
            py::arg("a"), py::arg("b"), py::arg("anchor"))
        .def(
            "step",
            [](World& w, double dt, int n) {
                py::gil_scoped_release release;
                for (int i = 0; i < n; ++i) w.step(dt);
            },
            py::arg("dt"), py::arg("n") = 1)
        .def("state_hash", &World::state_hash)
        .def("kinetic_energy", &World::kinetic_energy)
        .def("potential_energy", &World::potential_energy)
        .def("linear_momentum",
             [](const World& w) {
                 const Vec2 p = w.linear_momentum();
                 return std::make_tuple(p.x, p.y);
             })
        .def("angular_momentum", &World::angular_momentum)
        .def_property_readonly("body_count", [](const World& w) { return w.bodies.size(); })
        .def_property_readonly("pair_count", &World::pair_count)
        .def("positions",
             [](const World& w) {
                 return columns(w, 2, [](const Body& b, double* o) {
                     o[0] = b.x.x;
                     o[1] = b.x.y;
                 });
             })
        .def("rotations",
             [](const World& w) {
                 return columns(w, 2, [](const Body& b, double* o) {
                     o[0] = b.q.c;
                     o[1] = b.q.s;
                 });
             })
        .def("velocities",
             [](const World& w) {
                 return columns(w, 3, [](const Body& b, double* o) {
                     o[0] = b.v.x;
                     o[1] = b.v.y;
                     o[2] = b.w;
                 });
             })
        .def(
            "set_velocity",
            [](World& w, int i, double vx, double vy, double av) {
                Body& b = w.bodies.at(i);
                b.v = {vx, vy};
                b.w = av;
            },
            py::arg("i"), py::arg("vx"), py::arg("vy"), py::arg("w") = 0.0)
        .def("mass", [](const World& w, int i) { return w.bodies.at(i).mass; })
        .def("inertia", [](const World& w, int i) { return w.bodies.at(i).inertia; })
        .def("shape",
             [](const World& w, int i) -> py::object {
                 const Shape& s = w.bodies.at(i).shape;
                 if (s.type == ShapeType::Circle) return py::make_tuple("circle", s.radius);
                 py::list verts;
                 for (int k = 0; k < s.count; ++k) {
                     verts.append(py::make_tuple(s.vertices[k].x, s.vertices[k].y));
                 }
                 return py::make_tuple("polygon", verts);
             })
        .def("contacts",
             [](const World& w) {
                 py::list out;
                 for (const Manifold& m : w.manifolds()) {
                     for (int k = 0; k < m.count; ++k) {
                         const ContactPoint& cp = m.points[k];
                         out.append(py::dict(py::arg("a") = m.a, py::arg("b") = m.b,
                                             py::arg("normal") = py::make_tuple(m.normal.x, m.normal.y),
                                             py::arg("point") = py::make_tuple(cp.point.x, cp.point.y),
                                             py::arg("separation") = cp.separation,
                                             py::arg("normal_impulse") = cp.normal_impulse));
                     }
                 }
                 return out;
             })
        .def("timings",
             [](const World& w) {
                 const StepTimings& t = w.timings;
                 return py::dict(py::arg("broadphase_s") = t.broadphase_s,
                                 py::arg("narrowphase_s") = t.narrowphase_s,
                                 py::arg("solve_s") = t.solve_s,
                                 py::arg("integrate_s") = t.integrate_s,
                                 py::arg("steps") = t.steps);
             })
        .def("reset_timings", [](World& w) { w.timings = StepTimings{}; });

    m.def("pile", &scenes::pile, "Reference pile scene shared with the C++ tests and benchmarks.",
          py::arg("n"), py::arg("broadphase") = Broadphase::SweepAndPrune);
    m.attr("PILE100_HASH600") = scenes::kPile100Hash600;
    m.def("stack",&scenes::stack, "A vertical stack of n unit boxes.", py::arg("n"));
}
