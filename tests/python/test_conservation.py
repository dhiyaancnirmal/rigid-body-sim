"""Physical laws the solver must respect, checked on randomized scenes."""

import math

from hypothesis import given
from hypothesis import strategies as st

from rbsim.experiments import elastic_gas_energy_ratio, free_space

DT = 1.0 / 60.0

body = st.tuples(
    st.sampled_from(["circle", "box"]),
    st.floats(0.2, 0.8),  # size
    st.floats(0.5, 5.0),  # density
    st.floats(0.0, 1.0),  # friction
    st.floats(-2.0, 2.0),  # angular velocity
)


@given(bodies=st.lists(body, min_size=4, max_size=12), speed=st.floats(3.0, 30.0))
def test_contacts_and_rods_conserve_linear_and_angular_momentum(bodies, speed):
    # Bodies on a ring, all thrown at the centre, so they pile into each other
    # (at least two of them are free of the rod and meet there).
    # No gravity and no static bodies: nothing outside the system pushes on it.
    w = free_space()
    n = len(bodies)
    for i, (kind, size, density, friction, spin) in enumerate(bodies):
        a = 2.0 * math.pi * i / n
        x, y = 4.0 * math.cos(a), 4.0 * math.sin(a)
        kw = dict(vx=-speed * math.cos(a), vy=-speed * math.sin(a), w=spin,
                  density=density, friction=friction, restitution=0.3)
        if kind == "circle":
            w.add_circle(x, y, size, **kw)
        else:
            w.add_box(x, y, size, size, angle=a, **kw)
    # A rigid rod between the first two bodies' centres.
    p = w.positions()
    w.add_distance_joint(0, 1, tuple(p[0]), tuple(p[1]))

    p0, l0 = w.linear_momentum(), w.angular_momentum()
    scale = sum(w.mass(i) for i in range(n)) * speed * 5.0
    touched = False
    for _ in range(180):
        w.step(DT)
        touched = touched or bool(w.contacts())
    p1, l1 = w.linear_momentum(), w.angular_momentum()
    assert touched, "scene should produce contacts"
    assert abs(p1[0] - p0[0]) < 1e-12 * scale
    assert abs(p1[1] - p0[1]) < 1e-12 * scale
    assert abs(l1 - l0) < 1e-11 * scale


@given(
    r1=st.floats(0.1, 1.0), r2=st.floats(0.1, 1.0),
    d1=st.floats(0.2, 5.0), d2=st.floats(0.2, 5.0),
    v=st.floats(0.5, 40.0), offset=st.floats(-0.9, 0.9),
    w1=st.floats(-5.0, 5.0),
)
def test_elastic_circle_collision_conserves_energy(r1, r2, d1, d2, v, offset, w1):
    # Frictionless circles with restitution 1, including off-centre hits.
    w = free_space(restitution_threshold=0.0)
    w.add_circle(-3.0, 0.0, r1, vx=v, w=w1, density=d1, friction=0.0, restitution=1.0)
    w.add_circle(3.0, offset * (r1 + r2), r2, density=d2, friction=0.0, restitution=1.0)
    e0 = w.kinetic_energy()
    for _ in range(int(8.0 / (v * DT)) + 30):
        w.step(DT)
    v1 = w.velocities()
    assert v1[1, 0] > 0.0, "second circle should have been struck"
    assert abs(w.kinetic_energy() - e0) <= 1e-12 * e0


@given(v=st.floats(0.5, 200.0), angle=st.floats(-1.2, 1.2))
def test_elastic_bounce_off_static_wall_keeps_speed(v, angle):
    w = free_space(restitution_threshold=0.0)
    w.add_box(0.0, 0.0, 0.5, 50.0, static=True, friction=0.0, restitution=1.0)
    w.add_circle(-3.0, 0.0, 0.3, vx=v * math.cos(angle), vy=v * math.sin(angle),
                 friction=0.0, restitution=1.0)
    for _ in range(int(6.0 / (v * DT)) + 30):
        w.step(DT)
    vx, vy, _ = w.velocities()[1]
    assert vx < 0.0
    assert abs(math.hypot(vx, vy) - v) <= 1e-12 * v


def test_dense_elastic_gas_keeps_energy():
    # 100 elastic circles in an elastic box for 20 s, with many simultaneous
    # contacts. Regression: speculative contacts that switched on mid-solve
    # used to stop circles dead at the surface and lost 22-31% of the energy.
    for seed in range(3):
        assert abs(elastic_gas_energy_ratio(100, seed) - 1.0) < 1e-3
