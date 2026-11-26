"""Fast bodies must not pass through thin walls (speculative contacts)."""

from hypothesis import given
from hypothesis import strategies as st

import rbsim

WALL_HALF = 0.05  # a wall 10 cm thick


def shoot(speed, shape, size, angle, dt, speculative=True):
    s = rbsim.Settings()
    s.gravity = (0.0, 0.0)
    s.speculative = speculative
    w = rbsim.World(s)
    w.add_box(0.0, 0.0, WALL_HALF, 20.0, static=True)
    if shape == "circle":
        w.add_circle(-3.0, 0.3, size, vx=speed)
    else:
        w.add_box(-3.0, 0.3, size, size, angle=angle, vx=speed)
    w.step(dt, n=int(3.0 / (speed * dt)) + 60)
    return w.positions()[1, 0]


@given(
    speed=st.floats(1.0, 10_000.0),
    shape=st.sampled_from(["circle", "box"]),
    size=st.floats(0.05, 0.5),
    angle=st.floats(-0.7, 0.7),
    dt=st.sampled_from([1 / 30, 1 / 60, 1 / 240]),
)
def test_no_tunneling_through_static_wall_up_to_10_km_per_s(speed, shape, size, angle, dt):
    x = shoot(speed, shape, size, angle, dt)
    # Centre stays on the near side, outside the wall (allowing the slop).
    assert x < -WALL_HALF


def test_discrete_contacts_alone_do_tunnel():
    # The negative control: same shot with speculative contacts off passes
    # through once a step moves the body further than it is thick.
    assert shoot(10.0, "circle", 0.25, 0.0, 1 / 60, speculative=False) < -WALL_HALF
    assert shoot(100.0, "circle", 0.25, 0.0, 1 / 60, speculative=False) > WALL_HALF
