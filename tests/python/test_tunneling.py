"""Fast bodies must not pass through thin walls (speculative contacts)."""

from hypothesis import given
from hypothesis import strategies as st

from rbsim.experiments import WALL_HALF, shoot


@given(
    speed=st.floats(1.0, 10_000.0),
    shape=st.sampled_from(["circle", "box"]),
    size=st.floats(0.05, 0.5),
    angle=st.floats(-0.7, 0.7),
    dt=st.sampled_from([1 / 30, 1 / 60, 1 / 240]),
    restitution=st.sampled_from([0.0, 0.5, 1.0]),
)
def test_no_tunneling_through_static_wall_up_to_10_km_per_s(speed, shape, size, angle, dt,
                                                            restitution):
    # Centre stays on the near side, outside the 10 cm wall.
    assert shoot(speed, shape, size, angle, dt, restitution=restitution) < -WALL_HALF


def test_discrete_contacts_alone_do_tunnel():
    # Negative control: with speculative contacts off, the same shot passes
    # through once a step moves the body further than body plus wall.
    assert shoot(10.0, "circle", 0.25, 0.0, 1 / 60, speculative=False) < -WALL_HALF
    assert shoot(100.0, "circle", 0.25, 0.0, 1 / 60, speculative=False) > WALL_HALF
