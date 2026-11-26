"""Comparison against closed-form physics."""

import math

import pytest

from rbsim.experiments import (
    coulomb_slide_distance,
    exact_pendulum_period,
    simulated_pendulum_period,
    slide_distance,
)


@pytest.mark.parametrize("amplitude_deg", [10, 45, 90, 150])
def test_pendulum_period_matches_elliptic_integral(amplitude_deg):
    amp = math.radians(amplitude_deg)
    exact = exact_pendulum_period(1.0, amp)
    sim = simulated_pendulum_period(1.0, amp, 1.0 / 240.0)
    assert abs(sim - exact) / exact < 0.005


@pytest.mark.parametrize("slope_deg,mu", [(20, 0.5), (25, 0.5), (35, 0.8)])
def test_box_holds_on_incline_below_friction_angle(slope_deg, mu):
    assert math.tan(math.radians(slope_deg)) < mu
    assert abs(slide_distance(slope_deg, mu, 2.0)) < 1e-3


@pytest.mark.parametrize("slope_deg,mu", [(30, 0.2), (45, 0.5), (60, 0.3)])
def test_box_slides_with_coulomb_acceleration(slope_deg, mu):
    expected = coulomb_slide_distance(slope_deg, mu, 1.0)
    assert abs(slide_distance(slope_deg, mu, 1.0) - expected) / expected < 0.02
