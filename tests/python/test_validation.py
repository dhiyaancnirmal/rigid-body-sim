"""Comparison against closed-form physics."""

import math

import pytest

import rbsim

G = 9.81


def exact_pendulum_period(length, amplitude):
    # T = 4 K(k) / w0 with k = sin(amplitude / 2); K by the arithmetic-geometric mean.
    k = math.sin(amplitude / 2.0)
    a, b = 1.0, math.sqrt(1.0 - k * k)
    while abs(a - b) > 1e-15:
        a, b = 0.5 * (a + b), math.sqrt(a * b)
    return 4.0 * (math.pi / (2.0 * a)) / math.sqrt(G / length)


def simulated_pendulum_period(length, amplitude, dt, periods=4):
    w = rbsim.World()
    pivot = w.add_circle(0.0, 0.0, 0.01, static=True)
    x, y = length * math.sin(amplitude), -length * math.cos(amplitude)
    bob = w.add_circle(x, y, 0.05)
    w.add_distance_joint(pivot, bob, (0.0, 0.0), (x, y))
    crossings = []
    prev_x, t = x, 0.0
    t_end = (periods + 0.6) * exact_pendulum_period(length, amplitude)
    while t < t_end:
        w.step(dt)
        t += dt
        cur_x = w.positions()[bob, 0]
        if (prev_x > 0.0) != (cur_x > 0.0):
            crossings.append(t - dt * cur_x / (cur_x - prev_x))  # linear interpolation
        prev_x = cur_x
    halves = [b - a for a, b in zip(crossings, crossings[1:])]
    return 2.0 * sum(halves) / len(halves)


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
    th = math.radians(slope_deg)
    t = 1.0
    expected = 0.5 * G * (math.sin(th) - mu * math.cos(th)) * t * t
    assert abs(slide_distance(slope_deg, mu, t) - expected) / expected < 0.02


def slide_distance(slope_deg, mu, seconds, dt=1.0 / 240.0):
    th = math.radians(slope_deg)
    c, s = math.cos(th), math.sin(th)
    w = rbsim.World()
    w.add_box(0.0, 0.0, 20.0, 0.5, angle=th, static=True, friction=mu)
    # Box resting on the slope surface, 0.75 m along the slope normal.
    x0, y0 = -s * 0.75, c * 0.75
    w.add_box(x0, y0, 0.25, 0.25, angle=th, friction=mu)
    w.step(dt, n=round(seconds / dt))
    x, y = w.positions()[1]
    # Distance travelled down the slope.
    return -((x - x0) * c + (y - y0) * s)
