"""Measurement scenes shared by the tests and scripts/run_validation.py.

Each function builds a scene, runs it and returns a measured number. The
reference values they are compared against (exact pendulum period, Coulomb
sliding distance) come from closed-form physics, not from the simulator.
"""

import math

from ._rbsim import Settings, World

G = 9.81
WALL_HALF = 0.05  # thin wall used by the tunneling shots: 10 cm thick


def free_space(**overrides):
    """A world with no gravity, with optional Settings overrides."""
    s = Settings()
    s.gravity = (0.0, 0.0)
    for k, v in overrides.items():
        setattr(s, k, v)
    return World(s)


def exact_pendulum_period(length, amplitude):
    """T = 4 K(k) / w0 with k = sin(amplitude / 2); K by the arithmetic-geometric mean."""
    k = math.sin(amplitude / 2.0)
    a, b = 1.0, math.sqrt(1.0 - k * k)
    while abs(a - b) > 1e-15:
        a, b = 0.5 * (a + b), math.sqrt(a * b)
    return 4.0 * (math.pi / (2.0 * a)) / math.sqrt(G / length)


def simulated_pendulum_period(length, amplitude, dt, periods=4):
    """Period of a bob on a rigid rod, from interpolated zero crossings."""
    w = World()
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
            crossings.append(t - dt * cur_x / (cur_x - prev_x))
        prev_x = cur_x
    halves = [b - a for a, b in zip(crossings, crossings[1:])]
    return 2.0 * sum(halves) / len(halves)


def slide_distance(slope_deg, mu, seconds, dt=1.0 / 240.0):
    """Distance a box resting on an incline travels down it in `seconds`."""
    th = math.radians(slope_deg)
    c, s = math.cos(th), math.sin(th)
    w = World()
    w.add_box(0.0, 0.0, 20.0, 0.5, angle=th, static=True, friction=mu)
    x0, y0 = -s * 0.75, c * 0.75  # resting on the slope surface
    w.add_box(x0, y0, 0.25, 0.25, angle=th, friction=mu)
    w.step(dt, n=round(seconds / dt))
    x, y = w.positions()[1]
    return -((x - x0) * c + (y - y0) * s)


def coulomb_slide_distance(slope_deg, mu, seconds):
    th = math.radians(slope_deg)
    return 0.5 * G * (math.sin(th) - mu * math.cos(th)) * seconds * seconds


def shoot(speed, shape, size, angle, dt, speculative=True, restitution=0.0):
    """Fire a body at a thin static wall; return its final x (wall spans |x| < 0.05)."""
    w = free_space(speculative=speculative)
    w.add_box(0.0, 0.0, WALL_HALF, 20.0, static=True)
    if shape == "circle":
        w.add_circle(-3.0, 0.3, size, vx=speed, restitution=restitution)
    else:
        w.add_box(-3.0, 0.3, size, size, angle=angle, vx=speed, restitution=restitution)
    w.step(dt, n=int(3.0 / (speed * dt)) + 60)
    return w.positions()[1, 0]


def elastic_gas_energy_ratio(n, seed, seconds=20.0, dt=1.0 / 60.0):
    """Final/initial kinetic energy of n elastic circles in an elastic box."""
    import random

    rng = random.Random(seed)
    w = free_space(restitution_threshold=0.0)
    half = 10.0
    for x, y, hw, hh in [(0, -half - 0.5, half + 1, 0.5), (0, half + 0.5, half + 1, 0.5),
                         (-half - 0.5, 0, 0.5, half + 1), (half + 0.5, 0, 0.5, half + 1)]:
        w.add_box(x, y, hw, hh, static=True, friction=0.0, restitution=1.0)
    k = math.ceil(math.sqrt(n))
    gap = (2 * half - 2) / max(k - 1, 1)
    for i in range(n):
        r, c = divmod(i, k)
        w.add_circle(-half + 1 + c * gap, -half + 1 + r * gap, 0.4,
                     vx=rng.uniform(-5, 5), vy=rng.uniform(-5, 5), friction=0.0, restitution=1.0)
    e0 = w.kinetic_energy()
    w.step(dt, n=round(seconds / dt))
    return w.kinetic_energy() / e0


def stack_drift(n, seconds=10.0, dt=1.0 / 60.0, iterations=8):
    """Max sideways drift and top-box sink of an n-box stack after `seconds`."""
    w = World()
    w.settings.velocity_iterations = iterations
    w.add_box(0.0, -0.5, 20.0, 0.5, static=True)
    # Each box offset 1 cm alternately so the stack is not perfectly symmetric.
    for i in range(n):
        w.add_box(0.01 * (1 if i % 2 else -1), 0.5 + i, 0.5, 0.5, friction=0.6)
    w.step(dt, n=round(seconds / dt))
    p = w.positions()
    drift = max(abs(p[i, 0] - 0.01 * (1 if (i - 1) % 2 else -1)) for i in range(1, n + 1))
    sink = (n - 0.5) - p[n, 1]
    return drift, sink
