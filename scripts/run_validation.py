"""Measure the simulator against closed-form physics; write results/validation.json.

    PYTHONPATH=python python scripts/run_validation.py
"""

import json
import math
import os
import platform
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

import rbsim  # noqa: E402
from rbsim import experiments as ex  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "validation.json")


def pendulum():
    rows = []
    for amp in (10, 45, 90, 150):
        exact = ex.exact_pendulum_period(1.0, math.radians(amp))
        row = {"amplitude_deg": amp, "exact_period_s": exact}
        for hz in (60, 120, 240):
            sim = ex.simulated_pendulum_period(1.0, math.radians(amp), 1.0 / hz)
            row[f"rel_error_{hz}hz"] = (sim - exact) / exact
        rows.append(row)
    return rows


def incline():
    rows = []
    for slope, mu in ((20, 0.5), (25, 0.5), (35, 0.8), (30, 0.2), (45, 0.5), (60, 0.3)):
        holds = math.tan(math.radians(slope)) < mu
        sim = ex.slide_distance(slope, mu, 1.0)
        exact = 0.0 if holds else ex.coulomb_slide_distance(slope, mu, 1.0)
        rows.append({"slope_deg": slope, "mu": mu, "expected_m": exact, "simulated_m": sim})
    return rows


def tunneling():
    speeds = [5, 10, 15, 20, 22, 24, 26, 30, 50, 100, 1000, 10000]
    out = {"body": "circle r=0.25 m", "wall_thickness_m": 2 * ex.WALL_HALF, "dt_s": 1 / 60}
    for spec in (False, True):
        passed = [v for v in speeds if ex.shoot(v, "circle", 0.25, 0.0, 1 / 60, spec) > -ex.WALL_HALF]
        key = "speculative" if spec else "discrete_only"
        out[key] = {"tested_speeds_m_s": speeds, "tunneled_at_m_s": passed}
    return out


def elastic_gas():
    rows = []
    for n in (2, 10, 50, 100):
        ratios = [ex.elastic_gas_energy_ratio(n, seed) for seed in range(3)]
        rows.append({"circles": n, "energy_ratio_after_20s": ratios})
    return rows


STACK_HEIGHTS = (5, 10, 15, 20, 25, 30, 40)


def stacks():
    rows = []
    for iters in (8, 16, 32):
        best = 0
        for n in STACK_HEIGHTS:
            drift, sink = ex.stack_drift(n, iterations=iters)
            if drift < 0.01 and sink < 0.05:
                best = n
            else:
                break
        rows.append({"velocity_iterations": iters, "max_stable_boxes": best,
                     "tested_heights": list(STACK_HEIGHTS)})
    return rows


def momentum():
    # Same ring scene as the property test, at its worst-case size.
    worst_p = worst_l = 0.0
    for n in (4, 8, 12):
        for speed in (3.0, 10.0, 30.0):
            w = ex.free_space()
            for i in range(n):
                a = 2 * math.pi * i / n
                kw = dict(vx=-speed * math.cos(a), vy=-speed * math.sin(a), w=1.0, friction=0.5,
                          restitution=0.3)
                if i % 2:
                    w.add_circle(4 * math.cos(a), 4 * math.sin(a), 0.5, **kw)
                else:
                    w.add_box(4 * math.cos(a), 4 * math.sin(a), 0.5, 0.5, angle=a, **kw)
            p = w.positions()
            w.add_distance_joint(0, 1, tuple(p[0]), tuple(p[1]))  # a rigid rod
            p0, l0 = w.linear_momentum(), w.angular_momentum()
            scale = sum(w.mass(i) for i in range(n)) * speed * 5.0
            w.step(1 / 60, n=180)
            p1, l1 = w.linear_momentum(), w.angular_momentum()
            worst_p = max(worst_p, math.hypot(p1[0] - p0[0], p1[1] - p0[1]) / scale)
            worst_l = max(worst_l, abs(l1 - l0) / scale)
    return {"max_rel_linear_drift": worst_p, "max_rel_angular_drift": worst_l}


def main():
    results = {
        "python": platform.python_version(),
        "machine": platform.machine(),
        "pendulum_period": pendulum(),
        "incline_friction": incline(),
        "tunneling": tunneling(),
        "elastic_gas": elastic_gas(),
        "stack_stability": stacks(),
        "momentum": momentum(),
    }
    w = rbsim.pile(100)
    w.step(1 / 60, n=600)
    results["pile100_hash_after_600_steps"] = f"0x{w.state_hash():016x}"
    with open(OUT, "w") as f:
        json.dump(results, f, indent=2)
        f.write("\n")
    print(f"wrote {os.path.normpath(OUT)}")


if __name__ == "__main__":
    main()
