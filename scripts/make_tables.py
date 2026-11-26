"""Regenerate the README tables from results/*.json.

    python scripts/make_tables.py          # rewrite README.md in place
    python scripts/make_tables.py --check  # exit 1 if README.md is stale

Each table lives between <!-- table:NAME --> and <!-- /table:NAME --> markers.
"""

import json
import math
import os
import re
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
README = os.path.join(ROOT, "README.md")


def load(name):
    with open(os.path.join(ROOT, "results", name)) as f:
        return json.load(f)


def sci(x):
    if x == 0:
        return "0"
    e = math.floor(math.log10(abs(x)))
    m = x / 10**e
    return f"{m:.1f}e{e}"


def pct(x, digits=2):
    return f"{100 * x:+.{digits}f}%"


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def step_bench():
    d = load("step_bench.json")
    base = {(r["bodies"], r["broadphase"]): r for r in d["baseline"]["runs"]}
    cur = {(r["bodies"], r["broadphase"]): r for r in d["current"]["runs"]}
    rows = []
    for n in sorted({k[0] for k in cur}):
        brute = base[(n, "brute_force")]["step_ms"]
        sap = base[(n, "sweep_and_prune")]["step_ms"]
        packed = cur[(n, "sweep_and_prune")]["step_ms"]
        same = len({base[(n, b)]["final_hash"] for b in ("brute_force", "sweep_and_prune")}
                   | {cur[(n, b)]["final_hash"] for b in ("brute_force", "sweep_and_prune")}) == 1
        rows.append([n, f"{brute:.2f}", f"{sap:.2f}", f"{packed:.2f}", f"{brute / packed:.1f}x",
                     "yes" if same else "NO"])
    return table(["bodies", "brute force (ms)", "+ sweep and prune (ms)",
                  "+ packed solver bodies (ms)", "total speedup", "same final hash"], rows)


def phases():
    d = load("step_bench.json")
    rows = []
    for label, runs in (("baseline", d["baseline"]["runs"]), ("current", d["current"]["runs"])):
        for r in runs:
            if r["bodies"] != 4000:
                continue
            name = {"brute_force": "brute force", "sweep_and_prune": "sweep and prune"}[r["broadphase"]]
            name += ", packed solver" if label == "current" else ", Body structs"
            rows.append([name, f"{r['broadphase_ms']:.2f}", f"{r['narrowphase_ms']:.2f}",
                         f"{r['solve_ms']:.2f}", f"{r['step_ms']:.2f}"])
    return table(["4000 bodies", "broadphase (ms)", "narrowphase (ms)", "solver (ms)",
                  "whole step (ms)"], rows)


def bench_meta():
    m = load("step_bench.json")["meta"]
    return (f"Measured on {m['cpu']}, {m['compiler']}, {m['flags']}, pinned to one core, "
            f"best of 5 repeats of 120 steps after 120 warm-up steps "
            f"(load average {m['load_average_at_start']} at start).")


METHOD = {"explicit_euler": "explicit Euler", "semi_implicit_euler": "semi-implicit Euler",
          "velocity_verlet": "velocity Verlet", "rk4": "RK4"}
SYSTEM = {"harmonic_oscillator": "harmonic oscillator", "pendulum_115deg": "pendulum, 115 deg",
          "kepler_e0.5": "Kepler orbit, e = 0.5"}


def integrators():
    d = load("integrators.json")
    rows = []
    for r in d["runs"]:
        if r["evals_per_period"] != 800:
            continue
        e = "diverged" if r["blew_up"] else sci(r["max_rel_energy_error"])
        p = "diverged" if r["blew_up"] else sci(r["final_position_error"])
        rows.append([SYSTEM[r["system"]], METHOD[r["method"]], r["steps_per_period"], e, p,
                     f"{r['ns_per_step']:.1f}"])
    return table(["system", "method", "steps per period", "max energy error",
                  "position error after 100 periods", "ns per step"], rows)


def long_horizon():
    d = load("integrators.json")["long_horizon"]
    rows = [[METHOD[r["method"]], sci(r["p100"]), sci(r["p1000"]), sci(r["p10000"])]
            for r in d["max_rel_energy_error_by_periods"]]
    return table(["method (200 force evaluations per period)", "100 periods", "1,000 periods",
                  "10,000 periods"], rows)


def pendulum():
    d = load("validation.json")["pendulum_period"]
    rows = [[f"{r['amplitude_deg']} deg", f"{r['exact_period_s']:.4f}", pct(r["rel_error_60hz"], 3),
             pct(r["rel_error_120hz"], 3), pct(r["rel_error_240hz"], 3)] for r in d]
    return table(["amplitude", "exact period (s)", "error at 60 Hz", "error at 120 Hz",
                  "error at 240 Hz"], rows)


def incline():
    d = load("validation.json")["incline_friction"]
    rows = []
    for r in d:
        expect = "stays put" if r["expected_m"] == 0 else f"{r['expected_m']:.4f} m"
        rows.append([f"{r['slope_deg']} deg", r["mu"], expect, f"{r['simulated_m']:.4f} m"])
    return table(["slope", "friction", "Coulomb model, 1 s", "simulated, 1 s"], rows)


def robustness():
    v = load("validation.json")
    t = v["tunneling"]
    disc = t["discrete_only"]["tunneled_at_m_s"]
    spec = t["speculative"]["tunneled_at_m_s"]
    top = max(t["speculative"]["tested_speeds_m_s"])
    gas = {r["circles"]: r["energy_ratio_after_20s"] for r in v["elastic_gas"]}
    stack = {r["velocity_iterations"]: r["max_stable_boxes"] for r in v["stack_stability"]}
    det = load("determinism.json")
    same = [b["name"] for b in det["builds"] if b["hash"] == det["expected"]]
    diff = [b["name"] for b in det["builds"] if b["hash"] != det["expected"]]
    mom = v["momentum"]

    def gas_range(n):
        lo, hi = min(gas[n]), max(gas[n])
        return f"{pct(lo - 1, 1)} to {pct(hi - 1, 1)}"

    rows = [
        ["Tunneling, 0.25 m circle into 0.1 m wall at 60 Hz, discrete contacts",
         f"passes through from {min(disc)} m/s" if disc else "never"],
        ["Same, speculative contacts", "never, up to %s m/s" % f"{top:,}" if not spec else
         f"passes at {min(spec)} m/s"],
        ["Linear momentum drift, contact + friction + rod scenes",
         f"{sci(mom['max_rel_linear_drift'])} (relative)"],
        ["Angular momentum drift, same scenes", f"{sci(mom['max_rel_angular_drift'])} (relative)"],
        ["Elastic gas energy change after 20 s, 2 and 10 circles",
         f"{gas_range(2)}; {gas_range(10)}"],
        ["Same, 50 and 100 circles (simultaneous contacts)", f"{gas_range(50)}; {gas_range(100)}"],
        ["Tallest stable stack (1 cm offsets), 8 / 16 / 32 iterations",
         f"{stack[8]} / {stack[16]} / {stack[32]} boxes"],
        ["Bit-identical pile hash across builds", ", ".join(same)],
        ["Hash changes (negative control)", ", ".join(diff) or "none"],
    ]
    return table(["check", "result"], rows)


TABLES = {
    "step_bench": step_bench,
    "phases": phases,
    "bench_meta": bench_meta,
    "integrators": integrators,
    "long_horizon": long_horizon,
    "pendulum": pendulum,
    "incline": incline,
    "robustness": robustness,
}


def render(text):
    for name, fn in TABLES.items():
        pattern = re.compile(rf"(<!-- table:{name} -->\n)(?:.*?\n)?(<!-- /table:{name} -->)", re.S)
        if not pattern.search(text):
            raise SystemExit(f"README.md has no markers for table '{name}'")
        text = pattern.sub(lambda m: m.group(1) + fn() + "\n" + m.group(2), text)
    return text


def main():
    with open(README) as f:
        current = f.read()
    updated = render(current)
    if "--check" in sys.argv:
        if updated != current:
            print("README.md tables are stale; run python scripts/make_tables.py")
            sys.exit(1)
        print("README.md tables match results/")
        return
    with open(README, "w") as f:
        f.write(updated)
    print("README.md tables regenerated")


if __name__ == "__main__":
    main()
