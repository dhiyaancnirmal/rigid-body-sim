#!/usr/bin/env bash
# Run the step and integrator benchmarks pinned to one core; write
# results/step_bench.json and results/integrators.json.
#
# The step benchmark also runs at a baseline commit (before the solver used
# packed SolverBody arrays) built from a temporary git worktree, so the
# before/after numbers come from the same machine and session.
#
#   scripts/run_benchmarks.sh [baseline-rev]   (default bde00f7)
set -euo pipefail
cd "$(dirname "$0")/.."
baseline="${1:-bde00f7}"
cpu="${BENCH_CPU:-3}"
py="${PYTHON:-python3}"

cmake -S . -B build -DRBSIM_PYTHON=OFF > /dev/null
cmake --build build -j2 --target bench_step bench_integrators > /dev/null

wt="build-det/baseline-src"
rm -rf "$wt"
git worktree add --detach "$wt" "$baseline" > /dev/null 2>&1
cmake -S "$wt" -B build-det/baseline -DRBSIM_PYTHON=OFF > /dev/null
cmake --build build-det/baseline -j2 --target bench_step > /dev/null
git worktree remove --force "$wt"

meta="$(mktemp)"
{
  echo "{"
  echo "  \"cpu\": \"$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')\","
  echo "  \"compiler\": \"$(c++ --version | head -1)\","
  echo "  \"flags\": \"-O3 (CMake Release) -ffp-contract=off -fno-fast-math\","
  echo "  \"pinned_cpu\": $cpu,"
  echo "  \"load_average_at_start\": \"$(cut -d' ' -f1-3 /proc/loadavg)\","
  echo "  \"baseline_rev\": \"$baseline\""
  echo "}"
} > "$meta"

rounds="${ROUNDS:-3}"
tmp="$(mktemp -d)"
# Alternate baseline and current so load changes on a shared machine hit
# both equally; the summary keeps each configuration's fastest round.
for r in $(seq 1 "$rounds"); do
  taskset -c "$cpu" build-det/baseline/bench_step > "$tmp/baseline_$r.json"
  taskset -c "$cpu" build/bench_step > "$tmp/current_$r.json"
done
taskset -c "$cpu" build/bench_integrators > results/integrators.json

"$py" - "$meta" "$tmp" "$rounds" <<'EOF'
import json, sys
meta = json.load(open(sys.argv[1]))
tmp, rounds = sys.argv[2], int(sys.argv[3])
meta["load_average_at_end"] = " ".join(open("/proc/loadavg").read().split()[:3])
meta["rounds"] = rounds
out = {"meta": meta}
for label in ("baseline", "current"):
    best = {}
    for r in range(1, rounds + 1):
        d = json.load(open(f"{tmp}/{label}_{r}.json"))
        for run in d["runs"]:
            key = (run["bodies"], run["broadphase"])
            if key not in best or run["step_ms"] < best[key]["step_ms"]:
                best[key] = run
    out[label] = {k: d[k] for k in ("warmup_steps", "window_steps", "repeats")}
    out[label]["runs"] = [best[k] for k in sorted(best)]
with open("results/step_bench.json", "w") as f:
    json.dump(out, f, indent=2)
    f.write("\n")
print("wrote results/step_bench.json and results/integrators.json")
EOF
rm -rf "$meta" "$tmp"
