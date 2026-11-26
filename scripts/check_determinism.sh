#!/usr/bin/env bash
# Build the C++ core under several optimization levels and instruction sets
# and record the pile(100) hash each produces; write results/determinism.json.
# The last configuration is a negative control: it allows FMA contraction,
# which changes rounding, so its hash is expected to differ.
set -euo pipefail
cd "$(dirname "$0")/.."

configs=(
  "O0|-O0|-ffp-contract=off;-fno-fast-math"
  "O2|-O2|-ffp-contract=off;-fno-fast-math"
  "O3|-O3|-ffp-contract=off;-fno-fast-math"
  "O3_avx2_fma|-O3 -march=x86-64-v3|-ffp-contract=off;-fno-fast-math"
  "O3_avx2_fma_contract_fast|-O3 -march=x86-64-v3|-ffp-contract=fast"
)

out=results/determinism.json
compiler="$(${CXX:-c++} --version | head -1)"
{
  echo "{"
  echo "  \"compiler\": \"${compiler}\","
  echo "  \"expected\": \"$(grep -o '0x[0-9a-f]*ull' include/rbsim/scenes.hpp | head -1 | sed 's/ull//')\","
  echo "  \"builds\": ["
} > "$out"

first=1
for cfg in "${configs[@]}"; do
  IFS='|' read -r name cxxflags fpflags <<< "$cfg"
  dir="build-det/$name"
  cmake -S . -B "$dir" -DRBSIM_PYTHON=OFF -DCMAKE_BUILD_TYPE=None \
        -DCMAKE_CXX_FLAGS="$cxxflags" -DRBSIM_FP_FLAGS="$fpflags" > /dev/null
  cmake --build "$dir" -j2 --target rbsim_tests > /dev/null
  hash="$({ "$dir/rbsim_tests" || true; } | grep -o 'hash after 600 steps: 0x[0-9a-f]*' | awk '{print $NF}')"
  echo "$name $hash"
  [ $first -eq 1 ] || echo "," >> "$out"
  first=0
  printf '    {"name": "%s", "cxx_flags": "%s", "fp_flags": "%s", "hash": "%s"}' \
    "$name" "$cxxflags" "$fpflags" "$hash" >> "$out"
done
printf '\n  ]\n}\n' >> "$out"
echo "wrote $out"
