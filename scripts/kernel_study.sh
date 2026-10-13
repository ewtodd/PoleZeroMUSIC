#!/usr/bin/env bash
# What limits the pulse-history correction, and what would improve it.
# Fits every variant on the first half of a run and scores it on the second,
# so parameter counts are compared honestly. See README, "What the simulation
# says about the correction".
# Run from the project root (nix build is run by the script):  bash scripts/kernel_study.sh
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=output/kernel_study
mkdir -p "$OUT"
nix build
./result/bin/simulate --nEvents=200000 --rate_kHz=40 --outputFile="$OUT/sim_40kHz.root" > "$OUT/sim_40kHz.log" 2>&1
./result/bin/simulate --nEvents=8000 --rate_kHz=40 --seed=777 --outputFile="$OUT/sim_40kHz_small.root" > "$OUT/sim_small.log" 2>&1
./result/bin/kernelstudy "$OUT/sim_40kHz.root" "$OUT/sim_40kHz_small.root" 2>&1 \
  | grep -v "^Info in" | tee "$OUT/kernel_study.log"
