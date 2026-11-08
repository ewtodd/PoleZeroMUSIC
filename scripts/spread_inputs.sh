#!/usr/bin/env bash
# Regenerate the simulated runs behind doc/figs/spread_table.tex: two
# behaviours of the baseline restorer, each at 0.5 (isolated pulses), 10, 44
# and 80 kHz, with the preamplifier decay times of the spec (tauR = 10 us,
# tauL = 25 us) against the 20 us pole-zero. figures takes the decay times for
# the table's caption from these runs.
#   ideal - the true baseline is subtracted
#   wait  - conditional restart, the firmware's with EnergyFilterLFLimitation
#           on, as in the 37Cl runs: after the 7 us freeze, averaging resumes
#           only once the trapezoid has stayed within 80 ADC for 2 us
# Run from anywhere: bash scripts/spread_inputs.sh
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=output/spread
N=${NEVENTS:-300000}
mkdir -p "$OUT"
nix build

declare -A CASE=(
  [ideal]="--baselineIdeal=on"
  [wait]="--lfFilterOn=off --baselineRestartThrADC=80"
)
for reading in ideal wait; do
  for rate in 0.5 10 44 80; do
    n=$N
    # Without the smoothing a waiting baseline is rebuilt from 500 us before
    # every isolated pulse, which makes 0.5 kHz slow; 60k events suffice.
    if [ "$reading" = wait ] && [ "$rate" = 0.5 ]; then n=60000; fi
    # shellcheck disable=SC2086
    ./result/bin/simulate --nEvents="$n" --rate_kHz="$rate" --tauL_us=25 \
      --tauR_us=10 ${CASE[$reading]} --outputFile="$OUT/${reading}_r${rate}.root" \
      > "$OUT/${reading}_r${rate}.log" 2>&1 &
  done
done
wait
