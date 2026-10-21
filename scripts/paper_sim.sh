#!/usr/bin/env bash
# Simulation study for doc/kernel_note.tex: does the fitted pulse-history
# kernel work whatever the firmware's baseline does with the pole-zero tail?
# Three readings of the baseline, each at 0.5 (isolated pulses), 10, 44 and
# 80 kHz, with the run-97 R decay time (11 us); then the note's figures.
#   ideal  - the true baseline is subtracted (never takes in the tail)
#   wait   - frozen until the trapezoid stays within 80 ADC for 2 us
#   follow - averaging restarts as soon as the 7 us freeze ends
#   caen   - one single-sample check per trapezoid length at 20 ADC (the
#            reading that keeps the baseline stuck; ruled out by run 97)
# Run from the project root:  bash scripts/paper_sim.sh
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=output/paper_sim
N=${NEVENTS:-300000}
mkdir -p "$OUT"
nix build

declare -A CASE=(
  [ideal]="--baselineIdeal=on"
  [wait]="--lfFilterOn=off --baselineRestartThrADC=80"
  [follow]="--lfFilterOn=off --baselineRestartThrADC=0"
  [caen]="--lfFilterOn=off --baselineRefreeze=on --baselineRestartThrADC=20"
)
for c in ideal wait follow caen; do
  for r in 0.5 10 44 80; do
    n=$N
    # Without the smoothing a waiting baseline is rebuilt from 500 us before
    # every isolated pulse, which makes 0.5 kHz slow; 60k events suffice.
    if [ "$c" = wait ] && [ "$r" = 0.5 ]; then n=60000; fi
    # shellcheck disable=SC2086
    ./result/bin/simulate --nEvents="$n" --rate_kHz="$r" --tauR_us=11 ${CASE[$c]} \
      --outputFile="$OUT/${c}_r${r}.root" > "$OUT/${c}_r${r}.log" 2>&1 &
  done
done
wait

EXTRA=""
for c in ideal wait follow caen; do
  for r in 0.5 10 44 80; do EXTRA="$EXTRA,${c}_r${r}:$OUT/${c}_r${r}.root"; done
done
./result/bin/analyze --extra="${EXTRA#,}" --outdir="$OUT/analysis" 2>&1 \
  | grep -v "^Info in" > "$OUT/analysis.log"

for format in pdf png; do
  out=doc/figs
  [ "$format" = png ] && out=doc/figs/png
  ./result/bin/notefigures --log="$HOME/MUSIC/analysis/37Cl/pipeline_fused.log" \
    --sim="$OUT/analysis.log" --out="$out" --format=$format
done
