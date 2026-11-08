#!/usr/bin/env bash
# Simulated runs with reacted particles, to show the correction on them too.
# Their long-end deposits come from MUSIC's 37Cl simulation, (alpha,n) and
# (alpha,p) weighted by TALYS (scripts/reaction_pairs.C). The reaction rate is
# raised to 1 % for statistics; the physical rate is about 2e-5 per beam
# particle. Settings and seeds are those of scripts/spread_inputs.sh, so the
# beam is event for event that of the reaction-free runs in output/spread.
# Run from anywhere: bash scripts/reaction_inputs.sh
set -euo pipefail
cd "$(dirname "$0")/.."
MUSIC=${MUSIC:-$HOME/MUSIC}
SIMDIR=${SIMDIR:-$MUSIC/analysis/37Cl/sim_root_files}
OUT=output/reactions
FRACTION=${FRACTION:-0.01}
N=${NEVENTS:-300000}
mkdir -p "$OUT"
nix build
nix develop -c root -l -b -q "scripts/reaction_pairs.C(\"$SIMDIR\", \"$OUT/pairs.txt\")"

declare -A CASE=(
  [ideal]="--baselineIdeal=on"
  [wait]="--lfFilterOn=off --baselineRestartThrADC=80"
)
for run in ideal_r44 ideal_r80 wait_r80; do
  reading=${run%_r*}
  rate=${run#*_r}
  # shellcheck disable=SC2086
  ./result/bin/simulate --nEvents="$N" --rate_kHz="$rate" --tauL_us=25 \
    --tauR_us=10 ${CASE[$reading]} --reactionFraction="$FRACTION" \
    --reactionPairs="$OUT/pairs.txt" --outputFile="$OUT/$run.root" \
    > "$OUT/$run.log" 2>&1 &
done
wait
