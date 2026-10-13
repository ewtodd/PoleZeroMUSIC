#!/usr/bin/env bash
# Grid of rate x pole-zero mismatch for the applicability figure: a 10 us
# preamp on the 37Cl trapezoid (2 us rise, 3 us flat), the pole-zero set
# matched and mismatched both ways, each run corrected with its own kernel.
# Writes output/applicability/grid.txt (label rate_kHz tauPz_us rmsRaw rmsCorr
# rmsFloor) for mains/applicability.cpp.
# Run from the project root (nix build is run by the script):  bash scripts/applicability.sh
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=output/applicability
mkdir -p "$OUT"
nix build
N=${NEVENTS:-60000}
RATES="5 10 20 40 80"
PZS="8 10 12 15 20 30"
EXTRA=""
for r in $RATES; do
  for pz in $PZS; do
    label="r${r}_pz${pz}"
    ./result/bin/simulate --nEvents="$N" --rate_kHz="$r" --tauL_us=10 --tauR_us=10 \
      --trapDecay_us="$pz" --outputFile="$OUT/$label.root" > "$OUT/$label.log" 2>&1 &
    EXTRA="$EXTRA,$label:$OUT/$label.root"
  done
  wait
done
./result/bin/analyze --extra="${EXTRA#,}" --outdir="$OUT/analysis" 2>&1 \
  | grep -v "^Info in" > "$OUT/analysis.log"
# Summary rows look like: label  kHz  ch  N  mean  RMSraw  RMScorr  kernel ...
awk '$1 ~ /^r[0-9]+_pz[0-9]+$/ && $3 == "L" {
       split($1, a, "_"); r = substr(a[1], 2); pz = substr(a[2], 3);
       print $1, r, pz, $6, $7 }' "$OUT/analysis.log" | sort -k2,2n -k3,3n > "$OUT/grid.txt"
wc -l "$OUT/grid.txt"
cat "$OUT/grid.txt"
