#!/usr/bin/env bash
# Scan tau_R with a non-tracking baseline (High = 16384 samples) and compare
# the PulseHistory kernel PulseHistory measures on the simulation with the
# real 37Cl long-R kernel from pipeline_fused.log (run 97).
set -euo pipefail
cd "$(dirname "$0")/.."
# The results this script documents were produced with the earlier
# long-window defaults (262 us baseline, unconditional restart, width-derived
# noise). LEGACY reproduces them; the current defaults model the measured
# restorer instead (README, "Baseline restorer in production").
LEGACY="--baselineMeanSamples=32768 --baselineRestartThrADC=0 --lfFilterOn=off --noiseSigmaADC_L=960 --noiseSigmaADC_R=960"
N=${NEVENTS:-100000}
OUT=output/tau_scan
mkdir -p "$OUT"
EXTRA=""
for tau in 8 10 12 14 16; do
  ./result/bin/simulate $LEGACY --nEvents="$N" --rate_kHz=40 --tauR_us="$tau" --baselineMeanSamples=16384 \
    --outputFile="$OUT/sim_tauR${tau}.root" > "$OUT/sim_tauR${tau}.log" 2>&1 &
  EXTRA="$EXTRA,tauR${tau}:$OUT/sim_tauR${tau}.root"
done
wait
./result/bin/analyze --extra="${EXTRA#,}" --outdir="$OUT" 2>&1 | grep -v "^Info in <TCanvas::Print>" > "$OUT/analysis.log"
echo "real 37Cl long R (run97):   1us:-0.049 2us:-0.118 3us:-0.142 4us:-0.137 6us:-0.119 8us:-0.099 12us:-0.074 18us:-0.046 26us:-0.021 38us:-0.004 56us:+0.003 83us:+0.004"
grep -A1 "^    long R" "$OUT/analysis.log" | grep "kernel:" | paste <(grep "^=== " "$OUT/analysis.log" | awk '{print $2}') -
