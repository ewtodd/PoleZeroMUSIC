#!/usr/bin/env bash
# Extra studies beyond the original task, each analysed on its own kernel:
#   linearity   40 kHz with a 25 % beam spread, for the excess-vs-previous-
#               amplitude linearity check with a wide lever arm;
#   medium1024  the literal reading of "Medium" (1024 samples at 8 ns);
#   high16384   the "High" restorer.
# Run from the project root (nix build is run by the script):  bash scripts/extras.sh
set -euo pipefail
cd "$(dirname "$0")/.."
# The results this script documents were produced with the earlier
# long-window defaults (262 us baseline, unconditional restart, width-derived
# noise). LEGACY reproduces them; the current defaults model the measured
# restorer instead (README, "Baseline restorer in production").
LEGACY="--baselineMeanSamples=32768 --baselineRestartThrADC=0 --lfFilterOn=off --noiseSigmaADC_L=960 --noiseSigmaADC_R=960"
N=${NEVENTS:-200000}
OUT=output/extras
mkdir -p "$OUT"
nix build
./result/bin/simulate $LEGACY --nEvents="$N" --rate_kHz=40 --beamPeakSigmaFrac=0.25 --outputFile="$OUT/sim_linearity_40kHz.root" > "$OUT/sim_linearity.log" 2>&1 &
./result/bin/simulate $LEGACY --nEvents="$N" --rate_kHz=40 --baselineMeanSamples=1024 --outputFile="$OUT/sim_medium1024_40kHz.root" > "$OUT/sim_medium1024.log" 2>&1 &
./result/bin/simulate $LEGACY --nEvents="$N" --rate_kHz=40 --baselineMeanSamples=16384 --outputFile="$OUT/sim_high16384_40kHz.root" > "$OUT/sim_high16384.log" 2>&1 &
wait
./result/bin/analyze \
  --extra="linearity:$OUT/sim_linearity_40kHz.root,medium1024:$OUT/sim_medium1024_40kHz.root,high16384:$OUT/sim_high16384_40kHz.root" \
  --outdir="$OUT" 2>&1 | grep -v "^Info in <TCanvas::Print>" | tee "$OUT/analysis.log"
