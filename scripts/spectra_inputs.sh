#!/usr/bin/env bash
# Regenerate the simulated hits used by the pulse-height spectra. The left
# chain is matched to the 20 us pole-zero (tauL = 20 us); the right chain keeps
# tauR = 10 us. The baseline and noise settings reproduce the old spectra.
# Run from anywhere: bash scripts/spectra_inputs.sh
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=output/figures
N=${NEVENTS:-200000}
LEGACY="--baselineMeanSamples=32768 --baselineRestartThrADC=0 --lfFilterOn=off --noiseSigmaADC_L=960 --noiseSigmaADC_R=960"
mkdir -p "$OUT"
nix build
for rate in 1 40; do
  # shellcheck disable=SC2086
  ./result/bin/simulate $LEGACY --tauL_us=20 --nEvents="$N" --rate_kHz="$rate" \
    --outputFile="$OUT/long_${rate}kHz.root" > "$OUT/long_${rate}kHz.log" 2>&1 &
done
wait
