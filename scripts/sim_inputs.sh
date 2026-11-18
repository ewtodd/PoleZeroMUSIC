#!/usr/bin/env bash
# Regenerate the simulated runs behind the paper's simulation section: the
# spectra and kernel figures and the spread table. Every run uses the DPP-PHA
# settings of the paper's Table 1 (the constants of the simulator's Config), a
# 10 % Gaussian spread of the beam amplitudes, and the same per-sample noise on
# both chains (563 ADC, the right chain's run-70 traces), so the two chains
# differ only in their decay time: the left is matched to the 20 us pole-zero
# (tauL = 20 us), the right mismatched (tauR = 10 us).
# The baseline restorer is CAEN's description of EnergyFilterLFLimitation = On
# (support, 2026-10-19): after each trigger the baseline is frozen for 7 us;
# at the end of each freeze it resumes if the trapezoid is back close to the
# baseline, and otherwise freezes again for another 7 us. CAEN does not say how
# close or what decides the check. Both readings use 20 ADC, about 1 % of the
# beam pulse: the run-70 waveforms (scripts/run70_baseline.sh) show the
# baseline does not move between pulses on channels with an overshoot, which
# CAEN's description fits only with a threshold of 20 or less. Two readings of
# what decides the check, each at 10 and 50 kHz:
#   average - the running average of the trapezoid over the baseline window
#             (1024 samples), kept up during the freeze, against the held
#             baseline; noise barely moves it, so the baseline resumes only
#             once a tail has gone
#   single  - the last trapezoid sample alone, which the trapezoid noise can
#             carry into the band while a tail is still there
# Run from anywhere: bash scripts/sim_inputs.sh
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=output/sim
N=${NEVENTS:-200000}
COMMON="--tauL_us=20 --tauR_us=10 --noiseSigmaADC_L=563 --noiseSigmaADC_R=563
  --beamPeakSigmaFrac=0.10 --lfFilterOn=off --baselineRestartThrADC=20
  --baselineRefreeze=on"
declare -A CASE=(
  [average]="--baselineCheckSamples=1024"
  [single]="--baselineCheckSamples=1"
)
mkdir -p "$OUT"
nix build
for reading in average single; do
  for rate in 10 50; do
    # shellcheck disable=SC2086
    ./result/bin/simulate $COMMON ${CASE[$reading]} --nEvents="$N" \
      --rate_kHz="$rate" --outputFile="$OUT/${reading}_r${rate}.root" \
      > "$OUT/${reading}_r${rate}.log" 2>&1 &
  done
done
wait
