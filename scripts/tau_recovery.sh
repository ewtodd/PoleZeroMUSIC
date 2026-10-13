#!/usr/bin/env bash
# Does the energy-vs-dt fit return the preamp decay time it was given?
# The long ends run a 20 us pole-zero and the short ends a 10 us one, so
# simulating both settings with the same two preamps produces all four strip
# groups and the extraction can be scored against a known answer. Errors are
# bootstrapped over the hits, with the spread over the start of the fit added
# as a systematic.
# Run from the project root (nix build is run by the script):  bash scripts/tau_recovery.sh
set -euo pipefail
cd "$(dirname "$0")/.."
# The results this script documents were produced with the earlier
# long-window defaults (262 us baseline, unconditional restart, width-derived
# noise). LEGACY reproduces them; the current defaults model the measured
# restorer instead (README, "Baseline restorer in production").
LEGACY="--baselineMeanSamples=32768 --baselineRestartThrADC=0 --lfFilterOn=off --noiseSigmaADC_L=960 --noiseSigmaADC_R=960"
OUT=output/tau_recovery
mkdir -p "$OUT"
nix build
./result/bin/simulate $LEGACY --nEvents=200000 --rate_kHz=27 --tauL_us=25 --tauR_us=10 \
  --trapDecay_us=20 --outputFile="$OUT/sim_long.root" > "$OUT/sim_long.log" 2>&1 &
./result/bin/simulate $LEGACY --nEvents=200000 --rate_kHz=27 --tauL_us=25 --tauR_us=10 \
  --trapDecay_us=10 --outputFile="$OUT/sim_short.root" > "$OUT/sim_short.log" 2>&1 &
wait
./result/bin/taustudy "$OUT/sim_long.root" "$OUT/sim_short.root" 2>&1 \
  | grep -v "^Info in" | tee "$OUT/tau_recovery.log"
# The same with the measured restorer (the current defaults).
./result/bin/simulate --nEvents=200000 --rate_kHz=27 --tauL_us=25 --tauR_us=10 \
  --trapDecay_us=20 --outputFile="$OUT/sim_long_meas.root" > "$OUT/sim_long_meas.log" 2>&1 &
./result/bin/simulate --nEvents=200000 --rate_kHz=27 --tauL_us=25 --tauR_us=10 \
  --trapDecay_us=10 --outputFile="$OUT/sim_short_meas.root" > "$OUT/sim_short_meas.log" 2>&1 &
wait
./result/bin/taustudy "$OUT/sim_long_meas.root" "$OUT/sim_short_meas.root" 2>&1 \
  | grep -v "^Info in" | tee "$OUT/tau_recovery_meas.log"
