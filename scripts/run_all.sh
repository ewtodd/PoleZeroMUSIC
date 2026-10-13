#!/usr/bin/env bash
# Production: build, simulate the runs the study asks for (40, 10 and 1 kHz,
# the tau_R = 20 us control and the sign check), then analyse.
# Run from the project root (nix build is run by the script):  bash scripts/run_all.sh
set -euo pipefail
cd "$(dirname "$0")/.."
# The results this script documents were produced with the earlier
# long-window defaults (262 us baseline, unconditional restart, width-derived
# noise). LEGACY reproduces them; the current defaults model the measured
# restorer instead (README, "Baseline restorer in production").
LEGACY="--baselineMeanSamples=32768 --baselineRestartThrADC=0 --lfFilterOn=off --noiseSigmaADC_L=960 --noiseSigmaADC_R=960"
N=${NEVENTS:-200000}
OUT=output
mkdir -p "$OUT"
nix build

run() {
  local name=$1; shift
  echo "simulate $name: $*"
  ./result/bin/simulate $LEGACY --nEvents="$N" --outputFile="$OUT/sim_$name.root" "$@" > "$OUT/sim_$name.log" 2>&1
}

run 40kHz --rate_kHz=40 &
run 10kHz --rate_kHz=10 &
run 1kHz --rate_kHz=1 &
# Control: both preamps matched to the 20 us pole-zero (the L chain is a 25 us
# part in the real detector, so tauL is overridden here too), a true null test.
run control_40kHz --rate_kHz=40 --tauL_us=20 --tauR_us=20 &
run signcheck_40kHz --rate_kHz=40 --tauR_us=20 --trapDecay_us=10 &
wait
for f in "$OUT"/sim_*.log; do echo "--- $f"; tail -3 "$f"; done

./result/bin/analyze \
  --main="$OUT/sim_40kHz.root,$OUT/sim_10kHz.root,$OUT/sim_1kHz.root" \
  --control="$OUT/sim_control_40kHz.root" \
  --signcheck="$OUT/sim_signcheck_40kHz.root" \
  --outdir="$OUT" 2>&1 | grep -v "^Info in <TCanvas::Print>" | tee "$OUT/analysis.log"
