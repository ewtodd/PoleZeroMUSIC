#!/usr/bin/env bash
# The run-97 energy-vs-gap figure of doc/kernel_note.tex, before and after the
# pulse-history correction, as PDF (doc/figs) and PNG (doc/figs/png). Runs
# PulseHistory with the 37Cl configuration on the first run-97 file (about
# 15 GB of memory), built against the local ~/MUSIC tooling, as production is.
# Run from the project root:  bash scripts/run97_figure.sh
set -euo pipefail
cd "$(dirname "$0")/.."
SOL=${SOL:-/labdata/MUSIC/New-37Cl/data_raw/music_exp1915_097_00_66222_000.sol}
MUSIC=${MUSIC:-$HOME/MUSIC}
CONFIG=${CONFIG:-$MUSIC/analysis/37Cl/config/Constants.cpp}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

export MUSIC CONFIG SOL WORK
nix develop -c bash -c '
  g++ -O2 -std=c++17 $(root-config --cflags) -I"$MUSIC/tooling/include" \
    scripts/run97_pulse_history.cpp "$CONFIG" \
    "$MUSIC/tooling/src/PulseHistory.cpp" "$MUSIC/tooling/src/Constants.cpp" \
    $(root-config --glibs) -lMinuit -lanalysis-utils -o "$WORK/run97_pulse_history"
  for fmt in pdf png; do
    "$WORK/run97_pulse_history" $fmt "$WORK/$fmt" run97_c000 "$SOL" > "$WORK/$fmt.log"
  done'

src=$(dirname "$(find "$WORK/pdf" -name 'deviation_vs_dt_previous_R.pdf' | head -1)")
cp "$src/deviation_vs_dt_previous_R.pdf" doc/figs/note_run97_deviation_R.pdf
src=$(dirname "$(find "$WORK/png" -name 'deviation_vs_dt_previous_R.png' | head -1)")
mkdir -p doc/figs/png
cp "$src/deviation_vs_dt_previous_R.png" doc/figs/png/note_run97_deviation_R.png
grep -E "long R|long L|short L" "$WORK/pdf.log" | head -3
