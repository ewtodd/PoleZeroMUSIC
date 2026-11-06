#!/usr/bin/env bash
# The run-97 pulse-height-vs-gap figure of doc/kernel_note.tex, before and
# after the pulse-history correction, as PDF (doc/figs) and PNG (doc/figs/png).
# Runs PulseHistory with the 37Cl configuration on the first run-97 file
# (about 15 GB of memory), built against the local ~/MUSIC tooling, as
# production is, and keeps its profiles in output/run97/deviation.root for
# figures to draw.
# Run from the project root:  bash scripts/run97_figure.sh
set -euo pipefail
cd "$(dirname "$0")/.."
SOL=${SOL:-/labdata/MUSIC/New-37Cl/data_raw/music_exp1915_097_00_66222_000.sol}
MUSIC=${MUSIC:-$HOME/MUSIC}
CONFIG=${CONFIG:-$MUSIC/analysis/37Cl/config/Constants.cpp}
OUT=output/run97
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$OUT"
nix build

export MUSIC CONFIG SOL WORK OUT
nix develop -c bash -c '
  g++ -O2 -std=c++17 $(root-config --cflags) -I"$MUSIC/tooling/include" \
    scripts/run97_pulse_history.cpp "$CONFIG" \
    "$MUSIC/tooling/src/PulseHistory.cpp" "$MUSIC/tooling/src/Constants.cpp" \
    $(root-config --glibs) -lMinuit -lanalysis-utils -o "$WORK/run97_pulse_history"
  "$WORK/run97_pulse_history" "$OUT/deviation.root" run97_c000 "$SOL" > "$WORK/run.log"
  for format in pdf png; do
    out=doc/figs
    [ "$format" = png ] && out=doc/figs/png
    mkdir -p "$out"
    ./result/bin/figures --run97="$OUT/deviation.root" --out="$out" \
      --format=$format 2>&1 | { grep -v "^Info in" || true; }
  done'
grep -E "long R|long L|short L" "$WORK/run.log" | head -3
