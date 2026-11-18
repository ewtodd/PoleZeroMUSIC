#!/usr/bin/env bash
# What run 70's waveforms say about the DPP-PHA baseline restorer with
# EnergyFilterLFLimitation = On: the firmware baseline at every traced hit,
# its change between consecutive hits on a channel, and CAEN's description
# fitted to that (scripts/run70_baseline.C says how). Writes
# doc/figs/run70_baseline.tex, LaTeX macros for the paper's text, and
# output/run70/baseline_change.png.
# The extraction of the long-channel hits from the 2.9 GB converted traces is
# kept in output/run70/run70_hits.root and skipped when that file exists.
# Run from anywhere: bash scripts/run70_baseline.sh
set -euo pipefail
cd "$(dirname "$0")/.."
IN=${RUN70:-/labdata/MUSIC/New-37Cl/data_raw/run70.root}
OUT=output/run70
mkdir -p "$OUT/aclic" doc/figs
ACLIC="gSystem->SetBuildDir(\"$OUT/aclic\", kTRUE);"
if [ ! -f "$OUT/run70_hits.root" ]; then
  nix develop -c root -l -b -q -e "$ACLIC" \
    "scripts/run70_baseline_extract.C+(\"$IN\", \"$OUT/run70_hits.root\")"
fi
nix develop -c root -l -b -q -e "$ACLIC" \
  "scripts/run70_baseline.C+(\"$OUT/run70_hits.root\", \"doc/figs/run70_baseline.tex\", \"$OUT/baseline_change.png\")" |
  tee "$OUT/run70_baseline.log"
