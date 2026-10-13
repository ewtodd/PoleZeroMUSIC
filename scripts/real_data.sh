#!/usr/bin/env bash
# Real-data checks on the 37Cl SOL files (see README, "What the real data say"):
#   realwave  run-70 waveforms: preamp decay, trapezoid tails, firmware baseline;
#   realdt    energy-vs-dt profiles per strip group through the hit streams.
# Run from the project root (nix build is run by the script):  bash scripts/real_data.sh
set -euo pipefail
cd "$(dirname "$0")/.."
RAW=/labdata/MUSIC/New-37Cl/data_raw
nix build
mkdir -p output/real
./result/bin/realwave --out=output/real/waveforms/run70.root "$RAW"/music_exp1915_070_00_66222_0*.sol
for r in 072 074 080 090 091 097; do
  f=$(ls "$RAW"/music_exp1915_${r}_*_66222_000.sol | head -1)
  ./result/bin/realdt --label=run$r --out=output/real/profiles --maxhits=20000000 "$f"
done
