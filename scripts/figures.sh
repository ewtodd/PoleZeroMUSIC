#!/usr/bin/env bash
# Figures for doc/polezero.tex. Simulates the runs the note quotes, then draws
# the five plots into doc/figs as PDF.
# Run from the project root (nix build is run by the script):  bash scripts/figures.sh
set -euo pipefail
cd "$(dirname "$0")/.."
# The results this script documents were produced with the earlier
# long-window defaults (262 us baseline, unconditional restart, width-derived
# noise). LEGACY reproduces them; the current defaults model the measured
# restorer instead (README, "Baseline restorer in production").
LEGACY="--baselineMeanSamples=32768 --baselineRestartThrADC=0 --lfFilterOn=off --noiseSigmaADC_L=960 --noiseSigmaADC_R=960"
N=${NEVENTS:-200000}
SIM=output/figures
mkdir -p "$SIM" doc/figs
nix build

# Long ends (20 us pole-zero) at three rates, for the baseline-slope test.
for r in 40 10 1; do
  ./result/bin/simulate $LEGACY --nEvents="$N" --rate_kHz="$r" \
    --outputFile="$SIM/long_${r}kHz.root" > "$SIM/long_${r}kHz.log" 2>&1 &
done
# The two pole-zero settings at the production rate give all four strip groups.
./result/bin/simulate $LEGACY --nEvents="$N" --rate_kHz=27 --trapDecay_us=20 \
  --outputFile="$SIM/long_27kHz.root" > "$SIM/long_27kHz.log" 2>&1 &
./result/bin/simulate $LEGACY --nEvents="$N" --rate_kHz=27 --trapDecay_us=10 \
  --outputFile="$SIM/short_27kHz.root" > "$SIM/short_27kHz.log" 2>&1 &
wait

# Every figure twice through PlottingUtils: PDF for the note, PNG for
# everything else (the applicability figures need scripts/applicability.sh to
# have run for the grid; the plane figure needs nothing).
for format in pdf png; do
  out=doc/figs
  [ "$format" = png ] && out=doc/figs/png
  mkdir -p "$out"
  ./result/bin/figures --format=$format --out="$out" \
    --long40="$SIM/long_40kHz.root" --long10="$SIM/long_10kHz.root" \
    --long1="$SIM/long_1kHz.root" --long27="$SIM/long_27kHz.root" \
    --short27="$SIM/short_27kHz.root" 2>&1 | grep -v "^Info in"
  ./result/bin/applicability --format=$format --grid=output/applicability/grid.txt \
    --out="$out" 2>&1 | grep -v "^Info in"
  [ -f "$out/log_applicability_plane.$format" ] && \
    mv "$out/log_applicability_plane.$format" "$out/applicability_plane.$format"
done

# Crop the PDFs to their content (pdfcrop comes from the LaTeX flake in doc/).
(
  cd doc
  # doc/ is not tracked by git, so its flake is addressed as a path.
  nix develop path:. --command bash -c 'for f in figs/*.pdf; do pdfcrop "$f" "$f" > /dev/null 2>&1; done'
)
ls -la doc/figs doc/figs/png
