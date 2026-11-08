# MUSIC pole-zero figures

`figures` draws the figures and writes the tables of the pole-zero paper: the
filter responses, the two-pulse baseline figure, the rate-dependent
pulse-height spectra, the spread, reaction and amplitude-band tables, and the
run-97 results. `simulate` makes the hits the
spectra read. The inputs are simulated filter responses, simulated hits, and
the outputs of the 37Cl analysis.

## Build and run

```
nix build                       # -> result/bin/{figures,simulate}
bash scripts/spectra_inputs.sh  # regenerate the 1/40 kHz spectra inputs
bash scripts/spread_inputs.sh   # regenerate the spread-table runs
bash scripts/reaction_inputs.sh # runs with reacted particles, in output/reactions
./result/bin/figures            # PDFs and the table in doc/figs, PNGs in doc/figs/png
bash scripts/copy_figures.sh    # copy the PDFs and the table into the paper
```

Inside `nix develop`, `make -j` builds the same binaries in `bin/`. The
unmodified MUSIC `PulseHistory`, `Constants` and `Plots` are compiled from the
pinned `music` flake input. To build against a local checkout, use
`nix build --override-input music git+file://$HOME/MUSIC --no-write-lock-file`
(a `path:` input would copy the whole checkout into the store).
`src/SimDatasetConfig.cpp` follows the pinned tooling, so a newer checkout may
need it updated first.

```
figures [--spread=<dir>] [--reactions=<dir>]
        [--run97=<run97.root>]
        [--long1=<long_1kHz.root>] [--long40=<long_40kHz.root>]
        [--out=doc/figs] [--format=pdf|png]
```

With no arguments it reads

| Option | Default | Contains |
|---|---|---|
| `--spread` | `output/spread` | the spread-table runs, `<reading>_r<rate>.root` |
| `--reactions` | `output/reactions` | the reaction-table runs, `<reading>_r<rate>.root` |
| `--run97` | `output/run97/run97.root` | what `scripts/run97_figure.sh` measures on run 97: the kernels, the before/after profiles, the strip spectra and their FWHM, and the band spreads |
| `--long1` | `output/figures/long_1kHz.root` | the 1 kHz simulated hits, left matched, right mismatched |
| `--long40` | `output/figures/long_40kHz.root` | the 40 kHz simulated hits, left matched, right mismatched |

A figure whose input is missing is skipped, so the previous file is left in
place. `--format` writes one format straight into `--out`; the default writes
both, with the PNGs moved to `<out>/png/`.
## Figures

| File | Content |
|---|---|
| `trapezoid` | matched pole--zero output for one pulse |
| `polezero` | matched, under- and over-compensated tails |
| `single` | one right-chain pulse with the 37Cl settings |
| `pileup_pair` | a second pulse inside the baseline freeze, with the loss at its read point |
| `spectrum_1kHz` / `spectrum_40kHz` | raw pulse-height spectra at the two rates, matched left against mismatched right |
| `spectrum_1kHz_corrected` / `spectrum_40kHz_corrected` | the same spectra after the PulseHistory correction |
| `spread_table.tex` | a complete table float, not a figure, with its caption and the label `tab:spread`: rms spread before and after correction at 10, 44 and 80 kHz for the two baseline readings (ideal, and the conditional restart of the firmware with `EnergyFilterLFLimitation` on), for the paper to `\input` (needs `siunitx`, `booktabs` and `amsmath`) |
| `reaction_table.tex` | a table float like `spread_table.tex`, label `tab:reactions`: spread about the true amplitude, before and after correction, for beam particles, reacted particles and beam particles right after a reacted one |
| `kernel_run97` | the long-left, long-right and short-left kernels fitted on run 97 |
| `run97_deviation_R` | run 97 long-right profile before and after correction |
| `run97_spectra` / `run97_spectra_corrected` | run 97's strip DeltaE spectra (the segmented strips S1-S16, each in units of its beam peak) before and after the correction, as MUSIC's events summary draws them |
| `run97_fwhm.tex` | LaTeX macros for the text, not a float: the strips' beam-peak FWHM in percent of the peak, before and after the correction, for the strips whose long segment the right chain reads (`\FwhmLongRightBefore`, `\FwhmLongRightAfter`) and those whose long segment the left chain reads (`\FwhmLongLeft...`), each the mean over those strips, with `Min`/`Max` variants for the narrowest and widest |
| `band_table.tex` | a table float like `spread_table.tex`, label `tab:bands`: the spread before correction and after it with one kernel and with four amplitude bands, for the runs of the spread table and for run 97, in percent of the beam peak |

The filter-response figures use unit-amplitude pulses, so their vertical
scales are in arbitrary units. They use the DPP-PHA settings of the 37Cl runs
97-137 (board 66222): 8 ns samples, 2 us rise, 3 us flat top, 20 us pole-zero,
50 % peaking, 7 us baseline holdoff. They are constants at the top of
`mains/figures.cpp`. The spectra read the simulated `hits` trees and fill all
non-pileup hits; each run measures its own PulseHistory kernel before the
corrected spectrum is drawn. `scripts/spectra_inputs.sh` regenerates those
inputs with the left chain at `tauL = 20 us`, matched to the 20 us pole-zero.

`spread_table.tex` gives the rms of non-pileup pulses within 0.5-1.5 times the
beam peak; its isolated column is the corrected rms at 0.5 kHz. Each run
measures its own PulseHistory kernel, as the spectra do.
`scripts/spread_inputs.sh` makes the runs: two baseline behaviours at 0.5,
10, 44 and 80 kHz, with `tauR = 10 us` and `tauL = 25 us`. The caption takes
the decay times from the runs.

`scripts/reaction_inputs.sh` makes runs in which 1 % of the particles react
upstream of the two long ends (`simulate --reactionFraction --reactionPairs`).
Their deposits come from MUSIC's 37Cl simulation: `scripts/reaction_pairs.C`
takes the (alpha,n) and (alpha,p) events reacting in strips 2-8 and every
long-left / long-right neighbour pair at or after the reaction, as ratios to
the unreacted beam there, weighted by the TALYS cross sections. The physical
rate is about 2e-5 per beam particle; the 1 % is for statistics. Reactions
draw from their own random stream, so the beam is event for event that of the
reaction-free runs of `scripts/spread_inputs.sh`. `reaction_table.tex` scores
these runs against truth: the spread of each height about its chain's gain
times the true amplitude, the gain taken from the beam particles. Only
particles with no other particle, recorded or not, within 5 us enter; the
pileup flag knows only recorded triggers, so it lets merged double pulses
through.

`band_table.tex` compares the production kernel, one per channel group, with
one kernel per band of the previous pulse height (`PULSE_HISTORY_AMP_BINS =
4`: below 0.5, 0.5-1.5, 1.5-2.5 and above 2.5 times the beam peak). Its
simulated rows correct the spread table's runs again with four
bands; its run-97 row comes from `run97.root`. Both use the spread table's
measure, in percent of the beam peak, so they read side by side.

## Run-97 input

`scripts/run97_figure.sh` builds `output/run97/run97.root` with the local
MUSIC tooling and then redraws the figures. It runs PulseHistory on the first
run-97 file and keeps the group kernels it fits, the before/after profiles and
the strip spectra. Those are built as in MUSIC's events summary: events seeded
on the grid within the event window, strips 1-16 summing their two ends, each
strip in units of its own beam peak, before and after the correction. The
kernels do not come from the 37Cl pipeline log, so they follow the pinned
PulseHistory even before that pipeline is rerun.

It also keeps each strip's beam-peak FWHM, taken from the per-event values
(0.2 % bins smoothed over 1 %, not the spectra's 0.02 a.u. bins), and which
chain reads each strip's long segment, for `run97_fwhm.tex`. For
`band_table.tex` it corrects the file a second time with four amplitude
bands. It then keeps the long groups' spread before correction, with one band
and with four bands: hits with no other hit on the same channel within 5 us,
within 0.5-1.5 times the channel's beam peak, in percent of that peak. It
also keeps the long ends' mean trigger rate.

```
bash scripts/run97_figure.sh
```

## Layout

| Path | Contents |
|---|---|
| `mains/figures.cpp` | the figures, the spread table and their input parsing |
| `mains/simulate.cpp` | the hit simulation behind the spectra |
| `include/`, `src/` | the simulator, and the preamp, trapezoid and PulseHistory adapter `figures` uses |
| `scripts/` | the spectra, spread-table, reaction and run-97 input generators and the paper-figure copy |

`output/` holds the inputs and `doc/figs/` the generated figures and table; both
are gitignored.
