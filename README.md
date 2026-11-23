# MUSIC pole-zero figures

`figures` draws the figures and writes the table and text macros of the
pole-zero paper: the filter responses, the correction figure, the simulated
pulse-height spectra and kernels, the spread table, and the run-97 results.
`simulate` makes the hits the simulated figures read. The inputs are simulated
filter responses, simulated hits, and the outputs of the 37Cl analysis.

## Build and run

```
nix build                       # -> result/bin/{figures,simulate}
bash scripts/sim_inputs.sh      # regenerate the simulated runs, in output/sim
bash scripts/run70_baseline.sh  # the baseline restorer in run 70's waveforms
./result/bin/figures            # PDFs, table and macros in doc/figs, PNGs in doc/figs/png
bash scripts/copy_figures.sh    # copy the PDFs, table and macros into the paper
```

Inside `nix develop`, `make -j` builds the same binaries in `bin/`. The
unmodified MUSIC `PulseHistory`, `Constants` and `Plots` are compiled from the
pinned `music` flake input. To build against a local checkout, use
`nix build --override-input music git+file://$HOME/MUSIC --no-write-lock-file`
(a `path:` input would copy the whole checkout into the store).
`src/SimDatasetConfig.cpp` follows the pinned tooling, so a newer checkout may
need it updated first.

```
figures [--sim=<dir>] [--run97=<run97.root>]
        [--reaction-events=<reaction_events.root>]
        [--out=doc/figs] [--format=pdf|png]
```

With no arguments it reads

| Option | Default | Contains |
|---|---|---|
| `--sim` | `output/sim` | the simulated runs of `scripts/sim_inputs.sh`, `<reading>_r<rate>.root` |
| `--run97` | `output/run97/run97.root` | what `scripts/run97_figure.sh` measures on run 97: the kernels, the before/after profiles, the strip spectra, and the band spreads |
| `--reaction-events` | `output/reaction_events/reaction_events.root` | what `scripts/reaction_events.sh` collects from MUSIC's outputs: its tagged reaction events and pure-beam events, with and without the correction |

A figure whose input is missing is skipped, so the previous file is left in
place. `--format` writes one format straight into `--out`; the default writes
both, with the PNGs moved to `<out>/png/`.
## Figures

On the standard 1200x800 canvas, the paper's figures use three text sizes:
axis text at the style's 48 px; in-figure annotations, corner notes and panel
letters at 38 px, which print at the preprint's 12 pt body size; and legends at
34 px. Lines keep the style's own width. Colour is used only where it tells
curves apart: chains, thresholds, strips' parity, and beam against reactions.
`--annotation-px=` and `--legend-px=` draw the annotations or the legends at
another size, for comparison.

| File | Content |
|---|---|
| `trapezoid` | matched pole-zero output for one pulse, in black and grey |
| `polezero` | under- and over-compensated tails |
| `single` / `correction_train` | the two panels of one figure: (a) in black and grey, one right-chain pulse with the 37Cl settings, its read point, gain g and tail; (b) a four-pulse train on the right chain with an ideal baseline, the last pulse's shift the sum of the earlier pulses' tails at its read point, κ(Δt_i) h̃_i |
| `spectrum_10kHz` / `spectrum_50kHz` | the two panels of one figure: simulated pulse-height spectra at 10 and 50 kHz, matched (left chain) against mismatched (right chain), before (dashed) and after (solid) the PulseHistory correction, with the running-average check |
| `correction_kernel` | the kernels PulseHistory fits on the 50 kHz run with the running-average check, matched and mismatched chain, against the analytic kernel of Eq. (kerneltail) |
| `spread_table.tex` | the simulated rms spread before and after correction at 10 and 50 kHz, matched and mismatched chain, for the two readings of the baseline check, in percent of the chain's own mean height. Not a figure: a table body (the tabular) for the paper to `\input` inside its own table float, which holds the caption and the label `tab:spread` (needs `siunitx`, `booktabs` and `amsmath`); a comment line at its top gives what the caption needs from the runs |
| `band_macros.tex` | LaTeX macros for the text, not a float: the spread of run 97's Long Strips, Right before correction, with one kernel and with one kernel per bin of the previous pulse height (four bins), and the long ends' rate (`\BandCount`, `\BandRun...`) |
| `run70_baseline_change` | the mean change of the firmware baseline between consecutive pulses on the Long Strips, Right of run 70 (the paper's 37Cl diagnostic run), against the time between them, with CAEN's description at 10, 20 and 80 ADC and with no check; from `output/run70/run70_baseline.root` (`scripts/run70_baseline.sh`) |
| `kernel_run97` | the long-left, long-right and short-left kernels fitted on run 97 |
| `run97_deviation_R` | run 97 long-right profile before and after correction, to 40 us between pulses |
| `run97_spectra` / `run97_spectra_corrected` | run 97's strip DeltaE spectra (the segmented strips S1-S16, each in units of its beam peak) before and after the correction, on a log scale, odd strips (long segment left) blue and even strips (long segment right) red; the two panels (a) and (b) of one figure |
| `reaction_traces_before` / `reaction_traces_after` | ten events MUSIC tagged as (α,n) in strip 3 and ten of its pure-beam events, strip by strip, before and after the correction (the same events in both) |

The filter-response figures use unit-amplitude pulses, so their vertical
scales are in arbitrary units. They use the DPP-PHA settings of the 37Cl runs
97-137 (board 66222): 8 ns samples, 2 us rise, 3 us flat top, 20 us pole-zero,
50 % peaking, 7 us baseline holdoff. They are constants at the top of
`mains/figures.cpp`.

## Simulated runs

`scripts/sim_inputs.sh` makes every simulated run the paper uses, at 10 and
50 kHz, with the paper's Table 1 settings, a 10 % Gaussian spread of the beam
amplitudes, and the same per-sample noise on both chains (the right chain's
run-70 traces), so the matched left chain (`tauL = 20 us`) and the mismatched
right chain (`tauR = 10 us`) differ only in their decay time. The baseline
restorer is CAEN's description of `EnergyFilterLFLimitation = On`: after each
trigger the baseline is frozen for 7 us; at the end of each freeze it resumes
if the trapezoid is back within 20 ADC of it, and otherwise freezes again.
CAEN does not say how close, nor what decides the check, so two readings:

- `average`: the running average of the trapezoid over the baseline window
  (1024 samples), kept up during the freeze, against the held baseline. The
  figures use this one.
- `single`: the last trapezoid sample alone. The trapezoid noise then lets the
  baseline resume while a tail is still there, and at high rate it wanders.

The 20 ADC comes from run 70 (below). Each hit carries the frozen baseline it
was read against, relative to the ideal one (`baseline` branch). The spectra
fill all non-pileup hits; each run measures its own PulseHistory kernel.
`spread_table.tex` gives the rms of non-pileup pulses within 0.5-1.6 times the
beam peak, in percent of the chain's mean height there, so the mismatched
chain's gain does not enter.

`band_macros.tex` compares, on run 97, the production kernel, one per channel
group, with one kernel per band of the previous pulse height
(`PULSE_HISTORY_AMP_BINS = 4`: below 0.5, 0.5-1.5, 1.5-2.5 and above 2.5 times
the beam peak).

## Run-70 baseline

Run 70 recorded waveforms with the trapezoid (`EnergyFilter`) as the second
probe, before the preamp swap at run 74, at about 3.5 kHz per long channel:
the right group with a 10 us pole-zero (an overshoot after every pulse, its
tail decaying with ~23 us, so the 25 us preamps), the left with a 20 us one,
near matched. `scripts/run70_baseline.sh` infers the firmware's baseline at
each traced hit, the trapezoid at the reading point less the reported energy
over its scale, and its change between consecutive hits on a channel. It fits
CAEN's description to that change, for a range of thresholds and baseline
windows, prints that the two readings of the check fit it alike at the run's
own baseline average (16384 samples), and writes
`output/run70/run70_baseline.root` (the curves of `run70_baseline_change`) and
`output/run70/baseline_change.png`. The extraction of the long-channel hits
from `data_raw/run70.root` is kept in `output/run70/run70_hits.root`.

## Run-97 input

`scripts/run97_figure.sh` builds `output/run97/run97.root` with the local
MUSIC tooling and then redraws the figures. It runs PulseHistory on the first
run-97 file and keeps the group kernels it fits, the before/after profiles and
the strip spectra. Those are built as in MUSIC's events summary: events seeded
on the grid within the event window, strips 1-16 summing their two ends, each
strip in units of its own beam peak, before and after the correction. The
kernels do not come from the 37Cl pipeline log, so they follow the pinned
PulseHistory even before that pipeline is rerun.

For `band_macros.tex` it also corrects the file a second time with four amplitude
bands. It then keeps the long groups' spread before correction, with one band
and with four bands: hits with no other hit on the same channel within 5 us,
within 0.5-1.6 times the channel's beam peak, each channel in units of its
own peak, rms in percent of their mean as for the simulated runs. It
also keeps the long ends' mean trigger rate.

```
bash scripts/run97_figure.sh
```

## Reaction events

The reaction figure shows events MUSIC's 37Cl analysis tagged as reactions,
and events of its pure-beam sample, with and without the correction: the
same events in both, using MUSIC's own code throughout. Two steps, after the
production 37Cl pipeline (its strip-sum-scatter cache must carry run numbers,
cache v46):

1. `scripts/nocorrection.sh` runs MUSIC's pipeline and reaction search with
   the correction switched off, for the late epoch (runs 97-137), into
   `output/nocorrection` (`MUSIC_RESULTS_DIR`): events, calibration,
   summaries, plots, and a strip-sum-scatter cache with MUSIC's own tags on
   uncorrected data. The production results are not touched; the pipeline
   log it writes into the dataset directory is moved here and the production
   one put back. About as long as a production run and about 25 GB.
2. `scripts/reaction_events.sh` reads the production cache and the
   no-correction events and writes
   `output/reaction_events/reaction_events.root`: every tagged event's and
   every pure-beam event's calibrated deposits with the correction (from
   the cache) and without it (the same event, joined on run and grid seed
   timestamp, read with MUSIC's `EnergyView`). A couple of minutes.

Both build MUSIC's library from the local checkout and link their program to
it.

```
bash scripts/nocorrection.sh
bash scripts/reaction_events.sh
./result/bin/figures
```

## Layout

| Path | Contents |
|---|---|
| `mains/figures.cpp` | the figures, the spread table, the text macros and their input parsing |
| `mains/simulate.cpp` | the hit simulation behind the simulated figures and table |
| `include/`, `src/` | the simulator, and the preamp, trapezoid and PulseHistory adapter `figures` uses |
| `scripts/` | the simulated-run, run-70 and run-97 input generators, the reaction-event collection and the paper-figure copy |

`output/` holds the inputs and `doc/figs/` the generated figures and table; both
are gitignored.
