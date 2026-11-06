# MUSIC pole-zero figures

`figures` draws the figures of `doc/kernel_note.tex` and the two-pulse baseline
figure. Four are simulated filter responses; the other three read the outputs
of the 37Cl analysis.

## Build and run

```
nix build                       # -> result/bin/figures
./result/bin/figures            # PDFs in doc/figs, PNGs in doc/figs/png
bash scripts/copy_figures.sh    # copy doc/figs/*.pdf into the paper
```

Inside `nix develop`, `make -j` builds the same binary in `bin/`.

```
figures [--log=<pipeline_fused.log>] [--run=run97_c000]
        [--sim=<analysis.log>[,<analysis.log>...]]
        [--run97=<deviation.root>] [--out=doc/figs] [--format=pdf|png]
```

With no arguments it reads

| Option | Default | Contains |
|---|---|---|
| `--log` | `~/MUSIC/analysis/37Cl/pipeline_fused.log` | the `kernel:` lines of the 37Cl pipeline |
| `--sim` | `output/paper_sim/analysis.log` | `analyze` logs, one per baseline reading |
| `--run97` | `output/run97/deviation.root` | run 97's `before_long_right` and `after_long_right` profiles |

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
| `kernel_run97` | the long-left, long-right and short-left kernels fitted on run 97 |
| `sim_R` / `sim_L` | corrected and raw spread against rate for the three baseline readings |
| `run97_deviation_R` | run 97 long-right profile before and after correction |

The simulated figures use unit-amplitude pulses, so their vertical scales are
in arbitrary units. They use the DPP-PHA settings of the 37Cl runs 97-137
(board 66222): 8 ns samples, 2 us rise, 3 us flat top, 20 us pole-zero, 50 %
peaking, 7 us baseline holdoff. They are constants at the top of
`mains/figures.cpp`.

## Run-97 input

`scripts/run97_figure.sh` builds `output/run97/deviation.root` with the local
MUSIC tooling and then calls `figures` for that figure:

```
bash scripts/run97_figure.sh
```

## Layout

| Path | Contents |
|---|---|
| `mains/figures.cpp` | the figures and their input parsing |
| `include/`, `src/` | the preamp and trapezoid models used by the simulated figures |
| `scripts/` | the run-97 profile generator |

`output/` is gitignored and holds analysis inputs and figure output.
