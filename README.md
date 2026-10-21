# MUSIC pole-zero mismatch simulation

This project simulates two MUSIC anode channels (one long-L, one long-R)
through a CAEN DPP-PHA chain:

preamp -> 16-bit digitizer -> Jordanov-Knoll trapezoid (with pole-zero
cancellation) -> trapezoid baseline restorer -> energy readout.

The goal: test the hypothesis that the R preamp has a 10 us decay time under
a 20 us pole-zero setting, and that the `PulseHistory` correction in the
MUSIC analysis removes the resulting pulse-height dependence on the time
since the previous pulse.

## Build and run

```
nix build                   # analysis-utils (pinned), ROOT, and all binaries -> result/bin/
./result/bin/simulate --nEvents=100000 --rate_kHz=40 --outputFile=sim.root
./result/bin/analyze --main=sim.root --outdir=analysis
bash scripts/run_all.sh     # the five runs of the task (40/10/1 kHz, control, sign check) + analysis
bash scripts/extras.sh      # linearity run and the 1024 / 16384 restorer variants -> output/extras
bash scripts/tau_scan.sh    # tau_R scan against the real run-97 kernel -> output/tau_scan
bash scripts/real_data.sh   # run-70 waveforms and hit-stream profiles -> output/real
bash scripts/kernel_study.sh # what limits the correction -> output/kernel_study
bash scripts/tau_recovery.sh # does the fit return the input tau -> output/tau_recovery
bash scripts/applicability.sh # rate x mismatch grid for the applicability figure -> output/applicability
bash scripts/figures.sh     # the figures of doc/polezero.tex -> doc/figs (PDF) and doc/figs/png (PNG)
bash scripts/paper_sim.sh   # the simulation study and figures of doc/kernel_note.tex (PDF and PNG)
bash scripts/run97_figure.sh # the run-97 before/after figure of doc/kernel_note.tex (needs ~/MUSIC)
```

`nix build` outputs only the `result` symlink (`result/bin/simulate`, ...);
nothing lands in the source tree, so there is no `bin/` or `build/` to clean.
The `Makefile` remains for interactive work inside `nix develop`.

### Raw-waveform tools

For SOLARIS runs that saved traces (the 88Sr runs under
`/labdata/MUSIC/88Sr`; see `doc/polezero.tex` and the headers of the mains):

```
result/bin/wavecheck  --settings=<XSetting.dat> --out=<dir> <sol files>   # settings table, offline trapezoid, tails, spectra, aligned shapes
result/bin/deconvolve --in=<dir>/wavecheck.root --out=<dir> [--tau2=100]  # detector current from the 10 us strips, pushed through a slower preamp
result/bin/optfilter  --settings=<XSetting.dat> --template=<dir>/wavecheck.root --out=<dir> <sol files>   # matched / optimal filter against the firmware
result/bin/applicability [--grid=output/applicability/grid.txt] [--out=doc/figs]   # the two applicability figures
```

The channel map and per-channel settings parser these share are in
`include/SolWave.hpp`.

For the baseline-restorer question (see "Baseline restorer in production"):

```
result/bin/blkernel --label=run097 [--group=chN] [--slices=a,b --scoreout=f.root] \
    [--sim=label:file.root ...] <sol files>                    # how long a pulse shifts later pulses, data vs simulated baselines
result/bin/wavenoise [--nfiles=N] <run-70 sol files>          # per-sample ADC and trapezoid noise from the traces
```

`--sim` labels must not contain a colon.

**The MUSIC tooling is a pinned flake input** (`github:ewtodd/MUSIC`), so the
build never reads `~/MUSIC` directly. After changing a header under
`~/MUSIC/tooling`, build against the checkout with
`nix build --override-input music path:~/MUSIC`; commit the MUSIC change and
`nix flake update music` to pick it up from GitHub. The whole derivation is
rebuilt when the input rev moves, so objects can never end up split across two
`Constants.hpp` versions (a mix like that crashes at exit; it can still happen
with incremental `make` inside `nix develop`).

`output/` as committed holds only the figures of the five task runs and
`analysis.log` with the summary tables. The other scripts write to their own
subdirectories.

`result/bin/simulate --name=value ...` overrides any `Config` field (see
`include/Config.hpp`).

**The defaults changed on 2026-10-12** from the old long-window baseline
(a 262 us average that restarts right after each trapezoid, noise from the
pulse-height width) to the one that fits run 97 (an 8.2 us average that
stays frozen until the trapezoid is back within 80 ADC, a 16 us smoothing,
noise measured on traces); see "Baseline restorer in production". The scripts whose documented results came from the old
defaults (`run_all.sh`, `extras.sh`, `tau_scan.sh`, `figures.sh`, and the
first half of `tau_recovery.sh`) pass the old values explicitly as `$LEGACY`,
so they reproduce the tables below. `kernel_study.sh`, `applicability.sh` and
the second half of `tau_recovery.sh` use the new defaults.

The pinned MUSIC tooling ships every `PulseHistory` group disabled
(`PULSE_HISTORY_GROUPS`); `src/SimDatasetConfig.cpp` enables the long L and R
ends with the binned kernel. Without that, `analyze` reports "Measure failed"
and applies no correction.

`result/bin/analyze --main=f40,f10,f1 --control=f --signcheck=f
--extra=label:f,...` analyses the runs. It measures the `PulseHistory`
kernel on the first `--main` file and applies it to the other two. It
measures a kernel on every other file. Figures go to `output/<label>/` and
`output/pulse_history/<label>/` (the latter drawn by `PulseHistory::SavePlots`
itself). The summary tables go to stdout (`output/analysis.log`).

The unmodified `tooling/src/PulseHistory.cpp` and `Constants.cpp` of the
pinned MUSIC input are compiled into the binaries.
`src/SimDatasetConfig.cpp` plays the role of the per-dataset `Constants.cpp`
of the monorepo.
`src/PulseHistoryAdapter.cpp` turns the simulated hits into the `RawHit`
stream `PulseHistory` expects.

## Parameters and where they come from

Runs 97-137 of the 37Cl data set are on board 66222. Its per-run settings
(`/labdata/MUSIC/New-37Cl/data_raw/music_exp1915_<run>XSetting_66222.dat`)
are identical over all 41 runs.

The file named in the task, `Settings/setting_66226_PHA.dat`, belongs to
board 66226. That board only served runs 0-5 with different filter settings
(5 us rise, 1 us flat, 50 us pole-zero). So the 66222 values are used.

| Config field | Value | Source |
|---|---|---|
| `samplePeriod_ns` | 8 | V2740, ADC_SamplRate = 125 MS/s |
| `trapRise_ns` / `trapFlat_ns` | 2000 / 3000 | EnergyFilterRiseTimeT / FlatTopT, all long ends |
| `trapDecay_us` | 20 | EnergyFilterPoleZeroT = 20000 ns on every long end (the 16 short ends are set to 10 us, the Grid to 200 us) |
| `peakingDelay_ns` | 1500 | EnergyFilterPeakingPosition = 50 % of the flat top |
| `peakingAvgSamples` | 16 | EnergyFilterPeakingAvg = MediumAVG = 16 samples (SOLARIS_DAQ `DigiParameters.h`) |
| `baselineMeanSamples` | 1024 | EnergyFilterBaselineAvg = Medium = 1024 samples = 8.2 us at 8 ns; the 2740 PHA has no decimation (SOLARIS_DAQ `DigiParameters.h`: Fixed/VeryLow/Low/MediumLow/Medium/MediumHigh/High = 0/16/64/256/1024/4096/16384). Until 2026-10-12 the default was 32768 (262 us), which many results below were produced with; run 97 rules it out ("Baseline restorer in production") |
| `baselineHoldoff_ns` | 7000 | After each trigger the baseline stops averaging for the trapezoid length (2 rise + flat) plus EnergyFilterBaselineGuardT, which is 0 in every 37Cl run (confirmed by CAEN) |
| `baselineRestartThrADC` / `baselineRestartHoldSamples` | 80 / 250 | After that, the baseline stays frozen until the trapezoid has stayed within this many ADC of it for this many samples in a row (250 = 2 us). This models what CAEN describes for EnergyFilterLFLimitation = On; 80 ADC / 2 us fits run 97 and the ch-4 overshoot runs. 0 starts averaging right after the freeze, as with LFLimitation = Off |
| `baselineRefreeze` | off | On: the check is made once at the end of each freeze and, if the trapezoid is not back within `baselineRestartThrADC`, the baseline is frozen for another trapezoid length (CAEN's own description). Fits run 97 at 10-20 ADC up to 38 us, but holds the baseline too long beyond 80 us |
| `baselineIdeal` | off | On: subtract the true baseline (the trapezoid of the ADC offset alone) instead of the restorer's, so earlier pulses never enter it. Used to test the correction without assuming what the restorer does |
| `lfFilterOn` / `lfFilterTauUs` | on / 16.384 | A smoothing (one-pole low-pass) on the baseline. CAEN says the firmware has none; it is on only because it reproduces a small, unexplained effect of pulses 80-500 us earlier in run 97. 65 us or more is ruled out |
| `baselineDecimation` | 1 | Baseline averaged on a clock divided by N (32 was an old reading, ruled out) |
| `triggerHoldoff_ns` | 1000 | TimeFilterRetriggerGuardT = 0; the real kernels are populated from 1 us, so the time filter's re-arm is taken as 1 us |
| `adcBits` / `adcOffset` | 16 / 32768 | ADC_Nbit = 16, DCOffset = 50 % (the task's 14 bit does not match the board) |
| `beamPeakADC` | 1600 | long-end beam peak in Events_Run100_1_c000.root (L1 1583, L3 1644, L5 1691 ADC) |
| `noiseSigmaADC_L/R` | 210 / 563 | white noise per sample giving the trapezoid noise measured on the run-70 traces (19 / 51 in energy units; see below). The earlier default, 960, was derived from the pulse-height width (the measured L width is ~105 ADC, of which 32 ADC is the 2 % beam spread) |
| `tauL_us` / `tauR_us` | 25 / 10 | read off the hardware for runs 97+: the L chain is a 25 us preamp variant, the R chain a 10 us one, both under a 20 us pole-zero |

Runs 97+ carry no waveforms, so the default noise is fixed by the
pulse-height width. Run 70 did save ADC and trapezoid traces, and
`result/bin/wavenoise` measures the noise on isolated (200 us), single-trigger
traces, before the trigger:

| | ADC rms per sample | trace ADC / energy | trapezoid rms (energy units) |
|---|---|---|---|
| L long | 252 | 0.362 | 19 |
| R long | 275 | 0.360 | 51 (63 / 1.24) |

Most of the ~700 energy units of per-sample noise lie above the trapezoid's
band. The simulated trapezoid passes 0.0905 x a white per-sample noise, so
the measured trapezoid noise corresponds to `--noiseSigmaADC_L=210
--noiseSigmaADC_R=563`. This matters wherever a model compares the trapezoid
with a threshold (the return check below); the width-derived 960 gives a
trapezoid noise of ~87. The R value comes from run 70's R chain (the pre-swap
preamp under a 10 us pole-zero).

## Model notes

- The preamp impulse response is the closed form of a 1 us linear ramp
  convolved with the exponential. The waveform is generated sample by sample.
  Finished pulses are carried as one exponentially decaying state variable,
  so the sum over all past pulses is exact at O(1) per sample. Idle stretches
  are skipped analytically. The trapezoid (a finite-impulse response of
  rise + rise + flat samples) and the baseline window are rebuilt from a
  warm-up before the next pulse. This is exact in distribution. The warm-up is
  the filter and baseline window plus `warmupMarginSamples` (default 0, only
  defensive slack). A baseline that waits for the trapezoid to return can
  stay frozen longer than its window, so without the smoothing 500 us is
  added.

- Trigger holdoff: a particle inside 1 us of the previous accepted one gives
  no trigger. Its charge still lands in the accepted pulse. These merged
  double pulses form the blob at twice the beam peak (4 % of hits at 40 kHz).
  Means and RMS values are taken over 0.5-1.5 x beam peak, like the
  beam-like selection of `PulseHistory`.

- Pileup flag: set on both members of a pair closer than rise + flat = 5 us,
  as the DPP-PHA does. Spectra and profiles exclude flagged hits. The
  `PulseHistory` correction is run on every accepted hit, as on real data.

- The tree stores `timestamp` in picoseconds (`RawHit` convention: the SOL
  reader multiplies the 2740 nanosecond stamps by 1000), `energy` unrounded,
  `dtPrev`/`dtNext` in microseconds. The adapter rounds energies to `UShort_t`
  and adds a synthetic Grid hit (2000 ADC, inside the reference window) at
  every event. `PulseHistory::Measure` seeds beam-like events from the
  reference channel and requires every long end inside its beam window.

## What the simulation shows

Both long-end groups run a 20 us pole-zero. Neither preamp matches it.
The residual left on the trapezoid after a pulse scales as
(1 - tau_set/tau_real). So the two chains behave oppositely:

- **R, tau = 10 us:** the pole-zero over-compensates. The trapezoid
  undershoots by about 0.36 A exp(-t / tau_R). The next pulse on the channel
  reads LOW by an amount decaying with tau_R.

- **L, tau = 25 us:** the pole-zero under-compensates. The trapezoid
  overshoots. The next pulse reads HIGH. The factor (1 - 20/25) = 0.2 against
  (1 - 20/10) = -1.0 makes the L effect about five times smaller and of the
  opposite sign.

The baseline restorer decides how much of that residual survives into the
energy. While the baseline is frozen from the previous trigger (dt < 7 us),
the residual enters directly. After that, the energy is the residual at the
peaking point minus whatever the baseline window has averaged in. A baseline
that averages the tail in would cancel most of it and flip the sign at large
dt. The production data show no such sign change, so the profile follows the
residual itself. The old default reached that with a 262 us window. The
hardware instead keeps the tail out of a short window: the baseline stays
frozen until the trapezoid is back near the baseline, as CAEN confirms (see
"Baseline restorer in production"). The energy-vs-dt profile
is the same either way; quantities that depend on how long the baseline
remembers a pulse (the rate correlation below, the bias in the next section,
the reach of the kernel) are not.

At 40 kHz that gives an R deficit of about 250 ADC at the shortest dt,
recovering with tau_R. The L excess is roughly five times smaller. Both widen
the raw pulse-height distribution. The `PulseHistory` kernel removes both:
the R RMS drops from 148 to 103 ADC. The L RMS drops from 104 to 96 ADC.
A matched channel (the control run) is at 92 ADC.

Against the real 37Cl profiles the R side matches band by band. The L excess
is about twice the measured one: at dt 5-6 us the model gives +0.058 of the
beam peak and run 97 gives +0.031. An effective tau_L near 22 us would
reproduce that. The defaults keep the nominal parts. The decay seen at the
digitizer is set by the cabling and the loading as well as by the preamp,
and will not match the nominal value exactly.

Fitting a single exponential to the real profiles does not measure tau
cleanly. Besides the pole-zero tail they carry a rate correlation: a hit with
a long gap behind it sits in a quiet stretch. The long baseline average there
is less depressed by pile-up. So the profile turns positive at large dt (long
R, relative to pulses whose previous pulse was 120-150 us earlier, runs
-0.183 at 5-6 us through zero near 30 us to +0.014 at 38-60 us).
`result/bin/realdt` scans the start of the fit and shows the resulting instability.
Comparing the whole shape against the simulation is the meaningful test.
(This explanation assumes the old long baseline average. The baseline that
fits run 97 forgets a pulse mostly within ~100 us. The +0.014 at 38-60 us
matches in size the small effect of pulses 80-500 us earlier described in
"Baseline restorer in production", whose cause is not known.)

## Does the fit return the decay time it was given?

The long ends run a 20 us pole-zero and the short ends a 10 us one.
Simulating both settings with the same two preamps produces all four strip
groups with a known answer.

`result/bin/taustudy` (`scripts/tau_recovery.sh`) fits each one. It bootstraps the
error over the hits rather than trusting the fit's covariance. It adds the
spread over the start of the fit as a systematic.

| group | true | pole-zero | p0 + p1 exp(-dt/p2) | with + p3 dt |
|---|---|---|---|---|
| long L, 25 us preamp | 25 | 20 | 20.96 +- 0.49 +- 0.85 | 25.47 +- 1.59 +- 0.58 |
| long R, 10 us preamp | 10 | 20 | 8.59 +- 0.07 +- 1.00 | 9.85 +- 0.12 +- 0.19 |
| short L, 25 us preamp | 25 | 10 | 21.33 +- 0.14 +- 1.51 | 28.20 +- 0.51 +- 1.90 |
| short R, 10 us preamp | 10 | 10 | degenerate | degenerate |

The model the task specifies is biased about 15 percent low on every group
that has a mismatch, by four standard deviations or more.

The bias is the trapezoid baseline. It is an average over a few hundred
microseconds. A hit with a long gap behind it finds that average less
depressed by pile-up. So the profile is lifted roughly in proportion to the
gap. Fitting that lift with a linear term recovers the input: long L and
long R agree with truth inside one combined sigma. Short L is 1.6 sigma high.
`Analysis.cpp` now uses the four-parameter form. This is a deliberate
departure from the task's specification.

*Caveat:* this study ran with the 262 us window. Re-run with the measured
restorer (`--baselineMeanSamples=1024 --baselineRestartThrADC=80
--baselineRestartHoldSamples=250 --lfFilterOn=1 --lfFilterTauUs=16.384
--noiseSigmaADC_L=210 --noiseSigmaADC_R=563`, otherwise as
`scripts/tau_recovery.sh`), neither form recovers the input:

| group | true | p0 + p1 exp(-dt/p2) | with + p3 dt |
|---|---|---|---|
| long L | 25 | 11.00 +- 0.09 +- 2.73 | 17.58 +- 0.25 +- 3.55 |
| long R | 10 | 9.07 +- 0.05 +- 1.58 | 11.17 +- 0.10 +- 1.30 |
| short L | 25 | 24.98 +- 0.09 +- 1.25 | 38.42 +- 0.46 +- 1.32 |

The linear term is a property of the long window. On long L the overshoot
drops below the restart threshold part-way through its decay and the
baseline follows the rest, so both fits fail. Until the restorer is pinned
down, get tau by comparing the whole profile with a simulation of it (as
for tau_R = 11 us in "What the CAEN manual says").

Short R is matched (10 us into a 10 us pole-zero), so it has no tail and the
fit is correctly degenerate. That is worth holding on to: the short R ends of
run 97 behave the same way in the real data. That is independent support for
the R preamp being a 10 us part.

## What the simulation says about the correction

`result/bin/kernelstudy` (`scripts/kernel_study.sh`) measures what limits the
correction. The simulation knows the true amplitude of every pulse.
Each variant is fitted on the first half of a run and scored on the second.
A model with more parameters cannot flatter itself.
The log also scores every variant on PulseHistory's beam-like events only
(the other long end within `PULSE_HISTORY_BEAM_LO/HI` of its peak), and
`KERNELSTUDY_PH_REPORT=1` prints PulseHistory's own kernel report.

R channel, 40 kHz, RMS in ADC over the beam window:

| correction | parameters | 192k hits | 7.7k hits |
|---|---|---|---|
| none | 0 | 148.7 | 151.0 |
| PulseHistory at the time, 12 log bins to 100 us | 13 | 102.7 | 101.7 |
| 48 log bins to 100 us | 49 | 101.8 | 102.9 |
| 192 log bins to 100 us | 193 | 101.9 | 106.8 |
| **12 log bins to 316 us** | 13 | **100.1** | **100.1** |
| **15 log bins to 316 us** | 16 | **99.8** | **99.8** |
| 5 exponentials + freeze, to 100 us | 7 | 101.8 | 101.5 |
| 7 exponentials + freeze, to 300 us | 9 | 99.0 | 99.3 |
| the same, given noiseless previous amplitudes | 9 | 95.5 | 95.7 |
| matched pole-zero, nothing to correct | - | 91.7 | 91.7 |

**The reach is the binding limit, and it is two constants.**

Keeping the binned kernel exactly as it is and only letting it look further
back takes the R RMS from 102.7 to 100.1. Adding three bins to preserve the
bin width takes it to 99.8. On the L channel the same change gives
95.4 -> 93.3. The trapezoid baseline is an average over a few hundred
microseconds. Pulses beyond the 100 us the code keeps still move it.

**With the baseline that fits run 97 this no longer holds.** The premise is
the 262 us window. How much an earlier pulse still shifts the energy,
measured on the run-97 hit stream, is ~0.003 of a pulse at 40-80 us, falls to ~0.0005 by 200 us and to
the noise beyond ~330 us, about a tenth of what a 262 us window gives.
`scripts/kernel_study.sh` now runs with the current defaults (8.2 us window,
frozen until within 80 ADC, 16 us smoothing, measured noise); the numbers
below are in `output/kernel_study_meas/`:

| correction, 40 kHz | params | R 192k hits | R 7.7k hits | L 192k hits | L 7.7k hits |
|---|---|---|---|---|---|
| none | 0 | 130.3 | 132.7 | 60.9 | 60.3 |
| PulseHistory (pinned tooling: 30 log bins to 316 us), 1 amplitude band | 31 | **71.1** | **70.3** | **46.5** | **44.0** |
| 12 log bins to 100 us | 13 | 72.1 | 70.9 | 44.0 | 42.9 |
| 48 log bins to 100 us | 49 | 71.0 | 70.4 | 43.8 | 42.7 |
| 12 log bins to 316 us | 13 | 72.6 | 71.4 | 43.9 | 42.7 |
| 15 log bins to 316 us | 16 | 72.0 | 70.9 | 43.7 | 42.7 |
| 7 exponentials + freeze, to 300 us | 9 | 70.8 | 69.9 | 43.4 | 42.4 |
| the same, given noiseless previous amplitudes | 9 | 65.0 | 64.2 | 42.9 | 41.5 |
| matched pole-zero, nothing to correct | - | 60.4 | | 37.3 | |

(The absolute values are lower than in the table above because the noise is
now the measured one. Full output, including the amplitude-band rows below,
in `output/kernel_study_meas/kernel_study.log`.)

- **Reach past 100 us gains nothing on R (-0.1 to +0.5 ADC) and ~0.3 ADC
  on L.** There is little memory past 100 us to correct. The pinned tooling
  already reaches 316 us (30 bins, `kLogHi` = -3.5), which is harmless but
  not needed.
- The correction still removes most of the effect (R 130 -> 71 against a
  floor of 60).
- **The "noiseless predictor" row is not a clean noise test.** It feeds
  each earlier hit's true amplitude, which for a merged pile-up hit is the
  first particle's alone, so it also changes how pile-ups are weighted. Most
  of its gain is that, not noise (see "Amplitude bands" below). With pile-ups
  modelled separately the noise part is ~1.6 ADC (65.5 vs 63.9 on R).
- An exponential basis gains 1.3 ADC on R, still small.
- PulseHistory on L (46.5) scores worse than the simpler fits (44.0). Its L
  kernel is fitted only on events whose R end is within 0.8-1.3 x its peak
  (`PULSE_HISTORY_BEAM_LO/HI`); scored on those events alone it gives 42.8,
  in line with the rest. The difference is L hits from events the fit never
  sees.

### Amplitude bands (`PULSE_HISTORY_AMP_BINS`)

A kernel per band of the earlier pulse's amplitude ([0, 0.5), ~1, ~2,
>= 2.5 x the beam peak) would follow a shift that is not linear in that
amplitude, e.g. a baseline freeze that lasts longer after a bigger pulse. The amplitude mix is the same in the simulation and in run 97 (long
ends: ~85 % at ~1x, ~13 % at ~2x, ~1.5 % at >= 2.5x, 1 % below 0.5x;
`blkernel` and `kernelstudy` print it).

**In the simulation 4 bands help R noticeably:** PulseHistory 71.1 -> 65.8,
7 exponentials 70.8 -> 65.5 (floor 60.4). Two bands ([0, 0.5) and the rest)
give nothing. The gain is the same with the long window (70.2 -> 65.4) and
with the restart rule off (73.4 -> 69.0), so it is not the restorer. It
comes from the ~2x hits, which in the simulation are two particles merged
into one trigger: per unit of measured energy they move the baseline only
~0.6 x as much as a 1x hit (PulseHistory R kernel at 5 us: -0.119 vs
-0.182).

**On real data the bands barely matter.** PulseHistory run on all of run
97 with the 37Cl configuration (the 4-band result reproduces the production
log exactly):

| run 97, beam-like rms (ADC) | none | 1 band | 4 bands |
|---|---|---|---|
| long R | 186.4 | 110.0 | 109.1 |
| long L | 151.1 | 144.3 | 143.8 |

The production form fit gives the per-band tail amplitude on long R as
c = 0.305 / 0.319 / 0.314 for ~1x / ~2x / >= 2.5x: linear in amplitude
within 5 %. There is a hint of a longer tail after bigger pulses (10-40 %
more at 16-30 us), as the restart rule would give, but it is worth < 1 ADC.
The 37Cl configuration therefore uses a single band (changed 2026-10-12;
it was 4).

So the simulation forms ~2x hits differently from the V2740: it drops a
trigger inside the holdoff and folds the particle into the previous hit,
whereas the data have hits 1-3 us after the previous one, which the
simulation never produces. Simulated gains that hinge on pile-up hits should
be checked against data.

The rate x mismatch grid (`scripts/applicability.sh`, now on the new
defaults; numbers in `output/applicability_meas/grid.txt`, figure
`doc/figs/applicability_grid.pdf`) puts the corrected spread at 1.43 x the
matched floor at 40 kHz on a 10 us preamp under a 20 us pole-zero, against
1.07 with the long window. That is mostly the floor: with the measured noise
it is 37 ADC instead of 93. The residual the correction leaves,
sqrt(corr^2 - floor^2), is 38 ADC against 36 before (58 against 54 at a
30 us pole-zero, 84 against 80 at 80 kHz). The overshooting 8 us pole-zero
is now the weakest case at low rate (4-9 % above the floor at 5-20 kHz).

The analysis below, made with the old 262 us baseline, is kept for
reference. It was written against
the earlier tooling (12 bins to 100 us); the pinned tooling has since made
the reach change it proposes (30 bins to 316 us), so the tree-shape
consequences below have already been taken.

In `tooling/include/PulseHistory.hpp`:

```
const Int_t kNBins = 15;       // was 12
const Double_t kLogHi = -3.5;  // was -4.0, i.e. 100 us -> 316 us
```

And in the dataset config, `PULSE_HISTORY_APPLY_MAX_US` from 100.0 to 316.0.

Nothing else has to change. `BinOf`, `BinCentreUs`, the feature indices, the
SVD solve, the report and the plots are all written in terms of `kNBins` and
`kLogHi`.

Two consequences to be aware of:

1. The `pulse_history` tree written into the events file carries
   `Kernel[kMaxAmpBins * kNBins]` and `BinCentreUs[kNBins]`. Those arrays
   grow from 72 and 12 to 90 and 15. Old and new files differ in shape.

2. Each channel's deque holds about three times as many past pulses. The
   history loops in `Measure` and `Apply` cost proportionally more. At 46 kHz
   that is roughly 15 entries instead of 5.

The comment on `PULSE_HISTORY_APPLY_MAX_US` says the fitted bins beyond about
100 us are degenerate with the intercept. The degeneracy is real. The SVD
already handles it. The test that matters is out-of-sample: the score
improves and does not fall apart at low statistics.

**A finer binning is not the answer.** Twelve bins to 192 buys one ADC at
full statistics and costs five at low statistics, where the far bins are
fitted on too few pulses.

**An exponential basis is a further 0.3 to 1 ADC, not worth the surgery.**
Coefficients on exp(-dt/tau) over a geometric ladder of tau, plus an indicator
for the dt range where the baseline is still frozen, reach 99.0 with 9
parameters. That is better than 15 bins by about 0.8 ADC. It would need
`Feat`, `BinOf`, `Shift`, the report and the diagnostic plots all reworked.
Worth keeping in mind for the short ends, whose kernels are the least well
determined. But the constants above get almost all of it for free.

What is left after that is about 3.5 ADC that a noiseless predictor removes
and about 3.8 ADC that no history model of this form removes. The first was
read as the noise of the previous pulse's measured energy, but the
noiseless predictor also reweights merged pile-up hits (see "Amplitude
bands"), so it overstates the noise.

## For the paper: the kernel under three readings of the baseline

The fitted kernel needs only what is certain about the filter and the
firmware: the trapezoid is linear, the energy is read at a fixed time after
the trigger, and the baseline is an average of earlier trapezoid samples. How
the baseline treats the tail (CAEN's return check, whose threshold is not
documented) changes the kernel's shape but not its form, so the paper tests
the correction under three readings of it instead of modelling it
(`scripts/paper_sim.sh`; the draft is `doc/kernel_note.tex`, built with
`cd doc && nix develop path:. -c latexmk -pdf kernel_note.tex`):

- **never takes in the tail**: `--baselineIdeal=on` subtracts the true
  baseline, what the trapezoid gives for the ADC offset alone;
- **frozen until within 80 ADC for 2 us**: CAEN's description with a
  threshold that reproduces the run-97 energy-vs-gap profile;
- **averaging restarts after the 7 us freeze**: the baseline takes in as
  much of the tail as it can.

The first and last are the two extremes of how much tail the baseline can
take in. Spread of non-pile-up pulses near the beam peak (ADC), 44 kHz,
tau_R = 11 us, with the same `PulseHistory` code as the analysis:

| baseline | R before | R after | R isolated | L before | L after | L isolated |
|---|---|---|---|---|---|---|
| never takes in the tail | 122.9 | 68.4 | 58.2 | 62.4 | 43.5 | 37.6 |
| frozen until within 80 ADC | 127.2 | 77.4 | 65.9 | 62.6 | 49.0 | 39.3 |
| restarts after the freeze | 117.7 | 80.7 | 65.7 | 50.3 | 43.2 | 39.3 |

("isolated": the same at 0.5 kHz.) After correction every reading is within
10-25 % of isolated pulses; the excess left is similar in all three, so it is
a limit of the kernel itself, not of the baseline. The kernel explains 0.64,
0.59 and 0.51 of the long-R variance; on run 97, at the same rate, 0.65.

One reading of CAEN's description defeats the correction: a single-sample
check once per trapezoid length at 20 ADC (`--baselineRefreeze=on`) keeps
the baseline stuck for hundreds of us; R only falls from 191 to 165 ADC and
the kernel explains 0.27 of the variance. Run 97 rules it out.

`result/bin/notefigures` draws the note's figures from the analyze logs and
the 37Cl pipeline log (the run-97 kernels), as PDF (`doc/figs`) and PNG
(`doc/figs/png`). It prints its PDFs with ROOT's "Portrait" option:
`PlottingUtils::SaveFigure` puts a side-legend canvas, much wider than tall,
on a rotated page. `scripts/run97_figure.sh` makes the run-97 before/after
figure by running `PulseHistory` with the 37Cl configuration from `~/MUSIC`.

Only the small files these need are kept in `output/`:
`paper_sim/analysis.log` (the study's spreads), the two applicability
`grid.txt` tables and `kernel_study_meas/kernel_study.log` (numbers quoted
above). Everything else there is regenerated by the scripts.

## What CAEN says

### The manuals

UM7788 rev. 9, x2740 FELib PHA Parameters:

- `EnergyFilterBaselineAvg`: Fixed / VeryLow / Low / MediumLow / Medium /
  MediumHigh / High = 0 / 16 / 64 / 256 / 1024 / 4096 / 16384 samples.
  Medium, the 37Cl setting, is 1024 samples = 8.2 us.

- `EnergyFilterPeakingAvg`: OneShot / LowAVG / MediumAVG / HighAVG =
  1 / 4 / 16 / 64 samples.

- `EnergyFilterBaselineGuardT`: "freeze the baseline calculation beyond the
  trapezoid end". It is a time in ns (0-8000 in steps of 8, in the firmware's
  own parameter list and in SOLARIS_DAQ `DigiParameters.h`), not an on/off
  switch. It is 0 on all 64 channels in every board-66222 run (6-137),
  including all of production (97-137).

- `EnergyFilterLFLimitation`: "Allows to enable a low-frequency filter for
  the energy filter", and nothing more. It is On in every 37Cl channel and Off
  in the SOLARIS_DAQ defaults. The x725/x730 register manual (UM5678, reg
  0x1nA0 bit 29) describes the same bit as "optimization of the Baseline
  Restorer to avoid tails in the energy peaks".

- The event flags carry an "Energy filter baseline freeze" bit (0x3). Hits
  acquired while the baseline was frozen are identifiable in the raw data.

UM5960 sec. 4.2.4 (CoMPASS) describes the baseline as the average of a set
number of trapezoid samples before the trapezoid, frozen while the trapezoid
lasts: "once the trapezoid is returned to the baseline, the averaging
restarts to run". UM5960 and the CoMPASS QuickStart GD6300 (rev. 27) do not
mention the LF limitation.

The firmware image does not explain the LF limitation either (the
parameter-to-register mapping is in the FPGA gateware; see the local
`FIRMWARE_GUIDE.md`).

### Support ticket (October 2026)

We asked CAEN support what the baseline does after a pulse on this firmware
(V2740, CUP 2025102402, FPGA 1.0.112), and what `EnergyFilterLFLimitation`
does. The answers:

- After each trigger the baseline is frozen for rise + flat + rise (7 us
  here), plus `EnergyFilterBaselineGuardT`, and then averaging starts again
  (first reply, which did not mention the LF limitation).
- With `EnergyFilterLFLimitation` On, the baseline stays frozen for about one
  trapezoid length. The firmware then checks whether the signal has come back
  close to the baseline. If it has not, the baseline is frozen for about
  another trapezoid length, and so on. A correctly pole-zero-compensated
  signal is back by the first check, so the repeat only happens when the
  pole-zero is wrong.
- There is no low-pass filter. The running average is the only filtering of
  the baseline.

CAEN did not say how close "close to the baseline" is, or whether the check
looks at one sample or several. "Baseline restorer in production" below
compares this with run 97.

## What the real data say

Run 70 waveforms, hit streams of runs 72-97.

Run 70 saved a 32 us waveform for every event (ADC input and trapezoid
output probes). `result/bin/realwave` averages the traces of 230k isolated
pulses near the beam peak on the long ends and fits the preamp decay
directly:

| | tau from ADC input | trapezoid tail (pole-zero at run 70) |
|---|---|---|
| L long ends (8 channels) | 19.9 us, channels 19.3-20.3 | flat (PZ 20 us) |
| R long ends (8 channels) | 19.0 us, channels 18.5-20.0 | +0.35 overshoot, tau 17 us (PZ 10 us) |

At run 70 both preamps had a ~19-20 us decay time.

`result/bin/realdt` then follows the energy-vs-dt profile of each strip group
through the hit streams (no waveforms needed), relative to pulses whose
previous pulse was 120-150 us earlier:

| run | R long pole-zero | R long excess at 8-12 us | R/L peak |
|---|---|---|---|
| 72, 73 | 10 us | +0.20 (19 us preamp under-compensated) | 1.00 |
| 74-90 | 10 us | +0.02-0.03 (about matched) | 0.94-0.97 |
| 91-97 | 20 us | -0.12 to -0.14 (deficit) | 0.84 |

The R preamp was therefore replaced between runs 73 and 74 by a unit with a
decay time near 10-12 us. It stayed roughly matched to the 10 us pole-zero
until run 90. At run 91 the R long ends were set to 20 us (and their baseline
from High to Medium). That produced the deficit seen in runs 97+.

The L short ends went the other way at run 89 (pole-zero 20 -> 10 us on a
20 us preamp) and show the corresponding excess.

The simulation reproduces all three R configurations with tau_R = 19 us before
run 74 and 11 us after, a baseline that does not follow the undershoot
(High), and nothing else changed (sim files `output/realdt/sim_run*like.root`).

Two consequences for the 37Cl analysis:

1. The R-side kernel of runs 97-137 is the pole-zero undershoot of a
   ~10-12 us preamp under a 20 us setting. That is exactly the hypothesis.
   But it is a property of the replacement preamp, not of the original.

2. The R long-end pulse height at the peaking point is about 9 % low in
   that configuration (flat top decays). That is part of the R/L gain
   difference.


### R preamp decay time

The run-97 R kernel (what `PulseHistory` measures) falls more slowly than a
baseline that starts averaging again right after the 7 us freeze allows, and
it never changes sign. `scripts/tau_scan.sh` simulates tau_R = 8-16 us with
the High (16384-sample) baseline, which does not follow the undershoot
either, and compares:

| dt bin | real R (run 97) | tau_R 8 | tau_R 10 | tau_R 12 | tau_R 14 |
|---|---|---|---|---|---|
| 6 us | -0.119 | -0.182 | -0.130 | -0.091 | -0.061 |
| 8 us | -0.099 | -0.124 | -0.093 | -0.068 | -0.047 |
| 12 us | -0.074 | -0.073 | -0.061 | -0.048 | -0.034 |
| 18 us | -0.046 | -0.033 | -0.033 | -0.028 | -0.021 |
| 26 us | -0.021 | -0.005 | -0.009 | -0.010 | -0.009 |

The real kernel lies between tau_R = 10 and 12 us in size, and falls
somewhat more slowly than either. So the R preamp decays in roughly 11-12 us,
not 20 us.

That scan used 40 kHz. Repeated at the measured 44 kHz with the measured
noise, as the rms difference from run 97 over the six gap bands from 5 to
38 us (each band measured relative to pulses whose previous pulse was
38-60 us earlier):

| model | tau_R 9 | 10 | **11** | 12 | 13 | 14 |
|---|---|---|---|---|---|---|
| wait for return (80 ADC) + 16 us smoothing | 0.029 | 0.016 | **0.006** | 0.019 | 0.036 | 0.053 |
| fixed 27 us freeze + 16 us smoothing | 0.030 | 0.018 | **0.007** | 0.018 | 0.034 | 0.050 |
| High 16384, freeze 7 us | 0.020 | 0.010 | 0.014 | 0.028 | 0.042 | 0.056 |

At tau_R = 11 us the first model gives -0.203 / -0.177 / -0.136 / -0.086 /
-0.044 / -0.012 against the data's -0.197 / -0.184 / -0.142 / -0.087 /
-0.037 / -0.006. **tau_R = 11 us** whichever model is used. The High
model's profile is similar, but it remembers a pulse far too long (below).

### Baseline restorer in production (runs 91+)

No waveform run has a Medium baseline on a mismatched channel. Runs 37-41
and 70 have Medium only on the matched L ends and High on R, and the LF
limitation is On in every run. So the question is answered from the
production hit streams.

**Pulse height against the gap to the previous pulse.** R long ends, run 97
(`realdt`, first 30 M hits), in percent of the beam peak, relative to pulses
whose previous pulse was 120-150 us earlier. By then the earlier pulse no
longer matters: in the simulation that level is within 0.1 % of a pulse with
no earlier pulse at all. The simulations run at the run-97 rate (44 kHz),
tau_R = 11 us and the measured noise, 400k events each:

| gap [us] | 5-6 | 6-8 | 8-12 | 12-18 | 18-26 | 26-38 | 38-60 | 60-90 | 90-120 | rms |
|---|---|---|---|---|---|---|---|---|---|---|
| **run 97** | -18.3 | -17.0 | -12.9 | -7.4 | -2.4 | +0.8 | +1.4 | +0.6 | +0.1 | |
| averaging starts again right after the 7 us freeze | -12.0 | -9.6 | -2.3 | +6.0 | +4.9 | +2.0 | +0.6 | +0.1 | +0.0 | 7.0 |
| CAEN's description, one check per trapezoid length, 10 ADC | -18.8 | -16.2 | -12.1 | -7.0 | -3.2 | -0.7 | +0.3 | +0.7 | +0.4 | 0.8 |
| same, 20 ADC | -18.4 | -15.8 | -11.4 | -6.1 | -2.3 | -0.2 | +0.5 | +0.6 | +0.4 | 0.9 |
| same, 40 ADC | -17.9 | -15.3 | -10.7 | -5.1 | -1.7 | -0.0 | +0.1 | +0.2 | -0.1 | 1.3 |
| frozen until within 80 ADC for 2 us | -18.5 | -15.9 | -11.9 | -7.0 | -2.5 | +0.6 | +0.7 | +0.3 | +0.3 | 0.6 |
| same + 16 us smoothing (the default) | -18.5 | -15.9 | -11.9 | -7.0 | -2.8 | +0.3 | +1.4 | +0.7 | +0.2 | 0.6 |

(rms: difference from run 97 over the nine bands, in percent. The statistical
error on each simulated band is 0.2-0.5 %.)

Why the first row reads high from 12 us on: a pulse that follows another sits
on its undershoot, which makes it read low. But if the baseline starts
averaging again right after the freeze, it takes in the undershoot of the
last 8 us, which was deeper than the undershoot still under the pulse, and
subtracting that makes the pulse read high. While the undershoot is
recovering the second effect wins. Run 97 reads low at every gap, so its
baseline does not follow the undershoot. That is what CAEN describes for
`EnergyFilterLFLimitation` On: the baseline stays frozen until the trapezoid
is back close to the baseline. Every model that waits for that comes close
(rms 0.6-1.3 %, against 7 %). CAEN's own version (one check per trapezoid
length) fits best with a threshold of 10-20 ADC on a single sample; the
simulation's default rule, 80 ADC held for 2 us, fits as well. How strict the real check is
cannot be told apart from these data.

Even the largest guard time the firmware allows (`BaselineGuardT` = 8 us, a
15 us freeze) would not explain run 97: with it the pulse still reads high
beyond about 18 us.

**How long a pulse is remembered.** `result/bin/blkernel` takes the pulses
with no other pulse within 40 us before or 10 us after, and measures how much
each earlier pulse on the same channel shifts their energy, by how long before
it came. The values are the shift per earlier pulse of beam-peak size, as a
fraction of the beam peak (run 97, 2.0 M pulses on the eight R long ends; the
L long ends, whose pole-zero error has the opposite sign, are the control):

| age [us] | 40-80 | 80-130 | 130-200 | 200-260 | 260-330 | 330-450 |
|---|---|---|---|---|---|---|
| run 97, R | +0.0030 | +0.0009 | +0.0005 | +0.0003 | +0.0004 | +0.0001 |
| run 97, L | -0.0008 | -0.0004 | -0.0001 | -0.0001 | +0.0001 | +0.0003 |
| averaging starts again right after the 7 us freeze | +0.0013 | +0.0000 | -0.0000 | +0.0001 | -0.0000 | -0.0000 |
| CAEN's description, 10 ADC | +0.0039 | +0.0050 | +0.0046 | +0.0043 | +0.0035 | +0.0034 |
| CAEN's description, 20 ADC | +0.0052 | +0.0059 | +0.0052 | +0.0049 | +0.0038 | +0.0037 |
| frozen until within 80 ADC for 2 us | +0.0014 | +0.0001 | -0.0001 | +0.0001 | +0.0000 | +0.0000 |
| same + 16 us smoothing (the default) | +0.0020 | +0.0007 | +0.0005 | +0.0002 | -0.0000 | +0.0002 |

(Simulated errors are about 0.0001 per bin. The 40-80 us bin includes the
previous pulse's own undershoot.)

So:

- In run 97 a pulse 80-500 us earlier still shifts the energy, by +0.09 % of
  the beam peak per pulse at 80-130 us, falling to +0.01 % by 330-450 us.
  Added over all earlier pulses at 44 kHz this is about a +1 % shift, which
  is the +1.4 % at 38-60 us in the first table.
- No model without the smoothing reproduces this. Waiting for the trapezoid
  to return gives nothing beyond 80 us. CAEN's version with one check per
  trapezoid length keeps the baseline frozen for hundreds of us and gives
  5-10 times too much. The 16 us smoothing roughly reproduces it, but CAEN
  says the firmware has no such filter. **The cause of this small effect is
  not known.** The L control shows a smaller effect (a quarter to a half of
  the size) with the opposite sign at 40-130 us, so it is at least partly
  tied to the R chain rather than to the beam rate. A slow second decay
  component in the R preamp is one possibility. It has not been measured: the
  R preamp fitted at run 74 has never been recorded with waveforms.
- Long averaging windows are ruled out: a 262 us window (Medium counted on a
  ÷32 clock), or any smoothing of 65 us or more, predicts four to ten times
  the run-97 effect at 80-330 us. The firmware agrees: the 2740 PHA has no
  decimation, and its largest window (High) is 131 us.
- Waiting for the trapezoid to return is not linear in pulse height: a bigger
  pulse leaves a bigger undershoot, which keeps the baseline frozen longer.
  The `PulseHistory` correction, which is linear in the earlier pulse heights,
  is fitted on pulses near the beam peak and is exact only for pulses of that
  size.
- On ch 4 in runs 25-29 (High window, a 10 us pole-zero on the ~19 us
  preamp, so an overshoot, at 27 kHz) the default model also fits the effect
  beyond 80 us (-0.0016 falling to -0.0002 at 450-650 us), while averaging
  again right after the freeze gives 5-10 times too much. Runs 13-24 (Medium
  on ch 4) show no effect beyond 80 us, where a ÷32 clock would leave a
  262 us plateau.

**Simulation defaults.** The default stays "frozen until within 80 ADC for
2 us, plus a 16 us smoothing", because it is the closest to run 97 in both
tables. The smoothing is there only to stand in for the unexplained effect
above. For the firmware as CAEN describes it, use
`--lfFilterOn=off --baselineRefreeze=on --baselineRestartThrADC=20` (add
`--rate_kHz=44 --tauR_us=11` for run 97). A baseline that waits for the
trapezoid can stay frozen longer than its window, so without the smoothing
the simulation rebuilds the baseline from 500 us before each pulse that
follows a quiet stretch (with only the usual ~15 us it released the freeze
too early). The older long-window defaults give the same energy-vs-gap
profile but remember a pulse far too long; results that depended on that have
been redone (tau recovery, kernel reach, applicability grid) or are marked as
long-window results.

**For the analysis** none of this changes the correction: `PulseHistory` fits
its kernel from the data (30 bins out to 316 us). The ~1 % overall shift is
the same for every pulse at a given rate and goes into the calibration.

