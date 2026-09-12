# Empirical package (Ch. 7) — how to run it, and what it found so far

The package lives in its own `empirical/` subtree. It is **not** unpacked
over the repo root: the root already has `dgp/` and `estimators/` folders
holding the *simulation* machinery, and mixing the euro-area
data-construction code into them makes both harder to read. Results still
land in the repo-wide `results/`, next to the Monte Carlo output.

Every path inside the package is absolute, resolved from the package's own
location by `ea_paths.m`, and every entry point puts the repo on the path
itself. **You can run any of these from any working directory** — no
`cd`, no `addpath(genpath(pwd))` required.

```
empirical/
  ea_paths.m                 absolute locations (single source of truth)
  RUN_EMPIRICAL.m            the driver: steps A-F
  SMOKE_TEST_EMPIRICAL.m     two-second interface check (run before the long job)
  README_EMPIRICAL.md        this file
  docs/CH7_DESIGN.md         the design document (pre-registered expectations)
  data/
    build_shock_series.m     EA-EMPD events -> monthly surprise series
    fetch_outcome_data.m     download/verify the four outcome series
    assemble_dataset.m       -> data/ea_dataset.mat
    read_sdmx_csv.m          minimal SDMX-CSV reader
    ea_check_series.m        input validation (see "Data integrity" below)
    ea_extract_series.m      reduce an over-broad Eurostat export to 1 series
    ea_write_provenance.m    <file>.source.txt sidecars
    ea_relevance.m           instrument-relevance diagnostics
    raw/ea_empd_events.csv   shipped: 4832 events, 1999-01-07..2025-10-30
    raw/placeholder_rejected/  quarantined fake data — read its README
  dgp/simulate_fitted_bvar_dgp.m      the null bootstrap
  estimators/estimate_lp_lagaug.m     MOP lag-augmented LP benchmark
  montecarlo/run_null_calibration.m   Leg 2, the long job
  tests/make_synthetic_fixture.m      SIMULATED dataset for interface tests
  tests/test_empirical_pipeline.m     the package on that fixture (in run_all_tests)
```

## Status (2026-09-12)

**The data are in place on this machine.** The four outcome series were
downloaded on 2026-09-01 (three through the portal APIs, the
industrial-production file as a Data Browser export that had every
geography and three NACE aggregates stacked in one file; it was reduced in
place to `M.PRD.B-D.SCA.I21.EA20` with `ea_extract_series`, and the
original export is kept as `raw/ip_ea_full_export.csv`). All four pass
`ea_check_series`; their values match known history (12-month Euribor
5.38% in September 2008 and −0.26% in December 2019, the STOXX 50 trough
of 1994 in March 2009, industrial production 98.6 → 71.2 between February
and April 2020). Provenance sidecars name the series keys read from the
files' own headers. The raw files and the assembled datasets are
gitignored (third-party data), so a clean clone still has to run
`fetch_outcome_data`.

**The whole pipeline has run**, in MATLAB R2025b: `assemble_dataset`
(baseline and JK variants), `SMOKE_TEST_EMPIRICAL`, `RUN_EMPIRICAL`
steps A–F (about four minutes at `p = 12`, `H = 48`), and the null
calibrations (`run_null_calibration`, about 8 seconds per replication at `p = 12`).
The results and their reading are in `docs/CH7_RESULTS.md`; the headline
is that, read against a null calibrated on the fitted BVAR (`R = 500`),
two cells survive at the 5% level — the HICP own-lag block and the
interest-rate block in the industrial-production equation — while the
rest of the map is within what a correctly specified VAR(12) produces on
data like these. This file documents the three things the real data
forced us to decide, because they change the design as written on
2026-09-01.

### 1. The surprise is a weak instrument for the monthly-average 1-year rate

`ea_relevance` regresses the 1-year-rate innovation on the surprise
innovation, both from the BVAR residuals (that slope *is* the impact
response the 25 bp normalisation divides by). At `p = 12` the impact of a
1 pp surprise innovation is 0.24 pp with a robust t of 0.77 (F = 0.59);
the naive same-month regression of the rate change on the surprise gives
0.03 (t = 0.07). The reasons are visible in the data: the indicator is a
monthly *average*, so a mid-month event is spread over two months; and
the largest surprises of the sample sit in late 2008, when a positive
surprise (a smaller cut than expected) coincided with 90 bp monthly falls
in the rate. The design document budgeted an F of about 4 (AGKL, 1-month
OIS on both sides), so this is weaker than pre-registered.

What it does and does not affect. The τ diagnostic asks whether the
horizon-`h` LP coefficients of each block disagree with the VAR-implied
centre; that is a statement about the dynamic specification of the
five-variable system and does not involve the instrument's strength. The
*structural* reading of the IRFs as responses to a monetary policy shock,
and the 25 bp scale factor (`k = 1.06` at `p = 12`, `3.5` at `p = 2`), do
depend on it and must be presented with that caveat. Every run prints
these numbers and stores them (`rel` in `results/empirical_<tag>.mat`).

### 2. The published prior scale breaks on a white-noise variable (fixed by a floor)

FMAR's prior scale `psi_v(h)` is the Newey-West long-run variance of a
univariate `h`-step LP residual. For the surprise, a mean-reverting
near-white-noise series, that estimate is *one fifth of the residual
variance* at `h ≥ 20` (the Bartlett sum adds negative autocovariances).
Since the prior variance of every coefficient attached to lagged
surprises is `Sigma_ii lambda_h^2 / psi_mps(h)`, those 12 coefficients per
equation get an absurdly loose prior, and the marginal likelihood
compensates by driving the **global** `lambda_h` down to 0.01–0.03 for
`h ≤ 28`; from `h ≈ 27` the objective becomes bimodal (modes near 0.01 and
0.4 of almost equal height), the global search flips between them across
adjacent horizons, and the τ surface shows a vertical stripe at
`h ≈ 30` with the FMAR IRF jumping onto the LP there. This was
established with controlled variants of the objective (recorded in the
session log): flooring `psi` at the plain residual variance makes the
objective unimodal at every horizon and the `lambda_h` path smooth
(0.008 → 0.11 at `h = 12` → 0.4–0.5 at `h ≥ 24`); removing the surprise
block from the selection reproduces the four-variable level system's own
path (0.01 → 0.36 → 0.9); and a joint search wants the surprise-lag block
a thousand times tighter than the Minnesota scaling gives it.

The remedy adopted is the floor, `cfg.fmar.psi_floor = true`
(`priors/fmar_prior_scale.m`): a long-run variance is never allowed below
the variance it is a long-run version of. It never binds for the four
persistent series, and on the simulation designs it binds by a few
percent at `h = 1` only (largest IRF change 1.5e-4;
`tests/test_psi_floor.m`). The default remains `false`, so every stored
simulation result is what it was. `select_lambda_fmar` now counts the
local maxima of its objective and `estimate_blp_fmar` reports how many
horizons were multimodal (`blp.diag.lambda_multimodal`): 0 of 48 with the
floor, 28 of 48 without it at `p = 12`.

Even with the floor, the surprise block still drags `lambda_h` down
relative to the level system alone (0.11 against 0.36 at `h = 12`). The
null calibration absorbs this (it simulates the same pipeline), but the
*levels* of τ in the five-variable heatmap are inflated relative to the
four-variable companion exhibit (step F). The principled next step is
AGKL's own treatment of the instrument — zero restrictions on the
surprise block, i.e. lagged surprises out of the horizon regressions —
which removes the offending coefficients rather than rescaling them.

### 3. The instrument block cannot be a horseshoe block (held at τ = 1)

With every block free, the surprise block's τ collapsed to 0.03–0.4 in
every equation: eleven near-zero lag coefficients dominate the group
scale and drag the one coefficient that carries the identification — the
contemporaneous surprise — to the VAR centre. The adaptive IRFs then sat
on the BVAR at long horizons while the global BLP followed the LP, and
the coherence check of Sec. 5.3 failed for that reason alone (block closer
to the LP than FMAR in 31% of cells, no relation to τ). The block prior
assumes a common scale within a block; the instrument block violates it
by construction. It is therefore **held at τ = 1**
(`cfg.blocks.fixed_tau = 1`; `tests/test_fixed_tau_blocks.m`): it stays
under the global FMAR prior like the intercept, and the diagnostic is
read on the four level blocks of each equation. The surprise *equation*
remains the negative control (its four level blocks should be quiet).

### Output tags

Every output name carries the design: `p12` is the baseline (floor on,
instrument block held), `p12_nofloor` the published scale with every
block free (the literal FMAR-nested estimator), `p12_allblocks` the floor
with every block free. The level-system files (step F) are
`*_levels_p12*`. Null calibrations follow the same convention
(`null_calibration_p12.mat`, `..._p2.mat`, `..._p6.mat`,
`..._p12_levels.mat`).

## Run order

```matlab
build_shock_series();      % seconds. Check the printed stats against the
                           % reference numbers; all should match exactly.

fetch_outcome_data();      % must end with "all four outcome files ready
                           % and validated". If a download fails it prints
                           % the exact manual route for that series; an
                           % over-broad Eurostat export is reduced with
                           % ea_extract_series (see Data integrity).

assemble_dataset();        % seconds. T = 240, no NaN. Read the printed
                           % first stage: it is WEAK here (see Status 1).

SMOKE_TEST_EMPIRICAL       % two seconds. Verifies every estimator interface
                           % the long job depends on, on the real dataset.

RUN_EMPIRICAL              % steps A-F, about four minutes in MATLAB.
```

Then the null calibrations (each a separate job; the baseline one takes
about 85 minutes at `R = 500`):

```matlab
null = run_null_calibration(struct('null', struct('n_rep', 10)));   % TIMING
null = run_null_calibration(struct('null', struct('n_rep', 500)));  % baseline p = 12
null = run_null_calibration(struct('p', 2, 'null', struct('n_rep', 200)));   % dose-response
null = run_null_calibration(struct('p', 6, 'null', struct('n_rep', 200)));
null = run_null_calibration(struct('null', struct('n_rep', 200, ...           % level system
           'dataset', fullfile(ea_paths().data, 'ea_dataset_levels.mat'), ...
           'fixed_tau', [], 'stem', 'p12_levels')));
```

and finally the protocol tables:

```matlab
DO_A = false; DO_B = false; DO_D = false; RUN_EMPIRICAL   % steps C, E, F re-read the nulls
```

## Data integrity — read this before touching `ea_check_series`

`fetch_outcome_data` and `assemble_dataset` both validate the four outcome
series and **refuse** to proceed on data that fails. This is not
defensive boilerplate: the package originally shipped four *placeholder*
csvs under the real filenames, with the right columns and no gaps, whose
contents were generated rather than observed — HICP rose by exactly +0.13
index points every month for 323 months, the 1-year rate sat at exactly
−0.5000 for years, and neither industrial production nor the STOXX 50
showed the 2008 crisis or the 2020 collapse. The old `fetch_outcome_data`
reported *"all four outcome files ready"* for them. They are now in
`data/raw/placeholder_rejected/`; that folder's README documents them.

The general tripwire is the **roughness ratio** `sd(diff²y)/sd(diff y)`,
which is ≈1.4 for any real monthly macro series (exactly √2 for a random
walk) and was 0.00–0.18 for all four placeholders. Two related
fingerprints (modal-step share, longest flat run) catch clipped and
constant-increment series. Thresholds and reasoning are in the header of
`data/ea_check_series.m`. On the real files the ratio is 1.38 (HICP),
1.44 (IP), 1.24 (STOXX) and 0.78 (Euribor, a smooth monthly average).
Override with `opts.skip_plausibility = true` only after looking at the
series and deciding a fingerprint is a false alarm.

**On the URLs.** The three ECB/HICP keys are confirmed working. The
Eurostat industrial-production key was corrected after a live run: the
indicator code is **`PRD`** (not `PROD`) and the current index base is
**`I21`** (2021 = 100), giving
`sts_inpr_m/M.PRD.B-D.SCA.I21.EA20` — the one series the chapter uses
(the silent EA19 fallback is gone: a different geography). A different
index base year is harmless: the variables enter as `100*log(index)`, so
a rescaling is absorbed by the constant.

**Over-broad exports.** The Eurostat Data Browser exports every geo and
activity if you leave those filters open — one file with dozens of
stacked series, which parses but is not a series. `ea_check_series`
catches it (duplicated `TIME_PERIOD`). Rather than re-downloading:

```matlab
ea_extract_series('empirical/data/raw/ip_ea.csv');          % list dimensions
ea_extract_series('empirical/data/raw/ip_ea.csv', ...       % reduce in place
                  'empirical/data/raw/ip_ea.csv', ...
                  struct('geo', 'EA20', 'nace_r2', 'B-D'));
```

**Provenance.** Every accepted raw file gets a `<file>.source.txt`
sidecar naming the series and where it came from, and `assemble_dataset`
copies those lines into `ds.meta.sources`. Files placed by hand are
recorded as such.

## Repairs made along the way

* `assemble_dataset.m` once did not parse (a missing `end`); fixed on
  2026-09-08 and pinned by `tests/test_empirical_pipeline.m`.
* `SMOKE_TEST_EMPIRICAL.m` asserted `all()` of a *matrix*, which Octave
  reduces silently and MATLAB rejects; the first MATLAB run of the
  package stopped there. Fixed (`all(all(...))`); the rest of the repo
  was scanned for the same hazard and has none.
* The vectorised `cfg.fmar.isrw` prior mean (2026-09-01): the empirical
  system needs `isrw = [0 1 1 1 1]`, and `if cfg.fmar.isrw` on a vector is
  `all()` of it. `estimate_bvar_niw` builds the prior mean as
  `diag(isrw)`; scalars are bit-identical to before
  (`tests/test_isrw_vector.m`).
* The pooled estimator (`estimate_blp_blockpooled`) is now run alongside
  the independent one in every step and in the null calibration.

## Repo signatures (verified, not assumed)

`SMOKE_TEST_EMPIRICAL` asserts all of these on the real dataset, so they
are checked rather than documented:

* `bvar = estimate_bvar_niw(Y, cfg)` → `.B` (m×K, constant row first,
  then p lag blocks), `.c`, `.A` (K×K×p), `.Sigma`, `.F`, `.Psi`, `.b1n`,
  `.theta` (K×(H+1)), `.theta_lo/.theta_hi`, `.lambda`, `.max_eig`.
* `blpf = estimate_blp_fmar(Y, cfg, bvar)` → `.lambda` (K×H),
  `.theta_mean`, `.diag.lambda_multimodal`, `.diag.psi_floor_bound`.
* `blpb = estimate_blp_blockadaptive(Y, cfg, bvar, lambda_mat)` →
  `.theta_mean`, `.lo/.hi` (K×(H+1)), `.tau_mean` (K×G×H, G = K, index
  h = 1..H), `.tau_q`, `.p_tau_gt1`, `.theta_cond`, `.beta_block`,
  `.beta_mean` (m×K×H, the full posterior-mean coefficient vectors, which
  the coefficient-level coherence check of step E uses; `estimate_blp_fmar`
  returns the same field), `.fixed_tau_blocks`.
* `blpp = estimate_blp_blockpooled(Y, cfg, bvar, lambda_mat)` → the same
  interface plus `.kappa_mean`.
* `lp = estimate_lp_lagaug(Y, cfg, bvar)` → `.theta` (K×(H+1)).
* `rel = ea_relevance(Y, bvar, cfg)` → `.impact_b/.impact_t/.impact_F`,
  `.naive_b/.naive_t`, `.k25`.
* `rng(seed, 'twister')` works on Octave ≥ 8.

## Notes on the two drivers

* **`RUN_EMPIRICAL`** — every step saves its own `.mat`, so a run that
  dies part-way can be resumed by flipping the completed steps' `DO_`
  toggles to `false`; every toggle (and `PSI_FLOOR`, `INSTR_FIXED`,
  `DO_POOLED`) can be preset in the workspace before running. Steps B, D
  and E read step A's output and stop with an actionable message if it is
  not there; D refuses to compare against a step-A file written at a
  different `p` or `H`. Figures are wrapped in try/catch. The saved BVAR
  structs omit the per-draw companion matrices (14 MB at `p = 12`), which
  only the simulation-side sensitivity analysis uses.
* **`run_null_calibration`** — defaults to `p = 12`, `h1_mode = 'lp'`,
  `psi_floor = true`, `fixed_tau = 1`, matching `RUN_EMPIRICAL`. Thresholds
  calibrated under a different setting do not describe the null of the
  statistic they are read against, so the design is recorded in the
  output and in the checkpoint, and a checkpoint is resumed only when it
  matches. It records both estimators, the null distribution of
  `lambda_h` and how often the objective was multimodal.

The early-horizon window of the `tau_bar` statistic is `[2, min(12, H)]`
in both drivers, so they always compute the same statistic.

## Shipped data

`data/raw/ea_empd_events.csv` — extracted from the EA-EMPD workbook
(sheet EA-EMPD; 4832 events 1999-01-07 to 2025-10-30; one stray `W` event
dropped). Columns: `Date_time`, `Event_type`, `Days_until_next_GC`,
`Non_regular_trading_day`, `Outside_regular_trading_hours`,
`OIS_1M/2M/3M/6M/1Y/2Y`, `STOXX50E`.

The workbook's Notes sheet describes `Outside_regular_trading_hours` with
inverted wording; the column *name* is correct (1 = outside 09:00–18:00
CET), verified against the timestamps. `build_shock_series.m` relies on
the verified semantics.

## Reference numbers (2001m1–2019m12, printed by `build_shock_series`)

These match the shipped events file **exactly** (re-verified in MATLAB on
2026-09-12):

| statistic | value |
|---|---|
| GC_ME events | 221 |
| `mps_gc_1y` monthly std | 4.29 bp |
| zero months | 24 |
| JK kept | 115 / 221 |
| `mps_gc_jk` std | 3.51 bp |
| speeches kept | 1596 |
| `mps_all_1y` std | 5.12 bp |
| corr with `mps_gc_1y` | 0.856 |
| `mps_all_tgt` std | 4.90 bp |

The JK count is 113 events with a strictly opposite-sign STOXX move plus
2 with an exactly zero OIS_1Y change, which the documented tie rule
assigns to the policy series. A quarterly EA-EMPD update will shift these
slightly — the script warns, it does not error.
