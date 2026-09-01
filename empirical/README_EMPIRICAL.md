# Empirical package (Ch. 7) — how to run it

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
  RUN_EMPIRICAL.m            the driver: steps A-D
  SMOKE_TEST_EMPIRICAL.m     5-minute interface check (run before the long job)
  README_EMPIRICAL.md        this file
  docs/CH7_DESIGN.md         the design document
  data/
    build_shock_series.m     EA-EMPD events -> monthly surprise series
    fetch_outcome_data.m     download/verify the four outcome series
    assemble_dataset.m       -> data/ea_dataset.mat
    read_sdmx_csv.m          minimal SDMX-CSV reader
    ea_check_series.m        input validation (see "Data integrity" below)
    raw/ea_empd_events.csv   shipped: 4832 events, 1999-01-07..2025-10-30
    raw/placeholder_rejected/  quarantined fake data — read its README
  dgp/simulate_fitted_bvar_dgp.m      the null bootstrap
  estimators/estimate_lp_lagaug.m     MOP lag-augmented LP benchmark
  montecarlo/run_null_calibration.m   Leg 2, the long job
```

## Run order

```matlab
build_shock_series();      % seconds. Check the printed stats against the
                           % reference numbers; all should match exactly.

fetch_outcome_data();      % must end with "all four outcome files ready
                           % and validated". If a download fails it prints
                           % the exact manual route for that series; do the
                           % download, then rerun until it says ready.

assemble_dataset();        % seconds. Expect T = 240, no NaN, and a
                           % weak-ish first stage (meetings-only surprises).

SMOKE_TEST_EMPIRICAL       % ~1 minute. Verifies every estimator interface
                           % the long job depends on. If this passes,
                           % RUN_EMPIRICAL will not fail on a signature.

RUN_EMPIRICAL              % steps A-D, a few hours at p = 12, H = 48.
```

Then, as a **separate** job (do not queue it the same night):

```matlab
null = run_null_calibration(struct('null', struct('n_rep', 10)));  % TIMING
null = run_null_calibration();                                     % R = 200
```

Do the `n_rep = 10` timing run first and scale `n_rep` from the printed
seconds-per-rep.

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
`data/ea_check_series.m`. Override with `opts.skip_plausibility = true`
only after looking at the series and deciding a fingerprint is a false
alarm.

**The four outcome-series URLs are unverified.** They were written from
portal documentation and could not be tested from the environment this
package was prepared in (`ec.europa.eu` and `data-api.ecb.europa.eu` are
both unreachable from it). Several candidates are tried per series,
including the 2021 = 100 base Eurostat has been migrating STS datasets to.
A different index base year is harmless: the variables enter as
`100*log(index)`, so a rescaling is absorbed by the constant.

## The vectorised `isrw` prior mean (already applied)

The empirical system mixes a white-noise-ish surprise with persistent
levels, so `cfg.fmar.isrw` must be the vector `[0 1 1 1 1]`. The prior
mean in `estimators/estimate_bvar_niw.m` used to read the flag as
`if cfg.fmar.isrw`, which in MATLAB/Octave is `all()` of the elements: a
vector starting with 0 evaluates **false**, silently dropping the
random-walk centre for *every* variable rather than just the first. No
error, just a wrong prior.

The prior mean is now built as `b(2:K+1,:) = diag(isrw)`, with a scalar
expanded to `K` entries. **Scalar configurations are bit-identical to
before** (`false` → zeros, `true` → `eye(K)`); this was verified field by
field on the BVAR, the FMAR BLP and the block-adaptive posterior before
and after the change, and is pinned by `tests/test_isrw_vector.m`, which
is part of `run_all_tests`. This is the only edit the package needed
outside `empirical/`.

## Repo signatures (verified, not assumed)

`SMOKE_TEST_EMPIRICAL` asserts all of these on the real dataset, so they
are checked rather than documented:

* `bvar = estimate_bvar_niw(Y, cfg)` → `.B` (m×K, constant row first,
  then p lag blocks), `.c`, `.A` (K×K×p), `.Sigma`, `.F`, `.Psi`, `.b1n`,
  `.theta` (K×(H+1)), `.theta_lo/.theta_hi`, `.lambda`, `.max_eig`.
  `simulate_fitted_bvar_dgp` also asserts the `B` vs `(c, A)` layout
  numerically at run time, so a mismatch fails loudly.
* `blpf = estimate_blp_fmar(Y, cfg, bvar)` → `.lambda` (K×H),
  `.theta_mean`.
* `blpb = estimate_blp_blockadaptive(Y, cfg, bvar, lambda_mat)` →
  `.theta_mean`, `.lo/.hi` (K×(H+1)), `.tau_mean` (K×G×H, G = K, index
  h = 1..H).
* `lp = estimate_lp_lagaug(Y, cfg, bvar)` → `.theta` (K×(H+1)).
* utils: `build_lp_regressors(x, p)` with layout `[1, y_t, y_{t-1}, ...]`;
  `empirical_quantile(x, probs)`; `normal_quantile(p)`;
  `safe_chol_lower(S)`; `default_config()`.
* `rng(seed, 'twister')` works on Octave ≥ 8.

## Notes on the two drivers

* **`RUN_EMPIRICAL`** — every step saves its own `.mat`, so a run that
  dies part-way can be resumed by flipping the completed steps' `DO_`
  toggles to `false`. Steps B and D read step A's output and stop with an
  actionable message (not a cryptic load error) if it is not there, and D
  additionally refuses to compare against a step-A file written at a
  different `p` or `H`. Output filenames carry the lag order
  (`empirical_p12.mat`, `tau_heatmap_p12.csv`), so a run at another `p`
  cannot overwrite the baseline the write-up cites. Figures are wrapped in
  try/catch: on a headless machine you lose the pngs and keep every `.mat`
  and `.csv`, which are the real outputs.
* **`run_null_calibration`** — defaults to **`p = 12`**, matching
  `RUN_EMPIRICAL` step A. It used to inherit `default_config`'s `p = 2`,
  which belongs to the simulation DGPs; thresholds calibrated at the wrong
  lag order do not describe the null distribution of the statistic they
  are compared against, and nothing in the output would have said so.
  Override with `run_null_calibration(struct('p', 6))`. The checkpoint is
  only resumed when its stored `n_rep`, `p` and `H` all match.

The early-horizon window of the `tau_bar` statistic is `[2, min(12, H)]`
in both drivers, so they always compute the same statistic and a
reduced-`H` run cannot die on an out-of-bound index.

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

These match the shipped events file **exactly**:

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
