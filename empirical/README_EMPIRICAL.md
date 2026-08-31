# Empirical package (Ch. 7) — integration notes

Drop the contents of this folder into the repo root (`block_adaptive_blp/`).
Files land in the existing `data/`, `dgp/`, `estimators/`, `montecarlo/`
folders plus new `docs/` and `data/raw|derived`; there are no filename
collisions with the existing code.

## Run order

```
octave/matlab, from repo root, addpath(genpath(pwd)):
1. build_shock_series();                       % reads data/raw/ea_empd_events.csv (shipped)
2. fetch_outcome_data();                       % downloads or tells you what to download
3. assemble_dataset();                         % -> data/ea_dataset.mat  (baseline)
4. RUN_EMPIRICAL                               % steps A-D (edit toggles at top)
5. null = run_null_calibration(struct('null', struct('n_rep', 10)));   % TIMING RUN
6. null = run_null_calibration();              % real run (R=200; 500 for thesis)
```

## The one-line isrw patch (REQUIRED before step 4)

The empirical system mixes a white-noise-ish surprise (variable 1) with
persistent levels. `cfg.fmar.isrw` must accept a 1×K vector. In
`estimators/estimate_bvar_niw.m`, the prior-mean construction currently
uses the scalar flag (prior mean = isrw * eye(K) on the first own lag).
Replace that line with the diagonal form, e.g.

```matlab
% before:  b(2:K+1, :) = cfg.fmar.isrw * eye(K);      (or equivalent)
% after:
b(2:K+1, :) = diag(cfg.fmar.isrw(:));
```

(and analogously wherever `isrw` multiplies `eye(K)` — grep for `isrw`).
Scalar configs keep working if you wrap with
`if isscalar(cfg.fmar.isrw), cfg.fmar.isrw = cfg.fmar.isrw*ones(1,K); end`.
This was flagged in the FMAR-integration chat; it is the only repo edit
the package needs.

## Assumed repo signatures (verify once; all documented at call sites)

* `bvar = estimate_bvar_niw(Y, cfg)` with fields `.B` (m×K, constant row
  first, then p lag blocks), `.c`, `.A` (K×K×p), `.Sigma`, `.b1n`,
  `.theta` (K×(H+1)), `.lambda`, `.max_eig`.
  `simulate_fitted_bvar_dgp` asserts the B vs (c,A) layout numerically at
  run time, so a mismatch fails loudly, not silently.
* `blpf = estimate_blp_fmar(Y, cfg, bvar)` with `.lambda` (K×H).
* `blpb = estimate_blp_blockadaptive(Y, cfg, bvar, lambda_mat)` with
  `.theta_mean`, `.lo/.hi`, `.tau_mean` (K×G×H), G = K.
* utils: `build_lp_regressors(x, p)` layout `[1, y_t, y_{t-1}, ...]`;
  `empirical_quantile(x, probs)`; `normal_quantile(p)`;
  `safe_chol_lower(S)`; `default_config()`.
* `rng(seed, 'twister')` works on Octave ≥ 8 (as used in the repo tests).

If `estimate_blp_fmar`'s signature differs (e.g. it selects its own BVAR
internally), adapt the three call sites: `run_null_calibration.m` (one),
`RUN_EMPIRICAL.m` (two).

## What could not be verified from the build environment

* The four outcome-data URLs in `fetch_outcome_data.m` are written from
  portal documentation and flagged "verify"; the manual instructions in
  that file are the fallback and the parser only needs
  TIME_PERIOD/OBS_VALUE columns.
* Exact repo function signatures above (recovered from the development
  chats, not from the live repo).

## Shipped data

* `data/raw/ea_empd_events.csv` — extracted from your uploaded
  `EA-EMPD_en.xlsx` (sheet EA-EMPD; 4832 events 1999-01-07 to 2025-10-30;
  one stray `W` event dropped; columns: Date_time, Event_type,
  Days_until_next_GC, Non_regular_trading_day,
  Outside_regular_trading_hours, OIS_1M/2M/3M/6M/1Y/2Y, STOXX50E).
  NOTE: the workbook's Notes sheet describes `Outside_regular_trading_hours`
  with inverted wording; verified against timestamps, 1 = outside
  9:00–18:00 CET. `build_shock_series.m` relies on the verified semantics.

## Reference sanity numbers (2001m1–2019m12; printed by build_shock_series)

GC_ME events 221; mps_gc_1y std 4.29 bp, 23 zero months; JK keeps 114/221
(51.6%), std 3.52 bp; speeches kept 1596/2701; mps_all_1y std 5.12 bp,
corr 0.856 with GC-only; mps_all_tgt std 4.90 bp. A quarterly EA-EMPD
update shifts these slightly — the script warns, it does not error.
