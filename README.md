# Block-Adaptive Bayesian Local Projections

Companion code for the MSc thesis *"Where Does the VAR Prior Fail?
Block-Adaptive Bayesian Local Projections under Sparse Dynamic
Misspecification."*

Bayesian local projections (BLPs) shrink local-projection coefficients
towards the values an estimated VAR implies. This project gives every
coefficient *block* (the lags of one variable, in one equation, at one
horizon) its own escape scale `tau` with a grouped-horseshoe prior, so that
the part of the system where the VAR prior is wrong can leave it while the
rest stays disciplined. `tau = 1` recovers the published estimator of
Ferreira, Miranda-Agrippino and Ricco (REStat 2025) exactly. The repository
contains the estimators, a Monte Carlo study of when block adaptation
helps, and a euro-area monetary-policy application in which the posterior
`tau` map is read as a *calibrated specification diagnostic* against a null
simulated for the design at hand.

## Where to start

| you want to | read |
|---|---|
| understand the project in plain language, file by file | `docs/PROJECT_EXPLAINED.md` |
| the estimators, the fair-comparison conventions, the diagnostic and its null | `docs/METHOD.md` |
| the simulation evidence and what it does and does not support | `docs/RESULTS.md` |
| the euro-area application: design, data, results | `empirical/docs/DESIGN.md`, `empirical/docs/DATA.md`, `empirical/docs/RESULTS.md` |
| a critical assessment of the empirical chapter and what to do next | `docs/ASSESSMENT.md` |
| every command, with runtimes | `docs/REPRODUCE.md` |
| what each remaining approximation costs | `docs/APPROXIMATIONS.md` |
| how the project got here | `docs/CHANGELOG.md` |

## Status (2026-09-15)

* **Simulations (settled).** On four designs at `R = 500`/`250`, the
  adaptive estimators lower RMSE at early horizons and raise it at late
  ones; horizon pooling keeps the early gain and makes the integrated RMSE
  indistinguishable from the FMAR baseline. Localisation of a sparse
  prior–data conflict works (97% against a 37% null) and is honest when
  there is nothing to localise; the null is specification-dependent, and
  finding where the prior fails does not imply that escaping helps. The
  thesis is framed around the calibrated diagnostic, not around an
  improved estimator (`docs/RESULTS.md` §5).
* **Empirical chapter (current design).** Four-variable euro-area system
  (end-of-month 1-month OIS, industrial production, SA HICP, EURO STOXX 50;
  2000–2019), identified with the EA-EMPD 1-month OIS surprise as an
  external instrument. First-stage effective F 13.7 (LP-IV 19.4), with two
  thirds of the identifying covariance in 2009–2011. The `tau` map, read
  against its own nulls with Monte Carlo and family-wise p-values, flags one
  cell at `p = 12`, the short-rate block of the output equation, under both
  the iid-resampling and the heteroskedasticity-robust null; the escape is
  absent at `p = 4` and `6`, appears at `p = 2` only under the homoskedastic
  null, and carries no out-of-sample forecast content
  (`empirical/docs/RESULTS.md`).
* **Assessment.** `docs/ASSESSMENT.md` lists the weaknesses of the
  empirical application, what was investigated on 2026-09-15, and the
  prioritised next steps.

## Requirements

MATLAB R2018b or newer, **no toolboxes** (gamma / inverse-gamma / inverse-
Wishart draws, quantiles and normal quantiles are implemented in `utils/`).
Also runs under GNU Octave ≥ 8 (the simulation results were produced under
Octave 8.4; the test suite and the empirical chapter ran under MATLAB
R2025b). Parallelism, where used, is separate operating-system processes
(`scripts/`), never a toolbox. The euro-area outcome and indicator series
are third-party data and are not in the repository; `empirical/docs/DATA.md`
says how to obtain them.

## Quick start

```matlab
cd tests; run_all_tests           % 21 tests, ~12 min in MATLAB
RUN_FMAR_DEMO                     % one sparse dataset through the FMAR stack + exact nesting check
QUICK_DEMO = true; RUN_FMAR_DEMO  % ~2 minutes
report_montecarlo('results/simulation/mc_headline_sparse.mat')   % the headline tables
```

Empirical chapter (after the data are in place, `empirical/README.md`):

```matlab
addpath(genpath(pwd))
assemble_dataset_v2(struct('system', 'ois4'));
SYSTEM = 'ois4'; run(fullfile(ea_paths().empirical, 'RUN_EMPIRICAL_IV.m'))
```

## Repository layout

```
RUN_ME_FIRST.m, RUN_FMAR_DEMO.m   single-dataset demonstrations (prototype mode / FMAR mode)
config/        default_config.m: every setting (cfg.mode, cfg.fmar, cfg.blp.pool, cfg.dgp, ...)
dgp/           the four simulation designs + analytic true IRFs
estimators/    estimate_bvar_niw (GLP Bayesian VAR), estimate_blp_fmar (global baseline),
               estimate_blp_blockadaptive (independent tau), estimate_blp_blockpooled (pooled tau),
               estimate_lp / estimate_var / estimate_blp_global (prototype mode)
priors/        prior centre, Newey-West scales, marginal-likelihood tightness, block prior
samplers/      gibbs_block_horseshoe (one equation-horizon), gibbs_block_pooled_horizons (one equation, all h)
montecarlo/    run_montecarlo, summarize/report/export, tau_diagnostics, paired_comparison,
               presets, exact parallel chunking, the exploratory grid, sensitivity analysis
scripts/       shell drivers for the parallel Monte Carlo runs
plots/, utils/ figures; draws, quantiles, Cholesky, regressor builder, config save/load
tests/         assertion-based tests (run_all_tests)
empirical/     Chapter 7 (self-contained; see empirical/README.md)
  RUN_EMPIRICAL_IV.m          the driver of the current (external-instrument) design
  ea_paths, ea_design_key, ea_apply_protocol, ea_cross_p_table
  data/                       instrument construction, data fetchers, dataset assembly, identification, relevance
  estimators/ dgp/ montecarlo/ LP-IV, lag-augmented LP; the bootstrap null; null calibration; block ablation
  legacy_v1/                  the 2026-09-12 internal-instrument design, kept runnable
  tests/                      synthetic fixtures (refused by every driver) and package tests
  docs/                       DESIGN, DESIGN_V1_LEGACY, LITERATURE, DATA, RESULTS, RESULTS_LOG
docs/          PROJECT_EXPLAINED, METHOD, RESULTS, APPROXIMATIONS, REPRODUCE, ASSESSMENT, CHANGELOG
results/       simulation/  empirical/  empirical/legacy_v1/   (see results/README.md)
presentation/  status decks (LaTeX)
```

## Conventions

* Every entry point resolves its own absolute paths (`ea_paths.m` for the
  empirical package), so any script runs from any working directory. Under
  `matlab -batch`, call scripts by absolute path after `addpath(genpath(<repo>))`.
* Outputs never overwrite an earlier design: simulation runs carry a preset
  in their name, empirical outputs a system and lag-order tag, nulls a stem
  and the residual scheme (`_wild`), and a null is refused for a design
  whose key differs from its own.
* The empirical drivers refuse any dataset stamped `synthetic = true`. The
  fixtures exist only so the package can be tested on a machine without
  the licensed data; no empirical number may come from them.
* Which number to quote: the integrated RMSE over `h = 2..H` for the
  simulations; never a raw `tau`, always its null-relative ratio and
  p-values with the null named (`results/README.md`).
