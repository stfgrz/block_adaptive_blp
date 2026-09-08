Output directory for figures, .mat files and CSV exports (simulation AND
empirical).

WHAT LANDS HERE
---------------
mc_final_<dgp>.mat        headline Monte Carlo runs (R = 500 for the
                          sparse and correct DGPs, R = 250 for the
                          intermediate and dense ones).  Each contains
                          `mc` (raw replication output) and `s`
                          (summarize_montecarlo).  Read them with
                              report_montecarlo('results/mc_final_sparse.mat')
mc_final_<dgp>_*.csv      the same numbers as tidy CSV:
   *_metrics.csv          bias, variance, MSE, RMSE, coverage and
                          interval length per (estimator, response,
                          horizon), plus the posterior-quantile band
                          columns for the sampled estimators
   *_integrated.csv       both integrated-RMSE definitions (legacy
                          h = 1..H and the PREFERRED h = 2..H), the
                          early/late splits, integrated bias and variance
   *_tau.csv              posterior mean/median and P(tau > 1) for every
                          (equation, block, horizon)
   *_tau_localization.csv localisation, false-positive and chain
                          diagnostics
   *_meta.txt             the reproducibility record: DGP parameters and
                          seeds, R/T/p/H, estimator settings, sampler
                          sizes, preset name
grid/                     the exploratory experiment grid: one .mat and
                          one CSV set per cell, plus grid_summary.csv
                          with the per-cell verdict
chunks/                   intermediate per-process chunks of a parallel
                          run, and their logs.  Not tracked by git; the
                          merged mc_final_* files are what matters
sensitivity_approximations.csv
                          what the three remaining approximations cost
                          (see docs/APPROXIMATIONS.md)
fig*.png                  figures from RUN_ME_FIRST / RUN_FMAR_DEMO

EARLIER RESULTS
---------------
mc_*.mat files without the `final` prefix are from earlier runs and are
kept as-is.  They predate three changes and are NOT comparable
term-by-term with the new ones:
  * they report the legacy integrated RMSE over h = 1..H under the FMAR
    h = 1 convention (see README section 3a);
  * the adaptive point estimate was the average of the IRF draws rather
    than the Rao-Blackwellised posterior mean;
  * they have four estimators, not five.
Use mc_preset('legacy') to reproduce that configuration.

WHICH NUMBER TO QUOTE
---------------------
The integrated RMSE over h = 2..H (`s.irmse_h2`, column `irmse_h2` in
*_integrated.csv).  The h = 1..H version is kept for compatibility and
is not a fair estimator comparison under the FMAR h = 1 convention.
