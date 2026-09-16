# The project, explained as simply as possible

This note is the map you read when you come back to the project after a
break. It explains, in order: the question, the estimation procedure for
one dataset, the simulation study, the euro-area application, and finally
every file in the repository in one or two lines. Numbers are quoted only
where they change how you read something; the authoritative tables are in
`results/simulation/REPORT.txt` and `docs/RESULTS.md`.

---

## 1. The question in one paragraph

An impulse response function (IRF) says how variable `i` reacts `h` months
after a shock. Two textbook ways to estimate it:

* a **VAR** fits one dynamic system and iterates it forward. Efficient, but
  if the system is wrong anywhere, every horizon inherits the error;
* a **local projection (LP)** runs a separate regression of `y(t+h)` on
  today's information for every `h`. Robust, but noisy.

A **Bayesian local projection (BLP)** runs the LP regressions but shrinks
their coefficients toward the values a fitted VAR implies. The published
version, Ferreira, Miranda-Agrippino and Ricco (REStat 2025, "FMAR"),
uses one tightness `lambda_h` per horizon: every coefficient is pulled to
the VAR with the same force.

This project asks: **if the VAR is wrong in only one part of the system,
can we let exactly that part of the LP escape the VAR prior, and does that
help?** The device is a block-specific scale `tau_{i,g,h}` on the
deviation of each coefficient block from the VAR centre, with a half-Cauchy
(group horseshoe) prior. `tau = 1` recovers FMAR exactly; `tau > 1` is a
local escape; `tau < 1` is local tightening.

The answer the simulations gave (Section 3) is: as an *estimator* the
device is roughly RMSE-neutral once `tau` is smoothed across horizons; as a
*diagnostic* of where the VAR prior and the data disagree it works and can
be calibrated. The thesis is therefore framed around the diagnostic, with
the euro-area application (Section 4) as its demonstration.

---

## 2. The procedure for one dataset, step by step

Notation: `K` variables, `p` lags, `T` observations, horizons `h = 0..H`.
The regressor vector is always `z(t) = [1; y(t); y(t-1); ...; y(t-p+1)]`
(`m = 1 + K p` entries), built once by `utils/build_lp_regressors.m` so
that every estimator conditions on the same information.

**Step 1: fit the Bayesian VAR** (`estimators/estimate_bvar_niw.m`).
A VAR(p) with the Giannone-Lenza-Primiceri Minnesota / Normal-Inverse-
Wishart prior. The prior mean of each variable's first own lag is 1 for
variables flagged as random walks (`cfg.fmar.isrw`) and 0 otherwise. The
overall tightness `lambda` maximises the closed-form marginal likelihood
(`priors/niw_logml.m`) plus a Gamma hyperprior with mode 0.4. The BVAR
supplies four things to everything downstream:

1. the **prior centre** for the LP coefficients at every horizon: the
   companion matrix power `A^h` (the LP coefficients a VAR implies are
   exactly its iterated-forecast coefficients; `priors/var_implied_lp_prior.m`
   is the same mapping written out, and `tests/test_var_prior_mapping.m`
   checks it against brute-force iterated forecasts);
2. the **identification**: the shock is the recursive (Cholesky) innovation
   to variable 1, normalised to a unit impact, `b1n = B0(:,1)/B0(1,1)`;
3. the **deterministic trend** that is removed from the data before the
   horizon regressions (`utils/var_deterministic_trend.m`, FMAR's
   convention: the path the fitted VAR would follow without shocks);
4. the VAR's own IRF and credible bands, used as a comparator.

**Step 2: for every horizon `h >= 1`, run the horizon regression.** On the
detrended data `x`, regress `x(t+h)` on `z(t)`. Three prior ingredients are
built per horizon:

* the **centre** `mu_h = [0; (J A^h)']`, the constant centred at 0 because
  the data are detrended;
* the **scales** `psi(h)`: for each variable, the Newey-West long-run
  variance of the residuals of a univariate own-lag LP at horizon `h`
  (`priors/fmar_prior_scale.m`). The prior variance of a coefficient on
  any lag of variable `v` in equation `i` is `Sigma_ii * lambda_h^2 / psi_v(h)`,
  with no Minnesota lag decay (deviations from a VAR-implied centre have no
  reason to shrink with the lag);
* the **tightness** `lambda_h`, maximising the same NIW marginal likelihood
  plus a Gamma hyperprior whose standard deviation loosens with `h`
  (`priors/select_lambda_fmar.m`).

**Step 3a: the global FMAR estimator** (`estimators/estimate_blp_fmar.m`).
With the prior above the posterior mean has the closed form
`(Z'Z + W^-1)^-1 (Z'y + W^-1 mu)`, and the IRF is `theta_i(h) = b1n' * beta_ih(y(t) block)`.
No sampling is needed. Intervals at `h >= 2` are FMAR's "quasi-Bayesian"
Newey-West sandwich around the posterior mean (the horizon-`h` LP error is
MA(h-1), which the per-horizon Gaussian likelihood ignores; the sandwich
is the ex-post correction). `h = 0` is `b1n` for every estimator. At
`h = 1` the published convention reports the BVAR; the project's harmonised
setting `cfg.fmar.h1_mode = 'lp'` estimates `h = 1` as an LP too, so that
the adaptive estimators are compared like for like.

**Step 3b: the block-adaptive estimator** (`estimators/estimate_blp_blockadaptive.m`,
sampler `samplers/gibbs_block_horseshoe.m`). Same data, centre, scales and
`lambda_h` (it inherits `lambda_h` from step 3a so that `tau` is the only
difference). The coefficients of equation `i` at horizon `h` are split into
`G = K` blocks, block `g` = the `p` coefficients on the lags of variable
`g`. The prior on block `g` becomes
`beta_{i,g,h} - mu_{i,g,h} ~ N(0, sigma_i^2 lambda_h^2 tau_{i,g,h}^2 D_g)`,
`tau ~ half-Cauchy(0,1)`. A Gibbs sampler cycles through the coefficients
(Gaussian), the error variance (inverse-gamma), and each `tau_g^2` and its
auxiliary `nu_g` (both inverse-gamma, the Makalic-Schmidt representation).
Because the prior scales with `sigma_i^2`, fixing `tau = 1` makes the
conditional posterior mean of the coefficients exactly the closed form of
step 3a, for any `sigma^2` (`tests/test_fmar_nesting_exact.m` verifies this
to 3e-16 over 912 cells). The reported IRF is the **Rao-Blackwellised**
posterior mean (average of the conditional means over the draws), which
removes simulation noise a closed-form competitor never pays. Primary
bands are the same sandwich as step 3a around this mean; posterior
quantile bands are kept as a secondary output.

**Step 3c: the horizon-pooled estimator** (`estimators/estimate_blp_blockpooled.m`,
sampler `samplers/gibbs_block_pooled_horizons.m`). Identical to 3b except
that the `H` scales of one block are tied together by a random walk on
`log tau_{i,g,h}` across `h`, with innovation variance `kappa^2` sampled per
block. `kappa -> 0` gives one scale per block, `kappa -> infinity` gives
the independent estimator. All horizons of one equation are sampled in
one chain; only the `log tau` path needs Metropolis steps (a checkerboard
sweep over odd and even horizons plus a joint level move). Everything else
is conjugate. `tau = 1` still nests FMAR exactly.

**Step 4: the frequentist comparators.** Plain LP with Newey-West bands
(`estimators/estimate_lp.m`) and, for the empirical chapter, the lag-augmented
LP of Montiel Olea and Plagborg-Moller with Eicker-Huber-White bands
(`empirical/estimators/estimate_lp_lagaug.m`), which needs no HAC choice
and is valid under near-unit roots. Both hold `b1n` fixed like the
Bayesian estimators, so bands differ only through the dynamic coefficients.

**What `tau` means.** A large posterior `tau_{i,g,h}` says: at horizon `h`,
in the equation for `y_i`, the coefficients on the lags of `y_g` want to be
far from where the fitted VAR puts them, relative to the horizon's global
tightness. That is prior-data *disagreement*. Structural misspecification
of the VAR is one cause; sampling noise in the estimated prior centre, a
poorly chosen `lambda_h`, or a weakly identified horizon are others. This
is why every `tau` number in the project is read against a *null*: the
same statistic computed where nothing is misspecified.

---

## 3. The simulation study

**Designs** (`dgp/`). All at `K = 3`, `T = 200`, fitted `p = 2`, a unit
recursive shock to `y_1`, true IRFs computed analytically from the true
parameters and never from an estimator (`dgp/compute_true_irf.m`, checked
against shocked-minus-unshocked simulation in `tests/test_true_irf.m`):

| DGP | truth | what is wrong with the fitted VAR(2) prior |
|---|---|---|
| correct | stable VAR(2) | nothing |
| sparse | VAR(3) with one nonzero entry `A3(3,1) = -0.30` | the "lags of `y_1`" block of the `y_3` equation only |
| intermediate | VAR(3) with a nonzero first column of `A3` | the "lags of `y_1`" block, in every equation |
| dense | VARMA(2,1) with a dense MA matrix | every block of every equation |

**Estimators** compared in every replication on the same simulated data:
LP, BVAR, BLP-FMAR (global), BLP-block (independent `tau`), BLP-pooled.

**Metrics** (`montecarlo/summarize_montecarlo.m`). Bias, Monte Carlo
variance, RMSE, coverage and interval length per (estimator, variable,
horizon); integrated RMSE over `h = 2..H` (the fair convention; the legacy
`h = 1..H` version is kept but not quoted); early (`h = 2..6`) and late
windows. Because all estimators see the same data, differences are tested
with a **paired** Monte Carlo t-statistic (`montecarlo/paired_comparison.m`),
which is roughly ten times more precise than comparing two RMSE levels.

**Tau diagnostics** (`montecarlo/tau_diagnostics.m`). The statistic that
works is the *within-equation, within-horizon contrast*: in each equation,
which block has the largest early-horizon `tau`, and how often across
replications ("concentration"). Its null is the same statistic on the
correct DGP. The rule "some `(i,g)` has posterior `P(tau > 1) > 0.5`"
answers "is anything escaping?"; the argmax map answers "where?". Rules
that average `tau` over equations before comparing blocks do not
discriminate and are retired.

**Headline findings** (`R = 500` for correct and sparse, `250` for the
others; `results/simulation/REPORT.txt`, `docs/RESULTS.md`):

1. The adaptive estimators reliably *lower* RMSE at early horizons and
   reliably *raise* it at late ones; the integrated verdict depends on how
   many late horizons are averaged. Pooling keeps the early gain and halves
   the late penalty, so BLP-pooled is RMSE-neutral against FMAR on the
   sparse DGP and slightly better on the correct one (it can tighten harder
   than the global prior where the VAR fits).
2. Localisation works: on the sparse DGP the true (equation, block) cell
   is ranked first in 97% of replications at early horizons, against 37%
   on the correct DGP. The dense DGP gives the null value (nothing to
   localise, correctly).
3. Two caveats travel with the diagnostic. The null is not a constant (it
   depends on `K`, `T`, `p`, persistence), so it has to be simulated for
   the design at hand. And the map does not say whether escaping will
   *help*: on the intermediate DGP it localises perfectly while adaptation
   is worse at every horizon, because the group horseshoe needs the
   deviation from the centre to be sparse within an equation, not just
   confined to one block.
4. Sensitivity (`montecarlo/run_sensitivity_approximations.m`,
   `docs/APPROXIMATIONS.md`): ignoring cross-equation covariance in the
   sampler is immaterial; inheriting `lambda_h` is not an approximation at
   all (re-selection returns the identical value); holding the prior centre
   at the BVAR posterior mean is the one simplification that matters for
   coverage (16% of the band width), and it hits all estimators alike.

**How the runs were produced.** `montecarlo/run_montecarlo.m` runs one DGP;
replication `r` is seeded as `cfg.seed + 100000*dgp + r`, so the work can be
split into chunks (`run_mc_chunk.m`) across processes and merged exactly
(`merge_montecarlo.m`, verified bit-identical in `tests/test_chunk_merge.m`).
`montecarlo/mc_preset.m` names the configurations (`final`, `grid`,
`legacy`, `smoke`). The 13-cell exploratory grid (`experiment_grid_spec.m`,
`R = 40`) sweeps misspecification strength, sample size, lag order and
sparsity one factor at a time. `report_montecarlo.m` prints every table
from a stored result; `write_results_report.m` regenerates `results/simulation/REPORT.txt`.

---

## 4. The euro-area application (Chapter 7)

**Purpose.** A methods demonstration: run the estimators on real monthly
euro-area data and let the `tau` map answer "where does the VAR prior fail
here?", with the null distribution simulated for this exact design. Design:
`empirical/docs/DESIGN.md`; data: `empirical/docs/DATA.md`; results:
`empirical/docs/RESULTS.md`; critical assessment: `docs/ASSESSMENT.md`.

**Data** (`empirical/data/`). The headline system `ois4` has four monthly
series, 2000m1 to 2019m12 (`T = 240`; with `p = 12` the effective sample
starts in 2001m1):

| # | variable | source | transform |
|---|---|---|---|
| 1 | `ois1m` | 1-month EONIA swap (OIS) rate, end-of-month level (ECB Data Portal copy of the Refinitiv series, bid/ask mid) | percent |
| 2 | `ip` | industrial production, EA20, B-D, seasonally and calendar adjusted (Eurostat `sts_inpr_m`) | 100 log |
| 3 | `hicp_sa` | HICP all items, working-day and seasonally adjusted (ECB `ICP.M.U2.Y.000000.3.INX`) | 100 log |
| 4 | `stoxx` | EURO STOXX 50, monthly average (ECB) | 100 log |

The monetary-policy **instrument** does not enter the VAR. It is the
high-frequency surprise in the 1-month OIS around ECB Governing Council
meetings and Executive Board speeches (EA-EMPD; day-count adjusted as in
Altavilla et al. 2025), summed within the month (`z_gcs_1m_adj_sum`,
`build_instrument_series.m`). It is used once, at identification: the impact
vector is the covariance of the VAR innovations with the instrument,
normalised to a unit effect on the policy rate (`ea_identify_proxy.m`,
Stock–Watson 2018), and every estimator's IRF is its horizon-`h` coefficient
block times that vector. Fifty other instrument variants (meetings only,
3-month and 1-year surprises, Jarociński–Karadi sign splits, Kilian
day-weighted aggregation, orthogonalised on pre-event information, Jarociński's
published shocks) are carried in the dataset for the relevance table.

**Three legs** (`empirical/docs/DESIGN.md`):

* **Leg 1, the finding**: the `tau` heatmap at `p = 12`, `H = 48`, for both
  adaptive estimators (`RUN_EMPIRICAL_IV` steps A and B).
* **Leg 2, the null**: simulate `R` datasets from the BVAR fitted to the
  real data and push each through the whole pipeline
  (`empirical/montecarlo/run_null_calibration.m`). Two residual schemes:
  iid resampling (homoskedastic) and a wild bootstrap that keeps the
  data's volatility profile (the short-rate innovation is four times less
  volatile after 2012). The reading protocol (`ea_apply_protocol.m`)
  reports, per cell, the ratio to the null 95th percentile, a Monte Carlo
  p-value `(r + 1)/(R + 1)` and a family-wise p-value from the max
  statistic over the free cells; a null is refused for a different design
  (`ea_design_key.m`).
* **Leg 3, validation without a truth**: the same estimation at
  `p = 2, 4, 6, 12` against each design's own null (step D); an
  out-of-sample block ablation that forces an escaping cell back to
  `tau = 1` and compares forecasts (Clark–West, Diebold–Mariano; step E);
  the LP-IV benchmark with Anderson–Rubin sets (`estimate_lp_iv.m`).

**Identification strength.** `ea_relevance_iv.m` reports the first stage on
the VAR innovations (effective F against the Montiel Olea–Pflueger
thresholds), lead/lag placebos, predictability from pre-event information and
from the VAR's own lags, influence statistics, and where in the sample the
first-stage covariance comes from. On the headline instrument the effective F
is 13.7 (LP-IV 19.4), but two thirds of the covariance sits in 2009–2011 and
ten months carry half of it; every IRF statement in the chapter is
conditional on that. The `tau` diagnostic is a statement about the dynamic
specification of the four-variable system and does not involve the
instrument.

**Legacy v1.** Until 2026-09-13 the surprise sat inside the VAR as the first
variable. That design (weak first stage on a monthly-average 1-year rate, the
prior-scale floor `cfg.fmar.psi_floor`, the surprise block held at `tau = 1`
by `cfg.blocks.fixed_tau`) is kept runnable in `empirical/legacy_v1/`, its
outputs in `results/empirical/legacy_v1/`, its design in
`empirical/docs/DESIGN_V1_LEGACY.md`.

---

## 5. Every file, in one or two lines

**Top level**

* `RUN_ME_FIRST.m`: the original prototype demonstration (plain Gaussian
  prior, grid tightness, no FMAR machinery), kept as the verified reference.
* `RUN_FMAR_DEMO.m`: one sparse dataset through the FMAR stack, the
  `tau = 1` nesting check, and an optional small Monte Carlo with the full
  report.
* `config/default_config.m`: every setting, with comments. `cfg.mode`
  chooses `'prototype'` or `'fmar'`; `cfg.fmar.*` are the FMAR constants;
  `cfg.blp.pool.*` the horizon-pooling prior; `cfg.dgp.*` the design
  constants of the simulation DGPs.

**Estimation**

* `estimators/estimate_bvar_niw.m`: GLP Bayesian VAR (prior centre,
  identification, trend, comparator).
* `estimators/estimate_blp_fmar.m`: the global FMAR BLP, closed form, with
  the sandwich bands.
* `estimators/estimate_blp_blockadaptive.m`: the block-adaptive BLP,
  `tau` independent by horizon; works in both modes.
* `estimators/estimate_blp_blockpooled.m`: the horizon-pooled BLP.
* `estimators/estimate_lp.m`, `estimate_var.m`, `estimate_blp_global.m`:
  the prototype-mode LP, OLS VAR and global BLP.
* `samplers/gibbs_block_horseshoe.m`: one (equation, horizon) chain; every
  conditional is derived in its header. Used by the global prototype BLP
  (`tau` fixed) and the block-adaptive BLP.
* `samplers/gibbs_block_pooled_horizons.m`: all horizons of one equation
  jointly, with the random-walk prior on `log tau`.
* `priors/niw_logml.m`: closed-form NIW marginal likelihood (port of the
  FMAR/GLP code). `priors/select_lambda_fmar.m`: `lambda_h` by golden
  section on that objective plus the hyperprior. `priors/fmar_prior_scale.m`:
  the Newey-West long-run scales `psi(h)`. `priors/var_implied_lp_prior.m`:
  VAR to LP-coefficient mapping (prototype mode). `priors/build_block_prior.m`:
  block membership and the prototype scaling. `priors/select_global_lambda.m`:
  the prototype grid tightness.
* `utils/`: regressor builder, gamma / inverse-gamma / inverse-Wishart
  draws from `rand`/`randn` only, quantiles, safe Cholesky, the Gamma
  (mode, sd) to (shape, scale) conversion, the VAR trend, and the
  save/load helpers that turn the one anonymous function in `cfg` into text.

**Simulation**

* `dgp/simulate_var_dgp.m`, `simulate_sparse_misspec_dgp.m`,
  `simulate_intermediate_misspec_dgp.m`, `simulate_dense_misspec_dgp.m`:
  the four designs; `simulate_linear_dgp.m` is their shared engine;
  `compute_true_irf.m` the analytic truth.
* `montecarlo/run_montecarlo.m`: one DGP, all estimators, stores everything.
  `summarize_montecarlo.m`: the metrics. `tau_diagnostics.m`: localisation,
  false positives, chain health. `paired_comparison.m`: the paired test.
  `report_montecarlo.m`: prints the tables. `export_montecarlo_csv.m`: tidy
  CSVs plus the meta record. `mc_preset.m`: named configurations.
  `run_mc_chunk.m`, `merge_montecarlo.m`, `compact_mc.m`: exact parallel
  runs. `experiment_grid_spec.m`, `run_experiment_grid.m`,
  `collect_grid_results.m`: the exploratory grid.
  `run_sensitivity_approximations.m`: prices the remaining approximations.
  `write_results_report.m`: regenerates `results/simulation/REPORT.txt`.
* `scripts/run_grid_parallel.sh`, `run_final_parallel.sh`,
  `run_remaining_finals.sh`: shell drivers that launch several Octave
  processes and merge (written for the Linux machine that produced the
  stored runs; paths inside `run_remaining_finals.sh` are that machine's).
* `plots/plot_irfs.m`, `plot_block_scales.m`, `plot_rmse.m`, `can_plot.m`.

**Tests** (`tests/run_all_tests.m` runs all; about 12 minutes in MATLAB)

* mapping and truth: `test_var_prior_mapping`, `test_true_irf`,
  `test_estimator_consistency` (large-`T` convergence of VAR and LP to the
  analytic truth);
* samplers and nesting: `test_sampler_sanity`, `test_nesting` (prototype),
  `test_fmar_nesting` and `test_fmar_nesting_exact` (FMAR),
  `test_pooled_estimator`, `test_approximations`;
* fidelity: `test_fmar_port` (against the original FMAR functions; skipped
  unless `FMAR_PATH` points at their replication package),
  `test_isrw_vector`;
* accounting: `test_mc_metrics`, `test_tau_diagnostics`, `test_chunk_merge`;
* the two empirical-application options: `test_psi_floor` (the prior-scale
  floor is essentially inert on the simulations and keeps the nesting
  exact where it binds), `test_fixed_tau_blocks` (held blocks report
  `tau = 1`, all-held equals FMAR, the default is unchanged);
* the empirical package on a synthetic fixture: `empirical/tests/test_empirical_pipeline.m`.

**Empirical package** (`empirical/`; every entry point resolves its own
paths through `ea_paths.m`, so it runs from any directory)

* `RUN_EMPIRICAL_IV.m`: the driver of the current design (steps REL, A–E).
* `ea_paths.m`: the absolute locations. `ea_design_key.m`: the canonical
  design identifier a null is keyed to. `ea_apply_protocol.m`: the reading
  protocol (ratios, Monte Carlo and family-wise p-values, Holm, counts).
  `ea_cross_p_table.m`: the protocol across lag orders and nulls.
* `data/build_instrument_series.m`: EA-EMPD workbook extract to 51 monthly
  instruments. `data/import_ois_daily.m`: daily 1M OIS to monthly
  end-of-month / average levels. `data/ea_fetch_v2_series.m` and
  `data/fetch_outcome_data.m`: download and validate the ECB / Eurostat
  series. `data/assemble_dataset_v2.m`: the `T x K` dataset with the
  instruments alongside (`ds.Z`). `data/ea_identify_proxy.m`: the proxy
  impact vector. `data/ea_relevance_iv.m`: relevance, timing, influence and
  anatomy diagnostics. `data/ea_pre_event_info.m`,
  `data/import_external_instrument.m`, `data/read_sdmx_csv.m`,
  `data/ea_check_series.m` (coverage and synthetic-data tripwires),
  `data/ea_extract_series.m`, `data/ea_write_provenance.m`.
* `estimators/estimate_lp_iv.m`: LP-IV with HAC, lag-augmented EHW and
  Anderson–Rubin sets. `estimators/estimate_lp_lagaug.m`: the
  lag-augmented LP comparator.
* `dgp/simulate_fitted_bvar_dgp.m`: the bootstrap under the null (iid,
  wild, block, Gaussian residual schemes).
* `montecarlo/run_null_calibration.m`: Leg 2, checkpointed and resumable.
  `montecarlo/run_block_ablation.m`: the out-of-sample block ablation.
* `legacy_v1/`: the 2026-09-12 internal-instrument design
  (`RUN_EMPIRICAL_V1.m`, `SMOKE_TEST_EMPIRICAL_V1.m`, `assemble_dataset.m`,
  `build_shock_series.m`, `ea_relevance.m`), kept runnable.
* `tests/`: synthetic fixtures (`make_synthetic_fixture.m`,
  `make_synthetic_ois_daily.m`) and the package tests
  (`test_empirical_pipeline`, `test_empirical_iv`, `test_null_modularity`,
  `test_lp_iv`); every driver refuses a fixture.
* `docs/DESIGN.md` (current design), `docs/DESIGN_V1_LEGACY.md`,
  `docs/LITERATURE.md`, `docs/DATA.md`, `docs/RESULTS.md` (current results),
  `docs/RESULTS_LOG.md` (the chronological record). `README.md`: how to run.

**Documentation and results**

* `README.md`: purpose, status, how to run. `docs/RESULTS.md`: the headline
  numbers and what they support. `docs/APPROXIMATIONS.md`: the remaining
  simplifications and their measured cost. `docs/REPRODUCE.md`: every
  command with runtimes. This file.
* `results/simulation/`: `mc_headline_<dgp>.mat` and their CSV exports are
  the five-estimator headline runs; `grid/` the exploratory grid;
  `REPORT.txt` the regenerated full report; `mc_final_*` and `mc_fmar_*`
  are earlier four-estimator runs kept for reproducibility.
  `results/empirical/`: the current design's outputs (`iv_*`, nulls,
  ablations, `figures/`); `results/empirical/legacy_v1/`: the v1 outputs.
  See `results/README.md`.
* `presentation/thesis_pitch_2026-09-12.tex`: the eleven-slide status
  deck of 12 September 2026 (simulation evidence, the diagnostic framing,
  the v1 euro-area first pass; it predates the external-instrument
  redesign). `presentation/archive/thesis_pitch_2026-09-08.tex` is the
  earlier version, kept for the record.
