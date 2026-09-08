# Block-Adaptive Bayesian Local Projections — research prototype

Companion code for the MSc thesis project *"Where Does the VAR Prior Fail?
Block-Adaptive Bayesian Local Projections under Sparse Dynamic
Misspecification."*

## 1. Purpose

Local projections (LPs) estimate impulse responses flexibly but noisily;
VARs estimate them efficiently but can be misspecified. Bayesian Local
Projections (BLPs) split the difference by shrinking LP coefficients towards
the values implied by an estimated VAR. This project asks: **when the VAR is
wrong only in a small part of the system, can we let exactly that part of the
LP escape the VAR prior — and does that improve impulse-response
estimation?** The code simulates data where we know precisely where the VAR
prior fails, and compares estimators in a Monte Carlo experiment.

## 2. Core idea

- **Ordinary LP.** For each horizon `h`, regress `y(t+h)` on
  `z(t) = [1; y(t); y(t-1); ...; y(t-p+1)]` by OLS. Unbiased-ish, noisy.
- **Global VAR-centred BLP.** Same regression, plus a Gaussian prior centred
  on the LP coefficients *implied by an estimated VAR(p)*, with one global
  tightness `lambda_h`. Every coefficient is pulled towards the VAR equally.
- **Block-adaptive BLP (the new object).** Partition the coefficients into
  blocks `g = 1..G` (here: one block per variable, collecting that variable's
  lag coefficients) and give each block a local *escape* scale `tau_{g,h}`
  with a half-Cauchy (grouped-horseshoe) prior. Blocks where the VAR prior
  fits keep `tau` small and stay disciplined by the VAR; blocks where the VAR
  is misspecified can inflate `tau` and escape. **The shrinkage is on
  deviations from the VAR centre, not towards zero.**

## 3. Current status — read this first

**What the code contains.** Two operating modes, selected by `cfg.mode`:

- **`'prototype'`** (original pedagogical machinery, unchanged): plain
  Gaussian likelihood, std-ratio prior scales, grid marginal likelihood
  for the tightness, posterior-quantile bands. Kept as the verified
  reference so every earlier result stays reproducible.
- **`'fmar'`**: the published methodology of Ferreira,
  Miranda-Agrippino & Ricco (REStat 2025), ported from their
  replication code (`priorType = 'VAR'`, Cholesky identification) and
  verified function-by-function (`tests/test_fmar_port.m`). This is the
  mode all reported results use.

**Three estimators are kept separately available** in FMAR mode, and
they differ *only* in how the block escape scales `tau_{i,g,h}` are
treated:

| | estimator | `tau` |
|---|---|---|
| 1 | `estimate_blp_fmar` | fixed at 1 (the global FMAR baseline) |
| 2 | `estimate_blp_blockadaptive` | independent half-Cauchy at every horizon |
| 3 | `estimate_blp_blockpooled` | `log tau` smoothed across horizons by a random-walk (first-difference) prior |

Everything else — detrending, prior centre, Newey–West long-run scales,
`lambda_h`, band construction — is shared, and (2) and (3) inherit
`lambda_h` from (1) in the Monte Carlo so `tau` is the only difference.

### What is established

* **Implementation and port.** The FMAR machinery is reproduced
  function-by-function (`tests/test_fmar_port.m`, skipped without their
  replication package).
* **Exact nesting.** With `tau = 1`, both adaptive estimators reproduce
  the FMAR closed-form posterior mean **to machine precision**, not to a
  tolerance. `tests/test_fmar_nesting_exact.m` compares the samplers'
  *conditional* posterior mean — a deterministic quantity — against the
  closed form over **912 (equation, horizon) comparisons** spanning four
  DGPs, two seeds, every equation, every horizon, both `h = 1`
  conventions and both adaptive estimators: worst deviation **3.3e-16**
  against a `1e-10` tolerance.
* **Reproducible Monte Carlo evidence.** Every result carries a
  standardised record (DGP parameters and seeds, `R/T/p/H`, estimator
  settings, all metrics, `tau` statistics, sampler diagnostics) and is
  exported to tidy CSV. Runs split across processes are bit-identical to
  the serial run (`tests/test_chunk_merge.m`).
* **The benchmark comparison is now fair.** See below.

### What is promising but not settled

* Localisation of *sparse* prior–data conflict: on the sparse DGP the
  posterior `tau` does concentrate on the block the misspecification was
  put in, and the correct-DGP false-positive rates of the same rules are
  reported alongside so the detection numbers can be read.
* Early-horizon RMSE gains on the DGP the method was designed for.

### What is NOT established

* **Integrated RMSE superiority.** The adaptive estimators do not beat
  the global FMAR baseline on integrated RMSE in general.
* **Broad generalisation.** The experiment grid varies one factor at a
  time around a baseline and cannot identify interactions.
* **Empirical usefulness.** The Chapter 7 application has no completed
  dataset: the four outcome series are licensed and not in the
  repository. Nothing empirical has been estimated, and the package
  refuses to run on the synthetic fixture that makes it testable.

See `docs/RESULTS.md` for the numbers and `docs/APPROXIMATIONS.md` for
what the remaining simplifications cost.

## 3a. Comparing the estimators fairly

Two changes were needed before any RMSE comparison meant what it looked
like it meant.

**The `h = 1` mismatch.** FMAR report the Bayesian VAR at `h = 1`, while
the block-adaptive estimators run an `h = 1` local projection. Any
`h = 1` RMSE gap therefore mixed "VAR vs LP" into what was supposed to
be "global vs adaptive" — and the adaptive estimators were additionally
being handed the BVAR's *Minnesota* tightness as if it were an LP
tightness (0.099 vs 0.270 on a test dataset). Two fixes:

* the reported metric is now **integrated RMSE over `h = 2..H`**
  (`s.irmse_h2`), which is invariant to the `h = 1` convention. The
  legacy `h = 1..H` metric (`s.irmse`) is kept unchanged for
  compatibility and is *not* what the write-up should quote;
* `cfg.fmar.h1_mode = 'lp'` makes the global baseline run an `h = 1`
  local projection with the same FMAR machinery, so `h = 1` is
  like-for-like and the `tau = 1` nesting is exact from `h = 1` upwards.
  The FMAR convention remains the default (`'bvar'`).

**Simulation noise in the point estimate.** The FMAR baseline is a
closed form with no sampler; the adaptive estimators are sampled.
Charging them for avoidable Monte Carlo noise is not a fair comparison,
so their posterior mean is now **Rao-Blackwellised** by default: the
average of the conditional means `E[beta | tau, sigma2, y]` rather than
of the draws. Same estimand, roughly half the Monte Carlo error, and
exactly the closed form when `tau` is fixed.
`cfg.blp.point_estimate = 'draw_mean'` restores the old behaviour.

## 4. Requirements

- MATLAB R2018b or newer (uses only base MATLAB: `rand/randn/rng`, `chol`,
  backslash, basic plotting). **No toolboxes required** — gamma/inverse-gamma
  draws, quantiles and normal quantiles are implemented in `utils/`.
- The Parallel Computing Toolbox is *not* required (plain `for` loops).
- The code also runs under GNU Octave ≥ 8 (every result reported here was
  produced and verified under Octave 8.4), but MATLAB is the target
  platform. Two Octave-specific notes: the demonstration scripts skip
  their figures when no graphics toolkit is installed
  (`plots/can_plot.m`), and a saved configuration struct stores its one
  anonymous function as source text because Octave's `save` refuses
  function handles (`utils/cfg_to_savable.m` / `cfg_from_saved.m`).
- The parallel drivers in `scripts/` are plain shell: they launch several
  independent Octave processes and merge the results exactly. No
  toolbox, and the merged result is bit-identical to a serial run
  (`tests/test_chunk_merge.m`).

## 5. How to run

```matlab
cd block_adaptive_blp     % open MATLAB in the project root
RUN_ME_FIRST              % prototype mode demonstration (unchanged)
RUN_FMAR_DEMO             % FMAR-mode demonstration + exact nesting check + MC
QUICK_DEMO = true; RUN_FMAR_DEMO     % ~2 minutes
```

Test suite (set `FMAR_PATH` to the `demo_IRFs` folder of the FMAR
replication package to also run the fidelity test against their
original functions; it is skipped with a message otherwise):

```matlab
cd tests; run_all_tests
```

The Monte Carlo experiments, the exploratory grid, the sensitivity
analysis and the empirical pipeline each have their own command; they
are all listed, with expected runtimes and outputs, in
**`docs/REPRODUCE.md`**. The short version:

```bash
scripts/run_grid_parallel.sh 4                  # exploratory grid, ~1.5 h
scripts/run_final_parallel.sh sparse  500 4     # headline run, ~3 h
```

```matlab
report_montecarlo('results/mc_headline_sparse.mat')   % all headline tables
```

## 6. Folder structure

```
RUN_ME_FIRST.m      prototype-mode demonstration (read this to see the flow)
RUN_FMAR_DEMO.m     FMAR-mode demonstration + exact nesting check + MC
config/             default_config.m — every setting incl. cfg.mode, cfg.fmar,
                    cfg.blp.pool (horizon pooling), cfg.dgp (design constants)
dgp/                four simulation DGPs (correct / sparse / intermediate /
                    dense) + true-IRF computation
estimators/         estimate_lp / estimate_var / estimate_blp_global /
                    estimate_bvar_niw (GLP Bayesian VAR) /
                    estimate_blp_fmar (published FMAR global BLP) /
                    estimate_blp_blockadaptive (tau independent by horizon) /
                    estimate_blp_blockpooled  (log tau pooled across horizons)
priors/             var_implied_lp_prior, build_block_prior,
                    select_global_lambda, niw_logml, fmar_prior_scale,
                    select_lambda_fmar
samplers/           gibbs_block_horseshoe.m        one (equation, horizon)
                    gibbs_block_pooled_horizons.m  all horizons of one
                                                   equation jointly
montecarlo/         run_montecarlo / summarize_montecarlo / report_montecarlo
                    tau_diagnostics / export_montecarlo_csv
                    mc_preset (named configurations: final / grid / legacy)
                    run_mc_chunk + merge_montecarlo (exact parallel runs)
                    experiment_grid_spec / run_experiment_grid /
                    collect_grid_results
                    run_sensitivity_approximations
scripts/            run_grid_parallel.sh, run_final_parallel.sh
plots/              plot_irfs / plot_rmse / plot_block_scales
tests/              assertion-based tests + run_all_tests
utils/              regressor builder, gamma / inverse-gamma / quantile
                    draws, safe Cholesky, gamma_coef, draw_iw,
                    var_deterministic_trend, cfg_to_savable / cfg_from_saved
docs/               APPROXIMATIONS.md  what the remaining simplifications cost
                    REPRODUCE.md       every command, with runtimes
                    RESULTS.md         the headline numbers and how to read them
empirical/          Chapter 7 euro-area application, self-contained;
                    empirical/tests/ holds the synthetic fixture and the
                    pipeline test.  See empirical/README_EMPIRICAL.md
results/            figures, .mat and .csv output (simulation AND empirical)
```

## 7. Main mathematical specification

LP regression at horizon `h` (per equation, Gaussian likelihood):

```
y(t+h) = X_t' beta_h + u(t+h)
```

Prior centred on the VAR-implied coefficients `mu_h^VAR` with block-local
escape scales:

```
beta_{g,h} - mu_{g,h}^VAR  ~  N( 0,  lambda_h^2 * tau_{g,h}^2 * D_{g,h} )
tau_{g,h}                  ~  half-Cauchy(0, 1)
```

`tau_{g,h} = 1` for all `g` recovers the global VAR-centred prior
(`tests/test_nesting.m` verifies this numerically). The centre `mu_h^VAR` is
the exact conditional-mean mapping of the estimated VAR: the coefficient
block on `y(t)` equals the reduced-form IRF `Psi_h`, verified against
iterated forecasts in `tests/test_var_prior_mapping.m`. Posterior sampling
is a Gibbs sampler with the Makalic–Schmidt auxiliary representation of the
half-Cauchy; every conditional is derived in the header of
`samplers/gibbs_block_horseshoe.m`.

**FMAR mode.** With `cfg.mode = 'fmar'` the same block structure is
placed on top of the published FMAR prior. At horizon `h`, on detrended
data, equation `i`:

```
beta_ih | sigma_i², tau  ~  N( mu_ih^BVAR ,  sigma_i² · W_ih ),
W_ih = diag(w),  w_const = Vc = 1e5,
w_j  = lambda_h² · tau_{g(j),h}² / psi_v(h)   (coefficient on any lag of variable v),
sigma_i² ~ IG( 3/2, psi_i(h)/2 )              (exact IW-diagonal marginal),
```

with `psi(h)` the Newey–West long-run variances of univariate own-lag LP
residuals (Bartlett, truncation `h+1`) and `lambda_h` maximising the
closed-form NIW marginal likelihood plus the Gamma hyperprior with mode
0.4 and sd `0.1 + 0.4/(1+exp(-0.3(h-12)))`, on `[1e-4, 5]`. Because the
prior scales with `sigma_i²`, the `tau ≡ 1` conditional posterior mean is
`(Z'Z + W⁻¹)⁻¹(Z'y + W⁻¹ mu)` for any `sigma_i²` — the FMAR closed form,
hence exact nesting. Primary intervals at `h ≥ 2` are FMAR's
quasi-Bayesian Newey–West sandwich around the posterior-mean IRF with the
impact vector fixed; posterior-quantile bands are returned alongside.

**Horizon pooling of `tau`.** In the estimator above, each
`tau_{i,g,h}` is identified by a single horizon-`h` regression, i.e. by
`p_g = p` coefficient deviations. That is very little information: the
`tau` path is noisy, and the noise is passed straight into the posterior
mean of `beta_h`, adding variance at horizons where the VAR prior is in
fact fine. Since "how badly does the VAR prior fail for block `g`?"
should vary smoothly with the horizon, the third estimator puts a
first-difference shrinkage prior on the path, writing
`x_{i,g,h} = log tau_{i,g,h}`:

```
x_{i,g,1}  ~ half-Cauchy(0,1) on the tau scale   (same anchor as before)
x_{i,g,h}  = phi * x_{i,g,h-1} + e,   e ~ N(0, kappa_{i,g}^2),   h >= 2
kappa^2    ~ IG(a_kappa, b_kappa)
```

with `phi = 1` (a random walk) by default. This nests both extremes:
`kappa -> 0` gives one scale per block for all horizons, `kappa -> inf`
returns the independent estimator, and sampling `kappa` lets the data
choose per equation and block. All horizons of one equation are sampled
jointly (`samplers/gibbs_block_pooled_horizons.m`); the coefficient and
variance steps stay conjugate, and only the `log tau` path needs
Metropolis — two checkerboard half-sweeps (sites at odd horizons are
conditionally independent given those at even horizons) plus a joint
level move on the whole path, five times per Gibbs iteration. Fixing
`tau = 1` here reproduces the FMAR closed form exactly, as it does for
the independent estimator.

## 8. Simulation designs

1. **Correct** (`simulate_var_dgp`): stable VAR(2); the VAR prior is right.
2. **Sparse** (`simulate_sparse_misspec_dgp`): same VAR(2) plus one omitted
   delayed (lag-3) effect of the shock variable on variable 3 (true model is
   a VAR(3) with a single nonzero entry in `A3`, placed on the shock-
   transmission channel so the prior bias is large enough to detect). The
   fitted VAR(2) prior is wrong mainly for the "lags of variable 1" block in
   the `y_3` equation — the block-adaptive prior's target case.
3. **Intermediate** (`simulate_intermediate_misspec_dgp`): the lag-3
   matrix has a nonzero *column* — every equation omits a delayed effect
   of the shock variable. The target block `g* = 1` is still well
   defined, but the signal is spread over `K` equations instead of
   concentrated in one. This design separates "block-sparse in the
   coefficient space" from "block-sparse in the equation space".
4. **Dense** (`simulate_dense_misspec_dgp`): VARMA(2,1) with a dense MA
   matrix; the VAR(inf) representation has omitted dynamics in *every*
   block, so local escape has no sparse target.

The strength of each design is a configuration constant
(`cfg.dgp.sparse_a31`, `cfg.dgp.interm_scale`, `cfg.dgp.dense_scale`), so
the experiment grid can vary misspecification strength while holding
everything else fixed. Defaults reproduce the original designs exactly.

True IRFs are computed analytically from the true parameters (and
cross-checked by shocked-vs-unshocked simulation with common random
numbers); they never touch any estimator.

## 9. Expected output

`RUN_FMAR_DEMO` produces: an IRF comparison figure (truth, LP, BVAR,
global BLP-FMAR, block-adaptive, horizon-pooled) on the sparse DGP; the
`tau` surface under both adaptive estimators, with `P(tau > 1)` and the
posterior smoothing scale `kappa`; the exact `tau = 1` nesting check;
and, when the demonstration Monte Carlo is enabled, the full report from
`report_montecarlo` — both integrated-RMSE conventions, the
bias–variance decomposition by horizon window, interval coverage under
both band constructions, and the `tau` localisation and false-positive
tables.

`docs/RESULTS.md` holds the headline numbers from the `R = 500` runs and
says explicitly which conclusions they do and do not support.

## 10. Next research steps

1. **A non-centred parameterisation of the pooled sampler.** `log tau`
   mixes slowly (effective sample sizes of order 10–70 out of 700). It
   does not currently matter, because the reported IRF is
   Rao-Blackwellised and its Monte Carlo error is small relative to the
   posterior spread, but it would matter if the `tau` path itself became
   the object of inference.
2. **Joint selection of `lambda_h` and `tau`.** The adaptive estimators
   currently inherit `lambda_h` from the global BLP, which isolates the
   effect of `tau` but is not what a stand-alone adaptive estimator
   would do (`docs/APPROXIMATIONS.md`, item 2c).
3. **Propagating BVAR uncertainty.** The prior centre and the impact
   vector are held at posterior means, as in FMAR. The sensitivity
   analysis prices both; a fully Bayesian version would integrate over
   the BVAR posterior.
4. **A factorial rather than one-factor-at-a-time grid**, to identify
   interactions — above all between misspecification sparsity and sample
   size.
5. **The euro-area application.** The pipeline runs; what is missing is
   the four licensed outcome series and the `R >= 200` null calibration.
   See `empirical/README_EMPIRICAL.md`.
