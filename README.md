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

The code now has **two operating modes**, selected by `cfg.mode`:

- **`'prototype'`** (default; original pedagogical machinery, unchanged):
  plain Gaussian likelihood, simple std-ratio prior scales, grid marginal
  likelihood for the tightness, posterior-quantile bands. Kept as the
  verified reference and so every earlier result stays reproducible.
- **`'fmar'`**: the **published methodology of Ferreira,
  Miranda-Agrippino & Ricco (REStat 2025)** as the baseline, ported from
  their replication code (github.com/leonardonferreira/BLP,
  `priorType = 'VAR'`, Cholesky identification) and verified against it
  function-by-function (`tests/test_fmar_port.m`). Ingredients:
  1. GLP Bayesian VAR (`estimators/estimate_bvar_niw.m`): Minnesota NIW
     prior, tightness by marginal likelihood + Gamma(mode .4, sd .2)
     hyperprior; supplies `h ≤ 1` responses, identification, prior
     centres, and the deterministic trend;
  2. horizon regressions on **detrended** data
     (`utils/var_deterministic_trend.m`), LP constant centred at 0;
  3. NIW prior per horizon with **Newey–West long-run scales** from
     univariate own-lag LP residuals (`priors/fmar_prior_scale.m`), no
     Minnesota lag decay under the VAR-centred prior;
  4. per-horizon tightness `lambda_h` by closed-form marginal likelihood
     (`priors/niw_logml.m`) plus FMAR's horizon-dependent Gamma
     hyperprior (`priors/select_lambda_fmar.m`);
  5. **quasi-Bayesian Newey–West sandwich bands** at `h ≥ 2` around the
     posterior-mean IRF (`estimators/estimate_blp_fmar.m`) — FMAR's
     ex-post correction for the MA(h−1) LP errors.

  The **block-adaptive estimator is exactly nested** in this baseline:
  with the conjugate prior scaling (`prior_scales_with_sigma2` in
  `samplers/gibbs_block_horseshoe.m`) and `tau ≡ 1`, its posterior mean
  equals the FMAR closed form for any `sigma²`
  (`tests/test_fmar_nesting.m`).

Documented approximations of the FMAR port (see the file headers):
detrending trims `p` observations instead of FMAR's zero-padding; the
block estimator runs an LP at `h = 1` while the FMAR baseline reports
the BVAR there; the per-equation sampler ignores cross-equation `Sigma`
correlation (irrelevant for the `tau = 1` posterior mean, which is
`sigma²`-free, but not for posterior spread — hence the primary bands
use FMAR's sandwich construction in both estimators). Identification
remains recursive with the impact vector fixed in the `h ≥ 2` bands,
exactly as in the FMAR code.

## 4. Requirements

- MATLAB R2018b or newer (uses only base MATLAB: `rand/randn/rng`, `chol`,
  backslash, basic plotting). **No toolboxes required** — gamma/inverse-gamma
  draws, quantiles and normal quantiles are implemented in `utils/`.
- The Parallel Computing Toolbox is *not* required (plain `for` loops).
- The code also runs under GNU Octave ≥ 8 (used for automated verification
  of this prototype), but MATLAB is the target platform.

## 5. How to run

```matlab
cd block_adaptive_blp     % open MATLAB in the project root
RUN_ME_FIRST              % prototype mode demonstration (unchanged)
RUN_FMAR_DEMO             % FMAR-mode demonstration + nesting check + MC
```

Fast smoke test (small Monte Carlo, short chains):

```matlab
QUICK_DEMO = true; RUN_ME_FIRST      % or: QUICK_DEMO = true; RUN_FMAR_DEMO
```

Run the test suite (set `FMAR_PATH` to the `demo_IRFs` folder of the
FMAR replication package to also run the fidelity test against their
original functions; it is skipped with a message otherwise):

```matlab
cd tests; run_all_tests
```

The Chapter 7 euro-area application is a separate, self-contained
pipeline under `empirical/` — every entry point resolves its own paths,
so it runs from any working directory:

```matlab
build_shock_series();      % EA-EMPD events -> monthly surprises
fetch_outcome_data();      % the four outcome series (or manual route)
assemble_dataset();        % -> empirical/data/ea_dataset.mat
SMOKE_TEST_EMPIRICAL       % ~1 min interface check
RUN_EMPIRICAL              % the multi-hour job
```

Read `empirical/README_EMPIRICAL.md` first; it documents the run order,
the input-validation tripwires, and the one estimator edit the package
required.

## 6. Folder structure

```
RUN_ME_FIRST.m      prototype-mode demonstration (read this to see the flow)
RUN_FMAR_DEMO.m     FMAR-mode demonstration + nesting check + Monte Carlo
config/             default_config.m — every setting incl. cfg.mode, cfg.fmar
dgp/                three simulation DGPs + true-IRF computation
estimators/         estimate_lp / estimate_var / estimate_blp_global /
                    estimate_bvar_niw (GLP Bayesian VAR) /
                    estimate_blp_fmar (published FMAR global BLP) /
                    estimate_blp_blockadaptive (dual-mode)
priors/             var_implied_lp_prior (VAR -> LP prior centre),
                    build_block_prior (blocks and D),
                    select_global_lambda (prototype grid tightness),
                    niw_logml (closed-form NIW marginal likelihood),
                    fmar_prior_scale (NW long-run LP scales),
                    select_lambda_fmar (FMAR tightness selection)
samplers/           gibbs_block_horseshoe.m — all conditionals documented,
                    incl. the conjugate (sigma²-scaled) FMAR variant
montecarlo/         run_montecarlo / summarize_montecarlo (mode-aware)
plots/              plot_irfs / plot_rmse / plot_block_scales
tests/              assertion-based tests + run_all_tests
                    (test_fmar_port validates against the original FMAR
                    code; test_fmar_nesting checks the tau = 1 nesting;
                    test_isrw_vector guards the vectorised isrw prior mean)
utils/              small transparent helpers (regressor builder, gamma /
                    inverse-gamma / quantile draws, safe Cholesky,
                    gamma_coef, draw_iw, var_deterministic_trend)
empirical/          Chapter 7 euro-area application, self-contained:
                    data construction (EA-EMPD surprises + outcome series),
                    RUN_EMPIRICAL driver, SMOKE_TEST_EMPIRICAL,
                    null calibration, docs/CH7_DESIGN.md.
                    See empirical/README_EMPIRICAL.md
results/            figures and .mat output land here (simulation AND
                    empirical)
```

Deviations from the suggested layout: a `utils/` folder was added (shared
helpers used by everything), `dgp/simulate_linear_dgp.m` holds the common
simulation engine of the three DGPs, `priors/select_global_lambda.m`
isolates the tightness selection, and a fourth test (`test_true_irf.m`) plus
`run_all_tests.m` were added.

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

## 8. Simulation designs

1. **Correct** (`simulate_var_dgp`): stable VAR(2); the VAR prior is right.
2. **Sparse** (`simulate_sparse_misspec_dgp`): same VAR(2) plus one omitted
   delayed (lag-3) effect of the shock variable on variable 3 (true model is
   a VAR(3) with a single nonzero entry in `A3`, placed on the shock-
   transmission channel so the prior bias is large enough to detect). The
   fitted VAR(2) prior is wrong mainly for the "lags of variable 1" block in
   the `y_3` equation — the block-adaptive prior's target case.
3. **Dense** (`simulate_dense_misspec_dgp`): VARMA(2,1) with a dense MA
   matrix; the VAR(inf) representation has omitted dynamics in *every*
   block, so local escape has no sparse target.

True IRFs are computed analytically from the true parameters (and
cross-checked by shocked-vs-unshocked simulation with common random
numbers); they never touch any estimator.

## 9. Expected output

`RUN_ME_FIRST` produces: IRF comparison figures (truth, LP, VAR, global BLP,
block-adaptive BLP with credible band) for the correct and sparse DGPs; a
table and figure of posterior block scales `tau_{g,h}` (on the sparse DGP,
block 1 should visibly escape, above all in the equation for `y_3`); and a
Monte Carlo summary with per-horizon bias, RMSE, coverage, average interval
length, integrated RMSE, average block scales, and the probability that the
truly misspecified block receives the largest scale.

## 10. Next research steps

Items 1–4 of the original plan (reproduce FMAR, swap in their prior
scaling and tightness procedure, nest the block layer in that baseline)
are **done** — `cfg.mode = 'fmar'` with fidelity verified against their
replication code (`tests/test_fmar_port.m`) and exact `tau = 1` nesting
(`tests/test_fmar_nesting.m`). Remaining:

1. Re-run the full Monte Carlo (≥ 500 reps) in FMAR mode and compare
   with the prototype-mode results: interval coverage should improve
   markedly (sandwich bands), and the BVAR/BLP-FMAR baselines are much
   stronger comparators than the prototype BLP.
2. Uncertainty refinements: impact-vector (`b1n`) uncertainty is still
   ignored at `h ≥ 2` (as in FMAR); a fully Bayesian alternative would
   propagate BVAR `Bzero` draws. Also consider bands that reflect the
   adaptive prior (the posterior-quantile secondary bands are already
   returned as `.lo_post/.hi_post`).
3. Extend the experiment grid: more DGPs, sample sizes, block schemes
   (per-lag, per-variable-and-lag), and a horizon-pooled version of
   `tau_g` (random-walk prior on `log tau_{g,h}` across `h`).
4. Empirical application: euro-area monetary policy with the EA-EMPD
   high-frequency surprises (monthly aggregation, information-effect
   robustness). **Implemented in `empirical/`** — see
   `empirical/docs/CH7_DESIGN.md`. The K = 5 system mixes a
   white-noise-ish surprise with persistent levels, so it uses the
   per-variable form `cfg.fmar.isrw = [0 1 1 1 1]` rather than a scalar
   `true`; `estimate_bvar_niw` accepts either (scalars behave exactly as
   before, pinned by `tests/test_isrw_vector.m`). What remains is
   obtaining the four outcome series (`fetch_outcome_data` prints the
   manual route) and running the null calibration at `R >= 200`.
