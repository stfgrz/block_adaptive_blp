# The method

The estimator, the diagnostic and the conventions under which the three
estimators are compared. Numbers live in `RESULTS.md`; the plain-language
tour of the code is `PROJECT_EXPLAINED.md`; the remaining simplifications
and what they cost are in `APPROXIMATIONS.md`.

## 1. Core idea

Local projections (LPs) estimate impulse responses flexibly but noisily;
VARs estimate them efficiently but can be misspecified. Bayesian Local
Projections (BLPs) shrink LP coefficients towards the values implied by an
estimated VAR. The question of this project: **when the VAR is wrong only in
a small part of the system, can exactly that part of the LP be let out of
the VAR prior, and does that help?**

- **Ordinary LP.** For each horizon `h`, regress `y(t+h)` on
  `z(t) = [1; y(t); y(t-1); ...; y(t-p+1)]` by OLS.
- **Global VAR-centred BLP.** The same regression with a Gaussian prior
  centred on the LP coefficients *implied by an estimated VAR(p)* and one
  global tightness `lambda_h` per horizon: every coefficient is pulled
  towards the VAR equally. The published version is Ferreira,
  Miranda-Agrippino and Ricco (REStat 2025, "FMAR").
- **Block-adaptive BLP.** Partition the coefficients into blocks
  `g = 1..G` (one block per variable, collecting that variable's lag
  coefficients) and give each block a local *escape* scale `tau_{g,h}` with
  a half-Cauchy (grouped-horseshoe) prior. Blocks where the VAR prior fits
  keep `tau` small; blocks where it is misspecified can inflate `tau` and
  escape. **The shrinkage is on deviations from the VAR centre, not towards
  zero.** `tau = 1` recovers the global prior exactly.

## 2. The three estimators

All reported results use FMAR mode (`cfg.mode = 'fmar'`): the published
machinery ported from the FMAR replication code and verified function by
function (`tests/test_fmar_port.m`). The three estimators differ *only* in
how the block scales are treated:

| | estimator | `tau_{i,g,h}` |
|---|---|---|
| 1 | `estimate_blp_fmar` | fixed at 1 (the global FMAR baseline, closed form) |
| 2 | `estimate_blp_blockadaptive` | independent half-Cauchy at every horizon |
| 3 | `estimate_blp_blockpooled` | `log tau` smoothed across horizons by a random-walk prior |

Detrending, prior centre, Newey–West long-run scales, `lambda_h` and band
construction are shared, and (2) and (3) inherit `lambda_h` from (1), so
`tau` is the only difference. With `tau = 1` both adaptive estimators
reproduce the FMAR closed-form posterior mean to machine precision
(`tests/test_fmar_nesting_exact.m`: 912 (equation, horizon) comparisons,
worst deviation 3.3e-16). The prototype mode (`cfg.mode = 'prototype'`,
plain Gaussian likelihood, grid tightness) is kept as the verified
pedagogical reference and is not used for any reported number.

## 3. Specification

LP regression at horizon `h`, per equation `i`, Gaussian likelihood:

```
y(t+h) = X_t' beta_h + u(t+h)
```

Prior centred on the VAR-implied coefficients `mu_h^VAR` with block-local
escape scales:

```
beta_{g,h} - mu_{g,h}^VAR  ~  N( 0,  lambda_h^2 * tau_{g,h}^2 * D_{g,h} )
tau_{g,h}                  ~  half-Cauchy(0, 1)
```

The centre `mu_h^VAR` is the exact conditional-mean mapping of the estimated
VAR: the coefficient block on `y(t)` equals the reduced-form IRF `Psi_h`
(`tests/test_var_prior_mapping.m`). Posterior sampling is a Gibbs sampler
with the Makalic–Schmidt auxiliary representation of the half-Cauchy; every
conditional is derived in the header of `samplers/gibbs_block_horseshoe.m`.

**FMAR mode.** On detrended data, at horizon `h`, equation `i`:

```
beta_ih | sigma_i², tau  ~  N( mu_ih^BVAR ,  sigma_i² · W_ih ),
W_ih = diag(w),  w_const = Vc = 1e5,
w_j  = lambda_h² · tau_{g(j),h}² / psi_v(h)   (coefficient on any lag of variable v),
sigma_i² ~ IG( 3/2, psi_i(h)/2 )              (exact IW-diagonal marginal),
```

with `psi(h)` the Newey–West long-run variances of univariate own-lag LP
residuals (Bartlett, truncation `h+1`) and `lambda_h` maximising the
closed-form NIW marginal likelihood plus the Gamma hyperprior with mode 0.4
and sd `0.1 + 0.4/(1+exp(-0.3(h-12)))`, on `[1e-4, 5]`. Because the prior
scales with `sigma_i²`, the `tau ≡ 1` conditional posterior mean is
`(Z'Z + W⁻¹)⁻¹(Z'y + W⁻¹ mu)` for any `sigma_i²`: the FMAR closed form,
hence exact nesting. Primary intervals at `h ≥ 2` are FMAR's quasi-Bayesian
Newey–West sandwich around the posterior-mean IRF with the impact vector
fixed; posterior-quantile bands are returned alongside (and undercover; see
`RESULTS.md`).

**Horizon pooling of `tau`.** In estimator (2) each `tau_{i,g,h}` is
identified by a single horizon-`h` regression, i.e. by `p` coefficient
deviations. That is very little information: the `tau` path is noisy, and
the noise is passed straight into the posterior mean of `beta_h`. Since "how
badly does the VAR prior fail for block `g`?" should vary smoothly with the
horizon, estimator (3) puts a first-difference shrinkage prior on the path,
`x_{i,g,h} = log tau_{i,g,h}`:

```
x_{i,g,1}  ~ half-Cauchy(0,1) on the tau scale   (same anchor as before)
x_{i,g,h}  = phi * x_{i,g,h-1} + e,   e ~ N(0, kappa_{i,g}^2),   h >= 2
kappa^2    ~ IG(a_kappa, b_kappa)
```

with `phi = 1` (a random walk). `kappa -> 0` gives one scale per block for
all horizons, `kappa -> inf` returns the independent estimator, and sampling
`kappa` lets the data choose per equation and block. All horizons of one
equation are sampled jointly (`samplers/gibbs_block_pooled_horizons.m`); the
coefficient and variance steps stay conjugate, only the `log tau` path needs
Metropolis (two checkerboard half-sweeps plus a joint level move, five times
per Gibbs iteration). One consequence matters for the empirical protocol:
the pooled `tau` at `h ≤ 12` depends on the horizons above 12 through the
smoothing, so it must be compared with a null simulated at the same `H`.

## 4. Comparing the estimators fairly

Two changes were needed before any RMSE comparison meant what it looked like
it meant.

**The `h = 1` mismatch.** FMAR report the Bayesian VAR at `h = 1`, while the
block-adaptive estimators run an `h = 1` local projection, so any `h = 1`
RMSE gap mixed "VAR vs LP" into "global vs adaptive", and the adaptive
estimators were handed the BVAR's Minnesota tightness as if it were an LP
tightness. Two fixes: the reported metric is the **integrated RMSE over
`h = 2..H`** (`s.irmse_h2`), invariant to the `h = 1` convention; and
`cfg.fmar.h1_mode = 'lp'` makes the global baseline run an `h = 1` local
projection with the same machinery, so the `tau = 1` nesting is exact from
`h = 1` upwards. The FMAR convention remains the default (`'bvar'`).

**Simulation noise in the point estimate.** The FMAR baseline is a closed
form with no sampler; the adaptive estimators are sampled. Their posterior
mean is therefore **Rao-Blackwellised** by default: the average of the
conditional means `E[beta | tau, sigma2, y]` rather than of the draws. Same
estimand, roughly half the Monte Carlo error, and exactly the closed form
when `tau` is fixed. `cfg.blp.point_estimate = 'draw_mean'` restores the
old behaviour.

## 5. The diagnostic and its null

A large posterior `tau_{i,g,h}` says that the horizon-`h` LP coefficients
of block `g` in equation `i` disagree with the VAR-implied centre, relative
to the horizon's global tightness. That is **prior–data disagreement**, not
misspecification: sampling noise in the estimated centre, a badly chosen
`lambda_h`, weak identification of a horizon, a change in volatility or in
the relationship over the sample are other causes. The simulation study
found that the false-positive rate of any threshold rule depends on `K`,
`T`, the fitted lag order and the persistence of the system (the
false-positive baseline of the localisation statistic is 0.42 at `p = 2`
and 0.62–0.68 at `p = 3`/`p = 4`, where the fitted VAR nests the truth).
Two consequences shape everything downstream:

1. **No threshold transfers across specifications.** The null has to be
   simulated for the design at hand: datasets are generated from the VAR
   fitted to the data (so the VAR-centred prior is correct by
   construction), the whole pipeline is re-run on each, and the real-data
   statistic is read against the resulting distribution
   (`empirical/montecarlo/run_null_calibration.m`). Because the innovations
   of real data are not homoskedastic, the calibration offers an iid
   residual-resampling null and a wild-bootstrap null that keeps each
   month's own residual variance; an escape should survive both.
2. **The right statistic is a within-equation, within-horizon contrast.**
   `tau` rises with the horizon under every data-generating process. The
   statistics that work compare blocks within an equation (the argmax and
   its concentration across replications; a cell's early-horizon mean scale
   against its own null quantile), not thresholds on the level of `tau`
   averaged over equations, which do not discriminate at all.

In the empirical protocol (`empirical/ea_apply_protocol.m`) each cell
reports its early-horizon mean scale relative to the null 95th percentile,
a Monte Carlo p-value `(r + 1)/(R + 1)`, and a family-wise p-value from the
maximum of the ratio over the free cells in each null draw (single-step
max-T adjustment). An escape is *reported* when the family-wise p-value is
below 5% under the pre-specified null; it is *validated* only if forcing the
cell back to `tau = 1` costs out-of-sample forecast accuracy
(`empirical/montecarlo/run_block_ablation.m`). The simulation study showed
why the second step is not optional: on the intermediate DGP the map
localises the conflict perfectly while adaptation hurts at every horizon.

## 6. Identification in the empirical chapter

The euro-area system carries no instrument in its state vector. The
high-frequency policy surprise `z_t` enters once: the impact vector is
`b_z = Cov(u_t, z_t) / Cov(u_{s,t}, z_t)` with `u_t` the VAR innovations and
`s` the policy rate (Stock–Watson 2018), and every estimator's structural IRF
is its horizon-`h` coefficient block times `b_z` (`empirical/data/ea_identify_proxy.m`).
Under the maintained VAR(p) structure this equals the LP-IV estimand, which
`empirical/estimators/estimate_lp_iv.m` estimates directly with HAC bands and
Anderson–Rubin sets that stay valid when the instrument is weak. The `tau`
map does not depend on the instrument at all; the IRFs depend on it
entirely, and their credibility is the first stage's, which
`empirical/data/ea_relevance_iv.m` reports together with where in the sample
the identifying covariance comes from.
