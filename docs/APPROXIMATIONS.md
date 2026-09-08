# Remaining methodological approximations

This note lists every place where the block-adaptive estimator departs
from a fully Bayesian treatment of the model it claims to fit, says what
each departure could cost, and points at the code that measures it.
Nothing here is a defect that was discovered late: each item is a
deliberate choice, most of them inherited from the published FMAR
baseline so that the adaptive layer is compared against that baseline on
its own terms. The point of the note is that "inherited from FMAR" is
not by itself an argument that a simplification is harmless *for the
adaptive layer*, which asks a different question of the same posterior.

The numbers quoted below come from
`montecarlo/run_sensitivity_approximations.m`; re-run it with

```matlab
S = run_sensitivity_approximations(struct('mc_file', ...
        fullfile('results', 'mc_headline_sparse.mat')));
```

which writes `results/sensitivity_approximations.csv`.

## The numbers, in one table

Sparse DGP, 10 simulated datasets, exploratory preset. Every entry is
the **maximum over all (equation, horizon) cells** of a replication,
averaged over replications — a worst-cell measure, not a typical one.
Each variant is run from the **same seed** as the baseline, and the
first row is the yardstick: the identical estimator run from a
*different* seed. A modelling change that moves the IRF by less than
that has not been shown to matter at all.

| | max abs IRF change | ... / posterior sd | max abs change in log tau |
|---|---|---|---|
| **noise baseline** (same estimator, different seed) | 0.0125 | 0.198 | 0.558 |
| A. sigma2 fixed at the system NIW value | 0.0149 | 0.231 | 0.556 |
| B3. lambda re-selected instead of inherited | **0.0000** | **0.000** | — |

| | sd contributed | reported half band-width | ratio |
|---|---|---|---|
| B1. impact vector `b1n` held fixed | 0.0014 | 0.1344 | **0.010** |
| B2. prior centre held at the BVAR posterior mean | 0.0222 | 0.1344 | **0.165** |

**What this says.**

* **A is immaterial.** Fixing `sigma2` at a value that used the
  cross-equation information moves the IRF by 0.0149 against a
  pure-resampling baseline of 0.0125, and moves `log tau` by 0.556
  against a baseline of 0.558 — i.e. not at all. Combined with the exact
  result at `tau = 1` (the change is *identically zero*, verified in
  `tests/test_approximations.m`), the per-equation treatment of `Sigma`
  is not a live concern for the point estimate. It remains a real
  restriction on the posterior *spread*, which is why the primary bands
  are not posterior quantiles.
* **B3 is not an approximation.** The adaptive estimator's own
  marginal-likelihood selection returns *exactly* the global BLP's
  `lambda_h` (largest gap across datasets: 0.00e+00) — same objective,
  same data, both evaluated at `tau = 1`. Inheriting it changes nothing.
  What remains an approximation is that `lambda_h` is chosen at
  `tau = 1` and never re-selected jointly with `tau`; that is a
  different, larger model and is listed as future work.
* **B1 is negligible, B2 is the one that matters.** Impact-vector
  uncertainty contributes 1% of the reported half band-width;
  prior-centre uncertainty contributes 16.5%, sixteen times more. If any
  of the three fixed quantities is to be made stochastic, it is the
  prior centre. Note that both channels affect the global baseline and
  the adaptive estimators *identically*, so neither can explain an RMSE
  difference between them — they matter for coverage, not for the
  comparison.

---

## 1. The adaptive sampler ignores cross-equation covariance

**What is done.** `samplers/gibbs_block_horseshoe.m` runs one chain per
(equation, horizon). Equation *i* carries a scalar error variance
`sigma2_i` with an Inverse-Gamma prior; the off-diagonal elements of the
system covariance `Sigma` never enter. The FMAR baseline, by contrast,
is a genuine system NIW posterior in which `Sigma` is a full matrix.

**Why it is defensible.** The IG prior used per equation is not an
invention: `IG(3/2, psi_i/2)` is the *exact marginal* of the *i*-th
diagonal element of `Sigma ~ IW(diag(psi), K+2)`, so the per-equation
model matches the system prior in every respect a single equation can
see.

**What it provably cannot affect.** With `tau = 1` the conditional
posterior mean of `beta` is
`(Z'Z + W^{-1})^{-1} (Z'y + W^{-1} mu)`, which does not contain
`sigma2` at all — it cancels between the likelihood and the
`sigma2`-scaled prior. So the point estimate of the *global* baseline is
untouched by this approximation. That is not an argument, it is an
identity, and `tests/test_fmar_nesting_exact.m` verifies the resulting
equality numerically over 912 (equation, horizon) comparisons: worst
deviation **3.3e-16**.

**What it can affect.** Two things.
1. *Posterior spread.* The per-equation posterior of `beta_i` conditions
   on `sigma2_i` alone and cannot borrow information from the other
   equations' residuals, so the posterior-quantile bands are not the
   system posterior's bands. This is one reason the **primary** intervals
   are not posterior quantiles but FMAR's sandwich (item 3).
2. *The point estimate, once `tau` is sampled.* The `tau` conditional
   (3') divides the block deviation by `sigma2_i`. A different treatment
   of `sigma2_i` therefore moves the `tau` path, and through it the
   posterior mean. This channel is real and is the one that had to be
   measured rather than argued away.

**How it is measured.** `cfg.blp.sigma2_mode = 'fixed_niw'` re-runs the
estimator with `sigma2_i` held at the *i*-th diagonal of the system NIW
posterior mode at the same `lambda_h` — a value that *did* use the
cross-equation information. Part A of the sensitivity script reports the
resulting change in the IRF, both in IRF units and relative to the
posterior standard deviation of the IRF, and the change in `log tau`.

---

## 2. Prior centre, `lambda_h` and the impact vector are held fixed

**What is done.** At `h >= 2` the estimators condition on
* the BVAR **posterior mean** for the prior centre `mu_h = (J F^h)'`,
* a single **`lambda_h`**, inherited from the global BLP so the
  estimators differ only through `tau`,
* the impact vector **`b1n`** at its point estimate.

All three are exactly what FMAR do. All three understate uncertainty,
and the first two also fix a quantity the adaptive layer is
*interpreting*: `tau` measures the distance between the LP coefficients
and this particular centre, at this particular tightness.

**2a. Impact vector.** `theta_i(h) = b1n' beta_h(2:1+K)` is **linear**
in `b1n`, so this channel can be isolated exactly, with no
re-estimation: hold the coefficient posterior fixed and re-evaluate
`theta` at each NIW posterior draw of `b1n` (`bvar.b1n_draws`). Part B1
of the sensitivity script reports the induced standard deviation against
the reported half band-width. Note that this channel affects the
*global baseline and the adaptive estimators identically*, so it cannot
explain any RMSE difference between them; it matters for coverage, not
for the comparison.

**2b. Prior centre.** Part B2 re-runs the whole adaptive estimator with
the centre and `b1n` taken from BVAR posterior draws instead of the
posterior mean, and reports the spread of the resulting IRFs. This is
the honest version of "the prior centre is itself estimated": at
`T = 200` the centre is not sharp, and a `tau` that escapes may be
escaping from sampling noise in the centre rather than from a structural
defect. **This is the main reason a large `tau` is reported in this
project as evidence of prior-data disagreement, not of
misspecification.**

**2c. `lambda_h`.** In the Monte Carlo the adaptive estimators inherit
the global BLP's `lambda_h`. That is a deliberate design choice — it
makes `tau` the *only* difference between the three estimators — but it
is not what a stand-alone adaptive estimator would do. Part B3 compares
the inherited `lambda_h` against re-selecting it by marginal likelihood
inside the adaptive estimator, and reports both the `lambda` values and
the IRF change. A fully coherent treatment would select `lambda_h` and
`tau` jointly; that is a larger model and is listed as future work
rather than attempted here.

---

## 3. Quasi-Bayesian Newey–West bands

**What is done.** The primary intervals at `h >= 2` are FMAR's
frequentist sandwich `(Z'Z)^{-1} S (Z'Z)^{-1}` with Bartlett truncation
`L = h + 1`, evaluated at the posterior-mean residuals and centred at
the posterior-mean IRF, with `b1n` fixed.

**Why.** The per-horizon Gaussian likelihood is *wrong*: under correct
specification the `h`-step LP error is MA(`h-1`), so the posterior it
implies is overconfident at long horizons. FMAR's ex-post sandwich is
their correction, and using the same construction for all three
estimators means a coverage comparison reflects the point estimator, not
two different band recipes.

**What it costs.** The bands are no longer a Bayesian credible set: they
are a frequentist interval drawn around a Bayesian point estimate, and
they inherit none of the adaptive prior's information about `tau`. In
particular, in a replication where `tau` is large the posterior is wide
but the sandwich need not be. Both are computed and stored
(`.lo/.hi` are the sandwich, `.lo_post/.hi_post` the posterior
quantiles; the Monte Carlo keeps both and `summarize_montecarlo` reports
`coverage`/`avg_len` for the first and `coverage_post`/`avg_len_post`
for the second), so the price is measured rather than assumed — part C
of the sensitivity script.

---

## 4. Smaller items, listed for completeness

**Detrending trims `p` observations.** FMAR zero-pad the first `p`
detrended observations to preserve the sample length; this port drops
them instead (`utils/var_deterministic_trend.m`). Cleaner, at the cost
of `p` observations, and documented in that file's header.

**`tau^2` draws are clipped.** The horseshoe has Cauchy tails; draws are
clipped to `[1e-10, 1e8]` (`cfg.gibbs.tau2_min/max`), and every clip is
counted in `diag.n_tau_clip` so truncation is visible rather than
silent. In the runs reported here the count is ~0 per replication; if it
ever becomes large the reported `tau` summaries are truncated summaries
and must be labelled as such.

**Mixing of `log tau` in the pooled sampler.** The pooled estimator's
`log tau` path is updated by Metropolis, and it mixes slowly: effective
sample sizes of order 10–70 out of 700 retained draws, even with the
checkerboard sweep, the joint level move and five sweeps per Gibbs
iteration. This is the standard centred-parameterisation pathology of a
hierarchical smoothing prior. It matters much less than it looks,
because the reported quantity is not `tau` but the IRF, whose posterior
mean is **Rao-Blackwellised**: the estimator averages the conditional
means `E[beta | tau, sigma2, y]` instead of the draws. The diagnostic
that is actually asserted (`tests/test_pooled_estimator.m`) is the Monte
Carlo standard error of the reported mean relative to the posterior
standard deviation, which stays below 0.15. The `tau` *summaries* from
the pooled estimator are averaged over hundreds of Monte Carlo
replications, which removes what remains. A non-centred
parameterisation would fix the mixing properly and is the obvious next
step if the `tau` path itself, rather than the IRF, becomes the object
of interest.

**Rao-Blackwellisation changes reported numbers.** Averaging the
conditional means rather than the draws targets the *same* posterior
mean with strictly less simulation noise, and it is what makes the
comparison against the closed-form FMAR baseline fair — otherwise the
adaptive estimators are charged for simulation noise that a closed form
never pays. It does change the numbers relative to results produced
before this option existed; `cfg.blp.point_estimate = 'draw_mean'`
reproduces those.

**The experiment grid is one-factor-at-a-time.** A full factorial over
misspecification strength, sample size, lag order and sparsity is 144
cells and out of reach. The grid varies one factor at a time around a
baseline, so it **cannot identify interactions** — for instance whether
the sample size at which adaptation starts to pay depends on how sparse
the misspecification is. Every report of the grid says so.
