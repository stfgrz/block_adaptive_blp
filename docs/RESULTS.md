# Results

This note collects the headline numbers and, more importantly, says
which conclusions they support and which they do not. Every table here
is printed by `report_montecarlo` from a stored result and exported to
CSV alongside it; nothing is transcribed by hand.

> **Status of this file.** The sections below are filled in from the
> runs described in `docs/REPRODUCE.md`. Where a number is quoted, the
> file it came from is named. Where a run has not been done, the section
> says so rather than leaving an impression.

---

## How to read the numbers

**Which RMSE.** Quote **`irmse_h2`**, the integrated RMSE over
`h = 2..H`. The legacy `h = 1..H` version is kept for compatibility but
is not a fair estimator comparison under the FMAR `h = 1` convention
(see the README, §3a).

**Monte Carlo precision.** The standard error of an estimated RMSE is
about `rmse / sqrt(2R)`: **3.2%** at `R = 500`, **5%** at `R = 200`,
**11%** at `R = 40`. Differences *between* estimators are more precise
than that, because every estimator sees the same simulated samples in
every replication (common random numbers), but no difference smaller
than a few percent should be read as real at `R = 40`.

**Bias and variance.** `variance` is the Monte Carlo variance of the
estimator across replications, with the `1/R` convention that makes
`MSE = bias² + variance` exact (asserted in `summarize_montecarlo`).
This is the table that answers the question the project actually asks:
*does adaptation reduce bias where the VAR prior is wrong, and does it
pay for that with variance elsewhere?*

**Detection numbers are meaningless alone.** A `tau` localisation rate
on a misspecified DGP must be read next to the rate at which the same
rule fires on the **correctly specified** DGP, where every flag is a
false positive by construction. Both are reported, always, and the
diagnostic struct is stamped `false_positive = true` in the second case
so the two cannot be confused.

**`tau` is prior–data disagreement, not misspecification.** A large
`tau_{i,g,h}` says the horizon-`h` LP coefficients of block `g` in
equation `i` disagree with the VAR-implied centre. Structural
misspecification is one cause; sampling noise in the estimated centre,
a badly chosen `lambda_h`, and weak identification of that particular
horizon are others. The `p = 3` and `p = 4` grid cells exist precisely
to probe this: there the fitted VAR *nests* the truth of the sparse
DGP, so any escape is a false positive.

---

## What the adaptive layer actually does

Before any RMSE table, it helps to know what `tau` does in practice,
because "escape from the VAR prior" turns out to be only half of it.

On the `R = 500` sparse run, the posterior mean of `tau` in the
misspecified equation `y_3`, by block and horizon, is

| block | h=1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **1** (misspecified) | 0.28 | **2.16** | **1.49** | **1.57** | **1.58** | **1.52** | 1.42 | 1.33 | 1.23 | 1.17 | 1.14 | 1.11 |
| 2 | 0.28 | 0.59 | 0.72 | 0.77 | 0.81 | 0.87 | 0.93 | 0.99 | 1.04 | 1.06 | 1.08 | 1.09 |
| 3 | 0.28 | 0.53 | 0.64 | 0.70 | 0.75 | 0.81 | 0.84 | 0.88 | 0.91 | 0.96 | 0.98 | 1.01 |

and on the correctly specified run the same table has **no** contrast
between blocks (0.29 / 0.29 / 0.30 at `h = 1`, rising to ~1.05 at
`h = 12` in all three).

Three things follow, and they shape how everything below should be read.

**1. The signal is a WITHIN-equation, WITHIN-horizon contrast, not a
level.** `tau` rises with the horizon under *both* DGPs — that is a
property of how much the horizon-`h` regression has to say about the
deviation, not of misspecification. Comparing block 1 against blocks 2
and 3 *at the same horizon in the same equation* removes that common
profile, which is exactly what the argmax/concentration diagnostic does
and why it works while threshold rules on the level of `tau` do not.

**2. The adaptive layer TIGHTENS as often as it escapes.** Averaged over
all (equation, block, horizon) cells the posterior mean of `tau` is
**0.965** on the sparse DGP and **0.964** on the correct one. The
group-horseshoe is therefore not only an escape device: where the VAR
prior fits, it shrinks *harder* than the global prior, which buys
variance. Only the conflicted block goes above 1. This is why adaptation
can help even on the correctly specified DGP, which the exploratory grid
does show, and it is a more interesting mechanism than "let the wrong
block out".

**3. The `h = 1` value 0.28 was an artefact of the unfair comparison.**
Under the FMAR convention the adaptive estimators inherited the *BVAR's*
Minnesota tightness at `h = 1` (0.27) for what is an LP regression whose
own marginal-likelihood tightness is about 0.10. The prior was far too
loose, and `tau` shrank to compensate. With `cfg.fmar.h1_mode = 'lp'`
the `h = 1` regression gets an LP tightness and this collapses; it is
also why the legacy `h = 1..H` integrated RMSE moves with the `h = 1`
convention while `irmse_h2` does not.

---

## 1. Established: exact nesting

`tests/test_fmar_nesting_exact.m`. With `tau = 1`, both adaptive
estimators reproduce the FMAR closed-form posterior mean over **912
(equation, horizon) comparisons** — four DGPs, two seeds, every
equation, every horizon, both `h = 1` conventions, both adaptive
estimators.

| quantity | value |
|---|---|
| worst deviation | **3.3e-16** |
| tolerance | 1e-10 |

This is a tolerance-free check: with `tau` fixed and the conjugate
`sigma²`-scaled prior, the sampler's *conditional* posterior mean is a
deterministic function of the data, so the comparison is of two closed
forms rather than of a simulation against a closed form. The older
sampling-based test (`test_fmar_nesting.m`) is kept and now also reports
a maximum error of 0.0000, because the reported point estimate is
Rao-Blackwellised.

---

## 2. Exploratory grid — where does adaptation help?

`results/grid/grid_summary.csv`, produced by
`scripts/run_grid_parallel.sh` and `collect_grid_results`.
Thirteen cells at `R = 40`, `H = 12`.

<!-- RESULTS-GRID -->

---

## 3. Headline Monte Carlo runs

`results/mc_final_<dgp>.mat`, `R = 500` (sparse, correct) and `R = 250`
(intermediate, dense), at `K = 3`, `p = 2`, `T = 200`, `H = 20`, FMAR
mode with `h1_mode = 'lp'`.

<!-- RESULTS-FINAL -->

---

## 4. Sensitivity of the methodological approximations

`results/sensitivity_approximations.csv`; the reasoning is in
`docs/APPROXIMATIONS.md`.

<!-- RESULTS-SENSITIVITY -->
