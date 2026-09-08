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
