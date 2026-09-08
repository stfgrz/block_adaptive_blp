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
so the two cannot be confused. In particular the *nominal* chance level
`1/G` is a lower bound on the null of the concentration statistic, not
the null itself: that statistic is a maximum over `K*G` cells, and in
this design block 1 (the shock variable) is mildly favoured even under
correct specification. The correct-DGP run is the null.

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

### Why horizon pooling changes the picture

The independent estimator identifies each `tau_{i,g,h}` from `p_g = p`
coefficient deviations at a single horizon. The pooled estimator ties
the `H` horizons of one block together, so the same scale is informed by
roughly `H` times as much data. The consequence is visible in both
directions at once, on a single sparse dataset (`RUN_FMAR_DEMO`, `y_3`
equation, `h = 1..8`):

| | h=1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| `tau`, block 1, independent | 0.49 | 2.92 | 1.64 | 1.90 | 1.79 | 1.80 | 1.43 | 1.41 |
| `tau`, block 1, **pooled** | 1.57 | 1.70 | 1.69 | 1.65 | 1.64 | 1.60 | 1.53 | 1.43 |
| `tau`, block 3, independent | 0.52 | 0.55 | 0.55 | 0.67 | 0.62 | 0.64 | 0.57 | 0.70 |
| `tau`, block 3, **pooled** | 0.20 | 0.18 | 0.15 | 0.14 | 0.13 | 0.12 | 0.12 | 0.12 |
| `P(tau>1)`, block 1, independent | 0.12 | 0.98 | 0.74 | 0.77 | 0.78 | 0.71 | 0.57 | 0.54 |
| `P(tau>1)`, block 1, **pooled** | 0.85 | 0.96 | 0.96 | 0.96 | 0.97 | 0.94 | 0.91 | 0.83 |

Pooling does not merely smooth. It makes the escape of the conflicted
block *sustained and confident* (`P(tau > 1)` between 0.83 and 0.97 at
every horizon instead of a ragged 0.12–0.98) **and** it shrinks the
non-conflicted blocks far harder than the global prior does (0.12–0.20
instead of 0.5–0.7). That second half is where the variance saving comes
from, and it is why the pooled estimator can beat the global baseline
even on the correctly specified DGP, where there is nothing to escape
from at all.

### The null depends on the fitted specification

The `p = 3` grid cell is the sharpest test of "is a large `tau` evidence
of misspecification?". There the fitted VAR(3) **nests** the sparse
design's truth, so nothing is misspecified and every flag is a false
positive. Yet the concentration statistic reads **0.68**, against
**0.42** on the correct DGP at `p = 2`.

Nothing is wrong with the diagnostic; the reading is that a richer
fitted VAR estimates its prior centre less precisely at `T = 200`, and
the LP disagrees most with that centre exactly where the extra
coefficient lives. The diagnostic is detecting **prior-data
disagreement**, faithfully — and at `p = 3` the disagreement is
sampling noise in the centre, not structure.

Two consequences, and they matter more than any RMSE number here.

1. **A fixed threshold on `tau` cannot be calibrated once and reused.**
   The false-positive rate of any rule depends on `K`, `T`, the fitted
   lag order and the persistence of the system, because all of those
   determine how sharply the prior centre is estimated.
2. **The null has to be simulated for the specification at hand.** That
   is precisely what the Chapter 7 design already calls for (its "Leg 2"
   parametric bootstrap under the fitted BVAR), and this cell is the
   simulation-side argument for why that leg is not optional.

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

| case | axis | r_blk | r_pool | rE_blk | rE_pool | conc | FP? | verdict (block) |
|---|---|---|---|---|---|---|---|---|
| sparsity_correct | sparsity | 0.984 | **0.954** | 0.961 | 0.940 | 0.42 | yes | early horizons |
| sparsity_sparse | sparsity | 0.993 | **0.975** | 0.977 | 0.982 | 0.95 | no | early horizons |
| sparsity_interm | sparsity | 1.032 | 1.033 | 1.037 | 1.040 | 1.00 | no | diagnostic only |
| sparsity_dense | sparsity | 0.996 | 0.982 | 0.991 | 0.979 | 0.42 | n/l | no signal |
| strength_a015 | strength | 0.990 | 0.973 | 0.981 | 0.981 | 0.72 | no | diagnostic only |
| strength_a045 | strength | 0.983 | 0.958 | 0.962 | 0.960 | 1.00 | no | early horizons |
| strength_a060 | strength | 0.978 | **0.956** | 0.954 | 0.949 | 1.00 | no | improves RMSE |
| strength_dense_half | strength | 0.991 | 0.976 | 0.981 | 0.973 | 0.45 | n/l | no signal |
| T100 | sample size | 1.003 | 0.986 | 0.995 | 0.997 | 0.88 | no | diagnostic only |
| T400 | sample size | 0.987 | **0.958** | 0.950 | 0.940 | 1.00 | no | early horizons |
| p1 | lag order | 1.027 | 1.011 | 1.010 | 1.018 | 1.00 | no | diagnostic only |
| p3 | lag order | 0.985 | 0.953 | 0.979 | 0.965 | 0.68 | **yes** | early horizons |
| p4 | lag order | 0.988 | 0.964 | 0.986 | 0.977 | 0.62 | **yes** | no signal |

`r_blk`, `r_pool` = adaptive / global integrated RMSE over `h = 2..H`;
`rE_*` the same over the early window; `< 1` means adaptation helps.
`conc` = the concentration statistic (max over equation x block of the
early-horizon win rate). `FP?`: **yes** = nothing is misspecified in
that cell, so `conc` is a false-positive baseline; **n/l** = the cell IS
misspecified but has no localisable block, so `conc` is neither; **no** =
there is a true block and `conc` is a detection rate.

### What the grid shows

**1. Horizon pooling dominates.** The pooled estimator beats the
independent one in **12 of 13 cells** and beats the *global baseline* in
**12 of 13** (ratios 0.953–0.986). The single exception is the
intermediate DGP, where both adaptive estimators lose (1.03). At `R = 40`
no individual cell is decisive, but 12 of 13 falling the same way is not
noise.

**2. Adaptation pays more as the sample grows, not less.**
`T = 100 → 0.986`, `T = 200 → 0.975`, `T = 400 → 0.958` (pooled). If the
binding constraint were bias, the ordering would run the other way; that
it runs this way says the binding constraint is **variance** — exactly
the diagnosis the bias–variance table gives.

**3. Adaptation pays more as the conflict grows.** `a31 = −0.15 → 0.973`,
`−0.30 → 0.975`, `−0.45 → 0.958`, `−0.60 → 0.956`, and only the strongest
cell earns an unqualified "improves RMSE" verdict.

**4. Where it does NOT help: the intermediate design.** When the *same*
block is wrong in *every* equation (`A3(:,1) ≠ 0`), both adaptive
estimators lose by ~3%, even though the diagnostic localises perfectly
(conc = 1.00). Block adaptation needs the conflict confined to a block
*within an equation*; when it is spread across equations, a global
loosening would serve better and the per-block escape is pure cost. This
is the sharpest negative result in the study and it is a genuine
limitation of the method, not of the implementation.

**5. Where it does not help either: `p = 1`.** With a fitted VAR(1) the
prior is wrong in many places at once; the ratio is 1.027 (block) /
1.011 (pooled). Same lesson as (4) from a different direction.

**6. The diagnostic has a clean dose–response, and a
specification-dependent null.** Detection rises monotonically with
conflict strength (0.42 → 0.72 → 0.95 → 1.00) and correctly reports
nothing on the dense DGP (0.42, the same as the correct-DGP baseline).
But the false-positive baseline is **not a constant**: 0.42 at `p = 2`,
0.62–0.68 at `p = 3`/`p = 4` where the fitted VAR nests the truth. See
"The null depends on the fitted specification" above.

---

## 3. Headline Monte Carlo runs

`results/mc_headline_<dgp>.mat`, `R = 500` (sparse, correct) and `R = 250`
(intermediate, dense), at `K = 3`, `p = 2`, `T = 200`, `H = 20`, FMAR
mode with `h1_mode = 'lp'`.

### The single most informative table in the study

`results/mc_headline_sparse.mat`, `R = 500`, sparse DGP. Paired
comparison against the global FMAR baseline — every estimator sees the
same simulated sample in every replication, so the *difference* is
estimated far more precisely than either level (the per-replication MSEs
correlate at 0.96–0.98, which buys roughly an order of magnitude of
precision over an unpaired comparison). `< 1` means adaptation helps;
`t` is the Monte Carlo t-statistic of the paired MSE difference.

| estimator | early (h=2–6) | h=2–12 | **all (h=2–20)** | late (h=7–20) |
|---|---|---|---|---|
| BLP-block (independent τ) | **0.968** (t −14.3) | **0.988** (t −5.7) | 1.012 (t +3.9) | 1.040 (t +8.2) |
| BLP-pooled (τ pooled over h) | **0.973** (t −7.4) | **0.973** (t −8.3) | 1.003 (t **+0.7**) | 1.021 (t +3.1) |

**This answers the question the project was stuck on.** Block adaptation
does not "lose on average" in any interesting sense. It **reliably wins
at early horizons** (3% lower RMSE, `t = −14`) and **reliably loses at
late ones** (4% higher, `t = +8`), and the integrated verdict is decided
by nothing more than how many late horizons the average happens to
include. At `H = 20` the fourteen late horizons outvote the five early
ones; at `H = 12` they do not, and the sign flips. A single integrated
RMSE reported without this decomposition invites reading an arbitrary
choice of `H` as a property of the estimator.

**And it answers what horizon pooling buys.** Pooling *preserves* the
early-horizon correction (0.968 → 0.973, statistically the same) while
*halving* the late-horizon penalty (1.040 → 1.021). That takes the
integrated ratio from a reliably-worse 1.012 (`t = 3.9`) to
indistinguishable from the baseline: **1.003, `t = 0.67`**. The pooled
estimator buys the early-horizon bias reduction essentially for free.

### Bias and variance, sparse DGP, R = 500

| estimator | h=2–20 \|bias\| | var | h=2–6 \|bias\| | var | h=7–20 \|bias\| | var |
|---|---|---|---|---|---|---|
| BLP-FMAR | 0.0201 | 0.0040 | 0.0342 | 0.0053 | 0.0151 | 0.0036 |
| BLP-block | 0.0194 | 0.0042 | **0.0326** | 0.0051 | 0.0147 | 0.0039 |
| BLP-pooled | 0.0195 | 0.0041 | 0.0354 | **0.0049** | 0.0138 | 0.0038 |

The mechanism is exactly the one hypothesised: **bias falls where the
VAR prior is wrong** (in `y_3`, the misspecified response, the bias goes
0.0270 → 0.0244 → 0.0243) and **variance rises where it is not** (late
horizons, 0.0036 → 0.0039). Pooling attacks the second term without
giving up the first.

### Intervals, sparse DGP, R = 500 (nominal 90%)

| estimator | coverage (sandwich) | length | coverage (posterior quantiles) | length |
|---|---|---|---|---|
| LP | 0.836 | 0.278 | — | — |
| BVAR | 0.601 | 0.116 | — | — |
| BLP-FMAR | 0.937 | 0.284 | — | — |
| BLP-block | 0.933 | 0.284 | 0.810 | 0.181 |
| BLP-pooled | 0.934 | 0.286 | 0.691 | 0.138 |

This is the empirical price of approximation 3. The quasi-Bayesian
Newey–West sandwich covers at 0.93–0.94 against a nominal 0.90 (mildly
conservative); the **posterior-quantile bands undercover badly** — 0.81
for the independent estimator and **0.69** for the pooled one, whose
posterior is sharper precisely because pooling concentrates it. Keeping
the sandwich as the primary interval is doing real work, and the pooled
estimator's posterior spread should not be reported as a credible band.

### τ localisation, sparse DGP, R = 500

Concentration statistic (max over equation × block of the win rate):

| window | value | at | reading |
|---|---|---|---|
| all horizons | 0.80 | (eq 3, block 1) | the true cell |
| early (h ≤ 6) | **0.97** | (eq 3, block 1) | the true cell |

Per-equation × block early-horizon winner map:

| | block 1 | block 2 | block 3 |
|---|---|---|---|
| eq 1 | 0.31 | 0.37 | 0.31 |
| eq 2 | 0.35 | 0.33 | 0.32 |
| **eq 3** | **0.97** | 0.02 | 0.00 |

Equations 1 and 2 sit at chance, as they should — nothing is wrong
there. Equation 3 is found essentially every time. The aggregate flag
rules, by contrast, fire at 0.04–0.24 and are useless: averaging τ over
equations before comparing blocks dilutes a signal confined to one
equation almost to nothing.

---

## 4. Sensitivity of the methodological approximations

`results/sensitivity_approximations.csv`; the reasoning is in
`docs/APPROXIMATIONS.md`.

<!-- RESULTS-SENSITIVITY -->
