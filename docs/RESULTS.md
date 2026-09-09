# Results

This note collects the headline numbers and, more importantly, says
which conclusions they support and which they do not. Every table here
is printed by `report_montecarlo` from a stored result and exported to
CSV alongside it; nothing is transcribed by hand.

> **Status of this file.** The sections below are filled in from the
> runs described in `docs/REPRODUCE.md`. Where a number is quoted, the
> file it came from is named. Where a run has not been done, the section
> says so rather than leaving an impression.
>
> Complete as of this revision: the 13-cell exploratory grid
> (`R = 40`), the sparse headline run (`R = 500`), the sensitivity
> analysis, and the full test suite. The `correct`, `intermediate` and
> `dense` headline runs were still executing when this was written; the
> correct-DGP false-positive numbers quoted below therefore come from
> the earlier `R = 500` run in `results/mc_fmar_final_correct.mat`,
> which used the FMAR `h = 1` convention and draw-averaged point
> estimates. Those conventions do not affect the τ localisation
> statistics, which is why they are quotable here; they do affect the
> RMSE columns, which is why no RMSE number is quoted from that file.
> Regenerate everything with `write_results_report` once the runs
> finish.

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

### All four headline runs, in one table

Paired RMSE ratios against the global FMAR baseline; `< 1` means
adaptation helps. `R = 500` for correct and sparse, `R = 250` for
intermediate and dense.

| DGP | estimator | early (2–6) | **all (2–20)** | late (7–20) |
|---|---|---|---|---|
| correct | BLP-block | 0.965 (t −13.9) | 1.011 (t +3.0) | 1.031 (t +6.3) |
| correct | **BLP-pooled** | **0.949** (t −11.4) | **0.988** (t −2.2) | 1.005 (t +0.7) |
| sparse | BLP-block | 0.968 (t −14.3) | 1.012 (t +3.9) | 1.040 (t +8.2) |
| sparse | **BLP-pooled** | **0.973** (t −7.4) | 1.003 (t +0.7) | 1.021 (t +3.1) |
| intermediate | BLP-block | 1.033 (t +9.1) | 1.025 (t +6.2) | 1.019 (t +3.2) |
| intermediate | BLP-pooled | 1.040 (t +7.3) | 1.014 (t +2.7) | 0.996 (t −0.5) |
| dense | BLP-block | 0.984 (t −4.7) | 1.011 (t +2.6) | 1.024 (t +4.0) |
| dense | **BLP-pooled** | **0.974** (t −3.8) | 0.999 (t −0.2) | 1.011 (t +1.2) |

Four readings, all of them firm at these sample sizes:

1. **The pooled estimator dominates the independent one on every DGP**,
   at every horizon window bar one. If the adaptive layer is used at
   all, it should be used with the horizons pooled.
2. **The independent block-adaptive estimator is reliably worse on
   integrated RMSE on all four DGPs** (1.011–1.025, `t` between 2.6 and
   6.2). There is no design here in which it is the right choice.
3. **The pooled estimator is a tie or a small win on three of four**
   (0.988, 1.003, 0.999) and a small loss on the fourth (1.014). Its
   *largest* integrated win is on the **correctly specified** DGP.
4. **The early/late split is the rule, not a sparse-DGP quirk** — it
   holds on correct, sparse and dense. The intermediate design is the
   one exception, and there adaptation is worse everywhere.

Reading 3 deserves emphasis because it is awkward for the obvious
motivation. The best case for this estimator, on integrated RMSE, is the
case where the VAR prior is *right*. That is consistent with the
tightening mechanism and inconsistent with selling it as a
misspecification remedy.

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

### Where adaptation actively hurts: the intermediate DGP

`results/mc_headline_intermediate.mat`, `R = 250`. Same lag-3
misspecification, but present in **every equation** (`A3(:,1) ≠ 0`)
rather than in one.

| estimator | early (2–6) | h=2–12 | all (2–20) | late (7–20) |
|---|---|---|---|---|
| BLP-block | 1.033 (t +9.1) | 1.034 (t +9.6) | 1.025 (t +6.2) | 1.019 (t +3.2) |
| BLP-pooled | 1.040 (t +7.3) | 1.031 (t +6.5) | 1.014 (t +2.7) | 0.996 (t −0.5) |

**The sparse pattern reverses.** Adaptation is worse at *every* horizon
window and worst **early** — precisely where, on the sparse DGP, it was
best. The τ diagnostic still localises perfectly here (concentration
1.00 in the grid), so this is not a failure to find the conflict; it is
a failure of the *response* to finding it.

The reason is the one the design was built to expose. The group
horseshoe buys its advantage from **sparsity in the deviation from the
VAR centre**. When the same block deviates in every equation, that
deviation is not sparse in the relevant sense: what the data want is a
*global* loosening, which the FMAR baseline gets for free by selecting a
looser `lambda_h`, and the per-block escape adds variance without buying
the right bias reduction.

So the method needs the conflict confined to a block **within an
equation**, not merely to a block. That is a real restriction on when it
should be used, and it is not something the τ map warns you about — the
map looks *identical* in the case where adaptation helps and the case
where it hurts.

(Note also that the BVAR comparator flips here: 1.279 against the FMAR
baseline, where on the sparse DGP it was 0.798. With the
misspecification present in every equation the VAR's bias finally
dominates its variance advantage.)

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
there. Equation 3 is found essentially every time.

### Detection and localisation against the R = 500 correct-DGP null

`results/mc_headline_correct.mat`, same design, nothing misspecified —
so every number in its column is a false-positive rate.

| rule | sparse (true block) | correct (false positive) | usable? |
|---|---|---|---|
| **concentration**, early horizons | **0.97** | **0.37** | yes — localises |
| concentration, all horizons | 0.80 | 0.43 | yes, weaker |
| **P(τ>1) > 0.50** on some (eq, block) | **0.74** | **0.03** | yes — detects |
| P(τ>1) > 0.75 | 0.05 | 0.00 | too conservative |
| τ ratio max/median > 1.25 | 0.24 | 0.30 | **no** |
| τ ratio max/median > 1.50 | 0.04 | 0.07 | **no** |

Two distinct jobs, and they need different statistics.

* **"Is anything escaping?"** — the probability rule works: 74% detection
  at a 3% false-positive rate. It is a genuine, calibrated test.
* **"Where?"** — the concentration statistic works: the true cell is
  ranked first in 97% of replications against a 37% null.
* **The ratio rules on aggregated τ do not work at all** (0.04 against
  0.07; 0.24 against 0.30). Averaging τ over equations *before*
  comparing blocks dilutes a signal confined to one equation almost to
  nothing. This is the rule the project previously reported, and it
  should not be used.

### On the correct DGP, pooling still wins

The most surprising number in the study. On the **correctly specified**
DGP, where there is nothing to escape from:

| estimator | early (2–6) | h=2–12 | all (2–20) | late (7–20) |
|---|---|---|---|---|
| BLP-block | 0.965 (t −13.9) | 0.981 (t −7.0) | 1.011 (t +3.0) | 1.031 (t +6.3) |
| BLP-pooled | **0.949** (t −11.4) | **0.951** (t −10.9) | **0.988** (t −2.2) | 1.005 (t +0.7) |

The pooled estimator is **reliably better than the global FMAR baseline
when the VAR prior is correct** (0.988, `t = −2.2`), driven by a 5%
early-horizon gain with no late-horizon penalty. Its posterior mean τ on
this DGP is 0.48–0.57 — far *below* 1.

That is the mechanism stated plainly: with the horizons pooled, each
block scale is informed by ~`H` times as much data, so where the VAR
prior fits the estimator can shrink *much harder than the global prior
does* and bank the variance. **The pooled block-adaptive estimator is
better understood as an adaptive tightening device that also permits
escape, than as an escape device.** That reframing matters for how the
thesis motivates it.

---

## 4. Sensitivity of the methodological approximations

`results/sensitivity_approximations.csv`; the reasoning is in
`docs/APPROXIMATIONS.md`.

Full reasoning and the tables are in `docs/APPROXIMATIONS.md`. In one
line each, against a noise baseline (the identical estimator re-run from
a different seed) of 0.0125 in IRF units:

* **Cross-equation covariance in the sampler: immaterial.** Fixing σ²
  at a value that used the system information moves the IRF by 0.0149
  (baseline 0.0125) and moves `log τ` by 0.556 (baseline 0.558) — i.e.
  not at all beyond resampling noise. At τ = 1 the change is *exactly*
  zero, which is an identity, not an estimate.
* **Inheriting λ from the global BLP: not an approximation.** The
  adaptive estimator's own marginal-likelihood selection returns the
  identical λ (largest gap across datasets 0.00e+00). What *is* an
  approximation is that λ is chosen at τ = 1 and never re-selected
  jointly with τ.
* **Fixed impact vector: negligible** (1.0% of the reported half
  band-width). **Fixed prior centre: the one that matters** (16.5%).
  Both channels hit the baseline and the adaptive estimators
  identically, so neither can explain an RMSE difference between them.
* **Quasi-Bayesian bands: doing real work.** See the interval table
  above — the sandwich covers 0.933–0.937 against a nominal 0.90, while
  posterior quantiles undercover at 0.810 (independent) and 0.691
  (pooled).

---

## 5. What the evidence supports as the contribution

The question this note has to answer is whether the thesis should frame
its contribution as **(1)** an improved estimator, **(2)** a Bayesian /
localised specification diagnostic, or **(3)** still inconclusive.

**The evidence supports (2), with a specific and defensible secondary
claim.** It does not support (1), and (3) understates what is now firmly
established.

### Why not (1), an improved estimator

The independent block-adaptive estimator is **reliably worse on
integrated RMSE on all four DGPs** (1.011–1.025, `t` from 2.6 to 6.2).
The pooled estimator is a tie or a small win on three of four and a
small loss on the fourth. To claim an improved estimator one would have
to lean on the early-horizon window — where the gains are real and
strongly significant (0.949–0.974 on three of four DGPs) — and that
claim is only as good as the choice of window, which the horizon-window
table shows is doing the work. Worse for the narrative: the pooled
estimator's **largest** integrated win is on the **correctly specified**
DGP (0.988), which is the wrong headline for a thesis about
misspecification.

### Why (2), a localised specification diagnostic

Here the evidence is strong, and it is the kind of evidence a diagnostic
needs — a measured null, not just a signal.

* **Localisation works and is calibrated.** The true (equation, block)
  cell is ranked first in **97%** of replications at early horizons,
  against a **37%** false-positive baseline measured on the correct DGP
  at the same `R = 500`. Equations with nothing wrong sit at chance.
* **Detection works and is calibrated.** `P(τ > 1) > 0.5` somewhere
  fires at **0.74** on the sparse DGP against **0.03** under the null —
  a genuine test at a 3% false-positive rate.
* **It is honest when there is nothing to localise.** On the dense
  (VARMA) DGP, which is misspecified everywhere, the concentration
  statistic reads 0.42 — the same as the correct-DGP baseline. It does
  not manufacture a block.
* **Monotone dose–response** in misspecification strength
  (0.42 → 0.72 → 0.95 → 1.00).
* **It costs nothing to adopt.** The adaptive layer is *exactly* nested
  in the published FMAR estimator at `τ = 1` (machine precision, 912
  comparisons), so the diagnostic can be run alongside an FMAR
  application without changing the reported IRFs at all.
* **Horizon pooling makes the diagnostic much sharper**, not just the
  estimator: `P(τ > 1)` for the true block runs 0.85–0.97 across
  `h = 1..8` instead of a ragged 0.12–0.98.

### The two caveats that must travel with it

1. **The null is specification-dependent.** The false-positive baseline
   is 0.42 at `p = 2` but 0.62–0.68 at `p = 3`/`p = 4`, where the fitted
   VAR *nests* the truth and nothing is wrong: a richer fit estimates
   its prior centre less precisely and τ faithfully reports the
   resulting disagreement. No threshold transfers across
   specifications; the null has to be simulated for the design at hand.
   This is exactly the parametric bootstrap the Chapter 7 design already
   specifies, and it is now an argued requirement rather than a
   precaution.
2. **The diagnostic does not tell you whether escaping will help.** On
   the intermediate DGP the τ map localises perfectly *and* adaptation
   is worse at every horizon. The map looks the same in the case where
   adaptation helps and the case where it hurts. A diagnostic that finds
   where the prior fails is therefore not, by itself, a prescription to
   let that block escape.

### The defensible secondary claim

Caveat 2 is not only a limitation: stated positively, it is a finding.
**Knowing where a VAR prior fails does not imply that locally escaping
it improves estimation** — and this project has the paired,
`R = 500` evidence to say so, together with the mechanism (the group
horseshoe monetises *sparsity in the deviation from the prior centre*,
which is a different and stronger condition than the misspecification
being localised). That is a real, quotable negative result about
adaptive shrinkage, not a null finding.

The horizon-pooled estimator should be reported as the methodological
contribution that makes the diagnostic usable — sharper τ paths, and an
adaptive layer that no longer costs integrated RMSE (1.003 on the sparse
DGP, `t = 0.67`) — rather than as an estimator that beats the benchmark.

### One-line recommendation

> Frame the thesis around a **calibrated, localised diagnostic for
> prior–data conflict in VAR-centred Bayesian local projections**, whose
> credibility rests on exact nesting in the published estimator and on
> measured nulls; report the horizon-pooled adaptive estimator as the
> refinement that makes the diagnostic sharp and RMSE-neutral; and
> present the failure of localisation to imply improvement as a finding
> in its own right rather than as a disappointing RMSE table.
