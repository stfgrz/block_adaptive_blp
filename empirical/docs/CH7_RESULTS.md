# Chapter 7 — results of the euro-area application

Run of 2026-09-12 (MATLAB R2025b). Every number below is read from a
file in `results/`; the file is named next to each table. How the
pipeline works and why three design decisions were taken on the real data
is in `README_EMPIRICAL.md` (Status 1–3); the pre-registered design and
expectations are in `CH7_DESIGN.md`.

**Baseline** (tag `p12`): `T = 240` (2000m1–2019m12), `K = 5`
(`mps, i1y, ip, hicp, stoxx`), `p = 12`, `H = 48`, FMAR mode with
`h1_mode = 'lp'`, the prior-scale floor on, the surprise block held at
`tau = 1`, chains of 500 + 1500 draws, early-horizon window `h = 2..12`.
Variants: `p12_nofloor` (the published FMAR scale, every block free) and
`p12_allblocks` (floor on, every block free). The level system without
the surprise is `*_levels_p12*`.

---

## 1. Data, shock and instrument relevance (T7.1)

| # | variable | source | transform | mean | sd |
|---|---|---|---|---|---|
| 1 | `mps` | EA-EMPD GC monetary-event-window 1-year OIS surprise, monthly sum | pp | 0.001 | 0.044 |
| 2 | `i1y` | 12-month Euribor, monthly average (ECB `FM.M.U2.EUR.RT.MM.EURIBOR1YD_.HSTA`) | percent | 1.94 | 1.73 |
| 3 | `ip` | Eurostat `sts_inpr_m M.PRD.B-D.SCA.I21.EA20` | 100 log | 454.9 | 4.7 |
| 4 | `hicp` | Eurostat `prc_hicp_midx M.I15.CP00.EA` (NSA) | 100 log | 451.5 | 10.0 |
| 5 | `stoxx` | EURO STOXX 50, monthly average (ECB `FM.M.U2.EUR.DS.EI.DJES50I.HSTA`) | 100 log | 806.5 | 20.8 |

The shock series reproduces the reference statistics exactly (221
meeting events, monthly sd 4.29 bp, 24 zero months; `build_shock_series`).

**Relevance.** The instrument enters recursively (ordered first), so its
relevance is the regression of the 1-year-rate innovation on the surprise
innovation from the BVAR residuals — the slope that the 25 bp
normalisation divides by (`ea_relevance`, stored as `rel` in
`empirical_<tag>.mat` and `dose_response_<tag>.mat`):

| fitted `p` | impact of a 1 pp surprise innovation on `i1y` (pp) | robust t | F | 25 bp scale `k` |
|---|---|---|---|---|
| 2 | 0.071 | 0.22 | 0.05 | 3.50 |
| 6 | — | 0.52 | 0.27 | — |
| 12 | 0.236 | 0.77 | 0.59 | 1.06 |

The naive same-month regression of the rate change on the surprise gives
0.03 (t = 0.07); the two-month change gives −0.06. The largest surprises
of the sample sit in November–December 2008 (positive surprises of 14
and 9 bp against 90 bp monthly falls in the rate), and the monthly
*average* rate spreads any mid-month event over two months. Only the
information-type surprises (`mps_gc_info`, same-sign stock move) show a
positive first stage (0.92, t = 2.6). **The surprise is a weak instrument
for this indicator.** The τ diagnostic below does not depend on this (it
is about the dynamic specification of the system); the IRFs in Section 4
must be read with it.

---

## 2. The tightness path, and why the design changed

`lambda_h` selected by the FMAR marginal likelihood at `p = 12`
(`empirical_<tag>.mat`, `blpf.lambda(1,:)`; `fig_lambda_<tag>.png`):

| | h = 1 | 2 | 6 | 12 | 20 | 24 | 30 | 40 | 48 | horizons with a bimodal objective |
|---|---|---|---|---|---|---|---|---|---|---|
| published scale (`_nofloor`) | 0.008 | 0.021 | 0.030 | 0.017 | 0.008 | 0.008 | 0.420 | 0.428 | 0.255 | 28 of 48 (h = 21..48) |
| with the floor (baseline) | 0.008 | 0.024 | 0.054 | 0.106 | 0.281 | 0.398 | 0.436 | 0.442 | 0.292 | 0 of 48 |
| level system, no surprise | 0.012 | 0.049 | 0.184 | 0.361 | 0.676 | 0.744 | 0.876 | 0.898 | 0.702 | 0 of 48 |

Under the published scale the objective is bimodal from `h = 21` and the
selected value flips between 0.01 and 0.4 across adjacent horizons
(0.014 at h = 28, 0.41 at 29–30, 0.010 at 31, 0.46 at 32), which produces
a vertical stripe in the τ surface and a jump of the FMAR IRF onto the LP.
The floor removes the bimodality and gives a smooth, increasing path. The
level system's own path is higher still: the surprise block continues to
drag the global tightness down by a factor of about three at `h = 12`
even with the floor. Every τ in Section 3 is measured relative to the
baseline path in the second row.

---

## 3. Leg 1 — where the data disagree with the VAR prior (F7.3)

Early-horizon mean posterior scale, `tau_bar(i, g) = mean_{h=2..12}
tau_{i,g,h}`, rows = equations, columns = blocks. The `mps` block is
held at 1 by design (no signal). Source: `tau_heatmap_p12.csv`.

**Independent estimator (BLP-block)**

| equation \ block | mps (held) | i1y | ip | hicp | stoxx |
|---|---|---|---|---|---|
| mps | 1 | 0.52 | 0.63 | 0.60 | 0.69 |
| i1y | 1 | 2.02 | 1.34 | 0.86 | 1.72 |
| ip | 1 | **5.33** | 3.21 | 1.95 | 1.31 |
| hicp | 1 | 2.41 | 1.82 | **9.45** | 2.29 |
| stoxx | 1 | 2.58 | 1.93 | 0.63 | 0.70 |

**Horizon-pooled estimator (BLP-pooled)**

| equation \ block | mps (held) | i1y | ip | hicp | stoxx |
|---|---|---|---|---|---|
| mps | 1 | 0.07 | 0.10 | 0.05 | 0.09 |
| i1y | 1 | 4.57 | 2.41 | 0.95 | 2.70 |
| ip | 1 | **7.57** | 4.39 | 2.07 | 1.81 |
| hicp | 1 | 3.27 | 2.34 | **10.05** | 2.97 |
| stoxx | 1 | 4.39 | 3.35 | 0.38 | 0.25 |

Posterior `P(tau > 1)` averaged over the window, pooled estimator: above
0.90 for every cell of the `ip` and `hicp` equations except `hicp ← stoxx`
(0.90) and for `i1y ← i1y`, `i1y ← ip`, `i1y ← stoxx` and
`stoxx ← i1y`; at most 0.41 for `i1y ← hicp`; zero in the `mps` equation
and for `stoxx ← hicp`, `stoxx ← stoxx`.

What the map says, before the null:

1. **The negative control is quiet.** In the surprise equation every level
   block sits below 1 (0.5–0.7 independent, 0.05–0.13 pooled). The
   surprise is unforecastable from lagged levels, and the VAR prior says
   so correctly.
2. **The strongest disagreement is the own-lag block of HICP** (9.4 /
   10.0, `P(tau > 1) = 1.00`). HICP is not seasonally adjusted; the
   horizon-`h` projection of prices on their own last twelve months
   disagrees with the VAR(12)'s iterated dynamics at every `h`, and does
   so at `p = 2`, `6` and `12` alike (Section 5).
3. **The interest-rate block escapes in the real, nominal and financial
   equations** (`ip ← i1y` 5.3 / 7.6, `hicp ← i1y` 2.4 / 3.3,
   `stoxx ← i1y` 2.6 / 4.4) and in its own equation (2.0 / 4.6): the
   propagation of the policy indicator is where the LP and the VAR-implied
   centre part company. The `ip` block escapes in the `ip` and `stoxx`
   equations.
4. **The adaptive layer also tightens.** `stoxx ← hicp`, `stoxx ← stoxx`
   and `i1y ← hicp` are held below the global prior (0.25–0.95 pooled):
   there the VAR centre is, if anything, too loose.
5. **Pooling sharpens the same picture** rather than changing it: every
   argmax is the same under both estimators (`i1y` block in the `ip` and
   `stoxx` equations, own block in `i1y` and `hicp`), and the pooled
   `P(tau > 1)` is near 0 or near 1 almost everywhere.

Whether these cells exceed what a correctly specified VAR(12) would
produce on data like these is the question of Section 7.

---

## 4. IRFs to a 25 bp surprise (F7.2)

Peak response within `h ≤ 24` (25 bp normalisation, `k = 1.06`; units:
pp for `i1y`, 100 × log for the rest); source `empirical_p12.mat`.

| response | BVAR | BLP-FMAR | BLP-block | BLP-pooled | LP-LA | h of FMAR peak |
|---|---|---|---|---|---|---|
| `i1y` | −0.25 | −1.76 | −1.91 | −1.89 | −1.34 | 24 |
| `ip` | −1.51 | −7.33 | −7.23 | −6.99 | −5.56 | 22 |
| `hicp` | −0.95 | −1.94 | −0.47 | −0.43 | −3.06 | 19 |
| `stoxx` | −2.0 | −27.4 | −35.0 | −35.1 | 1.3 | 24 |

Read with Section 1 in mind: the impact response the scale divides by
has a t statistic of 0.77, the 90% bands (`fig_irf_p12.png`) include
zero almost everywhere, and the sign of the rate response at
`h = 20–30` is negative for every LP-based estimator — the pattern of a
weakly identified shock, not a transmission estimate. The exhibit's job
here is methodological: the two adaptive estimators now sit on the FMAR
path wherever FMAR leaves the BVAR (`i1y`, `ip`, `stoxx`), which is what
holding the instrument block at 1 was meant to achieve; with that block
free they sat on the BVAR instead (`empirical_p12_allblocks.mat`).

---

## 5. Leg 3 — dose–response in the lag order

Early-horizon `tau_bar`, independent estimator, at `p = 2, 6, 12` on
identical data (`dose_response_p12.csv`; the pooled estimator gives the
same ordering with larger values at `p = 12`):

| cell | p = 2 | p = 6 | p = 12 |
|---|---|---|---|
| `hicp ← hicp` | 4.83 | 5.23 | 9.52 |
| `ip ← i1y` | 2.09 | 1.42 | 5.64 |
| `hicp ← i1y` | 2.37 | 1.58 | 2.47 |
| `stoxx ← i1y` | 1.73 | 1.47 | 2.81 |
| `ip ← hicp` | 2.04 | 1.79 | 1.93 |
| `ip ← ip` | 1.27 | 1.14 | 3.33 |
| median over the 20 free cells (block / pooled) | 1.00 / 0.69 | 1.31 / 1.21 | 1.61 / 2.31 |
| `lambda_h` at h = 2 / 6 / 12 | 0.064 / 0.146 / 0.356 | 0.031 / 0.092 / 0.224 | 0.024 / 0.054 / 0.106 |

The pre-registered prediction was a map that lights up at `p = 2` and
fades with `p`. The *levels* do the opposite for the strongest cells. The
reason is not the data but the metric: `tau` is measured relative to
`lambda_h`, and `lambda_h` is three to six times tighter at `p = 12` than
at `p = 2` (0.024 vs 0.064 at h = 2; 0.106 vs 0.356 at h = 12), because
a VAR with 61 regressors is priced by the marginal likelihood very
differently from one with 11. In effective units, `lambda_h * tau`, the
HICP own-block deviation is 1.7 at `p = 2` and 1.0 at `p = 12`, i.e. it
does fade. The comparable statistic across `p` is therefore the number of
cells above their own null q95 at each `p`, reported next.

**Read through the protocol** (`dose_response_protocol_p12.csv`, nulls at
`p = 2` and `p = 6` with `R = 200`, `p = 12` with `R = 500`): the number of
free cells above their own null q95, and which they are.

| fitted `p` | independent | pooled |
|---|---|---|
| 2 | 1: hicp ← hicp (4.83 vs 4.24) | 1: hicp ← hicp (4.24 vs 4.05) |
| 6 | 1: hicp ← hicp (5.23 vs 4.90) | 0: none |
| 12 | 2: ip ← i1y (5.64 vs 4.98); hicp ← hicp (9.51 vs 6.03) | 2: ip ← i1y (7.35 vs 6.27); hicp ← hicp (10.08 vs 7.64) |

The calibrated dose–response is flat, not decreasing: a knowingly too
short monthly VAR(2) does **not** produce widespread calibrated escapes,
because its null is just as wide as its real-data map (the null
equation-wide q95 at `p = 2` runs 3.2–6.3 against 4.5–7.3 at `p = 12`). The
HICP own block is the one cell that is above its null at every lag order;
the interest-rate block in the `ip` equation joins it only at `p = 12`.
The null `lambda_h` path moves with `p` exactly as the real-data one does
(median null `lambda_12`: 0.30 at `p = 2`, 0.13 at `p = 6`, 0.07 at
`p = 12`), which is why the levels in the table above cannot be compared
across `p` and the counts can.

The negative control holds at every `p`: the largest level-block `tau_bar`
in the surprise equation is 0.78 / 0.65 / 0.65 (independent) and
0.07 / 0.09 / 0.16 (pooled) at `p = 2 / 6 / 12`.

---

## 6. Leg 3 — mechanical coherence, and the JK split (F7.5, T7.3)

The block prior acts on coefficient blocks, so coherence is checked
there (`coherence_p12.csv`): for every (equation, free block, horizon
≥ 2) the distance of the block's posterior-mean coefficients from their
unshrunk OLS value on the same detrended regressors, under the adaptive
estimator relative to under FMAR,
`ratio = ||b_adaptive − b_ols|| / ||b_fmar − b_ols||`. Escaping blocks must
have `ratio < 1`, tightened blocks `ratio > 1`. Over 940 cells:

| τ class | n | median ratio, BLP-block | share < 1 | median ratio, BLP-pooled | share < 1 |
|---|---|---|---|---|---|
| τ > 2 | 344 / 400 | 0.82 | 1.00 | 0.84 | 1.00 |
| 1 < τ ≤ 2 | 211 / 179 | 0.97 | 0.98 | 0.97 | 0.98 |
| τ ≤ 1 | 385 / 361 | 1.00 | 0.26 | 1.01 | 0.16 |

Spearman correlation between `log tau` and `log ratio`: −0.89 (block),
−0.87 (pooled). The escapes in Section 3 are therefore exactly where the
adaptive estimators release the prior, and the tightened cells are where
they pull harder than the global prior. An IRF-level version of this check
(block IRF closer to the LP than the FMAR IRF where τ escapes) is not
informative here: an IRF mixes all blocks' contemporaneous coefficients
through the impact vector, and with the surprise block free it failed for
that reason alone (31% of cells, no relation to τ).

**JK information-effect split** (`tau_jk_comparison_p12.csv`; independent
estimator, early-horizon `tau_bar`, baseline → policy-only surprises):
most level cells move down modestly (`ip ← i1y` 5.33 → 4.42,
`hicp ← hicp` 9.45 → 8.37, `hicp ← i1y` 2.41 → 1.98, `i1y ← stoxx`
1.72 → 1.18), while the `stoxx` equation moves up (`stoxx ← i1y`
2.58 → 2.73, `stoxx ← ip` 1.93 → 2.22). The pre-registered expectation of
*weaker* stoxx-related escapes under policy-only surprises is not
supported. The JK variant's first stage is no stronger (impact 0.28,
t = 0.78).

---

## 7. Leg 2 — the null calibration and the reading protocol (T7.2)

`null_calibration_p12.mat`, `null_thresholds_p12.csv`: `R = 500` datasets
simulated from the BVAR(12) fitted to the real data (residual resampling,
actual initial conditions), each pushed through the entire pipeline with
the baseline settings (floor on, surprise block held, `h1_mode = 'lp'`,
chains 300 + 700). The tightness objective was unimodal in every
replication; the median null `lambda_h` path runs 0.010 → 0.071 at
`h = 12`, against 0.106 on the real data. Band coverage of the
true-under-null IRF by the sandwich bands is 0.93–0.99 per equation at a
nominal 0.90.

**The null is wide.** Per-cell and per-equation 95th percentiles of the
early-horizon `tau_bar` under the null, with the real-data values
(`tau_protocol_p12.csv`):

**Independent estimator (BLP-block)** — entries are `tau_bar` (null q95); **bold** = above its own null q95.

| equation \ block | i1y | ip | hicp | stoxx | max over blocks (eq. q95) |
|---|---|---|---|---|---|
| mps | 0.52 (0.84) | 0.63 (0.89) | 0.60 (0.95) | 0.69 (0.95) | 0.69 (1.00) |
| i1y | 2.02 (5.75) | 1.34 (4.76) | 0.86 (3.07) | 1.72 (4.21) | 2.02 (6.30) |
| ip | **5.33 (4.98)** | 3.21 (4.34) | 1.95 (3.40) | 1.30 (3.62) | 5.33 (5.48) |
| hicp | 2.41 (6.08) | 1.82 (4.49) | **9.45 (6.03)** | 2.29 (4.29) | 9.45 (7.03) **escape** |
| stoxx | 2.58 (5.56) | 1.93 (5.93) | 0.63 (3.24) | 0.70 (4.69) | 2.58 (7.32) |

**Horizon-pooled estimator (BLP-pooled)** — entries are `tau_bar` (null q95); **bold** = above its own null q95.

| equation \ block | i1y | ip | hicp | stoxx | max over blocks (eq. q95) |
|---|---|---|---|---|---|
| mps | 0.07 (0.25) | 0.10 (0.29) | 0.05 (0.42) | 0.09 (0.32) | 0.10 (1.00) |
| i1y | 4.57 (7.27) | 2.41 (5.70) | 0.95 (3.64) | 2.70 (5.32) | 4.57 (7.70) |
| ip | **7.57 (6.27)** | 4.39 (5.69) | 2.07 (3.88) | 1.81 (4.49) | 7.57 (7.35) **escape** |
| hicp | 3.27 (7.59) | 2.34 (6.11) | **10.05 (7.64)** | 2.97 (5.62) | 10.05 (8.32) **escape** |
| stoxx | 4.39 (6.85) | 3.35 (7.18) | 0.38 (3.94) | 0.25 (6.09) | 4.39 (8.53) |

**Reading, applied mechanically.** Two cells exceed their own null q95
under both estimators or under the sharper one:

* the **HICP own-lag block** (9.45 against 7.03 independent; 10.05
  against 8.32 pooled), also the only cell that trips the equation-wide
  "anything escapes" flag under the independent estimator;
* the **interest-rate block in the industrial-production equation**
  (7.57 against 7.36 pooled; 5.33 independent, above its own cell's q95
  but just under the equation-wide 5.48).

Every other cell — including all the rate-block cells in the `i1y` and
`stoxx` equations that look dramatic in Section 3 — lies inside its null.
The surprise equation is quiet in the data and in the null alike.

Two features of the null deserve comment. First, the argmax statistic
does not attribute: under the null the free blocks of the `i1y`, `ip` and
`stoxx` equations win their equation's argmax with frequencies close to
the uniform 0.25, so the real-data argmax (`i1y` block in `ip` and
`stoxx`) is not evidence by itself; the level of the cell against its
q95 is.

| estimator | equation | i1y | ip | hicp | stoxx |
|---|---|---|---|---|---|
| block | i1y | i1y 0.30 | ip 0.19 | hicp 0.23 | stoxx 0.23 |
| block | ip | i1y 0.22 | ip 0.20 | hicp 0.31 | stoxx 0.23 |
| block | hicp | i1y 0.23 | ip 0.06 | hicp 0.64 | stoxx 0.07 |
| block | stoxx | i1y 0.20 | ip 0.28 | hicp 0.19 | stoxx 0.30 |
| pooled | i1y | i1y 0.28 | ip 0.15 | hicp 0.21 | stoxx 0.20 |
| pooled | ip | i1y 0.20 | ip 0.18 | hicp 0.32 | stoxx 0.19 |
| pooled | hicp | i1y 0.25 | ip 0.08 | hicp 0.59 | stoxx 0.06 |
| pooled | stoxx | i1y 0.18 | ip 0.26 | hicp 0.20 | stoxx 0.29 |

Second, in the `hicp` equation the HICP own block is the argmax in 59–64%
of null replications: a VAR(12) fitted to a seasonal NSA price index
estimates its own seasonal centre imprecisely even when it is the true
model, and τ reports that disagreement — the empirical counterpart of the
simulation study's `p = 3`/`p = 4` cells. The real-data value exceeds even
that inflated null, which is what makes it a finding rather than an
artefact.

Taken together with Section 2, the calibrated result is modest and
specific: at `p = 12` the FMAR prior is adequate for most of this system,
and the two places where the data disagree with it beyond what estimation
noise in the prior centre explains are the seasonal dynamics of prices and
the transmission of the policy rate to industrial production.

---

## 8. Companion exhibit — the four-variable level system

Same machinery on `i1y, ip, hicp, stoxx` alone (random-walk prior means,
recursive innovation to `i1y`, every block free; `tau_heatmap_levels_p12.csv`).
Early-horizon `tau_bar`, independent / pooled:

| equation \ block | i1y | ip | hicp | stoxx |
|---|---|---|---|---|
| i1y | 0.82 / 0.90 | 0.61 / 0.56 | 0.33 / 0.28 | 0.72 / 0.66 |
| ip | **1.80 / 1.88** | 1.16 / 1.20 | 0.68 / 0.63 | 0.58 / 0.44 |
| hicp | 0.85 / 0.80 | 0.67 / 0.63 | **3.31 / 3.31** | 0.87 / 0.87 |
| stoxx | 1.15 / 1.06 | 0.70 / 0.75 | 0.30 / 0.13 | 0.33 / 0.10 |

With a tightness path the level system chooses for itself (0.36 at
`h = 12`), the same two features stand out — the HICP own block and the
interest-rate block in the `ip` equation — at a third of the magnitude,
and everything else is at or below 1. Against its own null
(`null_calibration_p12_levels.mat`, `R = 200`; `tau_protocol_levels_p12.csv`)
the null is correspondingly narrower (equation-wide q95 of the largest
`tau_bar` 1.9–2.4), and **only the HICP own block survives** (3.31
against 2.16 independent, 2.30 pooled); the rate block in the `ip`
equation (1.80 against 2.11) does not. The comparison with Section 3 is a
direct measure of how much the instrument's presence inflates the
five-variable map, and the two calibrated readings agree on the one
robust feature.

A caution the level-system null exposes: the sandwich bands' coverage of
the true-under-null IRF is **0.48** at a nominal 0.90 for the rate's own
response (`i1y` equation; 0.82–0.92 for the other three), against
0.93–0.99 in the five-variable system, where the shock is the surprise
and the responses are small. A near-unit-root own response estimated by
local projections in levels is biased toward zero at intermediate
horizons, and a band centred on the shrunk estimate does not cover the
bias. Any IRF exhibit built on the level system needs the lag-augmented
LP or a bias correction alongside.

---

## 9. Against the pre-registered expectations (CH7_DESIGN Sec. 7)

| # | expectation | outcome |
|---|---|---|
| 1 | `p = 12`: few or no escapes; if any, in the stoxx equation at short h | Calibrated verdict (Section 7): **few escapes, yes** — two cells at q95, the HICP own block and the interest-rate block in the `ip` equation; **stoxx, no** — nothing in the stoxx equation survives its null. |
| 2 | `p = 2`: widespread escapes, HICP own block among the strongest | Not supported: read against its own null, `p = 2` shows one escape (the HICP own block), the same cell that survives at every `p`; nothing widespread (Section 5). |
| 3 | dose–response: `tau_bar` weakly decreasing in `p` | Not supported: in levels the path is confounded by the tightness; through the protocol the count is 1 / 1 / 2 at `p = 2 / 6 / 12` (independent), flat rather than decreasing. |
| 4 | negative control (`mps` equation) quiet at all `p` | **Supported** at every `p` and for both estimators. |
| 5 | JK variant: weaker stoxx-related escapes | Not supported: the stoxx-equation cells rise slightly under policy-only surprises while most other cells fall modestly (Section 6). |

Two things the design did not anticipate and the thesis must state: the
instrument is weak for a monthly-average indicator (Section 1), and the
published prior scale plus the block structure both need an adjustment
when a white-noise instrument sits inside the system (`README_EMPIRICAL.md`,
Status 2–3). Both are findings about applying VAR-centred local
projections to internal-instrument designs, not about the euro area.

---

## 10. What remains

* A stronger policy indicator: an end-of-month 1-year OIS level instead
  of the monthly-average Euribor would remove the two-month smearing and
  the bank-credit premium of 2008–2012.
* The instrument block: AGKL's zero restrictions (lagged surprises out of
  the horizon regressions) rather than holding a loose block at `tau = 1`;
  this also removes the remaining drag on `lambda_h` (Section 2).
* Robustness listed in the design and not yet run: `mps_gc_3m`,
  `mps_all_1y` (speeches), and a `K = 7` system with loans and spreads.
