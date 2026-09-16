# Chapter 7 — results of the euro-area application (current reading)

Headline system `ois4`, external-instrument design (`DESIGN.md`), read
with the protocol revised on 2026-09-15. Every number is read from a file
in `results/empirical/`, named next to each table. The chronological record
of all runs, including the v1 design and the 2026-09-14 first reading of
this system, is `RESULTS_LOG.md`; where the two disagree this file wins.

**Design.** `T = 240` (2000m1–2019m12; effective sample 2001m1–2019m12),
`K = 4` (`ois1m` end-of-month 1-month OIS, `ip`, `hicp_sa`, `stoxx`),
`p = 12`, `H = 48`, FMAR mode with `h1_mode = 'lp'`, prior-scale floor on,
every block free, real-data chains 500 + 1500, early-horizon window
`h = 2..12`. Headline instrument by the pre-specified rule:
`z_gcs_1m_adj_sum` (1-month OIS surprise, day-count adjusted, meetings and
speeches, within-month sum). Nulls: iid residual resampling (`R = 1000` at
`p = 12`, `R = 200` at `p = 2, 4, 6`) and the wild bootstrap that keeps each
month's own residual variance (`R = 500` at `p = 12`, `300` at `p = 2, 4, 6`).

---

## 1. Instrument relevance and where it comes from

First stage of the end-of-month 1M-OIS innovation on the instrument
(`iv_ois4_p12_relevance.csv`; effective F = HAC t², Montiel Olea–Pflueger
thresholds 12.05 / 15.06 / 23.11 / 37.42 for 30 / 20 / 10 / 5 % bias):

| instrument | F_eff | t | R² | LP-IV F | lead+1 t | F 2001–08 | F 2009–11 | F 2012–19 | F excl. 2008m9–09m6 |
|---|---|---|---|---|---|---|---|---|---|
| **`z_gcs_1m_adj_sum`** (headline) | **13.7** | 3.71 | 0.103 | 19.4 | −2.38 | 2.0 | 44.2 | 5.3 | 8.7 |
| `z_gc_1m_adj_sum` (meetings only) | 7.9 | 2.80 | 0.047 | 11.2 | −2.50 | 0.4 | 26.3 | 34.3 | 10.5 |
| `z_gcs_1m_adj_sum_bs` (orthogonalised on pre-event information) | 15.0 | 3.87 | 0.115 | 19.2 | −1.83 | 2.7 | 46.2 | 3.4 | 9.3 |
| `z_gcs_3m_sum` (3-month surprise) | 19.2 | 4.38 | 0.141 | — | −2.35 | 3.6 | 39.8 | 5.2 | 33.1 |
| `z_gcs_1m_adj_kilian` (day-weighted; wrong convention for an eom level) | 7.9 | 2.81 | 0.044 | 7.5 | −0.49 | 0.4 | 19.7 | 5.3 | 5.2 |
| `z_gc_1m_adj_jk_sum` (JK poor man's, policy events) | 2.7 | 1.65 | 0.021 | 7.2 | −2.58 | 0.0 | 33.7 | 12.6 | 5.0 |
| `z_ext_jk_mp_pm` (Jarociński MP shock) | 2.3 | 1.53 | 0.033 | 10.8 | −2.58 | 0.0 | 32.2 | 5.8 | 17.4 |
| legacy v1 (1-year surprise on the 12M Euribor) | 0.5 | 0.74 | 0.010 | — | — | — | — | — | — |

**Anatomy of the headline first stage** (`ea_relevance_iv`, block
`anatomy`; the signed share of Σ (z_t − z̄)(u_t − ū) by sub-period):

| | 2001–2008 | 2009–2011 | 2012–2019 | top 5 months | months carrying half the absolute mass |
|---|---|---|---|---|---|
| share of the first-stage covariance | 29 % | **63 %** | 7 % | 53 % | 10 of 228 |
| sd of the OIS innovation (pp) | 0.125 | 0.158 | 0.044 | | |

The five largest months are 2001m5 (15 %), 2001m4 (10 %), 2008m10 (10 %),
2011m2 (10 %) and 2011m11 (9 %). The identification is carried by the
2009–2011 sovereign-debt years and a few large surprises of 2001 and
2008; the 2012–2019 short end contributes nothing. The naive regression of
the monthly rate change on the surprise (AGKL's Figure 9) has F = 0.7
here, against 12.6 in AGKL's sample that includes the 2022–23 hiking
cycle: the innovation regression is what carries the first stage.

**Timing and exogeneity.** No lag placebo is significant; the lead
placebo `u_{s,t}` on `z_{t+1}` has t = −2.4 (a rate innovation this month
predicts a negative surprise next month), which the orthogonalised
instrument brings to −1.8. The instrument is unpredictable from the four
pre-event regressors (p = 0.39) but partly predictable from the VAR's own
48 lags on event months (R² 0.27, HAC Wald p ≈ 0.001, q = 48; the
one-regressor test on the BVAR-predicted rate change is insignificant,
t = −1.1). The proxy identification conditions on those lags, so this is
admissible, but it means the surprises are not news relative to the VAR's
information set.

## 2. IRFs

`iv_ois4_p12_irf.csv`, `iv_ois4_p12_lpiv.csv`; 25 bp normalisation on the
end-of-month 1M OIS; BLP-block point with 90 % sandwich band; LP-IV
Anderson–Rubin set (weak-instrument robust).

| variable | h | BLP-block | LP-IV AR set |
|---|---|---|---|
| ois1m | 6 | 0.37 [0.23, 0.50] | [−0.18, 0.32] |
| ois1m | 12 | 0.34 [0.16, 0.51] | [−0.54, 0.28] |
| ip | 3 | 1.63 [0.79, 2.46] | [−0.86, 1.63] |
| ip | 12 | 0.82 [0.00, 1.64] | [−3.13, 0.73] |
| ip | 24 | −0.98 [−1.81, −0.15] | [−0.83, 1.96] |
| hicp_sa | 12 | 0.15 [−0.03, 0.33] | [−1.23, 0.15] |
| stoxx | 3 | 4.33 [2.00, 6.65] | [2.23, 13.25] |
| stoxx | 12 | 1.95 [−1.69, 5.59] | [2.90, 16.58] |

A positive short-rate surprise is followed by higher stock prices and, in
the Bayesian LPs, higher output for a year with flat prices: the euro-area
information-effect pattern (Jarociński–Karadi, Kerssenfischer, Ricco–
Savini–Tuteja). The "pure policy" instruments that would remove it (JK
split, Jarociński MP) are too weak here (F 2–3) to overturn it, though
their impact vectors do flip the stock sign. Two cautions on top of the
first stage: with F_eff between the 30 % and 20 % bias thresholds the
Anderson–Rubin sets, not the sandwich bands, are the intervals to quote;
and under the null the sandwich bands cover the true rate response in
only **41 %** of replications at a nominal 90 % (`null_thresholds_iv_ois4_p12.csv`,
`eq_coverage_under_null`, iid null; 59 % for `ip`, 75 % for `hicp_sa`,
77 % for `stoxx`; 54 / 64 / 78 / 83 % under the wild null), the levels-LP
bias toward zero at intermediate horizons for a
near-unit-root own response. The rate's own IRF should be shown with the
lag-augmented LP or not at all.

## 3. The tau map against its nulls

Early-horizon mean scale `tau_bar(i, g)` over `h = 2..12`, read through
`ea_apply_protocol`: ratio to the cell's own null 95th percentile, Monte
Carlo p-value `p_mc = (r + 1)/(R + 1)`, and the family-wise p-value
`p_fwer` from the max statistic over the 16 free cells. Both estimators;
the pooled one re-run at `H = 12` for the protocol (its statistic moves by
up to 0.3 between `H = 48` and `H = 12`; the independent estimator's does
not move at all). Files: `iv_ois4_p12_tau_protocol.csv` (iid, `R = 1000`),
`iv_ois4_p12_tau_protocol_wild.csv` (wild, `R = 500`); the 2026-09-14
reading against the `R = 200` iid null is archived as `*_R200`.

| cell | estimator | iid null (R = 1000): ratio to q95 | p_mc | p_fwer | wild null (R = 500): ratio to q95 | p_mc | p_fwer |
|---|---|---|---|---|---|---|---|
| **ip ← ois1m** | block | 1.90 | 0.001 | **0.002** | 2.01 | 0.002 | **0.006** |
| **ip ← ois1m** | pooled | 2.00 | 0.002 | **0.006** | 2.12 | 0.002 | **0.006** |
| stoxx ← ois1m | block | 1.24 | 0.014 | 0.138 | 1.16 | 0.028 | 0.224 |
| stoxx ← ois1m | pooled | 1.49 | 0.004 | **0.039** | 1.32 | 0.018 | 0.098 |
| ip ← ip | block | 1.25 | 0.017 | 0.131 | 1.31 | 0.020 | 0.092 |
| ip ← ip | pooled | 1.23 | 0.016 | 0.167 | 1.29 | 0.020 | 0.122 |
| stoxx ← ip | block | 1.05 | 0.035 | 0.410 | 1.05 | 0.036 | 0.419 |
| stoxx ← ip | pooled | 1.20 | 0.016 | 0.201 | 1.24 | 0.008 | 0.170 |
| ois1m ← ois1m | block / pooled | 0.94 / 1.01 | 0.064 / 0.049 | 0.680 / 0.486 | below 0.9 | | |
| every other free cell | both | below 0.9 | | | below 0.9 | | |

Bold: family-wise p ≤ 0.05. With `R = 1000` Holm's step-down on `p_mc`
(attainable once `R ≥ 319`) rejects exactly the `ip ← ois1m` cell under
both nulls and both estimators, and nothing else.

**Reading.** One cell survives the family-wise 5 % test under both
nulls and both estimators: the short-rate block of the output equation
(`ip ← ois1m`), at twice its null 95th percentile. The stock-equation
short-rate cell clears 5 % only for the pooled estimator under the iid
null (0.039; 0.138 for the independent one) and for neither under the
wild null (0.098 / 0.224), whose 95th percentile for that cell is a fifth
higher once the short rate's own heteroskedasticity is respected; the
output equation's own-lag block and the rate's own block are inside both. Two of four
equations trip the equation-wide flag under either null. What the wild
null changes is where the short-rate innovation's variance matters: the
95th percentile of the short-rate block's scale rises by 24 % in the
rate's own equation and 21 % in the stock equation, and by nothing in
the output equation.

The `hicp_sa` equation is quiet under every null (ratios 0.5–0.7): with
seasonally adjusted prices the HICP own-lag escape of the NSA level system
(`RESULTS_LOG.md`, V2.3–V2.5) is gone, as the year-on-year variant
predicted. That escape, the only one in the project validated out of
sample (+56 % two-month MSFE when ablated), is residual seasonality that a
random-walk-centred VAR(12) prior under-fits: a mechanical
misspecification that the diagnostic finds correctly and that disappears
when the data are fixed.

## 4. Across lag orders

`iv_ois4_cross_p.csv`: escape counts per design and null, `p ∈ {2, 4, 6, 12}`.
Raw `tau_bar` is never compared across `p` (it is measured relative to a
tightness that changes with `p`); only null-relative quantities are.

| design | null (R) | estimator | cells above q95 | Holm | family-wise ≤ 5 % | min p_fwer | largest cell (ratio to q95) |
|---|---|---|---|---|---|---|---|
| p = 2 | iid (200) | block / pooled | 1 / 1 | 0 / 0 | 1 / 1 | 0.020 / 0.025 | ip ← ois1m (1.87 / 2.43) |
| p = 2 | wild (300) | block / pooled | 1 / 1 | 0 / 0 | 0 / 0 | 0.113 / 0.146 | ip ← ois1m (1.40 / 1.53) |
| p = 4 | iid (200) | block / pooled | 1 / 1 | 0 / 0 | 0 / 0 | 0.378 / 0.159 | stoxx ← hicp_sa (1.10 / 1.44) |
| p = 4 | wild (300) | block / pooled | 2 / 1 | 0 / 0 | 0 / 0 | 0.449 / 0.100 | stoxx ← hicp_sa (1.04 / 1.42) |
| p = 6 | iid (200) | block / pooled | 1 / 1 | 0 / 0 | 0 / 0 | 0.373 / 0.537 | ip ← hicp_sa (1.09 / 1.02) |
| p = 6 | wild (300) | block / pooled | 0 / 0 | 0 / 0 | 0 / 0 | 0.558 / 0.654 | ip ← hicp_sa (0.99 / 0.93) |
| p = 12 | iid (1000) | block / pooled | 4 / 5 | 1 / 1 | 1 / 2 | 0.002 / 0.006 | ip ← ois1m (1.88 / 2.00) |
| p = 12 | wild (500) | block / pooled | 4 / 4 | 1 / 1 | 1 / 1 | 0.006 / 0.006 | ip ← ois1m (1.98 / 2.12) |

The `p = 12` rows come from step D's own chains (seeded independently of
step A), which is why their ratios differ from Section 3 in the second
decimal; the reading is the same. At `p = 2, 4, 6` the nulls have `R = 200`
(iid) or `300` (wild), so Holm cannot reject there by construction and the
family-wise column is the one to read.

The output-equation short-rate cell is the only cell that ever clears
the family-wise test. It does so at `p = 12` under both nulls, at `p = 2`
under the homoskedastic null only (ratio 1.87 → 1.40 and family-wise
p 0.020 → 0.113 once the null keeps each month's own volatility), and
nowhere at `p = 4` or `6`, where the largest cells are price-block cells
at ratios 1.0–1.4 with family-wise p above 0.10 under every null. This
is not the monotone dose–response the v1 design pre-registered, and the
sharpness of the `p = 12` reading owes something to the fact that the
first stage is also strongest at `p = 12` only through the instrument,
which the tau map does not use: the two are independent by construction.

## 5. Out-of-sample block ablation

`ablation_iv_ois4_p12.csv` (60 expanding-window origins from 2010m1,
forecasts at `h = 2..12`, the whole pipeline re-estimated at each origin;
the escaping cell forced back to `tau = 1`; positive Clark–West = the
escape carries forecast content) and `ablation_iv_ois4_p2_from2005.csv`
(168 origins from 2005m1 at `p = 2`).

| cell | role | p, origins | h | MSFE ratio ablated/free | DM (HLN) | Clark–West |
|---|---|---|---|---|---|---|
| ip ← ois1m | escape | 12, 2010– | 2 | 0.99 | −0.45 | −0.38 |
| ip ← ois1m | escape | 12, 2010– | 12 | 0.84 | −1.42 | −0.69 |
| stoxx ← ois1m | escape | 12, 2010– | 2 | 1.00 | −0.36 | −0.25 |
| stoxx ← ois1m | escape | 12, 2010– | 12 | 0.92 | −0.75 | 0.01 |
| ois1m ← ois1m | escape | 12, 2010– | 12 | 0.92 | −2.30 | −2.68 |
| ois1m ← hicp_sa | control | 12, 2010– | 12 | 0.92 | −1.94 | −2.30 |
| ip ← ois1m | escape | 2, 2005– | 2 | 1.03 | 1.44 | 1.83 |
| ip ← ois1m | escape | 2, 2005– | 12 | 0.98 | −1.06 | −1.14 |

None of the `ois4` escapes buys forecast accuracy: releasing the
short-rate block is neutral at two months and slightly worse at twelve,
the control behaves the same way, and the one marginal gain (Clark–West
1.83 at `h = 2` over 2005–2019 at `p = 2`) is at the one-sided 5 % level
only. Two readings remain compatible with this: the disagreement lives in
the pre-2010 part of the sample that the evaluation window never has to
forecast, or the released coefficients are too noisy for the flexibility
to pay.

## 6. What the chapter can say

1. **On the diagnostic.** With seasonally adjusted data the FMAR prior is
   adequate for most of this four-variable system. The one place the data
   disagree with it beyond what a correctly specified VAR(12) produces is
   the short-rate block of the output equation at `p = 12`: at twice its
   null 95th percentile, family-wise p 0.006 under both the iid and the
   wild-bootstrap null, for both estimators. It is absent at `p = 4` and
   `6`, appears at `p = 2` only under the homoskedastic null, and carries
   no out-of-sample forecast content. The chapter should present it as a
   calibrated prior–data disagreement about the delayed pass-through of
   the short rate to output that a twelve-lag VAR prior misses, robust to
   the short rate's volatility regime but without predictive payoff, not
   as a transmission finding.
2. **On the positive control.** The NSA-price seasonal escape (level
   system) is the clean demonstration: found at every lag order, validated
   out of sample, gone when the data are seasonally adjusted.
3. **On identification.** The IRFs are those of the 2009–2011 surprises.
   The information-effect pattern is a conditional result about this
   instrument; every IRF statement carries the first stage (F 13.7, LP-IV
   19.4), the Anderson–Rubin sets and the anatomy of Section 1.
4. **On inference.** Holm's step-down over 16 cells needs at least 319
   null draws to be able to reject; the family-wise max-statistic
   p-value is the test to read, and the heteroskedasticity-robust null is
   the one to pre-specify next time.

## 7. Against the pre-registered expectations (`DESIGN_V1_LEGACY.md` Sec. 7, `DESIGN.md` Sec. 7)

| # | expectation | outcome |
|---|---|---|
| 1 | `p = 12`: few or no escapes | Few: one cell at the family-wise 5 % level under both nulls (the stock-equation cell only under the iid null, at the boundary). |
| 2 | `p = 2`: widespread escapes | Not supported: one cell (the same one) under the iid null, none under the wild null. |
| 3 | dose–response in `p` | Not supported at any lag order under the wild null; under the iid null the cell appears at `p = 2` and `12` and not at `4`, `6`. |
| 4 | negative control | Not applicable in the instrument-free design; the `hicp_sa` equation is quiet everywhere. |
| 5 | matched instrument strengthens the first stage | Supported (0.5 → 13.7), but the strength is a 2009–2011 phenomenon. |

## 8. What remains

* A null that also draws the fitted VAR's parameters (posterior draws
  rather than the posterior mean) to carry the prior-centre uncertainty
  that `APPROXIMATIONS.md` prices at 16 % of the band width.
* The sample to 2024 with a COVID volatility treatment, so that the
  2022–23 hiking cycle gives the 1-month instrument variance where it now
  has none; an end-of-month stock index; an information-purged instrument
  with a first stage.
* The lag-augmented LP as the interval for the rate's own response.
