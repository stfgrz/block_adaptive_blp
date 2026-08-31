# Chapter 7 — Empirical application design
## "Where does the VAR prior fail in euro-area monetary transmission?"

Status: first full draft (2026-09-01). Companion code: `data/`, `dgp/`,
`estimators/`, `montecarlo/`, `RUN_EMPIRICAL.m`. Integration notes:
`README_EMPIRICAL.md`.

---

## 0. Headline claim and scope

The chapter is a **methods demonstration**: the block-adaptive BLP runs on
real euro-area data, produces IRFs that are sane relative to three
benchmarks (BVAR, BLP-FMAR, lag-augmented LP), and — the new part — the
tau map delivers a *calibrated, interpretable* answer to "where does the
VAR prior fail here?". Every economic statement is conditional and
demonstrative; the chapter never claims a new stylized fact about euro-area
transmission. The upgrade path to a "diagnostic-as-product" paper
(headline (c)) is in Sec. 10.

Three legs:

* **Leg 1 (the finding):** the tau heatmap on the real data — which
  variable-blocks escape the VAR prior, in which equations, at which
  horizons.
* **Leg 2 (the inferential benchmark):** a parametric bootstrap under the
  null "the fitted BVAR *is* the DGP" calibrates the null distribution of
  the tau statistics for this exact design (K, T, p, persistence, trends).
  Real-data escapes are read against these thresholds, not against eyeball
  or the retired max/median rule. Code: `montecarlo/run_null_calibration.m`.
* **Leg 3 (validation without ground truth):** a dose–response in the lag
  order p ∈ {2, 6, 12}, a mechanical-coherence check (the block IRF should
  move away from the BVAR toward the LP exactly where tau escapes), and the
  JK information-effect re-run as a design-sensitivity dimension.

## 1. The shock series (Q4)

**Data.** EA-EMPD (Altavilla, Gürkaynak, Kind & Laeven, "Monetary
Transmission with Frequent Policy Events", ECB WP 3157 / CEPR DP 20196),
tick-cleaned event-window asset-price changes for all GC meetings and all
Executive-Board/President speeches since 1999. We use the **monetary event
window** (`GC_ME`: press release + conference, the union window), the
standard object in the euro-area HF literature since Altavilla et al.
(2019).

**Aggregation.** Sum of event surprises within the calendar month, zero in
months without events — exactly AGKL's construction for their own
transmission BVAR (their Sec. 4.3). Monthly panel produced by
`data/build_shock_series.m`; window statistics for 2001m1–2019m12: 221
GC_ME events in 206 of 228 months, monthly std 4.29 bp.

**Maturity — decision: OIS 1Y baseline, OIS 3M robustness.** The genuinely
contested choice; both sides, then the argument.

*The case for a short maturity (1M/3M, the "target" tradition).* Kuttner
(2001) logic: the shortest contract spanning exactly one meeting measures
the *policy-rate* surprise. AGKL themselves use the day-count-adjusted 1M
OIS as their surprise measure, and JK (2020, AEJ:Macro) use the 3M EONIA
swap for the euro area. A 1Y surprise mixes target news with
forward-guidance news, so the structural shock is a composite.

*The case for 1Y (the "policy indicator" tradition).* Gertler–Karadi
(2015) and Miranda-Agrippino–Ricco (2021) deliberately use a ~1-year rate
as the policy indicator: it reflects the effective stance *including
guidance*, which is the relevant stance concept when much of policy is
communication. Two sample-specific facts decide it here:

1. **The ELB.** Our sample (2001–2019) contains eight years in which the
   short end was pinned. Measured in the file: GC_ME surprise std in
   2012–2019 is **1.76 bp for OIS 3M vs 2.30 bp for OIS 1Y** (and 4.94 bp
   for 1Y pre-2009). The 1M is worse still. A short-maturity surprise
   series over 2012–2019 is close to degenerate → weak instrument →
   exactly the situation Bauer–Swanson (2021) warn about. AGKL can afford
   the 1M because their sample runs to 2025 and includes the hiking cycle;
   ours, by design, does not.
2. **Data completeness.** In the EA-EMPD file, GC_ME rows have **0 NaN in
   OIS_1Y from 1999**, vs 6 NaN in OIS_3M and 25 in OIS_1M (early years).

Also pragmatic: the AGKL day-count adjustment is defined for the 1M
contract (factor 30/(30−m)); for GC meetings m = 0 so the factor is 1, but
for the speeches extension the adjustment is material and noisy — the 1Y
sidesteps the mechanics entirely (accrual effects are second-order at that
maturity).

*Honest caveats to state in the thesis:* (i) the 1Y surprise identifies a
composite target+FG shock — acceptable with a single instrument, and the
composite is arguably the policy-relevant stance shock (GK tradition);
(ii) comparability with AGKL is at the design level, not the maturity
level — hence the OIS 3M robustness run; (iii) the Bauer–Swanson
predictability critique applies to all HF surprises and is out of scope
(AGKL's own predictability tests, their Tables 3–4, find surprises largely
unpredictable from macro news, with partial predictability from pre-event
stock moves on meeting days).

**Variants produced** (all by `build_shock_series.m`):
`mps_gc_1y` (baseline), `mps_gc_3m` (robustness), `mps_gc_jk` /
`mps_gc_info` (JK poor-man split, Sec. 6), `mps_all_1y` (meetings +
filtered speeches, extension), `mps_all_tgt` (AGKL-style day-count-adjusted
target composite; our documented cap: drop events with adjustment factor
> 6).

## 2. The system (Q6)

K = 5, ordering = identification (surprise first; recursive/internal
instrument as in AGKL, Stock–Watson 2018, Plagborg-Møller–Wolf 2021):

| # | var | transform | why |
|---|-----|-----------|-----|
| 1 | `mps` | pp (bp/100) | the instrument; zero-restricted block per AGKL eq. (7)-(8) |
| 2 | `i1y` | percent, level | stance level: needed for the +25 bp impact normalisation (Q5) and gives the diagnostic a "policy block" |
| 3 | `ip`  | 100·log | AGKL baseline outcome; the real block |
| 4 | `hicp`| 100·log | AGKL baseline outcome; the nominal block |
| 5 | `stoxx`|100·log | fast financial block: a-priori most likely locus of VAR-prior failure at short horizons; also ties to the JK narrative |

Rationale. AGKL's baseline is {mps, IP, inflation}; we add the stance
level (GK-style indicator, and without it the 25 bp normalisation has no
anchor) and one fast financial variable. Euro-area monthly VARs
conventionally include IP + HICP (e.g. the JK euro-area application and
the literature following it); the stock index doubles as the info-effect
classifier input. K = 5 keeps m = 1 + K·p = 61 regressors at p = 12
against T_eff = 228 − h, which is tight but exactly the regime BLP
shrinkage is for; adding more variables (credit spread, loans — AGKL's
extended set) is a listed extension, not baseline, because every extra
variable costs 12 regressors per equation and dilutes the block diagnostic.

**Normalisation (Q5).** IRFs are estimated in unit-surprise-innovation
scale and rescaled at reporting time by k = 0.25 / θ_bvar,i1y(0), i.e. a
shock normalised to raise the 1Y rate by 25 bp on impact. One common k
(from the BVAR impact) across estimators so their differences are not
confounded with normalisation differences.

**Deterministics/persistence.** FMAR machinery as in Ch. 5–6: GLP BVAR
with RW prior means for persistent variables. `isrw` must be the vector
[0 1 1 1 1] (white-noise centre for the surprise, RW for levels) — the
one-line patch in `estimate_bvar_niw.m` (README). HICP is used NSA
(Eurostat index); p = 12 absorbs seasonality at baseline, but note the
p = 2 dose-response leg deliberately does *not* — that is part of the
design (Sec. 5.2): residual seasonality is a known, real form of dynamic
misspecification at p = 2.

**1Y level source.** 12M Euribor (monthly average, from 1999) as the
practical default; it embeds a bank credit premium, material 2008–2012 —
disclose, and swap in a 1Y OIS level series if located in the ECB portal
(`fetch_outcome_data.m` documents both).

## 3. Sample (Q8)

Y spans **2000m1–2019m12** (T = 240); with p = 12 the first year is lag
initialisation, so the effective sample is **2001m1–2019m12** — matching
AGKL's January-2001 start and the user's requested window, and ending
before COVID. Y is held fixed across the p ∈ {2,6,12} runs (only the
prior's flexibility changes, not the data). Extension beyond 2019: Sec. 10.

## 4. Estimators

1. **BVAR** (`estimate_bvar_niw`) — the prior taken literally.
2. **BLP-FMAR** (`estimate_blp_fmar`) — scalar per-horizon escape
   (Ferreira–Miranda-Agrippino–Ricco, REStat 2025).
3. **BLP-block** (`estimate_blp_blockadaptive`) — the thesis estimator;
   tau map is Leg 1.
4. **LP-LA** (`estimators/estimate_lp_lagaug.m`, new) — Montiel
   Olea–Plagborg-Møller (Ecta 2021) lag-augmented LP with EHW bands: the
   frequentist benchmark an LP audience trusts, valid under the near-unit
   roots we certainly have. No HAC tuning to argue about.

Plus the plain LP of the repo where useful for continuity with Ch. 6.

## 5. The validation design

### 5.1 Leg 2 — null calibration (the key inferential exhibit)

On real data there is no ground truth, so the tau map alone is
descriptive. `run_null_calibration.m` simulates R datasets from the
BVAR fitted to the real data (residual resampling; initial conditions =
actual first p months) and pushes each through the *entire* pipeline
(BVAR re-estimation, FMAR λ re-selection, Gibbs). Outputs per
(equation, block):

* null mean and q90/95/99 of tau_bar = mean tau over h ∈ [2,12];
* null argmax frequencies (the empirical analogue of the Ch. 6
  correct-DGP false-positive benchmark, 0.330);
* per-equation q95 of max_g tau_bar (calibrated "anything escapes" flag —
  replaces the retired max/median > 1.5 rule, which Ch. 6 showed is
  miscalibrated in fmar mode);
* band coverage under the null (a free sanity check on the sandwich bands
  in the empirical design).

**Reading protocol (pre-registered):** a block escape is *reported* iff
tau_bar(i,g) > q95(i,g); the *attribution* statement uses the argmax
compared to its null frequency; the *pooled* flag uses maxstat_q95. All
statements in the thesis text follow this protocol mechanically.

### 5.2 Leg 3 — dose–response in p (and what it means)

p is the lag order of both the LP controls and the VAR that generates the
prior centre. At p = 12 a monthly 5-variable VAR is flexible: the prior
may simply be adequate, tau stays near/below 1 everywhere, and the heatmap
is quiet. **A quiet map at p = 12 is a finding, not a failure** — the
diagnostic is telling us the FMAR prior is fine for this design, which is
precisely what a specification diagnostic should say when nothing is
wrong. To show the diagnostic *would* speak if something were wrong, we
turn the dial we control: estimate at p = 2, 6, 12 on identical data. At
p = 2 a monthly VAR is knowingly too short (seasonality, delayed
transmission) — a real, not simulated, misspecification. Predictions:

* escapes light up at p = 2 (positive control on real data),
* fade monotonically as p grows (dose–response),
* and at p = 2 concentrate on blocks/equations whose deeper lags matter
  (e.g. own deep lags of hicp given NSA seasonality; deep real-block lags
  in the ip equation).

Together with Leg 2 this gives a falsifiable validation story with no
access to ground truth. (Development is unaffected by this subtlety: the
runs are the same code with three values of cfg.p.)

### 5.3 Mechanical coherence and negative control

* Coherence: at cells where tau escapes, |θ_block − θ_LP| < |θ_FMAR −
  θ_LP| should hold (the block estimator releases the prior exactly
  there); at quiet cells θ_block ≈ θ_FMAR. Scatter exhibit from saved
  IRFs.
* Negative control: the surprise equation (variable 1). The surprise is
  built from event windows; its dynamics are close to white noise and the
  VAR nests white noise, so its equation's blocks should be quiet. An
  escape there would flag a construction problem, not economics.

## 6. Information effects as a design dimension (Q10, option (b))

Rerun the whole Leg-1 exhibit with `mps_gc_jk` (JK "poor man's sign
restriction": keep the GC surprise only when the STOXX50E window move has
the opposite sign; 114/221 events survive, monthly std 3.52 bp). Question
asked of the diagnostic: *does cleaning information effects change where
the VAR prior fails?* Plausible pattern: the stoxx-block escapes in the
baseline shrink under JK surprises (information events are exactly where
stock–rate comovement violates the recursive VAR's implied dynamics). Any
outcome is reportable; this is a sensitivity dimension, not a hypothesis
test. `RUN_EMPIRICAL.m` step D emits the side-by-side early-h tau table.

## 7. Pre-registered expectations

Stated before running, to keep the demonstration honest:

1. p = 12 baseline: few or no escapes at q95; if any, most likely in the
   stoxx equation / fast-financial cells at h ≤ 6.
2. p = 2: widespread escapes; hicp own-block among the strongest (NSA
   seasonality at lag 12 truncated away).
3. Dose–response: median tau_bar weakly decreasing in p, cell-wise.
4. Negative control (mps equation): quiet at all p.
5. JK variant: weakly fewer/weaker escapes in stoxx-related cells.

Deviations from these are reported as such, not massaged.

## 8. Exhibits

* T7.1 Data and shock construction summary (+ first-stage relevance, with
  the AGKL F≈4 meetings-only benchmark quoted).
* F7.1 Monthly surprise series, baseline vs JK vs speeches-augmented.
* F7.2 IRFs at p = 12, four estimators, 25 bp normalisation, 90% bands.
* F7.3 Tau heatmaps per equation (K × h), p = 12, with q95 contours from
  Leg 2. **The headline exhibit.**
* T7.2 Null-calibration thresholds (from `null_thresholds.csv`).
* F7.4 Dose–response: early-h tau_bar by cell across p ∈ {2,6,12}.
* F7.5 Coherence scatter (Sec. 5.3).
* T7.3 JK comparison table (`tau_jk_comparison.csv`).

## 9. Run plan and budgets

1. `build_shock_series` + `fetch_outcome_data` + `assemble_dataset` —
   minutes; verify first-stage prints.
2. `RUN_EMPIRICAL` steps A–D — each single estimation is one BVAR + K·H
   λ-searches + K·H Gibbs chains; expect tens of minutes per estimation at
   p = 12, H = 48 with full chains. Steps A+C+D ≈ 5 estimations.
3. `run_null_calibration` — the long job. Defaults R = 200, H = 12,
   reduced chains; **do `n_rep = 10` first** and scale from the printed
   s/rep. Checkpointed every 10 reps, resume-safe. Thesis run R = 500 if
   the quick timing allows (user's machine did 3×500 Ch. 6 reps in ~1 h,
   but these reps are heavier: m = 61 vs 7 regressors, K = 5 vs 3 —
   budget overnight).

## 10. Extensions (kept out of the baseline deliberately)

* **Speeches dual-headline (c-option from the planning discussion):**
  rerun everything with `mps_all_1y` (all 228 months have events; corr
  0.856 with the baseline series; AGKL's headline result is precision
  gains — first-stage F from ≈4 to ≈13). The block diagnostic question
  becomes: does a stronger instrument change where the prior fails?
* **Post-2019 / COVID:** Lenza–Primiceri (2022)-style volatility
  rescaling: give COVID months (2020m3–2020m12, estimated decay) inflated
  residual variances; implementation = row weights 1/s_t applied to (Y, X)
  in both the BVAR and every horizon regression, plus the same weights
  inside the Newey–West prior-scale and sandwich steps. Alternative:
  Cascaldi-Garcia "pandemic priors" (dummy per pandemic month). Either
  extends the sample to 2024 without letting three months dominate every
  likelihood. This is a clean, self-contained Ch. 8 / future-paper module.
* **Headline upgrade path to (c):** once Legs 1–3 behave, the paper
  version inverts the framing: the tau map (with calibrated thresholds)
  becomes the product — "a specification diagnostic for VAR-prior local
  projections, with an application to euro-area transmission" — and the
  euro-area map becomes the evidence. Requires: the null calibration at
  R = 500+, the dose–response, at least one economically legible escape
  story told carefully, and the speeches variant as external validity.
* AGKL extended system (loans + lending spreads) as a K = 7 robustness.

## 11. Decisions log (traceability)

| Q | decision |
|---|----------|
| 1 | headline (a) methods demo; (c) as upgrade path (Sec. 10) |
| 2 | Legs 2+3 in scope |
| 3 | meetings-only baseline; speeches = extension (Sec. 10) |
| 4 | OIS 1Y baseline (ELB variance + completeness + GK indicator logic); 3M robustness; both-sides argument Sec. 1 |
| 5 | +25 bp impact normalisation on the 1Y rate, common k |
| 6 | K = 5: mps, i1y, ip, hicp, stoxx (Sec. 2) |
| 7 | dose–response explained Sec. 5.2; no impact on development |
| 8 | 2001–2019 effective sample; LP-COVID extension sketched Sec. 10 |
| 9 | merged protocol: calibrated thresholds + early-h argmax (Sec. 5.1) |
| 10 | JK info split as design dimension (Sec. 6) |
| 11 | MOP-M lag-augmented LP added (Sec. 4) |
| 12 | compute budgeted in Sec. 9 |
| 13 | data-construction + null-calibration modules delivered together |
