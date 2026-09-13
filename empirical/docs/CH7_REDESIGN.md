# Chapter 7 redesign (v2): external-instrument identification, 1-month OIS, and null-calibrated validation

Status: **design and code prepared on 2026-09-13; the 1-month OIS data arrive on Monday.**
Nothing in this document claims that the redesign succeeds. Everything that
can be run without the OIS level series has been run on the existing data and
is labelled as such in `CH7_RESULTS.md` (section "v2 preview"). The legacy
baseline of 2026-09-12 (`RUN_EMPIRICAL.m`, tag `p12`) is untouched and remains
the reference point.

Literature support for every material choice, with what was verified against
the papers' text and what is from memory, is in `CH7_LITERATURE.md`. This
document cites it as [L-n].

---

## 0. What changes, in one paragraph

The surprise leaves the VAR state vector. The system is the K-variable block
`{policy indicator, ip, hicp, stoxx, (loans, lending spread)}`; the Bayesian
VAR, the FMAR tightness `lambda_h`, the block scales `tau` and the null
calibration are all computed on that system alone. The EA-EMPD surprise
enters only once, at identification: the impact vector is the covariance of
the VAR innovations with the monthly instrument, normalised to a unit effect
on the policy indicator (Stock–Watson 2018 [L-6]), and the structural IRF at
horizon h is the Bayesian-LP reduced-form coefficient block times that vector.
Under the maintained VAR(p) structure this equals the LP-IV estimand [L-6,
L-7]. The instrument is rebuilt from the full EA-EMPD workbook as the
day-count-adjusted 1-month OIS surprise (AGKL eq. 3 [L-1]), with and without
speeches, aggregated to the month in a way that matches the policy indicator's
timing convention (Kilian 2024 [L-5]), and read against a proper relevance
module. The tau map is then compared across `p ∈ {2, 4, 6, 12}` through each
design's own null.

## 1. Identification: external instrument (proxy), not an internal lag block

**Legacy (v1).** `mps_t` ordered first in the VAR with lagged surprises in
every equation. AGKL do exactly this (their eqs. 7–8, an internal instrument
in the sense of Plagborg-Møller–Wolf 2021, with `mps_t` first and no
contemporaneous feedback) [L-1]. It is valid under non-invertibility [L-7].
The cost, found on our data: twelve near-zero lag coefficients per equation
that (i) violate the block prior's common-scale assumption and (ii) drag the
global tightness down (README_EMPIRICAL Status 2–3). The v1 result stays as
the legacy baseline.

**v2.** External-instrument (proxy) identification of the impact vector:

    u_t = VAR(p) innovations of Y (posterior-mean BVAR, K x 1),  z_t = monthly instrument
    b_z = Cov(u_t, z_t) / Cov(u_{s,t}, z_t)          (unit effect on indicator s; SW2018 [L-6])
    theta(h) = Gamma_h b_z,   Gamma_h = horizon-h coefficient block on y_t of the (Bayesian) LP

`Gamma_h` is whatever the estimator produces (BVAR companion power, FMAR
posterior mean, block-adaptive posterior mean); only `b_z` changes relative
to the recursive `b1n`. Implementation: `ea_identify_proxy.m` returns a copy
of the BVAR struct with `b1n := b_z` (and `theta`, `theta_lo/hi`, `b1n_draws`
recomputed from the NIW draws), so every estimator in the repo is
proxy-identified without modification.

*What is assumed.* Relevance `Cov(u_s, z) ≠ 0`, contemporaneous exogeneity,
and lead–lag exogeneity given the p lags (SW2018 [L-6]); invertibility of the
policy shock in the K-variable system (the price of the external-instrument
route, PMW 2021 [L-7]). The lead–lag condition is testable and is tested
(Sec. 5). Invertibility is not testable here; the v1 internal-instrument run
is the robustness check for it, and the two are reported side by side.

*Sign and scale.* `b_z(s) = 1`: IRFs are per 1 pp innovation in the policy
indicator induced by the surprise; reported scale ×0.25 (25 bp). No more
`k25 = 0.25/impact` amplification: the normalisation is the estimand itself
and its precision is the first stage.

*Bands.* For the BLP estimators the sandwich bands hold `b_z` fixed, exactly
as FMAR hold `b1n` fixed (documented approximation, `docs/APPROXIMATIONS.md`).
The frequentist benchmark `estimate_lp_iv.m` carries the instrument
uncertainty in full: 2SLS per horizon with HAC, lag-augmented EHW as a
variant, and an Anderson–Rubin confidence set that stays valid when the
instrument is weak (Montiel Olea–Stock–Watson 2021 [L-8]).

*Why not AGKL zero restrictions on lagged surprises.* Verified against the
paper: AGKL do **not** restrict lagged surprises; their "zero restrictions"
are the contemporaneous Cholesky ordering [L-1]. Our earlier note
("AGKL-style zero restrictions on the instrument block") was wrong and is
withdrawn. Dropping the lag block of the surprise while keeping `z_t` as a
contemporaneous regressor is a third option (Ramey 2016 "shock as regressor"),
equivalent to the proxy route under the VAR structure; it is not implemented.

## 2. The instrument: 1-month OIS surprise, adjusted, with speeches

Built by `build_instrument_series.m` from `raw/ea_empd_events_full.csv` (all
54 columns of the workbook, extracted 2026-09-13; the shipped 12-column
extract is kept for the legacy builder). Every variant is a monthly series in
basis points, zero in months without events (AGKL [L-1]).

| dimension | options | default and reason |
|---|---|---|
| event set | `gc` (GC_ME), `gcs` (GC_ME + EB/P speeches), `sp` (speeches only) | `gcs` for the headline instrument: AGKL's first stage rises from F 4.0 (meetings) to 12.6 (all events) [L-1]; `gc` reported alongside |
| surprise | `1m_adj` = 30/(30−m)·ΔOIS_1M, m = days until next GC (AGKL eq. 3 [L-1]; m = 0 for meetings); `1m_raw`; `3m`; `1y` | `1m_adj` (matches the 1M indicator, Kuttner logic [L-3]) |
| speech filters | regular trading day and inside 9–18 CEST (AGKL baseline [L-1]); events with m ≥ 30 dropped (formula undefined); factor cap `max_adj` (default 6) | AGKL's exclusion of windows containing a Bloomberg relevance-90 data release **cannot** be replicated (no release calendar): documented deviation |
| aggregation | `sum` (within-month sum) or `kilian` (day-weighted with carry-over) | tied to the indicator, Sec. 3 |
| information effects | `jk` (poor man's sign restriction on the GC event, STOXX50E sign [L-2]); `bs` (orthogonalised on pre-event public information, Bauer–Swanson 2023 [L-9]); `ext` (externally supplied series, Sec. 6) | headline uses none; all three are robustness dimensions |
| stance dimension | `path1y` = event-level residual of ΔOIS_1Y on the adjusted 1M surprise (GC only) | a rough forward-guidance proxy in the spirit of the EA-MPD rotation, in which factors 2–3 are orthogonal to the 1M OIS [L-4]; **not** the EA-MPD FG factor (no 5–10Y in our rotation) |

Column names: `z_<set>_<surprise>[_<info>]_<agg>`, e.g. `z_gcs_1m_adj_sum`,
`z_gc_1m_adj_jk_kilian`, `z_gc_path1y_sum`.

**What a 1M-OIS shock is, and is not.** ABGMR 2019 find that the press-release
window contains a single factor, Target, and that "the change in the
one-month OIS directly as the Target factor ... made no difference" [L-4].
A 1M-OIS instrument therefore identifies the **conventional short-rate**
dimension of policy. It is not the forward-guidance or QE dimension, and over
2012–2019 (short end pinned) its variance is small: GC_ME sd 1.8 bp against
3.3 bp in 2001–08 (computed from the file). The thesis text must say
"conventional policy shock", never "the monetary policy stance". The `path1y`
proxy and the legacy 1Y instrument are the documented complements.

## 3. Monthly aggregation must match the indicator's timing convention

Kilian (2024) [L-5]: summing the daily surprises within the month gives the
cumulative level shift, which is the right monthly proxy only if the VAR
variable is measured **end-of-month to end-of-month**. If the VAR variable is
a **monthly average**, a surprise on business day d of a month with D
business days raises that month's average by (D−d+1)/D of the surprise and
the *next* month's average by the remaining (d−1)/D; the conventional sum
misdates part of the surprise, which weakens the first stage. Gertler–Karadi
(2015) recognise the same point [L-10]. Hence:

| policy indicator | aggregation | code |
|---|---|---|
| end-of-month 1M OIS level (`ois1m_eom`) | within-month sum | `agg = 'sum'` |
| monthly average of same-day closes: 1M OIS (`ois1m_avg`), EONIA (`eonia_avg`) | day-weighted with carry-over, event day counted | `agg = 'kilian'` |
| monthly average of a fixing set before the event: 12M Euribor (`i1y`, legacy), 1M Euribor | day-weighted with carry-over, event day not counted (11:00 CET fixing precedes the 13:45 release) | `agg = 'kilianfix'` |

Business-day (default) or calendar-day weights in `build_instrument_series`.
All three conventions are always produced, so the relevance table can show
which one the data prefer — but the *headline* pairing is fixed a priori as
in the table above, not chosen by the F statistic.

The **legacy finding re-read**: v1 paired a summed surprise with a
monthly-average indicator (12M Euribor). By Kilian's argument that pairing is
misaligned, which is one candidate explanation of the t = 0.77 first stage.
The `kilian` aggregation on the legacy indicator is computed now (Sec. 8) so
that this explanation is tested rather than asserted.

## 4. Policy indicator and system

**Indicator.** 1-month EONIA-swap (OIS) level, Refinitiv `EUREON1M=`, daily,
imported by `import_ois_daily.m` into monthly `eom` (last quote of the
month), `avg` (mean of daily quotes), `first` and `chg_eom`. Reasons: same
underlying asset as the instrument (AGKL's own first stage: "monthly changes
in the one-month OIS on the high-frequency shocks of the same underlying
asset", F 4.0/14.1/12.6 [L-1]); no bank credit premium (unlike Euribor,
material 2007–2012); available since 1999. Known events: EONIA re-based to
€STR + 8.5 bp on 2 Oct 2019 (a level convention, continuous), discontinued 3
Jan 2022 (outside the sample). Fallbacks on disk or one fetch away: EONIA
monthly average (`FM.M.U2.EUR.4F.MM.EONIA.HSTA`, verified) and 1M Euribor
average (`FM.M.U2.EUR.RT.MM.EURIBOR1MD_.HSTA`, verified; credit premium).

**Systems** (`assemble_dataset_v2`, variable registry in the file):

| name | variables | sample | isrw |
|---|---|---|---|
| `lev4` (v2 core, runs today) | `i1y, ip, hicp, stoxx` | 2000m1–2019m12 | 1 1 1 1 |
| `lev4_yoy` | `i1y, ip, hicp_yoy, stoxx` | 2000m1–2019m12 | 1 1 1 1 |
| `ois4` (Monday) | `ois1m, ip, hicp_sa, stoxx` | 2000m1–2019m12 | 1 1 1 1 |
| `ois6` (Monday) | `ois1m, ip, hicp_sa, loans, lend_spread, stoxx` | 2003m1–2019m12 | 1 1 1 1 0 1 |

HICP: the seasonally adjusted overall index is the ECB's own series
`ICP.M.U2.Y.000000.3.INX` ("working day and seasonally adjusted", verified
on the Data Portal; the `S` adjustment code does not exist). `hicp_yoy` =
100·(log P_t − log P_{t−12}) from the NSA index needs no adjustment and is
available today. Both are offered because the literature does both: log
levels (JK 2020, GLP-style priors) and inflation rates (AGKL's "inflation");
the RW prior centre for `hicp_yoy` is a choice (euro-area YoY inflation has a
monthly AR(1) near 0.98) and is switchable per variable.

Loans: adjusted loans to NFCs, index of notional stocks, SA
(`BSI.M.U2.Y.U.A20T.A.I.U2.2240.Z01.E`, from 2003m1, verified). Lending
spread: composite cost of borrowing for NFCs
(`MIR.M.U2.B.A2I.AM.R.A.2240.EUR.N`, from 2003m1, verified) minus the 1M OIS
(registry option: minus 12M Euribor). AGKL's extended system has "loan volumes
and lending rate spreads" [L-1]; the exact spread definition is not stated in
the text we could verify, so ours is a documented choice. A spread gets a
white-noise centre (stationary), the practice of Bańbura–Giannone–Reichlin
2010 for stationary series [L-11, from memory]. The `ois6` sample starts in
2003m1 because MIR and the adjusted-loans index do; with p = 12 that is 192
effective months, which is thin for K = 6 (m = 73 regressors), so `ois6` is a
robustness system, not the headline.

## 5. Relevance, timing and influence diagnostics (`ea_relevance_iv.m`)

Reported for every (instrument, indicator) pair, printed and saved:

* first stage `u_{s,t} = a + b z_t + e_t` on the VAR innovations: `b`,
  EHW and HAC t, **effective F** (for one instrument the HAC-robust F; MOP
  2013 critical values 12.05 / 15.06 / 23.11 / 37.42 for 30/20/10/5 % Nagar
  bias, verified [L-8]; Staiger–Stock rule of thumb 10), partial R², N,
  number of zero months, sd(z) and sd(u_s) in bp, correlation;
* the AGKL-style naive first stage `Δ i_t` on `z_t` (their Fig. 9) [L-1];
* **lead/lag placebo**: `u_{s,t}` on `z_{t+k}`, k = −3..3. Only k = 0
  should be significant. A significant k = −1 is the signature of summed
  surprises against a monthly-average indicator (Sec. 3, Kilian); a
  significant lead means the instrument is dated after the indicator moves
  or lead–lag exogeneity fails (SW2018 [L-6]);
* **predictability** (Bauer–Swanson 2023 [L-9], AGKL Tables 3–4 [L-1]):
  `z_t` on pre-event public information available on disk — STOXX return
  over the previous 3 months, 12-month IP growth, 12-month inflation, change
  in the indicator over the previous 3 months — F, R², and the
  orthogonalised residual `z_..._bs`. AGKL find only the pre-event stock move
  predictive on meeting days [L-1];
* **autocorrelation** of `z_t` (MAR 2021 project on own lags [L-12]);
* **influence**: leave-one-month-out `b` and `F`, DFBETA, the five most
  influential months, and the first stage with 2008m9–2009m6 excluded and
  with `z` winsorised at 1/99 %. AGKL note that the 50 bp cut of November
  2008 produced a sizeable positive surprise and the January 2009 cut a
  negative one [L-1]; RST find non-linear information effects concentrated
  in a few crisis events explaining up to ~40 % of short-maturity price
  revisions [L-13];
* **subsamples**: 2001–08, 2009–11, 2012–19.

## 6. Information effects

Kept as robustness, not baseline, because the literature disagrees on what the
correction removes (JK: central-bank information [L-2]; Bauer–Swanson: the
central bank's response to public news, "not an information effect" [L-9]).

1. `jk`: poor man's sign restriction on GC events (opposite sign of the OIS
   and STOXX50E window moves; zeros assigned to policy) [L-2].
2. `bs`: orthogonalisation on lagged public information (Sec. 5).
3. `ext`: an externally supplied monthly series. Ready-made:
   Jarociński's updated ECB shocks, `shocks_ecb_mpd_me_m.csv` in
   `github.com/marekjarocinski/jkshocks_update_ecb` (MP_pm, CBI_pm,
   MP_median, CBI_median; first PC of OIS 1M–1Y scaled to the OIS-1Y sd;
   CC BY 4.0) [L-2]; the RST informationally robust factors (projection on
   Eurosystem staff projections, linear and non-linear) are not published
   with the paper and would have to be requested [L-13]. Adapter:
   `import_external_instrument.m` (any csv with year, month and a value
   column, or the JK file by column name).

RST's headline — that linear corrections barely move euro-area IRFs while
non-linear crisis-period information effects resolve the output/price/stock
puzzles [L-13] — is why the influence diagnostics of Sec. 5 (crisis months
in or out) are part of the headline exhibit rather than a footnote.

## 7. Reading the tau map across specifications (`ea_apply_protocol.m`, `ea_cross_p_table.m`)

`tau_bar` is measured relative to `lambda_h`, which differs by design, so raw
`tau` is never compared across p or across systems. Each design has its own
null (`run_null_calibration`, one per `(system, p, H, h_early, psi_floor,
fixed_tau, h1_mode, chains)`) and the comparison table reports, per cell:

* `tau_bar / q95` (1 = at the 5 % threshold),
* the **null percentile** of `tau_bar` (rank among the R null draws),
* the escape indicator at q95 and the per-equation max-statistic flag,
* and per design the **count** of escapes, the count under Holm's correction
  across the K·(K−1) free cells, and the count above the max-statistic
  threshold (the bootstrap analogue of Westfall–Young; the null draws already
  supply the max distribution).

`p ∈ {2, 4, 6, 12}`; p = 4 is AGKL's lag order [L-1].

**Null modularity.** `ea_design_key.m` hashes the design (variable names,
T, p, H, window, prior options, identification is irrelevant to tau and is
recorded but not part of the key). The null file stores it;
`ea_apply_protocol` refuses a null whose key differs from the run's. No
threshold can be read against the wrong design by accident.

## 8. Fallback validation: out-of-sample block ablation (`run_block_ablation.m`)

Without ground truth, an escape is credible if releasing that block buys
predictive accuracy. Design (Giannone–Lenza–Primiceri-style pseudo
out-of-sample [L-11]): expanding window from a chosen origin; at each origin
the whole pipeline (BVAR, `lambda_h`, block sampler) is re-estimated on data
up to the origin, and the horizon-h direct forecast of `y_{i,t+h}` is formed
(deterministic component + LP fit) under (a) the block-adaptive posterior and
(b) the same posterior with cell `(i, g)` **forced to tau = 1** (a new
per-cell `cfg.blocks.fixed_tau_mask`). Loss: MSFE over h in the early
window. The comparison is between nested Bayesian predictors, so the
Clark–West (2007) adjusted statistic is reported with the Diebold–Mariano
statistic (HAC with the h-step overlap, Harvey–Leybourne–Newbold small-sample
factor) [L-14]. Prediction: forcing an escaping cell back worsens the MSFE;
forcing a quiet control cell does not. Runs today on `lev4` for the HICP
own-lag cell (escaped in v1) and a quiet control.

## 9. Legacy baseline

`RUN_EMPIRICAL.m` and every `p12*` output are unchanged and labelled
"legacy baseline (v1, internal instrument)". The new driver is
`RUN_EMPIRICAL_IV.m`; its outputs carry the prefix `iv_` and a system tag, so
nothing collides.

## 10. What runs today and what waits for Monday

Today (real data on disk): `lev4` and `lev4_yoy` systems; instruments built
from the full workbook in all variants (1M adjusted, with speeches, 1Y,
3M, JK, BS, path) with both aggregations; proxy identification with the
**legacy indicator** (12M Euribor monthly average, paired with `kilian`
aggregation) and the relevance/placebo/influence module for every variant;
the block-ablation validation on `lev4`; nulls for `lev4` at p = 2, 4, 6
(p = 12 exists); the cross-p table. Monday: `ois4`/`ois6`, their nulls, the
headline first stage with the matched indicator. See
`MONDAY_DATA_CHECKLIST.md`.
