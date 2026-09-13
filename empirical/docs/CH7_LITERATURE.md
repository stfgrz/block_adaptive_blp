# Literature behind the v2 design — what each choice rests on

Compiled 2026-09-13. Tags: **[V]** verified in the paper's own text on this
date (page/equation given); **[S]** secondary source; **[M]** from memory,
not re-verified. Where papers disagree the disagreement is stated, not
resolved. Numbered [L-n] as cited in `CH7_REDESIGN.md`.

**[L-1] Altavilla, Gürkaynak, Kind & Laeven (2025), "Monetary transmission
with frequent policy events", ECB WP 3157 (EA-EMPD).**
[V] Surprise = change in the 1-month OIS in a 30-minute window; for
speeches the day-count adjustment eqs. (1)–(3), p. 8: with m accrual days at
the old rate and 30−m at the new one, `E_t i_n − E_{t−1} i_n = 30/(30−m)
(OIS_t − OIS_{t−1})`; example "2 bp with 20 days remaining ⇒ 3 bp". [V]
Speeches outside trading hours excluded in the baseline; events whose window
contains a data release with Bloomberg relevance index ≥ 90 excluded (we
cannot replicate the latter). [V] "To construct mps_t, we aggregate intraday
surprises from days in month t that feature Governing Council meetings or
speeches ... set to zero for months without" (p. 23). [V] BVAR eqs. (7)–(8):
`x_t = [mps_t; y_t]`, mps_t first, "internal instrument" (Stock–Watson 2018,
PMW 2021), zero restrictions = no contemporaneous feedback from y_t to mps_t;
no restriction on lagged surprises is stated. Baseline y_t = industrial
production and inflation; extended set adds loan volumes and lending-rate
spreads; sample January 2001–September 2025, four lags, Minnesota prior. [V]
First stage (Fig. 9): monthly change in the 1M OIS on mps_t of the same
asset — meetings: slope 0.46 (t 2.01), R² 1.53 %, F 4.03; speeches: 1.50
(3.76), 4.64 %, F 14.12; all events: 0.66 (3.55), 3.97 %, F 12.57. Whether
the monthly 1M OIS is an average or end-of-month value is not stated in the
text we could read. [V] Predictability: R² < 1 % on a macro-news index; on
meeting days the pre-event STOXX move is the only consistently significant
predictor (Tables 3–4). [V] "The 50 bp cut in November 2008 generated a
sizeable positive surprise, whereas a similar cut in January 2009 produced a
negative one." [V] Workbook assets: OIS 1W–20Y, DE/FR/IT/ES yields, FX,
STOXX50E, SX7E (Notes sheet of `EA-EMPD.en.xlsx`); no adjusted-surprise
column is shipped, the adjustment is applied by the user.

**[L-2] Jarociński & Karadi (2020, AEJ:Macro), "Deconstructing monetary
policy surprises".** [V] Sign restrictions: rate surprise and stock surprise
co-move negatively ⇒ monetary policy, positively ⇒ central-bank information
(Table 1); "poor man's" version uses the rate surprise only in months where
the stock surprise had the opposite sign, zero otherwise, proxy ordered first
with Cholesky. [V] Zero restrictions: other shocks have no contemporaneous
impact on the surprises, "plausible as long as the surprises are
unpredictable"; results unaffected by relaxing them (Appendix). [V] Euro-area
data: EONIA swaps 1M–2Y and Euro Stoxx 50 around the 13:45 CET press release.
[M] EA VAR variables: 1-year bund yield, log HICP, log IP (or interpolated
GDP), Euro Stoxx 50, BBB spread; monthly sums of surprises. [V] Updated ECB
series: `github.com/marekjarocinski/jkshocks_update_ecb`,
`shocks_ecb_mpd_me_m.csv` (monthly) and `_d.csv`; MP_pm, CBI_pm, MP_median,
CBI_median; interest-rate surprise = first PC of OIS 1M, 3M, 6M, 1Y scaled to
the OIS-1Y monetary-event-window sd; stock = 100·Δlog Euro Stoxx 50;
CC BY 4.0, cite the AEJ paper.

**[L-3] Kuttner (2001); Gürkaynak, Sack & Swanson (2005).** [M] The
contract spanning exactly the next meeting measures the target surprise; the
day-count scaling in [L-1] is the euro-area analogue of Kuttner's fed-funds
futures scaling. Matching the instrument's maturity to the indicator is the
reason AGKL regress the 1M OIS on 1M-OIS surprises.

**[L-4] Altavilla, Brugnolini, Gürkaynak, Motto & Ragusa (2019, JME),
"Measuring euro area monetary policy" (EA-MPD).** [V] Press-release window
contains a single factor, Target; press-conference window contains Timing,
Forward Guidance and (post-2014) QE. [V] Rotation: factors 2 and 3 do not
load on the one-month OIS, "the standard measure of the immediate policy
setting surprise"; using the change in the one-month OIS directly as the
Target factor "made no difference". ⇒ a 1M-OIS instrument identifies the
conventional-policy dimension only.

**[L-5] Kilian (2024, JEDC; Dallas Fed WP 2310, 2023), "How to construct
monthly VAR proxies based on daily surprises in futures markets".** [V]
Summing daily surprises within the month "would make sense if the price in
the VAR were measured as the percent change from the last day of the
preceding month to the last day of the current month", but VAR variables are
typically monthly averages; a surprise on day d of a month with T days
raises that month's average by (T−d+1)/T and creates "an additional shock in
the following month whose magnitude depends on the timing" (the remaining
(d−1)/T); the monthly-average proxy is only weakly correlated (42 % in his
application) with the conventional sum and improves the first stage. [V]
Gertler–Karadi (2015, fn. 11) recognise the timing issue for monthly-average
indicators.

**[L-6] Stock & Watson (2018, EJ), external instruments.** [V] LP-IV needs
relevance, contemporaneous exogeneity and "lead–lag exogeneity" (uncorrelated
with past and future shocks, at least after controls); the unit-effect
normalisation (`Θ_{0,11} = 1`) gives IRFs in native units; weak-instrument
robust inference for SVAR-IV in Montiel Olea–Stock–Watson. [M] Under a
VAR(p) DGP the LP-IV estimand at h equals Γ_h b_z with b_z ∝ E[u_t z_t]
(both are the projection of y_{t+h} on the instrumented innovation).

**[L-7] Plagborg-Møller & Wolf (2021, Ecta).** [V] "Structural estimation
with an instrument (proxy) can be carried out by ordering the instrument
first in a recursive VAR, even under noninvertibility"; the internal
instrument strategy is valid under non-invertibility "unlike the well-known
external instrument SVAR-IV approach"; the LP-IV estimand equals the
recursive VAR with the IV ordered first.

**[L-8] Montiel Olea & Pflueger (2013, JBES); Montiel Olea, Stock &
Watson (2021, J. Econometrics).** [V] Effective-F critical values for
K_eff = 1: 12.05 (τ = 30 %), 15.06 (20 %), 23.11 (10 %), 37.42 (5 %); "a
simple rule of thumb ... rejects when the effective F is greater than 23.1".
[M] With one instrument the effective F equals the HAC-robust first-stage F.
[V] MOSW: Anderson–Rubin set `{β : q(β) ≤ χ²_{1,1−a}}` valid for any
instrument strength; report it together with the robust first-stage F or
Wald statistic. [V] Ramey (2016) quotes the same thresholds (23 for 10 %,
37 for 5 %) and notes the Stock–Yogo 10 applies only with serially
uncorrelated first-stage errors.

**[L-9] Bauer & Swanson (2023, NBER Macro Annual; 2023 AER).** [V] Six
predictors: nonfarm-payrolls surprise, 12-month employment growth, S&P 500
return over the prior 3 months, change in the yield-curve slope, commodity
price change, Treasury implied skewness; surprises orthogonalised on them;
predictability attributed to the "Fed response to news" channel rather than
an information effect — a disagreement with [L-2] and [L-12] on
interpretation, not on the regression. Euro-area analogue with what is on
disk: prior-3-month STOXX return, 12-month IP growth, 12-month inflation,
prior-3-month change in the indicator.

**[L-10] Gertler & Karadi (2015, AEJ:Macro).** [V] Monthly-average policy
indicators; "a surprise at the end of a month can be expected to have a
smaller influence on the monthly average rate than a surprise at the
beginning"; monthly-average surprise construction (fn. 11, 31-day window
[M]).

**[L-11] Giannone, Lenza & Primiceri (2015); Bańbura, Giannone & Reichlin
(2010).** [M] Pseudo out-of-sample evaluation of shrinkage priors with
expanding windows and MSFE by horizon; random-walk centre for levels and
zero centre for stationary series.

**[L-12] Miranda-Agrippino & Ricco (2021, AEJ:Macro; BoE SWP 657).** [V]
Instrument = projection of market surprises on their own lags and on
Greenbook forecasts (the central bank's information set); high-frequency
surprises are autocorrelated and predictable. [M] Monthly sums; external
instrument in a BVAR.

**[L-13] Ricco, Savini & Tuteja (2024/2026), "Monetary policy, information
and country risk shocks in the euro area", CEPR DP 19679.** [V] Four
factors (target, forward guidance, QE/QT, country risk) from EA-MPD
press-release + press-conference responses summed; information correction
"survey-based" (projection on Eurosystem/ECB staff projections, published
after the meeting) versus "market-based" (JK, Cieslak–Schrimpf); linear
corrections "insufficient", non-linear information effects concentrated in
a few high-volatility crisis events explain up to ~40 % of short-maturity
price revisions and resolve the output/price/stock puzzles. [V] No public
release of the series is mentioned in the paper or on the authors' pages
checked (Savini's site lists papers only).

**[L-14] Diebold & Mariano (1995); Harvey, Leybourne & Newbold (1997);
Clark & West (2007).** [M] DM with HAC for h-step overlapping errors and the
HLN small-sample factor; Clark–West adjusts for the noise of the larger
nested model. Applied to the block ablation (nested by construction).

**Series keys verified on the ECB Data Portal API on 2026-09-13** (one
observation requested, nothing stored): `ICP.M.U2.Y.000000.3.INX` HICP
overall index, ECB, working-day and seasonally adjusted, 2015 = 100
(`...S...` does not exist); `BSI.M.U2.Y.U.A20T.A.I.U2.2240.Z01.E` adjusted
loans to NFCs, index of notional stocks, SA, from 2003-01;
`BSI.M.U2.N.A.A20.A.1.U2.2240.Z01.E` loans to NFCs, outstanding amounts,
NSA, from 1997-09; `MIR.M.U2.B.A2I.AM.R.A.2240.EUR.N` composite cost of
borrowing, NFCs, from 2003-01; `FM.M.U2.EUR.4F.MM.EONIA.HSTA` EONIA
monthly average 1994–2021; `FM.M.U2.EUR.RT.MM.EURIBOR1MD_.HSTA` 1M Euribor
monthly average from 1994; `FM.M.U2.EUR.RT.MM.EURIBOR1YD_.HSTA` (in use).
An end-of-month STOXX key (`...HSTE`) does not exist; the STOXX stays a
monthly average.
