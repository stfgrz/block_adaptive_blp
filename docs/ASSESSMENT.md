# Critical assessment of the empirical chapter (2026-09-15)

An outside reading of the euro-area application, written after
reorganising the repository, re-deriving the headline numbers from the
stored files and running the additional analyses described in Section 3.
The simulation study is not re-assessed here beyond what it implies for
the empirical protocol; its own verdict is in `RESULTS.md` §5.

## 1. What changed in the repository, and why

* **Layout.** `results/` was a flat directory of some 240 files mixing
  simulation output, two empirical designs and their nulls. It is now
  `results/simulation/`, `results/empirical/` (current design) and
  `results/empirical/legacy_v1/`, with figures in `figures/` subfolders and
  a `results/README.md` that says what each file is and which number to
  quote. The v1 empirical code moved to `empirical/legacy_v1/` and remains
  runnable; `ea_paths` carries the new locations and now points
  `dataset` at the headline system. Every code path, script and document
  reference was updated and the full test suite passes.
* **Documentation.** The root README was a 300-line status journal with
  dated banners; it is now a map. The method moved to `docs/METHOD.md`,
  the history to `docs/CHANGELOG.md`. The empirical documents were renamed
  to their permanent roles (`DESIGN`, `DESIGN_V1_LEGACY`, `LITERATURE`,
  `DATA`, `RESULTS`, `RESULTS_LOG`); the 830-line results journal is kept
  verbatim as the log, and a short `RESULTS.md` states the current reading.
* **Hygiene.** LaTeX build artefacts and a tracked 2.7 MB Refinitiv-sourced
  spreadsheet were untracked (the daily OIS csv built from it was already
  ignored as third-party data; the spreadsheet should not have been
  committed either). `.claude/` is ignored. Shell scripts no longer
  hard-code another machine's path.
* **Code.** Three methodological changes, each pinned by tests:
  1. `ea_apply_protocol` reports Monte Carlo p-values `(r + 1)/(R + 1)` and
     a family-wise max-statistic p-value, runs Holm on the former and says
     when `R` is too small for Holm to reject at all.
  2. `simulate_fitted_bvar_dgp` gained a wild-bootstrap (and a
     block-bootstrap) residual scheme; `run_null_calibration` names the
     scheme in its output stem; `RUN_EMPIRICAL_IV` reads every null it
     finds for a design and re-runs the pooled estimator at the null's
     horizon count.
  3. `ea_relevance_iv` reports where the first-stage covariance comes
     from and whether the instrument is predictable from the VAR's lags.

## 2. The main weaknesses of the empirical application

Ordered by how much they threaten the chapter's claims.

1. **Inference on the tau map was over-stated.** The 2026-09-14 reading
   called the short-rate escape in the output equation "Holm-significant".
   With `R = 200` null draws and 16 cells, Holm's first threshold
   (0.05/16 = 0.0031) lies below the smallest attainable Monte Carlo
   p-value (1/201 = 0.0050); the rejection came from treating an empirical
   p of exactly zero as zero. No Holm rejection was possible with that
   null. This is a correctable error, and it is corrected, but it is the
   kind a referee finds first.
2. **The null was homoskedastic in a sample that is not.** The standard
   deviation of the short-rate innovation is 0.13–0.16 pp in 2001–2011 and
   0.04 pp in 2012–2019. The iid resampling null spreads crisis-sized
   innovations over the whole sample. Under a wild-bootstrap null that keeps
   each month's own residual, the null distribution of the short-rate block's
   scale widens where the rate's own volatility matters: at `p = 2` the
   output-equation cell's 95th percentile rises from 1.84 to 2.47 and its
   family-wise p-value from 0.020 to 0.113; at `p = 12` the 95th percentile
   of the short-rate block rises by about a fifth in the rate's own and the
   stock equation (the stock-equation cell: family-wise p 0.138 / 0.039 for
   the independent / pooled estimator under the iid null, 0.224 / 0.098
   under the wild one) and is unchanged in the output equation, whose cell
   survives (ratio 1.9–2.1, family-wise p 0.002–0.006 under both nulls,
   both estimators). Heteroskedasticity
   therefore explains the peripheral escapes but not the central one; the
   diagnostic as specified still cannot tell a change in volatility from a
   change in dynamics, and the wild null has to be the pre-specified one.
3. **Identification is a 2009–2011 event.** The headline first stage
   (F_eff 13.7, LP-IV 19.4) clears the 30 % but not the 20 % Nagar-bias
   threshold, and 63 % of the identifying covariance comes from 2009–2011,
   53 % from five months, half of the absolute mass from ten months. The
   2012–2019 short end contributes 7 %. The IRFs are the responses to the
   sovereign-debt-crisis surprises, and their information-effect pattern
   (stocks and output up after a tightening surprise) is what one expects
   from that period. The pure-policy instruments that would remove it have
   F between 2 and 3.
4. **Exogeneity is imperfect.** A rate innovation this month predicts a
   negative surprise next month (t = −2.4), and the surprises are partly
   predictable from the VAR's 48 lags on event months (R² 0.27). The proxy
   route conditions on the lags, so the second is admissible, but both say
   the instrument is not news relative to the model's information set.
5. **The substantive finding is thin.** With seasonally adjusted prices,
   the calibrated escapes are the short-rate block in the output (and
   marginally the stock) equation; they carry no out-of-sample forecast
   content; they appear at `p = 2` and `p = 12` but not at `p = 4` or `6`,
   which is not the pre-registered dose–response and has only a post-hoc
   explanation. The one escape that is validated on every count (NSA HICP
   seasonality in the level system) is a mechanical misspecification of a
   random-walk-centred VAR(12) prior for a seasonal index: a textbook
   positive control, not a fact about monetary transmission.
6. **The IRF bands for the policy rate are not credible.** Under the null
   the sandwich bands cover the true rate response 41 % of the time at a
   nominal 90 % (59 % for output; 54 % and 64 % under the wild null); the
   levels-LP bias toward zero at
   intermediate horizons for a near-unit-root response is not something a
   band centred on the shrunk estimate can absorb. The lag-augmented LP is
   in the repository and should be the interval shown for that response.
7. **Design and data limits.** The stock index is a monthly average, so
   its impact response is smeared over two months; the sample stops in
   2019, before the only period (2022–23) in which the 1-month instrument
   has variance outside a crisis; `K = 4` has no credit or spread variable;
   the loans-and-spread system starts in 2003 and is too thin at `p = 12`.
8. **The invertibility assumption of the proxy route is untested.** The
   v1 internal-instrument design is called the robustness check, but it
   used a different instrument and indicator, so the two are not a
   like-for-like comparison.

## 3. What was investigated, and what was found

All on the stored `ois4` files and the real data, in MATLAB, 2026-09-15.

| question | method | finding |
|---|---|---|
| Does the Holm claim survive the standard Monte Carlo p-value? | recompute `p = (r + 1)/(R + 1)` from the stored `R = 200` null draws | No: zero Holm rejections for either estimator; 4 (block) / 5 (pooled) of 16 cells have unadjusted `p ≤ 0.05` against 0.8 expected. |
| What does a proper family-wise test say? | max over free cells of `tau_bar / q95` in each null draw (single-step max-T) | Under the `R = 200` iid null the output-equation short-rate cell has `p_fwer` 0.010 (block) / 0.020 (pooled), the stock-equation cell 0.050 / 0.025; nothing else. |
| Is the null's homoskedasticity doing work? | wild-bootstrap null (Rademacher signs on each month's residual), `R = 300` at `p = 2, 4, 6`, `R = 500` at `p = 12` | At `p = 2` the escape disappears (`p_fwer` 0.113 / 0.146); at `p = 12` the output-equation cell survives (ratio 2.01 / 2.12, `p_fwer` 0.006 / 0.006) and the stock-equation cell does not (0.224 / 0.098); at `p = 4, 6` nothing under either null. Coverage of the rate's own response under the wild null improves from 41 % to 54 %. |
| Does `R` matter for the iid reading? | iid null rerun with `R = 1000` at `p = 12` | The 95th percentiles move by a few hundredths; the output-equation cell keeps `p_fwer` 0.002 / 0.006 and is now also Holm-rejected (smallest attainable p 0.001 < 0.0031); the stock-equation cell's borderline 0.050 at `R = 200` becomes 0.138 (block) / 0.039 (pooled): the `R = 200` boundary reading was resolution noise. |
| Where does the first stage come from? | signed covariance shares by year, sub-period and month; predictability from the VAR lags | 63 % from 2009–2011, 29 % from 2001–2008 (two months of 2001 alone 25 %), 7 % from 2012–2019; ten months carry half the absolute mass; R² 0.27 on the lags. |
| Does the pooled estimator's protocol depend on `H`? | re-run at `H = 12` vs the stored `H = 48` | Cells move by up to 0.3 in `tau_bar` (e.g. 0.43 → 0.59, 2.04 → 1.74); the family-wise reading of the top two cells is unchanged. The independent estimator is `H`-invariant draw for draw. |
| Do the sandwich bands cover under the null? | `eq_coverage_under_null` in the null files | iid null: 41 % for the rate's own response, 59 % output, 75 % prices, 77 % stocks; wild null: 54 / 64 / 78 / 83 % (nominal 90 %). |

## 4. Is there a publication-worthy contribution?

**Not as an empirical paper about euro-area transmission.** The IRFs are
the well-documented information-effect pattern, identified from three
crisis years with an instrument of moderate strength, and the chapter is
right not to claim a transmission finding.

**Possibly, as a methods paper with a demonstration**, provided the
framing is the diagnostic and the demonstration is presented for what it
is. The defensible core is:

1. an adaptive shrinkage layer exactly nested in a published estimator,
   which makes the diagnostic free to add to any FMAR application;
2. a design-specific calibrated null with family-wise inference, and the
   demonstration that the null's assumptions matter: the same map reads
   differently under a homoskedastic and a heteroskedasticity-robust null,
   which is a useful warning for anyone reading posterior scales as
   evidence of misspecification;
3. a clean positive control on real data (seasonality in NSA prices:
   found at every lag order, validated out of sample, gone when fixed) and
   an instructive borderline case (the short-rate block of the output
   equation: a calibrated disagreement that survives both nulls at
   `p = 12`, is absent at `p = 4` and `6`, and carries no forecast
   content);
4. the simulation evidence that the diagnostic localises and that
   localisation does not imply improvement.

That is a coherent thesis chapter and a plausible field-journal
methods note. It is not yet a paper that a referee at a general-interest
journal would accept, because the empirical demonstration has no result
that is both novel and robust.

## 5. The most promising improvements

In order of expected return.

1. **Extend the sample to 2024 and treat 2020.** The 2022–23 hiking cycle
   is where a 1-month OIS surprise has variance outside a crisis (AGKL's
   first stage is built on it). Needs: the same ECB and Eurostat series to
   2024 (free; the daily OIS export already runs to 2026), EA-EMPD events
   (shipped to 2025), and a volatility treatment for 2020m3–m12 in the BVAR,
   the LPs and the Newey–West steps (Lenza–Primiceri rescaling or pandemic
   dummies). This is the single change most likely to turn the first stage
   from a 2009–2011 phenomenon into a sample-wide one.
2. **Make the wild-bootstrap null the pre-specified one** and add a null
   that draws the VAR parameters from their posterior, so that the
   prior-centre uncertainty (16 % of the band width) enters the null. Both
   are implementable with the existing code in a few hours.
3. **A stronger and cleaner instrument.** The 3-month surprise has F 19.2
   and the orthogonalised 1-month surprise F 15.0 with a quieter lead
   placebo; a factor across maturities would combine them. For the
   information effect, request the Ricco–Savini–Tuteja factors or build a
   Jarociński–Karadi split that keeps power (speeches included).
4. **Fix the stock-index timing** with daily EURO STOXX 50 closes (freely
   available) aggregated to end of month or to the event day, and add a
   credit variable (a euro corporate bond spread) so the system speaks to
   the transmission literature.
5. **Show the rate's own response with the lag-augmented LP** (already in
   the repository) and report Anderson–Rubin sets as the primary IRF
   intervals while F_eff is below 23.
6. **Test invertibility like for like**: run the internal-instrument design
   with the *same* 1-month instrument and indicator (`ois1m` system with
   the surprise ordered first) and compare the impact vectors.
7. **Validate escapes with a forecast design that covers the crisis**
   (rolling windows at `p = 6`, or density-forecast scores), since the
   current 2010–2019 window never has to forecast the period that carries
   the disagreement.

## 6. Changes implemented today and their effects

| change | effect on the reading |
|---|---|
| Monte Carlo p-values and the family-wise max test | The "Holm-significant" headline becomes "one cell at family-wise 5 % under the iid null"; the protocol now states the resolution of its null. |
| Wild-bootstrap null (`p = 2, 4, 6, 12`) | Removes the `p = 2` escape and the stock-equation cell at `p = 12`; leaves the output-equation cell at `p = 12` (family-wise p 0.006). The reading becomes one robust cell instead of two or three. |
| iid null at `R = 1000` (`p = 12`) | Holm becomes attainable and agrees with the max-statistic test on the one cell; the stock-equation cell's `R = 200` boundary reading turns out to be noise (0.138 block / 0.039 pooled). |
| Pooled estimator re-run at the null's `H` | Removes an approximation the docs had flagged; changes individual cells by up to 0.3 but not the top-two reading. |
| First-stage anatomy and lag predictability | Quantifies what the sub-sample F statistics only hinted at; belongs in the chapter's identification table. |
| QUICK mode tolerates a mismatched null key | The interface check runs end to end again (it errored once a null existed). |

## 7. Recommended next steps, by priority

1. Rewrite the chapter's results section from `empirical/docs/RESULTS.md`:
   the diagnostic finds the prior adequate once seasonality is handled;
   the one disagreement that survives both nulls (short rate → output at
   `p = 12`) has no forecast content and is absent at intermediate lag
   orders; the IRFs are conditional on a 2009–2011 instrument. Drop every
   Holm statement from the 2026-09-14 reading.
2. Pre-specify the wild-bootstrap null and rerun the four lag orders with
   `R ≥ 500` each (about 3 hours of machine time); regenerate the cross-p
   table.
3. Extend the sample to 2024 with a COVID treatment and rerun the first
   stage; if F_eff clears 23 on the full sample, the IRF exhibit becomes
   reportable.
4. Replace the rate's own IRF bands by the lag-augmented LP bands; move the
   Anderson–Rubin sets into the main IRF figure.
5. Fetch daily STOXX closes and build an end-of-month index; add a
   corporate spread; rebuild `ois4` / `ois6`.
6. Run the like-for-like invertibility check (internal instrument with the
   1-month surprise).
7. Add the posterior-draw null and a subsample stability check of the tau
   map (2001–2011 vs 2012–2019 at `p = 6`).
