# Project history

Dated record of the decisions that shaped the repository. Newest first.
Detailed numbers are in `RESULTS.md` (simulations) and
`../empirical/docs/RESULTS_LOG.md` (empirical chapter, chronological).

## 2026-09-15 — reorganisation, protocol correction, heteroskedasticity-robust null

* **Layout.** `results/` split into `simulation/` (Monte Carlo runs, grid,
  sensitivity, report, figures), `empirical/` (current design: `iv_*`,
  nulls, ablations, figures) and `empirical/legacy_v1/` (the 2026-09-12
  internal-instrument outputs). The v1 empirical code moved to
  `empirical/legacy_v1/` (`RUN_EMPIRICAL_V1`, `SMOKE_TEST_EMPIRICAL_V1`,
  `assemble_dataset`, `build_shock_series`, `ea_relevance`) and stays
  runnable. `ea_paths` gained `results_legacy`, `figures`, `dataset_legacy`;
  `ea_paths().dataset` is now the `ois4` headline dataset. Build artefacts
  (`presentation/build/`) and the Refinitiv-sourced daily OIS export were
  untracked; `.claude/` ignored. Shell scripts no longer hard-code a machine
  path. Documentation regrouped: `docs/METHOD.md`, this file,
  `docs/ASSESSMENT.md`; `empirical/README.md`; `empirical/docs/{DESIGN,
  DESIGN_V1_LEGACY, LITERATURE, DATA, RESULTS, RESULTS_LOG}.md`.
* **Protocol.** `ea_apply_protocol` reports Monte Carlo p-values
  `(r + 1)/(R + 1)` (Holm now runs on them and states when `R` is too small
  to reject at all) and a family-wise max-statistic p-value over the free
  cells; `ea_cross_p_table` carries the new counts and the null's residual
  scheme. The "Holm-significant" readings of 2026-09-14 rested on an
  `R = 200` null in which no Holm rejection was attainable; the headline
  iid null was rerun with `R = 1000`.
* **Null.** `simulate_fitted_bvar_dgp` gained `'wild'` (Rademacher signs on
  each month's own residual) and `'block'` schemes; `run_null_calibration`
  appends the scheme to the output stem. The short-rate innovation's
  standard deviation falls fourfold after 2012; under the wild null the
  short-rate block's 95th percentile rises by about a fifth in the rate's
  own and the stock equation and is unchanged in the output equation, so
  the stock-equation and the `p = 2` escapes disappear while the
  output-equation escape at `p = 12` survives (family-wise p 0.006); see
  `empirical/docs/RESULTS.md`.
* **Relevance.** `ea_relevance_iv` reports the anatomy of the first stage
  (covariance shares by year, sub-period, top months) and the instrument's
  predictability from the VAR lags. Two thirds of the headline first stage
  comes from 2009–2011.
* **Pooled protocol.** `RUN_EMPIRICAL_IV` re-runs the pooled estimator at
  `H = 12` for the protocol (the null is simulated at that `H`); the
  independent estimator's `tau` at `h ≤ 12` is `H`-invariant draw for draw.
* **QUICK mode** tolerates a mismatched null key (it is an interface
  check, and the key differs by construction).

## 2026-09-14 — the 1-month OIS data

The daily 1M OIS (ECB Data Portal copy of the Refinitiv series), five ECB
series (SA HICP, loans, cost of borrowing, EONIA, 1M Euribor) and
Jarociński's updated shocks were obtained. Headline system `ois4 = {ois1m
eom, ip, hicp_sa, stoxx}` with the matched instrument `z_gcs_1m_adj_sum`:
first stage F_eff 13.7 (LP-IV 19.4), strength concentrated in 2009–2011;
IRFs show the euro-area information-effect pattern. Nulls at
`p ∈ {2, 4, 6, 12}` with `R = 200`; ablations. Log: `RESULTS_LOG.md`, "Monday run".

## 2026-09-13 — external-instrument redesign (v2)

The surprise left the VAR state vector and became an external instrument
used only for the impact vector (`ea_identify_proxy`); instruments rebuilt
from the full EA-EMPD workbook (1M adjusted per AGKL, speeches, 3M, 1Y, JK,
path proxy; `sum` / `kilian` / `kilianfix` aggregations);
`ea_relevance_iv`, `estimate_lp_iv` (Anderson–Rubin), `ea_design_key` /
`ea_apply_protocol` / `ea_cross_p_table`, `run_block_ablation`. Literature
behind each choice verified in the papers' text (`LITERATURE.md`). Dress
rehearsal on the legacy indicator: the HICP own-lag escape of the NSA level
system carries out-of-sample content (+56% two-month MSFE when ablated).

## 2026-09-12 — first run on real data (v1)

Four outcome series validated (a set of placeholder files shipped earlier
was caught by `ea_check_series` and quarantined). `RUN_EMPIRICAL` steps A–F
and the null calibrations ran in MATLAB. Three findings forced by the data:
the 1-year surprise is a weak instrument for the monthly-average 12-month
Euribor (robust t 0.77); FMAR's Newey–West prior scale collapses on the
near-white-noise surprise (`cfg.fmar.psi_floor`); the surprise's lag block
violates the group prior's common-scale assumption (`cfg.blocks.fixed_tau`).
Calibrated escapes: the HICP own-lag block (seasonality) and the
interest-rate block in the IP equation.

## 2026-09-09 — headline Monte Carlo runs complete

`R = 500` (correct, sparse) and `R = 250` (intermediate, dense) with five
estimators; paired comparison; the early/late split; horizon pooling
RMSE-neutral; localisation 0.97 against a 0.37 null; the intermediate DGP
reversal. Contribution reframed as a calibrated localised diagnostic
(`RESULTS.md` §5).

## 2026-09-08 — fair benchmark, pooled estimator, diagnostics

`h1_mode = 'lp'`, integrated RMSE over `h = 2..H`, Rao-Blackwellised point
estimate, `estimate_blp_blockpooled`, `tau_diagnostics`, exact parallel
chunking, tidy CSV exports, the 13-cell exploratory grid.

## 2026-09-01 — Chapter 7 design and data pipeline

Design document with pre-registered expectations (`DESIGN_V1_LEGACY.md`),
EA-EMPD shock construction, outcome-data fetchers, null calibration, the
lag-augmented LP comparator.
