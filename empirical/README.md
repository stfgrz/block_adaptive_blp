# Empirical chapter — how to run it

The euro-area application of the block-adaptive BLP. Design:
`docs/DESIGN.md`; data sources and how to obtain them: `docs/DATA.md`;
literature behind each choice: `docs/LITERATURE.md`; current results:
`docs/RESULTS.md`; chronological record of every run: `docs/RESULTS_LOG.md`.

Every entry point puts the repository on the path itself and resolves its
paths through `ea_paths.m`, so any of them runs from any working directory:

```matlab
addpath(genpath('<repo>'))          % once per session (matlab -batch: absolute path)
```

## What the design is, in five lines

Four monthly euro-area series 2000m1–2019m12 (`ois4`: end-of-month 1-month
OIS rate, industrial production, SA HICP, EURO STOXX 50), a Bayesian VAR(12)
that supplies the prior centre, and the three BLP estimators. The
high-frequency 1-month OIS surprise around ECB meetings and speeches
(EA-EMPD, day-count adjusted, summed within the month) is an **external
instrument**: it never enters the VAR; it identifies the impact vector once
(`data/ea_identify_proxy.m`). The `tau` map of the two adaptive estimators
is read against nulls simulated from the fitted VAR for exactly this design
(`montecarlo/run_null_calibration.m`, `ea_apply_protocol.m`), compared
across lag orders (`ea_cross_p_table.m`) and validated out of sample by
forcing an escaping cell back to the global prior
(`montecarlo/run_block_ablation.m`).

## Files

```
RUN_EMPIRICAL_IV.m          driver: steps REL (relevance), A (estimation + IRFs), B (tau map + protocol),
                            C (LP-IV table), D (cross-p), E (ablation); QUICK = true is the interface check
ea_paths.m                  absolute locations (single source of truth)
ea_design_key.m             canonical design identifier a null is keyed to
ea_apply_protocol.m         the reading protocol: ratio to null q95, Monte Carlo p, family-wise p, Holm, counts
ea_cross_p_table.m          the protocol across designs / nulls
data/
  build_instrument_series.m   EA-EMPD workbook extract -> 51 monthly instruments (bp)
  import_ois_daily.m          daily 1M OIS -> monthly eom / avg / first
  fetch_outcome_data.m        Eurostat IP, NSA HICP, ECB STOXX, 12M Euribor (download + validation)
  ea_fetch_v2_series.m        ECB SA HICP, loans, cost of borrowing, EONIA, 1M Euribor
  import_external_instrument.m  any monthly csv (e.g. Jarocinski's shocks) as an extra instrument
  assemble_dataset_v2.m       -> data/ea_dataset_v2_<system>.mat (Y, instruments Z, pre-event info)
  ea_identify_proxy.m         proxy impact vector, HAC se, first stage; drops into every estimator
  ea_relevance_iv.m           first stage, placebos, predictability, influence, subsamples, ANATOMY
  ea_pre_event_info.m, ea_check_series.m, ea_extract_series.m, ea_write_provenance.m, read_sdmx_csv.m
estimators/  estimate_lp_iv.m (HAC, lag-augmented EHW, Anderson-Rubin sets), estimate_lp_lagaug.m
dgp/         simulate_fitted_bvar_dgp.m  the null: 'resample' | 'wild' | 'block' | 'gaussian' residuals
montecarlo/  run_null_calibration.m (checkpointed), run_block_ablation.m
legacy_v1/   RUN_EMPIRICAL_V1, SMOKE_TEST_EMPIRICAL_V1, assemble_dataset, build_shock_series, ea_relevance
tests/       make_synthetic_fixture, make_synthetic_ois_daily; test_empirical_pipeline, test_empirical_iv,
             test_null_modularity, test_lp_iv (all in tests/run_all_tests)
```

Outputs: `results/empirical/` (`results/README.md` lists every file).

## Run order

```matlab
build_instrument_series();                 % seconds; prints the reference statistics (221 GC events, sd 4.29 bp, ...)
fetch_outcome_data();  ea_fetch_v2_series();   % downloads + validates (or verifies files already on disk)
import_ois_daily();                        % needs raw/ois1m_ea_daily.csv (Date, Bid, Ask); see docs/DATA.md
jk = struct('path', fullfile(ea_paths().raw, 'shocks_ecb_mpd_me_m.csv'), 'column', 'MP_pm', 'name', 'jk_mp_pm', 'scale', 100);
assemble_dataset_v2(struct('system', 'ois4', 'external', {{jk}}));    % T = 240, K = 4, 106 instruments
QUICK = true; SYSTEM = 'ois4'; run(fullfile(ea_paths().empirical, 'RUN_EMPIRICAL_IV.m'))   % ~4 min interface check
SYSTEM = 'ois4'; run(fullfile(ea_paths().empirical, 'RUN_EMPIRICAL_IV.m'))                 % ~30 min; steps REL A B C D E
```

Nulls (separate jobs; ~6 s per replication at `p = 12`, 3 s at `p = 2`;
checkpointed every 10 replications and resumable):

```matlab
for p = [2 4 6 12]
    run_null_calibration(struct('p', p, 'null', struct('n_rep', 500, 'stem', sprintf('iv_ois4_p%d', p))));
    run_null_calibration(struct('p', p, 'null', struct('n_rep', 500, 'sim_method', 'wild', 'stem', sprintf('iv_ois4_p%d', p))));
end
SYSTEM = 'ois4'; DO_REL = false; DO_A = false; DO_C = false; DO_E = false; run(fullfile(ea_paths().empirical, 'RUN_EMPIRICAL_IV.m'))
```

Step B applies the protocol to every null found for the design
(`NULL_VARIANTS`, default the iid null and the `_wild` / `_block` ones) and
re-runs the pooled estimator at the null's `H` so the comparison is exact.
Holm's step-down over 16 cells needs `R >= 319` to be able to reject at all;
the headline iid null was run with `R = 1000`.

Other systems: `assemble_dataset_v2(struct('system', 'ois4_yoy', 'vars', {{'ois1m', 'ip', 'hicp_yoy', 'stoxx'}}))`,
`'ois4avg_nsa'` with `'ois1m_avg'` (the Kilian aggregation test), `'ois6'`
(loans and lending spread, 2003 start), `'lev4'` / `'lev4_yoy'` (the 12M
Euribor indicator). The headline instrument follows a pre-specified rule
(maturity matched to the indicator, meetings + speeches, aggregation matched
to the indicator's timing convention); `HEADLINE_Z` overrides it.

Toggles preset in the workspace before `run`: `DO_REL`, `DO_A` ... `DO_E`,
`DO_POOLED`, `P_BASE` (12), `H_IRF` (48), `PLIST` ([2 4 6 12]),
`HEADLINE_Z`, `ALT_Z`, `NULL_FILE`, `NULL_VARIANTS`, `POOLED_PROTOCOL_EXACT`,
`N_ORIGINS`, `ORIGIN0_YM`, `DATASET`, `QUICK`.

## Reading the output

* `iv_ois4_p12_relevance.csv`: one row per instrument. `F_eff` against the
  Montiel Olea–Pflueger thresholds 12.05 / 15.06 / 23.11 / 37.42; the lead
  placebos (`lead1_t` should be insignificant); `cov_share_2009_11`,
  `cov_share_top5`, `n_half_abs` say where the identifying covariance comes
  from; `lagpred_p` and `predchg_t` whether the instrument is predictable
  from the VAR's own past.
* `iv_ois4_p12_tau_protocol.csv` (iid null) and `_tau_protocol_wild.csv`:
  per cell `ratio_q95`, `p_mc`, `p_fwer`. Report a cell when `p_fwer <= 0.05`
  under the pre-specified null and say whether it holds under the other.
* `iv_ois4_cross_p.csv`: escape counts and minimum family-wise p per lag
  order and null.
* `ablation_iv_ois4_p12.csv`: forcing an escaping cell back to `tau = 1`;
  positive Clark–West means the escape carries forecast content.
* IRFs (`iv_ois4_p12_irf.csv`, `figures/fig_iv_ois4_p12_irf.png`) only
  together with the first stage; below the 10%-bias threshold report the
  Anderson–Rubin sets in `iv_ois4_p12_lpiv.csv`.

## Data integrity

`fetch_outcome_data`, `ea_fetch_v2_series` and `assemble_dataset_v2`
validate every series with `ea_check_series` and **refuse** to proceed on
data that fails: coverage of the estimation window, and three fingerprints
that catch generated series (the package once shipped four placeholder files
under the real names; they are quarantined in `data/raw/placeholder_rejected/`).
The tripwire is the roughness ratio `sd(diff²y)/sd(diff y)`, about 1.4 for a
real monthly macro series and 0.00–0.18 for the placeholders. Every accepted
raw file gets a `<file>.source.txt` sidecar naming the series key and the
download, and the dataset records those lines in `ds.meta.sources`.

## Legacy v1

The 2026-09-12 design (surprise inside the VAR as the first variable, 12M
Euribor indicator, prior-scale floor, surprise block held at `tau = 1`) is
in `legacy_v1/`, writes to `results/empirical/legacy_v1/`, and is described
in `docs/DESIGN_V1_LEGACY.md`. It remains the robustness check for
invertibility, which the external-instrument route assumes.
