# Reproducing the analysis

Every number in the write-up comes from one of the commands below.
Nothing depends on the state of a MATLAB/Octave workspace: each entry
point resolves its own paths and takes its settings from a **named
preset** (`montecarlo/mc_preset.m`), which is recorded inside every
result it produces (`mc.meta.preset`).

Requirements: MATLAB R2018b+ or GNU Octave ≥ 8. No toolboxes. The
simulations use plain `for` loops; parallelism, where used, is separate
operating-system processes, not a toolbox.

---

## 0. Test suite (run this first)

```matlab
cd tests
run_all_tests
```

About 12 minutes in MATLAB (25 in Octave); 21 tests. `test_fmar_port` is skipped unless `FMAR_PATH` points
at the `demo_IRFs` folder of the FMAR replication package.

---

## 1. Single-dataset demonstrations

```matlab
RUN_FMAR_DEMO                      % FMAR stack, all five estimators
QUICK_DEMO = true; RUN_FMAR_DEMO   % ~2 minutes
RUN_ME_FIRST                       % the original prototype stack, unchanged
```

`RUN_FMAR_DEMO` prints the selected `lambda` path, the posterior block
scales `tau` under both adaptive estimators, the exact `tau = 1` nesting
check, and (with `cfg.demo.run_montecarlo = true`) a small Monte Carlo
report.

**Side effect.** Both demonstrations write into `results/simulation/`:
`RUN_ME_FIRST` overwrites `mc_sparse.mat` and `figures/fig1`–`fig5`,
`RUN_FMAR_DEMO` overwrites `mc_fmar_{sparse,correct}.mat`,
`figures/fig_fmar_irf_sparse.png` and writes `mc_fmar_*` CSVs. Those
`.mat`/`.png` files are tracked as the earliest baselines, so after a demo
run restore them with
`git checkout -- results/simulation/mc_fmar_*.mat results/simulation/mc_sparse.mat results/simulation/figures/`
and delete the CSVs, unless you mean to replace them.

---

## 2. Exploratory experiment grid

Thirteen cells sweeping misspecification strength, sample size, fitted
lag order and sparsity, at `R = 40` (see `montecarlo/experiment_grid_spec.m`
for the design and for why it is one-factor-at-a-time).

```bash
scripts/run_grid_parallel.sh 4          # 4 worker processes, ~1.5 h
```

or, inside MATLAB/Octave:

```matlab
run_experiment_grid                     % all cells, serially
run_experiment_grid(3:5)                % a slice, for manual parallelism
collect_grid_results                    % -> results/simulation/grid/grid_summary.csv
```

Outputs land in `results/simulation/grid/`: one `.mat` and one set of CSVs per
cell, plus `grid_summary.csv` with the per-cell verdict (improves RMSE /
early horizons / diagnostic only / no signal).

**The grid is exploratory.** At `R = 40` the Monte Carlo standard error
of an RMSE is about 11% of it, so a cell locates a large systematic
effect and nothing finer.

---

## 3. Final Monte Carlo runs

`R = 500` at the headline design (`K = 3`, `p = 2`, `T = 200`, `H = 20`,
FMAR mode, `h1_mode = 'lp'`, all five estimators), split over four
processes and merged exactly:

```bash
scripts/run_final_parallel.sh sparse       500 4
scripts/run_final_parallel.sh correct      500 4
scripts/run_final_parallel.sh intermediate 250 4
scripts/run_final_parallel.sh dense        250 4
```

These take hours. Start them **detached**, or an interruption of the
launching shell takes the workers with it (which is exactly what
happened once here — the workers shared the shell's process group and
were killed mid-run):

```bash
setsid nohup scripts/run_remaining_finals.sh > finals.log 2>&1 < /dev/null &
```

Progress is easiest to follow through the chunk files appearing in
`results/simulation/chunks/`; Octave block-buffers stdout when it is redirected, so
the per-chunk logs lag well behind the actual work.

Each call writes `results/simulation/mc_headline_<dgp>.mat` plus the CSV exports.
Roughly 3 hours per `R = 500` DGP on 4 cores.

Serial equivalent (identical results, ~11 hours per DGP):

```matlab
cfg = mc_preset('final');
mc  = run_montecarlo(cfg, 'sparse');
s   = summarize_montecarlo(mc);
save(fullfile('results', 'simulation', 'mc_headline_sparse.mat'), 'mc', 's', '-v7');
```

The chunked and serial runs are **bit-identical**: replication `r`'s
seed is `cfg.seed + 100000*dgp_id + r` and nothing carries over between
replications. `tests/test_chunk_merge.m` asserts this, and
`merge_montecarlo` refuses overlapping chunks and marks incomplete
merges.

### Reading the results

```matlab
report_montecarlo('results/simulation/mc_headline_sparse.mat')
```

prints both integrated-RMSE conventions, the bias–variance
decomposition, interval performance under both band constructions, and
the `tau` diagnostics. The same numbers are in the CSVs:

| file | contents |
|---|---|
| `*_metrics.csv` | bias, variance, MSE, RMSE, coverage, length per (estimator, response, horizon) |
| `*_integrated.csv` | both IRMSE definitions, early/late splits, integrated bias and variance |
| `*_tau.csv` | posterior mean/median and `P(tau > 1)` per (equation, block, horizon) |
| `*_tau_localization.csv` | localisation, false-positive and chain diagnostics |
| `*_meta.txt` | DGP parameters, seeds, `R/T/p/H`, estimator settings, sampler sizes |

---

## 4. Sensitivity of the methodological approximations

```matlab
run_sensitivity_approximations(struct('mc_file', ...
    fullfile('results', 'simulation', 'mc_headline_sparse.mat')))
```

Writes `results/simulation/sensitivity_approximations.csv`; the interpretation is
in `docs/APPROXIMATIONS.md`.

---

## 5. Reproducing the earlier (pre-revision) numbers

The `legacy` preset restores the four-estimator stack, the FMAR `h = 1`
convention and the draw-averaged posterior mean:

```matlab
cfg = mc_preset('legacy');
mc  = run_montecarlo(cfg, 'sparse');
```

Two caveats. First, the sampler was made about three times faster
(hoisted loop invariants, vectorised block sums, and replacing Octave's
`assert` inside `draw_gamma`, which cost ~100 µs per call). The order in
which random numbers are consumed is unchanged, so results agree to
floating-point associativity — about `2e-16` on the IRFs — but not
bit-for-bit. Second, `cfg.blp.point_estimate` must be `'draw_mean'`
(the `legacy` preset sets it) to reproduce the old Monte Carlo noise in
the adaptive point estimates.

---

## 6. Chapter 7 empirical pipeline (current design: external instrument)

Full instructions, data sources and the design are in `empirical/README.md`,
`empirical/docs/DATA.md` and `empirical/docs/DESIGN.md`. The outcome and
indicator series are third-party data and are **not** in the repository
(they sit gitignored on the author's machine; the EA-EMPD event extracts and
the derived monthly instrument panel are shipped).

```matlab
addpath(genpath('<repo>'))
build_instrument_series();                       % EA-EMPD workbook extract -> 51 monthly instruments (seconds)
fetch_outcome_data();  ea_fetch_v2_series();     % Eurostat / ECB series (download + validation)
import_ois_daily();                              % daily 1M OIS export -> monthly eom / avg (see DATA.md)
jk = struct('path', fullfile(ea_paths().raw, 'shocks_ecb_mpd_me_m.csv'), 'column', 'MP_pm', 'name', 'jk_mp_pm', 'scale', 100);
assemble_dataset_v2(struct('system', 'ois4', 'external', {{jk}}));   % -> empirical/data/ea_dataset_v2_ois4.mat
QUICK = true; SYSTEM = 'ois4'; run(fullfile(ea_paths().empirical, 'RUN_EMPIRICAL_IV.m'))   % interface check, ~4 min, meaningless numbers
SYSTEM = 'ois4'; run(fullfile(ea_paths().empirical, 'RUN_EMPIRICAL_IV.m'))                 % steps REL A B C D E, ~30 min
```

The null calibrations, one per design (`p`) and per residual scheme; about
6 s per replication at `p = 12` for `K = 4` in MATLAB (3 s at `p = 2`):

```matlab
for p = [2 4 6 12]
    run_null_calibration(struct('p', p, 'null', struct('n_rep', 500, 'stem', sprintf('iv_ois4_p%d', p))));                        % iid residual resampling
    run_null_calibration(struct('p', p, 'null', struct('n_rep', 500, 'sim_method', 'wild', 'stem', sprintf('iv_ois4_p%d', p))));  % wild bootstrap (-> *_wild)
end
SYSTEM = 'ois4'; DO_REL = false; DO_A = false; DO_C = false; DO_E = false; run(fullfile(ea_paths().empirical, 'RUN_EMPIRICAL_IV.m'))   % protocol + cross-p against every null found
```

`run_null_calibration` defaults to the `ois4` dataset (`ea_paths().dataset`);
pass `'dataset'` for another system. The headline iid null at `p = 12` was run
with `n_rep = 1000` so that Holm's step-down has the resolution it needs over
16 cells (`R >= 319`); see `empirical/docs/RESULTS.md`.

Outputs land in `results/empirical/`: `iv_<system>_p<p>*.{mat,csv}`,
`null_calibration_iv_*.mat` (+ `.key.txt`, `null_thresholds_iv_*.csv`),
`ablation_iv_*`, `iv_<system>_cross_p.csv`, figures in `figures/`.

Without the data the package is still testable: `test_empirical_iv` and
`test_empirical_pipeline` exercise every interface on synthetic fixtures that
every driver refuses (`ds.synthetic = true`). **No empirical result may be
produced from a fixture.**

## 7. Legacy v1 pipeline (internal instrument, 2026-09-12)

Kept runnable in `empirical/legacy_v1/`; writes to
`results/empirical/legacy_v1/`. Design: `empirical/docs/DESIGN_V1_LEGACY.md`;
results: the first part of `empirical/docs/RESULTS_LOG.md`.

```matlab
build_shock_series(); fetch_outcome_data(); assemble_dataset();   % -> empirical/data/ea_dataset.mat
SMOKE_TEST_EMPIRICAL_V1                                            % interface check (2 s)
RUN_EMPIRICAL_V1                                                   % steps A-F, ~4 min at p = 12, H = 48
run_null_calibration(struct('null', struct('n_rep', 500, 'dataset', ea_paths().dataset_legacy)));   % ~85 min at K = 5
DO_A = false; DO_B = false; DO_D = false; RUN_EMPIRICAL_V1         % protocol tables
```
