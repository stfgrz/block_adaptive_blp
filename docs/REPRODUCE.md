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

About 25 minutes. `test_fmar_port` is skipped unless `FMAR_PATH` points
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
collect_grid_results                    % -> results/grid/grid_summary.csv
```

Outputs land in `results/grid/`: one `.mat` and one set of CSVs per
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

Each call writes `results/mc_headline_<dgp>.mat` plus the CSV exports.
Roughly 3 hours per `R = 500` DGP on 4 cores.

Serial equivalent (identical results, ~11 hours per DGP):

```matlab
cfg = mc_preset('final');
mc  = run_montecarlo(cfg, 'sparse');
s   = summarize_montecarlo(mc);
save(fullfile('results', 'mc_headline_sparse.mat'), 'mc', 's', '-v7');
```

The chunked and serial runs are **bit-identical**: replication `r`'s
seed is `cfg.seed + 100000*dgp_id + r` and nothing carries over between
replications. `tests/test_chunk_merge.m` asserts this, and
`merge_montecarlo` refuses overlapping chunks and marks incomplete
merges.

### Reading the results

```matlab
report_montecarlo('results/mc_headline_sparse.mat')
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
    fullfile('results', 'mc_headline_sparse.mat')))
```

Writes `results/sensitivity_approximations.csv`; the interpretation is
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

## 6. Chapter 7 empirical pipeline

The four outcome series are licensed and are **not** in the repository.

```matlab
build_shock_series();      % EA-EMPD events -> monthly surprises (ships with the repo)
fetch_outcome_data();      % downloads and VALIDATES the four series
assemble_dataset();        % -> empirical/data/ea_dataset.mat
SMOKE_TEST_EMPIRICAL       % interface check on the real dataset
RUN_EMPIRICAL              % the multi-hour job
```

Without the licensed data, the package is still testable:

```matlab
test_empirical_pipeline    % parses, builds shocks, runs all five
                           % estimators on a synthetic fixture
```

The fixture is simulated data stamped `ds.synthetic = true`, and both
`RUN_EMPIRICAL` and `SMOKE_TEST_EMPIRICAL` refuse a dataset carrying
that stamp. **No empirical result in this project may be produced from
it.** See `empirical/README_EMPIRICAL.md`.
