# Monday data checklist (v2 pipeline)

What to obtain, where to put it, what the file must contain, and the exact
commands. Nothing below is downloaded automatically except through the two
fetch helpers named here, and only when you run them.

## 1. Required: daily 1-month OIS level (Refinitiv/LSEG RIC `EUREON1M=`)

| item | value |
|---|---|
| file | `empirical/data/raw/ois1m_ea_daily.csv` (gitignored, like the other raw files) |
| frequency | daily (business days) |
| coverage needed | 1999-12-01 to 2019-12-31 at least (earlier and later is fine) |
| columns | one date column (header containing `Date`, `Time` or `Timestamp`) and one value column (header containing `Close`, `Last`, `Mid`, `Value`, `Price` or `EUREON`); percent per annum; other columns are ignored |
| date format | `yyyy-mm-dd`, `yyyy/mm/dd`, `dd/mm/yyyy`, `dd.mm.yyyy` or `dd-mmm-yyyy`; an appended time is ignored; **not** Excel serial numbers |
| field | mid or close of the fixed rate; if bid/ask only, export the mid |
| sidecar | after placing it, run the provenance line below so the dataset records the RIC and the export date |

What the importer does (`import_ois_daily`): validates coverage (no month
with fewer than 10 quotes, no gap over 10 days, values in [−2, 8]), builds
`derived/ois1m_ea_monthly.{csv,mat}` with `eom` (last quote of the month:
the END-OF-MONTH indicator, paired with `sum`-aggregated surprises), `avg`
(monthly average: paired with `kilian`-aggregated surprises), `first`,
`chg_eom`, `n_days`. The EONIA re-basing of 2 October 2019 (€STR + 8.5 bp)
is a level convention inside the last quarter of the sample and needs no
adjustment; EONIA ended on 3 January 2022, outside the sample.

Optional but useful on the same terminal: `EUREON1Y=` (1-year OIS, daily),
same layout, saved as `empirical/data/raw/ois1y_ea_daily.csv`. It would
give the legacy 1-year instrument a matched, credit-premium-free indicator
(the registry entry can be added in five lines; not wired yet).

## 2. Required for `ois4` / `ois6`: seasonally adjusted HICP (ECB Data Portal)

Run once (downloads and validates, writes provenance sidecars):

```matlab
addpath(genpath('/Users/stefanograziosi/Documents/GitHub/block_adaptive_blp')); ea_fetch_v2_series()
```

It fetches, with keys verified on 2026-09-13:

| file | series key | needed for |
|---|---|---|
| `hicp_sa_ea.csv` | `ICP.M.U2.Y.000000.3.INX` (HICP overall, ECB, working-day and seasonally adjusted, 2015 = 100) | `ois4`, `ois6` |
| `loans_nfc_ea.csv` | `BSI.M.U2.Y.U.A20T.A.I.U2.2240.Z01.E` (adjusted loans to NFCs, index of notional stocks, SA; from 2003-01) | `ois6` |
| `ccb_nfc_ea.csv` | `MIR.M.U2.B.A2I.AM.R.A.2240.EUR.N` (composite cost of borrowing, NFCs; from 2003-01) | `ois6` (`lend_spread` = ccb − indicator) |
| `eonia_ea.csv` | `FM.M.U2.EUR.4F.MM.EONIA.HSTA` (EONIA monthly average) | optional indicator `eonia_avg` |
| `euribor1m_ea.csv` | `FM.M.U2.EUR.RT.MM.EURIBOR1MD_.HSTA` (1M Euribor monthly average) | optional fallback indicator (credit premium) |

`ea_fetch_v2_series(struct('dry_run', true))` prints the URLs without
downloading, for a manual export (save as csvdata / SDMX-CSV under the file
name in the table).

## 3. Optional: an externally supplied information-robust shock

* Jarociński's updated ECB shocks (CC BY 4.0): download
  `shocks_ecb_mpd_me_m.csv` from `github.com/marekjarocinski/jkshocks_update_ecb`
  to `empirical/data/raw/` and pass it as an external instrument:
  ```matlab
  assemble_dataset_v2(struct('system', 'ois4', 'external', {{struct('path', ...
      '/Users/stefanograziosi/Documents/GitHub/block_adaptive_blp/empirical/data/raw/shocks_ecb_mpd_me_m.csv', ...
      'column', 'MP_pm', 'name', 'jk_mp_pm', 'scale', 100)}}))
  ```
  (`scale` = 100 if the file is in percentage points; check its README.
  The series is scaled to the 1-year-OIS surprise sd, so treat it as a
  1-year-type instrument: it is an alternative, not the headline for a 1M
  indicator.)
* Ricco–Savini–Tuteja informationally robust factors: not published with
  CEPR DP 19679; request from the authors. Any csv with `year, month,
  value` (or a `date` column) loads through the same `external` option.

## 4. Commands, in order (after the files are in place)

```matlab
addpath(genpath('/Users/stefanograziosi/Documents/GitHub/block_adaptive_blp'))
ea_write_provenance(fullfile(ea_paths().raw, 'ois1m_ea_daily.csv'), 'Refinitiv EUREON1M= (1M EONIA swap, daily mid)', 'LSEG Workspace export, <date>', true)
import_ois_daily();                                  % monthly eom / avg / first
assemble_dataset_v2(struct('system', 'ois4'));       % 2000m1-2019m12
assemble_dataset_v2(struct('system', 'ois6'));       % 2003m1-2019m12 (loans, spread)
SYSTEM = 'ois4'; run(fullfile(ea_paths().empirical, 'RUN_EMPIRICAL_IV.m'))
```

Then the nulls for the new design (about 2.5 s per replication at p = 2,
~5 s at p = 12 for K = 4; run them while reading the relevance table):

```matlab
for p = [2 4 6 12]
    run_null_calibration(struct('p', p, 'null', struct('n_rep', 200, ...
        'dataset', fullfile(ea_paths().data, 'ea_dataset_v2_ois4.mat'), 'stem', sprintf('iv_ois4_p%d', p))));
end
SYSTEM = 'ois4'; DO_REL = false; DO_A = false; DO_C = false; DO_E = false; run(fullfile(ea_paths().empirical, 'RUN_EMPIRICAL_IV.m'))
```

## 5. What to look at first

1. `results/iv_ois4_p12_relevance.csv`: the row `z_gcs_1m_adj_sum`
   (headline for the end-of-month indicator) and `z_gc_1m_adj_sum`; `F_eff`
   against 12.05 / 15.06 / 23.11 / 37.42 (Montiel Olea–Pflueger) and 10;
   `lead1_t`, `lead2_t` (should be insignificant); `F_excl_crisis`;
   `top_month`. AGKL's own numbers for the 1M OIS on 1M-OIS surprises are
   F = 4.0 (meetings), 14.1 (speeches), 12.6 (all events).
2. Whether `sum` beats `kilian` on the `eom` indicator and the reverse on
   `avg` (`assemble_dataset_v2(struct('system','ois4','vars',{{'ois1m_avg','ip','hicp_sa','stoxx'}}))`
   gives the average-indicator dataset). This is the Kilian test.
3. Only if `F_eff` is above the 10 %–bias threshold: the IRFs in
   `results/iv_ois4_p12_irf.csv` and the LP-IV Anderson–Rubin sets in
   `results/iv_ois4_p12_lpiv.csv`. Below it, report the AR sets only.
