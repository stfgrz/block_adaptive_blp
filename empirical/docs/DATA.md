# Data: sources, acquisition, validation, known issues

Every series the chapter uses, where it comes from, how it is obtained and
checked, and what is known to be imperfect about it. Raw files live in
`empirical/data/raw/` and are **gitignored** as third-party data, except the
EA-EMPD event extracts and the provenance sidecars; derived monthly panels
built from shipped inputs are tracked.

## 1. Series

| code | series | source and key | file | tracked |
|---|---|---|---|---|
| `ois1m` (`_eom`, `_avg`, `_first`) | 1-month EONIA swap (OIS) rate, business daily, bid and ask; mid taken by `import_ois_daily` | ECB Data Portal copy of the Refinitiv series `FM.B.U2.EUR.RT.SI.EUREON1M_` (.BID / .ASK), exported by the user 2026-09-14 | `raw/ois1m_ea_daily.csv` (Date, Bid, Ask); monthly `derived/ois1m_ea_monthly.{csv,mat}` | raw no; derived yes |
| `ip` | industrial production, EA20, industry excl. construction (B-D), seasonally and calendar adjusted, 2021 = 100 | Eurostat `sts_inpr_m`, `M.PRD.B-D.SCA.I21.EA20` (the 2026-09-01 Data Browser export covered 40 geographies × 3 aggregates and was reduced with `ea_extract_series`; the full export is kept locally) | `raw/ip_ea.csv` | no |
| `hicp` | HICP all items, NSA, 2015 = 100 | Eurostat `prc_hicp_midx`, `M.I15.CP00.EA` | `raw/hicp_ea.csv` | no |
| `hicp_sa` | HICP all items, working-day and seasonally adjusted, 2015 = 100 | ECB `ICP.M.U2.Y.000000.3.INX` | `raw/hicp_sa_ea.csv` | no |
| `hicp_yoy` | 100·(log P_t − log P_{t−12}) from the NSA index | derived in `assemble_dataset_v2` | — | — |
| `stoxx` | EURO STOXX 50, **monthly average** of daily closes | ECB `FM.M.U2.EUR.DS.EI.DJES50I.HSTA` | `raw/stoxx50_ea.csv` | no |
| `i1y` | 12-month Euribor, monthly average (legacy indicator) | ECB `FM.M.U2.EUR.RT.MM.EURIBOR1YD_.HSTA` | `raw/rate1y_ea.csv` | no |
| `euribor1m`, `eonia_avg` | 1-month Euribor and EONIA, monthly averages (fallback indicators) | ECB `FM.M.U2.EUR.RT.MM.EURIBOR1MD_.HSTA`, `FM.M.U2.EUR.4F.MM.EONIA.HSTA` | `raw/euribor1m_ea.csv`, `raw/eonia_ea.csv` | no |
| `loans` | adjusted loans to NFCs, index of notional stocks, SA (from 2003m1) | ECB `BSI.M.U2.Y.U.A20T.A.I.U2.2240.Z01.E` | `raw/loans_nfc_ea.csv` | no |
| `ccb`, `lend_spread` | composite cost of borrowing, NFCs (from 2003m1); spread = ccb − policy indicator | ECB `MIR.M.U2.B.A2I.AM.R.A.2240.EUR.N` | `raw/ccb_nfc_ea.csv` | no |
| instruments `z_*` | EA-EMPD event-window asset-price changes (OIS 1M/3M/1Y, STOXX50E) for all Governing Council meetings and Executive Board / President speeches since 1999 | Altavilla, Gürkaynak, Kind & Laeven (2025), ECB WP 3157, workbook `EA-EMPD.en.xlsx`; all 54 columns extracted 2026-09-13 | `raw/ea_empd_events_full.csv` (4833 events); legacy 12-column extract `raw/ea_empd_events.csv` | **yes** |
| `z_ext_jk_*` | Jarociński & Karadi updated ECB monetary-policy and information shocks (monthly and daily), CC BY 4.0 | `github.com/marekjarocinski/jkshocks_update_ecb`, `shocks_ecb_mpd_me_{m,d}.csv`, downloaded 2026-09-14 | `raw/shocks_ecb_mpd_me_{m,d}.csv` | no |

All keys were verified on the ECB Data Portal API on 2026-09-13; the
sidecars `raw/<file>.source.txt` (tracked) record series, URL and date for
every accepted file.

## 2. Acquisition

```matlab
addpath(genpath('<repo>'))
fetch_outcome_data();        % ip, hicp (NSA), stoxx, rate1y: downloads, validates, writes sidecars;
                             % prints the manual route for a series whose download fails
ea_fetch_v2_series();        % hicp_sa, loans, ccb, eonia, euribor1m (ECB portal, csvdata format)
ea_fetch_v2_series(struct('dry_run', true))   % just print the URLs
```

**1-month OIS (manual).** The daily series is not on the free portal API as
a single field. Export it from the ECB Data Portal (`FM.B.U2.EUR.RT.SI.EUREON1M_`,
fields BID and ASK; or from a Refinitiv/LSEG terminal as `EUREON1M=`) with
one date column and either a mid/close column or Bid and Ask columns, dates
as text (`yyyy-mm-dd` and the common European formats are recognised, Excel
serial numbers are not), coverage at least 1999-12-01 to 2019-12-31, save it
as `raw/ois1m_ea_daily.csv`, record its provenance and import it:

```matlab
ea_write_provenance(fullfile(ea_paths().raw, 'ois1m_ea_daily.csv'), '<series description>', '<where and when exported>', true)
import_ois_daily();          % checks coverage (no month < 10 quotes, no gap > 10 days, values in [-2, 8]) and
                             % writes derived/ois1m_ea_monthly.{csv,mat}: eom, avg, first, chg_eom, n_days
```

**Instruments.** `build_instrument_series()` reads the shipped EA-EMPD
extract and writes `derived/instruments_monthly.{csv,mat}` (51 monthly
series in basis points, zero in months without a contributing event) and
`derived/instruments_events.csv` (the event-level table). Every instrument is
also carried in a `_bs` version orthogonalised on pre-event public
information (`ea_pre_event_info`). `assemble_dataset_v2` attaches an
external series through `opts.external` (`import_external_instrument`; for
Jarociński's file `scale = 100` converts percentage points to basis points).

**Assembly.** `assemble_dataset_v2(struct('system', 'ois4'))` (or `'vars'`
explicitly) validates every input on the estimation window, builds `ds.Y`
(levels in percent, indices as 100·log), `ds.Z` (instruments), `ds.pred`
(pre-event information), `ds.isrw` (random-walk prior centre flags) and
`ds.meta.sources`, and saves `data/ea_dataset_v2_<system>.mat`. A synthetic
OIS fixture makes `ds.synthetic = true`, which every driver refuses.

## 3. Validation

`ea_check_series(name, ym, val, need0, need1)` runs on every raw series and
the builders **stop** when it fails:

* coverage of the required window, no duplicated `TIME_PERIOD` (an
  over-broad Eurostat export with several series stacked in one file is
  caught here; reduce it with `ea_extract_series`);
* the roughness ratio `sd(diff²y)/sd(diff y)`, ≈ 1.4 for any real monthly
  macro series (exactly √2 for a random walk); 0.00–0.18 flags a
  generated series (the quarantined placeholders in
  `raw/placeholder_rejected/` are the cautionary tale);
* the modal-step share and the longest flat run (clipped or
  constant-increment series).

On the real files the ratio is 1.38 (HICP), 1.44 (IP), 1.24 (STOXX), 0.78
(12M Euribor, a smooth monthly average). `opts.skip_plausibility = true`
overrides a fingerprint only after looking at the series.

Reference statistics of the shipped event file (window 2001m1–2019m12,
printed by `build_instrument_series` and pinned by `test_empirical_iv`):
221 GC monetary-event windows (216 with a 1-month quote), 1608 speeches
passing the AGKL filters (1021 with an adjusted 1-month surprise),
`z_gc_1y_sum` sd 4.29 bp with 24 zero months (identical to the legacy
`mps_gc_1y` to 5e-7 bp), `z_gcs_1m_adj_sum` sd 4.28 bp with 13 zero months.

## 4. Known issues and choices

* **The STOXX index is a monthly average.** A surprise on business day `d`
  of a month with `D` business days moves that month's average by
  `(D−d+1)/D` of its effect and the next month's by the rest (Kilian
  2024). The impact response of `stoxx` is therefore smeared over `h = 0`
  and `h = 1`; an end-of-month index level would fix it but the ECB key
  (`...HSTE`) does not exist. Daily STOXX 50 closes from another source
  (e.g. STOXX's own historical data, Stooq, Yahoo Finance) would allow an
  end-of-month or event-aligned series; not fetched.
* **The end-of-month 1M OIS pairs with the within-month sum** of the
  surprises; a monthly-average indicator would pair with the Kilian
  day-weighted aggregation (`kilian` / `kilianfix`). Both aggregations are
  always built so the relevance table can show which the data prefer,
  but the headline pairing is fixed a priori.
* **Sample period.** 2000m1–2019m12 (`T = 240`; with `p = 12` the
  effective sample is 2001–2019). The short rate is essentially pinned
  from 2012 (innovation standard deviation 0.04 pp against 0.13–0.16 in
  2001–2011), which is why the calibration offers a wild-bootstrap null.
  Post-2019 data would need the COVID months handled (Lenza–Primiceri
  volatility rescaling or pandemic dummies; not implemented).
* **1M OIS quotes are missing for 5 of the 221 meetings** in the EA-EMPD
  file (early years); the 1-year quote is complete.
* **EONIA re-basing** (€STR + 8.5 bp from 2 October 2019) is a level
  convention inside the last quarter of the sample; EONIA was discontinued
  on 3 January 2022, outside the sample.
* **12-month Euribor** embeds a bank credit premium, material 2007–2012;
  it is the legacy indicator only.
* **AGKL's exclusion of speech windows containing a high-relevance data
  release** cannot be replicated (no release calendar); speeches outside
  regular trading hours or on non-trading days are excluded, and events
  whose day-count factor exceeds 6 (or whose next meeting is 30 or more
  days away) are dropped from the adjusted 1-month surprise.
* **The instrument is predictable from the VAR's own lags** (event
  months: R² 0.27 on 48 lags, HAC Wald p ≈ 0.001; the one-regressor test
  on the BVAR-predicted rate change is insignificant). The proxy
  identification conditions on those lags, so this is admissible under
  Stock–Watson lead–lag exogeneity, but it is why the raw regression of the
  monthly rate change on the surprise (AGKL's Figure 9) shows nothing here
  (F 0.7) while the innovation regression gives 13.7.

## 5. Legacy v1 inputs

`legacy_v1/build_shock_series` builds the six monthly surprise series of the
v1 design (`mps_gc_1y` baseline, `mps_gc_3m`, JK split, speeches) from the
12-column extract `raw/ea_empd_events.csv` into `derived/shocks_monthly.{csv,mat}`
(tracked; `test_empirical_iv` checks that the new builder reproduces it).
`legacy_v1/assemble_dataset` builds the five-variable `data/ea_dataset.mat`.
