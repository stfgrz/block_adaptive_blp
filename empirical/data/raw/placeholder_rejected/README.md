# Quarantined placeholder files — NOT DATA

The four csv files in this folder are **not observed data**. They were
generated while the Ch. 7 pipeline was being wired up, and they sat in
`empirical/data/raw/` under the real filenames (`ip_ea.csv`,
`hicp_ea.csv`, `rate1y_ea.csv`, `stoxx50_ea.csv`) with the right columns,
the right monthly frequency and no gaps — so `fetch_outcome_data` parsed
them, reported *"all four outcome files ready"*, and the entire chapter
would have been estimated on fabricated numbers.

They are kept here rather than deleted so the diagnosis is checkable.
**Do not move them back.**

## How we know they are synthetic

| file | fingerprint |
|---|---|
| `hicp_ea.csv` | rises by **exactly +0.13 index points every month** for 323 consecutive months (roughness statistic 0.00) |
| `rate1y_ea.csv` | sits at **exactly −0.5000** for 87 consecutive months; 32% of its month-to-month changes are exactly identical; 12M Euribor was near 4% in 2023, not −0.5% |
| `ip_ea.csv` | **no COVID collapse** — April 2020 reads 98.18 against 98.09 in February, where euro-area industrial production actually fell about 28% |
| `stoxx50_ea.csv` | **no 2008 and no March-2020 crash**; a smooth monotone path throughout (roughness 0.09) |

The general tripwire is the *roughness ratio* `sd(diff²y)/sd(diff y)`,
which is ≈1.4 for any real monthly macro series (exactly √2 for a random
walk) and 0.00–0.17 for all four of these. `empirical/data/ea_check_series.m`
implements it and two related fingerprints, and both `fetch_outcome_data`
and `assemble_dataset` now refuse a series that trips them.

## What to do instead

Run `fetch_outcome_data()`. It will try the portal URLs and, if they
fail, print the exact manual download route for each series (portal,
dataset, filters, target filename). Once all four validate, it prints
*"all four outcome files ready and validated."*
