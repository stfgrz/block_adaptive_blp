# Thesis pitch deck

Two versions live here:

* `thesis_pitch_2026-09-12.tex` (eleven slides) records the state as of
  12 September 2026: the fair-benchmark headline runs with five estimators,
  the diagnostic framing, and the first pass at the euro-area application
  (the three design corrections the real data forced, the τ map, the
  calibrated reading). Numbers come from `../results/mc_headline_*.mat`
  and `../results/*_p12*.csv`.
* `archive/thesis_pitch_2026-09-08.tex` (nine slides) is the 8 September
  2026 version, kept for the record. It predates the pooled estimator, the
  fair benchmark and all empirical results; its roadmap items are done.

Both decks predate the external-instrument redesign of the empirical chapter
(13–15 September 2026); their euro-area slides describe the v1 design.

Build from this directory with TeX Live (`build/` is gitignored):

```sh
mkdir -p build
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build thesis_pitch_2026-09-12.tex
cp build/thesis_pitch_2026-09-12.pdf thesis_pitch_2026-09-12.pdf
```

The source is intentionally self-contained. The evidence chart uses `pgfplots`
with values extracted from `../results/simulation/mc_fmar_final_sparse.mat`, so the final
PDF does not depend on generated image assets.
