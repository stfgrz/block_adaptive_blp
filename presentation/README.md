# Thesis pitch deck

Two versions live here:

* `thesis_pitch_2026-09-12.tex` (eleven slides) records the state as of
  12 September 2026: the fair-benchmark headline runs with five estimators,
  the diagnostic framing, and the first pass at the euro-area application
  (the three design corrections the real data forced, the τ map, the
  calibrated reading). Numbers come from `../results/mc_headline_*.mat`
  and `../results/*_p12*.csv`.
* `thesis_pitch.tex` (nine slides) is the 8 September 2026 version, kept
  for the record. It predates the pooled estimator, the fair benchmark and
  all empirical results; its roadmap items are done.

Build either from this directory with TeX Live (replace the file name):

Build from this directory with TeX Live:

```sh
mkdir -p build
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build thesis_pitch.tex
cp build/thesis_pitch.pdf thesis_pitch.pdf
```

The source is intentionally self-contained. The evidence chart uses `pgfplots`
with values extracted from `../results/mc_fmar_final_sparse.mat`, so the final
PDF does not depend on generated image assets.
