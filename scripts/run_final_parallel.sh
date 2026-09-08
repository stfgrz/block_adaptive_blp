#!/usr/bin/env bash
# Run a FINAL (R = 500) Monte Carlo for one DGP across NCHUNK Octave
# processes, then merge the chunks and write the summary + csv exports.
#
# Replication r's seed is a pure function of r, so the chunks reassemble
# into exactly what a serial run would have produced
# (montecarlo/merge_montecarlo.m asserts the design matches and that the
# replication indices cover 1..R exactly once).
#
# OUTPUT NAMING.  The merged result is written to
#   results/mc_headline_<dgp>.mat
# NOT results/mc_<preset>_<dgp>.mat.  The repository already contains
# mc_final_<dgp>.mat files from earlier prototype-mode runs, and those
# are baselines that must be preserved; the "headline" prefix keeps the
# new five-estimator runs from overwriting them.  The per-process chunks
# still carry the preset name so a chunk directory is self-describing.
#
# Usage:  scripts/run_final_parallel.sh DGP [R] [NCHUNK] [PRESET]
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DGP="${1:?usage: run_final_parallel.sh DGP [R] [NCHUNK] [PRESET]}"
R="${2:-500}"
NCHUNK="${3:-4}"
PRESET="${4:-final}"
CH="$ROOT/results/chunks"
LOGDIR="$CH/logs"
mkdir -p "$CH" "$LOGDIR"

echo "final MC: dgp=$DGP R=$R preset=$PRESET over $NCHUNK chunks"
PER=$(( (R + NCHUNK - 1) / NCHUNK ))
FILES=""
for ((c = 0; c < NCHUNK; c++)); do
  R0=$(( c * PER + 1 ))
  R1=$(( (c + 1) * PER ))
  if (( R1 > R )); then R1=$R; fi
  if (( R0 > R )); then continue; fi
  F="$CH/${PRESET}_${DGP}_${R0}_${R1}.mat"
  FILES="$FILES '$F',"
  ( octave-cli --no-gui --quiet --eval \
      "addpath(genpath('$ROOT')); ov = struct('mc', struct('n_rep', $R)); \
       run_mc_chunk('$PRESET','$DGP',$R0,$R1,'$F',ov);" \
      > "$LOGDIR/${PRESET}_${DGP}_${R0}_${R1}.log" 2>&1 ) &
done
wait

# NOTE for anyone editing the --eval block below: bash removes the
# backslash-newlines, so the whole thing reaches Octave as ONE line.  A
# '%' comment inside it would comment out everything after it.  Keep the
# code comment-free and put explanations in named functions instead
# (that is why the compaction step is montecarlo/compact_mc.m).
octave-cli --no-gui --quiet --eval \
  "addpath(genpath('$ROOT')); \
   fl = {${FILES%,}}; ch = cell(1, numel(fl)); \
   for k = 1:numel(fl), L = load(fl{k}); ch{k} = L.mc; end; \
   mc = merge_montecarlo(ch); \
   s = summarize_montecarlo(mc); \
   mc = compact_mc(mc); \
   stem = sprintf('mc_headline_%s', '$DGP'); \
   out = fullfile('$ROOT', 'results', [stem '.mat']); \
   save(out, 'mc', 's', '-v7'); \
   export_montecarlo_csv(s, fullfile('$ROOT', 'results', stem)); \
   fprintf('merged %d replications -> %s\n', s.R, out);" \
  2>&1 | tee "$LOGDIR/merge_headline_${DGP}.log"
