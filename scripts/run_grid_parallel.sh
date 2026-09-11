#!/usr/bin/env bash
# Run the exploratory experiment grid across NPROC Octave processes.
#
# Each cell of montecarlo/experiment_grid_spec.m is independent (its own
# seeds, its own output file), so the cells are simply dealt out to the
# workers round-robin.  Nothing is shared and nothing is appended to, so
# the result does not depend on how many workers are used.
#
# Usage:  scripts/run_grid_parallel.sh [NPROC] [LOGDIR]
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
NPROC="${1:-4}"
LOGDIR="${2:-$ROOT/results/grid/logs}"
mkdir -p "$LOGDIR"

NCELL=$(octave-cli --no-gui --quiet --eval \
  "addpath(genpath('$ROOT')); c = experiment_grid_spec(); printf('%d\n', numel(c));" \
  2>/dev/null | tail -1)
echo "experiment grid: $NCELL cells over $NPROC workers; logs in $LOGDIR"

for ((w = 0; w < NPROC; w++)); do
  IDX=""
  for ((k = w + 1; k <= NCELL; k += NPROC)); do IDX="$IDX $k"; done
  # A worker with no assigned cells must NOT be launched: passing [] to
  # run_experiment_grid means "every cell", so an extra worker would
  # silently re-run the whole grid and overwrite the others' output.
  if [ -z "$IDX" ]; then
    echo "worker $w: no cells assigned ($NPROC workers for $NCELL cells); skipping"
    continue
  fi
  IDX="[$(echo "$IDX" | tr ' ' ',')]"
  echo "worker $w: cells $IDX"
  ( octave-cli --no-gui --quiet --eval \
      "addpath(genpath('$ROOT')); run_experiment_grid($IDX);" \
      > "$LOGDIR/worker_$w.log" 2>&1 ) &
done
wait
echo "grid finished; collecting"
octave-cli --no-gui --quiet --eval \
  "addpath(genpath('$ROOT')); collect_grid_results();" \
  2>&1 | tee "$LOGDIR/collect.log"
