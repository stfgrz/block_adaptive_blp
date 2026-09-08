#!/usr/bin/env bash
# Run the remaining headline Monte Carlo experiments, in priority order.
#
# Launched with setsid so the workers live in their own process group and
# survive an interruption of the shell that started them (an earlier run
# of these experiments was killed when the launching session paused).
set -u
cd /home/user/block_adaptive_blp
echo "started at $(date -u +%H:%M)"
./scripts/run_final_parallel.sh correct      500 4
echo "correct done at $(date -u +%H:%M)"
./scripts/run_final_parallel.sh intermediate 250 4
echo "intermediate done at $(date -u +%H:%M)"
./scripts/run_final_parallel.sh dense        250 4
echo "dense done at $(date -u +%H:%M)"
echo ALL_FINALS_DONE
