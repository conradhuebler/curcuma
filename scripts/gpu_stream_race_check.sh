#!/bin/bash
# Claude Generated (Oct 2026). Run-to-run check of the GPU GFN-FF MD for stream-ordering races.
# Usage: [GPU=rocm] scripts/gpu_stream_race_check.sh BINARY XYZ RUNS [SOLVER ...]   (default solvers: cholesky ldlt lu pcg ppcg batched)
# Runs the MD RUNS times per EEQ solver and prints how many distinct trajectories came out (physics columns only:
# the wall-clock columns differ between runs and are left out). One line per solver; "RUNS <hash>" = all identical.
# A race shows up much more often with a second GPU job on the same card, e.g.
#   while true; do curcuma -sp big.xyz -method gfn2 -gpu ${GPU:-cuda} -gpu_devices 0 -no_bmt -verbosity 0; done
# Note: `batched` drops the inter-fragment coupling, its numbers differ from the other solvers (by design).
B=$1; X=$2; N=$3; shift 3
SOLVERS=${@:-cholesky ldlt lu pcg ppcg batched}
D=$(mktemp -d); cp "$X" $D/s.xyz; cd $D || exit 1
for m in $SOLVERS; do
  for i in $(seq 1 $N); do
    rm -f s.topo.json
    $B -md s.xyz -method gfnff -gpu ${GPU:-cuda} -gpu_devices 0 -gpu_strict true -no_bmt -T 300 -dt 1 -seed 42 -max_time 60 \
       -print 10 -dump 100000000 -verbosity 1 -eeq_solver.solve_method $m 2>&1 | sed 's/\x1b\[[0-9;]*m//g' \
       | grep -E "^\s*[0-9]+\.[0-9]+\s+-?[0-9]+\.[0-9]+\s" | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | md5sum | cut -c1-8
  done | sort | uniq -c | sort -rn | tr '\n' ' '
  echo " <- $m"
done
rm -rf $D
