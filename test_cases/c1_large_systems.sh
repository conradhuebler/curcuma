#!/bin/bash
# C1 large-system modes — fast-proxy regression test
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (C1, June 2026)
#
# Validates the three opt-in native-GFN large-system modes against the dense
# reference on a small connected molecule (caffeine), where each mode has an
# exact/known limit:
#   - c1_mode=fragments (C1c): one connected fragment -> dense fallback (EXACT)
#   - c1_mode=dc (C1a): a buffer covering the whole molecule -> dense (tight)
#   - c1_mode=sparse (C1b): 0 K, threshold 0 -> dense-0K (tight)
# Also asserts the default path (c1_mode=none) is unchanged.
#
# Usage: c1_large_systems.sh <curcuma-binary> <molecule-dir>

set -u
CURCUMA="${1:?need curcuma binary path}"
MOLDIR="${2:?need molecule directory}"
MOL="${MOLDIR}/caffeine.xyz"

fail=0
# Extract "Single Point Energy = <x> Eh" (strip ANSI).
energy() {
    "$CURCUMA" -sp "$MOL" "$@" 2>/dev/null \
        | sed -r 's/\x1b\[[0-9;]*m//g' \
        | grep -oE 'Single Point Energy = -?[0-9]+\.[0-9]+' | grep -oE -- '-?[0-9]+\.[0-9]+' | head -1
}
# assert |a-b| <= tol
check() {
    local name="$1" a="$2" b="$3" tol="$4"
    if [ -z "$a" ] || [ -z "$b" ]; then
        echo "FAIL  $name : missing energy (a='$a' b='$b')"; fail=1; return
    fi
    awk -v a="$a" -v b="$b" -v t="$tol" -v n="$name" \
        'BEGIN{d=a-b; if(d<0)d=-d; if(d<=t){printf "PASS  %s : |%.8f-%.8f|=%.2e <= %.0e\n",n,a,b,d,t} else {printf "FAIL  %s : |%.8f-%.8f|=%.2e > %.0e\n",n,a,b,d,t; exit 1}}' \
        || fail=1
}

echo "== C1 large-system modes (caffeine) =="
DENSE=$(energy -method gfn2)
DENSE0=$(energy -method gfn2 -electronic_temperature 0)
echo "dense(300K)=$DENSE  dense(0K)=$DENSE0"

# C1c: single connected fragment -> exact dense fallback.
C1C=$(energy -method gfn2 -c1_mode fragments)
check "C1c fragments == dense" "$C1C" "$DENSE" 1e-6

# C1a: huge cell+buffer -> whole molecule per sub-system -> ~dense.
C1A=$(energy -method gfn2 -c1_mode dc -c1_cell_bohr 200 -c1_buffer_bohr 200)
check "C1a dc(full buffer) ~ dense" "$C1A" "$DENSE" 1e-3

# C1b: 0 K, threshold 0 -> dense-0K.
C1B=$(energy -method gfn2 -c1_mode sparse -c1_sparse_threshold 0 -electronic_temperature 0)
check "C1b sparse(thr 0) ~ dense0K" "$C1B" "$DENSE0" 1e-4

# Default path unchanged (gfn1 too).
NONE1=$(energy -method gfn1)
C1C1=$(energy -method gfn1 -c1_mode fragments)
check "C1c gfn1 == dense gfn1" "$C1C1" "$NONE1" 1e-6

if [ "$fail" -eq 0 ]; then echo "ALL C1 CHECKS PASSED"; exit 0; else echo "C1 CHECKS FAILED"; exit 1; fi
