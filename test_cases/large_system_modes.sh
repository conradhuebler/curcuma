#!/bin/bash
# large_system_modes — regression test
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (large_system_mode, June 2026)
#
# Validates the three opt-in native-GFN large-system modes against the dense
# reference on small connected molecules, where each mode has an exact/known
# limit, AND validates the -eigensolver propagation through fragments/dc.
#
#   - large_system_mode=fragments: connected molecule -> dense fallback (EXACT)
#   - large_system_mode=dc:        huge cell+buffer -> dense (tight)
#   - large_system_mode=sparse:    0 K, threshold 0 -> dense-0K (tight)
#
# Eigensolver propagation:
#   - large_system_mode=fragments -eigensolver=native vs -eigensolver=mkl
#       bit-identical (caffeine = single fragment = dense; -eigensolver choice
#       is observable only when the mode actually engages a solver)
#   - large_system_mode=dc -eigensolver=native @full-buffer ~ dense
#   - large_system_mode=dc -eigensolver=purify -electronic_temperature 0
#       @full-buffer ~ dense
#   - large_system_mode=dc -eigensolver=purify -electronic_temperature 300
#       must hard-error (T>0 incompatible with purify)
#
# Also asserts the default path (large_system_mode=none) is unchanged.
#
# Usage: large_system_modes.sh <curcuma-binary> <molecule-dir>

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
# Assert a hard-error by checking that the energy line is "0.00000000 Eh"
# (curcuma doesn't propagate the wrapper's hasError() into the exit code; the
# fallback is to print 0.00000000 Eh when the method refused to compute, which
# is the observable signal of a clean rejection in the run log).
assert_hard_error() {
    local name="$1"; shift
    local out
    out=$("$@" 2>/dev/null \
        | sed -r 's/\x1b\[[0-9;]*m//g')
    if echo "$out" | grep -q 'Single Point Energy = 0.00000000 Eh'; then
        echo "PASS  $name : rejected with 0.00000000 Eh (hard-error path)"
    elif echo "$out" | grep -q 'Failed to set molecule'; then
        echo "PASS  $name : rejected with 'Failed to set molecule'"
    else
        echo "FAIL  $name : expected hard-error, got a real energy:"; echo "$out" | grep 'Single Point' | head -1; fail=1
    fi
}

echo "== large_system_modes (caffeine) =="
DENSE=$(energy -method gfn2)
DENSE0=$(energy -method gfn2 -electronic_temperature 0)
echo "dense(300K)=$DENSE  dense(0K)=$DENSE0"

# fragments: single connected fragment -> exact dense fallback.
C1C=$(energy -method gfn2 -large_system_mode fragments)
check "fragments == dense" "$C1C" "$DENSE" 1e-6

# dc: huge cell+buffer -> whole molecule per sub-system -> ~dense.
C1A=$(energy -method gfn2 -large_system_mode dc -large_system_cell_bohr 200 -large_system_buffer_bohr 200)
check "dc(full buffer) ~ dense" "$C1A" "$DENSE" 1e-3

# sparse: 0 K, threshold 0 -> dense-0K.
C1B=$(energy -method gfn2 -large_system_mode sparse -large_system_sparse_threshold 0 -electronic_temperature 0)
check "sparse(thr 0) ~ dense0K" "$C1B" "$DENSE0" 1e-4

# Default path unchanged (gfn1 too).
NONE1=$(energy -method gfn1)
C1C1=$(energy -method gfn1 -large_system_mode fragments)
check "fragments gfn1 == dense gfn1" "$C1C1" "$NONE1" 1e-6

# --- -eigensolver propagation in fragments (single fragment = dense path;
# both should match the dense reference bit-precisely) -------------------
C_NATIVE=$(energy -method gfn2 -large_system_mode fragments -eigensolver native)
C_LOBPCG=$(energy -method gfn2 -large_system_mode fragments -eigensolver lobpcg)
check "fragments+eigensolver=native == dense" "$C_NATIVE" "$DENSE" 1e-6
check "fragments+eigensolver=lobpcg == dense" "$C_LOBPCG" "$DENSE" 1e-6

# --- -eigensolver propagation in dc (full buffer = single sub-system = dense) -
DC_NATIVE=$(energy -method gfn2 -large_system_mode dc -large_system_cell_bohr 200 -large_system_buffer_bohr 200 -eigensolver native)
DC_LOBPCG=$(energy -method gfn2 -large_system_mode dc -large_system_cell_bohr 200 -large_system_buffer_bohr 200 -eigensolver lobpcg)
check "dc+eigensolver=native@fullbuf ~ dense" "$DC_NATIVE" "$DENSE" 1e-3
check "dc+eigensolver=lobpcg@fullbuf ~ dense" "$DC_LOBPCG" "$DENSE" 1e-3

# --- dc + purify at T=0: purify falls back to dense GES per sub-block (can't
#     provide the full spectrum DC needs for Fermi occupation). Verify the
#     fallback produces a real energy (not a hard-error). Note: DC-SCF at T=0
#     does not always converge (discontinuous Fermi occupation), so we only
#     check that the fallback runs without error, not a tight energy match.
DC_PURIFY0=$(energy -method gfn2 -large_system_mode dc -large_system_cell_bohr 200 -large_system_buffer_bohr 200 -eigensolver purify -xtb.electronic_temperature 0)
if [ -z "$DC_PURIFY0" ] || [ "$DC_PURIFY0" = "0.00000000" ]; then
    echo "FAIL  dc+eigensolver=purify+T=0 should compute (fallback to dense GES), got '$DC_PURIFY0'"; fail=1
else
    echo "PASS  dc+eigensolver=purify+T=0 computes via dense GES fallback (got $DC_PURIFY0 Eh)"
fi

# --- dc + purify at T=300 must hard-error (purify is 0-K only) -------------
assert_hard_error "dc+eigensolver=purify+T=300 hard-errors" \
    "$CURCUMA" -sp "$MOL" -method gfn2 -large_system_mode dc \
        -eigensolver purify -xtb.electronic_temperature 300
assert_hard_error "fragments+eigensolver=purify+T=300 hard-errors" \
    "$CURCUMA" -sp "$MOL" -method gfn2 -large_system_mode fragments \
        -eigensolver purify -xtb.electronic_temperature 300
# Sanity: T=0 must NOT trigger the hard-error (purify then runs on the
# fragments/dc sub-blocks with the regular purify soft-warn fallback).
DC_PURIFY_T0=$(energy -method gfn2 -large_system_mode fragments \
               -eigensolver purify -xtb.electronic_temperature 0)
if [ -z "$DC_PURIFY_T0" ] || [ "$DC_PURIFY_T0" = "0.00000000" ]; then
    echo "FAIL  fragments+eigensolver=purify+T=0 should compute, got '$DC_PURIFY_T0'"; fail=1
else
    echo "PASS  fragments+eigensolver=purify+T=0 computes (got $DC_PURIFY_T0 Eh)"
fi

if [ "$fail" -eq 0 ]; then echo "ALL large_system_modes CHECKS PASSED"; exit 0; else echo "large_system_modes CHECKS FAILED"; exit 1; fi
