#!/bin/bash

# Test: RATTLE must not freeze inter-fragment motion (RMSD-MTD bias under constraints)
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Jul 2026) - regression guard for the RATTLE bug where SimpleMD::Rattle() called
# RemoveRotations() unconditionally on every constrained step. That zeroed EVERY fragment's linear
# and angular momentum each step, so the RMSD-MTD bias could not move two fragments relative to each
# other -- "the bias does not grip" whenever RATTLE was on (which ConfSearch does automatically for
# every cycle at T >= rattle_threshold_temp).
#
# Measured A/B on this exact input (3 ps, gfnff, T=500 K, seed 42, rattle 2):
#   with the RemoveRotations call: fragment-COM distance range 0.43 A,  44 hills deposited
#   without it (fixed):            fragment-COM distance range 2.17 A, 127 hills deposited
# The thresholds below sit between the two.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 15: RATTLE + RMSD-MTD inter-fragment motion"
TEST_DIR="$SCRIPT_DIR"

MIN_COM_RANGE=1.0    # Angstrom, spanned by the fragment-COM distance over the trajectory
MIN_HILLS=80         # deposited RMSD-MTD hills
MAX_CH_RANGE=0.01    # Angstrom, allowed variation of a RATTLE-constrained C-H bond

run_test() {
    cd "$TEST_DIR"
    cleanup_bmt_dirs
    rm -f stdout.log stderr.log
    timeout 200 "$CURCUMA" -md input.xyz -method gfnff -rmsd_mtd \
        -maxtime 3000 -temperature 500 -rattle 2 -md.seed 42 -threads 1 \
        > stdout.log 2> stderr.log
    RUN_EXIT=$?
    TRJ_FILE=$(find_output_file "input.trj.xyz")
}

validate_results() {
    local failed=0

    assert_exit_code $RUN_EXIT 0 "constrained MTD run completes"

    # 1. RATTLE really active (8 X-H bonds in the acetic acid dimer)
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "RATTLE: 8 constraints" stdout.log 2>/dev/null; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: RATTLE active with 8 X-H constraints"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: RATTLE constraint report missing (constraints not applied?)"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    if [ -z "$TRJ_FILE" ] || [ ! -s "$TRJ_FILE" ]; then
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: trajectory not found"
        TESTS_RUN=$((TESTS_RUN + 1)); TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    # 2. Fragment-COM distance must span more than MIN_COM_RANGE.
    # Atoms 1-8 are molecule A, atoms 9-16 molecule B (unweighted centroids are enough here).
    local com_range
    com_range=$(awk '
        NF==1 && $1==16 { n=0; getline; next }                 # atom count line, skip comment
        NF>=4 {
            n++
            if (n<=8)  { ax+=$2; ay+=$3; az+=$4 } else { bx+=$2; by+=$3; bz+=$4 }
            if (n==16) {
                dx=ax/8-bx/8; dy=ay/8-by/8; dz=az/8-bz/8
                d=sqrt(dx*dx+dy*dy+dz*dz)
                if (first==0) { mn=d; mx=d; first=1 }
                if (d<mn) mn=d; if (d>mx) mx=d
                ax=ay=az=bx=by=bz=0; n=0
            }
        }
        END { printf "%.4f", mx-mn }' "$TRJ_FILE")
    TESTS_RUN=$((TESTS_RUN + 1))
    if awk -v v="$com_range" -v t="$MIN_COM_RANGE" 'BEGIN{exit(v>t?0:1)}'; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: fragment-COM distance spans $com_range A (> $MIN_COM_RANGE A required)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: fragment-COM distance spans only $com_range A (<= $MIN_COM_RANGE A) -- inter-fragment motion frozen"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 3. The bias keeps depositing hills (exploration is alive)
    local hills
    hills=$(grep -oE "RMSD-MTD provenance: [0-9]+ deposits" stdout.log | grep -oE "[0-9]+" | head -1)
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "${hills:-0}" -ge "$MIN_HILLS" ]; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: ${hills} RMSD-MTD hills deposited (>= $MIN_HILLS required)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: only ${hills:-0} RMSD-MTD hills deposited (< $MIN_HILLS)"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 4. The constraint itself still holds: C2-H5 (atoms 2 and 5) must stay fixed
    local ch_range
    ch_range=$(awk '
        NF==1 && $1==16 { n=0; getline; next }
        NF>=4 {
            n++
            if (n==2) { cx=$2; cy=$3; cz=$4 }
            if (n==5) { hx=$2; hy=$3; hz=$4 }
            if (n==16) {
                dx=cx-hx; dy=cy-hy; dz=cz-hz
                d=sqrt(dx*dx+dy*dy+dz*dz)
                if (first==0) { mn=d; mx=d; first=1 }
                if (d<mn) mn=d; if (d>mx) mx=d
                n=0
            }
        }
        END { printf "%.5f", mx-mn }' "$TRJ_FILE")
    TESTS_RUN=$((TESTS_RUN + 1))
    if awk -v v="$ch_range" -v t="$MAX_CH_RANGE" 'BEGIN{exit(v<t?0:1)}'; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: constrained C-H bond varies by $ch_range A (< $MAX_CH_RANGE A)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: constrained C-H bond varies by $ch_range A (>= $MAX_CH_RANGE A) -- RATTLE not enforcing"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    return $failed
}

main() {
    test_header "$TEST_NAME"
    run_test
    validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
