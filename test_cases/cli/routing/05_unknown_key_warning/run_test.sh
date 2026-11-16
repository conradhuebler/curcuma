#!/bin/bash
# CLI Routing Test 05: a key that no module knows is reported, known keys are not
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026), docs/debt_2026-11/PARAMS.md M1
#
# `-sp mol.xyz -temprature 400` was accepted without a message. The run is not stopped:
# some capabilities read keys that have no PARAM.

set -e
export LC_NUMERIC=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="routing - 05: unknown key warning"
test_header "$TEST_NAME"

cd "$SCRIPT_DIR"

check() {
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "$2" = "0" ]; then
        echo -e "${GREEN}✓ PASS${NC}: $1"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: $1"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

"$CURCUMA" -sp water.xyz -method gfnff -temprature 400 -no_bmt > typo.log 2>&1 && rc=0 || rc=$?
check "misspelled flat flag still lets the run finish" $rc
grep -q "Parameter 'temprature' is unknown" typo.log && r=0 || r=1
check "-temprature is reported" $r

"$CURCUMA" -sp water.xyz -method gfnff -gfnff.cn_cutof 5 -no_bmt > typo2.log 2>&1 && rc=0 || rc=$?
grep -q "Parameter 'cn_cutof' is unknown" typo2.log && r=0 || r=1
check "misspelled -gfnff.cn_cutof is reported" $r

"$CURCUMA" -sp water.xyz -method gfnff -threads 1 -gpu none -dispersion false -gradient true -charge 0 -no_bmt > known.log 2>&1 && rc=0 || rc=$?
check "run with known flags finishes" $rc
grep -q "is unknown" known.log && r=1 || r=0
check "known flags (-threads, -gpu, -dispersion, -gradient, -charge) are not reported" $r

rm -f typo.log typo2.log known.log water*.topo.json curcuma_citations.bib water.sp.xyz

echo
echo "Tests run: $TESTS_RUN, passed: $TESTS_PASSED, failed: $TESTS_FAILED"
if [ "$TESTS_FAILED" -gt 0 ]; then
    exit 1
fi
exit 0
