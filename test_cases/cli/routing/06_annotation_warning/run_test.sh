#!/bin/bash
# CLI Routing Test 06: values outside the enum/min/max annotation of a PARAM are reported
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026), docs/debt_2026-11/PARAMS.md L4
#
# The annotations were only printed by -help-module. A value outside them gives a warning,
# the run is not stopped; enumerations compare case-insensitively.

set -e
export LC_NUMERIC=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="routing - 06: annotation warning"
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

# -sp does not read the simplemd scope, so the values only reach the check
"$CURCUMA" -sp water.xyz -method gfnff -md.thermostat banana -md.andersen_probability 3 -md.temperature -5 -no_bmt > bad.log 2>&1 && rc=0 || rc=$?
check "run with out-of-range values still finishes" $rc
grep -q "Value 'banana' of -thermostat (module simplemd) is not one of" bad.log && r=0 || r=1
check "enum violation is reported" $r
grep -q "of -andersen_probability (module simplemd) is above the maximum" bad.log && r=0 || r=1
check "maximum violation is reported" $r
grep -q "of -temperature (module simplemd) is below the minimum" bad.log && r=0 || r=1
check "minimum violation is reported" $r

"$CURCUMA" -sp water.xyz -method gfnff -md.thermostat CSVR -md.andersen_probability 0.5 -md.temperature 300 -no_bmt > ok.log 2>&1 && rc=0 || rc=$?
check "run with valid values finishes" $rc
grep -qE "is not one of|below the minimum|above the maximum" ok.log && r=1 || r=0
check "valid values (case-insensitive enum) are not reported" $r

rm -f bad.log ok.log water*.topo.json curcuma_citations.bib water.sp.xyz

echo
echo "Tests run: $TESTS_RUN, passed: $TESTS_PASSED, failed: $TESTS_FAILED"
if [ "$TESTS_FAILED" -gt 0 ]; then
    exit 1
fi
exit 0
