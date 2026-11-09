#!/bin/bash
# CLI Routing Test 04: numeric values of Bool flags are booleans
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026), docs/debt_2026-11/PARAMS.md H4
#
# H4: CLI2Json stores `-write_trajectory 0|1` as a number; `-opt` then threw a nlohmann
#     type_error and ended without optimising.

set -e
export LC_NUMERIC=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="routing - 04: numeric Bool flags"
test_header "$TEST_NAME"

cd "$SCRIPT_DIR"

check() {
    # check <label> <condition result 0|1>
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "$2" = "0" ]; then
        echo -e "${GREEN}✓ PASS${NC}: $1"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: $1"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

# --- H4 ---
for v in 0 1 true false; do
    rm -f water.opt.xyz
    "$CURCUMA" -opt water.xyz -method gfnff -write_trajectory $v -no_bmt > h4_$v.log 2>&1
    rc=$?
    [ "$rc" = "0" ] && [ -f water.opt.xyz ] && r=0 || r=1
    check "-opt -write_trajectory $v optimises (exit 0, water.opt.xyz written)" $r
done

"$CURCUMA" -opt water.xyz -method gfnff -write_trajectory 0 -export_run wt.json -no_bmt > /dev/null 2>&1 || true
got=$(python3 -c "import json; print(json.load(open('wt.json'))['opt']['write_trajectory'])" 2>/dev/null || echo PYTHON_ERROR)
[ "$got" = "False" ] && r=0 || r=1
check "-write_trajectory 0 is exported as boolean false (got '$got')" $r

rm -f h4_*.log wt.json water.opt.xyz water*.topo.json water.param.json curcuma_citations.bib
rm -rf water.opt.*

echo
echo "Tests run: $TESTS_RUN, passed: $TESTS_PASSED, failed: $TESTS_FAILED"
if [ "$TESTS_FAILED" -gt 0 ]; then
    exit 1
fi
exit 0
