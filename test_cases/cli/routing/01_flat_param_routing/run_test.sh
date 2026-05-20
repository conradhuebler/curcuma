#!/bin/bash
# CLI Routing Test 01: Flat CLI flags route to their owning module via ParameterRegistry
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Auto-routing + ambiguity + global passthrough verification

set -e
export LC_NUMERIC=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="routing - 01: Flat CLI flag routing"
test_header "$TEST_NAME"

cd "$SCRIPT_DIR"

# --- 1. Flat registered gfnff param routes to controller["gfnff"] ---
"$CURCUMA" -sp water.xyz -method gfnff -cn_cutoff_bohr 5.5 -export_run flat.json -silent > /dev/null 2>&1 || true
TESTS_RUN=$((TESTS_RUN + 1))
if [ ! -f flat.json ]; then
    echo -e "${RED}✗ FAIL${NC}: -export_run did not produce flat.json"
    TESTS_FAILED=$((TESTS_FAILED + 1))
else
    routed=$(python3 -c 'import json,sys; d=json.load(open("flat.json")); print(d.get("gfnff",{}).get("cn_cutoff_bohr"))')
    if [ "$routed" = "5.5" ]; then
        echo -e "${GREEN}✓ PASS${NC}: -cn_cutoff_bohr 5.5 routed to gfnff.cn_cutoff_bohr"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: gfnff.cn_cutoff_bohr=$routed (expected 5.5)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
fi

# --- 2. Dotted form still works (regression) ---
"$CURCUMA" -sp water.xyz -method gfnff -gfnff.cn_cutoff_bohr 4.5 -export_run dotted.json -silent > /dev/null 2>&1 || true
TESTS_RUN=$((TESTS_RUN + 1))
dotted=$(python3 -c 'import json; d=json.load(open("dotted.json")); print(d.get("gfnff",{}).get("cn_cutoff_bohr"))' 2>/dev/null || echo "MISSING")
if [ "$dotted" = "4.5" ]; then
    echo -e "${GREEN}✓ PASS${NC}: dotted -gfnff.cn_cutoff_bohr 4.5 still routes correctly"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    echo -e "${RED}✗ FAIL${NC}: dotted form broken, gfnff.cn_cutoff_bohr=$dotted (expected 4.5)"
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

# --- 3. Global params stay top-level (verbosity 1 must not be re-routed; values land at controller root) ---
"$CURCUMA" -sp water.xyz -method gfnff -verbosity 1 -threads 2 -export_run globals.json > /dev/null 2>&1 || true
TESTS_RUN=$((TESTS_RUN + 1))
result=$(python3 -c 'import json
d=json.load(open("globals.json"))
v=d.get("verbosity"); t=d.get("threads")
ok = (v in (1, 1.0, "1")) and (t in (2, 2.0, "2"))
print("PASS" if ok else f"FAIL v={v} t={t}")' 2>/dev/null || echo "PYTHON_ERROR")
if [ "$result" = "PASS" ]; then
    echo -e "${GREEN}✓ PASS${NC}: global params -verbosity, -threads remain at top level"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    echo -e "${RED}✗ FAIL${NC}: $result"
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

# --- 4. Unregistered flat flag stays in command's module ---
"$CURCUMA" -sp water.xyz -method gfnff -some_unknown_flag true -export_run unreg.json -silent > /dev/null 2>&1 || true
TESTS_RUN=$((TESTS_RUN + 1))
unreg=$(python3 -c 'import json; d=json.load(open("unreg.json")); print(d.get("opt",{}).get("some_unknown_flag"))' 2>/dev/null || echo "MISSING")
if [ "$unreg" = "True" ]; then
    echo -e "${GREEN}✓ PASS${NC}: unregistered -some_unknown_flag preserved in command module"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    echo -e "${RED}✗ FAIL${NC}: unregistered flag handling unexpected (got $unreg)"
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

# --- 5. Export includes full gfnff defaults block ---
"$CURCUMA" -sp water.xyz -method gfnff -export_run defaults.json -silent > /dev/null 2>&1 || true
TESTS_RUN=$((TESTS_RUN + 1))
n_gfnff_keys=$(python3 -c 'import json; d=json.load(open("defaults.json")); print(len(d.get("gfnff",{})))' 2>/dev/null || echo "0")
if [ "$n_gfnff_keys" -ge 10 ]; then
    echo -e "${GREEN}✓ PASS${NC}: exported gfnff block has $n_gfnff_keys params (full defaults present)"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    echo -e "${RED}✗ FAIL${NC}: gfnff block only has $n_gfnff_keys params (expected >= 10)"
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

# --- Summary ---
echo
echo "Tests run: $TESTS_RUN, passed: $TESTS_PASSED, failed: $TESTS_FAILED"
if [ "$TESTS_FAILED" -gt 0 ]; then
    exit 1
fi
exit 0
