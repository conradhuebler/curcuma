#!/bin/bash
# CLI Routing Test 03: ConfSearch parameters reach ConfSearch
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - regression guard for the ConfSearch registry migration (Jul 2026)
#
# Before ConfSearch was migrated into the ParameterRegistry it owned no parameter names, so
# main.cpp's auto-router silently moved its flags into whichever module DID own them:
#   -md_method / -opt_method  -> controller["polymerbuild"]
#   -thermostat / -coupling   -> controller["simplemd"]
#   -restart                  -> controller["confscan"]
# ConfSearch then fell back to its "gfnff" default and reported "Single-method mode" even
# though the user had asked for a dual-method run. The dotted escape hatch did not help
# either, because "-confsearch.opt_method" is stripped to "opt_method" before the router runs.
#
# This test pins BOTH directions: ConfSearch flags must stay in controller["confsearch"], and
# the SimpleMD-legacy names that no module registers must still reach "-md".

set -e
export LC_NUMERIC=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="routing - 03: ConfSearch parameter routing"
test_header "$TEST_NAME"

cd "$SCRIPT_DIR"

# This test asserts on the resolved controller written by -export_run, not on any search
# result -- but curcuma still executes the command, so keep every run trivially short.
FAST="-method gfnff -time 5 -repeat 1 -startT 300 -endT 300 -deltaT 100 -threads 1"

check_key() {
    # check_key <label> <json> <python-expr> <expected>
    local label="$1" file="$2" expr="$3" expected="$4"
    TESTS_RUN=$((TESTS_RUN + 1))
    local got
    got=$(python3 -c "import json; d=json.load(open('$file')); print($expr)" 2>/dev/null || echo "PYTHON_ERROR")
    if [ "$got" = "$expected" ]; then
        echo -e "${GREEN}✓ PASS${NC}: $label"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: $label (got '$got', expected '$expected')"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

# --- ConfSearch flags must land in controller["confsearch"] ---
"$CURCUMA" -confsearch water.xyz -method gfnff -md_method gfnff -opt_method gfn2 \
    -thermostat berendsen -coupling 25 -restart -charge -1 -spin 1 -rmsd 0.5 \
    -xtb.solvent water $FAST -export_run cs.json -no_bmt -silent > /dev/null 2>&1 || true

if [ ! -f cs.json ]; then
    echo -e "${RED}✗ FAIL${NC}: -export_run did not produce cs.json"
    TESTS_RUN=$((TESTS_RUN + 1))
    TESTS_FAILED=$((TESTS_FAILED + 1))
else
    check_key "-opt_method reaches confsearch (not polymerbuild)" cs.json "d['confsearch']['opt_method']" "gfn2"
    check_key "-md_method reaches confsearch (not polymerbuild)"  cs.json "d['confsearch']['md_method']"  "gfnff"
    check_key "-thermostat reaches confsearch (not simplemd)"     cs.json "d['confsearch']['thermostat']" "berendsen"
    check_key "-coupling reaches confsearch (not simplemd)"       cs.json "d['confsearch']['coupling']"   "25.0"
    check_key "-restart reaches confsearch (not confscan)"        cs.json "d['confsearch']['restart']"    "True"
    check_key "-rmsd reaches confsearch"                          cs.json "d['confsearch']['rmsd']"       "0.5"
    check_key "-charge reaches confsearch"                        cs.json "d['confsearch']['charge']"     "-1.0"
    check_key "-spin reaches confsearch"                          cs.json "d['confsearch']['spin']"       "1.0"
    check_key "method sub-scope preserved"                        cs.json "d['xtb']['solvent']"           "water"
    check_key "nothing leaked into polymerbuild"                  cs.json "'polymerbuild' in d"           "False"
    check_key "nothing leaked into simplemd.thermostat"           cs.json "d.get('simplemd',{}).get('thermostat')" "None"
    check_key "nothing leaked into confscan.restart"              cs.json "d.get('confscan',{}).get('restart')"    "None"
fi

# --- Aliases declared in the ConfSearch PARAM block resolve to their canonical names ---
"$CURCUMA" -confsearch water.xyz -dt 0.25 -velo 2 -dump 7 $FAST -export_run alias.json -no_bmt -silent > /dev/null 2>&1 || true
if [ -f alias.json ]; then
    check_key "-dt stays in confsearch (alias of time_step)"   alias.json "d['confsearch']['dt']"   "0.25"
    check_key "-velo stays in confsearch"                      alias.json "d['confsearch']['velo']" "2.0"
    check_key "-dump stays in confsearch"                      alias.json "d['confsearch']['dump']" "7.0"
fi

# --- Non-regression: SimpleMD-legacy names that NO module registers must still reach -md ---
# ConfSearch deliberately does not register these; claiming them would steal them from -md.
"$CURCUMA" -md water.xyz -unique true -respa 2 -dipole true -thermostat csvr -method gfnff -maxtime 5 \
    -export_run md.json -no_bmt -silent > /dev/null 2>&1 || true
if [ -f md.json ]; then
    check_key "-unique still reaches simplemd"     md.json "d['simplemd']['unique']"     "True"
    check_key "-respa still reaches simplemd"      md.json "d['simplemd']['respa']"      "2.0"
    check_key "-dipole still reaches simplemd"     md.json "d['simplemd']['dipole']"     "True"
    check_key "-thermostat still reaches simplemd" md.json "d['simplemd']['thermostat']" "csvr"
    check_key "not stolen by confsearch"           md.json "'unique' in d.get('confsearch',{})" "False"
fi

# --- Cleanup ---
rm -f cs.json alias.json md.json
rm -f water.bias*.xyz water.cumulative*.xyz water.input*.xyz water.mtd.xyz water.trj.xyz
rm -f water.r*.mtd.xyz water*.topo.json water.param.json water_md_params.json
rm -f water_child_params.json water.confsearch.restart.json curcuma_restart.json
rm -f curcuma_citations.bib water.md.restart.json
rm -rf water.snapshots water.r*.snapshots water.confsearch.* water.md.*

# --- Summary ---
echo
echo "Tests run: $TESTS_RUN, passed: $TESTS_PASSED, failed: $TESTS_FAILED"
if [ "$TESTS_FAILED" -gt 0 ]; then
    exit 1
fi
exit 0
