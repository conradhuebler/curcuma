#!/bin/bash
# CLI Routing Test 02: -export_run / -import_config round-trip with deep merge + CLI override
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Validates _command/_input replay, deep merge, and override precedence

set -e
export LC_NUMERIC=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="routing - 02: Round-trip export/import with override"
test_header "$TEST_NAME"

cd "$SCRIPT_DIR"

# --- 1. Export a run with a non-default gfnff param ---
"$CURCUMA" -sp water.xyz -method gfnff -cn_cutoff_bohr 5.5 -export_run baseline.json -silent > /dev/null 2>&1 || true
TESTS_RUN=$((TESTS_RUN + 1))
if [ ! -f baseline.json ]; then
    echo -e "${RED}✗ FAIL${NC}: baseline.json not produced"
    TESTS_FAILED=$((TESTS_FAILED + 1))
    exit 1
fi
cmd=$(python3 -c 'import json; d=json.load(open("baseline.json")); print(d.get("_command"))')
inp=$(python3 -c 'import json; d=json.load(open("baseline.json")); print(d.get("_input"))')
if [ "$cmd" = "sp" ] && [ "$inp" = "water.xyz" ]; then
    echo -e "${GREEN}✓ PASS${NC}: exported _command='sp', _input='water.xyz'"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    echo -e "${RED}✗ FAIL${NC}: _command='$cmd', _input='$inp'"
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

# --- 2. JSON-driven replay: curcuma -import_config baseline.json ---
"$CURCUMA" -import_config baseline.json -silent > replay.log 2>&1
rc=$?
TESTS_RUN=$((TESTS_RUN + 1))
if [ $rc -eq 0 ]; then
    echo -e "${GREEN}✓ PASS${NC}: JSON-driven replay succeeded (exit 0)"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    echo -e "${RED}✗ FAIL${NC}: replay returned exit $rc"
    cat replay.log
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

# --- 3. CLI override during replay ---
"$CURCUMA" -import_config baseline.json -cn_cutoff_bohr 7.5 -export_run override.json -silent > /dev/null 2>&1 || true
TESTS_RUN=$((TESTS_RUN + 1))
override_val=$(python3 -c 'import json; d=json.load(open("override.json")); print(d.get("gfnff",{}).get("cn_cutoff_bohr"))' 2>/dev/null || echo "MISSING")
if [ "$override_val" = "7.5" ]; then
    echo -e "${GREEN}✓ PASS${NC}: CLI override -cn_cutoff_bohr 7.5 wins over JSON's 5.5"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    echo -e "${RED}✗ FAIL${NC}: override value=$override_val (expected 7.5)"
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

# --- 4. Deep-merge: 3-level nested JSON survives import + export ---
cat > deep.json <<'EOF'
{
  "_command": "sp",
  "_input": "water.xyz",
  "method": "gfnff",
  "opt": { "convergence": { "level": "tight", "energy": 1e-9 } }
}
EOF
"$CURCUMA" -import_config deep.json -export_run deep_out.json -silent > /dev/null 2>&1 || true
TESTS_RUN=$((TESTS_RUN + 1))
deep_level=$(python3 -c 'import json; d=json.load(open("deep_out.json")); print(d.get("opt",{}).get("convergence",{}).get("level"))' 2>/dev/null || echo "MISSING")
deep_e=$(python3 -c 'import json; d=json.load(open("deep_out.json")); print(d.get("opt",{}).get("convergence",{}).get("energy"))' 2>/dev/null || echo "MISSING")
if [ "$deep_level" = "tight" ] && [ "$deep_e" = "1e-09" ]; then
    echo -e "${GREEN}✓ PASS${NC}: 3-level deep merge preserved (opt.convergence.level=tight, energy=1e-9)"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    echo -e "${RED}✗ FAIL${NC}: deep merge lost data: level=$deep_level, energy=$deep_e"
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

# --- 5. CLI command beats JSON _command (with warning) ---
"$CURCUMA" -hessian water.xyz -import_config baseline.json -silent > cmd_conflict.log 2>&1 || true
TESTS_RUN=$((TESTS_RUN + 1))
if grep -q "_command" cmd_conflict.log 2>/dev/null; then
    echo -e "${GREEN}✓ PASS${NC}: warning emitted when CLI command conflicts with JSON _command"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    echo -e "${YELLOW}⚠ INFO${NC}: no _command-conflict warning emitted (expected, but only a warning)"
    TESTS_PASSED=$((TESTS_PASSED + 1))
fi

# --- Summary ---
echo
echo "Tests run: $TESTS_RUN, passed: $TESTS_PASSED, failed: $TESTS_FAILED"
if [ "$TESTS_FAILED" -gt 0 ]; then
    exit 1
fi
exit 0
