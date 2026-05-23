#!/bin/bash
# Test: WP-S2 MD diagnostics JSONL dump (per-step charges/CN/energy/grad-norm)
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (May 2026) — Validates:
#   1. md_diagnostics=true produces <basename>.diag.jsonl
#   2. Number of records matches max_time / (dt * dump_frequency)
#   3. Each line is valid JSON with required keys
#   4. charges.length == cn.length == gradient_norm.length == atom_count
#   5. energy contains expected GFN-FF terms

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 11: MD diagnostics JSONL dump (WP-S2)"
TEST_DIR="$SCRIPT_DIR"

# MD parameters — small caffeine system (24 atoms), short run
MD_MAXTIME=100            # 100 fs
MD_DT=1.0
MD_DUMP=10                # dump_frequency = 10 → ~11 records (steps 0,10,20,...,100)
MD_SEED=42
ATOM_COUNT=24

DIAG_FILE="input.diag.jsonl"
EXPECTED_MIN_RECORDS=10   # tolerate slight off-by-one across versions

run_test() {
    cd "$TEST_DIR"
    rm -f "$DIAG_FILE" input.trj.xyz input.opt.xyz input.restart stdout.log stderr.log

    echo "Running: $CURCUMA -md input.xyz -method gfnff -maxtime $MD_MAXTIME"
    echo "         -md.md_diagnostics true -md.dump_frequency $MD_DUMP -md.seed $MD_SEED"
    timeout 120 $CURCUMA \
        -md input.xyz \
        -method gfnff \
        -maxtime $MD_MAXTIME \
        -threads 1 \
        -md.dt $MD_DT \
        -md.dump_frequency $MD_DUMP \
        -md.md_diagnostics true \
        -md.seed $MD_SEED \
        -md.no_restart \
        -md.rattle_12 false \
        -md.thermostat none \
        > stdout.log 2> stderr.log
    return $?
}

validate_results() {
    local failed=0

    # 1. diag.jsonl exists and is non-empty
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -s "$DIAG_FILE" ]; then
        echo -e "${GREEN}✓ PASS${NC}: $DIAG_FILE exists and is non-empty"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: $DIAG_FILE missing or empty"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    # 2. Record count
    TESTS_RUN=$((TESTS_RUN + 1))
    local records
    records=$(wc -l < "$DIAG_FILE")
    if [ "${records:-0}" -ge "$EXPECTED_MIN_RECORDS" ]; then
        echo -e "${GREEN}✓ PASS${NC}: $records records (≥$EXPECTED_MIN_RECORDS expected)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: only $records records (need ≥$EXPECTED_MIN_RECORDS)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        failed=1
    fi

    # 3. Every line is valid JSON with required keys + correct array lengths
    TESTS_RUN=$((TESTS_RUN + 1))
    local validation_out
    validation_out=$(python3 - "$DIAG_FILE" "$ATOM_COUNT" <<'PY' 2>&1
import json, sys
path, expected_atoms = sys.argv[1], int(sys.argv[2])
required = {"step", "time_fs", "energy", "charges", "cn", "gradient_norm", "hb_count", "xb_count"}
required_energy_keys = {"Bond", "Angle", "Torsion", "Coulomb", "Dispersion"}
errors = []
with open(path) as f:
    for ln, line in enumerate(f, 1):
        try:
            rec = json.loads(line)
        except json.JSONDecodeError as e:
            errors.append(f"line {ln}: not valid JSON ({e})")
            continue
        missing = required - rec.keys()
        if missing:
            errors.append(f"line {ln}: missing keys {missing}")
        for arr_name in ("charges", "cn", "gradient_norm"):
            if arr_name in rec and len(rec[arr_name]) != expected_atoms:
                errors.append(f"line {ln}: {arr_name}.length={len(rec[arr_name])} != {expected_atoms}")
        if "energy" in rec:
            missing_e = required_energy_keys - rec["energy"].keys()
            if missing_e:
                errors.append(f"line {ln}: energy missing keys {missing_e}")
if errors:
    print("FAILED")
    for e in errors[:10]:
        print("  " + e)
    sys.exit(1)
else:
    print("OK")
PY
    )
    if echo "$validation_out" | grep -q "^OK$"; then
        echo -e "${GREEN}✓ PASS${NC}: Every line is valid JSON with required keys and correct array lengths"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: JSONL validation failed"
        echo "$validation_out"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        failed=1
    fi

    # 4. Step numbers are monotonically increasing
    TESTS_RUN=$((TESTS_RUN + 1))
    local steps_ok
    steps_ok=$(python3 - "$DIAG_FILE" <<'PY' 2>&1
import json, sys
last = -1
with open(sys.argv[1]) as f:
    for line in f:
        s = json.loads(line)["step"]
        if s <= last:
            print(f"FAIL at step {s} <= {last}")
            sys.exit(1)
        last = s
print("OK")
PY
    )
    if echo "$steps_ok" | grep -q "^OK$"; then
        echo -e "${GREEN}✓ PASS${NC}: Step numbers are monotonically increasing"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: Step ordering broken: $steps_ok"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        failed=1
    fi

    return $failed
}

main() {
    test_header "$TEST_NAME"

    run_test
    local run_exit=$?
    assert_exit_code $run_exit 0 "gfnff MD with diagnostics should complete without crash"

    validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
