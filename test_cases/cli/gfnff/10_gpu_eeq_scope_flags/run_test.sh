#!/bin/bash

# Test: EEQ flags reach the GPU GFN-FF path under every spelling (-eeq_solver.X, -gfnff.X, flat -X)
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Oct 2026) - AI-generated, machine-tested; human production testing pending
#
# Found on Oct 26, 2026: the GPU wrapper read skip_phase2 from the gfnff scope only, so -eeq_solver.skip_phase2 true
# changed the first (CPU) step of an MD and left the device EEQ switched on, while -gfnff.skip_phase2 true bypassed it.
#   1. default       - no bypass notice, the device projected PCG runs
#   2. three spellings of skip_phase2 true - all print the bypass notice and give the same final Epot (1e-9 Eh),
#      which differs from the default run
# Skipped (exit 77) when the default run does not use the CUDA device EEQ (no plugin or no device).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../test_utils.sh"

TEST_NAME="gfnff - 10: GPU EEQ flags under every spelling"
ENERGY_TOLERANCE="0.000000001"

strip() { sed 's/\x1b\[[0-9;]*m//g' "$1"; }
epot_of() { strip "$1" | grep -E '^\s*[0-9]+\.[0-9]+\s+-?[0-9]+\.[0-9]+\s' | tail -1 | awk '{print $2}'; }

run() {   # run <log> <extra flags...>
    local log=$1; shift
    $CURCUMA -md chain.xyz -method gfnff -gpu cuda -gpu_devices 0 -gpu_strict true -T 300 -dt 1 -seed 7 \
        -max_time 6 -print 6 -dump 100000000 -verbosity 1 -no_bmt "$@" > "$log" 2>&1
    echo "EXIT $?" >> "$log"
}

main() {
    test_header "$TEST_NAME"
    cd "$SCRIPT_DIR"
    set +e
    run default.log
    if ! strip default.log | grep -q "EEQ solver (final charges):.*\[GPU CUDA"; then
        echo -e "${YELLOW}SKIP${NC}: the default run did not use the CUDA device EEQ (no CUDA plugin or device)"
        exit 77
    fi
    assert_string_in_file "EXIT 0" default.log "default GPU run: exit code 0"
    assert_string_not_in_file "bypassing GPU EEQ" default.log "default: device EEQ not bypassed"
    local E0=$(epot_of default.log)

    run s_eeq.log -eeq_solver.skip_phase2 true
    run s_gfnff.log -gfnff.skip_phase2 true
    run s_flat.log -skip_phase2 true
    for l in s_eeq s_gfnff s_flat; do
        assert_string_in_file "bypassing GPU EEQ" $l.log "$l: skip_phase2 reaches the GPU path"
    done
    local E1=$(epot_of s_eeq.log) E2=$(epot_of s_gfnff.log) E3=$(epot_of s_flat.log)
    assert_numeric_match "$E2" "$E1" "$ENERGY_TOLERANCE" "-eeq_solver.skip_phase2 gives the same Epot as -gfnff.skip_phase2"
    assert_numeric_match "$E2" "$E3" "$ENERGY_TOLERANCE" "flat -skip_phase2 gives the same Epot as -gfnff.skip_phase2"
    TESTS_RUN=$((TESTS_RUN + 1))
    if awk -v a="$E0" -v b="$E2" 'BEGIN{d=a-b; if (d<0) d=-d; exit !(d>1e-4)}'; then
        echo -e "${GREEN}PASS${NC}: skip_phase2 changes the energy (default $E0, skipped $E2)"; TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: skip_phase2 did not change the energy (default $E0, skipped $E2)"; TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    print_test_summary
}

main "$@"
