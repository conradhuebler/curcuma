#!/bin/bash
export PROJECT_ROOT='/home/conrad/src/claude_curcuma/curcuma/test_cases/cli/../..'

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    $CURCUMA -rmsd ref.xyz target.xyz > stdout.log 2> stderr.log
    assert_exit_code $? 0 "RMSD should succeed"
    assert_string_in_file "RMSD" stdout.log "RMSD output"
    return 0
}

validate_results() {
    local rmsd_line=$(grep "RMSD:" stdout.log | sed 's/\x1b\[[0-9;]*m//g' | head -1)
    if [ -n "$rmsd_line" ]; then
        local rmsd_value=$(echo "$rmsd_line" | grep -oP 'RMSD:\s+\K[0-9.]+' | head -1)
        if [ -n "$rmsd_value" ]; then
            echo -e "${BLUE}Info:${NC} RMSD value: $rmsd_value"
            TESTS_RUN=$((TESTS_RUN + 1))
            TESTS_PASSED=$((TESTS_PASSED + 1))
        fi
    fi
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; }

main() {
    test_header "RMSD Test"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
