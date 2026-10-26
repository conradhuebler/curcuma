#!/bin/bash
# Test: shared SCF flags reach the native QM engine (hf, hf-3c) with master's precedence
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Oct 27, 2026). AI-generated, machine-tested only.
#
# Native QM reads registry module "qm" and the shared module "scf" (scf_parameters.h). Before Oct 27, 2026 the flat
# -scf_threshold / -scf_mode stayed in the command module with a "registered in multiple modules" warning and never
# reached the SCF, -scf.* was not read at all, and an unconverged HF SCF returned its energy with exit code 0.
#   (a) flat -scf_threshold, -scf.scf_threshold and -qm.scf_threshold give the same (larger) iteration count, no warning
#   (b) precedence: -qm.<name> wins over -scf.<name>
#   (c) unconverged SCF is an error; -scf_allow_unconverged true accepts it with a warning
#   (d) an SCF mode of another solver (broyden) is an error for HF
#   (e) the old spelling -qm.diis_start equals -scf_diis_start and -qm.scf_diis_start (renamed Oct 27, 2026)
#   (f) at -sp the flat -max_iterations stays with the optimiser (module opt owns it) and does not stop the SCF
#   (g) hf-3c reads the same flags; all converged energies agree

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="routing - 04: shared SCF flags for native QM (hf, hf-3c)"
TEST_DIR="$SCRIPT_DIR"

strip()     { sed 's/\x1b\[[0-9;]*m//g' "$1"; }
energy_of() { strip "$1" | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1; }
iters_of()  { grep -c 'SCF iter' "$1"; }
exit_of()   { grep -oP 'EXIT \K[0-9]+$' "$1" | tail -1; }
run() { local LOG=$1; local M=$2; shift 2; $CURCUMA -sp water.xyz -method "$M" -threads 1 -verbosity 1 -no_bmt "$@" > "$LOG" 2>&1; echo "EXIT $?" >> "$LOG"; }
check() {  # check <condition result 0/1> <message>
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "$1" -eq 0 ]; then echo -e "${GREEN}PASS${NC}: $2"; TESTS_PASSED=$((TESTS_PASSED + 1));
    else echo -e "${RED}FAIL${NC}: $2"; TESTS_FAILED=$((TESTS_FAILED + 1)); fi
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    run base.log hf
    local N0=$(iters_of base.log); local E0=$(energy_of base.log)

    # (a) three spellings, same effect
    run flat.log hf -scf_threshold 1e-10
    run scf.log hf -scf.scf_threshold 1e-10
    run qm.log hf -qm.scf_threshold 1e-10
    local N1=$(iters_of flat.log)
    [ "$N1" -gt "$N0" ]; check $? "flat -scf_threshold 1e-10 reaches the SCF ($N0 -> $N1 iterations)"
    assert_equals "$N1" "$(iters_of scf.log)" "-scf.scf_threshold gives the iterations of the flat flag"
    assert_equals "$N1" "$(iters_of qm.log)" "-qm.scf_threshold gives the iterations of the flat flag"
    assert_string_not_in_file "registered in multiple modules" flat.log "flat -scf_threshold is routed to qm without the ambiguity warning"

    # (b) the method scope wins
    run prec.log hf -scf.scf_threshold 1e-10 -qm.scf_threshold 1e-4
    assert_equals "$N0" "$(iters_of prec.log)" "-qm.scf_threshold 1e-4 wins over -scf.scf_threshold 1e-10"

    # (c) unconverged SCF
    run unconv.log hf -scf_max_iterations 2
    [ "$(exit_of unconv.log)" != "0" ]; check $? "unconverged HF SCF ends with an error (exit $(exit_of unconv.log))"
    assert_string_in_file "scf_allow_unconverged" unconv.log "the error names -scf_allow_unconverged"
    run allow.log hf -scf_max_iterations 2 -scf_allow_unconverged true
    assert_equals "0" "$(exit_of allow.log)" "-scf_allow_unconverged true accepts the unconverged energy"

    # (d) mode of another solver
    run broyden.log hf -scf_mode broyden
    [ "$(exit_of broyden.log)" != "0" ]; check $? "-scf_mode broyden is an error for HF"
    assert_string_in_file "not available for native HF" broyden.log "the error says which modes HF takes"

    # (e) renamed DIIS start, old spelling still an alias
    run ds_old.log hf -qm.diis_start 6
    run ds_flat.log hf -scf_diis_start 6
    run ds_qm.log hf -qm.scf_diis_start 6
    local ND=$(iters_of ds_old.log)
    [ "$ND" -ne "$N0" ]; check $? "-qm.diis_start 6 changes the iterations ($N0 -> $ND)"
    assert_equals "$ND" "$(iters_of ds_flat.log)" "-scf_diis_start 6 equals the old spelling -qm.diis_start 6"
    assert_equals "$ND" "$(iters_of ds_qm.log)" "-qm.scf_diis_start 6 equals the old spelling -qm.diis_start 6"

    # (f) -max_iterations belongs to the optimiser at -sp
    run maxit.log hf -max_iterations 2
    assert_equals "$N0" "$(iters_of maxit.log)" "flat -max_iterations at -sp does not reach the SCF (module opt owns it)"

    # (g) hf-3c reads the same flags
    run h3_base.log hf-3c
    run h3_flat.log hf-3c -scf_threshold 1e-10
    [ "$(iters_of h3_flat.log)" -gt "$(iters_of h3_base.log)" ]; check $? "hf-3c: flat -scf_threshold reaches the SCF ($(iters_of h3_base.log) -> $(iters_of h3_flat.log))"
    assert_numeric_match "$(energy_of h3_base.log)" "$(energy_of h3_flat.log)" "0.00000001" "hf-3c: tight SCF gives the same energy"

    for L in flat scf qm prec ds_old ds_flat ds_qm maxit; do
        assert_numeric_match "$E0" "$(energy_of $L.log)" "0.00000001" "hf energy of $L.log equals the default run"
    done

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
