#!/bin/bash
# Test: per-method help, method overview, unread-parameter warning, flat routing of a -qm parameter
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Oct 31, 2026). AI-generated, machine-tested only.
#
# Help and parameters are organised per method (MethodCatalog, docs/METHODS_OVERVIEW.md).
#   (a) -help methods: exit 0, every MethodFactory name appears
#   (b) -help gfn2|gfnff|b2plyp|hf-3c|uff: mandatory sections; every parameter listed is registered in its module
#   (c) -help <capability> and -help <category> still work; the orca collision is reported; unknown names say so
#   (d) an unknown -qm.<name> and a parameter the method does not read give a warning, the run still succeeds
#   (e) flat -grid sg2 -method pbe equals -qm.grid sg2 (and differs from the default grid)
#   (f) a parameter the method reads gives no warning

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="routing - 05: per-method help and unread-parameter warning"
TEST_DIR="$SCRIPT_DIR"

strip()     { sed 's/\x1b\[[0-9;]*m//g' "$1"; }
energy_of() { strip "$1" | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1; }
exit_of()   { grep -oP 'EXIT \K[0-9]+$' "$1" | tail -1; }
run() { local LOG=$1; shift; $CURCUMA -sp water.xyz -threads 1 -verbosity 1 -no_bmt "$@" > "$LOG" 2>&1; echo "EXIT $?" >> "$LOG"; }
check() {  # check <condition result 0/1> <message>
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "$1" -eq 0 ]; then echo -e "${GREEN}PASS${NC}: $2"; TESTS_PASSED=$((TESTS_PASSED + 1));
    else echo -e "${RED}FAIL${NC}: $2"; TESTS_FAILED=$((TESTS_FAILED + 1)); fi
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    # (a) overview lists every factory name
    $CURCUMA -help methods > methods.log 2>&1; local RC=$?
    assert_equals "0" "$RC" "-help methods exits with 0"
    $CURCUMA -check_method_catalog routes 2>&1 | grep ' scope=' | awk '{print $1}' > routes.log
    local MISSING=""
    while read -r N; do
        grep -qE -- "(^|[ ,])$(printf '%s' "$N" | sed 's/[+.]/\\&/g')(,|$)" methods.log || MISSING="$MISSING $N"
    done < routes.log
    [ -z "$MISSING" ] && [ "$(wc -l < routes.log)" -gt 80 ]; check $? "-help methods names every MethodFactory name ($(wc -l < routes.log) names; missing:${MISSING:- none})"

    # (b) per-method help: sections and registered parameters
    for M in gfn2 gfnff b2plyp hf-3c uff; do
        $CURCUMA -help "$M" > "help_$M.log" 2>&1
        local OK=0
        for S in "Method: $M" "Summary:" "Implementation:" "Status:" "Gradient:" "GPU:" "Restrictions:" \
                 "Variants and aliases" "Parameters this method reads" "Shared SCF knobs"; do
            grep -qF -- "$S" "help_$M.log" || { OK=1; echo "  missing section '$S' in -help $M"; }
        done
        check $OK "-help $M has all mandatory sections"
        # every "(-module.name)" or leading "-module.name" in the parameter list must be registered
        local BAD="" COUNT=0
        while read -r MOD NAME; do
            COUNT=$((COUNT + 1))
            [ -f "mod_$MOD.log" ] || $CURCUMA -help_module "$MOD" > "mod_$MOD.log" 2>&1
            grep -qE "^  -$NAME <" "mod_$MOD.log" || BAD="$BAD $MOD.$NAME"
        done < <(grep -oP '^    -\K[a-z0-9_]+\.[a-z0-9_]+(?=  )|^    -[a-z0-9_]+  \(-\K[a-z0-9_]+\.[a-z0-9_]+' "help_$M.log" | tr '.' ' ')
        [ -z "$BAD" ] && [ "$COUNT" -gt 0 ]; check $? "-help $M: all $COUNT listed parameters are registered (unregistered:${BAD:- none})"
    done
    assert_string_in_file "scf_threshold" help_gfn2.log "-help gfn2 lists the shared -scf_threshold"
    assert_string_in_file "dh_pt2" help_b2plyp.log "-help b2plyp lists -qm.dh_pt2"
    assert_string_not_in_file "vwn_variant" help_b2plyp.log "-help b2plyp does not list vwn_variant (not read)"

    # (c) capability, category, collision, unknown
    $CURCUMA -help sp > help_sp.log 2>&1
    assert_string_in_file "Description:" help_sp.log "-help sp shows the capability"
    $CURCUMA -help analysis > help_analysis.log 2>&1
    assert_string_in_file "Analysis Capabilities" help_analysis.log "-help analysis still shows the category"
    $CURCUMA -help orca > help_orca.log 2>&1
    assert_string_in_file "also a method name" help_orca.log "-help orca reports the capability/method collision"
    $CURCUMA -help method:orca > help_method_orca.log 2>&1
    assert_string_in_file "Method: orca" help_method_orca.log "-help method:orca shows the method"
    $CURCUMA -help no_such_thing > help_unknown.log 2>&1
    assert_string_in_file "Unknown category" help_unknown.log "-help of an unknown name says so"

    # (d) warnings, run continues
    run unknown.log -method hf -qm.no_such_param 3
    assert_string_in_file "not a registered parameter of module qm" unknown.log "unknown -qm.no_such_param is reported"
    assert_equals "0" "$(exit_of unknown.log)" "the run with an unknown parameter still succeeds"
    run unread.log -method hf -grid sg2
    assert_string_in_file "is not read by method hf" unread.log "-grid with -method hf is reported as not read"
    assert_equals "0" "$(exit_of unread.log)" "the run with an unread parameter still succeeds"
    run unread_ff.log -method gfnff -scf_guess h0
    assert_string_in_file "not read by method gfnff" unread_ff.log "-scf_guess with -method gfnff is reported as not read"

    # (e) flat routing of a -qm parameter
    run pbe_base.log -method pbe
    run pbe_flat.log -method pbe -grid sg2
    run pbe_dot.log -method pbe -qm.grid sg2
    local EB=$(energy_of pbe_base.log); local EF=$(energy_of pbe_flat.log); local ED=$(energy_of pbe_dot.log)
    assert_numeric_match "$ED" "$EF" "0.0000000001" "flat -grid sg2 equals -qm.grid sg2 for pbe"
    python3 -c "import sys; sys.exit(0 if abs(float('$EF') - float('$EB')) > 1e-7 else 1)"; check $? "sg2 changes the pbe energy ($EB -> $EF)"

    # (f) no warning for a parameter the method reads
    assert_string_not_in_file "not read by method" pbe_flat.log "-grid with -method pbe gives no unread warning"
    run gfn2_guess.log -method gfn2 -scf_guess h0
    assert_string_not_in_file "not read by method" gfn2_guess.log "-scf_guess with -method gfn2 gives no unread warning"

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
