#!/bin/bash
# Test: native GFN2 SCF preconditioner (-scf_precond frag), docs/SCF_PRECOND.md
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Oct 2026)
#
# (a) default (no switch): energy equals the reference of the binary before the preconditioner existed
# (b) -scf_precond frag: same energy as (a) within 1e-7 Eh and fewer SCF iterations
# (c) error paths: a charged system and a cluster of isolated ions ignore the switch with a warning
#     and give the same energy and iteration count as without it. The isolated ions (4 ions, 5.5 A apart, near-degenerate
#     frontier orbitals) sit at the stability limit of the SCF at 300 K: the run may converge or end with the clear
#     "SCF NOT converged" error and exit 1, but the switch must not change which of the two happens, and an energy is
#     printed only if the SCF converged (changed Oct 26, 2026; the former absolute reference depended on rounding).
# (d) -scf_precond frag, 5 fs NVE MD on the unrelaxed 197-atom cluster (the guard fires there): runs without
#     error, and the first MD energy (same geometry, same fixed point) equals the run without the switch within
#     1e-7 Eh. No iteration condition: frag does not improve the cold start on this structure (measured n=1,
#     -sp -threads 8: none 24, frag 33 iterations; docs/SCF_PRECOND.md).
# (i) -scf_precond_segment_atoms: a 58-atom PEO cut segmented (restricted additive Schwarz): same energy,
#     fewer iterations, same result at -threads 1 and 2.
# (j) -scf_precond_jacobian response on the structure of (i): same energy as fd, iterations <= fd + 1.
# (k) GFN1 with segmentation and response Jacobian on the same structure: same energy as without M, iterations <= none.
# (l) Hessian with a segmented M shared by the workers: frequencies within 1 cm-1 of no M, fewer iterations, one build.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 12: GFN2 SCF preconditioner (scf_precond frag)"
TEST_DIR="$SCRIPT_DIR"

# Tolerance 1e-7 Eh: the preconditioner changes only the path to the SCF fixed point; the energy may
# move by the SCF threshold (1e-5 on max|dq|), which is far below this. Output resolution is 1e-8 Eh.
ENERGY_TOLERANCE="0.0000001"

# Reference values: curcuma release build of commit 1632241c (before -scf_precond existed),
# -method gfn2 -threads 1, measured 2026-10-10 (docs/SCF_PRECOND.md).
REF_UREAWATER="-83.04383863"       # urea3-water8.gfn2-relaxed, 11 SCF iterations there
REF_HYDRONIUM="-15.31425657"       # hydronium-water2, -charge 1
REF_NACL="-9.49239348"             # na2cl2-isolated-ions

energy_of() { sed 's/\x1b\[[0-9;]*m//g' "$1" | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1; }
iters_of()  { sed 's/\x1b\[[0-9;]*m//g' "$1" | grep -oP 'SCF converged in \K[0-9]+' | tail -1; }

run() {   # run <log> <xyz> [args...]
    local LOG=$1; local XYZ=$2; shift 2
    $CURCUMA -sp "$XYZ" -method gfn2 -threads 1 -verbosity 1 -no_bmt "$@" > "$LOG" 2>&1
    assert_exit_code $? 0 "curcuma -sp $XYZ $*"
}

# Since the SimpleMD output rework (origin/master, Oct 2026) the MD prints no calculator lines ("Final Energy", "SCF precond: ...")
# any more: the first MD energy is Epot of the table row at t = 0, the cluster SCFs of the preconditioner build show in SCFTRACE.
md_first_epot() { sed 's/\x1b\[[0-9;]*m//g' "$1" | grep -E '^ +0\.000000 ' | head -1 | awk '{print $2}'; }
# SCFTRACE goes to stderr, the iteration table to stdout: a start can be glued behind an unterminated stdout line, so count
# every occurrence, not only line starts (Oct 26, 2026: the line-start count missed one and gave 38 instead of 39 extra starts).
count_scf_starts() { sed 's/\x1b\[[0-9;]*m//g' "$1" | grep -o 'SCFTRACE iter 0 ' | wc -l; }

first_energy() { sed 's/\x1b\[[0-9;]*m//g' "$1" | grep -oP 'Final Energy: \K[-0-9.]+' | head -1; }

run_md() {   # run_md <name> [args...]: 5 fs NVE GFN2 MD, stdout and stderr to <name>.log
    local NAME=$1; shift
    CURCUMA_SCF_TRACE=1 $CURCUMA -md ureawater197.xyz -method gfn2 -max_time 5 -time_step 1 -temperature 300 \
        -thermostat none -seed 42 -threads 1 -verbosity 2 -no_bmt "$@" > "$NAME.log" 2>&1
    assert_exit_code $? 0 "curcuma -md ureawater197.xyz $*"
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    # (a) default path
    run none.log ureawater.xyz
    local E0=$(energy_of none.log); local N0=$(iters_of none.log)
    assert_numeric_match "$REF_UREAWATER" "$E0" "$ENERGY_TOLERANCE" "default: energy vs pre-change reference"
    assert_string_not_in_file "SCF precond" none.log "default: no preconditioner set up"

    # (b) fragment preconditioner
    run frag.log ureawater.xyz -scf_precond frag
    local E1=$(energy_of frag.log); local N1=$(iters_of frag.log)
    assert_numeric_match "$E0" "$E1" "$ENERGY_TOLERANCE" "frag: energy equals the default run"
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$N0" ] && [ -n "$N1" ] && [ "$N1" -lt "$N0" ]; then
        echo -e "${GREEN}PASS${NC}: frag: fewer SCF iterations ($N1 < $N0)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: frag: SCF iterations not reduced (frag '$N1', default '$N0')"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # (c1) charged system: switch ignored with a warning
    run ion_none.log hydronium.xyz -charge 1
    run ion_frag.log hydronium.xyz -charge 1 -scf_precond frag
    assert_string_in_file "scf_precond 'frag' ignored: total charge" ion_frag.log "charged: warning printed"
    assert_numeric_match "$REF_HYDRONIUM" "$(energy_of ion_frag.log)" "$ENERGY_TOLERANCE" "charged: energy unchanged"
    assert_equals "$(iters_of ion_none.log)" "$(iters_of ion_frag.log)" "charged: same SCF iterations as without the switch"

    # (c2) isolated ions (single atoms): no block qualifies, switch ignored with a warning. The SCF of this structure is at
    #      the stability limit (converges at 1000 K in 38 iterations, not at 300 K on all builds): same outcome with and
    #      without the switch, never an energy without convergence.
    $CURCUMA -sp nacl.xyz -method gfn2 -threads 1 -verbosity 1 -no_bmt > nacl_none.log 2>&1; local RC_NONE=$?
    $CURCUMA -sp nacl.xyz -method gfn2 -threads 1 -verbosity 1 -no_bmt -scf_precond frag > nacl_frag.log 2>&1; local RC_FRAG=$?
    assert_string_in_file "no fragment qualifies" nacl_frag.log "isolated ions: warning printed"
    assert_equals "$RC_NONE" "$RC_FRAG" "isolated ions: same exit code with and without the switch"
    if [ "$RC_FRAG" -eq 0 ]; then
        assert_numeric_match "$(energy_of nacl_none.log)" "$(energy_of nacl_frag.log)" "$ENERGY_TOLERANCE" "isolated ions: converged, energy unchanged"
        echo "   (isolated ions converged; reference of the former build: $REF_NACL Eh)"
    else
        assert_string_in_file "SCF NOT converged" nacl_frag.log "isolated ions: not converged -> clear error"
        assert_equals "" "$(energy_of nacl_frag.log)" "isolated ions: no energy printed without convergence"
    fi

    # (d) MD on a structure where the guard fires: no error, same first energy
    run_md md_none
    run_md md_frag -scf_precond frag
    local N0=$(count_scf_starts md_none.log); local N1=$(count_scf_starts md_frag.log)
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "$((N1 - N0))" -eq 39 ]; then
        echo -e "${GREEN}PASS${NC}: MD frag: preconditioner built (39 fragment SCFs more than without: $N1 against $N0 SCF starts)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: MD frag: expected 39 more SCF starts than without the switch (none $N0, frag $N1)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    local EF0=$(md_first_epot md_none.log); local EF1=$(md_first_epot md_frag.log)
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$EF0" ] && [ -n "$EF1" ]; then
        echo -e "${GREEN}PASS${NC}: MD frag: first energies found ($EF0, $EF1)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: MD frag: first energy missing (none '$EF0', frag '$EF1')"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    assert_numeric_match "$EF0" "$EF1" "$ENERGY_TOLERANCE" "MD frag: first MD energy equals the run without the switch"

    # (e) CLI guard and spin/parity (Oct 2026): "-spin 1" as ONE argument is an error, not a silently
    #     ignored key; a radical without -spin gets a parity warning; with -spin 1 the energy equals
    #     xtb 6.7.1 (Cl2-, doublet: -9.321669164 Eh, identical with and without --uhf 1 in xtb).
    local REF_CL2="-9.32166916"
    $CURCUMA -sp cl2.xyz -method gfn2 -threads 1 -verbosity 1 -no_bmt -charge -1 "-spin 1" > spin_token.log 2>&1
    assert_exit_code $? 1 "flag token with whitespace (\"-spin 1\" as one argument) is refused"
    assert_string_in_file "contains whitespace" spin_token.log "flag token: error message printed"
    run cl2_nospin.log cl2.xyz -charge -1
    assert_string_in_file "inconsistent parity" cl2_nospin.log "radical without -spin: parity warning printed"
    run cl2_spin.log cl2.xyz -charge -1 -spin 1
    assert_string_not_in_file "inconsistent parity" cl2_spin.log "radical with -spin 1: no parity warning"
    assert_numeric_match "$REF_CL2" "$(energy_of cl2_spin.log)" "$ENERGY_TOLERANCE" "Cl2- with -spin 1 equals xtb 6.7.1"

    # (f) Int parameters given on the CLI arrive as double (diis_start 30 -> 30.0); the xtb lookups used to drop
    #     them silently (is_number_integer). With -scf_mode diis the damped warmup length changes the iteration count
    #     (urea3-water8: default 5 -> 13 iterations, diis_start 30 -> 21; diis_start 1 happens to give 13 as well).
    run diis_default.log ureawater.xyz -scf_mode diis
    run diis_start30.log ureawater.xyz -scf_mode diis -xtb.diis_start 30
    local ND0=$(iters_of diis_default.log); local ND1=$(iters_of diis_start30.log)
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$ND0" ] && [ -n "$ND1" ] && [ "$ND0" != "$ND1" ]; then
        echo -e "${GREEN}PASS${NC}: -xtb.diis_start 30 reaches the SCF (iterations $ND1 vs default $ND0)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: -xtb.diis_start 30 had no effect (iterations '$ND1', default '$ND0')"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    # The scf_ spelling (flat, reaches the SCF at -sp where the optimiser owns -diis_start) and the scf scope give the same iterations.
    run diis_start30_flat.log ureawater.xyz -scf_mode diis -scf_diis_start 30
    run diis_start30_scf.log ureawater.xyz -scf_mode diis -scf.diis_start 30
    assert_equals "$ND1" "$(iters_of diis_start30_flat.log)" "-scf_diis_start 30 (flat) gives the iterations of -xtb.diis_start 30"
    assert_equals "$ND1" "$(iters_of diis_start30_scf.log)" "-scf.diis_start 30 (alias) gives the iterations of -xtb.diis_start 30"
    # A longer SCF path ends at another point inside the convergence threshold (2e-6 Eh measured), hence 1e-5 here.
    assert_numeric_match "$(energy_of diis_default.log)" "$(energy_of diis_start30.log)" "0.00001" "diis_start moves the energy by less than the SCF threshold scale (1e-5 Eh)"

    # (g) scf_extrapolation_order and (h) scf_xlbomd_correctors on the CLI (double -> int, as in (f)): total SCF
    #     iterations of a 15 fs MD (CURCUMA_SCF_TRACE counts them, -threads 1 is deterministic). Before the
    #     fix every order gave the same total; the correctors are a MINIMUM cycle count, so 12 exceeds the
    #     natural ~7 cycles while 1 does not. Measured (urea3-water8, aspc): default 118, order 0 122, order 6 120,
    #     correctors 1: 126, correctors 12: 178. Relations only, no absolute numbers.
    md_total() {   # md_total <name> [args...]: sets MD_TOTAL = number of SCF iterations
        local NAME=$1; shift
        CURCUMA_SCF_TRACE=1 $CURCUMA -md ureawater.xyz -method gfn2 -max_time 15 -time_step 1 -temperature 300 \
            -thermostat none -seed 42 -threads 1 -verbosity 1 -no_bmt -scf_extrapolation aspc "$@" > "$NAME.log" 2>&1
        assert_exit_code $? 0 "curcuma -md ureawater.xyz -scf_extrapolation aspc $*"
        MD_TOTAL=$(grep -c '^SCFTRACE iter' "$NAME.log")
    }
    md_total xo_default
    local T_DEF=$MD_TOTAL
    md_total xo_0 -xtb.scf_extrapolation_order 0
    local T_O0=$MD_TOTAL
    md_total xo_6 -xtb.scf_extrapolation_order 6
    local T_O6=$MD_TOTAL
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "$T_O0" -gt 0 ] && [ "$T_O0" != "$T_DEF" ] && [ "$T_O6" != "$T_DEF" ] && [ "$T_O0" != "$T_O6" ]; then
        echo -e "${GREEN}PASS${NC}: -xtb.scf_extrapolation_order reaches the SCF (iterations default $T_DEF, order 0 $T_O0, order 6 $T_O6)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: -xtb.scf_extrapolation_order had no effect (default $T_DEF, order 0 $T_O0, order 6 $T_O6)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    md_total xl_c1 -scf_extrapolation_apply xlbomd -xtb.scf_xlbomd_correctors 1
    local T_C1=$MD_TOTAL
    md_total xl_c12 -scf_extrapolation_apply xlbomd -xtb.scf_xlbomd_correctors 12
    local T_C12=$MD_TOTAL
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "$T_C1" -gt 0 ] && [ "$T_C12" -ge $((T_C1 + 20)) ]; then
        echo -e "${GREEN}PASS${NC}: -xtb.scf_xlbomd_correctors reaches the SCF (iterations correctors 1: $T_C1, 12: $T_C12)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: -xtb.scf_xlbomd_correctors had no effect (correctors 1: $T_C1, 12: $T_C12)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # (i) segmentation of a large molecule (restricted additive Schwarz, Oct 2026): PEO cut of 58 atoms, forced
    #     into segments by -scf_precond_max_fragment_atoms 20, 3 heavy atoms per segment -> 8 segments. Same
    #     energy as without M (|dE| < 1e-6 Eh), fewer iterations, and the same energy and iteration count with
    #     -threads 1 (one result slot per segment job). Measured: none 11, segmented 5 iterations, -83.89363037 Eh.
    seg_run() {   # seg_run <log> <threads> [args...]
        local LOG=$1; local NT=$2; shift 2
        $CURCUMA -sp peo58.xyz -method gfn2 -threads "$NT" -verbosity 2 -no_bmt "$@" > "$LOG" 2>&1
        assert_exit_code $? 0 "curcuma -sp peo58.xyz -threads $NT $*"
    }
    seg_run seg_none.log 2
    seg_run seg_t2.log 2 -scf_precond frag -scf_precond_max_fragment_atoms 20 -scf_precond_segment_atoms 3
    seg_run seg_t1.log 1 -scf_precond frag -scf_precond_max_fragment_atoms 20 -scf_precond_segment_atoms 3
    assert_string_in_file "8 segments" seg_t2.log "segmented: 8 segments built"
    local ES0=$(energy_of seg_none.log); local NS0=$(iters_of seg_none.log)
    local ES2=$(energy_of seg_t2.log);   local NS2=$(iters_of seg_t2.log)
    assert_numeric_match "$ES0" "$ES2" "0.000001" "segmented: energy equals the run without M"
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$NS0" ] && [ -n "$NS2" ] && [ "$NS2" -lt "$NS0" ]; then
        echo -e "${GREEN}PASS${NC}: segmented: fewer SCF iterations ($NS2 < $NS0)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: segmented: SCF iterations not reduced (segmented '$NS2', none '$NS0')"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    assert_numeric_match "$ES2" "$(energy_of seg_t1.log)" "0.00000001" "segmented: -threads 1 gives the same energy as -threads 2"
    assert_equals "$NS2" "$(iters_of seg_t1.log)" "segmented: -threads 1 gives the same iteration count as -threads 2"

    # (j) -scf_precond_jacobian response (Oct 2026): the cluster Jacobians by linear response instead of central
    #     differences. Same structure and switches as (i): same energy as the fd run (< 1e-8 Eh) and at most one
    #     iteration more (measured: identical iteration count; |J_resp - J_fd| / |J_fd| ~ 1e-9).
    seg_run seg_resp.log 2 -scf_precond frag -scf_precond_max_fragment_atoms 20 -scf_precond_segment_atoms 3 -scf_precond_jacobian response
    assert_string_in_file "Jacobian response" seg_resp.log "response: Jacobian mode reported"
    local ER=$(energy_of seg_resp.log); local NR=$(iters_of seg_resp.log)
    assert_numeric_match "$ES2" "$ER" "0.00000001" "response: energy equals the fd run"
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$NR" ] && [ -n "$NS2" ] && [ "$NR" -le $((NS2 + 1)) ]; then
        echo -e "${GREEN}PASS${NC}: response: iterations $NR <= fd $NS2 + 1"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: response: iterations '$NR' > fd '$NS2' + 1"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # (k) GFN1 (Oct 2026): -scf_precond frag with segmentation and the response Jacobian on the structure of (i).
    #     GFN1 mixes only the shell charges. Same energy as without M (< 1e-6 Eh), iterations <= without M.
    g1_run() {   # g1_run <log> [args...]
        local LOG=$1; shift
        $CURCUMA -sp peo58.xyz -method gfn1 -threads 2 -verbosity 2 -no_bmt "$@" > "$LOG" 2>&1
        assert_exit_code $? 0 "curcuma -sp peo58.xyz -method gfn1 $*"
    }
    g1_run g1_none.log
    g1_run g1_seg.log -scf_precond frag -scf_precond_max_fragment_atoms 20 -scf_precond_segment_atoms 3 -scf_precond_jacobian response
    assert_string_in_file "8 segments" g1_seg.log "GFN1: 8 segments built"
    assert_string_not_in_file "ignored" g1_seg.log "GFN1: switch not ignored"
    local EG0=$(energy_of g1_none.log); local NG0=$(iters_of g1_none.log)
    local EG1=$(energy_of g1_seg.log);  local NG1=$(iters_of g1_seg.log)
    assert_numeric_match "$EG0" "$EG1" "0.000001" "GFN1: energy equals the run without M"
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$NG0" ] && [ -n "$NG1" ] && [ "$NG1" -le "$NG0" ]; then
        echo -e "${GREEN}PASS${NC}: GFN1: iterations $NG1 <= $NG0 without M"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: GFN1: iterations '$NG1' > '$NG0' without M"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # (l) Hessian with a segmented M (Oct 2026): the Hessian builds M once at the reference geometry and its
    #     worker threads share it (XTB::adoptSharedPreconditioner). methyl-a-d-glucopyranoside (27 atoms), forced
    #     into segments by -scf_precond_max_fragment_atoms 10: frequencies equal the run without M within 1 cm-1
    #     (SCF-threshold noise, measured 0.5), fewer SCF iterations, M built once (4 segment clusters, not 4 per worker).
    hess_run() {   # hess_run <name> [args...]
        local NAME=$1; shift
        CURCUMA_SCF_TRACE=1 $CURCUMA -hessian glc.xyz -method gfn2 -threads 4 -verbosity 1 -no_bmt "$@" > "$NAME.out" 2> "$NAME.err"
        assert_exit_code $? 0 "curcuma -hessian glc.xyz $*"
    }
    hess_run h_none
    hess_run h_seg -scf_precond frag -scf_precond_max_fragment_atoms 10 -scf_precond_segment_atoms 3
    local HN=$(awk '$1=="SCFTRACE" && $2=="iter" && $5==27' h_none.err | wc -l)
    local HS=$(awk '$1=="SCFTRACE" && $2=="iter" && $5==27' h_seg.err | wc -l)
    local HB=$(awk '$1=="SCFTRACE" && $2=="iter" && $3==0 && $5!=27' h_seg.err | wc -l)
    local FDEV=$(python3 - h_none.out h_seg.out <<'PY'
import re, sys
def freqs(f):
    t = re.sub(r'\x1b\[[0-9;]*m', '', open(f, errors='replace').read())
    i = t.index('Vibrational Frequencies'); j = t.index('Finished after', i)
    return [float(m.group(1)) for m in re.finditer(r'([0-9.]+)(\(i\)|\(\*\))?', t[i + 30:j])]
a, b = freqs(sys.argv[1]), freqs(sys.argv[2])
print('%.3f' % max(abs(x - y) for x, y in zip(a, b)) if a and len(a) == len(b) else 'mismatch')
PY
)
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "$FDEV" != "mismatch" ] && python3 -c "import sys; sys.exit(0 if float('$FDEV') <= 1.0 else 1)" 2>/dev/null; then
        echo -e "${GREEN}PASS${NC}: Hessian: frequencies equal the run without M (max deviation $FDEV cm-1)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: Hessian: frequency deviation '$FDEV' cm-1"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "$HS" -gt 0 ] && [ "$HS" -lt "$HN" ] && [ "$HB" -eq 4 ]; then
        echo -e "${GREEN}PASS${NC}: Hessian: SCF iterations $HS < $HN, M built once ($HB cluster SCFs)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: Hessian: iterations $HS vs $HN, cluster SCFs $HB (expected 4)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
