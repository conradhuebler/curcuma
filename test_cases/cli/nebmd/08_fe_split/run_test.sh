#!/bin/bash
# NEB-MD CLI Test 08 - decomposition of the constrained mean force in fe.csv
# Copyright (C) 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
# Claude Generated 2026 (AI-generated, machine-tested)
#
# Ethane, GFN-FF, scheme 'independent' + umbrella + constrain, 6 windows, 2 ps per window, with a
# z wall at 0.005 A^2 (max z of the interior windows ~0.006 A^2, so the wall acts in a few % of the
# samples). Checks: exit 0; 4 estimator-I rows (the endpoint images are not sampled);
# fe.csv row header carries dGds_gradE_Eh,dGds_metric_Eh,dGds_wall_Eh; for every estimator-I row
# gradE + metric + wall = dGds_Eh within the 8-decimal rounding of dGds_Eh (6e-9 Eh); the wall column
# is non-zero in at least one row (active wall); every value finite.
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

pass() { TESTS_RUN=$((TESTS_RUN + 1)); TESTS_PASSED=$((TESTS_PASSED + 1)); echo -e "${GREEN}ok PASS${NC}: $1"; }
fail() { TESTS_RUN=$((TESTS_RUN + 1)); TESTS_FAILED=$((TESTS_FAILED + 1)); echo -e "${RED}xx FAIL${NC}: $1"; }

run_test() {
    cd "$TEST_DIR"
    $CURCUMA -nebmd start.xyz end.xyz \
        -nebmd.scheme independent -nebmd.umbrella true -nebmd.constrain true \
        -nebmd.nimages 6 -nebmd.z_wall 0.005 \
        -simplemd.method gfnff -simplemd.time_step 0.5 -simplemd.temperature 300 \
        -simplemd.thermostat csvr -simplemd.max_time 2000 -simplemd.seed 1 \
        > stdout.log 2> stderr.log
    assert_exit_code $? 0 "NEB-MD with z wall should succeed"
    return 0
}

validate_results() {
    cd "$TEST_DIR"
    local fe=$(find_output_file "start.neb.fe.csv")
    if [ -z "$fe" ] || [ ! -f "$fe" ]; then
        fail "start.neb.fe.csv not written"
        return 1
    fi
    grep -q "^estimator,s,dGds_Eh,G_kcal,samples,rejected_frac,dGds_gradE_Eh,dGds_metric_Eh,dGds_wall_Eh$" "$fe" \
        && pass "fe.csv row header with the three split columns" || fail "fe.csv row header without the split columns"
    local res=$(python3 - "$fe" <<'PYEOF'
import math, sys
rows = [l.strip().split(",") for l in open(sys.argv[1]) if l.startswith("I,")]
dmax, wmax, fin = 0.0, 0.0, 1
for r in rows:
    t, g, m, w = float(r[2]), float(r[6]), float(r[7]), float(r[8])
    fin &= all(math.isfinite(x) for x in (t, g, m, w))
    dmax = max(dmax, abs(g + m + w - t))
    wmax = max(wmax, abs(w))
print(len(rows), fin, "%.3e" % dmax, "%.3e" % wmax)
PYEOF
)
    echo "  rows / finite / max|gradE+metric+wall-dGds| / max|wall|: $res"
    local n=$(echo "$res" | awk '{print $1}')
    [ "$n" = "4" ] && pass "4 estimator-I rows (6 images, fixed endpoints)" || fail "$n estimator-I rows, expected 4"
    [ "$(echo "$res" | awk '{print $2}')" = "1" ] && pass "split columns finite" || fail "non-finite split value"
    python3 -c "import sys; sys.exit(0 if float('$(echo "$res" | awk '{print $3}')') <= 6e-9 else 1)" \
        && pass "gradE + metric + wall = dGds_Eh (<= 6e-9 Eh, 8-decimal rounding)" || fail "split does not add up to dGds_Eh"
    python3 -c "import sys; sys.exit(0 if float('$(echo "$res" | awk '{print $4}')') > 0.0 else 1)" \
        && pass "wall column non-zero (wall active)" || fail "wall column zero although z_wall < <z>"
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; rm -f start.neb.*; rm -rf start.nebmd.*; }

main() {
    test_header "NEB-MD fe.csv mean-force split (gfnff, ethane, z wall 0.005)"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
