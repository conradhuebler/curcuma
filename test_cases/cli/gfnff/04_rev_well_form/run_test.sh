#!/bin/bash

# Test: rev-gfnff stage 3a(iii) bond-well forms - identity of the default, liveness of the flag
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 18, 2026)
#
# -gfnff.rev_well_form selects the shape of the bond well: gauss (delivered), mg, erfmorse
# (docs/REV_GFNFF_STAGE3A.md). Three assertions, and the third is the one that keeps the other
# two honest:
#   1. IDENTITY - the DEFAULT and an explicit `mg` give the same energy to 12 digits, and an
#      explicit `gauss` still gives the recorded delivered-Gaussian value. A default flip must
#      not disturb the form it flipped away from.
#   2. gfnff is untouched - plain -method gfnff still gives its recorded 12-digit energy. The
#      whole rev path is gated on rev_enabled, so a well-form default can never reach it.
#   3. LIVENESS - gauss and erfmorse must each differ from the default by more than 1e-6 Eh, and
#      must differ from each other. Without this the test would pass on a build where the flag is
#      silently ignored, which is exactly how a switchable form dies.
#
# UPDATED Sep 19, 2026: the operator made `mg` the default (test_cases/revgfnff/_log/WORK_STATUS.md
# package 6b), so the roles of "default" and "gauss" are swapped relative to the Sep 18 version of
# this test. The PINNED value stays on `gauss`, the stable delivered form:
#   revgfnff gauss  -4.673521653477 Eh   gfnff  -4.672737068614 Eh
# The two new forms are AI-fitted per element pair (rev_well_table.h) and deliberately NOT
# pinned to a value here: their numbers will move when the table is refitted, while the identity
# and liveness statements above must survive that. Measured Sep 19 (binary 916847ff):
#   default = mg -4.546943898047, erfmorse -4.539136588545 (mg - gauss +0.126578 Eh).

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../test_utils.sh"

TEST_NAME="gfnff - 04: rev-gfnff bond-well form identity and liveness"
TEST_DIR="$SCRIPT_DIR"

REV_GAUSS_EH=-4.673521653477   # the DELIVERED Gaussian; the default is mg since Sep 19, 2026
GFNFF_EH=-4.672737068614
IDENTITY_TOL=1e-11
LIVENESS_MIN=1e-6

sp() {   # method topology extra... -> energy on stdout, 12+ digits
    # the printed "Single Point Energy" carries 8 decimals, which cannot express a 1e-11
    # identity statement; the -batch JSONL carries the full double.
    local method="$1" topo="$2"; shift 2
    local d
    d=$(mktemp -d)
    cp input.xyz "$d/s.xyz"
    ( cd "$d" && timeout 300 $CURCUMA -sp s.xyz -method "$method" -gfnff.topology_mode "$topo" \
        -gfnff.cache_topology false -no_bmt -threads 1 -verbosity 1 \
        -batch true -batch_out o.jsonl "$@" > o.log 2>&1 )
    python3 -c "import json,sys; print('%.12f' % json.loads(open('$d/o.jsonl').read().splitlines()[0])['energy_eh'])"
    rm -rf "$d"
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"

    local e_def e_gauss e_mg e_em e_gfnff
    e_def=$(sp revgfnff react)
    e_gauss=$(sp revgfnff react -gfnff.rev_well_form gauss)
    e_mg=$(sp revgfnff react -gfnff.rev_well_form mg)
    e_em=$(sp revgfnff react -gfnff.rev_well_form erfmorse)
    e_gfnff=$(sp gfnff static)

    echo "default   $e_def"
    echo "gauss     $e_gauss"
    echo "mg        $e_mg"
    echo "erfmorse  $e_em"
    echo "gfnff     $e_gfnff"

    TESTS_RUN=$((TESTS_RUN + 1))
    local py_rc
    set +e
    python3 - "$e_def" "$e_gauss" "$e_mg" "$e_em" "$e_gfnff" \
             "$REV_GAUSS_EH" "$GFNFF_EH" "$IDENTITY_TOL" "$LIVENESS_MIN" <<'PYEOF'
import sys
d, g, mg, em, ff, ref_gauss, ref_ff, tol, live = (float(x) for x in sys.argv[1:10])
ok, reasons = True, []
# the DEFAULT is mg since Sep 19, 2026 - assert that, and keep the pinned value on gauss
if abs(d - mg) > tol:
    ok = False; reasons.append(f"default {d:.12f} != explicit mg {mg:.12f} - mg is the default")
if abs(g - ref_gauss) > tol:
    ok = False; reasons.append(f"explicit gauss {g:.12f} != recorded {ref_gauss:.12f}")
if abs(ff - ref_ff) > tol:
    ok = False; reasons.append(f"gfnff {ff:.12f} != recorded {ref_ff:.12f}")
if abs(g - d) < live:
    ok = False; reasons.append(f"gauss differs from the default by only {abs(g-d):.2e} Eh - flag ignored?")
if abs(em - d) < live:
    ok = False; reasons.append(f"erfmorse differs from the default by only {abs(em-d):.2e} Eh - flag ignored?")
if abs(em - mg) < live:
    ok = False; reasons.append(f"erfmorse and mg differ by only {abs(em-mg):.2e} Eh")
print(f"mg - gauss = {(mg-g):+.6f} Eh, erfmorse - gauss = {(em-g):+.6f} Eh, "
      f"erfmorse - mg = {(em-mg):+.2e} Eh")
if not ok:
    print("FAIL reasons: " + "; ".join(reasons))
sys.exit(0 if ok else 1)
PYEOF
    py_rc=$?
    set -e
    if [ $py_rc -eq 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: default == mg, gauss pinned, gfnff untouched, all three forms live"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: well-form identity/liveness"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
