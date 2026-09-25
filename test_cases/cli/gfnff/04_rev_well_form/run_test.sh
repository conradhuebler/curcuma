#!/bin/bash

# Test: rev-gfnff stage 3a(iii) bond-well forms - identity of the default, liveness of the flag
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 18, 2026)
#
# -gfnff.rev_well_form selects the shape of the bond well: gauss (delivered), mg, erfmorse, mg2,
# mg3 (docs/REV_GFNFF_STAGE3A.md). Three assertions, and the third is the one that keeps the other
# two honest:
#   1. IDENTITY - the DEFAULT and an explicit `mg3` give the same energy to 12 digits, and an
#      explicit `gauss` still gives the recorded delivered-Gaussian value. A default flip must
#      not disturb the forms it flipped away from.
#   2. gfnff is untouched - plain -method gfnff still gives its recorded 12-digit energy. The
#      whole rev path is gated on rev_enabled, so a well-form default can never reach it.
#   3. LIVENESS - gauss, mg, erfmorse and mg2 must each differ from the default by more than
#      1e-6 Eh, and all five forms must be pairwise distinct. Without this the test would pass on
#      a build where the flag is silently ignored, which is exactly how a switchable form dies.
#
# UPDATED Sep 19, 2026: the operator made `mg` the default (test_cases/revgfnff/_log/WORK_STATUS.md
# package 6b), so the roles of "default" and "gauss" were swapped relative to the Sep 18 version.
# UPDATED Sep 22, 2026: the operator made `mg3` the default (WORK_STATUS package 12), so the
# identity clause now names `mg3` and the liveness clause was widened from two forms to four.
# No numeric threshold was changed: IDENTITY_TOL stays 1e-11 and LIVENESS_MIN stays 1e-6.
# The PINNED value stays on `gauss`, the stable delivered form:
#   revgfnff gauss  -4.673521653477 Eh   gfnff  -4.672737068614 Eh
# The four fitted forms are AI-fitted per element pair (rev_well_table.h / rev_well_table_v2.h)
# and deliberately NOT pinned to a value here: their numbers will move when the table is refitted,
# while the identity and liveness statements above must survive that. Measured Sep 22, 2026:
#   default = mg3 -4.789915106585, mg -4.546943898047, erfmorse -4.539136588545,
#   mg2 -4.550323010840 (mg3 - gauss -0.116393 Eh). The closest pair is mg/mg2 at 3.4e-3 Eh,
#   i.e. 3400x the 1e-6 liveness bound, so the widened clause is not a knife edge.
# RE-PINNED Sep 25, 2026 (multi-gpu merge, operator decision): feature/multi-gpu switched the
# reference-less bonded-triple ATM dispersion term off by default (-gfnff.dispersion_atm false).
# Both pins moved by +1.2195e-8 Eh, nothing else: gauss -4.673521653477 -> -4.673521641282,
# gfnff -4.672737068614 -> -4.672737056419; with -gfnff.dispersion_atm true the merged binary
# reproduces the old pins exactly. Thresholds unchanged.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../test_utils.sh"

TEST_NAME="gfnff - 04: rev-gfnff bond-well form identity and liveness"
TEST_DIR="$SCRIPT_DIR"

REV_GAUSS_EH=-4.673521641282   # the DELIVERED Gaussian; the default is mg3 since Sep 22, 2026
GFNFF_EH=-4.672737056419
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

    local e_def e_gauss e_mg e_em e_mg2 e_mg3 e_gfnff
    e_def=$(sp revgfnff react)
    e_gauss=$(sp revgfnff react -gfnff.rev_well_form gauss)
    e_mg=$(sp revgfnff react -gfnff.rev_well_form mg)
    e_em=$(sp revgfnff react -gfnff.rev_well_form erfmorse)
    e_mg2=$(sp revgfnff react -gfnff.rev_well_form mg2)
    e_mg3=$(sp revgfnff react -gfnff.rev_well_form mg3)
    e_gfnff=$(sp gfnff static)

    echo "default   $e_def"
    echo "gauss     $e_gauss"
    echo "mg        $e_mg"
    echo "erfmorse  $e_em"
    echo "mg2       $e_mg2"
    echo "mg3       $e_mg3"
    echo "gfnff     $e_gfnff"

    TESTS_RUN=$((TESTS_RUN + 1))
    local py_rc
    set +e
    python3 - "$e_def" "$e_gauss" "$e_mg" "$e_em" "$e_mg2" "$e_mg3" "$e_gfnff" \
             "$REV_GAUSS_EH" "$GFNFF_EH" "$IDENTITY_TOL" "$LIVENESS_MIN" <<'PYEOF'
import sys
d, g, mg, em, mg2, mg3, ff, ref_gauss, ref_ff, tol, live = (float(x) for x in sys.argv[1:12])
ok, reasons = True, []
# the DEFAULT is mg3 since Sep 22, 2026 - assert that, and keep the pinned value on gauss
if abs(d - mg3) > tol:
    ok = False; reasons.append(f"default {d:.12f} != explicit mg3 {mg3:.12f} - mg3 is the default")
if abs(g - ref_gauss) > tol:
    ok = False; reasons.append(f"explicit gauss {g:.12f} != recorded {ref_gauss:.12f}")
if abs(ff - ref_ff) > tol:
    ok = False; reasons.append(f"gfnff {ff:.12f} != recorded {ref_ff:.12f}")
# LIVENESS: every non-default form must move the energy, and all five must be pairwise distinct
for name, val in (("gauss", g), ("mg", mg), ("erfmorse", em), ("mg2", mg2)):
    if abs(val - d) < live:
        ok = False
        reasons.append(f"{name} differs from the default by only {abs(val-d):.2e} Eh - flag ignored?")
forms = [("gauss", g), ("mg", mg), ("erfmorse", em), ("mg2", mg2), ("mg3", mg3)]
for i in range(len(forms)):
    for j in range(i + 1, len(forms)):
        if abs(forms[i][1] - forms[j][1]) < live:
            ok = False
            reasons.append(f"{forms[i][0]} and {forms[j][0]} differ by only "
                           f"{abs(forms[i][1]-forms[j][1]):.2e} Eh")
print(f"mg - gauss = {(mg-g):+.6f} Eh, erfmorse - gauss = {(em-g):+.6f} Eh, "
      f"mg2 - gauss = {(mg2-g):+.6f} Eh, mg3 - gauss = {(mg3-g):+.6f} Eh, "
      f"closest pair = {min(abs(forms[i][1]-forms[j][1]) for i in range(5) for j in range(i+1, 5)):.2e} Eh")
if not ok:
    print("FAIL reasons: " + "; ".join(reasons))
sys.exit(0 if ok else 1)
PYEOF
    py_rc=$?
    set -e
    if [ $py_rc -eq 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: default == mg3, gauss pinned, gfnff untouched, all five forms live and distinct"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: well-form identity/liveness"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
