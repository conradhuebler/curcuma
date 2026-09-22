#!/bin/bash

# Test: the class-C radical-adduct falsifier of rev-gfnff stage 3a(ii)
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 18, 2026)
#
# A free hydrogen radical approaching a SATURATED carbon must be repelled: CH4 + H has no bond to
# make. The delivered valence share is a LEFT-OVER rule, f_i = clip((Val_i - sum_{k != j} w_ik) /
# w_ij), and with a term weight of ~1 on every partner one partner too many makes every pair of
# that atom see "nothing left" - so the carbon hands out 0 of its 4 valences and the growing
# budget papers over it, producing an artificial adduct far below the reference
# (FABLE_REVIEW_2 A.4, test_cases/revgfnff/_log/WORK_STATUS.md packages 1 and 3).
#
# The reference is the class-C r2SCAN-3c scan test_cases/revgfnff/ref/C/ch4_H (15 rigid points,
# d = 1.00 to 3.00 A). Both curves are referred to their own d = 3.00 A point, so the comparison
# is of the approach PROFILE and no absolute energy zero enters.
#
# Two arms, and the second is what makes this a falsifier rather than a golden value:
#   1. the DEFAULT (the valence-conserving share since Sep 19, 2026) must keep the whole curve
#      above -10 kcal/mol - measured min -1.5, at 1.0/1.2/1.3 A: +16.8 / +23.9 / +22.4. This is
#      the acceptance criterion for the conserving share; if it stops holding, that mode has
#      regressed.
#   2. -gfnff.rev_share_form delivered must still show the old, bad profile - measured min
#      deviation -90.2 kcal/mol - within 5 kcal/mol. That is a regression pin on the OLD
#      behaviour, kept so that the arm cannot silently become the good one and hide a
#      dispatch bug that routes both arms to the same code.
#
# UPDATED Sep 19, 2026 (WORK_STATUS package 6): the arms were swapped when `conserving` became
# the default, and the delivered pin moved -87.0 -> -89.4 because `mg` became the default WELL
# form in the same session (a deeper, wider well makes the artificial adduct slightly deeper).
# UPDATED Sep 22, 2026 (WORK_STATUS package 12): `mg3` became the default well form, so the
# delivered arm - which pins no well form of its own - is now evaluated with the mg3 well and its
# min deviation moved -89.4 -> -90.2 kcal/mol. Same cause as the Sep 19 move, and again a
# RE-MEASUREMENT of the old-bad behaviour, not a relaxation: the pin tracks the arm, the arm is
# still the bad one, and the conserving default arm still comes out at -1.53 (floor -10).
# The tolerance is unchanged at 5 kcal/mol in both updates.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../test_utils.sh"

TEST_NAME="gfnff - 03: rev-gfnff class-C radical-adduct falsifier (CH4 + H)"
TEST_DIR="$SCRIPT_DIR"

DELIVERED_MIN=-90.2      # kcal/mol, measured Sep 22, 2026 (mg3 well + delivered share; was -89.4 under mg)
DELIVERED_TOL=5.0
CONSERVING_FLOOR=-10.0

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    TESTS_RUN=$((TESTS_RUN + 1))
    local py_out py_rc
    set +e
    py_out=$(python3 - "$CURCUMA" "$PROJECT_ROOT" "$DELIVERED_MIN" "$DELIVERED_TOL" \
                      "$CONSERVING_FLOOR" <<'PYEOF'
import json, math, os, shutil, subprocess, sys, tempfile
from pathlib import Path

binary, root, delivered_min, delivered_tol, cons_floor = sys.argv[1], sys.argv[2], *map(float, sys.argv[3:6])
REF = Path(root) / "test_cases" / "revgfnff" / "ref" / "C" / "ch4_H"
H2K = 627.5094740631

info = json.loads((REF / "energies.json").read_text())
lines = (REF / "points.xyz").read_text().splitlines()
frames, i = [], 0
while i < len(lines):
    if not lines[i].strip():
        i += 1
        continue
    n = int(lines[i].split()[0])
    frames.append(lines[i:i + 2 + n])
    i += n + 2


def sp(frame, extra):
    d = tempfile.mkdtemp()
    try:
        (Path(d) / "s.xyz").write_text("\n".join(frame) + "\n")
        cmd = [binary, "-sp", "s.xyz", "-method", "revgfnff", "-gfnff.topology_mode", "react",
               "-gfnff.cache_topology", "false", "-charge", str(info["charge"]),
               "-no_bmt", "-threads", "1", "-verbosity", "1"] + extra
        r = subprocess.run(cmd, cwd=d, capture_output=True, text=True, timeout=300)
        for ln in r.stdout.splitlines():
            if "Single Point Energy" in ln:
                return float(ln.split("=")[-1].split()[0])
        return None
    finally:
        shutil.rmtree(d, ignore_errors=True)


Eref = [p["energy_eh"] for p in info["points"]]
dref = [(e - Eref[-1]) * H2K for e in Eref]
ok, reasons = True, []
# the DEFAULT arm passes no flag at all - that is the point of the swap: it must be the good one
for label, extra, expect in (("default", [], "floor"),
                             ("delivered", ["-gfnff.rev_share_form", "delivered"], None)):
    E = [sp(f, extra) for f in frames]
    if any(e is None for e in E):
        print(f"FAIL: {label} arm has {sum(e is None for e in E)} failed single points")
        sys.exit(1)
    dcur = [(e - E[-1]) * H2K for e in E]
    dev = [c - r for c, r in zip(dcur, dref)]
    print(f"{label:11s} dev min {min(dev):+8.1f}  at 1.0/1.2/1.3 A: "
          f"{dev[0]:+7.1f} {dev[2]:+7.1f} {dev[3]:+7.1f} kcal/mol")
    if expect is None:
        if abs(min(dev) - delivered_min) > delivered_tol:
            ok = False
            reasons.append(f"delivered arm min dev {min(dev):+.1f} is more than {delivered_tol} "
                           f"from its recorded {delivered_min:+.1f} kcal/mol")
    else:
        if min(dev) < cons_floor:
            ok = False
            reasons.append(f"conserving arm min dev {min(dev):+.1f} below the floor "
                           f"{cons_floor:+.1f} kcal/mol - the valence-conserving share regressed")
if not ok:
    print("FAIL reasons: " + "; ".join(reasons))
sys.exit(0 if ok else 1)
PYEOF
)
    py_rc=$?
    set -e
    echo "$py_out"
    if [ $py_rc -eq 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: the default is above the floor, the delivered arm still shows the old profile"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: class-C adduct falsifier"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
