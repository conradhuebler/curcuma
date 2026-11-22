#!/bin/bash

# Test: GFN-FF metal classification - one table (PeriodicTable::getMetalType) and the Sn/Pb/Bi small-CN rule
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026) - AI-generated, machine-tested; human production testing pending
#
# gfnff_par.h carried a second metal table that differed from the reference at Zn, Ga, Cd, Sn (1e-2..7e-2 Eh errors);
# it is removed. Sn/Pb/Bi with <= 4 neighbours count as non-metals in the EEQ ff terms and bond types
# (gfnff_ini.f90:234-235). Eight constructed systems (coordinates inline, not from the structure library),
# reference energies from xtb 6.7.1 --gfnff (pprcht/gfnff gave the same to 1e-8 Eh on Zn, Cd, Ga, Sn, Pb, Bi systems).

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../test_utils.sh"

TEST_NAME="gfnff - 11: metal type"

main() {
    test_header "$TEST_NAME"
    cd "$SCRIPT_DIR"
    set +e
    CURCUMA="$CURCUMA" python3 - <<'PYEOF'
import os, re, shutil, subprocess, sys, tempfile
B = os.environ["CURCUMA"]
CASES = [
    ("zncl4", -2, -2.28491114413, "Zn 0.000000 0.000000 0.000000;Cl 1.327906 1.327906 1.327906;Cl 1.327906 -1.327906 -1.327906;Cl -1.327906 1.327906 -1.327906;Cl -1.327906 -1.327906 1.327906"),
    ("znme2", 0, -1.103958004281, "Zn 0.000000 0.000000 0.000000;C 0.000000 0.000000 1.930000;H 1.027479 0.000000 2.293849;H -0.513740 0.889823 2.293849;H -0.513740 -0.889823 2.293849;C 0.000000 0.000000 -1.930000;H 0.513740 0.889823 -2.293849;H -1.027479 0.000000 -2.293849;H 0.513740 -0.889823 -2.293849"),
    ("cdme2", 0, -1.115222856133, "Cd 0.000000 0.000000 0.000000;C 0.000000 0.000000 2.110000;H 1.027479 0.000000 2.473849;H -0.513740 0.889823 2.473849;H -0.513740 -0.889823 2.473849;C 0.000000 0.000000 -2.110000;H 0.513740 0.889823 -2.473849;H -1.027479 0.000000 -2.473849;H 0.513740 -0.889823 -2.473849"),
    ("gacl3", 0, -0.169186765272, "Ga 0.000000 0.000000 0.000000;Cl 2.100000 0.000000 0.000000;Cl -1.050000 1.818653 0.000000;Cl -1.050000 -1.818653 0.000000"),
    ("sncl4", 0, -0.117496427496, "Sn 0.000000 0.000000 0.000000;Cl 1.316359 1.316359 1.316359;Cl 1.316359 -1.316359 -1.316359;Cl -1.316359 1.316359 -1.316359;Cl -1.316359 -1.316359 1.316359"),
    ("pbcl4", 0, -0.111489787533, "Pb 0.000000 0.000000 0.000000;Cl 1.402961 1.402961 1.402961;Cl 1.402961 -1.402961 -1.402961;Cl -1.402961 1.402961 -1.402961;Cl -1.402961 -1.402961 1.402961"),
    ("bicl3", 0, -0.108580290411, "Bi 0.000000 0.000000 0.000000;Cl 2.193689 0.000000 -1.156775;Cl -1.096844 1.899790 -1.156775;Cl -1.096844 -1.899790 -1.156775"),
    ("sncl6", -2, -2.373562507956, "Sn 0.000000 0.000000 0.000000;Cl 2.420000 0.000000 0.000000;Cl -2.420000 0.000000 0.000000;Cl 0.000000 2.420000 0.000000;Cl 0.000000 -2.420000 0.000000;Cl 0.000000 0.000000 2.420000;Cl 0.000000 0.000000 -2.420000"),
]
TOL = 6.0e-8
fails = 0
for name, q, ref, body in CASES:
    d = tempfile.mkdtemp()
    lines = body.split(";")
    open(os.path.join(d, name + ".xyz"), "w").write("%d\n%s\n%s\n" % (len(lines), name, "\n".join(lines)))
    p = subprocess.run([B, "-sp", name + ".xyz", "-method", "gfnff", "-charge", str(q), "-verbosity", "1"],
                       cwd=d, capture_output=True, text=True, stdin=subprocess.DEVNULL)
    m = re.search(r"Single Point Energy =\s*(-?\d+\.\d+)", p.stdout)
    shutil.rmtree(d, ignore_errors=True)
    e = float(m.group(1)) if m else None
    ok = e is not None and abs(e - ref) < TOL
    print("%s %-6s curcuma %s  xtb %.8f  diff %s" % ("PASS" if ok else "FAIL", name, e, ref, None if e is None else "%.2e" % (e - ref)))
    fails += 0 if ok else 1
sys.exit(1 if fails else 0)
PYEOF
    local rc=$?
    if [ $rc -eq 0 ]; then TESTS_PASSED=$((TESTS_PASSED + 1)); else TESTS_FAILED=$((TESTS_FAILED + 1)); fi
    TESTS_RUN=$((TESTS_RUN + 1))
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
