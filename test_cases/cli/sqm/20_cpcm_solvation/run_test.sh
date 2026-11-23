#!/bin/bash
# Test: native CPCM solvation (-xtb.solvent_model cpcm) for gfn1/gfn2 against tblite 0.4.0 --cpcm
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026). AI-generated, machine-tested only.
#
#   (a) single-point energies equal tblite 0.4.0 (commit 9ca469a) `tblite run --method gfnN --cpcm <solvent> --acc 0.0001
#       [--charge q] mol.xyz` (measured Nov 22, 2026 at the library geometries) within 5e-9 Eh; measured differences
#       <= 7.8e-10 Eh (caffeine gfn2), which equal the gas-phase differences of the native GFN methods to tblite
#   (b) a dielectric constant instead of a solvent name (-xtb.solvent 37.5, -xtb.solvent_epsilon 37.5) = tblite --cpcm 37.5;
#       a solvent without ALPB parameters (pyridine, tblite solvent table)
#   (c) -dump_gradient equals central finite differences of -sp energies (h = 2e-4 Angstrom) within 1e-6 Eh/Angstrom,
#       neutral (methanol, gfn2, water) and charged (NH4+, gfn1, hexane)
#   (d) errors: unknown solvent, a reference state other than gsolv (tblite --cpcm rejects it as well)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 20: native CPCM solvation (tblite model)"
TEST_DIR="$SCRIPT_DIR"

strip() { sed 's/\x1b\[[0-9;]*m//g' "$1"; }
exit_of() { grep -oP 'EXIT \K[0-9]+' "$1" | tail -1; }
# run <log> <xyz> <args...>: full-precision energy in <log>.grad
run() { local LOG=$1 XYZ=$2; shift 2; $CURCUMA -sp "$XYZ" -threads 1 -verbosity 1 -no_bmt -dump_gradient "${LOG%.log}.grad" "$@" > "$LOG" 2>&1; echo "EXIT $?" >> "$LOG"; }
eprec() { grep -oP '# energy \K[-0-9.e]+' "${1%.log}.grad" | head -1; }

expect_error() {      # log, text that must appear, description
    TESTS_RUN=$((TESTS_RUN + 1))
    local rc=$(exit_of "$1")
    if [ -n "$rc" ] && [ "$rc" != "0" ] && strip "$1" | grep -q -i "$2"; then
        echo -e "${GREEN}PASS${NC}: $3 (exit code $rc)"; TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: $3 (exit code $rc, expected non-zero and '$2')"; TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

# fd_check <tag> <xyz> <args...>: max |central FD - analytic| over all coordinates
fd_check() {
    local TAG=$1 XYZ=$2; shift 2
    run "${TAG}_ref.log" "$XYZ" -scf_threshold 1e-10 "$@"
    local NAT=$(head -1 "$XYZ")
    local i c
    for ((i = 0; i < NAT; i++)); do
        for c in 0 1 2; do
            for s in p m; do
                python3 - "$XYZ" "$i" "$c" "$s" > "${TAG}_${i}_${c}_${s}.xyz" <<'EOF'
import sys
path, i, c, s = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
lines = open(path).read().split("\n")
n = int(lines[0])
print(n)
print()
for k in range(n):
    p = lines[2 + k].split()
    x = [float(v) for v in p[1:4]]
    if k == i:
        x[c] += 2e-4 if s == "p" else -2e-4
    print(p[0], " ".join(f"{v:.12f}" for v in x))
EOF
                run "${TAG}_${i}_${c}_${s}.log" "${TAG}_${i}_${c}_${s}.xyz" -scf_threshold 1e-10 "$@"
            done
        done
    done
    python3 - "$TAG" "$NAT" <<'EOF'
import sys
tag, nat = sys.argv[1], int(sys.argv[2])
def energy(f):
    return float([l for l in open(f) if l.startswith("# energy")][0].split()[2])
g = [[float(v) for v in l.split()] for l in open(f"{tag}_ref.grad") if l.strip() and not l.startswith("#")]
d = max(abs((energy(f"{tag}_{i}_{c}_p.grad") - energy(f"{tag}_{i}_{c}_m.grad")) / 4e-4 - g[i][c])
        for i in range(nat) for c in range(3))
print(f"{d:.3e}")
EOF
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts
    rm -f *_ref.* fd*_*.xyz fd*_*.log fd*_*.grad

    # (a) energies against tblite 0.4.0 --cpcm (--acc 0.0001)
    check() {   # tag, xyz, tblite energy, args...
        local TAG=$1 XYZ=$2 REF=$3; shift 3
        run "$TAG.log" "$XYZ" -scf_threshold 1e-10 "$@"
        assert_numeric_match "$REF" "$(eprec $TAG.log)" "0.000000005" "$TAG = tblite 0.4.0 --cpcm"
    }
    check w_g2_water    water.xyz     -5.0683465147069  -method gfn2 -xtb.solvent water  -xtb.solvent_model cpcm
    check w_g1_hexane   water.xyz     -5.7672579257421  -method gfn1 -xtb.solvent hexane -xtb.solvent_model cpcm
    check m_g1_water    methanol.xyz  -8.9593915923690  -method gfn1 -xtb.solvent water  -xtb.solvent_model cpcm
    check m_g2_hexane   methanol.xyz  -8.2256826092880  -method gfn2 -xtb.solvent hexane -xtb.solvent_model cpcm
    check c_g2_water    caffeine.xyz -42.155082143933   -method gfn2 -xtb.solvent water  -xtb.solvent_model cpcm
    check c_g1_hexane   caffeine.xyz -44.512413147135   -method gfn1 -xtb.solvent hexane -xtb.solvent_model cpcm
    check n_g1_water    ammonium.xyz  -5.0077602836452  -method gfn1 -xtb.solvent water  -xtb.solvent_model cpcm -charge 1
    check n_g2_hexane   ammonium.xyz  -4.5429258193887  -method gfn2 -xtb.solvent hexane -xtb.solvent_model cpcm -charge 1
    check oh_g2_water   hydroxyl.xyz  -4.7780281172896  -method gfn2 -xtb.solvent water  -xtb.solvent_model cpcm -charge -1
    check oh_g1_hexane  hydroxyl.xyz  -5.4055810787050  -method gfn1 -xtb.solvent hexane -xtb.solvent_model cpcm -charge -1
    check s_g2_water    h2s.xyz       -4.2532945503669  -method gfn2 -xtb.solvent water  -xtb.solvent_model cpcm
    check s_g1_hexane   h2s.xyz       -4.7691127629204  -method gfn1 -xtb.solvent hexane -xtb.solvent_model cpcm

    # (b) dielectric constant instead of a name (tblite --cpcm 37.5), and a solvent of tblite's table without ALPB set
    check m_g2_eps_name methanol.xyz  -8.2252391162     -method gfn2 -xtb.solvent 37.5 -xtb.solvent_model cpcm
    check m_g2_eps_flag methanol.xyz  -8.2252391162     -method gfn2 -xtb.solvent_model cpcm -xtb.solvent_epsilon 37.5
    check m_g1_pyridine methanol.xyz  -8.9594935763     -method gfn1 -xtb.solvent pyridine -xtb.solvent_model cpcm

    # (c) analytic gradient against finite differences of the energy
    local D
    D=$(fd_check fdm methanol.xyz -method gfn2 -xtb.solvent water -xtb.solvent_model cpcm)
    assert_numeric_match "0" "$D" "0.000001" "gfn2 methanol water/CPCM: max |FD - analytic| (Eh/A)"
    D=$(fd_check fdn ammonium.xyz -method gfn1 -xtb.solvent hexane -xtb.solvent_model cpcm -charge 1)
    assert_numeric_match "0" "$D" "0.000001" "gfn1 NH4+ hexane/CPCM: max |FD - analytic| (Eh/A)"

    # (d) errors
    run e1.log water.xyz -method gfn2 -xtb.solvent watr -xtb.solvent_model cpcm
    expect_error e1.log "unknown solvent" "cpcm: unknown solvent is an error"
    run e2.log water.xyz -method gfn2 -xtb.solvent water -xtb.solvent_model cpcm -solvent_state bar1M
    expect_error e2.log "solvent_state" "cpcm: solvent_state bar1M is an error"

    rm -f fd*_*_*.xyz fd*_*_*.log fd*_*_*.grad
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
