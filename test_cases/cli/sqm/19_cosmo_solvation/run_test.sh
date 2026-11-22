#!/bin/bash
# Test: native ddCOSMO solvation (-xtb.solvent_model cosmo) for gfn1/gfn2 against xtb 6.7.1 --cosmo
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026). AI-generated, machine-tested only.
#
#   (a) single-point energies equal xtb 6.7.1 `xtb mol.xyz --gfn N --cosmo <solvent> --sp --acc 0.01 [--chrg q]`
#       (measured Nov 22, 2026 at the library geometries) within 1e-7 Eh; measured differences <= 2.8e-8 Eh, which is
#       the gas-phase difference of the native GFN methods to xtb (caffeine gas 3.1e-8 Eh)
#   (b) a dielectric constant instead of a solvent name (-xtb.solvent 37.5, -xtb.solvent_epsilon 37.5) = xtb --cosmo 37.5
#   (c) -dump_gradient equals central finite differences of -sp energies (h = 2e-4 Angstrom) within 1e-6 Eh/Angstrom,
#       neutral (water, gfn2) and charged (NH4+, gfn1)
#   (d) errors: unknown solvent, a reference state other than gsolv (xtb adds no state shift to COSMO)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 19: native COSMO (ddCOSMO) solvation"
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

    # (a) energies against xtb 6.7.1
    check() {   # tag, xyz, xtb energy, args...
        local TAG=$1 XYZ=$2 REF=$3; shift 3
        run "$TAG.log" "$XYZ" "$@"
        assert_numeric_match "$REF" "$(eprec $TAG.log)" "0.0000001" "$TAG = xtb 6.7.1 --cosmo"
    }
    check w_g2_water   water.xyz     -5.0714227341 -method gfn2 -xtb.solvent water   -xtb.solvent_model cosmo
    check w_g1_acetone water.xyz     -5.7790550387 -method gfn1 -xtb.solvent acetone -xtb.solvent_model cosmo
    check m_g2_acetone methanol.xyz  -8.2419654777 -method gfn2 -xtb.solvent acetone -xtb.solvent_model cosmo
    check m_g1_water   methanol.xyz  -8.9677720465 -method gfn1 -xtb.solvent water   -xtb.solvent_model cosmo
    check c_g2_water   caffeine.xyz -42.1640854347 -method gfn2 -xtb.solvent water   -xtb.solvent_model cosmo
    check n_g2_water   ammonium.xyz  -4.6324413624 -method gfn2 -xtb.solvent water   -xtb.solvent_model cosmo -charge 1
    check n_g1_acetone ammonium.xyz  -5.0400927154 -method gfn1 -xtb.solvent acetone -xtb.solvent_model cosmo -charge 1

    # (b) dielectric constant instead of a name (xtb --cosmo 37.5)
    check m_g2_eps_name methanol.xyz -8.230802012570 -method gfn2 -xtb.solvent 37.5 -xtb.solvent_model cosmo
    check m_g2_eps_flag methanol.xyz -8.230802012570 -method gfn2 -xtb.solvent_model cosmo -xtb.solvent_epsilon 37.5

    # (c) analytic gradient against finite differences of the energy
    local D
    D=$(fd_check fdw water.xyz -method gfn2 -xtb.solvent water -xtb.solvent_model cosmo)
    assert_numeric_match "0" "$D" "0.000001" "gfn2 water/COSMO: max |FD - analytic| (Eh/A)"
    D=$(fd_check fdn ammonium.xyz -method gfn1 -xtb.solvent acetone -xtb.solvent_model cosmo -charge 1)
    assert_numeric_match "0" "$D" "0.000001" "gfn1 NH4+ acetone/COSMO: max |FD - analytic| (Eh/A)"

    # (d) errors
    run e1.log water.xyz -method gfn2 -xtb.solvent watr -xtb.solvent_model cosmo
    expect_error e1.log "unknown solvent" "cosmo: unknown solvent is an error"
    run e2.log water.xyz -method gfn2 -xtb.solvent water -xtb.solvent_model cosmo -solvent_state bar1M
    expect_error e2.log "solvent_state" "cosmo: solvent_state bar1M is an error"

    rm -f fd*_*_*.xyz fd*_*_*.log fd*_*_*.grad
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
