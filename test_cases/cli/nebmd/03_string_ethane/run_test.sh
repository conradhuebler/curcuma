#!/bin/bash
# NEB-MD CLI Test 03 - scheme 'string' (finite-temperature string, stage 1)
# Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
# Claude Generated 2026 (AI-generated, machine-tested)
#
# Ethane (staggered -> eclipsed), GFN-FF, 6 images, 2 ps, mean-path update every 400 steps.
# Checks: exit 0; start.neb.string.xyz with 6 frames of 8 atoms, all finite; history with
# >= 4 updates and last delta_max < first; no reorder message; restraint never capped;
# consecutive frames are Kabsch-aligned (plain RMSD <= best-fit RMSD + 1e-6), which is what
# the stage-2 consumers (PathCV lambda, path spacing) assume. History has the trailing columns
# scale,n_dead with scale == 1 (no trust-radius scaling) and n_dead == 0 for every update.
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

pass() { TESTS_RUN=$((TESTS_RUN + 1)); TESTS_PASSED=$((TESTS_PASSED + 1)); echo -e "${GREEN}ok PASS${NC}: $1"; }
fail() { TESTS_RUN=$((TESTS_RUN + 1)); TESTS_FAILED=$((TESTS_FAILED + 1)); echo -e "${RED}xx FAIL${NC}: $1"; }

run_test() {
    cd "$TEST_DIR"
    $CURCUMA -nebmd start.xyz end.xyz \
        -nebmd.scheme string -nebmd.nimages 6 -nebmd.string_update 400 \
        -nebmd.string_max_time 2000 \
        -simplemd.method gfnff -simplemd.time_step 0.5 -simplemd.temperature 300 \
        -simplemd.thermostat csvr \
        > stdout.log 2> stderr.log
    assert_exit_code $? 0 "NEB-MD string should succeed"
    local sxyz=$(find_output_file "start.neb.string.xyz")
    assert_file_exists "$sxyz" "string path"
    local hist=$(find_output_file "start.neb.string_history.csv")
    assert_file_exists "$hist" "string history"
    return 0
}

validate_results() {
    local sxyz=$(find_output_file "start.neb.string.xyz")
    local hist=$(find_output_file "start.neb.string_history.csv")
    [ -z "$sxyz" ] || [ ! -f "$sxyz" ] && return 1
    [ -z "$hist" ] || [ ! -f "$hist" ] && return 1
    local clean=$(sed 's/\x1b\[[0-9;]*m//g' stdout.log)

    local nframes=$(grep -c "string update=" "$sxyz" || true)
    [ "${nframes:-0}" -eq 6 ] && pass "string.xyz has 6 frames" || fail "string.xyz has ${nframes:-0} frames, expected 6"

    local bad_nat=$(awk '/^[ \t]*[0-9]+[ \t]*$/ {if ($1 != 8) n++} END {print n+0}' "$sxyz")
    [ "$bad_nat" -eq 0 ] && pass "every frame has 8 atoms" || fail "$bad_nat frame headers with atom count != 8"

    local nonfinite=$(grep -ciE "nan|inf" "$sxyz" "$hist" | awk -F: '{s+=$NF} END{print s+0}')
    [ "$nonfinite" -eq 0 ] && pass "no NaN/Inf in string.xyz / history" || fail "$nonfinite NaN/Inf occurrences"

    local nupd=$(($(wc -l < "$hist") - 1))
    [ "$nupd" -ge 4 ] && pass "history has $nupd updates (>= 4)" || fail "history has $nupd updates, expected >= 4"

    # Claude Generated 2026: robustness columns appended at the END (columns 1-7 unchanged).
    # Ethane updates are <= 0.004 A, far below string_max_step = 0.05 -> scale 1 everywhere.
    local header=$(head -1 "$hist")
    [ "$header" = "update,step,delta_max,delta_rms,arclength,perp_max,below_tol,scale,n_dead,delta_perp_max,f_rms_max,f_se_max,f_sig_max,z_max,gain_max,prep_shift" ] \
        && pass "history header has the trailing columns scale,n_dead" || fail "history header '$header'"
    local nscale_bad=$(awk -F, 'NR>1 && ($8 == "" || $8 + 0 != 1) {n++} END {print n+0}' "$hist")
    [ "$nscale_bad" -eq 0 ] && pass "scale == 1 for all $nupd updates" || fail "$nscale_bad update(s) with scale != 1"
    local ndead_bad=$(awk -F, 'NR>1 && ($9 == "" || $9 + 0 != 0) {n++} END {print n+0}' "$hist")
    [ "$ndead_bad" -eq 0 ] && pass "n_dead == 0 for all updates" || fail "$ndead_bad update(s) with n_dead != 0"

    local first=$(awk -F, 'NR==2 {print $3}' "$hist")
    local last=$(awk -F, 'END {print $3}' "$hist")
    if awk -v a="$last" -v b="$first" 'BEGIN {exit !(a < b)}'; then
        pass "last delta_max $last < first $first"
    else
        fail "last delta_max $last not below first $first"
    fi

    local nreorder=$(echo "$clean" | grep -ci "reorder" || true)
    [ "${nreorder:-0}" -eq 0 ] && pass "no reorder message" || fail "$nreorder reorder message(s) in stdout"

    if echo "$clean" | grep -q "restraint force capped 0 times"; then
        pass "restraint force never capped"
    else
        fail "restraint cap counter not 0 (or line missing)"
    fi

    # Consecutive alignment: plain RMSD must equal the best-fit (unweighted Kabsch) RMSD.
    # Best fit via Horn's quaternion method, 4x4 Jacobi eigen solver - pure python, no numpy.
    local worst=$(python3 - "$sxyz" <<'PYEOF'
import math, sys
def frames(fn):
    L = [l.rstrip("\n") for l in open(fn)]
    i, out = 0, []
    while i < len(L):
        if not L[i].strip():
            i += 1; continue
        n = int(L[i].split()[0])
        out.append([[float(x) for x in L[i + 2 + k].split()[1:4]] for k in range(n)])
        i += n + 2
    return out
def jacobi_max(A):
    A = [r[:] for r in A]
    for _ in range(100):
        off = sum(A[p][q] ** 2 for p in range(4) for q in range(4) if p != q)
        if off < 1e-30:
            break
        for p in range(4):
            for q in range(p + 1, 4):
                if abs(A[p][q]) < 1e-300:
                    continue
                th = 0.5 * math.atan2(2 * A[p][q], A[q][q] - A[p][p])
                c, s = math.cos(th), math.sin(th)
                for k in range(4):
                    akp, akq = A[k][p], A[k][q]
                    A[k][p], A[k][q] = c * akp - s * akq, s * akp + c * akq
                for k in range(4):
                    apk, aqk = A[p][k], A[q][k]
                    A[p][k], A[q][k] = c * apk - s * aqk, s * apk + c * aqk
    return max(A[k][k] for k in range(4))
def rmsd(X, Y, fit):
    n = len(X)
    if fit:
        cx = [sum(x[c] for x in X) / n for c in range(3)]
        cy = [sum(y[c] for y in Y) / n for c in range(3)]
        X = [[x[c] - cx[c] for c in range(3)] for x in X]
        Y = [[y[c] - cy[c] for c in range(3)] for y in Y]
        S = [[sum(X[a][i] * Y[a][j] for a in range(n)) for j in range(3)] for i in range(3)]
        (xx, xy, xz), (yx, yy, yz), (zx, zy, zz) = S
        N = [[xx + yy + zz, yz - zy, zx - xz, xy - yx],
             [yz - zy, xx - yy - zz, xy + yx, zx + xz],
             [zx - xz, xy + yx, -xx + yy - zz, yz + zy],
             [xy - yx, zx + xz, yz + zy, -xx - yy + zz]]
        g = sum(v * v for x in X for v in x) + sum(v * v for y in Y for v in y)
        return math.sqrt(max(g - 2 * jacobi_max(N), 0.0) / n)
    return math.sqrt(sum((X[a][c] - Y[a][c]) ** 2 for a in range(n) for c in range(3)) / n)
F = frames(sys.argv[1])
print("%.3e" % max(rmsd(F[i], F[i - 1], False) - rmsd(F[i], F[i - 1], True) for i in range(1, len(F))))
PYEOF
)
    if [ -n "$worst" ] && awk -v w="$worst" 'BEGIN {exit !(w <= 1e-6)}'; then
        pass "consecutive plain RMSD <= best-fit RMSD + 1e-6 (max diff $worst A)"
    else
        fail "consecutive frames not aligned: max(plain - bestfit) = ${worst:-n/a}"
    fi
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; rm -f start.neb.*; }

main() {
    test_header "NEB-MD string scheme (gfnff, ethane)"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
