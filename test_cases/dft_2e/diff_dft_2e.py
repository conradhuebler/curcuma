#!/usr/bin/env python3
"""
WP2 4-centre-ERI validation orchestrator (Claude Generated, WP2, June 2026).

Pure-stdlib (no numpy) so it runs under any CMake-found Python3, matching the
project's validate_sqm.py / diff_dft_1e.py convention. Drives three independent
checks of the curcuma native-DFT 2-electron integrals:

  (a) KERNEL GATE (always run): element-wise ERI (cartesian 6d, chemists'
      (mu nu | lam sig)) from the curcuma dumper (dump_dft_2e) vs an independent
      Python McMurchie-Davidson witness (scripts/dft_2e_python_ints.py) that
      shares curcuma's exact cartesian AO order. Tolerance --tol (default 1e-10).
      This validates the ERI KERNEL directly, independent of any density.

  (b) INTERNAL CONSISTENCY (always run, on the cartesian dump):
        - nbf matches the witness
        - 8-fold ERI symmetry: max|ERI(a,b,c,d) - ERI(perm)| <= 1e-12 over all
          quartets and all 8 index permutations
        - J and K symmetric to 1e-12
        - Tr(P*S) = 2 (dummy closed-shell density, one spatial orbital)
        - Tr(P*J) == Tr(P*K) to 1e-9 (single-orbital identity: for P = 2 c c^T
          the standard J/K contractions satisfy Tr(PJ)==Tr(PK) for ANY 4-tensor,
          by dummy-index relabeling b<->c in the exchange contraction)
        - J and K element-wise vs the witness (same dummy P) to --tol, checking
          the C++ buildCoulomb/buildExchange index conventions

  (c) xcDFT CROSS-CHECK (optional): if --xcDFT-ref points at an He-VDZ ERI file
      (i j k l value, 1-indexed, unique elements, as produced by xcDFT's
      read_integrals), compare curcuma's He-VDZ ERI at --tol. Gated behind the
      flag so a checkout without the reference still runs (a)+(b).

Exit code 0 = all run checks pass; nonzero otherwise.

  diff_dft_2e.py --dump <dump_dft_2e> --python <witness.py> --xyz <file>
                 [--basis NAME] [--tol 1e-10] [--xcDFT-ref FILE] [--quiet]

Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>. GPL-3.0.
"""
import argparse, json, math, os, subprocess, sys


class Fail(Exception):
    pass


def run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise Fail("command failed (%d): %s\n%s" % (p.returncode, " ".join(cmd), p.stderr))
    return p.stdout


def mat(J):
    return [[float(x) for x in row] for row in J]


def check_symmetry(name, M, tol):
    n = len(M)
    d = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            a = abs(M[i][j] - M[j][i])
            if a > d:
                d = a
    if d > tol:
        raise Fail("%s not symmetric: max|M-MT|=%.3e > %.3e" % (name, d, tol))
    return d


def trace(A, B):
    """Tr(A B) = sum_ij A[i][j] B[j][i] (B symmetric here so == sum A_ij B_ij)."""
    n = len(A)
    s = 0.0
    for i in range(n):
        for j in range(n):
            s += A[i][j] * B[j][i]
    return s


# The 8 chemists'-notation permutations of (mu,nu,lam,sig) that leave the ERI
# invariant: swap within the bra {mu,nu}, within the ket {lam,sig}, and swap the
# two pairs. The canonical element itself is the first entry.
PERMS = [
    lambda a, b, c, d: (a, b, c, d),
    lambda a, b, c, d: (b, a, c, d),
    lambda a, b, c, d: (a, b, d, c),
    lambda a, b, c, d: (b, a, d, c),
    lambda a, b, c, d: (c, d, a, b),
    lambda a, b, c, d: (d, c, a, b),
    lambda a, b, c, d: (c, d, b, a),
    lambda a, b, c, d: (d, c, b, a),
]


def eight_fold_max(E, n):
    """max over all quartets and all 8 perms of |E(a,b,c,d) - E(perm)|."""
    def at(a, b, c, d):
        return ((a * n + b) * n + c) * n + d
    worst = 0.0
    for mu in range(n):
        for nu in range(n):
            for lam in range(n):
                for sig in range(n):
                    v = E[at(mu, nu, lam, sig)]
                    for p in PERMS:
                        a, b, c, d = p(mu, nu, lam, sig)
                        diff = abs(v - E[at(a, b, c, d)])
                        if diff > worst:
                            worst = diff
    return worst


def max_abs_flat(A, B):
    d = 0.0
    for i in range(len(A)):
        a = abs(A[i] - B[i])
        if a > d:
            d = a
    return d


def max_abs_mat(A, B):
    d = 0.0
    for i in range(len(A)):
        for j in range(len(A[0])):
            a = abs(A[i][j] - B[i][j])
            if a > d:
                d = a
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--python", required=True)
    ap.add_argument("--xyz", required=True)
    ap.add_argument("--basis", default="def2-SVP")
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--charge", type=float, default=0.0)
    ap.add_argument("--spin", type=int, default=0)
    ap.add_argument("--xcDFT-ref", default=None,
                   help="optional xcDFT He-VDZ ERI reference file (i j k l value, 1-indexed)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    verbose = not args.quiet
    failures = []

    def log(msg):
        if verbose:
            print(msg)

    stem = os.path.splitext(os.path.basename(args.xyz))[0]

    # ---- (a) kernel gate + (b) internal consistency ----
    try:
        curc = json.loads(run([args.dump, args.xyz, "--basis", args.basis,
                               "--charge", str(args.charge), "--spin", str(args.spin)]))
        py = json.loads(run([sys.executable, args.python, args.xyz,
                            "--basis", args.basis, "--charge", str(args.charge)]))

        n = curc["nbf"]
        if n != py["nbf"]:
            raise Fail("nbf mismatch: curcuma=%d witness=%d" % (n, py["nbf"]))
        if len(curc["ERI"]) != n ** 4 or len(py["ERI"]) != n ** 4:
            raise Fail("ERI length mismatch: curcuma=%d witness=%d (expected %d)"
                       % (len(curc["ERI"]), len(py["ERI"]), n ** 4))

        # (a) element-wise ERI kernel gate
        d_eri = max_abs_flat(curc["ERI"], py["ERI"])
        if d_eri > args.tol:
            raise Fail("kernel ERI max|curc-witness|=%.3e > %.3e" % (d_eri, args.tol))
        log("  kernel ERI max|curc-witness| = %.3e" % d_eri)

        E = curc["ERI"]

        # (b1) 8-fold symmetry
        sym_eri = eight_fold_max(E, n)
        if sym_eri > 1e-12:
            raise Fail("8-fold ERI symmetry max|diff|=%.3e > 1e-12" % sym_eri)
        log("  8-fold ERI symmetry max|diff| = %.3e" % sym_eri)

        # (b2) J / K symmetric
        S_c = mat(curc["S"]); P_c = mat(curc["P"]); J_c = mat(curc["J"]); K_c = mat(curc["K"])
        check_symmetry("J", J_c, 1e-12)
        check_symmetry("K", K_c, 1e-12)
        check_symmetry("S", S_c, 1e-12)
        log("  J / K / S symmetric to 1e-12 OK")

        # (b3) Tr(P S) = 2 (dummy closed-shell, one spatial orbital, 2 electrons)
        trPS = trace(P_c, S_c)
        if abs(trPS - 2.0) > 1e-9:
            raise Fail("Tr(P*S) != 2: %.12f" % trPS)
        log("  Tr(P*S) = %.12f" % trPS)

        # (b4) Tr(P J) == Tr(P K) (single-orbital identity)
        trPJ = trace(P_c, J_c)
        trPK = trace(P_c, K_c)
        dJK = abs(trPJ - trPK)
        if dJK > 1e-9:
            raise Fail("Tr(PJ) != Tr(PK): |%.12f - %.12f|=%.3e > 1e-9" % (trPJ, trPK, dJK))
        log("  Tr(PJ)=%.12f  Tr(PK)=%.12f  |diff|=%.3e" % (trPJ, trPK, dJK))

        # (b5) J / K element-wise vs witness (same dummy P -> validates the C++
        #      buildCoulomb/buildExchange index conventions)
        J_p = mat(py["J"]); K_p = mat(py["K"]); P_p = mat(py["P"])
        d_P = max_abs_mat(P_c, P_p)
        d_J = max_abs_mat(J_c, J_p)
        d_K = max_abs_mat(K_c, K_p)
        if d_P > 1e-10:
            raise Fail("dummy P max|curc-witness|=%.3e > 1e-10 (k mismatch?)" % d_P)
        if d_J > args.tol:
            raise Fail("J max|curc-witness|=%.3e > %.3e" % (d_J, args.tol))
        if d_K > args.tol:
            raise Fail("K max|curc-witness|=%.3e > %.3e" % (d_K, args.tol))
        log("  dummy P max|curc-witness| = %.3e" % d_P)
        log("  J max|curc-witness| = %.3e" % d_J)
        log("  K max|curc-witness| = %.3e" % d_K)

        log("[%s] KERNEL + INTERNAL OK (nbf=%d, tol=%.0e)" % (stem, n, args.tol))

    except Fail as e:
        failures.append("(a/b) %s" % e)

    # ---- (c) optional xcDFT He-VDZ cross-check ----
    if args.xcDFT_ref and os.path.exists(args.xcDFT_ref):
        try:
            # Parse the xcDFT reference: lines "i j k l value" (1-indexed,
            # unique elements). Build a flat n^4 tensor exploiting 8-fold
            # symmetry, then compare to curcuma's cartesian ERI.
            refs = {}
            n_ref = 0
            with open(args.xcDFT_ref) as f:
                for line in f:
                    t = line.split()
                    if len(t) < 5:
                        continue
                    i, j, k, l = int(t[0]), int(t[1]), int(t[2]), int(t[3])
                    v = float(t[4])
                    refs[(i, j, k, l)] = v
                    n_ref = max(n_ref, i, j, k, l)
            n_ref = n_ref  # 1-indexed size
            Eref = [0.0] * (n_ref ** 4)
            def atref(a, b, c, d):
                return ((a * n_ref + b) * n_ref + c) * n_ref + d
            for (i, j, k, l), v in refs.items():
                a, b, c, d = i - 1, j - 1, k - 1, l - 1
                for p in PERMS:
                    pa, pb, pc, pd = p(a, b, c, d)
                    Eref[atref(pa, pb, pc, pd)] = v
            if curc["nbf"] != n_ref:
                raise Fail("xcDFT ref nbf=%d != curcuma nbf=%d" % (n_ref, curc["nbf"]))
            d_xc = max_abs_flat(curc["ERI"], Eref)
            if d_xc > args.tol:
                raise Fail("xcDFT ERI max|curc-xcDFT|=%.3e > %.3e" % (d_xc, args.tol))
            log("  xcDFT He-VDZ ERI max|curc-xcDFT| = %.3e" % d_xc)
            log("[%s] xcDFT CROSS-CHECK OK (tol=%.0e)" % (stem, args.tol))
        except Fail as e:
            failures.append("(c) %s" % e)
    elif verbose and args.xcDFT_ref:
        log("[%s] xcDFT ref not found -> cross-check skipped" % stem)

    if failures:
        for f in failures:
            sys.stderr.write("FAIL %s: %s\n" % (stem, f))
        sys.exit(1)
    if not verbose:
        print("PASS %s" % stem)
    sys.exit(0)


if __name__ == "__main__":
    main()