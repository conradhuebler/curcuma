#!/usr/bin/env python3
"""
WP1 1e-integral validation orchestrator (Claude Generated, WP1, June 2026).

Pure-stdlib (no numpy) so it runs under any CMake-found Python3, matching the
project's validate_sqm.py convention. Drives three independent checks of the
curcuma native-QM 1e integrals:

  (a) KERNEL GATE (always run): element-wise S/T/V from the curcuma dumper
      (--cartesian_d, 6 cartesian d) vs an independent Python Obara-Saika
      witness (qm_1e_python_ints.py) that shares curcuma's exact AO order.
      Tolerance --tol (default 1e-10). This validates the integral KERNEL
      directly, independent of the spherical transform.

  (b) INTERNAL CONSISTENCY (always run, on the cartesian dump):
        - S, T, V, H symmetric to 1e-12
        - H = T + V element-wise to 1e-12
        - S_ii = 1 (renormalized contracted AOs) to 1e-9
        - nbf matches the witness
        - Tr(P*S) = 1 for a dummy density built from an S^{-1/2} column
          (idempotency proxy) to 1e-9

  (c) ORCA REFERENCE (optional): if <stem>.orca_ref.json sits next to the
      xyz (produced by scripts/qm_1e_reference.py via ORCA), compare curcuma's
      smallest SPHERICAL overlap eigenvalue to ORCA's reported "Smallest
      eigenvalue" (from the overlap diagonalisation). ORCA runs def2-SVP in the
      spherical 5d convention, so the spherical dump (without --cartesian_d) is
      compared, NOT the cartesian 6d set used in (a/b) (the 6d->5d projection is
      non-square, so the spectra differ). This is a WP1-valid 1e cross-check:
      AO-order-invariant and extremely basis-sensitive, it confirms that
      curcuma's def2-SVP.dat == ORCA's internal def2-SVP and that curcuma's
      cartesian->spherical d transform is correct. The full MO spectrum vs
      ORCA is a WP3+ deliverable (ORCA's orbital energies are the SCF Fock
      spectrum Hc+2J-K, not the 1e Hcore spectrum curcuma has at WP1), so
      diff_qm_1e.py does NOT assert it.

Exit code 0 = all run checks pass; nonzero otherwise.

  diff_qm_1e.py --dump <dump_qm_1e> --python <witness.py> --xyz <file>
                 [--basis NAME] [--tol 1e-10] [--tol-orca 1e-4] [--quiet]

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
    """JSON nested list -> list of lists of float."""
    return [[float(x) for x in row] for row in J]


def max_abs(A):
    m = 0.0
    for row in A:
        for x in row:
            a = abs(x)
            if a > m:
                m = a
    return m


def transpose(A):
    return [list(col) for col in zip(*A)]


def matmul(A, B):
    n, k, p = len(A), len(A[0]), len(B[0])
    BT = transpose(B)
    return [[sum(A[i][l] * BT[j][l] for l in range(k)) for j in range(p)] for i in range(n)]


def sym_matmul(A, B):
    return matmul(A, matmul(B, transpose(A)))


def jacobi_eig(A, max_sweeps=200, eps=1e-30):
    """Symmetric eigendecomposition A = V diag(w) V^T via cyclic Jacobi.

    Converges to machine precision: the off-diagonal Frobenius^2 is driven
    below eps (1e-30, i.e. ~1e-15 per element), with a generous sweep cap so the
    S^{-1/2} used for the Tr(P*S) idempotency proxy is accurate to ~1e-12.
    """
    n = len(A)
    V = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    w = [A[i][i] for i in range(n)]
    prev_off = float("inf")
    for _ in range(max_sweeps):
        off = 0.0
        for p in range(n):
            for q in range(p + 1, n):
                off += A[p][q] * A[p][q]
        if off < eps or off >= prev_off:
            break
        prev_off = off
        for p in range(n):
            for q in range(p + 1, n):
                apq = A[p][q]
                if abs(apq) < 1e-300:
                    continue
                app = A[p][p]
                aqq = A[q][q]
                tau = (aqq - app) / (2.0 * apq)
                if tau >= 0:
                    t = 1.0 / (tau + math.sqrt(1.0 + tau * tau))
                else:
                    t = -1.0 / (-tau + math.sqrt(1.0 + tau * tau))
                c = 1.0 / math.sqrt(1.0 + t * t)
                s = t * c
                for i in range(n):
                    aip = A[i][p]
                    aiq = A[i][q]
                    A[i][p] = c * aip - s * aiq
                    A[i][q] = s * aip + c * aiq
                for i in range(n):
                    api = A[p][i]
                    aqi = A[q][i]
                    A[p][i] = c * api - s * aqi
                    A[q][i] = s * api + c * aqi
                for i in range(n):
                    vip = V[i][p]
                    viq = V[i][q]
                    V[i][p] = c * vip - s * viq
                    V[i][q] = s * vip + c * viq
    w = [A[i][i] for i in range(n)]
    return w, V


def sinv_half(S):
    w, V = jacobi_eig([row[:] for row in S])
    wmin = min(w)
    if wmin < 1e-12:
        raise Fail("S has near-zero eigenvalue: %.3e" % wmin)
    D = [[1.0 / math.sqrt(w[i]) if i == j else 0.0 for j in range(len(w))] for i in range(len(w))]
    # S^{-1/2} = V D V^T
    return sym_matmul(V, D)


def mo_energies(H, S):
    """Generalized eigenvalue H C = S C eps via symmetric orthogonalization."""
    X = sinv_half(S)
    Hp = sym_matmul(X, H)  # X H X^T (X symmetric)
    eps, _ = jacobi_eig([row[:] for row in Hp])
    return sorted(eps)


def dummy_density_tracePS(S):
    """Tr(P*S) for P = s s^T, s = S^{-1/2}[:,k] (column of largest norm)."""
    X = sinv_half(S)
    k = max(range(len(X)), key=lambda c: sum(X[i][c] ** 2 for i in range(len(X))))
    n = len(S)
    tr = 0.0
    for i in range(n):
        for j in range(n):
            tr += X[i][k] * X[j][k] * S[i][j]
    return tr


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--python", required=True)
    ap.add_argument("--xyz", required=True)
    ap.add_argument("--basis", default="def2-SVP")
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--tol-orca", type=float, default=1e-4,
                   help="ORCA smallest-S-eigenvalue tolerance (ORCA prints ~4 sig figs)")
    ap.add_argument("--charge", type=float, default=0.0)
    ap.add_argument("--spin", type=int, default=0)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    verbose = not args.quiet
    failures = []

    def log(msg):
        if verbose:
            print(msg)

    stem = os.path.splitext(os.path.basename(args.xyz))[0]
    xyz_dir = os.path.dirname(os.path.abspath(args.xyz))

    # ---- (a) kernel gate: cartesian curcuma vs python witness ----
    try:
        curc_cart = json.loads(run([args.dump, args.xyz, "--basis", args.basis,
                                    "--cartesian_d", "--charge", str(args.charge),
                                    "--spin", str(args.spin)]))
        py_json = json.loads(run([sys.executable, args.python, args.xyz,
                                 "--basis", args.basis, "--charge", str(args.charge)]))

        S_c = mat(curc_cart["S"]); T_c = mat(curc_cart["T"])
        V_c = mat(curc_cart["V"]); H_c = mat(curc_cart["H"])
        S_p = mat(py_json["S"]); T_p = mat(py_json["T"])
        V_p = mat(py_json["V"]); H_p = mat(py_json["H"])

        if curc_cart["nbf"] != py_json["nbf"]:
            raise Fail("nbf mismatch: curcuma=%d witness=%d" % (curc_cart["nbf"], py_json["nbf"]))
        if len(S_c) != len(S_p) or (S_c and len(S_c[0]) != len(S_p[0])):
            raise Fail("shape mismatch: curcuma %sx%s vs witness %sx%s"
                       % (len(S_c), len(S_c[0]) if S_c else 0, len(S_p), len(S_p[0]) if S_p else 0))

        for name, Mc, Mp in (("S", S_c, S_p), ("T", T_c, T_p),
                             ("V", V_c, V_p), ("H", H_c, H_p)):
            d = 0.0
            for i in range(len(Mc)):
                for j in range(len(Mc[0])):
                    a = abs(Mc[i][j] - Mp[i][j])
                    if a > d:
                        d = a
            if d > args.tol:
                raise Fail("kernel %s max|curc-witness|=%.3e > %.3e" % (name, d, args.tol))
            log("  kernel %-2s max|curc-witness| = %.3e" % (name, d))

        # ---- (b) internal consistency (on the cartesian dump) ----
        sym_tol = 1e-12
        for name, M in (("S", S_c), ("T", T_c), ("V", V_c), ("H", H_c)):
            check_symmetry(name, M, sym_tol)
        log("  symmetry S/T/V/H max|M-MT| <= %.0e OK" % sym_tol)

        dHTV = 0.0
        n = len(H_c)
        for i in range(n):
            for j in range(n):
                a = abs(H_c[i][j] - (T_c[i][j] + V_c[i][j]))
                if a > dHTV:
                    dHTV = a
        if dHTV > 1e-12:
            raise Fail("H != T+V: max|H-(T+V)|=%.3e" % dHTV)
        log("  H = T + V max|diff| = %.3e" % dHTV)

        dSii = max(abs(S_c[i][i] - 1.0) for i in range(n))
        if dSii > 1e-9:
            raise Fail("S_ii != 1: max|S_ii-1|=%.3e" % dSii)
        log("  S_ii = 1 max|S_ii-1| = %.3e" % dSii)

        trPS = dummy_density_tracePS(S_c)
        if abs(trPS - 1.0) > 1e-9:
            raise Fail("Tr(PS) != 1: %.12f" % trPS)
        log("  Tr(P*S) = %.12f (dummy idempotent density)" % trPS)

        log("[%s] KERNEL + INTERNAL OK (nbf=%d, tol=%.0e)" % (stem, curc_cart["nbf"], args.tol))

    except Fail as e:
        failures.append("(a/b) %s" % e)

    # ---- (c) optional ORCA smallest-overlap-eigenvalue cross-check ----
    # ORCA runs def2-SVP in the SPHERICAL 5d convention; the cartesian 6d set
    # carries an extra s-like (dxx+dyy+dzz) combination that is near-linearly-
    # dependent with the s block, so the cartesian smallest S eigenvalue is
    # systematically smaller than the spherical one. The 6d->5d transform is a
    # 5x6 projection, NOT a square congruence, so the spectra genuinely differ.
    # The cross-check MUST therefore use curcuma's spherical S (dump without
    # --cartesian_d), which is what ORCA diagonalises. This jointly validates
    # (i) curcuma's def2-SVP.dat == ORCA's internal def2-SVP (basis-sensitive,
    # AO-invariant) and (ii) curcuma's cartesian->spherical d transform.
    orca_ref = os.path.join(xyz_dir, stem + ".orca_ref.json")
    if os.path.exists(orca_ref):
        try:
            ref = json.load(open(orca_ref))
            s_min_orca = float(ref["s_min_eigenvalue"])
            curc_sph = json.loads(run([args.dump, args.xyz, "--basis", args.basis,
                                       "--charge", str(args.charge),
                                       "--spin", str(args.spin)]))
            S_sph = mat(curc_sph["S"])
            w_s, _ = jacobi_eig([row[:] for row in S_sph])
            s_min_curc = min(w_s)
            d = abs(s_min_curc - s_min_orca)
            if d > args.tol_orca:
                raise Fail("smallest S eigenvalue (spherical): curcuma=%.6e orca=%.6e diff=%.3e > %.3e"
                           % (s_min_curc, s_min_orca, d, args.tol_orca))
            log("  ORCA smallest S eigenvalue (spherical): curc=%.6e orca=%.6e diff=%.3e"
                % (s_min_curc, s_min_orca, d))
            log("[%s] ORCA S-EIGENVALUE OK (tol=%.0e)" % (stem, args.tol_orca))
        except Fail as e:
            failures.append("(c) %s" % e)
        except (KeyError, ValueError) as e:
            failures.append("(c) malformed orca_ref.json: %s" % e)
    elif verbose:
        log("[%s] no .orca_ref.json -> ORCA comparison skipped (run scripts/qm_1e_reference.py)" % stem)

    if failures:
        for f in failures:
            sys.stderr.write("FAIL %s: %s\n" % (stem, f))
        sys.exit(1)
    if not verbose:
        print("PASS %s" % stem)
    sys.exit(0)


if __name__ == "__main__":
    main()