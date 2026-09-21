#!/usr/bin/env python3
"""
Independent Python witness for the curcuma native-DFT 4-centre ERI (WP2).

Reads an XYZ and the SAME def2-SVP.dat curcuma uses, builds the contracted
cartesian GTO basis in the SAME AO order as curcuma's
BasisSetParser::createGTOFromBasis (per atom in input order; per shell in file
order; p -> [px,py,pz]; d -> [dxx,dyy,dzz,dxy,dxz,dyz]), and computes the
4-centre electron-repulsion integral tensor (mu nu | lam sig) in CHEMISTS'
notation via an independently written McMurchie-Davidson implementation, plus
the Coulomb matrix J and exchange matrix K built from a dummy closed-shell
density P.

This mirrors scripts/dft_1e_python_ints.py (same basis parser, same
hermite_coeffs / boys_array primitives, same AO order) so the ERI can be
compared element by element against dump_dft_2e. It is pure-stdlib (no numpy /
scipy / pyscf) to run under any CMake-found Python3, matching the project's
validate_sqm.py / diff_dft_1e.py convention.

McMurchie-Davidson ERI (Helgaker, Molecular Electronic-Structure Theory, ch. 9.9;
McMurchie & Davidson, J. Comput. Phys. 26, 218 (1977)):

  (ab|cd) = (2 pi^(5/2))/(p q sqrt(p+q))
            * sum_{t,u,v} sum_{tau,ups,om}
              E^ab_{t,u,v} E^cd_{tau,ups,om} (-1)^(tau+ups+om)
              * R^0_{t+tau, u+ups, v+om}

  p = a+b, q = c+d, rho = pq/(p+q), Boys arg T = rho*|P-Q|^2.
  R^n_{0,0,0} = (-2 rho)^n F_n(T); R recurrence with the (P-Q) displacement
  (bra recurrence, plus sign) -- gamma-independent.

Output JSON schema (matches dump_dft_2e.cpp):
  { "molecule":{...}, "basis","cartesian_d":true, "nbf","num_electrons",
    "nuclear_repulsion", "S":[[...]], "eri_order":"mu_nu_lam_sig",
    "ERI":[...flat n^4...], "P":[[...]], "J":[[...]], "K":[[...]],
    "dummy_density":"2*c*c^T, c=ones normalized so c^T S c=1 (rank-1, no eig)" }

  dft_2e_python_ints.py <input.xyz> [--basis NAME] [--dat FILE] [--charge Q]

Claude Generated (WP2). GPL-3.0.
"""
import argparse, json, math, os, sys

PI = 3.14159265358979323846264338327950288
AA_TO_BOHR = 1.0 / 0.529177210903

SYM_TO_Z = {"H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10}
Z_TO_SYM = {v:k for k,v in SYM_TO_Z.items()}

SHELL_COMPONENTS = {
    "S": [(0,0,0)],
    "P": [(1,0,0),(0,1,0),(0,0,1)],
    "D": [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)],
}

def dfact(n):
    r = 1.0; k = 1
    while k <= n:
        r *= (2*k-1); k += 1
    return r if n > 0 else 1.0

def primitive_norm(alpha, l, m, n):
    L = l+m+n
    denom = math.sqrt(dfact(l)*dfact(m)*dfact(n)) or 1.0
    return (2.0*alpha/PI)**0.75 * (4.0*alpha)**(L/2.0) / denom

def overlap1d_table(imax, jmax, PA, PB, gamma):
    g2 = 0.5/gamma
    S = [[0.0]*(jmax+1) for _ in range(imax+1)]
    S[0][0] = 1.0
    for i in range(1, imax+1):
        prev = S[i-2][0] if i >= 2 else 0.0
        S[i][0] = PA*S[i-1][0] + g2*(i-1)*prev
    for i in range(imax+1):
        for j in range(1, jmax+1):
            si1 = S[i-1][j-1] if i >= 1 else 0.0
            sj1 = S[i][j-2]   if j >= 2 else 0.0
            S[i][j] = PB*S[i][j-1] + g2*(i*si1 + (j-1)*sj1)
    return S

def gaussian_product(a1, a2, A, B):
    gamma = a1+a2
    P = [(a1*A[k]+a2*B[k])/gamma for k in range(3)]
    zeta = a1*a2/gamma
    R2 = sum((A[k]-B[k])**2 for k in range(3))
    return gamma, P, math.exp(-zeta*R2)

def prim_overlap(la,ma,na, lb,mb,nb, a1,a2, A, B):
    gamma, P, K = gaussian_product(a1,a2,A,B)
    if K == 0.0: return 0.0
    pref = (PI/gamma)**1.5 * K
    Sx = overlap1d_table(max(la,lb), max(la,lb), P[0]-A[0], P[0]-B[0], gamma)
    Sy = overlap1d_table(max(ma,mb), max(ma,mb), P[1]-A[1], P[1]-B[1], gamma)
    Sz = overlap1d_table(max(na,nb), max(na,nb), P[2]-A[2], P[2]-B[2], gamma)
    return pref * Sx[la][lb] * Sy[ma][mb] * Sz[na][nb]

def hermite_coeffs(iA, iB, PA, PB, gamma, K):
    """E[t][i][j] Hermite expansion coefficients (K folded into E[0][0][0])."""
    tmax = iA+iB
    g2 = 0.5/gamma
    E = [[[0.0]*(iB+1) for _ in range(iA+1)] for _ in range(tmax+1)]
    if K == 0.0: return E
    E[0][0][0] = K
    for i in range(1, iA+1):
        prev = E[0][i-2][0] if i >= 2 else 0.0
        E[0][i][0] = PA*E[0][i-1][0] + g2*(i-1)*prev
    for j in range(1, iB+1):
        for i in range(iA+1):
            ei = E[0][i-1][j-1] if i >= 1 else 0.0
            ej = E[0][i][j-2]   if j >= 2 else 0.0
            E[0][i][j] = PB*E[0][i][j-1] + g2*(i*ei + (j-1)*ej)
    for t in range(1, tmax+1):
        for i in range(iA+1):
            for j in range(iB+1):
                ei = E[t-1][i-1][j] if i >= 1 else 0.0
                ej = E[t-1][i][j-1] if j >= 1 else 0.0
                et = E[t-2][i][j]   if t >= 2 else 0.0
                E[t][i][j] = PA*E[t-1][i][j] + g2*(i*ei + j*ej + t*et)
    return E

def boys_array(maxN, T):
    F = [0.0]*(maxN+1)
    if T < 1e-14:
        for n in range(maxN+1): F[n] = 1.0/(2*n+1)
        return F
    M = maxN + 25
    G = [0.0]*(M+1)
    eT = math.exp(-T)
    for n in range(M-1, -1, -1):
        G[n] = (2.0*T*G[n+1] + eT)/(2*n+1)
    for n in range(maxN+1): F[n] = G[n]
    return F

def build_R_eri(tmax, umax, vmax, Nmax, W, rho, boys):
    """ERI Coulomb auxiliary R[t][u][v][N] for the combined (bra+ket) Hermite
    orders. Base R^n_{0,0,0} = (-2*rho)^n * F_n(T); recurrence with the (P-Q)
    displacement W = P - Q (PLUS sign, bra recurrence)."""
    Wx, Wy, Wz = W
    R = [[[[0.0]*(Nmax+1) for _ in range(vmax+1)] for _ in range(umax+1)] for _ in range(tmax+1)]
    b = -2.0*rho; bp = 1.0
    for N in range(Nmax+1):
        R[0][0][0][N] = bp*boys[N]; bp *= b
    for t in range(tmax+1):
        if t >= 1:
            for N in range(Nmax):
                low = (t-1)*R[t-2][0][0][N+1] if t >= 2 else 0.0
                R[t][0][0][N] = low + Wx*R[t-1][0][0][N+1]
        for u in range(1, umax+1):
            for N in range(Nmax):
                low = (u-1)*R[t][u-2][0][N+1] if u >= 2 else 0.0
                R[t][u][0][N] = low + Wy*R[t][u-1][0][N+1]
        for u in range(umax+1):
            for v in range(1, vmax+1):
                for N in range(Nmax):
                    low = (v-1)*R[t][u][v-2][N+1] if v >= 2 else 0.0
                    R[t][u][v][N] = low + Wz*R[t][u][v-1][N+1]
    return R

def prim_eri(la,ma,na, lb,mb,nb, lc,mc,nc, ld,md,nd,
             a1,a2,a3,a4, A,B,C,D):
    """Primitive (ga gb | gc gd) in chemists' notation."""
    p, P, Kab = gaussian_product(a1, a2, A, B)
    if Kab == 0.0: return 0.0
    q, Q, Kcd = gaussian_product(a3, a4, C, D)
    if Kcd == 0.0: return 0.0
    Ex = hermite_coeffs(la, lb, P[0]-A[0], P[0]-B[0], p, Kab)
    Eu = hermite_coeffs(ma, mb, P[1]-A[1], P[1]-B[1], p, 1.0)
    Ev = hermite_coeffs(na, nb, P[2]-A[2], P[2]-B[2], p, 1.0)
    Ex2 = hermite_coeffs(lc, ld, Q[0]-C[0], Q[0]-D[0], q, Kcd)
    Eu2 = hermite_coeffs(mc, md, Q[1]-C[1], Q[1]-D[1], q, 1.0)
    Ev2 = hermite_coeffs(nc, nd, Q[2]-C[2], Q[2]-D[2], q, 1.0)
    tmax,umax,vmax = la+lb, ma+mb, na+nb
    taumax,upsmax,ommax = lc+ld, mc+md, nc+nd
    Ttot, Utot, Vtot = tmax+taumax, umax+upsmax, vmax+ommax
    Nmax = Ttot+Utot+Vtot
    rho = p*q/(p+q)
    T = rho*sum((P[k]-Q[k])**2 for k in range(3))
    F = boys_array(Nmax, T)
    R = build_R_eri(Ttot, Utot, Vtot, Nmax, [P[0]-Q[0],P[1]-Q[1],P[2]-Q[2]], rho, F)
    pref = 2.0*PI**2.5 / (p*q*math.sqrt(p+q))
    val = 0.0
    for t in range(tmax+1):
        et = Ex[t][la][lb]
        if et == 0.0: continue
        for u in range(umax+1):
            eu = Eu[u][ma][mb]
            if eu == 0.0: continue
            for v in range(vmax+1):
                ev = Ev[v][na][nb]
                if ev == 0.0: continue
                bra = et*eu*ev
                for tau in range(taumax+1):
                    et2 = Ex2[tau][lc][ld]
                    if et2 == 0.0: continue
                    for ups in range(upsmax+1):
                        eu2 = Eu2[ups][mc][md]
                        if eu2 == 0.0: continue
                        for om in range(ommax+1):
                            ev2 = Ev2[om][nc][nd]
                            if ev2 == 0.0: continue
                            sgn = -1.0 if (tau+ups+om) & 1 else 1.0
                            val += bra*et2*eu2*ev2*sgn*R[t+tau][u+ups][v+om][0]
    return pref*val

def contracted_eri(a, b, c, d):
    la,ma,na = a["lmn"]; lb,mb,nb = b["lmn"]; lc,mc,nc = c["lmn"]; ld,md,nd = d["lmn"]
    A,B,C,D = a["center"],b["center"],c["center"],d["center"]
    s = 0.0
    for ia in range(len(a["exps"])):
        ca = a["coeff"][ia]
        for ib in range(len(b["exps"])):
            cb = b["coeff"][ib]
            for ic in range(len(c["exps"])):
                cc = c["coeff"][ic]
                for idd in range(len(d["exps"])):
                    s += ca*cb*cc*d["coeff"][idd]*prim_eri(
                        la,ma,na, lb,mb,nb, lc,mc,nc, ld,md,nd,
                        a["exps"][ia], b["exps"][ib], c["exps"][ic], d["exps"][idd],
                        A,B,C,D)
    return s

def build_eri(orbs):
    """Full n^4 ERI tensor in chemists' (mu,nu|lam,sig), flat row-major
    index = ((mu*n+nu)*n+lam)*n+sig. Canonical-quartet loop with 8-fold fill."""
    n = len(orbs)
    E = [0.0]*(n*n*n*n)
    def at(mu,nu,lam,sig):
        # row-major 4-index, identical to C++ ERITensor::operator():
        # ((mu*n+nu)*n + lam)*n + sig = (mu*n+nu)*n*n + lam*n + sig.
        return ((mu*n+nu)*n + lam)*n + sig
    for mu in range(n):
        for nu in range(mu, n):
            p1 = mu*n+nu
            for lam in range(n):
                for sig in range(lam, n):
                    p2 = lam*n+sig
                    if p1 > p2: continue
                    v = contracted_eri(orbs[mu],orbs[nu],orbs[lam],orbs[sig])
                    # 8-fold chemists' symmetry
                    E[at(mu,nu,lam,sig)] = v
                    E[at(nu,mu,lam,sig)] = v
                    E[at(mu,nu,sig,lam)] = v
                    E[at(nu,mu,sig,lam)] = v
                    E[at(lam,sig,mu,nu)] = v
                    E[at(sig,lam,mu,nu)] = v
                    E[at(lam,sig,nu,mu)] = v
                    E[at(sig,lam,nu,mu)] = v
    return E, n

def mat_S(orbs):
    n = len(orbs)
    S = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            a, b = orbs[i], orbs[j]
            v = 0.0
            for ia in range(len(a["exps"])):
                for ib in range(len(b["exps"])):
                    v += a["coeff"][ia]*b["coeff"][ib]*prim_overlap(
                        *a["lmn"], *b["lmn"], a["exps"][ia], b["exps"][ib],
                        a["center"], b["center"])
            S[i][j] = S[j][i] = v
    return S

def dummy_density(S):
    """Closed-shell single-orbital dummy density P = 2 c c^T, where c is the
    ones-vector S-orthonormalized by the SCALAR c^T S c = 1 (no eigendecomposition).
    This is solver-independent and robust to degenerate S eigensubspaces (which
    made an argmax-column-of-S^-1/2 pick ambiguous on symmetric molecules such as
    linear BeH2). P is rank-1, so Tr(P J) == Tr(P K) holds for ANY ERI by
    dummy-index relabeling; Tr(P S) = 2 c^T S c = 2 (one spatial orbital, 2
    electrons). The fixed ones-vector exercises every AO index in J/K."""
    n = len(S)
    c = [1.0] * n
    q = 0.0
    for i in range(n):
        for j in range(n):
            q += c[i] * S[i][j] * c[j]
    inv = 1.0 / math.sqrt(q)
    c = [x * inv for x in c]
    return [[2.0 * c[i] * c[j] for j in range(n)] for i in range(n)]

def build_JK(E, n, P):
    """J_munu = sum P_lamsig (mu,nu|lam,sig); K_munu = sum P_lamsig (mu,lam|nu,sig)."""
    def eri(mu,nu,lam,sig):
        return ((mu*n+nu)*n + lam)*n + sig
    J = [[0.0]*n for _ in range(n)]
    K = [[0.0]*n for _ in range(n)]
    for mu in range(n):
        for nu in range(n):
            js = 0.0; ks = 0.0
            for lam in range(n):
                for sig in range(n):
                    Pls = P[lam][sig]
                    js += Pls*E[eri(mu,nu,lam,sig)]
                    ks += Pls*E[eri(mu,lam,nu,sig)]
            J[mu][nu] = js; K[mu][nu] = ks
    return J, K

def parse_basis_dat(path):
    with open(path) as f:
        lines = f.readlines()
    name_to_sym = {"HYDROGEN":"H","HELIUM":"He","LITHIUM":"Li","BERYLLIUM":"Be",
                   "BORON":"B","CARBON":"C","NITROGEN":"N","OXYGEN":"O",
                   "FLUORINE":"F","NEON":"Ne"}
    in_data = False
    i = 0
    basis = {}
    n = len(lines)
    while i < n and not lines[i].strip().startswith("$DATA"):
        i += 1
    i += 1
    cur = None
    while i < n:
        s = lines[i].strip()
        i += 1
        if not s: continue
        if s.startswith("$END"): break
        if s.startswith("!"): continue
        toks = s.split()
        first = toks[0]
        if first.isalpha() and first.upper() in name_to_sym and (len(toks)==1 or not toks[1][0].isdigit()):
            cur = name_to_sym[first.upper()]
            basis[cur] = []
            continue
        letter = first[0].upper()
        nprim = int(toks[1])
        ncontr = int(toks[2]) if len(toks) > 2 else 1
        exps = []
        coeffs = [[] for _ in range(ncontr)]
        for _ in range(nprim):
            p = lines[i].split()
            i += 1
            exps.append(float(p[1]))
            for c in range(ncontr):
                coeffs[c].append(float(p[2+c]))
        basis[cur].append({"letter":letter,"exps":exps,"coeffs":coeffs})
    return basis

def build_basis(basis_map, atoms, coord_bohr):
    orbs = []
    for iat, Z in enumerate(atoms):
        sym = Z_TO_SYM[Z]
        if sym not in basis_map:
            raise RuntimeError(f"no basis for {sym}")
        C = coord_bohr[iat]
        for shell in basis_map[sym]:
            letter = shell["letter"]
            comps = SHELL_COMPONENTS[letter]
            for contr in range(len(shell["coeffs"])):
                raw = shell["coeffs"][contr]
                exps = shell["exps"]
                for (l,m,n) in comps:
                    coeff = [raw[a]*primitive_norm(exps[a], l, m, n) for a in range(len(exps))]
                    sii = 0.0
                    for a in range(len(exps)):
                        for b in range(len(exps)):
                            sii += coeff[a]*coeff[b]*prim_overlap(l,m,n,l,m,n,
                                exps[a],exps[b], C, C)
                    if sii > 1e-15:
                        inv = 1.0/math.sqrt(sii)
                        coeff = [c*inv for c in coeff]
                    orbs.append({"lmn":(l,m,n),"exps":exps,"coeff":coeff,"center":C})
    return orbs

def read_xyz(path):
    with open(path) as f:
        nat = int(f.readline())
        name = f.readline().strip()
        atoms, coord_ang = [], []
        for _ in range(nat):
            p = f.readline().split()
            atoms.append(SYM_TO_Z[p[0]])
            coord_ang.append([float(p[1]),float(p[2]),float(p[3])])
    return atoms, coord_ang, name

def default_dat(basis):
    here = os.path.dirname(os.path.abspath(__file__))
    src = os.path.join(here, "..", "src", "core", "energy_calculators", "qm_methods")
    for d in (os.environ.get("CURCUMA_DFT_BASIS"), os.environ.get("CURCUMA_DATA"),
              src, "."):
        if not d: continue
        p = os.path.join(d, basis+".dat") if not basis.endswith(".dat") else (d if d.endswith(".dat") else os.path.join(d, os.path.basename(basis)))
        if os.path.exists(p): return p
    return os.path.join(src, basis+".dat")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xyz")
    ap.add_argument("--basis", default="def2-SVP")
    ap.add_argument("--dat", default=None)
    ap.add_argument("--charge", type=float, default=0.0)
    args = ap.parse_args()
    dat = args.dat or default_dat(args.basis)
    atoms, coord_ang, name = read_xyz(args.xyz)
    coord_bohr = [[c*AA_TO_BOHR for c in row] for row in coord_ang]
    basis_map = parse_basis_dat(dat)
    orbs = build_basis(basis_map, atoms, coord_bohr)
    n = len(orbs)
    S = mat_S(orbs)
    E, _ = build_eri(orbs)
    P = dummy_density(S)
    J, K = build_JK(E, n, P)
    e_nn = 0.0
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            r = math.sqrt(sum((coord_bohr[i][k]-coord_bohr[j][k])**2 for k in range(3)))
            if r < 1e-12: continue
            e_nn += atoms[i]*atoms[j]/r
    nelec = sum(atoms) - int(round(args.charge))
    jat = [{"z":atoms[i],"x":coord_bohr[i][0],"y":coord_bohr[i][1],"z":coord_bohr[i][2]} for i in range(len(atoms))]
    out = {"molecule":{"name":name,"natoms":len(atoms),"atoms":jat},
           "basis":args.basis,"cartesian_d":True,"nbf":n,"num_electrons":nelec,
           "nuclear_repulsion":e_nn, "eri_order":"mu_nu_lam_sig",
           "S":S, "ERI":E, "P":P, "J":J, "K":K,
           "dummy_density":"2*c*c^T, c=ones normalized so c^T S c=1 (rank-1, no eig)"}
    print(json.dumps(out))

if __name__ == "__main__":
    main()