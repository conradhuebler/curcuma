#!/usr/bin/env python3
"""
Independent Python witness for the curcuma native-DFT 1e integrals (WP1).

Reads an XYZ and the SAME def2-SVP.dat curcuma uses, builds the contracted
cartesian GTO basis in the SAME AO order as curcuma's
BasisSetParser::createGTOFromBasis (per atom in input order; per shell in file
order; p -> [px,py,pz]; d -> [dxx,dyy,dzz,dxy,dxz,dyz]), and computes the
overlap S, kinetic T (gradient identity) and nuclear-attraction V (McMurchie-
Davidson + Boys) via an independently written Obara-Saika implementation.

Output JSON schema (matches dump_dft_1e.cpp):
  { "molecule":{...}, "basis","cartesian_d":true, "nbf","num_electrons",
    "nuclear_repulsion", "S":[[...]],"T":[[...]],"V":[[...]],"H":[[...]] }

This is the element-wise kernel witness: because curcuma and this script share
the exact AO ordering (cartesian 6d), the matrices can be compared element by
element. ORCA (which only exposes MOs, not the AO integral matrices, via
orca_2json) is validated separately via the generalized-eigenvalue spectrum
(see diff_dft_1e.py).

  dft_1e_python_ints.py <input.xyz> [--basis NAME] [--dat FILE] [--charge Q]

Claude Generated (WP1). GPL-3.0.
"""
import argparse, json, math, os, sys

PI = 3.14159265358979323846264338327950288
AA_TO_BOHR = 1.0 / 0.529177210903

# Element symbol -> Z (WP1 scope H-Ne).
SYM_TO_Z = {"H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10}
Z_TO_SYM = {v:k for k,v in SYM_TO_Z.items()}

# Cartesian (l,m,n) components per shell-letter, in curcuma emission order.
SHELL_COMPONENTS = {
    "S": [(0,0,0)],
    "P": [(1,0,0),(0,1,0),(0,0,1)],
    "D": [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)],
}

def dfact(n):
    r = 1.0
    k = 1
    while k <= n:
        r *= (2*k-1); k += 1
    return r if n > 0 else 1.0  # (-1)!! = 1, (-2)!! = 1

def primitive_norm(alpha, l, m, n):
    L = l+m+n
    denom = math.sqrt(dfact(l)*dfact(m)*dfact(n)) or 1.0
    return (2.0*alpha/PI)**0.75 * (4.0*alpha)**(L/2.0) / denom

def overlap1d_table(imax, jmax, PA, PB, gamma):
    """OS 1D overlap recurrence S[i][j]."""
    g2 = 0.5/gamma
    S = [[0.0]*(jmax+1) for _ in range(imax+1)]
    S[0][0] = 1.0
    for i in range(1, imax+1):
        prev = S[i-2][0] if i >= 2 else 0.0
        S[i][0] = PA*S[i-1][0] + g2*(i-1)*prev
    for i in range(imax+1):
        for j in range(1, jmax+1):
            si1 = S[i-1][j-1] if i >= 1 else 0.0   # i*S(i-1,j-1)
            sj1 = S[i][j-2]   if j >= 2 else 0.0   # (j-1)*S(i,j-2)
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

def prim_overlap_shifted(la,ma,na, lb,mb,nb, gamma, K, PAx,PBx, PAy,PBy, PAz,PBz):
    """Overlap with arbitrary angular momenta (used by the kinetic kernel)."""
    if K == 0.0: return 0.0
    pref = (PI/gamma)**1.5 * K
    Sx = overlap1d_table(max(la,lb), max(la,lb), PAx, PBx, gamma)
    Sy = overlap1d_table(max(ma,mb), max(ma,mb), PAy, PBy, gamma)
    Sz = overlap1d_table(max(na,nb), max(na,nb), PAz, PBz, gamma)
    return pref * Sx[la][lb] * Sy[ma][mb] * Sz[na][nb]

def hermite_coeffs(iA, iB, PA, PB, gamma, K):
    """E[t][i][j] Hermite expansion coefficients (K folded into E[0][0][0]).

    Standard McMurchie-Davidson forward recursion (Helgaker 9.5.5/9.5.6): raise
    the polynomial powers i then j, filling every Hermite order t at once.
      E^{i+1,j}_t = (1/(2p)) E^{i,j}_{t-1} + PA E^{i,j}_t + (t+1) E^{i,j}_{t+1}
      E^{i,j+1}_t = (1/(2p)) E^{i,j}_{t-1} + PB E^{i,j}_t + (t+1) E^{i,j}_{t+1}
    (The earlier "raise t" variant matched only t<=1: E[2][0][0] came out 0.5
    instead of 0.0 and E[2][1][1] 0.25 instead of 0.0625.)
    """
    tmax = iA+iB
    g2 = 0.5/gamma
    E = [[[0.0]*(iB+1) for _ in range(iA+1)] for _ in range(tmax+1)]
    if K == 0.0: return E
    E[0][0][0] = K

    def at(t, i, j):
        if t < 0 or t > tmax or i < 0 or j < 0 or i > iA or j > iB: return 0.0
        return E[t][i][j]
    for i in range(1, iA+1):
        for t in range(tmax+1):
            E[t][i][0] = g2*at(t-1, i-1, 0) + PA*at(t, i-1, 0) + (t+1)*at(t+1, i-1, 0)
    for j in range(1, iB+1):
        for i in range(iA+1):
            for t in range(tmax+1):
                E[t][i][j] = g2*at(t-1, i, j-1) + PB*at(t, i, j-1) + (t+1)*at(t+1, i, j-1)
    return E

def boys_array(maxN, T):
    F = [0.0]*(maxN+1)
    if T < 1e-14:
        for n in range(maxN+1): F[n] = 1.0/(2*n+1)
        return F
    # Large T: the downward recursion needs a start index far above T; a fixed
    # maxN+25 start collapses for T >~ 15 (F_0(37) 50x too small, F_0(T>=50)
    # exactly 0). Use the closed form F_0(T) = 0.5 sqrt(pi/T) erf(sqrt(T)) and
    # recur UPWARD -- the stable direction once e^{-T} is small (its subtraction
    # cancels at small T, hence the split at T = 1).
    if T >= 1.0:
        F[0] = 0.5*math.sqrt(math.pi/T)*math.erf(math.sqrt(T))
        eT = math.exp(-T)
        for n in range(maxN):
            F[n+1] = ((2*n+1)*F[n] - eT)/(2*T)
        return F
    M = maxN + 25
    G = [0.0]*(M+1)
    eT = math.exp(-T)
    for n in range(M-1, -1, -1):
        G[n] = (2.0*T*G[n+1] + eT)/(2*n+1)
    for n in range(maxN+1): F[n] = G[n]
    return F

def build_R(tmax, umax, vmax, Nmax, PC, gamma, boys):
    """R[t][u][v][N] Coulomb auxiliary. Returns nested lists."""
    Px,Py,Pz = PC
    R = [[[[0.0]*(Nmax+1) for _ in range(vmax+1)] for _ in range(umax+1)] for _ in range(tmax+1)]
    pref = 2.0*PI/gamma
    # base R^N_{000} = (2 pi / gamma) (-2 gamma)^N F_N(T): the (-2 gamma)^N factor is
    # required by the t >= 2 recurrence (whose R^{n+1} carries one more (-2 gamma)).
    # It was missing, so R^0_{200} came out +(2 pi/gamma) F_1 instead of
    # -(2 pi/gamma) 2 gamma F_1 -- wrong magnitude AND sign.
    neg2p = 1.0
    for N in range(Nmax+1):
        R[0][0][0][N] = pref*neg2p*boys[N]
        neg2p *= -2.0*gamma
    # Displacement term carries a PLUS sign: this is the bra recurrence
    # (differentiation w.r.t. P), same form as the ERI's build_R_eri. The former
    # minus reproduced every on-centre integral (P-C == 0) but inverted the sign of
    # off-centre elements.
    for t in range(tmax+1):
        if t >= 1:
            for N in range(Nmax):
                low = (t-1)*R[t-2][0][0][N+1] if t >= 2 else 0.0
                R[t][0][0][N] = low + Px*R[t-1][0][0][N+1]
        for u in range(1, umax+1):
            for N in range(Nmax):
                low = (u-1)*R[t][u-2][0][N+1] if u >= 2 else 0.0
                R[t][u][0][N] = low + Py*R[t][u-1][0][N+1]
        for u in range(umax+1):
            for v in range(1, vmax+1):
                for N in range(Nmax):
                    low = (v-1)*R[t][u][v-2][N+1] if v >= 2 else 0.0
                    R[t][u][v][N] = low + Pz*R[t][u][v-1][N+1]
    return R

def prim_nuclear(la,ma,na, lb,mb,nb, a1,a2, A, B, C):
    gamma, P, K = gaussian_product(a1,a2,A,B)
    if K == 0.0: return 0.0
    Ex = hermite_coeffs(la, lb, P[0]-A[0], P[0]-B[0], gamma, K)
    Eu = hermite_coeffs(ma, mb, P[1]-A[1], P[1]-B[1], gamma, 1.0)
    Ev = hermite_coeffs(na, nb, P[2]-A[2], P[2]-B[2], gamma, 1.0)
    tmax,umax,vmax = la+lb, ma+mb, na+nb
    Nmax = tmax+umax+vmax
    T = gamma*sum((P[k]-C[k])**2 for k in range(3))
    F = boys_array(Nmax, T)
    R = build_R(tmax,umax,vmax,Nmax, [P[k]-C[k] for k in range(3)], gamma, F)
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
                val += et*eu*ev*R[t][u][v][0]
    return val

def parse_basis_dat(path):
    """Parse the curcuma $DATA basis format -> {symbol: [shells]}."""
    with open(path) as f:
        lines = f.readlines()
    name_to_sym = {"HYDROGEN":"H","HELIUM":"He","LITHIUM":"Li","BERYLLIUM":"Be",
                   "BORON":"B","CARBON":"C","NITROGEN":"N","OXYGEN":"O",
                   "FLUORINE":"F","NEON":"Ne"}
    in_data = False
    i = 0
    basis = {}
    n = len(lines)
    # skip to $DATA
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
        # element header?
        if first.isalpha() and first.upper() in name_to_sym and (len(toks)==1 or not toks[1][0].isdigit()):
            cur = name_to_sym[first.upper()]
            basis[cur] = []
            continue
        # shell line: <S|P|D> nprim [ncontr]
        letter = first[0].upper()
        nprim = int(toks[1])
        ncontr = int(toks[2]) if len(toks) > 2 else 1
        exps = []
        coeffs = [[] for _ in range(ncontr)]
        for _ in range(nprim):
            p = lines[i].split()
            i += 1
            # p = [idx, exp, c0, c1, ...]
            exps.append(float(p[1]))
            for c in range(ncontr):
                coeffs[c].append(float(p[2+c]))
        basis[cur].append({"letter":letter,"exps":exps,"coeffs":coeffs})
    return basis

def build_basis(basis_map, atoms, coord_bohr):
    """Build flat contracted cartesian basis in curcuma AO order, normalized."""
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
                    # renormalize contracted shell so S_ii = 1
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

def mat_T(orbs):
    n = len(orbs)
    T = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            a, b = orbs[i], orbs[j]
            A, B = a["center"], b["center"]
            v = 0.0
            for ia in range(len(a["exps"])):
                al = a["exps"][ia]
                for ib in range(len(b["exps"])):
                    be = b["exps"][ib]
                    gamma, P, K = gaussian_product(al, be, A, B)
                    if K == 0.0: continue
                    PAx,PAy,PAz = P[0]-A[0],P[1]-A[1],P[2]-A[2]
                    PBx,PBy,PBz = P[0]-B[0],P[1]-B[1],P[2]-B[2]
                    la,ma,na = a["lmn"]; lb,mb,nb = b["lmn"]
                    def ov(LA,MA,NA,LB,MB,NB):
                        return prim_overlap_shifted(LA,MA,NA,LB,MB,NB,gamma,K,PAx,PBx,PAy,PBy,PAz,PBz)
                    tot = 0.0
                    # x
                    tx = 0.0
                    if la>0 and lb>0: tx += la*lb*ov(la-1,ma,na,lb-1,mb,nb)
                    if la>0:         tx -= 2*be*la*ov(la-1,ma,na,lb+1,mb,nb)
                    if lb>0:         tx -= 2*al*lb*ov(la+1,ma,na,lb-1,mb,nb)
                    tx += 4*al*be*ov(la+1,ma,na,lb+1,mb,nb)
                    # y
                    ty = 0.0
                    if ma>0 and mb>0: ty += ma*mb*ov(la,ma-1,na,lb,mb-1,nb)
                    if ma>0:         ty -= 2*be*ma*ov(la,ma-1,na,lb,mb+1,nb)
                    if mb>0:         ty -= 2*al*mb*ov(la,ma+1,na,lb,mb-1,nb)
                    ty += 4*al*be*ov(la,ma+1,na,lb,mb+1,nb)
                    # z
                    tz = 0.0
                    if na>0 and nb>0: tz += na*nb*ov(la,ma,na-1,lb,mb,nb-1)
                    if na>0:         tz -= 2*be*na*ov(la,ma,na-1,lb,mb,nb+1)
                    if nb>0:         tz -= 2*al*nb*ov(la,ma,na+1,lb,mb,nb-1)
                    tz += 4*al*be*ov(la,ma,na+1,lb,mb,nb+1)
                    tot = tx+ty+tz
                    v += 0.5*a["coeff"][ia]*b["coeff"][ib]*tot
            T[i][j] = T[j][i] = v
    return T

def mat_V(orbs, atoms, coord_bohr):
    n = len(orbs)
    V = [[0.0]*n for _ in range(n)]
    for C in range(len(atoms)):
        Z = atoms[C]
        if Z == 0: continue
        Cpos = coord_bohr[C]
        for i in range(n):
            for j in range(i, n):
                a, b = orbs[i], orbs[j]
                val = 0.0
                for ia in range(len(a["exps"])):
                    for ib in range(len(b["exps"])):
                        val += a["coeff"][ia]*b["coeff"][ib]*prim_nuclear(
                            *a["lmn"], *b["lmn"], a["exps"][ia], b["exps"][ib],
                            a["center"], b["center"], Cpos)
                val = -Z*val
                V[i][j] += val
                if i != j: V[j][i] += val
    return V

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
    S = mat_S(orbs); T = mat_T(orbs); V = mat_V(orbs, atoms, coord_bohr)
    H = [[T[i][j]+V[i][j] for j in range(n)] for i in range(n)]
    # nuclear repulsion
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
           "nuclear_repulsion":e_nn,
           "S":S,"T":T,"V":V,"H":H}
    print(json.dumps(out))

if __name__ == "__main__":
    main()