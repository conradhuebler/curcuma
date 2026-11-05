#!/usr/bin/env python3
"""Claude Generated (Nov 2026). Independent COSX reference for test_cases/qm_cosx/reference.json: PySCF RHF/RKS whose K is the seminumerical exchange on a GIVEN grid
(points+weights file written by curcuma with CURCUMA_QM_COSX_DUMP), AO values from PySCF eval_gto and
the potential integrals from libcint int1e_grids. Prints exact and COSX total energies.
usage: cosx_pyscf_reference.py xyz basis method gridfile [fit 1|0] [xcgrid nrad nang]"""
import sys, json, time
import numpy as np, scipy.linalg
from pyscf import gto, scf, dft, lib

xyz, basis, method, gridfile = sys.argv[1:5]
fit = int(sys.argv[5]) if len(sys.argv) > 5 else 1
nrad = int(sys.argv[6]) if len(sys.argv) > 6 else 50
nang = int(sys.argv[7]) if len(sys.argv) > 7 else 194
lines = open(xyz).read().splitlines()
nat = int(lines[0]); atoms = "; ".join(l for l in lines[2:2 + nat])
mol = gto.M(atom=atoms, basis=basis, unit="Angstrom", verbose=0, cart=False)
g = np.loadtxt(gridfile); coords, w = g[:, :3], g[:, 3]
ao = mol.eval_gto("GTOval", coords)
S = mol.intor("int1e_ovlp")
proj = scipy.linalg.solve(ao.T @ (w[:, None] * ao), S) if fit else np.eye(mol.nao)


def cosx_k(dm):
    Pp = proj @ dm
    F = (w[:, None] * ao) @ Pp
    K = np.zeros_like(dm)
    bs = max(1, int(4e8 / 8 / mol.nao ** 2))
    for p0 in range(0, len(w), bs):
        A = mol.intor("int1e_grids", grids=coords[p0:p0 + bs])  # (ng, n, n) of 1/|r - r_g|
        G = np.einsum("gml,gl->gm", A, F[p0:p0 + bs])
        K += ao[p0:p0 + bs].T @ G
    return 0.5 * (K + K.T)


def make(m):
    if method == "hf":
        mf = scf.RHF(m)
    else:
        mf = dft.RKS(m); mf.xc = "b3lyp5"
        mf.grids.atom_grid = (nrad, nang); mf.grids.prune = None
    mf.conv_tol = 1e-12; mf.conv_tol_grad = 1e-8; mf.max_cycle = 200
    return mf

t = time.time(); mf = make(mol); e_exact = mf.kernel(); t_exact = time.time() - t
mf2 = make(mol)
base_jk = mf2.get_jk
def get_jk(mol=None, dm=None, hermi=1, with_j=True, with_k=True, omega=None):
    vj, _ = base_jk(mol, dm, hermi, with_j=True, with_k=False)
    vk = cosx_k(np.asarray(dm)) if with_k else None
    return vj, vk
mf2.get_jk = get_jk
t = time.time(); e_cosx = mf2.kernel(dm0=mf.make_rdm1()); t_cosx = time.time() - t
print(json.dumps({"exact": e_exact, "cosx": e_cosx, "diff": e_cosx - e_exact, "npts": len(w), "nao": mol.nao,
                  "converged": bool(mf.converged and mf2.converged)}))
