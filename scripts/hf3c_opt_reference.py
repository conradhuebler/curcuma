#!/usr/bin/env python3
"""Reference HF-3c minimum for test_cases/qm_grad/*_hf3c_min_reference.xyz: PySCF RHF/MINIX
+ simple-dftd3 D3(BJ) + gCP/SRB, minimised with scipy BFGS (gtol 1e-7).
Usage: hf3c_opt_reference.py <start.xyz> <out.xyz>.  Claude Generated (Sep 2026)."""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hf3c_reference import load_gamess_basis, read_xyz, MINIX, AATOAU
from pyscf import gto, scf
from dftd3.interface import GeometricCounterpoise, RationalDampingParam, DispersionModel
from scipy.optimize import minimize
B=load_gamess_basis(MINIX)
atoms=read_xyz(sys.argv[1]); sym=[s for s,_ in atoms]; x0=np.array([x for _,x in atoms]).ravel()*AATOAU
num=np.array([gto.charge(s) for s in sym])
dm=[None]
def f(x):
    pos=x.reshape(-1,3)
    mol=gto.M(atom=[[s,p] for s,p in zip(sym,pos)],basis={s:B[s] for s in sym},unit='Bohr',verbose=0)
    mf=scf.RHF(mol); mf.conv_tol=1e-12; e=mf.kernel(dm0=dm[0]); dm[0]=mf.make_rdm1()
    g=mf.nuc_grad_method().kernel()
    d=DispersionModel(num,pos).get_dispersion(RationalDampingParam(method="hf3c"),grad=True)
    c=GeometricCounterpoise(num,pos,method="hf3c").get_counterpoise(grad=True)
    return e+float(d['energy'])+float(c['energy']), (g+d['gradient']+c['gradient']).ravel()
r=minimize(f,x0,jac=True,method='BFGS',options={'gtol':1e-7,'maxiter':500})
print("E_min %.10f  |g|max %.2e  iters %d"%(r.fun,np.abs(r.jac).max(),r.nit))
with open(sys.argv[2],'w') as o:
    o.write("%d\nreference HF-3c minimum E=%.10f\n"%(len(sym),r.fun))
    for s,p in zip(sym,r.x.reshape(-1,3)/AATOAU): o.write("%s %.8f %.8f %.8f\n"%(s,*p))
