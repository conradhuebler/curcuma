#!/usr/bin/env python3
"""Reference witness for the HF-3c gCP correction (native-port WIP).

Port of the published gCP model (H. Kruse, S. Grimme, J. Chem. Phys. 136, 154101
(2012)) in the form implemented by the reference codes:
  - dftd3/simple-dftd3  src/dftd3/gcp.f90  +  src/dftd3/gcp/param.f90
    (SPDX LGPL-3.0-or-later) -- the algorithm transcribed here,
  - grimme-lab/gcp      (the mctc-gcp driver).
The parameters are the hf/minix (i.e. HF-3c) set: sigma/alpha/beta/eta =
0.1290/1.1549/1.1763/1.1526, the per-element `emiss` and `nbas` tables, and the
Slater exponents `eta * slater_exp(Z)`.

WHAT IS VERIFIED (against /opt/orca_6_1/otool_gcp, the gCP tool ORCA ships):
  - H2 reproduces otool_gcp to 4.6e-10 Eh, which pins sigma/alpha/beta/eta/emiss
    for H and the whole 1s-1s Slater-overlap path.
  - All 10 dft_1e molecules are within 2.5e-7 Eh, six of them within 1e-8.
  - The A_n(x) = int_1^inf t^n e^{-xt} dt and B_n(x) = int_-1^1 t^n e^{-xt} dt
    auxiliaries match numerical quadrature to ~1e-15 for every n used.
  - Fortran is case-insensitive: the ZA/ZB in the source ARE za/zb, so the <1s|2s>
    norm is built from the exponents *after* the <2s|1s> swap. Getting this wrong
    is worth ~8x on the overlap.
  - B_n must use its closed form, NOT the truncated (i<=12) `bint` series, for the
    different-exponent branch; using the series there costs ~1e-4 Eh on BeH2/BH.

OPEN: the residual ~1e-7 on BeH2/BH/HF/NH3/H2O is not explained by the A/B
auxiliaries (verified) nor by the constants (H2 pins them). The leading
hypothesis is that otool_gcp v1.06 (Sep 2014) carries older per-element tables for
the heavier elements than the current param.f90 -- to be settled before the C++
port is trusted. 12-2026.

Usage:  gcp_reference_witness.py <file.xyz> [<file.xyz> ...]
"""
import math
# gCP, ported from dftd3/simple-dftd3 src/dftd3/gcp.f90 + gcp/param.f90 (LGPL-3.0-or-later).
# NOTE: Fortran is case-insensitive -- ZA/ZB in the source ARE za/zb, so the norm uses
# the exponents *after* the <2s|1s> swap.
ALPHA,BETA,ETA,SIGMA = 1.1549,1.1763,1.1526,0.1290
EMISS={1:0.04240,2:0.02832,3:0.17787,4:0.17160,5:0.22424,6:0.27995,7:0.35791,
       8:0.47901,9:0.63852,10:0.83235}
NBAS ={1:1,2:1,3:5,4:5,5:5,6:5,7:5,8:5,9:5,10:5}
SL_S=[1.2000,1.6469,0.6534,1.0365,1.3990,1.7210,2.0348,2.2399,2.5644,2.8812]
SL_P=[0.0000,0.0000,0.5305,0.8994,1.2685,1.6105,1.9398,2.0477,2.4022,2.7421]
def slater_exp(z): return SL_S[z-1] if z<=2 else 0.5*(SL_S[z-1]+SL_P[z-1])

def P(n,x):                      # coefficients n!/k!  (from the reference B0..B4)
    return sum(math.factorial(n)/math.factorial(k)*x**k for k in range(n+1))
def A(n,x):
    s=0.0; t=1.0
    for k in range(n+1):
        if k: t*=(n-k+1)
        s+=t*x**(n-k)
    return math.exp(-x)*s/x**(n+1)
def Bn(n,x):
    if abs(x)<1e-12:                    # exact limit at x=0
        return 2.0/(n+1) if n%2==0 else 0.0
    return (P(n,-x)*math.exp(x)-P(n,x)*math.exp(-x))/x**(n+1)

def ssovl(r,na,nb,za,zb):
    shell=lambda z: 1 if z<=2 else 2
    ii=shell(na)*shell(nb); R05=r*0.5
    ax=(za+zb)*R05; bx=(zb-za)*R05
    if ii==1:
        n=0.25*math.sqrt((za*zb*r*r)**3)
        return n*(A(2,ax)*Bn(0,bx)-Bn(2,bx)*A(0,ax))
    if ii==2:                            # <1s|2s>
        if shell(na)>=shell(nb): za,zb=zb,za; ax=(za+zb)*R05; bx=(zb-za)*R05
        n=math.sqrt((za**3)*(zb**5))*(r**4)*0.125
        return math.sqrt(1/3.)*n*(A(3,ax)*Bn(0,bx)-Bn(3,bx)*A(0,ax)+A(2,ax)*Bn(1,bx)-Bn(2,bx)*A(1,ax))
    if ii==4:                            # <2s|2s>
        n=math.sqrt((za*zb)**5)*(r**5)*0.0625
        return n*(A(4,ax)*Bn(0,bx)+Bn(4,bx)*A(0,ax)-2.0*A(2,ax)*Bn(2,bx))/3.0
    raise ValueError("shell pair %d outside H-Ne scope"%ii)

def gcp(atoms,coords_ang):
    A2B=1.0/0.529177210903
    xyz=[[c*A2B for c in p] for p in coords_ang]; n=len(atoms); en=[0.0]*n
    xv=[(1.0/math.sqrt(NBAS[z]-0.5*z) if NBAS[z]-0.5*z>=0.5 else 0.0) for z in atoms]
    for iat in range(n):
        for jat in range(iat+1):
            izp,jzp=atoms[iat],atoms[jat]
            emi=EMISS[izp]*xv[jat]*SIGMA; emj=EMISS[jzp]*xv[iat]*SIGMA
            r1=math.dist(xyz[iat],xyz[jat])
            if r1<1e-14 or r1>60.0: continue
            sij=ssovl(r1,izp,jzp,ETA*slater_exp(izp),ETA*slater_exp(jzp))
            dE=math.exp(-ALPHA*r1**BETA)/math.sqrt(sij)
            en[iat]+=emi*dE
            if iat!=jat: en[jat]+=emj*dE
    return sum(en), en
