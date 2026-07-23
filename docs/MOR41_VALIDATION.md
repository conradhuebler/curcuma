# MOR41: Native GFN-FF / GFN1 / GFN2 vs xtb 6.6.1 — Validation

Status: AI-generated, machine-tested (Jul 2026). NOT human production-tested.
The AI does not assign ✅ TESTED / ✅ APPROVED.

## Update (Jul 16, 2026): native GFN1 + GFN2 transition metals FIXED

The original negative result below was traced to concrete bugs; **native GFN1 and GFN2
now reproduce tblite for transition metals** (residual ~1e-3 Eh, down from up to 15 Eh).
Three bugs, all in `xtb_native.cpp` / `STO_CGTO.hpp`:

1. **Shell-vs-angular indexing** — `reference_occ`, `p_kcn`, `p_shpoly` are indexed by
   **angular momentum** in tblite (`gfn2.f90` / `gfn1.f90`: `(0:2, elem)`), but were read
   by **shell index**. Main-group is unaffected (`ang_sh==ish`); transition metals order
   shells `[d,s,p]`, so shell-indexing scrambled the d/s/p occupations and CN/shpoly →
   wrong `n0`/`dq` → multi-Eh error. Fix: index by the shell's angular momentum.
   - GFN2: all three by angular momentum (no shell shares an l).
   - GFN1: valence-aware — `reference_occ` goes to the **valence** shell of each l
     (first shell of that l), polarisation shells get 0 (H has two s-shells); `p_shpoly`
     by angular momentum for all shells; **`p_kcn` stays shell-indexed** (tblite gfn1
     keeps it `max_shell`).
2. **6s/6p STO-NG** — the Stewart STO-NG tables in `STO_CGTO.hpp` stopped at n=5, so the
   6s/6p shells of 5d metals (W/Ir/Pt) used the wrong Gaussians (6s→2p, 6p→3d slots) →
   ~0.1–0.5 Eh per 5d atom. Added tblite's dedicated 6s/6p STO-6G arrays.

Result (native vs tblite, identical geometries):
- **GFN2**: 3d metals (Cr/Mn/Fe/Co/Ni) **exact (1e-8)**; 4d/5d **~1e-3 Eh**. MOR41 reaction
  MAD native-GFN2-vs-xtb **507 → 0.14 kcal/mol** (max 5735 → 1.24). curcuma-GFN2 vs DLPNO
  MAD ~15, comparable to xtb-GFN2's own 11.8.
- **GFN1**: **72/95 structures exact (1e-8), 88/95 < 1e-3 Eh** (were off by 100s of
  kcal/mol); 5d (Ir PR22) 4e-8 after the 6s/6p fix.
- **No regression**: all main-group sqm molecules (incl. d-shell H2S/PH3/SiH4/HCl and H)
  stay 1e-8 vs tblite for both GFN1 and GFN2; gfn2 ecomp ctests pass.

The remaining residual was then fully localized to **dispersion** (overlap and all electronic
params are bit-exact vs tblite) and closed by four more fixes: the **5p STO-4G transcription
error** (period-5 overlap), the **D4 `r4/r2` table truncation at Z≤36** (heavy-element D4
under-binding), the **D3/D4 CN-Gaussian weight-normalization threshold** (collapsed C6 for
sparse-reference-CN metals → GFN1 D3 under-binding), and the **missing D4 `sscale` entry for
reference-system Na** (refsys=11) which left the 4th high-CN reference of Sc/Ti/V/Zr/Nb/Hf/Ta
with C6 ~2.4× too large (the gfn2 PR40/Ti over-binding).

**Final result vs tblite: native GFN1 and GFN2 reproduce the reference for ALL 95 MOR41
structures to <1e-6 Eh (95/95 both methods; ~85/95 at the 1e-8 print-precision floor).**
MOR41 gfn2 reaction MAD vs xtb: **507 → 0.03 kcal/mol.** No main-group/3d/GFN-FF regression;
all energy-validation ctests pass.

The only remaining transition-metal gap is **GFN-FF** (a force field, separate from the GFN1/2
tight-binding stack). Its metal **bond term** was ~1.95× too strong because the Fortran
metal-bond branch (`btyp>=5`) was not implemented; **that branch has since been implemented**
(`53d6aeb`, refined by `14fc648`) and, together with the four-list neighbour port (`5afbb02`),
brings MOR41 GFN-FF per-structure MAD to **7.30** kcal/mol (max **37.80** = ED07, 39/95 within
1 kcal/mol). See [GFNFF_METAL_BOND_ANALYSIS.md](GFNFF_METAL_BOND_ANALYSIS.md) — note that file
is the **pre-fix** diagnostic and is flagged as superseded.

**Two GFN-FF metal fixes (Jul 2026)** took the per-structure MAD from 7.295 to **5.288**
(max 37.80 → **35.08**, within-1 39 → **45/95**), with **zero structures net worse**:

1. **ED07 bond-list bug** (−37.80 → +6.34): the bond-term generators re-derived a 70-bond list
   from an old distance heuristic while the ported getnb criterion (and xtb) said 68. Both
   generators now consume the getnb list. See [GFNFF_NEIGHBOR_LISTS.md](GFNFF_NEIGHBOR_LISTS.md).
2. **Metal-H `mtyp` + `fsrb2`** (fixed a whole class of hydride complexes: PR07 −27→+3.5,
   PR09 −33→−2.1, PR08 −24→−0.9, ED17, PR34/35/33/17, PR06, ED14, PR10): hydrogen was wrongly
   excluded from `mtyp=1` (it is group 1), which mis-set `fqq`/`fcn` on metal-H bonds; and the
   `fsrb2` EN-scaling was keyed on `mtyp>0` instead of an actual metal bond, so fixing `mtyp(H)`
   also required gating `fsrb2` on `imetal>0` to keep organic C-H/O-H bonds unchanged. Both
   ground-truthed per-bond against the in-tree Fortran analyzer.

**eta (`btyp=6`) + TM-TM `mchar` wired (Jul 2026):** eta bonds now promote to `bstren[6]=0.78`
+ the `eta_shift*nb(20,metal)` r0 term (detection matches the Fortran analyzer bond-for-bond:
PR15 4/4, ED13 7/7, PR13 9/9); TM-TM bonds get `bstrength *= 1-min(2·mchar_i+2·mchar_j, 0.5)`
(PR22 Ir-Ir prefactor now exact vs Fortran). Per-structure MAD 5.288 → **4.658**, within-1 45 →
**48/95**; per-reaction MAD → **4.067**.

**Phase-1 π-metal `qa` charge fixed (Jul 2026):** the eta `fqq` error (PR15 1.0235 vs analyzer
1.7496) was a **directed-graph bug** in the sparse topological-distance builder. `neighbor_lists`
(== Fortran `nbdum`) is asymmetric for eta bonds — the metal lists the eta ligand but the ligand
omits the metal (`nbm`, gfnff_ini2.f90:199) — and the Dijkstra adjacency was directed, so the eta
ligand was disconnected from its metal in the Phase-1 graph. Fortran (and this file's
Floyd-Warshall variant) symmetrize; the sparse variant now does too. PR15 `topo%qa` now matches
the analyzer bit-for-bit (Ni 0.289140, eta-C −0.070091), eta `fqq` → 1.749625, PR15 +40.18 →
+0.09. Per-structure MAD → **4.226**, within-1 → **49/95**; per-reaction MAD → **3.470**, max
31.83 → **17.54**. Non-eta systems byte-identical.

Session arc for native GFN-FF on MOR41 (per-reaction MAD): 57 (start) → 8.17 (bond-list) →
4.07 (metal-H + eta/mchar) → **3.47** (π-metal qa). Remaining residuals are no longer dominated
by any single mechanism identified so far.

## Update (Jul 23, 2026): FT-HMO π-occupation fixes + the reference-implementation split

Status: AI-generated, machine-tested. NOT human production-tested.

### Two Hückel (FT-HMO) occupation fixes — `huckel_solver.cpp`

The metal-coordinated aromatic rings (Cp / N-heteroaromatic ligands) had **odd** π-electron
counts that curcuma's closed-shell FT-HMO occupation handled wrongly. Two faithful ports of the
Fortran `gfnff_qm.f90` / `gfnff_ini.f90` closed the gap:

1. **Open-shell `occu` split for odd nelpi** (was: `nel_alpha=nel/2`, one smear doubled — placed
   only `nel-1` electrons and ran the biradical test at the wrong orbital). Now mirrors
   `scc_core.f90 occu`: `ihomoa=nel/2+1`, `ihomob=nel/2`, two separate `fermismear` passes,
   biradical test at the true alpha HOMO. Even nelpi is bit-identical to before (benzene / COT /
   real antiaromatics untouched). Fixes the 5-electron Cp rings (ED36, PR04, PR10).
2. **`pisip>0.40` "wrong pi occupation" fallback** (`gfnff_ini.f90:1082-1099`, xtb variant — NOT
   gated on the print flag, unlike the pprcht standalone which nests it in `if(pr2)`; redo at
   et=4000). If the converged HOMO energy > 0.40, redo the final solve with `nelpi-1`. Fixes the
   5-atom/7-electron N-heteroaromatic rings (ED21/PR16/ED16a), which without it put all 7
   electrons into an antibonding orbital → +180 kcal too weak.

Effect vs the port reference (**pprcht/gfnff**, standalone, all 95 structures):
per-structure MAD **7.27 → 1.61** (max 82.8 → 20.9, within-1 52 → 63/95). 71/71 runnable gfnff
ctests still pass; 64/95 MOR41 structures bit-identical (only odd-nelpi π-systems change).

### Two carbene fixes — the central amidine C (ED16b, `gfnff_method.cpp`)

Ground-truthed per-bond/per-angle against the standalone analyzer, ED16b (an organic
N-heterocyclic diamine, NO metal) had a −20.9 kcal residual = bond −12 + angle −9, both from its
2-coordinate central amidine carbon (N-C-N, a carbene):

1. **itag reaches the FT-HMO π-count** — the `GFNFF::calculatePiBondOrders` call site hard-coded
   an all-zeros `itag` (`gfnff_method.cpp:8065`), so the carbene C (correctly tagged `itag=1` by
   `determineHybridizationFortran`) was counted as +1 π-electron → nelpi 7 not 6 → wrong pibo +
   a spurious `pisip>0.40` fallback. Now passes `topo_info.itag` through. Bond term → 1e-7 Eh.
2. **carbene angle equilibrium** (`gfnff_ini.f90:1573-1577`) — a 2-coordinate group-4 center with
   `itag=1` uses θ0=145° (C) / 90° (heavier), overriding the hyb/ring rule. curcuma lacked it, so
   the carbene kept the 5-ring θ0=109° → too little angle strain. Added before the metal-center
   branch; `fbsmall`/force-constant follow from θ0 automatically.

**ED16b: −20.85 → 0.0000 kcal** vs pprcht **and** xtb (they agree here). Only ED16b changes in
MOR41 (surgical, no regression). Overall per-structure MAD vs pprcht **1.61 → 1.39**
(max 20.9 → 12.7, within-1 63 → 64/95).

### Halogen-bond B-atom topological-distance filter (`gfnff_method.cpp:9310`)

`detectHalogenBondsNative` only excluded a base B directly bonded to the X donor (topological
distance 1), but Fortran `gfnff_ini.f90:872` requires `bpair(B,X) > 3` — B must be **more than 3
bonds** from X (A…B, not X-B). So B atoms 2-3 bonds from X — e.g. the phosphine P and ring C
reachable through the Ru center from an S donor in PR34 — were wrongly admitted: **114 X-bonds,
−0.0201 Eh vs the reference 22 / −0.0018 Eh** (S is a valid chalcogen X donor, so it is not fully
spurious, just over-counted ~11×). Now filters on `topo_info.bpair[X][B] <= 3` (999 = beyond BFS
depth / other fragment, correctly kept). PR34 **−10.48 → +0.99 kcal**.

Broad effect (31 structures touched, all with S/P/metal X donors): per-structure MAD vs pprcht
**1.39 → 0.975**, within-1 64 → **71/95**; 71/71 runnable gfnff ctests pass. Verified faithful —
curcuma's X-bond energy now matches the analyzer bit-for-bit (PR34 −0.00179 vs −0.00183; ED18/PR23
both exactly 0). Several structures that *appeared* to worsen (ED18 −0.8→+3.2, PR07, ED07) had a
spurious negative X-bond masking a real residual in another term — now correctly exposed for the
next pass.

### EEQ `gam`/`chi` heavy-element array corruption (`gfnff_par.h`)

ED07 (the previously "unexplained" W complex) traced to a **wrong EEQ hardness**: `gam_eeq[73]`
(W) was **+0.064240** vs the Fortran `gam_angewChem2020(74) = −0.003724`. A full array diff showed
**`gam_eeq` was corrupt for all of Z=56-86** (Ba + lanthanides + 5d/6p — placeholder/interpolated
values; the Z≥87 slots even repeated the Z=55-56 chunk), and **`chi_eeq` for Z=57-71** (La-Lu). Z≤55
and cnf/alp were correct. Replaced Z=56-86 `gam` and Z=57-71 `chi` with the verbatim Fortran
`*_angewChem2020` arrays.

Wrong W hardness → wrong W EEQ charge (**qa 0.326 → 0.351**, the reference value) → wrong Coulomb
(+2.5 kcal) and metal-bond `fqq` (W-P 1.266 → 1.366, bond +6.5 kcal). Fixing it resolved every 5d-metal
(W/Ir/Pt) structure at once: **ED07 +8.54 → +0.09, PR07 +6.65 → +0.10, ED18 +3.15 → +0.17, PR22
+2.88 → +0.23, ED22 +2.15 → 0.00, PR31, PR08, ED08, ED29, …** (14 improved, 1 marginal). 3d/4d
metals (Z≤55) unchanged. Per-structure MAD vs pprcht **0.975 → 0.630**, within-1 71 → **82/95**;
71/71 runnable gfnff ctests pass. (No MOR41 lanthanides, so the `chi` Z=57-71 fix is correctness-only.)

**Session arc vs pprcht/gfnff (per-structure MAD):** 7.27 (start) → 1.61 (FT-HMO occu + pisip) →
1.39 (carbene itag + angle) → 0.975 (halogen bpair) → 0.630 (EEQ gam array) → 0.543
(eta-aware X-bond bpair) → 0.469 (SP3-specials torsion order) → **0.29** (bond fcn neighbour
count, PR27), within-1 52 → 86/95. **Reaction-level MAD vs pprcht: 0.255 kcal/mol** (max 1.46 =
rxn 23, the `.CHRG` quirk where curcuma is correct) — the per-structure residuals largely cancel
in the reaction energies, so the port reproduces its reference at the level that matters.

### SP3-specials torsion order — ED30/PR30 (`gfnff_torsions.cpp`)

ED30 (an organic aminophosphine, no metal) and PR30 (its Pt complex) had a ~−3.5 kcal residual
entirely in the torsion term (curcuma ~half the reference). The aminophosphine P-N torsions got
curcuma's phi0=180/f1~0.006 default instead of the reference phi0=60/f1=3.0. Root cause: the
Fortran **SP3-specials** block (`gfnff_ini.f90:1746`) — a bond between two RAW-sp3 group-5 atoms
(N-N, P-P, N-P) → nrot=3, phi0=60, f1=3.0, using the raw hyb — sits AFTER the pi-sp3 case (:1733)
and **overrides** it. curcuma applied its equivalent SP3-specials BEFORE its pi-sp3 override (which
was relocated to the end of `getGFNFFTorsionParameters`), so for a P-N bond with the ring N in a
pi-system the pi-sp3 override (`!in_ring && k_in_pi && !j_in_pi && hyb_j==3`) reset it back to
phi0=180/f1=0.2. Moved the raw-hyb SP3-specials override to just before the force-constant assembly
(after the pi-sp3 block), matching the Fortran order; the pi f2 0.55 scaling is re-applied for the
rare conjugated (pibo>0) case. **ED30 −3.44 → +0.00, PR30 −3.72 → +0.08** (only these two change).
Per-structure MAD **0.543 → 0.469**, within-1 83 → **85/95**; 71/71 runnable gfnff ctests pass.

### Bond fcn neighbour count — PR27 (`getGFNFFBondParameters`, `gfnff_method.cpp:~4754`)

PR27 (a CpRu(PR3)2Cl complex) had a **+17.7 kcal** per-structure residual, an order of magnitude
above everything else and NOT a metal-bond fine-precision case (`vs_ppr != vs_xtb`). Term-by-term
against the analyzer, the entire gap was in **one bond**: the spurious cis P…P bond (R=2.645 Å,
both phosphines on the Ru) had curcuma `fc=-0.0002` vs reference `-0.046`, losing **−18.5 kcal**.

Root cause: the heavy-heavy bond weakener `fcn = 1/(1+0.007·nb(20,i)²)/(1+0.007·nb(20,j)²)`
(`gfnff_ini.f90:1181-1183`) uses `topo%nb(20,i)` — the **bonded-neighbour COUNT** (slot 20 of the
`nb` array holds the degree, a Fortran convention; `nb(1:19,i)` are the neighbour indices). curcuma
called `countNeighborsWithin20Bohr()`, a literal 20-Bohr distance sphere (~40 atoms in a compact
48-atom complex), collapsing fcn to ~0.007 and nuking every **non-metal** heavy-heavy bond (P-P,
P-S, S-S, …). Metal bonds were unaffected — the metal branch (`gfnff_ini.f90:1254-1259`) already
used the bonded degree via `topo.neighbor_lists`. Now the normal path uses the same
`topo.neighbor_lists[atom].size()`.

**PR27 +17.67 → −0.76**, P-P bond E −0.00015 → −0.02961 Eh (ref −0.02963). Only PR27 changes in
MOR41 (sole non-metal heavy-heavy bond); 35/35 golden-value gfnff ctests byte-identical, 71/71
runnable gfnff ctests pass. Per-structure MAD **0.469 → 0.29**, within-1 85 → **86/95**. Note this
is a general correctness fix — it affects any molecule with a non-metal heavy-heavy bond (disulfides,
phosphine pairs, P-S, …), not just PR27.

### Eta-aware X-bond bpair — ED33 (`detectHalogenBondsNative`)

ED33's residual was the halogen filter dropping a valid X-bond (curcuma X-bond 0 vs the reference
−0.0101 Eh). Root cause via the instrumented analyzer: the reference accepts A=Ru, X=P, B=eta-C at
`bpair=5`, but curcuma computes `bpair(P,C)=2`. The two ring carbons are **eta-coordinated** to Ru;
the Fortran `bpair` (nbondmat/pairsbond) records a distance only for **symmetric** reachability
(`dai .and. daj`), and eta bonds are stored asymmetrically (the metal lists the eta-C, the eta-C
omits the metal — `gfnff_ini2.f90:199`), so an eta bond can never bridge a ≥2-bond path. curcuma's
`bpair` is a plain BFS that bridges through the eta Ru-C, wrongly shortcutting P…C_eta to 2. The
X-bond filter now rebuilds the distance on an **eta-free adjacency** (metal↔`itag=−1` edges removed;
normal Ru-P/Ru-S kept), so eta bonds no longer bridge. **ED33 +6.76 → +0.45** (the residual is now a
smaller bond/BATM term); **PR34 unchanged** (S donor, no eta); other eta structures (ED36/PR13/PR15)
byte-identical. Per-structure MAD **0.630 → 0.543**, within-1 82 → **83/95**; 71/71 runnable gfnff
ctests pass.

### The reference-implementation split (OPEN PORTING QUESTION)

Building the standalone **pprcht/gfnff** analyzer (`external/gfnff`, `-Dbuild_exe=ON`, the code
curcuma was ported from — "adapted from the xtb code", Spicher/Grimme 2020) revealed that the
**two GFN-FF implementations themselves disagree** on transition-metal complexes:

| comparison (per-structure, n=95) | MAD | max |
|----------------------------------|----:|----:|
| curcuma (both fixes) vs **pprcht/gfnff** (the port reference) | **1.61** | 20.9 |
| **pprcht/gfnff** vs **xtb 6.7.1** (two GFN-FF impls) | **11.6** | **178** |
| curcuma vs xtb 6.7.1 | 12.7 | 178 |

So the large "curcuma vs xtb" residual on TMs (PR10 −178, PR04 −164, ED10 −111, PR28 −103, …) is
**overwhelmingly the pprcht-vs-xtb divergence, not a curcuma port error** — curcuma now tracks its
reference (pprcht) to MAD 1.6. xtb 6.7.1's GFN-FF has evidently changed on TMs since the pprcht
snapshot was extracted (cf. Moradi et al., *J. Comput. Chem.* 2026, "Extensions to Extended
Tight-Binding Methods for Transition-Metal Containing Systems").

**RESOLVED against DLPNO-CCSD(T) (Jul 23, 2026).** The 41 reaction energies were assembled for all
three engines and compared to the DLPNO-CCSD(T)/CBS reference (`reactions.dat` col 2, Table S1;
`scripts`-style driver in the session scratchpad):

| engine (reaction dE vs DLPNO, n=41) | MAD | RMSD | max | MD |
|-------------------------------------|----:|-----:|----:|----:|
| **pprcht/gfnff** | **62.6** | 80.7 | 279 | +56.4 |
| curcuma native (all fixes) | 63.1 | 81.0 | 279 | +56.9 |
| xtb 6.7.1 `--gfnff` | 71.7 | 99.9 | 419 | +61.4 |

Two conclusions:
1. **GFN-FF is fundamentally unsuitable for MOR41 reaction thermochemistry** — all three impls are
   ~60–70 kcal/mol from DLPNO (systematic +56 MD: they underestimate the exothermicity of the
   metal-ligand bond formation/breaking these reactions involve). This is a **force-field
   limitation, not a port bug** — a FF has no electronic bond-dissociation energy. (Tight-binding
   GFN2 reaches ~12 MAD on the same set; GFN-FF is not built for this property.)
2. **pprcht is closer to the true QM than xtb** (MAD 62.6 vs 71.7; RMSD 81 vs 100), and on the 8
   reactions where pprcht and xtb diverge most it is closer in **6/8** (only R30/R31 favour xtb).
   So curcuma tracking pprcht (its port source) follows the marginally better reference — there is
   **no QM justification to re-target xtb.**

Caveat: the session fixes made curcuma *faithful to pprcht* (MAD 0.47), the correct port goal; a
pre-fix curcuma happened to sit ~53 vs DLPNO because some of its bugs pointed toward QM by luck —
that was not a correctness advantage. The GFN-FF-vs-QM gap is the method's, not the port's.

### Remaining pure-port residuals (curcuma vs pprcht, where pprcht==xtb so not the split)

After the PR27 fcn fix the worst per-structure residuals are PR25 +2.87, ED25 +2.77, PR21 +2.51,
ED21 +2.27, PR16 +1.95, PR23 +1.63, PR35 +1.42 — all in the **fine-precision** class and all
verified (high-precision analyzer print, `6f14.8`) to have **no single-parameter smoking gun**:
the mismatches are sub-0.0016 in the bond force constant (from 3rd-decimal FT-HMO π-bond-order
differences on conjugated C=C/C=N bonds) and sub-0.01 Å in the metal-bond equilibrium r0 (the
CN-dependent `gfnffrab` quantity, scattered across element pairs). These correlated reactant/product
residuals **cancel at the reaction level** (reaction MAD 0.255). Irreducible fine-precision, not a
bug. Ground-truth with the `external/gfnff/_build/gfnff` analyzer (per-bond `pibo`/`fqq`/force-constant
print via `pr=.true.`; temporarily raise the `6f8.3` bond format for sub-mEh comparison).

The sections below are the original (pre-fix) diagnostic and remain valid for GFN-FF.

## What this is

MOR41 (Dohm, Hansen, Steinmetz, Grimme, Checinski, *J. Chem. Theory Comput.* 2018,
14, 2596) is a set of **41 realistic closed-shell metal-organic reactions**. Every
one of the 95 structures contains a transition metal (Co, Cr, Fe, Ir, Mn, Mo, Ni,
Pd, Pt, Rh, Ru, Ti, W). It is the transition-metal analog of the S30L check in
[S30L_GFNNF_VALIDATION.md](S30L_GFNNF_VALIDATION.md): the goal is to verify whether
curcuma's **native** GFN-FF / GFN1 / GFN2 reproduce the xtb 6.6.1 reference
implementation on this set. The DLPNO-CCSD(T)/CBS reference energies (SI Table S1)
are carried along only as accuracy context.

All structures are neutral closed-shell singlets: the SI title says "closed-shell",
SI Table S14 shows large singlet-triplet gaps for every molecule, and the electron
count of all 95 `mol.xyz` is even. So every single point is run at **charge 0,
multiplicity 1**, no `.CHRG`/`.UHF`.

## Method

`scripts/mor41_validation.py` runs a single point on each structure's cartesian
`mol.xyz` with both engines and the matching method
(`curcuma -sp -method {gfnff,gfn1,gfn2}` vs `xtb --sp {--gfnff,--gfn 1,--gfn 2}`),
caches every energy in `_run/energies.json`, then assembles the 41 reaction energies
`dE = sum(coeff*E)` from `test_cases/MOR41-testset/reactions.dat` (SI Table S1).
Two views are produced:

- **Per-structure** (`_run/mor41_per_structure.md`) — `E_curcuma - E_xtb` per
  molecule. This is the cleanest implementation check: it isolates exactly which
  molecules native curcuma fails to reproduce, and is immune to the fact that a
  reaction energy mixes several structures.
- **Per-reaction** (`_run/mor41_results*.md/.csv`) — `dE_cur`, `dE_xtb`, `dE_ref`
  and their differences.

## Outcome

### Per-structure agreement, curcuma − xtb (kcal/mol), n=95

| method | within 1 kcal/mol | MAD | max |
|--------|:-----------------:|----:|----:|
| gfnff  | **21 / 95** |   96 |   621 |
| gfn1   | **25 / 95** |  585 | 42736 |
| gfn2   | **26 / 95** | 1667 |  9825 |

### Per-reaction (kcal/mol), n=41

| method | curcuma vs xtb (MAD / max) | curcuma vs DLPNO (MAD) | xtb vs DLPNO (MAD / max) |
|--------|:--------------------------:|:----------------------:|:------------------------:|
| gfnff  |  57 / 282  |  53 |  62 / 279 |
| gfn1   | 1102 / 42981 | 54 | 1063 / 43080 |
| gfn2   |  507 / 5735 | 509 | **11.8 / 37** |

## Key findings

1. **Native curcuma reproduces xtb only for the metal-free molecules.** The
   structures that agree with xtb to < 1 kcal/mol are exactly the main-group
   ligands and gases (CO, CO2, H2, C2H4/C2H6/C3H8, Bz, COD, I2, MeCN, MeI/AcI/AcCl,
   MeOH, PhOH/PhSH/PhSeH, PMe3/PCy3, …) — e.g. CO gfn2 matches to 1e-8 Eh. This is
   consistent with the validated `sqm_reference` main-group set.

2. **Every transition-metal complex diverges, often massively.** Per-structure
   errors run to hundreds of kcal/mol and, for gfn2, up to **15.6 Eh** (ED33, a
   Ru/P complex: curcuma −77.14 Eh vs xtb −63.74 Eh). The SCF *converges* in
   curcuma (e.g. PR05/Mn: 23 Broyden iterations) — the energies are simply wrong,
   not unconverged. This sharpens the CLAUDE.md note "transition metals enabled but
   **unvalidated**" to: **native GFN1/GFN2/GFN-FF transition-metal energies do not
   match xtb and must not be used for TM systems.**

3. **gfn2 is the clean diagnosis.** xtb-gfn2 vs the DLPNO reference has MAD 11.8
   kcal/mol (max 37) — sane and in line with the paper — so xtb-gfn2 is a valid
   reference and the fault is entirely on the curcuma-native side (MAD 507 vs xtb).
   Heaviest failures are the 4d/5d complexes (Ru, Mo, Pd, Rh, Ir, Pt) and Se.

4. **xtb-gfn1 is itself numerically unstable on this set.** For PR40 (22 atoms)
   xtb `--gfn 1` reports "convergence criteria satisfied after 53 iterations" but at
   a spurious −93.6 Eh (ED40a/ED40b are −11.8/−13.1 Eh; curcuma-gfn1 gives the
   physical −25.5 Eh). This single xtb pathology dominates the gfn1 statistics
   (max 42736 kcal/mol). Hence the gfn1 comparison is doubly unreliable — both
   engines struggle with GFN1 on transition metals — and gfn1 curcuma-vs-xtb MAD
   should not be read as a curcuma error alone.

## What was NOT done / caveats

- No attempt was made to *fix* the native TM parametrization — this is a diagnostic
  run. The result is a negative one: native GFN methods are not TM-ready.
- xtb SCF settings left at defaults; no accelerated/level-shift options tried to
  rescue the xtb-gfn1 PR40 convergence (it reports success, so it would not be
  retried automatically anyway).
- DLPNO comparison is context only; GFN-family methods are not expected to reach
  the 2 kcal/mol reference accuracy.
- Geometries are the published MOR41 minima; no re-optimization.

## Reproduce

```bash
cd release && make -j4                      # ensure release/curcuma is current
python scripts/mor41_validation.py --only 5 # smoke test one small reaction
python scripts/mor41_validation.py          # full sweep, all three methods
# inspect:
#   test_cases/MOR41-testset/_run/mor41_per_structure.md   (per-molecule cur-xtb)
#   test_cases/MOR41-testset/_run/mor41_results.md         (per-reaction, all methods)
```

Energies are cached in `test_cases/MOR41-testset/_run/energies.json`; re-runs and
`--only <ids>` subsets are served from cache. Requires the xtb 6.6.1 binary at
`~/Downloads/xtb-6.6.1/bin/xtb`.
