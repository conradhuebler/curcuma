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

Still open in the GFN-FF metal path: `btyp=6` (eta) promotion is not wired, and the TM-TM
`mchar` attenuation is omitted. Neither is exercised by the largest remaining residuals.

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
