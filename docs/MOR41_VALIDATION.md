# MOR41: Native GFN-FF / GFN1 / GFN2 vs xtb 6.6.1 — Validation

Status: AI-generated, machine-tested (Jul 2026). NOT human production-tested.
The AI does not assign ✅ TESTED / ✅ APPROVED.

## Update (Jul 16, 2026): GFN2 transition-metal shell-indexing bug FIXED

The original negative result below was traced to a concrete bug and **GFN2 is now
fixed for 3d transition metals**. Three GFN2 parameter tables — `reference_occ`,
`p_kcn`, `p_shpoly` — are indexed by **angular momentum** in tblite
(`gfn2.f90`: `(0:2, elem)`), but `xtb_native.cpp` read them by **shell index**.
Main-group is unaffected (shells stored `[s,p,d]`, so `ang_sh==ish`); transition
metals order shells `[d,s,p]`, so shell-indexing scrambled the d/s/p reference
occupations and CN/shpoly → wrong `n0`/`dq` → multi-Eh errors. Fix: index those three
by the shell's angular momentum for GFN2 (`xtb_native.cpp` buildReferenceOccupations +
H0 setup).

Result (native GFN2 vs tblite-GFN2, identical geometries):
- **3d metals (Cr, Mn, Fe, Co, Ni): EXACT to tblite (1e-8)** — e.g. Fe(CO)4 was 4.93 Eh
  off, now 0.0.
- 4d (Mo, Ru, Rh, Pd): ~1e-3 Eh residual; 5d (W, Ir, Pt): ~0.1–0.5 Eh residual — a
  separate, unresolved heavy-element band-energy issue (Rh's H0 params match tblite
  exactly, so it is not another indexing bug).
- MOR41 reaction energies, native GFN2 vs xtb 6.6.1: **MAD 507 → 4.29 kcal/mol**
  (max 5735 → 35; the residual is entirely the 4d/5d reactions). curcuma-GFN2 vs DLPNO
  MAD is now 15.4, comparable to xtb-GFN2's own 11.8.
- **No regression**: all main-group sqm molecules (incl. d-shell H2S/PH3/SiH4/HCl) stay
  1e-8 vs tblite for both GFN1 and GFN2; gfn2 ecomp ctests pass.

**GFN1 is NOT fixed** and still shows the TM error. GFN1's basis has valence+polarisation
shells that share an angular momentum (e.g. H = two s-shells with distinct `p_kcn`, and
tblite `gfn1.f90` keeps `p_kcn` shell-indexed), so the same angular-momentum remap would
collapse those shells (verified: it breaks even H2O/SiH4). A GFN1 TM fix needs
valence/polarisation-aware shell mapping and is deferred. GFN-FF TM deviations are a
separate force-field-parameter matter (also open).

The sections below are the original (pre-fix) diagnostic and remain valid for GFN1 and
GFN-FF.

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
