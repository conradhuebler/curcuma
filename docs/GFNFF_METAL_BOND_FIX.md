# GFN-FF Metal Bond-Stretch Fix — Implementation Report

Status: 🤖 AI-generated, ⚙️ machine-tested (2026-07-17). Human production testing pending.
Scope: native GFN-FF (`-method gfnff`) bond-stretch force constant for bonds touching a
metal. Metal-free / organic bonds are untouched (gated) and remain bit-identical.

Root-cause analysis this implements: [GFNFF_METAL_BOND_ANALYSIS.md](GFNFF_METAL_BOND_ANALYSIS.md).

## 1. What changed

Single file: `src/core/energy_calculators/ff_methods/gfnff_method.cpp`,
`GFNFF::getGFNFFBondParameters()`.

Added a **metal bond-stretch branch** (new block at ~line 4856, immediately after the
`mtyp1`/`mtyp2` classification) that ports the mutually-exclusive Fortran `else` branch
`gfnff_ini.f90:1188-1264` (bbtyp>=5). It is gated on `imetal1 > 0 || imetal2 > 0`, so
metal-free bonds skip it entirely and take the unchanged normal-bond path.

Inside the branch, for a metal bond it **overrides** the three factors the old code left on
the normal-bond path:

1. **`bstrength = bstren(bbtyp)`** (Fortran `gfnff_ini.f90:1190-1197`): M-X = `bstren(5)=1.00`,
   TM-TM = `bstren(7/8)=3.40` with the periodic-row correction. Replaces the hybridization
   `bsmat` value (which was 1.4842 for Ni-C, 1.3234 for Fe-C, etc.). This also fixes `alpha`,
   since `alpha` uses `srb3*bstrength` (self-corrects, no separate change needed).
2. **metal `fqq`** (Fortran `gfnff_ini.f90:1209-1211`): charge factor `25.0` (not the
   normal-bond `70.0`) and `qfacbm(mtyp)` averaged over the two atoms
   (`qfacbm = {1.0, -0.2, -0.2, 0.70, 0.50}` for mtyp 0..4). Minor for the carbonyls
   (~1.00) but correct in general.
3. **metal `fcn`** (Fortran `gfnff_ini.f90:1254-1259`): `1/(1+c·nb²)` with `c=0.100`
   (group 1/2), `0.030` (main-group metal), **`0.036` (transition metal)** — applied per
   metal atom. Uses the **bonded-neighbour count** `topo.neighbor_lists[i].size()`
   (== Fortran `topo%nb(20,i)`, the bonded degree), **not** `countNeighborsWithin20Bohr()`
   (a distance-sphere count that would be wrong). The branch resets `fcn = 1.0` first, so the
   erroneous normal-path `fcn` (which fires when both atoms have Z>10, e.g. Rh-Cl) is discarded.

The existing metal `fheavy` / `fpi` / `metal_shift` / `fsrb2` add-ons (lines ~4857-4966) were
already correct and are **not** duplicated. The bond r0/shift handling is unchanged (the
analysis confirmed native r0 already matches xtb for these systems).

### Not ported (residuals, no target exercises them)
- **eta bonds (bbtyp=6)**: need `itag(atom)==-1 && piadr(atom)>0` for the non-metal partner,
  not exposed at this call site; an eta bond falls back to M-X. No eta (π-coordinated) ligand
  in the validation targets.
- **TM-TM mchar "metallic" scaling** (`gfnff_ini.f90:1195-1197`, `1-min(2·mchar_i+2·mchar_j,0.5)`):
  `mchar` is not available in curcuma. The `bstren` selection for TM-TM is ported; the mchar
  attenuation is omitted. No M-M bond in the validation targets.

## 2. Results — metal complexes vs xtb 6.6.1 `--gfnff`

Command (native): `curcuma -sp <xyz> -method gfnff -charge 0 -no_bmt`.
Reference: `xtb 6.6.1 mol.xyz --sp --gfnff`.

### Bond energy (Eh)

| System            | native BEFORE | native AFTER | xtb 6.6.1   | Δ bond AFTER |
|-------------------|---------------|--------------|-------------|--------------|
| ED03 Ni(CO)3      | -1.5723589    | -1.3628209   | -1.3627544  | **-6.6e-5**  |
| ED02 Fe(CO)4      | -2.1682121    | -1.8300341   | -1.8299549  | **-7.9e-5**  |
| ED39 RhCl(CO)2    | -1.1565627    | -0.9836324   | -0.9965633  | +1.29e-2     |
| PR01 Cr(CO)6      | (n/a)         | -2.6395568   | -2.6404317  | +8.7e-4      |
| PR05 Mn(CO)5H     | (n/a)         | -2.6505557   | -2.6512630  | +7.1e-4      |

### Total energy (Eh)

| System            | native BEFORE | native AFTER | xtb 6.6.1   | Δ total AFTER |
|-------------------|---------------|--------------|-------------|---------------|
| ED03 Ni(CO)3      | -1.6206979    | -1.4111599   | -1.4114817  | +3.2e-4       |
| ED02 Fe(CO)4      | -2.2210729    | -1.8828948   | -1.8761418  | -6.75e-3      |
| ED39 RhCl(CO)2    | -1.2133981    | -1.0404678   | -1.0533869  | +1.29e-2      |
| PR01 Cr(CO)6      | (n/a)         | -2.7083509   | -2.7117825  | +3.4e-3       |
| PR05 Mn(CO)5H     | (n/a)         | -2.7257818   | -2.7231681  | -2.6e-3       |

**The bond-stretch over-binding is fixed.** ED03 and ED02 (Ni/Fe carbonyls, the primary
targets) match xtb bond energy to <1e-4 Eh; the pure carbonyls PR01/PR05 to <1e-3 Eh. The
per-complex bond error dropped from -131 / -212 / -100 kcal/mol to <0.05 / <0.05 / +8.1 kcal/mol.

### Residuals (honest)
- **ED39 RhCl(CO)2 (+12.9 mEh bond, +12.9 mEh total)**: the Rh-C(O) and C-O ligand bonds match
  the ED03 pattern (which is exact); the residual localizes to the **Rh-Cl halide bond**
  (`fheavy=1.30` TM-halogen). The dominant fix reduced this from -160 mEh to +12.9 mEh
  (native now slightly *under*-binds). Likely a small r0/shift or halide-specific difference,
  not the force-constant factors ported here.
- **ED02 Fe(CO)4 total (-6.75 mEh)** despite the bond matching to 7.9e-5: this residual is in a
  **non-bond term** (Fe angle +0.0045 Eh, plus small Coulomb/repulsion deltas), pre-existing and
  outside the scope of this bond-stretch fix.

## 3. No regression on metal-free systems

The branch is gated on `imetal1>0 || imetal2>0`; metal-free bonds are structurally untouched.

- `cd release && ctest -R gfnff` → **75/76 pass**. The neutral `gfnff` validation label
  (21 tests: caffeine, triose, acetic_acid_dimer, …) and the `validation` label (37 tests) all
  pass against their fixed golden energies — i.e. metal-free energies are unchanged.
- The single failure, `gfnff_gpu_vs_cpu_h2`, aborts with `File not found:
  ../../external/gfnff/test/h2.xyz` (a missing test-fixture file, before any calculation).
  This is a **pre-existing infrastructure failure unrelated to this change** (H2 is metal-free;
  the file is genuinely absent from the tree).
- Spot checks unchanged: benzene (`Bz`) total -2.3633652, C2H6 -1.0718002, CO2 -0.6874420.

## 4. Files
- `src/core/energy_calculators/ff_methods/gfnff_method.cpp` — metal bond-stretch branch in
  `getGFNFFBondParameters()` (~line 4856). No other source file modified.
