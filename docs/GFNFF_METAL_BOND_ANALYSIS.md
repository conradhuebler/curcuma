# GFN-FF Metal Bond-Stretch Overbinding — Root-Cause Analysis

Status: 🤖 AI-generated diagnostic (2026-07-16). No source code modified.
Scope: native GFN-FF (`-method gfnff`) bond-stretch term for bonds involving a
transition metal. Organic/metal-free bonds are unaffected.

## 1. Confirmed measurement (ED03 = Ni(CO)3, 7 atoms, neutral)

Command:
`./release/curcuma -sp test_cases/MOR41-testset/ED03/mol.xyz -method gfnff -charge 0 -verbosity 3 -no_bmt`
vs `xtb 6.6.1 mol.xyz --sp --gfnff`.

| Component      | native (Eh) | xtb --gfnff (Eh) | delta (Eh) |
|----------------|-------------|------------------|------------|
| **bond**       | **-1.5723589** | **-1.3627544** | **-0.209604** |
| angle          | +3.9e-8     | +4.4e-8          | ~0 |
| repulsion      | +0.106413   | +0.106336        | +8e-5 |
| dispersion     | -0.0034168  | -0.0034321       | +1.5e-5 |
| electrostatic  | -0.1517677  | -0.1521842       | +4e-4 |
| bonded ATM     | +0.000553   | +0.000553        | ~0 |

Every term except **bond** matches to ~1e-4 Eh. The entire discrepancy
(-0.2096 Eh = -131.5 kcal/mol) is in the bond-stretch term.

## 2. Per-bond breakdown (from `BOND_CSV` / `BOND_FACTORS`, verbosity 3)

`E_bond = k_b * exp(-alpha*(r-r0)^2)`, `k_b` stored negative.

| Bond (3x each) | rij (Bohr) | r0 | k_b (fc) | alpha | E per bond | 3-bond sum |
|----------------|-----------|-----|----------|-------|-----------|-----------|
| Ni-C (Z28-Z6)  | 3.4084 | 2.9518 | -0.159531 | 0.502982 | -0.143646 | **-0.430938** |
| C-O  (Z6-Z8)   | 2.1704 | 1.8289 | -0.409752 | 0.635556 | -0.380471 | -1.141413 |
|                |        |        |          |        | **total** | **-1.572351** |

The C-O ligand bonds are ordinary covalent bonds and are correct. Subtracting the
native C-O sum from the xtb bond total isolates the xtb Ni-C contribution:

`xtb Ni-C = -1.3627544 - (-1.141413) = -0.221341 Eh`  vs  native Ni-C `-0.430938 Eh`.

**Native Ni-C bonds are 1.95x too strong. The full -0.2096 Eh error is the Ni-C term.**

## 3. Root cause: the Fortran metal-bond branch is not implemented

### Fortran reference (`external/xtb/src/gfnff/gfnff_ini.f90`)

Bond parameter setup splits on bond type at line 1067. Every bond touching a metal
gets `btyp=5` (metal), 6 (eta), or 7 (M-M) at lines **1061-1064**, and then takes a
**mutually exclusive `else` branch (lines 1129-1206)** that recomputes the force
constant from scratch — it does NOT reuse the normal-bond `bstrength`, `fqq`, or `fcn`:

- **line 1131**: `bstrength = gen%bstren(bbtyp)`  → for M-X (`bbtyp=5`) this is
  `gen%bstren(5) = 1.00` (`gfnff_param.f90:654`). The hybridization/`bsmat` value is
  discarded.
- **lines 1150-1152**: metal `fqq` uses charge factor **25.0** (not 70.0) and
  `gen%qfacbm(mtyp)` (not `qfacbm0=0.047`).
- **lines 1195-1200**: metal `fcn` is applied per metal atom regardless of the ligand
  being light:
  `if(mtyp==4) fcn = fcn/(1+0.036*nb(20,i)^2)` (TM),
  `1+0.100*...` (group 1/2), `1+0.030*...` (main-group metal).
  Here `topo%nb(20,i)` is the **count of bonded neighbours**, not a distance count.
- lines 1167-1186 (M-CO/M-CN `fpi=1.5`/`shift=-0.45`), 1187-1200 (metal shifts),
  1201-1205 (`fsrb2` TM sign flip) — these curcuma already has.
- final assembly, line **1226**:
  `vbond(3) = -bond(ia)*bond(ja)*ringf*bstrength*fqq*fheavy*fpi*fxh*fcn`.

`gen%bstren` (`gfnff_param.f90:650-656`): 1=1.00, 2=1.24, 3=1.98, 4=1.22,
**5(M-X)=1.00**, 6(M-eta)=0.78, 7/8(M-M)=3.40.

### What curcuma does (`src/core/energy_calculators/ff_methods/gfnff_method.cpp`)

`getGFNFFBondParameters()` (line 4184) has **no metal branch**. It always follows the
normal-bond path and applies the metal corrections as multiplicative *add-ons*:

- **bstrength** computed only from hybridization/`bsmat` at lines **4464-4526**; there
  is no `if (metal bond) bstrength = bstren(5)` override. For Ni-C this yields
  **1.4842** (should be **1.00**).
- **fqq** always uses the normal formula (factor 70.0, `qfacbm0=0.047`), lines
  **4576-4579**. (Numerically ~1.0 for ED03, so minor here, but wrong in general.)
- **fcn** applied only `if (z1>10 && z2>10)` with factor 0.007, lines **4697-4709**.
  For Ni-C the ligand C is light (Z=6), so the guard fails and **fcn=1.0**
  (should be `1/(1+0.036*nb^2)`).
- metal `fheavy`/`fpi`/`metal_shift` add-ons at lines 4834-4924 (these are correct).
- fc assembled at line **4935**:
  `-(bond_i*bond_j*bstrength*fqq*ringf*fheavy*fpi*fxh*fcn)`.
- `alpha` (lines 4947-4978) contains `srb3*bstrength`, so the inflated bstrength also
  makes alpha too steep (0.5030 vs Fortran 0.4571) — a secondary consequence.

### Numerical reconciliation for Ni-C

Fortran-correct factors: `bstrength=1.00`, `fpi=1.5`, `fheavy=1.0`, `fqq≈1.0`,
`fcn = 1/(1+0.036*nb^2)` with `nb(20,Ni)=3` bonded C → `fcn=1/1.324=0.7553`.

- fc_correct = -(0.186003 * 0.385248 * 1.00 * 1.5 * 0.7553) = **-0.081196**
  (vs native -0.159531 → factor 0.509).
- alpha_correct = 0.3731*(1 - 0.06976*dEN^2 + 0.2538*1.00) = **0.4571**
- E_correct = -0.081196 * exp(-0.4571*0.4566^2) = -0.073815/bond, x3 = **-0.221446 Eh**

This matches the xtb-implied Ni-C sum (-0.221341 Eh) to ~1e-4 Eh, i.e. the two missing
pieces (bstrength 1.4842→1.00 and fcn 1.0→0.7553) fully account for the -0.2096 Eh error.
The two factors contribute roughly equally (0.674 and 0.755 → combined ~0.509, "~half").

## 4. Metal-general (Fe, Rh) — same mechanism

| System | native bond (Eh) | xtb bond (Eh) | delta (Eh / kcal) | metal factors (native) |
|--------|------------------|---------------|-------------------|------------------------|
| ED03 Ni(CO)3 | -1.572359 | -1.362754 | -0.2096 / -131.5 | bstr 1.4842, fcn 1.0, nb(Ni)=3 |
| ED02 Fe(CO)4 | -2.168212 | -1.829955 | -0.3383 / -212.3 | bstr 1.3234, fcn 1.0, nb(Fe)=4 |
| ED39 Rh...   | -1.156563 | -0.996563 | -0.1600 / -100.4 | (metal M-C/M-Cl bonds) |

For Fe-C the correct factors give `(1.00/1.3234)*[1/(1+0.036*4^2)] = 0.756*0.6345 = 0.480`
per bond (≈2.08x too strong), and Fe has 4 CO ligands, so its absolute error is larger.
The overbinding scales with metal coordination number, as expected from the missing
`fcn` and `bstren(5)` terms.

## 5. Localized causes (file:line)

1. **Missing metal-bond branch / wrong `bstrength`** — dominant.
   `gfnff_method.cpp:4464-4526` (bstrength from `bsmat` only; no metal override).
   Must mirror Fortran `gfnff_ini.f90:1129-1131` — for a metal bond set
   `bstrength = bstren(bbtyp)` (M-X=1.00, eta=0.78, M-M=3.40 with the row/mchar
   corrections at 1132-1139).

2. **Missing metal `fcn`** — co-dominant.
   `gfnff_method.cpp:4697-4709` (only fires when both atoms Z>10, factor 0.007).
   Must add the per-metal-atom factors of Fortran `gfnff_ini.f90:1195-1200`
   (0.100 group1/2, 0.030 main-group, **0.036 TM**).
   NOTE: this requires the **bonded-neighbour count** `nb(20,i)`, but curcuma's
   `countNeighborsWithin20Bohr()` (`gfnff_method.cpp:7320-7349`) counts every atom
   within a 20-Bohr sphere, which for a small isolated complex is ~(N-1), not the
   bonded degree. Fixing the metal `fcn` must use the actual bonded-neighbour count
   (e.g. `topo.neighbor_lists[i].size()`), not this function.

3. **Wrong metal `fqq`** — minor for these systems, wrong in general.
   `gfnff_method.cpp:4576-4579` uses the normal formula. Metal bonds need
   `gfnff_ini.f90:1150-1152` (factor 25.0, `qfacbm(mtyp)` averaged).

4. **`alpha` steepness inflated** — secondary consequence of (1).
   `gfnff_method.cpp:4976` (`srb3*bstrength`); self-corrects once bstrength is fixed.

Bond energy is evaluated in `forcefieldthread.cpp:921` (`CalculateGFNFFBondContribution`,
formula lines 977-998) — that code is correct; the error is entirely in the parameters
fed to it.

## 6. Recommendation

Add the Fortran metal branch to `getGFNFFBondParameters()`: when a bond has
`imetal1>0 || imetal2>0` (bbtyp>=5), stop using the `bsmat` hybridization value and
instead:
- set `bstrength = bstren(bbtyp)` (M-X=1.00, M-eta=0.78, M-M=3.40 + the TM-TM
  `mchar`/row corrections, `gfnff_ini.f90:1131-1139`),
- compute `fqq` with the metal formula (25.0 / `qfacbm(mtyp)`, lines 1150-1152),
- apply the metal `fcn` factors (0.100/0.030/0.036) using the **bonded-neighbour
  count**, lines 1195-1200.
Keep the existing `fheavy`/`fpi`/`metal_shift`/`fsrb2` metal add-ons (already correct).

Validate against the MOR41 metal-carbonyl subset (ED02/ED03/ED39) targeting the xtb
`:: bond energy` line to ~1e-4 Eh, and re-run the S30L / neutral validation set to
confirm zero regression on metal-free systems (the branch is gated on `imetal>0`, so
organics are untouched).
