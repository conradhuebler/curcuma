# QMDFF in Curcuma

**Status: 🤖 AI-generated, ⚙️ machine-tested — human production testing pending.**

QMDFF (Quantum Mechanically Derived Force Field) was introduced by S. Grimme,
*J. Chem. Theory Comput.* **2014**, *10*, 4497–4514 (DOI
[10.1021/ct500573f](https://doi.org/10.1021/ct500573f)).

Curcuma's `-method qmdff` implements the **QMDFF energy expression** as a port of
the Fortran evaluator that used to ship with xtb as `src/qmdff.f90` (module
`xtb_qmdff`). That file was deleted from xtb in commit `b7dbd36` ("Refactoring of
external drivers", PR #568); the last version is recoverable from the vendored
repository at `external/xtb`:

```bash
git -C external/xtb show b7dbd36^:src/qmdff.f90
```

All line references in the source comments (`qmdff.f90:NNN`) point into that file.

---

## What is implemented

| Term | Formula | Reference |
|------|---------|-----------|
| Bond (LJ) | `E = k [1 + (r₀/r)^a − 2 (r₀/r)^(a/2)]` | `qmdff.f90:329-334` |
| Bond (Morse) | `E = k [1 − exp(−α(r−r₀))]²`, `α = a/(2r₀)` | `qmdff.f90:335-338` |
| Angle | `E = k (cosθ − cosθ₀)² · damp`, harmonic in θ for linear θ₀ | `qmdff.f90:357-402` |
| Torsion | Fourier sum, erf-switched between the ±φ₀ branches, damped | `qmdff.f90:405-464` |
| Out-of-plane | `E = v(1 − cos(ω−ω₀))` or `v(cosω − cosω₀)²`, damped | `qmdff.f90:466-504` |
| Dispersion | **D4** with Becke–Johnson damping, `s8 = 2.0`, `a1 = 0.58`, `a2 = 4.80` (see below) | `qmdff.f90:541-553` |
| Electrostatics | `E = ε₁(nk) qᵢqⱼ / r` | `qmdff.f90:561-569` |
| Repulsion | `E = ε₂(nk) ZᵢZⱼ exp(−α r)/r`, `α = 16.5/R0ab^1.5` | `qmdff.f90:571-584` |
| Hydrogen bond | `E = −DA(q) · damp_long · angle_term / r_AB³`, analytic gradient | `qmdff.f90:766-881` |
| Halogen bond | `E = −c_A · damp_long · angle_term / r_XB²`, numerical gradient | `qmdff.f90:887-986` |

Supporting pieces ported verbatim: the bend/torsion dissociation damping
`abdamp` (`qmdff.f90:633-648`), the charge-dependent HB/XB scaling `hbpara`
(`qmdff.f90:753-760`), the effective valence-electron table `valelff`
(`qmdff.f90:1084-1158`), the LJ→Morse promotion `setmorse` (`qmdff.f90:994-1019`),
and the D3 `R0ab` table from `xtb/src/disp/dftd3.f`.

### Source layout

| File | Contents |
|------|----------|
| `src/core/energy_calculators/ff_methods/qmdff_terms.h` | All energy/gradient kernels + the term structs. Self-contained, atomic units. |
| `src/core/energy_calculators/ff_methods/qmdff_data.{h,cpp}` | Element tables: D3 `R0ab` (4465 packed values), atomic radii, valence electrons, `scalehb`/`scalexb`, `eps1`/`eps2`. |
| `src/core/energy_calculators/ff_methods/ff_workspace_uff.cpp` | `executeQMDFF` + the five `calcQMDFF*` drivers (threaded, partitioned). |
| `src/core/energy_calculators/ff_methods/forcefieldgenerator.cpp` | `setQMDFFTerms()` — builds the torsion / NCI / HB lists. |
| `test_cases/test_qmdff_gradients.cpp` | Analytic vs. finite-difference gradient test (ctest `qmdff_gradients`). |

### Units

The kernels in `qmdff_terms.h` work in **atomic units** (Bohr, Hartree), like the
Fortran reference. `FFWorkspace` stores the UFF/QMDFF geometry in **Ångström**, so
each `calcQMDFF*` driver scales positions by `m_au = 1.889726125` on the way in and
the gradient by the same factor on the way out. The generic `Bond`/`Angle` structs
keep curcuma's Ångström convention (`Bond::r0_ij`); the QMDFF-specific structs
(`QMDFFNonCovalent::r0_bj`, `alpha`) are in atomic units.

---

## What is *not* the reference: the parameters

This is the most important caveat.

The published QMDFF is a **fit to a QM Hessian** carried out by Grimme's separate
QMDFF program, which writes a parameter file that `qmdff.f90` then merely
evaluates.

Curcuma now has that fit (`-qmdfffit`, see below), but it is **not** what a bare
`-method qmdff` gives you. Without a fit, curcuma derives the parameters from its own
UFF-style topology perception — that is what this table describes:

| List | Origin in curcuma |
|------|-------------------|
| Bond `r₀` | the input structure's bond distance |
| Bond `k` | fixed default (0.01) — meant to be refined by `-qmdfffit` |
| Bond `a` | `ka(Zi)·ka(Zj) + kEN·ΔEN²`, scaled (`qmdff_par.h`) |
| Angle `θ₀` | the input structure's angle |
| Angle `k` | UFF force constant, refined by `-qmdfffit` |
| Torsion | curcuma's UFF torsion terms, re-expressed exactly in the QMDFF Fourier form (see below) |
| Out-of-plane | UFF improper, re-expressed in the QMDFF single-minimum form; **the barrier magnitude is UFF's, not a fitted QMDFF value** |
| NCI (dispersion/repulsion) | essentially parameter-free: **C6 from D4** (charge- and CN-dependent, using the same QM charges as the electrostatics), R0 from `r2r4`, Z from `valelff`, α from D3 `R0ab` |
| NCI (electrostatics) | needs QM partial charges — see below |
| HB/XB | topological detection + the reference's `hbpara` charge scaling |

**So without `-qmdfffit`: the functional form is the reference one, the parameters are
curcuma's.** Do not compare absolute energies against published QMDFF numbers. After a
fit, the bond/angle force constants and torsion/out-of-plane barriers come from a QM
Hessian — the topology, r₀/θ₀ and the bond exponent still do not.

Because `r₀` and `θ₀` are taken from the input structure, a single point at that
structure has *zero* bond and angle energy by construction — QMDFF is a local
force field around a reference geometry.

### Torsion mapping (exact)

Curcuma's UFF torsion `E = V/2 · (1 − cos(nφ₀)cos(nφ))` is written as two QMDFF
Fourier terms with the QMDFF reference angle set to π. At φ₀ = π both erf branches
coincide (`dphi1 == dphi2 == φ − π`), so

```
v₁(1 + cos nφ) + v₂(1 − cos nφ) = (v₁+v₂) + (v₁−v₂) cos nφ
v₁ = V/4 (1 − cos nφ₀),   v₂ = V/4 (1 + cos nφ₀)
```

reproduces the UFF term exactly for any φ₀, while gaining QMDFF's dissociation
damping. The phases `nπ` and `nπ+π` undo the `φ → φ−π` shift.

### Partial charges

Electrostatics, hydrogen bonds and halogen bonds all need QM-derived atomic
charges. Curcuma has none for a bare single point, so those three terms are
**switched off** and a warning is printed at verbosity ≥ 1:

```
[WARN] QMDFF: no partial charges available — electrostatics and HB/XB terms are
       switched off. Use -qmdfffit to derive charges from a QM calculation.
```

`-qmdfffit` computes GFN2 charges, stores them in `scf.json`, and writes a
parameter set that carries `qmdff_ncis` with populated `qq` values.

---

## Usage

```bash
# Single point / optimisation with generated parameters
curcuma -sp  molecule.xyz -method qmdff
curcuma -opt molecule.xyz -method qmdff

# Derive charges + fit force constants against a GFN2 Hessian
curcuma -qmdfffit molecule.xyz -threads 4
```

Relevant parameters (registry module `forcefield`):

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `qmdff_morsethr` | 99.0 | Well depth above which a 1,2 stretch switches from Lennard-Jones to Morse (`setmorse`) |

---

## Validation

`ctest -R qmdff_gradients` compares the analytic gradient against a central
finite difference of the same parameter set, on a displaced geometry (parameters
generated at the reference structure, energy evaluated away from it — the
situation an optimiser or MD run produces).

| Molecule | Terms exercised | max &#124;Δg&#124; (Eh/Å) |
|----------|-----------------|-----------------|
| H2O | bond, angle | 4.0e-11 |
| CH4 | bond, angle | 8.7e-12 |
| CH3OH | + torsion, electrostatics | 1.4e-11 |
| CH3OCH3 | + several torsions | 1.3e-11 |
| C6H6 | + ring torsions, out-of-plane | 6.6e-11 |
| H2O_dimer | + NCI, hydrogen bond | 4.2e-11 |
| CH3Cl···H2O (manual) | + halogen bond | 4.8e-11 |
| CH3Br···H2O (manual) | + halogen bond | 4.9e-11 |

**What this establishes**: the analytic derivatives are consistent with the energy
expressions, for every term group.

**What this does not establish**: agreement with the Fortran reference or with
published QMDFF results. No QMDFF reference implementation is available to compare
against — `qmdff.f90` was removed from xtb, and it required parameter files that
curcuma does not produce. The port was done by line-by-line transcription and is
**not numerically validated against the original**.

### Deliberately kept from the reference

The halogen-bond gradient is a **central finite difference with step 1e-6 Bohr**,
exactly as in `eabxag` (`qmdff.f90:933-986`). An independently derived analytic
gradient would be a silent source of disagreement with the reference; the XB term
contributes little to the total force, so the cost is irrelevant.

The dispersion coefficients C6 are computed **once**, from the coordination
numbers of the structure the parameters are generated on (`ff_ini`,
`qmdff.f90:146-155`), and treated as constants by `ff_nonb`. Neither the reference
nor this port carries a dC6/dCN response term in the gradient.

---

## Bugs found and fixed on the way

Two defects in the pre-port QMDFF code, both of which produced plausible but wrong
energies:

1. **Bond exponent** — `calcQMDFFBonds` used `2·(r₀/r)^(0.75a)` where the QMDFF
   form is `2·(r₀/r)^(a/2)`. The attractive branch was therefore too steep.
2. **Angle reference unit** — `theta0_ijk` is stored in radians (the generator
   writes `acos(costheta)`, and GFN-FF reads it as radians), but both QMDFF angle
   routines converted it as if it were degrees (`cos(θ₀·π/180)`). Every reference
   angle was wrong.

Consequence: `-method qmdff` energies change. For CH3OCH3 the single point goes
from 2.58 Eh (nonsense — larger than a covalent bond) to 0.032 Eh, next to UFF's
0.048 Eh on the same structure.

---

## Parametrisation: `-qmdfffit` (August 2026)

QMDFF's defining feature is that its force constants come from a **QM Hessian**.
Curcuma now does that, in `src/core/energy_calculators/ff_methods/qmdff_parametrisation.{h,cpp}`.

The key observation is that the QMDFF energy is **exactly linear** in every constant
worth determining — bond `fc`, angle `fc`, torsion `scale`, out-of-plane `scale`:

```
H_FF(x0) = sum_p  k_p * U_p(x0)
```

where `U_p` is the Hessian of term `p` with its own constant set to 1. Fitting is
therefore a single **linear least-squares solve**, and no force-field Hessian is
needed at all. Two stages:

1. **Starting guess** — Seminario projection for bonds and angles
   (Seminario, *IJQC* **60** (1996) 1271; symmetrised per Allen/Payne/Cole,
   *JCTC* **14** (2018) 274). For torsions and out-of-plane terms there is *no*
   defensible Seminario analogue — a torsion constant is a barrier height, which a
   Hessian at one geometry cannot see — so the Rayleigh quotient
   `<U_p,H>/<U_p,U_p>` is used, which is the exact single-parameter solution and is
   well defined for every term kind.
2. **Global least squares** — normal equations `(A + λ diag A) k = b + λ diag A k0`
   with `A_pq = <U_p,U_q>`, assembled only over terms sharing an atom (O(N)). λ is
   chosen by a ladder sweep; non-negativity via a bounded Lawson–Hanson active set.

The `U_p` are built by central differences of the **already FD-validated analytic
term gradients** in `qmdff_terms.h`, over only the term's own 3·k coordinates:
12–24 tiny kernel calls per term, O(1) each.

### Cost

| | 9 atoms | 50 atoms |
|---|---|---|
| old Levenberg–Marquardt fit | 37 260 FF-gradient calls per step | 1.33·10⁶ per step |
| new linear fit | **0** | **0** |
| remaining cost | the QM Hessian (6N gradient calls) | same |

### Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `method` | gfn2 | QM method for the reference Hessian and charges |
| `lambda` | −1 (auto) | Tikhonov strength for bonds/angles |
| `lambda_torsion` | 1.0 | separate, much stronger prior on torsion/out-of-plane barriers |
| `nonnegative` | true | constrain constants to physical lower bounds |
| `tie_torsions` | true | one shared scale per central bond |
| `fd_step` | 1e-4 Bohr | step for the unit Hessians |
| `verify` | true | recompute one FF Hessian and compare |
| `write_param` | true | also write `<input>.param.json` |

### Two non-obvious findings, both fixed

- **Do not project translations/rotations out of the target.** A term's Hessian
  annihilates rotations only where that term's *gradient* vanishes. Bonds and angles
  sit at their reference by construction, but a QMDFF torsion generally does not, so
  its unit Hessian legitimately carries a rotational component. Projecting the target
  while the model keeps it makes the least-squares problem inconsistent — on CH3OCH3
  it turned an exactly recoverable synthetic system into a 3.5e-3 residual with
  constants off by 22%. `project_tr` is therefore off by default; the information is
  reported as `tr_content` instead.
- **Torsion barriers need their own prior.** Fitted on the same footing as bonds and
  angles, every torsion scale is driven to zero — the Hessian carries almost no
  barrier information and the overlapping angle terms absorb what little there is.
  A force field with free internal rotation is not a useful outcome, hence
  `lambda_torsion = 1`.

### Validation of the fit

`ctest -R qmdff_hessian_fit` — a **synthetic round trip**: known force constants →
Hessian built through the *real* FFWorkspace path → fit back → require recovery.
Independent of the fit's own machinery, so an indexing slip, an Å-vs-Bohr factor or a
mass-weighted target cannot pass.

| check | CH3OCH3 | C6H6 | H2O_dimer |
|---|---|---|---|
| model with `k_true` reproduces the synthetic Hessian | 5.6e-8 | 1.3e-7 | 6.8e-8 |
| fitted Hessian reproduces it | 2.1e-8 | 3.0e-8 | 4.2e-8 |
| force constants recovered (max rel. error) | 1.5e-6 | 1.6e-6 | 7.8e-7 |

Plus closed-form checks (`U_bond == a²/(2r₀²)·bbᵀ`, LJ ≡ Morse, to 4e-9) and
translational invariance of every unit Hessian.

Each real run additionally verifies itself: one force-field Hessian on the fitted
set is compared against the fit's own prediction (`4.0e-08` on CH3OCH3), which proves
in one number that the parameters reached `ForceField`, that the unit-Hessian
indexing and units are right, and that the linearity assumption holds.

### Multi-basin fitting (`-basins N`)

Every residual type is linear in the same constants, so several fitting points stack into
one system — `A = Sum_b A(x_b)`, `b = Sum_b b(x_b)` — giving ONE parameter set for the
whole ensemble, with no basin-assignment problem and no discontinuity between basins:

```
Hessian  at basin b :  Sum_p k_p U_p(x_b)           ~ H_QM(x_b) - H_nb(x_b)
energy   at basin b :  Sum_p k_p [e_p(b) - e_p(0)]  ~ dE_QM(b,0) - dE_nb(b,0)
gradient at basin b :  Sum_p k_p g_p(x_b)           ~ g_QM(x_b) - g_nb(x_b)
```

Each block is normalised by its own target norm, so `weight_hessian` / `weight_energy` /
`weight_gradient` are dimensionless.

Measured on the 107-atom peptide, 4 basins:

| configuration | rel. residual | R² | constants at floor |
|---|---|---|---|
| 1 basin, Hessian only | 0.1344 | +0.982 | 17 / 364 |
| 4 basins, Hessian only | **0.1339** | **+0.982** | 17 / 364 |
| 4 basins, Hessian + energy | 0.1832 | +0.966 | 105 / 364 |
| 4 basins, all three blocks | 1.0147 | **−0.030** | 120 / 364 |

Two things follow.

**Multi-basin Hessians are free.** Four basins fit as well as one (0.1339 vs 0.1344), so
the force constants really are transferable between basins — the curvatures do not
conflict. But by the same token this cannot fix conformer ranking: it barely changes the
solution.

**The gradient block must stay off** (`weight_gradient` default 0) unless every basin
shares the reference internal coordinates. `r0`/`theta0` belong to the structure the
parameter set was generated on; at any other basin the bonded terms are displaced, so the
force field has a large gradient there while the reference gradient at an optimised
conformer is about 0. With non-negative constants the only way to satisfy that is to
switch terms off — visible as 120 of 364 constants pinned at their floor and R² negative.

The energy block sits in between: it costs little in fit quality but already pins 105
constants, i.e. it too is fighting a reference geometry that belongs to a different basin.
**The remaining obstacle is not the force constants — it is that `r0`/`theta0`/`phi0` are
basin-specific and enter non-linearly.**

### Fitting GFN-FF's terms instead (`-potential gfnff`)

Every GFN-FF bonded term is linear in exactly one constant, so the same engine applies:

| term | form | constant |
|---|---|---|
| bond | `fc·exp(−α·dr²)` | `Bond::fc` (**negative** — it is a well depth) |
| angle | `fc·(cosθ−cosθ₀)²·damp` | `Angle::fc` |
| torsion | `V·(1+cos(n(φ−φ₀)+π))·damp` | `Dihedral::V` |
| inversion | `V·(1−cosω)·damp` | `Inversion::fc` |

Motivation: keep GFN-FF's non-covalent block — which carries ~80 % of the conformational
spread and beat QMDFF's on the peptide — and fit only the bonded constants to the QM
Hessian. Kernels in `gfnff_unit_terms.h`, selected by `FitOptions::source`.

**Status: verified end to end.** The frozen-parameter model reproduces the production
GFN-FF Hessian at the fitted constants to **5.9e-09** (CH3OCH3, 9 atoms) and **6.4e-09**
(WEKLQ peptide, 107 atoms, 601 constants). ⚙️ Machine-tested; not yet human production tested.

Five real defects were found on the way there, each invisible in the energy alone:

1. **Sign convention.** GFN-FF's bond constant is negative. The shared non-negativity
   constraint pinned all 8 CH3OCH3 bonds at `+1e-4`. The solver now works internally with
   `s = sign·k >= 0` and undoes the flip on write-back.
2. **Torsion tying.** Tying torsions about a shared central bond is only valid when they
   share a constant. QMDFF's `scale` is 1 for all of them; GFN-FF gives every dihedral its
   own `V`, so collapsing a group to one number corrupted both the guess and the remainder
   (1.9e-2 residual). Tying is now disabled for GFN-FF.
3. **Hydrogen-bridge bonds.** `a *= (1 - 0.1*hb_cn_H)` for bonds with `nr_hb >= 1`, and
   `hb_cn_H` is computed inside `FFWorkspace` at evaluation time — it is 0 in the exported
   parameter set with no accessor. Those bonds (16 of 108 on the peptide) are excluded from
   the fit rather than fitted with the wrong exponent; the remainder still accounts for them.
4. **Gradient units.** `GFNFFComputationalMethod::getGradient()` returns Hartree/**Bohr**,
   unlike the Hartree/Angstrom the rest of curcuma passes around. Differencing it over an
   Angstrom step left every reference Hessian a uniform factor 1.8897 too small.
5. **The coordination numbers never reached the term build** (below).

**Production-path injection.** `GFNFF::setExternalParameterSet()` (forwarded by
`GFNFFComputationalMethod::setExternalGFNFFParameters`) replaces parameter generation with a
supplied set, so a fitted force field is evaluated through the real path — FFWorkspace and
all the state `GFNFFMethod` provides. Both the remainder `H_nb = H_FF(k_current) -
sum k_current*U_p` and the verification use it; a bare `ForceField` + `setGFNFFParameters`
builds no workspace and was measured 90 % off.

#### The dynamic bond r0 makes the bond unit Hessian non-local

GFN-FF's bond reference distance follows the coordination numbers,

```
r0 = (r0_base_i + cnfak_i*cn_i + r0_base_j + cnfak_j*cn_j + rabshift) * ff
```

and the force field carries that through to the gradient as `dE/dcn * dcn/dx`. Since CN
reaches every neighbour, the exact unit Hessian of a bond is **not** confined to its two
atoms. Writing `u = dr/dx` and `w = ff*(cnfak_i*dcn_i/dx + cnfak_j*dcn_j/dx)` (both frozen
at the reference geometry), the bond gradient is exactly `dE/dr * (u - w)`, so

```
U_bond = (local 2-atom block)  -  1/2 * d2E/dr2 * (u w^T + w u^T)
```

Everything is still proportional to the force constant, so the fit remains **one linear
solve** — the correction is stored as the rank-2 factorisation (`TermHessian::extra_*`),
never as a dense 3N x 3N matrix. `w` comes from the frozen `CNDerivStore`, densified once
per run and handed over via `setCNDerivatives()`.

Omitting this term costs 1.5 % on the bond kernel and made the end-to-end verification fail
at 7.7e-04.

#### Why the CN had to be plumbed twice

`setCoordinationNumbers()` and `setCNDerivatives()` **rebuild the term data**. The
constructor builds it first, at which point neither exists, and a GFN-FF bond then falls
back silently to the static `Bond::r0_ij`. That is only ~5e-3 Bohr off the dynamic r0 —
invisible in the longitudinal curvature (0.3 %) but a **2 % error in the transverse
component**, which is proportional to `r - r0` and which GFN-FF's off-minimum bonds
(`dr/r ~ 0.10`) make clearly visible. It was the entire residual after the CN response was
neutralised.

For a multi-basin fit, set `BasinData::cn` / `BasinData::dcn` per basin; both are geometry
dependent and otherwise basin 0's values are reused.

#### Self-check

`gfnff_kernel_check_<kind>` scales one term kind by 5 % and compares the real Hessian change
against the model's. This isolates a single kernel and does not depend on the remainder or
on the fitted constants. Measured:

| kind | CH3OCH3 | peptide (107 atoms) |
|---|---|---|
| bond | 5.3e-08 (8) | 5.2e-08 (92) |
| angle | 1.2e-07 (13) | 9.6e-08 (191) |
| torsion | 1.4e-07 (6) | 1.3e-07 (293) |
| out-of-plane | — | 1.2e-07 (25) |

All at the finite-difference floor. The end-to-end verification (one production Hessian at
the fitted constants vs the fit's own prediction) is 5.9e-09 and 6.4e-09.

### What it does NOT establish

The fit reproduces *a* Hessian; it does not make QMDFF a good force field. On
CH3OCH3 the fit reaches R² = 0.981, i.e. the QMDFF functional form cannot represent
~2% of the GFN2 Hessian's variance — that part is absorbed by whichever term
overlaps most, usually the angles, which makes those constants non-transferable.
See [QMDFF_CONFORMER_USECASE.md](QMDFF_CONFORMER_USECASE.md) for what the fitted
force field is and is not good at on a real conformer ensemble.

## Known limitations and open points

- **No reference-number validation of the energy expression** (see above). This is
  the main open item.
- **Out-of-plane barrier** starts from UFF's improper force constant before the fit,
  not a QMDFF value.
- **`nk` classification** of non-covalent pairs is derived from the bond-graph
  distance (1,4 → class 3; 1,5 → class 4; beyond → class 5). A real QMDFF
  parameter file stores this list explicitly; the interpretation of class 6
  (electrostatics off, dispersion on) has no counterpart in curcuma's generator.
- **`quff`** shares the QMDFF code path (`m_ff_type == 2`) and therefore inherits
  all of the above.
- Legacy parameter files written before this port carry no `qmdff_torsions` /
  `qmdff_ncis` lists; those fall back to the UFF torsion/inversion/LJ terms so
  that old `*.param.json` files keep working.

## Related pre-existing issue (not part of this work)

The **UFF** analytic gradient disagrees with its own finite difference by ~1e-2
Eh/Å on every molecule tested, including H2O and CH4 which have only bond and
angle terms (`test_qmdff_gradients --method=uff`). This is independent of QMDFF
and was found incidentally; it has not been investigated.
