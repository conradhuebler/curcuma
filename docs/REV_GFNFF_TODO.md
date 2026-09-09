# rev-gfnff — where curcuma could be better than the reference

⚠️ AI-generated, machine-tested only — not human production tested.

A running list of places where curcuma's GFN-FF **deliberately does, or could do, something
different from the reference** because the reference is demonstrably wrong on the physics.
Everything here is separate from port fidelity: the default is to reproduce
`pprcht/gfnff` (and, where they agree, `xtb`), because that is what the GFN-FF parameters
were fitted against. This file collects the cases where that costs accuracy, so a future
revised force field can pick them up deliberately rather than rediscover them.

**The distinction that decides each case**: were the GFN-FF parameters fitted *against* the
quirk? If yes, "fixing" it invalidates the parametrisation and must stay opt-in. If the
quantity is a hard-coded ab-initio constant or a plain coding slip, curcuma may lead.

## Already deviating from the reference by default

### 1. `nh_linear_fix` — GEODEP angle rule creates artefact minima at N-H centres
Started on branch `confsearch` (`51830efa`), ported and then refined here (see
[GFNFF_STATUS.md](GFNFF_STATUS.md)). The reference's "input angle > 160 deg -> sp" fallback
makes a thermally stretched =N-H its own equilibrium, so an optimisation walks into the
artefact. **External verdict** on the same two geometries (formamidine, only the imine
C-N-H angle varied), ΔE(179° − 119°):

| method | kcal/mol |
|---|---:|
| r²SCAN-3c | +24.1 |
| GFN2 | +24.0 |
| GFN-FF, guard on (curcuma default) | +14.0 |
| GFN-FF, guard off (= xtb/pprcht) | **−66.4** |

The reference inverts the sign; curcuma keeps it. `-gfnff.nh_linear_fix false` restores
bit-faithfulness. Status: **done, default on.**

### 2. `storsion_reference_loop_bug` — sTors evaluates only its last entry
`gfnff_engrad.F90:494-503` loops `do i=1,m` but calls `sTors_eg(m, ...)` with the array
SIZE, so pprcht *and* xtb evaluate one entry m times and drop the rest. Its `erefhalf` is a
DLPNO-CCSD(T) diphenylacetylene value, **not a fitted parameter**, so summing correctly
cannot invalidate anything. curcuma sums; `-gfnff.storsion_reference_loop_bug true`
reproduces the reference. Status: **done, default correct.**

## Open — physics known to be wrong in BOTH implementations

### 3. EEQ over-rewards charge delocalisation (charged radicals, 2c-3e bonds)
Measured on GMTKN55 `G21EA/EA_25`, the dichlorine radical anion Cl₂⁻ at 2.73 Å.
E(Cl₂⁻) − E(Cl⁻) − E(Cl):

| method | kcal/mol |
|---|---:|
| r²SCAN-3c (ORCA 6, UHF) | **−41.5** |
| GFN2 (xtb, UHF) | −33.9 |
| GFN-FF, nfrag = 2 (reference behaviour, charge pinned on one Cl) | −6.6 |
| GFN-FF, nfrag = 1 (charge delocalised over both) | −106.3 |

Experiment puts the dissociation near −30 kcal/mol. **Neither GFN-FF setting is right**, and
the fragment count is the only lever — a discrete switch that brackets the truth without
hitting it. Term decomposition (curcuma):

| | bond | Coulomb | total |
|---|---:|---:|---:|
| nfrag = 2 | −11.5 | +5.3 | −6.6 |
| nfrag = 1 | −11.2 | **−94.8** | −106.3 |

The bond term is essentially unchanged; the entire difference is the EEQ. Going from
q = (−1, 0) to (−0.5, −0.5) is worth **−94.8 kcal/mol** to GFN-FF's EEQ self-energy
(which scales as q², so splitting one unit charge over two centres halves the penalty),
where the true stabilisation is ~−35. A revised model needs a delocalisation term that is
not just the EEQ self-energy — a **method change**, not a parameter or perception fix.

Curcuma currently reproduces the reference here (Known Issue #17: the second q-loop pass
keeps pass 1's fragmentation, `gfnff_ini.f90:467`), because that is both faithful and the
smaller error. Affects the whole charged-complex family: AHB21, BH76 SN2 transition states,
G21EA, SIE4x4, CHB6.

### 3b. Carbene angle: GFN-FF bends CH2- 30 deg wrong either way
Bending scan of the methylene anion CH2- (C-H fixed at 1.129 A, H-C-H varied), relative
energies in kcal/mol:

| H-C-H | GFN-FF carbene (θ0=145) | GFN-FF reference (θ0=120) | GFN2 | r²SCAN-3c |
|---:|---:|---:|---:|---:|
| 90 | 13.77 | 10.42 | 2.56 | 1.58 |
| **100** | 8.87 | 5.72 | **0.00** | **0.00** |
| 110 | 5.27 | 2.62 | 1.66 | 1.14 |
| 120 | 2.77 | 0.81 | 6.65 | 4.57 |
| 130 | 1.17 | **0.00** | 14.27 | 9.72 |
| 145 | **0.00** | 0.02 | 28.81 | 19.02 |

GFN2 and r²SCAN-3c agree on a minimum at **100°** (experiment ~102°). The reference-faithful
treatment minimises at 130°, the carbene one at 145° — so restoring port fidelity (Known
Issue #18) also moved GFN-FF *towards* the truth, unlike case 3. But 130° vs 100° is still a
30° error: GFN-FF has no equilibrium angle that fits a bent carbanion. Unlike case 3 there is
no conflict here, only a ceiling.

### 4. GFN-FF cannot do MOR41 reaction thermochemistry at all
Against DLPNO-CCSD(T) (Table S1, 41 reactions): pprcht MAD 62.6, curcuma 63.1, xtb 71.7
kcal/mol, where GFN2 reaches ~12. A method limitation, not a port issue — recorded so nobody
re-opens it as a bug. See [MOR41_VALIDATION.md](MOR41_VALIDATION.md).

### 5. Bond perception in charged species is chemically right but energetically worse
Same root as (3), stated separately because it may deserve its own treatment: curcuma's
charge-shrunk radii correctly perceive Cl₂⁻ as **bonded** (it is — a real radical anion at
~2.7 Å), and the reference's own second pass agrees (`#bonds : 1`). But GFN-FF's bond
parameters are fitted for ordinary 2-electron bonds and overbind a 2-centre-3-electron bond.
Getting the perception right therefore does not help until the bond term can express a
half bond.

### 6. The angle term is ~3x too soft at an SN2 transition state
Umbrella scan of the F...CH3...F- transition state (C-F and C-H fixed, the three F-C-H
angles varied together), energy relative to the D3h point:

| F-C-H | GFN-FF | GFN2 | r²SCAN-3c |
|---:|---:|---:|---:|
| 80° | **2.13** | 6.03 | **6.43** |
| 85° | 0.48 | 1.49 | 1.61 |
| 90° | 0.00 | 0.00 | 0.00 |

GFN2 reproduces r²SCAN-3c almost exactly; GFN-FF is a factor of three too soft. Known Issue
#19 restored port fidelity by removing a non-reference `|qa| < 1` guard around the angle
`fqq`, which **softens** the F-C-H constants further (0.2631 -> 0.224) — i.e. away from the
physics, like case 3. It was still the right call: the accidental stiffening only ever
applied to angles at a full-unit-charge fragment, so it was not a principled correction.
A revised model would need a stiffer bend at a hypercoordinate centre generally.

*Caveat on the numbers*: the scan geometry is symmetric, so its perception differs slightly
from the GMTKN55 structure (there one fluoride is its own fragment). The GFN-FF column is
therefore the reference-faithful `fqq` regime throughout — the comparison measures the term,
not the guard.

## How to add to this list

An entry belongs here when an **external** reference (r²SCAN-3c, GFN2, DLPNO, experiment)
has been consulted and shows the reference implementation to be wrong — not merely different
from curcuma. Record the numbers, not the impression.
