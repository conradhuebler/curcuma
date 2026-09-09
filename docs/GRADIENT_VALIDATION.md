# Gradient validation: curcuma vs xtb, set-wide

⚠️ AI-generated, machine-tested only — not human production tested.

Until September 2026 curcuma's **energies** were validated over thousands of structures
(GMTKN55, MOR41, S30L-CI) while its **gradients** were only spot-checked: three small
molecules in `ctest` plus ad-hoc finite differences. This document closes that gap.

Runners:
* `scripts/gradient_compare.py` — set-wide `max |g_curcuma − g_xtb|` per structure.
* `scripts/gradient_arbitrate.py` — for one structure, decides **which** code is right by
  finite-differencing each program's *own* total energy.
* `test_cases/check_gradient_units.py` — the permanent `gradient_unit_contract` ctest.

```bash
python scripts/gradient_compare.py --set gmtkn55 --method all
python scripts/gradient_arbitrate.py test_cases/GMTKN55-testset/BH76/hcnts/struc.xyz gfn1
```

## Units — the thing to get right

`EnergyCalculator::Gradient()` returns **Eh/Angstrom** for every method. `SimpleMD` feeds
it straight into a Verlet integrator whose coordinates are Angstrom, the LBFGS interfaces
take it raw, and `Hessian` differences it against Angstrom displacements. Native xTB
converts explicitly (`xtb_native.cpp:1480`, `m_gradient /= au`); GFN-FF does so since
Sep 2026 (`gfnff_method.cpp`). xtb's Turbomole `gradient` file is Eh/Bohr, so the compare
script multiplies curcuma by `au`. Everything below is Eh/Bohr.

## Result (GMTKN55, 2462 structures, xtb 6.7.1)

Metric: the largest single Cartesian component of `|g_curcuma − g_xtb|` per structure.

| method | n | median | max | ≤1e-6 | ≤1e-5 |
|---|---:|---:|---:|---:|---:|
| gfn1  | 2462 | 3.0e-07 | 4.8e-02 | 86.1 % | 95.2 % |
| gfn2  | 2458 | 3.4e-07 | 4.8e-02 | 82.6 % | 95.3 % |
| gfnff | 2451 | 3.7e-08 | 6.5e-02 | 76.2 % | 81.4 % |

(gfn2 was 4.3e-07 / 79.7 % before the SCF-convergence fix of Known Issue #29.)

(gfn2: 4 structures where one code produced no gradient. gfnff: 9 such, plus
`PX13/hf_4_ts` and `hf_6_ts` where **xtb** returns NaN and curcuma returns finite values —
the same two structures Known Issue #24 records as NaN in pprcht and xtb.)

The distribution is bimodal: most structures agree at the numerical floor, a minority
disagree by a lot. **A deviation from xtb is not by itself evidence against curcuma**, so
every large one below was arbitrated with finite differences.

## Arbitration: who is right where they differ

A correct analytic gradient must reproduce the central finite difference of the energy it
belongs to. That referee needs no third code.

### Every outlier was FD-checked, not just a sample

For **all** gfn1/gfn2 structures deviating from xtb by more than 1e-4 Eh/Bohr, curcuma's
analytic gradient was compared against the finite difference of curcuma's **own** energy at
the largest gradient component. Since the energies match xtb (MAD 0.000 kcal/mol), that is
the port criterion. First pass: **14 of 202 failed — all H2/H2+ with gfn2**, a real curcuma
bug (Known Issue #29: the SCF convergence test ignored the multipole moments the mixer was
already mixing). After the fix: **0 of 188 fail**, worst relative deviation 3.4e-03, which
is FD truncation. Extrapolating from the five hand-arbitrated cases below would have missed
it — they all happened to be xtb-side.

### gfn1 / gfn2 — curcuma matches the FD of *both* energies, xtb does not

Five independent structures, at 2–3 step sizes each:

| structure | method | curcuma analytic | xtb analytic | FD of curcuma E | FD of xtb E |
|---|---|---:|---:|---:|---:|
| `ACONF/B_G` atom 0 x | gfn2 | +0.001670 | −0.005970 | +0.001670 | +0.001669 |
| `ACONF/H_ttt` atom 0 x | gfn2 | +0.002099 | −0.006336 | +0.002099 | +0.002100 |
| `WCPT18/ts6` atom 0 x | gfn1 | +0.030772 | −0.017099 | +0.030772 | +0.030771 |
| `BH76/hcnts` atom 1 x | gfn1 | +0.004213 | +0.022244 | +0.004213 | +0.004213 |
| `PNICO23/2` atom 0 x | gfn1 | +0.012554 | −0.002066 | +0.012553 | +0.012554 |

The two codes' **energies agree** (their finite differences coincide), and curcuma's
analytic gradient reproduces them; xtb's does not. The affected geometries share a
signature: **atoms with exactly equal Cartesian coordinates** — `BH76/hcnts` has C and N at
identical x, the ACONF alkanes several carbons at identical x — which is common in
idealised benchmark geometries and rare in relaxed ones. Consistently, at xtb's *own
optimised* geometry the two agree again (H2O gfn2: |g| 5.40e-5 vs 5.47e-5 Eh/Bohr, largest
component +0.000027 both). Not root-caused inside xtb; the claim is scoped to these tests.

### gfnff — both codes are right about their own energy; the energies differ

| structure | curcuma analytic | xtb analytic | FD of curcuma E | FD of xtb E |
|---|---:|---:|---:|---:|
| `AL2X6/al2me6` atom 7 z | −0.005748 | −0.070667 | −0.005748 | −0.070667 |
| `MB16-43/35` atom 14 y | +0.010352 | −0.052262 | +0.010352 | −0.052262 |
| `INV24/NCl3_TS` atom 3 y | +0.024390 | +0.059154 | +0.024390 | +0.059154 |

Each analytic gradient reproduces the finite difference of its **own** energy exactly, and
the two finite differences differ. So the GFN-FF gradient deviation is the *energy* split
between the pprcht port source and xtb — `AL2X6/al2me6` and `MB16-43` are exactly the
structures Known Issues #21/#22 name as that split — not a gradient defect.

## Two bugs this validation found (both fixed, Sep 2026)

See CLAUDE.md Known Issue #28 for the full account.

1. **GFN-FF returned Eh/Bohr where the interface contract is Eh/Angstrom**, so MD forces
   were a factor 1/au = 1.8897 too small. Found by NVE energy conservation, which
   separates the causes cleanly — an integration error shrinks as dt², a force/energy
   inconsistency does not:

   | dt / fs | GFN-FF before | GFN-FF after | gfn2 (control) |
   |---|---:|---:|---:|
   | 0.5 | 3.4e-03 | 7.6e-05 | 7.7e-05 |
   | 0.25 | 3.5e-03 | 1.9e-05 | 2.2e-05 |
   | 0.125 | 3.5e-03 | 5.0e-06 | 6.0e-06 |

   (Etot drift, CH4, 200 fs, thermostat off.) The existing `gfnff_numgrad_builtin` test
   could not see it: it compares GFN-FF's *internal* gradient against an *internal* finite
   difference, both in Bohr.

2. **The finite-difference Hessian was assembled in Eh/Ang² and consumed as Eh/Bohr²**, with
   an empirical `+ 47.349 cm⁻¹` offset masking part of it. Every frequency was too high —
   gfn1/gfn2 by 1/au, gfnff by sqrt(1/au) (one factor cancelled against bug 1). After both
   fixes, H2O against xtb 6.7.1:

   | method | curcuma | xtb 6.7.1 | max dev |
   |---|---|---|---:|
   | gfn1 | 1490.4 / 3515.4 / 3621.4 | 1490.57 / 3515.28 / 3626.07 | 0.13 % |
   | gfn2 | 1574.9 / 3486.3 / 3487.0 | 1574.99 / 3486.05 / 3491.64 | 0.13 % |
   | gfnff | 1631.9 / 3635.0 / 3637.6 | 1632.01 / 3634.77 / 3637.31 | 0.01 % |

Energies and optimised geometries were never affected by either.

## Are MD and geometry optimisation safe now?

Measured, not asserted.

**MD — yes, for all three methods.** NVE total-energy drift on CH4 / H2O / NH3, thermostat
off, 200 fs:

| molecule | method | dt 0.5 fs | dt 0.25 | dt 0.125 | drift ratio (dt^2 => ~4) |
|---|---|---:|---:|---:|---|
| CH4 | gfn1 | 8.4e-05 | 2.5e-05 | 6e-06 | 3.4 / 4.2 |
| CH4 | gfn2 | 7.7e-05 | 2.2e-05 | 6e-06 | 3.5 / 3.7 |
| CH4 | gfnff | 7.6e-05 | 1.9e-05 | 5e-06 | 4.0 / 3.8 |
| H2O | gfn1 / gfn2 / gfnff | 1.5-1.8e-04 | 3-5e-05 | 0.9-2.1e-05 | 3.3-4.7 / 2.5-3.4 |
| NH3 | gfn1 / gfn2 / gfnff | 2.6-2.8e-04 | 5.8-6.6e-05 | 1.5-1.8e-05 | 4.0-4.9 / 3.2-4.4 |

All nine combinations scale as dt^2, i.e. the forces are consistent with the energy. Before
the fix GFN-FF was flat in dt (3.4e-03 at every step size).

**Optimisation — yes, with one caveat that is not about gradients.** Final minima against
`xtb --opt tight` on the sqm_reference molecules, all three methods: energies agree to
**<= 0.0085 kcal/mol** and geometries to <= 0.0007 A RMSD for rigid molecules
(`acetic_acid_dimer` differs by 0.04-0.19 A on a flat, floppy surface while its energy still
agrees to 0.0085 kcal/mol; that RMSD has no rotational fit and is an upper bound).
Robustness over 30 randomly drawn GMTKN55 structures x 3 methods: **89/90 converge**; the
single failure is `DIPCS10/nh3_2+` with GFN-FF, a dissociative dication.

The caveat: the default **`auto` optimiser can abort with "Energy rise exceeded maximum
allowed"** on a high-symmetry system whose only active mode is the totally symmetric one.
Reproducible on perfect-Td CH4 with gfn1 (initial |g| = 9.0e-03 Eh/Ang, so genuinely not
converged) — and it is **pre-existing**, a pre-fix binary fails identically, so it is a
step-control issue and not a gradient issue. `-opt.optimizer lbfgs` optimises the same case
correctly: E = -4.27427100 vs xtb's -4.27427080 Eh, C-H 1.08609 vs 1.08580 A.

## What was tested / not tested

- **Tested**: analytic gradients of gfn1/gfn2/gfnff against xtb 6.7.1 on all 2462 GMTKN55
  geometries at identical charge/multiplicity; the largest disagreements arbitrated against
  finite differences of both codes' energies; the gradient unit contract and H2O
  frequencies as a permanent ctest.
- **Not tested**: MOR41 and S30L-CI gradients (the script supports `--set mor41`/`s30lci`,
  they have not been run); solvated gradients beyond the existing `xtb_solvation_numgrad` /
  `gfnff_numgrad_alpb_water` tests; GPU gradient paths; gradients of any method other than
  these three (UFF, EHT and QMDFF behave differently in a first probe and were not pursued
  — see the note below).
- **Open**: a first unit probe on `-method uff`, `eht` and `qmdff` gave gradient/FD ratios
  of 10.8, ~0 and 1.02 on H2O, i.e. none of them clearly honours the Eh/Angstrom contract.
  Not investigated; the `gradient_unit_contract` ctest deliberately covers only the three
  validated methods.

Human production testing pending until the operator removes this note.
