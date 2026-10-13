# WP: Periodic non-bonded interactions (minimum image)

**Status: OPEN.** The periodic container (`wall_potential = pbc`, Sep 2026) relocates
molecules geometrically but the energy is still non-periodic. The clearance test and
the elastic-reflection fallback described below are a **stopgap**, not the intended
solution, and this work package is what replaces them.

## What exists today

`SimpleMD::wrapIntoContainer()` moves a whole fragment back into the container when
its centre of mass leaves it: per axis for a box, to the antipodal point for a
sphere. It exerts no force, so unlike the harmonic and logfermi walls it does not
heat the system and does not do work on an escaping molecule.

## Why it needs a stopgap

The GFN-FF energy has **no minimum-image convention**. `FFWorkspace`, the CN
calculator and the EEQ solver all use plain `r_ij = |x_i - x_j|`; only the unused
legacy `ForceFieldThread` has MIC code. So a wrapped molecule does not continue its
interactions across the boundary — it arrives at coordinates that, as far as the
force field is concerned, are simply somewhere else in the same box.

That has one immediate consequence: **the destination may be occupied.** In a true
periodic system this cannot bite, because the arriving molecule was already
interacting with those neighbours through their periodic images. Here it is a hard
overlap. Measured on 2 N2 + 6 H2 at 3000 K in a 4.5 A sphere: a wrap dropped a
molecule 0.024 A from another atom and the run turned NaN one step later.

The stopgap in `wrapIntoContainer()`:

1. Before wrapping, test the destination against every atom outside the fragment and
   refuse the wrap when anything is within `kWrapClearance` (1.8 A).
2. When the wrap is refused, reflect the fragment elastically instead: reverse the
   outward component of its centre-of-mass velocity. That conserves kinetic energy
   exactly, so it does not heat either, and it keeps the molecule in the container.

It is contained and it does not heat, but it is **not periodic**: a molecule near the
boundary feels no neighbours on the other side, and whether it wraps or bounces
depends on what happens to be standing at the destination.

## What "done" means

Minimum-image distances in the non-bonded terms, so that the container is a real
periodic cell and the wrap becomes a pure bookkeeping operation.

1. **Cell definition.** One cell carried through `Mol`/`Molecule` (`m_unit_cell`,
   `m_has_pbc` already exist) into `GFNFF` (which already copies them at
   `gfnff_method.cpp:757-760`) and on into `FFWorkspace`. The spherical container has
   no lattice and therefore stays outside this: periodic non-bonded interactions are
   a **box** feature, and the sphere keeps the container semantics above.
2. **Minimum image in the pair loops.** Dispersion (D4), repulsion, Coulomb and the
   CN calculation all need the same `r_ij` correction. One shared helper, applied
   where the pair distance is formed, rather than four copies.
3. **Cutoff versus cell size.** MIC requires the interaction cutoff to stay below
   `L/2`. `docs/GFNFF-CPU.md:189` already records this as out of scope; it becomes a
   hard check here, with a clear error when the cell is too small for the cutoffs in
   use.
4. **EEQ.** The charge equilibration is global and its matrix is built from all
   pairs; a periodic treatment needs at least MIC there too, and honestly an Ewald
   sum for the long-range part. Scope this deliberately — MIC-only EEQ is an
   approximation that should be measured against a non-periodic reference before it
   is trusted.
5. **Reactive topology.** `detectReactiveBondChanges()` uses plain distances as well,
   so bonds cannot form across the boundary today. It needs the same helper.
6. **Then remove the stopgap**: with MIC in place the destination is never "occupied"
   in a meaningful sense, so the clearance test and the reflection fallback come out
   of `wrapIntoContainer()` and it reduces to the wrap itself.

## Validation the finished version needs

- Energy of a molecule sitting on the cell boundary equals its energy in the cell
  centre (translation invariance under the periodic image).
- A dimer split across the boundary has the same energy as the unsplit dimer.
- NVE drift over a few ps with a molecule crossing the boundary repeatedly.
- The non-periodic path stays bit-identical when no cell is set.
