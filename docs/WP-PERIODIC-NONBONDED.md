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

## A0 survey (2026-11-01, Claude Generated, read-only)

Done by an agent on `master` (`c3f570fb`); full table with lines in the session scratchpad. Re-checked on
`origin/devel/methods` (`1bd051a5`) only for items 1, 2 and 5. **Line numbers of the other items are from `master` and must be
re-read before coding.** Distances in Bohr unless noted.

1. **CPU reads the cell nowhere.** `GFNFF::InitialiseMolecule` copies `m_has_pbc`/`m_unit_cell` (`ff_methods/gfnff_method.cpp:880-882`);
   only accessors in `gfnff.h:2559-2562` and a copy in `gfnff_frag_charge.cpp:259` use them. `FFWorkspace`, CN, EEQ and D4 have no
   cell member. This is why `unitcell 9 9 9` leaves the CPU energy unchanged.
2. **`SimpleMD::applyPeriodicBoundaryConditions()` wraps every atom separately** into [0,1) (`simplemd.cpp:3447`, called at 3753 and
   4121) once the XYZ carries a cell. That tears molecules, and the bonded terms use raw differences. Needs fragment-wise wrap for a
   box, or MIC in the bonded terms. Open: does the topology recomputation (0.5 Bohr jump) replace the workspace bond list?
3. **Topology assumes a connected geometry** (dense N x N distance matrix, hybridisation, neighbour lists, linearity test); 240 MB at
   5475 atoms. Fragments and phase-1 EEQ depend only on the bond graph.
4. **Six pair lists are built from non-periodic distances:** D4, non-bonded repulsion, explicit Coulomb, HB, XB, CN neighbours.
   `SpatialCellList` (`tools/spatial_cell_list.h`) has no wrap-around; one periodic cell list would cover all six.
5. **Lattice convention mismatch:** `buildLatticeVectors` writes vectors as rows (`tools/pbc_utils.h:60-95`, verified), the MIC helper
   uses `cell_inv * r` (columns). Identical for orthorhombic cells, wrong for triclinic. Scope: orthorhombic only until fixed.
6. **GPU MIC is incomplete:** present in dispersion, repulsion, Coulomb, bonds, angles, HB/XB energies, CN; absent in dihedrals,
   inversions, sTors, ATM, mixed XB, the pair-list builders, HB detection and all EEQ matrix kernels (so EEQ charges are
   computed without MIC and inserted into a Coulomb with MIC). ROCm likewise. (From code reading, not measured.)
7. **Cutoffs (master):** CN 10; CN derivative range max(40, eeq_cut, disp_cut); D4 pair build 60 / evaluation 50 (skin 10); ATM off
   (`dispersion_atm` false); repulsion 20 (+2 skin); Coulomb and EEQ **no cutoff** (`coulomb_r_cut` 0, `eeq_distance_cutoff` 0);
   HB 15.8 / 21.2; XB 21.2; reactive scan 1.6 / 2.6 x covalent radii sum, every 5 calls.
   Against L/2 for 1825 H2O (L = 37.9 A, L/2 = 35.8 Bohr): D4 evaluation (50 Bohr) and the CN derivative range (40) exceed it;
   `dispersion_cutoff_bohr` must be <= 35. The remainder beyond 35 Bohr is not quantified.
8. **D3 vs D4:** D3 in native GFN1 (pair loop without cutoff, CN 25, s9 = 0) and uff-d3 (cutoff 100); the GFN-FF D3 variant aborts with
   an error. D4 in native GFN2 (`D4Evaluator`, pair 50, ATM 25) and in GFN-FF (own variant).
9. **Native GFN1/GFN2 have no MIC** (xtb_h0, xtb_gradient, xtb_response, xtb_precond, xtb_multipole); not part of this work.
10. **Reusable:** `PBCUtils::applyMinimumImage` / `calculateDistancePBC` (header-only). Proposal: one inline orthorhombic MIC in Bohr
    with the same `rint` arithmetic as `applyMIC` on the GPU, so CPU and GPU stay comparable.

## A1 foundation (2026-11-01, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

Done on `feature/bulk` (base `1bd051a5`). **Still no minimum image in any energy term and no Ewald**: a GFN-FF energy
with a cell is the non-periodic energy (ctest `pbc_foundation` checks this bit for bit).

- **Shared MIC** in `src/tools/pbc_utils.h` (MIC inline; check and wrap in `pbc_utils.cpp`, part of `curcuma_core`). Chosen
  over `ff_methods/` because SimpleMD (`capabilities/`) and `Molecule` use the same header: `OrthoCell`, `makeOrthoCell()`,
  `minimumImageOrtho()`, `applyMinimumImageOrtho()`. Arithmetic of the GPU `applyMIC`
  (`s = d*inv_L; s -= rint(s); d = L*s`), `inv_L` from the Eigen inverse like the GPU upload. Unit-agnostic; Bohr in GFN-FF.
- **Cell into `FFWorkspace`**: `GFNFF::configurePeriodicCell()` (start of `InitialiseMolecule()`) builds the Bohr cell,
  `FFWorkspace::setPeriodicCell()` stores it (member only, nothing reads it yet). A non-orthorhombic cell is an error;
  the message names the row/column mismatch of `buildLatticeVectors()` vs `applyMinimumImage()` (not repaired).
- **Cutoff check** (`GFNFF::activeCutoffsBohr()`, `PBCUtils::checkCutoffsAgainstCell()`): every listed radius must be
  `< L_min/2`, else `std::runtime_error` naming each cutoff, its value and `L_min/2` in Bohr. Checked: CN list
  (`cn_cutoff_bohr` 10), CN-derivative range `max(40, eeq_distance_cutoff, dispersion_cutoff_bohr)`, D4 evaluation 50 and
  build 60 (both capped by `dispersion_cutoff_bohr`), repulsion 20 (+ `nonbonded_skin_bohr` from `nb_cell_list_min_atoms`
  on), HB `sqrt(hbthr1)` 15.8 and `sqrt(hbthr2)` 21.2 (also the XB radius), `coulomb_r_cut` if set, reactive scan
  (`topology_mode react`). Coulomb/EEQ without cutoff are not checked (Ewald). Without a cell nothing changes.
- **Consequence**: with the defaults the smallest accepted cube is L > 120 Bohr (63.5 A, D4 build radius). Even with
  `dispersion_cutoff_bohr <= 35` the CN-derivative floor of 40 Bohr needs L > 80 Bohr (42.3 A), so the 1825-water box of
  A0 (37.9 A) is refused. Whether that 40-Bohr floor is needed at all is open (the CN derivative itself ends at 10 Bohr);
  answered in A2: not needed, with a cell it is now the CN range.
- **SimpleMD**: `applyPeriodicBoundaryConditions()` wraps whole fragments (`Molecule::GetFragments()`, covalent-radius
  graph, MIC distances with a cell) by their centre of mass into fractional [0,1) via `PBCUtils::wrapFragmentsIntoCell()`;
  O(N) gate before the O(N^2) fragment search. Systems with CG beads keep the per-bead wrap (2 A bead radius would merge a
  melt into one fragment). A fragment already torn in the input is moved as it is, not reassembled.
- **Not done in A1**: MIC in the pair loops and pair-list builders (A2), EEQ/Coulomb periodicity, triclinic cells, GPU
  parity of the check (the GPU wrapper calls the same `GFNFF` set-up, not verified on a GPU), removal of the
  `wrapIntoContainer()` stopgap. Human production testing pending.

## A2 periodic coordination number (2026-11-01, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `273ff2cd`). **Only the CN is periodic.** Dispersion, repulsion, Coulomb/EEQ, HB/XB, bonded terms,
topology bonds and the D4 generator's own set-up CN still use plain distances, so a periodic GFN-FF energy is not yet
periodic. Without a cell every path is unchanged (checked bit for bit, see below).

- **CN value**: `CNCalculator::calculateGFNFFCN` (both overloads) and `calculateGFNFFCNWithNeighbors` take an optional
  `const PBCUtils::OrthoCell*`; active cell: every pair vector is the minimum image (`applyMinimumImageOrtho`). The
  neighbour-list pass uses a periodic `SpatialCellList` (`calculateGFNFFCNWithNeighborsPeriodic`), neighbours sorted
  ascending (the non-periodic order), thread-count independent. All GFN-FF CPU call sites pass `&m_cell_bohr_ortho`:
  per-step CN and its N^2 fallback, `updateDynamicState`, topology-init CN, angle CN, `mchar` (incl. its pair loop),
  `calculateTopologyCharges`, `NumGradFixedCharges`.
- **CN derivative**: `calculateCoordinationNumberDerivatives` (every branch: threaded/serial, `GFNFF_FAST_EXP` or not,
  neighbour list or threshold, the cn_raw recompute) and the fused path (`GFNFF_CN_DCN_FUSION`, off in `release/`, not
  compiled or tested) take r_ij = x_j - x_i as minimum image, so both pair entries and the diagonal refer to the nearest image.
- **40-Bohr floor**: the derivative support is the CN support (dCN_i/dx_k = 0 unless k is a CN partner of i). With the
  neighbour list (default `cn_cutoff_bohr` 10) the threshold is not read at all; the "stencil must cover the dispersion
  range" comment at that site is wrong (the reach of j enters dE/dCN_i, not dCN_i/dx_k). With a cell the range is now
  `cnDerivativeRangeBohr()` = `cn_cutoff_bohr` (40 for `cn_cutoff_bohr` 0, matching the 1600 Bohr^2 per-step CN there);
  without a cell `max(40, eeq_distance_cutoff, dispersion_cutoff_bohr)` as before. The L/2 check lists the periodic value.
- **Smallest accepted cube with defaults** is still L > 120 Bohr (D4 build 60); next blockers repulsion 22, HB/XB 21.2.
- **Topology-init CN** keeps its 40-Bohr threshold with the minimum image; for L/2 < 40 pairs beyond L/2 count with their
  nearest image only (their erf contribution is below 1e-20 for the pairs measured; not checked for heavy elements).
- **Tests** (`ctest -R pbc_`): `pbc_cn` (cell list vs brute-force MIC on 4 grids; GFN-FF CN/dCN of the water dimer with one
  molecule shifted by L and/or torn across the boundary equal to the intact dimer without cell to 4.4e-16 / 3.6e-15;
  small (8,8,9) Bohr cell with cut 3.9 equal to a 27-image brute-force sum, 1.3e-15, images change CN by 0.81;
  2x1x1/1x1x2 supercells equal per atom; 384-atom 4x4x4 replica equal to the single cell and bitwise for 1 and 4 threads;
  analytic dCN vs central FD (Angstrom) for the torn dimer 4.7e-8 at h = 1e-4 A, quadratic in h). `pbc_foundation`'s
  100 A cell check is now 1e-12 Eh, not bitwise (MIC arithmetic L*(d/L) is not exact).
- **Not tested**: total energy or gradient with a cell (not periodic yet), CN with charge groups/react mode under a cell,
  GPU paths. Without a cell: -sp energy (15 decimals) and full gradient dump (17 digits, md5) identical before/after for
  benzene, caffeine, methyl-glucopyranoside, water8, hydronium-water2, urea-water197, peo201-water500 at 1 and 4 threads.
- **GPU (code reading only)**: `k_cn_compute`/`k_cn_chainrule` use `applyMIC` with the same `rint` arithmetic (the
  off-diagonal products are exact zeros), so CPU and GPU can agree to rounding; differences: GPU pair test `r2 > thr`
  (inclusive) vs CPU strict `<`, radii `d_rcov_d3 * 4/3` vs `gfnff_cn_rcov_bohr()`, fixed-point CN sums. Not run.
- Human production testing pending.

