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
  `cnDerivativeRangeBohr()` = `cn_cutoff_bohr` (for `cn_cutoff_bohr` 0 sqrt(cnthr), see the next item);
  without a cell `max(40, eeq_distance_cutoff, dispersion_cutoff_bohr)` as before. The L/2 check lists the periodic value.
- **`cn_cutoff_bohr` 0 unified (resolved)**: the per-step N^2 CN used a fixed 1600 Bohr^2, `updateDynamicState` cnthr = 100 - 50 log10(`cn_accuracy`); all now use cnthr via `CNCalculator::gfnffAccuracyThresholdSqBohr` / `GFNFF::cnThresholdSqBohr` (reference gfnff_param.f90:551, gfnff_engrad.F90:299; the reference accuracy is 0.1, i.e. 150 Bohr^2, our `cn_accuracy` default 1 gives 100); default path bitwise unchanged on the 7 structures, `cn_cutoff_bohr` 0 bitwise unchanged on caffeine, water8, urea197.
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


## A3b periodic D4 dispersion pair term (2026-11-01, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `674d5a27`). **CN (A2) and the D4 pair term are periodic; nothing else.** Repulsion, Coulomb/EEQ,
HB/XB, bonded terms, topology bonds and ATM still use plain distances; the long-range tail followed in A3c (below).

- **Kernel** `FFWorkspace::calcD4Dispersion`: pair vector as minimum image from `m_cell_bohr` (no-op without a cell);
  energy, force and the dC6/dCN terms use it. **Generator** `D4ParameterGenerator::setPeriodicCell()`: CN with minimum image,
  pair list `i < j` at the minimum-image distance, each pair once, never images.
- **Radii, one source**: `GFNFF::dispersionRadiiBohr()` (eval, build, capped) is read by `applyDispersionCutoff`,
  `dispersionSkinBohr`, `activeCutoffsBohr`, the generator and the verbosity-2 table. Without a cell: 50 / 60, capped by
  `dispersion_cutoff_bohr` as before. With a cell and no user cutoff the evaluation radius becomes `nextafter(L_min/2, 0)`
  when 50 is not below it, the list radius eval + 10; verbosity 1 prints "D4 evaluation cutoff reduced from X to Y Bohr
  (L/2); long-range tail correction not yet applied". An explicit cutoff >= L_min/2 is still an error.
- **Why the build radius may exceed L/2** (comment in `activeCutoffsBohr`): the list stores pairs, the kernel forms the
  minimum image at every call; with eval < L/2 only the nearest image can lie within eval, and no self-image. A pair absent
  from the list cannot reach eval before an atom moved (build - eval)/2, because the minimum-image distance is 1-Lipschitz.
  With a cell only the evaluation radius is checked against L/2.
- **Verlet trigger**: `maxAtomDisplacement(a, b, cell)` takes each displacement as minimum image, so a fragment wrap by
  SimpleMD is no motion. The plain version was safe (it over-estimates) but rebuilt the list after every wrap.
  Valid while no atom moves more than L/2 - skin/2 between two calls.
- **Smallest cube with defaults** is now set by HB/XB 21.2 Bohr (L > 42.4 Bohr; water dimer: 23.0 A accepted, 22.0 A
  refused), with the repulsion cell list (`nb_cell_list_min_atoms`, large systems) by repulsion 22 (L > 44); D4 no longer blocks.
- **Tests** (`ctest -R pbc_dispersion`, water molecules from `water-dimer.roo2.8`, intact, moved rigidly): dimer shifted
  by L / contact across the boundary = no cell, 6.5e-19 Eh; 21 random waters (63 atoms, L 24 Bohr, r_c 11.5, hbond and
  repulsion off) = 27-image brute force 3.5e-18 (images -6.6e-3 Eh), no pair missing, threads 1/4 3e-17; 2x1x1 supercell
  = 2 x cell 1.4e-17; analytic gradient (incl. dE/dCN, max 1.6e-5 of max 1.6e-4 Eh/A) vs central FD for 6 atoms 0.06-0.4 A from a face
  1.1e-11 Eh/A at h = 1e-4 A, 4.3e-11 at 2e-4 (quadratic); default cutoff reduced to 12 Bohr and = brute force 1.7e-18;
  explicit 12.5 / 50 refused; wrap by +L: no rebuild, dE 0; small moves = fresh instance; 5.7 Bohr move rebuilds.
- **Without a cell**: -sp energy (16 digits) and `-dump_gradient` md5 identical before/after on the 7 structures of A2 at 1
  and 4 threads.
- **GPU (code reading, not changed)**: `k_dispersion` applies `applyMIC` with the same arithmetic and an inclusive
  `rij <= r_cut` (CPU skips `rij > r_cut`: same). Missing: the device pair builder `k_disp_pairs_build` uses plain
  distances, the 60 Bohr constant and a hard-coded `r_cut = 50` (neither the cap nor the minimum image), and
  `k_check_displacement` is not wrap-safe. With `gpu_disp_pairs_on_device` (default) a cell therefore gives a wrong
  dispersion energy, now also for small cells that A1 used to refuse. The host-built list (`gpu_disp_pairs_on_device false`)
  carries the MIC pairs and the capped `r_cut`, not run.
- **GPU: fallback to the CPU with a cell until the port** (Nov 1, 2026): `MethodFactory::periodicCellNeedsCpu()`, called in
  `EnergyCalculator::setMolecule()` before the device setup (all backends, also the multi-GPU batch) and for the deltaml base:
  warning "periodic cell set: GPU path has no minimum image yet, using CPU", counted fallback, `-gpu_strict true` exits 3.
  Native gfn1/gfn2 untouched (CPU and GPU both ignore the cell). ctest `pbc_gpu_guard` (stand-in wrapper, no device).
- **ATM (off by default)**: triples and kernel without minimum image; three independent minimum images close the triangle
  only if r_ij + r_jk < L/2, so with the 25 Bohr ATM cutoff a cell needs L > 100 Bohr or a per-triple image choice.
  The bonded ATM (batm, 1,4 pairs) is a bonded term and shares their state (fine for intact molecules, wrong for torn ones).
- Human production testing pending.

### Cutoff dependence measured in A3a (clusters, no cell; input for A3c)

E(r_c) - E(50) in kcal/mol, `-sp gfnff`, E_tot change = E_disp change; r_c = 50/60/0 bit-identical (eval radius min(50, r_c)).

| system (N) | E_disp(50) Eh | 20 | 25 | 30 | 35 | 40 |
|---|---|---|---|---|---|---|
| urea197 (197) | -0.1612275 | 0.158 | 0.0081 | 6e-5 | 0 | 0 |
| peo500 (2910) | -2.1199138 | 14.08 | 5.29 | 2.20 | 0.943 | 0.390 |
| urea1000 (6200) | -6.4752113 | 52.61 | 21.36 | 9.53 | 4.371 | 1.909 |
| peo1500 (7320) | -4.7189039 | 33.34 | 13.19 | 5.81 | 2.674 | 1.188 |

- log-log slope of E(r_c) - E(inf) over 30-40 Bohr: -4.7 (urea1000), -4.3 (peo1500), not the bulk -3, because cluster
  shells beyond 30 Bohr are only partly filled; effective <C6 zeta> per pair is distance-independent (urea1000 6.25-6.33 au).
- C8 share of the tail: 6.6 % (20-25 Bohr shell), 3.1 % (30-35), 2.3 % (35-40).
- Extrapolation to the 1825-H2O box (r_c 35.8 Bohr), **not measured**: E_tail ~ -14.7 kcal/mol (~0.4 % of E_disp,
  <C6 zeta> = 6.3 au assumed), pressure ~ -37 bar.
- Recommendation for A3c: variant (ii), tail per step from the current CN through the per-atom factorisation
  sum_i sum_j zeta_i zeta_j C6_ij (O(N K 7), no pairs), dE_tail/dCN into `acc.dEdcn`, virial dE/dV = -E_tail/V.

## A3c homogeneous long-range D4 tail (2026-11-01, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `4801aab1`), CPU only, active only with a cell (`-gfnff.dispersion_tail`, default true, no effect
without a cell). Variant (ii) of A3a: per call from the current CN, force through dE/dCN.

- **Formula** (`D4ParameterGenerator::homogeneousDispersionTail`, derivation in the comment there): beyond the evaluation
  radius r_c (`dispersionRadiiBohr().eval`, < L_min/2) g(r) = 1, i.e. every atom j of the cell has images with density 1/V:
  E_tail = (1/2) sum_i sum_j (1/V) int_{r_c}^inf 4 pi r^2 e_ij dr = -(2 pi/V) [S6/(3 r_c^3) + 2 S8/(5 r_c^5)],
  S6 = sum_ij zeta_i zeta_j C6_ij, S8 = sum_ij zeta_i zeta_j r4r2ij C6_ij over all N^2 ordered pairs incl. i = j (self images).
  s6 = 1, s8 = 2 as in the kernel. Factorised per partner element and reference (O(N K 49), no pairs); C6 from the current
  Gaussian weights, zeta fixed (topology charges); dE/dCN_k = 2 zeta_k sum_j f_kj dC6_kj/dCN_k into the CN chain rule.
- **Damping (large-r form, R0 dropped)**: over-binds the C6 part by at most (1/3)(R0/r_c)^6, the C8 part by (5/13)(R0/r_c)^8
  (R0 O-O 7.41 Bohr): 3.0e-5 at r_c 35, 2.1e-2 at 11.5 (measured on a single species: 2.9e-5 and 2.1e-2 of the tail).
  At small r_c this is the dominant error: keeping the damping in the same g(r) = 1 integral (diagnostic in the test, per
  element pair by quadrature, same cost) reduces the r_c dependence below from 8.61 / 0.88 to 0.99 / 0.003 kcal/mol.
- **Reported** as "Dispersion tail" (verbosity 2 table, only when non-zero), `FFEnergyComponents::disp_tail`,
  `GFNFF::DispersionTailEnergy()`, JSON key `DispersionTail`; `DispersionEnergy()` stays the pair term. The A3b warning now
  appears only with the tail off. Not implemented: virial/pressure, dE_tail/dV = -E_tail/V at fixed r_c (needed for NPT).
- **Tests** (`ctest -R pbc_dispersion_tail`, waters from `water-dimer.roo2.8`): prefactor vs radial quadrature 6e-16 /
  2e-15; factorised S6/S8 vs the explicit O(N^2) sum over `getChargeWeightedC6` 5.4e-14 / 7.4e-15 (63 atoms); central FD of
  pair + tail (6 atoms 0.05-0.26 A from a face) 3.1e-11 -> 7.8e-12 Eh/A at h 2e-4 -> 1e-4 A, tail alone 3.2e-12 -> 7.8e-13
  (tail force max 1.5e-6 of max dispersion force 1.8e-4 Eh/A); no cell or tail off: tail 0, energy and gradient bitwise
  equal; 2x1x1 supercell E_tail = 2 x cell 2.2e-15. `pbc_dispersion` and `pbc_foundation` now set `dispersion_tail false`
  (pair-term checks; the dimer tail in the 100 A cell is -4e-10 Eh).
- **r_c dependence** (100 waters at 0.92 x liquid density, random sequential placement, L 28 Bohr, NOT a liquid):
  E_pair + E_tail at r_c 8 / 10 / 12 Bohr relative to 12: -8.61 / -0.88 / 0 kcal/mol, E_pair alone moves 35.1 kcal/mol.
- **Replica check** (`test_pbc_dispersion_tail --measure`): small cell r_c = L/2 + tail vs its 2x2x2 replica r_c = L + tail,
  per small cell. Random packing (300 atoms, L 28 Bohr): explicit shell L/2..L -8.114 kcal/mol, tail difference -8.149, so
  E_pair + E_tail differ by 0.035 kcal/mol (1.2e-4 per atom; tail -9.22). Cubic water lattice 5^3 (spacing 3.1 A, L 29.3
  Bohr): shell -8.392 vs tail -9.631, difference 1.24 kcal/mol (a lattice has no g(r) = 1 at 14.6 Bohr). No periodic liquid
  water box exists in the repository (structures/bulk are clusters without a cell); none was built.
- **1825-H2O extrapolation of A3a checked on a lattice**: 12^3 = 1728 waters on the 3.1 A lattice (liquid density, r_c 35.15
  Bohr): E_tail -10.14 kcal/mol (-1.96e-3 per atom, 0.37 % of E_pair, C8 share 1.6 %), against the extrapolated -14.7. The
  difference is the assumed <zeta^2 C6> = 6.3 au (urea/water); water gives 4.26 au. At fixed density E_tail at r_c = L/2 does
  not depend on N, so 1728 and 1825 compare directly.
- **Without a cell**: -sp energy (16 digits) and `-dump_gradient` md5 identical before/after on the 7 structures of A2 at 1
  and 4 threads. **GPU**: not changed; with a cell the A3b guard runs GFN-FF on the CPU, so the tail is always the CPU one.
- **Not tested**: real liquids, other elements than H/O, charged systems (zeta from topology charges), MD energy
  conservation with the tail, the P1a CN-cache threshold and static-CN mode (the tail follows the same stale weights as the
  pair term), frag_charge_model variants with a cell. Human production testing pending.

## A4a periodic repulsion and bonded terms (2026-11-01, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `34657480`), CPU only. **Periodic now: CN (A2), D4 pair term and tail (A3b/c), bonded and non-bonded
repulsion, bonds, angles, torsions (incl. extra torsions), inversions, sTors.** Still plain distances: HB/XB (incl. the HB
modulation of X-H bonds, `computeHBCoordinationNumbers`), Coulomb/EEQ, topology generation (distance matrix), reactive
topology, BATM/ATM; these are A4b and later.

- **Repulsion** (`calcBondedRepulsion`, `calcNonbondedRepulsion`): pair vector as minimum image, as in `calcD4Dispersion`.
  The distance-filtered list (`nb_cell_list_min_atoms`) is built with the periodic `SpatialCellList` (pairs i < j, each once,
  never images); the skin trigger uses `maxAtomDisplacement(a, b, cell)`, so a fragment wrap is no motion. The all-pairs
  list below the threshold needs no change.
- **Bonded terms**: with a cell each term is evaluated on a local copy of its atoms (`placeMinimumImage`,
  `placeChainOfFour` in `ff_workspace_gfnff.cpp`): bond j next to i; angle i, k next to the centre j; torsion and sTors
  i, k next to j, l next to k; inversion the three neighbours next to the centre. The angles therefore come from
  minimum-image bond vectors. The copy differs from the true positions by whole lattice vectors, so the derivatives are
  unchanged and go to the same atom rows with the same sign (checked per term in the test). Without a cell the kernels read
  `m_geometry` as before; one branch per term, outside the arithmetic.
- **Cutoff check**: with a cell only the repulsion evaluation radius 20 Bohr is listed (argument of A3b: the list stores
  pairs, the kernel images at every call, the trigger is wrap-safe), not 20 + skin. Smallest accepted cube with defaults
  stays HB/XB-bound, L > 42.4 Bohr; with `hbond false` L > 40 Bohr.
- **Torn input**: the terms are correct for a molecule torn across a face, but the topology is still built from plain
  distances, so such a geometry can only be handed to an instance initialised intact (`UpdateMolecule`, constant topology).
- **Tests** (`ctest -R pbc_bonded_repulsion`, hbond off, 44 Bohr cube; water dimer, caffeine rotated and distorted by up to
  0.05 A so the inversions are non-zero): shifted by (L,-L,2L), centred on the corner, per-atom wrapped (torn) and one atom
  by L vs the cell centre: every term <= 2.7e-15 Eh, term gradients per atom <= 5.8e-15 Eh/Bohr (without a cell the same
  geometries differ by up to 2.0 Eh). 40 random waters (120 atoms, L 41 Bohr): E_nbrep = 27-image brute force 3.5e-17 for
  the all-pairs list and the periodic cell list (build radius 22 > L/2), images 7.1e-4 Eh, threads 1/4 7.6e-17; L 39.8
  refused. 2x1x1 supercell: water terms 1.4e-14, caffeine 3.3e-16. Central FD (torn caffeine, 8 atoms 0.08-0.32 A from a
  face) at h 2e-4 -> 1e-4 A: angle 1.2e-8 -> 2.9e-9, torsion + inversion 4.0e-9 -> 1.0e-9, repulsion 3.6e-8 -> 9.1e-9 Eh/A
  (max gradients 2.6e-2 / 9.0e-3 / 1.1e-1). Repulsion list: wrap by +L no rebuild, dE 0; 1 A move rebuilds, = fresh instance.
- **Without a cell**: -sp energy (16 digits) and `-dump_gradient` md5 identical before/after on the 7 structures of A2 at 1
  and 4 threads; peo201-water500 -sp wall time (median of 3) 0.350 -> 0.350 s at 4 threads, 0.690 -> 0.691 s at 1 thread.
- **Not tested**: sTors and extra torsions (no library structure has them under a cell; code identical to the torsion path),
  bonds with dynamic r0 through a torn bond only via the term energies, MD with a cell, GPU (unchanged; the A3b guard runs a
  cell on the CPU). Human production testing pending.

## A4b periodic hydrogen and halogen bonds (2026-11-01, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `196f3a50`), CPU only. **Periodic now additionally: HB and XB terms, their detection, the HB
coordination number (modulation of the X-H bond term) and the re-detection trigger.** Still plain distances: Coulomb/EEQ,
topology generation, reactive topology, BATM/ATM, and the charge ensemble `fragEnsembleBlend` (see below).

- **Kernels** (`calcHydrogenBonds`, `calcHalogenBonds`): each stored triple is evaluated in a local image frame anchored on H
  (on X for XB): A and B are the minimum images next to H, the neighbours of B and the acceptor parent C next to B, the
  neighbours D of C next to C (case 3), so distances, out-of-line terms, the lone-pair vector (case 4), the H...B=C angle and
  the D-B-C-H torsions come from one consistent minimum-image geometry; gradient rows unchanged (lattice-vector argument of A4a).
- **Detection** (`detectHydrogenBondsNative`, `detectHalogenBondsNative`): periodic `SpatialCellList` in the cell-list branch,
  minimum-image distances in both legacy branches, in the case-1 pre-test (cached hydrogens and j as images next to i), in
  `estimateHBStrengthCase1` / `estimateXBStrength` (kernel frame) and in the carbene angle. Lists store atom triples, never images.
- **HB coordination number** (`computeHBCoordinationNumbers`): H...B vector as minimum image. **Trigger** (`shouldUpdateHBXB`):
  per-atom displacement as minimum image, so a wrap by a lattice vector is no motion (valid while no atom moves > L/2).
- **Which radii are checked against L_min/2** (comment in `activeCutoffsBohr`): the search radii sqrt(hbthr1) 15.8 and
  sqrt(hbthr2) 21.2 (the detection must see the one image within them). With r_AB^2 + r_AH^2 + r_BH^2 < hbthr2 the three
  minimum-image vectors of a case-1 triple close a triangle (sum < sqrt(3 hbthr2) = 36.7 Bohr < L_min), case 2-4 likewise
  (r_HB <= r_AB + r_AH), so the accepted triple is the one the kernel evaluates. New: the XB kernel radius r_XB 20 Bohr (image
  of B within it unique). Not checked: the HB kernel radius r_cut 50 (acts on r_AB inside the frame, selects no image); a
  triple whose r_HB grows past L/2 between re-detections switches to the other image of B (in the water box the largest
  r_HB of any triple is 17.5 Bohr). HB CN: erf term exactly 0 beyond ~2 rcov.
- **Smallest accepted cube**: unchanged, defaults L > 42.43 Bohr (sqrt(hbthr2)); `hbond false` L > 40 Bohr (repulsion 20).
- **Tests** (`ctest -R pbc_hbond`; library structures hydronium-water2 (+1), acetic-acid-dimer.cyclic-hb, fi-pyridine, new
  registry names): (a) 44 Bohr cube, one instance initialised in the centre, `hb_update_force_every 1` (lists re-detected on each
  variant): shifted by (L,-L,2L), each fragment by a different lattice vector, centred on the corner, torn per atom: E_HB,
  E_XB, E_bond <= 1.3e-15 Eh, HB/XB/bond gradients <= 3.3e-15 Eh/Bohr, triple lists identical (21 / 48 / 550 HB, 3 XB; cases
  1/2/3 present, case 4 in none); without a cell the fragment-shifted geometry loses E_HB entirely (7.7e-3 Eh). (c) 40 random
  waters, L 43 Bohr, `hb_min_pair_energy_eh 0`: 1371 triples = independent 27-image enumeration, E_HB = reference 1.7e-18 Eh,
  656 triples through images (-7.1e-4 Eh); 1/4 threads 4e-18; cell-list = legacy detection; L 42.3 refused, with `hbond false`
  accepted. (d) 2x1x1 supercell: E_HB, E_bond = 2 x cell 1.4e-14, 2742 = 2 x 1371 triples. (e) torn composite (acetic acid
  dimer + F-I...pyridine, 10 atoms within 0.01 A of a face): HB FD 3.4e-11 -> 8.5e-12 Eh/A at h 2e-4 -> 1e-4 A (max |g| 4.6e-4).
  (f) wrap by +L: no re-detection, dE 0; RMSD 0.33 Bohr move re-detects, = instance re-detecting every call.
- **Pre-existing XB gradient defect (not periodic, not fixed)**: in `calcHalogenBonds` the short-range damping derivative has the
  wrong sign and the long-range one uses r/c instead of 1/r (`cuda/gfnff_kernels.cu` the same); FD residual 1.2e-4 of 5.4e-4
  Eh/A with and without a cell. Both corrections together give 4e-12 (tried, reverted: outside A4b). The test checks only that
  the cell changes nothing (residual torn = residual intact without cell, 4.8e-21). Energy unaffected.
- **Charge ensemble not periodic**: hydronium-water2 torn per atom changed E_HB by 4.4e-3 and E_bond by 0.32 Eh with the
  default `frag_charge_model ensemble` (net charge, several fragments): `fragEnsembleBlend` re-fragments with plain distances.
  The test uses `frag_charge_model reference`; open item (TODO.md).
- **XB list near L/2**: a system wider than L/2 sees its own image within the 21.2 Bohr XB search radius (composite: one extra
  triple F...O image, r_XB > 20, energy 0). Correct for a periodic system, but a cell-vs-no-cell list comparison differs.
- **Without a cell**: -sp energy (16 digits) and `-dump_gradient` md5 identical before/after on the 7 structures of A2 at 1 and
  4 threads; peo201-water500 wall time (median of 3) 0.347 -> 0.346 s at 4 threads, 0.691 -> 0.692 s at 1 thread.
- **Not tested**: HB case 4, XB in a halogen-bonded crystal (only fi-pyridine exists, its N is bonded to I in the topology, so
  its XB triples go to other acceptors), MD with a cell, `nb_cell_list_min_atoms` XB cell-list branch with a cell (HB one
  tested), GPU (unchanged; the A3b guard runs a cell on the CPU, `prepareHBDetectInputs` still plain). Human production testing pending.
