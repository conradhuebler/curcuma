# WP: Periodic non-bonded interactions (minimum image)

**Status (2026-11-03, after A6): native CPU GFN-FF with an orthorhombic unit cell is periodic** (MIC in CN, D4,
repulsion, bonded, HB/XB, topology; Ewald Coulomb and EEQ with `coulomb_ewald`, A1-B3, P1). With such a cell
`wall_potential = pbc` in a box is the pure cell wrap, without clearance test or reflection (A6). The stopgap
described below **stays** for everything that has no periodic energy: no cell, the spherical container, every other
method (gfn1/gfn2, external xtb/tblite/Ulysses, UFF, QMDFF, CG) and GFN-FF with `coulomb_ewald false`; for those a
cell is not wrapped and the COM removal is the vacuum one (A6b, matrix there). Open: GPU paths with a cell (CPU
fallback), triclinic cells, NPT.

## What exists today

`SimpleMD::wrapIntoContainer()` moves a whole fragment back into the container when
its centre of mass leaves it: per axis for a box, to the antipodal point for a
sphere. It exerts no force, so unlike the harmonic and logfermi walls it does not
heat the system and does not do work on an escaping molecule. Since A6 it is skipped when the method reports a
periodic energy and a cell is set (box container only); then `applyPeriodicBoundaryConditions()` alone wraps.

## Why it needs a stopgap (still true without a periodic energy)

The text of this section describes the state before A2-B3 and still describes every case without a periodic
energy (list above). The GFN-FF energy without a cell has **no minimum-image convention**. `FFWorkspace`, the CN
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
   of `wrapIntoContainer()` and it reduces to the wrap itself. Done in A6 for the periodic case only (the container
   wrap is skipped and the cell wrap remains); the stopgap stays for all methods without a periodic energy.

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

## A4c periodic topology generation and reactive topology (2026-11-01, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `af776bf8`), CPU only. **Periodic now additionally: the bond criterion, the whole topology build
(hybridisation, neighbour sets, rings, pi systems, angle/torsion/inversion/sTors generation, HB/XB set-up detection), the
reactive bond scan, the charge-ensemble contact window (`fragEnsembleBlend`) and the topology cache key.** A structure with a
cell whose molecules are torn in the input (per-atom wrap, SimpleMD fragment wrap) is initialised directly. Still plain
distances: Coulomb and the per-step Phase-2 EEQ (Ewald is block B), ATM/BATM, native GFN1/GFN2; GPU unchanged (A3b guard).

- **Bond criterion** (`getCachedBondList`): with a cell the pair distance is the minimum image; the distance functor is
  chosen once per call (generic lambda), the plain instantiation is the old loop. The same for the reactive scan
  (`detectReactiveBondChanges`, form/break thresholds) and the contact window `fragWindowEdges`; the bond-graph split and
  the window gradient of `fragEnsembleBlend` use `pairVectorBohr()` (one branch per pair, only in the ensemble path).
- **Reassembled geometry instead of hundreds of MIC sites** (`unwrapFragmentsForTopology`, `TopologyGeometryScope` at
  the start of `generateGFNFFParameterSet`): every connected component of the minimum-image bond list is placed as one
  intact molecule (lowest-index atom keeps its position, neighbours at the minimum image next to the atom they were
  reached from, breadth first); `m_geometry_bohr`/`m_geometry` are swapped for the build and restored afterwards (also on
  an exception). The copy differs from the input by whole lattice vectors per atom, so every minimum-image quantity (CN,
  D4, repulsion, HB/XB detection, bonded terms) is unchanged, and the plain-distance topology code (`bangl` in the
  hybridisation, `nnNearestNoM`, the icase-2/3 neighbour predicate, angle/torsion/inversion generation, `sdum` sorting,
  the Phase-1 EEQ matrix) sees intact molecules. A pair accepted by the icase-2/3 predicate is a bond (its threshold is
  the icase-1 one without `fm >= 1`), so its plain distance in the copy is the minimum image; a plain distance is never
  shorter than the minimum image, so no spurious contact appears. The bond list is detected on the input first (it
  defines the placement) and once more inside the build (same set: lattice shifts leave the minimum image unchanged;
  one extra O(N^2) pass per build, cell only). The reference geometries recorded during the build (pair-list skins, HB
  trigger, `m_last_eeq_geometry`) hold the reassembled copy; the skin and HB triggers take minimum-image displacements
  (A3b/A4b) and see the restored input as unmoved, the EEQ reference differs and costs one Phase-2 solve on the first call.
- **Refused**: a component whose bonds close around the cell (a chain bonded to its own periodic image, mismatch > 1e-6
  Bohr on a second path to an atom): `std::runtime_error` naming the bond. Such a system has no intact placement; the
  bonded terms alone would be right (A4a), the topology around the closing bond not. **Reported, not repaired**: a bond
  that only the charge-shrunk pass 2 of the q-loop adds between two separately placed fragments and that spans the cell
  (warning per bond; its angles/torsions are built from the plain vector). Not covered: `topology_mode auto` re-derives
  the topology from `Calculation()` on the input geometry as is (the scope is only in `generateGFNFFParameterSet`).
- **Dense N x N distance matrix removed** (`calculateTopologyInfoOnce` phase 2.1): its only reader was `nnNearestNoM`
  (bonded pairs), which forms the distance itself with the same expression (bitwise). 240 MB at 5475 atoms; on
  peo201-water500 (2910 atoms, 4 threads) topology 105.8 -> 90.5 ms, parameter generation 330 -> 281 ms, max RSS
  559 -> 494 MB, `-sp` wall time (median of 3) 0.348 -> 0.329 s at 4 threads, 0.691 -> 0.667 s at 1 thread.
  With a 100 A cube: topology 95.5 ms, parameter generation 296 ms, RSS 501 MB, wall 0.346 s.
- **Phase-1 EEQ / topology charges** (corrected in B3): the Phase-1 matrix uses topological (bond-graph) distances
  only and is zero between fragments (`calculateTopologyChargesMultiRHS`), so it reads no geometry across fragments
  and needs no cell; the earlier statement here ("plain Coulomb over the reassembled geometry, non-periodic between
  fragments") was wrong. Until B3, Coulomb and the per-step Phase-2 EEQ (`computeSharedDistances`, `prepareCNAndEEQ`)
  read the input geometry with plain distances: a torn input gave a different Coulomb energy than the intact one
  (caffeine torn -0.1622 vs intact -0.1744 Eh). Since B3 both are periodic by Ewald summation (`coulomb_ewald`).
- **Topology cache** (`computeTopologyFingerprint`): `|cell=Lx,Ly,Lz` (Bohr, 12 digits) is appended only when a cell is
  active, so the fingerprint and every existing cache file without a cell are unchanged (checked: caffeine, urea197
  fingerprints and charges before/after identical); a cache written with a cell is not read without one, nor with another
  cell, and vice versa ("fingerprint mismatch - recalculating").
- **Reactive topology**: minimum-image scan, scan trigger `react_check_disp_bohr` on minimum-image displacements (a
  fragment wrap triggers no scan); a rebuild runs the full `generateGFNFFParameterSet` with the reassembled geometry, so
  two molecules that bond across a face become one intact component. The displacement test of the react scan and the
  cutoff list are unchanged otherwise.
- **Smallest accepted cell**: the topology thresholds (bond criterion < 6 Bohr, reactive scan 2.6 x rcov sum, topology-init
  CN 40 Bohr with minimum image, A2) add no blocker; defaults L > 42.43 Bohr (sqrt(hbthr2)), `hbond false` L > 40 Bohr
  (repulsion 20; water dimer 40.5 accepted, 39.9 refused).
- **Tests** (`ctest -R pbc_topology`, labels pbc/gfnff/react): (a) caffeine, hydronium-water2 (+1, default charge
  ensemble), acetic acid dimer, 44 Bohr cube, torn per atom at the corner and initialised directly: bond list, angle/
  dihedral/inversion counts, hybridisation, fragment count identical to the intact instance, topology charges identical
  (anchored and fully torn variants), every periodic term <= 1.5e-15 Eh (hydronium E_HB -8.31e-3 Eh equal under the
  ensemble model); with cell == without cell for the intact molecule; torn without a cell has fewer bonds (not vacuous).
  (b) 2x1x1 supercell: counts exactly doubled, periodic terms = 2 x cell to 7.1e-15 Eh (Coulomb not). (c) 60 random
  waters, L 43.5 Bohr, wrapped per atom (8 torn): 120 bonds = the O-H pairs = an independent 27-image enumeration of the
  charge-shrunk criterion with the instance's Phase-1 charges; 110 without a cell; threads 1/4 identical. (d) react,
  H + H from 4.0 to 1.2 Bohr in the centre and across the x face: formation at the same step (20), E_after, E_final and
  dE_jump (-0.1471063 Eh) equal to 1e-12. (e) cache fingerprints: cell != plain, same cell reproduced, other cell
  different. (f) smallest cell as above.
- **Without a cell**: -sp energy (16 digits) and `-dump_gradient` md5 identical before/after on the 7 structures of A2 at 1
  and 4 threads (benzene, caffeine, methyl-glucopyranoside, water8, hydronium-water2, urea-water197, peo201-water500).
- **Not tested**: polymers wider than the cell (the refusal path only by reading), metals and eta fragments under a
  cell, `nb_cell_list_min_atoms` >= 800-atom neighbour-set path with torn input (peo201 with a 100 A cube runs, numbers
  not compared), charge groups / ionpairs with a cell, MD with a cell and react, `topology_mode auto`, GPU. Human
  production testing pending.

## B0 stand-alone Ewald sum for point charges (2026-11-01, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

`src/core/energy_calculators/ff_methods/ewald.h/.cpp` (namespace `curcuma::ewald`, class `EwaldSum`): direct Ewald
sum (no PME) in an orthorhombic `PBCUtils::OrthoCell`, atomic units (Bohr, e, Eh, E = q_i q_j / r). **Not connected to
GFN-FF or the EEQ**; nothing calls it outside `test_ewald`, so the default path is untouched.

- **Terms**: real space `sum_{i<j} q_i q_j erfc(a r)/r` (minimum image, r < r_cut <= L_min/2, linked cells when every
  axis has >= 3 cells, else all pairs); reciprocal `(1/(2V)) sum_{k!=0} (4 pi/k^2) exp(-k^2/(4a^2)) |S(k)|^2` over the
  half space (n and -n, weight doubled), phase factors per atom and axis by recursion; self `-(a/sqrt pi) sum q^2`;
  background `-pi Q^2/(2 a^2 V)`. Boundary condition tin-foil: no dipole term, k = 0 dropped.
- **Interface**: `energy`, `gradient` (dE/dr, Eh/Bohr), `potentials` (phi_i = dE/dq_i = (A q)_i without A),
  `evaluate` (all in one pass, with the four terms), `matrix` (dense A, tests only). `realSpaceKernel()` is the hook for
  B2 (Gaussian-charge correction -erfc(g_ij r)/r), not built.
- **Parameters** (`chooseEwaldParameters(cell, eps, r_cut)`, r_cut default L_min/2): `a = sqrt(-ln eps)/r_cut` from
  exp(-(a r_cut)^2) = eps, `k_max = 2 a sqrt(-ln eps)` from exp(-k_max^2/(4a^2)) = eps. A per-term criterion, not a bound
  on the total error: at eps = 1e-6 (N = 1000 random charges) the energy is off by 2.6e-6 relative, the gradient by
  8.5e-9 of its maximum.
- **Literature**: Ewald (1921), Allen & Tildesley (2017), Frenkel & Smit, de Leeuw, Perram & Smith (1980) cited from
  memory, not checked against the pages; the formulas were re-derived and are checked by the tests below.

Tests (`ctest -R ewald_`, `test_cases/unit/test_ewald.cpp`; `ewald_correctness` 0.4 s, `ewald_performance` 6.4 s, label `perf`):

| Test | Result |
|---|---|
| (a) NaCl, 8-ion cubic cell, Madelung M | Ewald 1.747564594633 (3 alphas, spread 2e-13); Evjen (neutral cubes, n = 40/80, Richardson) 1.7475645946, rel 8e-14 (independent of Ewald); literature 1.747565 (one web search: Wikipedia, Spektrum lexicon; 1.747564594633 quoted from Borwein et al. in the search summary, page not opened) |
| (b) alpha independence, 64 neutral random charges, L 30 Bohr, eps 1e-12, r_cut 15/12/9 (9 = linked-cell path) | E 6.3e-13 rel, gradient 5.2e-13, potentials 1.1e-13 (of max) |
| (c) naive long-double Ewald (all images, full k cube, no recursion; self-consistency, not independent), 12 charges, 9 x 10 x 11 Bohr | E rel error 4.6e-7 / 4.4e-9 / 5.1e-11 / 4.9e-13 / 1.6e-13 at eps 1e-6 / 1e-8 / 1e-10 / 1e-12 / 1e-14 |
| (d) FD, 16 charges, charged cell | gradient 1.7e-10 Eh/Bohr at h = 1e-4, 6.7e-10 at 2e-4 (ratio 4.04, quadratic); phi vs dE/dq 1.9e-11 rel; E = 1/2 q.phi 2e-14 |
| (e) translation by an arbitrary vector and by L; single atoms shifted by L | 6e-15 and 4e-16 (E and gradient) |
| (f) 2x1x1 supercell; dense A | E/N 1.7e-15, gradient, potentials 1e-15; A q = phi 5e-14, A symmetric exactly |
| (g) one charge in a uniform background, cubic L = 10 Bohr, E = -xi q^2/(2L) | Ewald xi = 2.837297479 (3 alphas, spread 1.1e-12); independent neutral-cube sum (charge + background cube, Bethe shift +pi/(6L), m = 40/80, Richardson) 2.8372974957, rel 5.7e-9; the value 2.837297 is from the task text, not found by the search |

(h) Timing, eps 1e-6, energy + gradient + potentials, one thread, uniform random positions at the atom density of
water (0.1003 / A^3), scan of 8 r_cut values from 8 Bohr to L/2:

| N | L (Bohr) | fastest r_cut (Bohr) | alpha | N_k (half space) | time (s) |
|---|---|---|---|---|---|
| 1000 | 40.7 | 15.1 | 0.247 | 3537 | 0.016 |
| 3000 | 58.7 | 17.1 | 0.217 | 7121 | 0.086 |
| 6000 | 73.9 | 20.4 | 0.182 | 8439 | 0.250 |

Scaling exponent 1.54 (1000 to 6000 and 3000 to 6000), as expected for O(N^1.5) with optimal alpha. The k loop in
`EwaldSum::evaluate` (N_k x N complex products) is the place to parallelise (over k with per-thread S, or over atoms).

- **Virial (derived here, not implemented, not tested)**: diagonal strain of an orthorhombic box, W_aa = -dE/d eps_aa:
  real `-sum_pairs q_i q_j K'(r) d_a^2 / r`; reciprocal per k `E_k [1 - 2 k_a^2 (1/k^2 + 1/(4a^2))]`; self 0;
  background `E_bg`. Needed for block C (NPT).
- **Not implemented**: PME/particle mesh, triclinic cells, virial/pressure, vacuum (dipole) boundary term, slab
  correction, threading, cached geometry between calls.
- **Not tested**: very small cells with r_cut far below the nearest pair distances, a few close pairs (r << 1 Bohr),
  more than 6000 charges, eps below 1e-14. Human production testing pending.

## B2 periodic GFN-FF Coulomb pair term by Ewald summation (2026-11-01, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

Switch `-gfnff.coulomb_ewald` (default **false**; cell only). The EEQ charges are still the non-periodic ones (B3
pending): with the switch on, the energy is not variational in the charges, and the forces are the fixed-charge
derivatives (the finite-difference check below is therefore at fixed charges). B3 sets the default to true once the
EEQ matrix is periodic. Without a cell, or with the switch off, the path is bitwise what it was (test (g)).

### Step 0: what the Coulomb code actually does (checked against the source before anything was changed)

- (i) Pair term `q_i q_j erf(gamma_ij r)/r`, `gamma_ij = 1/sqrt(alp_i + alp_j)` with `alp_i = alpeeq_i =
  (alpha_i + ff q_i^(1))^2` (charge-corrected, Phase-1 charge): `ff_workspace_gfnff.cpp:1090-1093` (implicit loop),
  `:1132-1146` (stored list), `gfnff_method.cpp:12572` (`c.gamma_ij`), `eeq_solver.cpp:1231` (alpeeq). The hypothesis
  "gamma = 1/sqrt(alpha_i^2 + alpha_j^2)" holds with `alpha_i^2 = alp_i`. Far-field shortcut: `erf = 1` for
  `gamma r >= 6.5` (`:1035`), exact in FP64.
- (ii) Local terms, unchanged by B2 (`ff_workspace.cpp:549-552`, TERM 2+3 in `postProcess`): `-q_i chi_i(CN)` and
  `1/2 q_i^2 (gam_i + sqrt(2/pi)/sqrt(alp_i))`. The second part equals the Gaussian self-energy
  `q^2/sqrt(2 pi alp)` = `1/2 q^2 lim erf(gamma_ii r)/r`, so the GFN-FF energy is "Gaussian clouds minus their
  self-interaction", which is what the periodic sum has to reproduce.
- (iii) No pair exceptions: every `i < j` pair, bonded or not, same or different fragment
  (`gfnff_method.cpp:12635-12639`, implicit loop `ff_workspace_gfnff.cpp:1060-1104`). Only cutoffs exist:
  `coulomb_r_cut` (default none) and `eeq_distance_cutoff` (explicit list only); pairs with a NaN charge or
  `alp <= 0` are skipped.
- (iv) `alpha_eeq` (`gfnff_par.h:99`, Bohr): H 0.585, C 0.903, N 1.278, O 0.905, Cl 2.057, Br 2.002, I 2.178,
  Ne 2.942 (largest), He 0.432 (smallest); `alp = alpha^2` from 0.19 (He) to 8.65 (Ne), with the charge correction
  `ff` up to 0.5 (group 7). `gamma_min = 1/sqrt(2 max alp)`: 0.34/Bohr for Cl-Cl, 0.24 for Ne-Ne. The correction
  `erfc(gamma r)/r` is below 2e-17 (one ulp of the 1/r it corrects, `erfc(6)`) beyond `r_g = 6 sqrt(2 max alp)`:
  7.7 Bohr for a water box (O with charge correction: alp 0.84), 17.4 Bohr with Cl, 25 Bohr with Ne.
  `GAUSSIAN_CORRECTION_XMAX = 6` in `ewald.h`; erfc(5.04) would be 1e-12.
- (v) dE/dq does not exist in the workspace. The charges enter through the EEQ matrix
  `A_ij = erf(gamma_ij r)/r` (`eeq_solver.cpp:1475-1483`, `:3391`), the gradient treats them as fixed
  (Hellmann-Feynman by the EEQ stationarity, TERM 1b for chi(CN)); `NumGradFixedCharges` (`gfnff_method.cpp:14105`)
  runs `m_workspace->calculate(false)` at the stored charges. B2 adds `FFWorkspace::coulombPotentials()`.

### Derivation (Gaussian correction and the background of a charged cell)

Target: `E_target = 1/2 sum'_{(i,j,n)} q_i q_j erf(g_ij |r_ij + n|)/|r_ij + n|` (prime: no `(i,i,0)`), with the
tin-foil / "k = 0 dropped" convention of B0. Let `F_G = (1/2V) sum_{k != 0} (4 pi/k^2) |rho_G(k)|^2` with
`rho_G(k) = sum_i q_i exp(-k^2 alp_i/4) e^{i k r_i}` (cloud `q (pi alp)^{-3/2} exp(-r^2/alp)`, pair potential
`erf(r/sqrt(alp_i + alp_j))/r`), which contains the Gaussian self-energies `S_G = sum_i q_i^2/sqrt(2 pi alp_i)`;
`E_target = F_G - S_G`. With the point-charge structure factor `S(k)`:
`F_G - F_pt = sum_{k != 0} (4 pi/(2V k^2)) (|rho_G|^2 - |S|^2)`. The absolutely convergent lattice sum of
`h_ij(r) = -erfc(g_ij r)/r` over all `(i,j,n)` equals, by Poisson summation, the same k sum *including* k = 0, where
`h_ij(k = 0) = -pi/g_ij^2 = -pi (alp_i + alp_j)`, so the k = 0 part is `-(pi Q/V) sum_i q_i alp_i`. Its `(i,i,0)`
term is `S_G - (point self)`, and B0's `E_pt = F_pt - (point self)`. Hence

    E_target = E_pt(real + recip + self + background)
               - sum'_{(i,j,n)} 1/2 q_i q_j erfc(g_ij |r_ij + n|)/|r_ij + n|
               + (pi Q / V) sum_i q_i alp_i

The last term is the interaction of the neutralising background with the difference between a cloud and a point
(space average of `-q erfc(r/sqrt alp)/r` is `-pi q alp`); it vanishes for Q = 0, is alpha_E-independent, and is
therefore invisible to the alpha-independence test. It is checked against an independent pure reciprocal-space sum
of the clouds (test e2). Limit `L -> inf`, `r_g << L`: `E_pt -> sum_{i<j} q_i q_j/r`, the correction
`-> -sum_{i<j} q_i q_j erfc(g r)/r`, sum `= sum_{i<j} q_i q_j erf(g r)/r`; the images `(i,i,n != 0)` carry
`erfc(g L) < 1e-16` and are not summed (enforced: `max(ewald_r_cut, r_g) < L_min/2`, also in the A1 check).

Implementation: one real-space kernel `[erfc(a r) - erfc(g_ij r)]/r` per minimum-image pair within
`r_cut = max(ewald_r_cut, r_g)` (`EwaldSum::realSpaceRange`, correction only for `g r < 6`), reciprocal part, self,
background and Gaussian background from the B0 code (`reciprocalRange`, `selfAndBackground`; `evaluate()` is now
the composition of these, B0 tests unchanged, 24/24). In `FFWorkspace`: `calcCoulombEwald(p)` takes the linear
slice `p` of the per-step pair list (linked cells) and of the half-space k list, every partition accumulates its own
energy, gradient and potentials (`coulomb_potentials`), `reduce()` sums in partition order, `postProcess()` adds the
O(N) terms; the same `CxxThreadPool` pattern as the other kernels, reproducible at a fixed thread count. Parameters:
`ewald_tolerance` (1e-8), `ewald_r_cut` (0 = `min(L/2, max(15 sqrt(L/40), r_g))`, from the B0 timing optimum
~ sqrt(L)), `ewald_alpha` (0 = from the tolerance; set for the alpha tests). Refused: `eeq_distance_cutoff > 0`,
`coulomb_r_cut > 0`, `r_cut >= L_min/2`.

### Tests (`ctest -R pbc_coulomb_ewald`, `test_cases/unit/test_pbc_coulomb_ewald.cpp`, 41 checks, 10 s)

| Test | Result |
|---|---|
| (a) alpha independence, 20 H2O in L = 44 Bohr, r_cut 19, alpha 0.25/0.32/0.42 (N_k 10k/21k/47k), tol 1e-10 | E_coul 1.3e-12 rel, gradient 1.1e-10, potentials 4.1e-12 (of max); EEQ charges bitwise equal across alphas |
| (b) large-cell limit, periodic E_coul (fixed charges) minus the plain path, L = 150/200/300 Bohr | H2O (dipole 0.692 e Bohr): -2.97e-7 / -1.25e-7 / -3.72e-8 Eh = exactly -2 pi M^2/(3V), exponent 3.00; water dimer (registry geometry has no net dipole) 7.4e-9 / 1.7e-9 / 2.3e-10, exponent 5.00; benzene 3.3e-10 / 7.8e-11 / 1.0e-11, exponent 5.00 |
| (c) fixed charges and fixed alp/gam/chi (the non-periodic Phase-1 EEQ gives a torn or replicated box slightly different alpeeq): translation, atom-wise wrap, 2x1x1 supercell | energy 0 / 1.3e-16 / 3.1e-15 rel, Coulomb gradient 6.6e-15 / 6.4e-15 / 1.4e-14; 4 threads vs 1: E 7.6e-16, gradient 6.6e-14, potentials 5.1e-15; 4 threads run-to-run bitwise |
| (d) FD at fixed charges, h = 1e-4 / 2e-4 A, dispersion_tail off | Coulomb component vs `gradientCoulomb`: 1.36e-9 / 5.42e-9 Eh/Bohr (ratio 4.00, quadratic; 1.7e-7 of max |g|); potentials vs FD dE_pair/dq 4.8e-12; E_pair = 1/2 q.phi 1.3e-16. The total vs `NumGradFixedCharges` differs by 2.2e-6 Eh/Bohr independent of h, the same with the plain Coulomb in the cell and 1.5e-6 without a cell (`--fd-diagnose`): pre-existing, not Coulomb, not located (TODO.md) |
| (e) charged cell, H3O+ 2 H2O + 18 H2O (Q = +1), L = 44, 3 alphas | E_coul 1.1e-11 rel, gradient 7.4e-11, potentials 2.4e-12; background -2.95e-4 .. -1.05e-4 (alpha dependent), Gaussian background -2.29e-4 (alpha independent) |
| (e2) `EwaldSum` with Gaussian widths vs an independent pure k-space sum of the clouds (long double, 121^3 k, 10 charges, Q = 1.39, cell 12 x 13 x 14, alp 0.08-0.34 Bohr^2 so that r_g = 4.9 < L/2) | energy 1.5e-13 rel, potentials 1.3e-13; without the Gaussian background the error is 8.2e-4 (the test discriminates it); 3 r_cut spread 5e-14; gradient vs FD 3.5e-9; dense Gaussian A q = phi 7.3e-14 |
| (f) alp = 1e-12 (point limit) vs B0 point-charge Ewald | energy 3.7e-15, gradient 0, potentials 1.2e-14 |
| (g) no cell: `coulomb_ewald true` vs false; CLI `-sp` + md5 of `-dump_gradient` on benzene, caffeine, methyl-a-d-glucopyranoside, water8-cluster, hydronium-water2 (+1), urea-water-cluster197, peo201-water500 at 1 and 4 threads, binary before vs after | bitwise; all 14 md5 identical (energies -2.36477550, -4.67273706, -3.96002715, -2.61818827, 0.29809137, -28.60465251, -367.39771445 Eh) |
| (i) `ctest -R 'gfnff|frag|eeq|simplemd|pbc|ewald|react'` list before/after | identical plus `pbc_coulomb_ewald` (116 -> 117) |

(h) Timing (`test_pbc_coulomb_ewald --perf <xyz> <L>` and `--perf-box <nwat> <L>`, tol 1e-8, energy + gradient,
median of 3 `calculate(true)`; "Coulomb" = CPU-sum of the Coulomb term over partitions incl. pair-list/phase setup,
"wall" = the whole workspace step):

| System | Path | 1 thread Coulomb / wall (ms) | 4 threads Coulomb CPU-sum / wall (ms) |
|---|---|---|---|
| peo201-water500 (2910 atoms, droplet) in L = 110 Bohr | plain O(N^2), no cell (4.2 M pairs, far-field shortcut) | 18.6 / 54 | 19.1 / 15.4 |
| | Ewald auto r_cut 24.9 (N_k 36k, 547k pairs) | 233 / 276 | 241 / 79 |
| | Ewald r_cut 40 (N_k 8.8k, 1.6 M pairs) | 152 / 197 | 151 / 60 |
| random water box, 970 H2O, L = 58 Bohr (water density) | plain, no cell | 18.8 / 86 | 19.4 / 28 |
| | Ewald auto r_cut 18.1 (N_k 14k, 534k pairs) | 119 / 204 | 122 / 66 |
| | Ewald r_cut 20 (N_k 10k, 726k pairs) | 109 / 194 | 112 / 63 |

The Ewald Coulomb costs 6-12 x the all-pairs path (which is 4.4 ns per pair thanks to the far-field shortcut; the
k loop is N_k x N complex products, no SIMD, no PME); the workspace step with Ewald scales 3.1 x from 1 to 4 threads.
The E_coul of the droplet in the cell (-53.2373) differs from the free droplet (-53.2347) by the images; the water
box shows -82.501 (periodic) against -82.737 (the same coordinates without images).

- **Not implemented**: periodic EEQ (B3: the matrix-vector product is `coulombPotentials()`, the dense matrix
  `EwaldSum::matrix(pos, &alp)`), GPU (cell -> CPU fallback), PME, virial, triclinic cells, a Verlet skin for the
  Ewald pair list (rebuilt every step, 3-11 ms at 2910 atoms).
- **Not tested**: MD energy conservation with Ewald, systems with Ne/heavy halogens near the `r_g < L/2` limit,
  more than 2910 atoms, `topology_mode react` with Ewald, the ALPB solvation combined with a cell.
- Human production testing pending.

## B3 periodic EEQ with the Ewald matrix (2026-11-01, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `31da9f67`), CPU only. **With a cell, `coulomb_ewald` is now on by default**: the Phase-2 EEQ
charges are solved with the Ewald matrix of the Gaussian charges, the Coulomb energy is the B2 Ewald sum with the same
parameters, and the total gradient is consistent with the energy (finite-difference gate below). `coulomb_ewald false`
is the emergency exit (plain distances in Coulomb and EEQ, not periodic). Without a cell nothing changes (bitwise, (g)).

### Step 0: what the EEQ code does (checked against the source before the change)

- **Matrix build per solver path** (`eeq_solver.cpp`): every path starts from the same dense augmented matrix
  `m_phase2_A` ((N+nfrag)^2, row-major): `calculateFinalCharges` fills the N x N block in 64 x 64 tiles of the lower
  triangle from the packed distance buffer (`m_shared_srab`, plain distances, `GFNFF::computeSharedDistances`),
  `A_ij = erf(r/sqrt(alp_i + alp_j))/r` (far field `1/r` for `r^2 >= 42.25 (alp_i + alp_j)`), diagonal
  `gam_i + sqrt(2/pi)/sqrt(alp_i)`, then the constraint rows/columns (0/1 per fragment; left unwritten for ppcg,
  `fillDeferredConstraintBlock`) and `dispatchSolve`: `cholesky`/`ldlt` = Schur complement over A_nn with the
  constraints as right-hand sides, `lu` = the augmented system, `ppcg` = projected CG on A_nn with the atom -> fragment
  map and the `localinv` preconditioner built from the entries of A_nn within a radius, `pcg` = nfrag+1 block-Jacobi
  solves, `batched` = per-fragment LU without inter-fragment coupling. `eeq_matrix_rebuild_eps_bohr` caches the N x N
  block, `eeq_refactor_eps_bohr` the Cholesky factor. `alp_i` = `alpeeq_i = (alpha_i + ff q_i^(1))^2` from the
  topology (Phase-1 charge, fixed per topology).
- **chi(CN)**: Phase 2 right-hand side `x_i = -chi_i + dxi_i (+ amide-H -0.02) + cnf_i sqrt(CN_i)` with the
  per-step CN (`calculateFinalCharges`, "4. Setup RHS"); the workspace self-energy uses the same
  `chi_base_i + cnf_i sqrt(CN_i)` (`postProcess`, TERM 2+3) and its derivative is TERM 1b
  (`qtmp_i = q_i cnf_i/(2 sqrt(CN_i))` into the CN chain rule `m_dcn.applyAdd`). `dxi` is topology-based; for
  2-coordinate carbon it reads a bond angle from the plain geometry (`detectHybridization`, carbene test).
- **Phase 1** (`calculateTopologyChargesMultiRHS`): Floyd-Warshall/Dijkstra topological distances over the bond
  graph, `A_ij = 0` between fragments, chi with the integer neighbour count. **It reads no geometry across
  fragments**, so the A4c statement "non-periodic between fragments" was wrong: with the periodic bond graph (A4c)
  Phase 1 needs no cell and is unchanged by B3. Its only geometry use is the dxi bond angle above.
- **Charge consumers in the gradient**: the Phase-2 charges `q` enter only the Coulomb energy (`calcCoulomb` /
  `calcCoulombEwald`, TERM 2+3, TERM 1b). HB/XB strengths, `zetac6` of D4, the bond `fqq` and `alpeeq` use the
  Phase-1 topology charges, constant per topology. Hence `dE/dr = dE/dr|_q + dE/dq . dq/dr` with `dE/dq = 0` on the
  constraint surface (EEQ stationarity): the fixed-charge Coulomb gradient (dA/dr at fixed q, B2) plus TERM 1b
  (dchi/dCN) is the complete gradient. **No bilinear term `z^T (dA/dr) q` and no polarisation identity are needed**,
  provided the charges minimise the same functional the workspace evaluates; the B2 gradient is exactly the gradient
  of the quadratic form `1/2 q^T A_ew(r) q` at fixed q (its FD test d). What remained inconsistent in B2 was only that
  the charges came from the non-periodic matrix.
- **Pre-existing**: the h-independent 1.5-2.2e-6 Eh/Bohr discrepancy of `NumGradFixedCharges` against the analytic
  gradient (B2 test d, TODO.md) is not a gradient error: `NumGradFixedCharges` evaluates `calculate(false)` with the
  CN of the plain 40-Bohr threshold formula and without refreshing the D4 C6 weights, i.e. a different function than
  the one the gradient belongs to. The full energy with re-solved charges (test (e) below) agrees to 2.3e-8 Eh/A.

### Formulas

Phase-2 functional with the pair sum over all images (`gfnff_engrad.F90` goed_gfnff, images added):

    E(q) = sum_i [ -chi_i(CN) q_i + 1/2 (gam_i + sqrt(2/pi)/sqrt(alp_i)) q_i^2 ] + 1/2 q^T A_ew q
    A_ew,ij = erfc(a r)/r - erfc(g_ij r)/r [MIC pairs, r < r_cut, g r < 6] + (4 pi/V) sum_{k!=0} e^{-k^2/4a^2}/k^2 cos(k.r_ij)
              - (2a/sqrt(pi)) delta_ij - pi/(a^2 V) + (pi/V)(alp_i + alp_j)

(ewald.h; `E_pair = 1/2 q^T A_ew q` holds for the B2 sum incl. the Gaussian correction and both backgrounds, test
e2/f of B2). The system `(A_ew + D) q + C^T lambda = chi`, `C q = q_frag` has the same chi, D and constraints as
without a cell. The Gaussian background of a charged cell, `(pi Q/V) sum_i q_i alp_i`, is kept as the symmetric
quadratic form `1/2 q^T [(pi/V)(alp_i + alp_j)] q` (identical for every q, not only on the constraint surface);
treating it as the linear shift `chi_i -> chi_i - pi Q alp_i/V` gives the same charges, because the two gradients
differ by `(pi/V)(sum_j q_j alp_j)` times the all-ones vector, which lies in the span of the constraint rows and is
absorbed by the multipliers. The diagonal of `A_ew` is `sum_k 2 w_k - 2a/sqrt(pi) - pi/(a^2 V) + 2 pi alp_i/V`
(`EwaldSum::diagonal`), the reciprocal part being the same for every atom (`|e^{ik r}|^2 = 1`): the Jacobi
preconditioner needs no matrix.

### Implementation

- `EEQSolver::setPeriodicEwald(cell, parameters, dense_max)` (from `GFNFF::configurePeriodicEwald()`, called in
  `InitialiseMolecule` after the topology-cache reset, before the topology build whose Phase-2 solve is the first
  periodic one); `calculateFinalCharges` branches to `calculateFinalChargesPeriodic` after chi/gam/alpeeq are known.
  Two paths: **dense** (`EwaldSum::matrix`, now k vectors in blocks with symmetric rank updates and the phase-table
  recursion, O(N^2 N_k) flops as GEMM; then `dispatchSolve` with the direct method fixed, cholesky / lu / ldlt) and
  **matrix-free pcg** (`PeriodicEEQOperator::apply` = `realSpaceRange` + `reciprocalRange` split over the pool,
  reduced in partition order, + `selfAndBackground` + D x; projected CG on the constraint null space with the
  constraint-preconditioner `z = M^-1 r - M^-1 C^T (C M^-1 C^T)^-1 C M^-1 r`, Jacobi M, warm start from the previous
  charges with the constraints re-enforced, stop `|P r| <= eeq_ppcg_tol (|b|+1)`, cap `eeq_ppcg_max_iter`; a
  breakdown or no convergence falls back to the dense solve for N <= 4 x dense_max, else error). `auto` = dense
  up to `eeq_ewald_dense_max_atoms` (1000), pcg above; `ppcg` -> pcg with one warning; `batched` refused;
  reaction field, `use_iterative_refinement`, `eeq_distance_cutoff`, `coulomb_r_cut` refused. The periodic path has
  no matrix cache and no cached factor. `GFNFF::Calculation` skips `computeSharedDistances()` (O(N^2) buffer) when
  the EEQ is periodic. One parameter set `m_ewald_params` (chosen from the element-table alpha before the topology
  exists) feeds both the solver and `FFWorkspace::setEwald`; at workspace creation the Gaussian radius from the
  charge-corrected alpeeq is re-checked: when it exceeds `r_cut` (iodine: 19.03 vs 18.46 Bohr from the element
  table) the parameters are re-derived with the larger radius, handed to the solver again and the set-up charges
  are discarded so the first `Calculation()` re-solves them with the final matrix (verbosity 1 reports the raise;
  refused when r_g does not fit below L_min/2).
- A lattice-vector shift of any atom leaves `A_ew` unchanged (MIC pairs, `k . nL = 2 pi n`), so the Phase-2 solve
  inside the topology build (A4c reassembled geometry) and the per-step solve on the input geometry see the same
  matrix; the torn-input checks (d) confirm this.

### Tests (`ctest -R pbc_eeq_ewald`, `test_cases/unit/test_pbc_eeq_ewald.cpp`, 59 checks, 44 s)

Water boxes: 20 H2O random rigid packing in L = 44 Bohr (60 atoms, several atoms within 0.03 A of a face), the same
with H3O+(H2O)2 (Q = +1, 64 atoms, default ensemble charge model); solver tolerance `eeq_ppcg_tol` 1e-12 (1e-14 in
(e)), Ewald `ewald_tolerance` 1e-10 / 1e-12.

| Test | Result |
|---|---|
| (a) dense cholesky vs lu vs matrix-free pcg (tol 1e-13) | charges 3.0e-15 / 1.6e-12 e, energy 5.3e-15 Eh, gradient 2.7e-13 Eh/Bohr; `ppcg` redirected, `batched` refused, `auto` switches at `eeq_ewald_dense_max_atoms`; pcg 4 threads vs 1 thread 3.6e-15 e, 4 threads run-to-run bitwise |
| (b) alpha_E 0.25 / 0.32 / 0.42, r_cut 19 (N_k 10k / 21k / 47k) | neutral: charges 3.8e-12 e, E_tot 3.6e-13 rel, total gradient 1.3e-10 of max; Q = +1: 5.1e-12 e, 9.2e-13, 7.6e-11 |
| (c) large cell L = 200 / 300 Bohr vs the non-periodic path (charges, E_tot, gradient) | water dimer (no dipole): 8.1e-10 e, 2.3e-10 Eh, 8.1e-11 Eh/Bohr, exponent 5.00; benzene 8.6e-11 e, 1.0e-11 Eh, exponent 5.00; 8-water cluster (M = 2.96 e Bohr): E(pbc) - E = -2.294e-6 / -6.796e-7 Eh against the tin-foil dipole term -2 pi M^2/(3V) = -2.293e-6 / -6.795e-7 (rest 3e-11), charges 4.6e-6 / 1.4e-6 e (exponent 2.99), gradient 1.7e-7 Eh/Bohr at L = 300 |
| (d) translation, atom-wise wrap (1 / 28 torn coordinates), 2x1x1 supercell (charges per atom, E, gradient), `hb_min_pair_energy_eh 0` | neutral: 1.2e-14 e / 0 / 3.4e-15; 1.3e-14 / 4.1e-16 / 3.6e-15; 1.4e-14 e / 9.5e-16 rel / 6.9e-15 Eh/Bohr. Q = +1: 1.0e-14 / 1.6e-16 / 7.0e-15; 1.1e-14 / 1.6e-16 / 6.3e-15; 1.1e-14 / 1.6e-16 / 1.2e-14. With the default HB candidate filter (1e-9 Eh) the Q = +1 supercell differed by 1.8e-9 Eh (one dropped candidate), independent of the Ewald tolerance |
| (e) **gate**: analytic total gradient vs central FD of the total energy with re-solved charges, h = 1e-4 / 2e-4 A, `hb_update_force_every 1`, `hb_min_pair_energy_eh 0` | neutral, dense: 2.33e-8 / 9.31e-8 Eh/A (ratio 4.00), Richardson 8.9e-11; neutral, pcg tol 1e-14: 2.33e-8 / 9.31e-8, Richardson 5.8e-11; Q = +1, dense: 2.35e-8 / 9.38e-8 (ratio 4.00), Richardson 1.1e-10; max |g| 1.6e-2 / 7.2e-2 Eh/A. Thresholds 1e-7 and 1e-8 met |
| (f) MIC-only EEQ (minimum-image erf(g r)/r matrix with the instance's fragment constraints, built in the test) vs Ewald EEQ, 50 H2O in L = 44 Bohr (0.014 atoms/A^3, dilute) | max charge deviation 4.8e-3 e (mean 1.5e-3, max q 0.63); E_EEQ(MIC, q_mic) - E_EEQ(Ewald, q_ew) = -3.6e-3 Eh = E_tot(MIC) - E_tot(Ewald) (the workspace Coulomb component is E_EEQ(q) exactly, every other term uses Phase-1 charges); the MIC charges in the Ewald functional lie 1.0e-4 Eh above its minimum. The test-built dense Ewald system reproduces the instance charges to 7e-15 e |
| (f, liquid density, `--mic-vs-ewald 970 58`, not in ctest) 970 H2O in L = 58 Bohr, 2910 atoms, 0.101 atoms/A^3, random rigid packing at 1.7 A minimum distance (947 fragments: some waters merge) - **not an equilibrated liquid** | max charge deviation 2.0e-2 e (mean 3.7e-3, max q 0.70); E_EEQ(MIC) - E_EEQ(Ewald) = E_tot(MIC) - E_tot(Ewald) = -3.4e-2 Eh (-21 kcal/mol, 1.2e-5 Eh per atom); MIC charges in the Ewald functional 1.3e-2 Eh above the minimum; instance charges (pcg) reproduced by the test-built dense system to 7.5e-14 e. This is the B1 measurement: a MIC-only EEQ is not an acceptable approximation of the periodic one |
| (g0) / (g) no cell | `coulomb_ewald` default vs false bitwise (energy, charges, gradient); CLI `-sp` + md5 of `-dump_gradient` on benzene, caffeine, methyl-a-d-glucopyranoside, water8-cluster, hydronium-water2 (+1), urea-water-cluster197, peo201-water500 at 1 and 4 threads: 14/14 md5 identical to the binary before B3 |
| (i) `ctest -R 'gfnff|frag|eeq|simplemd|pbc|ewald|react'` list | identical plus `pbc_eeq_ewald` (117 -> 118) |

(h) Timing (`test_pbc_eeq_ewald --perf-box <nwat> <L> <threads>`; median of 3 `Calculation(true)` at perturbed geometries
after set-up, `eeq_ppcg_tol` 1e-12, `ewald_tolerance` 1e-8; random water packing at liquid density, not a liquid; the
no-cell EEQ (ppcg + localinv) is not timed separately):

| System | Path | EEQ solve (ms, iterations) | Coulomb CPU-sum (ms) | Calculation(true) wall (ms) |
|---|---|---|---|---|
| 970 H2O, 2910 atoms, L = 58 Bohr, 4 threads | no cell | - | 19.1 | 42.2 |
| | cell, pcg, r_cut auto 18.06 (N_k 14008) | 386.7 (15) | 123.2 | 459.2 |
| | cell, pcg, r_cut 25 (N_k 5317) | 477.3 (15) | 128.7 | 554.3 |
| 970 H2O, 1 thread | no cell | - | 18.9 | 120.2 |
| | cell, pcg, r_cut auto | 1364.3 (15) | 119.0 | 1582.2 |
| | cell, pcg, r_cut 25 | 1606.1 (15) | 124.7 | 1828.1 |
| 300 H2O, 900 atoms, L = 44 Bohr (dilute, the smallest accepted cube), 4 threads | cell, dense cholesky (N_k 9174) | 85.0 | 22.1 | 98.7 |
| 500 H2O, 1500 atoms, L = 46.5 Bohr, 4 threads | no cell | - | 5.3 | 15.3 |
| | cell, pcg (N_k 10002) | 138.6 (15) | 44.5 | 168.3 |
| | cell, dense cholesky | 238.1 | 44.2 | 267.6 |

With a cell the step costs 11 x the all-pairs path at 2910 atoms (13 x at 1 thread): the EEQ solve dominates (15 CG
iterations with warm start, each a real-space pass over 0.5 M pairs plus N x N_k complex products; the pair list and
the phase table are built twice per step, once for the solver and once for the workspace). The dense path is slower
than pcg already at 1500 atoms (238 vs 139 ms), so `eeq_ewald_dense_max_atoms` = 1000 is a two-point estimate, not
tuned. CLI `-md` on the same box, 20 steps of 0.5 fs, 4 threads: 14.5 s with the cell (20-21 iterations per step, cold
start every step is not the case: the warm start is the previous step) vs 1.4 s without (both including set-up).

(j) NVE energy conservation (**finding, n = 1 per row**): 60 H2O random packing in L = 44 Bohr (0.014 atoms/A^3, a
dilute gas, not a liquid), dt 0.5 fs, CSVR 300 K for 1 ps (coupling 100 fs, seed 7) then 2 ps NVE from the restart
(velocities included), 4 threads, defaults unless noted; drift = linear fit of E_tot(t) over the 41 records, relative
to the mean kinetic energy of the run. The system is still relaxing (T rises from 300 to 350-430 K in every run).

| run | E_tot drift (Eh/ps) | relative (<E_kin>/ps) | E_tot std / max-min (Eh) |
|---|---|---|---|
| no cell (same molecules as a cluster in vacuum) | +1.1e-4 | +3.7e-4 | 1.8e-4 / 9.1e-4 |
| cell, defaults (Ewald Coulomb + EEQ, dense cholesky, tol 1e-8) | -4.0e-3 | -1.1e-2 | 2.4e-3 / 8.6e-3 |
| cell, `coulomb_ewald false` (plain Coulomb/EEQ with wrap, the state before B3) | +2.2e-2 | +6.4e-2 | 1.3e-2 / 4.2e-2 |
| cell, `hb_update_force_every 1`, `hb_min_pair_energy_eh 0` | -5.1e-3 | -1.5e-2 | 3.1e-3 / 1.0e-2 |
| cell, `ewald_tolerance 1e-10` | -8.5e-3 | -2.3e-2 | 5.1e-3 / 1.5e-2 |
| cell, `ewald_r_cut 21` (N_k smaller) | -6.1e-3 | -1.7e-2 | 3.6e-3 / 1.3e-2 |
| cell, dt 0.25 fs (1 ps) | -9.4e-3 | -2.8e-2 | 2.9e-3 / 1.1e-2 |
| cell, `dispersion_tail false` | -7.5e-3 | -2.1e-2 | 4.6e-3 / 1.6e-2 |
| cell, `hbond false` | -2.5e-3 | -8.5e-3 | 1.6e-3 / 5.0e-3 |
| cell, `dispersion false` | -4.4e-3 | -1.3e-2 | 2.7e-3 / 8.1e-3 |

With the cell the total energy is lost at 2.5-9e-3 Eh/ps (1-3 % of E_kin per ps), 20-80 x the cluster, in every
variant; the periodic Coulomb/EEQ makes it 5 x better than the plain-distance Coulomb with wrapping, but the loss is
not removed by a tighter Ewald tolerance, a smaller time step, the HB lists, the D4 tail or dispersion as a whole,
so it is not located in a single term by these runs (trajectories diverge, so the variants are different n = 1
samples, not a controlled series). The static gradient (test e) is consistent to 2e-8 Eh/A, so the candidate is a
discontinuity between steps (pair-list, minimum-image or wrap events), not the force at a fixed geometry. A 0.3 ps
run with per-step diagnostics (`md_diagnostics`, dump every step) shows: per-step |dE_tot| noise 2.4e-4 Eh (the same
without the cell: the 0.5 fs step at 400 K), the loss comes in plateaus and drops (E_tot - E_tot(0) at 0.1 ps
intervals: -1.6e-3, -0.9e-3, -0.9e-3, -1.1e-3, -3.7e-3, ... -7.9e-3), 11 fragment wraps in 600 steps with dE_tot at the
wrap steps within the noise (-3.4e-4 .. +3.3e-4), HB-list changes every 10 steps likewise; `-coulomb false` runs are
unphysical with and without the cell (-0.2 Eh/ps, discarded). **Open (TODO.md): not located; next a longer per-step
series with a work balance (dE_pot vs -F.dx per step) to find the steps where energy and force disagree.** **Resolved in D1 below (2026-11-02): the loss was the COM/rotation removal, and these runs were fresh starts because `-restart_file` had rejected the restart file.**

## D1 NVE energy loss with a cell: cause and fix (2026-11-02, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `7b4beca1`), CPU only. The B3 finding (j) is reproduced and root-caused: **the energy was
removed by SimpleMD's periodic COM/rotation removal, not by the periodic energy surface.** Two bookkeeping defects, both
in `src/capabilities/simplemd.cpp`:

1. **`RemoveRotation()` (remove_com_mode 1, the default, every `remove_com_motion` = 100 fs) subtracts the rigid
   rotation `omega = I^-1 L` of the whole system together with the COM velocity.** In a cluster in vacuum the angular
   momentum L is a constant of motion, so after the first call nothing is taken out. With a cell L about the COM of the
   wrapped coordinates is not conserved (no rotational symmetry; every wrap moves a molecule by a lattice vector), so
   every call removed `1/2 L^T I^-1 L` of kinetic energy: the loss is a staircase with one step per 200 MD steps. Fix:
   with a cell (`m_has_pbc`) the translation only is removed, and the DOF count subtracts 3 instead of 6 (the reported
   T of a periodic run changes by that ratio, e.g. 534 -> 537 DOF for 180 atoms). Modes 2/3 (per-fragment rotations,
   opt-in) are unchanged.
2. **`-restart_file` rejected curcuma's own restart files.** The files SimpleMD writes wrap the state as
   `{"MD": {...}}`; the `-restart_file` path handed the outer object to the validator, which failed on the missing
   `method` and **silently started a fresh run** (`[ERROR] Restart file validation failed ... Starting fresh simulation
   instead`, exit 0). Every B3 "NVE after CSVR" run therefore started from the random packing with fresh 300 K
   velocities and no thermostat, which is why T rose to 350-430 K. Fix: the wrapped layout is unwrapped. Note that a
   restart file carries `thermostat`, `dT`, `MaxTime`, `T` and overrides the CLI; the protocol below edits the file.

### Evidence (pre-fix binary saved, 2 ps NVE from a fresh start at 300 K, dt 0.5 fs, 3 threads; drift = linear fit of
E_tot over the last 1.5 ps; `max-min` of E_tot over the run; box60 = the B3 packing, box60b = a second packing, seed 2)

| experiment | run | drift (Eh/ps) | max-min (Eh) | comment |
|---|---|---|---|---|
| 1 reproduce, cell, default removal (n = 3) | box60 seed 7 / seed 11 / box60b seed 7 | -5.9e-3 / -9.6e-3 / -6.5e-3 | 1.2e-2 / 1.7e-2 / 1.8e-2 | B3 (j): -4.0e-3 / 8.6e-3 |
| 1 no cell, default removal (n = 3) | same three | -1.0e-5 / -1.7e-5 / -1.5e-5 | 1.3e-3 / 1.6e-3 / 1.5e-3 | reference level |
| discriminator: cell, `-remove_com_motion 0` (n = 3) | same three | -4.2e-5 / -3.6e-5 / -6.8e-5 | 1.28e-3 / 1.56e-3 / 1.52e-3 | **equal to no cell** (max-min to 2 digits) |
| no cell, `-remove_com_motion 0` (n = 3) | same three | -1.9e-5 / -2.5e-5 / -2.4e-5 | 1.26e-3 / 1.56e-3 / 1.52e-3 | removal is inert in vacuum |
| A giant cell L = 300 Bohr, cluster centred, no wrap, no image interaction | default removal / off | -5.9e-6 / -1.4e-5 | 1.40e-3 / 1.40e-3 | cell plumbing alone does not drift: L is conserved without wraps |
| B ladder, L = 44 Bohr, 3 ps, default removal | 1 / 2 / 8 H2O | -2.8e-7 / -1.8e-6 / **-3.0e-3** | 1.9e-5 / 9.8e-5 / 5.8e-3 | drift appears with collisions + wraps |
| B ladder, removal off | 1 / 2 / 8 H2O | -2.9e-7 / +1.4e-8 / -2.3e-6 | 1.9e-5 / 8.8e-5 / 1.5e-4 | 8 H2O: 1300 x smaller |
| C dt scaling, cell, default removal, box60 seed 7 | dt 0.25 / 0.5 / 1.0 fs | -5.6e-3 / -5.9e-3 / -9.7e-3 | 1.0e-2 / 1.2e-2 / 1.9e-2 | loss per time, not per step: an event every 100 fs |
| per-step series (B3 diagnostics run, 600 steps, cell) | dE_tot at steps 200 and 400 | -6.7e-4 / -6.9e-4 Eh | | the two largest of 600 steps (per-step noise 2.4e-4); no cell: rank 70 and 143. Sum -1.36e-3 Eh / 0.3 ps = -4.5e-3 Eh/ps = the B3 drift |

(D, the per-term path integral along the trajectory, was not needed: with the removal off the periodic run conserves
energy as well as the cluster, so no term of the energy surface jumps between steps at the 1e-5 Eh/ps level resolved
here. The B3 variants (Ewald tolerance, HB lists, dispersion, tail, dt) all lost energy because the removal ran in all
of them.) `-coulomb_ewald false` (B3: +2.2e-2 Eh/ps) was not re-measured; that path is non-periodic with wrapping and
is not expected to conserve energy.

### After the fix (release binary, same protocol)

| run | drift (Eh/ps) | max-min (Eh) |
|---|---|---|
| cell, default removal, box60 seed 7 / seed 11 / box60b seed 7 | -3.6e-5 / -1.7e-5 / -2.6e-5 | 1.28e-3 / 1.56e-3 / 1.53e-3 |
| no cell, box60 seed 7 (bitwise the same run as before the fix) | -1.0e-5 | 1.26e-3 |
| 8 H2O in L = 44 Bohr, default removal, 2 ps | -3.1e-6 | 1.5e-4 |
| **B3 protocol** (CSVR 300 K 1 ps, then NVE 2 ps from the restart file with `thermostat` edited to none; `-restart_file` now loads it, verified: no "fresh simulation" line, periodic EEQ line present): cell seed 7 / seed 11 | -5.3e-5 / -8.5e-5 | 5.0e-4 / 5.8e-4 |
| B3 protocol, no cell, seed 7 | -4.0e-6 | 4.1e-4 |
| **dense box**: 420 H2O in L = 44 Bohr (1260 atoms, 0.100 atoms/A^3, random packing at 1.7 A, 378 fragments after CSVR, i.e. merged waters, not a liquid; 300 L-BFGS steps, not converged; CSVR 1 ps raised E_tot by 1.57 Eh while the packing relaxed), NVE 2 ps dt 0.5 fs, pcg EEQ 23 iterations, default removal, n = 1 | -6.4e-4 (-3.7e-4 <E_kin>/ps; <E_tot> per 0.5 ps window -144.17989 / .18020 / .18054 / .18082, monotone) | 3.0e-3 (std 5.4e-4) |
| dense box, same restart, `-hb_update_force_every 1 -hb_min_pair_energy_eh 0` (HB/XB lists every step), dt 0.5 fs, n = 1 | **-2.9e-5** (windows -144.17976 / .17974 / .17983 / .17973, flat) | 2.0e-3 (std 3.4e-4) |
| dense box, same restart, dt 0.25 fs, default HB refresh (every 10 steps = 2.5 fs), n = 1 | -3.5e-4 (windows monotone) | 9.3e-4 (std 2.3e-4) |

**Second, smaller source at liquid density: the HB/XB list refresh interval.** With the default `hb_update_force_every 10` the
HB/XB pair lists are stale for 10 steps and the energy surface jumps at every refresh; at 1260 atoms with many hydrogen bonds
this costs -6.4e-4 Eh/ps (-3.7e-4 <E_kin>/ps) and goes away with the lists rebuilt every step (-2.9e-5 Eh/ps, 22 x less, n = 1
each from the same restart). Halving dt halves the loss (refresh every 2.5 fs instead of 5 fs), it does not quarter it, so it
is not the Verlet dt^2 error. This is not a cell effect (the same refresh runs without a cell, where B3 measured it on 60 H2O as
below the noise) but a density effect that the periodic box makes reachable; the default is a speed trade-off (TODO.md,
section GPU: 10 -> 20 saves 1.1 ms/step at 7320 atoms) and is left to the operator. The dilute protocol residual of
-5e-5..-8.5e-5 Eh/ps above is of the same order as this run with per-step lists and was not pursued further.

The cell run now conserves energy at the level of the cluster (fluctuation amplitude max-min equal to 2 digits in every pair; the
fitted drift of the protocol runs, -5e-5..-8e-5 Eh/ps against -4e-6 without the cell, is n = 2 against n = 1 of a still relaxing
system (T 300 -> 367-379 K with the cell, 308 K without, because the periodic gas condenses) and is 50-100 x below the B3 loss;
whether it is a residual cell effect (HB-list refresh every 10 steps, D4/repulsion cutoffs) is not resolved here). The per-cause
balance: the rotation removal accounts for the whole B3 loss within the resolution of these runs; the restart defect changed the
starting state (fresh 300 K velocities on the unrelaxed packing instead of the CSVR state) but not the loss rate (fresh-start
runs lose -5.9e-3..-9.6e-3 Eh/ps with the pre-fix binary, like the B3 runs that were, unknowingly, fresh starts too).

### Checks

- `test_cases/unit/test_pbc_md_com_removal.cpp` (ctest `pbc_md_com_removal`, 0.04 s): 8 registry waters in L = 44 Bohr, a restart
  velocity field thermal + COM drift + rigid rotation, one step of 1e-6 fs with the removal due at step 0; with a cell the
  velocities after the step equal v_in - v_com (8e-6, tolerance 3e-5 from the 6-decimal restart format) and E_kin equals
  E_kin(v_in - v_com); without a cell v_in - v_com - omega x (r - r_com) as before; the `{"MD": {...}}` file is accepted.
  Before the fix 5 of its 8 checks failed (cell: 2.4e-2 velocity residual, E_kin 0.0384 instead of 0.0816 Eh; wrapped restart
  ignored, 1.1e-1), after it 8/8 pass.
- No-cell path bitwise: `-sp -dump_gradient` md5 and printed energy on benzene, caffeine, methyl-a-d-glucopyranoside, water8,
  hydronium-water2 (+1), urea-water197, peo201-water500 at 1 and 4 threads: 14/14 identical to the pre-fix binary; the no-cell
  NVE run (box60 seed 7) is bitwise the same trajectory before and after (drift -1.00e-5, max-min 1.26e-3, E(end)-E(0) +6.65e-4).
- `ctest -j1 -R 'gfnff|frag|eeq|simplemd|pbc|ewald|react'`: 118 passed, 1 skipped (GPU), 0 failed (119 incl. the new test).
- `scripts/check_docs.py`: the 9 known errors only.

Not covered: `-gpu` runs (the wrap and the removal are CPU code in SimpleMD, so the same fix applies, not run), modes 2/3 under a
cell (per-fragment rotation removal takes rotational energy out of every molecule, with or without a cell; opt-in, unchanged),
the spherical container (`wall_type spheric` has no lattice, `m_has_pbc` false: the full removal still runs there).

## P1 CPU performance of the periodic EEQ and Ewald Coulomb (2026-11-03, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `a4587768`), CPU only. Harness: `test_pbc_eeq_ewald --perf-box <nwat> <L> <threads> [json]`
(median of 3 `Calculation(true)` after a 1e-3 A move of one atom, warm start), `--perf-md <nwat> <L> <threads> <nsteps> <dt_fs>
[json]` (velocity-Verlet NVE on the same random packing from 300 K, median over steps 6..20; **the packing is unrelaxed, so charges
change by ~2e-2 e per 0.5 fs step: a hot, non-equilibrium system, n = 1 box per size**), `--tol-table`, `--write-box`. Boxes:
970 H2O in L = 58 Bohr (2910 atoms) and 1916 H2O in L = 72.8 Bohr (5748 atoms), liquid density, 4 threads unless noted,
nothing else running on the machine.

### Step 0: profile before (B3 code, 2910 atoms, 4 threads, `Calculation(true)` 460 ms)

EEQ pcg 388 ms = 16 operator calls (15 iterations + start residual) of 23.5 ms wall each; per call (1-thread run) real space
31 ms (erfc/exp of 534k pairs every call) and reciprocal 53 ms (std::complex products of 14008 k vectors x 2910 atoms); pair
list and phase table 11 ms, built twice per step (solver and workspace). Workspace Coulomb 123 ms CPU (real 33, recip 78,
setup 11). Rest of the step ~45 ms. Thread scaling 1/2/4/8: 1576/848/460/318 ms.

### What was changed

1. **Shared geometry data** (`EwaldGeometryCache`, `EwaldSum::buildGeometryCache`): pair list, kernel K_ij = [erfc(a r) -
   erfc(g_ij r)]/r and dK/dr per pair, phase table, built once per geometry by whoever sees it first (the EEQ solve) and reused
   by the workspace (key: positions, alp, alpha/r_cut/k_max, compared bit for bit; GFNFF hands both the same object). The EEQ
   matrix-vector product is a gather over the stored kernel. Bitwise the same numbers as the per-call kernel (one non-inlined
   `pairKernel`; ctest (p1)). Cap `ewald_kernel_cache_mb` (1024 MB, 16 bytes per pair), above it the per-call kernel.
   Memory: 21.7 MB at 2910 atoms (850k pairs), 94.9 MB at 5748 (2.36M pairs); peak RSS of the harness process 800 vs 731 MB
   (no cell) at 2910, 2110 vs 2071 MB at 5748.
2. **Reciprocal pass vectorised**: structure-of-arrays phase table (n-major cos/sin), real arithmetic, `omp simd` loops
   (S(k) reduction reassociated, so bits differ from B3 at 1e-16 relative). 53 -> 15 ms per call at N_k 14008 (1 thread).
3. **Pair search** with cells of edge r_cut/2 (125 neighbours) and threads over (cx, cy) columns, concatenated in order (same
   list for every thread count): 5.7 -> 3.7 ms at 4 threads.
4. **r_cut heuristic** (`GFNFF::ewaldParameters`): 15 -> 17.5 sqrt(L_min/40) (scan below).
5. **Block preconditioner** (`eeq_ewald_precond block`, default): per fragment (chunks <= 32 atoms) the inverse of the
   non-periodic Gaussian EEQ block erf(r/sqrt(alp_i+alp_j))/r at minimum-image distance (derivation in `eeq_solver.cpp`,
   `PeriodicBlockPreconditioner`); `jacobi` restores B3. 15 -> 7 iterations (perf-box), 20.5 -> 11 (MD).
6. **Own tolerance** `eeq_ewald_tol` 1e-9 (table below); `eeq_ppcg_tol` (1e-12) now only for the non-periodic ppcg.
7. **Start vector** `eeq_ewald_extrapolation aspc`, order 1 (`curcuma::extrapolation::weights`, as `eeq_extrapolation` of the
   non-periodic pcg): 7 -> 6 iterations in MD.

### Time per step (ms, `Calculation(true)` wall, 4 threads; EEQ ms and pcg iterations in brackets)

| stage | 2910 perf-box | 2910 MD | 5748 perf-box | 5748 MD |
|---|---|---|---|---|
| before (B3) | 460.9 (388.9, 15 it) | not measured in the harness; CLI `-md` 640 ms/step | 1205.7 (1012.0, 14 it) | - |
| 1-3: shared data, stored kernel, SIMD reciprocal, pair search (Jacobi, 1e-12, r_cut 18.06 / 20.24) | 118.2 (79.7, 15) | - | 305.6 (208.1, 14) | - |
| 4: r_cut 21.07 / 23.61 | 106.1 (69.2, 15) | 119.8 (82.0, 20.5) | 256.6 (165.3, 14) | 327.7 (225.9, 21) |
| 5: block preconditioner | 80.4 (43.8, 7) | 92.8 (55.4, 11) | 216.8 (123.7, 8) | 237.8 (145.3, 11) |
| 6: `eeq_ewald_tol` 1e-9 | 70.0 (33.1, 4) | 80.7 (42.6, 7) | 183.1 (90.3, 4) | 204.7 (114.2, 7) |
| 7: aspc start (final defaults) | - | 77.5 (39.5, 6) | - | 200.9 (106.7, 6) |
| no cell (reference) | 46.4 | - | 133.6 | - |

CLI `-md` on the 2910-atom box (4 threads, dt 0.5 fs, difference of 20- and 40-step runs, median of 3): 640 -> 114 ms/step
(no cell: 47 ms/step). The CLI step with a cell is ~37 ms longer than `Calculation(true)` in the harness (77.5 ms); that part is
SimpleMD, not located here. Thread scaling after (MD, 2910 atoms) 1/2/4/8: 248/133/79/59 ms (perf-box: 230/123/70/48).

### Stage 2 scan: r_cut with the stored kernel (MD, 6 pcg iterations per step, 4 threads)

| 2910 atoms, L 58: r_cut | 15 | 18.06 | 21 | 24 |
|---|---|---|---|---|
| step (EEQ) ms | 102.2 (58.7) | 82.2 (42.1) | 77.4 (40.5) | 83.4 (47.0) |

| 5748 atoms, L 72.8: r_cut | 17 | 20.24 | 23 | 26 |
|---|---|---|---|---|
| step (EEQ) ms | 270.4 (163.9) | 219.3 (118.8) | 199.8 (107.3) | 208.8 (114.9) |

The optimum moves with the iteration count (more iterations favour larger r_cut: with 15 Jacobi iterations it was flat over
21-27 Bohr at 2910 atoms); 17.5 sqrt(L/40) gives 21.07 / 23.61. Energies move by up to 1e-7 Eh over the scan (2910 atoms;
`ewald_tolerance` 1e-8 is a per-term criterion), i.e. 3e-11 Eh per atom.

### Tolerance table (pcg from a cold start, `eeq_ewald_tol` vs 1e-14; criterion fixed beforehand: max |dG| <= 1e-8 Eh/Bohr, |dE| <= 1e-9 Eh per atom)

| tol | 50 H2O L 44: max abs dq / dE per atom / max abs dG | Q = +1 (64 atoms) | 970 H2O L 58 |
|---|---|---|---|
| 1e-6 | 5.8e-6 / 8.6e-13 / 1.2e-6 | 6.0e-6 / 3.6e-13 / 8.8e-7 | 7.1e-6 / 3.2e-13 / 1.2e-6 |
| 1e-8 | 2.1e-9 / 5e-17 / 3.7e-10 | 2.3e-8 / 4e-17 / 4.2e-9 | 1.1e-7 / 3e-16 / **1.9e-8** |
| 1e-9 | 2.1e-9 / 5e-17 / 3.7e-10 | 2.9e-10 / 4e-17 / 7.4e-11 | 1.2e-8 / 4e-16 / 1.9e-9 |
| 1e-10 | 5.0e-11 / 2e-17 / 8.3e-12 | 2.9e-10 / 4e-17 / 7.4e-11 | 1.4e-9 / 8e-16 / 2.4e-10 |
| 1e-12 | 9.1e-13 / 7e-17 / 1.7e-13 | 5.4e-12 / 1e-17 / 1.0e-12 | 3.2e-12 / 1e-16 / 4.9e-13 |

(e, Eh, Eh/Bohr.) 1e-8 fails the gradient criterion on the 2910-atom box, 1e-9 is the largest tolerance that meets both on all
three (n = 3 systems; the energy error is quadratic in the charge error and far below the limit everywhere). Not checked at
5748 atoms (tol relative to |b| ~ sqrt(N)).

### Discarded or not done

- Larger extrapolation orders: aspc order 3 gave 7 -> 5-5.5 iterations in MD at 0.5 fs but 6.05 -> 6.65 on an irregular
  steepest-descent sequence (n = 1 each); gauss order 2 7 -> 6-7. Order 1 was never worse in these samples (MD at 1 fs: 7 -> 6.9).
  An earlier measurement with 1 iteration per step came from a harness bug (time step 1000x too small) and is void.
- A preconditioner coupling neighbouring fragments (cluster blocks, localinv-like) was not tried: ~1 ms per pcg iteration left.
- Faster erfc in the kernel build (39 of 56 ms set-up at 5748 atoms, 4 threads) and a Verlet skin for the pair list (17 ms)
  were not done; the set-up is now the largest Ewald item (perf-box 5748: set-up 56 ms, 5 operator calls 38 ms wall = 7.7 ms
  each, of which reciprocal 26.5 ms CPU per call; workspace Coulomb 43 ms CPU).
- PME: not justified by the profile after stages 1-6 (reciprocal ~50 of 183 ms wall at 5748 atoms, set-up 56); kept as an
  option for larger cells: O(N log N) with B-spline charge assignment, its own interpolation error to validate on forces,
  estimated 1-2 weeks including tests.
- `eeq_ewald_dense_max_atoms` (1000) unchanged: the matrix-free pcg is 7.2 vs 84 ms (dense) at 900 atoms and 1.5 vs 16 ms at
  300 atoms, but the small-box ctests run the exact dense path at 1e-14-level thresholds; decision left to the operator.
- **Found, not caused by P1, not fixed**: an explicit `-solve_method pcg` changes the energy of the 2910-atom box by 3.75e-6 Eh
  against `auto` (which also runs the periodic pcg there), and cholesky vs pcg at 900 atoms differ by 5.4e-7 Eh, with angle,
  dihedral and HB terms differing (1.5e-7, 5e-8, 4e-8 Eh): the solver name also selects a non-periodic EEQ solve inside the
  topology set-up whose charges feed parameters. The same in the B3 binary (TODO.md).

### Acceptance

- (i) before vs after with defaults (CLI `-sp`, 4 threads): 50 H2O L 44 (dense) dE 2.6e-9 Eh (1.8e-11 per atom), max |dG|
  3.4e-10 Eh/Bohr, max |dq| 9e-10 e; Q = +1 (64 atoms, dense) 8.6e-10 / 1.5e-10 / 6e-10; 970 H2O (pcg) 2.5e-8 Eh (8.6e-12 per
  atom) / 1.8e-9 / 1.2e-8. With the old r_cut forced and `eeq_ewald_tol` 1e-14 the solver alone accounts for <= 1.5e-12 Eh/Bohr
  (2910 atoms, before the r_cut change; the forced r_cut differs from the automatic one in the 4th decimal, which moves the
  dense boxes by 1e-9 Eh).
- (ii, vii) `ctest -j1 -R 'gfnff|frag|eeq|simplemd|pbc|ewald|react'`: 118 passed, 1 skipped (GPU), 0 failed, the same 119 tests.
  Thresholds unchanged; `pbc_eeq_ewald` (a) and (e) set `eeq_ewald_tol` instead of `eeq_ppcg_tol` (the pcg they test now reads
  it), and gained (p1) (stored == per-call kernel bitwise, shared data reused, 4-thread pair list == 1-thread list, block vs
  Jacobi 1.6e-12 e, aspc vs plain start 1.6e-12 e) and an FD gate with the P1 defaults.
- (iii) FD gate (20 H2O): dense Richardson 9.4e-11 Eh/A, pcg tol 1e-14 5.1e-11, **pcg with the P1 defaults (tol 1e-9) 7.5e-11**,
  Q = +1 dense 9.3e-11 (limit 1e-8).
- (iv) no cell: `-sp` energy (16 digits in the dump) and `-dump_gradient` md5 of the 7 structures at 1 and 4 threads: 14/14
  identical to the binary before P1.
- (v) 4 threads run to run: bitwise (energy, gradient, charges) on all three boxes.
- (vi) NVE, 60 H2O in L = 44 Bohr, 1 ps from 300 K, dt 0.5 fs, 3 threads, fit 0.25-1 ps: -5.2e-5 Eh/ps (default dense and forced
  pcg) against -4.0e-5 with the binary before P1 (n = 1 each; trajectories differ through the r_cut change); D1 level met.
- Not tested: charged or heavy-element boxes at the new r_cut beyond the ctest set, more than 5748 atoms, `ewald_kernel_cache_mb`
  fallback at production size (only the 60-atom bitwise check), relaxed liquids. Human production testing pending.

## A6 wall_potential pbc with a periodic energy: the stopgap is skipped (2026-11-03, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `40461124`), CPU only.

### Step 0: how the two wraps interacted (read from the code before the change)

- Two independent mechanisms, both called after every integration step in this order: `applyPeriodicBoundaryConditions()`
  (only with a cell, `m_has_pbc` from `Molecule::hasPBC()`; fragment COM into fractional [0,1), no clearance test, every
  method) and then `wrapIntoContainer()` (only `wall_potential pbc`, `m_wall_wrap`; box per axis or sphere antipodal, with
  clearance test `kWrapClearance` 1.8 A and reflection), at the end of `SimpleMD::Verlet()` and `SimpleMD::Rattle()`.
- Box equal to the cell `[0, L)`: the cell wrap runs first, so the container wrap never finds a COM outside; the
  stopgap could not fire (it still ran the O(N^2) fragment search whenever an atom of a straddling molecule was outside).
- Box different from the cell (e.g. the auto-sized box, or `[-L/2, L/2)`): both wrap against each other. Measured in
  ctest `pbc_md_wrap_pure` (d), GFN-FF with `coulomb_ewald false`, dimer across the face, 200 steps: 398 cell wraps and
  400 container wraps, i.e. one fragment moved by +L and -L in every step.
- Periodic energy: only the native CPU GFN-FF with a cell and `coulomb_ewald` (periodic Ewald EEQ active); the GPU
  GFN-FF with a cell is replaced by the CPU method in `EnergyCalculator::setMolecule()` (`periodicCellNeedsCpu`); no
  other method reads the cell (`grep` for MIC/`m_unit_cell` in `energy_calculators/`: GFN-FF files only, plus a copy
  in `xtb_fragment_scf.cpp`).

### Decision and change

- **One source of truth**: `ComputationalMethod::hasPeriodicEnergy()` (default false), overridden by
  `GFNFFComputationalMethod` (`hasPBC() && periodicEEQActive()`), forwarded by `EnergyCalculator::hasPeriodicEnergy()`.
  A virtual on the method rather than a name comparison in SimpleMD, like `gfnffInstance()`/`requestsCpuFallback()`; it is
  asked after `setMolecule()`, so the GPU-to-CPU swap has happened. A new virtual changes the `ComputationalMethod`
  vtable: GPU plugins must be rebuilt from the same tree (always the case for in-tree builds).
- **Periodic case** (`InitialiseWalls()`): box container + `wall_potential pbc` + cell + periodic energy sets
  `m_wall_wrap_by_cell`; `wrapIntoContainer()` then returns at once, and the wrap is the existing cell wrap (one lattice
  vector per fragment, no clearance test, no reflection). Chosen over a second wrap implementation in the container path:
  one function, and the cell (not the container bounds) defines where the periodic images are. Bounds that differ from
  `[0, L)` produce a warning that they are ignored. Log at verbosity 2 names the method.
- **Everything else unchanged bit for bit**: no cell, a sphere, a method without periodic energy, `coulomb_ewald false`.
  New: one warning at verbosity 1 per run (`wall_potential pbc: energy of method X is not periodic (...); wrap uses the
  clearance/reflection fallback`), and with a cell but no container and a non-periodic method a warning that fragments are
  wrapped without minimum image in the energy. Counters (`cellWrapCount()`, `containerWrapCount()`,
  `containerRefusedCount()`, `containerReflectionCount()`, `wallWrapByCell()`) for tests; they do not touch the numbers.

### Tests (`ctest -R pbc_md_wrap_pure`, `test_cases/unit/test_pbc_md_wrap_pure.cpp`, < 1 s)

| part | setup | result |
|---|---|---|
| (a) source of truth | gfnff + cell / + `coulomb_ewald false` / no cell / uff + cell | true / false / false / false |
| (b) pure wrap | registry water dimer split across the +x face of L = 44 Bohr (MIC O-O 2.8 A), drift +x, NVE 2 ps, dt 0.5 fs, COM removal off, box = cell | 2 cell wraps, 0 container wraps, 0 refused, 0 reflections; drift 9.4e-8 Eh/ps (fit over 2 ps), max-min 1.0e-5 Eh; max abs(P - P0) 1.7e-14 (abs(P0) 0.72) |
| (c) fallback | uff / gfnff, no cell, water landing on a water | refused 20 / 1, reflections 4 / 1 |
| (d) box != cell | `[-L/2, L/2)`, periodic / `coulomb_ewald false` | 1 cell wrap, no container action / 398 + 400 wraps (pre-A6 behaviour, kept) |

(b) does **not** show that the old clearance test would have refused there: the smallest raw contact right after a wrap
was 2.78 A (> 1.8 A). In the box = cell geometry the old code would not have refused either (cell wrap first, see Step 0);
what changes there is that the container path is no longer entered. (n = 1 trajectory per part.)

### Acceptance

- No cell, no container: `-sp` energy (16 digits in the dump header) and `-dump_gradient` md5 of benzene, caffeine,
  methyl-a-d-glucopyranoside, water8, hydronium-water2 (+1), urea-water197, peo201-water500 at 1 and 4 threads: 14/14
  identical to the binary before A6.
- `ctest -j1 -R 'gfnff|frag|eeq|simplemd|pbc|ewald|react|restart|md_'` (124 tests incl. the new one): 3 full runs; runs 2
  and 3 (the last after the final build) 0 failed, 1 skipped (GPU); run 1 had one failure, `cli_gfnff_09_farfield_shortcut`
  (MD part "two waters ... Epot series identical"), reproduced 1 of 15 in isolation. **Cause (not A6)**: the test runs
  two MD jobs (`farfield_shortcut` true/false) without `-md.seed`; the default `seed -1` takes the wall clock in seconds
  (`simplemd.cpp` `if (m_seed == -1)`, two sites), so two runs that straddle a second boundary start from different
  velocities. `cli_simplemd_16_md_seed` asserts "runs without a seed are identical" with the same default and has the same
  flake. Fix (not applied, test change): pass a fixed `-md.seed` in both tests, or drop that assertion of the seed test.
- Not tested: `-gpu` with a cell (falls back to CPU, then periodic; not run), the RATTLE integrator with a cell,
  CG systems with a cell, triclinic cells (refused by GFN-FF). Human production testing pending.

### Options not implemented (behaviour changes, operator decision)

1. **Box becomes the cell**: `wall_type rect` + `wall_potential pbc` without a cell sets `m_unit_cell` from the box (origin
   shifted to 0) before the energy set-up, so GFN-FF runs periodically. Effort ~0.5 day (Molecule/Mol cell before
   `EnergyCalculator::setMolecule`, origin shift of the coordinates; the A1 cutoff check then applies, e.g. the HB/XB
   radius 21.2 Bohr needs L > 42.4 Bohr (22.4 A) with the defaults, not re-checked for every cutoff after A2-A3c). Risk: every existing `wall_potential pbc` run changes
   energy (images interact), and small auto-sized boxes turn into a hard error.
2. **Cell + non-periodic method**: skip the container wrap whenever a cell is set (removes the +L/-L ping-pong of (d)),
   or refuse a cell for methods without periodic energy. Effort < 1 h each; changes trajectories of such runs.
   **Done differently in A6b**: the cell wrap is skipped instead (operator decision).

## A6b a unit cell with a non-periodic method: no cell wrap, vacuum COM removal (2026-11-03, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `15e91ab0`), CPU only. Closes option 2 of A6 (operator decision): with a cell and a method whose
energy does not read the cell (`hasPeriodicEnergy()` false: native gfn1/gfn2, xtb/tblite/Ulysses, uff, qmdff, cg, GFN-FF
with `coulomb_ewald false`) `applyPeriodicBoundaryConditions()` is skipped. Reason: such an energy is translation invariant,
so the wrap changes nothing physical, but it moves a fragment by L without a clearance test (overlap risk) and, with a
container box different from the cell, fought the container wrap (+L/-L every step, A6 (d)).

- **One flag**: `SimpleMD::m_periodic_energy = m_has_pbc && EnergyCalculator::hasPeriodicEnergy()`, set once in
  `SimpleMD::Initialise()` right after `setMolecule()` (GPU-to-CPU swap done), read by `InitialiseWalls()` (A6 decision, no
  second query), `InitConstrainedBonds()` (DOF), `applyPeriodicBoundaryConditions()` and `RemoveRotation()`.
- **Callers of the cell wrap** (`grep`): only `SimpleMD::Verlet()` and `SimpleMD::Rattle()` (end of step, before
  `wrapIntoContainer()`); no call at initialisation or restart, so the gate in the function covers every path.
- **D1 consequence**: the D1 translation-only removal and -3 DOF were needed because wraps break angular momentum
  conservation. Without a wrap the run is the vacuum run in a cell that nobody reads, L is conserved, so the vacuum
  removal (translation + rotation, -6 DOF; -5 for 2 atoms, -3 for 1) applies again. Measured: uff/gfn2 + cell give the
  E_tot series bitwise equal to the same run without a cell (below).
- **Warning** at verbosity 1, once per run (`InitialiseWalls()`): `unit cell set but the energy of method X is not periodic:
  coordinates are not wrapped into the cell` (replaces the A6 warning "fragments are wrapped ... without minimum image",
  now also printed together with the `wall_potential pbc` fallback warning). The container wrap is unchanged.
- Behaviour change also for pure-CG systems with a cell: the per-bead wrap is gone (the cg force field does not read the
  cell either, `grep` for cell/MIC in `energy_calculators/`).

| cell | periodic energy | cell wrap | container wrap (`wall_potential pbc`) | COM removal (mode 1) | DOF removed |
|---|---|---|---|---|---|
| no | (always false) | no | clearance/reflection fallback | translation + rotation | 6 |
| yes | no (all but CPU GFN-FF + `coulomb_ewald`) | **no (A6b; before: yes)** | clearance/reflection fallback | **translation + rotation (A6b; before: translation)** | **6 (before: 3)** |
| yes | yes | yes (fragment COM into [0,1)) | skipped in a box (A6), fallback in a sphere | translation only (D1) | 3 |

### Tests

- New ctest `pbc_md_nonperiodic_cell` (`test_cases/unit/test_pbc_md_nonperiodic_cell.cpp`, 0.1 s): (a) registry water dimer at
  the origin, cell L = 44 Bohr ([0, L), so the dimer starts outside it), `wall_potential pbc` box [-L/2, L/2), fresh 300 K start
  seed 7, NVE 200 steps of 0.5 fs, removal every 10 fs: uff and gfn2 each 0 cell wraps, 0 container actions, DOF 12 = 3N-6, and the
  E_tot series **bitwise equal** to the run without a cell, so the drift is that of the vacuum reference (gfn2 9.6e-6 Eh/ps,
  max-min 6.6e-5 Eh, n = 1); (b) uff, box = cell, dimer drifting through +x: 0 cell wraps, 2 container wraps; (c) gfnff + cell DOF
  15 = 3N-3 (periodic), gfnff + cell + `coulomb_ewald false` DOF 12.
- `pbc_md_wrap_pure` (d) fallback (`coulomb_ewald false`, box [-L/2, L/2)): 0 cell wraps, 2 container wraps in 200 steps
  (A6: 398 + 400); two checks added. `pbc_md_com_removal` (gfnff, periodic) unchanged and green.
- **Finding, not A6b**: uff NVE on the water dimer is not energy conserving at dt 0.5 fs, with or without a cell and with the
  binary before A6b as well (CLI `-md -method uff`, `water-dimer.roo2.8`, seed 7: E_tot 0.009 -> 0.39 Eh in 100 fs, identical
  table old/new). The uff part of the test therefore checks bitwise equality with vacuum, not a drift bound. Not root-caused.
- No cell, bitwise vs the binary before A6b: `-sp -dump_gradient` (file md5 incl. the 16-digit energy line) on benzene, caffeine,
  methyl-a-d-glucopyranoside, water8, hydronium-water2 (+1), urea-water197, peo201-water500 at 1 and 4 threads 14/14 identical;
  `-md` gfnff 200 fs seed 7, 4 threads, caffeine and water8: energy table (without the wall-clock column) and trajectory (without
  the version string) identical.
- `ctest -j1 -R 'gfnff|frag|eeq|simplemd|pbc|ewald|react|restart|md_'`: 125 tests, 0 failed, 1 skipped (GPU), one run; the two
  seed flakes of A6 did not occur, so their tests were left unchanged. `scripts/check_docs.py`: the 9 known errors only.
- Not tested: `-gpu` with a cell, RATTLE with a cell, modes 2/3, a restart file with a cell (restart keeps `m_dof = 3N` as
  before, independent of A6b), CG with a cell. Human production testing pending.

## S0 GPU safety net before the periodic GPU port: plugin ABI version, wrapper refuses a cell (2026-11-03, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `89945e40`). Stage S0 of the GPU port plan (G0, scratchpad). Problem: A6 put
`ComputationalMethod::hasPeriodicEnergy()` in the middle of the vtable, so every GPU plugin built before A6 calls wrong
slots in a newer executable, and the loader checked only the Eigen alignment.

- **ABI version**: `CURCUMA_GPU_PLUGIN_ABI_VERSION` 1 in `src/core/energy_calculators/gpu_plugin_abi.h` (rules for raising
  it are in the header); CUDA, ROCm and Vulkan entry files export `curcuma_<backend>_abi_version()`; the loader checks it
  right after `dlopen` (`gpu_plugin::checkPluginCompatibility`), before any plugin object exists. Missing symbol or other
  value: warning naming path, found and expected version and "rebuild the GPU plugin from the same source tree as the
  curcuma executable", counted fallback (`GPU plugin built from another source state (ABI/variant check)`), CPU;
  `-gpu_strict true` exits with 3. Details: [GPU_PLUGIN_STARTUP.md](GPU_PLUGIN_STARTUP.md).
- **Wrapper** `GFNFFGpuMethodImpl<Backend>` (CUDA and ROCm): `hasPeriodicEnergy()` explicitly false (until S4);
  `setMolecule()` with `m_has_pbc` returns false with an error message, as a second line behind
  `MethodFactory::periodicCellNeedsCpu()`; the multi-GPU Coulomb/EEQ split is cleared with a cell (third line; not
  reachable today because `setMolecule()` refuses first). All batch paths (`GpuDevicePool` workers) go through the same
  `setMolecule()`.
- **Native gfn1/gfn2 on the GPU** (`XtbGpuComputationalMethod`, `xtb_gpu_adapter.h`): no cell handling, like the CPU
  `NativeXtbMethod` (reads no `m_has_pbc`, TODO "Periodisch, kleinere Punkte"); unchanged, a cell is silently ignored on
  both. Vulkan has no GFN-FF.
- **Accessor for S3**: `GFNFF::activeEwaldParameters()` (public, const, returns a copy of `m_ewald_params`; zeros without
  Ewald). No caller yet.

Tests:
- `gpu_plugin_abi` (new, no GPU, stand-in plugins from `test_cases/unit/gpu_plugin_abi_stub/`): without the symbol and
  with version 2 refused with path, found/expected version and the rebuild hint; version 1 accepted; Eigen alignment 8
  refused; absent plugin is no rejection; 23/23 checks. `gpu_plugin_abi_cli` (new): stand-in `libcurcuma_cuda.so`
  without version through `curcuma -sp -gpu cuda`: message, counted fallback, energy equal to the CPU run
  (-2.618188272291070 Eh, water8); `-gpu_strict true` exit 3; right version: loaded, then "missing entry point" fallback.
- CUDA (RTX 5080, CUDA 13.3, `-DCMAKE_CUDA_ARCHITECTURES=120`, fresh build directory outside the repo, `make -j6 curcuma
  curcuma_cuda` 448 s, incremental 176 s): water8 without a cell GPU -2.618188272291102 vs CPU -2.618188272291070 Eh
  (dE 3.2e-14), max |dG| 5.6e-15 Eh/A, two GPU runs md5-identical; with `unitcell 200 200 200`: warning, 1 counted
  fallback, result md5-identical to the CPU run with the cell, `-gpu_strict true` exit 3; the real plugin with the
  version symbol renamed in the binary: refused, CPU energy, strict exit 3; `gfnff_gpu_pbc_refuse` (new, CUDA only:
  wrapper constructed directly, guard bypassed): cell refused with the message, no cell sets up, `hasPeriodicEnergy()`
  false; `gfnff_gpu_builtin_molecules`, `pbc_gpu_guard`, `gpu_plugin_abi(_cli)` green in the CUDA build.
- No cell, CPU: `-sp -dump_gradient` 7 structures x 1/4 threads 14/14 md5-identical to the binary before S0;
  `ctest -j1 -R 'gfnff|frag|eeq|simplemd|pbc|ewald|react|restart|md_|gpu|plugin'` 129 tests, 0 failed, 2 skipped (GPU
  tests without plugin, as before); the A6 seed flakes did not occur, so they were left unchanged.
- Not tested: ROCm (no SDK here; the entry file and the wrapper change are compiled only in the CUDA build, the ROCm
  entry file not at all), Vulkan plugin build, Windows, multi-GPU split with a cell (one GPU), an old CUDA plugin
  actually built from a pre-A6 tree (the renamed-symbol copy stands in for it), GPU bitwise against the plugin before S0.
  `libcurcuma_cuda_mgpu.so` (loaded by the CUDA plugin) has no version. Human production testing pending.

## S1 local GFN-FF terms with a cell on the CUDA device (2026-11-04, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `77eac176`). Stage S1 of the GPU plan (G0): every local term, the device pair builders, the
device H-bond detection and the D4 tail apply the minimum image of an orthorhombic cell on the device. **Released
only behind `-gfnff.gpu_periodic true` (PARAM, default false, fallback in `GFNFF::GFNFF(const json&)`) and only with
`coulomb_ewald false`**: the device Coulomb is still the minimum-image pair sum of April 2026 (not Ewald) and the device
EEQ non-periodic (S3/S4). `MethodFactory::periodicCellNeedsCpu()` keeps the CPU for `gpu_periodic false` (message as
before) and for `gpu_periodic true` with `coulomb_ewald true` ("periodic Coulomb/EEQ on the GPU is not available yet
(S3/S4)", counted fallback); the wrapper's `setMolecule()` refuses the same cases and every non-CUDA backend
(`Backend::has_periodic_local_terms`, true for CUDA, false for ROCm). `hasPeriodicEnergy()` stays false, so SimpleMD
keeps the A6b vacuum treatment (no cell wrap) with the GPU.

- **Kernels (`cuda/gfnff_kernels.cu`)**: `micPlace()` / `micChainOfFour()` are the arithmetic of the CPU
  `placeMinimumImage()` / `imageNear()` (ref + mic(x - ref)); without a cell they copy the position, so the non-periodic
  kernels keep their expressions. `k_angles` (i, k next to j; the angle and its gradient were still plain), `k_dihedrals`
  and `k_storsions` (chain of four), `k_inversions` (three neighbours next to the centre), `k_hbonds` (frame on H: A, B
  next to H, neighbours of B and the acceptor parent C next to B, D next to C for the D-B-C-H torsions, lone-pair vector
  from the frame), `k_xbonds` (frame on X), `k_check_displacement` (minimum-image displacement, so a wrap triggers no
  topology or CN-list rebuild), `k_disp_pairs_count/_build` (minimum-image distance, per-pair `r_cut` =
  `GFNFF::dispersionRadiiBohr().eval` instead of the literal 50, build radius `.build`), `repPairWanted()` (minimum
  image), `d_hb_dist2()` / `d_hb_case1_strength()` (minimum image, kernel frame as `estimateHBStrengthCase1`); the FP32
  pre-filters of `k_hb_detect` / `d_hb_strength_decision_f` carry no image and are skipped with a cell (the double test
  decides, as on the host). The FP32 `k_xbonds_mixed` has no frame: with a cell the FP64 kernel runs. Host side:
  `prepareHBDetectInputs()` takes the carbene angle as minimum image (A4b said "still plain").
- **D4 tail (A3c) on the device**: `k_d4_tail_W` (one thread per element slot and reference, atoms in index order:
  deterministic) and `k_d4_tail_atoms` (per atom, the CPU loop order; energy per atom in fixed point into the new slot
  `E_DISP_TAIL`, dE/dCN by `fix_add` into `d_dEdcn` before the CN chain rule), launched in
  `launchChargeDependentAndFinish()` on the main stream after the Gaussian weights of the step. Tables from
  `FFWorkspaceGPU::setDispersionTail()` (zeta from `D4ParameterGenerator::atomZeta()`, radius and volume from the host
  `GFNFF`), set in `initGPUWorkspace()`. Differs from the CPU in the summation order of S6/S8 only.
- **Cell upload**: a `cudaMemcpyToSymbol` through the `extern __constant__` declarations of `ff_workspace_gpu.cu` failed
  silently (the error surfaced as `cudaErrorInvalidDevice` in the next checked call of the set-up, and `d_has_pbc` stayed
  0): the symbols are now written by `gfnffUploadCell()` in `gfnff_kernels.cu` with the status checked. The symbols are
  one state per device for the process, so `uploadCellToDevice()` sets them (also "no cell") before the set-up pair
  build and at every step; unchanged values are not re-copied.
- **Parity fixes of diagnostic switches found on the way**: `coulomb false` on the CPU kept the self term and TERM 1b
  (water8: +0.706 Eh, gradient 4.5e-3 Eh/A against the GPU, which skips both); now both are off with the switch
  (`FFWorkspace::postProcess`). `hbond false` on the GPU left the XB kernel running (the CPU switches HB and XB
  together); now gated alike. Still asymmetric: `hbond false` on the GPU keeps the HB-CN modulation of the X-H bonds,
  the CPU drops it with the detection (urea197: 6.5e-3 Eh/A; TODO.md). Defaults are untouched.
- **ABI**: `GFNFF` gained `m_gpu_periodic`, so `CURCUMA_GPU_PLUGIN_ABI_VERSION` is 2 (plugins must be rebuilt).

Tests (`gfnff_gpu_pbc_terms`, CUDA; returns 77 = skipped without CUDA; `coulomb false`, `coulomb_ewald false`,
`gpu_periodic true`, `hb_update_force_every 1`, 44 Bohr cube; thresholds fixed beforehand: 1e-11 Eh per term, 1e-12
Eh/A for the gradient): water dimer, caffeine, hydronium-water2 (+1, `frag_charge_model reference`), acetic acid dimer,
acetic dimer + F-I...pyridine composite, 40 and 60 random waters (the 60-water box also with
`nb_cell_list_min_atoms 0`: device repulsion build, cell-list HB branch); states cell centre, shifted by (L,-L,2L),
fragments by different lattice vectors, corner, torn per atom (also initialised torn, A4c), 2x1x1 supercell.
- GPU vs CPU over all 8 systems x 5 states (+ dispersion off / repulsion off at centre and torn): worst term 2.2e-13 Eh
  (disp_tail; bond 7.1e-14, angle 4.2e-14, bonded repulsion 1.0e-13, non-bonded 8.6e-14, dispersion 7.8e-15, HB 3.1e-14,
  XB 2.6e-15, batm 1.9e-14, dihedral 3.1e-15, inversion 7.7e-16; no structure has sTors), worst gradient 1.5e-13 Eh/A;
  isolated dispersion / repulsion gradients <= 4.3e-15; HB triple counts equal (caffeine 358, composite 550, 60 waters 1468).
- GPU invariance between the states: terms <= 1e-11, gradient <= 1e-12 (BATM has no minimum image on either side, A3b,
  so its term and, where present, the gradient are left out of the invariance check; GPU == CPU holds for them).
- 2x1x1 supercell (intact molecules, neutral systems): every term = 2 x cell within 5.7e-14 Eh, per-atom gradient
  3.3e-14 Eh/A, HB counts doubled. (Hydronium is excluded: the reference charge rule puts the net charge on the first fragment.)
- FD of the GPU energy on the torn geometries (8 atoms nearest a face, h 2e-4 / 1e-4 A): caffeine 8.5e-8 -> 2.1e-8,
  Richardson 1.9e-10; composite 3.1e-8 -> 8.0e-9, 2.2e-10; 40 waters (4 atoms) 9.8e-8 -> 2.4e-8, 1.7e-9 Eh/A.
- Verlet: a wrap of one molecule by +L rebuilds no D4 / repulsion list (device-built, cell list), dE 3.6e-15, dG 2.9e-15;
  a 1.5 A move rebuilds. Two device instances on the torn 60-water box and the composite: energy and gradient bitwise.
- CLI (water8 in a 25 A cube): `-gpu cuda` alone -> CPU as before; with `gpu_periodic` alone -> S3/S4 fallback, result
  identical to the CPU run; with `coulomb_ewald false`: `coulomb false` GPU vs CPU dE 5.1e-14 Eh, full GPU run twice and
  once next to a concurrent GPU MD: md5-identical. `gfnff_gpu_pbc_refuse` extended ((d) refusal with the S3/S4 text,
  (e) accepted with `coulomb_ewald false`, 200 A cube = no cell), `pbc_gpu_guard` extended (g, h).
- CUDA build: `ctest -R 'gpu|plugin|pbc'` 168 of 168 pass (3 skipped: split, scope-flags, deltaml). The 17
  `gfnff_gpu_vs_cpu_*` had failed with "File not found: ../../test_cases/structures/..." because this build directory
  lies outside the repository (their reference JSON files carry repository-relative paths); `test_gfnff_gpu` now falls
  back to `CURCUMA_SOURCE_DIR` (CMake) when the relative form does not exist (caffeine 2.3e-14 Eh, 1.5e-13 Eh/Bohr).
  CPU build: `ctest -j1 -R 'gfnff|frag|eeq|simplemd|pbc|ewald|react|restart|md_|gpu|plugin'` 130 of 130 (2 skipped);
  `pbc_gpu_guard` has two more checks.
- Without a cell: CPU `-sp -dump_gradient` 7 structures x 1/4 threads 14/14 md5-identical to the binary before S1. GPU
  (plugin + executable saved before S1, each structure run twice): energies (18 digits) and gradient dumps md5-identical
  14/14. **Bisection of a first attempt** (the kernels read the frame through local arrays, plain copies without a
  cell): urea197 and peo201-water500 differed on one N atom by one fixed-point quantum (4.2e-16 Eh/A), with HB,
  dispersion, repulsion and Coulomb switched off as well, so in the bonded kernels. `cuobjdump -sass` of the two plugins
  showed every kernel not touched by S1 byte-identical and the rewritten bonded kernels with another instruction mix
  (k_inversions DFMA 248 -> 284, DMUL 347 -> 365, DADD 78 -> 96; k_dihedrals 262/360/99 -> 298/378/125; k_angles
  198/159/54 -> 194/157/66; k_storsions 204/197/45 -> 240/215/72): the compiler contracted the array-based expressions
  into other fused multiply-adds, which rounds differently. Fix: the six local-term kernels are templates on PBC; the
  `false` instantiation keeps the former statements verbatim (`if constexpr` at each divergence point, launch through
  `launch_k_*()` in `gfnff_kernels.cu`, because nvcc gives the host stubs of template kernels hidden visibility), the
  `true` instantiation builds the frames. The SASS of `k_angles/k_dihedrals/k_inversions/k_storsions/k_xbonds/k_hbonds<false>`
  is instruction-for-instruction identical to the pre-S1 kernels (6626/12498/11762/7330/7250/33714 instructions); the
  kernels that still differ in SASS are the ones whose cell-less behaviour is a decision, not a sum (`k_hb_detect`,
  `k_hb_ab_pairs`, `k_rep_pairs_*`, `k_disp_pairs_*`, `k_check_displacement`).
- Not tested: timing (none measured), MD with a cell on the GPU (SimpleMD treats the GPU run as non-periodic, A6b),
  ROCm (refuses; kernels untouched), Vulkan, multi-GPU split (one card), charged multi-fragment systems under the
  ensemble charge model (CPU fallback as before), sTors and extra torsions under a cell (no structure), case-4 H-bonds,
  `topology_mode auto` / react on the device with a cell, systems with more than one element class in the tail beyond
  H/C/N/O/F/I. Human production testing pending.
