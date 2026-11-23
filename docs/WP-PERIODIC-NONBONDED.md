# WP: Periodic non-bonded interactions (minimum image)

**Status (2026-11-05, after N1): native CPU GFN-FF with an orthorhombic unit cell is periodic** (MIC in CN, D4,
repulsion, bonded, HB/XB, topology; Ewald Coulomb and EEQ with `coulomb_ewald`, A1-B3, P1) **and has an analytic
atomic virial / instantaneous pressure** (N1, `-md -compute_pressure true`). With such a cell
`wall_potential = pbc` in a box is the pure cell wrap, without clearance test or reflection (A6). The stopgap
described below **stays** for everything that has no periodic energy: no cell, the spherical container, every other
method (gfn1/gfn2, external xtb/tblite/Ulysses, UFF, QMDFF, CG) and GFN-FF with `coulomb_ewald false`; for those a
cell is not wrapped and the COM removal is the vacuum one (A6b, matrix there). Open: GPU paths with a cell (S2 hybrid
opt-in, no virial there), triclinic cells, the barostat (N2). **Since Nov 14, 2026 the native gfn1/gfn2 with an
orthorhombic cell are periodic energies too (energy and gradient, no virial): [WP-PERIODIC-XTB.md](WP-PERIODIC-XTB.md);
the A6b lists naming gfn1/gfn2 as non-periodic are superseded there.**

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

- **Virial (derived here; implemented and FD-tested in N1, 2026-11-05)**: diagonal strain of an orthorhombic box, W_aa = -dE/d eps_aa:
  real `-sum_pairs q_i q_j K'(r) d_a^2 / r`; reciprocal per k `E_k [1 - 2 k_a^2 (1/k^2 + 1/(4a^2))]`; self 0;
  background `E_bg`. The full tensor and its check are in section N1.
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

## S2 CUDA periodic hybrid: local terms on the device, periodic EEQ and Ewald Coulomb on the host (2026-11-04, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `8ee4f638`). Stage S2 = variant A' of the GPU plan (G0): the first release of `-gpu cuda` with a
cell and the default `coulomb_ewald true`, behind `-gfnff.gpu_periodic true` (PARAM default false). The result is the
periodic energy of the CPU path (`hasPeriodicEnergy()` true, so SimpleMD wraps into the cell, A6/D1). Device Ewald
(S3) and device periodic EEQ (S4) are not part of it; `coulomb_ewald false` keeps the S1 behaviour (minimum-image
Coulomb on the device, non-periodic EEQ, `hasPeriodicEnergy()` false); ROCm/Vulkan refuse as before.

### Architecture of the seam

- **Charges**: the device computes the CN (S1 minimum image, `k_cn_compute`); the wrapper takes the existing
  `eeqCpuFragmentRoute()` (CN download, `GFNFF::prepareCNAndEEQ(gradient, gpu_only, &cn, skip_eeq=false)`, charges
  up): with a cell and `coulomb_ewald` the host `EEQSolver` is already periodic (`calculateFinalChargesPeriodic`,
  B3/P1: dense below `eeq_ewald_dense_max_atoms`, block-preconditioned pcg with aspc warm start and `eeq_ewald_tol`
  above). Nothing new on this path; the device EEQ solvers are never entered with a cell (also not the
  first-step `isFragmentTopoValid()` branch).
- **Coulomb**: `FFWorkspaceGPU::setCoulombPairOnDevice(false)` skips the device pair kernels (`k_coulomb_tiles` /
  `k_coulomb_implicit` / `k_coulomb`) only; `k_coulomb_postprocess` (self term `-q chi(CN) + 1/2 q^2 (gam +
  sqrt(2/pi)/sqrt(alp))` = TERM 2+3, and TERM 1b `qtmp` into the CN chain rule) stays on the device, because the CPU
  keeps exactly these in `postProcess()` while `calcCoulombEwald()` holds only the pair part (`ff_workspace.cpp`
  `postProcess`, `ff_workspace_gfnff.cpp:1173`): the split is the same on both sides, nothing is counted twice or
  missed. The host side is the new `FFWorkspace::evaluateEwaldCoulomb(geometry_bohr, charges, gradient, &grad, &terms)`:
  the B2 code (`calcCoulombEwald` per partition over the pool, reduced in partition order, `selfAndBackground`), on
  the P1 geometry cache that the EEQ solve of the same step just built (same positions, widths and parameters bit for
  bit, so the pair list, kernel and phase table are reused, no second build). Energy in Eh and gradient in Eh/Bohr are
  added to the device result in the wrapper (`m_last_energy`, `m_cached_gradient` before the `/au` of the accessors).
  The Ewald pair part is CN-independent at fixed charges and the charges are stationary (B3), so no dE/dCN or dE/dq
  term crosses the seam.
- **Gating**: `MethodFactory::periodicCellNeedsCpu()` keeps the GPU for every `gpu_periodic true` run (the S1 "S3/S4"
  fallback is gone); the wrapper's `setMolecule()` refuses a cell without the opt-in and on a backend without
  `has_periodic_local_terms`, checks after `InitialiseMolecule()` that the host periodic EEQ and the workspace Ewald
  are set up, and switches the device pair Coulomb off (set in `initGPUWorkspace()`, so a react-topology rebuild of
  the workspace keeps it off; react with a cell on the device is otherwise untested). `hasPeriodicEnergy()` = hybrid
  active. No change to
  `ComputationalMethod`, the plugin entry points or a shared class layout (`FFWorkspace` gained a method, no member;
  the new flag lives in the plugin-only `FFWorkspaceGPU`): ABI version stays 2.
- Diagnostics: `CURCUMA_GFNFF_STEPPHASE=1` prints `ewald=` (host Ewald ms) next to `eeq=`; `getLastPrepTiming()`
  has `gpu_host_ewald`; `getEnergyDecomposition()` adds `Coulomb_host_ewald(_real/_reciprocal)`; the verbosity-2 report
  shows the Coulomb term as device self term + host Ewald.

### Acceptance (thresholds fixed beforehand; ctest `gfnff_gpu_pbc_hybrid`, CUDA only, 77 = skipped without CUDA)

Boxes as in `pbc_eeq_ewald` (rigid random packings): 50 H2O L 44 (dense), H3O+(H2O)2 + 18 H2O Q = +1 L 44 (dense,
`frag_charge_model reference` on both sides, else the wrapper asks for the CPU), 970 H2O L 58 (2910 atoms, pcg, 4
threads); each at the input geometry (both sides reuse the charges of the set-up solve) and at a geometry with every
atom moved by up to 0.005 A (per-step solve, device CN on the GPU).

| check | result |
|---|---|
| (i) GPU hybrid vs CPU, default tolerance | charges <= 8.9e-16 e (bitwise at the input geometry), energy 50 H2O 3.2e-13 / 7.8e-14 Eh, Q = +1 7.0e-14 / 0, 970 H2O 2.7e-12 / 2.4e-12 Eh (limits 1e-12 e, 1e-9 Eh) |
| (i) `eeq_ewald_tol` 1e-14 on both sides | energy <= 2.6e-12 Eh, gradient <= 4.1e-14 Eh/A, host Ewald pair part <= 3.6e-15 Eh (limits 1e-11 Eh, 1e-12 Eh/A) |
| (ii) FD gate, hybrid total gradient vs central FD of the total energy with re-solved charges, h 2e-4 / 1e-4 A, 6 atoms nearest a face | 20 H2O: 8.3e-8 -> 2.1e-8 (ratio 3.95), Richardson 6.0e-10 Eh/A; Q = +1: 7.9e-8 -> 2.0e-8, Richardson 7.9e-10 (limit 1e-8; CPU B3 (e): 5.8e-11..1.1e-10) |
| (iii) invariances of the hybrid (20 H2O, Q = +1): translation (L,-L,2L), atom-wise wrap (on the same instance and initialised torn), 2x1x1 supercell (neutral) | charges <= 7.8e-15 e, energy <= 3.2e-14 Eh (supercell = 2 x cell), per-atom gradient <= 1.9e-14 Eh/A |
| (iv) two hybrid instances, two steps each (second from the aspc warm start), 1 and 4 threads | energy, gradient and charges bit for bit |
| CLI, 50 H2O L 44, 4 threads | hybrid vs CPU periodic dE 1.4e-13 Eh, max dG 1.7e-14 Eh/A; two hybrid runs and one next to a concurrent GPU MD: `-dump_gradient` md5 identical; `-gpu cuda` without the opt-in: CPU fallback, md5 equal to the CPU run; `coulomb_ewald false` (S1): dE 1.0e-3 Eh (minimum-image Coulomb, as expected) |
| (vi) without a cell, plugin and executable before S2 vs after | GPU `-sp -dump_gradient` 7 structures x 2 runs 14/14 md5 identical; CPU 7 structures x 1/4 threads 14/14 md5 identical (no kernel touched) |
| (vii) ctest | CUDA build `-R 'gpu|plugin|pbc'` 169 of 169 pass (3 skipped as in S1: split, scope-flags, deltaml), before and after the HB-CN fix below; CPU `release/` `-j1 -R 'gfnff|frag|eeq|simplemd|pbc|ewald|react|restart|md_|gpu|plugin'` 131 of 131 (4 skipped: GPU) |
| (viii) guards | `pbc_gpu_guard` (g) and `gfnff_gpu_pbc_refuse` (d) changed on purpose: `gpu_periodic` with the default `coulomb_ewald` is accepted (hybrid, `hasPeriodicEnergy()` true) instead of the S1 "S3/S4" refusal; the default (no opt-in) and ROCm paths are unchanged |

(v) NVE, 60 H2O in L 44 Bohr (the D1 box60), dt 0.5 fs, CSVR 300 K 1 ps then NVE 2 ps from the edited restart file
(D1 protocol), 4 threads, n = 2 seeds each (drift = linear fit over the last 1.5 ps):

| run | drift (Eh/ps) | E_tot max-min (Eh) | E_tot std (Eh) | E(end)-E(0) (Eh) | <T> (K) |
|---|---|---|---|---|---|
| CPU periodic, seed 7 / 11 | -7.9e-5 / -1.1e-4 | 4.5e-4 / 5.2e-4 | 8.2e-5 / 9.4e-5 | -1.1e-4 / -1.1e-4 | 372 / 383 |
| GPU hybrid, seed 7 / 11 (after the HB-CN fix below) | -3.9e-5 / -7.5e-5 | 4.1e-4 / 5.2e-4 | 8.3e-5 / 9.5e-5 | -5.5e-5 / -1.1e-4 | 371 / 378 |
| GPU hybrid, seed 7 / 11 (before the fix, kept as evidence) | -2.7e-5 / -1.1e-4 | 1.8e-3 / 2.2e-3 | 1.1e-4 / 1.5e-4 | -1.5e-3 / -2.1e-3 | 377 / 385 |

The hybrid conserves energy at the CPU level (D1: -5e-5..-8.5e-5 Eh/ps with this protocol); the CSVR phases of the CPU
and GPU runs start identically and separate chaotically (end states differ by up to 0.16 A/fs in the velocities after
1 ps), so the NVE phases are n = 2 independent samples per side, not the same trajectory. The hybrid runs wrapped into
the cell (`hasPeriodicEnergy()` true, A6/D1 path) with no GPU fallback.

### Time per MD step (ms, CLI `-md`, dt 0.5 fs, 4 host threads, difference of 40- and 20-step runs, 3 repetitions; RTX 5080, nothing else running)

| box | CPU periodic | GPU hybrid | gain | hybrid `calculateEnergy()` median (ms): total / host EEQ / host Ewald / device finish |
|---|---|---|---|---|
| 970 H2O, 2910 atoms, L 58 | 113.8 / 114.9 / 116.3 | 76.0 / 78.7 / 80.0 | 1.46x | 47.2 / 42.2 (pcg) / 4.0 / 0.5 |
| 1916 H2O, 5748 atoms, L 72.8 | 324.5 / 324.9 / 323.7 | 240.2 / 240.4 / 243.6 | 1.35x | 121.2 / 107.1 / 11.1 / 1.0 |

The G0 estimate (~1.5x) holds, for the reason G0 gave: the step is the host EEQ (89 % of the energy call), the host
Ewald pair term (8 %) and the device work that is hidden behind the CN download and the host solve. The CPU numbers
here are the CLI step (P1 measured 114 ms/step CLI at 2910 atoms, 77.5 ms for `Calculation(true)` in the harness);
the hybrid energy call is 47 / 121 ms, the CLI step 78 / 241 ms, so SimpleMD spends 31 ms (2910) and 120 ms (5748)
per step outside the energy call, growing like ~N^1.7 (P1 noted 37 ms at 2910 atoms, not located; with the hybrid it is
now the second-largest item and the first thing to look at before S3/S4).

### Finding on the way (GPU path, not S2-specific, fixed)

The first GPU NVE runs (restart protocol) lost 1.6e-3 / 1.8e-3 Eh between step 0 and step 1 while the CPU runs did not
(E_tot max-min 1.8e-3 / 2.2e-3 against 4.5e-4 / 5.2e-4 Eh, drift fits of the last 1.5 ps fine). Term diff of a two-step
restart (verbosity 2): only the bond term, +1.56e-3 Eh at step 0 (GPU -16.0736371, CPU -16.0752005 Eh), equal from step
1; a fresh `-sp` on the restart geometry agrees to 1e-13 Eh; the same offset with `coulomb_ewald false` (S1 path), so
not the seam. Cause: the wrapper computed the per-bond HB coordination numbers (`computeBondHBCN`, Step 2b) BEFORE the
HB/XB refresh of the step (Step 2c, which uploads the new HB alpha pairs and bond metadata); on every refresh step the
X-H bond modulation used the previous list (the CPU computes it in `calculate()` after `updateHBXBIfNeeded()`). On a
restart the first geometry is 1 ps away from the set-up lists, hence the visible step. Fixed by moving Step 2b behind
Step 2c; a single point is unchanged (no refresh at step 0: no-cell 14/14 md5 identical, CUDA ctest 169/169 before
and after), GPU MD steps on which a refresh fires now use the current list, as the CPU does.

### Not done / not tested

- Not part of S2: Ewald and periodic EEQ on the device (S3/S4), ROCm/Vulkan (refuse), the multi-GPU split (off with a
  cell), a Verlet skin for the host Ewald pair list, the SimpleMD per-step overhead above.
- Not tested: charged boxes under the ensemble charge model (CPU fallback as before), `topology_mode auto`/react with a
  cell on the device, ALPB with a cell (refused by B3), heavy-element boxes near the `r_g < L/2` limit, more than 5748
  atoms, `ewald_kernel_cache_mb` fallback at production size, an energy-only MD (`calculateEnergy(false)` is covered
  only by the FD gate), NPT. Human production testing pending.

## M1 SimpleMD step profile with a cell; RATTLE with a cell (2026-11-05, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `5adb1ca1`). Locates the per-step SimpleMD cost outside the energy call that S2 measured (31 / 120 ms at
2910 / 5748 atoms) and checks RATTLE with a cell.

### Profile method

`CURCUMA_MD_STEPPHASE=1` (new, environment only, read once in `LoadControlJson()`; off = one bool test per item, no clock read,
no number changes) prints one `MDPHASE` line per step on stderr: `com` (COM/rotation removal), `kick1` (first half kick, drift,
NaN checks), `geom` (`EnergyCalculator::updateGeometry`), `energy` (`CalculateEnergy(true)`), `grad` (gradient copy, NaN check),
`bias` (walls, external potentials, Rg bias, restraints), `kick2` (second half kick), `thermo` (thermostat + `EKin()`), `wrap`
(cell wrap without the fragment search), `frag` (`setGeometry` + `GetFragments()` inside the wrap), `container`
(`wrapIntoContainer()`), `after` (averages, aborts, output, restart, status line, CheckStop). Boxes of P1/S2 (970 H2O L 58 Bohr,
1916 H2O L 72.8 Bohr), the same coordinates without a cell as vacuum reference, CLI `-md`, dt 0.5 fs, 4 threads, seed 7, 40
steps, median over steps 6..40, median of 3 repetitions; GPU = S2 hybrid on the RTX 5080 (`-gpu cuda -gfnff.gpu_periodic true`).

### Time per step (ms; all items not listed are <= 0.03 ms)

| run | energy | frag | rest of the step | step total |
|---|---|---|---|---|
| 2910 cell CPU, before | 77.6 | **29.5** | 0.1 | 107.8 |
| 2910 cell CPU, after | 79.4 | 1.2 | 0.1 | 80.6 |
| 2910 no cell CPU (before = after) | 45.5 | 0 | 0.05 | 45.5 |
| 2910 cell GPU hybrid, after | 48.8 | 1.2 | 0.1 | 50.1 |
| 5748 cell CPU, before | 198.0 | **113.5** | 0.1 | 313.9 |
| 5748 cell CPU, after | 200.4 | 2.3 | 0.1 | 202.9 |
| 5748 no cell CPU (before = after) | 147.6 | 0 | 0.05 | 147.7 |
| 5748 cell GPU hybrid, after | 122.9 | 2.3 | 0.1 | 125.5 |

The whole cell-specific SimpleMD cost was `Molecule::GetFragments()` in `applyPeriodicBoundaryConditions()`: in a liquid some
molecule always straddles a face, so the O(N) gate never skipped it, and the search was an all-pairs breadth-first scan over a
`std::vector` with `erase` (O(N^2) minimum-image distances, 113.5 / 29.5 = 3.85 for 1.98x the atoms). Verlet, thermostat, COM
removal (O(N), not due in 40 steps), copies, output: each <= 0.03 ms. Wall clock of the CLI (difference of 40- and 20-step
runs, median of 3; GPU "before" = executable before M1 with the current plugin):

| box | CPU before / after | GPU hybrid before / after |
|---|---|---|
| 2910 | 113.9 / 85.9 | 80.8 / 53.7 |
| 5748 | 324.6 / 212.0 | 246.7 / 134.6 |

The wall clock equals the mean (not the median) of `total` over the same steps (86.5 / 214.0 / 53.5 / 134.7 ms): the 3-9 ms
between table and wall clock are steps with a slower energy call (HB/XB list refresh), not SimpleMD. SimpleMD outside the energy
call is now 1.3 / 2.4 ms. The GPU hybrid gain against the CPU rises from 1.41x / 1.32x to 1.60x / 1.58x.

### What was changed (no number moves)

- `Molecule::GetFragments()` takes from `Molecule::fragmentCellListMinAtoms()` atoms on (256, setter for tests) the candidate
  partners from a cell grid (`Molecule::fragmentsByCellList()`): wrapped fractional coordinates with a cell, Cartesian bins without,
  bin edge >= 1.001 x 2 max(r_cov) x scaling, bin count capped at 8N+64 (larger bins), fewer than 3 bins on a periodic axis =
  the old scan. Same seed order, same member order, the accepted partners of one member sorted by index, the same
  `CalculateDistance(PBC)(u, v)` call and threshold expression, the same mass sum order (the seed mass is still left out of the
  ordering key), so fragments, their order and `FragmentMass()` are bitwise those of the scan (completeness argument in the code
  comment). It speeds up every `GetFragments()` caller above 256 atoms (wrap, container wrap, topo check, DOF set-up).
- `CURCUMA_MD_STEPPHASE` profile (above). `Molecule::getScaling()` accessor.

### Checks

- ctest `molecule_fragments_celllist` (new, 0.1 s): grid == scan bitwise (fragments, order, masses) on 300-water boxes at liquid
  density with molecules shifted by lattice vectors and single atoms torn across faces (10 seeds, with and without the cell),
  two triclinic cells, a dense random gas of caffeine atoms (merged chains, 3 seeds, with and without the cell), a 5 A cell (scan
  fallback), 300 waters over 5000 A and one atom 3.7e5 A away (bin cap); 970 waters: scan 35.9 ms, grid 1.15 ms.
- No cell, bitwise against the executable before M1: `-sp -dump_gradient` 7 structures x 1/4 threads 14/14 md5 identical;
  gfnff `-md` 200 steps dt 0.5 fs seed 7, 4 threads, caffeine and water8: energy table (13 columns, without the two timing
  columns) and trajectory (without the version line) md5 identical, last frame identical.
- With a cell, bitwise against the executable before M1: `-sp -dump_gradient` of the 2910-atom box and box60 md5 identical;
  `-md` box60 200 steps (180 atoms, old path), box60 with `-rattle 2`, and the 2910-atom box 100 steps (new grid path, 939
  fragments): table, trajectory and last frame md5 identical (so the energies agree exactly, not only to 1e-12 Eh).
- CPU `ctest -R 'gfnff|frag|eeq|simplemd|pbc|ewald|react|restart|md_|gpu|plugin' -j1`: 133 of 133 passed (4 skipped: GPU),
  incl. the two new tests; the A6 seed flakes did not occur, so their tests were left unchanged. CUDA build (`build_cuda_bulk`,
  full `make`) `ctest -R 'gpu|plugin|pbc'`: 170 of 170 passed (3 skipped as in S1/S2: split, scope-flags, deltaml). ABI version
  unchanged (no layout of a plugin-facing class touched; `make` relinked the plugin).

### Not changed, recommended (each changes numbers or behaviour)

- The remaining 1.2 / 2.3 ms are the grid search every step. Re-using the fragments while no atom moved more than a skin since the
  last search would remove it, but makes the wrap depend on a history (a reactive merge within the skin would be seen later).
- `InitConstrainedBonds()` calls `Molecule::DistanceMatrix()`: two dense N x N matrices (2 x 68 MB at 2910 atoms, 2 x 264 MB at
  5748) and an O(N^2) loop at set-up; a neighbour-list bond search would give the same list (not measured, set-up only).

### RATTLE with a cell (CPU periodic GFN-FF, `coulomb_ewald` on)

(a) **Constraint pairs across a wrap.** The constraint list is built from raw distances (`simplemd.cpp:1234` `DistanceMatrix()`,
raw norm, `molecule.cpp:2284` `<=`; reference length `simplemd.cpp:1246` `CalculateDistance`, raw) and `Rattle()` evaluates raw
coordinate differences (`simplemd.cpp:3951`/`3960` 1-2 position stage, `4009`/`4018` 1-3, `4153`/`4159` and `4188`/`4194` velocity
stage). The cell wrap moves whole fragments, and a constrained pair is always in one fragment (its length stays below the
fragment threshold), so the raw differences stay valid: 60 H2O L 44 Bohr, `-rattle 2`, dt 1 fs, 2 ps, 165 / 167 molecule wraps
(seeds 7 / 11, counted every 10 steps): max |d^2 - d0^2| over all O-H pairs and frames 1.02e-4 A^2 (= `rattle_tol_12` 1e-4), equal to
the vacuum run; longest raw O-H 0.9635 A. **Limitation found**: a bond that already crosses a face in the input (atom wrapped alone)
is bonded for GFN-FF and the fragment search (minimum image) but not for RATTLE, so it is silently unconstrained for the whole run
and the DOF count treats it as free. Now detected and warned at set-up (`rattleCrossingBonds()`, verbosity 1); the list is unchanged.
The P1/S2 boxes and box60 have no such bond (checked: O-H count raw = minimum image).

(b) **NVE**, box60 (60 H2O, L 44 Bohr), fresh 300 K start, no thermostat, 2 ps, drift = linear fit of E_tot over the last 1.5 ps,
4 threads, n = 2 seeds (7 / 11):

| run | drift (Eh/ps) | E_tot max-min (Eh) |
|---|---|---|
| Verlet dt 0.5, cell / vacuum | -2.8e-5, -2.1e-5 / -1.7e-5, -3.6e-5 | 1.4e-3, 1.6e-3 / 1.4e-3, 1.6e-3 |
| Verlet dt 1.0, cell / vacuum | -1.1e-4, -1.3e-4 / -7.3e-5, -1.1e-4 | 6.2e-3, 7.2e-3 / 6.2e-3, 7.3e-3 |
| `-rattle 2` (= `-rattle 1` for water: same 120 O-H constraints, bitwise the same run) dt 1.0, cell / vacuum | -2.1e-5, -5.5e-5 / -3.1e-5, +1.3e-5 | 6.2e-4, 4.7e-4 / 4.9e-4, 4.1e-4 |
| same, GPU hybrid, cell | -2.2e-5, -5.9e-5 | 6.1e-4, 4.9e-4 |
| rigid water (`-rattle 1 -rattle_13 true`, atoms reordered H-O-H, 180 constraints) dt 1.0, cell / vacuum | -4.3e-4, -5.2e-4 / -2.6e-4, -2.5e-4 | 1.3e-3, 1.3e-3 / 8.6e-4, 8.7e-4 |
| rigid water, `rattle_tol_12` = `rattle_tol_13` = 1e-7, seed 7, cell / vacuum | +2.6e-6 / +1.1e-5 | 1.6e-4 / 1.2e-4 |

RATTLE on O-H at 1 fs conserves energy as well as Verlet at 0.5 fs, with and without the cell (drift within the n = 2 spread of
the vacuum runs). The rigid-water loss is the default tolerance (`rattle_tol_13` 1e-3 A^2), not the cell: with 1e-7 it drops 100x
(n = 1). **Second limitation found**: the 1-3 list (`simplemd.cpp:1251`, `for (int k = 0; k < j; ++k)` inside `j < i`) only finds
chains i-j-k with i > j > k, so H-O-H with O first in the file gets 0 of 60 angle constraints (`-rattle_13 true` silently does
nothing on the P1/S2 boxes and the registry water); a fix changes every `rattle_13` run (TODO.md).

(c) **DOF**: 3N - N_c - 3 with the cell (translation only, D1) and 3N - N_c - 6 without: box60 417 / 414, rigid water 357; ctest
checks 53 / 50 on 8 waters.

(d) **GPU hybrid** (`-gpu cuda -gfnff.gpu_periodic true -rattle 2`, hybrid confirmed by the `EEQ GPU Phase 2: periodic hybrid` line
and `STEPPHASE`): same constraint result through the wraps (max |d^2 - d0^2| 1.02e-4 A^2, 166 / 172 wraps), DOF 417, NVE in the table
above. RATTLE is host code, so nothing GPU-specific is expected beyond the energy.

ctest `pbc_md_rattle` (new, 0.4 s): 8 registry waters in L 44 Bohr, 6 placed 0.02 A inside a face, gfnff periodic, `-rattle 2`, dt 1
fs, 400 steps: 6 wraps, O-H constraints within 2 x `rattle_tol_12` in raw coordinates in every step, drift -9.7e-6 / -1.4e-5 Eh/ps and
max-min 7.1e-5 / 6.9e-5 Eh (cell / vacuum; gate |drift| < 2e-4, max-min < 5 x vacuum + 2e-5), DOF 53 / 50; rattle 1 + rattle_13:
DOF follows the applied constraints (16 + 0, the defect above); torn input: 15 constraints, 1 crossing bond, DOF 54.

Not tested: triclinic cells with RATTLE (GFN-FF refuses them), RATTLE with `topology_mode react` (refused by the PARAM text),
restart files with RATTLE and a cell, molecules other than water, more than 180 atoms with RATTLE, ROCm. Human production testing
pending.

## N1 atomic virial and instantaneous pressure of the periodic GFN-FF energy (2026-11-05, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `71761fad`), CPU only. Prerequisite for a barostat (N2, not part of this); no cell changes here.

### Definition and formulas

- **Atomic virial** `Xi_ab = -dE/d(eps_ab)` under an affine strain of the cell AND of every atom at fixed fractional
  coordinates (molecules are deformed with the cell, no molecular/COM correction). Isotropic pressure
  `P = (2 K + tr Xi) / (3 V)` with K the kinetic energy of all atoms and V the orthorhombic cell volume; atomic units
  Eh, Bohr, Eh/Bohr^3; `1 Eh/Bohr^3 = 2.94210157e13 Pa = 2.94210157e8 bar` (CODATA 2018, `CurcumaUnit::Pressure`,
  checked in the test to 1e-8). The full 3x3 tensor is accumulated (as cheap as the trace); SimpleMD uses the trace.
- **Rule for every kernel** (`ff_workspace_gfnff.cpp`, `VirialTermRows`, `virialPair`): a term evaluated at the positions
  `x_a^frame` of its local minimum-image frame (A4a/A4b; plain positions without a cell) gives
  `Xi_ab -= sum_a x_{a,a}^frame G_{a,b}` with the gradient rows G the term wrote; the row sum of one term is zero, so the
  frame origin drops out (positions are taken relative to the first atom of the term). Pair terms:
  `Xi_ab -= d_a G_{i,b}` with `d = x_i - x_j` as the kernel formed it (minimum image). Applied to bonds (incl. the
  HB-CN modulation, frame B next to H), angles, torsions (incl. extra torsions), inversions, sTors, D4 pair term, bonded
  and non-bonded repulsion, HB (frame on H: A, B, neighbours of B, C and D as the kernel images them; an atom in two
  roles is counted once), XB (frame on X), BATM and ATM (plain positions, as their energies), plain Coulomb (cluster).
- **CN chain rule** (`FFWorkspace::postProcessVirial`): the gradient `G_i += sum_j v_j dCN_j/dx_i` from the
  `CNDerivStore` pair entries (i, j, c) gives `Xi_ab -= sum_entries v_j d_a c_b` with `d = x_i - x_j` (minimum image, the
  vector the derivative was built with); the diagonal entry is the negative row sum and drops out. `v` is the same
  `dEdcn_combined` as in the gradient, split as the gradient components: bond (dr0/dCN), dispersion (dC6/dCN of the D4
  pair term and the D4 tail), Coulomb (TERM 1b, chi(CN), `-qtmp`).
- **Coulomb by Ewald** (`EwaldSum::realSpaceVirialRange`, `reciprocalVirialRange`, `backgroundVirial`; derivation in
  `ewald.cpp`): strain at fixed fractional coordinates scales every pair vector by (1 + eps), V by (1 + tr eps) and
  k by (1 - eps) (k.r invariant). Real space `Xi_ab -= q_i q_j K'(r) d_a d_b / r` with the stored P1 kernel derivative
  (point kernel minus the Gaussian correction). Reciprocal: `E_k = w_k |S(k)|^2`, `w_k = (4 pi/V) exp(-k^2/4a^2)/k^2`,
  `d(k^2)/d(eps_ab) = -2 k_a k_b`, `dV/d(eps_ab) = V delta_ab`, hence
  `Xi_ab += sum_k E_k [delta_ab - 2 k_a k_b (1/k^2 + 1/(4 a^2))]`, trace `sum_k E_k [1 - k^2/(2 a^2)]` (the B0 note per
  axis is this tensor's diagonal; sign and factors confirmed by the FD gate below). Self term `-(a/sqrt pi) sum q^2`: 0.
  Background `-pi Q^2/(2 a^2 V)` and Gaussian background `(pi Q/V) sum q alp` are ~ 1/V: `Xi_ab = E delta_ab`, trace 3E.
  Charges are fixed (Hellmann-Feynman: the Phase-2 charges are stationary with the constraints, B3), alpha is fixed.
- **D4 tail** (A3c) at fixed absolute r_c: `E_tail ~ 1/V`, so `Xi_ab = E_tail delta_ab` (trace 3 E_tail), plus its dE/dCN
  through the chain rule. Note: with the default `dispersion_cutoff_bohr` 0 the evaluation radius is `nextafter(L_min/2)`
  and would follow a strain; the energy the virial belongs to is the one with r_c held fixed (the test sets an absolute
  cutoff; a barostat must either hold r_c fixed or differentiate r_c(L) as well).
  The reported pressure adds the crossing-pair correction 3 E6 + 5 E8 since V-LRC (section below, `pressure_lrc`).
- **Not in the virial**: the ensemble charge blend (`frag_charge_model ensemble`, geometry-dependent window weights),
  ALPB, ATM (no minimum image): `GFNFF::setWantVirial(true)` refuses them with the reason; without a cell or with
  `coulomb_ewald false` likewise. The GPU hybrid (S2) evaluates the local terms on the device and produces no virial
  (`virialValid()` false).

### Implementation

- Every GFN-FF kernel of `FFWorkspace` is a `template <bool VIRIAL>` (`calc*Impl`); the non-template name dispatches on
  `setWantVirial()` outside the loops. The `false` instantiation is the former kernel statement for statement (`if
  constexpr` blocks discarded, the row recorder is an empty struct there), the `true` one records the accumulator rows
  of the term's atoms before the writes and contracts the difference afterwards. Per-partition `FFVirialComponents`
  (per term, plus `cn_bond/cn_dispersion/cn_coulomb`, `coulomb_bg`, `disp_tail`), reduced in partition order; the Ewald
  virial as separate passes over the same pair/k ranges (the energy/gradient loops of P1 are untouched).
- `GFNFF::setWantVirial()`, `virial()` (3x3, Eh), `virialComponents()`, `virialValid()`, `cellVolumeBohr3()`;
  `FFWorkspace::evaluateEwaldCoulomb()` (hybrid) clears `virialValid()`. Plugin ABI version 3 (GFNFF, FFWorkspace and
  FFAccumulator gained members; CUDA plugins must be rebuilt, not done here).
- **SimpleMD** `-compute_pressure true` (PARAM, default false, no effect without a cell): after each `FastEnergy()`
  `readVirial()` takes tr Xi from `gfnffInstance()`; `AverageQuantities()` forms `P = (2 Ekin + tr Xi)/(3 V)` in bar;
  the status line gains the columns `P`, `<P>` (bar), `V` (A^3), `rho` (g/cm^3, from the atom masses). Refused at
  set-up with RATTLE ("pressure with constraints not implemented"), with `cleanenergy`, with a non-periodic energy; a
  GPU-evaluated energy throws at the first step ("the energy of this step carries no atomic virial"). With
  `-gpu cuda` and a cell without `gpu_periodic` the A3b/S0 guards already run the CPU method, so the pressure works
  (counted fallback as before).

### Tests (`ctest -R pbc_virial`, `test_cases/unit/test_pbc_virial.cpp`)

Protocol fixed before the first run: Ewald `alpha` 0.32, `r_cut` 19, `ewald_tolerance` 1e-14 (k set converged far below
the FD noise, no parameter follows L), `dispersion_cutoff_bohr` 20 (absolute r_c), `hb_update_force_every` 1,
`hb_min_pair_energy_eh` 0, `eeq_ewald_tol` 1e-14 (dense EEQ below 1000 atoms anyway), `frag_charge_model reference`,
L = 44 Bohr. FD: `Xi_FD = -[E(+eps) - E(-eps)]/(2 eps)` with `L -> L(1+eps)`, `x -> x(1+eps)` for every atom, a fresh
instance per strained geometry (topology, lists and charges re-derived), eps = 1e-5 and 2e-5 (Richardson). Criterion
`|Xi_an - Xi_FD(1e-5)| / max(|Xi_FD|, 1e-6 Eh) <= 1e-6` per term and total. Per term means per energy component
(there are no bond/angle/torsion switches in GFN-FF; dispersion/repulsion/hbond/coulomb switches exist but are not
needed: each component energy is a function of the geometry on its own, and the Coulomb component is the EEQ
functional, so its FD with re-solved charges is the fixed-charge derivative). Cutoff events: the number of pairs within
20 Bohr (D4 and repulsion radius), HB and XB triples is counted at every point of the stencil; in all five systems
they were constant (736/826/276/390 pairs, 163/255/358/550 triples, 3 XB), so no cutoff-jump contribution entered
any FD value here (a pair crossing r_c = 20 Bohr under eps = 2e-5 would need a pair within 4e-4 Bohr of it; the D4
pair energy there is ~1e-7 Eh, i.e. a jump of ~5e-3 Eh in Xi_FD if it happened).

| term | 20 H2O (60 atoms) | 60 H2O (180, 4 threads) | H3O+(H2O)2 + 18 H2O, Q=+1 | caffeine | acetic dimer + F-I...pyridine |
|---|---|---|---|---|---|
| bond (+ dr0/dCN, HB-CN) | 4.9e-10 (Xi -2.62) | 4.6e-10 (-7.86) | 5.2e-10 (-2.74) | 9.1e-10 (-2.97) | 4.0e-10 (-3.22) |
| angle | 1.8e-11 (0.022) | 2.0e-11 (0.067) | 4.2e-10 (0.036) | 2.1e-10 (0.062) | 6.9e-11 (0.122) |
| dihedral (+ extra) | 0 | 0 | 0 | 1.7e-9 (8.6e-3) | 1.8e-9 (0.077) |
| inversion | 0 | 0 | 0 | 0 (none) | 1.4e-9 (0.020) |
| sTors | 0 (none in any system) | 0 | 0 | 0 | 0 |
| dispersion pair + tail (+ dC6/dCN) | 4.3e-10 (-4.1e-3) | 2.5e-10 (-0.059) | 4.3e-10 (-5.0e-3) | 2.2e-9 (-0.019) | 2.7e-10 (-0.020) |
| repulsion (bonded + non-bonded) | 1.5e-9 (5.91) | 1.5e-9 (17.8) | 1.5e-9 (6.50) | 1.8e-9 (3.18) | 1.7e-9 (3.38) |
| Coulomb (Ewald pair + bg + self + TERM 1b) | 1.7e-10 (-2.78) | 3.0e-10 (-8.41) | 6.9e-11 (-2.93) | 2.0e-9 (-0.144) | 4.2e-10 (-0.397) |
| hbond | 3.3e-10 (-2.6e-3) | 5.2e-10 (-0.035) | 1.1e-9 (0.046) | 7.1e-11 (-2.6e-7) | 1.8e-9 (0.021) |
| xbond | 0 | 0 | 0 | 0 | 2.7e-9 (-2.2e-3) |
| batm | 0 | 0 | 0 | 1.8e-9 (-0.010) | 1.8e-9 (-0.038) |
| **total** (eps 1e-5 / Richardson) | 1.3e-8 / 8.6e-11 (0.527) | 1.4e-8 / 8.0e-10 (1.479) | 8.9e-9 / 2.7e-10 (0.905) | 8.6e-8 / 3.0e-9 (0.102) | 1.5e-7 / 6.8e-9 (-0.050) |

Relative deviation at eps = 1e-5 (tr Xi in Eh in brackets). The per-term values are 1e-9-level: the FD truncation
of a ~3 Eh term at eps = 1e-5 (Richardson brings every term to <= 1.1e-9). The totals are worse in relative terms
because tr Xi_total is a small difference of bond/repulsion/Coulomb traces of 3-18 Eh (composite: -0.050 from +-3.4),
and Richardson reduces them 20-30x, so the residual is truncation, not a missing term. Every tensor is symmetric to
<= 9.4e-16 Eh. (b) 20 H2O: translation by (L,-L,2L) plus an arbitrary shift 1.2e-13 rel (components 2.6e-14 abs),
atom-wise wrap (5 coordinates moved) 2.5e-14, 2x1x1 supercell = 2 x cell 1.9e-13 rel (components 9.4e-14), 4 vs 1
threads 9.9e-14. (c) with a cell the gradient and the energy with the virial requested are bitwise those without;
without a cell and with `coulomb_ewald false` the request is refused with the reason. (d) 60 H2O random packing in
L = 44 Bohr (0.142 g/cm^3, not equilibrated, not a liquid): tr Xi = 1.479 Eh, virial pressure 1702.5 bar, kinetic
part at 300 K (3N-3 DOF) 587.4 bar, sum 2290 bar; per term: repulsion +20465, bond -9051, Coulomb -9679, angle +77,
dispersion -68, HB -41 bar. 20 H2O: 606 bar (virial part). Numbers reported, not assessed.

### CLI, cost, bitwise checks

- `-md -compute_pressure true` on the D1 box60 (60 H2O, L 44 Bohr, 0.142 g/cm^3, 4 threads, dt 0.5 fs): columns `P <P> V rho`
  = `2375.75 2375.75 12622.966 0.14219` at step 0 (virial part 1769 bar at step 2, kinetic 2 K/(3V) ~ 600 bar at 298 K),
  swinging between -1580 and +1777 bar over the first 10 fs of an unrelaxed packing. The P1/S2 box of 970 H2O in
  L 58 Bohr (2910 atoms, 1.004 g/cm^3, random rigid packing at 1.7 A, NOT a liquid): P = 45361 bar at step 0, -7700
  at 5 fs, +30592 at 10 fs (T 298 -> 728 K: the packing is relaxing). Numbers reported, not assessed.
- Cost: CLI MD step on the 2910-atom box (4 threads, difference of 40- and 20-step runs, 3 repetitions, median):
  87.8 ms without and 92.9 ms with `compute_pressure` (+5.1 ms = the Ewald virial passes (real + reciprocal, the
  reciprocal one re-forms S(k) per k) and the row recorders of the local terms; the pair/k loops of the energy are
  untouched). With `compute_pressure false` the kernels run their non-virial instantiations: zero cost.
- `-gpu cuda` with a cell and `compute_pressure`: the A3b/S0 guards keep the CPU method (no plugin present here:
  "GPU plugin not present, calculation on the CPU", counted fallback), pressure columns as on the CPU (1769 bar at
  step 2, equal to the CPU run). The hybrid path (`gpu_periodic`, plugin present) was NOT run (no CUDA plugin built
  from this tree); by construction it throws at the first step ("the energy of this step carries no atomic virial").
  `-rattle 2 -compute_pressure true`: `[ERROR] compute_pressure: pressure with constraints not implemented (rattle)`,
  exit 1. No cell: a verbosity-1 warning, the run continues with the 15 standard columns.
- **Bitwise** (binary saved before N1): no cell, `-sp -dump_gradient` md5 of benzene, caffeine,
  methyl-a-d-glucopyranoside, water8, hydronium-water2 (+1), urea-water197, peo201-water500 at 1 and 4 threads
  14/14 identical (the `false` template instantiations); with a cell and the default `compute_pressure false`:
  `-sp -dump_gradient` of box60 and the 2910-atom box md5 identical, `-md` 200 steps box60 (also with `-rattle 2`)
  and 100 steps 2910-atom box: energy table, trajectory and last frame md5 identical; no-cell MD caffeine/water8
  200 steps identical.
- `scripts/check_docs.py`: the 9 known errors only. CPU `ctest -j1 -R 'gfnff|frag|eeq|simplemd|pbc|ewald|react|restart|md_|gpu|plugin|molecule'`:
  135 tests, 131 passed, 4 skipped (GPU), 0 failed; `pbc_virial` 1.3 s (dense EEQ, 5 systems).
- **Not tested**: the GPU hybrid with the virial request (see above), sTors under a cell (no structure), HB case 4,
  charged boxes under the ensemble model (refused), the pressure during a long MD (only 20 steps here), triclinic
  cells (refused by GFN-FF), `topology_mode react` with the virial (the react rebuild swaps the workspace lists; the
  virial flag survives `rebuildInteractionLists`, not run), the kinetic contribution with COM removal (K includes all
  atoms, the 3 removed DOF are not corrected). Human production testing pending.

## N2 isotropic barostat for the periodic GFN-FF: Berendsen and stochastic cell rescaling (2026-11-05, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `a26b82da`), CPU only. Goal: NPT runs of water boxes (W1 project, `Versuchsplan.md`). User page:
[NPT.md](NPT.md). Parameters (module `simplemd`): `barostat` none | berendsen | crescale, `pressure` 1.01325 bar, `tau_p`
1000 fs, `compressibility` 4.5e-5 1/bar, `pcoupl_interval` 10. A barostat switches `compute_pressure` on; refused with
RATTLE, `cleanenergy`, `-gpu`, a wall potential, a non-periodic energy or no cell; crescale without a thermostat warns.

### Equations (code comment at `SimpleMD::applyBarostat`)

- Coupling every N = `pcoupl_interval` steps with Dt = N dt from the full instantaneous pressure of the step just finished,
  P = (2 K + tr Xi)/(3 V) (K after the thermostat, Xi at x_{n+1}); cell edges and EVERY atom position scaled by mu about
  the origin (fixed fractional coordinates, the strain of the N1 virial).
- Berendsen (J. Chem. Phys. 81, 3684 (1984)): mu^3 = 1 - (Dt/tau_p) beta (P0 - P); velocities unchanged.
- C-rescale (Bernetti, Bussi, J. Chem. Phys. 153, 114107 (2020), arXiv:2006.09250): d ln V = -(beta/tau_p)(P0 - P) dt +
  sqrt(2 kB T beta/(V tau_p)) dW, one Euler-Maruyama step per coupling (dW = sqrt(Dt) xi, xi ~ N(0,1) from the seeded MD
  generator `m_rng`, drawn after the thermostat's numbers of the same step: reproducible with `-seed`; T = thermostat
  target), positions x exp(eps/3), velocities x exp(-eps/3). **The isotropic equations are those of the work order and
  of memory, NOT checked against the paper's text: the one web search allowed returned the paper (journal, volume,
  article, arXiv id) but not its equations.** Stationary density (Ito noise ~ 1/V): V exp(-(H + P0 V)/kT) in (r, p, eps),
  i.e. exp(-(H + P0 V)/kT) in (r, p, V). Ideal gas with dof thermostatted momentum components: rho(V) ~ V^(dof/3)
  exp(-P0 V/kT), <V> = N_eff kT/P0, <dV^2>/<V>^2 = 1/N_eff with N_eff = dof/3 + 1 (the paper's N + 1 for dof = 3N;
  curcuma's periodic DOF count is 3N - 3, D1, so N_eff = N here).
- Placement in velocity Verlet: kick - drift - force(x_{n+1}) - kick - thermostat - barostat. The forces kept for the
  next first half-kick are those of the unscaled x_{n+1} (as GROMACS does; |mu - 1| ~ 1e-4 per coupling at the
  defaults): no extra energy evaluation, and without a barostat the step is bitwise the N1 step (checks below).

### Cell update (`EnergyCalculator::setUnitCell` -> `GFNFF::updateCell`)

Between two energy calls, in this order: (1) new `PBCUtils::OrthoCell` from the cell in A (orthorhombic check); (2) every
active cutoff against the new L_min/2 (`activeCutoffsBohr()`, which now reports the Ewald r_cut fixed at the set-up
instead of re-deriving one from L); a failure restores the old cell and throws; (3) `m_unit_cell` (A), `m_cell_bohr_ortho`,
`FFWorkspace::setPeriodicCell`; (4) `FFWorkspace::setEwald(m_ewald_params)` and `EEQSolver::updatePeriodicCell(cell)`: new
`EwaldSum` objects (lattice vectors, volume, k vectors with n_max per axis from the FIXED k_max and the new L; alpha and
r_cut fixed), the shared `EwaldGeometryCache` invalidated (pair list, stored kernel, phase table belong to the old
lattice), the EEQ warm-start history (`m_periodic_q_hist`, aspc start values) kept; (5) `D4ParameterGenerator::
setPeriodicCell` with the (frozen) build radius; the D4 tail reads V and the frozen r_c per call; (6) `m_lists_rebuild_
pending`: the D4 list, the distance-filtered repulsion list, the cutoff Coulomb list, the HB/XB detection and the
reactive scan are rebuilt at the next `Calculation()` (the Verlet-skin argument of A3b/A4a bounds pair-distance changes
by atom displacements at a fixed lattice; a changed lattice moves image distances by n (L' - L), which no skin covers, so
the simple and correct choice is the rebuild; cost below); `m_last_eeq_geometry` cleared, `m_hb_reference->needs_update`.
The topology, its Phase-1 charges and the `.topo.json` fingerprint (start cell) are untouched. `cellUpdateCount()`.

### Frozen cutoffs (`GFNFF::freezeCutoffsForCellChanges`, called by SimpleMD at the barostat set-up)

The N1 virial is the strain derivative at fixed radii. Frozen: the D4 evaluation radius at min(current,
nextafter(L_min/2 (1.15)^(-1/3))), list radius eval + 10 (so that the density may rise by 15 % before the D4 radius
alone would stop the run; the tail covers the part beyond r_c, A3c), `dispersionRadiiBohr()` then returns the frozen
values and the D4 list is rebuilt once (the per-pair `r_cut` stamp). Absolute already: Ewald r_cut (`m_ewald_params`),
repulsion 20, HB/XB search radii sqrt(hbthr1/2) = 15.8 / 21.21 Bohr, XB 20. W1 boxes (432 H2O, L0 = 44.39 Bohr): D4
radius 21.19 Bohr (< 21.21), the HB/XB search radius stays the limiter: L_min > 42.43 Bohr, density +12.9 %. A
violated check at a coupling: `[ERROR] NPT: cell shrank below the cutoff limit L=... Bohr (needs > ... Bohr, <cutoff>) at
step ... (t, mu, P); restart written: <basename>.npt_limit.json. Reason: ...`, the unscaled state is restored and saved,
`StopReason::CellLimit`, `finalizeRun()` runs. An explicit `dispersion_cutoff_bohr` is kept as given.

### Output, restart, trajectory

Status line: N1's P, <P>, V, rho plus `L` (shortest edge, A) with a barostat. Restart: `cell` (9 numbers, A; also with
`compute_pressure` alone), `barostat`, `pressure`, `tau_p`, `compressibility`, `pcoupl_interval`, `pcoupl_calls`,
`pcoupl_events`, `average_P`; a loaded restart sets the molecule's cell before the energy set-up and the barostat state
(the file overrides the CLI as for the thermostat). Trajectory frames: `... unitcell a b c 90.0 90.0 90.0` appended to the
comment line, only with a barostat (without one the trajectory stays md5-identical to N1; operator decision whether every
periodic MD should write it); a frame read back gives the cell (`-sp` on the last frame of the crescale run: periodic cell
44.5195 Bohr = 23.5587 A).

### Tests (`ctest -R pbc_npt`)

- **(c) `pbc_npt_cell_update`** (`test_cases/unit/test_pbc_npt.cpp`, 20 waters random-packed in L = 44 Bohr, `ewald_alpha`
  0.32, `ewald_r_cut` 19, `ewald_tolerance` 1e-10, dense EEQ): frozen D4 radius 20.998587 Bohr = 22 (1.15)^(-1/3);
  instance A at the start == fresh instance with `dispersion_cutoff_bohr` = that value: dE 0, max|dG| 0, d tr Xi 0.
  A scaled by mu = 1.001 (`updateCell` + `UpdateMolecule`) vs a fresh instance on the scaled structure (same alpha, r_cut,
  k_max): **dE 0.0 Eh, max|dG| 0.0 Eh/Bohr, d tr Xi 0.0 Eh** (E -6.544991750852 vs start -6.544081495695; tr Xi 0.876717245);
  round trip to mu = 1: 0 against the start; a cell of 39.6 Bohr refused (four cutoffs named) with the instance unchanged;
  4 threads after the update vs 1 thread: dE 6e-15, max|dG| 1e-15.
- **(d) limit abort** (same test, 20 waters, Berendsen P0 = 2e5 bar, tau_p 100 fs, every step, dt 0.5: mu = 0.985 per
  step): stop at step 2 with `NPT: cell shrank below the cutoff limit L=42.0454 Bohr (needs > 42.4264 Bohr, HB A-H-B / XB
  A-B search radius sqrt(hbthr2))`, `StopReason::CellLimit`, kept cell 42.6832 Bohr (the last valid one),
  `pbc_npt_limit.snapshots/pbc_npt_limit.npt_limit.json` with `cell` and `barostat`.
- **(e) restart** (same test, crescale, every 2nd step): 40 steps = 20 cell changes (L 23.2838 -> 23.319076 A); a new
  SimpleMD from the saved state carries that cell and `pcoupl_events` 20, 20 more steps give 30 and a moved cell.
- **(a) `pbc_npt_idealgas_crescale` / `_berendsen`** (`test_pbc_npt_idealgas.cpp`, 64 Ne on a 4^3 lattice in L = 200
  Bohr, CSVR 300 K tau 100 fs, dt 1 fs, coupling every step, compressibility 1/P0, P0 = 65 kT/V0 = 2.2710 bar,
  `ewald_r_cut` 60, `ewald_tolerance` 1e-4, matrix-free pcg, `dispersion_cutoff_bohr` 40; E_pot of the lattice -1.1e-7 Eh;
  200 ps, 5 ps equilibration, 99 blocks of 2 ps; dof = 189, N_eff = 64): crescale (tau_p 50 fs) **<V> = 64.29 +- 0.28
  kT/P0** (block SE 0.44 %; N_eff = 64: 1.0 SE; the paper's N + 1 = 65 would be 2.5 SE off, the Berendsen value 63
  4.6 SE), **relative variance 0.01502 = 0.961/N_eff**, <P_virial> 5.8e-3 bar = 0.25 % of P0 (collisions: second virial
  coefficient; at L = 120 Bohr it was 1.1 %, and a 3.4 sigma volume excursion there hit the D4 limit L > 100 Bohr after
  75 ps: a legitimate NPT fluctuation, hence the larger box with fixed radii). Berendsen (tau_p 100 fs): <V> = 63.55 +-
  0.28 kT/P0 (expected (dof/3)(1 + var) = 63.46, 0.3 SE; <P> = 2.2707 vs P0 2.2710 bar), relative variance 0.00736 =
  0.471/N_eff (what is left is the thermal noise of the kinetic pressure integrated over tau_p: ~ tau_corr/(2 tau_p));
  with tau_p 200 fs (the ctest setting) <V> = 63.56 +- 0.28 kT/P0 (expected 63.34, 0.8 SE), relative variance 0.00543 =
  0.348/N_eff, <P> 2.2706 bar. 205000 steps take ~60 s per test (0.3 ms per step; n = 1 run per setting, seed 7).
- **(f) bitwise** (binary saved before N2): no cell, `-sp -dump_gradient` energy (16 digits) and gradient md5 on benzene,
  caffeine, methyl-a-d-glucopyranoside, water8, hydronium-water2 (+1), urea-water197, peo201-water500 at 1 and 4 threads:
  14/14 identical. With a cell, no barostat: `-sp -dump_gradient` of box60 (180 atoms, L 44 Bohr) and box970 (2910 atoms)
  identical; `-md` 100 steps (dt 0.5, seed 7, 4 threads): energy table (13 columns), trajectory (comment lines removed:
  they carry the git tag) and last frame md5 identical for both boxes.
- **(g)** CPU `ctest -R 'gfnff|frag|eeq|simplemd|pbc|ewald|react|restart|md_|gpu|plugin|molecule|xyz|parser|curcumaopt' -j1`:
  147 run, 100 % passed (148 listed: 4 skipped GPU, 1 disabled), incl. the 3 new `pbc_npt_*`, `xyz_comment_unitcell`,
  `cli_curcumaopt_08_not_converged_report`; `scripts/check_docs.py`: the 9 known errors only (`docs/NPT.md` marked
  intent-to-add, `git add -N`, nothing committed).
- XYZ parser and `-opt` display fix (separate parts of the order): ctests `xyz_comment_unitcell`,
  `cli_curcumaopt_08_not_converged_report` (AIChangelog.md).

### (b) 432-water box (W1 `box432_s1.xyz`, 1296 atoms, L0 23.4905 A = 44.39 Bohr, 0.997 g/cm^3)

1 ps CSVR heat-up (NVT, T 298.15, tau 100 fs, dt 0.5, seed 7; the TIP3P-like packing relaxes: T 298 -> 474 K at 0.2 ps,
351 K at 1 ps, P 23874 bar at step 0, -5432 at 0.2 ps, 1360 at 1 ps), then 2 ps from its restart file: NVT with
`compute_pressure`, Berendsen NPT, C-rescale NPT (P0 1.01325 bar, tau_p 1000 fs, beta 4.5e-5, every 10 steps), 4 threads.
Both NPT runs finish without an error; the cell moves smoothly: Berendsen L 23.5224 -> 23.5927 A (max |dL| per 100 steps
0.026 A), rho 0.9929 -> 0.9841 g/cm^3, <P> over the 2 ps 6 bar (std 545, last ps 129); C-rescale L 23.4932 -> 23.5587 A
(max |dL| per 100 steps 0.085 A), rho 0.9967 -> 0.9884, <P> 395 bar (std 1452, last ps 128); NVT reference <P> of the
same 2 ps -208 bar (instantaneous -1440 .. +1270). E_pot still falls in all three runs (-146.87 -> -147.44 / -147.60 /
-147.56 Eh) and T is at 300-320 K: **the box is relaxing, not equilibrated; the densities are not a result** (the W1
protocol equilibrates 10 ps before 20 ps production). Excerpt (t ps, Epot Eh, T K, P bar, rho g/cm^3, L A):

| run | 1.0 | 1.7 | 2.1 | 2.5 | 2.9-3.0 |
|---|---|---|---|---|---|
| NVT | -146.872, 343, 1270, 0.9970 | -147.26 (1.6), 309, -678 | -147.44 (2.0), 300, -584 | -147.47, 309, -373 | -147.61 (2.8), 308, -1277 |
| Berendsen | same start | -147.351, 307, 55, 0.9865, 23.574 | -147.468, 304, 107, 0.9912, 23.536 | -147.52 (2.6), 314, -1071, 0.9875, 23.566 | -147.600, 302, 183, 0.9849, 23.587 |
| C-rescale | same start | -147.257, 308, 2358, 0.9955, 23.502 | -147.409, 320, -1412, 0.9890, 23.554 | -147.456, 302, 1350, 0.9904, 23.542 | -147.557, 309, 2170, 0.9884, 23.559 |

Time per step (4 threads, difference of the 4000- and 2000-step runs from the same restart, 3 repetitions):

| run (3 runs in parallel on a 32-core machine, with other work) | per step (ms), 3 repetitions | median |
|---|---|---|
| NVT + `compute_pressure` | 44.2 / 43.2 / 40.5 | 43.2 |
| Berendsen, every 10 steps | 46.7 / 46.4 / 43.4 | 46.4 |
| C-rescale, every 10 steps | 46.8 / 46.0 / 43.3 | 46.0 |

The barostat costs +3 ms per step (+7 %): the forced rebuild of the D4 list, the repulsion list and the HB/XB detection
at every 10th step (the HB/XB refresh runs every 10 steps anyway, `hb_update_force_every`), the new Ewald lattice (k
vectors, cache rebuild: the EEQ pcg starts from the kept history) and the scaling itself (O(N)). The three repetitions
vary by 10 % (shared machine), the +3 ms is consistent across them (n = 3 per mode).

### Not done / not tested

- The C-rescale equations against the paper's text (above); anisotropic or semi-isotropic coupling; a molecular
  (COM-based) virial and scaling; the constraint virial (RATTLE refused, so NPT production runs without constraints);
  the GPU hybrid (refused); `topology_mode react` with a barostat (scan forced at every cell change, not run); charged
  boxes; a liquid at equilibrium (only 2 ps of a relaxing packing); boxes above 1296 atoms with a barostat; restart
  files with a barostat and RATTLE. The forced list rebuild at every coupling could be replaced by a skin that budgets
  the strain (|mu - 1| L sqrt(3)/2 per coupling) if the rebuild cost matters at `pcoupl_interval` 1. Plugin ABI 4: the
  CUDA plugin must be rebuilt (not built here). Human production testing pending.

## V-LRC long-range pressure correction of the D4 tail (2026-11-09, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `0351a245`), CPU only. Trigger: the P1 cutoff scan of the W1 water project (project folder
`03_auswertung/cutoff_scan/BERICHT.md`): with the tail the pressure still followed `dP(r_c) = c + A r_c^-3.15`, with an
amplitude about equal to the tail's own pressure; without the tail about twice that.

### Derivation (code comment at `FFWorkspace::postProcessVirial`)

One pair class `u(r) = -C_n r^-n`, N atoms, rho = N/V, g(r) = 1 beyond the sharp radius r_c:
- energy tail `E_n = (N rho/2) int_rc^inf 4 pi r^2 u dr = -2 pi N rho C_n / ((n-3) r_c^(n-3))` (the A3c formula),
- virial tail `W_n = -(N rho/2) int 4 pi r^2 r u'(r) dr = n E_n` (since `r u' = -n u`), pressure `P_n = (n/3) E_n / V`.
This is the classical long-range (tail) correction of the virial pressure for a truncated pair potential (as in Allen and
Tildesley, Computer Simulation of Liquids, long-range corrections; checked here against their Lennard-Jones forms: the
r^-6 parts of E_lrc and P_lrc follow from the two lines above with C_6 = 4 eps sigma^6).
The N1 virial differentiates the model energy at a fixed pair set: `E_n ~ 1/V`, tr Xi = 3 E_n, i.e. `P = E_n/V`, factor 1
instead of n/3. The difference `(n-3)/3 E_n/V = -(2 pi/3) rho^2 r_c^3 u(r_c)` is the impulsive term of the pairs that cross
r_c under a strain (energy jump u(r_c) each); a finite-difference strain at constant pair count (the N1 gate) cannot see it.
With the two D4 tail parts (C6, n = 6; C8, n = 8): **tr Xi_lrc = 3 E6_tail + 5 E8_tail** (the order's derivation is
confirmed, factors 2 and 8/3 of the full tail pressure instead of 1). The kernel is a hard cutoff (`rij > disp.r_cut`:
skip), so the crossing term exists. Approximations: g(r) = 1 beyond r_c and the large-r form of the damping (as E_tail,
error <= (1/3)(R0/r_c)^6 relative); a statistical estimator, not a derivative of the model energy.

### Implementation

- `FFVirialComponents::tail_lrc = (3 E6 + 5 E8)/3 * I`; `total()` unchanged (model derivative), `pressureTotal()` =
  `total() + tail_lrc`; `FFWorkspace::virial()` (and with it `GFNFF::virial()`, SimpleMD's P column, the barostat) returns
  `pressureTotal()`, `FFWorkspace::modelVirial()` the model part. E6/E8 from the same sums S6/S8 as the tail energy
  (`homogeneousTailFromSums`), passed with `setDispersionTail(e, dEdcn, e6, e8, pressure_lrc)`; `dispersionTailE6/E8()`.
- PARAM `gfnff.pressure_lrc` (Bool, default true, category Dispersion; constructor fallback true). false: tail_lrc = 0,
  pressure as N1/N2 (bitwise). Only with `dispersion_tail` on; no effect on energies, forces or without a cell. The C8 part,
  dE_tail/dCN and the CN chain rule are unchanged. Plugin ABI 5 -> 6 (`gpu_plugin_abi.h`; the GPU hybrid has no virial).

### Numbers

- Estimate: water at 1.026 g/cm^3 (432 box, rho_mol 5.08e-3 Bohr^-3), <zeta_i zeta_j C6> per molecule pair 40 au (from
  E_tail -0.01992 Eh at r_c 21.00 Bohr, V 84990 Bohr^3): E6/V = -(2 pi/3) rho^2 <C6>/r_c^3 = -69.0 / -34.7 / -17.1 / -8.7 bar
  at r_c 21.0 / 26.4 / 33.4 / 41.9 Bohr; P1 "tail off" +67.2 / +33.2 / +16.4 / +8.2 bar and the P1 fit P(inf) - P(ref)
  -70.0 / -34.0 / -16.4 / -8.0 bar. Measured shift of P at the matrix radius (new binary, lrc true - false, 5 frames):
  -70.6 / -34.3 / -16.8 bar (432 / 864 / 1728); C8 share of E_tail 4-5 % (water boxes in the test).
- Density effect of the old estimator (linear, kappa_T 4.5e-5 /bar, barostat setting): NPT densities with the matrix radii
  were 0.32 / 0.15 / 0.08 / 0.04 % too low (432 / 864 / 1728 / 3456; 3456 from the P1 fit).

### Tests (thresholds fixed before the runs)

- (a) `pbc_virial` (FD gate against `total()`) unchanged and passing with the default `pressure_lrc true`.
- (b) ctest `pbc_virial_lrc`: 20/60 H2O in L = 44 Bohr (r_c 20 and default L/2): E6 + E8 = E_tail (0 rel), tr tail_lrc =
  3 E6 + 5 E8 (0 rel), isotropic, `virial()` = `total()` + tail_lrc (8e-17), pressure_lrc false: tail_lrc 0 and virial()
  == total() bitwise, model components bitwise equal; E6/E8 from E_tail at r_c 15 and 20 (two-radius fit) 4e-16 / 4e-15.
- (c) same test, 800 Ne (random, d_min 5 Bohr, L 70 Bohr, rho 2.33e-3 Bohr^-3, 5 configurations, fixed Ewald): paired
  dP(r_c) - dP(30) model +6.17 / +2.02 / +0.61 bar, with lrc -0.12 / -0.013 / -0.004 bar (r_c 15 / 20 / 25, SE 0.02 /
  0.004 / 0.001); ratio at 20 Bohr 0.006 (criterion <= 0.2); max |q| 4e-14 e, |E_coul| 5e-13 Eh. The -0.12 bar at 15 Bohr
  is of the size of the damping and d_min-correlation error (not separated).
- (d) cutoff scan, project folder `03_auswertung/cutoff_scan_lrc/` (script `cutoff_scan_lrc.py`, P1 protocol, 5 frames,
  paired against the matrix radius): dP over 21..27 Bohr with lrc 0.15 bar (864), 0.14 bar (1728) against 37.7 / 38.1 bar
  without (ratio 0.004, criterion <= 0.30). 432 (L/2 22 Bohr: 24 and 27 Bohr not admissible, cell check aborts; 1728 run
  instead as the second box): over 15..21.48 3.4 bar against 137 bar (0.025). With lrc the residual is -3.4 to -3.6 bar
  at 15 Bohr and <= 0.4 bar from 18 Bohr on in all three boxes. `pressure_lrc false` reproduces the old scan exactly
  (P 0.000 bar, E_tot 0 Eh, 85 runs).
- (e) without `compute_pressure`: `-sp` energy and `-dump_gradient` md5 of the 7 structures at 1 and 4 threads 14/14
  identical to the HEAD binary; `-md` box60 100 fs with a cell: table (timing columns excluded) and trajectory coordinates
  identical; with `compute_pressure` / Berendsen (432 box, 40 steps) and `pressure_lrc false` identical incl. P, V and the
  final restart file; with lrc true P (and the NPT cell) differ as intended.
- (f) CPU `ctest -R 'gfnff|frag|eeq|simplemd|pbc|ewald|react|restart|md_|gpu|plugin|molecule|xyz|parser|curcumaopt' -j1`:
  "100 % tests passed out of 149" (4 GPU tests skipped, 1 disabled), incl. `pbc_virial` and `pbc_virial_lrc` (16 s);
  `scripts/check_docs.py`: the 9 known errors only.
- **Not tested**: liquids other than water and the Ne gas, charged boxes, NPT production runs with the correction (the
  density effect above is the linear estimate), the GPU hybrid (no virial), the CUDA plugin at ABI 6 (not built). No
  crossing term exists for the other sharp radii (HB/XB 21.2, repulsion 20 Bohr). Human production testing pending.

## F-MD four SimpleMD defects from the matrix and pilot runs: print_frequency, restart seed, restart DOF, seed flakes (2026-11-09, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/bulk` (base `46c35811`), CPU only, all in `src/capabilities/simplemd.cpp/.h`. New ctests
`simplemd_print_frequency` and `simplemd_restart_state` (`test_cases/unit/test_simplemd_fmd.cpp`), both red with the
sources of `46c35811` (SIGFPE, exit 136; 5 of 16 checks) and green after the fix.

1. **`-print_frequency 0.5` ended with SIGFPE.** The PARAM was `Int` and read with `get<int>`: 0.5 became 0, used as
   the modulus in `static_cast<int>(m_step * m_dT) % m_print` (old line ~3244). print_frequency is a time in fs (as the
   test always treated it; the help said "steps"). Fix: PARAM `Double`, `m_print` double, test in `printDue()`; a value
   <= 0 prints every step, a fractional one prints when `floor(step*dt / freq)` changes, an integer one keeps the old
   test unchanged.
2. **Repeated table rows at dt < 1 fs: not fixed.** The test truncates the time, so with dt 0.5 fs the steps 100 and
   101 both print the 50 fs row (and step 1 a row at 0.5 fs). A variant that kept only the first step per truncated time
   was built and ran: rows and values otherwise identical, but `gfnff_d4_fused_c6` (`check_d4_fused_c6.py`) counts at
   least 20 rows for 20 steps of 0.5 fs with `print_frequency 1`, i.e. it relies on one row per step, and failed. Per the
   rule "existing tests must not change" integer frequencies keep the old test exactly; the project parsers read the
   50 fs grid and tolerate the repeats. A fix needs that test adapted (operator decision).
3. **After a restart the generator was not seeded.** `seedRandomEngine()` was only called from `InitVelocities()`, which
   a restart skips: every restarted run drew the default-constructed `mt19937` sequence (CSVR, Andersen, C-rescale noise
   identical for every repetition; seeds 202 and 303 bit-identical). Fix (`Initialise()`, after the restart load): an
   explicit `-seed` (in the `simplemd`/`global` scope or flat, detected in the constructor) wins over the file's `seed`;
   without `-seed` the file's value is used (note: the stored value is the start seed plus one per CSVR call, so a
   continuation draws a new stream); `-seed -1` (clock) and `0` resolve as before; then `seedRandomEngine()` (seed and
   restart geometry). Test: CSVR, 20 steps after the restart: 202 vs 303 differ (max dEkin 2.3e-2 Eh), 202 twice and
   "no -seed" twice identical, "no -seed" == `-seed <file value>`.
4. **After a restart m_dof was 3N.** `InitConstrainedBonds()` (the only DOF count) was skipped on restart. The count is
   now `countDegreesOfFreedom()` (split off unchanged) and also runs after a restart, with the constraint list loaded from
   the file (not rebuilt; RATTLE in the file without a list warns). Test, water dimer, DOF fresh / after restart: vacuum
   uff 12/12 (before 18), cell + GFN-FF 15/15 (before 18), RATTLE 8/8 (before 18), cell + uff 12/12 (before 18).
5. **Seed flakes**: `cli_gfnff_09_farfield_shortcut` (MD part) and `cli_simplemd_16_md_seed` ("no seed" pair) ran
   without `-md.seed`; they now use `-md.seed 1` and `-md.seed 0`. The default of `-md.seed` (-1) is unchanged.

**Runs made before this fix** (frozen matrix binaries keep all four defects): fresh starts are not affected by 3 and
4. A run continued with `-restart_file` and **no thermostat or barostat** (the `e_nve` restarts) has the same
trajectory as after the fix (DOF and the generator only enter thermostats, barostat noise and `adaptive_step`); only
the printed T and <T> are too low by the factor (3N - n_removed - n_constraints)/3N (box: 3/3N, e.g. -0.08 %), P is
not affected (Ekin enters directly). A restarted run **with** CSVR/Berendsen/Nose-Hoover regulated to a kinetic energy
too high by 3N/(3N - removed - constraints) (box +0.08 %, droplet about +1 %), and with CSVR/Andersen/C-rescale its
noise was the same for every repetition and every seed: such repetitions are not independent samples.

Checks: without restart and with integer print frequency bitwise unchanged: `-sp -dump_gradient` (16-digit energy line
and md5) on the 7 structures at 1 and 4 threads 14/14; `-md` caffeine and water8 gfnff, 200 steps of 1 fs, seed 7:
trajectory identical apart from the version string, final restart file byte-identical, table rows identical (with
`-print_frequency 50` and with the default 1000); dt 0.5 fs with `-print_frequency 2`: the same rows as before.
CLI restart (water8, CSVR, MaxTime edited to 40 fs): `-seed 202` twice and `-md.seed 202` identical, `-seed 303`
different, no `-seed` twice identical; DOF after the restart 66 = 3N-6. `ctest -R 'gfnff|frag|eeq|simplemd|pbc|ewald|
react|restart|md_|gpu|plugin|molecule|xyz|parser|curcumaopt' -j1`: "100 % tests passed out of 152" (4 GPU skipped, 1
disabled); `cli_gfnff_09_farfield_shortcut` and `cli_simplemd_16_md_seed` 30 times each, 0 failures; `check_docs.py`:
the 9 known errors only.
- **Trap met on the way**: the first build declared `PARAM(print_frequency, Double, 1000, ...)`; the registry stores the
  literal as `int`, `getDefaultJson()` skips a default whose `std::any` cast fails, and every MD without an explicit
  `-print_frequency` threw `ConfigManager: Parameter 'print_frequency' not found` (23 ctests red). A `Double` PARAM needs
  `1000.0`. The new unit test sets the flag explicitly and could not see it; the CLI ctests did.
- **Not changed**: the other restart fields still override the CLI (`thermostat`, `dT`, `MaxTime`, `T`; TODO "MD-Restart:
  Datei schlaegt CLI"); a RATTLE restart takes its constraint list from the file; `dump_frequency 0` and the other step
  moduli are not guarded. Not tested: Nose-Hoover after a restart, remove_com_mode 2/3 after a restart, GPU MD. Human
  production testing pending.

## A1-GPU atomic virial and NPT of the CUDA periodic hybrid (2026-11-23, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/gfnff-gpu-virial` (base `94303824`, devel/methods). Item A1 of [GPU_PERIODIC_GAPS.md](GPU_PERIODIC_GAPS.md):
`-md -compute_pressure true` and `-barostat berendsen|crescale` run with `-gpu cuda -gfnff.gpu_periodic true` (S2
hybrid). The virial is the N1 one (atomic, fixed fractional coordinates, every term).

### Where each part comes from

| part (FFVirialComponents) | evaluated by |
|---|---|
| bond (incl. the HB-CN modulation of the X-H bonds), D4 pair, bonded and non-bonded repulsion | separate device kernels `k_bonds_virial`, `k_hb_alpha_virial`, `k_dispersion_virial`, `k_repulsion_virial` (pair vector and dE/dr formed as in the energy kernels, which are untouched) |
| angle, dihedral (incl. extra torsions), inversion, sTors, HB (frame on H, every gradient row with its frame position: A, B, neighbours of B, C, D, LP neighbours), XB (frame on X), BATM (plain positions, as the energy) | VIR instantiation of the kernel body: `k_*_body<PBC, VIR>` is a `__forceinline__` device function, `k_*<PBC>` (former signature) calls it with VIR false, `k_*_vir` with VIR true |
| CN chain rule (bond dr0/dCN, dC6/dCN incl. the D4 tail, TERM 1b chi(CN)) | `k_cn_chainrule_virial` with the final device dE/dCN: one tensor, stored in `cn_bond` (`cn_dispersion`, `cn_coulomb` zero; the device keeps one dE/dCN vector) |
| Coulomb Ewald real + reciprocal pair part, background terms | host, `FFWorkspace::evaluateEwaldCoulomb` (the CPU passes `realSpaceVirialRange`, `reciprocalVirialRange`, `backgroundVirial`) |
| D4 tail `E_tail I`, V-LRC `(3 E6 + 5 E8)/3 I` | device tail energy; E8 from a second `k_d4_tail_atoms` launch with pref6 = 0 (no dE/dCN) |

Reduction: per thread nine doubles, rounded to 2^-48 Eh and summed as integers (warp, block, one atomic per block and
component; `blockReduceAddVirial`), so the result does not depend on thread order (docs/GPU_REPRODUCIBILITY.md). Slot
layout `cuda/gfnff_virial_slots.h`. Assembly: `GFNFFGpuMethodImpl::assembleHybridVirial()` ->
`FFWorkspace::setHybridVirial()`, read through `GFNFF::virial()` as on the CPU. With the virial on, the CUDA graphs are
not used; with it off every launch is the former one.

Cell changes: `EnergyCalculator::setUnitCell` -> `GFNFF::updateCell` on the host instance (N2). The wrapper
(`followHostCellChanges()`) sees the changed `cellUpdateCount()` or frozen D4 radius at the next call: device tail
tables (1/V, r_c), CN-derivative pair list, graphs; the cell symbols are uploaded every step anyway; the host list
decisions (HB/XB refresh, repulsion, D4) see the rebuild flag and the device rebuilds follow; the wrapper clears the
flag after them (`GFNFF::acknowledgeListsRebuilt()`, the CPU does it in `Calculation`). `SimpleMD::setupBarostat` accepts
`-gpu` when the method is GFN-FF with a periodic energy (the hybrid); `coulomb_ewald false` (S1, not periodic) stays
refused. Plugin ABI 25 -> 26.

### Finding on the way: an added kernel parameter changes the last bits

A first version added the virial pointer as a last parameter of the six templated kernels. Without a cell the GPU
gradients then differed from the base build in the last digits on 5 of 6 structures (energy 1e-14 Eh). The SASS of the
`false` instantiations had another instruction mix (k_dihedrals DADD 96 -> 75, k_hbonds DMUL 717 -> 714,
k_storsions DFMA 202 -> 200) although no statement of that path changed. Now the kernels keep their former signature and
call a shared body: `cuobjdump -sass` of the plugin, all 305 functions of the base plugin instruction-for-instruction
identical (12 new functions).

### Acceptance (thresholds fixed beforehand)

| check | result |
|---|---|
| ctest `gfnff_gpu_pbc_virial` (a): GPU vs CPU per term (14 parts + total + pressure virial) and tensor component; 20 / 60 H2O, H3O+(H2O)2 + 18 H2O (Q +1), caffeine, acetic dimer + F-I...pyridine (XB) in L 44 Bohr | worst 5.0e-13 Eh (limit 1e-7) |
| same, W1 frames and boxes (4 threads, production settings): 432 H2O (equilibrated frame and `box432_s1`), 864, 1728 H2O, urea400 + water1000 with an assigned cell L = 42 A (not relaxed) | 3.6e-12 / 4.2e-12 / 3.4e-12 / 2.7e-11 / 2.8e-11 Eh; total and pressure virial <= 4.8e-11 Eh; dE <= 2.7e-11 Eh, dG <= 8.7e-14 Eh/A |
| (b) FD gate of the hybrid energy, fresh GPU instance per strain, eps 1e-5 / 2e-5 Richardson, isotropic and axial eps_xx, eps_yy, eps_zz; 20 H2O and the XB composite | <= 5.6e-9 Eh (limit 1e-7) |
| (c) second GPU instance (60 H2O, 4 threads): energy, gradient, every tensor; energy and gradient with the virial request vs without | bit for bit (both) |
| (d) `updateCell` (x1.001) on the GPU instance vs a fresh GPU instance on the scaled box | energy, gradient and virial: difference 0 (limits 1e-10 Eh, 1e-10 Eh/A, 1e-9 Eh) |
| (e) SimpleMD Berendsen 20 steps, 60 atoms, CPU vs `-gpu cuda` (no fallback counted) | 10 couplings each, cell edge equal to 1e-9 A, P 647.0860 bar both |
| CLI NPT 432 H2O (W1 frame, 20 steps, CSVR 300 K, Berendsen every 2 steps, 4 threads) | 21 status lines: Epot (6 decimals), P (2), V (3), rho (5), L (5) identical CPU vs GPU |
| no cell, 6 structures x gfnff/gfn1/gfn2 x 1/4 threads + 2 MD, `-sp -dump_gradient` against `merge-devel` builds of the same options | CPU 38/38 lines identical; CUDA 38/38 identical |

Note on the reference builds: the first comparison failed for every method including gfn2, because the build
directories of this branch had been configured with the options as one quoted string (zsh does not split `$OPTS`),
i.e. without MKL/AVX/march native. With the options as separate words the comparison above is clean.

### Time per NPT step (ms; CLI, 4 host threads, Berendsen every 2 steps, difference of 60- and 20-step runs, n = 3; RTX 5080)

| box (W1 frame) | CPU | GPU hybrid | gain |
|---|---|---|---|
| 432 H2O (1296 atoms) | 48.0 / 44.9 / 48.0 | 26.3 / 28.2 / 30.6 | 1.7x |
| 864 H2O (2592) | 117.9 / 119.5 / 119.7 | 69.2 / 71.1 / 67.8 | 1.7x |
| 1728 H2O (5184) | 316.2 / 317.1 / 319.2 | 189.2 / 188.5 / 186.5 | 1.7x |

216 H2O (`water216-box.sc-s1`, L 18.63 A = 35.2 Bohr) cannot be run with a cell: the HB/XB search radius needs L > 42.43
Bohr (refused by the cell check, CPU and GPU alike). The host EEQ and Ewald (S3/S4, item A2) remain the main cost.

### Not done / not tested

- ROCm (B3; the wrapper code is guarded by `has_periodic_local_terms`, not compiled here), multi-GPU (off with a
  cell), FP32 mixed precision with the virial (BATM takes the FP64 kernel, repulsion keeps its FP32 gradient but an FP64
  virial), the CN chain-rule split per term on the device, `topology_mode react` with the virial, charged boxes under
  the ensemble charge model (CPU fallback), a long NPT run on the GPU (only 20 steps here), C-rescale on the GPU (shares
  the code path, not run). Human production testing pending.

## A3-P pressure terms of the sharp GFN-FF radii, measured (2026-11-23, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

Item A3 of [GPU_PERIODIC_GAPS.md](GPU_PERIODIC_GAPS.md); CPU path (the GPU hybrid has the same virial, A1-GPU). Tool:
`test_pbc_pressure_cutoffs --xyz FRAME --strain --hb f1,f2,...` (ctest `pbc_pressure_cutoffs` checks only the per-term
split). Frames: W1 cutoff-scan frames (`03_auswertung/cutoff_scan/runs/frames_*`, equilibrated Berendsen NPT, 300 K),
D4 radius frozen as an NPT run freezes it, production settings otherwise, 2 threads.

**Method.** The N1 virial is -dE/d(eps) with every pair and triple list fixed. A list with a sharp radius gains or loses
members under a strain; a member whose energy at the radius is not zero then adds an impulsive term to the pressure of
the truncated model (V-LRC derives it for the D4 radius). Measured per term: `P_thermo = -(dE_t/d eps)/(3V)` from a
least-squares quadratic fit of E_t over 9 isotropic strains eps = 0, +-1, 2, 4, 8 x 1e-4, a fresh instance (lists
re-derived) at every strain, against `P_vir = tr Xi_t/(3V)` at eps = 0; crossing = P_thermo - P_vir. Fit residuals
<= 1.8e-6 Eh, i.e. a slope uncertainty below 1 bar.

| frame | HB crossing | D4 crossing | V-LRC (in the reported P) | repulsion | bond (CN) | Coulomb (CN) | total crossing |
|---|---|---|---|---|---|---|---|
| 432 H2O, t 3.9 ps | -37.7 | -71.0 | -71.2 | +1.31 | -0.29 | -0.17 | -107.9 |
| 432 H2O, t 11.9 ps | -31.2 | -70.3 | -71.0 | +1.30 | -0.28 | -0.17 | -100.6 |
| 432 H2O, t 19.9 ps | -30.9 | -68.0 | -66.3 | +1.28 | -0.27 | -0.17 | -98.0 |
| 864 H2O, t 19.9 ps | -30.5 | -33.2 | -33.6 | +1.29 | -0.28 | -0.17 | -62.8 |

(bar; n = 4 frames, 2 box sizes.) The same 432 frame with the |E_HB| candidate filter off (`hb_min_pair_energy_eh 0`,
933 431 instead of 152 135 triples): HB crossing -30.8 bar, HB virial pressure -0.2 bar, so the filter's own crossing
term is about 0.1 bar.

**HB search radii (thr1, thr2 scaled together, fixed-list virial pressure relative to the default, 1728 H2O frame
t 19.9 ps, L 69.8 Bohr):** x0.8 +87.0, x1.2 -29.5, x1.5 -44.9, x2.0 -50.7, x2.5 -51.1 bar (P1 measured +87 / -29 on
other frames); E_HB changes by -2.3e-6 Eh per atom from x1 to x2.5: the radius is part of the model (reference
value hb_accuracy 0.1), not a numerical parameter. 432 frame: x0.8 +81.1, x1.05 -9.2 bar (larger factors exceed L/2).

**Findings.**
- V-LRC is confirmed by an independent measurement: the measured D4 crossing term equals the V-LRC estimate to
  0.4-1.7 bar (3 % at most).
- The HB search radius has a crossing term of -31 to -38 bar that the reported pressure does not contain: the barostat
  sees a pressure about 31 bar too high for the model it simulates; NPT densities are then too low by about 0.2 %
  (linear estimate with kappa_T ~ 6e-5 /bar, a rough number, not the measured compressibility of GFN-FF).
- Against the converged search radius (x2 and beyond) the default-radius pressure is 51 bar too high (density
  about -0.3 % by the same rough estimate); the crossing term covers 60 % of that. Which of the two is wanted is a
  model decision: the crossing term makes the pressure consistent with the reference model as it is (its truncated
  HB energy), the radius limit changes the model.
- Repulsion (+1.3 bar), CN-dependent bond and Coulomb parts (-0.45 bar together) are negligible.

**Not changed (operator decision).** A per-step HB crossing estimator would need the energy density of the triples at
two different sharp surfaces (r_AB^2 < hbthr1; case 1 also r_AB^2 + r_AH^2 + r_BH^2 < hbthr2) plus the energy filter,
estimated from one side of the surface each step, and the device-detected HB lists on the GPU; it is not a cheap,
exact term like V-LRC (whose g(r) = 1 tail is analytic). Enlarging the radii only for the pressure would make the
pressure inconsistent with the forces and roughly triple the HB work (1728 H2O: 0.41 / 0.89 / 2.0 million triples at
x0.8 / x1.2 / x2; the default lies between, about 0.6 million by the 432 count of 152 135 per 432 waters).
Options for the operator: (a) leave as is and correct NPT densities by about +0.2 % in the analysis, (b) a statistical
HB crossing estimator in the virial (CPU and GPU), (c) a larger hb_accuracy-derived radius as a deliberate model
deviation (switch, port fidelity stays the default). Not measured: other liquids, XB (no halogen in the frames),
charged boxes, the GPU-detected lists (the CPU lists are used here).

## S3 device Ewald of the CUDA periodic hybrid: Ewald sum and pcg operator on the device (2026-11-23, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

On `feature/gfnff-gpu-ewald` (base `f7229704`, A1). Item A2 of [GPU_PERIODIC_GAPS.md](GPU_PERIODIC_GAPS.md), stages S3a
(Ewald pair sum at given charges) and S3b (operator of the periodic pcg), built and gated together. Default in the hybrid
(`-gpu cuda -gfnff.gpu_periodic true`, `coulomb_ewald` true); `-gfnff.gpu_ewald_device false` restores the S2/A1 host Ewald
in the same binary.

### What runs where

| part | where |
|---|---|
| pair list (full list per atom, ascending j), kernel K_ij = erfc(a r)/r - erfc(g_ij r)/r and dK/dr, phase table | device, once per geometry (`GfnffEwaldGPU::prepare`, key: positions, alp, cell, parameters bit for bit) |
| operator of the periodic pcg: real-space gather + reciprocal (structure factors, per-atom pass) | device (`potentials`), one call per CG iteration |
| Ewald pair energy, gradient (Eh/Bohr), virial (real + reciprocal) at the final charges | device (`evaluate`) inside `FFWorkspace::evaluateEwaldCoulomb` |
| self, background, Gaussian background (O(N)), diag_local; CG loop (start vector, aspc, block preconditioner, projection, dots, stopping rule) | host, the CPU code unchanged (`calculateFinalChargesPeriodic`) |
| dense EEQ path (N <= `eeq_ewald_dense_max_atoms`, 1000) | host, as before; the Ewald energy then still runs on the device |

Interface: `curcuma::ewald::EwaldDeviceOps` (ewald.h, pure virtual), implemented by `GfnffEwaldGPU`
(`cuda/ewald_gpu.cu`), owned by the wrapper, handed over by `GFNFF::setEwaldDeviceOps` to `EEQSolver` and `FFWorkspace`
(not owned there). Every sum has a fixed order (warp shuffles, block trees, slices summed in order, one-block final sums):
no atomic additions. A device failure (any call returns false) is a counted fallback (`reportGpuFallback`,
`gfnff_periodic_ewald_device`; `-gpu_strict` stops): the pcg restarts on the host cache, the Ewald energy takes the host
partitions. Fault injection `CURCUMA_GFNFF_EWALD_DEVICE_FAIL=prepare|potentials|evaluate`; profile lines
`CURCUMA_GFNFF_EWALD_PROFILE=1` (EWALDPROF) and `CURCUMA_GFNFF_STEPPHASE=1` (STEPPHASE_EEQ: pcg split). Plugin ABI 26 -> 27.
Multi-GPU with a cell stays off (unchanged); ROCm has no device Ewald (`createEwaldDevice()` returns nullptr, the cell is
refused there anyway).

### Acceptance (thresholds fixed beforehand; ctest `gfnff_gpu_pbc_ewald`, CUDA only, 77 without CUDA)

| check | result |
|---|---|
| (a) device Ewald vs CPU and vs the host-Ewald GPU hybrid, input geometry and every atom moved by <= 0.005 A (per-step solve); 20 H2O (FD settings), 60 H2O, H3O+(H2O)2 + 18 H2O (Q = +1), 970 H2O L 58 (pcg), 4 threads | charges <= 3.8e-15 e, energy <= 5.2e-14 Eh, gradient <= 7.9e-15 Eh/A, model virial <= 9.2e-14 Eh (limits 1e-9 e, 1e-9 Eh, 1e-9 Eh/A, 1e-7 Eh); device used (7 operator calls per solve at 2910 atoms), no fallback counted |
| (a) W1 boxes 432 / 864 / 1728 H2O (`box*_s1`, 4 threads, production settings) | vs CPU: charges <= 7.7e-15 e, energy <= 7.7e-12 Eh, gradient <= 4.2e-14 Eh/A, virial <= 6.5e-11 Eh; vs host-Ewald GPU: energy <= 9.1e-13, gradient <= 1.0e-15, virial <= 9.3e-13 |
| (a) ionic mixture Na2Cl2 + 380 H2O, L 48 Bohr (1144 atoms, pcg) | vs host-Ewald GPU: energy 4.3e-14 Eh, gradient 1.1e-15 Eh/A, virial 5.7e-14 Eh; vs CPU 3.4e-9 Eh / 3.5e-8 Eh/A / 1.1e-7 Eh (charges 5e-11 e): **pre-existing**, the same numbers with the A1 base build (D4 pair term 3.5e-9 Eh and its CN chain rule for the isolated ions; not Ewald), test limit there 1e-8 / 1e-7 / 2e-7 (TODO.md) |
| (b) FD gate, device Ewald: gradient vs Richardson FD of the total energy (h 2e-4 / 1e-4 A, 6 atoms near a face) and virial vs FD under isotropic and axial strain; 20 H2O and Q = +1, dense EEQ and forced pcg (`eeq_ewald_dense_max_atoms 0`, 423 / 558 device operator calls) | gradient <= 9.4e-10 Eh/A (limit 1e-8), virial <= 6.3e-9 Eh (limit 1e-7) |
| (c) second GPU instance, 970 H2O, two steps, 4 threads | energy, gradient, charges, virial bit for bit |
| (d) fault injection prepare / potentials / evaluate | energy equal to the host-Ewald hybrid within 0 / 2.8e-13 / 1.4e-13 Eh, fallback counted 2 / 1 / 1 times |
| no cell, 6 structures x gfnff/gfn1/gfn2 x 1/4 threads + 2 MD, `-sp -dump_gradient` against the A1 build (`gfnff-gpu-virial`, same options) | CPU 38/38 lines identical, CUDA 38/38 identical; SASS (`cuobjdump -sass`) of every function of the base plugin identical (names in anonymous namespaces carry the source path hash; the cub scan kernels now have a second instantiation) |
| CPU with a cell, 432 / 864 H2O `-sp` (pcg inside the set-up), new vs A1 build | energy and gradient md5 identical (the restructured CPU pcg is bitwise the old one) |
| CLI NPT 432 H2O (`box432_s1`, 20 steps, CSVR 300 K, Berendsen every 2 steps, 4 threads): CPU, A1 hybrid, new | 15 status lines: Epot (6 decimals), P (2), V (3), rho (5), L (5) identical in all three |
| existing ctests `gfnff_gpu_pbc_hybrid`, `gfnff_gpu_pbc_virial` (now on the device Ewald by default) | pass |
| full `ctest -j1` | `release/` 494 of 496 run, `release_cuda/` 778 of 780; the two failures in both are the known `cli_confsearch_05_seed_cold_bias`, `cli_deltaml_04_label_consistency` (red in the base build too) |

### Time per NPT step (ms; CLI, 4 host threads, CSVR, Berendsen every 2 steps, difference of 60- and 20-step runs, n = 3, median; RTX 5080; load average 0.6-3.3 from the runs themselves, nothing else running)

| box (`box*_s1`) | atoms | CPU | A1 hybrid (host EEQ + Ewald) | S3 (device Ewald + operator) | S3 vs A1 | S3 vs CPU |
|---|---|---|---|---|---|---|
| 432 H2O | 1296 | 47.0 (45.9-49.2) | 23.6 (21.4-23.7) | 12.6 (12.4-12.8) | 1.9x | 3.7x |
| 864 H2O | 2592 | 119.6 (118.0-123.8) | 63.0 (61.3-64.2) | 32.6 (32.5-32.7) | 1.9x | 3.7x |
| 1728 H2O | 5184 | 308.8 (308.1-309.1) | 175.8 (170.5-176.1) | 93.4 (93.2-93.4) | 1.9x | 3.3x |
| 3456 H2O | 10368 | 871.0 (870.0-871.5) | 514.2 (512.4-520.4) | 278.4 (278.1-279.4) | 1.8x | 3.1x |

Split at 3456 H2O (STEPPHASE, step without a barostat event): energy call 96 ms, of it EEQ 73 ms = device pair build 17 ms
(4.6 ms alone, the rest is contention with the device terms on other streams) + 7 operator calls 48 ms (6.8 ms each: structure
factors 3.1, per-atom pass 3.5, real space 0.16; FP64-bound, about half the FP64 peak of the card) + CN download 5.7 ms;
Ewald energy/gradient/virial 9 ms; the host share of the CG is about 1 ms. A step with a barostat event costs 396 ms: the
forced list rebuilds after the cell change (HB/XB 101 ms, other host lists ~200 ms; A1/N2 behaviour, unchanged here) are now
the largest item of an NPT step (about half of the 278 ms average).

### S4 (device-resident pcg): not done, by measurement

The CG loop itself (vectors, projection, block preconditioner, dots) costs about 1 ms of the 66 ms solve at 10368 atoms
(0.2 ms at 5184); moving it to the device would save that and seven 80 KB round trips. The solve is the device operator.
Levers that remain, each measured or estimated above: (1) the forced list rebuilds on barostat steps (~150 ms per step
averaged at 3456 H2O); (2) a larger real-space radius moves work out of the FP64-bound reciprocal sum (the device real space
costs 0.16 ms for 11 million entries): `-gfnff.ewald_r_cut` up to L/2 changes the energy at the 1e-7 Eh level (ewald_tolerance
is per term) and must be set for CPU comparisons too; at L/2 instead of 17.5 sqrt(L/40) the k-vector count falls by about
3.5x at 3456 H2O (estimate, not run); (3) a pair list with a skin instead of the per-step rebuild (4.6 ms alone).

### Not done / not tested

- ROCm (no device Ewald; B2/B3), multi-GPU with a cell, triclinic cells, energy-only MD on the device path (covered only
  by the FD energies), `ewald_kernel_cache_mb` (the device keeps every entry; memory 20 bytes per entry, 225 MB at
  10368 atoms), more than 10368 atoms, heavy elements near the r_g < L/2 limit, charged boxes under the ensemble charge
  model (CPU fallback as before), a long NPT run (20 steps here), C-rescale (same code path). Human production testing
  pending.
