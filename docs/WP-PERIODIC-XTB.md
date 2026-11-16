# WP: Periodic native GFN1/GFN2 (orthorhombic cell, Gamma point, CPU)

**Status (2026-11-15): P0 plan, P1 periodic classical terms, P2a/P2b periodic GFN1/GFN2 energy, P3 periodic gradient,
P4 virial/pressure/barostat, P5 integration and T threading done; G (CUDA port) in progress, see the G section.**
`-method gfn1|gfn2` with an orthorhombic unit cell gives the periodic Gamma-point energy, gradient and atomic virial
(gfn1: 3-7e-11 Eh, 8.6e-7 Eh/Bohr and 1.8e-6 Eh virial from tblite on the water dimer; gfn2: tblite's periodic multipole
Ewald is not alpha-consistent, curcuma uses the consistent form and differs from tblite by up to 1e-3 Eh at 8 A, see P2b;
`-xtb.pbc_multipole_tblite_factors true` reproduces tblite), so `-opt`, NVT/NVE and NPT `-md` with a cell run (SimpleMD
wraps fragments, -3 DOF, `-compute_pressure`, `-barostat berendsen|crescale`). Without a cell nothing changes (bitwise
gate in every stage).
Branch `feature/periodic-xtb`, base `feature/bulk` 957afab7. All work here is 🤖 AI-generated; stages are ⚙️ machine-tested
only when the stage section says so. No human production testing anywhere in this document.

Goal: the same periodic treatment GFN-FF got in [WP-PERIODIC-NONBONDED.md](WP-PERIODIC-NONBONDED.md) (A0-A6b, B0-B3, N1/N2,
V-LRC) for the native GFN1/GFN2, so that a water box without surface can be run with GFN2 (density, later NPT). Motive:
GFN2 water droplets in vacuum come out about 40-45 % too dense (to be verified; `GFN-FF-Wasser/06_gfn2/BERICHT.md`).

Reference: **tblite 0.4.0** (`external/tblite`, built in `release_tblite/`), which does periodic GFN1/GFN2 at the Gamma
point only (no k-point sampling; the C API takes `lattice` and `periodic`). Checked in P0: `tblite_periodic_reference`
(`test_cases/sqm_reference/tblite_periodic_reference.c`, build line in its header, linked against the static libraries of
`release_tblite/_deps`, no network) gives for the registry water dimer GFN2 vacuum -10.1346829493237554 Eh, cubic cell 30 A
-10.1346783862059819 (difference 4.6e-6 Eh: tin-foil boundary of a dipolar cell plus the Klopman-Ohno switch, see P0.4),
cell 8 A -10.1324596572071250; GFN1 cell 8 A -11.5308428959022216, gradient and virial are returned. tblite is therefore a
usable reference for energy, gradient and (P4) virial.

## P0 survey: where the geometry enters the native xTB (read 2026-11-14 on 957afab7)

All in `src/core/energy_calculators/qm_methods/`. Every pair loop below uses plain `r_ij = x_i - x_j`; none reads
`Mol::m_has_pbc`/`m_unit_cell` (`NativeXtbMethod::setMolecule` copies the `Mol` but passes no cell to `XTB`).

| Term | Function | Geometry dependence | tblite 0.4.0 periodic treatment (reference) |
|---|---|---|---|
| CN for H0 self-energies | `cn_exp` (GFN1) / `cn_gfn` (GFN2), `parameters/xtb_params_extra.hpp`, called from `XTB::computeCoordinationNumbers` (`xtb_h0.cpp`) | all pairs, cutoff 25 Bohr | sum over lattice images (j, T) != (i, 0) within 25 Bohr; self images count (`ncoord/type.f90:167-182`) |
| Overlap S, H0 | `XTB::getHamiltonianH0` (`xtb_h0.cpp:184`), dense nao x nao, every shell pair, no cutoff; on-atom block `h_factor = 1`, `S_diag = 1` | distance in `shpoly`, `hscale`, CGTO overlap | S(Gamma) = sum_T S(T), same for H0 (`xtb/h0.f90:196-290`): images of the same atom (T != 0) use the off-site formula with Z_i = Z_j; the T = 0 on-atom block separately (`:290-336`); adjacency list within `min(sqrt(2 intcut/min_alpha), 40)` Bohr, intcut 25 at accuracy 1 (`basis/type.f90:166-199`) |
| GFN2 dipole/quadrupole AO integrals | `XTB::setupMultipole` (`xtb_multipole.cpp:137`), origin at the column atom, traceless | same pairs as S | same image sum, `shift_operator` per image (`h0.f90`) |
| Isotropic Coulomb gamma (Klopman-Ohno) | `coulomb::build_gamma_matrix` (`xtb_coulomb.hpp:98`) via `XTB::buildGammaMatrix` | nsh x nsh, `1/(r^2 + gam^-2)^(1/2)` for every shell pair on different atoms | Ewald: real space over images of `fcut*KO + (1-fcut)/r - erf(a r)/r`, `fcut` switches KO to 1/r between 9 and 10 Bohr (`rcut = 10`, `fsmooth`), reciprocal `(4pi/V) sum_k exp(-k^2/4a^2)/k^2 cos(k.r)`, self `-2a/sqrt(pi)` on the diagonal, on-atom hardness added as in 0d; **no** background term `-pi/(a^2 V)`; Wigner-Seitz image weights (`coulomb/charge/effective.f90:300-447`) |
| Third order | `xtb_thirdorder.cpp`, `xtb_coulomb.hpp:188-279` | none (on-site) | unchanged |
| GFN2 multipole electrostatics (AES) | `XTB::buildMultipoleInteractionMatrices` (`xtb_multipole.cpp:68`): `amat_sd`, `amat_dd`, `amat_sq` nat x nat, damped `1/r^3`, `1/r^5` kernels, `mrad` from the GFN2 CN | all atom pairs, no cutoff | Ewald of the bare multipole kernels plus the damped-minus-bare remainder in real space, reciprocal `sin`/`cos` sums with k-vector tensors, self terms `-2/3 a^3/sqrt(pi)` (dipole) and `4/9 a^3/sqrt(pi)` (quadrupole trace, zero for traceless), conv = 100 eps (`coulomb/multipole.f90:478-640`) |
| Repulsion | `XTB::calcRepulsionEnergy` (`xtb_h0.cpp:369`), all pairs | `Z_i Z_j / r^rexp exp(-sqrt(a_i a_j) r^kexp)` | images within 25 Bohr incl. self images, each (i, j, T) once with 1/2 per atom (`repulsion/effective.f90:159-178`) |
| Halogen bond (GFN1) | `calcHalogenBondEnergy`, `buildHalogenBondList` (`xtb_h0.cpp:476-579`): X = Cl/Br/I/At, A = N/O/P/S, B = nearest neighbour of X | X-A pairs within 20 Bohr, B over the whole molecule | X-A pairs over images within 20 Bohr, B over atoms and images (`classical/halogen.f90:149-245`) |
| D3 (GFN1) | `D3ParameterGenerator::prepareForEnergyGradient` / `getEnergyAndGradient` (`ff_methods/d3param_generator.*`), own CN (exp counting, 25 Bohr), pair sum without cutoff, s9 = 0 | all pairs | s-dftd3 through tblite: CN 25, pair cutoff 50 Bohr (`disp/d3.f90:59`), images |
| D4 (GFN2) | `D4ParameterGenerator::GenerateParameters` (reference weights from the EN-weighted covalent CN, `dispersion/d4_ncoord.*`, 30 Bohr), `D4Evaluator::computeEnergyAndGradient` (pairs i < j, charge-scaled C6 from the SCF Mulliken charges), `computeATM` (triples, 25 Bohr), `addD4CovalentCNGradient`; in-SCF potential `computeD4PotentialDedq` | all pairs / triples | dftd4 through tblite: CN 30, pair cutoff 50, ATM 25 Bohr (`disp/d4.f90:82`), images in CN, pairs and triples |
| Gradient | `XTB::calculateGradient` (`xtb_gradient.cpp:51`): sections 1 repulsion, 1b XB, 2a/2b H0/Pulay (atom pairs iat < jat, `dS`, `shpoly`, multipole integral Pulay), 3 Coulomb (`dgamma/dr` per shell pair), 5 multipoles, 3b D4 (cached), 4 CN chain rule | every section is a plain pair loop | tblite `get_hamiltonian_gradient` over the adjacency list with `trans`, Ewald gradients (`get_multipole_gradient_3d`), dftd4/s-dftd3 image gradients |
| SCF guess / mixer / eigensolver / preconditioner / CPSCF | `seedEEQGuess` (non-periodic EEQ, guess only), Broyden, `solveEigen`, `xtb_precond.cpp`, `xtb_response.cpp` | guess and preconditioner only change the path, not the fixed point; CPSCF only with `d4_charge_source cpscf` (not default) | unchanged |
| Not combinable with a cell | implicit solvation (`m_solvation`), `-large_system_mode fragments|dc|sparse` (`FragmentScfDriver`), GPU wrappers (`XtbGpuComputationalMethod`, HIP, Vulkan) | | error (solvation, large-system) or counted CPU fallback (GPU, P5) |

## P0 decisions (physics and numerics)

1. **Gamma point only, orthorhombic cell, CPU.** The Gamma-point approximation sums the Bloch phases with k = 0: every AO
   interacts with every image of every other AO with phase 1. For a cell that is small against the electronic coherence
   length (a few molecules) the band dispersion is sampled by one point and the energy per cell carries an error that
   decays with the cell size; the test is the supercell convergence E(n x n x n)/n^3. tblite has the same limitation, so the
   agreement test (<= 1e-6 Eh) checks the implementation, the supercell test checks the physics. Triclinic cells are
   refused like in GFN-FF (`PBCUtils::makeOrthoCell`).
2. **Image sums, not minimum image, for every short-range term** (S/H0, CN, repulsion, XB, D3/D4 pairs and triples, the
   Klopman-Ohno remainder, the multipole remainder): the lattice vectors T with `|r_ij + T| < R_c` are summed explicitly,
   so `R_c < L/2` is **not** required; a self image (i, i, T != 0) counts. The cutoffs follow tblite 0.4.0 so that the
   reference comparison is a check of the code, not of a model choice: integral cutoff `min(sqrt(50/min_alpha), 40)` Bohr
   (`-xtb.pbc_integral_cutoff_bohr`, 0 = this rule), CN 25 (H0), repulsion 25, XB 20, D3 CN 25 / pairs 50, D4 CN 30 /
   pairs 50 / ATM 25 Bohr. Without a cell every one of these loops keeps its present (cutoff-free) form.
3. **Long-range Coulomb by Ewald summation**, reusing `curcuma::ewald::EwaldSum` (B0/B2): the point-charge Ewald matrix
   `A_pt(i, j)` (real space erfc(a r)/r under the minimum image with `r_cut <= L_min/2`, reciprocal, self, background) with
   `alpha`, `k_max` from `chooseEwaldParameters(cell, pbc_ewald_tolerance = 1e-10, r_cut)`. The Ewald real-space part is the
   one place where the minimum image is used; the tolerance makes the neglected images (erfc(a r_cut) < 1e-10) irrelevant
   and the result alpha-independent (test).
4. **Klopman-Ohno is not Gaussian: exact split into Ewald and a short-range remainder, following tblite.**
   `KO(r) = 1/sqrt(r^2 + gam^-2) = 1/r + [KO(r) - 1/r]`; the remainder behaves as `-gam^-2/(2 r^3)`, and a lattice sum of
   `r^-3` diverges logarithmically, so the periodic Klopman-Ohno **needs a cutoff by construction**. tblite switches the
   kernel to 1/r with the quintic `fsmooth` between 9 and 10 Bohr (`rcut = 10`); the remainder `fcut (KO - 1/r)` is then a
   short-range image sum (within 10 Bohr) and `gamma_pbc(is, js) = A_pt(at(is), at(js)) + sum_T fcut(r) [KO_gam(r) - 1/r]
   + on-atom terms as in 0d`. Consequence: the periodic energy of a cell larger than the molecule equals the vacuum energy
   only when every pair is below 9 Bohr (water dimer: yes); for larger molecules the KO -> 1/r switch beyond 9 Bohr is a
   model difference to the vacuum GFN2 of a few 1e-5 Eh per pair at 10 Bohr (3 % of 1/r). This is tblite's definition of
   periodic GFN-xTB and is kept as `-xtb.pbc_coulomb_switch_bohr` (default 10; change = no tblite agreement).
   **It is a model constant, not a convergence parameter**: energy and pressure of a liquid grow with ln(rc) without
   limit (measured in "Klopman-Ohno switch radius" below).
   The background term `-pi/(a^2 V)` (missing in tblite) is included so that charged cells are alpha-independent; a
   charged cell therefore differs from tblite by `-pi Q^2/(2 a^2 V)` (documented, tested neutral only).
5. **GFN2 multipoles with Ewald** (tblite's split: damped kernel minus Ewald-screened bare kernel in real space, bare
   kernel in reciprocal space, dipole self term). The self-consistent D4 potential and the third order are local in the
   charges and stay as they are. If the multipole Ewald does not fit into P2, GFN2 with a cell aborts with an error, never
   silently.
6. **No-cell path bit for bit**: every periodic branch is a separate function (`*Periodic`) chosen by `if (m_pbc)` at the
   call site, outside the arithmetic; the FMA-contraction trap of WP-PERIODIC-NONBONDED S1/N1 applies. Gate: `-sp` energy
   (16 digits) and md5 of `-dump_gradient` for gfn1 and gfn2 on several registry structures at 1 and 4 threads against
   the binary built from 957afab7 (kept outside the tree).
7. **Energy per term is kept**: `getEnergyDecomposition()` reports the same containers; the P1 tests read them.

## Stages and effort (estimates of 2026-11-14)

| Stage | Content | Estimate |
|---|---|---|
| P0 | survey, tblite check, plan (this section) | done |
| P1 | lattice-image helper in `XTB`; periodic CN, repulsion, XB, D3, D4 (CN, pairs, ATM) energies; tests: supercell scaling, vacuum limit, tblite classical terms where separable | 1 day; D4/ATM images are the largest part |
| P2a | GFN1 electronic: S/H0 image sums, gamma by Ewald + KO remainder, SCF unchanged; tblite <= 1e-6 Eh | 1 day |
| P2b | GFN2: periodic AO multipole integrals (image sum with origin shift), multipole Ewald matrices, D4 potential | 1-2 days |
| P3 | gradients of every term incl. images and Ewald (point-charge Ewald gradient of `EwaldSum`, KO remainder, multipole Ewald gradient port), FD gate of the assembled surface, tblite forces | 1-2 days |
| P4 | atomic virial (strain derivative at fixed fractional coordinates) for every term, `-compute_pressure`, barostat hookup; without it SimpleMD refuses gfn1/gfn2 + barostat | 1-2 days |
| P5 | `hasPeriodicEnergy()` for gfn1/gfn2 with a cell, SimpleMD wrap/DOF, GPU guard (counted CPU fallback), tests | 0.5 day |

Cost note for the water-box target (432 H2O, 1296 atoms, GFN2 nsh = 1728, nao = 2592): the dense Gamma-point matrices are
those of the vacuum calculation; the image sums add a factor of the number of images per pair within the cutoffs
(a few at L ~ 25 A); the Ewald matrix via `EwaldSum::matrix` is nat x nat (13 MB); the multipole reciprocal sums must be
factorised over k vectors (rank updates), a per-pair k loop would be O(nat^2 N_k). Measured numbers follow in P2/P3.

## P1 periodic classical terms (2026-11-14, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

**Scope.** `XTB::setPeriodicCell()` (orthorhombic, Angstrom in, Bohr inside; refuses triclinic cells and, for GFN2,
`d4_charge_source` other than `mulliken`) and the image-sum twins of the classical kernels in
`qm_methods/xtb_periodic.cpp`, each reached by `if (m_pbc)` at the entry of the cell-free kernel:
`computeCoordinationNumbersPeriodic` (H0 CN, 25 Bohr), `calcRepulsionEnergyPeriodic` (25), `calcHalogenBondEnergyPeriodic`
(20; X-A over images, B the nearest atom or image of X), `calcD3DispersionPeriodic` (GFN1: CN 25, pairs 50, CN chain rule over
images; `D3ParameterGenerator::prepareForEnergyGradientWithCN`, `getEnergyAndGradientPeriodic`), `calcD4DispersionPeriodic`
(GFN2: `computeD4CovalentCNPeriodic` 30 pushed into the generator, `D4Evaluator::computeEnergyAndGradientPeriodic` 50,
`computeATMPeriodic` 25, `addD4CovalentCNGradientPeriodic`), and the in-SCF `computeD4PotentialDedq` with images. Shared
helper `PBCUtils::LatticeImageList` / `forEachImage` in `src/tools/pbc_utils.h` (minimum image first, then the translations
`|n_a| <= (R + L_a/2)/L_a`; a self image (i, i, T != 0) counts; convention `d = x_i - x_j - T` as in tblite).
**Not yet wired to the CLI**: `NativeXtbMethod` does not forward the cell until the electronic part is periodic (P2), so
`-method gfn1|gfn2` with a unit cell is still the vacuum energy; the P1 code is reached only through the `XTB` object.

- **Counting rules** (dftd4/s-dftd3/tblite): pairs j <= i, every image, `dE/2` to i and (i != j) to j, so a self-image
  pair (i, i, T) counts half per translation (T and -T are both visited) and has no direct Cartesian gradient;
  `dE/dCN_i -= dC6_ii/dCN_i(first argument) edisp` once per translation. The D3 generator's `computeDC6DCN()` now also
  fills the diagonal (dC6_ii/dCN_i), read only by the periodic self-image sum; the cell-free loops use i < j.
- **Finding, dftd4 periodic ATM gradient**: `get_atm_dispersion_derivs` (dftd4 3.x, `damping/atm.f90`) applies the
  `triple_scale` to the energy and the virial but not to `dGij/dGik/dGjk`. In a molecule every triple has three distinct
  atoms (scale 1), with lattice images the self-image triples (iat = jat, 1/2; iat = jat = kat, 1/6) are double counted
  in the force. Found by the FD gate below (periodic ATM gradient 10 % off FD on the 8 A water-dimer cell, pair term
  exact); `computeATMPeriodic` scales the gradient with `triple`. Code reading of dftd4 plus the FD here; not confirmed
  against a dftd4 run. Energies are unaffected, so tblite energies remain the reference; tblite periodic *forces* would
  carry this error for small cells (to be checked with `tblite_periodic_reference --gradient` against its own FD in P3).

Tests (`ctest -R pbc_xtb_classical`, `test_cases/unit/test_pbc_xtb_classical.cpp`, registry water dimer, GFN1 and GFN2, test
charges O -0.6 / H +0.3 for the D4 zeta weighting):

| Check | Result |
|---|---|
| vacuum limit: 60 A cube vs no cell (repulsion, XB, CN, dispersion, dispersion gradient) | 0 / 0 / 0 / 0 / 7e-21 |
| 2x1x1 supercell (16 x 8 x 8 A) vs 8 A cube: repulsion and dispersion per cell, CN and dispersion gradient per atom | 1.4e-16, 2e-17, 1.1e-15, 3e-20 (gfn1), 1.4e-16, 2e-17, 1.1e-15, 1e-19 (gfn2) |
| second molecule shifted by (L, -L, 2L); dimer moved and wrapped atom-wise into [0, L) | repulsion 1.8e-16 / 4e-17, dispersion 1e-19 / 0, CN 4e-16 / 0, gradient 3e-20 / 2e-20 |
| repulsion in the 8 A cube vs a 13^3-image brute-force sum written in the test | 1.4e-17 (gfn1), 4.2e-17 (gfn2) |
| dispersion gradient (pairs + ATM + CN chain rule) vs central FD, h 1e-4 A, fixed charges | no cell 4.3e-13 / 5.4e-13 Eh/Bohr; 8 A cube 5.6e-13 / 5.9e-13 (max gradient 7.6e-6 / 2.9e-5) |
| images are seen in the 8 A cube | E_disp -7.254e-4 vs free -5.475e-4 (gfn1), -9.006e-4 vs -7.384e-4 (gfn2); CN(O) 1.98897 vs 1.98881; GFN1 repulsion identical to 12 digits (exp(-a r^1.5) at 15 Bohr), GFN2 2e-8 |

- **Without a cell**: `-sp` energy and `-dump_gradient` md5 on benzene, caffeine, water8, hydronium-water2 (+1), methanol,
  urea3-water8 for gfn1 and gfn2 at 1 and 4 threads: 24/24 identical to the binary built from 957afab7.
- **Not tested**: the halogen bond with a cell (no registry structure with Br/I; code mirrors tblite `get_xbond_list` and the
  vacuum kernel, vacuum limit on water is trivially 0), elements beyond H/O, charged cells, the self-image half weights for
  cells so small that `L < 2 r_cov` (bond to the own image), threads (the periodic kernels are serial), tblite agreement of
  the individual terms (not separable through the C API; total energies in P2).

## P2a periodic GFN1 energy: S/H0 image sums and Ewald Coulomb (2026-11-14, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

- **S and H0 at the Gamma point** (`XTB::getHamiltonianH0Periodic`, chosen by `if (m_pbc)` at the entry of
  `getHamiltonianH0`): every shell pair over the lattice images within `pbcIntegralCutoffBohr()` = tblite's
  `min(sqrt(2 * 25 / min_alpha), 40)` Bohr (`-xtb.pbc_integral_cutoff_bohr`, 0 = this rule); the image block uses the
  off-site `hscale * shpoly(r)` of the image distance (also for the self images, Z_i = Z_j); the T = 0 on-atom block is the
  vacuum one (`h_factor` 1, exact diagonal). Same `parallelStripes` over the first shell as the vacuum kernel.
- **gamma by Ewald summation** (`XTB::buildGammaMatrixPeriodic`): `gamma(is, js) = A_pt(at_i, at_j) + sum_T fcut(r)
  [KO_gam(r) - 1/r] + on-atom hardness`, `A_pt = EwaldSum::matrix(positions)` (point charges: erfc(a r)/r under the
  minimum image, reciprocal, self, background `-pi/(a^2 V)`), `alpha`/`k_max` from `chooseEwaldParameters(cell,
  pbc_ewald_tolerance = 1e-10, r_cut = L_min/2)`; the Klopman-Ohno remainder over images within
  `pbc_coulomb_switch_bohr` = 10 with tblite's quintic `fsmooth` over the last Bohr (P0 decision 4). Reported at
  verbosity 3 as alpha and the number of k vectors.
- **Wrapper** (`NativeXtbMethod::setMolecule`): `Mol::m_has_pbc` -> `XTB::setPeriodicCell(m_unit_cell)`; errors (never a
  silent vacuum energy) for gfn2 (P2b pending), implicit solvation, a non-orthorhombic cell, `d4_charge_source` other
  than `mulliken`; `-large_system_mode` keeps its own driver (no cell, documented as not periodic). `hasPeriodicEnergy()`
  is true with a cell. A gradient request with a cell is refused until P3 (`calculateEnergy(true)` sets an error).
- **tblite agreement** (water dimer `water-dimer.roo2.8`, `tblite_periodic_reference gfn1 ... L L L`, tblite 0.4.0
  accuracy 1, curcuma `-sp -method gfn1`, cubic cells):

| L (A) | E curcuma (Eh) | E tblite (Eh) | dE (Eh) |
|---|---|---|---|
| 8 | -11.530842895860 | -11.530842895902 | +4.2e-11 |
| 10 | -11.531571559511 | -11.531571559541 | +3.0e-11 |
| 12 | -11.531722031160 | -11.531722031225 | +6.6e-11 |
| 30 | -11.531825251411 | -11.531825251444 | +3.3e-11 |

  tblite's own run-to-run scatter is 1e-14 (OpenMP reductions); tblite converges the energy to 1e-6 * accuracy and the
  charges to 2e-5 * accuracy (`xtb/singlepoint.f90`), at `--accuracy 0.001` its L = 8 energy moves by 5e-11, so the
  3-7e-11 above is the SCF convergence of both codes. Charges: curcuma at `scf_threshold` 1e-9 against tblite at
  accuracy 0.001 agree to 2.5e-10 e (at accuracy 1 tblite's charges are 1e-5 off their converged values). The 8 A cell
  differs from vacuum by 1.0e-3 Eh, the 30 A cell by 1.9e-5 (tin-foil dipole term of the dipolar cell; the Klopman-Ohno
  switch does not act below 9 Bohr).
- **Without a cell**: 24/24 `-sp`/`-dump_gradient` md5 identical to the 957afab7 binary (6 structures x gfn1/gfn2 x 1/4
  threads).
- ctest `pbc_xtb_tblite` (`test_cases/unit/test_pbc_xtb_tblite.cpp`): the three energies 8/10/12 A against the tblite
  values above (tolerance 1e-9 Eh), charges against tblite (1e-6 e), `hasPeriodicEnergy()`, gfn2 + cell refused,
  gradient request refused.
- **Not tested**: charged cells (tblite has no background term, see P0.4), open shell, d-shell elements, Ewald r_cut other
  than L_min/2 (`-xtb.pbc_ewald_r_cut_bohr`), threads > 1 with a cell (the Ewald matrix and the KO remainder are serial),
  cells so small that S(Gamma) becomes ill-conditioned (Cholesky fallback path), anything beyond H/O.

## P2b periodic GFN2 energy: multipole electrostatics by Ewald summation (2026-11-14, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

- **AO dipole/quadrupole integrals** (`XTB::setupMultipolePeriodic`, chosen by `if (m_pbc)` in `setupMultipole`): the
  same image sum as S/H0 (integral cutoff); every image block is shifted to the image position of its column atom
  (origin at atom(nu) + T, the periodic form of tblite's `shift_operator`), the traceless transform is applied to the
  sum. The damping radii `mrad` use the periodic GFN2 CN (`computeCoordinationNumbers()`).
- **Interaction matrices** (`buildMultipoleInteractionMatricesPeriodic`, port of tblite `get_multipole_matrix_3d` with
  tblite's Ewald parameters: `search_alpha` with the multipole terms, `search_cutoff` at conv = 100 eps, box-shaped
  translation lists; `alpha` 0.156 at L = 8 and 10 A, 0.141 at 12 A, 0.070 at 30 A): per atom pair the direct sum
  of damped kernel minus Ewald-screened bare kernel over the translation box, the reciprocal sums, the self terms.
  The direct box also truncates the damping remainder `(fdmp - 1) r^-3` (~ r^-6), which is why tblite's cutoffs are
  kept (a different direct cutoff is a different, converging, value; conv 1e-10 instead of 1.5e-6 moves the L = 12
  energy by 6e-8 Eh). Matrices are stored (anti)symmetric by construction (pairs j <= i mirrored).
- **Finding: tblite 0.4.0's periodic multipole matrices are not Ewald-consistent.** In `get_amat_sdq_rec_3d`
  (`coulomb/multipole.f90`) the charge-dipole reciprocal term is `2 k sin(k.r) w_k` and the charge-quadrupole term
  `+k_a k_b cos(k.r) w_k` over the full k set. The kernels are derivatives of the Coulomb kernel: charge-dipole
  `v/r^3 = -grad(1/r)`, dipole-dipole `delta/r^3 - 3 v v/r^5 = -grad grad(1/r)`, charge-quadrupole (traceless)
  `v_a v_b/r^5 = (1/3) d_a d_b (1/r)`. With `phi_rec = (4 pi/V) sum_k e^{-k^2/4a^2}/k^2 cos(k.r)` the reciprocal
  parts are `+k sin w` (factor 1), `+k_a k_b cos w` (as in tblite) and `-(1/3) k_a k_b cos w`; tblite's direct
  parts are exactly the matching `-grad erfc(ar)/r`, `-grad grad`, `(1/3) grad grad` (checked term by term, the
  `e1`/`e2` expressions). Measured with the probe `CURCUMA_PBC_MP_DEBUG` (rebuild at 1.5 alpha, L = 12): with
  tblite's factors the charge-dipole block changes by 1.1e-2 (of 7.3e-2) and the charge-quadrupole block by 6.6e-3
  (of 3.3e-2); with the factors 1 and -1/3 all three blocks change by <= 3.3e-7 (the conv level). Energies: with
  tblite's factors the port reproduces tblite at L = 8 to 2.1e-10 Eh (so the port itself is exact), with the
  consistent factors the total energy is alpha independent (0.14 vs 0.3: 5.5e-8 Eh at L = 12) and converges
  smoothly to vacuum; tblite's own L series does not (table). Decision: **the consistent factors are used**
  (`CURCUMA_PBC_MP_SD_FACTOR`, `CURCUMA_PBC_MP_SQ_FACTOR`, `CURCUMA_PBC_MP_ALPHA`, `CURCUMA_PBC_MP_CONV` are kept as
  diagnostic environment switches, defaults 1, -1/3, tblite alpha, 100 eps). Consequence: **periodic GFN2 energies
  differ from tblite 0.4.0 by up to 1e-3 Eh at L = 8 A** (3.6e-5 at 16 A, 5e-5 at 20 A, 2.9e-6 at 30 A) and the
  stop rule of the task (1e-5 against tblite) is formally violated with the cause located in the reference; not
  confirmed against a corrected tblite or an independent periodic GFN2 code. The dipole-dipole constant `fac/6` (a
  `2 pi/(3V) |M|^2` term, tblite convention) and the self terms are kept as in tblite.

| L (A) | E curcuma consistent (Eh) | E tblite 0.4.0 (Eh) | curcuma - vacuum | tblite - vacuum | GFN1 periodic - vacuum (tblite = curcuma) |
|---|---|---|---|---|---|
| 8 | -10.133473316 | -10.132459657 | 1.21e-3 | 2.22e-3 | 9.8e-4 |
| 10 | -10.134288580 | -10.133111595 | 3.9e-4 | 1.57e-3 | 2.5e-4 |
| 12 | -10.134522049 | -10.133359079 | 1.61e-4 | 1.32e-3 | 1.0e-4 |
| 16 | -10.134644660 | -10.134562075 | 3.8e-5 | 1.2e-4 | |
| 20 | -10.134670320 | -10.134618108 | 1.3e-5 | 6.5e-5 | 8.4e-6 |
| 30 | -10.134681289 | -10.134678386 | 1.7e-6 | 4.6e-6 | 1.1e-6 |

  Vacuum: GFN2 -10.134682949 (curcuma = tblite -10.1346829493), GFN1 -11.531826377. The periodic-minus-vacuum
  difference of a dipolar cell is dominated by the tin-foil dipole term and vanishes with the cell; the consistent
  GFN2 series follows GFN1's, tblite's GFN2 series does not.
- **Tests** (`ctest -R pbc_xtb_tblite`): gfn1 as in P2a; gfn2 (a) tblite factors reproduce tblite at L = 8 (1e-9),
  (b) alpha 0.14 vs 0.3 at L = 12 within 5e-7, (c) L = 30 within 5e-6 of vacuum, (d) regression values L = 8 / 12,
  charge neutrality. `pbc_xtb_classical` unchanged. Without a cell 24/24 bitwise (6 structures x gfn1/gfn2 x 1/4
  threads).
- **Cost**: the reciprocal part is a per-pair k loop, O(nat^2 N_k) (342 k vectors at L = 8-12 A; 1296 atoms x 342 is
  5.7e8 trig evaluations per geometry, a few seconds); a structure-factor factorisation is open (P-perf).
- **Not tested**: d-shell elements with a cell (the d block uses `sphericalMultipoleBlock` per image), charged and
  open-shell cells, cells smaller than 8 A, more than two molecules, threads > 1 with a cell (matrices built with
  `parallelStripes`, not compared across thread counts).

## P3 periodic gradient and FD gate (2026-11-14, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

- `XTB::calculateGradientPeriodic()` (chosen by `if (m_pbc)` at the entry of `calculateGradient`), every section of the
  vacuum gradient over the lattice images: repulsion (pairs j < i, images within 25 Bohr), GFN1 halogen bond
  (`addHalogenBondGradientPeriodic`, B over atoms and images), H0/Pulay per image with the image's own overlap block
  (the shpoly and dE/dCN terms need `S_T`, not the Gamma-point sum) and the multipole-integral Pulay term at the image
  position (origin shift `dR = x_j + T - x_i`), the self images (i, i, T != 0) only through dE/dCN (their Cartesian
  derivative is zero), isotropic Coulomb as `EwaldSum::gradient(positions, q_at)` (point charges, the same
  alpha/r_cut/k_max as the energy) plus the Klopman-Ohno remainder `d/dr [fcut (KO - 1/r)]` over images within the
  switch radius, GFN2 multipoles (`addMultipoleGradientPeriodic`: tblite `get_multipole_gradient_3d` with the P2b
  reciprocal factors 1 and -1/3, direct part verbatim, `dE/d(mrad)` into dE/dCN), the cached periodic dispersion
  gradient (D3/D4 of P1, CN chain rule already folded), and the CN chain rule of the H0 coordination number over
  images. `NativeXtbMethod` no longer refuses a gradient with a cell.
- **Not implemented with a cell**: d-shell elements (`setPeriodicCell` refuses them: the periodic gradient of the d
  multipole integrals would need the per-image spherical blocks); threads > 1 only in the H0/Pulay stripes.
- **Gate** (`ctest -R pbc_xtb_gradient`, `test_cases/unit/test_pbc_xtb_gradient.cpp`; water dimer, `scf_threshold`
  1e-9, central differences with h = 1e-4 A in Angstrom through `EnergyCalculator`; CLI gate
  `-dump_gradient` vs `-sp` differences with the same h): max |analytic - FD| 8.1e-9 (gfn1) / 8.9e-9 (gfn2) Eh/A at
  L = 8 A, 8.1e-9 / 8.9e-9 at 12 A, against max |g| 3.1e-2 / 2.8e-2; forces sum to zero (6e-16); GFN1 forces at L = 8
  within 8.6e-7 Eh/Bohr of tblite's (tblite accuracy 1, its charges converged to 2e-5). CLI gate gfn1 8.1e-9, gfn2
  8.9e-9 Eh/A. Without a cell 24/24 bitwise (6 structures x gfn1/gfn2 x 1/4 threads) to the 957afab7 binary.
- **Not compared**: GFN2 forces with tblite (tblite's periodic multipole gradient carries the factor-2 / +1 terms of
  P2b and the dftd4 ATM `triple_scale` issue of P1); forces of charged, open-shell or d-shell cells; thread-count
  independence of the periodic gradient.

## P4 atomic virial, pressure and barostat (2026-11-15, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

- **Definition**: Xi_ab = -dE/d(eps_ab) at fixed fractional coordinates (positions and cell scaled together), as
  `GFNFF::virial()` (N1); P = (2 K + tr Xi)/(3 V). Computed in `XTB::calculateGradientPeriodic()` when
  `setWantVirial(true)` (one gradient call; `virialValid()`), symmetrised.
- **Rule used everywhere**: a term that depends on an image pair vector d = x_i - x_j - T contributes -(dE/dd) (x) d,
  i.e. every per-image gradient vector of P3 contracted with its vector; **the self images (i, i, T != 0) now
  enter** (|T| scales with the cell) although their Cartesian derivative is zero: repulsion and D3/D4 pairs with weight
  1/2 per translation, the H0/Pulay self-image blocks with weight 1/2 (the loop visits both AO orderings; P3 had added
  their dE/dCN with weight 1, corrected here, 1e-7-level), the Klopman-Ohno remainder of on-site shell pairs (1/2 for
  is = js), the CN chain rules of the H0 CN, the D3 CN and the D4 CN (dE/dCN_i f'(r) d (x) d / r once per translation),
  the ATM through its three edge vectors (dG_ij (x) v_ij + ...), the halogen bond through its three vectors.
  Point-charge Ewald: `EwaldSum::realSpaceVirialRange`, `reciprocalVirialRange`, `backgroundVirial` (N1) at q_at.
  Multipole Ewald (`addMultipoleGradientPeriodic`): direct images -g (x) vec (self 1/2); reciprocal with k_a ->
  k_a - eps_ba k_b, w_k -> w_k [1 - tr eps + k_a k_b eps_ab (2/k^2 + 1/(2 alpha^2))], k.r invariant: charge-dipole
  E = sd_f sin w A (A' = -(mu_i,a q_j - mu_j,a q_i) k_b), dipole-dipole E = DDw cos w (k.mu_i)(k.mu_j), charge-quadrupole
  E = sq_f cos w B (B' = -2 [q_i (theta_j k)_a + q_j (theta_i k)_a] k_b), the 1/V constant fac/6 of the dipole-dipole
  block (-(fac/6) mu_i.mu_j delta_ab), alpha held fixed (the energy is alpha independent to conv). The strain
  derivative of the AO multipole integrals is the P3 Pulay term at the image position contracted with d (translation
  invariance: d/dd = d/dA at fixed B). Third order and on-site terms: no strain dependence.
- **Interfaces**: `XTB::setWantVirial/virial/virialValid`, `XTB::updateCell` (= `setPeriodicCell`: image lists, Ewald
  parameters and boxes are rebuilt from the cell at every call, nothing frozen, no cutoff limit);
  `EnergyCalculator::setWantVirial/virialValid/virial/largestActiveCutoffBohr` dispatch to GFN-FF (`gfnffInstance()`)
  or the native xTB engine through the new virtual `ComputationalMethod::xtbInstance()` (G1: a `dynamic_cast` to
  `NativeXtbMethod` does not resolve for the GPU adapters compiled into the plugins; plugin ABI 7 -> 8), `setUnitCell`
  accepts gfn1/gfn2; SimpleMD `-compute_pressure` and `-barostat berendsen|crescale` accept a periodic gfn1/gfn2
  (`readVirial` via `EnergyCalculator`), the GFN-FF cutoff freeze stays GFN-FF only.
- **Cell changes (barostat)**: the image lists adapt, the point-charge Ewald uses r_cut = L_min/2 with alpha, k_max
  from the tolerance (smooth in L), the multipole Ewald uses tblite's alpha search and box counts (`ceil(cut/L)`), which
  jump with L at the conv level (1.5e-6 in matrix elements, ~1e-7 Eh): continuous enough for MD, visible in a
  cell-scaling FD at h = 1e-5 only through the energy, not through the analytic virial (which is at fixed alpha).
- **Gate** (`ctest -R pbc_xtb_virial`, `test_cases/unit/test_pbc_xtb_virial.cpp`; water dimer, `scf_threshold` 1e-9,
  central differences of the energy under eps_xx, eps_yy, eps_zz and the isotropic strain, h = 1e-5, positions and
  cell scaled): gfn1 L = 8 A: 4e-11 / 1.2e-10 / 3e-11 Eh, trace 4.9e-10; gfn2 L = 8 A: 2.5e-10 / 9e-11 / 2.7e-10,
  trace 1.2e-9; gfn2 L = 12 A: 1.6e-10 / 1.9e-10 / 1.2e-11, trace 6.5e-10 (target 1e-7). Term isolation (probes
  `dispersionVirialNow`, `multipoleEnergyProbe`/`multipoleVirialProbe` at fixed moments, radii and alpha): D3/D4 incl.
  ATM 1e-13, dipole-dipole and charge-quadrupole 1e-13, charge-dipole 1e-13 after a sign fix found by this probe (the
  reciprocal charge-dipole energy is +sd_f sin(k.v0) w A with v0 = x_j - x_i). gfn1 virial against tblite's sigma
  (`virial_eh`, Xi = -sigma): max 1.8e-6 Eh (tblite accuracy 1); gfn2 not compared (P2b). Without a cell 24/24 bitwise.
- **Smoke**: `-md w27 (27 H2O lattice, 9.3 A) -method gfn1 -barostat berendsen -compute_pressure true -thermostat
  berendsen`, 40 steps: runs, pressure and density columns, cell updated (n = 1, no physics claim).
- **tblite factors as a registered switch** (`-xtb.pbc_multipole_tblite_factors`, Bool, default false): true sets the
  reciprocal charge-dipole / charge-quadrupole factors to tblite 0.4.0's 2 / +1 (energy, gradient and virial; the env
  variables `CURCUMA_PBC_MP_SD_FACTOR`/`SQ_FACTOR` still override both); `pbc_xtb_tblite` check (a) uses it. The
  default stays the Ewald-consistent 1 / -1/3 (operator decision pending, P2b).
- **Not tested**: off-diagonal strains (a sheared cell is not orthorhombic and is refused), the virial of charged or
  open-shell cells, d-shell elements (refused), thread-count independence of the virial, NPT densities (a water box
  run is the operator's campaign), the barostat cell-limit path (no limit for xTB).

## P5 integration (2026-11-14, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

- `NativeXtbMethod::hasPeriodicEnergy()` is true once the engine holds a cell; SimpleMD then applies the A6/A6b rules
  (fragment wrap into the cell, translation-only COM removal, -3 DOF). Smoke run: `-md dimer_L8.xyz -method gfn1
  -MaxTime 10 -dt 0.5 -thermostat none`: `COM/rot removal (mode 1): -3 DOF | effective DOF = 15`, 20 steps, no error
  (energy conservation not evaluated, n = 1).
- GPU with a cell: `EnergyCalculator::setMolecule()` swaps a `-gpu cuda|rocm|vulkan` gfn1/gfn2 request with a cell to
  the CPU method with a counted GPU fallback (`native xTB: periodic cell set, GPU path has no lattice images,
  calculation on the CPU`, same pattern as the GFN-FF A3b guard); `XTB::Calculation()` additionally refuses a cell
  with an attached device backend (`setHardError`). Checked only in a build without the plugin (the plugin fallback
  fires first, the periodic guard second, energy = CPU periodic energy, test in `pbc_xtb_tblite`); not run on a GPU.
- External methods (`xtb-*`, `tblite-*`, Ulysses) unchanged (a cell is still ignored there, A6b rules).

## Water-box cost (2026-11-14, one shared 32-core machine, `-threads 4`, `-sp` with `-dump_gradient`, wall incl. start-up)

Cubic water lattices at 0.997 g/cm^3 (scratch `make_water_box.py`, random orientations, 3.10 A spacing; not liquids,
timing only), energy + gradient:

| box | atoms | L (A) | gfn1 (s) | gfn2 (s) | gfn2 split (verbosity 3) |
|---|---|---|---|---|---|
| 3^3 = 27 H2O | 81 | 9.3 | 0.15 | 0.38 | |
| 4^3 = 64 H2O | 192 | 12.4 | 0.30 | 1.1 | |
| 6^3 = 216 H2O | 648 | 18.6 | 2.5 | 11.0 | setup 0.45 (gamma Ewald 0.12, multipoles 0.25), SCF 10 iterations 1.3, post-SCF 8.0 of which periodic D4 (pairs + ATM over images, serial) 7.98 |

The periodic D4 (`computeEnergyAndGradientPeriodic`, `computeATMPeriodic`) is serial and dominates; for 432 H2O
(1296 atoms, L ~ 23.5 A) a factor 4 on the N^2 parts gives an estimate of 35-45 s per energy+gradient at 4 threads,
i.e. a GFN2 NVT step of about half a minute until the periodic D4 is threaded and the ATM gets a neighbour list (the
vacuum `computeATM` has both). GFN1 (D3, no ATM) is 4-5x cheaper. NPT needs P4.

## Follow-up stages (operator decision of 2026-11-14): T threading, P4 virial, G CUDA port

Interim status 2026-11-15: **T done** (below), **P4 done** (section above), **G0 design written, G1-G7 see the G
status below**. The tblite multipole finding stays open; the factors remain switchable
(`-xtb.pbc_multipole_tblite_factors`, env `CURCUMA_PBC_MP_SD_FACTOR`/`SQ_FACTOR`, default consistent = Ewald-correct).

### T threading of the periodic loops (2026-11-15, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

- **Scheme** (deterministic reduction): the periodic D4 pair term, the periodic ATM (`D4Evaluator::runBlocks`) and the
  periodic multipole gradient (`addMultipoleGradientPeriodic`) split their row loop into `NB = min(nat, 32)` fixed blocks
  (`iat % NB == blk`), each block has its own energy/gradient/dE/dCN/dE/dq accumulators, and the blocks are reduced in
  block order. The block count depends only on the atom count, so these terms give the same bits for 1, 4 or 8
  threads (and run to run); worker t takes the blocks t, t + T, ... on the `CxxThreadPool` of the evaluator /
  `parallelStripes`, gated by `effectiveIntraThreads` (serial under molecule-level parallelism and for small systems,
  docs/SQM_THREADING.md). The in-SCF periodic dE/dq (`computeD4PotentialDedq`) uses the same blocks.
- **ATM neighbour list**: the images of every partner are collected once per `iat` (O(N) per atom) instead of once per
  (iat, jat, kat) triple; the triple loop runs over the partners with images in the same (iat, jat, kat) and image
  order, so the ATM sum is unchanged bit for bit (the plain O(N^3) image collection was the hot spot: 8.0 s of 11.0 s
  at 216 H2O).
- **Unchanged (still serial or per-thread partials)**: periodic CN and D3 (0.26 s at 216 H2O/8 threads, GFN1 only),
  the Klopman-Ohno remainder of gamma, the AO multipole-integral image sum (`parallelStripes`, disjoint rows), the
  H0/Pulay gradient (`parallelStripes` with per-thread partials as in vacuum: thread-count dependent at the last bit),
  the reciprocal multipole matrices (`parallelStripes` over `iat`, every element written once: deterministic).
- **Remaining thread dependence** of the total: energy 2e-12 Eh, gradient 1e-7 Eh/A between 1 and 4/8 threads (w64
  gfn2, w216 gfn1), the same for gfn1 (no D4 blocks involved), i.e. from the SCF (MKL eigensolver threads) and the
  per-thread partials of the H0/Pulay gradient as in vacuum (root CLAUDE.md: bit-identical only at a fixed thread count).
- **Timing** (water lattices as above, `-sp` + `-dump_gradient`, wall incl. start-up, s):

| box | method | threads | before T | after T |
|---|---|---|---|---|
| 64 H2O | gfn1 | 1 / 4 / 8 | 0.66 / 0.28 / 0.22 | 0.70 / 0.33 / 0.28 |
| 64 H2O | gfn2 | 1 / 4 / 8 | 1.42 / 1.10 / 1.07 | 1.45 / 0.60 / 0.64 |
| 216 H2O | gfn1 | 1 / 4 / 8 | 6.6 / 2.45 / 1.88 | 6.9 / 2.65 / 2.06 |
| 216 H2O | gfn2 | 1 / 4 / 8 | 13.6 / 11.0 / 10.7 | 7.9 / 3.7 / 2.5 |

  216 H2O gfn2 at 8 threads: setup 0.37, SCF (10 it.) 1.14, post-SCF 0.71 (D4 pairs 0.016, ATM 0.61, CN fold 0.04),
  total 2.2 s of the 2.5 s wall. Estimate for 432 H2O: 6-10 s per energy+gradient at 8 threads (the SCF eigensolve
  grows as nao^3). The gfn1 numbers are the start-up and SCF (the 0.1-0.2 s differences are machine noise, shared node).
- Gates: `pbc_xtb_classical`, `pbc_xtb_gradient` (FD 8e-9 Eh/A) and `pbc_xtb_tblite` green; the 216 H2O gfn2 gradient
  before/after the ATM restructuring agrees to 7e-13 Eh and 6e-8 Eh/A (both runs at 8 threads, i.e. within the
  run-to-run scatter of the threaded SCF; a 1-thread comparison of the old binary was not kept); without a cell 24/24
  bitwise to 957afab7.

### G0 design: CUDA port of the periodic gfn1/gfn2 paths (2026-11-15, Claude Generated; design only, nothing ported)

What runs on the device today (CUDA plugin `libcurcuma_cuda.so`, `XtbGpuComputationalMethod` in `xtb_gpu_adapter.h`,
`XtbGpuContext` in `qm_methods/cuda/`, stages 0-6 of [SQM_GPU.md](SQM_GPU.md)): S/H0/CN/gamma/multipole integrals
(`beginComputed`), the resident GFN1 and GFN2 SCF loops (Fock, eigensolve, populations, Broyden, in-SCF D4 and
multipole potentials), the nuclear gradient (repulsion, H0/Pulay, Coulomb, multipole Pulay), post-SCF D4 (pairs, ATM,
CN fold) and D3; all **without lattice images or Ewald sums**. With a cell the P5 guard swaps to the CPU (counted
fallback) and `XTB::Calculation()` refuses an attached backend; the GPU is used by the periodic path nowhere.

| G sub-stage | device work | host today | estimate | note |
|---|---|---|---|---|
| G1 | image-summed S/H0 and multipole AO integrals: the device integral kernels get a per-pair translation list (`LatticeImageList` upload) and accumulate `sum_T` into the same nao x nao matrices; the on-atom T = 0 block as today | `getHamiltonianH0Periodic`, `setupMultipolePeriodic` (0.03 + 0.17 s at 216 H2O) | 2-3 days | the SCF loop then runs unchanged on the device (Gamma point = same dense matrices) |
| G2 | periodic gamma: point-charge Ewald matrix on the device (real pairs under MIC, reciprocal as structure-factor rank updates `C W C^T + S W S^T`, self/background) plus the Klopman-Ohno remainder over images | `buildGammaMatrixPeriodic` (0.13 s) | 1-2 days | reuse the B0/B2 formulas; `EwaldSum::matrix` is the reference |
| G3 | periodic multipole interaction matrices (direct box + reciprocal with the P2b factors, tblite alpha/cutoffs) | `buildMultipoleInteractionMatricesPeriodic` (0.25 s, O(nat^2 N_k)) | 1-2 days | the k loop maps well to one thread per (pair, k block) |
| G4 | periodic D4/D3 (pairs, ATM, CN fold) with images; the existing device D4 kernels take a translation list | `calcD{3,4}DispersionPeriodic` (0.6 s at 8 threads) | 1-2 days | self-image weights 1/2 and the ATM `triple_scale` in the gradient as on the CPU |
| G5 | periodic gradient: image loops of H0/Pulay (incl. multipole Pulay at the image position), repulsion, KO remainder, Ewald gradient (device), multipole gradient, CN chain | `calculateGradientPeriodic` (~1 s) | 2-3 days | the self-image dE/dCN terms and the 1/2 weights of P4 |
| G6 | virial (P4 formulas: -G (x) d per image, Ewald/multipole strain terms) | `m_virial` accumulation | 1 day | after G5 |
| G7 | lift the P5 guard for CUDA only when G1-G6 exist; ROCm/Vulkan keep the CPU fallback (not ported, documented); every new `GpuScfBackend` virtual raises `CURCUMA_GPU_PLUGIN_ABI_VERSION` (S0 pattern) | | 0.5 day | parity gate: CPU vs CUDA energy/gradient/virial <= 1e-8 Eh resp. Eh/A on the dimer cells and a water lattice; bitwise CUDA run to run ([GPU_REPRODUCIBILITY.md](GPU_REPRODUCIBILITY.md)) |

Order of benefit: G1 (enables the resident SCF with a cell, the SCF is 1.1 of 2.2 s at 216 H2O and grows as nao^3)
before G4/G5 (the CPU image loops are already threaded), G2/G3 last (cheap on the CPU). Build: `release_cuda/` of the
main checkout holds a CUDA build (`-DUSE_CUDA_XTB=ON -DCMAKE_CUDA_ARCHITECTURES=120`), the RTX 5080 was idle (0 %)
on 2026-11-15; no GPU run was made in this session. **Status: G0 written; G1 done in the reduced form below
(host-built periodic matrices, device-resident SCF); of G2-G6 done: G3 (multipole matrices), G4a (ATM), G4b (in-SCF
D4 dE/dq), G5a (multipole gradient + virial pair loop), see the sections below; G1 device integral kernels, G2, the rest
of G4/G5/G6 open; G7 lifted for CUDA in the G1 form.**

### G1 (reduced): device-resident SCF with host-built periodic matrices, CUDA only (2026-11-15, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

Revision of the G0 design after reading the code: the `GpuScfBackend` already has an upload path next to the
device integral build (`begin(H0, S, L)`, `beginMultipole(dp_int, qp_int)`, Stage 2 of [SQM_GPU.md](SQM_GPU.md)) and a
host-driven loop in which the host assembles the potential and the device does Fock/eigensolve/density/moments. With a
cell this path needs no new device kernel: the host builds the periodic S/H0/L (image sums), the Ewald gamma and the
periodic multipole integrals/matrices as on the CPU and uploads them once per geometry; the device-resident SCF
(`solve`/`density`/`atomicCharges`/`multipoleMoments`) is linear in the uploaded matrices and therefore periodic by
construction (Gamma point = dense matrices). Everything that has no lattice images on the device stays on the host
and is gated with `!m_pbc` in `XTB::Calculation`: the device integral build (`beginComputed`), the device CN, the
device potential (`beginPotential*`, the fully resident loops for GFN2 and GFN1), the device multipole potential
(`beginMultipolePotential`), the in-SCF device D4 dE/dq (`dispersionDedq`), the device D3 and post-SCF D4, and the
device gradient (`calculateGradientGpu`); the P3/P4 host gradient and virial run on the downloaded density. The CUDA
wrapper alone enables this (`XTB::setPeriodicDeviceScf(true)` in `XtbGpuComputationalMethod`); ROCm/Vulkan keep the
P5 swap to the CPU (`EnergyCalculator::setMolecule`, `gpu != cuda`) and `XTB::Calculation` refuses a cell with any
other backend (hard error). `CURCUMA_GPU_PLUGIN_ABI_VERSION` 7 -> 8 (`XTB` gained data members in P4 and here; the
class is shared with the wrappers compiled into the plugins).

What this buys: the SCF eigensolve (the nao^3 part, 1.1 of 2.2 s at 216 H2O on 8 CPU threads, the dominant term for
432 H2O) runs on the device; the host keeps the O(nat^2) image sums, the Ewald/multipole matrix builds and the
periodic D4 (threaded, stage T). Not ported (G1 device kernels, G2-G6): image-summed device integrals, device Ewald
gamma, device multipole Ewald matrices, device periodic D3/D4, device periodic gradient and virial; the host-driven
loop transfers v_ao and the charges/moments per iteration (the fully resident loop is off with a cell).

Gate (`ctest -R pbc_xtb_gpu_parity`, `test_cases/unit/test_pbc_xtb_gpu_parity.cpp`; water dimer 8 and 12 A, gfn1 and
gfn2, `scf_threshold` 1e-9, `-gpu cuda` through `EnergyCalculator` against the CPU; the test SKIPs when the GPU
request fell back; RTX 5080, CUDA plugin of `release_cuda/` in the worktree, 2026-11-15): |dE| 5e-15 / 0 / 5e-15 /
5e-15 Eh, max |dG| 3.9e-11 / 4.2e-11 / 1.3e-10 / 9.0e-11 Eh/A, max |dXi| 5.2e-11 / 5.8e-11 / 2.3e-10 / 1.2e-10 Eh for
gfn1 8 A / gfn1 12 A / gfn2 8 A / gfn2 12 A; the second CUDA run is bitwise identical in energy, gradient and
virial in all four cases (n = 1 repeat each). Verbosity 2 shows `SCF: device-resident GFN2 path (H0/S/L + multipole
integrals uploaded)` with the cell. `pbc_xtb_tblite` keeps its `-gpu cuda` check (1e-8 instead of 1e-10: device
eigensolver rounding) and passes in the CUDA build. Without a cell the CUDA path is untouched (every new gate is
`!m_pbc`); the CPU binary without a cell is 24/24 bitwise to 957afab7 (6 structures x gfn1/gfn2 x 1/4 threads).

Water lattice 216 H2O (648 atoms, 18.6 A, `-sp` + `-dump_gradient`, wall incl. start-up, same binary, one run each):

| method | threads | SCF threshold | CPU (s) | `-gpu cuda` (s) | E (8 printed decimals) | max dG cuda - cpu (Eh/A) |
|---|---|---|---|---|---|---|
| gfn2 | 8 | default | 2.35 | 2.14 | identical | 9.2e-7 |
| gfn1 | 8 | default | 1.85 | 1.52 | identical | 1.1e-6 |
| gfn2 | 1 | 1e-9 | 10.53 | 3.13 | identical | 9.1e-11 |

  The 1e-6 Eh/A at 8 threads with the default threshold is the SCF convergence level plus the per-thread partials of
  the host H0/Pulay gradient (stage T: 1e-7 between thread counts); with a tight SCF at one thread the two paths agree
  to 9e-11. At 8 host threads the device eigensolve saves little (the host still builds every periodic matrix and the
  threaded MKL eigensolve is already fast at nao = 2592); at one host thread the GPU run is 3.4x faster. For 432 H2O
  (nao 5184, 8x the eigensolve work) the device eigensolve is the lever; the host matrix builds (O(nat^2) image sums,
  multipole k loop) then dominate the GPU run and are the next items (G1 device integrals, G3).

SimpleMD: `setupBarostat()` refused every `-gpu` request (N2 rule for the GFN-FF hybrid); it now refuses only when the
method has no native xTB engine (`xtbInstance()`), so `-md ... -gpu cuda -barostat berendsen -compute_pressure true`
runs with gfn1/gfn2 (virial and cell update on the host). Smoke (27 H2O lattice 9.31 A, gfn1, 20 steps of 0.5 fs,
Berendsen thermostat and barostat, `-threads 2`): CUDA and CPU both run, P/V/rho/L columns present, cell moves
(9.3105 -> 9.3102 A on CUDA, -> 9.3119 A on the CPU; the two runs start from different unseeded velocities, so their
pressures differ from step 0 on; n = 1, no physics claim).

Not tested: MD with a cell on the GPU beyond the 20-step smoke (no energy conservation or seeded CPU/CUDA
trajectory comparison), `-gpu cuda` with a cell and threads > 1 beyond the 8-thread lattice single points, open-shell
or charged cells on the device, d-shell elements (refused with a cell anyway), the ROCm/Vulkan refusal on a machine
with those plugins.

### G2-G6 profile: host time per step in CUDA mode (2026-11-15, Claude Generated)

Water boxes of `build_waterbox.py` (operator project GFN-FF-Wasser; 216 H2O built with `--n 216 --seed 1`, L 18.64 A;
432 H2O = `box432_s1.xyz`, L 23.49 A, nao 5184), gfn2, `-sp` + `-dump_gradient`, `-gpu cuda`, `-threads 8`, RTX 5080,
binary of G1 plus verbosity-3 timers (no numerical change), one run each on a shared machine (repeat scatter of single
items 10-20 %):

| item (ms) | 216 H2O | 432 H2O | stage |
|---|---|---|---|
| CN / overlap + H0 / orthonormalizer | 17 / 29 / 22 | 56 / 78 / 143 | G1 integrals |
| gamma: Ewald matrix / Klopman-Ohno remainder | 118 / 15 | 396 / 59 | G2 |
| multipoles: AO dp/qp image sums / interaction matrices | 81 / 104 | 334 / 728 | G1 / G3 |
| pre-SCF upload (H0/S/L, multipole ints) | 55 | 143 | - |
| SCF 11 it: total / device eigensolve / host D4 dE/dq | 540 / 284 / 220 | 1758 / 1062 / 532 | G4 (dE/dq) |
| post-SCF D4: pairs / ATM / CN fold | 21 / 768 / 41 | 42 / 1415 / 102 | G4 |
| gradient: W / repulsion+XB / H0-Pulay / Ewald / KO / multipole / CN chain | 58 / 23 / 63 / 9 / 16 / 124 / 37 | 331 / 70 / 136 / 34 / 59 / 894 / 133 | G5 |
| TOTAL (`Timing breakdown`) / wall | 2284 / 2710 | 7284 / 7610 | |

Order of benefit at 432 H2O: ATM (1.4 s), multipole gradient (0.9 s), multipole matrices (0.7 s), in-SCF D4 dE/dq
(0.5 s), Ewald gamma (0.4 s), AO multipole image sums and W (0.33 s each). The virial is not part of this table (`-sp`
requests none); its extra cost is in the MD measurements below. `CURCUMA_PBC_GPU_HOST=<tokens>` (diagnostic env,
`atm`, `mpmat`, `mpgrad`, `d4dedq`, `all`) keeps a ported stage on the host in the same binary.

### G4a periodic D4 ATM on the device (2026-11-15, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

- `GpuScfBackend::dispersionATMPeriodic` (CUDA: `XtbGpuContext::dispersionATMPeriodic`, kernels `k_atm_pbc_nb<FILL>`,
  `k_d4_atm_pbc_w`): image neighbour list of every atom on the device (candidates atom x translation in ascending
  order, warp-ballot compaction, two passes), then a warp-per-atom gather over the unordered pairs of list points;
  every translation class of triples appears once per member, so E = 1/3 sum e_atom, the gradient and dE/dCN of an
  atom take the derivative with respect to the centre only, and the virial is 1/3 of the per-atom sum (the
  `triple_scale` weights of `computeATMPeriodic` are the multiplicities of this gather). No atomics; the host adds
  the per-atom results in atom order. The q = 0 C6 matrices come from `buildAtmC6Flat` (host, 10 ms at 432 H2O).
  Per list entry the kernel stores r, 1/r^2, r0_aX, sqrt(C6_aX) and dC6_aX/dCN_a / C6_aX, so the triple loop has one
  sqrt, two divisions and one cbrt (FP64 sqrt/div are slow on consumer GPUs; 490 -> 240 ms at 432 H2O).
- G1 repair found here: with a cell `beginBasis` was skipped, so the device lacked the CSR offsets of the fixed-order
  charge/moment sums and used the atomic kernels; the CUDA run of 216 H2O was not bitwise reproducible (dE 2e-13 Eh
  between two runs, the dimer cells of the G1 gate were). `XTB::Calculation` now uploads the (geometry-free) basis with
  a cell as well; 216 H2O gfn2 and gfn1, 1 and 8 threads: second run md5-identical.
- Plugin ABI 8 -> 9 (new `GpuScfBackend` virtual). Diagnostic `CURCUMA_PBC_GPU_HOST=atm` keeps the host loop.
- Gates: `pbc_xtb_gpu_parity` (now also the water trimer in a 6.5 A cube; virial limit 1e-7 as specified) gfn1/gfn2:
  |dE| <= 5.3e-15 Eh, max |dG| <= 1.3e-10 Eh/A, max |dXi| <= 2.3e-10 Eh, second CUDA run bitwise in 6/6; 216 H2O gfn2
  (`-scf_threshold 1e-9`, 8 threads) device ATM vs `CURCUMA_PBC_GPU_HOST=atm` in the same binary: dE 2.3e-13 Eh,
  max dG 4e-15 Eh/A; `pbc_xtb_*` 5/5 in `release/` and `release_cuda/`; CPU without a cell 24/24 md5 to 957afab7,
  CUDA without a cell 24/24 md5 to a 957afab7 build with the `release_cuda/` configuration (AVX2/AVX512/MKL/march native).

### G3 periodic multipole interaction matrices on the device (2026-11-15, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

- `GpuScfBackend::multipoleMatricesPeriodic` (CUDA kernel `k_mp_amat_pbc`): the pair body of
  `buildMultipoleInteractionMatricesPeriodic` statement for statement (host k list with weights and P2b factors, direct
  box translations), one thread per pair jat <= iat writing both mirrored elements (each element once, deterministic);
  the 18 nat x nat blocks are downloaded into the host matrices (the potential and energy stay on the host), the
  diagonal self terms are added by the host as before. `CURCUMA_PBC_GPU_HOST=mpmat` keeps the host loop. Plugin ABI 10.
- Time (`interaction mats`, 8 threads, one run each, same binary with and without the env switch): 216 H2O 95 -> 22 ms,
  432 H2O 719 -> 108 ms (the 18 x 13 MB download is most of the device time).
- Gates: `pbc_xtb_gpu_parity` 6/6 unchanged within 1e-15 (gfn2 trimer |dE| 4.4e-15 Eh), second CUDA run bitwise in 6/6;
  device vs host matrices in the same binary: 216 H2O dE 2.3e-13 Eh (1 ulp of E), max dG 3.8e-15 Eh/A; 432 H2O dE 0,
  max dG 4.9e-15 Eh/A; `pbc_xtb_*` 5/5 in both builds; without a cell CPU 24/24 and CUDA 24/24 md5 to 957afab7.

### G5a periodic multipole gradient on the device (2026-11-15, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

- `GpuScfBackend::multipoleGradientPeriodic` (CUDA kernel `k_mp_grad_pbc<WANT_V>`): the pair body (`pairTerms`) of
  `addMultipoleGradientPeriodic` statement for statement, one thread per pair jat <= iat; gradient and dE/d(mrad)
  added with `gfix_add` (64-bit fixed point, 2^-52, order independent, the GFN1/GFN2 gradient scheme of
  [GPU_REPRODUCIBILITY.md](GPU_REPRODUCIBILITY.md)), the virial (-dE/d(eps)) as a fixed-order block tree per 256
  threads then one fixed-point add per block; a non-finite or out-of-range contribution (`g_gfix_bad`) fails the call
  and the host blocks run. The dE/d(mrad) -> dE/dCN fold stays on the host. `CURCUMA_PBC_GPU_HOST=mpgrad`. ABI 11.
- Time (`multipole` in `periodic gradient`, 8 threads, without virial): 216 H2O 135 -> 18 ms, 432 H2O 901 -> 110 ms.
- Gates: `pbc_xtb_gpu_parity` 6/6 unchanged (incl. the virial, max |dXi| 2.3e-10 Eh, set by the SCF), bitwise run to
  run; device vs host in the same binary: 216 H2O dE 0, max dG 1.1e-14 Eh/A, 432 H2O dE 0, max dG 1.8e-14 Eh/A; 216 H2O
  second device run md5-identical; MD 216 H2O gfn2 20 steps of 0.5 fs, `-seed 7 -compute_pressure true`, all stages
  on the device vs `CURCUMA_PBC_GPU_HOST=all`: after 20 steps E_pot identical to 6 decimals, T 458.811166 / 458.811149 K,
  P 15331.03 / 15331.03 bar (lattice start, not equilibrated); `pbc_xtb_*` 5/5 in both builds; without a cell CPU 24/24
  and CUDA 24/24 md5 to 957afab7.

### G4b periodic in-SCF D4 dE/dq on the device (2026-11-15, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

- The BJ factor of a pair does not depend on the charges, so `GpuScfBackend::beginDispersionImages` (after the existing
  `beginDispersion` upload, once per geometry) builds D(i, j) = sum_T [s6/(r^6 + R0^6) + s8 r4r2_ij/(r^8 + R0^8)] over the
  images of the 50 Bohr pair cutoff (kernel `k_d4_pbc_dmat`, one thread per element); per SCF iteration
  `dispersionDedqPeriodic` (kernel `k_d4_dedq_pbc`, one thread per atom, fixed j order) contracts the host-built
  W/dWq with D, the terms of `computeEnergyAndGradientPeriodic`'s dE/dq. `XTB::m_pbc_d4_dev_ready` follows
  `m_d4_prepared`; `CURCUMA_PBC_GPU_HOST=d4dedq`. Plugin ABI 12. The post-SCF D4 pair term (42 ms at 432 H2O) stays
  on the host.
- Time (`potential build - of which D4`, 8 threads): 216 H2O 13.6 -> 2.5 ms per iteration, 432 H2O 45.9 -> 10.2 ms per
  iteration (the rest is the host `buildRefWFlat`).
- Gates: `pbc_xtb_gpu_parity` 6/6 (gfn2 |dE| <= 8.9e-15 Eh, gradient and virial as before), bitwise run to run; device
  vs host in the same binary: 216 H2O dE 0, max dG 4.6e-15 Eh/A, 432 H2O dE 0, max dG 4.2e-15 Eh/A; second 216 H2O
  device run md5-identical; `pbc_xtb_*` 5/5 in both builds; without a cell CPU 24/24 and CUDA 24/24 md5 to 957afab7.

### G2-G6 status and timing after G4a/G3/G5a/G4b (2026-11-15, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

Same binary (`release_cuda/`, ABI 12), RTX 5080, `-threads 8`, gfn2; "before" = `CURCUMA_PBC_GPU_HOST=all` (every
ported stage on the host, i.e. the G1 split), "after" = default. `-sp` wall with `-dump_gradient` (no virial); MD:
`-md -time_step 0.5 -max_time 10` (20 steps) `-seed 7 -compute_pressure true` (gradient + virial every step), mean of
the per-step wall column over the 20 steps; one run each, machine shared (load 1-7):

| box | `-sp` before (s) | `-sp` after (s) | MD step before (s) | MD step after (s) |
|---|---|---|---|---|
| 216 H2O (648 atoms) | 2.18 | 1.38 | 1.75 | 0.99 |
| 432 H2O (1296 atoms, nao 5184) | 6.39 | 4.09 | 6.90 | 3.73 |

One 432 H2O MD step with the virial after G4b (verbosity 3, ms): TOTAL 3753 = setup 918 (gamma: Ewald matrix 348 +
KO remainder 55, AO multipole image sums 175, multipole matrices 104 on the device, orthonormalizer 115, S/H0 68) +
pre-SCF 120 + SCF 8 it 1007 (device eigensolve 786, D4 dE/dq 107) + post-SCF 656 (ATM 309 on the device, CN fold 94,
D4 pairs 48) + gradient 1048 (W 287, multipole + virial 293 on the device, H0/Pulay 143, CN chain 126, Ewald 67, KO 56,
repulsion 68). Not ported (largest first): Ewald gamma matrix (G2, 0.35-0.40 s; reciprocal part is a rank update
C C^T + S S^T over 6602 k vectors, a cuBLAS dsyrk candidate), W = C f eps C^T (0.29 s, could use the device C), AO
multipole image sums and S/H0 image sums (G1 device integral kernels, 0.24 s), H0/Pulay and CN-chain image loops
(G5 rest, 0.27 s), KO remainder (0.11 s), periodic D3 (GFN1), Ewald gradient/virial (0.07 s); G6 (virial) is ported
only where its pair loop is (ATM, multipole), the other virial terms stay on the host.

Regression after G4b (2026-11-15, test selection of the regression run below, `-j1`): `release/` 367 of 368 passed,
`release_cuda/` 638 of 639 passed; the only failure in both is the pre-existing `xtb_cpscf`; 3 tests `Not Run
(Disabled)` as before.

## Regression run (2026-11-14, `release/`, `ctest -R 'gfn|xtb|sqm|gfnff|pbc|ewald|simplemd|restart|md_|gpu|plugin|molecule|xyz|parser|curcumaopt|structinv' -j1`)

369 tests selected, 2 failed before the fix below: `pbc_md_nonperiodic_cell` (its check (a) assumed gfn2 + cell to be
non-periodic; now uses `eht` there and checks gfn2 + cell as periodic with DOF 3N-3 in a new (d), green) and `xtb_cpscf`
(checks (C), 8e-5 .. 3e-4 against 1e-6 .. 1.5e-4; fails identically with the executable of the main `release/` build of
Nov 9, i.e. pre-existing and unrelated to this work; not changed). 3 tests `Not Run (Disabled)` as before. The
WP-PERIODIC-NONBONDED A6b statement "native gfn1/gfn2 have no periodic energy" is superseded by P5.

## Klopman-Ohno switch radius: model constant, not a convergence parameter (2026-11-16, Claude Generated; 🤖 AI-generated, ⚙️ machine-tested)

Question from the GFN2 water box (density 1.31 g/cm^3 at rc = 10 Bohr, pressure shifts of 10^4 bar with rc): is the
rc dependence of `-xtb.pbc_coulomb_switch_bohr` a bug or the model? Answer: the model (tblite's definition); no code
change. Measurements and scripts: project folder GFN-FF-Wasser, subfolder 06_gfn2/switch (operator's scratch, not in the repo).

- **Code = tblite statement by statement**: kernel `f KO + (1 - f)/r - erf(a r)/r` per image (`effective.f90:404-412`),
  quintic `fsmooth` over [rc - 1, rc] (`:449-466`), `rcut = 10` hard-coded (`:114`); `buildGammaMatrixPeriodic` and
  the gradient/virial loop (`xtb_periodic.cpp`, `koSwitch`) are the same. The CUDA path uses the host gamma.
- **Why it cannot converge**: `KO - 1/r = -1/(2 gbar^2 r^3) + O(r^-5)` with `gbar` the shell-pair hardness average.
  For two neutral molecules the r^-3 term does not cancel because `gbar` differs per shell pair: point model
  `u(R) = -C/R^3`, `C = sum_{a in A, b in B} q_a q_b / (2 gbar_ab^2)`. GFN1 (harmonic mean) gives `C = D_A D_B / 4`,
  `D = sum_a q_a / g_a`, i.e. orientation-independent and attractive between equal molecules; GFN2 (arithmetic mean)
  is not separable but behaves alike. GFN2 water shell charges (64-H2O box image, rc 10): O s +0.24, O p -0.90,
  H +0.32 e, giving `C` = 0.135 Eh Bohr^3. A lattice sum of R^-3 in 3D diverges logarithmically:
  `dE/dln rc = -2 pi rho C N`, `dP/dln rc = -2 pi rho^2 C` (-10200 bar per e-fold at 1.29 g/cm^3).
- **Measured** (1-step MD, CPU, `-scf_threshold 1e-9`, n = 1 image each; 64 H2O GFN2 box end image, rho 1.29; 432 H2O
  GFN-FF box image, rho 1.02):

| rc (Bohr) | 6 | 8 | 10 | 12 | 15 | 20 | 25 | 30 | 40 | 50 |
|---|---|---|---|---|---|---|---|---|---|---|
| 64: E - E(10) (Eh) | -0.066 | +0.212 | 0 | +0.042 | -0.049 | -0.150 | -0.232 | -0.295 | -0.404 | -0.490 |
| 64: P_vir - P_vir(10) (bar) | -8690 | +9216 | 0 | -4988 | -1720 | -11909 | -10535 | -13623 | -20552 | -23269 |
| 432: E - E(10) (Eh) | +0.075 | +0.616 | 0 | -0.246 | -0.646 | -1.124 | -1.514 | -1.843 | -2.362 | |
| 432: P_vir - P_vir(10) (bar) | -14940 | +4679 | 0 | +948 | -1448 | -4031 | -6354 | -7940 | -10753 | |

  For rc >= 8 E(rc) - E(10) is reproduced by the remainder sum at fixed rc = 10 shell charges to within 6.2 % on both
  images (e.g. 64: +0.223 / +0.044 / -0.469 Eh at 8 / 12 / 50); the rest is the SCF response (larger at rc = 6:
  16 % and 49 %, where the switch reaches the first neighbour shell). Beyond 15 Bohr E falls
  linearly in ln rc (64: -0.37 Eh per e-fold, point model -0.35; 432: -1.79, point model -1.69). P falls by
  -10^4 bar per e-fold (64) and -6000 to -10000 bar (432); below 15 Bohr it oscillates by +-10^4 bar with the
  pair structure in the 1-Bohr switch window (the fixed-charge estimate of dP has the measured sign at every rc of both images; magnitude mostly within a factor 2, at 432/15 Bohr a factor 10).
- **Consistency**: `fsmooth` has zero first and second derivatives at both ends, so E is C2 in the geometry and in rc.
  Cell-scaling FD (central, isotropic, h = 1e-4 / 2e-4, scf_threshold 1e-10, 64 image) minus analytic virial:
  +31.2 / +20.7 bar at rc 8, 10 and 12 alike (spread 0.07 bar), so the KO term is virial-consistent. The offset is the
  hard D4 ATM cutoff (25 Bohr, energy discontinuous, FD h-dependent): with `-xtb.d4_atm_cutoff 40` -1.2 / -5.1 bar,
  60: -7.9 bar (h = 1e-4); remaining scatter from the other hard cutoffs (D4 pairs 50, CN). Not a code error.
- **Consequence**: a periodic GFN2 (and GFN1) liquid density is defined only together with rc; there is no rc -> inf
  limit, and the vacuum GFN2 of a cluster corresponds to an rc of the cluster size (expected denser; not measured).
  A larger default would not be "more converged". Proposal, not implemented: keep 10 (tblite parity) and report rc with
  every periodic GFN density; a larger value only as an explicitly named model variant. NPT 64-H2O runs at rc = 8 and
  12 are queued in the project folder (density response, pending).
- Not tested: GFN1 numbers (only the analytic form), other liquids, cells below 2 rc (image sums are exact there too).

## Not implemented / not tested (updated per stage)

- k-point sampling, non-orthorhombic cells, slab/1D periodicity, the vacuum (dipole) boundary term, off-diagonal
  strains (a sheared cell is refused), threading of the periodic CN/D3, Klopman-Ohno and Ewald matrix builds, the
  structure-factor factorisation of the multipole k sums, d-shell elements with a cell (refused), implicit solvation
  with a cell (refused), `-large_system_mode` with a cell (its driver ignores the cell: documented as non-periodic, not
  refused), GPU with a cell: CUDA in the G1 split plus the device stages G3/G4a/G4b/G5a (the rest on the host, see G2-G6
  status), ROCm/Vulkan CPU fallback (not planned).
- Not tested: open shell, charged cells (background term differs from tblite), halogen bond with a cell (no Br/I
  structure), thread-count independence of the periodic energy and gradient, MD energy conservation with a cell,
  the GPU fallback on a machine with a GPU plugin, anything beyond H/O, real liquid boxes (only a water dimer and
  water lattices were run). Human production testing pending for everything in this document.
