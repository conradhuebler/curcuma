# WP: Periodic native GFN1/GFN2 (orthorhombic cell, Gamma point, CPU)

**Status (interim, 2026-11-14): P0 plan, P1 periodic classical terms and P2a periodic GFN1 energy done; P2b (GFN2
multipole Ewald), P3 (gradient), P4 (virial), P5 (integration) open.** `-method gfn1` with an orthorhombic unit cell
gives the periodic Gamma-point energy (tblite agreement 3-7e-11 Eh on the water dimer, see P2a); its gradient is not
periodic yet, so `-opt`/`-md`/`-hessian` with a cell are refused until P3; `-method gfn2` with a cell is refused until
P2b. Without a cell nothing changes (bitwise gate below).
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

## Not implemented / not tested (updated per stage)

- k-point sampling, non-orthorhombic cells, slab/1D periodicity, the vacuum (dipole) boundary term.
- Not yet: everything (P0). Open shell, charged cells, d-shell elements, GPU with a cell, implicit solvation with a cell,
  `-large_system_mode` with a cell will be listed here with their status at each stage.
