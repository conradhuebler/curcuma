# Analytic RHF Hessian (`-hessian -method hf -hessian.analytic true`, Nov 2026)

🤖 AI-generated, ⚙️ machine-tested only. Human production testing is pending.

Second nuclear derivatives of the closed-shell native HF energy, with the orbital response from the
coupled-perturbed HF (CPHF) equations. Code: `src/core/energy_calculators/qm_methods/qm_hessian.{h,cpp}`
(derivation in the file header), entry `QMEngine::computeHessian`, capability hook
`Hessian::CalculateHessianAnalytic` in `src/capabilities/hessian.cpp`.

## Use

```sh
curcuma -hessian water.xyz -method hf -qm.basis def2-SVP -hessian.analytic true
```

| parameter | default | meaning |
|---|---|---|
| `-hessian.analytic` | `false` | analytic Hessian for `-method hf`; any other method, or a refused case, runs the numerical Hessian after a warning |
| `-hessian.cphf_threshold` | `1e-10` | largest \|residual\| (Eh) of each of the 3N CPHF columns |
| `-hessian.cphf_max_iterations` | `100` | conjugate-gradient iterations per column; not converged = fallback |

The numerical (finite-difference) Hessian stays the default and the reference. Refused, with a warning and the
numerical path: open shell, Kohn-Sham functionals, `hf-3c` and every other method, ECP atoms, `-qm.ri j|jk`,
g or higher shells, an unconverged SCF, a non-positive CPHF operator, a non-converged CPHF.

## Formula

    H_xy = sum P h^xy + 1/2 sum D (mn|ls)^xy - sum W S^xy + E_nn^xy + Tr(P^y F^x) - Tr(W^y S^x)
    F^x  = h^x + G^x[P],   P^y = -1/2 P S^y P + 2 (C_v U^y C_o^T + h.c.),   W^y = d/dy (1/2 P F P)
    (eps_a - eps_i) U_ai + [C_v^T G[D(U)] C_o]_ai = -(F^y_ai - eps_i S^y_ai + [C_v^T G[-1/2 P S^y P] C_o]_ai)

Pople, Krishnan, Schlegel, Binkley, Int. J. Quantum Chem. Symp. 13, 225 (1979); Yamaguchi, Osamura, Goddard,
Schaefer, *A New Dimension to Quantum Chemistry* (1994); Helgaker, Jorgensen, Olsen, ch. 9 (integrals).
Derivative integrals: McMurchie-Davidson with shifted Gaussian powers, centre derivatives once or twice per centre,
the last centre (nucleus of V, fourth shell of an ERI) from translational invariance; cartesian basis, densities and
matrices mapped through the spherical transform as in the WP8 gradient. The 2e part runs over canonical shell
quartets, integral-direct (no derivative tensor is stored), without screening. CPHF: preconditioned conjugate
gradient on the virtual x occupied space, all 3N columns in lockstep; with the stored ERI tensor the J/K of all
response densities of one iteration are two GEMMs over the tensor, with `-qm.scf_direct on` (or `-gpu`) the engine's
`buildJK` runs once per density. Units: Eh/Bohr^2 from the engine, Eh/Angstrom^2 inside the capability.

## Validation (Nov 5, 2026; release build, 4 threads, this container: 4 vCPU, no GPU)

`ctest -L qm_hess` (10 tests, def2-SVP; 10/10 pass, 218 s), reference `test_cases/qm_hessian/reference.json` from
`scripts/hf_hessian_reference.py` (PySCF 2.14.0 `pyscf.hessian.rhf`, conv_tol 1e-13, conv_tol_cpscf 1e-12, plus the
4-point finite difference of PySCF's analytic gradient; regenerated with PySCF 2.14.0 on Nov 5: same Hessians to
<= 3.5e-12). The other rows are manual runs of `dump_qm_hessian --fd --stencil 4 --h 5e-4` against a PySCF analytic
Hessian computed at the same geometry. All in Eh/Bohr^2 (max |H| 0.5..1.1).

| molecule | basis | vs FD of own gradient | vs FD of PySCF gradient | vs PySCF analytic | asymmetry | translation / rotation |
|---|---|---|---|---|---|---|
| H2O | def2-SVP | 3.7e-10 | 4.7e-9 | 3.3e-8 | 2.1e-11 | 2.1e-11 / 3.1e-11 |
| NH3 | def2-SVP | 2.0e-9 | 6.8e-9 | 1.7e-8 | 1.6e-11 | 1.3e-11 / 3.2e-11 |
| CH4 | def2-SVP | 6.6e-11 | 4.4e-9 | 1.1e-8 | 7.6e-12 | 5.5e-12 / 9.2e-12 |
| HF | def2-SVP | 6.1e-12 | 9.1e-9 | 4.1e-8 | 4.4e-16 | 4.1e-15 / 2.0e-13 |
| H2O distorted | def2-SVP | 1.2e-10 | 4.1e-9 | 8.4e-8 | 1.8e-11 | 1.2e-11 / 3.3e-11 |
| formaldehyde distorted, non-planar | def2-SVP | 1.1e-9 | 7.6e-9 | 8.9e-8 | 2.8e-11 | 1.3e-11 / 4.4e-11 |
| formaldehyde distorted, `-qm.scf_direct on` | def2-SVP | 9.6e-9 | 7.6e-9 | 8.9e-8 | 3.2e-11 | 1.2e-11 / 4.5e-11 |
| formaldehyde (HF-3c minimum) | def2-SVP | 9.6e-9 (3.2e-9 at h = 1e-3) | - | 7.8e-8 | 2.8e-11 | 2.4e-11 / 4.3e-11 |
| H2O distorted, `-qm.cartesian_d true` | def2-SVP | 1.3e-10 | - | 9.0e-8 (PySCF cart) | 2.1e-11 | 1.2e-11 / 3.5e-11 |
| HCl (library `hydrogen-chloride`) | def2-SVP | 1.1e-8 | - | 7.2e-9 | 2.0e-12 | 1.1e-14 / 2.3e-11 |
| H2S (library `hydrogen-sulfide.d-shell-validation`) | def2-SVP | 1.1e-9 | - | 8.3e-9 | 1.7e-11 | 9.2e-12 / 3.4e-11 |
| H2O | def2-TZVP | 4.4e-10 | - | 1.8e-8 | 1.3e-11 | 6.2e-12 / 1.5e-11 |
| NH3 | def2-TZVP | 9.1e-10 | - | 1.2e-8 | 4.4e-11 | 4.2e-11 / 6.1e-11 |
| CH4 | def2-TZVP | 4.1e-10 | - | 2.5e-8 | 2.1e-11 | 3.0e-11 / 3.9e-11 |
| HF | def2-TZVP | 9.3e-11 | - | 5.0e-9 | 2.7e-15 | 4.6e-15 / 1.0e-12 |
| H2O distorted | def2-TZVP | 4.0e-11 | - | 1.8e-8 | 1.1e-11 | 4.7e-12 / 8.1e-12 |
| formaldehyde distorted | def2-TZVP | 2.8e-10 | - | 3.2e-8 | 5.8e-11 | 4.9e-11 / 8.6e-11 |
| formaldehyde (HF-3c minimum) | def2-TZVP | 8.5e-9 | - | 1.7e-8 | 4.4e-11 | 3.5e-11 / 7.0e-11 |

- FD: 4-point central difference of the fully assembled `QMMethod::getGradient()` (Eh/Angstrom by contract,
  converted), h = 5e-4 Bohr, SCF orbital gradient 1e-10. The larger FD rows (~1e-8) are gradient noise of the
  displaced SCFs, not truncation: formaldehyde/SVP drops from 9.6e-9 to 3.2e-9 when h is doubled.
- PySCF's analytic Hessian is itself 7e-9..9e-8 off the FD of its own gradient on the def2-SVP cases (not changed by
  conv_tol_cpscf 1e-12 or direct_scf_tol 1e-16, cause not investigated); curcuma agrees with PySCF's FD to <= 9.1e-9.
- Rotational invariance is checked as H (L r) = L g (L the rotation generator), which holds at every geometry.
- Fixed-density pieces (`qm_hess_pieces_H2O_distorted`, def2-SVP): S^x, (T+V)^x, derivative Fock G^x, the 1e and 2e
  skeleton second derivatives vs finite differences of `buildOverlap`, `buildCoreHamiltonian`,
  `buildCoulomb - buildExchange/2`, `gradientOneElectron`, `gradientTwoElectron` at fixed P: relative <= 7.4e-9 (FD
  truncation); the quartet digest reproduces the stored-tensor J - K/2 to 8e-15.
- Newton minimisation with the analytic gradient and Hessian (`dump_qm_hessian --newton`) converges quadratically,
  e.g. H2O/def2-SVP |g| 8.8e-2, 8.0e-3, 2.9e-4, 7.2e-8, 8.5e-14 Eh/Bohr; at the minima the six smallest |eigenvalues|
  are <= 3.2e-11 (non-linear; five for linear HF).

**Frequencies vs PySCF** at those Newton minima (PySCF gradient there <= 1.7e-10 Eh/Bohr; PySCF analytic Hessian at
the same geometry, both with PySCF's isotope masses, translations and rotations projected): max |diff| per molecule
(cm^-1)

| basis | H2O | NH3 | CH4 | HF | H2CO |
|---|---|---|---|---|---|
| def2-SVP | 1.2e-4 | 3.0e-5 | 1.8e-5 | 1.3e-4 | 3.9e-5 |
| def2-TZVP | 4.9e-5 | 1.1e-4 | 9.2e-5 | 1.4e-5 | 1.6e-5 |

`qm_hess_cli_frequencies_H2O` checks analytic vs numerical (step 5e-4, tight SCF) through the `-hessian` CLI
(1477.00/4624.00/4787.40 cm^-1 both, at the non-stationary `qm_1e/H2O.xyz`). No ORCA comparison was made.

## Timing (`curcuma -hessian -method hf -threads 4 -no_bmt`, wall time incl. start-up and SCF, Nov 5 2026, 4 vCPU)

Numerical = the default `-hessian` path (threaded central differences of the gradient, step 5e-3, default SCF
thresholds); analytic = `-hessian.analytic true`. Frequencies of the two differ by <= 0.4 cm^-1 (up to the first 12 modes compared per molecule).

| molecule | basis | n_bf | numerical (s) | analytic (s) |
|---|---|---|---|---|
| H2O (`qm_1e`) | def2-SVP | 24 | 0.6 | 0.8 |
| H2O | def2-TZVP | 43 | 3.9 | 2.7 |
| formaldehyde (HF-3c minimum) | def2-SVP | 38 | 3.8 | 3.4 |
| formaldehyde | def2-TZVP | 74 | 56.8 | 10.3 |
| methanol (library) | def2-SVP | 48 | 7.3 | 1.8 |
| dimethyl ether (library) | def2-SVP | 72 | 52.5 | 31.4 |
| dimethyl ether | def2-TZVP | 129 | killed (SIGKILL after 53 s, 15 GB machine) | 173.9 |
| dimethyl ether, numerical with `-qm.scf_direct on` | def2-TZVP | 129 | 676 | 155 (second run, `-verbosity 2`) |

Dimethyl ether/TZVP analytic at `-verbosity 2`: 1e 0.7 s, 2e derivative integrals 99.6 s, CPHF 40.5 s (431 J/K
builds, <= 16 CG iterations), total 154 s in this run: the unscreened 2e derivative pass dominates. The numerical
default was killed there; likely out of memory (four workers, each n^4 * 8 B = 2.2 GB of ERI tensor plus
scratch; not measured), a property of the numerical path.

## Tested / not tested / not implemented

- Tested: H, C, N, O, F, S, Cl; closed shell, neutral; def2-SVP (d) and def2-TZVP (d, f); spherical and cartesian d
  (SVP only); stored tensor and integral-direct SCF; 1 to 4 atoms for accuracy, methanol and dimethyl ether
  (6 and 9 atoms) for timing and CLI frequencies only; 4 threads.
- Not tested: charged molecules, K-Kr (all-electron, not refused), molecules above 9 atoms, near-degenerate
  HOMO-LUMO, transition states, MINIX and other def2 sets, f shells with cartesian d, thread counts other than 4 for
  accuracy; the `-hessian.read` path; `-gpu` (the CPHF J/K then go through `buildJK`, never run here).
- Not implemented: UHF/ROHF, Kohn-Sham (XC kernel and grid response), hf-3c (D3/gCP/SRB second derivatives),
  ECP second derivatives, RI-J/RI-JK Hessian, g shells, GPU kernels, screening of the derivative quartets, a Hessian
  file output (the capability prints frequencies only, as before), dipole derivatives / IR intensities.
- Cost: the 2e derivative pass is O(N^4) over all multi-centre quartets without screening; the batched CPHF holds
  3N n x n response matrices next to the stored n^4 ERI tensor.
- Human production testing is pending.
