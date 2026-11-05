# Native QM SCF: level shift and wavefunction stability analysis (Nov 4, 2026)

> 🤖 AI-generated, ⚙️ machine-tested. One machine (4 vCPU container, no GPU), n = 1 per number, not reviewed by the
> operator. Human production testing is pending.

Native ab-initio SCF (`QMEngine`: `hf`, `hf-3c`, `lda`/`pbe`/`b3lyp`/`tpss`/`scan`, closed shell `qm_scf.cpp`, open
shell `qm_uscf.cpp`). Both options are opt-in; with them off the SCF is bitwise the one before (checked below).

## Level shift (`-qm.scf_level_shift b`, `-qm.scf_level_shift_off g`, `-qm.scf_mode level-shift`)

- Diagonalised matrix F + b (S - S P_s S) (closed shell S - 1/2 S P S; Saunders/Hillier 1973), applied after the
  DIIS/ADIIS extrapolation (as PySCF). DIIS error, energy and the three convergence tests use the unshifted F.
- When a shifted iteration meets every test, the unshifted (extrapolated) F of that iteration is diagonalised once more
  and its orbitals, orbital energies and density are returned (the rule of the run without shift). A first version
  continued the loop unshifted instead: undamped Roothaan steps (ozone) left the converged point again.
- `scf_level_shift_off g > 0` removes the shift once ||g||_F < g. `scf_mode level-shift` = diis + b = 0.2 Eh.
- Host loops only. With `-gpu` the device-resident RHF loop is skipped and counted (`qm_scf`); not run (no GPU here).
- Code: `ScfLevelShift` in `qm_scf_internal.h`, used by `runSCF` and `runUSCF`.

Result does not depend on b (dump_qm_energy, scf_threshold 1e-10, def2-SVP; iterations from 2-thread probes, the
differences from `ctest -L qm_stability` at 1 thread):

| case | b (Eh) | iterations | E - E(b = 0) (Eh) | max dg (Eh/Bohr) |
|---|---|---|---|---|
| water HF, gradient | 0.2 / 0.5 / 1.0 | 15 -> 19 / 26 / 39 | 0 (b = 0.5) | 1.3e-12 (b = 0.5) |
| H2 2.0 A B3LYP sg1, gradient | 0.5 | - | 1e-15 | 2.4e-16 |
| water B3LYP sg1, gradient | 0.5 | 12 -> 34 | - | - |
| CH3 UHF, gradient | 0.5 | 18 -> 27 | -3e-14 | 8.0e-14 |
| ozone HF | 0.3 / 1.0 | 24 -> 29 / 51 | -1.2e-13 / -3.6e-13 | - |
| ozone HF, plain undamped (fails without) | 0.3 | 76 | vs PySCF -1.1e-13 | - |
| N2 2.0 A HF, plain undamped (fails without) | 1.0 | 51 | vs PySCF -2.6e-13 | - |

Where it helps and where not (iterations to scf_threshold 1e-10, "fail" = not converged in the limit):

| case | without | b = 0.3 | b = 0.5 | b = 1.0 |
|---|---|---|---|---|
| ozone HF, `scf_mode plain`, `scf_damping 1` (undamped Roothaan), 150 it. | fail | 76 | - | 127 |
| N2 at 2.0 A HF, plain, damping 1, 150 it. | fail | fail | - | 51 |
| water, O-H 1.8 A HF, plain, damping 1, 200 it. | fail | - | fail | fail |
| ozone UHF (`-qm.reference uhf -qm.uhf_guess_mix true`), DIIS, 200 it. | 60 | 99 | 56 / 73 / 70 (off 1e-1 / 1e-2 / 1e-3) | fail |
| ozone RHF, DIIS | 24 | 29 | 24 / 24 / 26 (off 1e-1 / 1e-2 / 1e-3) | 51 |

With the default DIIS no tested case failed without the shift, and the shift always cost iterations unless it was
removed early (`scf_level_shift_off`). It rescues oscillating undamped Roothaan iterations. No DIIS failure case was found
in this test set, so the benefit for DIIS runs is not demonstrated.

## Stability analysis (`-qm.scf_stability off|check|follow`, `-qm.scf_stability_max_follow 3`)

- After a converged SCF (`QMEngine::Calculation`): lowest eigenvalue of the real orbital Hessian (A + B) by a
  Davidson solve (`qm_stability.cpp`, start: three lowest diagonal unit vectors plus one fixed vector with all
  components, correction r / (theta - diag), residual 1e-6, collapse at 40 vectors). The product (A + B) x is one J/K
  build of D = C_v x C_o^T + transpose (UHF: one two-density build); the MO Fock blocks make it exact for
  non-canonical orbitals. Operators: RHF -> RHF (singlet, `rhf_rhf`), RHF -> UHF (triplet, `rhf_uhf`), UHF -> UHF
  (`uhf_uhf`). Eigenvalue < -1e-5 Eh: warning (verbosity >= 1).
- KS response: singlet exact (`qmxc::buildFxc`); triplet and UKS a central difference of the spin V_xc with
  t max|D| = 1e-4 (no spin-resolved f_xc exists in the code). Meta-GGA triplet/UKS: not available (info line).
- `follow`: occupied orbitals rotated by exp(K) with the unit-norm eigenvector (closed form through the SVD), SCF re-run
  from that density; kept only if the energy is lower, else the previous solution is restored. A singlet instability
  is followed before a triplet one. RHF -> UHF switches the geometry to UHF/UKS (alpha +v, beta -v; it stays
  unrestricted for later SCFs of the same geometry, `m_broken_symmetry_serial`); refused with `-qm.reference rhf`, for
  closed-shell-only wrappers (mp2, ri-mp2, lmp2, ccsd, double hybrids) and meta-GGAs.
- Report: energy decomposition `stability` (`rhf_rhf`, `rhf_uhf`, `uhf_uhf`, `stable`, `followed`), printed by
  `dump_qm_energy` as `stability <key> <value>`. Convention: eigenvalues of (A + B) in Eh; the RHF energy Hessian is
  4 (A + B). PySCF's `stability()` prints its own scaling; the comparison below uses a dense matrix of this convention.

Reference: `scripts/qm_stability_reference.py` (PySCF 2.14, dense (A + B) from PySCF's UHF response function with
exact f_xc, `numpy.linalg.eigvalsh`; followed solution from PySCF `stability()` + SCF), `test_cases/qm_stability/reference.json`.
def2-SVP, B3LYP = B3LYP5 on curcuma's sg1 grid. d = curcuma - PySCF (Eh).

| case | E | rhf_rhf / uhf_uhf | d | rhf_uhf | d |
|---|---|---|---|---|---|
| water HF | -75.961008957568 | 0.3523628716 | -5e-11 | 0.2815853005 | -3e-11 |
| H2 2.0 A HF | -0.921734205941 | 0.3889758583 | -3e-14 | -0.2295873842 | 1.6e-13 |
| H2 2.0 A UHF (restricted solution) | same | uhf -0.2295873842 | 4.9e-13 | | |
| water O-H 1.8 A HF | -75.576343840390 | 0.0731864542 | 1.9e-10 | -0.2452823596 | 3.8e-10 |
| ozone HF | -224.069232446976 | 0.0831601761 | -1.1e-11 | -0.2012087798 | 1.6e-10 |
| N2 2.0 A HF (atom guess) | -108.214505731375 | -0.1918336066 | -1.0e-11 | -0.4117922189 | -1.7e-11 |
| CH3 UHF | -39.532955448646 | uhf 0.2817846007 | 9e-11 | | |
| H2 2.0 A B3LYP | -0.992223815209 (PySCF -...261) | 0.4640989352 | -7.3e-9 | -0.0984866833 | 1.1e-8 |
| H2 2.0 A UKS B3LYP | same | uhf -0.0984866991 | -4.8e-9 | | |
| ozone B3LYP | -225.139087849779 | 0.0929297024 | -1.5e-10 | -0.0164562085 | 1.0e-11 |

Followed solutions (`follow`) vs PySCF: H2 2.0 A HF -1.002773459950 Eh (dE 2e-15, <S^2> 0.90437376 both), H2 with
`-qm.reference uhf` (UHF internal) the same, water 1.8 A HF -75.725117786072 (dE -7e-14, <S^2> 1.67223202), H2
B3LYP -1.012098168967 (dE 5.7e-11, <S^2> 0.69845775). FD gate of the assembled surface after following (water 1.8 A,
HF, UHF after the RHF -> UHF switch, 4-point central difference of `-sp` energies with follow, 9 coordinates): max
3.2e-11 Eh/Bohr; level shift CH3 UHF b = 0.5: max 2.0e-10 Eh/Bohr.

Other observations: N2 at 2.0 A from the atom guess follows singlet -> triplet -> UHF internal to -108.546016 Eh
(<S^2> 2.685, three follows, final analysis stable with a zero eigenvalue, the rotation about the bond axis); PySCF
following only the external mode reaches another UHF minimum, -108.562699 Eh. Which solution following reaches depends
on the route; this case is not a ctest. RI-JK (H2 2.0 A): eigenvalues of the RI operator (0.388971 / -0.229505 vs
exact 0.388976 / -0.229587); no external RI reference. Integral-direct J/K gives the stored-tensor values to 4e-16.
Ozone B3LYP at 1 thread: 4.3 min for the two operators (22 products; each V_xc build takes ~3 s at 1 thread on this
grid, 2.5 s for the whole run at 4 threads), so the KS ctests run with 4 threads.

## Verification

- Bitwise with both options off: 13 runs (hf/b3lyp/pbe RI-JK, energy-only and gradient, plain mode, integral-direct,
  UHF/UKS, `-qm.reference uhf` singlet, hf-3c; 1 thread) through `dump_qm_energy` with `CURCUMA_QM_SCF_TRACE=1`:
  energies and gradients at 17 digits plus iteration counts, md5 `24099b0562288ceb312959214c9c6167` for the baseline
  binary (029a72e), stage A and stage B. Densities themselves were not dumped (the gradient at 17 digits stands in).
- `ctest -L qm_stability` (21 tests), `qm_1e|qm_2e|qm_grad|qm|qm_uhf|qm_uks|qm_b3lyp|qm_lda|qm_pbe|qm_dft_grad` after
  stage A: 273/277, the 4 failures were timeouts (`qm_lmp2_exact_threads1_acetic_acid_dimer_fc` 900 s,
  `qm_dft_opt_formaldehyde_distorted_{lda,pbe,b3lyp}` 1800 s) while another job loaded the 4-core machine; those code
  paths run with both options off, where the md5 check shows identical results. After stage B (all QM test binaries
  rebuilt): `ctest -L "qm_1e|qm_2e|qm_grad|qm_uhf|qm_uks|qm_b3lyp|qm_stability"` 151/151 and `method_catalog_*` 2/2.
- The three new library structures add 15 generated `structinv_*`/`structcomp_*` tests (gfn1, gfn2, gfnff): 15/15. In
  this CPU build no structinv test runs a native ab-initio method.

## Not tested / not implemented

- Not tested: the device loops (no GPU in the container: the counted `qm_scf` fallback with level shift was never
  run), ROCm/Vulkan, molecules beyond 3 heavy atoms, other bases than def2-SVP/MINIX, charged systems, RI-JK against
  an external reference, meta-GGA singlet operator against PySCF, `-opt`/`-md` with `follow` (each geometry would be
  analysed and followed again), level shift with RI mixed precision or integral-direct incremental builds.
- Not implemented: real -> complex and RHF/UHF -> GHF instabilities, a spin-resolved analytic f_xc (UKS/triplet use a
  central difference), meta-GGA triplet/UKS operators, stability analysis on the device, a second-order (Newton) SCF
  (stage C, assessment only, see below).
- Human production testing is pending.

## Stage C: second-order SCF (assessment, not implemented)

The stability code already provides what a second-order SCF needs: the orbital gradient (the DIIS error), the
Hessian-vector product (A + B) x at one J/K build (UHF: one two-density build) on every J/K path, exact KS singlet response, the
exp(K) rotation and the SCF restart from a density. An augmented-Hessian / trust-region Newton step (PySCF `newton()`
style: solve [[0, g^T], [g, H]] for the lowest root with the existing Davidson, scale to a trust radius, accept by the
energy, switch from DIIS once ||g|| < 1e-2) is about 200-300 lines in `qm_stability.cpp`-like code plus a loop driver.
Costs: 3-8 J/K builds per macro iteration, so for well-behaved molecules it is slower than DIIS (12-15 iterations);
the benefit is for cases DIIS cannot converge, of which this test set had none (only undamped Roothaan failed).
Open points before implementing: UKS needs the spin f_xc (the central difference is too noisy for a Newton solver
near 1e-8), the incremental integral-direct builds and RI mixed precision assume the DIIS loop, the device loop would
stay out. Realistic for HF/RKS on the CPU; not started.
