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
- Report: energy decomposition `stability`, printed by `dump_qm_energy` as `stability <key> <value>`: `rhf_rhf`,
  `rhf_uhf`, `uhf_uhf` and `<key>_converged` (Davidson) of the FINAL solution, `initial_<key>` of the first analysis,
  `stable`, `followed` (Nov 5 fix: values of an earlier, followed solution used to stay next to `stable`). A Davidson
  solve that does not converge (60 iterations, residual 1e-6) is a warning; its value is an upper bound. Convention: eigenvalues of (A + B) in Eh; the RHF energy Hessian is
  4 (A + B). PySCF's `stability()` prints its own scaling; the comparison below uses a dense matrix of this convention.

Reference: `scripts/qm_stability_reference.py` (PySCF 2.14, dense (A + B) from PySCF's UHF response function with
exact f_xc, `numpy.linalg.eigh`; followed solution from PySCF `stability()` + SCF, for RHF -> UHF also
`followed_unit_*` from curcuma's start, unit-norm exp(K)), `test_cases/qm_stability/reference.json`. Independent of
that dense matrix: PySCF's `newton_ah.gen_g_hop_rhf` Hessian (dense) has lowest eigenvalue 2 x `rhf_rhf` on H2, water,
water 1.8 A, ozone, N2 (e.g. N2 -0.3836672132 / 2 = -0.1918336066).
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
B3LYP -1.012098168967 (dE 5.7e-11, <S^2> 0.69845775); ozone HF -224.142171083472 (dE -5.7e-14, <S^2> 0.89289086),
ozone B3LYP -225.139930993513 (UKS, <S^2> 0.19430292, 8.4e-4 Eh below RKS; PySCF from the same unit exp(K) start
-225.139930993529, dE 1.6e-11). PySCF's own `stability()` + UKS returns to the RKS solution for ozone/B3LYP (its rotation
is too small; a PySCF scan of the start amplitude t: t = 0.1 stays RKS, t = 0.3 ... 1.5 reach -225.139930993529), so its
`followed_energy` there is the RKS energy and the ctest uses `followed_unit_energy`. FD gate of the assembled surface after following (water 1.8 A,
HF, UHF after the RHF -> UHF switch, 4-point central difference of `-sp` energies with follow, 9 coordinates): max
3.2e-11 Eh/Bohr; level shift CH3 UHF b = 0.5: max 2.0e-10 Eh/Bohr.

Other observations: N2 at 2.0 A from the atom guess follows singlet -> triplet -> UHF internal to -108.546016 Eh
(<S^2> 2.685, three follows, final analysis stable with a zero eigenvalue, the rotation about the bond axis); PySCF
following only the external mode reaches another UHF minimum, -108.562699 Eh. Which solution following reaches depends
on the route; this case is not a ctest. RI-JK (H2 2.0 A): eigenvalues of the RI operator (0.388971 / -0.229505 vs
exact 0.388976 / -0.229587); no external RI reference. Integral-direct J/K gives the stored-tensor values to 4e-16.
Ozone B3LYP `check`: 1.3-1.8 s for the whole run at 4 threads (Nov 5 2026, re-measured; an earlier 1-thread figure of
4.3 min was not re-measured). Level shift: the commutator F'PS - SPF' of the shifted matrix equals that of F for an
idempotent P (S - 1/2 S P S commutes with P), so the shift changes the DIIS extrapolation and the diagonalisation only.

## Verification

Re-run Nov 5 2026 by a second session (the first session's logs were lost in a container restart; its numbers are
not reproduced here). 4-vCPU container, 15 GB RAM, no GPU, CPU build `release/` (USE_CUDA off).

- Bitwise with both options off: 11 runs through `dump_qm_energy --gradient --verbosity 1 --threads 1` with
  `CURCUMA_QM_SCF_TRACE=1` (water hf, b3lyp, hf integral-direct, hf `scf_mode plain`, hf RI-JK, hf-3c; CH3 UHF and UKS;
  stretched H2 2.0 A, N2 2.0 A, water 1.8 A hf): energy and gradient at 17 digits, every SCF iteration line (E_elec,
  dE, |g|, dP) and the iteration count; md5 of each file identical between a 029a72e build (separate worktree) and
  feature/qm with the fixes below, also with `-qm.scf_level_shift 0 -qm.scf_stability off` given explicitly.
  Densities are not dumped; the per-iteration dP and the gradient stand in.
- Level shift on vs off (same 11 runs): b = 0.5: |dE| <= 5.7e-14 Eh, max |dg| <= 5.3e-12 Eh/Bohr, iterations up
  (water 15 -> 26, water 1.8 A 17 -> 41, N2 2.0 A 14 -> 22, H2 2.0 A 6 -> 20, plain 64 -> 102); b = 0.2 with
  `scf_level_shift_off 1e-2`: |dE| <= 5.7e-14, max |dg| <= 7.9e-13, iterations unchanged except plain 64 -> 66, CH3 18 -> 19.
- FD gate (central difference at +-h and +-2h, h = 1e-3 A, of `dump_qm_energy` energies vs `--gradient`, scf_threshold 1e-10):
  water 1.8 A HF with `follow` (UHF after the switch) 8.7e-11 Eh/Bohr; CH3 UHF b = 0.5 4.6e-11; ozone B3LYP with
  `follow` (UKS) 1.1e-10.
- `ctest -L "qm_1e|qm_2e|qm_grad|qm_stability"` 64/64 (23 qm_stability), `ctest -R method_catalog` 2/2.
- PySCF 2.14 (venv): `scripts/qm_stability_reference.py` reproduced every stored value of the first session bitwise
  (1 thread); the regenerated file differs only in the new `followed_unit_*` keys and in the last digits (<= 1.1e-14)
  from `eigh` instead of `eigvalsh`.
- A full `make` in `release/` does not fit this container (about 40 GB usable disk, every test executable is
  0.2-0.65 GB): the build stopped at the link step of unrelated tests (`No space left on device`, later an OOM kill
  of cc1plus); `curcuma_core` and the targets above were built with exit status 0.

## Not tested / not implemented

- Not tested: the device loops (no GPU in the container: the counted `qm_scf` fallback with level shift was never
  run; with `--gpu cuda` and no plugin the J/K fallback `qm_jk` fires first and the device loop is never reached;
  the `qm_scf` count is also emitted for KS runs without device XC, where the device loop would not have run either),
  ctest labels qm_uhf/qm_uks/qm_b3lyp/qm_dft_grad/qm_lmp2 in the Nov 5 re-run (their executables were not rebuilt), ROCm/Vulkan, molecules beyond 3 heavy atoms, other bases than def2-SVP/MINIX, charged systems, RI-JK against
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

Hook points (checked Nov 5 2026 against the code): in `runSCF` (`qm_scf.cpp`) after `conv.observe(...)`, where the
level shift sits now: once ||g|| is below a switch threshold, replace `solveFock(fock_use)` + `buildDensity` by a
micro-iteration on the current orbitals (`vstop.C_prev`) with the operator of `stabilityLowest` kind 0 (needs C and F
in the MO basis of the current density, and a gradient vector g = C_v^T F C_o, not the AO commutator), then
`rotateOrbitals` and the density from it. Same in `runUSCF` with kind 2. `ScfVerifiedStop`, `DirectFockCadence` (incremental builds; Newton products
need full builds of D, which `buildJK` already does) and the device loop (`runSCFDevice`, skipped as for the level
shift) are the places that would have to opt out. The Davidson in `qm_stability.cpp` would need a linear-solve or
augmented-Hessian variant; the operator itself is reusable as is.
