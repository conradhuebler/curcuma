# Native ab-initio QM Implementation Guide (HF / HF-3c / KS-DFT)

**Document Purpose**: Implementation guide for the native Gaussian-basis QM engine in
Curcuma (HF, HF-3c; LDA / PBE / B3LYP pending), ported/extended from the xcDFT TCC winter
school 2019 code with ORCA 6.1 as the validation reference.

**Naming (Sep 2026)**: the engine was renamed from `DFT` to **`QMEngine`** because HF and
HF-3c, not only DFT functionals, run on it. Files `dft*.{h,cpp}` -> `qm_engine`, `qm_scf`,
`qm_integrals`, `qm_method`; namespace `dft1e` -> `qmint`; `DFTFunctional` -> `QMFunctional`;
parameter module `dft` -> **`qm`** (`-qm.basis`, `-qm.scf_*`, `-qm.threads`,
`-qm.eri_screening`; `-dft.*` is still merged as a legacy scope); `CURCUMA_DFT_BASIS` ->
`CURCUMA_QM_BASIS` (old name still read); tests `dft_1e`/`dft_2e` -> `qm_1e`/`qm_2e`.
The historical sections below keep the names that were current at the time.
Performance work and the GPU plan: [QM_GPU_ROADMAP.md](QM_GPU_ROADMAP.md).

**Status**: WP1 + WP2 + **WP3 working** (July 2026). The closed-shell HF SCF now
converges and reproduces ORCA 6.1 HF/def2-SVP on all 10 validation molecules —
9 of them to ≤4e-9 Eh with ORCA's default guess, BH on the same SCF branch (see
"What was wrong until July 2026" below; three real kernel defects were found and
fixed). `ctest -L qm_1e` + `ctest -L qm_2e` 20/20. The analytic HF gradient (WP8)
exists since Sep 2026 (see the WP8 section; `-opt` works for `hf` and `hf-3c`).
LDA/PBE/B3LYP (WP5-WP7) are still open: those methods intentionally return the
nuclear repulsion only. ⚠️ AI-generated / ⚙️ machine-tested only — no human
production test, and the validation set is 10 molecules in a H-Ne basis.

## WP3 -- HF-SCF vs ORCA 6.1 (July 2026)

`curcuma -sp <mol>.xyz -method hf` (def2-SVP, default spherical 5d, default SAD
guess), against `orca ! HF def2-SVP TightSCF`:

| molecule | curcuma | ORCA | diff |
|---|---|---|---|
| H2 | -1.12890672 | -1.128906716206 | -3.8e-09 |
| He | -2.85516048 | -2.855160479342 | -6.6e-10 |
| LiH | -7.97866220 | -7.978662196021 | -4.0e-09 |
| BeH2 | -15.75349906 | -15.753499058370 | -1.6e-09 |
| CH4 | -40.16917039 | -40.169170388185 | -1.8e-09 |
| NH3 | -56.14009755 | -56.140097551803 | +1.8e-09 |
| H2O | -75.95751380 | -75.957513799019 | -9.8e-10 |
| HF | -99.93249914 | -99.932499138781 | -1.2e-09 |
| Ne | -128.37640681 | -128.376406809861 | -1.4e-10 |
| BH | -25.09917485 | -25.099174849888 | -1.1e-10 |

The residual (≈1e-11 relative) does not respond to a tighter SCF threshold
(`-dft.scf_threshold 1e-11` leaves it unchanged), and ORCA's own
TightSCF↔VeryTightSCF shift is only 0-2e-10, so it is **curcuma's**, not ORCA's
convergence: ≈4e-9 Eh is the honest agreement figure, and where exactly it comes
from is not yet pinned down (accumulated roundoff in the Fock/eigensolve path, or
the Boys/Hermite kernels at the 1e-11 level, are the candidates).

**BH and the initial guess.** BH has more than one closed-shell RHF stationary
point: ORCA's default model-potential guess finds -25.099174849888 and ORCA with
`! HCore` finds -24.871630934. The fix is the guess, not the SCF: the default is
now **`-dft.scf_guess sad`** (superposition of atomic densities), which reaches
ORCA's lower solution reproducibly (-25.09917485, 1.1e-10) and leaves every other
molecule bit-unchanged. Implemented in `buildAtomicGuess()`: per atom, diagonalize
that atom's block of the core Hamiltonian in its own atomic basis and fill the
atom's electrons into the resulting atomic orbitals (aufbau, fractional occupation
of the last when Z is odd), then sum. `Tr(P S) = N` by construction (rescaled for
charged systems).

`-dft.scf_guess h0` restores the bare-core start and usually reaches the *higher*
solution (-24.871630934, i.e. ORCA's `! HCore` branch, matching to 5e-9), but it is
**not a reliable way to select a BH solution**: the zero-density start sits near a
bifurcation there, and both outcomes were observed for the *identical* command and
binary -- 8/8 runs at -24.871630934, and separately a 10-run streak at
-25.09917485, at every threshold, mixer and thread count tried (the flag itself was
verified to arrive: verbosity 3 prints `SCF guess: h0 (Tr(PS) = 0.000000)`). SAD
reached the lower solution in **every** observation, so use `sad` and treat the h0
branch's BH basin as unreproducible. No other molecule in the set shows this: they
have a single solution and are stable run to run (H2O 5/5, all ten at ≤4e-9). The
bistability is a BH property plus floating-point detail; whether MKL's dynamic
thread count is the trigger was not established.

**Note (fixed July 2026)**: `DFTMethod` used to merge only the controller's *top*
level into its defaults, so **every `-dft.*` flag was silently ignored** --
`-dft.scf_mode`, `-dft.scf_threshold` and `-dft.scf_guess` all had no effect. The
CLI routes them into `controller["dft"]` (the `dft` module scope), which is now
merged explicitly, top level first, matching the other method wrappers.

### What was wrong until July 2026 (three kernel defects)

The WP3 SCF loop was never the problem: an independent Python RHF written on the
WP1/WP2 witness integrals reproduced curcuma's iteration trace *bit for bit*
(including its period-2 oscillation), so the defect had to be in the integrals.
All three defects lived in `qm_integrals.cpp` and were invisible to the WP1/WP2
gates because the "independent" Python witness shared the same algorithms.

1. **`boysArray` -- the dominant error.** The downward recursion started at a fixed
   `M = maxN + 25`, which is only enough for small T. For T >~ 15 it collapses:
   `F_0(37)` came out 50x too small and `F_0(T >= 50)` exactly 0. Every integral
   whose Gaussian pair sits away from the nucleus was therefore wrong -- i.e. every
   real molecule (a single atom has T = 0, which is why the atom-only cases looked
   fine). Fix: `F_0(T) = 0.5 sqrt(pi/T) erf(sqrt(T))` plus the **upward** recursion
   for T >= 1 (stable there; its subtraction cancels at small T). Worst relative
   error over T in [1e-14, 800], n <= 5: 4.7e-14 (was 1.0).
2. **`hermiteCoeffs` -- wrong from t = 2 on.** The "raise t" recursion reproduced
   only the t <= 1 coefficients: `E[2][0][0]` came out 0.5 instead of 0.0 and
   `E[2][1][1]` 0.25 instead of 0.0625. Fix: the standard McMurchie-Davidson
   forward recursion over i then j (Helgaker 9.5.5/9.5.6). Shared by the 1e
   nuclear attraction and the 4-centre ERI -- so the WP2 ERI was affected too.
3. **1e nuclear-attraction R auxiliary.** The base was missing its `(-2 gamma)^N`
   factor, and the displacement term had the wrong sign. On-centre both defects
   are invisible (P - C = 0, base index 0); off-centre they gave the right
   magnitude with an inverted sign, e.g. `<pz_A|1/r_B|s_A>` came out -0.238696
   instead of +0.238695. Fix: base `(2 pi/gamma)(-2 gamma)^N F_N(T)` and the PLUS
   sign of the bra recurrence (matching `buildRblockERI`).

Evidence for the fixes: exact closed forms (`<px|1/r|px>` = pi/6, formerly pi/3;
`<ss|1/r|ss>` = pi unchanged), deterministic 3D quadrature for off-centre
elements, and ORCA's own core Hamiltonian (He `h_pp` = +0.372308 Eh, and the
generalised (H,S) spectrum for BH/CH4 matching to ORCA's 6-decimal printout).
Before the fixes the HF energies were off by up to 19 Eh (H2O) and the SCF had no
fixed point at all.

**Lesson for the gates**: a witness that shares the algorithm cannot catch a
conceptual error in that algorithm. The WP1/WP2 gates validated the kernels
against a second implementation of the *same* recursions, and the one genuinely
independent reference they carried (ORCA, via the S eigenvalue and MO energies)
only covered quantities that were already right. A total-energy comparison
against ORCA belongs in `qm_1e` from the start.

### WP3 limits / not implemented

Only closed-shell (even electron count); no open-shell, no charged systems beyond
what the closed-shell code path supports. Two guesses (`sad`, `h0`) and DIIS or
plain damping -- no SOSCF, and no verification that a converged solution is a
minimum (a secondary solution could still be reported). No dispersion (D3/D4)
coupling, no solvation, no gradient at the time (added Sep 2026, see WP8). Basis
scope is H-Ne (`def2-SVP.dat`).
The `qm_1e`/`qm_2e` graders' ORCA reference JSONs carry MO energies only, so the
new total-energy agreement is currently checked by hand, not in CI.

---

## HF-3c (native, Sep 2026)

`HF-3c` (R. Sure, S. Grimme, J. Comput. Chem. 34, 1672 (2013)) is native since Sep 2026:
`-method hf-3c` -> `HF3CMethod` (`hf3c_method.{h,cpp}`). The external ORCA version is
reachable as `-method orca-hf-3c`.

    E(HF-3c) = E(HF/MINIX) + E_D3(BJ) + E_gCP + E_SRB

| term | implementation | parameters |
|---|---|---|
| HF/MINIX | `QMEngine(HF)` with the basis forced to `MINIX.dat` (ORCA export, H-Ne) | -- |
| D3(BJ) | `D3ParameterGenerator::createForHF3C()` | s6=1, s8=0.8777, a1=0.4171, a2=2.9149, **s9=0** (two-body only; ORCA prints only E6/E8) |
| gCP | `gcp.{h,cpp}`, ported from simple-dftd3 `gcp.f90` + `gcp/param.f90` (hf/minix) | sigma=0.1290, eta=1.1526, alpha=1.1549, beta=1.1763, `emiss`/`nbas` tables |
| SRB ("bas") | same file (`base` term) | `-qscal sum (Z_A Z_B)^1.5 exp(-rscal R0_AB^0.75 R_AB)`, qscal=0.03, rscal=0.7, R0 = D3 pair radii |

**Validation** (`ctest -L qm_hf3c`, 17 molecules: the 10 `qm_1e` ones plus benzene,
water dimer, formaldehyde, HCN, BF3, LiF, He...CH4): each term against programs that share
no code with curcuma -- PySCF (RHF/MINIX, basis parsed from the same file) and the
simple-dftd3 Python API (`DispersionModel` + `RationalDampingParam(method="hf3c")`,
`GeometricCounterpoise(method="hf3c")`), generated by `scripts/hf3c_reference.py`. Largest
deviations: HF 2.0e-11, D3 3.6e-11, gCP 4.3e-13, SRB 1.6e-12, total 2.6e-11 Eh. H2O
against ORCA 6.1 `! HF-3c`: -75.501489461789 vs -75.501489461360 (4.3e-10; ORCA's own
decomposition: HF -75.48881664772478, D3 -0.002646402235, gCP+bas -0.010026411). gCP/SRB
alone over all 55 H-Ne element pairs: energy <= 2.3e-12, gradient <= 2.3e-12 Eh/Bohr.

**Traps in the gCP port**, all in the source comments:
- The reference writes the Slater exponents as `[real(wp) :: 1.2000, ...]` **without
  `_wp`**, so Fortran reads them as single precision (2.5644 -> 2.5643999576...). Using the
  exact decimals shifts E_gCP by ~2e-8 relative (6e-10 Eh on F2); the port uses `float`.
- The ~1e-7 residual the earlier Python witness (`scripts/gcp_reference_witness.py`) showed
  against `otool_gcp` was **its own rounded `emiss` table** (5 digits instead of 6), not
  older tables in otool_gcp; with the exact values it matches simple-dftd3 to 7e-10.
- The like-exponent branch (|zeta_A - zeta_B| < 0.1 -- like-element pairs and He-C) uses
  the truncated 12-term B_n series, the other branch the closed form; He...CH4 exercises
  the He-C case.

**Not implemented / not tested**: closed shell only; H-Ne only (MINIX file and gCP tables);
no charged/open-shell reference comparison; not tested against ORCA beyond H2O; `-md` with
`hf-3c` technically possible (gradient exists) but not tested. Human production testing pending.

---

## WP8 -- Analytic RHF gradient, and `-opt` with `hf` / `hf-3c` (Sep 2026)

`-method hf` and `-method hf-3c` now have an analytic nuclear gradient, so `-opt` works.
🤖 AI-generated / ⚙️ machine-tested only.

**Formula** (closed-shell RHF, spin-summed density P; Pople, Krishnan, Schlegel, Binkley,
Int. J. Quantum Chem. S13, 225 (1979); Helgaker/Jørgensen/Olsen ch. 9):

    dE/dA = Σ P_mn d(T+V)_mn/dA − Σ W_mn dS_mn/dA + ½ Σ D_mnls d(mn|ls)/dA + dE_nn/dA
    W_mn  = 2 Σ_i^occ ε_i C_mi C_ni              (energy-weighted density)
    D_mnls = P_mn P_ls − ¼ (P_ml P_ns + P_ms P_nl)

Every derivative integral is built from the centre derivative of a cartesian Gaussian,
`d/dA_x g_l = 2α g_{l+1} − l g_{l−1}`, i.e. from ordinary integrals with a shifted power;
only derivatives with respect to the first function's centre are needed (symmetry of P and
of D), and the nuclear-attraction operator term follows from translational invariance.
The derivatives are taken in the cartesian basis and P/W are mapped back through the
(geometry-independent) 6d→5d transform. Code: `gradientOneElectron` /
`gradientTwoElectron` (WP8 block in `qm_integrals.cpp`), `QMEngine::computeGradient()`.
The 2e part uses the shell-blocked kernel with a two-step McMurchie-Davidson contraction.
HF-3c adds the D3 gradient (with the CN chain rule, same call pattern as the GFN1 D3) and
the gCP/SRB gradient. `getGradient()` returns Eh/Å (the `ComputationalMethod` contract).

**Validation** (`ctest -L qm_grad`, 19 tests):

| check | cases | result |
|---|---|---|
| analytic vs PySCF analytic RHF gradient (def2-SVP, with d shells) | 8 molecules incl. 2 strongly distorted | ≤ 1.9e-10 Eh/Bohr |
| analytic vs PySCF RHF/MINIX + simple-dftd3 D3 + gCP gradients (hf-3c) | 9 molecules (H2O…benzene, BF3, LiF, He···CH4) | ≤ 1.5e-10 Eh/Bohr |
| analytic vs central FD of curcuma's own energy (h = 1e-4 Bohr) | all 17 | ~3e-9 (FD truncation) |
| `-opt -method hf-3c` vs PySCF + simple-dftd3 minimum (scipy BFGS, |g| < 1e-7) | distorted H2O, distorted formaldehyde | ΔE ≤ 3e-12 Eh, RMSD ≤ 1e-6 Å |

Also checked by hand (not a ctest): `-opt -method hf` (def2-SVP) on distorted formaldehyde
against the PySCF minimum: RMSD 1.0e-5 Å, energy equal to 8 decimals.

**Cost**: the 2e gradient runs over the canonical shell quartets only (as the ERI build),
with derivatives on the centres of A, B and C and the one on D from translational invariance;
one-centre quartets are skipped. Benzene HF-3c: ERI 1.14 s / 2e gradient 4.5 s on 1 thread,
0.35 s / 1.24 s on 4 threads (first version: 8.5 s / 2.2 s, every ordered bra pair). Benzene
HF/def2-SVP (d shells, 4 threads): ERI 6.2 s, 2e gradient 12.3 s; there the gradient agrees
with PySCF to 3e-11 Eh/Bohr at `-qm.scf_threshold 1e-10` but only 5e-9 at the default
1e-6 -- the gradient is first order in the SCF error (irrelevant for `-opt`, whose gradient
threshold is 5e-4). Density-weighted screening is still open.
Whole optimisations with this gradient, the SCF-DIIS below and the warm start (4 threads):
benzene HF-3c 42.3 -> 26.2 s, formaldehyde HF/def2-SVP 6.25 -> 4.19 s, same final energies.

**Found on the way**: `-qm.threads` was silently ignored -- the constructor accepted only
`is_number_integer()`, the CLI/registry deliver `1.0`, and QMDriver's default of 4 threads
stayed in place for every run. Earlier "4 threads" timings were therefore real; the default
is now the registry's 1 (or the global `-threads`).

**Not tested / not implemented**: open shell; elements beyond Ne; MD with `hf`/`hf-3c`
(energy conservation not checked); gradients of `lda`/`pbe`/`b3lyp` (no V_xc yet, they
report `hasGradient() == false`). Every optimisation step still rebuilds the ERI tensor.

**SCF warm start** (Sep 2026, `-qm.scf_warm_start`, default on): the SCF of a new geometry
of the same molecule starts from the previous converged occupied orbitals, re-orthonormalised
in the new overlap metric, `C' = C (C^T S C)^{-1/2}`, `P = 2 C' C'^T` (exact electron count
and idempotency). Falls back to `-qm.scf_guess` for the first geometry, a different molecule,
or a near-singular `C^T S C`. Measured: water stepped 5x by 0.01 A, 61 vs 90 SCF iterations,
energies identical to 7e-14 Eh; `-opt` formaldehyde HF/def2-SVP 246 -> 149 iterations,
benzene HF-3c 110 -> 86, same final energies. **Wall time barely moves** (6.45 -> 6.25 s,
42.9 -> 42.3 s on 4 threads): per step the ERI rebuild and the gradient dominate, not the SCF.
**SCF-DIIS (Sep 2026)**: the QM SCF has its own accelerator (`QMScfAccelerator` in
`qm_scf.cpp`; the shared `DIISAccelerator` stays as it is for native xTB). Pulay DIIS now
starts at the second Fock matrix (`-qm.diis_start 1`, was 3), keeps 8 entries
(`-qm.diis_subspace`), uses the commutator in the orthonormal basis `X^T(FPS-SPF)X`
(Pulay 1982) and drops the oldest entries when the B matrix is near-singular. Over 33
molecule/basis cases (H-Ne, def2-SVP and MINIX): 357 -> 325 SCF iterations at the default
threshold, 503 -> 413 at 1e-9, every energy unchanged to 8 decimals. `-qm.scf_mode adiis`
adds ADIIS (Hu & Yang 2010) far from convergence, blended into DIIS (Garza & Scuseria 2012);
on this set it has **no measured benefit** (342 / 427 iterations; with the bare-core start
348 DIIS vs 369 ADIIS, and BH lands on the same higher RHF solution with both), so it is
opt-in and only checked for reaching the same energies.

**Integral-direct SCF (Sep 2026)**: `-qm.scf_direct on|off|auto` (default `auto`).
`qmint::DirectJK` keeps only the shell-pair tables and Schwarz factors and recomputes every
screened canonical shell quartet per Fock build, digesting it straight into J and K
(Almloef, Faegri, Korsell 1982; one weight `v*deg/8` per canonical quartet, the eight index
orderings recovered by symmetrising J and K). Screening is density-weighted: a quartet is
skipped when `Q_AB Q_CD max|P|` over the six shell blocks it touches is below
`-qm.eri_screening` (Haeser & Ahlrichs 1989). The SCF builds incrementally,
`J(P_n) = J(P_ref) + J(P_n - P_ref)`, with a full build every 8 iterations and always a
full build for the final energy. Memory O(n^2) instead of n^4 doubles. `auto` stays with
the stored tensor (fastest J/K) as long as it fits in `-qm.eri_max_memory_mb` (4000 MB,
about 150 functions). Measured, 4 threads:

| system (def2-SVP) | functions | stored | direct |
|---|---:|---:|---:|
| benzene, energy | 114 | 5.0 s, peak 1320 MB | 28.4 s, peak 34 MB |
| naphthalene (auto -> direct) | 180 | would need 8.4 GB | 195 s, E = PySCF to 8 decimals (print precision) |

J/K with screening off equals the stored path to 4e-15 (`qm_direct_jk`); SCF energies agree
to 1e-13 Eh, including a warm-started geometry sequence, and `-opt` formaldehyde gives
the identical geometry. Density screening skipped 15 % (benzene) / 42 % (naphthalene) of
the quartets; loosening it to 1e-10 did not pay (more DIIS iterations). The per-build cost
is the integral kernel itself, so direct mode is for systems the stored tensor cannot hold,
not a speed-up. The 2e gradient was already direct.

**Faster ERI kernel (Sep 2026)**: four changes to `shellQuartet`, each checked element-wise
against the previous tensor (max |diff| 5e-15 on formaldehyde and water dimer, def2-SVP):
(1) the Boys function comes from a grid (spacing 0.05) with a 7-term Taylor step and the
downward recurrence instead of long-double `expl`/`erfl`; vs an independent series it is
accurate to 2.8e-15 relative (the old routine: 5.9e-14), `ctest -R qm_boys`;
(2) the ket Hermite sum runs one Cartesian direction at a time and is accumulated over the
ket primitives before the bra side is touched, with an s-type-ket shortcut;
(3) primitive screening -- primitive pairs are sorted by their own Schwarz factor and a
primitive quartet below `eri_screening * 1e-3` is skipped (off with `eri_screening 0`);
(4) the pair with more components goes on the bra side, `(ab|cd) = (cd|ab)`.
Profiling note: callgrind (instruction counts) pointed at the Hermite loops, but the
wall-clock gain came mostly from the Boys function and from the contraction order -- an
x87 `erfl` is cheap in instructions and expensive in cycles. Measured, def2-SVP:

| | before | after |
|---|---:|---:|
| benzene, one full J/K build, 1 / 4 threads (no screening) | 10.5 / 2.6 s | 3.4 / 0.94 s |
| benzene HF energy, stored tensor (ERI build) | 5.0 s (3.4 s) | 3.6 s (1.9 s) |
| benzene HF energy, integral-direct | 28.4 s | 8.6 s |
| naphthalene HF energy, integral-direct | 195 s | 58 s |
| `-opt` benzene HF-3c / formaldehyde HF | 26.2 / 4.19 s | 19.2 / 3.33 s |

All energies and optimised geometries unchanged at print precision. The 2e gradient kernel
only picked up the Boys change; the same reordering is open there.

Also fixed: `setMolecule()` with a *different* molecule on a reused method object kept the
old molecule's integrals (HF-3c water on an object last used for BH: -17.44 instead of
-75.50 Eh); `resetForNewMolecule()` now runs first. Both covered by `qm_update_geometry`.


## Attribution and Provenance

**Origin**: xcDFT -- TCCM winter school 2019: DFT (Fortran, RKS, LDA+HF, Euler-Maclaurin x
Lebedev grid, DIIS, Gaussian basis). Located at `/home/conrad/src/claude_curcuma/xcDFT`.

Directly ported modules carry the marker
`// Ported from xcDFT (TCCM winter school 2019: DFT), src/<file>.f90`.

**Curcuma-native extensions** (not present in xcDFT, to be implemented in later WPs):
- Becke multi-centre atomic grid (the xcDFT grid is single-centre / atom-only)
- GGA and hybrid functionals (xcDFT ships only LDA + HF stubs)
- The 4-centre ERI engine (xcDFT reads pre-computed integrals)
- Analytic nuclear gradient
- VWN5 correction of the (known-broken) xcDFT LDA correlation formula

**Second reference**: ORCA 6.1 (`/opt/orca_6_1`, reached via `ORCA_PATH=/opt/orca_6_1/orca`
through the existing curcuma-ORCA interface) with the def2-SVP basis, used to validate every
WP against an independent ab-initio implementation.

**Copyright / license header convention** (every new DFT source file):
```cpp
/*
 * <Native KS-DFT ...>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Ported from xcDFT (TCCM winter school 2019: DFT), src/<file>.f90
 * Second reference: ORCA 6.1 (def2-SVP).
 *
 * Claude Generated: ...
 *
 * This program is free software under GPL-3.0
 */
```

---

## Method Names (no umbrella)

There is **no** `-method dft`. Each functional is its own method name, exactly like
`gfn1` / `gfn2` / `pm3`:

| Method name | Functional | Rung | WP |
|------------|------------|------|----|
| `hf`       | Hartree-Fock (exact exchange only) | 666 | WP3 |
| `lda`       | Slater-Dirac + VWN5               | 1   | WP5 |
| `pbe`       | PBE-GGA (+ grad rho)             | 2   | WP6 |
| `b3lyp`     | B3LYP hybrid (+ exact exchange)   | 4   | WP7 |

| `hf-3c`    | HF/MINIX + D3(BJ) + gCP + SRB (composite, `HF3CMethod`) | -- | Sep 2026 |

The functional is fixed by the method name and mapped in `MethodFactory::create()` to
`QMMethod(QMFunctional::..., config)` (or `HF3CMethod`). It is **not** a parameter.

Shared parameters live in the `qm` parameter module (`-qm.basis`, `-qm.grid`,
`-qm.scf_*`, `-qm.threads`, `-qm.eri_screening`); see `qm_engine.h`
(`BEGIN_PARAMETER_DEFINITION(qm)`). The pre-Sep-2026 `-dft.*` scope is still merged
(below `-qm.*`).

---

## Roadmap

Full plan: [`docs/DFT_ROADMAP/`](DFT_ROADMAP/) (README + WP0..WP9) and the master plan file
`/home/conrad/.claude/plans/wir-wollen-dft-in-clever-gizmo.md`.

| WP | Title | Status |
|----|-------|--------|
| WP0 | Geruest, Quellenangabe, Setup | WP0 scaffold done (this doc) |
| WP1 | GTO 1e integrals (S/T/V) | done — `ctest -L qm_1e` 10/10 (kernels corrected Jul 2026) |
| WP2 | 4-centre ERI (McMurchie-Davidson) | done — `ctest -L qm_2e` 10/10 (kernels corrected Jul 2026) |
| WP3 | HF-SCF (rung 666, hard ERI gate) | **done (Jul 2026)** — 10/10 molecules ≤4e-9 Eh vs ORCA, see above |
| WP4 | DFT grid (Euler-Maclaurin + Lebedev + Becke) | open |
| WP5 | LDA (Slater-Dirac + VWN5) | open |
| WP6 | PBE-GGA (+ grad rho) | open |
| WP7 | B3LYP hybrid (+ exact exchange) | open |
| WP8 | Analytic gradient | **done for `hf`/`hf-3c` (Sep 2026)** — vs PySCF ≤1.9e-10 Eh/Bohr, `-opt` works (`ctest -L qm_grad`) |
| WP9 | Integration, CLI, parameters, docs, CTest | open |

---

## WP0 -- What was tested / not tested / not implemented

> AI-generated, machine-tested only. No `TESTED` / `APPROVED` label (human-only per CLAUDE.md).

**Implemented (WP0)**:
- `dft.h` / `dft.cpp` -- `class DFT : public QMDriver` with `DFTFunctional` enum (HF/LDA/PBE/B3LYP);
  `Calculation()` returns the nuclear repulsion energy only (electronic energy == 0) and logs
  "native DFT -- nur Geruest". `MakeOverlap` / `MakeH` are stubs.
- `dft_method.h` / `dft_method.cpp` -- `DFTMethod : public ComputationalMethod` wrapper that
  delegates to `DFT`; `hasGradient() == false` until WP8.
- `BEGIN_PARAMETER_DEFINITION(dft)` -- basis / grid / threads / scf_max_iterations /
  scf_threshold / scf_mode (functional deliberately absent).
- `MethodFactory` registration of `hf` / `lda` / `pbe` / `b3lyp` (always available, native).
- `kEnergyCalcMethodScopes` gains `"dft"`; CMakeLists.txt gains `dft.cpp` + `dft_method.cpp`.

**Tested**: `curcuma -sp <He.xyz> -method {hf,lda,pbe,b3lyp}` runs and prints the scaffold
notice with `E = E_nn` (0.0 Eh for a single He atom); `make -j4` and `make GenerateParams`
warning-free; `curcuma -methods` lists the four under Quantum Methods; full `ctest` green.

**Not tested**: any molecule with more than one atom at a physically meaningful energy
(scaffold returns only `E_nn`); comparison against ORCA/xcDFT (no electronic energy exists yet).

**Not implemented** (later WPs): GTO 1e integrals (T, V), 4-centre ERIs, the SCF loop, the XC
functionals, the quadrature grid, the analytic gradient, dispersion (D3/D4) coupling, solvent.

---

## WP2 -- What was tested / not tested / not implemented

> AI-generated, machine-tested only. No `TESTED` / `APPROVED` label (human-only per CLAUDE.md).

**Implemented (WP2)**:
- `qm_integrals.{hpp,cpp}` -- `ERITensor` (flat n^4, chemists' (mu nu | lam sig),
  `operator()` read / `at()` write / `set8()` 8-fold fill); `buildERI` via
  McMurchie-Davidson (reuses WP1 `hermiteCoeffs` + `boysArray`; new ERI R-auxiliary
  `R^n_{000}=(-2 rho)^n F_n(T)` with the (P-Q) displacement, plus the
  Hermite-Hermite Coulomb contraction over the ket pair with the (-1)^(tau+ups+om)
  sign; canonical quartet loop mu<=nu, lam<=sig, pair(mu nu)<=pair(lam sig), one
  compute + 8-fold fill); `buildCoulomb` J_munu = sum P_lamsig (mu nu|lam sig);
  `buildExchange` K_munu = sum P_lamsig (mu lam|nu sig); `applySphericalTransformERI`
  (4-index Q-transform for the WP3 SCF spherical path).
- `dft.h` / `dft.cpp` -- lazy `DFT::cartesianERI()` builds `m_eri_cart` on first
  request and caches it; NOT called by `Calculation()` (scaffold `-sp` path
  unchanged: E_nn + 0, electronic energy = 0). `eriReady()` query added.
- `test_cases/qm_2e/` -- `dump_qm_2e.cpp` (emits ERI/J/K + dummy P as JSON),
  `diff_qm_2e.py` (pure-stdlib orchestrator, three gates), `CMakeLists.txt`
  (registers `ctest -L qm_2e` for the 10 qm_1e molecules, tol 1e-10).
- `scripts/qm_2e_python_ints.py` -- independent pure-stdlib (no numpy/scipy/pyscf)
  MD ERI witness sharing curcuma's exact cartesian AO order, plus J/K from the
  same dummy density.

**Tested**: `ctest -L qm_2e` 10/10 on H2, He, LiH, BeH2, BH, CH4, NH3, H2O, HF, Ne
(def2-SVP, cartesian 6d). Gates per molecule: (a) element-wise ERI vs the Python MD
witness max|curc-witness| = 2.1e-14 (worst CH4/NH3/H2O; He/LiH exact 0); (b) 8-fold
ERI symmetry exact (0.0), J/K symmetric to 1e-12, Tr(P*S)=2, Tr(P*J)==Tr(P*K) to
1.8e-13 (worst Ne), J/K element-wise vs witness to 7e-14; (ss|ss) primitive spot-check
vs the closed form `2 pi^{5/2}/(p q sqrt(p+q)) K_AB K_CD F_0(T)` = 2.2e-16 (incl. T=0
same-centre). `ctest -L qm_1e` regression still 10/10; `make -j4` and
`make GenerateParams` warning-free; scaffold `curcuma -sp <mol>.xyz -method hf`
unchanged (E_nn only, ERI not built).

**Not tested**: atoms beyond Ne (def2-SVP.dat scope is H-Ne, matching WP1);
f-shells and general contraction (code-correct but untested); spherical-5d ERI
(`applySphericalTransformERI` provided, not exercised by the kernel gate which
runs cartesian 6d); open-shell / anionic / charged densities (dummy P is closed-shell
rank-1); the SCF loop and any real SCF density (WP3); the analytic gradient (WP8);
comparison against ORCA / xcDFT 2e output (the optional `--xcDFT-ref` He-VDZ gate is
prepared but the reference file is not in the repo, so it auto-skips); numerical
stability for large basis sets (>~31 basis functions, where the dense n^4 tensor and
the O(n^4 * primitives^4) MD contraction become expensive -- fine for WP2's ≤31 bf
target, not a production ERI engine).

**Not implemented** (later WPs): the SCF loop (WP3, will use J/K with a real
density and the closed-shell Fock H + 2J - K); the XC functional / quadrature grid
(WP4-WP7); the analytic gradient (WP8); dispersion (D3/D4) coupling; solvent. The
exchange factor for the closed-shell Fock (-K, or -0.5*K depending on convention) is
applied by the WP3 SCF, NOT by `buildExchange` (which returns the bare K with
E_exchange = 0.25*Tr(P*K)).

**Conservatism note**: agreement with the Python witness on 10 def2-SVP molecules
validates the MD kernel against an INDEPENDENT implementation (different code, same
algorithm) to ~1e-14, and the (ss|ss) spot-check anchors the prefactor/Boys base to
the analytic formula. This is strong but not absolute: both implementations share the
MD algorithm, so a conceptual MD error would pass both. The 8-fold symmetry and the
Tr(PJ)==Tr(PK) identity are internal-consistency checks (symmetry of the build;
rank-1 dummy-index relabeling), not independent correctness proofs. Human production
testing against an external 2e reference (ORCA/xcDFT/libint) is pending.

---

## Files

| File | Role |
|------|------|
| `src/core/energy_calculators/qm_methods/qm_engine.h/.cpp` | `QMFunctional` enum, `qm` PARAM block, `QMEngine : public QMDriver` (basis, 1e integrals, lazy ERI, geometry updates) |
| `src/core/energy_calculators/qm_methods/qm_scf.cpp` | closed-shell RHF SCF (SAD/H0 guess, DIIS), Fock build, active-basis ERI |
| `src/core/energy_calculators/qm_methods/qm_integrals.hpp/.cpp` | 1e kernels, shell-quartet-blocked MD ERI (Schwarz, OpenMP, on-the-fly 5d), J/K, spherical transforms (`namespace qmint`) |
| `src/core/energy_calculators/qm_methods/qm_method.h/.cpp` | `QMMethod : public ComputationalMethod` wrapper, `engineConfig()` scope merge |
| `src/core/energy_calculators/qm_methods/hf3c_method.h/.cpp` | `HF3CMethod`: HF/MINIX + D3(BJ) + gCP + SRB |
| `src/core/energy_calculators/qm_methods/gcp.h/.cpp` | gCP + SRB energy and gradient (simple-dftd3 port, H-Ne) |
| `src/core/energy_calculators/qm_methods/MINIX.dat`, `def2-SVP.dat` | basis sets (H-Ne) |
| `src/core/energy_calculators/method_factory.cpp` | `hf`/`lda`/`pbe`/`b3lyp`/`hf-3c` dispatch, `qm` scope |
| `test_cases/qm_1e/`, `test_cases/qm_2e/` | 1e / ERI gates vs Python witnesses (+ `bench_qm_integrals`, timing only) |
| `test_cases/qm_hf3c/` | HF-3c gate vs PySCF + simple-dftd3 + ORCA, `qm_update_geometry` regression test |
| `test_cases/qm_grad/` | gradient gate (analytic vs FD + PySCF/simple-dftd3) and `-opt` vs reference minimum |
| `scripts/qm_1e_python_ints.py`, `scripts/qm_2e_python_ints.py` | independent Python integral witnesses |
| `scripts/hf3c_reference.py`, `scripts/gcp_reference_witness.py` | HF-3c reference generator, standalone gCP witness |
| `scripts/hf_gradient_reference.py`, `scripts/hf3c_opt_reference.py` | reference gradients (PySCF + simple-dftd3), reference HF-3c minima |
