# Native KS-DFT Implementation Guide

**Document Purpose**: Implementation guide for the native Kohn-Sham DFT engine in Curcuma
(HF / LDA / PBE / B3LYP), ported/extended from the xcDFT TCC winter school 2019 code with
ORCA 6.1 as the validation reference.

**Status**: WP1 + WP2 + **WP3 working** (July 2026). The closed-shell HF SCF now
converges and reproduces ORCA 6.1 HF/def2-SVP on all 10 validation molecules —
9 of them to ≤4e-9 Eh with ORCA's default guess, BH on the same SCF branch (see
"What was wrong until July 2026" below; three real kernel defects were found and
fixed). `ctest -L dft_1e` + `ctest -L dft_2e` 20/20. LDA/PBE/B3LYP (WP5-WP7) and
the analytic gradient (WP8) are still open: those methods intentionally return the
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
All three defects lived in `dft_integrals.cpp` and were invisible to the WP1/WP2
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
against ORCA belongs in `dft_1e` from the start.

### WP3 limits / not implemented

Only closed-shell (even electron count); no open-shell, no charged systems beyond
what the closed-shell code path supports. Two guesses (`sad`, `h0`) and DIIS or
plain damping -- no SOSCF, and no verification that a converged solution is a
minimum (a secondary solution could still be reported). No dispersion (D3/D4)
coupling, no solvation, no gradient (`hasGradient() == false` until WP8). Basis
scope is H-Ne (`def2-SVP.dat`).
The `dft_1e`/`dft_2e` graders' ORCA reference JSONs carry MO energies only, so the
new total-energy agreement is currently checked by hand, not in CI.

---


---

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

The functional is fixed by the method name and mapped in `MethodFactory::create()` to
`DFTMethod(DFTFunctional::..., config)`. It is **not** a parameter.

Shared parameters live in the `dft` parameter module (`-dft.basis`, `-dft.grid`,
`-dft.scf_*`); see `dft.h` (`BEGIN_PARAMETER_DEFINITION(dft)`).

---

## Roadmap

Full plan: [`docs/DFT_ROADMAP/`](DFT_ROADMAP/) (README + WP0..WP9) and the master plan file
`/home/conrad/.claude/plans/wir-wollen-dft-in-clever-gizmo.md`.

| WP | Title | Status |
|----|-------|--------|
| WP0 | Geruest, Quellenangabe, Setup | WP0 scaffold done (this doc) |
| WP1 | GTO 1e integrals (S/T/V) | done — `ctest -L dft_1e` 10/10 (kernels corrected Jul 2026) |
| WP2 | 4-centre ERI (McMurchie-Davidson) | done — `ctest -L dft_2e` 10/10 (kernels corrected Jul 2026) |
| WP3 | HF-SCF (rung 666, hard ERI gate) | **done (Jul 2026)** — 10/10 molecules ≤4e-9 Eh vs ORCA, see above |
| WP4 | DFT grid (Euler-Maclaurin + Lebedev + Becke) | open |
| WP5 | LDA (Slater-Dirac + VWN5) | open |
| WP6 | PBE-GGA (+ grad rho) | open |
| WP7 | B3LYP hybrid (+ exact exchange) | open |
| WP8 | Analytic gradient | open |
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
- `dft_integrals.{hpp,cpp}` -- `ERITensor` (flat n^4, chemists' (mu nu | lam sig),
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
- `test_cases/dft_2e/` -- `dump_dft_2e.cpp` (emits ERI/J/K + dummy P as JSON),
  `diff_dft_2e.py` (pure-stdlib orchestrator, three gates), `CMakeLists.txt`
  (registers `ctest -L dft_2e` for the 10 dft_1e molecules, tol 1e-10).
- `scripts/dft_2e_python_ints.py` -- independent pure-stdlib (no numpy/scipy/pyscf)
  MD ERI witness sharing curcuma's exact cartesian AO order, plus J/K from the
  same dummy density.

**Tested**: `ctest -L dft_2e` 10/10 on H2, He, LiH, BeH2, BH, CH4, NH3, H2O, HF, Ne
(def2-SVP, cartesian 6d). Gates per molecule: (a) element-wise ERI vs the Python MD
witness max|curc-witness| = 2.1e-14 (worst CH4/NH3/H2O; He/LiH exact 0); (b) 8-fold
ERI symmetry exact (0.0), J/K symmetric to 1e-12, Tr(P*S)=2, Tr(P*J)==Tr(P*K) to
1.8e-13 (worst Ne), J/K element-wise vs witness to 7e-14; (ss|ss) primitive spot-check
vs the closed form `2 pi^{5/2}/(p q sqrt(p+q)) K_AB K_CD F_0(T)` = 2.2e-16 (incl. T=0
same-centre). `ctest -L dft_1e` regression still 10/10; `make -j4` and
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
| `src/core/energy_calculators/qm_methods/dft.h` | `DFTFunctional` enum, `dft` PARAM block, `DFT : public QMDriver` |
| `src/core/energy_calculators/qm_methods/dft.cpp` | Engine: nuclear repulsion + scaffold notice |
| `src/core/energy_calculators/qm_methods/dft_method.h` | `DFTMethod : public ComputationalMethod` wrapper |
| `src/core/energy_calculators/qm_methods/dft_method.cpp` | Wrapper delegation + `getDefaultConfig` |
| `src/core/energy_calculators/method_factory.cpp` | `hf`/`lda`/`pbe`/`b3lyp` dispatch + listings |
| `src/core/energy_calculators/qm_methods/dft_integrals.hpp` | WP1 1e + WP2 ERI/J/K + spherical-transform API (`namespace dft1e`) |
| `src/core/energy_calculators/qm_methods/dft_integrals.cpp` | WP1 1e kernels + WP2 MD ERI / Coulomb / exchange / 4-index spherical |
| `test_cases/dft_2e/dump_dft_2e.cpp` | WP2 dumper: ERI (cartesian 6d) + J/K + dummy P as JSON |
| `test_cases/dft_2e/diff_dft_2e.py` | WP2 orchestrator: kernel gate + internal consistency + optional xcDFT |
| `test_cases/dft_2e/CMakeLists.txt` | `ctest -L dft_2e` registration (10 mols, reuses dft_1e .xyz) |
| `scripts/dft_2e_python_ints.py` | Independent pure-stdlib MD ERI + J/K witness (shared AO order) |