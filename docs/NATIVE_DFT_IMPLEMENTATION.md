# Native KS-DFT Implementation Guide

**Document Purpose**: Implementation guide for the native Kohn-Sham DFT engine in Curcuma
(HF / LDA / PBE / B3LYP), ported/extended from the xcDFT TCC winter school 2019 code with
ORCA 6.1 as the validation reference.

**Status**: WP0 scaffold done. **WP1 done** (June 2026): real 1e GTO integrals
(S/T/V over contracted cartesian Gaussian, Obara-Saika 1986 + McMurchie-Davidson/Boys,
`def2-SVP.dat` for H-Ne); `MakeOverlap`→S, `MakeH`→Hc=T+V; `Calculation()` still
returns E_nn+0 (no SCF/XC yet). Validation: `ctest -L dft_1e` 10/10 (kernel vs
independent Python witness ≤1e-14, internal consistency, ORCA smallest spherical
S eigenvalue ≤1e-4). **WP2 done** (June 2026): 4-centre ERI engine
(McMurchie-Davidson, chemists' (μν|λσ), 8-fold symmetry) + Coulomb J / HF exchange K
from a dummy closed-shell density; `DFT::cartesianERI()` built lazily (NOT on the
scaffold `-sp` path). Validation: `ctest -L dft_2e` 10/10 (kernel vs independent
Python MD witness ≤2.1e-14, 8-fold symmetry exact, Tr(PJ)==Tr(PK) ≤1.8e-13,
(ss|ss) closed form 2.2e-16). All later work packages (WP3-WP8) fill in the SCF
loop, the XC functionals, and the gradient. ⚠️ AI-generated / ⚙️ machine-tested only.

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
| WP1 | GTO 1e integrals (S/T/V) | done — `ctest -L dft_1e` 10/10 |
| WP2 | 4-centre ERI (McMurchie-Davidson) | done — `ctest -L dft_2e` 10/10 |
| WP3 | HF-SCF (rung 666, hard ERI gate) | open |
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