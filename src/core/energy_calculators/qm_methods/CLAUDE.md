# CLAUDE.md - Quantum Methods Directory

## Overview

The qm_methods directory contains the quantum mechanical method implementations and interfaces for Curcuma. This module provides a flexible, extensible framework for various quantum chemical methods including semi-empirical, tight-binding, and ab initio approaches.

## Structure

```
qm_methods/
├── interface/
│   ├── abstract_interface.h  # Base QMInterface class
│   ├── ulysses.cpp/h         # Ulysses semi-empirical interface
├── qm_driver.cpp/h           # Base driver for matrix-based QM methods
├── eht.cpp/h                 # Extended Hückel Theory implementation
├── eht_parameters.cpp/h      # EHT parameter database
├── gfnff.cpp/h               # Native GFN-FF implementation (gfnff)
├── gfnff_advanced.cpp/h      # Advanced GFN-FF features
├── xtbinterface.cpp/h        # XTB method interface
├── tbliteinterface.cpp/h     # TBLite method interface
├── dftd3interface.cpp/h      # DFT-D3 dispersion corrections
├── dftd4interface.cpp/h      # DFT-D4 dispersion corrections
├── STOIntegrals.hpp          # Slater-type orbital integrals
├── GTOIntegrals.hpp          # Gaussian-type orbital integrals
├── ParallelEigenSolver.hpp   # Parallel matrix diagonalization
├── basissetparser.hpp        # Basis set file parsing
└── QM_ARCHITECTURE.md        # Detailed architecture documentation
```

## Architecture Overview

### Interface Layer (`QMInterface`)
Unified polymorphic interface for all quantum methods:
```cpp
class QMInterface {
    virtual bool InitialiseMolecule() = 0;
    virtual double Calculation(bool gradient, bool verbose) = 0;
    virtual Vector Charges() const = 0;
    virtual Vector BondOrders() const = 0;
    virtual Geometry Gradient() const = 0;
};
```

### Driver Layer (`QMDriver`)
Base class for matrix-based quantum methods providing:
- Common matrix storage (Hamiltonian, overlap, MO coefficients)
- Threading support and parallel computation
- Template method pattern for customizable calculation steps

## Method Implementations

### Native Methods
- **Extended Hückel Theory (EHT)**: Complete semi-empirical implementation
- **Native GFN-FF (gfnff)**: Curcuma's own GFN-FF implementation (WORK IN PROGRESS)

### External Interfaces
- **XTB Interface**: Extended tight-binding methods (GFN-FF, GFN1, GFN2)
- **TBLite Interface**: Tight-binding DFT methods (GFN1, GFN2, iPEA1)
- **Ulysses Interface**: Various semi-empirical methods (PM3, AM1, MNDO)
- **ORCA Interface** (Jun 2026): External composite methods (HF-3c, B97-3c, r2SCAN-3c, PBEh-3c, custom via `orca`) — see ORCA Notes below
- **DFT-D3/D4**: Dispersion correction interfaces

### Integral Support
- **STO Integrals**: Analytical Slater-type orbital overlap integrals
- **GTO Integrals**: Primitive Gaussian overlap integrals
- **Parallel Solver**: Optimized eigenvalue solver with threading

## EnergyCalculator Integration

Methods are created by the polymorphic **`MethodFactory`** (`method_factory.cpp`) —
the old `SwitchMethod()` case-switch was removed in the Jan-2025 refactor. Each name
maps to a `ComputationalMethod` subclass: `gfn1`/`gfn2` -> `NativeXtbMethod`,
`gfnff` -> `GFNFFComputationalMethod`, `eht`/`pm3` -> native, `uff`/`qmdff` ->
`ForceFieldMethod`, `xtb-*`/`tblite-*`/`ipea1`/`ugfn*` -> external interfaces. See the
top-level CLAUDE.md "Supported Method Hierarchies".

## Configuration System

JSON-based configuration with common defaults:
- Threading control
- SCF convergence parameters
- Method-specific settings
- Debugging and verbosity options

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Quantum method development priorities, theoretical enhancements, and performance optimization goals to be defined by operator/programmer*

## Verbosity System Standards

### Universal Verbosity Levels (All QM/MM Methods)

#### Level 0: Silent Mode
- **ABSOLUTELY NO OUTPUT** - Critical for optimization and MD
- Only internal calculations, zero console output
- Exception: Critical errors that terminate calculation

#### Level 1: Minimal Results
- Final energy and convergence status only
- Brief method identification
- Essential warnings (convergence failures)

#### Level 2: Properties and Analysis  
- **SCF/Iteration Progress**: Energy per cycle with iteration count
- **Orbital Properties**: HOMO, LUMO, band gap (formatted tables)
- **Molecular Properties**: Charges, bond orders, dipole moment
- **Method Parameters**: Key settings and convergence criteria

#### Level 3: Complete Analysis
- **Full Orbital Lists**: All energies in structured tables
- **Coefficients**: MO contributions, population analysis  
- **Debug Details**: SCF parameters, convergence mechanisms
- **Method Internals**: Basis sets, grid details, algorithm specifics

### Implementation Guidelines
```cpp
// Standard pattern for all QM methods
if (CurcumaLogger::get_verbosity() >= 1) {
    CurcumaLogger::energy_abs(final_energy, "SCF Energy");
}
if (CurcumaLogger::get_verbosity() >= 2) {
    // Orbital properties, molecular analysis
}
if (CurcumaLogger::get_verbosity() >= 3) {
    // Full coefficient matrices, debug output
}
```

## Variable Section

### AI Implementation Status

> ⚠️ **All native QM methods: AI-implemented, machine-tested only — not human production tested.**

| Method | Test Status | Notes |
|--------|-------------|-------|
| EHT | not systematically tested | qualitative only |
| Native GFN1-xTB | vs tblite: **10/12 SQM molecules at 1e-8** (only He2 ~1.5e-8 floor, complex ~2.6e-7 open) | `gfn1` routes to `curcuma::xtb::XTB(GFN1)` (canonical; the `ngfn1` alias was removed). D3(BJ) dispersion computed (`createForGFN1`, s9=0). **WP2 — three fixes brought GFN1 from "1e-8 on none" to "10/12":** (1) **scrambled D3 reference-CN table** (`d3_reference_cn.cpp`) regenerated from s-dftd3 → C6 exact. (2) **electronic residual** — GFN1 third-order potential was double-counted: `addThirdOrderPotential` (`xtb_thirdorder.cpp`) added `v_at(i)=q_i²·Γ_i` to **both** `pot.v_at` and (broadcast) `pot.v_sh`, but `expand_potential` sums `v_ao=v_sh+v_at` → 2× in Fock → wrong SCF fixed point; same bug in CPSCF response (`xtb_response.cpp`). Fix: keep in `v_at` only. (3) **D3 C8 tail** — `getR6` (`d3param_generator.cpp`) used empirical `C8/C6=10.72·r4r2²` (+0.06%/pair); replaced with exact s-dftd3 `3·√(½·⟨r⁴⟩/⟨r²⟩·√Z)²`. triose now bit-matches tblite (electronic AND dispersion). **Validation is vs tblite/s-dftd3 via GFN1 only**; the shared D3 `getR6` fix benefits UFF-D3/`gfnff-d3` but those stay unvalidated. See [docs/SQM_WP2_gfn1_accuracy.md](../../../../docs/SQM_WP2_gfn1_accuracy.md). xfails reduced to He2+complex. |
| Native GFN2-xTB | revalidation pending | **Broyden charge mixing is the new default SCF (2026-05-29)** — `complex` (231 atoms) divergence fixed. Default `-scf_mode broyden` (`broyden_mixer.h`, modified Broyden / Johnson 1988) mixes the SCC charge vector like tblite; `-method gfn2` converges `complex` from bare-H0 with no flags (34 it → −329.52707823 Eh), energy-identical to `-scf_mode diis` on the small set. Other modes: `-scf_mode diis\|plain\|level-shift` (`applyLevelShift` in `xtb_scf.cpp`), `-scf_guess h0\|eeq` (`seedEEQGuess` in `xtb_native.cpp`). SCF params now wired (`scf_damping`/`diis_start`/`diis_subspace`/`level_shift`, `xtb` scope). Full native-GFN ctest suite green (gradient/CPSCF/D4/ngfn1+ngfn2-baseline, 7/7). See [docs/SCF_MODES.md](../../../../docs/SCF_MODES.md). **D4 dispersion integrated (AP7, 2026-05)** via `curcuma::dispersion::D4Evaluator` with Caldeweyher 2019 BJ params (s6=1, s8=2.7, a1=0.52, a2=5.0). Energy + analytical gradient (FD-validated < 5e-5 Eh/Å on H₂O/CH₄/NH₃/C₆H₆). **q-response now analytic for `d4_charge_source="eeq"` (default)**: ∂E_D4/∂q (Phase 1) × ∂q/∂x from single-shot dftd4 EEQ (`d4_charge_model`, Phase 2), validated by `test_d4_dedq`. `d4_charge_source="mulliken"` (CPSCF on GFN2 SCF) now meets the <5e-5 D4-gradient target for H₂O/HCN after a sign fix (`RHS_SIGN=+1`; the earlier −1 was a sign-blind-metric artifact that inverted the response). Multipole charge-response terms gated off pending re-derivation; raw ∂q/∂x is ~12-25% on polar molecules (immaterial to the D4 gradient). Gated by `test_xtb_cpscf` Tests C/D. **Alignment-Roadmap vs tblite (Komponenten-Audit, SCF-Status, D4-Tiefe):** [docs/GFN2_NATIVE_ROADMAP.md](../../../../docs/GFN2_NATIVE_ROADMAP.md). See [docs/PHASE3B4](../../../../docs/PHASE3B4_MULLIKEN_RESPONSE_WP.md) / [PHASE3B5](../../../../docs/PHASE3B5_MULTIPOLE_RESPONSE_WP.md) |
| PM3/AM1/MNDO | 21/21 vs Ulysses (< 4 µEh) | most complete |
| PM6 | not tested | parameters present, untested |

### Current Development Status
- **✅ Universal Verbosity**: All QM methods integrated with CurcumaLogger
- **✅ ConfigManager Integration**: All QM interfaces accept ConfigManager (Phases 3A-3C)
- **✅ EHT Implementation**: Functional with orbital analysis and verbosity control
- **✅ XTB/TBLite Interfaces**: Native library verbosity synchronized with CurcumaLogger
- **✅ Ulysses Interface**: Complete CurcumaLogger integration with SCF progress
- **✅ Native GFN1/GFN2 Gradients (AP4)**: `xtb_gradient.cpp` — repulsion + H0/Pulay + Coulomb + CN; `-opt` converges on H₂O, CH₄, NH₃
- **✅ GFN2 Multipole Gradient Schritt 1 (AP5)**: Section 5 in `xtb_gradient.cpp` — SD/DD/SQ direct interaction gradient + mrad/CN chain-rule (port of `get_multipole_gradient_0d`)
- **✅ GFN2 Multipole Integral Pulay (AP5b)**: `xtb_gradient.cpp` lines ~316–390 — d(dp_iat)/dR and d(qp_iat)/dR terms via origin-shift correction from `dD_dA`; `cgto_multipole_grad_transformed()` in `xtb_multipole_ints.hpp`
- **✅ Gradient unit fix**: `m_gradient /= au` (was `*= au`) in `xtb_native.cpp` — caused au² error (~72% wrong) at non-equilibrium geometries
- **⚙️ SCF performance (2026-06-01)**: deep single-core retiming + 5 fixes — energy-neutral, 45/45 native-GFN ctests green. complex(231) energy+gradient: **gfn1 2982→1221 ms (beats tblite 2562 & xtb 1367); gfn2 2944→1364 ms (beats tblite 1567; 1.39× xtb 979)**. Fixes: (1) **EEQ initial guess default** (`scf_guess=eeq`) — gfn1 35→16 it, gfn2 34→22 it; (2) **Cholesky reduction** of the generalized eigenproblem — cache L=chol(S), per-iter `dsygst`+triangular back-transform replaces the dense `S^{-1/2}` double-GEMM (`buildOrthonormalizer`/`solveEigen`; `m_X` is now column-major `Eigen::MatrixXd` for Fortran dsygst); (3) occupied-only **density GEMM** (`leftCols(ncol)`); (4) **GFN2 D4 q-response routed to analytic EEQ** (`d4_charge_source=eeq` now actually wired) — post-SCF 653→83 ms, replaces the 574 ms Mulliken CPSCF; (5) **scf_threshold 1e-6→1e-5 default** (energy bit-identical, fewer iters; MD/opt may tighten). Per-iter `-verbosity 3` breakdown + `scripts/sqm_bench.sh` harness added. See [docs/SQM_PERFORMANCE.md](../../../../docs/SQM_PERFORMANCE.md)
- **⚙️ SCF intra-molecule multi-threading (2026-06-01)**: a single large molecule (`-sp`/`-opt`/MD) fans the SCF over `-threads N` via the project `CxxThreadPool`; default serial + bit-identical. Auto-gated (`src/core/intra_parallel_context.h`): molecule-level batch workers suppress it (no N×N oversubscribe), size-guarded. Parallelised: overlap+H0, GFN2 multipole integrals, gradient (thread-local reduce), `buildFock`, the GFN2 **D4 ATM 3-body** (`computeATM`, 2026-06-05), AND the per-iteration generalized **eigensolve** (`MklThreadScope` bumps MKL threads only around `solveEigen`; the build links threaded `mkl_gnu_thread`). **complex/231 @ `-threads 8` (2026-06-05): gfn2 1262→488 ms (2.58×), gfn1 1148→345 ms (3.33×); solve-eigen gfn2 607→206 (2.95×), gfn1 731→212 (3.45×).** The eigensolve saturates at ~8 threads (16/32 regress — D&C `dsyevd` memory-bound) and is at MKL's practical ceiling; partial-diag/GPU are dead-ends/neutral. Energy+gradient bit-identical t1/t8 (serial path byte-unchanged). See [docs/SQM_THREADING.md](../../../../docs/SQM_THREADING.md)
- **⚙️ large_system_mode (2026-06-02)**: opt-in `-large_system_mode fragments|dc|sparse` scale the native GFN SCF past ~1000 atoms by locality; default `none` (dense) byte-unchanged. All reuse the dense `XTB` + tblite-validated `evaluateComponentsAtFixedDensity`. **`fragments`** (`xtb_fragment_scf.cpp`): bond-fragment SCF, summed E + block-diagonal gradient (only mode with a gradient); `-eigensolver` propagates per fragment. **`dc`** (`xtb_dc.cpp`): Yang divide-and-conquer (global Fock + sub-block diag + shared-µ + core-projection + Broyden mixing); energy-only; `-eigensolver` propagates per sub-block, but `purify`/`lobpcg` fall back to dense GES (need full spectrum for Fermi occupation); complex 291→0.1 mEh / polymer-1410 90.8→0.10 mEh monotonic vs buffer. **`sparse`** (`xtb_sparse.cpp`): non-orthogonal Palser-Manolopoulos density purification (0 K, gapped); `-eigensolver` ignored (sparse IS the eigensolver). **Hard rule:** `-large_system_mode=fragments|dc + -eigensolver=purify` requires `-electronic_temperature 0` (wrapper rejects in `setMolecule`). ctest `large_system_modes` (label `large_system`) green. See [docs/SQM_LARGE_SYSTEMS.md](../../../../docs/SQM_LARGE_SYSTEMS.md)
- **⚙️ SQM perf WP AP2+AP1-GPU (2026-06-04)**: **AP2** — `weightedC6Gfn2` split into `buildAtomRefW`+`contractC6Gfn2`, `D4Evaluator` hoists the per-atom D4 reference weights out of the O(N²) pair loop + ATM `c6` fill (bit-identical; GFN2 post-SCF E 73→66 ms on complex, ~10%). **AP1-GPU** — opt-in partial diagonalization `-scf_gpu_partial_diag` (default OFF): the GPU resident eigensolve computes only the lowest `nocc+~5%` eigenpairs (`cusolverDnDsyevdx`/`Ssyevdx`), sentinel-pads eps, auto-widens on a tiny-gap tail; energies identical to the full solve, but **measured net-neutral on an RTX 5080** (tridiagonalization dominates the dense eigensolve, same as the CPU `dsyevr` dead-end) so the default GPU path stays the full `cusolverDnDsyevd`. See [docs/SQM_PERF_OPT_WP.md](../../../../docs/SQM_PERF_OPT_WP.md)
- **⚙️ SQM GPU Stage 6 — fully device-resident GFN2 SCF loop (2026-06-05)**: the whole per-iteration loop body (potential→Fock→eigensolve→occupation→density→charges/moments→SCC energy→**device Broyden mixer**) runs on the GPU; the host polls only `max|Δq_sh|`+4 energy scalars (O(1))/step — eps/occ/pop_ao/moments/q_sh never cross the bus. New kernels `k_occupations`/`k_qsh_scatter`/`k_d4_build_refw`/`k_energy_*`/`k_broyden_*`/`k_maxabsdiff`; seam `GpuScfBackend::{supportsResidentLoop,beginResidentLoop,residentScfStep,residentLoopCharges,…}`; gated GFN2 device-potential+Broyden path, host loop is byte-unchanged fallback (GFN1 unaffected). Energy bit-stable vs CPU (complex/231 −329.52707822); `gpu_gfn2_validation` 12/12; gradient ~1e-10 Eh/Å at `-scf_threshold 1e-8` (device Broyden lands at a different point in the loose 1e-5 ball — energy stationary→1e-8, gradient needs tight SCF); `ctest -L gpu_scf` 31/31 (occupation/q_sh/d4refw/scc_energy/broyden component tests — re-authored 2026-06-05 after the 5 files were lost in a WIP restore; residuals ~0..1e-14); sanitizer 0 errors; non-CUDA `release/` 24/24. **Honest: residency/correctness milestone, NOT a measured `-sp` speed-up** (complex SCF 341 ms GPU vs 417 ms CPU but ≈ Stage-5 ~515 ms TOTAL; the eigensolve + residual O(1) host-pointer cuBLAS syncs ~5/iter dominate, not the removed transfers). See [docs/SQM_GPU.md](../../../../docs/SQM_GPU.md)
- **⚙️ Multi-step SCC extrapolation (2026-06-05)**: opt-in generalisation of the 1-step warm-start — predicts the new-step SCC vector (GFN1 `q_sh`, GFN2 `[q_sh; dp_at; qp_at]`) from a history (`m_scf_history`) of converged steps. `-scf_extrapolation none|aspc|gauss` (ASPC Kolafa-2004 binomial coeffs / least-squares polynomial; both `Σw_j=1` charge-conserving), `-scf_extrapolation_order`, `-scf_extrapolation_apply guess|xlbomd`. `guess` keeps the full SCF (safe, same fixpoint within `scf_threshold`). `xlbomd` (Phase 2, EXPERIMENTAL) = extended-Lagrangian MD: the SCC density is a time-reversibly propagated auxiliary (Verlet + Niklasson dissipation, `xlbomdCoefficients` K=3..7) seeding a **converged** corrector (a naive bare-map/no-convergence corrector is non-contractive for tight-binding and diverges — verified), so the energy is exact; value is low MD energy drift (unvalidated). `m_xlbomd_aux` deque; bootstrap seeds K+1 converged states, then propagate from D[P_n]. `scf_xlbomd_correctors` = min corrector iters (default 1, no effect). Helpers `packSccState`/`unpackSccState`/`extrapolationWeights` in `xtb_native.cpp`; push in `UpdateMolecule`, predict in `Calculation` guess block, reset in `InitialiseMolecule`. Default `none` byte-unchanged. **Caffeine smooth-trajectory (guess, thr 1e-8):** SCF iters gfn1 170→79(aspc)/72(gauss), gfn2 215→90/92 (~55% fewer), final E bit-identical; xlbomd 170→141/215→180, energy bit-identical. Cites `aspc`(Kolafa 2004)/`density_extrapolation`(Pulay-Fogarasi 2004)/`xlbomd`(Niklasson 2008) on use. `ctest -L scf_extrapolation` green; default-path 29/29. Caveat: ideal-case numbers (uniform steps); real opt smaller, MD closer. See [docs/SQM_SCF_EXTRAPOLATION.md](../../../../docs/SQM_SCF_EXTRAPOLATION.md)
- **🔧 Native GFN-FF (gfnff)**: Architecture complete, parameter debugging in progress
- **🤖 AI-generated ORCA Interface** (Jun 2026): Wrapper for ORCA composite methods via external process (popen). Supports HF-3c, B97-3c, r2SCAN-3c, PBEh-3c and custom input. Status: AI-generated, pending human production test (no ✅ TESTED label).

### Verbosity Integration Status ✅
- **✅ EHT**: Complete integration with `printOrbitalAnalysisVerbose()` for Level 2
- **✅ XTBInterface**: XTB native verbosity synchronized (XTB_VERBOSITY_MUTED/MINIMAL/FULL)
- **✅ TBLiteInterface**: TBLite context verbosity synchronized with CurcumaLogger levels
- **✅ UlyssesInterface**: Full SCF progress, orbital analysis, and molecular properties
- **✅ All Methods**: Silent mode (Level 0) for optimization/MD, debug mode (Level 3) available

### Active Issues

#### ConfigManager Migration
- **DFT-D3/D4 UpdateParameters**: 🟡 LOW PRIORITY - `UpdateParameters()` methods still accept JSON instead of ConfigManager
  - Affects: `dftd3interface.h/cpp`, `dftd4interface.cpp`, `dispersion_method.cpp`
  - Impact: Seldom used, low priority for migration
  - Solution: Add ConfigManager overloads with delegating JSON constructors (Pattern established in Phase 3B)

#### Native GFN-FF (gfnff)
- **Parameter Generation**: JSON serialization creates null values for some parameters
- **Placeholder Parameters**: Missing real GFN-FF force field parameters from literature
- **Validation**: Incomplete parameter consistency checks with external GFN-FF reference

**GFN-FF Implementation Architecture** (for complete checklist see `../ff_methods/CLAUDE.md`):

- ✅ **Parameter Generation** (GFNFF class in gfnff.cpp)
  - CN-dependent radii, EEQ charges, topology detection, hybridization
  - Methods: generateTopologyAwareBonds(), generateGFNFFDispersionPairs(), etc.

- ✅ **Term Calculation** (ForceFieldThread in ../ff_methods/forcefieldthread.cpp)
  - Multi-threaded energy/gradient calculations
  - Methods: CalculateGFNFFBondContribution(), CalculateGFNFFDispersionContribution(), etc.

- ✅ **All 7 Terms Implemented**: Bond, Angle, Torsion, Inversion, Dispersion, Repulsion, Coulomb

**GFN-FF Components** (topology-aware parameter generation):
- ✅ CN-Berechnung (Coordination Numbers) - D3 method in `gfnff.cpp:calculateCoordinationNumbers()`
- ✅ CN-dependent radii - `r0_gfnff[86]` and `cnfak_gfnff[86]` in `gfnff.cpp`
- ✅ Row-dependent electronegativity - `p_enpoly[6][2]` for 6 periods
- ✅ Hybridization detection - Topology-based in `gfnff.cpp:determineHybridization()`
- ✅ Hybridization bond-strength matrix - `bsmat[4][4]` for mixed hybridizations
- ✅ EEQ charges - `gfnff.cpp:calculateEEQCharges()` with angewChem2020 parameters
- ✅ EEQ charge-dependent corrections (fqq) - Sigmoid function for bond strength
- ✅ Ring strain (ringf) - `findSmallestRings()`, fringbo=0.020 for 3-6 membered rings
- ✅ XH bond corrections (fxh) - Element-specific: BH(+10%), NH(+6%), OH(-7%)
- ✅ CN-dependent heavy atom (fcn) - Coordination-based bond weakening for Z>10

#### Other Issues
- **Memory Optimization**: Needed for large basis sets (>1000 atoms)
- **Ulysses D3H4X/D3H+ Corrections**: Corrections calculated internally but not extractable via getter methods - energies remain identical with/without corrections

### Recent Major Achievements (January-October 2025)
- **🎯 Universal Verbosity System**: All QM methods support consistent 4-level output control
- **🔧 ConfigManager Integration**: Complete type-safe parameter system (Phases 3A-3C)
  - 240 parameters with PARAM macro definitions
  - All 8 QM/FF interfaces accept ConfigManager constructors
  - Type-safe parameter access: `config.get<int>("accuracy")`
  - End-to-end parameter flow from CLI to external libraries
- **🔧 Native Library Integration**: XTB and TBLite verbosity controlled by CurcumaLogger
- **📊 Scientific Output**: HOMO/LUMO analysis, orbital properties, molecular analysis at Level 2
- **🚀 Performance**: Zero overhead silent mode for iterative calculations
- **🏗️ Enhanced Error Handling**: All methods use CurcumaLogger for consistent error reporting
- **✅ Complete Ulysses Integration**: All 27 semi-empirical methods (9 base × 3 correction modes) functional in Curcuma
- **🔧 MethodFactory Enhancement**: Fixed AM1/PM3 method recognition and universal method calculation support

### Performance Optimizations
- Threading support in matrix operations
- Efficient integral calculation algorithms
- Memory-optimized basis set handling
- **Silent Mode**: Zero-overhead Level 0 for iterative calculations

### ORCA Interface Notes (Jun 2026)
- **🤖 AI-generated** — pending human production test (no ✅ TESTED label).
- **NOT thread-safe**: each thread needs its own instance with unique orca_basename; `OrcaMethod::isThreadSafe() == false`.
- **CG atoms (element 226)**: rejected by default; opt-in via `orca_allow_cg=true`.
- **JSON schema**: validated against orca_2json (ORCA 5.0.4); unknown schemas fall back to text parsing.
- **Parser supports ORCA 5.x and 6.x** output formats (header variants: `The cartesian gradient`, `CARTESIAN GRADIENT`, `Gradient (Eh/Bohr)`).
- **Legacy CLI** (`curcuma -orca <input>`): preserved via `runExistingInput()`, which delegates to the new code path.
- **Tested**: unit-tested via `test_orca_interface` (8 sections, 30+ assertions). No integration test against a real ORCA binary in CI.

---

*This documentation covers all quantum mechanical methods and computational infrastructure. See QM_ARCHITECTURE.md for detailed technical specifications.*