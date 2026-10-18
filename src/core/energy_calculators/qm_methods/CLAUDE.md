# CLAUDE.md - Quantum Methods Directory

Quantum-mechanical and semi-empirical method implementations and the `ComputationalMethod` wrappers around them.
The former long version (status log of 2025 and 2026, per-feature performance records, old verbosity text) is kept in
[docs/archive/QM_METHODS_NOTES_2026-10.md](../../../../docs/archive/QM_METHODS_NOTES_2026-10.md), with a list of the statements in it that are out of date.
Detailed architecture: `QM_ARCHITECTURE.md` in this directory and [docs/SQM_VALIDATION.md](../../../../docs/SQM_VALIDATION.md).

> ⚠️ **All native methods here are 🤖 AI-implemented and ⚙️ machine-tested only. None is ✅ TESTED.**

## Layout

- **Wrappers** (`ComputationalMethod`): `native_xtb_method.*` (gfn1/gfn2), `xtb_method.*` (external XTB), `xtb_gpu_method.*`, `xtb_hip_method.*`, `xtb_vulkan_method.*`, `gfnff_method.*` (adapter to the native GFN-FF in `../ff_methods/`), `eht_method.*`, `nddo_method.*`, `orca_method.*`, `dispersion_method.*`. External: `tbliteinterface.*`, `xtbinterface.*`, `ulyssesinterface.*`, `orcainterface.*`, `dftd3interface.*`, `dftd4interface.*`.
- **Native xTB** (`xtb_native.cpp/h` and the `xtb_*.cpp` group: `xtb_h0`, `xtb_scf`, `xtb_gradient`, `xtb_coulomb`, `xtb_multipole`, `xtb_thirdorder`, `xtb_response`, `xtb_dc`, `xtb_fragment_scf`, `xtb_sparse`, `xtb_precond`): `curcuma::xtb::XTB` implements `QMInterface` directly and does not derive from `QMDriver`. The only GFN1/GFN2 parameter tables are `parameters/gfn1_params.hpp` and `parameters/gfn2_params.hpp`.
- **Driver** (`qm_driver.cpp/h`): base class for the STO/GTO-basis methods EHT and NDDO (matrix storage, `MakeOverlap/MakeH` hooks). NDDO (`nddo.*`, `nddo_params.*`) covers MNDO/AM1/PM3/PM6.
- **Integrals and solvers**: `STOIntegrals.hpp`, `GTOIntegrals.hpp`, `STO_CGTO.hpp`, `xtb_ao_utils.hpp`, `xtb_multipole_ints.hpp`, `native_eigensolver.*`, `ParallelEigenSolver.hpp`, `basissetparser.hpp`, GPU contexts under `cuda/` and the plugin loader `../gpu_plugin.*`.
- **Interface** (`interface/abstract_interface.h`, `interface/ulysses.*`): `QMInterface` with `InitialiseMolecule()`, `Calculation(gradient, verbose)`, `Charges()`, `BondOrders()`, `Gradient()`.

Methods are created by `MethodFactory` (`../method_factory.cpp`, table-driven, see `../CLAUDE.md`); `gfn1`/`gfn2` map to `NativeXtbMethod`, `xtb-*`/`tblite-*`/`ipea1`/`ugfn*` to the external providers.

## Status per method (what was compared with what)

| Method | Compared with | Result |
|---|---|---|
| `gfn1`, `gfn2` native | tblite, xtb 6.7.1 | root `CLAUDE.md` "Current Capabilities", [docs/SQM_VALIDATION.md](../../../../docs/SQM_VALIDATION.md), [docs/GMTKN55_VALIDATION.md](../../../../docs/GMTKN55_VALIDATION.md), [docs/MOR41_VALIDATION.md](../../../../docs/MOR41_VALIDATION.md) |
| `pm3`/`am1`/`mndo` | Ulysses | 21 of 21 tests, below 4 µEh |
| Native ab-initio QM (hf/hf-3c/lda/pbe/b3lyp) | `hf` vs ORCA 6.1 10/10 ≤4e-9 Eh; `hf-3c` terms vs PySCF+simple-dftd3 ≤4e-11 Eh; gradient vs PySCF ≤2e-10 Eh/Bohr; LDA/PBE/B3LYP open | `QMEngine` (`qm_engine.*`, `qm_scf.cpp`, `qm_integrals.*`), `HF3CMethod`; closed shell, H-Ne; AI-generated, machine-tested. Details: [docs/NATIVE_QM_IMPLEMENTATION.md](../../../../docs/NATIVE_QM_IMPLEMENTATION.md#status-summary), GPU: [docs/QM_GPU_ROADMAP.md](../../../../docs/QM_GPU_ROADMAP.md) |
| `ccsd` native (WP13) | PySCF `cc.CCSD`, ORCA 6.1 | 7 molecules H2O..benzene, full + fc, def2-SVP: <= 3.5e-10 / 1.5e-9 Eh; energy only, CPU, H-Ne; not tested: other bases, charged or open-shell systems, > 114 basis functions |
| `lmp2` native (WP15) | curcuma RI-MP2 (same aux), ORCA 6.1 DLPNO-MP2 (qualitative) | complete spaces = RI-MP2 to <= 2e-12 Eh; defaults 3.7e-6..2.4e-5 Eh/atom (alkanes C4-C40, (H2O)4-20, acenes to anthracene), tetracene/pentacene need `-qm.lmp2_pair_cut 25`; energy only, CPU, H/C/O tested; not tested: other elements and bases, charged systems; human production test pending |
| `pm6` | none found | parameters present, no test in `ctest -N` |
| `eht` | none | not systematically tested, qualitative only |
| `orca` | ORCA 5.x/6.x output formats | see below |

## Invariants and traps

- The gradient contract of every method is Eh/Angstrom; native xTB converts with `m_gradient /= au` (`xtb_native.cpp`) (Known Issue #28).
- A NaN, unconverged or implausible SCF must be an error, never a result (`-scf_allow_unconverged true` opts out); the convergence test covers everything the mixer mixes, including the GFN2 multipole part (#9, #29).
- The tblite reference is the target for gfn1/gfn2; the two documented places where xtb differs are opt-in switches (`-xtb.sto6g_legacy_4sp`, `-xtb.d4_atm_cutoff`) (#27, #30).
- Element-indexed tables: transition-metal shells are ordered `[d,s,p]`, so parameters indexed by angular momentum must not be read by shell index (#5).
- Silent mode (verbosity 0) must produce no output except errors; verbosity levels and logger functions: [docs/LOGGING_SYSTEM.md](../../../../docs/LOGGING_SYSTEM.md). No numerical path may depend on the print level.
- GPU backends are runtime plugins; fallbacks to the CPU are counted and reported (`gpu_fallback.h`).
- Native QM on the GPU (Q-G1/Q-G2, `cuda/qm_gpu_{jk,eri}.cu`, seam `qm_gpu_backend.h`): the device ERI kernel reads only `qmint::buildEriGpuTables` output, so any change to the CPU shell-pair tables/screening in `qm_integrals.cpp` must keep that function in step; `ctest -L qm_gpu` compares the tensor element-wise (AI-generated, machine-tested; [docs/QM_GPU_ROADMAP.md](../../../../docs/QM_GPU_ROADMAP.md)).
- Integral-direct J/K on the GPU (Q-G3, `cuda/qm_gpu_direct.cu`, `-qm.gpu_direct`): kernels shared with Q-G2 in `cuda/qm_gpu_eri_kernels.cuh`; its screen must stay the expression of `qmint::DirectJK::build` (same order), `ctest -L qm_gpu` checks J/K and the quartet count against the CPU; J/K last bits are schedule-dependent (atomics) (AI-generated, machine-tested).
- Integral-direct incremental J/K builds have an absolute screening noise (up to `eri_screening` per quartet, independent of dP): ||dP|| floats at 1e-8..3e-7 for benzene..pyrene; `DirectFockCadence` (`qm_scf.cpp`, shared by host and device loop) switches to full builds after a stall and then drops the DIIS entries from incremental builds; without it a direct SCF asked for less than the floor never converges (AI-generated, machine-tested; [docs/QM_GPU_ROADMAP.md](../../../../docs/QM_GPU_ROADMAP.md) section 2).
- Device-resident RHF SCF loop (Q-G4, `cuda/qm_gpu_scf.cu`, seam `qm_gpu_scf_backend.h`, `-qm.gpu_scf`): the driver stays in `QMEngine::runSCF`/`runSCFDevice`, host and device share `QMScfAccelerator` (only inner products differ), so any change to the host loop (Fock cadence, DIIS, convergence test) must be made in both; `ctest -L qm_gpu` checks E and iteration count (AI-generated, machine-tested).
- Two-electron RHF gradient on the GPU (Q-G7, `cuda/qm_gpu_grad.cu`, seam `qm_gpu_grad_backend.h`, `-qm.gpu_grad`): reads `qmint::buildEriGradGpuTables` (cartesian, Hermite tables raised by one on both shells) and must stay the sum of `qmint::gradientTwoElectron` (same D, canonical quartets, translational invariance), so a change on one side goes into both; quartet screen at `eri_screening`/100; `ctest -L qm_gpu` (`qm_gpu_grad_*`) compares with the CPU 2e gradient (AI-generated, machine-tested).
- One-electron RHF gradient on the GPU (Q-G8, `cuda/qm_gpu_grad1e.cu`, seam `qm_gpu_grad1e_backend.h`, `-qm.gpu_grad1e`): reads the same `buildEriGradGpuTables` output (built once per gradient for both parts) and must stay the quantity of `qmint::gradientOneElectron` (P dT + P dV incl. Hellmann-Feynman - W dS, no E_nn); `qm_gpu_grad1e_*` compare with it at 1e-10 (AI-generated, machine-tested).
- Density fitting (WP10, `qm_ri.h`/`qm_ri.cpp`, `-qm.ri off|j|jk`, default off) is a method variant, not a replacement: `buildJK` dispatches, the exact J/K stays in `buildJKExact`; `cuda/qm_gpu_ri.cu` (`-qm.gpu_ri`) only runs the contractions; aux files `def2-universal-{JFIT,JKFIT}.turbomole` (H-Ne); `ctest -L qm_ri` vs PySCF/ORCA (AI-generated, machine-tested; [WP10 doc](../../../../docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md)).
- RI 3-index tensor is PACKED (Oct 2026): row Q of `RIJK::tensor()` (and the gradient density G) holds the lower triangle, column `qmri::pairIndex(m,k)`, n(n+1)/2 columns; contract with `packContractionWeights`, full blocks via `unpackRow`/`unpackSymmetric`; the device seam `QMGpuRI` takes the same layout (entry `curcuma_cuda_create_qm_ri_packed`).
- RI-JK speed (Oct 2026): heavy RI steps are our OpenMP over single-threaded level-3 BLAS calls (MKL sequential); CPU `RIJK::build` forms J and K in one pass (dsymm on the packed triangle + SYRK; device: unpack kernel + GEMM + SYRK), so its J equals `coulomb(P)` only up to rounding; `buildJK` passes sqrt(2) C_occ when P is bitwise the density `runSCF` built from it; `cuda/qm_gpu_ri_fit.cu` (`-qm.gpu_ri_fit`, default off) is the setup solve B = L^-1 T.
- RI gradient (`qm_ri.cpp` section 6): derivative of the RI energy (exact for the Cholesky metric, not for eigen with dropped eigenvalues); RI-J adds `qmint::gradientTwoElectronExchange` (declared in `qm_ri.h`, `qm_integrals.hpp` untouched); `cuda/qm_gpu_ri_grad.cu` (`-qm.gpu_ri_grad`) must stay the quantity of `RIJK::threeIndexGradient`; PySCF references need `conv_tol_grad` <= 1e-11; `ctest -L qm_ri` (AI-generated, machine-tested).
- Canonical RMP2 (WP11, `qm_mp2.*`, `mp2_method.*`, `cuda/qm_gpu_mp2.cu`, seam `qm_gpu_mp2_backend.h`): CPU and device run the same four transformations and the same [i][a][j][b] layout, so a change on one side goes into both; it needs an n^4 AO tensor (no integral-direct MP2); `ctest -L qm_mp2` vs PySCF at 1e-8 (AI-generated, machine-tested; [WP11](../../../../docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md)).
- MP2 gradient (`qm_mp2_grad.cpp`, derivation at its top; mp2/scs/sos, frozen core): n^4 AO tensor plus n^4 AO 2PDM on the host; Z-vector J/K use the tensor rebuilt after a direct SCF (direct-builder noise stalls it); `qmint::gradientTwoElectronMP2` is declared in `qm_mp2.h` and shares `gradientTwoElectronDriver` with the RHF gradient.
- `cuda/qm_gpu_mp2_grad.cu` (`-qm.gpu_grad`, fallback `qm_mp2_grad`) is a COPY of the Q-G7 kernel with another quartet weight: a change of the derivative scheme goes into both .cu files. PySCF's analytic MP2 gradient is ~1e-7 off its own FD, so the ctest reference is PySCF's 4-point FD; RI SCF + MP2 gradient is an error (AI-generated, machine-tested; [WP11](../../../../docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md)).
- RI-MP2 (`qm_rimp2.cpp`, `cuda/qm_gpu_rimp2.cu`, seam `qm_gpu_rimp2_backend.h`, methods `ri-mp2`/`scs-ri-mp2`/`sos-ri-mp2`): aux `-qm.ri_mp2_aux` (default def2-universal-JKFIT, no def2-SVP/C: a method variant, the output label says so); reuses the SCF RI tensor only for the same aux file; CPU and device use the same B[i][a][P] layout and energy expression; own `qmint::gemmABt` because the linked BLAS is single-threaded MKL; no gradient; `ctest -L qm_mp2` (`qm_rimp2_*`) vs PySCF at 1e-8 (AI-generated, machine-tested; [WP10 doc](../../../../docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md)).
- Canonical RCCSD (WP13, `qm_ccsd.*`, `ccsd_method.*`, method `ccsd`): CPU only, exact integrals, full MO tensor incl. v^4 (memory check `-qm.cc_max_memory_mb`); contractions keep the index strings of `scripts/local_ccsd/ccsd_canonical.py`, change both together; no gradient, `-qm.ri` is an error; `ctest -L qm_ccsd` vs PySCF at 1e-8 on H, C, N, O only, thread scaling not measured on an idle machine (AI-generated, machine-tested; [WP13](../../../../docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md)).
- Local PNO-RI-MP2 (WP15, `qm_lmp2.*`, `lmp2_method.*`, method `lmp2`): follows `scripts/local_mp2/lmp2_common.py` step by step (change both together); RI tensor via `QMEngine::correlationRI` (shared with RI-MP2); the DIIS error block must stay normalised (unscaled it stalls at max|R| ~1e-9); `ctest -L qm_lmp2`: exactness vs RI-MP2 at 1e-9 plus a self-generated regression reference (AI-generated, machine-tested; [WP15](../../../../docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md)).
- The Broyden seed hook (`BroydenMixer::setStepHook`, `-scf_precond frag`, `xtb_precond.h`) may change only the path, never the fixed point; with `-scf_precond none` nothing is set; the guard (`-scf_precond_guard`) only removes the hook (and with `best` restarts the scalar step from the best earlier input) (AI-generated, machine-tested; [docs/SCF_PRECOND.md](../../../../docs/SCF_PRECOND.md)).
- The device copy of M (resident CUDA loop, `XtbGpuContext::broydenSetPreconditioner`) uses the device packing of the SCC vector (moments atom-major, the host `packSCC` is component-major); `xtb_native.cpp` remaps the indices before the upload. Without M the device Broyden runs its unchanged kernels (AI-generated, machine-tested).

## ORCA interface (🤖 AI-generated, pending human production test)

- Not thread-safe: `OrcaMethod::isThreadSafe()` is false; each thread needs its own instance with a unique `orca_basename`.
- Element 226 (coarse-grained beads) is rejected unless `orca_allow_cg=true`.
- Output parsing for ORCA 5.x and 6.x; JSON read validated against orca_2json 5.0.4, unknown schemas fall back to text parsing.
- The legacy `curcuma -orca <input>` path is kept through `runExistingInput()`.
- Unit test `test_orca_interface` exists and **fails** in the 2026-10-01 run (cause not investigated); there is no integration test against a real ORCA binary.

## Open items

- `DFT-D3/D4` `UpdateParameters()` still takes JSON instead of `ConfigManager` (low priority; `dftd3interface.*`, `dftd4interface.cpp`, `dispersion_method.cpp`).
- Memory for large basis sets (> 1000 atoms): see `-large_system_mode` in [docs/SQM_LARGE_SYSTEMS.md](../../../../docs/SQM_LARGE_SYSTEMS.md).
- Ulysses D3H4X/D3H+ corrections are computed internally but not exposed through getters; energies are identical with and without them.
- Open GPU and gradient items: root `CLAUDE.md` "Open Items".

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Quantum method development priorities, theoretical enhancements and performance goals, to be defined by the operator.*
