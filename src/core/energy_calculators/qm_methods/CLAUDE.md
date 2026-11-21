# CLAUDE.md - Quantum Methods Directory

Quantum-mechanical and semi-empirical method implementations and the `ComputationalMethod` wrappers around them.
The former long version (status log of 2025 and 2026, per-feature performance records, old verbosity text) is kept in
[docs/archive/QM_METHODS_NOTES_2026-10.md](../../../../docs/archive/QM_METHODS_NOTES_2026-10.md), with a list of the statements in it that are out of date.
Detailed architecture: `QM_ARCHITECTURE.md` in this directory and [docs/SQM_VALIDATION.md](../../../../docs/SQM_VALIDATION.md).

> ⚠️ **All native methods here are 🤖 AI-implemented and ⚙️ machine-tested only. None is ✅ TESTED.**

## Layout

- **Wrappers** (`ComputationalMethod`): `native_xtb_method.*` (gfn1/gfn2), `xtb_method.*` (external XTB), `xtb_gpu_method.*`, `xtb_hip_method.*`, `xtb_vulkan_method.*`, `gfnff_method.*` (adapter to the native GFN-FF in `../ff_methods/`), `eht_method.*`, `nddo_method.*`, `orca_method.*`, `dispersion_method.*`. External: `tbliteinterface.*`, `xtbinterface.*`, `ulyssesinterface.*`, `orcainterface.*`, `dftd3interface.*`, `dftd4interface.*`.
- **Native xTB** (`xtb_native.cpp/h` and the `xtb_*.cpp` group: `xtb_h0`, `xtb_scf`, `xtb_gradient`, `xtb_coulomb`, `xtb_multipole`, `xtb_thirdorder`, `xtb_response`, `xtb_dc`, `xtb_fragment_scf`, `xtb_sparse`, `xtb_precond`, `xtb_periodic`): `curcuma::xtb::XTB` implements `QMInterface` directly and does not derive from `QMDriver`. The only GFN1/GFN2 parameter tables are `parameters/gfn1_params.hpp` and `parameters/gfn2_params.hpp`.
- **Periodic cell (WP-PERIODIC-XTB, Nov 2026, 🤖 ⚙️)**: `XTB::setPeriodicCell()` (orthorhombic, Gamma point) switches every kernel to its image-sum twin in `xtb_periodic.cpp` (`if (m_pbc)` at the kernel entry, cell-free statements verbatim; tblite cutoffs `kPbc*`): CN, repulsion, XB, D3, D4 (pairs, ATM, dE/dq), S/H0, gamma = point-charge `EwaldSum::matrix` + Klopman-Ohno remainder switched to 1/r at 10 Bohr (tblite). [docs/WP-PERIODIC-XTB.md](../../../../docs/WP-PERIODIC-XTB.md).
- Periodic cell, wrapper: `NativeXtbMethod` forwards `Mol::m_has_pbc`, `hasPeriodicEnergy()` true with a cell; refusals (error, never a vacuum energy) for solvation + cell, d-shell elements + cell, non-orthorhombic cells; `-gpu rocm|vulkan` + cell swaps to the CPU (counted fallback in `EnergyCalculator::setMolecule`, hard error in `XTB::Calculation`); gfn1 energy 3-7e-11 Eh and forces 8.6e-7 Eh/Bohr from tblite on the water dimer.
- Periodic cell on CUDA (G1): `XtbGpuComputationalMethod` sets `XTB::setPeriodicDeviceScf(true)`; the host builds the periodic S/H0/L, gamma and multipole integrals and uploads them (`begin`/`beginMultipole`), the device runs the resident Fock/eigensolve/density/moments; every image-free device kernel (integral build, CN, potential, resident loops, D3/D4, gradient) is gated `!m_pbc`; parity gate ctest `pbc_xtb_gpu_parity` (1e-8), plugin ABI 8.
- Periodic cell on CUDA (G2-G6, Nov 2026): ported stages run as self-contained `GpuScfBackend` calls from the periodic host kernels (`dispersionATMPeriodic` G4a, `multipoleMatricesPeriodic` G3, `multipoleGradientPeriodic` G5a, `beginDispersionImages`/`dispersionDedqPeriodic` G4b, `ewaldGammaPeriodic` G2, `energyWeightedDensity`/`cnChainGradientPeriodic`/`repulsionPeriodic` G5b-d), host loop on any failure, `CURCUMA_PBC_GPU_HOST=<tokens>` keeps a stage on the host; with a cell `beginBasis` still runs (CSR offsets for bitwise charge/moment sums).
- Periodic cell on CUDA (G1b): `beginComputedPeriodic` builds S/H0/L and the GFN2 AO multipole integrals on the device; host copies deferred (`m_pbc_sx_deferred`, gamma-only consumers use `ensureHostGamma`); `CURCUMA_PBC_GPU_HOST=ints` keeps the host build. G1c: `setHostLoopPseudo` allows the device pseudo-diagonalisation in the host-driven periodic loop (`pseudo` token off).
- Periodic cell on CUDA (E1/A1/M1/M2/C1): `setHostLoopDensityComplement` (P = 2 (S^-1 - C_v C_v^T), token `dcompl`), ATM per triple class (`triview` = old kernel), multipole gradient/virial k sums by structure factors (`mppair` = old pair kernel)., `d4CnFoldPeriodic` (D4 CN chain rule, token `cnfold`). H1: `h0PulayGradientPeriodic` (H0/Pulay gradient + virial over images, W built on the device, token `h0grad`). H2/H3: `ewaldGradientPeriodic`, `koGradientPeriodic`, `cnPeriodic` (tokens `ewgrad`, `kograd`, `cn`). H4/H5: `d3DispersionPeriodic`, `d4PairsPeriodic` (tokens `d3`, `d4pair`). D1: in-SCF D4 dE/dq one block per atom (token `dedq1`).
- Periodic gradient (`calculateGradientPeriodic`, P3): every vacuum section over the images, self images only through dE/dCN, `EwaldSum::gradient` for the point-charge part; FD gate 8e-9 Eh/A (ctest `pbc_xtb_gradient`). Periodic D4/ATM/dE/dq and the multipole gradient run on fixed row blocks (`D4Evaluator::runBlocks`, thread-count independent bits, stage T).
- Periodic virial (P4, `XTB::setWantVirial/virial/updateCell`): -(dE/dd)(x)d per image incl. self images (weights 1/2), Ewald/multipole strain terms; cell-scaling FD gate 1e-9 Eh (ctest `pbc_xtb_virial`); barostat and `compute_pressure` accepted through `EnergyCalculator`.
- Periodic GFN2 multipoles (`setupMultipolePeriodic`, `buildMultipoleInteractionMatricesPeriodic`): tblite's alpha/cutoffs, but the reciprocal charge-dipole factor 1 (tblite 2) and charge-quadrupole -1/3 (tblite +1); tblite 0.4.0's periodic multipole matrices are alpha dependent, so periodic gfn2 differs from tblite by up to 1e-3 Eh (8 A cell). Diagnostic env `CURCUMA_PBC_MP_*` ([docs/WP-PERIODIC-XTB.md](../../../../docs/WP-PERIODIC-XTB.md) P2b).
- **Driver** (`qm_driver.cpp/h`): base class for the STO/GTO-basis methods EHT and NDDO (matrix storage, `MakeOverlap/MakeH` hooks). NDDO (`nddo.*`, `nddo_params.*`) covers MNDO/AM1/PM3/PM6.
- **Integrals and solvers**: `STOIntegrals.hpp`, `GTOIntegrals.hpp`, `STO_CGTO.hpp`, `xtb_ao_utils.hpp`, `xtb_multipole_ints.hpp`, `native_eigensolver.*`, `ParallelEigenSolver.hpp`, GPU contexts under `cuda/` and the plugin loader `../gpu_plugin.*`.
- **GFN-FF GPU wrapper**: template `gfnff_gpu_method_impl.h` with traits `gfnff_gpu_method.h` (CUDA) and `gfnff_hip_method.h` (ROCm); ROCm-only switches (`rocm_dc6dcn_tables`, `rocm_device_schur`, WP9; `rocm_fixed_point`, WP10) reach the HIP workspace and EEQ solver through `GFNFFHipComputationalMethod` and `rocmGfnffOptionsSlot()` (set around `setMolecule`/`calculateEnergy`), not through the template.
- **Interface** (`interface/abstract_interface.h`, `interface/ulysses.*`): `QMInterface` with `InitialiseMolecule()`, `Calculation(gradient, verbose)`, `Charges()`, `BondOrders()`, `Gradient()`.

Methods are created by `MethodFactory` (`../method_factory.cpp`, table-driven, see `../CLAUDE.md`); `gfn1`/`gfn2` map to `NativeXtbMethod`, `xtb-*`/`tblite-*`/`ipea1`/`ugfn*` to the external providers.

## Status per method (what was compared with what)

| Method | Compared with | Result |
|---|---|---|
| `gfn1`, `gfn2` native | tblite, xtb 6.7.1 | root `CLAUDE.md` "Current Capabilities", [docs/SQM_VALIDATION.md](../../../../docs/SQM_VALIDATION.md), [docs/GMTKN55_VALIDATION.md](../../../../docs/GMTKN55_VALIDATION.md), [docs/MOR41_VALIDATION.md](../../../../docs/MOR41_VALIDATION.md) |
| `pm3`/`am1`/`mndo` | Ulysses | 21 of 21 tests, below 4 µEh |
| `pm6` | none found | parameters present, no test in `ctest -N` |
| `eht` | none | not systematically tested, qualitative only |
| `orca` | ORCA 5.x/6.x output formats | see below |

## Invariants and traps

- The gradient contract of every method is Eh/Angstrom; native xTB converts with `m_gradient /= au` (`xtb_native.cpp`) (Known Issue #28).
- A NaN, unconverged or implausible SCF must be an error, never a result (`-scf_allow_unconverged true` opts out); the convergence test covers everything the mixer mixes, including the GFN2 multipole part (#9, #29).
- Generic SCF knobs live in registry module `scf` (`scf_parameters.h`, defaults are "method default" markers); per-method defaults and precedence (`-xtb.*` > `scf` > default): [docs/SCF_MODES.md](../../../../docs/SCF_MODES.md). NDDO keeps its inline DIIS (DIISAccelerator moved energies by 1.8e-7 Eh).
- The tblite reference is the target for gfn1/gfn2; the two documented places where xtb differs are opt-in switches (`-xtb.sto6g_legacy_4sp`, `-xtb.d4_atm_cutoff`) (#27, #30).
- Element-indexed tables: transition-metal shells are ordered `[d,s,p]`, so parameters indexed by angular momentum must not be read by shell index (#5).
- Silent mode (verbosity 0) must produce no output except errors; verbosity levels and logger functions: [docs/LOGGING_SYSTEM.md](../../../../docs/LOGGING_SYSTEM.md). No numerical path may depend on the print level.
- GPU backends are runtime plugins; fallbacks to the CPU are counted and reported (`gpu_fallback.h`).
- ROCm GFN2 builds the SCF potential of its host loop on the device (`-xtb.rocm_device_potential`, default true; `HipScfBackend` in `xtb_hip_method.cpp`, `resSolvePotential` in `rocm/xtb_hip_context.hip`); a failed step returns to the host potential via `GpuScfBackend::devicePotentialFallback` (AI-generated, machine-tested; [docs/ROCM_RESIDENT_LOOP_WP.md](../../../../docs/ROCM_RESIDENT_LOOP_WP.md)).
- ROCm GFN2 and GFN1 run the whole SCF iteration on the device (`-xtb.rocm_resident_loop`, default true; `resStep` in `rocm/xtb_hip_context.hip`, GFN1 via `resBeginLoopGfn1`: SCC vector = q_sh, atom-resolved third order, no WP1 potential needed): SCC vector in the host `packSCC` order and D4 weights at the mixed input charges, so it follows the host-driven loop to rounding; a failed step restarts that SCF on the host-driven loop via `GpuScfBackend::residentLoopFallback` (AI-generated, machine-tested; [docs/ROCM_RESIDENT_LOOP_WP.md](../../../../docs/ROCM_RESIDENT_LOOP_WP.md)).
- ROCm open shell (n_alpha != n_beta) stays on the resident loop (`-xtb.rocm_resident_open_shell`, default true; `HipScfBackend::setResidentOpenShell`, `k_res_occupations_uhf`): two spin fillings over the same eigenvectors, summed, Newton sums in the host order (AI-generated, machine-tested; [docs/ROCM_RESIDENT_LOOP_WP.md](../../../../docs/ROCM_RESIDENT_LOOP_WP.md)).
- ROCm `-scf_precond frag` (GFN1 and GFN2) runs on the resident loop (`-xtb.rocm_resident_precond`, default true; `resBroyden`/`pcSet` in `rocm/xtb_hip_context.hip`): `HipScfBackend::residentSetPreconditioner` maps the CUDA atom-major indices that `xtb_native.cpp` produces back to the host `packSCC` order (AI-generated, machine-tested; [docs/ROCM_RESIDENT_LOOP_WP.md](../../../../docs/ROCM_RESIDENT_LOOP_WP.md)).
- ROCm native ALPB/GBSA in the SCF potential on the device (`-xtb.rocm_device_solvation`, default true; `solvApply` in `rocm/xtb_hip_context.hip`): GFN2 v_at += B q_at, GFN1 also + B cm5 (host passes it via `GpuScfBackend::beginSolvationGfn1`); solvation energy and gradient stay on the host (AI-generated, machine-tested; [docs/ROCM_RESIDENT_LOOP_WP.md](../../../../docs/ROCM_RESIDENT_LOOP_WP.md)).
- ROCm pseudo-diagonalisation on the resident loop (generic `-scf_pseudo_diag`, ROCm default off, unlike CUDA; `eigensolveStep` in `rocm/xtb_hip_context.hip`): FP32 (and with `_fp64` FP64) steps rotate the previous eigenvectors, the gradient takes W from G_oo after a pseudo step (AI-generated, machine-tested; [docs/ROCM_RESIDENT_LOOP_WP.md](../../../../docs/ROCM_RESIDENT_LOOP_WP.md) WP7).
- ROCm post-SCF dispersion and CN on the device (`-xtb.rocm_d4_postscf_device`, `-xtb.rocm_d3_device`, `-xtb.rocm_cn_device`, default true; `d4PostScf`/`d3Dispersion` in `rocm/xtb_hip_context.hip`): GFN2 D4 needs the resident-loop reference tables of the geometry, else the host code runs (AI-generated, machine-tested; [docs/ROCM_RESIDENT_LOOP_WP.md](../../../../docs/ROCM_RESIDENT_LOOP_WP.md) WP8a).
- ROCm gradient path after the resident loop (`-xtb.rocm_grad_resident`, `-xtb.rocm_grad_tail_device`, `-xtb.rocm_defer_host_integrals`, default true; `convPotentials`/`gradTail` in `rocm/xtb_hip_context.hip`): device gradient potentials, tail and repulsion; with `GpuScfBackend::deferHostCopies` the host copies of S/H0, the multipole integrals and matrices are built only on demand (AI-generated, machine-tested; [docs/ROCM_RESIDENT_LOOP_WP.md](../../../../docs/ROCM_RESIDENT_LOOP_WP.md) WP8b).
- ROCm fixed-point gradient accumulation (`-xtb.rocm_fixed_point`, default true; `xfix_add` in `rocm/xtb_hip_context.hip`): the four gradient kernels and the host-loop moments sum order-independently, an out-of-range value fails the device call; the SCF loop has no FP atomics (AI-generated, machine-tested; [docs/ROCM_RESIDENT_LOOP_WP.md](../../../../docs/ROCM_RESIDENT_LOOP_WP.md) WP10).
- ROCm GFN2 on-the-fly multipole interaction (`-xtb.rocm_otf_multipole`, default true; `potBeginOtf`/`k_mp_potential_otf` in `rocm/xtb_hip_context.hip`): with `-gpu_multipole_otf auto` only above the 4 GB stored limit (`GpuScfBackend::preferStoredMultipole`), `on` forces it; bit-identical to the stored kernels (AI-generated, machine-tested; [docs/ROCM_RESIDENT_LOOP_WP.md](../../../../docs/ROCM_RESIDENT_LOOP_WP.md) WP11).
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
