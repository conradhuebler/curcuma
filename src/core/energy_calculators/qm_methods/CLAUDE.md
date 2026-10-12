# CLAUDE.md - Quantum Methods Directory

Quantum-mechanical and semi-empirical method implementations and the `ComputationalMethod` wrappers around them.
The former long version (status log of 2025 and 2026, per-feature performance records, old verbosity text) is kept in
[docs/archive/QM_METHODS_NOTES_2026-10.md](../../../../docs/archive/QM_METHODS_NOTES_2026-10.md), with a list of the statements in it that are out of date.
Detailed architecture: `QM_ARCHITECTURE.md` in this directory and [docs/SQM_VALIDATION.md](../../../../docs/SQM_VALIDATION.md).

> ⚠️ **All native methods here are 🤖 AI-implemented and ⚙️ machine-tested only. None is ✅ TESTED.**

## Layout

- **Wrappers** (`ComputationalMethod`): `native_xtb_method.*` (gfn1/gfn2), `xtb_method.*` (external XTB), `xtb_gpu_method.*`, `xtb_hip_method.*`, `xtb_vulkan_method.*`, `gfnff_method.*` (adapter to the native GFN-FF in `../ff_methods/`), `eht_method.*`, `nddo_method.*`, `orca_method.*`, `dispersion_method.*`. External: `tbliteinterface.*`, `xtbinterface.*`, `ulyssesinterface.*`, `orcainterface.*`, `dftd3interface.*`, `dftd4interface.*`.
- **Native xTB** (`xtb_native.cpp/h` and the `xtb_*.cpp` group: `xtb_h0`, `xtb_scf`, `xtb_gradient`, `xtb_coulomb`, `xtb_multipole`, `xtb_thirdorder`, `xtb_response`, `xtb_dc`, `xtb_fragment_scf`, `xtb_sparse`): `curcuma::xtb::XTB` implements `QMInterface` directly and does not derive from `QMDriver`. The only GFN1/GFN2 parameter tables are `parameters/gfn1_params.hpp` and `parameters/gfn2_params.hpp`.
- **Driver** (`qm_driver.cpp/h`): base class for the STO/GTO-basis methods EHT and NDDO (matrix storage, `MakeOverlap/MakeH` hooks). NDDO (`nddo.*`, `nddo_params.*`) covers MNDO/AM1/PM3/PM6.
- **Integrals and solvers**: `STOIntegrals.hpp`, `GTOIntegrals.hpp`, `STO_CGTO.hpp`, `xtb_ao_utils.hpp`, `xtb_multipole_ints.hpp`, `native_eigensolver.*`, `ParallelEigenSolver.hpp`, `basissetparser.hpp`, GPU contexts under `cuda/` and the plugin loader `../gpu_plugin.*`.
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
- The tblite reference is the target for gfn1/gfn2; the two documented places where xtb differs are opt-in switches (`-xtb.sto6g_legacy_4sp`, `-xtb.d4_atm_cutoff`) (#27, #30).
- Element-indexed tables: transition-metal shells are ordered `[d,s,p]`, so parameters indexed by angular momentum must not be read by shell index (#5).
- Silent mode (verbosity 0) must produce no output except errors; verbosity levels and logger functions: [docs/LOGGING_SYSTEM.md](../../../../docs/LOGGING_SYSTEM.md). No numerical path may depend on the print level.
- GPU backends are runtime plugins; fallbacks to the CPU are counted and reported (`gpu_fallback.h`).

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
