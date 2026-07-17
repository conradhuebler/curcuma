[![CodeFactor](https://www.codefactor.io/repository/github/conradhuebler/curcuma/badge)](https://www.codefactor.io/repository/github/conradhuebler/curcuma) [![Build](https://github.com/conradhuebler/curcuma/workflows/AutomaticBuild/badge.svg)](https://github.com/conradhuebler/curcuma/actions)  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4302722.svg)](https://doi.org/10.5281/zenodo.4302722)

![curcuma Logo](https://github.com/conradhuebler/curcuma/raw/master/misc/curcuma_II.png)

# Curcuma

A simple Open Source molecular modelling tool.

## Download and requirements
Dependencies are fetched automatically via CMake FetchContent (no manual submodule init required).
- [LBFGSpp](https://github.com/conradhuebler/LBFGSpp) a fork of [yixuan/LBFGSpp](https://github.com/yixuan/LBFGSpp/) provides LBFGS optimiser, the fork allows performing single step optimisation without resetting any calculated optimsation history
- [XTB](https://github.com/grimme-lab/xtb) the official xtb program — optional external backend for `xtb-gfn1`, `xtb-gfn2`, `xtb-gfnff` (USE_XTB build flag). Not required for the native GFN1/GFN2/GFN-FF implementations.
- [tblite](https://github.com/tblite/tblite) optional external backend for `tblite-gfn1`, `tblite-gfn2`, `xtb-gfn1`, `xtb-gfn2` (USE_TBLITE build flag). The canonical `gfn1`/`gfn2` methods use the native curcuma implementation instead.
- [simple-d3](https://github.com/dftd3/simple-dftd3) Dispersion Correction D3
- [cpp-d4](https://github.com/conradhuebler/cpp-d4) Fork of the cpp-d4 repository for Dispersion Correction D4
- [CxxThreadPool](https://github.com/conradhuebler/CxxThreadPool) - C++ Thread Pool for parallel calculation
- [eigen](https://gitlab.com/libeigen/eigen) provides eigen C++ library for linear algebra. Eigen is not downloaded automatically, but will be fetched and updated if the build scripts in the **scripts** subdirectory are used.
- [fmt](https://github.com/fmtlib/fmt) formatted console output
- [plumped](https://github.com/plumed/plumed2) Support for Metadynamics, must be compiled manually and enabled manually (Option USE_Plumed)
- [ulysses](https://gitlab.com/siriius/ulysses) Support for several semiemprical models via ulysses

Additionally, [nlohmann/json](https://github.com/nlohmann/json) is obtained via cmake.

A C++/Eigen implementation of the Munkres Algorithmus (Hungarian Method) based on [the workshop here](https://brc2.com/the-algorithm-workshop/) is included.

### GPU acceleration (optional)

The GPU backends for `gfn1`/`gfn2`/`gfnff` are **off by default** and each needs extra
system dependencies (one backend per build dir: `release_cuda/`, `release_rocm/`,
`release_vulkan/`). The default `release/` build needs none of these.

- **CUDA** (`-gpu cuda`, `-DUSE_CUDA=ON`): NVIDIA CUDA toolkit — `nvcc`,
  cuSOLVER, cuBLAS, cudart (Arch: `cuda`). See [docs/SQM_GPU.md](docs/SQM_GPU.md).
- **ROCm / HIP** (`-gpu rocm`, `-DUSE_ROCM=ON -DCMAKE_PREFIX_PATH=/opt/rocm`):
  `hip-runtime-amd`, `rocm-llvm`, `rocm-device-libs`, `rocminfo`; **`rocblas` + `rocsolver`**
  for the GPU eigensolver (xTB) / EEQ solve (GFN-FF). Set `-DROCM_GPU_ARCH` to your GPU's `gfx`
  (e.g. `gfx1150`; `rocminfo | grep gfx`). `USE_ROCM` enables gfn1/gfn2, `USE_ROCM`
  enables gfnff. See [docs/SQM_ROCM.md](docs/SQM_ROCM.md).
- **Vulkan** (`-gpu vulkan`, `-DUSE_VULKAN=ON`): `vulkan-icd-loader` +
  `vulkan-headers` + an FP64-capable driver (AMD `vulkan-radeon`/RADV — **no ROCm needed**,
  NVIDIA `nvidia-utils`, Intel `vulkan-intel`); `shaderc`/`glslang` only to regenerate the
  (committed) SPIR-V. Needs a device with `shaderFloat64`. See [docs/SQM_VULKAN.md](docs/SQM_VULKAN.md).

## Validation Status Labels

Curcuma contains a mix of production-tested and AI-generated code. The following labels appear throughout this README and the internal `CLAUDE.md` documentation to indicate the confidence level of each feature:

| Label | Meaning | Who sets it |
|-------|---------|-------------|
| 🤖 AI-generated | Code written by AI, not yet reviewed by a human | AI |
| ⚙️ Machine-tested | Passes automated tests (CI, ctest) | AI |
| 👁️ Human-reviewed | Human has read and understood the code | Human only |
| ✅ TESTED | Human has run it on real problems and verified correct behaviour | **Human only** |
| ✅ APPROVED | Human confirms correctness, ready for production use | **Human only** |

**Important**: passing automated tests does not imply physical correctness for all inputs. AI-generated scientific code can produce plausible but wrong results in untested regimes. Features marked only 🤖/⚙️ should be cross-checked against an established reference before use in research.

### UFF, xTB, GFN-FF and Dispersion Correction
Curcuma has an interface to tblite, xtb as well simple-d3 and cpp-d4, enabling semiempirical calculations or combinations of UFF with D3, D4 and H4 (no parameters are adjusted yet). To use one of the methods, please add **-method methodname** to your arguments:

UFF (default)
- uff : Universal Force Field

Native force field (no external dependency required):
- **gfnff** : Native C++ GFN-FF — full energy and gradient, validated against Fortran reference (see status below)
- **gfnff** + `-gpu cuda` : CUDA-accelerated variant; topology cached, charges on CPU, all kernels on GPU
- **xtb-gfnff** : GFN-FF via the xtb Fortran library (USE_GFNFF build flag)

Native GFN methods (no external dependency required, canonical backends since AP3 2026-04-25):
- **gfn1** : Native GFN1-xTB — 10/12 validation molecules at 1e-8 vs tblite
- **gfn2** : Native GFN2-xTB — 11/12 validation molecules at 1e-8 vs tblite (only `complex` open at 6.95e-5)

> Native GFN1/GFN2 are validated against tblite to a 1e-8 Eh target — see [docs/SQM_VALIDATION.md](docs/SQM_VALIDATION.md). For explicit tblite or xtb backends use `tblite-gfn1`/`tblite-gfn2` or `xtb-gfn1`/`xtb-gfn2`.

> **d-shell elements (X-I1, June 2026):** native GFN1/GFN2 now handle d-shell basis functions (S, P, Cl, Si and other main-group d elements), matching tblite to ≤1e-8 Eh; analytic gradients FD-validated. CPU only — on `-gpu` a d-shell system falls back to the CPU integral/SCF path. Transition metals are enabled but not yet validated. See [docs/SQM_DSHELL_WP.md](docs/SQM_DSHELL_WP.md).

> Native GFN1/GFN2 can use multiple cores **within one calculation** of a single large molecule: pass `-threads N` to a `-sp`/`-opt`/MD run (default is serial and bit-identical). Integral setup, gradient and Fock build scale ~3–5×; see [docs/SQM_THREADING.md](docs/SQM_THREADING.md).

> Opt-in **MKL-free / GPU-portable eigensolve kernels** are available for the native GFN SCF (MKL stays the default): `-eigensolver native` (own Householder + Cuppen divide-and-conquer), `-eigensolver purify` (0 K density-matrix purification, GEMM-only, no diagonalization), `-eigensolver lobpcg` (seeded block LOBPCG, experimental), and `CURCUMA_EIG_TRED2=blocked` (BLAS-3 blocked tridiagonalization). See [docs/SQM_EIGENSOLVE_GPU.md](docs/SQM_EIGENSOLVE_GPU.md).

> Opt-in **CUDA GPU path** for the native GFN1/GFN2 solver: `-method gfn1|gfn2 -gpu cuda` (build `release_cuda/` with `-DUSE_CUDA=ON`). Staged cuSOLVER/cuBLAS port (the CPU path is unchanged and `#ifdef`-free); both **GFN1** and **GFN2** run a device-resident SCF under the default Broyden mixing, and **Stage 3 builds the integrals (CN/S/H0/L/γ/multipole) on the device and Stage 4 the nuclear gradient — so `-opt`/`-md` are fully device-resident** (only xyz up, gradient+energy down per step; every device kernel matches the CPU elementwise to ~1e-15). 🤖 AI-generated / ⚙️ machine-tested only. See [docs/SQM_GPU.md](docs/SQM_GPU.md).

> **AMD/ROCm** (`-gpu rocm`, build `release_rocm/` with `-DUSE_ROCM=ON`) and **Vulkan compute** (`-gpu vulkan`, hand-written SPIR-V, `-DUSE_VULKAN=ON`) backends for the same `gfn1`/`gfn2`/`gfnff` methods. `-gpu auto` picks the first compiled backend (cuda > rocm > vulkan), else CPU. 🤖 **Vulkan: GFN1 = Stage 2** (device-resident SCF — Fock/eigensolve/density/populations/band on the GPU via a device-built Löwdin S⁻¹ᐟ²; only `v_ao`/`occ` up and `eps`/`pop`/`band` down per iteration), **GFN2 = Stage 1** (per-iteration eigensolve on GPU). gfn1/gfn2 single-point + opt match the CPU bit-for-bit on the validation set (AMD 890M/RADV); integrals/gradient still CPU. **ROCm: GFN1 = Stage 4 (fully device-resident)** — the integral build (CN/S/H0/L/γ), the SCF (Fock/density/eigensolve via HIP kernels + rocBLAS + rocSOLVER) and the nuclear gradient (repulsion/Pulay/Coulomb HIP kernels) all run on the GPU; only the dispersion gradient + CN chain-rule on the host. **GFN2** uses the device integrals + rocSOLVER eigensolver (gradient on host). gfn1/gfn2 single-point + opt match the CPU bit-for-bit, incl. the full `-opt` trajectory (AMD 890M, needs `rocsolver`+`rocblas`). **ROCm GFN-FF** (`-DUSE_ROCM=ON`, June 2026): the full energy + nuclear-gradient kernel stack runs on the GPU (single-TU hipify of the CUDA gfnff kernels; EEQ via rocSOLVER `dpotrf`/`dgetrf` + host CPU-Schur); single-point energy and gradient match CPU ≤1e-7 on water/CH4/caffeine/231-atom complex. **Two opt-in CUDA-only GFN-FF GPU flags (default OFF, ROCm mirrors pending):** `-gfnff.eeq_mixed_precision` (FP32-factor + FP64-refine EEQ solve) and `-gfnff.gpu_disp_pairs_on_device` (on-device D4 pair build) — bit-identical to the host but not a measured speedup (residency milestones). See [docs/SQM_ROCM.md](docs/SQM_ROCM.md) / [docs/SQM_VULKAN.md](docs/SQM_VULKAN.md) / [docs/GFNFF_PERFORMANCE_LEVERS.md](docs/GFNFF_PERFORMANCE_LEVERS.md).

> Opt-in **approximate large-system modes** scale the native GFN SCF beyond ~1000 atoms by exploiting locality (default is the exact dense path): `-large_system_mode fragments` (disconnected-fragment SCF, energy+gradient, `-eigensolver` propagates per fragment), `-large_system_mode dc` (divide-and-conquer, energy-only, `-eigensolver` propagates per sub-block, `-large_system_buffer_bohr` accuracy knob), `-large_system_mode sparse` (non-orthogonal density purification, 0 K gapped, `-eigensolver` ignored, `-large_system_sparse_threshold` knob). Each converges to the dense energy as its knob tightens; combining `-large_system_mode=fragments|dc` with `-eigensolver=purify` requires `-electronic_temperature 0` (hard error otherwise). See [docs/SQM_LARGE_SYSTEMS.md](docs/SQM_LARGE_SYSTEMS.md).

> Opt-in **multi-step SCC extrapolation** for the native GFN SCF cuts SCF iterations across geometry steps in `-opt`/`-md` by predicting the next charge state from several past converged steps (generalises the 1-step warm-start; default `none` is unchanged). `-scf_extrapolation aspc` (Kolafa ASPC, best for fixed-timestep MD) or `-scf_extrapolation gauss` (least-squares, better for irregular opt steps), with `-scf_extrapolation_order`. The safe default `guess` coupling still converges the SCF fully; `-scf_extrapolation_apply xlbomd` is an experimental extended-Lagrangian Born-Oppenheimer mode (time-reversible auxiliary density + converged corrector, for low MD energy drift). On a smooth caffeine trajectory, `aspc`/`gauss` roughly halve SCF iterations (gfn2 215→90, gfn1 170→79) with bit-identical converged energy. 🤖 AI-generated / ⚙️ machine-tested only. See [docs/SQM_SCF_EXTRAPOLATION.md](docs/SQM_SCF_EXTRAPOLATION.md).

xtb methods:
- xtb-gfnff : GFN-FF via the xtb library
- xtb-gfn1
- xtb-gfn2

Using only **d3** or **d4** should be possible.

Native GFN2 includes an analytic D4 dispersion charge-response gradient
(∂E_D4/∂q · ∂q/∂x). The zeta charges default to a single-shot dftd4 EEQ model
(`-d4_charge_source eeq`, analytic ∂q/∂x); `-d4_charge_source mulliken` feeds the
GFN2 SCF charges (energy + ∂E/∂q; the CPSCF gradient response is still pending —
see [docs/D4_Q_RESPONSE.md](docs/D4_Q_RESPONSE.md)). Current alignment vs tblite:
11/12 at 1e-8 — only `complex` (231 atoms, 6.95e-5 Eh residual) remains open.
Status tracked in [docs/GFN2_NATIVE_ROADMAP.md](docs/GFN2_NATIVE_ROADMAP.md) and
[docs/GFN2_D4_STATUS.md](docs/GFN2_D4_STATUS.md); `ctest -L d4_diag`.

The native GFN SCF defaults to `broyden` mixing — a modified-Broyden quasi-Newton
scheme on the SCC charge vector, the same mixer tblite/xtb use — which converges
large polar systems that the old Fock-DIIS diverged on (e.g. the 231-atom
`complex` now converges from the bare guess with plain `-method gfn2`). Other
modes remain selectable: `-scf_mode diis|plain|level-shift` and `-scf_guess
h0|eeq` (plus `-scf_damping`, `-diis_start`, `-level_shift`). See
[docs/SCF_MODES.md](docs/SCF_MODES.md).

Please cite xtb, tblite etc if external methods are used within curcuma! The most recent information can be found at the respective github pages, some are listed below.

UFF
- J. Am. Chem. Soc. (1992) 114(25) p. 10024-10035,
- with the H4 hydrogen bond correction (J. Chem. Theory Comput. 8, 141-151 (2012)) included (same parameters as applied in case of PM6-D3 for now).

GFN-FF (native C++ implementation):
- S. Spicher and S. Grimme, Angew. Chem. Int. Ed. 2020, 59, 15665. DOI: 10.1002/anie.202004239

### Native GFN-FF Status (April 2026)

The native `gfnff` implementation is **AI-implemented and machine-tested** — human production testing is pending.

**What works (validated by automated tests):**
- All energy terms: bonds, angles, torsions, inversions, repulsion, dispersion (D4), Coulomb (EEQ), hydrogen bonds, halogen bonds, triple-bond torsions, BATM, ATM
- Analytical gradients for all terms; GPU (CUDA) analytical gradients correct
- 20 validation molecules (H₂ to a 1280-atom polymer) — energy vs. Fortran reference within tolerances
- CUDA acceleration: topology caching, async CPU/GPU overlap, shared-memory reduction
- Geometry optimization and MD using gradients
- **GFN-FF ALPB solvation**: self-consistent Born reaction field coupled into EEQ (`A_eeq += B`); validated against xtb 6.7.1 (`--gfnff --alpb`) to ≤1e-8 Eh (7 molecules × 4 solvents, `ctest -L gfnff_solvation`; June 2026). Gradient FD-validated at frozen solvated charges (same approximation as the Fortran reference). `-gfnff.solvent_model gbsa` maps to ALPB (GFN-FF has no separate GBSA model; warns at runtime). See [docs/SQM_SOLVATION_WP.md](docs/SQM_SOLVATION_WP.md).

**Not validated / not implemented:**
- **Periodic boundary conditions**: Not implemented
- **Organometallics / transition metals**: No test molecule with metal center; parameter quality unknown
- **Gradient accuracy for large systems**: Dispersion gradients show √N accumulation error (expected for O(N²) terms, scientifically acceptable for MD/opt)
- **GPU energy for polymer (1280 atoms)**: 8.9 µEh vs. 1 µEh tolerance — pre-existing, under investigation

**Known differences from Fortran reference** (see [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md)):
- Sub-mEh agreement for most small/medium molecules
- EEQ charge environment corrections (dxi) partially implemented
- Metal-specific EEQ corrections (fqq) not implemented

Do not use for production on untested system classes without cross-checking against `xtb-gfnff`.

D3:
- J. Chem. Phys. 132, 154104 (2010); https://doi.org/10.1063/1.3382344

D4:
- E. Caldeweyher, C. Bannwarth and S. Grimme, J. Chem. Phys., 2017, 147, 034112. DOI: 10.1063/1.4993215
- E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth and S. Grimme, J. Chem. Phys., 2019, 150, 154122. DOI: 10.1063/1.5090222

Dispersion correction parameters are yet complicated to change, this will be improved sooner than later.

## Compiling
To compile Curcuma you will need [CMake](https://cmake.org/download/) 3.15 or newer and a C++17-capable compiler, both gcc and icc (quite recent version) work. One possible option is MinGW. For Windows, it is further necessary to add the bin-folder in the MinGW installation to the path (Edit the system environment variables > Environment Variables > under "System Variables" select "Path" > Edit > New > paste path, for example "C:\MinGW\bin").

To obtain the most recent version
```sh
git clone --recursive https://github.com/conradhuebler/curcuma
```
For Windows: you need to make sure to navigate to a folder outside the Windows System before clone to avoid conflicts of usage rights.

Compile it as follows on Unix Platform:
```sh
cd curcuma 
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```
For Windows: you need to tell the CMD which compiler to use to avoid errors if using the wrong compiler. Using MinGW as example, the "cmake ..." command would become
```sh
cmake .. -DCMAKE_BUILD_TYPE=Release -G "MinGW Makefiles"
```

### Windows with OpenMP support (w64devkit)

Standard MinGW does not ship with OpenMP. To enable OpenMP on Windows, use [w64devkit](https://github.com/skeeto/w64devkit), which provides a GCC toolchain with OpenMP support out of the box.

1. Install [MinGW](https://sourceforge.net/projects/mingw/) and [CMake](https://cmake.org/download/), and add their `bin` directories to the System environment variables (`PATH`).
2. Download and extract [w64devkit](https://github.com/skeeto/w64devkit/releases).
3. Open the w64devkit terminal from the extracted folder (run `w64devkit.exe`).
4. Inside that terminal, build Curcuma:

```sh
git clone --recursive https://github.com/conradhuebler/curcuma
cd curcuma
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -G "MinGW Makefiles"
mingw32-make
```

The w64devkit environment provides the correct `libgomp` runtime, so `-DUSE_OpenMP=ON` will be picked up automatically by CMake.

## Modern Parameter System (October 2025)

Curcuma features an **automated parameter registry system** for all molecular modeling capabilities:

- **Auto-generated help** directly from source code annotations
- **Type-safe** parameter definitions with compile-time validation
- **JSON export/import** for reproducible computational workflows
- **Alias support** for multiple parameter names
- **Build-time validation** detects duplicate or conflicting parameters

### For Users

**Export default configuration:**
```sh
./curcuma -export-config analysis > my_analysis.json
```

**Modify and run with custom config:**
```sh
# Edit my_analysis.json with your preferred settings
./curcuma -analysis input.xyz -import-config my_analysis.json
```

**List available modules:**
```sh
./curcuma -list-modules
```

**Capture and replay a full run (2026):** any registered parameter is reachable by its flat CLI name, and `-export_run` / `-import_config` form a round-trip. See [docs/CLI_ROUND_TRIP.md](docs/CLI_ROUND_TRIP.md).
```sh
./curcuma -sp water.xyz -method gfnff -cn_cutoff_bohr 5.5 -export_run run.json
./curcuma -import_config run.json                       # replay
./curcuma -import_config run.json -cn_cutoff_bohr 7.0   # replay with override
```

### For Developers

All new capabilities must use the Parameter Registry System. See:
- **Technical Documentation**: [docs/PARAMETER_SYSTEM.md](docs/PARAMETER_SYSTEM.md)
- **Migration Guide**: [docs/PARAMETER_MIGRATION_GUIDE.md](docs/PARAMETER_MIGRATION_GUIDE.md)
- **Reference Implementation**: `src/capabilities/analysis.h`

**Quick Start for New Capabilities:**
```cpp
// In your capability header:
#include "src/core/parameter_macros.h"

class MyCapability : public CurcumaMethod {
private:
    BEGIN_PARAMETER_DEFINITION(my_capability)

    PARAM(max_iterations, Int, 100,
          "Maximum number of iterations",
          "Algorithm",
          {"max_iter"})

    PARAM(output_file, String, "",
          "Optional output file path",
          "Output",
          {"out"})

    END_PARAMETER_DEFINITION
};
```

Build system automatically extracts parameters and generates unified registry with validation and help.

# Usage

## General
curcuma catches the Ctrl-C signals from the console if used on Linux platform. It will then create an empty file called "stop". Some methods, like **confscan**, regularly check for that file and finalise the current task. If Ctrl-C is signaled and a "stop" file already exists, curcuma will stop immediately.

## RMSD Calculator 
```sh
curcuma -rmsd file1.xyz file2.xyz
```
Computes RMSD. If the two structures are ordered differently, curcuma will automatically reorder the atom list. To force reordering use
```sh
curcuma -rmsd file1.xyz file2.xyz -reorder
```

Two basic approaches are currently implemented, one incremental (testing many, but not all possible orders, parallelised) method without any requirements and another Kuhn-Munkres like approach. 
The Kuhn-Munkres like approach needs a prior alignment (only cude) of the structures, which may be obtained using a smaller substructure (template). One way is to use a molecule in a supramolecular structure . Use
```sh
-method template -fragment 1
```
to use the second structure, eg the one bound non-covalently by the first structure, as template. With this approach a much faster reordering is obtained. Omitting -fragment, curcuma tries using the smallest fragment.

Alternatively, generate templates using the incremental approach. For this, all atoms of one element (nitrogen is the default) are used. After the templates are generated, the structures oriented according the templates are reordered using the Kuhn-Munkres like approach. Use
```sh
-method hybrid -element 8
```
for taking oxygen as template.
 
Add
```sh
-heavy
```
to perform calculation only on non-proton atoms.

A new, and all previous published (and not yet published) methods was propsed Feb 2023 by Vásquez-Pérez and coworkers. It outperforms all methods so far (and sometimes the methods implemented in curcuma). 
It can be obtained at [Github](https://github.com/qcuaeh/molalignlib) and included in curcuma for RMSD calculation conformational filtering.
```sh
-reorder -method molalign -molalignbin /anypath/molalign
```

If the method was used, please cite the authors!
[J. Chem. Inf. Model. 2023, 63, 4, 1157–1165](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.2c01187)
The method can be applied during confscan, however, problems with the random numbers occur during the test runs for larger problems, making the less ideal than the natively in curcuma implemented methods.
```sh
-rmsdmethod molalign -molalignbin /anypath/molalign
```

```json
{ "reorder", false },
{ "check", false },
{ "heavy", false },
{ "fragment", -1 },
{ "fragment_reference", -1 },
{ "fragment_target", -1 },
{ "init", -1 },
{ "pt", 0 },
{ "silent", false },
{ "storage", 1.0 },
{ "method", "incr" },
{ "noreorder", false },
{ "threads", 1 },
{ "Element", 7 },
{ "DynamicCenter", false },
{ "order", "" },
{ "check", false },
{ "topo", 0 },
{ "write", 0 },
{ "moi", false },
{ "update-rotation", false },
{ "damping", 0.8 },
{ "split", false },
{ "nomunkres", false },
{ "dmix", -1 },
{ "molalignbin", "molalign" }
```


## Docking tool
Some docking can be performed (WIP).

Use
```sh
curcuma -dock -host A.xyz -guest B.xyz
```
to perform docking of B as guest and A as host molecule, use
```sh
curcuma -dock -complex AB.xyz
```
to perform docking on a complex or use
```sh
curcuma -dock -complex AB.xyz -guest C.xyz
```
to replace substrat structures in complex AB with C.


Use
```sh
curcuma -dock -host A.xyz -guest B.xyz -Step_X X  -Step_Y Y -Step_Z Z
```
with X, Y and Z being the steps of rotation. With X = 10, 10 rotations with 360/X ° will be performed.

Use
```sh
curcuma -dock -host A.xyz -guest B.xyz -Pos_X X  -Pos_X Y -Pos_X Z
```
with {X, Y, Z} being the initial anchor position for the substrat.

After docking a PseudoFF optimisation of the docking position will be performed, where the Lennard-Jones-Potential between both structures is calculated. XTB is then used to preoptimise the unique docking structures and the results are filtered using ConfScan and the template based reordering approach.

```json
{ "Pos_X", 0.0 },
{ "Pos_Y", 0.0 },
{ "Pos_Z", 0.0 },
{ "AutoPos", true },
{ "Filter", true },
{ "PostOpt", true },
{ "Step_X", 10 },
{ "Step_Y", 10 },
{ "Step_z", 10 },
{ "Host", "none" },
{ "Guest", "none" },
{ "Complex", "none" },
{ "scaling", 1.5 },
{ "NoOpt", false },
{ "CentroidMaxDistance", 1e5 },
{ "CentroidTolDis", 1e-1 },
{ "RotationTolDis", 1e-1 },
{ "Threads", 1 },
{ "DockingThreads", 1 },
{ "Charge", 0 },
{ "Cycles", 1 },
{ "RMSDMethod", "incr" },
{ "RMSDThreads", 1 },
{ "RMSDElement", 7 }
```



## Conformation Filter
Curcuma has some conformation filter based on energy, rmsd, rotation constants and rank limitation. As some structures may be identic, yet can not be aligned due to atom ordering, depending on the difference of the energy and rotational constants and the rmsd, automatic reordering in rmsd calculation will be performed.

Use
```sh
curcuma -confscan conformation.xyz
```
to simple filter the conformation. The results will be stored in an additional xyz file.

Use
```sh
curcuma -confscan conformation.xyz -MaxHTopoDiff 0
```
to ensure, that even the rmsd is smaller than the threshold, the second molecule is only rejected if there is no difference in the hydrogen bond pattern. Set to ***-MaxHTopoDiff 2*** if two two changes are allowed. A detailed description will follow some time.
```sh
curcuma -confscan conformation.xyz -reorder
```
to force reordering for every rmsd calculation.

Add
```sh
-heavy
```
to perform rmsd calculation and reordering only on non-proton atoms. Reordering is much faster then!

Adding
```sh
-RMSDmethod template
```
the template based reordering is used, **hybrid** is also working.

With
```sh
-Useorders X
```
up to X reorder results will be reused from the last reorder calculation. Template based approaches result in only one reorder rule, however X is set to 0 in automatically, but can be changed with that argument.

By adding the argument
```sh
-accepted accepted.xyz
```
a file with already accepted structures can be passed to curcuma. Molecules in that file (accepted.xyz) will be rejected if they appear in the conformation.xyz, thus several files with conformation can be joined.

Confscan will write a restart file, finalise and quit, if a file called "stop" is found in the working directory. Such file will be generated if Ctrl-C is hit or if it is created using for example the **touch stop** command.
Within a restart file, the last energy difference and the atom indicies from reordering are stored. A restart file will automatically be read upon the start of curcuma. The content of the restart file will be used to speed up the 2nd step of the conformation filtering procedure.

Confscan supports the molalign tool. However, as too often reordering with molalign is not working, it can efficiently be used if the RMSD is only slightly above the threshold. 
```sh
-domolalign 1.1
```
Sets the threshold to 1.1*RMSDthreshold. If the molecule was accepted as to different, but the RMSD is blow 1.1*RMSDthreshold molalign will check too.

Confscan write a statistic file, where for each rejected molecule the reference alongside the energy difference and the RMSD is printed out. Furthermore, the reordered indices are given, if available. Molalign does not return the reordered indices, hence they are empty or marked **0,0** if the reordered was finally performed using molalign in a standard run.

```json
{ "noname", true },
{ "restart", true },
{ "heavy", false },
{ "rmsd", -1 },
{ "rank", -1 },
{ "writeXYZ", false },
{ "forceReorder", false },
{ "check", false },
{ "energy", 1.0 },
{ "maxenergy", -1.0 },
{ "preventreorder", false },
{ "scaleLoose", 1.5 },
{ "scaleTight", 0.1 },
{ "scaleLooseEnergy", 1.2 },
{ "scaleTightEnergy", 0.1 },
{ "scaleLooseRotational", 1.2 },
{ "scaleTightRotational", 0.1 },
{ "scaleLooseRipser", 1.2 },
{ "scaleTightRipser", 0.1 },
{ "skip", 0 },
{ "allxyz", false },
{ "update", false },
{ "MaxParam", -1 },
{ "UseOrders", -1 },
{ "RMSDMethod", "hybrid" },
{ "MaxHTopoDiff", -1 },
{ "threads", 1 },
{ "RMSDElement", 7 },
{ "accepted", "" },
{ "method", "" },
{ "lastdE", -1 },
{ "fewerFile", false },
{ "dothird", true },
{ "skipfirst", false },
{ "ignoreRotation", false },
{ "ignoreBarCode", false },
{ "skipless", false },
{ "looseThresh", 7 },
{ "tightThresh", 3 },
{ "update-rotation", false },
{ "damping", 0.8 },
{ "split", false },
{ "writefiles", false },
{ "nomunkres", false },
{ "molalignbin", "molalign" },
{ "ripser_xmax", 4 },
{ "ripser_xmin", 0 },
{ "ripser_ymax", 4 },
{ "ripser_ymin", 0 },
{ "ripser_bins", 10 },
{ "ripser_scaling", 0.1 },
{ "ripser_stdx", 10 },
{ "ripser_stdy", 10 },
{ "ripser_ratio", 1 },
{ "ripser_dimension", 2 },
{ "domolalign", -1 }
```

```cpp
/* rotational = 1
 * ripser     = 2
 * energy     = 4 */
int looseThresh = 1 * (diff_rot < m_diff_rot_threshold_loose) + 2 * (diff < m_diff_ripser_threshold_loose) + 4 * (std::abs(mol1->Energy() - mol2->Energy()) * 2625.5 < m_diff_energy_threshold_loose);
if ((looseThresh & m_looseThresh) == m_looseThresh) 
{

}
```
## Find unique structures in trajectories
xyz and trj are handled equally.
```sh
curcuma -rmsdtraj XXX.trj -writeUnique -rmsd 1.5
```


```json
{ "writeUnique", false },
{ "writeAligned", false },
{ "rmsd", 1.5 },
{ "fragment", -1 },
{ "reference", "none" },
{ "second", "none" },
{ "heavy", false },
{ "pcafile", false },
{ "allxyz", false },
{ "RefFirst", false },
{ "noreorder", true },
{ "opt", false },
{ "filter", false },
{ "writeRMSD", true },
{ "offset", 0 }
```

## Geometry optimisation (batch mode possible)
Geometry optimisation can be performed with curcuma using 
```sh
curcuma -opt XXX.xyz
```
A file called XXX.opt.xyz with the optimised structures will be written. The individual steps are stored in XXX.trj.xyz. The number of threads can be controlled with
```sh
-threads X
```

```json
{ "writeXYZ", true },
{ "printOutput", true },
{ "dE", 0.1 },
{ "dRMSD", 0.01 },
{ "method", "uff" },
{ "MaxIter", 5000 },
{ "LBFGS_eps", 1e-5 },
{ "StoreIntermediate", 2000 },
{ "SingleStep", 20 },
{ "ConvCount", 11 },
{ "GradNorm", 0.001 },
{ "Threads", 1 },
{ "Charge", 0 },
{ "Spin", 0 },
{ "SinglePoint", false },
{ "optH", false },
{ "serial", false }
```


```cpp
/*
 * Energy = 1
 * RMSD = 2
 * LBFGS Conv = 4
 * Gradient Norm = 8
 * */
converged = 1 * (abs(fun.m_energy - final_energy) * 2625.5 < dE)
    + 2 * (driver->RMSD() < dRMSD)
    + 4 * (solver.isConverged())
    + 8 * (solver.final_grad_norm() < GradNorm);
perform_optimisation = (converged != ConvCount) && (fun.isError() == 0);
}

## Reorder and Align trajectories
To reorder trajectory files with dissordered atomic indicies, for example after merging several minimum energy path files from NEB calculation, use
```sh
curcuma -rmsdtraj XXX.xyz -writeAligned
```
Reordering will be done with respect to the previouse structure in the trajectory. If the first structure should be used, add ***-reffirst*** as additional argument. The new trajectory is called XXX_aligned.xyz.

Using ***-rmsdtraj*** argument, a file **XXX_rmsd.dat** will be written, where the rmsd is stored.

## Distance and angle calculation

```sh
curcuma -distance XXX.trj atom1 atom2
```

```sh
curcuma -angle XXX.trj atom1 atom2 atom3
```

The index starts with 1. Using grep and sed via ***|grep '::' |sed 's/:://g'*** omitts unused output.

## Compare two RDG vs sign(λ<sub>2</sub>)ρ plots 
Using 
```sh
curcuma -nci file1.dat file2.dat
```
one can ''remove'' RDG vs sign(λ<sub>2</sub>)ρ points which occur in both plots (file1.dat and file2.dat). The similarity of two points is set to true, if the distance is below a threshold distance, which is defined by the averaged distance of two adjacent points.

## Molecular Dynamics and Metadynamics
Curcuma has now a Molecular Dynamics modul, which can be used with:
```sh
curcuma -md input.xyz
```

### Possible options
```json
{ "writeXYZ", true },
{ "printOutput", true },
{ "MaxTime", 5000 },
{ "T", 298.15 },
{ "dt", 1 }, // single step in fs
{ "rm_COM", 100 }, // remove translation and rotation every x fs
{ "charge", 0 },
{ "Spin", 0 },
{ "rmrottrans", 0 },
{ "nocenter", false },
{ "dump", 50 },
{ "print", 1000 },
{ "unique", false },
{ "rmsd", 1.5 },
{ "opt", false },
{ "hmass", 1 },
{ "velo", 1 },
{ "rescue", false },
{ "coupling", 10 },
{ "MaxTopoDiff", 15 },
{ "impuls", 0 },
{ "method", "uff" },
{ "impuls_scaling", 0.75 },
{ "writeinit", false },
{ "initfile", "none" },
{ "norestart", false },
{ "writerestart", 1000 },
{ "rattle", false },
{ "rattle_tolerance", 1e-6 },
{ "rattle_maxiter", 10 },
{ "thermostat", "csvr" },
{ "respa", 1 },
{ "dipole", false },
{ "seed", 1 },
{ "cleanenergy", false },
{ "wall", "none" }, // can be spheric or rect
{ "wall_type", "logfermi" }, // can be logfermi or harmonic
{ "wall_spheric_radius", 0 },
{ "wall_xl", 0 },
{ "wall_yl", 0 },
{ "wall_zl", 0 },
{ "wall_x_min", 0 },
{ "wall_x_max", 0 },
{ "wall_y_min", 0 },
{ "wall_y_max", 0 },
{ "wall_z_min", 0 },
{ "wall_z_max", 0 },
{ "wall_temp", 298.15 },
{ "wall_beta", 6 },
{ "mtd", false },
{ "plumed", "plumed.dat" }
```

For example, using 
```sh
curcuma -md input.xyz -method gfnff
``` 
the GFN-FF approach will be used.

```sh
curcuma -md input.xyz -method gfnff -T 500  -berendson 200 -dt 1 -hmass 1 -thermostat_steps 400 -velo 4 -maxtime 2e4 -dt 0.5 -impuls 500 -impuls_scaling 0.75
``` 
will perform some kind of conformational search using GFN-FF. Results are stored in **input.unique.xyz**! Repeating it will result in other conformations and the previous results stored in **input.unique.xyz** will be overwritten. Bonds may break from time to time ...

Rattle can be used to constrain (currently) all bonds, allowing larger time steps for integration. Up to 8 fs might be possible.
```sh
curcuma -md input.xyz -rattle -dt 4
``` 

The MD implementation integrates well into curcuma, hence calculation can be stopped with Ctrl-C (or a "stop" file) and will be resumed (velocities and geometries are stored) if a restart file is found.

The thermostat target temperature can follow a multi-stage **ramp** and individual atom subsets can be thermostatted as separate **regions**:
```sh
curcuma -md input.xyz -method gfnff -temperature 300 -temp_ramp true -temp_schedule "600:steps:5000;300:reach:10"
```
See [docs/TEMPERATURE_RAMP.md](docs/TEMPERATURE_RAMP.md) for the schedule grammar (`steps`/`reach`), the `temp_regions` JSON array, live temperature control, and per-thermostat support.

With
```sh
curcuma -md input.xyz -mtd
``` 
a metadynamics simulation can be performed using plumed. It is a ***plumed.dat*** expected, or can be set with
```sh
curcuma -md input.xyz -mtd -plumed plumed.dat
```

See [docs/PLUMED_HELP.md](docs/PLUMED_HELP.md) for the full PLUMED integration guide (unit conversions, output files, available CVs, thermal equilibration gate, internal RMSD-MTD).

## Conformational Search (dual-method)

The MD-driven conformational search (`-confsearch`) can explore with a cheap method and refine/rank the discovered conformers with a more accurate one:

```sh
curcuma -confsearch input.xyz -md_method gfnff -opt_method gfn2
```

`-md_method` runs the MD exploration and the pre-optimization; `-opt_method` runs the per-cycle accurate re-optimization and the final ranking. Both fall back to `-method` when unset, so `curcuma -confsearch input.xyz -method gfnff` keeps the single-method behaviour. See [docs/CONFSEARCH_DUAL_METHOD.md](docs/CONFSEARCH_DUAL_METHOD.md).

The search is **restartable** with `-restart`: a self-contained checkpoint (bias pool, cumulative conformers, seeds, energies, schedule position) is written after every MD phase and every temperature cycle, into the BMT dir and copied back to the start directory. Re-running the same command with `-restart` resumes from it (kill the process to interrupt; the checkpoint persists). See [docs/CONFSEARCH_RESTART.md](docs/CONFSEARCH_RESTART.md).

## Output Directory System (BMT)

By default, all curcuma commands create a **Basename.Method.Timestamp** directory for their output files. For example:

```sh
curcuma -md water.xyz -method gfnff
# Output goes to: water.md.20260609_143052/
```

The BMT directory contains all trajectory files, restart data, and a `metadata.json` file with calculation details (basename, method, timestamp, input file). This keeps the working directory clean and makes it easy to compare runs.

To copy specific files back to the working directory after the calculation finishes, use the `-bak` flag:

```sh
curcuma -opt water.xyz -method gfnff -bak water.opt.xyz
# water.opt.xyz is copied from the BMT directory to the working directory
```

Multiple files can be specified: `-bak water.opt.xyz -bak water.trj.xyz`.

To disable BMT and write output to the working directory (legacy behavior):

```sh
curcuma -md input.xyz -method uff -no_bmt
```

# Funding
The development of curcuma is funded by:

-   2026 Stiftung Innovation in der Hochschullehre

![STIL Logo](https://github.com/conradhuebler/curcuma/raw/master/STIL_Funding.jpg)


# Citation
Please cite the software package if you obtain results:
[conradhuebler/curcuma: Curcuma Zenodo Citation](https://doi.org/10.5281/zenodo.4302722)

Have a lot of fun!
