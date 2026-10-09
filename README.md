[![CodeFactor](https://www.codefactor.io/repository/github/conradhuebler/curcuma/badge)](https://www.codefactor.io/repository/github/conradhuebler/curcuma) [![Build](https://github.com/conradhuebler/curcuma/workflows/AutomaticBuild/badge.svg)](https://github.com/conradhuebler/curcuma/actions)  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4302722.svg)](https://doi.org/10.5281/zenodo.4302722)

![curcuma Logo](https://github.com/conradhuebler/curcuma/raw/master/misc/curcuma_II.png)

# Curcuma

An open source toolkit for molecular modelling: energies, gradients, geometry optimisation, molecular
dynamics, conformer search and structure analysis with native implementations of GFN-FF, GFN1-xTB,
GFN2-xTB, PM3/AM1/MNDO and UFF. The native methods need no external program.

## What to know before you install

**Validation status.** A large part of curcuma, including every native quantum-chemical and force-field
method, was written with AI assistance and checked by automated tests and by comparison with reference
programs. It has **not been tested by a human on production problems**. Passing tests does not imply
physical correctness for every input; cross-check results against an established reference before using
them in research. The labels used in this README and in the developer documentation:

| Label | Meaning | Who sets it |
|-------|---------|-------------|
| 🤖 AI-generated | Code written by AI, not reviewed by a human | AI |
| ⚙️ Machine-tested | Passes automated tests (CI, ctest) | AI |
| 👁️ Human-reviewed | A human has read and understood the code | Human only |
| ✅ TESTED | A human has run it on real problems and it behaves correctly | **Human only** |
| ✅ APPROVED | Human confirms correctness, ready for production | **Human only** |

**Methods** (`-method NAME`; `curcuma -methods` lists what your build provides). Native methods are 🤖 ⚙️:

| Method | Kind | Use | What it was compared with |
|---|---|---|---|
| `gfnff` | GFN-FF force field, native | **default of every capability**; fast (milliseconds per gradient) | pprcht/gfnff (the port source), every structure of MOR41 and GMTKN55: max 0.07 kcal/mol; S30L-CI: max 0.58, caused by one deliberate difference in the triple-bond torsion (0.05 with `-gfnff.storsion_reference_loop_bug true`). Against xtb it differs where xtb and the port source differ. |
| `gfn2` | GFN2-xTB, native | the accurate choice, about 100x slower than `gfnff` | tblite and xtb 6.7.1: 11 of 12 reference molecules within 1e-8 Eh, GMTKN55 within 0.017 kcal/mol |
| `gfn1` | GFN1-xTB, native | as `gfn2` | tblite: 14 of 16 reference molecules within 1e-8 Eh; GMTKN55 max 0.011 kcal/mol |
| `pm3`, `am1`, `mndo`, `pm6` | NDDO, native | semi-empirical | PM3/AM1/MNDO: 21 of 21 tests against Ulysses |
| `eht`, `uff`, `qmdff` | Hückel, force fields | | machine-tested only |
| `xtb-*`, `tblite-*`, `ugfn2`, `ipea1` | external backends | need `USE_XTB`, `USE_TBLITE`, `USE_ULYSSES` | production-quality interfaces |

Use `gfnff` to explore and `gfn2` to decide. Details per method and per data set:
[docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md), [docs/SQM_VALIDATION.md](docs/SQM_VALIDATION.md),
[docs/GMTKN55_VALIDATION.md](docs/GMTKN55_VALIDATION.md), [docs/MOR41_VALIDATION.md](docs/MOR41_VALIDATION.md).
The figures above are from September 2026 and are not updated automatically; the scripts in `scripts/`
re-measure them.

**Limits you should know.**
- GFN-FF is a force field. On the 41 MOR41 reactions it is 60 to 70 kcal/mol away from DLPNO-CCSD(T) in mean absolute deviation (GFN2: about 12). Port fidelity to the reference is high; the reference itself is not accurate for reaction energies of metal complexes.
- GFN-FF has no periodic boundary conditions.
- Open-shell systems work for `gfn1`/`gfn2` (`-spin N`, N = number of unpaired electrons).
- GPU backends (`-gpu cuda|rocm|vulkan`) are optional and off by default. CUDA has been run on an NVIDIA H200, RTX 5090, RTX 5080, RTX A4500 and GTX 1660; ROCm on a Radeon 890M (gfx1150) only, so the ROCm plugin build is unverified on other setups; Vulkan is opt-in and brings no speed-up. A GPU run that falls back to the CPU says so after the run (`-gpu_strict true` stops instead).
- Results of the native methods can differ in the last digits from tblite/xtb in documented places (`-xtb.sto6g_legacy_4sp`, `-xtb.d4_atm_cutoff 40.0` reproduce xtb's choices).
- Molecular dynamics times recorded with versions before the fix of 2026-09 (`ef462fcf`) run 1.9516 times too long; multiply them by that factor.
- GFN-FF writes its perceived topology to `<basename>.topo.json` next to the input and reads it again. Delete the file after changing the geometry file in a way that alters bonding.

**Capabilities** (each prints its options with `curcuma -<capability>` without arguments): single point
`-sp`, optimisation `-opt`, Hessian and frequencies `-hessian`, molecular dynamics `-md` (thermostats,
RATTLE, temperature ramps, PLUMED metadynamics), conformer search `-confsearch` and filtering `-confscan`,
RMSD with atom reordering `-rmsd`, docking `-dock`, interaction energies `-interaction`, implicit solvation
(ALPB/GBSA for native GFN1/GFN2, ALPB for GFN-FF), trajectory and structure analysis `-analysis`.

## Download and requirements

Prebuilt binaries are on the [Releases page](https://github.com/conradhuebler/curcuma/releases): a Linux x86_64 AppImage
(`curcuma-<version>-x86_64-Linux.AppImage`) and a Windows x86_64 archive (`curcuma-<version>-x86_64-Windows.zip`). Every release
is marked as a pre-release and is built by the CI from one commit; the releases named `Curcuma CI (feature/...)` are builds
of development branches.

To build from source:

```sh
git clone --recursive https://github.com/conradhuebler/curcuma
```

You need [CMake](https://cmake.org/download/) 3.18 or newer and a C++17 compiler (gcc, clang, icc, MinGW).
Dependencies are fetched by CMake (FetchContent); nothing needs to be initialised by hand.

- [LBFGSpp](https://github.com/conradhuebler/LBFGSpp), a fork of [yixuan/LBFGSpp](https://github.com/yixuan/LBFGSpp/) (the fork allows single-step optimisation without resetting the history)
- [simple-d3](https://github.com/dftd3/simple-dftd3) (D3) and [cpp-d4](https://github.com/conradhuebler/cpp-d4) (D4, fork)
- [CxxThreadPool](https://github.com/conradhuebler/CxxThreadPool), [fmt](https://github.com/fmtlib/fmt), [nlohmann/json](https://github.com/nlohmann/json)
- [Eigen](https://gitlab.com/libeigen/eigen) is not downloaded automatically; the build scripts in `scripts/` fetch and update it
- Optional: [xtb](https://github.com/grimme-lab/xtb) (`USE_XTB`), [tblite](https://github.com/tblite/tblite) (`USE_TBLITE`), [Ulysses](https://gitlab.com/siriius/ulysses) (`USE_ULYSSES`), [PLUMED](https://github.com/plumed/plumed2) (`USE_Plumed`, must be built manually). The native `gfn1`, `gfn2` and `gfnff` do not need xtb or tblite.

A C++/Eigen implementation of the Munkres algorithm (Hungarian method) based on [this workshop](https://brc2.com/the-algorithm-workshop/) is included.

### Build (Linux, macOS)

```sh
cd curcuma
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

### Build (Windows)

Clone outside the Windows system folders. Add the `bin` directory of your MinGW installation to `PATH`
and select the generator explicitly:

```sh
cmake .. -DCMAKE_BUILD_TYPE=Release -G "MinGW Makefiles"
```

Standard MinGW has no OpenMP. [w64devkit](https://github.com/skeeto/w64devkit) provides a GCC with
OpenMP: install [CMake](https://cmake.org/download/), extract w64devkit, start `w64devkit.exe`, and
inside that terminal run `git clone --recursive ...`, `mkdir build`, `cd build`,
`cmake .. -DCMAKE_BUILD_TYPE=Release -G "MinGW Makefiles"`, `mingw32-make`. CMake then picks up
`-DUSE_OpenMP=ON` with the right `libgomp`. The Windows nightly build uses `-DUSE_PORTABLE_MATH=ON`,
see [docs/PORTABLE_ERF.md](docs/PORTABLE_ERF.md); native Windows runs have not been compared with Linux
results on a Windows machine.

### GPU acceleration (optional)

One backend per build directory (`release_cuda/`, `release_rocm/`, `release_vulkan/`); the default build
needs none of these. The backends are loaded at run time as plugins (`libcurcuma_{cuda,rocm,vulkan}.so`),
so a CPU-only run does not touch them.

- **CUDA** (`-DUSE_CUDA=ON`, run with `-gpu cuda`): NVIDIA CUDA toolkit (`nvcc`, cuSOLVER, cuBLAS, cudart). [docs/SQM_GPU.md](docs/SQM_GPU.md)
- **ROCm / HIP** (`-DUSE_ROCM=ON -DCMAKE_PREFIX_PATH=/opt/rocm -DROCM_GPU_ARCH=gfxNNNN`): `hip-runtime-amd`, `rocm-llvm`, `rocm-device-libs`, `rocminfo`, `rocblas`, `rocsolver`. [docs/SQM_ROCM.md](docs/SQM_ROCM.md)
- **Vulkan** (`-DUSE_VULKAN=ON`): `vulkan-icd-loader`, `vulkan-headers` and a driver with `shaderFloat64`; no ROCm needed on AMD. [docs/SQM_VULKAN.md](docs/SQM_VULKAN.md)
- **Several GPUs**: `-gpu_device N` pins a run; batch runs spread over the visible devices (`-gpu_devices`, `-gpu_workers_per_device`); one large molecule can be split over several GPUs. [docs/MULTI_GPU.md](docs/MULTI_GPU.md), [docs/GPU_TUNING.md](docs/GPU_TUNING.md)
- **Tuning**: every performance knob is a CLI flag; `python scripts/tuning_sweep.py mol.xyz --method gfn2 [--gpu cuda]` measures them on your machine and checks that no setting changes the energy.

## Using curcuma

```sh
curcuma -sp water.xyz                       # single point, default method gfnff
curcuma -sp water.xyz -method gfn2          # accurate native GFN2-xTB
curcuma -opt water.xyz -method gfn2         # geometry optimisation
curcuma -hessian water.xyz -method gfn2     # frequencies
curcuma -md water.xyz -method gfnff         # molecular dynamics
curcuma -rmsd a.xyz b.xyz                   # RMSD (reorders atoms if needed)
curcuma -confsearch mol.xyz -md_method gfnff -opt_method gfn2
```

### Getting help

| Command | Output |
|---|---|
| `curcuma -help` | all capabilities by category |
| `curcuma -help <category>` | detailed help of one category (e.g. `optimization`, `dynamics`) |
| `curcuma -methods` | methods available in this build, GPU devices |
| `curcuma -list-modules` | modules with their parameter counts |
| `curcuma -help-module <module>` | every parameter of a module with type, default and description (e.g. `gfnff`) |
| `curcuma -export-config <module>` | the module's defaults as JSON |

The help text and the exported defaults come from the parameter definitions in the source, so they
describe the installed binary. The per-tool notes for RMSD, docking, ConfScan, trajectories,
optimisation, MD and ConfSearch (examples, options, convergence rules) are in
[docs/USAGE_TOOLS.md](docs/USAGE_TOOLS.md); its parameter lists date from 2025.

### Parameters, JSON, reproducible runs

Any parameter can be given as a flat flag (`-cn_cutoff_bohr 5.5`) or scoped (`-gfnff.cn_cutoff_bohr 5.5`).
`-export-config <module>` prints a module's defaults as JSON to the standard output (the program banner
precedes the JSON in the build checked on 2026-10-01, so remove the lines before the first `{` when saving it).

A complete run, with the resolved parameters, can be captured and replayed
([docs/CLI_ROUND_TRIP.md](docs/CLI_ROUND_TRIP.md)):

```sh
curcuma -sp water.xyz -method gfnff -cn_cutoff_bohr 5.5 -export_run run.json
curcuma -import_config run.json                       # replay
curcuma -import_config run.json -cn_cutoff_bohr 7.0   # replay with override
```

### Output directories (BMT)

Every command writes into a directory named `Basename.Method.Timestamp`, for example
`water.opt.20261001_091303/`, with the output files and a `metadata.json`. `-bak FILE` copies a file back
to the working directory, `-no_bmt` writes into the working directory.

```sh
curcuma -opt water.xyz -method gfnff -bak water.opt.xyz
```

### Stopping a run

On Linux, Ctrl-C makes curcuma create an empty file `stop`. ConfScan and MD check for it,
write a restart file and finish. A second Ctrl-C while `stop` exists ends the program at once.

### More

[docs/](docs/) holds one document per feature: the SCF modes ([docs/SCF_MODES.md](docs/SCF_MODES.md)),
solvation ([docs/SOLVATION.md](docs/SOLVATION.md)), temperature ramps
([docs/TEMPERATURE_RAMP.md](docs/TEMPERATURE_RAMP.md)), large-system MD
([docs/MD_LARGE_SYSTEMS.md](docs/MD_LARGE_SYSTEMS.md)), PLUMED ([docs/PLUMED_HELP.md](docs/PLUMED_HELP.md)),
reactive GFN-FF MD ([docs/GFNFF_REACT_TOPOLOGY.md](docs/GFNFF_REACT_TOPOLOGY.md), experimental), benchmark
sets ([docs/TESTSET_RETRIEVAL.md](docs/TESTSET_RETRIEVAL.md)). Method status paragraphs that used to stand in
this README are in [docs/README_METHOD_NOTES_2026.md](docs/README_METHOD_NOTES_2026.md).

## For developers

New capabilities define their parameters with the Parameter Registry (`BEGIN_PARAMETER_DEFINITION` /
`PARAM` in the capability header); the build extracts them into the help, the JSON export and the
validation. See [docs/PARAMETER_SYSTEM.md](docs/PARAMETER_SYSTEM.md),
[docs/archive/PARAMETER_MIGRATION_GUIDE.md](docs/archive/PARAMETER_MIGRATION_GUIDE.md) and the reference
implementation `src/capabilities/analysis.h`. Each source directory has a `CLAUDE.md` with its design notes;
the root [CLAUDE.md](CLAUDE.md) holds the project rules, the open items and the validation traps, and
[docs/KNOWN_ISSUES_ARCHIVE.md](docs/KNOWN_ISSUES_ARCHIVE.md) the bug and validation history.

```sh
cd build && ctest --output-on-failure    # test suite
```

## Please cite

Software: [conradhuebler/curcuma, Zenodo](https://doi.org/10.5281/zenodo.4302722).
For the conformer filter (`-confscan`): C. Hübler, *A conformational filter protocol for structures with topological symmetry*, ChemRxiv preprint (2026), [doi:10.26434/chemrxiv.15009180/v1](https://chemrxiv.org/doi/full/10.26434/chemrxiv.15009180/v1).
Please also cite the methods and programs you use; the most recent information is on their project pages.

- UFF: A. K. Rappe et al., J. Am. Chem. Soc. 114 (1992) 10024-10035. With the H4 hydrogen-bond correction (J. Chem. Theory Comput. 8, 141-151 (2012)), same parameters as for PM6-D3 for now.
- GFN-FF: S. Spicher and S. Grimme, Angew. Chem. Int. Ed. 59 (2020) 15665. DOI: 10.1002/anie.202004239
- GFN1-xTB: S. Grimme, C. Bannwarth and P. Shushkov, J. Chem. Theory Comput. 13 (2017) 1989.
- GFN2-xTB: C. Bannwarth, S. Ehlert and S. Grimme, J. Chem. Theory Comput. 15 (2019) 1652.
- D3: J. Chem. Phys. 132, 154104 (2010), https://doi.org/10.1063/1.3382344
- D4: E. Caldeweyher et al., J. Chem. Phys. 147 (2017) 034112, DOI: 10.1063/1.4993215; J. Chem. Phys. 150 (2019) 154122, DOI: 10.1063/1.5090222
- Atom reordering with molalign (optional, `-reorder -method molalign -molalignbin PATH`): [J. Chem. Inf. Model. 2023, 63, 4, 1157-1165](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.2c01187)

Curcuma prints the references for the methods used in a run and writes them to `curcuma_citations.bib`.

## Funding

The development of curcuma is funded by 2026 Stiftung Innovation in der Hochschullehre.

![STIL Logo](https://github.com/conradhuebler/curcuma/raw/master/STIL_Funding.jpg)
