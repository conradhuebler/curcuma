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

Native ab-initio QM engine (**`hf` and `hf-3c` usable, LDA/PBE/B3LYP not yet** — no external dependency required; Gaussian-basis 1e integrals + 4-centre ERIs (McMurchie-Davidson, shell-quartet blocked, Schwarz-screened, OpenMP) and a closed-shell RHF SCF, integral-direct above ~150 basis functions (`-qm.scf_direct`); analytic gradient, so `-opt` works; H-Ne):
- **hf** : closed-shell RHF, matches ORCA 6.1 HF/def2-SVP on 10/10 validation molecules within 4e-9 Eh; engine settings under `-qm.*` (e.g. `-qm.basis MINIX -qm.threads 4`; the old `-dft.*` scope still works); since Oct 25, 2026 a faster CPU integral kernel (Hermite-matrix McMurchie-Davidson, energies unchanged to 10 decimals): C80H162/def2-SVP 712 -> 213 s on 24 threads (ORCA 6.1: 471 s), anthracene 219 -> 48 s on one thread (ORCA: 43 s, at the tighter convergence test of the comparison), see [docs/QM_TIMING_ORCA_PYSCF_2026-10.md](docs/QM_TIMING_ORCA_PYSCF_2026-10.md)
- **hf-3c** : native HF-3c = HF/MINIX + D3(BJ) + gCP + SRB (Sure & Grimme 2013); every term matches PySCF + simple-dftd3 to ≤4e-11 Eh on 17 molecules, H2O total within 4.3e-10 Eh of ORCA 6.1; analytic gradient ≤2e-10 Eh/Bohr vs PySCF + simple-dftd3, `-opt -method hf-3c` reaches their minimum to 1e-6 Å. The external ORCA version is `orca-hf-3c`
- **GPU** (`-gpu cuda`): ERI tensor and J/K on the device; when the tensor does not fit, integral-direct J/K on the device (`-qm.gpu_direct auto|on|off`); from 200 basis functions the whole SCF loop stays on the device (`-qm.gpu_scf auto|on|off`, small gain, J/K dominate); both parts of the analytic gradient run on the device too (`-qm.gpu_grad auto|on|off` two-electron, `-qm.gpu_grad1e auto|on|off` one-electron; benzene HF/def2-SVP `-opt` 10 s vs 77 s on 24 CPU threads); exact integral-direct SCF of C80H162 (1930 functions) in 130 s and C120H242 (2890) in 306 s on one RTX 5090 (Oct 19, 2026: sorted shell-pair lists, `-qm.scf_guess atom` start density, `-qm.scf_early_screening`, `-qm.eri_prim_density`; energies unchanged to 10 decimals); since Oct 20, 2026 thread-per-quartet McMurchie-Davidson kernels for all integral classes (SCF loop 4.2x faster, C120H242 116 s, energies unchanged to 10 decimals) and FP32 integrals in the early SCF phase (`-qm.gpu_eri_precision fp64|fp32|mixed|auto`, `-qm.gpu_fp32_switch`; auto = mixed on GeForce-class cards from 150 functions; convergence and the final energy from FP64 builds, within 2e-10 Eh of FP64; another 10-23 % of the SCF loop, C120H242 112 s); AI-generated, machine-tested on one RTX 5090, see the GPU roadmap below
- **Density fitting** (`-qm.ri j|jk`, default off): RI-J (Coulomb fitted, exchange exact) or RI-JK (both fitted) with def2/J resp. def2/JK (H-Ne); a method variant (RI error 5e-6 to 4e-5 Eh per atom on the test molecules), energies match PySCF to 1e-10 Eh and ORCA 6.1 to 2e-9 Eh; RI-JK makes HF several times faster on the CPU (naphthalene 25 vs 406 ms per Fock build on 16 threads) with O(n^2 N_aux) memory (3-index tensor stored as packed n(n+1)/2 pairs); analytic RI gradient for `-opt` (vs PySCF <= 6.1e-11 Eh/Bohr, anthracene RI-JK 2e gradient 0.6 s vs 43 s exact on 16 threads); contractions on the GPU with `-gpu cuda` (metric solve optionally too, `-qm.gpu_ri_fit on`); heptacene (510 basis functions) 19 s on 16 threads, 12 s with the GPU, C40H82 (970) 62 s with the tensor on a 32 GB card; resource-saving modes since Oct 19 2026: `-qm.ri_storage batched` (automatic above `-qm.ri_max_memory_mb`; recomputed 3-index integrals, K from half-transformed integrals, same energies to 1e-10 Eh, several times slower per build), `-qm.gpu_ri_stream` (tensor in host memory, streamed to the card), `-qm.ri_mixed_precision` (FP32 tensor, FP64 corrections, final energy FP64) and `-qm.ri_fit_radius` (local exchange fit, an approximation with a measured error); AI-generated, machine-tested ([WP10](docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md))
- **ri-mp2**, **scs-ri-mp2**, **sos-ri-mp2** : density-fitted closed-shell MP2 on the native RHF with def2-universal-JKFIT (default) or -JFIT as correlation fitting basis (`-qm.ri_mp2_aux`; no def2-SVP/C shipped, so a method variant, not RI-MP2/def2-SVP/C); within 1.7e-9 Eh of PySCF (same aux files) and 1.8e-9 Eh of ORCA 6.1 (AuxC def2/JK); RI error vs canonical 3e-6 to 3e-5 Eh per atom; anthracene correlation part 0.75 s vs 3.7 s canonical (16 threads) with <1 GB instead of 31 GB; `-qm.ri jk` for an RI SCF, `-gpu cuda` for the energy step; energy only; AI-generated, machine-tested ([WP10/WP11 doc](docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md))
- **lmp2** (alias `local-ri-mp2`) : local PNO-based RI-MP2 on the native RHF (Pipek-Mezey or Boys orbitals, Boughton-Pulay + 7 Bohr PAO domains, T_CutPNO 1e-8 with semicanonical correction, 15 Bohr pair cut with dipole estimate; `-qm.lmp2_*`, aux `-qm.ri_mp2_aux`); with complete spaces identical to RI-MP2 (2e-12 Eh); default thresholds 3.7e-6 to 2.4e-5 Eh per atom against RI-MP2 on alkanes C4-C40, water clusters up to 20 and acenes up to anthracene (tetracene/pentacene need `-qm.lmp2_pair_cut 25`); faster than RI-MP2 from about C16H34 (C40H82 12 s vs 696 s on 24 threads), since Oct 18 2026 local density fitting with integral-direct 3-index integrals by default (`-qm.lmp2_fit local|global`, `-qm.lmp2_fit_radius` 8 Bohr, `-qm.lmp2_fit_domain orbital|pair`, `-qm.lmp2_ao_cut`; adds <= 1e-6 Eh per atom on the test sets), no stored RI tensor: C80H162 (1930 basis functions) LMP2 in 39 s and 7 GB with an exact GPU direct SCF (`-qm.ri off -gpu cuda`, 1521 s), before 93 s on a 133 GB RI-JK tensor (the SCF limits the size); localisation sweeps `-qm.lmp2_localizer_algo fast|classic` (fast default since Oct 19 2026, same energies, C80H162 7.8 s -> 0.4 s); energy only, closed shell; AI-generated, machine-tested, see [WP15](docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md)
- **ccsd** : canonical closed-shell CCSD (RCCSD) on the native RHF, exact integrals, CPU (`-qm.frozen_core`, `-qm.cc_*`); prints E_corr(CCSD), the MP2 energy of iteration 0, iterations and the T1 diagnostic; def2-SVP correlation energies within 3.5e-10 Eh of PySCF and 1.5e-9 Eh of ORCA 6.1 on 7 molecules up to benzene (full and frozen core); energy only, no RI, needs memory for the v^4 integrals; AI-generated, machine-tested, see [WP13](docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md)
- **mp2**, **scs-mp2**, **sos-mp2** : canonical closed-shell MP2 on the native RHF (no RI, `-qm.frozen_core true` for the chemical frozen core, `-gpu cuda` runs the transformation on the device); def2-SVP correlation energies within 1.8e-9 Eh of PySCF and 2.0e-9 Eh of ORCA 6.1 on 8 molecules (full and frozen core); analytic gradient (relaxed density via Z-vector, frozen core, `-sp -gradient`, `-opt`; 2e derivative part also on CUDA) within 1.4e-9 Eh/Bohr of the finite difference of PySCF's MP2 energy and 7.7e-8 of ORCA 6.1; AI-generated, machine-tested, see [WP11](docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md)
- **lda**, **pbe**, **b3lyp** : method names reserved, currently the nuclear repulsion only (no V_xc yet)
- see [docs/NATIVE_QM_IMPLEMENTATION.md](docs/NATIVE_QM_IMPLEMENTATION.md), the [DFT roadmap](docs/DFT_ROADMAP/) and the [performance/GPU roadmap](docs/QM_GPU_ROADMAP.md)

> Native GFN1/GFN2 are validated against tblite to a 1e-8 Eh target — see [docs/SQM_VALIDATION.md](docs/SQM_VALIDATION.md). For explicit tblite or xtb backends use `tblite-gfn1`/`tblite-gfn2` or `xtb-gfn1`/`xtb-gfn2`.
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
- **GFN-FF EEQ solver**: `-eeq_solver.solve_method auto|cholesky|ppcg|pcg|lu|ldlt|batched` (alias `-gfnff.solve`) picks the charge solver; `auto` switches to `ppcg` at 500 atoms. The names mean the same on the CPU and with `-gpu cuda`; `-gpu rocm` accepts only `auto` and `cholesky`. `batched` drops cross-fragment Coulomb (inexact). [docs/EEQ_SOLVER_CONTROL.md](docs/EEQ_SOLVER_CONTROL.md)
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

**Container walls**: `-md.wall_potential harmonic|logfermi|pbc`. The first two push an
escaping atom back and do work on it; the harmonic force constant is `wall_temp * kB`, far
too soft at its 298 K default to hold a hot gas. `pbc` exerts no force and instead moves a
whole molecule that left the container back in on the opposite side, so it neither heats nor
lets anything out. It is a container, not a periodic cell: the energy has no minimum image,
so interactions are not continued across the boundary (see
[docs/WP-PERIODIC-NONBONDED.md](docs/WP-PERIODIC-NONBONDED.md)).

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
the opt-in GFN1/GFN2 SCF preconditioner `-scf_precond frag`, with segmentation of large molecules ([docs/SCF_PRECOND.md](docs/SCF_PRECOND.md), experimental),
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
