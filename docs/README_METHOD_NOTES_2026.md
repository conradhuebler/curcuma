# Method notes moved from the README (as of 2026-10-01)

Per-method status paragraphs that stood in the README under "UFF, xTB, GFN-FF and Dispersion Correction"
and "Native GFN-FF Status (April 2026)", moved unchanged. Dated statements are as old as their date; the
current status is in [GFNFF_STATUS.md](GFNFF_STATUS.md), [SQM_VALIDATION.md](SQM_VALIDATION.md) and
[KNOWN_ISSUES_ARCHIVE.md](KNOWN_ISSUES_ARCHIVE.md).

Statements in this text that later work contradicts:
- "Organometallics / transition metals: No test molecule with metal center" and "Metal-specific EEQ
  corrections (fqq) not implemented": MOR41 (95 structures, transition metals) is the main GFN-FF
  validation set since July 2026 (Known Issues #5, #6); GFN-FF itself remains a force field, 60-70 kcal/mol
  from DLPNO on MOR41 reactions.
- "20 validation molecules": the GFN-FF validation now covers MOR41, GMTKN55 and S30L-CI.
- "15/16 validation molecules at 1e-8" for gfn2 and "14/16" for gfn1: `docs/SQM_VALIDATION.md` carries the
  current table (root CLAUDE.md: gfn2 11/12 with `complex` open).

### UFF, xTB, GFN-FF and Dispersion Correction
Curcuma has an interface to tblite, xtb as well simple-d3 and cpp-d4, enabling semiempirical calculations or combinations of UFF with D3, D4 and H4 (no parameters are adjusted yet). To use one of the methods, please add **-method methodname** to your arguments:

Classical force field:
- uff : Universal Force Field (no longer any capability's default since Sep 2026)

Native force field (no external dependency required, **the default for every capability**):
- **gfnff** : Native C++ GFN-FF — full energy and gradient, validated against Fortran reference (see status below)
- **gfnff** + `-gpu cuda` : CUDA-accelerated variant; topology cached, charges on CPU, all kernels on GPU
- **xtb-gfnff** : GFN-FF via the xtb Fortran library (USE_GFNFF build flag)

**Which method should I use?** `gfnff` is the **fast** one and the default for every
capability (single point, optimisation, MD, Hessian, conformer search): a native GFN-FF
force field, no external dependency, milliseconds per gradient. `gfn2` is the **accurate**
one: native GFN2-xTB, semi-empirical QM, roughly two orders of magnitude slower but with
real electronic structure (charges, orbitals, bond breaking). Use `gfnff` to explore and
`gfn2` to decide.

> **Charged multi-fragment systems (GFN-FF, Sep 2026):** the net charge now goes to the chemically right fragment instead of the one that happens to contain atom 1 (`-gfnff.frag_charge_model ensemble`, the new default; `reference` restores the old rule). Energies of charged species that GFN-FF sees as several fragments can therefore differ from earlier releases and from xtb (GMTKN55 WATER27 reaction MAD 58.6 -> 21.4 kcal/mol). See [docs/FRAG_CHARGE_MODEL.md](FRAG_CHARGE_MODEL.md).

Native GFN methods (no external dependency required, canonical backends since AP3 2026-04-25):
- **gfn1** : Native GFN1-xTB — 14/16 validation molecules at 1e-8 vs tblite; includes the GFN1-only halogen-bond correction (B–X···A, added Sep 2026)
- **gfn2** : Native GFN2-xTB — 15/16 validation molecules at 1e-8 vs tblite (only `complex` open at 7.3e-8)

> Native GFN1/GFN2 are validated against tblite to a 1e-8 Eh target — see [docs/SQM_VALIDATION.md](SQM_VALIDATION.md). For explicit tblite or xtb backends use `tblite-gfn1`/`tblite-gfn2` or `xtb-gfn1`/`xtb-gfn2`.

> **Halogen bonds (GFN1, Sep 2026):** GFN1 carries a classical B–X···A correction (X = Cl/Br/I/At, acceptor = N/O/P/S) that GFN2 does not. It was previously unimplemented; with it, all 2462 GMTKN55 structures reproduce xtb 6.7.1 to MAD 0.00007 / max 0.011 kcal/mol (was 0.041 / 11.93, and every deviation above 0.1 kcal was a halogen-bonded `HAL59` structure). See [docs/GMTKN55_VALIDATION.md](GMTKN55_VALIDATION.md).

> **4th-period elements (GFN1):** what is left of that 0.011 kcal/mol is a genuine xtb-vs-tblite disagreement, not a curcuma error — the two references carry different STO-6G 4s/4p tables (Z = 19–36; GFN2 uses STO-4G there and is unaffected). curcuma follows tblite, whose expansion fits the exact Slater function 3–5× better. `-xtb.sto6g_legacy_4sp true` switches to xtb's tables and reproduces the binary bit-for-bit. Details in [docs/GMTKN55_VALIDATION.md](GMTKN55_VALIDATION.md).

> **Gradients (Sep 2026):** analytic gradients are now validated set-wide against xtb 6.7.1 on all 2462 GMTKN55 geometries (median deviation 3e-7 / 4e-7 / 4e-8 Eh/Bohr for gfn1 / gfn2 / gfnff), with the outliers arbitrated by finite differences of each code's own energy. That sweep found and fixed two unit bugs — GFN-FF MD forces were a factor 1.89 too small, and vibrational frequencies were too high for every method — see [docs/GRADIENT_VALIDATION.md](GRADIENT_VALIDATION.md). Frequencies now match xtb to ≤0.13 % on H2O for all three methods.

> **GFN-FF MD/optimisation correctness (Sep 2026):** several non-bonded pair lists were built once at setup and never revisited during a run — a non-bonded repulsion pair starting more than 20 Bohr apart could diffuse to near-zero distance with **zero** repulsive force (root-caused from a real crash on an H200), and D4 dispersion's C6 values were frozen at the setup geometry for every GFN-FF MD/optimisation run, not only large ones, while the gradient used the current CN — energy and gradient belonged to different functions after the first geometry change. Both fixed with a periodic geometry-triggered refresh; MOR41 (95/95) and GMTKN55 gfnff (2462/2462) bit-identical before/after. See [docs/GFNFF_PAIR_LIST_REFRESH.md](GFNFF_PAIR_LIST_REFRESH.md).

> **Speed:** on a 231-atom complex (single core, energy+gradient) native `gfn1` runs in ~1.02 s and `gfn2` in ~1.08 s, versus xtb 6.7.1 at 1.37 s / 0.98 s — i.e. gfn1 is faster than xtb and gfn2 within ~11%. See [docs/SQM_PERFORMANCE.md](SQM_PERFORMANCE.md) for the single-core record and [docs/SQM_THREADING.md](SQM_THREADING.md) for `-threads N` scaling.

> **d-shell elements (X-I1, June 2026):** native GFN1/GFN2 now handle d-shell basis functions (S, P, Cl, Si and other main-group d elements), matching tblite to ≤1e-8 Eh; analytic gradients FD-validated. CPU only — on `-gpu` a d-shell system falls back to the CPU integral/SCF path. Transition metals: after the Jul 2026 fixes (shell-vs-angular parameter indexing + 6s/6p STO-6G expansion), **native GFN1 and GFN2 reproduce tblite for transition metals** — GFN2 3d exact (1e-8), GFN1 72/95 MOR41 structures exact; both leave a small **~1e-3 Eh** residual for 4d/5d (heavy-element band/multipole/D4, still open). **GFN-FF transition metals are not yet validated.** See [docs/SQM_DSHELL_WP.md](SQM_DSHELL_WP.md) and [docs/MOR41_VALIDATION.md](MOR41_VALIDATION.md).

> Native GFN1/GFN2 can use multiple cores **within one calculation** of a single large molecule: pass `-threads N` to a `-sp`/`-opt`/MD run (default is serial and bit-identical). Integral setup, gradient and Fock build scale ~3–5×; see [docs/SQM_THREADING.md](SQM_THREADING.md).

> **Benchmark test sets on demand:** `python scripts/fetch_testset.py fetch mor41` downloads the Grimme-group MOR41/GMTKN55/S30L benchmark sets into the layout the validation scripts expect (S30L's Supporting Information is paywalled and must be placed by hand; instructions are printed). `scripts/testset_perf.py` then times CPU/threading/GPU performance on whatever set is fetched. See [docs/TESTSET_RETRIEVAL.md](TESTSET_RETRIEVAL.md).

> Opt-in **MKL-free / GPU-portable eigensolve kernels** are available for the native GFN SCF (MKL stays the default): `-eigensolver native` (own Householder + Cuppen divide-and-conquer), `-eigensolver purify` (0 K density-matrix purification, GEMM-only, no diagonalization), `-eigensolver lobpcg` (seeded block LOBPCG, experimental), and `CURCUMA_EIG_TRED2=blocked` (BLAS-3 blocked tridiagonalization). See [docs/SQM_EIGENSOLVE_GPU.md](SQM_EIGENSOLVE_GPU.md).

> Opt-in **CUDA GPU path** for the native GFN1/GFN2 solver: `-method gfn1|gfn2 -gpu cuda` (build `release_cuda/` with `-DUSE_CUDA=ON`). Staged cuSOLVER/cuBLAS port (the CPU path is unchanged and `#ifdef`-free); both **GFN1** and **GFN2** run a device-resident SCF under the default Broyden mixing, and **Stage 3 builds the integrals (CN/S/H0/L/γ/multipole) on the device and Stage 4 the nuclear gradient — so `-opt`/`-md` are fully device-resident** (only xyz up, gradient+energy down per step; every device kernel matches the CPU elementwise to ~1e-15). 🤖 AI-generated / ⚙️ machine-tested only. See [docs/SQM_GPU.md](SQM_GPU.md).

> **AMD/ROCm** (`-gpu rocm`, build `release_rocm/` with `-DUSE_ROCM=ON`) and **Vulkan compute** (`-gpu vulkan`, hand-written SPIR-V, `-DUSE_VULKAN=ON`) backends for the same `gfn1`/`gfn2`/`gfnff` methods. `-gpu auto` picks the first compiled backend (cuda > rocm > vulkan), else CPU. 🤖 **Vulkan: GFN1 = Stage 2** (device-resident SCF — Fock/eigensolve/density/populations/band on the GPU via a device-built Löwdin S⁻¹ᐟ²; only `v_ao`/`occ` up and `eps`/`pop`/`band` down per iteration), **GFN2 = Stage 1** (per-iteration eigensolve on GPU). gfn1/gfn2 single-point + opt match the CPU bit-for-bit on the validation set (AMD 890M/RADV); integrals/gradient still CPU. **ROCm: GFN1 = Stage 4 (fully device-resident)** — the integral build (CN/S/H0/L/γ), the SCF (Fock/density/eigensolve via HIP kernels + rocBLAS + rocSOLVER) and the nuclear gradient (repulsion/Pulay/Coulomb HIP kernels) all run on the GPU; only the dispersion gradient + CN chain-rule on the host. **GFN2** uses the device integrals + rocSOLVER eigensolver (gradient on host). gfn1/gfn2 single-point + opt match the CPU bit-for-bit, incl. the full `-opt` trajectory (AMD 890M, needs `rocsolver`+`rocblas`). **ROCm GFN-FF** (`-DUSE_ROCM=ON`, June 2026): the full energy + nuclear-gradient kernel stack runs on the GPU (single-TU hipify of the CUDA gfnff kernels; EEQ via rocSOLVER `dpotrf`/`dgetrf` + host CPU-Schur); single-point energy and gradient match CPU ≤1e-7 on water/CH4/caffeine/231-atom complex. **Two opt-in CUDA-only GFN-FF GPU flags (default OFF, ROCm mirrors pending):** `-gfnff.eeq_mixed_precision` (FP32-factor + FP64-refine EEQ solve) and `-gfnff.gpu_disp_pairs_on_device` (on-device D4 pair build) — bit-identical to the host but not a measured speedup (residency milestones). See [docs/SQM_ROCM.md](SQM_ROCM.md) / [docs/SQM_VULKAN.md](SQM_VULKAN.md) / [docs/GFNFF_PERFORMANCE_LEVERS.md](GFNFF_PERFORMANCE_LEVERS.md).

> Opt-in **approximate large-system modes** scale the native GFN SCF beyond ~1000 atoms by exploiting locality (default is the exact dense path): `-large_system_mode fragments` (disconnected-fragment SCF, energy+gradient, `-eigensolver` propagates per fragment), `-large_system_mode dc` (divide-and-conquer, energy-only, `-eigensolver` propagates per sub-block, `-large_system_buffer_bohr` accuracy knob), `-large_system_mode sparse` (non-orthogonal density purification, 0 K gapped, `-eigensolver` ignored, `-large_system_sparse_threshold` knob). Each converges to the dense energy as its knob tightens; combining `-large_system_mode=fragments|dc` with `-eigensolver=purify` requires `-electronic_temperature 0` (hard error otherwise). See [docs/SQM_LARGE_SYSTEMS.md](SQM_LARGE_SYSTEMS.md).

> Opt-in **multi-step SCC extrapolation** for the native GFN SCF cuts SCF iterations across geometry steps in `-opt`/`-md` by predicting the next charge state from several past converged steps (generalises the 1-step warm-start; default `none` is unchanged). `-scf_extrapolation aspc` (Kolafa ASPC, best for fixed-timestep MD) or `-scf_extrapolation gauss` (least-squares, better for irregular opt steps), with `-scf_extrapolation_order`. The safe default `guess` coupling still converges the SCF fully; `-scf_extrapolation_apply xlbomd` is an experimental extended-Lagrangian Born-Oppenheimer mode (time-reversible auxiliary density + converged corrector, for low MD energy drift). On a smooth caffeine trajectory, `aspc`/`gauss` roughly halve SCF iterations (gfn2 215→90, gfn1 170→79) with bit-identical converged energy. 🤖 AI-generated / ⚙️ machine-tested only. See [docs/SQM_SCF_EXTRAPOLATION.md](SQM_SCF_EXTRAPOLATION.md).

xtb methods:
- xtb-gfnff : GFN-FF via the xtb library
- xtb-gfn1
- xtb-gfn2

Using only **d3** or **d4** should be possible.

Native GFN2 includes an analytic D4 dispersion charge-response gradient
(∂E_D4/∂q · ∂q/∂x). The zeta charges default to a single-shot dftd4 EEQ model
(`-d4_charge_source eeq`, analytic ∂q/∂x); `-d4_charge_source mulliken` feeds the
GFN2 SCF charges (energy + ∂E/∂q; the CPSCF gradient response is still pending —
see [docs/D4_Q_RESPONSE.md](D4_Q_RESPONSE.md)). Current alignment vs tblite:
11/12 at 1e-8 — only `complex` (231 atoms, 6.95e-5 Eh residual) remains open.
Status tracked in [docs/GFN2_NATIVE_ROADMAP.md](GFN2_NATIVE_ROADMAP.md) and
[docs/GFN2_D4_STATUS.md](GFN2_D4_STATUS.md); `ctest -L d4_diag`.

The native GFN SCF defaults to `broyden` mixing — a modified-Broyden quasi-Newton
scheme on the SCC charge vector, the same mixer tblite/xtb use — which converges
large polar systems that the old Fock-DIIS diverged on (e.g. the 231-atom
`complex` now converges from the bare guess with plain `-method gfn2`). Other
modes remain selectable: `-scf_mode diis|plain|level-shift` and `-scf_guess
h0|eeq` (plus `-scf_damping`, `-diis_start`, `-level_shift`). See
[docs/SCF_MODES.md](SCF_MODES.md).

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
- **GFN-FF ALPB solvation**: self-consistent Born reaction field coupled into EEQ (`A_eeq += B`); validated against xtb 6.7.1 (`--gfnff --alpb`) to ≤1e-8 Eh (7 molecules × 4 solvents, `ctest -L gfnff_solvation`; June 2026). Gradient FD-validated at frozen solvated charges (same approximation as the Fortran reference). `-gfnff.solvent_model gbsa` maps to ALPB (GFN-FF has no separate GBSA model; warns at runtime). See [docs/SQM_SOLVATION_WP.md](SQM_SOLVATION_WP.md).

**Not validated / not implemented:**
- **Periodic boundary conditions**: Not implemented
- **Organometallics / transition metals**: No test molecule with metal center; parameter quality unknown

**Large-system precision (re-verified Sep 2026, superseding an older note)**: the two caveats
previously listed here — "dispersion gradients show √N accumulation error on large systems"
and "GPU energy for polymer (1280 atoms): 8.9 µEh vs. 1 µEh tolerance" — no longer reproduce.
`test_gfnff_validation` on the current 1410-atom `polymer.xyz` (`ctest -R gfnff_val_polymer`):
dispersion GradComp max_err 9.9e-9 Eh/Bohr (tol 1e-4, was ~4.9e-4 in Mar 2026 — likely fixed
incidentally by later D3/D4 precision work, e.g. docs/KNOWN_ISSUES_ARCHIVE.md Known Issues #5). CPU-vs-GPU
single-point energy on the same molecule, ROCm (gfx1150): 0.33 µEh (well under the 1 µEh
target; CUDA hardware was not available to re-check that backend directly).

**Reactive MD (experimental)**: `-gfnff.topology_mode react` lets bonds form and break during MD (hysteresis re-detection + bonded-term rebuild, NVT-only) — see [docs/GFNFF_REACT_TOPOLOGY.md](GFNFF_REACT_TOPOLOGY.md).

**Cross-platform determinism (`-DUSE_PORTABLE_MATH=ON`)**: Wine and native Windows can round `erf`/`acos`/`exp`/`log` differently in the last bit (different CRT-DLL reimplementations), which can flip a GFN-FF classification threshold into a different bond term. Vendored fdlibm-derived replacements close this; off by default, on for the Windows nightly build — see [docs/PORTABLE_ERF.md](PORTABLE_ERF.md).

**One unit system (`-DUSE_LEGACY_UNIT_CONSTANTS=ON` to revert)**: every Bohr/Ångström and Hartree conversion uses CODATA 2018 (`src/core/units.h`). Before Sep 2026 seven different Bohr radii were in use, which put a systematic 5e-7 Eh between CPU and GPU GFN-FF on a 7320-atom system; the legacy build restores the old per-site values bit for bit — see [docs/UNIT_CONSTANTS.md](UNIT_CONSTANTS.md).

**Known differences from Fortran reference** (see [docs/GFNFF_STATUS.md](GFNFF_STATUS.md)):
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
