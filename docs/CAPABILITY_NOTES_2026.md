# Capability notes 2026 (long form)

The detailed capability list and the "Completed Developments (2026)" list as they stood in the root
`CLAUDE.md` until 2026-10-01, moved here unchanged. The root `CLAUDE.md` now carries a short table.
Measured numbers in this text are dated by the entry they come from and are not kept current; the
campaign history behind them is in [KNOWN_ISSUES_ARCHIVE.md](KNOWN_ISSUES_ARCHIVE.md).

Known inconsistency inside this text: the d-shell bullet (June 2026) calls transition metals
"unvalidated"; Known Issue #5 (Jul 16, 2026) later reproduced tblite for all 95 MOR41 structures.

## Current Capabilities

### 1. Quantum Mechanical Methods

#### Native Implementations (Educational, No External Dependencies)

> **Default method is `gfnff` for every capability** (`-sp`, `-opt`, `-md`, `-hessian`,
> `-confsearch`, `-casino`): it is the **fast** general-purpose choice. **`gfn2` is the
> accurate one** — native GFN2-xTB, ~100x slower, but with electronic structure. `uff`
> remains available but is no longer any capability's default (Sep 2026). The one place
> `uff` stays the default is `ForceFieldGenerator`'s own `method` parameter, which selects
> which *force field* the FF engine builds (uff | uff-d3 | qmdff | cg) — GFN-FF does not
> go through that class, and its `ff_type = 4` branch there is dead code since the legacy
> `ForceField` GFN-FF path was removed.

> ⚠️ **All native QM methods are AI-implemented and machine-tested only — not human production tested.**
> Results should be validated against external references (TBLite, Ulysses, XTB) before use in research.

- ⚠️ **Extended Hückel Theory (EHT)** - AI-implemented, machine-tested
- ⚠️ **GFN2-xTB (Native)** - AI-implemented, machine-tested; canonical `gfn2` backend; `-opt` works; Broyden SCF default (`-scf_mode diis|plain|level-shift`, `-scf_guess h0|eeq|fragments`, the last opt-in for clusters of small molecules); a NaN or unconverged SCF is an error since Sep 27, 2026 (`-scf_allow_unconverged true` accepts the latter); 11/12 sqm_reference @1e-8 vs tblite (only 231-atom `complex` open) — [docs/SQM_VALIDATION.md](SQM_VALIDATION.md), [docs/SCF_MODES.md](SCF_MODES.md)
  - **d-shell support (X-I1, June 2026)**: S/P/Cl/Si/… (main-group d) now compute via the cartesian→spherical transform, ≤1e-8 Eh vs tblite. CPU + **CUDA GPU** (validated on a GTX 1660: energy bit-identical to CPU, gradient ~1e-16) + **ROCm GPU** (B6 port, validated on a Radeon 890M/gfx1150: energy bit-identical to CPU at 8 dp, `-opt` tracks the CPU trajectory); **Vulkan d still falls back to CPU** (GLSL shaders pending). Transition metals enabled but **unvalidated**. See [docs/SQM_DSHELL_WP.md](SQM_DSHELL_WP.md)
  - **Threading**: intra-molecule `-threads N` (setup 4×, gradient 3.6×), which since Sep 2026 drives the eigensolve too (`-eigensolver_max_threads N` caps it separately) — [docs/SQM_THREADING.md](SQM_THREADING.md)
  - **Tuning knobs are CLI flags** (Sep 2026): `-scf_reduce`, `-eigensolver_max_threads`, `-scf_fp32_stall_patience`, `-scf_fp32_false_fixpoint_factor`, `-gpu_multipole_otf` — the former environment variables still work and still win. `scripts/tuning_sweep.py` scans them on the machine it runs on, checks every result's energy against the baseline and prints a command line — [docs/GPU_TUNING.md](GPU_TUNING.md) sections 4-5
  - **Integral setup** (Jul 2026): shell-pair-blocked overlap/multipole/gradient kernels — setup 209→91 ms, gfn2 total 1199→1083 ms on complex/231 — [docs/SQM_PERFORMANCE.md](SQM_PERFORMANCE.md)
  - **Integral numerics caveat**: those kernels are algebraically exact but ~1 ulp off the old values (GCC FMA contraction); energies bit-identical, gradients ≤1.7e-14 Eh/Bohr
  - **Eigensolvers**: opt-in `-eigensolver native|purify|lobpcg`, `CURCUMA_EIG_TRED2=blocked` (MKL-free / GPU-portable) — [docs/SQM_EIGENSOLVE_GPU.md](SQM_EIGENSOLVE_GPU.md)
  - **Large systems**: `-large_system_mode fragments|dc|sparse` scales SCF past ~1000 atoms — [docs/SQM_LARGE_SYSTEMS.md](SQM_LARGE_SYSTEMS.md)
  - **SCF extrapolation**: `-scf_extrapolation aspc|gauss` cuts SCF iters in opt/MD (caffeine gfn2 215→90); experimental `xlbomd` = extended-Lagrangian MD — [docs/SQM_SCF_EXTRAPOLATION.md](SQM_SCF_EXTRAPOLATION.md)
  - **GPU backends** (`-gpu cuda|rocm|vulkan|auto`): all three device-resident through `-opt`/`-md`; ROCm has FP32 mixed-precision ON by default (real win), Vulkan is opt-in only (eigensolve-bound, no speedup). Detail in [docs/SQM_GPU.md](SQM_GPU.md) / [docs/SQM_ROCM.md](SQM_ROCM.md) / [docs/SQM_VULKAN.md](SQM_VULKAN.md); roadmap in [docs/SQM_GPU_ROADMAP.md](SQM_GPU_ROADMAP.md)
  - **GPU fallbacks are counted** (Sep 28, 2026, 🤖): summary after every run, `-gpu_strict true` = exit 3 at the first one (`src/core/gpu_fallback.h`, [docs/GPU_TUNING.md](GPU_TUNING.md) section 1)
  - **Multi-GPU** (Sep 2026, 🤖 machine-tested): `-gpu_device N`; batch workers (`-sp`/`-opt` multi-XYZ, ConfSearch, Hessian) leased over a device pool (`-gpu_devices`, `-gpu_workers_per_device`); 7k-atom GFN2 fits one 20 GB card (screened integrals, on-the-fly multipoles); tuning options for manual tweaking in [docs/GPU_TUNING.md](GPU_TUNING.md); multi-GPU eigensolve + density **on by default above 4000 basis functions when several GPUs are visible** (`-gpu_eigensolver_devices`/`-gpu_density_devices ... |none`; cuSOLVERMp/NCCL, cusolverMg FP64 fallback; polymer_2x 419 -> 194 s on 4 GPUs) — [docs/MULTI_GPU.md](MULTI_GPU.md)
  - **All GPU backends are runtime `dlopen` plugins** (`libcurcuma_{cuda,rocm,vulkan}.so`; CUDA Jul 2026, ROCm/Vulkan Sep 2026): the main binary has no backend symbol or `#ifdef`, CPU-only runs start in ~9 ms; `-gpu <backend>` loads the plugin lazily and warns + falls back to CPU when it is absent; `-methods` lists the plugins found. ROCm plugin build unverified (no SDK here). See [docs/GPU_PLUGIN_STARTUP.md](GPU_PLUGIN_STARTUP.md)
- ⚠️ **GFN1-xTB (Native)** - AI-implemented, machine-tested; canonical `gfn1` backend; `-opt` works; Broyden SCF default; 14/16 sqm_reference @1e-8 vs tblite (He2 + `complex` remain) — [docs/SQM_WP2_gfn1_accuracy.md](SQM_WP2_gfn1_accuracy.md); **halogen-bond correction** (B-X···A, GFN1-only) implemented Sep 2026, GMTKN55 MAD 0.041→0.00007 kcal/mol (Known Issue #26); the residual max 0.011 is the xtb-vs-tblite STO-6G 4s/4p split for Z=19-36 (`-xtb.sto6g_legacy_4sp true` reproduces xtb bit-for-bit; Known Issue #27)
  - **GFN2 GPU pseudo-diagonalisation** (Sep 29, 2026, 🤖 machine-tested, opt-in): `-scf_pseudo_diag` replaces most FP32 eigensolves by one occupied-virtual rotation + S-re-orthonormalisation (polymer_2x 4 GPUs 113 -> 89 s; MOR41/GMTKN55 on the GPU identical); `-scf_pseudo_diag_fp64` does it in FP64 too, with an exact W for the gradient - correct but slower on FP64-weak cards, meant for full-rate-FP64 GPUs (unmeasured) — [docs/GFN2_GPU_COST_PLAN.md](GFN2_GPU_COST_PLAN.md)
  - **GFN-FF, one molecule on several GPUs** (Sep 29, 2026, 🤖 machine-tested): Coulomb tiles + projected-PCG EEQ split over the devices (`-gfnff.gpu_split_devices auto|list|none`, from 1000 / 4000 atoms); each-pair-once Coulomb kernel; single fragments take projected PCG like the CPU. polymer_2x MD step 0.377 s (1 GPU, Sep 28) -> 0.076 s (4 GPUs) — [docs/GPU_TUNING.md](GPU_TUNING.md) section 3
- ⚠️ **PM3/AM1/MNDO (Native NDDO)** - AI-implemented, machine-tested; 21/21 tests vs Ulysses reference (< 4 µEh)
- ⚠️ **Native GFN-FF** - AI-implemented, machine-tested; see [docs/GFNFF_STATUS.md](GFNFF_STATUS.md); S30L vs xtb 6.6.1 (Jul 2026): 29/30 within ~0.5 kcal/mol after F1/F2/CLI/F3/ipis + Jul 8 torsion fixes (MAD 433→1.20); residuals: 23 (.CHRG/harness quirk, curcuma judged correct), 30/AB (bond pibo) — see [docs/S30L_GFNNF_VALIDATION.md](S30L_GFNNF_VALIDATION.md). **Fragment-charge placement for charged multi-fragment species changed default (Sep 2026)**: the carrier is now chosen by chemistry, not atom index (`-gfnff.frag_charge_model` default `reference`→`ensemble`) — fixes a real index bug (GMTKN55 WATER27 reaction MAD 58.6→21.4); a continuous cross-threshold window stays opt-in (`frag_charge_s_max`, default `1.0` = off) since it is a net loss for plain GFN-FF alone — see Known Issue #31, [docs/FRAG_CHARGE_MODEL.md](FRAG_CHARGE_MODEL.md). A reused calculator's energy-only calls now also see the current CN and D4 C6 (Known Issue #32).

#### External Interfaces (Production Quality, Requires Compilation)
- **TBLite Interface** - Tight-binding DFT methods (GFN1, GFN2, iPEA1) + **Solvation** (CPCM, GB, ALPB)
- **XTB Interface** - Extended tight-binding methods (GFN-FF, GFN1, GFN2)
- **Ulysses Interface** - Semi-empirical methods (PM3, PM6, AM1, MNDO, RM1, etc.) + **Solvation** (GBSA)
- **Native GFN-FF** - Curcuma's own implementation (`gfnff`) - ✅ **IMPLEMENTED**

### 2. Force Field Methods
- **Universal Force Field (UFF)** - General-purpose molecular mechanics
- **GFN-FF** (`gfnff`) - ✅ **FULLY IMPLEMENTED** - See [docs/GFNFF_STATUS.md](GFNFF_STATUS.md)
- ⚠️ **Coarse-grained beads** (`cg`, Sep 2026, AI/machine-tested) - LJ spheres/ellipsoids on the workspace engine, input via `-load_ff_json FILE` (`cg_default`, `cg_per_atom`, `pair_interactions`); analytic sphere gradient - see [docs/CLEANUP_2026_09.md](CLEANUP_2026_09.md)
- **QMDFF** - Quantum Mechanically Derived Force Fields
- **Universal Parameter Caching** - Automatic save/load for all FF methods

### 3. Solvation Models (Implicit Solvent)
- ✅ **TBLite Solvation** - CPCM, GB (Generalized Born), ALPB for GFN methods. ⚠️ Pending in the `USE_TBLITE=OFF` dev build (see [docs/SOLVATION.md](SOLVATION.md)) - unaffected: native ALPB/GBSA below.
- ✅ **Ulysses Solvation** - GBSA (Generalized Born + SA) for GFN/MNDO methods
- ⚠️ **Native GFN1/GFN2 ALPB + GBSA** (June 2026, AI/machine-tested) - self-consistent
  ALPB (`-xtb.solvent_model alpb`, P16 kernel) and GBSA (`-xtb.solvent_model gbsa`, Still kernel)
  in the native xTB SCF, matching tblite total ΔG (Born + CDS + shift; CM5 for gfn1) to
  ≤1e-8 Eh on the validation set (CPU + GPU); `-method gfn2 -xtb.solvent water -xtb.solvent_model gbsa`
  (legacy numeric codes 3/2 still accepted). CPCM native solvation still pending.
  See [docs/SQM_SOLVATION_WP.md](SQM_SOLVATION_WP.md)
- ⚠️ **Native GFN-FF ALPB** (June 2026, AI/machine-tested) - self-consistent: the Born
  reaction field couples into the EEQ solve (`A_eeq += B`), so charges polarize in the solvent.
  `-method gfnff -gfnff.solvent water -gfnff.solvent_model alpb` matches **xtb 6.7.1** (`--gfnff
  --alpb`) to **≤1e-8 Eh** (7 mol × 4 solvents); analytic gradient FD-validated. GFN-FF has
  no separate GBSA (reference uses ALPB), so `-gfnff.solvent_model gbsa` maps to ALPB. See
  [docs/SQM_SOLVATION_WP.md](SQM_SOLVATION_WP.md) WP5
- **25+ Solvents** - water, methanol, DMSO, acetone, benzene, etc.
- **Auto-Activation** - Specify `-solvent water` to enable
- **Documentation** - See [docs/SOLVATION.md](SOLVATION.md) for details

### 4. Dispersion and Non-Covalent Corrections
- **DFT-D3** - Grimme's D3 dispersion correction
- **DFT-D4** - Next-generation D4 dispersion correction
- **H4 Correction** - Hydrogen bonding and halogen bonding corrections

### 5. Geometry Optimization
- **LBFGS Optimizer** - Limited-memory Broyden-Fletcher-Goldfarb-Shanno
- **Multiple Convergence Criteria** - Energy, gradient, RMSD-based
- **Constrained Optimization** - Distance, angle, and dihedral constraints

### 6. Conformational Analysis ✅ REFACTORED 2025
- **ConfSearch** - Systematic conformational searching (unified trajectory framework); supports **dual-method** runs (`-md_method` explore + pre-opt, `-opt_method` refine + rank; both fall back to `-method`) — see [docs/CONFSEARCH_DUAL_METHOD.md](CONFSEARCH_DUAL_METHOD.md); **restartable** via `-restart` (self-contained checkpoint: bias pool + cumulative + seeds + schedule, written to CWD + BMT) — see [docs/CONFSEARCH_RESTART.md](CONFSEARCH_RESTART.md); **registry-backed since Jul 2026** (67 PARAMs) so its flags are no longer auto-routed away, and every child computation (MD / 4 opt sites / 2 ConfScan passes) shares one `ChildConfig()` carrying charge, spin, gpu and the method sub-scopes; **RMSD-MTD bias speedup (Jul 2026)**: a rigorous Gaussian-cutoff screen (`-rmsd_mtd_screen`, default ON, physics-preserving) skips far hills before the Kabsch, plus an enforced pool cap (`-rmsd_mtd_max_gaussians`) — see [docs/CONFSEARCH_MTD_SCREEN.md](CONFSEARCH_MTD_SCREEN.md)
- **ConfScan** - Conformational scanning along reaction coordinates
- **RMSD Analysis** - Structure comparison and alignment
- **Energy-based Filtering** - Automatic conformer ranking
- **Refactored Geometry Commands** - TrajectoryWriter for JSON format (Phase 5)

### 7. Molecular Dynamics
- **SimpleMD** - Basic molecular dynamics simulation
- ⚠️ **Temperature ramps / live T / thermal regions** (Jun 2026, AI/machine-tested) - `setTargetTemperature()` live setpoint; multi-stage `temp_ramp`/`temp_schedule` (`steps`/`reach` modes); per-atom-subset `temp_regions` (Berendsen/CSVR/Andersen; NH falls back to global). No-region path byte-identical to legacy. See [docs/TEMPERATURE_RAMP.md](TEMPERATURE_RAMP.md)
- ⚠️ **Step-rejecting integrator** (`-adaptive_step`, Sep 2026, AI/machine-tested) - redoes a step subdivided when it violates conservation; two self-calibrating channels, the global energy drift and (`-adaptive_step_local`, on with it) the **hottest atom over the per-atom mean**, which keeps its contrast at 7320 atoms where the global one loses it. **Off by default**, explicit `false` and `-adaptive_step_local false` bit-identical. polymer_2x 300 fs NVE dt=0.5: **+67.59 -> +0.53 Eh, 2175 -> 246 K**, with 92 of 104 rejections from the local channel alone. Measurements in [docs/MD_LARGE_SYSTEMS.md](MD_LARGE_SYSTEMS.md)
- ⚠️ **MD time step FIXED** (`ef462fcf`, Sep 2026, AI/machine-tested) - SimpleMD advanced **1.9516 fs per requested fs** (the step multiplied a velocity in sqrt(Eh/amu) as if it were A/fs), so every MD in curcuma - all methods, ConfSearch exploration included - ran with a step ~2x too large and `-MaxTime` stretched by the same factor. Energies, forces, temperature and energy conservation were unaffected and correct; only the clock was wrong, which is why nothing caught it. Now converted in `Verlet()`/`Rattle()`/`NoseHover()` via `CurcumaUnit::Constants::FS_TO_MD_TIME`; new ctest `md_time_axis` ties the reported time to a Hessian frequency. **Every MD time recorded before this commit must be multiplied by 1.9516.**
- **Large-system MD** - why a GFN-FF MD of 7320 atoms heats: **one water molecule collapses** (H-O-H 103.4 -> 3.6 deg, that molecule alone +10.96 Eh, 12 atoms holding 99 % of the kinetic energy). Not size, not truncation - the water and polymer halves each conserve, and the earlier "dt^2 truncation" reading is corrected in [docs/MD_LARGE_SYSTEMS.md](MD_LARGE_SYSTEMS.md)
- **NEB Docking** - Nudged elastic band for transition states
- **Trajectory Analysis** - Analysis of MD trajectories
- **PLUMED Metadynamics** - Enhanced sampling via PLUMED2 plugin (`-mtd` flag) — see [docs/PLUMED_HELP.md](PLUMED_HELP.md)

### 8. Analysis Tools
- **✅ Parallel Analysis** - Frame-level parallelization with CxxThreadPool (3-8x speedup, January 2026)

### 9. Output Directory System
- **🤖 BMT (Basename.Method.Timestamp)** - Default output directory for all commands — see `src/tools/CLAUDE.md`
- **`-bak` flag** - Copy specified files from BMT directory back to CWD
- **`-no_bmt`** - Disable BMT, write output to CWD (legacy behavior)
- **✅ TrajectoryWriter** - Unified output system for Human/CSV/JSON/DAT formats
- **✅ Scattering Analysis** - P(q)/S(q) with logarithmic q-spacing and automatic gnuplot visualization (2026)
- **RMSD Calculations** - Root-mean-square deviation analysis
- **Persistent Diagram** - Topological data analysis
- **Hessian Analysis** - Second derivative calculations
- **Orbital Analysis** - Molecular orbital visualization and analysis

### 10. Core Computational Libraries
- ✅ **MNDO Integrals** - Dewar-Thiel multipole expansion for semi-empirical 2e⁻ integrals, see [docs/MNDO_INTEGRALS.md](MNDO_INTEGRALS.md)


## Completed Developments (2026)

✅ **GFN-FF / GFN1 / GFN2 cleanup + speedup** (Sep 2026, AI/machine-tested) - one FF engine
(`FFWorkspace`; legacy `ForceField`/`ForceFieldThread` GFN-FF path and its per-step feed
removed), native xTB decoupled from `QMDriver`, duplicate parameter headers / old NDDO classes /
dead EEQ paths deleted (~20k lines), table-driven `MethodFactory` registry + one shared
sub-scope list, exact hot-path fixes (LAPACK scratch, HB-gradient index, EEQ views). Energies
identical to the last digit; GFN-FF SP 1.35-1.9x, GFN1 polymer 1.45x. Numbers, method and
open items in [docs/CLEANUP_2026_09.md](CLEANUP_2026_09.md)

> Older 2025 work (parameter registry, polymorphic EnergyCalculator, native GFN2/GFN1/PM3, MNDO integrals, GFN-FF full implementation, scattering, analysis parallelization, dependency gating) is in `AIChangelog.md` + git history.

✅ **`-interaction` capability** (June 2026) - supramolecular interaction energy `E(AB)−E(A)−E(B)` for the S30L host-guest set; modes: S30L A/B/AB dir (+`.CHRG`), batch vs `reference_s30l` (MAD/RMSD), explicit `-fragA/-fragB`, single-AB auto-split
✅ **GFN-FF aromatic ring torsions fixed** (June 2026) - acyclic-only pi-sp3 rules were wrongly applied to ring torsions; gated on `!in_ring`; S30L host A now bit-identical to Fortran, validation 18/18 — see [docs/GFNFF_STATUS.md](GFNFF_STATUS.md)
✅ **GFN-FF GPU HB-freeze resolved + per-frame gradient diagnostic** (June 2026) - the GPU HB-charge freeze is correct (matches CPU+Fortran); `test_gfnff_grad_traj` is the clean force metric (MD heat-exchange is not) — see [docs/GPU_GFNNF_DISCREPANCIES.md](GPU_GFNNF_DISCREPANCIES.md)

