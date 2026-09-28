# Multi-GPU gap analysis: GFN2/GFN1 and GFN-FF (inventory, Sep 28, 2026)

Inventory only - nothing here was changed in the code. Basis: `feature/multi-gpu` at `a2fdaceb`.
Method: two read-only code inventories (one per method family) with file:line evidence, then the
claims that matter most re-checked by hand, plus one test on the GPU. Evidence labels:

- **[R]** re-checked by hand in the code (Sep 28, 2026)
- **[C]** found in the code by the inventory, file:line given, not re-checked by hand
- **[M]** measured on Sep 28, 2026
- **[D]** stated only in docs; quoted with source and date

Paths below are relative to `src/core/energy_calculators/` unless they start with `src/`.

## Summary

- **One large molecule on several GPUs exists only for GFN2/GFN1 on CUDA**, and only for the
  eigensolve, the screened-pattern density and the gradient's W. Everything else stays on the
  calculation's device, which caps the gain: polymer_2x (7320 atoms) GFN2 on 1/2/4 A4500 cards
  249/208/144 s, i.e. 1.73x on 4 [D: lab journal "MD-Stabilitaet grosser Systeme", Sep 21, 2026].
- **GFN-FF has no multi-GPU path for one molecule at all.** Multi-GPU means one device per batch
  worker. The single-GPU step (182 ms on polymer_2x, Sep 27) is dominated by two O(N^2) device
  phases, EEQ ~79 ms and Coulomb ~93 ms [D: TODO.md:779-780].
- **Batch distribution works and scales** (8 copies of polymer, GFN2: 1 GPU 100 s, 4 GPUs 28.5 s)
  [D: MULTI_GPU.md:228-235, Sep 17], but several consumers do not use it and one flag combination
  makes every worker grab every GPU.
- **Two correctness defects on the GPU path that no doc mentions** (F-1 ROCm without Coulomb,
  F-2 rejected EEQ charges used anyway), plus one stale-C6 case (F-11).
- **Several degradations are silent** at default verbosity (G2-7, G2-11, F-19); `gpu_strict`
  from the plan does not exist.
- **The big open measurements**: NVLink, a per-phase profile on the target hardware, any GFN-FF
  batch timing, and the full GFN-FF MD step after the Sep 27 changes.

## 1. GFN2/GFN1 (`-method gfn2|gfn1 -gpu ...`)

### What is distributed today (CUDA)

| stage | where | evidence |
|---|---|---|
| S, H0, L, gamma, multipole integrals | calculation device | `qm_methods/xtb_native.cpp:521-615` [C] |
| potential, Fock build, occupations, charges/multipoles, Broyden mixing, SCC energy | calculation device (GFN2 resident loop), GFN1 potential on the host | `xtb_native.cpp:1041-1219`, `cuda/xtb_gpu_context.cu:3142` [C] |
| **eigensolve** | **several devices** above 4000 basis functions (FP64: full generalized solve via cuSOLVERMp/cuBLASMp; FP32: `syevd` only) | `cuda/xtb_gpu_context.cu:2935-3129` [C] |
| **density** | **several devices**, GFN2 resident loop with screened storage only | `cuda/xtb_gpu_context.cu:3209-3366` [C] |
| D4 (+ATM), repulsion, halogen/solvation | calculation device or host | `xtb_native.cpp:1741-1842, 2730-2803` [C] |
| gradient | calculation device and host; **only W** is split | `cuda/xtb_gpu_context.cu:5347-5556` [C] |

The part that stays on one device was ~98 of 194 s on 4 A4500 [D: GPU_TUNING.md:89-95, Sep 17];
on the optimised geometry a constant ~30 s against 249 s single-GPU [D: lab journal, Sep 21].

### Gaps

**Correctness / wrong settings**

- **G2-13 ROCm forces mixed precision on** *(open, documentation only - operator decision Sep 28)* and overrides a user's `-scf_mixed_precision false`:
  `qm_methods/xtb_hip_method.cpp:109` calls `setMixedPrecision(true)` unconditionally [R]. The CUDA
  wrapper had the same bug and fixed it for itself (`xtb_gpu_method.cpp:695-705`) [C].
- **G2-17 (suspected, refuted on one molecule)**: GFN1 with screened storage might read an unfilled
  pattern density in the gradient. Caffeine, `-gpu_sparse_integrals on` (screened storage active,
  100 % of pairs), gradient identical to the dense GPU run; both differ from the CPU by 1.9e-6
  Eh/A (loose default SCF threshold) [M]. One molecule only.

**Device selection and batch coverage**

- **G2-2 The split ignores `-gpu_devices` and `-gpu_device`.** Default and `all` mean every visible
  device, device 0 included (`xtb_gpu_method.cpp:556-570`) [R]. A large single run on a shared node
  places ranks on busy GPUs; workaround: an explicit list or `CUDA_VISIBLE_DEVICES`.
- **G2-3 An explicit `-gpu_eigensolver_devices`/`-gpu_density_devices` inside a batch worker
  bypasses the lease guard** (`xtb_gpu_method.cpp:556-558, 607`) [R] and `all` then expands to every
  device - every worker builds NCCL communicators over all GPUs. GPU_TUNING.md:40 says "batch
  workers keep it off", true only for the default. There is no API to lease k devices for one worker.
- **G2-16 Batch coverage is partial** [C]: `-opt` multi-XYZ at `-threads <= 1` runs sequentially
  (`src/capabilities/optimizer_factory.cpp:407-417`); ConfScan's energy recomputation is sequential
  (`src/capabilities/confscan.cpp:545, 608`); the Hessian defaults to 1 thread, so 1 GPU
  (`src/capabilities/hessian.cpp:283`); only the `-sp` batch raises its worker count to the GPU slots
  (`src/main.cpp:1757-1760`). The pool balances by lease count, not free memory or foreign load
  (`src/core/gpu_device_pool.cpp:164-185`).
- **G2-14 ConfSearch forwards only `gpu`**, not `gpu_device`, `gpu_sparse_integrals`,
  `gpu_memory_check` or `gpu_eigensolver_*` (`src/capabilities/confsearch.cpp:1078-1097`) [C].
- **G2-12 Vulkan counts devices it cannot use** (`vulkan/vk_context.cpp:209-237`); the pool can
  lease them and that worker runs on the CPU [C].

**Scaling of one molecule**

- **G2-4 Integrals, Fock, potential, mixing, D4 and most of the gradient are not distributed** -
  the ceiling above. Distributing them (column-block ownership) is the planned next step; whether it
  pays depends on the per-phase profile on the target hardware, which does not exist (X-2).
- **G2-9 Distributed buffers and the distributed L are freed after every SCF**
  (`cuda/xtb_gpu_context.cu:5079`, `xtb_distributed_eigensolver.cpp:370-391`) [C], so every MD/opt
  step re-allocates and re-scatters them. Cost unmeasured.
- **G2-10 Dense host downloads every geometry** (S, H0, L, gamma; P and C unless the deferred path
  is active; `xtb_native.cpp:543-550, 1674-1679`) [C] - 1.9 GB per n^2 matrix at nao 15444.
- **G2-5 Device 0 keeps every dense n^2 matrix**; adding GPUs only moves the eigensolver workspace
  off it: 17.9 -> 13.5 GB on 4 GPUs, 12.3 GB with device 0 out of the solver [D: GPU_TUNING.md:72-76,
  Sep 17]. The memory estimate is single-device (`cuda/xtb_gpu_context.cu:3517-3552`), so a system
  that would fit only when split can be refused. `CudaBuffer` uses `int` sizes
  (`ff_methods/cuda/gfnff_soa.h:48`): dense n^2 buffers stop at nao ~46340 whatever the GPU count [C].
- **G2-6 The distributed density exists only for CUDA + GFN2 resident loop + screened storage**
  [C]. GFN1 gets none and prints no status; with dense storage the status wrongly says "nao below
  min" (`cuda/xtb_gpu_context.cu:3643-3644` vs `3212`).

**Silent or misleading behaviour**

- **G2-7** [C]: cusolverMg chosen instead of cuSOLVERMp is info-level only
  (`xtb_gpu_method.cpp:441-442`); without cuBLASMp the FP64 generalized path quietly becomes a
  distributed `syevd` (`xtb_distributed_eigensolver.cpp:155-166`); a failed verification-buffer
  allocation is not a warning (`cuda/xtb_gpu_context.cu:403-425`); `-scf_gpu_partial_diag` disables
  the split without a message.
- **G2-8** [C]: Mg FP32 is not excluded up front - only the per-solve check rejects it; with
  `-gpu_eigensolver_verify false` it is accepted unchecked (`cuda/xtb_gpu_context.cu:394, 646-666`),
  while GPU_TUNING.md:44/99 and MULTI_GPU.md:66 say "FP64 only".
- **G2-11** [C]: batch workers run at verbosity 0 (`src/main.cpp:1766, 1783`), so a CPU fallback
  of a worker is invisible; the printed "device N" is the lease, not proof the GPU ran.
- **G2-15** [C]: `gpu_strict`, `gpu_multi_mode`, `gpu_profile`, `gpu_count` from the multi-GPU plan
  are not implemented (no hits in `src/`).

**Backend parity**

- **G2-1 ROCm and Vulkan cannot split one molecule**: `setDistributed*` is called only from
  `qm_methods/xtb_gpu_method.cpp:597` [R]. ROCm/Vulkan device selection is compile-unverified here
  [D: MULTI_GPU.md:80].

## 2. GFN-FF (`-method gfnff -gpu ...`)

### What runs where in one MD energy+gradient call (CUDA)

| stage | where | evidence |
|---|---|---|
| topology, hybridisation, pi, fragments | host at setup; full two-pass recomputation can fire during MD when the 0.5 Bohr displacement flag is set | `ff_methods/gfnff_method.cpp:1053-1180` [C] |
| CN, dCN pair list | device (atomics); dCN list rebuilt every gradient step while the flag stays set | `ff_methods/cuda/ff_workspace_gpu.cu:2930-3034`, `qm_methods/gfnff_gpu_method_impl.h:739-771` [C] |
| Phase-2 EEQ | device: nfrag=1 dense Cholesky every step; nfrag>1 and N>=500 projected PCG (WP7-E) on a dense N x N matrix | `qm_methods/gfnff_gpu_method_impl.h:1119` [R], `ff_methods/cuda/eeq_solver_gpu.cu:1925-2142` [C] |
| D4 pair list | host rebuild (60/50 Bohr skin) + re-upload; device build opt-in | `ff_methods/gfnff_method.cpp:1963-2019` [C] |
| repulsion list | host cell list, rebuilt when the 2 Bohr skin is exceeded | `ff_methods/gfnff_method.cpp:1880-1944` [C] |
| HB/XB lists | host, forced every 10 calls | `ff_methods/gfnff_method.cpp:1799-1873` [C] |
| bonded, repulsion, BATM, HB/XB terms | device, 3 streams | `ff_methods/cuda/ff_workspace_gpu.cu:1995-2314` [C] |
| Coulomb | device, implicit all-pairs gather | `ff_methods/cuda/gfnff_kernels.cu:382-428` [C] |
| reduction, output | global `atomicAdd`; D2H of energies, gradient, dEdCN, charges | `ff_methods/cuda/gfnff_kernels.cu:30-36` [C] |

Measured: 182 ms per energy call on polymer_2x on one A4500, EEQ ~79 ms, Coulomb phase ~93 ms
[D: TODO.md:779-780, GPU_TUNING.md:162, Sep 27]. Host list rebuilds, when they fire: repulsion
~240 ms, HB/XB ~100 ms, D4 host ~670 ms at 7320 atoms [D: GFNFF_PAIR_LIST_REFRESH.md:58, 98, 159].

### Gaps

**Correctness**

- **F-1 ROCm computes no Coulomb term by default - since `ab6e3f5e` (Sep 17, 2026).** *Open,
  documentation only (operator decision Sep 28: no ROCm code changes without ROCm hardware); workaround
  in TODO.md and CLAUDE.md Known Issue #35.* The shared
  wrapper sets implicit Coulomb pairs by default (`qm_methods/gfnff_gpu_method_impl.h:403-406`), which
  clears the host Coulomb list (`ff_methods/gfnff_method.cpp:3757-3761`). The HIP workspace has no
  implicit branch: it takes the self-energy parameters only from a non-empty pair list
  (`ff_methods/rocm/gfnff_rocm.hip:4948-4975`) and launches Coulomb only if `coulomb.n > 0`
  (`:6337`) [R]. Not run (no ROCm SDK here). Workaround `-gfnff.gpu_coulomb_implicit false`; the
  per-atom self-energy fix of CLAUDE.md Known Issue #8 is not mirrored on ROCm either [C].
- **F-2 A rejected GPU EEQ solution is still used.** *FIXED Sep 28, 2026 (stage A): the D2D copy now
  follows the validation, `setEEQCharges` clears the pending-device flag, and the log names what is used.
  Measured with the new test hook `CURCUMA_EEQ_GPU_FORCE_REJECT=1` (triose NVE MD, frame at 5 fs): before,
  Epot = -9.873827 = the CPU Phase-2 single point at that geometry, i.e. the rejected solution; after,
  -9.873197, neither Phase 1 (-9.846663) nor Phase 2 (-9.873546). Correction to the text below: the
  fallback keeps the charges of the last ACCEPTED solve (as the CPU does, `gfnff_method.cpp`, "keeping
  previous m_charges"), not the Phase-1 topology charges the log used to claim.* [M]
  Original finding: In the single-fragment GPU path the charges are
  copied device-to-device (`qm_methods/gfnff_gpu_method_impl.h:1328`) **before** the NaN/|q|>50 check
  (`:1330-1354`). The fallback only sets the host vector (`:1477`); `setEEQCharges` does not clear
  `m_device_charges_ready` (`ff_methods/cuda/ff_workspace_gpu.cu:1255-1266`), so the charge upload is
  skipped (`:2429-2433`) and Coulomb runs on the rejected charges while the log says "Phase 1
  topology charges" [R]. Code reading, not reproduced.
- **F-11 Energy-only GPU calls use stale C6** - *FIXED Sep 28, 2026 (stage A): in energy mode the C6
  refresh is enqueued before Phase 1 from `d_cn_final`. New ctest `gfnff_gpu_energy_only_geometry_*`
  (energy-only call at a displaced geometry vs a gradient call there): triose 2.1e-4 -> 2e-15 Eh, caffeine
  3.0e-5 -> 0 [M].* Original finding: the refresh is gated on `gradient`
  (`qm_methods/gfnff_gpu_method_impl.h:909`) [R]; the CPU refreshes in energy-only calls too
  (`ff_methods/gfnff_method.cpp:1633-1644`). Affects line searches and energy-only scans.
- **F-12 Possible unordered cross-stream read** - *addressed Sep 28, 2026 (stage A): WP7-E now orders
  its stream after an event recorded on the legacy default stream (all curcuma streams are blocking
  streams). Never observed to fail, so there is no before/after number; polymer_2x GPU MD energies are
  unchanged.* Original finding: `d_rhs_atoms` is written on the workspace stream
  (`ff_methods/cuda/ff_workspace_gpu.cu:1957-1966`) and read by WP7-E on the solver stream
  (`ff_methods/cuda/eeq_solver_gpu.cu:1966`) without an event wait; the sync the comment relies on was
  removed by WP5-D [C]. Not observed to fail.
- **F-13 The CPU EEQ inside a CUDA process segfaults on polymer_2x** - *the CUDA wrapper no longer
  reads `eeq_rocm_cpu_fragment_threshold` (stage A, Sep 28); the segfault itself (F-Q9) stays open.* (TECHNICAL_DEBT.md F-Q9,
  Sep 24, unresolved). The default CUDA path avoids it, but `eeq_rocm_cpu_fragment_threshold`
  (help: "ROCm only") is also read on CUDA (`qm_methods/gfnff_gpu_method_impl.h:552`) [C].

- **F-20 (found Sep 28, 2026 while testing F-11) - the CPU had the same defect class in the Coulomb
  term.** The self-energy uses `chi = chi_base + cnf*sqrt(CN)` with the workspace CN
  (`ff_methods/ff_workspace.cpp:522`), which only the gradient path set (`setCNDerivatives`). An
  energy-only call after a geometry change therefore used the CN of the last gradient call, or the
  static chi of the start geometry: triose displaced by up to 0.08 A, 4.6e-3 Eh, caffeine 2.8e-4 Eh,
  all in the Coulomb term. Affects CPU line searches, energy-only scans and the energy-FD Hessian.
  *FIXED (stage A): `FFWorkspace::setCN()` in the energy-only path; ctest `gfnff_energy_only_geometry_*`
  now exact; MOR41 (95) and GMTKN55 (2462) gfnff single points bit-identical to the previous binary.* [M]

**Scaling of the single-GPU step (prerequisite for any multi-GPU gain)**

- **F-4 A single fragment always takes the dense O(N^3) GPU Cholesky every step**
  (`qm_methods/gfnff_gpu_method_impl.h:1119`) [R]; the CPU switches to projected PCG from 500 atoms
  (`ff_methods/eeq_solver.cpp:1583-1586`). Contradicts GPU_TUNING.md:161 ("select the solver alike").
- **F-5 EEQ is O(N^2) per step**: the dense N x N matrix is rebuilt each step, every PCG iteration is a
  dense `dsymv` plus three blocking host syncs (`ff_methods/cuda/eeq_solver_gpu.cu:1994-2034, 2084`)
  [C]. The only open alternative is a matrix-free/FMM matvec.
- **F-6 Coulomb evaluates all N^2 pairs, each twice** (`ff_methods/cuda/gfnff_kernels.cu:401-421`) [C].
- **F-7 List maintenance is serial on the host** and re-uploads with blocking `cudaMemcpy`
  (`ff_methods/cuda/gfnff_soa.h:85-92`); the CPU `FFWorkspace` is kept alive and updated on each
  rebuild [C]. The device D4 build exists but is opt-in (630 -> 32 ms) [D: GPU_TUNING.md:158].
- **F-8 The two-pass topology recomputation can fire during MD** [C]; its cost has never been measured.
- **F-9 Per-step overhead**: a blocking displacement sync, three CN downloads nobody reads, WP5-C
  skip-check kernels whose flag has no reader (`ff_methods/cuda/ff_workspace_gpu.cu:3000-3016, 2719`) [C].
- **F-10 CUDA graphs are permanently off** (graph capture failed on Blackwell;
  `ff_methods/cuda/ff_workspace_gpu.cu:1036-1042`) [R]; the comment at `:2394` claims the opposite.
- **F-15 Size wall at N >= 46341**: `int` products `N*N` overflow (`ff_methods/cuda/eeq_solver_gpu.cu:273, 357`,
  `ff_methods/cuda/ff_workspace_gpu.cu:1543, 1637`, `ff_methods/cuda/gpu_utils.cpp:62`); the memory
  estimate ignores the pair lists [C].
- **F-18** FP32 kernels exist but are off; CUDA still accumulates the gradient with per-pair atomics
  where ROCm has gather kernels [C, D: SQM_GFNFF_GPU_WP.md].

**Multi-GPU**

- **F-3 No path spreads one GFN-FF molecule over several GPUs** - no NCCL, peer access or
  multi-device logic in the GFN-FF GPU code [C]. A standalone `-md` uses one device
  (`src/capabilities/simplemd.cpp:852-858`); batch consumers (ConfSearch MD/opt, `-sp`/`-opt`
  batch, CurcumaOpt threads, Hessian) lease one device per worker [C]. A split would need at each
  step: coordinate broadcast, allreduce of CN and dEdCN, a row-block distributed EEQ matrix with an
  allreduce per PCG dot product, a charge broadcast for Coulomb and a gradient reduction - with the
  host stages F-7/F-8 still serial.
- **F-16 Batch workers leak** the host GFNFF object and the full parameter set on destruction
  (`qm_methods/gfnff_gpu_method_impl.h:368-370, 500-505`); the comment says "~100 KB", but the set
  holds every pair list [C]. Size unmeasured.
- No GFN-FF batch timing exists; all pool numbers are GFN2 [D: MULTI_GPU.md:226-237].

**Determinism and visibility**

- **F-14** GPU results are not run-to-run reproducible (`atomicAdd` gradient, energy, CN, fragment
  sums, pair order) [C]; known, deferred by the operator [D: TODO.md:772-781, Sep 27].
- **F-19** "All EEQ paths exhausted -> Phase-1 charges" warns only at verbosity >= 1
  (`qm_methods/gfnff_gpu_method_impl.h:1474-1479`), so MD children at verbosity 0 never show it [C].

**ROCm parity**

- **F-17** Every device EEQ variant is a stub on ROCm (GPU-Schur, WP7-A/C/E, batched;
  `ff_methods/rocm/eeq_solver_hip.hiph:406-445`); below 16 fragments a dense rocSOLVER solve with host
  Schur step, above it a host EEQ. Energy reductions hard-code wave32 (`ff_methods/rocm/gfnff_rocm.hip:86-128`),
  wrong on CDNA/MI (wave64). The September mirrors are uncompiled [C].

## 3. Cross-cutting

- **X-1 GPU/EEQ tuning parameters missing from the registry** [M]: the build reports 26
  "Malformed PARAM" warnings (18 in `ff_methods/eeq_solver.h`, 8 in `ff_methods/gfnff.h`) - multi-line
  PARAM macros that the extractor drops. Among the dropped names: `max_pcg_iterations`,
  `pcg_tolerance`, `pcg_large_threshold`, `gpu_block_jacobi_max_frag_atoms`,
  `gpu_block_jacobi_max_nfrag`, `eeq_pcg_nfrag_threshold`, `eeq_contact_prefer_exact`,
  `eeq_distance_cutoff`, `eeq_extrapolation`, `solvent_model`; `solve_method` has no registry entry
  either. The values still reach the code in dotted form, but they are missing from `-help`, flat-flag
  routing and `-export_run`. Known trap (single-line PARAM rule), not yet cleaned up.
- **X-2 Missing measurements**: no NVLink run; no per-phase GPU profile on the target (H200)
  hardware [D: MULTI_GPU.md:355-357]; no GFN-FF batch timing; the full GFN-FF MD step was last timed
  on Sep 21 (1.475 s/step with one A4500 vs 1.25 s on 16 CPU threads, MD_LARGE_SYSTEMS.md:117-121),
  before WP7-E and the skin default - the Sep 27 figure is the energy call only; the cost of G2-9
  and F-8 inside MD is unknown.
- **X-3 Doc drift** [C]: MD_LARGE_SYSTEMS.md:117-131 ("GPU buys nothing") predates the Sep 27 changes;
  MULTI_GPU.md:28 (CPU/GPU 5.7e-7 Eh) is superseded by CLAUDE.md #33; TODO.md:1011-1018 says ConfSearch
  switches the GPU off at `-threads > 1` - true only while the device pool is inactive (one
  visible GPU); with an active pool the GPU stays on (`src/capabilities/confsearch.cpp:161-170`) [R]; GPU_TUNING.md:36 ("off by default") contradicts :40 and
  the code; GPU_TUNING.md:161 contradicts F-4; GFNFF_PAIR_LIST_REFRESH.md gives the repulsion rebuild
  as 52 ms (:116) and ~240 ms (:159); comments in `qm_methods/gfnff_gpu_method_impl.h:34-39, 959-965`
  describe the old CPU-EEQ pipeline.

## 4. What to tackle, in order

Ordered by risk first, then by what limits scaling. No numbers here are promises; each item needs
its own measurement.

1. **Correctness on the GPU path**: F-1 (ROCm Coulomb), F-2 (rejected charges), F-11 (stale C6 in
   energy-only calls), G2-13 (ROCm overrides the user's precision), F-13 (CPU EEQ reachable in a CUDA
   process). Each is small and local.
2. **Visibility and device ownership**: a `gpu_strict` / loud-fallback mode (G2-7, G2-11, F-19,
   G2-15); the split respects `-gpu_devices`/`-gpu_device` and the lease (G2-2, G2-3); ConfSearch
   forwards the GPU keys (G2-14); Vulkan filters unusable devices (G2-12); fix the dropped PARAMs (X-1).
3. **Measure before building** (X-2): per-phase profile of GFN2 and of one GFN-FF MD step on the
   target hardware, one NVLink run, one GFN-FF batch run. This decides items 4 and 5.
4. **GFN2 one-molecule scaling**: keep distributed buffers and L across MD/opt steps (G2-9);
   multi-GPU-aware memory estimate (G2-5); then, if the profile says so, distribute integrals and the
   Fock build (G2-4).
5. **GFN-FF single-GPU step** (before any GFN-FF multi-GPU work, since the host stages cap it):
   projected PCG on the GPU for single fragments (F-4); a matrix-free EEQ/Coulomb (F-5, F-6); list
   maintenance on the device, starting by making the existing device D4 build the default (F-7);
   remove the dead per-step syncs (F-9).
6. **GFN-FF multi-GPU** (F-3), only if 5 leaves the O(N^2) phases dominant; the coupling points are
   listed under F-3.
7. **ROCm parity** (F-17, wave64) and the GFN1 distributed density (G2-6).

Human review pending; this is an AI-generated inventory. Evidence marked [C] was not re-checked by hand.
