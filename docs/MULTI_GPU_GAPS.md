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

- **G2-2 The split ignores `-gpu_devices` and `-gpu_device`.** *FIXED Sep 28, 2026 (stage B): `all`/default
  means the `-gpu_devices` list, a single calculation with a pool runs on the first listed device
  (EnergyCalculator, key `gpu_device_auto`), an explicit `-gpu_device N` keeps the default split off.
  Measured by sampling `nvidia-smi` during a polymer GFN2 run with `-gpu_devices 2,3` (split thresholds
  lowered to 1000): old binary on devices 0,1,2,3, new on 2,3 only; energy identical (-2088.25340678).* [M]
  Original finding: Default and `all` mean every visible
  device, device 0 included (`xtb_gpu_method.cpp:556-570`) [R]. A large single run on a shared node
  places ranks on busy GPUs; workaround: an explicit list or `CUDA_VISIBLE_DEVICES`.
- **G2-3 An explicit `-gpu_eigensolver_devices`/`-gpu_density_devices` inside a batch worker
  bypasses the lease guard** *- FIXED Sep 28 (stage B): ignored in a leased worker with a warning and a
  counted fallback (4-structure `-sp` batch with `-gpu_eigensolver_devices all`: summary "4 x ...
  ignored in a batch worker"). A k-device lease is still not implemented.* [M] (`xtb_gpu_method.cpp:556-558, 607`) [R] and `all` then expands to every
  device - every worker builds NCCL communicators over all GPUs. GPU_TUNING.md:40 says "batch
  workers keep it off", true only for the default. There is no API to lease k devices for one worker.
- **G2-16 Batch coverage is partial** *- partly FIXED Sep 28 (stage B): the `-opt` multi-XYZ batch and the
  threaded Hessian run at least one worker per GPU slot (4 caffeine, gfn2: opt 5.28 -> 2.40 s; Hessian
  1 -> 4 devices, 4.72 -> 2.49 s). ConfScan's recomputation stays sequential. Found on the way (H-1,
  fixed): with more than one Hessian worker the frequencies were computed but never printed, because
  the workers' logger save/restore interleaved and left the level at 0 (CLAUDE.md Known Issue #3,
  pre-existing, also on the CPU with `-threads 4`). Measured side effects, not fixed: the FD Hessian
  depends on the worker split at the default SCF threshold (GPU 1 vs 4 devices: max 0.3 cm-1; CPU 1 vs 4
  threads: 0.1 = print precision) and CPU vs GPU differ by up to 4.5 cm-1 (pre-existing, same with the
  old binary); `-scf_threshold`/`-hessian.scf_threshold` did NOT reach the Hessian workers (H-2: `executeHessian`
  forwarded only `controller["hessian"]`, and ConfigManager drops object values). FIXED Sep 29: scopes are
  forwarded and re-attached; with `-scf_threshold 1e-9` CPU, 1 GPU and 4 GPUs give identical frequencies,
  i.e. the 4.5 cm-1 CPU/GPU gap was loose-SCF noise.* [M] Original finding [C]: `-opt` multi-XYZ at `-threads <= 1` runs sequentially
  (`src/capabilities/optimizer_factory.cpp:407-417`); ConfScan's energy recomputation is sequential
  (`src/capabilities/confscan.cpp:545, 608`); the Hessian defaults to 1 thread, so 1 GPU
  (`src/capabilities/hessian.cpp:283`); only the `-sp` batch raises its worker count to the GPU slots
  (`src/main.cpp:1757-1760`). The pool balances by lease count, not free memory or foreign load
  (`src/core/gpu_device_pool.cpp:164-185`).
- **G2-14 ConfSearch forwards only `gpu`** *- FIXED Sep 28 (stage B): `ChildConfig()` forwards every global
  GPU key present in the controller (code change; not exercised by a run)*, not `gpu_device`, `gpu_sparse_integrals`,
  `gpu_memory_check` or `gpu_eigensolver_*` (`src/capabilities/confsearch.cpp:1078-1097`) [C].
- **G2-12 Vulkan counts devices it cannot use** *(open - the Vulkan plugin is not built here, so a filter
  could not even be compiled; documented only)* (`vulkan/vk_context.cpp:209-237`); the pool can
  lease them and that worker runs on the CPU [C].

**Scaling of one molecule**

- **G2-4 Integrals, Fock, potential, mixing, D4 and most of the gradient are not distributed** -
  the ceiling above. Distributing them (column-block ownership) is the planned next step; whether it
  pays depends on the per-phase profile on the target hardware, which does not exist (X-2).
  **Part 1 (Sep 29, 2026)**: the dense S/H0 build skips atom pairs outside the screening cutoff
  (`k_overlap_h0` 11.5 -> 1.3 s at polymer_2x, on every GPU count - a single-device fix made a split
  unnecessary) and the Cholesky of S runs as cuSOLVERMp `potrf` on the eigensolver's devices
  (3.5 -> 1.8 s on 4 GPUs); see MULTI_GPU.md "Setup of one large molecule".
- **G2-9 Distributed buffers and the distributed L are freed after every SCF**
  (`cuda/xtb_gpu_context.cu:5079`, `xtb_distributed_eigensolver.cpp:370-391`) [C], so every MD/opt
  step re-allocates and re-scatters them. Cost unmeasured.
- **G2-10 Dense host downloads every geometry** (S, H0, L, gamma; P and C unless the deferred path
  is active; `xtb_native.cpp:543-550, 1674-1679`) [C] - 1.9 GB per n^2 matrix at nao 15444.
  Measured Sep 29, 2026 (`-verbosity 3` setup split): 3.7 s per geometry at polymer_2x, the same on
  1/2/4 GPUs - now the largest part of the setup after the device build (5.0 s).
- **G2-5 Device 0 keeps every dense n^2 matrix**; adding GPUs only moves the eigensolver workspace
  off it: 17.9 -> 13.5 GB on 4 GPUs, 12.3 GB with device 0 out of the solver [D: GPU_TUNING.md:72-76,
  Sep 17]. The memory estimate is single-device (`cuda/xtb_gpu_context.cu:3517-3552`), so a system
  that would fit only when split can be refused. `CudaBuffer` uses `int` sizes
  (`ff_methods/cuda/gfnff_soa.h:48`): dense n^2 buffers stop at nao ~46340 whatever the GPU count [C].
- **G2-6 The distributed density exists only for CUDA + GFN2 resident loop + screened storage**
  [C]. *Status message FIXED Sep 28 (stage B): it now names the real cause (dense storage, nao below the
  threshold, or "never reached: GFN1 has no such path"); the GFN1 path itself is still missing.* GFN1 gets none and prints no status; with dense storage the status wrongly says "nao below
  min" (`cuda/xtb_gpu_context.cu:3643-3644` vs `3212`).

**Silent or misleading behaviour**

- **G2-7** *- FIXED Sep 28 (stage B): every such status is a warning and a counted fallback (see G2-15);
  the partial-diag and no-cuBLASMp cases got their own status (`distributedEigensolverDegradation`).
  Not triggered on this machine (cuSOLVERMp + cuBLASMp present), so exercised only by code reading.* [C]: cusolverMg chosen instead of cuSOLVERMp is info-level only
  (`xtb_gpu_method.cpp:441-442`); without cuBLASMp the FP64 generalized path quietly becomes a
  distributed `syevd` (`xtb_distributed_eigensolver.cpp:155-166`); a failed verification-buffer
  allocation is not a warning (`cuda/xtb_gpu_context.cu:403-425`); `-scf_gpu_partial_diag` disables
  the split without a message.
- **G2-8** *- FIXED Sep 28 (stage B): cusolverMg never takes an FP32 solve (`distReady`).* [C]: Mg FP32 is not excluded up front - only the per-solve check rejects it; with
  `-gpu_eigensolver_verify false` it is accepted unchecked (`cuda/xtb_gpu_context.cu:394, 646-666`),
  while GPU_TUNING.md:44/99 and MULTI_GPU.md:66 say "FP64 only".
- **G2-11** *- FIXED Sep 28 (stage B) by the fallback summary (G2-15).* [C]: batch workers run at verbosity 0 (`src/main.cpp:1766, 1783`), so a CPU fallback
  of a worker is invisible; the printed "device N" is the lease, not proof the GPU ran.
- **G2-15** [C]: `gpu_strict`, `gpu_multi_mode`, `gpu_profile`, `gpu_count` from the multi-GPU plan
  are not implemented (no hits in `src/`). *`gpu_strict` IMPLEMENTED Sep 28 (stage B):
  `src/core/gpu_fallback.{h,cpp}` counts every GPU degradation (worker CPU fallback, plugin missing or
  declined, device init failure, eigensolver degradations, GFN-FF EEQ fallbacks, split flags ignored in
  a worker); `main` prints a summary after the run at warning level whatever the verbosity; `-gpu_strict
  true` ends the run at the first one with exit code 3. Measured with the EEQ test hook: triose MD at
  `-verbosity 0` prints "4 x GFN-FF EEQ: no valid device solution", `-gpu_strict true` exits 3, a clean
  run with `-gpu_strict true` exits 0 with no report. `gpu_strict` is a global key like the other gpu_*
  keys (main.cpp `global_params`), not a registry PARAM. The other three names stay unimplemented.* [M]

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

- **F-4 A single fragment always takes the dense O(N^3) GPU Cholesky every step** *- FIXED Sep 29, 2026:
  one fragment takes projected PCG from eeq_ppcg_min_atoms on, as on the CPU (WP5-A stays the fallback;
  the per-fragment dense block is no longer allocated for nfrag = 1). polymer (1410, one fragment): EEQ
  8.6 -> 7.6-8.1 ms, energies identical; larger single fragments gain more (O(N^2 it) vs O(N^3)), not
  measured - no such test system here.* [M]
  (`qm_methods/gfnff_gpu_method_impl.h:1119`) [R]; the CPU switches to projected PCG from 500 atoms
  (`ff_methods/eeq_solver.cpp:1583-1586`). Contradicts GPU_TUNING.md:161 ("select the solver alike").
- **F-5 EEQ is O(N^2) per step**: the dense N x N matrix is rebuilt each step, every PCG iteration is a
  dense `dsymv` plus three blocking host syncs (`ff_methods/cuda/eeq_solver_gpu.cu:1994-2034, 2084`)
  [C]. The only open alternative is a matrix-free/FMM matvec.
- **F-6 Coulomb evaluates all N^2 pairs, each twice** *- FIXED Sep 29, 2026: `k_coulomb_tiles`, each pair
  once, j tile in shared memory: polymer_2x 93.7 -> 31.6 ms, polymer 17.4 -> 2.7 ms; gradient equal to 9e-16,
  energy to 3e-13 Eh.* [M] Original: (`ff_methods/cuda/gfnff_kernels.cu:401-421`) [C].
- **F-7 List maintenance is serial on the host** and re-uploads with blocking `cudaMemcpy`
  (`ff_methods/cuda/gfnff_soa.h:85-92`); the CPU `FFWorkspace` is kept alive and updated on each
  rebuild [C]. The device D4 build exists but is opt-in (630 -> 32 ms) [D: GPU_TUNING.md:158].
- **F-8 The two-pass topology recomputation can fire during MD** [C]. *Measured and FIXED Sep 29, 2026:
  with `topology_mode auto` every forced HB/XB refresh (every 10 calls) re-ran the full topology (CPU EEQ
  N=7320, Hueckel, BATM; ~1.65 s per firing, ~160 ms/step on polymer_2x) because an atom had moved > 0.5
  Bohr, and every repulsion rebuild re-derived the bond list (26.8 M pair tests, ~75 % of ~265 ms). The
  reference never recomputes the topology in MD (gfnff_hbset reads only the setup topology). Default is
  now `topology_mode constant`, and the bond list is frozen with it. MD results unchanged: triose CPU 300
  steps bit-identical auto vs constant (energy terms, charges, gradients), polymer_2x GPU 45 steps <= 4e-12
  Eh (GPU run-to-run noise). polymer_2x MD step: GPU 0.377 -> 0.200 s, CPU (16 thr) 1.36 -> 1.08 s.* [M]
- **F-9 Per-step overhead**: a blocking displacement sync, three CN downloads nobody reads, WP5-C
  skip-check kernels whose flag has no reader (`ff_methods/cuda/ff_workspace_gpu.cu:3000-3016, 2719`) [C].
- **F-10 CUDA graphs are permanently off** (graph capture failed on Blackwell;
  `ff_methods/cuda/ff_workspace_gpu.cu:1036-1042`) [R]. *Correction Sep 28, 2026: the comment near
  `:2394` does not claim the opposite - it is conditional ("if Phase-1 ran as a graph"); it now says
  that this never happens today.*
- **F-15 Size wall at N >= 46341**: `int` products `N*N` overflow (`ff_methods/cuda/eeq_solver_gpu.cu:273, 357`,
  `ff_methods/cuda/ff_workspace_gpu.cu:1543, 1637`, `ff_methods/cuda/gpu_utils.cpp:62`); the memory
  estimate ignores the pair lists [C].
- **F-18** FP32 kernels exist but are off; CUDA still accumulates the gradient with per-pair atomics
  where ROCm has gather kernels [C, D: SQM_GFNFF_GPU_WP.md].

**Multi-GPU**

- **F-3 No path spreads one GFN-FF molecule over several GPUs** *- IMPLEMENTED Sep 29, 2026 (CUDA):
  Coulomb tiles + projected-PCG EEQ (column blocks + matvec) over the devices, peer copies, ctest
  `gfnff_gpu_split_equals_single`; polymer_2x MD step 1/2/4 GPUs 0.139/0.099/0.076 s. Numbers and options:
  GPU_TUNING.md section 3. Remaining EEQ cost at 4 GPUs is ~70 PCG iterations x ~0.5 ms. The loop now
  syncs once per iteration instead of three times (alpha/beta on the device; 4 GPUs 0.076 -> 0.073 s/step,
  charges within 3e-11 e of the exact CPU solve). Tried and NOT kept (measured, polymer_2x): an exact block
  inverse per small fragment (1500 waters) as preconditioner - iterations 75/68/68 -> 73/66/68; a linear
  extrapolation of the start vector - 1-2 of ~68 iterations. The conditioning of the inter-fragment
  Coulomb coupling sets the iteration count; a global preconditioner (deflation, multigrid) would be the
  next lever. `CURCUMA_PPCG_ITERS=1` prints the iterations of every solve.* [M] Original: - no NCCL, peer access or
  multi-device logic in the GFN-FF GPU code [C]. A standalone `-md` uses one device
  (`src/capabilities/simplemd.cpp:852-858`); batch consumers (ConfSearch MD/opt, `-sp`/`-opt`
  batch, CurcumaOpt threads, Hessian) lease one device per worker [C]. A split would need at each
  step: coordinate broadcast, allreduce of CN and dEdCN, a row-block distributed EEQ matrix with an
  allreduce per PCG dot product, a charge broadcast for Coulomb and a gradient reduction - with the
  host stages F-7/F-8 still serial.
- **F-16 Batch workers leak** *- FIXED Sep 29, 2026. Root cause was NOT CUDA: the nvcc-compiled .cu TUs
  of the CUDA plugin lacked `-DEIGEN_MAX_ALIGN_BYTES=64` (hipcc got it in June), so Eigen's inline
  allocation code existed in two ABIs inside the plugin; freeing the host GFNFF object aborted with
  `free(): invalid pointer` (reproduced deterministically, gdb: `GFNFF::~GFNFF`). The March 2026 "heap
  corruption" workaround leaked ~194 MB per structure (54 MB parameter set + ~140 MB GFNFF object,
  LD_PRELOAD tracker). Now `CMAKE_CUDA_FLAGS` carries the same pin, the wrapper frees both and keeps only
  bonds + bond_hb_data. Peak RSS 8/32/64 x polymer: 1.95/1.90/1.93 GB (was 2.7/7.3/13.3); 10/10 runs of 64
  structures clean under `GLIBC_TUNABLES=glibc.malloc.check=3`; energies identical.
  `CURCUMA_GFNFF_GPU_LEAK=1` restores the old leak. Whether F-Q9 (TECHNICAL_DEBT) had the same cause is
  plausible (its ASan build had vectorization off, i.e. no ABI split) but not proven: F-Q9 no longer
  reproduces with the Sep 27 binary (0/5).* [M] Original: the host GFNFF object and the full parameter set on destruction
  (`qm_methods/gfnff_gpu_method_impl.h:368-370, 500-505`); the comment says "~100 KB", but the set
  holds every pair list [C]. Size unmeasured.
- No GFN-FF batch timing exists; all pool numbers are GFN2 [D: MULTI_GPU.md:226-237].

**Determinism and visibility**

- **F-14** GPU results are not run-to-run reproducible (`atomicAdd` gradient, energy, CN, fragment
  sums, pair order) [C]; known, deferred by the operator [D: TODO.md:772-781, Sep 27].
- **F-19** *- FIXED Sep 28 (stage B): one ungated warning with the correct wording plus a counted fallback
  (only when the device solve was applicable; the first-step topology route is by design).* "All EEQ paths exhausted -> Phase-1 charges" warns only at verbosity >= 1
  (`qm_methods/gfnff_gpu_method_impl.h:1474-1479`), so MD children at verbosity 0 never show it [C].

**ROCm parity**

- **F-17** Every device EEQ variant is a stub on ROCm (GPU-Schur, WP7-A/C/E, batched;
  `ff_methods/rocm/eeq_solver_hip.hiph:406-445`); below 16 fragments a dense rocSOLVER solve with host
  Schur step, above it a host EEQ. Energy reductions hard-code wave32 (`ff_methods/rocm/gfnff_rocm.hip:86-128`),
  wrong on CDNA/MI (wave64). The September mirrors are uncompiled [C].

## 3. Cross-cutting

- **X-1 GPU/EEQ tuning parameters missing from the registry** *- FIXED Sep 28, 2026 (stage C): the parser
  tokenizes PARAM macros (strings, comments, adjacent literals) instead of matching a regex; 0 malformed
  warnings, 725 instead of 705 definitions. The 26 warnings were 20 dropped names plus fragments of PARAMs
  that were registered anyway. Old and new parser on the same inputs: the 705 common entries are
  identical (one, `eeq_refactor_force_every`, differs only by joined adjacent literals - the same runtime
  string). Two newly registered defaults were set to the values the code actually used, because the PARAM
  text had never been effective: `max_pcg_iterations` 200 -> 100 and `pcg_large_system_iterations`
  5000 -> 100 (CPU fallbacks in `eeq_solver.cpp`; the accuracy profiles in `accuracy_profile.cpp` that set
  other values are never called - dead code, noted). Routing: `-export_run` for 9 fixed commands differs
  only by the new `gfnff/solvent` + `gfnff/solvent_model` default entries and by the flat EEQ flags moving
  from the command module to `eeq_solver`; all energies identical. Behaviour change, intended: a flat
  `-solve_method pcg -max_pcg_iterations 300` was ignored by the CPU solver before (log: `solve_method=
  cholesky`) and is applied now (`solve_method=pcg`, 300 iterations); the GPU reads these keys from the
  `eeq_solver` scope first as well.* [M] Original finding: the build reports 26
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
- **X-3 Doc drift** *- FIXED Sep 28 (stage B), each place corrected with a dated note; the repulsion
  52 vs ~240 ms item turned out to be two different quantities (CPU rebuild vs GPU call saving), left
  as an open measurement.* [C]: MD_LARGE_SYSTEMS.md:117-131 ("GPU buys nothing") predates the Sep 27 changes;
  MULTI_GPU.md:28 (CPU/GPU 5.7e-7 Eh) is superseded by CLAUDE.md #33; TODO.md:1011-1018 says ConfSearch
  switches the GPU off at `-threads > 1` - true only while the device pool is inactive (one
  visible GPU); with an active pool the GPU stays on (`src/capabilities/confsearch.cpp:161-170`) [R]; GPU_TUNING.md:36 ("off by default") contradicts :40 and
  the code; GPU_TUNING.md:161 contradicts F-4; GFNFF_PAIR_LIST_REFRESH.md gives the repulsion rebuild
  as 52 ms (:116) and ~240 ms (:159); comments in `qm_methods/gfnff_gpu_method_impl.h:34-39, 959-965`
  describe the old CPU-EEQ pipeline.

## Measurements Sep 28, 2026 (stage D, after the fixes of stages A-C)

4x RTX A4500 (20 GB, no NVLink), i9-10980XE (18 cores); an unrelated ollama process held 7-9 GB on
every card during the runs (0 % load). Lab journal: vault `Labor/curcuma MD-Stabilität großer
Systeme.md`, entry Sep 28. GFN-FF MD: `-dt 1 -thermostat csvr -T 300 -seed 7 -threads 16`.

| what | result | note |
|---|---|---|
| GFN-FF MD step, polymer_2x (7320) | **GPU 0.38-0.39 s, CPU 16 threads 1.33 s** | difference of 5/35 and 5/95-step runs; Sep 21: 1.475 vs 1.25 |
| of that: energy call (GPU) | 171 ms = EEQ 77 + Coulomb/finish 94 | verbosity-2 report, late step |
| of that: SimpleMD `step_total` | 213 ms median (dump steps) | **~0.18 s/step not attributed**; no perf/nsys here |
| EEQ, one fragment, polymer (1410) | GPU dense Cholesky 8.6 ms; CPU CN+EEQ 36 ms (projected PCG) / 52 ms (exact) | F-4 pays only for large single-fragment systems; none available here |
| GFN-FF `-sp` batch, 8x polymer | **1 GPU 3.2 s, 4 GPUs 1.3 s, CPU 16 threads 0.98 s** | GFN-FF single points belong on the CPU |
| GPU batch peak RSS, 8/16/32x polymer | **2.7 / 4.2 / 7.3 GB** (CPU 1.3 / 2.3 / 2.3) | F-16 leak confirmed: ~190 MB per structure, linear |
| GFN2 polymer, split thresholds lowered to 1000 | 1/2/4 GPUs 8.95 / 11.7 / 12.8 s; 1 GPU: setup 1.94 s, SCF 4.0 s, post-SCF 2.34 s of which D4 1.86 s | proxy only; the split is a loss at this size (gate 4000 is right) |
| GFN2 polymer_2x on 1/2/4 GPUs (D1, D2 = G2-9) | **not measured** | needs 13.5-17.9 GB on device 0, ~11.5 GB free |

Measurement artefact found: `-md_diagnostics_timing` reports host prep times on the GPU path (`eeq_solve`
769 ms, same as the CPU run), and with `-dump 1` the diagnostics writing itself slows the step to ~1.5 s.

## Profiling Sep 29, 2026 (nsys + perf, agents; all timed runs serialised on one lock)

**GFN-FF MD step, polymer_2x, 1 A4500, 16 host threads** (before the F-8 fix: 0.38 s/step; after: 0.20 s):

| part | ms/step | note |
|---|---:|---|
| `k_coulomb_implicit` | 85 | all pairs, each twice |
| EEQ projected PCG (`symv`) | 53 | ~69 iterations x 0.76 ms, bandwidth-bound dense N x N |
| `k_eeq_build_matrix` | 14 | dense N x N every step |
| dispersion + dC6/dCN | 14 | |
| CN pair-list regeneration | 6 | fires on 16/25 steps |
| GPU idle (syncs, 208 `cudaStreamSynchronize`/step) | 5-6 | GPU busy 97 % of a normal step |
| forced HB/XB refresh incl. full topology | 160 | F-8, fixed |
| repulsion rebuild incl. bond list | 34 | bond list part fixed with F-8 |
| outside the energy call | <= 2 | |

**GFN2 polymer_2x SP (nao 15444), 1/2/4 A4500:** wall 287 / 227 / 143 s, SCF iterations 15 / 14 / 12
(part of the "speedup" is FP32 rounding luck). Per FP32 iteration only 1.78x on 4 GPUs. Not scaling at 4
GPUs (~55 of 140 s): setup 20.8 s (`k_overlap_h0` 11.6 s and Cholesky of S 3.5 s on device 0), post-SCF
5 s (D4 ATM 3.8 s on device 0), the FP64 back-transform `trsm` (10.6 / 20.8 / 10.9 s - does not scale;
**done Sep 29, 2026**: column split with a full L per device, 10.6 / 5.6 / 2.7 s),
~1.5 s/iteration device-0 share (potential 0.33, SCC energy 0.17, FP32 copy/back-transform 0.94 s).
FP32 `syevd` inside cuSOLVERMp is the largest bucket (4.8 s/it at 4 GPUs, NCCL 22-30 % of busy time).
G2-9 (buffer re-allocation between opt steps) measured <= 0.2 s per step - not worth it. Device-0 peak
15.3 / 18.5 / 15.8 GB for 1 / 2 / 4 GPUs (2 GPUs needs more than 1). `-opt` on 4 GPUs ran 5 SCFs for 3
printed steps (two extra at the start geometry, ~216 s) - under investigation.

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
