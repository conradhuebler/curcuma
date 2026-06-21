# GPU backend roadmap — remaining work packages (ROCm & Vulkan)

> Status: 🤖 AI-authored plan. These are **proposed** work packages, not implemented
> features. CUDA is the complete reference; ROCm and Vulkan are staged ports of it. Every
> AP below ports proven CUDA math through the existing `GpuScfBackend` seam — no new
> algorithms, so each step stays bit-checkable against the CPU/CUDA path.

## Where each backend stands

The native xTB SCF offloads through the abstract `GpuScfBackend` seam (`xtb_native.h`); any
hook returning `false` falls back to the CPU, so each backend is brought up stage-by-stage
and is correct at every step. Stage numbers below match the CUDA stages in
[SQM_GPU.md](SQM_GPU.md).

| Stage | CUDA | ROCm | Vulkan |
|-------|:----:|:----:|:------:|
| 0 build + dispatch + device handshake | ✅ | ✅ | ✅ |
| 1 GPU eigensolver (per-iteration `F C = S C ε`) | ✅ | ✅ (rocSOLVER `dsygvd`) | ✅ (FP64 Jacobi) |
| 2a GFN1 device-resident SCF | ✅ | ✅ | ✅ |
| 3 integrals on device — isotropic (CN/S/H0/L/γ) | ✅ | ✅ | ✅ |
| 4 nuclear gradient on device — GFN1 | ✅ | ✅ | ✅ (V-AP1, gather kernels) |
| 3m integrals on device — GFN2 multipole (dp/qp) | ✅ | ✅ (R-AP1) | ✅ (V-AP2) |
| 2b GFN2 device-resident multipole SCF | ✅ | ✅ (R-AP2) | ✅ (V-AP3) |
| 4m nuclear gradient on device — GFN2 multipole | ✅ | ❌ **R-AP3** | ❌ **V-AP4** |
| 5/6 device GFN2 potential + fully resident loop | ✅ | ❌ **R-AP4** (opt) | ❌ **V-AP5** (opt) |
| GFN-FF on the GPU | ✅ (full) | ✅ energy+grad (June 2026); EEQ = rocSOLVER `dpotrf`/`dgetrf` + host CPU-Schur (device GPU-Schur/PCG not ported) | ❌ **V-AP6** |
| in-SCF solvation (ALPB/GBSA) on device | ✅ | ❌ **X-AP1** | ❌ **X-AP1** |

**Summary:** Vulkan now matches ROCm/CUDA on the full GFN1 stack — integrals, resident
SCF **and** the on-device nuclear gradient (V-AP1 done), so GFN1 `-opt`/`-md` are fully
device-resident. Both backends now also build the **GFN2 multipole integrals
(dp_int/qp_int) on the device** (R-AP1 / V-AP2 done) and download them so the host GFN2
SCF skips its O(nao²) integral loop. They still trail CUDA on the GFN2 **resident**
multipole SCF + multipole gradient, GFN-FF, and device solvation. The s/p-only limit
(H/C/N/O…, no d shells) of the CPU native path applies to all GPU ports.

> **Vulkan note (no FP64 atomics):** Vulkan/GLSL has no `atomicAdd(double)` (unlike
> CUDA/ROCm), so the four scatter kernels are restructured as **per-atom GATHER** (one
> thread owns atom A's gradient slot and sums over all partners). The pair forces are
> antisymmetric (`f_ij = −f_ji`) and every H0-Pulay scalar term is pair-symmetric, so the
> gather reproduces the scatter exactly; the overlap derivative is evaluated with A's
> shell as the first argument (`dS/dR_A`) and `dxij = x_A − x_foreign`, carrying the sign.
> Cost is ~2× the overlap-gradient evaluations (each unordered AO pair seen from both
> atoms) — acceptable on an iGPU where this is a residency/correctness milestone.

## Port mechanics (shared)

- **Vulkan**: each kernel is a hand-written FP64 GLSL compute shader in
  `qm_methods/vulkan/shaders/*.comp`, compiled to SPIR-V (`compile_shaders.sh`, dev-time
  only) and embedded as `*.spv.inc` via `shaders/spirv_kernels.h`. No BLAS/LAPACK on
  Vulkan — GEMM/TRSM/eigensolve are also hand-written shaders (already present).
- **ROCm**: each kernel is a HIP `__global__`/`__device__` function in
  `qm_methods/rocm/xtb_hip_integrals.hiph` + `xtb_hip_context.hip`, compiled by `hipcc -c`
  into a plain relocatable object linked with g++ as an `EXTERNAL_OBJECT` (no
  `enable_language(HIP)`, no `--offload-arch`/`--hip-link` on the GNU link). rocBLAS/
  rocSOLVER available.
- **Port source (CUDA)**: device math lives in `cuda/xtb_gpu_integrals_device.cuh`
  (`cgto_multipole`, `d_cgto_multipole_grad_transformed`, the Obara-Saika overlap helpers)
  and the kernels/host glue in `cuda/xtb_gpu_context.cu` (`k_multipole_ints`,
  `solveMultipole`, `multipoleMoments`). ROCm already ports the *isotropic* subset
  verbatim into `xtb_hip_integrals.hiph`; the multipole helpers extend the same file.
- **Validation pattern (mandatory, per [CLAUDE.md](../CLAUDE.md) accuracy rules)**: prove
  each kernel elementwise vs the CPU build (component test, target ~1e-13…1e-15), then the
  end-to-end check — `-sp` energy bit-identical to CPU at 8 dp over the 12-molecule
  `test_cases/sqm_reference` set (incl. 231-atom `complex`), and `-opt` matching the CPU
  trajectory **step-by-step in energy AND gradient norm** (a wrong gradient diverges
  immediately, so this is the decisive gradient check). Add a ctest label mirroring CUDA's
  `gpu_integrals` / `gpu_gradient`.

---

## Vulkan work packages

> **Tag note:** `V-AP*` here is the **staging sequence** (mirrors the CUDA/ROCm stages). The
> Vulkan *performance* work — the Householder eigensolver and the workgroup-per-atom gradient —
> is a separate track tagged **`V-PERF-*`** (`V-PERF-1`/`V-PERF-2`), even though their commits
> `be159a1`/`61994ce` were authored as `V-AP5`/`V-AP6` before the retag. The `V-AP5`/`V-AP6`
> below are the planned device-potential-build / GFN-FF stages, unrelated to those commits. See
> `docs/SQM_VULKAN_EIGENSOLVER_WP.md`.

### V-AP1 — On-device nuclear gradient, GFN1 (Stage 4) — ✅ DONE (2026-06)
- **Goal**: bring Vulkan to ROCm/CUDA parity for GFN1; make GFN1 `-opt`/`-md`
  fully device-resident (only xyz up, gradient+energy down per step).
- **Done**: three SPIR-V gather kernels — `grad_rep` (section 1), `grad_coulomb`
  (section 3), `grad_pulay` (sections 2a on-site CN + 2b H0/Pulay off-site, with the
  inline Obara-Saika `cgto_overlap_grad` ported into the shader). The energy-weighted
  density `W = C·diag(2ε)·Cᵀ` is built on device via the existing `scale_cols`+`gemm`
  from the resident `rC` (Jacobi order; `rd[ridx[k]] = 2ε[k]`). `XtbVulkanContext::gradient`
  + `VulkanScfBackend::{supportsGradient→true, gradient(...)}` over the resident
  density/MO coefficients (`pc_resident`, no P/C upload). Dispersion gradient + CN
  chain-rule stay on the host (as CUDA/ROCm). No FP64 atomics (see Vulkan note above).
- **Validated** (AMD 890M / RADV): `-sp` GFN1 energy + gradient norm bit-identical to CPU
  over the full 12-molecule `sqm_reference` set incl. 231-atom `complex`; `-opt` caffeine
  35 steps step-by-step identical in energy AND gradient norm. ctest
  `cli_gpu_gradient_01_vulkan_gfn1_gradient` (labels `gpu_gradient;gpu;vulkan;cli`).
- **Unlocks**: closed the single biggest Vulkan gap; reuses the resident density from
  Stage 2a. GFN2 still falls back to the CPU gradient (host SCF — `use_gpu_resident` is
  false for GFN2 on Vulkan, so the device gradient is only ever taken for GFN1).

### V-AP2 — GFN2 multipole integrals on device (Stage 3m) — ✅ DONE (2026-06)
- **Goal**: build `dp_int`/`qp_int` (AO dipole/quadrupole) on the GPU so the host GFN2 SCF
  skips its O(nao²) `setupMultipole` integral loop.
- **Done**: SPIR-V `multipole_ints` shader (one thread per AO pair: global-origin
  `cgto_multipole` then the per-column origin shift + traceless transform using the
  resident overlap `rS`). `XtbVulkanContext::{beginMultipoleComputed, downloadMultipole}`
  + `gDpInt`(3·nao²)/`gQpInt`(6·nao²) buffers; no new basis fields (Stage 3/4 already
  upload ao2sh/ao2at/iao_sh/ang_sh/prim_*). Host: `GpuScfBackend::downloadMultipoleInts`
  hook + `setupMultipole(bool integrals_on_device)` skips steps 1-2; wired in
  `xtb_native.cpp` for `!supportsMultipole()` backends (CUDA's resident loop excluded).
  Since Vulkan GFN2 runs the host SCF, the device integrals feed the host Fock — bit-
  identical GFN2 energy proves them.
- **Validated** (AMD 890M / RADV): gfn2 `-sp` energy bit-identical to CPU (8 dp) over the
  full 12-molecule `sqm_reference` set incl. 231-atom `complex`; gfn2 `-opt` (NH3) identical
  trajectory (integrals rebuilt each geometry). ctest
  `cli_gpu_multipole_01_vulkan_gfn2_multipole`.
- **Depends on**: nothing (independent of V-AP1). Next: V-AP3 (resident multipole SCF).

### V-AP3 — GFN2 device-resident multipole SCF (Stage 2b) — ✅ DONE (2026-06)
- **Goal**: flip `supportsMultipole()→true` so GFN2 enters the resident loop instead of
  falling back to the Stage-1 eigensolver; keep the anisotropic channel on the device.
- **Done**: SPIR-V `fock_multipole` (adds `F −= ½·Σ_k(dp_int[k]·v_dp + qp_int[k]·v_qp)` to
  the resident isotropic Fock, one thread per AO pair) + `multipole_moments` (atomic
  dp_at/qp_at — per-atom GATHER, no FP64 atomics). `XtbVulkanContext::{solveMultipole,
  multipoleMoments, beginMultipoleComputed}` reuse the GFN1 resident Löwdin+Jacobi solve
  with the extra Fock term; v_dp/v_qp/dp_at/qp_at buffers + sets in begBasis. The SCF loop
  in `xtb_native.cpp` is method-agnostic (already drove the CUDA resident multipole path),
  so no loop change. Only `v_dp`(3·nat)/`v_qp`(6·nat) up, `dp_at`/`qp_at` + eps down per
  iteration (no nao² eigensolver transfer).
- **Gradient gotcha (fixed)**: making GFN2 resident sets `use_gpu_resident=true`, which
  routed GFN2 through `calculateGradientGpu` → the GFN1-only device `gradient()` (missing
  the multipole-integral Pulay). Fix: the device `gradient()` returns false for GFN2
  (non-empty v_dp/v_qp) → falls back to the full **host** `calculateGradient` (V-AP4 pending).
- **Validated** (AMD 890M / RADV): gfn2 `-sp` energy bit-identical (8 dp) over the
  12-molecule set incl. 231-atom `complex`; gfn2 `-sp` gradient norm matches CPU; gfn2
  `-opt` (NH3) step-by-step identical. **Honest: NOT a speed-up on this iGPU** — the FP64
  Jacobi eigensolve dominates (GFN1, already resident, is also ~1.3× the CPU on `complex`);
  the removed per-iteration transfers are cheap shared-memory copies on an iGPU, not PCIe.
  V-AP3 is the correct architecture + the prerequisite for V-AP4 and the discrete-GPU win;
  the iGPU speed lever is an **FP32 eigensolve** (see "FP32 mixed precision" below).
- **Depends on**: V-AP2 (resident `dp_int`/`qp_int`). Next: V-AP4 (GFN2 device gradient).

### V-AP4 — GFN2 nuclear gradient incl. multipole (Stage 4m)
- **Goal**: device-resident GFN2 `-opt`/`-md`.
- **Port from**: `d_cgto_multipole_grad_transformed` (the multipole-integral Pulay) + the
  SD/DD/SQ direct-interaction gradient (`get_multipole_gradient_0d` port already in the CPU
  path / CUDA Stage 4).
- **Depends on**: V-AP1 (the GFN1 gradient skeleton) + V-AP3 (resident moments).

### V-AP5 — (optional) device GFN2 potential build + fully resident loop (Stage 5/6)
- **Goal**: move the in-SCF D4/EEQ + isotropic+multipole potential and the
  occupation/populations/SCC-energy/Broyden loop body onto the GPU.
- **Honest note**: per the CUDA Stage-6 measurement this is a **residency/correctness
  milestone, not a measured `-sp` speed-up** (the eigensolve + O(1) host syncs dominate, not
  the removed transfers). Low priority on an iGPU where FP64 is already slow. Defer until a
  discrete GPU makes it worthwhile.

### V-AP6 — GFN-FF on Vulkan
- **Goal**: EEQ solve + GFN-FF energy/gradient terms on the GPU (CN, bonds, angles,
  torsions, dispersion, repulsion, Coulomb).
- **Note**: separate track from the GFN1/GFN2 tight-binding path; sizeable. Scope after the
  GFN2 stack lands. CUDA's GFN-FF stack is only partial, so this is partly net-new, not a
  pure port.

---

## ROCm work packages

ROCm already has GFN1 fully device-resident (Stage 4). The remaining gap is the **GFN2
anisotropic stack** (then GFN-FF). The HIP ports are lower-risk than Vulkan because
rocBLAS/rocSOLVER cover the dense linear algebra and the isotropic `.hiph` already exists.

### R-AP1 — GFN2 multipole integrals on device (Stage 3m) — ✅ DONE (2026-06)
- **Goal**: build `dp_int`/`qp_int` on the GPU so the host GFN2 SCF skips its integral loop.
- **Done**: `d_moment1d`/`d_type_to_cart`/`d_primitive_multipole`/`d_cgto_multipole` ported
  verbatim from the CUDA `.cuh` into `xtb_hip_integrals.hiph` (the `__device__` math
  compiles unchanged under hipcc) + the `k_multipole_ints` kernel (origin shift + traceless
  transform with the resident `dS0`) in `xtb_hip_context.hip`. `dDpInt`/`dQpInt` buffers +
  `XtbHipContext::{beginMultipoleComputed, downloadMultipoleInts}`; no new `XtbHipBasisData`
  fields (Stage 4 already uploads ao2sh/ao2at). Shared host plumbing with V-AP2.
- **Validated** (Radeon 890M/gfx1150): gfn2 `-sp` energy bit-identical to CPU (8 dp) over
  the full 12-molecule `sqm_reference` set incl. 231-atom `complex`; gfn2 `-opt` (NH3)
  identical trajectory. (CLI-validated like the other ROCm stages; no ctest harness.)
- **Depends on**: nothing. Next: R-AP2 (resident multipole SCF).

### R-AP2 — GFN2 device-resident multipole SCF (Stage 2b) — ✅ DONE (2026-06)
- **Goal**: GFN2 enters a device-resident loop instead of the Stage-1 eigensolver.
- **Done**: `k_add_fock_multipole` + `k_multipole_moments` (verbatim CUDA ports — HIP has
  `atomicAdd(double)`, so the moment scatter is direct) added to `xtb_hip_context.hip`;
  `XtbHipContext::{solveMultipole, multipoleMoments, beginMultipoleComputed}` reuse the
  GFN1 resident `k_fock` + rocSOLVER `dsygvd` with the extra Fock term. Same shared host
  plumbing + gradient fix as V-AP3 (GFN2 gradient stayed on the host at R-AP2; R-AP3 below
  moves it onto the device). Validated on gfx1150: gfn2 `-sp` energy bit-identical + gradient
  norm matches CPU; gfn2 `-opt` (NH3) step-by-step identical. Same iGPU honest note as V-AP3
  (eigensolve-bound, not a speed-up; FP32 is the lever).
- **Depends on**: R-AP1. Next: R-AP3 (GFN2 device gradient).

### R-AP3 — GFN2 nuclear gradient incl. multipole-integral Pulay (Stage 4m) — ✅ DONE (2026-06-16)
- **Goal**: move the GFN2 multipole-integral Pulay term onto the device so the GFN2 nuclear
  gradient is device-resident (it joins the GFN1 repulsion/Pulay/Coulomb kernels).
- **Done**: `d_cgto_multipole_grad_transformed` + helpers (`d_mpg_dGdA`,
  `d_primitive_multipole_grad`) ported verbatim from the CUDA `.cuh` into
  `xtb_hip_integrals.hiph`; `k_grad_h0_pulay` extended with the dp/qp-integral derivative term
  contracted against the converged `v_dp`/`v_qp` (uploaded once per gradient to `dVdp`/`dVqp`);
  `XtbHipContext::gradient` + `HipScfBackend` thread `v_dp`/`v_qp` through, and the GFN2 dispatch
  in `xtb_hip_method.cpp` now runs the device gradient (was: return false → host fallback). The
  multipole SD/DD/SQ **interaction** gradient (§5) + dispersion + CN chain-rule stay on the host.
- **Validated** (Radeon 890M/gfx1150): with `-scf_mixed_precision false` gfn2 `-sp` gradient norm
  bit-identical to CPU over `sqm_reference` incl. 231-atom `complex` (reldiff = 0); the
  ROCm-routed FD gradient check (`test_xtb_gradient`) is byte-identical to the CPU path per
  molecule (proving device analytic == FD-validated CPU AP5 gradient); gfn2 `-opt` (triose)
  tracks CPU then diverges in the LBFGS tail exactly like the known-good GFN1 control
  (eigensolver-backend noise, not R-AP3). `ctest -L native_xtb` 10/10. At the default FP32-mixed
  path the gradient differs ~1e-7 (the X-AP3 lever).
- **Depends on**: R-AP1/R-AP2 (the resident `dp_int`/`qp_int` + `v_dp`/`v_qp`).

### R-AP4 — (optional) device GFN2 potential + fully resident loop (Stage 5/6)
- Same correctness-not-speed caveat as V-AP5. Defer.

### R-AP5 — GFN-FF on ROCm
- Same scope as V-AP6; the HIP route can reuse rocBLAS for the dense pieces.

---

## Cross-cutting

### X-AP1 — In-SCF solvation (ALPB/GBSA) on the device
- **Goal**: fold the Born reaction field `v_at += B·q_solute` into the device potential so
  `-xtb.solvent_model alpb|gbsa` works on `-gpu rocm`/`vulkan` (CUDA already does — the
  reaction field is folded into the uploaded `v_ao`; the fully-resident loop is
  solvent-aware).
- **Note**: for the host-driven SCF (current ROCm/Vulkan GFN2) this is the cheap interim
  path — add `B·q` to the host `v_ao` before upload. Only the fully-resident loop (AP5)
  needs the Born matrix on the device.

### X-AP2 — Validation harness parity
- **Goal**: add `ctest` labels for ROCm/Vulkan mirroring CUDA's `gpu_integrals` /
  `gpu_gradient` / `gpu_gfn{1,2}_validation`, plus per-kernel component tests
  (device-vs-CPU elementwise). Run `compute-sanitizer`-equivalent (`rocgdb`/RADV
  validation layers) where available.

### X-AP3 — FP32 mixed-precision eigensolve (the iGPU speed lever) — ✅ DONE (2026-06)
- **Why**: residency (V-AP3/R-AP2) is correct but does **not** beat the CPU on the AMD 890M
  because the **FP64** eigensolve dominates. Consumer/iGPU FP64 runs at ~1/16 of FP32, so
  an FP32 eigensolve far from convergence (FP64 only for the final iterations) was the
  hypothesised speed-up.
- **Plumbing**: `-scf_mixed_precision` (+ `-scf_fp32_threshold`) sets `m_eig_fp32` per
  iteration; the SCF loop passes it into `solve()`/`solveMultipole()`; the wrappers forward
  it to the backend contexts. Convergence is never accepted on an FP32 step, so the final
  density is FP64 and the energy matches FP64 within `scf_threshold`.
- **ROCm — real win, default ON**: `rocsolver_ssygvd` (FP32) on FP32 copies of F/S, the
  eigenvectors/eigenvalues cast back to FP64 (`k_d2f`/`k_f2d`); `dS0` (resident overlap)
  stays intact. Measured on a Radeon 890M, `complex` (231 atoms): **GFN1 1844→1281 ms
  (1.44×), GFN2 1924→1510 ms (1.27×)** over the resident FP64 path (which is already ~6-8×
  the CPU because rocSOLVER is a mature kernel). Energies bit-identical to CPU. Mixed
  precision is **defaulted ON for `-gpu rocm`** (matches CUDA). Caveat: at the loose default
  `scf_threshold` (1e-5) the FP32 path reaches a slightly different fixed point, so the
  **gradient norm differs ~1e-7 from CPU** (immaterial to `-opt` convergence, but not
  bit-identical); `-scf_threshold 1e-8` restores bit-identical gradients. Opt-out:
  `-scf_mixed_precision false`.
- **Vulkan — no win, opt-in only**: FP32 variants of the two-sided Jacobi shaders
  (`angles_f32`/`col_f32`/`row_f32`/`vec_f32`) + FP32 work buffers (`rAtil32`/`rCtil32`/
  `rCs32`); Ã cast FP64→FP32 before the FP32 Jacobi, C̃ back to FP64. **Measured net-neutral
  to net-negative**: the hand-written cyclic Jacobi is dispatch/barrier/bandwidth-bound, not
  FP64-arithmetic-bound, so FP32 is ~equal per iteration (`complex` GFN2: FP32 828 vs FP64
  844 ms/iter) and the perturbed early iterations occasionally cost an extra SCF cycle
  (GFN1 14→15), making it slower overall. The real Vulkan lever is the **eigensolve
  algorithm** (the naive Jacobi), not its precision. So mixed precision is **NOT defaulted on
  for `-gpu vulkan`** — the FP32 path stays as a correct, documented opt-in
  (`-scf_mixed_precision true`). Energies still match CPU.
- **Lesson**: the iGPU bottleneck is backend-specific. ROCm's rocSOLVER benefits from lower
  precision; Vulkan's bottleneck is the algorithm, which precision cannot fix. A faster
  Vulkan eigensolve (blocked/one-sided Jacobi, or a divide-and-conquer port) is the open lever.
- **Depends on**: nothing (orthogonal to the GFN2 stack). Independent per backend.

### X-AP4 — In-SCF GFN2 D4 dE_D4/dq on device — ✅ DONE (2026-06, both backends)
- **Why**: with `supportsDevicePotential()` still false on ROCm/Vulkan, the GFN2 SCF runs the
  host potential build every iteration, dominated by the O(N²) D4 dE_D4/dq
  (`XTB::addDispersionPotential` → `computeD4PotentialDedq`). CUDA alone routed this to the
  device (`supportsDeviceDispersion()`); this WP ports that one hook to the other two backends.
- **Done**: `supportsDeviceDispersion()`→true + `beginDispersion`/`dispersionDedq` on both
  `HipScfBackend`/`XtbHipContext` (HIP `k_d4_dedq`, verbatim CUDA port) and
  `VulkanScfBackend`/`XtbVulkanContext` (FP64 SPIR-V `d4_dedq.comp`, 9 SSBO + `uint nat` push,
  per-atom GATHER, no atomics). The host plumbing (`xtb_native.cpp` `use_device_disp`:
  geometry-fixed reference upload once + per-iteration `buildRefWFlat` → `dispersionDedq`) was
  already present, so only the three overrides + the kernel were added. Element C6 block
  uploaded once/process; Z/nref/√r4r2/xyz per geometry; W/dWq per iteration.
- **Validated** (Radeon 890M / RADV + gfx1150): GFN2 `-sp` energy bit-identical to CPU over
  `sqm_reference` incl. 231-atom `complex`; GFN2 `-opt` (NH3) same minimum; device path
  confirmed active by per-iteration invocation trace (not a silent host fallback); GFN1 + CPU
  byte-unchanged; Vulkan `ctest -L gpu` 22/22.
- **Honest**: **not a per-iteration speed-up at tested sizes** — the threaded host D4 evaluator
  is ~3.2 ms/it on `complex`/231 and the device round-trip (W/dWq up + kernel + dEdq down,
  sync/submit-bound on the iGPU) is ~4.0 ms/it on both backends. The benefit scales with N
  (O(N²) host cost) and only lands fully when fused into a device-resident potential/loop
  (V-AP5/R-AP4, Stage 5/6) so dEdq never round-trips.
- **2-body D4 nuclear gradient — DONE (2026-06, ROCm + Vulkan)**: the post-SCF 2-body D4
  (energy + grad + dE/dCN + dE/dq) runs on the device in one gather (`k_d4_grad` HIP /
  `d4_grad.comp` SPIR-V, a superset of the dEdq kernel reusing beginDispersion's resident
  reference data). Seam `GpuScfBackend::dispersionGradient`; `calcDispersionEnergy` calls it
  behind `supportsDeviceDispersion()`, then **ATM + the CN-distribution + the q-response stay
  host** and sum on top. Validated: GFN2 `-opt` (triose) Vulkan-FP64 trajectory bit-identical
  to CPU (every column incl. Energy-Change + gradient norm); ROCm 6-dp (FP32-SCF delta only);
  device path confirmed active (trace). **Still not a speed-up** (once/geometry) — pure
  architectural residency, the prerequisite for a fully device-resident GFN2 `-opt`/`-md`.
- **ATM 3-body D4 gradient — DONE (2026-06, ROCm + Vulkan)**: the post-SCF ATM (energy + grad +
  dE/dCN) runs on the device in a per-atom gather over pairs (`k_d4_atm` HIP / `d4_atm.comp`
  SPIR-V). The triple-edge force is antisymmetric (dG_ba=-dG_ab, dang symmetric in the edge
  endpoints) → no atomics, each triple's energy counted once per member (E=⅓·Σ). The host does the
  cheap O(N²) q=0 reference C6/∂C6∂CN build (`buildAtmC6Flat`) + uploads the nat² matrices; the
  device does the O(N³) triple loop reusing the resident geometry+√r4r2. Vulkan computes the
  non-integer damping power `(r0/r1)^(alp/3)` as `cbrt(b^alp)` via Newton (no `dlog`). Seam
  `GpuScfBackend::dispersionATM`. Validated: GFN2 `-opt` (triose) Vulkan-FP64 bit-identical to CPU
  (every column); ROCm 6-dp; energy 8-dp; device path confirmed (trace). Still no speed-up
  (once/geometry).
- **EEQ q-response ∂q/∂x — DONE on ROCm (2026-06) → fully device-resident ROCm GFN2 D4 gradient**:
  the 6 `k_d4eeq_*` kernels (CN / augmented matrix / RHS / adjoint / response gather) ported from
  CUDA + the (N+1) LU via rocSOLVER `dgetrf`/`dgetrs` (factor reused for the adjoint). The
  `GpuScfBackend::eeqCharges/eeqChargeResponse` seam + host integration already existed; only the
  ROCm backend was added. `supportsDeviceEeq()`→true also routes the scf_guess=eeq through the
  device. Validated: GFN2 `-opt` (triose) FP64 bit-identical to CPU (every column), `-sp` energy
  8-dp + same SCF iters, trace-confirmed active, `ctest` 34/34. **ROCm now does the WHOLE GFN2 D4
  gradient on device** (in-SCF dEdq + 2-body + ATM + q-response).
- **EEQ q-response ∂q/∂x — DONE on Vulkan (2026-06-18) → fully device-resident Vulkan GFN2 D4
  gradient**: the missing device dense (N+1) solve is a hand-written single-workgroup FP64
  **Gaussian elimination without pivoting** (`d4eeq_solve.comp`) — the screened-Coulomb hardness
  block + the charge-constraint border is non-singular without row swaps, so the no-pivot factor
  agrees with the host/ROCm partial-pivoted LU to rounding. The other three kernels (`d4eeq_cn`
  raw CN, `d4eeq_build` augmented matrix, `d4eeq_resp` ∂q/∂x gather) port the CUDA/HIP math
  verbatim; GLSL-fp64 has no `erf`, so a Numerical-Recipes `derf` (erfcc, ~1.2e-7) was added to
  `vk_integrals.glsl` (ample — q-response is a small gradient term). The O(N) log-compression /
  RHS / u-weight stay on the host (mapped buffers), as ROCm does on the device. `eeqCharges`/
  `eeqChargeResponseGradient` on the context + the `supportsDeviceEeq()→true` wrapper override
  reuse the existing CUDA-shared `GpuScfBackend` seam + host integration; `scf_guess=eeq` also
  routes through the device. Validated (AMD 890M/RADV): GFN2 `-sp` E+gradnorm match CPU (H2O/HCN/
  NH3, 7 dp); device EEQ charges = host model exactly; per-iter SCF trace bit-identical from iter 1,
  final E identical, iter-0 ±1e-10 (the `derf` signature confirms the device path runs, not a
  silent host fallback); GFN2 `-opt` (HCN) trajectory bit-identical to CPU; `ctest -L vulkan`
  22/22. **Vulkan now does the WHOLE GFN2 D4 gradient on device, on par with ROCm.**
- **CUDA — TODO, now the laggard for the D4 nuclear gradient (inversion)**: the 2-body + ATM
  D4-gradient kernels (`k_d4_grad`/`d4_grad.comp`, `k_d4_atm`/`d4_atm.comp`) + the
  `dispersionGradient`/`dispersionATM` seam are **net-new for ROCm/Vulkan — they do NOT exist in
  the CUDA `.cu`** (CUDA only has the in-SCF `dispersionDedq` + the EEQ `eeqChargeResponse`). So
  CUDA still computes the **2-body + ATM D4 gradient on the host** (`disp_grad_device=false` →
  host `computeEnergyAndGradient`+`computeATM`). Result correct, but CUDA's "fully device-resident
  `-opt`/`-md`" story now has a D4-gradient host hold-out that **ROCm does not**. To catch up, port
  `k_d4_grad`/`k_d4_atm` into `cuda/xtb_gpu_context.cu` + override `dispersionGradient`/
  `dispersionATM` in `xtb_gpu_method.cpp` (straight port of the HIP kernels). **D4-gradient
  residency by backend:** ROCm = full (dedq+2body+ATM+q-response); Vulkan = full
  (dedq+2body+ATM+q-response, since 2026-06-18); CUDA = only dedq + q-response (2-body + ATM still host).
- **Depends on**: nothing (orthogonal). Independent per backend.

## Suggested order

1. ~~**V-AP1** (Vulkan GFN1 gradient) — closes the only Vulkan↔ROCm asymmetry.~~ ✅ DONE (2026-06).
2. ~~**R-AP1 + V-AP2** (GFN2 multipole integrals, both backends) — shared port.~~ ✅ DONE (2026-06).
3. ~~**R-AP2 + V-AP3** (resident multipole SCF) — correct but iGPU-eigensolve-bound.~~ ✅ DONE (2026-06).
4. ~~**X-AP3** (FP32 mixed-precision eigensolve) — the actual iGPU speed lever.~~ ✅ DONE (2026-06):
   real 1.27-1.44× on ROCm (rocsolver_ssygvd, default ON); net-neutral on Vulkan (Jacobi is
   dispatch-bound, not FP64-bound — opt-in only). Vulkan eigensolve algorithm is the open lever.
5. **R-AP3 + V-AP4** (GFN2 device gradient) → device-resident GFN2 `-opt`/`-md`.
6. **X-AP1** (device solvation), **GFN-FF** (R-AP5 / V-AP6), **Stage 5/6** last.
