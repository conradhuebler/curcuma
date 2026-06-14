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
| 3m integrals on device — GFN2 multipole (dp/qp) | ✅ | ❌ **R-AP1** | ❌ **V-AP2** |
| 2b GFN2 device-resident multipole SCF | ✅ | ❌ **R-AP2** | ❌ **V-AP3** |
| 4m nuclear gradient on device — GFN2 multipole | ✅ | ❌ **R-AP3** | ❌ **V-AP4** |
| 5/6 device GFN2 potential + fully resident loop | ✅ | ❌ **R-AP4** (opt) | ❌ **V-AP5** (opt) |
| GFN-FF on the GPU | partial | ❌ **R-AP5 / V-AP6** | ❌ |
| in-SCF solvation (ALPB/GBSA) on device | ✅ | ❌ **X-AP1** | ❌ **X-AP1** |

**Summary:** Vulkan now matches ROCm/CUDA on the full GFN1 stack — integrals, resident
SCF **and** the on-device nuclear gradient (V-AP1 done), so GFN1 `-opt`/`-md` are fully
device-resident. Both backends still trail CUDA on the entire GFN2 anisotropic
(dipole/quadrupole) stack, GFN-FF, and device solvation. The s/p-only limit (H/C/N/O…,
no d shells) of the CPU native path applies to all GPU ports.

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

### V-AP2 — GFN2 multipole integrals on device (Stage 3m)
- **Goal**: build `dp_int`/`qp_int` (AO dipole/quadrupole) on the GPU so GFN2 stops
  uploading the 9 nao² matrices per geometry.
- **Port from**: `cgto_multipole` (`xtb_gpu_integrals_device.cuh`) + `k_multipole_ints`.
- **New/changed**: SPIR-V `multipole_ints` shader; `XtbVulkanContext::beginMultipoleComputed()`
  + a `downloadMultipole`/resident-buffer path; extend `XtbVulkanBasisData` with the
  multipole basis fields.
- **Depends on**: nothing (independent of V-AP1).

### V-AP3 — GFN2 device-resident multipole SCF (Stage 2b)
- **Goal**: flip `supportsMultipole()→true` so GFN2 enters the resident loop instead of
  falling back to the Stage-1 eigensolver; keep the anisotropic channel on the device.
- **Port from**: CUDA `solveMultipole` (isotropic Fock **+** `F −= ½·Σ_k(dp_int[k]·v_dp +
  qp_int[k]·v_qp)`) and `multipoleMoments` (atomic dp/qp moments via per-AO `atomicAdd`
  scatter). The multipole *potential* (`v_dp`/`v_qp` from the moments) stays on the host as
  in CUDA — only `v_dp`(3·nat)/`v_qp`(6·nat) up, `dp_at`/`qp_at`(9·nat) down per iteration.
- **New/changed**: SPIR-V `fock_multipole` + `multipole_moments` shaders;
  `beginMultiple`/`beginMultipoleComputed` resident hooks in `VulkanScfBackend`.
- **Depends on**: V-AP2 (needs the resident `dp_int`/`qp_int`).

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

### R-AP1 — GFN2 multipole integrals on device (Stage 3m)
- **Goal**: build `dp_int`/`qp_int` on the GPU (`xtb_hip_integrals.hiph` currently states
  "no multipole").
- **Port from**: `cgto_multipole` + `k_multipole_ints` (CUDA) → HIP `__device__`/`k_`
  equivalents in `xtb_hip_integrals.hiph`; extend `XtbHipBasisData` + `beginComputed`.
- **Depends on**: nothing.

### R-AP2 — GFN2 device-resident multipole SCF (Stage 2b)
- **Goal**: GFN2 enters a device-resident loop (today it uses device integrals + the
  Stage-1 rocSOLVER eigensolver per iteration with the host assembling the multipole Fock).
- **Port from**: CUDA `solveMultipole` + `multipoleMoments` (HIP kernels + rocBLAS).
- **Depends on**: R-AP1.

### R-AP3 — GFN2 nuclear gradient incl. multipole (Stage 4m)
- **Goal**: device-resident GFN2 `-opt`/`-md` (GFN2 gradient is on the host today).
- **Port from**: `d_cgto_multipole_grad_transformed` + SD/DD/SQ direct gradient → HIP, reusing
  the existing `k_grad_*` skeleton.
- **Depends on**: R-AP1 (the resident `dp_int`/`qp_int`).

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

## Suggested order

1. ~~**V-AP1** (Vulkan GFN1 gradient) — closes the only Vulkan↔ROCm asymmetry.~~ ✅ DONE
   (2026-06). **Next up:** step 2.
2. **R-AP1 + V-AP2** (GFN2 multipole integrals, both backends) — shared port, independent.
3. **R-AP2 + V-AP3** (resident multipole SCF) → **R-AP3 + V-AP4** (GFN2 gradient).
4. **X-AP1** (device solvation, interim host-`v_ao` form is cheap) alongside step 3.
5. **GFN-FF** (R-AP5 / V-AP6) and the **Stage 5/6** APs last (largest / lowest measured ROI).
