# Native GFN1/GFN2 on the GPU (CUDA)

> Status: 🤖 AI-generated, ⚙️ machine-tested. **Not** human production tested.
> Validation is vs the committed **tblite** references (1e-8 Eh gate) and vs the
> CPU native path; bit-identical GPU↔CPU is impossible (GEMM reduction order ≠
> MKL), the converged fixed point agrees to FP64 rounding (≪ 1e-8).

## What it is

`-method gfn1|gfn2 -gpu cuda` runs the native xTB tight-binding solver with its
heavy linear algebra on an NVIDIA GPU (cuSOLVER / cuBLAS), mirroring the existing
GFN-FF CUDA stack. `-gpu auto` uses the GPU when the build and a device are
available, else falls back to the CPU. The CPU path (`XTB` class, all `xtb_*.cpp`)
is unchanged and `#ifdef`-free — the GPU variant is a sibling, selected by the
factory, that **owns** a validated CPU `XTB` (host setup + properties) and a
device context (`XtbGpuContext`, cuSOLVER/cuBLAS handles + one stream).

This is a staged port. Each stage is independently buildable and validated at the
1e-8 Eh gate before the next begins.

| Stage | Scope | Status |
|------|-------|--------|
| 0 | Scaffolding: context, method wrapper, CMake sub-option, factory branch (pass-through to CPU) | ✅ done |
| 1 | Per-iteration generalized eigensolve on the GPU (cuSOLVER), rest on CPU | ✅ done |
| 2a | **GFN1 device-resident SCF**: H0/S/L stay on the device for the whole loop; per-iter only length-nao vectors cross | ✅ done |
| 2b | **GFN2 device-resident SCF**: + multipole integrals resident, anisotropic Fock + atomic moments on the device (isotropic potential incl. in-SCF D4 stays on host) | ✅ done |
| 3 | Integrals on the GPU (S, H0, CN, γ, multipole) | ⏳ |
| 4 | Gradients on the GPU → device-resident `-opt`/`-md` | ⏳ |
| 5 | Mixed precision (FP32 bulk / FP64 refinement) + tuning | ⏳ |

## Build

A dedicated build directory (mirrors `release_tblite/`; **not** the canonical
`release/`):

```bash
cmake -S . -B release_cuda -DCMAKE_BUILD_TYPE=Release -DC17=ON \
  -DUSE_MKL=ON -DUSE_BLAS=ON -DMKL_ROOT=/opt/intel/oneapi/mkl/latest \
  -DUSE_AVX2=ON -DUSE_AVX512=ON -DUSE_MARCH_NATIVE=ON -DUSE_PCH=ON \
  -DUSE_CUDA=ON -DUSE_CUDA_XTB=ON -DCMAKE_CUDA_ARCHITECTURES=120
cmake --build release_cuda -j8
```

- `USE_CUDA_XTB` (defaults to `USE_CUDA`) selects which `.cu` sources are
  compiled — it lets you build the GFN-FF GPU stack without the xTB kernels.
  The code gate is the single macro `USE_CUDA`; `USE_CUDA_XTB` is only a
  feature-availability flag for the factory.
- `-DCMAKE_CUDA_ARCHITECTURES=120` (Blackwell / RTX 5080) keeps compile fast.
- Reconfiguring an existing build dir offline: add
  `-DFETCHCONTENT_FULLY_DISCONNECTED=ON`. Do not run a bare `cmake .`.

## How Stage 2a works (GFN1 device-resident SCF)

Only GFN1 with the default Broyden charge mixing takes the resident path; GFN2
and the non-Broyden SCF modes keep the Stage-1 per-iteration eigensolver hook
(still a GPU eigensolve, CPU everything else). The design keeps every subtle piece
of physics and convergence on the **validated CPU path** and moves only the
O(nao²/nao³) kernels to the device:

- **Resident on the device for the whole loop** (uploaded once per geometry):
  `H0`, overlap `S`, lower Cholesky factor `L` (S = L·Lᵀ), and the per-iteration
  density `P` and MO coefficients `C`.
- **Per iteration, only length-nao vectors cross the bus**: the AO potential
  `v_ao` up; eigenvalues `eps`, AO populations `pop_ao` and the band energy
  scalar down. The nao×nao matrices never move during the loop (Stage 1 still
  re-uploaded the Fock every iteration).
- **On the device** (`XtbGpuContext::resident{Begin,Solve,Density,Finalize}`,
  `qm_methods/cuda/xtb_gpu_context.cu`):
  - Fock `F = H0 − ½·S·(v_ao⊕v_ao)` (one kernel, fused into the solve),
  - generalized eigensolve via the cached `L` (2× `cublasDtrsm` + `cusolverDnDsyevd`),
  - density `P = C·diag(occ)·Cᵀ` (column scale + `cublasDgemm`),
  - Mulliken AO populations `pop_ao(μ) = Σ_ν P_μν·S_μν` (reduction kernel),
  - band energy `Σ_μν P_μν·H0_μν` (`cublasDdot`).
- **On the host** (unchanged `XTB::Calculation`): potential build (Coulomb
  `γ·q_sh`, third-order), occupation / Fermi smearing, Mulliken shell/atom
  aggregation, the Coulomb + third-order energies, **Broyden mixing**, and
  convergence. After convergence `finalize()` downloads `P`/`C` once for the
  post-SCF energies and the (still CPU) gradient.

The seam is the abstract `GpuScfBackend` interface (`xtb_native.h`), which speaks
only project types (`Matrix`/`Vector`/`Eigen`) so the core `XTB` class stays free
of every CUDA header. Set to `nullptr` (default, and always in a non-CUDA build)
it is completely inert — the CPU SCF is byte-unchanged.

## How Stage 2b works (GFN2 device-resident SCF)

GFN2 adds an anisotropic (dipole/quadrupole) channel on top of the isotropic
GFN1 loop. Stage 2b keeps that channel's heavy parts on the device too, so the
whole GFN2 SCF (Broyden) is device-resident; only the cheap O(nat²) bits and the
mixing stay on the host.

- **Additionally resident** (uploaded once per geometry, `beginMultipole`): the
  3 dipole + 6 quadrupole AO integral matrices (`m_dp_int`, `m_qp_int`, each
  n×n) and the AO→atom map.
- **Device, per iteration:**
  - `solveMultipole`: the isotropic Fock kernel **plus** the multipole
    contribution `F −= ½·Σ_k(dp_int[k]·v_dp + qp_int[k]·v_qp)` (the GFN2 branch
    of `buildFock`), then the same cached-L eigensolve as GFN1.
  - `multipoleMoments`: the atomic dipole/quadrupole moments
    `dp_at(k,iat) = −Σ_{μ∈iat}Σ_ν P_νμ·dp_int[k]_νμ` (and quadrupole), one thread
    per AO with an `atomicAdd` scatter into the owning atom — the multipole block
    of `updatePopulations`.
- **Host, per iteration (unchanged CPU code):** the multipole *potential*
  `v_dp`/`v_qp`/`v_at`-shift (`addMultipolePotential`, O(nat²)), the shell
  third-order, the **in-SCF D4** charge coupling (`addDispersionPotential` folds
  `dE_D4/dq` into `v_at` from the host-resident charges — no device work, no
  special callback), the Fermi occupation, the Mulliken shell/atom aggregation,
  the Coulomb/third-order/multipole energies, the Broyden mix of the larger SCC
  vector (`q_sh` + `dp_at` + `qp_at`), and convergence.

So GFN2's only extra device traffic per iteration is `v_dp` (3·nat) and `v_qp`
(6·nat) up, and `dp_at`/`qp_at` (9·nat) down — all small. The n×n multipole
integrals never move during the loop.

## Validation

```bash
cd release_cuda
ctest -L "gpu_gfn1_validation"   # GFN1, -gpu cuda, 1e-8 vs tblite (resident path)
ctest -L "gpu_gfn2_validation"   # GFN2, -gpu cuda, 1e-8 vs tblite (Stage-1 path)
ctest -L "gfn1_validation" -L "gfn2_validation"   # CPU must stay green too
```

Measured (2026-06-03, RTX 5080, CUDA 13.1):

- `gpu_gfn1_validation` (resident GFN1, Stage 2a) **12/12** and
  `gpu_gfn2_validation` (resident GFN2, Stage 2b) **12/12** — vs tblite at 1e-8,
  the documented xfails (GFN1 He2 ~1.5e-8 floor + complex; GFN2 complex)
  reproduced exactly. CPU `gfn{1,2}_validation` **12/12**.
- GPU-vs-CPU total energy matches to all printed decimals on H2O, caffeine,
  triose, acetic_acid_dimer and **complex (231 atoms)** for both GFN1 and GFN2 —
  ≤ 1e-8 Eh.
- `compute-sanitizer --tool memcheck`: **0 errors** on the GFN1 and the GFN2
  (atomicAdd-scatter) resident paths.
- `release/` (canonical, no CUDA) rebuilds with the core edits and passes the CPU
  `native_xtb` + `gfn{1,2}_validation` suites — the seam is `#ifdef`-free.

## Stage 3 — integrals on the device (2026-06-03, DONE + validated)

The one-time-per-geometry integral build now runs on the GPU, so the resident SCF
**computes** S/H0/L instead of receiving them — no nao²-sized matrix is uploaded
per geometry. New device header `cuda/xtb_gpu_integrals_device.cuh` carries the
`__device__` ports; `XTB::exportGpuBasis` flattens the basis/H0 (post-ortho
primitives, GFN1 valence flags, per-shell Coulomb hardness, ao2at/ao2sh) into
CUDA-free bundles; `XtbGpuContext::beginBasis` uploads them once and
`computeIntegrals` runs `k_cn → k_self_energy → k_overlap_h0 → cusolverDnDpotrf(L)
→ k_gamma (+k_multipole_ints for GFN2)`. The seam prefers the device-computed path
and falls back to the Stage-2 upload. Each kernel matches the CPU **elementwise**:
CN/self-energy bit-identical, S/H0/L ~1e-15, γ ~5e-17, GFN2 dp_int/qp_int ~1e-13
(`ctest -L gpu_integrals`, 42 component tests).

## Stage 4 — nuclear gradient on the device (2026-06-03, DONE + validated)

The native gradient runs on the GPU when the SCF is device-resident, so
`-opt`/`-md` are **fully device-resident** (integrals + SCF + gradient; only xyz
up, gradient+energy down per step). GFN1: repulsion / on-site CN / H0-Pulay
(reusing the validated `d_cgto_overlap_grad`) / isotropic Coulomb on the GPU
(`k_grad_*`), energy-weighted density `W=2·Σ ε_i c_i c_iᵀ` built on the device.
GFN2 adds the multipole-integral Pulay (`d_cgto_multipole_grad_transformed`) on the
GPU; the direct SD/DD/SQ interaction + mrad/CN chain (section 5) and the dispersion
gradient + CN chain-rule stay on the host. Full gradient device-vs-CPU matches to
**~1e-15 Eh/Å** on H2O/CH4/NH3/C6H6/HCN/triose (both methods, `ctest -L
gpu_gradient`); GPU `-opt` converges to the same geometry/energy as CPU;
compute-sanitizer 0 errors. The `γ`/CN host-consumer choice was option (a) — γ is
computed on the device and downloaded once per geometry (nsh², cheap); moving the
Coulomb gemv on-device (option b) would *add* per-iteration transfers until the
whole isotropic potential moves to the device.

## Notes / limits

- GeForce FP64 is ~1/64 of FP32, so a pure-FP64 resident GFN1 SCF is **not**
  necessarily faster than the 8-core MKL CPU on small/medium systems; Stage 2a is
  a correctness + residency milestone. Device-resident `-opt`/`-md` (no host
  round-trip per step) landed with Stage 4.
- **Mixed-precision eigensolve (default ON for `-gpu`)**: far-from-convergence SCF
  iterations solve in FP32 (`cusolverDnSsyevd`), reverting to FP64 near convergence
  (max|dq| < `fp32_threshold`, default 1e-3) so the converged energy is FP64 and
  the 1e-8 gate holds. Convergence is never accepted on an FP32 step. This closed
  the gap to the threaded CPU: `complex`/231 gfn2 `-sp -threads 8` eigensolve
  519→278 ms, TOTAL ~515 ms — **~8% faster than CPU `-threads 8`** (~554 ms). gfn1
  is on par (both bounded by the host D3 finite-difference dispersion gradient).
  Use `-threads N` so the host setup/potential build parallelise alongside the GPU.
- Memory headroom (16 GB): `complex` (231 atoms) is ~tens of MB; the dense path's
  practical ceiling is ~`polymer` (1410 atoms, ~3–4 GB). Larger → large-system
  modes (CPU-orchestrated; the GPU dense kernel is the per-block solver).
- Large-system modes (`-large_system_mode fragments|dc|sparse`) are **not** wired
  to `-gpu` yet; with `-gpu` GFN1/Broyden they run the CPU fragment/DC driver.
