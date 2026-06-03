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
| 2b | GFN2 device-resident SCF (multipole Fock + multipole Mulliken + in-SCF D4 host callback) | ⏳ next |
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

## Validation

```bash
cd release_cuda
ctest -L "gpu_gfn1_validation"   # GFN1, -gpu cuda, 1e-8 vs tblite (resident path)
ctest -L "gpu_gfn2_validation"   # GFN2, -gpu cuda, 1e-8 vs tblite (Stage-1 path)
ctest -L "gfn1_validation" -L "gfn2_validation"   # CPU must stay green too
```

Measured (2026-06-03, RTX 5080, CUDA 13.1):

- `gpu_gfn1_validation` (resident GFN1) **12/12** — 10 molecules ≤ 1e-8 vs tblite,
  the two documented xfails (He2 ~1.5e-8 floor, complex) reproduced exactly.
- `gpu_gfn2_validation` **12/12** (unchanged), CPU `gfn{1,2}_validation` **12/12**.
- GFN1 GPU-vs-CPU total energy matches to all printed decimals on H2O, caffeine,
  triose and **complex (231 atoms)** — ≤ 1e-8 Eh.
- `compute-sanitizer --tool memcheck`: **0 errors** on the resident path.
- `release/` (canonical, no CUDA) rebuilds with the core edits and passes the CPU
  `native_xtb` + `gfn{1,2}_validation` suites — the seam is `#ifdef`-free.

## Notes / limits

- GeForce FP64 is ~1/64 of FP32, so a pure-FP64 resident GFN1 SCF is **not**
  necessarily faster than the 8-core MKL CPU on small/medium systems; Stage 2a is
  a correctness + residency milestone. The throughput win arrives with Stage 5
  (FP32 bulk) and device-resident `-opt`/`-md` (Stage 4, no host round-trip per
  step).
- Memory headroom (16 GB): `complex` (231 atoms) is ~tens of MB; the dense path's
  practical ceiling is ~`polymer` (1410 atoms, ~3–4 GB). Larger → large-system
  modes (CPU-orchestrated; the GPU dense kernel is the per-block solver).
- Large-system modes (`-large_system_mode fragments|dc|sparse`) are **not** wired
  to `-gpu` yet; with `-gpu` GFN1/Broyden they run the CPU fragment/DC driver.
