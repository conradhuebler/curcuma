# Work Package: Stage 3 — native xTB integrals on the GPU

> 🤖 AI-prepared work package (2026-06-03), not yet implemented. Companion to the
> staged GPU plan in `~/.claude/plans/ich-w-rde-gerne-eine-wiggly-scroll.md` and
> the status doc [docs/SQM_GPU.md](SQM_GPU.md). Stages 0–2 are committed
> (`6956f16`, `20bfb29`, `72c50e8`); this WP specifies Stage 3.

## 0. Where we are / what Stage 3 is

Stages 0–2 made the **SCF loop** device-resident for GFN1 (2a) and GFN2 (2b): H0,
S, the Cholesky factor L, the density P, MO coefficients C (and for GFN2 the
multipole integrals dp_int/qp_int) live on the GPU for the whole loop. **But the
integrals themselves are still built on the host** (`getHamiltonianH0`,
`computeCoordinationNumbers`, `buildGammaMatrix`, `setupMultipole`) and uploaded
once per geometry by `GpuScfBackend::begin` / `beginMultipole`.

Stage 3 moves that one-time-per-geometry integral build onto the GPU: CN, S, H0,
γ, and the GFN2 multipole integrals. After Stage 3, `begin()` receives the
geometry + basis + parameters and **computes** S/H0/L on the device instead of
receiving them; nothing nao²-sized is uploaded per geometry.

### Honest value assessment (read before prioritising)

- The integrals are built **once per geometry**, not per SCF iteration. For `-sp`
  they are setup cost (~49 ms threaded for `complex`/231; the SCF loop dominates).
- The real payoff is **`-opt`/`-md`**: every geometry step currently does a CPU
  integral build + a host→device upload of S/H0/L (+ 9 multipole matrices for
  GFN2). Stage 3 removes both per step.
- **BUT** `-opt`/`-md` are not device-resident until **Stage 4** (gradients):
  today the gradient runs on the CPU and pulls P/C/integrals back to the host
  every step. So Stage 3's throughput win is only fully realised together with
  Stage 4. Stage 4's dS/dR, dH0/dR kernels are *derivatives of the Stage-3
  integral kernels*, so Stage 3 is the natural prerequisite — but if the goal is
  "fastest path to a useful GPU -opt", consider whether Stage 3+4 should be
  planned as one push.

## 1. CPU code to mirror (exact references)

| Component | CPU entry | Kernel math |
|---|---|---|
| CN | `XTB::computeCoordinationNumbers` `xtb_h0.cpp:54`; `cn_exp` `xtb_params_extra.hpp:154` (GFN1), `cn_gfn` `xtb_params_extra.hpp:183` (GFN2) | O(N²) pairwise counting function (exp/erf damping) |
| self-energies | `XTB::getSelfEnergies` `xtb_h0.cpp:79` | `ε_sh(s) = selfenergy[s] − kcn[s]·CN(atom(s))` |
| S + H0 | `XTB::getHamiltonianH0` `xtb_h0.cpp:95-232`; `CGTO::cgto_overlap` `STO_CGTO.hpp:242` | shell-pair overlap × `avg_eps·h_factor` |
| γ | `XTB::buildGammaMatrix` `xtb_coulomb.cpp:21`; `coulomb::build_gamma_matrix` `xtb_coulomb.hpp:98` | O(nsh²) Klopman–Ohno |
| multipole ints | `XTB::setupMultipole` `xtb_multipole.cpp:68`; `multipole_ints::cgto_multipole` `xtb_multipole_ints.hpp:145`, `moment1d` `:47` | dipole(3)/quad(6) AO ints + origin shift + traceless transform |
| params | `xtb_params_extra.hpp` `namespace gfn1` (`:90`: `kpair :111`, `kshell :96`, `enscale`) / `gfn2` (`:134`: `kpair`, `kshell`, `wexp`); `pauling_en[]`, `atomic_rad_au()`; per-element tables in `parameters/gfn{1,2}_params.hpp` (hubbard, shell_hubbard, gexp, p_vcn, p_rad, p_dkernel, p_qkernel, mp_* constants) |

### S/H0 `h_factor` logic (the heart of the hard kernel)

For shell pair (a,b) with `avg_eps = ½(ε_sh(a)+ε_sh(b))`, `H0_μν = avg_eps·h_factor·S_μν`:
- **On-atom** (iat==jat): `h_factor = 1.0` (and the same-orbital diagonal S=1).
- **Off-atom**: `h_factor = hs · π_ij` where
  - `π_ij = (1 + shpoly[a]·rr)(1 + shpoly[b]·rr)`, `rr = sqrt(sqrt(R²)/(rad(zi)+rad(zj)))` (`atomic_rad_au`).
  - **GFN1** `hs`: valence flags (first shell of each ℓ per atom is valence,
    `xtb_h0.cpp:113-124`). `vi&&vj`: `kpair(zi,zj)·kshell(la,lb)·(1+enscale·ΔEN²)`;
    `vi^vj`: `0.5(kshell(l,l)+kdiff)`; else `kdiff`. `ΔEN = pauling_en[zi-1]−pauling_en[zj-1]`.
  - **GFN2** `hs`: `zij·km`, `km = kpair(zi,zj)·kshell(la,lb)·(1+enscale·ΔEN²)`,
    `zij = (2√(ζa·ζb)/(ζa+ζb))^wexp` (ζ = `cgto.slater_exp`).

### Overlap kernel facts (from `STO_CGTO.hpp`)

- **s/p only** (`ao_to_type` returns −1 for d; `cgto_overlap` only handles types
  0=s,1=px,2=py,3=pz). Bit-for-bit identical to the CPU coverage.
- p-ordering inside a shell is tblite's `{py=2, pz=3, px=1}` (`ao_to_type`
  `xtb_h0.cpp:41`). Preserve it.
- Primitive overlaps are **closed form** (no recursion): `primitive_ss/sp/pp_overlap`
  (`STO_CGTO.hpp:187-227`). Only `pow`/`exp` needed — both exist as CUDA device math.
- **The basis is already orthogonalised on the host** (`CGTO::orthogonalize`,
  `STO_CGTO.hpp:313`, appends extra primitives to 2nd same-ℓ shells). So the device
  just consumes the final `m_basis.cgto[ish]` primitives (alpha/coeff) — do NOT
  re-run orthogonalisation on the device; upload the post-ortho primitives.

## 2. Device design

### 2.1 Flattened basis (host → device, once per molecule)

`CGTO::Shell` uses `std::vector` (host-only). Flatten `m_basis` into device arrays
(upload once in `begin` / a new `beginBasis`):

```
int   nsh, nao, nat;
int   sh_ang[nsh], sh_nprim[nsh], sh_prim_off[nsh];   // prim_off = exclusive prefix sum of nprim
int   sh2at[nsh], iao_sh[nsh], nao_sh[nsh], ish_at[nat], nsh_at[nat];
int   ao2sh[nao], ao2at[nao];
double prim_alpha[Σnprim], prim_coeff[Σnprim];        // post-orthogonalisation
double sh_zeta[nsh];                                  // cgto.slater_exp (GFN2 zij)
double shpoly[nsh], selfenergy[nsh], kcn[nsh];        // from m_h0 (H0Data)
int    z[nat]; double xyz_bohr[3*nat];                // geometry per Calculation
```

Geometry (`xyz_bohr`) and the geometry-dependent self-energies change per step;
the rest is molecule-constant. Keep the split so `-opt`/`-md` re-upload only xyz.

### 2.2 `__device__` integral functions

Port verbatim into a device header (e.g. `cuda/xtb_gpu_integrals_device.cuh`):
- `d_primitive_ss/sp/pp_overlap` (trivial).
- `d_cgto_overlap(...)` operating on flattened prim arrays + an explicit `nprim`/
  `offset` instead of the `Shell` struct.
- For multipole: `d_moment1d` + `d_cgto_multipole` (port `xtb_multipole_ints.hpp`),
  then the origin-shift + traceless transform inline (mirror `setupMultipole`
  `xtb_multipole.cpp:120-159`).

### 2.3 Parameter tables as `__constant__`

Element tables (`hubbard_parameter[86]`, `shell_hubbard[86][3]`, `pauling_en[86]`,
`atomic_rad`, `p_vcn/p_rad/p_dkernel/p_qkernel[86]`, repulsion tables) → `__constant__`
device arrays, copied with `cudaMemcpyToSymbol` once. `kpair`/`kshell`/`enscale`/`wexp`/
`gexp`/`mp_*` are small scalars/tiny tables — pass as kernel args or `__constant__`.

### 2.4 Kernels

| Kernel | Grid | Notes |
|---|---|---|
| `k_cn` | 1 thread / atom, inner loop over atoms | mirror `cn_exp`/`cn_gfn`; reuse the GFN-FF `k_cn_compute` structure (`ff_methods/cuda/gfnff_kernels.cu`) |
| `k_self_energy` | 1 thread / shell | `se[s] = selfenergy[s] − kcn[s]·CN[sh2at[s]]` |
| `k_overlap_h0` | **one block per shell-pair** (or 1 thread/AO-pair) | compute shell-pair `h_factor` once, broadcast to AO pairs; write S and H0. Watch the on-atom same-orbital S=1 special case (`xtb_h0.cpp:214`) |
| `k_gamma` | 1 thread / shell-pair | Klopman–Ohno (`build_gamma_matrix`); diagonal = raw hardness, on-atom off-diag = averaged, cross-atom = damped |
| `k_multipole_ints` | one block per shell-pair (like overlap) | dipole+quad + origin shift + traceless transform; needs S (so run after `k_overlap_h0`) |

All column-major n×n outputs (S, H0 symmetric → layout-safe; dp_int/qp_int NOT
symmetric — match the CPU `(μ,ν)` indexing used by `k_add_fock_multipole` /
`k_multipole_moments` from Stage 2b).

### 2.5 Integration into the resident path

Rework `XtbGpuContext::residentBegin` (and `residentBeginMultipole`) to **compute**
rather than receive:
- New `beginBasis(flattened basis + params)` — once per molecule.
- `residentBegin(xyz, z)` per geometry: `k_cn` → `k_self_energy` → `k_overlap_h0`
  → **L = chol(S) on device via `cusolverDnDpotrf`** (replaces host Eigen LLT in
  `buildOrthonormalizer`; keep lower triangle, zero the upper or rely on the
  trsm fill mode) → `k_gamma`. S/H0/L/γ now resident; **nothing nao² uploaded**.
- GFN2: `residentBeginMultipole` runs `k_multipole_ints` (no upload).
- γ is currently consumed **host-side** (`addCoulombShellPotential` does
  `v_sh += γ·q_sh` on the host). Options: (a) download γ once per geometry (nsh²,
  cheap) so the host potential build is unchanged — simplest, recommended first;
  (b) move the Coulomb `gemv` to the device too (folds into the resident loop,
  later refinement). Same choice for CN: the host still needs CN for the in-SCF
  D4 (`addDispersionPotential`) and the multipole potential `mrad`/scalar shift —
  download CN (length nat) once per geometry.

> Net per-geometry transfers after Stage 3: **up** = xyz (3·nat) + the molecule-
> constant basis (once ever); **down** = γ (nsh²) + CN (nat) for the host
> consumers (until those move to the device too). The nao²-sized S/H0/L/dp/qp no
> longer cross the bus.

## 3. Decomposition (each sub-stage buildable + validated before the next)

The consumers of CN and γ are still on the host, so the *valuable* unit is "all
integrals on device + begin() computes". But implement and validate bottom-up:

- **3a — CN** (`k_cn`, `k_self_energy`): validate CN vs `computeCoordinationNumbers`
  elementwise @1e-9 (GFN1 `cn_exp` and GFN2 `cn_gfn`). Lowest risk; reuses GFN-FF
  kernel structure.
- **3b — S + H0** (`k_overlap_h0` + device `cgto_overlap` + `__constant__` params +
  device chol): the hard core. Validate S and H0 elementwise @1e-9 vs
  `getHamiltonianH0` on H2O/CH4/NH3/C6H6/HCN/caffeine/complex, GFN1 **and** GFN2
  (different `h_factor`). Then rewire `begin()` to compute S/H0/L on device and run
  the full `gpu_gfn1_validation` @1e-8.
- **3c — γ** (`k_gamma`): validate vs `build_gamma_matrix` @1e-9; download for the
  host potential (option a). Folds into GFN1 + GFN2 begin.
- **3d — multipole ints** (`k_multipole_ints`): validate dp_int/qp_int @1e-9 vs
  `setupMultipole`; rewire `residentBeginMultipole` to compute; full
  `gpu_gfn2_validation` @1e-8.

## 4. Validation strategy

- **Component (per sub-stage):** new ctest comparing the device integral to the CPU
  one elementwise @1e-9 (max-abs over the matrix). Add a small executable or extend
  the `sqm_*` harness; gate behind `if(USE_CUDA)`, runtime device-skip.
- **System:** `gpu_gfn{1,2}_validation` @1e-8 vs tblite (existing labels) must stay
  green with the documented xfails (He2, complex). `compute-sanitizer` on ≥1
  molecule per sub-stage (the overlap kernel has data-dependent inner loops — check
  for races/OOB).
- **No-CUDA + CPU regression:** rebuild canonical `release/`; CPU
  `gfn{1,2}_validation` + `native_xtb` must stay green (the host integral code is
  unchanged — Stage 3 only adds a device path behind `GpuScfBackend`).
- Reminder: GPU≠CPU bit-identical (GEMM/reduction order). The 1e-8 gate is vs
  tblite; component checks are @1e-9 elementwise (the integrals are not reductions,
  so they should match to ~1e-12, but allow 1e-9 for `pow`/`exp` ULP differences).

## 5. Risks / gotchas

- **`std::pow`/`std::exp` in device code:** use the CUDA device overloads (`pow`,
  `exp`); they differ from libm by a few ULP → keep component tolerance at 1e-9,
  not 1e-12. This is also why the converged SCF stays within the 1e-8 gate.
- **Orthogonalised basis:** upload the *post*-`orthogonalize` primitives; never
  re-orthogonalise on device (it would double-apply).
- **p-orbital ordering** `{py,pz,px}` and the on-atom **S=1 diagonal** special case
  (`xtb_h0.cpp:214`) must be reproduced exactly or the SCF shifts.
- **GFN2 traceless quadrupole transform** (`xtb_multipole.cpp:149-156`) and the
  per-column origin shift (origin = atom of the *column* AO) are easy to get wrong;
  mirror the CPU index-for-index. dp_int/qp_int are **not symmetric**.
- **Device Cholesky:** `cusolverDnDpotrf` writes one triangle; the Stage-1/2 trsm
  path uses `CUBLAS_FILL_MODE_LOWER` — make sure the lower L matches the host
  Eigen `matrixL()` convention (zero/garbage in the upper is fine for trsm with
  fill=lower).
- **Param tables `__constant__`:** 64 KB constant limit — the 86-element tables fit
  easily; the per-prim arrays are global, not constant.
- **Determinism:** integrals are pairwise (no cross-thread reduction), so the
  device result is deterministic and should match the CPU to ULP — easier to
  validate than the SCF.

## 6. Files

- New: `cuda/xtb_gpu_integrals_device.cuh` (device `cgto_overlap`/`cgto_multipole`/
  CN), kernels + `__constant__` tables in `cuda/xtb_gpu_context.cu` (or a new
  `cuda/xtb_gpu_integrals.cu` sharing the pimpl — note the Stage-2 buffers live in
  `xtb_gpu_context.cu`'s `Impl`, so adding there is simplest first).
- Extend: `cuda/xtb_gpu_context.{h,cu}` (`beginBasis`, compute-on-device
  `residentBegin`/`residentBeginMultipole`, device chol).
- Extend: `xtb_gpu_method.cpp` (flatten `m_basis`/`m_h0` and pass to `beginBasis`);
  the `GpuScfBackend` seam may need a `beginBasis(const BasisMap&, const H0Data&, …)`
  method (project types → keep core CUDA-free).
- `XTB::Calculation`: when the device computes the integrals, the host
  `getHamiltonianH0`/`buildGammaMatrix`/`setupMultipole`/`buildOrthonormalizer`
  calls can be skipped on the GPU-resident path (guard, don't delete — the CPU
  path stays the default and the validation reference).
- Tests: `test_cases/sqm_reference/CMakeLists.txt` — new component tests behind
  `if(USE_CUDA)`.

## 7. Suggested order for the next session

1. 3a CN (fast win, low risk, unblocks self-energies).
2. 3b S/H0 + device chol — the bulk; rewire GFN1 `begin()` and prove
   `gpu_gfn1_validation` 12/12.
3. 3c γ (download for host potential).
4. 3d multipole ints; prove `gpu_gfn2_validation` 12/12.
5. Update `docs/SQM_GPU.md` (Stage 3 → done), CLAUDE.md/README/AIChangelog, memory.

Then Stage 4 (gradients) becomes the path to fully device-resident `-opt`/`-md`,
reusing the Stage-3 integral kernels as the basis for dS/dR, dH0/dR.
