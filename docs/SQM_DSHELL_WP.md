# WP: native GFN1/GFN2 d-shell integrals (X-I1)

Status: **CPU implemented + machine-validated (2026-06-27)** — B1–B5 done, B0 guard removed.
B6 = GPU device d kernels: **CUDA + ROCm done**, Vulkan pending. 🤖 AI/machine-tested only,
NOT ✅ TESTED.

## What was done / tested (2026-06-27, Claude Generated)
- **B1** `xtb_ao_utils.hpp` (new): single source of truth for `as_cgto_shell`/`ao_to_type`
  (X-I5 resolved) + the cartesian(6)->spherical(5) `dtrafo` (verbatim from
  tblite `integral/trafo.f90`), general `cartMoment1d` (la,lb<=3) and the d overlap /
  multipole / gradient blocks. **Per-shell-pair rule:** only d-touching pairs use the
  transform path; pure s/p pairs keep the existing scalar kernels -> s/p molecules stay
  **bit-identical** (verified: H2O/CH4/NH3/C6H6/HCN/caffeine dE = 0 vs tblite).
- **B2** d overlap: `max|ΔS|` vs tblite **1.5e-10** on H2S/PH3/SiH4/HCl (gfn1+gfn2),
  `test_xtb_overlap` extended to the d block.
- **B3** d in H0 (GFN1): energy vs tblite **<=1e-8 Eh** (4/4 d molecules).
- **B4** GFN2 d multipole integrals (`setupMultipole` step 1b): energy vs tblite
  **<=1e-8 Eh** (4/4 d molecules).
- **B5** d gradients: `test_xtb_gradient` FD — GFN1 ~1e-6, GFN2 ~1-4e-4 Eh/Å (the same
  SCF/FD floor as the s/p GFN2 baseline H2O 1.7e-4 / NH3 5.0e-4), all PASS @5e-4;
  `-opt` on H2S converges.
- **B7** B0 guard removed (`xtb_native.cpp`); d molecules added to the `sqm_val_*` ctest
  suite (1e-8 gate); docs updated.

## GPU device d kernels (B6)
- **CUDA: implemented + validated (2026-06-27).** The device integral (overlap/H0/multipole)
  + SCF + gradient kernels handle d via the same dtrafo (per-spherical-AO-element device
  functions in `cuda/xtb_gpu_integrals_device.cuh`; 4 kernels in `xtb_gpu_context.cu`
  branched on `la>=2||lb>=2`). Gated by `GpuScfBackend::supportsDshell()` (CUDA overrides
  true). **Validated on a GTX 1660:** `-gpu cuda` energy **bit-identical to CPU** (H2S/PH3/
  SiH4/HCl, gfn1+gfn2), device gradient vs CPU **~1e-16** (`test_xtb_cuda_gradient`), `-opt`
  converges, s/p `gpu_gfn{1,2}_validation` 32/32 unchanged.
- **ROCm: implemented + validated (2026-06-28, B6 port).** The d-shell device helpers were
  mirrored verbatim from `cuda/xtb_gpu_integrals_device.cuh` into
  `rocm/xtb_hip_integrals.hiph`, and the same `la>=2||lb>=2` branch added to the three HIP
  kernels (`k_overlap_h0`/`k_multipole_ints`/`k_grad_h0_pulay` in `xtb_hip_context.hip`);
  `HipScfBackend::supportsDshell()` now returns true, so d systems run fully device-resident.
  **Validated on an AMD Radeon 890M (gfx1150):** `-gpu rocm` energy **bit-identical to CPU at
  8 dp** (H2S/PH3/SiH4/HCl, gfn1+gfn2; `sqm_val_rocm_*` 32/32 incl. s/p non-regression), and
  `-opt` on H2S (gfn2) tracks the CPU trajectory step-by-step in energy and gradient norm
  (gradient floor ~1e-6 from the default FP32 mixed-precision eigensolve; `-scf_threshold 1e-8`
  or `-scf_mixed_precision false` for tighter). s/p path byte-unchanged (the `!dpair` branch).
- **Vulkan: CPU fallback (device d kernels pending).** `supportsDshell()` stays false for the
  Vulkan backend; d systems use the validated CPU path. Porting needs GLSL/SPIR-V shader work
  (the HIP/CUDA `__device__` helpers do not carry over verbatim) — deferred.
- **Transition metals.** Enabled (guard removed) but **NOT validated** — no TM reference
  set; only main-group d (S/P/Cl/Si) was checked. Open-d SCF convergence is untested.
- Precision floor: the validate_sqm gate parses the 8-dp printed energy (~1e-8 resolution),
  same as the s/p suite.

---

## Original brief (for reference)

🤖 AI-authored brief, machine-context only. Nothing here is ✅ TESTED.

## Why this is its own session
The integral kernels (overlap, H0, GFN2 multipole, and all their gradients) are **s/p-only**.
GFN uses **spherical** harmonics; for s/p the cartesian and spherical bases coincide (1 and 3
functions), so the current code uses direct cartesian "type codes". For **d** there are 6
cartesian but 5 spherical functions, so d requires a cartesian(6)→spherical(5) transform
applied consistently to S, H0, the GFN2 dipole/quadrupole integrals, and every gradient.
Getting this wrong yields *plausible but wrong* energies — exactly the failure mode CLAUDE.md
warns about — so it must be implemented incrementally and validated to µEh vs tblite, not
folded into an unrelated batch.

## Current state (verified 2026-06-26)
- The basis builder already **emits** d shells: `XTB::buildBasis()`
  (`src/core/energy_calculators/qm_methods/xtb_native.cpp`) pushes `nao = 2*ang+1 = 5` AOs for
  any shell with `ang==2`. So the AO space is allocated; only the **integrals** over those AOs
  are missing (they become zero rows → singular overlap → eigensolver failure).
- Which elements have a d shell: `gfn2_params.hpp` / `gfn1_params.hpp` `ang_shell[86][3]` —
  **Na, Mg, Al, Si, P, S, Cl, Ar and every transition metal** (`{0,1,2}` or `{2,0,1}`).
  d Slater exponents exist (`slater_exponent[idx][2]`); **confirm** d-shell H0 self-energy /
  `kshell(l=2)` / level params are present (tblite carries them) before B3.
- **B0 guard to remove when validated:** `XTB::InitialiseMolecule()` in `xtb_native.cpp`, the
  loop `for sh: if (m_basis.ang_sh[sh] > 1) setHardError(...); return false;` (search "B0 safety
  scaffold"). Also the wrapper message in `native_xtb_method.cpp` is generic-surfaced; nothing
  else to revert.
- Empirical proof of the bug (pre-guard): `curcuma -sp h2s.xyz -method gfn2` → singular overlap
  → eigensolver fail → was `0.0 Eh` exit 0. Post-guard: clear error, exit 1.

## The s/p-only kernels to extend (~13 files)
CPU:
- `STO_CGTO.hpp` — `cgto_overlap` handles type ∈ {0=s,1=px,2=py,3=pz} only (comment at ~:361);
  `cgto_overlap_grad` returns 0 for `type > 3` (~:383). The 1-D `moment1d` machinery already
  takes arbitrary cartesian `lx,ly,lz`, so the d primitive integrals are reachable — the gap is
  the AO mapping + the spherical transform.
- `xtb_multipole_ints.hpp` — `type_to_cart` maps only 0..3; `cgto_multipole` (GFN2 dip/quad).
- `xtb_h0.cpp` — `ao_to_type(ang,local_ao)` returns −1 for `ang==2` (~:42); the S/H0 loops
  `if (t_a < 0) continue;` skip d (~:308).
- `xtb_multipole.cpp`, `xtb_response.cpp` (CPSCF H0/hscale + `ao_to_type` copies),
  `xtb_gradient.cpp` (`as_cgto_shell_g`/`ao_to_type_g` `_g`-suffixed copies; Pulay).

GPU (mirror the CPU kernels; **gate off for d in B6 first**, implement later):
- `cuda/xtb_gpu_context.cu`, `cuda/xtb_gpu_integrals_device.cuh`
- `rocm/xtb_hip_context.hip`, `rocm/xtb_hip_integrals.hiph`
- `vulkan/shaders/{overlap_h0,multipole_ints,grad_pulay}.comp`, `vulkan/vk_integrals.glsl`,
  `vulkan/xtb_vulkan_context.h`

Note the duplication debt (X-I5/X-G1): `ao_to_type` / `as_cgto_shell` are copied ~3×; consider
extracting a shared `xtb_ao_utils.hpp` as part of B1 so d is added in one place.

## Plan
- **B1 — AO mapping + spherical transform.** Define the 5 spherical-d AO ordering matching
  tblite, and the cartesian(6)→spherical(5) `dtrafo` matrix. Extend `ao_to_type`/`type_to_cart`
  (and the `_g` copies) to d. Single source of truth (shared header).
- **B2 — d overlap.** Extend `cgto_overlap` (STO_CGTO) to `l=2` via the cartesian
  Hermite/Obara-Saika `moment1d` path + `dtrafo`. Validate the **S matrix** vs tblite first
  (largest, easiest signal).
- **B3 — d in H0.** `kshell(l=2)`, d self-energy/level (confirm params), shpoly. Validate H0
  then the SCF energy vs tblite.
- **B4 — GFN2 multipole integrals for d.** `cgto_multipole` + consumers (`xtb_multipole.cpp`,
  `xtb_response.cpp`).
- **B5 — gradients for d.** `cgto_overlap_grad` (drop the `type>3 → 0`), H0/Pulay
  (`xtb_gradient.cpp`), multipole-integral Pulay. FD-validate.
- **B6 — GPU.** Initially detect a d shell in the GPU wrappers and fall back to CPU with a clear
  one-time message (do NOT silently run the broken device kernel). Implement device d as a
  follow-up.
- **B7 — validation + remove guard.** Compare native vs `tblite-gfn2` / `dump_tblite_reference`
  on **H2S, PH3, HCl, SiH4** (+ one transition-metal case) to ≤1e-6 Eh; FD-validate the
  d-containing gradients. Add ctest reference data alongside the native-GFN suite. **Only then**
  remove the B0 guard and update CLAUDE.md / TECHNICAL_DEBT.md (X-I1 → resolved).

## Validation harness
- Energies: `curcuma -sp <mol> -method gfn2` vs `-method tblite-gfn2` (needs `USE_TBLITE=ON`),
  or the `dump_tblite_reference` tool already used for solvation refs.
- Build/test in `release/` (`make -j4`, `ctest`). Keep all existing s/p molecules
  bit-identical (H2O/CH4/NH3/caffeine/complex) — d support must not perturb the s/p path.
- Suggested first targets: H2S, PH3, HCl, SiH4 (small, one heavy d-element each), then a TM.

## Cross-refs
- Audit row: `docs/TECHNICAL_DEBT.md` X-I1 (and the 2026-06-26 resolved block: B0 guard).
- Integral background: `src/core/energy_calculators/qm_methods/QM_ARCHITECTURE.md`.
