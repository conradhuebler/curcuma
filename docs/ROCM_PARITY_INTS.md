# ROCm INTS: shell-pair kernel for the periodic S/H0 build (Nov 25, 2026)

🤖 AI-generated, ⚙️ machine-tested on one RX 9060 XT (gfx1200). Human production testing pending.
Worktree start `c781caa7` (branch `feature/rocm-gfn1-ints`); changes uncommitted.

## What was changed

Follow-up of G1b ([ROCM_PARITY_G1B.md](ROCM_PARITY_G1B.md)): on this GPU the AO-pair kernel `k_overlap_h0_pbc` made the
device build of the periodic GFN1 integrals slower than the host build for many images (216 H2O, 125 images).

- `rocm/xtb_hip_periodic_kernels.hiph`: new kernels `pbc_dev::k_overlap_h0_pbc_sh` and `pbc_dev::k_pbc_kind_pairs`;
  `k_overlap_h0_pbc` is unchanged (ISA of all 150 former functions of `xtb_hip_context.hip.o` identical to the base).
  - One work item per shell pair (a, b); the (2la+1)(2lb+1) AO elements share d0, the image loop with its cutoff test,
    rr, pi_ij, hs and, per primitive pair and image, g, P and S00. Templates per (la, lb) in {0, 1}^2; a d shell writes 0
    as the AO-pair kernel does (a cell refuses d shells before this point).
  - Work items walk the shells in an l-sorted order (stable), so a wave of 32 mostly runs one (la, lb) branch.
  - `pow(pi/g, 1.5)` and `-ai aj / g` are tabulated once per pair of shell kinds (`k_pbc_kind_pairs`; a kind = shells
    with the same primitive count and bitwise equal exponents; more than 256 kinds: computed inline as before).
  - If the translation list is the `PBCUtils::makeLatticeImageList` grid (checked on the host, bitwise `T = n * L`),
    only the lattice cells that can lie within the cutoff per axis are visited, in list order; otherwise the whole
    list is walked. The exact r2 test is unchanged.
  - Element arithmetic and summation order (primitive pairs i, j inside an image, images in list order, on-atom block
    first) are those of `k_overlap_h0_pbc` / `d_cgto_overlap`: S and H0 are bitwise equal (measured below).
- `rocm/xtb_hip_periodic_host.hiph`, `XtbHipContext::computeIntegralsPeriodic`: launches the shell-pair kernel by
  default (block 32 x 8); `CURCUMA_PBC_GPU_INTS_LEGACY=1` (any value not starting with `0`) launches the G1b AO-pair
  kernel. Builds the shell order, kind map and grid check on the host (three small device-to-host copies of basis
  arrays per call). `PbcState` gained `dPbcIntShOrd`, `dPbcIntShKind`, `dPbcIntKindRep`, `dPbcIntKindTab`.
- `k_multipole_ints_pbc` (GFN2 dipole/quadrupole integrals) stays one thread per AO pair.

## What was measured

Build `release/` of the worktree (`-DUSE_ROCM=ON -DUSE_MKL=ON -DUSE_PCH=ON -DUSE_AVX2=ON`); base = plugin built from
`c781caa7` before any edit (`CURCUMA_PLUGIN_DIR`, same core). Every GPU run under `flock /tmp/claude-1000/rocm-gpu.lock`,
one job at a time; `-sp -scf_threshold 1e-9 -dump_gradient`. Cells (orthorhombic only; the native periodic xTB refuses
other cells): water dimer 8 A (`wd8`) and 8 x 10 x 12 A (`wdnc`), caffeine 14 A (`caf14`), 14 x 16 x 18 A (`cafnc`),
30 A (`caf30`), 216 H2O 18.6258 A (`w216`, 4 threads, others 1 thread).

| Gate | Result |
|---|---|
| S/H0 element by element, shell-pair vs AO-pair kernel in the same process (temporary diagnostic, removed), 6 cells x gfn1/gfn2 | 0 differing bits in all 12 cases (up to 2985984 elements, w216 gfn1) |
| Same with the full-list walk and without the kind table forced (temporary switches, removed), wd8/caf14/w216 x gfn1/gfn2 | dumps md5-identical to the default path, 18/18 |
| ROCm vs CPU, CLI, 6 cells x gfn1/gfn2 | worst dE 4.8e-12 Eh (w216 gfn2), dG 3.0e-10 Eh/A (cafnc gfn2) |
| Shell-pair vs legacy switch, 12 cases | dumps md5-identical 12/12 |
| Legacy switch vs base plugin, 12 cases | dumps md5-identical 12/12 |
| Second run, separate process, 12 cases | md5-identical 12/12 |
| No cell (caffeine, water dimer x gfn1/gfn2), new vs base | md5-identical 4/4 |
| `pbc_xtb_gpu_parity_rocm` (10 cases), default and legacy switch | passed; worst dE 4.8e-12 Eh, dG 2.6e-10 Eh/A, dXi 5.1e-10 Eh; 10/10 second run bitwise; both outputs line-identical |
| ctest `pbc_xtb_gpu_fallback_rocm`, `pbc_xtb_gpu_combinations_rocm`, `pbc_xtb_gradient`, `pbc_xtb_virial`, `pbc_xtb_cell_refusals` | passed |
| ISA of `xtb_hip_context.hip.o`, base vs new (scratch `a2/isa.py`) | 150/150 former functions identical, 2 new kernels |

No FD check was run: the integrals and every dump are bitwise equal to the G1b kernel, whose FD gates are in
ROCM_PARITY_G1B.md.

Kernel time (rocprofv3 kernel trace, w216 gfn1, n = 1): `k_overlap_h0_pbc` 620 ms, `k_overlap_h0_pbc_sh` 153 ms,
`k_pbc_kind_pairs` 0.01 ms. Intermediate steps of the same trace: shell pairs with a per-thread primitive table 325 ms,
plus the l-sorted order 208 ms.

Timing, `CURCUMA_GPU_PROFILE=1 -verbosity 3`, n = 3 per variant, variants alternating; median [range]. host = G1b host
build (`CURCUMA_PBC_GPU_HOST=ints`), leg = `CURCUMA_PBC_GPU_INTS_LEGACY=1`, new = default. The "overlap + H0" bucket of
the device variants includes the device Cholesky (and for GFN2 the multipole kernel); that of the host variant does not.

| Cell | Method | nao / images | overlap + H0 ms host / leg / new | wall s host / leg / new |
|---|---|---|---|---|
| w216 | gfn1 | 1728 / 125 | 273 [267-274] / 568 [566-568] / 190 [189-191] | 4.43 [4.42-4.52] / 4.67 [4.67-4.68] / 4.29 [4.28-4.29] |
| w216 | gfn2 | 1296 / 27 | 48 [44-48] / 161 [160-162] / 116 [116-117] | 3.15 [3.14-3.15] / 2.99 [2.98-2.99] / 2.94 [2.93-2.94] |
| caf14 | gfn1 | 76 / 125 | 3 / 7 / 5 | 0.39 / 0.43 / 0.43 |
| caf14 | gfn2 | 66 / 125 | 1 / 5 / 4 | 0.42 / 0.42 / 0.42 |
| cafnc | gfn1 | 76 / 125 | 2 / 6 / 4 | 0.40 / 0.44 / 0.43 |
| cafnc | gfn2 | 66 / 125 | 1 / 5 / 4 | 0.43 / 0.42 / 0.42 |
| caf30 | gfn1 | 76 / 27 | 2 / 5 / 4 | 0.39 / 0.43 / 0.43 |
| caf30 | gfn2 | 66 / 27 | 1 / 4 / 4 | 0.42 / 0.42 / 0.42 |
| wd8 | gfn1 | 16 / 125 | 1 / 14 / 10 | 0.35 / 0.40 / 0.40 |
| wd8 | gfn2 | 12 / 125 | 1 / 8 / 6 | 0.39 / 0.39 / 0.40 |

Ranges of the small cells are within 0.01 s (wall) and 5 ms (bucket) of the median. The setup bucket of w216 gfn1:
host 405, leg 616, new 235 ms. On every cell the new kernel is no slower than the legacy one. On the small cells the
device variants (leg and new alike) are 0.03-0.05 s slower in wall time than the host build for gfn1; this was already
so with the G1b kernel and is not in the integral bucket (not traced further).

## Not tested

- Any other GPU than one gfx1200; CUDA is not touched (the CUDA kernel stays per AO pair).
- Cells larger than 216 H2O; elements other than H/C/N/O; shells with more primitives than the GFN1/GFN2 H/C/N/O basis.
- The inline fallback for more than 256 shell kinds (only forced with a temporary switch, see above).
- MD/opt/NPT with a cell (repeated geometry updates through the new kernel); `-gpu_sparse_integrals` with a cell.
- Thread counts other than 1 and 4.

## Not implemented

- The GFN2 multipole kernel `k_multipole_ints_pbc` keeps one thread per AO pair (w216 gfn2: overlap + H0 bucket
  116 ms device vs 48 ms host build of S/H0 alone; wall time favours the device).
- The CUDA `k_overlap_h0_pbc` has no shell-pair counterpart.
- No CLI parameter for the kernel choice; only the environment switch.

## Proposed lines (central files, for the coordinator)

- TODO.md: replace the G1b line "GFN1 mit vielen Bildern ... langsamer als der Host-Aufbau" by "ROCm periodisch: S/H0
  per Schalenpaar-Kernel (docs/ROCM_PARITY_INTS.md), 216 H2O GFN1 Setup 235 statt 405 ms (Host); offen: GFN2-Multipol-
  Kernel pro AO-Paar, CUDA-Gegenstück".
- AIChangelog.md: "2026-11-25 ROCm INTS: shell-pair kernel `k_overlap_h0_pbc_sh` for the periodic S/H0 (l-sorted
  waves, kind-pair primitive table, lattice-range image walk), bitwise equal to the AO-pair kernel; 216 H2O GFN1 kernel
  620 -> 153 ms, wall 4.67 -> 4.29 s (host 4.43 s); `CURCUMA_PBC_GPU_INTS_LEGACY=1` restores the G1b kernel ->
  docs/ROCM_PARITY_INTS.md"
- GPU_PARITY_CUDA_ROCM.md / `gpu_parity_known.json`: `k_overlap_h0_pbc_sh`, `k_pbc_kind_pairs` ROCm-only (performance
  variant, same values as `k_overlap_h0_pbc`).
- GPU_PERIODIC_GAPS.md B4: "ROCm G1b: S/H0 by shell pair since Nov 25, 2026 (docs/ROCM_PARITY_INTS.md)".

## Plugin ABI change needed

No. Only plugin-internal code changed (`XtbHipContext::computeIntegralsPeriodic`, `PbcState`, new kernels); no
interface between core and plugin was touched.
