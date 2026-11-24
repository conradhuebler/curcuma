# ROCm parity N: cause of the ROCm SCF run-to-run differences (Nov 24, 2026)

Item N of [ROCM_PARITY_PLAN_2026-11.md](ROCM_PARITY_PLAN_2026-11.md). 🤖 AI-generated, ⚙️ machine-tested. Machine: one
RX 9060 XT (gfx1200), ROCm 7.2.4, Linux 6.18.55. Worktree HEAD `ebd1b357` (feature/rocm-parity), `release/` built with
`-DUSE_ROCM=ON -DUSE_MKL=ON -DUSE_PCH=ON -DUSE_AVX2=ON`.

## Result

The differences reported in [ROCM_PARITY_A3.md](ROCM_PARITY_A3.md) ("Open") are reproduced by **a second process
using rocSOLVER on the same GPU at the same time**, and not by a preceding CPU calculation in the process. Which
concurrent process was on the GPU during the A3 runs is not known (not recorded). With the GPU to itself, every ROCm series below is bitwise identical, with or without a preceding CPU calculation; with one concurrent
GPU process running the SCF (or only `rocsolver_dsygvd`), the same series differ, and `test_pbc_xtb_gpu_parity rocm`
fails on exactly the A3 cases (216 H2O box, gfn1 and gfn2). This is the rocSOLVER concurrency fault already described in
[ROCM_GFN_POTENTIAL_2026-10-14.md](ROCM_GFN_POTENTIAL_2026-10-14.md#wrong-rocsolver-results-with-several-processes-on-one-gpu-oct-1617-2026):
`rocsolver_dpotrf`/`dsygvd`/`ssygvd` return wrong results (max |diff| 3e-4 to 1.3e-3 on the factor or eigenvectors),
and the converged energies of the affected curcuma runs differ by 1e-14 to 4e-12 Eh (that the SCF iterates a wrong
early step away is the likely reading, not traced per iteration).

## Measurements

Reproducer: scratch program (`EnergyCalculator`, registry molecules, `scf_threshold` 1e-9, verbosity 0), one optional
calculation before the series ("pre"), then n calculations from fresh calculators in the same process; hash = FNV-1a
over the energy and gradient bytes (plus virial for the cell). Every run under the shared GPU lock; "load" = a second
process started inside the same lock hold (repro caffeine gfn2 in a loop, or the routine named).

| Variable changed (one per row) | Cell | n | distinct hashes |
|---|---|---|---|
| none (baseline) | caffeine gfn1, pre cpu, threads 1 | 4 proc x 10 | 1 |
| pre: none instead of cpu | caffeine gfn1, threads 1 | 4 x 10 | 1 |
| threads 4 (CPU pre and GPU series) | caffeine gfn1 | 4 x 10 | 1 (same hash as threads 1) |
| method gfn2 | caffeine, pre cpu | 4 x 10 | 1 |
| cell, pre cpu with 4 threads, GPU series 1 thread (A3 probe setup) | 216 H2O, L 18.63 A, gfn1 | 1 x 10 | 1 |
| CPU load: 12 busy loops on 10 cores | caffeine gfn1 / gfn2 / 216 H2O | 3x10 / 3x10 / 1x10 | 1 / 1 / 1 |
| separate `curcuma -sp -gpu rocm` processes, md5 of `-dump_gradient` | caffeine gfn1 t1 / gfn1 t4 / gfn2 t1 | 10 / 10 / 10 | 1 / 1 / 1 |
| `test_pbc_xtb_gpu_parity rocm`, GPU to itself | full test | 3 runs | 0 DIFFERS, 3/3 pass |
| **one concurrent GPU process (gfn2 SCF loop)** | caffeine gfn1, pre cpu | 20 | **5** (16 + 4 singles; dE up to 1.4e-14 Eh) |
| same, repeat | caffeine gfn1, pre cpu / pre none | 20 / 20 | 2 / 1 |
| same, FP64 only (`scf_mixed_precision false`) | caffeine gfn1 | 60 | 2 |
| **one concurrent GPU process** | 216 H2O cell, A3 probe setup | 10 | **9** (E spread 3.8e-12 Eh) |
| **one concurrent GPU process** | `test_pbc_xtb_gpu_parity rocm` | 2 runs | **2/2 fail**: 216 H2O gfn1 and gfn2 DIFFERS |
| the load processes themselves (caffeine gfn2 loop, no other variable) | while the measured job ran | 36 to 449 runs per cell (7 cells) | 2 to 38 (a 16-run cell: 1) |
| concurrent load = GFN-FF on ROCm (no rocSOLVER) | caffeine gfn1 measured / GFN-FF load | 20 / 338 | 1 / 1 |

Library level, no curcuma (`bitcheck`: same input every iteration, count of outputs not bitwise equal to iteration 0,
n = 66, the caffeine GFN2 basis size):

| Measured routine | load: none | load: `dsygvd` loop | load: `dgemm` loop | load: own FMA kernel |
|---|---|---|---|---|
| `rocsolver_dpotrf` | 0 / 2000, 0 / 20000 | 7 / 2000, 1 / 2000, **39 / 20000** (max 7.0e-4) | 0 / 2000, 0 / 20000 | 0 / 2000, 0 / 20000 |
| `rocsolver_dsygvd` | 0 / 2000 (twice) | 3 / 2000, 1 / 2000 (max 1.3e-3) | 0 / 2000 | 0 / 2000 |
| `rocsolver_ssygvd` | 0 / 2000 | 3 / 2000 (max 9.4e-4) | | |
| `rocsolver_dsyevd`, `ssyevd` | | 0 / 2000 each | | |
| `rocblas_dgemm` | | 0 / 2000 | | |
| own FMA kernel | 0 / 2000 | 0 / 2000 | | |

So the fault needs rocSOLVER (here `sygvd`) in the *other* process; a concurrent `dgemm` or plain kernel load does not
trigger it (0 of 20000). The `syevd` zeros (n = 2000) do not exclude a lower rate.

## Excluded (each with the numbers above)

- Preceding CPU calculation in the process (pre cpu vs none, solo): 1 hash each.
- FP control state: MXCSR 0x1f80 (FTZ 0, DAZ 0, round-to-nearest) and x87 CW 0x037f at start, after the CPU
  calculation and at the end, in all 31 reproducer processes.
- Host thread count (`-threads` 1 vs 4, also CPU pre 4 / GPU 1): same hash; a CPU part of the ROCm run that depends on
  the thread count would change it.
- CPU load / scheduling of host threads (12 busy loops): 1 hash.
- Process-level state (separate processes, md5 of `-dump_gradient`): 1 hash per cell, 3 cells x 10.
- Unsynchronised device-to-host copies: every `hipMemcpyAsync(..., DeviceToHost, stream)` in `xtb_hip_context.hip` and
  `xtb_hip_periodic_host.hiph` is followed by a `hipStreamSynchronize` before the host reads (read, not measured).
- Mixed precision: the effect stays with FP64 only.

Not separated: rocBLAS atomics mode (`rocblas_set_atomics_mode`) was not switched; the dgemm load never triggered a
fault and the measured `dgemm` never differed, so it was not pursued. Resident vs host-driven loop, cached Cholesky,
explicit Linv and `rocm_density_syrk` under load were not measured (all paths call `potrf` or `sygvd`, which both fail
under load). `-xtb.rocm_fixed_point false` not measured.

## Change

None (no source or test file changed). Proposed, not applied: `RESOURCE_LOCK "gpu"` for the 32 `sqm_val_rocm_*`
tests in `test_cases/sqm_reference/CMakeLists.txt`, the only ROCm-labelled tests without it (32 of 128), so that
`ctest -j` cannot run them next to other GPU tests. It was not applied because its gate could not be met: without the
lock, `ctest -j8 -R "sqm_val_rocm|pbc_xtb_gpu_parity_rocm"` passed 2 of 2 times (0 DIFFERS), i.e. the short sqm tests
did not overlap the 216 H2O cases in that run, so no before/after difference was measurable. The A3 failures therefore
stay unattributed to a specific concurrent process; what is shown is that a concurrent rocSOLVER process reproduces
them (2 of 2) and that none occurs without one (3 of 3 passes, every in-process series 1 hash).

## Limits

- The lock only serialises tests inside one `ctest` call. Two `ctest` runs, an agent's run outside the shared
  `flock`, or a user job on the same GPU still produce the fault; with four processes the Oct 2026 measurement also
  found converged energies off by 1e-8 to 9e-3 Eh. A GPU job per process stays the only safe use on this card.
- The CUDA `sqm_val_gpu_*` tests have no `RESOURCE_LOCK` either; CUDA was not tested for this fault, left unchanged.
- Cause inside rocSOLVER, the driver or the hardware not separated; one card, one ROCm version; no report to AMD.

Reproducer sources (scratch, not in the repository): `repro.cpp` (curcuma series) and `bitcheck.cpp` (library level,
build `hipcc -O2 --offload-arch=gfx1200 bitcheck.cpp -lrocsolver -lrocblas`) in the agent N scratch folder.

## Not tested

Other GPUs, other ROCm versions, n other than 66 in `bitcheck`, more than one load process at the library level,
GFN-FF ROCm measured under a rocSOLVER load (it uses `potrf`/`getrf` only in some EEQ paths).

## Proposed lines (central files, for the coordinator)

- TODO.md, replace the A3 line "ROCm: SCF nicht reproduzierbar im selben Prozess nach einer CPU-Rechnung ..." by: "ROCm:
  bitweise Abweichungen der SCF (1e-14 bis 4e-12 Eh) entstehen nur, wenn ein zweiter Prozess gleichzeitig rocSOLVER
  auf derselben GPU nutzt (nicht durch eine vorige CPU-Rechnung; docs/ROCM_PARITY_N.md); gehört zum bestehenden Eintrag
  'rocSOLVER liefert bei mehreren Prozessen ... falsche Ergebnisse'."
- TODO.md, existing rocSOLVER entry, append: "Nov 24: auch `ssygvd`; tritt nur auf, wenn der andere Prozess rocSOLVER
  rechnet (potrf 39/20000 bei sygvd-Last, 0/20000 bei dgemm- oder Kernel-Last), docs/ROCM_PARITY_N.md."
- AIChangelog.md: "2026-11-24 ROCm parity N: the in-process ROCm SCF differences of A3 are reproduced by the rocSOLVER
  multi-process fault (solo 1 hash in every cell, one concurrent rocSOLVER process 9 of 10 distinct on 216 H2O) and
  not by a preceding CPU calculation -> docs/ROCM_PARITY_N.md"
- test_cases/sqm_reference/CMakeLists.txt (coordinator decision): `RESOURCE_LOCK "gpu"` for `sqm_val_rocm_*`.
- GPU_REPRODUCIBILITY.md "Not covered" (ROCm): "rocSOLVER results are not reproducible while another process uses
  rocSOLVER on the same GPU (docs/ROCM_PARITY_N.md)".
- GPU_PERIODIC_GAPS.md B4 "open: in-process reproducibility" -> "closed: not a periodic-path defect, concurrent
  rocSOLVER (docs/ROCM_PARITY_N.md)".

## Plugin ABI

Change needed: no (no source change).
