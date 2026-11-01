# Delta-ML labels on the GPU: measured, not implemented (Oct 27, 2026)

Status: AI-generated, machine-tested measurement. No source change. Human production testing pending.

Question: should `-activelearn` compute its base and reference labels (`LabelThread`,
`src/capabilities/activelearn.cpp`) on the GPU? Decision rule set beforehand: implement an
opt-in `label_gpu` only if the best GPU configuration reaches >= 1.5x the best CPU configuration
in structures per minute on at least one molecule, with label differences below 1e-6 Eh and
1e-5 Eh/Angstrom.

## Setup

- Branch `feature/delta-ml`, HEAD `be71113b`, build `build_deltaml_cuda` (CUDA plugin, same
  `activelearn.cpp` as HEAD). AMD Ryzen 9 9950X3D (16 cores, 32 SMT threads), one RTX 5080
  (16 GB). No other jobs; load average 0.2 before, 4.1 after the matrix.
- Inputs: peo8 (58 atoms, 32 frames of a 6.5 ps GFN-FF MD at 400 K) and the 90-atom host
  (`test_cases/structures/ensembles/aaa-host-frames353.xyz`, 32 evenly spaced frames). Each set
  repeated 4x to 64x so that the label time is long enough to read.
- Runs: the real label path, `-activelearn in.xyz -method gfnff -ref_method gfn2 -threads N
  [-gpu cuda]`; the time is activelearn's own "labelled N structures in X s" (base and reference,
  energy and gradient each, fresh calculators per structure). "Threads" = N in-process label
  workers; "processes" = P separate curcuma processes with `-threads 1`, each on 1/P of the frames.
- Finding during the setup: `-gpu cuda` already reaches the label calculators today. `LabelThread`
  copies the full controller, which carries the global `gpu` key, so `-activelearn -gpu cuda`
  (meant for the descriptor/fit, [DELTA_ML_GPU_2026-10.md](DELTA_ML_GPU_2026-10.md)) also labels on
  the GPU, all workers on one device without a device lease. Confirmed by label differences (below)
  and 0.6-0.9 GB device memory per worker.

## Throughput (structures per minute, gfnff + gfn2 with gradients; two repeats, identical within 3 %)

| Configuration | peo8 (58) | host (90) | GPU memory |
|---|---|---|---|
| CPU, 8 threads | 15 750 | 4 900 | - |
| CPU, 16 threads | 29 260 | 6 000 | - |
| CPU, 32 threads (best CPU) | 32 340 | 6 000 | - |
| GPU, 1 worker | 1 500 | 1 105 | 0.93 GB |
| GPU, 2 / 4 / 8 threads, one process (best GPU) | 2 000 | 1 500 | 1.3 / 2.0 / 3.4 GB |
| GPU, 2 / 4 / 8 processes | 1 460 / 1 260 / 1 230 | 1 060 / 960 / 910 | 1.9 / 3.7 / 7.2 GB |
| Best GPU / best CPU | 0.06x | 0.25x | |

Base label alone (gfnff for base and reference): CPU 1 thread 25 600 (peo8) and 6 000 (host),
CPU 32 threads 137 000-154 000 and 77 000-82 000; GPU 1 worker 6 000 and 4 800, GPU 8 threads
6 000 and 6 000. GFN-FF on the GPU is slower than one CPU core at both sizes.

One GFN2 label costs about 28 ms on one CPU core at 58 atoms; a fresh GPU calculator per structure
costs about 40 ms at 58 and 54 ms at 90 atoms, so device setup and launch latency dominate. The
host saturates the CPU at 16 threads (16 physical cores).

## Label differences, GPU vs CPU (256 structures each)

| | base dE (Eh) | base dG (Eh/A) | ref dE (Eh) | ref dG (Eh/A) |
|---|---|---|---|---|
| peo8 | 1.3e-13 | 4.7e-14 | 1.5e-11 | 1.7e-06 |
| host | 2.0e-13 | 5.7e-14 | 3.0e-11 | 3.7e-06 |

GPU labels are run-to-run identical (two repeats, 0.0) and independent of 1 vs 8 workers. The
GFN2 gradient difference (<= 3.7e-6 Eh/A) is below the 1e-5 gate and three orders below the
held-out Delta-ML gradient errors (about 1e-2 Eh/A) and the fit's `sigma_force` (1e-3); it does
not matter for the fit.

## Decision

Accuracy passes, throughput fails by a factor of 25 (peo8) and 6 (host) against the 1.5x gate.
No `label_gpu` switch was implemented. Consequence for users: `-activelearn -gpu cuda` today makes
labelling 4x (host) to 16x (peo8) slower than the CPU labels of the same run with `-threads 32`;
whether the labels should be pinned to the CPU when `-gpu` is given is an open item in `TODO.md`.

## Not tested

ROCm and Vulkan; multi-GPU; molecules above 90 atoms (the GPU break-even for single structures
lies in the hundreds of atoms, [GPU_TUNING.md](GPU_TUNING.md)); ORCA/DFT references (no GPU path,
stay CPU); a warm per-worker GPU calculator (ruled out by the fresh-calculator rule in
`LabelThread`). Scripts: scratch only (`bench.py`, `cmpds.py`), not committed.
