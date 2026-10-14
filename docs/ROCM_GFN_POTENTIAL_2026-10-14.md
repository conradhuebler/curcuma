# ROCm GFN1/GFN2: where the time goes and what could be gained (Oct 14, 2026)

Status: measurement and analysis, AI-generated, no code changed. Machine: RX 9060 XT (gfx1200, FP64 about 1/17 of FP32),
Core Ultra 5 225, `-threads 10`, branch `feature/rocm-step-pauses-mirror`. All wall times are single runs (n = 1) unless noted.

## CPU against ROCm, `-sp`, same energies to 8 decimals

| System | Atoms | gfn2 CPU / ROCm | gfn1 CPU / ROCm |
|---|---|---|---|
| caffeine | 24 | 0.02 / 0.33 s | 0.02 / 0.31 s |
| urea-water-cluster197 | 197 | 0.44 / 1.55 s | 0.31 / 1.02 s |
| pep20 | 309 | 0.80 / 1.81 s | 0.75 / 1.54 s |
| water300 | 300 | 0.35 / 0.93 s | 0.40 / 1.02 s |
| water1002 | 1002 | 4.51 / 4.42 s | 5.90 / 7.17 s |
| peo201-chain | 1410 | 18.98 / 15.28 s | 24.78 / 19.06 s |

`-gpu auto` selects ROCm for GFN1/GFN2 at every size (checked on caffeine: "using native xTB on GPU (rocm)"); the
`gfnff.gpu_auto_min_atoms` threshold covers GFN-FF only. Break-even is near 1000 atoms for gfn2, above 1000 for gfn1.
ROCm `-opt` of a 114-atom and a 66-atom molecule is 3x slower than the CPU (n = 3, `ROCM_EVALUATION_2026-10-13.md`).

## Where the ROCm time goes (gfn2, peo201-chain, 1410 atoms, nao 3222, 15.0 s)

SCF 12.4 s (82 %), of which the "solve eigen" timer shows 11.7 s: 7 FP32 iterations at 0.64 s and 5 FP64 iterations at
1.55 s (switch at `scf_fp32_threshold` 1e-3); the solver itself is 9.7 s of that (7 x 428 + 5 x 1337 ms), the remaining
2.0 s is the FP64 density GEMM (see below). Setup 1.2 s (`k_overlap_h0` 0.46 s), D4 ATM 3-body 1.2 s (one FP64 kernel,
`k_d4_atm`). peo201-water500 (2910 atoms, nao 6222, 79 s): "solve eigen" timer 86 % (FP32 3.4 s and FP64 9.4 s per iteration, density GEMM included; at standalone dgemm speed 1.1 s per iteration it would be about 13 s of the 68 s, estimate).
rocprofv3 on the 1410-atom run: GPU busy 13.0 of about 15 s; FP64 GEMM kernels (inside sygst/syevd) 5.5 s,
`latrd` (tridiagonalisation, gemv) 1.8 s.

## Standalone rocSOLVER timing (subagent, median of 3, ms, FP32 / FP64)

| n | sygvd | syevd | sygst | potrf + trsm | gemm TFLOPS |
|---|---|---|---|---|---|
| 3222 | 410 / 1327 | 240 / 966 | 148 / 247 | 23 / 149 | 6.8 / 0.43 |
| 6222 | 1992 / 7998 | 1534 / 6212 | 365 / 994 | 88 / 885 | 7.2 / 0.43 |

In curcuma the solver call is as fast as standalone: temporary timers in `resSolveMP` (not committed) gave 416 to 435 ms
(ssygvd, hipEvent 417 to 426 ms, n = 7 calls) and 1330 to 1396 ms (dsygvd, n = 5 calls) at n = 3222, so the earlier
impression of "1.2 to 1.7 times slower" came from comparing the SCF iteration time (640 / 1545 ms) and the
"solve eigen" timer of `-verbosity 3` (11.7 s = 7 x 428 + 5 x 1337 ms + about 2.0 s) with the standalone number.
The difference is the density matrix: `resDensity` computes P = (C diag(occ)) C^T with a full n x n x n `rocblas_dgemm`
(`xtb_hip_context.hip`, `resDensity`), although columns beyond `ncol` are zero. rocprofv3 shows one FP64 GEMM kernel of
165 ms per iteration (n = 3222; 0.43 TFLOPS peak FP64 gives 2 n^3 / 0.43e12 = 160 ms), and thread count does not matter
(`-threads` 1 / 4 / 10: 710 / 655 / 660 ms per FP32 iteration, n = 1). Jacobi, QR and bisection variants are 10 to 50 times slower than
syevd (`syevdj` 80 s at n = 6222). FP32 purification (2 GEMM + geam per iteration, 17.3 ms at n = 3222, 120 ms at 6222)
only beats sygvd below about 10 to 16 iterations (arithmetic from the table, not a converging purification).

## Density GEMM (new finding, same day)

P is built with a full n^3 FP64 GEMM. Using only the `ncol` occupied columns and a symmetric rank-k update
(`rocblas_dsyrk` on C diag(sqrt(occ)), valid because occupations are non-negative) would cut the flops by about a factor 4
(E: 165 -> about 40 ms per iteration at n = 3222, 12 iterations: about 1.5 s of the 15 s run; at n = 6222 about 10 s of
79 s). Not implemented; results change at rounding level only. In FP32 iterations the same product could run in FP32.

## Levers (estimates marked E are arithmetic on the numbers above, not measured in curcuma)

0. **Density GEMM** (previous section).
1. **Precomputed L instead of potrf + sygst per iteration.** The device already holds L (S is constant per geometry);
   F' = L^-1 F L^-T needs two trsm (7 ms each FP32, 85 ms FP64 at n = 3222). E: FP32 solve 410 -> 265 ms (-35 %),
   FP64 1327 -> 1190 ms (-10 %).
2. **Later FP64 switch.** Measured on peo201-chain: `-scf_fp32_threshold` 1e-3 / 1e-4 / 1e-5 / 1e-6 gives 15.1 / 13.3 / 12.2 /
   15.0 s wall (12 / 12 / 13 / 17 iterations), energy -2088.25340678 Eh in all four (n = 1 structure). 1e-6 stalls in
   FP32 and loses. `TODO.md` asks not to change `scf_fp32_stall_patience` until the fixed-point question is settled.
3. **Pseudo-diagonalisation on ROCm.** CUDA has `-scf_pseudo_diag` (113 -> 89 s on polymer_2x, 4 GPUs, MOR41 95/95
   identical, GMTKN55 2459/2460); ROCm has none. A pseudo step costs three GEMMs; E: FP32 about 40 ms against 640 ms,
   FP64 about 470 ms against 1550 ms at n = 3222. Largest single lever, largest effort.
4. **FP32 eigenvectors plus FP64 correction instead of the FP64 solves** (CUDA plan item 0.2, `GFN2_GPU_COST_PLAN.md`).
5. **D4 ATM 3-body (1.2 s) and integrals (1.2 s)**: FP64-bound single kernels; cutoff and tiling untested here.
6. **Small systems**: plugin start 0.2 s, per-call overhead, host-side potential and gradient parts
   (`ROCM_EVALUATION_2026-10-13.md`). Cheapest fix: `-gpu auto` uses the CPU for GFN1/GFN2 below about 1000 atoms (the
   GFN-FF rule applied to xTB); the threshold needs a measurement on more system types than water and peo.

## Not measured

MD or `-opt` per-step cost of GFN2 on ROCm (only whole optimisations, n = 3); gradient phases on the device; gfn1
eigensolve split; memory; gfx1150; systems with transition metals; CPU large-system modes (`-large_system_mode`) as the
comparison at 1410 atoms (the CPU time above is the default mode).
