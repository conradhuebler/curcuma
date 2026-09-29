# GFN2 on several GPUs: where the time goes and what to change

Status (Sep 29, 2026): stage 0 measured, stage 1a done, 1b tried and rejected, 1c deferred;
stages 2-3 open (design below, awaiting the operator's go).
All numbers: 4x RTX A4500 (20 GB, PCIe, no NVLink), polymer_2x (7320 atoms, nao 15444, 19252
electrons, nocc 9626), commit `bfbdc569`, `CURCUMA_GPU_PROFILE=1`.

## Baseline: one single point, 116 s on 4 GPUs, 12 SCF iterations

| item | time | share | runs on |
|---|---:|---:|---|
| FP32 eigensolve (11 x, cuSOLVERMp `syevd`) | 51.9 s | 45 % | 4 GPUs |
| final FP64 eigensolve (1 x: `sygst` 4.2 + `syevd` 19.5 + back-transform 2.7) | 26.9 s | 23 % | 4 GPUs |
| density P + populations (12 x, `k_density_sp`, 0.74 s per call) | 10.2 s | 9 % | device 0 |
| FP32 reduction L^-1 F L^-T (11 x, two `strsm`) | 6.3 s | 5 % | device 0 |
| setup (see 0.3) | ~15 s | 13 % | mixed |
| potential, SCC energy, post-SCF D4 | ~6 s | 5 % | device 0 / host |

## Stage 0: measurements (done)

- **0.1 Process grid.** `test_cases/cuda/bench_syevd_mp.cpp` now takes the grid rows as 5th
  argument. n = 15444, 4 GPUs: FP32 1x4 nb128 5.28 s, 2x2 nb64/128/256 5.59/5.32/5.23 s; FP64 1x4
  20.54 s, 2x2 21.62/21.00/20.98 s. **No gain** - the 1 x ndev layout stays. The block size was
  already tuned (GPU_TUNING.md).
- **0.2 Can FP32 eigenvectors plus FP64 correction replace the final FP64 solve?**
  `CURCUMA_DUMP_FS=<file>` (xtb_scf.cpp, zero cost when unset) writes the last iteration's F and S;
  a standalone Eigen probe compared, on polymer (nao 3222, CPU run, `-scf_threshold 1e-9`), against
  an FP64 reduction + FP64 eigensolve:

  | variant | max\|dP\| | Tr(A dP) (Eh) | HOMO error (Eh) |
  |---|---:|---:|---:|
  | FP32 vectors as returned | 5.3e-5 | -9.4e-3 | 4.5e-7 |
  | + FP64 re-orthonormalisation | 1.9e-6 | 1.0e-8 | 9.5e-12 |
  | + 1 FP64 rotation step | **7.9e-12** | 4.2e-12 | 1.5e-11 |
  | + 2 FP64 rotation steps | 5.7e-14 | 4.2e-12 | 1.5e-11 |

  The rotation step is the first-order occupied-virtual Jacobi update
  theta_ai = G_ai / (G_ii - G_aa), G = C^T A C, followed by re-orthonormalisation
  (Stewart, Csaszar, Pulay, J. Comput. Chem. 3 (1982) 227). The occupied-virtual block shrinks
  7e-7 -> 6e-12 -> 2e-16 (quadratic). **Stage 3 is viable**; one probe, one system - the gate in
  the real code is the SCF energy and gradient against the current FP64 path.
- **0.3 Setup of polymer_2x (~13 s, nsys).** Kernels 3.4 s (overlap/H0 1.3, distributed Cholesky
  1.7, multipole/gamma 0.4); first NCCL/cuSOLVERMp use ~1.5 s (pinned allocations,
  `cuMemSetAccess`; on polymer the distributed Cholesky takes 1.72 s of which 0.44 s factorise,
  vs 41 ms on one GPU); host download + copy of S/H0/L/gamma 3.7 s; host multipole interaction
  matrices 1.5 s; pre-SCF EEQ guess on device 0 ~1.5 s.

## Stage 1 results

- **1a done.** S/H0 are downloaded column-major straight into `m_S`/`m_H0` and transposed in
  place (exact element swaps), L and gamma straight into `m_X`/`m_gamma`, and the screened-storage
  scatter runs on up to 8 host threads. polymer_2x: download + host copy **3.67 -> 1.07 s**
  (downloads 1.3 s, host copies 2.36 s before). Checked element-wise against the old path on
  caffeine (dense storage) and polymer (screened storage): all four matrices identical. A full
  deferral (no download at all on the resident path) would save the remaining ~1 s and 5.7 GB of
  host memory; not done, it touches every host consumer of these matrices.
- **1b rejected (measured).** The density as an FP32 `syrk` per GPU plus a gather of the stored
  pairs cut the density step from 0.85 to 0.29 s, but every FP32 iteration's SCC energy came out
  **2.7 mEh too high** (-11857.3256 vs -11857.3287 Eh) - the FP32 resolution of an 11857 Eh
  electronic energy (2e-7 relative). The SCF then cannot meet its energy criterion in FP32 and
  switched to FP64 three times: 16 instead of 12 iterations, 178 instead of 113 s. The FP32
  eigenvectors themselves do not cause this because their density is accumulated in FP64. Kept
  out of the code; any faster density must accumulate in FP64 (e.g. a tiled FP64 kernel that
  reuses C rows - the current kernel reads C once per stored pair and runs at ~93 GFLOP/s, i.e.
  memory-bound).
- **1c deferred.** The FP32 reduction already runs at ~13 TFLOP/s on device 0 (two `strsm`,
  7.4 TFLOP in 0.57 s, near the A4500's FP32 peak), so four GPUs could save at most ~0.4 s per
  iteration. Stage 2 below works in the AO basis and needs no reduction in its iterations, which
  leaves the reduction only in the few full-diagonalisation steps.

## Stages 2-3: design after stage 0/1

Pseudo-diagonalisation in the **AO basis** (C^T S C = I), orbitals split over the GPUs by
virtual columns, FP32 like today's mixed-precision iterations, density still accumulated in FP64:

1. Broadcast F (FP32, 0.95 GB) to all GPUs; S is geometry-constant and stored once per GPU.
2. GPU k: T_k = F C_virt[:,J_k], G_ov[:,J_k] = C_occ^T T_k (C_occ replicated, 0.6 GB FP32).
3. theta = G_ai / (G_ii - G_aa) (eigenvalue estimates: G_ii = diag C^T F C, also split).
4. C_occ += C_virt theta (partial per GPU, allreduce), C_virt -= C_occ theta^T (local).
5. Re-orthonormalise C_occ against S (C_occ^T S C_occ, Cholesky, trsm), split by columns.

Estimate ~0.8 s per iteration on 4 A4500 against 4.7 s (eigensolve) + 0.57 s (reduction) today.
Full diagonalisation stays on the first iteration, whenever the SCF stalls or the HOMO-LUMO gap
is small against kT, and every k-th iteration as a safety net. Stage 3 is the same step in FP64
once, after the last FP32 full solve, replacing the FP64 `sygst + syevd + trsm` (26.9 s; est.
8-9 s, and no FP64 reduction needed).

## Stages 1-3 (planned)

1. **Low risk, results unchanged.** (a) Do not download S/H0/L/gamma on the resident path; build
   the host copies on demand (as for the multipole integrals, `xtb_native.cpp:599`): -3.7 s per
   geometry, -5.7 GB host memory. (b) Density P over all GPUs, FP32 in the FP32 iterations:
   10.2 s -> est. 2-3 s. (c) Reduction on all GPUs: L^-1 once per geometry, each GPU forms its own
   block-cyclic columns with two GEMMs, which also removes the scatter: 6.3 s -> est. 1.5 s.
2. **Pseudo-diagonalisation** in the middle SCF iterations (the rotation step of 0.2 in FP32,
   ~10 TFLOP per iteration, GEMM only): 11 x 4.7 s -> est. 3 full + 8 x ~0.3 s. Full solve on the
   first iteration, before convergence, when the HOMO-LUMO gap is small against kT, and on stalls.
   Intermediate iterations are no longer bit-identical; the gate is the converged energy.
3. **Final FP64 solve replaced** by the last FP32 solve + FP64 reduction + one FP64 rotation step
   (0.2): 26.9 s -> est. 8-10 s.

Estimates are back-of-envelope; each stage gets its own before/after measurement. FP64 is ~1/32
of FP32 on the A4500 and near 1:1 on an H200, so stage 3's balance differs there (not measurable
here).
