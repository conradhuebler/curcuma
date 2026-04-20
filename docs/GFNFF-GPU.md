# GPU GFN-FF Implementation — Status & Analysis

Updated 2026-04-12. All optimizations from the roadmap (G1a–G2a) are implemented.

---

## Completed Optimizations

| ID | Optimization | Status | Result |
|----|-------------|--------|--------|
| G1a | Dispersion pair sorting by idx_i | ✅ Done | GPU 1.50× (20.0 vs 30.0 ms/step) |
| G1b | `__ldg()` for SoA pair parameters | ✅ Done | Combined with G1a measurement |
| G1c | Coulomb + Repulsion pair sorting | ✅ Done | No additional effect at 1410 atoms |
| G2a | `__launch_bounds__(256,2)` for bonded kernels | ✅ Done | GPU 1.58× cumulative |
| G2b | Nsight Compute profiling | ⏭️ Skipped | Not installed; used ptxas -v instead |

### G2a Register Pressure Results (ptxas -v, sm_75)

| Kernel | Before (512,2) | After (256,2) | Stack Before | Stack After | Spills After |
|--------|----------------|---------------|-------------|-------------|-------------|
| k_angles | 64 regs | 82 regs | 88 B | 40 B | 0 |
| k_storsions | 64 regs | 106 regs | 256 B | 40 B | 0 |
| k_inversions | 64 regs | 116 regs | 264 B | 40 B | 0 |
| k_dihedrals | 64 regs | 128 regs | 432 B | 88 B | 36/36 B |
| k_hbonds | 64 regs | 128 regs | 632 B | 312 B | 548/780 B |
| k_dispersion | 43 regs | 43 regs | 0 | 0 | 0 |
| k_repulsion | 40 regs | 40 regs | 0 | 0 | 0 |
| k_coulomb | — | — | — | — | — |

k_angles, k_storsions, k_inversions have **zero spills** after G2a. k_dihedrals reduced from 432→88 bytes stack. k_hbonds still spills significantly — would require kernel splitting per case type.

---

## Current Architecture

### Memory Layout

- **Coordinates**: SoA (separate cx[N], cy[N], cz[N]) — fully coalesced warp reads
- **Pair indices**: Separate idx_i[n], idx_j[n] — gather access (coalescing depends on sort order, now sorted by idx_i for dispersion/coulomb/repulsion)
- **Pair parameters**: SoA layout (C6[tid], r4r2ij[tid], etc.) — coalesced sequential reads
- **Constant memory**: d_rcov_d3[87], d_refcn_const[826], d_refn_const[118], d_lattice[9], d_lattice_inv[9], d_has_pbc, d_cov_radii[100] — ~8 KB total

### Launch Configuration

- **Lightweight kernels** (dispersion, repulsion, Coulomb, CN, etc.): `__launch_bounds__(512, 2)` — 64 regs/thread budget
- **Heavy kernels** (angles, dihedrals, inversions, storsions, hbonds): `__launch_bounds__(256, 2)` — 128 regs/thread budget, eliminates most register spilling
- **Adaptive block sizing**: `getLaunchConfig(n, maxBlockSize)` — 32/128/256/512 threads depending on problem size

### Stream Topology

Three concurrent streams during Phase 1 (charge-independent):
- **Stream A** (pairwise): k_dispersion, k_repulsion (bonded + non-bonded), k_hb_alpha_chainrule
- **Stream B** (bonded): k_bonds, k_angles, k_dihedrals, k_inversions, k_storsions
- **Stream C** (three-body): k_batm, k_atm, k_xbonds, k_hbonds

Phase 2 (charge-dependent, after EEQ):
- **Stream A**: k_dispersion (gradient, deferred), k_coulomb
- **Main stream**: k_coulomb_postprocess, k_cn_chainrule

### Energy Reduction

`blockReduceAddEnergy()` uses warp shuffle (`__shfl_down_sync`) + shared memory aggregation. One `atomicAdd` per block instead of per thread. 256 bytes shared memory per block.

---

## Remaining Bottlenecks

### 1. k_hbonds Register Spilling (548/780 bytes)

The most complex kernel with 4 case types, neighbor loops, and nested `dihedral_angle_grad` calls. At 128 regs/thread, still spills 548 bytes. Options:
- Split into per-case-type kernels (complex, requires SoA restructuring)
- Use shared memory for intermediate arrays (limited benefit, 256 bytes already used for energy reduction)
- Accept spilling (HB is a rare term, marginal overall impact)

### 2. Coordinate Gather (Uncoalesced)

Pair kernels access `cx[i], cx[j]` where i,j are atom indices. Even with idx_i sorting, j indices are still random. This is inherent to the algorithm — texture memory or shared memory tiling could help but adds complexity.

### 3. Atomic Gradient Accumulation

6 `atomicAdd` per pair per kernel (grad[3*atom + dim]). High-coordination atoms (e.g., central carbon) receive many concurrent writes. Block-level gradient reduction would reduce contention but adds shared memory complexity.

### 4. No Further GPU Optimizations Recommended

- **Texture memory**: Legacy API, `__ldg()` already provides equivalent functionality
- **Cooperative groups**: No benefit over existing warp shuffle reduction
- **CUTLASS/Tensor Cores**: GFN-FF is not a matrix problem
- **Further kernel splitting**: Marginal benefit for HB (rare term)

---

## Performance Summary (2026-04-12)

**System:** AMD Ryzen 9 9950X3D, NVIDIA GeForce RTX 5080, 4 Threads
**Testfall:** Polymer, 1410 Atome, 1000 MD-Schritte

| Methode | ms/Schritt | Faktor vs. Baseline |
|---------|-----------|---------------------|
| gfnff-gpu (nach G1a+b+G2a) | 19.0 ms | 1.58× |
| gfnff-gpu (Baseline) | 30.0 ms | 1.00× |
| gfnff CPU (4 Threads) | 148.8 ms | 0.20× |

GPU is now **7.8× faster** than CPU (was 5.3× at baseline).

---

## Kernel Inventory (22 kernels)

| Kernel | Launch Bounds | Registers | Stack | Purpose |
|--------|--------------|-----------|-------|---------|
| k_dispersion | 512,2 | 43 | 0 | D4 dispersion |
| k_repulsion | 512,2 | 40 | 0 | Bonded + non-bonded repulsion |
| k_repulsion_mixed | 512,2 | 30 | 0 | FP32 intermediates variant |
| k_coulomb | 512,2 | — | — | Coulomb TERM 1 |
| k_coulomb_self | 512,2 | 28 | 0 | Coulomb TERM 2+3 |
| k_coulomb_postprocess | 512,2 | 35 | 0 | Fused self-energy + qtmp |
| k_subtract_qtmp | 512,2 | 25 | 0 | Coulomb TERM 1b |
| k_cn_chainrule | 512,2 | 38 | 0 | CN gradient chain rule |
| k_hb_alpha_chainrule | 512,2 | 37 | 0 | HB alpha gradient |
| k_bonds | 512,2 | 44+ | 256+ | Bond stretching |
| k_angles | 256,2 | 82 | 40 | Angle bending |
| k_dihedrals | 256,2 | 128 | 88 | Dihedral torsion |
| k_inversions | 256,2 | 116 | 40 | Out-of-plane inversion |
| k_storsions | 256,2 | 106 | 40 | Triple bond torsion |
| k_batm | 512,2 | 64 | 40 | Bonded ATM |
| k_batm_mixed | 512,2 | 46 | 0 | FP32 ATM variant |
| k_atm | 512,2 | 64 | 40 | 3-body dispersion |
| k_xbonds | 512,2 | 64 | 72 | Halogen bonds |
| k_xbonds_mixed | 512,2 | 54 | 0 | FP32 XB variant |
| k_hbonds | 256,2 | 128 | 312 | Hydrogen bonds |
| k_gaussian_weights | 512,2 | 48 | 0 | D4 Gaussian weight computation |
| k_dc6dcn_per_pair | 512,2 | 50 | 0 | D4 dc6/dcn per dispersion pair |
| k_cn_compute | 512,2 | 36 | 0 | Coordination number |
| k_check_displacement | 512,2 | 6 | 0 | Topology cache displacement check |

---

*Original analysis: March-April 2026. Updated: April 2026 after G1a+b+G1c+G2a implementation.*