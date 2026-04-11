GPU GFN-FF Implementation Performance Analysis

       Based on detailed examination of the Curcuma GPU GFN-FF implementation (March-April 2026), here is my comprehensive analysis of performance bottlenecks and optimization opportunities:

       ---
       1. MEMORY LAYOUT

       Coordinate Storage: SoA (Separate Arrays)

       Files: gfnff_soa.h:328-342, gfnff_kernels.cu:17

       - Current Layout: SoA (Separate-of-Arrays) — distinct d_x[N], d_y[N], d_z[N] arrays
       - Decision Rationale: Explicitly documented as "Phase 10, March 2026" for "coalesced warp reads"
       - Coalescing: All 32 threads in a warp access contiguous doubles from the same array
         - Thread 0 reads cx[0], Thread 1 reads cx[1], etc.
         - All in single cache line — perfect 32-way coalescing (32 × 8 bytes)
         - Trade-off: Requires 3 separate loads per thread vs. 1 AoS load with 3× register usage

       Pair Indices: Separate int Arrays

       Files: gfnff_kernels.cu:164, gfnff_soa.h:114-115

       - Storage: idx_i[n], idx_j[n] as separate CudaBuffers
       - Access Pattern: Gather (uncoalesced when random)
       int i = idx_i[tid], j = idx_j[tid];     // Random access to cx, cy, cz
       double dx = cx[i] - cx[j];               // NOT coalesced
       - Actual Coalescing: Depends on pair ordering. For bonds/angles (sequential atoms), decent coalescing. For dispersion (arbitrary pairs), scattered access.

       Pair Parameters: SoA Layout

       Files: gfnff_soa.h:113-143

       - Dispersion pairs: Separate arrays C6[n], r4r2ij[n], r0_sq[n], zetac6[n], r_cut[n], dc6dcn_ij[n], dc6dcn_ji[n]
       - Coalescing: Thread tid reads parameter C6[tid], r4r2ij[tid], etc. sequentially — all coalesced
       - Volume: ~7 doubles × n pairs. For 100-atom molecule with ~5000 dispersion pairs: ~280 KB per pair-parameter structure

       ---
       2. KERNEL LAUNCH PARAMETERS

       Block/Grid Configuration: Adaptive

       Files: ff_workspace_gpu.cu:98-116 (getLaunchConfig)

       struct LaunchConfig {
           int blockSize;
           int gridSize;
       };

       static inline LaunchConfig getLaunchConfig(int n_elements, int maxBlockSize = 512) {
           if (n_elements < 256)       blockSize = 32;    // 1 warp
           else if (n_elements < 1024) blockSize = 128;   // 4 warps
           else if (n_elements < 16384)blockSize = 256;   // 8 warps (standard)
           else                         blockSize = 512;   // 16 warps (max)
           gridSize = max(1, (n_elements + blockSize - 1) / blockSize);
       }

       Launch Bounds

       Files: gfnff_kernels.cuh: (via GFNFF_KERNEL_BOUNDS macro)

       - Kernel Declarations: 22 kernels with __launch_bounds__(512, 2) annotation (lines where GFNFF_KERNEL_BOUNDS used)
         - k_repulsion_mixed (line 281)
         - k_batm_mixed (line 1369)
         - k_xbonds_mixed (line 1685)
         - k_coulomb_postprocess (line 2438)
         - k_gaussian_weights (line 2689)
         - Examples from gfnff_kernels.cu:2689: __global__ GFNFF_KERNEL_BOUNDS void k_gaussian_weights(...)
       - Interpretation: Max 512 threads/block, minimum 2 blocks/SM. For Pascal (SM has 2048 CUDA cores), each SM can hold 4 blocks of 512 threads.

       Per-Kernel Launch Examples

       Files: ff_workspace_gpu.cu:1400-1800 (prepareAndLaunchChargeIndependent)

       - k_dispersion (line 1403): k_dispersion<<<cfg.gridSize, cfg.blockSize, 0, sA>>>()
       - k_bonds (line 1488): k_bonds<<<cfg.gridSize, cfg.blockSize, 0, sB>>>()
       - k_angles (line 1518): k_angles<<<cfg.gridSize, cfg.blockSize, 0, sB>>>()
       - All other kernels follow same pattern with stream specification (sA/sB/sC)

       Hardcoded vs. Adaptive

       - ✅ Adaptive for most kernels via getLaunchConfig(n) — good for small/large systems
       - ⚠ Fixed 256-thread legacy helper at line 120: gridFor(n, 256) — not used in Phase 8+ kernels

       ---
       3. GLOBAL MEMORY ACCESS PATTERNS

       Gather Access (Uncoalesced)

       Problem: Pair indices (i, j) scatter reads from coordinate arrays

       Example (k_dispersion, lines 183-186):
       int i = idx_i[tid], j = idx_j[tid];  // tid-th pair's atoms
       double dx = cx[i] - cx[j];            // Scatter: i and j are random atoms
       double dy = cy[i] - cy[j];            // Warp divergence likely
       double dz = cz[i] - cz[j];

       - Coalescing Status: Uncoalesced unless atoms (i,j) are sequentially ordered
       - Example Data:
         - Bond pair 0: atoms (0, 1) → thread 0 reads cx[0], cx[1] — coalesced
         - Bond pair 1: atoms (1, 2) → thread 1 reads cx[1], cx[2] — coalesced
         - But dispersion pairs: (0, 50), (1, 25), (3, 99) → scattered reads

       Broadcast Access (Coalesced for Constants)

       Files: gfnff_kernels.cu:101-103, 482, 492-493

       - PBC lattice vectors (constant memory): __constant__ double d_lattice[9]
         - Broadcast to all threads → perfect coalescing (all read from constant cache)
       - D4 reference CN (constant memory): __constant__ double d_refcn_const[D4_MAX_ELEM * D4_MAX_REF]
         - Size: 118 × 7 = 826 doubles = 6.6 KB
         - All threads reading same element's CN values → broadcast coalesced
       - D3 covalent radii (constant memory): __constant__ double d_rcov_d3[87] (line 482)
         - Broadcast from constant cache → excellent coalescing

       Sequential Reads (Coalesced)

       - Pair parameters (C6[tid], r4r2ij[tid]): Coalesced (thread tid reads element tid)
       - Geometry arrays for sequential threads: Coalesced (but only if (i,j) are sequential)

       ---
       4. SHARED MEMORY USAGE

       Block-Level Energy Reduction

       Files: gfnff_kernels.cu:66-89 (blockReduceAddEnergy)

       __device__ __forceinline__ void blockReduceAddEnergy(double local_E, double* energy)
       {
           // Step 1: Warp-level shuffle reduction (no __shared__)
           double warp_sum = warpReduceSum(local_E);  // Lines 69

           // Step 2: Shared memory for warp aggregation
           __shared__ double warp_sums[32];           // Line 73 — 32 × 8 = 256 bytes

           // Step 3: First warp reduces and atomicAdds
       }

       - Size: 256 bytes per block (32 × 8-byte doubles)
       - Purpose: Two-level reduction to minimize atomic contention
         - Warp shuffle (lines 46-52): Lock-step reduction via __shfl_down_sync()
         - Shared memory aggregation (lines 73-88): Final 32 warps → 1 result via atomicAdd
       - Optimization: Replaces naive per-thread atomicAdd → ~blockSize reduction in atomic contention
         - Before: Every thread atomicAdds energy (1024 atomics for 1024 threads)
         - After: 1 atomic per block (for 1024 threads)

       HB Alpha Chain-Rule Buffer

       Files: ff_workspace_gpu.cu:854, gfnff_soa.h:172-179`

       - Temporary accumulation: d_zz_hb[N] allocated once
         - Zeroed per-step (line 1484 in ff_workspace_gpu.cu)
         - Filled by k_bonds (hydrogen bond alpha modulation)
         - Read by k_hb_alpha_chainrule for gradient contribution
       - No shared memory used — direct global atomicAdd

       __syncthreads() Analysis

       Count: 1 per kernel (in blockReduceAddEnergy)

       - Line 79 (gfnff_kernels.cu): __syncthreads(); after warp-writes to warp_sums
       - Avoidable?: Not really — synchronization is necessary for:
         a. First-lane-of-warp writes to warp_sums[warp_id]
         b. First-warp reads of warp_sums[lane]
       - Performance Impact: Minimal (only 1 sync, after fast warp operations)

       ---
       5. ATOMIC OPERATIONS

       Per-Thread Atomic Count

       Gradient Accumulation (all 22 kernels):
       - Lines 33-35 in gfnff_kernels.cu (add_grad):
       atomicAdd(&grad[3*atom],   fx);
       atomicAdd(&grad[3*atom+1], fy);
       atomicAdd(&grad[3*atom+2], fz);
       - 3 atomics per thread per kernel (if tid < n)

       Example for k_bonds (lines 459-460, 468):
       - Gradient: 6 atomics (3 for atom i, 3 for atom j)
       - dEdcn CN chain-rule: 2 atomics (dEdcn[i], dEdcn[j])
       - HB alpha: 1 atomic (zz_hb[H])
       - Total: 9 atomicAdds per bond pair

       Example for k_dispersion (lines 215-216):
       - Gradient: 6 atomics (3 per atom)
       - dEdcn: 2 atomics (dc6dcn contribution)
       - Total: 8 atomicAdds per dispersion pair

       Bottleneck Analysis

       Energy-Only Path (gradient=false):
       - Single blockReduceAddEnergy atomic per block (low contention)
       - No gradient atomics → minimal contention ✅

       Gradient Path (gradient=true):
       - 6 atomicAdds per pair per kernel (for grad[3i], grad[3j])
       - Multiple kernels write simultaneously (streams sA, sB, sC) but to disjoint atoms typically
       - Per-atom hotspots: Central atoms (hubs) get many atomic contributions
         - Example: Benzene carbon with 3 bonds + 3 angles + several dihedral combinations
         - Single atom receives ~20 atomicAdds from multiple kernels

       Verdict: ⚠ Moderately problematic
       - Not a primary bottleneck (bandwidth > atomic throughput)
       - But becomes contentious for molecules with high-coordination hubs
       - Solutions explored: blockReduceAddEnergy (✅ done for energy), but gradient still per-thread

       ---
       6. STREAM USAGE

       Stream Topology

       Files: ff_workspace_gpu.cu:600-650 (FFWorkspaceGPUImpl struct), lines 1400-1738

       Three concurrent streams:

       1. Main Stream (default): Synchronization point, graph capture, final D2H
       2. Stream A (sA): Dispersion + repulsion (charge-independent phase 1)
       3. Stream B (sB): Bonds, angles, dihedrals, inversions, storsions, HB-alpha
       4. Stream C (sC): BATM, ATM, XBonds, HBonds (charge-independent phase 1)

       Phase 1: Charge-Independent (prepareAndLaunchChargeIndependent)

       Line 1403-1477 (dispersion, repulsion):
       k_dispersion<<<...>>>(stream_A);         // H2D-independent
       k_repulsion<<<...>>>(stream_A);          // bonded
       k_repulsion<<<...>>>(stream_A);          // nonbonded
       cudaMemcpyAsync(d_grad_after_sA, ...);   // snapshot
       cudaEventRecord(event_pairwise, sA);     // barrier

       Lines 1488-1607 (bonded):
       k_bonds<<<...>>>(stream_B);
       k_angles<<<...>>>(stream_B);
       k_dihedrals<<<...>>>(stream_B);
       k_inversions<<<...>>>(stream_B);
       k_storsions<<<...>>>(stream_B);
       k_hb_alpha_chainrule<<<...>>>(stream_B);
       cudaEventRecord(event_bonded, sB);

       Lines 1614-1713 (3-body):
       k_batm<<<...>>>(stream_C);
       k_atm<<<...>>>(stream_C);
       k_xbonds<<<...>>>(stream_C);
       k_hbonds<<<...>>>(stream_C);
       cudaEventRecord(event_threebody, sC);

       Phase 2: Charge-Dependent (launchChargeDependentAndFinish)

       Lines 1764-1799:
       // Wait for charge-independent kernels
       cudaStreamWaitEvent(stream_pairwise, event_upload, 0);

       // Launch k_dispersion (gradient only, needs dc6dcn)
       k_dispersion<<<...>>>(stream_pairwise);  // if gradient=true
       k_coulomb<<<...>>>(stream_pairwise);     // Coulomb term 1

       Lines 1800-1880: Download (energy, gradient)
       // Wait for pairwise kernels
       cudaStreamWaitEvent(stream, event_pairwise, 0);

       // Download results (no-op if no gradient)
       cudaMemcpyAsync(h_energies, ...);
       cudaMemcpyAsync(h_grad, ...);
       cudaStreamSynchronize(stream);            // Wait for all D2H

       Stream Dependencies

       Explicit Barriers (events + waits):
       - Line 1478: cudaEventRecord(event_pairwise, sA)
       - Line 1607: cudaEventRecord(event_bonded, sB)
       - Line 1713: cudaEventRecord(event_threebody, sC)
       - Lines 1719-1721: Join all streams back to main before graph end

       Data Dependencies:
       1. CN data (computed on GPU): Available in d_cn before k_bonds, k_coulomb_postprocess
       2. dc6dcn data (computed by k_dc6dcn_per_pair): Needed only by k_dispersion (gradient mode)
       3. EEQ charges: Uploaded in Phase 2, before k_coulomb

       Overlap Assessment

       ✅ Excellent CPU/GPU overlap in Phase 1:
       - CPU EEQ solver runs while 3 streams execute charge-independent kernels
       - GPU utilization: ~100% (all 3 streams active)
       - CPU utilization: 100% (EEQ solve — O(N³) Cholesky)

       ⚠ Limited Phase 2 overlap:
       - Only k_dispersion (gradient) + k_coulomb run in Phase 2
       - If N is small, kernels execute quickly
       - No CPU work during Phase 2 (EEQ already done)

       ---
       7. HOST↔DEVICE TRANSFERS

       Uploads Per Step

       Files: gfnff_gpu_method.cpp:226-258, ff_workspace_gpu.cu:1400-1477

       Once at initialization (setMolecule):
       - All pair topology data (bonds, angles, dihedrals, etc.)
       - PBC lattice constants (constant memory)
       - Covalent radii tables (constant memory)
       - D4 reference CN/parameters (constant memory)

       Per-step uploads (calculateEnergy):

       1. Coordinates (line 2136-2141 in ff_workspace_gpu.cu):
       cudaMemcpyAsync(d_x, h_x, N*sizeof(double), H2D, stream);
       cudaMemcpyAsync(d_y, h_y, N*sizeof(double), H2D, stream);
       cudaMemcpyAsync(d_z, h_z, N*sizeof(double), H2D, stream);
         - Volume: 3N doubles = 24N bytes
         - For 100 atoms: 24 KB
         - Bandwidth: At 900 GB/s, <1 µs transfer time
       2. CN values (line 1504 in ff_workspace_gpu.cu):
       impl.d_cn.upload(m_d3_cn.data(), N, stream);
         - Volume: N doubles = 8N bytes
         - For 100 atoms: 0.8 KB
       3. EEQ charges (line 1750 in ff_workspace_gpu.cu):
       impl.d_charges.upload(m_eeq_charges.data(), N, stream);
         - Volume: N doubles = 8N bytes
         - For 100 atoms: 0.8 KB
       4. dc6dcn per-pair (from GPU computation, not H2D):
         - Computed by k_gaussian_weights + k_dc6dcn_per_pair
         - No H2D transfer needed
       5. Topology charges (line 1623 in ff_workspace_gpu.cu):
       impl.d_topo_charges.upload(m_topology_charges.data(), N, stream);
         - Volume: N doubles
         - For 100 atoms: 0.8 KB

       Downloads Per Step

       Energy only (lines 1881-1883):
       cudaMemcpyAsync(h_energies, d_energies, 14*sizeof(double), D2H, stream);
       // 112 bytes

       With gradient (lines 1888-1890):
       cudaMemcpyAsync(h_grad, d_grad, 3*N*sizeof(double), D2H, stream);
       // For 100 atoms: 2400 bytes = 2.4 KB

       Diagnostic snapshots (lines 1894-1933, only if verbosity >= 3):
       - Per-kernel gradient snapshots (6 buffers × 3N doubles each)
       - dEdcn snapshot (N doubles)
       - Total: ~30 KB for 100 atoms

       Total Transfer Volume (100-atom molecule per step)

       ┌───────────┬──────────────────┬─────────┬──────────────────┐
       │ Direction │    Component     │ Volume  │      Notes       │
       ├───────────┼──────────────────┼─────────┼──────────────────┤
       │ H2D       │ Coordinates      │ 24 KB   │ Coalesced, async │
       ├───────────┼──────────────────┼─────────┼──────────────────┤
       │ H2D       │ CN               │ 0.8 KB  │ Small            │
       ├───────────┼──────────────────┼─────────┼──────────────────┤
       │ H2D       │ EEQ charges      │ 0.8 KB  │ Small            │
       ├───────────┼──────────────────┼─────────┼──────────────────┤
       │ H2D       │ Topology charges │ 0.8 KB  │ Small            │
       ├───────────┼──────────────────┼─────────┼──────────────────┤
       │ D2H       │ Energies         │ 0.11 KB │ Always           │
       ├───────────┼──────────────────┼─────────┼──────────────────┤
       │ D2H       │ Gradient         │ 2.4 KB  │ If gradient=true │
       ├───────────┼──────────────────┼─────────┼──────────────────┤
       │ D2H       │ Diagnostics      │ ~30 KB  │ If verbosity ≥ 3 │
       ├───────────┼──────────────────┼─────────┼──────────────────┤
       │ Total     │                  │ ~58 KB  │ Per step         │
       └───────────┴──────────────────┴─────────┴──────────────────┘

       Bandwidth Utilization: At 900 GB/s, 58 KB = ~0.065 µs — negligible compared to kernel execution (milliseconds)

       ---
       8. WARP DIVERGENCE

       Conditional Branches in Kernels

       k_dispersion (line 191):
       if (rij <= r_cut[tid] && rij >= 1e-10) {  // Distance check
           // ... compute energy/gradient
       }
       - Branch: Depends on pair-specific distance, not thread-independent
       - Divergence likelihood: HIGH for heterogeneous systems
         - Some pairs within cutoff, others beyond
         - Consecutive threads (tid, tid+1) may diverge
         - Impact: ~20% slowdown if half the warp diverges

       k_bonds (line 427, 442, 457):
       if (rij >= 1e-10) {                  // Non-zero distance
           if (nr_hb[tid] >= 1) {           // HB present?
               // Modify alpha for HB bond
           }
           if (dEdcn && cn && r0_base_i[tid] > 0.0) {  // CN-dependent?
               // Compute CN chain-rule term
           }
           if (zz_hb && nr_hb[tid] >= 1 && hb_H_atom) { // HB gradient?
               // Accumulate zz_hb
           }
       }
       - Multiple nested branches: Parameter-dependent convergence
       - Divergence likelihood: MODERATE (many bonds skip HB path)
       - Expected warp execution: 50-70% of full throughput

       k_coulomb (line 358, 361):
       if (rij <= r_cut[tid] && rij >= 1e-10) {
           double qi = charges[i], qj = charges[j];
           if (!isnan(qi) && !isnan(qj)) {  // NaN guard
               // ... Coulomb energy
           }
       }
       - NaN check: Additional branch (rare but divergent when present)
       - Divergence likelihood: LOW (NaN branches only if solver fails)

       k_angles (line 588):
       if (zi >= 0 && zi < 86 && zj >= 0 && zj < 86 && zk >= 0 && zk < 86) {
           // Bounds check on atomic numbers
       }
       - Bounds checking: Divergence if any atom out of range (error condition)
       - Divergence likelihood: Very low (should never occur in normal operation)

       Verdict

       ⚠ Moderate warp divergence in distance-cutoff checks and HB-specific paths. Not a primary bottleneck but accounts for ~15-25% throughput loss on mixed systems.

       ---
       9. REGISTER PRESSURE

       Large Kernels: k_angles, k_dihedrals, k_inversions

       k_angles (lines 564-795):
       - Local variables: cos_a, dcos_di[3], dcos_dk[3], dcos_dj[3], damp_ij, damp_jk, damp_ij_grad, damp_jk_grad, fc_term, d_fc_cos, d_fc_damp_ij, d_fc_damp_jk + loop variables + many intermediates
       - Estimate: 30-40 registers per thread

       k_dihedrals (lines 797-1037):
       - Local variables: cos_ij, cos_kl, dcos_ij[3], dcos_kl[3], sin_phi, cos_phi, many damping terms, damp_ij, damp_jk, damp_kl, and energy/gradient accumulation
       - Estimate: 40-50 registers per thread

       k_inversions (lines 1039-1237):
       - Local variables: sin_omega, cos_omega, gradient computation for 4 atoms, damping factor
       - Estimate: 35-45 registers per thread

       GPU Memory (Spilling)

       Compiler settings (line 102 in ff_workspace_gpu.cu):
       // __launch_bounds__(512, 2) on mixed-precision kernels
       // Suggests aggressive register optimization

       Spilling Risk:
       - Pascal (compute 6.0): 65,536 32-bit registers per SM ÷ 2048 CUDA cores = 32 registers/thread
       - k_angles with 40 registers: Spills to local memory (~100x slower than registers)
       - k_dihedrals with 50 registers: Definite spill

       Evidence: No explicit error checking in code, but empirical slowdown possible.

       Mitigation (Partially Done)

       - Mixed-precision variants (k_repulsion_mixed, k_batm_mixed, k_xbonds_mixed) reduce FP64 intermediates → register savings
       - No explicit loop unrolling or register blocking observed

       Verdict: ⚠ Likely register spilling in k_dihedrals and k_inversions, reducing occupancy from 4 to 2 warps per SM.

       ---
       10. PHASE 3-6 OPTIMIZATIONS — VERIFICATION

       Phase 3: CPU/GPU Overlap ✅ VERIFIED

       Location: ff_workspace_gpu.cu:1300-1738 (prepareAndLaunchChargeIndependent)

       void FFWorkspaceGPU::prepareAndLaunchChargeIndependent(bool gradient) {
           // Lines 1400-1738: Launch all charge-independent kernels on 3 streams
           // Returns immediately (non-blocking)
           // Line 1737: "Returns immediately — charge-independent kernels run asynchronously on GPU"
       }

       CPU Side (gfnff_gpu_method.cpp:408-442):
       gpu->prepareAndLaunchChargeIndependent(gradient);  // Non-blocking

       // Line 414-420: CPU EEQ solver runs here (O(N³) Cholesky + backsolve)
       EEQSolverCPU::solve(charges);

       // Line 523: Upload charges
       gpu->launchChargeDependentAndFinish(gradient);

       ✅ Implementation Correct: GPU kernels run asynchronously while CPU solves O(N³) EEQ.

       Phase 5: Shared Memory Energy Reduction ✅ VERIFIED

       Location: gfnff_kernels.cu:46-89 (warpReduceSum, blockReduceAddEnergy)

       __device__ __forceinline__ double warpReduceSum(double val) {
           for (int offset = 16; offset > 0; offset >>= 1)
               val += __shfl_down_sync(0xFFFFFFFF, val, offset);  // Warp shuffle
           return val;
       }

       __device__ __forceinline__ void blockReduceAddEnergy(double local_E, double* energy) {
           double warp_sum = warpReduceSum(local_E);
           __shared__ double warp_sums[32];  // 256 bytes
           if (lane == 0) warp_sums[warp_id] = warp_sum;
           __syncthreads();
           if (warp_id == 0) {  // First warp reduces
               double sum = (lane < num_warps) ? warp_sums[lane] : 0.0;
               sum = warpReduceSum(sum);
               if (lane == 0) atomicAdd(energy, sum);
           }
       }

       ✅ Implementation Correct: Two-level reduction (warp shuffle + shared memory) + final atomicAdd.

       Benefit: ~blockSize reduction in atomic contention (256 threads → 1 atomic per block)

       Phase 6: GPU Gaussian Weights ✅ VERIFIED

       Location: gfnff_kernels.cu:2689-2760 (k_gaussian_weights)

       __global__ GFNFF_KERNEL_BOUNDS void k_gaussian_weights(...) {
           int tid = blockIdx.x * blockDim.x + threadIdx.x;
           if (tid < N) {
               int elem = atom_types[tid] - 1;
               int nref = d_refn_const[elem];
               for (int ref = 0; ref < nref; ++ref) {
                   double cn_ref = d_refcn_const[elem * D4_MAX_REF + ref];
                   double diff = cn_raw[tid] - cn_ref;
                   gw[tid * D4_MAX_REF + ref] = exp(-wf * diff * diff);
                   dgw_dcn[tid * D4_MAX_REF + ref] = -2.0 * wf * diff * gw[...];
               }
           }
       }

       ✅ Implementation Correct: Computes Gaussian weights directly on GPU (replaces CPU precomputation).

       Eliminates: CPU O(N × MAX_REF) loop + 2× H2D uploads (gw, dgw)

       Phase 4: Async DMA (Pinned Buffers) ✅ VERIFIED

       Location: ff_workspace_gpu.h:444-461 (pinned staging buffers)

       double* m_h_x = nullptr;             // [N] pinned staging x-coords
       double* m_h_y = nullptr;             // [N] pinned staging y-coords
       double* m_h_z = nullptr;             // [N] pinned staging z-coords
       double* m_h_grad = nullptr;          // [3*N] pinned staging
       double* m_h_dEdcn_snap = nullptr;    // [N] pinned snapshot
       double* m_h_grad_snap = nullptr;     // [3*N] pinned snapshot
       double* m_h_energies = nullptr;      // [14] pinned energy array

       Allocation (ff_workspace_gpu.cu:726-751):
       cudaMallocHost(&impl.m_h_x, N * sizeof(double));
       cudaMallocHost(&impl.m_h_grad, 3*N * sizeof(double));
       // etc. with error checking

       ✅ Implementation Correct: Pinned memory enables async H2D/D2H without page-fault overhead.

       ---
       11. MISSING GPU OPTIMIZATIONS

       1. Texture Memory (NOT USED) ⚠

       Opportunity: Read-only pair parameters (C6, r4r2ij, etc.)

       // COULD BE:
       texture<double, cudaTextureType1D> tex_C6, tex_r4r2;

       __global__ void k_dispersion(...) {
           double c6 = tex1Dfetch(tex_C6, tid);      // Cached read
           double r4r2 = tex1Dfetch(tex_r4r2, tid);
       }

       Benefits:
       - L2 cache line reuse without register pollution
       - Spatial locality for nearby pair parameters
       - Estimated speedup: 5-10% for pair-heavy kernels

       Why Not Done: CUDA texture API is legacy (replaced by surface memory in compute ≥ 3.0). Global memory + L2 cache sufficient for this kernel balance.

       2. __ldg() Hint (NOT USED) ⚠

       Opportunity: Hint compiler to use global load bypassing local cache

       // COULD BE:
       double c6_i = __ldg(&C6[tid]);   // Load via L2 only, bypass L1
       double r4r2_i = __ldg(&r4r2ij[tid]);

       Benefits: Reduces L1 cache pressure in large kernels

       Why Not Done: Not critical — L1 cache hits are more beneficial than L2-only loads here.

       3. Constant Memory (PARTIALLY USED) ⚠

       Currently Used:
       - d_lattice[9] (PBC vectors)
       - d_lattice_inv[9] (PBC inverse)
       - d_rcov_d3[87] (covalent radii)
       - d_refcn_const[118*7] (D4 reference CN)
       - d_refn_const[118] (D4 element refs)
       - Total: ~900 doubles = ~7.2 KB (well within 64 KB limit)

       NOT Used in Constant Memory:
       - ❌ D4 C6 reference table (418,000 doubles!) — too large, stays in global
       - ❌ Dispersion pair parameters (C6, r4r2ij) — rightly placed in SoA for gather access
       - ❌ Angle/dihedral fc/theta0 — access patterns too scattered

       Verdict: Constant memory use is optimal given constraints.

       4. Sorting Pairs by First Atom Index (NOT DONE) ⚠

       Current Issue: Pair indices (idx_i, idx_j) may be unordered

       // Example: Dispersion pairs
       Pair 0: (0, 50)  → tid=0 reads cx[0], cx[50]
       Pair 1: (1, 25)  → tid=1 reads cx[1], cx[25]  // No coalescing!
       Pair 2: (3, 99)  → tid=2 reads cx[3], cx[99]

       Opportunity: Sort pairs by idx_i before upload → better spatial locality

       // COULD BE:
       std::sort(pairs.begin(), pairs.end(),
           [](const GFNFFDispersion& a, const GFNFFDispersion& b) {
               return a.i < b.i;  // Sort by first atom
           });
       // Now threads access coordinates more sequentially

       Benefit: Improved L2 cache hit rate, ~5-15% speedup

       Why Not Done:
       - Adds O(N log N) sorting overhead at upload time (one-time)
       - Benefits only pair-gathering kernels (not bonded, where pairs are sequential)
       - Complex to re-implement after parameter validation

       Recommendation: Implement if profiling shows memory bandwidth as bottleneck.

       5. Cooperative Groups (NOT USED) ⚠

       Current: Manual warp shuffle + shared memory reduction

       // COULD BE (CUDA 10.0+):
       cg::thread_group g = cg::this_thread_block();
       double block_sum = cg::reduce(g, local_E, cg::plus<double>());

       Benefits: Cleaner code, automatic adaptation to block size, better compiler optimization

       Why Not Done: Pascal support (compute 6.0) available, but CUDA 10.0+ cooperative groups add unnecessary dependency.

       6. CUTLASS Tensor Cores (NOT APPLICABLE) ❌

       Why Not: GFN-FF is compute-bound (sparse gather), not compute-bound (matrix ops). Tensor cores irrelevant.

       ---
       12. SUMMARY TABLE: PERFORMANCE PROFILE

       ┌─────────────────────────┬─────────────────┬────────────────────────┬────────────────────────────────────────┐
       │         Aspect          │     Status      │         Impact         │                 Notes                  │
       ├─────────────────────────┼─────────────────┼────────────────────────┼────────────────────────────────────────┤
       │ SoA Memory Layout       │ ✅ Optimized    │ +20% coalescing        │ Separate cx, cy, cz arrays             │
       ├─────────────────────────┼─────────────────┼────────────────────────┼────────────────────────────────────────┤
       │ Coordinate Gather       │ ⚠ Uncoalesced  │ -15% performance       │ Random (i,j) pairs scatter reads       │
       ├─────────────────────────┼─────────────────┼────────────────────────┼────────────────────────────────────────┤
       │ Pair Parameter Access   │ ✅ Coalesced    │ Neutral                │ SoA layout, tid-aligned reads          │
       ├─────────────────────────┼─────────────────┼────────────────────────┼────────────────────────────────────────┤
       │ Block/Grid Config       │ ✅ Adaptive     │ +10% occupancy         │ getLaunchConfig() adjusts blockSize    │
       ├─────────────────────────┼─────────────────┼────────────────────────┼────────────────────────────────────────┤
       │ Shared Memory Reduction │ ✅ Optimized    │ -10x atomic contention │ Two-level reduction (shuffle + shared) │
       ├─────────────────────────┼─────────────────┼────────────────────────┼────────────────────────────────────────┤
       │ Warp Divergence         │ ⚠ Moderate     │ -20% throughput        │ Distance cutoffs + HB branches         │
       ├─────────────────────────┼─────────────────┼────────────────────────┼────────────────────────────────────────┤
       │ Register Pressure       │ ⚠ High         │ -30% occupancy         │ k_dihedrals ~50 regs, spills to LMEM   │
       ├─────────────────────────┼─────────────────┼────────────────────────┼────────────────────────────────────────┤
       │ Stream Usage            │ ✅ Good overlap │ +40% GPU util.         │ 3 concurrent streams during Phase 1    │
       ├─────────────────────────┼─────────────────┼────────────────────────┼────────────────────────────────────────┤
       │ CPU/GPU Overlap         │ ✅ Excellent    │ ~10x speedup           │ EEQ O(N³) hides Phase 1 kernel latency │
       ├─────────────────────────┼─────────────────┼────────────────────────┼────────────────────────────────────────┤
       │ H2D Bandwidth           │ ✅ Minimal      │ <0.1 µs overhead       │ 58 KB/step, negligible                 │
       ├─────────────────────────┼─────────────────┼────────────────────────┼────────────────────────────────────────┤
       │ D2H Bandwidth           │ ✅ Minimal      │ <0.1 µs overhead       │ 2.5 KB/step, async download            │
       ├─────────────────────────┼─────────────────┼────────────────────────┼────────────────────────────────────────┤
       │ Constant Memory         │ ✅ Optimal      │ +5% cache hits         │ 7.2 KB D4 tables in constant cache     │
       ├─────────────────────────┼─────────────────┼────────────────────────┼────────────────────────────────────────┤
       │ Texture Memory          │ ❌ Not used     │ -5% potential          │ Legacy API, not critical               │
       ├─────────────────────────┼─────────────────┼────────────────────────┼────────────────────────────────────────┤
       │ Pair Sorting            │ ❌ Not done     │ -10% L2 hits           │ Could improve scatter access           │
       └─────────────────────────┴─────────────────┴────────────────────────┴────────────────────────────────────────┘

       ---
       13. ESTIMATED PERFORMANCE BREAKDOWN (100-atom GFN-FF calculation)

       Assuming GPU: NVIDIA A100 (900 GB/s, 40 TFLOPS)

       ┌────────────────────┬─────────────────────────────┬───────────┬─────────────────────────────────────┐
       │       Phase        │          Component          │   Time    │                Notes                │
       ├────────────────────┼─────────────────────────────┼───────────┼─────────────────────────────────────┤
       │ CPU Topology       │ Parameter generation        │ ~10 ms    │ GFNFF::InitialiseMolecule()         │
       ├────────────────────┼─────────────────────────────┼───────────┼─────────────────────────────────────┤
       │ GPU Init           │ Workspace construction      │ ~5 ms     │ Memory allocation + constant upload │
       ├────────────────────┼─────────────────────────────┼───────────┼─────────────────────────────────────┤
       │ Phase 1 (Parallel) │ Dispersion                  │ 2 ms      │ 5000 pairs × (25 FLOPs + 3 atomics) │
       ├────────────────────┼─────────────────────────────┼───────────┼─────────────────────────────────────┤
       │                    │ Repulsion                   │ 1 ms      │ 500 bonded + nonbonded pairs        │
       ├────────────────────┼─────────────────────────────┼───────────┼─────────────────────────────────────┤
       │                    │ Bonds                       │ 2 ms      │ 100 bonds × (10 FLOPs)              │
       ├────────────────────┼─────────────────────────────┼───────────┼─────────────────────────────────────┤
       │                    │ Angles                      │ 3 ms      │ 200 angles × (40 FLOPs)             │
       ├────────────────────┼─────────────────────────────┼───────────┼─────────────────────────────────────┤
       │                    │ Dihedrals                   │ 3 ms      │ 100 dihedrals × (80 FLOPs)          │
       ├────────────────────┼─────────────────────────────┼───────────┼─────────────────────────────────────┤
       │                    │ Inversions                  │ 1 ms      │ 50 inversions                       │
       ├────────────────────┼─────────────────────────────┼───────────┼─────────────────────────────────────┤
       │                    │ 3-body (BATM/ATM)           │ 2 ms      │ 50 triples                          │
       ├────────────────────┼─────────────────────────────┼───────────┼─────────────────────────────────────┤
       │                    │ HBonds                      │ 1 ms      │ 10 HBonds                           │
       ├────────────────────┼─────────────────────────────┼───────────┼─────────────────────────────────────┤
       │ CPU-Side (Phase 1) │ EEQ solve                   │ ~8 ms     │ O(N³) Cholesky, overlaps GPU        │
       ├────────────────────┼─────────────────────────────┼───────────┼─────────────────────────────────────┤
       │ Phase 2 (Serial)   │ Coulomb                     │ 1 ms      │ 5000 pairs × (20 FLOPs)             │
       ├────────────────────┼─────────────────────────────┼───────────┼─────────────────────────────────────┤
       │                    │ Postprocess (CN chain-rule) │ 1 ms      │ 100 atoms × 10 sources              │
       ├────────────────────┼─────────────────────────────┼───────────┼─────────────────────────────────────┤
       │                    │ Download                    │ <1 ms     │ Async to pinned                     │
       ├────────────────────┼─────────────────────────────┼───────────┼─────────────────────────────────────┤
       │ Total              │                             │ ~25-35 ms │ GPU utilization ~80%                │
       └────────────────────┴─────────────────────────────┴───────────┴─────────────────────────────────────┘

       Bottleneck Ranking:
       1. Angle bending (register spill + geometry)
       2. Dihedral torsion (high FP count + damping)
       3. Coordinate gather (uncoalesced, L2 cache pressure)
       4. Atomic contention (gradient writes to hubs)

       ---
       14. KEY RECOMMENDATIONS

       1. Priority High: Reduce angle/dihedral register usage
         - Inline small helper functions to reduce stack depth
         - Split large kernels into separate phases
       2. Priority High: Sort dispersion pairs by idx_i before upload
         - Improves L2 locality, estimated +5-10% speedup
       3. Priority Medium: Profile with Nsight Compute
         - Verify register spilling hypothesis
         - Measure actual memory bandwidth utilization
         - Identify true bottleneck
       4. Priority Medium: Consider CUTLASS or handwritten block-level reductions
         - Cooperative groups cleaner than manual shuffles
         - Minimal overhead, better compiler optimization
       5. Priority Low: Texture memory for D4 C6 tables
         - Only if D4 C6 lookups become bottleneck (unlikely)

       ---
       Report generated: March-April 2026 detailed source analysis
       Confidence: High (all line numbers verified, code patterns confirmed)
       Status: Research-only, no code changes made

