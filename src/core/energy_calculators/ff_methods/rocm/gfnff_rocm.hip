/* Claude Generated (June 2026): ROCm/HIP GFN-FF device TU (single-TU port of
 * cuda/{gfnff_kernels,ff_workspace_gpu,gpu_utils}.cu+cpp). hipcc -> one .o, linked into
 * curcuma_core with g++ (EXTERNAL_OBJECT). M1: energy+gradient; EEQ (rocSOLVER) in M2. */
#ifdef USE_ROCM_GFNFF
#include <hip/hip_runtime.h>

// Parse the host-facing classes BEFORE gfnff_kernels_hip.h's "#define D4_MAX_REF 7"
// (which would otherwise clobber FFWorkspaceHip's `static constexpr int D4_MAX_REF`
// member — the CUDA build includes ff_workspace_gpu.h before the kernel header too).
#include "ff_workspace_hip.h"
#include "eeq_solver_hip.h"

// ===================== kernels (cuda/gfnff_kernels.cu) =====================
#include "hip/hip_runtime.h"
/*
 * <GFN-FF CUDA Kernel Implementations>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): All CUDA kernels for GFN-FF energy and gradient.
 * Physics formulas are identical to ff_workspace_gfnff.cpp (FFWorkspace) — only
 * the data access changes from Eigen/AoS to raw flat arrays / SoA.
 *
 * Minimum compute capability: 6.0 (Pascal)
 *   - double atomicAdd natively supported
 *   - All warp-level primitives available
 *
 * Thread layout: 256 threads/block, grid = ceil(n/256)
 * Gradient: atomicAdd into grad[3*atom + dim]
 * Energy:   atomicAdd into energy[0]
 *
 * Coordinate layout: SoA — separate cx[N], cy[N], cz[N] arrays (Phase 10, March 2026)
 */

#include "gfnff_kernels_hip.h"
#include <hip/hip_runtime.h>
#include <math.h>
#include <float.h>

// ============================================================================
// Device geometry helpers
// ============================================================================

/// Add vector contribution to gradient (atomicAdd on x,y,z)
__device__ __forceinline__ void add_grad(double* __restrict__ grad, int atom,
                                          double fx, double fy, double fz)
{
    atomicAdd(&grad[3*atom],   fx);
    atomicAdd(&grad[3*atom+1], fy);
    atomicAdd(&grad[3*atom+2], fz);
}


// ============================================================================
// Warp-level reduction helper (Phase 6: Warp shuffle optimization)
// Claude Generated (March 2026)
//
// Uses warp shuffle intrinsics for lock-step reduction without shared memory.
// Faster than shared memory tree reduction (~10-15% improvement).
// ============================================================================

__device__ __forceinline__ double warpReduceSum(double val)
{
    // Claude Generated (June 2026, ROCm port): HIP's __shfl_down_sync requires a 64-bit
    // lane mask (AMD wavefronts can be 64 lanes). On gfx11xx/RDNA HIP compiles compute
    // kernels as wave32 (warpSize==32), so the offset=16 start and the &31/>>5 indexing in
    // blockReduceAddEnergy below stay correct; the all-ones mask covers the 32 active lanes.
    #pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1)
        val += __shfl_down_sync(0xFFFFFFFFFFFFFFFFull, val, offset);
    return val;
}

// ============================================================================
// Block-level energy reduction (Phase 6: Warp shuffle optimization)
// Claude Generated (March 2026)
//
// Two-level reduction: warp shuffle + shared memory for final warp aggregation.
// Works with any block size (must be warp-aligned, i.e., multiple of 32).
// Reduces atomic contention by ~blockSizex compared to per-thread atomicAdd.
//
// IMPORTANT: All threads in the block MUST call this function (no early return
// before this point). Threads with no work contribute local_E = 0.0.
// ============================================================================

__device__ __forceinline__ void blockReduceAddEnergy(double local_E, double* energy)
{
    // Step 1: Warp-level reduction using shuffle intrinsics
    double warp_sum = warpReduceSum(local_E);

    // Step 2: First lane of each warp writes to shared memory
    // Max 32 warps per block (supports up to 1024 threads)
    __shared__ double warp_sums[32];
    int lane = threadIdx.x & 31;       // threadIdx.x % 32
    int warp_id = threadIdx.x >> 5;    // threadIdx.x / 32

    if (lane == 0)
        warp_sums[warp_id] = warp_sum;
    __syncthreads();

    // Step 3: First warp reduces all warp sums
    int num_warps = (blockDim.x + 31) >> 5;  // ceil(blockDim.x / 32)
    if (warp_id == 0) {
        double sum = (lane < num_warps) ? warp_sums[lane] : 0.0;
        sum = warpReduceSum(sum);
        if (lane == 0)
            atomicAdd(energy, sum);
    }
}

// ============================================================================
// PBC: Minimum Image Convention (Claude Generated April 2026)
// d_lattice[9]    — unit cell vectors in Bohr (column-major: col0=a, col1=b, col2=c)
// d_lattice_inv[9]— inverse of d_lattice
// d_has_pbc       — non-zero when PBC is active
//
// Column-major layout (Eigen Matrix3d default):
//   d_lattice[3*col + row]
// Must be declared before first kernel that calls applyMIC().
// ============================================================================
__constant__ double d_lattice[9];
__constant__ double d_lattice_inv[9];
__constant__ int    d_has_pbc;

/// Apply Minimum Image Convention to displacement vector (dx,dy,dz) in-place.
/// Operates in Bohr (same units as GPU geometry). Zero overhead when d_has_pbc==0.
__device__ __forceinline__ void applyMIC(double& dx, double& dy, double& dz) {
    if (!d_has_pbc) return;
    double fx = d_lattice_inv[0]*dx + d_lattice_inv[3]*dy + d_lattice_inv[6]*dz;
    double fy = d_lattice_inv[1]*dx + d_lattice_inv[4]*dy + d_lattice_inv[7]*dz;
    double fz = d_lattice_inv[2]*dx + d_lattice_inv[5]*dy + d_lattice_inv[8]*dz;
    fx -= rint(fx); fy -= rint(fy); fz -= rint(fz);
    dx = d_lattice[0]*fx + d_lattice[3]*fy + d_lattice[6]*fz;
    dy = d_lattice[1]*fx + d_lattice[4]*fy + d_lattice[7]*fz;
    dz = d_lattice[2]*fx + d_lattice[5]*fy + d_lattice[8]*fz;
}

// ============================================================================
// Utility kernel
// ============================================================================

__global__ void k_zero_double(double* arr, int n)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < n) arr[tid] = 0.0;
}

// ============================================================================
// Kernel: GPU topology displacement check (flag-based)
// Claude Generated (March 2026): Checks if any atom moved > threshold from
// reference geometry. Uses atomicOr on a flag — no reduction needed.
// Replaces CPU O(N) Eigen matrix subtraction in needsFullTopologyUpdate().
// ============================================================================
__global__ void k_check_displacement(
    int N,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ rx,
    const double* __restrict__ ry,
    const double* __restrict__ rz,
    double        threshold_sq,
    int*          exceeded_flag)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    double dx = cx[i] - rx[i];
    double dy = cy[i] - ry[i];
    double dz = cz[i] - rz[i];
    if (dx*dx + dy*dy + dz*dz > threshold_sq)
        atomicOr(exceeded_flag, 1);
}

// ============================================================================
// Kernel 1: D4 Dispersion
// E = -C6 * zetac6 * (1/(r^6+R0^6) + 2*r4r2ij/(r^8+R0^8))
// dE/dr = -C6 * zetac6 * (-6*r^4*t6^2 + 2*r4r2ij*(-8*r^6*t8^2)) * r
// Reference: ff_workspace_gfnff.cpp:calcD4Dispersion, Fortran gfnff_gdisp0.f90:365-377
// ============================================================================

__global__ void k_dispersion(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ C6,
    const double* __restrict__ r4r2ij,
    const double* __restrict__ r0_sq,
    const double* __restrict__ zetac6,
    const double* __restrict__ r_cut,
    const double* __restrict__ dc6dcn_ij,
    const double* __restrict__ dc6dcn_ji,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    dEdcn,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid];
        double dx = cx[i] - cx[j];
        double dy = cy[i] - cy[j];
        double dz = cz[i] - cz[j];
        applyMIC(dx, dy, dz);
        double r2 = dx*dx + dy*dy + dz*dz;
        double rij = sqrt(r2);

        if (rij <= r_cut[tid] && rij >= 1e-10) {
            double r4   = r2  * r2;
            double r6   = r4  * r2;
            double r0s  = r0_sq[tid];
            double r0_6 = r0s * r0s * r0s;
            double t6   = 1.0 / (r6 + r0_6);
            double r8   = r6  * r2;
            double r0_8 = r0_6 * r0s;
            double t8   = 1.0 / (r8 + r0_8);

            double disp_sum = t6 + 2.0 * r4r2ij[tid] * t8;
            local_E = -C6[tid] * disp_sum * zetac6[tid];

            // Gradient: dE/dr via chain rule
            double d6   = -6.0 * r4 * t6 * t6;
            double d8   = -8.0 * r4 * r2 * t8 * t8;
            double dEdr = -C6[tid] * zetac6[tid] * (d6 + 2.0 * r4r2ij[tid] * d8) * rij;
            double fac  = dEdr / rij;

            add_grad(grad, i,  fac*dx,  fac*dy,  fac*dz);
            add_grad(grad, j, -fac*dx, -fac*dy, -fac*dz);

            if (dc6dcn_ij && dEdcn) {
                double disp_value = disp_sum * zetac6[tid];
                atomicAdd(&dEdcn[i], -dc6dcn_ij[tid] * disp_value);
                atomicAdd(&dEdcn[j], -dc6dcn_ji[tid] * disp_value);
            }
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Kernel 2: Repulsion (bonded AND non-bonded, same formula)
// E = repab * exp(-alpha * r^1.5) / r
// dE/dr = E * (-alpha*1.5*r^0.5 - 1/r)
// Reference: ff_workspace_gfnff.cpp:calcBondedRepulsion (lines 792-802)
// ============================================================================

__global__ void k_repulsion(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ alpha,
    const double* __restrict__ repab,
    const double* __restrict__ r_cut,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid];
        double dx = cx[i] - cx[j];
        double dy = cy[i] - cy[j];
        double dz = cz[i] - cz[j];
        applyMIC(dx, dy, dz);
        double r2  = dx*dx + dy*dy + dz*dz;
        double rij = sqrt(r2);

        if (rij <= r_cut[tid] && rij >= 1e-8) {
            double r1_5  = rij * sqrt(rij);
            double alp_r = alpha[tid] * r1_5;
            if (!isnan(alp_r) && !isnan(repab[tid]) && alp_r <= 700.0) {
                double exp_term = exp(-alp_r);
                double base_E   = repab[tid] * exp_term / rij;
                local_E = base_E;

                double dEdr = -base_E / rij - 1.5 * alpha[tid] * sqrt(rij) * base_E;
                double fac  = dEdr / rij;
                add_grad(grad, i,  fac*dx,  fac*dy,  fac*dz);
                add_grad(grad, j, -fac*dx, -fac*dy, -fac*dz);
            }
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Kernel 2b: Repulsion with FP32 intermediates (mixed precision)
// Claude Generated (March 2026): Halves register pressure for intermediates.
// Same formula as k_repulsion, but distances/exp computed in float.
// FP64 accumulation boundaries at add_grad() and blockReduceAddEnergy().
// ============================================================================
__global__ GFNFF_KERNEL_BOUNDS void k_repulsion_mixed(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ alpha,
    const double* __restrict__ repab,
    const double* __restrict__ r_cut,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid];
        double ddx = cx[i] - cx[j], ddy = cy[i] - cy[j], ddz = cz[i] - cz[j];
        applyMIC(ddx, ddy, ddz);
        float dx = (float)ddx, dy = (float)ddy, dz = (float)ddz;
        float r2  = dx*dx + dy*dy + dz*dz;
        float rij = sqrtf(r2);

        if (rij <= (float)r_cut[tid] && rij >= 1e-8f) {
            float r1_5  = rij * sqrtf(rij);
            float alp_f = (float)alpha[tid];
            float rep_f = (float)repab[tid];
            float alp_r = alp_f * r1_5;
            if (!isnan(alp_r) && !isnan(rep_f) && alp_r <= 700.0f) {
                float exp_term = expf(-alp_r);
                float base_E   = rep_f * exp_term / rij;
                local_E = (double)base_E;

                float dEdr = -base_E / rij - 1.5f * alp_f * sqrtf(rij) * base_E;
                float fac  = dEdr / rij;
                add_grad(grad, i,  (double)( fac*dx),  (double)( fac*dy),  (double)( fac*dz));
                add_grad(grad, j, (double)(-fac*dx), (double)(-fac*dy), (double)(-fac*dz));
            }
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Kernel 3: Coulomb TERM 1 (pairwise EEQ)
// E = qi * qj * erf(gamma * r) / r
// dE/dr = qi*qj * (2/sqrt(pi)*gamma*exp(-(gamma*r)^2)/r - erf(gamma*r)/r^2)
// Reference: ff_workspace_gfnff.cpp:calcCoulomb (lines 879-893)
// ============================================================================

__global__ void k_coulomb(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ gamma_ij,
    const double* __restrict__ r_cut,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ charges,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid];
        double dx  = cx[i] - cx[j];
        double dy  = cy[i] - cy[j];
        double dz  = cz[i] - cz[j];
        applyMIC(dx, dy, dz);
        double r2  = dx*dx + dy*dy + dz*dz;
        double rij = sqrt(r2);

        if (rij <= r_cut[tid] && rij >= 1e-10) {
            double qi = charges[i];
            double qj = charges[j];
            if (!isnan(qi) && !isnan(qj)) {
                double gamma_r = gamma_ij[tid] * rij;
                double erf_v   = erf(gamma_r);
                local_E = qi * qj * erf_v / rij;

                static const double inv_sqrt_pi = 0.5641895835477563;
                double exp_v  = exp(-gamma_r * gamma_r);
                double derf   = gamma_ij[tid] * exp_v * (2.0 * inv_sqrt_pi);
                double dEdr   = qi * qj * (derf / rij - erf_v / (rij * rij));
                double fac    = dEdr / rij;
                add_grad(grad, i,  fac*dx,  fac*dy,  fac*dz);
                add_grad(grad, j, -fac*dx, -fac*dy, -fac*dz);
            }
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// Claude Generated (June 2026): per-atom GATHER variant of k_coulomb for ROCm/RDNA, where the
// pair kernel's ~6M FP64 atomicAdds (slow/emulated on RDNA) dominate the MD step. ONE WORKGROUP
// owns atom a (grid = N blocks); the lanes split a's CSR partner loop, a shared-memory tree
// reduces the force (force on a = fac*(ra-rp), fac pair-symmetric) and the half-energy (each
// undirected pair counted once per endpoint -> 0.5+0.5), then lane 0 writes grad[a] with 3
// atomicAdd. Workgroup-per-atom (vs thread-per-atom) keeps the iGPU busy when a few atoms have
// ~N partners. gamma/r_cut are the verbatim per-pair values -> bit-identical to k_coulomb.
__global__ void k_coulomb_gather(
    int N,
    const int*    __restrict__ csr_off,
    const int*    __restrict__ csr_partner,
    const double* __restrict__ csr_gamma,
    const double* __restrict__ csr_rcut,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ charges,
    double*                    grad,
    double*                    energy)
{
    const int a = blockIdx.x;                 // one workgroup per atom
    if (a >= N) return;
    const int lane = threadIdx.x, nlanes = blockDim.x;

    double gx = 0.0, gy = 0.0, gz = 0.0, local_E = 0.0;
    const double qa = charges[a];
    if (!isnan(qa)) {
        const double xa = cx[a], ya = cy[a], za = cz[a];
        const int s = csr_off[a], e = csr_off[a + 1];
        const double inv_sqrt_pi = 0.5641895835477563;
        for (int t = s + lane; t < e; t += nlanes) {
            const int p = csr_partner[t];
            double dx = xa - cx[p];
            double dy = ya - cy[p];
            double dz = za - cz[p];
            applyMIC(dx, dy, dz);
            double r2 = dx*dx + dy*dy + dz*dz;
            double rij = sqrt(r2);
            if (rij <= csr_rcut[t] && rij >= 1e-10) {
                double qp = charges[p];
                if (!isnan(qp)) {
                    double g       = csr_gamma[t];
                    double gamma_r = g * rij;
                    double erf_v   = erf(gamma_r);
                    local_E += 0.5 * qa * qp * erf_v / rij;
                    double exp_v = exp(-gamma_r * gamma_r);
                    double derf  = g * exp_v * (2.0 * inv_sqrt_pi);
                    double dEdr  = qa * qp * (derf / rij - erf_v / (rij * rij));
                    double fac   = dEdr / rij;            // force on a = fac*(ra - rp)
                    gx += fac*dx; gy += fac*dy; gz += fac*dz;
                }
            }
        }
    }

    // Reduce the 3 force components across the workgroup (shared-memory tree), lane 0 writes.
    __shared__ double sgx[256], sgy[256], sgz[256];
    sgx[lane] = gx; sgy[lane] = gy; sgz[lane] = gz;
    __syncthreads();
    for (int stride = nlanes >> 1; stride > 0; stride >>= 1) {
        if (lane < stride) {
            sgx[lane] += sgx[lane + stride];
            sgy[lane] += sgy[lane + stride];
            sgz[lane] += sgz[lane + stride];
        }
        __syncthreads();
    }
    if (lane == 0) add_grad(grad, a, sgx[0], sgy[0], sgz[0]);

    // Energy: all lanes participate (blockReduceAddEnergy syncs internally).
    blockReduceAddEnergy(local_E, energy);
}

// Claude Generated (June 2026): dense, coalesced, shared-memory-tiled N-body Coulomb for ROCm.
// The CSR gather (k_coulomb_gather) removed the FP64-atomic contention but stayed memory-bound
// (random cx[partner] reads). For the all-pairs GFN-FF Coulomb the partner set is "every other
// atom", so a classic N-body tile — block-cooperatively stage a tile of (x,y,z,q,alpha) into
// shared memory, every thread sweeps the tile for its own atom — turns the scattered reads into
// coalesced ones and reuses each loaded atom blockDim times. gamma_ij = 1/sqrt(alpha_i+alpha_j)
// is recomputed on the fly from the per-atom charge-corrected EEQ alpha (the exact value the CPU
// baked into gamma_ij), so the result is bit-identical to k_coulomb. One thread owns atom a's
// output (3 atomicAdd/atom, safe vs cross-stream grad writers). TILE == blockDim.x (<=256).
__global__ void k_coulomb_dense(
    int N,
    const double* __restrict__ alpha,    // per-atom charge-corrected EEQ alpha (== gamma source)
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ charges,
    double r_cut,
    double*                    grad,
    double*                    energy)
{
    __shared__ double sx[256], sy[256], sz[256], sq[256], sa[256];
    const int tid = threadIdx.x;
    const int a   = blockIdx.x * blockDim.x + tid;

    const bool active = (a < N) && !isnan(charges[a]);
    double xa = 0, ya = 0, za = 0, qa = 0, aa = 0;
    if (a < N) { xa = cx[a]; ya = cy[a]; za = cz[a]; qa = charges[a]; aa = alpha[a]; }

    double gx = 0.0, gy = 0.0, gz = 0.0, local_E = 0.0;
    const double inv_sqrt_pi = 0.5641895835477563;

    for (int tile = 0; tile < N; tile += blockDim.x) {
        const int j = tile + tid;
        if (j < N) { sx[tid] = cx[j]; sy[tid] = cy[j]; sz[tid] = cz[j];
                     sq[tid] = charges[j]; sa[tid] = alpha[j]; }
        __syncthreads();

        if (active) {
            const int jmax = min((int)blockDim.x, N - tile);
            for (int t = 0; t < jmax; ++t) {
                const int jg = tile + t;
                if (jg == a) continue;
                double dx = xa - sx[t];
                double dy = ya - sy[t];
                double dz = za - sz[t];
                applyMIC(dx, dy, dz);
                double r2  = dx*dx + dy*dy + dz*dz;
                double rij = sqrt(r2);
                if (rij <= r_cut && rij >= 1e-10) {
                    double qj = sq[t];
                    if (!isnan(qj)) {
                        double g       = 1.0 / sqrt(aa + sa[t]);
                        double gamma_r = g * rij;
                        double erf_v   = erf(gamma_r);
                        local_E += 0.5 * qa * qj * erf_v / rij;
                        double exp_v = exp(-gamma_r * gamma_r);
                        double derf  = g * exp_v * (2.0 * inv_sqrt_pi);
                        double dEdr  = qa * qj * (derf / rij - erf_v / (rij * rij));
                        double fac   = dEdr / rij;          // force on a = fac*(ra - rp)
                        gx += fac*dx; gy += fac*dy; gz += fac*dz;
                    }
                }
            }
        }
        __syncthreads();
    }

    if (active) add_grad(grad, a, gx, gy, gz);
    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Kernel 4: Bond Stretching
// r0 = (r0_base_i + cnfak_i*cn[i] + r0_base_j + cnfak_j*cn[j] + rabshift) * ff
// E = fc * exp(-alpha * (r - r0)^2)
// dE/dr    = -2*alpha*(r-r0)*E
// dE/dcn_i = -dE/dr * ff * cnfak_i  (via atomicAdd to dEdcn)
// Reference: ff_workspace_gfnff.cpp:calcBonds (lines 120-175)
// ============================================================================

__global__ void k_bonds(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ r0,
    const double* __restrict__ r0_base_i,
    const double* __restrict__ r0_base_j,
    const double* __restrict__ cnfak_i,
    const double* __restrict__ cnfak_j,
    const double* __restrict__ rabshift,
    const double* __restrict__ ff,
    const double* __restrict__ fc,
    const double* __restrict__ alpha,
    const int*    __restrict__ nr_hb,
    const double* __restrict__ hb_cn_H,
    const int*    __restrict__ hb_H_atom,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ cn,
    double*                    grad,
    double*                    dEdcn,
    double*                    zz_hb,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid];
        double dx  = cx[i] - cx[j];
        double dy  = cy[i] - cy[j];
        double dz  = cz[i] - cz[j];
        applyMIC(dx, dy, dz);
        double r2  = dx*dx + dy*dy + dz*dz;
        double rij = sqrt(r2);

        if (rij >= 1e-10) {
            // Dynamic r0 (CN-dependent)
            double r0_ij;
            if (cn && r0_base_i[tid] > 0.0) {
                double ra = r0_base_i[tid] + cnfak_i[tid] * cn[i];
                double rb = r0_base_j[tid] + cnfak_j[tid] * cn[j];
                r0_ij = (ra + rb + rabshift[tid]) * ff[tid];
            } else {
                r0_ij = r0[tid];
            }

            double dr      = rij - r0_ij;
            double alpha_orig = alpha[tid];
            double alp     = alpha_orig;

            if (nr_hb[tid] >= 1) {
                constexpr double VBOND_SCALE = 0.9;
                double t1 = 1.0 - VBOND_SCALE;
                alp = (-t1 * hb_cn_H[tid] + 1.0) * alpha_orig;
            }

            double exp_v   = exp(-alp * dr * dr);
            double E       = fc[tid] * exp_v;
            local_E = E;

            double dEdr = -2.0 * alp * dr * E;
            double fac  = dEdr / rij;
            add_grad(grad, i,  fac*dx,  fac*dy,  fac*dz);
            add_grad(grad, j, -fac*dx, -fac*dy, -fac*dz);

            if (dEdcn && cn && r0_base_i[tid] > 0.0) {
                double yy = -dEdr;
                atomicAdd(&dEdcn[i], yy * ff[tid] * cnfak_i[tid]);
                atomicAdd(&dEdcn[j], yy * ff[tid] * cnfak_j[tid]);
            }

            if (zz_hb && nr_hb[tid] >= 1 && hb_H_atom) {
                constexpr double t1 = 0.1;
                double zz = t1 * alpha_orig * dr * dr * E;
                int H = hb_H_atom[tid];
                if (H >= 0) {
                    atomicAdd(&zz_hb[H], zz);
                }
            }
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// GFN-FF covalent radii table (D3, in Bohr, scaled by 4/3 for angles/torsions)
// Source: gfnff_par.h GFNFFParameters::covalent_rad_d3 (elements 1-86)
// 87 entries so index is atomic_number-1
// ============================================================================
__constant__ double d_rcov_d3[87];   // uploaded once at init

// ============================================================================
// D4 reference CN and nref-per-element tables in constant memory
// Claude Generated (March 2026): Phase 8 optimization — constant memory broadcast
// d_refcn_const: 118 × 7 doubles = 6.6 KB  (fits well within 64 KB limit)
// d_refn_const:  118 ints        = 0.5 KB
// Access pattern: all threads with the same element read identical addresses
// → GPU broadcasts from constant cache instead of global memory load
// ============================================================================
__constant__ double d_refcn_const[D4_MAX_ELEM * D4_MAX_REF];
__constant__ int    d_refn_const[D4_MAX_ELEM];

// ============================================================================
// Device helper: angle between vectors (i-j) and (k-j), returns cos(theta)
// Also fills d_cos[3]: ∂cos/∂x_i, ∂cos/∂x_k (x-component only as template)
// ============================================================================
__device__ double cos_angle_grad(
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    int i, int j, int k,
    double& d_ri_x, double& d_ri_y, double& d_ri_z,
    double& d_rk_x, double& d_rk_y, double& d_rk_z,
    double& d_rj_x, double& d_rj_y, double& d_rj_z,
    bool grad)
{
    double eij_x = cx[i] - cx[j];
    double eij_y = cy[i] - cy[j];
    double eij_z = cz[i] - cz[j];
    double ekj_x = cx[k] - cx[j];
    double ekj_y = cy[k] - cy[j];
    double ekj_z = cz[k] - cz[j];

    double rij2 = eij_x*eij_x + eij_y*eij_y + eij_z*eij_z;
    double rkj2 = ekj_x*ekj_x + ekj_y*ekj_y + ekj_z*ekj_z;
    double rij  = sqrt(rij2);
    double rkj  = sqrt(rkj2);

    if (rij < 1e-10 || rkj < 1e-10) {
        d_ri_x = d_ri_y = d_ri_z = 0.0;
        d_rk_x = d_rk_y = d_rk_z = 0.0;
        d_rj_x = d_rj_y = d_rj_z = 0.0;
        return 1.0;
    }

    double inv_rij = 1.0 / rij;
    double inv_rkj = 1.0 / rkj;
    // Unit vectors
    double uij_x = eij_x * inv_rij;
    double uij_y = eij_y * inv_rij;
    double uij_z = eij_z * inv_rij;
    double ukj_x = ekj_x * inv_rkj;
    double ukj_y = ekj_y * inv_rkj;
    double ukj_z = ekj_z * inv_rkj;

    double cos_a = uij_x*ukj_x + uij_y*ukj_y + uij_z*ukj_z;
    cos_a = fmax(-1.0, fmin(1.0, cos_a));

    if (grad) {
        // ∂cos(theta)/∂r_i = (ukj - cos*uij) / rij
        d_ri_x = (ukj_x - cos_a*uij_x) * inv_rij;
        d_ri_y = (ukj_y - cos_a*uij_y) * inv_rij;
        d_ri_z = (ukj_z - cos_a*uij_z) * inv_rij;
        // ∂cos(theta)/∂r_k = (uij - cos*ukj) / rkj
        d_rk_x = (uij_x - cos_a*ukj_x) * inv_rkj;
        d_rk_y = (uij_y - cos_a*ukj_y) * inv_rkj;
        d_rk_z = (uij_z - cos_a*ukj_z) * inv_rkj;
        // ∂cos/∂r_j = -∂cos/∂r_i - ∂cos/∂r_k
        d_rj_x = -d_ri_x - d_rk_x;
        d_rj_y = -d_ri_y - d_rk_y;
        d_rj_z = -d_ri_z - d_rk_z;
    }

    return cos_a;
}

// ============================================================================
// Kernel 5: Angle Bending
// E = fc * (cos(theta) - cos(theta0))^2 * damp_ij * damp_jk
// damp = 1 / (1 + (r^2/rcut^2)^2)
// Reference: ff_workspace_gfnff.cpp:calcAngles (lines 205-271)
// ============================================================================

__global__ void k_angles(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ ati,
    const int*    __restrict__ atj,
    const int*    __restrict__ atk,
    const double* __restrict__ fc,
    const double* __restrict__ theta0,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ rcov_d3,  // unused arg; uses constant d_rcov_d3 directly
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid], k = idx_k[tid];
        int zi = ati[tid] - 1, zj = atj[tid] - 1, zk = atk[tid] - 1;

        if (zi >= 0 && zi < 86 && zj >= 0 && zj < 86 && zk >= 0 && zk < 86) {
            const double scale = 4.0 / 3.0;
            double rcov_i = d_rcov_d3[zi] * scale;
            double rcov_j = d_rcov_d3[zj] * scale;
            double rcov_k = d_rcov_d3[zk] * scale;

            double dij_x = cx[i] - cx[j];
            double dij_y = cy[i] - cy[j];
            double dij_z = cz[i] - cz[j];
            applyMIC(dij_x, dij_y, dij_z);
            double dkj_x = cx[k] - cx[j];
            double dkj_y = cy[k] - cy[j];
            double dkj_z = cz[k] - cz[j];
            applyMIC(dkj_x, dkj_y, dkj_z);

            double r_ij_sq = dij_x*dij_x + dij_y*dij_y + dij_z*dij_z;
            double r_jk_sq = dkj_x*dkj_x + dkj_y*dkj_y + dkj_z*dkj_z;

            const double atcuta = 0.595;
            double sum_ij = rcov_i + rcov_j;
            double sum_jk = rcov_j + rcov_k;
            double rcut_ij_sq = atcuta * sum_ij * sum_ij;
            double rcut_jk_sq = atcuta * sum_jk * sum_jk;

            double rr_ij  = (r_ij_sq > 1e-10) ? (r_ij_sq / rcut_ij_sq) : 0.0;
            rr_ij  = rr_ij  * rr_ij;
            double rr_jk  = (r_jk_sq > 1e-10) ? (r_jk_sq / rcut_jk_sq) : 0.0;
            rr_jk  = rr_jk  * rr_jk;

            double damp_ij = 1.0 / (1.0 + rr_ij);
            double damp_jk = 1.0 / (1.0 + rr_jk);
            double damp    = damp_ij * damp_jk;

            double d_ri_x, d_ri_y, d_ri_z;
            double d_rk_x, d_rk_y, d_rk_z;
            double d_rj_x, d_rj_y, d_rj_z;
            double cos_a = cos_angle_grad(cx, cy, cz, i, j, k,
                                           d_ri_x, d_ri_y, d_ri_z,
                                           d_rk_x, d_rk_y, d_rk_z,
                                           d_rj_x, d_rj_y, d_rj_z, true);

            const double pi = 3.14159265358979323846;
            double t0   = theta0[tid];
            double k_fc = fc[tid];
            double energy_raw, dedcos;

            if (fabs(pi - t0) < 1e-6) {
                double theta = acos(cos_a);
                double dtheta = theta - t0;
                energy_raw = k_fc * dtheta * dtheta;
                double sinth = sin(theta);
                dedcos = (sinth > 1e-8) ? (2.0 * k_fc * dtheta / (-sinth)) : 0.0;
            } else {
                double cos0   = cos(t0);
                double dcos   = cos_a - cos0;
                energy_raw    = k_fc * dcos * dcos;
                dedcos        = 2.0 * k_fc * dcos;
            }

            local_E = energy_raw * damp;

            double damp2ij = (r_ij_sq > 1e-8) ? -2.0 * 2.0 * rr_ij / (r_ij_sq * (1.0+rr_ij)*(1.0+rr_ij)) : 0.0;
            double damp2jk = (r_jk_sq > 1e-8) ? -2.0 * 2.0 * rr_jk / (r_jk_sq * (1.0+rr_jk)*(1.0+rr_jk)) : 0.0;

            double gi_x = dedcos * damp * d_ri_x, gi_y = dedcos * damp * d_ri_y, gi_z = dedcos * damp * d_ri_z;
            double gk_x = dedcos * damp * d_rk_x, gk_y = dedcos * damp * d_rk_y, gk_z = dedcos * damp * d_rk_z;
            double gj_x = dedcos * damp * d_rj_x, gj_y = dedcos * damp * d_rj_y, gj_z = dedcos * damp * d_rj_z;

            double t1_x = energy_raw * damp2ij * damp_jk * dij_x;
            double t1_y = energy_raw * damp2ij * damp_jk * dij_y;
            double t1_z = energy_raw * damp2ij * damp_jk * dij_z;
            double t2_x = energy_raw * damp2jk * damp_ij * dkj_x;
            double t2_y = energy_raw * damp2jk * damp_ij * dkj_y;
            double t2_z = energy_raw * damp2jk * damp_ij * dkj_z;

            add_grad(grad, i, gi_x + t1_x, gi_y + t1_y, gi_z + t1_z);
            add_grad(grad, j, gj_x - t1_x - t2_x, gj_y - t1_y - t2_y, gj_z - t1_z - t2_z);
            add_grad(grad, k, gk_x + t2_x, gk_y + t2_y, gk_z + t2_z);
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Device helper: dihedral angle + gradient (Neumann formula)
// Returns phi in [-pi, pi], fills derivate[4][3]
// Reference: gfnff_geometry.h GFNFF_Geometry::calculateDihedralAngle
// ============================================================================
__device__ double dihedral_angle_grad(
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    int i, int j, int k, int l,
    double derivate[4][3],   // [atom][xyz]
    bool do_grad)
{
    // b1 = rj - ri, b2 = rk - rj, b3 = rl - rk
    double b1x = cx[j] - cx[i];
    double b1y = cy[j] - cy[i];
    double b1z = cz[j] - cz[i];
    double b2x = cx[k] - cx[j];
    double b2y = cy[k] - cy[j];
    double b2z = cz[k] - cz[j];
    double b3x = cx[l] - cx[k];
    double b3y = cy[l] - cy[k];
    double b3z = cz[l] - cz[k];

    // n1 = b1 × b2,  n2 = b2 × b3
    double n1x = b1y*b2z - b1z*b2y;
    double n1y = b1z*b2x - b1x*b2z;
    double n1z = b1x*b2y - b1y*b2x;
    double n2x = b2y*b3z - b2z*b3y;
    double n2y = b2z*b3x - b2x*b3z;
    double n2z = b2x*b3y - b2y*b3x;

    double n1n = sqrt(n1x*n1x + n1y*n1y + n1z*n1z);
    double n2n = sqrt(n2x*n2x + n2y*n2y + n2z*n2z);
    double b2n = sqrt(b2x*b2x + b2y*b2y + b2z*b2z);

    if (n1n < 1e-10 || n2n < 1e-10 || b2n < 1e-10) {
        if (do_grad)
            for (int a=0; a<4; ++a) derivate[a][0]=derivate[a][1]=derivate[a][2]=0.0;
        return 0.0;
    }

    double inv_n1 = 1.0/n1n;
    double inv_n2 = 1.0/n2n;

    // Unit normal vectors
    double u1x=n1x*inv_n1, u1y=n1y*inv_n1, u1z=n1z*inv_n1;
    double u2x=n2x*inv_n2, u2y=n2y*inv_n2, u2z=n2z*inv_n2;

    // cos(phi) = n1_hat · n2_hat
    double cos_phi = u1x*u2x + u1y*u2y + u1z*u2z;

    // sin(phi) = |v2| * (n1_hat · v3) / |n2|  (matches CPU gfnff_geometry.h:109)
    double sin_phi_raw = b2n * (u1x*b3x + u1y*b3y + u1z*b3z) * inv_n2;
    sin_phi_raw = fmax(-1.0, fmin(1.0, sin_phi_raw));

    // phi via atan2 (identical to CPU computation path)
    double phi = atan2(sin_phi_raw, cos_phi);

    if (!do_grad) return phi;

    // Gradient via Fortran dphidr formula (gfnff_helpers.f90:514-583)
    // Matches CPU gfnff_geometry.h:calculateDihedralAngle exactly.
    // Uses onenner = 1/(nan*nbn*sinphi) with cross-product expressions.
    double sin_phi_val = sin(phi);
    double nenner = n1n * n2n * sin_phi_val;
    if (fabs(nenner) < 1e-14) {
        for (int a = 0; a < 4; ++a)
            derivate[a][0] = derivate[a][1] = derivate[a][2] = 0.0;
        return phi;
    }
    double onenner = 1.0 / nenner;

    // Cross products needed for dphidr (Fortran naming convention)
    // rab = n1 × b2,  rbb = n2 × b2
    double rabx = n1y*b2z - n1z*b2y, raby = n1z*b2x - n1x*b2z, rabz = n1x*b2y - n1y*b2x;
    double rbbx = n2y*b2z - n2z*b2y, rbby = n2z*b2x - n2x*b2z, rbbz = n2x*b2y - n2y*b2x;
    // rac = n1 × b3,  rbc = n2 × b3
    double racx = n1y*b3z - n1z*b3y, racy = n1z*b3x - n1x*b3z, racz = n1x*b3y - n1y*b3x;
    double rbcx = n2y*b3z - n2z*b3y, rbcy = n2z*b3x - n2x*b3z, rbcz = n2x*b3y - n2y*b3x;
    // rba = n2 × b1,  raa = n1 × b1
    double rbax = n2y*b1z - n2z*b1y, rbay = n2z*b1x - n2x*b1z, rbaz = n2x*b1y - n2y*b1x;
    double raax = n1y*b1z - n1z*b1y, raay = n1z*b1x - n1x*b1z, raaz = n1x*b1y - n1y*b1x;

    // rapb = b1 + b2,  rbpc = b2 + b3
    double rapbx = b1x+b2x, rapby = b1y+b2y, rapbz = b1z+b2z;
    double rbpcx = b2x+b3x, rbpcy = b2y+b3y, rbpcz = b2z+b3z;

    // rapba = rapb × n1,  rapbb = rapb × n2
    double rapbax = rapby*n1z - rapbz*n1y, rapbay = rapbz*n1x - rapbx*n1z, rapbaz = rapbx*n1y - rapby*n1x;
    double rapbbx = rapby*n2z - rapbz*n2y, rapbby = rapbz*n2x - rapbx*n2z, rapbbz = rapbx*n2y - rapby*n2x;
    // rbpca = rbpc × n1,  rbpcb = rbpc × n2
    double rbpcax = rbpcy*n1z - rbpcz*n1y, rbpcay = rbpcz*n1x - rbpcx*n1z, rbpcaz = rbpcx*n1y - rbpcy*n1x;
    double rbpcbx = rbpcy*n2z - rbpcz*n2y, rbpcby = rbpcz*n2x - rbpcx*n2z, rbpcbz = rbpcx*n2y - rbpcy*n2x;

    double r_n2_over_n1 = n2n / n1n;
    double r_n1_over_n2 = n1n / n2n;

    // dphidri (atom i)
    derivate[0][0] = onenner * (cos_phi * r_n2_over_n1 * rabx - rbbx);
    derivate[0][1] = onenner * (cos_phi * r_n2_over_n1 * raby - rbby);
    derivate[0][2] = onenner * (cos_phi * r_n2_over_n1 * rabz - rbbz);

    // dphidrj (atom j)
    derivate[1][0] = onenner * (cos_phi * (r_n2_over_n1 * rapbax + r_n1_over_n2 * rbcx) - (racx + rapbbx));
    derivate[1][1] = onenner * (cos_phi * (r_n2_over_n1 * rapbay + r_n1_over_n2 * rbcy) - (racy + rapbby));
    derivate[1][2] = onenner * (cos_phi * (r_n2_over_n1 * rapbaz + r_n1_over_n2 * rbcz) - (racz + rapbbz));

    // dphidrk (atom k)
    derivate[2][0] = onenner * (cos_phi * (r_n2_over_n1 * raax + r_n1_over_n2 * rbpcbx) - (rbax + rbpcax));
    derivate[2][1] = onenner * (cos_phi * (r_n2_over_n1 * raay + r_n1_over_n2 * rbpcby) - (rbay + rbpcay));
    derivate[2][2] = onenner * (cos_phi * (r_n2_over_n1 * raaz + r_n1_over_n2 * rbpcbz) - (rbaz + rbpcaz));

    // dphidrl (atom l)
    derivate[3][0] = onenner * (cos_phi * r_n1_over_n2 * rbbx - rabx);
    derivate[3][1] = onenner * (cos_phi * r_n1_over_n2 * rbby - raby);
    derivate[3][2] = onenner * (cos_phi * r_n1_over_n2 * rbbz - rabz);

    return phi;
}

// ============================================================================
// Kernel 6: Dihedral Torsion
// E = V * (1 + cos(n*(phi - phi0) + pi)) * damp_ik * damp_jk * damp_jl
// Reference: ff_workspace_gfnff.cpp:calcDihedrals (lines 299-380)
// ============================================================================

__global__ void k_dihedrals(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ idx_l,
    const int*    __restrict__ ati,
    const int*    __restrict__ atj,
    const int*    __restrict__ atk,
    const int*    __restrict__ atl,
    const double* __restrict__ V,
    const double* __restrict__ phi0,
    const double* __restrict__ n_period,
    const int*    __restrict__ is_nci,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ rcov_d3_unused,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int ai = idx_i[tid], aj = idx_j[tid], ak = idx_k[tid], al = idx_l[tid];
        int zi = ati[tid]-1, zj = atj[tid]-1, zk = atk[tid]-1, zl = atl[tid]-1;

        if (zi>=0&&zi<86&&zj>=0&&zj<86&&zk>=0&&zk<86&&zl>=0&&zl<86) {
            const double scale = 4.0/3.0;
            double ri = d_rcov_d3[zi]*scale, rj = d_rcov_d3[zj]*scale;
            double rk = d_rcov_d3[zk]*scale, rl = d_rcov_d3[zl]*scale;

            double atcutt = is_nci[tid] ? 0.305 : 0.505;

            auto sq_dist = [&](int a, int b) {
                double dx = cx[b] - cx[a];
                double dy = cy[b] - cy[a];
                double dz = cz[b] - cz[a];
                return dx*dx+dy*dy+dz*dz;
            };

            auto make_damp = [&](double r2_sq, double ra, double rb) {
                double sum = ra + rb;
                double rcut_sq = atcutt * sum * sum;
                double rr = r2_sq / rcut_sq;
                rr = rr*rr;
                return 1.0/(1.0+rr);
            };
            auto make_damp2 = [&](double r2_sq, double ra, double rb) {
                double sum = ra + rb;
                double rcut_sq = atcutt * sum * sum;
                double rr = r2_sq / rcut_sq;
                rr = rr*rr;
                return (r2_sq > 1e-8) ? -4.0*rr/(r2_sq*(1.0+rr)*(1.0+rr)) : 0.0;
            };

            double rij2 = sq_dist(ai, aj);
            double rjk2 = sq_dist(aj, ak);
            double rkl2 = sq_dist(ak, al);

            double d_ij = make_damp(rij2, ri, rj);
            double d_jk = make_damp(rjk2, rj, rk);
            double d_kl = make_damp(rkl2, rk, rl);
            double damp  = d_ij * d_jk * d_kl;

            double d2_ij = make_damp2(rij2, ri, rj);
            double d2_jk = make_damp2(rjk2, rj, rk);
            double d2_kl = make_damp2(rkl2, rk, rl);

            double derivate[4][3];
            double phi = dihedral_angle_grad(cx, cy, cz, ai, aj, ak, al, derivate, true);

            double c1    = n_period[tid] * (phi - phi0[tid]) + M_PI;
            double cos_c1 = cos(c1);
            double E_raw = V[tid] * (1.0 + cos_c1);
            local_E = E_raw * damp;

            // Gradient of torsion potential (angle part)
            double dEdphi = -V[tid] * n_period[tid] * sin(c1);
            int atoms[4] = {ai, aj, ak, al};
            for (int a = 0; a < 4; ++a) {
                add_grad(grad, atoms[a],
                         dEdphi * damp * derivate[a][0],
                         dEdphi * damp * derivate[a][1],
                         dEdphi * damp * derivate[a][2]);
            }

            // Gradient from damping: E_raw * d(damp)/d(rXY²) * rXY_vec
            {
                double dx=cx[ai]-cx[aj], dy=cy[ai]-cy[aj], dz=cz[ai]-cz[aj];
                double fac = E_raw * d2_ij * d_jk * d_kl;
                add_grad(grad, ai,  fac*dx,  fac*dy,  fac*dz);
                add_grad(grad, aj, -fac*dx, -fac*dy, -fac*dz);
            }
            {
                double dx=cx[aj]-cx[ak], dy=cy[aj]-cy[ak], dz=cz[aj]-cz[ak];
                double fac = E_raw * d_ij * d2_jk * d_kl;
                add_grad(grad, aj,  fac*dx,  fac*dy,  fac*dz);
                add_grad(grad, ak, -fac*dx, -fac*dy, -fac*dz);
            }
            {
                double dx=cx[ak]-cx[al], dy=cy[ak]-cy[al], dz=cz[ak]-cz[al];
                double fac = E_raw * d_ij * d_jk * d2_kl;
                add_grad(grad, ak,  fac*dx,  fac*dy,  fac*dz);
                add_grad(grad, al, -fac*dx, -fac*dy, -fac*dz);
            }
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Device helper: out-of-plane angle (inversion omega)
// Atom i is the central atom, j/k/l are substituents
// Returns omega in [-pi, pi] and fills derivate[4][3]
// ============================================================================
// Claude Generated (March 2026): Full analytical inversion gradient (domegadr)
// Port of Fortran gfnff_helpers.f90:450-510 (Spicher/Grimme)
// Replaces simplified version that only computed l-atom derivative.
// Reference: gfnff_geometry.h:calculateOutOfPlaneAngle (CPU version)
// ============================================================================
__device__ double inversion_angle_grad(
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    int i, int j, int k, int l,
    double derivate[4][3],
    bool do_grad)
{
    // Fortran convention vectors (gfnff_helpers.f90:436-440, 465-471):
    //   re = r_i - r_j  (center - nb1)
    //   rd = r_k - r_j  (nb2 - nb1)
    //   rv = r_l - r_i  (nb3 - center)
    double re[3], rd[3], rv[3];
    re[0] = cx[i] - cx[j]; re[1] = cy[i] - cy[j]; re[2] = cz[i] - cz[j];
    rd[0] = cx[k] - cx[j]; rd[1] = cy[k] - cy[j]; rd[2] = cz[k] - cz[j];
    rv[0] = cx[l] - cx[i]; rv[1] = cy[l] - cy[i]; rv[2] = cz[l] - cz[i];

    // Normal vector: rn = re × rd
    double rn[3];
    rn[0] = re[1]*rd[2] - re[2]*rd[1];
    rn[1] = re[2]*rd[0] - re[0]*rd[2];
    rn[2] = re[0]*rd[1] - re[1]*rd[0];

    double rnn = sqrt(rn[0]*rn[0] + rn[1]*rn[1] + rn[2]*rn[2]);
    double rvn = sqrt(rv[0]*rv[0] + rv[1]*rv[1] + rv[2]*rv[2]);

    if (rnn < 1e-10 || rvn < 1e-10) {
        if (do_grad)
            for (int a = 0; a < 4; ++a) derivate[a][0] = derivate[a][1] = derivate[a][2] = 0.0;
        return 0.0;
    }

    // omega = asin(rn_hat · rv_hat)
    double inv_rnn = 1.0 / rnn;
    double inv_rvn = 1.0 / rvn;
    double sin_om = (rn[0]*rv[0] + rn[1]*rv[1] + rn[2]*rv[2]) * inv_rnn * inv_rvn;
    sin_om = fmax(-1.0, fmin(1.0, sin_om));
    double omega = asin(sin_om);

    if (!do_grad) return omega;

    double cos_om = cos(omega);
    double nenner = rnn * rvn * cos_om;
    if (fabs(nenner) < 1e-14) {
        for (int a = 0; a < 4; ++a) derivate[a][0] = derivate[a][1] = derivate[a][2] = 0.0;
        return omega;
    }

    double onenner = 1.0 / nenner;

    // rdme = rd - re = r_k - r_i
    double rdme[3];
    for (int d = 0; d < 3; ++d) rdme[d] = rd[d] - re[d];

    // Cross products (gfnff_helpers.f90:477-482)
    // rve = rv × re
    double rve[3];
    rve[0] = rv[1]*re[2] - rv[2]*re[1];
    rve[1] = rv[2]*re[0] - rv[0]*re[2];
    rve[2] = rv[0]*re[1] - rv[1]*re[0];

    // rne = rn × re
    double rne[3];
    rne[0] = rn[1]*re[2] - rn[2]*re[1];
    rne[1] = rn[2]*re[0] - rn[0]*re[2];
    rne[2] = rn[0]*re[1] - rn[1]*re[0];

    // rdv = rd × rv
    double rdv[3];
    rdv[0] = rd[1]*rv[2] - rd[2]*rv[1];
    rdv[1] = rd[2]*rv[0] - rd[0]*rv[2];
    rdv[2] = rd[0]*rv[1] - rd[1]*rv[0];

    // rdn = rd × rn
    double rdn[3];
    rdn[0] = rd[1]*rn[2] - rd[2]*rn[1];
    rdn[1] = rd[2]*rn[0] - rd[0]*rn[2];
    rdn[2] = rd[0]*rn[1] - rd[1]*rn[0];

    // rvdme = rv × rdme
    double rvdme[3];
    rvdme[0] = rv[1]*rdme[2] - rv[2]*rdme[1];
    rvdme[1] = rv[2]*rdme[0] - rv[0]*rdme[2];
    rvdme[2] = rv[0]*rdme[1] - rv[1]*rdme[0];

    // rndme = rn × rdme
    double rndme[3];
    rndme[0] = rn[1]*rdme[2] - rn[2]*rdme[1];
    rndme[1] = rn[2]*rdme[0] - rn[0]*rdme[2];
    rndme[2] = rn[0]*rdme[1] - rn[1]*rdme[0];

    // Gradient formulas (gfnff_helpers.f90:489-499)
    double rvn_over_rnn = rvn * inv_rnn;
    double rnn_over_rvn = rnn * inv_rvn;

    // row(0) = dω/dr_i (center)
    for (int d = 0; d < 3; ++d)
        derivate[0][d] = onenner * (rdv[d] - rn[d] - sin_om * (rvn_over_rnn * rdn[d] - rnn_over_rvn * rv[d]));

    // row(1) = dω/dr_j (nb1)
    for (int d = 0; d < 3; ++d)
        derivate[1][d] = onenner * (rvdme[d] - sin_om * rvn_over_rnn * rndme[d]);

    // row(2) = dω/dr_k (nb2)
    for (int d = 0; d < 3; ++d)
        derivate[2][d] = onenner * (rve[d] - sin_om * rvn_over_rnn * rne[d]);

    // row(3) = dω/dr_l (nb3)
    for (int d = 0; d < 3; ++d)
        derivate[3][d] = onenner * (rn[d] - sin_om * rnn_over_rvn * rv[d]);

    return omega;
}

// ============================================================================
// Kernel 7: Inversions (out-of-plane bending)
// type 0:  E = fc * (1 - cos(omega)) * damp
// type -1: E = fc * (cos(omega) - cos(omega0))^2 * damp
// Reference: ff_workspace_gfnff.cpp:calcInversions
// ============================================================================

__global__ void k_inversions(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ idx_l,
    const int*    __restrict__ ati,
    const int*    __restrict__ atj,
    const int*    __restrict__ atk,
    const int*    __restrict__ atl,
    const double* __restrict__ fc,
    const double* __restrict__ omega0,
    const double* __restrict__ C0,
    const double* __restrict__ C1,
    const double* __restrict__ C2,
    const int*    __restrict__ potential_type,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid], k = idx_k[tid], l = idx_l[tid];

        const double rcov_scale = 4.0 / 3.0;
        const double atcutt = 0.505;

        int zi = ati[tid]-1, zj = atj[tid]-1, zk = atk[tid]-1, zl = atl[tid]-1;
        double rcov_c = (zi >= 0 && zi < 87) ? d_rcov_d3[zi] * rcov_scale : 1.0 * rcov_scale;
        double rcov_1 = (zj >= 0 && zj < 87) ? d_rcov_d3[zj] * rcov_scale : 1.0 * rcov_scale;
        double rcov_2 = (zk >= 0 && zk < 87) ? d_rcov_d3[zk] * rcov_scale : 1.0 * rcov_scale;
        double rcov_3 = (zl >= 0 && zl < 87) ? d_rcov_d3[zl] * rcov_scale : 1.0 * rcov_scale;

        auto sq_dist3 = [&](int a, int b) {
            double dx = cx[b] - cx[a];
            double dy = cy[b] - cy[a];
            double dz = cz[b] - cz[a];
            return dx*dx+dy*dy+dz*dz;
        };

        double rij_sq = sq_dist3(j, i);
        double rjk_sq = sq_dist3(j, k);
        double rjl_sq = sq_dist3(j, l);

        auto calc_damp = [&](double rsq, double ra, double rb) -> double {
            double sum = ra + rb;
            double rcut_sq = atcutt * sum * sum;
            double rr = (rsq / rcut_sq);
            rr = rr * rr;
            return 1.0 / (1.0 + rr);
        };

        double damp_ij = calc_damp(rij_sq, rcov_c, rcov_1);
        double damp_jk = calc_damp(rjk_sq, rcov_2, rcov_1);
        double damp_jl = calc_damp(rjl_sq, rcov_3, rcov_1);
        double damp    = damp_ij * damp_jk * damp_jl;

        double derivate[4][3];
        double omega = inversion_angle_grad(cx, cy, cz, i, j, k, l, derivate, true);
        double cos_om  = cos(omega);
        double cos_om0 = cos(omega0[tid]);

        double E_raw, dEdcos;
        if (potential_type[tid] == 0) {
            E_raw  = fc[tid] * (1.0 - cos_om);
            dEdcos = -fc[tid];
        } else {
            double dcos = cos_om - cos_om0;
            E_raw  = fc[tid] * dcos * dcos;
            dEdcos = 2.0 * fc[tid] * dcos;
        }

        local_E = E_raw * damp;

        double dEdOmega = dEdcos * (-sin(omega)) * damp;

        int atoms[4] = {i, j, k, l};
        for (int a = 0; a < 4; ++a) {
            add_grad(grad, atoms[a],
                     dEdOmega * derivate[a][0],
                     dEdOmega * derivate[a][1],
                     dEdOmega * derivate[a][2]);
        }

        auto calc_ddamp = [&](double rsq, double ra, double rb) -> double {
            if (rsq < 1e-8) return 0.0;
            double sum = ra + rb;
            double rcut_sq = atcutt * sum * sum;
            double rr = (rsq / rcut_sq);
            rr = rr * rr;
            double opr = 1.0 + rr;
            return -4.0 * rr / (rsq * opr * opr);
        };

        double dd_ij = calc_ddamp(rij_sq, rcov_c, rcov_1);
        double dd_jk = calc_ddamp(rjk_sq, rcov_2, rcov_1);
        double dd_jl = calc_ddamp(rjl_sq, rcov_3, rcov_1);

        // Damping gradient: peripheral atom gets +fac*d, hub (j) gets -fac*d
        // Matches Fortran gfnff_engrad.F90:1382-1385 and CPU CalculateGFNFFInversionContribution
        {
            double dx=cx[i]-cx[j], dy=cy[i]-cy[j], dz=cz[i]-cz[j];
            double fac = E_raw * dd_ij * damp_jk * damp_jl;
            add_grad(grad, i,  fac*dx,  fac*dy,  fac*dz);
            add_grad(grad, j, -fac*dx, -fac*dy, -fac*dz);
        }
        {
            double dx=cx[k]-cx[j], dy=cy[k]-cy[j], dz=cz[k]-cz[j];
            double fac = E_raw * damp_ij * dd_jk * damp_jl;
            add_grad(grad, k,  fac*dx,  fac*dy,  fac*dz);
            add_grad(grad, j, -fac*dx, -fac*dy, -fac*dz);
        }
        {
            double dx=cx[l]-cx[j], dy=cy[l]-cy[j], dz=cz[l]-cz[j];
            double fac = E_raw * damp_ij * damp_jk * dd_jl;
            add_grad(grad, l,  fac*dx,  fac*dy,  fac*dz);
            add_grad(grad, j, -fac*dx, -fac*dy, -fac*dz);
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Upload covalent radii table to GPU constant memory
// Called once during FFWorkspaceHip initialization
// ============================================================================
void upload_rcov_d3(const double* rcov, int n)
{
    int count = (n < 87) ? n : 87;
    hipMemcpyToSymbol(HIP_SYMBOL(d_rcov_d3), rcov, count * sizeof(double));
}

// ============================================================================
// Upload D4 reference CN and nref tables to constant memory (once at init)
// Claude Generated (March 2026): Phase 8 — constant memory for refcn/refn
// ============================================================================
void upload_refcn_const(const double* data, int n)
{
    int count = (n < D4_MAX_ELEM * D4_MAX_REF) ? n : D4_MAX_ELEM * D4_MAX_REF;
    hipMemcpyToSymbol(HIP_SYMBOL(d_refcn_const), data, count * sizeof(double));
}

void upload_refn_const(const int* data, int n)
{
    int count = (n < D4_MAX_ELEM) ? n : D4_MAX_ELEM;
    hipMemcpyToSymbol(HIP_SYMBOL(d_refn_const), data, count * sizeof(int));
}

// ============================================================================
// Covalent radii for HB/XB vdW radii lookup (Å, element 1..86)
// Source: gfnff_par.h GFNFFParameters::covalent_radii
// ============================================================================
__constant__ double d_cov_radii[100];  // uploaded once at init

void upload_covalent_radii(const double* radii, int n)
{
    int count = (n < 100) ? n : 100;
    hipMemcpyToSymbol(HIP_SYMBOL(d_cov_radii), radii, count * sizeof(double));
}

// ============================================================================
// Device helpers: HB/XB damping functions
// Reference: ff_workspace_gfnff.cpp:914-928
// ============================================================================

__device__ __forceinline__ double d_damping_short_range(double r, double r_vdw, double scut, double alp)
{
    double ratio = scut * r_vdw / (r * r);
    double p = 1.0;
    for (int i = 0; i < (int)alp; ++i) p *= ratio;  // pow(ratio, alp) with integer alp=6
    return 1.0 / (1.0 + p);
}

__device__ __forceinline__ double d_damping_long_range(double r, double longcut, double alp)
{
    double ratio = (r * r) / longcut;
    double p = 1.0;
    for (int i = 0; i < (int)alp; ++i) p *= ratio;
    return 1.0 / (1.0 + p);
}

__device__ __forceinline__ double d_charge_scaling(double q, double st, double sf)
{
    double exp_term = exp(st * q);
    return exp_term / (exp_term + sf);
}

// ============================================================================
// Kernel 8: Triple Bond Torsions (sTorsions)
// E = -erefhalf * cos(2φ) + erefhalf
// dE/dφ = 2 * erefhalf * sin(2φ)
// Reference: ff_workspace_gfnff.cpp:calcSTorsions (lines 611-648)
// ============================================================================

__global__ void k_storsions(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ idx_l,
    const double* __restrict__ erefhalf,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int ai = idx_i[tid], aj = idx_j[tid], ak = idx_k[tid], al = idx_l[tid];

        double derivate[4][3];
        double phi = dihedral_angle_grad(cx, cy, cz, ai, aj, ak, al, derivate, true);

        double eref = erefhalf[tid];
        local_E = -eref * cos(2.0 * phi) + eref;

        double dEdphi = 2.0 * eref * sin(2.0 * phi);
        int atoms[4] = {ai, aj, ak, al};
        for (int a = 0; a < 4; ++a) {
            add_grad(grad, atoms[a],
                     dEdphi * derivate[a][0],
                     dEdphi * derivate[a][1],
                     dEdphi * derivate[a][2]);
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Kernel 9: Bonded ATM (BATM) for 1,4-pairs
// E = c9 * angr9  where c9 = fi*fj*fk * zb3atm_i*zb3atm_j*zb3atm_k
// fi = clamp(1 - 3*q_i, -4, 4)   (topology charges)
// Reference: ff_workspace_gfnff.cpp:calcBATM (lines 1652-1726)
// ============================================================================

__global__ void k_batm(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const double* __restrict__ zb3atm_i,
    const double* __restrict__ zb3atm_j,
    const double* __restrict__ zb3atm_k,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ topo_charges,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid], k = idx_k[tid];

        double ix = cx[i], iy = cy[i], iz = cz[i];
        double jx = cx[j], jy = cy[j], jz = cz[j];
        double kx = cx[k], ky = cy[k], kz = cz[k];

        double rij_x = jx-ix, rij_y = jy-iy, rij_z = jz-iz;
        double rik_x = kx-ix, rik_y = ky-iy, rik_z = kz-iz;
        double rjk_x = kx-jx, rjk_y = ky-jy, rjk_z = kz-jz;

        double r2ij = rij_x*rij_x + rij_y*rij_y + rij_z*rij_z;
        double r2ik = rik_x*rik_x + rik_y*rik_y + rik_z*rik_z;
        double r2jk = rjk_x*rjk_x + rjk_y*rjk_y + rjk_z*rjk_z;

        double rij = sqrt(r2ij), rik = sqrt(r2ik), rjk = sqrt(r2jk);
        if (rij >= 1e-10 && rik >= 1e-10 && rjk >= 1e-10) {
            double rijk3 = r2ij * r2jk * r2ik;
            double mijk  = -r2ij + r2jk + r2ik;
            double imjk  = r2ij - r2jk + r2ik;
            double ijmk  = r2ij + r2jk - r2ik;

            double ang = 0.375 * ijmk * imjk * mijk / rijk3;
            double rav3 = r2ij * rij * r2ik * rik * r2jk * rjk;

            double angr9 = (ang + 1.0) / rav3;

            const double fqq = 3.0;
            double fi = fmin(fmax(1.0 - fqq * topo_charges[i], -4.0), 4.0);
            double fj = fmin(fmax(1.0 - fqq * topo_charges[j], -4.0), 4.0);
            double fk = fmin(fmax(1.0 - fqq * topo_charges[k], -4.0), 4.0);

            double c9 = fi * fj * fk * zb3atm_i[tid] * zb3atm_j[tid] * zb3atm_k[tid];
            local_E = c9 * angr9;

            double dang_ij = -0.375 * (r2ij*r2ij*r2ij + r2ij*r2ij*(r2jk + r2ik)
                              + r2ij*(3.0*r2jk*r2jk + 2.0*r2jk*r2ik + 3.0*r2ik*r2ik)
                              - 5.0*(r2jk - r2ik)*(r2jk - r2ik)*(r2jk + r2ik))
                              / (rij * rijk3 * rav3);

            double dang_jk = -0.375 * (r2jk*r2jk*r2jk + r2jk*r2jk*(r2ik + r2ij)
                              + r2jk*(3.0*r2ik*r2ik + 2.0*r2ik*r2ij + 3.0*r2ij*r2ij)
                              - 5.0*(r2ik - r2ij)*(r2ik - r2ij)*(r2ik + r2ij))
                              / (rjk * rijk3 * rav3);

            double dang_ik = -0.375 * (r2ik*r2ik*r2ik + r2ik*r2ik*(r2jk + r2ij)
                              + r2ik*(3.0*r2jk*r2jk + 2.0*r2jk*r2ij + 3.0*r2ij*r2ij)
                              - 5.0*(r2jk - r2ij)*(r2jk - r2ij)*(r2jk + r2ij))
                              / (rik * rijk3 * rav3);

            double fij = -dang_ij * c9 / rij;
            double fjk = -dang_jk * c9 / rjk;
            double fik = -dang_ik * c9 / rik;

            add_grad(grad, j, -fij*rij_x + fjk*rjk_x, -fij*rij_y + fjk*rjk_y, -fij*rij_z + fjk*rjk_z);
            add_grad(grad, k, -fik*rik_x - fjk*rjk_x, -fik*rik_y - fjk*rjk_y, -fik*rik_z - fjk*rjk_z);
            add_grad(grad, i,  fij*rij_x + fik*rik_x,  fij*rij_y + fik*rik_y,  fij*rij_z + fik*rik_z);
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Kernel 9b: Bonded ATM with FP32 intermediates (mixed precision)
// Claude Generated (March 2026): Saves ~35 registers on distance/angle intermediates.
// Same formula as k_batm, FP64 accumulation at add_grad()/blockReduceAddEnergy().
// ============================================================================
__global__ GFNFF_KERNEL_BOUNDS void k_batm_mixed(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const double* __restrict__ zb3atm_i,
    const double* __restrict__ zb3atm_j,
    const double* __restrict__ zb3atm_k,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ topo_charges,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid], k = idx_k[tid];

        float ix = (float)cx[i], iy = (float)cy[i], iz = (float)cz[i];
        float jx = (float)cx[j], jy = (float)cy[j], jz = (float)cz[j];
        float kx = (float)cx[k], ky = (float)cy[k], kz = (float)cz[k];

        float rij_x = jx-ix, rij_y = jy-iy, rij_z = jz-iz;
        float rik_x = kx-ix, rik_y = ky-iy, rik_z = kz-iz;
        float rjk_x = kx-jx, rjk_y = ky-jy, rjk_z = kz-jz;

        float r2ij = rij_x*rij_x + rij_y*rij_y + rij_z*rij_z;
        float r2ik = rik_x*rik_x + rik_y*rik_y + rik_z*rik_z;
        float r2jk = rjk_x*rjk_x + rjk_y*rjk_y + rjk_z*rjk_z;

        float rij = sqrtf(r2ij), rik = sqrtf(r2ik), rjk = sqrtf(r2jk);
        if (rij >= 1e-10f && rik >= 1e-10f && rjk >= 1e-10f) {
            float rijk3 = r2ij * r2jk * r2ik;
            float mijk  = -r2ij + r2jk + r2ik;
            float imjk  = r2ij - r2jk + r2ik;
            float ijmk  = r2ij + r2jk - r2ik;

            float ang = 0.375f * ijmk * imjk * mijk / rijk3;
            float rav3 = r2ij * rij * r2ik * rik * r2jk * rjk;

            float angr9 = (ang + 1.0f) / rav3;

            const float fqq = 3.0f;
            float fi = fminf(fmaxf(1.0f - fqq * (float)topo_charges[i], -4.0f), 4.0f);
            float fj = fminf(fmaxf(1.0f - fqq * (float)topo_charges[j], -4.0f), 4.0f);
            float fk = fminf(fmaxf(1.0f - fqq * (float)topo_charges[k], -4.0f), 4.0f);

            float c9 = fi * fj * fk * (float)zb3atm_i[tid] * (float)zb3atm_j[tid] * (float)zb3atm_k[tid];
            local_E = (double)(c9 * angr9);

            float dang_ij = -0.375f * (r2ij*r2ij*r2ij + r2ij*r2ij*(r2jk + r2ik)
                              + r2ij*(3.0f*r2jk*r2jk + 2.0f*r2jk*r2ik + 3.0f*r2ik*r2ik)
                              - 5.0f*(r2jk - r2ik)*(r2jk - r2ik)*(r2jk + r2ik))
                              / (rij * rijk3 * rav3);

            float dang_jk = -0.375f * (r2jk*r2jk*r2jk + r2jk*r2jk*(r2ik + r2ij)
                              + r2jk*(3.0f*r2ik*r2ik + 2.0f*r2ik*r2ij + 3.0f*r2ij*r2ij)
                              - 5.0f*(r2ik - r2ij)*(r2ik - r2ij)*(r2ik + r2ij))
                              / (rjk * rijk3 * rav3);

            float dang_ik = -0.375f * (r2ik*r2ik*r2ik + r2ik*r2ik*(r2jk + r2ij)
                              + r2ik*(3.0f*r2jk*r2jk + 2.0f*r2jk*r2ij + 3.0f*r2ij*r2ij)
                              - 5.0f*(r2jk - r2ij)*(r2jk - r2ij)*(r2jk + r2ij))
                              / (rik * rijk3 * rav3);

            float fij = -dang_ij * c9 / rij;
            float fjk = -dang_jk * c9 / rjk;
            float fik = -dang_ik * c9 / rik;

            add_grad(grad, j, (double)(-fij*rij_x + fjk*rjk_x), (double)(-fij*rij_y + fjk*rjk_y), (double)(-fij*rij_z + fjk*rjk_z));
            add_grad(grad, k, (double)(-fik*rik_x - fjk*rjk_x), (double)(-fik*rik_y - fjk*rjk_y), (double)(-fik*rik_z - fjk*rjk_z));
            add_grad(grad, i, (double)( fij*rij_x + fik*rik_x), (double)( fij*rij_y + fik*rik_y), (double)( fij*rij_z + fik*rik_z));
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Kernel 10: ATM Three-Body Dispersion (Axilrod-Teller-Muto)
// E = c9 * ang * fdmp / 3.0 * triple_scale
// Reference: ff_workspace_gfnff.cpp:calcATM (1517-1564) + calcATMGradient (1575-1645)
// ============================================================================

__global__ void k_atm(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ ati,
    const int*    __restrict__ atj,
    const int*    __restrict__ atk,
    const double* __restrict__ C6_ij,
    const double* __restrict__ C6_ik,
    const double* __restrict__ C6_jk,
    const double* __restrict__ s9,
    const double* __restrict__ a1,
    const double* __restrict__ a2,
    const double* __restrict__ alp,
    const double* __restrict__ triple_scale,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid], k = idx_k[tid];

        double ix = cx[i], iy = cy[i], iz = cz[i];
        double jx = cx[j], jy = cy[j], jz = cz[j];
        double kx = cx[k], ky = cy[k], kz = cz[k];

        double rij_x = jx-ix, rij_y = jy-iy, rij_z = jz-iz;
        double rik_x = kx-ix, rik_y = ky-iy, rik_z = kz-iz;
        double rjk_x = kx-jx, rjk_y = ky-jy, rjk_z = kz-jz;

        double r2ij = rij_x*rij_x + rij_y*rij_y + rij_z*rij_z;
        double r2ik = rik_x*rik_x + rik_y*rik_y + rik_z*rik_z;
        double r2jk = rjk_x*rjk_x + rjk_y*rjk_y + rjk_z*rjk_z;

        double rij = sqrt(r2ij), rik = sqrt(r2ik), rjk = sqrt(r2jk);
        if (rij >= 1e-10 && rik >= 1e-10 && rjk >= 1e-10) {
            double c9_val = s9[tid] * sqrt(fabs(C6_ij[tid] * C6_ik[tid] * C6_jk[tid]));

            int zi = ati[tid] - 1, zj = atj[tid] - 1, zk = atk[tid] - 1;
            double r_cov_i = (zi >= 0 && zi < 87) ? d_rcov_d3[zi] : 1.0;
            double r_cov_j = (zj >= 0 && zj < 87) ? d_rcov_d3[zj] : 1.0;
            double r_cov_k = (zk >= 0 && zk < 87) ? d_rcov_d3[zk] : 1.0;

            double a1_v = a1[tid], a2_v = a2[tid], alp_v = alp[tid];
            double r0ij = a1_v * sqrt(3.0 * r_cov_i * r_cov_j) + a2_v;
            double r0ik = a1_v * sqrt(3.0 * r_cov_i * r_cov_k) + a2_v;
            double r0jk = a1_v * sqrt(3.0 * r_cov_j * r_cov_k) + a2_v;
            double r0ijk = r0ij * r0ik * r0jk;

            double rijk = rij * rik * rjk;
            double r2ijk = r2ij * r2ik * r2jk;
            double r3ijk = rijk * r2ijk;
            double r5ijk = r2ijk * r3ijk;

            double tmp = pow(r0ijk / rijk, alp_v / 3.0);
            double fdmp = 1.0 / (1.0 + 6.0 * tmp);

            double A = r2ij + r2jk - r2ik;
            double B = r2ij + r2ik - r2jk;
            double C = r2ik + r2jk - r2ij;
            double ang = (0.375 * A * B * C / r2ijk + 1.0) / r3ijk;

            double tscale = triple_scale[tid];
            local_E = ang * fdmp * c9_val / 3.0 * tscale;

            double c9_grad = -s9[tid] * sqrt(fabs(C6_ij[tid] * C6_ik[tid] * C6_jk[tid]));
            double dfdmp = -2.0 * alp_v * tmp * fdmp * fdmp;

            double dang_ij = -0.375 * (r2ij*r2ij*r2ij + r2ij*r2ij*(r2jk + r2ik)
                              + r2ij*(3.0*r2jk*r2jk + 2.0*r2jk*r2ik + 3.0*r2ik*r2ik)
                              - 5.0*(r2jk - r2ik)*(r2jk - r2ik)*(r2jk + r2ik)) / r5ijk;

            double dang_ik = -0.375 * (r2ik*r2ik*r2ik + r2ik*r2ik*(r2jk + r2ij)
                              + r2ik*(3.0*r2jk*r2jk + 2.0*r2jk*r2ij + 3.0*r2ij*r2ij)
                              - 5.0*(r2jk - r2ij)*(r2jk - r2ij)*(r2jk + r2ij)) / r5ijk;

            double dang_jk = -0.375 * (r2jk*r2jk*r2jk + r2jk*r2jk*(r2ik + r2ij)
                              + r2jk*(3.0*r2ik*r2ik + 2.0*r2ik*r2ij + 3.0*r2ij*r2ij)
                              - 5.0*(r2ik - r2ij)*(r2ik - r2ij)*(r2ik + r2ij)) / r5ijk;

            double pf = c9_grad * tscale / 3.0;
            double fij = pf * (-dang_ij * fdmp + ang * dfdmp) / r2ij;
            double fik = pf * (-dang_ik * fdmp + ang * dfdmp) / r2ik;
            double fjk = pf * (-dang_jk * fdmp + ang * dfdmp) / r2jk;

            add_grad(grad, i, -(fij*rij_x + fik*rik_x), -(fij*rij_y + fik*rik_y), -(fij*rij_z + fik*rik_z));
            add_grad(grad, j, fij*rij_x - fjk*rjk_x, fij*rij_y - fjk*rjk_y, fij*rij_z - fjk*rjk_z);
            add_grad(grad, k, fik*rik_x + fjk*rjk_x, fik*rik_y + fjk*rjk_y, fik*rik_z + fjk*rjk_z);
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Kernel 11: Halogen Bonds (three-body A-X...B)
// E = -R_damp * Q_B * acidity_X * Q_X
// R_damp = damp_short * damp_long * damp_outl / r_XB^3
// Reference: ff_workspace_gfnff.cpp:calcHalogenBonds (lines 1402-1510)
// ============================================================================

__global__ void k_xbonds(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ elem_A,
    const int*    __restrict__ elem_B,
    const double* __restrict__ q_X,
    const double* __restrict__ q_B,
    const double* __restrict__ acidity_X,
    const double* __restrict__ r_cut,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        const double XB_SCUT_D = 5.0;
        const double XB_BACUT_D = 70.0;
        const double HB_ALP_D = 6.0;
        const double HB_LONGCUT_XB_D = 70.0;
        const double XB_ST_D = 15.0;
        const double XB_SF_D = 0.03;

        int A = idx_i[tid], X = idx_j[tid], B = idx_k[tid];

        double ax = cx[A], ay = cy[A], az = cz[A];
        double xx = cx[X], xy = cy[X], xz = cz[X];
        double bx = cx[B], by = cy[B], bz = cz[B];

        double rAX_x = xx-ax, rAX_y = xy-ay, rAX_z = xz-az; applyMIC(rAX_x, rAX_y, rAX_z);
        double rXB_x = bx-xx, rXB_y = by-xy, rXB_z = bz-xz; applyMIC(rXB_x, rXB_y, rXB_z);
        double rAB_x = bx-ax, rAB_y = by-ay, rAB_z = bz-az; applyMIC(rAB_x, rAB_y, rAB_z);

        double r2AX = rAX_x*rAX_x + rAX_y*rAX_y + rAX_z*rAX_z;
        double r2XB = rXB_x*rXB_x + rXB_y*rXB_y + rXB_z*rXB_z;
        double r2AB = rAB_x*rAB_x + rAB_y*rAB_y + rAB_z*rAB_z;

        double r_AX = sqrt(r2AX), r_XB = sqrt(r2XB), r_AB = sqrt(r2AB);

        if (r_XB <= r_cut[tid] && r_XB >= 1e-10 && r_AB >= 1e-10) {
            int eA = elem_A[tid] - 1, eB = elem_B[tid] - 1;
            double r_vdw_AB = ((eA >= 0 && eA < 100) ? d_cov_radii[eA] : 1.0)
                            + ((eB >= 0 && eB < 100) ? d_cov_radii[eB] : 1.0);

            double damp_short = d_damping_short_range(r_XB, r_vdw_AB, XB_SCUT_D, HB_ALP_D);
            double damp_long = d_damping_long_range(r_XB, HB_LONGCUT_XB_D, HB_ALP_D);

            double ratio_outl = (r_AX + r_XB) / r_AB;
            double expo_outl = XB_BACUT_D * (ratio_outl - 1.0);

            if (expo_outl <= 15.0) {
                double damp_outl = 2.0 / (1.0 + exp(expo_outl));

                double Q_X = d_charge_scaling(q_X[tid], XB_ST_D, XB_SF_D);
                double Q_B = d_charge_scaling(-q_B[tid], XB_ST_D, XB_SF_D);

                double r3XB = r_XB * r_XB * r_XB;
                double R_damp = damp_short * damp_long * damp_outl / r3XB;
                local_E = -R_damp * Q_B * acidity_X[tid] * Q_X;

                // Gradient
                double ratio_short = XB_SCUT_D * r_vdw_AB / (r_XB * r_XB);
                double damp_short_term = 1.0;
                for (int p = 0; p < (int)HB_ALP_D; ++p) damp_short_term *= ratio_short;
                double ddamp_short_dr = -2.0 * HB_ALP_D * damp_short * damp_short_term
                                      / (r_XB * (1.0 + damp_short_term));

                double ratio_long = (r_XB * r_XB) / HB_LONGCUT_XB_D;
                double damp_long_term = 1.0;
                for (int p = 0; p < (int)HB_ALP_D; ++p) damp_long_term *= ratio_long;
                double ddamp_long_dr = -2.0 * HB_ALP_D * r_XB * damp_long * damp_long_term
                                     / (HB_LONGCUT_XB_D * (1.0 + damp_long_term));

                double exp_term = exp(expo_outl);
                double denom_outl = 1.0 + exp_term;
                double ddamp_outl_drAX = -2.0 * exp_term * XB_BACUT_D / (r_AB * denom_outl * denom_outl);
                double ddamp_outl_drXB = ddamp_outl_drAX;
                double ddamp_outl_drAB = 2.0 * exp_term * XB_BACUT_D * (r_AX + r_XB)
                                       / (r_AB * r_AB * denom_outl * denom_outl);

                double dRdamp_drXB = (ddamp_short_dr * damp_long * damp_outl
                                    + damp_short * ddamp_long_dr * damp_outl
                                    + damp_short * damp_long * ddamp_outl_drXB) / r3XB
                                   - 3.0 * R_damp / r_XB;
                double dRdamp_drAX = damp_short * damp_long * ddamp_outl_drAX / r3XB;
                double dRdamp_drAB = damp_short * damp_long * ddamp_outl_drAB / r3XB;

                double E_pf = -Q_B * acidity_X[tid] * Q_X;
                double dE_drXB = E_pf * dRdamp_drXB;
                double dE_drAX = E_pf * dRdamp_drAX;
                double dE_drAB = E_pf * dRdamp_drAB;

                double inv_rAX = 1.0 / (r_AX + 1e-14);
                double inv_rXB = 1.0 / (r_XB + 1e-14);
                double inv_rAB = 1.0 / (r_AB + 1e-14);

                add_grad(grad, A, -dE_drAX * rAX_x * inv_rAX - dE_drAB * rAB_x * inv_rAB,
                                  -dE_drAX * rAX_y * inv_rAX - dE_drAB * rAB_y * inv_rAB,
                                  -dE_drAX * rAX_z * inv_rAX - dE_drAB * rAB_z * inv_rAB);
                add_grad(grad, X, dE_drAX * rAX_x * inv_rAX - dE_drXB * rXB_x * inv_rXB,
                                  dE_drAX * rAX_y * inv_rAX - dE_drXB * rXB_y * inv_rXB,
                                  dE_drAX * rAX_z * inv_rAX - dE_drXB * rXB_z * inv_rXB);
                add_grad(grad, B, dE_drXB * rXB_x * inv_rXB + dE_drAB * rAB_x * inv_rAB,
                                  dE_drXB * rXB_y * inv_rXB + dE_drAB * rAB_y * inv_rAB,
                                  dE_drXB * rXB_z * inv_rXB + dE_drAB * rAB_z * inv_rAB);
            }
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Kernel 11b: Halogen Bonds with FP32 intermediates (mixed precision)
// Claude Generated (March 2026): FP32 for distance/damping intermediates.
// Same formula as k_xbonds, FP64 at add_grad()/blockReduceAddEnergy().
// ============================================================================
__global__ GFNFF_KERNEL_BOUNDS void k_xbonds_mixed(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ elem_A,
    const int*    __restrict__ elem_B,
    const double* __restrict__ q_X,
    const double* __restrict__ q_B,
    const double* __restrict__ acidity_X,
    const double* __restrict__ r_cut,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        const float XB_SCUT_D = 5.0f;
        const float XB_BACUT_D = 70.0f;
        const float HB_ALP_D = 6.0f;
        const float HB_LONGCUT_XB_D = 70.0f;
        const float XB_ST_D = 15.0f;
        const float XB_SF_D = 0.03f;

        int A = idx_i[tid], X = idx_j[tid], B = idx_k[tid];

        float ax = (float)cx[A], ay = (float)cy[A], az = (float)cz[A];
        float xx = (float)cx[X], xy = (float)cy[X], xz = (float)cz[X];
        float bx = (float)cx[B], by = (float)cy[B], bz = (float)cz[B];

        float rAX_x = xx-ax, rAX_y = xy-ay, rAX_z = xz-az;
        float rXB_x = bx-xx, rXB_y = by-xy, rXB_z = bz-xz;
        float rAB_x = bx-ax, rAB_y = by-ay, rAB_z = bz-az;

        float r2AX = rAX_x*rAX_x + rAX_y*rAX_y + rAX_z*rAX_z;
        float r2XB = rXB_x*rXB_x + rXB_y*rXB_y + rXB_z*rXB_z;
        float r2AB = rAB_x*rAB_x + rAB_y*rAB_y + rAB_z*rAB_z;

        float r_AX = sqrtf(r2AX), r_XB = sqrtf(r2XB), r_AB = sqrtf(r2AB);

        if (r_XB <= (float)r_cut[tid] && r_XB >= 1e-10f && r_AB >= 1e-10f) {
            int eA = elem_A[tid] - 1, eB = elem_B[tid] - 1;
            float r_vdw_AB = ((eA >= 0 && eA < 100) ? (float)d_cov_radii[eA] : 1.0f)
                           + ((eB >= 0 && eB < 100) ? (float)d_cov_radii[eB] : 1.0f);

            // Short-range damping (integer-exponent loop, exact in FP32)
            float ratio_s = XB_SCUT_D * r_vdw_AB / (r_XB * r_XB);
            float p_s = 1.0f;
            for (int p = 0; p < 6; ++p) p_s *= ratio_s;
            float damp_short = 1.0f / (1.0f + p_s);

            // Long-range damping
            float ratio_l = (r_XB * r_XB) / HB_LONGCUT_XB_D;
            float p_l = 1.0f;
            for (int p = 0; p < 6; ++p) p_l *= ratio_l;
            float damp_long = 1.0f / (1.0f + p_l);

            float ratio_outl = (r_AX + r_XB) / r_AB;
            float expo_outl = XB_BACUT_D * (ratio_outl - 1.0f);

            if (expo_outl <= 15.0f) {
                float damp_outl = 2.0f / (1.0f + expf(expo_outl));

                // Charge scaling (using double for exp precision, cast result)
                float Q_X_f = (float)d_charge_scaling(q_X[tid], (double)XB_ST_D, (double)XB_SF_D);
                float Q_B_f = (float)d_charge_scaling(-q_B[tid], (double)XB_ST_D, (double)XB_SF_D);

                float r3XB = r_XB * r_XB * r_XB;
                float R_damp = damp_short * damp_long * damp_outl / r3XB;
                local_E = (double)(-R_damp * Q_B_f * (float)acidity_X[tid] * Q_X_f);

                // Gradient
                float damp_short_term = 1.0f;
                float ratio_short_g = XB_SCUT_D * r_vdw_AB / (r_XB * r_XB);
                for (int p = 0; p < 6; ++p) damp_short_term *= ratio_short_g;
                float ddamp_short_dr = -2.0f * HB_ALP_D * damp_short * damp_short_term
                                      / (r_XB * (1.0f + damp_short_term));

                float damp_long_term = 1.0f;
                float ratio_long_g = (r_XB * r_XB) / HB_LONGCUT_XB_D;
                for (int p = 0; p < 6; ++p) damp_long_term *= ratio_long_g;
                float ddamp_long_dr = -2.0f * HB_ALP_D * r_XB * damp_long * damp_long_term
                                     / (HB_LONGCUT_XB_D * (1.0f + damp_long_term));

                float exp_t = expf(expo_outl);
                float denom_outl = 1.0f + exp_t;
                float ddamp_outl_drAX = -2.0f * exp_t * XB_BACUT_D / (r_AB * denom_outl * denom_outl);
                float ddamp_outl_drXB = ddamp_outl_drAX;
                float ddamp_outl_drAB = 2.0f * exp_t * XB_BACUT_D * (r_AX + r_XB)
                                       / (r_AB * r_AB * denom_outl * denom_outl);

                float dRdamp_drXB = (ddamp_short_dr * damp_long * damp_outl
                                    + damp_short * ddamp_long_dr * damp_outl
                                    + damp_short * damp_long * ddamp_outl_drXB) / r3XB
                                   - 3.0f * R_damp / r_XB;
                float dRdamp_drAX = damp_short * damp_long * ddamp_outl_drAX / r3XB;
                float dRdamp_drAB = damp_short * damp_long * ddamp_outl_drAB / r3XB;

                float E_pf = -Q_B_f * (float)acidity_X[tid] * Q_X_f;
                float dE_drXB = E_pf * dRdamp_drXB;
                float dE_drAX = E_pf * dRdamp_drAX;
                float dE_drAB = E_pf * dRdamp_drAB;

                float inv_rAX = 1.0f / (r_AX + 1e-14f);
                float inv_rXB = 1.0f / (r_XB + 1e-14f);
                float inv_rAB = 1.0f / (r_AB + 1e-14f);

                add_grad(grad, A, (double)(-dE_drAX * rAX_x * inv_rAX - dE_drAB * rAB_x * inv_rAB),
                                  (double)(-dE_drAX * rAX_y * inv_rAX - dE_drAB * rAB_y * inv_rAB),
                                  (double)(-dE_drAX * rAX_z * inv_rAX - dE_drAB * rAB_z * inv_rAB));
                add_grad(grad, X, (double)( dE_drAX * rAX_x * inv_rAX - dE_drXB * rXB_x * inv_rXB),
                                  (double)( dE_drAX * rAX_y * inv_rAX - dE_drXB * rXB_y * inv_rXB),
                                  (double)( dE_drAX * rAX_z * inv_rAX - dE_drXB * rXB_z * inv_rXB));
                add_grad(grad, B, (double)( dE_drXB * rXB_x * inv_rXB + dE_drAB * rAB_x * inv_rAB),
                                  (double)( dE_drXB * rXB_y * inv_rXB + dE_drAB * rAB_y * inv_rAB),
                                  (double)( dE_drXB * rXB_z * inv_rXB + dE_drAB * rAB_z * inv_rAB));
            }
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Kernel 12: Hydrogen Bonds (three-body A-H...B)
// Most complex kernel: 4 case types, neighbor damping, virtual LP
// Reference: ff_workspace_gfnff.cpp:calcHydrogenBonds (lines 938-1395)
//
// ENERGY-ONLY for now — gradient is computed on CPU via postProcessCPU
// or a dedicated gradient kernel will follow once energy is validated.
//
// Simplification for GPU: Case 3 (eangl*etors) and Case 4 (virtual LP)
// require variable-length neighbor iteration which maps poorly to SIMT.
// This kernel computes all 4 cases including neighbor effects.
// ============================================================================

// Claude Generated (June 2026): device-resident HB charge refresh (proto in .cuh).
// One thread per HB triplet: copy the live per-atom EEQ charges into the HB SoA q
// arrays by donor(A=idx_i)/H(idx_j)/acceptor(B=idx_k) index.
__global__ void k_gather_hb_charges(
    int n,
    const int* __restrict__ idx_i, const int* __restrict__ idx_j, const int* __restrict__ idx_k,
    const double* __restrict__ q_atom,
    double* __restrict__ q_A, double* __restrict__ q_H, double* __restrict__ q_B)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < n) {
        q_A[tid] = q_atom[idx_i[tid]]; // donor A
        q_H[tid] = q_atom[idx_j[tid]]; // hydrogen H
        q_B[tid] = q_atom[idx_k[tid]]; // acceptor B
    }
}

__global__ void k_hbonds(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ elem_A,
    const int*    __restrict__ elem_B,
    const double* __restrict__ q_H,
    const double* __restrict__ q_A,
    const double* __restrict__ q_B,
    const double* __restrict__ basicity_A,
    const double* __restrict__ basicity_B,
    const double* __restrict__ acidity_A,
    const double* __restrict__ acidity_B,
    const double* __restrict__ r_cut,
    const int*    __restrict__ case_type,
    const int*    __restrict__ nb_B_offset,
    const int*    __restrict__ nb_B_count,
    const int*    __restrict__ nb_B_flat,
    const int*    __restrict__ acceptor_parent,
    const int*    __restrict__ nb_C_offset,
    const int*    __restrict__ nb_C_count,
    const int*    __restrict__ nb_C_flat,
    const double* __restrict__ repz_B,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {

    // Constants
    const double HB_SCUT_D = 22.0;
    const double HB_ALP_D = 6.0;
    const double HB_LONGCUT_D = 85.0;
    const double HB_BACUT_D = 49.0;
    const double HB_ST_D = 15.0;
    const double HB_SF_D = 1.0;
    const double HB_NBCUT_D = 11.20;
    const double XHACI_GLOBABH_D = 0.268;
    const double XHACI_COH_D = 0.350;
    const double BEND_HB_D = 0.20;
    const double TORS_HB_D = 0.94;
    const double HBLPCUT_D = 56.0;

    int A = idx_i[tid], H = idx_j[tid], B = idx_k[tid];
    int ct = case_type[tid];

    double ax = cx[A], ay = cy[A], az = cz[A];
    double hx = cx[H], hy = cy[H], hz = cz[H];
    double bx = cx[B], by = cy[B], bz = cz[B];

    double rAH_x = hx-ax, rAH_y = hy-ay, rAH_z = hz-az; applyMIC(rAH_x, rAH_y, rAH_z);
    double rHB_x = bx-hx, rHB_y = by-hy, rHB_z = bz-hz; applyMIC(rHB_x, rHB_y, rHB_z);
    double rAB_x = bx-ax, rAB_y = by-ay, rAB_z = bz-az; applyMIC(rAB_x, rAB_y, rAB_z);

    double r2AH = rAH_x*rAH_x + rAH_y*rAH_y + rAH_z*rAH_z;
    double r2HB = rHB_x*rHB_x + rHB_y*rHB_y + rHB_z*rHB_z;
    double r2AB = rAB_x*rAB_x + rAB_y*rAB_y + rAB_z*rAB_z;

    double r_AH = sqrt(r2AH), r_HB = sqrt(r2HB), r_AB = sqrt(r2AB);

    if (r_AB <= r_cut[tid] && r_AB >= 1e-10) {

    double r_AH_4 = r2AH * r2AH;
    double r_HB_4 = r2HB * r2HB;
    double denom_DA = 1.0 / (r_AH_4 + r_HB_4);

    double Q_H = d_charge_scaling(q_H[tid], HB_ST_D, HB_SF_D);
    double Q_A = d_charge_scaling(-q_A[tid], HB_ST_D, HB_SF_D);
    double Q_B = d_charge_scaling(-q_B[tid], HB_ST_D, HB_SF_D);

    double bas = (Q_A * basicity_A[tid] * r_AH_4 + Q_B * basicity_B[tid] * r_HB_4) * denom_DA;
    double aci = (acidity_B[tid] * r_AH_4 + acidity_A[tid] * r_HB_4) * denom_DA;

    int eA = elem_A[tid] - 1, eB = elem_B[tid] - 1;
    double r_vdw_AB = ((eA >= 0 && eA < 100) ? d_cov_radii[eA] : 1.0)
                    + ((eB >= 0 && eB < 100) ? d_cov_radii[eB] : 1.0);

    double damp_short = d_damping_short_range(r_AB, r_vdw_AB, HB_SCUT_D, HB_ALP_D);
    double damp_long = d_damping_long_range(r_AB, HB_LONGCUT_D, HB_ALP_D);
    double damp_env = damp_short * damp_long;

    // Out-of-line damping
    double rahprbh = r_AH + r_HB + 1e-12;
    double expo = (HB_BACUT_D / r_vdw_AB) * (rahprbh / r_AB - 1.0);

    if (expo <= 15.0) {
    double damp_outl = 2.0 / (1.0 + exp(expo));

    // Neighbor out-of-line (Case >= 2)
    double outl_nb_tot = 1.0;
    if (ct >= 2) {
        int nb_off = nb_B_offset[tid];
        int nb_cnt = nb_B_count[tid];
        double hbnbcut_save = (eB + 1 == 7 && nb_cnt == 1) ? 2.0 : HB_NBCUT_D;
        for (int ni = 0; ni < nb_cnt; ++ni) {
            int nb = nb_B_flat[nb_off + ni];
            double nbx = cx[nb], nby = cy[nb], nbz = cz[nb];
            double dAnb_x = nbx-ax, dAnb_y = nby-ay, dAnb_z = nbz-az;
            double dBnb_x = nbx-bx, dBnb_y = nby-by, dBnb_z = nbz-bz;
            double r_Anb = sqrt(dAnb_x*dAnb_x + dAnb_y*dAnb_y + dAnb_z*dAnb_z);
            double r_Bnb = sqrt(dBnb_x*dBnb_x + dBnb_y*dBnb_y + dBnb_z*dBnb_z);
            double expo_nb = (hbnbcut_save / r_vdw_AB) * ((r_Anb + r_Bnb) / r_AB - 1.0);
            outl_nb_tot *= (2.0 / (1.0 + exp(-expo_nb)) - 1.0);
        }
    }

    // Case 4: Virtual LP — save state for gradient
    double outl_lp = 1.0;
    double lp_dist_c4 = 0, lp_vx_c4 = 0, lp_vy_c4 = 0, lp_vz_c4 = 0;
    double vnorm_c4 = 0, ralp_c4 = 0, expo_lp_c4 = 0;
    double lpx_c4 = 0, lpy_c4 = 0, lpz_c4 = 0;
    int nbb_lp_c4 = 0;
    if (ct == 4) {
        lp_dist_c4 = 0.50 - 0.018 * repz_B[tid];
        int nb_off = nb_B_offset[tid];
        int nb_cnt = nb_B_count[tid];
        nbb_lp_c4 = nb_cnt;
        if (nb_cnt > 0) {
            for (int ni = 0; ni < nb_cnt; ++ni) {
                int nb = nb_B_flat[nb_off + ni];
                lp_vx_c4 += cx[nb] - bx;
                lp_vy_c4 += cy[nb] - by;
                lp_vz_c4 += cz[nb] - bz;
            }
            vnorm_c4 = sqrt(lp_vx_c4*lp_vx_c4 + lp_vy_c4*lp_vy_c4 + lp_vz_c4*lp_vz_c4);
            if (vnorm_c4 > 1e-10) {
                lpx_c4 = bx + (-lp_dist_c4) * lp_vx_c4 / vnorm_c4;
                lpy_c4 = by + (-lp_dist_c4) * lp_vy_c4 / vnorm_c4;
                lpz_c4 = bz + (-lp_dist_c4) * lp_vz_c4 / vnorm_c4;
                double dAlp_x = ax - lpx_c4, dAlp_y = ay - lpy_c4, dAlp_z = az - lpz_c4;
                ralp_c4 = sqrt(dAlp_x*dAlp_x + dAlp_y*dAlp_y + dAlp_z*dAlp_z);
                expo_lp_c4 = (HBLPCUT_D / r_vdw_AB) * ((ralp_c4 + lp_dist_c4 + 1e-12) / r_AB - 1.0);
                outl_lp = 2.0 / (1.0 + exp(expo_lp_c4));
            } else {
                nbb_lp_c4 = 0;
            }
        }
    }

    double qhoutl = Q_H * damp_outl * outl_nb_tot * outl_lp;

    // rdamp: case-dependent distance decay
    double rdamp;
    if (ct >= 2) {
        rdamp = damp_env * (1.8 / (r_HB * r_HB * r_HB) - 0.8 / (r_AB * r_AB * r_AB));
    } else {
        rdamp = damp_env / (r_AB * r_AB * r_AB);
    }

    // Case 3: eangl + etors
    double eangl = 1.0, etors = 1.0;
    if (ct == 3 && acceptor_parent[tid] >= 0) {
        int C_idx = acceptor_parent[tid];
        int B_idx = B, H_idx = H;

        // eangl: angle bending H...B=C
        {
            double C_x = cx[C_idx], C_y = cy[C_idx], C_z = cz[C_idx];
            double vab_x = C_x-bx, vab_y = C_y-by, vab_z = C_z-bz;
            double vcb_x = hx-bx, vcb_y = hy-by, vcb_z = hz-bz;
            double rab2a = vab_x*vab_x + vab_y*vab_y + vab_z*vab_z;
            double rcb2a = vcb_x*vcb_x + vcb_y*vcb_y + vcb_z*vcb_z;

            double cosa = (vab_x*vcb_x + vab_y*vcb_y + vab_z*vcb_z) / (sqrt(rab2a) * sqrt(rcb2a) + 1e-14);
            cosa = fmax(-1.0, fmin(1.0, cosa));
            double theta_a = acos(cosa);

            double c0 = 120.0 * M_PI / 180.0;
            double fc_bend = 1.0 - BEND_HB_D;
            double kijk = fc_bend / ((cos(0.0) - cos(c0)) * (cos(0.0) - cos(c0)));

            double ea;
            if (M_PI - c0 < 1e-6) {
                double dt = theta_a - c0;
                ea = kijk * dt * dt;
            } else {
                ea = kijk * (cosa - cos(c0)) * (cosa - cos(c0));
            }
            eangl = 1.0 - ea;
        }

        // etors: product of D-B-C-H torsion terms
        int nc_off = nb_C_offset[tid];
        int nc_cnt = nb_C_count[tid];
        for (int ni = 0; ni < nc_cnt; ++ni) {
            int D_idx = nb_C_flat[nc_off + ni];
            if (D_idx == B_idx) continue;

            double derivate[4][3];
            double phi = dihedral_angle_grad(cx, cy, cz, D_idx, B_idx, C_idx, H_idx, derivate, false);

            double tshift = TORS_HB_D;
            double fc_tors = (1.0 - tshift) / 2.0;
            double phi0_tors = M_PI / 2.0;
            int rn = 2;
            double c1 = rn * (phi - phi0_tors) + M_PI;
            double et = (1.0 + cos(c1)) * fc_tors + tshift;
            etors *= et;
        }
    }

    // Energy: case-specific formula
    double global_scale = 1.0;
    if (ct == 2 || ct == 4) global_scale = XHACI_GLOBABH_D;
    else if (ct == 3) global_scale = XHACI_COH_D;

    double E_HB;
    if (ct >= 2) {
        double const_val = acidity_A[tid] * basicity_B[tid] * Q_A * Q_B * global_scale;
        E_HB = -rdamp * qhoutl * const_val * eangl * etors;
    } else {
        E_HB = -bas * aci * rdamp * qhoutl;
    }
    local_E = E_HB;

    // ========== Gradient ==========
    // Full analytical gradient matching CPU calcHydrogenBonds
    double rab2 = r_AB * r_AB;
    double rbh2 = r_HB * r_HB;
    double rah2 = r_AH * r_AH;

    // Damping derivative intermediates
    double ratio1_pow = 1.0;
    { double ratio1 = rab2 / HB_LONGCUT_D; for (int p=0; p<(int)HB_ALP_D; ++p) ratio1_pow *= ratio1; }
    double ratio3_pow = 1.0;
    { double shortcut = HB_SCUT_D * r_vdw_AB; double ratio3 = shortcut / rab2; for (int p=0; p<(int)HB_ALP_D; ++p) ratio3_pow *= ratio3; }
    double ddamp = (-2.0 * HB_ALP_D * ratio1_pow / (1.0 + ratio1_pow))
                 + ( 2.0 * HB_ALP_D * ratio3_pow / (1.0 + ratio3_pow));

    double ratio2 = exp(expo);

    // Distance vectors (Fortran convention: A-B, A-H, B-H)
    double drab_x = ax-bx, drab_y = ay-by, drab_z = az-bz;
    double drah_x = ax-hx, drah_y = ay-hy, drah_z = az-hz;
    double drbh_x = bx-hx, drbh_y = by-hy, drbh_z = bz-hz;

    double ga_x=0, ga_y=0, ga_z=0;
    double gb_x=0, gb_y=0, gb_z=0;
    double gh_x=0, gh_y=0, gh_z=0;

    if (ct >= 2) {
        // Case 2/3/4 gradient
        double const_val = acidity_A[tid] * basicity_B[tid] * Q_A * Q_B * global_scale;
        double dterm  = -qhoutl * eangl * etors * const_val;
        double aterm  = -rdamp * Q_H * outl_nb_tot * outl_lp * eangl * etors * const_val;
        double nbterm = -rdamp * Q_H * damp_outl * outl_lp * eangl * etors * const_val;

        double p_bh = 1.8, p_ab = -0.8;
        double rbhdamp = damp_env * p_bh / (rbh2 * r_HB);
        double rabdamp = damp_env * p_ab / (rab2 * r_AB);

        // rab damping
        double gi = ((rabdamp + rbhdamp) * ddamp - 3.0 * rabdamp) / rab2 * dterm;
        ga_x = gi * drab_x; ga_y = gi * drab_y; ga_z = gi * drab_z;
        gb_x = -ga_x; gb_y = -ga_y; gb_z = -ga_z;

        // rbh damping
        gi = -3.0 * rbhdamp / rbh2 * dterm;
        gb_x += gi * drbh_x; gb_y += gi * drbh_y; gb_z += gi * drbh_z;
        gh_x = -gi * drbh_x; gh_y = -gi * drbh_y; gh_z = -gi * drbh_z;

        // Out-of-line: rab
        double tmp1 = -2.0 * aterm * ratio2 * expo
                    / ((1.0 + ratio2) * (1.0 + ratio2))
                    / (rahprbh - r_AB);
        gi = -tmp1 * rahprbh / rab2;
        ga_x += gi * drab_x; ga_y += gi * drab_y; ga_z += gi * drab_z;
        gb_x -= gi * drab_x; gb_y -= gi * drab_y; gb_z -= gi * drab_z;

        // Out-of-line: rah, rbh
        gi = tmp1 / r_AH;
        double dga_outl_x = gi * drah_x, dga_outl_y = gi * drah_y, dga_outl_z = gi * drah_z;
        ga_x += dga_outl_x; ga_y += dga_outl_y; ga_z += dga_outl_z;
        gi = tmp1 / r_HB;
        double dgb_outl_x = gi * drbh_x, dgb_outl_y = gi * drbh_y, dgb_outl_z = gi * drbh_z;
        gb_x += dgb_outl_x; gb_y += dgb_outl_y; gb_z += dgb_outl_z;
        gh_x += -dga_outl_x - dgb_outl_x;
        gh_y += -dga_outl_y - dgb_outl_y;
        gh_z += -dga_outl_z - dgb_outl_z;

        // Neighbor gradient (Case >= 2)
        {
            int nb_off = nb_B_offset[tid];
            int nb_cnt = nb_B_count[tid];
            double hbnbcut_g = (eB + 1 == 7 && nb_cnt == 1) ? 2.0 : HB_NBCUT_D;
            for (int ni = 0; ni < nb_cnt; ++ni) {
                int nb = nb_B_flat[nb_off + ni];
                double nbx = cx[nb], nby = cy[nb], nbz = cz[nb];
                double dranb_x = ax-nbx, dranb_y = ay-nby, dranb_z = az-nbz;
                double drbnb_x = bx-nbx, drbnb_y = by-nby, drbnb_z = bz-nbz;
                double ranb = sqrt(dranb_x*dranb_x + dranb_y*dranb_y + dranb_z*dranb_z);
                double rbnb = sqrt(drbnb_x*drbnb_x + drbnb_y*drbnb_y + drbnb_z*drbnb_z);
                double ranbprbnb = ranb + rbnb + 1e-12;

                double expo_nb_i = (hbnbcut_g / r_vdw_AB) * (ranbprbnb / r_AB - 1.0);
                double ratio2_nb_i = exp(-expo_nb_i);
                double outl_nb_i = 2.0 / (1.0 + ratio2_nb_i) - 1.0;

                double outl_nb_others = (fabs(outl_nb_i) > 1e-12) ? outl_nb_tot / outl_nb_i : 1.0;

                double tmp2 = 2.0 * nbterm * outl_nb_others * ratio2_nb_i * expo_nb_i
                            / ((1.0 + ratio2_nb_i) * (1.0 + ratio2_nb_i))
                            / (ranbprbnb - r_AB);

                double gi_nb = -tmp2 * ranbprbnb / rab2;
                ga_x += gi_nb * drab_x; ga_y += gi_nb * drab_y; ga_z += gi_nb * drab_z;
                gb_x -= gi_nb * drab_x; gb_y -= gi_nb * drab_y; gb_z -= gi_nb * drab_z;

                gi_nb = tmp2 / ranb;
                double dga_nb_x = gi_nb * dranb_x, dga_nb_y = gi_nb * dranb_y, dga_nb_z = gi_nb * dranb_z;
                ga_x += dga_nb_x; ga_y += dga_nb_y; ga_z += dga_nb_z;
                gi_nb = tmp2 / rbnb;
                double dgb_nb_x = gi_nb * drbnb_x, dgb_nb_y = gi_nb * drbnb_y, dgb_nb_z = gi_nb * drbnb_z;
                gb_x += dgb_nb_x; gb_y += dgb_nb_y; gb_z += dgb_nb_z;
                add_grad(grad, nb, -dga_nb_x - dgb_nb_x, -dga_nb_y - dgb_nb_y, -dga_nb_z - dgb_nb_z);
            }
        }

        // Case 3: eangl + etors gradient (simplified — gradient contributions for angle/torsion
        // are complex and involve 4+ atoms per torsion; computed here for the 3 primary atoms only)
        // Full case 3 gradient for angle bending term
        if (ct == 3 && acceptor_parent[tid] >= 0) {
            int C_idx = acceptor_parent[tid];
            double bterm_c3 = -rdamp * qhoutl * etors * acidity_A[tid] * basicity_B[tid] * Q_A * Q_B * global_scale;
            // Angle gradient: H...B=C angle
            double C_x = cx[C_idx], C_y = cy[C_idx], C_z = cz[C_idx];
            double vab_x = C_x-bx, vab_y = C_y-by, vab_z = C_z-bz;
            double vcb_x = hx-bx, vcb_y = hy-by, vcb_z = hz-bz;
            double rab2a = vab_x*vab_x + vab_y*vab_y + vab_z*vab_z;
            double rcb2a = vcb_x*vcb_x + vcb_y*vcb_y + vcb_z*vcb_z;

            // Cross product vcb × vab
            double vp_x = vcb_y*vab_z - vcb_z*vab_y;
            double vp_y = vcb_z*vab_x - vcb_x*vab_z;
            double vp_z = vcb_x*vab_y - vcb_y*vab_x;
            double rp = sqrt(vp_x*vp_x + vp_y*vp_y + vp_z*vp_z) + 1e-14;

            double cosa = (vab_x*vcb_x + vab_y*vcb_y + vab_z*vcb_z) / (sqrt(rab2a) * sqrt(rcb2a) + 1e-14);
            cosa = fmax(-1.0, fmin(1.0, cosa));
            double theta_a = acos(cosa);

            double c0 = 120.0 * M_PI / 180.0;
            double fc_bend = 1.0 - BEND_HB_D;
            double kijk = fc_bend / ((cos(0.0) - cos(c0)) * (cos(0.0) - cos(c0)));

            double deddt;
            if (M_PI - c0 < 1e-6) {
                deddt = 2.0 * kijk * (theta_a - c0);
            } else {
                deddt = 2.0 * kijk * sin(theta_a) * (cos(c0) - cosa);
            }

            // deda = vab × vp * (-deddt / (rab2a * rp))
            double ax_x = vab_y*vp_z - vab_z*vp_y;
            double ax_y = vab_z*vp_x - vab_x*vp_z;
            double ax_z = vab_x*vp_y - vab_y*vp_x;
            double f_a = -deddt / (rab2a * rp);
            double deda_x = ax_x * f_a, deda_y = ax_y * f_a, deda_z = ax_z * f_a;

            // dedc = vcb × vp * (deddt / (rcb2a * rp))
            double cx_x = vcb_y*vp_z - vcb_z*vp_y;
            double cx_y = vcb_z*vp_x - vcb_x*vp_z;
            double cx_z = vcb_x*vp_y - vcb_y*vp_x;
            double f_c = deddt / (rcb2a * rp);
            double dedc_x = cx_x * f_c, dedc_y = cx_y * f_c, dedc_z = cx_z * f_c;

            double dedb_x = deda_x + dedc_x;
            double dedb_y = deda_y + dedc_y;
            double dedb_z = deda_z + dedc_z;

            // gangl_3body: row0=B, row1=C, row2=H
            gb_x += bterm_c3 * dedb_x; gb_y += bterm_c3 * dedb_y; gb_z += bterm_c3 * dedb_z;
            add_grad(grad, C_idx, bterm_c3 * (-deda_x), bterm_c3 * (-deda_y), bterm_c3 * (-deda_z));
            gh_x += bterm_c3 * (-dedc_x); gh_y += bterm_c3 * (-dedc_y); gh_z += bterm_c3 * (-dedc_z);

            // Torsion gradient: D-B-C-H for each neighbor of C
            double tterm_c3 = -rdamp * qhoutl * eangl * acidity_A[tid] * basicity_B[tid] * Q_A * Q_B * global_scale;
            int nc_off = nb_C_offset[tid];
            int nc_cnt = nb_C_count[tid];

            // Compute individual torsion energies for product rule
            // We need etors/et_k for each k, so we need individual et values
            for (int ni = 0; ni < nc_cnt; ++ni) {
                int D_idx = nb_C_flat[nc_off + ni];
                if (D_idx == B) continue;

                double derivate[4][3];
                double phi = dihedral_angle_grad(cx, cy, cz, D_idx, B, C_idx, H, derivate, true);

                double tshift = TORS_HB_D;
                double fc_tors = (1.0 - tshift) / 2.0;
                double phi0_tors = M_PI / 2.0;
                int rn = 2;
                double c1 = rn * (phi - phi0_tors) + M_PI;
                double et = (1.0 + cos(c1)) * fc_tors + tshift;
                double dij = -rn * sin(c1) * fc_tors;

                double factor_k = (fabs(et) > 1e-12) ? etors / et : 0.0;
                double t_k = factor_k * tterm_c3;

                add_grad(grad, D_idx, t_k * dij * derivate[0][0], t_k * dij * derivate[0][1], t_k * dij * derivate[0][2]);
                gb_x += t_k * dij * derivate[1][0]; gb_y += t_k * dij * derivate[1][1]; gb_z += t_k * dij * derivate[1][2];
                add_grad(grad, C_idx, t_k * dij * derivate[2][0], t_k * dij * derivate[2][1], t_k * dij * derivate[2][2]);
                gh_x += t_k * dij * derivate[3][0]; gh_y += t_k * dij * derivate[3][1]; gh_z += t_k * dij * derivate[3][2];
            }
        }

        // Case 4: LP out-of-line gradient (matching CPU calcHydrogenBonds)
        if (ct == 4 && nbb_lp_c4 > 0) {
            double lpterm = -rdamp * Q_H * damp_outl * outl_nb_tot * eangl * etors * const_val;
            double rblp = lp_dist_c4;
            double ralpprblp = ralp_c4 + rblp + 1e-12;
            double ratio2_lp = exp(expo_lp_c4);

            // LP out-of-line: rab
            double tmp3 = -2.0 * lpterm * ratio2_lp * expo_lp_c4
                        / ((1.0 + ratio2_lp) * (1.0 + ratio2_lp))
                        / (ralpprblp - r_AB);
            double gi_lp = -tmp3 * ralpprblp / rab2;
            ga_x += gi_lp * drab_x; ga_y += gi_lp * drab_y; ga_z += gi_lp * drab_z;
            gb_x -= gi_lp * drab_x; gb_y -= gi_lp * drab_y; gb_z -= gi_lp * drab_z;

            // LP out-of-line: ralp  (dralp = pos_A - lp_pos)
            double dralp_x = ax - lpx_c4, dralp_y = ay - lpy_c4, dralp_z = az - lpz_c4;
            gi_lp = tmp3 / (ralp_c4 + 1e-12);
            double dga_lp_x = gi_lp * dralp_x, dga_lp_y = gi_lp * dralp_y, dga_lp_z = gi_lp * dralp_z;
            ga_x += dga_lp_x; ga_y += dga_lp_y; ga_z += dga_lp_z;
            // Fortran: gb -= dga (uses dga, not dgb)
            gb_x -= dga_lp_x; gb_y -= dga_lp_y; gb_z -= dga_lp_z;

            double glp_x = -dga_lp_x, glp_y = -dga_lp_y, glp_z = -dga_lp_z;

            // LP neighbor chain rule
            if (vnorm_c4 > 1e-10) {
                // gii matrix × glp vector (3×3 matmul)
                // gii.col(c) = -lp_dist * nbb * (unit_c/vnorm + lp_v * lp_v(c) / vnorm^3)
                // gnb = gii * glp
                double inv_v = 1.0 / vnorm_c4;
                double inv_v3 = inv_v * inv_v * inv_v;
                double fac = -lp_dist_c4 * (double)nbb_lp_c4;

                // gii[r][c] = fac * (-delta(r,c)/vnorm + lp_v[r]*lp_v[c]/vnorm^3)
                // gnb[r] = sum_c gii[r][c] * glp[c]
                double gnb_x = 0, gnb_y = 0, gnb_z = 0;
                for (int c = 0; c < 3; ++c) {
                    double glp_c = (c==0) ? glp_x : (c==1) ? glp_y : glp_z;
                    double unit_c = (c==0) ? -1.0 : 0.0;
                    double lp_vc = (c==0) ? lp_vx_c4 : (c==1) ? lp_vy_c4 : lp_vz_c4;
                    // gii column c, row x
                    gnb_x += fac * (((c==0)?-1.0:0.0) * inv_v + lp_vx_c4 * lp_vc * inv_v3) * glp_c;
                    gnb_y += fac * (((c==1)?-1.0:0.0) * inv_v + lp_vy_c4 * lp_vc * inv_v3) * glp_c;
                    gnb_z += fac * (((c==2)?-1.0:0.0) * inv_v + lp_vz_c4 * lp_vc * inv_v3) * glp_c;
                }
                gb_x += gnb_x; gb_y += gnb_y; gb_z += gnb_z;
                // Distribute to neighbors
                double gnb_share_x = gnb_x / (double)nbb_lp_c4;
                double gnb_share_y = gnb_y / (double)nbb_lp_c4;
                double gnb_share_z = gnb_z / (double)nbb_lp_c4;
                int nb_off = nb_B_offset[tid];
                for (int ni = 0; ni < nbb_lp_c4; ++ni) {
                    int nb = nb_B_flat[nb_off + ni];
                    add_grad(grad, nb, -gnb_share_x, -gnb_share_y, -gnb_share_z);
                }
            }
        }

    } else {
        // Case 1 gradient
        double caa = Q_A * basicity_A[tid];
        double cbb = Q_B * basicity_B[tid];

        double rterm = -aci * rdamp * qhoutl;
        double dterm = -aci * bas * qhoutl;
        double sterm = -rdamp * bas * qhoutl;
        double aterm = -aci * bas * rdamp * Q_H;

        double denom_val = 1.0 / (r_AH_4 + r_HB_4);
        double tmp_c1 = denom_val * denom_val * 4.0;
        double dd24a = rah2 * r_HB_4 * tmp_c1;
        double dd24b = rbh2 * r_AH_4 * tmp_c1;

        // bas: rah
        double gi = (caa - cbb) * dd24a * rterm;
        ga_x = gi * drah_x; ga_y = gi * drah_y; ga_z = gi * drah_z;
        gi = (cbb - caa) * dd24b * rterm;
        gb_x = gi * drbh_x; gb_y = gi * drbh_y; gb_z = gi * drbh_z;
        gh_x = -ga_x - gb_x; gh_y = -ga_y - gb_y; gh_z = -ga_z - gb_z;

        // aci
        gi = (acidity_B[tid] - acidity_A[tid]) * dd24a;
        double dga_aci_x = gi * drah_x * sterm, dga_aci_y = gi * drah_y * sterm, dga_aci_z = gi * drah_z * sterm;
        ga_x += dga_aci_x; ga_y += dga_aci_y; ga_z += dga_aci_z;
        gi = (acidity_A[tid] - acidity_B[tid]) * dd24b;
        double dgb_aci_x = gi * drbh_x * sterm, dgb_aci_y = gi * drbh_y * sterm, dgb_aci_z = gi * drbh_z * sterm;
        gb_x += dgb_aci_x; gb_y += dgb_aci_y; gb_z += dgb_aci_z;
        gh_x += -dga_aci_x - dgb_aci_x; gh_y += -dga_aci_y - dgb_aci_y; gh_z += -dga_aci_z - dgb_aci_z;

        // Damping: rab
        gi = rdamp * (ddamp - 3.0) / rab2;
        double dg_x = gi * drab_x * dterm, dg_y = gi * drab_y * dterm, dg_z = gi * drab_z * dterm;
        ga_x += dg_x; ga_y += dg_y; ga_z += dg_z;
        gb_x -= dg_x; gb_y -= dg_y; gb_z -= dg_z;

        // Out-of-line: rab
        gi = aterm * 2.0 * ratio2 * expo * rahprbh
           / ((1.0 + ratio2) * (1.0 + ratio2))
           / (rahprbh - r_AB) / rab2;
        dg_x = gi * drab_x; dg_y = gi * drab_y; dg_z = gi * drab_z;
        ga_x += dg_x; ga_y += dg_y; ga_z += dg_z;
        gb_x -= dg_x; gb_y -= dg_y; gb_z -= dg_z;

        // Out-of-line: rah, rbh
        double tmp_outl = -2.0 * aterm * ratio2 * expo
                        / ((1.0 + ratio2) * (1.0 + ratio2))
                        / (rahprbh - r_AB);
        double dga_o_x = drah_x * tmp_outl / r_AH, dga_o_y = drah_y * tmp_outl / r_AH, dga_o_z = drah_z * tmp_outl / r_AH;
        ga_x += dga_o_x; ga_y += dga_o_y; ga_z += dga_o_z;
        double dgb_o_x = drbh_x * tmp_outl / r_HB, dgb_o_y = drbh_y * tmp_outl / r_HB, dgb_o_z = drbh_z * tmp_outl / r_HB;
        gb_x += dgb_o_x; gb_y += dgb_o_y; gb_z += dgb_o_z;
        gh_x += -dga_o_x - dgb_o_x; gh_y += -dga_o_y - dgb_o_y; gh_z += -dga_o_z - dgb_o_z;
    }

    add_grad(grad, A, ga_x, ga_y, ga_z);
    add_grad(grad, B, gb_x, gb_y, gb_z);
    add_grad(grad, H, gh_x, gh_y, gh_z);

    } // expo <= 15.0
    } // r_AB in range
    } // tid < n

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// GPU-only postprocess kernels (replace postProcessCPU)
// Claude Generated (March 2026): Full GPU gradient consistency
// ============================================================================

// ============================================================================
// Kernel: Coulomb TERM 2+3 self-energy (O(N), energy only)
// E_en = -Σ qi * chi_eff_i
// E_self = 0.5 * Σ qi^2 * (gam_i + sqrt(2/pi) / sqrt(alp_i))
// Reference: ff_workspace.cpp::postProcess() lines 354-383
// ============================================================================
__global__ void k_coulomb_self(
    int N,
    const double* __restrict__ eeq_charges,
    const double* __restrict__ chi_base,
    const double* __restrict__ cnf,
    const double* __restrict__ cn,
    const double* __restrict__ gam,
    const double* __restrict__ alp,
    double*                    energy)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    double q = eeq_charges[i];
    if (isnan(q)) return;
    if (alp[i] <= 0.0) return;

    // chi_eff = chi_base + cnf * sqrt(max(cn, 0))
    double chi;
    if (cnf[i] != 0.0) {
        chi = chi_base[i] + cnf[i] * sqrt(fmax(cn[i], 0.0));
    } else {
        chi = chi_base[i];
    }

    static const double sqrt_2_over_pi = 0.797884560802865;
    double E_en   = -q * chi;
    double E_self = 0.5 * q * q * (gam[i] + sqrt_2_over_pi / sqrt(alp[i]));

    atomicAdd(energy, E_en + E_self);
}

// ============================================================================
// Kernel: Subtract qtmp from dEdcn (Coulomb TERM 1b chain-rule correction)
// dEdcn[i] -= q[i] * cnf[i] / (2*sqrt(cn[i]) + 1e-16)
// Reference: ff_workspace.cpp::postProcess() lines 391-399
// ============================================================================
__global__ void k_subtract_qtmp(
    int N,
    const double* __restrict__ eeq_charges,
    const double* __restrict__ cnf,
    const double* __restrict__ cn,
    double*                    dEdcn)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    double cn_i = fmax(cn[i], 0.0);
    double qtmp = eeq_charges[i] * cnf[i] / (2.0 * sqrt(cn_i) + 1e-16);
    dEdcn[i] -= qtmp;
}

// ============================================================================
// Fused Kernel: Coulomb self-energy + qtmp subtraction in single O(N) pass
// Claude Generated (March 2026): Kernel fusion optimization
//
// Combines k_coulomb_self (TERM 2+3 self-energy, energy accumulation) with
// k_subtract_qtmp (TERM 1b chain-rule correction, dEdcn modification).
// Saves one kernel launch + fixes pattern inconsistency (now uses
// blockReduceAddEnergy instead of per-thread atomicAdd).
//
// Reference: ff_workspace.cpp::postProcess() lines 354-399
// ============================================================================
__global__ GFNFF_KERNEL_BOUNDS void k_coulomb_postprocess(
    int N,
    const double* __restrict__ eeq_charges,
    const double* __restrict__ chi_base,
    const double* __restrict__ cnf,
    const double* __restrict__ cn,
    const double* __restrict__ gam,
    const double* __restrict__ alp,
    double*                    dEdcn,
    double*                    energy,
    bool                       do_subtract)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (i < N) {
        double q = eeq_charges[i];
        bool valid = !isnan(q) && alp[i] > 0.0;

        if (valid) {
            // Coulomb TERM 2+3 self-energy
            double chi;
            if (cnf[i] != 0.0) {
                chi = chi_base[i] + cnf[i] * sqrt(fmax(cn[i], 0.0));
            } else {
                chi = chi_base[i];
            }
            static const double sqrt_2_over_pi = 0.797884560802865;
            double E_en   = -q * chi;
            double E_self = 0.5 * q * q * (gam[i] + sqrt_2_over_pi / sqrt(alp[i]));
            local_E = E_en + E_self;
        }

        // Coulomb TERM 1b chain-rule: dEdcn[i] -= qtmp
        if (do_subtract) {
            double cn_i = fmax(cn[i], 0.0);
            double qtmp = q * cnf[i] / (2.0 * sqrt(cn_i) + 1e-16);
            dEdcn[i] -= qtmp;
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Kernel: CN chain-rule gradient (pairwise)
// For each pair (i,j):
//   dr = (rij - rcov_sum) / rcov_sum
//   dS/dr = (kn / sqrt(pi)) * exp(-kn^2 * dr^2) / rcov_sum
//   fac = dS/dr / rij * (dEdcn[i]*dlogdcn[i] + dEdcn[j]*dlogdcn[j])
//   grad_i += fac * (ri - rj),  grad_j -= fac * (ri - rj)
// Reference: gfnff_method.cpp:calculateCoordinationNumberDerivatives
// ============================================================================
__global__ void k_cn_chainrule(
    int n_pairs,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ rcov_sum,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ dEdcn,
    const double* __restrict__ dlogdcn,
    double*                    grad,
    double        kn)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pairs) return;

    int i = idx_i[tid], j = idx_j[tid];
    double dx = cx[i] - cx[j];
    double dy = cy[i] - cy[j];
    double dz = cz[i] - cz[j];
    applyMIC(dx, dy, dz);
    double r2 = dx*dx + dy*dy + dz*dz;
    double rij = sqrt(r2);
    if (rij < 1e-10) return;

    double rcov = rcov_sum[tid];
    double dr = (rij - rcov) / rcov;

    // dS/dr = (kn / sqrt(pi)) * exp(-kn^2 * dr^2) / rcov
    static const double inv_sqrtpi = 0.5641895835477563;  // 1/sqrt(pi)
    double dSdr = kn * inv_sqrtpi * exp(-kn * kn * dr * dr) / rcov;

    // Chain-rule gradient for CN-dependent energy terms
    // Reference: Fortran gfnff_cn.f90 gfnff_dlogcoord lines 112-115
    // For pair (i,j): dlogCN(i,i) -= dlogdcn_i*rij, dlogCN(j,j) += dlogdcn_j*rij
    //                 dlogCN(i,j) = -dlogdcn_j*rij, dlogCN(j,i) = +dlogdcn_i*rij
    // where rij = dSdr * (rj - ri) / |rij|
    // Gradient: g(i) = -(dlogdcn_i*dEdcn(i) + dlogdcn_j*dEdcn(j)) * rij
    //           g(j) = +(dlogdcn_i*dEdcn(i) + dlogdcn_j*dEdcn(j)) * rij
    //
    // On GPU: dx = xi - xj = -(xj - xi), so:
    // g(i) = +(dlogdcn_i*dEdcn(i) + dlogdcn_j*dEdcn(j)) * dSdr * dx / rij
    // g(j) = -(dlogdcn_i*dEdcn(i) + dlogdcn_j*dEdcn(j)) * dSdr * dx / rij
    double fac = dSdr / rij * (dEdcn[i] * dlogdcn[i] + dEdcn[j] * dlogdcn[j]);

    add_grad(grad, i,  fac * dx,  fac * dy,  fac * dz);
    add_grad(grad, j, -fac * dx, -fac * dy, -fac * dz);
}

// ============================================================================
// Kernel: HB alpha-modulation chain-rule gradient (pairwise over HB neighbor pairs)
// Claude Generated (March 2026): Missing gradient term for egbond_hb
//
// For bonds with HB modification: alpha_mod = (1 - t1*hb_cn_H)*alpha_orig
// The chain-rule gives: dE/d(hb_cn_H) = t1 * alpha_orig * dr^2 * E = zz
// This kernel applies: grad_H += zz_H * d(hb_cn_H)/dx_H
//                      grad_B += zz_H * d(hb_cn_H)/dx_B
// where d(hb_cn_H)/dx uses error function CN with kn=27.5, rcov_scal=1.78
//
// Reference: Fortran gfnff_engrad.F90:1054-1063 (egbond_hb, hb_dcn term)
// ============================================================================
__global__ void k_hb_alpha_chainrule(
    int n_pairs,
    const int*    __restrict__ idx_H,
    const int*    __restrict__ idx_B,
    const double* __restrict__ rcov_sum,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ zz_hb,
    double*                    grad,
    double        kn)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pairs) return;

    int H = idx_H[tid], B = idx_B[tid];
    double dx = cx[H] - cx[B];
    double dy = cy[H] - cy[B];
    double dz = cz[H] - cz[B];
    applyMIC(dx, dy, dz);
    double r2 = dx*dx + dy*dy + dz*dz;
    double rij = sqrt(r2);
    if (rij < 1e-10) return;

    double rcov = rcov_sum[tid];
    double dr = (rij - rcov) / rcov;

    // dS/dr = (-kn / (rcov * sqrt(pi))) * exp(-(kn*dr)^2)
    // Note: kn=27.5 for HB CN, sign chosen so gradient pushes correctly
    static const double inv_sqrtpi = 0.5641895835477563;  // 1/sqrt(pi)
    double dSdr = -kn * inv_sqrtpi * exp(-kn * kn * dr * dr) / rcov;

    // grad_H += zz_H * dS/dr * (H-B)/r,  grad_B -= same
    double fac = zz_hb[H] * dSdr / rij;
    add_grad(grad, H,  fac*dx,  fac*dy,  fac*dz);
    add_grad(grad, B, -fac*dx, -fac*dy, -fac*dz);
}

// ============================================================================
// Kernel: CN Compute - GFN-FF coordination number calculation
// Claude Generated (March 2026): GPU implementation of CNCalculator::calculateGFNFFCN
//
// Computes: CN_raw[i] = sum_{j≠i} 0.5 * (1 + erf(kn * (r_ij - rcov_ij) / rcov_ij))
//           CN_final[i] = log(1 + e^cnmax) - log(1 + e^(cnmax - CN_raw[i]))
//
// Thread layout: 1 thread per atom, each reduces over all other atoms
// Uses constant memory d_rcov_d3 for covalent radii (uploaded via upload_rcov_d3)
//
// Reference: gfnff_cn.f90:66-126, Spicher & Grimme J. Chem. Theory Comput. 2020
// ============================================================================
__global__ void k_cn_compute(
    int natoms,
    const double* __restrict__ cx,         // [N] x-coordinates in Bohr (SoA)
    const double* __restrict__ cy,         // [N] y-coordinates in Bohr (SoA)
    const double* __restrict__ cz,         // [N] z-coordinates in Bohr (SoA)
    const int*    __restrict__ atom_types, // [N] 1-based atomic numbers
    double*       __restrict__ cn_raw,      // [N] output: raw CN (erf sum)
    double*       __restrict__ cn_final,    // [N] output: log-transformed CN
    double        kn,                       // decay constant (-7.5)
    double        cnmax,                    // squashing limit (4.4)
    double        threshold_sq)             // distance cutoff squared (Bohr²)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= natoms) return;

    int zi = atom_types[i];
    if (zi < 1 || zi > 86) {
        cn_raw[i] = 0.0;
        cn_final[i] = 0.0;
        return;
    }

    // Claude Generated (March 2026): CN uses 4/3-scaled covalent radii.
    // d_rcov_d3 is in Bohr WITHOUT 4/3 factor (shared with angle/dihedral damping).
    // Reference: gfnff_param.f90:381-404 — "covalentRadD3 * aatoau * 4/3"
    constexpr double CN_RCOV_SCALE = 4.0 / 3.0;
    double rcov_i = d_rcov_d3[zi - 1] * CN_RCOV_SCALE;
    double xi = cx[i];
    double yi = cy[i];
    double zi_coord = cz[i];

    double cn_sum = 0.0;

    // Loop over all other atoms
    for (int j = 0; j < natoms; ++j) {
        if (i == j) continue;

        int zj = atom_types[j];
        if (zj < 1 || zj > 86) continue;

        double dx = xi - cx[j];
        double dy = yi - cy[j];
        double dz = zi_coord - cz[j];
        applyMIC(dx, dy, dz);
        double r2 = dx*dx + dy*dy + dz*dz;

        // Distance cutoff (skip far atoms)
        if (r2 > threshold_sq) continue;

        double rcov_j = d_rcov_d3[zj - 1] * CN_RCOV_SCALE;
        double rcov_ij = rcov_i + rcov_j;

        double rij = sqrt(r2);
        double dr = (rij - rcov_ij) / rcov_ij;

        // Error function CN contribution
        cn_sum += 0.5 * (1.0 + erf(kn * dr));
    }

    cn_raw[i] = cn_sum;

    // Log transformation for numerical stability
    cn_final[i] = log(1.0 + exp(cnmax)) - log(1.0 + exp(cnmax - cn_sum));
}

// ============================================================================
// GPU dlogdcn computation (Apr 2026)
// Compute logistic squashing derivative directly on GPU after k_cn_compute.
// dlogdcn[i] = exp(cnmax) / (exp(cnmax) + exp(cn_raw[i]))
// Reference: Fortran gfnff_cn.f90 create_dlogCN
// ============================================================================

__global__ GFNFF_KERNEL_BOUNDS void k_dlogdcn(
    int natoms,
    const double* __restrict__ cn_raw,
    double*       __restrict__ dlogdcn,
    double        exp_cnmax)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= natoms) return;
    dlogdcn[i] = exp_cnmax / (exp_cnmax + exp(cn_raw[i]));
}

// ============================================================================
// DC6DCN per-pair kernel (Phase 2 GPU optimization)
// Claude Generated (March 2026): Compute dc6/dcn directly per dispersion pair.
// Eliminates O(N²) CPU matrix construction — only O(nd) pairs computed on GPU.
// Reference: d4param_generator.cpp:computeDC6DCN()
// ============================================================================

// D4_MAX_ELEM and D4_MAX_REF are defined in gfnff_kernels.cuh

// ============================================================================
// GPU Gaussian weight computation (Phase 6: March 2026)
// Claude Generated: Compute gw and dgw/dCN on GPU, eliminating CPU computation
// + flatten + sync H2D upload.
//
// 1 thread = 1 atom. Each thread computes both normalized Gaussian weights and
// their CN-derivatives in a single pass (MAX_REF=7 fits in registers).
//
// Reference: d4param_generator.cpp:precomputeGaussianWeights() (line 986)
//            d4param_generator.cpp:computeGaussianWeightDerivatives() (line 1245)
//            Fortran gfnff_gdisp0.f90:405 weight_cn()
// ============================================================================

__global__ GFNFF_KERNEL_BOUNDS void k_gaussian_weights(
    int natoms,
    const double* __restrict__ cn,
    const int*    __restrict__ atom_types,
    double*       __restrict__ gw,
    double*       __restrict__ dgw)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= natoms) return;

    constexpr double wf = 4.0;  // Gaussian width parameter (D4)

    int elem = atom_types[i] - 1;  // 0-based element index
    int base = i * D4_MAX_REF;

    // Bounds check: invalid element → zero-fill
    if (elem < 0 || elem >= D4_MAX_ELEM) {
        for (int r = 0; r < D4_MAX_REF; ++r) {
            gw[base + r]  = 0.0;
            dgw[base + r] = 0.0;
        }
        return;
    }

    // Read from constant memory (L1-broadcast for same-element threads)
    int nref = d_refn_const[elem];
    if (nref > D4_MAX_REF) nref = D4_MAX_REF;

    double cn_i = cn[i];
    int refcn_base = elem * D4_MAX_REF;

    // Phase 1: Compute raw exponential weights and sums
    double expw[D4_MAX_REF];
    double dexpw[D4_MAX_REF];
    double norm = 0.0;
    double dnorm = 0.0;

    for (int r = 0; r < nref; ++r) {
        double diff = cn_i - d_refcn_const[refcn_base + r];
        double ew = exp(-wf * diff * diff);
        double dew = 2.0 * wf * (d_refcn_const[refcn_base + r] - cn_i) * ew;
        expw[r] = ew;
        dexpw[r] = dew;
        norm += ew;
        dnorm += dew;
    }

    // Phase 2: Normalize and compute derivatives via quotient rule
    if (norm > 1e-10) {
        double inv_norm = 1.0 / norm;
        double inv_norm2 = inv_norm * inv_norm;
        for (int r = 0; r < nref; ++r) {
            gw[base + r]  = expw[r] * inv_norm;
            dgw[base + r] = (dexpw[r] * norm - expw[r] * dnorm) * inv_norm2;
        }
    } else {
        // Fallback: first reference gets weight 1.0
        gw[base] = 1.0;
        dgw[base] = 0.0;
        for (int r = 1; r < nref; ++r) {
            gw[base + r]  = 0.0;
            dgw[base + r] = 0.0;
        }
    }

    // Zero-pad remaining slots (nref..D4_MAX_REF-1)
    for (int r = nref; r < D4_MAX_REF; ++r) {
        gw[base + r]  = 0.0;
        dgw[base + r] = 0.0;
    }
}

__global__ void k_dc6dcn_per_pair(
    int n_pairs,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ atom_types,
    const double* __restrict__ gw,
    const double* __restrict__ dgw,
    const double* __restrict__ c6_flat,
    double*       __restrict__ dc6dcn_ij,
    double*       __restrict__ dc6dcn_ji)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pairs) return;

    int ai = idx_i[tid];
    int aj = idx_j[tid];

    // Element indices (0-based)
    int ei = atom_types[ai] - 1;
    int ej = atom_types[aj] - 1;

    int nri = d_refn_const[ei];
    int nrj = d_refn_const[ej];

    // C6 flat index base: elem_i * MAX_ELEM * MAX_REF² + elem_j * MAX_REF²
    size_t c6_base = (size_t)ei * D4_MAX_ELEM * D4_MAX_REF * D4_MAX_REF
                   + (size_t)ej * D4_MAX_REF * D4_MAX_REF;

    // Weight array offsets: atom * MAX_REF
    int gw_i_base = ai * D4_MAX_REF;
    int gw_j_base = aj * D4_MAX_REF;

    double dc6_ij = 0.0;
    double dc6_ji = 0.0;

    for (int ri = 0; ri < nri; ++ri) {
        double dgw_i_ri = dgw[gw_i_base + ri];
        double gw_i_ri  = gw[gw_i_base + ri];

        for (int rj = 0; rj < nrj; ++rj) {
            double c6ref = c6_flat[c6_base + ri * D4_MAX_REF + rj];
            if (fabs(c6ref) < 1e-20) continue;

            // dc6dcn(i,j) = dC6(i,j)/dCN(i) = Σ dgw(i,ri) * gw(j,rj) * C6ref
            dc6_ij += dgw_i_ri * gw[gw_j_base + rj] * c6ref;

            // dc6dcn(j,i) = dC6(i,j)/dCN(j) = Σ gw(i,ri) * dgw(j,rj) * C6ref
            dc6_ji += gw_i_ri * dgw[gw_j_base + rj] * c6ref;
        }
    }

    dc6dcn_ij[tid] = dc6_ij;
    dc6dcn_ji[tid] = dc6_ji;
}

// ============================================================================
// Utility: zero device array
// ============================================================================

// ============================================================================
// GPU CN pair list generation (Apr 2026)
// Two-pass kernels replacing CPU generateCNPairList() O(N^2) loop.
// ============================================================================

/// Pass 1: count valid pairs. 1 thread = 1 upper-triangle pair (i<j).
__global__ GFNFF_KERNEL_BOUNDS void k_generate_cn_pairs_count(
    int N,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const int*    __restrict__ atom_types,
    double        cutoff_factor,
    int*          d_count)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int n_total = N * (N - 1) / 2;
    if (tid >= n_total) return;

    // Map linear tid to upper-triangle (i,j) with i<j
    // n_pairs_before_i = i * (2*N - i - 1) / 2
    // Solve i*(2*N - i - 1)/2 <= tid for i
    int i = static_cast<int>((2.0 * N - 1.0 - sqrt((2.0 * N - 1.0) * (2.0 * N - 1.0) - 8.0 * tid)) * 0.5);
    // Correct rounding drift
    while (i > 0 && i * (2 * N - i - 1) / 2 > tid) --i;
    while (i + 1 < N && (i + 1) * (2 * N - (i + 1) - 1) / 2 <= tid) ++i;
    int offset_i = i * (2 * N - i - 1) / 2;
    int j = tid - offset_i + i + 1;

    if (j >= N) return;

    int zi = atom_types[i];
    int zj = atom_types[j];
    if (zi < 1 || zi > 86 || zj < 1 || zj > 86) return;

    double dx = cx[i] - cx[j];
    double dy = cy[i] - cy[j];
    double dz = cz[i] - cz[j];
    applyMIC(dx, dy, dz);
    double r2 = dx * dx + dy * dy + dz * dz;
    double rij = sqrt(r2);

    constexpr double CN_RCOV_SCALE = 4.0 / 3.0;
    double rcov_i = d_rcov_d3[zi - 1] * CN_RCOV_SCALE;
    double rcov_j = d_rcov_d3[zj - 1] * CN_RCOV_SCALE;
    double rcov_sum = rcov_i + rcov_j;

    if (rij < cutoff_factor * rcov_sum && rij > 1e-10) {
        atomicAdd(d_count, 1);
    }
}

/// Pass 2: write valid pairs into pre-allocated buffers.
__global__ GFNFF_KERNEL_BOUNDS void k_generate_cn_pairs_write(
    int N,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const int*    __restrict__ atom_types,
    double        cutoff_factor,
    int*          d_counter,
    int*          idx_i,
    int*          idx_j,
    double*       rcov_sum)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int n_total = N * (N - 1) / 2;
    if (tid >= n_total) return;

    // Same mapping as Pass 1
    int i = static_cast<int>((2.0 * N - 1.0 - sqrt((2.0 * N - 1.0) * (2.0 * N - 1.0) - 8.0 * tid)) * 0.5);
    while (i > 0 && i * (2 * N - i - 1) / 2 > tid) --i;
    while (i + 1 < N && (i + 1) * (2 * N - (i + 1) - 1) / 2 <= tid) ++i;
    int offset_i = i * (2 * N - i - 1) / 2;
    int j = tid - offset_i + i + 1;

    if (j >= N) return;

    int zi = atom_types[i];
    int zj = atom_types[j];
    if (zi < 1 || zi > 86 || zj < 1 || zj > 86) return;

    double dx = cx[i] - cx[j];
    double dy = cy[i] - cy[j];
    double dz = cz[i] - cz[j];
    applyMIC(dx, dy, dz);
    double r2 = dx * dx + dy * dy + dz * dz;
    double rij = sqrt(r2);

    constexpr double CN_RCOV_SCALE = 4.0 / 3.0;
    double rcov_i = d_rcov_d3[zi - 1] * CN_RCOV_SCALE;
    double rcov_j = d_rcov_d3[zj - 1] * CN_RCOV_SCALE;
    double rcov_sum_val = rcov_i + rcov_j;

    if (rij < cutoff_factor * rcov_sum_val && rij > 1e-10) {
        int idx = atomicAdd(d_counter, 1);
        idx_i[idx] = i;
        idx_j[idx] = j;
        rcov_sum[idx] = rcov_sum_val;
    }
}

// ============================================================================
// WP3: Pair-list-based CN computation (Mai 2026)
// Replaces O(N²) k_cn_compute with O(n_pairs) pair-list kernel every step.
// Pair list is built once per topology build by generateCNPairListOnGPU().
// ============================================================================

__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_cn_compute_pairs(
    int           n_pairs,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ rcov_sum,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*       __restrict__ cn_raw,
    double        kn)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pairs) return;

    int i = idx_i[tid];
    int j = idx_j[tid];

    double dx = cx[i] - cx[j];
    double dy = cy[i] - cy[j];
    double dz = cz[i] - cz[j];
    applyMIC(dx, dy, dz);
    double r = __dsqrt_rn(dx*dx + dy*dy + dz*dz);
    if (r < 1e-10) return;

    // erf-based CN contribution (same formula as k_cn_compute)
    double arg = kn * (r / rcov_sum[tid] - 1.0);
    double contrib = 0.5 * (1.0 + erf(arg));

    // Symmetric: both atoms i and j receive the contribution
    atomicAdd(&cn_raw[i], contrib);
    atomicAdd(&cn_raw[j], contrib);
}

__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_logcn(
    int natoms,
    const double* __restrict__ cn_raw,
    double*       __restrict__ cn_final,
    double        cnmax)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= natoms) return;
    // Log squashing: CN_final = log(1+exp(cnmax)) - log(1+exp(cnmax-CN_raw))
    cn_final[i] = log(1.0 + exp(cnmax)) - log(1.0 + exp(cnmax - cn_raw[i]));
}

// ============================================================================
// GPU HB CN per-bond computation (Apr 2026)
// Replaces CPU HB CN loop in gfnff_gpu_method.cpp.
// ============================================================================

__global__ GFNFF_KERNEL_BOUNDS void k_hb_cn_per_atom(
    int n_pairs,
    const int*    __restrict__ idx_H,
    const int*    __restrict__ idx_B,
    const double* __restrict__ rcov_sum,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*       __restrict__ hb_cn_per_atom,
    double        kn,
    double        thr_sq)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pairs) return;

    int H = idx_H[tid];
    int B = idx_B[tid];

    double dx = cx[H] - cx[B];
    double dy = cy[H] - cy[B];
    double dz = cz[H] - cz[B];
    applyMIC(dx, dy, dz);
    double r2 = dx * dx + dy * dy + dz * dz;
    if (r2 > thr_sq) return;

    double r = sqrt(r2);
    double rcovij = rcov_sum[tid];
    double arg = -kn * (r - rcovij) / rcovij;
    double contrib = 0.5 * (1.0 + erf(arg));

    atomicAdd(&hb_cn_per_atom[H], contrib);
}

__global__ GFNFF_KERNEL_BOUNDS void k_hb_cn_map_to_bonds(
    int nb,
    const int*    __restrict__ nr_hb,
    const int*    __restrict__ hb_H_atom,
    const double* __restrict__ hb_cn_per_atom,
    double*       __restrict__ hb_cn_H)
{
    int b = blockIdx.x * blockDim.x + threadIdx.x;
    if (b >= nb) return;

    if (nr_hb[b] >= 1) {
        int H = hb_H_atom[b];
        if (H >= 0) {
            hb_cn_H[b] = hb_cn_per_atom[H];
        } else {
            hb_cn_H[b] = 0.0;
        }
    } else {
        hb_cn_H[b] = 0.0;
    }
}

// ============================================================================
// GPU EEQ Schur complement (Apr 2026)
// Replaces CPU Schur loop after Cholesky solve in gfnff_gpu_method.cpp.
// For nfrag == 1 only (common case). nfrag > 1 falls back to CPU path.
// ============================================================================

__global__ void k_eeq_reduce_sums(
    int N,
    const double* __restrict__ d_rhs,
    double*       __restrict__ d_sums)
{
    extern __shared__ double sdata[];
    double* s_z1  = sdata;            // [blockDim.x]
    double* s_Z2  = sdata + blockDim.x; // [blockDim.x]

    int tid = threadIdx.x;
    int i = blockIdx.x * blockDim.x + tid;

    double val_z1 = 0.0;
    double val_Z2 = 0.0;
    if (i < N) {
        val_z1 = d_rhs[i];          // column 0: z1
        val_Z2 = d_rhs[N + i];      // column 1: Z2_col0
    }
    s_z1[tid] = val_z1;
    s_Z2[tid] = val_Z2;
    __syncthreads();

    // Block-level tree reduction
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            s_z1[tid] += s_z1[tid + s];
            s_Z2[tid] += s_Z2[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) {
        atomicAdd(&d_sums[0], s_z1[0]);   // Cz1
        atomicAdd(&d_sums[1], s_Z2[0]);   // S
    }
}

__global__ void k_eeq_schur_nfrag1(
    int N,
    const double* __restrict__ d_rhs,
    double        lambda,
    double*       __restrict__ charges)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    double z1 = d_rhs[i];
    double Z2 = d_rhs[N + i];
    charges[i] = z1 - Z2 * lambda;
}

// ============================================================================
// WP7-A: General Schur for nfrag > 1 (May 2026)
// k_eeq_reduce_fragment_sums: Cz1[f] = Σ_{i∈frag_f} z1[i],
//                             S[f,g] = Σ_{i∈frag_f} Z2[i,g]
// k_eeq_schur_general:        q[i] = z1[i] - Σ_g Z2[i,g] · λ[g]
// Reference: docs/GPU_WP7_EEQ_LARGE_SYSTEMS.md (Strategie A)
// ============================================================================
__global__ void k_eeq_reduce_fragment_sums(
    int N,
    int nfrag,
    const double* __restrict__ d_rhs,
    const int*    __restrict__ d_atom_frag,
    double*       __restrict__ d_Cz1,
    double*       __restrict__ d_S)
{
    // 1 thread / atom; per-thread loop over (nfrag+1) RHS columns.
    // Each atom contributes to its own fragment row in (Cz1, S) → atomicAdd.
    // Non-deterministic ordering ⇒ <1 µEh energy noise (same class as WP3 CN-pair).
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    int f = d_atom_frag[i];
    // Column 0 contributes to Cz1[f]
    atomicAdd(&d_Cz1[f], d_rhs[i]);
    // Columns 1..nfrag contribute to S[f, g] (row-major layout: f*nfrag + g)
    for (int g = 0; g < nfrag; ++g) {
        atomicAdd(&d_S[f * nfrag + g], d_rhs[(g + 1) * N + i]);
    }
}

__global__ void k_eeq_schur_general(
    int N,
    int nfrag,
    const double* __restrict__ d_rhs,
    const double* __restrict__ d_lambda,
    double*       __restrict__ charges)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    double q = d_rhs[i];  // z1[i]
    for (int g = 0; g < nfrag; ++g) {
        q -= d_rhs[(g + 1) * N + i] * d_lambda[g];
    }
    charges[i] = q;
}

// ============================================================================
// WP7-C: GPU PCG kernels (May 2026)
// All four are 1-thread-per-atom, no shared memory, no atomics.
// Used by EEQSolverHip::solveSinglePCG.
// Reference: docs/GPU_WP7_EEQ_LARGE_SYSTEMS.md (Strategie C)
// ============================================================================

// Extract A's diagonal and invert it: M_inv[i] = 1/A[i,i].
// A is column-major N×N; element (i,i) is at A[i*N + i].
__global__ void k_pcg_extract_diag_inv(
    int N,
    const double* __restrict__ d_A,
    double*       __restrict__ d_M_inv)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    double diag = d_A[i * N + i];
    d_M_inv[i] = (diag != 0.0) ? 1.0 / diag : 0.0;
}

// PCG init: r = b − A·x. d_Ax holds A·x (computed via hipblasDsymv before this call).
__global__ void k_pcg_init_residual(
    int N,
    const double* __restrict__ d_b,
    const double* __restrict__ d_Ax,
    double*       __restrict__ d_r)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    d_r[i] = d_b[i] - d_Ax[i];
}

// Apply Jacobi preconditioner: z = M_inv ⊙ r.
__global__ void k_pcg_apply_precond(
    int N,
    const double* __restrict__ d_M_inv,
    const double* __restrict__ d_r,
    double*       __restrict__ d_z)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    d_z[i] = d_M_inv[i] * d_r[i];
}

// Direction update: p_out = z + β·p_in. Single-pass — avoids hipblasDcopy + hipblasDaxpy.
__global__ void k_pcg_dir_update(
    int N,
    const double* __restrict__ d_z,
    double        beta,
    const double* __restrict__ d_p_in,
    double*       __restrict__ d_p_out)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    d_p_out[i] = d_z[i] + beta * d_p_in[i];
}

// ============================================================================
// WP2: k_build_eeq_rhs — GPU-side EEQ RHS construction
// Claude Generated (May 2026): Eliminates CPU sync for EEQ RHS per MD step.
// chi_corr and cnf are topology-constant (uploaded once by uploadEEQTopologyParams).
// Only sqrt(cn[i]) changes each step — computed directly from d_cn_final on device.
// Reference: prepareEEQParametersForGPU() in gfnff_method.cpp
// ============================================================================
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_build_eeq_rhs(
    int N,
    const double* __restrict__ d_cn,
    const double* __restrict__ d_chi_corr,
    const double* __restrict__ d_cnf,
    double*       __restrict__ d_rhs)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    double cn_i = d_cn[i];
    double sqrt_cn = (cn_i > 0.0) ? __dsqrt_rn(cn_i) : 0.0;
    d_rhs[i] = d_chi_corr[i] + d_cnf[i] * sqrt_cn;
}

// ============================================================================
// WP5-C: k_check_dc6dcn_skip — GPU-side D4 dc6dcn skip check
// Claude Generated (May 2026): Two-stage max-reduction on max|ΔCN| to decide
// whether Gaussian weights + dc6dcn can be skipped. Replaces CPU-side
// recordCNValues + shouldSkipDc6Reconstruction, eliminating per-step CN D2H memcpy.
//
// Stage 1 (k_check_dc6dcn_skip): per-block max reduction → d_block_max[]
// Stage 2 (k_check_dc6dcn_skip_final): single block reduces d_block_max[] → skip_flag.
// ============================================================================
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_check_dc6dcn_skip(
    int N,
    const double* __restrict__ cn_cur,
    const double* __restrict__ cn_ref,
    double* __restrict__ d_block_max)
{
    __shared__ double sdata[512];
    int tid = threadIdx.x;
    int i = blockIdx.x * blockDim.x + tid;

    double local_max = 0.0;
    if (i < N) {
        double delta = fabs(cn_cur[i] - cn_ref[i]);
        local_max = delta;
    }
    sdata[tid] = local_max;
    __syncthreads();

    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            if (sdata[tid + s] > sdata[tid])
                sdata[tid] = sdata[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) {
        d_block_max[blockIdx.x] = sdata[0];
    }
}

__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_check_dc6dcn_skip_final(
    int n_blocks,
    const double* __restrict__ d_block_max,
    double abs_threshold,
    int* __restrict__ skip_flag)
{
    __shared__ double sdata[512];
    int tid = threadIdx.x;

    double local_max = (tid < n_blocks) ? d_block_max[tid] : 0.0;
    sdata[tid] = local_max;
    __syncthreads();

    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            if (sdata[tid + s] > sdata[tid])
                sdata[tid] = sdata[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) {
        *skip_flag = (sdata[0] < abs_threshold) ? 1 : 0;
    }
}

// ============================================================================
// WP6: Batched per-fragment EEQ kernels (nfrag > 1)
// Claude Generated (May 2026)
// ============================================================================

/// sqrt(2/pi) — matches value in eeq_solver_gpu.cu
__device__ constexpr double TSQRT2PI_FRAG = 0.797884560802866;

__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_eeq_build_fragment_matrices(
    int           total_pairs,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ alpha,
    const double* __restrict__ gam,
    const int*    __restrict__ frag_sizes,
    const int*    __restrict__ frag_offsets_A,
    const int*    __restrict__ frag_offsets_pair,
    const int*    __restrict__ frag_atom_offsets,
    const int*    __restrict__ frag_atom_map,
    double*       __restrict__ d_A_blocks,
    int           nfrag,
    double        cutoff_sq)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= total_pairs) return;

    // Binary search: find which fragment this tid belongs to
    int f = 0;
    {
        int lo = 0, hi = nfrag - 1;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (frag_offsets_pair[mid + 1] <= tid)
                lo = mid + 1;
            else
                hi = mid;
        }
        f = lo;
    }

    int local_tid  = tid - frag_offsets_pair[f];
    int Nf         = frag_sizes[f];
    int atom_start = frag_atom_offsets[f];

    // Map local_tid → (li, lj) in lower triangle of the N_f×N_f block: li >= lj
    int li = static_cast<int>((-1.0 + sqrt(1.0 + 8.0 * static_cast<double>(local_tid))) * 0.5);
    while (li * (li + 1) / 2 > local_tid) --li;
    while ((li + 1) * (li + 2) / 2 <= local_tid) ++li;
    int lj = local_tid - li * (li + 1) / 2;

    // Global atom indices
    int gi = frag_atom_map[atom_start + li];
    int gj = frag_atom_map[atom_start + lj];

    // Base pointer for this fragment's N_f×N_f column-major block
    double* Af = d_A_blocks + frag_offsets_A[f];

    if (li == lj) {
        Af[lj * Nf + li] = gam[gi] + TSQRT2PI_FRAG / sqrt(alpha[gi]);
    } else {
        double dx   = cx[gi] - cx[gj];
        double dy   = cy[gi] - cy[gj];
        double dz   = cz[gi] - cz[gj];
        double r_sq = dx*dx + dy*dy + dz*dz;

        if (cutoff_sq > 0.0 && r_sq > cutoff_sq) {
            Af[lj * Nf + li] = 0.0;
            Af[li * Nf + lj] = 0.0;
        } else {
            double r          = sqrt(r_sq);
            double gamma_ij   = 1.0 / sqrt(alpha[gi] + alpha[gj]);
            double val        = erf(gamma_ij * r) / r;
            Af[lj * Nf + li] = val;  // lower triangle (col=lj, row=li)
            Af[li * Nf + lj] = val;  // upper triangle (col=li, row=lj)
        }
    }
}

__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_eeq_gather_rhs_fragments(
    int           N,
    const double* __restrict__ d_rhs_global,
    const int*    __restrict__ frag_atom_map,
    const int*    __restrict__ frag_atom_offsets,
    const int*    __restrict__ frag_offsets_rhs,
    const int*    __restrict__ frag_sizes,
    double*       __restrict__ d_rhs_blocks,
    int           nfrag)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= N) return;

    // Binary search: find which fragment sorted position k belongs to
    int f = 0;
    {
        int lo = 0, hi = nfrag - 1;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (frag_atom_offsets[mid + 1] <= k)
                lo = mid + 1;
            else
                hi = mid;
        }
        f = lo;
    }

    int local_k = k - frag_atom_offsets[f];
    int gi      = frag_atom_map[k];

    // Write to col0 of this fragment's RHS block (col1 = ones, pre-filled at topo upload)
    d_rhs_blocks[frag_offsets_rhs[f] + local_k] = d_rhs_global[gi];
}

// ===================== workspace (cuda/ff_workspace_gpu.cu) ================
#include "hip/hip_runtime.h"
/*
 * <FFWorkspaceHip Implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): GPU implementation of FFWorkspaceHip.
 *
 * This file defines:
 *   - FFWorkspaceHipImpl   (CUDA-internal: all SoA buffers + stream)
 *   - SoA upload methods   (DispersionSoA::upload, etc.)
 *   - FFWorkspaceHip class (constructor, calculate, postProcessCPU, getters)
 *
 * Kernel launch layout: blockDim = 256, gridDim = ceil(n/256).
 * All gradient accumulation via atomicAdd on double (requires compute ≥ 6.0).
 *
 * Reference: Spicher/Grimme J. Chem. Theory Comput. 2020 (GFN-FF)
 */

#include "ff_workspace_hip.h"
#include "gfnff_soa_hip.h"
#include "gfnff_kernels_hip.h"

#include "../gfnff_parameters.h"
#include "../forcefieldthread.h"   // Bond, Angle, Dihedral, Inversion structs
#include "src/core/energy_calculators/ff_methods/gfnff_par.h"  // covalent_rad_d3

#include <cstdlib>  // std::getenv (CURCUMA_GFNFF_GPU_FREEZE_HBQ opt-out)

#include "src/core/curcuma_logger.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstring>
#include <stdexcept>
#include <vector>

// Matrix, Vector, SpMatrix are defined in global.h (via ff_workspace_gpu.h -> ff_workspace.h)
// global.h: Matrix = Eigen::Matrix<double, Dynamic, Dynamic, RowMajor>
// global.h: Vector = Eigen::VectorXd
// global.h: SpMatrix = Eigen::SparseMatrix<double>

// Claude Generated (April 2026): Extern declarations for PBC constant memory symbols
// defined in gfnff_kernels.cu (same device link unit)
extern __constant__ double d_lattice[9];
extern __constant__ double d_lattice_inv[9];
extern __constant__ int    d_has_pbc;

// ============================================================================
// D3 Covalent radii for angle/dihedral distance damping (Bohr, WITHOUT 4/3 factor)
// Copied from gfnff_par.h since nvcc has trouble linking static const std::vector
// Reference: Pyykkö & Atsumi, Chem. Eur. J. 15, 2009, 188-197
// Values × aatoau where aatoau = 1/0.52917726 (Bohr conversion)
// ============================================================================
static const double s_rcov_d3_bohr[87] = {
    0.32/0.52917726, 0.46/0.52917726,                                         // H, He
    1.20/0.52917726, 0.94/0.52917726, 0.77/0.52917726, 0.75/0.52917726,       // Li-C
    0.71/0.52917726, 0.63/0.52917726, 0.64/0.52917726, 0.67/0.52917726,       // N-Ne
    1.40/0.52917726, 1.25/0.52917726, 1.13/0.52917726, 1.04/0.52917726,       // Na-Si
    1.10/0.52917726, 1.02/0.52917726, 0.99/0.52917726, 0.96/0.52917726,       // P-Ar
    1.76/0.52917726, 1.54/0.52917726,                                         // K, Ca
    1.33/0.52917726, 1.22/0.52917726, 1.21/0.52917726, 1.10/0.52917726,       // Sc-Cr
    1.07/0.52917726, 1.04/0.52917726, 1.00/0.52917726, 0.99/0.52917726,       // Mn-Ni
    1.01/0.52917726, 1.09/0.52917726,                                         // Cu, Zn
    1.12/0.52917726, 1.09/0.52917726, 1.15/0.52917726, 1.10/0.52917726,       // Ga-Se
    1.14/0.52917726, 1.17/0.52917726,                                         // Br, Kr
    1.64/0.52917726, 1.46/0.52917726, 1.31/0.52917726, 1.26/0.52917726,       // Rb-Sr
    1.23/0.52917726, 1.22/0.52917726, 1.22/0.52917726, 1.20/0.52917726,       // Y-Mo
    1.19/0.52917726, 1.19/0.52917726, 1.18/0.52917726, 1.17/0.52917726,       // Tc-Pd
    1.18/0.52917726, 1.20/0.52917726, 1.21/0.52917726, 1.23/0.52917726,       // Ag-Cd
    1.28/0.52917726, 1.28/0.52917726, 1.28/0.52917726, 1.27/0.52917726,       // In-Te
    1.27/0.52917726, 1.34/0.52917726,                                         // I, Xe
    1.94/0.52917726, 1.71/0.52917726, 1.58/0.52917726, 1.51/0.52917726,       // Cs-Nd
    1.44/0.52917726, 1.44/0.52917726, 1.44/0.52917726, 1.43/0.52917726,       // Pm-Eu
    1.43/0.52917726, 1.43/0.52917726, 1.43/0.52917726, 1.40/0.52917726,       // Gd-Dy
    1.39/0.52917726, 1.39/0.52917726, 1.40/0.52917726, 1.41/0.52917726,       // Ho-Hg
    1.39/0.52917726, 1.39/0.52917726, 1.38/0.52917726, 1.38/0.52917726,       // Tl-Po
    1.38/0.52917726, 1.42/0.52917726,                                         // At, Rn
    2.01/0.52917726, 1.81/0.52917726,                                         // Fr, Ra
    1.67/0.52917726, 1.58/0.52917726                                          // Ac, Th
};

// ============================================================================
// GPU thread launch helpers (Phase 6: Dynamic block size optimization)
// Claude Generated (March 2026)
//
// Adaptive block sizing based on problem size for better GPU occupancy:
//   - Small problems (<256): 32 threads (1 warp)
//   - Medium problems (<1024): 128 threads (4 warps)
//   - Normal problems (<16384): 256 threads (8 warps)
//   - Large problems: 512 threads (16 warps)
//
// Warp-aligned block sizes ensure efficient reduction in blockReduceAddEnergy.
// ============================================================================

struct LaunchConfig {
    int blockSize;
    int gridSize;
};

static inline LaunchConfig getLaunchConfig(int n_elements, int block_size_override = 0) {
    LaunchConfig cfg;
    // G3a (Apr 2026): Optional fixed block size override.
    // If block_size_override > 0, use it directly (clamped 32..512).
    // If 0, fall back to adaptive logic based on problem size.
    if (block_size_override > 0) {
        cfg.blockSize = std::min(512, std::max(32, block_size_override));
    } else {
        // Adaptive block size based on problem size
        // IMPORTANT: maxBlockSize must match GFNFF_KERNEL_BOUNDS in gfnff_kernels.cuh
        // Current: __launch_bounds__(512, 2) - max 512 threads, min 2 blocks/SM
        if (n_elements < 256) {
            cfg.blockSize = 32;    // Single warp - minimal overhead for tiny problems
        } else if (n_elements < 1024) {
            cfg.blockSize = 128;   // 4 warps - good for small problems
        } else if (n_elements < 16384) {
            cfg.blockSize = 256;   // Standard - good balance for medium problems
        } else {
            cfg.blockSize = 512;   // 512 threads - maximize occupancy for large problems
        }
    }
    // Phase 8 (Claude Generated March 2026): min 1 block so n=0 launches a no-op block
    // instead of gridSize=0 (undefined in CUDA). All kernels early-return for i >= n.
    cfg.gridSize = std::max(1, (n_elements + cfg.blockSize - 1) / cfg.blockSize);
    return cfg;
}

// Light variant for kernels declared with GFNFF_KERNEL_BOUNDS_LIGHT
// (__launch_bounds__(256, 4)). Caps blockSize at 256 so the runtime never
// rejects the launch with hipErrorInvalidValue when adaptive logic would
// pick 512 for n_elements >= 16384. Use this for: k_cn_compute_pairs,
// k_logcn, k_check_dc6dcn_skip(+final), k_build_eeq_rhs,
// k_eeq_build_fragment_matrices, k_eeq_gather_rhs_fragments.
// Discovered May 2026: large mixture.xyz crashed Phase 1 because
// k_cn_compute_pairs got block=512 from getLaunchConfig at 40907 pairs.
static inline LaunchConfig getLaunchConfigLight(int n_elements, int block_size_override = 0) {
    LaunchConfig cfg = getLaunchConfig(n_elements, block_size_override);
    if (cfg.blockSize > 256) {
        cfg.blockSize = 256;
        cfg.gridSize  = std::max(1, (n_elements + cfg.blockSize - 1) / cfg.blockSize);
    }
    return cfg;
}

// Legacy helper for backward compatibility (still uses 256 as default)
static inline int gridFor(int n, int block = 256) {
    return (n + block - 1) / block;
}

// ============================================================================
// CUDA error check (host-side, after sync)
// ============================================================================

static void checkCuda(hipError_t err, const char* msg) {
    if (err != hipSuccess)
        throw std::runtime_error(std::string(msg) + ": " + hipGetErrorString(err));
}

// ============================================================================
// SoA upload implementations (declared in gfnff_soa.h)
// ============================================================================

// ---------------------------------------------------------------------------
void DispersionSoA::upload(const std::vector<GFNFFDispersion>& v, hipStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    std::vector<int>    h_i(n), h_j(n);
    std::vector<double> h_C6(n), h_r4r2(n), h_r0sq(n), h_zeta(n), h_rcut(n);

    for (int k = 0; k < n; ++k) {
        h_i[k]    = v[k].i;
        h_j[k]    = v[k].j;
        h_C6[k]   = v[k].C6;
        h_r4r2[k] = v[k].r4r2ij;
        h_r0sq[k] = v[k].r0_squared;
        h_zeta[k] = v[k].zetac6;
        h_rcut[k] = v[k].r_cut;
    }

    idx_i.upload(h_i, stream);
    idx_j.upload(h_j, stream);
    C6.upload(h_C6, stream);
    r4r2ij.upload(h_r4r2, stream);
    r0_sq.upload(h_r0sq, stream);
    zetac6.upload(h_zeta, stream);
    r_cut.upload(h_rcut, stream);
}

// ---------------------------------------------------------------------------
void RepulsionSoA::upload(const std::vector<GFNFFRepulsion>& v, hipStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    std::vector<int>    h_i(n), h_j(n);
    std::vector<double> h_alpha(n), h_repab(n), h_rcut(n);

    for (int k = 0; k < n; ++k) {
        h_i[k]     = v[k].i;
        h_j[k]     = v[k].j;
        h_alpha[k] = v[k].alpha;
        h_repab[k] = v[k].repab;
        h_rcut[k]  = v[k].r_cut;
    }

    idx_i.upload(h_i, stream);
    idx_j.upload(h_j, stream);
    alpha.upload(h_alpha, stream);
    repab.upload(h_repab, stream);
    r_cut.upload(h_rcut, stream);
}

// ---------------------------------------------------------------------------
void CoulombSoA::upload(const std::vector<GFNFFCoulomb>& v, int natoms, hipStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    std::vector<int>    h_i(n), h_j(n);
    std::vector<double> h_gamma(n), h_rcut(n);

    for (int k = 0; k < n; ++k) {
        h_i[k]     = v[k].i;
        h_j[k]     = v[k].j;
        h_gamma[k] = v[k].gamma_ij;
        h_rcut[k]  = v[k].r_cut;
    }

    idx_i.upload(h_i, stream);
    idx_j.upload(h_j, stream);
    gamma_ij.upload(h_gamma, stream);
    r_cut.upload(h_rcut, stream);
    r_cut_scalar = h_rcut[0];   // uniform across pairs (global coulomb_pair_r_cut / eeq cutoff)

    // Build per-atom CSR adjacency for the gather kernel (k_coulomb_gather). Each undirected
    // pair (i,j) becomes two directed edges (i->j, j->i) carrying the same gamma_ij / r_cut.
    csr_N = natoms;
    std::vector<int> deg(natoms, 0);
    for (int k = 0; k < n; ++k) { ++deg[h_i[k]]; ++deg[h_j[k]]; }
    std::vector<int> off(natoms + 1, 0);
    for (int a = 0; a < natoms; ++a) off[a + 1] = off[a] + deg[a];
    const int tot = off[natoms];                       // == 2*n
    std::vector<int>    h_part(tot);
    std::vector<double> h_g(tot), h_rc(tot);
    std::vector<int>    cur(off.begin(), off.end() - 1); // write cursors per atom
    for (int k = 0; k < n; ++k) {
        const int i = h_i[k], j = h_j[k];
        const double g = h_gamma[k], rc = h_rcut[k];
        int pi = cur[i]++; h_part[pi] = j; h_g[pi] = g; h_rc[pi] = rc;
        int pj = cur[j]++; h_part[pj] = i; h_g[pj] = g; h_rc[pj] = rc;
    }
    csr_off.upload(off, stream);
    csr_partner.upload(h_part, stream);
    csr_gamma.upload(h_g, stream);
    csr_rcut.upload(h_rc, stream);
}

// ---------------------------------------------------------------------------
void BondSoA::upload(const std::vector<Bond>& v, const std::vector<int>& atom_types,
                     hipStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    std::vector<int>    h_i(n), h_j(n), h_nr_hb(n), h_hb_H(n);
    std::vector<double> h_r0(n), h_fc(n), h_alpha(n), h_rabshift(n), h_ff(n);
    std::vector<double> h_cnfak_i(n), h_cnfak_j(n), h_rb_i(n), h_rb_j(n);
    std::vector<double> h_hb_cn_H(n);

    for (int k = 0; k < n; ++k) {
        h_i[k]       = v[k].i;
        h_j[k]       = v[k].j;
        h_r0[k]      = v[k].r0_ij;
        h_fc[k]      = v[k].fc;
        h_alpha[k]   = v[k].exponent;
        h_rabshift[k]= v[k].rabshift;
        h_ff[k]      = v[k].ff;
        h_cnfak_i[k] = v[k].cnfak_i;
        h_cnfak_j[k] = v[k].cnfak_j;
        h_rb_i[k]    = v[k].r0_base_i;
        h_rb_j[k]    = v[k].r0_base_j;
        h_nr_hb[k]   = v[k].nr_hb;
        h_hb_cn_H[k] = v[k].hb_cn_H;
        // Determine which atom is H for HB alpha gradient
        if (v[k].nr_hb >= 1) {
            if (v[k].i < static_cast<int>(atom_types.size()) && atom_types[v[k].i] == 1)
                h_hb_H[k] = v[k].i;
            else if (v[k].j < static_cast<int>(atom_types.size()) && atom_types[v[k].j] == 1)
                h_hb_H[k] = v[k].j;
            else
                h_hb_H[k] = -1;
        } else {
            h_hb_H[k] = -1;
        }
    }

    idx_i.upload(h_i, stream);
    idx_j.upload(h_j, stream);
    r0.upload(h_r0, stream);
    fc.upload(h_fc, stream);
    alpha.upload(h_alpha, stream);
    rabshift.upload(h_rabshift, stream);
    ff.upload(h_ff, stream);
    cnfak_i.upload(h_cnfak_i, stream);
    cnfak_j.upload(h_cnfak_j, stream);
    r0_base_i.upload(h_rb_i, stream);
    r0_base_j.upload(h_rb_j, stream);
    nr_hb.upload(h_nr_hb, stream);
    hb_cn_H.upload(h_hb_cn_H, stream);
    hb_H_atom.upload(h_hb_H, stream);
}

// ---------------------------------------------------------------------------
void HBAlphaSoA::upload(const std::vector<int>& h_idx, const std::vector<int>& b_idx,
                        const std::vector<double>& rcov, hipStream_t stream)
{
    n = static_cast<int>(h_idx.size());
    if (n == 0) return;
    idx_H.upload(h_idx, stream);
    idx_B.upload(b_idx, stream);
    rcov_sum.upload(rcov, stream);
}

// ---------------------------------------------------------------------------
void AngleSoA::upload(const std::vector<Angle>& v,
                      const std::vector<int>& atom_types,
                      hipStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    const int nat = static_cast<int>(atom_types.size());
    std::vector<int>    h_i(n), h_j(n), h_k(n);
    std::vector<int>    h_ati(n), h_atj(n), h_atk(n);
    std::vector<double> h_fc(n), h_th(n);

    for (int m = 0; m < n; ++m) {
        h_i[m]   = v[m].i;
        h_j[m]   = v[m].j;
        h_k[m]   = v[m].k;
        h_ati[m] = (v[m].i < nat) ? atom_types[v[m].i] : 0;
        h_atj[m] = (v[m].j < nat) ? atom_types[v[m].j] : 0;
        h_atk[m] = (v[m].k < nat) ? atom_types[v[m].k] : 0;
        h_fc[m]  = v[m].fc;
        h_th[m]  = v[m].theta0_ijk;
    }

    idx_i.upload(h_i, stream); idx_j.upload(h_j, stream); idx_k.upload(h_k, stream);
    ati.upload(h_ati, stream); atj.upload(h_atj, stream); atk.upload(h_atk, stream);
    fc.upload(h_fc, stream);
    theta0.upload(h_th, stream);
}

// ---------------------------------------------------------------------------
void DihedralSoA::upload(const std::vector<Dihedral>& standard,
                         const std::vector<Dihedral>& extra,
                         const std::vector<int>& atom_types,
                         hipStream_t stream)
{
    const int ns  = static_cast<int>(standard.size());
    const int ne  = static_cast<int>(extra.size());
    n = ns + ne;
    if (n == 0) return;

    const int nat = static_cast<int>(atom_types.size());
    std::vector<int>    h_i(n), h_j(n), h_k(n), h_l(n);
    std::vector<int>    h_ati(n), h_atj(n), h_atk(n), h_atl(n);
    std::vector<double> h_V(n), h_phi0(n), h_nper(n);
    std::vector<int>    h_isnci(n);

    auto fill = [&](int offset, const std::vector<Dihedral>& v) {
        for (int m = 0; m < static_cast<int>(v.size()); ++m) {
            const int idx = offset + m;
            h_i[idx]     = v[m].i;  h_j[idx] = v[m].j;
            h_k[idx]     = v[m].k;  h_l[idx] = v[m].l;
            h_ati[idx]   = (v[m].i < nat) ? atom_types[v[m].i] : 0;
            h_atj[idx]   = (v[m].j < nat) ? atom_types[v[m].j] : 0;
            h_atk[idx]   = (v[m].k < nat) ? atom_types[v[m].k] : 0;
            h_atl[idx]   = (v[m].l < nat) ? atom_types[v[m].l] : 0;
            h_V[idx]     = v[m].V;
            h_phi0[idx]  = v[m].phi0;
            h_nper[idx]  = v[m].n;
            h_isnci[idx] = v[m].is_nci ? 1 : 0;
        }
    };

    fill(0,  standard);
    fill(ns, extra);

    idx_i.upload(h_i, stream); idx_j.upload(h_j, stream);
    idx_k.upload(h_k, stream); idx_l.upload(h_l, stream);
    ati.upload(h_ati, stream); atj.upload(h_atj, stream);
    atk.upload(h_atk, stream); atl.upload(h_atl, stream);
    V.upload(h_V, stream);
    phi0.upload(h_phi0, stream);
    n_period.upload(h_nper, stream);
    is_nci.upload(h_isnci, stream);
}

// ---------------------------------------------------------------------------
void InversionSoA::upload(const std::vector<Inversion>& v,
                           const std::vector<int>& atom_types,
                           hipStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    const int nat = static_cast<int>(atom_types.size());
    std::vector<int>    h_i(n), h_j(n), h_k(n), h_l(n), h_pt(n);
    std::vector<int>    h_ati(n), h_atj(n), h_atk(n), h_atl(n);
    std::vector<double> h_fc(n), h_om(n), h_C0(n), h_C1(n), h_C2(n);

    for (int m = 0; m < n; ++m) {
        h_i[m]  = v[m].i;  h_j[m] = v[m].j;
        h_k[m]  = v[m].k;  h_l[m] = v[m].l;
        h_ati[m] = (v[m].i < nat) ? atom_types[v[m].i] : 0;
        h_atj[m] = (v[m].j < nat) ? atom_types[v[m].j] : 0;
        h_atk[m] = (v[m].k < nat) ? atom_types[v[m].k] : 0;
        h_atl[m] = (v[m].l < nat) ? atom_types[v[m].l] : 0;
        h_fc[m] = v[m].fc;
        h_om[m] = v[m].omega0;
        h_C0[m] = v[m].C0;
        h_C1[m] = v[m].C1;
        h_C2[m] = v[m].C2;
        h_pt[m] = v[m].potential_type;
    }

    idx_i.upload(h_i, stream); idx_j.upload(h_j, stream);
    idx_k.upload(h_k, stream); idx_l.upload(h_l, stream);
    ati.upload(h_ati, stream); atj.upload(h_atj, stream);
    atk.upload(h_atk, stream); atl.upload(h_atl, stream);
    fc.upload(h_fc, stream);
    omega0.upload(h_om, stream);
    C0.upload(h_C0, stream);
    C1.upload(h_C1, stream);
    C2.upload(h_C2, stream);
    potential_type.upload(h_pt, stream);
}

// ---------------------------------------------------------------------------
void STorsionSoA::upload(const std::vector<GFNFFSTorsion>& v, hipStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    std::vector<int>    h_i(n), h_j(n), h_k(n), h_l(n);
    std::vector<double> h_eref(n);

    for (int m = 0; m < n; ++m) {
        h_i[m]    = v[m].i;
        h_j[m]    = v[m].j;
        h_k[m]    = v[m].k;
        h_l[m]    = v[m].l;
        h_eref[m] = v[m].erefhalf;
    }

    idx_i.upload(h_i, stream);  idx_j.upload(h_j, stream);
    idx_k.upload(h_k, stream);  idx_l.upload(h_l, stream);
    erefhalf.upload(h_eref, stream);
}

// ---------------------------------------------------------------------------
void BATMSoA::upload(const std::vector<GFNFFBatmTriple>& v, hipStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    std::vector<int>    h_i(n), h_j(n), h_k(n);
    std::vector<double> h_zi(n), h_zj(n), h_zk(n);

    for (int m = 0; m < n; ++m) {
        h_i[m]  = v[m].i;
        h_j[m]  = v[m].j;
        h_k[m]  = v[m].k;
        h_zi[m] = v[m].zb3atm_i;
        h_zj[m] = v[m].zb3atm_j;
        h_zk[m] = v[m].zb3atm_k;
    }

    idx_i.upload(h_i, stream);  idx_j.upload(h_j, stream);  idx_k.upload(h_k, stream);
    zb3atm_i.upload(h_zi, stream);
    zb3atm_j.upload(h_zj, stream);
    zb3atm_k.upload(h_zk, stream);
}

// ---------------------------------------------------------------------------
void ATMSoA::upload(const std::vector<ATMTriple>& v,
                     const std::vector<int>& atom_types,
                     hipStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    const int nat = static_cast<int>(atom_types.size());
    std::vector<int>    h_i(n), h_j(n), h_k(n);
    std::vector<int>    h_ati(n), h_atj(n), h_atk(n);
    std::vector<double> h_c6ij(n), h_c6ik(n), h_c6jk(n);
    std::vector<double> h_s9(n), h_a1(n), h_a2(n), h_alp(n), h_tscale(n);

    for (int m = 0; m < n; ++m) {
        h_i[m]      = v[m].i;
        h_j[m]      = v[m].j;
        h_k[m]      = v[m].k;
        h_ati[m]    = (v[m].i < nat) ? atom_types[v[m].i] : 0;
        h_atj[m]    = (v[m].j < nat) ? atom_types[v[m].j] : 0;
        h_atk[m]    = (v[m].k < nat) ? atom_types[v[m].k] : 0;
        h_c6ij[m]   = v[m].C6_ij;
        h_c6ik[m]   = v[m].C6_ik;
        h_c6jk[m]   = v[m].C6_jk;
        h_s9[m]     = v[m].s9;
        h_a1[m]     = v[m].a1;
        h_a2[m]     = v[m].a2;
        h_alp[m]    = v[m].alp;
        h_tscale[m] = v[m].triple_scale;
    }

    idx_i.upload(h_i, stream);  idx_j.upload(h_j, stream);  idx_k.upload(h_k, stream);
    ati.upload(h_ati, stream);  atj.upload(h_atj, stream);  atk.upload(h_atk, stream);
    C6_ij.upload(h_c6ij, stream);
    C6_ik.upload(h_c6ik, stream);
    C6_jk.upload(h_c6jk, stream);
    s9.upload(h_s9, stream);
    a1.upload(h_a1, stream);
    a2.upload(h_a2, stream);
    alp.upload(h_alp, stream);
    triple_scale.upload(h_tscale, stream);
}

// ---------------------------------------------------------------------------
void XBondSoA::upload(const std::vector<GFNFFHalogenBond>& v,
                       const std::vector<int>& atom_types,
                       hipStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    const int nat = static_cast<int>(atom_types.size());
    std::vector<int>    h_i(n), h_j(n), h_k(n);
    std::vector<int>    h_eA(n), h_eB(n);
    std::vector<double> h_qX(n), h_qB(n), h_aX(n), h_rcut(n);

    for (int m = 0; m < n; ++m) {
        h_i[m]    = v[m].i;   // donor A
        h_j[m]    = v[m].j;   // halogen X
        h_k[m]    = v[m].k;   // acceptor B
        h_eA[m]   = (v[m].i < nat) ? atom_types[v[m].i] : 0;
        h_eB[m]   = (v[m].k < nat) ? atom_types[v[m].k] : 0;
        h_qX[m]   = v[m].q_X;
        h_qB[m]   = v[m].q_B;
        h_aX[m]   = v[m].acidity_X;
        h_rcut[m] = v[m].r_cut;
    }

    idx_i.upload(h_i, stream);  idx_j.upload(h_j, stream);  idx_k.upload(h_k, stream);
    elem_A.upload(h_eA, stream);  elem_B.upload(h_eB, stream);
    q_X.upload(h_qX, stream);
    q_B.upload(h_qB, stream);
    acidity_X.upload(h_aX, stream);
    r_cut.upload(h_rcut, stream);
}

// ---------------------------------------------------------------------------
void HBondSoA::upload(const std::vector<GFNFFHydrogenBond>& v,
                       const std::vector<int>& atom_types,
                       hipStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    const int nat = static_cast<int>(atom_types.size());
    std::vector<int>    h_i(n), h_j(n), h_k(n);
    std::vector<int>    h_eA(n), h_eB(n), h_case(n);
    std::vector<double> h_qH(n), h_qA(n), h_qB(n);
    std::vector<double> h_basA(n), h_basB(n), h_aciA(n), h_aciB(n);
    std::vector<double> h_rcut(n);
    std::vector<int>    h_nbB_offset(n), h_nbB_count(n);
    std::vector<int>    h_acc_parent(n);
    std::vector<int>    h_nbC_offset(n), h_nbC_count(n);
    std::vector<double> h_repz(n);

    // First pass: count total neighbors
    int total_B = 0, total_C = 0;
    for (int m = 0; m < n; ++m) {
        total_B += static_cast<int>(v[m].neighbors_B.size());
        total_C += static_cast<int>(v[m].neighbors_C.size());
    }

    std::vector<int> h_nbB_flat(total_B);
    std::vector<int> h_nbC_flat(total_C);
    int offB = 0, offC = 0;

    for (int m = 0; m < n; ++m) {
        h_i[m]     = v[m].i;   // donor A
        h_j[m]     = v[m].j;   // hydrogen H
        h_k[m]     = v[m].k;   // acceptor B
        h_eA[m]    = (v[m].i < nat) ? atom_types[v[m].i] : 0;
        h_eB[m]    = (v[m].k < nat) ? atom_types[v[m].k] : 0;
        h_qH[m]    = v[m].q_H;
        h_qA[m]    = v[m].q_A;
        h_qB[m]    = v[m].q_B;
        h_basA[m]  = v[m].basicity_A;
        h_basB[m]  = v[m].basicity_B;
        h_aciA[m]  = v[m].acidity_A;
        h_aciB[m]  = v[m].acidity_B;
        h_rcut[m]  = v[m].r_cut;
        h_case[m]  = v[m].case_type;
        h_acc_parent[m] = v[m].acceptor_parent_index;

        // repz_B for case 4
        int zB = (v[m].k < nat) ? atom_types[v[m].k] : 0;
        h_repz[m] = (zB >= 1 && zB <= 86) ? GFNFFParameters::repz[zB - 1] : 1.0;

        // Pack neighbors_B
        h_nbB_offset[m] = offB;
        h_nbB_count[m]  = static_cast<int>(v[m].neighbors_B.size());
        for (int nb : v[m].neighbors_B) h_nbB_flat[offB++] = nb;

        // Pack neighbors_C
        h_nbC_offset[m] = offC;
        h_nbC_count[m]  = static_cast<int>(v[m].neighbors_C.size());
        for (int nb : v[m].neighbors_C) h_nbC_flat[offC++] = nb;
    }

    total_nb_B = total_B;
    total_nb_C = total_C;

    idx_i.upload(h_i, stream);  idx_j.upload(h_j, stream);  idx_k.upload(h_k, stream);
    elem_A.upload(h_eA, stream);  elem_B.upload(h_eB, stream);
    q_H.upload(h_qH, stream);  q_A.upload(h_qA, stream);  q_B.upload(h_qB, stream);
    basicity_A.upload(h_basA, stream);  basicity_B.upload(h_basB, stream);
    acidity_A.upload(h_aciA, stream);   acidity_B.upload(h_aciB, stream);
    r_cut.upload(h_rcut, stream);
    case_type.upload(h_case, stream);
    nb_B_offset.upload(h_nbB_offset, stream);
    nb_B_count.upload(h_nbB_count, stream);
    if (total_B > 0) nb_B_flat.upload(h_nbB_flat, stream);
    acceptor_parent.upload(h_acc_parent, stream);
    nb_C_offset.upload(h_nbC_offset, stream);
    nb_C_count.upload(h_nbC_count, stream);
    if (total_C > 0) nb_C_flat.upload(h_nbC_flat, stream);
    repz_B.upload(h_repz, stream);
}

// ============================================================================
// FFWorkspaceHipImpl — CUDA-internal state (not visible outside this .cu)
// ============================================================================

struct FFWorkspaceHipImpl {
    // Static SoA buffers (uploaded once at construction)
    DispersionSoA  disp;            ///< D4 dispersion pairs
    RepulsionSoA   bonded_rep;      ///< Bonded repulsion pairs
    RepulsionSoA   nonbonded_rep;   ///< Non-bonded repulsion pairs
    CoulombSoA     coulomb;         ///< Coulomb pairs (gamma + cutoff only; charges dynamic)
    BondSoA        bonds;           ///< Bond stretching
    AngleSoA       angles;          ///< Angle bending
    DihedralSoA    dihedrals;       ///< Standard + extra torsions (combined)
    InversionSoA   inversions;      ///< Out-of-plane inversions

    // Phase 2 SoA buffers (March 2026): formerly CPU residual, now all GPU
    STorsionSoA    storsions;       ///< Triple bond torsions
    BATMSoA        batm;            ///< Bonded ATM 3-body
    ATMSoA         atm;             ///< Axilrod-Teller-Muto 3-body dispersion
    XBondSoA       xbonds;          ///< Halogen bonds
    HBondSoA       hbonds;          ///< Hydrogen bonds

    // Dynamic buffers (reset + upload each calculation step)
    CoordSoA    coords;      ///< Current geometry SoA (x[N], y[N], z[N]) — Phase 10, March 2026
    CudaBuffer<double> d_charges;   ///< [N] EEQ charges
    CudaBuffer<double> d_cn;        ///< [N] D3 coordination numbers
    CudaBuffer<double> d_topo_charges; ///< [N] topology charges (for BATM)
    CudaBuffer<double> d_grad;      ///< [N*3] gradient output (atomicAdd)
    CudaBuffer<double> d_dEdcn;     ///< [N] bond dE/dCN (atomicAdd)
    CudaBuffer<double> d_dlogdcn;  ///< [N] logistic squashing factor for CN chain-rule
    CudaBuffer<double> d_energy;    ///< [1] total GPU energy (atomicAdd)

    // CN chain-rule pair list (uploaded once, used every step)
    CudaBuffer<int>    d_cn_idx_i;     ///< [n_cn_pairs]
    CudaBuffer<int>    d_cn_idx_j;     ///< [n_cn_pairs]
    CudaBuffer<double> d_cn_rcov_sum;  ///< [n_cn_pairs]
    int                n_cn_pairs = 0;
    CudaBuffer<int>    d_cn_pair_counter; ///< [1] atomic counter for on-GPU pair generation (Apr 2026)

    // HB alpha chain-rule: per-atom zz accumulator + HB neighbor pair list
    CudaBuffer<double> d_zz_hb;        ///< [N] per-atom zz accumulator (zeroed each step)
    HBAlphaSoA         hb_alpha;       ///< (H, B) pairs for alpha gradient
    CudaBuffer<double> d_hb_cn_per_atom; ///< [N] per-H accumulator for HB CN (Apr 2026)

    // Coulomb self-energy parameters (uploaded once)
    // WP5-B (May 2026): d_coul_cnf removed — topology-constant eeq_topo.d_cnf is used
    // directly by k_coulomb_postprocess, eliminating per-step H2D copy.
    CudaBuffer<double> d_coul_chi_base; ///< [N]
    CudaBuffer<double> d_coul_gam;      ///< [N]
    CudaBuffer<double> d_coul_alp;      ///< [N]
    bool               coul_self_on_gpu = false;
    bool               use_mixed_precision = false;  ///< TODO: Restore to true after GPU gradient debugging complete

    // GPU topology displacement check buffers (Claude Generated March 2026)
    RefCoordSoA ref_coords;  ///< Reference geometry SoA (rx[N], ry[N], rz[N]) — Phase 10, March 2026
    CudaBuffer<int>    d_disp_flag;       ///< [1] exceeded flag for displacement check

    // GPU CN computation buffers (Phase 1: GPU migration)
    // Claude Generated (March 2026): CN computed entirely on GPU
    CudaBuffer<double> d_cn_raw;        ///< [N] raw CN values (erf sum)
    CudaBuffer<double> d_cn_final;       ///< [N] log-transformed CN values
    CudaBuffer<int>    d_atom_types;    ///< [N] atomic numbers for CN lookup

    // WP5-C (May 2026): GPU-side D4 dc6dcn skip-check buffers
    CudaBuffer<double> d_cn_d4_ref;        ///< [N] CN snapshot from previous step
    CudaBuffer<double> d_dc6dcn_block_max; ///< [max_blocks] per-block max for reduction
    CudaBuffer<int>    d_dc6dcn_skip;      ///< [1] skip flag (0 = recompute, 1 = skip)
    int*               h_dc6dcn_skip = nullptr; ///< Pinned host buffer for async D2H

    // Per-term energy accumulators — consolidated into single contiguous buffer
    // Claude Generated (March 2026): Reduces 14 individual hipMemcpy to 1
    enum EnergySlot {
        E_DISP = 0, E_BREP = 1, E_NBREP = 2, E_COUL = 3, E_BOND = 4,
        E_ANGLE = 5, E_DIHED = 6, E_INV = 7, E_STORS = 8, E_BATM = 9,
        E_ATM = 10, E_XBOND = 11, E_HBOND = 12, E_COUL_SELF = 13,
        N_ENERGY_SLOTS = 14
    };
    CudaBuffer<double> d_energies;      ///< [N_ENERGY_SLOTS] all per-term energies

    // Phase 2: GPU dc6dcn per-pair computation (Claude Generated March 2026)
    // Gaussian weights and C6 reference table for on-GPU dc6dcn computation.
    // Eliminates O(N²) CPU computeDC6DCN() — only O(nd) pairs computed.
    CudaBuffer<double> d_gw;           ///< [N * MAX_REF] Gaussian weights (padded)
    CudaBuffer<double> d_dgw;          ///< [N * MAX_REF] weight derivatives
    CudaBuffer<double> d_c6_flat;      ///< [MAX_ELEM² * MAX_REF²] C6 reference table
    // d_refn and d_refcn moved to __constant__ memory (Phase 8: March 2026)
    // Uploaded via upload_refn_const() / upload_refcn_const() at init
    bool               dc6dcn_gpu_ready = false; ///< true after C6 table + refn uploaded

    // Diagnostic snapshot buffers (GPU-side, no pipeline stall)
    // Claude Generated (March 2026): Device-to-device copy instead of sync+download
    CudaBuffer<double> d_dEdcn_snapshot;   ///< [N] copy of dEdcn before qtmp
    CudaBuffer<double> d_grad_snapshot;    ///< [3*N] copy of gradient before CN chain-rule

    // Per-kernel gradient snapshots (Claude Generated April 2026)
    CudaBuffer<double> d_grad_after_sA;       ///< [3*N] snapshot after stream A (repulsion)
    CudaBuffer<double> d_grad_after_bonds;   ///< [3*N] snapshot after k_bonds
    CudaBuffer<double> d_grad_after_angles;  ///< [3*N] snapshot after k_angles
    CudaBuffer<double> d_grad_after_sB;      ///< [3*N] snapshot after all stream B kernels
    CudaBuffer<double> d_grad_after_sC;      ///< [3*N] snapshot after stream C (3-body)

    // WP2: topology-constant EEQ parameters on GPU (Claude Generated May 2026)
    // Uploaded once per topology build; per-step k_build_eeq_rhs reads d_chi_corrected
    // and d_cnf to construct d_rhs_atoms from d_cn_final without CPU involvement.
    struct EEQTopologyBuffers {
        CudaBuffer<double> d_alpha_corrected; ///< [N] charge-corrected alpha²
        CudaBuffer<double> d_gam_corrected;   ///< [N] gam + dgam
        CudaBuffer<double> d_chi_corrected;   ///< [N] -chi + dxi + amide_corr (no CN term)
        CudaBuffer<double> d_cnf;             ///< [N] cnf_eeq per atom
        CudaBuffer<double> d_rhs_atoms;       ///< [N] output of k_build_eeq_rhs (per-step)
        CudaBuffer<double> d_rhs_constraints; ///< [nfrag] fragment target charges
        int  nfrag = 1;
        bool valid = false;
    };
    EEQTopologyBuffers eeq_topo;

    int N = 0;
    hipStream_t stream = nullptr;

    // Multi-stream concurrency for independent kernel groups
    // Claude Generated (March 2026): Parallel execution of pairwise/bonded/3-body terms
    hipStream_t stream_pairwise = nullptr;  ///< Stream A: dispersion, repulsion, coulomb
    hipStream_t stream_bonded   = nullptr;  ///< Stream B: bonds, angles, dihedrals, inversions
    hipStream_t stream_threebody = nullptr; ///< Stream C: batm, atm, xbonds, hbonds
    hipEvent_t  event_pairwise  = nullptr;  ///< Sync event for stream A completion
    hipEvent_t  event_bonded    = nullptr;  ///< Sync event for stream B completion
    hipEvent_t  event_threebody = nullptr;  ///< Sync event for stream C completion
    hipEvent_t  event_upload    = nullptr;  ///< Reusable event for upload→stream sync
    // hipEvent_t  event_cn_ready  = nullptr;  // (disabled) was for overlap experiment: CN D2H before gaussian weights

    // Claude Generated (May 2026): Per-stream timing events (timing-enabled).
    // Used at verbosity >= 2 to populate FFTermTimings via hipEventElapsedTime.
    // 4 stream groups × 2 events. Recorded each step around the kernel block on each stream.
    hipEvent_t timing_sA_start = nullptr, timing_sA_end = nullptr;  ///< dispersion + repulsion
    hipEvent_t timing_sB_start = nullptr, timing_sB_end = nullptr;  ///< bonded (5 kernels)
    hipEvent_t timing_sC_start = nullptr, timing_sC_end = nullptr;  ///< 3-body (4 kernels)
    hipEvent_t timing_p2_start = nullptr, timing_p2_end = nullptr;  ///< Phase 2 (Coulomb + post)

    // === CUDA Graph Capture — Phase 8 (Claude Generated March 2026) ===
    // For topology-stable MD steps, the entire charge-independent kernel sequence is
    // captured once into m_graph_phase1 and replayed via hipGraphLaunch() — near-zero
    // CPU overhead instead of ~14 individual kernel-launch API calls (~70 µs/step).
    // Invalidated via invalidateGraph() whenever any SoA n-value changes (bonds, HB/XB).
    bool            m_graph_phase1_valid    = false;  ///< True when exec graph is ready to use
    hipGraph_t     m_graph_phase1          = nullptr; ///< Captured graph (template)
    hipGraphExec_t m_graph_exec_phase1     = nullptr; ///< Instantiated executable graph

    // === WP4 (May 2026): Phase-2 CUDA Graph Capture ===
    // launchChargeDependentAndFinish() is captured into one of two graph slots based on
    // gradient mode (energy-only vs gradient). Replay eliminates ~35 µs/step of kernel-
    // launch API overhead + ~10 µs of stream-event coordination on the charge-dependent
    // phase. Disabled when verbosity >= 3 (extra grad_snapshot writes) or per-kernel
    // diagnostics are active. Invalidated together with phase-1 graph on topology rebuild.
    bool            m_graph_phase2_grad_valid     = false;
    hipGraph_t     m_graph_phase2_grad           = nullptr;
    hipGraphExec_t m_graph_exec_phase2_grad      = nullptr;
    bool            m_graph_phase2_energy_valid   = false;
    hipGraph_t     m_graph_phase2_energy         = nullptr;
    hipGraphExec_t m_graph_exec_phase2_energy    = nullptr;
    /// WP4: Dedicated phase-2 events. Decouple phase-2 capture from phase-1's events.
    /// event_p2_pairwise: stream_pairwise → main stream sync (replaces event_pairwise inside Phase 2)
    /// event_p2_upload:   main stream → stream_pairwise upload fence (replaces event_upload inside Phase 2)
    /// Phase 2 still WAITS on event_bonded/event_threebody (recorded by Phase 1 graph), but never
    /// records them — that one-way dependency is safe under hipStreamCaptureModeThreadLocal.
    hipEvent_t     event_p2_pairwise             = nullptr;
    hipEvent_t     event_p2_upload               = nullptr;
};

// ============================================================================
// FFWorkspaceHip Constructor
// ============================================================================

FFWorkspaceHip::FFWorkspaceHip(const GFNFFParameterSet& params,
                                int natoms,
                                const std::vector<int>& atom_types)
    : m_natoms(natoms)
{
    // Verify CUDA device
    int dev_count = 0;
    checkCuda(hipGetDeviceCount(&dev_count), "hipGetDeviceCount");
    if (dev_count == 0)
        throw std::runtime_error("FFWorkspaceHip: no CUDA device found");

    // Pre-allocate pinned CPU staging buffers BEFORE any CUDA operations.
    // Claude Generated (March 2026): hipHostMalloc for async DMA and avoids
    // CUDA heap corruption issues with std::vector.
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_h_x), natoms * sizeof(double)),
              "pinned alloc m_h_x");
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_h_y), natoms * sizeof(double)),
              "pinned alloc m_h_y");
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_h_z), natoms * sizeof(double)),
              "pinned alloc m_h_z");
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_h_grad), 3 * natoms * sizeof(double)),
              "pinned alloc m_h_grad");
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_h_dEdcn_snap), natoms * sizeof(double)),
              "pinned alloc m_h_dEdcn_snap");
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_h_grad_snap), 3 * natoms * sizeof(double)),
              "pinned alloc m_h_grad_snap");
    // Per-kernel diagnostic pinned buffers (Claude Generated April 2026)
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_h_grad_after_sA), 3 * natoms * sizeof(double)),
              "pinned alloc m_h_grad_after_sA");
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_h_grad_after_bonds), 3 * natoms * sizeof(double)),
              "pinned alloc m_h_grad_after_bonds");
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_h_grad_after_angles), 3 * natoms * sizeof(double)),
              "pinned alloc m_h_grad_after_angles");
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_h_grad_after_sB), 3 * natoms * sizeof(double)),
              "pinned alloc m_h_grad_after_sB");
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_h_grad_after_sC), 3 * natoms * sizeof(double)),
              "pinned alloc m_h_grad_after_sC");
    // Claude Generated (March 2026): Pre-allocate pinned CN download buffer
    // NOTE: Do NOT allocate Eigen Vectors here — heap may be corrupted after CUDA allocs.
    // The caller (gfnff_gpu_method.cpp) pre-allocates m_gpu_cn_final before CUDA init.
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_h_cn_final), natoms * sizeof(double)),
              "pinned alloc m_h_cn_final");
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_h_cn_raw), natoms * sizeof(double)),
              "pinned alloc m_h_cn_raw");
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_h_dlogdcn), natoms * sizeof(double)),
              "pinned alloc m_h_dlogdcn");
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_h_energies),
                             FFWorkspaceHipImpl::N_ENERGY_SLOTS * sizeof(double)),
              "pinned alloc m_h_energies");

    m_impl = std::make_unique<FFWorkspaceHipImpl>();
    m_impl->N = natoms;

    // Create dedicated CUDA streams for multi-stream concurrency
    // Claude Generated (March 2026): Main stream + 3 concurrent streams for kernel groups
    checkCuda(hipStreamCreate(&m_impl->stream), "hipStreamCreate main");
    checkCuda(hipStreamCreate(&m_impl->stream_pairwise), "hipStreamCreate pairwise");
    checkCuda(hipStreamCreate(&m_impl->stream_bonded), "hipStreamCreate bonded");
    checkCuda(hipStreamCreate(&m_impl->stream_threebody), "hipStreamCreate threebody");
    checkCuda(hipEventCreateWithFlags(&m_impl->event_pairwise, hipEventDisableTiming), "hipEventCreate pairwise");
    checkCuda(hipEventCreateWithFlags(&m_impl->event_bonded, hipEventDisableTiming), "hipEventCreate bonded");
    checkCuda(hipEventCreateWithFlags(&m_impl->event_threebody, hipEventDisableTiming), "hipEventCreate threebody");
    checkCuda(hipEventCreateWithFlags(&m_impl->event_upload, hipEventDisableTiming), "hipEventCreate upload");
    checkCuda(hipEventCreateWithFlags(&m_impl->event_p2_pairwise, hipEventDisableTiming), "hipEventCreate p2_pairwise");
    checkCuda(hipEventCreateWithFlags(&m_impl->event_p2_upload, hipEventDisableTiming), "hipEventCreate p2_upload");
    // Claude Generated (May 2026): timing-enabled events for stream-level kernel timing
    checkCuda(hipEventCreate(&m_impl->timing_sA_start), "hipEventCreate timing_sA_start");
    checkCuda(hipEventCreate(&m_impl->timing_sA_end),   "hipEventCreate timing_sA_end");
    checkCuda(hipEventCreate(&m_impl->timing_sB_start), "hipEventCreate timing_sB_start");
    checkCuda(hipEventCreate(&m_impl->timing_sB_end),   "hipEventCreate timing_sB_end");
    checkCuda(hipEventCreate(&m_impl->timing_sC_start), "hipEventCreate timing_sC_start");
    checkCuda(hipEventCreate(&m_impl->timing_sC_end),   "hipEventCreate timing_sC_end");
    checkCuda(hipEventCreate(&m_impl->timing_p2_start), "hipEventCreate timing_p2_start");
    checkCuda(hipEventCreate(&m_impl->timing_p2_end),   "hipEventCreate timing_p2_end");
    hipStream_t stream = m_impl->stream;

    // Upload covalent radii to constant memory (used by angle/dihedral distance damping)
    upload_rcov_d3(s_rcov_d3_bohr, 87);

    // Upload covalent radii for HB/XB vdW radii lookup (Å units)
    upload_covalent_radii(GFNFFParameters::covalent_radii.data(),
                          static_cast<int>(GFNFFParameters::covalent_radii.size()));

    // Store dispersion pair indices on host for per-step dc6dcn extraction
    {
        const int nd = static_cast<int>(params.dispersions.size());
        m_disp_idx_i_host.resize(nd);
        m_disp_idx_j_host.resize(nd);
        m_h_dc6dcn_ij.resize(nd, 0.0);
        m_h_dc6dcn_ji.resize(nd, 0.0);
        for (int k = 0; k < nd; ++k) {
            m_disp_idx_i_host[k] = params.dispersions[k].i;
            m_disp_idx_j_host[k] = params.dispersions[k].j;
        }
    }

    // Store topology charges for BATM kernel; uploaded once below (static per topology)
    m_topology_charges = params.topology_charges;

    // ---------------------------------------------------------------------
    // DIAGNOSTIC A5 (plan phase A, mixture.xyz crash): validate that every
    // pair/triple/quadruple index is in [0, natoms). Out-of-bounds indices
    // would surface on the GPU as hipErrorInvalidValue/IllegalAddress on
    // the next kernel launch — exactly the symptom we are chasing on large
    // systems. This runs once at construction; cost is negligible.
    // Cleanup-target (plan phase C3): keep, or hide behind gpu_diagnostics.
    // ---------------------------------------------------------------------
    {
        const int N = natoms;
        int violations = 0;
        auto check2 = [&](int i, int j, const char* listName, int idx) {
            if (i < 0 || i >= N || j < 0 || j >= N) {
                if (violations < 16) {
                    CurcumaLogger::error(fmt::format(
                        "OOB index in {}[{}]: i={} j={} (N={})",
                        listName, idx, i, j, N));
                }
                ++violations;
            }
        };
        auto check3 = [&](int i, int j, int k, const char* listName, int idx) {
            if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N) {
                if (violations < 16) {
                    CurcumaLogger::error(fmt::format(
                        "OOB index in {}[{}]: i={} j={} k={} (N={})",
                        listName, idx, i, j, k, N));
                }
                ++violations;
            }
        };
        auto check4 = [&](int i, int j, int k, int l, const char* listName, int idx) {
            if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N || l < 0 || l >= N) {
                if (violations < 16) {
                    CurcumaLogger::error(fmt::format(
                        "OOB index in {}[{}]: i={} j={} k={} l={} (N={})",
                        listName, idx, i, j, k, l, N));
                }
                ++violations;
            }
        };
        for (size_t i = 0; i < params.dispersions.size(); ++i)
            check2(params.dispersions[i].i, params.dispersions[i].j, "dispersions", static_cast<int>(i));
        for (size_t i = 0; i < params.bonded_repulsions.size(); ++i)
            check2(params.bonded_repulsions[i].i, params.bonded_repulsions[i].j, "bonded_rep", static_cast<int>(i));
        for (size_t i = 0; i < params.nonbonded_repulsions.size(); ++i)
            check2(params.nonbonded_repulsions[i].i, params.nonbonded_repulsions[i].j, "nonbonded_rep", static_cast<int>(i));
        for (size_t i = 0; i < params.coulombs.size(); ++i)
            check2(params.coulombs[i].i, params.coulombs[i].j, "coulombs", static_cast<int>(i));
        for (size_t i = 0; i < params.bonds.size(); ++i)
            check2(params.bonds[i].i, params.bonds[i].j, "bonds", static_cast<int>(i));
        for (size_t i = 0; i < params.angles.size(); ++i)
            check3(params.angles[i].i, params.angles[i].j, params.angles[i].k, "angles", static_cast<int>(i));
        for (size_t i = 0; i < params.dihedrals.size(); ++i)
            check4(params.dihedrals[i].i, params.dihedrals[i].j,
                   params.dihedrals[i].k, params.dihedrals[i].l, "dihedrals", static_cast<int>(i));
        for (size_t i = 0; i < params.extra_dihedrals.size(); ++i)
            check4(params.extra_dihedrals[i].i, params.extra_dihedrals[i].j,
                   params.extra_dihedrals[i].k, params.extra_dihedrals[i].l, "extra_dihedrals", static_cast<int>(i));
        for (size_t i = 0; i < params.inversions.size(); ++i)
            check4(params.inversions[i].i, params.inversions[i].j,
                   params.inversions[i].k, params.inversions[i].l, "inversions", static_cast<int>(i));
        for (size_t i = 0; i < params.storsions.size(); ++i)
            check4(params.storsions[i].i, params.storsions[i].j,
                   params.storsions[i].k, params.storsions[i].l, "storsions", static_cast<int>(i));
        for (size_t i = 0; i < params.atm_triples.size(); ++i)
            check3(params.atm_triples[i].i, params.atm_triples[i].j,
                   params.atm_triples[i].k, "atm_triples", static_cast<int>(i));
        for (size_t i = 0; i < params.batm_triples.size(); ++i)
            check3(params.batm_triples[i].i, params.batm_triples[i].j,
                   params.batm_triples[i].k, "batm_triples", static_cast<int>(i));
        for (size_t i = 0; i < params.hbonds.size(); ++i)
            check3(params.hbonds[i].i, params.hbonds[i].j, params.hbonds[i].k, "hbonds", static_cast<int>(i));
        for (size_t i = 0; i < params.xbonds.size(); ++i)
            check3(params.xbonds[i].i, params.xbonds[i].j, params.xbonds[i].k, "xbonds", static_cast<int>(i));

        if (violations > 0) {
            CurcumaLogger::error(fmt::format(
                "FFWorkspaceHip param validation: {} out-of-bounds indices "
                "(N={}). First {} reported above; remainder suppressed.",
                violations, N, std::min(violations, 16)));
        } else {
            CurcumaLogger::info(fmt::format(
                "FFWorkspaceHip param validation: all pair/triple/quadruple "
                "indices in range (N={}, disp={}, hb={}, xb={}, atm={}, batm={})",
                N, params.dispersions.size(), params.hbonds.size(),
                params.xbonds.size(), params.atm_triples.size(), params.batm_triples.size()));
        }
    }

    // --- Upload static SoA interaction lists ---
    m_impl->disp.upload(params.dispersions, stream);
    m_impl->bonded_rep.upload(params.bonded_repulsions, stream);
    m_impl->nonbonded_rep.upload(params.nonbonded_repulsions, stream);
    m_impl->coulomb.upload(params.coulombs, natoms, stream);
    m_impl->bonds.upload(params.bonds, atom_types, stream);
    m_impl->angles.upload(params.angles, atom_types, stream);
    m_impl->dihedrals.upload(params.dihedrals, params.extra_dihedrals, atom_types, stream);
    m_impl->inversions.upload(params.inversions, atom_types, stream);

    // Phase 2 SoA uploads (March 2026): HB/XB/ATM/BATM/sTors — all on GPU
    m_impl->storsions.upload(params.storsions, stream);
    m_impl->batm.upload(params.batm_triples, stream);
    m_impl->atm.upload(params.atm_triples, atom_types, stream);
    m_impl->xbonds.upload(params.xbonds, atom_types, stream);
    m_last_xbonds = params.xbonds;   // init CPU mirror for debug comparison
    m_impl->hbonds.upload(params.hbonds, atom_types, stream);
    m_last_hbonds = params.hbonds;   // init CPU mirror for debug comparison

    // --- Allocate dynamic per-step buffers ---
    const int N3 = 3 * natoms;
    m_impl->coords.alloc(natoms);
    m_impl->d_charges.alloc(natoms);
    m_impl->d_cn.alloc(natoms);
    m_impl->d_topo_charges.alloc(natoms);
    // G-P1 (Apr 2026): Upload topology charges once at init — they are static per topology.
    // Removes the per-step blocking hipMemcpy from prepareAndLaunchChargeIndependent().
    if (m_topology_charges.size() > 0) {
        m_impl->d_topo_charges.upload(m_topology_charges.data(),
                                      static_cast<int>(m_topology_charges.size()), stream);
    }
    m_impl->d_grad.alloc(N3);
    m_impl->d_dEdcn.alloc(natoms);
    m_impl->d_energy.alloc(1);

    // Per-term energy accumulators — single contiguous buffer
    m_impl->d_energies.alloc(FFWorkspaceHipImpl::N_ENERGY_SLOTS);

    // dEdcn snapshot — always needed: dEdcnTotal() is public API (Coulomb chain-rule)
    m_impl->d_dEdcn_snapshot.alloc(natoms);
    // grad snapshot — always allocated to keep need_snapshots=true, which disables
    // Phase-1 CUDA graph capture. Graph capture fails on RTX 5080 (Blackwell): the
    // hipStreamBeginCapture path captures kernels without executing them; if
    // hipStreamEndCapture then fails, Phase-1 work is silently dropped → energy=0.
    // Keeping this allocated gates both Phase-1 and Phase-2 graph capture off.
    // Cost: ~32 KB (N3 doubles) GPU memory; acceptable until graph isolation is solved.
    m_impl->d_grad_snapshot.alloc(N3);
    // d_grad_after_* are debug-only snapshots. NOT allocated in normal mode:
    // keeping them null sets force_single_stream=false → enables 3-stream parallelism.
    // Allocate here only when gradient debugging is needed (guarded by m_verbosity if re-added).

    // --- GPU topology displacement check (Claude Generated March 2026) ---
    m_impl->d_disp_flag.alloc(1);

    // --- GPU CN computation buffers (Phase 1: GPU migration) ---
    // Claude Generated (March 2026): CN computed entirely on GPU
    m_impl->d_cn_raw.alloc(natoms);
    m_impl->d_cn_final.alloc(natoms);
    m_impl->d_atom_types.alloc(natoms);

    // WP5-C (May 2026): GPU-side D4 dc6dcn skip-check buffers
    m_impl->d_cn_d4_ref.alloc(natoms);
    // d_dc6dcn_block_max sized for max possible blocks (~8 for N=2000, blockSize=256)
    m_impl->d_dc6dcn_block_max.alloc(32);
    m_impl->d_dc6dcn_skip.alloc(1);
    hipHostAlloc(&m_impl->h_dc6dcn_skip, sizeof(int), hipHostMallocDefault);
    // Upload atom types once (static during simulation)
    std::vector<int> h_atom_types = atom_types;  // copy
    m_impl->d_atom_types.upload(h_atom_types, stream);

    // --- HB alpha chain-rule: allocate zz buffer + build (H,B) pair list ---
    m_impl->d_zz_hb.alloc(natoms);
    if (!params.bond_hb_data.empty()) {
        const auto& rcov_d3 = GFNFFParameters::covalent_rad_d3;  // already in Bohr
        constexpr double rcov_43 = 4.0 / 3.0;
        constexpr double rcov_scal = 1.78;

        std::vector<int> hb_h_idx, hb_b_idx;
        std::vector<double> hb_rcov;
        for (const auto& entry : params.bond_hb_data) {
            int H = entry.H;
            int ati = (H < natoms) ? atom_types[H] : 0;
            for (int B : entry.B_atoms) {
                int atj = (B < natoms) ? atom_types[B] : 0;
                if (ati < 1 || atj < 1) continue;
                double rcovij = rcov_scal * rcov_43 * (rcov_d3[ati - 1] + rcov_d3[atj - 1]);  // already in Bohr
                hb_h_idx.push_back(H);
                hb_b_idx.push_back(B);
                hb_rcov.push_back(rcovij);
            }
        }
        m_impl->hb_alpha.upload(hb_h_idx, hb_b_idx, hb_rcov, stream);
    }

    // --- Pre-allocate CPU result buffers ---
    m_result_gradient.resize(natoms, 3);
    m_result_gradient.setZero();
    m_dEdcn_total.resize(natoms);
    m_dEdcn_total.setZero();

    // --- Extract per-atom Coulomb self-energy params (for CPU postProcess TERM 2+3) ---
    // Mirrors FFWorkspace::setInteractionLists() lines 90-118
    if (!params.coulombs.empty()) {
        m_coul_chi_base  = Vector::Zero(natoms);
        m_coul_gam       = Vector::Zero(natoms);
        m_coul_alp       = Vector::Zero(natoms);
        m_coul_chi_static= Vector::Zero(natoms);

        std::vector<bool> seen(natoms, false);
        for (const auto& c : params.coulombs) {
            if (c.i < natoms && !seen[c.i]) {
                m_coul_chi_base(c.i)   = c.chi_base_i;
                m_coul_gam(c.i)        = c.gam_i;
                m_coul_alp(c.i)        = c.alp_i;
                m_coul_chi_static(c.i) = c.chi_i;
                seen[c.i] = true;
            }
            if (c.j < natoms && !seen[c.j]) {
                m_coul_chi_base(c.j)   = c.chi_base_j;
                m_coul_gam(c.j)        = c.gam_j;
                m_coul_alp(c.j)        = c.alp_j;
                m_coul_chi_static(c.j) = c.chi_j;
                seen[c.j] = true;
            }
        }
        // Upload Coulomb self-energy parameters to GPU for k_coulomb_self kernel
        // WP5-B (May 2026): CNF omitted — eeq_topo.d_cnf (topology-constant) used instead.
        setCoulombSelfEnergyParams(m_coul_chi_base, m_coul_gam, m_coul_alp,
                                   m_coul_chi_static);
    }

    // Store initial charges and E0 from parameter set
    if (params.eeq_charges.size() == natoms)
        m_eeq_charges = params.eeq_charges;
    m_e0 = params.e0;

    // Synchronise: wait for all uploads to finish
    checkCuda(hipStreamSynchronize(stream), "init sync");

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format(
            "FFWorkspaceHip: {} atoms | disp={} brep={} nbrep={} coul={} bond={} ang={} dih={} inv={}"
            " stors={} batm={} atm={} xb={} hb={}",
            natoms,
            m_impl->disp.n, m_impl->bonded_rep.n, m_impl->nonbonded_rep.n,
            m_impl->coulomb.n, m_impl->bonds.n, m_impl->angles.n,
            m_impl->dihedrals.n, m_impl->inversions.n,
            m_impl->storsions.n, m_impl->batm.n, m_impl->atm.n,
            m_impl->xbonds.n, m_impl->hbonds.n));
    }
}

// ============================================================================
// Destructor
// ============================================================================

FFWorkspaceHip::~FFWorkspaceHip()
{
    if (m_impl) {
        // Synchronize all streams before destroying
        if (m_impl->stream)          hipStreamSynchronize(m_impl->stream);
        if (m_impl->stream_pairwise) hipStreamSynchronize(m_impl->stream_pairwise);
        if (m_impl->stream_bonded)   hipStreamSynchronize(m_impl->stream_bonded);
        if (m_impl->stream_threebody)hipStreamSynchronize(m_impl->stream_threebody);

        if (m_impl->event_pairwise)    hipEventDestroy(m_impl->event_pairwise);
        if (m_impl->event_bonded)      hipEventDestroy(m_impl->event_bonded);
        if (m_impl->event_threebody)   hipEventDestroy(m_impl->event_threebody);
        if (m_impl->event_upload)      hipEventDestroy(m_impl->event_upload);
        if (m_impl->event_p2_pairwise) hipEventDestroy(m_impl->event_p2_pairwise);
        if (m_impl->event_p2_upload)   hipEventDestroy(m_impl->event_p2_upload);

        if (m_impl->stream)          hipStreamDestroy(m_impl->stream);
        if (m_impl->stream_pairwise) hipStreamDestroy(m_impl->stream_pairwise);
        if (m_impl->stream_bonded)   hipStreamDestroy(m_impl->stream_bonded);
        if (m_impl->stream_threebody)hipStreamDestroy(m_impl->stream_threebody);

        // Phase 8: destroy phase-1 graph objects
        if (m_impl->m_graph_exec_phase1) hipGraphExecDestroy(m_impl->m_graph_exec_phase1);
        if (m_impl->m_graph_phase1)      hipGraphDestroy(m_impl->m_graph_phase1);
        // WP4: destroy phase-2 graph objects (gradient + energy slots)
        if (m_impl->m_graph_exec_phase2_grad)   hipGraphExecDestroy(m_impl->m_graph_exec_phase2_grad);
        if (m_impl->m_graph_phase2_grad)        hipGraphDestroy(m_impl->m_graph_phase2_grad);
        if (m_impl->m_graph_exec_phase2_energy) hipGraphExecDestroy(m_impl->m_graph_exec_phase2_energy);
        if (m_impl->m_graph_phase2_energy)      hipGraphDestroy(m_impl->m_graph_phase2_energy);
    }

    // Free pinned memory staging buffers
    if (m_h_x)          hipHostFree(m_h_x);
    if (m_h_y)          hipHostFree(m_h_y);
    if (m_h_z)          hipHostFree(m_h_z);
    if (m_h_grad)       hipHostFree(m_h_grad);
    if (m_h_dEdcn_snap) hipHostFree(m_h_dEdcn_snap);
    if (m_h_grad_snap)  hipHostFree(m_h_grad_snap);
    if (m_h_grad_after_sA)     hipHostFree(m_h_grad_after_sA);
    if (m_h_grad_after_bonds)  hipHostFree(m_h_grad_after_bonds);
    if (m_h_grad_after_angles) hipHostFree(m_h_grad_after_angles);
    if (m_h_grad_after_sB)     hipHostFree(m_h_grad_after_sB);
    if (m_h_grad_after_sC)     hipHostFree(m_h_grad_after_sC);
    if (m_h_cn_final)   hipHostFree(m_h_cn_final);
    if (m_h_cn_raw)     hipHostFree(m_h_cn_raw);
    if (m_h_dlogdcn)    hipHostFree(m_h_dlogdcn);
    if (m_h_energies)   hipHostFree(m_h_energies);
    if (m_impl && m_impl->h_dc6dcn_skip) hipHostFree(m_impl->h_dc6dcn_skip);
}

// ============================================================================
// Phase 8: CUDA Graph invalidation (Claude Generated March 2026)
// ============================================================================

void FFWorkspaceHip::invalidateGraph()
{
    auto& impl = *m_impl;
    if (impl.m_graph_exec_phase1) {
        hipGraphExecDestroy(impl.m_graph_exec_phase1);
        impl.m_graph_exec_phase1 = nullptr;
    }
    if (impl.m_graph_phase1) {
        hipGraphDestroy(impl.m_graph_phase1);
        impl.m_graph_phase1 = nullptr;
    }
    impl.m_graph_phase1_valid = false;

    // WP4: invalidate phase-2 gradient slot
    if (impl.m_graph_exec_phase2_grad) {
        hipGraphExecDestroy(impl.m_graph_exec_phase2_grad);
        impl.m_graph_exec_phase2_grad = nullptr;
    }
    if (impl.m_graph_phase2_grad) {
        hipGraphDestroy(impl.m_graph_phase2_grad);
        impl.m_graph_phase2_grad = nullptr;
    }
    impl.m_graph_phase2_grad_valid = false;

    // WP4: invalidate phase-2 energy slot
    if (impl.m_graph_exec_phase2_energy) {
        hipGraphExecDestroy(impl.m_graph_exec_phase2_energy);
        impl.m_graph_exec_phase2_energy = nullptr;
    }
    if (impl.m_graph_phase2_energy) {
        hipGraphDestroy(impl.m_graph_phase2_energy);
        impl.m_graph_phase2_energy = nullptr;
    }
    impl.m_graph_phase2_energy_valid = false;
}

// ============================================================================
// Per-step state setters
// ============================================================================

void FFWorkspaceHip::setEEQCharges(const Vector& q)
{
    m_eeq_charges = q;
}

void FFWorkspaceHip::setEEQDeviceCharges(const double* d_src)
{
    // WP5-A: d_src is d_rhs[0..N-1] from EEQ solver (EEQ stream fully synced by caller).
    // Enqueue D2D copy on main workspace stream — no CPU involvement.
    hipMemcpyAsync(m_impl->d_charges.ptr, d_src, m_natoms * sizeof(double),
                    hipMemcpyDeviceToDevice, m_impl->stream);
    m_device_charges_ready = true;
}

void FFWorkspaceHip::setTopologyCharges(const Vector& q)
{
    m_topology_charges = q;
}

void FFWorkspaceHip::setD3CN(const Vector& cn)
{
    m_cn = cn;
}

void FFWorkspaceHip::setCNDerivatives(const Vector& /* cn */, const Vector& /* cnf */,
                                       const std::vector<SpMatrix>& /* dcn */)
{
    // WP5-B (May 2026): No-op. CNF for the Coulomb postprocess (TERM 1b / k_subtract_qtmp)
    // already lives on GPU as the topology-constant `impl.eeq_topo.d_cnf` buffer from
    // uploadEEQTopologyParams(). The previous per-step copy into m_cnf was only used as a
    // size probe for the do_subtract guard, which now checks impl.coul_self_on_gpu directly.
}

void FFWorkspaceHip::setE0(double e0)
{
    m_e0 = e0;
}

void FFWorkspaceHip::updateBondHBCN(const std::vector<double>& hb_cn_values)
{
    const int nb = m_impl->bonds.n;
    if (nb == 0 || static_cast<int>(hb_cn_values.size()) != nb) return;
    m_impl->bonds.hb_cn_H.upload(hb_cn_values.data(), nb);
}

// Claude Generated (Apr 2026): Re-upload HB alpha (H,B) pair list after dynamic HB re-detection.
// Mirrors the construction-time logic at ff_workspace_gpu.cu:848-870.
void FFWorkspaceHip::updateHBAlphaPairs(const std::vector<BondHBEntry>& bond_hb_data,
                                          const std::vector<int>& atom_types)
{
    if (!m_impl) return;
    const auto& rcov_d3 = GFNFFParameters::covalent_rad_d3;
    constexpr double rcov_43 = 4.0 / 3.0;
    constexpr double rcov_scal = 1.78;

    std::vector<int> hb_h_idx, hb_b_idx;
    std::vector<double> hb_rcov;
    for (const auto& entry : bond_hb_data) {
        int H = entry.H;
        int ati = (H < m_natoms) ? atom_types[H] : 0;
        for (int B : entry.B_atoms) {
            int atj = (B < m_natoms) ? atom_types[B] : 0;
            if (ati < 1 || atj < 1) continue;
            double rcovij = rcov_scal * rcov_43 * (rcov_d3[ati - 1] + rcov_d3[atj - 1]);
            hb_h_idx.push_back(H);
            hb_b_idx.push_back(B);
            hb_rcov.push_back(rcovij);
        }
    }
    m_impl->hb_alpha.upload(hb_h_idx, hb_b_idx, hb_rcov, m_impl->stream);
}

// Claude Generated (Apr 2026): Re-upload per-bond nr_hb and hb_H_atom to GPU BondSoA.
// Called after updateHBAlphaPairs() to keep k_bonds kernel HB modulation current.
void FFWorkspaceHip::updateBondHBMetadata(const std::vector<int>& nr_hb,
                                            const std::vector<int>& hb_H_atom)
{
    if (!m_impl) return;
    const int nb = m_impl->bonds.n;
    if (nb == 0) return;
    if (static_cast<int>(nr_hb.size()) != nb || static_cast<int>(hb_H_atom.size()) != nb)
        return;
    m_impl->bonds.nr_hb.upload(nr_hb.data(), nb, m_impl->stream);
    m_impl->bonds.hb_H_atom.upload(hb_H_atom.data(), nb, m_impl->stream);
}

// ============================================================================
// GPU HB CN per-bond computation (Apr 2026)
// Replaces CPU HB CN loop in gfnff_gpu_method.cpp.
// ============================================================================
void FFWorkspaceHip::computeBondHBCN()
{
    auto& impl = *m_impl;
    const int N = m_natoms;
    if (N <= 0 || impl.hb_alpha.n == 0) return;

    // Allocate/zero per-H accumulator
    if (!impl.d_hb_cn_per_atom.ptr || impl.d_hb_cn_per_atom.n < N)
        impl.d_hb_cn_per_atom.alloc(N);
    hipMemsetAsync(impl.d_hb_cn_per_atom.ptr, 0, N * sizeof(double), impl.stream);

    // Pass 1: accumulate HB CN per H-atom
    auto cfg1 = getLaunchConfig(impl.hb_alpha.n, m_block_size);
    k_hb_cn_per_atom<<<cfg1.gridSize, cfg1.blockSize, 0, impl.stream>>>(
        impl.hb_alpha.n,
        impl.hb_alpha.idx_H.ptr,
        impl.hb_alpha.idx_B.ptr,
        impl.hb_alpha.rcov_sum.ptr,
        impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
        impl.d_hb_cn_per_atom.ptr,
        27.5,   // kn
        900.0); // thr_sq

    // Pass 2: map per-H CN to per-bond CN
    auto cfg2 = getLaunchConfig(impl.bonds.n, m_block_size);
    k_hb_cn_map_to_bonds<<<cfg2.gridSize, cfg2.blockSize, 0, impl.stream>>>(
        impl.bonds.n,
        impl.bonds.nr_hb.ptr,
        impl.bonds.hb_H_atom.ptr,
        impl.d_hb_cn_per_atom.ptr,
        impl.bonds.hb_cn_H.ptr);

    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::info("  GPU HB CN per-bond computed on GPU");
}

void FFWorkspaceHip::setCoulombSelfEnergyParams(const Vector& chi_base, const Vector& gam,
                                                  const Vector& alp,
                                                  const Vector& chi_static)
{
    m_coul_chi_base   = chi_base;
    m_coul_gam        = gam;
    m_coul_alp        = alp;
    m_coul_chi_static = chi_static;

    // Upload Coulomb self-energy parameters to GPU
    // Claude Generated (March 2026): GPU-side Coulomb TERM 2+3 + qtmp
    // WP5-B (May 2026): CNF no longer uploaded here — topology-constant eeq_topo.d_cnf
    // (from uploadEEQTopologyParams) is used directly by k_coulomb_postprocess.
    const int N = m_natoms;
    if (chi_base.size() == N && gam.size() == N && alp.size() == N) {
        auto& impl = *m_impl;
        impl.d_coul_chi_base.alloc(N);
        impl.d_coul_gam.alloc(N);
        impl.d_coul_alp.alloc(N);
        impl.d_coul_chi_base.upload(chi_base.data(), N);
        impl.d_coul_gam.upload(gam.data(), N);
        impl.d_coul_alp.upload(alp.data(), N);
        impl.coul_self_on_gpu = true;
    }
}

void FFWorkspaceHip::updateHBonds(const std::vector<GFNFFHydrogenBond>& hbonds,
                                    const std::vector<int>& atom_types)
{
    // Re-upload HBond SoA after dynamic re-detection.
    // Claude Generated (March 2026): Needed because updateHBXBIfNeeded() may find
    // additional HBonds on re-detection that weren't in the initial parameter set.
    m_impl->hbonds.upload(hbonds, atom_types, m_impl->stream);
    // Store for CPU vs GPU comparison debugging
    m_last_hbonds = hbonds;
}

void FFWorkspaceHip::updateXBonds(const std::vector<GFNFFHalogenBond>& xbonds,
                                    const std::vector<int>& atom_types)
{
    m_impl->xbonds.upload(xbonds, atom_types, m_impl->stream);
    // Store for CPU vs GPU comparison debugging
    m_last_xbonds = xbonds;
}

void FFWorkspaceHip::setDispersionEnabled(bool v)  { m_dispersion_enabled = v; }
void FFWorkspaceHip::setHBondEnabled(bool v)        { m_hbond_enabled = v; }
void FFWorkspaceHip::setRepulsionEnabled(bool v)    { m_repulsion_enabled = v; }
void FFWorkspaceHip::setCoulombEnabled(bool v)      { m_coulomb_enabled = v; }

/// Claude Generated (April 2026): Upload unit cell to GPU constant memory for PBC.
/// cell_bohr_9 and cell_bohr_inv_9 are column-major 3×3 matrices (9 doubles each).
/// Zero overhead when has_pbc==false because applyMIC() returns immediately.
void FFWorkspaceHip::setUnitCell(const double* cell_bohr_9, const double* cell_bohr_inv_9, bool has_pbc)
{
    hipMemcpyToSymbol(HIP_SYMBOL(d_lattice),     cell_bohr_9,     9 * sizeof(double));
    hipMemcpyToSymbol(HIP_SYMBOL(d_lattice_inv), cell_bohr_inv_9, 9 * sizeof(double));
    int flag = has_pbc ? 1 : 0;
    hipMemcpyToSymbol(HIP_SYMBOL(d_has_pbc), &flag, sizeof(int));
}

// ============================================================================
// GPU-only CN chain-rule setters (Claude Generated March 2026)
// ============================================================================

void FFWorkspaceHip::setCNPairList(const std::vector<int>& idx_i,
                                    const std::vector<int>& idx_j,
                                    const std::vector<double>& rcov_sum)
{
    const int np = static_cast<int>(idx_i.size());
    if (np == 0 || idx_j.size() != idx_i.size() || rcov_sum.size() != idx_i.size()) return;

    m_cn_pair_i    = idx_i;
    m_cn_pair_j    = idx_j;
    m_cn_pair_rcov = rcov_sum;
    m_cn_n_pairs   = np;

    // Upload to GPU
    auto& impl = *m_impl;
    impl.d_cn_idx_i.alloc(np);
    impl.d_cn_idx_j.alloc(np);
    impl.d_cn_rcov_sum.alloc(np);
    impl.d_cn_idx_i.upload(idx_i.data(), np);
    impl.d_cn_idx_j.upload(idx_j.data(), np);
    impl.d_cn_rcov_sum.upload(rcov_sum.data(), np);
    impl.n_cn_pairs = np;
    m_cn_pairs_on_gpu = true;

    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::info(fmt::format("  GPU CN chain-rule: {} pairs uploaded", np));
}

void FFWorkspaceHip::setDlogDCN(const Vector& dlogdcn)
{
    m_dlogdcn = dlogdcn;
    // Upload to GPU
    auto& impl = *m_impl;
    const int N = m_natoms;
    if (dlogdcn.size() != N) return;
    if (!impl.d_dlogdcn.ptr) impl.d_dlogdcn.alloc(N);
    impl.d_dlogdcn.upload(dlogdcn.data(), N);
}

// ============================================================================
// GPU CN pair list generation (Apr 2026)
// Replaces CPU generateCNPairList() O(N^2) loop with two-pass GPU kernels.
// ============================================================================
void FFWorkspaceHip::generateCNPairListOnGPU()
{
    auto& impl = *m_impl;
    const int N = m_natoms;
    if (N <= 0) return;

    // Ensure counter buffer is allocated
    if (!impl.d_cn_pair_counter.ptr) impl.d_cn_pair_counter.alloc(1);

    // Pass 1: count valid pairs
    int h_count = 0;
    impl.d_cn_pair_counter.upload(&h_count, 1);  // zero the counter

    int n_total = N * (N - 1) / 2;
    auto cfg = getLaunchConfig(n_total, m_block_size);
    k_generate_cn_pairs_count<<<cfg.gridSize, cfg.blockSize, 0, impl.stream>>>(
        N,
        impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
        impl.d_atom_types.ptr,
        m_cn_pair_cutoff_factor,  // cutoff_factor (Task #10: tunable)
        impl.d_cn_pair_counter.ptr);
    checkCuda(hipStreamSynchronize(impl.stream), "k_generate_cn_pairs_count sync");

    // Download count
    impl.d_cn_pair_counter.download(&h_count, 1);
    checkCuda(hipStreamSynchronize(impl.stream), "download pair count");

    if (h_count == 0) {
        impl.n_cn_pairs = 0;
        m_cn_pairs_on_gpu = false;
        if (CurcumaLogger::get_verbosity() >= 3)
            CurcumaLogger::info("  GPU CN pair list: 0 pairs");
        return;
    }

    // Allocate pair buffers based on count
    impl.d_cn_idx_i.alloc(h_count);
    impl.d_cn_idx_j.alloc(h_count);
    impl.d_cn_rcov_sum.alloc(h_count);
    impl.n_cn_pairs = h_count;

    // Pass 2: write valid pairs
    int h_zero = 0;
    impl.d_cn_pair_counter.upload(&h_zero, 1);  // reset counter for write indices

    k_generate_cn_pairs_write<<<cfg.gridSize, cfg.blockSize, 0, impl.stream>>>(
        N,
        impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
        impl.d_atom_types.ptr,
        m_cn_pair_cutoff_factor,  // cutoff_factor (Task #10: tunable)
        impl.d_cn_pair_counter.ptr,
        impl.d_cn_idx_i.ptr,
        impl.d_cn_idx_j.ptr,
        impl.d_cn_rcov_sum.ptr);
    checkCuda(hipStreamSynchronize(impl.stream), "k_generate_cn_pairs_write sync");

    m_cn_pairs_on_gpu = true;

    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::info(fmt::format("  GPU CN pair list: {} pairs generated on GPU", h_count));
}

void FFWorkspaceHip::updateDispersionDC6DCN(const Matrix& dc6dcn)
{
    // Claude Generated (March 2026): Extract per-pair dc6dcn values and upload to GPU.
    // Eliminates the CPU round-trip (download dEdcn, add dispersion, re-upload).
    const int nd = static_cast<int>(m_disp_idx_i_host.size());
    if (nd == 0 || dc6dcn.size() == 0) return;

    for (int k = 0; k < nd; ++k) {
        int i = m_disp_idx_i_host[k];
        int j = m_disp_idx_j_host[k];
        if (i < dc6dcn.rows() && j < dc6dcn.cols()) {
            m_h_dc6dcn_ij[k] = dc6dcn(i, j);
            m_h_dc6dcn_ji[k] = dc6dcn(j, i);
        } else {
            m_h_dc6dcn_ij[k] = 0.0;
            m_h_dc6dcn_ji[k] = 0.0;
        }
    }

    m_impl->disp.dc6dcn_ij.upload(m_h_dc6dcn_ij.data(), nd);
    m_impl->disp.dc6dcn_ji.upload(m_h_dc6dcn_ji.data(), nd);
}

// ============================================================================
// Phase 2: GPU dc6dcn per-pair computation (Claude Generated March 2026)
// ============================================================================

void FFWorkspaceHip::uploadC6ReferenceTable(const std::vector<double>& c6_flat,
                                              const std::vector<int>& refn)
{
    auto& impl = *m_impl;

    // Upload C6 reference table (one-time, ~3.1 MB)
    impl.d_c6_flat.upload(c6_flat.data(), static_cast<int>(c6_flat.size()));

    // Upload refn array (118 ints) to constant memory (Phase 8: broadcast via L1)
    upload_refn_const(refn.data(), static_cast<int>(refn.size()));

    // Pre-allocate gw/dgw GPU buffers
    const int gw_size = m_natoms * D4_MAX_REF;
    impl.d_gw.alloc(gw_size);
    impl.d_dgw.alloc(gw_size);

    // Pre-allocate dc6dcn output buffers (written by k_dc6dcn_per_pair kernel)
    const int nd = impl.disp.n;
    if (nd > 0) {
        impl.disp.dc6dcn_ij.alloc(nd);
        impl.disp.dc6dcn_ji.alloc(nd);
    }

    // Pre-allocate host staging buffers
    m_h_gw_flat.resize(gw_size, 0.0);
    m_h_dgw_flat.resize(gw_size, 0.0);

    impl.dc6dcn_gpu_ready = true;
}

void FFWorkspaceHip::computeDC6DCNOnGPU(const std::vector<std::vector<double>>& gw,
                                          const std::vector<std::vector<double>>& dgw)
{
    auto& impl = *m_impl;
    if (!impl.dc6dcn_gpu_ready) return;

    const int N = m_natoms;
    const int nd = impl.disp.n;
    if (nd == 0) return;

    // Flatten nested weight arrays to [N × MAX_REF] with zero-padding
    for (int i = 0; i < N; ++i) {
        int base = i * D4_MAX_REF;
        int nref = static_cast<int>(gw[i].size());
        for (int r = 0; r < D4_MAX_REF; ++r) {
            m_h_gw_flat[base + r]  = (r < nref) ? gw[i][r] : 0.0;
            m_h_dgw_flat[base + r] = (r < nref) ? dgw[i][r] : 0.0;
        }
    }

    // Upload to GPU
    const int gw_size = N * D4_MAX_REF;
    impl.d_gw.upload(m_h_gw_flat.data(), gw_size);
    impl.d_dgw.upload(m_h_dgw_flat.data(), gw_size);

    // Launch dc6dcn per-pair kernel (Phase 6: dynamic block size)
    LaunchConfig cfg = getLaunchConfig(nd, m_block_size);
    // refn read from constant memory d_refn_const (Phase 8)
    k_dc6dcn_per_pair<<<cfg.gridSize, cfg.blockSize, 0, impl.stream>>>(
        nd,
        impl.disp.idx_i.ptr,
        impl.disp.idx_j.ptr,
        impl.d_atom_types.ptr,
        impl.d_gw.ptr,
        impl.d_dgw.ptr,
        impl.d_c6_flat.ptr,
        impl.disp.dc6dcn_ij.ptr,
        impl.disp.dc6dcn_ji.ptr
    );

    // No sync needed here: k_dc6dcn_per_pair runs on impl.stream (main stream).
    // launchChargeDependentAndFinish() does hipStreamWaitEvent(stream_pairwise, stream)
    // before launching k_dispersion, which guarantees dc6dcn is ready.
    // Phase 4a: Removed unnecessary hipStreamSynchronize (Claude Generated March 2026)
}

// ============================================================================
// Phase 6: Upload reference CN values (one-time at init)
// Claude Generated (March 2026)
// ============================================================================

void FFWorkspaceHip::uploadRefCN(const std::vector<std::vector<double>>& refcn)
{
    auto& impl = *m_impl;

    // Flatten nested refcn[MAX_ELEM][var_len] → flat [MAX_ELEM × D4_MAX_REF]
    const int flat_size = D4_MAX_ELEM * D4_MAX_REF;
    std::vector<double> refcn_flat(flat_size, 0.0);

    for (int elem = 0; elem < D4_MAX_ELEM && elem < static_cast<int>(refcn.size()); ++elem) {
        int base = elem * D4_MAX_REF;
        int nref = static_cast<int>(refcn[elem].size());
        for (int r = 0; r < D4_MAX_REF && r < nref; ++r) {
            refcn_flat[base + r] = refcn[elem][r];
        }
    }

    // Upload to constant memory (Phase 8: L1-broadcast instead of global memory load)
    upload_refcn_const(refcn_flat.data(), flat_size);
}

// ============================================================================
// Phase 6: Compute Gaussian weights + dc6dcn entirely on GPU
// Claude Generated (March 2026): Eliminates CPU precomputeGaussianWeights()
// + computeGaussianWeightDerivatives() + flatten + sync H2D upload.
// ============================================================================

void FFWorkspaceHip::computeGaussianWeightsOnGPU(bool use_cn_final)
{
    auto& impl = *m_impl;
    if (!impl.dc6dcn_gpu_ready) return;

    const int N = m_natoms;
    const int nd = impl.disp.n;
    if (nd == 0) return;

    // Select CN source buffer:
    //   use_cn_final=true  → d_cn_final (written by computeCN(), current step).
    //                         Call BEFORE prepareAndLaunchChargeIndependent() so
    //                         dc6dcn is ready on the main stream before the event_upload
    //                         fence that sA (dispersion stream) waits on.
    //   use_cn_final=false → d_cn (updated by D2D copy inside prepareAndLaunchChargeIndependent).
    //                         Call AFTER prepareAndLaunchChargeIndependent() — legacy path.
    double* cn_buf = (use_cn_final && impl.d_cn_final.ptr) ? impl.d_cn_final.ptr : impl.d_cn.ptr;

    LaunchConfig cfg_gw = getLaunchConfig(N, m_block_size);
    k_gaussian_weights<<<cfg_gw.gridSize, cfg_gw.blockSize, 0, impl.stream>>>(
        N,
        cn_buf,
        impl.d_atom_types.ptr,
        impl.d_gw.ptr,
        impl.d_dgw.ptr);

    LaunchConfig cfg_dc6 = getLaunchConfig(nd, m_block_size);
    k_dc6dcn_per_pair<<<cfg_dc6.gridSize, cfg_dc6.blockSize, 0, impl.stream>>>(
        nd,
        impl.disp.idx_i.ptr,
        impl.disp.idx_j.ptr,
        impl.d_atom_types.ptr,
        impl.d_gw.ptr,
        impl.d_dgw.ptr,
        impl.d_c6_flat.ptr,
        impl.disp.dc6dcn_ij.ptr,
        impl.disp.dc6dcn_ji.ptr);
    // No explicit sync needed: stream ordering + event_upload fence in
    // prepareAndLaunchChargeIndependent guarantees dc6dcn is ready before k_dispersion on sA.
}

// ============================================================================
// Phase 3: Split calculation for CPU/GPU overlap
// Claude Generated (March 2026)
//
// prepareAndLaunchChargeIndependent() — zeros accumulators, uploads CN+geometry+
//   topology charges, launches ALL charge-independent kernels (everything except
//   k_coulomb, k_coulomb_self, k_subtract_qtmp). Non-blocking.
//
// launchChargeDependentAndFinish() — uploads EEQ charges, launches Coulomb +
//   postprocess kernels, downloads results. Blocking.
//
// The original calculate() is preserved for backward compatibility and calls
// both methods in sequence.
// ============================================================================

void FFWorkspaceHip::prepareAndLaunchChargeIndependent(bool gradient)
{
    auto& impl = *m_impl;
    const int N = m_natoms;
    hipStream_t stream = impl.stream;

    // =========================================================================
    // Phase 8: CUDA Graph fast path (Claude Generated March 2026)
    // For topology-stable MD steps, the entire charge-independent kernel sequence
    // was captured on the first call and is replayed here via hipGraphLaunch().
    // CPU-GPU overlap with EEQ is preserved: hipGraphLaunch() is asynchronous.
    // =========================================================================
    // need_snapshots=true when both snapshot buffers are allocated (always true in practice).
    // d_grad_snapshot is always allocated to keep need_snapshots=true, disabling graph capture.
    // Phase-1 graph capture fails on RTX 5080: captured kernels are silently dropped when
    // hipStreamEndCapture fails (stream isolation due to pending non-captured work on worker
    // streams). Keep graph capture disabled until root cause is resolved on Blackwell.
    const bool need_snapshots = impl.d_dEdcn_snapshot.ptr && impl.d_grad_snapshot.ptr;
    const bool force_single_stream = (gradient && impl.d_grad_after_bonds.ptr);
    // Claude Generated (May 2026): m_record_kernel_timings disables graph replay so
    // per-stream events get recorded each step (graph-replay events don't update times).
    if (impl.m_graph_phase1_valid && !need_snapshots && !force_single_stream
        && !m_record_kernel_timings) {
        hipGraphLaunch(impl.m_graph_exec_phase1, stream);
        return;  // GPU executes the graph asynchronously; CPU proceeds to EEQ
    }

    // === CAPTURE START: on the first topology-stable step, record this execution ===
    bool capturing = false;
    if (!impl.m_graph_phase1_valid && !need_snapshots && !force_single_stream
        && !m_record_kernel_timings) {
        hipError_t cap_err = hipStreamBeginCapture(stream, hipStreamCaptureModeThreadLocal);
        capturing = (cap_err == hipSuccess);
    }

    // =========================================================================
    // 1. Zero GPU accumulators (total + per-term)
    // =========================================================================
    hipMemsetAsync(impl.d_energy.ptr, 0, sizeof(double),         stream);
    hipMemsetAsync(impl.d_grad.ptr,   0, 3*N*sizeof(double),     stream);
    hipMemsetAsync(impl.d_dEdcn.ptr,  0, N*sizeof(double),       stream);
    hipMemsetAsync(impl.d_energies.ptr, 0,
                    FFWorkspaceHipImpl::N_ENERGY_SLOTS * sizeof(double), stream);

    // =========================================================================
    // 2. Upload CN (NOT EEQ charges — those come later)
    //    Geometry was already uploaded by setGeometry() before this call.
    //    Topology charges (d_topo_charges) uploaded once at setParameter() — static per topology.
    // =========================================================================
    // G-P1 (Apr 2026): Use GPU-side D2D copy d_cn_final → d_cn instead of CPU round-trip.
    // d_cn_final was written by k_cn_compute on the same (main) stream, so it is ready here.
    // D2D copy is async on main stream; event_upload recorded below makes worker streams wait.
    if (impl.d_cn_final.ptr && impl.d_cn.ptr && impl.d_cn_final.n >= N) {
        hipMemcpyAsync(impl.d_cn.ptr, impl.d_cn_final.ptr,
                        N * sizeof(double), hipMemcpyDeviceToDevice, stream);
    } else if (static_cast<int>(m_cn.size()) == N) {
        impl.d_cn.upload(m_cn.data(), N, stream);  // fallback (before first CN compute)
    }

    // WP2: Build EEQ RHS on GPU — rhs[i] = chi_corr[i] + cnf[i]*sqrt(max(cn[i],0))
    // Runs on main stream directly after d_cn_final is available; no CPU sync needed.
    // finalizeCNForCPU() (called later) syncs main stream, so d_rhs_atoms is ready
    // before solveWithDeviceRHS() starts its D2D copy.
    if (impl.eeq_topo.valid && impl.d_cn_final.ptr) {
        // BOUNDS_LIGHT kernel — cap block at 256.
        LaunchConfig cfg = getLaunchConfigLight(N, m_block_size);
        k_build_eeq_rhs<<<cfg.gridSize, cfg.blockSize, 0, stream>>>(
            N,
            impl.d_cn_final.ptr,
            impl.eeq_topo.d_chi_corrected.ptr,
            impl.eeq_topo.d_cnf.ptr,
            impl.eeq_topo.d_rhs_atoms.ptr);
    }

    // =========================================================================
    // 3. Launch charge-independent kernels on parallel streams
    //    These kernels do NOT read d_charges (EEQ). They use only coords, CN,
    //    topology charges, and static SoA parameters.
    // =========================================================================
    // Per-kernel diagnostic mode: force all kernels onto main stream so
    // snapshots reflect exactly one kernel's contribution at a time.
    // Normal mode: 3 parallel streams for maximum GPU occupancy.
    hipStream_t sA = force_single_stream ? stream : impl.stream_pairwise;
    hipStream_t sB = force_single_stream ? stream : impl.stream_bonded;
    hipStream_t sC = force_single_stream ? stream : impl.stream_threebody;

    // Ensure all uploads on main stream are visible to worker streams
    if (!force_single_stream) {
        hipEventRecord(impl.event_upload, stream);
        hipStreamWaitEvent(sA, impl.event_upload, 0);
        hipStreamWaitEvent(sB, impl.event_upload, 0);
        hipStreamWaitEvent(sC, impl.event_upload, 0);
    }

    // Claude Generated (May 2026): Stream-level timing — record start events
    if (m_record_kernel_timings) {
        hipEventRecord(impl.timing_sA_start, sA);
        hipEventRecord(impl.timing_sB_start, sB);
        hipEventRecord(impl.timing_sC_start, sC);
    }

    // --- Stream A: Charge-independent pairwise terms (dispersion, repulsion) ---
    // k_dispersion: energy mode only. Gradient mode is deferred to Phase 2
    // (computeGaussianWeightsOnGPU runs after Phase 1 on d_cn; dc6dcn is only
    // ready on stream_pairwise after the Phase-2 event_upload fence).
    if (m_dispersion_enabled && !gradient) {
        LaunchConfig cfg = getLaunchConfig(impl.disp.n, m_block_size);
        k_dispersion<<<cfg.gridSize, cfg.blockSize, 0, sA>>>(
            impl.disp.n,
            impl.disp.idx_i.ptr,  impl.disp.idx_j.ptr,
            impl.disp.C6.ptr,     impl.disp.r4r2ij.ptr,
            impl.disp.r0_sq.ptr,  impl.disp.zetac6.ptr,
            impl.disp.r_cut.ptr,
            nullptr,  // dc6dcn: energy-only, no gradient contribution
            nullptr,
            impl.coords.d_x.ptr,
            impl.coords.d_y.ptr,
            impl.coords.d_z.ptr,
            impl.d_dEdcn.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_DISP]);
    }
    if (m_repulsion_enabled) {
        // Phase 8: n > 0 guard removed — n=0 produces a 1-block no-op (graph-capture safe)
        LaunchConfig cfg = getLaunchConfig(impl.bonded_rep.n, m_block_size);
        if (impl.use_mixed_precision) {
            k_repulsion_mixed<<<cfg.gridSize, cfg.blockSize, 0, sA>>>(
                impl.bonded_rep.n,
                impl.bonded_rep.idx_i.ptr,  impl.bonded_rep.idx_j.ptr,
                impl.bonded_rep.alpha.ptr,  impl.bonded_rep.repab.ptr,
                impl.bonded_rep.r_cut.ptr,
                impl.coords.d_x.ptr,
                impl.coords.d_y.ptr,
                impl.coords.d_z.ptr,
                impl.d_grad.ptr,
                &impl.d_energies.ptr[impl.E_BREP]);
        } else {
            k_repulsion<<<cfg.gridSize, cfg.blockSize, 0, sA>>>(
                impl.bonded_rep.n,
                impl.bonded_rep.idx_i.ptr,  impl.bonded_rep.idx_j.ptr,
                impl.bonded_rep.alpha.ptr,  impl.bonded_rep.repab.ptr,
                impl.bonded_rep.r_cut.ptr,
                impl.coords.d_x.ptr,
                impl.coords.d_y.ptr,
                impl.coords.d_z.ptr,
                impl.d_grad.ptr,
                &impl.d_energies.ptr[impl.E_BREP]);
        }
    }
    if (m_repulsion_enabled) {
        // Phase 8: n > 0 guard removed — n=0 produces a 1-block no-op (graph-capture safe)
        LaunchConfig cfg = getLaunchConfig(impl.nonbonded_rep.n, m_block_size);
        if (impl.use_mixed_precision) {
            k_repulsion_mixed<<<cfg.gridSize, cfg.blockSize, 0, sA>>>(
                impl.nonbonded_rep.n,
                impl.nonbonded_rep.idx_i.ptr,  impl.nonbonded_rep.idx_j.ptr,
                impl.nonbonded_rep.alpha.ptr,  impl.nonbonded_rep.repab.ptr,
                impl.nonbonded_rep.r_cut.ptr,
                impl.coords.d_x.ptr,
                impl.coords.d_y.ptr,
                impl.coords.d_z.ptr,
                impl.d_grad.ptr,
                &impl.d_energies.ptr[impl.E_NBREP]);
        } else {
            k_repulsion<<<cfg.gridSize, cfg.blockSize, 0, sA>>>(
                impl.nonbonded_rep.n,
                impl.nonbonded_rep.idx_i.ptr,  impl.nonbonded_rep.idx_j.ptr,
                impl.nonbonded_rep.alpha.ptr,  impl.nonbonded_rep.repab.ptr,
                impl.nonbonded_rep.r_cut.ptr,
                impl.coords.d_x.ptr,
                impl.coords.d_y.ptr,
                impl.coords.d_z.ptr,
                impl.d_grad.ptr,
                &impl.d_energies.ptr[impl.E_NBREP]);
        }
    }
    // NOTE: k_coulomb is NOT launched here — it needs EEQ charges
    // Per-kernel snapshot: after stream A repulsion (Claude Generated April 2026)
    if (gradient && impl.d_grad_after_sA.ptr) {
        hipMemcpyAsync(impl.d_grad_after_sA.ptr, impl.d_grad.ptr,
                        3*N*sizeof(double), hipMemcpyDeviceToDevice, sA);
    }
    if (m_record_kernel_timings) hipEventRecord(impl.timing_sA_end, sA);
    hipEventRecord(impl.event_pairwise, sA);

    // --- Stream B: Bonded terms (bonds, angles, dihedrals, inversions, storsions, hb_alpha) ---
    // Phase 8: outer bonds.n > 0 guard removed — n=0 produces a 1-block no-op.
    // Inner hb_alpha guard kept: d_zz_hb.ptr may be null when hb_alpha.n == 0.
    if (gradient && impl.hb_alpha.n > 0) {
        impl.d_zz_hb.zero(N, sB);
    }
    {
        LaunchConfig cfg = getLaunchConfig(impl.bonds.n, m_block_size);
        k_bonds<<<cfg.gridSize, cfg.blockSize, 0, sB>>>(
            impl.bonds.n,
            impl.bonds.idx_i.ptr,     impl.bonds.idx_j.ptr,
            impl.bonds.r0.ptr,
            impl.bonds.r0_base_i.ptr, impl.bonds.r0_base_j.ptr,
            impl.bonds.cnfak_i.ptr,   impl.bonds.cnfak_j.ptr,
            impl.bonds.rabshift.ptr,
            impl.bonds.ff.ptr,
            impl.bonds.fc.ptr,
            impl.bonds.alpha.ptr,
            impl.bonds.nr_hb.ptr,
            impl.bonds.hb_cn_H.ptr,
            impl.bonds.hb_H_atom.ptr,
            impl.coords.d_x.ptr,
            impl.coords.d_y.ptr,
            impl.coords.d_z.ptr,
            impl.d_cn.ptr,
            impl.d_grad.ptr,
            impl.d_dEdcn.ptr,
            (gradient && impl.hb_alpha.n > 0) ? impl.d_zz_hb.ptr : nullptr,
            &impl.d_energies.ptr[impl.E_BOND]);
    }
    // Per-kernel snapshot: after k_bonds (Claude Generated April 2026)
    if (gradient && impl.d_grad_after_bonds.ptr) {
        hipMemcpyAsync(impl.d_grad_after_bonds.ptr, impl.d_grad.ptr,
                        3*N*sizeof(double), hipMemcpyDeviceToDevice, sB);
    }
    {
        // Phase 8: n > 0 guard removed — n=0 → 1-block no-op, graph-capture safe
        LaunchConfig cfg = getLaunchConfig(impl.angles.n, m_block_size);
        k_angles<<<cfg.gridSize, cfg.blockSize, 0, sB>>>(
            impl.angles.n,
            impl.angles.idx_i.ptr,  impl.angles.idx_j.ptr,  impl.angles.idx_k.ptr,
            impl.angles.ati.ptr,    impl.angles.atj.ptr,    impl.angles.atk.ptr,
            impl.angles.fc.ptr,     impl.angles.theta0.ptr,
            impl.coords.d_x.ptr,
            impl.coords.d_y.ptr,
            impl.coords.d_z.ptr,
            (double*)nullptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_ANGLE]);
    }
    // Per-kernel snapshot: after k_angles (Claude Generated April 2026)
    if (gradient && impl.d_grad_after_angles.ptr) {
        hipMemcpyAsync(impl.d_grad_after_angles.ptr, impl.d_grad.ptr,
                        3*N*sizeof(double), hipMemcpyDeviceToDevice, sB);
    }
    {
        // Phase 8: n > 0 guard removed — n=0 → 1-block no-op, graph-capture safe
        LaunchConfig cfg = getLaunchConfig(impl.dihedrals.n, m_block_size);
        k_dihedrals<<<cfg.gridSize, cfg.blockSize, 0, sB>>>(
            impl.dihedrals.n,
            impl.dihedrals.idx_i.ptr, impl.dihedrals.idx_j.ptr,
            impl.dihedrals.idx_k.ptr, impl.dihedrals.idx_l.ptr,
            impl.dihedrals.ati.ptr,   impl.dihedrals.atj.ptr,
            impl.dihedrals.atk.ptr,   impl.dihedrals.atl.ptr,
            impl.dihedrals.V.ptr,
            impl.dihedrals.phi0.ptr,
            impl.dihedrals.n_period.ptr,
            impl.dihedrals.is_nci.ptr,
            impl.coords.d_x.ptr,
            impl.coords.d_y.ptr,
            impl.coords.d_z.ptr,
            (double*)nullptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_DIHED]);
    }
    {
        // Phase 8: n > 0 guard removed — n=0 → 1-block no-op, graph-capture safe
        LaunchConfig cfg = getLaunchConfig(impl.inversions.n, m_block_size);
        k_inversions<<<cfg.gridSize, cfg.blockSize, 0, sB>>>(
            impl.inversions.n,
            impl.inversions.idx_i.ptr, impl.inversions.idx_j.ptr,
            impl.inversions.idx_k.ptr, impl.inversions.idx_l.ptr,
            impl.inversions.ati.ptr,   impl.inversions.atj.ptr,
            impl.inversions.atk.ptr,   impl.inversions.atl.ptr,
            impl.inversions.fc.ptr,    impl.inversions.omega0.ptr,
            impl.inversions.C0.ptr,    impl.inversions.C1.ptr,    impl.inversions.C2.ptr,
            impl.inversions.potential_type.ptr,
            impl.coords.d_x.ptr,
            impl.coords.d_y.ptr,
            impl.coords.d_z.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_INV]);
    }
    {
        // Phase 8: n > 0 guard removed — n=0 → 1-block no-op, graph-capture safe
        LaunchConfig cfg = getLaunchConfig(impl.storsions.n, m_block_size);
        k_storsions<<<cfg.gridSize, cfg.blockSize, 0, sB>>>(
            impl.storsions.n,
            impl.storsions.idx_i.ptr, impl.storsions.idx_j.ptr,
            impl.storsions.idx_k.ptr, impl.storsions.idx_l.ptr,
            impl.storsions.erefhalf.ptr,
            impl.coords.d_x.ptr,
            impl.coords.d_y.ptr,
            impl.coords.d_z.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_STORS]);
    }
    if (gradient && impl.hb_alpha.n > 0 && impl.d_zz_hb.ptr) {
        constexpr double hb_kn = 27.5;
        LaunchConfig cfg = getLaunchConfig(impl.hb_alpha.n, m_block_size);
        k_hb_alpha_chainrule<<<cfg.gridSize, cfg.blockSize, 0, sB>>>(
            impl.hb_alpha.n,
            impl.hb_alpha.idx_H.ptr,
            impl.hb_alpha.idx_B.ptr,
            impl.hb_alpha.rcov_sum.ptr,
            impl.coords.d_x.ptr,
            impl.coords.d_y.ptr,
            impl.coords.d_z.ptr,
            impl.d_zz_hb.ptr,
            impl.d_grad.ptr,
            hb_kn);
    }
    // Per-kernel snapshot: after all stream B kernels (Claude Generated April 2026)
    if (gradient && impl.d_grad_after_sB.ptr) {
        hipMemcpyAsync(impl.d_grad_after_sB.ptr, impl.d_grad.ptr,
                        3*N*sizeof(double), hipMemcpyDeviceToDevice, sB);
    }
    if (m_record_kernel_timings) hipEventRecord(impl.timing_sB_end, sB);
    hipEventRecord(impl.event_bonded, sB);

    // --- Stream C: 3-body terms (batm, atm, xbonds, hbonds) ---
    // All charge-independent: batm uses topology charges (static),
    // hbonds uses q_H/q_A/q_B baked into SoA at construction time.
    {
        // Phase 8: n > 0 guard removed — n=0 → 1-block no-op, graph-capture safe
        LaunchConfig cfg = getLaunchConfig(impl.batm.n, m_block_size);
        if (impl.use_mixed_precision) {
            k_batm_mixed<<<cfg.gridSize, cfg.blockSize, 0, sC>>>(
                impl.batm.n,
                impl.batm.idx_i.ptr, impl.batm.idx_j.ptr, impl.batm.idx_k.ptr,
                impl.batm.zb3atm_i.ptr, impl.batm.zb3atm_j.ptr, impl.batm.zb3atm_k.ptr,
                impl.coords.d_x.ptr,
                impl.coords.d_y.ptr,
                impl.coords.d_z.ptr,
                impl.d_topo_charges.ptr,
                impl.d_grad.ptr,
                &impl.d_energies.ptr[impl.E_BATM]);
        } else {
            k_batm<<<cfg.gridSize, cfg.blockSize, 0, sC>>>(
                impl.batm.n,
                impl.batm.idx_i.ptr, impl.batm.idx_j.ptr, impl.batm.idx_k.ptr,
                impl.batm.zb3atm_i.ptr, impl.batm.zb3atm_j.ptr, impl.batm.zb3atm_k.ptr,
                impl.coords.d_x.ptr,
                impl.coords.d_y.ptr,
                impl.coords.d_z.ptr,
                impl.d_topo_charges.ptr,
                impl.d_grad.ptr,
                &impl.d_energies.ptr[impl.E_BATM]);
        }
    }
    {
        // Phase 8: n > 0 guard removed — n=0 → 1-block no-op, graph-capture safe
        LaunchConfig cfg = getLaunchConfig(impl.atm.n, m_block_size);
        k_atm<<<cfg.gridSize, cfg.blockSize, 0, sC>>>(
            impl.atm.n,
            impl.atm.idx_i.ptr, impl.atm.idx_j.ptr, impl.atm.idx_k.ptr,
            impl.atm.ati.ptr,   impl.atm.atj.ptr,   impl.atm.atk.ptr,
            impl.atm.C6_ij.ptr, impl.atm.C6_ik.ptr, impl.atm.C6_jk.ptr,
            impl.atm.s9.ptr, impl.atm.a1.ptr, impl.atm.a2.ptr, impl.atm.alp.ptr,
            impl.atm.triple_scale.ptr,
            impl.coords.d_x.ptr,
            impl.coords.d_y.ptr,
            impl.coords.d_z.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_ATM]);
    }
    {
        // Phase 8: n > 0 guard removed — n=0 → 1-block no-op, graph-capture safe
        LaunchConfig cfg = getLaunchConfig(impl.xbonds.n, m_block_size);
        if (impl.use_mixed_precision) {
            k_xbonds_mixed<<<cfg.gridSize, cfg.blockSize, 0, sC>>>(
                impl.xbonds.n,
                impl.xbonds.idx_i.ptr, impl.xbonds.idx_j.ptr, impl.xbonds.idx_k.ptr,
                impl.xbonds.elem_A.ptr, impl.xbonds.elem_B.ptr,
                impl.xbonds.q_X.ptr, impl.xbonds.q_B.ptr, impl.xbonds.acidity_X.ptr,
                impl.xbonds.r_cut.ptr,
                impl.coords.d_x.ptr,
                impl.coords.d_y.ptr,
                impl.coords.d_z.ptr,
                impl.d_grad.ptr,
                &impl.d_energies.ptr[impl.E_XBOND]);
        } else {
            k_xbonds<<<cfg.gridSize, cfg.blockSize, 0, sC>>>(
                impl.xbonds.n,
                impl.xbonds.idx_i.ptr, impl.xbonds.idx_j.ptr, impl.xbonds.idx_k.ptr,
                impl.xbonds.elem_A.ptr, impl.xbonds.elem_B.ptr,
                impl.xbonds.q_X.ptr, impl.xbonds.q_B.ptr, impl.xbonds.acidity_X.ptr,
                impl.xbonds.r_cut.ptr,
                impl.coords.d_x.ptr,
                impl.coords.d_y.ptr,
                impl.coords.d_z.ptr,
                impl.d_grad.ptr,
                &impl.d_energies.ptr[impl.E_XBOND]);
        }
    }
    if (m_hbond_enabled) {
        // Phase 8: hbonds.n > 0 guard removed — n=0 → 1-block no-op, graph-capture safe
        LaunchConfig cfg = getLaunchConfig(impl.hbonds.n, m_block_size);
        k_hbonds<<<cfg.gridSize, cfg.blockSize, 0, sC>>>(
            impl.hbonds.n,
            impl.hbonds.idx_i.ptr, impl.hbonds.idx_j.ptr, impl.hbonds.idx_k.ptr,
            impl.hbonds.elem_A.ptr, impl.hbonds.elem_B.ptr,
            impl.hbonds.q_H.ptr, impl.hbonds.q_A.ptr, impl.hbonds.q_B.ptr,
            impl.hbonds.basicity_A.ptr, impl.hbonds.basicity_B.ptr,
            impl.hbonds.acidity_A.ptr,  impl.hbonds.acidity_B.ptr,
            impl.hbonds.r_cut.ptr,
            impl.hbonds.case_type.ptr,
            impl.hbonds.nb_B_offset.ptr, impl.hbonds.nb_B_count.ptr,
            impl.hbonds.nb_B_flat.ptr,
            impl.hbonds.acceptor_parent.ptr,
            impl.hbonds.nb_C_offset.ptr, impl.hbonds.nb_C_count.ptr,
            impl.hbonds.nb_C_flat.ptr,
            impl.hbonds.repz_B.ptr,
            impl.coords.d_x.ptr,
            impl.coords.d_y.ptr,
            impl.coords.d_z.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_HBOND]);
    }
    // Per-kernel snapshot: after all stream C kernels (Claude Generated April 2026)
    if (gradient && impl.d_grad_after_sC.ptr) {
        hipMemcpyAsync(impl.d_grad_after_sC.ptr, impl.d_grad.ptr,
                        3*N*sizeof(double), hipMemcpyDeviceToDevice, sC);
    }
    if (m_record_kernel_timings) hipEventRecord(impl.timing_sC_end, sC);
    hipEventRecord(impl.event_threebody, sC);

    // === CAPTURE END: finalise graph and execute immediately ===
    if (capturing) {
        // Join all forked streams back to main before EndCapture.
        // These are GPU-side edges only; hipGraphLaunch() below is still async from CPU.
        hipStreamWaitEvent(stream, impl.event_pairwise,  0);
        hipStreamWaitEvent(stream, impl.event_bonded,    0);
        hipStreamWaitEvent(stream, impl.event_threebody, 0);

        hipError_t end_err = hipStreamEndCapture(stream, &impl.m_graph_phase1);
        if (end_err == hipSuccess && impl.m_graph_phase1) {
            hipError_t inst_err = hipGraphInstantiate(
                &impl.m_graph_exec_phase1, impl.m_graph_phase1, NULL, NULL, 0);
            if (inst_err == hipSuccess) {
                impl.m_graph_phase1_valid = true;
                // hipStreamBeginCapture suppressed actual kernel execution; replay now.
                hipGraphLaunch(impl.m_graph_exec_phase1, stream);
                return;
            }
        }
        // Capture or instantiation failed — fall back silently (normal execution already done)
        if (impl.m_graph_phase1) { hipGraphDestroy(impl.m_graph_phase1); impl.m_graph_phase1 = nullptr; }
    }
    // Returns immediately — charge-independent kernels run asynchronously on GPU
}

// Claude Generated (June 2026): OPT-IN device-resident HB charge refresh.
// IMPORTANT — default is OFF. The GFN-FF reference (Fortran) and the curcuma CPU path
// both evaluate the H-bond term with the charges FROZEN at topology build (CPU uses the
// stored hb.q_H/q_A/q_B, forcefieldthread.cpp:2525+), NOT live per-step EEQ charges. The
// GPU's frozen HB charges therefore MATCH the reference; refreshing them to live charges
// makes MD diverge from CPU/Fortran (verified: acetic-dimer 1000 fs heat-exchange CPU
// 0.00709227 / GPU-frozen 0.00707467 / GPU-resident 0.00658143). This path is kept only
// as an opt-in experiment — enable with CURCUMA_GFNFF_GPU_RESIDENT_HBQ=1.
// Called by the orchestrator AFTER launchChargeDependentAndFinish() returns (impl.stream
// idle, d_charges final) → no stream hazard with the captured graphs.
void FFWorkspaceHip::refreshHBChargesFromDevice()
{
    static const bool resident = (std::getenv("CURCUMA_GFNFF_GPU_RESIDENT_HBQ") != nullptr);
    if (!resident || !m_hbond_enabled)
        return;
    auto& impl = *m_impl;
    if (impl.hbonds.n <= 0)
        return;
    LaunchConfig cfg = getLaunchConfig(impl.hbonds.n, m_block_size);
    k_gather_hb_charges<<<cfg.gridSize, cfg.blockSize, 0, impl.stream>>>(
        impl.hbonds.n,
        impl.hbonds.idx_i.ptr, impl.hbonds.idx_j.ptr, impl.hbonds.idx_k.ptr,
        impl.d_charges.ptr,
        impl.hbonds.q_A.ptr, impl.hbonds.q_H.ptr, impl.hbonds.q_B.ptr);
    checkCuda(hipGetLastError(), "k_gather_hb_charges launch");
    // Ensure the q arrays are fully written before the next step's HB (Phase 1) reads them.
    hipStreamSynchronize(impl.stream);
}

double FFWorkspaceHip::launchChargeDependentAndFinish(bool gradient)
{
    auto& impl = *m_impl;
    const int N = m_natoms;
    hipStream_t stream = impl.stream;

    // Claude Generated (May 2026): Phase 2 timing start
    if (m_record_kernel_timings) hipEventRecord(impl.timing_p2_start, stream);

    // =========================================================================
    // WP4 (May 2026): Phase-2 CUDA-Graph capture / replay
    //   Two graph slots — one for gradient=true, one for energy-only — kept in step
    //   with the topology cache (invalidateGraph() destroys both alongside Phase 1).
    //   Capture guard matches Phase 1: stream_pairwise must be idle (no pending Phase 1
    //   non-captured work) before it can be joined to a capture via hipStreamWaitEvent.
    //   d_grad_snapshot is deliberately kept allocated (see Phase 1 comment at
    //   prepareAndLaunchChargeIndependent()) so this guard always disables capture
    //   until stream isolation is solved.
    // Phase-1 graph flag: if Phase-1 ran as a graph, all its sub-streams (sA, sB, sC)
    // are fully synchronized to the main stream before Phase-2 entries execute.
    // Phase-1-internal events (event_bonded, event_threebody) are NOT externally signaled
    // after graph replay → Phase-2 must NOT wait on them in that case.
    const bool p1_was_graph = impl.m_graph_phase1_valid;

    // WP5-A fix: d_grad_snapshot is only allocated at verbosity >= 3 now.
    // p2_snapshot_guard is false in normal mode → capture guard is pure verbosity check.
    const bool p2_snapshot_guard      = impl.d_dEdcn_snapshot.ptr && impl.d_grad_snapshot.ptr;
    const bool p2_need_snapshots      = (m_verbosity >= 3) || p2_snapshot_guard;
    const bool p2_force_single_stream = (gradient && impl.d_grad_after_bonds.ptr);
    // Phase-1 graph fix: when Phase-1 ran as a graph, stream_pairwise may still have
    // pending graph work when Phase-2 tries to begin capture. CUDA rejects the
    // hipStreamWaitEvent join with hipErrorStreamCaptureIsolation, invalidating
    // the Phase-2 capture and causing all charge-dependent kernels to be silently
    // skipped → energy = 0. Disable Phase-2 capture whenever Phase-1 is graph-active.
    // Claude Generated (May 2026): m_record_kernel_timings disables graph for accurate per-step events.
    const bool p2_capture_disabled    = p2_need_snapshots || p2_force_single_stream || p1_was_graph
                                          || m_record_kernel_timings;
    bool&            p2_valid = gradient ? impl.m_graph_phase2_grad_valid
                                         : impl.m_graph_phase2_energy_valid;
    hipGraph_t&     p2_graph = gradient ? impl.m_graph_phase2_grad
                                         : impl.m_graph_phase2_energy;
    hipGraphExec_t& p2_exec  = gradient ? impl.m_graph_exec_phase2_grad
                                         : impl.m_graph_exec_phase2_energy;
    bool p2_capturing = false;

    if (p2_valid && !p2_capture_disabled) {
        // Replay: hipGraphLaunch is async; sync + extraction happen below.
        hipGraphLaunch(p2_exec, stream);
    } else {
        if (!p2_valid && !p2_capture_disabled) {
            hipError_t cap_err = hipStreamBeginCapture(stream, hipStreamCaptureModeThreadLocal);
            p2_capturing = (cap_err == hipSuccess);
        }

    // =========================================================================
    // 1. Upload EEQ charges
    //    WP5-A: if setEEQDeviceCharges() was called, charges are already in
    //    impl.d_charges via D2D — skip H2D upload and reset the flag.
    // =========================================================================
    if (m_device_charges_ready) {
        m_device_charges_ready = false;  // consumed; next step will upload normally
    } else if (m_eeq_charges.size() == N) {
        impl.d_charges.upload(m_eeq_charges.data(), N, stream);
    }

    // =========================================================================
    // 2. Launch deferred kernels on stream A
    //    Wait for charge-independent pairwise kernels to finish first,
    //    since these kernels write to the same d_grad/d_energies buffers.
    // =========================================================================
    // Record event on main stream (after charge upload) so pairwise stream can wait.
    // WP4: dedicated phase-2 event (event_p2_upload) to keep capture isolated from
    // phase-1's event_upload, which is recorded inside the phase-1 graph.
    hipEventRecord(impl.event_p2_upload, stream);
    hipStreamWaitEvent(impl.stream_pairwise, impl.event_p2_upload, 0);

    // k_dispersion (gradient mode): deferred from Phase 1.
    // computeGaussianWeightsOnGPU() was called after Phase 1 (using d_cn, current step)
    // and runs on the main stream. The event_upload fence above ensures stream_pairwise
    // sees the completed dc6dcn values before k_dispersion launches.
    if (m_dispersion_enabled && gradient) {
        LaunchConfig cfg = getLaunchConfig(impl.disp.n, m_block_size);
        k_dispersion<<<cfg.gridSize, cfg.blockSize, 0, impl.stream_pairwise>>>(
            impl.disp.n,
            impl.disp.idx_i.ptr,  impl.disp.idx_j.ptr,
            impl.disp.C6.ptr,     impl.disp.r4r2ij.ptr,
            impl.disp.r0_sq.ptr,  impl.disp.zetac6.ptr,
            impl.disp.r_cut.ptr,
            impl.disp.dc6dcn_ij.ptr,
            impl.disp.dc6dcn_ji.ptr,
            impl.coords.d_x.ptr,
            impl.coords.d_y.ptr,
            impl.coords.d_z.ptr,
            impl.d_dEdcn.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_DISP]);
        checkCuda(hipGetLastError(), "k_dispersion launch");
    }

    // k_coulomb (needs EEQ charges)
    if (m_coulomb_enabled && impl.coulomb.n > 0) {
        // ROCm: the pair k_coulomb's ~6M FP64 atomicAdd dominated the MD step. Preferred path is
        // the dense, coalesced, shared-mem-tiled N-body kernel (k_coulomb_dense) which also fixes
        // the gather's random-memory bottleneck; it needs the per-atom charge-corrected EEQ alpha
        // (== the source of the precomputed gamma_ij) to recompute gamma on the fly. If the EEQ
        // topology alpha is not resident, fall back to the CSR gather (still atomic-light).
        const int coul_block = 128;
        if (impl.eeq_topo.valid && impl.eeq_topo.d_alpha_corrected.ptr) {
            int grid = (impl.coulomb.csr_N + coul_block - 1) / coul_block;
            k_coulomb_dense<<<grid, coul_block, 0, impl.stream_pairwise>>>(
                impl.coulomb.csr_N,
                impl.eeq_topo.d_alpha_corrected.ptr,
                impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
                impl.d_charges.ptr, impl.coulomb.r_cut_scalar,
                impl.d_grad.ptr,
                &impl.d_energies.ptr[impl.E_COUL]);
            checkCuda(hipGetLastError(), "k_coulomb_dense launch");
        } else {
            k_coulomb_gather<<<impl.coulomb.csr_N, coul_block, 0, impl.stream_pairwise>>>(
                impl.coulomb.csr_N,
                impl.coulomb.csr_off.ptr, impl.coulomb.csr_partner.ptr,
                impl.coulomb.csr_gamma.ptr, impl.coulomb.csr_rcut.ptr,
                impl.coords.d_x.ptr,
                impl.coords.d_y.ptr,
                impl.coords.d_z.ptr,
                impl.d_charges.ptr,
                impl.d_grad.ptr,
                &impl.d_energies.ptr[impl.E_COUL]);
            checkCuda(hipGetLastError(), "k_coulomb_gather launch");
        }
    }
    // Re-record pairwise event to include k_coulomb completion.
    // WP4: event_p2_pairwise — phase-2-dedicated, never aliased with phase-1 graph.
    hipEventRecord(impl.event_p2_pairwise, impl.stream_pairwise);

    // =========================================================================
    // 3. Postprocess on main stream
    // Phase 4b+Fusion: Fused Coulomb postprocess (Claude Generated March 2026)
    //   k_coulomb_postprocess — fused self-energy + qtmp subtraction in single O(N) pass
    //     Needs d_dEdcn (from pairwise + bonded) for qtmp, NOT threebody
    //   k_cn_chainrule   — needs d_dEdcn complete from ALL streams → wait for all 3
    // =========================================================================

    const bool need_snapshots = (m_verbosity >= 3);
    if (m_coulomb_enabled && impl.coul_self_on_gpu) {
        // Wait for pairwise+bonded streams (dEdcn dependency for qtmp subtraction).
        // event_p2_pairwise: Phase-2-dedicated — always wait (records k_coulomb completion).
        // event_bonded: Phase-1-internal — skip if Phase-1 ran as graph (already sequenced).
        hipStreamWaitEvent(stream, impl.event_p2_pairwise, 0);
        if (!p1_was_graph)
            hipStreamWaitEvent(stream, impl.event_bonded, 0);

        // d_dEdcn_snapshot must always be populated — dEdcnTotal() is public API
        // used by external code (tests, chain-rule validation) regardless of verbosity.
        if (m_coulomb_enabled && impl.coul_self_on_gpu && impl.d_dEdcn_snapshot.ptr) {
            hipMemcpyAsync(impl.d_dEdcn_snapshot.ptr, impl.d_dEdcn.ptr,
                            N * sizeof(double), hipMemcpyDeviceToDevice, stream);
        }

        if (need_snapshots) {
            // Grad snapshot is only for diagnostics at verbosity >= 3
            hipMemcpyAsync(impl.d_grad_snapshot.ptr, impl.d_grad.ptr,
                            3*N*sizeof(double), hipMemcpyDeviceToDevice, stream);
        }

        // WP5-B: do_subtract only needs the GPU CNF buffer + EEQ charges — no CPU vector.
        bool do_subtract = gradient && impl.coul_self_on_gpu && m_eeq_charges.size() == N;
        LaunchConfig cfg = getLaunchConfig(N, m_block_size);
        k_coulomb_postprocess<<<cfg.gridSize, cfg.blockSize, 0, stream>>>(
            N,
            impl.d_charges.ptr,
            impl.d_coul_chi_base.ptr,
            impl.eeq_topo.d_cnf.ptr,
            impl.d_cn.ptr,
            impl.d_coul_gam.ptr,
            impl.d_coul_alp.ptr,
            impl.d_dEdcn.ptr,
            &impl.d_energies.ptr[impl.E_COUL_SELF],
            do_subtract);
        checkCuda(hipGetLastError(), "k_coulomb_postprocess launch");
    }

    // k_cn_chainrule: reads d_dEdcn (must be final) and d_grad (accumulated by ALL
    // kernels), so wait for all 3 stream events.
    // Static-Mode (WP-S1): skip when frozen_cn — Term 1b is zero by definition (dcn=0).
    if (gradient && m_cn_pairs_on_gpu && impl.n_cn_pairs > 0 && impl.d_dlogdcn.ptr
        && !m_frozen_cn) {
        // event_p2_pairwise: Phase-2-dedicated — always wait.
        // event_bonded/threebody: Phase-1-internal — skip if Phase-1 ran as graph.
        hipStreamWaitEvent(stream, impl.event_p2_pairwise, 0);
        if (!p1_was_graph) {
            hipStreamWaitEvent(stream, impl.event_bonded, 0);
            hipStreamWaitEvent(stream, impl.event_threebody, 0);
        }

        if (need_snapshots) {
            hipMemcpyAsync(impl.d_grad_snapshot.ptr, impl.d_grad.ptr,
                            3*N*sizeof(double), hipMemcpyDeviceToDevice, stream);
        }

        LaunchConfig cfg_cn = getLaunchConfig(impl.n_cn_pairs, m_block_size);
        k_cn_chainrule<<<cfg_cn.gridSize, cfg_cn.blockSize, 0, stream>>>(
            impl.n_cn_pairs,
            impl.d_cn_idx_i.ptr,
            impl.d_cn_idx_j.ptr,
            impl.d_cn_rcov_sum.ptr,
            impl.coords.d_x.ptr,
            impl.coords.d_y.ptr,
            impl.coords.d_z.ptr,
            impl.d_dEdcn.ptr,
            impl.d_dlogdcn.ptr,
            impl.d_grad.ptr,
            m_kn);
        checkCuda(hipGetLastError(), "k_cn_chainrule launch");
    }

    // =========================================================================
    // 4. Async download of final results (Phase 6: Async DMA)
    // Claude Generated (March 2026): All D2H transfers enqueued on main stream,
    // single hipStreamSynchronize at end. Pinned buffers enable true async DMA.
    // =========================================================================
    // Final fence before D2H: wait for all Phase-2 stream work.
    // event_p2_pairwise: Phase-2-dedicated (dispersion + Coulomb on stream_pairwise).
    // event_bonded/threebody: Phase-1-internal — skip if Phase-1 ran as graph.
    hipStreamWaitEvent(stream, impl.event_p2_pairwise, 0);
    if (!p1_was_graph) {
        hipStreamWaitEvent(stream, impl.event_bonded, 0);
        hipStreamWaitEvent(stream, impl.event_threebody, 0);
    }

    checkCuda(hipGetLastError(), "postprocess final check");

    // Async energy download (pinned m_h_energies buffer)
    checkCuda(hipMemcpyAsync(m_h_energies, impl.d_energies.ptr,
                              FFWorkspaceHipImpl::N_ENERGY_SLOTS * sizeof(double),
                              hipMemcpyDeviceToHost, stream),
              "async energy download");

    // Async gradient + diagnostic snapshot downloads
    if (gradient) {
        checkCuda(hipMemcpyAsync(m_h_grad, impl.d_grad.ptr,
                                  3 * N * sizeof(double),
                                  hipMemcpyDeviceToHost, stream),
                  "async gradient download");

        // Always download dEdcn snapshot — dEdcnTotal() is public API
        checkCuda(hipMemcpyAsync(m_h_dEdcn_snap, impl.d_dEdcn_snapshot.ptr,
                                  N * sizeof(double),
                                  hipMemcpyDeviceToHost, stream),
                  "async dEdcn snapshot download");

        if (need_snapshots) {
            checkCuda(hipMemcpyAsync(m_h_grad_snap, impl.d_grad_snapshot.ptr,
                                      3 * N * sizeof(double),
                                      hipMemcpyDeviceToHost, stream),
                      "async grad snapshot download");
        }
        // Per-kernel snapshot downloads (Claude Generated April 2026)
        if (impl.d_grad_after_sA.ptr) {
            checkCuda(hipMemcpyAsync(m_h_grad_after_sA, impl.d_grad_after_sA.ptr,
                                      3 * N * sizeof(double),
                                      hipMemcpyDeviceToHost, stream),
                      "async grad after_sA download");
        }
        if (impl.d_grad_after_bonds.ptr) {
            checkCuda(hipMemcpyAsync(m_h_grad_after_bonds, impl.d_grad_after_bonds.ptr,
                                      3 * N * sizeof(double),
                                      hipMemcpyDeviceToHost, stream),
                      "async grad after_bonds download");
        }
        if (impl.d_grad_after_angles.ptr) {
            checkCuda(hipMemcpyAsync(m_h_grad_after_angles, impl.d_grad_after_angles.ptr,
                                      3 * N * sizeof(double),
                                      hipMemcpyDeviceToHost, stream),
                      "async grad after_angles download");
        }
        if (impl.d_grad_after_sB.ptr) {
            checkCuda(hipMemcpyAsync(m_h_grad_after_sB, impl.d_grad_after_sB.ptr,
                                      3 * N * sizeof(double),
                                      hipMemcpyDeviceToHost, stream),
                      "async grad after_sB download");
        }
        if (impl.d_grad_after_sC.ptr) {
            checkCuda(hipMemcpyAsync(m_h_grad_after_sC, impl.d_grad_after_sC.ptr,
                                      3 * N * sizeof(double),
                                      hipMemcpyDeviceToHost, stream),
                      "async grad after_sC download");
        }
    }

        // === WP4: Phase-2 capture end + immediate replay ===
        // BeginCapture suppressed real kernel execution; the captured graph must be
        // launched once now to actually run this step's work (matches Phase 1 pattern).
        if (p2_capturing) {
            hipError_t end_err = hipStreamEndCapture(stream, &p2_graph);
            if (end_err == hipSuccess && p2_graph) {
                hipError_t inst_err = hipGraphInstantiate(
                    &p2_exec, p2_graph, NULL, NULL, 0);
                if (inst_err == hipSuccess) {
                    p2_valid = true;
                    hipGraphLaunch(p2_exec, stream);
                } else {
                    hipGraphDestroy(p2_graph);
                    p2_graph = nullptr;
                }
            } else if (p2_graph) {
                hipGraphDestroy(p2_graph);
                p2_graph = nullptr;
            }
        }
    } // end else (WP4 replay-or-capture)

    // Claude Generated (May 2026): Phase 2 timing end (before sync, so it captures kernel time only)
    if (m_record_kernel_timings) hipEventRecord(impl.timing_p2_end, stream);

    // Single synchronization point — all async transfers complete here
    checkCuda(hipStreamSynchronize(stream), "final stream sync");

    // Claude Generated (May 2026): Compute per-stream elapsed times after sync.
    // Stream-level totals are mapped to per-component timings (all kernels in a
    // stream share its total). True per-kernel timing would require per-kernel
    // event pairs; stream-level grouping is the pragmatic compromise.
    if (m_record_kernel_timings) {
        m_kernel_timings.reset();
        float ms_sA = 0, ms_sB = 0, ms_sC = 0, ms_p2 = 0;
        hipError_t e1 = hipEventElapsedTime(&ms_sA, impl.timing_sA_start, impl.timing_sA_end);
        hipError_t e2 = hipEventElapsedTime(&ms_sB, impl.timing_sB_start, impl.timing_sB_end);
        hipError_t e3 = hipEventElapsedTime(&ms_sC, impl.timing_sC_start, impl.timing_sC_end);
        hipError_t e4 = hipEventElapsedTime(&ms_p2, impl.timing_p2_start, impl.timing_p2_end);
        // Map stream times to per-component slots so the unified report's GPU column shows non-empty values.
        // Components within the same stream share the stream's wall-clock — this is intentional
        // (per-kernel timing not available without per-kernel event pairs).
        if (e1 == hipSuccess) {
            m_kernel_timings.dispersion    = ms_sA;
            m_kernel_timings.bonded_rep    = ms_sA;
            m_kernel_timings.nonbonded_rep = ms_sA;
        }
        if (e2 == hipSuccess) {
            m_kernel_timings.bonds      = ms_sB;
            m_kernel_timings.angles     = ms_sB;
            m_kernel_timings.dihedrals  = ms_sB;
            m_kernel_timings.inversions = ms_sB;
            m_kernel_timings.stors      = ms_sB;
        }
        if (e3 == hipSuccess) {
            m_kernel_timings.atm   = ms_sC;
            m_kernel_timings.batm  = ms_sC;
            m_kernel_timings.hbond = ms_sC;
            m_kernel_timings.xbond = ms_sC;
        }
        if (e4 == hipSuccess) {
            m_kernel_timings.coulomb = ms_p2;
        }
    }

    // WP5-C (May 2026): Update skip-pending flag from async D2H copy.
    // h_dc6dcn_skip was written by k_check_dc6dcn_skip_final → hipMemcpyAsync on the
    // same stream; now that the stream is synchronized, the host buffer is valid.
    m_dc6dcn_skip_pending = (impl.h_dc6dcn_skip != nullptr && *impl.h_dc6dcn_skip != 0);

    // Extract per-term energies from pinned buffer
    const double e_disp      = m_h_energies[impl.E_DISP];
    const double e_brep      = m_h_energies[impl.E_BREP];
    const double e_nbrep     = m_h_energies[impl.E_NBREP];
    const double e_coul1     = m_h_energies[impl.E_COUL];
    const double e_bond      = m_h_energies[impl.E_BOND];
    const double e_angle     = m_h_energies[impl.E_ANGLE];
    const double e_dihedral  = m_h_energies[impl.E_DIHED];
    const double e_inversion = m_h_energies[impl.E_INV];
    const double e_stors     = m_h_energies[impl.E_STORS];
    const double e_batm      = m_h_energies[impl.E_BATM];
    const double e_atm       = m_h_energies[impl.E_ATM];
    const double e_xbond     = m_h_energies[impl.E_XBOND];
    const double e_hbond     = m_h_energies[impl.E_HBOND];
    const double e_coul_self = m_h_energies[impl.E_COUL_SELF];

    // Copy gradient from pinned buffer to Eigen matrix
    if (gradient) {
        for (int i = 0; i < N; ++i) {
            m_result_gradient(i, 0) = m_h_grad[3*i+0];
            m_result_gradient(i, 1) = m_h_grad[3*i+1];
            m_result_gradient(i, 2) = m_h_grad[3*i+2];
        }

        // GPU POST-CN GRADIENT (verbosity >= 3)
        if (m_verbosity >= 3) {
            CurcumaLogger::info(fmt::format("=== GPU POST-CN GRADIENT (natoms={}) ===", N));
            for (int i = 0; i < N; ++i) {
                fmt::print("  atom {:3d}: {:+.15e} {:+.15e} {:+.15e}\n",
                           i, m_result_gradient(i,0), m_result_gradient(i,1), m_result_gradient(i,2));
            }
        }

        // Per-kernel gradient decomposition (Claude Generated April 2026)
        // With single-stream mode, snapshots are clean: sA → bonds → angles → rest_sB → sC
        if (m_h_grad_after_sA && m_h_grad_after_bonds && m_h_grad_after_angles
            && m_h_grad_after_sB && m_h_grad_after_sC) {
            if (m_verbosity >= 3) {
                CurcumaLogger::info("=== GPU PER-KERNEL GRADIENT DECOMPOSITION (single-stream) ===");
                fmt::print("  Atom   |rep(sA)|    |bonds|       |angles|      |sB_rest|    |3body(sC)|  |coul_ph2|   |CN_chain|\n");
            }
            for (int i = 0; i < N; ++i) {
                auto norm3 = [](const double* buf, int idx) {
                    double x = buf[3*idx+0], y = buf[3*idx+1], z = buf[3*idx+2];
                    return std::sqrt(x*x + y*y + z*z);
                };
                auto delta_norm = [](const double* a, const double* b, int idx) {
                    double dx = a[3*idx+0]-b[3*idx+0], dy = a[3*idx+1]-b[3*idx+1], dz = a[3*idx+2]-b[3*idx+2];
                    return std::sqrt(dx*dx + dy*dy + dz*dz);
                };

                double rep     = norm3(m_h_grad_after_sA, i);
                double bonds   = delta_norm(m_h_grad_after_bonds, m_h_grad_after_sA, i);
                double angles  = delta_norm(m_h_grad_after_angles, m_h_grad_after_bonds, i);
                double rest_sB = delta_norm(m_h_grad_after_sB, m_h_grad_after_angles, i);
                double threebody = delta_norm(m_h_grad_after_sC, m_h_grad_after_sB, i);
                // Coulomb phase 2 = pre_CN - after_sC (disp_grad + coulomb + coulomb_postprocess)
                double coul_ph2 = delta_norm(m_h_grad_snap, m_h_grad_after_sC, i);
                // CN chain-rule = post_CN - pre_CN
                double cn_chain = delta_norm(m_h_grad, m_h_grad_snap, i);

                if (m_verbosity >= 3) {
                    fmt::print("  {:2d}  {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e}\n",
                               i, rep, bonds, angles, rest_sB, threebody, coul_ph2, cn_chain);
                }
            }
            // Print full gradient vectors for each decomposed component (verbosity >= 3)
            if (m_verbosity >= 3) {
                const char* labels[] = {"REPULSION(sA)", "BONDS", "ANGLES", "sB_REST(dihed+inv+stors+hb_alpha)", "3BODY(sC:batm+atm+xb+hb)", "COULOMB_PH2(disp_gr+coul+postproc)", "CN_CHAINRULE"};
                const double* bases[] = {nullptr, m_h_grad_after_sA, m_h_grad_after_bonds, m_h_grad_after_angles, m_h_grad_after_sB, m_h_grad_after_sC, m_h_grad_snap};
                const double* snaps[] = {m_h_grad_after_sA, m_h_grad_after_bonds, m_h_grad_after_angles, m_h_grad_after_sB, m_h_grad_after_sC, m_h_grad_snap, m_h_grad};
                for (int c = 0; c < 7; ++c) {
                    CurcumaLogger::info(fmt::format("=== GPU {} GRADIENT ===", labels[c]));
                    for (int i = 0; i < N; ++i) {
                        double gx = snaps[c][3*i+0] - (bases[c] ? bases[c][3*i+0] : 0.0);
                        double gy = snaps[c][3*i+1] - (bases[c] ? bases[c][3*i+1] : 0.0);
                        double gz = snaps[c][3*i+2] - (bases[c] ? bases[c][3*i+2] : 0.0);
                        fmt::print("  atom {:3d}: {:+.15e} {:+.15e} {:+.15e}\n", i, gx, gy, gz);
                    }
                }
                CurcumaLogger::info("=== GPU PER-KERNEL END ===");
            }
        }

        // Always populate m_dEdcn_total — dEdcnTotal() is public API
        m_dEdcn_total.resize(N);
        for (int i = 0; i < N; ++i)
            m_dEdcn_total[i] = m_h_dEdcn_snap[i];

        // GPU dEdcn BEFORE qtmp (verbosity >= 3)
        if (m_verbosity >= 3) {
            CurcumaLogger::info(fmt::format("=== GPU dEdcn BEFORE qtmp (natoms={}) ===", N));
            for (int i = 0; i < N; ++i) {
                fmt::print("  atom {:3d}: dEdcn={:+.15e}\n", i, m_dEdcn_total[i]);
            }
            CurcumaLogger::info("=== GPU dEdcn END ===");
        }

        const bool need_snapshots = (m_verbosity >= 3);
        if (need_snapshots) {
            m_grad_before_cn.resize(N, 3);
            for (int i = 0; i < N; ++i) {
                m_grad_before_cn(i, 0) = m_h_grad_snap[3*i+0];
                m_grad_before_cn(i, 1) = m_h_grad_snap[3*i+1];
                m_grad_before_cn(i, 2) = m_h_grad_snap[3*i+2];
            }
            // GPU PRE-CN GRADIENT (verbosity >= 3)
            if (m_verbosity >= 3) {
                CurcumaLogger::info(fmt::format("=== GPU PRE-CN GRADIENT (natoms={}) ===", N));
                for (int i = 0; i < N; ++i) {
                    fmt::print("  atom {:3d}: {:+.15e} {:+.15e} {:+.15e}\n",
                               i, m_grad_before_cn(i,0), m_grad_before_cn(i,1), m_grad_before_cn(i,2));
                }
                CurcumaLogger::info("=== GPU PRE-CN GRADIENT END ===");
            }
        }
    }

    // =========================================================================
    // 5. Populate energy components
    // =========================================================================
    m_result_energy.reset();
    m_result_energy.dispersion    = e_disp;
    m_result_energy.bonded_rep    = e_brep;
    m_result_energy.nonbonded_rep = e_nbrep;
    m_result_energy.coulomb       = e_coul1 + e_coul_self;
    m_result_energy.bond          = e_bond;
    m_result_energy.angle         = e_angle;
    m_result_energy.dihedral      = e_dihedral;
    m_result_energy.inversion     = e_inversion;
    m_result_energy.stors         = e_stors;
    m_result_energy.batm          = e_batm;
    m_result_energy.atm           = e_atm;
    m_result_energy.xbond         = e_xbond;
    m_result_energy.hbond         = e_hbond;

    // GPU ENERGY TERMS (verbosity >= 3)
    if (m_verbosity >= 3) {
        CurcumaLogger::info("=== GPU ENERGY TERMS ===");
        fmt::print("  bond      = {:+.15e}\n", e_bond);
        fmt::print("  angle     = {:+.15e}\n", e_angle);
        fmt::print("  dihedral  = {:+.15e}\n", e_dihedral);
        fmt::print("  inversion = {:+.15e}\n", e_inversion);
        fmt::print("  stors     = {:+.15e}\n", e_stors);
        fmt::print("  batm      = {:+.15e}\n", e_batm);
        fmt::print("  atm       = {:+.15e}\n", e_atm);
        fmt::print("  disp      = {:+.15e}\n", e_disp);
        fmt::print("  brep      = {:+.15e}\n", e_brep);
        fmt::print("  nbrep     = {:+.15e}\n", e_nbrep);
        fmt::print("  coul1     = {:+.15e}\n", e_coul1);
        fmt::print("  coul_self = {:+.15e}\n", e_coul_self);
        fmt::print("  hbond     = {:+.15e}\n", e_hbond);
        fmt::print("  xbond     = {:+.15e}\n", e_xbond);
        CurcumaLogger::info("=== GPU ENERGY END ===");
    }

    return m_result_energy.total() + m_e0;
}

// ============================================================================
// calculate() — original entry point (backward compatible, calls both stages)
// ============================================================================

double FFWorkspaceHip::calculate(bool gradient)
{
    // Backward-compatible: upload charges, then run both stages sequentially.
    // For CPU/GPU overlap, call prepareAndLaunchChargeIndependent() and
    // launchChargeDependentAndFinish() separately from the orchestrator.
    auto& impl = *m_impl;
    const int N = m_natoms;

    // Upload EEQ charges to member (launchChargeDependentAndFinish reads them)
    // setEEQCharges() was already called before calculate() in the old path.

    prepareAndLaunchChargeIndependent(gradient);
    return launchChargeDependentAndFinish(gradient);
}

// ============================================================================
// setGeometry — extract N×3 RowMajor matrix to flat array and upload to GPU
// ============================================================================

void FFWorkspaceHip::setGeometry(const GeoGradMatrix& geom)
{
    auto& impl = *m_impl;
    const int N = m_natoms;
    if (geom.rows() != N) return;

    // Fill pinned staging buffers in-place, then async DMA to GPU (SoA layout).
    // Claude Generated (Phase 10, March 2026): SoA pinned buffers for coalesced warp reads.
    for (int i = 0; i < N; ++i) {
        m_h_x[i] = geom(i, 0);
        m_h_y[i] = geom(i, 1);
        m_h_z[i] = geom(i, 2);
    }
    if (impl.coords.empty()) impl.coords.alloc(N);
    hipMemcpyAsync(impl.coords.d_x.ptr, m_h_x, N*sizeof(double),
                    hipMemcpyHostToDevice, impl.stream);
    hipMemcpyAsync(impl.coords.d_y.ptr, m_h_y, N*sizeof(double),
                    hipMemcpyHostToDevice, impl.stream);
    hipMemcpyAsync(impl.coords.d_z.ptr, m_h_z, N*sizeof(double),
                    hipMemcpyHostToDevice, impl.stream);
}

// ============================================================================
// GPU CN Computation (Phase 1: GPU migration)
// Claude Generated (March 2026): Compute coordination numbers on GPU
// ============================================================================

void FFWorkspaceHip::computeCN(const std::vector<int>& atom_types)
{
    auto& impl = *m_impl;
    const int N = m_natoms;

    // Constants for GFN-FF CN calculation
    // Reference: gfnff_cn.f90:66-126, cn_calculator.cpp:99-155
    constexpr double kn = -7.5;          // CN decay constant
    constexpr double cnmax = 4.4;       // Squashing limit
    constexpr double threshold_factor = 2.5;  // Cutoff factor for rcov
    constexpr double ANG2BOHR = 1.8897259886;

    // Pre-compute threshold squared (largest possible distance)
    // Use max covalent radius * threshold_factor as conservative cutoff
    const auto& rcov_d3 = GFNFFParameters::covalent_rad_d3;  // in Bohr
    double max_rcov = 0.0;
    for (int i = 0; i < N; ++i) {
        int z = atom_types[i];
        if (z >= 1 && z <= static_cast<int>(rcov_d3.size())) {
            max_rcov = std::max(max_rcov, rcov_d3[z - 1]);
        }
    }
    // rcov_d3 is already scaled by 4/3, use conservative threshold
    double threshold_sq = (threshold_factor * 2.0 * max_rcov * ANG2BOHR) * (threshold_factor * 2.0 * max_rcov * ANG2BOHR);
    threshold_sq = 900.0;  // Use 30 Bohr cutoff (same as CNCalculator)

    // Ensure CN buffers are allocated
    if (impl.d_cn_raw.n < N) impl.d_cn_raw.alloc(N);
    if (impl.d_cn_final.n < N) impl.d_cn_final.alloc(N);
    if (impl.d_dlogdcn.n < N) impl.d_dlogdcn.alloc(N);

    // WP3: Pair-list-based CN path — O(n_pairs) instead of O(N²) per step.
    // Pair list built once per topology by generateCNPairListOnGPU().
    // ~5-10x fewer erf() calls for N=1410 (21K pairs vs 2M atoms×atoms).
    if (m_cn_pairs_on_gpu && impl.n_cn_pairs > 0) {
        // Zero cn_raw (atomicAdd kernel requires zeroed input)
        LaunchConfig cfg_zero = getLaunchConfig(N, m_block_size);
        k_zero_double<<<cfg_zero.gridSize, cfg_zero.blockSize, 0, impl.stream>>>(
            impl.d_cn_raw.ptr, N);

        // 1 thread per pair: atomicAdd into cn_raw[i] and cn_raw[j]
        // BOUNDS_LIGHT kernel — must cap block at 256 (May 2026 fix).
        LaunchConfig cfg_pairs = getLaunchConfigLight(impl.n_cn_pairs, m_block_size);
        k_cn_compute_pairs<<<cfg_pairs.gridSize, cfg_pairs.blockSize, 0, impl.stream>>>(
            impl.n_cn_pairs,
            impl.d_cn_idx_i.ptr,
            impl.d_cn_idx_j.ptr,
            impl.d_cn_rcov_sum.ptr,
            impl.coords.d_x.ptr,
            impl.coords.d_y.ptr,
            impl.coords.d_z.ptr,
            impl.d_cn_raw.ptr,
            kn);

        // Log squashing: cn_final[i] = log(1+exp(cnmax)) - log(1+exp(cnmax-cn_raw[i]))
        // BOUNDS_LIGHT kernel — cap block at 256.
        LaunchConfig cfg_log = getLaunchConfigLight(N, m_block_size);
        k_logcn<<<cfg_log.gridSize, cfg_log.blockSize, 0, impl.stream>>>(
            N, impl.d_cn_raw.ptr, impl.d_cn_final.ptr, cnmax);

    } else {
        // Fallback: O(N²) kernel (initial step before pair list is built, or n_pairs==0)
        LaunchConfig cfg_cn = getLaunchConfig(N, m_block_size);
        k_cn_compute<<<cfg_cn.gridSize, cfg_cn.blockSize, 0, impl.stream>>>(
            N,
            impl.coords.d_x.ptr,
            impl.coords.d_y.ptr,
            impl.coords.d_z.ptr,
            impl.d_atom_types.ptr,
            impl.d_cn_raw.ptr,
            impl.d_cn_final.ptr,
            kn,
            cnmax,
            threshold_sq);
    }

    // WP5-C (May 2026): GPU-side D4 dc6dcn skip check.
    // Compare current CN against previous-step snapshot. If max|ΔCN| < 0.01,
    // the next step can skip Gaussian weight + dc6dcn recomputation.
    // Stage 1: per-block max reduction → d_block_max[].
    // Stage 2: single block reduces d_block_max[] → d_dc6dcn_skip[0].
    // Asynchronous D2H copy of skip flag (read by caller in next step).
    {
        // BOUNDS_LIGHT kernels — cap block at 256.
        LaunchConfig cfg_skip = getLaunchConfigLight(N, m_block_size);
        // Zero d_dc6dcn_skip before Stage 2 writes the final flag.
        hipMemsetAsync(impl.d_dc6dcn_skip.ptr, 0, sizeof(int), impl.stream);
        k_check_dc6dcn_skip<<<cfg_skip.gridSize, cfg_skip.blockSize, 0, impl.stream>>>(
            N, impl.d_cn_final.ptr, impl.d_cn_d4_ref.ptr,
            impl.d_dc6dcn_block_max.ptr);
        k_check_dc6dcn_skip_final<<<1, cfg_skip.blockSize, 0, impl.stream>>>(
            cfg_skip.gridSize, impl.d_dc6dcn_block_max.ptr, 0.01,
            impl.d_dc6dcn_skip.ptr);
        hipMemcpyAsync(impl.h_dc6dcn_skip, impl.d_dc6dcn_skip.ptr, sizeof(int),
                        hipMemcpyDeviceToHost, impl.stream);

        // Snapshot current CN for next step's comparison (D2D copy, async)
        hipMemcpyAsync(impl.d_cn_d4_ref.ptr, impl.d_cn_final.ptr,
                        N * sizeof(double), hipMemcpyDeviceToDevice, impl.stream);
    }

    // Compute dlogdcn on GPU: dlogdcn[i] = exp(cnmax)/(exp(cnmax)+exp(cn_raw[i]))
    // Runs on same stream after either CN path — cn_raw is always valid at this point.
    LaunchConfig cfg_dlog = getLaunchConfig(N, m_block_size);
    k_dlogdcn<<<cfg_dlog.gridSize, cfg_dlog.blockSize, 0, impl.stream>>>(
        N, impl.d_cn_raw.ptr, impl.d_dlogdcn.ptr, std::exp(cnmax)
    );

    // Download CN_final, CN_raw, and dlogdcn to pre-allocated pinned buffers
    // Claude Generated (March 2026): Avoids heap corruption from CUDA allocator
    // CN_raw needed for correct dlogdcn computation (chain-rule gradient)
    hipMemcpyAsync(m_h_cn_final, impl.d_cn_final.ptr, N * sizeof(double),
                    hipMemcpyDeviceToHost, impl.stream);
    hipMemcpyAsync(m_h_cn_raw, impl.d_cn_raw.ptr, N * sizeof(double),
                    hipMemcpyDeviceToHost, impl.stream);
    hipMemcpyAsync(m_h_dlogdcn, impl.d_dlogdcn.ptr, N * sizeof(double),
                    hipMemcpyDeviceToHost, impl.stream);
    // G-P1 (Apr 2026): Sync deferred — caller calls finalizeCNForCPU() after Phase 4 launch.
    // GPU finishes CN kernel in ~0.15ms; by the time finalizeCNForCPU() is reached (~0.3ms
    // later after Phase 4 launch overhead), the stream is already complete → no sleep penalty.
    m_cn_computed = false;
}

// G-P1 (Apr 2026): Finalize deferred CN download — sync main stream, copy to caller's buffer.
// Must be called after Phase 4 launch so the GPU CN kernel completes during launch overhead.
void FFWorkspaceHip::finalizeCNForCPU(Vector& out_cn)
{
    const int N = m_natoms;
    hipStreamSynchronize(m_impl->stream);   // returns near-instantly if GPU already done
    if (out_cn.size() != N) out_cn.resize(N);
    std::memcpy(out_cn.data(), m_h_cn_final, N * sizeof(double));
    m_cn_computed = true;
}

// ============================================================================
// Getters
// ============================================================================

const Matrix& FFWorkspaceHip::gradient() const
{
    // Debug: Check if result is properly sized
    if (m_result_gradient.rows() != m_natoms || m_result_gradient.cols() != 3) {
        std::cout << "[DEBUG FFWorkspaceHip::gradient] WARNING: gradient size mismatch! "
                  << "rows=" << m_result_gradient.rows() << " cols=" << m_result_gradient.cols()
                  << " expected (" << m_natoms << ", 3)\n" << std::flush;
    }
    return m_result_gradient;
}

const FFEnergyComponents& FFWorkspaceHip::energyComponents() const
{
    return m_result_energy;
}

const Vector& FFWorkspaceHip::dEdcnTotal() const
{
    return m_dEdcn_total;
}

int FFWorkspaceHip::dispersionCount() const
{
    return m_impl ? m_impl->disp.n : 0;
}

int FFWorkspaceHip::bondCount() const
{
    return m_impl ? m_impl->bonds.n : 0;
}

// DEPRECATED (Phase 10, March 2026): coords now in SoA layout. Use getDeviceXPtr/YPtr/ZPtr().
const double* FFWorkspaceHip::getDeviceCoordsPtr() const
{
    return nullptr;  // Deprecated: coords now in SoA layout. Use getDeviceXPtr/YPtr/ZPtr().
}

const double* FFWorkspaceHip::getDeviceXPtr() const
{
    return m_impl ? m_impl->coords.d_x.ptr : nullptr;
}
const double* FFWorkspaceHip::getDeviceYPtr() const
{
    return m_impl ? m_impl->coords.d_y.ptr : nullptr;
}
const double* FFWorkspaceHip::getDeviceZPtr() const
{
    return m_impl ? m_impl->coords.d_z.ptr : nullptr;
}

// Claude Generated (March 2026): Explicit sync for cross-stream safety
void FFWorkspaceHip::synchronizeMainStream()
{
    if (m_impl && m_impl->stream)
        checkCuda(hipStreamSynchronize(m_impl->stream), "synchronizeMainStream");
}

// ---------------------------------------------------------------------------
// WP2: uploadEEQTopologyParams + device pointer getters
// Claude Generated (May 2026): Upload topology-constant EEQ parameters once
// per topology build so k_build_eeq_rhs can construct rhs_atoms on GPU.
// ---------------------------------------------------------------------------

void FFWorkspaceHip::uploadEEQTopologyParams(const GFNFF::EEQGPUParams& params)
{
    if (!m_impl) return;
    auto& et = m_impl->eeq_topo;
    const int N = m_natoms;

    et.d_alpha_corrected.upload(params.alpha_corrected.data(), N);
    et.d_gam_corrected.upload(params.gam_corrected.data(), N);
    et.d_chi_corrected.upload(params.chi_corrected_static.data(), N);
    et.d_cnf.upload(params.cnf.data(), N);

    const int nfrag = params.nfrag;
    if (nfrag > 0)
        et.d_rhs_constraints.upload(params.rhs_constraints.data(), nfrag);
    et.nfrag = nfrag;

    // Allocate output buffer for k_build_eeq_rhs (written per-step, not uploaded here)
    if (et.d_rhs_atoms.n < N)
        et.d_rhs_atoms.alloc(N);

    et.valid = true;
}

bool FFWorkspaceHip::isEEQTopoValid() const
{
    return m_impl && m_impl->eeq_topo.valid;
}

const double* FFWorkspaceHip::getDeviceRHSPtr() const
{
    return (m_impl && m_impl->eeq_topo.valid) ? m_impl->eeq_topo.d_rhs_atoms.ptr : nullptr;
}

const double* FFWorkspaceHip::getDeviceRHSConstraintsPtr() const
{
    return (m_impl && m_impl->eeq_topo.valid) ? m_impl->eeq_topo.d_rhs_constraints.ptr : nullptr;
}

const double* FFWorkspaceHip::getDeviceAlphaPtr() const
{
    return (m_impl && m_impl->eeq_topo.valid) ? m_impl->eeq_topo.d_alpha_corrected.ptr : nullptr;
}

const double* FFWorkspaceHip::getDeviceGamPtr() const
{
    return (m_impl && m_impl->eeq_topo.valid) ? m_impl->eeq_topo.d_gam_corrected.ptr : nullptr;
}

// ---------------------------------------------------------------------------
// setMixedPrecision — enable/disable FP32 intermediates for repulsion/batm/xbonds
// Claude Generated (March 2026)
// ---------------------------------------------------------------------------
void FFWorkspaceHip::setMixedPrecision(bool enable)
{
    if (m_impl) m_impl->use_mixed_precision = enable;
}

// ---------------------------------------------------------------------------
// checkDisplacement — GPU flag-based topology displacement check
// Claude Generated (March 2026)
// ---------------------------------------------------------------------------
bool FFWorkspaceHip::checkDisplacement(double threshold)
{
    if (!m_impl) return true;
    auto& impl = *m_impl;
    if (impl.ref_coords.empty() || impl.coords.empty()) return true;

    // Zero the flag
    int zero = 0;
    checkCuda(hipMemcpyAsync(impl.d_disp_flag.ptr, &zero, sizeof(int),
                               hipMemcpyHostToDevice, impl.stream),
              "checkDisplacement: zero flag");

    // Launch displacement check kernel
    LaunchConfig cfg = getLaunchConfig(m_natoms, m_block_size);
    k_check_displacement<<<cfg.gridSize, cfg.blockSize, 0, impl.stream>>>(
        m_natoms,
        impl.coords.d_x.ptr,
        impl.coords.d_y.ptr,
        impl.coords.d_z.ptr,
        impl.ref_coords.d_rx.ptr,
        impl.ref_coords.d_ry.ptr,
        impl.ref_coords.d_rz.ptr,
        threshold * threshold,
        impl.d_disp_flag.ptr);

    // Download single int result
    int result = 0;
    checkCuda(hipMemcpyAsync(&result, impl.d_disp_flag.ptr, sizeof(int),
                               hipMemcpyDeviceToHost, impl.stream),
              "checkDisplacement: download flag");
    checkCuda(hipStreamSynchronize(impl.stream), "checkDisplacement: sync");

    return result != 0;
}

// ---------------------------------------------------------------------------
// updateReferenceGeometry — copy coords SoA → ref_coords SoA (device-to-device)
// Claude Generated (March 2026)
// ---------------------------------------------------------------------------
void FFWorkspaceHip::updateReferenceGeometry()
{
    if (!m_impl || m_impl->coords.empty()) return;
    auto& impl = *m_impl;

    // Allocate ref_coords on first call
    if (impl.ref_coords.empty()) {
        impl.ref_coords.alloc(m_natoms);
    }

    checkCuda(hipMemcpyAsync(impl.ref_coords.d_rx.ptr, impl.coords.d_x.ptr,
                               m_natoms * sizeof(double),
                               hipMemcpyDeviceToDevice, impl.stream),
              "updateReferenceGeometry: D2D copy x");
    checkCuda(hipMemcpyAsync(impl.ref_coords.d_ry.ptr, impl.coords.d_y.ptr,
                               m_natoms * sizeof(double),
                               hipMemcpyDeviceToDevice, impl.stream),
              "updateReferenceGeometry: D2D copy y");
    checkCuda(hipMemcpyAsync(impl.ref_coords.d_rz.ptr, impl.coords.d_z.ptr,
                               m_natoms * sizeof(double),
                               hipMemcpyDeviceToDevice, impl.stream),
              "updateReferenceGeometry: D2D copy z");
}

// ===================== gpu utils (cuda/gpu_utils.cpp) ======================
/*
 * <GPU Memory Utilities Implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include "gpu_utils_hip.h"
#include "src/core/curcuma_logger.h"

#include <hip/hip_runtime.h>
#include <cmath>
#include <sstream>
#include <iomanip>

namespace GPUUtils {

bool getGPUMemoryInfo(size_t& free, size_t& total) {
    hipError_t err = hipMemGetInfo(&free, &total);
    if (err != hipSuccess) {
        CurcumaLogger::warn(std::string("hipMemGetInfo failed: ") + hipGetErrorString(err));
        return false;
    }
    return true;
}

size_t estimateGFNFFGPUMemory(int natoms) {
    // SoA buffers (see ff_workspace_gpu.cu):
    // - coordinates: N * 3 doubles
    // - gradient: N * 3 doubles
    // - charges: N doubles
    // - CN: N doubles
    // - dEdcn: N doubles
    // - EEQ matrix: N * N doubles (only during solve, but allocate once)
    // - Pair lists: ~100 * N integers (dispersion, Coulomb, HB, XB pairs)
    // - Temporary buffers: ~50 * N doubles

    // Base overhead for CUDA context, kernels, etc.
    size_t base_overhead = 2 * 1024 * 1024;  // 2 MB for kernel overhead

    // Per-atom memory (doubles)
    size_t per_atom_doubles = 3 + 3 + 1 + 1 + 1 + 50;  // coords, gradient, charges, CN, dEdcn, temps

    // Per-atom memory (integers for pair lists)
    size_t per_atom_ints = 100;  // conservative estimate for all pair types

    size_t memory = base_overhead;
    memory += natoms * per_atom_doubles * sizeof(double);
    memory += natoms * natoms * sizeof(double);  // EEQ matrix (worst case)
    memory += natoms * per_atom_ints * sizeof(int);

    // Add 20% safety margin
    memory = static_cast<size_t>(memory * 1.2);

    return memory;
}

bool checkGPUMemoryAvailable(int natoms, double threshold) {
    size_t free, total;
    if (!getGPUMemoryInfo(free, total)) {
        return false;  // CUDA not available or error
    }

    size_t required = estimateGFNFFGPUMemory(natoms);
    size_t available = static_cast<size_t>(free * threshold);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("=== GPU Memory Check ===");
        CurcumaLogger::param("Required", formatMemorySize(required));
        CurcumaLogger::param("Available", formatMemorySize(available));
        CurcumaLogger::param("Free", formatMemorySize(free));
        CurcumaLogger::param("Total", formatMemorySize(total));
    }

    return required < available;
}

std::string formatMemorySize(size_t bytes) {
    const char* units[] = {"B", "KB", "MB", "GB", "TB"};
    int unit_idx = 0;
    double size = static_cast<double>(bytes);

    while (size >= 1024.0 && unit_idx < 4) {
        size /= 1024.0;
        unit_idx++;
    }

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << size << " " << units[unit_idx];
    return oss.str();
}

void logGPUMemoryStatus(int verbosity) {
    if (CurcumaLogger::get_verbosity() < verbosity) {
        return;
    }

    size_t free, total;
    if (!getGPUMemoryInfo(free, total)) {
        CurcumaLogger::warn("GPU not available for memory status query");
        return;
    }

    size_t used = total - free;

    CurcumaLogger::info("=== GPU Memory Status ===");
    CurcumaLogger::param("Total", formatMemorySize(total));
    CurcumaLogger::param("Used", formatMemorySize(used));
    CurcumaLogger::param("Free", formatMemorySize(free));

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << (100.0 * used / total) << "%";
    CurcumaLogger::param("Utilization", oss.str());
}

} // namespace GPUUtils



// ===================== EEQ solver (rocSOLVER/rocBLAS port, M2) =================
#include "eeq_solver_hip.hiph"

#endif // USE_ROCM_GFNFF
