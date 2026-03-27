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
 * Coordinate layout: coords[3*i + 0] = x_i, [3*i+1] = y_i, [3*i+2] = z_i
 */

#include "gfnff_kernels.cuh"
#include <cuda_runtime.h>
#include <math.h>
#include <float.h>

// ============================================================================
// Device geometry helpers
// ============================================================================

/// Squared distance between atoms i and j
__device__ __forceinline__ double dist_sq(const double* __restrict__ c, int i, int j)
{
    double dx = c[3*j]   - c[3*i];
    double dy = c[3*j+1] - c[3*i+1];
    double dz = c[3*j+2] - c[3*i+2];
    return dx*dx + dy*dy + dz*dz;
}

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
    #pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1)
        val += __shfl_down_sync(0xFFFFFFFF, val, offset);
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
    const double* __restrict__ current_coords,
    const double* __restrict__ ref_coords,
    double        threshold_sq,
    int*          exceeded_flag)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    double dx = current_coords[3*i]   - ref_coords[3*i];
    double dy = current_coords[3*i+1] - ref_coords[3*i+1];
    double dz = current_coords[3*i+2] - ref_coords[3*i+2];
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
    const double* __restrict__ coords,
    double*                    dEdcn,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid];
        double dx = coords[3*i]   - coords[3*j];
        double dy = coords[3*i+1] - coords[3*j+1];
        double dz = coords[3*i+2] - coords[3*j+2];
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

            // dEdcn chain-rule: dc6/dcn contribution for CN gradient
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
    const double* __restrict__ coords,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid];
        double dx = coords[3*i]   - coords[3*j];
        double dy = coords[3*i+1] - coords[3*j+1];
        double dz = coords[3*i+2] - coords[3*j+2];
        double r2  = dx*dx + dy*dy + dz*dz;
        double rij = sqrt(r2);

        if (rij <= r_cut[tid] && rij >= 1e-8) {
            double r1_5      = rij * sqrt(rij);
            double alp_r     = alpha[tid] * r1_5;
            if (!isnan(alp_r) && !isnan(repab[tid]) && alp_r <= 700.0) {
                double exp_term  = exp(-alp_r);
                double base_E    = repab[tid] * exp_term / rij;
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
    const double* __restrict__ coords,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid];
        float dx = (float)(coords[3*i]   - coords[3*j]);
        float dy = (float)(coords[3*i+1] - coords[3*j+1]);
        float dz = (float)(coords[3*i+2] - coords[3*j+2]);
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
                add_grad(grad, i,  (double)( fac*dx), (double)( fac*dy), (double)( fac*dz));
                add_grad(grad, j,  (double)(-fac*dx), (double)(-fac*dy), (double)(-fac*dz));
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
    const double* __restrict__ coords,
    const double* __restrict__ charges,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid];
        double dx  = coords[3*i]   - coords[3*j];
        double dy  = coords[3*i+1] - coords[3*j+1];
        double dz  = coords[3*i+2] - coords[3*j+2];
        double r2  = dx*dx + dy*dy + dz*dz;
        double rij = sqrt(r2);

        if (rij <= r_cut[tid] && rij >= 1e-10) {
            double qi = charges[i];
            double qj = charges[j];
            if (!isnan(qi) && !isnan(qj)) {
                double gamma_r = gamma_ij[tid] * rij;
                double erf_v   = erf(gamma_r);
                local_E = qi * qj * erf_v / rij;

                // Gradient
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
    const double* __restrict__ coords,
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
        double dx  = coords[3*i]   - coords[3*j];
        double dy  = coords[3*i+1] - coords[3*j+1];
        double dz  = coords[3*i+2] - coords[3*j+2];
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
// Device helper: angle between vectors (i-j) and (k-j), returns cos(theta)
// Also fills d_cos[3]: ∂cos/∂x_i, ∂cos/∂x_k (x-component only as template)
// ============================================================================
__device__ double cos_angle_grad(
    const double* __restrict__ c, int i, int j, int k,
    double& d_ri_x, double& d_ri_y, double& d_ri_z,
    double& d_rk_x, double& d_rk_y, double& d_rk_z,
    double& d_rj_x, double& d_rj_y, double& d_rj_z,
    bool grad)
{
    double eij_x = c[3*i]   - c[3*j];
    double eij_y = c[3*i+1] - c[3*j+1];
    double eij_z = c[3*i+2] - c[3*j+2];
    double ekj_x = c[3*k]   - c[3*j];
    double ekj_y = c[3*k+1] - c[3*j+1];
    double ekj_z = c[3*k+2] - c[3*j+2];

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
    const double* __restrict__ coords,
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

            double dij_x = coords[3*i]   - coords[3*j];
            double dij_y = coords[3*i+1] - coords[3*j+1];
            double dij_z = coords[3*i+2] - coords[3*j+2];
            double dkj_x = coords[3*k]   - coords[3*j];
            double dkj_y = coords[3*k+1] - coords[3*j+1];
            double dkj_z = coords[3*k+2] - coords[3*j+2];

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
            double cos_a = cos_angle_grad(coords, i, j, k,
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
    const double* __restrict__ c,
    int i, int j, int k, int l,
    double derivate[4][3],   // [atom][xyz]
    bool do_grad)
{
    // b1 = rj - ri, b2 = rk - rj, b3 = rl - rk
    double b1x = c[3*j]   - c[3*i];
    double b1y = c[3*j+1] - c[3*i+1];
    double b1z = c[3*j+2] - c[3*i+2];
    double b2x = c[3*k]   - c[3*j];
    double b2y = c[3*k+1] - c[3*j+1];
    double b2z = c[3*k+2] - c[3*j+2];
    double b3x = c[3*l]   - c[3*k];
    double b3y = c[3*l+1] - c[3*k+1];
    double b3z = c[3*l+2] - c[3*k+2];

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
    const double* __restrict__ coords,
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
                double dx=coords[3*b]-coords[3*a];
                double dy=coords[3*b+1]-coords[3*a+1];
                double dz=coords[3*b+2]-coords[3*a+2];
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
            double phi = dihedral_angle_grad(coords, ai, aj, ak, al, derivate, true);

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
                double dx=coords[3*ai]-coords[3*aj], dy=coords[3*ai+1]-coords[3*aj+1], dz=coords[3*ai+2]-coords[3*aj+2];
                double fac = E_raw * d2_ij * d_jk * d_kl;
                add_grad(grad, ai,  fac*dx,  fac*dy,  fac*dz);
                add_grad(grad, aj, -fac*dx, -fac*dy, -fac*dz);
            }
            {
                double dx=coords[3*aj]-coords[3*ak], dy=coords[3*aj+1]-coords[3*ak+1], dz=coords[3*aj+2]-coords[3*ak+2];
                double fac = E_raw * d_ij * d2_jk * d_kl;
                add_grad(grad, aj,  fac*dx,  fac*dy,  fac*dz);
                add_grad(grad, ak, -fac*dx, -fac*dy, -fac*dz);
            }
            {
                double dx=coords[3*ak]-coords[3*al], dy=coords[3*ak+1]-coords[3*al+1], dz=coords[3*ak+2]-coords[3*al+2];
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
    const double* __restrict__ c,
    int i, int j, int k, int l,
    double derivate[4][3],
    bool do_grad)
{
    // Fortran convention vectors (gfnff_helpers.f90:436-440, 465-471):
    //   re = r_i - r_j  (center - nb1)
    //   rd = r_k - r_j  (nb2 - nb1)
    //   rv = r_l - r_i  (nb3 - center)
    double re[3], rd[3], rv[3];
    for (int d = 0; d < 3; ++d) {
        re[d] = c[3*i+d] - c[3*j+d];
        rd[d] = c[3*k+d] - c[3*j+d];
        rv[d] = c[3*l+d] - c[3*i+d];
    }

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
    const double* __restrict__ coords,
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
            double dx=coords[3*b]-coords[3*a];
            double dy=coords[3*b+1]-coords[3*a+1];
            double dz=coords[3*b+2]-coords[3*a+2];
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
        double omega = inversion_angle_grad(coords, i, j, k, l, derivate, true);
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

        {
            double dx=coords[3*i]-coords[3*j], dy=coords[3*i+1]-coords[3*j+1], dz=coords[3*i+2]-coords[3*j+2];
            double fac = E_raw * dd_ij * damp_jk * damp_jl;
            add_grad(grad, j,  fac*dx,  fac*dy,  fac*dz);
            add_grad(grad, i, -fac*dx, -fac*dy, -fac*dz);
        }
        {
            double dx=coords[3*k]-coords[3*j], dy=coords[3*k+1]-coords[3*j+1], dz=coords[3*k+2]-coords[3*j+2];
            double fac = E_raw * damp_ij * dd_jk * damp_jl;
            add_grad(grad, j,  fac*dx,  fac*dy,  fac*dz);
            add_grad(grad, k, -fac*dx, -fac*dy, -fac*dz);
        }
        {
            double dx=coords[3*l]-coords[3*j], dy=coords[3*l+1]-coords[3*j+1], dz=coords[3*l+2]-coords[3*j+2];
            double fac = E_raw * damp_ij * damp_jk * dd_jl;
            add_grad(grad, j,  fac*dx,  fac*dy,  fac*dz);
            add_grad(grad, l, -fac*dx, -fac*dy, -fac*dz);
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Upload covalent radii table to GPU constant memory
// Called once during FFWorkspaceGPU initialization
// ============================================================================
void upload_rcov_d3(const double* rcov, int n)
{
    int count = (n < 87) ? n : 87;
    cudaMemcpyToSymbol(d_rcov_d3, rcov, count * sizeof(double));
}

// ============================================================================
// Covalent radii for HB/XB vdW radii lookup (Å, element 1..86)
// Source: gfnff_par.h GFNFFParameters::covalent_radii
// ============================================================================
__constant__ double d_cov_radii[100];  // uploaded once at init

void upload_covalent_radii(const double* radii, int n)
{
    int count = (n < 100) ? n : 100;
    cudaMemcpyToSymbol(d_cov_radii, radii, count * sizeof(double));
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
    const double* __restrict__ coords,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int ai = idx_i[tid], aj = idx_j[tid], ak = idx_k[tid], al = idx_l[tid];

        double derivate[4][3];
        double phi = dihedral_angle_grad(coords, ai, aj, ak, al, derivate, true);

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
    const double* __restrict__ coords,
    const double* __restrict__ topo_charges,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid], k = idx_k[tid];

        double ix = coords[3*i], iy = coords[3*i+1], iz = coords[3*i+2];
        double jx = coords[3*j], jy = coords[3*j+1], jz = coords[3*j+2];
        double kx = coords[3*k], ky = coords[3*k+1], kz = coords[3*k+2];

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
    const double* __restrict__ coords,
    const double* __restrict__ topo_charges,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid], k = idx_k[tid];

        float ix = (float)coords[3*i], iy = (float)coords[3*i+1], iz = (float)coords[3*i+2];
        float jx = (float)coords[3*j], jy = (float)coords[3*j+1], jz = (float)coords[3*j+2];
        float kx = (float)coords[3*k], ky = (float)coords[3*k+1], kz = (float)coords[3*k+2];

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
    const double* __restrict__ coords,
    double*                    grad,
    double*                    energy)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid], k = idx_k[tid];

        double ix = coords[3*i], iy = coords[3*i+1], iz = coords[3*i+2];
        double jx = coords[3*j], jy = coords[3*j+1], jz = coords[3*j+2];
        double kx = coords[3*k], ky = coords[3*k+1], kz = coords[3*k+2];

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
    const double* __restrict__ coords,
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

        double ax = coords[3*A], ay = coords[3*A+1], az = coords[3*A+2];
        double xx = coords[3*X], xy = coords[3*X+1], xz = coords[3*X+2];
        double bx = coords[3*B], by = coords[3*B+1], bz = coords[3*B+2];

        double rAX_x = xx-ax, rAX_y = xy-ay, rAX_z = xz-az;
        double rXB_x = bx-xx, rXB_y = by-xy, rXB_z = bz-xz;
        double rAB_x = bx-ax, rAB_y = by-ay, rAB_z = bz-az;

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
    const double* __restrict__ coords,
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

        float ax = (float)coords[3*A], ay = (float)coords[3*A+1], az = (float)coords[3*A+2];
        float xx = (float)coords[3*X], xy = (float)coords[3*X+1], xz = (float)coords[3*X+2];
        float bx = (float)coords[3*B], by = (float)coords[3*B+1], bz = (float)coords[3*B+2];

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
    const double* __restrict__ coords,
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

    double ax = coords[3*A], ay = coords[3*A+1], az = coords[3*A+2];
    double hx = coords[3*H], hy = coords[3*H+1], hz = coords[3*H+2];
    double bx = coords[3*B], by = coords[3*B+1], bz = coords[3*B+2];

    double rAH_x = hx-ax, rAH_y = hy-ay, rAH_z = hz-az;
    double rHB_x = bx-hx, rHB_y = by-hy, rHB_z = bz-hz;
    double rAB_x = bx-ax, rAB_y = by-ay, rAB_z = bz-az;

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
            double nbx = coords[3*nb], nby = coords[3*nb+1], nbz = coords[3*nb+2];
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
                lp_vx_c4 += coords[3*nb]   - bx;
                lp_vy_c4 += coords[3*nb+1] - by;
                lp_vz_c4 += coords[3*nb+2] - bz;
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
            double cx = coords[3*C_idx], cy = coords[3*C_idx+1], cz = coords[3*C_idx+2];
            double vab_x = cx-bx, vab_y = cy-by, vab_z = cz-bz;
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
            double phi = dihedral_angle_grad(coords, D_idx, B_idx, C_idx, H_idx, derivate, false);

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
                double nbx = coords[3*nb], nby = coords[3*nb+1], nbz = coords[3*nb+2];
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
            double cx = coords[3*C_idx], cy = coords[3*C_idx+1], cz = coords[3*C_idx+2];
            double vab_x = cx-bx, vab_y = cy-by, vab_z = cz-bz;
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
                double phi = dihedral_angle_grad(coords, D_idx, B, C_idx, H, derivate, true);

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
    const double* __restrict__ coords,
    const double* __restrict__ dEdcn,
    const double* __restrict__ dlogdcn,
    double*                    grad,
    double        kn)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pairs) return;

    int i = idx_i[tid], j = idx_j[tid];
    double dx = coords[3*i]   - coords[3*j];
    double dy = coords[3*i+1] - coords[3*j+1];
    double dz = coords[3*i+2] - coords[3*j+2];
    double r2 = dx*dx + dy*dy + dz*dz;
    double rij = sqrt(r2);
    if (rij < 1e-10) return;

    double rcov = rcov_sum[tid];
    double dr = (rij - rcov) / rcov;

    // dS/dr = (kn / sqrt(pi)) * exp(-kn^2 * dr^2) / rcov
    static const double inv_sqrtpi = 0.5641895835477563;  // 1/sqrt(pi)
    double dSdr = kn * inv_sqrtpi * exp(-kn * kn * dr * dr) / rcov;

    // Combined chain-rule factor
    double fac = dSdr / rij * (dEdcn[i] * dlogdcn[i] + dEdcn[j] * dlogdcn[j]);

    add_grad(grad, i,  fac*dx,  fac*dy,  fac*dz);
    add_grad(grad, j, -fac*dx, -fac*dy, -fac*dz);
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
    const double* __restrict__ coords,
    const double* __restrict__ zz_hb,
    double*                    grad,
    double        kn)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pairs) return;

    int H = idx_H[tid], B = idx_B[tid];
    double dx = coords[3*H]   - coords[3*B];
    double dy = coords[3*H+1] - coords[3*B+1];
    double dz = coords[3*H+2] - coords[3*B+2];
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
    const double* __restrict__ coords,    // [3*N] in Bohr
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
    double xi = coords[3*i];
    double yi = coords[3*i + 1];
    double zi_coord = coords[3*i + 2];

    double cn_sum = 0.0;

    // Loop over all other atoms
    for (int j = 0; j < natoms; ++j) {
        if (i == j) continue;

        int zj = atom_types[j];
        if (zj < 1 || zj > 86) continue;

        double dx = xi - coords[3*j];
        double dy = yi - coords[3*j + 1];
        double dz = zi_coord - coords[3*j + 2];
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
    const double* __restrict__ refcn,
    const int*    __restrict__ refn,
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

    int nref = refn[elem];
    if (nref > D4_MAX_REF) nref = D4_MAX_REF;

    double cn_i = cn[i];
    int refcn_base = elem * D4_MAX_REF;

    // Phase 1: Compute raw exponential weights and sums
    double expw[D4_MAX_REF];
    double dexpw[D4_MAX_REF];
    double norm = 0.0;
    double dnorm = 0.0;

    for (int r = 0; r < nref; ++r) {
        double diff = cn_i - refcn[refcn_base + r];
        double ew = exp(-wf * diff * diff);
        double dew = 2.0 * wf * (refcn[refcn_base + r] - cn_i) * ew;
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
    const int*    __restrict__ refn,
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

    int nri = refn[ei];
    int nrj = refn[ej];

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
