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

#include <thrust/execution_policy.h>
#include <thrust/scan.h>
#include "gfnff_kernels.cuh"
#include <cuda_runtime.h>
#include <math.h>
#include <float.h>

// ============================================================================
// Device geometry helpers
// ============================================================================

// ============================================================================
// Order-independent accumulation (Claude Generated, Oct 2026; docs/GPU_REPRODUCIBILITY.md)
// FP64 atomicAdd depends on the arrival order of the threads, so two runs differ in the last bits and
// the difference is amplified by the dynamics. Integer addition is associative: the contribution is
// rounded to a multiple of 2^-FIX_BITS, added to the same 8 bytes as a 64-bit integer (two's
// complement wraps harmlessly, only the final sum has to fit) and converted back in place by
// k_fix_to_double once the accumulation is over. The buffer must be zeroed first (0 bits = 0 and 0.0).
// ============================================================================
/// Claude Generated (Oct 2026, CUDA/ROCm parity item 12): set (sticky) when a contribution was not finite or alone
/// outside the accumulator range; read by the host after the step (FFWorkspaceGPU, gfnffFixBadReadAsync), which then
/// returns NaN as the FP64 path would. Before, a NaN/Inf contribution became llrint(NaN) = LLONG_MIN (a silent wrong
/// number) and an Inf force looped for ever in the chunk loop.
__device__ int d_gfnff_fix_bad;

/// llrint(x * scale) with the range check of the ROCm port: outside the int64 range (or not finite) the flag is set
/// and 0 is returned. Identical to llrint for every value that fits.
__device__ __forceinline__ long long fix_round(double x, double scale)
{
    const double y = x * scale;
    if (!(fabs(y) < 9.2e18)) {
        atomicOr(&d_gfnff_fix_bad, 1);
        return 0ll;
    }
    return llrint(y);
}

__device__ __forceinline__ void fix_add_u64(unsigned long long* a, double v, double scale)
{
    const double av = fabs(v);
    if (!(av * scale < 9.2e18)) {                      // beyond the int64 range (or NaN/Inf): the sum is wrong anyway
        atomicOr(&d_gfnff_fix_bad, 1);
        if (!(av < 5.7e17)) return;                    // not finite or beyond 2^59: add nothing usable
    }
    if (av < 64.0) {                                   // the common case: one rounded value
        atomicAdd(a, static_cast<unsigned long long>(llrint(v * scale)));
        return;
    }
    // |v| >= 64: the former loop added copysign(64, v) n times and then the remainder. The n chunks are folded into one
    // integer, n * llrint(64 scale) + llrint(r scale) with r = v - n * 64 (exact): modulo 2^64 the same accumulator
    // bits as the loop, with one atomic and a bounded cost (the loop was unbounded for a huge v).
    const double n = trunc(av * (1.0 / 64.0));
    const double r = v - copysign(n * 64.0, v);
    const unsigned long long chunk = static_cast<unsigned long long>(llrint(copysign(64.0, v) * scale));
    atomicAdd(a, static_cast<unsigned long long>(n) * chunk + static_cast<unsigned long long>(llrint(r * scale)));
}

__device__ __forceinline__ void fix_add(double* addr, double v, double scale)
{
    fix_add_u64(reinterpret_cast<unsigned long long*>(addr), v, scale);
}

/// Copy of a fixed-point accumulator as doubles (the accumulator itself stays untouched; diagnostic snapshots).
__global__ void k_fix_copy(int n, const double* __restrict__ src, double* __restrict__ dst, double inv_scale)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    dst[i] = static_cast<double>(reinterpret_cast<const long long*>(src)[i]) * inv_scale;
}

/// In-place conversion of a fixed-point accumulator (see fix_add) into doubles.
__global__ void k_fix_to_double(int n, double* __restrict__ buf, double inv_scale)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const long long x = reinterpret_cast<long long*>(buf)[i];
    buf[i] = static_cast<double>(x) * inv_scale;
}

/// Add vector contribution to gradient (atomicAdd on x,y,z)
__device__ __forceinline__ void add_grad(double* __restrict__ grad, int atom,
                                          double fx, double fy, double fz)
{
    // Fixed-point (order-independent) accumulation: grad is an integer accumulator until
    // k_fix_to_double converts it at the end of the step (FFWorkspaceGPU::calculate*, docs/GPU_REPRODUCIBILITY.md).
    fix_add(&grad[3*atom],   fx, GFNFF_FIX_FORCE);
    fix_add(&grad[3*atom+1], fy, GFNFF_FIX_FORCE);
    fix_add(&grad[3*atom+2], fz, GFNFF_FIX_FORCE);
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

// Integer version (Oct 2026): thread values are fixed-point integers (multiples of 2^-48 Eh, dispersion 2^-52), so the warp, block and
// atomic sums are associative and do not depend on which pair or atom a thread got - the energy is the same bit for
// bit however an atomically built list was ordered. The slot is converted by k_fix_to_double at the end of the step.
__device__ __forceinline__ void blockReduceAddEnergyFix(long long local_E, double* energy)
{
    long long warp_sum = local_E;
    for (int off = 16; off > 0; off >>= 1) warp_sum += __shfl_down_sync(0xffffffffu, warp_sum, off);

    __shared__ long long warp_sums_fx[32];
    int lane = threadIdx.x & 31;
    int warp_id = threadIdx.x >> 5;

    if (lane == 0)
        warp_sums_fx[warp_id] = warp_sum;
    __syncthreads();

    int num_warps = (blockDim.x + 31) >> 5;
    if (warp_id == 0) {
        long long sum = (lane < num_warps) ? warp_sums_fx[lane] : 0ll;
        for (int off = 16; off > 0; off >>= 1) sum += __shfl_down_sync(0xffffffffu, sum, off);
        if (lane == 0)
            atomicAdd(reinterpret_cast<unsigned long long*>(energy), static_cast<unsigned long long>(sum));
    }
}

__device__ __forceinline__ void blockReduceAddEnergy(double local_E, double* energy)
{
    blockReduceAddEnergyFix(fix_round(local_E, GFNFF_FIX_ENERGY), energy);
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

// Claude Generated (Nov 2026, WP-PERIODIC-NONBONDED S1): local image frame of a bonded / HB / XB term, the
// arithmetic of the CPU placeMinimumImage() / imageNear() (ff_workspace_gfnff.cpp): the image of atom `atom`
// next to the reference point (rx, ry, rz) is ref + mic(x_atom - ref). Without a cell the position is copied
// as it is (no "ref + (x - ref)" round trip), so the non-periodic kernels keep their bits.
__device__ __forceinline__ void micPlace(const double* __restrict__ cx, const double* __restrict__ cy,
                                         const double* __restrict__ cz, int atom,
                                         double rx, double ry, double rz, double& x, double& y, double& z)
{
    if (!d_has_pbc) { x = cx[atom]; y = cy[atom]; z = cz[atom]; return; }
    double dx = cx[atom] - rx, dy = cy[atom] - ry, dz = cz[atom] - rz;
    applyMIC(dx, dy, dz);
    x = rx + dx; y = ry + dy; z = rz + dz;
}

/// Four-atom chain a-b-c-d in slots 0..3: b at its position, a and c next to b, d next to c (placeChainOfFour()).
__device__ __forceinline__ void micChainOfFour(const double* __restrict__ cx, const double* __restrict__ cy,
                                               const double* __restrict__ cz, int a, int b, int c, int d,
                                               double* lx, double* ly, double* lz)
{
    lx[1] = cx[b]; ly[1] = cy[b]; lz[1] = cz[b];
    micPlace(cx, cy, cz, a, lx[1], ly[1], lz[1], lx[0], ly[0], lz[0]);
    micPlace(cx, cy, cz, c, lx[1], ly[1], lz[1], lx[2], ly[2], lz[2]);
    micPlace(cx, cy, cz, d, lx[2], ly[2], lz[2], lx[3], ly[3], lz[3]);
}

// ============================================================================
// Atomic virial on the device (Claude Generated, Nov 2026, GPU_PERIODIC_GAPS A1; definition at FFVirialComponents)
// Xi_ab -= p_a g_b for every gradient row g a term writes, p the position of that atom in the frame the term used,
// taken relative to one atom of the term (the rows of a term sum to zero, so the origin drops out). Row-major
// 3a+b in a per-thread array; blockReduceAddVirial rounds each of the nine values to a multiple of 2^-48 Eh and sums
// them as integers (order-independent, docs/GPU_REPRODUCIBILITY.md) into the nine slots of the term. Only the
// VIRIAL instantiations of the kernels and the separate *_virial kernels call these; the other kernels are unchanged.
// ============================================================================
__device__ __forceinline__ void vir_outer(double* v, double px, double py, double pz, double gx, double gy, double gz)
{
    v[0] -= px * gx; v[1] -= px * gy; v[2] -= px * gz;
    v[3] -= py * gx; v[4] -= py * gy; v[5] -= py * gz;
    v[6] -= pz * gx; v[7] -= pz * gy; v[8] -= pz * gz;
}

/// All threads of the block must call it (no early return before); slot = 9 fixed-point accumulators.
__device__ __forceinline__ void blockReduceAddVirial(const double* v, double* slot)
{
    __shared__ long long vir_sums_fx[9][32];
    const int lane = threadIdx.x & 31;
    const int warp_id = threadIdx.x >> 5;
    #pragma unroll
    for (int c = 0; c < 9; ++c) {
        long long x = fix_round(v[c], GFNFF_FIX_ENERGY);
        for (int off = 16; off > 0; off >>= 1) x += __shfl_down_sync(0xffffffffu, x, off);
        if (lane == 0) vir_sums_fx[c][warp_id] = x;
    }
    __syncthreads();
    const int num_warps = (blockDim.x + 31) >> 5;
    if (warp_id == 0) {
        #pragma unroll
        for (int c = 0; c < 9; ++c) {
            long long s = (lane < num_warps) ? vir_sums_fx[c][lane] : 0ll;
            for (int off = 16; off > 0; off >>= 1) s += __shfl_down_sync(0xffffffffu, s, off);
            if (lane == 0 && s != 0ll)
                atomicAdd(reinterpret_cast<unsigned long long*>(&slot[c]), static_cast<unsigned long long>(s));
        }
    }
}

/// Pair term with gradient fac*d on the first atom and -fac*d on the second (d = x_i - x_j as the kernel formed it)
__device__ __forceinline__ void vir_pair(double* v, double fac, double dx, double dy, double dz)
{
    vir_outer(v, dx, dy, dz, fac * dx, fac * dy, fac * dz);
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
    applyMIC(dx, dy, dz);   // Nov 2026 (S1): a wrap by a lattice vector is no motion (CPU maxAtomDisplacement, A3b/A4a)
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
    long long local_Efx = 0;

    // K2 (Oct 2026): grid-stride loop. With the full grid (gridDim*blockDim >= n) every thread
    // handles at most one pair, exactly as before; a capped grid (background launch) loops.
    const int stride = gridDim.x * blockDim.x;
    for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < n; tid += stride) {
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
            local_Efx += fix_round(-C6[tid] * disp_sum * zetac6[tid], GFNFF_FIX_ENERGY_DISP);   // per pair: list order does not matter

            // Gradient: dE/dr via chain rule
            double d6   = -6.0 * r4 * t6 * t6;
            double d8   = -8.0 * r4 * r2 * t8 * t8;
            double dEdr = -C6[tid] * zetac6[tid] * (d6 + 2.0 * r4r2ij[tid] * d8) * rij;
            double fac  = dEdr / rij;

            add_grad(grad, i,  fac*dx,  fac*dy,  fac*dz);
            add_grad(grad, j, -fac*dx, -fac*dy, -fac*dz);

            if (dc6dcn_ij && dEdcn) {
                double disp_value = disp_sum * zetac6[tid];
                fix_add(&dEdcn[i], -dc6dcn_ij[tid] * disp_value, GFNFF_FIX_FORCE);
                fix_add(&dEdcn[j], -dc6dcn_ji[tid] * disp_value, GFNFF_FIX_FORCE);
            }
        }
    }

    blockReduceAddEnergyFix(local_Efx, energy);
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
                const double g = gamma_ij[tid];
                double fac;
                if (g * g * r2 >= GFNFF_COUL_FARFIELD_X2) {
                    // Far field (Claude Generated, Oct 2026): erf == 1, exp term < 2.3e-16, see
                    // GFNFF_COUL_FARFIELD_X2. E = qq/r, fac = dEdr/r = -qq/r^3. sqrt and division
                    // (no rsqrt): this path carries a per-pair cutoff, tested on rij above.
                    local_E = qi * qj / rij;
                    fac     = -(qi * qj) / (r2 * rij);
                } else {
                    double gamma_r = g * rij;
                    double erf_v   = erf(gamma_r);
                    local_E = qi * qj * erf_v / rij;

                    static const double inv_sqrt_pi = 0.5641895835477563;
                    double exp_v  = exp(-gamma_r * gamma_r);
                    double derf   = g * exp_v * (2.0 * inv_sqrt_pi);
                    double dEdr   = qi * qj * (derf / rij - erf_v / (rij * rij));
                    fac           = dEdr / rij;
                }
                add_grad(grad, i,  fac*dx,  fac*dy,  fac*dz);
                add_grad(grad, j, -fac*dx, -fac*dy, -fac*dz);
            }
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// Implicit all-pairs Coulomb, see gfnff_kernels.cuh. Pair arithmetic identical to k_coulomb.
__global__ void k_coulomb_implicit(
    int natoms,
    const double* __restrict__ alp,
    double                     r_cut,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ charges,
    double*                    grad,
    double*                    energy)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;

    if (i < natoms) {
        const double qi = charges[i];
        if (!isnan(qi)) {
            static const double inv_sqrt_pi = 0.5641895835477563;
            const bool no_cut = isinf(r_cut);   // default: no Coulomb cutoff (coulomb_r_cut)
            double gx = 0.0, gy = 0.0, gz = 0.0, e_half = 0.0;
            for (int j = 0; j < natoms; ++j) {
                if (j == i) continue;
                double dx  = cx[i] - cx[j];
                double dy  = cy[i] - cy[j];
                double dz  = cz[i] - cz[j];
                applyMIC(dx, dy, dz);
                double r2  = dx*dx + dy*dy + dz*dz;
                const double asum = alp[i] + alp[j];
                if (no_cut && r2 >= GFNFF_COUL_FARFIELD_X2 * asum) {
                    // Far field, no cutoff (Claude Generated, Oct 2026): same arithmetic as
                    // k_coulomb_tiles (erf == 1, exp term < 2.3e-16, rsqrt).
                    const double qj = charges[j];
                    if (isnan(qj)) continue;
                    const double ri = rsqrt(r2);
                    const double e  = qi * qj * ri;
                    e_half += 0.5 * e;
                    const double f = -e * ri * ri;
                    gx += f*dx; gy += f*dy; gz += f*dz;
                    continue;
                }
                double rij = sqrt(r2);
                if (rij > r_cut || rij < 1e-10) continue;
                double qj = charges[j];
                if (isnan(qj)) continue;
                double fac;
                if (r2 >= GFNFF_COUL_FARFIELD_X2 * asum) {
                    // Far field inside a finite cutoff: sqrt and division (Claude Generated, Oct 2026).
                    e_half += 0.5 * (qi * qj / rij);
                    fac     = -(qi * qj) / (r2 * rij);
                } else {
                    double gamma_ij = 1.0 / sqrt(asum);
                    double gamma_r = gamma_ij * rij;
                    double erf_v   = erf(gamma_r);
                    e_half += 0.5 * (qi * qj * erf_v / rij);
                    double exp_v  = exp(-gamma_r * gamma_r);
                    double derf   = gamma_ij * exp_v * (2.0 * inv_sqrt_pi);
                    double dEdr   = qi * qj * (derf / rij - erf_v / (rij * rij));
                    fac           = dEdr / rij;
                }
                gx += fac*dx; gy += fac*dy; gz += fac*dz;
            }
            local_E = e_half;
            add_grad(grad, i, gx, gy, gz);   // atomic: other kernels may touch atom i concurrently
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Kernel 3b: tiled implicit Coulomb, each pair once (Claude Generated, Sep 2026)
//
// E = sum_{i<j} q_i q_j erf(gamma_ij r_ij) / r_ij,   gamma_ij = 1/sqrt(alp_i + alp_j)
// dE/dr = q_i q_j (2 gamma/sqrt(pi) exp(-gamma^2 r^2) / r - erf(gamma r) / r^2)
// Identical to k_coulomb_implicit, but every pair is evaluated once: that kernel pays the
// erf + exp of each pair twice (once from each atom's row), which dominated the GFN-FF MD step
// on large systems (polymer_2x, 7320 atoms: 85 ms of ~180 ms per energy call, nsys Sep 29).
// Thread t of tile (bi, bj) owns atom i = bi*T + t; it walks the j tile in a rotated order
// (jj = (t + k) mod T) so that concurrent shared-memory atomics on the j gradients rarely collide.
// Far field (Oct 2026): pairs with r^2 >= 36 (alp_i + alp_j) skip erf/exp (GFNFF_COUL_FARFIELD_X2);
// polymer_2x: 11.2 -> 2.4 ms per call on an RTX 5080 (CUPTI, n = 22 calls), energy unchanged at 1e-13 Eh.
// ============================================================================

__global__ void k_coulomb_tiles(
    int natoms,
    const int* __restrict__    tiles,
    const double* __restrict__ alp,
    double                     r_cut,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ charges,
    double*                    grad,
    double*                    energy)
{
    constexpr int T = GFNFF_COUL_TILE;
    __shared__ double sx[T], sy[T], sz[T], sq[T], sa[T];
    __shared__ unsigned long long sgx[T], sgy[T], sgz[T];   // fixed-point partner gradients (order-independent)

    const int bi = tiles[2 * blockIdx.x];
    const int bj = tiles[2 * blockIdx.x + 1];
    const int t  = threadIdx.x;
    const int i  = bi * T + t;
    const int jl = bj * T + t;

    if (jl < natoms) {
        sx[t] = cx[jl]; sy[t] = cy[jl]; sz[t] = cz[jl];
        sq[t] = charges[jl]; sa[t] = alp[jl];
    } else {
        sx[t] = sy[t] = sz[t] = 0.0; sq[t] = nan(""); sa[t] = 1.0;
    }
    sgx[t] = 0ull; sgy[t] = 0ull; sgz[t] = 0ull;
    __syncthreads();

    double local_E = 0.0;
    double gx = 0.0, gy = 0.0, gz = 0.0;
    const bool active = (i < natoms);
    const double xi = active ? cx[i] : 0.0;
    const double yi = active ? cy[i] : 0.0;
    const double zi = active ? cz[i] : 0.0;
    const double qi = active ? charges[i] : nan("");
    const double ai = active ? alp[i] : 1.0;
    const bool diag = (bi == bj);
    const bool no_cut = isinf(r_cut);   // default: no Coulomb cutoff (coulomb_r_cut)
    static const double two_inv_sqrt_pi = 1.1283791670955126;

    if (active && !isnan(qi)) {
        for (int k = 0; k < T; ++k) {
            const int jj = (t + k) & (T - 1);
            const int j  = bj * T + jj;
            if (j >= natoms || (diag && j <= i)) continue;
            const double qj = sq[jj];
            if (isnan(qj)) continue;
            double dx = xi - sx[jj];
            double dy = yi - sy[jj];
            double dz = zi - sz[jj];
            applyMIC(dx, dy, dz);
            const double r2  = dx * dx + dy * dy + dz * dz;
            const double asum = ai + sa[jj];
            const double qq   = qi * qj;
            double fac;
            if (no_cut && r2 >= GFNFF_COUL_FARFIELD_X2 * asum) {
                // Far field, no cutoff (Claude Generated, Oct 2026): erf(gamma r) == 1 and the exp
                // term is < 2.3e-16 (GFNFF_COUL_FARFIELD_X2). rsqrt (1 ulp) replaces sqrt and two
                // divisions. r2 >= 36 asum > 0, so the r < 1e-10 guard below cannot apply here.
                const double ri = rsqrt(r2);
                const double e  = qq * ri;
                local_E += e;
                fac = -e * ri * ri;
                const double fx = fac * dx, fy = fac * dy, fz = fac * dz;
                gx += fx; gy += fy; gz += fz;
                fix_add_u64(&sgx[jj], -fx, GFNFF_FIX_FORCE);
                fix_add_u64(&sgy[jj], -fy, GFNFF_FIX_FORCE);
                fix_add_u64(&sgz[jj], -fz, GFNFF_FIX_FORCE);
                continue;
            }
            const double rij = sqrt(r2);
            if (rij > r_cut || rij < 1e-10) continue;
            if (r2 >= GFNFF_COUL_FARFIELD_X2 * asum) {
                // Far field inside a finite cutoff: same shortcut with sqrt and division,
                // taken after applyMIC and the unchanged cutoff test.
                local_E += qq / rij;
                fac = -qq / (r2 * rij);
            } else {
                const double gamma_ij = 1.0 / sqrt(asum);
                const double gamma_r  = gamma_ij * rij;
                const double erf_v    = erf(gamma_r);
                local_E += qq * erf_v / rij;
                const double derf = gamma_ij * exp(-gamma_r * gamma_r) * two_inv_sqrt_pi;
                fac = qq * (derf / rij - erf_v / r2) / rij;
            }
            const double fx = fac * dx, fy = fac * dy, fz = fac * dz;
            gx += fx; gy += fy; gz += fz;
            fix_add_u64(&sgx[jj], -fx, GFNFF_FIX_FORCE);
            fix_add_u64(&sgy[jj], -fy, GFNFF_FIX_FORCE);
            fix_add_u64(&sgz[jj], -fz, GFNFF_FIX_FORCE);
        }
        add_grad(grad, i, gx, gy, gz);
    }
    __syncthreads();
    if (jl < natoms && (sgx[t] != 0ull || sgy[t] != 0ull || sgz[t] != 0ull)) {
        const double inv = 1.0 / GFNFF_FIX_FORCE;   // exact: the value is a multiple of 2^-52 and is re-added unchanged
        add_grad(grad, jl, static_cast<double>(static_cast<long long>(sgx[t])) * inv,
                 static_cast<double>(static_cast<long long>(sgy[t])) * inv,
                 static_cast<double>(static_cast<long long>(sgz[t])) * inv);
    }

    blockReduceAddEnergy(local_E, energy);
}

__global__ void k_add_coulomb_parts(int n3, int nparts, const double* __restrict__ parts,
                                    const double* __restrict__ part_E,
                                    double* grad, double* energy)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n3) {
        // the helper devices' partial gradients are fixed-point accumulators (see fix_add): integer sum
        unsigned long long s = 0ull;
        for (int p = 0; p < nparts; ++p) s += reinterpret_cast<const unsigned long long*>(parts)[static_cast<size_t>(p) * n3 + k];
        atomicAdd(reinterpret_cast<unsigned long long*>(&grad[k]), s);   // other streams add to d_grad concurrently
    }
    if (k == 0) {
        unsigned long long e = 0ull;   // helper energies are fixed-point accumulators (blockReduceAddEnergy): integer sum
        for (int p = 0; p < nparts; ++p) e += reinterpret_cast<const unsigned long long*>(part_E)[p];
        atomicAdd(reinterpret_cast<unsigned long long*>(energy), e);
    }
}

__global__ void k_add_disp_parts(int N, int nparts, const double* __restrict__ gparts,
                                 const double* __restrict__ dparts, const double* __restrict__ part_E,
                                 double* grad, double* dEdcn, double* energy)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < 3 * N) {
        unsigned long long s = 0ull;   // fixed-point partial gradients: integer sum
        for (int p = 0; p < nparts; ++p) s += reinterpret_cast<const unsigned long long*>(gparts)[static_cast<size_t>(p) * 3 * N + k];
        atomicAdd(reinterpret_cast<unsigned long long*>(&grad[k]), s);   // other streams add to d_grad concurrently
    }
    if (k < N) {
        unsigned long long s = 0ull;
        for (int p = 0; p < nparts; ++p) s += reinterpret_cast<const unsigned long long*>(dparts)[static_cast<size_t>(p) * N + k];
        atomicAdd(reinterpret_cast<unsigned long long*>(&dEdcn[k]), s);
    }
    if (k == 0) {
        unsigned long long e = 0ull;   // helper energies are fixed-point accumulators (blockReduceAddEnergy): integer sum
        for (int p = 0; p < nparts; ++p) e += reinterpret_cast<const unsigned long long*>(part_E)[p];
        atomicAdd(reinterpret_cast<unsigned long long*>(energy), e);
    }
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
                fix_add(&dEdcn[i], yy * ff[tid] * cnfak_i[tid], GFNFF_FIX_FORCE);
                fix_add(&dEdcn[j], yy * ff[tid] * cnfak_j[tid], GFNFF_FIX_FORCE);
            }

            if (zz_hb && nr_hb[tid] >= 1 && hb_H_atom) {
                constexpr double t1 = 0.1;
                double zz = t1 * alpha_orig * dr * dr * E;
                int H = hb_H_atom[tid];
                if (H >= 0) {
                    fix_add(&zz_hb[H], zz, GFNFF_FIX_FORCE);
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

// Nov 2026 (S1): the local-term kernels are templates on PBC. The false instantiation keeps the former statements
// verbatim (a rewrite through local arrays changed the FMA contraction of the bonded kernels by one fixed-point
// quantum on urea197/peo201; the SASS of k_*<false> is compared against the pre-S1 plugin), the true instantiation
// builds the CPU image frame. FFWorkspaceGPU launches k_*<true> only with an active cell.
// A1 (Nov 2026): the body is a device function; the kernel k_angles<PBC> keeps its former signature and statements
// (an added kernel parameter changed the register/FMA choices and the last bits without a cell); k_angles_vir is
// the virial instantiation (cell only).
template <bool PBC, bool VIR>
__device__ __forceinline__ void k_angles_body(
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
    double*                    energy,
    double*                    vir)      // A1 (Nov 2026)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;
    [[maybe_unused]] double vacc[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid], k = idx_k[tid];
        int zi = ati[tid] - 1, zj = atj[tid] - 1, zk = atk[tid] - 1;

        if (zi >= 0 && zi < 86 && zj >= 0 && zj < 86 && zk >= 0 && zk < 86) {
            const double scale = 4.0 / 3.0;
            double rcov_i = d_rcov_d3[zi] * scale;
            double rcov_j = d_rcov_d3[zj] * scale;
            double rcov_k = d_rcov_d3[zk] * scale;

            // Nov 2026 (S1, PBC): local image frame, i and k as the images next to the centre j (slots 0..2 =
            // i, j, k; CPU calcAngles); the angle and its gradient read the frame too. Else: the former statements.
            [[maybe_unused]] double lx[3], ly[3], lz[3];
            double dij_x, dij_y, dij_z, dkj_x, dkj_y, dkj_z;
            if constexpr (PBC) {
                lx[1] = cx[j]; ly[1] = cy[j]; lz[1] = cz[j];
                micPlace(cx, cy, cz, i, lx[1], ly[1], lz[1], lx[0], ly[0], lz[0]);
                micPlace(cx, cy, cz, k, lx[1], ly[1], lz[1], lx[2], ly[2], lz[2]);
                dij_x = lx[0] - lx[1];
                dij_y = ly[0] - ly[1];
                dij_z = lz[0] - lz[1];
                dkj_x = lx[2] - lx[1];
                dkj_y = ly[2] - ly[1];
                dkj_z = lz[2] - lz[1];
            } else {
                dij_x = cx[i] - cx[j];
                dij_y = cy[i] - cy[j];
                dij_z = cz[i] - cz[j];
                applyMIC(dij_x, dij_y, dij_z);
                dkj_x = cx[k] - cx[j];
                dkj_y = cy[k] - cy[j];
                dkj_z = cz[k] - cz[j];
                applyMIC(dkj_x, dkj_y, dkj_z);
            }

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
            double cos_a;
            if constexpr (PBC)
                cos_a = cos_angle_grad(lx, ly, lz, 0, 1, 2,
                                       d_ri_x, d_ri_y, d_ri_z,
                                       d_rk_x, d_rk_y, d_rk_z,
                                       d_rj_x, d_rj_y, d_rj_z, true);
            else
                cos_a = cos_angle_grad(cx, cy, cz, i, j, k,
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
            if constexpr (VIR) {   // A1: frame positions relative to the centre j (slot 1)
                vir_outer(vacc, dij_x, dij_y, dij_z, gi_x + t1_x, gi_y + t1_y, gi_z + t1_z);
                vir_outer(vacc, dkj_x, dkj_y, dkj_z, gk_x + t2_x, gk_y + t2_y, gk_z + t2_z);
            }
        }
    }

    blockReduceAddEnergy(local_E, energy);
    if constexpr (VIR) blockReduceAddVirial(vacc, vir);
}

template <bool PBC>
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
    k_angles_body<PBC, false>(n, idx_i, idx_j, idx_k, ati, atj, atk, fc, theta0, cx, cy, cz, rcov_d3, grad, energy, nullptr);
}

__global__ __launch_bounds__(512) void k_angles_vir(
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
    double*                    energy,
    double*                    vir)      // A1 (Nov 2026)
{
    k_angles_body<true, true>(n, idx_i, idx_j, idx_k, ati, atj, atk, fc, theta0, cx, cy, cz, rcov_d3, grad, energy, vir);
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

// A1 (Nov 2026): the body is a device function; the kernel k_dihedrals<PBC> keeps its former signature and statements
// (an added kernel parameter changed the register/FMA choices and the last bits without a cell); k_dihedrals_vir is
// the virial instantiation (cell only).
template <bool PBC, bool VIR>
__device__ __forceinline__ void k_dihedrals_body(
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
    double*                    energy,
    double*                    vir)      // A1 (Nov 2026)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;
    [[maybe_unused]] double vacc[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    if (tid < n) {
        int ai = idx_i[tid], aj = idx_j[tid], ak = idx_k[tid], al = idx_l[tid];
        int zi = ati[tid]-1, zj = atj[tid]-1, zk = atk[tid]-1, zl = atl[tid]-1;

        if (zi>=0&&zi<86&&zj>=0&&zj<86&&zk>=0&&zk<86&&zl>=0&&zl<86) {
            const double scale = 4.0/3.0;
            double ri = d_rcov_d3[zi]*scale, rj = d_rcov_d3[zj]*scale;
            double rk = d_rcov_d3[zk]*scale, rl = d_rcov_d3[zl]*scale;

            double atcutt = is_nci[tid] ? 0.305 : 0.505;

            // Nov 2026 (S1, PBC): local image frame (slots 0..3 = i, j, k, l; placeChainOfFour of the CPU);
            // every distance and vector below reads the frame. Else: the former statements on cx/cy/cz.
            [[maybe_unused]] double lx[4], ly[4], lz[4];
            if constexpr (PBC)
                micChainOfFour(cx, cy, cz, ai, aj, ak, al, lx, ly, lz);

            auto sq_dist = [&](int a, int b) {
                if constexpr (PBC) {
                    double dx = lx[b] - lx[a];
                    double dy = ly[b] - ly[a];
                    double dz = lz[b] - lz[a];
                    return dx*dx+dy*dy+dz*dz;
                } else {
                    double dx = cx[b] - cx[a];
                    double dy = cy[b] - cy[a];
                    double dz = cz[b] - cz[a];
                    return dx*dx+dy*dy+dz*dz;
                }
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

            double rij2, rjk2, rkl2;
            if constexpr (PBC) {
                rij2 = sq_dist(0, 1);
                rjk2 = sq_dist(1, 2);
                rkl2 = sq_dist(2, 3);
            } else {
                rij2 = sq_dist(ai, aj);
                rjk2 = sq_dist(aj, ak);
                rkl2 = sq_dist(ak, al);
            }

            double d_ij = make_damp(rij2, ri, rj);
            double d_jk = make_damp(rjk2, rj, rk);
            double d_kl = make_damp(rkl2, rk, rl);
            double damp  = d_ij * d_jk * d_kl;

            double d2_ij = make_damp2(rij2, ri, rj);
            double d2_jk = make_damp2(rjk2, rj, rk);
            double d2_kl = make_damp2(rkl2, rk, rl);

            double derivate[4][3];
            double phi;
            if constexpr (PBC)
                phi = dihedral_angle_grad(lx, ly, lz, 0, 1, 2, 3, derivate, true);
            else
                phi = dihedral_angle_grad(cx, cy, cz, ai, aj, ak, al, derivate, true);

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
            if constexpr (VIR) {   // A1: frame slots relative to slot 1 (j)
                for (int a = 0; a < 4; ++a)
                    vir_outer(vacc, lx[a] - lx[1], ly[a] - ly[1], lz[a] - lz[1],
                              dEdphi * damp * derivate[a][0], dEdphi * damp * derivate[a][1], dEdphi * damp * derivate[a][2]);
            }

            // Gradient from damping: E_raw * d(damp)/d(rXY²) * rXY_vec
            {
                double dx, dy, dz;
                if constexpr (PBC) { dx=lx[0]-lx[1]; dy=ly[0]-ly[1]; dz=lz[0]-lz[1]; }
                else               { dx=cx[ai]-cx[aj]; dy=cy[ai]-cy[aj]; dz=cz[ai]-cz[aj]; }
                double fac = E_raw * d2_ij * d_jk * d_kl;
                add_grad(grad, ai,  fac*dx,  fac*dy,  fac*dz);
                add_grad(grad, aj, -fac*dx, -fac*dy, -fac*dz);
                if constexpr (VIR) vir_pair(vacc, fac, dx, dy, dz);
            }
            {
                double dx, dy, dz;
                if constexpr (PBC) { dx=lx[1]-lx[2]; dy=ly[1]-ly[2]; dz=lz[1]-lz[2]; }
                else               { dx=cx[aj]-cx[ak]; dy=cy[aj]-cy[ak]; dz=cz[aj]-cz[ak]; }
                double fac = E_raw * d_ij * d2_jk * d_kl;
                add_grad(grad, aj,  fac*dx,  fac*dy,  fac*dz);
                add_grad(grad, ak, -fac*dx, -fac*dy, -fac*dz);
                if constexpr (VIR) vir_pair(vacc, fac, dx, dy, dz);
            }
            {
                double dx, dy, dz;
                if constexpr (PBC) { dx=lx[2]-lx[3]; dy=ly[2]-ly[3]; dz=lz[2]-lz[3]; }
                else               { dx=cx[ak]-cx[al]; dy=cy[ak]-cy[al]; dz=cz[ak]-cz[al]; }
                double fac = E_raw * d_ij * d_jk * d2_kl;
                add_grad(grad, ak,  fac*dx,  fac*dy,  fac*dz);
                add_grad(grad, al, -fac*dx, -fac*dy, -fac*dz);
                if constexpr (VIR) vir_pair(vacc, fac, dx, dy, dz);
            }
        }
    }

    blockReduceAddEnergy(local_E, energy);
    if constexpr (VIR) blockReduceAddVirial(vacc, vir);
}

template <bool PBC>
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
    k_dihedrals_body<PBC, false>(n, idx_i, idx_j, idx_k, idx_l, ati, atj, atk, atl, V, phi0, n_period, is_nci, cx, cy, cz, rcov_d3_unused, grad, energy, nullptr);
}

__global__ __launch_bounds__(512) void k_dihedrals_vir(
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
    double*                    energy,
    double*                    vir)      // A1 (Nov 2026)
{
    k_dihedrals_body<true, true>(n, idx_i, idx_j, idx_k, idx_l, ati, atj, atk, atl, V, phi0, n_period, is_nci, cx, cy, cz, rcov_d3_unused, grad, energy, vir);
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

// A1 (Nov 2026): the body is a device function; the kernel k_inversions<PBC> keeps its former signature and statements
// (an added kernel parameter changed the register/FMA choices and the last bits without a cell); k_inversions_vir is
// the virial instantiation (cell only).
template <bool PBC, bool VIR>
__device__ __forceinline__ void k_inversions_body(
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
    double*                    energy,
    double*                    vir)      // A1 (Nov 2026)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;
    [[maybe_unused]] double vacc[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    if (tid < n) {
        int i = idx_i[tid], j = idx_j[tid], k = idx_k[tid], l = idx_l[tid];

        const double rcov_scale = 4.0 / 3.0;
        const double atcutt = 0.505;

        int zi = ati[tid]-1, zj = atj[tid]-1, zk = atk[tid]-1, zl = atl[tid]-1;
        double rcov_c = (zi >= 0 && zi < 87) ? d_rcov_d3[zi] * rcov_scale : 1.0 * rcov_scale;
        double rcov_1 = (zj >= 0 && zj < 87) ? d_rcov_d3[zj] * rcov_scale : 1.0 * rcov_scale;
        double rcov_2 = (zk >= 0 && zk < 87) ? d_rcov_d3[zk] * rcov_scale : 1.0 * rcov_scale;
        double rcov_3 = (zl >= 0 && zl < 87) ? d_rcov_d3[zl] * rcov_scale : 1.0 * rcov_scale;

        // Nov 2026 (S1, PBC): local image frame, the centre i at its position and its three neighbours as the
        // images next to it (slots 0..3 = i, j, k, l; CPU calcInversions). Else: the former statements.
        [[maybe_unused]] double lx[4], ly[4], lz[4];
        if constexpr (PBC) {
            lx[0] = cx[i]; ly[0] = cy[i]; lz[0] = cz[i];
            micPlace(cx, cy, cz, j, lx[0], ly[0], lz[0], lx[1], ly[1], lz[1]);
            micPlace(cx, cy, cz, k, lx[0], ly[0], lz[0], lx[2], ly[2], lz[2]);
            micPlace(cx, cy, cz, l, lx[0], ly[0], lz[0], lx[3], ly[3], lz[3]);
        }

        auto sq_dist3 = [&](int a, int b) {
            if constexpr (PBC) {
                double dx = lx[b] - lx[a];
                double dy = ly[b] - ly[a];
                double dz = lz[b] - lz[a];
                return dx*dx+dy*dy+dz*dz;
            } else {
                double dx = cx[b] - cx[a];
                double dy = cy[b] - cy[a];
                double dz = cz[b] - cz[a];
                return dx*dx+dy*dy+dz*dz;
            }
        };

        double rij_sq, rjk_sq, rjl_sq;
        if constexpr (PBC) {
            rij_sq = sq_dist3(1, 0);
            rjk_sq = sq_dist3(1, 2);
            rjl_sq = sq_dist3(1, 3);
        } else {
            rij_sq = sq_dist3(j, i);
            rjk_sq = sq_dist3(j, k);
            rjl_sq = sq_dist3(j, l);
        }

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
        double omega;
        if constexpr (PBC)
            omega = inversion_angle_grad(lx, ly, lz, 0, 1, 2, 3, derivate, true);
        else
            omega = inversion_angle_grad(cx, cy, cz, i, j, k, l, derivate, true);
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
        if constexpr (VIR) {   // A1: frame slots relative to the centre (slot 0)
            for (int a = 1; a < 4; ++a)
                vir_outer(vacc, lx[a] - lx[0], ly[a] - ly[0], lz[a] - lz[0],
                          dEdOmega * derivate[a][0], dEdOmega * derivate[a][1], dEdOmega * derivate[a][2]);
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
            double dx, dy, dz;
            if constexpr (PBC) { dx=lx[0]-lx[1]; dy=ly[0]-ly[1]; dz=lz[0]-lz[1]; }
            else               { dx=cx[i]-cx[j]; dy=cy[i]-cy[j]; dz=cz[i]-cz[j]; }
            double fac = E_raw * dd_ij * damp_jk * damp_jl;
            add_grad(grad, i,  fac*dx,  fac*dy,  fac*dz);
            add_grad(grad, j, -fac*dx, -fac*dy, -fac*dz);
            if constexpr (VIR) vir_pair(vacc, fac, dx, dy, dz);
        }
        {
            double dx, dy, dz;
            if constexpr (PBC) { dx=lx[2]-lx[1]; dy=ly[2]-ly[1]; dz=lz[2]-lz[1]; }
            else               { dx=cx[k]-cx[j]; dy=cy[k]-cy[j]; dz=cz[k]-cz[j]; }
            double fac = E_raw * damp_ij * dd_jk * damp_jl;
            add_grad(grad, k,  fac*dx,  fac*dy,  fac*dz);
            add_grad(grad, j, -fac*dx, -fac*dy, -fac*dz);
            if constexpr (VIR) vir_pair(vacc, fac, dx, dy, dz);
        }
        {
            double dx, dy, dz;
            if constexpr (PBC) { dx=lx[3]-lx[1]; dy=ly[3]-ly[1]; dz=lz[3]-lz[1]; }
            else               { dx=cx[l]-cx[j]; dy=cy[l]-cy[j]; dz=cz[l]-cz[j]; }
            double fac = E_raw * damp_ij * damp_jk * dd_jl;
            add_grad(grad, l,  fac*dx,  fac*dy,  fac*dz);
            add_grad(grad, j, -fac*dx, -fac*dy, -fac*dz);
            if constexpr (VIR) vir_pair(vacc, fac, dx, dy, dz);
        }
    }

    blockReduceAddEnergy(local_E, energy);
    if constexpr (VIR) blockReduceAddVirial(vacc, vir);
}

template <bool PBC>
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
    k_inversions_body<PBC, false>(n, idx_i, idx_j, idx_k, idx_l, ati, atj, atk, atl, fc, omega0, C0, C1, C2, potential_type, cx, cy, cz, grad, energy, nullptr);
}

__global__ __launch_bounds__(512) void k_inversions_vir(
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
    double*                    energy,
    double*                    vir)      // A1 (Nov 2026)
{
    k_inversions_body<true, true>(n, idx_i, idx_j, idx_k, idx_l, ati, atj, atk, atl, fc, omega0, C0, C1, C2, potential_type, cx, cy, cz, grad, energy, vir);
}

// Claude Generated (Nov 2026, S1): the cell symbols live in this translation unit; a cudaMemcpyToSymbol on an
// `extern __constant__` declaration in another unit failed (the error then surfaced in the next checked call of
// the workspace set-up as cudaErrorInvalidDevice), so the upload is done here and its status returned.
cudaError_t gfnffUploadCell(const double* lattice_9, const double* lattice_inv_9, int has_pbc)
{
    if (has_pbc) {
        cudaError_t e = cudaMemcpyToSymbol(d_lattice, lattice_9, 9 * sizeof(double));
        if (e != cudaSuccess) return e;
        e = cudaMemcpyToSymbol(d_lattice_inv, lattice_inv_9, 9 * sizeof(double));
        if (e != cudaSuccess) return e;
    }
    return cudaMemcpyToSymbol(d_has_pbc, &has_pbc, sizeof(int));
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
// Upload D4 reference CN and nref tables to constant memory (once at init)
// Claude Generated (March 2026): Phase 8 — constant memory for refcn/refn
// ============================================================================
void upload_refcn_const(const double* data, int n)
{
    int count = (n < D4_MAX_ELEM * D4_MAX_REF) ? n : D4_MAX_ELEM * D4_MAX_REF;
    cudaMemcpyToSymbol(d_refcn_const, data, count * sizeof(double));
}

void upload_refn_const(const int* data, int n)
{
    int count = (n < D4_MAX_ELEM) ? n : D4_MAX_ELEM;
    cudaMemcpyToSymbol(d_refn_const, data, count * sizeof(int));
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

// A1 (Nov 2026): the body is a device function; the kernel k_storsions<PBC> keeps its former signature and statements
// (an added kernel parameter changed the register/FMA choices and the last bits without a cell); k_storsions_vir is
// the virial instantiation (cell only).
template <bool PBC, bool VIR>
__device__ __forceinline__ void k_storsions_body(
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
    double*                    energy,
    double*                    vir)      // A1 (Nov 2026)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;
    [[maybe_unused]] double vacc[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    [[maybe_unused]] double vpx[4], vpy[4], vpz[4];   // A1: frame positions relative to slot 1 (VIR only)

    if (tid < n) {
        int ai = idx_i[tid], aj = idx_j[tid], ak = idx_k[tid], al = idx_l[tid];

        // Nov 2026 (S1, PBC): local image frame as in k_dihedrals (CPU calcSTorsions, placeChainOfFour)
        double derivate[4][3];
        double phi;
        if constexpr (PBC) {
            double lx[4], ly[4], lz[4];
            micChainOfFour(cx, cy, cz, ai, aj, ak, al, lx, ly, lz);
            phi = dihedral_angle_grad(lx, ly, lz, 0, 1, 2, 3, derivate, true);
            if constexpr (VIR)
                for (int a = 0; a < 4; ++a) { vpx[a] = lx[a] - lx[1]; vpy[a] = ly[a] - ly[1]; vpz[a] = lz[a] - lz[1]; }
        } else {
            phi = dihedral_angle_grad(cx, cy, cz, ai, aj, ak, al, derivate, true);
        }

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
        if constexpr (VIR)
            for (int a = 0; a < 4; ++a)
                vir_outer(vacc, vpx[a], vpy[a], vpz[a], dEdphi * derivate[a][0], dEdphi * derivate[a][1], dEdphi * derivate[a][2]);
    }

    blockReduceAddEnergy(local_E, energy);
    if constexpr (VIR) blockReduceAddVirial(vacc, vir);
}

template <bool PBC>
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
    k_storsions_body<PBC, false>(n, idx_i, idx_j, idx_k, idx_l, erefhalf, cx, cy, cz, grad, energy, nullptr);
}

__global__ __launch_bounds__(512) void k_storsions_vir(
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
    double*                    energy,
    double*                    vir)      // A1 (Nov 2026)
{
    k_storsions_body<true, true>(n, idx_i, idx_j, idx_k, idx_l, erefhalf, cx, cy, cz, grad, energy, vir);
}

// ============================================================================
// Kernel 9: Bonded ATM (BATM) for 1,4-pairs
// E = c9 * angr9  where c9 = fi*fj*fk * zb3atm_i*zb3atm_j*zb3atm_k
// fi = clamp(1 - 3*q_i, -4, 4)   (topology charges)
// Reference: ff_workspace_gfnff.cpp:calcBATM (lines 1652-1726)
// ============================================================================

// A1 (Nov 2026): body as a device function, k_batm keeps its former signature (bounds from the .cuh declaration)
template <bool VIR>
__device__ __forceinline__ void k_batm_body(
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
    double*                    energy,
    double*                    vir)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;
    [[maybe_unused]] double vacc[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

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
            if constexpr (VIR) {   // positions relative to i
                vir_outer(vacc, rij_x, rij_y, rij_z, -fij*rij_x + fjk*rjk_x, -fij*rij_y + fjk*rjk_y, -fij*rij_z + fjk*rjk_z);
                vir_outer(vacc, rik_x, rik_y, rik_z, -fik*rik_x - fjk*rjk_x, -fik*rik_y - fjk*rjk_y, -fik*rik_z - fjk*rjk_z);
            }
        }
    }

    blockReduceAddEnergy(local_E, energy);
    if constexpr (VIR) blockReduceAddVirial(vacc, vir);
}

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
    k_batm_body<false>(n, idx_i, idx_j, idx_k, zb3atm_i, zb3atm_j, zb3atm_k, cx, cy, cz, topo_charges, grad, energy, nullptr);
}

__global__ GFNFF_KERNEL_BOUNDS void k_batm_vir(
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
    double*                    energy,
    double*                    vir)
{
    k_batm_body<true>(n, idx_i, idx_j, idx_k, zb3atm_i, zb3atm_j, zb3atm_k, cx, cy, cz, topo_charges, grad, energy, vir);
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

// A1 (Nov 2026): the body is a device function; the kernel k_xbonds<PBC> keeps its former signature and statements
// (an added kernel parameter changed the register/FMA choices and the last bits without a cell); k_xbonds_vir is
// the virial instantiation (cell only).
template <bool PBC, bool VIR>
__device__ __forceinline__ void k_xbonds_body(
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
    double*                    energy,
    double*                    vir)      // A1 (Nov 2026)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;
    [[maybe_unused]] double vacc[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    if (tid < n) {
        const double XB_SCUT_D = 5.0;
        const double XB_BACUT_D = 70.0;
        const double HB_ALP_D = 6.0;
        const double HB_LONGCUT_XB_D = 70.0;
        const double XB_ST_D = 15.0;
        const double XB_SF_D = 0.03;

        int A = idx_i[tid], X = idx_j[tid], B = idx_k[tid];

        // Nov 2026 (S1, PBC): local image frame anchored on X, A and B as the images next to X (CPU
        // calcHalogenBonds, A4b). Else: the former statements (three minimum-image vectors).
        double ax, ay, az, xx, xy, xz, bx, by, bz;
        double rAX_x, rAX_y, rAX_z, rXB_x, rXB_y, rXB_z, rAB_x, rAB_y, rAB_z;
        if constexpr (PBC) {
            xx = cx[X]; xy = cy[X]; xz = cz[X];
            micPlace(cx, cy, cz, A, xx, xy, xz, ax, ay, az);
            micPlace(cx, cy, cz, B, xx, xy, xz, bx, by, bz);
            rAX_x = xx-ax; rAX_y = xy-ay; rAX_z = xz-az;
            rXB_x = bx-xx; rXB_y = by-xy; rXB_z = bz-xz;
            rAB_x = bx-ax; rAB_y = by-ay; rAB_z = bz-az;
        } else {
            ax = cx[A]; ay = cy[A]; az = cz[A];
            xx = cx[X]; xy = cy[X]; xz = cz[X];
            bx = cx[B]; by = cy[B]; bz = cz[B];
            rAX_x = xx-ax; rAX_y = xy-ay; rAX_z = xz-az; applyMIC(rAX_x, rAX_y, rAX_z);
            rXB_x = bx-xx; rXB_y = by-xy; rXB_z = bz-xz; applyMIC(rXB_x, rXB_y, rXB_z);
            rAB_x = bx-ax; rAB_y = by-ay; rAB_z = bz-az; applyMIC(rAB_x, rAB_y, rAB_z);
        }

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
                // Known Issue #46: d(damp_s)/dr = +2*alp*damp_s*t/(r*(1+t)), d(damp_l)/dr = -2*alp*damp_l*u/(r*(1+u))
                // (Fortran rbxgfnff_eg); the sign of the first and 1/r of the second were wrong
                double ddamp_short_dr = 2.0 * HB_ALP_D * damp_short * damp_short_term
                                      / (r_XB * (1.0 + damp_short_term));

                double ratio_long = (r_XB * r_XB) / HB_LONGCUT_XB_D;
                double damp_long_term = 1.0;
                for (int p = 0; p < (int)HB_ALP_D; ++p) damp_long_term *= ratio_long;
                double ddamp_long_dr = -2.0 * HB_ALP_D * damp_long * damp_long_term
                                     / (r_XB * (1.0 + damp_long_term));

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
                if constexpr (VIR) {   // A1: frame anchored on X
                    vir_outer(vacc, ax - xx, ay - xy, az - xz,
                              -dE_drAX * rAX_x * inv_rAX - dE_drAB * rAB_x * inv_rAB,
                              -dE_drAX * rAX_y * inv_rAX - dE_drAB * rAB_y * inv_rAB,
                              -dE_drAX * rAX_z * inv_rAX - dE_drAB * rAB_z * inv_rAB);
                    vir_outer(vacc, bx - xx, by - xy, bz - xz,
                              dE_drXB * rXB_x * inv_rXB + dE_drAB * rAB_x * inv_rAB,
                              dE_drXB * rXB_y * inv_rXB + dE_drAB * rAB_y * inv_rAB,
                              dE_drXB * rXB_z * inv_rXB + dE_drAB * rAB_z * inv_rAB);
                }
            }
        }
    }

    blockReduceAddEnergy(local_E, energy);
    if constexpr (VIR) blockReduceAddVirial(vacc, vir);
}

template <bool PBC>
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
    k_xbonds_body<PBC, false>(n, idx_i, idx_j, idx_k, elem_A, elem_B, q_X, q_B, acidity_X, r_cut, cx, cy, cz, grad, energy, nullptr);
}

__global__ __launch_bounds__(512) void k_xbonds_vir(
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
    double*                    energy,
    double*                    vir)      // A1 (Nov 2026)
{
    k_xbonds_body<true, true>(n, idx_i, idx_j, idx_k, elem_A, elem_B, q_X, q_B, acidity_X, r_cut, cx, cy, cz, grad, energy, vir);
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
                // Known Issue #46 (see the FP64 kernel above)
                float ddamp_short_dr = 2.0f * HB_ALP_D * damp_short * damp_short_term
                                      / (r_XB * (1.0f + damp_short_term));

                float damp_long_term = 1.0f;
                float ratio_long_g = (r_XB * r_XB) / HB_LONGCUT_XB_D;
                for (int p = 0; p < 6; ++p) damp_long_term *= ratio_long_g;
                float ddamp_long_dr = -2.0f * HB_ALP_D * damp_long * damp_long_term
                                     / (r_XB * (1.0f + damp_long_term));

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

// A1 (Nov 2026): the body is a device function; the kernel k_hbonds<PBC> keeps its former signature and statements
// (an added kernel parameter changed the register/FMA choices and the last bits without a cell); k_hbonds_vir is
// the virial instantiation (cell only).
template <bool PBC, bool VIR>
__device__ __forceinline__ void k_hbonds_body(
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
    double*                    energy,
    double*       __restrict__ e_triple,
    double*                    vir)      // A1 (Nov 2026)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;
    // A1: frame positions relative to H (the frame anchor), contracted with every gradient row written (VIR only)
    [[maybe_unused]] double vacc[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

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

    // Nov 2026 (S1, PBC): local image frame anchored on H (CPU calcHydrogenBonds, A4b): A and B as the images
    // next to H, below the neighbours of B and the acceptor parent C next to B, the neighbours D of C next to C.
    // Else: the former statements (three minimum-image vectors, neighbours read as they are).
    double ax, ay, az, hx, hy, hz, bx, by, bz;
    double rAH_x, rAH_y, rAH_z, rHB_x, rHB_y, rHB_z, rAB_x, rAB_y, rAB_z;
    if constexpr (PBC) {
        hx = cx[H]; hy = cy[H]; hz = cz[H];
        micPlace(cx, cy, cz, A, hx, hy, hz, ax, ay, az);
        micPlace(cx, cy, cz, B, hx, hy, hz, bx, by, bz);
        rAH_x = hx-ax; rAH_y = hy-ay; rAH_z = hz-az;
        rHB_x = bx-hx; rHB_y = by-hy; rHB_z = bz-hz;
        rAB_x = bx-ax; rAB_y = by-ay; rAB_z = bz-az;
    } else {
        ax = cx[A]; ay = cy[A]; az = cz[A];
        hx = cx[H]; hy = cy[H]; hz = cz[H];
        bx = cx[B]; by = cy[B]; bz = cz[B];
        rAH_x = hx-ax; rAH_y = hy-ay; rAH_z = hz-az; applyMIC(rAH_x, rAH_y, rAH_z);
        rHB_x = bx-hx; rHB_y = by-hy; rHB_z = bz-hz; applyMIC(rHB_x, rHB_y, rHB_z);
        rAB_x = bx-ax; rAB_y = by-ay; rAB_z = bz-az; applyMIC(rAB_x, rAB_y, rAB_z);
    }

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
            double nbx, nby, nbz;
            if constexpr (PBC) micPlace(cx, cy, cz, nb, bx, by, bz, nbx, nby, nbz);   // neighbour next to B
            else { nbx = cx[nb]; nby = cy[nb]; nbz = cz[nb]; }
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
                if constexpr (PBC) {
                    double nbx, nby, nbz;
                    micPlace(cx, cy, cz, nb, bx, by, bz, nbx, nby, nbz);   // neighbour next to B
                    lp_vx_c4 += nbx - bx;
                    lp_vy_c4 += nby - by;
                    lp_vz_c4 += nbz - bz;
                } else {
                    lp_vx_c4 += cx[nb] - bx;
                    lp_vy_c4 += cy[nb] - by;
                    lp_vz_c4 += cz[nb] - bz;
                }
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

        // PBC: C next to B (frame); the D-B-C-H chain below reads slots 0..3 = D, B, C, H with D next to C
        double C_x, C_y, C_z;
        if constexpr (PBC) micPlace(cx, cy, cz, C_idx, bx, by, bz, C_x, C_y, C_z);
        else { C_x = cx[C_idx]; C_y = cy[C_idx]; C_z = cz[C_idx]; }

        // eangl: angle bending H...B=C
        {
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
            if (D_idx == B) continue;

            double derivate[4][3];
            double phi;
            if constexpr (PBC) {
                double tx[4] = { 0.0, bx, C_x, hx }, ty[4] = { 0.0, by, C_y, hy }, tz[4] = { 0.0, bz, C_z, hz };
                micPlace(cx, cy, cz, D_idx, C_x, C_y, C_z, tx[0], ty[0], tz[0]);   // D next to C
                phi = dihedral_angle_grad(tx, ty, tz, 0, 1, 2, 3, derivate, false);
            } else {
                phi = dihedral_angle_grad(cx, cy, cz, D_idx, B, C_idx, H, derivate, false);
            }

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
                double nbx, nby, nbz;
                if constexpr (PBC) micPlace(cx, cy, cz, nb, bx, by, bz, nbx, nby, nbz);   // neighbour next to B (as above)
                else { nbx = cx[nb]; nby = cy[nb]; nbz = cz[nb]; }
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
                if constexpr (VIR)
                    vir_outer(vacc, nbx - hx, nby - hy, nbz - hz, -dga_nb_x - dgb_nb_x, -dga_nb_y - dgb_nb_y, -dga_nb_z - dgb_nb_z);
            }
        }

        // Case 3: eangl + etors gradient (simplified — gradient contributions for angle/torsion
        // are complex and involve 4+ atoms per torsion; computed here for the 3 primary atoms only)
        // Full case 3 gradient for angle bending term
        if (ct == 3 && acceptor_parent[tid] >= 0) {
            int C_idx = acceptor_parent[tid];
            double bterm_c3 = -rdamp * qhoutl * etors * acidity_A[tid] * basicity_B[tid] * Q_A * Q_B * global_scale;
            // Angle gradient: H...B=C angle (PBC: C next to B, the frame of the energy part)
            double C_x, C_y, C_z;
            if constexpr (PBC) micPlace(cx, cy, cz, C_idx, bx, by, bz, C_x, C_y, C_z);
            else { C_x = cx[C_idx]; C_y = cy[C_idx]; C_z = cz[C_idx]; }
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
            if constexpr (VIR)
                vir_outer(vacc, C_x - hx, C_y - hy, C_z - hz, bterm_c3 * (-deda_x), bterm_c3 * (-deda_y), bterm_c3 * (-deda_z));
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
                double phi;
                [[maybe_unused]] double Dx = 0.0, Dy = 0.0, Dz = 0.0;   // A1: D in the frame (VIR only)
                if constexpr (PBC) {
                    double tx[4] = { 0.0, bx, C_x, hx }, ty[4] = { 0.0, by, C_y, hy }, tz[4] = { 0.0, bz, C_z, hz };
                    micPlace(cx, cy, cz, D_idx, C_x, C_y, C_z, tx[0], ty[0], tz[0]);   // D next to C
                    phi = dihedral_angle_grad(tx, ty, tz, 0, 1, 2, 3, derivate, true);
                    if constexpr (VIR) { Dx = tx[0]; Dy = ty[0]; Dz = tz[0]; }
                } else {
                    phi = dihedral_angle_grad(cx, cy, cz, D_idx, B, C_idx, H, derivate, true);
                }

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
                if constexpr (VIR) {
                    vir_outer(vacc, Dx - hx, Dy - hy, Dz - hz, t_k * dij * derivate[0][0], t_k * dij * derivate[0][1], t_k * dij * derivate[0][2]);
                    vir_outer(vacc, C_x - hx, C_y - hy, C_z - hz, t_k * dij * derivate[2][0], t_k * dij * derivate[2][1], t_k * dij * derivate[2][2]);
                }
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
                    if constexpr (VIR) {
                        double nbx, nby, nbz;
                        micPlace(cx, cy, cz, nb, bx, by, bz, nbx, nby, nbz);   // neighbour next to B, as the LP vector
                        vir_outer(vacc, nbx - hx, nby - hy, nbz - hz, -gnb_share_x, -gnb_share_y, -gnb_share_z);
                    }
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
    if constexpr (VIR) {   // H is the origin of the frame: its row drops out
        vir_outer(vacc, ax - hx, ay - hy, az - hz, ga_x, ga_y, ga_z);
        vir_outer(vacc, bx - hx, by - hy, bz - hz, gb_x, gb_y, gb_z);
    }

    } // expo <= 15.0
    } // r_AB in range

    // Claude Generated (Oct 2026): per-triple energy on request (FFWorkspaceGPU::setRecordHBondEnergies).
    // local_E is E_HB of this triple, 0 when it is outside r_cut or the out-of-line exponent limit -
    // the same value the CPU workspace records. Written after the term, so the sum is untouched.
    if (e_triple)
        e_triple[tid] = local_E;
    } // tid < n

    blockReduceAddEnergy(local_E, energy);
    if constexpr (VIR) blockReduceAddVirial(vacc, vir);
}

template <bool PBC>
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
    double*                    energy,
    double*       __restrict__ e_triple)
{
    k_hbonds_body<PBC, false>(n, idx_i, idx_j, idx_k, elem_A, elem_B, q_H, q_A, q_B, basicity_A, basicity_B, acidity_A, acidity_B, r_cut, case_type, nb_B_offset, nb_B_count, nb_B_flat, acceptor_parent, nb_C_offset, nb_C_count, nb_C_flat, repz_B, cx, cy, cz, grad, energy, e_triple, nullptr);
}

__global__ __launch_bounds__(512) void k_hbonds_vir(
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
    double*                    energy,
    double*       __restrict__ e_triple,
    double*                    vir)      // A1 (Nov 2026)
{
    k_hbonds_body<true, true>(n, idx_i, idx_j, idx_k, elem_A, elem_B, q_H, q_A, q_B, basicity_A, basicity_B, acidity_A, acidity_B, r_cut, case_type, nb_B_offset, nb_B_count, nb_B_flat, acceptor_parent, nb_C_offset, nb_C_count, nb_C_flat, repz_B, cx, cy, cz, grad, energy, e_triple, vir);
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

    fix_add(energy, E_en + E_self, GFNFF_FIX_ENERGY);
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
    bool                       do_subtract,
    bool                       cn_zero_guard)
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
            // Claude Generated (Oct 2026): gfnff.cn_zero_guard, same exact test as FFWorkspace::postProcess
            // (CPU): at CN_i == 0.0 the energy -q*cnf*sqrt(max(CN,0)) is flat, so TERM 1b is 0 (no threshold).
            if (!(cn_zero_guard && cn_i == 0.0)) {
                double qtmp = q * cnf[i] / (2.0 * sqrt(cn_i) + 1e-16);
                dEdcn[i] -= qtmp;
            }
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
    double        kn,
    double        r2_max)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pairs) return;

    int i = idx_i[tid], j = idx_j[tid];
    double dx = cx[i] - cx[j];
    double dy = cy[i] - cy[j];
    double dz = cz[i] - cz[j];
    applyMIC(dx, dy, dz);
    double r2 = dx*dx + dy*dy + dz*dz;
    // A16 (Nov 2026, Claude Generated): the CPU CN range (cn_cutoff_bohr / cnthr), see k_cn_compute_pairs.
    if (r2 > r2_max) return;
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

    // Phase 2: Normalize and compute derivatives via quotient rule.
    // Claude Generated (Oct 2026): same rule as D4ParameterGenerator::precomputeGaussianWeights()
    // and computeGaussianWeightDerivatives(): always divide by the finite sum (1/sum is finite for
    // any nonzero sum, also for atoms whose reference CNs are far from the actual CN, e.g. Mn and Ru
    // in carbonyl complexes) and guard each value against NaN/Inf from a true underflow, where the
    // reference with the maximal CN gets weight 1 and the derivative 0. The former `norm > 1e-10`
    // switch collapsed such atoms onto the first reference (ROCm MOR41 gfnff: PR05, ED05, PR26-28,
    // ED26 differed from the CPU by up to 1.6 kcal/mol in the dispersion term).
    {
        double max_cn_ref = d_refcn_const[refcn_base];
        for (int r = 1; r < nref; ++r)
            max_cn_ref = fmax(max_cn_ref, d_refcn_const[refcn_base + r]);
        const double inv_norm = 1.0 / norm;
        const double inv_norm2 = 1.0 / (norm * norm);
        for (int r = 0; r < nref; ++r) {
            double g = expw[r] * inv_norm;
            if (!isfinite(g))
                g = (d_refcn_const[refcn_base + r] == max_cn_ref) ? 1.0 : 0.0;
            double dg = (dexpw[r] * norm - expw[r] * dnorm) * inv_norm2;
            if (!isfinite(dg))
                dg = 0.0;
            gw[base + r]  = g;
            dgw[base + r] = dg;
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
    double*       __restrict__ dc6dcn_ji,
    double*       __restrict__ c6_out)
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
    double c6     = 0.0;   // Claude Generated (Sep 2026): C6 itself, same double sum

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

            // C6(i,j) = Σ gw(i,ri) * gw(j,rj) * C6ref (Caldeweyher 2019 CN-weighted C6)
            c6 += gw_i_ri * gw[gw_j_base + rj] * c6ref;
        }
    }

    dc6dcn_ij[tid] = dc6_ij;
    dc6dcn_ji[tid] = dc6_ji;
    if (c6_out) c6_out[tid] = c6;
}

// ---------------------------------------------------------------------------
// Tabulated dc6dcn (Claude Generated, Oct 2026)
// k_dc6dcn_per_pair contracts the reference table for every pair: 3 sums over nref(i) x nref(j) terms.
// Both factors of that contraction depend on one atom and the partner ELEMENT only, so they are
// tabulated once per atom and element type:
//   A[a][t][r'] = sum_r gw_a[r] * C6ref(Z_a, Z_t; r, r')      (partner reference r' of element t)
//   W[a][t][r'] = sum_r gw_a[r] * C6ref(Z_t, Z_a; r', r)
// and the pair needs nref terms only:
//   C6(i,j)      = sum_rj gw_j[rj]  A[i][t_j][rj]
//   dC6/dCN(j)   = sum_rj dgw_j[rj] A[i][t_j][rj]
//   dC6/dCN(i)   = sum_ri dgw_i[ri] W[j][t_i][ri]
// identical to the double loop up to the order of summation (|C6ref| < 1e-20 skipped as there).
// ---------------------------------------------------------------------------
__global__ void k_d4_partner_tables(int N, int nt, const int* __restrict__ type_z, const int* __restrict__ atom_types,
                                    const double* __restrict__ gw, const double* __restrict__ c6_flat,
                                    double* __restrict__ A, double* __restrict__ W)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N * nt) return;
    const int a = tid / nt, t = tid % nt;
    const int za = atom_types[a] - 1, zt = type_z[t] - 1;
    const int nra = d_refn_const[za], nrt = d_refn_const[zt];
    const size_t base_at = (static_cast<size_t>(za) * D4_MAX_ELEM + zt) * D4_MAX_REF * D4_MAX_REF;
    const size_t base_ta = (static_cast<size_t>(zt) * D4_MAX_ELEM + za) * D4_MAX_REF * D4_MAX_REF;
    for (int rp = 0; rp < D4_MAX_REF; ++rp) {
        double sa = 0.0, sw = 0.0;
        if (rp < nrt) {
            for (int r = 0; r < nra; ++r) {
                const double g = gw[a * D4_MAX_REF + r];
                const double c_at = c6_flat[base_at + r * D4_MAX_REF + rp];
                if (fabs(c_at) >= 1e-20) sa += g * c_at;
                const double c_ta = c6_flat[base_ta + rp * D4_MAX_REF + r];
                if (fabs(c_ta) >= 1e-20) sw += g * c_ta;
            }
        }
        A[(static_cast<size_t>(a) * nt + t) * D4_MAX_REF + rp] = sa;
        W[(static_cast<size_t>(a) * nt + t) * D4_MAX_REF + rp] = sw;
    }
}

__global__ void k_dc6dcn_per_pair_tab(int n_pairs, const int* __restrict__ idx_i, const int* __restrict__ idx_j,
                                      const int* __restrict__ atom_types, const int* __restrict__ atype_idx, int nt,
                                      const double* __restrict__ gw, const double* __restrict__ dgw,
                                      const double* __restrict__ A, const double* __restrict__ W,
                                      double* __restrict__ dc6dcn_ij, double* __restrict__ dc6dcn_ji,
                                      double* __restrict__ c6_out)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pairs) return;
    const int ai = idx_i[tid], aj = idx_j[tid];
    const int nri = d_refn_const[atom_types[ai] - 1], nrj = d_refn_const[atom_types[aj] - 1];
    const double* Aij = A + (static_cast<size_t>(ai) * nt + atype_idx[aj]) * D4_MAX_REF;
    const double* Wji = W + (static_cast<size_t>(aj) * nt + atype_idx[ai]) * D4_MAX_REF;
    double c6 = 0.0, dji = 0.0, dij = 0.0;
    for (int rj = 0; rj < nrj; ++rj) {
        c6  += gw[aj * D4_MAX_REF + rj] * Aij[rj];
        dji += dgw[aj * D4_MAX_REF + rj] * Aij[rj];
    }
    for (int ri = 0; ri < nri; ++ri)
        dij += dgw[ai * D4_MAX_REF + ri] * Wji[ri];
    dc6dcn_ij[tid] = dij;
    dc6dcn_ji[tid] = dji;
    if (c6_out) c6_out[tid] = c6;
}

void launch_dc6dcn_tabulated(cudaStream_t s, int N, int nt, int n_pairs, const int* idx_i, const int* idx_j,
                             const int* atom_types, const int* atype_idx, const int* type_z,
                             const double* gw, const double* dgw, const double* c6_flat, double* A, double* W,
                             double* dc6ij, double* dc6ji, double* c6_out, int block)
{
    k_d4_partner_tables<<<(N * nt + 127) / 128, 128, 0, s>>>(N, nt, type_z, atom_types, gw, c6_flat, A, W);
    if (n_pairs <= 0) return;
    k_dc6dcn_per_pair_tab<<<(n_pairs + block - 1) / block, block, 0, s>>>(
        n_pairs, idx_i, idx_j, atom_types, atype_idx, nt, gw, dgw, A, W, dc6ij, dc6ji, c6_out);
}

// Claude Generated (Oct 2026): cost of one dispersion pair for the multi-GPU split. k_dc6dcn_per_pair
// loops over nref(i) x nref(j) reference pairs (1..49), k_dispersion costs about the same as an average
// dc6dcn pair, so cost = nref_i * nref_j + 10. Inclusive prefix sum over the pair list.
__global__ void k_disp_pair_cost(int n, const int* __restrict__ idx_i, const int* __restrict__ idx_j,
                                 const int* __restrict__ atom_types, long long* __restrict__ out)
{
    const int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    out[t] = static_cast<long long>(d_refn_const[atom_types[idx_i[t]] - 1])
           * d_refn_const[atom_types[idx_j[t]] - 1] + 10;
}

void disp_pair_cost_scan(int n, const int* idx_i, const int* idx_j, const int* atom_types,
                         long long* out, cudaStream_t stream)
{
    if (n <= 0) return;
    k_disp_pair_cost<<<(n + 255) / 256, 256, 0, stream>>>(n, idx_i, idx_j, atom_types, out);
    thrust::inclusive_scan(thrust::cuda::par.on(stream), out, out + n, out);
}

void exclusive_scan_int(const int* in, int* out, int n, cudaStream_t stream)
{
    if (n <= 0) return;
    thrust::exclusive_scan(thrust::cuda::par.on(stream), in, in + n, out);
}

// ============================================================================
// Utility: zero device array
// ============================================================================

// ============================================================================
// GPU CN pair list generation (Apr 2026)
// Two-pass kernels replacing the former host O(N^2) CN pair loop.
// ============================================================================

/// Pass 1: count valid pairs. 1 thread = 1 upper-triangle pair (i<j).
__global__ GFNFF_KERNEL_BOUNDS void k_generate_cn_pairs_count(
    int N,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const int*    __restrict__ atom_types,
    double        cutoff_factor,
    double        extra_bohr,
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

    if (rij < cutoff_factor * rcov_sum + extra_bohr && rij > 1e-10) {
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
    double        extra_bohr,
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

    if (rij < cutoff_factor * rcov_sum_val + extra_bohr && rij > 1e-10) {
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
    double        kn,
    double        r2_max)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pairs) return;

    int i = idx_i[tid];
    int j = idx_j[tid];

    double dx = cx[i] - cx[j];
    double dy = cy[i] - cy[j];
    double dz = cz[i] - cz[j];
    applyMIC(dx, dy, dz);
    // A16 (Nov 2026, Claude Generated): the CN range of the CPU (CNCalculator::calculateGFNFFCN*):
    // r2_max = nextafter(cn_cutoff_bohr^2, 0) (strict r^2 < rc^2) or cnthr(cn_accuracy) with
    // cn_cutoff_bohr = 0. The pair list itself reaches f*(rcov_i+rcov_j) + 2 skin, which exceeds
    // 10 Bohr for heavy atoms (Na-Cl 17 Bohr); before A16 those pairs entered the GPU CN.
    const double r2 = dx*dx + dy*dy + dz*dz;
    if (r2 > r2_max) return;
    double r = __dsqrt_rn(r2);
    if (r < 1e-10) return;

    // erf-based CN contribution (same formula as k_cn_compute)
    double arg = kn * (r / rcov_sum[tid] - 1.0);
    double contrib = 0.5 * (1.0 + erf(arg));

    // Symmetric: both atoms i and j receive the contribution
    fix_add(&cn_raw[i], contrib, GFNFF_FIX_CN);   // order-independent: the CN decides the charges
    fix_add(&cn_raw[j], contrib, GFNFF_FIX_CN);
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
// WP-A (Jun 2026): on-device D4 dispersion pair-list build
// Reference (bit-identical target):
//   d4param_generator.cpp:GenerateDispersionPairsNative (line 1789)
//   getChargeWeightedC6 (line 1223) — Σ gw_i·gw_j·c6ref (CN-only weighting)
//   GFNFFParameters::zetaChargeScale (gfnff_par.h:933)
// ============================================================================

/// Device port of GFNFFParameters::zetaChargeScale (gfnff_par.h:933).
/// Caller passes the per-element zeff/c; Z out of [1,86] -> 1.0 at the call site.
__device__ inline double zetaChargeScaleDev(double zeff, double c, double q)
{
    double qmod = zeff + q;
    if (qmod < 0.0) return exp(3.0);            // highly anionic -> max zeta
    return exp(3.0 * (1.0 - exp(c * (1.0 - zeff / qmod))));
}

/// Map a linear upper-triangle index tid -> (i,j), i<j (same arithmetic as the
/// CN pair kernels so the surviving pair *set* matches).
__device__ inline void dispPairFromTid(int tid, int N, int& i, int& j)
{
    i = static_cast<int>((2.0 * N - 1.0 - sqrt((2.0 * N - 1.0) * (2.0 * N - 1.0) - 8.0 * tid)) * 0.5);
    while (i > 0 && i * (2 * N - i - 1) / 2 > tid) --i;
    while (i + 1 < N && (i + 1) * (2 * N - (i + 1) - 1) / 2 <= tid) ++i;
    j = tid - i * (2 * N - i - 1) / 2 + i + 1;
}

__global__ GFNFF_KERNEL_BOUNDS void k_disp_pairs_count(
    int N,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const int*    __restrict__ atom_types,
    double        cutoff_sq,
    int*          d_count)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int n_total = N * (N - 1) / 2;
    if (tid >= n_total) return;

    int i, j;
    dispPairFromTid(tid, N, i, j);
    if (j >= N) return;

    int zi = atom_types[i];
    int zj = atom_types[j];
    if (zi < 1 || zi > D4_MAX_ELEM || zj < 1 || zj > D4_MAX_ELEM) return;

    double dx = cx[i] - cx[j];
    double dy = cy[i] - cy[j];
    double dz = cz[i] - cz[j];
    applyMIC(dx, dy, dz);   // Nov 2026 (S1): pair once at its minimum-image distance (host generator, A3b)
    double r2 = dx * dx + dy * dy + dz * dz;
    if (r2 > cutoff_sq) return;

    atomicAdd(&d_count[i], 1);   // per-row count (Oct 2026): rows are placed in order by the build pass
}

__global__ GFNFF_KERNEL_BOUNDS void k_disp_pairs_build(
    int N,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const int*    __restrict__ atom_types,
    double        cutoff_sq,
    const double* __restrict__ gw,
    const double* __restrict__ c6_flat,
    const double* __restrict__ sqrtzr4r2,
    const double* __restrict__ topo_q,
    const double* __restrict__ zeta_zeff,
    const double* __restrict__ zeta_c,
    double        a1,
    double        a2,
    double        r_cut_eval,
    const int*    __restrict__ row_off,
    int*          d_counter,
    int*          out_i,
    int*          out_j,
    double*       out_C6,
    double*       out_r4r2,
    double*       out_r0sq,
    double*       out_zeta,
    double*       out_rcut)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int n_total = N * (N - 1) / 2;
    if (tid >= n_total) return;

    int i, j;
    dispPairFromTid(tid, N, i, j);
    if (j >= N) return;

    int zi = atom_types[i];
    int zj = atom_types[j];
    if (zi < 1 || zi > D4_MAX_ELEM || zj < 1 || zj > D4_MAX_ELEM) return;

    double dx = cx[i] - cx[j];
    double dy = cy[i] - cy[j];
    double dz = cz[i] - cz[j];
    applyMIC(dx, dy, dz);   // Nov 2026 (S1): same test as k_disp_pairs_count
    double r2 = dx * dx + dy * dy + dz * dz;
    if (r2 > cutoff_sq) return;

    // C6 = Σ_{ri,rj} gw_i[ri]·gw_j[rj]·c6ref  (== getChargeWeightedC6, CN-only weighting)
    int ei = zi - 1, ej = zj - 1;
    int nri = d_refn_const[ei];
    int nrj = d_refn_const[ej];
    size_t c6_base = (size_t)ei * D4_MAX_ELEM * D4_MAX_REF * D4_MAX_REF
                   + (size_t)ej * D4_MAX_REF * D4_MAX_REF;
    int gi = i * D4_MAX_REF;
    int gj = j * D4_MAX_REF;
    double c6 = 0.0;
    for (int ri = 0; ri < nri; ++ri) {
        double wi = gw[gi + ri];
        size_t b = c6_base + (size_t)ri * D4_MAX_REF;
        for (int rj = 0; rj < nrj; ++rj) {
            double c6ref = c6_flat[b + rj];
            if (fabs(c6ref) < 1e-20) continue;
            c6 += wi * gw[gj + rj] * c6ref;
        }
    }

    double si = sqrtzr4r2[ei];
    double sj = sqrtzr4r2[ej];
    double r4r2ij = 3.0 * si * sj;
    double t = a1 * sqrt(r4r2ij) + a2;
    double r0_sq = t * t;

    double qi = topo_q[i];
    double qj = topo_q[j];
    double zsi = (zi >= 1 && zi <= 86) ? zetaChargeScaleDev(zeta_zeff[zi - 1], zeta_c[zi - 1], qi) : 1.0;
    double zsj = (zj >= 1 && zj <= 86) ? zetaChargeScaleDev(zeta_zeff[zj - 1], zeta_c[zj - 1], qj) : 1.0;

    // Oct 2026: grouped by row i (row_off from an exclusive scan of the row counts), so that the list is
    // ordered like the host list; an arrival-ordered atomic slot scattered the pairs and made
    // k_dc6dcn_per_pair 1.9x and k_dispersion 1.6x slower (non-coalesced gw/coordinate loads).
    int slot = row_off[i] + atomicAdd(&d_counter[i], 1);
    out_i[slot]    = i;
    out_j[slot]    = j;
    out_C6[slot]   = c6;
    out_r4r2[slot] = r4r2ij;
    out_r0sq[slot] = r0_sq;
    out_zeta[slot] = zsi * zsj;
    out_rcut[slot] = r_cut_eval;  // D4 pair-energy cutoff: GFNFF::dispersionRadiiBohr().eval (50 without a cell; S1)
}

// ---------------------------------------------------------------------------
// L1 (Nov 2026, Claude Generated): ordered D4 pair build, one block per row i. The same pair test and the same
// per-pair statements as k_disp_pairs_count / k_disp_pairs_build, but each block walks j = i+1 .. N-1 in tiles and
// places the surviving pairs by a block prefix sum: the list comes out in (i, j) order without the key sort, and
// the row count needs no atomics. Selected by FFWorkspaceGPU::setListFastRebuild (default on).
// ---------------------------------------------------------------------------
#include <cub/block/block_reduce.cuh>
#include <cub/block/block_scan.cuh>

/// L1 (Nov 2026, Claude Generated): float minimum image of one orthorhombic component (pre-filter only). A difference
/// beyond 2 L returns NaN: every float test then fails to reject and the double path decides (keeps the float error
/// below 6e-8 (|d| + L) per component, see hb_detect_body).
__device__ __forceinline__ float hbWrapF(float d, float L, float iL)
{
    return fabsf(d) <= 2.f * L ? d - L * rintf(d * iL) : __int_as_float(0x7fc00000);
}

constexpr int kDispRowThreads = 256;

/// Float cell of the row kernels: orthorhombic lengths (off-diagonal noise < 1e-9 Bohr) up to 400 Bohr, else off.
struct DispRowCell { bool ok; float Lx, Ly, Lz, iLx, iLy, iLz; };
__device__ __forceinline__ DispRowCell dispRowCell()
{
    DispRowCell c;
    c.ok = d_has_pbc && d_lattice[0] <= 400.0 && d_lattice[4] <= 400.0 && d_lattice[8] <= 400.0
           && fabs(d_lattice[1]) < 1e-9 && fabs(d_lattice[2]) < 1e-9 && fabs(d_lattice[3]) < 1e-9
           && fabs(d_lattice[5]) < 1e-9 && fabs(d_lattice[6]) < 1e-9 && fabs(d_lattice[7]) < 1e-9;
    c.Lx = static_cast<float>(d_lattice[0]); c.Ly = static_cast<float>(d_lattice[4]); c.Lz = static_cast<float>(d_lattice[8]);
    c.iLx = 1.f / c.Lx; c.iLy = 1.f / c.Ly; c.iLz = 1.f / c.Lz;
    return c;
}

/// The pair test of k_disp_pairs_count (r^2 <= cutoff^2 at the double minimum image). With a cell a float minimum
/// image of the double differences decides first where it is safe: its r^2 is within about 2e-3 Bohr^2 of the double
/// value (|d| <= 2 L, L <= 400 Bohr; hbWrapF), so outside the band cutoff^2 (1 +- 1e-4) (0.1 Bohr^2 at 31 Bohr) it
/// gives the double answer; inside the band, or with NaN, the double test runs as before.
__device__ __forceinline__ bool dispRowPairIn(const double* __restrict__ cx, const double* __restrict__ cy,
    const double* __restrict__ cz, const int* __restrict__ atom_types, int i, int j, double cutoff_sq,
    const DispRowCell& fc, float cut_lo, float cut_hi)
{
    int zj = atom_types[j];
    if (zj < 1 || zj > D4_MAX_ELEM) return false;
    double dx = cx[i] - cx[j];
    double dy = cy[i] - cy[j];
    double dz = cz[i] - cz[j];
    if (fc.ok) {
        const float fx = hbWrapF(static_cast<float>(dx), fc.Lx, fc.iLx);
        const float fy = hbWrapF(static_cast<float>(dy), fc.Ly, fc.iLy);
        const float fz = hbWrapF(static_cast<float>(dz), fc.Lz, fc.iLz);
        const float r2f = fx * fx + fy * fy + fz * fz;
        if (r2f > cut_hi) return false;
        if (r2f < cut_lo) return true;
    }
    applyMIC(dx, dy, dz);
    double r2 = dx * dx + dy * dy + dz * dz;
    return !(r2 > cutoff_sq);
}

/// Per-atom charge scaling of the D4 pairs (the zsi / zsj of k_disp_pairs_build, same function and arguments).
__global__ void k_disp_atom_zeta(int N, const int* __restrict__ atom_types, const double* __restrict__ topo_q,
    const double* __restrict__ zeta_zeff, const double* __restrict__ zeta_c, double* __restrict__ zeta_atom)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    const int zi = atom_types[i];
    const double qi = topo_q[i];
    zeta_atom[i] = (zi >= 1 && zi <= 86) ? zetaChargeScaleDev(zeta_zeff[zi - 1], zeta_c[zi - 1], qi) : 1.0;
}

__global__ __launch_bounds__(kDispRowThreads) void k_disp_rows_count(
    int N, const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const int* __restrict__ atom_types, double cutoff_sq, int* __restrict__ row_count)
{
    using Reduce = cub::BlockReduce<int, kDispRowThreads>;
    __shared__ typename Reduce::TempStorage tmp;
    const int i = blockIdx.x;
    if (i >= N) return;
    const int zi = atom_types[i];
    const bool row_ok = !(zi < 1 || zi > D4_MAX_ELEM);
    const DispRowCell fc = dispRowCell();
    const float cut_lo = static_cast<float>(cutoff_sq) * (1.f - 1e-4f), cut_hi = static_cast<float>(cutoff_sq) * (1.f + 1e-4f);
    int c = 0;
    if (row_ok)
        for (int j = i + 1 + threadIdx.x; j < N; j += kDispRowThreads)
            if (dispRowPairIn(cx, cy, cz, atom_types, i, j, cutoff_sq, fc, cut_lo, cut_hi)) ++c;
    const int total = Reduce(tmp).Sum(c);
    if (threadIdx.x == 0) row_count[i] = total;
}

__global__ __launch_bounds__(kDispRowThreads) void k_disp_rows_build(
    int N, const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const int* __restrict__ atom_types, double cutoff_sq,
    const double* __restrict__ gw, const double* __restrict__ c6_flat, const double* __restrict__ sqrtzr4r2,
    const double* __restrict__ zeta_atom,
    double a1, double a2, double r_cut_eval, const int* __restrict__ row_off,
    int* out_i, int* out_j, double* out_C6, double* out_r4r2, double* out_r0sq, double* out_zeta, double* out_rcut)
{
    using Scan = cub::BlockScan<int, kDispRowThreads>;
    __shared__ typename Scan::TempStorage tmp;
    const int i = blockIdx.x;
    if (i >= N) return;
    const int zi = atom_types[i];
    if (zi < 1 || zi > D4_MAX_ELEM) return;   // uniform per block: no thread reaches the scan
    const DispRowCell fc = dispRowCell();
    const float cut_lo = static_cast<float>(cutoff_sq) * (1.f - 1e-4f), cut_hi = static_cast<float>(cutoff_sq) * (1.f + 1e-4f);
    int base = row_off[i];
    for (int j0 = i + 1; j0 < N; j0 += kDispRowThreads) {
        const int j = j0 + threadIdx.x;
        const int flag = (j < N && dispRowPairIn(cx, cy, cz, atom_types, i, j, cutoff_sq, fc, cut_lo, cut_hi)) ? 1 : 0;
        int pos = 0, tile_total = 0;
        Scan(tmp).ExclusiveSum(flag, pos, tile_total);
        if (flag) {
            const int zj = atom_types[j];
            // statements of k_disp_pairs_build (C6 contraction, r4r2ij, R0^2, zeta), same order
            int ei = zi - 1, ej = zj - 1;
            int nri = d_refn_const[ei];
            int nrj = d_refn_const[ej];
            size_t c6_base = (size_t)ei * D4_MAX_ELEM * D4_MAX_REF * D4_MAX_REF
                           + (size_t)ej * D4_MAX_REF * D4_MAX_REF;
            int gi = i * D4_MAX_REF;
            int gj = j * D4_MAX_REF;
            double c6 = 0.0;
            for (int ri = 0; ri < nri; ++ri) {
                double wi = gw[gi + ri];
                size_t b = c6_base + (size_t)ri * D4_MAX_REF;
                for (int rj = 0; rj < nrj; ++rj) {
                    double c6ref = c6_flat[b + rj];
                    if (fabs(c6ref) < 1e-20) continue;
                    c6 += wi * gw[gj + rj] * c6ref;
                }
            }
            double si = sqrtzr4r2[ei];
            double sj = sqrtzr4r2[ej];
            double r4r2ij = 3.0 * si * sj;
            double t = a1 * sqrt(r4r2ij) + a2;
            double r0_sq = t * t;
            double zsi = zeta_atom[i];   // k_disp_atom_zeta: the same zetaChargeScaleDev values
            double zsj = zeta_atom[j];
            const int slot = base + pos;
            out_i[slot]    = i;
            out_j[slot]    = j;
            out_C6[slot]   = c6;
            out_r4r2[slot] = r4r2ij;
            out_r0sq[slot] = r0_sq;
            out_zeta[slot] = zsi * zsj;
            out_rcut[slot] = r_cut_eval;
        }
        base += tile_total;
        __syncthreads();   // the scan storage is reused by the next tile
    }
}

// ---------------------------------------------------------------------------
// Non-bonded repulsion pair list on the device (Claude Generated, Oct 2026)
// ---------------------------------------------------------------------------
__device__ inline int repKeyFind(const long long* keys, int n, long long key)
{
    int lo = 0, hi = n - 1;
    while (lo <= hi) {
        const int mid = (lo + hi) >> 1;
        const long long v = keys[mid];
        if (v == key) return mid;
        if (v < key) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}

/// Pair survives the distance test, has valid elements and is not an excluded bond (i < j).
__device__ inline bool repPairWanted(int N, int i, int j, const double* cx, const double* cy, const double* cz,
                                     double cutoff_sq, const RepDevTables& t)
{
    if (!t.valid[i] || !t.valid[j]) return false;
    double dx = cx[i] - cx[j], dy = cy[i] - cy[j], dz = cz[i] - cz[j];
    applyMIC(dx, dy, dz);   // Nov 2026 (S1): minimum-image build distance (host periodic SpatialCellList, A4a)
    if (dx * dx + dy * dy + dz * dz > cutoff_sq) return false;
    return repKeyFind(t.bond_keys, t.n_bond_keys, static_cast<long long>(i) * N + j) < 0;
}

__global__ GFNFF_KERNEL_BOUNDS void k_rep_pairs_count(
    int N, const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    double cutoff_sq, RepDevTables t, int* d_count)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N * (N - 1) / 2) return;
    int i, j;
    dispPairFromTid(tid, N, i, j);
    if (j >= N) return;
    if (repPairWanted(N, i, j, cx, cy, cz, cutoff_sq, t)) atomicAdd(&d_count[i], 1);
}

__global__ GFNFF_KERNEL_BOUNDS void k_rep_pairs_build(
    int N, const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const int* __restrict__ atom_types, double cutoff_sq, RepDevTables t,
    const int* __restrict__ row_off, int* d_counter,
    int* out_i, int* out_j, double* out_alpha, double* out_repab, double* out_rcut)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N * (N - 1) / 2) return;
    int i, j;
    dispPairFromTid(tid, N, i, j);
    if (j >= N) return;
    if (!repPairWanted(N, i, j, cx, cy, cz, cutoff_sq, t)) return;

    // element-pair factor, same branch order as make_nb_rep()
    const int zi = atom_types[i], zj = atom_types[j];
    double ff = 1.0;
    if (zi == 1 && zj == 1) {
        ff = t.hhfac;
        const int k = repKeyFind(t.hh_keys, t.n_hh_keys, static_cast<long long>(i) * N + j);
        if (k >= 0) ff *= (t.hh_dist[k] == 2) ? t.hh13 : t.hh14;
    } else if ((zi == 1 && t.metal[j]) || (zj == 1 && t.metal[i])) {
        ff = 0.85;
    } else if ((zi == 1 && zj == 6) || (zj == 1 && zi == 6)) {
        ff = 0.91;
    } else if ((zi == 1 && zj == 8) || (zj == 1 && zi == 8)) {
        ff = 1.04;
    }
    const int slot = row_off[i] + atomicAdd(&d_counter[i], 1);
    out_i[slot] = i;
    out_j[slot] = j;
    out_alpha[slot] = sqrt(t.dum[i] * t.dum[j]) * ff;
    out_repab[slot] = t.repz[i] * t.repz[j] * t.repscaln;
    out_rcut[slot] = t.r_cut;
}

// L1 (Nov 2026, Claude Generated): ordered repulsion pair build with a cell, one block per row (as k_disp_rows_*):
// the pair test of repPairWanted with the float band pre-filter of dispRowPairIn, the pair values of
// k_rep_pairs_build in the same statements; the list comes out in (i, j) order (the two-pass build places the
// pairs of a row in atomic arrival order; the energy sums are fixed-point, so the order does not change a bit).
__device__ __forceinline__ bool repRowPairWanted(int N, int i, int j, const double* cx, const double* cy, const double* cz,
                                                 double cutoff_sq, const RepDevTables& t, const DispRowCell& fc,
                                                 float cut_lo, float cut_hi)
{
    if (!t.valid[i] || !t.valid[j]) return false;
    double dx = cx[i] - cx[j], dy = cy[i] - cy[j], dz = cz[i] - cz[j];
    int decided = 0;   // +1 in, -1 out, 0 double test
    if (fc.ok) {
        const float fx = hbWrapF(static_cast<float>(dx), fc.Lx, fc.iLx);
        const float fy = hbWrapF(static_cast<float>(dy), fc.Ly, fc.iLy);
        const float fz = hbWrapF(static_cast<float>(dz), fc.Lz, fc.iLz);
        const float r2f = fx * fx + fy * fy + fz * fz;
        if (r2f > cut_hi) decided = -1;
        else if (r2f < cut_lo) decided = 1;
    }
    if (decided < 0) return false;
    if (decided == 0) {
        applyMIC(dx, dy, dz);
        if (dx * dx + dy * dy + dz * dz > cutoff_sq) return false;
    }
    return repKeyFind(t.bond_keys, t.n_bond_keys, static_cast<long long>(i) * N + j) < 0;
}

__global__ __launch_bounds__(kDispRowThreads) void k_rep_rows_count(
    int N, const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    double cutoff_sq, RepDevTables t, int* __restrict__ row_count)
{
    using Reduce = cub::BlockReduce<int, kDispRowThreads>;
    __shared__ typename Reduce::TempStorage tmp;
    const int i = blockIdx.x;
    if (i >= N) return;
    const DispRowCell fc = dispRowCell();
    const float cut_lo = static_cast<float>(cutoff_sq) * (1.f - 1e-4f), cut_hi = static_cast<float>(cutoff_sq) * (1.f + 1e-4f);
    int c = 0;
    for (int j = i + 1 + threadIdx.x; j < N; j += kDispRowThreads)
        if (repRowPairWanted(N, i, j, cx, cy, cz, cutoff_sq, t, fc, cut_lo, cut_hi)) ++c;
    const int total = Reduce(tmp).Sum(c);
    if (threadIdx.x == 0) row_count[i] = total;
}

__global__ __launch_bounds__(kDispRowThreads) void k_rep_rows_build(
    int N, const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const int* __restrict__ atom_types, double cutoff_sq, RepDevTables t, const int* __restrict__ row_off,
    int* out_i, int* out_j, double* out_alpha, double* out_repab, double* out_rcut)
{
    using Scan = cub::BlockScan<int, kDispRowThreads>;
    __shared__ typename Scan::TempStorage tmp;
    const int i = blockIdx.x;
    if (i >= N) return;
    const DispRowCell fc = dispRowCell();
    const float cut_lo = static_cast<float>(cutoff_sq) * (1.f - 1e-4f), cut_hi = static_cast<float>(cutoff_sq) * (1.f + 1e-4f);
    int base = row_off[i];
    for (int j0 = i + 1; j0 < N; j0 += kDispRowThreads) {
        const int j = j0 + threadIdx.x;
        const int flag = (j < N && repRowPairWanted(N, i, j, cx, cy, cz, cutoff_sq, t, fc, cut_lo, cut_hi)) ? 1 : 0;
        int pos = 0, tile_total = 0;
        Scan(tmp).ExclusiveSum(flag, pos, tile_total);
        if (flag) {
            // element-pair factor, same branch order as make_nb_rep() (statements of k_rep_pairs_build)
            const int zi = atom_types[i], zj = atom_types[j];
            double ff = 1.0;
            if (zi == 1 && zj == 1) {
                ff = t.hhfac;
                const int k = repKeyFind(t.hh_keys, t.n_hh_keys, static_cast<long long>(i) * N + j);
                if (k >= 0) ff *= (t.hh_dist[k] == 2) ? t.hh13 : t.hh14;
            } else if ((zi == 1 && t.metal[j]) || (zj == 1 && t.metal[i])) {
                ff = 0.85;
            } else if ((zi == 1 && zj == 6) || (zj == 1 && zi == 6)) {
                ff = 0.91;
            } else if ((zi == 1 && zj == 8) || (zj == 1 && zi == 8)) {
                ff = 1.04;
            }
            const int slot = base + pos;
            out_i[slot] = i;
            out_j[slot] = j;
            out_alpha[slot] = sqrt(t.dum[i] * t.dum[j]) * ff;
            out_repab[slot] = t.repz[i] * t.repz[j] * t.repscaln;
            out_rcut[slot] = t.r_cut;
        }
        base += tile_total;
        __syncthreads();
    }
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

    fix_add(&hb_cn_per_atom[H], contrib, GFNFF_FIX_CN);
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

// Single-block twin of k_eeq_reduce_sums (Oct 2026): every thread sums a fixed strided slice, a fixed tree adds the
// slices - the same two sums in the same order in every run (the block-wise atomicAdd added in arrival order).
__global__ void k_eeq_reduce_sums_det(int N, const double* __restrict__ d_rhs, double* __restrict__ d_sums)
{
    extern __shared__ double sdata[];
    double* s_z1 = sdata;
    double* s_Z2 = sdata + blockDim.x;
    const int tid = threadIdx.x;
    double a = 0.0, b = 0.0;
    for (int i = tid; i < N; i += blockDim.x) { a += d_rhs[i]; b += d_rhs[N + i]; }
    s_z1[tid] = a; s_Z2[tid] = b;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) { s_z1[tid] += s_z1[tid + s]; s_Z2[tid] += s_Z2[tid + s]; }
        __syncthreads();
    }
    if (tid == 0) { d_sums[0] = s_z1[0]; d_sums[1] = s_Z2[0]; }
}

// Fragment sums over a CSR list of the atoms of each fragment (ascending atom index): one block per fragment,
// fixed strided slices and a fixed tree -> deterministic (k_frag_sum added in arrival order).
__global__ void k_frag_sum_csr(int nfrag, const double* __restrict__ v, const int* __restrict__ off,
                               const int* __restrict__ idx, double* __restrict__ out)
{
    extern __shared__ double sh[];
    const int f = blockIdx.x;
    if (f >= nfrag) return;
    const int tid = threadIdx.x;
    double a = 0.0;
    for (int p = off[f] + tid; p < off[f + 1]; p += blockDim.x) a += v[idx[p]];
    sh[tid] = a;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) sh[tid] += sh[tid + s];
        __syncthreads();
    }
    if (tid == 0) out[f] = sh[0];
}

// Warp per fragment (Oct 2026): for many small fragments (water boxes: 1500 blocks of 256 threads with 8 barriers each
// cost 21 us per call, twice per PCG iteration). Lane l sums atoms l, l+32, ... of the fragment in ascending order,
// a fixed shuffle tree adds the 32 partials: the same result in every run, independent of the launch shape.
__global__ void k_frag_sum_warp(int nfrag, const double* __restrict__ v, const int* __restrict__ off,
                                const int* __restrict__ idx, double* __restrict__ out)
{
    const int lane = threadIdx.x & 31;
    const int f = blockIdx.x * (blockDim.x >> 5) + (threadIdx.x >> 5);
    if (f >= nfrag) return;
    double a = 0.0;
    for (int p = off[f] + lane; p < off[f + 1]; p += 32) a += v[idx[p]];
    for (int o = 16; o > 0; o >>= 1) a += __shfl_down_sync(0xffffffffu, a, o);
    if (lane == 0) out[f] = a;
}

// Deterministic twin of k_eeq_reduce_fragment_sums: one thread per (fragment, column), atoms in ascending order.
// Column 0 -> Cz1[f], column g + 1 -> S[f * nfrag + g].
__global__ void k_eeq_reduce_fragment_sums_csr(int N, int nfrag, const double* __restrict__ d_rhs,
                                               const int* __restrict__ off, const int* __restrict__ idx,
                                               double* __restrict__ d_Cz1, double* __restrict__ d_S)
{
    const long long t = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    const long long ncol = static_cast<long long>(nfrag) + 1;
    if (t >= static_cast<long long>(nfrag) * ncol) return;
    const int f = static_cast<int>(t / ncol), c = static_cast<int>(t % ncol);
    double acc = 0.0;
    for (int p = off[f]; p < off[f + 1]; ++p) acc += d_rhs[static_cast<size_t>(c) * N + idx[p]];
    if (c == 0) d_Cz1[f] = acc; else d_S[static_cast<size_t>(f) * nfrag + (c - 1)] = acc;
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
// Used by EEQSolverGPU::solveSinglePCG.
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

// PCG init: r = b − A·x. d_Ax holds A·x (computed via cublasDsymv before this call).
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

// Direction update: p_out = z + β·p_in. Single-pass — avoids cublasDcopy + cublasDaxpy.
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
// WP7-D: GPU block-Jacobi preconditioner (Jun 2026, Claude Generated)
// ============================================================================

// Mirror the lower triangle of each per-fragment block into the upper triangle.
// Column-major: element (row, col) at col*Nf + row; potri(LOWER) fills row >= col.
__global__ void k_eeq_symmetrize_blocks(
    int           nfrag,
    double*       __restrict__ d_A_blocks,
    const int*    __restrict__ frag_sizes,
    const int*    __restrict__ frag_offsets_A)
{
    int f = blockIdx.x;
    if (f >= nfrag) return;
    int Nf = frag_sizes[f];
    if (Nf <= 0) return;
    double* A = d_A_blocks + frag_offsets_A[f];
    for (int idx = threadIdx.x; idx < Nf * Nf; idx += blockDim.x) {
        int col = idx / Nf;
        int row = idx % Nf;
        if (row < col)                       // upper entry (garbage after potri LOWER)
            A[col * Nf + row] = A[row * Nf + col];  // = lower partner (col,row)
    }
}

// z = blockdiag(A_ff^-1)·r. One thread-block per fragment; gather/scatter via the
// fragment-sorted atom map so r/z stay in global atom order (matches the CPU
// BlockJacobiPC::apply contract). Dynamic shared mem stages r_f (max_frag_N doubles).
__global__ void k_eeq_block_jacobi_apply(
    int           nfrag,
    const double* __restrict__ d_Ainv_blocks,
    const int*    __restrict__ frag_sizes,
    const int*    __restrict__ frag_offsets_A,
    const int*    __restrict__ frag_atom_offsets,
    const int*    __restrict__ frag_atom_map,
    const double* __restrict__ d_r,
    double*       __restrict__ d_z)
{
    extern __shared__ double s_rf[];   // [Nf] staged residual for this fragment
    int f = blockIdx.x;
    if (f >= nfrag) return;
    int Nf = frag_sizes[f];
    if (Nf <= 0) return;
    int atom_start    = frag_atom_offsets[f];
    const double* Ainv = d_Ainv_blocks + frag_offsets_A[f];

    // Stage r_f into shared memory (gather global → local).
    for (int li = threadIdx.x; li < Nf; li += blockDim.x)
        s_rf[li] = d_r[frag_atom_map[atom_start + li]];
    __syncthreads();

    // z_f = Ainv_f · r_f  (symmetric, column-major: Ainv[(col)*Nf + row]).
    for (int li = threadIdx.x; li < Nf; li += blockDim.x) {
        double acc = 0.0;
        for (int lj = 0; lj < Nf; ++lj)
            acc += Ainv[lj * Nf + li] * s_rf[lj];
        d_z[frag_atom_map[atom_start + li]] = acc;   // scatter local → global
    }
}

// ── WP7-E: GPU projected PCG kernels (Sep 2026) ──────────────────────────────
__global__ void k_frag_sum(int N,
                            const double* __restrict__ v,
                            const int*    __restrict__ atom_frag,
                            double*       __restrict__ frag_sum)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    atomicAdd(&frag_sum[atom_frag[i]], v[i]);
}

__global__ void k_frag_scatter_add(int N,
                                    double*       __restrict__ v,
                                    const int*    __restrict__ atom_frag,
                                    const double* __restrict__ frag_delta)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    v[i] += frag_delta[atom_frag[i]];
}

__global__ void k_frag_project_delta(int nfrag,
                                      const double* __restrict__ frag_sum,
                                      const double* __restrict__ inv_frag_atoms,
                                      double*       __restrict__ frag_delta)
{
    int f = blockIdx.x * blockDim.x + threadIdx.x;
    if (f >= nfrag) return;
    frag_delta[f] = -frag_sum[f] * inv_frag_atoms[f];
}

__global__ void k_frag_feasibility_delta(int nfrag,
                                          const double* __restrict__ frag_sum,
                                          const double* __restrict__ rhs_constraints,
                                          const double* __restrict__ inv_frag_atoms,
                                          double*       __restrict__ frag_delta)
{
    int f = blockIdx.x * blockDim.x + threadIdx.x;
    if (f >= nfrag) return;
    frag_delta[f] = (rhs_constraints[f] - frag_sum[f]) * inv_frag_atoms[f];
}

__global__ void k_frag_precond_correct(int N,
                                        const double* __restrict__ Minv,
                                        const int*    __restrict__ atom_frag,
                                        const double* __restrict__ frag_sum,
                                        const double* __restrict__ frag_sM,
                                        double*       __restrict__ z)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    int f = atom_frag[i];
    double sM = frag_sM[f];
    if (sM > 0.0) z[i] -= Minv[i] * frag_sum[f] / sM;
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
            const double asum = alpha[gi] + alpha[gj];
            // Far field (Claude Generated, Oct 2026): erf == 1.0 exactly, value bit-identical
            double val        = (r_sq >= GFNFF_COUL_FARFIELD_X2 * asum)
                                    ? 1.0 / r
                                    : erf((1.0 / sqrt(asum)) * r) / r;
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

// ============================================================================
// Device-side H-bond detection (Claude Generated, Oct 2026)
// Host reference: GFNFF::detectHydrogenBondsNative(), cell-list branch (process_pair / make_nhb2 and
// case1StrengthDecision / estimateHBStrengthCase1). Constants as in k_hbonds.
// ============================================================================

__device__ __forceinline__ bool d_hb_bonded(int a, int b, const int* __restrict__ off, const int* __restrict__ idx)
{
    for (int p = off[a]; p < off[a + 1]; ++p)
        if (idx[p] == b) return true;
    return false;
}

__device__ __forceinline__ double d_hb_dist2(const double* __restrict__ cx, const double* __restrict__ cy,
                                              const double* __restrict__ cz, int a, int b)
{
    double dx = cx[a] - cx[b], dy = cy[a] - cy[b], dz = cz[a] - cz[b];
    applyMIC(dx, dy, dz);   // Nov 2026 (S1): minimum-image search distances (host detection, A4b)
    return dx * dx + dy * dy + dz * dz;
}

/// |E_HB| of a case-1 triple, the formula of estimateHBStrengthCase1() (and of the energy kernel).
__device__ double d_hb_case1_strength(const double* __restrict__ cx, const double* __restrict__ cy,
                                       const double* __restrict__ cz, const int* __restrict__ Z,
                                       const double* __restrict__ bas, const double* __restrict__ aci,
                                       const double* __restrict__ q, int A, int H, int B, double thr)
{
    // Nov 2026 (S1): the kernel frame of the host estimate, A and B as the images next to H
    double ax, ay, az, bx, by, bz;
    micPlace(cx, cy, cz, A, cx[H], cy[H], cz[H], ax, ay, az);
    micPlace(cx, cy, cz, B, cx[H], cy[H], cz[H], bx, by, bz);
    const double r_AH = sqrt((cx[H] - ax) * (cx[H] - ax) + (cy[H] - ay) * (cy[H] - ay) + (cz[H] - az) * (cz[H] - az));
    const double r_HB = sqrt((bx - cx[H]) * (bx - cx[H]) + (by - cy[H]) * (by - cy[H]) + (bz - cz[H]) * (bz - cz[H]));
    const double r_AB = sqrt((bx - ax) * (bx - ax) + (by - ay) * (by - ay) + (bz - az) * (bz - az));
    const double r_vdw = d_cov_radii[Z[A] - 1] + d_cov_radii[Z[B] - 1];
    // Upper bound first (case1StrengthDecision() of the host): no pow, rejects most candidates. Valid for
    // non-negative basicities and acidities; returns 0 when the bound is below the threshold.
    {
        const double expo_b = (49.0 / r_vdw) * ((r_AH + r_HB) / r_AB - 1.0);
        if (expo_b > 15.0) return 0.0;
        const double bA = bas[A], bB = bas[B], aA = aci[A], aB = aci[B];
        if (bA >= 0.0 && bB >= 0.0 && aA >= 0.0 && aB >= 0.0) {
            const double t = r_AB * r_AB / 85.0;
            const double t2 = t * t;
            const double outl_bound = (expo_b > 0.0) ? 2.0 * exp(-expo_b) : 2.0;
            const double bound = fmax(bA, bB) * fmax(aA, aB) / (1.0 + t2 * t2 * t2) / (r_AB * r_AB * r_AB) * outl_bound;
            if (isfinite(bound) && bound < thr * (1.0 - 1e-6)) return 0.0;
        }
    }
    const double r_AH_4 = r_AH * r_AH * r_AH * r_AH;
    const double r_HB_4 = r_HB * r_HB * r_HB * r_HB;
    const double denom_DA = 1.0 / (r_AH_4 + r_HB_4);
    const double Q_H = d_charge_scaling(q[H], 15.0, 1.0);
    const double Q_A = d_charge_scaling(-q[A], 15.0, 1.0);
    const double Q_B = d_charge_scaling(-q[B], 15.0, 1.0);
    const double b = (Q_A * bas[A] * r_AH_4 + Q_B * bas[B] * r_HB_4) * denom_DA;
    const double a = (aci[B] * r_AH_4 + aci[A] * r_HB_4) * denom_DA;
    const double ratio_s = 22.0 * r_vdw / (r_AB * r_AB);
    const double damp_short = 1.0 / (1.0 + pow(ratio_s, 6.0));
    const double damp_long = 1.0 / (1.0 + pow(r_AB * r_AB / 85.0, 6.0));
    const double expo = (49.0 / r_vdw) * ((r_AH + r_HB) / r_AB - 1.0);
    const double damp_outl = (expo > 15.0) ? 0.0 : 2.0 / (1.0 + exp(expo));
    const double rdamp = damp_short * damp_long / (r_AB * r_AB * r_AB);
    return fabs(b * a * rdamp * (Q_H * damp_outl));
}

__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_hb_to_float(
    int N, const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    float* fx, float* fy, float* fz)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= N) return;
    fx[a] = static_cast<float>(cx[a]); fy[a] = static_cast<float>(cy[a]); fz[a] = static_cast<float>(cz[a]);
}

/// FP32 evaluation of the case-1 strength of d_hb_case1_strength(): returns -1 (surely below thr), +1 (surely at
/// or above thr) or 0 (within a 1 % band around thr or a borderline out-of-line exponent: the caller runs the
/// double path). The float evaluation is good to about 1e-5 relative, the 1 % band covers it with a wide margin.
__device__ __forceinline__ int d_hb_strength_decision_f(
    const float* __restrict__ fx, const float* __restrict__ fy, const float* __restrict__ fz,
    const int* __restrict__ Z, const double* __restrict__ bas, const double* __restrict__ aci,
    const double* __restrict__ q, int A, int H, int B, float thr)
{
    if (d_has_pbc) return 0;                   // Nov 2026 (S1): plain float coordinates carry no image; the double path decides
    const float ax = fx[H] - fx[A], ay = fy[H] - fy[A], az = fz[H] - fz[A];
    const float hx = fx[B] - fx[H], hy = fy[B] - fy[H], hz = fz[B] - fz[H];
    const float bx = fx[B] - fx[A], by = fy[B] - fy[A], bz = fz[B] - fz[A];
    const float r_AH2 = ax * ax + ay * ay + az * az;
    const float r_HB2 = hx * hx + hy * hy + hz * hz;
    const float r_AB2 = bx * bx + by * by + bz * bz;
    const float r_AH = sqrtf(r_AH2), r_HB = sqrtf(r_HB2), r_AB = sqrtf(r_AB2);
    const float r_vdw = static_cast<float>(d_cov_radii[Z[A] - 1] + d_cov_radii[Z[B] - 1]);
    const float expo = (49.f / r_vdw) * ((r_AH + r_HB) / r_AB - 1.f);
    if (expo > 15.01f) return -1;              // damping is exactly 0 beyond 15
    if (expo > 14.99f) return 0;               // borderline: the double path decides
    const float bA = static_cast<float>(bas[A]), bB = static_cast<float>(bas[B]);
    const float aA = static_cast<float>(aci[A]), aB = static_cast<float>(aci[B]);
    const float r_AH4 = r_AH2 * r_AH2, r_HB4 = r_HB2 * r_HB2;
    const float denom = 1.f / (r_AH4 + r_HB4);
    auto cs = [](float qq) { const float e = expf(15.f * qq); return e / (e + 1.f); };
    const float Q_H = cs(static_cast<float>(q[H]));
    const float Q_A = cs(-static_cast<float>(q[A]));
    const float Q_B = cs(-static_cast<float>(q[B]));
    const float b = (Q_A * bA * r_AH4 + Q_B * bB * r_HB4) * denom;
    const float a = (aB * r_AH4 + aA * r_HB4) * denom;
    const float rs = 22.f * r_vdw / r_AB2;
    const float rs2 = rs * rs;
    const float damp_short = 1.f / (1.f + rs2 * rs2 * rs2);
    const float t = r_AB2 / 85.f;
    const float t2 = t * t;
    const float damp_long = 1.f / (1.f + t2 * t2 * t2);
    const float damp_outl = 2.f / (1.f + expf(expo));
    const float est = fabsf(b * a * (damp_short * damp_long / (r_AB2 * r_AB)) * (Q_H * damp_outl));
    if (!isfinite(est)) return 0;
    if (est < 0.99f * thr) return -1;
    if (est > 1.01f * thr) return 1;
    return 0;
}

/// L1 (Nov 2026, Claude Generated): d_hb_strength_decision_f with the three displacements passed in (H - A, B - H,
/// B - A as float minimum images of the double differences, periodic cell); the statements after them are those
/// of d_hb_strength_decision_f, so the decision band (1 %) and the expo test are the same.
__device__ __forceinline__ int d_hb_strength_decision_fv(
    float ax, float ay, float az, float hx, float hy, float hz, float bx, float by, float bz,
    const int* __restrict__ Z, const double* __restrict__ bas, const double* __restrict__ aci,
    const double* __restrict__ q, int A, int H, int B, float thr)
{
    const float r_AH2 = ax * ax + ay * ay + az * az;
    const float r_HB2 = hx * hx + hy * hy + hz * hz;
    const float r_AB2 = bx * bx + by * by + bz * bz;
    const float r_AH = sqrtf(r_AH2), r_HB = sqrtf(r_HB2), r_AB = sqrtf(r_AB2);
    const float r_vdw = static_cast<float>(d_cov_radii[Z[A] - 1] + d_cov_radii[Z[B] - 1]);
    const float expo = (49.f / r_vdw) * ((r_AH + r_HB) / r_AB - 1.f);
    if (expo > 15.01f) return -1;              // damping is exactly 0 beyond 15
    if (expo > 14.99f) return 0;               // borderline: the double path decides
    const float bA = static_cast<float>(bas[A]), bB = static_cast<float>(bas[B]);
    const float aA = static_cast<float>(aci[A]), aB = static_cast<float>(aci[B]);
    const float r_AH4 = r_AH2 * r_AH2, r_HB4 = r_HB2 * r_HB2;
    const float denom = 1.f / (r_AH4 + r_HB4);
    auto cs = [](float qq) { const float e = expf(15.f * qq); return e / (e + 1.f); };
    const float Q_H = cs(static_cast<float>(q[H]));
    const float Q_A = cs(-static_cast<float>(q[A]));
    const float Q_B = cs(-static_cast<float>(q[B]));
    const float b = (Q_A * bA * r_AH4 + Q_B * bB * r_HB4) * denom;
    const float a = (aB * r_AH4 + aA * r_HB4) * denom;
    const float rs = 22.f * r_vdw / r_AB2;
    const float rs2 = rs * rs;
    const float damp_short = 1.f / (1.f + rs2 * rs2 * rs2);
    const float t = r_AB2 / 85.f;
    const float t2 = t * t;
    const float damp_long = 1.f / (1.f + t2 * t2 * t2);
    const float damp_outl = 2.f / (1.f + expf(expo));
    const float est = fabsf(b * a * (damp_short * damp_long / (r_AB2 * r_AB)) * (Q_H * damp_outl));
    if (!isfinite(est)) return 0;
    if (est < 0.99f * thr) return -1;
    if (est > 1.01f * thr) return 1;
    return 0;
}

__global__ GFNFF_KERNEL_BOUNDS void k_hb_ab_pairs(
    int N, const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const double* __restrict__ bas, const double* __restrict__ aci, const char* __restrict__ ab_ok,
    double thr1, int cap, int* d_counter, int* ab_i, int* ab_j)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int n_total = N * (N - 1) / 2;
    if (tid >= n_total) return;
    int i, j;
    dispPairFromTid(tid, N, i, j);
    if (j >= N) return;
    if (!ab_ok[i] || !ab_ok[j]) return;
    if (d_hb_dist2(cx, cy, cz, i, j) > thr1) return;
    const double s1 = bas[i] * aci[j], s2 = bas[j] * aci[i];
    if (s1 < 1e-6 && s2 < 1e-6) return;
    const int pos = atomicAdd(d_counter, 1);
    if (pos < cap) { ab_i[pos] = i; ab_j[pos] = j; }
}


/// Body of the H-bond detection (one warp per AB pair). PBCF = false is the former kernel statement for statement
/// (k_hb_detect); PBCF = true adds a float minimum-image pre-filter for a periodic cell (k_hb_detect_pbc). The
/// pre-filter only rejects: a triple that passes the double tests (r_iH^2 <= thr2, r_AB^2 + r_iH^2 + r_jH^2 < thr2)
/// has a float estimate within about 5e-3 Bohr^2 of the double value (double differences rounded once to float,
/// |d| <= 2 L, L <= 400 Bohr), inside the 1e-4 thr2 margin (0.045 Bohr^2 for water), and a near-minimal float
/// image differs from the double one only where both are near L/2, beyond the search radius (< L_min/2). So the
/// emitted set, and after the key sort the list, equal the double-only test. The case-1 strength filter uses the float
/// decision on the same minimum-image displacements (1 % band, double path inside it), as without a cell.
template <bool PBCF>
__device__ __forceinline__ void hb_detect_body(
    int n_ab, const int* __restrict__ ab_i, const int* __restrict__ ab_j,
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const float* __restrict__ fx, const float* __restrict__ fy, const float* __restrict__ fz, int N,
    const int* __restrict__ Z, const double* __restrict__ bas, const double* __restrict__ aci,
    const double* __restrict__ q,
    const int* __restrict__ hyd_list, int n_hyd,
    const int* __restrict__ bond_off, const int* __restrict__ bond_idx,
    const int* __restrict__ hyd_off, const int* __restrict__ hyd_idx,
    const int* __restrict__ nbr_off, const int* __restrict__ nbr_idx,
    double thr1, double thr2, double min_pair_energy,
    int cap, int* d_counter, long long* out_key, int* out_ct, int* out_par,
    const int* __restrict__ hnb_off, const int* __restrict__ hnb_idx)
{
    const int warp = (blockIdx.x * blockDim.x + threadIdx.x) >> 5;
    const int lane = threadIdx.x & 31;
    if (warp >= n_ab) return;
    const int i = ab_i[warp], j = ab_j[warp];
    const double r_AB_sq = d_hb_dist2(cx, cy, cz, i, j);

    const long long NN = N;
    auto emit = [&](int A, int H, int B, int ct, int par) {
        const int pos = atomicAdd(d_counter, 1);
        if (pos < cap) {
            out_key[pos] = (static_cast<long long>(A) * NN + H) * NN + B;
            out_ct[pos] = ct;
            out_par[pos] = par;
        }
    };

    // nhb2: hydrogen bonded to the donor (make_nhb2): case 2, 3 (carbonyl-like acceptor O) or 4 (two-coordinate N)
    if (lane == 0 && !d_hb_bonded(i, j, bond_off, bond_idx)) {
        auto emit_nhb2 = [&](int donor, int H, int acc) {
            int ct = 2, par = -1;
            const int deg = nbr_off[acc + 1] - nbr_off[acc];
            if (Z[acc] == 8 && deg == 1) {
                const int parent = nbr_idx[nbr_off[acc]];
                if (Z[parent] == 6 || Z[parent] == 7) { ct = 3; par = parent; }
            } else if (Z[acc] == 7 && deg == 2) {
                ct = 4;
            }
            emit(donor, H, acc, ct, par);
        };
        for (int p = hyd_off[i]; p < hyd_off[i + 1]; ++p) emit_nhb2(i, hyd_idx[p], j);
        for (int p = hyd_off[j]; p < hyd_off[j + 1]; ++p) {
            const int H = hyd_idx[p];
            if (d_hb_bonded(i, H, bond_off, bond_idx)) continue;  // bonded to i takes precedence
            emit_nhb2(j, H, i);
        }
    }

    // case 1: hydrogen bonded to neither, sum-of-squared-distances criterion (gfnff_ini2.f90:742). The lanes
    // split the HB hydrogens; a float pre-filter (relative margin 1e-4, far above the float rounding of the
    // squared distances) rejects most of them before the double test, which decides as on the host.
    const float fxi = fx[i], fyi = fy[i], fzi = fz[i], fxj = fx[j], fyj = fy[j], fzj = fz[j];
    const float f_rab = static_cast<float>(r_AB_sq);
    const float f_reject = static_cast<float>(thr2) * 1.0001f;
    // L1 (Nov 2026): orthorhombic cell lengths for the float pre-filter of the PBC instantiation; off for a
    // non-diagonal lattice or L > 400 Bohr (the double test alone decides, as before).
    [[maybe_unused]] bool pbc_filter = false;
    [[maybe_unused]] float Lx = 0.f, Ly = 0.f, Lz = 0.f, iLx = 0.f, iLy = 0.f, iLz = 0.f;
    [[maybe_unused]] double cxi = 0.0, cyi = 0.0, czi = 0.0, cxj = 0.0, cyj = 0.0, czj = 0.0;
    if constexpr (PBCF) {
        // off-diagonal entries are rounding noise of the cell conversion (5e-15 Bohr seen); above 1e-9 Bohr the
        // cell is treated as non-orthorhombic and the filter stays off
        pbc_filter = d_has_pbc && d_lattice[0] <= 400.0 && d_lattice[4] <= 400.0 && d_lattice[8] <= 400.0
                     && fabs(d_lattice[1]) < 1e-9 && fabs(d_lattice[2]) < 1e-9 && fabs(d_lattice[3]) < 1e-9
                     && fabs(d_lattice[5]) < 1e-9 && fabs(d_lattice[6]) < 1e-9 && fabs(d_lattice[7]) < 1e-9;
        Lx = static_cast<float>(d_lattice[0]); Ly = static_cast<float>(d_lattice[4]); Lz = static_cast<float>(d_lattice[8]);
        iLx = 1.f / Lx; iLy = 1.f / Ly; iLz = 1.f / Lz;
        cxi = cx[i]; cyi = cy[i]; czi = cz[i]; cxj = cx[j]; cyj = cy[j]; czj = cz[j];
    }
    // L1 (Nov 2026): the PBC instantiation walks the hydrogens near i (k_hb_hnb_*: a superset of every H with
    // r_iH^2 <= thr2, the first test below); otherwise every HB hydrogen as before.
    int h_begin = 0, h_end = n_hyd;
    if constexpr (PBCF) {
        if (hnb_off) { h_begin = hnb_off[i]; h_end = hnb_off[i + 1]; }
    }
    for (int hq = h_begin + lane; hq < h_end; hq += 32) {
        int H;
        if constexpr (PBCF) H = hnb_off ? hnb_idx[hq] : hyd_list[hq];
        else H = hyd_list[hq];
        if (!d_has_pbc) {   // Nov 2026 (S1): the float pre-filter has no image; with a cell the double test alone decides
            const float hx = fx[H], hy = fy[H], hz = fz[H];
            const float ax = hx - fxi, ay = hy - fyi, az = hz - fzi;
            const float f_ih = ax * ax + ay * ay + az * az;
            if (f_ih > f_reject) continue;
            const float bx = hx - fxj, by = hy - fyj, bz = hz - fzj;
            if (f_rab + f_ih + (bx * bx + by * by + bz * bz) > f_reject) continue;
        }
        if constexpr (PBCF) {   // L1 (Nov 2026): float minimum-image pre-filter with a cell, see k_hb_detect_pbc
            if (pbc_filter) {
                const float ax = hbWrapF(static_cast<float>(cx[H] - cxi), Lx, iLx);
                const float ay = hbWrapF(static_cast<float>(cy[H] - cyi), Ly, iLy);
                const float az = hbWrapF(static_cast<float>(cz[H] - czi), Lz, iLz);
                const float f_ih = ax * ax + ay * ay + az * az;
                if (f_ih > f_reject) continue;
                const float bx = hbWrapF(static_cast<float>(cx[H] - cxj), Lx, iLx);
                const float by = hbWrapF(static_cast<float>(cy[H] - cyj), Ly, iLy);
                const float bz = hbWrapF(static_cast<float>(cz[H] - czj), Lz, iLz);
                if (f_rab + f_ih + (bx * bx + by * by + bz * bz) > f_reject) continue;
            }
        }
        if (H == i) continue;                                  // forEachNeighbor excludes the query atom
        const double r_iH_sq = d_hb_dist2(cx, cy, cz, i, H);
        if (r_iH_sq > thr2) continue;                          // cell query radius sqrt(thr2)
        const double r_jH_sq = d_hb_dist2(cx, cy, cz, j, H);
        if (!(r_AB_sq + r_iH_sq + r_jH_sq < thr2)) continue;
        if (d_hb_bonded(i, H, bond_off, bond_idx) || d_hb_bonded(j, H, bond_off, bond_idx)) continue;
        if (min_pair_energy > 0.0) {
            int dec;
            if constexpr (PBCF) {   // L1: the float decision on minimum-image displacements (d_hb_strength_decision_f returns 0 with a cell)
                if (pbc_filter) {
                    const float ax = hbWrapF(static_cast<float>(cx[H] - cxi), Lx, iLx);
                    const float ay = hbWrapF(static_cast<float>(cy[H] - cyi), Ly, iLy);
                    const float az = hbWrapF(static_cast<float>(cz[H] - czi), Lz, iLz);
                    const float hx = hbWrapF(static_cast<float>(cxj - cx[H]), Lx, iLx);
                    const float hy = hbWrapF(static_cast<float>(cyj - cy[H]), Ly, iLy);
                    const float hz = hbWrapF(static_cast<float>(czj - cz[H]), Lz, iLz);
                    const float bx = hbWrapF(static_cast<float>(cxj - cxi), Lx, iLx);
                    const float by = hbWrapF(static_cast<float>(cyj - cyi), Ly, iLy);
                    const float bz = hbWrapF(static_cast<float>(czj - czi), Lz, iLz);
                    dec = d_hb_strength_decision_fv(ax, ay, az, hx, hy, hz, bx, by, bz, Z, bas, aci, q, i, H, j,
                                                    static_cast<float>(min_pair_energy));
                } else {
                    dec = d_hb_strength_decision_f(fx, fy, fz, Z, bas, aci, q, i, H, j, static_cast<float>(min_pair_energy));
                }
            } else {
                dec = d_hb_strength_decision_f(fx, fy, fz, Z, bas, aci, q, i, H, j, static_cast<float>(min_pair_energy));
            }
            if (dec < 0) continue;
            if (dec == 0) {                                    // within the float margin: the exact (double) value decides
                const double est = d_hb_case1_strength(cx, cy, cz, Z, bas, aci, q, i, H, j, min_pair_energy);
                if (est < min_pair_energy) continue;           // NaN keeps the entry, as on the host
            }
        }
        emit(i, H, j, 1, -1);
    }
}

__global__ GFNFF_KERNEL_BOUNDS void k_hb_detect(
    int n_ab, const int* __restrict__ ab_i, const int* __restrict__ ab_j,
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const float* __restrict__ fx, const float* __restrict__ fy, const float* __restrict__ fz, int N,
    const int* __restrict__ Z, const double* __restrict__ bas, const double* __restrict__ aci,
    const double* __restrict__ q,
    const int* __restrict__ hyd_list, int n_hyd,
    const int* __restrict__ bond_off, const int* __restrict__ bond_idx,
    const int* __restrict__ hyd_off, const int* __restrict__ hyd_idx,
    const int* __restrict__ nbr_off, const int* __restrict__ nbr_idx,
    double thr1, double thr2, double min_pair_energy,
    int cap, int* d_counter, long long* out_key, int* out_ct, int* out_par)
{
    hb_detect_body<false>(n_ab, ab_i, ab_j, cx, cy, cz, fx, fy, fz, N, Z, bas, aci, q, hyd_list, n_hyd, bond_off, bond_idx,
                          hyd_off, hyd_idx, nbr_off, nbr_idx, thr1, thr2, min_pair_energy, cap, d_counter, out_key, out_ct, out_par,
                          nullptr, nullptr);
}

/// L1 (Nov 2026, Claude Generated): k_hb_detect with the float minimum-image pre-filter (periodic cell); same list.
__global__ GFNFF_KERNEL_BOUNDS void k_hb_detect_pbc(
    int n_ab, const int* __restrict__ ab_i, const int* __restrict__ ab_j,
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const float* __restrict__ fx, const float* __restrict__ fy, const float* __restrict__ fz, int N,
    const int* __restrict__ Z, const double* __restrict__ bas, const double* __restrict__ aci,
    const double* __restrict__ q,
    const int* __restrict__ hyd_list, int n_hyd,
    const int* __restrict__ bond_off, const int* __restrict__ bond_idx,
    const int* __restrict__ hyd_off, const int* __restrict__ hyd_idx,
    const int* __restrict__ nbr_off, const int* __restrict__ nbr_idx,
    double thr1, double thr2, double min_pair_energy,
    int cap, int* d_counter, long long* out_key, int* out_ct, int* out_par,
    const int* __restrict__ hnb_off, const int* __restrict__ hnb_idx)
{
    hb_detect_body<true>(n_ab, ab_i, ab_j, cx, cy, cz, fx, fy, fz, N, Z, bas, aci, q, hyd_list, n_hyd, bond_off, bond_idx,
                         hyd_off, hyd_idx, nbr_off, nbr_idx, thr1, thr2, min_pair_energy, cap, d_counter, out_key, out_ct, out_par,
                         hnb_off, hnb_idx);
}

/// L1 (Nov 2026, Claude Generated): hydrogens near each AB atom for k_hb_detect_pbc. Atom a (ab_ok) keeps every HB
/// hydrogen whose float minimum-image r^2 is not above thr2 (1 + 1e-4) (NaN keeps it), in hyd_list order: a superset
/// of the hydrogens with double r_aH^2 <= thr2 (float error <= 5e-3 Bohr^2, see hb_detect_body). Without the float
/// cell (non-orthorhombic, L > 400 Bohr) every hydrogen is kept. count = nullptr: write pass at off.
__global__ void k_hb_hnb(int N, const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const char* __restrict__ ab_ok, const int* __restrict__ hyd_list, int n_hyd, double thr2,
    int* __restrict__ count, const int* __restrict__ off, int* __restrict__ out)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= N) return;
    if (!ab_ok[a]) { if (count) count[a] = 0; return; }
    const DispRowCell fc = dispRowCell();
    const float thr_hi = static_cast<float>(thr2) * 1.0001f;
    const double xa = cx[a], ya = cy[a], za = cz[a];
    int c = 0;
    const int base = count ? 0 : off[a];
    for (int hq = 0; hq < n_hyd; ++hq) {
        const int H = hyd_list[hq];
        bool keep = true;
        if (fc.ok) {
            const float dx = hbWrapF(static_cast<float>(cx[H] - xa), fc.Lx, fc.iLx);
            const float dy = hbWrapF(static_cast<float>(cy[H] - ya), fc.Ly, fc.iLy);
            const float dz = hbWrapF(static_cast<float>(cz[H] - za), fc.Lz, fc.iLz);
            keep = !(dx * dx + dy * dy + dz * dz > thr_hi);
        }
        if (!keep) continue;
        if (!count) out[base + c] = H;
        ++c;
    }
    if (count) count[a] = c;
}

__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_hb_fill_columns(
    int n, int N, const long long* __restrict__ key_sorted, const int* __restrict__ perm,
    const int* __restrict__ in_ct, const int* __restrict__ in_par, const int* __restrict__ Z,
    const double* __restrict__ bas, const double* __restrict__ aci, const double* __restrict__ q,
    const double* __restrict__ repz, const int* __restrict__ nbr_off, const int* __restrict__ nbr_idx,
    int* idx_i, int* idx_j, int* idx_k, int* elem_A, int* elem_B,
    double* q_H, double* q_A, double* q_B, double* bas_A, double* bas_B, double* aci_A, double* aci_B,
    double* r_cut, int* case_type, int* acc_parent, double* repz_B, int* nbB_count, int* nbC_count)
{
    const int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    const long long key = key_sorted[t];
    const long long NN = N;
    const int B = static_cast<int>(key % NN);
    const int H = static_cast<int>((key / NN) % NN);
    const int A = static_cast<int>(key / (NN * NN));
    const int src = perm[t];
    const int ct = in_ct[src];
    const int par = in_par[src];
    idx_i[t] = A; idx_j[t] = H; idx_k[t] = B;
    elem_A[t] = Z[A]; elem_B[t] = Z[B];
    q_H[t] = q[H]; q_A[t] = q[A]; q_B[t] = q[B];
    bas_A[t] = bas[A]; bas_B[t] = bas[B];
    aci_A[t] = aci[A]; aci_B[t] = aci[B];
    r_cut[t] = 50.0;
    case_type[t] = ct;
    acc_parent[t] = par;
    const int zB = Z[B];
    repz_B[t] = (zB >= 1 && zB <= 86) ? repz[zB - 1] : 1.0;
    nbB_count[t] = (ct >= 2) ? (nbr_off[B + 1] - nbr_off[B]) : 0;
    int cC = 0;
    if (ct == 3) {
        for (int p = nbr_off[par]; p < nbr_off[par + 1]; ++p)
            if (nbr_idx[p] != B) ++cC;
    }
    nbC_count[t] = cC;
}

__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_hb_fill_neighbors(
    int n, const int* __restrict__ idx_k, const int* __restrict__ case_type, const int* __restrict__ acc_parent,
    const int* __restrict__ nbr_off, const int* __restrict__ nbr_idx,
    const int* __restrict__ nbB_offset, const int* __restrict__ nbC_offset, int* nbB_flat, int* nbC_flat)
{
    const int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    const int ct = case_type[t];
    if (ct < 2) return;
    const int B = idx_k[t];
    int w = nbB_offset[t];
    for (int p = nbr_off[B]; p < nbr_off[B + 1]; ++p) nbB_flat[w++] = nbr_idx[p];
    if (ct == 3) {
        const int par = acc_parent[t];
        w = nbC_offset[t];
        for (int p = nbr_off[par]; p < nbr_off[par + 1]; ++p)
            if (nbr_idx[p] != B) nbC_flat[w++] = nbr_idx[p];
    }
}

// ============================================================================
// Homogeneous long-range D4 tail on the device (WP-PERIODIC-NONBONDED A3c, S1) - Claude Generated (Nov 2026)
// Mirror of D4ParameterGenerator::homogeneousDispersionTail(): beyond the evaluation radius r_c every atom j
// has images with density 1/V, so
//   E_tail = pref6 * S6 + pref8 * S8,  pref6 = -(2 pi/V)/(3 r_c^3),  pref8 = -(2 pi/V) 2/(5 r_c^5),
//   S6 = sum_ij zeta_i zeta_j C6_ij,  S8 = sum_ij zeta_i zeta_j (3 sqz_i sqz_j) C6_ij   (all ordered pairs, i = j too),
// factorised per partner element slot s and reference b: W[s][b] = sum_{j in s} zeta_j w_j[b] (k_d4_tail_W, one
// thread per (s, b), atoms in index order: deterministic), then per atom i (k_d4_tail_atoms, same loop order as
// the CPU): c6sum_i = sum_s sum_a w_i[a] sum_b C6ref(Z_i,a; Z_s,b) W[s][b], dE/dCN_i = 2 zeta_i (pref6 d6 + pref8
// r4 d8) with the weight derivatives in place of w_i[a]. The energy is summed per atom in fixed point
// (blockReduceAddEnergy, slot E_DISP_TAIL), dE/dCN by fix_add into d_dEdcn before the CN chain rule. Differs
// from the CPU by the summation order of S6/S8 only (per-atom partials instead of two running sums).
// ============================================================================
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_d4_tail_W(
    int N, int K, const int* __restrict__ slot_of_atom, const int* __restrict__ nref_slot,
    const double* __restrict__ zeta, const double* __restrict__ gw, double* __restrict__ W)
{
    const int t = blockIdx.x * blockDim.x + threadIdx.x;   // t = s * D4_MAX_REF + b
    if (t >= K * D4_MAX_REF) return;
    const int s = t / D4_MAX_REF, b = t % D4_MAX_REF;
    double acc = 0.0;
    if (b < nref_slot[s])
        for (int j = 0; j < N; ++j)
            if (slot_of_atom[j] == s) acc += zeta[j] * gw[j * D4_MAX_REF + b];
    W[t] = acc;
}

__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_d4_tail_atoms(
    int N, int K, const int* __restrict__ slot_of_atom, const int* __restrict__ slot_elem,
    const int* __restrict__ nref_slot, const double* __restrict__ sqz_slot, const double* __restrict__ zeta,
    const double* __restrict__ gw, const double* __restrict__ dgw, const double* __restrict__ c6_flat,
    const double* __restrict__ W, double pref6, double pref8, int want_grad,
    double* __restrict__ dEdcn, double* __restrict__ energy)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;
    if (i < N) {
        const int si = slot_of_atom[i];
        const int ei = slot_elem[si] - 1;
        const int na = nref_slot[si];
        double c6sum = 0.0, c8sum = 0.0, d6sum = 0.0, d8sum = 0.0;
        for (int s = 0; s < K; ++s) {
            const size_t base = static_cast<size_t>(ei) * D4_MAX_ELEM * D4_MAX_REF * D4_MAX_REF
                              + static_cast<size_t>(slot_elem[s] - 1) * D4_MAX_REF * D4_MAX_REF;
            const double* Ws = W + s * D4_MAX_REF;
            const int nb = nref_slot[s];
            double c = 0.0, d = 0.0;
            for (int a = 0; a < na; ++a) {
                double h = 0.0;
                const double* row = c6_flat + base + static_cast<size_t>(a) * D4_MAX_REF;
                for (int b = 0; b < nb; ++b) h += row[b] * Ws[b];
                c += gw[i * D4_MAX_REF + a] * h;
                if (want_grad) d += dgw[i * D4_MAX_REF + a] * h;
            }
            c6sum += c;
            c8sum += sqz_slot[s] * c;
            d6sum += d;
            d8sum += sqz_slot[s] * d;
        }
        const double zi = zeta[i];
        const double r4 = 3.0 * sqz_slot[si];
        local_E = pref6 * (zi * c6sum) + pref8 * (zi * r4 * c8sum);
        if (want_grad)
            fix_add(&dEdcn[i], 2.0 * zi * (pref6 * d6sum + pref8 * r4 * d8sum), GFNFF_FIX_FORCE);
    }
    blockReduceAddEnergy(local_E, energy);
}

// ============================================================================
// Virial-only pair kernels (Claude Generated, Nov 2026, GPU_PERIODIC_GAPS A1)
// One thread per pair, the pair vector and the factor fac = (dE/dr)/r exactly as the energy/gradient kernel forms
// them (k_bonds, k_dispersion, k_repulsion, k_cn_chainrule, k_hb_alpha_chainrule; those kernels are not touched),
// Xi -= fac d d^T into the nine slots `vir`. Launched after the step's gradient kernels on the main stream, only
// when the virial is requested; they read the same inputs (CN, dc6dcn, dEdcn after the conversion, zz_hb).
// ============================================================================
__global__ __launch_bounds__(512) void k_bonds_virial(
    int n, const int* __restrict__ idx_i, const int* __restrict__ idx_j, const double* __restrict__ r0,
    const double* __restrict__ r0_base_i, const double* __restrict__ r0_base_j, const double* __restrict__ cnfak_i,
    const double* __restrict__ cnfak_j, const double* __restrict__ rabshift, const double* __restrict__ ff,
    const double* __restrict__ fc, const double* __restrict__ alpha, const int* __restrict__ nr_hb,
    const double* __restrict__ hb_cn_H, const double* __restrict__ cx, const double* __restrict__ cy,
    const double* __restrict__ cz, const double* __restrict__ cn, double* vir)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double v[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    if (tid < n) {
        const int i = idx_i[tid], j = idx_j[tid];
        double dx = cx[i] - cx[j], dy = cy[i] - cy[j], dz = cz[i] - cz[j];
        applyMIC(dx, dy, dz);
        const double rij = sqrt(dx * dx + dy * dy + dz * dz);
        if (rij >= 1e-10) {
            double r0_ij;
            if (cn && r0_base_i[tid] > 0.0) {
                const double ra = r0_base_i[tid] + cnfak_i[tid] * cn[i];
                const double rb = r0_base_j[tid] + cnfak_j[tid] * cn[j];
                r0_ij = (ra + rb + rabshift[tid]) * ff[tid];
            } else {
                r0_ij = r0[tid];
            }
            const double dr = rij - r0_ij;
            double alp = alpha[tid];
            if (nr_hb[tid] >= 1)
                alp = (-(1.0 - 0.9) * hb_cn_H[tid] + 1.0) * alpha[tid];
            const double E = fc[tid] * exp(-alp * dr * dr);
            const double dEdr = -2.0 * alp * dr * E;
            vir_pair(v, dEdr / rij, dx, dy, dz);
        }
    }
    blockReduceAddVirial(v, vir);
}

__global__ __launch_bounds__(512) void k_dispersion_virial(
    int n, const int* __restrict__ idx_i, const int* __restrict__ idx_j, const double* __restrict__ C6,
    const double* __restrict__ r4r2ij, const double* __restrict__ r0_sq, const double* __restrict__ zetac6,
    const double* __restrict__ r_cut, const double* __restrict__ cx, const double* __restrict__ cy,
    const double* __restrict__ cz, double* vir)
{
    double v[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    const int stride = gridDim.x * blockDim.x;
    for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < n; tid += stride) {
        const int i = idx_i[tid], j = idx_j[tid];
        double dx = cx[i] - cx[j], dy = cy[i] - cy[j], dz = cz[i] - cz[j];
        applyMIC(dx, dy, dz);
        const double r2 = dx * dx + dy * dy + dz * dz;
        const double rij = sqrt(r2);
        if (rij <= r_cut[tid] && rij >= 1e-10) {
            const double r4 = r2 * r2, r6 = r4 * r2, r0s = r0_sq[tid];
            const double r0_6 = r0s * r0s * r0s;
            const double t6 = 1.0 / (r6 + r0_6);
            const double t8 = 1.0 / (r6 * r2 + r0_6 * r0s);
            const double d6 = -6.0 * r4 * t6 * t6;
            const double d8 = -8.0 * r4 * r2 * t8 * t8;
            const double dEdr = -C6[tid] * zetac6[tid] * (d6 + 2.0 * r4r2ij[tid] * d8) * rij;
            vir_pair(v, dEdr / rij, dx, dy, dz);
        }
    }
    blockReduceAddVirial(v, vir);
}

__global__ __launch_bounds__(512) void k_repulsion_virial(
    int n, const int* __restrict__ idx_i, const int* __restrict__ idx_j, const double* __restrict__ alpha,
    const double* __restrict__ repab, const double* __restrict__ r_cut, const double* __restrict__ cx,
    const double* __restrict__ cy, const double* __restrict__ cz, double* vir)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double v[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    if (tid < n) {
        const int i = idx_i[tid], j = idx_j[tid];
        double dx = cx[i] - cx[j], dy = cy[i] - cy[j], dz = cz[i] - cz[j];
        applyMIC(dx, dy, dz);
        const double rij = sqrt(dx * dx + dy * dy + dz * dz);
        if (rij <= r_cut[tid] && rij >= 1e-8) {
            const double alp_r = alpha[tid] * (rij * sqrt(rij));
            if (!isnan(alp_r) && !isnan(repab[tid]) && alp_r <= 700.0) {
                const double base_E = repab[tid] * exp(-alp_r) / rij;
                const double dEdr = -base_E / rij - 1.5 * alpha[tid] * sqrt(rij) * base_E;
                vir_pair(v, dEdr / rij, dx, dy, dz);
            }
        }
    }
    blockReduceAddVirial(v, vir);
}

/// CN chain rule (k_cn_chainrule): fac = dS/dr / r (dEdcn_i dlogdcn_i + dEdcn_j dlogdcn_j), all CN-dependent terms
__global__ __launch_bounds__(512) void k_cn_chainrule_virial(
    int n_pairs, const int* __restrict__ idx_i, const int* __restrict__ idx_j, const double* __restrict__ rcov_sum,
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const double* __restrict__ dEdcn, const double* __restrict__ dlogdcn, double kn, double r2_max, double* vir)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double v[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    if (tid < n_pairs) {
        const int i = idx_i[tid], j = idx_j[tid];
        double dx = cx[i] - cx[j], dy = cy[i] - cy[j], dz = cz[i] - cz[j];
        applyMIC(dx, dy, dz);
        const double r2 = dx * dx + dy * dy + dz * dz;
        const double rij = sqrt(r2);
        if (rij >= 1e-10 && r2 <= r2_max) {   // A16: the CN range of k_cn_chainrule
            const double rcov = rcov_sum[tid];
            const double dr = (rij - rcov) / rcov;
            const double dSdr = kn * 0.5641895835477563 * exp(-kn * kn * dr * dr) / rcov;
            vir_pair(v, dSdr / rij * (dEdcn[i] * dlogdcn[i] + dEdcn[j] * dlogdcn[j]), dx, dy, dz);
        }
    }
    blockReduceAddVirial(v, vir);
}

/// HB-CN modulation of the X-H bonds (k_hb_alpha_chainrule): fac = zz_H dS/dr / r; counted with the bond term (CPU)
__global__ __launch_bounds__(512) void k_hb_alpha_virial(
    int n_pairs, const int* __restrict__ idx_H, const int* __restrict__ idx_B, const double* __restrict__ rcov_sum,
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const double* __restrict__ zz_hb, double kn, double* vir)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double v[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    if (tid < n_pairs) {
        const int H = idx_H[tid], B = idx_B[tid];
        double dx = cx[H] - cx[B], dy = cy[H] - cy[B], dz = cz[H] - cz[B];
        applyMIC(dx, dy, dz);
        const double rij = sqrt(dx * dx + dy * dy + dz * dz);
        if (rij >= 1e-10) {
            const double rcov = rcov_sum[tid];
            const double dr = (rij - rcov) / rcov;
            const double dSdr = -kn * 0.5641895835477563 * exp(-kn * kn * dr * dr) / rcov;
            vir_pair(v, zz_hb[H] * dSdr / rij, dx, dy, dz);
        }
    }
    blockReduceAddVirial(v, vir);
}

// Claude Generated (Oct 2026): host access to d_gfnff_fix_bad, which lives in this translation unit only.
void gfnffFixBadReadAsync(int* h_pinned, cudaStream_t stream)
{
    cudaMemcpyFromSymbolAsync(h_pinned, d_gfnff_fix_bad, sizeof(int), 0, cudaMemcpyDeviceToHost, stream);
}

void gfnffFixBadReset()
{
    const int zero = 0;
    cudaMemcpyToSymbol(d_gfnff_fix_bad, &zero, sizeof(int));
}



// ============================================================================
// Claude Generated (Nov 2026, S1): host launch wrappers of the PBC-templated local-term kernels. The templates are
// instantiated here (nvcc gives the host stubs of template kernels hidden visibility, so another translation unit
// cannot launch them); pbc false launches the former kernel, true the image-frame variant.
// ============================================================================
void launch_k_angles(bool pbc, dim3 grid, dim3 block, cudaStream_t stream, int n, const int* idx_i, const int* idx_j, const int* idx_k, const int* ati, const int* atj, const int* atk, const double* fc, const double* theta0, const double* cx, const double* cy, const double* cz, const double* rcov_d3, double* grad, double* energy, double* vir)
{
    // A1 (Nov 2026): vir != nullptr (only with a cell) launches the virial instantiation
    if (pbc && vir) k_angles_vir<<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, ati, atj, atk, fc, theta0, cx, cy, cz, rcov_d3, grad, energy, vir);
    else if (pbc) k_angles<true><<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, ati, atj, atk, fc, theta0, cx, cy, cz, rcov_d3, grad, energy);
    else     k_angles<false><<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, ati, atj, atk, fc, theta0, cx, cy, cz, rcov_d3, grad, energy);
}

void launch_k_dihedrals(bool pbc, dim3 grid, dim3 block, cudaStream_t stream, int n, const int* idx_i, const int* idx_j, const int* idx_k, const int* idx_l, const int* ati, const int* atj, const int* atk, const int* atl, const double* V, const double* phi0, const double* n_period, const int* is_nci, const double* cx, const double* cy, const double* cz, const double* rcov_d3_unused, double* grad, double* energy, double* vir)
{
    // A1 (Nov 2026): vir != nullptr (only with a cell) launches the virial instantiation
    if (pbc && vir) k_dihedrals_vir<<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, idx_l, ati, atj, atk, atl, V, phi0, n_period, is_nci, cx, cy, cz, rcov_d3_unused, grad, energy, vir);
    else if (pbc) k_dihedrals<true><<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, idx_l, ati, atj, atk, atl, V, phi0, n_period, is_nci, cx, cy, cz, rcov_d3_unused, grad, energy);
    else     k_dihedrals<false><<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, idx_l, ati, atj, atk, atl, V, phi0, n_period, is_nci, cx, cy, cz, rcov_d3_unused, grad, energy);
}

void launch_k_inversions(bool pbc, dim3 grid, dim3 block, cudaStream_t stream, int n, const int* idx_i, const int* idx_j, const int* idx_k, const int* idx_l, const int* ati, const int* atj, const int* atk, const int* atl, const double* fc, const double* omega0, const double* C0, const double* C1, const double* C2, const int* potential_type, const double* cx, const double* cy, const double* cz, double* grad, double* energy, double* vir)
{
    // A1 (Nov 2026): vir != nullptr (only with a cell) launches the virial instantiation
    if (pbc && vir) k_inversions_vir<<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, idx_l, ati, atj, atk, atl, fc, omega0, C0, C1, C2, potential_type, cx, cy, cz, grad, energy, vir);
    else if (pbc) k_inversions<true><<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, idx_l, ati, atj, atk, atl, fc, omega0, C0, C1, C2, potential_type, cx, cy, cz, grad, energy);
    else     k_inversions<false><<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, idx_l, ati, atj, atk, atl, fc, omega0, C0, C1, C2, potential_type, cx, cy, cz, grad, energy);
}

void launch_k_storsions(bool pbc, dim3 grid, dim3 block, cudaStream_t stream, int n, const int* idx_i, const int* idx_j, const int* idx_k, const int* idx_l, const double* erefhalf, const double* cx, const double* cy, const double* cz, double* grad, double* energy, double* vir)
{
    // A1 (Nov 2026): vir != nullptr (only with a cell) launches the virial instantiation
    if (pbc && vir) k_storsions_vir<<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, idx_l, erefhalf, cx, cy, cz, grad, energy, vir);
    else if (pbc) k_storsions<true><<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, idx_l, erefhalf, cx, cy, cz, grad, energy);
    else     k_storsions<false><<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, idx_l, erefhalf, cx, cy, cz, grad, energy);
}

void launch_k_xbonds(bool pbc, dim3 grid, dim3 block, cudaStream_t stream, int n, const int* idx_i, const int* idx_j, const int* idx_k, const int* elem_A, const int* elem_B, const double* q_X, const double* q_B, const double* acidity_X, const double* r_cut, const double* cx, const double* cy, const double* cz, double* grad, double* energy, double* vir)
{
    // A1 (Nov 2026): vir != nullptr (only with a cell) launches the virial instantiation
    if (pbc && vir) k_xbonds_vir<<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, elem_A, elem_B, q_X, q_B, acidity_X, r_cut, cx, cy, cz, grad, energy, vir);
    else if (pbc) k_xbonds<true><<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, elem_A, elem_B, q_X, q_B, acidity_X, r_cut, cx, cy, cz, grad, energy);
    else     k_xbonds<false><<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, elem_A, elem_B, q_X, q_B, acidity_X, r_cut, cx, cy, cz, grad, energy);
}

void launch_k_hbonds(bool pbc, dim3 grid, dim3 block, cudaStream_t stream, int n, const int* idx_i, const int* idx_j, const int* idx_k, const int* elem_A, const int* elem_B, const double* q_H, const double* q_A, const double* q_B, const double* basicity_A, const double* basicity_B, const double* acidity_A, const double* acidity_B, const double* r_cut, const int* case_type, const int* nb_B_offset, const int* nb_B_count, const int* nb_B_flat, const int* acceptor_parent, const int* nb_C_offset, const int* nb_C_count, const int* nb_C_flat, const double* repz_B, const double* cx, const double* cy, const double* cz, double* grad, double* energy, double* e_triple, double* vir)
{
    // A1 (Nov 2026): vir != nullptr (only with a cell) launches the virial instantiation
    if (pbc && vir) k_hbonds_vir<<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, elem_A, elem_B, q_H, q_A, q_B, basicity_A, basicity_B, acidity_A, acidity_B, r_cut, case_type, nb_B_offset, nb_B_count, nb_B_flat, acceptor_parent, nb_C_offset, nb_C_count, nb_C_flat, repz_B, cx, cy, cz, grad, energy, e_triple, vir);
    else if (pbc) k_hbonds<true><<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, elem_A, elem_B, q_H, q_A, q_B, basicity_A, basicity_B, acidity_A, acidity_B, r_cut, case_type, nb_B_offset, nb_B_count, nb_B_flat, acceptor_parent, nb_C_offset, nb_C_count, nb_C_flat, repz_B, cx, cy, cz, grad, energy, e_triple);
    else     k_hbonds<false><<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, elem_A, elem_B, q_H, q_A, q_B, basicity_A, basicity_B, acidity_A, acidity_B, r_cut, case_type, nb_B_offset, nb_B_count, nb_B_flat, acceptor_parent, nb_C_offset, nb_C_count, nb_C_flat, repz_B, cx, cy, cz, grad, energy, e_triple);
}

// A1 (Nov 2026): BATM launch (template on VIR, host stub visibility as above)
void launch_k_batm(double* vir, dim3 grid, dim3 block, cudaStream_t stream, int n, const int* idx_i, const int* idx_j, const int* idx_k, const double* zb3atm_i, const double* zb3atm_j, const double* zb3atm_k, const double* cx, const double* cy, const double* cz, const double* topo_charges, double* grad, double* energy)
{
    if (vir) k_batm_vir<<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, zb3atm_i, zb3atm_j, zb3atm_k, cx, cy, cz, topo_charges, grad, energy, vir);
    else     k_batm<<<grid, block, 0, stream>>>(n, idx_i, idx_j, idx_k, zb3atm_i, zb3atm_j, zb3atm_k, cx, cy, cz, topo_charges, grad, energy);
}
