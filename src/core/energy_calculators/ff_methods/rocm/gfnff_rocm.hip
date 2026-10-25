/* Claude Generated (June 2026): ROCm/HIP GFN-FF device TU (single-TU port of
 * cuda/{gfnff_kernels,ff_workspace_gpu,gpu_utils}.cu+cpp). hipcc -> one .o, linked into
 * curcuma_core with g++ (EXTERNAL_OBJECT). M1: energy+gradient; EEQ (rocSOLVER) in M2. */
#ifdef USE_ROCM
#include <hip/hip_runtime.h>

// Parse the host-facing classes BEFORE gfnff_kernels_hip.h's "#define D4_MAX_REF 7"
// (which would otherwise clobber FFWorkspaceHip's `static constexpr int D4_MAX_REF`
// member — the CUDA build includes ff_workspace_gpu.h before the kernel header too).
#include "ff_workspace_hip.h"
#include "eeq_solver_hip.h"

#include <csignal>
// Claude Generated (June 2026): GPU-wedge hardening. SIGKILL (kill -9) can't be caught, but
// SIGINT (Ctrl-C), SIGTERM (kill / `timeout`) and SIGABRT CAN — and if the process dies while a
// HIP kernel is in flight without tearing down the device, the RDNA compute queue can wedge
// (every later submission then hangs at init until a GPU reset/reboot). We install a one-shot
// handler that best-effort resets the device, restores the default disposition, and re-raises,
// so an interrupted run cleans up the GPU instead of leaving it hung. (hipDeviceReset from a
// handler is not strictly async-signal-safe, but a clean GPU on abort beats a wedged one; the
// only true non-recoverable case is kill -9, which the operator must simply avoid.)
namespace {
volatile std::sig_atomic_t g_hip_sig_installed = 0;
void hipCleanupSignalHandler(int sig) {
    hipDeviceReset();
    std::signal(sig, SIG_DFL);
    std::raise(sig);
}
void installHipSignalHandlers() {
    if (g_hip_sig_installed) return;
    g_hip_sig_installed = 1;
    std::signal(SIGINT,  hipCleanupSignalHandler);
    std::signal(SIGTERM, hipCleanupSignalHandler);
    std::signal(SIGABRT, hipCleanupSignalHandler);
}
} // namespace

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

// ============================================================================
// Order-independent (fixed-point) accumulation, WP10 (Claude Generated, Oct 2026; port of cuda/gfnff_kernels.cu
// fix_add, docs/GPU_REPRODUCIBILITY.md). FP64 atomicAdd depends on the arrival order of the threads, so two
// runs differ in the last bits. Integer addition is associative: with d_gfnff_fix != 0 a contribution is
// rounded to a multiple of 2^-52 (forces, dE/dCN, zz_hb; energy slot 0 = dispersion), 2^-48 (other energy
// slots) or 2^-56 (CN sums), added as a 64-bit integer into the same 8 bytes of the double buffer, and the
// buffer is converted in place by k_fix_to_double once the accumulation is over. Same scales as CUDA.
// d_gfnff_fix == 0 keeps the FP64 atomicAdd path unchanged (switch -gfnff.rocm_fixed_point).
// Range (as CUDA): +-2048 per gradient/dE/dCN component and for the dispersion energy, +-32768 Eh per other
// energy slot, +-128 for CN sums; a final sum beyond that wraps silently.
// ============================================================================
constexpr double GFNFF_FIX_CN          = 72057594037927936.0;   // 2^56
constexpr double GFNFF_FIX_ENERGY      = 281474976710656.0;     // 2^48
constexpr double GFNFF_FIX_ENERGY_DISP = 4503599627370496.0;    // 2^52
constexpr double GFNFF_FIX_FORCE       = 4503599627370496.0;    // 2^52

/// 1: accumulators of the GFN-FF step are fixed-point integers (set by setGfnffFixModeHip, host side).
__constant__ int d_gfnff_fix;

/// Set (sticky) when a contribution was not finite or alone outside the accumulator range; read by the host after
/// the step (FFWorkspaceHip::launchChargeDependentAndFinish), which then returns NaN as the FP64 path would.
__device__ int d_gfnff_fix_bad;

/// The CUDA version adds |v| >= 64 in chunks of 64 with one atomic each (unbounded for a huge or infinite v).
/// Here the n chunks are folded into one integer, n * llrint(64 scale) + llrint(r scale), r = v - n * 64 (exact for
/// |v| < 2^59): modulo 2^64 the same accumulator bits as the chunk loop, with one atomic and a bounded cost.
__device__ __forceinline__ void fix_add_u64(unsigned long long* a, double v, double scale)
{
    const double av = fabs(v);
    if (!(av * scale < 9.2e18)) {                     // beyond int64 range (or NaN/Inf): the sum is wrong anyway
        if (!(av < 5.7e17)) {                         // 2^59 > |v| or not finite: mark, add nothing usable
            atomicOr(&d_gfnff_fix_bad, 1);
            return;
        }
        atomicOr(&d_gfnff_fix_bad, 1);
    }
    if (av < 64.0) {                                  // the common case: one rounded value
        atomicAdd(a, static_cast<unsigned long long>(llrint(v * scale)));
        return;
    }
    const double n = trunc(av * (1.0 / 64.0));        // number of 64-chunks of the CUDA loop
    const double r = v - copysign(n * 64.0, v);       // exact: the loop's remainder
    const unsigned long long chunk = static_cast<unsigned long long>(llrint(copysign(64.0, v) * scale));
    const unsigned long long sum = static_cast<unsigned long long>(n) * chunk
                                 + static_cast<unsigned long long>(llrint(r * scale));
    atomicAdd(a, sum);
}

__device__ __forceinline__ void fix_add(double* addr, double v, double scale)
{
    fix_add_u64(reinterpret_cast<unsigned long long*>(addr), v, scale);
}

/// Claude Generated (Oct 2026, WP10): one accumulation into a step accumulator: fixed point with the given
/// scale when d_gfnff_fix is set, otherwise the FP64 atomicAdd of the original kernels.
__device__ __forceinline__ void acc_add(double* addr, double v, double scale)
{
    if (d_gfnff_fix) fix_add(addr, v, scale);
    else             atomicAdd(addr, v);
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

/// Add vector contribution to gradient (FP64 atomicAdd on x,y,z, or fixed point with d_gfnff_fix)
__device__ __forceinline__ void add_grad(double* __restrict__ grad, int atom,
                                          double fx, double fy, double fz)
{
    if (d_gfnff_fix) {
        fix_add(&grad[3*atom],   fx, GFNFF_FIX_FORCE);
        fix_add(&grad[3*atom+1], fy, GFNFF_FIX_FORCE);
        fix_add(&grad[3*atom+2], fz, GFNFF_FIX_FORCE);
        return;
    }
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

/// Claude Generated (Oct 2026, WP10): integer twin of blockReduceAddEnergy (port of CUDA blockReduceAddEnergyFix).
/// Thread values are fixed-point integers, so warp, block and atomic sums are associative. Wave32 as below.
__device__ __forceinline__ void blockReduceAddEnergyFix(long long local_E, double* energy)
{
    long long warp_sum = local_E;
    for (int off = 16; off > 0; off >>= 1) warp_sum += __shfl_down_sync(0xFFFFFFFFFFFFFFFFull, warp_sum, off);

    __shared__ long long warp_sums_fx[32];
    int lane = threadIdx.x & 31;
    int warp_id = threadIdx.x >> 5;

    if (lane == 0)
        warp_sums_fx[warp_id] = warp_sum;
    __syncthreads();

    int num_warps = (blockDim.x + 31) >> 5;
    if (warp_id == 0) {
        long long sum = (lane < num_warps) ? warp_sums_fx[lane] : 0ll;
        for (int off = 16; off > 0; off >>= 1) sum += __shfl_down_sync(0xFFFFFFFFFFFFFFFFull, sum, off);
        if (lane == 0)
            atomicAdd(reinterpret_cast<unsigned long long*>(energy), static_cast<unsigned long long>(sum));
    }
}

__device__ __forceinline__ void blockReduceAddEnergyFP(double local_E, double* energy);

/// Block energy reduction into one slot; with d_gfnff_fix every thread value is rounded to a multiple of
/// 1/fix_scale first (GFNFF_FIX_ENERGY, dispersion kernels GFNFF_FIX_ENERGY_DISP, as CUDA).
__device__ __forceinline__ void blockReduceAddEnergy(double local_E, double* energy, double fix_scale = GFNFF_FIX_ENERGY)
{
    if (d_gfnff_fix) {
        if (!(fabs(local_E) * fix_scale < 9.2e18)) {   // not finite or beyond the int64 range: mark, contribute 0
            atomicOr(&d_gfnff_fix_bad, 1);
            local_E = 0.0;
        }
        blockReduceAddEnergyFix(static_cast<long long>(llrint(local_E * fix_scale)), energy);
    } else {
        blockReduceAddEnergyFP(local_E, energy);
    }
}

__device__ __forceinline__ void blockReduceAddEnergyFP(double local_E, double* energy)
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

/// Claude Generated (Oct 2026): both displacement tests of a step in one single-workgroup launch,
/// without a zeroed flag buffer: flags[0] = some atom moved more than sqrt(thr_sq) from (rx, ry, rz)
/// (topology reference, checkDisplacement), flags[1] the same against (qx, qy, qz) with qthr_sq (CN
/// pair-list reference, cnPairListStale); qx == nullptr skips the second test (flags[1] = 0).
/// Same per-atom test as k_check_displacement. Launch: one workgroup of 1024 threads.
__global__ __launch_bounds__(1024)
void k_check_displacement_pair(
    int N,
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const double* __restrict__ rx, const double* __restrict__ ry, const double* __restrict__ rz,
    double thr_sq,
    const double* __restrict__ qx, const double* __restrict__ qy, const double* __restrict__ qz,
    double qthr_sq,
    int* __restrict__ flags)
{
    int f0 = 0, f1 = 0;
    for (int i = threadIdx.x; i < N; i += blockDim.x) {
        if (rx) {
            const double dx = cx[i] - rx[i];
            const double dy = cy[i] - ry[i];
            const double dz = cz[i] - rz[i];
            if (dx*dx + dy*dy + dz*dz > thr_sq) f0 = 1;
        }
        if (qx) {
            const double dx = cx[i] - qx[i];
            const double dy = cy[i] - qy[i];
            const double dz = cz[i] - qz[i];
            if (dx*dx + dy*dy + dz*dz > qthr_sq) f1 = 1;
        }
    }
    f0 = __syncthreads_or(f0);
    f1 = __syncthreads_or(f1);
    if (threadIdx.x == 0) { flags[0] = f0; flags[1] = f1; }
}

/// Claude Generated (Oct 2026): zero the per-step accumulators of Phase 1 (energy, gradient, dE/dCN,
/// per-term energies) in one launch instead of four hipMemsetAsync blits.
__global__ void k_zero_step_accumulators(double* __restrict__ energy, double* __restrict__ grad, int n3,
                                         double* __restrict__ dEdcn, int n, double* __restrict__ energies, int nslots)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n3) grad[i] = 0.0;
    if (i < n) dEdcn[i] = 0.0;
    if (i < nslots) energies[i] = 0.0;
    if (i == 0) *energy = 0.0;
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
                acc_add(&dEdcn[i], -dc6dcn_ij[tid] * disp_value, GFNFF_FIX_FORCE);
                acc_add(&dEdcn[j], -dc6dcn_ji[tid] * disp_value, GFNFF_FIX_FORCE);
            }
        }
    }

    blockReduceAddEnergy(local_E, energy, GFNFF_FIX_ENERGY_DISP);   // one pair per thread: rounded per pair
}

// Claude Generated (June 2026): per-atom GATHER variant of k_dispersion for ROCm/RDNA. In
// gradient mode the deferred pair k_dispersion is the dominant MD cost (~113 ms on polymer):
// ~6M grad + ~2M dEdcn FP64 atomicAdd over 943k pairs. One thread owns atom a, sweeps a's CSR
// edges, accumulates force (force on a = fac*(ra-rp), fac pair-symmetric), half-energy and the
// dEdcn term locally, then writes grad[a] (3 atomicAdd) + dEdcn[a] (1 atomicAdd). The dynamic
// per-pair dc6dcn stays pair-indexed (read via csr_pairidx; csr_is_i picks ij vs ji) and needs no
// per-step CSR rebuild. Oct 2026: pair arithmetic without sqrt and with one division (equal to
// k_dispersion up to FP rounding); the default gradient path is k_dispersion_gather_uv below.
__global__ void k_dispersion_gather(
    int N,
    const int*    __restrict__ csr_off,
    const int*    __restrict__ csr_partner,
    const int*    __restrict__ csr_pairidx,
    const int*    __restrict__ csr_is_i,
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
    // Claude Generated (Oct 2026): one workgroup per atom, lanes stride the atom's CSR row. The former
    // one-thread-per-atom loop made neighbouring lanes read CSR rows ~N_partner entries apart (no coalescing);
    // 38 ms per step on 7320 atoms (rocprofv3), the largest kernel after the Coulomb far-field branch.
    // Same pair terms; the per-atom sums are added in a different order (FP reordering only).
    const int a = blockIdx.x;
    const int lane = threadIdx.x, nlanes = blockDim.x;
    double local_E = 0.0;
    double gx = 0.0, gy = 0.0, gz = 0.0, dEdcn_a = 0.0;

    if (a < N) {
        const double xa = cx[a], ya = cy[a], za = cz[a];
        const bool has_dcn = (dc6dcn_ij != nullptr && dEdcn != nullptr);
        const int s = csr_off[a], e = csr_off[a + 1];
        for (int t = s + lane; t < e; t += nlanes) {
            const int p   = csr_partner[t];
            const int tid = csr_pairidx[t];
            double dx = xa - cx[p];
            double dy = ya - cy[p];
            double dz = za - cz[p];
            applyMIC(dx, dy, dz);
            double r2  = dx*dx + dy*dy + dz*dz;
            // Claude Generated (Oct 2026): cutoff tested on r^2 (no sqrt), t6 and t8 from one FP64
            // reciprocal, and fac = (dE/dr)/r without the r factor that cancels. Same terms as
            // k_dispersion up to FP rounding. r <= 60 Bohr keeps (r6+R0^6)(r8+R0^8) below ~1e25.
            const double rc = r_cut[tid];
            if (r2 <= rc * rc && r2 >= 1e-20) {
                double r4   = r2 * r2;
                double r6   = r4 * r2;
                double r0s  = r0_sq[tid];
                double r0_6 = r0s * r0s * r0s;
                double r8   = r6 * r2;
                double r0_8 = r0_6 * r0s;
                double den6 = r6 + r0_6;
                double den8 = r8 + r0_8;
                double inv  = 1.0 / (den6 * den8);
                double t6   = den8 * inv;
                double t8   = den6 * inv;
                double rr   = r4r2ij[tid];
                double disp_sum = t6 + 2.0 * rr * t8;
                local_E += 0.5 * (-C6[tid] * disp_sum * zetac6[tid]);
                double d6   = -6.0 * r4 * t6 * t6;
                double d8   = -8.0 * r4 * r2 * t8 * t8;
                double fac  = -C6[tid] * zetac6[tid] * (d6 + 2.0 * rr * d8);   // force on a = fac*(ra - rp)
                gx += fac*dx; gy += fac*dy; gz += fac*dz;
                if (has_dcn) {
                    double disp_value = disp_sum * zetac6[tid];
                    double dcn_a = csr_is_i[t] ? dc6dcn_ij[tid] : dc6dcn_ji[tid];
                    dEdcn_a += -dcn_a * disp_value;
                }
            }
        }
    }

    // Reduce the force components and dEdcn over the workgroup (shared-memory tree), lane 0 writes.
    __shared__ double sgx[256], sgy[256], sgz[256], sdc[256];
    sgx[lane] = gx; sgy[lane] = gy; sgz[lane] = gz; sdc[lane] = dEdcn_a;
    __syncthreads();
    for (int stride = nlanes >> 1; stride > 0; stride >>= 1) {
        if (lane < stride) {
            sgx[lane] += sgx[lane + stride];
            sgy[lane] += sgy[lane + stride];
            sgz[lane] += sgz[lane + stride];
            sdc[lane] += sdc[lane + stride];
        }
        __syncthreads();
    }
    if (lane == 0 && a < N) {
        add_grad(grad, a, sgx[0], sgy[0], sgz[0]);
        if (dc6dcn_ij != nullptr && dEdcn != nullptr) acc_add(&dEdcn[a], sdc[0], GFNFF_FIX_FORCE);
    }

    blockReduceAddEnergy(local_E, energy, GFNFF_FIX_ENERGY_DISP);
}

// Claude Generated (Oct 2026): k_dispersion_gather with C6 and dC6/dCN computed in the edge loop from the
// factorised D4 sums (k_disp_uv) instead of being read from the per-pair arrays C6 / dc6dcn_ij / dc6dcn_ji.
// From atom a's side, for any partner p: C6(a,p) = sum_r v(a,Zp,r) gw(p,r) and dC6/dCN_a = sum_r u(a,Zp,r) gw(p,r).
// For a = idx_j this uses C6ref(Za,Zp,ra,rp) = C6ref(Zp,Za,rp,ra); the host enables this kernel only when the
// uploaded reference table is bitwise symmetric (c6_table_symmetric). Each pair's C6 is now formed twice (once
// per CSR direction) with different rounding; k_dc6dcn_per_pair_uv and its 3 x n per-pair writes are skipped.
// Same energy/force formula as k_dispersion_gather.
/// @param atom_ctype [N] compact element type per atom
/// @param uv         [N * ntypes * 2*D4_MAX_REF] output of k_disp_uv for the current gw
/// @param gw         [N * D4_MAX_REF] Gaussian weights k_disp_uv was built from (zero-padded)
/// @param dEdcn      [N] CN chain-rule accumulator (nullptr: no dC6/dCN term)
__global__ void k_dispersion_gather_uv(
    int N, int ntypes,
    const int*    __restrict__ csr_off,
    const int*    __restrict__ csr_partner,
    const int*    __restrict__ csr_pairidx,
    const int*    __restrict__ atom_ctype,
    const double* __restrict__ uv,
    const double* __restrict__ gw,
    const double* __restrict__ r4r2ij,
    const double* __restrict__ r0_sq,
    const double* __restrict__ zetac6,
    const double* __restrict__ r_cut,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    dEdcn,
    double*                    grad,
    double*                    energy)
{
    const int a = blockIdx.x;
    const int lane = threadIdx.x, nlanes = blockDim.x;
    double local_E = 0.0;
    double gx = 0.0, gy = 0.0, gz = 0.0, dEdcn_a = 0.0;

    if (a < N) {
        const double xa = cx[a], ya = cy[a], za = cz[a];
        const bool has_dcn = (dEdcn != nullptr);
        const double* uv_a = uv + (size_t)a * ntypes * (2 * D4_MAX_REF);
        const int s = csr_off[a], e = csr_off[a + 1];
        for (int t = s + lane; t < e; t += nlanes) {
            const int p   = csr_partner[t];
            const int tid = csr_pairidx[t];
            double dx = xa - cx[p];
            double dy = ya - cy[p];
            double dz = za - cz[p];
            applyMIC(dx, dy, dz);
            double r2  = dx*dx + dy*dy + dz*dz;
            const double rc = r_cut[tid];
            if (r2 <= rc * rc && r2 >= 1e-20) {
                const double* u  = uv_a + atom_ctype[p] * (2 * D4_MAX_REF);
                const double* v  = u + D4_MAX_REF;
                const double* gp = gw + p * D4_MAX_REF;
                double c6 = 0.0, dc6 = 0.0;
#pragma unroll
                for (int r = 0; r < D4_MAX_REF; ++r) {
                    const double g = gp[r];
                    c6  += v[r] * g;
                    dc6 += u[r] * g;
                }
                double r4   = r2 * r2;
                double r6   = r4 * r2;
                double r0s  = r0_sq[tid];
                double r0_6 = r0s * r0s * r0s;
                double r8   = r6 * r2;
                double r0_8 = r0_6 * r0s;
                double den6 = r6 + r0_6;
                double den8 = r8 + r0_8;
                double inv  = 1.0 / (den6 * den8);
                double t6   = den8 * inv;
                double t8   = den6 * inv;
                double rr   = r4r2ij[tid];
                double zc   = zetac6[tid];
                double disp_sum = t6 + 2.0 * rr * t8;
                local_E += 0.5 * (-c6 * disp_sum * zc);
                double d6   = -6.0 * r4 * t6 * t6;
                double d8   = -8.0 * r4 * r2 * t8 * t8;
                double fac  = -c6 * zc * (d6 + 2.0 * rr * d8);   // force on a = fac*(ra - rp)
                gx += fac*dx; gy += fac*dy; gz += fac*dz;
                if (has_dcn) dEdcn_a += -dc6 * (disp_sum * zc);
            }
        }
    }

    __shared__ double sgx[256], sgy[256], sgz[256], sdc[256];
    sgx[lane] = gx; sgy[lane] = gy; sgz[lane] = gz; sdc[lane] = dEdcn_a;
    __syncthreads();
    for (int stride = nlanes >> 1; stride > 0; stride >>= 1) {
        if (lane < stride) {
            sgx[lane] += sgx[lane + stride];
            sgy[lane] += sgy[lane + stride];
            sgz[lane] += sgz[lane + stride];
            sdc[lane] += sdc[lane + stride];
        }
        __syncthreads();
    }
    if (lane == 0 && a < N) {
        add_grad(grad, a, sgx[0], sgy[0], sgz[0]);
        if (dEdcn != nullptr) acc_add(&dEdcn[a], sdc[0], GFNFF_FIX_FORCE);
    }

    blockReduceAddEnergy(local_E, energy, GFNFF_FIX_ENERGY_DISP);
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

// Claude Generated (Oct 2026): per-atom gather variant of k_repulsion for the non-bonded list on ROCm.
// On gfx1200 every FP64 atomicAdd is a compare-and-swap loop (global_atomic_cmpswap_b64), and the pair
// kernel issues 6 of them per pair; the list is grouped by the lower atom index, so the lanes of a
// wavefront hit the same grad[3*i] address and serialise (8.3 ms per step for 1.47M pairs, 6200 atoms).
// Here one workgroup owns atom a and sweeps a's CSR row (each pair appears once per direction), sums the
// force on a and half the pair energy in registers and writes grad[a] with 3 atomicAdd per atom.
// Same pair formula as k_repulsion (E = repab exp(-alpha r^1.5) / r), with one reciprocal of r instead
// of three divisions: equal to k_repulsion up to FP rounding. Force on a = fac * (r_a - r_p).
/// @param N           number of atoms (grid size; one workgroup per atom, blockDim a power of two <= 256)
/// @param csr_off     [N+1] row offsets of the directed-edge CSR
/// @param csr_partner [2n] partner atom per directed edge
/// @param csr_pairidx [2n] index of the pair in the RepulsionSoA arrays
/// @param alpha,repab,r_cut [n] per-pair constants (RepulsionSoA)
/// @param grad        [3N] gradient (3 atomicAdd per atom)
/// @param energy      [1] energy slot (block reduction + one atomicAdd per block)
__global__ void k_repulsion_gather(
    int N,
    const int*    __restrict__ csr_off,
    const int*    __restrict__ csr_partner,
    const int*    __restrict__ csr_pairidx,
    const double* __restrict__ alpha,
    const double* __restrict__ repab,
    const double* __restrict__ r_cut,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy)
{
    const int a = blockIdx.x;
    const int lane = threadIdx.x, nlanes = blockDim.x;
    double local_E = 0.0;
    double gx = 0.0, gy = 0.0, gz = 0.0;

    if (a < N) {
        const double xa = cx[a], ya = cy[a], za = cz[a];
        const int s = csr_off[a], e = csr_off[a + 1];
        for (int t = s + lane; t < e; t += nlanes) {
            const int p = csr_partner[t];
            const int k = csr_pairidx[t];
            double dx = xa - cx[p];
            double dy = ya - cy[p];
            double dz = za - cz[p];
            applyMIC(dx, dy, dz);
            double r2  = dx*dx + dy*dy + dz*dz;
            double rij = sqrt(r2);
            if (rij <= r_cut[k] && rij >= 1e-8) {
                const double sr    = sqrt(rij);
                const double alp   = alpha[k];
                const double rep   = repab[k];
                const double alp_r = alp * (rij * sr);
                if (!isnan(alp_r) && !isnan(rep) && alp_r <= 700.0) {
                    const double inv_r  = 1.0 / rij;
                    const double base_E = rep * exp(-alp_r) * inv_r;
                    local_E += 0.5 * base_E;
                    // dE/dr = -E/r - 1.5 alpha sqrt(r) E; fac = (dE/dr) / r
                    const double fac = -(base_E * inv_r + 1.5 * alp * sr * base_E) * inv_r;
                    gx += fac * dx; gy += fac * dy; gz += fac * dz;
                }
            }
        }
    }

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
    if (lane == 0 && a < N)
        add_grad(grad, a, sgx[0], sgy[0], sgz[0]);

    blockReduceAddEnergy(local_E, energy);
}

// Claude Generated (Oct 2026): per-atom sum of per-slot gradient contributions written by a scatter-free
// term kernel (k_hbonds in slot mode). slot_g holds one 3-vector per slot, atom_slots lists the slots that
// belong to atom a (CSR by atom). One workgroup per atom, shared-memory tree reduction, 3 atomicAdd per
// atom (other streams add to grad concurrently).
/// @param N          number of atoms (grid size; blockDim a power of two <= 256)
/// @param atom_off   [N+1] row offsets into atom_slots
/// @param atom_slots [nslots] slot indices grouped by owning atom
/// @param slot_g     [3*nslots] gradient contribution per slot (x, y, z)
/// @param grad       [3N] gradient (atomicAdd)
__global__ void k_slot_gather(
    int N,
    const int*    __restrict__ atom_off,
    const int*    __restrict__ atom_slots,
    const double* __restrict__ slot_g,
    double*                    grad)
{
    const int a = blockIdx.x;
    const int lane = threadIdx.x, nlanes = blockDim.x;
    double gx = 0.0, gy = 0.0, gz = 0.0;
    int s = 0, e = 0;
    if (a < N) { s = atom_off[a]; e = atom_off[a + 1]; }
    for (int t = s + lane; t < e; t += nlanes) {
        const int sl = atom_slots[t];
        gx += slot_g[3*sl];
        gy += slot_g[3*sl + 1];
        gz += slot_g[3*sl + 2];
    }
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
    if (lane == 0 && a < N && e > s)
        add_grad(grad, a, sgx[0], sgy[0], sgz[0]);
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
                double erf_v, dEdr;
                if (gamma_r >= GFNFF_COUL_FARFIELD_GAMMA_R) {   // far field: erf == 1, Gaussian < 2.3e-16
                    erf_v = 1.0;
                    dEdr  = qi * qj * (-1.0 / (rij * rij));
                } else {
                    erf_v = erf(gamma_r);
                    static const double inv_sqrt_pi = 0.5641895835477563;
                    double exp_v  = exp(-gamma_r * gamma_r);
                    double derf   = gamma_ij[tid] * exp_v * (2.0 * inv_sqrt_pi);
                    dEdr = qi * qj * (derf / rij - erf_v / (rij * rij));
                }
                local_E = qi * qj * erf_v / rij;
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
                    double erf_v, dEdr;
                    if (gamma_r >= GFNFF_COUL_FARFIELD_GAMMA_R) {   // far field, see k_coulomb
                        erf_v = 1.0;
                        dEdr  = qa * qp * (-1.0 / (rij * rij));
                    } else {
                        erf_v = erf(gamma_r);
                        double exp_v = exp(-gamma_r * gamma_r);
                        double derf  = g * exp_v * (2.0 * inv_sqrt_pi);
                        dEdr = qa * qp * (derf / rij - erf_v / (rij * rij));
                    }
                    local_E += 0.5 * qa * qp * erf_v / rij;
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

// Claude Generated (Oct 2026): pair-once tiled Coulomb for ROCm (port of the idea of CUDA k_coulomb_tiles).
//
// E = sum_{i<j} q_i q_j erf(gamma_ij r_ij) / r_ij,   gamma_ij = 1/sqrt(alp_i + alp_j)
// dE/dr = q_i q_j (2 gamma/sqrt(pi) exp(-gamma^2 r^2) / r - erf(gamma r) / r^2)
// k_coulomb_dense evaluated every pair twice (once from each atom's row, 0.5 * E each). Here
// block (bj, bi) with bi <= bj handles tile pair (bi, bj) of GFNFF_COUL_TILE_HIP atoms each, and
// every pair i < j once: thread t owns atom i = bi*T + t, the j tile sits in shared memory, and the
// j gradients are accumulated in shared memory (rotated walk jj = (t + k) mod T so concurrent shared
// atomics rarely hit the same slot) and added to the global gradient once per block.
// Per pair the arithmetic is that of k_coulomb_dense (sqrt and division, no rsqrt; far field
// r^2 >= 36 (alp_i + alp_j), i.e. gamma r >= 6, with erf == 1 exactly; r_cut / 1e-10 test before the
// far-field test; pairs with a NaN charge skipped; applyMIC), so per-pair energies and forces are
// bit-identical to the old kernel and only the summation order changes.
// Launch: grid dim3(nb, nb), nb = ceil(N / T), block T; blocks with bj < bi return at once.
// Formerly k_coulomb_dense (June 2026: coalesced shared-memory N-body sweep over all ordered pairs).
constexpr int GFNFF_COUL_TILE_HIP = 128;   ///< default atoms per tile == block size (power of two)

/// Claude Generated (Oct 2026): below this atom count the per-step kernels run on the main stream only
/// (FFWorkspaceHipImpl::single_stream). Every switch between hardware queues costs the GPU tens of
/// microseconds on the RX 9060 XT, more than small kernels gain from running side by side
/// (MD ms/step, worker streams -> main stream only: 24 atoms 0.77 -> 0.56, 66: 0.96 -> 0.68,
/// 114: 1.17 -> 0.94, 300: 1.81 -> 1.71; but 600: 2.43 -> 2.70, 1002: 3.04 -> 3.24).
constexpr int kRocmSingleStreamMaxAtoms = 500;

/// Claude Generated (Oct 2026): streams of the single-stream workspace whose step is being computed on this
/// host thread (set in prepareAndLaunchChargeIndependent, nullptr otherwise). The small-system device EEQ
/// (eeq_solver_hip.hiph, same translation unit) records g_rocm_small_eeq_event on the main stream and runs
/// on the auxiliary stream (stream_threebody, idle in this mode) behind it. Measured on the RX 9060 XT
/// (MD ms/step, 24 atoms): on the main stream 0.62, on the null stream 0.46 or 0.95 depending on how many
/// streams the process had created (their spread over the four hardware queues), on the auxiliary
/// workspace stream 0.44 with either stream count.
static thread_local hipStream_t g_rocm_small_eeq_stream = nullptr;
static thread_local hipStream_t g_rocm_small_eeq_aux = nullptr;
static thread_local hipEvent_t  g_rocm_small_eeq_event = nullptr;

// Claude Generated (Oct 2026): i tile (block size, TI >= 32 for the wave reduction) and j tile (TJ)
// are separate template parameters. Block (bj, bi) handles i in [bi*TI, bi*TI+TI) and j in
// [bj*TJ, bj*TJ+TJ), every pair with j > i once; blocks whose whole j range lies at or below their
// smallest i return at once. TI == TJ == 128 is the former kernel (same pairs per block, same walk
// order, bit-identical). Small systems take smaller tiles: at 66 atoms the 128 tile is a single
// workgroup on one compute unit and the FP64 erf/exp of all 2145 pairs ran there (~140 us on the
// RX 9060 XT); with 32 x 8 tiles 27 workgroups share the work. Per-pair arithmetic is unchanged;
// only the summation order of energy and gradient differs between tile shapes.
template <int TI, int TJ>
__global__ __launch_bounds__(TI)
void k_coulomb_tiles_hip(
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
    static_assert((TJ & (TJ - 1)) == 0 && TJ <= TI && TI % 32 == 0, "tile shapes: TJ power of two <= TI, TI wave multiple");
    const int bj = blockIdx.x;
    const int bi = blockIdx.y;
    if (bj * TJ + TJ - 1 <= bi * TI) return;   // block-uniform: no pair with j > i in this block

    __shared__ double sx[TJ], sy[TJ], sz[TJ], sq[TJ], sa[TJ];
    __shared__ double sgx[TJ], sgy[TJ], sgz[TJ];
    const int t  = threadIdx.x;
    const int i  = bi * TI + t;

    if (t < TJ) {
        const int jl = bj * TJ + t;
        if (jl < N) { sx[t] = cx[jl]; sy[t] = cy[jl]; sz[t] = cz[jl]; sq[t] = charges[jl]; sa[t] = alpha[jl]; }
        else        { sx[t] = sy[t] = sz[t] = 0.0; sq[t] = 0.0; sa[t] = 1.0; }
        sgx[t] = 0.0; sgy[t] = 0.0; sgz[t] = 0.0;
    }
    __syncthreads();

    const bool active = (i < N) && !isnan(charges[i]);
    double xi = 0, yi = 0, zi = 0, qi = 0, ai = 0;
    if (active) { xi = cx[i]; yi = cy[i]; zi = cz[i]; qi = charges[i]; ai = alpha[i]; }
    const double inv_sqrt_pi = 0.5641895835477563;
    double gx = 0.0, gy = 0.0, gz = 0.0, local_E = 0.0;
    const bool fixmode = (d_gfnff_fix != 0);   // 0.0 and integer 0 have the same bits: the zeroing above serves both
    // WP10 (Claude Generated, Oct 2026): with TI == TJ the rotated walk maps the TI threads one-to-one onto the TJ
    // slots at a fixed k, so with one barrier per k every slot receives its partner forces in k order without
    // atomics (deterministic FP sums; 6200 atoms: cheaper than three fixed-point shared atomics per pair). The
    // per-pair arithmetic is a verbatim copy of the loop below. Smaller j tiles (TJ < TI) keep shared atomics.
    const bool lockstep = fixmode && (TI == TJ);

    if (lockstep) {
        for (int k = 0; k < TJ; ++k) {
            const int jj = (t + k) & (TJ - 1);
            const int j  = bj * TJ + jj;
            double fx = 0.0, fy = 0.0, fz = 0.0;
            bool hit = false;
            if (active && j < N && j > i) {
                double dx = xi - sx[jj];
                double dy = yi - sy[jj];
                double dz = zi - sz[jj];
                applyMIC(dx, dy, dz);
                double r2  = dx*dx + dy*dy + dz*dz;
                double rij = sqrt(r2);
                const double qj = sq[jj];
                if ((rij <= r_cut && rij >= 1e-10) && !isnan(qj)) {
                    const double asum = ai + sa[jj];
                    double erf_v, dEdr;
                    if (r2 >= 36.0 * asum) {
                        erf_v = 1.0;
                        dEdr  = qi * qj * (-1.0 / (rij * rij));
                    } else {
                        double g       = 1.0 / sqrt(asum);
                        double gamma_r = g * rij;
                        erf_v = erf(gamma_r);
                        double exp_v = exp(-gamma_r * gamma_r);
                        double derf  = g * exp_v * (2.0 * inv_sqrt_pi);
                        dEdr = qi * qj * (derf / rij - erf_v / (rij * rij));
                    }
                    local_E += qi * qj * erf_v / rij;
                    const double fac = dEdr / rij;
                    fx = fac*dx; fy = fac*dy; fz = fac*dz;
                    gx += fx; gy += fy; gz += fz;
                    hit = true;
                }
            }
            if (hit) { sgx[jj] -= fx; sgy[jj] -= fy; sgz[jj] -= fz; }   // slot jj belongs to this thread at this k
            __syncthreads();
        }
        if (active) add_grad(grad, i, gx, gy, gz);
    } else if (active) {
        for (int k = 0; k < TJ; ++k) {
            const int jj = (t + k) & (TJ - 1);
            const int j  = bj * TJ + jj;
            if (j >= N || j <= i) continue;
            double dx = xi - sx[jj];
            double dy = yi - sy[jj];
            double dz = zi - sz[jj];
            applyMIC(dx, dy, dz);
            double r2  = dx*dx + dy*dy + dz*dz;
            double rij = sqrt(r2);
            if (!(rij <= r_cut && rij >= 1e-10)) continue;
            const double qj = sq[jj];
            if (isnan(qj)) continue;
            const double asum = ai + sa[jj];
            double erf_v, dEdr;
            if (r2 >= 36.0 * asum) {   // far field (gamma r >= GFNFF_COUL_FARFIELD_GAMMA_R): erf == 1, no erf/exp
                erf_v = 1.0;
                dEdr  = qi * qj * (-1.0 / (rij * rij));
            } else {
                double g       = 1.0 / sqrt(asum);
                double gamma_r = g * rij;
                erf_v = erf(gamma_r);
                double exp_v = exp(-gamma_r * gamma_r);
                double derf  = g * exp_v * (2.0 * inv_sqrt_pi);
                dEdr = qi * qj * (derf / rij - erf_v / (rij * rij));
            }
            local_E += qi * qj * erf_v / rij;
            const double fac = dEdr / rij;   // force on i = fac*(ri - rj), on j the negative
            const double fx = fac*dx, fy = fac*dy, fz = fac*dz;
            gx += fx; gy += fy; gz += fz;
            if (fixmode) {   // WP10: fixed-point partner gradients in shared memory (order-independent, as CUDA)
                fix_add_u64(reinterpret_cast<unsigned long long*>(&sgx[jj]), -fx, GFNFF_FIX_FORCE);
                fix_add_u64(reinterpret_cast<unsigned long long*>(&sgy[jj]), -fy, GFNFF_FIX_FORCE);
                fix_add_u64(reinterpret_cast<unsigned long long*>(&sgz[jj]), -fz, GFNFF_FIX_FORCE);
            } else {
                atomicAdd(&sgx[jj], -fx);
                atomicAdd(&sgy[jj], -fy);
                atomicAdd(&sgz[jj], -fz);
            }
        }
        if (gx != 0.0 || gy != 0.0 || gz != 0.0 || TI == TJ)   // TI == TJ: former behaviour (always add)
            add_grad(grad, i, gx, gy, gz);
    }
    __syncthreads();
    if (t < TJ) {
        const int jl = bj * TJ + t;
        if (fixmode && !lockstep) {
            // integer partner sums go to the integer global accumulator unchanged (exact, no second rounding)
            const unsigned long long* ux = reinterpret_cast<const unsigned long long*>(&sgx[t]);
            const unsigned long long* uy = reinterpret_cast<const unsigned long long*>(&sgy[t]);
            const unsigned long long* uz = reinterpret_cast<const unsigned long long*>(&sgz[t]);
            if (jl < N && (*ux != 0ull || *uy != 0ull || *uz != 0ull)) {
                atomicAdd(reinterpret_cast<unsigned long long*>(&grad[3*jl]),   *ux);
                atomicAdd(reinterpret_cast<unsigned long long*>(&grad[3*jl+1]), *uy);
                atomicAdd(reinterpret_cast<unsigned long long*>(&grad[3*jl+2]), *uz);
            }
        } else if (jl < N && (sgx[t] != 0.0 || sgy[t] != 0.0 || sgz[t] != 0.0)) {
            add_grad(grad, jl, sgx[t], sgy[t], sgz[t]);
        }
    }

    blockReduceAddEnergy(local_E, energy);
}

/// Claude Generated (Oct 2026): tile shape of k_coulomb_tiles_hip for N atoms. From 2000 atoms on the
/// former 128 x 128 kernel runs unchanged (bit-identical results; at 2910 atoms 32 x 16 was within 5 %,
/// at 6200 atoms slower). Below, j tiles of 16 / 8 / 4 atoms spread the pairs over more workgroups.
/// Kernel time per step on the RX 9060 XT (gfx1200, 128x128 -> chosen shape, 101 calls each):
/// 24 atoms 52 -> 24 us (32x8), 66: 141 -> 27 us, 114: 269 -> 36 us, 300: 182 -> 56 us,
/// 600: 216 -> 131 us (32x16), 1002: 404 -> 233 us (32x16).
/// CURCUMA_ROCM_COUL_TILE=TIxTJ forces a shape (128x128 = former kernel; 64x64, 32x32, 32x16, 32x8, 32x4).
static void coulombTileShapeHip(int N, int& ti, int& tj)
{
    static const int forced = [] {
        const char* e = std::getenv("CURCUMA_ROCM_COUL_TILE");
        if (!e || !*e) return 0;
        const std::string v(e);
        if (v == "128x128") return 1;
        if (v == "64x64")   return 2;
        if (v == "32x32")   return 3;
        if (v == "32x16")   return 4;
        if (v == "32x8")    return 5;
        if (v == "32x4")    return 6;
        return 0;
    }();
    int sel = forced;
    if (sel == 0) {
        if      (N >= 2000) sel = 1;
        else if (N >= 400)  sel = 4;
        else if (N >= 48)   sel = 5;
        else                sel = 6;
    }
    static const int TIs[7] = {128, 128, 64, 32, 32, 32, 32}, TJs[7] = {128, 128, 64, 32, 16, 8, 4};
    ti = TIs[sel]; tj = TJs[sel];
}

/// Claude Generated (Oct 2026): launch k_coulomb_tiles_hip with the shape of coulombTileShapeHip().
static void launchCoulombTilesHip(int N, const double* alpha, const double* cx, const double* cy, const double* cz,
                                  const double* charges, double r_cut, double* grad, double* energy, hipStream_t stream)
{
    int ti = GFNFF_COUL_TILE_HIP, tj = GFNFF_COUL_TILE_HIP;
    coulombTileShapeHip(N, ti, tj);
    const dim3 grid((N + tj - 1) / tj, (N + ti - 1) / ti);
#define COUL_TILES_LAUNCH(A, B) k_coulomb_tiles_hip<A, B><<<grid, A, 0, stream>>>(N, alpha, cx, cy, cz, charges, r_cut, grad, energy)
    if      (ti == 128) COUL_TILES_LAUNCH(128, 128);
    else if (ti == 64)  COUL_TILES_LAUNCH(64, 64);
    else if (tj == 32)  COUL_TILES_LAUNCH(32, 32);
    else if (tj == 16)  COUL_TILES_LAUNCH(32, 16);
    else if (tj == 8)   COUL_TILES_LAUNCH(32, 8);
    else                COUL_TILES_LAUNCH(32, 4);
#undef COUL_TILES_LAUNCH
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
                acc_add(&dEdcn[i], yy * ff[tid] * cnfak_i[tid], GFNFF_FIX_FORCE);
                acc_add(&dEdcn[j], yy * ff[tid] * cnfak_j[tid], GFNFF_FIX_FORCE);
            }

            if (zz_hb && nr_hb[tid] >= 1 && hb_H_atom) {
                constexpr double t1 = 0.1;
                double zz = t1 * alpha_orig * dr * dr * E;
                int H = hb_H_atom[tid];
                if (H >= 0) {
                    acc_add(&zz_hb[H], zz, GFNFF_FIX_FORCE);
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

// Claude Generated (Oct 2026): per-triple charge scaling factors Q_H, Q_A, Q_B of k_hbonds, formed once per
// charge upload (the HB charges are frozen at topology build, or refreshed by k_gather_hb_charges) instead of
// three exp per triple and step. Same function and arguments as in k_hbonds.
/// @param n          number of HB triples
/// @param q_H,q_A,q_B [n] HB charges per triple
/// @param Qs_H,Qs_A,Qs_B [n] output: d_charge_scaling(q_H), d_charge_scaling(-q_A), d_charge_scaling(-q_B)
__global__ void k_hb_charge_scaling(
    int n,
    const double* __restrict__ q_H, const double* __restrict__ q_A, const double* __restrict__ q_B,
    double* __restrict__ Qs_H, double* __restrict__ Qs_A, double* __restrict__ Qs_B)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < n) {
        Qs_H[tid] = d_charge_scaling(q_H[tid], 15.0, 1.0);    // HB_ST_D, HB_SF_D of k_hbonds
        Qs_A[tid] = d_charge_scaling(-q_A[tid], 15.0, 1.0);
        Qs_B[tid] = d_charge_scaling(-q_B[tid], 15.0, 1.0);
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
    const int*    __restrict__ slot_off,
    double*       __restrict__ slot_g,
    const double* __restrict__ Qs_H,
    const double* __restrict__ Qs_A,
    const double* __restrict__ Qs_B,
    double*                    grad,
    double*                    energy)
{
    // Claude Generated (Oct 2026): slot mode (slot_g != nullptr). Instead of FP64 atomicAdd into grad (a CAS
    // loop on gfx1200; the triples of one H/A sit next to each other, so the lanes of a wavefront collide),
    // every gradient contribution is stored in this triple's own slots [slot_off[tid], slot_off[tid+1]):
    // A, H, B, neighbours of B (case >= 2), then C and the neighbours of C (case 3 with a parent). Each slot
    // is written by this thread only; k_slot_gather sums them per atom. Inactive triples zero their slots.
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double local_E = 0.0;
    bool slots_written = false;

    if (tid < n) {
    const int sbase = slot_g ? slot_off[tid] : 0;

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

    // Claude Generated (Oct 2026): Qs_* (k_hb_charge_scaling, same function) when given: 3 exp per triple fewer.
    double Q_H = Qs_H ? Qs_H[tid] : d_charge_scaling(q_H[tid], HB_ST_D, HB_SF_D);
    double Q_A = Qs_A ? Qs_A[tid] : d_charge_scaling(-q_A[tid], HB_ST_D, HB_SF_D);
    double Q_B = Qs_B ? Qs_B[tid] : d_charge_scaling(-q_B[tid], HB_ST_D, HB_SF_D);

    double bas = (Q_A * basicity_A[tid] * r_AH_4 + Q_B * basicity_B[tid] * r_HB_4) * denom_DA;
    double aci = (acidity_B[tid] * r_AH_4 + acidity_A[tid] * r_HB_4) * denom_DA;

    int eA = elem_A[tid] - 1, eB = elem_B[tid] - 1;
    double r_vdw_AB = ((eA >= 0 && eA < 100) ? d_cov_radii[eA] : 1.0)
                    + ((eB >= 0 && eB < 100) ? d_cov_radii[eB] : 1.0);

    // Claude Generated (Oct 2026): d_damping_short_range / d_damping_long_range written out so that the
    // gradient reuses (scut r_vdw / r^2)^6 and (r^2 / longcut)^6 instead of recomputing them (same operations,
    // same values bitwise).
    double p_short = 1.0, p_long = 1.0;
    {
        const double rr = r_AB * r_AB;
        const double ratio_s = HB_SCUT_D * r_vdw_AB / rr;
        const double ratio_l = rr / HB_LONGCUT_D;
        for (int i = 0; i < (int)HB_ALP_D; ++i) { p_short *= ratio_s; p_long *= ratio_l; }
    }
    double damp_short = 1.0 / (1.0 + p_short);
    double damp_long = 1.0 / (1.0 + p_long);
    double damp_env = damp_short * damp_long;

    // Out-of-line damping
    double rahprbh = r_AH + r_HB + 1e-12;
    double expo = (HB_BACUT_D / r_vdw_AB) * (rahprbh / r_AB - 1.0);

    if (expo <= 15.0) {
    const double exp_expo = exp(expo);   // also ratio2 of the gradient
    double damp_outl = 2.0 / (1.0 + exp_expo);

    // Neighbor out-of-line (Case >= 2)
    // Claude Generated (Oct 2026): |A - nb| and |B - nb| of the first HB_NB_CACHE neighbours are kept for the
    // gradient loop below (same values bitwise: the difference vectors there are the negated ones).
    constexpr int HB_NB_CACHE = 4;
    double c_ranb[HB_NB_CACHE], c_rbnb[HB_NB_CACHE];
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
            if (ni < HB_NB_CACHE) { c_ranb[ni] = r_Anb; c_rbnb[ni] = r_Bnb; }
            double expo_nb = (hbnbcut_save / r_vdw_AB) * ((r_Anb + r_Bnb) / r_AB - 1.0);
            outl_nb_tot *= (2.0 / (1.0 + exp(-expo_nb)) - 1.0);
        }
    }

    // Case 4: Virtual LP — save state for gradient
    double outl_lp = 1.0;
    double lp_dist_c4 = 0, lp_vx_c4 = 0, lp_vy_c4 = 0, lp_vz_c4 = 0;
    double vnorm_c4 = 0, ralp_c4 = 0, expo_lp_c4 = 0, exp_lp_c4 = 0;   // exp_lp_c4: Oct 2026, reused by the gradient
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
                exp_lp_c4 = exp(expo_lp_c4);
                outl_lp = 2.0 / (1.0 + exp_lp_c4);
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

    // Damping derivative intermediates (Oct 2026: the powers and exp(expo) come from the energy part)
    const double ratio1_pow = p_long;    // (rab2 / longcut)^6
    const double ratio3_pow = p_short;   // (scut r_vdw / rab2)^6
    double ddamp = (-2.0 * HB_ALP_D * ratio1_pow / (1.0 + ratio1_pow))
                 + ( 2.0 * HB_ALP_D * ratio3_pow / (1.0 + ratio3_pow));

    const double ratio2 = exp_expo;

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
                double ranb, rbnb;
                if (ni < HB_NB_CACHE) { ranb = c_ranb[ni]; rbnb = c_rbnb[ni]; }
                else {
                    ranb = sqrt(dranb_x*dranb_x + dranb_y*dranb_y + dranb_z*dranb_z);
                    rbnb = sqrt(drbnb_x*drbnb_x + drbnb_y*drbnb_y + drbnb_z*drbnb_z);
                }
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
                if (slot_g) {
                    double* g = slot_g + 3 * (sbase + 3 + ni);   // first write of this slot
                    g[0] = -dga_nb_x - dgb_nb_x; g[1] = -dga_nb_y - dgb_nb_y; g[2] = -dga_nb_z - dgb_nb_z;
                } else {
                    add_grad(grad, nb, -dga_nb_x - dgb_nb_x, -dga_nb_y - dgb_nb_y, -dga_nb_z - dgb_nb_z);
                }
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
            // Slot mode: C sits after the nb_B slots, its neighbours D follow (same order as nb_C_flat).
            const int sC = sbase + 3 + nb_B_count[tid];
            if (slot_g) {
                double* g = slot_g + 3 * sC;   // first write of the C slot
                g[0] = bterm_c3 * (-deda_x); g[1] = bterm_c3 * (-deda_y); g[2] = bterm_c3 * (-deda_z);
            } else {
                add_grad(grad, C_idx, bterm_c3 * (-deda_x), bterm_c3 * (-deda_y), bterm_c3 * (-deda_z));
            }
            gh_x += bterm_c3 * (-dedc_x); gh_y += bterm_c3 * (-dedc_y); gh_z += bterm_c3 * (-dedc_z);

            // Torsion gradient: D-B-C-H for each neighbor of C
            double tterm_c3 = -rdamp * qhoutl * eangl * acidity_A[tid] * basicity_B[tid] * Q_A * Q_B * global_scale;
            int nc_off = nb_C_offset[tid];
            int nc_cnt = nb_C_count[tid];

            // Compute individual torsion energies for product rule
            // We need etors/et_k for each k, so we need individual et values
            for (int ni = 0; ni < nc_cnt; ++ni) {
                int D_idx = nb_C_flat[nc_off + ni];
                if (D_idx == B) {
                    if (slot_g) {
                        double* g = slot_g + 3 * (sC + 1 + ni);
                        g[0] = 0.0; g[1] = 0.0; g[2] = 0.0;
                    }
                    continue;
                }

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

                if (slot_g) {
                    double* gD = slot_g + 3 * (sC + 1 + ni);
                    gD[0] = t_k * dij * derivate[0][0]; gD[1] = t_k * dij * derivate[0][1]; gD[2] = t_k * dij * derivate[0][2];
                    double* gC = slot_g + 3 * sC;
                    gC[0] += t_k * dij * derivate[2][0]; gC[1] += t_k * dij * derivate[2][1]; gC[2] += t_k * dij * derivate[2][2];
                } else {
                    add_grad(grad, D_idx, t_k * dij * derivate[0][0], t_k * dij * derivate[0][1], t_k * dij * derivate[0][2]);
                    add_grad(grad, C_idx, t_k * dij * derivate[2][0], t_k * dij * derivate[2][1], t_k * dij * derivate[2][2]);
                }
                gb_x += t_k * dij * derivate[1][0]; gb_y += t_k * dij * derivate[1][1]; gb_z += t_k * dij * derivate[1][2];
                gh_x += t_k * dij * derivate[3][0]; gh_y += t_k * dij * derivate[3][1]; gh_z += t_k * dij * derivate[3][2];
            }
        }

        // Case 4: LP out-of-line gradient (matching CPU calcHydrogenBonds)
        if (ct == 4 && nbb_lp_c4 > 0) {
            double lpterm = -rdamp * Q_H * damp_outl * outl_nb_tot * eangl * etors * const_val;
            double rblp = lp_dist_c4;
            double ralpprblp = ralp_c4 + rblp + 1e-12;
            double ratio2_lp = exp_lp_c4;   // nbb_lp_c4 > 0 implies it was set in the energy part

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
                    if (slot_g) {
                        double* g = slot_g + 3 * (sbase + 3 + ni);   // already written by the out-of-line loop
                        g[0] += -gnb_share_x; g[1] += -gnb_share_y; g[2] += -gnb_share_z;
                    } else {
                        int nb = nb_B_flat[nb_off + ni];
                        add_grad(grad, nb, -gnb_share_x, -gnb_share_y, -gnb_share_z);
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

    if (slot_g) {
        double* g = slot_g + 3 * sbase;   // slots 0, 1, 2 = A, H, B
        g[0] = ga_x; g[1] = ga_y; g[2] = ga_z;
        g[3] = gh_x; g[4] = gh_y; g[5] = gh_z;
        g[6] = gb_x; g[7] = gb_y; g[8] = gb_z;
        slots_written = true;
    } else {
        add_grad(grad, A, ga_x, ga_y, ga_z);
        add_grad(grad, B, gb_x, gb_y, gb_z);
        add_grad(grad, H, gh_x, gh_y, gh_z);
    }

    } // expo <= 15.0
    } // r_AB in range

    // Slot mode: a triple outside the range contributes nothing; clear its slots.
    if (slot_g && !slots_written) {
        for (int s = sbase; s < slot_off[tid + 1]; ++s) {
            slot_g[3*s] = 0.0; slot_g[3*s + 1] = 0.0; slot_g[3*s + 2] = 0.0;
        }
    }
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

    acc_add(energy, E_en + E_self, GFNFF_FIX_ENERGY);
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

// ============================================================================
// D4 C6 and dC6/dCN per pair, factorised over the partner element.
// Claude Generated (Oct 2026). The double sum
//   C6(i,j)        = sum_ri sum_rj gw(i,ri)  gw(j,rj)  C6ref(Zi,Zj,ri,rj)
//   dC6(i,j)/dCN_i = sum_ri sum_rj dgw(i,ri) gw(j,rj)  C6ref(Zi,Zj,ri,rj)
//   dC6(i,j)/dCN_j = sum_ri sum_rj gw(i,ri)  dgw(j,rj) C6ref(Zi,Zj,ri,rj)
// is split into a per-atom half that depends only on atom i and the ELEMENT of j,
//   u(i,e,rj) = sum_ri dgw(i,ri) C6ref(Zi,e,ri,rj),   v(i,e,rj) = sum_ri gw(i,ri) C6ref(Zi,e,ri,rj),
// and a per-pair half C6 = sum_rj v(i,Zj,rj) gw(j,rj), dC6/dCN_i = sum_rj u(i,Zj,rj) gw(j,rj),
// dC6/dCN_j = sum_rj v(i,Zj,rj) dgw(j,rj): 21 FP64 FMA per pair instead of up to 3 x 49.
// Same terms as the double sum (reference entries with |C6ref| < 1e-20 are dropped in both);
// the products are summed in a different order (FP rounding only).
// Elements are compacted to the types present (ctype_elem); u/v layout
// uv[(atom * ntypes + type) * 2*D4_MAX_REF + {0..6: u, 7..13: v}], zero beyond nref.
// Reference: d4param_generator.cpp:computeDC6DCN(), Fortran gfnff_gdisp0.f90:262-305
// ============================================================================

/// Per-atom half of the factorised D4 C6 / dC6dCN sums: one thread per (atom, element type).
/// @param natoms      number of atoms
/// @param ntypes      number of element types present in the system
/// @param atom_types  [N] atomic numbers (1-based)
/// @param ctype_elem  [ntypes] 0-based element index of each compact type
/// @param gw, dgw     [N * D4_MAX_REF] Gaussian weights and their CN derivatives (zero-padded)
/// @param c6_flat     [MAX_ELEM^2 * MAX_REF^2] C6 reference table
/// @param uv          [N * ntypes * 2 * D4_MAX_REF] output u (dgw-contracted) and v (gw-contracted)
__global__ void k_disp_uv(
    int natoms, int ntypes,
    const int*    __restrict__ atom_types,
    const int*    __restrict__ ctype_elem,
    const double* __restrict__ gw,
    const double* __restrict__ dgw,
    const double* __restrict__ c6_flat,
    double*       __restrict__ uv)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= natoms * ntypes) return;
    const int a = k / ntypes;
    const int e = k - a * ntypes;
    double* out = uv + (size_t)k * (2 * D4_MAX_REF);

    const int ea = atom_types[a] - 1;
    const int eb = ctype_elem[e];
    if (ea < 0 || ea >= D4_MAX_ELEM || eb < 0 || eb >= D4_MAX_ELEM) {
        for (int r = 0; r < 2 * D4_MAX_REF; ++r) out[r] = 0.0;
        return;
    }
    int nra = d_refn_const[ea]; if (nra > D4_MAX_REF) nra = D4_MAX_REF;
    int nrb = d_refn_const[eb]; if (nrb > D4_MAX_REF) nrb = D4_MAX_REF;
    const size_t c6_base = (size_t)ea * D4_MAX_ELEM * D4_MAX_REF * D4_MAX_REF
                         + (size_t)eb * D4_MAX_REF * D4_MAX_REF;
    double g[D4_MAX_REF], dg[D4_MAX_REF];
    for (int r = 0; r < D4_MAX_REF; ++r) { g[r] = gw[a * D4_MAX_REF + r]; dg[r] = dgw[a * D4_MAX_REF + r]; }

    for (int rb = 0; rb < D4_MAX_REF; ++rb) {
        double u = 0.0, v = 0.0;
        if (rb < nrb) {
            for (int ra = 0; ra < nra; ++ra) {
                const double c6ref = c6_flat[c6_base + ra * D4_MAX_REF + rb];
                if (fabs(c6ref) < 1e-20) continue;
                u += dg[ra] * c6ref;
                v += g[ra] * c6ref;
            }
        }
        out[rb]              = u;
        out[D4_MAX_REF + rb] = v;
    }
}

/// Per-pair half of the factorised D4 sums (see the block comment above): one thread per pair.
/// @param uv         output of k_disp_uv for the same gw/dgw
/// @param atom_ctype [N] compact element type of each atom
/// @param c6_out     [n] C6(i,j) at the current CN, or nullptr (C6 kept)
__global__ void k_dc6dcn_per_pair_uv(
    int n_pairs, int ntypes,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ atom_ctype,
    const double* __restrict__ uv,
    const double* __restrict__ gw,
    const double* __restrict__ dgw,
    double*       __restrict__ dc6dcn_ij,
    double*       __restrict__ dc6dcn_ji,
    double*       __restrict__ c6_out)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pairs) return;
    const int ai = idx_i[tid];
    const int aj = idx_j[tid];
    const double* u  = uv + ((size_t)ai * ntypes + atom_ctype[aj]) * (2 * D4_MAX_REF);
    const double* v  = u + D4_MAX_REF;
    const double* gj = gw  + aj * D4_MAX_REF;
    const double* dj = dgw + aj * D4_MAX_REF;
    double dc6_ij = 0.0, dc6_ji = 0.0, c6 = 0.0;
#pragma unroll
    for (int r = 0; r < D4_MAX_REF; ++r) {
        const double g = gj[r];
        dc6_ij += u[r] * g;
        dc6_ji += v[r] * dj[r];
        c6     += v[r] * g;
    }
    dc6dcn_ij[tid] = dc6_ij;
    dc6dcn_ji[tid] = dc6_ji;
    if (c6_out) c6_out[tid] = c6;
}

// ---------------------------------------------------------------------------
// WP9-B (Claude Generated, Oct 2026): tabulated dc6dcn, port of the CUDA k_d4_partner_tables /
// k_dc6dcn_per_pair_tab (cuda/gfnff_kernels.cu). Per atom a and partner element type t, from gw only:
//   A[a][t][r'] = sum_r gw_a[r] C6ref(Z_a, Z_t; r, r'),   W[a][t][r'] = sum_r gw_a[r] C6ref(Z_t, Z_a; r', r)
// and per pair C6 = sum gw_j A[i][t_j], dC6/dCN(j) = sum dgw_j A[i][t_j], dC6/dCN(i) = sum dgw_i W[j][t_i].
// Against k_disp_uv / k_dc6dcn_per_pair_uv: C6 and dC6/dCN(j) are the same sums (A == v), dC6/dCN(i)
// contracts the other index first (differs by rounding). The compact types are ROCm's ctype_elem /
// atom_ctype (0-based element per type) instead of CUDA's type_z / atype_idx.
// ---------------------------------------------------------------------------
__global__ void k_d4_partner_tables(int N, int nt, const int* __restrict__ ctype_elem, const int* __restrict__ atom_types,
                                    const double* __restrict__ gw, const double* __restrict__ c6_flat,
                                    double* __restrict__ A, double* __restrict__ W)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N * nt) return;
    const int a = tid / nt, t = tid % nt;
    const int za = atom_types[a] - 1, zt = ctype_elem[t];
    double* outA = A + (static_cast<size_t>(a) * nt + t) * D4_MAX_REF;
    double* outW = W + (static_cast<size_t>(a) * nt + t) * D4_MAX_REF;
    if (za < 0 || za >= D4_MAX_ELEM || zt < 0 || zt >= D4_MAX_ELEM) {
        for (int r = 0; r < D4_MAX_REF; ++r) { outA[r] = 0.0; outW[r] = 0.0; }
        return;
    }
    int nra = d_refn_const[za]; if (nra > D4_MAX_REF) nra = D4_MAX_REF;
    int nrt = d_refn_const[zt]; if (nrt > D4_MAX_REF) nrt = D4_MAX_REF;
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
        outA[rp] = sa;
        outW[rp] = sw;
    }
}

/// Claude Generated (Oct 2026, WP9-B): per-pair half of the tabulated dc6dcn (see the block above).
/// @param c6_out [n] C6(i,j) at the current CN, or nullptr (C6 kept)
__global__ void k_dc6dcn_per_pair_tab(int n_pairs, const int* __restrict__ idx_i, const int* __restrict__ idx_j,
                                      const int* __restrict__ atom_types, const int* __restrict__ atom_ctype, int nt,
                                      const double* __restrict__ gw, const double* __restrict__ dgw,
                                      const double* __restrict__ A, const double* __restrict__ W,
                                      double* __restrict__ dc6dcn_ij, double* __restrict__ dc6dcn_ji,
                                      double* __restrict__ c6_out)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pairs) return;
    const int ai = idx_i[tid], aj = idx_j[tid];
    int nri = d_refn_const[atom_types[ai] - 1]; if (nri > D4_MAX_REF) nri = D4_MAX_REF;
    int nrj = d_refn_const[atom_types[aj] - 1]; if (nrj > D4_MAX_REF) nrj = D4_MAX_REF;
    const double* Aij = A + (static_cast<size_t>(ai) * nt + atom_ctype[aj]) * D4_MAX_REF;
    const double* Wji = W + (static_cast<size_t>(aj) * nt + atom_ctype[ai]) * D4_MAX_REF;
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

// ============================================================================
// WP-A (Jun 2026): on-device D4 dispersion pair-list build (HIP port of
// cuda/gfnff_kernels.cu k_disp_pairs_count / k_disp_pairs_build). Reference target
// (bit-identical up to FP order of the host Gaussian weights):
//   d4param_generator.cpp:GenerateDispersionPairsNative + getChargeWeightedC6
//   GFNFFParameters::zetaChargeScale (gfnff_par.h)
// ============================================================================

/// Device port of GFNFFParameters::zetaChargeScale. Caller passes per-element zeff/c;
/// Z out of [1,86] -> 1.0 at the call site.
__device__ inline double zetaChargeScaleDev(double zeff, double c, double q)
{
    double qmod = zeff + q;
    if (qmod < 0.0) return exp(3.0);            // highly anionic -> max zeta
    return exp(3.0 * (1.0 - exp(c * (1.0 - zeff / qmod))));
}

/// Map a linear upper-triangle index tid -> (i,j), i<j.
__device__ inline void dispPairFromTid(int tid, int N, int& i, int& j)
{
    i = static_cast<int>((2.0 * N - 1.0 - sqrt((2.0 * N - 1.0) * (2.0 * N - 1.0) - 8.0 * tid)) * 0.5);
    while (i > 0 && i * (2 * N - i - 1) / 2 > tid) --i;
    while (i + 1 < N && (i + 1) * (2 * N - (i + 1) - 1) / 2 <= tid) ++i;
    j = tid - i * (2 * N - i - 1) / 2 + i + 1;
}

/// Claude Generated (Oct 2026): per-pair static D4 data of the pair (i, j), shared by
/// k_disp_pairs_build and k_disp_pairs_static (same expressions as before the split).
__device__ inline void dispPairStaticData(
    int i, int j, int zi, int zj,
    const double* __restrict__ gw, const double* __restrict__ c6_flat,
    const double* __restrict__ sqrtzr4r2, const double* __restrict__ topo_q,
    const double* __restrict__ zeta_zeff, const double* __restrict__ zeta_c,
    double a1, double a2,
    double& c6_out, double& r4r2_out, double& r0sq_out, double& zeta_out)
{
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

    c6_out   = c6;
    r4r2_out = r4r2ij;
    r0sq_out = r0_sq;
    zeta_out = zsi * zsj;
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
    double r2 = dx * dx + dy * dy + dz * dz;
    if (r2 > cutoff_sq) return;

    atomicAdd(d_count, 1);
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
    double r2 = dx * dx + dy * dy + dz * dz;
    if (r2 > cutoff_sq) return;

    double c6, r4r2ij, r0_sq, zeta;
    dispPairStaticData(i, j, zi, zj, gw, c6_flat, sqrtzr4r2, topo_q, zeta_zeff, zeta_c, a1, a2,
                       c6, r4r2ij, r0_sq, zeta);

    int slot = atomicAdd(d_counter, 1);
    out_i[slot]    = i;
    out_j[slot]    = j;
    out_C6[slot]   = c6;
    out_r4r2[slot] = r4r2ij;
    out_r0sq[slot] = r0_sq;
    out_zeta[slot] = zeta;
    out_rcut[slot] = 50.0;  // D4 pair-energy cutoff (struct default, matches host)
}

// ============================================================================
// Claude Generated (Oct 2026): ordered device build of the D4 pair list and of its per-atom
// gather CSR. k_disp_pairs_build above takes its output slot from a global atomic counter, so
// the pair order depended on the arrival order of the wavefronts; the CSR then had to be built
// on the host from a downloaded copy of the indices (13.7 M pairs at 6200 atoms: 137 ms of host
// loops, 13-26 ms download, 10 ms upload per rebuild). Here one workgroup owns one atom row:
//   k_disp_rows_count : lo[a] = #partners p < a, hi[a] = #partners p > a in the list
//   (host)            : exclusive scans -> pair row offsets rowoff (over hi), csr_off (over lo+hi)
//   k_disp_rows_build : partners in ascending order (wavefront ballot + per-tile prefix); pair
//                       (a, p > a) goes to slot rowoff[a] + rank, i.e. the host order of
//                       GenerateDispersionPairsNative (i ascending, j ascending = key i*N+j);
//                       CSR row of a = all partners ascending, is_i = (p > a)
//   k_disp_pairs_static: C6, r4r2ij, R0^2, zetac6, r_cut per pair (one thread per pair)
//   k_disp_csr_lower  : csr_pairidx of the p < a half by binary search of a in the
//                       (ascending) idx_j row of p
// The CSR rows are those DispersionSoA::upload builds from the same pair list: there each row
// receives its entries in pair order, which for a list sorted by i*N+j is ascending partner order.
// ============================================================================

/// Same predicate as k_disp_pairs_count / k_disp_pairs_build for the pair (i, j), i < j.
__device__ inline bool dispPairInList(int i, int j,
                                      const double* __restrict__ cx, const double* __restrict__ cy,
                                      const double* __restrict__ cz, const int* __restrict__ atom_types,
                                      double cutoff_sq)
{
    const int zi = atom_types[i];
    const int zj = atom_types[j];
    if (zi < 1 || zi > D4_MAX_ELEM || zj < 1 || zj > D4_MAX_ELEM) return false;
    const double dx = cx[i] - cx[j];
    const double dy = cy[i] - cy[j];
    const double dz = cz[i] - cz[j];
    const double r2 = dx * dx + dy * dy + dz * dz;
    return !(r2 > cutoff_sq);
}

constexpr int DISP_ROW_BLOCK = 256;

/// Claude Generated (Oct 2026): per-atom partner counts below (lo) and above (hi) the atom index.
__global__ __launch_bounds__(DISP_ROW_BLOCK) void k_disp_rows_count(
    int N,
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const int* __restrict__ atom_types, double cutoff_sq,
    int* __restrict__ lo, int* __restrict__ hi)
{
    __shared__ int s_lo[DISP_ROW_BLOCK];
    __shared__ int s_hi[DISP_ROW_BLOCK];
    const int a = blockIdx.x;
    if (a >= N) return;
    int nlo = 0, nhi = 0;
    for (int p = threadIdx.x; p < N; p += blockDim.x) {
        if (p == a) continue;
        const int i = p < a ? p : a;
        const int j = p < a ? a : p;
        if (dispPairInList(i, j, cx, cy, cz, atom_types, cutoff_sq)) {
            if (p < a) ++nlo; else ++nhi;
        }
    }
    s_lo[threadIdx.x] = nlo;
    s_hi[threadIdx.x] = nhi;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (static_cast<int>(threadIdx.x) < s) {
            s_lo[threadIdx.x] += s_lo[threadIdx.x + s];
            s_hi[threadIdx.x] += s_hi[threadIdx.x + s];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) { lo[a] = s_lo[0]; hi[a] = s_hi[0]; }
}

/// Claude Generated (Oct 2026): writes the pair indices (a, p > a) in ascending p at rowoff[a] and
/// the CSR row of atom a (all partners ascending). csr_pairidx of the p < a half: k_disp_csr_lower;
/// per-pair static data: k_disp_pairs_static (one thread per pair: here most lanes of a tile have
/// no pair, and the C6 contraction took 17-27 of 20-31 ms at 13.7 M pairs when done in this kernel).
__global__ __launch_bounds__(DISP_ROW_BLOCK) void k_disp_rows_build(
    int N,
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const int* __restrict__ atom_types, double cutoff_sq,
    const int* __restrict__ lo, const int* __restrict__ rowoff, const int* __restrict__ csr_off,
    int* __restrict__ out_i, int* __restrict__ out_j,
    int* __restrict__ csr_partner, int* __restrict__ csr_pairidx, int* __restrict__ csr_is_i)
{
    __shared__ int s_scan[DISP_ROW_BLOCK / 32];   // per-wavefront counts of the current tile
    const int a = blockIdx.x;
    if (a >= N) return;
    const int t = threadIdx.x;
    const int row_csr = csr_off[a];
    const int row_pair = rowoff[a];
    const int n_lo = lo[a];
    int base = 0;   // partners of a found in the previous tiles
    for (int p0 = 0; p0 < N; p0 += DISP_ROW_BLOCK) {
        const int p = p0 + t;
        bool in = false;
        if (p < N && p != a) {
            const int i = p < a ? p : a;
            const int j = p < a ? a : p;
            in = dispPairInList(i, j, cx, cy, cz, atom_types, cutoff_sq);
        }
        // tile prefix count of the flags: ballot + popcount within the wavefront, then the
        // wavefront totals (warpSize 32 or 64, at most DISP_ROW_BLOCK / 32 of them)
        const int lane = t % warpSize;
        const int wave = t / warpSize;
        const int nwave = DISP_ROW_BLOCK / warpSize;
        const unsigned long long mask = __ballot(in ? 1 : 0);
        const int in_wave = __popcll(mask & ((1ull << lane) - 1ull));
        if (lane == 0) s_scan[wave] = __popcll(mask);
        __syncthreads();
        int before = 0, tile_total = 0;
        for (int w = 0; w < nwave; ++w) {
            if (w < wave) before += s_scan[w];
            tile_total += s_scan[w];
        }
        if (in) {
            const int rank = base + before + in_wave;       // rank of p among the partners of a
            const int e = row_csr + rank;
            csr_partner[e] = p;
            csr_is_i[e]    = (p > a) ? 1 : 0;
            if (p > a) {
                const int k = row_pair + rank - n_lo;
                out_i[k]    = a;
                out_j[k]    = p;
                csr_pairidx[e] = k;
            }
        }
        base += tile_total;
        __syncthreads();   // s_scan is rewritten by the next tile
    }
}

/// Claude Generated (Oct 2026): per-pair static D4 data (C6 contraction, r4r2ij, R0^2, zetac6,
/// r_cut) for the ordered list, one thread per pair.
__global__ GFNFF_KERNEL_BOUNDS void k_disp_pairs_static(
    int n, const int* __restrict__ idx_i, const int* __restrict__ idx_j,
    const int* __restrict__ atom_types,
    const double* __restrict__ gw, const double* __restrict__ c6_flat,
    const double* __restrict__ sqrtzr4r2, const double* __restrict__ topo_q,
    const double* __restrict__ zeta_zeff, const double* __restrict__ zeta_c,
    double a1, double a2,
    double* __restrict__ out_C6, double* __restrict__ out_r4r2, double* __restrict__ out_r0sq,
    double* __restrict__ out_zeta, double* __restrict__ out_rcut)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= n) return;
    const int i = idx_i[k];
    const int j = idx_j[k];
    double c6, r4r2ij, r0_sq, zeta;
    dispPairStaticData(i, j, atom_types[i], atom_types[j], gw, c6_flat, sqrtzr4r2, topo_q,
                       zeta_zeff, zeta_c, a1, a2, c6, r4r2ij, r0_sq, zeta);
    out_C6[k]   = c6;
    out_r4r2[k] = r4r2ij;
    out_r0sq[k] = r0_sq;
    out_zeta[k] = zeta;
    out_rcut[k] = 50.0;  // D4 pair-energy cutoff (struct default, matches host)
}

/// Claude Generated (Oct 2026): pair index of the CSR entries (a, p < a): the pair (p, a) sits
/// in the row of p, whose idx_j are ascending.
__global__ __launch_bounds__(DISP_ROW_BLOCK) void k_disp_csr_lower(
    int N, const int* __restrict__ lo, const int* __restrict__ rowoff, const int* __restrict__ csr_off,
    const int* __restrict__ idx_j, const int* __restrict__ csr_partner, int* __restrict__ csr_pairidx)
{
    const int a = blockIdx.x;
    if (a >= N) return;
    const int s = csr_off[a];
    const int n_lo = lo[a];
    for (int t = threadIdx.x; t < n_lo; t += blockDim.x) {
        const int p = csr_partner[s + t];
        int l = rowoff[p], r = rowoff[p + 1];        // search a in idx_j[l, r)
        while (l < r) {
            const int m = l + ((r - l) >> 1);
            if (idx_j[m] < a) l = m + 1; else r = m;
        }
        csr_pairidx[s + t] = l;                       // idx_j[l] == a by construction
    }
}

// ============================================================================
// WP9-A (Claude Generated, Oct 2026): device rebuild of the non-bonded repulsion list, the ordered
// row scheme of the D4 list above applied to the pair predicate of GFNFF::make_nb_rep() (CUDA:
// k_rep_pairs_count / k_rep_pairs_build, arrival order). Pairs (i < j) come out in ascending i*N+j,
// and the gather CSR row of atom a lists all its partners ascending, which is what
// buildRepulsionGatherCSR() makes of that list. Tables: GFNFF::RepulsionDeviceInputs.
// ============================================================================

/// Device copies of GFNFF::RepulsionDeviceInputs (pointers into the workspace buffers).
struct RepDevTablesHip {
    const int* valid = nullptr;
    const int* metal = nullptr;
    const double* dum = nullptr;
    const double* repz = nullptr;
    const long long* bond_keys = nullptr;
    int n_bond_keys = 0;
    const long long* hh_keys = nullptr;
    const int* hh_dist = nullptr;
    int n_hh_keys = 0;
    double repscaln = 0.0, hhfac = 0.0, hh13 = 0.0, hh14 = 0.0, r_cut = 20.0;
};

/// Index of key in the sorted array keys[0, n), or -1.
__device__ inline int repKeyFindHip(const long long* __restrict__ keys, int n, long long key)
{
    int lo = 0, hi = n - 1;
    while (lo <= hi) {
        const int mid = lo + ((hi - lo) >> 1);
        const long long v = keys[mid];
        if (v == key) return mid;
        if (v < key) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}

/// Pair (i < j) is in the non-bonded list: valid elements, r^2 <= cutoff^2, not an (i, j) bond-list entry.
__device__ inline bool repPairInListHip(int N, int i, int j, const double* __restrict__ cx,
                                        const double* __restrict__ cy, const double* __restrict__ cz,
                                        double cutoff_sq, const RepDevTablesHip& t)
{
    if (!t.valid[i] || !t.valid[j]) return false;
    const double dx = cx[i] - cx[j], dy = cy[i] - cy[j], dz = cz[i] - cz[j];
    if (dx * dx + dy * dy + dz * dz > cutoff_sq) return false;
    return repKeyFindHip(t.bond_keys, t.n_bond_keys, static_cast<long long>(i) * N + j) < 0;
}

/// Claude Generated (Oct 2026): per-atom partner counts below (lo) and above (hi) the atom index.
__global__ __launch_bounds__(DISP_ROW_BLOCK) void k_rep_rows_count(
    int N, const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    double cutoff_sq, RepDevTablesHip t, int* __restrict__ lo, int* __restrict__ hi)
{
    __shared__ int s_lo[DISP_ROW_BLOCK];
    __shared__ int s_hi[DISP_ROW_BLOCK];
    const int a = blockIdx.x;
    if (a >= N) return;
    int nlo = 0, nhi = 0;
    for (int p = threadIdx.x; p < N; p += blockDim.x) {
        if (p == a) continue;
        const int i = p < a ? p : a;
        const int j = p < a ? a : p;
        if (repPairInListHip(N, i, j, cx, cy, cz, cutoff_sq, t)) {
            if (p < a) ++nlo; else ++nhi;
        }
    }
    s_lo[threadIdx.x] = nlo;
    s_hi[threadIdx.x] = nhi;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (static_cast<int>(threadIdx.x) < s) {
            s_lo[threadIdx.x] += s_lo[threadIdx.x + s];
            s_hi[threadIdx.x] += s_hi[threadIdx.x + s];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) { lo[a] = s_lo[0]; hi[a] = s_hi[0]; }
}

/// Claude Generated (Oct 2026): pair indices (a, p > a) in ascending p at rowoff[a]; with csr_partner set
/// also the gather CSR row of a (all partners ascending; pairidx of the p < a half: k_disp_csr_lower).
__global__ __launch_bounds__(DISP_ROW_BLOCK) void k_rep_rows_build(
    int N, const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    double cutoff_sq, RepDevTablesHip t,
    const int* __restrict__ lo, const int* __restrict__ rowoff, const int* __restrict__ csr_off,
    int* __restrict__ out_i, int* __restrict__ out_j, int* __restrict__ csr_partner, int* __restrict__ csr_pairidx)
{
    __shared__ int s_scan[DISP_ROW_BLOCK / 32];   // per-wavefront counts of the current tile (warpSize 32 or 64)
    const int a = blockIdx.x;
    if (a >= N) return;
    const int tt = threadIdx.x;
    const int row_pair = rowoff[a];
    const int n_lo = lo[a];
    const int row_csr = csr_partner ? csr_off[a] : 0;
    int base = 0;
    for (int p0 = 0; p0 < N; p0 += DISP_ROW_BLOCK) {
        const int p = p0 + tt;
        bool in = false;
        if (p < N && p != a) {
            const int i = p < a ? p : a;
            const int j = p < a ? a : p;
            in = repPairInListHip(N, i, j, cx, cy, cz, cutoff_sq, t);
        }
        const int lane = tt % warpSize;
        const int wave = tt / warpSize;
        const int nwave = DISP_ROW_BLOCK / warpSize;
        const unsigned long long mask = __ballot(in ? 1 : 0);
        const int in_wave = __popcll(mask & ((1ull << lane) - 1ull));
        if (lane == 0) s_scan[wave] = __popcll(mask);
        __syncthreads();
        int before = 0, tile_total = 0;
        for (int w = 0; w < nwave; ++w) {
            if (w < wave) before += s_scan[w];
            tile_total += s_scan[w];
        }
        if (in) {
            const int rank = base + before + in_wave;
            if (csr_partner) csr_partner[row_csr + rank] = p;
            if (p > a) {
                const int k = row_pair + rank - n_lo;
                out_i[k] = a;
                out_j[k] = p;
                if (csr_partner) csr_pairidx[row_csr + rank] = k;
            }
        }
        base += tile_total;
        __syncthreads();
    }
}

/// Claude Generated (Oct 2026): alpha, repab, r_cut of every pair, same expressions and branch order as
/// make_nb_rep() (alpha = sqrt(dum_i dum_j) ff, repab = repz_i repz_j REPSCALN).
__global__ GFNFF_KERNEL_BOUNDS void k_rep_pairs_static(
    int n, int N, const int* __restrict__ idx_i, const int* __restrict__ idx_j,
    const int* __restrict__ atom_types, RepDevTablesHip t,
    double* __restrict__ out_alpha, double* __restrict__ out_repab, double* __restrict__ out_rcut)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= n) return;
    const int i = idx_i[k], j = idx_j[k];
    const int zi = atom_types[i], zj = atom_types[j];
    double ff = 1.0;
    if (zi == 1 && zj == 1) {
        ff = t.hhfac;
        const int h = repKeyFindHip(t.hh_keys, t.n_hh_keys, static_cast<long long>(i) * N + j);
        if (h >= 0) {
            if (t.hh_dist[h] == 2) ff *= t.hh13;
            else if (t.hh_dist[h] == 3) ff *= t.hh14;
        }
    } else if ((zi == 1 && t.metal[j]) || (zj == 1 && t.metal[i])) {
        ff = 0.85;
    } else if ((zi == 1 && zj == 6) || (zj == 1 && zi == 6)) {
        ff = 0.91;
    } else if ((zi == 1 && zj == 8) || (zj == 1 && zi == 8)) {
        ff = 1.04;
    }
    out_alpha[k] = sqrt(t.dum[i] * t.dum[j]) * ff;
    out_repab[k] = t.repz[i] * t.repz[j] * t.repscaln;
    out_rcut[k]  = t.r_cut;
}

// ============================================================================
// Device-side H-bond detection (Claude Generated, Oct 2026): port of the CUDA kernels of
// cuda/gfnff_kernels.cu (d_hb_* helpers and k_hb_to_float unchanged; AB pairs per row (k_hb_ab_rows), k_hb_detect with
// row counting; ordering by rows instead of a sort of 64-bit keys, see FFWorkspaceHip::detectHBondsOnDevice()).
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
    const double dx = cx[a] - cx[b], dy = cy[a] - cy[b], dz = cz[a] - cz[b];
    return dx * dx + dy * dy + dz * dz;
}

/// |E_HB| of a case-1 triple, the formula of estimateHBStrengthCase1() (and of the energy kernel).
__device__ double d_hb_case1_strength(const double* __restrict__ cx, const double* __restrict__ cy,
                                       const double* __restrict__ cz, const int* __restrict__ Z,
                                       const double* __restrict__ bas, const double* __restrict__ aci,
                                       const double* __restrict__ q, int A, int H, int B, double thr)
{
    const double r_AH = sqrt(d_hb_dist2(cx, cy, cz, A, H));
    const double r_HB = sqrt(d_hb_dist2(cx, cy, cz, H, B));
    const double r_AB = sqrt(d_hb_dist2(cx, cy, cz, A, B));
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

/// Claude Generated (Oct 2026): float copies of the HB hydrogen positions in hyd_list order, packed (x, y, z, unused)
/// for the case-1 loop of k_hb_detect_rows: one coalesced 16-byte load per candidate instead of the index load
/// followed by three dependent gathers (k_hb_detect_rows 10.8 -> 8.8 ms at 6200 atoms on the RX 9060 XT).
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_hb_pack_hyd(int n_hyd, const int* __restrict__ hyd_list,
                                                        const float* __restrict__ fx, const float* __restrict__ fy,
                                                        const float* __restrict__ fz, float4* __restrict__ hpos)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= n_hyd) return;
    const int H = hyd_list[k];
    hpos[k] = make_float4(fx[H], fy[H], fz[H], 0.f);
}

/// FP32 evaluation of the case-1 strength of d_hb_case1_strength(): returns -1 (surely below thr), +1 (surely at
/// or above thr) or 0 (within a 1 % band around thr or a borderline out-of-line exponent: the caller runs the
/// double path). The float evaluation is good to about 1e-5 relative, the 1 % band covers it with a wide margin.
__device__ __forceinline__ int d_hb_strength_decision_f(
    const float* __restrict__ fx, const float* __restrict__ fy, const float* __restrict__ fz,
    const int* __restrict__ Z, const double* __restrict__ bas, const double* __restrict__ aci,
    const double* __restrict__ q, int A, int H, int B, float thr)
{
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

/// Claude Generated (Oct 2026): AB pairs of the H-bond detection, same predicate as the CUDA k_hb_ab_pairs, but one
/// workgroup per atom i (partners j > i over the threads) instead of one thread per upper-triangle index: no
/// double-precision sqrt for the index mapping and no work at all for atoms that cannot be donor or acceptor
/// (one thread per pair: 2.9 ms, per row: 0.3 ms on the RX 9060 XT at 6200 atoms). Arrival order (atomic slot).
__global__ __launch_bounds__(256) void k_hb_ab_rows(
    int N, const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const double* __restrict__ bas, const double* __restrict__ aci, const char* __restrict__ ab_ok,
    double thr1, int cap, int* d_counter, int* ab_i, int* ab_j)
{
    const int i = blockIdx.x;
    if (i >= N || !ab_ok[i]) return;
    const double bi = bas[i], ai = aci[i];
    for (int j = i + 1 + static_cast<int>(threadIdx.x); j < N; j += blockDim.x) {
        if (!ab_ok[j]) continue;
        if (d_hb_dist2(cx, cy, cz, i, j) > thr1) continue;
        const double s1 = bi * aci[j], s2 = bas[j] * ai;
        if (s1 < 1e-6 && s2 < 1e-6) continue;
        const int pos = atomicAdd(d_counter, 1);
        if (pos < cap) { ab_i[pos] = i; ab_j[pos] = j; }
    }
}

/// Claude Generated (Oct 2026): threads per AB pair (one workgroup per pair) in k_hb_detect_rows and capacity of its
/// shared-memory candidate queue. The kernel synchronises with __syncthreads and shared-memory atomics only (no
/// shuffle, ballot or other wavefront intrinsic), so it is correct for wave32 (gfx11xx, gfx12xx) and wave64 alike.
constexpr int HB_DETECT_THREADS = 128;
constexpr int HB_DETECT_QUEUE = 1024;
constexpr int HB_DETECT_PER_THREAD = 4;   ///< hydrogens per thread and pre-filter round (queue holds two rounds)

/// Claude Generated (Oct 2026): the H-bond triples of the CUDA k_hb_detect (cuda/gfnff_kernels.cu, same tests in the
/// same order), one workgroup per AB pair. Thread 0 emits the nhb2 triples (cases 2-4). For case 1 the threads first
/// apply the float pre-filter to the HB hydrogens and queue the survivors in shared memory (about 1 % of the
/// hydrogens at 6200 atoms); the queue is then worked off by all threads (exact double tests, bonded tests, strength
/// decision), so the slow path does not run on one lane while the others of its wavefront wait (8.8 -> 4.2 ms at
/// 6200 atoms, 76k AB pairs x 3600 hydrogens, RX 9060 XT). Emission: the
/// workgroup reserves its slots with one global atomic per queue round. Each triple goes to an unordered slot together
/// with its row, row = (case - 1) * N + A (case_rows) or A, whose length is counted in row_count;
/// FFWorkspaceHip::detectHBondsOnDevice() orders the rows afterwards. Counting continues past `cap` (the caller reruns).
__global__ __launch_bounds__(HB_DETECT_THREADS) void k_hb_detect_rows(
    int n_ab, const int* __restrict__ ab_i, const int* __restrict__ ab_j,
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const float* __restrict__ fx, const float* __restrict__ fy, const float* __restrict__ fz, int N,
    const int* __restrict__ Z, const double* __restrict__ bas, const double* __restrict__ aci,
    const double* __restrict__ q,
    const int* __restrict__ hyd_list, const float4* __restrict__ hpos, int n_hyd,
    const int* __restrict__ bond_off, const int* __restrict__ bond_idx,
    const int* __restrict__ hyd_off, const int* __restrict__ hyd_idx,
    const int* __restrict__ nbr_off, const int* __restrict__ nbr_idx,
    double thr2, double min_pair_energy, int case_rows,
    int cap, int* d_counter, int* row_count,
    int* e_A, int* e_H, int* e_B, int* e_ct, int* e_par)
{
    __shared__ int s_queue[HB_DETECT_QUEUE];
    __shared__ int s_qn, s_emit, s_base;
    const int pair = blockIdx.x;
    if (pair >= n_ab) return;   // uniform per workgroup
    const int tid = threadIdx.x;
    const int i = ab_i[pair], j = ab_j[pair];
    const double r_AB_sq = d_hb_dist2(cx, cy, cz, i, j);

    auto store = [&](int pos, int A, int H, int B, int ct, int par) {
        atomicAdd(row_count + (case_rows ? (ct - 1) * N + A : A), 1);
        if (pos < cap) {
            e_A[pos] = A; e_H[pos] = H; e_B[pos] = B; e_ct[pos] = ct; e_par[pos] = par;
        }
    };

    // nhb2: hydrogen bonded to the donor (make_nhb2): case 2, 3 (carbonyl-like acceptor O) or 4 (two-coordinate N)
    if (tid == 0) {
        s_qn = 0; s_emit = 0;
        if (!d_hb_bonded(i, j, bond_off, bond_idx)) {
            auto emit_nhb2 = [&](int donor, int H, int acc) {
                int ct = 2, par = -1;
                const int deg = nbr_off[acc + 1] - nbr_off[acc];
                if (Z[acc] == 8 && deg == 1) {
                    const int parent = nbr_idx[nbr_off[acc]];
                    if (Z[parent] == 6 || Z[parent] == 7) { ct = 3; par = parent; }
                } else if (Z[acc] == 7 && deg == 2) {
                    ct = 4;
                }
                store(atomicAdd(d_counter, 1), donor, H, acc, ct, par);
            };
            for (int p = hyd_off[i]; p < hyd_off[i + 1]; ++p) emit_nhb2(i, hyd_idx[p], j);
            for (int p = hyd_off[j]; p < hyd_off[j + 1]; ++p) {
                const int H = hyd_idx[p];
                if (d_hb_bonded(i, H, bond_off, bond_idx)) continue;  // bonded to i takes precedence
                emit_nhb2(j, H, i);
            }
        }
    }
    __syncthreads();

    // case 1: hydrogen bonded to neither, sum-of-squared-distances criterion (gfnff_ini2.f90:742). Float pre-filter
    // (relative margin 1e-4) into the queue, then the double test that decides as on the host.
    const float fxi = fx[i], fyi = fy[i], fzi = fz[i], fxj = fx[j], fyj = fy[j], fzj = fz[j];
    const float f_rab = static_cast<float>(r_AB_sq);
    const float f_reject = static_cast<float>(thr2) * 1.0001f;
    // exact case-1 test of queued hydrogen hq (true: emit (i, H, j))
    auto exact = [&](int hq) -> bool {
        const int H = hyd_list[hq];
        if (H == i) return false;                              // forEachNeighbor excludes the query atom
        const double r_iH_sq = d_hb_dist2(cx, cy, cz, i, H);
        if (r_iH_sq > thr2) return false;                      // cell query radius sqrt(thr2)
        const double r_jH_sq = d_hb_dist2(cx, cy, cz, j, H);
        if (!(r_AB_sq + r_iH_sq + r_jH_sq < thr2)) return false;
        if (d_hb_bonded(i, H, bond_off, bond_idx) || d_hb_bonded(j, H, bond_off, bond_idx)) return false;
        if (min_pair_energy > 0.0) {
            const int dec = d_hb_strength_decision_f(fx, fy, fz, Z, bas, aci, q, i, H, j, static_cast<float>(min_pair_energy));
            if (dec < 0) return false;
            if (dec == 0) {                                    // within the float margin: the exact (double) value decides
                const double est = d_hb_case1_strength(cx, cy, cz, Z, bas, aci, q, i, H, j, min_pair_energy);
                if (est < min_pair_energy) return false;       // NaN keeps the entry, as on the host
            }
        }
        return true;
    };
    constexpr int round = HB_DETECT_THREADS * HB_DETECT_PER_THREAD;
    for (int h0 = 0; h0 < n_hyd; h0 += round) {
        for (int u = 0; u < HB_DETECT_PER_THREAD; ++u) {
            const int hq = h0 + u * HB_DETECT_THREADS + tid;
            if (hq >= n_hyd) break;
            const float4 hp = hpos[hq];
            const float ax = hp.x - fxi, ay = hp.y - fyi, az = hp.z - fzi;
            const float f_ih = ax * ax + ay * ay + az * az;
            const float bx = hp.x - fxj, by = hp.y - fyj, bz = hp.z - fzj;
            if (!(f_ih > f_reject) && !(f_rab + f_ih + (bx * bx + by * by + bz * bz) > f_reject))
                s_queue[atomicAdd(&s_qn, 1)] = hq;
        }
        __syncthreads();
        const int qn = s_qn;   // uniform after the barrier
        __syncthreads();       // every thread has read s_qn before the next round appends to it
        if (qn + round <= HB_DETECT_QUEUE && h0 + round < n_hyd) continue;   // room for one more round
        // work off the queue in rounds of HB_DETECT_THREADS; one global atomic per round reserves the slots
        for (int r0 = 0; r0 < qn; r0 += HB_DETECT_THREADS) {
            const int x = r0 + tid;
            const int hq_x = (x < qn) ? s_queue[x] : -1;
            const bool emit = (hq_x >= 0) && exact(hq_x);
            const int local = emit ? atomicAdd(&s_emit, 1) : -1;
            __syncthreads();
            if (tid == 0) { s_base = (s_emit > 0) ? atomicAdd(d_counter, s_emit) : 0; }
            __syncthreads();
            if (emit) store(s_base + local, i, hyd_list[hq_x], j, 1, -1);
            __syncthreads();
            if (tid == 0) s_emit = 0;
            __syncthreads();
        }
        if (tid == 0) s_qn = 0;
        __syncthreads();
    }
}
// ---------------------------------------------------------------------------
// Claude Generated (Oct 2026): exclusive scan of an int array on the device (out has n + 1 entries, out[n] is the
// total; *total_out receives it as well when given). Three launches: tiles of HB_SCAN_TILE elements, the tile sums
// in one workgroup, the tile offsets added. Shared-memory scan with __syncthreads only (no wavefront intrinsic).
// ---------------------------------------------------------------------------
constexpr int HB_SCAN_THREADS = 256;
constexpr int HB_SCAN_ITEMS = 4;
constexpr int HB_SCAN_TILE = HB_SCAN_THREADS * HB_SCAN_ITEMS;

/// Exclusive prefix of v over the workgroup (HB_SCAN_THREADS threads); *total receives the sum.
__device__ inline int hbBlockExclusiveScan(int v, int* s, int* total)
{
    const int t = threadIdx.x;
    s[t] = v;
    __syncthreads();
    for (int off = 1; off < HB_SCAN_THREADS; off <<= 1) {
        const int x = (t >= off) ? s[t - off] : 0;
        __syncthreads();
        s[t] += x;
        __syncthreads();
    }
    const int incl = s[t];
    *total = s[HB_SCAN_THREADS - 1];
    __syncthreads();   // s is reused by the caller
    return incl - v;
}

__global__ __launch_bounds__(HB_SCAN_THREADS) void k_hb_scan_tiles(int n, const int* __restrict__ in,
                                                                   int* __restrict__ out, int* __restrict__ tile_sums)
{
    __shared__ int s[HB_SCAN_THREADS];
    const int base = blockIdx.x * HB_SCAN_TILE + threadIdx.x * HB_SCAN_ITEMS;
    int v[HB_SCAN_ITEMS];
    int sum = 0;
    for (int k = 0; k < HB_SCAN_ITEMS; ++k) {
        v[k] = (base + k < n) ? in[base + k] : 0;
        sum += v[k];
    }
    int total = 0;
    int run = hbBlockExclusiveScan(sum, s, &total);
    for (int k = 0; k < HB_SCAN_ITEMS; ++k) {
        if (base + k < n) out[base + k] = run;
        run += v[k];
    }
    if (threadIdx.x == 0) tile_sums[blockIdx.x] = total;
}

__global__ __launch_bounds__(HB_SCAN_THREADS) void k_hb_scan_tile_sums(int ntiles, int* __restrict__ tile_sums,
                                                                       int* __restrict__ out, int n,
                                                                       int* __restrict__ total_out)
{
    __shared__ int s[HB_SCAN_THREADS];
    int carry = 0;
    for (int b0 = 0; b0 < ntiles; b0 += HB_SCAN_THREADS) {
        const int b = b0 + static_cast<int>(threadIdx.x);
        const int v = (b < ntiles) ? tile_sums[b] : 0;
        int total = 0;
        const int ex = hbBlockExclusiveScan(v, s, &total);
        if (b < ntiles) tile_sums[b] = carry + ex;
        carry += total;
    }
    if (threadIdx.x == 0) {
        out[n] = carry;
        if (total_out) *total_out = carry;
    }
}

__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_hb_scan_add(int n, int* __restrict__ out, const int* __restrict__ tile_sums)
{
    const int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) out[t] += tile_sums[t / HB_SCAN_TILE];
}

// ---------------------------------------------------------------------------
// Claude Generated (Oct 2026): ordering of rows. The scatter puts each element into its row in arrival order;
// k_hb_row_sort then sorts each row (one workgroup per row) by its key: (H, B) for the triples (key H*N+B), the
// slot index for the per-atom slot lists. Keys within a row are unique, so the result does not depend on the
// arrival order. Bitonic network in its "flip" form (every comparator puts the smaller key at the lower index), so a
// row of length L that is not a power of two is sorted by skipping every comparator whose upper index is >= L
// (those positions act as +infinity and never move). Rows up to the launch's shared capacity (the next power of two
// of the longest row, at most HB_SORT_SHM entries) are sorted in shared memory,
// longer ones in place in global memory with the keys recomputed from the values.
// ---------------------------------------------------------------------------
constexpr int HB_SORT_THREADS = 256;
constexpr int HB_SORT_SHM = 4096;   ///< largest row sorted in shared memory (48 KB for the triples)

/// Element e of a row: atomic slot within its row (arrival order).
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_hb_row_scatter(int n, const int* __restrict__ row_of_a,
                                                           const int* __restrict__ row_of_ct, int N, int case_rows,
                                                           const int* __restrict__ row_off, int* __restrict__ row_fill,
                                                           int* __restrict__ out)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= n) return;
    // triples: row = (case - 1) * N + A (case_rows) or A; slots (row_of_ct == nullptr): row = atom
    const int a = row_of_a[e];
    if (a < 0 || a >= N) return;   // slots of atoms outside [0, N) are not gathered (buildHBSlotTables)
    const int row = (row_of_ct && case_rows) ? (row_of_ct[e] - 1) * N + a : a;
    out[row_off[row] + atomicAdd(row_fill + row, 1)] = e;
}

/// Key of value v: H*N+B of triple v (eH != nullptr) or v itself (slot index).
__device__ __forceinline__ long long hbRowKey(int v, const int* __restrict__ eH, const int* __restrict__ eB, int N)
{
    return eH ? static_cast<long long>(eH[v]) * N + eB[v] : static_cast<long long>(v);
}

/// Sorts n_rows rows of vals (row r is vals[row_off[r] .. row_off[r+1])) ascending by hbRowKey.
__global__ __launch_bounds__(HB_SORT_THREADS) void k_hb_row_sort(int n_rows, const int* __restrict__ row_off,
                                                                 int* __restrict__ vals, const int* __restrict__ eH,
                                                                 const int* __restrict__ eB, int N, int shm_cap)
{
    // dynamic shared memory of the launch: shm_cap keys (8 bytes) and values (4 bytes) for the triples (eH given),
    // shm_cap values only for the slot lists (the value is the key)
    extern __shared__ long long s_dyn[];
    const int r = blockIdx.x;
    if (r >= n_rows) return;
    const int off = row_off[r];
    const int L = row_off[r + 1] - off;
    if (L <= 1) return;
    int* v = vals + off;
    const bool with_key = (eH != nullptr);
    long long* s_key = s_dyn;
    int* s_val = with_key ? reinterpret_cast<int*>(s_dyn + shm_cap) : reinterpret_cast<int*>(s_dyn);
    int P = 1;
    while (P < L) P <<= 1;
    const bool shm = (L <= shm_cap);
    if (shm) {
        for (int x = threadIdx.x; x < L; x += blockDim.x) {
            s_val[x] = v[x];
            if (with_key) s_key[x] = hbRowKey(v[x], eH, eB, N);
        }
        __syncthreads();
    }
    auto cmpswap = [&](int a, int b) {
        if (shm) {
            if (with_key) {
                if (s_key[a] > s_key[b]) {
                    const long long tk = s_key[a]; s_key[a] = s_key[b]; s_key[b] = tk;
                    const int tv = s_val[a]; s_val[a] = s_val[b]; s_val[b] = tv;
                }
            } else if (s_val[a] > s_val[b]) {
                const int tv = s_val[a]; s_val[a] = s_val[b]; s_val[b] = tv;
            }
        } else if (hbRowKey(v[a], eH, eB, N) > hbRowKey(v[b], eH, eB, N)) {
            const int tv = v[a]; v[a] = v[b]; v[b] = tv;
        }
    };
    const int half_P = P >> 1;
    for (int k = 2; k <= P; k <<= 1) {
        const int hk = k >> 1;
        for (int idx = threadIdx.x; idx < half_P; idx += blockDim.x) {   // flip: i with its mirror in the k-block
            const int a = (idx / hk) * k + (idx % hk);
            const int b = a ^ (k - 1);
            if (b < L) cmpswap(a, b);
        }
        __syncthreads();
        for (int d = k >> 2; d > 0; d >>= 1) {                           // half-cleaners
            for (int idx = threadIdx.x; idx < half_P; idx += blockDim.x) {
                const int a = (idx / d) * 2 * d + (idx % d);
                const int b = a + d;
                if (b < L) cmpswap(a, b);
            }
            __syncthreads();
        }
    }
    if (shm)
        for (int x = threadIdx.x; x < L; x += blockDim.x) v[x] = s_val[x];
}

/// Claude Generated (Oct 2026): *out = max over cnt[0..n) (one workgroup of HB_SORT_THREADS threads).
__global__ __launch_bounds__(HB_SORT_THREADS) void k_hb_max(int n, const int* __restrict__ cnt, int* __restrict__ out)
{
    __shared__ int s[HB_SORT_THREADS];
    int m = 0;
    for (int x = threadIdx.x; x < n; x += blockDim.x) m = max(m, cnt[x]);
    s[threadIdx.x] = m;
    __syncthreads();
    for (int h = HB_SORT_THREADS / 2; h > 0; h >>= 1) {
        if (static_cast<int>(threadIdx.x) < h) s[threadIdx.x] = max(s[threadIdx.x], s[threadIdx.x + h]);
        __syncthreads();
    }
    if (threadIdx.x == 0) *out = s[0];
}

/// Claude Generated (Oct 2026): port of k_hb_fill_columns; triple t of the final order is element order[t] of the
/// unordered emission. Also writes the gradient-slot count of each triple (k_hbonds slot mode, buildHBSlotTables).
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_hb_fill_columns_rows(
    int n, const int* __restrict__ order,
    const int* __restrict__ e_A, const int* __restrict__ e_H, const int* __restrict__ e_B,
    const int* __restrict__ in_ct, const int* __restrict__ in_par, const int* __restrict__ Z,
    const double* __restrict__ bas, const double* __restrict__ aci, const double* __restrict__ q,
    const double* __restrict__ repz, const int* __restrict__ nbr_off, const int* __restrict__ nbr_idx,
    int* idx_i, int* idx_j, int* idx_k, int* elem_A, int* elem_B,
    double* q_H, double* q_A, double* q_B, double* bas_A, double* bas_B, double* aci_A, double* aci_B,
    double* r_cut, int* case_type, int* acc_parent, double* repz_B, int* nbB_count, int* nbC_count,
    int* slot_count)
{
    const int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    const int src = order[t];
    const int A = e_A[src], H = e_H[src], B = e_B[src];
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
    const int cB = (ct >= 2) ? (nbr_off[B + 1] - nbr_off[B]) : 0;
    int cC = 0;
    if (ct == 3) {
        for (int p = nbr_off[par]; p < nbr_off[par + 1]; ++p)
            if (nbr_idx[p] != B) ++cC;
    }
    nbB_count[t] = cB;
    nbC_count[t] = cC;
    if (slot_count)
        slot_count[t] = 3 + (ct >= 2 ? cB : 0) + ((ct == 3 && par >= 0) ? 1 + cC : 0);
}

/// Claude Generated (Oct 2026): port of k_hb_fill_neighbors plus the owning atom of every gradient slot in the
/// order of buildHBSlotTables (A, H, B, nb(B) for case >= 2, C and nb(C) for case 3 with a parent) and the
/// per-atom slot counts (atom_cnt, atoms outside [0, N) not counted).
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_hb_fill_neighbors_slots(
    int n, int N, const int* __restrict__ idx_i, const int* __restrict__ idx_j, const int* __restrict__ idx_k,
    const int* __restrict__ case_type, const int* __restrict__ acc_parent,
    const int* __restrict__ nbr_off, const int* __restrict__ nbr_idx,
    const int* __restrict__ nbB_offset, const int* __restrict__ nbC_offset, int* nbB_flat, int* nbC_flat,
    const int* __restrict__ slot_off, int* slot_atom, int* atom_cnt)
{
    const int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    const int ct = case_type[t];
    const int B = idx_k[t];
    const int par = acc_parent[t];
    if (ct >= 2) {
        int w = nbB_offset[t];
        for (int p = nbr_off[B]; p < nbr_off[B + 1]; ++p) nbB_flat[w++] = nbr_idx[p];
        if (ct == 3) {
            w = nbC_offset[t];
            for (int p = nbr_off[par]; p < nbr_off[par + 1]; ++p)
                if (nbr_idx[p] != B) nbC_flat[w++] = nbr_idx[p];
        }
    }
    if (!slot_off) return;
    int s = slot_off[t];
    auto put = [&](int a) {
        slot_atom[s++] = a;
        if (a >= 0 && a < N) atomicAdd(atom_cnt + a, 1);
    };
    put(idx_i[t]); put(idx_j[t]); put(B);
    if (ct >= 2)
        for (int p = nbr_off[B]; p < nbr_off[B + 1]; ++p) put(nbr_idx[p]);
    if (ct == 3 && par >= 0) {
        put(par);
        for (int p = nbr_off[par]; p < nbr_off[par + 1]; ++p)
            if (nbr_idx[p] != B) put(nbr_idx[p]);
    }
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
    double        extra_bohr,    // K1 (Oct 2026, Claude Generated): 2*skin of the Verlet list, 0 = legacy
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
    double        extra_bohr,    // K1 (Oct 2026, Claude Generated): 2*skin of the Verlet list, 0 = legacy
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
    acc_add(&cn_raw[i], contrib, GFNFF_FIX_CN);   // WP10: order-independent with d_gfnff_fix (the CN decides the charges)
    acc_add(&cn_raw[j], contrib, GFNFF_FIX_CN);
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

    acc_add(&hb_cn_per_atom[H], contrib, GFNFF_FIX_CN);
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
            const double asum = alpha[gi] + alpha[gj];
            double val;
            if (r_sq >= 36.0 * asum) val = 1.0 / r;   // far field: erf == 1.0, bit-identical
            else                     val = erf((1.0 / sqrt(asum)) * r) / r;
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
#include "../ff_terms.h"   // Bond, Angle, Dihedral, Inversion structs
#include "src/core/energy_calculators/ff_methods/gfnff_par.h"  // covalent_rad_d3

#include <cstdlib>  // std::getenv (CURCUMA_GFNFF_GPU_FREEZE_HBQ opt-out)
#include <thread>
#include <mutex>   // WP10: setGfnffFixModeHip
#include <limits>  // WP10: quiet_NaN

#include "src/core/curcuma_logger.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <cstdio>
#include "src/core/gpu_fallback.h"   // WP9: counted fallbacks of the device repulsion rebuild

/// Claude Generated (Oct 2026, WP9): thread-local slot of the ROCm GFN-FF switches (ff_workspace_hip.h).
const RocmGfnffOptions*& rocmGfnffOptionsSlot()
{
    static thread_local const RocmGfnffOptions* slot = nullptr;
    return slot;
}

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
// GFN-FF covalent radii for the CN and the angle/dihedral distance damping are built at
// init directly from GFNFFParameters::covalent_rad_d3 (see the upload site), so they can
// never drift from the CPU. (Jul 2026) The previous hand-copied s_rcov_d3_bohr[] had
// drifted to the *Grimme D3* radii for Z>=37 (metals: Pt 1.38 vs GFN-FF 1.12 Å), making
// the GPU CN — and thus the metal bond-r0 and angle terms — diverge from the CPU.

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

/// Claude Generated (Oct 2026, WP10): set the device switch d_gfnff_fix (fixed-point step accumulators) on the
/// current device. The value is per process and device, so it is uploaded only when it changes (normally once,
/// in the workspace constructor); a change waits for the device to be idle so no running step sees a mixed
/// mode. Called again before every step entry point, outside graph capture, for processes with several
/// workspaces of different settings.
static void setGfnffFixModeHip(bool on)
{
    static std::mutex mtx;
    static int uploaded_dev = -1, uploaded_val = -1;
    std::lock_guard<std::mutex> lock(mtx);
    int dev = 0;
    checkCuda(hipGetDevice(&dev), "hipGetDevice (fixed-point switch)");
    const int v = on ? 1 : 0;
    if (dev == uploaded_dev && v == uploaded_val) return;
    checkCuda(hipDeviceSynchronize(), "device sync before fixed-point switch");
    checkCuda(hipMemcpyToSymbol(HIP_SYMBOL(d_gfnff_fix), &v, sizeof(int)), "upload d_gfnff_fix");
    uploaded_dev = dev;
    uploaded_val = v;
}

/// Claude Generated (Oct 2026, WP10): in-place fixed-point -> double conversion of n doubles on a stream.
static void fixToDoubleHip(double* buf, int n, double scale, hipStream_t s)
{
    if (!buf || n <= 0) return;
    k_fix_to_double<<<(n + 255) / 256, 256, 0, s>>>(n, buf, 1.0 / scale);
}

/// Claude Generated (Oct 2026, WP10): snapshot copy of a step accumulator: the doubles themselves, or the
/// converted copy of a fixed-point accumulator (the accumulator keeps its integers).
static void snapshotAccHip(bool fixed, double* dst, const double* src, int n, hipStream_t s)
{
    if (fixed) k_fix_copy<<<(n + 255) / 256, 256, 0, s>>>(n, src, dst, 1.0 / GFNFF_FIX_FORCE);
    else       (void)hipMemcpyAsync(dst, src, n * sizeof(double), hipMemcpyDeviceToDevice, s);
}

// Claude Generated (Oct 2026): atomic-free gather tables of the ROCm backend (FP64 atomicAdd is a CAS loop on gfx1200).
/// Directed-edge CSR of the non-bonded repulsion list for k_repulsion_gather (rebuilt with the list).
struct HipPairGatherCSR {
    CudaBuffer<int> off;      ///< [N+1] row offsets
    CudaBuffer<int> partner;  ///< [2n] partner atom per directed edge
    CudaBuffer<int> pairidx;  ///< [2n] pair index into the RepulsionSoA arrays
    int N = 0;
};
/// Claude Generated (Oct 2026): buffers of the device-side H-bond detection (FFWorkspaceHip::detectHBondsOnDevice()).
/// Buffers only grow (ensure); e_* and order are the candidate buffers of `cap` entries (24 bytes each).
struct HipHBDetectDev {
    CudaBuffer<double> bas, aci, q, repz;
    CudaBuffer<char>   ab_ok;
    CudaBuffer<int>    hyd_list, bond_off, bond_idx, hyd_off, hyd_idx, nbr_off, nbr_idx;
    CudaBuffer<int>    counter, ab_counter, ab_i, ab_j;
    CudaBuffer<float>  fx, fy, fz;
    CudaBuffer<float4> hpos;                         ///< [n_hyd] packed HB hydrogen positions (k_hb_pack_hyd)
    CudaBuffer<int>    e_A, e_H, e_B, e_ct, e_par;   ///< [cap] triples in arrival order
    CudaBuffer<int>    order;                        ///< [cap] arrival index of triple t of the final order
    CudaBuffer<int>    row_cnt, row_off, row_fill;   ///< rows (case, A) or A: counts, offsets [R+1], fill cursors
    CudaBuffer<int>    slot_cnt, slot_atom, atom_cnt; ///< gradient slots per triple, owning atom per slot, per atom
    CudaBuffer<int>    scan_tiles, totals;           ///< scan scratch; [5] totals (nb_B, nb_C, slots, longest rows)
    int ab_cap = 0, cap = 0;
    int*   pin = nullptr;                           ///< pinned host stage of the (A, H, B) download
    size_t pin_n = 0;
    HipHBDetectDev() = default;
    HipHBDetectDev(const HipHBDetectDev&) = delete;
    HipHBDetectDev& operator=(const HipHBDetectDev&) = delete;
    ~HipHBDetectDev() { if (pin) (void)hipHostFree(pin); }
    /// Pinned stage of at least m ints (grows with 1/8 headroom); nullptr if no pinned memory could be had.
    int* pinned(size_t m)
    {
        if (m <= pin_n) return pin;
        if (pin) (void)hipHostFree(pin);
        pin = nullptr; pin_n = 0;
        const size_t want = m + m / 8;
        if (hipHostMalloc(reinterpret_cast<void**>(&pin), want * sizeof(int)) != hipSuccess) { pin = nullptr; return nullptr; }
        pin_n = want;
        return pin;
    }
};

/// Claude Generated (Oct 2026): host arrays of one HB list upload (fillHBondSoA, buildHBSlotTables), kept between
/// refreshes so that a refresh neither zero-fills nor page-faults ~100 MB of fresh memory at 582k triples. Every
/// element that is uploaded is written in each fill.
struct HipHBHostStage {
    std::vector<int>    i, j, k, eA, eB, cs, nbB_off, nbB_cnt, acc_parent, nbC_off, nbC_cnt, nbB_flat, nbC_flat;
    std::vector<double> qH, qA, qB, basA, basB, aciA, aciB, rcut, repz;
    std::vector<int>    slot_off, slot_atom, atom_off, atom_slots;
    std::vector<std::vector<int>> hist;
};

/// Claude Generated (Oct 2026): resize a staging array, with 1/8 headroom on growth (the HB list grows and shrinks
/// by about 1 % between refreshes, so it is reallocated rarely).
template <typename V>
static void hbStageResize(V& a, size_t m)
{
    if (a.capacity() < m) a.reserve(m + m / 8);
    a.resize(m);
}

/// Claude Generated (Oct 2026): pinned staging arena for the H-bond list uploads of the single-stream mode.
/// CudaBuffer::upload is a blocking hipMemcpy from pageable memory (~30 us each on the RX 9060 XT, 24 per
/// refresh at 66 atoms); here each array is copied into pinned memory and sent with hipMemcpyAsync on the
/// main stream. put() waits for the stream before it reuses or replaces the block.
struct HipUploadArena {
    bool   enabled = false;
    char*  base = nullptr;
    size_t cap = 0, used = 0;
    HipUploadArena() = default;
    HipUploadArena(const HipUploadArena&) = delete;
    HipUploadArena& operator=(const HipUploadArena&) = delete;
    ~HipUploadArena() { if (base) (void)hipHostFree(base); }
    /// Pinned copy of bytes at src (256-byte aligned); nullptr if no pinned memory could be had.
    void* put(const void* src, size_t bytes, hipStream_t stream)
    {
        const size_t need = (bytes + 255) & ~static_cast<size_t>(255);
        if (used + need > cap) {
            (void)hipStreamSynchronize(stream);   // pending copies read the old block
            if (need > cap) {
                if (base) (void)hipHostFree(base);
                base = nullptr;
                cap = 0;
                const size_t ncap = std::max<size_t>(need, std::max<size_t>(2 * cap, 1 << 16));
                if (hipHostMalloc(reinterpret_cast<void**>(&base), ncap) != hipSuccess) { base = nullptr; return nullptr; }
                cap = ncap;
            }
            used = 0;
        }
        char* dst = base + used;
        std::memcpy(dst, src, bytes);
        used += need;
        return dst;
    }
};

/// Claude Generated (Oct 2026): upload count elements; through the arena (async from pinned memory,
/// capacity with 25 % headroom so a list that grows by a few entries is not reallocated) when it is
/// enabled, else CudaBuffer::upload as before.
template <typename T>
static void stagedUploadHip(CudaBuffer<T>& buf, const T* src, int count, hipStream_t stream, HipUploadArena* A)
{
    if (count <= 0) return;
    if (!A || !A->enabled) { buf.upload(src, count, stream); return; }
    if (buf.n < count) buf.alloc(count + count / 4 + 16);
    void* pinned = A->put(src, static_cast<size_t>(count) * sizeof(T), stream);
    if (!pinned) { buf.upload(src, count, stream); return; }
    checkCuda(hipMemcpyAsync(buf.ptr, pinned, static_cast<size_t>(count) * sizeof(T), hipMemcpyHostToDevice, stream),
              "staged list upload");
}
template <typename T>
static void stagedUploadHip(CudaBuffer<T>& buf, const std::vector<T>& v, hipStream_t stream, HipUploadArena* A)
{
    stagedUploadHip(buf, v.data(), static_cast<int>(v.size()), stream, A);
}

/// Claude Generated (Oct 2026): pinned host stage of the repulsion list refresh (stagedRepulsionUploadHip), kept
/// between refreshes. Grows with 1/8 headroom (the non-bonded list changes by well under 1 % per refresh), so a
/// refresh neither allocates nor page-faults host memory. cnt holds the per-thread atom degrees of the gather CSR.
struct HipRepStage {
    char*  pin = nullptr;
    size_t cap = 0;
    std::vector<int> cnt;
    HipRepStage() = default;
    HipRepStage(const HipRepStage&) = delete;
    HipRepStage& operator=(const HipRepStage&) = delete;
    ~HipRepStage() { if (pin) (void)hipHostFree(pin); }
    /// Pinned block of at least bytes; nullptr if no pinned memory could be had. Callers synchronise the
    /// stream that read the previous contents before calling.
    char* get(size_t bytes)
    {
        if (bytes <= cap) return pin;
        if (pin) (void)hipHostFree(pin);
        pin = nullptr;
        cap = 0;
        const size_t want = bytes + bytes / 8;
        if (hipHostMalloc(reinterpret_cast<void**>(&pin), want) != hipSuccess) { pin = nullptr; return nullptr; }
        cap = want;
        return pin;
    }
};

/// Per-triple gradient slots of the H-bond list for k_hbonds slot mode + k_slot_gather (rebuilt with the list).
struct HipHBSlotTables {
    HipHBHostStage stage;           ///< Oct 2026: host arrays reused across HB list refreshes
    HipUploadArena arena;           ///< Oct 2026: staged uploads (single-stream mode)
    CudaBuffer<int>    slot_off;    ///< [n+1] slot range per triple: A, H, B, nb(B) (case >= 2), C + nb(C) (case 3)
    CudaBuffer<int>    atom_off;    ///< [N+1] row offsets into atom_slots
    CudaBuffer<int>    atom_slots;  ///< [nslots] slot indices grouped by atom
    CudaBuffer<double> slot_g;      ///< [3*nslots] per-slot gradient written by k_hbonds
    CudaBuffer<double> Qs_H, Qs_A, Qs_B;  ///< [n] charge scaling per triple (k_hb_charge_scaling)
    int nslots = 0;
    int N = 0;
    int nq = 0;                     ///< triples covered by Qs_* (0: k_hbonds computes the scaling itself)
    hipStream_t upload_stream = nullptr;   ///< Oct 2026: stream of the staged uploads (set by fillHBondSoA)
};

// ============================================================================
// SoA upload implementations (declared in gfnff_soa.h)
// ============================================================================

// ---------------------------------------------------------------------------
void DispersionSoA::upload(const std::vector<GFNFFDispersion>& v, int natoms, hipStream_t stream)
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

    // Build per-atom CSR adjacency for k_dispersion_gather. Each pair t = (i,j) becomes two
    // directed edges carrying the pair index t and an is_i flag (so the gather picks
    // dc6dcn_ij vs dc6dcn_ji for the owning atom). The dynamic dc6dcn arrays stay pair-indexed.
    csr_N = natoms;
    std::vector<int> deg(natoms, 0);
    for (int k = 0; k < n; ++k) { ++deg[h_i[k]]; ++deg[h_j[k]]; }
    std::vector<int> off(natoms + 1, 0);
    for (int a = 0; a < natoms; ++a) off[a + 1] = off[a] + deg[a];
    const int tot = off[natoms];                       // == 2*n
    std::vector<int> h_part(tot), h_pidx(tot), h_isi(tot);
    std::vector<int> cur(off.begin(), off.end() - 1);
    for (int k = 0; k < n; ++k) {
        const int i = h_i[k], j = h_j[k];
        int pi = cur[i]++; h_part[pi] = j; h_pidx[pi] = k; h_isi[pi] = 1;
        int pj = cur[j]++; h_part[pj] = i; h_pidx[pj] = k; h_isi[pj] = 0;
    }
    csr_off.upload(off, stream);
    csr_partner.upload(h_part, stream);
    csr_pairidx.upload(h_pidx, stream);
    csr_is_i.upload(h_isi, stream);
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
/// Claude Generated (Oct 2026): run f(begin, end) over [0, n) split into `threads` contiguous ranges, the first
/// range on the calling thread. Every caller writes disjoint index ranges, so the result does not depend on the
/// thread count. Below `min_per_thread` items per thread fewer threads are used.
template <typename F>
static void hipHostParallelFor(int n, int threads, F&& f, int min_per_thread = 16384)
{
    if (n <= 0) return;
    int t = std::max(1, std::min(threads, n / std::max(1, min_per_thread)));
    if (t <= 1) { f(0, n); return; }
    std::vector<std::thread> pool;
    pool.reserve(t - 1);
    for (int c = 1; c < t; ++c) {
        const int b = static_cast<int>(static_cast<long long>(n) * c / t);
        const int e = static_cast<int>(static_cast<long long>(n) * (c + 1) / t);
        pool.emplace_back([&f, b, e]() { f(b, e); });
    }
    f(0, static_cast<int>(static_cast<long long>(n) / t));
    for (auto& th : pool) th.join();
}

// ---------------------------------------------------------------------------
/// Claude Generated (Oct 2026): per-triple slot ranges and the per-atom slot CSR for k_hbonds slot mode, built
/// from the host SoA arrays of fillHBondSoA (device order). Slot order per triple must match k_hbonds: A, H, B,
/// then nb(B) if case >= 2, then C and nb(C) if case == 3 and the parent index is >= 0. Slots whose atom is
/// outside [0, natoms) are not gathered. The per-atom slot lists hold the slots in ascending order (k_slot_gather
/// sums them in that order); the counting sort is split over `threads` ranges of slots with one histogram per
/// range, which places every slot exactly where the serial counting sort does.
static void buildHBSlotTables(int n, const std::vector<int>& h_i, const std::vector<int>& h_j, const std::vector<int>& h_k,
                              const std::vector<int>& h_case,
                              const std::vector<int>& nbB_off, const std::vector<int>& nbB_cnt, const std::vector<int>& nbB_flat,
                              const std::vector<int>& acc_parent,
                              const std::vector<int>& nbC_off, const std::vector<int>& nbC_cnt, const std::vector<int>& nbC_flat,
                              int natoms, HipHBSlotTables& T, int threads)
{
    T.N = natoms;
    std::vector<int>& slot_off = T.stage.slot_off;
    hbStageResize(slot_off, n + 1);
    int tot = 0;
    for (int m = 0; m < n; ++m) {
        slot_off[m] = tot;
        tot += 3 + (h_case[m] >= 2 ? nbB_cnt[m] : 0)
                 + ((h_case[m] == 3 && acc_parent[m] >= 0) ? 1 + nbC_cnt[m] : 0);
    }
    slot_off[n] = tot;
    T.nslots = tot;
    std::vector<int>& slot_atom = T.stage.slot_atom;
    hbStageResize(slot_atom, tot);
    hipHostParallelFor(n, threads, [&](int m0, int m1) {
        for (int m = m0; m < m1; ++m) {
            int q = slot_off[m];
            slot_atom[q++] = h_i[m]; slot_atom[q++] = h_j[m]; slot_atom[q++] = h_k[m];
            if (h_case[m] >= 2)
                for (int t = 0; t < nbB_cnt[m]; ++t) slot_atom[q++] = nbB_flat[nbB_off[m] + t];
            if (h_case[m] == 3 && acc_parent[m] >= 0) {
                slot_atom[q++] = acc_parent[m];
                for (int t = 0; t < nbC_cnt[m]; ++t) slot_atom[q++] = nbC_flat[nbC_off[m] + t];
            }
        }
    });
    // Stable counting sort of the slots by atom: range r of slots counts into hist[r], the start of range r in the
    // list of atom a is atom_off[a] + sum of hist[r' < r][a].
    const int nr = std::max(1, std::min(threads, tot / 65536));
    std::vector<std::vector<int>>& hist = T.stage.hist;
    hist.resize(nr);
    for (auto& h : hist) h.assign(natoms, 0);
    auto range_begin = [&](int r) { return static_cast<int>(static_cast<long long>(tot) * r / nr); };
    hipHostParallelFor(nr, nr, [&](int r0, int r1) {
        for (int r = r0; r < r1; ++r) {
            auto& h = hist[r];
            for (int q = range_begin(r); q < range_begin(r + 1); ++q) {
                const int a = slot_atom[q];
                if (a >= 0 && a < natoms) ++h[a];
            }
        }
    }, 1);
    std::vector<int>& atom_off = T.stage.atom_off;
    atom_off.assign(natoms + 1, 0);
    for (int a = 0; a < natoms; ++a) {
        int c = 0;
        for (int r = 0; r < nr; ++r) { const int h = hist[r][a]; hist[r][a] = atom_off[a] + c; c += h; }
        atom_off[a + 1] = atom_off[a] + c;
    }
    std::vector<int>& atom_slots = T.stage.atom_slots;
    hbStageResize(atom_slots, atom_off[natoms]);
    hipHostParallelFor(nr, nr, [&](int r0, int r1) {
        for (int r = r0; r < r1; ++r) {
            auto& cur = hist[r];
            for (int q = range_begin(r); q < range_begin(r + 1); ++q) {
                const int a = slot_atom[q];
                if (a >= 0 && a < natoms) atom_slots[cur[a]++] = q;
            }
        }
    }, 1);
    HipUploadArena* A = &T.arena;
    stagedUploadHip(T.slot_off, slot_off.data(), n + 1, T.upload_stream, A);
    stagedUploadHip(T.atom_off, atom_off.data(), natoms + 1, T.upload_stream, A);
    if (atom_off[natoms] > 0) stagedUploadHip(T.atom_slots, atom_slots.data(), atom_off[natoms], T.upload_stream, A);
    if (T.slot_g.n < 3 * tot) T.slot_g.alloc(3 * (tot + tot / 8));   // headroom: the list grows and shrinks between refreshes
}

// Claude Generated (Oct 2026): body of HBondSoA::upload with an optional permutation: device triple m is
// v[perm[m]] (perm == nullptr: identity). Used to group the triples by case type (fillHBondSoA callers).
static void fillHBondSoA(HBondSoA& s, const std::vector<GFNFFHydrogenBond>& vin, const int* perm,
                         const std::vector<int>& atom_types, hipStream_t stream,
                         HipHBSlotTables* T = nullptr, int natoms = 0, int threads = 1);

void HBondSoA::upload(const std::vector<GFNFFHydrogenBond>& v,
                       const std::vector<int>& atom_types,
                       hipStream_t stream)
{
    fillHBondSoA(*this, v, nullptr, atom_types, stream);
}

// Oct 2026: the arrays are written directly in device order (no separate permutation pass) by `threads` host
// threads over disjoint triple ranges; the neighbour lists nb_B_flat / nb_C_flat are stored in device order too
// (each triple carries its own offset and count, so k_hbonds reads the same values).
static void fillHBondSoA(HBondSoA& s, const std::vector<GFNFFHydrogenBond>& vin, const int* perm,
                         const std::vector<int>& atom_types, hipStream_t stream,
                         HipHBSlotTables* T, int natoms, int threads)
{
    int& n = s.n;
    n = static_cast<int>(vin.size());
    if (T) { T->nslots = 0; T->N = natoms; }   // Oct 2026: an empty list must not gather stale slots
    if (n == 0) return;
    const std::vector<GFNFFHydrogenBond>& v = vin;

    const int nat = static_cast<int>(atom_types.size());
    HipHBHostStage local_stage;
    HipHBHostStage& st = T ? T->stage : local_stage;
    for (auto* a : { &st.i, &st.j, &st.k, &st.eA, &st.eB, &st.cs, &st.nbB_cnt, &st.acc_parent, &st.nbC_cnt })
        hbStageResize(*a, n);
    for (auto* a : { &st.qH, &st.qA, &st.qB, &st.basA, &st.basB, &st.aciA, &st.aciB, &st.rcut, &st.repz })
        hbStageResize(*a, n);
    hbStageResize(st.nbB_off, n + 1);
    hbStageResize(st.nbC_off, n + 1);
    std::vector<int>    &h_i = st.i, &h_j = st.j, &h_k = st.k;
    std::vector<int>    &h_eA = st.eA, &h_eB = st.eB, &h_case = st.cs;
    std::vector<double> &h_qH = st.qH, &h_qA = st.qA, &h_qB = st.qB;
    std::vector<double> &h_basA = st.basA, &h_basB = st.basB, &h_aciA = st.aciA, &h_aciB = st.aciB;
    std::vector<double> &h_rcut = st.rcut;
    std::vector<int>    &h_nbB_offset = st.nbB_off, &h_nbB_count = st.nbB_cnt;
    std::vector<int>    &h_acc_parent = st.acc_parent;
    std::vector<int>    &h_nbC_offset = st.nbC_off, &h_nbC_count = st.nbC_cnt;
    std::vector<double> &h_repz = st.repz;

    // Pass 1: per-triple values in device order
    hipHostParallelFor(n, threads, [&](int m0, int m1) {
        for (int m = m0; m < m1; ++m) {
            const GFNFFHydrogenBond& hb = v[perm ? perm[m] : m];
            h_i[m]     = hb.i;   // donor A
            h_j[m]     = hb.j;   // hydrogen H
            h_k[m]     = hb.k;   // acceptor B
            h_eA[m]    = (hb.i < nat) ? atom_types[hb.i] : 0;
            h_eB[m]    = (hb.k < nat) ? atom_types[hb.k] : 0;
            h_qH[m]    = hb.q_H;
            h_qA[m]    = hb.q_A;
            h_qB[m]    = hb.q_B;
            h_basA[m]  = hb.basicity_A;
            h_basB[m]  = hb.basicity_B;
            h_aciA[m]  = hb.acidity_A;
            h_aciB[m]  = hb.acidity_B;
            h_rcut[m]  = hb.r_cut;
            h_case[m]  = hb.case_type;
            h_acc_parent[m] = hb.acceptor_parent_index;
            // repz_B for case 4
            const int zB = (hb.k < nat) ? atom_types[hb.k] : 0;
            h_repz[m] = (zB >= 1 && zB <= 86) ? GFNFFParameters::repz[zB - 1] : 1.0;
            h_nbB_count[m] = static_cast<int>(hb.neighbors_B.size());
            h_nbC_count[m] = static_cast<int>(hb.neighbors_C.size());
        }
    });
    // Offsets of the neighbour lists (device order)
    h_nbB_offset[0] = 0; h_nbC_offset[0] = 0;
    for (int m = 0; m < n; ++m) {
        h_nbB_offset[m + 1] = h_nbB_offset[m] + h_nbB_count[m];
        h_nbC_offset[m + 1] = h_nbC_offset[m] + h_nbC_count[m];
    }
    const int total_B = h_nbB_offset[n], total_C = h_nbC_offset[n];
    std::vector<int>& h_nbB_flat = st.nbB_flat;
    std::vector<int>& h_nbC_flat = st.nbC_flat;
    hbStageResize(h_nbB_flat, total_B);
    hbStageResize(h_nbC_flat, total_C);
    // Pass 2: neighbour lists
    hipHostParallelFor(n, threads, [&](int m0, int m1) {
        for (int m = m0; m < m1; ++m) {
            const GFNFFHydrogenBond& hb = v[perm ? perm[m] : m];
            int ob = h_nbB_offset[m];
            for (int nb : hb.neighbors_B) h_nbB_flat[ob++] = nb;
            int oc = h_nbC_offset[m];
            for (int nb : hb.neighbors_C) h_nbC_flat[oc++] = nb;
        }
    });

    s.total_nb_B = total_B;
    s.total_nb_C = total_C;

    // Oct 2026 (Claude Generated): staged uploads in the single-stream mode (T->arena.enabled), see
    // HipUploadArena; otherwise the blocking CudaBuffer::upload as before. The stage vectors are resized
    // with headroom, so the element counts are passed explicitly.
    HipUploadArena* A = T ? &T->arena : nullptr;
    if (T) T->upload_stream = stream;
    stagedUploadHip(s.idx_i, h_i.data(), n, stream, A);  stagedUploadHip(s.idx_j, h_j.data(), n, stream, A);
    stagedUploadHip(s.idx_k, h_k.data(), n, stream, A);
    stagedUploadHip(s.elem_A, h_eA.data(), n, stream, A);  stagedUploadHip(s.elem_B, h_eB.data(), n, stream, A);
    stagedUploadHip(s.q_H, h_qH.data(), n, stream, A);  stagedUploadHip(s.q_A, h_qA.data(), n, stream, A);
    stagedUploadHip(s.q_B, h_qB.data(), n, stream, A);
    stagedUploadHip(s.basicity_A, h_basA.data(), n, stream, A);  stagedUploadHip(s.basicity_B, h_basB.data(), n, stream, A);
    stagedUploadHip(s.acidity_A, h_aciA.data(), n, stream, A);   stagedUploadHip(s.acidity_B, h_aciB.data(), n, stream, A);
    stagedUploadHip(s.r_cut, h_rcut.data(), n, stream, A);
    stagedUploadHip(s.case_type, h_case.data(), n, stream, A);
    stagedUploadHip(s.nb_B_offset, h_nbB_offset.data(), n, stream, A);
    stagedUploadHip(s.nb_B_count, h_nbB_count.data(), n, stream, A);
    if (total_B > 0) stagedUploadHip(s.nb_B_flat, h_nbB_flat.data(), total_B, stream, A);
    stagedUploadHip(s.acceptor_parent, h_acc_parent.data(), n, stream, A);
    stagedUploadHip(s.nb_C_offset, h_nbC_offset.data(), n, stream, A);
    stagedUploadHip(s.nb_C_count, h_nbC_count.data(), n, stream, A);
    if (total_C > 0) stagedUploadHip(s.nb_C_flat, h_nbC_flat.data(), total_C, stream, A);
    stagedUploadHip(s.repz_B, h_repz.data(), n, stream, A);

    if (T) buildHBSlotTables(n, h_i, h_j, h_k, h_case, h_nbB_offset, h_nbB_count, h_nbB_flat,
                             h_acc_parent, h_nbC_offset, h_nbC_count, h_nbC_flat, natoms, *T, threads);
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
    // Claude Generated (Sep 2026): implicit all-pairs Coulomb (no pair list), the HIP counterpart of
    // the CUDA path. k_coulomb_tiles_hip enumerates every pair from the per-atom alpeeq (d_coul_alp).
    bool           coulomb_implicit = false;
    double         coulomb_rcut = 100.0;
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

    // Claude Generated (Oct 2026): atomic-free gather tables (FP64 atomicAdd is a CAS loop on gfx1200).
    using PairGatherCSR = HipPairGatherCSR;
    using HBSlotTables  = HipHBSlotTables;
    PairGatherCSR nb_rep_csr;
    HBSlotTables hb_slots;
    HipHBDetectDev hbdet;            ///< Claude Generated (Oct 2026): device-side H-bond detection buffers
    HipRepStage    rep_stage;        ///< Claude Generated (Oct 2026): pinned stage of the repulsion list refresh
    bool rep_gather_enabled = true;  ///< CURCUMA_ROCM_REP_GATHER=0 restores the pair kernel k_repulsion (non-bonded list)
    // WP9-A (Claude Generated, Oct 2026): tables and scratch of the device repulsion rebuild
    // (generateNonbondedRepulsionOnGPU); topology-constant, uploaded at every rebuild (small).
    CudaBuffer<int>       rep_valid, rep_metal, rep_hh_dist, rep_lo, rep_hi, rep_rowoff;
    CudaBuffer<double>    rep_dum, rep_repz;
    CudaBuffer<long long> rep_bond_keys, rep_hh_keys;
    long rep_device_builds = 0;      ///< device rebuilds done (verbosity-2 summary at destruction)
    long rep_host_fallbacks = 0;     ///< device rebuilds that failed and were replaced by the host build
    bool rep_device_off = false;     ///< set after a failure: later rebuilds take the host build directly
    bool hb_gather_enabled  = true;  ///< CURCUMA_ROCM_HB_GATHER=0 restores k_hbonds with atomicAdd
    /// WP10 (Claude Generated, Oct 2026): fixed-point step accumulators (-gfnff.rocm_fixed_point, see fix_add);
    /// d_grad, d_dEdcn, d_energies, d_zz_hb, d_cn_raw and d_hb_cn_per_atom hold integers until k_fix_to_double.
    bool fixed_point = false;
    int* h_fix_bad = nullptr;        ///< pinned read-back of d_gfnff_fix_bad (fixed_point only)
    long fix_bad_steps = 0;          ///< steps whose fixed-point accumulation saw a non-finite / out-of-range value

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
    RefCoordSoA        cn_ref_coords;     ///< K1 (Oct 2026, Claude Generated): geometry of the last CN pair-list build
    CudaBuffer<int>    d_cn_disp_flag;    ///< K1 (Oct 2026, Claude Generated): [1] flag for cnPairListStale()
    // Claude Generated (Oct 2026): lean per-step bookkeeping (CURCUMA_ROCM_LEAN_STEP=0 restores the former
    // sequence): both displacement tests in one launch with one pinned read-back, no dc6dcn skip test
    // (its flag has no reader), CN host copies only when finalizeCNForCPU() asks, one zeroing kernel.
    bool               lean_step = true;
    CudaBuffer<int>    d_disp_flags2;     ///< [2] flags of k_check_displacement_pair
    int*               h_disp_flags = nullptr;   ///< pinned [2] read-back of d_disp_flags2
    double             cn_stale_thr = -1.0;      ///< threshold of the last cnPairListStale() call (Bohr)
    /// Claude Generated (Oct 2026): all per-step kernels on the main stream (no worker streams, no
    /// cross-stream events) below kRocmSingleStreamMaxAtoms; CURCUMA_ROCM_SINGLE_STREAM=1/0 forces it.
    bool               single_stream = false;
    /// Claude Generated (Oct 2026): host mirrors of the two displacement references (x, y, z blocks of N),
    /// copied from the pinned geometry that setGeometry() uploaded. With single_stream the displacement
    /// tests run on the host (same per-atom test), which saves a device round trip per step.
    std::vector<double> h_ref_xyz, h_cnref_xyz;
    bool               h_geom_ok = false, h_ref_ok = false, h_cnref_ok = false;
    /// Claude Generated (Oct 2026): inputs of the last HB/XB list uploads (single_stream mode only). The
    /// forced HB re-detection every hb_update_force_every steps mostly finds the same lists; an identical
    /// input leaves the device arrays as they are instead of ~30 blocking uploads.
    std::vector<GFNFFHydrogenBond> prev_hb_input;
    bool               prev_hb_valid = false;
    std::vector<int>    prev_alpha_h, prev_alpha_b;
    std::vector<double> prev_alpha_rcov;
    bool               prev_alpha_valid = false;
    std::vector<int>    prev_nr_hb, prev_hb_H;
    bool               prev_bond_hb_valid = false;
    bool               prev_xb_valid = false;
    /// Claude Generated (Oct 2026): CN host copies for finalizeCNForCPU() (single-stream mode): event after
    /// k_dlogdcn, copies on the idle stream_threebody, so the host EEQ waits for the CN only, not for
    /// Phase 1. No stream of its own: one more stream doubled the 66-atom step (0.59 -> 1.17 ms), the
    /// hardware queues are then shared differently.
    hipEvent_t         event_cn_ready = nullptr;
    bool               cn_event_recorded = false;
    hipEvent_t         event_eeq_order = nullptr;   ///< Oct 2026: main -> auxiliary stream order of the small EEQ
    /// Claude Generated (Oct 2026): K2 on the single-stream path - the gradient-mode dispersion gather
    /// was launched by launchDispersionGradientEarly() this step (Phase 2 skips it).
    bool               disp_launched_early = false;
    double*            h_charges_stage = nullptr;   ///< pinned [N] staging for the host-EEQ charge upload
    int                cn_stale_cached = -1;     ///< flag computed ahead by checkDisplacement(), -1 = none

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

    // Claude Generated (Oct 2026): factorised D4 C6 / dC6dCN sums (k_disp_uv, k_dc6dcn_per_pair_uv).
    int                disp_ntypes = 0;        ///< element types present (compact index)
    CudaBuffer<int>    d_atom_ctype;           ///< [N] compact element type per atom
    CudaBuffer<int>    d_ctype_elem;           ///< [disp_ntypes] 0-based element per compact type
    CudaBuffer<double> d_disp_uv;              ///< [N * disp_ntypes * 2*D4_MAX_REF] u | v per (atom, type)
    bool disp_uv_valid        = false; ///< d_disp_uv matches d_gw/d_dgw (cleared when the pair list is replaced)
    bool disp_pair_c6_current = true;  ///< disp.C6 / dc6dcn_* hold the values of the last d_disp_uv
    bool c6_table_symmetric   = false; ///< C6ref(Zi,Zj,ri,rj) == C6ref(Zj,Zi,rj,ri) bitwise (fused gather allowed)
    bool disp_fused_enabled   = true;  ///< CURCUMA_ROCM_DISP_FUSED=0 switches the fused gather off
    // WP9-B (Claude Generated, Oct 2026): tabulated dc6dcn (-gfnff.rocm_dc6dcn_tables, CURCUMA_DC6DCN_TABLES
    // overrides): per-pair C6/dc6dcn from k_d4_partner_tables + k_dc6dcn_per_pair_tab; the fused gather is off.
    bool dc6_tab_enabled      = false;
    CudaBuffer<double> d_tabA, d_tabW;          ///< [N * disp_ntypes * D4_MAX_REF] each

    // WP-A (Jun 2026): on-device D4 dispersion pair-list build buffers
    CudaBuffer<double> d_sqrtzr4r2;        ///< [118] sqrt(Z·r4/r2) per element
    CudaBuffer<double> d_topo_q;           ///< [N] Phase-1 topology charges (for zetac6)
    CudaBuffer<double> d_zeta_zeff;        ///< [86] zeta zeff per element
    CudaBuffer<double> d_zeta_c;           ///< [86] zeta c per element
    CudaBuffer<int>    d_disp_pair_counter; ///< [1] atomic counter for the two-pass build
    bool               disp_built_on_gpu = false; ///< true when disp SoA was device-built

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

/// Claude Generated (Oct 2026): directed-edge CSR of a repulsion pair list for k_repulsion_gather.
/// Pair k = (i, j) becomes the edges i -> j and j -> i, both carrying k; rows keep the list order.
/// Pairs with an atom index outside [0, natoms) are left out (the pair kernel would read out of bounds).
static void buildRepulsionGatherCSR(const std::vector<GFNFFRepulsion>& v, int natoms,
                                    FFWorkspaceHipImpl::PairGatherCSR& c)
{
    c.N = natoms;
    const int n = static_cast<int>(v.size());
    auto ok = [natoms](int a) { return a >= 0 && a < natoms; };
    std::vector<int> off(natoms + 1, 0);
    for (int k = 0; k < n; ++k)
        if (ok(v[k].i) && ok(v[k].j)) { ++off[v[k].i + 1]; ++off[v[k].j + 1]; }
    for (int a = 0; a < natoms; ++a) off[a + 1] += off[a];
    const int tot = off[natoms];
    std::vector<int> part(tot), pidx(tot);
    std::vector<int> cur(off.begin(), off.end() - 1);
    for (int k = 0; k < n; ++k) {
        const int i = v[k].i, j = v[k].j;
        if (!ok(i) || !ok(j)) continue;
        int pi = cur[i]++; part[pi] = j; pidx[pi] = k;
        int pj = cur[j]++; part[pj] = i; pidx[pj] = k;
    }
    c.off.upload(off);
    if (tot > 0) { c.partner.upload(part); c.pairidx.upload(pidx); }
}

/// Claude Generated (Oct 2026): repulsion list refresh through the pinned stage impl.rep_stage; the device receives
/// the same bytes as RepulsionSoA::upload() for both lists plus buildRepulsionGatherCSR() for the non-bonded one
/// (when gather is true). The non-bonded list is packed by up to `threads` host threads over contiguous ranges of
/// pairs; each range counts its own atom degrees, and its directed edges are placed after those of the lower ranges,
/// so every CSR row keeps the ascending pair order of the serial build. Copies are hipMemcpyAsync on impl.stream (the
/// SoA columns go out while the CSR is scattered); returns after the stream is synchronised, since the repulsion
/// kernels run on another stream. Device buffers grow with 1/8 headroom. urea400-water1000 (1.47 M pairs): the
/// serial path spent ~25 ms per refresh, mostly in fresh host vectors and the serial CSR build, ~2 ms in the copies.
/// @return false if no pinned memory could be had (nothing was uploaded; the caller takes the serial path)
static bool stagedRepulsionUploadHip(FFWorkspaceHipImpl& impl, const std::vector<GFNFFRepulsion>& bonded,
                                     const std::vector<GFNFFRepulsion>& nonbonded, bool gather, int natoms,
                                     int threads)
{
    auto al = [](size_t b) { return (b + 255) & ~static_cast<size_t>(255); };
    auto soaBytes = [&al](size_t m) { return 3 * al(m * sizeof(double)) + 2 * al(m * sizeof(int)); };
    const size_t nb = bonded.size(), nn = nonbonded.size();
    const size_t bytes = soaBytes(nb) + soaBytes(nn)
                       + (gather ? al((natoms + 1) * sizeof(int)) + 2 * al(2 * nn * sizeof(int)) : 0);
    char* base = impl.rep_stage.get(bytes);
    if (!base) return false;

    struct Cols { double *alpha, *repab, *rcut; int *ii, *jj; };
    size_t used = 0;
    auto carve = [&](size_t b) { char* p = base + used; used += al(b); return p; };
    auto cols = [&](size_t m) {
        Cols c;
        c.alpha = reinterpret_cast<double*>(carve(m * sizeof(double)));
        c.repab = reinterpret_cast<double*>(carve(m * sizeof(double)));
        c.rcut  = reinterpret_cast<double*>(carve(m * sizeof(double)));
        c.ii    = reinterpret_cast<int*>(carve(m * sizeof(int)));
        c.jj    = reinterpret_cast<int*>(carve(m * sizeof(int)));
        return c;
    };
    const Cols cb = cols(nb), cn = cols(nn);
    auto pack = [](const std::vector<GFNFFRepulsion>& v, const Cols& c, size_t a, size_t e) {
        for (size_t k = a; k < e; ++k) {
            c.alpha[k] = v[k].alpha;
            c.repab[k] = v[k].repab;
            c.rcut[k]  = v[k].r_cut;
            c.ii[k]    = v[k].i;
            c.jj[k]    = v[k].j;
        }
    };
    auto ok = [natoms](int a) { return a >= 0 && a < natoms; };

    // Ranges of the non-bonded list: identical in both passes (the CSR placement depends on them).
    const int T = (nn < 32768) ? 1 : std::max(1, std::min(threads, static_cast<int>(nn / 16384)));
    auto rangeOf = [nn, T](int t, size_t& a, size_t& e) { a = nn * t / T; e = nn * (t + 1) / T; };
    auto runT = [T](auto&& fn) {
        if (T == 1) { fn(0); return; }
        std::vector<std::thread> pool;
        pool.reserve(T - 1);
        for (int t = 1; t < T; ++t) pool.emplace_back([&fn, t]() { fn(t); });
        fn(0);
        for (auto& th : pool) th.join();
    };
    std::vector<int>& cnt = impl.rep_stage.cnt;
    if (gather) cnt.assign(static_cast<size_t>(T) * natoms, 0);

    // Pass 1: SoA columns and, for the gather CSR, per-range atom degrees.
    pack(bonded, cb, 0, nb);
    runT([&](int t) {
        size_t a, e;
        rangeOf(t, a, e);
        pack(nonbonded, cn, a, e);
        if (!gather) return;
        int* deg = cnt.data() + static_cast<size_t>(t) * natoms;
        for (size_t k = a; k < e; ++k) {
            const int i = cn.ii[k], j = cn.jj[k];
            if (ok(i) && ok(j)) { ++deg[i]; ++deg[j]; }
        }
    });

    hipStream_t s = impl.stream;
    auto put = [s](auto& buf, const auto* src, size_t m) {
        if (m == 0) return;
        const int mi = static_cast<int>(m);
        if (buf.n < mi) buf.alloc(mi + mi / 8 + 16);
        checkCuda(hipMemcpyAsync(buf.ptr, src, m * sizeof(*src), hipMemcpyHostToDevice, s), "staged repulsion upload");
    };
    auto putSoA = [&put](RepulsionSoA& soa, const Cols& c, size_t m) {
        soa.n = static_cast<int>(m);
        put(soa.idx_i, c.ii, m);
        put(soa.idx_j, c.jj, m);
        put(soa.alpha, c.alpha, m);
        put(soa.repab, c.repab, m);
        put(soa.r_cut, c.rcut, m);
    };
    putSoA(impl.bonded_rep, cb, nb);
    putSoA(impl.nonbonded_rep, cn, nn);

    // Pass 2: gather CSR. Row offsets from the summed degrees, then the start of range t in row a is the row
    // offset plus the degrees of the lower ranges in that row.
    if (gather) {
        int* off = reinterpret_cast<int*>(carve((natoms + 1) * sizeof(int)));
        int* part = reinterpret_cast<int*>(carve(2 * nn * sizeof(int)));
        int* pidx = reinterpret_cast<int*>(carve(2 * nn * sizeof(int)));
        off[0] = 0;
        for (int a = 0; a < natoms; ++a) {
            int run = off[a];
            for (int t = 0; t < T; ++t) {
                int& c = cnt[static_cast<size_t>(t) * natoms + a];
                const int d = c;
                c = run;   // becomes the write cursor of range t in row a
                run += d;
            }
            off[a + 1] = run;
        }
        const int tot = off[natoms];
        runT([&](int t) {
            size_t a, e;
            rangeOf(t, a, e);
            int* cur = cnt.data() + static_cast<size_t>(t) * natoms;
            for (size_t k = a; k < e; ++k) {
                const int i = cn.ii[k], j = cn.jj[k];
                if (!ok(i) || !ok(j)) continue;
                const int pi = cur[i]++; part[pi] = j; pidx[pi] = static_cast<int>(k);
                const int pj = cur[j]++; part[pj] = i; pidx[pj] = static_cast<int>(k);
            }
        });
        FFWorkspaceHipImpl::PairGatherCSR& c = impl.nb_rep_csr;
        c.N = natoms;
        put(c.off, off, static_cast<size_t>(natoms) + 1);
        put(c.partner, part, static_cast<size_t>(tot));
        put(c.pairidx, pidx, static_cast<size_t>(tot));
    }
    checkCuda(hipStreamSynchronize(s), "staged repulsion upload sync");
    return true;
}

/// Claude Generated (Oct 2026): device order of the H-bond triples, grouped by case type (stable within a case).
/// k_hbonds branches on the case (case 1 gradient, neighbour loops for case >= 2, angle and torsions for case 3);
/// in list order the cases alternate every ~4 triples (urea400-water1000), so most wavefronts ran every branch.
static std::vector<int> hbCasePermutation(const std::vector<GFNFFHydrogenBond>& v)
{
    // Counting sort on the case type (one sequential pass to count, one to place): stable, O(n). A comparison
    // sort over the ~180-byte structs cost tens of ms per HB list refresh at 582k triples.
    const int n = static_cast<int>(v.size());
    std::vector<int> ct(n);
    int ct_max = 0;
    for (int m = 0; m < n; ++m) { ct[m] = std::max(0, v[m].case_type); ct_max = std::max(ct_max, ct[m]); }
    std::vector<int> start(ct_max + 2, 0);
    for (int m = 0; m < n; ++m) ++start[ct[m] + 1];
    for (int c = 0; c <= ct_max; ++c) start[c + 1] += start[c];
    std::vector<int> perm(n);
    for (int m = 0; m < n; ++m) perm[start[ct[m]]++] = m;
    return perm;
}

/// Claude Generated (Oct 2026): fill the per-triple charge scaling factors (k_hb_charge_scaling) for the current
/// HB charges; synchronises impl.stream so that the next k_hbonds launch (stream C) sees them.
static void updateHBChargeScaling(FFWorkspaceHipImpl& impl)
{
    auto& T = impl.hb_slots;
    T.nq = 0;
    if (!impl.hb_gather_enabled || impl.hbonds.n <= 0)
        return;
    const int n = impl.hbonds.n;
    T.Qs_H.ensure(n); T.Qs_A.ensure(n); T.Qs_B.ensure(n);
    k_hb_charge_scaling<<<(n + 255) / 256, 256, 0, impl.stream>>>(
        n, impl.hbonds.q_H.ptr, impl.hbonds.q_A.ptr, impl.hbonds.q_B.ptr,
        T.Qs_H.ptr, T.Qs_A.ptr, T.Qs_B.ptr);
    checkCuda(hipGetLastError(), "k_hb_charge_scaling launch");
    // Oct 2026 (Claude Generated): k_hbonds reads Qs_* later on the same stream; the single-stream mode
    // leaves out this host wait.
    if (!impl.single_stream)
        checkCuda(hipStreamSynchronize(impl.stream), "k_hb_charge_scaling sync");
    T.nq = n;
}

// ============================================================================
// FFWorkspaceHip Constructor
// ============================================================================

FFWorkspaceHip::FFWorkspaceHip(const GFNFFParameterSet& params,
                                int natoms,
                                const std::vector<int>& atom_types)
    : m_natoms(natoms)
{
    // GPU-wedge hardening: reset the device on Ctrl-C / kill / abort (see top of file).
    installHipSignalHandlers();

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
    checkCuda(hipEventCreateWithFlags(&m_impl->event_cn_ready, hipEventDisableTiming), "hipEventCreate cn_ready");
    checkCuda(hipEventCreateWithFlags(&m_impl->event_eeq_order, hipEventDisableTiming), "hipEventCreate eeq_order");
    checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_impl->h_charges_stage), std::max(1, natoms) * sizeof(double)),
              "pinned alloc h_charges_stage");
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

    // Upload covalent radii to constant memory (used by the GPU CN and the angle/dihedral
    // distance damping). Taken directly from GFNFFParameters::covalent_rad_d3 (already in
    // Bohr) — the same table the CPU CN and damping use — so GPU and CPU CN agree for metals.
    upload_rcov_d3(GFNFFParameters::covalent_rad_d3.data(), 87);

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
    m_impl->disp.upload(params.dispersions, natoms, stream);
    m_impl->disp_uv_valid = false; m_impl->disp_pair_c6_current = true;  // Oct 2026: new pair list, d_disp_uv stale
    m_impl->bonded_rep.upload(params.bonded_repulsions, stream);
    m_impl->nonbonded_rep.upload(params.nonbonded_repulsions, stream);
    // Claude Generated (Oct 2026): atomic-free gathers for the non-bonded repulsion and the H bonds.
    {
        const char* env_rep = std::getenv("CURCUMA_ROCM_REP_GATHER");
        const char* env_hb  = std::getenv("CURCUMA_ROCM_HB_GATHER");
        m_impl->rep_gather_enabled = !(env_rep && std::string(env_rep) == "0");
        m_impl->hb_gather_enabled  = !(env_hb  && std::string(env_hb)  == "0");
    }
    {   // WP10 (Claude Generated, Oct 2026): fixed-point step accumulators, -gfnff.rocm_fixed_point (method options)
        const RocmGfnffOptions* opt = rocmGfnffOptionsSlot();
        m_impl->fixed_point = opt ? opt->fixed_point : RocmGfnffOptions{}.fixed_point;
        setGfnffFixModeHip(m_impl->fixed_point);
        if (m_impl->fixed_point) {
            checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_impl->h_fix_bad), sizeof(int)), "pinned alloc h_fix_bad");
            *m_impl->h_fix_bad = 0;
        }
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info(fmt::format("ROCm GFN-FF accumulation: {}", m_impl->fixed_point
                ? "fixed point (order-independent, -gfnff.rocm_fixed_point true)"
                : "FP64 atomics (-gfnff.rocm_fixed_point false)"));
    }
    if (m_impl->rep_gather_enabled)
        buildRepulsionGatherCSR(params.nonbonded_repulsions, natoms, m_impl->nb_rep_csr);
    m_impl->coulomb.upload(params.coulombs, natoms, stream);
    if (params.coulomb_implicit && params.coul_self_alp.size() == natoms) {
        m_impl->coulomb_implicit = true;
        m_impl->coulomb_rcut = params.coulomb_implicit_rcut;
    }
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
    if (m_impl->hb_gather_enabled) {
        // Claude Generated (Oct 2026): triples grouped by case type, slot tables in the same order
        const std::vector<int> perm = hbCasePermutation(params.hbonds);
        fillHBondSoA(m_impl->hbonds, params.hbonds, perm.data(), atom_types, stream, &m_impl->hb_slots, natoms);
    } else {
        m_impl->hbonds.upload(params.hbonds, atom_types, stream);
    }
    updateHBChargeScaling(*m_impl);   // Claude Generated (Oct 2026)
    if (CurcumaLogger::get_verbosity() >= 3)   // Oct 2026: debug mirror only where it is read (see updateHBonds)
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
    if (m_impl->h_dc6dcn_skip) *m_impl->h_dc6dcn_skip = 0;
    {   // Claude Generated (Oct 2026): lean per-step bookkeeping, see FFWorkspaceHipImpl::lean_step
        const char* env = std::getenv("CURCUMA_ROCM_LEAN_STEP");
        m_impl->lean_step = !(env && std::string(env) == "0");
        m_impl->d_disp_flags2.alloc(2);
        if (!m_impl->h_disp_flags)
            checkCuda(hipHostMalloc(reinterpret_cast<void**>(&m_impl->h_disp_flags), 2 * sizeof(int)),
                      "pinned alloc h_disp_flags");
        m_impl->h_disp_flags[0] = m_impl->h_disp_flags[1] = 0;
        const char* env_ss = std::getenv("CURCUMA_ROCM_SINGLE_STREAM");
        if (env_ss && *env_ss) m_impl->single_stream = (*env_ss != '0');
        else                   m_impl->single_stream = m_impl->lean_step && natoms < kRocmSingleStreamMaxAtoms;
        m_impl->hb_slots.arena.enabled = m_impl->single_stream;   // staged HB list uploads (HipUploadArena)
    }
    // Upload atom types once (static during simulation)
    std::vector<int> h_atom_types = atom_types;  // copy
    m_impl->d_atom_types.upload(h_atom_types, stream);

    // Claude Generated (Oct 2026): compact element types for the factorised D4 C6 sums (k_disp_uv).
    {
        std::vector<int> type_of_elem(D4_MAX_ELEM, -1), ctype_elem, atom_ctype(natoms, 0);
        for (int a = 0; a < natoms; ++a) {
            const int e = atom_types[a] - 1;
            if (e < 0 || e >= D4_MAX_ELEM) continue;      // k_disp_uv writes zeros for such atoms
            if (type_of_elem[e] < 0) { type_of_elem[e] = static_cast<int>(ctype_elem.size()); ctype_elem.push_back(e); }
            atom_ctype[a] = type_of_elem[e];
        }
        if (ctype_elem.empty()) ctype_elem.push_back(-1);
        m_impl->disp_ntypes = static_cast<int>(ctype_elem.size());
        m_impl->d_ctype_elem.upload(ctype_elem, stream);
        m_impl->d_atom_ctype.upload(atom_ctype, stream);
        m_impl->d_disp_uv.alloc(natoms * m_impl->disp_ntypes * 2 * D4_MAX_REF);
        // WP9-B (Claude Generated, Oct 2026): tabulated dc6dcn. Switch from the method options, the CUDA
        // environment name CURCUMA_DC6DCN_TABLES (0/1) overrides it; CURCUMA_ROCM_FAIL_DC6TAB=1 makes the table
        // allocation fail (fault injection). A failed allocation keeps the k_disp_uv path, counted as a fallback.
        {
            const RocmGfnffOptions* opt = rocmGfnffOptionsSlot();
            bool tab = opt ? opt->dc6dcn_tables : RocmGfnffOptions{}.dc6dcn_tables;
            const char* e = std::getenv("CURCUMA_DC6DCN_TABLES");
            if (e && *e) tab = std::string(e) != "0";
            if (tab) {
                try {
                    const char* f = std::getenv("CURCUMA_ROCM_FAIL_DC6TAB");
                    if (f && *f && std::string(f) != "0")
                        throw std::runtime_error("fault injection CURCUMA_ROCM_FAIL_DC6TAB");
                    const int sz = natoms * m_impl->disp_ntypes * D4_MAX_REF;
                    m_impl->d_tabA.alloc(sz);
                    m_impl->d_tabW.alloc(sz);
                    m_impl->dc6_tab_enabled = true;
                } catch (const std::exception& ex) {
                    m_impl->d_tabA.free();
                    m_impl->d_tabW.free();
                    CurcumaLogger::warn(fmt::format("ROCm: tabulated dc6dcn unavailable ({}); k_disp_uv path kept", ex.what()));
                    curcuma::reportGpuFallback("GFN-FF dc6dcn: tables unavailable, factorised uv path", ex.what());
                }
            }
            if (CurcumaLogger::get_verbosity() >= 2)
                CurcumaLogger::info(fmt::format("ROCm dc6dcn: {}", m_impl->dc6_tab_enabled
                    ? "tabulated (k_d4_partner_tables + k_dc6dcn_per_pair_tab, fused gather off)"
                    : "factorised uv (k_disp_uv + per-pair or fused gather)"));
        }
    }

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
    if (params.coulomb_implicit && params.coul_self_alp.size() == natoms) {
        // No pair list to scan: the per-atom self-energy inputs come straight from the
        // parameter set (GFNFF::generateCoulombSelfEnergyNative).
        m_coul_chi_base   = params.coul_self_chi_base;
        m_coul_gam        = params.coul_self_gam;
        m_coul_alp        = params.coul_self_alp;
        m_coul_chi_static = params.coul_self_chi_static;
        setCoulombSelfEnergyParams(m_coul_chi_base, m_coul_gam, m_coul_alp, m_coul_chi_static);
    } else if (!params.coulombs.empty()) {
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
    if (m_impl && g_rocm_small_eeq_stream == m_impl->stream) {
        g_rocm_small_eeq_stream = nullptr;
        g_rocm_small_eeq_aux = nullptr;
        g_rocm_small_eeq_event = nullptr;
    }
    if (m_impl && (m_impl->rep_device_builds > 0 || m_impl->rep_host_fallbacks > 0)
        && CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::info(fmt::format("ROCm non-bonded repulsion list: {} device rebuilds, {} host fallbacks",
                                        m_impl->rep_device_builds, m_impl->rep_host_fallbacks));
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
        if (m_impl->event_cn_ready)  (void)hipEventDestroy(m_impl->event_cn_ready);
        if (m_impl->event_eeq_order) (void)hipEventDestroy(m_impl->event_eeq_order);
        if (m_impl->h_charges_stage) (void)hipHostFree(m_impl->h_charges_stage);

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
    if (m_impl && m_impl->h_disp_flags) (void)hipHostFree(m_impl->h_disp_flags);
    if (m_impl && m_impl->h_fix_bad) (void)hipHostFree(m_impl->h_fix_bad);   // WP10
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
    // Claude Generated (Oct 2026): unchanged pair list -> device arrays unchanged (single_stream mode only).
    if (m_impl->single_stream) {
        if (m_impl->prev_alpha_valid && hb_h_idx == m_impl->prev_alpha_h && hb_b_idx == m_impl->prev_alpha_b
            && hb_rcov == m_impl->prev_alpha_rcov)
            return;
        m_impl->prev_alpha_h = hb_h_idx;
        m_impl->prev_alpha_b = hb_b_idx;
        m_impl->prev_alpha_rcov = hb_rcov;
        m_impl->prev_alpha_valid = true;
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
    // Claude Generated (Oct 2026): unchanged metadata -> device arrays unchanged (single_stream mode only).
    if (m_impl->single_stream) {
        if (m_impl->prev_bond_hb_valid && nr_hb == m_impl->prev_nr_hb && hb_H_atom == m_impl->prev_hb_H) return;
        m_impl->prev_nr_hb = nr_hb;
        m_impl->prev_hb_H = hb_H_atom;
        m_impl->prev_bond_hb_valid = true;
    }
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
    setGfnffFixModeHip(impl.fixed_point);   // WP10

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
    if (impl.fixed_point)   // WP10: integer sums -> doubles before the map reads them
        fixToDoubleHip(impl.d_hb_cn_per_atom.ptr, N, GFNFF_FIX_CN, impl.stream);

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

/// Claude Generated (Oct 2026): true when two H-bond lists agree in every field fillHBondSoA() uploads
/// (and hbCasePermutation() reads), so the device arrays built from them are identical.
static bool sameHBondInput(const std::vector<GFNFFHydrogenBond>& a, const std::vector<GFNFFHydrogenBond>& b)
{
    if (a.size() != b.size()) return false;
    for (size_t m = 0; m < a.size(); ++m) {
        const GFNFFHydrogenBond& x = a[m];
        const GFNFFHydrogenBond& y = b[m];
        if (x.i != y.i || x.j != y.j || x.k != y.k || x.case_type != y.case_type
            || x.acceptor_parent_index != y.acceptor_parent_index
            || x.q_H != y.q_H || x.q_A != y.q_A || x.q_B != y.q_B
            || x.basicity_A != y.basicity_A || x.basicity_B != y.basicity_B
            || x.acidity_A != y.acidity_A || x.acidity_B != y.acidity_B || x.r_cut != y.r_cut
            || x.neighbors_A != y.neighbors_A || x.neighbors_B != y.neighbors_B || x.neighbors_C != y.neighbors_C)
            return false;
    }
    return true;
}

/// Claude Generated (Oct 2026): the same for halogen bonds (fields of XBondSoA::upload).
static bool sameXBondInput(const std::vector<GFNFFHalogenBond>& a, const std::vector<GFNFFHalogenBond>& b)
{
    if (a.size() != b.size()) return false;
    for (size_t m = 0; m < a.size(); ++m) {
        const GFNFFHalogenBond& x = a[m];
        const GFNFFHalogenBond& y = b[m];
        if (x.i != y.i || x.j != y.j || x.k != y.k || x.basicity_B != y.basicity_B || x.acidity_X != y.acidity_X
            || x.q_X != y.q_X || x.q_B != y.q_B || x.r_cut != y.r_cut)
            return false;
    }
    return true;
}

void FFWorkspaceHip::updateHBonds(const std::vector<GFNFFHydrogenBond>& hbonds,
                                    const std::vector<int>& atom_types)
{
    // Claude Generated (Oct 2026): unchanged list -> device arrays unchanged (single_stream mode only).
    if (m_impl->single_stream) {
        if (m_impl->prev_hb_valid && sameHBondInput(hbonds, m_impl->prev_hb_input)) return;
        m_impl->prev_hb_input = hbonds;
        m_impl->prev_hb_valid = true;
    }
    // Re-upload HBond SoA after dynamic re-detection.
    // Claude Generated (March 2026): Needed because updateHBXBIfNeeded() may find
    // additional HBonds on re-detection that weren't in the initial parameter set.
    if (m_impl->hb_gather_enabled) {
        // Claude Generated (Oct 2026): triples grouped by case type, slot tables in the same order; packed by
        // the -threads host threads (m_host_threads, set through setStagedListUpload() by the wrapper)
        const std::vector<int> perm = hbCasePermutation(hbonds);
        fillHBondSoA(m_impl->hbonds, hbonds, perm.data(), atom_types, m_impl->stream, &m_impl->hb_slots, m_natoms,
                     m_host_threads);
    } else {
        m_impl->hbonds.upload(hbonds, atom_types, m_impl->stream);
    }
    updateHBChargeScaling(*m_impl);   // Claude Generated (Oct 2026)
    // Host copy of the list, read only by the wrapper's CPU-vs-GPU list comparison at verbosity >= 3.
    // Claude Generated (Oct 2026): not kept below that level (13-30 ms per refresh at 582k H bonds).
    if (CurcumaLogger::get_verbosity() >= 3)
        m_last_hbonds = hbonds;
    else
        std::vector<GFNFFHydrogenBond>().swap(m_last_hbonds);
}

/// Claude Generated (Oct 2026): grow a device buffer with 1/8 headroom (the H-bond list changes by about 1 % between
/// refreshes; ensure() without headroom reallocated on every growth, about 1 ms each at 590k triples).
template <typename T>
static void hbEnsure(CudaBuffer<T>& b, int count)
{
    if (count > b.n) b.alloc(count + count / 8 + 16);
}

/// Claude Generated (Oct 2026): k_hb_row_sort with the shared memory of the longest row (next power of two of
/// max_row, at most HB_SORT_SHM entries): 12 bytes per entry for the triples (eH given), 4 for the slot lists. Less
/// shared memory than a fixed 4096-entry buffer lets more workgroups run per compute unit.
static void hbLaunchRowSort(int n_rows, const int* row_off, int* vals, const int* eH, const int* eB, int N, int max_row,
                            hipStream_t stream)
{
    if (n_rows <= 0 || max_row <= 1) return;
    int cap = 1;
    while (cap < max_row && cap < HB_SORT_SHM) cap <<= 1;
    const size_t bytes = static_cast<size_t>(cap) * (eH ? sizeof(long long) + sizeof(int) : sizeof(int));
    const int threads = std::min(HB_SORT_THREADS, std::max(64, cap / 2));
    k_hb_row_sort<<<n_rows, threads, bytes, stream>>>(n_rows, row_off, vals, eH, eB, N, cap);
    checkCuda(hipGetLastError(), "k_hb_row_sort launch");
}

/// Claude Generated (Oct 2026): exclusive scan of n ints on the device (k_hb_scan_*); out holds n + 1 entries,
/// out[n] = total, also written to *total_out (device pointer) when given.
static void hbDeviceExclusiveScan(const int* in, int* out, int n, CudaBuffer<int>& tiles, int* total_out,
                                  hipStream_t stream)
{
    const int ntiles = std::max(1, (n + HB_SCAN_TILE - 1) / HB_SCAN_TILE);
    tiles.ensure(ntiles);
    k_hb_scan_tiles<<<ntiles, HB_SCAN_THREADS, 0, stream>>>(n, in, out, tiles.ptr);
    k_hb_scan_tile_sums<<<1, HB_SCAN_THREADS, 0, stream>>>(ntiles, tiles.ptr, out, n, total_out);
    if (n > 0) {
        const LaunchConfig c = getLaunchConfigLight(n);
        k_hb_scan_add<<<c.gridSize, c.blockSize, 0, stream>>>(n, out, tiles.ptr);
    }
    checkCuda(hipGetLastError(), "hb scan launch");
}

/// Claude Generated (Oct 2026): ROCm port of FFWorkspaceGPU::detectHBondsOnDevice() (cuda/ff_workspace_gpu.cu).
/// Same tests as the CUDA kernels (here k_hb_ab_rows, k_hb_detect_rows); the order is built without a sort of 64-bit keys: each
/// triple is counted into its row (case, A) while it is emitted, the rows get their offsets by a device scan, the
/// triples are scattered into their rows and k_hb_row_sort orders each row by (H, B). The result is the list
/// sorted by (A, H, B) (the CUDA order) and then grouped by case type, stably, i.e. exactly what fillHBondSoA()
/// with hbCasePermutation() uploads from that list. With CURCUMA_ROCM_HB_GATHER=0 the rows are the A alone (no
/// case grouping, no slot tables), as in updateHBonds(). The gradient-slot tables of the gather mode
/// (buildHBSlotTables()) are built on the device as well: slot ranges by a scan, the per-atom slot lists by
/// counting, scatter and a per-atom row sort of the slot indices (ascending, as on the host).
int FFWorkspaceHip::detectHBondsOnDevice(const GFNFFHBDetectInputs& in, std::vector<int>& hb_i,
                                         std::vector<int>& hb_j, std::vector<int>& hb_k,
                                         std::vector<int>* hb_case, std::vector<int>* hb_parent)
{
    auto& impl = *m_impl;
    HipHBDetectDev& D = impl.hbdet;
    HBondSoA& soa = impl.hbonds;
    HipHBSlotTables& T = impl.hb_slots;
    const int N = m_natoms;
    hipStream_t stream = impl.stream;
    const bool gather = impl.hb_gather_enabled;   // case-grouped order and slot tables (fillHBondSoA with T)

    // the device list replaces whatever updateHBonds() uploaded last
    impl.prev_hb_valid = false;
    std::vector<GFNFFHydrogenBond>().swap(impl.prev_hb_input);
    std::vector<GFNFFHydrogenBond>().swap(m_last_hbonds);

    // inputs (a few arrays of N and the CSR lists)
    D.bas.upload(in.basicity.data(), N);
    D.aci.upload(in.acidity.data(), N);
    D.q.upload(in.charge.data(), N);
    D.ab_ok.upload(in.ab_ok.data(), N);
    D.hyd_list.upload(in.hb_hyd_list.data(), static_cast<int>(in.hb_hyd_list.size()));
    D.bond_off.upload(in.bond_off.data(), N + 1);
    D.bond_idx.upload(in.bond_idx.data(), static_cast<int>(in.bond_idx.size()));
    D.hyd_off.upload(in.hyd_off.data(), N + 1);
    D.hyd_idx.upload(in.hyd_idx.data(), static_cast<int>(in.hyd_idx.size()));
    D.nbr_off.upload(in.nbr_off.data(), N + 1);
    D.nbr_idx.upload(in.nbr_idx.data(), static_cast<int>(in.nbr_idx.size()));
    if (!D.repz.ptr) {
        std::vector<double> repz86(86, 1.0);
        for (size_t z = 0; z < repz86.size() && z < GFNFFParameters::repz.size(); ++z) repz86[z] = GFNFFParameters::repz[z];
        D.repz.upload(repz86.data(), 86);
    }
    D.counter.ensure(1);
    D.ab_counter.ensure(1);
    D.totals.ensure(5);

    // pass 1: AB pairs (compacted, arrival order); pass 2: one workgroup per AB pair. Both
    // count past the capacity on overflow; the pass is then rerun with a larger buffer.
    D.fx.ensure(N); D.fy.ensure(N); D.fz.ensure(N);
    {
        const LaunchConfig c = getLaunchConfigLight(N, m_block_size);
        k_hb_to_float<<<c.gridSize, c.blockSize, 0, stream>>>(N, impl.coords.d_x.ptr, impl.coords.d_y.ptr,
            impl.coords.d_z.ptr, D.fx.ptr, D.fy.ptr, D.fz.ptr);
        checkCuda(hipGetLastError(), "k_hb_to_float launch");
        const int n_hyd = static_cast<int>(in.hb_hyd_list.size());
        D.hpos.ensure(std::max(1, n_hyd));
        if (n_hyd > 0) {
            const LaunchConfig ch = getLaunchConfigLight(n_hyd, m_block_size);
            k_hb_pack_hyd<<<ch.gridSize, ch.blockSize, 0, stream>>>(n_hyd, D.hyd_list.ptr, D.fx.ptr, D.fy.ptr, D.fz.ptr,
                                                                    D.hpos.ptr);
            checkCuda(hipGetLastError(), "k_hb_pack_hyd launch");
        }
    }
    int n_ab = 0;
    int ab_cap = std::max(D.ab_cap, 1 << 17);
    for (int attempt = 0; attempt < 4; ++attempt) {
        D.ab_i.ensure(ab_cap); D.ab_j.ensure(ab_cap);
        checkCuda(hipMemsetAsync(D.ab_counter.ptr, 0, sizeof(int), stream), "hb ab counter reset");
        k_hb_ab_rows<<<std::max(1, N), 256, 0, stream>>>(
            N, impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr, D.bas.ptr, D.aci.ptr, D.ab_ok.ptr,
            in.hbthr1, ab_cap, D.ab_counter.ptr, D.ab_i.ptr, D.ab_j.ptr);
        checkCuda(hipGetLastError(), "k_hb_ab_rows launch");
        D.ab_counter.download(&n_ab, 1, stream);
        checkCuda(hipStreamSynchronize(stream), "hb ab count download");
        if (n_ab <= ab_cap) break;
        ab_cap = n_ab + n_ab / 4;
    }
    D.ab_cap = std::max(D.ab_cap, ab_cap);
    if (n_ab > ab_cap) throw std::runtime_error("ROCm H-bond detection: AB pair count exceeds the buffer after reruns");

    const int R = gather ? 4 * N : N;   // rows: (case - 1) * N + A, or A
    D.row_cnt.ensure(R); D.row_off.ensure(R + 1); D.row_fill.ensure(std::max(R, N));
    int cap = std::max(D.cap, 1 << 20);
    int count = 0, max_row = 0;
    for (int attempt = 0; attempt < 4 && n_ab > 0; ++attempt) {
        D.e_A.ensure(cap); D.e_H.ensure(cap); D.e_B.ensure(cap); D.e_ct.ensure(cap); D.e_par.ensure(cap);
        checkCuda(hipMemsetAsync(D.counter.ptr, 0, sizeof(int), stream), "hb counter reset");
        checkCuda(hipMemsetAsync(D.row_cnt.ptr, 0, sizeof(int) * R, stream), "hb row count reset");
        k_hb_detect_rows<<<n_ab, HB_DETECT_THREADS, 0, stream>>>(
            n_ab, D.ab_i.ptr, D.ab_j.ptr, impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
            D.fx.ptr, D.fy.ptr, D.fz.ptr, N, impl.d_atom_types.ptr, D.bas.ptr, D.aci.ptr, D.q.ptr,
            D.hyd_list.ptr, D.hpos.ptr, static_cast<int>(in.hb_hyd_list.size()),
            D.bond_off.ptr, D.bond_idx.ptr, D.hyd_off.ptr, D.hyd_idx.ptr, D.nbr_off.ptr, D.nbr_idx.ptr,
            in.hbthr2, in.min_pair_energy, gather ? 1 : 0, cap, D.counter.ptr, D.row_cnt.ptr,
            D.e_A.ptr, D.e_H.ptr, D.e_B.ptr, D.e_ct.ptr, D.e_par.ptr);
        checkCuda(hipGetLastError(), "k_hb_detect_rows launch");
        k_hb_max<<<1, HB_SORT_THREADS, 0, stream>>>(R, D.row_cnt.ptr, D.totals.ptr + 3);
        D.counter.download(&count, 1, stream);
        checkCuda(hipMemcpyAsync(&max_row, D.totals.ptr + 3, sizeof(int), hipMemcpyDeviceToHost, stream), "hb max row");
        checkCuda(hipStreamSynchronize(stream), "hb count download");
        if (count <= cap) break;
        cap = count + count / 4;
    }
    D.cap = std::max(D.cap, cap);
    if (count > cap) throw std::runtime_error("ROCm H-bond detection: entry count exceeds the buffer after reruns");

    const int n = count;
    soa.n = n;
    hb_i.clear(); hb_j.clear(); hb_k.clear();
    if (hb_case) hb_case->clear();
    if (hb_parent) hb_parent->clear();
    if (gather) { T.nslots = 0; T.N = N; }   // an empty list must not gather stale slots (fillHBondSoA)
    if (n == 0) {
        soa.total_nb_B = 0; soa.total_nb_C = 0;
        updateHBChargeScaling(impl);
        return 0;
    }

    // order: row offsets, arrival-order scatter into the rows, then each row sorted by (H, B)
    hbDeviceExclusiveScan(D.row_cnt.ptr, D.row_off.ptr, R, D.scan_tiles, nullptr, stream);
    checkCuda(hipMemsetAsync(D.row_fill.ptr, 0, sizeof(int) * R, stream), "hb row fill reset");
    hbEnsure(D.order, n);
    {
        const LaunchConfig c = getLaunchConfigLight(n, m_block_size);
        k_hb_row_scatter<<<c.gridSize, c.blockSize, 0, stream>>>(n, D.e_A.ptr, D.e_ct.ptr, N, gather ? 1 : 0,
            D.row_off.ptr, D.row_fill.ptr, D.order.ptr);
        hbLaunchRowSort(R, D.row_off.ptr, D.order.ptr, D.e_H.ptr, D.e_B.ptr, N, max_row, stream);
        checkCuda(hipGetLastError(), "hb row order launch");
    }

    // SoA columns (the arrays of HBondSoA::upload), neighbour offsets and slot ranges by scans
    CudaBuffer<int>* icols[] = { &soa.idx_i, &soa.idx_j, &soa.idx_k, &soa.elem_A, &soa.elem_B, &soa.case_type,
                                 &soa.nb_B_count, &soa.acceptor_parent, &soa.nb_C_count };
    CudaBuffer<double>* dcols[] = { &soa.q_H, &soa.q_A, &soa.q_B, &soa.basicity_A, &soa.basicity_B,
                                    &soa.acidity_A, &soa.acidity_B, &soa.r_cut, &soa.repz_B };
    for (auto* b : icols) hbEnsure(*b, n);
    for (auto* b : dcols) hbEnsure(*b, n);
    hbEnsure(soa.nb_B_offset, n + 1);
    hbEnsure(soa.nb_C_offset, n + 1);
    if (gather) { hbEnsure(D.slot_cnt, n); hbEnsure(T.slot_off, n + 1); }
    const LaunchConfig cfg = getLaunchConfigLight(n, m_block_size);
    k_hb_fill_columns_rows<<<cfg.gridSize, cfg.blockSize, 0, stream>>>(
        n, D.order.ptr, D.e_A.ptr, D.e_H.ptr, D.e_B.ptr, D.e_ct.ptr, D.e_par.ptr, impl.d_atom_types.ptr,
        D.bas.ptr, D.aci.ptr, D.q.ptr, D.repz.ptr, D.nbr_off.ptr, D.nbr_idx.ptr,
        soa.idx_i.ptr, soa.idx_j.ptr, soa.idx_k.ptr, soa.elem_A.ptr, soa.elem_B.ptr,
        soa.q_H.ptr, soa.q_A.ptr, soa.q_B.ptr, soa.basicity_A.ptr, soa.basicity_B.ptr,
        soa.acidity_A.ptr, soa.acidity_B.ptr, soa.r_cut.ptr, soa.case_type.ptr, soa.acceptor_parent.ptr,
        soa.repz_B.ptr, soa.nb_B_count.ptr, soa.nb_C_count.ptr, gather ? D.slot_cnt.ptr : nullptr);
    checkCuda(hipGetLastError(), "k_hb_fill_columns_rows launch");
    hbDeviceExclusiveScan(soa.nb_B_count.ptr, soa.nb_B_offset.ptr, n, D.scan_tiles, D.totals.ptr + 0, stream);
    hbDeviceExclusiveScan(soa.nb_C_count.ptr, soa.nb_C_offset.ptr, n, D.scan_tiles, D.totals.ptr + 1, stream);
    if (gather)
        hbDeviceExclusiveScan(D.slot_cnt.ptr, T.slot_off.ptr, n, D.scan_tiles, D.totals.ptr + 2, stream);
    int totals[3] = { 0, 0, 0 };
    D.totals.download(totals, gather ? 3 : 2, stream);
    checkCuda(hipStreamSynchronize(stream), "hb totals download");
    soa.total_nb_B = totals[0];
    soa.total_nb_C = totals[1];
    const int nslots = gather ? totals[2] : 0;
    if (soa.total_nb_B > 0) hbEnsure(soa.nb_B_flat, soa.total_nb_B);
    if (soa.total_nb_C > 0) hbEnsure(soa.nb_C_flat, soa.total_nb_C);
    if (gather) {
        hbEnsure(D.slot_atom, nslots);
        D.atom_cnt.ensure(N + 1);
        checkCuda(hipMemsetAsync(D.atom_cnt.ptr, 0, sizeof(int) * N, stream), "hb atom count reset");
    }
    k_hb_fill_neighbors_slots<<<cfg.gridSize, cfg.blockSize, 0, stream>>>(
        n, N, soa.idx_i.ptr, soa.idx_j.ptr, soa.idx_k.ptr, soa.case_type.ptr, soa.acceptor_parent.ptr,
        D.nbr_off.ptr, D.nbr_idx.ptr, soa.nb_B_offset.ptr, soa.nb_C_offset.ptr, soa.nb_B_flat.ptr, soa.nb_C_flat.ptr,
        gather ? T.slot_off.ptr : nullptr, D.slot_atom.ptr, D.atom_cnt.ptr);
    checkCuda(hipGetLastError(), "k_hb_fill_neighbors_slots launch");

    if (gather) {
        // per-atom slot lists (buildHBSlotTables): counts -> offsets, scatter, slots ascending within each atom
        T.atom_off.ensure(N + 1);
        hbDeviceExclusiveScan(D.atom_cnt.ptr, T.atom_off.ptr, N, D.scan_tiles, nullptr, stream);
        hbEnsure(T.atom_slots, std::max(1, nslots));
        checkCuda(hipMemsetAsync(D.row_fill.ptr, 0, sizeof(int) * N, stream), "hb slot fill reset");
        const LaunchConfig cs = getLaunchConfigLight(nslots, m_block_size);
        k_hb_row_scatter<<<cs.gridSize, cs.blockSize, 0, stream>>>(nslots, D.slot_atom.ptr, nullptr, N, 0,
            T.atom_off.ptr, D.row_fill.ptr, T.atom_slots.ptr);
        k_hb_max<<<1, HB_SORT_THREADS, 0, stream>>>(N, D.atom_cnt.ptr, D.totals.ptr + 4);
        int max_slots = 0;
        checkCuda(hipMemcpyAsync(&max_slots, D.totals.ptr + 4, sizeof(int), hipMemcpyDeviceToHost, stream), "hb max slots");
        checkCuda(hipStreamSynchronize(stream), "hb max slots sync");
        hbLaunchRowSort(N, T.atom_off.ptr, T.atom_slots.ptr, nullptr, nullptr, N, max_slots, stream);
        checkCuda(hipGetLastError(), "hb slot lists launch");
        T.nslots = nslots;
        T.N = N;
        if (T.slot_g.n < 3 * nslots) T.slot_g.alloc(3 * (nslots + nslots / 8));   // headroom as buildHBSlotTables
    }
    updateHBChargeScaling(impl);

    // (A, H, B) to the host through a pinned stage (pageable copies ran at about 4 GB/s: 1.6 ms at 590k triples)
    int* pin = D.pinned(3 * static_cast<size_t>(n));
    if (pin) {
        const size_t bytes = sizeof(int) * static_cast<size_t>(n);
        checkCuda(hipMemcpyAsync(pin, soa.idx_i.ptr, bytes, hipMemcpyDeviceToHost, stream), "hb triple download");
        checkCuda(hipMemcpyAsync(pin + n, soa.idx_j.ptr, bytes, hipMemcpyDeviceToHost, stream), "hb triple download");
        checkCuda(hipMemcpyAsync(pin + 2 * static_cast<size_t>(n), soa.idx_k.ptr, bytes, hipMemcpyDeviceToHost, stream),
                  "hb triple download");
    } else {
        hb_i.resize(n); hb_j.resize(n); hb_k.resize(n);
        soa.idx_i.download(hb_i.data(), n, stream);
        soa.idx_j.download(hb_j.data(), n, stream);
        soa.idx_k.download(hb_k.data(), n, stream);
    }
    if (hb_case) { hb_case->resize(n); soa.case_type.download(hb_case->data(), n, stream); }
    if (hb_parent) { hb_parent->resize(n); soa.acceptor_parent.download(hb_parent->data(), n, stream); }
    checkCuda(hipStreamSynchronize(stream), "hb triple download sync");
    if (pin) {
        hb_i.assign(pin, pin + n);
        hb_j.assign(pin + n, pin + 2 * static_cast<size_t>(n));
        hb_k.assign(pin + 2 * static_cast<size_t>(n), pin + 3 * static_cast<size_t>(n));
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        const double cand_mb = (5.0 * D.e_A.n + D.order.n) * sizeof(int) / 1048576.0;
        const double row_mb = (3.0 * std::max(R, N) + 1.0) * sizeof(int) / 1048576.0;
        const double slot_mb = (static_cast<double>(D.slot_atom.n) + D.slot_cnt.n) * sizeof(int) / 1048576.0;
        CurcumaLogger::info(fmt::format("  [DEBUG] ROCm HB detection: {} AB pairs, {} triples, {} slots; device "
                                        "buffers: candidates {:.1f} MB (cap {}), rows {:.2f} MB, slot scratch {:.1f} MB",
                                        n_ab, n, nslots, cand_mb, D.cap, row_mb, slot_mb));
    }
    return n;
}

void FFWorkspaceHip::updateXBonds(const std::vector<GFNFFHalogenBond>& xbonds,
                                    const std::vector<int>& atom_types)
{
    // Claude Generated (Oct 2026): unchanged list -> device arrays unchanged (single_stream mode only).
    if (m_impl->single_stream && m_impl->prev_xb_valid && sameXBondInput(xbonds, m_last_xbonds)) return;
    m_impl->xbonds.upload(xbonds, atom_types, m_impl->stream);
    // Store for CPU vs GPU comparison debugging
    m_last_xbonds = xbonds;
    m_impl->prev_xb_valid = true;
}

static void dumpRepulsionListHip(const char* source, std::vector<GFNFFRepulsion> v);   // below (WP9-A)

void FFWorkspaceHip::updateRepulsion(const std::vector<GFNFFRepulsion>& bonded_reps,
                                       const std::vector<GFNFFRepulsion>& nonbonded_reps)
{
    // Claude Generated (Sep 2026): mirror of FFWorkspaceGPU::updateRepulsion (ff_workspace_gpu.cu)
    // for the ROCm plugin. See GFNFF::updateNonbondedRepulsionIfNeeded().
    // Claude Generated (Oct 2026): -gfnff.gpu_staged_list_upload (default true) takes the pinned staged path with the
    // -threads host threads, same device bytes; false, or no pinned memory, keeps the serial path below.
    dumpRepulsionListHip("host", nonbonded_reps);   // Claude Generated (Oct 2026): CURCUMA_ROCM_REP_DUMP only
    if (m_staged_list_upload
        && stagedRepulsionUploadHip(*m_impl, bonded_reps, nonbonded_reps, m_impl->rep_gather_enabled, m_natoms,
                                    m_host_threads))
        return;
    m_impl->bonded_rep.upload(bonded_reps, m_impl->stream);
    m_impl->nonbonded_rep.upload(nonbonded_reps, m_impl->stream);
    if (m_impl->rep_gather_enabled)
        buildRepulsionGatherCSR(nonbonded_reps, m_natoms, m_impl->nb_rep_csr);   // Claude Generated (Oct 2026)
}

/// Claude Generated (Oct 2026, WP9-A): test dump of a non-bonded repulsion list, sorted by (i, j), to the file
/// named by CURCUMA_ROCM_REP_DUMP (appended, one block per refresh: "# refresh <n> <source> <pairs>", then
/// "i j alpha repab r_cut" with 17 digits). Lets the device list be compared with the host list of the same geometry.
static void dumpRepulsionListHip(const char* source, std::vector<GFNFFRepulsion> v)
{
    static const char* path = std::getenv("CURCUMA_ROCM_REP_DUMP");
    if (!path || !*path) return;
    static long refresh = 0;
    std::sort(v.begin(), v.end(), [](const GFNFFRepulsion& a, const GFNFFRepulsion& b) {
        return a.i != b.i ? a.i < b.i : a.j < b.j;
    });
    FILE* f = std::fopen(path, "a");
    if (!f) return;
    std::fprintf(f, "# refresh %ld %s %zu\n", ++refresh, source, v.size());
    for (const auto& r : v)
        std::fprintf(f, "%d %d %.17g %.17g %.17g\n", r.i, r.j, r.alpha, r.repab, r.r_cut);
    std::fclose(f);
}

/// Claude Generated (Oct 2026, WP9-A): host build of the non-bonded list from the device tables and the host
/// coordinate mirror, pairs in ascending (i, j): the fallback of a failed device rebuild. Same predicate and
/// pair expressions as k_rep_rows_* / k_rep_pairs_static.
static std::vector<GFNFFRepulsion> repulsionHostBuildFromTables(const GFNFF::RepulsionDeviceInputs& in, int N,
                                                                const double* x, const double* y, const double* z,
                                                                const std::vector<int>& Z)
{
    std::vector<GFNFFRepulsion> out;
    const double cut_sq = in.build_radius_bohr * in.build_radius_bohr;
    auto find = [](const std::vector<long long>& k, long long key) -> long {
        auto it = std::lower_bound(k.begin(), k.end(), key);
        return (it != k.end() && *it == key) ? static_cast<long>(it - k.begin()) : -1;
    };
    for (int i = 0; i < N; ++i) {
        if (!in.valid[i]) continue;
        for (int j = i + 1; j < N; ++j) {
            if (!in.valid[j]) continue;
            const double dx = x[i] - x[j], dy = y[i] - y[j], dz = z[i] - z[j];
            if (dx * dx + dy * dy + dz * dz > cut_sq) continue;
            const long long key = static_cast<long long>(i) * N + j;
            if (find(in.bond_keys, key) >= 0) continue;
            double ff = 1.0;
            if (Z[i] == 1 && Z[j] == 1) {
                ff = in.hhfac;
                const long h = find(in.hh_keys, key);
                if (h >= 0) {
                    if (in.hh_dist[h] == 2) ff *= in.hh13;
                    else if (in.hh_dist[h] == 3) ff *= in.hh14;
                }
            } else if ((Z[i] == 1 && in.metal[j]) || (Z[j] == 1 && in.metal[i])) {
                ff = 0.85;
            } else if ((Z[i] == 1 && Z[j] == 6) || (Z[j] == 1 && Z[i] == 6)) {
                ff = 0.91;
            } else if ((Z[i] == 1 && Z[j] == 8) || (Z[j] == 1 && Z[i] == 8)) {
                ff = 1.04;
            }
            GFNFFRepulsion r;
            r.i = i;
            r.j = j;
            r.alpha = std::sqrt(in.dum[i] * in.dum[j]) * ff;
            r.repab = in.repz[i] * in.repz[j] * in.repscaln;
            r.r_cut = in.r_cut_bohr;
            out.push_back(r);
        }
    }
    return out;
}

/**
 * @brief WP9-A: rebuild the non-bonded repulsion list on the device from the current device coordinates.
 *
 * Claude Generated (Oct 2026). Port of FFWorkspaceGPU::generateNonbondedRepulsionOnGPU with the ordered
 * row build of the ROCm D4 list (k_rep_rows_count, two row-offset scans on the host, k_rep_rows_build,
 * k_rep_pairs_static, k_disp_csr_lower for the lower CSR half): pairs in ascending (i, j) and the gather CSR
 * built on the device (no host list, no index download). The bonded list is not touched (topology only).
 * A failure (or CURCUMA_ROCM_FAIL_REPDEV=1) builds the list on the host from the same tables and the host
 * coordinate mirror, uploads it as updateRepulsion() does, counts a GPU fallback and keeps the host build
 * for the rest of the run. Selected by -gfnff.gpu_rep_pairs_on_device (wrapper, trait has_device_repulsion_build).
 */
void FFWorkspaceHip::generateNonbondedRepulsionOnGPU(const GFNFF::RepulsionDeviceInputs& in)
{
    auto& impl = *m_impl;
    const int N = m_natoms;
    if (N <= 1) return;
    static const bool fail_inject = [] {
        const char* e = std::getenv("CURCUMA_ROCM_FAIL_REPDEV");
        return e && *e && std::string(e) != "0";
    }();
    static const bool dump = [] { const char* e = std::getenv("CURCUMA_ROCM_REP_DUMP"); return e && *e; }();
    const bool gather = impl.rep_gather_enabled;
    bool ok = !impl.rep_device_off;
    int n = 0;
    try {
        if (ok && fail_inject) throw std::runtime_error("fault injection CURCUMA_ROCM_FAIL_REPDEV");
        if (ok) {
            impl.rep_valid.upload(in.valid);
            impl.rep_metal.upload(in.metal);
            impl.rep_dum.upload(in.dum);
            impl.rep_repz.upload(in.repz);
            if (!in.bond_keys.empty()) impl.rep_bond_keys.upload(in.bond_keys);
            if (!in.hh_keys.empty()) { impl.rep_hh_keys.upload(in.hh_keys); impl.rep_hh_dist.upload(in.hh_dist); }
            RepDevTablesHip t;
            t.valid = impl.rep_valid.ptr;  t.metal = impl.rep_metal.ptr;
            t.dum = impl.rep_dum.ptr;      t.repz = impl.rep_repz.ptr;
            t.bond_keys = impl.rep_bond_keys.ptr;  t.n_bond_keys = static_cast<int>(in.bond_keys.size());
            t.hh_keys = impl.rep_hh_keys.ptr;      t.hh_dist = impl.rep_hh_dist.ptr;
            t.n_hh_keys = static_cast<int>(in.hh_keys.size());
            t.repscaln = in.repscaln;  t.hhfac = in.hhfac;  t.hh13 = in.hh13;  t.hh14 = in.hh14;
            t.r_cut = in.r_cut_bohr;
            const double cut_sq = in.build_radius_bohr * in.build_radius_bohr;

            impl.rep_lo.ensure(N);
            impl.rep_hi.ensure(N);
            impl.rep_rowoff.ensure(N + 1);
            k_rep_rows_count<<<N, DISP_ROW_BLOCK, 0, impl.stream>>>(
                N, impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr, cut_sq, t,
                impl.rep_lo.ptr, impl.rep_hi.ptr);
            checkCuda(hipGetLastError(), "k_rep_rows_count launch");
            std::vector<int> h_lo(N), h_hi(N);
            impl.rep_lo.download(h_lo.data(), N, impl.stream);
            impl.rep_hi.download(h_hi.data(), N, impl.stream);
            checkCuda(hipStreamSynchronize(impl.stream), "k_rep_rows_count sync");
            std::vector<int> rowoff(N + 1, 0), off(N + 1, 0);
            for (int a = 0; a < N; ++a) {
                rowoff[a + 1] = rowoff[a] + h_hi[a];
                off[a + 1]    = off[a] + h_lo[a] + h_hi[a];
            }
            n = rowoff[N];
            if (off[N] != 2 * n)
                throw std::runtime_error(fmt::format("CSR size {} != 2 x {} pairs", off[N], n));

            auto& nb = impl.nonbonded_rep;
            if (nb.idx_i.n < n) {   // grow with headroom; the SoA keeps its capacity across rebuilds
                const int cap = n + n / 8 + 16;
                nb.idx_i.alloc(cap); nb.idx_j.alloc(cap); nb.alpha.alloc(cap); nb.repab.alloc(cap); nb.r_cut.alloc(cap);
            }
            nb.n = n;
            if (gather) {
                auto& c = impl.nb_rep_csr;
                c.N = N;
                c.off.upload(off);
                if (c.partner.n < 2 * n) { const int cap = 2 * n + n / 4 + 16; c.partner.alloc(cap); c.pairidx.alloc(cap); }
            }
            if (n > 0) {
                impl.rep_rowoff.upload(rowoff);
                k_rep_rows_build<<<N, DISP_ROW_BLOCK, 0, impl.stream>>>(
                    N, impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr, cut_sq, t,
                    impl.rep_lo.ptr, impl.rep_rowoff.ptr, gather ? impl.nb_rep_csr.off.ptr : nullptr,
                    nb.idx_i.ptr, nb.idx_j.ptr,
                    gather ? impl.nb_rep_csr.partner.ptr : nullptr, gather ? impl.nb_rep_csr.pairidx.ptr : nullptr);
                checkCuda(hipGetLastError(), "k_rep_rows_build launch");
                LaunchConfig cfg = getLaunchConfig(n, m_block_size);
                k_rep_pairs_static<<<cfg.gridSize, cfg.blockSize, 0, impl.stream>>>(
                    n, N, nb.idx_i.ptr, nb.idx_j.ptr, impl.d_atom_types.ptr, t,
                    nb.alpha.ptr, nb.repab.ptr, nb.r_cut.ptr);
                checkCuda(hipGetLastError(), "k_rep_pairs_static launch");
                if (gather) {
                    k_disp_csr_lower<<<N, DISP_ROW_BLOCK, 0, impl.stream>>>(
                        N, impl.rep_lo.ptr, impl.rep_rowoff.ptr, impl.nb_rep_csr.off.ptr,
                        nb.idx_j.ptr, impl.nb_rep_csr.partner.ptr, impl.nb_rep_csr.pairidx.ptr);
                    checkCuda(hipGetLastError(), "k_disp_csr_lower (repulsion) launch");
                }
            }
            // The repulsion kernels run on another stream: finish the build here.
            checkCuda(hipStreamSynchronize(impl.stream), "k_rep_rows_build sync");
            ++impl.rep_device_builds;
            if (dump) {
                std::vector<GFNFFRepulsion> v(n);
                std::vector<int> hi(n), hj(n);
                std::vector<double> ha(n), hr(n), hc(n);
                nb.idx_i.download(hi.data(), n, impl.stream);
                nb.idx_j.download(hj.data(), n, impl.stream);
                nb.alpha.download(ha.data(), n, impl.stream);
                nb.repab.download(hr.data(), n, impl.stream);
                nb.r_cut.download(hc.data(), n, impl.stream);
                checkCuda(hipStreamSynchronize(impl.stream), "repulsion dump sync");
                for (int k = 0; k < n; ++k) v[k] = GFNFFRepulsion{hi[k], hj[k], ha[k], hr[k], hc[k]};
                dumpRepulsionListHip("device", std::move(v));
            }
            if (CurcumaLogger::get_verbosity() >= 3)
                CurcumaLogger::info(fmt::format("  ROCm non-bonded repulsion list: {} pairs built on the device{}",
                                                n, gather ? ", gather CSR on the device" : ""));
            return;
        }
    } catch (const std::exception& e) {
        (void)hipGetLastError();   // clear a sticky launch error before the host path
        impl.rep_device_off = true;
        CurcumaLogger::warn(fmt::format("ROCm: device repulsion list rebuild failed ({}); host build from now on", e.what()));
        curcuma::reportGpuFallback("GFN-FF repulsion list: device rebuild failed, host build", e.what());
    }
    // Host build from the same tables (also every later rebuild once the device path failed).
    if (!impl.h_geom_ok) {
        CurcumaLogger::error("ROCm: repulsion list host build without host coordinates; list kept");
        return;
    }
    std::vector<int> Z(N);
    impl.d_atom_types.download(Z.data(), N, impl.stream);
    checkCuda(hipStreamSynchronize(impl.stream), "atom types download");
    std::vector<GFNFFRepulsion> v = repulsionHostBuildFromTables(in, N, m_h_x, m_h_y, m_h_z, Z);
    ++impl.rep_host_fallbacks;
    impl.nonbonded_rep.upload(v, impl.stream);
    if (gather) buildRepulsionGatherCSR(v, N, impl.nb_rep_csr);
    if (dump) dumpRepulsionListHip("host-fallback", v);
}

void FFWorkspaceHip::updateDispersion(const std::vector<GFNFFDispersion>& pairs)
{
    // Claude Generated (Sep 2026): mirror of FFWorkspaceGPU::updateDispersion (ff_workspace_gpu.cu);
    // the ROCm DispersionSoA::upload() additionally rebuilds the per-atom gather CSR from the new
    // pairs. Mechanical mirror, unverified (no ROCm SDK here).
    auto& impl = *m_impl;
    impl.disp.upload(pairs, m_natoms, impl.stream);
    { impl.disp_uv_valid = false; impl.disp_pair_c6_current = true; }  // Oct 2026: new pair list, d_disp_uv stale
    const int nd = impl.disp.n;
    if (impl.dc6dcn_gpu_ready && nd > 0) {
        impl.disp.dc6dcn_ij.ensure(nd);
        impl.disp.dc6dcn_ji.ensure(nd);
    }
    m_disp_idx_i_host.resize(nd);
    m_disp_idx_j_host.resize(nd);
    m_h_dc6dcn_ij.assign(nd, 0.0);
    m_h_dc6dcn_ji.assign(nd, 0.0);
    for (int k = 0; k < nd; ++k) {
        m_disp_idx_i_host[k] = pairs[k].i;
        m_disp_idx_j_host[k] = pairs[k].j;
    }
}

void FFWorkspaceHip::updateCoulombPairs(const std::vector<GFNFFCoulomb>& pairs)
{
    // Claude Generated (Sep 2026): mirror of FFWorkspaceGPU::updateCoulombPairs; the ROCm
    // CoulombSoA::upload() also rebuilds its gather CSR. Unverified (no ROCm SDK here).
    m_impl->coulomb.upload(pairs, m_natoms, m_impl->stream);
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
/// Claude Generated (Oct 2026): host copy of the geometry last uploaded by setGeometry() (the values the
/// device coordinates hold), for the host displacement tests of the single-stream mode.
static bool copyHostGeometryMirror(bool geom_ok, int N, const double* hx, const double* hy, const double* hz,
                                   std::vector<double>& dst)
{
    if (!geom_ok || !hx || !hy || !hz) return false;
    dst.resize(3 * static_cast<size_t>(N));
    std::memcpy(dst.data(), hx, N * sizeof(double));
    std::memcpy(dst.data() + N, hy, N * sizeof(double));
    std::memcpy(dst.data() + 2 * N, hz, N * sizeof(double));
    return true;
}

/// Claude Generated (Oct 2026): true when some atom of (hx, hy, hz) is farther than sqrt(thr_sq) from the
/// reference ref (x, y, z blocks of N); the per-atom test of k_check_displacement.
static bool hostDisplacementExceeded(int N, const double* hx, const double* hy, const double* hz,
                                     const std::vector<double>& ref, double thr_sq)
{
    const double* rx = ref.data();
    const double* ry = rx + N;
    const double* rz = ry + N;
    for (int i = 0; i < N; ++i) {
        const double dx = hx[i] - rx[i];
        const double dy = hy[i] - ry[i];
        const double dz = hz[i] - rz[i];
        if (dx*dx + dy*dy + dz*dz > thr_sq) return true;
    }
    return false;
}

void FFWorkspaceHip::generateCNPairListOnGPU()
{
    auto& impl = *m_impl;
    const int N = m_natoms;
    if (N <= 0) return;

    impl.cn_stale_cached = -1;   // Oct 2026: the pair-list reference changes below

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
        2.0 * m_cn_pair_skin_bohr, // K1 (Oct 2026): Verlet skin, 0 = legacy list
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
        2.0 * m_cn_pair_skin_bohr, // K1 (Oct 2026): Verlet skin, 0 = legacy list
        impl.d_cn_pair_counter.ptr,
        impl.d_cn_idx_i.ptr,
        impl.d_cn_idx_j.ptr,
        impl.d_cn_rcov_sum.ptr);
    checkCuda(hipStreamSynchronize(impl.stream), "k_generate_cn_pairs_write sync");

    m_cn_pairs_on_gpu = true;
    ++m_cn_pair_builds;

    // K1 (Oct 2026, Claude Generated): remember the geometry this list was built from, so
    // cnPairListStale() can measure the largest atom displacement since this build.
    if (impl.cn_ref_coords.empty()) impl.cn_ref_coords.alloc(N);
    impl.h_cnref_ok = copyHostGeometryMirror(impl.h_geom_ok, N, m_h_x, m_h_y, m_h_z, impl.h_cnref_xyz);
    checkCuda(hipMemcpyAsync(impl.cn_ref_coords.d_rx.ptr, impl.coords.d_x.ptr, N * sizeof(double),
                             hipMemcpyDeviceToDevice, impl.stream), "CN pair ref: D2D x");
    checkCuda(hipMemcpyAsync(impl.cn_ref_coords.d_ry.ptr, impl.coords.d_y.ptr, N * sizeof(double),
                             hipMemcpyDeviceToDevice, impl.stream), "CN pair ref: D2D y");
    checkCuda(hipMemcpyAsync(impl.cn_ref_coords.d_rz.ptr, impl.coords.d_z.ptr, N * sizeof(double),
                             hipMemcpyDeviceToDevice, impl.stream), "CN pair ref: D2D z");

    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::info(fmt::format("  GPU CN pair list: {} pairs generated on GPU (build {})",
                                        h_count, m_cn_pair_builds));
}

// ============================================================================
// K1 (Oct 2026, Claude Generated): staleness test of the CN pair list (port of the CUDA
// FFWorkspaceGPU::cnPairListStale). True when some atom moved more than max_disp_bohr
// since the last generateCNPairListOnGPU() (or when no list reference exists yet). Uses
// the flag kernel of checkDisplacement(), but against the pair-list reference geometry,
// which is independent of the topology reference used by checkDisplacement().
// ============================================================================
bool FFWorkspaceHip::cnPairListStale(double max_disp_bohr)
{
    if (!m_impl) return true;
    auto& impl = *m_impl;
    if (impl.cn_ref_coords.empty() || impl.coords.empty()) return true;
    // Claude Generated (Oct 2026): lean path - the flag was computed by checkDisplacement() for this
    // geometry and threshold; otherwise compute it here (one launch, pinned read-back).
    if (impl.lean_step && impl.h_disp_flags) {
        if (impl.cn_stale_cached >= 0 && impl.cn_stale_thr == max_disp_bohr) {
            const int r = impl.cn_stale_cached;
            impl.cn_stale_cached = -1;
            return r != 0;
        }
        impl.cn_stale_thr = max_disp_bohr;
        if (impl.single_stream && impl.h_cnref_ok && impl.h_geom_ok)   // host test (see checkDisplacement)
            return hostDisplacementExceeded(m_natoms, m_h_x, m_h_y, m_h_z, impl.h_cnref_xyz,
                                            max_disp_bohr * max_disp_bohr);
        k_check_displacement_pair<<<1, 1024, 0, impl.stream>>>(
            m_natoms, impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
            nullptr, nullptr, nullptr, 0.0,
            impl.cn_ref_coords.d_rx.ptr, impl.cn_ref_coords.d_ry.ptr, impl.cn_ref_coords.d_rz.ptr,
            max_disp_bohr * max_disp_bohr, impl.d_disp_flags2.ptr);
        checkCuda(hipMemcpyAsync(impl.h_disp_flags, impl.d_disp_flags2.ptr, 2 * sizeof(int),
                                 hipMemcpyDeviceToHost, impl.stream), "cnPairListStale: download flags");
        checkCuda(hipStreamSynchronize(impl.stream), "cnPairListStale: sync");
        return impl.h_disp_flags[1] != 0;
    }
    if (!impl.d_cn_disp_flag.ptr) impl.d_cn_disp_flag.alloc(1);

    int zero = 0;
    checkCuda(hipMemcpyAsync(impl.d_cn_disp_flag.ptr, &zero, sizeof(int),
                             hipMemcpyHostToDevice, impl.stream), "cnPairListStale: zero flag");
    LaunchConfig cfg = getLaunchConfig(m_natoms, m_block_size);
    k_check_displacement<<<cfg.gridSize, cfg.blockSize, 0, impl.stream>>>(
        m_natoms,
        impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
        impl.cn_ref_coords.d_rx.ptr, impl.cn_ref_coords.d_ry.ptr, impl.cn_ref_coords.d_rz.ptr,
        max_disp_bohr * max_disp_bohr,
        impl.d_cn_disp_flag.ptr);
    int result = 0;
    checkCuda(hipMemcpyAsync(&result, impl.d_cn_disp_flag.ptr, sizeof(int),
                             hipMemcpyDeviceToHost, impl.stream), "cnPairListStale: download flag");
    checkCuda(hipStreamSynchronize(impl.stream), "cnPairListStale: sync");
    return result != 0;
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

    // Claude Generated (Oct 2026): k_dispersion_gather_uv evaluates C6(a,p) from either atom's side,
    // which needs C6ref(Zi,Zj,ri,rj) == C6ref(Zj,Zi,rj,ri). The D4 generator stores both halves with the
    // same value (d4param_generator.cpp); check it bitwise instead of relying on it.
    {
        const size_t R = D4_MAX_REF, E = D4_MAX_ELEM;
        bool sym = (c6_flat.size() >= E * E * R * R);
        for (size_t ei = 0; sym && ei < E; ++ei)
            for (size_t ej = ei; sym && ej < E; ++ej)
                for (size_t ri = 0; sym && ri < R; ++ri)
                    for (size_t rj = 0; rj < R; ++rj)
                        if (c6_flat[((ei * E + ej) * R + ri) * R + rj] != c6_flat[((ej * E + ei) * R + rj) * R + ri]) { sym = false; break; }
        impl.c6_table_symmetric = sym;
        const char* env = std::getenv("CURCUMA_ROCM_DISP_FUSED");
        impl.disp_fused_enabled = !(env && std::string(env) == "0");
        if (CurcumaLogger::get_verbosity() >= 3)
            CurcumaLogger::info(fmt::format("  D4 C6 table symmetric: {}, fused dispersion gather: {}",
                                            sym ? "yes" : "no", (sym && impl.disp_fused_enabled) ? "on" : "off"));
    }

    // Upload refn array (118 ints) to constant memory (Phase 8: broadcast via L1)
    upload_refn_const(refn.data(), static_cast<int>(refn.size()));

    // Pre-allocate gw/dgw GPU buffers
    const int gw_size = m_natoms * D4_MAX_REF;
    impl.d_gw.alloc(gw_size);
    impl.d_dgw.alloc(gw_size);

    // Pre-allocate dc6dcn output buffers (written by k_dc6dcn_per_pair_uv)
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

/// Claude Generated (Oct 2026): factorised C6 / dC6dCN per pair from the current d_gw / d_dgw:
/// k_disp_uv (per atom and element type) then k_dc6dcn_per_pair_uv (per pair). Writes C6 only
/// when c6_update is set (per-step C6 refresh), as the former k_dc6dcn_per_pair did.
/// With with_uv false only the per-pair half runs (d_disp_uv must be current).
static void launchDispC6PerPair(FFWorkspaceHipImpl& impl, int natoms, bool c6_update,
                                int block_size, hipStream_t stream, bool with_uv = true)
{
    const int nd = impl.disp.n;
    if (nd == 0 || natoms <= 0 || impl.disp_ntypes <= 0) return;
    if (with_uv && impl.dc6_tab_enabled) {
        // WP9-B (Claude Generated, Oct 2026): tabulated dc6dcn; d_disp_uv is not written on this path.
        LaunchConfig cfg_t = getLaunchConfig(natoms * impl.disp_ntypes, block_size);
        k_d4_partner_tables<<<cfg_t.gridSize, cfg_t.blockSize, 0, stream>>>(
            natoms, impl.disp_ntypes, impl.d_ctype_elem.ptr, impl.d_atom_types.ptr,
            impl.d_gw.ptr, impl.d_c6_flat.ptr, impl.d_tabA.ptr, impl.d_tabW.ptr);
        impl.disp_uv_valid = false;
        impl.disp_pair_c6_current = true;
        LaunchConfig cfg_p = getLaunchConfig(nd, block_size);
        k_dc6dcn_per_pair_tab<<<cfg_p.gridSize, cfg_p.blockSize, 0, stream>>>(
            nd, impl.disp.idx_i.ptr, impl.disp.idx_j.ptr, impl.d_atom_types.ptr, impl.d_atom_ctype.ptr,
            impl.disp_ntypes, impl.d_gw.ptr, impl.d_dgw.ptr, impl.d_tabA.ptr, impl.d_tabW.ptr,
            impl.disp.dc6dcn_ij.ptr, impl.disp.dc6dcn_ji.ptr, c6_update ? impl.disp.C6.ptr : nullptr);
        return;
    }
    if (with_uv) {
        LaunchConfig cfg_uv = getLaunchConfig(natoms * impl.disp_ntypes, block_size);
        k_disp_uv<<<cfg_uv.gridSize, cfg_uv.blockSize, 0, stream>>>(
            natoms, impl.disp_ntypes, impl.d_atom_types.ptr, impl.d_ctype_elem.ptr,
            impl.d_gw.ptr, impl.d_dgw.ptr, impl.d_c6_flat.ptr, impl.d_disp_uv.ptr);
        impl.disp_uv_valid = true;
    }
    impl.disp_pair_c6_current = true;
    LaunchConfig cfg = getLaunchConfig(nd, block_size);
    k_dc6dcn_per_pair_uv<<<cfg.gridSize, cfg.blockSize, 0, stream>>>(
        nd, impl.disp_ntypes,
        impl.disp.idx_i.ptr, impl.disp.idx_j.ptr, impl.d_atom_ctype.ptr,
        impl.d_disp_uv.ptr, impl.d_gw.ptr, impl.d_dgw.ptr,
        impl.disp.dc6dcn_ij.ptr, impl.disp.dc6dcn_ji.ptr,
        c6_update ? impl.disp.C6.ptr : nullptr);
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

    // Factorised C6 / dc6dcn per pair (Oct 2026; refn read from constant memory d_refn_const)
    launchDispC6PerPair(impl, N, m_disp_c6_update, m_block_size, impl.stream);

    // No sync needed here: k_disp_uv + k_dc6dcn_per_pair_uv run on impl.stream (main stream).
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
// WP-A (Jun 2026): on-device D4 dispersion pair-list build. HIP mirror of
// FFWorkspaceGPU::generateDispersionPairListOnGPU (cuda/ff_workspace_gpu.cu).
// Replaces the host O(N^2) GenerateDispersionPairsNative loop + the per-build H2D
// upload of the pair SoA. ROCm-specific extra step: the per-atom CSR adjacency
// (k_dispersion_gather). Default (Oct 2026): ordered build k_disp_rows_* with the CSR on
// the device; CURCUMA_GPU_DISP_SORT=0: atomic-slot build, CSR rebuilt on the host.
// ============================================================================

void FFWorkspaceHip::generateDispersionPairListOnGPU(
    const std::vector<double>& gw_flat,
    const std::vector<double>& sqrtzr4r2,
    const std::vector<double>& topo_charges,
    double a1, double a2, double cutoff_bohr)
{
    auto& impl = *m_impl;
    const int N = m_natoms;
    if (N <= 0) return;
    if (!impl.dc6dcn_gpu_ready) {
        if (CurcumaLogger::get_verbosity() >= 1)
            CurcumaLogger::warn("GPU disp pair build skipped: C6 reference table not uploaded");
        return;
    }

    // Use the HOST Gaussian weights (the same m_gaussian_weights getChargeWeightedC6
    // uses) for the C6 contraction so the result is bit-identical to the host list.
    // gw_flat is laid out [atom*D4_MAX_REF + ref], zero-padded to D4_MAX_REF.
    if (static_cast<int>(gw_flat.size()) >= N * D4_MAX_REF)
        impl.d_gw.upload(gw_flat.data(), N * D4_MAX_REF);

    // Upload the element/atom tables the build kernel needs.
    impl.d_sqrtzr4r2.upload(sqrtzr4r2.data(), static_cast<int>(sqrtzr4r2.size()));
    impl.d_zeta_zeff.upload(GFNFFParameters::zeta_zeff.data(),
                            static_cast<int>(GFNFFParameters::zeta_zeff.size()));
    impl.d_zeta_c.upload(GFNFFParameters::zeta_c.data(),
                         static_cast<int>(GFNFFParameters::zeta_c.size()));
    {
        std::vector<double> q(N, 0.0);
        for (int i = 0; i < N && i < static_cast<int>(topo_charges.size()); ++i)
            q[i] = topo_charges[i];
        impl.d_topo_q.upload(q.data(), N);
    }

    if (!impl.d_disp_pair_counter.ptr) impl.d_disp_pair_counter.alloc(1);

    // Claude Generated (Oct 2026): ordered build (k_disp_rows_*): pairs in the host order i*N+j and
    // the gather CSR built on the device, no index download. CURCUMA_GPU_DISP_SORT=0 keeps the
    // arrival-order build below (atomic slot, host CSR rebuild), as on CUDA.
    const char* env_sort = std::getenv("CURCUMA_GPU_DISP_SORT");
    if (!(env_sort && env_sort[0] == '0')) {
        const double cut_sq = cutoff_bohr * cutoff_bohr;
        CudaBuffer<int> d_lo, d_hi, d_rowoff;
        d_lo.alloc(N);
        d_hi.alloc(N);
        k_disp_rows_count<<<N, DISP_ROW_BLOCK, 0, impl.stream>>>(
            N, impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
            impl.d_atom_types.ptr, cut_sq, d_lo.ptr, d_hi.ptr);
        std::vector<int> h_lo(N), h_hi(N);
        d_lo.download(h_lo.data(), N, impl.stream);
        d_hi.download(h_hi.data(), N, impl.stream);
        checkCuda(hipStreamSynchronize(impl.stream), "k_disp_rows_count sync");
        std::vector<int> rowoff(N + 1, 0), off(N + 1, 0);
        for (int a = 0; a < N; ++a) {
            rowoff[a + 1] = rowoff[a] + h_hi[a];
            off[a + 1]    = off[a] + h_lo[a] + h_hi[a];
        }
        const int n = rowoff[N];
        if (off[N] != 2 * n)
            throw std::runtime_error(fmt::format("GPU D4 pair build: CSR size {} != 2 x {} pairs", off[N], n));

        { impl.disp_uv_valid = false; impl.disp_pair_c6_current = true; }  // new pair list, d_disp_uv stale
        impl.disp.n = n;
        impl.disp.csr_N = N;
        m_disp_idx_i_host.clear();   // no reader on this path (the host CSR rebuild below is skipped)
        m_disp_idx_j_host.clear();
        // capacity only: the kernels use impl.disp.n, so the buffers are not shrunk (no free + malloc per rebuild)
        impl.disp.idx_i.ensure(n);
        impl.disp.idx_j.ensure(n);
        impl.disp.C6.ensure(n);
        impl.disp.r4r2ij.ensure(n);
        impl.disp.r0_sq.ensure(n);
        impl.disp.zetac6.ensure(n);
        impl.disp.r_cut.ensure(n);
        impl.disp.dc6dcn_ij.ensure(n);
        impl.disp.dc6dcn_ji.ensure(n);
        impl.disp.csr_partner.ensure(2 * n);
        impl.disp.csr_pairidx.ensure(2 * n);
        impl.disp.csr_is_i.ensure(2 * n);
        impl.disp.csr_off.upload(off, impl.stream);
        if (n > 0) {
            d_rowoff.upload(rowoff, impl.stream);
            k_disp_rows_build<<<N, DISP_ROW_BLOCK, 0, impl.stream>>>(
                N, impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
                impl.d_atom_types.ptr, cut_sq,
                d_lo.ptr, d_rowoff.ptr, impl.disp.csr_off.ptr,
                impl.disp.idx_i.ptr, impl.disp.idx_j.ptr,
                impl.disp.csr_partner.ptr, impl.disp.csr_pairidx.ptr, impl.disp.csr_is_i.ptr);
            LaunchConfig cfg_st = getLaunchConfig(n, m_block_size);
            k_disp_pairs_static<<<cfg_st.gridSize, cfg_st.blockSize, 0, impl.stream>>>(
                n, impl.disp.idx_i.ptr, impl.disp.idx_j.ptr, impl.d_atom_types.ptr,
                impl.d_gw.ptr, impl.d_c6_flat.ptr,
                impl.d_sqrtzr4r2.ptr, impl.d_topo_q.ptr,
                impl.d_zeta_zeff.ptr, impl.d_zeta_c.ptr,
                a1, a2,
                impl.disp.C6.ptr, impl.disp.r4r2ij.ptr, impl.disp.r0_sq.ptr,
                impl.disp.zetac6.ptr, impl.disp.r_cut.ptr);
            k_disp_csr_lower<<<N, DISP_ROW_BLOCK, 0, impl.stream>>>(
                N, d_lo.ptr, d_rowoff.ptr, impl.disp.csr_off.ptr,
                impl.disp.idx_j.ptr, impl.disp.csr_partner.ptr, impl.disp.csr_pairidx.ptr);
            checkCuda(hipGetLastError(), "k_disp_rows_build launch");
        }
        // d_lo / d_rowoff are freed on return: finish the kernels first.
        checkCuda(hipStreamSynchronize(impl.stream), "k_disp_rows_build sync");
        impl.disp_built_on_gpu = true;
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info(fmt::format(
                "  GPU D4 dispersion pair list: {} pairs built on GPU in host order, gather CSR on GPU", n));
        return;
    }

    const double cutoff_sq = cutoff_bohr * cutoff_bohr;
    const int n_total = N * (N - 1) / 2;
    LaunchConfig cfg = getLaunchConfig(n_total, m_block_size);

    // Pass 1: count surviving pairs.
    int h_count = 0;
    impl.d_disp_pair_counter.upload(&h_count, 1);
    k_disp_pairs_count<<<cfg.gridSize, cfg.blockSize, 0, impl.stream>>>(
        N, impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
        impl.d_atom_types.ptr, cutoff_sq, impl.d_disp_pair_counter.ptr);
    checkCuda(hipStreamSynchronize(impl.stream), "k_disp_pairs_count sync");
    impl.d_disp_pair_counter.download(&h_count, 1);
    checkCuda(hipStreamSynchronize(impl.stream), "download disp pair count");

    { impl.disp_uv_valid = false; impl.disp_pair_c6_current = true; }  // Oct 2026: new pair list, d_disp_uv stale
    // (Re)allocate the DispersionSoA to the device-counted pair count.
    impl.disp.n = h_count;
    impl.disp.idx_i.alloc(h_count);
    impl.disp.idx_j.alloc(h_count);
    impl.disp.C6.alloc(h_count);
    impl.disp.r4r2ij.alloc(h_count);
    impl.disp.r0_sq.alloc(h_count);
    impl.disp.zetac6.alloc(h_count);
    impl.disp.r_cut.alloc(h_count);
    impl.disp.dc6dcn_ij.alloc(h_count);
    impl.disp.dc6dcn_ji.alloc(h_count);

    if (h_count == 0) {
        impl.disp_built_on_gpu = true;
        m_disp_idx_i_host.clear();
        m_disp_idx_j_host.clear();
        impl.disp.csr_N = N;
        return;
    }

    // Pass 2: write idx + per-pair static data (C6 contraction, r4r2ij, R0², zetac6).
    int h_zero = 0;
    impl.d_disp_pair_counter.upload(&h_zero, 1);
    k_disp_pairs_build<<<cfg.gridSize, cfg.blockSize, 0, impl.stream>>>(
        N, impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
        impl.d_atom_types.ptr, cutoff_sq,
        impl.d_gw.ptr, impl.d_c6_flat.ptr,
        impl.d_sqrtzr4r2.ptr, impl.d_topo_q.ptr,
        impl.d_zeta_zeff.ptr, impl.d_zeta_c.ptr,
        a1, a2,
        impl.d_disp_pair_counter.ptr,
        impl.disp.idx_i.ptr, impl.disp.idx_j.ptr,
        impl.disp.C6.ptr, impl.disp.r4r2ij.ptr, impl.disp.r0_sq.ptr,
        impl.disp.zetac6.ptr, impl.disp.r_cut.ptr);
    checkCuda(hipStreamSynchronize(impl.stream), "k_disp_pairs_build sync");

    // Keep the host idx mirror consistent (CPU dc6dcn fallback / displacement consumers
    // AND the CSR rebuild below). One download per build (topology-change frequency).
    m_disp_idx_i_host.resize(h_count);
    m_disp_idx_j_host.resize(h_count);
    impl.disp.idx_i.download(m_disp_idx_i_host.data(), h_count);
    impl.disp.idx_j.download(m_disp_idx_j_host.data(), h_count);
    checkCuda(hipStreamSynchronize(impl.stream), "download disp idx mirror");

    // ROCm-specific: rebuild the per-atom CSR adjacency for k_dispersion_gather from the
    // device-built pair list — identical arithmetic to DispersionSoA::upload. The dynamic
    // dc6dcn arrays stay pair-indexed; csr_pairidx/csr_is_i index them per directed edge.
    {
        const int n = h_count;
        const std::vector<int>& hi = m_disp_idx_i_host;
        const std::vector<int>& hj = m_disp_idx_j_host;
        impl.disp.csr_N = N;
        std::vector<int> deg(N, 0);
        for (int k = 0; k < n; ++k) { ++deg[hi[k]]; ++deg[hj[k]]; }
        std::vector<int> off(N + 1, 0);
        for (int a = 0; a < N; ++a) off[a + 1] = off[a] + deg[a];
        const int tot = off[N];                            // == 2*n
        std::vector<int> h_part(tot), h_pidx(tot), h_isi(tot);
        std::vector<int> cur(off.begin(), off.end() - 1);
        for (int k = 0; k < n; ++k) {
            const int i = hi[k], j = hj[k];
            int pi = cur[i]++; h_part[pi] = j; h_pidx[pi] = k; h_isi[pi] = 1;
            int pj = cur[j]++; h_part[pj] = i; h_pidx[pj] = k; h_isi[pj] = 0;
        }
        impl.disp.csr_off.upload(off, impl.stream);
        impl.disp.csr_partner.upload(h_part, impl.stream);
        impl.disp.csr_pairidx.upload(h_pidx, impl.stream);
        impl.disp.csr_is_i.upload(h_isi, impl.stream);
    }

    impl.disp_built_on_gpu = true;

    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::info(fmt::format(
            "  GPU D4 dispersion pair list: {} pairs built on GPU + CSR rebuilt for gather", h_count));
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

    // Claude Generated (Oct 2026): factorised C6 sums. In gradient mode (use_cn_final false) with the
    // per-step C6 refresh on and a symmetric reference table, only the per-atom half runs here and
    // k_dispersion_gather_uv forms C6 and dC6/dCN per edge; the per-pair arrays are then stale and are
    // filled before the next energy-only k_dispersion (prepareAndLaunchChargeIndependent).
    const bool fused = !use_cn_final && m_disp_c6_update && impl.c6_table_symmetric && impl.disp_fused_enabled
                       && !impl.dc6_tab_enabled
                       && impl.disp.csr_N == N;
    if (fused) {
        LaunchConfig cfg_uv = getLaunchConfig(N * impl.disp_ntypes, m_block_size);
        k_disp_uv<<<cfg_uv.gridSize, cfg_uv.blockSize, 0, impl.stream>>>(
            N, impl.disp_ntypes, impl.d_atom_types.ptr, impl.d_ctype_elem.ptr,
            impl.d_gw.ptr, impl.d_dgw.ptr, impl.d_c6_flat.ptr, impl.d_disp_uv.ptr);
        impl.disp_uv_valid = true;
        impl.disp_pair_c6_current = false;
    } else {
        launchDispC6PerPair(impl, N, m_disp_c6_update, m_block_size, impl.stream);
    }
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
    impl.disp_launched_early = false;   // Oct 2026: set again by launchDispersionGradientEarly() this step
    setGfnffFixModeHip(impl.fixed_point);   // WP10: before any capture or launch of this step
    g_rocm_small_eeq_stream = impl.single_stream ? stream : nullptr;   // Oct 2026, see its declaration
    g_rocm_small_eeq_aux = impl.single_stream ? impl.stream_threebody : nullptr;
    g_rocm_small_eeq_event = impl.single_stream ? impl.event_eeq_order : nullptr;

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
    if (impl.lean_step) {   // Claude Generated (Oct 2026): one kernel instead of four blits
        const int nz = std::max(3 * N, static_cast<int>(FFWorkspaceHipImpl::N_ENERGY_SLOTS));
        k_zero_step_accumulators<<<(nz + 255) / 256, 256, 0, stream>>>(
            impl.d_energy.ptr, impl.d_grad.ptr, 3 * N, impl.d_dEdcn.ptr, N,
            impl.d_energies.ptr, FFWorkspaceHipImpl::N_ENERGY_SLOTS);
    } else {
    hipMemsetAsync(impl.d_energy.ptr, 0, sizeof(double),         stream);
    hipMemsetAsync(impl.d_grad.ptr,   0, 3*N*sizeof(double),     stream);
    hipMemsetAsync(impl.d_dEdcn.ptr,  0, N*sizeof(double),       stream);
    hipMemsetAsync(impl.d_energies.ptr, 0,
                    FFWorkspaceHipImpl::N_ENERGY_SLOTS * sizeof(double), stream);
    }

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
    // Oct 2026 (Claude Generated): impl.single_stream also keeps everything on the main stream.
    const bool one_stream = force_single_stream || impl.single_stream;
    hipStream_t sA = one_stream ? stream : impl.stream_pairwise;
    hipStream_t sB = one_stream ? stream : impl.stream_bonded;
    hipStream_t sC = one_stream ? stream : impl.stream_threebody;

    // Ensure all uploads on main stream are visible to worker streams
    if (!one_stream) {
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
        // Claude Generated (Oct 2026): after fused gradient steps the per-pair C6 is stale (frozen-CN
        // energy calls reach this without a weight refresh); fill it from the current d_disp_uv first.
        if (!impl.disp_pair_c6_current && impl.disp_uv_valid)
            launchDispC6PerPair(impl, m_natoms, m_disp_c6_update, m_block_size, sA, /*with_uv=*/false);
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
        } else if (impl.rep_gather_enabled && impl.nb_rep_csr.N == N) {
            // Claude Generated (Oct 2026): per-atom gather, one workgroup of 64 lanes per atom (power of two
            // for the tree reduction); no per-pair atomicAdd (a CAS loop on gfx1200).
            k_repulsion_gather<<<std::max(1, N), 64, 0, sA>>>(
                impl.nb_rep_csr.N,
                impl.nb_rep_csr.off.ptr, impl.nb_rep_csr.partner.ptr, impl.nb_rep_csr.pairidx.ptr,
                impl.nonbonded_rep.alpha.ptr, impl.nonbonded_rep.repab.ptr,
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
        snapshotAccHip(impl.fixed_point, impl.d_grad_after_sA.ptr, impl.d_grad.ptr, 3*N, sA);
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
        snapshotAccHip(impl.fixed_point, impl.d_grad_after_bonds.ptr, impl.d_grad.ptr, 3*N, sB);
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
        snapshotAccHip(impl.fixed_point, impl.d_grad_after_angles.ptr, impl.d_grad.ptr, 3*N, sB);
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
        if (impl.fixed_point)   // WP10: k_bonds is the only producer (same stream): integers -> doubles
            fixToDoubleHip(impl.d_zz_hb.ptr, N, GFNFF_FIX_FORCE, sB);
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
        snapshotAccHip(impl.fixed_point, impl.d_grad_after_sB.ptr, impl.d_grad.ptr, 3*N, sB);
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
        // Claude Generated (Oct 2026): slot mode (per-triple slots + per-atom gather) unless disabled.
        const bool hb_slots = impl.hb_gather_enabled && impl.hb_slots.N == N && impl.hb_slots.nslots > 0;
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
            hb_slots ? impl.hb_slots.slot_off.ptr : nullptr,
            hb_slots ? impl.hb_slots.slot_g.ptr : nullptr,
            (impl.hb_slots.nq == impl.hbonds.n) ? impl.hb_slots.Qs_H.ptr : nullptr,
            (impl.hb_slots.nq == impl.hbonds.n) ? impl.hb_slots.Qs_A.ptr : nullptr,
            (impl.hb_slots.nq == impl.hbonds.n) ? impl.hb_slots.Qs_B.ptr : nullptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_HBOND]);
        if (hb_slots) {
            k_slot_gather<<<std::max(1, N), 64, 0, sC>>>(
                impl.hb_slots.N, impl.hb_slots.atom_off.ptr, impl.hb_slots.atom_slots.ptr,
                impl.hb_slots.slot_g.ptr, impl.d_grad.ptr);
        }
    }
    // Per-kernel snapshot: after all stream C kernels (Claude Generated April 2026)
    if (gradient && impl.d_grad_after_sC.ptr) {
        snapshotAccHip(impl.fixed_point, impl.d_grad_after_sC.ptr, impl.d_grad.ptr, 3*N, sC);
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
    updateHBChargeScaling(impl);   // Claude Generated (Oct 2026): Q_* follow the refreshed charges (synchronises)
    // Ensure the q arrays are fully written before the next step's HB (Phase 1) reads them.
    hipStreamSynchronize(impl.stream);
}

/// Claude Generated (Oct 2026): the gradient-mode dispersion gather of Phase 2 (former inline block of
/// launchChargeDependentAndFinish), also launched early by launchDispersionGradientEarly().
static void launchDispersionGatherHip(FFWorkspaceHipImpl& impl, hipStream_t sP)
{
    {
        // ROCm: per-atom gather (grad + dEdcn accumulated locally) instead of the pair kernel —
        // the deferred gradient-mode k_dispersion's ~8M FP64 atomicAdd was the dominant MD cost.
        // Fused path (Oct 2026, default): C6 and dC6/dCN per edge from d_disp_uv, no per-pair C6/dc6dcn
        // arrays. Otherwise (CURCUMA_ROCM_DISP_FUSED=0, no per-step C6 refresh, asymmetric C6 table,
        // stale d_disp_uv) the per-pair arrays are read via csr_pairidx. Both equal k_dispersion up to FP rounding.
        const int disp_block = 64;   // one workgroup per atom (power of two for the tree reduction)
        int grid = impl.disp.csr_N;
        if (impl.disp_uv_valid && !impl.disp_pair_c6_current) {
            // Claude Generated (Oct 2026): C6 and dC6/dCN formed per edge from d_disp_uv (fused path).
            k_dispersion_gather_uv<<<grid, disp_block, 0, sP>>>(
                impl.disp.csr_N, impl.disp_ntypes,
                impl.disp.csr_off.ptr, impl.disp.csr_partner.ptr, impl.disp.csr_pairidx.ptr,
                impl.d_atom_ctype.ptr, impl.d_disp_uv.ptr, impl.d_gw.ptr,
                impl.disp.r4r2ij.ptr, impl.disp.r0_sq.ptr, impl.disp.zetac6.ptr, impl.disp.r_cut.ptr,
                impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
                impl.d_dEdcn.ptr, impl.d_grad.ptr, &impl.d_energies.ptr[impl.E_DISP]);
        } else {
        k_dispersion_gather<<<grid, disp_block, 0, sP>>>(
            impl.disp.csr_N,
            impl.disp.csr_off.ptr, impl.disp.csr_partner.ptr,
            impl.disp.csr_pairidx.ptr, impl.disp.csr_is_i.ptr,
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
        }
        checkCuda(hipGetLastError(), "k_dispersion_gather launch");
    }
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
        if (impl.lean_step && impl.h_charges_stage && impl.d_charges.n >= N) {
            // Claude Generated (Oct 2026): pinned stage + async copy instead of a blocking pageable upload.
            std::memcpy(impl.h_charges_stage, m_eeq_charges.data(), N * sizeof(double));
            checkCuda(hipMemcpyAsync(impl.d_charges.ptr, impl.h_charges_stage, N * sizeof(double), hipMemcpyHostToDevice,
                                     stream), "staged EEQ charge upload");
        } else {
            impl.d_charges.upload(m_eeq_charges.data(), N, stream);
        }
    }

    // =========================================================================
    // 2. Launch deferred kernels on stream A
    //    Wait for charge-independent pairwise kernels to finish first,
    //    since these kernels write to the same d_grad/d_energies buffers.
    // =========================================================================
    // Record event on main stream (after charge upload) so pairwise stream can wait.
    // WP4: dedicated phase-2 event (event_p2_upload) to keep capture isolated from
    // phase-1's event_upload, which is recorded inside the phase-1 graph.
    // Oct 2026 (Claude Generated): with impl.single_stream the pairwise kernels run on the main stream
    // and the cross-stream events of this phase are left out (stream order is the dependency).
    const bool one_stream = impl.single_stream;
    hipStream_t sP = one_stream ? stream : impl.stream_pairwise;
    if (!one_stream) {
        hipEventRecord(impl.event_p2_upload, stream);
        hipStreamWaitEvent(impl.stream_pairwise, impl.event_p2_upload, 0);
    }

    // k_dispersion (gradient mode): deferred from Phase 1.
    // computeGaussianWeightsOnGPU() was called after Phase 1 (using d_cn, current step)
    // and runs on the main stream. The event_upload fence above ensures stream_pairwise
    // sees the completed dc6dcn values before k_dispersion launches.
    if (m_dispersion_enabled && gradient && !impl.disp_launched_early)
        launchDispersionGatherHip(impl, sP);
    impl.disp_launched_early = false;

    // k_coulomb (needs EEQ charges)
    if (m_coulomb_enabled && impl.coulomb_implicit) {
        launchCoulombTilesHip(m_natoms, impl.d_coul_alp.ptr,
            impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
            impl.d_charges.ptr, impl.coulomb_rcut,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_COUL], sP);
        checkCuda(hipGetLastError(), "k_coulomb_tiles_hip (implicit) launch");
    } else if (m_coulomb_enabled && impl.coulomb.n > 0) {
        // ROCm: the pair k_coulomb's ~6M FP64 atomicAdd dominated the MD step. Preferred path is
        // the shared-memory tiled all-pairs kernel (k_coulomb_tiles_hip, each pair once) which also fixes
        // the gather's random-memory bottleneck; it needs the per-atom charge-corrected EEQ alpha
        // (== the source of the precomputed gamma_ij) to recompute gamma on the fly. If the EEQ
        // topology alpha is not resident, fall back to the CSR gather (still atomic-light).
        const int coul_block = 128;
        // The dense kernel loops ALL j with a uniform runtime cutoff, so it is only equivalent to
        // the pair list when that list IS all-pairs (default eeq_distance_cutoff=0). When a spatial
        // EEQ cutoff is active the list is a FIXED param-gen-geometry subset, and a runtime all-j
        // loop would diverge from CPU during MD/opt (newly-close pairs) — fall back to the gather
        // (exact CSR list) there. csr_N*(csr_N-1)/2 == n iff the list is the full upper triangle.
        const long all_pairs = static_cast<long>(impl.coulomb.csr_N) * (impl.coulomb.csr_N - 1) / 2;
        const bool coulomb_all_pairs = (static_cast<long>(impl.coulomb.n) == all_pairs);
        if (impl.eeq_topo.valid && impl.eeq_topo.d_alpha_corrected.ptr && coulomb_all_pairs) {
            launchCoulombTilesHip(impl.coulomb.csr_N,
                impl.eeq_topo.d_alpha_corrected.ptr,
                impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
                impl.d_charges.ptr, impl.coulomb.r_cut_scalar,
                impl.d_grad.ptr,
                &impl.d_energies.ptr[impl.E_COUL], sP);
            checkCuda(hipGetLastError(), "k_coulomb_tiles_hip launch");
        } else {
            k_coulomb_gather<<<impl.coulomb.csr_N, coul_block, 0, sP>>>(
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
    if (!one_stream) hipEventRecord(impl.event_p2_pairwise, impl.stream_pairwise);

    // =========================================================================
    // 3. Postprocess on main stream
    // Phase 4b+Fusion: Fused Coulomb postprocess (Claude Generated March 2026)
    //   k_coulomb_postprocess — fused self-energy + qtmp subtraction in single O(N) pass
    //     Needs d_dEdcn (from pairwise + bonded) for qtmp, NOT threebody
    //   k_cn_chainrule   — needs d_dEdcn complete from ALL streams → wait for all 3
    // =========================================================================

    const bool need_snapshots = (m_verbosity >= 3);
    // WP10: with fixed_point, d_dEdcn / d_grad / d_energies are integer accumulators (as CUDA): dEdcn is converted
    // once every producer is done (before its first double consumer), d_grad after the CN chain rule (the last
    // accumulation), the energy slots before the download.
    bool dEdcn_converted = !impl.fixed_point;
    auto convert_dEdcn = [&] {
        if (dEdcn_converted) return;
        fixToDoubleHip(impl.d_dEdcn.ptr, N, GFNFF_FIX_FORCE, stream);
        dEdcn_converted = true;
    };
    if (m_coulomb_enabled && impl.coul_self_on_gpu) {
        // Wait for pairwise+bonded streams (dEdcn dependency for qtmp subtraction).
        // event_p2_pairwise: Phase-2-dedicated — always wait (records k_coulomb completion).
        // event_bonded: Phase-1-internal — skip if Phase-1 ran as graph (already sequenced).
        if (!one_stream) {
            hipStreamWaitEvent(stream, impl.event_p2_pairwise, 0);
            if (!p1_was_graph)
                hipStreamWaitEvent(stream, impl.event_bonded, 0);
        }

        convert_dEdcn();
        // d_dEdcn_snapshot must always be populated — dEdcnTotal() is public API
        // used by external code (tests, chain-rule validation) regardless of verbosity.
        if (m_coulomb_enabled && impl.coul_self_on_gpu && impl.d_dEdcn_snapshot.ptr) {
            hipMemcpyAsync(impl.d_dEdcn_snapshot.ptr, impl.d_dEdcn.ptr,
                            N * sizeof(double), hipMemcpyDeviceToDevice, stream);
        }

        if (need_snapshots) {
            // Grad snapshot is only for diagnostics at verbosity >= 3
            snapshotAccHip(impl.fixed_point, impl.d_grad_snapshot.ptr, impl.d_grad.ptr, 3*N, stream);
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
        if (!one_stream) {
            hipStreamWaitEvent(stream, impl.event_p2_pairwise, 0);
            if (!p1_was_graph) {
                hipStreamWaitEvent(stream, impl.event_bonded, 0);
                hipStreamWaitEvent(stream, impl.event_threebody, 0);
            }
        }

        convert_dEdcn();
        if (need_snapshots) {
            snapshotAccHip(impl.fixed_point, impl.d_grad_snapshot.ptr, impl.d_grad.ptr, 3*N, stream);
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
    if (!one_stream) {
        hipStreamWaitEvent(stream, impl.event_p2_pairwise, 0);
        if (!p1_was_graph) {
            hipStreamWaitEvent(stream, impl.event_bonded, 0);
            hipStreamWaitEvent(stream, impl.event_threebody, 0);
        }
    }

    if (impl.fixed_point) {
        convert_dEdcn();   // no-op unless neither the Coulomb postprocess nor the chain rule ran
        if (gradient)      // last accumulation is done: fixed point -> double before the download
            fixToDoubleHip(impl.d_grad.ptr, 3 * N, GFNFF_FIX_FORCE, stream);
        static_assert(FFWorkspaceHipImpl::E_DISP == 0, "slot 0 is the dispersion slot (own fixed-point scale)");
        fixToDoubleHip(impl.d_energies.ptr, 1, GFNFF_FIX_ENERGY_DISP, stream);
        fixToDoubleHip(impl.d_energies.ptr + 1, FFWorkspaceHipImpl::N_ENERGY_SLOTS - 1, GFNFF_FIX_ENERGY, stream);
        if (impl.h_fix_bad)   // range / finiteness flag of the accumulation, read after the final sync
            checkCuda(hipMemcpyFromSymbolAsync(impl.h_fix_bad, HIP_SYMBOL(d_gfnff_fix_bad), sizeof(int), 0,
                                               hipMemcpyDeviceToHost, stream), "fixed-point flag read-back");
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
    m_dc6dcn_skip_pending = !impl.lean_step && (impl.h_dc6dcn_skip != nullptr && *impl.h_dc6dcn_skip != 0);

    // WP10: a non-finite or out-of-range contribution cannot be represented in the fixed-point accumulators; the
    // FP64 path would return NaN or a huge value there, so the step returns NaN (energy and gradient) instead of a
    // wrapped number. The flag is sticky on the device and reset here.
    bool fix_bad = false;
    if (impl.fixed_point && impl.h_fix_bad && *impl.h_fix_bad != 0) {
        fix_bad = true;
        *impl.h_fix_bad = 0;
        const int zero = 0;
        checkCuda(hipMemcpyToSymbol(HIP_SYMBOL(d_gfnff_fix_bad), &zero, sizeof(int)), "reset fixed-point flag");
        if (impl.fix_bad_steps++ == 0)
            CurcumaLogger::warn("ROCm GFN-FF: a force, energy or CN contribution was not finite or outside the "
                                "fixed-point range (-gfnff.rocm_fixed_point); the step returns NaN");
    }

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

    if (fix_bad) {
        if (gradient) m_result_gradient.setConstant(std::numeric_limits<double>::quiet_NaN());
        return std::numeric_limits<double>::quiet_NaN();
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
    impl.cn_stale_cached = -1;   // Claude Generated (Oct 2026): flag computed ahead belongs to the old geometry

    // Fill pinned staging buffers in-place, then async DMA to GPU (SoA layout).
    // Claude Generated (Phase 10, March 2026): SoA pinned buffers for coalesced warp reads.
    for (int i = 0; i < N; ++i) {
        m_h_x[i] = geom(i, 0);
        m_h_y[i] = geom(i, 1);
        m_h_z[i] = geom(i, 2);
    }
    impl.h_geom_ok = true;
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
    // Distance cutoff of the O(N^2) CN sum: 30 Bohr. Sep 2026: dead rcov-based estimate
    // (overwritten right after, with a CODATA-1986 factor on radii already in Bohr) removed,
    // mirroring the CUDA file - Claude Generated; uncompiled here (no ROCm SDK).
    const double threshold_sq = 900.0;

    setGfnffFixModeHip(impl.fixed_point);   // WP10
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
        if (impl.fixed_point)   // WP10: integer CN sums -> doubles before the log squashing
            fixToDoubleHip(impl.d_cn_raw.ptr, N, GFNFF_FIX_CN, impl.stream);

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
    // Oct 2026 (Claude Generated): skipped on the lean path - nothing reads the flag
    // (FFWorkspaceGPU::shouldSkipDc6dcn() has no caller); CURCUMA_ROCM_LEAN_STEP=0 restores it.
    if (!impl.lean_step) {
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
    if (impl.single_stream && impl.event_cn_ready) {   // Claude Generated (Oct 2026), see finalizeCNForCPU()
        checkCuda(hipEventRecord(impl.event_cn_ready, impl.stream), "record cn_ready");
        impl.cn_event_recorded = true;
    }

    // Download CN_final, CN_raw, and dlogdcn to pre-allocated pinned buffers
    // Claude Generated (March 2026): Avoids heap corruption from CUDA allocator
    // CN_raw needed for correct dlogdcn computation (chain-rule gradient)
    // Oct 2026 (Claude Generated): on the lean path finalizeCNForCPU() enqueues these copies, the only
    // place that reads the host buffers (the normal device-EEQ step never does).
    if (!impl.lean_step) {
        hipMemcpyAsync(m_h_cn_final, impl.d_cn_final.ptr, N * sizeof(double),
                        hipMemcpyDeviceToHost, impl.stream);
        hipMemcpyAsync(m_h_cn_raw, impl.d_cn_raw.ptr, N * sizeof(double),
                        hipMemcpyDeviceToHost, impl.stream);
        hipMemcpyAsync(m_h_dlogdcn, impl.d_dlogdcn.ptr, N * sizeof(double),
                        hipMemcpyDeviceToHost, impl.stream);
    }
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
    auto& impl = *m_impl;
    if (impl.lean_step && impl.d_cn_final.n >= N && impl.d_cn_raw.n >= N && impl.d_dlogdcn.n >= N) {
        // Claude Generated (Oct 2026): the copies computeCN() enqueues on the former path. With the event of
        // computeCN() (single-stream mode) they run on the otherwise idle stream_threebody right behind
        // k_dlogdcn, so the host (EEQ on the CPU) does not wait for the Phase-1 kernels queued after it;
        // d_cn_final, d_cn_raw and d_dlogdcn are not written again before the next computeCN().
        const bool early = impl.cn_event_recorded && impl.single_stream && impl.stream_threebody;
        hipStream_t s = early ? impl.stream_threebody : impl.stream;
        if (early) checkCuda(hipStreamWaitEvent(s, impl.event_cn_ready, 0), "finalizeCNForCPU: wait cn_ready");
        checkCuda(hipMemcpyAsync(m_h_cn_final, impl.d_cn_final.ptr, N * sizeof(double), hipMemcpyDeviceToHost, s),
                  "finalizeCNForCPU: cn_final");
        checkCuda(hipMemcpyAsync(m_h_cn_raw, impl.d_cn_raw.ptr, N * sizeof(double), hipMemcpyDeviceToHost, s),
                  "finalizeCNForCPU: cn_raw");
        checkCuda(hipMemcpyAsync(m_h_dlogdcn, impl.d_dlogdcn.ptr, N * sizeof(double), hipMemcpyDeviceToHost, s),
                  "finalizeCNForCPU: dlogdcn");
        checkCuda(hipStreamSynchronize(s), "finalizeCNForCPU: sync");
        if (out_cn.size() != N) out_cn.resize(N);
        std::memcpy(out_cn.data(), m_h_cn_final, N * sizeof(double));
        m_cn_computed = true;
        return;
    }
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

/// K3 is not used on ROCm (Claude Generated, Oct 2026): running the device EEQ on its own stream behind an
/// input event was measured slower (hardware-queue switches); no event, the solve keeps its null-stream order.
void FFWorkspaceHip::recordEEQInputsReady() {}
void* FFWorkspaceHip::eeqInputsReadyEvent() const { return nullptr; }

/// K2 on ROCm (Claude Generated, Oct 2026), single-stream mode only: the gradient-mode dispersion gather
/// reads geometry, dc6dcn / d_disp_uv and the Phase-1 topology charges, never the EEQ charges, so it is
/// enqueued right after the Gaussian weights; while the host solves the EEQ (many-fragment route) the GPU
/// works on it. Same stream, same order relative to the other gradient writers as in Phase 2, so the
/// result is unchanged. Returns false (Phase 2 launches it) with worker streams or without dispersion.
bool FFWorkspaceHip::launchDispersionGradientEarly()
{
    if (!m_impl || !m_impl->single_stream || !m_dispersion_enabled) return false;
    launchDispersionGatherHip(*m_impl, m_impl->stream);
    m_impl->disp_launched_early = true;
    return true;
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

    // Claude Generated (Oct 2026): lean path - one launch for this test and, when cnPairListStale() has
    // been used before, its test against the CN pair-list reference; one pinned read-back.
    if (impl.lean_step && impl.single_stream && impl.h_ref_ok && impl.h_geom_ok) {
        // Claude Generated (Oct 2026): host test against the host mirrors, no device round trip.
        const bool want_cn = impl.cn_stale_thr >= 0.0 && !impl.cn_ref_coords.empty() && impl.h_cnref_ok;
        impl.cn_stale_cached = want_cn
            ? (hostDisplacementExceeded(m_natoms, m_h_x, m_h_y, m_h_z, impl.h_cnref_xyz,
                                        impl.cn_stale_thr * impl.cn_stale_thr) ? 1 : 0)
            : -1;
        return hostDisplacementExceeded(m_natoms, m_h_x, m_h_y, m_h_z, impl.h_ref_xyz, threshold * threshold);
    }
    if (impl.lean_step && impl.h_disp_flags) {
        const bool want_cn = impl.cn_stale_thr >= 0.0 && !impl.cn_ref_coords.empty();
        k_check_displacement_pair<<<1, 1024, 0, impl.stream>>>(
            m_natoms, impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
            impl.ref_coords.d_rx.ptr, impl.ref_coords.d_ry.ptr, impl.ref_coords.d_rz.ptr, threshold * threshold,
            want_cn ? impl.cn_ref_coords.d_rx.ptr : nullptr, want_cn ? impl.cn_ref_coords.d_ry.ptr : nullptr,
            want_cn ? impl.cn_ref_coords.d_rz.ptr : nullptr, impl.cn_stale_thr * impl.cn_stale_thr,
            impl.d_disp_flags2.ptr);
        checkCuda(hipMemcpyAsync(impl.h_disp_flags, impl.d_disp_flags2.ptr, 2 * sizeof(int),
                                 hipMemcpyDeviceToHost, impl.stream), "checkDisplacement: download flags");
        checkCuda(hipStreamSynchronize(impl.stream), "checkDisplacement: sync");
        impl.cn_stale_cached = want_cn ? impl.h_disp_flags[1] : -1;
        return impl.h_disp_flags[0] != 0;
    }

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
    impl.h_ref_ok = copyHostGeometryMirror(impl.h_geom_ok, m_natoms, m_h_x, m_h_y, m_h_z, impl.h_ref_xyz);

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

#endif // USE_ROCM
