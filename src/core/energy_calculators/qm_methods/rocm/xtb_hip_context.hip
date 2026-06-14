/*
 * <Native xTB ROCm/HIP Context — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): native GFN1/GFN2 ROCm device engine.
 *   Stage 0: device handshake (query + stream).
 *   Stage 1: solveGeneralized — the per-iteration F C = S C ε via rocSOLVER dsygvd.
 *   Stage 2: device-resident GFN1 SCF (residentBegin/Solve/Density/Finalize) — Fock
 *            build + populations via the __global__ kernels below, density via rocBLAS,
 *            eigensolve via rocSOLVER; H0/S stay resident.
 * Compiled by hipcc into a plain object (CMake add_custom_command), then linked into
 * curcuma_core with g++ — so the device kernels build, but the HIP-only --offload-arch /
 * --hip-link flags never reach the GNU link (which would drop libgomp). Links libamdhip64
 * + librocsolver + librocblas.
 */

#ifdef USE_ROCM_XTB

#include "xtb_hip_context.h"

#include <hip/hip_runtime.h>
#ifdef HAVE_ROCSOLVER
#include <rocblas/rocblas.h>
#include <rocsolver/rocsolver.h>
#include "xtb_hip_integrals.hiph"                 // CN/overlap/H0/gamma device helpers
#include "../parameters/xtb_params_extra.hpp"     // covalent_rad_d3_au / pauling_en / atomic_rad_au
#endif

#include <memory>
#include <string>
#include <vector>

namespace curcuma {
namespace xtb {
namespace gpu {

#ifdef HAVE_ROCSOLVER
namespace {
// GFN1 isotropic Fock: F[i,j] = H0[i,j] − ½·S[i,j]·(v[i] + v[j]). Column-major n×n.
__global__ void k_fock(const double* H0, const double* S, const double* v, double* F, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;  // row
    const int j = blockIdx.y * blockDim.y + threadIdx.y;  // col
    if (i >= n || j >= n) return;
    const int idx = i + j * n;
    F[idx] = H0[idx] - 0.5 * S[idx] * (v[i] + v[j]);
}
// Column scaling: B[i,j] = A[i,j]·d[j] (forms Cw = C·diag(occ)).
__global__ void k_scale_cols(const double* A, const double* d, double* B, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= n || j >= n) return;
    const int idx = i + j * n;
    B[idx] = A[idx] * d[j];
}
// Mulliken populations + band-energy partials from the resident density P (symmetric):
//   out[i]   = Σ_j P[i,j]·S[i,j],   out[n+i] = Σ_j P[i,j]·H0[i,j].
__global__ void k_popband(const double* P, const double* S, const double* H0, double* out, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double pop = 0.0, bnd = 0.0;
    for (int j = 0; j < n; ++j) {
        const double p = P[i + j * n];
        pop += p * S[i + j * n];
        bnd += p * H0[i + j * n];
    }
    out[i] = pop; out[n + i] = bnd;
}

// ---- Stage 3 integral build kernels (ported from cuda/xtb_gpu_context.cu) ----
// Element parameter tables (covrad/pauling/arad, length 86) are passed as device
// arrays rather than __constant__ symbols.

// Coordination numbers: one thread per atom (cutoff 25 bohr).
__global__ void k_cn(int nat, const double* __restrict__ xyz, const int* __restrict__ z,
                     int is_gfn2, const double* __restrict__ covrad, double* __restrict__ cn)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    const double xi = xyz[3 * i + 0], yi = xyz[3 * i + 1], zi = xyz[3 * i + 2];
    const double rci = covrad[z[i] - 1];
    const double cutoff2 = 25.0 * 25.0;
    double sum = 0.0;
    for (int j = 0; j < nat; ++j) {
        if (j == i) continue;
        const double dx = xi - xyz[3 * j + 0], dy = yi - xyz[3 * j + 1], dz = zi - xyz[3 * j + 2];
        const double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 > cutoff2 || r2 < 1.0e-12) continue;
        const double r  = sqrt(r2);
        const double rc = rci + covrad[z[j] - 1];
        sum += d_cn_pair(r, rc, is_gfn2 != 0);
    }
    cn[i] = sum;
}

// CN-shifted shell self-energies: se[s] = selfenergy[s] − kcn[s]·CN[sh2at[s]].
__global__ void k_self_energy(int nsh, const double* __restrict__ selfenergy,
                              const double* __restrict__ kcn, const int* __restrict__ sh2at,
                              const double* __restrict__ cn, double* __restrict__ se)
{
    const int s = blockIdx.x * blockDim.x + threadIdx.x;
    if (s >= nsh) return;
    se[s] = selfenergy[s] - kcn[s] * cn[sh2at[s]];
}

// Overlap S + bare Hamiltonian H0: one thread per shell-pair (a,b).
__global__ void k_overlap_h0(
    int nsh, int nao, int is_gfn2,
    const int* __restrict__ sh2at, const int* __restrict__ ang_sh,
    const int* __restrict__ iao_sh, const int* __restrict__ nao_sh,
    const int* __restrict__ sh_nprim, const int* __restrict__ sh_prim_off,
    const double* __restrict__ prim_alpha, const double* __restrict__ prim_coeff,
    const double* __restrict__ sh_zeta, const double* __restrict__ shpoly,
    const double* __restrict__ se, const int* __restrict__ z,
    const int* __restrict__ valence, const double* __restrict__ xyz,
    const double* __restrict__ atomic_rad, const double* __restrict__ pauling,
    double* __restrict__ S, double* __restrict__ H0)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    const int b = blockIdx.y * blockDim.y + threadIdx.y;
    if (a >= nsh || b >= nsh) return;

    const int iat = sh2at[a], jat = sh2at[b];
    const int la = ang_sh[a], lb = ang_sh[b];
    const double xa = xyz[3 * iat + 0], ya = xyz[3 * iat + 1], za = xyz[3 * iat + 2];
    const double xb = xyz[3 * jat + 0], yb = xyz[3 * jat + 1], zb = xyz[3 * jat + 2];
    const double avg_eps = 0.5 * (se[a] + se[b]);

    double h_factor;
    if (iat == jat) {
        h_factor = 1.0;
    } else {
        const int zi = z[iat], zj = z[jat];
        const double dx = xa - xb, dy = ya - yb, dz = za - zb;
        const double r2 = dx * dx + dy * dy + dz * dz;
        const double rr = sqrt(sqrt(r2) / (atomic_rad[zi - 1] + atomic_rad[zj - 1]));
        const double pi_ij = (1.0 + shpoly[a] * rr) * (1.0 + shpoly[b] * rr);
        double hs;
        if (is_gfn2 == 0) {
            const bool vi = valence[a] != 0, vj = valence[b] != 0;
            if (vi && vj) {
                double den = pauling[zi - 1] - pauling[zj - 1]; den *= den;
                hs = d_kpair_gfn1(zi, zj) * d_kshell_gfn1(la, lb) * (1.0 + (-7.0e-3) * den);
            } else if (vi && !vj) {
                hs = 0.5 * (d_kshell_gfn1(la, la) + 2.85);
            } else if (!vi && vj) {
                hs = 0.5 * (d_kshell_gfn1(lb, lb) + 2.85);
            } else {
                hs = 2.85;
            }
        } else {
            double den = pauling[zi - 1] - pauling[zj - 1]; den *= den;
            const double enp = 1.0 + 2.0e-2 * den;
            const double km  = d_kshell_gfn2(la, lb) * enp;
            const double za_ = sh_zeta[a], zb_ = sh_zeta[b];
            const double zij = pow(2.0 * sqrt(za_ * zb_) / (za_ + zb_), 0.5);
            hs = zij * km;
        }
        h_factor = hs * pi_ij;
    }

    const int ia_start = iao_sh[a], ia_nao = nao_sh[a];
    const int jb_start = iao_sh[b], jb_nao = nao_sh[b];
    const double* pa_alpha = prim_alpha + sh_prim_off[a];
    const double* pa_coeff = prim_coeff + sh_prim_off[a];
    const int npa = sh_nprim[a];
    const double* pb_alpha = prim_alpha + sh_prim_off[b];
    const double* pb_coeff = prim_coeff + sh_prim_off[b];
    const int npb = sh_nprim[b];

    for (int ia = 0; ia < ia_nao; ++ia) {
        const int mu = ia_start + ia;
        const int ta = d_ao_to_type(la, ia);
        if (ta < 0) continue;
        for (int jb = 0; jb < jb_nao; ++jb) {
            const int nu = jb_start + jb;
            const int tb = d_ao_to_type(lb, jb);
            if (tb < 0) continue;
            const double s_ab = (iat == jat && a == b && ta == tb)
                ? 1.0
                : d_cgto_overlap(pa_alpha, pa_coeff, npa, pb_alpha, pb_coeff, npb,
                                 xa, ya, za, xb, yb, zb, ta, tb);
            const size_t idx = static_cast<size_t>(mu) + static_cast<size_t>(nu) * nao;
            S[idx]  = s_ab;
            H0[idx] = avg_eps * h_factor * s_ab;
        }
    }
}

// Coulomb γ matrix: one thread per shell-pair (gexp=2; per-shell hardness g host-built).
__global__ void k_gamma(int nsh, int is_gfn2, const int* __restrict__ sh2at,
                        const double* __restrict__ g, const double* __restrict__ xyz,
                        double* __restrict__ gamma)
{
    const int is = blockIdx.x * blockDim.x + threadIdx.x;
    const int js = blockIdx.y * blockDim.y + threadIdx.y;
    if (is >= nsh || js >= nsh) return;
    const int iat = sh2at[is], jat = sh2at[js];
    double val;
    if (is == js) {
        val = g[is];
    } else if (iat == jat) {
        val = d_coulomb_average(g[is], g[js], is_gfn2 != 0);
    } else {
        const double dx = xyz[3 * iat + 0] - xyz[3 * jat + 0];
        const double dy = xyz[3 * iat + 1] - xyz[3 * jat + 1];
        const double dz = xyz[3 * iat + 2] - xyz[3 * jat + 2];
        const double r1g = dx * dx + dy * dy + dz * dz;          // gexp=2 → R²
        const double gam = d_coulomb_average(g[is], g[js], is_gfn2 != 0);
        val = pow(r1g + pow(gam, -2.0), -0.5);
    }
    gamma[static_cast<size_t>(is) + static_cast<size_t>(js) * nsh] = val;
}
} // namespace
#endif // HAVE_ROCSOLVER

struct XtbHipContext::Impl {
    bool        ok = false;
    int         device = -1;
    std::string name;
    hipStream_t stream = nullptr;

#ifdef HAVE_ROCSOLVER
    rocblas_handle handle = nullptr;   // rocSOLVER reuses the rocBLAS handle
    bool    handle_ok = false;
    int     solver_n  = 0;             // cached device-buffer size
    double* dF   = nullptr;            // F → eigenvectors (n×n)
    double* dS   = nullptr;            // S (overwritten by its Cholesky factor)
    double* dEps = nullptr;            // eigenvalues (n)
    double* dE   = nullptr;            // rocSOLVER workspace (n)
    rocblas_int* dInfo = nullptr;      // convergence flag

    bool ensureSolver(int n) {
        if (!handle_ok) {
            if (rocblas_create_handle(&handle) != rocblas_status_success) return false;
            rocblas_set_stream(handle, stream);
            handle_ok = true;
        }
        if (n == solver_n && dF) return true;
        freeSolver();
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMalloc(&dF, nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dS, nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dEps, n * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dE, n * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dInfo, sizeof(rocblas_int)) != hipSuccess) return false;
        solver_n = n;
        return true;
    }
    void freeSolver() {
        for (void* p : { (void*)dF, (void*)dS, (void*)dEps, (void*)dE, (void*)dInfo })
            if (p) hipFree(p);
        dF = dS = dEps = dE = nullptr; dInfo = nullptr; solver_n = 0;
    }

    // ----- resident GFN1 SCF (Stage 2) ------------------------------------
    // Reuses dF (eigenvectors C), dS (Cholesky work copy), dEps/dE/dInfo from the
    // solver cache; adds the geometry-resident H0/S + density work buffers.
    int     res_n = 0;
    double* dH0 = nullptr;   // resident H0 (n×n)
    double* dS0 = nullptr;   // resident S  (n×n, untouched; dS is the work copy)
    double* dvao = nullptr;  // v_ao (n)
    double* dP  = nullptr;   // density (n×n)
    double* dCw = nullptr;   // C·diag(occ) (n×n)
    double* dOcc = nullptr;  // occupations (n)
    double* dPB = nullptr;   // pop + bandpart (2n)

    bool ensureResident(int n) {
        if (!ensureSolver(n)) return false;
        if (n == res_n && dH0) return true;
        freeResident();
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMalloc(&dH0,  nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dS0,  nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dvao, n  * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dP,   nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dCw,  nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dOcc, n  * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dPB,  2 * n * sizeof(double)) != hipSuccess) return false;
        res_n = n;
        return true;
    }
    void freeResident() {
        for (void* p : { (void*)dH0, (void*)dS0, (void*)dvao, (void*)dP, (void*)dCw, (void*)dOcc, (void*)dPB })
            if (p) hipFree(p);
        dH0 = dS0 = dvao = dP = dCw = dOcc = dPB = nullptr; res_n = 0;
    }

    bool resBegin(const double* H0, const double* S, int n) {
        if (!ok || n <= 0 || !ensureResident(n)) return false;
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMemcpy(dH0, H0, nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dS0, S,  nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        return true;
    }
    bool resSolve(const double* v_ao, double* eps_out) {
        const int n = res_n; if (n <= 0) return false;
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMemcpy(dvao, v_ao, n * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        dim3 blk(16, 16), grd((n + 15) / 16, (n + 15) / 16);
        hipLaunchKernelGGL(k_fock, grd, blk, 0, stream, dH0, dS0, dvao, dF, n);   // F → dF
        if (hipMemcpyAsync(dS, dS0, nn * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess) return false;
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        // F C = S C ε on the device; dF ← eigenvectors (ascending), dS destroyed.
        if (rocsolver_dsygvd(handle, rocblas_eform_ax, rocblas_evect_original, rocblas_fill_lower,
                             n, dF, n, dS, n, dEps, dE, dInfo) != rocblas_status_success) return false;
        if (hipDeviceSynchronize() != hipSuccess) return false;
        rocblas_int info = 0;
        if (hipMemcpy(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (info != 0) return false;
        if (hipMemcpy(eps_out, dEps, n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return true;   // dF holds C resident, ascending order
    }
    bool resDensity(const double* occ, int ncol, double* pop_ao, double* band) {
        const int n = res_n; if (n <= 0 || ncol < 0 || ncol > n) return false;
        std::vector<double> occ_full(n, 0.0);          // ascending order matches C columns
        for (int r = 0; r < ncol; ++r) occ_full[r] = occ[r];
        if (hipMemcpy(dOcc, occ_full.data(), n * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        dim3 blk(16, 16), grd((n + 15) / 16, (n + 15) / 16);
        hipLaunchKernelGGL(k_scale_cols, grd, blk, 0, stream, dF, dOcc, dCw, n);   // Cw = C·diag(occ)
        const double one = 1.0, zero = 0.0;
        // P = Cw · Cᵀ
        if (rocblas_dgemm(handle, rocblas_operation_none, rocblas_operation_transpose,
                          n, n, n, &one, dCw, n, dF, n, &zero, dP, n) != rocblas_status_success) return false;
        dim3 blk1(64), grd1((n + 63) / 64);
        hipLaunchKernelGGL(k_popband, grd1, blk1, 0, stream, dP, dS0, dH0, dPB, n);
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        std::vector<double> pb(static_cast<size_t>(2) * n);
        if (hipMemcpy(pb.data(), dPB, 2 * n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        double b = 0.0;
        for (int i = 0; i < n; ++i) { pop_ao[i] = pb[i]; b += pb[n + i]; }
        *band = b;
        return true;
    }
    bool resFinalize(double* P, double* C) {
        const int n = res_n; if (n <= 0) return false;
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMemcpy(P, dP, nn * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (hipMemcpy(C, dF, nn * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return true;
    }

    // ----- Stage 3: device integral build ---------------------------------
    bool basis_ready = false;
    int  bnat = 0, bnsh = 0, bnao = 0, bnprim = 0, bis_gfn2 = 0, bgeom_nat = 0;
    int*    dz = nullptr; int* dsh2at = nullptr; int* dang_sh = nullptr; int* diao_sh = nullptr;
    int*    dnao_sh = nullptr; int* dsh_nprim = nullptr; int* dsh_prim_off = nullptr; int* dvalence = nullptr;
    double* dprim_alpha = nullptr; double* dprim_coeff = nullptr; double* dsh_zeta = nullptr;
    double* dselfenergy = nullptr; double* dkcn = nullptr; double* dshpoly = nullptr; double* dghard = nullptr;
    double* dcovrad = nullptr; double* dpauling = nullptr; double* darad = nullptr;   // element tables (86)
    double* dxyz = nullptr; double* dcn = nullptr; double* dse = nullptr; double* dgamma = nullptr; double* dL = nullptr;

    template <class T> bool up(T*& dptr, const T* host, size_t count) {
        if (dptr) hipFree(dptr);
        dptr = nullptr;
        if (hipMalloc(&dptr, count * sizeof(T)) != hipSuccess) return false;
        return hipMemcpy(dptr, host, count * sizeof(T), hipMemcpyHostToDevice) == hipSuccess;
    }

    bool begBasis(const XtbHipBasisData& bd) {
        if (!ok || bd.nao <= 0 || !ensureResident(bd.nao)) return false;  // allocates dH0/dS0/...
        bnat = bd.nat; bnsh = bd.nsh; bnao = bd.nao; bnprim = bd.nprim_total; bis_gfn2 = bd.is_gfn2;
        bool g = true;
        g &= up(dz, bd.z, bnat);
        g &= up(dsh2at, bd.sh2at, bnsh);    g &= up(dang_sh, bd.ang_sh, bnsh);
        g &= up(diao_sh, bd.iao_sh, bnsh);  g &= up(dnao_sh, bd.nao_sh, bnsh);
        g &= up(dsh_nprim, bd.sh_nprim, bnsh); g &= up(dsh_prim_off, bd.sh_prim_off, bnsh);
        if (bd.valence) g &= up(dvalence, bd.valence, bnsh);   // GFN1 only; null for GFN2
        else { if (dvalence) hipFree(dvalence); dvalence = nullptr; }
        g &= up(dprim_alpha, bd.prim_alpha, bnprim); g &= up(dprim_coeff, bd.prim_coeff, bnprim);
        g &= up(dsh_zeta, bd.sh_zeta, bnsh); g &= up(dselfenergy, bd.selfenergy, bnsh);
        g &= up(dkcn, bd.kcn, bnsh);         g &= up(dshpoly, bd.shpoly, bnsh);
        g &= up(dghard, bd.shell_hardness, bnsh);
        double covrad[86], pauling[86], arad[86];
        for (int z = 1; z <= 86; ++z) {
            covrad[z - 1]  = curcuma::xtb::covalent_rad_d3_au(z);
            pauling[z - 1] = curcuma::xtb::pauling_en[z - 1];
            arad[z - 1]    = curcuma::xtb::atomic_rad_au(z);
        }
        g &= up(dcovrad, covrad, 86); g &= up(dpauling, pauling, 86); g &= up(darad, arad, 86);
        if (bnat != bgeom_nat || !dgamma) {
            for (void* p : { (void*)dxyz, (void*)dcn, (void*)dse, (void*)dgamma, (void*)dL }) if (p) hipFree(p);
            dxyz = dcn = dse = dgamma = dL = nullptr;
            g &= hipMalloc(&dxyz, 3 * static_cast<size_t>(bnat) * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dcn, static_cast<size_t>(bnat) * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dse, static_cast<size_t>(bnsh) * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dgamma, static_cast<size_t>(bnsh) * bnsh * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dL, static_cast<size_t>(bnao) * bnao * sizeof(double)) == hipSuccess;
            bgeom_nat = bnat;
        }
        basis_ready = g;
        return g;
    }

    bool computeIntegrals(const double* xyz) {
        if (!basis_ready || bnao <= 0) return false;
        if (hipMemcpy(dxyz, xyz, 3 * static_cast<size_t>(bnat) * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        { dim3 b(64), gr((bnat + 63) / 64);
          hipLaunchKernelGGL(k_cn, gr, b, 0, stream, bnat, dxyz, dz, bis_gfn2, dcovrad, dcn); }
        { dim3 b(64), gr((bnsh + 63) / 64);
          hipLaunchKernelGGL(k_self_energy, gr, b, 0, stream, bnsh, dselfenergy, dkcn, dsh2at, dcn, dse); }
        { dim3 b(16, 16), gr((bnsh + 15) / 16, (bnsh + 15) / 16);
          hipLaunchKernelGGL(k_overlap_h0, gr, b, 0, stream, bnsh, bnao, bis_gfn2,
                             dsh2at, dang_sh, diao_sh, dnao_sh, dsh_nprim, dsh_prim_off,
                             dprim_alpha, dprim_coeff, dsh_zeta, dshpoly, dse, dz, dvalence, dxyz,
                             darad, dpauling, dS0, dH0); }   // S → dS0, H0 → dH0 (resident)
        { dim3 b(16, 16), gr((bnsh + 15) / 16, (bnsh + 15) / 16);
          hipLaunchKernelGGL(k_gamma, gr, b, 0, stream, bnsh, bis_gfn2, dsh2at, dghard, dxyz, dgamma); }
        // L = chol(S): work on a copy of S0 (resident S must survive).
        if (hipMemcpyAsync(dL, dS0, static_cast<size_t>(bnao) * bnao * sizeof(double),
                           hipMemcpyDeviceToDevice, stream) != hipSuccess) return false;
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        if (rocsolver_dpotrf(handle, rocblas_fill_lower, bnao, dL, bnao, dInfo) != rocblas_status_success) return false;
        if (hipDeviceSynchronize() != hipSuccess) return false;
        rocblas_int info = 0;
        if (hipMemcpy(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return info == 0;   // dS0/dH0 are now the device-built integrals (resident)
    }

    bool dlMat(double* out, const double* src, size_t count) {
        return hipMemcpy(out, src, count * sizeof(double), hipMemcpyDeviceToHost) == hipSuccess;
    }

    void freeBasis() {
        for (void* p : { (void*)dz,(void*)dsh2at,(void*)dang_sh,(void*)diao_sh,(void*)dnao_sh,
                         (void*)dsh_nprim,(void*)dsh_prim_off,(void*)dvalence,(void*)dprim_alpha,
                         (void*)dprim_coeff,(void*)dsh_zeta,(void*)dselfenergy,(void*)dkcn,
                         (void*)dshpoly,(void*)dghard,(void*)dcovrad,(void*)dpauling,(void*)darad,
                         (void*)dxyz,(void*)dcn,(void*)dse,(void*)dgamma,(void*)dL })
            if (p) hipFree(p);
        dz = dsh2at = dang_sh = diao_sh = dnao_sh = dsh_nprim = dsh_prim_off = dvalence = nullptr;
        dprim_alpha = dprim_coeff = dsh_zeta = dselfenergy = dkcn = dshpoly = dghard = nullptr;
        dcovrad = dpauling = darad = dxyz = dcn = dse = dgamma = dL = nullptr;
        basis_ready = false; bgeom_nat = 0;
    }
#endif // HAVE_ROCSOLVER
};

XtbHipContext::XtbHipContext()
    : m_impl(std::make_unique<Impl>())
{
    int count = 0;
    if (hipGetDeviceCount(&count) != hipSuccess || count <= 0)
        return;
    if (hipSetDevice(0) != hipSuccess)
        return;

    hipDeviceProp_t prop;
    if (hipGetDeviceProperties(&prop, 0) != hipSuccess)
        return;
    if (hipStreamCreate(&m_impl->stream) != hipSuccess)
        return;

    m_impl->device = 0;
    m_impl->name   = prop.name;
    m_impl->ok     = true;
}

XtbHipContext::~XtbHipContext()
{
    if (!m_impl) return;
#ifdef HAVE_ROCSOLVER
    m_impl->freeBasis();
    m_impl->freeResident();
    m_impl->freeSolver();
    if (m_impl->handle_ok) rocblas_destroy_handle(m_impl->handle);
#endif
    if (m_impl->stream)
        hipStreamDestroy(m_impl->stream);
}

bool XtbHipContext::ok() const { return m_impl && m_impl->ok; }

std::string XtbHipContext::deviceName() const
{
    return m_impl ? m_impl->name : std::string();
}

int XtbHipContext::deviceId() const { return m_impl ? m_impl->device : -1; }

bool XtbHipContext::deviceAvailable()
{
    int count = 0;
    return hipGetDeviceCount(&count) == hipSuccess && count > 0;
}

bool XtbHipContext::solveGeneralized(const double* F, const double* S, int n,
                                     double* eps_out, double* C_colmajor)
{
#ifndef HAVE_ROCSOLVER
    (void)F; (void)S; (void)n; (void)eps_out; (void)C_colmajor;
    return false;   // built without rocSOLVER → caller uses the CPU eigensolver
#else
    if (!m_impl || !m_impl->ok || n <= 0 || !F || !S || !eps_out || !C_colmajor)
        return false;
    Impl& d = *m_impl;
    if (!d.ensureSolver(n)) return false;
    const size_t nn = static_cast<size_t>(n) * n;

    // F C = S C ε  (itype = AX), eigenvectors of the ORIGINAL problem (CᵀSC = I).
    // dF is overwritten with the eigenvectors, dS with its Cholesky factor.
    if (hipMemcpy(d.dF, F, nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
    if (hipMemcpy(d.dS, S, nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;

    const rocblas_status st = rocsolver_dsygvd(
        d.handle, rocblas_eform_ax, rocblas_evect_original, rocblas_fill_lower,
        n, d.dF, n, d.dS, n, d.dEps, d.dE, d.dInfo);
    if (st != rocblas_status_success) return false;
    if (hipDeviceSynchronize() != hipSuccess) return false;

    rocblas_int info = 0;
    if (hipMemcpy(&info, d.dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
    if (info != 0) return false;   // S not positive-definite or no convergence

    if (hipMemcpy(C_colmajor, d.dF, nn * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
    if (hipMemcpy(eps_out, d.dEps, n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
    return true;
#endif
}

bool XtbHipContext::residentBegin(const double* H0, const double* S, int n)
{
#ifndef HAVE_ROCSOLVER
    (void)H0; (void)S; (void)n; return false;
#else
    return m_impl && m_impl->ok && H0 && S && n > 0 && m_impl->resBegin(H0, S, n);
#endif
}
bool XtbHipContext::residentSolve(const double* v_ao, double* eps_out)
{
#ifndef HAVE_ROCSOLVER
    (void)v_ao; (void)eps_out; return false;
#else
    return m_impl && m_impl->ok && v_ao && eps_out && m_impl->resSolve(v_ao, eps_out);
#endif
}
bool XtbHipContext::residentDensity(const double* occ, int ncol, double* pop_ao, double* band)
{
#ifndef HAVE_ROCSOLVER
    (void)occ; (void)ncol; (void)pop_ao; (void)band; return false;
#else
    return m_impl && m_impl->ok && occ && pop_ao && band && m_impl->resDensity(occ, ncol, pop_ao, band);
#endif
}
bool XtbHipContext::residentFinalize(double* P_colmajor, double* C_colmajor)
{
#ifndef HAVE_ROCSOLVER
    (void)P_colmajor; (void)C_colmajor; return false;
#else
    return m_impl && m_impl->ok && P_colmajor && C_colmajor && m_impl->resFinalize(P_colmajor, C_colmajor);
#endif
}

bool XtbHipContext::beginBasis(const XtbHipBasisData& basis)
{
#ifndef HAVE_ROCSOLVER
    (void)basis; return false;
#else
    return m_impl && m_impl->ok && m_impl->begBasis(basis);
#endif
}
bool XtbHipContext::beginComputed(const double* xyz_bohr)
{
#ifndef HAVE_ROCSOLVER
    (void)xyz_bohr; return false;
#else
    return m_impl && m_impl->ok && xyz_bohr && m_impl->computeIntegrals(xyz_bohr);
#endif
}
bool XtbHipContext::downloadOverlap(double* S_colmajor)
{
#ifndef HAVE_ROCSOLVER
    (void)S_colmajor; return false;
#else
    return m_impl && m_impl->ok && S_colmajor && m_impl->basis_ready
        && m_impl->dlMat(S_colmajor, m_impl->dS0, static_cast<size_t>(m_impl->bnao) * m_impl->bnao);
#endif
}
bool XtbHipContext::downloadH0(double* H0_colmajor)
{
#ifndef HAVE_ROCSOLVER
    (void)H0_colmajor; return false;
#else
    return m_impl && m_impl->ok && H0_colmajor && m_impl->basis_ready
        && m_impl->dlMat(H0_colmajor, m_impl->dH0, static_cast<size_t>(m_impl->bnao) * m_impl->bnao);
#endif
}
bool XtbHipContext::downloadCholesky(double* L_colmajor)
{
#ifndef HAVE_ROCSOLVER
    (void)L_colmajor; return false;
#else
    return m_impl && m_impl->ok && L_colmajor && m_impl->basis_ready
        && m_impl->dlMat(L_colmajor, m_impl->dL, static_cast<size_t>(m_impl->bnao) * m_impl->bnao);
#endif
}
bool XtbHipContext::downloadGamma(double* gamma_colmajor)
{
#ifndef HAVE_ROCSOLVER
    (void)gamma_colmajor; return false;
#else
    return m_impl && m_impl->ok && gamma_colmajor && m_impl->basis_ready
        && m_impl->dlMat(gamma_colmajor, m_impl->dgamma, static_cast<size_t>(m_impl->bnsh) * m_impl->bnsh);
#endif
}

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_ROCM_XTB
