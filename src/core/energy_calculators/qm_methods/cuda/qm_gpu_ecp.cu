/*
 * < CUDA ECP matrix and ECP gradient of the native QM engine (semi-local def2-ECP) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Nov 2026, docs/QM_GPU_ROADMAP.md section 4.6). AI-generated, machine-tested only.
//
// Device counterpart of the heavy parts of qm_ecp.cpp (formulas and literature: qm_ecp.h). Per ECP centre C:
//
//   local      V_ab += sum_k c_k <a| exp(-zeta_k r_C^2) |b>      analytic three-centre overlap (kern::threeCentreOverlap)
//   semi-local V_ab += sum_l sum_k w_k r_k^2 dU_l(r_k) sum_m F_a,lm(r_k) F_b,lm(r_k)
//              F_f,lm(r_k) = sum_terms c sum_ijk B_ijk r_k^(i+j+k) sum_lambda e^(-alpha (r_k - d)^2) [e^-x i_lambda(x)]
//                            T_ijk,lambda,lm,   x = 2 alpha d r_k  (Flores-Moreno half-numerical scheme)
//   gradient   G(a, d) = sum_b P_ab <d_d a| U_C |b> with the derivative functions of the setup (rows nbasis + 3a + d)
//
// The host (qm_ecp.cpp) keeps the quadrature-order loop, the acceptance test and the angular tables T; this file
// receives T and computes at one order:
//   1. ecpRadialKernel   R[slot][lambda][k] = e^(-alpha (r_k-d)^2) e^-x i_lambda(x) per unique (atom, exponent) slot
//   2. ecpProjKernel     F[(lm, k)][f] (column-major, f fastest), one thread per (function, point), 25 accumulators
//   3. ecpOmegaKernel / ecpWeightKernel   FW = F diag(w_k r_k^2 dU_l(r_k)) on the basis rows
//   4. cublasDgemm       V += FW F^T (energy) or X = P FW (gradient), then ecpGradDotKernel G(a,d) += dF . X
// Radial points are processed in chunks whose size depends only on the problem dimensions (no free-memory
// dependence), each thread sums in a fixed order (no atomics): the result is reproducible run to run. The building
// blocks (Bessel series, binomials, Gaussian moments, overlap, radial point, monomial index) are those of
// qm_ecp_kernels.h, shared with the CPU. Limits: l_proj <= 5 (s-g projectors, the lanthanides' h local channel),
// cartesian powers <= 6 (g shells + one derivative), FP64 throughout.

#include "../qm_gpu_ecp_backend.h"
#include "../qm_ecp_kernels.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace {

using qmecp::kern::Prim;

constexpr int kMaxLproj = 5;            // projectors l < 5 -> 25 lm
constexpr int kMaxNlm = kMaxLproj * kMaxLproj;
constexpr int kMaxPow = 6;              // g shell + one derivative = 5; one spare
constexpr int kMaxLam = kMaxPow + kMaxLproj;
constexpr int kBlock = 128;
// Chunk budget in doubles (F, FW, X and R of one chunk of radial points): 2^26 = 512 MB. Fixed, so the chunking
// (and the summation order) depends only on the problem dimensions.
constexpr long long kChunkBudget = 1LL << 26;

void setWhyEcp(char* why, int len, const std::string& s)
{
    if (why && len > 0) {
        std::strncpy(why, s.c_str(), (size_t)len - 1);
        why[len - 1] = '\0';
    }
}

struct DevFns {
    int nbasis, nfn, natoms, nslot;
    const int* fn_atom;
    const int* fn_term;
    const double* fn_xyz;
    const double* term_c;
    const double* term_alpha;
    const int* term_p;
    const int* term_slot;
    const int* slot_atom;
    const double* slot_alpha;
};

struct DevCentre {
    double C[3];
    int nlocal, lproj, nlm, lam_max_all;
    const double* local_coef;
    const double* local_zeta;
    const int* semi_off;
    const double* semi_coef;
    const double* semi_zeta;
    double rcut;
    const double* atom_D;
    const double* atom_d;
    const long long* ang_off;
    const int* ang_nmax;
    const double* ang_T;
};

// --- 1. radial factors per slot and point --------------------------------------------------------------------
__global__ void ecpRadialKernel(DevFns f, DevCentre c, int npts, int k0, int kc, double* R)
{
    const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= (long long)f.nslot * kc) return;
    const int k = (int)(idx % kc), s = (int)(idx / kc);
    const int a = f.slot_atom[s];
    const int L1 = c.lam_max_all + 1;
    double* out = R + (size_t)s * L1 * kc;
    for (int lam = 0; lam < L1; ++lam) out[(size_t)lam * kc + k] = 0.0;
    if (c.ang_off[a] < 0) return;
    double r, w;
    qmecp::kern::radialPoint(k0 + k + 1, npts, r, w);
    if (r > c.rcut) return;
    const double d = c.atom_d[a], alpha = f.slot_alpha[s];
    const double ex = alpha * (r - d) * (r - d);
    if (ex > qmecp::kern::kExpCut) return;
    const bool on_centre = d < 1.0e-12;
    const int lam_max = on_centre ? 0 : c.ang_nmax[a] + c.lproj;
    double bes[kMaxLam + 1];
    qmecp::kern::scaledBessel(2.0 * alpha * d * r, lam_max, bes);
    const double e = exp(-ex);
    for (int lam = 0; lam <= lam_max; ++lam) out[(size_t)lam * kc + k] = e * bes[lam];
}

// Binomial coefficients n over k, n <= 6 (exact integers, equal to kern::binom)
__constant__ double kBinom[kMaxPow + 1][kMaxPow + 1] = {
    { 1, 0, 0, 0, 0, 0, 0 },  { 1, 1, 0, 0, 0, 0, 0 },   { 1, 2, 1, 0, 0, 0, 0 },    { 1, 3, 3, 1, 0, 0, 0 },
    { 1, 4, 6, 4, 1, 0, 0 },  { 1, 5, 10, 10, 5, 1, 0 }, { 1, 6, 15, 20, 15, 6, 1 } };

// --- 2. angular projections F[(lm * kc + k) * nfn + f] -------------------------------------------------------
__global__ void ecpProjKernel(DevFns f, DevCentre c, int npts, int k0, int kc, const double* __restrict__ R,
                              double* __restrict__ F)
{
    const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= (long long)f.nfn * kc) return;
    const int fi = (int)(idx % f.nfn), k = (int)(idx / f.nfn);
    const int a = f.fn_atom[fi];
    const int nlm = c.nlm, L1 = c.lam_max_all + 1;
    double acc[kMaxNlm];
    for (int lm = 0; lm < nlm; ++lm) acc[lm] = 0.0;
    const long long off = c.ang_off[a];
    if (off >= 0) {
        double r, w;
        qmecp::kern::radialPoint(k0 + k + 1, npts, r, w);
        const double d = c.atom_d[a];
        const bool on_centre = d < 1.0e-12;
        const int lam_max = c.ang_nmax[a] + c.lproj;  // table dimension of this atom
        const double* T = c.ang_T + off;
        const double D0 = c.atom_D[3 * a], D1 = c.atom_D[3 * a + 1], D2 = c.atom_D[3 * a + 2];
        double rp[kMaxPow + 1], px[kMaxPow + 1], py[kMaxPow + 1], pz[kMaxPow + 1];
        rp[0] = px[0] = py[0] = pz[0] = 1.0;
        for (int n = 1; n <= kMaxPow; ++n) {
            rp[n] = rp[n - 1] * r;
            // (-D)^e with pow, as the CPU (repeated products differ by a few ulp, amplified by the cancellation of the
            // binomial shift); once per thread, not per term
            px[n] = pow(-D0, (double)n);
            py[n] = pow(-D1, (double)n);
            pz[n] = pow(-D2, (double)n);
        }
        for (int t = f.fn_term[fi]; t < f.fn_term[fi + 1]; ++t) {
            const double alpha = f.term_alpha[t];
            if (d > c.rcut && alpha * (d - c.rcut) * (d - c.rcut) > qmecp::kern::kExpCut) continue;
            const double* Rs = R + (size_t)f.term_slot[t] * L1 * kc + k;
            const int p0 = f.term_p[3 * t], p1 = f.term_p[3 * t + 1], p2 = f.term_p[3 * t + 2];
            const double ct = f.term_c[t];
            for (int i = 0; i <= p0; ++i)
                for (int j = 0; j <= p1; ++j)
                    for (int q = 0; q <= p2; ++q) {
                        const double B = kBinom[p0][i] * px[p0 - i] * kBinom[p1][j] * py[p1 - j] * kBinom[p2][q] * pz[p2 - q];
                        if (B == 0.0) continue;
                        const int n = i + j + q, mi = qmecp::kern::monoIndex(i, j, q);
                        const double base = ct * B * rp[n];
                        const double* Tm = T + (size_t)mi * (lam_max + 1) * nlm;
                        for (int l = 0; l < c.lproj; ++l) {
                            const int lam_top = on_centre ? 0 : min(n + l, lam_max);
                            for (int lam = (n + l) & 1; lam <= lam_top; lam += 2) {
                                const double v = base * Rs[(size_t)lam * kc];
                                if (v == 0.0) continue;
                                const double* Tl = Tm + (size_t)lam * nlm + l * l;
                                for (int m = 0; m < 2 * l + 1; ++m) acc[l * l + m] += v * Tl[m];
                            }
                        }
                    }
        }
    }
    for (int lm = 0; lm < nlm; ++lm) F[((size_t)lm * kc + k) * f.nfn + fi] = acc[lm];
}

// --- 3. FW = F diag(omega) on the basis rows, omega_l(k) = w_k r_k^2 sum_t c_t exp(-zeta_t r_k^2) -------------
__global__ void ecpOmegaKernel(DevCentre c, int npts, int k0, int kc, double* __restrict__ om)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= c.lproj * kc) return;
    const int k = idx % kc, l = idx / kc;
    double r, w;
    qmecp::kern::radialPoint(k0 + k + 1, npts, r, w);
    double u = 0.0;
    for (int t = c.semi_off[l]; t < c.semi_off[l + 1]; ++t) u += c.semi_coef[t] * exp(-c.semi_zeta[t] * r * r);
    om[idx] = w * r * r * u;
}

__global__ void ecpWeightKernel(DevFns f, DevCentre c, int kc, const double* __restrict__ om, const double* __restrict__ F,
                                double* __restrict__ FW)
{
    const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    const long long total = (long long)f.nbasis * c.nlm * kc;
    if (idx >= total) return;
    const int fi = (int)(idx % f.nbasis);
    const long long col = idx / f.nbasis;
    const int k = (int)(col % kc), lm = (int)(col / kc);
    int l = 0;
    while ((l + 1) * (l + 1) <= lm) ++l;
    FW[(size_t)col * f.nbasis + fi] = F[(size_t)col * f.nfn + fi] * om[(size_t)l * kc + k];
}

// --- 4b. G(a, d) += sum_col dF(nbasis + 3a + d, col) X(a, col) ------------------------------------------------
// Two stages with a fixed summation order (no atomics): block (row tile, column split) sums its columns into
// part[split][row], ecpGradSumKernel adds the splits in ascending order. Threads tx run over 32 consecutive rows
// (coalesced F reads), ty over the columns of the split with stride kDotTy.
constexpr int kDotRows = 32, kDotTy = 8, kDotSplits = 64;
__global__ void ecpGradDotKernel(DevFns f, long long ncol, const double* __restrict__ F, const double* __restrict__ X,
                                 double* __restrict__ part)
{
    __shared__ double sh[kDotTy][kDotRows];
    const int tx = threadIdx.x, ty = threadIdx.y;
    const int r = blockIdx.x * kDotRows + tx;  // derivative row 0 .. 3 nbasis - 1
    const long long per = (ncol + kDotSplits - 1) / kDotSplits;
    const long long c0 = (long long)blockIdx.y * per, c1 = min(ncol, c0 + per);
    double s = 0.0;
    if (r < 3 * f.nbasis) {
        const int a = r / 3, row = f.nbasis + r;
        for (long long col = c0 + ty; col < c1; col += kDotTy) s += F[(size_t)col * f.nfn + row] * X[(size_t)col * f.nbasis + a];
    }
    sh[ty][tx] = s;
    __syncthreads();
    if (ty == 0 && r < 3 * f.nbasis) {
        double t = 0.0;
        for (int k = 0; k < kDotTy; ++k) t += sh[k][tx];
        part[(size_t)blockIdx.y * 3 * f.nbasis + r] = t;
    }
}

__global__ void ecpGradSumKernel(int nrow, const double* __restrict__ part, double* __restrict__ G)
{
    const int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= nrow) return;
    double t = 0.0;
    for (int k = 0; k < kDotSplits; ++k) t += part[(size_t)k * nrow + r];
    G[r] += t;
}

// --- local part -------------------------------------------------------------------------------------------------
__device__ double localElementDev(const DevFns& f, const DevCentre& c, int fa, int fb)
{
    const double A[3] = { f.fn_xyz[3 * fa], f.fn_xyz[3 * fa + 1], f.fn_xyz[3 * fa + 2] };
    const double B[3] = { f.fn_xyz[3 * fb], f.fn_xyz[3 * fb + 1], f.fn_xyz[3 * fb + 2] };
    double v = 0.0;
    for (int ta = f.fn_term[fa]; ta < f.fn_term[fa + 1]; ++ta) {
        const Prim pa{ f.term_c[ta], f.term_alpha[ta], { f.term_p[3 * ta], f.term_p[3 * ta + 1], f.term_p[3 * ta + 2] } };
        for (int tb = f.fn_term[fb]; tb < f.fn_term[fb + 1]; ++tb) {
            const Prim pb{ f.term_c[tb], f.term_alpha[tb], { f.term_p[3 * tb], f.term_p[3 * tb + 1], f.term_p[3 * tb + 2] } };
            for (int k = 0; k < c.nlocal; ++k)
                v += pa.c * pb.c * c.local_coef[k] * qmecp::kern::threeCentreOverlap(pa, A, pb, B, c.local_zeta[k], c.C);
        }
    }
    return v;
}

__global__ void ecpLocalKernel(DevFns f, DevCentre c, double* V)
{
    const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    const long long n = f.nbasis;
    if (idx >= n * n) return;
    const int i = (int)(idx / n), j = (int)(idx % n);
    if (j > i) return;
    const double v = localElementDev(f, c, i, j);
    V[(size_t)i * n + j] = v;
    V[(size_t)j * n + i] = v;
}

// E[(3i + d) * nbasis + j] = P_ij <d_d a_i| U_L |b_j>, then G(i, d) = sum_j E in ascending j (ecpRowSumKernel)
__global__ void ecpLocalGradKernel(DevFns f, DevCentre c, const double* __restrict__ P, double* __restrict__ E)
{
    const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    const long long nb = f.nbasis;
    if (idx >= 3 * nb * nb) return;
    const int j = (int)(idx % nb), r = (int)(idx / nb), i = r / 3, d = r % 3;
    const double p = P[(size_t)i * nb + j];
    E[idx] = p != 0.0 ? p * localElementDev(f, c, f.nbasis + 3 * i + d, j) : 0.0;
}

__global__ void ecpRowSumKernel(int nrow, int ncol, const double* __restrict__ E, double* __restrict__ G)
{
    const int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= nrow) return;
    double s = 0.0;
    for (int j = 0; j < ncol; ++j) s += E[(size_t)r * ncol + j];
    G[r] = s;
}

// One cuBLAS handle per device for the life of the process (creating one costs up to ~0.2 s, an ECP object is made
// per geometry). Never destroyed: the plugin stays loaded until exit.
cublasHandle_t cachedCublas(int device)
{
    static std::mutex mtx;
    static std::map<int, cublasHandle_t> handles;
    std::lock_guard<std::mutex> lock(mtx);
    auto it = handles.find(device);
    if (it != handles.end()) return it->second;
    cublasHandle_t h = nullptr;
    if (cublasCreate(&h) != CUBLAS_STATUS_SUCCESS) return nullptr;
    handles[device] = h;
    return h;
}

unsigned blocksFor(long long n) { return (unsigned)std::max<long long>(1, (n + kBlock - 1) / kBlock); }

class QmGpuEcpCuda final : public qmecp::QMGpuECP {
public:
    ~QmGpuEcpCuda() override
    {
        freeCentre();
        freeWork();
        for (void* p : m_owned) cudaFree(p);
        for (auto& e : m_ev)
            if (e) cudaEventDestroy(e);
    }

    bool init(const qmecp::QMGpuEcpSetup* s, int device, std::string& why)
    {
        if (s->nbasis <= 0 || (s->nfn != s->nbasis && s->nfn != 4 * s->nbasis)) {
            why = "ECP setup: nfn must be nbasis or 4 nbasis";
            return false;
        }
        for (int t = 0; t < s->nterm; ++t)
            if (s->term_p[3 * t] + s->term_p[3 * t + 1] + s->term_p[3 * t + 2] > kMaxPow) {
                why = "ECP setup: cartesian power above 6";
                return false;
            }
        m_f.nbasis = s->nbasis;
        m_f.nfn = s->nfn;
        m_f.natoms = s->natoms;
        m_f.nslot = s->nslot;
        int *fa = nullptr, *ft = nullptr, *tp = nullptr, *ts = nullptr, *sa = nullptr;
        double *fx = nullptr, *tc = nullptr, *tal = nullptr, *sal = nullptr;
        const bool ok = up(&fa, s->fn_atom, s->nfn) && up(&ft, s->fn_term, s->nfn + 1) && up(&fx, s->fn_xyz, 3 * (size_t)s->nfn)
            && up(&tc, s->term_c, s->nterm) && up(&tal, s->term_alpha, s->nterm) && up(&tp, s->term_p, 3 * (size_t)s->nterm)
            && up(&ts, s->term_slot, s->nterm) && up(&sa, s->slot_atom, s->nslot) && up(&sal, s->slot_alpha, s->nslot);
        if (!ok) {
            why = "device allocation / upload of the ECP functions failed";
            return false;
        }
        m_f.fn_atom = fa;
        m_f.fn_term = ft;
        m_f.fn_xyz = fx;
        m_f.term_c = tc;
        m_f.term_alpha = tal;
        m_f.term_p = tp;
        m_f.term_slot = ts;
        m_f.slot_atom = sa;
        m_f.slot_alpha = sal;
        int dev = 0;
        cudaGetDevice(&dev);
        m_cublas = cachedCublas(dev);
        if (!m_cublas) {
            why = "cublasCreate failed";
            return false;
        }
        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device >= 0 ? device : 0);
        m_name = std::string("cuda:") + prop.name;
        // CURCUMA_QM_GPU_ECP_TIMING=1: per-kernel event timing in stats() (adds a synchronisation per chunk)
        if (const char* t = std::getenv("CURCUMA_QM_GPU_ECP_TIMING"); t && t[0] == '1') {
            m_timing = true;
            for (auto& e : m_ev) cudaEventCreate(&e);
        }
        return cudaGetLastError() == cudaSuccess;
    }

    std::string name() const override { return m_name; }

    bool setCentre(const qmecp::QMGpuEcpCentre& c) override
    {
        if (c.bytes != (int)sizeof(qmecp::QMGpuEcpCentre)) return false;
        if (c.lproj < 0 || c.lproj > kMaxLproj) return false;
        freeCentre();
        DevCentre d{};
        for (int k = 0; k < 3; ++k) d.C[k] = c.C[k];
        d.nlocal = c.nlocal;
        d.lproj = c.lproj;
        d.nlm = c.lproj * c.lproj;
        d.rcut = c.rcut;
        int lam_all = 0;
        for (int a = 0; a < m_f.natoms; ++a)
            if (c.ang_off[a] >= 0) lam_all = std::max(lam_all, c.ang_nmax[a] + c.lproj);
        if (lam_all > kMaxLam) return false;
        d.lam_max_all = lam_all;
        double *lc = nullptr, *lz = nullptr, *sc = nullptr, *sz = nullptr, *aD = nullptr, *ad = nullptr, *T = nullptr;
        int *so = nullptr, *an = nullptr;
        long long* ao = nullptr;
        const int nsemi = c.lproj > 0 ? c.semi_off[c.lproj] : 0;
        const bool ok = upc(&lc, c.local_coef, c.nlocal) && upc(&lz, c.local_zeta, c.nlocal)
            && upc(&so, c.semi_off, c.lproj + 1) && upc(&sc, c.semi_coef, nsemi) && upc(&sz, c.semi_zeta, nsemi)
            && upc(&aD, c.atom_D, 3 * (size_t)m_f.natoms) && upc(&ad, c.atom_d, m_f.natoms)
            && upc(&ao, c.ang_off, m_f.natoms) && upc(&an, c.ang_nmax, m_f.natoms) && upc(&T, c.ang_T, (size_t)c.ang_size);
        if (!ok) return false;
        d.local_coef = lc;
        d.local_zeta = lz;
        d.semi_off = so;
        d.semi_coef = sc;
        d.semi_zeta = sz;
        d.atom_D = aD;
        d.atom_d = ad;
        d.ang_off = ao;
        d.ang_nmax = an;
        d.ang_T = T;
        m_c = d;
        m_have_centre = true;
        return true;
    }

    bool local(double* V) override
    {
        if (!m_have_centre) return false;
        const auto t0 = std::chrono::steady_clock::now();
        const size_t n2 = (size_t)m_f.nbasis * m_f.nbasis;
        if (!ensure(&m_V, m_V_cap, n2)) return false;
        cudaMemset(m_V, 0, n2 * sizeof(double));
        if (m_c.nlocal > 0) ecpLocalKernel<<<blocksFor((long long)n2), kBlock>>>(m_f, m_c, m_V);
        if (cudaGetLastError() != cudaSuccess) return false;
        const bool ok = cudaMemcpy(V, m_V, n2 * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess;
        m_local_ms += ms(t0);
        return ok;
    }

    bool semilocal(int npts, double* V) override
    {
        if (!m_have_centre || npts <= 0) return false;
        const auto t0 = std::chrono::steady_clock::now();
        const int nb = m_f.nbasis;
        const size_t n2 = (size_t)nb * nb;
        if (!ensure(&m_V, m_V_cap, n2)) return false;
        cudaMemset(m_V, 0, n2 * sizeof(double));
        bool ok = true;
        if (m_c.nlm > 0) {
            const int kc = chunk(npts, false);
            const double one = 1.0;
            for (int k0 = 0; ok && k0 < npts; k0 += kc) {
                const int kn = std::min(kc, npts - k0);
                ok = project(npts, k0, kn, false);
                const long long ncol = (long long)m_c.nlm * kn;
                ok = ok && cublasDgemm(m_cublas, CUBLAS_OP_N, CUBLAS_OP_T, nb, nb, (int)ncol, &one, m_FW, nb, m_F, m_f.nfn,
                                       &one, m_V, nb) == CUBLAS_STATUS_SUCCESS;
            }
        }
        // column-major V on the device = V^T row-major; V is symmetric up to rounding, transpose on the host to keep
        // the GEMM orientation of the CPU (V = PiW Pi^T, row a of PiW)
        std::vector<double> h(n2);
        ok = ok && cudaMemcpy(h.data(), m_V, n2 * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess
            && cudaGetLastError() == cudaSuccess;
        if (ok)
            for (int i = 0; i < nb; ++i)
                for (int j = 0; j < nb; ++j) V[(size_t)i * nb + j] = h[(size_t)j * nb + i];
        m_semi_ms += ms(t0);
        ++m_semi_calls;
        m_points += npts;
        return ok;
    }

    bool gradient(int npts, const double* P, double* G) override
    {
        if (!m_have_centre || m_f.nfn != 4 * m_f.nbasis || npts <= 0) return false;
        const auto t0 = std::chrono::steady_clock::now();
        const int nb = m_f.nbasis;
        const size_t n2 = (size_t)nb * nb;
        if (!ensure(&m_P, m_P_cap, n2) || !ensure(&m_G, m_G_cap, 3 * (size_t)nb)) return false;
        if (cudaMemcpy(m_P, P, n2 * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        cudaMemset(m_G, 0, 3 * (size_t)nb * sizeof(double));
        if (m_c.nlocal > 0) {
            if (!ensure(&m_E, m_E_cap, 3 * n2)) return false;
            ecpLocalGradKernel<<<blocksFor(3LL * (long long)n2), kBlock>>>(m_f, m_c, m_P, m_E);
            ecpRowSumKernel<<<blocksFor(3LL * nb), kBlock>>>(3 * nb, nb, m_E, m_G);
        }
        bool ok = cudaGetLastError() == cudaSuccess;
        if (ok && m_c.nlm > 0) {
            const int kc = chunk(npts, true);
            const double one = 1.0, zero = 0.0;
            for (int k0 = 0; ok && k0 < npts; k0 += kc) {
                const int kn = std::min(kc, npts - k0);
                ok = project(npts, k0, kn, true);
                const long long ncol = (long long)m_c.nlm * kn;
                // X = P FW: P is row-major, i.e. P^T in the column-major view; OP_T restores P
                ok = ok && cublasDgemm(m_cublas, CUBLAS_OP_T, CUBLAS_OP_N, nb, (int)ncol, nb, &one, m_P, nb, m_FW, nb, &zero,
                                       m_X, nb) == CUBLAS_STATUS_SUCCESS;
                if (ok && ensure(&m_part, m_part_cap, (size_t)kDotSplits * 3 * nb)) {
                    const dim3 grid((unsigned)((3 * nb + kDotRows - 1) / kDotRows), kDotSplits), block(kDotRows, kDotTy);
                    ecpGradDotKernel<<<grid, block>>>(m_f, ncol, m_F, m_X, m_part);
                    ecpGradSumKernel<<<blocksFor(3LL * nb), kBlock>>>(3 * nb, m_part, m_G);
                } else
                    ok = false;
                ok = ok && cudaGetLastError() == cudaSuccess;
            }
        }
        ok = ok && cudaMemcpy(G, m_G, 3 * (size_t)nb * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess
            && cudaGetLastError() == cudaSuccess;
        m_grad_ms += ms(t0);
        m_points += npts;
        return ok;
    }

    std::string stats() const override
    {
        char buf[240];
        std::snprintf(buf, sizeof buf, "%d functions, local %.1f ms, semi-local %d orders %.1f ms, gradient %.1f ms, %lld radial points",
                      m_f.nfn, m_local_ms, m_semi_calls, m_semi_ms, m_grad_ms, m_points);
        std::string out = buf;
        if (m_timing) {
            std::snprintf(buf, sizeof buf, "; kernels radial %.1f proj %.1f weight %.1f ms", m_t_radial, m_t_proj, m_t_weight);
            out += buf;
        }
        return out;
    }

private:
    static double ms(std::chrono::steady_clock::time_point t0)
    {
        return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    }

    // Radial points per chunk from the dimensions only (kChunkBudget doubles for R, F, FW and X).
    int chunk(int npts, bool grad) const
    {
        const long long per = (long long)m_f.nslot * (m_c.lam_max_all + 1)
            + (long long)m_c.nlm * (m_f.nfn + (grad ? 2LL : 1LL) * m_f.nbasis);
        return (int)std::max<long long>(1, std::min<long long>(npts, kChunkBudget / std::max<long long>(1, per)));
    }

    // R, F and FW of the points k0 .. k0 + kn - 1 of the npts-point order (and the X buffer for a gradient).
    bool project(int npts, int k0, int kn, bool grad)
    {
        const size_t nR = (size_t)m_f.nslot * (m_c.lam_max_all + 1) * kn;
        const size_t nF = (size_t)m_f.nfn * m_c.nlm * kn, nW = (size_t)m_f.nbasis * m_c.nlm * kn;
        if (!ensure(&m_R, m_R_cap, nR) || !ensure(&m_F, m_F_cap, nF) || !ensure(&m_FW, m_FW_cap, nW)) return false;
        if (grad && !ensure(&m_X, m_X_cap, nW)) return false;
        if (m_timing) cudaEventRecord(m_ev[0]);
        ecpRadialKernel<<<blocksFor((long long)m_f.nslot * kn), kBlock>>>(m_f, m_c, npts, k0, kn, m_R);
        if (m_timing) cudaEventRecord(m_ev[1]);
        ecpProjKernel<<<blocksFor((long long)m_f.nfn * kn), kBlock>>>(m_f, m_c, npts, k0, kn, m_R, m_F);
        if (m_timing) cudaEventRecord(m_ev[2]);
        if (!ensure(&m_om, m_om_cap, (size_t)m_c.lproj * kn)) return false;
        ecpOmegaKernel<<<blocksFor((long long)m_c.lproj * kn), kBlock>>>(m_c, npts, k0, kn, m_om);
        ecpWeightKernel<<<blocksFor((long long)nW), kBlock>>>(m_f, m_c, kn, m_om, m_F, m_FW);
        if (m_timing) {
            cudaEventRecord(m_ev[3]);
            cudaEventSynchronize(m_ev[3]);
            float a = 0, b = 0, c = 0;
            cudaEventElapsedTime(&a, m_ev[0], m_ev[1]);
            cudaEventElapsedTime(&b, m_ev[1], m_ev[2]);
            cudaEventElapsedTime(&c, m_ev[2], m_ev[3]);
            m_t_radial += a;
            m_t_proj += b;
            m_t_weight += c;
        }
        return cudaGetLastError() == cudaSuccess;
    }

    bool ensure(double** p, size_t& cap, size_t count)
    {
        if (count <= cap && *p) return true;
        if (*p) cudaFree(*p);
        *p = nullptr;
        cap = 0;
        if (cudaMalloc(p, std::max<size_t>(count, 1) * sizeof(double)) != cudaSuccess) {
            *p = nullptr;
            cudaGetLastError();
            return false;
        }
        cap = count;
        return true;
    }

    template <typename T>
    bool up(T** d, const T* h, size_t n)
    {
        if (cudaMalloc(d, std::max<size_t>(n, 1) * sizeof(T)) != cudaSuccess) { *d = nullptr; return false; }
        m_owned.push_back(*d);
        return n == 0 || cudaMemcpy(*d, h, n * sizeof(T), cudaMemcpyHostToDevice) == cudaSuccess;
    }
    template <typename T>
    bool upc(T** d, const T* h, size_t n)
    {
        if (cudaMalloc(d, std::max<size_t>(n, 1) * sizeof(T)) != cudaSuccess) { *d = nullptr; return false; }
        m_centre_owned.push_back(*d);
        return n == 0 || cudaMemcpy(*d, h, n * sizeof(T), cudaMemcpyHostToDevice) == cudaSuccess;
    }

    void freeCentre()
    {
        for (void* p : m_centre_owned) cudaFree(p);
        m_centre_owned.clear();
        m_have_centre = false;
    }
    void freeWork()
    {
        for (double* p : { m_R, m_F, m_FW, m_X, m_V, m_P, m_G, m_E, m_part, m_om })
            if (p) cudaFree(p);
        m_R = m_F = m_FW = m_X = m_V = m_P = m_G = m_E = m_part = m_om = nullptr;
    }

    DevFns m_f{};
    DevCentre m_c{};
    bool m_have_centre = false;
    std::vector<void*> m_owned, m_centre_owned;
    double *m_R = nullptr, *m_F = nullptr, *m_FW = nullptr, *m_X = nullptr, *m_V = nullptr, *m_P = nullptr, *m_G = nullptr;
    double *m_E = nullptr, *m_part = nullptr, *m_om = nullptr;
    size_t m_R_cap = 0, m_F_cap = 0, m_FW_cap = 0, m_X_cap = 0, m_V_cap = 0, m_P_cap = 0, m_G_cap = 0, m_E_cap = 0, m_part_cap = 0,
           m_om_cap = 0;
    cublasHandle_t m_cublas = nullptr;
    std::string m_name = "cuda";
    double m_local_ms = 0.0, m_semi_ms = 0.0, m_grad_ms = 0.0;
    int m_semi_calls = 0;
    bool m_timing = false;
    cudaEvent_t m_ev[4] = { nullptr, nullptr, nullptr, nullptr };
    double m_t_radial = 0.0, m_t_proj = 0.0, m_t_weight = 0.0;
    long long m_points = 0;
};

}  // namespace

// ECP device object (QMGpuECP*) or nullptr with a reason in `why`. `setup_bytes` must equal sizeof(QMGpuEcpSetup)
// of this plugin (layout check); the _v1 in the name changes with the vtable. Claude Generated (Nov 2026).
extern "C" void* curcuma_cuda_create_qm_ecp_v1(const qmecp::QMGpuEcpSetup* s, int setup_bytes, int device, char* why,
                                               int whylen)
{
    if (!s || setup_bytes != (int)sizeof(qmecp::QMGpuEcpSetup) || s->bytes != (int)sizeof(qmecp::QMGpuEcpSetup)) {
        setWhyEcp(why, whylen, "QMGpuEcpSetup layout differs between core and plugin");
        return nullptr;
    }
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess) {
        setWhyEcp(why, whylen, "cudaSetDevice failed");
        return nullptr;
    }
    cudaFree(nullptr);  // context creation
    auto b = std::make_unique<QmGpuEcpCuda>();
    std::string err;
    if (!b->init(s, device, err)) {
        cudaGetLastError();
        setWhyEcp(why, whylen, err);
        return nullptr;
    }
    return static_cast<qmecp::QMGpuECP*>(b.release());
}
