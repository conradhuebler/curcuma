/*
 * <Canonical closed-shell CCSD on top of the native RHF (WP13) -- implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Equations, references and cost: see qm_ccsd.h. This file holds, in reading order:
 *   1. ccsdContract / addPermuted   the two tensor helpers (an einsum restricted to "pairwise
 *                                   contraction" = permute + one GEMM, and a permuted addition)
 *   2. ccsdIntegrals                AO -> MO transformation and the integral blocks
 *   3. ccsdEnergy, ccsdResidual     the CCSD energy and amplitude equations
 *   4. ccsdSolve                    Jacobi + DIIS iterations
 *   5. QMEngine::computeCCSD        checks, AO tensor, memory limit, driver
 *
 * Every contraction in ccsdResidual carries the index string of the Python reference
 * scripts/local_ccsd/ccsd_canonical.py (validated against PySCF and ORCA), so the C++ can be
 * checked line by line against it.
 *
 * BLAS: the GEMMs are Eigen products, which call the linked BLAS (single-threaded MKL on the
 * development machine, see the WP10 doc "BLAS (Korrektur)"); the threading is our OpenMP over
 * row or column blocks of each product, one BLAS call per block.
 *
 * Claude Generated (Oct 2026, WP13). AI-generated, machine-tested.
 *
 * This program is free software under GPL-3.0
 */

#include "qm_ccsd.h"
#include "qm_engine.h"
#include "qm_mp2.h"

#include "src/core/curcuma_logger.h"
#include "src/core/gpu_fallback.h"

#include <fmt/format.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using RowMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using RowMap = Eigen::Map<RowMat>;
using ConstRowMap = Eigen::Map<const RowMat>;

double msSince(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

// ---------------------------------------------------------------------------------------------
// 1. Tensor helpers
// ---------------------------------------------------------------------------------------------

/// dst[dl] (+)= alpha * src[sl]: copy (or add) a row-major tensor with index labels `sl` into the
/// index order `dl` (a permutation of sl). When the last destination axis is not the last source
/// axis, the two axes are walked in 32 x 32 tiles, so both reads and writes stay in cache (a plain
/// loop reads the source with a large stride). Claude Generated (Oct 2026, WP13).
void permuteInto(const double* src, const std::vector<int>& sdims, const std::string& sl, const std::string& dl,
                 double* dst, double alpha, bool accumulate, int threads)
{
    const int r = (int)dl.size();
    if (r == 0) {
        dst[0] = (accumulate ? dst[0] : 0.0) + alpha * src[0];
        return;
    }
    std::vector<size_t> sstride(sl.size(), 1);
    for (int k = (int)sl.size() - 2; k >= 0; --k) sstride[k] = sstride[k + 1] * (size_t)sdims[k + 1];
    std::vector<int> ddims(r);
    std::vector<size_t> map(r), dstride(r, 1);  // source stride / destination stride of each destination axis
    int dsrc_last = r - 1;                       // destination axis that is the last (stride-1) source axis
    for (int d = 0; d < r; ++d) {
        const size_t s = sl.find(dl[d]);
        ddims[d] = sdims[s];
        map[d] = sstride[s];
        if ((int)s == (int)sl.size() - 1) dsrc_last = d;
    }
    for (int d = r - 2; d >= 0; --d) dstride[d] = dstride[d + 1] * (size_t)ddims[d + 1];
    size_t total = 1;
    for (int d = 0; d < r; ++d) total *= (size_t)ddims[d];
    const size_t inner = (size_t)ddims[r - 1], istr = map[r - 1];

    if (dsrc_last == r - 1) {
        // last axes coincide: contiguous rows on both sides
        const size_t outer = total / inner;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (total > 200000)
#endif
        for (long long row = 0; row < (long long)outer; ++row) {
            size_t rest = (size_t)row, base = 0;  // decompose the row index over the leading axes
            for (int d = r - 2; d >= 0; --d) {
                base += (rest % (size_t)ddims[d]) * map[d];
                rest /= (size_t)ddims[d];
            }
            double* out = dst + (size_t)row * inner;
            const double* in = src + base;
            if (accumulate)
                for (size_t k = 0; k < inner; ++k) out[k] += alpha * in[k * istr];
            else
                for (size_t k = 0; k < inner; ++k) out[k] = alpha * in[k * istr];
        }
        return;
    }
    // tiled 2D transpose over (destination axis x = dsrc_last, destination axis k = last) for every
    // combination of the remaining axes
    const int x = dsrc_last;
    const size_t nx = (size_t)ddims[x], dxs = dstride[x];
    const size_t outer = total / (nx * inner);
    constexpr size_t TB = 32;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (total > 200000)
#endif
    for (long long row = 0; row < (long long)outer; ++row) {
        size_t rest = (size_t)row, sbase = 0, dbase = 0;
        for (int d = r - 2; d >= 0; --d) {
            if (d == x) continue;
            const size_t id = rest % (size_t)ddims[d];
            rest /= (size_t)ddims[d];
            sbase += id * map[d];
            dbase += id * dstride[d];
        }
        for (size_t x0 = 0; x0 < nx; x0 += TB)
            for (size_t k0 = 0; k0 < inner; k0 += TB) {
                const size_t x1 = std::min(nx, x0 + TB), k1 = std::min(inner, k0 + TB);
                for (size_t k = k0; k < k1; ++k)
                    for (size_t ix = x0; ix < x1; ++ix) {  // source contiguous in ix
                        double& o = dst[dbase + ix * dxs + k];
                        const double val = alpha * src[sbase + ix + k * istr];
                        o = accumulate ? o + val : val;
                    }
            }
    }
}

/// C (M x N) += alpha * op(A) op(B), op(A) M x K, op(B) K x N, row-major contiguous storage:
/// A is stored M x K (transA = false) or K x M (transA = true), B is K x N or N x K. Transposed
/// operands go to BLAS as transpose flags, no copy. OpenMP over row (or column) blocks of C, one
/// block per thread and one BLAS GEMM per block, so each element is summed over k in one call
/// (every block re-packs the whole other operand inside BLAS; 4 blocks per thread were slower).
/// Claude Generated (Oct 2026, WP13).
void gemmAccumulate(const double* A, bool transA, const double* B, bool transB, double* C, size_t M, size_t N,
                    size_t K, double alpha, int threads)
{
    if (M == 0 || N == 0 || K == 0) return;
    using Idx = Eigen::Index;
    RowMap Cm(C, (Idx)M, (Idx)N);
    const ConstRowMap As(A, transA ? (Idx)K : (Idx)M, transA ? (Idx)M : (Idx)K);
    const ConstRowMap Bs(B, transB ? (Idx)N : (Idx)K, transB ? (Idx)K : (Idx)N);
    // rows s..s+l of op(A), columns s..s+l of op(B), as blocks of the stored matrices
    auto prod = [&](Idx rs, Idx rl, Idx cs, Idx cl) {
        auto Cb = Cm.block(rs, cs, rl, cl);
        if (!transA && !transB)
            Cb.noalias() += alpha * (As.middleRows(rs, rl) * Bs.middleCols(cs, cl));
        else if (transA && !transB)
            Cb.noalias() += alpha * (As.middleCols(rs, rl).transpose() * Bs.middleCols(cs, cl));
        else if (!transA && transB)
            Cb.noalias() += alpha * (As.middleRows(rs, rl) * Bs.middleRows(cs, cl).transpose());
        else
            Cb.noalias() += alpha * (As.middleCols(rs, rl).transpose() * Bs.middleRows(cs, cl).transpose());
    };
    const double flops = 2.0 * (double)M * (double)N * (double)K;
    if (threads <= 1 || flops < 4.0e6) {
        prod(0, (Idx)M, 0, (Idx)N);
        return;
    }
    const bool by_rows = M >= N;
    const Idx len = by_rows ? (Idx)M : (Idx)N;
    const Idx chunk = std::max<Idx>(8, (len + threads - 1) / threads);
    const Idx nchunk = (len + chunk - 1) / chunk;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static, 1)
#endif
    for (Idx c = 0; c < nchunk; ++c) {
        const Idx st = c * chunk, l = std::min(chunk, len - st);
        if (by_rows)
            prod(st, l, 0, (Idx)N);
        else
            prod(0, (Idx)M, st, l);
    }
}
}  // namespace

namespace qmint {

CCTensor::~CCTensor() = default;
CCSDIntegrals::~CCSDIntegrals() = default;

/// @brief out += alpha * einsum(spec, A, B) for a pairwise contraction, e.g. "kcld,ilcd->ki".
/// Every index of A and B appears either in the output (a free index of exactly one operand) or in
/// both operands (summed); traces and batch indices are rejected. Done as: A -> [free_A | summed],
/// B -> [summed | free_B] (a transposed operand goes to BLAS as a flag, a permuted copy is made
/// only for other orders), one GEMM, and the result written in the output order. An empty output ("ia,ia->") is a full contraction.
/// Claude Generated (Oct 2026, WP13).
void ccsdContract(const char* spec_c, double alpha, const CCTensor& A, const CCTensor& B, CCTensor& out, int threads)
{
    const std::string spec(spec_c);
    const size_t comma = spec.find(','), arrow = spec.find("->");
    if (comma == std::string::npos || arrow == std::string::npos)
        throw std::runtime_error("ccsdContract: bad spec " + spec);
    const std::string la = spec.substr(0, comma), lb = spec.substr(comma + 1, arrow - comma - 1),
                      lc = spec.substr(arrow + 2);
    if (la.size() != A.dims.size() || lb.size() != B.dims.size() || lc.size() != out.dims.size())
        throw std::runtime_error("ccsdContract: rank mismatch in " + spec);
    int dim[128] = { 0 };
    auto setDim = [&](char ch, int d) {
        if (dim[(int)ch] && dim[(int)ch] != d) throw std::runtime_error("ccsdContract: dimension mismatch in " + spec);
        dim[(int)ch] = d;
    };
    for (size_t k = 0; k < la.size(); ++k) setDim(la[k], A.dims[k]);
    for (size_t k = 0; k < lb.size(); ++k) setDim(lb[k], B.dims[k]);
    for (size_t k = 0; k < lc.size(); ++k) setDim(lc[k], out.dims[k]);

    std::string fa, fb, sum;  // free indices of A and B (in output order), summed indices (in A order)
    for (char ch : lc) {
        const bool ina = la.find(ch) != std::string::npos, inb = lb.find(ch) != std::string::npos;
        if (ina == inb) throw std::runtime_error("ccsdContract: output index must come from exactly one operand: " + spec);
        (ina ? fa : fb) += ch;
    }
    for (char ch : la) {
        const bool inb = lb.find(ch) != std::string::npos, inc = lc.find(ch) != std::string::npos;
        if (!inb && !inc) throw std::runtime_error("ccsdContract: trace index not supported: " + spec);
        if (inb && !inc) sum += ch;
    }
    for (char ch : lb)
        if (la.find(ch) == std::string::npos && lc.find(ch) == std::string::npos)
            throw std::runtime_error("ccsdContract: trace index not supported: " + spec);

    size_t M = 1, N = 1, K = 1;
    for (char ch : fa) M *= (size_t)dim[(int)ch];
    for (char ch : fb) N *= (size_t)dim[(int)ch];
    for (char ch : sum) K *= (size_t)dim[(int)ch];

    // Layout choice (no copy where possible): A is usable as [fa|sum] or, transposed, as [sum|fa];
    // B as [sum|fb] or [fb|sum]. The order of the summed indices is taken from A, or from B when A
    // would need a copy anyway. The output is written directly as [fa|fb], or as [fb|fa] by forming
    // the transposed product op(B)^T op(A)^T; only other output orders go through a temporary.
    auto inOrder = [](const std::string& labels, const std::string& set) {
        std::string r;
        for (char ch : labels)
            if (set.find(ch) != std::string::npos) r += ch;
        return r;
    };
    const std::string sumA = inOrder(la, sum), sumB = inOrder(lb, sum);
    const bool aDirect = (la == fa + sumA || la == sumA + fa);
    const bool bDirectWithB = (lb == sumB + fb || lb == fb + sumB);
    const std::string sk = (!aDirect && bDirectWithB) ? sumB : sumA;

    std::vector<double> bufA, bufB, bufC;
    const double* pa = A.data.data();
    const double* pb = B.data.data();
    bool tA = false, tB = false;
    if (la == fa + sk) {
    } else if (la == sk + fa) {
        tA = true;
    } else {
        bufA.resize(A.size());
        permuteInto(A.data.data(), A.dims, la, fa + sk, bufA.data(), 1.0, false, threads);
        pa = bufA.data();
    }
    if (lb == sk + fb) {
    } else if (lb == fb + sk) {
        tB = true;
    } else {
        bufB.resize(B.size());
        permuteInto(B.data.data(), B.dims, lb, sk + fb, bufB.data(), 1.0, false, threads);
        pb = bufB.data();
    }
    if (lc == fa + fb) {
        gemmAccumulate(pa, tA, pb, tB, out.data.data(), M, N, K, alpha, threads);
    } else if (lc == fb + fa) {
        // out^T = op(B)^T op(A)^T  (N x M)
        gemmAccumulate(pb, !tB, pa, !tA, out.data.data(), N, M, K, alpha, threads);
    } else {
        bufC.assign(M * N, 0.0);
        gemmAccumulate(pa, tA, pb, tB, bufC.data(), M, N, K, 1.0, threads);
        std::vector<int> cdims;
        for (char ch : fa + fb) cdims.push_back(dim[(int)ch]);
        permuteInto(bufC.data(), cdims, fa + fb, lc, out.data.data(), alpha, true, threads);
    }
}

/// @brief out += alpha * A with A's indices relabelled, spec "ijab->jiba" (the Python
/// `A.transpose(1, 0, 3, 2)`). Claude Generated (Oct 2026, WP13).
void addPermuted(const char* spec_c, double alpha, const CCTensor& A, CCTensor& out, int threads)
{
    const std::string spec(spec_c);
    const size_t arrow = spec.find("->");
    permuteInto(A.data.data(), A.dims, spec.substr(0, arrow), spec.substr(arrow + 2), out.data.data(), alpha, true,
                threads);
}

/// Copy of A permuted into new index order (`A.transpose(...)` as a new tensor). Claude Generated.
CCTensor permuted(const char* spec_c, const CCTensor& A, int threads)
{
    const std::string spec(spec_c);
    const size_t arrow = spec.find("->");
    const std::string sl = spec.substr(0, arrow), dl = spec.substr(arrow + 2);
    std::vector<int> d;
    for (char ch : dl) d.push_back(A.dims[sl.find(ch)]);
    CCTensor r(d);
    permuteInto(A.data.data(), A.dims, sl, dl, r.data.data(), 1.0, false, threads);
    return r;
}

// ---------------------------------------------------------------------------------------------
// 2. MO integrals
// ---------------------------------------------------------------------------------------------

double ccsdMemoryMB(int n, int n_frozen, int n_occ, int diis)
{
    const double m = n - n_frozen, o = n_occ - n_frozen, v = n - n_occ, N = n;
    const double v4 = v * v * v * v, ov3 = o * v * v * v, o2v2 = o * o * v * v;
    // transformation: buffers n^3 m and n^2 m^2, then the blocks (v^4 + ov^3 + small ones)
    const double transform = N * N * N * m + N * N * m * m + v4 + ov3 + 4.0 * o2v2;
    // iterations: v^4 + ovvv and up to three permuted copies, DIIS (2 x diis vectors), ~30 o^2v^2 work tensors
    const double iterate = v4 + 4.0 * ov3 + (2.0 * std::max(diis, 1) + 30.0) * (o2v2 + o * v);
    return std::max(transform, iterate) * 8.0 / 1.0e6;
}

namespace {
/// One quarter transformation: out (R x m) = in^T (R x d0) * C (d0 x m), in row-major (d0 x R).
/// The transformed index moves from the front to the back: In[p][q r s] -> Out[q r s][p'].
void quarterTransform(const double* in, size_t d0, size_t R, const RowMat& C, double* out, int threads)
{
    const Eigen::Index m = C.cols();
    const ConstRowMap In(in, (Eigen::Index)d0, (Eigen::Index)R);
    RowMap Out(out, (Eigen::Index)R, m);
    const Eigen::Index rows = (Eigen::Index)R;
    const Eigen::Index chunk = std::max<Eigen::Index>(256, (rows + 8 * threads - 1) / (8 * threads));
    const Eigen::Index nchunk = (rows + chunk - 1) / chunk;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
    for (Eigen::Index c = 0; c < nchunk; ++c) {
        const Eigen::Index r0 = c * chunk, len = std::min(chunk, rows - r0);
        Out.middleRows(r0, len).noalias() = In.middleCols(r0, len).transpose() * C;
    }
}
}  // namespace

// Claude Generated (Oct 2026, WP13): full transformation to the active MO space, then the blocks.
void ccsdIntegrals(const double* eri, int n, const Matrix& C, const Vector& eps, int n_frozen, int n_occ,
                   int threads, CCSDIntegrals& I)
{
    threads = std::max(1, threads);
    const int m = n - n_frozen, o = n_occ - n_frozen, v = n - n_occ;
    I.o = o;
    I.v = v;
    const RowMat Cact = C.block(0, n_frozen, n, m);
    const size_t N = n, M = m;

    // (pq|rs) -> four quarter steps, each moves the transformed index to the back:
    //   [p q r s] -> [q r s p'] -> [r s p' q'] -> [s p' q' r'] -> [p' q' r' s']
    std::vector<double> buf1(N * N * N * M), buf2(N * N * M * M);
    quarterTransform(eri, N, N * N * N, Cact, buf1.data(), threads);
    quarterTransform(buf1.data(), N, N * N * M, Cact, buf2.data(), threads);
    quarterTransform(buf2.data(), N, N * M * M, Cact, buf1.data(), threads);
    quarterTransform(buf1.data(), N, M * M * M, Cact, buf2.data(), threads);
    std::vector<double>().swap(buf1);
    const double* g = buf2.data();  // (pq|rs) in the active MO basis, occupied 0..o-1, virtual o..m-1
    auto G = [&](int p, int q, int r, int s) { return g[(((size_t)p * M + q) * M + r) * M + s]; };

    I.oooo = CCTensor({ o, o, o, o });
    I.ovoo = CCTensor({ o, v, o, o });
    I.ovov = CCTensor({ o, v, o, v });
    I.oovv = CCTensor({ o, o, v, v });
    I.ovvo = CCTensor({ o, v, v, o });
    I.ovvv = CCTensor({ o, v, v, v });
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
    for (int i = 0; i < o; ++i) {
        for (int j = 0; j < o; ++j)
            for (int k = 0; k < o; ++k)
                for (int l = 0; l < o; ++l) I.oooo.data[(((size_t)i * o + j) * o + k) * o + l] = G(i, j, k, l);
        for (int a = 0; a < v; ++a) {
            for (int j = 0; j < o; ++j)
                for (int k = 0; k < o; ++k) I.ovoo.data[(((size_t)i * v + a) * o + j) * o + k] = G(i, o + a, j, k);
            for (int j = 0; j < o; ++j)
                for (int b = 0; b < v; ++b) I.ovov.data[(((size_t)i * v + a) * o + j) * v + b] = G(i, o + a, j, o + b);
            for (int b = 0; b < v; ++b) {
                for (int j = 0; j < o; ++j) I.ovvo.data[(((size_t)i * v + a) * v + b) * o + j] = G(i, o + a, o + b, j);
                for (int c = 0; c < v; ++c) I.ovvv.data[(((size_t)i * v + a) * v + b) * v + c] = G(i, o + a, o + b, o + c);
            }
        }
        for (int j = 0; j < o; ++j)
            for (int a = 0; a < v; ++a)
                for (int b = 0; b < v; ++b) I.oovv.data[(((size_t)i * o + j) * v + a) * v + b] = G(i, j, o + a, o + b);
    }
    // vvvv as the ladder matrix V[(a,b)][(c,d)] = (ac|bd)
    const size_t V = v;
    I.vvvv.assign(V * V * V * V, 0.0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
    for (int a = 0; a < v; ++a)
        for (int b = 0; b < v; ++b)
            for (int c = 0; c < v; ++c)
                for (int d = 0; d < v; ++d)
                    I.vvvv[(((size_t)a * V + b) * V + c) * V + d] = G(o + a, o + c, o + b, o + d);

    // Canonical Fock blocks: f_pq = e_p delta_pq in the RHF MO basis
    I.foo = CCTensor({ o, o });
    I.fvv = CCTensor({ v, v });
    I.fov = CCTensor({ o, v });
    for (int i = 0; i < o; ++i) I.foo.data[(size_t)i * o + i] = eps(n_frozen + i);
    for (int a = 0; a < v; ++a) I.fvv.data[(size_t)a * v + a] = eps(n_occ + a);
}

// ---------------------------------------------------------------------------------------------
// 3. Energy and residual
// ---------------------------------------------------------------------------------------------

namespace {
/// tau_ijab = t2_ijab + t1_ia t1_jb
CCTensor makeTau(const CCTensor& t1, const CCTensor& t2, int threads)
{
    CCTensor tau = t2;
    ccsdContract("ia,jb->ijab", 1.0, t1, t1, tau, threads);
    return tau;
}

/// e_ia = f_ii - f_aa (canonical: orbital energy differences, negative)
double eia(const CCSDIntegrals& I, int i, int a)
{
    return I.foo.data[(size_t)i * I.o + i] - I.fvv.data[(size_t)a * I.v + a];
}

/// Removes the diagonal of f (canonical part, moved to the denominator) from X.
void subtractDiagonal(CCTensor& X, const CCTensor& f)
{
    const int d = X.dims[0];
    for (int p = 0; p < d; ++p) X.data[(size_t)p * d + p] -= f.data[(size_t)p * d + p];
}
}  // namespace

double ccsdEnergy(const CCTensor& t1, const CCTensor& t2, const CCSDIntegrals& I, int threads)
{
    // E = 2 sum f_ia t_ia + 2 sum tau_ijab (ia|jb) - sum tau_ijab (ib|ja)
    const CCTensor tau = makeTau(t1, t2, threads);
    CCTensor e(std::vector<int>{});  // rank 0: one number
    ccsdContract("ia,ia->", 2.0, I.fov, t1, e, threads);
    ccsdContract("ijab,iajb->", 2.0, tau, I.ovov, e, threads);
    ccsdContract("ijab,ibja->", -1.0, tau, I.ovov, e, threads);
    return e.data[0];
}

// Claude Generated (Oct 2026, WP13): port of residual() in scripts/local_ccsd/ccsd_canonical.py,
// which transcribes PySCF rintermediates.py / rccsd.update_amps (Hirata et al. 2004, Eqs. 35-45).
void ccsdResidual(const CCTensor& t1, const CCTensor& t2, const CCSDIntegrals& I, int threads, CCTensor& R1,
                  CCTensor& R2)
{
    const int o = I.o, v = I.v, T = threads;
    const CCTensor &ovov = I.ovov, &oovv = I.oovv, &ovvo = I.ovvo, &ovoo = I.ovoo, &ovvv = I.ovvv;
    const CCTensor &foo = I.foo, &fvv = I.fvv, &fov = I.fov;

    // ---- kappa intermediates F (Eqs. 37-39) ----
    // W1_kc = 2 sum_ld (kc|ld) t_ld - sum_ld (kd|lc) t_ld ;  Fov = f_ov + W1
    CCTensor W1({ o, v });
    ccsdContract("kcld,ld->kc", 2.0, ovov, t1, W1, T);
    ccsdContract("kdlc,ld->kc", -1.0, ovov, t1, W1, T);
    CCTensor Foo = foo, Fvv = fvv, Fov = fov;
    ccsdContract("kcld,ilcd->ki", 2.0, ovov, t2, Foo, T);
    ccsdContract("kdlc,ilcd->ki", -1.0, ovov, t2, Foo, T);
    ccsdContract("kc,ic->ki", 1.0, W1, t1, Foo, T);   // = 2 (kcld,ic,ld) - (kdlc,ic,ld)
    ccsdContract("kcld,klad->ac", -2.0, ovov, t2, Fvv, T);
    ccsdContract("kdlc,klad->ac", 1.0, ovov, t2, Fvv, T);
    ccsdContract("kc,ka->ac", -1.0, W1, t1, Fvv, T);  // = -2 (kcld,ka,ld) + (kdlc,ka,ld)
    for (size_t x = 0; x < Fov.size(); ++x) Fov.data[x] += W1.data[x];

    // ---- lambda intermediates L (Eqs. 40-41), from F before its diagonal is removed ----
    // S_ki = 2 sum_lc (lc|ki) t_lc - (kc|li) t_lc ;  Rv_ac = 2 sum_kd (kd|ac) t_kd - (kc|ad) t_kd
    CCTensor S({ o, o }), Rv({ v, v });
    ccsdContract("lcki,lc->ki", 2.0, ovoo, t1, S, T);
    ccsdContract("kcli,lc->ki", -1.0, ovoo, t1, S, T);
    ccsdContract("kdac,kd->ac", 2.0, ovvv, t1, Rv, T);
    ccsdContract("kcad,kd->ac", -1.0, ovvv, t1, Rv, T);
    CCTensor Loo = Foo, Lvv = Fvv;
    ccsdContract("kc,ic->ki", 1.0, fov, t1, Loo, T);
    for (size_t x = 0; x < Loo.size(); ++x) Loo.data[x] += S.data[x];
    ccsdContract("kc,ka->ac", -1.0, fov, t1, Lvv, T);
    for (size_t x = 0; x < Lvv.size(); ++x) Lvv.data[x] += Rv.data[x];
    subtractDiagonal(Foo, foo);
    subtractDiagonal(Fvv, fvv);
    subtractDiagonal(Loo, foo);
    subtractDiagonal(Lvv, fvv);

    // ---- T1 numerator N1 ----
    CCTensor N1 = fov;
    {
        CCTensor Q({ v, v });  // -2 (kc,ka,ic->ia) = -2 sum_c [sum_k f_kc t_ka] t_ic
        ccsdContract("kc,ka->ca", 1.0, fov, t1, Q, T);
        ccsdContract("ca,ic->ia", -2.0, Q, t1, N1, T);
    }
    ccsdContract("ac,ic->ia", 1.0, Fvv, t1, N1, T);
    ccsdContract("ki,ka->ia", -1.0, Foo, t1, N1, T);
    ccsdContract("kc,kica->ia", 2.0, Fov, t2, N1, T);
    ccsdContract("kc,ikca->ia", -1.0, Fov, t2, N1, T);
    {
        CCTensor P({ o, o });  // (kc,ic,ka->ia) with Fov
        ccsdContract("kc,ic->ki", 1.0, Fov, t1, P, T);
        ccsdContract("ki,ka->ia", 1.0, P, t1, N1, T);
    }
    ccsdContract("kcai,kc->ia", 2.0, ovvo, t1, N1, T);
    ccsdContract("kiac,kc->ia", -1.0, oovv, t1, N1, T);
    ccsdContract("kdac,ikcd->ia", 2.0, ovvv, t2, N1, T);
    ccsdContract("kcad,ikcd->ia", -1.0, ovvv, t2, N1, T);
    ccsdContract("ac,ic->ia", 1.0, Rv, t1, N1, T);   // = 2 (kdac,kd,ic) - (kcad,kd,ic)
    ccsdContract("lcki,klac->ia", -2.0, ovoo, t2, N1, T);
    ccsdContract("kcli,klac->ia", 1.0, ovoo, t2, N1, T);
    ccsdContract("ki,ka->ia", -1.0, S, t1, N1, T);   // = -2 (lcki,lc,ka) + (kcli,lc,ka)

    // ---- T2 numerator N2; terms of the form X_ijab + X_jiba are collected in Sym ----
    CCTensor N2 = permuted("iajb->ijab", ovov, T);  // (ia|jb)
    CCTensor Sym({ o, o, v, v });
    {
        // tmp2_abic = (ia|cb) - sum_k (ki|bc) t_ka   (Python: ovvv.transpose(1,3,0,2) + einsum("kibc,ka->abic", oovv, -t1))
        CCTensor tmp2 = permuted("iacb->abic", ovvv, T);
        ccsdContract("kibc,ka->abic", -1.0, oovv, t1, tmp2, T);
        ccsdContract("abic,jc->ijab", 1.0, tmp2, t1, Sym, T);
    }
    {
        // tmp2_akij = (ia|jk) + sum_c (kc|ai) t_jc   (Python: ovoo.transpose(1,3,0,2) + einsum("kcai,jc->akij", ovvo, t1))
        CCTensor tmp2 = permuted("iajk->akij", ovoo, T);
        ccsdContract("kcai,jc->akij", 1.0, ovvo, t1, tmp2, T);
        ccsdContract("akij,kb->ijab", -1.0, tmp2, t1, Sym, T);
    }

    // ---- chi intermediates W (Eqs. 42-45) ----
    CCTensor Woooo = permuted("kilj->klij", I.oooo, T);  // (ki|lj)
    ccsdContract("lcki,jc->klij", 1.0, ovoo, t1, Woooo, T);
    ccsdContract("kclj,ic->klij", 1.0, ovoo, t1, Woooo, T);
    ccsdContract("kcld,ijcd->klij", 1.0, ovov, t2, Woooo, T);
    {
        CCTensor X({ o, o, o, v });  // (kcld,ic,jd->klij) in two steps
        ccsdContract("kcld,ic->kild", 1.0, ovov, t1, X, T);
        ccsdContract("kild,jd->klij", 1.0, X, t1, Woooo, T);
    }
    CCTensor Wvoov = permuted("kcai->akic", ovvo, T);  // (kc|ai)
    ccsdContract("kcad,id->akic", 1.0, ovvv, t1, Wvoov, T);
    ccsdContract("kcli,la->akic", -1.0, ovoo, t1, Wvoov, T);
    ccsdContract("ldkc,ilda->akic", -0.5, ovov, t2, Wvoov, T);
    ccsdContract("lckd,ilad->akic", -0.5, ovov, t2, Wvoov, T);
    ccsdContract("ldkc,ilad->akic", 1.0, ovov, t2, Wvoov, T);
    {
        CCTensor X({ o, o, v, o });  // -(ldkc,id,la->akic)
        ccsdContract("ldkc,id->lkci", 1.0, ovov, t1, X, T);
        ccsdContract("lkci,la->akic", -1.0, X, t1, Wvoov, T);
    }
    CCTensor Wvovo = permuted("kiac->akci", oovv, T);  // (ki|ac)
    ccsdContract("kdac,id->akci", 1.0, ovvv, t1, Wvovo, T);
    ccsdContract("lcki,la->akci", -1.0, ovoo, t1, Wvovo, T);
    ccsdContract("lckd,ilda->akci", -0.5, ovov, t2, Wvovo, T);
    {
        CCTensor X({ o, v, o, o });  // -(lckd,id,la->akci)
        ccsdContract("lckd,id->lcki", 1.0, ovov, t1, X, T);
        ccsdContract("lcki,la->akci", -1.0, X, t1, Wvovo, T);
    }

    // ---- remaining T2 terms ----
    const CCTensor tau = makeTau(t1, t2, T);
    ccsdContract("klij,klab->ijab", 1.0, Woooo, tau, N2, T);
    // particle-particle ladder sum_cd (ac|bd) tau_ijcd: one (o^2 x v^2)(v^2 x v^2) GEMM, the O(o^2 v^4) step
    gemmAccumulate(tau.data.data(), false, I.vvvv.data(), false, N2.data.data(), (size_t)o * o, (size_t)v * v, (size_t)v * v, 1.0, T);
    {
        // t1 dressing of the ladder: -sum_k t_kb sum_cd (kd|ac) tau_ijcd  (o^3 v^3, no second v^4 array)
        CCTensor Y({ o, o, o, v });
        ccsdContract("kdac,ijcd->ijka", 1.0, ovvv, tau, Y, T);
        ccsdContract("kb,ijka->ijab", -1.0, t1, Y, Sym, T);
    }
    ccsdContract("ac,ijcb->ijab", 1.0, Lvv, t2, Sym, T);
    ccsdContract("ki,kjab->ijab", -1.0, Loo, t2, Sym, T);
    ccsdContract("akic,kjcb->ijab", 2.0, Wvoov, t2, Sym, T);
    ccsdContract("akci,kjcb->ijab", -1.0, Wvovo, t2, Sym, T);
    ccsdContract("akic,kjbc->ijab", -1.0, Wvoov, t2, Sym, T);
    ccsdContract("bkci,kjac->ijab", -1.0, Wvovo, t2, Sym, T);
    for (size_t x = 0; x < N2.size(); ++x) N2.data[x] += Sym.data[x];
    addPermuted("ijab->jiba", 1.0, Sym, N2, T);

    // ---- residuals R = N - e t ----
    R1 = N1;
    for (int i = 0; i < o; ++i)
        for (int a = 0; a < v; ++a) R1.data[(size_t)i * v + a] -= eia(I, i, a) * t1.data[(size_t)i * v + a];
    R2 = N2;
#ifdef _OPENMP
#pragma omp parallel for num_threads(T) schedule(static)
#endif
    for (int i = 0; i < o; ++i)
        for (int j = 0; j < o; ++j)
            for (int a = 0; a < v; ++a)
                for (int b = 0; b < v; ++b) {
                    const size_t x = (((size_t)i * o + j) * v + a) * v + b;
                    R2.data[x] -= (eia(I, i, a) + eia(I, j, b)) * t2.data[x];
                }
}

// ---------------------------------------------------------------------------------------------
// 4. Solver: Jacobi step + Pulay DIIS
// ---------------------------------------------------------------------------------------------

namespace {
/// Pulay DIIS on flat amplitude vectors, error vector = Jacobi step (scripts/local_mp2/lmp2_common.py DIIS).
struct AmplitudeDIIS {
    int nmax;
    std::vector<std::vector<double>> x, e;
    explicit AmplitudeDIIS(int n) : nmax(n) {}
    std::vector<double> update(const std::vector<double>& xn, const std::vector<double>& en)
    {
        if (nmax <= 1) return xn;
        x.push_back(xn);
        e.push_back(en);
        if ((int)x.size() > nmax) {
            x.erase(x.begin());
            e.erase(e.begin());
        }
        const int n = (int)x.size();
        if (n < 2) return xn;
        Eigen::MatrixXd B = Eigen::MatrixXd::Constant(n + 1, n + 1, -1.0);
        B(n, n) = 0.0;
        for (int a = 0; a < n; ++a)
            for (int b = a; b < n; ++b) {
                double s = 0.0;
                for (size_t k = 0; k < e[a].size(); ++k) s += e[a][k] * e[b][k];
                B(a, b) = B(b, a) = s;
            }
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n + 1);
        rhs(n) = -1.0;
        const Eigen::VectorXd c = B.fullPivLu().solve(rhs);
        if (!c.allFinite()) return xn;
        std::vector<double> r(xn.size(), 0.0);
        for (int a = 0; a < n; ++a)
            for (size_t k = 0; k < r.size(); ++k) r[k] += c(a) * x[a][k];
        return r;
    }
};
}  // namespace

// Claude Generated (Oct 2026, WP13): port of solve_ccsd() in scripts/local_ccsd/ccsd_canonical.py.
void ccsdSolve(const CCSDIntegrals& I, const CCSDOptions& opt, int threads, CCSDResult& res)
{
    const auto t0 = std::chrono::steady_clock::now();
    const int o = I.o, v = I.v;
    const size_t n1 = (size_t)o * v, n2 = n1 * n1;
    // start: t1 = 0, t2 = (ia|jb) / e_ijab (MP2 amplitudes) -> iteration-0 energy = MP2
    CCTensor t1({ o, v }), t2({ o, o, v, v });
    for (int i = 0; i < o; ++i)
        for (int j = 0; j < o; ++j)
            for (int a = 0; a < v; ++a)
                for (int b = 0; b < v; ++b)
                    t2.data[(((size_t)i * o + j) * v + a) * v + b] =
                        I.ovov.data[(((size_t)i * v + a) * o + j) * v + b] / (eia(I, i, a) + eia(I, j, b));
    double e = ccsdEnergy(t1, t2, I, threads);
    res.e_mp2 = e;
    double e_old = e;
    AmplitudeDIIS diis(opt.diis);
    CCTensor R1, R2;
    std::vector<double> xn(n1 + n2), dn(n1 + n2);
    const bool print = CurcumaLogger::get_verbosity() >= 3;
    if (print)
        CurcumaLogger::info(fmt::format("CCSD iteration   0  E = {:.12f} (MP2 amplitudes)", e));
    res.converged = false;
    int it = 0;
    for (it = 1; it <= opt.max_iter; ++it) {
        ccsdResidual(t1, t2, I, threads, R1, R2);
        double rmax = 0.0;
        // Jacobi step d = R / e (e < 0, so t + d moves against the residual); x = t + d
        for (int i = 0; i < o; ++i)
            for (int a = 0; a < v; ++a) {
                const size_t k = (size_t)i * v + a;
                rmax = std::max(rmax, std::abs(R1.data[k]));
                dn[k] = R1.data[k] / eia(I, i, a);
                xn[k] = t1.data[k] + dn[k];
            }
        for (int i = 0; i < o; ++i)
            for (int j = 0; j < o; ++j)
                for (int a = 0; a < v; ++a)
                    for (int b = 0; b < v; ++b) {
                        const size_t k = (((size_t)i * o + j) * v + a) * v + b;
                        rmax = std::max(rmax, std::abs(R2.data[k]));
                        dn[n1 + k] = R2.data[k] / (eia(I, i, a) + eia(I, j, b));
                        xn[n1 + k] = t2.data[k] + dn[n1 + k];
                    }
        const std::vector<double> x = diis.update(xn, dn);
        std::copy(x.begin(), x.begin() + n1, t1.data.begin());
        std::copy(x.begin() + n1, x.end(), t2.data.begin());
        e = ccsdEnergy(t1, t2, I, threads);
        res.max_residual = rmax;
        res.last_de = e - e_old;
        if (print)
            CurcumaLogger::info(fmt::format("CCSD iteration {:3d}  E = {:.12f}  dE = {:.2e}  max|R| = {:.2e}  {:.1f} ms",
                                            it, e, e - e_old, rmax, msSince(t0)));
        if (!std::isfinite(e)) break;
        if (std::abs(e - e_old) < opt.e_threshold && rmax < opt.r_threshold) {
            res.converged = true;
            break;
        }
        e_old = e;
    }
    res.iterations = std::min(it, opt.max_iter);
    res.e_corr = e;
    double s = 0.0;
    for (double t : t1.data) s += t * t;
    res.t1_diagnostic = o > 0 ? std::sqrt(s / (2.0 * o)) : 0.0;
    res.t_iterations_ms = msSince(t0);
}

}  // namespace qmint

// =================================================================================
// 5. QMEngine::computeCCSD -- checks, AO tensor, memory, driver
// =================================================================================

// Claude Generated (Oct 2026, WP13). Needs a converged closed-shell RHF with exact integrals.
bool QMEngine::computeCCSD(const qmint::CCSDOptions& opt, qmint::CCSDResult& res)
{
    res = qmint::CCSDResult();
    if (m_functional != QMFunctional::HF) {
        CurcumaLogger::error("CCSD: needs the HF reference (QMFunctional::HF)");
        return false;
    }
    if (riActive()) {
        // The CCSD integrals are exact; an RI-HF reference would make a mixed method that no other
        // program defines. Not combined.
        CurcumaLogger::error("CCSD: -qm.ri " + m_ri_mode + " is not supported with ccsd (exact-integral RHF reference only; use -qm.ri off)");
        return false;
    }
    if (!m_scf_converged) {
        CurcumaLogger::error("CCSD: the RHF SCF is not converged -- no CCSD energy");
        return false;
    }
    const int n = m_nbf;
    const int n_occ = m_num_electrons / 2;
    if (m_num_electrons % 2 != 0) {
        CurcumaLogger::error("CCSD: closed shell only (odd number of electrons)");
        return false;
    }
    int n_frozen = 0;
    if (opt.frozen_core)
        for (int Z : m_atoms) n_frozen += qmint::frozenCoreOrbitals(Z);
    if (n_frozen >= n_occ && n_occ > 0) {
        CurcumaLogger::error(fmt::format("CCSD: frozen core ({} orbitals) leaves no active occupied orbital ({} occupied)",
                                         n_frozen, n_occ));
        return false;
    }
    if (m_mo.rows() != n || m_mo.cols() != n || m_energies.size() != n) {
        CurcumaLogger::error("CCSD: no MO coefficients / orbital energies from the SCF");
        return false;
    }
    res.n_frozen = n_frozen;
    res.n_occ_active = n_occ - n_frozen;
    res.n_virt = n - n_occ;
    res.n_basis = n;
    if (res.n_virt == 0 || res.n_occ_active == 0) {
        res.converged = true;  // nothing to correlate: E_corr = 0
        return true;
    }

    res.memory_mb = qmint::ccsdMemoryMB(n, n_frozen, n_occ, opt.diis);
    if (res.memory_mb > opt.max_memory_mb) {
        CurcumaLogger::error(fmt::format(
            "CCSD: needs about {:.0f} MB (MO integrals incl. v^4 = {:.0f} MB, transformation buffers, amplitudes; "
            "n = {}, o = {}, v = {}), above -qm.cc_max_memory_mb {:.0f}; raise the limit or use a smaller system",
            res.memory_mb, std::pow((double)res.n_virt, 4) * 8.0 / 1.0e6, n, res.n_occ_active, res.n_virt,
            opt.max_memory_mb));
        return false;
    }
    if (gpuRequested())
        curcuma::reportGpuFallback("qm_ccsd", "CCSD has no device implementation; transformation and iterations on the CPU");

    // AO tensor: the same source as MP2 (stored SCF tensor, device download, or a rebuild within
    // -qm.mp2_max_memory_mb)
    qmint::MP2Options mopt;
    mopt.max_memory_mb = opt.ao_max_memory_mb;
    qmint::MP2Result mres;
    qmint::ERITensor local;
    std::string src;
    const double* eri = mp2HostTensor(mopt, mres, local, src);
    if (!eri) return false;
    res.t_integrals_ms = mres.t_integrals_ms;

    const int threads = std::max(1, m_threads);
    qmint::CCSDIntegrals I;
    const auto t0 = std::chrono::steady_clock::now();
    qmint::ccsdIntegrals(eri, n, m_mo, m_energies, n_frozen, n_occ, threads, I);
    res.t_transform_ms = msSince(t0);
    local = qmint::ERITensor();  // a rebuilt tensor is not needed any more

    qmint::ccsdSolve(I, opt, threads, res);
    res.info = fmt::format("{}; CPU, {} threads, estimated peak {:.0f} MB", src, threads, res.memory_mb);
    if (!std::isfinite(res.e_corr)) {
        CurcumaLogger::error("CCSD: non-finite correlation energy (diverged amplitudes)");
        return false;
    }
    if (!res.converged) {
        CurcumaLogger::error(fmt::format("CCSD: not converged after {} iterations (last dE = {:.2e}, max|R| = {:.2e}); "
                                         "raise -qm.cc_max_iter", res.iterations, res.last_de, res.max_residual));
        return false;
    }
    return true;
}
