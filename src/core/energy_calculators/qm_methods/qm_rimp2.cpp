/*
 * <RI-MP2: density-fitted closed-shell MP2 on top of the native RHF>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Resolution-of-the-identity MP2 (Feyereisen, Fitzgerald, Komornicki, Chem. Phys. Lett. 208, 359
 * (1993); Weigend, Haeser, Theor. Chem. Acc. 97, 331 (1997)). The MO integrals of canonical MP2
 * (qm_mp2.h) are replaced by the Coulomb-metric fit
 *
 *   (ia|jb) ~ sum_P B_ia^P B_jb^P,   B_ia^P = sum_Q [L^-1]_PQ (Q|ia),   J = L L^T,  J_PQ = (P|Q),
 *   (Q|ia)  = sum_mn C_mi (Q|mn) C_na,
 *
 * and the energy expression is the canonical one (E_OS, E_SS, SCS, SOS; same frozen-core rule).
 * Steps here:
 *   1. RI tensor: the SCF's own one when -qm.ri j|jk ran with the same auxiliary file, else a
 *      separate qmri::RIJK setup in mode j (raw (Q|mn) + metric factor, no n^2-column forward
 *      substitution) with -qm.ri_mp2_aux.
 *   2. (Q|ia) = Co^T T^Q Cv per auxiliary row (two gemmABt calls, OpenMP over Q), then the metric
 *      factor on the o*v columns (RIJK::fitRows), then the transpose to B[i][a][P].
 *   3. Pair energies: CPU qmint::riMp2EnergyHost (pairs i >= j, one v x v GEMM each) or the device
 *      (cuda/qm_gpu_rimp2.cu, batches over i, cuBLAS DGEMM), CPU fallback counted as qm_ri_mp2.
 *
 * Auxiliary basis caveat: curcuma ships no def2-SVP/C file; the default correlation fitting basis
 * is def2-universal-JKFIT (or -JFIT). That is a method variant, not RI-MP2/def2-SVP/C as in other
 * programs. Energy only (no RI-MP2 gradient).
 *
 * Claude Generated (Oct 2026). AI-generated, machine-tested.
 *
 * This program is free software under GPL-3.0
 */

#include "qm_engine.h"
#include "qm_gpu_rimp2_backend.h"
#include "qm_mp2.h"
#include "qm_ri.h"

#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/gpu_plugin.h"
#include "src/core/gpu_fallback.h"
#include "src/core/units.h"

#include <fmt/format.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <vector>

namespace {
double msSince(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

/// File name without directory and extension (label of the auxiliary set).
std::string auxLabel(const std::string& path)
{
    std::string s = path;
    const size_t slash = s.find_last_of('/');
    if (slash != std::string::npos) s = s.substr(slash + 1);
    const size_t dot = s.rfind(".turbomole");
    if (dot != std::string::npos) s = s.substr(0, dot);
    return s;
}
}  // namespace

namespace qmint {

// Claude Generated (Oct 2026, RI-MP2). C = A B^T with 4 x 4 register blocks of dot products along
// the contiguous k dimension (the compiler vectorises the q loop with -march=native).
void gemmABt(const double* A, int m, const double* B, int n, int k, double* C)
{
    int i = 0;
    for (; i + 4 <= m; i += 4) {
        const double* a0 = A + (size_t)i * k;
        const double *a1 = a0 + k, *a2 = a1 + k, *a3 = a2 + k;
        int j = 0;
        for (; j + 4 <= n; j += 4) {
            const double* b0 = B + (size_t)j * k;
            const double *b1 = b0 + k, *b2 = b1 + k, *b3 = b2 + k;
            double s00 = 0, s01 = 0, s02 = 0, s03 = 0, s10 = 0, s11 = 0, s12 = 0, s13 = 0;
            double s20 = 0, s21 = 0, s22 = 0, s23 = 0, s30 = 0, s31 = 0, s32 = 0, s33 = 0;
#pragma omp simd reduction(+ : s00, s01, s02, s03, s10, s11, s12, s13, s20, s21, s22, s23, s30, s31, s32, s33)
            for (int q = 0; q < k; ++q) {
                const double x0 = a0[q], x1 = a1[q], x2 = a2[q], x3 = a3[q];
                const double y0 = b0[q], y1 = b1[q], y2 = b2[q], y3 = b3[q];
                s00 += x0 * y0; s01 += x0 * y1; s02 += x0 * y2; s03 += x0 * y3;
                s10 += x1 * y0; s11 += x1 * y1; s12 += x1 * y2; s13 += x1 * y3;
                s20 += x2 * y0; s21 += x2 * y1; s22 += x2 * y2; s23 += x2 * y3;
                s30 += x3 * y0; s31 += x3 * y1; s32 += x3 * y2; s33 += x3 * y3;
            }
            double* c = C + (size_t)i * n + j;
            c[0] = s00; c[1] = s01; c[2] = s02; c[3] = s03;
            c += n;
            c[0] = s10; c[1] = s11; c[2] = s12; c[3] = s13;
            c += n;
            c[0] = s20; c[1] = s21; c[2] = s22; c[3] = s23;
            c += n;
            c[0] = s30; c[1] = s31; c[2] = s32; c[3] = s33;
        }
        for (; j < n; ++j)
            for (int r = 0; r < 4; ++r) {
                const double* a = A + (size_t)(i + r) * k;
                const double* b = B + (size_t)j * k;
                double s = 0.0;
                for (int q = 0; q < k; ++q) s += a[q] * b[q];
                C[(size_t)(i + r) * n + j] = s;
            }
    }
    for (; i < m; ++i)
        for (int j = 0; j < n; ++j) {
            const double* a = A + (size_t)i * k;
            const double* b = B + (size_t)j * k;
            double s = 0.0;
            for (int q = 0; q < k; ++q) s += a[q] * b[q];
            C[(size_t)i * n + j] = s;
        }
}

// Claude Generated (Oct 2026, RI-MP2). Pair (i,j), i >= j: K_ab = (ia|jb) = sum_P B_ia^P B_jb^P,
// (ib|ja) = K_ba. The (j,i) pair gives K^T with the same denominators, hence the factor 2 for i != j.
void riMp2EnergyHost(const double* B, int o, int v, int naux, const double* eo, const double* ev, int threads,
                     MP2Result& res)
{
    res.e_os = res.e_ss = 0.0;
    if (o <= 0 || v <= 0) return;
    threads = std::max(1, threads);
    const int npairs = o * (o + 1) / 2;
    std::vector<int> pi(npairs), pj(npairs);
    for (int i = 0, p = 0; i < o; ++i)
        for (int j = 0; j <= i; ++j, ++p) {
            pi[p] = i;
            pj[p] = j;
        }
    std::vector<double> pos(npairs, 0.0), pss(npairs, 0.0);
    const auto t0 = std::chrono::steady_clock::now();
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
    {
        std::vector<double> K((size_t)v * v);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (int p = 0; p < npairs; ++p) {
            const int i = pi[p], j = pj[p];
            gemmABt(B + (size_t)i * v * naux, v, B + (size_t)j * v * naux, v, naux, K.data());
            const double eij = eo[i] + eo[j];
            double os = 0.0, ss = 0.0;
            for (int a = 0; a < v; ++a) {
                const double eija = eij - ev[a];
                const double* row = K.data() + (size_t)a * v;
                for (int b = 0; b < v; ++b) {
                    const double t = row[b];
                    const double tx = K[(size_t)b * v + a];
                    const double inv_d = 1.0 / (eija - ev[b]);
                    os += t * t * inv_d;
                    ss += t * (t - tx) * inv_d;
                }
            }
            const double f = (i == j) ? 1.0 : 2.0;
            pos[p] = f * os;
            pss[p] = f * ss;
        }
    }
    for (int p = 0; p < npairs; ++p) {  // fixed order: independent of the thread count
        res.e_os += pos[p];
        res.e_ss += pss[p];
    }
    res.t_energy_ms += msSince(t0);
}

}  // namespace qmint

// =================================================================================
// QMEngine::correlationRI (shared by RI-MP2 and LMP2)
// =================================================================================

// Claude Generated (Oct 2026): moved out of computeRIMP2 so that LMP2 (qm_lmp2.cpp) uses the same tensor.
std::shared_ptr<qmri::RIJK> QMEngine::correlationRI(const std::string& file, std::string& src, double& setup_ms,
                                                    std::string& why)
{
    setup_ms = 0.0;
    if (riActive() && m_ri && m_ri_file == file) {
        src = "RI tensor of the SCF reused (-qm.ri " + m_ri_mode + ")";
        return m_ri;
    }
    const auto t0 = std::chrono::steady_clock::now();
    std::shared_ptr<qmri::RIJK> ri;
    try {
        const qmri::AuxBasisMap aux = qmri::parseTurbomoleBasis(file);
        Matrix pos = Matrix::Zero(m_atomcount, 3);
        for (int i = 0; i < m_atomcount; ++i)
            for (int k = 0; k < 3; ++k) pos(i, k) = CurcumaUnit::Length::angstrom_to_bohr(m_geometry(i, k));
        const std::vector<qmri::Shell> sh = qmri::auxiliaryShells(aux, m_atoms, pos);
        qmri::RIOptions ro;
        ro.mode = "j";  // raw (Q|mn) + metric factor; the factor goes onto the transformed columns (fitRows)
        ro.metric = m_ri_metric;
        ro.metric_cutoff = m_ri_metric_cutoff;
        ro.screening = m_eri_screening;
        ro.max_memory_mb = m_ri_max_memory_mb;
        ro.spherical_aux = !m_cartesian_d;
        ro.threads = m_threads;
        auto r = std::make_shared<qmri::RIJK>();
        if (r->setup(m_gto_basis, m_Q.size() != 0 ? &m_Q : nullptr, sh, ro, why)) ri = r;
    } catch (const std::exception& e) {
        why = e.what();
    }
    if (!ri) return nullptr;
    setup_ms = msSince(t0);
    src = "own RI setup: " + ri->stats();
    return ri;
}

// =================================================================================
// QMEngine::computeRIMP2
// =================================================================================

// Claude Generated (Oct 2026, RI-MP2). Needs a converged closed-shell RHF (Calculation() first).
bool QMEngine::computeRIMP2(const qmint::MP2Options& opt, qmint::MP2Result& res)
{
    res = qmint::MP2Result();
    res.ri = true;
    if (m_functional != QMFunctional::HF) {
        CurcumaLogger::error("RI-MP2: needs the HF reference (QMFunctional::HF)");
        return false;
    }
    if (!m_scf_converged) {
        CurcumaLogger::error("RI-MP2: the RHF SCF is not converged -- no MP2 energy");
        return false;
    }
    const int n = m_nbf;
    const int n_occ = m_num_electrons / 2;
    int n_frozen = 0;
    if (opt.frozen_core)
        for (int Z : m_atoms) n_frozen += qmint::frozenCoreOrbitals(Z);
    if (n_frozen >= n_occ && n_occ > 0) {
        CurcumaLogger::error(fmt::format("RI-MP2: frozen core ({} orbitals) leaves no active occupied orbital ({} occupied)",
                                         n_frozen, n_occ));
        return false;
    }
    if (m_mo.rows() != n || m_mo.cols() != n || m_energies.size() != n) {
        CurcumaLogger::error("RI-MP2: no MO coefficients / orbital energies from the SCF");
        return false;
    }
    const int o = n_occ - n_frozen, v = n - n_occ;
    res.n_frozen = n_frozen;
    res.n_occ_active = o;
    res.n_virt = v;
    res.n_basis = n;
    if (riActive() && m_ri) res.scf_ri = fmt::format("RI-{} {}", m_ri_mode == "j" ? "J" : "JK", auxLabel(m_ri_file));

    // ---- 1. RI tensor: reuse the SCF's one for the same auxiliary file, else a mode-j setup ----
    std::string file;
    try {
        file = qmri::resolveAuxBasisFile(opt.ri_aux, "jk");
    } catch (const std::exception& e) {
        CurcumaLogger::error(std::string("RI-MP2: ") + e.what());
        return false;
    }
    res.aux_name = auxLabel(file);
    if (o == 0 || v == 0) return true;  // nothing to correlate

    std::string ri_src, why;
    double setup_ms = 0.0;
    const std::shared_ptr<qmri::RIJK> ri = correlationRI(file, ri_src, setup_ms, why);
    if (!ri) {
        CurcumaLogger::error("RI-MP2: density fitting setup with " + file + " not possible: " + why);
        return false;
    }
    res.t_ri_setup_ms = setup_ms;

    // ---- 2. B[i][a][P] ----
    const int nrows = (int)ri->tensor().rows();
    const size_t ov = (size_t)o * v;
    const double b_mb = (double)ov * nrows * 8.0 / 1.0e6;
    if (2.0 * b_mb > m_ri_max_memory_mb) {
        CurcumaLogger::error(fmt::format("RI-MP2: the (ia|P) tensor needs 2 x {:.0f} MB (o v naux = {} x {} x {}, plus one "
                                         "transposed copy), above -qm.ri_max_memory_mb {:g}",
                                         b_mb, o, v, nrows, m_ri_max_memory_mb));
        return false;
    }
    const auto t1 = std::chrono::steady_clock::now();
    const int threads = std::max(1, m_threads);
    Matrix X(nrows, (Eigen::Index)ov);  // (Q|ia), row Q = [i][a]
    {
        const Matrix Co = m_mo.block(0, n_frozen, n, o);
        const Matrix Cv = m_mo.block(0, n_occ, n, v);
        const Matrix CoT = Co.transpose(), CvT = Cv.transpose();  // row-major o x n, v x n
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
        {
            std::vector<double> Y((size_t)o * n), TQ((size_t)n * n);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
            for (int Q = 0; Q < nrows; ++Q) {
                ri->unpackRow(Q, TQ.data());  // the tensor stores the packed lower triangle per row (qm_ri.h)
                // T^Q is symmetric: (Co^T T^Q)[i][m] = sum_k CoT[i][k] T^Q[m][k]
                qmint::gemmABt(CoT.data(), o, TQ.data(), n, n, Y.data());
                // (Q|ia) = sum_m Y[i][m] Cv[m][a] = sum_m Y[i][m] CvT[a][m]
                qmint::gemmABt(Y.data(), o, CvT.data(), v, n, X.data() + (size_t)Q * ov);
            }
        }
    }
    ri->fitRows(X);  // metric factor (no-op when the SCF tensor already is B)
    const int naux = (int)X.rows();
    std::vector<double> B(ov * naux);
    {
        const size_t blk = 64;
        const size_t nbq = ((size_t)naux + blk - 1) / blk, nbc = (ov + blk - 1) / blk;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) collapse(2)
#endif
        for (size_t bq = 0; bq < nbq; ++bq)
            for (size_t bc = 0; bc < nbc; ++bc)
                for (size_t q = bq * blk; q < std::min((size_t)naux, (bq + 1) * blk); ++q)
                    for (size_t c = bc * blk; c < std::min(ov, (bc + 1) * blk); ++c)
                        B[c * naux + q] = X(q, c);
    }
    Matrix().swap(X);
    res.t_transform_ms = msSince(t1);
    res.n_aux = naux;

    std::vector<double> eo(o), ev(v);
    for (int i = 0; i < o; ++i) eo[i] = m_energies(n_frozen + i);
    for (int a = 0; a < v; ++a) ev[a] = m_energies(n_occ + a);

    // ---- 3. pair energies: device (-qm.gpu_mp2) or CPU ----
    bool done = false;
    const std::string mode = opt.gpu_mp2;
    if (mode != "auto" && mode != "on" && mode != "off" && mode != "true" && mode != "false"
        && CurcumaLogger::get_verbosity() >= 1)
        CurcumaLogger::warn(fmt::format("RI-MP2: unknown -qm.gpu_mp2 '{}' (use auto|on|off); using auto", mode));
    if (gpuRequested() && mode != "off" && mode != "false") {
        const std::string backend = (m_gpu == "auto") ? gpu_plugin::firstAvailable() : m_gpu;
        if (backend == "none" || !gpu_plugin::available(backend)) {
            curcuma::reportGpuFallback("qm_ri_mp2", "no GPU plugin for -gpu " + m_gpu);
        } else {
            qmint::QMGpuRIMP2Input in;
            in.o = o;
            in.v = v;
            in.naux = naux;
            in.B = B.data();
            in.eo = eo.data();
            in.ev = ev.data();
            in.device = m_gpu_device;
            in.max_device_mb = m_gpu_max_memory_mb;
            in.batch_mb = opt.batch_mb;
            qmint::QMGpuRIMP2Output out;
            std::string why;
            if (gpu_plugin::qmRiMp2Energy(backend, in, out, why)) {
                res.e_os = out.e_os;
                res.e_ss = out.e_ss;
                res.t_energy_ms = out.t_upload_ms + out.t_gemm_ms + out.t_energy_ms;
                res.on_device = true;
                res.info = fmt::format("{}; energy on {} ({} batch(es) of <= {} occupied orbitals; upload {:.1f} ms, "
                                       "GEMM {:.1f} ms, sum {:.1f} ms)",
                                       ri_src, backend, out.nbatches, out.batch, out.t_upload_ms, out.t_gemm_ms,
                                       out.t_energy_ms);
                done = true;
            } else {
                curcuma::reportGpuFallback("qm_ri_mp2", "device RI-MP2 failed: " + why);
            }
        }
    }
    if (!done) {
        qmint::riMp2EnergyHost(B.data(), o, v, naux, eo.data(), ev.data(), threads, res);
        res.info = fmt::format("{}; energy on CPU, {} threads, {} pairs", ri_src, threads, o * (o + 1) / 2);
    }
    if (!std::isfinite(res.e_os) || !std::isfinite(res.e_ss)) {
        CurcumaLogger::error("RI-MP2: non-finite correlation energy (degenerate HOMO/LUMO?)");
        return false;
    }
    return true;
}
