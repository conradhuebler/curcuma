/*
 * <Canonical closed-shell MP2 on top of the native RHF (WP11) -- implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Formulas and the transformation scheme: see qm_mp2.h. This file holds
 *   - qmint::frozenCoreOrbitals   chemical frozen-core count per element
 *   - qmint::mp2EnergyHost        the O(n^5) transformation + pair-energy sum on the CPU
 *   - QMEngine::computeMP2        picks the AO tensor (stored host tensor, device tensor,
 *                                 download, or a fresh build) and the CPU or device path
 *
 * Claude Generated (Oct 2026, WP11). AI-generated, machine-tested.
 *
 * This program is free software under GPL-3.0
 */

#include "qm_mp2.h"
#include "qm_engine.h"
#include "qm_gpu_mp2_backend.h"

#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/gpu_plugin.h"
#include "src/core/gpu_fallback.h"

#include <fmt/format.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <vector>

namespace {
using RowMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using RowMap = Eigen::Map<RowMat>;
using ConstRowMap = Eigen::Map<const RowMat>;

double msSince(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}
}  // namespace

namespace qmint {

int frozenCoreOrbitals(int Z)
{
    if (Z <= 2) return 0;    // H, He
    if (Z <= 10) return 1;   // Li-Ne: 1s (also Li, Be, where ORCA/PySCF freeze nothing; kept, H-Ne results unchanged)
    // Z > 10 (Nov 2026): the ORCA 6.1 / PySCF chemcore table (all-electron core, before an ECP is subtracted, see
    // QMEngine::frozenCoreCount). ORCA 6.1 "Freezing NCore" checked for Na, Mg, Al, K, Ca, Sc, Zn, Ga, Kr, Rb, Sr,
    // Y, Cd, In, I, Cs, Ba, La, Hf, Hg, Tl, Pb with def2-SVP; PySCF pyscf.data.elements.chemcore_atm has the same.
    if (Z <= 12) return 1;   // Na, Mg: 1s
    if (Z <= 30) return 5;   // Al-Zn: [Ne]
    if (Z <= 38) return 9;   // Ga-Sr: [Ar]
    if (Z <= 48) return 14;  // Y-Cd: [Ar]3d
    if (Z <= 70) return 18;  // In-Yb: [Kr]
    if (Z <= 80) return 23;  // Lu-Hg: [Kr]4d
    if (Z <= 102) return 34; // Tl-No: [Xe]4f
    return 50;
}

// Claude Generated (Oct 2026, WP11): CPU canonical MP2, batches over the active occupied index.
void mp2EnergyHost(const double* eri, int n, const Matrix& C, const Vector& eps, int n_frozen, int n_occ,
                   int threads, double batch_mb, MP2Result& res)
{
    const int o = n_occ - n_frozen;  // active occupied
    const int v = n - n_occ;         // virtuals
    res.e_os = res.e_ss = 0.0;
    res.t_transform_ms = res.t_energy_ms = 0.0;
    if (o <= 0 || v <= 0) return;
    threads = std::max(1, threads);

    // MO blocks: Co (n x o, active occupied), Cv (n x v); transposed copies for the left GEMMs.
    const RowMat Co = C.block(0, n_frozen, n, o);
    const RowMat Cv = C.block(0, n_occ, n, v);
    const RowMat CoT = Co.transpose(), CvT = Cv.transpose();
    const size_t n2 = (size_t)n * n, n3 = n2 * n;

    // Batch size over i: (i nu|lam sig) [n^3] plus the finished (ia|jb) block [v o v] per i.
    const double per_i_mb = (double)(n3 + (size_t)o * v * v) * 8.0 / 1.0e6;
    const int batch = std::max(1, std::min(o, (int)std::floor(batch_mb / per_i_mb)));
    std::vector<double> X1((size_t)batch * n3);
    std::vector<double> K((size_t)batch * v * o * v);  // [i][a][j][b] = (ia|jb)
    std::vector<double> os_i(o, 0.0), ss_i(o, 0.0);

    // The same memory read column-major: E_cm(nu lam sig, mu) = (mu nu|lam sig), an (n^3 x n) matrix.
    // Step 1 as the tall-skinny product E_cm (n^3 x n) * Co (n x ib) runs at BLAS speed and threads
    // over the long dimension; the equivalent short-wide row-major form (ib x n)(n x n^3) ran at
    // ~2.5 GFLOPS on one thread whatever the thread count (benzene: 2.8 s).
    const Eigen::Map<const Eigen::MatrixXd> Ecm(eri, (Eigen::Index)n3, n);
    const Eigen::MatrixXd CoCm = Co;  // column-major copy (n x o)
    int nbatches = 0;
    double t_step1 = 0.0;
    for (int i0 = 0; i0 < o; i0 += batch) {
        const int ib = std::min(batch, o - i0);
        ++nbatches;
        const auto t0 = std::chrono::steady_clock::now();
        // 1. (i nu|lam sig) = sum_mu C_mu,i (mu nu|lam sig), as X1_cm (n^3 x ib) = E_cm (n^3 x n) * Co.
        //    X1_cm(rest, i) sits at rest + i*n^3, i.e. row-major X1[i][nu][lam][sig]. OpenMP over row
        //    blocks of the long dimension, one GEMM per block: the BLAS linked here may be single-threaded
        //    (this machine: MKL sequential, see WP10 doc "BLAS (Korrektur)"), so threading inside one dgemm cannot be relied on.
        {
            const Eigen::Index rows = (Eigen::Index)n3;
            // NOTE (Oct 31 2026): the block height depends on the thread count (mp2GradientHost uses fixed 8192 rows
            // since then); energies were bitwise over 1/4/8 threads on H2O and NH3, larger systems not checked. Left
            // unchanged on purpose: changing it would change the MP2 energy bits.
            const Eigen::Index chunk = std::max<Eigen::Index>(1024, (rows + 8 * threads - 1) / (8 * threads));
            const Eigen::Index nchunk = (rows + chunk - 1) / chunk;
            Eigen::Map<Eigen::MatrixXd> X(X1.data(), rows, ib);
            const Eigen::MatrixXd Cb = CoCm.middleCols(i0, ib);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
            for (Eigen::Index c = 0; c < nchunk; ++c) {
                const Eigen::Index r0 = c * chunk, len = std::min(chunk, rows - r0);
                X.middleRows(r0, len).noalias() = Ecm.middleRows(r0, len) * Cb;
            }
        }
        t_step1 += msSince(t0);
        // 2.-4. per i (OpenMP over i, single-threaded GEMMs inside):
        //   W[a][lam sig] = sum_nu Cv[nu][a] X1_i[nu][lam sig]       (v x n)(n x n^2)
        //   K_i[a][j][b]  = sum_lam,sig Co[lam][j] W_a[lam][sig] Cv[sig][b]   (o x n)(n x n)(n x v)
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
        {
            RowMat W(v, (Eigen::Index)n2), T(o, n);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
            for (int ii = 0; ii < ib; ++ii) {
                const ConstRowMap A(X1.data() + (size_t)ii * n3, n, (Eigen::Index)n2);
                W.noalias() = CvT * A;
                double* Ki = K.data() + (size_t)ii * v * o * v;
                for (int a = 0; a < v; ++a) {
                    const ConstRowMap M(W.data() + (size_t)a * n2, n, n);
                    T.noalias() = CoT * M;
                    RowMap Ka(Ki + (size_t)a * o * v, o, v);
                    Ka.noalias() = T * Cv;
                }
            }
        }
        res.t_transform_ms += msSince(t0);

        // Pair energies of this batch; (ib|ja) = K_i[b][j][a] lives in the same i block.
        const auto t1 = std::chrono::steady_clock::now();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
        for (int ii = 0; ii < ib; ++ii) {
            const int i = i0 + ii;
            const double ei = eps(n_frozen + i);
            const double* Ki = K.data() + (size_t)ii * v * o * v;
            double os = 0.0, ss = 0.0;
            for (int a = 0; a < v; ++a) {
                const double ea = eps(n_occ + a);
                for (int j = 0; j < o; ++j) {
                    const double eij_a = ei + eps(n_frozen + j) - ea;
                    const double* row = Ki + ((size_t)a * o + j) * v;  // (ia|jb), b = 0..v-1
                    for (int b = 0; b < v; ++b) {
                        const double t = row[b];
                        const double tx = Ki[((size_t)b * o + j) * v + a];  // (ib|ja)
                        const double inv_d = 1.0 / (eij_a - eps(n_occ + b));
                        os += t * t * inv_d;
                        ss += t * (t - tx) * inv_d;
                    }
                }
            }
            os_i[i] = os;
            ss_i[i] = ss;
        }
        res.t_energy_ms += msSince(t1);
    }
    // Fixed summation order over i: the result does not depend on the thread count.
    for (int i = 0; i < o; ++i) {
        res.e_os += os_i[i];
        res.e_ss += ss_i[i];
    }
    res.info += fmt::format("CPU, {} threads, {} batch(es) of <= {} occupied orbitals ({:.1f} MB per batch); "
                            "first quarter {:.1f} ms",
                            threads, nbatches, batch, batch * per_i_mb, t_step1);
}

}  // namespace qmint

// =================================================================================
// QMEngine::computeMP2 -- tensor selection and CPU / device dispatch
// =================================================================================

// Claude Generated (Oct 2026, WP11). Needs a converged closed-shell RHF (Calculation() first).
// AO tensor, in this order:
//   device path: the stored device tensor of the SCF (QMGpuJK::deviceERI, Q-G1/Q-G2), else the
//                host tensor below uploaded once;
//   host tensor: the cached active-basis tensor of the SCF; else a download of the device tensor;
//                else a fresh qmint::buildERI (same screening, same 5d transform) if
//                n^4 * 8 bytes <= opt.max_memory_mb -- an integral-direct SCF stores no tensor and
//                there is no integral-direct MP2 yet, so above that limit MP2 stops with an error.
// Claude Generated (Oct 2026): host AO tensor for MP2 (energy and gradient), formerly a lambda in computeMP2.
const double* QMEngine::mp2HostTensor(const qmint::MP2Options& opt, qmint::MP2Result& res, qmint::ERITensor& local,
                                      std::string& tensor_src)
{
    const int n = m_nbf;
    const double tensor_mb = std::pow((double)n, 4) * 8.0 / 1.0e6;
    const auto t0 = std::chrono::steady_clock::now();
    const bool cached = (m_Q.size() == 0) ? m_eri_ready : m_eri_active_ready;
    if (cached) {
        tensor_src = "stored host tensor of the SCF";
        return eriActive().data();
    }
    if (tensor_mb > opt.max_memory_mb) {
        CurcumaLogger::error(fmt::format(
            "MP2: the AO tensor ({:.1f} MB) is not stored (integral-direct SCF) and exceeds -qm.mp2_max_memory_mb "
            "{:.1f}; integral-direct MP2 is not implemented", tensor_mb, opt.max_memory_mb));
        return nullptr;
    }
    if (m_gpu_jk && m_gpu_jk->n() == n) {
        local = qmint::ERITensor(n);
        if (m_gpu_jk->downloadERI(local.data())) {
            res.t_integrals_ms += msSince(t0);
            tensor_src = fmt::format("device tensor downloaded ({:.1f} MB)", tensor_mb);
            return local.data();
        }
    }
    local = qmint::buildERI(m_gto_basis, m_threads, m_eri_screening, m_Q.size() != 0 ? &m_Q : nullptr);
    res.t_integrals_ms += msSince(t0);
    tensor_src = fmt::format("AO tensor rebuilt for MP2 ({:.1f} MB, {:.1f} ms)", tensor_mb, res.t_integrals_ms);
    return local.data();
}

bool QMEngine::computeMP2(const qmint::MP2Options& opt, qmint::MP2Result& res)
{
    res = qmint::MP2Result();
    // WP12 stage 1 (Oct 2026): the double hybrids evaluate the same PT2 expression on their KS orbitals and
    // orbital energies (the HF path is unchanged)
    if (m_functional != QMFunctional::HF && m_functional != QMFunctional::DoubleHybrid) {
        CurcumaLogger::error("MP2: needs the HF or double-hybrid KS reference");
        return false;
    }
    if (!m_scf_converged) {
        CurcumaLogger::error("MP2: the RHF SCF is not converged -- no MP2 energy");
        return false;
    }
    const int n = m_nbf;
    const int n_occ = m_num_electrons / 2;
    int n_frozen = 0;
    if (opt.frozen_core)
        n_frozen += frozenCoreCount();
    if (n_frozen >= n_occ && n_occ > 0) {
        CurcumaLogger::error(fmt::format("MP2: frozen core ({} orbitals) leaves no active occupied orbital ({} occupied)",
                                         n_frozen, n_occ));
        return false;
    }
    if (m_mo.rows() != n || m_mo.cols() != n || m_energies.size() != n) {
        CurcumaLogger::error("MP2: no MO coefficients / orbital energies from the SCF");
        return false;
    }
    res.n_frozen = n_frozen;
    res.n_occ_active = n_occ - n_frozen;
    res.n_virt = n - n_occ;
    res.n_basis = n;
    if (res.n_virt == 0 || res.n_occ_active == 0) return true;  // nothing to correlate: E_corr = 0

    qmint::ERITensor local;  // temporary host tensor (download or rebuild), freed on return
    std::string tensor_src;
    // at most one acquisition per call (a device fallback must not rebuild the tensor)
    const double* host_eri = nullptr;
    bool host_tried = false;
    auto hostTensor = [&]() -> const double* {
        if (!host_tried) {
            host_tried = true;
            host_eri = mp2HostTensor(opt, res, local, tensor_src);
        }
        return host_eri;
    };

    // ---- device path (WP11 step 1, GPU variant) ----
    bool done = false;
    const std::string mode = opt.gpu_mp2;
    if (mode != "auto" && mode != "on" && mode != "off" && mode != "true" && mode != "false"
        && CurcumaLogger::get_verbosity() >= 1)
        CurcumaLogger::warn(fmt::format("MP2: unknown -qm.gpu_mp2 '{}' (use auto|on|off); using auto", mode));
    if (gpuRequested() && mode != "off" && mode != "false") {
        const std::string backend = m_gpu;
        if (backend == "none" || !gpu_plugin::available(backend)) {
            curcuma::reportGpuFallback("qm_mp2", "no GPU plugin for -gpu " + m_gpu);
        } else {
            qmint::QMGpuMP2Input in;
            in.n = n;
            in.n_frozen = n_frozen;
            in.n_occ = n_occ;
            const Matrix Cc = m_mo;  // row-major contiguous
            const Vector ec = m_energies;
            in.C = Cc.data();
            in.eps = ec.data();
            in.device = m_gpu_device;
            in.max_device_mb = m_gpu_max_memory_mb;
            in.batch_mb = opt.batch_mb;
            in.eri_device = (m_gpu_jk && m_gpu_jk->n() == n) ? m_gpu_jk->deviceERI() : nullptr;
            if (in.eri_device)
                tensor_src = "stored device tensor of the SCF (" + m_gpu_jk->name() + ")";
            else
                in.eri_host = hostTensor();
            std::string why;
            qmint::QMGpuMP2Output out;
            if ((in.eri_device || in.eri_host) && gpu_plugin::qmMp2Energy(backend, in, out, why)) {
                res.e_os = out.e_os;
                res.e_ss = out.e_ss;
                res.t_integrals_ms += out.t_upload_ms;
                res.t_transform_ms = out.t_transform_ms;
                res.t_energy_ms = out.t_energy_ms;
                res.on_device = true;
                res.info = fmt::format("{} ({}; {} batch(es) of <= {} occupied orbitals; upload {:.1f} ms)",
                                       backend, tensor_src, out.nbatches, out.batch, out.t_upload_ms);
                done = true;
            } else if (in.eri_device || in.eri_host) {
                curcuma::reportGpuFallback("qm_mp2", "device MP2 failed: " + why);
            } else {
                return false;  // no tensor at all (error printed by hostTensor)
            }
        }
    }

    // ---- CPU path ----
    if (!done) {
        const double* eri = hostTensor();
        if (!eri) return false;
        res.info = tensor_src + "; ";
        qmint::mp2EnergyHost(eri, n, m_mo, m_energies, n_frozen, n_occ, m_threads, opt.batch_mb, res);
    }
    if (!std::isfinite(res.e_os) || !std::isfinite(res.e_ss)) {
        CurcumaLogger::error("MP2: non-finite correlation energy (degenerate HOMO/LUMO?)");
        return false;
    }
    return true;
}
