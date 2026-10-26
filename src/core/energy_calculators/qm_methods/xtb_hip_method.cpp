/*
 * <Native xTB ROCm/HIP Method Wrapper — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06, restructured Sep 2026): everything the three GPU wrappers
 * share — the ComputationalMethod forwarding, the device handshake/logging/CPU fallback,
 * and the GpuScfBackend marshalling — now lives in xtb_gpu_adapter.h. This file keeps
 * ONLY what is specific to ROCm:
 *   - the eigensolver hook (rocSOLVER dsygvd solves the GENERALIZED problem directly, so
 *     the host rebuilds S = L·Lᵀ; Vulkan instead reduces to a standard problem on the host),
 *   - two GpuScfBackend overrides (the device-built Cholesky download, the d-shell flag),
 *   - the mixed-precision policy (ON for ROCm),
 *   - the stage banner + the no-rocSOLVER Stage-0 note.
 * Host-compiled (no device code) — the HIP kernels live in the .hip TU.
 */

#ifdef USE_ROCM

#include "xtb_hip_method.h"

namespace {
using curcuma::xtb::MethodType;
using curcuma::xtb::gpu::XtbHipContext;

// WP8b (Claude Generated, Oct 2026): defaults of -xtb.rocm_grad_resident, -xtb.rocm_grad_tail_device and
// -xtb.rocm_defer_host_integrals (keep equal to the PARAM defaults in native_xtb_method.h).
constexpr bool kRocmGradResidentDefault = true;
constexpr bool kRocmGradTailDefault = true;
constexpr bool kRocmDeferHostDefault = true;
// WP10 (Claude Generated, Oct 2026): default of -xtb.rocm_fixed_point (keep equal to the PARAM in native_xtb_method.h).
constexpr bool kRocmFixedPointDefault = true;
// WP11 (Claude Generated, Oct 2026): default of -xtb.rocm_otf_multipole (keep equal to the PARAM in native_xtb_method.h).
constexpr bool kRocmOtfMultipoleDefault = true;

#ifdef HAVE_ROCSOLVER
/**
 * @brief Device-resident GFN1/GFN2 SCF backend over an XtbHipContext.
 *
 * All ~40 GpuScfBackend virtuals are forwarded by the shared XtbGpuResidentBackend
 * template (xtb_gpu_adapter.h). Only the two genuinely ROCm-specific pieces are here.
 * Claude Generated (Stage 2/R-AP1/R-AP2/R-AP3; restructured Sep 2026).
 */
class HipScfBackend
    : public curcuma::xtb::gpu::XtbGpuResidentBackend<XtbHipContext,
                                                     curcuma::xtb::gpu::XtbHipBasisData> {
public:
    explicit HipScfBackend(XtbHipContext* ctx)
        : XtbGpuResidentBackend(ctx) {}

    // ---- BACKEND-SPECIFIC 1: the Cholesky factor -------------------------
    // ROCm builds L = chol(S) on the device (rocSOLVER) and downloads it. Vulkan has no
    // Cholesky shader and computes it host-side from the device S, so this cannot live
    // in the shared template.
    bool downloadCholesky(Eigen::MatrixXd& L_out) override
    {
        if (!m_ctx || m_n <= 0) return false;
        L_out.resize(m_n, m_n);
        return m_ctx->downloadCholesky(L_out.data());
    }

    // ---- BACKEND-SPECIFIC 2: d shells ------------------------------------
    // X-I1 B6: the HIP integral/SCF/gradient kernels handle d shells (cartesian->spherical
    // dtrafo in rocm/xtb_hip_integrals.hiph + the dpair branch in k_overlap_h0 /
    // k_multipole_ints / k_grad_h0_pulay). The Vulkan GLSL shaders do not yet, so they
    // keep the base false (d systems route to the CPU path). Claude Generated.
    bool supportsDshell() const override { return true; }

    // ---- BACKEND-SPECIFIC 3: multipole potential of the host-driven loop --------------
    // Claude Generated (Oct 2026), -xtb.rocm_device_multipole (default false). The GFN2
    // multipole potential and energy (XTB::addMultipolePotential / energyMultipole) run on the
    // device from the host interaction matrices, uploaded once per setupMultipole(). Vulkan has
    // no such kernels, so this is not in the shared template.
    void setDeviceMultipole(bool on) { m_device_mp = on; }
    bool supportsDeviceMultipolePotential() const override { return m_device_mp; }

    /// Device limit for the 15 nat^2 doubles of the matrices (4 GB, nat up to about 5770).
    static constexpr double kMpMaxBytes = 4.0e9;

    bool beginMultipolePotential(const std::array<Eigen::MatrixXd, 3>& sd,
                                 const std::array<std::array<Eigen::MatrixXd, 3>, 3>& dd,
                                 const std::array<Eigen::MatrixXd, 6>& sq,
                                 const std::vector<double>& dkernel,
                                 const std::vector<double>& qkernel,
                                 unsigned long long version) override
    {
        if (!m_ctx || m_nat <= 0) return false;
        const Eigen::Index n = m_nat;
        auto okm = [n](const Eigen::MatrixXd& m) { return m.rows() == n && m.cols() == n; };
        const double* mats[15] = {
            sd[0].data(), sd[1].data(), sd[2].data(),
            dd[0][0].data(), dd[0][1].data(), dd[0][2].data(), dd[1][1].data(), dd[1][2].data(), dd[2][2].data(),
            sq[0].data(), sq[1].data(), sq[2].data(), sq[3].data(), sq[4].data(), sq[5].data() };
        bool shapes = static_cast<Eigen::Index>(dkernel.size()) == n
                   && static_cast<Eigen::Index>(qkernel.size()) == n;
        for (int k = 0; k < 3; ++k) shapes = shapes && okm(sd[k]);
        for (int a = 0; a < 3; ++a) for (int b = a; b < 3; ++b) shapes = shapes && okm(dd[a][b]);
        for (int k = 0; k < 6; ++k) shapes = shapes && okm(sq[k]);
        if (!shapes) { m_mp_error = "matrix shapes do not match the atom count"; return false; }
        if (!m_ctx->beginMultipolePotential(m_nat, mats, dkernel.data(), qkernel.data(), version,
                                            kMpMaxBytes)) {
            m_mp_error = m_ctx->multipolePotentialError();
            return false;
        }
        return true;
    }
    bool multipolePotential(const Vector& q_at, const Eigen::MatrixXd& dp_at,
                            const Eigen::MatrixXd& qp_at, Vector& v_at,
                            Eigen::MatrixXd& v_dp, Eigen::MatrixXd& v_qp) override
    {
        if (!m_ctx || m_nat <= 0 || q_at.size() != m_nat
            || dp_at.rows() != 3 || dp_at.cols() != m_nat
            || qp_at.rows() != 6 || qp_at.cols() != m_nat) return false;
        v_at.resize(m_nat);
        v_dp.resize(3, m_nat);
        v_qp.resize(6, m_nat);
        return m_ctx->multipolePotential(q_at.data(), dp_at.data(), qp_at.data(),
                                         v_at.data(), v_dp.data(), v_qp.data());
    }
    bool multipoleEnergy(const Vector& q_at, const Eigen::MatrixXd& dp_at,
                         const Eigen::MatrixXd& qp_at, double& e) override
    {
        if (!m_ctx || m_nat <= 0 || q_at.size() != m_nat
            || dp_at.rows() != 3 || dp_at.cols() != m_nat
            || qp_at.rows() != 6 || qp_at.cols() != m_nat) return false;
        return m_ctx->multipoleEnergy(q_at.data(), dp_at.data(), qp_at.data(), &e);
    }
    bool solveMultipoleResident(const Eigen::VectorXd& v_ao, Vector& eps, bool fp32) override
    {
        if (!m_ctx || static_cast<int>(v_ao.size()) != m_n) return false;
        eps.resize(m_n);
        return m_ctx->solveMultipoleResident(v_ao.data(), eps.data(), fp32);
    }
    void noteDeviceMultipoleFallback(const std::string& why) override
    {
        if (!m_mp_fallback.empty()) return;   // keep the first event
        m_mp_fallback = why + (m_mp_error.empty() ? std::string() : " (" + m_mp_error + ")");
    }
    /// First device-multipole fallback of this backend ("" = none). Claude Generated (Oct 2026).
    const std::string& deviceMultipoleFallback() const { return m_mp_fallback; }

    // ---- BACKEND-SPECIFIC 4: GFN2 SCF potential on the device (WP1) -------------------
    // Claude Generated (Oct 2026), -xtb.rocm_device_potential. The host SCF loop uploads the mixed
    // q_sh/dp_at/qp_at and the host D4 reference weights; the device builds the whole potential
    // (gamma q_sh, third order, multipole, D4 dE/dq), the Fock and the eigensolve
    // (XtbHipContext::solvePotential). The interaction matrices are built on the device
    // (beginPotentialStored) or, if that fails, uploaded from the host (beginPotential). A failed
    // setup or SCF step moves this object to the host potential build for the rest of the run
    // (one notice, devicePotentialFallback()). Vulkan has no such kernels.
    void setDevicePotential(bool on) { m_device_pot = on; }
    bool supportsDevicePotential() const override { return m_device_pot && m_pot_fallback.empty(); }
    bool supportsStoredMultipoleOnDevice() const override { return supportsDevicePotential(); }

    bool beginBasis(const curcuma::xtb::GpuBasisFlat& bf, const curcuma::xtb::GpuH0Flat& hf) override
    {
        m_sh2at = bf.sh2at;   // host copy for q_at = sum of q_sh per atom in solvePotential
        return XtbGpuResidentBackend::beginBasis(bf, hf);
    }
    bool beginPotential(int nat, int nsh, const double* amat_sd, const double* amat_dd,
                        const double* amat_sq, const double* dkernel, const double* qkernel,
                        const double* gamma3) override
    {
        m_res_error.clear();
        if (!m_ctx || !supportsDevicePotential()) return false;
        return notePotentialBegin(m_ctx->beginPotential(nat, nsh, amat_sd, amat_dd, amat_sq, dkernel,
                                                        qkernel, gamma3, kMpMaxBytes));
    }
    bool beginPotentialStored(int nat, int nsh, const double* xyz_bohr, const double* mrad,
                              double dmp3, double dmp5, const double* dkernel, const double* qkernel,
                              const double* gamma3) override
    {
        m_res_error.clear();
        if (!m_ctx || !supportsDevicePotential()) return false;
        return notePotentialBegin(m_ctx->beginPotentialStored(nat, nsh, xyz_bohr, mrad, dmp3, dmp5,
                                                              dkernel, qkernel, gamma3, kMpMaxBytes));
    }
    // WP11 (Claude Generated, Oct 2026), -xtb.rocm_otf_multipole: the device potential without the stored interaction
    // matrices (XtbHipContext::beginPotentialOnTheFly). In the automatic choice of the host (-gpu_multipole_otf auto)
    // preferStoredMultipole keeps the stored path while its 15 nat^2 doubles fit kMpMaxBytes, so only systems above
    // the limit (about 5770 atoms), which took the host potential before, change path. A failed begin takes the
    // host potential build like a failed stored begin (notePotentialBegin, devicePotentialFallback).
    void setOnTheFlyMultipole(bool on) { m_otf_mp = on; }
    bool supportsOnTheFlyMultipole() const override { return m_otf_mp && supportsDevicePotential(); }
    bool preferStoredMultipole(int nat) const override
    {
        return 15.0 * static_cast<double>(nat) * static_cast<double>(nat) * sizeof(double) <= kMpMaxBytes;
    }
    bool beginPotentialOnTheFly(int nat, int nsh, const double* xyz_bohr, const double* mrad, double dmp3,
                                double dmp5, const double* dkernel, const double* qkernel,
                                const double* gamma3) override
    {
        m_res_error.clear();
        if (!m_ctx || !supportsOnTheFlyMultipole()) return false;
        const bool ok = notePotentialBegin(m_ctx->beginPotentialOnTheFly(nat, nsh, xyz_bohr, mrad, dmp3, dmp5,
                                                                         dkernel, qkernel, gamma3));
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info(ok ? fmt::format("ROCm: GFN2 multipole interaction rebuilt on the fly in every step "
                                                 "({} atoms, no stored matrices, -xtb.rocm_otf_multipole)", nat)
                                   : "ROCm: on-the-fly multipole setup failed (" + m_pot_begin_failed + ")");
        return ok;
    }
    bool solvePotential(const Vector& q_sh, const Eigen::MatrixXd& dp_at, const Eigen::MatrixXd& qp_at,
                        const std::vector<double>& W, const std::vector<double>& dWq, Vector& eps,
                        bool fp32 = false, int n_eig = 0) override
    {
        (void)n_eig;   // full spectrum (the caller disables the partial window on this path)
        m_res_error.clear();
        if (!m_ctx || m_nat <= 0 || q_sh.size() != m_nsh || static_cast<int>(m_sh2at.size()) != m_nsh
            || dp_at.rows() != 3 || dp_at.cols() != m_nat || qp_at.rows() != 6 || qp_at.cols() != m_nat
            || W.size() < 7 * static_cast<size_t>(m_nat) || dWq.size() < 7 * static_cast<size_t>(m_nat))
            return false;
        // q_at in the order of the host unpackSCC (xtb_native.cpp), so both paths see the same charges.
        std::vector<double> q_at(m_nat, 0.0);
        for (int s = 0; s < m_nsh; ++s) q_at[m_sh2at[s]] += q_sh(s);
        eps.resize(m_n);
        return m_ctx->solvePotential(q_sh.data(), q_at.data(), dp_at.data(), qp_at.data(), W.data(),
                                     dWq.data(), eps.data(), fp32);
    }
    bool devicePotentialFallback(const std::string& why) override
    {
        if (!m_device_pot) return false;
        if (noteSolvationStepFailure(why)) return true;   // WP6: only the solvation goes back to the host
        if (m_pot_fallback.empty()) {
            const std::string err = m_ctx ? m_ctx->devicePotentialError() : std::string();
            m_pot_fallback = why + (err.empty() ? std::string() : " (" + err + ")");
        }
        return true;
    }
    std::string lastError() const override
    {
        if (!m_res_error.empty()) return m_res_error;   // the last resident-loop call failed (WP2)
        return m_ctx ? m_ctx->devicePotentialError() : std::string();
    }

    // ---- BACKEND-SPECIFIC 5: fully device-resident GFN2 loop, closed shell (WP2) --------------
    // Claude Generated (Oct 2026), -xtb.rocm_resident_loop. The host polls residentScfStep (one fused
    // iteration on the device, XtbHipContext::residentScfStep) instead of driving the loop itself. Needs the
    // WP1 device potential of the same geometry. Open shell: setResidentOpenShell (WP4) below; the Broyden
    // seed preconditioner: residentSetPreconditioner (WP5) below; solvent: beginSolvation (WP6) below. The device
    // gradient potentials (prepareConvergedPotentials) keep their host path. A failed setup or step turns the loop
    // off for the rest of the run with one notice (residentLoopFallbackReason()).
    void setResidentLoop(bool on) { m_res_loop = on; }
    bool supportsResidentLoop() const override
    {
        return m_res_loop && supportsDevicePotential() && m_res_fallback.empty();
    }
    bool beginDispersionWeights(const std::vector<double>& cn, const std::vector<double>& gi,
                                const std::vector<double>& zeff, const std::vector<double>& refcn,
                                const std::vector<double>& refcovcn, const std::vector<double>& refq,
                                const std::vector<int>& nref) override
    {
        m_res_error.clear();
        if (!m_ctx || !supportsResidentLoop()) return false;
        const size_t n = cn.size(), w = 7 * n;
        if (n == 0 || static_cast<int>(n) != m_nat || gi.size() != n || zeff.size() != n || nref.size() != n
            || refcn.size() != w || refcovcn.size() != w || refq.size() != w)
            return noteResidentBegin(false, "D4 reference tables do not match the atom count");
        return noteResidentBegin(m_ctx->beginDispersionWeights(m_nat, cn.data(), gi.data(), zeff.data(),
                                                               refcn.data(), refcovcn.data(), refq.data(),
                                                               nref.data()), nullptr);
    }
    bool beginResidentLoop(const Vector& q_sh0, const Eigen::MatrixXd& dp_at0, const Eigen::MatrixXd& qp_at0,
                           const Vector& q_at0, const Vector& n0_sh, const Vector& n0_at, double Tele,
                           double n_elec, int nocc_pairs, double alpha, int max_hist, double w0) override
    {
        m_res_error.clear();
        if (!m_ctx || !supportsResidentLoop()) return false;
        if (q_sh0.size() != m_nsh || n0_sh.size() != m_nsh || q_at0.size() != m_nat || n0_at.size() != m_nat
            || dp_at0.rows() != 3 || dp_at0.cols() != m_nat || qp_at0.rows() != 6 || qp_at0.cols() != m_nat)
            return noteResidentBegin(false, "initial guess does not match the basis");
        return noteResidentBegin(m_ctx->beginResidentLoop(m_nsh, m_nat, m_n, Tele, n_elec, nocc_pairs,
                                                          q_sh0.data(), dp_at0.data(), qp_at0.data(),
                                                          q_at0.data(), n0_sh.data(), n0_at.data(), alpha,
                                                          max_hist, w0), nullptr);
    }
    // ---- BACKEND-SPECIFIC 6: fully device-resident GFN1 loop, closed shell (WP3) --------------
    // Claude Generated (Oct 2026), -xtb.rocm_resident_loop (and the generic -scf_gpu_gfn1_resident, read by the
    // host). SCC vector = shell charges, v_sh = gamma q_sh, v_at = Gamma_A q_A^2; does not need the GFN2 device
    // potential. residentScfStep, residentLoopCharges and the fallback below serve both loops.
    void setResidentLoopGfn1(bool on) { m_res_loop_gfn1 = on; }
    bool supportsResidentLoopGfn1() const override { return m_res_loop_gfn1 && m_res_fallback.empty(); }
    bool beginResidentLoopGfn1(const Vector& q_sh0, const Vector& q_at0, const Vector& n0_sh, const Vector& n0_at,
                               const std::vector<double>& gamma3_at, double Tele, double n_elec, int nocc_pairs,
                               double alpha, int max_hist, double w0) override
    {
        m_res_error.clear();
        if (!m_ctx || !supportsResidentLoopGfn1()) return false;
        if (q_sh0.size() != m_nsh || n0_sh.size() != m_nsh || q_at0.size() != m_nat || n0_at.size() != m_nat
            || static_cast<int>(gamma3_at.size()) != m_nat)
            return noteResidentBegin(false, "initial guess does not match the basis");
        return noteResidentBegin(m_ctx->beginResidentLoopGfn1(m_nsh, m_nat, m_n, Tele, n_elec, nocc_pairs,
                                                              q_sh0.data(), q_at0.data(), n0_sh.data(),
                                                              n0_at.data(), gamma3_at.data(), alpha, max_hist,
                                                              w0), nullptr);
    }
    // ---- BACKEND-SPECIFIC 7: open shell on the resident loop (WP4) -------------------------------
    // Claude Generated (Oct 2026), -xtb.rocm_resident_open_shell. The host calls this after a successful
    // beginResidentLoop / beginResidentLoopGfn1 for n_alpha != n_beta; true keeps the resident loop, whose
    // occupation step then fills the same eigenvectors per spin channel (k_res_occupations_uhf). false (switch
    // off) leaves the host-driven loop, as before WP4.
    void setResidentOpenShellAllowed(bool on) { m_res_open = on; }
    bool setResidentOpenShell(double n_alpha, double n_beta) override
    {
        if (!m_ctx || !m_res_open) return false;
        if (m_ctx->setResidentOpenShell(n_alpha, n_beta)) return true;
        m_res_error = m_ctx->residentLoopError();
        return false;
    }
    bool residentScfStep(bool fp32, double& dq, double& e_band, double& e_coulomb, double& e_third,
                         double& e_multipole) override
    {
        m_res_error.clear();
        if (!m_ctx) return false;
        const bool ok = m_ctx->residentScfStep(fp32, &dq, &e_band, &e_coulomb, &e_third, &e_multipole);
        if (!ok) m_res_error = m_ctx->residentLoopError();
        return ok;
    }
    bool residentLoopCharges(Vector& q_sh, Vector& q_at, Eigen::MatrixXd& dp_at, Eigen::MatrixXd& qp_at,
                             Vector& eps) override
    {
        if (!m_ctx || q_sh.size() != m_nsh || q_at.size() != m_nat || dp_at.rows() != 3 || dp_at.cols() != m_nat
            || qp_at.rows() != 6 || qp_at.cols() != m_nat || eps.size() != m_n) return false;
        return m_ctx->residentLoopCharges(q_sh.data(), q_at.data(), dp_at.data(), qp_at.data(), eps.data());
    }
    bool residentLoopFallback(const std::string& why) override
    {
        if (!m_res_loop && !m_res_loop_gfn1) return false;
        if (noteSolvationStepFailure(why)) { m_res_error.clear(); return true; }   // WP6, see devicePotentialFallback
        if (m_res_fallback.empty()) m_res_fallback = why;
        m_res_error.clear();
        return true;
    }
    // ---- BACKEND-SPECIFIC 8: Broyden seed preconditioner on the resident loop (WP5) -----------------
    // Claude Generated (Oct 2026), -scf_precond frag with -xtb.rocm_resident_precond (default true). The host builds
    // M (xtb_precond.cpp), flattens it and remaps the indices to the CUDA device packing (moments atom-major,
    // xtb_native.cpp); the ROCm loop keeps the host packSCC order, so the indices are mapped back here before the
    // upload (XtbHipContext::broydenSetPreconditioner). GFN1 (SCC vector = q_sh) needs no remap. No multi-GPU
    // split on ROCm (base multiGpuScfSplit false).
    void setResidentPreconditioner(bool on) { m_res_pc = on; }
    bool supportsResidentPreconditioner() const override { return m_res_pc; }
    bool supportsResidentPreconditionerGfn1() const override { return m_res_pc; }
    unsigned long long residentPreconditionerTag() const override
    {
        return m_ctx ? m_ctx->broydenPreconditionerTag() : 0;
    }
    bool residentSetPreconditioner(const curcuma::xtb::GpuPrecondFlat* data, unsigned long long tag,
                                   curcuma::xtb::GpuBroydenGuardFn guard, bool keep, bool best) override
    {
        m_res_error.clear();
        if (!m_ctx || !m_res_pc) return false;
        auto done = [this](bool ok) {
            if (!ok) m_res_error = m_ctx->residentLoopError();
            return ok;
        };
        if (!data) return done(m_ctx->broydenSetPreconditioner(nullptr, tag, std::move(guard), keep, best));
        const int N = data->N, nsh = m_nsh, nat = m_nat;
        const bool gfn2 = (N == nsh + 9 * nat);
        if ((N != nsh && !gfn2) || static_cast<int>(data->fallback.size()) != N) {
            m_res_error = "preconditioner size does not match the SCC vector";
            return false;
        }
        // Inverse of the host remap: device dp at nsh + 3*atom + k -> host nsh + k*nat + atom, device qp at
        // nsh + 3*nat + 6*atom + k -> host nsh + 3*nat + k*nat + atom.
        auto to_host = [nsh, nat](int j) {
            if (j < nsh) return j;
            if (j < nsh + 3 * nat) { const int r = j - nsh; return nsh + (r % 3) * nat + r / 3; }
            const int r = j - nsh - 3 * nat;
            return nsh + 3 * nat + (r % 6) * nat + r / 6;
        };
        std::vector<int> idx(data->idx);
        std::vector<double> fb(data->fallback);
        if (gfn2) {
            for (int& j : idx) j = to_host(j);
            for (int j = 0; j < N; ++j) fb[to_host(j)] = data->fallback[j];
        }
        XtbHipContext::BroydenPrecond pc;
        pc.N = N;
        pc.nblocks = data->nblocks();
        pc.max_dim = data->max_dim;
        pc.idx_off = data->idx_off.data();
        pc.idx = idx.data();
        pc.b_off = data->b_off.data();
        pc.B = data->B.data();
        pc.fallback = fb.data();
        return done(m_ctx->broydenSetPreconditioner(&pc, tag, std::move(guard), keep, best));
    }
    // ---- BACKEND-SPECIFIC 9: implicit solvation in the device potential (WP6) ----------------------
    // Claude Generated (Oct 2026), -xtb.rocm_device_solvation. Port of the CUDA supportsDeviceSolvation /
    // beginSolvation (xtb_gpu_method.cpp): the Born matrix B of ALPB/GBSA goes up once per geometry and the device
    // potential adds v_at += B q_at (GFN2: WP1 step and WP2 loop; GFN1: WP3 loop, plus the CM5 shift B cm5). The
    // solvation energy and its nuclear gradient stay on the host (post-SCF, from the converged charges), as for
    // CUDA. A failed upload or step returns the solvation to the host for the rest of the run (one notice); the
    // device potential and the resident loop stay available for gas-phase calculations.
    void setDeviceSolvation(bool on) { m_device_solv = on; }
    bool supportsDeviceSolvation() const override
    {
        return m_device_solv && m_solv_fallback.empty() && supportsDevicePotential();
    }
    bool supportsDeviceSolvationGfn1() const override
    {
        return m_device_solv && m_solv_fallback.empty() && supportsResidentLoopGfn1();
    }
    bool beginSolvation(int nat, const double* born_mat) override
    {
        if (!m_ctx || !supportsDeviceSolvation()) return false;
        return noteSolvationBegin(m_ctx->beginSolvation(nat, born_mat, nullptr, solvationMaxBytes()));
    }
    bool beginSolvationGfn1(int nat, const double* born_mat, const double* v_shift) override
    {
        if (!m_ctx || !supportsDeviceSolvationGfn1()) return false;
        return noteSolvationBegin(m_ctx->beginSolvation(nat, born_mat, v_shift, solvationMaxBytes()));
    }
    /**
     * @brief Called after each calculation: a failed solvation upload of this calculation becomes a permanent
     *        fallback (later solvated calculations keep the solvation on the host) with one notice.
     * @return the first device-solvation fallback of this backend ("" = none)
     */
    const std::string& solvationFallbackReason()
    {
        if (m_solv_fallback.empty() && !m_solv_begin_failed.empty())
            m_solv_fallback = "setup (" + m_solv_begin_failed + ")";
        m_solv_begin_failed.clear();
        return m_solv_fallback;
    }

    /**
     * @brief Called after each calculation: a setup whose begin failed (the host then took the
     *        host-driven loop) becomes a permanent fallback with one notice.
     * @return the first resident-loop fallback of this backend ("" = none)
     */
    const std::string& residentLoopFallbackReason()
    {
        if (m_res_fallback.empty() && !m_res_begin_failed.empty())
            m_res_fallback = "setup (" + m_res_begin_failed + ")";
        m_res_begin_failed.clear();
        return m_res_fallback;
    }
    /**
     * @brief Called after each calculation: a setup whose last begin attempt failed (the host then
     *        built the potential itself) turns into a permanent fallback with one notice.
     * @return the first device-potential fallback of this backend ("" = none)
     */
    const std::string& devicePotentialFallbackReason()
    {
        if (m_pot_fallback.empty() && !m_pot_begin_failed.empty())
            m_pot_fallback = "setup (" + m_pot_begin_failed + ")";
        m_pot_begin_failed.clear();
        return m_pot_fallback;
    }

    // ---- BACKEND-SPECIFIC 10: post-SCF dispersion and CN on the device (WP8a) ----------------------------
    // Claude Generated (Oct 2026). Port of the CUDA overrides (xtb_gpu_method.cpp), same contracts; the host
    // (xtb_native.cpp calcDispersionEnergy / primeD4CnFromDevice / the integral block) runs its own code when
    // these return false. Switches -xtb.rocm_d4_postscf_device (GFN2 D4), -xtb.rocm_d3_device (GFN1 D3) and
    // -xtb.rocm_cn_device (Hamiltonian CN download and D4 covalent CN). A declined call (precondition, e.g. no
    // device D4 reference tables after a host-driven SCF) only takes the host code; a failed call turns that
    // path off for the rest of the run with one notice (postScfFallbackReason).
    enum PostScfKind { kPsD4 = 0, kPsD3 = 1, kPsCn = 2 };
    void setPostScfDevice(bool d4, bool d3, bool cn) { m_ps_on[kPsD4] = d4; m_ps_on[kPsD3] = d3; m_ps_on[kPsCn] = cn; }
    bool d4CovalentCnDevice(int nat, const int* Z, const double* xyz_bohr, const double* rad, const double* en,
                            double kn, double k4, double k5, double k6, double cutoff, double* cn_out) override
    {
        if (!postScfAllowed(kPsCn)) return false;
        return notePostScf(kPsCn, "D4 covalent CN",
                           m_ctx->d4CovalentCnDevice(nat, Z, xyz_bohr, rad, en, kn, k4, k5, k6, cutoff, cn_out));
    }
    bool downloadCn(double* cn_out, int nat) override
    {
        if (!postScfAllowed(kPsCn) || !cn_out || nat != m_nat) return false;
        return notePostScf(kPsCn, "Hamiltonian CN", m_ctx->downloadCn(cn_out));
    }
    bool d3DispersionDevice(int nat, const int* Z, const double* xyz_ang, int nt, const int* atom_type,
                            const double* c6tab, const double* r6tab, const int* nref_tab,
                            const double* refcn_tab, const double* rcov_ang, int nrcov,
                            double s6, double s8, double a1, double a2, double k1, double k2,
                            double cn_cut_ang, double bohr_per_ang, double fold_unit,
                            bool need_gradient, double* energy_out, double* grad_out) override
    {
        if (!postScfAllowed(kPsD3)) return false;
        return notePostScf(kPsD3, "D3(BJ) dispersion",
                           m_ctx->d3DispersionDevice(nat, Z, xyz_ang, nt, atom_type, c6tab, r6tab, nref_tab,
                                                     refcn_tab, rcov_ang, nrcov, s6, s8, a1, a2, k1, k2,
                                                     cn_cut_ang, bohr_per_ang, fold_unit, need_gradient,
                                                     energy_out, grad_out));
    }
    bool d4PostScfDevice(int nat, const double* q_at, bool with_atm, double s9, double a1, double a2, double alp,
                         double atm_cutoff, double cn_cutoff, const double* rad, const double* en, double kn,
                         double k4, double k5, double k6, double* energy_out, double* grad_out,
                         double* dEdq_out) override
    {
        if (!postScfAllowed(kPsD4)) return false;
        const bool ok = notePostScf(kPsD4, "post-SCF D4 (weights, 2-body, ATM, CN fold)",
                                    m_ctx->d4PostScfDevice(nat, q_at, with_atm, s9, a1, a2, alp, atm_cutoff,
                                                           cn_cutoff, rad, en, kn, k4, k5, k6, energy_out,
                                                           grad_out, dEdq_out));
        if (ok && CurcumaLogger::get_verbosity() >= 3 && m_ctx->postScfAtmMs() >= 0.0)
            CurcumaLogger::info(fmt::format("  D4 ATM kernel (ROCm, FP64): {:.2f} ms", m_ctx->postScfAtmMs()));
        return ok;
    }
    /// First failure of a post-SCF device path ("" = none); the caller prints it once.
    const std::string& postScfFallbackReason(int kind) const { return m_ps_fallback[kind]; }

    // ---- BACKEND-SPECIFIC 11: gradient path and lazy host copies (WP8b) ------------------------------------
    // Claude Generated (Oct 2026). Port of the CUDA overrides prepareConvergedPotentials / gradientResident (switch
    // -xtb.rocm_grad_resident), gradientTailDevice / repulsionEnergyDevice (-xtb.rocm_grad_tail_device), plus the lazy
    // host copies (-xtb.rocm_defer_host_integrals): deferHostCopies (the host skips the multipole-integral download
    // and the post-SCF potential rebuild) and downloadOverlapInto / downloadH0Into (ensureHostIntegrals downloads
    // the device S/H0 instead of rebuilding them on the host). A false return runs the existing host code; a device
    // failure turns that path off for the rest of the run with one notice (gradPathFallbackReason).
    enum GradKind { kGrRes = 0, kGrTail = 1, kGrDefer = 2 };
    void setGradPath(bool resident, bool tail, bool defer)
    {
        m_gr_on[kGrRes] = resident; m_gr_on[kGrTail] = tail; m_gr_on[kGrDefer] = defer;
    }
    bool prepareConvergedPotentials(bool gfn2) override
    {
        if (!gradAllowed(kGrRes)) return false;
        return noteGrad(kGrRes, "converged gradient potentials", m_ctx->convergedPotentials(gfn2));
    }
    bool gradientResident(const Vector& eps, int nocc_orbs, Matrix& grad_out, Vector& dEdcn_out) override
    {
        if (!gradAllowed(kGrRes) || m_nat <= 0 || static_cast<int>(eps.size()) < nocc_orbs) return false;
        std::vector<double> grad(3 * static_cast<size_t>(m_nat), 0.0), dEdcn(m_nat, 0.0);
        if (!noteGrad(kGrRes, "gradient with the device potentials",
                      m_ctx->gradientResident(eps.data(), nocc_orbs, grad.data(), dEdcn.data())))
            return false;
        grad_out.resize(m_nat, 3);
        for (int i = 0; i < m_nat; ++i)
            for (int k = 0; k < 3; ++k) grad_out(i, k) = grad[3 * i + k];
        dEdcn_out.resize(m_nat);
        for (int i = 0; i < m_nat; ++i) dEdcn_out(i) = dEdcn[i];
        return true;
    }
    bool gradientTailDevice(int nat, bool is_gfn2, bool with_mp, const double* xyz_bohr, const double* dEdcn,
                            const double* rcov_au, const double* q_at, const double* dp_at, const double* qp_at,
                            const double* mrad, const double* cn, const double* rad, const double* vcn,
                            double dmp3, double dmp5, double shift, double kexp, double rmax,
                            double* grad_out) override
    {
        if (!gradAllowed(kGrTail) || nat != m_nat) return false;
        return noteGrad(kGrTail, "gradient tail (multipole interaction, CN chain rule)",
                        m_ctx->gradientTailDevice(nat, is_gfn2, with_mp, xyz_bohr, dEdcn, rcov_au, q_at, dp_at, qp_at,
                                                  mrad, cn, rad, vcn, dmp3, dmp5, shift, kexp, rmax, grad_out));
    }
    bool repulsionEnergyDevice(int nat, bool is_gfn2, const int* z, const double* xyz_bohr, const double* alpha,
                               const double* zeff, double kexp, double rexp, double kexp_light,
                               double* energy_out) override
    {
        if (!gradAllowed(kGrTail) || nat != m_nat) return false;
        return noteGrad(kGrTail, "repulsion energy",
                        m_ctx->repulsionEnergyDevice(nat, is_gfn2, z, xyz_bohr, alpha, zeff, kexp, rexp, kexp_light,
                                                     energy_out));
    }
    bool deferHostCopies() const override { return gradAllowed(kGrDefer); }
    bool downloadOverlapInto(double* S_colmajor) override
    {
        if (!gradAllowed(kGrDefer) || m_n <= 0) return false;
        return noteGrad(kGrDefer, "host S/H0 copies", deferFault() ? false
                        : m_ctx->downloadOverlap(S_colmajor));
    }
    bool downloadH0Into(double* H0_colmajor) override
    {
        if (!gradAllowed(kGrDefer) || m_n <= 0) return false;
        return noteGrad(kGrDefer, "host S/H0 copies", m_ctx->downloadH0(H0_colmajor));
    }
    /// First failure of a WP8b path ("" = none); the caller prints it once.
    const std::string& gradPathFallbackReason(int kind) const { return m_gr_fallback[kind]; }

private:
    bool gradAllowed(int kind) const { return m_ctx && m_gr_on[kind] && m_gr_fallback[kind].empty(); }
    /// CURCUMA_ROCM_FAIL_GRAD=defer: the lazy download of the host S/H0 copies fails (the host then builds them).
    static bool deferFault()
    {
        static const bool f = [] { const char* e = std::getenv("CURCUMA_ROCM_FAIL_GRAD"); return e && std::string(e) == "defer"; }();
        return f;
    }
    bool noteGrad(int kind, const char* what, bool ok)
    {
        if (ok) {
            if (CurcumaLogger::get_verbosity() >= 2)
                CurcumaLogger::info(kind == kGrDefer ? fmt::format("ROCm: {} from the device integrals", what)
                                                     : fmt::format("ROCm: {} on the device", what));
            return true;
        }
        const std::string err = kind != kGrDefer ? m_ctx->gradientPathError()
                              : (deferFault() ? "forced by CURCUMA_ROCM_FAIL_GRAD=defer" : "device download failed");
        static const char* const label[3] = { "gradient potentials (-xtb.rocm_grad_resident)",
                                              "gradient tail and repulsion (-xtb.rocm_grad_tail_device)",
                                              "lazy host S/H0 copies (-xtb.rocm_defer_host_integrals)" };
        if (m_gr_fallback[kind].empty())
            m_gr_fallback[kind] = std::string(label[kind]) + " (failed at the " + what
                                  + (err.empty() ? std::string() : ": " + err) + ")";
        return false;
    }
    bool        m_gr_on[3] = {false, false, false};   ///< WP8b switches (resident, tail, defer)
    std::string m_gr_fallback[3];                     ///< first failure per WP8b path ("" = none)

    bool postScfAllowed(int kind) const { return m_ctx && m_ps_on[kind] && m_ps_fallback[kind].empty(); }
    bool notePostScf(int kind, const char* what, bool ok)
    {
        if (ok) {
            if (CurcumaLogger::get_verbosity() >= 2)
                CurcumaLogger::info(fmt::format("ROCm: {} on the device", what));
            return true;
        }
        const std::string err = m_ctx->postScfError();
        if (m_ctx->postScfUnavailable()) {
            if (CurcumaLogger::get_verbosity() >= 2)
                CurcumaLogger::info(fmt::format("ROCm: {} on the host ({})", what, err));
            return false;
        }
        static const char* const label[3] = { "post-SCF D4", "D3(BJ) dispersion", "Hamiltonian and D4 covalent CN" };
        if (m_ps_fallback[kind].empty())
            m_ps_fallback[kind] = std::string(label[kind]) + " (failed at the " + what
                                  + (err.empty() ? std::string() : ": " + err) + ")";
        return false;
    }
    bool        m_ps_on[3] = {false, false, false};   ///< WP8a switches (d4, d3, cn), set by the constructor
    std::string m_ps_fallback[3];                     ///< first failure per WP8a path ("" = none)

    bool notePotentialBegin(bool ok)
    {
        if (ok) m_pot_begin_failed.clear();
        else {
            m_pot_begin_failed = m_ctx->devicePotentialError();
            if (m_pot_begin_failed.empty()) m_pot_begin_failed = "unknown error";
        }
        return ok;
    }
    bool noteResidentBegin(bool ok, const char* why)
    {
        if (ok) { m_res_begin_failed.clear(); return true; }
        m_res_error = why ? std::string(why) : (m_ctx ? m_ctx->residentLoopError() : std::string());
        m_res_begin_failed = m_res_error.empty() ? std::string("unknown error") : m_res_error;
        return false;
    }
    /// WP6: limit for the device Born matrix (nat^2 doubles, 2 GB = about 15800 atoms); the developer variable
    /// CURCUMA_ROCM_SOLV_MAX_MB lowers it to exercise the size fallback.
    static double solvationMaxBytes()
    {
        static const double lim = [] {
            const char* e = std::getenv("CURCUMA_ROCM_SOLV_MAX_MB");
            const double mb = e ? std::atof(e) : 0.0;
            return mb > 0.0 ? mb * 1.0e6 : 2.0e9;
        }();
        return lim;
    }
    bool noteSolvationBegin(bool ok)
    {
        if (ok) { m_solv_begin_failed.clear(); return true; }
        m_solv_begin_failed = m_ctx ? m_ctx->solvationError() : std::string();
        if (m_solv_begin_failed.empty()) m_solv_begin_failed = "unknown error";
        return false;
    }
    /// A step failed: if the solvation part failed, record it as the solvation fallback and return true (the host
    /// then rebuilds that SCF's potential on the host, incl. the solvation, and later solvated calculations keep
    /// the host solvation because supportsDeviceSolvation() turns false).
    bool noteSolvationStepFailure(const std::string& why)
    {
        if (!m_ctx || !m_ctx->solvationFailed()) return false;
        if (m_solv_fallback.empty()) {
            const std::string err = m_ctx->solvationError();
            m_solv_fallback = why + (err.empty() ? std::string() : " (" + err + ")");
        }
        return true;
    }
    bool        m_device_solv = false;     ///< WP6: -xtb.rocm_device_solvation
    std::string m_solv_begin_failed;       ///< error of the last failed solvation upload of the current calculation
    std::string m_solv_fallback;           ///< first permanent device-solvation fallback ("" = none)

    bool        m_res_loop = false;
    bool        m_res_loop_gfn1 = false;   ///< WP3: GFN1 resident loop allowed
    bool        m_res_open = false;        ///< WP4: open shell on the resident loop allowed
    bool        m_res_pc = false;          ///< WP5: Broyden seed preconditioner on the resident loop allowed
    std::string m_res_error;         ///< reason of the last failed resident-loop call (lastError)
    std::string m_res_begin_failed;   ///< error of the last failed begin of the current calculation
    std::string m_res_fallback;       ///< first permanent resident-loop fallback ("" = none)

    bool             m_device_pot = false;
    bool             m_otf_mp = false;     ///< WP11: -xtb.rocm_otf_multipole
    std::vector<int> m_sh2at;             ///< shell -> atom (host copy)
    std::string      m_pot_begin_failed;   ///< error of the last failed begin of the current calculation
    std::string      m_pot_fallback;       ///< first permanent device-potential fallback ("" = none)

    bool        m_device_mp = false;
    std::string m_mp_error;      ///< reason of the last failed upload
    std::string m_mp_fallback;   ///< first fallback, read by reportSolverFallbacks
};
#endif // HAVE_ROCSOLVER
}  // namespace

XtbHipComputationalMethod::XtbHipComputationalMethod(MethodType method, const json& config)
    : XtbGpuAdapter(method, config, "ROCm", "ROCm/HIP device")
{
    if (!gpuActive()) return;   // the adapter already warned; the CPU path stands

#ifdef HAVE_ROCSOLVER
    curcuma::xtb::XTB* xtb = cpuSolver();
    if (!xtb) return;
    XtbHipContext* ctx = context();

    // Stage 1: install the GPU eigensolver. solveEigen() delegates the per-iteration
    // generalized eigenproblem (F, S = L·Lᵀ) → (C, eps) to rocSOLVER (dsygvd) on the
    // device; everything else (integrals, Fock, density, gradient) stays on the CPU.
    // ROCm-SPECIFIC: rocSOLVER solves the GENERALIZED problem directly, so the host
    // reconstructs S from the Cholesky factor and hands both matrices over — no host
    // reduction/back-transform (that is the Vulkan route). Works for GFN1 and GFN2 alike
    // (the host hands it the complete Fock).
    xtb->setExternalEigensolver(
        [ctx](const Matrix& F, const Eigen::MatrixXd& L,
              Matrix& C, Vector& eps) -> bool {
            const int n = static_cast<int>(F.rows());
            if (n <= 0 || L.rows() != n || L.cols() != n) return false;
            Eigen::MatrixXd Fcm = F;                                   // column-major
            Eigen::MatrixXd Ll  = L.triangularView<Eigen::Lower>();    // lower Cholesky
            Eigen::MatrixXd Scm = Ll * Ll.transpose();                 // S = L·Lᵀ
            eps.resize(n);
            Eigen::MatrixXd Ccm(n, n);
            if (!ctx->solveGeneralized(Fcm.data(), Scm.data(), n, eps.data(), Ccm.data()))
                return false;
            C = Ccm;
            return true;
        });

    // Stage 2: install the device-resident SCF backend. Under the default Broyden mixing,
    // XTB::Calculation keeps H0/S (and the per-iteration density and MO coefficients)
    // resident on the device — the Fock build + populations are HIP kernels, the density a
    // rocBLAS GEMM, the eigensolve rocSOLVER. GFN1 runs the isotropic loop, GFN2 the
    // multipole loop (R-AP2).
    bool device_multipole = false;   // -xtb.rocm_device_multipole (Claude Generated, Oct 2026)
    bool device_potential = false;   // -xtb.rocm_device_potential (WP1, Claude Generated, Oct 2026)
    bool resident_loop = false;      // -xtb.rocm_resident_loop (WP2, Claude Generated, Oct 2026)
    bool resident_loop_gfn1 = false; // the same switch for GFN1 (WP3, Claude Generated, Oct 2026)
    bool resident_open = false;      // -xtb.rocm_resident_open_shell (WP4, Claude Generated, Oct 2026)
    bool resident_precond = false;   // -xtb.rocm_resident_precond (WP5, Claude Generated, Oct 2026)
    bool device_solvation = false;   // -xtb.rocm_device_solvation (WP6, Claude Generated, Oct 2026)
    bool ps_d4 = false, ps_d3 = false, ps_cn = false;   // WP8a switches (Claude Generated, Oct 2026)
    bool gr_res = false, gr_tail = false, gr_defer = false;   // WP8b switches (Claude Generated, Oct 2026)
    bool otf_mp = false;             // -xtb.rocm_otf_multipole (WP11, Claude Generated, Oct 2026)
    // Resident-solver switches -xtb.rocm_cached_cholesky (default false),
    // -xtb.rocm_density_syrk (default false) and -xtb.rocm_explicit_linv (default false); the fallbacks here are the defaults of the
    // PARAM block in native_xtb_method.h. Claude Generated (Oct 2026).
    {
        auto cfgBool = [&config](const char* key, bool def) -> bool {
            const json* x = nullptr;
            if (config.contains("xtb") && config["xtb"].is_object() && config["xtb"].contains(key))
                x = &config["xtb"][key];
            else if (config.contains(key))
                x = &config[key];
            if (!x) return def;
            if (x->is_boolean()) return x->get<bool>();
            if (x->is_number()) return x->get<double>() != 0.0;
            if (x->is_string()) {
                const std::string v = x->get<std::string>();
                if (v == "true" || v == "1" || v == "on") return true;
                if (v == "false" || v == "0" || v == "off") return false;
            }
            return def;
        };
        const bool cached_chol  = cfgBool("rocm_cached_cholesky", false);
        const bool density_syrk = cfgBool("rocm_density_syrk", false);
        const bool explicit_linv = cfgBool("rocm_explicit_linv", false);    // Claude Generated (Oct 2026), default false since Oct 17
        ctx->setSolverOptions(cached_chol, density_syrk, explicit_linv);
        // WP10 (Claude Generated, Oct 2026): -xtb.rocm_fixed_point, order-independent gradient accumulation.
        ctx->setFixedPoint(cfgBool("rocm_fixed_point", kRocmFixedPointDefault));
        // Claude Generated (Oct 2026): -xtb.rocm_device_multipole, default false (PARAM default
        // in native_xtb_method.h; this fallback is what takes effect without the key).
        device_multipole = cfgBool("rocm_device_multipole", false);
        // Claude Generated (Oct 2026): -xtb.rocm_device_potential (WP1), default true; this fallback
        // is the default that takes effect (PARAM in native_xtb_method.h).
        device_potential = cfgBool("rocm_device_potential", true);
        // Claude Generated (Oct 2026): -xtb.rocm_resident_loop (WP2), default true; this fallback is the
        // default that takes effect (PARAM in native_xtb_method.h). The device Broyden restarts its
        // history every SCF, so a CURCUMA_SCF_PRECOND_FILE and -keep_diis true keep the host-driven loop;
        // -scf_precond frag runs on the device since WP5 (rocm_resident_precond below).
        // WP3 (Claude Generated, Oct 2026): GFN1 uses the same switch and exclusions but not the GFN2 device
        // potential (its loop builds the potential itself).
        resident_loop = cfgBool("rocm_resident_loop", true);
        // Claude Generated (Oct 2026): -xtb.rocm_resident_open_shell (WP4), default true; this fallback is the
        // default that takes effect (PARAM in native_xtb_method.h). Only acts with the resident loop on.
        resident_open = cfgBool("rocm_resident_open_shell", true);
        // Claude Generated (Oct 2026): -xtb.rocm_resident_precond (WP5), default true; this fallback is the
        // default that takes effect (PARAM in native_xtb_method.h). true: -scf_precond frag runs on the resident
        // loop (device Broyden with M); false: -scf_precond other than none keeps the host-driven loop (before WP5).
        resident_precond = cfgBool("rocm_resident_precond", true);
        // Claude Generated (Oct 2026): -xtb.rocm_device_solvation (WP6); this fallback is the default that takes
        // effect (PARAM in native_xtb_method.h). Acts with the GFN2 device potential and the GFN1 resident loop.
        device_solvation = cfgBool("rocm_device_solvation", true);
        // Claude Generated (Oct 2026, WP8a): -xtb.rocm_d4_postscf_device, -xtb.rocm_d3_device, -xtb.rocm_cn_device;
        // these fallbacks are the defaults that take effect (PARAM in native_xtb_method.h).
        ps_d4 = cfgBool("rocm_d4_postscf_device", true);
        ps_d3 = cfgBool("rocm_d3_device", true);
        ps_cn = cfgBool("rocm_cn_device", true);
        // Claude Generated (Oct 2026, WP8b): -xtb.rocm_grad_resident, -xtb.rocm_grad_tail_device,
        // -xtb.rocm_defer_host_integrals; these fallbacks are the defaults that take effect (PARAM in native_xtb_method.h).
        gr_res = cfgBool("rocm_grad_resident", kRocmGradResidentDefault);
        gr_tail = cfgBool("rocm_grad_tail_device", kRocmGradTailDefault);
        gr_defer = cfgBool("rocm_defer_host_integrals", kRocmDeferHostDefault);
        // WP11 (Claude Generated, Oct 2026): -xtb.rocm_otf_multipole; this fallback is the default that takes effect
        // (PARAM in native_xtb_method.h). Acts with the GFN2 device potential.
        otf_mp = cfgBool("rocm_otf_multipole", kRocmOtfMultipoleDefault);
        if (CurcumaLogger::get_verbosity() >= 2 && (gr_res || gr_tail || gr_defer))
            CurcumaLogger::info(fmt::format(
                "{}: ROCm gradient path: converged potentials {}, gradient tail and repulsion {}, lazy host "
                "integral copies {} (-xtb.rocm_grad_resident, -xtb.rocm_grad_tail_device, "
                "-xtb.rocm_defer_host_integrals)", getMethodName(), gr_res ? "on the device" : "on the host",
                gr_tail ? "on the device" : "on the host", gr_defer ? "on" : "off"));
        if (CurcumaLogger::get_verbosity() >= 2 && (ps_d4 || ps_d3 || ps_cn))
            CurcumaLogger::info(fmt::format(
                "{}: ROCm post-SCF on the device: D4 {}, D3 {}, CN {} (-xtb.rocm_d4_postscf_device, "
                "-xtb.rocm_d3_device, -xtb.rocm_cn_device)", getMethodName(), ps_d4 ? "on" : "off",
                ps_d3 ? "on" : "off", ps_cn ? "on" : "off"));
        // Claude Generated (Oct 2026, WP7): the generic -scf_pseudo_diag / _max_dq / _max_steps / _fp64 on the
        // ROCm resident loop (XtbHipContext::setPseudoDiagonalisation, port of the CUDA pseudo step). The PARAM
        // default true is the CUDA default; on ROCm the switch is off unless -scf_pseudo_diag is given (this
        // fallback is the ROCm default, docs/ROCM_RESIDENT_LOOP_WP.md WP7). The gates use the PARAM defaults.
        {
            auto cfgNum = [&config](const char* key, double def) -> double {
                const json* x = nullptr;
                if (config.contains("xtb") && config["xtb"].is_object() && config["xtb"].contains(key))
                    x = &config["xtb"][key];
                else if (config.contains(key))
                    x = &config[key];
                if (!x) return def;
                if (x->is_number()) return x->get<double>();
                if (x->is_string()) { try { return std::stod(x->get<std::string>()); } catch (...) {} }
                return def;
            };
            constexpr bool kRocmPseudoDiagDefault = false;
            const bool pseudo = cfgBool("scf_pseudo_diag", kRocmPseudoDiagDefault);
            const bool pseudo64 = cfgBool("scf_pseudo_diag_fp64", false);
            const double pdq = cfgNum("scf_pseudo_diag_max_dq", 0.3);
            const int pst = static_cast<int>(std::lround(cfgNum("scf_pseudo_diag_max_steps", 12)));
            ctx->setPseudoDiagonalisation(pseudo, pdq, pst, pseudo64);
            m_pseudo_diag = pseudo;
            if (CurcumaLogger::get_verbosity() >= 2 && pseudo)
                CurcumaLogger::info(fmt::format(
                    "{}: ROCm pseudo-diagonalisation on the resident loop (-scf_pseudo_diag): {} steps, previous "
                    "max|dq| < {}, at most {} in a row", getMethodName(), pseudo64 ? "FP32 and FP64" : "FP32",
                    pdq, pst));
        }
        if (resident_loop) {
            std::string precond = "none";
            const json* x = nullptr;
            if (config.contains("xtb") && config["xtb"].is_object() && config["xtb"].contains("scf_precond"))
                x = &config["xtb"]["scf_precond"];
            else if (config.contains("scf_precond"))
                x = &config["scf_precond"];
            if (x && x->is_string()) precond = x->get<std::string>();
            // The dense CURCUMA_SCF_PRECOND_FILE matrix and the history kept across SCFs (-keep_diis true) stay on
            // the host-driven loop (the device Broyden restarts its history at every begin).
            if ((!resident_precond && precond != "none" && !precond.empty())
                || std::getenv("CURCUMA_SCF_PRECOND_FILE") || cfgBool("keep_diis", false))
                resident_loop = false;
        }
        resident_loop_gfn1 = resident_loop && getMethodName() == "gfn1";
        resident_loop = resident_loop && device_potential;
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info(fmt::format(
                "{}: ROCm solver: {} eigensolve, density via {}", getMethodName(),
                !cached_chol ? "sygvd"
                             : (explicit_linv ? "cached-Cholesky with explicit L^-1 (trmm/syevd/trmm)"
                                              : "cached-Cholesky (sygst/syevd/trsm)"),
                density_syrk ? "dsyrk" : "dgemm"));
        if (CurcumaLogger::get_verbosity() >= 2 && device_multipole && getMethodName() == "gfn2")
            CurcumaLogger::info(fmt::format(
                "{}: ROCm multipole potential and energy on the device (-xtb.rocm_device_multipole)",
                getMethodName()));
        if (CurcumaLogger::get_verbosity() >= 2 && device_potential && getMethodName() == "gfn2")
            CurcumaLogger::info(fmt::format(
                "{}: ROCm SCF potential (Coulomb, third order, multipole, D4 dE/dq) built on the device "
                "(-xtb.rocm_device_potential)", getMethodName()));
        if (CurcumaLogger::get_verbosity() >= 2 && device_solvation
            && ((device_potential && getMethodName() == "gfn2") || resident_loop_gfn1))
            CurcumaLogger::info(fmt::format(
                "{}: ROCm implicit-solvation reaction field (ALPB/GBSA Born matrix) in the device SCF potential "
                "(-xtb.rocm_device_solvation)", getMethodName()));
        if (CurcumaLogger::get_verbosity() >= 2
            && ((resident_loop && getMethodName() == "gfn2") || resident_loop_gfn1))
            CurcumaLogger::info(fmt::format(
                "{}: ROCm SCF iterations fully on the device (-xtb.rocm_resident_loop), {}",
                getMethodName(), resident_open ? "closed and open shell" : "closed shell only"));
    }

    {
        auto backend = std::make_unique<HipScfBackend>(ctx);
        backend->setDeviceMultipole(device_multipole);
        backend->setDevicePotential(device_potential);
        backend->setResidentLoop(resident_loop);
        backend->setResidentLoopGfn1(resident_loop_gfn1);
        backend->setResidentOpenShellAllowed(resident_open);
        backend->setResidentPreconditioner(resident_precond);
        backend->setDeviceSolvation(device_solvation);
        backend->setPostScfDevice(ps_d4, ps_d3, ps_cn);
        backend->setGradPath(gr_res, gr_tail, gr_defer);
        backend->setOnTheFlyMultipole(otf_mp);
        m_scf_backend = std::move(backend);
    }
    xtb->setGpuScfBackend(m_scf_backend.get());

    // ROCm-SPECIFIC: FP64 is ~1/16 of FP32 on this iGPU and the per-iteration eigensolve
    // dominates the resident SCF, so default mixed precision ON for the GPU path:
    // far-from-convergence iterations solve in FP32 (rocsolver_ssygvd), reverting to FP64
    // once max|dq| < scf_fp32_threshold so the converged fixed point and energy stay FP64.
    // Matches the CUDA backend (X-AP3); Vulkan deliberately opts out.
    xtb->setMixedPrecision(true);

    if (CurcumaLogger::get_verbosity() >= 2) {
        // The format string must be a compile-time constant (fmt >= 9 checks it in a
        // consteval context), so the method-dependent part is selected as plain text.
        // Claude Generated (Oct 2026): text corrected after a code reading (docs/ROCM_HOST_PARTS_2026-10-17.md):
        // the SCF loop runs on the host with a device eigensolver; potentials, occupations and the mixer are host work.
        const char* detail = (getMethodName() == "gfn2" && resident_loop)
            ? "integral build (CN/S/H0/L/gamma + dp/qp multipole and interaction matrices), SCF potential, "
              "Fock, eigensolve, D4 reference weights, occupations, density, charges and moments, SCC "
              "energies, Broyden mixer and nuclear gradient incl. multipole-integral Pulay and D4 pair terms "
              "on the GPU"
            : (getMethodName() == "gfn2" && device_potential)
            ? "integral build (CN/S/H0/L/gamma + dp/qp multipole and interaction matrices), SCF potential "
              "(Coulomb, third order, multipole, D4 pair contraction), Fock, eigensolve, density and nuclear "
              "gradient incl. multipole-integral Pulay and D4 pair terms on the GPU; D4 reference weights, "
              "occupations, mixer, SCC energies, multipole interaction gradient and CN chain-rule on the CPU"
            : getMethodName() == "gfn2"
            ? "integral build (CN/S/H0/L/gamma + dp/qp multipole), Fock, eigensolve, density and "
              "nuclear gradient incl. multipole-integral Pulay and D4 pair terms on the GPU; "
              "SCF potentials (Coulomb, third order, multipole, D4 weights), occupations, mixer, "
              "multipole interaction gradient and CN chain-rule on the CPU"
            : resident_loop_gfn1
            ? "integral build (CN/S/H0/L/gamma), SCF potential (Coulomb, third order), Fock, eigensolve, "
              "occupations, density, charges, SCC energies, Broyden mixer and nuclear gradient on the GPU; "
              "halogen bond on the CPU"
            : "integral build (CN/S/H0/L/gamma), Fock, eigensolve, density and nuclear gradient on the GPU; "
              "SCF potentials, occupations, mixer, D3 dispersion, halogen bond and CN chain-rule on the CPU";
        // Claude Generated (Oct 2026, WP2/WP3): with -xtb.rocm_resident_loop the closed-shell iterations run on the device.
        const char* loop = ((getMethodName() == "gfn2" && resident_loop) || resident_loop_gfn1)
            ? (resident_open ? "a device-resident SCF loop for closed and open shell"
                             : "a device-resident SCF loop for closed shell (open shell: host SCF loop)")
            : "a host SCF loop";
        // WP8b (Claude Generated, Oct 2026): the gradient-path items of the resident loop follow their switches.
        std::string detail_s = detail;
        if ((getMethodName() == "gfn2" && resident_loop) || resident_loop_gfn1) {
            std::string host = gr_res ? "" : "post-SCF potentials for the gradient";
            if (!gr_tail)
                host += std::string(host.empty() ? "" : ", ")
                      + (getMethodName() == "gfn2" ? "multipole interaction gradient, CN chain-rule and repulsion energy"
                                                   : "CN chain-rule and repulsion energy");
            if (!host.empty()) detail_s += "; " + host + " on the CPU";
        }
        CurcumaLogger::info(fmt::format(
            "{}: ROCm device kernels (rocSOLVER + HIP) with {}: {}",
            getMethodName(), loop, detail_s));
    }
#else
    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::info(fmt::format(
            "{}: ROCm Stage 0 (device handshake); SCF/integrals/gradient on CPU "
            "(build without rocSOLVER)", getMethodName()));
#endif
}

// The destructor reports what an optimiser or MD driver kept silent: they run the energy
// calls at logger verbosity 0, so a fallback seen there is still pending here and is printed
// at level 1, as printGpuFallbackSummary() does for the GPU fallback summary.
// Claude Generated (Oct 2026).
XtbHipComputationalMethod::~XtbHipComputationalMethod() { reportSolverFallbacks(true); }

// Claude Generated (Oct 2026): one notice per kind of solver fallback and method object.
double XtbHipComputationalMethod::calculateEnergy(bool gradient)
{
    const double e = XtbGpuAdapter::calculateEnergy(gradient);
    reportSolverFallbacks();
    // WP7 (Claude Generated, Oct 2026): pseudo-diagonalisation counts of the last resident SCF.
    if (m_pseudo_diag && CurcumaLogger::get_verbosity() >= 2) {
        if (XtbHipContext* ctx = context()) {
            int taken = 0, rejected = 0;
            ctx->pseudoDiagonalisationCounts(taken, rejected);
            CurcumaLogger::info(fmt::format("{}: ROCm pseudo-diagonalisation (last resident SCF): {} steps taken, "
                                            "{} rejected by the gate", getMethodName(), taken, rejected));
        }
    }
    return e;
}

// Claude Generated (Oct 2026). The resident solver (xtb_hip_context.hip) keeps running on a
// slower device path when the explicit L^-1 or the cached Cholesky factor cannot be formed;
// the results stay the same up to rounding. The plugin context has no logger, so the events
// are read here. These are GPU-internal path changes, not CPU fallbacks, so they are not
// counted in the gpu_fallback registry (which -gpu_strict turns into an exit).
void XtbHipComputationalMethod::reportSolverFallbacks(bool final_report)
{
    XtbHipContext* ctx = context();
    if (!ctx) return;
#ifdef HAVE_ROCSOLVER
    // Claude Generated (Oct 2026): -xtb.rocm_device_multipole fell back to the host potential.
    // Same once-per-object rule as the solver notices below.
    // Claude Generated (Oct 2026): -xtb.rocm_device_potential fell back to the host potential build.
    if (!m_reported_pot_fallback && m_scf_backend) {
        const std::string& why = static_cast<HipScfBackend*>(m_scf_backend.get())->devicePotentialFallbackReason();
        const int v0 = CurcumaLogger::get_verbosity();
        if (!why.empty() && (v0 >= 1 || final_report)) {
            if (v0 < 1) CurcumaLogger::set_verbosity(1);
            m_reported_pot_fallback = true;
            CurcumaLogger::warn(fmt::format(
                "{}: ROCm device SCF potential failed at the {}; the SCF potential is built on the host "
                "for the rest of the run (reported once)", getMethodName(), why));
            if (v0 < 1) CurcumaLogger::set_verbosity(v0);
        }
    }
    // Claude Generated (Oct 2026): -xtb.rocm_resident_loop fell back to the host-driven loop.
    if (!m_reported_res_fallback && m_scf_backend) {
        const std::string& why = static_cast<HipScfBackend*>(m_scf_backend.get())->residentLoopFallbackReason();
        const int v0 = CurcumaLogger::get_verbosity();
        if (!why.empty() && (v0 >= 1 || final_report)) {
            if (v0 < 1) CurcumaLogger::set_verbosity(1);
            m_reported_res_fallback = true;
            CurcumaLogger::warn(fmt::format(
                "{}: ROCm device-resident SCF loop failed at the {}; the SCF runs on the host-driven loop "
                "for the rest of the run (reported once)", getMethodName(), why));
            if (v0 < 1) CurcumaLogger::set_verbosity(v0);
        }
    }
    // Claude Generated (Oct 2026, WP6): -xtb.rocm_device_solvation fell back to the host solvation.
    if (!m_reported_solv_fallback && m_scf_backend) {
        const std::string& why = static_cast<HipScfBackend*>(m_scf_backend.get())->solvationFallbackReason();
        const int v0 = CurcumaLogger::get_verbosity();
        if (!why.empty() && (v0 >= 1 || final_report)) {
            if (v0 < 1) CurcumaLogger::set_verbosity(1);
            m_reported_solv_fallback = true;
            CurcumaLogger::warn(fmt::format(
                "{}: ROCm device solvation failed at the {}; the solvation potential is built on the host "
                "for the rest of the run (reported once)", getMethodName(), why));
            if (v0 < 1) CurcumaLogger::set_verbosity(v0);
        }
    }
    // WP8a (Claude Generated, Oct 2026): a post-SCF device path (D4, D3, CN) failed and returned to the host code.
    for (int k = 0; k < 3 && m_scf_backend; ++k) {
        if (m_reported_ps_fallback[k]) continue;
        const std::string& why = static_cast<HipScfBackend*>(m_scf_backend.get())->postScfFallbackReason(k);
        const int v0 = CurcumaLogger::get_verbosity();
        if (!why.empty() && (v0 >= 1 || final_report)) {
            if (v0 < 1) CurcumaLogger::set_verbosity(1);
            m_reported_ps_fallback[k] = true;
            CurcumaLogger::warn(fmt::format(
                "{}: ROCm device {} failed; the host code computes it for the rest of the run (reported once)",
                getMethodName(), why));
            if (v0 < 1) CurcumaLogger::set_verbosity(v0);
        }
    }
    // WP8b (Claude Generated, Oct 2026): a gradient-path device step failed and returned to the host code.
    for (int k = 0; k < 3 && m_scf_backend; ++k) {
        if (m_reported_gr_fallback[k]) continue;
        const std::string& why = static_cast<HipScfBackend*>(m_scf_backend.get())->gradPathFallbackReason(k);
        const int v0 = CurcumaLogger::get_verbosity();
        if (!why.empty() && (v0 >= 1 || final_report)) {
            if (v0 < 1) CurcumaLogger::set_verbosity(1);
            m_reported_gr_fallback[k] = true;
            CurcumaLogger::warn(fmt::format(
                "{}: ROCm device {} failed; the host code runs it for the rest of the run (reported once)",
                getMethodName(), why));
            if (v0 < 1) CurcumaLogger::set_verbosity(v0);
        }
    }
    // WP7 (Claude Generated, Oct 2026): the pseudo-diagonalisation was turned off (allocation or step failure).
    if (!m_reported_pseudo_fallback && m_pseudo_diag) {
        const std::string why = ctx->pseudoDiagonalisationNotice();
        const int v0 = CurcumaLogger::get_verbosity();
        if (!why.empty() && (v0 >= 1 || final_report)) {
            if (v0 < 1) CurcumaLogger::set_verbosity(1);
            m_reported_pseudo_fallback = true;
            CurcumaLogger::warn(fmt::format(
                "{}: ROCm pseudo-diagonalisation turned off: {}; every SCF step runs the full eigensolve for the "
                "rest of the run (reported once)", getMethodName(), why));
            if (v0 < 1) CurcumaLogger::set_verbosity(v0);
        }
    }
    if (!m_reported_mp_fallback && m_scf_backend) {
        const std::string& why = static_cast<HipScfBackend*>(m_scf_backend.get())->deviceMultipoleFallback();
        const int v0 = CurcumaLogger::get_verbosity();
        if (!why.empty() && (v0 >= 1 || final_report)) {
            if (v0 < 1) CurcumaLogger::set_verbosity(1);
            m_reported_mp_fallback = true;
            CurcumaLogger::warn(fmt::format(
                "{}: ROCm device multipole potential failed at the {}; the multipole potential and "
                "energy are computed on the host for the rest of the run (reported once)",
                getMethodName(), why));
            if (v0 < 1) CurcumaLogger::set_verbosity(v0);
        }
    }
#endif
    const unsigned ev = ctx->solverFallbackEvents();
    const unsigned fresh = ev & ~m_reported_fallbacks
        & (XtbHipContext::kFallbackLinvAlloc | XtbHipContext::kFallbackLinvTrtri
           | XtbHipContext::kFallbackCholFailed | XtbHipContext::kFallbackFactorRetried);
    if (!fresh) return;
    // At verbosity 0 the warning would be swallowed: keep the bits pending for a later call;
    // the final report (destructor) raises the level to 1 for these lines.
    const int saved_verbosity = CurcumaLogger::get_verbosity();
    if (saved_verbosity < 1 && !final_report) return;
    if (saved_verbosity < 1) CurcumaLogger::set_verbosity(1);
    m_reported_fallbacks |= fresh;
    const char* src = (ev & XtbHipContext::kFallbackInjected)
        ? " [forced by CURCUMA_ROCM_FAIL_LINV]" : "";
    if (fresh & XtbHipContext::kFallbackLinvAlloc)
        CurcumaLogger::warn(fmt::format(
            "{}: ROCm solver: the explicit L^-1 buffer could not be allocated; the SCF uses the "
            "cached-Cholesky path (sygst/syevd/trsm) for the rest of the run{} (reported once)",
            getMethodName(), src));
    if (fresh & XtbHipContext::kFallbackLinvTrtri)
        CurcumaLogger::warn(fmt::format(
            "{}: ROCm solver: trtri could not invert the Cholesky factor of S; this geometry uses "
            "the cached-Cholesky path (sygst/syevd/trsm), later geometries try L^-1 again{} "
            "(reported once)", getMethodName(), src));
    if (fresh & XtbHipContext::kFallbackCholFailed)
        CurcumaLogger::warn(fmt::format(
            "{}: ROCm solver: the cached Cholesky factor of S could not be formed (allocation "
            "failed or potrf reported S not positive definite); the affected geometries use "
            "rocsolver sygvd instead{} (reported once)", getMethodName(), src));
    if (fresh & XtbHipContext::kFallbackFactorRetried)
        CurcumaLogger::warn(fmt::format(
            "{}: ROCm solver: rocSOLVER returned a Cholesky factor or inverse that failed the "
            "L L^T v = S v probe and it was recomputed{} (reported once)", getMethodName(), src));
    if (saved_verbosity < 1) CurcumaLogger::set_verbosity(saved_verbosity);
}

#endif // USE_ROCM
