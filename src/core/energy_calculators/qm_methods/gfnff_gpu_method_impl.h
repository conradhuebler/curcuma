/*
 * <GFNFFGpuMethodImpl — shared CUDA/HIP host wrapper for the GFN-FF GPU path>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): ComputationalMethod adapter for the gfnff GPU path.
 * Claude Generated (Sep 2026): the CUDA (gfnff_gpu_method) and HIP (gfnff_hip_method)
 * wrappers were the same ~1400 lines twice, differing only in the device class names,
 * a few log strings and ONE functional detail (ROCm has no device GPU-Schur EEQ, so it
 * routes many-fragment systems to the exact CPU PCG).  They are now one class template
 * parameterised on a small backend traits struct:
 *
 *   struct Backend {
 *       using Workspace = FFWorkspaceGPU;      // or FFWorkspaceHip
 *       using EEQSolver = EEQSolverGPU;        // or EEQSolverHip
 *       static constexpr const char* name           = "CUDA";           // log prose
 *       static constexpr const char* workspace_name = "FFWorkspaceGPU"; // error text
 *       static constexpr bool has_device_schur      = true;   // false on ROCm
 *       static constexpr bool has_device_ppcg       = true;   // WP7-E projected PCG (ROCm too, Oct 2026)
 *       static constexpr bool has_localinv_upload   = true;   // host-built localinv upload (ROCm too, Oct 2026)
 *       static constexpr bool has_localinv_device_build = true; // device-built localinv (ROCm too, Oct 2026)
 *       static constexpr bool has_cn_pair_skin      = true;   // K1 CN pair-list skin (ROCm too, Oct 2026)
 *       static constexpr bool has_stream_overlap    = true;   // false on ROCm (K2/K3, Oct 2026)
 *       static constexpr int  default_eeq_cpu_fragment_threshold = 0;   // 16 on ROCm
 *       static void downloadDoubles(double* host, const double* device, int n);
 *   };
 *
 * The traits are defined in gfnff_gpu_method.h / gfnff_hip_method.h, which also derive
 * the concrete GFNFFGPUComputationalMethod / GFNFFHipComputationalMethod from this
 * template and declare the explicit instantiation.  The definitions are instantiated in
 * exactly ONE translation unit per backend (gfnff_gpu_method.cpp compiled by nvcc,
 * gfnff_hip_method.cpp compiled by g++), so the generated code is the same as before.
 *
 * Architecture (unchanged):
 *   - m_gfnff          : GFNFF instance (CPU topology + EEQ charges + CN)
 *   - m_gpu_workspace  : device workspace (bonds, angles, dihedrals, inversions,
 *                        dispersion, repulsion, Coulomb, HB, XB, ATM, BATM, sTors)
 *   - every energy term runs on the device (no host residual workspace since Jun 2026)
 *
 * Per-step calculation (orchestrated here, NOT delegated to GFNFF::Calculation);
 * updated Sep 28, 2026 - the list below used to describe the March 2026 CPU-EEQ pipeline:
 *   1. host list maintenance (HB/XB, D4 pair skin, repulsion skin) + GPU computeCN()
 *   2. GPU: Phase 1 = charge-independent kernels (prepareAndLaunchChargeIndependent),
 *      D4 Gaussian weights / C6 refresh (energy mode before Phase 1, gradient mode after)
 *   3. Phase-2 EEQ on the device (nfrag == 1: dense Schur-Cholesky WP5-A; nfrag > 1: WP7-E
 *      projected PCG above the ppcg thresholds, else WP7-A); host EEQ only for skip_phase2,
 *      frozen charges, or while the fragment topology is not on the device yet
 *   4. GPU: Coulomb + post-processing + download (launchChargeDependentAndFinish)
 *
 * Usage:
 *   ./curcuma -sp mol.xyz -method gfnff -gpu cuda    # CUDA backend
 *   ./curcuma -sp mol.xyz -method gfnff -gpu rocm    # ROCm backend
 *   ./curcuma -sp mol.xyz -method gfnff              # CPU (default)
 */

#pragma once

#if defined(USE_CUDA) || defined(USE_ROCM)

#include "../computational_method.h"
#include "../ff_methods/gfnff.h"
#include "../ff_methods/eeq_solver.h"  // EEQSolveMethod enum
#include "src/core/curcuma_logger.h"
#include "src/core/gpu_fallback.h"
#include "src/core/gpu_device_pool.h"
#include "src/core/citation_registry.h"
#include "src/core/energy_calculators/ff_methods/gfnff_par.h"
#include "src/core/energy_calculators/ff_methods/cuda/gpu_utils.h"  // backend-neutral
#include "src/core/energy_calculators/ff_methods/cn_calculator.h"
#include "src/core/energy_calculators/ff_methods/forcefield.h"
#include "src/core/energy_calculators/ff_methods/ff_terms.h"

#include <atomic>
#include <chrono>
#include <type_traits>
#include <sstream>
#include <cmath>
#include <cstring>
#include <memory>

/// Claude Generated (Oct 2026): true when the backend has a static downloadEEQCharges(host, device, n)
/// (ROCm: charges of its small-system device EEQ come from a host copy made with the status read-back).
/// A backend without the member (CUDA) gets false and keeps Backend::downloadDoubles unchanged.
template <class B, class = void>
struct gfnffBackendEEQChargeDownload : std::false_type {};
template <class B>
struct gfnffBackendEEQChargeDownload<B, std::void_t<decltype(&B::downloadEEQCharges)>> : std::true_type {};

/// Claude Generated (Oct 2026): true iff the backend traits declare `has_device_repulsion_build = true` (CUDA).
template <class B, class = void>
struct gfnffBackendRepulsionBuild : std::false_type {};
template <class B>
struct gfnffBackendRepulsionBuild<B, std::void_t<decltype(B::has_device_repulsion_build)>>
    : std::bool_constant<B::has_device_repulsion_build> {};

template <class Backend>
class GFNFFGpuMethodImpl : public ComputationalMethod {
public:
    /// Backend-provided device types (FFWorkspaceGPU/EEQSolverGPU resp. ...Hip)
    using GpuWorkspace = typename Backend::Workspace;
    using GpuEEQSolver = typename Backend::EEQSolver;

    explicit GFNFFGpuMethodImpl(const std::string& method_name, const json& config);
    ~GFNFFGpuMethodImpl();

    // === ComputationalMethod interface ===

    bool setMolecule(const Mol& mol) override;
    double calculateEnergy(bool gradient = false) override;
    Matrix getGradient() const override;
    void copyGradientTo(Matrix& target) const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override;
    Position getDipole() const override;
    bool updateGeometry(const Matrix& geometry) override;
    bool hasGradient() const override { return true; }
    bool isThreadSafe() const override { return false; }
    std::string getMethodName() const override { return m_method_name; }

    // === Configuration ===

    void setThreadCount(int threads) override;
    void setParameters(const json& params) override;
    json getParameters() const override;

    // === Error handling ===

    bool hasError() const override;
    void clearError() override;
    std::string getErrorMessage() const override;

    // === Energy decomposition ===

    json getEnergyDecomposition() const override;

    /// Expose GPU dEdcn for diagnostics (valid after calculateEnergy with gradient)
    const Vector& getGPUdEdcn() const;

    /// Expose GPU workspace for gradient diagnostics (e.g. gradientBeforeCN)
    GpuWorkspace* getGPUWorkspace() const { return m_gpu_workspace.get(); }

    /// Get GPU CN result (valid after calculateEnergy)
    const Vector& getGPUCN() const { return m_gpu_cn_final; }

    // WP-P1 (May 2026): per-phase and per-stream timings for MD diagnostics
    json getLastPrepTiming() const override;
    json getStreamTimings() const override;
    void setForcePhaseTiming(bool on) override;

    /// CPU-side GFNFF instance (topology, EEQ charges, react bond list and its
    /// events). Overrides the virtual hook rather than being a plain inline
    /// accessor: this class lives in a dlopen plugin and derives from an
    /// `extern template` base, so only a vtable dispatch is resolvable from outside
    /// the plugin. Claude Generated (Sep 2026).
    GFNFF* getGFNFF() const { return m_gfnff.get(); }
    GFNFF* gfnffInstance() const override { return m_gfnff.get(); }
    bool requestsCpuFallback() const override { return m_cpu_fallback_requested; }
    std::string cpuFallbackReason() const override { return m_cpu_fallback_reason; }
    /// Device index this method is pinned to (-1 = runtime default). Claude Generated (Sep 2026).
    int gpuDevice() const { return m_gpu_device; }

private:
    /// Make the pinned device current on the calling thread (no-op when unpinned).
    void bindDevice() const { if (m_gpu_device >= 0) Backend::setDevice(m_gpu_device); }

    /**
     * @brief Initialize GPU workspace from GFNFF parameter set.
     * Called after m_gfnff->InitialiseMolecule() succeeds.
     * @return true on success; sets m_has_error on CUDA failure.
     */
    bool initGPUWorkspace();

    std::unique_ptr<GFNFF>          m_gfnff;
    // Claude Generated (Sep 2026, docs/MULTI_GPU_GAPS.md F-16): the part of the host parameter
    // set the GPU path still reads after the workspace is built (bonds + bond_hb_data for the
    // HB bond cross-reference). It used to be the FULL set, deliberately leaked as a guard against
    // a March 2026 heap corruption whose cause (async uploads from destroyed temporaries) was fixed
    // in the same commit; measured Sep 2026: 54 MB per molecule of 1410 atoms, plus ~140 MB for
    // the leaked GFNFF object - 194 MB per structure in a GPU batch. Now freed normally;
    // CURCUMA_GFNFF_GPU_LEAK=1 restores the old leak (safety valve for TECHNICAL_DEBT F-Q9).
    std::unique_ptr<GFNFFParameterSet> m_gpu_params;
    std::unique_ptr<GpuWorkspace>   m_gpu_workspace;

    // Claude Generated (March 2026): GPU EEQ solver (cuSOLVER Cholesky)
    std::unique_ptr<GpuEEQSolver>   m_eeq_gpu;

    // Pre-allocated buffers for GPU EEQ Schur complement (avoid per-step heap allocs)
    std::vector<double> m_eeq_z1;           ///< [N] A⁻¹ · b_atoms
    std::vector<double> m_eeq_Z2;           ///< [N*nfrag] A⁻¹ · C^T (column-major)
    std::vector<double> m_eeq_charges_gpu;  ///< [N] final charges from GPU path
    std::vector<double> m_schur_workspace;  ///< [nfrag*(nfrag+2)] CPU Schur: S matrix + rhs + lambda (no heap after CUDA init)

    // WP2: cached topology-constant EEQ data (avoid prepareEEQParametersForGPU per step)
    std::vector<int>    m_eeq_fraglist;          ///< [N] fragment IDs (topology-constant)
    std::vector<double> m_eeq_rhs_constraints;   ///< [nfrag] target charges (topology-constant)
    std::string         m_cpu_fallback_reason;             ///< text for the counted fallback (empty: charge-rule text)
    bool                m_cpu_fallback_requested = false;  ///< charged multi-fragment system: CPU method wanted (see requestsCpuFallback())
    int                 m_eeq_nfrag = 1;         ///< number of fragments (topology-constant)

    json             m_parameters;
    std::string      m_method_name;
    std::vector<int> m_atom_types;   ///< Element numbers (Z) — stored from setMolecule()
    bool             m_allow_unconverged_charges = false; ///< Respect allow_unconverged_charges for GPU→CPU fallback (Claude Generated Apr 2026)
    bool             m_skip_phase2 = false;               ///< Skip Phase 2 EEQ refinement, use Phase 1 topology charges (Claude Generated Apr 2026)
    bool             m_initialized   = false;
    bool             m_has_error     = false;
    int              m_gpu_device    = -1;   ///< pinned device (-1 = runtime default)
    std::string      m_error_message;
    double           m_last_energy   = 0.0;
    double           m_gpu_upload_time_ms = 0.0;
    Matrix           m_cached_gradient; ///< Cached gradient (copied from GPU workspace after calculate)

    // Claude Generated (Sep 2026): GPU-path per-phase wall-clock breakdown, so
    // getLastPrepTiming() reports real numbers instead of the CPU PrepTiming's zeros
    // (prepareCNAndEEQ is always called with skip_eeq=true on this path — see
    // calculateEnergy() below — so its own m_last_prep_timing.eeq_solve etc. stay 0).
    // These std::chrono timestamps were already taken unconditionally on every step
    // (used only for the verbosity>=2 printout below); storing them costs nothing
    // extra and the fields are populated every call — zero-cost is enforced at the
    // CALLER (simplemd.cpp only invokes LastPrepTiming() when md_diagnostics_timing
    // is set), not here.
    struct GpuPrepTiming {
        double cn = 0.0;              ///< GPU CN kernel launch (host wall-clock)
        double dlogdcn = 0.0;         ///< dlogdcn / CN-pair-list setup (gradient only)
        double hb_cn = 0.0;           ///< per-bond HB coordination number kernel
        double hbxb = 0.0;            ///< dynamic HB/XB re-detection
        // Claude Generated (Oct 2026): host-side list refreshes that were inside "total" but in no
        // named phase (step-time spikes of 35 ms / 180 ms on urea400-water1000, see Labor journal).
        double pre = 0.0;             ///< setGeometry + react topology + displacement check, before the CN launch
        double rep_update = 0.0;      ///< non-bonded repulsion pair-list refresh (step 2d)
        // Claude Generated (Oct 2026): split of hbxb (host detection / bond cross reference / device uploads)
        double hb_detect = 0.0;       ///< updateHBXBIfNeeded(): trigger test + host detection of the HB/XB lists
        double hb_rebuild = 0.0;      ///< rebuildBondHBData(): bond-HB cross reference (host)
        double hb_upload = 0.0;       ///< updateHBonds/XBonds/HBAlphaPairs/BondHBMetadata: device uploads
        double rep_upload = 0.0;      ///< part of rep_update spent in the device upload (rest: host list generation)
        double rep_pairs = 0.0;       ///< non-bonded repulsion pairs after the last refresh (count, not ms)
        double disp_update = 0.0;     ///< D4 dispersion pair-list refresh incl. device rebuild (step 2e)
        double coul_update = 0.0;     ///< Coulomb pair-list refresh (step 2e)
        double launch = 0.0;          ///< charge-independent GPU kernel launch overhead
        double eeq_solve = 0.0;       ///< EEQ Phase 2 solve — GPU Schur/PCG or CPU fallback
        double coulomb_postprocess = 0.0; ///< Coulomb kernel + postprocess + D2H download
        double total = 0.0;           ///< calculateEnergy() wall time (excl. topology/paramgen)
    } m_last_gpu_prep_timing;

    // CN chain-rule pair list: built on the device (generateCNPairListOnGPU()); this flag
    // records whether the current list is valid.
    bool                m_cn_pairs_generated = false;

    // Task #10 (Jun 2026): CN-derivative pair-list regen controls (from PARAMs)
    bool                m_cn_pair_regen = true;       ///< regen on topology displacement
    int                 m_cn_pair_regen_every = 0;    ///< force regen every N grad steps (0=off)
    double              m_cn_pair_cutoff_factor = 2.5;///< cutoff = factor*(rcov_i+rcov_j)
    long                m_cn_pair_grad_steps = 0;     ///< gradient-step counter for regen_every
    // K1 (Oct 2026, Claude Generated): Verlet skin of the CN pair list (gpu_cn_pair_skin_bohr).
    // < 0: legacy (rebuild whenever the 0.5 Bohr topology check fires, i.e. every MD step).
    double              m_cn_pair_skin = 1.0;         ///< skin in Bohr, list cutoff + 2*skin
    double              m_cn_pair_max_disp = 0.0;     ///< allowed displacement since the last build (Bohr)
    // K2/K3 (Oct 2026, Claude Generated): stream overlap switches (gpu_eeq_overlap, gpu_dispersion_overlap).
    bool                m_eeq_overlap = true;         ///< K3: device ppcg waits on the EEQ-input event, not on dc6dcn
    bool                m_disp_overlap = true;        ///< K2: gradient k_dispersion right after dc6dcn
    bool                m_eeq_priority_set = false;   ///< K2/K3: EEQ stream priority raised once

    // GPU CN result (always used — GPU CN is the default path)
    Vector              m_gpu_cn_final;        ///< Cached GPU CN result

    // RMSD-based EEQ lazy refactorization (Apr 2026)
    // Caller tracks geometry RMSD from last full Cholesky refactorization.
    // m_eeq_rmsd_threshold == 0.0 → always refactorize (default, matches reference).
    // m_eeq_rmsd_threshold >  0.0 → lazy: reuse L when per-atom RMSD < threshold (Bohr).
    double m_eeq_rmsd_threshold = 0.0;  ///< Bohr; from eeq_rmsd_threshold param
    Matrix m_eeq_ref_geom;              ///< geometry (Bohr) at last EEQ refactorization
    bool   m_eeq_has_ref_geom = false;

    // EEQ Coulomb-matrix distance cutoff (Bohr). Default 0 = no cutoff, matches Fortran
    // goed_gfnff. Set non-zero only for performance experiments — produces HF-inconsistent
    // gradients vs. the un-truncated Coulomb energy and degrades MD energy conservation.
    double m_eeq_distance_cutoff = 0.0;  ///< Bohr; from eeq_distance_cutoff param

    // WP7-B (May 2026): GPU EEQ solver strategy for nfrag>1 systems.
    // SchurCholesky → WP5-A (nfrag=1) / WP7-A (nfrag>1) — exact, full N×N Cholesky.
    // Batched      → WP7-B (nfrag>1)                   — per-fragment Cholesky, drops cross-fragment Coulomb.
    // PCG          → WP7-C (nfrag>1)                   — iterative O(k·N²), warm-started.
    // LU           → dense LU of the augmented system (CUDA only, Oct 2026).
    // LDLT         → Bunch-Kaufman LDL^T of A_nn + Schur (CUDA only, Oct 2026).
    // Auto         → EEQSolver::resolveMethod (ppcg above the size thresholds, else cholesky).
    EEQSolveMethod m_eeq_strategy = EEQSolveMethod::Auto;
    /// Claude Generated (Oct 2026): energy criterion of the batched solver on the GPU path. Rebuilds
    /// the residual r = A q_b - x of the full EEQ system on the host (one O(N^2) pass, once per
    /// topology) from the device alpha/gam/rhs and the accepted batched charges, then
    /// EEQSolver::estimateBatchedError(); warns above eeq_batched_max_error_eh.
    void estimateBatchedErrorOnHost(int N, double cutoff_sq)
    {
        const auto t0 = std::chrono::steady_clock::now();
        std::vector<double> alpha(N), gam(N), rhs(N);
        Backend::downloadDoubles(alpha.data(), m_gpu_workspace->getDeviceAlphaPtr(), N);
        Backend::downloadDoubles(gam.data(), m_gpu_workspace->getDeviceGamPtr(), N);
        Backend::downloadDoubles(rhs.data(), m_gpu_workspace->getDeviceRHSPtr(), N);
        const GeoGradMatrix& g = m_gfnff->getGeometryBohr();
        const std::vector<double>& q = m_eeq_charges_gpu;
        constexpr double tsqrt2pi = 0.797884560802866;   // as in k_eeq_build_matrix
        auto elem = [&](int i, int j) -> double {
            if (i == j) return gam[i] + tsqrt2pi / std::sqrt(alpha[i]);
            const double dx = g(i, 0) - g(j, 0), dy = g(i, 1) - g(j, 1), dz = g(i, 2) - g(j, 2);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (cutoff_sq > 0.0 && r2 > cutoff_sq) return 0.0;
            const double r = std::sqrt(r2);
            return std::erf(r / std::sqrt(alpha[i] + alpha[j])) / r;
        };
        Vector res(N);
        for (int i = 0; i < N; ++i) res(i) = elem(i, i) * q[i] - rhs[i];
        for (int i = 1; i < N; ++i)
            for (int j = 0; j < i; ++j) {
                const double a = elem(i, j);
                res(i) += a * q[j];
                res(j) += a * q[i];
            }
        std::vector<std::vector<int>> frag_atoms(m_eeq_nfrag);
        for (int i = 0; i < N && i < static_cast<int>(m_eeq_fraglist.size()); ++i) {
            const int f = m_eeq_fraglist[i] - 1;
            if (f >= 0 && f < m_eeq_nfrag) frag_atoms[f].push_back(i);
        }
        const double dE = EEQSolver::estimateBatchedError(res, frag_atoms, elem);
        const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        if (dE > m_eeq_batched_max_error)
            CurcumaLogger::warn(fmt::format(
                "EEQ batched (final charges, GPU): estimated energy error from the dropped inter-fragment coupling "
                "{:.2e} Eh > eeq_batched_max_error_eh {:.1e} Eh; cholesky or ppcg give the exact charges",
                dE, m_eeq_batched_max_error));
        else if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info(fmt::format("EEQ batched (final charges, GPU): estimated energy error {:.2e} Eh ({:.1f} ms)", dE, ms));
        if (CurcumaLogger::get_verbosity() >= 3)
            CurcumaLogger::info(fmt::format("EEQ batched error estimate (GPU path, host) took {:.1f} ms", ms));
    }

    // Oct 2026: energy criterion for the explicit batched solver (replaces eeq_batched_min_distance).
    double m_eeq_batched_max_error = 1e-6;   ///< eeq_batched_max_error_eh; 0 = no estimate
    bool   m_eeq_batched_check_pending = true; ///< estimate after the next batched solve (set per topology)


    // WP7-C (May 2026): GPU PCG parameters (mirror CPU EEQSolver defaults).
    int    m_eeq_pcg_max_iter   = 200;     ///< per-PCG-call iteration cap
    double m_eeq_pcg_tolerance  = 1e-10;   ///< convergence tolerance on |r|
    // Sep 2026 (curcuma EEQ-Löser-Benchmark): block-Jacobi preconditioner size gates — see
    // eeq_solver.h PARAM gpu_block_jacobi_max_frag_atoms / gpu_block_jacobi_max_nfrag.
    int    m_eeq_block_jacobi_max_frag_atoms = 300;
    int    m_eeq_block_jacobi_max_nfrag      = 400;

    // WP7-E (Sep 2026): GPU projected PCG — single-solve replacement for WP7-C's
    // nfrag+1-solves loop on many-fragment, large systems (curcuma EEQ-Löser-Benchmark).
    // Reuses the CPU EEQSolver's own ppcg PARAMs so CPU/GPU auto-selection agree.
    int    m_eeq_ppcg_max_iter   = 500;
    double m_eeq_ppcg_tolerance  = 1e-10;  ///< GPU tol is absolute on |r| (unlike CPU's relative eeq_ppcg_tol)
    int    m_eeq_ppcg_min_nfrag  = 1;      ///< Auto strategy: prefer ProjectedPCG over PCG when nfrag>=this
    int    m_eeq_ppcg_min_atoms  = 500;    ///< ... and N>=this
    // Claude Generated (Oct 2026): local-inverse preconditioner of WP7-E (eeq_ppcg_precond localinv).
    // CUDA and ROCm (Oct 2026): built on the device (requestLocalInverseDeviceBuild, default,
    // has_localinv_device_build) or, with CURCUMA_LOCALINV_HOST_BUILD=1 or the CUDA multi-GPU split,
    // on the host by EEQSolver::buildLocalInverse and uploaded as CSR (has_localinv_upload).
    bool   m_eeq_localinv = false;
    double m_eeq_localinv_radius = 8.0;
    double m_eeq_localinv_ext = 10.0;
    double m_eeq_localinv_rebuild = 1.0;
    int    m_eeq_localinv_threads = 1;
    int    m_host_threads = 1;   ///< -threads: host threads of the staged list upload (Oct 2026)
    std::shared_ptr<LocalInversePreconditioner> m_eeq_li;   ///< last build (host: G, H, S; device: geometry only)
    long   m_eeq_li_builds = 0, m_eeq_li_failures = 0;
    double m_eeq_li_build_ms = 0.0;
    long   m_eeq_li_seen_dev_builds = 0, m_eeq_li_seen_dev_failures = 0, m_eeq_li_seen_runtime = 0;
    static constexpr int kLocalInvMaxNfrag = 4096;   ///< host build: dense S^-1 on the device, nfrag^2 doubles (128 MB here)
    /// Builds (or keeps) the local inverse for the current geometry and uploads it; on any failure the
    /// device solver runs with the diagonal Jacobi. Same rebuild rule as the CPU solver.
    void prepareLocalInverseGPU(int N, double cutoff_sq);
    /// After a WP7-E solve: books device builds/rejections and reports every Jacobi re-solve as a
    /// GPU fallback (curcuma::reportGpuFallback, -gpu_strict). Claude Generated (Oct 2026).
    void finishLocalInverseGPU();
    /// true when the local inverse for N atoms is built on the device (backend has the device build,
    /// CURCUMA_LOCALINV_HOST_BUILD unset, and the solver can build it, e.g. no multi-GPU split).
    /// Claude Generated (Oct 2026).
    bool localInverseDeviceBuildSelected(int N) const
    {
        if constexpr (Backend::has_localinv_device_build) {
            static const bool host_build_env = [] { const char* e = std::getenv("CURCUMA_LOCALINV_HOST_BUILD"); return e && *e && *e != '0'; }();
            return !host_build_env && m_eeq_gpu && m_eeq_gpu->localInverseDeviceBuildPossible(N);
        } else {
            (void)N;
            return false;
        }
    }

    // Deliverable 3 (Jun 2026): fragment count at/above which the device EEQ solve is
    // replaced by the exact CPU PCG/block-Jacobi solver. The device path does a dense
    // N x N Cholesky for nfrag>1 (O(N^3), intractable for solvent boxes); the CPU PCG is
    // O(N^2 k) and exact. 0 = always use the device solve (CUDA default, which has the
    // device Schur paths). From eeq_rocm_cpu_fragment_threshold; the traits supply the
    // per-backend default (0 for CUDA, 16 for ROCm).
    int    m_eeq_cpu_fragment_threshold = Backend::default_eeq_cpu_fragment_threshold;
    /// Claude Generated (Oct 2026), ROCm: eeq_rocm_cpu_fragment_threshold given explicitly. Then
    /// nfrag >= threshold always takes the CPU EEQ (the pre-Oct-2026 route), also when the solver
    /// resolves to ppcg; by default the device projected PCG runs there.
    bool   m_eeq_cpu_route_forced = false;

    bool m_split_configured = false;  ///< Claude Generated (Sep 2026): multi-GPU split set up
    bool m_split_reported = false;
    int  m_calc_count = 0;  ///< counts calculateEnergy() calls; first 5 always print timing

    // Claude Generated (Sep 2026): setForcePhaseTiming(true) (from -md_diagnostics_timing)
    // used to be silently overridden every step by the unconditional
    // `setRecordKernelTimings(verbosity >= 2)` at the top of calculateEnergy() below, so a
    // GFN-FF MD run at the normal (< 2) MD verbosity never recorded GPU kernel timings even
    // though it had asked to. Mirrors GFNFF::m_force_phase_timing's do_timing OR-pattern.
    bool m_force_phase_timing = false;

    /// Claude Generated (Sep 2026): WP-A on-device D4 pair-list build from the host generator's
    /// current CN + Gaussian weights (factored out of initGPUWorkspace() so the skin-triggered
    /// rebuild of GFNFF::updateDispersionPairsIfNeeded() can repeat it).
    void buildDispersionPairsOnDevice();
    /// Claude Generated (Oct 2026): compare the device H-bond list with detectHydrogenBondsNative() (sets, case 1
    /// with A < B); prints a summary to stderr and warns on a difference (-gfnff.gpu_hb_detect_verify).
    void verifyDeviceHBList(const std::vector<int>& di, const std::vector<int>& dj, const std::vector<int>& dk,
                            const std::vector<int>& dcase, const std::vector<int>& dparent);
    bool m_disp_pairs_on_device = false;  ///< WP-A mode active (gpu_disp_pairs_on_device)
    bool m_disp_pairs_requested = false;  ///< the D4 pair list is built on the device (setting, backend default, ROCm env)
    bool m_hb_detect_on_device = false;   ///< Oct 2026: H-bond triples enumerated on the device (gpu_hb_detect_on_device)
    bool m_hb_detect_verify = false;      ///< compare the device list with the host detection at every refresh
    bool m_rep_pairs_on_device = false;   ///< Oct 2026: non-bonded repulsion list rebuilt on the device (gpu_rep_pairs_on_device)
};

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

/// Claude Generated (Sep 2026, F-16): CURCUMA_GFNFF_GPU_LEAK=1 restores the old deliberate leak
/// of the host GFNFF object and parameter set (safety valve for the unresolved F-Q9 heap issue).
inline bool gfnffGpuLeakRequested()
{
    static const bool leak = [] {
        const char* e = std::getenv("CURCUMA_GFNFF_GPU_LEAK");
        return e && *e && std::string(e) != "0";
    }();
    return leak;
}

template <class Backend>
GFNFFGpuMethodImpl<Backend>::GFNFFGpuMethodImpl(const std::string& method_name,
                                                      const json& config)
    : m_parameters(config)
    , m_method_name(method_name)
{
    // Claude Generated (April 2026): Read accuracy-related flags for GPU fallback behavior
    // config.value("gfnff", config) falls back to config itself so top-level CLI params work.
    json gfnff_cfg = config.value("gfnff", config);
    m_allow_unconverged_charges = gfnff_cfg.value("allow_unconverged_charges", false);   // re-read below with the eeq_solver scope first
    m_skip_phase2 = gfnff_cfg.value("skip_phase2", false);                               // re-read below with the eeq_solver scope first
    // eeq_rmsd_threshold: per-atom RMSD (Bohr) below which the Cholesky factor is reused.
    // 0.0 (default) = always refactorize, matches reference CPU behavior exactly.
    // Suggested MD value: 0.05–0.10 Bohr (saves ~12 ms/step when geometry barely changes).
    m_eeq_rmsd_threshold = gfnff_cfg.value("eeq_rmsd_threshold", 0.0);

    // EEQ Coulomb-matrix distance cutoff (Bohr). Default 0 = no cutoff, matches Fortran
    // goed_gfnff (gfnff_engrad.F90:1274-1391) and the CPU EEQSolver default. Non-zero
    // values violate Hellmann-Feynman vs. the full Coulomb energy → MD energy drift.
    // Read below through eeq_value() (eeq_solver scope first).

    // WP7-B/C (May 2026): EEQ solver strategy for nfrag>1.
    //   "cholesky" / "schur_cholesky" → WP5-A/WP7-A (exact, default).
    //   "batched"                     → WP7-B per-fragment Cholesky (drops cross-fragment Coulomb).
    //   "pcg"                         → WP7-C iterative PCG (warm-started).
    //   "auto"                        → EEQSolver::resolveMethod (ppcg by eeq_ppcg_min_atoms/min_nfrag, else cholesky).
    //   "lu" / "ldlt"                 → device LU / LDL^T (CUDA only; other backends reject them, EEQSolver::gpuSupports).
    // WP7-E: reuses the CPU EEQSolver's ppcg PARAMs (eeq_solver.h) for both the iteration
    // cap and the Auto-strategy nfrag/atom thresholds, so CPU and GPU agree on when to
    // prefer projected PCG. Note m_eeq_ppcg_tolerance is an ABSOLUTE tolerance on |r| (GPU
    // convention, like m_eeq_pcg_tolerance above) whereas the CPU's eeq_ppcg_tol is relative
    // (tol*(|b|+1)) — same PARAM value, different meaning; read as pcg_tolerance's GPU sibling.
    // Claude Generated (Sep 2026): these PARAMs belong to eeq_solver (eeq_solver.h) and the CPU
    // EEQSolver reads them from that scope. Read them from there first, so the same flag
    // (-eeq_solver.eeq_ppcg_min_nfrag 0 = exact solve) steers CPU and GPU alike.
    auto eeq_value = [&](const char* key, auto dflt) {
        const json* scopes[2] = {
            gfnff_cfg.contains("eeq_solver") ? &gfnff_cfg.at("eeq_solver") : nullptr,
            config.contains("eeq_solver") ? &config.at("eeq_solver") : nullptr };
        for (const json* scope : scopes) {
            if (scope && scope->is_object() && scope->contains(key))
                return scope->value(key, dflt);
        }
        return gfnff_cfg.value(key, dflt);
    };
    // Claude Generated (Oct 26, 2026): -eeq_solver.skip_phase2 and -eeq_solver.allow_unconverged_charges steered the CPU
    // EEQSolver but not this wrapper, which read the gfnff scope only (GPU MD of peo201-chain: -eeq_solver.skip_phase2 true
    // changed the first CPU step and left the device EEQ on, -gfnff.skip_phase2 true bypassed it).
    m_allow_unconverged_charges = eeq_value("allow_unconverged_charges", false);
    m_skip_phase2 = eeq_value("skip_phase2", false);
    // Claude Generated (Sep 28, 2026, docs/MULTI_GPU_GAPS.md X-1): the solver keys below are
    // eeq_solver PARAMs since the parameter parser registers multi-line PARAMs, so a flat flag
    // (-pcg_tolerance 1e-12) is routed to the eeq_solver scope. Read that scope first, like the
    // CPU solver, then gfnff; the GPU keeps its own fallback defaults.
    {
        std::string strategy_str = eeq_value("solve_method", gfnff_cfg.value("solve", std::string("auto")));
        m_eeq_strategy = EEQSolver::parseSolveMethod(strategy_str);
    }
    m_eeq_distance_cutoff = eeq_value("eeq_distance_cutoff", 0.0);   // 0 = no cutoff (see above)
    m_eeq_batched_max_error = eeq_value("eeq_batched_max_error_eh", 1e-6);
    m_eeq_pcg_max_iter   = eeq_value("max_pcg_iterations", 200);
    m_eeq_pcg_tolerance  = eeq_value("pcg_tolerance", 1e-10);
    m_eeq_block_jacobi_max_frag_atoms = eeq_value("gpu_block_jacobi_max_frag_atoms", 300);
    m_eeq_block_jacobi_max_nfrag      = eeq_value("gpu_block_jacobi_max_nfrag", 400);
    m_eeq_ppcg_max_iter  = eeq_value("eeq_ppcg_max_iter", 500);
    m_eeq_ppcg_tolerance = eeq_value("pcg_tolerance", 1e-10);
    m_eeq_ppcg_min_nfrag = eeq_value("eeq_ppcg_min_nfrag", 1);
    m_eeq_ppcg_min_atoms = eeq_value("eeq_ppcg_min_atoms", 500);
    {   // Oct 2026: local-inverse preconditioner (same PARAMs and validation as the CPU EEQSolver)
        const std::string pc = eeq_value("eeq_ppcg_precond", std::string("localinv"));
        if (pc == "localinv") m_eeq_localinv = true;
        else if (pc != "diag")
            throw std::invalid_argument(fmt::format("eeq_ppcg_precond '{}' unknown (valid: diag, localinv)", pc));
        m_eeq_localinv_radius  = eeq_value("eeq_localinv_radius_bohr", 8.0);
        m_eeq_localinv_ext     = std::max(eeq_value("eeq_localinv_ext_bohr", 10.0), m_eeq_localinv_radius);
        m_eeq_localinv_rebuild = eeq_value("eeq_localinv_rebuild_bohr", 1.0);
        if (config.contains("threads") && config["threads"].is_number())
            m_eeq_localinv_threads = std::max(1, static_cast<int>(config["threads"].get<double>()));
        m_host_threads = m_eeq_localinv_threads;
    }

    // Claude Generated (Sep 2026, multi-GPU): `gpu_device` (global CLI key, set per worker by
    // the batch capabilities) pins this method to one device. Everything below that allocates
    // on the device (workspace, EEQ solver, constant tables) is created after this bind.
    // An invalid index throws; the plugin entry point turns that into nullptr -> CPU fallback.
    if (config.contains("gpu_device")) {
        const auto& v = config["gpu_device"];
        if (v.is_number()) m_gpu_device = static_cast<int>(v.get<double>());   // CLI stores 3.0
        else if (v.is_string() && !v.get<std::string>().empty()) m_gpu_device = std::stoi(v.get<std::string>());
    }
    if (m_gpu_device >= 0) {
        const int count = Backend::deviceCount();
        if (m_gpu_device >= count || !Backend::setDevice(m_gpu_device))
            throw std::runtime_error(std::string(Backend::name) + " device " + std::to_string(m_gpu_device)
                                     + " is not available (" + std::to_string(count) + " visible)");
    }

    m_gfnff = std::make_unique<GFNFF>(config);
}

// ---------------------------------------------------------------------------
// Destructor — controlled teardown to survive CUDA heap corruption
// ---------------------------------------------------------------------------

template <class Backend>
GFNFFGpuMethodImpl<Backend>::~GFNFFGpuMethodImpl()
{
    // Destroy GPU workspace first (holds CUDA resources) - on the device that owns them.
    bindDevice();
    m_gpu_workspace.reset();

    // Claude Generated (Sep 2026, F-16): freed normally now (it used to be leaked, ~140 MB per
    // molecule of 1410 atoms, see m_gpu_params). CURCUMA_GFNFF_GPU_LEAK=1 keeps the old leak.
    if (gfnffGpuLeakRequested()) {
        (void)m_gpu_params.release();
        (void)m_gfnff.release();
    } else {
        m_gpu_params.reset();
        m_gfnff.reset();
    }
}

// ---------------------------------------------------------------------------
// setMolecule — init topology (CPU), then upload parameters to GPU
// ---------------------------------------------------------------------------

template <class Backend>
bool GFNFFGpuMethodImpl<Backend>::setMolecule(const Mol& mol)
{
    bindDevice();
    if (!m_gfnff) {
        m_has_error = true;
        m_error_message = "GFNFFGPUMethod: m_gfnff is nullptr";
        CurcumaLogger::error(m_error_message);
        return false;
    }

    // Store atom types for GPU workspace (from mol, before InitialiseMolecule)
    m_atom_types = mol.m_atoms;

    auto t_init_start = std::chrono::high_resolution_clock::now();

    // WP-A (Jun 2026): if the GPU builds the D4 pair list on device, tell the host
    // generator to skip its O(N^2) pair loop (keep only CN + gw the device reuses).
    // Must be set before InitialiseMolecule triggers host parameter generation.
    {
        // Claude Generated (Oct 14, 2026): default on for CUDA and ROCm (device build in the host pair order; checked on
        // 27 structures on CUDA and 16 on ROCm, docs/GPU_STEP_PAUSES_2026-10-14.md). The gfnff scope wins over the top level.
        bool disp_dev = Backend::default_disp_pairs_on_device;
        if (m_parameters.contains("gfnff") && m_parameters["gfnff"].is_object()
            && m_parameters["gfnff"].contains("gpu_disp_pairs_on_device"))
            disp_dev = m_parameters["gfnff"].value("gpu_disp_pairs_on_device", disp_dev);
        else
            disp_dev = m_parameters.value("gpu_disp_pairs_on_device", disp_dev);
        m_disp_pairs_requested = disp_dev;
        m_gfnff->setSkipHostDispPairs(disp_dev);
        // Claude Generated (Sep 2026): the device enumerates the Coulomb pairs itself;
        // `-gpu_coulomb_implicit false` restores the host pair list.
        bool coul_implicit = m_parameters.value("gpu_coulomb_implicit", true);
        if (m_parameters.contains("gfnff") && m_parameters["gfnff"].is_object())
            coul_implicit = m_parameters["gfnff"].value("gpu_coulomb_implicit", coul_implicit);
        m_gfnff->setImplicitCoulombPairs(coul_implicit);
        m_gfnff->setKeepFullParameterSet(true);  // consumeCachedParameterSet() needs every pair list
        // Claude Generated (Oct 2026): the wrapper never runs the CPU workspace of m_gfnff (no
        // GFNFF::Calculation()), so its HB/XB lists need no update on a re-detection (every backend
        // since the CUDA measurement of Oct 14; first ROCm only).
        m_gfnff->setSkipHostHBWorkspaceUpdate(true);
        // Claude Generated (Oct 2026): the repulsion lists of the wrapper come from getLast*Repulsions();
        // the CPU workspace never runs, so a refresh need not copy them into it (every backend).
        m_gfnff->setSkipHostRepulsionWorkspaceUpdate(true);
    }

    // CPU topology + parameter generation (same as CPU gfnff)
    if (!m_gfnff->InitialiseMolecule(mol)) {
        m_has_error = true;
        m_error_message = "GFNFFGPUMethod: InitialiseMolecule failed";
        CurcumaLogger::error(m_error_message);
        return false;
    }
    auto t_after_topo = std::chrono::high_resolution_clock::now();

    // Upload parameters to GPU + create CPU residual
    if (!initGPUWorkspace()) {
        return false;
    }

    // Oct 2026: from here on the non-bonded repulsion list is rebuilt on the device (default on for CUDA,
    // `-gfnff.gpu_rep_pairs_on_device false` keeps the host build). Set after the setup build on purpose.
    if constexpr (gfnffBackendRepulsionBuild<Backend>::value) {
        bool rep_dev = true;
        if (m_parameters.contains("gpu_rep_pairs_on_device")) rep_dev = m_parameters.value("gpu_rep_pairs_on_device", rep_dev);
        if (m_parameters.contains("gfnff") && m_parameters["gfnff"].is_object())
            rep_dev = m_parameters["gfnff"].value("gpu_rep_pairs_on_device", rep_dev);
        typename GFNFF::RepulsionDeviceInputs probe;
        if (rep_dev && m_gfnff->repulsionDeviceInputs(probe)) {
            m_rep_pairs_on_device = true;
            m_gfnff->setSkipHostRepulsionPairs(true);
        }
    }
    auto t_init_end = std::chrono::high_resolution_clock::now();

    m_gpu_upload_time_ms = std::chrono::duration<double, std::milli>(t_init_end - t_after_topo).count();

    m_initialized = true;

    // Claude Generated (Oct 2026): the GPU wrapper never calls GFNFF::Calculation(), so the
    // frag_charge_model ensemble (CPU default) does not act here: the reference rule would put the whole
    // net charge on the fragment holding atom 1 (GMTKN55 gfnff, ROCm against CPU: 25 of 2462 structures
    // differ, up to 131 kcal/mol, all of them charged). A charged molecule with several fragments therefore
    // asks EnergyCalculator for the CPU method (requestsCpuFallback()); -gfnff.frag_charge_model reference
    // keeps the GPU and makes the CPU use the same rule.
    m_cpu_fallback_requested = (m_eeq_nfrag > 1 && mol.m_charge != 0 && m_gfnff->fragEnsembleRequested());
    // Claude Generated (Oct 2026): the ALPB solvation term exists only in GFNFF::Calculation(), which this
    // wrapper never calls; before this check -gfnff.solvent water on -gpu rocm returned the gas-phase
    // energy without a message (caffeine, 1 run: -4.67273706 against -4.68976245 Eh on the CPU).
    if (m_gfnff->solvationActive()) {
        m_cpu_fallback_requested = true;
        m_cpu_fallback_reason = "GFN-FF: implicit solvent (ALPB) is not applied by the GPU path, calculation on the CPU";
    }

    CitationRegistry::cite("gfnff");
    CitationRegistry::cite("d4", "gfnff");

    return true;
}

// ---------------------------------------------------------------------------
// initGPUWorkspace — split params → GPU workspace + CPU residual workspace
// ---------------------------------------------------------------------------

template <class Backend>
bool GFNFFGpuMethodImpl<Backend>::initGPUWorkspace()
{
    try {
        const std::vector<int>& atom_types = m_atom_types;
        const int natoms = static_cast<int>(atom_types.size());
        if (natoms == 0) {
            m_has_error = true;
            m_error_message = "GFNFFGPUMethod: zero atoms after InitialiseMolecule";
            CurcumaLogger::error(m_error_message);
            return false;
        }

        // Claude Generated (March 2026): GPU memory check before allocation
        // Check available GPU memory before proceeding
        if (!GPUUtils::checkGPUMemoryAvailable(natoms, 0.8)) {
            m_has_error = true;
            m_error_message = fmt::format(
                "Insufficient GPU memory for {} atoms. Use CPU fallback: -method gfnff -gpu none",
                natoms);
            CurcumaLogger::error(m_error_message);
            return false;
        }

        // Log GPU memory status at verbosity >= 2
        GPUUtils::logGPUMemoryStatus(2);

        // Consume pre-generated params from initializeForceField().
        // CRITICAL: Do NOT call generateGFNFFParameterSet() again — causes heap corruption.
        std::unique_ptr<GFNFFParameterSet> pending = m_gfnff->consumeCachedParameterSet();
        if (!pending) {
            m_has_error = true;
            m_error_message = "GFNFFGPUMethod: no cached params (was InitialiseMolecule called?)";
            CurcumaLogger::error(m_error_message);
            return false;
        }

        // Pre-allocate ALL Eigen Vectors and staging buffers BEFORE CUDA init.
        // After the GPU workspace is constructed, the driver may corrupt adjacent heap
        // metadata, making Eigen Vector resize/assign crash. All per-step buffers must
        // be allocated now and reused via memcpy in the hot path.
        m_gpu_cn_final = Vector::Zero(natoms);
        m_cached_gradient = Matrix::Zero(natoms, 3);
        m_gfnff->preAllocateForGPUPath(natoms);

        // Pre-allocate EEQ staging buffers before CUDA init (heap safety).
        // After the GPU workspace is constructed, the driver may corrupt adjacent heap
        // metadata, making any malloc/new crash. All per-step buffers must be allocated NOW.
        {
            int nfrag_actual = m_gfnff->getTopologyInfo().nfrag;
            int nf = std::max(nfrag_actual, 1);
            m_eeq_z1.resize(natoms, 0.0);
            m_eeq_Z2.resize(natoms * nf, 0.0);
            m_eeq_charges_gpu.resize(natoms, 0.0);
            // CPU Schur workspace: [S: nf×nf | rhs: nf | lambda: nf]
            // Avoids Eigen MatrixXd/VectorXd heap allocs in the nfrag>1 hot path.
            m_schur_workspace.assign(nf * (nf + 2), 0.0);
        }

        // === 1. Create GPU workspace from full parameter set ===
        // All energy terms (including HB, XB, ATM, BATM, sTors) run on GPU.
        m_gpu_workspace = std::make_unique<GpuWorkspace>(*pending, natoms, atom_types);

        // Claude Generated (Sep 2026, F-16): keep only what the GPU path reads later (bonds and
        // the bond-HB cross-reference, see calculateEnergy); the rest (D4 pairs, repulsion
        // pairs, ...) is on the device now and is freed here. Old behaviour (leak the whole set)
        // with CURCUMA_GFNFF_GPU_LEAK=1.
        if (gfnffGpuLeakRequested()) {
            (void)m_gpu_params.release();
            m_gpu_params.reset(pending.release());
        } else {
            auto keep = std::make_unique<GFNFFParameterSet>();
            keep->bonds = std::move(pending->bonds);
            keep->bond_hb_data = std::move(pending->bond_hb_data);
            pending.reset();
            m_gpu_params = std::move(keep);
        }

        // Pass verbosity to GPU workspace for conditional diagnostics
        m_gpu_workspace->setVerbosity(CurcumaLogger::get_verbosity());

        // Claude Generated (April 2026): Forward term enable/disable flags to GPU workspace
        // Without this, GPU always computes all terms regardless of config flags.
        // Flags live under m_parameters["gfnff"] (same layout as CPU ForceField path).
        json gfnff_cfg;
        if (m_parameters.contains("gfnff") && m_parameters["gfnff"].is_object())
            gfnff_cfg = m_parameters["gfnff"];
        if (gfnff_cfg.contains("dispersion"))
            m_gpu_workspace->setDispersionEnabled(gfnff_cfg.value("dispersion", true));
        if (gfnff_cfg.contains("hbond"))
            m_gpu_workspace->setHBondEnabled(gfnff_cfg.value("hbond", true));
        if (gfnff_cfg.contains("repulsion"))
            m_gpu_workspace->setRepulsionEnabled(gfnff_cfg.value("repulsion", true));
        if (gfnff_cfg.contains("coulomb"))
            m_gpu_workspace->setCoulombEnabled(gfnff_cfg.value("coulomb", true));
        // G3a (Apr 2026): Optional fixed GPU kernel block size (0 = adaptive default)
        if (gfnff_cfg.contains("gpu_block_size"))
            m_gpu_workspace->setBlockSize(gfnff_cfg.value("gpu_block_size", 0));

        // Task #10 (Jun 2026): CN-derivative pair-list regen + cutoff knobs.
        // Prefer the nested gfnff-scoped value, fall back to top-level m_parameters.
        auto cfg_get_bool = [&](const char* k, bool d) {
            if (gfnff_cfg.contains(k)) return gfnff_cfg.value(k, d);
            return m_parameters.value(k, d);
        };
        auto cfg_get_int = [&](const char* k, int d) {
            if (gfnff_cfg.contains(k)) return gfnff_cfg.value(k, d);
            return m_parameters.value(k, d);
        };
        auto cfg_get_dbl = [&](const char* k, double d) {
            if (gfnff_cfg.contains(k)) return gfnff_cfg.value(k, d);
            return m_parameters.value(k, d);
        };
        m_cn_pair_regen         = cfg_get_bool("gpu_cn_pair_regen", true);
        m_cn_pair_regen_every   = cfg_get_int("gpu_cn_pair_regen_every", 0);
        m_cn_pair_cutoff_factor = cfg_get_dbl("gpu_cn_pair_cutoff_factor", 2.5);
        m_gpu_workspace->setCNPairCutoffFactor(m_cn_pair_cutoff_factor);
        // K1 (Oct 2026, Claude Generated): CN pair list as a Verlet list with skin s.
        // The list holds every pair with r < f*rs + 2s (rs = rcov_i + rcov_j, 4/3-scaled D3 radii,
        // f = gpu_cn_pair_cutoff_factor). A pair NOT in the list had r_build >= f*rs + 2s; after every
        // atom moved by at most D since the build, r >= f*rs + 2s - 2D. The list is rebuilt once
        // D > D_allow = s + rs_min/2 * (f - 1 - x_tail/|kn|), rs_min = smallest rs among the elements
        // present. Then every missing pair has r >= rs*(1 + x_tail/|kn|), so its erf-CN term is at most
        // 0.5*erfc(x_tail) and its dCN/dr at most |kn|/(rs*sqrt(pi))*exp(-x_tail^2). With kn = -7.5 and
        // x_tail = 8: 5.6e-30 per pair (CN) and 4.2e-28 / Bohr per pair (derivative, H-H, rs = 1.61 Bohr).
        // f = 2.5, H-H: rs_min = 1.61 Bohr, D_allow = s + 0.35 Bohr. Legacy (s < 0): rebuild on the
        // topology check only, list cutoff f*rs exactly as before.
        m_cn_pair_skin = cfg_get_dbl("gpu_cn_pair_skin_bohr", 1.0);
        if constexpr (Backend::has_cn_pair_skin) {
            if (m_cn_pair_skin >= 0.0) {
                constexpr double kn_abs = 7.5;   // |kn| of the GFN-FF erf CN (computeCN, k_cn_compute_pairs)
                constexpr double x_tail = 8.0;   // erf argument the missing pairs are guaranteed to exceed
                double rcov_min = 1e30;
                for (int z : m_atom_types) {
                    if (z >= 1 && z <= static_cast<int>(GFNFFParameters::covalent_rad_d3.size()))
                        rcov_min = std::min(rcov_min, GFNFFParameters::covalent_rad_d3[z - 1] * (4.0 / 3.0));
                }
                const double rs_min = 2.0 * rcov_min;
                // For f < 1 + x_tail/|kn| (2.07) even a fresh list misses larger terms; then the margin is
                // 0 and a stale list misses nothing that a fresh list keeps (r >= f*rs for every rs).
                const double c_tail = m_cn_pair_cutoff_factor - 1.0 - x_tail / kn_abs;
                const double margin = c_tail > 0.0 ? 0.5 * rs_min * c_tail : 0.0;
                m_cn_pair_max_disp = std::max(0.0, m_cn_pair_skin + margin);
                m_gpu_workspace->setCNPairSkin(m_cn_pair_skin);
                if (CurcumaLogger::get_verbosity() >= 3)
                    CurcumaLogger::info(fmt::format("  GPU CN pair list: skin {:.3f} Bohr, rebuild above {:.3f} Bohr displacement",
                                                    m_cn_pair_skin, m_cn_pair_max_disp));
            }
        }
        // K2/K3 (Oct 2026, Claude Generated): stream overlap; fallbacks are the PARAM defaults.
        m_eeq_overlap  = cfg_get_bool("gpu_eeq_overlap", true);
        m_disp_overlap = cfg_get_bool("gpu_dispersion_overlap", true);
        // Oct 2026 (Claude Generated): rebuilt D4/repulsion lists reach the device through the
        // staged pinned upload (same bytes, parallel conversion); fallback is the PARAM default.
        if constexpr (Backend::has_staged_list_upload)
            m_gpu_workspace->setStagedListUpload(cfg_get_bool("gpu_staged_list_upload", true), m_host_threads);
        // Claude Generated (Oct 2026): H-bond triples enumerated on the device (CUDA and ROCm, default on).
        if constexpr (Backend::has_device_hb_detect) {
            m_hb_detect_on_device = cfg_get_bool("gpu_hb_detect_on_device", true);
            m_hb_detect_verify = m_hb_detect_on_device && cfg_get_bool("gpu_hb_detect_verify", false);
        }
        // Claude Generated (Sep 2026): per-step D4 C6 refresh on the device, same switch as the
        // CPU path (GFNFF::refreshDispersionC6()).
        m_gpu_workspace->setDispersionC6Update(cfg_get_bool("dispersion_c6_update", true));

        // Deliverable 3 (Jun 2026): route high-fragment EEQ to the exact CPU PCG.
        // Backend default: 0 on CUDA (device Schur handles it), 16 on ROCm.
        // Claude Generated (Sep 2026, MULTI_GPU_GAPS F-13): read only where the backend has
        // no device Schur paths (ROCm). The key is documented as "ROCm only", and on CUDA the
        // host EEQ it selects is the path that segfaults inside a CUDA process on polymer_2x
        // (TECHNICAL_DEBT.md F-Q9), so the CUDA wrapper keeps its trait default (0).
        if constexpr (!Backend::has_device_schur) {
            m_eeq_cpu_fragment_threshold = cfg_get_int("eeq_rocm_cpu_fragment_threshold",
                                                       Backend::default_eeq_cpu_fragment_threshold);
            m_eeq_cpu_route_forced = gfnff_cfg.contains("eeq_rocm_cpu_fragment_threshold")
                                  || m_parameters.contains("eeq_rocm_cpu_fragment_threshold");
        }

        // Phase 2: Upload C6 reference table for GPU dc6dcn per-pair computation
        D4ParameterGenerator* d4 = m_gfnff->getD4Generator();
        if (d4 && !d4->getC6FlatCache().empty()) {
            m_gpu_workspace->uploadC6ReferenceTable(d4->getC6FlatCache(), d4->getRefN());
            // Phase 6: Upload reference CN values for GPU Gaussian weight kernel
            m_gpu_workspace->uploadRefCN(d4->getRefCN());

            // WP-A (Jun 2026): opt-in on-device D4 dispersion pair-list build. Replaces
            // the host GenerateDispersionPairsNative O(N^2) loop + the per-build H2D
            // upload of the pair arrays (the host-uploaded SoA from the ctor is
            // overwritten here). Needs CN + gw on the device first. Default OFF: the
            // proven host build stays the reference; the ROCm workspace additionally rebuilds
            // the gather CSR. See docs/GFNFF_PERFORMANCE_LEVERS.md.
            if (m_disp_pairs_requested) {
                m_disp_pairs_on_device = true;
                buildDispersionPairsOnDevice();
                if (CurcumaLogger::get_verbosity() >= 1)
                    CurcumaLogger::info("GFN-FF GPU: D4 dispersion pair list built on device (gpu_disp_pairs_on_device)");
            }
        }

        // Claude Generated (March 2026): GPU EEQ solver (cuSOLVER/rocSOLVER Cholesky)
        // EEQ staging buffers pre-allocated above (before device init). Create solver here.
        // Lazy refactorization is RMSD-based: force_refactor flag passed per solve() call.
        m_eeq_gpu = std::make_unique<GpuEEQSolver>(natoms);
        m_eeq_priority_set = false;   // K2/K3: new solver stream has default priority

        // Claude Generated (Sep 2026, docs/MULTI_GPU_GAPS.md F-3): split ONE molecule over several
        // GPUs - the Coulomb tiles and the projected-PCG EEQ (matrix column blocks + matvec).
        // gpu_split_devices: auto (default) = the -gpu_devices pool, else every visible device,
        // unless this is a leased batch worker or the run is pinned with an explicit -gpu_device;
        // a comma list; none/off. Used from gpu_split_min_atoms on.
        if constexpr (Backend::has_multi_gpu_split) {
            std::string spec = "auto";
            for (const json* scope : { &gfnff_cfg, &m_parameters })
                if (scope->contains("gpu_split_devices") && (*scope)["gpu_split_devices"].is_string()) {
                    spec = (*scope)["gpu_split_devices"].get<std::string>();
                    break;
                }
            const int min_atoms = cfg_get_int("gpu_split_min_atoms", 1000);
            const int leased = curcuma::leasedGpuDevice();
            const bool auto_dev = m_parameters.value("gpu_device_auto", false);
            const bool pinned = leased < 0 && m_parameters.contains("gpu_device") && !auto_dev;
            std::vector<int> devs;
            if (spec == "auto" || spec == "all") {
                if (leased < 0 && (!pinned || spec == "all")) {
                    const auto& pool = curcuma::GpuDevicePool::instance();
                    if (pool.active()) devs = pool.devices();
                    else for (int d = 0; d < Backend::deviceCount(); ++d) devs.push_back(d);
                }
            } else if (spec != "none" && spec != "off" && !spec.empty()) {
                std::stringstream ss(spec);
                for (std::string tok; std::getline(ss, tok, ',');) {
                    try { devs.push_back(std::stoi(tok)); } catch (...) {
                        CurcumaLogger::warn("gpu_split_devices: cannot parse '" + tok + "'");
                    }
                }
                if (leased >= 0 && !devs.empty()) {
                    CurcumaLogger::warn("-gpu_split_devices ignored in a batch worker (leased device "
                                        + std::to_string(leased) + ")");
                    curcuma::reportGpuFallback("-gpu_split_devices ignored in a batch worker");
                    devs.clear();
                }
            }
            if (devs.size() > 1) {
                m_gpu_workspace->setCoulombSplitDevices(devs, min_atoms);
                m_eeq_gpu->setSplitMixedPrecision(cfg_get_bool("gpu_eeq_mixed_precision", true),
                                                  cfg_get_dbl("gpu_eeq_mixed_rho", 1e-5));
                m_eeq_gpu->setSplitDevices(devs, std::max(min_atoms, cfg_get_int("gpu_split_eeq_min_atoms", 4000)));
                m_split_configured = true;
            }
        }

        // WP-B (Jun 2026): FP32-factor + FP64 iterative-refinement EEQ solve. Opt-in on
        // both backends (default off until per-card tuned, mirroring the xTB
        // scf_mixed_precision X-AP3 default-off-on-CUDA policy). Reads
        // gfnff.eeq_mixed_precision, falling back to the global scf_mixed_precision intent;
        // refine steps via eeq_mixed_precision_iters (default 2, min 1). Targets the
        // factor-dominated few-fragment paths; the nfrag>1 general path stays FP64 (see
        // eeq_solver_gpu.cu / rocm/eeq_solver_hip.hiph).
        {
            bool mp = cfg_get_bool("eeq_mixed_precision",
                                   cfg_get_bool("scf_mixed_precision", false));
            int  mp_iters = cfg_get_int("eeq_mixed_precision_iters", 2);
            m_eeq_gpu->setMixedPrecision(mp, mp_iters);
            if (mp && CurcumaLogger::get_verbosity() >= 1)
                CurcumaLogger::info(fmt::format(
                    "GFN-FF GPU EEQ: FP32 factor + FP64 iterative refinement ({} step(s))",
                    (mp_iters < 1 ? 1 : mp_iters)));
        }

        // WP2: Upload topology-constant EEQ parameters to GPU.
        // cn=Zero because chi_corrected_static and cnf are CN-independent.
        {
            Vector dummy_cn = Vector::Zero(natoms);
            GFNFF::EEQGPUParams topo_params = m_gfnff->prepareEEQParametersForGPU(dummy_cn);
            m_gpu_workspace->uploadEEQTopologyParams(topo_params);
            m_eeq_fraglist         = topo_params.fraglist;
            m_eeq_rhs_constraints  = topo_params.rhs_constraints;
            m_eeq_nfrag            = topo_params.nfrag;

            // WP6: upload fragment topology for batched Cholesky (nfrag > 1), and since Sep 2026
            // also for one fragment on the device-Schur backend, which then takes the projected
            // PCG from eeq_ppcg_min_atoms on like the CPU (F-4; see singleFragmentPPCG()).
            if (m_eeq_nfrag > 1 || (m_eeq_nfrag == 1 && (Backend::has_device_schur || Backend::has_device_ppcg))) {
                m_eeq_gpu->uploadFragmentTopology(m_eeq_nfrag, m_eeq_fraglist, natoms);
                // WP7-B: cache min inter-fragment distance for batched-solver warning.
                m_eeq_batched_check_pending = true;  // batched error estimate once per topology
            }
        }

        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::success(fmt::format(
                "gfnff-gpu: GPU workspace ready ({} atoms, {} bonds, {} disp pairs, GPU EEQ enabled)",
                natoms, m_gpu_workspace->bondCount(), m_gpu_workspace->dispersionCount()));
        }

        return true;

    } catch (const std::exception& e) {
        m_has_error = true;
        m_error_message = std::string(Backend::workspace_name) + " init failed: " + e.what();
        CurcumaLogger::error(m_error_message);
        return false;
    }
}

// ---------------------------------------------------------------------------
// calculateEnergy — orchestrate CN/EEQ/state/calculate (no delegation to GFNFF::Calculation)
// ---------------------------------------------------------------------------

template <class Backend>
double GFNFFGpuMethodImpl<Backend>::calculateEnergy(bool gradient)
{
    bindDevice();
    CitationRegistry::cite("gfnff");
    CitationRegistry::cite("d4", "gfnff");
    CitationRegistry::cite("eeq", "gfnff");
    CitationRegistry::cite("pyykko", "gfnff");
    CitationRegistry::cite("sanderson", "gfnff");
    CitationRegistry::cite("ghosh_islam", "gfnff");
    CitationRegistry::cite("atm", "d3");
    CitationRegistry::cite("bj", "d3");
    CitationRegistry::cite("casimir_polder", "d4");
    if (!m_initialized || !m_gfnff || !m_gpu_workspace) {
        CurcumaLogger::error("GFNFFGPUMethod: not initialized");
        return 0.0;
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== GFNFFGPUMethod::calculateEnergy() START (GPU path) ===");
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    // CURCUMA_GFNFF_STEPPHASE=1: GAP lines name the host sections that no phase timer covers (> 0.5 ms)
    static const bool gap_log = [] { const char* e = std::getenv("CURCUMA_GFNFF_STEPPHASE"); return e && *e && std::string(e) != "0"; }();
    auto gap_t = t0;
    auto gap_lap = [&](const char* what) {
        if (!gap_log) return;
        const auto now = std::chrono::high_resolution_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(now - gap_t).count();
        if (ms > 0.5) fmt::print(stderr, "GAP {} {:.2f}\n", what, ms);
        gap_t = now;
    };

    // Claude Generated (May 2026): Enable per-stream kernel timing at verbosity >= 2.
    // This disables CUDA Graph replay so events get recorded each step.
    // Claude Generated (Sep 2026): OR with m_force_phase_timing (set by
    // setForcePhaseTiming(), driven by -md_diagnostics_timing) — this used to be an
    // unconditional overwrite, so a diagnostics run below verbosity 2 silently lost
    // GPU kernel timing every step even though it had asked to force it on.
    if (m_gpu_workspace) {
        m_gpu_workspace->setRecordKernelTimings(
            CurcumaLogger::get_verbosity() >= 2 || m_force_phase_timing);
        // Static-Mode (WP-S1, May 2026): propagate frozen-state flags so GPU kernels for
        // CN / Gaussian-weights / dc6dcn / EEQ-charge-upload are skipped this step.
        m_gpu_workspace->setStaticFlags(
            m_gfnff->staticCNFrozen(),
            m_gfnff->staticChargesFrozen());
    }

    try {

    // === Phase 3: CPU/GPU Overlap Architecture (Claude Generated March 2026) ===
    //
    // Pipeline:
    //   1. GPU: computeCN → sync (need CN for EEQ)
    //   2. CPU: prepareCNAndEEQ (O(N³) EEQ solver)
    //      GPU: charge-independent kernels run in parallel (dispersion, repulsion,
    //           bonds, angles, dihedrals, inversions, storsions, hbonds, batm, atm, xbonds)
    //   3. CPU done → upload EEQ charges → launch Coulomb + postprocess → download
    //
    // This overlaps the expensive EEQ solver with ~12 GPU kernels that don't need charges.

    const int N = static_cast<int>(m_atom_types.size());
    const GeoGradMatrix& geom_bohr = m_gfnff->getGeometryBohr();

    // Claude Generated (April 2026): Upload PBC unit cell to GPU constant memory
    if (m_gfnff->hasPBC()) {
        Eigen::Matrix3d cell_bohr = m_gfnff->getUnitCellBohr();
        Eigen::Matrix3d cell_bohr_inv = cell_bohr.inverse();
        m_gpu_workspace->setUnitCell(cell_bohr.data(), cell_bohr_inv.data(), true);
    }

    // --- Step 0: GPU topology displacement check (Claude Generated March 2026) ---
    // Runs on GPU where coords already live. Result fed to GFNFF::needsFullTopologyUpdate()
    // to skip CPU O(N) Eigen matrix subtraction.
    const auto t_pre_start = std::chrono::high_resolution_clock::now();
    m_gpu_workspace->setGeometry(geom_bohr);

    // React topology mode (Claude Generated Aug 2026): run the CPU hysteresis scan;
    // when it rebuilt the bonded terms, the whole device workspace is stale (bond/angle/
    // torsion SoAs, repulsion partition, Coulomb params, EEQ fragment topology). Rebuild
    // it from the refreshed cached parameter set by re-running initGPUWorkspace().
    // (Since Sep 2026 the previous parameter set is freed here; it used to be leaked, F-16.)
    m_gfnff->updateReactiveTopologyIfNeeded();
    if (m_gfnff->consumeReactRebuild()) {
        m_gpu_workspace.reset();
        m_eeq_gpu.reset();
        if (!initGPUWorkspace()) {
            m_has_error = true;
            m_error_message = "react topology rebuild: GPU workspace reconstruction failed";
            CurcumaLogger::error(m_error_message);
            return 0.0;
        }
        m_gpu_workspace->setGeometry(geom_bohr);
        m_cn_pairs_generated = false;
        m_eeq_has_ref_geom = false;
    }
    gap_lap("after_setgeometry_react");

    // K1 (Oct 2026, Claude Generated): true when the CN pair list is a Verlet list with its own
    // reference geometry (see initGPUWorkspace() for the threshold derivation).
    bool cn_pair_skin_mode = false;
    if constexpr (Backend::has_cn_pair_skin)
        cn_pair_skin_mode = (m_cn_pair_skin >= 0.0);

    {
        bool needs_topo_update = m_gpu_workspace->checkDisplacement(0.5);
        m_gfnff->setExternalTopologyDecision(needs_topo_update);
        if (needs_topo_update && !cn_pair_skin_mode) {
            // Reference geometry updated after topology recalculation (triggered inside prepareCNAndEEQ)
            // — updateReferenceGeometry() called below after prepareCNAndEEQ.
            // Task #10 (Jun 2026): the CN-derivative pair list is built from the current
            // geometry; when atoms have moved past the displacement threshold it is stale,
            // so drop the latch to force a rebuild in the gradient block below.
            if (m_cn_pair_regen)
                m_cn_pairs_generated = false;
        }
    }
    gap_lap("after_checkDisplacement");

    // K1 (Oct 2026, Claude Generated): rebuild an existing CN pair list BEFORE computeCN when an
    // atom moved more than m_cn_pair_max_disp since its build, so this step's CN and chain rule
    // both see a list that misses no pair with a contribution above 0.5*erfc(8). The first list
    // is still built after computeCN (below), exactly as before, so a single point is unchanged.
    if constexpr (Backend::has_cn_pair_skin) {
        if (cn_pair_skin_mode && m_cn_pair_regen && m_cn_pairs_generated
            && m_gpu_workspace->cnPairListStale(m_cn_pair_max_disp)) {
            m_gpu_workspace->generateCNPairListOnGPU();
        }
    }
    gap_lap("after_K1_cnlist_rebuild");

    // --- Step 1: GPU CN computation (G-P1: truly async, no sync here) ---
    // k_cn_compute + k_dlogdcn + 3 async D2H launches on main stream.
    // d_cn_final / d_dlogdcn are already on GPU and correct after this call.
    // finalizeCNForCPU() is called after Phase 4 launch to sync without sleeping.
    // Static-Mode (WP-S1): skip when frozen_cn — d_cn_final/d_dlogdcn retain previous-step values.
    auto t_cn_start = std::chrono::high_resolution_clock::now();
    m_last_gpu_prep_timing.pre = std::chrono::duration<double, std::milli>(t_cn_start - t_pre_start).count();
    if (!m_gpu_workspace->frozenCN()) {
        m_gpu_workspace->computeCN(m_atom_types);
    }
    // G-P1: removed immediate sync + memcpy + setD3CN here.
    // d_cn_final → d_cn D2D copy happens inside prepareAndLaunchChargeIndependent().
    // setD3CN() no longer needed: d_cn is populated GPU-side from d_cn_final.
    auto t_cn_end = std::chrono::high_resolution_clock::now();

    // --- Step 2a: Set up gradient state (CN pairs) ---
    // G-P1: dlogdcn is already on GPU from k_dlogdcn; setDlogDCN() removed (redundant).
    // d_dlogdcn is used directly by k_cn_chainrule via the main stream ordering.
    if (gradient) {
        // Task #10 (Jun 2026): optional periodic forced regen as an MD safety net for
        // slow drift that never trips the 0.5 Bohr displacement check above.
        if (m_cn_pair_regen_every > 0 && (m_cn_pair_grad_steps % m_cn_pair_regen_every == 0))
            m_cn_pairs_generated = false;
        ++m_cn_pair_grad_steps;

        if (!m_cn_pairs_generated) {
            m_gpu_workspace->generateCNPairListOnGPU();
            m_cn_pairs_generated = true;
        }
        // G-P1: removed: memcpy dlogdcn from pinned + setDlogDCN()
        // k_dlogdcn already wrote d_dlogdcn on main stream before prepareAndLaunch.
    }
    auto t_dlogdcn_end = std::chrono::high_resolution_clock::now();

    // --- Step 2b: Compute per-bond HB coordination numbers ---
    // Claude Generated (Apr 2026): HB CN computed entirely on GPU via two-pass kernel.
    // Replaces CPU loop over bond_hb_data (O(n_hb_pairs) erf() evaluations).
    auto t_hb_cn_start = std::chrono::high_resolution_clock::now();
    if (m_gpu_params && !m_gpu_params->bond_hb_data.empty()) {
        m_gpu_workspace->computeBondHBCN();
    }
    auto t_hb_cn_end = std::chrono::high_resolution_clock::now();

    // --- Step 2c: Dynamic HB/XB re-detection (must happen before kernel launch) ---
    auto t_hbxb_start = std::chrono::high_resolution_clock::now();
    // CURCUMA_GFNFF_STEPPHASE=1: sub-timers of a fired refresh (stderr), Claude Generated (Oct 2026)
    static const bool hb_phase_log = [] { const char* e = std::getenv("CURCUMA_GFNFF_STEPPHASE"); return e && *e && std::string(e) != "0"; }();
    auto hb_t = std::chrono::steady_clock::now();
    auto hb_lap = [&](const char* what) {
        if (!hb_phase_log) return;
        const auto now = std::chrono::steady_clock::now();
        fmt::print(stderr, "HBXBSUB {} {:.2f}\n", what, std::chrono::duration<double, std::milli>(now - hb_t).count());
        hb_t = now;
    };
    bool hbxb_updated = false;
    bool hb_on_device = false;                       // the H-bond list of this refresh was built on the device
    std::vector<int> dev_hb_i, dev_hb_j, dev_hb_k;   // its (A, H, B) triples on the host
    if constexpr (Backend::has_device_hb_detect) {
        if (m_hb_detect_on_device) {
            hb_on_device = true;
            if (m_gfnff->hbxbRefreshDue()) {
                const Vector& topo_q = m_gfnff->getTopologyInfo().topology_charges;
                const auto in = m_gfnff->prepareHBDetectInputs(topo_q);
                std::vector<int> dev_case, dev_parent;
                const int n_hb = m_gpu_workspace->detectHBondsOnDevice(in, dev_hb_i, dev_hb_j, dev_hb_k,
                    m_hb_detect_verify ? &dev_case : nullptr, m_hb_detect_verify ? &dev_parent : nullptr);
                if (m_hb_detect_verify)
                    verifyDeviceHBList(dev_hb_i, dev_hb_j, dev_hb_k, dev_case, dev_parent);
                m_gfnff->commitHBXBDeviceUpdate(n_hb, m_gfnff->detectXBForDevice());
                hbxb_updated = m_gfnff->consumeHBXBUpdate();
            }
        }
    }
    if (!hb_on_device) {
        m_gfnff->updateHBXBIfNeeded(nullptr);
        hbxb_updated = m_gfnff->consumeHBXBUpdate();
    }
    const auto t_hb_det_end = std::chrono::high_resolution_clock::now();
    m_last_gpu_prep_timing.hb_detect = std::chrono::duration<double, std::milli>(t_hb_det_end - t_hbxb_start).count();
    m_last_gpu_prep_timing.hb_rebuild = 0.0;
    m_last_gpu_prep_timing.hb_upload = 0.0;
    if (hbxb_updated) hb_lap("detect");
    if (hbxb_updated) {
        if (!hb_on_device)
            m_gpu_workspace->updateHBonds(m_gfnff->getLastHBonds(), m_atom_types);
        m_gpu_workspace->updateXBonds(m_gfnff->getLastXBonds(), m_atom_types);
        const auto t_hb_up1 = std::chrono::high_resolution_clock::now();
        m_last_gpu_prep_timing.hb_upload = std::chrono::duration<double, std::milli>(t_hb_up1 - t_hb_det_end).count();
        hb_lap("upload_hb_xb");

        // Claude Generated (Apr 2026): Propagate HB re-detection to bond metadata and alpha pairs.
        // Without this, the bond SoA's nr_hb/hb_H_atom and the HB alpha pair list become stale,
        // causing gradient errors at H-bond atoms during optimization/MD.
        // Reference: Fortran gfnff_ini2.f90:1008-1060 (bond_hb_AHB_set0/set1)
        if (m_gpu_params) {
            const auto t_hb_rb0 = std::chrono::high_resolution_clock::now();
            auto hb_update = hb_on_device
                ? m_gfnff->rebuildBondHBDataFromTriples(dev_hb_i, dev_hb_j, dev_hb_k, m_gpu_params->bonds)
                : m_gfnff->rebuildBondHBData(m_gfnff->getLastHBonds(), m_gpu_params->bonds);
            const auto t_hb_rb1 = std::chrono::high_resolution_clock::now();
            m_last_gpu_prep_timing.hb_rebuild = std::chrono::duration<double, std::milli>(t_hb_rb1 - t_hb_rb0).count();
            hb_lap("rebuildBondHBData");

            // Update HB alpha (H,B) pair list on GPU
            m_gpu_workspace->updateHBAlphaPairs(hb_update.bond_hb_data, m_atom_types);
            hb_lap("updateHBAlphaPairs");

            // Update per-bond nr_hb and hb_H_atom metadata on GPU
            m_gpu_workspace->updateBondHBMetadata(hb_update.bond_nr_hb, hb_update.bond_hb_H_atom);
            m_last_gpu_prep_timing.hb_upload += std::chrono::duration<double, std::milli>(
                std::chrono::high_resolution_clock::now() - t_hb_rb1).count();
            hb_lap("updateBondHBMetadata");

            // Update leaked params' bond_hb_data for HB CN computation (Step 2b)
            m_gpu_params->bond_hb_data = std::move(hb_update.bond_hb_data);

            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format(
                    "  [DEBUG] HB re-detection: updated bond_hb_data ({} entries), "
                    "hb_alpha ({} pairs), bond nr_hb for {} bonds",
                    static_cast<int>(m_gpu_params->bond_hb_data.size()),
                    [&]() { int n = 0; for (const auto& e : m_gpu_params->bond_hb_data) n += static_cast<int>(e.B_atoms.size()); return n; }(),
                    static_cast<int>(hb_update.bond_nr_hb.size())));
            }
        }

        // Phase 8: HB/XB SoA n-values changed → captured graph is stale
        m_gpu_workspace->invalidateGraph();
    }
    auto t_hbxb_end = std::chrono::high_resolution_clock::now();
    gap_t = t_hbxb_end;

    // --- Step 2d: Periodic non-bonded repulsion pair-list refresh ---
    // Claude Generated (Sep 2026): the repulsion pair list is otherwise built ONCE from a
    // hard 20 Bohr cutoff at setMolecule() time and never revisited (identical bug on
    // CPU and GPU, since both go through GFNFF::InitialiseMolecule() -> initGPUWorkspace()
    // once). See GFNFF::updateNonbondedRepulsionIfNeeded() for the full rationale.
    const auto t_rep_start = std::chrono::high_resolution_clock::now();
    m_gfnff->updateNonbondedRepulsionIfNeeded(nullptr);
    m_last_gpu_prep_timing.rep_upload = 0.0;
    gap_lap("repulsion_host_build");
    if constexpr (gfnffBackendRepulsionBuild<Backend>::value) {
        if (m_rep_pairs_on_device && m_gfnff->consumeNonbondedRepulsionUpdate()) {
            // The skin rule fired on the host; the list is rebuilt from the device coordinates. The
            // tables depend on the topology only and are cheap to collect, so no cache to invalidate.
            typename GFNFF::RepulsionDeviceInputs rep_in;
            if (m_gfnff->repulsionDeviceInputs(rep_in)) {
                m_gpu_workspace->generateNonbondedRepulsionOnGPU(rep_in);   // device coordinates set at the start of the call
                m_gpu_workspace->invalidateGraph();
            } else {
                m_has_error = true;
                m_error_message = "GFN-FF GPU: repulsion tables unavailable after the device build was enabled";
                CurcumaLogger::error(m_error_message);
                return 0.0;
            }
        }
    }
    if (m_gfnff->consumeNonbondedRepulsionUpdate()) {
        const auto t_rep_up = std::chrono::high_resolution_clock::now();
        m_gpu_workspace->updateRepulsion(m_gfnff->getLastBondedRepulsions(),
                                          m_gfnff->getLastNonbondedRepulsions());
        m_last_gpu_prep_timing.rep_upload = std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - t_rep_up).count();
        m_last_gpu_prep_timing.rep_pairs = static_cast<double>(m_gfnff->getLastNonbondedRepulsions().size());
        // Repulsion SoA n-values changed → captured graph is stale (same reasoning as HB/XB).
        m_gpu_workspace->invalidateGraph();
    }
    gap_lap("repulsion_upload");

    // --- Step 2e: D4 dispersion pair-list skin check + explicit Coulomb list refresh ---
    // Claude Generated (Sep 2026): same one-shot-list bug as the repulsion list above; see
    // GFNFF::updateDispersionPairsIfNeeded() / updateCoulombPairsIfNeeded(). Must precede
    // prepareAndLaunchChargeIndependent() and computeGaussianWeightsOnGPU(), which read the
    // dispersion SoA (k_dc6dcn_per_pair, k_dispersion). Both lists change their SoA n ->
    // the captured CUDA graph is stale and must be invalidated, exactly like HB/XB/repulsion.
    const auto t_disp_start = std::chrono::high_resolution_clock::now();
    m_last_gpu_prep_timing.rep_update = std::chrono::duration<double, std::milli>(t_disp_start - t_rep_start).count();
    gap_lap("before_dispersion_update");
    m_gfnff->updateDispersionPairsIfNeeded(nullptr);
    gap_lap("dispersion_host_generator");
    if (m_gfnff->consumeDispersionPairsUpdate()) {
        if (m_disp_pairs_on_device) {
            buildDispersionPairsOnDevice();   // host generator refreshed CN + weights; device builds pairs
            gap_lap("dispersion_device_build");
        } else if (m_gfnff->getWorkspace()) {
            m_gpu_workspace->updateDispersion(m_gfnff->getWorkspace()->d4Dispersions());
        }
        m_gpu_workspace->invalidateGraph();
    }
    const auto t_coul_start = std::chrono::high_resolution_clock::now();
    m_last_gpu_prep_timing.disp_update = std::chrono::duration<double, std::milli>(t_coul_start - t_disp_start).count();
    gap_lap("dispersion_update");
    m_gfnff->updateCoulombPairsIfNeeded(nullptr);
    if (m_gfnff->consumeCoulombPairsUpdate() && m_gfnff->getWorkspace()) {
        m_gpu_workspace->updateCoulombPairs(m_gfnff->getWorkspace()->coulombPairs());
        m_gpu_workspace->invalidateGraph();
    }
    m_last_gpu_prep_timing.coul_update = std::chrono::duration<double, std::milli>(
        std::chrono::high_resolution_clock::now() - t_coul_start).count();

    // HB/XB pair list consistency: CPU vs GPU (verbosity >= 3)
    if (CurcumaLogger::get_verbosity() >= 3 && !hb_on_device) {
        const auto& hb_cpu = m_gfnff->getLastHBonds();
        const auto& hb_gpu = m_gpu_workspace->getLastHBonds();
        if (hb_cpu.size() != hb_gpu.size()) {
            CurcumaLogger::warn(fmt::format(
                "  HB pair count mismatch: CPU={} GPU={}",
                hb_cpu.size(), hb_gpu.size()));
        } else {
            CurcumaLogger::info(fmt::format("  HB pairs: {} (CPU==GPU)", hb_cpu.size()));
        }
        size_t n_compare = std::min(hb_cpu.size(), hb_gpu.size());
        for (size_t p = 0; p < n_compare && p < 20; ++p) {
            bool same = (hb_cpu[p].i == hb_gpu[p].i &&
                         hb_cpu[p].j == hb_gpu[p].j &&
                         hb_cpu[p].k == hb_gpu[p].k &&
                         hb_cpu[p].case_type == hb_gpu[p].case_type);
            if (!same) {
                CurcumaLogger::warn(fmt::format(
                    "  HB pair[{}] diff: CPU(A={} H={} B={} case={}) GPU(A={} H={} B={} case={})",
                    p, hb_cpu[p].i, hb_cpu[p].j, hb_cpu[p].k, hb_cpu[p].case_type,
                    hb_gpu[p].i, hb_gpu[p].j, hb_gpu[p].k, hb_gpu[p].case_type));
            }
        }

        const auto& xb_cpu = m_gfnff->getLastXBonds();
        const auto& xb_gpu = m_gpu_workspace->getLastXBonds();
        if (xb_cpu.size() != xb_gpu.size()) {
            CurcumaLogger::warn(fmt::format(
                "  XB pair count mismatch: CPU={} GPU={}",
                xb_cpu.size(), xb_gpu.size()));
        } else {
            CurcumaLogger::info(fmt::format("  XB pairs: {} (CPU==GPU)", xb_cpu.size()));
        }
    }

    gap_lap("coulomb_update_and_checks");
    auto t_prep_end = std::chrono::high_resolution_clock::now();

    // === OVERLAP: Launch charge-independent GPU kernels ===
    // These run asynchronously while CPU computes EEQ parameters below.
    // k_dispersion in gradient mode is not in Phase 1 (needs dc6dcn from
    // computeGaussianWeightsOnGPU, which runs AFTER Phase 1 on d_cn — current step); it starts
    // right after dc6dcn (K2) or in Phase 2 (gpu_dispersion_overlap false).
    // G-P1: d_cn_final → d_cn D2D copy is enqueued inside prepareAndLaunchChargeIndependent().
    //
    // Claude Generated (Sep 2026): energy-only calls need the C6 of THIS geometry too (line
    // searches, energy-only scans), as on the CPU (gfnff_method.cpp, dispersionC6Stale). In
    // energy mode k_dispersion is launched inside Phase 1 on stream A, so the refresh has to be
    // enqueued BEFORE it, from d_cn_final (current step, computeCN above); the event_upload
    // fence in Phase 1 orders stream A after it.
    if (!gradient && !m_gpu_workspace->frozenCN()) {
        m_gpu_workspace->computeGaussianWeightsOnGPU(/*use_cn_final=*/true);
    }
    m_gpu_workspace->prepareAndLaunchChargeIndependent(gradient);

    // K3 (Oct 2026, Claude Generated): mark the main-stream point where coords + rhs_atoms are
    // written, BEFORE dc6dcn is enqueued; the device ppcg waits on it instead of on everything.
    bool eeq_input_event = false;
    (void)eeq_input_event;   // unused on backends without stream overlap
    if constexpr (Backend::has_stream_overlap) {
        if (m_eeq_overlap) {
            m_gpu_workspace->recordEEQInputsReady();
            eeq_input_event = true;
        }
    }

    // Gaussian weights + dc6dcn: must run AFTER Phase 1 so d_cn holds the current step's
    // values (Phase 1 D2D-copies d_cn_final → d_cn on main stream before sA fences).
    // Using d_cn guarantees consistency with the gradient-mode dispersion kernel.
    // Static-Mode (WP-S1): skip when frozen_cn — gw/dgw/dc6dcn buffers retain previous values.
    if (gradient && !m_gpu_workspace->frozenCN()) {
        m_gpu_workspace->computeGaussianWeightsOnGPU(/*use_cn_final=*/false);
    }
    // K2 (Oct 2026, Claude Generated): k_dispersion reads no EEQ charge (zetac6 holds the Phase-1
    // topology charges), so it runs on the pairwise stream right after dc6dcn, parallel to the EEQ
    // solve, instead of in Phase 2 after it. Off: Phase 2 launches it as before.
    if constexpr (Backend::has_stream_overlap) {
        if (gradient && m_disp_overlap)
            m_gpu_workspace->launchDispersionGradientEarly();
    }
    auto t_launch_end = std::chrono::high_resolution_clock::now();

    // WP5-D (May 2026): In the normal GPU path (nfrag==1, not skip_phase2) no CPU
    // consumer of CN remains after WP5-B+C. finalizeCNForCPU + prepareCNAndEEQ are
    // kept only for fallback branches that still need CPU-side CN/EEQ.
    // WP6 (May 2026): nfrag > 1 with valid fragment topology also skips CPU path.
    if (m_skip_phase2 || (m_eeq_nfrag != 1 && !m_eeq_gpu->isFragmentTopoValid())) {
        m_gpu_workspace->finalizeCNForCPU(m_gpu_cn_final);
        m_gfnff->prepareCNAndEEQ(gradient, /*gpu_only=*/true, &m_gpu_cn_final, /*skip_eeq=*/true);
    }

    // Update GPU reference geometry ONLY when a full topology recalculation happened.
    // Must NOT update every step — that would make the displacement check always return 0
    // and prevent topology updates during MD (instability bug).
    if (m_gfnff->consumeFullTopologyUpdate()) {
        m_gpu_workspace->updateReferenceGeometry();
        // Phase 8: all SoAs rebuilt after topology change → captured graph is stale
        m_gpu_workspace->invalidateGraph();
        // Topology rebuild means large geometry change → reset EEQ reference geometry
        // so the next step always gets a fresh Cholesky factorization.
        m_eeq_has_ref_geom = false;

        // WP2: Topology-constant EEQ parameters changed — re-upload to GPU
        Vector dummy_cn = Vector::Zero(N);
        GFNFF::EEQGPUParams topo_params = m_gfnff->prepareEEQParametersForGPU(dummy_cn);
        m_gpu_workspace->uploadEEQTopologyParams(topo_params);
        m_eeq_fraglist        = topo_params.fraglist;
        m_eeq_rhs_constraints = topo_params.rhs_constraints;
        m_eeq_nfrag           = topo_params.nfrag;

        // WP6: upload fragment topology for batched Cholesky (nfrag > 1); one fragment too (F-4)
        if (m_eeq_nfrag > 1 || (m_eeq_nfrag == 1 && (Backend::has_device_schur || Backend::has_device_ppcg))) {
            m_eeq_gpu->uploadFragmentTopology(m_eeq_nfrag, m_eeq_fraglist, N);
            // WP7-B: refresh min inter-fragment distance after topology rebuild.
            m_eeq_batched_check_pending = true;  // batched error estimate once per topology
        }
    }

    // WP5-B (May 2026): setCNDerivatives no-op — CNF lives on GPU as eeq_topo.d_cnf.
    // WP5-C (May 2026): D4 skip-check moved to GPU (k_check_dc6dcn_skip in computeCN).
    // WP5-D (May 2026): finalizeCNForCPU + prepareCNAndEEQ removed from normal path.
    // All CN consumers are now GPU-side; no CPU round-trip needed.

    // Claude Generated (April 2026): EEQ charge selection, in priority order (updated Sep 28, 2026;
    // the old list predated WP6/WP7 and the device Schur paths):
    //   1. frozen charges (static mode) or skip_phase2: no Phase-2 solve
    //   2. host EEQ: charges already current for this geometry, the ROCm many-fragment route,
    //      or nfrag > 1 before the fragment topology is on the device
    //   3. device EEQ: nfrag == 1 -> dense Schur-Cholesky (WP5-A; always dense, see
    //      docs/MULTI_GPU_GAPS.md F-4); nfrag > 1 -> WP7-E projected PCG / WP7-A / WP7-B / WP7-C
    //   4. a rejected or failed device solve keeps the last accepted charges (reported as a
    //      GPU fallback, src/core/gpu_fallback.h)
    Vector charges(N);
    // Deliverable 3 (Jun 2026) / Oct 2026, ROCm only (threshold 0 on CUDA): exact CPU EEQ
    // (prepareCNAndEEQ with the configured solver; auto = projected PCG for many fragments) for
    // high fragment counts where the device has no usable solve. The normal device path leaves CN
    // GPU-resident (WP5-D), so m_gpu_cn_final is synced to the current geometry first.
    auto eeqCpuFragmentRoute = [&]() {
        m_gpu_workspace->finalizeCNForCPU(m_gpu_cn_final);
        m_gfnff->prepareCNAndEEQ(gradient, /*gpu_only=*/true, &m_gpu_cn_final, /*skip_eeq=*/false);
        charges = m_gfnff->getLastCharges();
        m_gpu_workspace->setEEQCharges(charges);
    };
    // Claude Generated (Oct 2026): a ROCm many-fragment system whose solver resolves to ppcg takes
    // the device projected PCG (WP7-E) instead of the CPU route, unless the threshold was given
    // explicitly (m_eeq_cpu_route_forced). Same resolution rule as the dispatch below.
    const bool many_frag_cpu_route = m_eeq_cpu_fragment_threshold > 0
        && m_eeq_nfrag >= m_eeq_cpu_fragment_threshold;
    bool many_frag_device_ppcg = false;
    if constexpr (Backend::has_device_ppcg) {
        many_frag_device_ppcg = many_frag_cpu_route && !m_eeq_cpu_route_forced
            && m_gpu_workspace->isEEQTopoValid() && m_eeq_gpu->isFragmentTopoValid()
            && EEQSolver::resolveMethod(m_eeq_strategy, N, m_eeq_nfrag, m_eeq_ppcg_min_atoms,
                                        m_eeq_ppcg_min_nfrag) == EEQSolveMethod::ProjectedPCG;
    }
    bool cpu_route_after_device = false;   ///< device ppcg failed on a many-fragment ROCm system
    if (m_gpu_workspace->frozenCharges()) {
        // Static-Mode (WP-S1): GPU d_charges retain previous-step values, no upload needed.
        // CPU-side charges fetched for any consumer that still needs them.
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("EEQ GPU Phase 2: frozen charges (static_charges mode) - skipping solve+upload");
        }
        charges = m_gfnff->getLastCharges();
    } else if (m_skip_phase2) {
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::warn("GPU path: skip_phase2=true - bypassing GPU EEQ, using CPU Phase 1 topology charges");
        }
        m_gfnff->prepareCNAndEEQ(gradient, /*gpu_only=*/true, &m_gpu_cn_final, /*skip_eeq=*/false);
        charges = m_gfnff->getLastCharges();
    } else if (m_gfnff->areEEQChargesCurrent(m_gfnff->getGeometryBohr())) {
        // Step A (Claude Generated, May 2026): Geometry unchanged since last EEQ solve.
        // Reuse the CPU charges (which the topology-stage CPU EEQ already computed,
        // including LU fallback for indefinite matrices) instead of running a fresh
        // GPU Cholesky. Mirrors the CPU "Skipping redundant Phase-2 EEQ" branch in
        // prepareCNAndEEQ (gfnff_method.cpp:761-765) and avoids the indefinite-matrix
        // divergence on systems like polymer.xyz where dpotrf cannot factor A.
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info(
                "EEQ GPU Phase 2: reusing cached CPU charges (geometry unchanged)");
        }
        charges = m_gfnff->getLastCharges();
        m_gpu_workspace->setEEQCharges(charges);
    } else if (many_frag_cpu_route && !many_frag_device_ppcg) {
        // Deliverable 3 (Jun 2026), ROCm only (threshold 0 on CUDA): high fragment count
        // (solvent box) and the solver does not resolve to ppcg (or the threshold was given
        // explicitly). The device EEQ would be a dense N x N Cholesky plus a host Schur step over
        // nfrag+1 columns, intractable here; the device PCG/batched/general-Schur variants are not
        // ported on ROCm. Route to the CPU EEQ, which keeps the charges identical to the CPU path.
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info(fmt::format(
                "EEQ GPU Phase 2: nfrag={} >= {} - routing to the exact CPU EEQ (device dense Cholesky avoided)",
                m_eeq_nfrag, m_eeq_cpu_fragment_threshold));
        }
        eeqCpuFragmentRoute();
    } else {
        // === GPU EEQ: Build Coulomb matrix + Cholesky solve on GPU ===
        // WP2: use topology-constant alpha/gam/chi/cnf already on GPU;
        // rhs_atoms built by k_build_eeq_rhs (launched in prepareAndLaunchChargeIndependent).
        // fraglist and nfrag are cached from last topology build (m_eeq_fraglist/m_eeq_nfrag).
        const bool use_device_rhs = m_gpu_workspace->isEEQTopoValid();

        // Guard: z1, Z2 and charges buffers must be large enough
        if (static_cast<int>(m_eeq_z1.size()) < N)
            m_eeq_z1.resize(N, 0.0);
        if (static_cast<int>(m_eeq_Z2.size()) < N * m_eeq_nfrag)
            m_eeq_Z2.resize(N * m_eeq_nfrag, 0.0);
        if (static_cast<int>(m_eeq_charges_gpu.size()) < N)
            m_eeq_charges_gpu.resize(N, 0.0);

        // GPU EEQ paths:
        //   nfrag == 1: WP5-A single-fragment Cholesky (always applicable)
        //   nfrag >  1: WP6 batched per-fragment Cholesky (applicable once topo uploaded)
        const bool gpu_eeq_applicable = (m_eeq_nfrag == 1)
            || (m_eeq_nfrag > 1 && m_eeq_gpu->isFragmentTopoValid());

        bool used_gpu_eeq = false;
        bool gpu_charges_accepted = true; // cleared if downloaded charges are NaN/|q|>50
        const double rhs_c0 = m_eeq_rhs_constraints.empty() ? 0.0 : m_eeq_rhs_constraints[0];

        // EEQ distance cutoff (matches CPU EEQSolver behaviour at eeq_solver.cpp:2641-2643).
        // CPU truncates Coulomb-matrix off-diagonals at 30 Bohr only when natoms > 200.
        // Below that threshold both paths must use 0.0 to keep small-molecule results
        // bit-identical to the previous behaviour.
        const double eeq_cutoff_sq = (N > 200 && m_eeq_distance_cutoff > 0.0)
                                       ? m_eeq_distance_cutoff * m_eeq_distance_cutoff
                                       : 0.0;

        // nfrag > 1: which solver. Default (cholesky) -> WP7-E above the CPU's ppcg thresholds,
        // WP7-A below them. "batched" -> WP7-B. "pcg" -> WP7-C (explicit only). "ppcg" -> WP7-E.
        // Claude Generated (Sep 2026): resolved ONCE here for both the verbosity-2 path report
        // and the dispatch below - they were two copies, and the report missed the new default.
        //
        // The default does what the CPU EEQSolver does with it: its SchurCholesky dispatch
        // switches to projected PCG for a large many-fragment system (eeq_solver.cpp
        // ppcg_auto: nfrag >= eeq_ppcg_min_nfrag && N >= eeq_ppcg_min_atoms). The GPU used to
        // keep the exact dense WP7-A there: a different algorithm than the CPU and 17x slower
        // on polymer_2x (1325 vs 79 ms per solve). The algorithm difference is NOT what makes
        // CPU and GPU MD trajectories part: on polymer_2x (2 ps, CSVR) two CPU runs with
        // identical settings separate by 1e-5 A at ~730 fs, CPU vs GPU WP7-A at ~700 fs, CPU vs
        // GPU WP7-E at ~710 fs - chaotic growth of rounding differences either way.
        // -eeq_solver.eeq_ppcg_min_nfrag 0 forces the exact solve, on both sides.
        // Oct 2026: the SAME rule as the CPU (EEQSolver::resolveMethod): an explicit method runs as
        // requested, auto picks projected PCG / Cholesky by eeq_ppcg_min_atoms / eeq_ppcg_min_nfrag.
        // The GPU used to carry its own copy of this logic, including an Auto branch of its own.
        const EEQSolveMethod resolved_method = EEQSolver::resolveMethod(
            m_eeq_strategy, N, m_eeq_nfrag, m_eeq_ppcg_min_atoms, m_eeq_ppcg_min_nfrag);
        const EEQSolveMethod many_frag_method = resolved_method;

        // One fragment: projected PCG (WP7-E) when it is the resolved method; WP5-A stays the fallback.
        // It is also what lets the multi-GPU split (gpu_split_devices) reach single-fragment molecules.
        const bool single_frag_ppcg = Backend::has_device_ppcg && m_eeq_nfrag == 1
            && resolved_method == EEQSolveMethod::ProjectedPCG
            && m_eeq_gpu->isFragmentTopoValid();
        const bool single_frag_pcg = Backend::has_device_schur && m_eeq_nfrag == 1
            && resolved_method == EEQSolveMethod::PCG
            && m_eeq_gpu->isFragmentTopoValid();

        if (CurcumaLogger::get_verbosity() >= 1) {
            const char* path = "CPU-fallback";
            if (gpu_eeq_applicable) {
                if constexpr (!Backend::has_device_schur) {
                    // ROCm: no device pcg/batched/Schur/lu/ldlt (EEQSolver::gpuSupports allows auto, cholesky
                    // and ppcg). Oct 2026: projected PCG with the diagonal Jacobi or the host-built local
                    // inverse (eeq_ppcg_precond); every other solve is the dense device Cholesky of WP2
                    // plus the host Schur step.
                    if (resolved_method == EEQSolveMethod::ProjectedPCG && m_eeq_gpu->isFragmentTopoValid()) {
                        if (m_eeq_localinv && localInverseDeviceBuildSelected(N))
                            path = m_eeq_nfrag == 1 ? "WP7-E projected PCG, one fragment, precond=localinv gpu build=device"
                                                    : "WP7-E projected PCG, precond=localinv gpu build=device";
                        else if (m_eeq_localinv && m_eeq_nfrag <= kLocalInvMaxNfrag)
                            path = m_eeq_nfrag == 1 ? "WP7-E projected PCG, one fragment, precond=localinv gpu build=host"
                                                    : "WP7-E projected PCG, precond=localinv gpu build=host";
                        else
                            path = m_eeq_nfrag == 1 ? "WP7-E projected PCG, one fragment, diagonal Jacobi"
                                                    : "WP7-E projected PCG, diagonal Jacobi";
                    } else
                        path = "WP2 dense Cholesky + CPU Schur";
                } else if (resolved_method == EEQSolveMethod::LU) {
                    path = "augmented LU (getrf)";
                } else if (resolved_method == EEQSolveMethod::LDLT) {
                    path = "LDL^T (sytrf) + Schur";
                } else if (m_eeq_nfrag == 1) {
                    path = single_frag_ppcg ? "WP7-E projected PCG, one fragment"
                         : single_frag_pcg  ? "WP7-C PCG, one fragment"
                         : resolved_method == EEQSolveMethod::Batched ? "WP5-A Cholesky (batched is exact for one fragment)"
                                            : "WP5-A Cholesky";
                } else {
                    if      (resolved_method == EEQSolveMethod::Batched)      path = "WP7-B batched";
                    else if (resolved_method == EEQSolveMethod::ProjectedPCG) path = "WP7-E projected PCG";
                    else if (resolved_method == EEQSolveMethod::PCG)          path = "WP7-C PCG";
                    else                                                       path = "WP7-A Cholesky";
                }
            }
            // Once per run (and again if the path changes), not every MD/opt step.
            // On ROCm only ppcg and the dense Cholesky run on the device; anything else is shown as cholesky.
            const EEQSolveMethod shown_method = (Backend::has_device_schur
                                                 || resolved_method == EEQSolveMethod::ProjectedPCG)
                ? resolved_method : EEQSolveMethod::SchurCholesky;
            const std::string line = fmt::format("EEQ solver (final charges): {}{} [GPU {}, {}], N={}, nfrag={}",
                EEQSolver::methodName(shown_method),
                m_eeq_strategy == EEQSolveMethod::Auto ? " (auto)" : "",
                Backend::name, path, N, m_eeq_nfrag);
            EEQSolver::reportSolverLine(line);
        }

        // K3 (Oct 2026): one-shot input event for this step's device solve (cleared below).
        if constexpr (Backend::has_stream_overlap) {
            m_eeq_gpu->setInputReadyEvent(eeq_input_event ? m_gpu_workspace->eeqInputsReadyEvent() : nullptr);
            // The CG is the critical path; dc6dcn/dispersion on default-priority streams fill its gaps.
            if (!m_eeq_priority_set && (m_eeq_overlap || m_disp_overlap)) {
                m_eeq_priority_set = true;
                if (!m_eeq_gpu->setStreamPriorityHigh(true) && CurcumaLogger::get_verbosity() >= 2)
                    CurcumaLogger::info("EEQ GPU: high stream priority not available, default priority kept");
            }
        }
        if (gpu_eeq_applicable) {
            // RMSD-based EEQ lazy refactorization:
            // Compute per-atom RMSD (Bohr) from geometry at last full Cholesky build.
            // If RMSD < threshold: reuse cached L (dpotrs only, saves ~12 ms/step).
            // If RMSD >= threshold or threshold==0: full refactorization.
            const GeoGradMatrix& cur_geom = m_gfnff->getGeometryBohr();
            bool force_refactor = true;
            if (m_eeq_rmsd_threshold > 0.0 && m_eeq_has_ref_geom
                    && m_eeq_ref_geom.rows() == cur_geom.rows()) {
                double rmsd_sq = 0.0;
                for (int i = 0; i < N; ++i) {
                    double dx = cur_geom(i,0) - m_eeq_ref_geom(i,0);
                    double dy = cur_geom(i,1) - m_eeq_ref_geom(i,1);
                    double dz = cur_geom(i,2) - m_eeq_ref_geom(i,2);
                    rmsd_sq += dx*dx + dy*dy + dz*dz;
                }
                force_refactor = (std::sqrt(rmsd_sq / N) >= m_eeq_rmsd_threshold);
            }

            bool eeq_ok = false;
            bool used_gpu_schur = false;

            if (use_device_rhs) {
                // Claude Generated (Oct 2026): solve_method lu -> dense LU of the augmented
                // (N+nfrag) system on the device, any nfrag (same semantics as the CPU lu).
                // Only the CUDA backend has it (EEQSolver::gpuSupports rejects lu elsewhere);
                // the discarded branch is not instantiated for the ROCm solver type.
                if constexpr (Backend::has_device_schur) {
                    if (resolved_method == EEQSolveMethod::LU && m_eeq_gpu->isFragmentTopoValid()) {
                        eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUAugmentedLU(
                            N, m_eeq_nfrag,
                            m_gpu_workspace->getDeviceXPtr(),
                            m_gpu_workspace->getDeviceYPtr(),
                            m_gpu_workspace->getDeviceZPtr(),
                            m_gpu_workspace->getDeviceAlphaPtr(),
                            m_gpu_workspace->getDeviceGamPtr(),
                            m_gpu_workspace->getDeviceRHSPtr(),
                            m_eeq_rhs_constraints,
                            eeq_cutoff_sq);
                        if (eeq_ok) {
                            used_gpu_schur = true;
                        } else {
                            CurcumaLogger::warn("EEQ GPU: augmented LU failed, using the dense Cholesky path for this step");
                            curcuma::reportGpuFallback("GFN-FF EEQ: lu failed, dense Cholesky step");
                        }
                    }
                    // Oct 2026: solve_method ldlt -> Bunch-Kaufman LDL^T of A_nn + Schur (CUDA only).
                    if (resolved_method == EEQSolveMethod::LDLT && m_eeq_gpu->isFragmentTopoValid()) {
                        eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUSchurLDLT(
                            N, m_eeq_nfrag,
                            m_gpu_workspace->getDeviceXPtr(),
                            m_gpu_workspace->getDeviceYPtr(),
                            m_gpu_workspace->getDeviceZPtr(),
                            m_gpu_workspace->getDeviceAlphaPtr(),
                            m_gpu_workspace->getDeviceGamPtr(),
                            m_gpu_workspace->getDeviceRHSPtr(),
                            m_eeq_rhs_constraints,
                            eeq_cutoff_sq);
                        if (eeq_ok) {
                            used_gpu_schur = true;
                        } else {
                            CurcumaLogger::warn("EEQ GPU: LDL^T failed, using the dense Cholesky path for this step");
                            curcuma::reportGpuFallback("GFN-FF EEQ: ldlt failed, dense Cholesky step");
                        }
                    }
                }
                if (eeq_ok) {
                    // lu solved above
                } else if (m_eeq_nfrag == 1) {
                  if (single_frag_ppcg) {
                    // F-4 (Sep 2026): one fragment, projected PCG as on the CPU; WP5-A on failure.
                    prepareLocalInverseGPU(N, eeq_cutoff_sq);
                    eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUProjectedPCG(
                        N, 1,
                        m_gpu_workspace->getDeviceXPtr(),
                        m_gpu_workspace->getDeviceYPtr(),
                        m_gpu_workspace->getDeviceZPtr(),
                        m_gpu_workspace->getDeviceAlphaPtr(),
                        m_gpu_workspace->getDeviceGamPtr(),
                        m_gpu_workspace->getDeviceRHSPtr(),
                        m_eeq_fraglist,
                        m_eeq_rhs_constraints,
                        m_eeq_ppcg_max_iter,
                        m_eeq_ppcg_tolerance,
                        eeq_cutoff_sq,
                        force_refactor);
                    finishLocalInverseGPU();
                    if (eeq_ok) used_gpu_schur = true;
                    else {
                        CurcumaLogger::warn("EEQ GPU: projected PCG (one fragment) stalled, using WP5-A");
                        curcuma::reportGpuFallback("GFN-FF EEQ: projected PCG stalled, dense Cholesky step");
                    }
                  }
                  if (!eeq_ok && single_frag_pcg) {
                    // Explicit pcg on one fragment: WP7-C with nfrag = 1 (diagonal Jacobi).
                    eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUPCG(
                        N, 1,
                        m_gpu_workspace->getDeviceXPtr(),
                        m_gpu_workspace->getDeviceYPtr(),
                        m_gpu_workspace->getDeviceZPtr(),
                        m_gpu_workspace->getDeviceAlphaPtr(),
                        m_gpu_workspace->getDeviceGamPtr(),
                        m_gpu_workspace->getDeviceRHSPtr(),
                        m_eeq_fraglist,
                        m_eeq_rhs_constraints,
                        m_eeq_pcg_max_iter,
                        m_eeq_pcg_tolerance,
                        eeq_cutoff_sq,
                        force_refactor,
                        m_eeq_block_jacobi_max_frag_atoms,
                        m_eeq_block_jacobi_max_nfrag);
                    if (eeq_ok) used_gpu_schur = true;
                    else {
                        CurcumaLogger::warn("EEQ GPU: PCG (one fragment) stalled, using WP5-A");
                        curcuma::reportGpuFallback("GFN-FF EEQ: PCG stalled, dense Cholesky step");
                    }
                  }
                  if (!eeq_ok) {
                    // WP5-A: single-fragment GPU Schur (fast path).
                    eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUSchur(
                        N, 1,
                        m_gpu_workspace->getDeviceXPtr(),
                        m_gpu_workspace->getDeviceYPtr(),
                        m_gpu_workspace->getDeviceZPtr(),
                        m_gpu_workspace->getDeviceAlphaPtr(),
                        m_gpu_workspace->getDeviceGamPtr(),
                        m_gpu_workspace->getDeviceRHSPtr(),
                        m_eeq_fraglist,
                        rhs_c0,
                        eeq_cutoff_sq,
                        force_refactor);
                    if (eeq_ok) {
                        used_gpu_schur = true;
                    } else {
                        if constexpr (Backend::has_device_schur) {
                            CurcumaLogger::warn("EEQ GPU: WP5-A failed, falling back to WP2 + CPU Schur");
                            curcuma::reportGpuFallback("GFN-FF EEQ: device Schur failed, host Schur step");
                        } else if (CurcumaLogger::get_verbosity() >= 3) {
                            CurcumaLogger::info("EEQ ROCm: device GPU-Schur (nfrag=1) not ported; using WP2 GPU solve + CPU Schur");
                        }
                        // Cholesky failed (not SPD): fall back to WP2 + CPU Schur
                        eeq_ok = m_eeq_gpu->solveWithDeviceRHS(
                            N, 1,
                            m_gpu_workspace->getDeviceXPtr(),
                            m_gpu_workspace->getDeviceYPtr(),
                            m_gpu_workspace->getDeviceZPtr(),
                            m_gpu_workspace->getDeviceAlphaPtr(),
                            m_gpu_workspace->getDeviceGamPtr(),
                            m_gpu_workspace->getDeviceRHSPtr(),
                            m_gpu_workspace->getDeviceRHSConstraintsPtr(),
                            m_eeq_fraglist,
                            m_eeq_z1.data(),
                            m_eeq_Z2.data(),
                            eeq_cutoff_sq,
                            force_refactor);
                    }
                  }
                } else {
                    // nfrag > 1: strategy resolved above (many_frag_method). All fall back to
                    // WP2 + CPU-Schur on failure.
                    const EEQSolveMethod resolved = many_frag_method;

                    if (resolved == EEQSolveMethod::ProjectedPCG && m_eeq_gpu->isFragmentTopoValid()) {
                        // WP7-E: single projected-PCG solve (no nfrag+1 loop, no block-Jacobi).
                        prepareLocalInverseGPU(N, eeq_cutoff_sq);
                        eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUProjectedPCG(
                            N, m_eeq_nfrag,
                            m_gpu_workspace->getDeviceXPtr(),
                            m_gpu_workspace->getDeviceYPtr(),
                            m_gpu_workspace->getDeviceZPtr(),
                            m_gpu_workspace->getDeviceAlphaPtr(),
                            m_gpu_workspace->getDeviceGamPtr(),
                            m_gpu_workspace->getDeviceRHSPtr(),
                            m_eeq_fraglist,
                            m_eeq_rhs_constraints,
                            m_eeq_ppcg_max_iter,
                            m_eeq_ppcg_tolerance,
                            eeq_cutoff_sq,
                            force_refactor);
                        finishLocalInverseGPU();
                        if (eeq_ok) {
                            used_gpu_schur = true;
                        } else if (many_frag_device_ppcg) {
                            // ROCm many-fragment system (Claude Generated, Oct 2026): no dense fallback
                            // (nfrag+1 right-hand sides plus a host Schur step), the exact CPU EEQ instead.
                            cpu_route_after_device = true;
                            CurcumaLogger::warn("EEQ GPU: WP7-E projected PCG stalled, using the CPU EEQ for this step");
                            curcuma::reportGpuFallback("GFN-FF EEQ: projected PCG stalled, CPU EEQ step");
                        } else {
                            CurcumaLogger::warn("EEQ GPU: WP7-E projected PCG stalled, falling through to WP7-A");
                            curcuma::reportGpuFallback("GFN-FF EEQ: projected PCG stalled, dense Cholesky step");
                        }
                        // On stall: silently fall through to WP7-A cholesky below.
                    }
                    if (!eeq_ok && resolved == EEQSolveMethod::PCG && m_eeq_gpu->isFragmentTopoValid()) {
                        // WP7-C: iterative PCG with warm-start (explicit solve_method=pcg only —
                        // Auto prefers WP7-E above; see PARAM gpu_block_jacobi_max_nfrag doc).
                        eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUPCG(
                            N, m_eeq_nfrag,
                            m_gpu_workspace->getDeviceXPtr(),
                            m_gpu_workspace->getDeviceYPtr(),
                            m_gpu_workspace->getDeviceZPtr(),
                            m_gpu_workspace->getDeviceAlphaPtr(),
                            m_gpu_workspace->getDeviceGamPtr(),
                            m_gpu_workspace->getDeviceRHSPtr(),
                            m_eeq_fraglist,
                            m_eeq_rhs_constraints,
                            m_eeq_pcg_max_iter,
                            m_eeq_pcg_tolerance,
                            eeq_cutoff_sq,
                            force_refactor,
                            m_eeq_block_jacobi_max_frag_atoms,
                            m_eeq_block_jacobi_max_nfrag);
                        if (eeq_ok) {
                            used_gpu_schur = true;
                        } else {
                            CurcumaLogger::warn("EEQ GPU: WP7-C PCG stalled, falling through to WP7-A");
                            curcuma::reportGpuFallback("GFN-FF EEQ: PCG stalled, dense Cholesky step");
                        }
                        // On stall: silently fall through to WP7-A cholesky below.
                    }
                    if (!eeq_ok && resolved == EEQSolveMethod::Batched
                            && m_eeq_gpu->isFragmentTopoValid()) {
                        // WP7-B: per-fragment Cholesky (drops cross-fragment Coulomb).
                        // Approximate; the error is estimated after the solve (energy criterion).
                        eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUSchurBatched(
                            N, m_eeq_nfrag,
                            m_gpu_workspace->getDeviceXPtr(),
                            m_gpu_workspace->getDeviceYPtr(),
                            m_gpu_workspace->getDeviceZPtr(),
                            m_gpu_workspace->getDeviceAlphaPtr(),
                            m_gpu_workspace->getDeviceGamPtr(),
                            m_gpu_workspace->getDeviceRHSPtr(),
                            m_eeq_rhs_constraints,
                            eeq_cutoff_sq,
                            force_refactor);
                        if (eeq_ok) {
                            used_gpu_schur = true;
                        } else {
                            CurcumaLogger::warn("EEQ GPU: WP7-B batched failed, falling through to WP7-A");
                            curcuma::reportGpuFallback("GFN-FF EEQ: batched failed, dense Cholesky step");
                        }
                    }
                    // WP7-A: full N×N Cholesky + GPU Schur complement for nfrag > 1.
                    // Replaces D2H of z1+Z2 (~N·(1+nfrag) doubles) with D2H of
                    // Cz1+S (nfrag + nfrag² doubles) and a tiny CPU Gauss-elim.
                    // Cross-fragment Coulomb correctly retained (single N×N Cholesky).
                    if (!eeq_ok && !cpu_route_after_device && m_eeq_gpu->isFragmentTopoValid()) {
                        eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUSchurGeneral(
                            N, m_eeq_nfrag,
                            m_gpu_workspace->getDeviceXPtr(),
                            m_gpu_workspace->getDeviceYPtr(),
                            m_gpu_workspace->getDeviceZPtr(),
                            m_gpu_workspace->getDeviceAlphaPtr(),
                            m_gpu_workspace->getDeviceGamPtr(),
                            m_gpu_workspace->getDeviceRHSPtr(),
                            m_eeq_fraglist,
                            m_eeq_rhs_constraints,
                            eeq_cutoff_sq,
                            force_refactor);
                        if (eeq_ok) {
                            used_gpu_schur = true;
                        } else {
                            if constexpr (Backend::has_device_schur) {
                                CurcumaLogger::warn("EEQ GPU: WP7-A cholesky failed, falling back to WP2 + CPU Schur");
                                curcuma::reportGpuFallback("GFN-FF EEQ: device Schur failed, host Schur step");
                            } else if (CurcumaLogger::get_verbosity() >= 3) {
                                CurcumaLogger::info("EEQ ROCm: device GPU-Schur not ported; using WP2 GPU solve + CPU Schur");
                            }
                        }
                    }
                    if (!eeq_ok && !cpu_route_after_device) {
                        // Fallback: WP2 exact GPU solve + CPU Schur complement.
                        // Triggers when fragment topo invalid OR Cholesky failed in WP7-A.
                        if constexpr (Backend::has_device_schur) {
                            CurcumaLogger::warn("EEQ GPU: all GPU Schur paths failed, falling back to WP2 + CPU Schur");
                            curcuma::reportGpuFallback("GFN-FF EEQ: device Schur failed, host Schur step");
                        } else if (CurcumaLogger::get_verbosity() >= 3) {
                            CurcumaLogger::info("EEQ ROCm: GPU-Schur not ported; using WP2 GPU solve + CPU Schur complement");
                        }
                        eeq_ok = m_eeq_gpu->solveWithDeviceRHS(
                            N, m_eeq_nfrag,
                            m_gpu_workspace->getDeviceXPtr(),
                            m_gpu_workspace->getDeviceYPtr(),
                            m_gpu_workspace->getDeviceZPtr(),
                            m_gpu_workspace->getDeviceAlphaPtr(),
                            m_gpu_workspace->getDeviceGamPtr(),
                            m_gpu_workspace->getDeviceRHSPtr(),
                            nullptr,  // d_rhs_constraints not used by WP2
                            m_eeq_fraglist,
                            m_eeq_z1.data(),
                            m_eeq_Z2.data(),
                            eeq_cutoff_sq,
                            force_refactor);
                        // On further failure: fall through to CPU EEQ at line 622
                    }
                }
            } else {
                // Fallback: EEQ topo not yet uploaded to GPU (should not happen after init)
                CurcumaLogger::warn("EEQ GPU: device RHS topo not valid, falling back to CPU param upload");
                const Vector& cn = m_gfnff->getLastCN();
                GFNFF::EEQGPUParams eeq_params = m_gfnff->prepareEEQParametersForGPU(cn);
                eeq_ok = m_eeq_gpu->solve(
                    N, eeq_params.nfrag,
                    m_gpu_workspace->getDeviceXPtr(),
                    m_gpu_workspace->getDeviceYPtr(),
                    m_gpu_workspace->getDeviceZPtr(),
                    eeq_params.alpha_corrected.data(),
                    eeq_params.gam_corrected.data(),
                    eeq_params.fraglist,
                    eeq_params.rhs_atoms.data(),
                    eeq_params.rhs_constraints.data(),
                    m_eeq_z1.data(),
                    m_eeq_Z2.data(),
                    eeq_cutoff_sq,
                    force_refactor);
                // update cached topology data in case it changed
                m_eeq_fraglist        = eeq_params.fraglist;
                m_eeq_rhs_constraints = eeq_params.rhs_constraints;
                m_eeq_nfrag           = eeq_params.nfrag;
            }

            // On successful refactorization: cache reference geometry for next RMSD check.
            if (eeq_ok && force_refactor) {
                m_eeq_ref_geom = cur_geom;
                m_eeq_has_ref_geom = true;
            }

            if (eeq_ok) {
                if (used_gpu_schur) {
                    // Download for storeChargesFromGPU() — EEQ stream already synced.
                    if constexpr (gfnffBackendEEQChargeDownload<Backend>::value)
                        Backend::downloadEEQCharges(m_eeq_charges_gpu.data(),
                                                    m_eeq_gpu->getDeviceChargesPtr(), N);
                    else
                        Backend::downloadDoubles(m_eeq_charges_gpu.data(),
                                                 m_eeq_gpu->getDeviceChargesPtr(), N);
                    // Test hook (Claude Generated, Sep 2026): CURCUMA_EEQ_GPU_FORCE_REJECT=1
                    // rejects every GPU EEQ solution, so the fallback below can be exercised on
                    // a well-conditioned molecule. Read once; zero cost when unset.
                    static const bool force_reject = [] {
                        const char* e = std::getenv("CURCUMA_EEQ_GPU_FORCE_REJECT");
                        return e && *e && std::string(e) != "0";
                    }();
                    if (force_reject) gpu_charges_accepted = false;
                    // Validate — mirrors gfnff_method.cpp:864-880 CPU check.
                    // GPU Cholesky can return info=0 with NaN/huge charges on ill-conditioned systems.
                    for (int i = 0; i < N && gpu_charges_accepted; ++i)
                        if (std::isnan(m_eeq_charges_gpu[i]) || std::isinf(m_eeq_charges_gpu[i]))
                            gpu_charges_accepted = false;
                    if (gpu_charges_accepted) {
                        double max_q = 0.0, min_q = 1e300, mean_q = 0.0;
                        for (int i = 0; i < N; ++i) {
                            double q = m_eeq_charges_gpu[i];
                            max_q = std::max(max_q, std::abs(q));
                            min_q = std::min(min_q, q);
                            mean_q += q;
                        }
                        mean_q /= N;
                        if (max_q > 50.0) gpu_charges_accepted = false;
                        CurcumaLogger::info(fmt::format(
                            "EEQ GPU charges stats: min={:.6f} max_abs={:.6f} mean={:.6f}",
                            min_q, max_q, mean_q));
                    }
                    if (gpu_charges_accepted) {
                        m_gfnff->storeChargesFromGPU(m_eeq_charges_gpu.data(), N);
                        // WP5-A: D2D copy into impl.d_charges — skips the H2D upload in
                        // Phase 2. Only AFTER validation (Claude Generated, Sep 2026): this copy
                        // used to precede the check, so a rejected solution stayed on the device
                        // and Coulomb ran on it while the log reported the Phase-1 fallback.
                        m_gpu_workspace->setEEQDeviceCharges(m_eeq_gpu->getDeviceChargesPtr());
                        if (resolved_method == EEQSolveMethod::Batched && m_eeq_nfrag > 1
                            && m_eeq_batched_check_pending && m_eeq_batched_max_error > 0.0
                            && CurcumaLogger::get_verbosity() >= 1) {   // diagnostic only, see the CPU path
                            m_eeq_batched_check_pending = false;
                            estimateBatchedErrorOnHost(N, eeq_cutoff_sq);
                        }
                    } else {
                        CurcumaLogger::warn("EEQ GPU: invalid charges (NaN/Inf or |q|>50), solution rejected");
                    }
                } else {
                    // CPU Schur complement: q = z1 - Z2 * λ
                    if (m_eeq_nfrag == 1) {
                        double S = 0.0, Cz1 = 0.0;
                        for (int i = 0; i < N; ++i) {
                            S   += m_eeq_Z2[i];
                            Cz1 += m_eeq_z1[i];
                        }
                        double lambda = (Cz1 - rhs_c0) / S;
                        for (int i = 0; i < N; ++i)
                            m_eeq_charges_gpu[i] = m_eeq_z1[i] - m_eeq_Z2[i] * lambda;
                    } else {
                        // Multi-fragment Schur: S·λ = C·z1 - qfrag, q = z1 - Z2·λ
                        // Uses pre-allocated m_schur_workspace to avoid any heap allocation
                        // after device init (the GPU workspace corrupts adjacent heap metadata).
                        const int nf = m_eeq_nfrag;
                        // Grow workspace only if nfrag increased (rare topology change)
                        if (nf * (nf + 2) > static_cast<int>(m_schur_workspace.size()))
                            m_schur_workspace.assign(nf * (nf + 2), 0.0);
                        double* Sw   = m_schur_workspace.data();           // [nf×nf] row-major
                        double* rhsw = m_schur_workspace.data() + nf * nf; // [nf]
                        double* lamw = rhsw + nf;                           // [nf]

                        // Build S(f,g) = Σ_{i∈frag_f} Z2[i + g*N]
                        for (int f = 0; f < nf; ++f) {
                            for (int g = 0; g < nf; ++g) {
                                double sum = 0.0;
                                for (int i = 0; i < N; ++i)
                                    if (m_eeq_fraglist[i] == f + 1)
                                        sum += m_eeq_Z2[i + g * N];
                                Sw[f * nf + g] = sum;
                            }
                        }
                        // Build rhs(f) = Σ_{i∈frag_f} z1[i] - qfrag[f]
                        for (int f = 0; f < nf; ++f) {
                            double sum = 0.0;
                            for (int i = 0; i < N; ++i)
                                if (m_eeq_fraglist[i] == f + 1)
                                    sum += m_eeq_z1[i];
                            double qfrag_f = (f < static_cast<int>(m_eeq_rhs_constraints.size()))
                                                 ? m_eeq_rhs_constraints[f] : 0.0;
                            rhsw[f] = sum - qfrag_f;
                        }
                        // Gaussian elimination with partial pivoting (in-place on Sw/rhsw)
                        for (int col = 0; col < nf; ++col) {
                            int pivot = col;
                            double maxv = std::abs(Sw[col * nf + col]);
                            for (int row = col + 1; row < nf; ++row) {
                                double v = std::abs(Sw[row * nf + col]);
                                if (v > maxv) { maxv = v; pivot = row; }
                            }
                            if (pivot != col) {
                                for (int k = 0; k < nf; ++k)
                                    std::swap(Sw[col * nf + k], Sw[pivot * nf + k]);
                                std::swap(rhsw[col], rhsw[pivot]);
                            }
                            double inv_d = (std::abs(Sw[col * nf + col]) > 1e-15)
                                               ? 1.0 / Sw[col * nf + col] : 0.0;
                            for (int row = col + 1; row < nf; ++row) {
                                double fac = Sw[row * nf + col] * inv_d;
                                for (int k = col + 1; k < nf; ++k)
                                    Sw[row * nf + k] -= fac * Sw[col * nf + k];
                                rhsw[row] -= fac * rhsw[col];
                                Sw[row * nf + col] = 0.0;
                            }
                        }
                        // Back substitution → lamw
                        for (int i = nf - 1; i >= 0; --i) {
                            lamw[i] = rhsw[i];
                            for (int j = i + 1; j < nf; ++j)
                                lamw[i] -= Sw[i * nf + j] * lamw[j];
                            double diag = Sw[i * nf + i];
                            lamw[i] = (std::abs(diag) > 1e-15) ? lamw[i] / diag : 0.0;
                        }
                        // Apply q[i] = z1[i] - Σ_f Z2[i+f*N] * λ[f]
                        for (int i = 0; i < N; ++i) {
                            double q = m_eeq_z1[i];
                            for (int f = 0; f < nf; ++f)
                                q -= m_eeq_Z2[i + f * N] * lamw[f];
                            m_eeq_charges_gpu[i] = q;
                        }
                    }
                    charges = Eigen::Map<const Vector>(m_eeq_charges_gpu.data(), N);
                    m_gfnff->storeChargesFromGPU(m_eeq_charges_gpu.data(), N);
                    m_gpu_workspace->setEEQCharges(charges);
                }
                if (CurcumaLogger::get_verbosity() >= 3) {
                    std::string frag_info;
                    for (int f = 0; f < m_eeq_nfrag; ++f) {
                        double qsum = 0.0;
                        for (int i = 0; i < N; ++i)
                            if (m_eeq_fraglist[i] == f + 1) qsum += m_eeq_charges_gpu[i];
                        frag_info += fmt::format(" frag{}={:+.4f}", f + 1, qsum);
                    }
                    CurcumaLogger::success(fmt::format("EEQ GPU done:{}", frag_info));
                }
                // Step B (May 2026): warn once when the indefinite-matrix LU
                // fallback kicks in. Mirrors the CPU dispatcher line
                // "EEQ PCG: ill-conditioned matrix ... using LU directly".
                if (m_eeq_gpu->isUsingLUFallback() && force_refactor
                        && CurcumaLogger::get_verbosity() >= 1) {
                    CurcumaLogger::warn(fmt::format(
                        "EEQ GPU: matrix indefinite (Cholesky info={}), using LU fallback",
                        m_eeq_gpu->getLastCholInfo()));
                }
                if (gpu_charges_accepted) used_gpu_eeq = true;
            } else if (cpu_route_after_device) {
                eeqCpuFragmentRoute();   // reported above (projected PCG stalled)
                used_gpu_eeq = true;     // charges are the exact CPU solution, not a fallback to old ones
            } else {
                CurcumaLogger::warn("GPU EEQ Cholesky failed (not SPD), falling back to CPU solver");
                curcuma::reportGpuFallback("GFN-FF EEQ: device Cholesky failed (not SPD)");
            }
        } else {
            // Fragment topo not yet uploaded (first MD step before topology build).
            // CPU fallback below will populate m_eeq_fraglist for next step.
        }
        // K3 (Oct 2026): the input event belongs to this step only.
        if constexpr (Backend::has_stream_overlap) {
            m_eeq_gpu->setInputReadyEvent(nullptr);
        }

        if (!used_gpu_eeq) {
            // Fallback: skip the EEQ solve and keep the charges of the last ACCEPTED solve
            // (getLastCharges() is not updated with skip_eeq). Same as the CPU path, which keeps
            // the previous m_charges on an invalid solution (gfnff_method.cpp, "keeping previous
            // m_charges"). CN derivatives needed for gradient are still computed by
            // prepareCNAndEEQ. Claude Generated (Sep 2026): the message used to say "Phase 1
            // topology charges", which is not what is used (measured on triose: the fallback
            // energy matches neither the Phase-1 nor the current Phase-2 single point).
            CurcumaLogger::warn("GPU EEQ: no valid solution - keeping the charges of the last accepted EEQ solve");
            // Counted only when the device solve was applicable and got rejected; the first-step
            // "fragment topology not yet uploaded" route is by design (F-19, Sep 2026).
            if (gpu_eeq_applicable)
                curcuma::reportGpuFallback("GFN-FF EEQ: no valid device solution, previous charges kept");
            m_gfnff->prepareCNAndEEQ(gradient, /*gpu_only=*/true, &m_gpu_cn_final, /*skip_eeq=*/true);
            charges = m_gfnff->getLastCharges();
            m_gpu_workspace->setEEQCharges(charges);
        }
    }

    auto t_eeq_end = std::chrono::high_resolution_clock::now();

    // === Finish: Coulomb + postprocess + download ===
    // Note: setEEQCharges / setEEQDeviceCharges already called above per path.
    m_last_energy = m_gpu_workspace->launchChargeDependentAndFinish(gradient);

    // Claude Generated (Sep 2026): say once which parts of the multi-GPU split ran (or why not).
    if constexpr (Backend::has_multi_gpu_split) {
        if (m_split_configured && !m_split_reported && m_calc_count >= 1) {
            m_split_reported = true;
            const std::string c = m_gpu_workspace->coulombSplitStatus();
            const std::string e = m_eeq_gpu ? m_eeq_gpu->splitStatus() : std::string();
            for (const std::string& st : { c, e }) {
                if (st.rfind("unavailable", 0) == 0) {
                    CurcumaLogger::warn("GFN-FF multi-GPU: " + st);
                    curcuma::reportGpuFallback("GFN-FF multi-GPU split unavailable", st);
                } else if (!st.empty() && CurcumaLogger::get_verbosity() >= 2) {
                    CurcumaLogger::info("GFN-FF multi-GPU: " + st);
                }
            }
        }
    }

    // Claude Generated (June 2026): OPT-IN device-resident HB charges (default OFF).
    // The reference (Fortran) and CPU both freeze the H-bond charges at topology build,
    // so the GPU's frozen HB charges already match; refreshing to live charges makes MD
    // diverge (verified). No-op unless CURCUMA_GFNFF_GPU_RESIDENT_HBQ=1. See
    // GpuWorkspace::refreshHBChargesFromDevice().
    m_gpu_workspace->refreshHBChargesFromDevice();
    auto t4 = std::chrono::high_resolution_clock::now();

    // Claude Generated (Sep 2026): store the per-phase breakdown for getLastPrepTiming()
    // (WP-P1 MD diagnostics). These are the SAME std::chrono timestamps the verbosity>=2
    // printout below already computed durations from on every step — storing them here is
    // a handful of subtractions, not a new measurement, and costs nothing beyond what this
    // function already paid. "eeq_solve" is the actual EEQ Phase 2 solve, wherever it ran
    // (GPU Schur/PCG in the normal case, CPU PCG/Cholesky fallback for the high-fragment or
    // skip_phase2 branches above) — this is what m_gfnff->getLastPrepTiming().eeq_solve
    // cannot report on this path (prepareCNAndEEQ is only ever called here with
    // skip_eeq=true, so its own eeq_solve timer never fires).
    m_last_gpu_prep_timing.cn        = std::chrono::duration<double, std::milli>(t_cn_end - t_cn_start).count();
    m_last_gpu_prep_timing.dlogdcn   = std::chrono::duration<double, std::milli>(t_dlogdcn_end - t_cn_end).count();
    m_last_gpu_prep_timing.hb_cn     = std::chrono::duration<double, std::milli>(t_hb_cn_end - t_dlogdcn_end).count();
    m_last_gpu_prep_timing.hbxb      = std::chrono::duration<double, std::milli>(t_hbxb_end - t_hbxb_start).count();
    m_last_gpu_prep_timing.launch    = std::chrono::duration<double, std::milli>(t_launch_end - t_prep_end).count();
    m_last_gpu_prep_timing.eeq_solve = std::chrono::duration<double, std::milli>(t_eeq_end - t_launch_end).count();
    m_last_gpu_prep_timing.coulomb_postprocess = std::chrono::duration<double, std::milli>(t4 - t_eeq_end).count();
    m_last_gpu_prep_timing.total     = std::chrono::duration<double, std::milli>(t4 - t0).count();

    ++m_calc_count;

    // Claude Generated (Oct 2026): CURCUMA_GFNFF_STEPPHASE=1 prints one line per calculator call
    // (ms): cn, dlogdcn, hb_cn, hbxb (fired flag), launch, eeq_solve, coulomb+postprocess, total.
    {
        static const bool phase_log = [] {
            const char* e = std::getenv("CURCUMA_GFNFF_STEPPHASE");
            return e && *e && std::string(e) != "0";
        }();
        if (phase_log) {
            const auto& T = m_last_gpu_prep_timing;
            fmt::print(stderr, "STEPPHASE {} grad={} cn={:.2f} dlogdcn={:.2f} hb_cn={:.2f} hbxb={:.2f}{} launch={:.2f} eeq={:.2f} coul={:.2f} total={:.2f}\n",
                       m_calc_count, gradient ? 1 : 0, T.cn, T.dlogdcn, T.hb_cn, T.hbxb, hbxb_updated ? "[F]" : "",
                       T.launch, T.eeq_solve, T.coulomb_postprocess, T.total);
        }
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        // Claude Generated (Apr 2026): Granular GPU timing + wallclock-normalized breakdown
        // Fix (May 2026): t_hbxb now measures actual updateHBXBIfNeeded duration (was ≈0).
        double t_gpu_cn    = std::chrono::duration<double, std::milli>(t_cn_end - t_cn_start).count();
        double t_dlogdcn   = std::chrono::duration<double, std::milli>(t_dlogdcn_end - t_cn_end).count();
        double t_hb_cn     = std::chrono::duration<double, std::milli>(t_hb_cn_end - t_dlogdcn_end).count();
        double t_hbxb      = std::chrono::duration<double, std::milli>(t_hbxb_end - t_hbxb_start).count();
        double t_launch    = std::chrono::duration<double, std::milli>(t_launch_end - t_prep_end).count();
        double t_cpu_eeq   = std::chrono::duration<double, std::milli>(t_eeq_end - t_launch_end).count();
        double t_coulomb   = std::chrono::duration<double, std::milli>(t4 - t_eeq_end).count();
        double t_calc      = std::chrono::duration<double, std::milli>(t4 - t0).count();
        double t_init      = m_gfnff->getTopologyTimeMs() + m_gfnff->getParamGenTimeMs();
        double t_wall      = t_calc + t_init;

        auto pct = [&](double t) { return t_wall > 0 ? 100.0 * t / t_wall : 0.0; };

        CurcumaLogger::result("\n[RESULT] GFN-FF (GPU) Timing Breakdown:");
        CurcumaLogger::result("  Phase                          Time (ms)   %Wall");
        CurcumaLogger::result("  --------------------------------------------------");
        CurcumaLogger::result(fmt::format("  GPU CN compute                 {:>10.1f}  {:>5.1f}%", t_gpu_cn, pct(t_gpu_cn)));
        if (gradient)
            CurcumaLogger::result(fmt::format("  dlogdcn compute                {:>10.1f}  {:>5.1f}%", t_dlogdcn, pct(t_dlogdcn)));
        CurcumaLogger::result(fmt::format("  HB CN per-bond loop            {:>10.1f}  {:>5.1f}%", t_hb_cn, pct(t_hb_cn)));
        CurcumaLogger::result(fmt::format("  HB/XB detection{}               {:>10.1f}  {:>5.1f}%",
            hbxb_updated ? "[fired]" : "       ", t_hbxb, pct(t_hbxb)));
        CurcumaLogger::result(fmt::format("  GPU kernel launch overhead     {:>10.1f}  {:>5.1f}%", t_launch, pct(t_launch)));
        CurcumaLogger::result(fmt::format("  CPU EEQ (parallel to GPU)      {:>10.1f}  {:>5.1f}%", t_cpu_eeq, pct(t_cpu_eeq)));
        CurcumaLogger::result(fmt::format("  Coulomb + postprocess + DMA    {:>10.1f}  {:>5.1f}%", t_coulomb, pct(t_coulomb)));
        CurcumaLogger::result("  --------------------------------------------------");
        CurcumaLogger::result(fmt::format("  Calculation Total              {:>10.1f}  {:>5.1f}%", t_calc, pct(t_calc)));

        double seq_time = t_cpu_eeq + t_coulomb;
        double overlap_time = std::max(0.0, t_cpu_eeq - t_launch);
        CurcumaLogger::result("\n[RESULT] GPU Pipeline Analysis:");
        CurcumaLogger::result(fmt::format(
            "  Sequential work (EEQ + finish):             {:.1f} ms ({:.1f}% of calc)",
            seq_time, t_calc > 0 ? 100.0 * seq_time / t_calc : 0.0));
        CurcumaLogger::result(fmt::format(
            "  CPU/GPU overlap (EEQ || charge-indep):    {:.1f} ms ({:.1f}% of calc)",
            overlap_time, t_calc > 0 ? 100.0 * overlap_time / t_calc : 0.0));
    }

    // Claude Generated (May 2026): Unified verbosity-2 GFN-FF report (GPU path).
    // Same struct + format as CPU path. GPU kernel timings come from CUDA events
    // (stream-level grouping — see GpuWorkspace::kernelTimings()).
    if (CurcumaLogger::get_verbosity() >= 2) {
        const auto& comp = m_gpu_workspace->energyComponents();
        GFNFFEnergyReport rep;
        rep.is_gpu = true;

        rep.bond          = comp.bond;
        rep.angle         = comp.angle;
        rep.dihedral      = comp.dihedral;
        rep.inversion     = comp.inversion;
        rep.stors         = comp.stors;
        rep.dispersion    = comp.dispersion;
        rep.bonded_rep    = comp.bonded_rep;
        rep.nonbonded_rep = comp.nonbonded_rep;
        rep.coulomb       = comp.coulomb;
        rep.hbond         = comp.hbond;
        rep.xbond         = comp.xbond;
        rep.atm           = comp.atm;
        rep.batm          = comp.batm;
        rep.total         = m_last_energy;

        // Per-component GPU timings (stream-level: components in same stream share the time).
        // Stream A: dispersion + repulsion;  Stream B: bonded;  Stream C: 3-body;  P2: coulomb
        const auto& kt = m_gpu_workspace->kernelTimings();
        rep.t_bond.gpu          = kt.bonds;
        rep.t_angle.gpu         = kt.angles;
        rep.t_dihedral.gpu      = kt.dihedrals;
        rep.t_inversion.gpu     = kt.inversions;
        rep.t_stors.gpu         = kt.stors;
        rep.t_dispersion.gpu    = kt.dispersion;
        rep.t_bonded_rep.gpu    = kt.bonded_rep;
        rep.t_nonbonded_rep.gpu = kt.nonbonded_rep;
        rep.t_coulomb.gpu       = kt.coulomb;
        rep.t_hbond.gpu         = kt.hbond;
        rep.t_xbond.gpu         = kt.xbond;
        rep.t_atm.gpu           = kt.atm;
        rep.t_batm.gpu          = kt.batm;

        // Phase summary — repurpose the locally-measured GPU phase wall-clocks
        double t_gpu_cn_ms   = std::chrono::duration<double, std::milli>(t_cn_end - t_cn_start).count();
        double t_cpu_eeq_ms  = std::chrono::duration<double, std::milli>(t_eeq_end - t_launch_end).count();
        double t_phase2_ms   = std::chrono::duration<double, std::milli>(t4 - t_eeq_end).count();
        double t_calc_ms     = std::chrono::duration<double, std::milli>(t4 - t0).count();

        rep.t_gpu_cn           = t_gpu_cn_ms;
        rep.t_cpu_eeq_gpu_path = t_cpu_eeq_ms;
        rep.t_gpu_phase2       = t_phase2_ms;
        rep.t_wall             = t_calc_ms;
        rep.t_topology         = m_gfnff->getTopologyTimeMs();
        rep.t_param_gen        = m_gfnff->getParamGenTimeMs();
        rep.t_gpu_upload       = m_gpu_upload_time_ms;

        // Gradient. The gradient Matrix is allocated in the device TU; its Eigen allocator
        // is pinned to the same EIGEN_MAX_ALIGN_BYTES as this host build (CMakeLists), so
        // the aligned packet loads in .norm() are safe (a 16-vs-64 alignment ABI mismatch
        // used to crash here on heap addresses that were only 16-byte aligned).
        if (gradient && m_gpu_workspace) {
            const Matrix& gpu_grad = m_gpu_workspace->gradient();
            if (gpu_grad.size() > 0) {
                rep.gradient_norm = gpu_grad.norm();
                rep.t_gradient.gpu = kt.coulomb;  // gradient mostly computed in Phase 2
            }
        }

        printGFNFFEnergyReport(rep);
    }

    // Cache gradient immediately — use raw memcpy into pre-allocated buffer
    // to avoid Eigen heap allocation (CUDA corrupts adjacent heap metadata)
    if (gradient && m_gpu_workspace) {
        const Matrix& gpu_grad = m_gpu_workspace->gradient();
        if (gpu_grad.rows() == m_cached_gradient.rows() &&
            gpu_grad.cols() == m_cached_gradient.cols()) {
            std::memcpy(m_cached_gradient.data(), gpu_grad.data(),
                        gpu_grad.rows() * gpu_grad.cols() * sizeof(double));
        }
    }

    return m_last_energy;

    } catch (const std::exception& e) {
        CurcumaLogger::error(fmt::format(
            "GFNFFGPU calculateEnergy EXCEPTION: {} — rethrowing", e.what()));
        throw;
    } catch (...) {
        CurcumaLogger::error("GFNFFGPU calculateEnergy: unknown exception — rethrowing");
        throw;
    }
}

// ---------------------------------------------------------------------------
// Remaining ComputationalMethod interface — all delegate to m_gfnff
// ---------------------------------------------------------------------------

// Claude Generated (Oct 2026): host build + upload of the local-inverse preconditioner (see
// LocalInversePreconditioner in eeq_solver.h). The block entries are recomputed from the device
// alpha/gam with the formula of k_eeq_build_matrix; a preconditioner only has to be close to A^-1,
// so rounding differences against the device matrix do not matter (the CG still converges to the
// same tolerance on the true device matrix).
template <class Backend>
void GFNFFGpuMethodImpl<Backend>::prepareLocalInverseGPU(int N, double cutoff_sq)
{
    if constexpr (!Backend::has_localinv_upload) {
        (void)N; (void)cutoff_sq;   // backend without an upload entry point
    } else {
        if (!m_eeq_localinv) return;
        const int nfrag = m_eeq_nfrag;
        // Oct 2026: the default is the device build (CUDA, and ROCm since Oct 2026); the host build
        // stays for the CUDA multi-GPU split and as comparison (CURCUMA_LOCALINV_HOST_BUILD=1).
        const bool device_build = localInverseDeviceBuildSelected(N);
        if (!device_build && nfrag > kLocalInvMaxNfrag) {
            if (!m_eeq_li) {
                CurcumaLogger::warn(fmt::format("EEQ GPU ppcg localinv: nfrag = {} > {} (dense S^-1 on the device), "
                                                "using diagonal Jacobi", nfrag, kLocalInvMaxNfrag));
                m_eeq_li = std::make_shared<LocalInversePreconditioner>();
                m_eeq_li->natoms = N;
                m_eeq_li->nfrag = nfrag;
            }
            return;
        }
        const GeoGradMatrix& g = m_gfnff->getGeometryBohr();
        Matrix X(N, 3);
        for (int i = 0; i < N; ++i)
            for (int d = 0; d < 3; ++d) X(i, d) = g(i, d);
        std::vector<int> frag(N, 0);
        for (int i = 0; i < N && i < static_cast<int>(m_eeq_fraglist.size()); ++i) frag[i] = m_eeq_fraglist[i] - 1;

        // A device-side breakdown dropped the uploaded M: treat it like a rejected build.
        if (m_eeq_li && m_eeq_li->ok && !m_eeq_gpu->localInversePreconditionerActive()) {
            m_eeq_li->ok = false;
            m_eeq_li->why = "indefinite in the GPU CG";
        }
        bool rebuild = !m_eeq_li || m_eeq_li->natoms != N || m_eeq_li->nfrag != nfrag || m_eeq_li->frag != frag
                    || m_eeq_localinv_rebuild <= 0.0;
        if (!rebuild) {
            const double max_dr = (X - m_eeq_li->geometry).rowwise().norm().maxCoeff();
            rebuild = !(max_dr <= m_eeq_localinv_rebuild);
        }
        if (!rebuild) return;

        if constexpr (Backend::has_localinv_device_build) {
            if (device_build) {
                // The device builds G inside the next WP7-E solve, from the matrix assembled there; this
                // record keeps the geometry for the rebuild rule, finishLocalInverseGPU books the result.
                auto rec = std::make_shared<LocalInversePreconditioner>();
                rec->natoms = N;
                rec->nfrag = nfrag;
                rec->frag = frag;
                rec->geometry = X;
                rec->ok = true;
                m_eeq_li = rec;
                m_eeq_gpu->requestLocalInverseDeviceBuild(m_eeq_localinv_radius, m_eeq_localinv_ext);
                return;
            }
        }

        std::vector<double> alpha(N), gam(N);
        Backend::downloadDoubles(alpha.data(), m_gpu_workspace->getDeviceAlphaPtr(), N);
        Backend::downloadDoubles(gam.data(), m_gpu_workspace->getDeviceGamPtr(), N);
        constexpr double tsqrt2pi = 0.797884560802866;   // as in k_eeq_build_matrix
        auto fill = [&](const std::vector<int>& idx, double* out) {
            const int n = static_cast<int>(idx.size());
            for (int c = 0; c < n; ++c) {
                const int j = idx[c];
                for (int rr = c; rr < n; ++rr) {
                    const int i = idx[rr];
                    double v;
                    if (i == j) v = gam[i] + tsqrt2pi / std::sqrt(alpha[i]);
                    else {
                        const double dx = X(i, 0) - X(j, 0), dy = X(i, 1) - X(j, 1), dz = X(i, 2) - X(j, 2);
                        const double r2 = dx * dx + dy * dy + dz * dz;
                        if (cutoff_sq > 0.0 && r2 > cutoff_sq) v = 0.0;
                        else { const double r = std::sqrt(r2); v = std::erf(r / std::sqrt(alpha[i] + alpha[j])) / r; }
                    }
                    out[static_cast<size_t>(c) * n + rr] = v;
                }
            }
        };
        m_eeq_li = EEQSolver::buildLocalInverse(X, frag, nfrag, fill, m_eeq_localinv_radius, m_eeq_localinv_ext,
                                                nullptr, m_eeq_localinv_threads, /*want_dense_S_inverse=*/true);
        m_eeq_li_build_ms += m_eeq_li->build_ms;
        bool uploaded = false;
        double upload_ms = 0.0;
        if (m_eeq_li->ok) {
            const auto tu = std::chrono::steady_clock::now();
            auto& G = m_eeq_li->G;
            auto& H = m_eeq_li->H;
            G.makeCompressed();
            H.makeCompressed();
            auto csr = [](const Eigen::SparseMatrix<double, Eigen::RowMajor>& M, std::vector<int>& rp,
                          std::vector<int>& col, std::vector<double>& val) {
                rp.assign(M.outerIndexPtr(), M.outerIndexPtr() + M.rows() + 1);
                col.assign(M.innerIndexPtr(), M.innerIndexPtr() + M.nonZeros());
                val.assign(M.valuePtr(), M.valuePtr() + M.nonZeros());
            };
            std::vector<int> grp, gcol, hrp, hcol;
            std::vector<double> gval, hval;
            csr(G, grp, gcol, gval);
            csr(H, hrp, hcol, hval);
            std::vector<double> sinv(m_eeq_li->S_inv.data(), m_eeq_li->S_inv.data() + m_eeq_li->S_inv.size());
            uploaded = m_eeq_gpu->setLocalInversePreconditioner(N, nfrag, grp, gcol, gval, hrp, hcol, hval, sinv);
            upload_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tu).count();
            if (!uploaded) { m_eeq_li->ok = false; m_eeq_li->why = "device upload failed"; }
        }
        if (m_eeq_li->ok) {
            ++m_eeq_li_builds;
        } else {
            m_eeq_gpu->clearLocalInversePreconditioner();
            ++m_eeq_li_failures;
            CurcumaLogger::warn(fmt::format("EEQ GPU ppcg localinv: preconditioner rejected ({}), using diagonal Jacobi "
                                            "until the next rebuild (failures so far: {})", m_eeq_li->why, m_eeq_li_failures));
        }
        static const bool print_iters = [] { const char* e = std::getenv("CURCUMA_PPCG_ITERS"); return e && *e && *e != '0'; }();
        if (print_iters)
            fprintf(stderr, "LOCALINV GPU build ok=%d build_ms=%.1f blocks_ms=%.1f check_ms=%.1f sinv_ms=%.1f upload_ms=%.1f nnz_row=%.1f\n",
                    m_eeq_li->ok ? 1 : 0, m_eeq_li->build_ms, m_eeq_li->blocks_ms, m_eeq_li->check_ms, m_eeq_li->sinv_ms, upload_ms,
                    m_eeq_li->nnz_per_row);
    }
}

template <class Backend>
void GFNFFGpuMethodImpl<Backend>::finishLocalInverseGPU()
{
    if constexpr (Backend::has_localinv_device_build) {
        if (!m_eeq_localinv) return;
        const auto st = m_eeq_gpu->localInverseDeviceStats();
        if (st.builds > m_eeq_li_seen_dev_builds) {
            m_eeq_li_builds += st.builds - m_eeq_li_seen_dev_builds;
            m_eeq_li_build_ms += st.last_build_ms;
        }
        if (st.build_failures > m_eeq_li_seen_dev_failures) {
            // A rejected build runs Jacobi until the next rebuild (as on the host): exact charges, so
            // it is a warning, not a GPU fallback.
            m_eeq_li_failures += st.build_failures - m_eeq_li_seen_dev_failures;
            if (m_eeq_li) { m_eeq_li->ok = false; m_eeq_li->why = st.last_why; }
            CurcumaLogger::warn(fmt::format("EEQ GPU ppcg localinv: preconditioner rejected ({}), using diagonal Jacobi "
                                            "until the next rebuild (failures so far: {})", st.last_why, m_eeq_li_failures));
        }
        for (long k = m_eeq_li_seen_runtime; k < st.runtime_fallbacks; ++k)
            curcuma::reportGpuFallback("GFN-FF EEQ: localinv re-solve with Jacobi", st.last_why);
        m_eeq_li_seen_dev_builds = st.builds;
        m_eeq_li_seen_dev_failures = st.build_failures;
        m_eeq_li_seen_runtime = st.runtime_fallbacks;
    } else if constexpr (Backend::has_localinv_upload) {
        // ROCm (Claude Generated, Oct 2026): host build only, so builds and rejections are booked in
        // prepareLocalInverseGPU; here only the Jacobi re-solves, as GPU fallbacks like on CUDA.
        if (!m_eeq_localinv) return;
        const long rt = m_eeq_gpu->localInverseRuntimeFallbacks();
        for (long k = m_eeq_li_seen_runtime; k < rt; ++k)
            curcuma::reportGpuFallback("GFN-FF EEQ: localinv re-solve with Jacobi",
                                       "breakdown (rz <= 0 / pAp <= 0) or stall");
        m_eeq_li_seen_runtime = rt;
    }
}

template <class Backend>
bool GFNFFGpuMethodImpl<Backend>::updateGeometry(const Matrix& geometry)
{
    bindDevice();
    if (!m_gfnff) return false;
    return m_gfnff->UpdateMolecule(geometry);
}

template <class Backend>
Matrix GFNFFGpuMethodImpl<Backend>::getGradient() const
{
    // Return cached gradient (copied immediately after calculate() to avoid
    // heap corruption issues when copying from GPU workspace later).
    //
    // Claude Generated (Sep 2026): the ComputationalMethod::getGradient() contract is
    // Eh/Angstrom (Known Issue #28). The CPU GFN-FF wrapper got the Bohr->Angstrom
    // conversion in that fix, this GPU wrapper (CUDA + ROCm) did not: the workspace works in
    // Bohr like the Fortran reference, so GPU MD forces were 1/au = 1.89x too small and every
    // gfnff_gpu_vs_cpu ctest reported a gradient norm ratio of 0.529. m_cached_gradient stays
    // in Eh/Bohr (internal); only the two public accessors convert.
    return m_cached_gradient / au;
}

template <class Backend>
void GFNFFGpuMethodImpl<Backend>::copyGradientTo(Matrix& target) const
{
    // Claude Generated (March 2026): Direct copy into pre-allocated target.
    // Avoids temporary Matrix creation (heap alloc) that crashes on CUDA-corrupted heap.
    // Sep 2026: converted to Eh/Angstrom in place, same contract as getGradient().
    if (m_cached_gradient.size() > 0 && target.rows() == m_cached_gradient.rows()
        && target.cols() == m_cached_gradient.cols()) {
        const double* src = m_cached_gradient.data();
        double* dst = target.data();
        const Eigen::Index n = m_cached_gradient.size();
        for (Eigen::Index k = 0; k < n; ++k)
            dst[k] = src[k] / au;
    } else {
        target = m_cached_gradient / au;
    }
}

template <class Backend>
Vector GFNFFGpuMethodImpl<Backend>::getCharges() const
{
    return m_gfnff->Charges();
}

template <class Backend>
Vector GFNFFGpuMethodImpl<Backend>::getBondOrders() const
{
    return m_gfnff->BondOrders();
}

template <class Backend>
Position GFNFFGpuMethodImpl<Backend>::getDipole() const
{
    return Position{0.0, 0.0, 0.0};
}

template <class Backend>
void GFNFFGpuMethodImpl<Backend>::setThreadCount(int threads)
{
    (void)threads;  // GPU path does not use CPU threads for bonded terms
}

template <class Backend>
void GFNFFGpuMethodImpl<Backend>::setParameters(const json& params)
{
    m_parameters = params;
    if (m_gfnff)
        m_gfnff->setParameters(params);
}

template <class Backend>
json GFNFFGpuMethodImpl<Backend>::getParameters() const
{
    return m_parameters;
}

template <class Backend>
bool GFNFFGpuMethodImpl<Backend>::hasError() const
{
    return m_has_error;
}

template <class Backend>
void GFNFFGpuMethodImpl<Backend>::clearError()
{
    m_has_error = false;
    m_error_message.clear();
}

template <class Backend>
std::string GFNFFGpuMethodImpl<Backend>::getErrorMessage() const
{
    return m_error_message;
}

template <class Backend>
json GFNFFGpuMethodImpl<Backend>::getEnergyDecomposition() const
{
    json energy_json;
    energy_json["GPU_Total"]  = m_last_energy;
    return energy_json;
}

template <class Backend>
const Vector& GFNFFGpuMethodImpl<Backend>::getGPUdEdcn() const
{
    return m_gpu_workspace->dEdcnTotal();
}

// WP-P1 (May 2026): per-phase CPU timings (CN/EEQ on host) for diagnostics
template <class Backend>
json GFNFFGpuMethodImpl<Backend>::getLastPrepTiming() const
{
    if (!m_gfnff) return {};
    // Claude Generated (Sep 2026): the GPU path calls GFNFF::prepareCNAndEEQ ONLY with
    // skip_eeq=true (CN passed in externally too), so m_gfnff's own PrepTiming has real
    // numbers for "eeq_topo" (topology-input caching) and "cnf" (per-atom CNF fill, still
    // done CPU-side for the gradient) but "cn"/"eeq_solve"/"dcn"/"d4_gw" are structurally
    // always 0 there — those steps happen on the GPU (or, for eeq_solve, sometimes on the
    // CPU fallback path) and are measured in m_last_gpu_prep_timing instead. Report the
    // real number for each field regardless of where it was measured.
    const auto& t = m_gfnff->getLastPrepTiming();
    const auto& g = m_last_gpu_prep_timing;
    return {
        {"cn",          g.cn},              // GPU CN kernel (host wall-clock around the launch)
        {"eeq_topo",    t.eeq_topo},         // CPU: topology-input cache (real)
        {"cnf",         t.cnf},              // CPU: per-atom CNF fill (real, gradient only)
        {"dcn",         t.dcn},              // 0 on this path: CN derivatives live GPU-side (WP5-B)
        {"d4_gw",       t.d4_gw},            // 0 on this path: Gaussian weights computed on GPU
        {"eeq_solve",   g.eeq_solve},        // EEQ Phase 2 solve, GPU Schur/PCG or CPU fallback
        {"charge_dist", t.charge_dist},      // 0 on this path: charges set directly on GPU
        {"total",       g.total},            // full calculateEnergy() wall time (GPU step == "prep")
        // GPU-only sub-phases with no CPU-path analogue (additive; harmless for any
        // consumer that only reads the fields above).
        {"gpu_dlogdcn",             g.dlogdcn},
        {"gpu_hb_cn",               g.hb_cn},
        {"gpu_hbxb_update",         g.hbxb},
        {"gpu_pre",                 g.pre},
        {"gpu_hb_detect",           g.hb_detect},
        {"gpu_hb_rebuild",          g.hb_rebuild},
        {"gpu_hb_upload",           g.hb_upload},
        {"gpu_rep_update",          g.rep_update},
        {"gpu_rep_upload",          g.rep_upload},
        {"gpu_rep_pairs",           g.rep_pairs},
        {"gpu_disp_update",         g.disp_update},
        {"gpu_coul_update",         g.coul_update},
        {"gpu_kernel_launch",       g.launch},
        {"gpu_coulomb_postprocess", g.coulomb_postprocess},
    };
}

template <class Backend>
void GFNFFGpuMethodImpl<Backend>::setForcePhaseTiming(bool on)
{
    if (m_gfnff) m_gfnff->setForcePhaseTiming(on);
    m_force_phase_timing = on;
    if (m_gpu_workspace) m_gpu_workspace->setRecordKernelTimings(on);
}

// WP-P1 (May 2026): per-stream / per-kernel-category GPU timings (FFTermTimings)
template <class Backend>
json GFNFFGpuMethodImpl<Backend>::getStreamTimings() const
{
    if (!m_gpu_workspace) return {};
    const auto& kt = m_gpu_workspace->kernelTimings();
    json out;
    auto add = [&](const char* name, double v) { if (v >= 0.0) out[name] = v; };
    add("bonds",          kt.bonds);
    add("angles",         kt.angles);
    add("dihedrals",      kt.dihedrals);
    add("inversions",     kt.inversions);
    add("stors",          kt.stors);
    add("dispersion",     kt.dispersion);
    add("bonded_rep",     kt.bonded_rep);
    add("nonbonded_rep",  kt.nonbonded_rep);
    add("coulomb",        kt.coulomb);
    add("hbond",          kt.hbond);
    add("xbond",          kt.xbond);
    add("atm",            kt.atm);
    add("batm",           kt.batm);
    return out;
}

template <class Backend>
void GFNFFGpuMethodImpl<Backend>::verifyDeviceHBList(const std::vector<int>& di, const std::vector<int>& dj,
                                                     const std::vector<int>& dk, const std::vector<int>& dcase,
                                                     const std::vector<int>& dparent)
{
    using Entry = std::array<int, 5>;   // A, H, B, case, parent
    auto canon = [](int A, int H, int B, int ct, int par) {
        if (ct == 1 && A > B) std::swap(A, B);   // case 1: the host orientation follows its cell order
        return Entry{A, H, B, ct, par};
    };
    std::vector<Entry> dev, host;
    dev.reserve(di.size());
    for (size_t m = 0; m < di.size(); ++m) dev.push_back(canon(di[m], dj[m], dk[m], dcase[m], dparent[m]));
    const auto hb = m_gfnff->detectHBHostReference();
    host.reserve(hb.size());
    for (const auto& x : hb) host.push_back(canon(x.i, x.j, x.k, x.case_type, x.acceptor_parent_index));
    std::sort(dev.begin(), dev.end());
    std::sort(host.begin(), host.end());
    std::vector<Entry> only_dev, only_host;
    std::set_difference(dev.begin(), dev.end(), host.begin(), host.end(), std::back_inserter(only_dev));
    std::set_difference(host.begin(), host.end(), dev.begin(), dev.end(), std::back_inserter(only_host));
    long long dup_dev = 0;
    for (size_t m = 1; m < dev.size(); ++m) if (dev[m] == dev[m - 1]) ++dup_dev;
    fmt::print(stderr, "HBVERIFY device {} host {} only-device {} only-host {} duplicates-device {}\n",
               dev.size(), host.size(), only_dev.size(), only_host.size(), dup_dev);
    auto show = [](const char* tag, const std::vector<Entry>& v) {
        for (size_t m = 0; m < v.size() && m < 5; ++m)
            fmt::print(stderr, "HBVERIFY   {} A={} H={} B={} case={} parent={}\n", tag, v[m][0], v[m][1], v[m][2], v[m][3], v[m][4]);
    };
    show("only-device", only_dev);
    show("only-host", only_host);
    if (!only_dev.empty() || !only_host.empty())
        CurcumaLogger::warn(fmt::format("GPU H-bond detection differs from the host list: {} only on the device, {} only on the host",
                                         only_dev.size(), only_host.size()));
}

template <class Backend>
void GFNFFGpuMethodImpl<Backend>::buildDispersionPairsOnDevice()
{
    // WP-A (Jun 2026) device build, moved here unchanged from initGPUWorkspace() (Sep 2026).
    D4ParameterGenerator* d4 = m_gfnff->getD4Generator();
    if (!d4) return;
    // Oct 2026: the setup-time call comes from initGPUWorkspace() before the first setGeometry(): the device
    // coordinates were zero, every pair counted as "within 60 Bohr" and the list held all N(N-1)/2 pairs
    // (polymer_2x: 26.8 instead of 9.4 M; energies unaffected since k_dispersion skips pairs beyond r_cut,
    // but dc6dcn and dispersion ran 1.9x / 1.6x slower until the first rebuild). Upload the geometry first.
    m_gpu_workspace->setGeometry(m_gfnff->getGeometryBohr());
    const int Nd = static_cast<int>(m_atom_types.size());
    // Host Gaussian weights (same as getChargeWeightedC6 uses) -> flat [N*MAX_REF].
    const auto& hgw = d4->getGaussianWeights();
    std::vector<double> gw_flat(static_cast<size_t>(Nd) * D4ParameterGenerator::MAX_REF, 0.0);
    for (int a = 0; a < Nd && a < static_cast<int>(hgw.size()); ++a) {
        int nref = std::min<int>(static_cast<int>(hgw[a].size()), D4ParameterGenerator::MAX_REF);
        for (int r = 0; r < nref; ++r) gw_flat[static_cast<size_t>(a) * D4ParameterGenerator::MAX_REF + r] = hgw[a][r];
    }
    std::vector<double> sqrtzr4r2(118, 0.0);
    for (int z = 1; z <= 118; ++z) sqrtzr4r2[z - 1] = d4->getSqrtZr4r2(z);
    const Vector& tc = m_gfnff->getTopologyInfo().topology_charges;
    std::vector<double> topo_q(tc.data(), tc.data() + tc.size());
    // GFN-FF D4: a1=0.58, a2=4.80, 60 Bohr cutoff (matches the host generator).
    m_gpu_workspace->generateDispersionPairListOnGPU(gw_flat, sqrtzr4r2, topo_q, 0.58, 4.80,
                                                     D4ParameterGenerator::PAIR_BUILD_CUTOFF_BOHR);
}

#endif // USE_CUDA || USE_ROCM
