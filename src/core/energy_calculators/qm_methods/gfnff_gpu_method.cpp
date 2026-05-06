/*
 * <GFNFFGPUComputationalMethod — GPU-accelerated GFN-FF wrapper>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): ComputationalMethod adapter for gfnff GPU path.
 * Clean GPU/CPU separation: GFNFF has no GPU knowledge, this class orchestrates.
 *
 * Usage:
 *   ./curcuma -sp mol.xyz -method gfnff -gpu cuda
 *   ./curcuma -sp mol.xyz -method gfnff -gpu auto   # GPU if available
 */

#ifdef USE_CUDA

#include "gfnff_gpu_method.h"
#include "src/core/curcuma_logger.h"
#include "src/core/citation_registry.h"
#include "src/core/energy_calculators/ff_methods/gfnff_par.h"
#include "src/core/energy_calculators/ff_methods/cuda/gpu_utils.h"
#include "src/core/energy_calculators/ff_methods/cn_calculator.h"
#include "src/core/energy_calculators/ff_methods/forcefield.h"
#include "src/core/energy_calculators/ff_methods/forcefieldthread.h"

#include <chrono>
#include <cmath>
#include <cstring>

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

GFNFFGPUComputationalMethod::GFNFFGPUComputationalMethod(const std::string& method_name,
                                                      const json& config)
    : m_parameters(config)
    , m_method_name(method_name)
{
    // Claude Generated (April 2026): Read accuracy-related flags for GPU fallback behavior
    // config.value("gfnff", config) falls back to config itself so top-level CLI params work.
    json gfnff_cfg = config.value("gfnff", config);
    m_allow_unconverged_charges = gfnff_cfg.value("allow_unconverged_charges", false);
    m_skip_phase2 = gfnff_cfg.value("skip_phase2", false);
    // eeq_rmsd_threshold: per-atom RMSD (Bohr) below which the Cholesky factor is reused.
    // 0.0 (default) = always refactorize, matches reference CPU behavior exactly.
    // Suggested MD value: 0.05–0.10 Bohr (saves ~12 ms/step when geometry barely changes).
    m_eeq_rmsd_threshold = gfnff_cfg.value("eeq_rmsd_threshold", 0.0);

    // EEQ Coulomb-matrix distance cutoff (Bohr). Default 0 = no cutoff, matches Fortran
    // goed_gfnff (gfnff_engrad.F90:1274-1391) and the CPU EEQSolver default. Non-zero
    // values violate Hellmann-Feynman vs. the full Coulomb energy → MD energy drift.
    m_eeq_distance_cutoff = gfnff_cfg.value("eeq_distance_cutoff", 0.0);

    m_gfnff = std::make_unique<GFNFF>(config);
}

// ---------------------------------------------------------------------------
// Destructor — controlled teardown to survive CUDA heap corruption
// ---------------------------------------------------------------------------

GFNFFGPUComputationalMethod::~GFNFFGPUComputationalMethod()
{
    // Destroy GPU workspace first (holds CUDA resources)
    m_gpu_workspace.reset();

    // Leak GFNFF instance — its Eigen member destructors (m_last_cn, m_charges, etc.)
    // crash on CUDA-corrupted heap metadata.  Cost: ~10 KB, process is ending anyway.
    m_gfnff.release();
}

// ---------------------------------------------------------------------------
// setMolecule — init topology (CPU), then upload parameters to GPU
// ---------------------------------------------------------------------------

bool GFNFFGPUComputationalMethod::setMolecule(const Mol& mol)
{
    if (!m_gfnff) {
        m_has_error = true;
        m_error_message = "GFNFFGPUMethod: m_gfnff is nullptr";
        CurcumaLogger::error(m_error_message);
        return false;
    }

    // Store atom types for GPU workspace (from mol, before InitialiseMolecule)
    m_atom_types = mol.m_atoms;

    auto t_init_start = std::chrono::high_resolution_clock::now();

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
    auto t_init_end = std::chrono::high_resolution_clock::now();

    if (CurcumaLogger::get_verbosity() >= 1) {
        double ms_topo = std::chrono::duration<double, std::milli>(t_after_topo - t_init_start).count();
        double ms_gpu  = std::chrono::duration<double, std::milli>(t_init_end - t_after_topo).count();
        CurcumaLogger::result_fmt("gfnff-gpu init: topology={:.1f}ms GPU upload={:.1f}ms", ms_topo, ms_gpu);
    }

    m_initialized = true;

    CitationRegistry::cite("gfnff");
    CitationRegistry::cite("d4", "gfnff");

    return true;
}

// ---------------------------------------------------------------------------
// initGPUWorkspace — split params → GPU workspace + CPU residual workspace
// ---------------------------------------------------------------------------

bool GFNFFGPUComputationalMethod::initGPUWorkspace()
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
        // After FFWorkspaceGPU construction, CUDA may corrupt adjacent heap metadata,
        // making Eigen Vector resize/assign crash. All per-step buffers must be
        // allocated now and reused via memcpy in the hot path.
        m_gpu_cn_final = Vector::Zero(natoms);
        m_cached_gradient = Matrix::Zero(natoms, 3);
        m_gfnff->preAllocateForGPUPath(natoms);

        // Pre-allocate EEQ staging buffers before CUDA init (heap safety).
        // m_eeq_Z2 is A^{-1} * C^T — size N × nfrag. Use actual nfrag from topology;
        // the "max 8 fragments" hardcode was too small for multi-fragment systems.
        {
            int nfrag_actual = m_gfnff->getTopologyInfo().nfrag;
            m_eeq_z1.resize(natoms, 0.0);
            m_eeq_Z2.resize(natoms * std::max(nfrag_actual, 1), 0.0);
            m_eeq_charges_gpu.resize(natoms, 0.0);
        }

        // === 1. Create GPU workspace from full parameter set ===
        // All energy terms (including HB, XB, ATM, BATM, sTors) run on GPU.
        m_gpu_workspace = std::make_unique<FFWorkspaceGPU>(*pending, natoms, atom_types);

        // WORKAROUND: Leak the parameter set — FFWorkspaceGPU's CUDA allocations corrupt
        // adjacent heap metadata, making the GFNFFParameterSet unfreeable.
        // Cost: ~100 KB one-time.  TODO: investigate CUDA root cause.
        m_gpu_params_leaked = pending.release();

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

        // Phase 2: Upload C6 reference table for GPU dc6dcn per-pair computation
        D4ParameterGenerator* d4 = m_gfnff->getD4Generator();
        if (d4 && !d4->getC6FlatCache().empty()) {
            m_gpu_workspace->uploadC6ReferenceTable(d4->getC6FlatCache(), d4->getRefN());
            // Phase 6: Upload reference CN values for GPU Gaussian weight kernel
            m_gpu_workspace->uploadRefCN(d4->getRefCN());
        }

        // Claude Generated (March 2026): GPU EEQ solver (cuSOLVER Cholesky)
        // EEQ staging buffers pre-allocated above (before CUDA init). Create solver here.
        // Lazy refactorization is RMSD-based: force_refactor flag passed per solve() call.
        m_eeq_gpu = std::make_unique<EEQSolverGPU>(natoms);

        // WP2: Upload topology-constant EEQ parameters to GPU.
        // cn=Zero because chi_corrected_static and cnf are CN-independent.
        {
            Vector dummy_cn = Vector::Zero(natoms);
            GFNFF::EEQGPUParams topo_params = m_gfnff->prepareEEQParametersForGPU(dummy_cn);
            m_gpu_workspace->uploadEEQTopologyParams(topo_params);
            m_eeq_fraglist         = topo_params.fraglist;
            m_eeq_rhs_constraints  = topo_params.rhs_constraints;
            m_eeq_nfrag            = topo_params.nfrag;
        }

        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::success(fmt::format(
                "gfnff-gpu: GPU workspace ready ({} atoms, {} bonds, {} disp pairs, GPU EEQ enabled)",
                natoms, m_gpu_workspace->bondCount(), m_gpu_workspace->dispersionCount()));
        }

        return true;

    } catch (const std::exception& e) {
        m_has_error = true;
        m_error_message = std::string("FFWorkspaceGPU init failed: ") + e.what();
        CurcumaLogger::error(m_error_message);
        return false;
    }
}

// ---------------------------------------------------------------------------
// calculateEnergy — orchestrate CN/EEQ/state/calculate (no delegation to GFNFF::Calculation)
// ---------------------------------------------------------------------------

double GFNFFGPUComputationalMethod::calculateEnergy(bool gradient)
{
    if (!m_initialized || !m_gfnff || !m_gpu_workspace) {
        CurcumaLogger::error("GFNFFGPUMethod: not initialized");
        return 0.0;
    }

    CitationRegistry::cite("gfnff");
    CitationRegistry::cite("d4", "gfnff");

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== GFNFFGPUMethod::calculateEnergy() START (GPU path) ===");
    }

    auto t0 = std::chrono::high_resolution_clock::now();

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
    const Matrix& geom_bohr = m_gfnff->getGeometryBohr();

    // Claude Generated (April 2026): Upload PBC unit cell to GPU constant memory
    if (m_gfnff->hasPBC()) {
        Eigen::Matrix3d cell_bohr = m_gfnff->getUnitCellBohr();
        Eigen::Matrix3d cell_bohr_inv = cell_bohr.inverse();
        m_gpu_workspace->setUnitCell(cell_bohr.data(), cell_bohr_inv.data(), true);
    }

    // --- Step 0: GPU topology displacement check (Claude Generated March 2026) ---
    // Runs on GPU where coords already live. Result fed to GFNFF::needsFullTopologyUpdate()
    // to skip CPU O(N) Eigen matrix subtraction.
    m_gpu_workspace->setGeometry(geom_bohr);
    {
        bool needs_topo_update = m_gpu_workspace->checkDisplacement(0.5);
        m_gfnff->setExternalTopologyDecision(needs_topo_update);
        if (needs_topo_update) {
            // Reference geometry updated after topology recalculation (triggered inside prepareCNAndEEQ)
            // — updateReferenceGeometry() called below after prepareCNAndEEQ.
        }
    }

    // --- Step 1: GPU CN computation (G-P1: truly async, no sync here) ---
    // k_cn_compute + k_dlogdcn + 3 async D2H launches on main stream.
    // d_cn_final / d_dlogdcn are already on GPU and correct after this call.
    // finalizeCNForCPU() is called after Phase 4 launch to sync without sleeping.
    auto t_cn_start = std::chrono::high_resolution_clock::now();
    m_gpu_workspace->computeCN(m_atom_types);
    // G-P1: removed immediate sync + memcpy + setD3CN here.
    // d_cn_final → d_cn D2D copy happens inside prepareAndLaunchChargeIndependent().
    // setD3CN() no longer needed: d_cn is populated GPU-side from d_cn_final.
    auto t_cn_end = std::chrono::high_resolution_clock::now();

    // --- Step 2a: Set up gradient state (CN pairs) ---
    // G-P1: dlogdcn is already on GPU from k_dlogdcn; setDlogDCN() removed (redundant).
    // d_dlogdcn is used directly by k_cn_chainrule via the main stream ordering.
    if (gradient) {
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
    if (m_gpu_params_leaked && !m_gpu_params_leaked->bond_hb_data.empty()) {
        m_gpu_workspace->computeBondHBCN();
    }
    auto t_hb_cn_end = std::chrono::high_resolution_clock::now();

    // --- Step 2c: Dynamic HB/XB re-detection (must happen before kernel launch) ---
    auto t_hbxb_start = std::chrono::high_resolution_clock::now();
    m_gfnff->updateHBXBIfNeeded(nullptr);
    bool hbxb_updated = m_gfnff->consumeHBXBUpdate();
    if (hbxb_updated) {
        m_gpu_workspace->updateHBonds(m_gfnff->getLastHBonds(), m_atom_types);
        m_gpu_workspace->updateXBonds(m_gfnff->getLastXBonds(), m_atom_types);

        // Claude Generated (Apr 2026): Propagate HB re-detection to bond metadata and alpha pairs.
        // Without this, the bond SoA's nr_hb/hb_H_atom and the HB alpha pair list become stale,
        // causing gradient errors at H-bond atoms during optimization/MD.
        // Reference: Fortran gfnff_ini2.f90:1008-1060 (bond_hb_AHB_set0/set1)
        if (m_gpu_params_leaked) {
            auto hb_update = m_gfnff->rebuildBondHBData(
                m_gfnff->getLastHBonds(), m_gpu_params_leaked->bonds);

            // Update HB alpha (H,B) pair list on GPU
            m_gpu_workspace->updateHBAlphaPairs(hb_update.bond_hb_data, m_atom_types);

            // Update per-bond nr_hb and hb_H_atom metadata on GPU
            m_gpu_workspace->updateBondHBMetadata(hb_update.bond_nr_hb, hb_update.bond_hb_H_atom);

            // Update leaked params' bond_hb_data for HB CN computation (Step 2b)
            m_gpu_params_leaked->bond_hb_data = std::move(hb_update.bond_hb_data);

            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format(
                    "  [DEBUG] HB re-detection: updated bond_hb_data ({} entries), "
                    "hb_alpha ({} pairs), bond nr_hb for {} bonds",
                    static_cast<int>(m_gpu_params_leaked->bond_hb_data.size()),
                    [&]() { int n = 0; for (const auto& e : m_gpu_params_leaked->bond_hb_data) n += static_cast<int>(e.B_atoms.size()); return n; }(),
                    static_cast<int>(hb_update.bond_nr_hb.size())));
            }
        }

        // Phase 8: HB/XB SoA n-values changed → captured graph is stale
        m_gpu_workspace->invalidateGraph();
    }
    auto t_hbxb_end = std::chrono::high_resolution_clock::now();

    // HB/XB pair list consistency: CPU vs GPU (verbosity >= 3)
    if (CurcumaLogger::get_verbosity() >= 3) {
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

    auto t_prep_end = std::chrono::high_resolution_clock::now();

    // === OVERLAP: Launch charge-independent GPU kernels ===
    // These run asynchronously while CPU computes EEQ parameters below.
    // k_dispersion in gradient mode is deferred to Phase 2 (needs dc6dcn from
    // computeGaussianWeightsOnGPU, which runs AFTER Phase 1 on d_cn — current step).
    // G-P1: d_cn_final → d_cn D2D copy is enqueued inside prepareAndLaunchChargeIndependent().
    m_gpu_workspace->prepareAndLaunchChargeIndependent(gradient);

    // Gaussian weights + dc6dcn: must run AFTER Phase 1 so d_cn holds the current step's
    // values (Phase 1 D2D-copies d_cn_final → d_cn on main stream before sA fences).
    // Using d_cn guarantees consistency with the dispersion kernel launched in Phase 2.
    if (gradient) {
        m_gpu_workspace->computeGaussianWeightsOnGPU(/*use_cn_final=*/false);
    }
    auto t_launch_end = std::chrono::high_resolution_clock::now();

    // WP5-D (May 2026): In the normal GPU path (nfrag==1, not skip_phase2) no CPU
    // consumer of CN remains after WP5-B+C. finalizeCNForCPU + prepareCNAndEEQ are
    // kept only for fallback branches that still need CPU-side CN/EEQ.
    if (m_skip_phase2 || m_eeq_nfrag != 1) {
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
    }

    // WP5-B (May 2026): setCNDerivatives no-op — CNF lives on GPU as eeq_topo.d_cnf.
    // WP5-C (May 2026): D4 skip-check moved to GPU (k_check_dc6dcn_skip in computeCN).
    // WP5-D (May 2026): finalizeCNForCPU + prepareCNAndEEQ removed from normal path.
    // All CN consumers are now GPU-side; no CPU round-trip needed.

    // Claude Generated (April 2026): EEQ charge selection
    // Three paths, in priority order:
    //   1. skip_phase2: bypass entirely, use CPU Phase 1 topology charges
    //   2. GPU EEQ Cholesky (nfrag==1 only): O(N²) matrix build + O(N³/6) GPU Cholesky
    //      For nfrag>1 the system is block-diagonal; dense Cholesky is O(N³) on the
    //      full matrix — inefficient. TODO: batched-Cholesky per fragment block.
    //   3. CPU Phase 2 EEQ: always valid fallback (already cached from InitialiseMolecule)
    Vector charges(N);
    if (m_skip_phase2) {
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::warn("GPU path: skip_phase2=true — bypassing GPU EEQ, using CPU Phase 1 topology charges");
        }
        m_gfnff->prepareCNAndEEQ(gradient, /*gpu_only=*/true, &m_gpu_cn_final, /*skip_eeq=*/false);
        charges = m_gfnff->getLastCharges();
    } else {
        // === GPU EEQ: Build Coulomb matrix + Cholesky solve on GPU ===
        // WP2: use topology-constant alpha/gam/chi/cnf already on GPU;
        // rhs_atoms built by k_build_eeq_rhs (launched in prepareAndLaunchChargeIndependent).
        // fraglist and nfrag are cached from last topology build (m_eeq_fraglist/m_eeq_nfrag).
        const bool use_device_rhs = m_gpu_workspace->isEEQTopoValid();

        // Guard: Z2 buffer must hold N × nfrag elements
        if (static_cast<int>(m_eeq_Z2.size()) < N * m_eeq_nfrag) {
            m_eeq_Z2.resize(N * m_eeq_nfrag, 0.0);
        }

        // GPU EEQ is only applicable for nfrag==1 (single connected system).
        // For nfrag>1 (multi-fragment box): use CPU Phase 2 charges directly.
        // TODO: replace with batched per-fragment Cholesky for multi-fragment GPU path.
        const bool gpu_eeq_applicable = (m_eeq_nfrag == 1);

        bool used_gpu_eeq = false;
        const double rhs_c0 = m_eeq_rhs_constraints.empty() ? 0.0 : m_eeq_rhs_constraints[0];

        // EEQ distance cutoff (matches CPU EEQSolver behaviour at eeq_solver.cpp:2641-2643).
        // CPU truncates Coulomb-matrix off-diagonals at 30 Bohr only when natoms > 200.
        // Below that threshold both paths must use 0.0 to keep small-molecule results
        // bit-identical to the previous behaviour.
        const double eeq_cutoff_sq = (N > 200 && m_eeq_distance_cutoff > 0.0)
                                       ? m_eeq_distance_cutoff * m_eeq_distance_cutoff
                                       : 0.0;

        if (gpu_eeq_applicable) {
            // RMSD-based EEQ lazy refactorization:
            // Compute per-atom RMSD (Bohr) from geometry at last full Cholesky build.
            // If RMSD < threshold: reuse cached L (dpotrs only, saves ~12 ms/step).
            // If RMSD >= threshold or threshold==0: full refactorization.
            const Matrix& cur_geom = m_gfnff->getGeometryBohr();
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
                // WP5-A: try GPU Schur path first (nfrag==1 only).
                // Eliminates 22 KB D2H (z1+Z2) + CPU O(N) Schur loop + 11 KB H2D.
                // Falls back to WP2 CPU-Schur for nfrag>1 or Cholesky failure.
                eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUSchur(
                    N, m_eeq_nfrag,
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
                    used_gpu_schur = true;  // restore original
                } else {
                    // nfrag>1 or Cholesky failed: fall back to WP2 solve + CPU Schur
                    eeq_ok = m_eeq_gpu->solveWithDeviceRHS(
                        N, m_eeq_nfrag,
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
            } else {
                // Fallback: EEQ topo not yet uploaded to GPU (should not happen after init)
                CurcumaLogger::warn("WP2: EEQ topo not valid — falling back to CPU param upload");
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
                    // WP5-A: charges already on GPU (d_rhs[0..N-1] in EEQ solver).
                    // D2D copy into impl.d_charges — skips H2D upload in Phase-2.
                    m_gpu_workspace->setEEQDeviceCharges(m_eeq_gpu->getDeviceChargesPtr());
                    // Download for storeChargesFromGPU() — EEQ stream already synced.
                    cudaMemcpy(m_eeq_charges_gpu.data(), m_eeq_gpu->getDeviceChargesPtr(),
                               N * sizeof(double), cudaMemcpyDeviceToHost);
                    m_gfnff->storeChargesFromGPU(m_eeq_charges_gpu.data(), N);
                } else {
                    // CPU Schur complement: q = z1 - Z2 * λ  (nfrag>1 or WP2 fallback)
                    double S = 0.0, Cz1 = 0.0;
                    for (int i = 0; i < N; ++i) {
                        S   += m_eeq_Z2[i];
                        Cz1 += m_eeq_z1[i];
                    }
                    double lambda = (Cz1 - rhs_c0) / S;
                    for (int i = 0; i < N; ++i)
                        m_eeq_charges_gpu[i] = m_eeq_z1[i] - m_eeq_Z2[i] * lambda;
                    charges = Eigen::Map<const Vector>(m_eeq_charges_gpu.data(), N);
                    m_gfnff->storeChargesFromGPU(m_eeq_charges_gpu.data(), N);
                    m_gpu_workspace->setEEQCharges(charges);
                }
                used_gpu_eeq = true;
            } else {
                CurcumaLogger::warn("GPU EEQ Cholesky failed (not SPD), falling back to CPU solver");
            }
        } else {
            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::info(fmt::format(
                    "GPU EEQ: nfrag={} — using CPU Phase 2 charges (TODO: batched-Cholesky per fragment)",
                    m_eeq_nfrag));
            }
        }

        if (!used_gpu_eeq) {
            // CPU Phase 2 EEQ fallback (cached from InitialiseMolecule on first step)
            m_gfnff->prepareCNAndEEQ(gradient, /*gpu_only=*/true, &m_gpu_cn_final, /*skip_eeq=*/false);
            charges = m_gfnff->getLastCharges();
            m_gpu_workspace->setEEQCharges(charges);
        }
    }

    auto t_eeq_end = std::chrono::high_resolution_clock::now();

    // === Finish: Coulomb + postprocess + download ===
    // Note: setEEQCharges / setEEQDeviceCharges already called above per path.
    m_last_energy = m_gpu_workspace->launchChargeDependentAndFinish(gradient);
    auto t4 = std::chrono::high_resolution_clock::now();

    ++m_calc_count;

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::energy_abs(m_last_energy, "GFN-FF (GPU) Energy");

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
        CurcumaLogger::result(fmt::format("  Initialization (not above)     {:>10.1f}  {:>5.1f}%", t_init, pct(t_init)));
        CurcumaLogger::result(fmt::format("  Wallclock Total                {:>10.1f}  100.0%", t_wall));

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

    // Energy decomposition at verbosity 2+
    if (CurcumaLogger::get_verbosity() >= 2) {
        const auto& comp = m_gpu_workspace->energyComponents();
        CurcumaLogger::info("GPU Energy Decomposition:");
        CurcumaLogger::param("  Bond", fmt::format("{:.10f} Eh", comp.bond));
        CurcumaLogger::param("  Angle", fmt::format("{:.10f} Eh", comp.angle));
        CurcumaLogger::param("  Dihedral", fmt::format("{:.10f} Eh", comp.dihedral));
        CurcumaLogger::param("  Inversion", fmt::format("{:.10f} Eh", comp.inversion));
        CurcumaLogger::param("  Dispersion", fmt::format("{:.10f} Eh", comp.dispersion));
        CurcumaLogger::param("  BondedRep", fmt::format("{:.10f} Eh", comp.bonded_rep));
        CurcumaLogger::param("  NonbondedRep", fmt::format("{:.10f} Eh", comp.nonbonded_rep));
        CurcumaLogger::param("  Coulomb", fmt::format("{:.10f} Eh", comp.coulomb));
        CurcumaLogger::param("  HBond", fmt::format("{:.10f} Eh", comp.hbond));
        CurcumaLogger::param("  XBond", fmt::format("{:.10f} Eh", comp.xbond));
        CurcumaLogger::param("  ATM", fmt::format("{:.10f} Eh", comp.atm));
        CurcumaLogger::param("  BATM", fmt::format("{:.10f} Eh", comp.batm));
        CurcumaLogger::param("  sTors", fmt::format("{:.10f} Eh", comp.stors));

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
}

// ---------------------------------------------------------------------------
// Remaining ComputationalMethod interface — all delegate to m_gfnff
// ---------------------------------------------------------------------------

bool GFNFFGPUComputationalMethod::updateGeometry(const Matrix& geometry)
{
    if (!m_gfnff) return false;
    return m_gfnff->UpdateMolecule(geometry);
}

Matrix GFNFFGPUComputationalMethod::getGradient() const
{
    // Return cached gradient (copied immediately after calculate() to avoid
    // heap corruption issues when copying from GPU workspace later)
    return m_cached_gradient;
}

void GFNFFGPUComputationalMethod::copyGradientTo(Matrix& target) const
{
    // Claude Generated (March 2026): Direct memcpy into pre-allocated target.
    // Avoids temporary Matrix creation (heap alloc) that crashes on CUDA-corrupted heap.
    if (m_cached_gradient.size() > 0 && target.rows() == m_cached_gradient.rows()
        && target.cols() == m_cached_gradient.cols()) {
        std::memcpy(target.data(), m_cached_gradient.data(),
                     m_cached_gradient.size() * sizeof(double));
    } else {
        target = m_cached_gradient;
    }
}

Vector GFNFFGPUComputationalMethod::getCharges() const
{
    return m_gfnff->Charges();
}

Vector GFNFFGPUComputationalMethod::getBondOrders() const
{
    return m_gfnff->BondOrders();
}

Position GFNFFGPUComputationalMethod::getDipole() const
{
    return Position{0.0, 0.0, 0.0};
}

void GFNFFGPUComputationalMethod::setThreadCount(int threads)
{
    (void)threads;  // GPU path does not use CPU threads for bonded terms
}

void GFNFFGPUComputationalMethod::setParameters(const json& params)
{
    m_parameters = params;
    if (m_gfnff)
        m_gfnff->setParameters(params);
}

json GFNFFGPUComputationalMethod::getParameters() const
{
    return m_parameters;
}

bool GFNFFGPUComputationalMethod::hasError() const
{
    return m_has_error;
}

void GFNFFGPUComputationalMethod::clearError()
{
    m_has_error = false;
    m_error_message.clear();
}

std::string GFNFFGPUComputationalMethod::getErrorMessage() const
{
    return m_error_message;
}

json GFNFFGPUComputationalMethod::getEnergyDecomposition() const
{
    json energy_json;
    energy_json["GPU_Total"]  = m_last_energy;
    return energy_json;
}

const Vector& GFNFFGPUComputationalMethod::getGPUdEdcn() const
{
    return m_gpu_workspace->dEdcnTotal();
}

// ---------------------------------------------------------------------------
// Generate CN pair list for GPU CN chain-rule kernel
// Claude Generated (March 2026): Replaces sparse dcn matrices
// ---------------------------------------------------------------------------

void GFNFFGPUComputationalMethod::generateCNPairList(const Matrix& geom_bohr)
{
    const int N = static_cast<int>(geom_bohr.rows());
    const auto& rcov_d3 = GFNFFParameters::covalent_rad_d3;
    constexpr double rcov_scale = 4.0 / 3.0;

    // Cutoff: where exp(-kn^2 * dr^2) < 1e-12 (kn=-7.5) → |dr| > 1.2 → r > 2.2 * rcov_sum
    constexpr double cutoff_factor = 2.5;

    m_cn_pair_i.clear();
    m_cn_pair_j.clear();
    m_cn_pair_rcov.clear();

    for (int i = 0; i < N; ++i) {
        int zi = m_atom_types[i];
        if (zi < 1 || zi > static_cast<int>(rcov_d3.size())) continue;
        double rcov_i = rcov_d3[zi - 1] * rcov_scale;

        for (int j = i + 1; j < N; ++j) {
            int zj = m_atom_types[j];
            if (zj < 1 || zj > static_cast<int>(rcov_d3.size())) continue;
            double rcov_j = rcov_d3[zj - 1] * rcov_scale;
            double rcov_sum = rcov_i + rcov_j;

            double dx = geom_bohr(i, 0) - geom_bohr(j, 0);
            double dy = geom_bohr(i, 1) - geom_bohr(j, 1);
            double dz = geom_bohr(i, 2) - geom_bohr(j, 2);
            double r2 = dx*dx + dy*dy + dz*dz;
            double rij = std::sqrt(r2);

            if (rij < cutoff_factor * rcov_sum && rij > 1e-10) {
                m_cn_pair_i.push_back(i);
                m_cn_pair_j.push_back(j);
                m_cn_pair_rcov.push_back(rcov_sum);
            }
        }
    }

    m_cn_pairs_generated = true;

    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::info(fmt::format("  CN pair list: {} pairs (N={})", m_cn_pair_i.size(), N));
}

#endif // USE_CUDA
