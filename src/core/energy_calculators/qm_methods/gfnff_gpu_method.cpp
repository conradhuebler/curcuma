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
#include "src/core/energy_calculators/ff_methods/gfnff_par.h"
#include "src/core/energy_calculators/ff_methods/cuda/gpu_utils.h"

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

        // Pre-allocate HB-CN staging buffers (before CUDA init corrupts heap)
        m_hb_cn_values.resize(pending->bonds.size(), 0.0);
        m_hb_cn_per_atom.resize(natoms, 0.0);

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

        // Phase 2: Upload C6 reference table for GPU dc6dcn per-pair computation
        D4ParameterGenerator* d4 = m_gfnff->getD4Generator();
        if (d4 && !d4->getC6FlatCache().empty()) {
            m_gpu_workspace->uploadC6ReferenceTable(d4->getC6FlatCache(), d4->getRefN());
            // Phase 6: Upload reference CN values for GPU Gaussian weight kernel
            m_gpu_workspace->uploadRefCN(d4->getRefCN());
        }

        // Claude Generated (March 2026): GPU EEQ solver (cuSOLVER Cholesky)
        // Pre-allocate EEQ buffers before creating the solver (heap safety)
        m_eeq_z1.resize(natoms, 0.0);
        m_eeq_Z2.resize(natoms * 8, 0.0);  // max 8 fragments
        m_eeq_charges_gpu.resize(natoms, 0.0);
        m_eeq_gpu = std::make_unique<EEQSolverGPU>(natoms);

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

    // --- Step 1: GPU CN computation ---
    m_gpu_workspace->computeCN(m_atom_types);

    // Copy CN from pinned buffers into pre-allocated Vectors (no heap allocs)
    std::memcpy(m_gpu_cn_final.data(), m_gpu_workspace->getCNPinnedBuffer(), N * sizeof(double));

    // Pass GPU CN to prepareCNAndEEQ — but DON'T call EEQ yet, first launch GPU kernels
    // prepareCNAndEEQ does: Gaussian weights, EEQ solve, weight derivatives, dc6dcn
    // We need CN uploaded to GPU before launching charge-independent kernels.
    m_gpu_workspace->setD3CN(m_gpu_cn_final);

    // --- Step 2a: Set up gradient state (CN pairs, dlogdcn, dc6dcn) ---
    // These must happen before launching charge-independent kernels because
    // k_dispersion needs dc6dcn and k_bonds needs d_cn.
    if (gradient) {
        if (!m_cn_pairs_generated) {
            generateCNPairList(geom_bohr);
            m_gpu_workspace->setCNPairList(m_cn_pair_i, m_cn_pair_j, m_cn_pair_rcov);
        }

        // Compute dlogdcn (logistic squashing derivative) — needed for CN chain-rule
        // Reference: Fortran gfnff_cn.f90 create_dlogCN
        // dlogdcn[i] = exp(cnmax) / (exp(cnmax) + exp(cn_raw[i]))
        // CRITICAL: Must use cn_raw (erf sum), NOT cn_final (log-transformed).
        // The dlogdcn formula is the derivative of log-squashing w.r.t. raw CN.
        constexpr double cnmax = 4.4;
        const double exp_cnmax = std::exp(cnmax);
        Vector cn_raw(N);
        std::memcpy(cn_raw.data(), m_gpu_workspace->getCNRawPinnedBuffer(), N * sizeof(double));
        Vector dlogdcn(N);
        for (int i = 0; i < N; ++i) {
            dlogdcn[i] = exp_cnmax / (exp_cnmax + std::exp(cn_raw[i]));
        }
        m_gpu_workspace->setDlogDCN(dlogdcn);
    }

    // --- Step 2b: Compute per-bond HB coordination numbers ---
    if (m_gpu_params_leaked && !m_gpu_params_leaked->bond_hb_data.empty()) {
        static const std::vector<double>& rcov_base = GFNFFParameters::covalent_rad_d3;
        constexpr double rcov_43 = 4.0 / 3.0;
        constexpr double kn = 27.5;
        constexpr double rcov_scal = 1.78;
        constexpr double thr = 900.0;

        std::fill(m_hb_cn_per_atom.begin(), m_hb_cn_per_atom.end(), 0.0);
        for (const auto& entry : m_gpu_params_leaked->bond_hb_data) {
            int H = entry.H;
            int ati = m_atom_types[H];
            for (int B : entry.B_atoms) {
                int atj = m_atom_types[B];
                double dx = geom_bohr(B, 0) - geom_bohr(H, 0);
                double dy = geom_bohr(B, 1) - geom_bohr(H, 1);
                double dz = geom_bohr(B, 2) - geom_bohr(H, 2);
                double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > thr) continue;
                double r = std::sqrt(r2);
                double rcovij = rcov_scal * rcov_43 * (rcov_base[ati - 1] + rcov_base[atj - 1]);
                double arg = -kn * (r - rcovij) / rcovij;
                m_hb_cn_per_atom[H] += 0.5 * (1.0 + std::erf(arg));
            }
        }

        const auto& bonds = m_gpu_params_leaked->bonds;
        const int nb = static_cast<int>(bonds.size());
        std::fill(m_hb_cn_values.begin(), m_hb_cn_values.end(), 0.0);
        for (int b = 0; b < nb; ++b) {
            if (bonds[b].nr_hb < 1) continue;
            int H = -1;
            if (m_atom_types[bonds[b].i] == 1) H = bonds[b].i;
            else if (m_atom_types[bonds[b].j] == 1) H = bonds[b].j;
            if (H >= 0) {
                m_hb_cn_values[b] = m_hb_cn_per_atom[H];
            }
        }
        m_gpu_workspace->updateBondHBCN(m_hb_cn_values);
    }

    // --- Step 2c: Dynamic HB/XB re-detection (must happen before kernel launch) ---
    m_gfnff->updateHBXBIfNeeded(nullptr);
    if (m_gfnff->consumeHBXBUpdate()) {
        m_gpu_workspace->updateHBonds(m_gfnff->getLastHBonds(), m_atom_types);
        m_gpu_workspace->updateXBonds(m_gfnff->getLastXBonds(), m_atom_types);
        // Phase 8: HB/XB SoA n-values changed → captured graph is stale
        m_gpu_workspace->invalidateGraph();
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    // === OVERLAP: Launch charge-independent GPU kernels ===
    // These run asynchronously while CPU computes EEQ parameters below.
    m_gpu_workspace->prepareAndLaunchChargeIndependent(gradient);

    // === CPU: CN distribution + EEQ parameter extraction (skip CPU EEQ solve) ===
    // prepareCNAndEEQ with skip_eeq=true: does CN, CNF, dcn setup but NO matrix build/solve.
    // Inside, getCachedTopology() may trigger full topology recalculation if displacement check flagged it.
    m_gfnff->prepareCNAndEEQ(gradient, /*gpu_only=*/true, &m_gpu_cn_final, /*skip_eeq=*/true);

    // Update GPU reference geometry ONLY when a full topology recalculation happened.
    // Must NOT update every step — that would make the displacement check always return 0
    // and prevent topology updates during MD (instability bug).
    if (m_gfnff->consumeFullTopologyUpdate()) {
        m_gpu_workspace->updateReferenceGeometry();
        // Phase 8: all SoAs rebuilt after topology change → captured graph is stale
        m_gpu_workspace->invalidateGraph();
    }

    const Vector& cn = m_gfnff->getLastCN();

    // Upload CN derivatives and compute dc6dcn for dispersion gradient
    if (gradient) {
        const Vector& cnf = m_gfnff->getLastCNF();
        m_gpu_workspace->setCNDerivatives(cn, cnf, {});

        // Phase 6: Gaussian weights + dc6dcn computed entirely on GPU
        // (k_gaussian_weights + k_dc6dcn_per_pair, no CPU gw/dgw needed)
        m_gpu_workspace->computeGaussianWeightsOnGPU();
    }

    // === GPU EEQ: Build Coulomb matrix + Cholesky solve on GPU ===
    // O(N) CPU parameter extraction + O(N²) GPU matrix build + O(N³/6) GPU Cholesky
    GFNFF::EEQGPUParams eeq_params = m_gfnff->prepareEEQParametersForGPU(cn);

    // Ensure coord upload on main stream is complete before EEQ stream reads SoA buffers.
    // Currently implicit (computeCN syncs main stream), but explicit for safety.
    m_gpu_workspace->synchronizeMainStream();
    bool eeq_ok = m_eeq_gpu->solve(
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
        m_eeq_Z2.data());

    // CPU Schur complement: q = z1 - Z2 * S^{-1} * (C*z1 - d)
    // For nfrag=1 (common case): scalar division
    Vector charges(N);
    if (eeq_ok) {
        const int nfrag = eeq_params.nfrag;
        if (nfrag == 1) {
            // Fast path: single fragment → scalar Schur complement
            // S = C · Z2 (1×1), where C is all-ones row
            double S = 0.0, Cz1 = 0.0;
            for (int i = 0; i < N; ++i) {
                S += m_eeq_Z2[i];
                Cz1 += m_eeq_z1[i];
            }
            double lambda = (Cz1 - eeq_params.rhs_constraints[0]) / S;
            for (int i = 0; i < N; ++i) {
                m_eeq_charges_gpu[i] = m_eeq_z1[i] - m_eeq_Z2[i] * lambda;
            }
        } else {
            // Multi-fragment: dense Schur complement (nfrag typically small)
            // S = C * Z2 (nfrag × nfrag)
            Eigen::MatrixXd S_mat(nfrag, nfrag);
            Eigen::VectorXd Cz1_vec(nfrag);
            for (int f = 0; f < nfrag; ++f) {
                Cz1_vec(f) = 0.0;
                for (int i = 0; i < N; ++i) {
                    bool in_frag = eeq_params.fraglist.empty()
                                     ? (f == 0)
                                     : (eeq_params.fraglist[i] == f + 1);
                    if (in_frag) Cz1_vec(f) += m_eeq_z1[i];
                }
                for (int g = 0; g < nfrag; ++g) {
                    S_mat(f, g) = 0.0;
                    const double* Z2_col = m_eeq_Z2.data() + g * N;
                    for (int i = 0; i < N; ++i) {
                        bool in_frag = eeq_params.fraglist.empty()
                                         ? (f == 0)
                                         : (eeq_params.fraglist[i] == f + 1);
                        if (in_frag) S_mat(f, g) += Z2_col[i];
                    }
                }
            }
            Eigen::VectorXd d_vec = Eigen::Map<const Eigen::VectorXd>(
                eeq_params.rhs_constraints.data(), nfrag);
            Eigen::VectorXd lambda_vec = S_mat.partialPivLu().solve(Cz1_vec - d_vec);
            for (int i = 0; i < N; ++i) {
                double correction = 0.0;
                for (int f = 0; f < nfrag; ++f) {
                    correction += m_eeq_Z2[f * N + i] * lambda_vec(f);
                }
                m_eeq_charges_gpu[i] = m_eeq_z1[i] - correction;
            }
        }
        charges = Eigen::Map<const Vector>(m_eeq_charges_gpu.data(), N);

        // Store charges in GFNFF for consistency (raw memcpy, no Eigen heap alloc)
        m_gfnff->storeChargesFromGPU(m_eeq_charges_gpu.data(), N);
    } else {
        // GPU Cholesky failed (not SPD) — fall back to CPU EEQ
        CurcumaLogger::warn("GPU EEQ Cholesky failed, falling back to CPU solver");
        m_gfnff->prepareCNAndEEQ(gradient, /*gpu_only=*/true, &m_gpu_cn_final, /*skip_eeq=*/false);
        charges = m_gfnff->getLastCharges();
    }

    auto t2 = std::chrono::high_resolution_clock::now();

    // === Upload EEQ charges and finish: Coulomb + postprocess + download ===
    m_gpu_workspace->setEEQCharges(charges);
    m_last_energy = m_gpu_workspace->launchChargeDependentAndFinish(gradient);
    auto t4 = std::chrono::high_resolution_clock::now();

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::energy_abs(m_last_energy, "GFN-FF (GPU) Energy");

        double ms_prep    = std::chrono::duration<double, std::milli>(t1 - t0).count();
        double ms_overlap = std::chrono::duration<double, std::milli>(t2 - t1).count();
        double ms_finish  = std::chrono::duration<double, std::milli>(t4 - t2).count();
        double ms_total   = std::chrono::duration<double, std::milli>(t4 - t0).count();
        CurcumaLogger::result_fmt("GFN-FF (GPU) Timing: CN+Prep={:.1f}ms EEQ||GPU={:.1f}ms Coul+Fin={:.1f}ms Total={:.1f}ms",
                                   ms_prep, ms_overlap, ms_finish, ms_total);
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
