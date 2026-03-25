/*
 * <GGFNFFComputationalMethod — GPU-accelerated GFN-FF wrapper>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): ComputationalMethod adapter for ggfnff.
 * Clean GPU/CPU separation: GFNFF has no GPU knowledge, this class orchestrates.
 */

#ifdef USE_CUDA

#include "ggfnff_method.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/gfnff_par.h"

#include <chrono>
#include <cmath>
#include <cstring>
#include <unordered_map>

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

GGFNFFComputationalMethod::GGFNFFComputationalMethod(const std::string& method_name,
                                                      const json& config)
    : m_parameters(config)
    , m_method_name(method_name)
{
    m_gfnff = std::make_unique<GFNFF>(config);
}

// ---------------------------------------------------------------------------
// Destructor — controlled teardown to survive CUDA heap corruption
// ---------------------------------------------------------------------------

GGFNFFComputationalMethod::~GGFNFFComputationalMethod()
{
    // Destroy GPU workspace first (holds CUDA resources)
    m_gpu_workspace.reset();
}

// ---------------------------------------------------------------------------
// setMolecule — init topology (CPU), then upload parameters to GPU
// ---------------------------------------------------------------------------

bool GGFNFFComputationalMethod::setMolecule(const Mol& mol)
{
    if (!m_gfnff) {
        m_has_error = true;
        m_error_message = "GGFNFFComputationalMethod: m_gfnff is nullptr";
        CurcumaLogger::error(m_error_message);
        return false;
    }

    // Store atom types for GPU workspace (from mol, before InitialiseMolecule)
    m_atom_types = mol.m_atoms;

    auto t_init_start = std::chrono::high_resolution_clock::now();

    // CPU topology + parameter generation (same as CPU gfnff)
    if (!m_gfnff->InitialiseMolecule(mol)) {
        m_has_error = true;
        m_error_message = "GGFNFFComputationalMethod: InitialiseMolecule failed";
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
        CurcumaLogger::result_fmt("ggfnff init: topology={:.1f}ms GPU upload={:.1f}ms", ms_topo, ms_gpu);
    }

    m_initialized = true;
    return true;
}

// ---------------------------------------------------------------------------
// initGPUWorkspace — split params → GPU workspace + CPU residual workspace
// ---------------------------------------------------------------------------

bool GGFNFFComputationalMethod::initGPUWorkspace()
{
    try {
        const std::vector<int>& atom_types = m_atom_types;
        const int natoms = static_cast<int>(atom_types.size());
        if (natoms == 0) {
            m_has_error = true;
            m_error_message = "GGFNFFComputationalMethod: zero atoms after InitialiseMolecule";
            CurcumaLogger::error(m_error_message);
            return false;
        }

        // Consume pre-generated params from initializeForceField().
        // CRITICAL: Do NOT call generateGFNFFParameterSet() again — causes heap corruption.
        std::unique_ptr<GFNFFParameterSet> pending = m_gfnff->consumeCachedParameterSet();
        if (!pending) {
            m_has_error = true;
            m_error_message = "GGFNFFComputationalMethod: no cached params (was InitialiseMolecule called?)";
            CurcumaLogger::error(m_error_message);
            return false;
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

        // Pre-allocate cached gradient BEFORE any CUDA operations that may corrupt heap
        m_cached_gradient = Matrix::Zero(natoms, 3);

        // Pre-allocate HB-CN staging buffers (reused every step without heap allocs)
        if (m_gpu_params_leaked) {
            m_hb_cn_values.resize(m_gpu_params_leaked->bonds.size(), 0.0);
        }

        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::success(fmt::format(
                "ggfnff: GPU workspace ready ({} atoms, {} bonds, {} disp pairs)",
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

double GGFNFFComputationalMethod::calculateEnergy(bool gradient)
{
    if (!m_initialized || !m_gfnff || !m_gpu_workspace) {
        CurcumaLogger::error("GGFNFFComputationalMethod: not initialized");
        return 0.0;
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== GGFNFFComputationalMethod::calculateEnergy() START (GPU path) ===");
    }

    auto t0 = std::chrono::high_resolution_clock::now();

    // === Step 1: CN + EEQ on CPU (via GFNFF helper, gpu_only=true) ===
    // gpu_only: skip sparse dcn matrices and CPU forcefield/workspace distribution
    m_gfnff->prepareCNAndEEQ(gradient, /*gpu_only=*/true);
    auto t1 = std::chrono::high_resolution_clock::now();

    // === Step 2: Distribute state to GPU workspace ===
    const Matrix& geom_bohr = m_gfnff->getGeometryBohr();
    const Vector& cn = m_gfnff->getLastCN();
    const Vector& charges = m_gfnff->getLastCharges();

    m_gpu_workspace->setGeometry(geom_bohr);
    m_gpu_workspace->setD3CN(cn);
    m_gpu_workspace->setEEQCharges(charges);

    if (gradient) {
        const Vector& cnf = m_gfnff->getLastCNF();

        // GPU CN chain-rule: pair list + dlogdcn (replaces sparse dcn matrices)
        // Generate CN pair list once (first call with gradient).
        // The cutoff (2.5× rcov) is generous enough for short MD runs.
        if (!m_cn_pairs_generated) {
            generateCNPairList(geom_bohr);
            m_gpu_workspace->setCNPairList(m_cn_pair_i, m_cn_pair_j, m_cn_pair_rcov);
        }

        // Compute dlogdcn = dCN_squashed/dCN_raw (logistic squashing derivative)
        constexpr double cnmax = 4.4;
        const double lncnmax = std::log(1.0 + std::exp(cnmax));
        const int N = static_cast<int>(cn.size());
        Vector dlogdcn(N);
        for (int i = 0; i < N; ++i) {
            double expval = std::exp(lncnmax - cn[i]);
            dlogdcn[i] = (expval - 1.0) / expval;
        }
        m_gpu_workspace->setDlogDCN(dlogdcn);
        m_gpu_workspace->setCNDerivatives(cn, cnf, {});

        // Upload per-pair dc6/dcn values to GPU for dispersion dEdcn kernel
        const Matrix* dc6dcn = m_gfnff->getDC6DCNPtr();
        if (dc6dcn) {
            m_gpu_workspace->updateDispersionDC6DCN(*dc6dcn);
        }
    }

    // === Step 2b: Compute per-bond HB coordination numbers ===
    // Reference: ff_workspace_gfnff.cpp:46-97 (computeHBCoordinationNumbers)
    // egbond_hb: Modified exponent for X-H bonds in hydrogen bonding
    if (m_gpu_params_leaked && !m_gpu_params_leaked->bond_hb_data.empty()) {
        static const std::vector<double>& rcov_base = GFNFFParameters::covalent_rad_d3;
        constexpr double rcov_43 = 4.0 / 3.0;
        constexpr double kn = 27.5;
        constexpr double rcov_scal = 1.78;
        constexpr double thr = 900.0;

        // Reuse pre-allocated map and vector (no heap allocs on corrupted heap)
        m_hb_cn_map.clear();
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
                m_hb_cn_map[H] += 0.5 * (1.0 + std::erf(arg));
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
                auto it = m_hb_cn_map.find(H);
                m_hb_cn_values[b] = (it != m_hb_cn_map.end()) ? it->second : 0.0;
            }
        }
        m_gpu_workspace->updateBondHBCN(m_hb_cn_values);
    }

    auto t2 = std::chrono::high_resolution_clock::now();

    // === Step 3: Dynamic HB/XB re-detection ===
    m_gfnff->updateHBXBIfNeeded(nullptr);

    // If HB/XB lists changed, re-upload SoA to GPU workspace
    if (m_gfnff->consumeHBXBUpdate()) {
        m_gpu_workspace->updateHBonds(m_gfnff->getLastHBonds(), m_atom_types);
        m_gpu_workspace->updateXBonds(m_gfnff->getLastXBonds(), m_atom_types);
    }
    auto t3 = std::chrono::high_resolution_clock::now();

    // === Step 4: Calculate — all terms on GPU ===
    m_last_energy = m_gpu_workspace->calculate(gradient);
    auto t4 = std::chrono::high_resolution_clock::now();

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::energy_abs(m_last_energy, "GFN-FF (GPU) Energy");

        double ms_cneeq  = std::chrono::duration<double, std::milli>(t1 - t0).count();
        double ms_upload  = std::chrono::duration<double, std::milli>(t2 - t1).count();
        double ms_hbxb   = std::chrono::duration<double, std::milli>(t3 - t2).count();
        double ms_gpu    = std::chrono::duration<double, std::milli>(t4 - t3).count();
        double ms_total  = std::chrono::duration<double, std::milli>(t4 - t0).count();
        CurcumaLogger::result_fmt("GFN-FF (GPU) Timing: CN+EEQ={:.1f}ms Upload+HBCN={:.1f}ms HB/XB={:.1f}ms GPU={:.1f}ms Total={:.1f}ms",
                                   ms_cneeq, ms_upload, ms_hbxb, ms_gpu, ms_total);
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

bool GGFNFFComputationalMethod::updateGeometry(const Matrix& geometry)
{
    if (!m_gfnff) return false;
    return m_gfnff->UpdateMolecule(geometry);
}

Matrix GGFNFFComputationalMethod::getGradient() const
{
    // Return cached gradient (copied immediately after calculate() to avoid
    // heap corruption issues when copying from GPU workspace later)
    return m_cached_gradient;
}

void GGFNFFComputationalMethod::copyGradientTo(Matrix& target) const
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

Vector GGFNFFComputationalMethod::getCharges() const
{
    return m_gfnff->Charges();
}

Vector GGFNFFComputationalMethod::getBondOrders() const
{
    return m_gfnff->BondOrders();
}

Position GGFNFFComputationalMethod::getDipole() const
{
    return Position{0.0, 0.0, 0.0};
}

void GGFNFFComputationalMethod::setThreadCount(int threads)
{
    (void)threads;  // GPU path does not use CPU threads for bonded terms
}

void GGFNFFComputationalMethod::setParameters(const json& params)
{
    m_parameters = params;
    if (m_gfnff)
        m_gfnff->setParameters(params);
}

json GGFNFFComputationalMethod::getParameters() const
{
    return m_parameters;
}

bool GGFNFFComputationalMethod::hasError() const
{
    return m_has_error;
}

void GGFNFFComputationalMethod::clearError()
{
    m_has_error = false;
    m_error_message.clear();
}

std::string GGFNFFComputationalMethod::getErrorMessage() const
{
    return m_error_message;
}

json GGFNFFComputationalMethod::getEnergyDecomposition() const
{
    json energy_json;
    energy_json["GPU_Total"]  = m_last_energy;
    return energy_json;
}

const Vector& GGFNFFComputationalMethod::getGPUdEdcn() const
{
    return m_gpu_workspace->dEdcnTotal();
}

// ---------------------------------------------------------------------------
// Generate CN pair list for GPU CN chain-rule kernel
// Claude Generated (March 2026): Replaces sparse dcn matrices
// ---------------------------------------------------------------------------

void GGFNFFComputationalMethod::generateCNPairList(const Matrix& geom_bohr)
{
    const int N = static_cast<int>(geom_bohr.rows());
    const auto& rcov_d3 = GFNFFParameters::covalent_rad_d3;
    constexpr double rcov_scale = 4.0 / 3.0;
    constexpr double kn = -7.5;

    // Cutoff: where exp(-kn^2 * dr^2) < 1e-12 → |dr| > 1.2 → r > 2.2 * rcov_sum
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
