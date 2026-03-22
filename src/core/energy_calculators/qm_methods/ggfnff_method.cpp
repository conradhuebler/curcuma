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
    // CUDA allocations corrupt adjacent heap metadata. To prevent the
    // FFWorkspace destructor from crashing on corrupted free-list entries,
    // disconnect the CPU residual BEFORE destroying the GPU workspace,
    // then leak the CPU residual (same pattern as m_gpu_params_leaked).
    if (m_gpu_workspace) {
        m_gpu_workspace->setCPUResidualWorkspace(nullptr);
    }
    m_gpu_workspace.reset();

    // Leak the CPU residual — its heap metadata may be corrupted by CUDA.
    // Cost: proportional to HB/XB/ATM/BATM list sizes (typically <50 KB).
    if (m_cpu_residual) {
        m_cpu_residual.release();
    }
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
        m_gpu_workspace = std::make_unique<FFWorkspaceGPU>(*pending, natoms, atom_types);

        // === 2. Create CPU residual workspace for HB/XB/ATM/BATM/sTors ===
        // These terms remain on CPU because they are complex/sparse.
        // GPU workspace's calculate() automatically calls m_cpu_residual->calculate().
        {
            GFNFFParameterSet residual_params;
            // Copy only CPU-residual terms
            residual_params.hbonds = pending->hbonds;
            residual_params.xbonds = pending->xbonds;
            residual_params.atm_triples = pending->atm_triples;
            residual_params.batm_triples = pending->batm_triples;
            residual_params.storsions = pending->storsions;
            residual_params.bond_hb_data = pending->bond_hb_data;
            // Charges needed for BATM and HB/XB energy calculations
            residual_params.eeq_charges = pending->eeq_charges;
            residual_params.topology_charges = pending->topology_charges;
            residual_params.method_type = pending->method_type;

            m_cpu_residual = std::make_unique<FFWorkspace>(1);  // single-threaded
            m_cpu_residual->setAtomTypes(atom_types);
            m_cpu_residual->setInteractionLists(std::move(residual_params));
            m_cpu_residual->setTopologyCharges(pending->topology_charges);
            m_cpu_residual->partition();
        }

        // === 3. Connect CPU residual to GPU workspace ===
        // FFWorkspaceGPU::calculate() will automatically add CPU residual energy/gradient.
        // State setters (setGeometry, setD3CN, setEEQCharges, etc.) forward to m_cpu_residual.
        m_gpu_workspace->setCPUResidualWorkspace(m_cpu_residual.get());

        // WORKAROUND: Leak the parameter set — FFWorkspaceGPU's CUDA allocations corrupt
        // adjacent heap metadata, making the GFNFFParameterSet unfreeable.
        // Cost: ~100 KB one-time.  TODO: investigate CUDA root cause.
        m_gpu_params_leaked = pending.release();

        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::success(fmt::format(
                "ggfnff: GPU workspace ready ({} atoms, {} bonds, {} disp pairs, "
                "{} HB, {} XB on CPU residual)",
                natoms, m_gpu_workspace->bondCount(), m_gpu_workspace->dispersionCount(),
                m_cpu_residual->getHBondCount(), m_cpu_residual->getXBondCount()));
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

    // === Step 1: CN + EEQ on CPU (via GFNFF helper) ===
    m_gfnff->prepareCNAndEEQ(gradient);
    auto t1 = std::chrono::high_resolution_clock::now();

    // === Step 2: Distribute state to GPU workspace ===
    // FFWorkspaceGPU setters automatically forward to m_cpu_residual.
    const Matrix& geom_bohr = m_gfnff->getGeometryBohr();
    const Vector& cn = m_gfnff->getLastCN();
    const Vector& charges = m_gfnff->getLastCharges();

    m_gpu_workspace->setGeometry(geom_bohr);
    m_gpu_workspace->setD3CN(cn);
    m_gpu_workspace->setEEQCharges(charges);

    if (gradient) {
        const Vector& cnf = m_gfnff->getLastCNF();
        const auto& dcn = m_gfnff->getLastCNDerivatives();
        m_gpu_workspace->setCNDerivatives(cn, cnf, dcn);

        const Matrix* dc6dcn = m_gfnff->getDC6DCNPtr();
        if (dc6dcn) {
            m_gpu_workspace->setDC6DCNPtr(dc6dcn);
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

        std::unordered_map<int, double> hb_cn_map;
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
                hb_cn_map[H] += 0.5 * (1.0 + std::erf(arg));
            }
        }

        const auto& bonds = m_gpu_params_leaked->bonds;
        const int nb = static_cast<int>(bonds.size());
        std::vector<double> hb_cn_values(nb, 0.0);
        for (int b = 0; b < nb; ++b) {
            if (bonds[b].nr_hb < 1) continue;
            int H = -1;
            if (m_atom_types[bonds[b].i] == 1) H = bonds[b].i;
            else if (m_atom_types[bonds[b].j] == 1) H = bonds[b].j;
            if (H >= 0) {
                auto it = hb_cn_map.find(H);
                hb_cn_values[b] = (it != hb_cn_map.end()) ? it->second : 0.0;
            }
        }
        m_gpu_workspace->updateBondHBCN(hb_cn_values);
    }

    auto t2 = std::chrono::high_resolution_clock::now();

    // === Step 3: Dynamic HB/XB re-detection ===
    // Updates m_cpu_residual (and internal GFNFF m_workspace/m_forcefield)
    m_gfnff->updateHBXBIfNeeded(m_cpu_residual.get());
    auto t3 = std::chrono::high_resolution_clock::now();

    // === Step 4: Calculate — GPU terms + CPU residual (HB/XB/ATM/BATM) ===
    m_last_energy = m_gpu_workspace->calculate(gradient);
    auto t4 = std::chrono::high_resolution_clock::now();

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::energy_abs(m_last_energy, "GFN-FF (GPU) Energy");

        double ms_cneeq  = std::chrono::duration<double, std::milli>(t1 - t0).count();
        double ms_upload  = std::chrono::duration<double, std::milli>(t2 - t1).count();
        double ms_hbxb   = std::chrono::duration<double, std::milli>(t3 - t2).count();
        double ms_gpu    = std::chrono::duration<double, std::milli>(t4 - t3).count();
        double ms_total  = std::chrono::duration<double, std::milli>(t4 - t0).count();
        CurcumaLogger::result_fmt("GFN-FF (GPU) Timing: CN+EEQ={:.1f}ms Upload+HBCN={:.1f}ms HB/XB={:.1f}ms GPU+Residual={:.1f}ms Total={:.1f}ms",
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
    // Gradient comes from GPU workspace (includes CPU residual contribution)
    if (m_gpu_workspace) {
        return m_gpu_workspace->gradient();
    }
    return Matrix();
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

#endif // USE_CUDA
