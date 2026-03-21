/*
 * <GGFNFFComputationalMethod — GPU-accelerated GFN-FF wrapper>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): ComputationalMethod adapter for ggfnff.
 */

#ifdef USE_CUDA

#include "ggfnff_method.h"
#include "src/core/curcuma_logger.h"

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

    // CPU topology + parameter generation (same as CPU gfnff)
    if (!m_gfnff->InitialiseMolecule(mol)) {
        m_has_error = true;
        m_error_message = "GGFNFFComputationalMethod: InitialiseMolecule failed";
        CurcumaLogger::error(m_error_message);
        return false;
    }

    // Upload parameters to GPU
    if (!initGPUWorkspace()) {
        // Fall back gracefully — error already set
        return false;
    }

    m_initialized = true;
    return true;
}

// ---------------------------------------------------------------------------
// initGPUWorkspace — extract GFNFFParameterSet → FFWorkspaceGPU
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

        // Claude Generated (March 2026): Consume pre-generated params from initializeForceField().
        // CRITICAL: Do NOT call generateGFNFFParameterSet() again — causes heap corruption.
        std::unique_ptr<GFNFFParameterSet> pending = m_gfnff->consumePendingGPUParams();
        if (!pending) {
            m_has_error = true;
            m_error_message = "GGFNFFComputationalMethod: no pending GPU params (was InitialiseMolecule called?)";
            CurcumaLogger::error(m_error_message);
            return false;
        }

        // Create GPU workspace — uploads all static SoA data
        m_gpu_workspace = std::make_unique<FFWorkspaceGPU>(*pending, natoms, atom_types);

        // WORKAROUND: Leak the parameter set — FFWorkspaceGPU's CUDA allocations corrupt
        // adjacent heap metadata, making the GFNFFParameterSet unfreeable ("double free or
        // corruption (out)").  Cost: ~100 KB one-time.  TODO: investigate CUDA root cause.
        m_gpu_params_leaked = pending.release();

        // Inject non-owning pointer into GFNFF so Calculation() routes to GPU
        m_gfnff->setGPUWorkspace(m_gpu_workspace.get());

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
// calculateEnergy — delegates entirely to GFNFF::Calculation()
//   (which internally routes through FFWorkspaceGPU when m_gpu_workspace is set)
// ---------------------------------------------------------------------------

double GGFNFFComputationalMethod::calculateEnergy(bool gradient)
{
    if (!m_initialized || !m_gfnff) {
        CurcumaLogger::error("GGFNFFComputationalMethod: not initialized");
        return 0.0;
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== GGFNFFComputationalMethod::calculateEnergy() START (GPU path) ===");
    }

    m_last_energy = m_gfnff->Calculation(gradient);

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::energy_abs(m_last_energy, "GFN-FF (GPU) Energy");
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
    return m_gfnff->Gradient();
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
    // GPU workspace combines all terms into a single scalar;
    // individual components not available per-term from GPU in Phase 1.
    energy_json["Bond"]       = 0.0;
    energy_json["Angle"]      = 0.0;
    energy_json["Torsion"]    = 0.0;
    energy_json["Inversion"]  = 0.0;
    energy_json["Dispersion"] = 0.0;
    energy_json["Coulomb"]    = 0.0;
    energy_json["HBond"]      = m_gfnff ? m_gfnff->HydrogenBondEnergy() : 0.0;
    energy_json["XBond"]      = m_gfnff ? m_gfnff->HalogenBondEnergy()  : 0.0;
    energy_json["ATM"]        = m_gfnff ? m_gfnff->ATMEnergy()           : 0.0;
    energy_json["BATM"]       = m_gfnff ? m_gfnff->BatmEnergy()          : 0.0;
    energy_json["GPU_Total"]  = m_last_energy;
    return energy_json;
}

#endif // USE_CUDA
