/*
 * < GFN-FF Computational Method Wrapper >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "gfnff_method.h"
#include "src/core/citation_registry.h"
#include "src/tools/general.h"

GFNFFComputationalMethod::GFNFFComputationalMethod(const std::string& method_name, const json& config)
    : m_parameters(config)
    , m_method_name(method_name)
{
    m_gfnff = std::make_unique<GFNFF>(config);
}

bool GFNFFComputationalMethod::setMolecule(const Mol& mol) {
    if (!m_gfnff) {
        CurcumaLogger::error("GFNFFComputationalMethod: m_gfnff is nullptr!");
        return false;
    }

    // Cache check is now done inside GFNFF::initializeForceField()
    // where it has access to the ForceField instance

    if (!m_gfnff->InitialiseMolecule(mol)) {
        CurcumaLogger::error("GFNFFComputationalMethod: InitialiseMolecule failed");
        return false;
    }

    // F-Q4 (Claude Generated): refuse a topology built on EEQ placeholder charges
    // instead of letting a wrong Coulomb energy/gradient propagate silently.
    if (m_gfnff->eeqSolveFailed()) {
        m_has_error = true;
        m_error_message = "GFN-FF EEQ solver fell back to placeholder charges during "
                          "topology setup; refusing to return an energy built on wrong charges";
        CurcumaLogger::error(m_error_message);
        return false;
    }

    return true;
}

bool GFNFFComputationalMethod::updateGeometry(const Matrix& geometry) {
    if (!m_gfnff) {
        CurcumaLogger::error("GFNFFComputationalMethod: m_gfnff is nullptr!");
        return false;
    }

    return m_gfnff->UpdateMolecule(geometry);
}

double GFNFFComputationalMethod::calculateEnergy(bool gradient) {
    CitationRegistry::cite("gfnff");
    CitationRegistry::cite("d4", "gfnff");
    CitationRegistry::cite("eeq", "gfnff");
    CitationRegistry::cite("pyykko", "gfnff");
    CitationRegistry::cite("sanderson", "gfnff");
    CitationRegistry::cite("ghosh_islam", "gfnff");
    CitationRegistry::cite("atm", "d3");
    CitationRegistry::cite("bj", "d3");
    CitationRegistry::cite("casimir_polder", "d4");
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== GFNFFComputationalMethod::calculateEnergy() START ===");
    }

    m_last_energy = m_gfnff->Calculation(gradient);

    // F-Q4 (Claude Generated): a per-step EEQ fallback (opt/MD) means the charges — and
    // hence the Coulomb energy/gradient — are wrong. Flag it so EnergyCalculator refuses
    // the result instead of feeding E to an optimiser/integrator.
    if (m_gfnff->eeqSolveFailed()) {
        m_has_error = true;
        m_error_message = "GFN-FF EEQ solver fell back to placeholder charges; the Coulomb "
                          "energy/gradient is unreliable";
        CurcumaLogger::error(m_error_message);
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("GFNFFComputationalMethod::calculateEnergy complete");
        CurcumaLogger::energy_abs(m_last_energy, "GFN-FF Energy");
    }

    return m_last_energy;
}

Matrix GFNFFComputationalMethod::getGradient() const {
    return m_gfnff->Gradient();
}

Vector GFNFFComputationalMethod::getCharges() const {
    return m_gfnff->Charges();
}

Vector GFNFFComputationalMethod::getBondOrders() const {
    return m_gfnff->BondOrders();
}

Position GFNFFComputationalMethod::getDipole() const {
    return Position{0.0, 0.0, 0.0};
}

// Claude Generated 2026: Return GFN-FF coordination numbers from last calculateEnergy() call
Vector GFNFFComputationalMethod::getCN() const {
    return m_gfnff ? m_gfnff->getLastCN() : Vector();
}

void GFNFFComputationalMethod::setThreadCount(int threads) {
    // Claude Generated (WP1, May 2026): forward to GFNFF::setThreadCount, which keeps
    // m_parameters["threads"] and the cached m_threads member in sync. EnergyCalculator
    // calls this after MethodFactory::create, so the post-construction path now works.
    if (m_gfnff) {
        m_gfnff->setThreadCount(threads);
    }
    m_parameters["threads"] = (threads > 0 ? threads : 1);
}

void GFNFFComputationalMethod::setParameters(const json& params) {
    m_parameters = params;
    if (m_gfnff) {
        m_gfnff->setParameters(params);
    }
}

json GFNFFComputationalMethod::getParameters() const {
    return m_parameters;
}

bool GFNFFComputationalMethod::hasError() const {
    return m_has_error;
}

void GFNFFComputationalMethod::clearError() {
    m_has_error = false;
    m_error_message.clear();
}

json GFNFFComputationalMethod::getEnergyDecomposition() const {
    json energy_json;

    if (!m_gfnff) {
        // Return zero JSON if GFNFF not initialized
        energy_json = {
            {"Bond", 0.0},
            {"Angle", 0.0},
            {"Torsion", 0.0},
            {"Inversion", 0.0},
            {"Dispersion", 0.0},
            {"Coulomb", 0.0},
            {"HBond", 0.0},
            {"XBond", 0.0},
            {"ATM", 0.0},
            {"BATM", 0.0}
        };
        return energy_json;
    }

    // Get all energy components from GFNFF
    energy_json["Bond"] = m_gfnff->BondEnergy();
    energy_json["Angle"] = m_gfnff->AngleEnergy();
    energy_json["Torsion"] = m_gfnff->DihedralEnergy();
    energy_json["Inversion"] = m_gfnff->InversionEnergy();
    energy_json["Dispersion"] = m_gfnff->DispersionEnergy();
    energy_json["Coulomb"] = m_gfnff->CoulombEnergy();
    energy_json["HBond"] = m_gfnff->HydrogenBondEnergy();
    energy_json["XBond"] = m_gfnff->HalogenBondEnergy();
    energy_json["ATM"] = m_gfnff->ATMEnergy();
    energy_json["BATM"] = m_gfnff->BatmEnergy();

    return energy_json;
}

std::string GFNFFComputationalMethod::getErrorMessage() const {
    return m_error_message;
}

// WP-S2 (May 2026): per-step diagnostics hooks for MDDiagnosticsWriter
int GFNFFComputationalMethod::getHBCount() const
{
    return m_gfnff ? static_cast<int>(m_gfnff->getLastHBonds().size()) : 0;
}

int GFNFFComputationalMethod::getXBCount() const
{
    return m_gfnff ? static_cast<int>(m_gfnff->getLastXBonds().size()) : 0;
}

// WP-P1 (May 2026): expose the cached PrepTiming as JSON for MDDiagnosticsWriter
json GFNFFComputationalMethod::getLastPrepTiming() const
{
    if (!m_gfnff) return {};
    const auto& t = m_gfnff->getLastPrepTiming();
    return {
        {"cn",          t.cn},
        {"eeq_topo",    t.eeq_topo},
        {"cnf",         t.cnf},
        {"dcn",         t.dcn},
        {"d4_gw",       t.d4_gw},
        {"eeq_solve",   t.eeq_solve},
        {"charge_dist", t.charge_dist},
        {"total",       t.total},
    };
}

void GFNFFComputationalMethod::setForcePhaseTiming(bool on)
{
    if (m_gfnff) m_gfnff->setForcePhaseTiming(on);
}