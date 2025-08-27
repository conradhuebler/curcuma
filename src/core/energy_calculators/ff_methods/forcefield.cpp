/*
 * < Generic force field class for curcuma . >
 * Copyright (C) 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "forcefieldderivaties.h"
#include "qmdff_par.h"
#include "src/core/curcuma_logger.h"
#include "uff_par.h"

#include "forcefieldfunctions.h"

#include <fmt/core.h>
#include <fmt/format.h>

#include "forcefield.h"
#include "forcefieldthread.h"

ForceField::ForceField(const json& controller)
{
    json parameter = MergeJson(UFFParameterJson, controller);

    m_threadpool = new CxxThreadPool();
    m_threadpool->setProgressBar(CxxThreadPool::ProgressBarType::None);
    m_threads = parameter["threads"];
    m_gradient_type = parameter["gradient"];

    // Auto-detect parameter file based on geometry file
    if (parameter.contains("geometry_file")) {
        std::string geom_file = parameter["geometry_file"];
        m_auto_param_file = generateParameterFileName(geom_file);
    }
}

ForceField::~ForceField()
{
    delete m_threadpool;
}

// Claude Generated: Temporary method for EnergyCalculator compatibility
// TODO: Eventually merge QMInterface and ForceField into unified interface
void ForceField::setMolecule(const Mol& mol)
{
    // Extract basic molecular information
    m_natoms = mol.m_number_atoms;
    m_geometry = mol.m_geometry;

    // Set atom types from atomic numbers
    std::vector<int> atom_types;
    atom_types.reserve(m_natoms);
    for (int i = 0; i < m_natoms; ++i) {
        atom_types.push_back(mol.m_atoms[i]);
    }
    setAtomTypes(atom_types);
}

void ForceField::UpdateGeometry(const Matrix& geometry)
{
    m_geometry = geometry;
}

void ForceField::UpdateGeometry(const double* coord)
{
#pragma message("replace with raw data")
    for (int i = 0; i < m_natoms; ++i) {
        m_geometry(i, 0) = coord[3 * i + 0];
        m_geometry(i, 1) = coord[3 * i + 1];
        m_geometry(i, 2) = coord[3 * i + 2];
    }
}

void ForceField::UpdateGeometry(const std::vector<std::array<double, 3>>& geometry)
{
#pragma message("replace with raw data")
    for (int i = 0; i < m_natoms; ++i) {
        m_geometry(i, 0) = geometry[i][0];
        m_geometry(i, 1) = geometry[i][1];
        m_geometry(i, 2) = geometry[i][2];
    }
}

void ForceField::setParameter(const json& parameters)
{
    static bool in_setParameter = false;
    if (in_setParameter) {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::warn("Recursive setParameter call detected - preventing infinite loop");
        }
        return;
    }
    in_setParameter = true;

    std::string method_name = "unknown";
    if (parameters.contains("method") && !parameters["method"].is_null()) {
        method_name = parameters["method"].get<std::string>();
    }

    // Level 2+: Force field initialization
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Initializing force field parameters");
        CurcumaLogger::param("method", method_name);
    }

    bool loaded_from_cache = false;

    // Try loading cached parameters if caching enabled and same method
    if (m_enable_caching && parameters.contains("method")) {
        std::string method = parameters["method"];
        loaded_from_cache = tryLoadAutoParameters(method);
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("loaded_from_cache", loaded_from_cache ? "true" : "false");
        }
    }

    if (!loaded_from_cache) {
        // Set new parameters (generation or explicit)
        if (parameters.contains("bonds"))
            setBonds(parameters["bonds"]);
        if (parameters.contains("angles"))
            setAngles(parameters["angles"]);
        if (parameters.contains("dihedrals"))
            setDihedrals(parameters["dihedrals"]);
        if (parameters.contains("inversions"))
            setInversions(parameters["inversions"]);
        if (parameters.contains("vdws"))
            setvdWs(parameters["vdws"]);
        m_parameters = parameters;
        m_method = m_parameters["method"];
        if (m_parameters.contains("e0"))
            m_e0 = m_parameters["e0"];

        AutoRanges();

        // Auto-save new parameters (only if caching enabled)
        if (m_enable_caching) {
            autoSaveParameters();
        }
    }

    // Set method type for ForceFieldThread (regardless of cache)
    int method_type = 1; // default UFF
    if (m_method == "qmdff" || m_method == "quff") {
        method_type = 2;
    } else if (m_method == "gfnff") {
        method_type = 3; // GFN-FF
    }

    // Claude Generated: Print parameter summary after setting parameters
    printParameterSummary();

    in_setParameter = false; // Reset the recursive guard
}

void ForceField::setParameterFile(const std::string& file)
{
    if (!loadParametersFromFile(file)) {
        CurcumaLogger::warn(fmt::format("Failed to load parameter file: {}", file));
    }
}

void ForceField::setBonds(const json& bonds)
{
    m_bonds.clear();
    for (int i = 0; i < bonds.size(); ++i) {
        json bond = bonds[i].get<json>();
        Bond b;
        b.type = bond["type"];
        b.i = bond["i"];
        b.j = bond["j"];
        b.k = bond["k"];
        b.distance = bond["distance"];
        b.exponent = bond["exponent"];

        b.r0_ij = bond["r0_ij"];
        b.r0_ik = bond["r0_ik"];

        b.fc = bond["fc"];
        m_bonds.push_back(b);
    }
}

void ForceField::setAngles(const json& angles)
{
    m_angles.clear();
    for (int i = 0; i < angles.size(); ++i) {
        json angle = angles[i].get<json>();
        Angle a;

        a.type = angle["type"];

        a.i = angle["i"];
        a.j = angle["j"];
        a.k = angle["k"];
        a.C0 = angle["C0"];
        a.C1 = angle["C1"];
        a.C2 = angle["C2"];
        a.fc = angle["fc"];
        a.r0_ij = angle["r0_ij"];
        a.r0_ik = angle["r0_ik"];
        a.theta0_ijk = angle["theta0_ijk"];
        m_angles.push_back(a);
    }
}

void ForceField::setDihedrals(const json& dihedrals)
{
    m_dihedrals.clear();
    for (int i = 0; i < dihedrals.size(); ++i) {
        json dihedral = dihedrals[i].get<json>();
        Dihedral d;
        d.type = dihedral["type"];

        d.i = dihedral["i"];
        d.j = dihedral["j"];
        d.k = dihedral["k"];
        d.l = dihedral["l"];
        d.V = dihedral["V"];
        d.n = dihedral["n"];
        d.phi0 = dihedral["phi0"];
        m_dihedrals.push_back(d);
    }
}

void ForceField::setInversions(const json& inversions)
{
    m_inversions.clear();
    for (int i = 0; i < inversions.size(); ++i) {
        json inversion = inversions[i].get<json>();
        Inversion inv;
        inv.type = inversion["type"];

        inv.i = inversion["i"];
        inv.j = inversion["j"];
        inv.k = inversion["k"];
        inv.l = inversion["l"];
        inv.fc = inversion["fc"];
        inv.C0 = inversion["C0"];
        inv.C1 = inversion["C1"];
        inv.C2 = inversion["C2"];

        m_inversions.push_back(inv);
    }
}

void ForceField::setvdWs(const json& vdws)
{
    m_vdWs.clear();
    for (int i = 0; i < vdws.size(); ++i) {
        json vdw = vdws[i].get<json>();
        vdW v;
        v.type = vdw["type"];

        v.i = vdw["i"];
        v.j = vdw["j"];
        v.C_ij = vdw["C_ij"];
        v.r0_ij = vdw["r0_ij"];

        m_vdWs.push_back(v);
    }
}

void ForceField::setESPs(const json& esps)
{
    m_EQs.clear();
    for (int i = 0; i < esps.size(); ++i) {
        json esp = esps[i].get<json>();
        EQ v;
        v.type = esp["type"];

        v.i = esp["i"];
        v.j = esp["j"];
        v.q_i = esp["q_i"];
        v.q_j = esp["q_j"];
        v.epsilon = esp["epsilon"];

        m_EQs.push_back(v);
    }
}

void ForceField::AutoRanges()
{
    // Level 3+: AutoRanges debug info
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Setting up force field calculation ranges");
        CurcumaLogger::param("method", m_method);
        CurcumaLogger::param("bonds", static_cast<int>(m_bonds.size()));
        CurcumaLogger::param("angles", static_cast<int>(m_angles.size()));
        CurcumaLogger::param("threads", m_threads);
    }

    int free_threads = m_threads;
#pragma message("revert")
    int d3 = false; // m_parameters["d3"];
    if (d3) {
        if (free_threads > 1)
            free_threads--;
        D3Thread* thread = new D3Thread(m_threads - 1, free_threads);
        thread->setParamater(m_parameters);
        thread->Initialise(m_atom_types);
        m_threadpool->addThread(thread);
        m_stored_threads.push_back(thread);
    }
    int h4 = false; // m_parameters["h4"];
    if (h4) {
        if (free_threads > 1)
            free_threads--;
        H4Thread* thread = new H4Thread(m_threads - 1, free_threads);
        thread->setParamater(m_parameters);
        thread->Initialise(m_atom_types);

        m_threadpool->addThread(thread);
        m_stored_threads.push_back(thread);
    }
    for (int i = 0; i < free_threads; ++i) {
        ForceFieldThread* thread = new ForceFieldThread(i, free_threads);
        thread->setGeometry(m_geometry, false);
        m_threadpool->addThread(thread);
        m_stored_threads.push_back(thread);
        if (std::find(m_uff_methods.begin(), m_uff_methods.end(), m_method) != m_uff_methods.end()) {
            thread->setMethod(1);
        } else if (std::find(m_qmdff_methods.begin(), m_qmdff_methods.end(), m_method) != m_qmdff_methods.end()) {
            thread->setMethod(2);
        } else if (m_method == "gfnff") {
            thread->setMethod(3); // GFN-FF
        }
        for (int j = int(i * m_bonds.size() / double(free_threads)); j < int((i + 1) * m_bonds.size() / double(free_threads)); ++j) {
            if (m_method == "gfnff") {
                thread->addGFNFFBond(m_bonds[j]);
            } else {
                thread->addBond(m_bonds[j]);
            }
        }

        for (int j = int(i * m_angles.size() / double(free_threads)); j < int((i + 1) * m_angles.size() / double(free_threads)); ++j) {
            if (m_method == "gfnff") {
                thread->addGFNFFAngle(m_angles[j]);
            } else {
                thread->addAngle(m_angles[j]);
            }
        }

        for (int j = int(i * m_dihedrals.size() / double(free_threads)); j < int((i + 1) * m_dihedrals.size() / double(free_threads)); ++j) {
            if (m_method == "gfnff") {
                thread->addGFNFFDihedral(m_dihedrals[j]);
            } else {
                thread->addDihedral(m_dihedrals[j]);
            }
        }

        for (int j = int(i * m_inversions.size() / double(free_threads)); j < int((i + 1) * m_inversions.size() / double(free_threads)); ++j) {
            if (m_method == "gfnff") {
                thread->addGFNFFInversion(m_inversions[j]);
            } else {
                thread->addInversion(m_inversions[j]);
            }
        }

        for (int j = int(i * m_vdWs.size() / double(free_threads)); j < int((i + 1) * m_vdWs.size() / double(free_threads)); ++j) {
            if (m_method == "gfnff") {
                thread->addGFNFFvdW(m_vdWs[j]);
            } else {
                thread->addvdW(m_vdWs[j]);
            }
        }

        for (int j = int(i * m_EQs.size() / double(free_threads)); j < int((i + 1) * m_EQs.size() / double(free_threads)); ++j)
            thread->addEQ(m_EQs[j]);
    }

    // Level 3+: AutoRanges completion info
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("created_threads", free_threads);
        CurcumaLogger::param("total_stored_threads", static_cast<int>(m_stored_threads.size()));
        CurcumaLogger::info("Force field calculation ranges setup completed");
    }
}

Eigen::MatrixXd ForceField::NumGrad()
{
    Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(m_natoms, 3);

    double dx = 1e-6; // m_d;
    // bool g = m_CalculateGradient;
    // m_CalculateGradient = false;
    double E1, E2;
    for (int i = 0; i < m_natoms; ++i) {
        for (int j = 0; j < 3; ++j) {
            m_geometry(i, j) += dx;
            E1 = Calculate(false);
            m_geometry(i, j) -= 2 * dx;
            E2 = Calculate(false);
            gradient(i, j) = (E1 - E2) / (2 * dx);
            m_geometry(i, j) += dx;
        }
    }
    // m_CalculateGradient = g;
    return gradient;
}

bool ForceField::saveParametersToFile(const std::string& filename) const
{
    try {
        json output = exportCurrentParameters();

        std::ofstream file(filename);
        if (!file.is_open()) {
            CurcumaLogger::error(fmt::format("Cannot open file {} for writing", filename));
            return false;
        }

        file << output.dump(4); // Pretty print with 4 spaces
        file.close();

        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::success(fmt::format("Force field parameters saved to: {}", filename));
        }
        return true;

    } catch (const std::exception& e) {
        CurcumaLogger::error(fmt::format("Error saving parameters: {}", e.what()));
        return false;
    }
}

bool ForceField::loadParametersFromFile(const std::string& filename)
{
    try {
        std::ifstream file(filename);
        if (!file.is_open()) {
            CurcumaLogger::error(fmt::format("Cannot open file {} for reading", filename));
            return false;
        }

        json loaded_params;
        file >> loaded_params;
        file.close();

        // Validate that this is a force field parameter file
        if (!loaded_params.contains("method") || !loaded_params.contains("bonds")) {
            CurcumaLogger::error("Invalid force field parameter file format");
            return false;
        }

        // Apply loaded parameters
        setParameter(loaded_params);

        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::success(fmt::format("Force field parameters loaded from: {}", filename));
            CurcumaLogger::param("method", loaded_params["method"].get<std::string>());
            CurcumaLogger::param("bonds", static_cast<int>(loaded_params["bonds"].size()));
            CurcumaLogger::param("angles", static_cast<int>(loaded_params["angles"].size()));
        }

        return true;

    } catch (const std::exception& e) {
        CurcumaLogger::error(fmt::format("Error loading parameters: {}", e.what()));
        return false;
    }
}

json ForceField::exportCurrentParameters() const
{
    json output;

    // Method identification
    output["method"] = m_method;
    output["natoms"] = m_natoms;
    output["e0"] = m_e0;

    // Export bonds
    json bonds = json::array();
    for (const auto& bond : m_bonds) {
        json b;
        b["type"] = bond.type;
        b["i"] = bond.i;
        b["j"] = bond.j;
        b["k"] = bond.k;
        b["distance"] = bond.distance;
        b["fc"] = bond.fc;
        b["exponent"] = bond.exponent;
        b["r0_ij"] = bond.r0_ij;
        b["r0_ik"] = bond.r0_ik;
        bonds.push_back(b);
    }
    output["bonds"] = bonds;

    // Export angles
    json angles = json::array();
    for (const auto& angle : m_angles) {
        json a;
        a["type"] = angle.type;
        a["i"] = angle.i;
        a["j"] = angle.j;
        a["k"] = angle.k;
        a["fc"] = angle.fc;
        a["r0_ij"] = angle.r0_ij;
        a["r0_ik"] = angle.r0_ik;
        a["theta0_ijk"] = angle.theta0_ijk;
        a["C0"] = angle.C0;
        a["C1"] = angle.C1;
        a["C2"] = angle.C2;
        angles.push_back(a);
    }
    output["angles"] = angles;

    // Export dihedrals
    json dihedrals = json::array();
    for (const auto& dihedral : m_dihedrals) {
        json d;
        d["type"] = dihedral.type;
        d["i"] = dihedral.i;
        d["j"] = dihedral.j;
        d["k"] = dihedral.k;
        d["l"] = dihedral.l;
        d["V"] = dihedral.V;
        d["n"] = dihedral.n;
        d["phi0"] = dihedral.phi0;
        dihedrals.push_back(d);
    }
    output["dihedrals"] = dihedrals;

    // Export inversions
    json inversions = json::array();
    for (const auto& inversion : m_inversions) {
        json inv;
        inv["type"] = inversion.type;
        inv["i"] = inversion.i;
        inv["j"] = inversion.j;
        inv["k"] = inversion.k;
        inv["l"] = inversion.l;
        inv["fc"] = inversion.fc;
        inv["C0"] = inversion.C0;
        inv["C1"] = inversion.C1;
        inv["C2"] = inversion.C2;
        inversions.push_back(inv);
    }
    output["inversions"] = inversions;

    // Export vdW terms
    json vdws = json::array();
    for (const auto& vdw : m_vdWs) {
        json v;
        v["type"] = vdw.type;
        v["i"] = vdw.i;
        v["j"] = vdw.j;
        v["C_ij"] = vdw.C_ij;
        v["r0_ij"] = vdw.r0_ij;
        vdws.push_back(v);
    }
    output["vdws"] = vdws;

    // Export electrostatic terms
    json eqs = json::array();
    for (const auto& eq : m_EQs) {
        json e;
        e["type"] = eq.type;
        e["i"] = eq.i;
        e["j"] = eq.j;
        e["q_i"] = eq.q_i;
        e["q_j"] = eq.q_j;
        e["epsilon"] = eq.epsilon;
        eqs.push_back(e);
    }
    output["electrostatics"] = eqs;

    // Add metadata
    output["generated_by"] = "curcuma_forcefield";
    output["timestamp"] = std::chrono::system_clock::now().time_since_epoch().count();

    return output;
}

std::string ForceField::generateParameterFileName(const std::string& geometry_file)
{
    // input.xyz -> input.param.json
    // path/to/molecule.xyz -> path/to/molecule.param.json

    size_t last_dot = geometry_file.find_last_of('.');
    if (last_dot == std::string::npos) {
        return geometry_file + ".param.json";
    }

    std::string base = geometry_file.substr(0, last_dot);
    return base + ".param.json";
}

bool ForceField::tryLoadAutoParameters(const std::string& method)
{
    if (m_auto_param_file.empty()) {
        return false; // No auto-file detected
    }

    // Check if parameter file exists
    std::ifstream test_file(m_auto_param_file);
    if (!test_file.good()) {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("cache_status", fmt::format("No cached parameters found at: {}", m_auto_param_file));
        }
        return false;
    }
    test_file.close();

    // Try to load parameters
    if (!loadParametersFromFile(m_auto_param_file)) {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::warn(fmt::format("Failed to load parameters from: {}", m_auto_param_file));
        }
        return false;
    }

    // Check if method matches
    if (m_parameters.contains("method") && !m_parameters["method"].is_null()) {
        std::string cached_method = m_parameters["method"].get<std::string>();
        if (cached_method == method) {
            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::success(fmt::format("Loaded cached {} parameters from: {}", method, m_auto_param_file));
            }
            return true;
        } else {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::warn(fmt::format("Method mismatch in cached parameters (found: {}, expected: {})",
                    cached_method, method));
            }
            m_parameters.clear();
            return false;
        }
    } else {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::warn("No method field found in cached parameters or method is null");
        }
        m_parameters.clear();
        return false;
    }
}

bool ForceField::autoSaveParameters() const
{
    if (m_auto_param_file.empty()) {
        return false;
    }

    bool success = saveParametersToFile(m_auto_param_file);
    if (success) {
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::success(fmt::format("Auto-saved parameters to: {}", m_auto_param_file));
        }
    }
    return success;
}

double ForceField::Calculate(bool gradient)
{
    // Level 3+: Calculation debug info
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Starting force field calculation");
        CurcumaLogger::param("gradient", gradient ? "analytical" : "none");
        CurcumaLogger::param("stored_threads", static_cast<int>(m_stored_threads.size()));
    }

    m_gradient = Eigen::MatrixXd::Zero(m_geometry.rows(), 3);
    double energy = 0.0;
    double d4_energy = 0;
    double d3_energy = 0;
    double bond_energy = 0.0;
    double angle_energy = 0.0;
    double dihedral_energy = 0.0;
    double inversion_energy = 0.0;
    double vdw_energy = 0.0;
    double rep_energy = 0.0;
    double eq_energy = 0.0;
    double h4_energy = 0.0;
    double hh_energy = 0.0;

    for (int i = 0; i < m_stored_threads.size(); ++i) {
        m_stored_threads[i]->UpdateGeometry(m_geometry, gradient);
    }

    m_threadpool->Reset();
    m_threadpool->setActiveThreadCount(m_threads);

    m_threadpool->StartAndWait();
    // m_threadpool->setWakeUp(m_threadpool->WakeUp() / 2);

    for (int i = 0; i < m_stored_threads.size(); ++i) {
        bond_energy += m_stored_threads[i]->BondEnergy();
        angle_energy += m_stored_threads[i]->AngleEnergy();
        dihedral_energy += m_stored_threads[i]->DihedralEnergy();
        inversion_energy += m_stored_threads[i]->InversionEnergy();
        if (m_stored_threads[i]->Type() != 3)
            vdw_energy += m_stored_threads[i]->VdWEnergy();
        else
            h4_energy += m_stored_threads[i]->VdWEnergy();
        if (m_stored_threads[i]->Type() != 3)
            rep_energy += m_stored_threads[i]->RepEnergy();
        else
            hh_energy += m_stored_threads[i]->RepEnergy();
        // eq_energy += m_stored_threads[i]->RepEnergy();

        m_gradient += m_stored_threads[i]->Gradient();
    }

    energy = m_e0 + bond_energy + angle_energy + dihedral_energy + inversion_energy + vdw_energy + rep_energy + eq_energy + h4_energy + hh_energy;
    // Level 1+: Final energy result
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::energy_abs(energy, "Force Field Energy");
    }

    // Level 2+: Energy decomposition
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Force field energy decomposition:");
        if (m_e0 != 0.0) {
            CurcumaLogger::param("E0_baseline", fmt::format("{:.6f} Eh", m_e0));
        }
        CurcumaLogger::param("bond_energy", fmt::format("{:.6f} Eh", bond_energy));
        CurcumaLogger::param("angle_energy", fmt::format("{:.6f} Eh", angle_energy));
        CurcumaLogger::param("dihedral_energy", fmt::format("{:.6f} Eh", dihedral_energy));
        CurcumaLogger::param("inversion_energy", fmt::format("{:.6f} Eh", inversion_energy));
        CurcumaLogger::param("nonbonded_energy", fmt::format("{:.6f} Eh", vdw_energy + rep_energy));
        if (d3_energy != 0.0) {
            CurcumaLogger::param("D3_energy", fmt::format("{:.6f} Eh", d3_energy));
        }
        if (d4_energy != 0.0) {
            CurcumaLogger::param("D4_energy", fmt::format("{:.6f} Eh", d4_energy));
        }
        if (h4_energy != 0.0) {
            CurcumaLogger::param("HBond_correction", fmt::format("{:.6f} Eh", h4_energy));
        }
        if (hh_energy != 0.0) {
            CurcumaLogger::param("HH_repulsion", fmt::format("{:.6f} Eh", hh_energy));
        }

        if (gradient) {
            double grad_norm = m_gradient.norm();
            CurcumaLogger::param("gradient_norm", fmt::format("{:.6f} Eh/Bohr", grad_norm));
        }
    }

    // Level 3+: Thread and calculation details
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("final_energy", fmt::format("{:.6f} Eh", energy));
        CurcumaLogger::info("Force field calculation completed");
    }
    return energy;
}

// Claude Generated: Print comprehensive parameter summary
void ForceField::printParameterSummary() const
{
    if (m_parameters.empty()) {
        return;
    }

    try {
        // Count different parameter types - parameters are stored as arrays
        int bonds = 0, angles = 0, dihedrals = 0, inversions = 0, vdws = 0, esps = 0;

        if (m_parameters.contains("bonds") && m_parameters["bonds"].is_array()) {
            bonds = m_parameters["bonds"].size();
        }
        if (m_parameters.contains("angles") && m_parameters["angles"].is_array()) {
            angles = m_parameters["angles"].size();
        }
        if (m_parameters.contains("dihedrals") && m_parameters["dihedrals"].is_array()) {
            dihedrals = m_parameters["dihedrals"].size();
        }
        if (m_parameters.contains("inversions") && m_parameters["inversions"].is_array()) {
            inversions = m_parameters["inversions"].size();
        }
        if (m_parameters.contains("vdws") && m_parameters["vdws"].is_array()) {
            vdws = m_parameters["vdws"].size();
        }
        if (m_parameters.contains("esps") && m_parameters["esps"].is_array()) {
            esps = m_parameters["esps"].size();
        }

        // Print summary in consistent format
        fmt::print("  • {} bond terms\n", bonds);
        fmt::print("  • {} angle terms\n", angles);
        fmt::print("  • {} dihedral terms\n", dihedrals);
        fmt::print("  • {} inversion terms\n", inversions);
        fmt::print("  • {} van der Waals pairs\n", vdws);
        fmt::print("  • {} electrostatic pairs\n", esps);

        // Print scaling factors
        fmt::print("  Scaling factors:\n");
        if (m_parameters.contains("vdw_scaling")) {
            fmt::print("    vdW: {:.3f}\n", m_parameters["vdw_scaling"].get<double>());
        }
        if (m_parameters.contains("bond_scaling")) {
            fmt::print("    bond: {:.3f}\n", m_parameters["bond_scaling"].get<double>());
        }
        if (m_parameters.contains("angle_scaling")) {
            fmt::print("    angle: {:.3f}\n", m_parameters["angle_scaling"].get<double>());
        }
        if (m_parameters.contains("dihedral_scaling")) {
            fmt::print("    dihedral: {:.3f}\n", m_parameters["dihedral_scaling"].get<double>());
        }
        if (m_parameters.contains("inversion_scaling")) {
            fmt::print("    inversion: {:.3f}\n", m_parameters["inversion_scaling"].get<double>());
        }
        if (m_parameters.contains("coulomb_scaling")) {
            fmt::print("    coulomb: {:.3f}\n", m_parameters["coulomb_scaling"].get<double>());
        }
        if (m_parameters.contains("rep_scaling")) {
            fmt::print("    repulsion: {:.3f}\n", m_parameters["rep_scaling"].get<double>());
        }

        // Print dispersion and hydrogen bonding flags
        fmt::print("  Force field flags:\n");
        if (m_parameters.contains("d3") && m_parameters["d3"].get<double>() != 0) {
            fmt::print("    D3 dispersion: enabled (s6={:.3f}, s8={:.3f})\n",
                m_parameters.value("d3_s6", 0.0), m_parameters.value("d3_s8", 0.0));
        }
        if (m_parameters.contains("d4") && m_parameters["d4"].get<double>() != 0) {
            fmt::print("    D4 dispersion: enabled\n");
        }
        if (m_parameters.contains("h4") && m_parameters["h4"].get<double>() != 0) {
            fmt::print("    H4 hydrogen bonding: enabled (scaling={:.3f})\n",
                m_parameters.value("h4_scaling", 1.0));
            if (m_parameters.contains("h4_nh_o")) {
                fmt::print("      NH···O: {:.3f}\n", m_parameters["h4_nh_o"].get<double>());
            }
            if (m_parameters.contains("h4_oh_n")) {
                fmt::print("      OH···N: {:.3f}\n", m_parameters["h4_oh_n"].get<double>());
            }
        }

    } catch (const std::exception& e) {
        fmt::print("Warning: Could not display parameter summary: {}\n", e.what());
    }
}
