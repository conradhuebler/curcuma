/*
 * < Generic force field class for curcuma . >
 * Copyright (C) 2024 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "forcefield.h"
#include "src/core/elements.h"

#include "ff_workspace.h"
#include "uff_par.h"  // UFFParameterJson: controller defaults

#include "src/core/curcuma_logger.h"
#include "src/core/global.h"

#include <fmt/core.h>
#include <fmt/format.h>

#include <chrono>  // exportCurrentParameters() timestamp
#include <cmath>
#include <fstream>
#include <memory>

ForceField::ForceField(const json& controller)
{
    json parameter = MergeJson(UFFParameterJson, controller);

    m_threadpool = new CxxThreadPool();
    m_threadpool->setProgressBar(CxxThreadPool::ProgressBarType::None);
    m_threads = parameter["threads"];

    // Auto-detect parameter file based on geometry file
    if (parameter.contains("geometry_file")) {
        std::string geom_file = parameter["geometry_file"];
        m_auto_param_file = generateParameterFileName(geom_file);
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::param("auto_param_file", m_auto_param_file);
        }
    } else {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::warn("No geometry_file in ForceField parameter - automatic caching disabled");
        }
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
    std::string method_name = "unknown";
    if (parameters.contains("method") && !parameters["method"].is_null()) {
        method_name = parameters["method"].get<std::string>();
    }

    bool loaded_from_cache = false;

    // Claude Generated (December 2025): Prevent infinite recursion during cache loading
    // Skip cache loading if we're already being called from within loadParametersFromFile
    static thread_local bool loading_from_cache = false;
    // Coarse-grained parameters come from the user's JSON (cg_default, ...), not from the
    // UFF/QMDFF parameter cache: never load or save a .param.json for them (Sep 2026).
    const bool is_cg_method = (method_name == "cg" || method_name == "cg-lj");

    if (!is_cg_method && !loading_from_cache && !m_in_setParameter && m_enable_caching && parameters.contains("method")) {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Attempting to load cached parameters");
        }
        std::string method = parameters["method"];
        loading_from_cache = true;
        loaded_from_cache = tryLoadAutoParameters(method);
        loading_from_cache = false;

        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::param("loaded_from_cache", loaded_from_cache ? "true" : "false");
        }

        // If we loaded from cache successfully, we're done (tryLoadAutoParameters already called setParameter recursively)
        if (loaded_from_cache) {
            return;
        }
    } else {
        // Claude Generated (Dec 2025): Diagnostic warnings for caching issues
        if (!m_enable_caching && CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::warn("Parameter caching is disabled - parameters will be regenerated");
        }
        if (m_auto_param_file.empty() && m_enable_caching && CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::warn("No geometry_file provided - automatic caching disabled");
        }
    }

    // Claude Generated: Recursion guard - prevents infinite loops from nested setParameter calls
    if (m_in_setParameter) {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::warn("Recursive setParameter call detected - preventing infinite loop");
        }
        return;
    }
    m_in_setParameter = true;

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Initializing force field parameters");
        CurcumaLogger::param("method", method_name);
    }

    if (!loaded_from_cache) {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Cache miss - generating new force field parameters");
        }

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

        // uff-d3: pairwise D3 terms from D3ParameterGenerator (forwarded by ForceFieldGenerator)
        if (parameters.contains("d3_dispersion_pairs"))
            setD3DispersionPairs(parameters["d3_dispersion_pairs"]);

        m_parameters = parameters;
        m_method = m_parameters["method"];
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("method_selected", m_method);
        }
        if (m_method == "cg" || m_method == "cg-lj")
            generateCGParameters(parameters);
        if (m_parameters.contains("e0"))
            m_e0 = m_parameters["e0"];

        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Parameter generation complete");
        }

        // Auto-save new parameters (only if caching enabled)
        if (m_enable_caching && !is_cg_method) {
            autoSaveParameters();
        }
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        const bool is_qmdff = (m_method == "qmdff" || m_method == "quff");
        CurcumaLogger::param("method_type", fmt::format("{} ({})", is_qmdff ? 2 : 1, is_qmdff ? "QMDFF" : "UFF"));
    }

    // Claude Generated: Print parameter summary after setting parameters
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Generating parameter summary");
    }

    printParameterSummary();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Parameter summary completed");
    }

    // Claude Generated (March 2026): FFWorkspace evaluates the bonded + vdW (+ D3) terms.
    buildWorkspace();

    m_in_setParameter = false; // Reset the recursive guard - FIX: use member variable

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success("ForceField::setParameter() complete");
    }
}

void ForceField::setParameterFile(const std::string& file)
{
    if (!loadParametersFromFile(file)) {
        CurcumaLogger::warn(fmt::format("Failed to load parameter file: {}", file));
    }
}

// Claude Generated (Sep 2026): FFWorkspace construction, shared by setParameter() and the
// Calculate() safety fallback. qmdff/quff select the QMDFF term functions, everything else
// (uff, uff-d3) the UFF ones. uff-d3 additionally hands the D3 pair list to the workspace.
// ---------------------------------------------------------------------------
// Coarse-grained parameters (Sep 2026: restored from the removed thread engine)
// One vdW entry of type 3 per CG-CG pair; evaluated by FFWorkspace::calcCGPairs.
// ---------------------------------------------------------------------------
void ForceField::generateCGParameters(const json& cg_config)
{
    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::info("Generating CG parameters for coarse-grained simulation");
    if (!cg_config.contains("cg_default")) {
        CurcumaLogger::error("CG config missing required 'cg_default' section (use -load_ff_json FILE)");
        throw std::invalid_argument("CG configuration missing 'cg_default' section");
    }
    const auto& cg_default = cg_config["cg_default"];
    if (cg_default.contains("shape_vector")) {
        auto shape = cg_default["shape_vector"];
        if (shape.size() != 3) {
            CurcumaLogger::error("CG shape_vector must have exactly 3 elements");
            throw std::invalid_argument("Invalid shape_vector: expected 3 elements");
        }
        if (shape[0] <= 0 || shape[1] <= 0 || shape[2] <= 0)
            CurcumaLogger::warn("CG shape_vector contains non-positive values");
    }
    if (cg_default.contains("epsilon") && cg_default["epsilon"].get<double>() < 0)
        CurcumaLogger::warn("CG epsilon parameter is negative (unusual for attractive interaction)");

    m_vdWs.clear();
    if (cg_config.contains("bonds")) {
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info("Loading CG bond parameters");
        setBonds(cg_config["bonds"]);
    }
    json pair_overrides;
    if (cg_config.contains("pair_interactions"))
        pair_overrides = cg_config["pair_interactions"];

    int pair_count = 0;
    for (int i = 0; i < m_natoms; ++i) {
        for (int j = i + 1; j < m_natoms; ++j) {
            if (m_atom_types[i] != CG_ELEMENT || m_atom_types[j] != CG_ELEMENT) continue;
            vdW cg_pair;
            cg_pair.type = 3;
            cg_pair.i = i;
            cg_pair.j = j;
            cg_pair.shape_i = getCGShapeForAtom(i, cg_config);
            cg_pair.shape_j = getCGShapeForAtom(j, cg_config);
            cg_pair.orient_i = getCGOrientationForAtom(i, cg_config);
            cg_pair.orient_j = getCGOrientationForAtom(j, cg_config);
            const std::string pair_key = fmt::format("{}-{}", i, j);
            if (!pair_overrides.empty() && pair_overrides.contains(pair_key)) {
                const auto& pp = pair_overrides[pair_key];
                cg_pair.sigma = pp.value("sigma", cg_default.value("sigma", 4.0));
                cg_pair.epsilon = pp.value("epsilon", cg_default.value("epsilon", 0.0));
                cg_pair.cg_potential_type = pp.value("potential_type", cg_default.value("potential_type", 1));
            } else {
                cg_pair.sigma = cg_default.value("sigma", 4.0);
                cg_pair.epsilon = cg_default.value("epsilon", 0.0);
                cg_pair.cg_potential_type = cg_default.value("potential_type", 1);
            }
            m_vdWs.push_back(cg_pair);
            ++pair_count;
        }
    }
    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::success(fmt::format("Generated {} CG pair interactions", pair_count));
}

Eigen::Vector3d ForceField::getCGShapeForAtom(int atom_index, const json& config) const
{
    if (config.contains("cg_per_atom")) {
        const std::string key = std::to_string(atom_index);
        if (config["cg_per_atom"].contains(key) && config["cg_per_atom"][key].contains("shape_vector")) {
            const auto& s = config["cg_per_atom"][key]["shape_vector"];
            return Eigen::Vector3d(s[0], s[1], s[2]);
        }
    }
    if (config.contains("cg_default") && config["cg_default"].contains("shape_vector")) {
        const auto& s = config["cg_default"]["shape_vector"];
        return Eigen::Vector3d(s[0], s[1], s[2]);
    }
    return Eigen::Vector3d(2.0, 2.0, 2.0);
}

Eigen::Vector3d ForceField::getCGOrientationForAtom(int atom_index, const json& config) const
{
    if (config.contains("cg_per_atom")) {
        const std::string key = std::to_string(atom_index);
        if (config["cg_per_atom"].contains(key) && config["cg_per_atom"][key].contains("orientation")) {
            const auto& o = config["cg_per_atom"][key]["orientation"];
            return Eigen::Vector3d(o[0], o[1], o[2]);
        }
    }
    if (config.contains("cg_default") && config["cg_default"].contains("orientation")) {
        const auto& o = config["cg_default"]["orientation"];
        return Eigen::Vector3d(o[0], o[1], o[2]);
    }
    return Eigen::Vector3d(0.0, 0.0, 0.0);
}

void ForceField::buildWorkspace()
{
    const bool is_qmdff = (m_method == "qmdff" || m_method == "quff");
    const bool is_cg = (m_method == "cg" || m_method == "cg-lj");
    if (!is_qmdff && !is_cg && m_method != "uff" && m_method != "uff-d3") {
        CurcumaLogger::warn(fmt::format(
            "ForceField: method '{}' is not uff/uff-d3/qmdff/cg - evaluating with the UFF term functions", m_method));
    }

    ForceFieldParameterSet ws_params;
    ws_params.method_type = is_qmdff ? FFMethodType::QMDFF : (is_cg ? FFMethodType::CG : FFMethodType::UFF);
    ws_params.bonds      = m_bonds;
    ws_params.angles     = m_angles;
    ws_params.dihedrals  = m_dihedrals;
    ws_params.inversions = m_inversions;
    ws_params.vdws       = m_vdWs;
    if (m_method == "uff-d3" && !m_d3_dispersions.empty()) {
        ws_params.dispersions = m_d3_dispersions;
        ws_params.dispersion_method = "d3";
        ws_params.dispersion_enabled = true;
    } else {
        ws_params.dispersion_enabled = false;
    }
    ws_params.hbond_enabled = false;
    ws_params.repulsion_enabled = false;
    ws_params.coulomb_enabled = false;

    m_workspace = std::make_unique<FFWorkspace>(m_threads);
    m_workspace->setPool(m_threadpool);  // Claude Generated (March 2026): Required for multi-thread execution
    m_workspace->setAtomTypes(m_atom_types);
    m_workspace->setInteractionLists(std::move(ws_params));
    m_workspace->partition();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success(fmt::format("FFWorkspace (UFF/QMDFF): {} bonds, {} angles, {} vdWs, T={}",
            m_bonds.size(), m_angles.size(), m_vdWs.size(), m_threads));
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
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("angles_processing", fmt::format("Processing {} angle parameters", angles.size()));
    }
    m_angles.clear();
    for (int i = 0; i < angles.size(); ++i) {
        json angle = angles[i].get<json>();
        Angle a;

        a.type = angle["type"];

        a.i = angle["i"];
        a.j = angle["j"];
        a.k = angle["k"];

        // Claude Generated: Optional Fourier coefficients (UFF/QMDFF may use Fourier expansion)
        a.C0 = angle.value("C0", 0.0);
        a.C1 = angle.value("C1", 0.0);
        a.C2 = angle.value("C2", 0.0);

        a.fc = angle["fc"];
        a.r0_ij = angle["r0_ij"];
        a.r0_ik = angle["r0_ik"];
        a.theta0_ijk = angle["theta0_ijk"];
        m_angles.push_back(a);
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("angles_processed", fmt::format("{}", m_angles.size()));
    }
}

void ForceField::setDihedrals(const json& dihedrals)
{
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("dihedrals_processing", fmt::format("Processing {} dihedral parameters", dihedrals.size()));
    }
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

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("dihedrals_processed", fmt::format("{}", m_dihedrals.size()));
    }
}

void ForceField::setInversions(const json& inversions)
{
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("inversions_processing", fmt::format("Processing {} inversion parameters", inversions.size()));
    }
    m_inversions.clear();
    for (int i = 0; i < inversions.size(); ++i) {
        json inversion = inversions[i].get<json>();
        Inversion inv;
        inv.type = inversion["type"];

        inv.i = inversion["i"];
        inv.j = inversion["j"];
        inv.k = inversion["k"];
        inv.l = inversion["l"];

        // UFF/QMDFF style: Fourier coefficients
        inv.fc = inversion["fc"];
        inv.C0 = inversion["C0"];
        inv.C1 = inversion["C1"];
        inv.C2 = inversion["C2"];

        m_inversions.push_back(inv);
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("inversions_processed", fmt::format("{}", m_inversions.size()));
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

// Claude Generated (Sep 2026): was setGFNFFDispersions(); only the uff-d3 D3 pair list
// still arrives here. The BJ radius r0^2 = (a1*sqrt(r4r2ij) + a2)^2 with r4r2ij = C8/C6 is
// precomputed here exactly as before so FFWorkspace evaluates the same pair energy.
void ForceField::setD3DispersionPairs(const json& pairs)
{
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("setD3DispersionPairs: Loading {} dispersion pairs", pairs.size()));
    }

    m_d3_dispersions.clear();

    for (int i = 0; i < pairs.size(); ++i) {
        json disp_json = pairs[i].get<json>();
        GFNFFDispersion disp;

        disp.i = disp_json["i"];
        disp.j = disp_json["j"];
        disp.C6 = disp_json["C6"];
        disp.r_cut = disp_json["r_cut"];

        // r4r2ij = 3 * sqrtZr4r2_i * sqrtZr4r2_j (implicit C8/C6 factor)
        // r0_squared = (a1*sqrt(r4r2ij) + a2)^2
        if (disp_json.contains("r4r2ij") && disp_json.contains("r0_squared")) {
            disp.r4r2ij = disp_json["r4r2ij"];
            disp.r0_squared = disp_json["r0_squared"];
        } else {
            // D3ParameterGenerator output: derive from C8/a1/a2
            double a1 = disp_json.contains("a1") ? disp_json["a1"].get<double>() : 0.58;
            double a2 = disp_json.contains("a2") ? disp_json["a2"].get<double>() : 4.80;

            if (disp_json.contains("C8") && disp.C6 > 1e-10) {
                double c8 = disp_json["C8"].get<double>();
                disp.r4r2ij = c8 / disp.C6;
            } else {
                disp.r4r2ij = 1.0;
            }
            disp.r0_squared = std::pow(a1 * std::sqrt(disp.r4r2ij) + a2, 2);
        }

        // Zeta charge scaling (GFN-FF only; 1.0 for D3 pairs)
        disp.zetac6 = disp_json.contains("zetac6") ? disp_json["zetac6"].get<double>() : 1.0;

        m_d3_dispersions.push_back(disp);
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success(fmt::format("Loaded {} D3 dispersion pairs", m_d3_dispersions.size()));
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

Eigen::MatrixXd ForceField::NumGrad()
{
    Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(m_natoms, 3);

    double dx = 1e-6; // m_d;
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
            if (loaded_params.contains("d3_dispersion_pairs")) {
                CurcumaLogger::param("d3_dispersion_pairs", static_cast<int>(loaded_params["d3_dispersion_pairs"].size()));
            }
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

    // Claude Generated: Export D3 dispersion parameters if they exist (uff-d3).
    // The generator output is stored as-is so a cache reload re-parses exactly what a
    // fresh generation would have parsed.
    if (m_parameters.contains("d3_dispersion_pairs")) {
        output["d3_dispersion_pairs"] = m_parameters["d3_dispersion_pairs"];
    }
    if (m_parameters.contains("d3_damping")) {
        output["d3_damping"] = m_parameters["d3_damping"];
    }
    if (m_parameters.contains("d3_enabled")) {
        output["d3_enabled"] = m_parameters["d3_enabled"];
    }

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
            // Claude Generated (Dec 2025): Show cache success at verbosity ≥1 (important user info)
            if (CurcumaLogger::get_verbosity() >= 1) {
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
    // Claude Generated (March 2026): Safety fallback — the workspace is built in setParameter();
    // rebuild it here if interaction lists exist without one.
    if (!m_workspace && !m_bonds.empty()) {
        CurcumaLogger::warn(fmt::format(
            "UFF/QMDFF workspace not initialized for method '{}' — rebuilding", m_method));
        buildWorkspace();
    }

    if (!m_workspace) {
        CurcumaLogger::error("ForceField::Calculate() - ForceField not initialized! No parameters available.");
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("This usually means setParameter() was not called or parameter generation failed.");
        }
        return 0.0;
    }

    m_workspace->setGeometry(m_geometry);
    double energy = m_workspace->calculate(gradient);
    if (gradient)
        m_gradient = m_workspace->gradient();
    const auto& e = m_workspace->energyComponents();
    m_bond_energy        = e.bond;
    m_angle_energy       = e.angle;
    m_dihedral_energy    = e.dihedral;
    m_inversion_energy   = e.inversion;
    m_vdw_energy         = e.vdw;
    m_rep_energy         = e.rep;
    m_dispersion_energy  = e.dispersion;
    m_d3_energy          = e.dispersion;  // D3 dispersion stored here for uff-d3
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
        int bonds = 0, angles = 0, dihedrals = 0, inversions = 0, vdws = 0, esps = 0, d3_pairs = 0;

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
        if (m_parameters.contains("d3_dispersion_pairs") && m_parameters["d3_dispersion_pairs"].is_array()) {
            d3_pairs = m_parameters["d3_dispersion_pairs"].size();
        }

        if (CurcumaLogger::get_verbosity() >= 2) {
            // Basic force field topology
            CurcumaLogger::info("Force field topology summary:");
            if (bonds > 0) CurcumaLogger::param(fmt::format("bonds_count", bonds), fmt::format("{} bond terms", bonds));
            if (angles > 0) CurcumaLogger::param(fmt::format("angles_count", angles), fmt::format("{} angle terms", angles));
            if (dihedrals > 0) CurcumaLogger::param(fmt::format("dihedrals_count", dihedrals), fmt::format("{} dihedral terms", dihedrals));
            if (inversions > 0) CurcumaLogger::param(fmt::format("inversions_count", inversions), fmt::format("{} inversion terms", inversions));
            if (vdws > 0) CurcumaLogger::param(fmt::format("vdw_count", vdws), fmt::format("{} van der Waals pairs", vdws));
            if (esps > 0) CurcumaLogger::param(fmt::format("electrostatic_count", esps), fmt::format("{} electrostatic pairs", esps));
            if (d3_pairs > 0) CurcumaLogger::param("d3_dispersion_pairs", fmt::format("{}", d3_pairs));

            // Print scaling factors
            CurcumaLogger::info("Scaling factors:");
            if (m_parameters.contains("vdw_scaling")) {
                CurcumaLogger::param("vdw_scaling", fmt::format("{:.3f}", m_parameters["vdw_scaling"].get<double>()));
            }
            if (m_parameters.contains("bond_scaling")) {
                CurcumaLogger::param("bond_scaling", fmt::format("{:.3f}", m_parameters["bond_scaling"].get<double>()));
            }
            if (m_parameters.contains("angle_scaling")) {
                CurcumaLogger::param("angle_scaling", fmt::format("{:.3f}", m_parameters["angle_scaling"].get<double>()));
            }
            if (m_parameters.contains("dihedral_scaling")) {
                CurcumaLogger::param("dihedral_scaling", fmt::format("{:.3f}", m_parameters["dihedral_scaling"].get<double>()));
            }
            if (m_parameters.contains("inversion_scaling")) {
                CurcumaLogger::param("inversion_scaling", fmt::format("{:.3f}", m_parameters["inversion_scaling"].get<double>()));
            }
            if (m_parameters.contains("coulomb_scaling")) {
                CurcumaLogger::param("coulomb_scaling", fmt::format("{:.3f}", m_parameters["coulomb_scaling"].get<double>()));
            }
            if (m_parameters.contains("rep_scaling")) {
                CurcumaLogger::param("repulsion_scaling", fmt::format("{:.3f}", m_parameters["rep_scaling"].get<double>()));
            }
        }

        if (CurcumaLogger::get_verbosity() >= 2) {
            // Print dispersion and hydrogen bonding flags
            CurcumaLogger::info("Force field flags:");
            if (m_parameters.contains("d3") && m_parameters["d3"].get<double>() != 0) {
                CurcumaLogger::info("D3 dispersion enabled");
                CurcumaLogger::param("d3_s6", fmt::format("{:.3f}", m_parameters.value("d3_s6", 0.0)));
                CurcumaLogger::param("d3_s8", fmt::format("{:.3f}", m_parameters.value("d3_s8", 0.0)));
            }
            if (m_parameters.contains("d4") && m_parameters["d4"].get<double>() != 0) {
                CurcumaLogger::info("D4 dispersion enabled");
            }
            if (m_parameters.contains("h4") && m_parameters["h4"].get<double>() != 0) {
                CurcumaLogger::info("H4 hydrogen bonding enabled");
                CurcumaLogger::param("h4_scaling", fmt::format("{:.3f}", m_parameters.value("h4_scaling", 1.0)));
                if (m_parameters.contains("h4_nh_o")) {
                    CurcumaLogger::param("h4_nh_o_scaling", fmt::format("{:.3f}", m_parameters["h4_nh_o"].get<double>()));
                }
                if (m_parameters.contains("h4_oh_n")) {
                    CurcumaLogger::param("h4_oh_n_scaling", fmt::format("{:.3f}", m_parameters["h4_oh_n"].get<double>()));
                }
            }
        }

    } catch (const std::exception& e) {
        CurcumaLogger::warn(fmt::format("Could not display parameter summary: {}", e.what()));
    }
}
