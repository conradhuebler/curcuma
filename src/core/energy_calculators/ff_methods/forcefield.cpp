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
#include "src/core/global.h"
#include "uff_par.h"

#include "cg_potentials.h"
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

void ForceField::distributeEEQCharges(const Vector& charges)
{
    // Phase 5A: Distribute EEQ charges to all threads for fqq calculation
    // Claude Generated (Nov 2025)
    for (int i = 0; i < m_stored_threads.size(); ++i) {
        m_stored_threads[i]->setEEQCharges(charges);
    }
}

void ForceField::setParameter(const json& parameters)
{
    // Claude Generated: Fix recursive guard - use member variable, not static
    if (m_in_setParameter) {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::warn("Recursive setParameter call detected - preventing infinite loop");
        }
        std::cout << "DEBUG: setParameter() RECURSION - returning early!" << std::endl;
        return;
    }
    m_in_setParameter = true;

    // Claude Generated: Debug - always log setParameter entry
    // std::cout << "DEBUG: ForceField::setParameter() ENTRY (set in_setParameter=true)" << std::endl;

    std::string method_name = "unknown";
    if (parameters.contains("method") && !parameters["method"].is_null()) {
        method_name = parameters["method"].get<std::string>();
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Initializing force field parameters");
        CurcumaLogger::param("method", method_name);
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Force field setup: parameters validated");
    }

    bool loaded_from_cache = false;

    // Try loading cached parameters if caching enabled and same method
    if (m_enable_caching && parameters.contains("method")) {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Attempting to load cached parameters");
        }
        std::string method = parameters["method"];
        loaded_from_cache = tryLoadAutoParameters(method);
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::param("loaded_from_cache", loaded_from_cache ? "true" : "false");
        }
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

        // Phase 4.2: GFN-FF pairwise non-bonded parameters (Claude Generated 2025)
        // Support both "gfnff_dispersions" (from GFNFF) and "d3_dispersion_pairs" (from D3ParameterGenerator)
        if (parameters.contains("gfnff_dispersions"))
            setGFNFFDispersions(parameters["gfnff_dispersions"]);
        else if (parameters.contains("d3_dispersion_pairs"))
            setGFNFFDispersions(parameters["d3_dispersion_pairs"]);

        if (parameters.contains("gfnff_repulsions"))
            setGFNFFRepulsions(parameters["gfnff_repulsions"]);
        if (parameters.contains("gfnff_coulombs"))
            setGFNFFCoulombs(parameters["gfnff_coulombs"]);

        // Phase 3: GFN-FF hydrogen bond and halogen bond setters (Claude Generated 2025)
        if (parameters.contains("gfnff_hbonds"))
            setGFNFFHydrogenBonds(parameters["gfnff_hbonds"]);
        if (parameters.contains("gfnff_xbonds"))
            setGFNFFHalogenBonds(parameters["gfnff_xbonds"]);

        m_parameters = parameters;
        m_method = m_parameters["method"];
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("method_selected", m_method);
        }
        if (m_parameters.contains("e0"))
            m_e0 = m_parameters["e0"];

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Calculating parameter ranges");
        }
        std::cout << "DEBUG: About to call AutoRanges()..." << std::endl;
        AutoRanges();
        std::cout << "DEBUG: AutoRanges() completed" << std::endl;

        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Parameter generation complete");
        }

        // Auto-save new parameters (only if caching enabled)
        if (m_enable_caching) {
            autoSaveParameters();
        }
    }

    // Set method type for ForceFieldThread (regardless of cache)
    int method_type = 1; // default UFF
    if (m_method == "qmdff" || m_method == "quff") {
        method_type = 2;
    } else if (m_method == "gfnff" || m_method == "cgfnff") {  // Claude Generated (Dec 2025): Accept both gfnff and cgfnff
        method_type = 3; // GFN-FF
    } else if (m_method == "cg" || m_method == "cg-lj") {
        method_type = 4; // CG methods
        // Generate CG parameters if not loaded from cache
        if (m_enable_caching && !loaded_from_cache) {
            generateCGParameters(parameters);
        }
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        std::string method_name = "UFF"; // default
        if (method_type == 2) method_name = "QMDFF";
        else if (method_type == 3) method_name = "GFN-FF";
        else if (method_type == 4) method_name = "CG";
        CurcumaLogger::param("method_type", fmt::format("{} ({})", method_type, method_name));
    }

    // Claude Generated: Print parameter summary after setting parameters
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Generating parameter summary");
    }

    printParameterSummary();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Parameter summary completed");
    }

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

// Claude Generated: CG parameter generation for sphere-based coarse graining
// Prepared for future ellipsoid extension but currently implements spheres only
void ForceField::generateCGParameters(const json& cg_config)
{
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Generating CG parameters for coarse-grained simulation");
    }

    // Validate required CG config sections
    if (!cg_config.contains("cg_default")) {
        CurcumaLogger::error("CG config missing required 'cg_default' section");
        throw std::invalid_argument("CG configuration missing 'cg_default' section");
    }

    const auto& cg_default = cg_config["cg_default"];

    // Validate shape vector if present
    if (cg_default.contains("shape_vector")) {
        auto shape = cg_default["shape_vector"];
        if (shape.size() != 3) {
            CurcumaLogger::error("CG shape_vector must have exactly 3 elements");
            throw std::invalid_argument("Invalid shape_vector: expected 3 elements");
        }
        if (shape[0] <= 0 || shape[1] <= 0 || shape[2] <= 0) {
            CurcumaLogger::warn("CG shape_vector contains non-positive values");
        }
    }

    // Validate epsilon (warn if negative)
    if (cg_default.contains("epsilon") && cg_default["epsilon"].get<double>() < 0) {
        CurcumaLogger::warn("CG epsilon parameter is negative (unusual for attractive interaction)");
    }

    // Clear existing vdW interactions to prepare for CG pairs
    m_vdWs.clear();

    // Load bond parameters if present (for constraint bonds between CG particles)
    if (cg_config.contains("bonds")) {
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Loading CG bond parameters");
        }
        setBonds(cg_config["bonds"]);
    }

    // Get pair interaction overrides if present
    json pair_overrides;
    if (cg_config.contains("pair_interactions")) {
        pair_overrides = cg_config["pair_interactions"];
    }

    // Generate individual vdW structs for each CG atom pair
    int pair_count = 0;
    for (int i = 0; i < m_natoms; ++i) {
        for (int j = i + 1; j < m_natoms; ++j) {

            // Only CG-CG interactions (both atoms must be CG_ELEMENT)
            if (m_atom_types[i] == CG_ELEMENT && m_atom_types[j] == CG_ELEMENT) {

                vdW cg_pair;
                cg_pair.type = 3; // CG type
                cg_pair.i = i;
                cg_pair.j = j;

                // Load shape parameters - currently sphere setup (pragmatic approach)
                cg_pair.shape_i = getCGShapeForAtom(i, cg_config);
                cg_pair.shape_j = getCGShapeForAtom(j, cg_config);

                // Orientation parameters (prepared for ellipsoids, unused for spheres)
                cg_pair.orient_i = getCGOrientationForAtom(i, cg_config);
                cg_pair.orient_j = getCGOrientationForAtom(j, cg_config);

                // Check for pair-specific interaction parameter overrides
                std::string pair_key = fmt::format("{}-{}", i, j);
                if (!pair_overrides.empty() && pair_overrides.contains(pair_key)) {
                    // Load custom pair parameters
                    const auto& pair_params = pair_overrides[pair_key];
                    cg_pair.sigma = pair_params.value("sigma", cg_default.value("sigma", 4.0));
                    cg_pair.epsilon = pair_params.value("epsilon", cg_default.value("epsilon", 0.0));
                    cg_pair.cg_potential_type = pair_params.value("potential_type", cg_default.value("potential_type", 1));

                    if (CurcumaLogger::get_verbosity() >= 2) {
                        CurcumaLogger::info(fmt::format("Applied custom parameters to pair {}-{}", i, j));
                    }
                } else {
                    // Use default parameters
                    cg_pair.sigma = cg_default.value("sigma", 4.0);
                    cg_pair.epsilon = cg_default.value("epsilon", 0.0);
                    cg_pair.cg_potential_type = cg_default.value("potential_type", 1);
                }

                m_vdWs.push_back(cg_pair);
                pair_count++;

                if (CurcumaLogger::get_verbosity() >= 3) {
                    CurcumaLogger::param(fmt::format("CG pair {}-{}", i, j),
                        fmt::format("sigma={:.2f}, epsilon={:.4f}, type={}",
                            cg_pair.sigma, cg_pair.epsilon, cg_pair.cg_potential_type));
                }
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success(fmt::format("Generated {} CG pair interactions", pair_count));
    }
}

// Claude Generated: Get CG shape parameters for specific atom (sphere setup)
Eigen::Vector3d ForceField::getCGShapeForAtom(int atom_index, const json& config)
{
    // Check for per-atom parameters
    if (config.contains("cg_per_atom")) {
        std::string atom_key = std::to_string(atom_index);
        if (config["cg_per_atom"].contains(atom_key)) {
            auto& shape = config["cg_per_atom"][atom_key]["shape_vector"];
            return Eigen::Vector3d(shape[0], shape[1], shape[2]);
        }
    }

    // Fall back to default - pragmatic sphere setup
    if (config.contains("cg_default")) {
        auto& shape = config["cg_default"]["shape_vector"];
        return Eigen::Vector3d(shape[0], shape[1], shape[2]);
    }

    // Ultimate fallback: sphere with radius 2.0
    return Eigen::Vector3d(2.0, 2.0, 2.0);
}

// Claude Generated: Get CG orientation parameters for specific atom (prepared for ellipsoids)
Eigen::Vector3d ForceField::getCGOrientationForAtom(int atom_index, const json& config)
{
    // Check for per-atom orientation parameters
    if (config.contains("cg_per_atom")) {
        std::string atom_key = std::to_string(atom_index);
        if (config["cg_per_atom"].contains(atom_key) && config["cg_per_atom"][atom_key].contains("orientation")) {
            auto& orient = config["cg_per_atom"][atom_key]["orientation"];
            return Eigen::Vector3d(orient[0], orient[1], orient[2]);
        }
    }

    // Fall back to default orientation
    if (config.contains("cg_default") && config["cg_default"].contains("orientation")) {
        auto& orient = config["cg_default"]["orientation"];
        return Eigen::Vector3d(orient[0], orient[1], orient[2]);
    }

    // No rotation (spheres don't need orientation)
    return Eigen::Vector3d(0.0, 0.0, 0.0);
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

        // Claude Generated: Optional Fourier coefficients (not used by GFN-FF)
        // GFN-FF uses simple angle bending, UFF/QMDFF may use Fourier expansion
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

        // GFN-FF uses different inversion parameters (barrier, omega0) vs UFF (fc, C0, C1, C2)
        if (inversion.contains("barrier")) {
            // GFN-FF style: use barrier and omega0
            inv.fc = inversion["barrier"];
            inv.C0 = inversion.value("omega0", 0.0);
            inv.C1 = 0.0;
            inv.C2 = 0.0;
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::param("inversion_format", "GFN-FF (barrier/omega0)");
            }
        } else {
            // UFF/QMDFF style: use Fourier coefficients
            inv.fc = inversion["fc"];
            inv.C0 = inversion["C0"];
            inv.C1 = inversion["C1"];
            inv.C2 = inversion["C2"];
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::param("inversion_format", "UFF (Fourier)");
            }
        }

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

// Phase 4.2: GFN-FF pairwise non-bonded parameter setters (Claude Generated 2025)

void ForceField::setGFNFFDispersions(const json& dispersions)
{
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("setGFNFFDispersions: Loading {} dispersions", dispersions.size()));
    }

    m_gfnff_dispersions.clear();
    for (int i = 0; i < dispersions.size(); ++i) {
        json disp_json = dispersions[i].get<json>();
        GFNFFDispersion disp;

        disp.i = disp_json["i"];
        disp.j = disp_json["j"];
        disp.C6 = disp_json["C6"];
        disp.C8 = disp_json["C8"];
        disp.s6 = disp_json["s6"];
        disp.s8 = disp_json["s8"];
        disp.a1 = disp_json["a1"];
        disp.a2 = disp_json["a2"];
        disp.r_cut = disp_json["r_cut"];

        m_gfnff_dispersions.push_back(disp);
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success(fmt::format("Loaded {} GFN-FF dispersion pairs", m_gfnff_dispersions.size()));
    }
}

void ForceField::setGFNFFRepulsions(const json& repulsions)
{
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("setGFNFFRepulsions: Loading {} repulsions", repulsions.size()));
    }

    m_gfnff_repulsions.clear();
    for (int i = 0; i < repulsions.size(); ++i) {
        json rep_json = repulsions[i].get<json>();
        GFNFFRepulsion rep;

        rep.i = rep_json["i"];
        rep.j = rep_json["j"];
        rep.alpha = rep_json["alpha"];
        rep.repab = rep_json["repab"];
        rep.r_cut = rep_json["r_cut"];

        m_gfnff_repulsions.push_back(rep);
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success(fmt::format("Loaded {} GFN-FF repulsion pairs", m_gfnff_repulsions.size()));
    }
}

void ForceField::setGFNFFCoulombs(const json& coulombs)
{
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("setGFNFFCoulombs: Loading {} Coulomb pairs", coulombs.size()));
    }

    m_gfnff_coulombs.clear();
    for (int i = 0; i < coulombs.size(); ++i) {
        json coul_json = coulombs[i].get<json>();
        GFNFFCoulomb coul;

        coul.i = coul_json["i"];
        coul.j = coul_json["j"];
        coul.q_i = coul_json["q_i"];
        coul.q_j = coul_json["q_j"];
        coul.gamma_ij = coul_json["gamma_ij"];
        coul.chi_i = coul_json.value("chi_i", 0.0);      // Default to 0 if missing (backward compat)
        coul.chi_j = coul_json.value("chi_j", 0.0);
        coul.gam_i = coul_json.value("gam_i", 0.0);      // Chemical hardness
        coul.gam_j = coul_json.value("gam_j", 0.0);
        coul.alp_i = coul_json.value("alp_i", 0.0);
        coul.alp_j = coul_json.value("alp_j", 0.0);
        coul.r_cut = coul_json.value("r_cut", 50.0);

        m_gfnff_coulombs.push_back(coul);
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success(fmt::format("Loaded {} GFN-FF Coulomb pairs", m_gfnff_coulombs.size()));
    }
}

void ForceField::setGFNFFHydrogenBonds(const json& hbonds)
{
    // Claude Generated (2025): Phase 3 - HB Parameter Setter
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("setGFNFFHydrogenBonds: Loading {} hydrogen bonds", hbonds.size()));
    }

    m_gfnff_hbonds.clear();
    for (const auto& hb : hbonds) {
        GFNFFHydrogenBond bond;

        bond.i = hb["i"];  // Donor atom A
        bond.j = hb["j"];  // Hydrogen
        bond.k = hb["k"];  // Acceptor atom B

        bond.basicity_A = hb["basicity_A"];
        bond.basicity_B = hb["basicity_B"];
        bond.acidity_A = hb["acidity_A"];

        bond.q_H = hb["q_H"];
        bond.q_A = hb["q_A"];
        bond.q_B = hb["q_B"];

        bond.case_type = hb.value("case", 1);

        if (bond.case_type == 2 && hb.contains("neighbors_B")) {
            bond.neighbors_B = hb["neighbors_B"].get<std::vector<int>>();
        }

        m_gfnff_hbonds.push_back(bond);
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success(fmt::format("Loaded {} GFN-FF hydrogen bonds", hbonds.size()));
    }
}

void ForceField::setGFNFFHalogenBonds(const json& xbonds)
{
    // Claude Generated (2025): Phase 3 - XB Parameter Setter
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("setGFNFFHalogenBonds: Loading {} halogen bonds", xbonds.size()));
    }

    m_gfnff_xbonds.clear();
    for (const auto& xb : xbonds) {
        GFNFFHalogenBond bond;

        bond.i = xb["i"];  // Donor atom A
        bond.j = xb["j"];  // Halogen X
        bond.k = xb["k"];  // Acceptor atom B

        bond.basicity_B = xb["basicity_B"];
        bond.acidity_X = xb["acidity_X"];

        bond.q_X = xb["q_X"];
        bond.q_B = xb["q_B"];

        m_gfnff_xbonds.push_back(bond);
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success(fmt::format("Loaded {} GFN-FF halogen bonds", xbonds.size()));
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

    // Claude Generated: Clear stored threads before creating new ones
    // This prevents thread accumulation in AutoRanges() if called multiple times
    if (!m_stored_threads.empty()) {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::warn("Clearing " + std::to_string(m_stored_threads.size()) + " existing threads before creating new ones");
        }
        for (auto t : m_stored_threads) {
            delete t;
        }
        m_stored_threads.clear();
    }

    int free_threads = m_threads;
#pragma message("revert")
    int d3 = false; // m_parameters["d3"];
    #ifdef USE_D3
    if (d3) {
        if (free_threads > 1)
            free_threads--;
        D3Thread* thread = new D3Thread(m_threads - 1, free_threads);
        thread->setParamater(m_parameters);
        thread->Initialise(m_atom_types);
        m_threadpool->addThread(thread);
        m_stored_threads.push_back(thread);
    }
    #endif
    // TEMPORARILY DISABLED - H4Thread has build errors
    // Claude Generated Comment - 2025-11-30
    /*
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
    */

    // Claude Generated: Ensure at least 1 thread is created for energy calculations
    // (free_threads could be 0 if all threads are reserved for D3/H4 corrections)
    int thread_count = (free_threads > 0) ? free_threads : 1;

    for (int i = 0; i < thread_count; ++i) {
        ForceFieldThread* thread = new ForceFieldThread(i, thread_count);
        thread->setGeometry(m_geometry, false);
        thread->Initialise(m_atom_types);  // Phase 3: Initialize atom types for covalent radius calculations
        m_threadpool->addThread(thread);
        m_stored_threads.push_back(thread);

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("created_thread", i);
        }
        if (std::find(m_uff_methods.begin(), m_uff_methods.end(), m_method) != m_uff_methods.end()) {
            thread->setMethod(1);
        } else if (std::find(m_qmdff_methods.begin(), m_qmdff_methods.end(), m_method) != m_qmdff_methods.end()) {
            thread->setMethod(2);
        } else if (m_method == "gfnff" || m_method == "cgfnff") { // Claude Generated (2025-12-13): Support both method names
            thread->setMethod(3); // GFN-FF

            // Phase 2: Configure GFN-FF parameter flags (Claude Generated Dec 2025)
            // These control which energy terms are calculated to save CPU time
            bool dispersion = m_parameters.value("dispersion", true);
            bool hbond = m_parameters.value("hbond", true);
            bool repulsion = m_parameters.value("repulsion", true);
            bool coulomb = m_parameters.value("coulomb", true);

            thread->setDispersionEnabled(dispersion);
            thread->setHydrogenBondEnabled(hbond);
            thread->setRepulsionEnabled(repulsion);
            thread->setCoulombEnabled(coulomb);
        }
        for (int j = int(i * m_bonds.size() / double(free_threads)); j < int((i + 1) * m_bonds.size() / double(free_threads)); ++j) {
            if (m_method == "gfnff" || m_method == "cgfnff") { // Claude Generated (2025-12-13): Support both method names
                thread->addGFNFFBond(m_bonds[j]);
            } else {
                thread->addBond(m_bonds[j]);
            }
        }

        for (int j = int(i * m_angles.size() / double(free_threads)); j < int((i + 1) * m_angles.size() / double(free_threads)); ++j) {
            if (m_method == "gfnff" || m_method == "cgfnff") {
                thread->addGFNFFAngle(m_angles[j]);
            } else {
                thread->addAngle(m_angles[j]);
            }
        }

        for (int j = int(i * m_dihedrals.size() / double(free_threads)); j < int((i + 1) * m_dihedrals.size() / double(free_threads)); ++j) {
            if (m_method == "gfnff" || m_method == "cgfnff") {
                thread->addGFNFFDihedral(m_dihedrals[j]);
            } else {
                thread->addDihedral(m_dihedrals[j]);
            }
        }

        for (int j = int(i * m_inversions.size() / double(free_threads)); j < int((i + 1) * m_inversions.size() / double(free_threads)); ++j) {
            if (m_method == "gfnff" || m_method == "cgfnff") {
                thread->addGFNFFInversion(m_inversions[j]);
            } else {
                thread->addInversion(m_inversions[j]);
            }
        }

        for (int j = int(i * m_vdWs.size() / double(free_threads)); j < int((i + 1) * m_vdWs.size() / double(free_threads)); ++j) {
            if (m_method == "gfnff" || m_method == "cgfnff") {
                thread->addGFNFFvdW(m_vdWs[j]);
            } else {
                thread->addvdW(m_vdWs[j]);
            }
        }

        // Phase 4.2: Distribute GFN-FF pairwise non-bonded interactions (Claude Generated 2025)
        if (m_method == "gfnff" || m_method == "cgfnff") { // Claude Generated (2025-12-13): Support both method names
            for (int j = int(i * m_gfnff_dispersions.size() / double(free_threads)); j < int((i + 1) * m_gfnff_dispersions.size() / double(free_threads)); ++j) {
                thread->addGFNFFDispersion(m_gfnff_dispersions[j]);
            }

            for (int j = int(i * m_gfnff_repulsions.size() / double(free_threads)); j < int((i + 1) * m_gfnff_repulsions.size() / double(free_threads)); ++j) {
                thread->addGFNFFRepulsion(m_gfnff_repulsions[j]);
            }

            for (int j = int(i * m_gfnff_coulombs.size() / double(free_threads)); j < int((i + 1) * m_gfnff_coulombs.size() / double(free_threads)); ++j) {
                thread->addGFNFFCoulomb(m_gfnff_coulombs[j]);
            }

            // Phase 6: Distribute atoms for self-energy calculation (Claude Generated Dec 2025)
            // CRITICAL: Each atom assigned to EXACTLY ONE thread to avoid duplicate self-energy
            // Distribute all atoms (0 to natoms-1) across threads
            int start_atom = int(i * m_natoms / double(free_threads));
            int end_atom = int((i + 1) * m_natoms / double(free_threads));
            std::vector<int> assigned_atoms;
            for (int atom_id = start_atom; atom_id < end_atom; ++atom_id) {
                assigned_atoms.push_back(atom_id);
            }
            thread->assignAtomsForSelfEnergy(assigned_atoms);

            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::param(fmt::format("thread_{}_atoms_for_self_energy", i),
                    fmt::format("[{}, {}) = {} atoms", start_atom, end_atom, assigned_atoms.size()));
            }
        }

        for (int j = int(i * m_EQs.size() / double(free_threads)); j < int((i + 1) * m_EQs.size() / double(free_threads)); ++j)
            thread->addEQ(m_EQs[j]);
    }

    // Level 3+: AutoRanges completion info
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("created_threads", thread_count);
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
    // Claude Generated (2025-12-13): Check for uninitialized ForceField
    if (m_stored_threads.empty()) {
        CurcumaLogger::error("ForceField::Calculate() - ForceField not initialized! No threads available.");
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("This usually means setParameter() was not called or AutoRanges() failed.");
            CurcumaLogger::info("Check if parameter generation succeeded and produced valid JSON arrays.");
        }
        return 0.0;
    }

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
    // Claude Generated: Reset all energy components for regression testing (Nov 2025)
    m_bond_energy = 0.0;
    m_angle_energy = 0.0;
    m_dihedral_energy = 0.0;
    m_inversion_energy = 0.0;
    m_vdw_energy = 0.0;
    m_rep_energy = 0.0;
    m_eq_energy = 0.0;
    m_dispersion_energy = 0.0;
    m_coulomb_energy = 0.0;
    m_energy_hbond = 0.0;    // Claude Generated (2025): Phase 5 - Reset HB energy
    m_energy_xbond = 0.0;    // Claude Generated (2025): Phase 5 - Reset XB energy
    m_hh_energy = 0.0;       // Claude Generated (Dec 2025, Session 9): Reset GFN-FF HH repulsion

    double h4_energy = 0.0;
    double cg_energy = 0.0; // Claude Generated: CG pair interaction energy

    // Claude Generated: GFN-FF specific non-bonded energies (Phase 4.4 fix)
    // NOTE: These are now accumulated into member variables m_dispersion_energy and m_coulomb_energy
    // (local variables removed - using m_dispersion_energy and m_coulomb_energy member variables)

    for (int i = 0; i < m_stored_threads.size(); ++i) {
        m_stored_threads[i]->UpdateGeometry(m_geometry, gradient);
    }

    m_threadpool->Reset();
    m_threadpool->setActiveThreadCount(m_threads);

    m_threadpool->StartAndWait();
    // m_threadpool->setWakeUp(m_threadpool->WakeUp() / 2);

    for (int i = 0; i < m_stored_threads.size(); ++i) {
        m_bond_energy += m_stored_threads[i]->BondEnergy();
        m_angle_energy += m_stored_threads[i]->AngleEnergy();
        m_dihedral_energy += m_stored_threads[i]->DihedralEnergy();
        m_inversion_energy += m_stored_threads[i]->InversionEnergy();

        // Claude Generated (2025): Debug thread type
        int thread_type = m_stored_threads[i]->Type();
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param(fmt::format("thread_{}_type", i), std::to_string(thread_type));
        }

        if (thread_type != 3) {
            m_vdw_energy += m_stored_threads[i]->VdWEnergy();
            m_rep_energy += m_stored_threads[i]->RepEnergy();
        } else {
            // GFN-FF (Type == 3) uses different energy components
            // Claude Generated: Store GFN-FF energies in both old and new variables for API compatibility
            h4_energy += m_stored_threads[i]->VdWEnergy();
            m_hh_energy += m_stored_threads[i]->RepEnergy();  // Claude Generated (Session 9): Use member variable
            // CRITICAL FIX (Nov 2025): Do NOT also add to m_rep_energy for GFN-FF!
            // For GFN-FF (method==3), repulsion goes ONLY into m_hh_energy, not m_rep_energy
            // This was causing double-counting: H2 showed 0.266 Eh instead of 0.050 Eh
            // m_rep_energy is for UFF/QMDFF only (method != 3)

            // CRITICAL FIX: Also collect GFN-FF dispersion and Coulomb energies!
            double thread_disp = m_stored_threads[i]->DispersionEnergy();
            double thread_coul = m_stored_threads[i]->CoulombEnergy();
            m_dispersion_energy += thread_disp;
            m_coulomb_energy += thread_coul;

            // Claude Generated (2025): Phase 5 - Collect HB/XB energies
            double thread_hb = m_stored_threads[i]->HydrogenBondEnergy();
            double thread_xb = m_stored_threads[i]->HalogenBondEnergy();
            m_energy_hbond += thread_hb;
            m_energy_xbond += thread_xb;

            // Claude Generated (2025): Debug individual energy components
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::param(fmt::format("thread_{}_dispersion", i),
                    fmt::format("{:.6f} Eh", thread_disp));
                CurcumaLogger::param(fmt::format("thread_{}_coulomb", i),
                    fmt::format("{:.6f} Eh", thread_coul));
                if (thread_hb != 0.0) {
                    CurcumaLogger::param(fmt::format("thread_{}_hbond", i),
                        fmt::format("{:.6f} Eh", thread_hb));
                }
                if (thread_xb != 0.0) {
                    CurcumaLogger::param(fmt::format("thread_{}_xbond", i),
                        fmt::format("{:.6f} Eh", thread_xb));
                }
            }
        }

        m_gradient += m_stored_threads[i]->Gradient();
    }

    // Claude Generated: CG pair interaction calculations (spherical implementation)
    // Only calculated for CG methods (method type 4)
    if (m_method == "cg" || m_method == "cg-lj") {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Calculating CG pair interactions");
        }

        for (const auto& pair : m_vdWs) {
            if (pair.type == 3) { // CG interaction
                Vector3d pos_i = m_geometry.row(pair.i);
                Vector3d pos_j = m_geometry.row(pair.j);

                double pair_energy = CGPotentials::calculateCGPairEnergy(pair, pos_i, pos_j);
                cg_energy += pair_energy;

                if (CurcumaLogger::get_verbosity() >= 3) {
                    CurcumaLogger::param(fmt::format("CG_pair_{}-{}", pair.i, pair.j),
                        fmt::format("{:.6f} Eh", pair_energy));
                }

                // Calculate gradients if requested (numerical differentiation)
                if (gradient) {
                    // TODO: Implement analytical gradients for CG potentials
                    // For now, gradients handled by numerical differentiation in EnergyCalculator
                }
            }
        }

        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::param("total_cg_energy", fmt::format("{:.6f} Eh", cg_energy));
        }
    }

    // Claude Generated: Add GFN-FF dispersion, Coulomb, HB, and XB energies to total
    energy = m_e0 + m_bond_energy + m_angle_energy + m_dihedral_energy + m_inversion_energy + m_vdw_energy + m_rep_energy + m_eq_energy + h4_energy + m_hh_energy + cg_energy + m_dispersion_energy + m_coulomb_energy + m_energy_hbond + m_energy_xbond;  // Claude Generated (Session 9): Use m_hh_energy

    // Claude Generated (2025): Debug total GFN-FF energies
    if (CurcumaLogger::get_verbosity() >= 3 && (m_dispersion_energy != 0.0 || m_coulomb_energy != 0.0)) {
        CurcumaLogger::param("total_gfnff_dispersion", fmt::format("{:.6f} Eh", m_dispersion_energy));
        CurcumaLogger::param("total_gfnff_coulomb", fmt::format("{:.6f} Eh", m_coulomb_energy));
        CurcumaLogger::param("total_before_gfnff", fmt::format("{:.6f} Eh",
            m_e0 + m_bond_energy + m_angle_energy + m_dihedral_energy + m_inversion_energy + m_vdw_energy + m_rep_energy + m_eq_energy + h4_energy + m_hh_energy + cg_energy));  // Claude Generated (Session 9): Use m_hh_energy
    }

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
        CurcumaLogger::param("bond_energy", fmt::format("{:.6f} Eh", m_bond_energy));
        CurcumaLogger::param("angle_energy", fmt::format("{:.6f} Eh", m_angle_energy));
        CurcumaLogger::param("dihedral_energy", fmt::format("{:.6f} Eh", m_dihedral_energy));
        CurcumaLogger::param("inversion_energy", fmt::format("{:.6f} Eh", m_inversion_energy));
        CurcumaLogger::param("nonbonded_energy", fmt::format("{:.6f} Eh", m_vdw_energy + m_rep_energy));
        if (d3_energy != 0.0) {
            CurcumaLogger::param("D3_energy", fmt::format("{:.6f} Eh", d3_energy));
        }
        if (d4_energy != 0.0) {
            CurcumaLogger::param("D4_energy", fmt::format("{:.6f} Eh", d4_energy));
        }
        if (h4_energy != 0.0) {
            CurcumaLogger::param("HBond_correction", fmt::format("{:.6f} Eh", h4_energy));
        }
        if (m_hh_energy != 0.0) {  // Claude Generated (Session 9): Use m_hh_energy
            CurcumaLogger::param("HH_repulsion", fmt::format("{:.6f} Eh", m_hh_energy));
        }
        if (m_dispersion_energy != 0.0) {
            CurcumaLogger::param("GFNFF_dispersion", fmt::format("{:.6f} Eh", m_dispersion_energy));
        }
        if (m_coulomb_energy != 0.0) {
            CurcumaLogger::param("GFNFF_coulomb", fmt::format("{:.6f} Eh", m_coulomb_energy));
        }
        if (m_energy_hbond != 0.0) {
            CurcumaLogger::param("GFNFF_hydrogen_bond", fmt::format("{:.6f} Eh", m_energy_hbond));
        }
        if (m_energy_xbond != 0.0) {
            CurcumaLogger::param("GFNFF_halogen_bond", fmt::format("{:.6f} Eh", m_energy_xbond));
        }
        if (cg_energy != 0.0) {
            CurcumaLogger::param("CG_interactions", fmt::format("{:.6f} Eh", cg_energy));
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
        int gfnff_dispersions = 0, gfnff_repulsions = 0, gfnff_coulombs = 0;

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

        // Claude Generated (2025): Check for GFN-FF specific pairwise parameters
        if (m_parameters.contains("gfnff_dispersions") && m_parameters["gfnff_dispersions"].is_array()) {
            gfnff_dispersions = m_parameters["gfnff_dispersions"].size();
        }
        if (m_parameters.contains("gfnff_repulsions") && m_parameters["gfnff_repulsions"].is_array()) {
            gfnff_repulsions = m_parameters["gfnff_repulsions"].size();
        }
        if (m_parameters.contains("gfnff_coulombs") && m_parameters["gfnff_coulombs"].is_array()) {
            gfnff_coulombs = m_parameters["gfnff_coulombs"].size();
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

            // Print GFN-FF specific parameters if present
            if (gfnff_dispersions > 0 || gfnff_repulsions > 0 || gfnff_coulombs > 0) {
                CurcumaLogger::info("GFN-FF pairwise interactions:");
                if (gfnff_dispersions > 0) CurcumaLogger::param("dispersion_pairs", fmt::format("{}", gfnff_dispersions));
                if (gfnff_repulsions > 0) CurcumaLogger::param("repulsion_pairs", fmt::format("{}", gfnff_repulsions));
                if (gfnff_coulombs > 0) CurcumaLogger::param("coulomb_pairs", fmt::format("{}", gfnff_coulombs));
            }

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
