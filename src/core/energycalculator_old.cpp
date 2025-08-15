/*
 * < General Energy and Gradient Calculator >
 * Copyright (C) 2022 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#ifdef USE_TBLITE
#include "src/core/tbliteinterface.h"
#endif

#ifdef USE_XTB
#include "src/core/xtbinterface.h"
#endif

#ifndef _WIN32
#include <filesystem>
namespace fs = std::filesystem;
#endif

#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>

#include <fmt/core.h>
#include <fmt/format.h>

#include "forcefieldgenerator.h"
#include "src/core/elements.h"

#include "energycalculator.h"

EnergyCalculator::EnergyCalculator(const std::string& method, const json& controller)
    : m_method(method)
    , m_basename("")
{
    initializeCommon(controller);
}

// Claude Generated: Constructor with basename for parameter caching
EnergyCalculator::EnergyCalculator(const std::string& method, const json& controller, const std::string& basename)
    : m_method(method)
    , m_basename(basename)
{
    initializeCommon(controller);
}

void EnergyCalculator::initializeCommon(const json& controller)
{
    std::transform(m_method.begin(), m_method.end(), m_method.begin(), [](unsigned char c) { return std::tolower(c); });

    m_controller = MergeJson(EnergyCalculatorJson ,controller);

    if (controller.contains("param_file")) {
        m_param_file = controller["param_file"];
    }

    if (controller.contains("write_param")) {
        m_writeparam = controller["write_param"];
    }

    if (controller.contains("geometry_file")) {
        m_geometry_file = controller["geometry_file"];
    }

    m_bonds = []() {
        return std::vector<std::vector<double>>{ {} };
    };

    m_mult = m_controller["multi"];
    m_SCFmaxiter = m_controller["SCFmaxiter"];
    m_solvent = m_controller["solvent"];
    m_Tele = m_controller["Tele"];

    switch (SwitchMethod(m_method)) {
    case 8:
        /*
             m_qmdff = new QMDFF(controller);
             if (m_parameter.size())
                 m_qmdff->setParameter(m_parameter);
             m_ecengine = [this](bool gradient, bool verbose) {
                 this->CalculateQMDFF(gradient, verbose);
             };
             */
        break;
    case 7:
        /*
            m_uff = new eigenUFF(controller);
            m_ecengine = [this](bool gradient, bool verbose) {
                this->CalculateUFF(gradient, verbose);
            };*/
        break;

    case 6:
        m_qminterface = new EHT();
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateQMInterface(gradient, verbose);
        };
        break;

    case 5:
        // m_d4 = new DFTD4Interface(controller);
        //  m_ecengine = [this](bool gradient, bool verbose) {
        //      this->CalculateD4(gradient, verbose);
        //  };
        break;

    case 4:
#ifdef USE_D3
        m_qminterface = new DFTD3Interface(controller);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateD3(gradient, verbose);
        };
#else
        std::cout << "D3 was not included ..." << std::endl;
        exit(1);
#endif

        break;

    case 3:
#ifdef USE_ULYSSES
        m_qminterface = new UlyssesInterface(controller);
        m_qminterface->setMult(m_mult);

        m_qminterface->setMethod(m_method);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateUlysses(gradient, verbose);
        };
#else
        std::cout << "Ulysses was not included ..." << std::endl;
        exit(1);
#endif
        break;

    case 2:
#ifdef USE_XTB
        m_qminterface = new XTBInterface(m_controller);
        m_qminterface->setMult(m_mult);
        m_qminterface->setMethod(m_method);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateXTB(gradient, verbose);
        };
#else
        std::cout << "XTB was not included ..." << std::endl;
        exit(1);
#endif
    break;

    case 1:
#ifdef USE_TBLITE
        m_qminterface = new TBLiteInterface(m_controller);
        m_qminterface->setMult(m_mult);

        m_qminterface->setMethod(m_method);

        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateTBlite(gradient, verbose);
            m_error = this->m_qminterface->Error();
        };
#else
        std::cout << "TBLite was not included ..." << std::endl;
        exit(1);
#endif
        /*
        m_bonds = [this]() {
            return this->m_tblite->BondOrders();
        };
      */
        break;

    case 0:
    default:
        m_forcefield = new ForceField(controller);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateFF(gradient, verbose);
        };
        break;
    }
}

EnergyCalculator::~EnergyCalculator()
{
    // Claude Generated: Use centralized dispatch for cleanup
    dispatchMethodAction("cleanup");
    // delete[] m_coord;
    // delete[] m_grad;
}

void EnergyCalculator::setMolecule(const Mol& mol)
{
    m_mol = mol;
    m_atoms = mol.m_number_atoms;
    m_geometry = mol.m_geometry;
    m_gradient = Eigen::MatrixXd::Zero(m_atoms, 3);
    m_xtb_gradient = Eigen::VectorXd::Zero(3 * m_atoms);

    /*
    std::vector<int> atoms = molecule.Atoms();
    m_coord = new double[3 * m_atoms];
    m_grad = new double[3 * m_atoms];
    m_gradient = Eigen::MatrixXd::Zero(m_atoms, 3);
    m_geometry = molecule.getGeometry();
    */

    switch (SwitchMethod(m_method)) {
    case 8:
        // m_qmdff->setMolecule(m_atoms, m_geometry);
        // m_qmdff->Initialise();
        break;
    case 7:

        break;

    case 6:
    case 3:
    case 2:
    case 1:
        // Claude Generated: Use centralized dispatch for standard QM methods
        dispatchMethodAction("setmolecule", &mol);
        if (SwitchMethod(m_method) == 2 || SwitchMethod(m_method) == 1) {
            m_qminterface->setMethod(m_method);
        }
        break;

    case 5:
#ifdef USE_D4
        static_cast<DFTD4Interface*>(m_qminterface)->InitialiseMolecule(mol, 1 / au);
#endif
        break;

    case 4:
#ifdef USE_D3
        static_cast<DFTD3Interface*>(m_qminterface)->InitialiseMolecule(mol.m_atoms);
#endif
        break;

    case 0:
    default:
        if (m_parameter.size() == 0) {
            // Claude Generated: Try loading cached parameters first
            bool loaded_from_cache = tryLoadAutoParameters(mol);

            if (!loaded_from_cache) {
                if (!std::filesystem::exists(m_param_file)) {
                    ForceFieldGenerator ff(m_controller);
                    ff.setMolecule(mol);
                    ff.Generate();
                    m_parameter = ff.getParameter();

                    // Auto-save parameters with standardized naming
                    saveAutoParameters(mol, m_parameter);

                    // if (m_writeparam) {
                    //     std::ofstream parameterfile("ff_param.json");
                    //     parameterfile << m_parameter;
                    // }
                }

                // Claude Generated: Always show parameter analysis after loading/generating

            } // end if (!loaded_from_cache)
        } // end if (m_parameter.size() == 0)
        // Claude Generated: Use centralized dispatch for ForceField
        dispatchMethodAction("setmolecule", &mol);
        m_forcefield->setParameter(m_parameter);
        break;
    }
    m_initialised = true;
}

int EnergyCalculator::SwitchMethod(const std::string& method)
{
    // Claude Generated: Helper function to get module name from case number
    auto getModuleName = [](int case_num) -> std::string {
        switch (case_num) {
        case 0:
            return "ForceField";
        case 1:
            return "TBLite";
        case 2:
            return "XTB";
        case 3:
            return "Ulysses";
        case 4:
            return "DFT-D3";
        case 5:
            return "DFT-D4";
        case 6:
            return "EHT";
        case 7:
            return "QMDFF (legacy)";
        case 8:
            return "UFF (legacy)";
        default:
            return "Unknown";
        }
    };

    // Claude Generated: Debug logger to show available modules
    static bool logged = false;
    if (!logged) {
        std::cout << "=== Available Curcuma Modules ===" << std::endl;
        std::cout << "ForceField: YES (always available)" << std::endl;
#ifdef USE_TBLITE
        std::cout << "TBLite: YES" << std::endl;
#else
        std::cout << "TBLite: NO" << std::endl;
#endif
#ifdef USE_XTB
        std::cout << "XTB: YES" << std::endl;
#else
        std::cout << "XTB: NO" << std::endl;
#endif
#ifdef USE_ULYSSES
        std::cout << "Ulysses: YES" << std::endl;
#else
        std::cout << "Ulysses: NO" << std::endl;
#endif
#ifdef USE_D3
        std::cout << "DFT-D3: YES" << std::endl;
#else
        std::cout << "DFT-D3: NO" << std::endl;
#endif
#ifdef USE_D4
        std::cout << "DFT-D4: YES" << std::endl;
#else
        std::cout << "DFT-D4: NO" << std::endl;
#endif
        std::cout << "EHT: YES (always available)" << std::endl;
        std::cout << "===================================" << std::endl;
        logged = true;
    }

    // Claude Generated: Lookup-table driven method resolution with priorities and fallbacks

    // 1. Explicit method mappings (highest priority, no fallbacks)
    static const std::map<std::string, int> explicit_methods = {
        { "cgfnff", 0 }, // Native GFN-FF
        { "ugfn2", 3 }, // Explicit Ulysses GFN2
        { "eht", 6 }, // Extended Hückel Theory
        { "d3", 4 }, // DFT-D3
        { "d4", 5 } // DFT-D4
    };

    // 2. Shared methods with priority fallbacks (TBLite > Ulysses > XTB)
    struct MethodPriority {
        std::string method_name;
        std::vector<std::pair<int, std::string>> priorities; // {case_number, compilation_flag}
    };

    static const std::vector<MethodPriority> shared_methods = {
        { "gfn2", { { 1, "USE_TBLITE" }, { 3, "USE_ULYSSES" }, { 2, "USE_XTB" } } },
        { "gfn1", { { 1, "USE_TBLITE" }, { 2, "USE_XTB" } } },
        { "ipea1", { { 1, "USE_TBLITE" } } },
        { "gfnff", { { 0, "" } } } // Alias for cgfnff → ForceField
    };

    // 3. Library-specific methods (exact match from lists)
    struct LibraryMethods {
        const StringList* method_list;
        int case_number;
        std::string compilation_flag;
    };

    static const std::vector<LibraryMethods> library_methods = {
        { &m_ff_methods, 0, "" },
        { &m_ulysses_methods, 3, "USE_ULYSSES" },
        { &m_xtb_methods, 2, "USE_XTB" },
        { &m_qmdff_method, 7, "" },
        { &m_uff_methods, 8, "" } // Legacy
    };

    // Method resolution logic:

    // Check explicit methods first
    auto explicit_it = explicit_methods.find(method);
    if (explicit_it != explicit_methods.end()) {
        std::cout << "Method '" << method << "' resolved to: " << getModuleName(explicit_it->second) << " (explicit mapping)" << std::endl;
        return explicit_it->second;
    }

    // Check shared methods with priority fallbacks
    for (const auto& shared : shared_methods) {
        if (method == shared.method_name) {
            for (const auto& priority : shared.priorities) {
                if (priority.second.empty() || isCompiled(priority.second)) {
                    std::cout << "Method '" << method << "' resolved to: " << getModuleName(priority.first) << " (priority fallback)" << std::endl;
                    return priority.first;
                }
            }
        }
    }

    // Check library-specific methods
    for (const auto& lib : library_methods) {
        for (const auto& lib_method : *lib.method_list) {
            if (method.find(lib_method) != std::string::npos) {
                if (lib.compilation_flag.empty() || isCompiled(lib.compilation_flag)) {
                    std::cout << "Method '" << method << "' resolved to: " << getModuleName(lib.case_number) << " (library-specific)" << std::endl;
                    return lib.case_number;
                } else {
                    std::cerr << "Library not available for " << method << std::endl;
                }
            }
        }
    }

    // Default fallback
    std::cerr << "Unknown method: " << method << ", using default ForceField" << std::endl;
    std::cout << "Method '" << method << "' resolved to: " << getModuleName(0) << " (default fallback)" << std::endl;
    return 0;
}

// Claude Generated: Helper function to check compilation flags
bool EnergyCalculator::isCompiled(const std::string& flag) const
{
    if (flag == "USE_TBLITE") {
#ifdef USE_TBLITE
        return true;
#else
        return false;
#endif
    }
    if (flag == "USE_ULYSSES") {
#ifdef USE_ULYSSES
        return true;
#else
        return false;
#endif
    }
    if (flag == "USE_XTB") {
#ifdef USE_XTB
        return true;
#else
        return false;
#endif
    }
    // Add other flags as needed
    return true; // Empty flag means always available
}

// Claude Generated: Centralized method dispatch to reduce repetitive switch statements
void EnergyCalculator::dispatchMethodAction(const std::string& action, const Mol* mol)
{
    int method_case = SwitchMethod(m_method);

    if (action == "cleanup") {
        // Destructor logic - centralized cleanup
        switch (method_case) {
        case 6:
        case 4:
        case 3:
        case 2:
        case 1:
            delete m_qminterface;
            break;
        case 0:
            delete m_forcefield;
            break;
            // cases 8, 7, 5 are commented out legacy methods
        }
    } else if (action == "setmolecule" && mol != nullptr) {
        // setMolecule logic - centralized molecule setting
        switch (method_case) {
        case 6:
        case 4:
        case 3:
        case 2:
        case 1:
            m_qminterface->InitialiseMolecule(*mol);
            break;
        case 0:
            m_forcefield->setMolecule(*mol);
            break;
            // cases 8, 7, 5 are commented out legacy methods
        }
    }
    // Note: Constructor logic remains in switch statement due to complexity
}
void EnergyCalculator::updateGeometry(const Eigen::VectorXd& geometry)
{
#pragma message("Eigen::VectorXd ....")
    for (int i = 0; i < m_atoms; ++i) {
        m_geometry(i, 0) = geometry[3 * i + 0];
        m_geometry(i, 1) = geometry[3 * i + 1];
        m_geometry(i, 2) = geometry[3 * i + 2];
    }

    // m_containsNaN = std::isnan(m_geometry[m_atoms - 1][0]);
}

void EnergyCalculator::updateGeometry(const double* coord)
{

    for (int i = 0; i < m_atoms; ++i) {
        m_geometry(i, 0) = coord[3 * i + 0];
        m_geometry(i, 1) = coord[3 * i + 1];
        m_geometry(i, 2) = coord[3 * i + 2];
    }
    // m_containsNaN = std::isnan(m_geometry[m_atoms - 1][0]);
}

void EnergyCalculator::updateGeometry(const std::vector<double>& geometry)
{
    for (int i = 0; i < m_atoms; ++i) {
        m_geometry(i, 0) = geometry[3 * i + 0];
        m_geometry(i, 1) = geometry[3 * i + 1];
        m_geometry(i, 2) = geometry[3 * i + 2];
    }
    // m_containsNaN = std::isnan(m_geometry[m_atoms - 1][0]);
}

void EnergyCalculator::updateGeometry(const Matrix& geometry)
{
    m_geometry = geometry;
}

double EnergyCalculator::CalculateEnergy(bool gradient, bool verbose)
{
    m_ecengine(gradient, verbose);
    return m_energy;
}

Eigen::MatrixXd EnergyCalculator::NumGrad()
{
    Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(m_atoms, 3);

    double dx = 1e-4; // m_d;
    double E1, E2;
    for (int i = 0; i < m_atoms; ++i) {
        for (int j = 0; j < 3; ++j) {
            m_geometry(i, j) += dx;
            E1 = CalculateEnergy(false, false);
            m_geometry(i, j) -= 2 * dx;
            E2 = CalculateEnergy(false, false);
            //   std::cout << E1 << " " << E2 << " " << m_energy << std::endl;
            gradient(i, j) = (E1 - E2) / (2 * dx);
            m_geometry(i, j) += dx;
        }
    }
    return gradient;
}

void EnergyCalculator::CalculateUFF(bool gradient, bool verbose)
{
    /*
    m_uff->UpdateGeometry(m_geometry);
    m_energy = m_uff->Calculate(gradient, verbose);
    if (gradient) {
        m_gradient = m_uff->Gradient();
        // m_gradient = m_uff->NumGrad();
    }
    */
}

void EnergyCalculator::CalculateTBlite(bool gradient, bool verbose)
{
#ifdef USE_TBLITE
    m_qminterface->UpdateMolecule(m_geometry);
    m_energy = m_qminterface->Calculation(gradient, verbose);
    if (gradient)
        m_gradient = m_qminterface->Gradient();

#endif
}

void EnergyCalculator::CalculateXTB(bool gradient, bool verbose)
{
#ifdef USE_XTB
    m_qminterface->UpdateMolecule(m_geometry);
    m_energy = m_qminterface->Calculation(gradient, verbose);
    if (gradient)
        m_gradient = m_qminterface->Gradient();

#endif
}

void EnergyCalculator::CalculateUlysses(bool gradient, bool verbose)
{
    m_qminterface->UpdateMolecule(m_geometry);
    m_energy = m_qminterface->Calculation(gradient, verbose);
    if (gradient) {
        m_gradient = m_qminterface->Gradient();
    }
}

void EnergyCalculator::CalculateQMInterface(bool gradient, bool verbose)
{
    m_qminterface->UpdateMolecule(m_geometry);
    m_energy = m_qminterface->Calculation(gradient, verbose);
    if (gradient) {
        if (m_qminterface->hasGradient())
            m_gradient = m_qminterface->Gradient();
        else
            m_gradient = NumGrad();
    }
}

void EnergyCalculator::CalculateD3(bool gradient, bool verbose)
{
#ifdef USE_D3
    for (int i = 0; i < m_atoms; ++i) {
        static_cast<DFTD3Interface*>(m_qminterface)->UpdateAtom(i, m_geometry(i, 0), m_geometry(i, 1), m_geometry(i, 2));
    }
    if (gradient) {
        m_energy = m_qminterface->Calculation(m_xtb_gradient.data());
        for (int i = 0; i < m_atoms; ++i) {
            m_gradient(i, 0) = m_xtb_gradient[3 * i + 0] * au;
            m_gradient(i, 1) = m_xtb_gradient[3 * i + 1] * au;
            m_gradient(i, 2) = m_xtb_gradient[3 * i + 2] * au;
        }
    } else
        m_energy = m_qminterface->Calculation(false, verbose);
#endif
}

void EnergyCalculator::CalculateD4(bool gradient, bool verbose)
{
#ifdef USE_D4
    for (int i = 0; i < m_atoms; ++i) {
        static_cast<DFTD4Interface*>(m_qminterface)->UpdateAtom(i, m_geometry(i, 0) / au, m_geometry(i, 1) / au, m_geometry(i, 2) / au);
    }
    if (gradient) {
        m_energy = m_qminterface->Calculation(m_grad);
        for (int i = 0; i < m_atoms; ++i) {
            m_gradient(i, 0) = m_grad[3 * i + 0] * au;
            m_gradient(i, 1) = m_grad[3 * i + 1] * au;
            m_gradient(i, 2) = m_grad[3 * i + 2] * au;
        }
    } else
        m_energy = m_qminterface->Calculation();
#endif
}

void EnergyCalculator::CalculateQMDFF(bool gradient, bool verbose)
{
    /*
    m_qmdff->UpdateGeometry(m_geometry);
    m_energy = m_qmdff->Calculate(gradient, verbose);
    if (gradient) {
        m_gradient = m_qmdff->Gradient();
    }*/
}

void EnergyCalculator::CalculateFF(bool gradient, bool verbose)
{
    m_forcefield->UpdateGeometry(m_geometry);
    m_energy = m_forcefield->Calculate(gradient, verbose);
    if (gradient) {
        m_gradient = m_forcefield->Gradient();
    }
}

Vector EnergyCalculator::Charges() const
{
    if (m_qminterface == nullptr)
        return Vector{};
    return m_qminterface->Charges();
}

Position EnergyCalculator::Dipole() const
{
    if (m_qminterface == nullptr)
        return Position{};
    return m_qminterface->Dipole();
}

std::vector<std::vector<double>> EnergyCalculator::BondOrders() const
{
    return m_bonds();
}

// Claude Generated: Auto-save force field parameters with intelligent naming
void EnergyCalculator::saveAutoParameters(const Mol& mol, const json& parameters)
{
    std::string param_filename;

    // Priority 1: Use geometry file from controller
    if (!m_geometry_file.empty()) {
        param_filename = ForceField::generateParameterFileName(m_geometry_file);
    }
    // Priority 2: Use basename from CurcumaMethod
    else if (!m_basename.empty()) {
        param_filename = m_basename + ".param.json";
    }
    // Priority 3: Fallback to molecular formula
    else {
        std::string formula = mol.m_formula;
        param_filename = fmt::format("{}_{}_{}.param.json", m_method, formula, mol.m_number_atoms);
    }

    try {
        std::ofstream param_file(param_filename);
        if (param_file.is_open()) {
            // Add metadata for validation
            json output_params = parameters;
            output_params["generated_by"] = "curcuma_energycalculator";
            output_params["timestamp"] = std::chrono::system_clock::now().time_since_epoch().count();
            output_params["atoms"] = mol.m_number_atoms;

            param_file << std::setw(2) << output_params << std::endl;
            param_file.close();

            fmt::print("Auto-saved force field parameters to: {}\n", param_filename);
        } else {
            fmt::print("Warning: Could not save parameters to {}\n", param_filename);
        }
    } catch (const std::exception& e) {
        fmt::print("Error saving parameters: {}\n", e.what());
    }
}

// Claude Generated: Try loading cached force field parameters with intelligent naming
bool EnergyCalculator::tryLoadAutoParameters(const Mol& mol)
{
    std::string param_filename;

    // Priority 1: Use geometry file from controller
    if (!m_geometry_file.empty()) {
        param_filename = ForceField::generateParameterFileName(m_geometry_file);
    }
    // Priority 2: Use basename from CurcumaMethod
    else if (!m_basename.empty()) {
        param_filename = m_basename + ".param.json";
    }
    // Priority 3: Fallback to molecular formula
    else {
        std::string formula = mol.m_formula;
        param_filename = fmt::format("{}_{}_{}.param.json", m_method, formula, mol.m_number_atoms);
    }

    // Check if parameter file exists
    if (!std::filesystem::exists(param_filename)) {
        return false;
    }

    try {
        std::ifstream param_file(param_filename);
        if (!param_file.is_open()) {
            return false;
        }

        json loaded_params;
        param_file >> loaded_params;
        param_file.close();

        // Validate parameters
        if (!loaded_params.contains("method") || !loaded_params.contains("atoms")) {
            fmt::print("Warning: Invalid cached parameter file: {}\n", param_filename);
            return false;
        }

        // Check method compatibility
        std::string cached_method = loaded_params["method"];
        if (cached_method != m_method) {
            fmt::print("Method mismatch in cached parameters (found: {}, expected: {})\n",
                cached_method, m_method);
            return false;
        }

        // Check atom count compatibility
        int cached_atoms = loaded_params["atoms"];
        if (cached_atoms != mol.m_number_atoms) {
            fmt::print("Atom count mismatch in cached parameters (found: {}, expected: {})\n",
                cached_atoms, mol.m_number_atoms);
            return false;
        }

        // Parameters are valid, use them
        m_parameter = loaded_params;

        fmt::print("Loaded cached {} parameters from: {}\n", m_method, param_filename);
        return true;

    } catch (const std::exception& e) {
        fmt::print("Error loading cached parameters from {}: {}\n", param_filename, e.what());
        return false;
    }
}
