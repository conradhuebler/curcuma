//
// Created by gerd on 11.12.24.
// Extended by Claude for ComputationalMethod integration (June 2026)
//

#include "orcainterface.h"
#include "src/core/citation_registry.h"
#include "src/core/curcuma_logger.h"
#include "src/core/units.h"

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <filesystem>
#include <iomanip>

#ifndef _WIN32
#include <sys/wait.h>
#endif

namespace fs = std::filesystem;

// Simple element symbol lookup for ORCA input generation.
// Claude Generated - June 2026
static const char* const kElementSymbols[] = {
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc",
    "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc",
    "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
    "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
};

static std::string elementSymbol(int z)
{
    if (z >= 1 && z <= 118)
        return kElementSymbols[z];
    return "X";
}

OrcaInterface::OrcaInterface(const ConfigManager& config)
    : m_config(config)
{
    inputFilePath = "orca.inp";
    outputFilePath = "orca.out";
    m_input_path = m_config.get<std::string>("orca_basename", "orca_calc") + ".inp";
    m_output_path = m_config.get<std::string>("orca_basename", "orca_calc") + ".out";
}

OrcaInterface::OrcaInterface()
    : m_config("orca", json{})
{
    inputFilePath = "orca.inp";
    outputFilePath = "orca.out";
}

OrcaInterface::~OrcaInterface()
{
    // Optional cleanup
}

void OrcaInterface::setInputFile(const std::string& inputFile)
{
    inputFilePath = inputFile;
    m_input_path = inputFile;
}

bool OrcaInterface::createInputFile(const std::string& content)
{
    std::ofstream outFile(inputFilePath);
    if (!outFile) {
        std::cerr << "Fehler beim Erstellen der Eingabedatei!" << std::endl;
        return false;
    }
    outFile << content;
    outFile.close();
    return true;
}

bool OrcaInterface::executeOrcaProcess()
{
    std::stringstream command;
    command << "orca " << inputFilePath << " > " << outputFilePath;
    int result = std::system(command.str().c_str());
    return (result == 0);
}

bool OrcaInterface::runOrca()
{
    CitationRegistry::cite("orca");
    std::cout << "Starte ORCA..." << std::endl;
    if (executeOrcaProcess()) {
        std::cout << "ORCA abgeschlossen!" << std::endl;
        return true;
    } else {
        std::cerr << "Fehler beim Ausführen von ORCA!" << std::endl;
        return false;
    }
}

void OrcaInterface::readOrcaJSON()
{
    std::ifstream property(inputFilePath + ".property.json");
    property >> OrcaJSON;
}

bool OrcaInterface::getOrcaJSON()
{
    std::stringstream command;
    command << "orca_2json " << inputFilePath << " -property >> " << outputFilePath;
    const int result = std::system(command.str().c_str());
    return (result == 0);
}

// =================================================================================
// ComputationalMethod support (Claude Generated - June 2026)
// =================================================================================

std::string OrcaInterface::findOrcaExecutable() const
{
    // Priority 1: Config parameter orca_executable
    std::string config_path = m_config.get<std::string>("orca_executable", "");
    if (!config_path.empty() && fs::exists(config_path)) {
        return config_path;
    }

    // Priority 2: Environment variable ORCA_PATH (full path to binary)
    const char* orca_path_env = std::getenv("ORCA_PATH");
    if (orca_path_env && fs::exists(orca_path_env)) {
        return std::string(orca_path_env);
    }

    // Priority 3: Environment variable ORCA_ROOT (directory, append /orca)
    const char* orca_root_env = std::getenv("ORCA_ROOT");
    if (orca_root_env) {
        std::string root_bin = std::string(orca_root_env) + "/orca";
        if (fs::exists(root_bin)) {
            return root_bin;
        }
    }

    // Priority 4: Fallback to "orca" via system PATH.
    // ORCA parallel runs require the ABSOLUTE path, so we resolve it here.
#ifndef _WIN32
    FILE* which_pipe = popen("which orca 2> /dev/null", "r");
    if (which_pipe) {
        char path_buf[1024];
        if (fgets(path_buf, sizeof(path_buf), which_pipe) != nullptr) {
            std::string path(path_buf);
            // Trim trailing newline
            if (!path.empty() && path.back() == '\n') path.pop_back();
            pclose(which_pipe);
            if (!path.empty() && fs::exists(path)) {
                return path;
            }
        } else {
            pclose(which_pipe);
        }
    }
#endif

    return "";
}

bool OrcaInterface::checkOrcaExecutable()
{
    // Quick static check without config — sufficient for isAvailable()
#ifdef _WIN32
    int result = std::system("where orca > NUL 2>>1");
#else
    int result = std::system("which orca > /dev/null 2>>1");
#endif
    if (result == 0) return true;

    // Also check environment variables (no config needed)
    const char* orca_path_env = std::getenv("ORCA_PATH");
    if (orca_path_env && fs::exists(orca_path_env)) return true;

    const char* orca_root_env = std::getenv("ORCA_ROOT");
    if (orca_root_env && fs::exists(std::string(orca_root_env) + "/orca")) return true;

    return false;
}

std::string OrcaInterface::buildGeometryBlock(const Mol& mol) const
{
    const int natoms = mol.m_number_atoms;
    const int charge = mol.m_charge;
    // Multiplicity: 2S+1. Default singlet (1); use mol.m_spin if non-zero.
    int multiplicity = 1;
    if (mol.m_spin > 0) {
        multiplicity = static_cast<int>(2.0 * mol.m_spin + 1.0 + 0.5);
    }

    std::ostringstream block;
    block << "*xyz " << charge << " " << multiplicity << "\n";

    for (int i = 0; i < natoms; ++i) {
        int element = mol.m_atoms[i];
        block << "  " << elementSymbol(element) << "  "
              << std::fixed << std::setprecision(8)
              << mol.m_geometry(i, 0) << "  "
              << mol.m_geometry(i, 1) << "  "
              << mol.m_geometry(i, 2) << "\n";
    }
    block << "*\n";
    return block.str();
}

bool OrcaInterface::generateInput(const Mol& mol, const std::string& method_keyword)
{
    clearError();

    // Determine paths from config
    std::string basename = m_config.get<std::string>("orca_basename", "orca_calc");
    m_input_path = basename + ".inp";
    m_output_path = basename + ".out";

    // Determine method keyword
    std::string method = method_keyword.empty()
        ? m_config.get<std::string>("orca_method", "HF-3c")
        : method_keyword;

    std::string custom_kw = m_config.get<std::string>("orca_keywords", "");
    std::string extra_kw = m_config.get<std::string>("orca_extra_keywords", "EnGrad TightSCF");
    int nprocs = m_config.get<int>("orca_nprocs", 1);
    int maxcore = m_config.get<int>("orca_maxcore", 500);

    std::ostringstream content;

    if (!custom_kw.empty()) {
        // Custom keywords mode: user supplies everything after !
        content << "! " << custom_kw << "\n";
    } else {
        // Auto-generated input
        content << "! " << method;
        if (!extra_kw.empty())
            content << " " << extra_kw;
        content << "\n";
    }

    // Parallel settings
    if (nprocs > 1) {
        content << "%pal nprocs " << nprocs << " end\n";
    }
    if (maxcore > 0) {
        content << "%maxcore " << maxcore << "\n";
    }

    // Solvation
    std::string cpcm_solvent = m_config.get<std::string>("cpcm_solvent", "none");
    if (cpcm_solvent != "none" && !cpcm_solvent.empty()) {
        content << "%cpcm\n  solvent \"" << cpcm_solvent << "\"\nend\n";
    }

    std::string alpb_solvent = m_config.get<std::string>("alpb_solvent", "none");
    if (alpb_solvent != "none" && !alpb_solvent.empty()) {
        content << "%alpb\n  solvent \"" << alpb_solvent << "\"\nend\n";
    }

    content << "\n";
    content << buildGeometryBlock(mol);

    std::ofstream out(m_input_path);
    if (!out) {
        m_has_error = true;
        m_error_message = "Failed to write ORCA input: " + m_input_path;
        return false;
    }
    out << content.str();
    out.close();

    // Update original paths for backward compatibility
    inputFilePath = m_input_path;
    outputFilePath = m_output_path;

    return true;
}

bool OrcaInterface::updateGeometryInInput(const Mol& mol)
{
    clearError();

    std::ifstream in(m_input_path);
    if (!in) {
        m_has_error = true;
        m_error_message = "ORCA input file not found for geometry update: " + m_input_path;
        return false;
    }

    std::string content((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());
    in.close();

    // Find and replace geometry block
    size_t xyz_pos = content.find("*xyz");
    if (xyz_pos == std::string::npos)
        xyz_pos = content.find("* XYZ");

    if (xyz_pos == std::string::npos) {
        m_has_error = true;
        m_error_message = "No *xyz block found in ORCA input: " + m_input_path;
        return false;
    }

    size_t end_pos = content.find("*\n", xyz_pos + 1);
    if (end_pos == std::string::npos)
        end_pos = content.find("* \n", xyz_pos + 1);

    if (end_pos == std::string::npos) {
        m_has_error = true;
        m_error_message = "No closing * found for xyz block in: " + m_input_path;
        return false;
    }

    std::string new_content = content.substr(0, xyz_pos);
    new_content += buildGeometryBlock(mol);
    new_content += content.substr(end_pos + 2);

    std::ofstream out(m_input_path);
    if (!out) {
        m_has_error = true;
        m_error_message = "Failed to write updated ORCA input: " + m_input_path;
        return false;
    }
    out << new_content;
    out.close();

    return true;
}

bool OrcaInterface::produceAndLoadJSON()
{
    // Run orca_2json to produce .property.json
    std::string json_cmd = "orca_2json " + m_input_path + " -property > /dev/null 2>>1";
    int result = std::system(json_cmd.c_str());
    if (result != 0) {
        return false;
    }
    return loadPropertyJSON();
}

bool OrcaInterface::loadPropertyJSON()
{
    std::string json_path = m_input_path + ".property.json";
    std::ifstream in(json_path);
    if (!in) {
        return false;
    }
    try {
        in >> m_property_json;
        return true;
    } catch (...) {
        return false;
    }
}

double OrcaInterface::parseEnergy() const
{
    // Try JSON first
    if (!m_property_json.is_null() && m_property_json.contains("Calculation")) {
        try {
            auto calc = m_property_json["Calculation"];
            if (calc.contains("CalculationType") && calc["CalculationType"].contains("Energy")) {
                auto energy = calc["CalculationType"]["Energy"];
                if (energy.contains("TotalEnergy") && energy["TotalEnergy"].contains("Value")) {
                    return energy["TotalEnergy"]["Value"].get<double>();
                }
            }
        } catch (...) {
            // Fall through to text parsing
        }
    }

    // Text fallback: search for "FINAL SINGLE POINT ENERGY"
    std::ifstream out(m_output_path);
    if (!out) return 0.0;

    std::string line;
    while (std::getline(out, line)) {
        size_t pos = line.find("FINAL SINGLE POINT ENERGY");
        if (pos != std::string::npos) {
            std::istringstream iss(line.substr(pos + 25));
            double energy = 0.0;
            if (iss >> energy) {
                return energy;
            }
        }
    }
    return 0.0;
}

Matrix OrcaInterface::parseGradient(int natoms) const
{
    // Try JSON first
    if (!m_property_json.is_null() && m_property_json.contains("Calculation")) {
        try {
            auto calc = m_property_json["Calculation"];
            if (calc.contains("CalculationType") &&
                calc["CalculationType"].contains("Properties") &&
                calc["CalculationType"]["Properties"].contains("Gradient")) {
                auto grad = calc["CalculationType"]["Properties"]["Gradient"];
                if (grad.is_array() && grad.size() == static_cast<size_t>(natoms * 3)) {
                    Matrix g(natoms, 3);
                    for (int i = 0; i < natoms; ++i) {
                        for (int j = 0; j < 3; ++j) {
                            g(i, j) = grad[i * 3 + j].get<double>();
                        }
                    }
                    return g;
                }
            }
        } catch (...) {
            // Fall through to text parsing
        }
    }

    return parseGradientFromText(natoms);
}

Matrix OrcaInterface::parseGradientFromText(int natoms) const
{
    std::ifstream out(m_output_path);
    if (!out) return Matrix::Zero(natoms, 3);

    std::string line;
    bool found_gradient = false;

    while (std::getline(out, line)) {
        // Look for "CARTESIAN GRADIENT" or "The cartesian gradient:"
        if (line.find("The cartesian gradient:") != std::string::npos ||
            line.find("CARTESIAN GRADIENT") != std::string::npos) {
            found_gradient = true;
            break;
        }
    }

    if (!found_gradient) {
        return Matrix::Zero(natoms, 3);
    }

    // Skip blank lines
    while (std::getline(out, line)) {
        if (!line.empty()) break;
    }

    // Read natoms lines of gradient data
    Matrix g(natoms, 3);
    for (int i = 0; i < natoms; ++i) {
        if (!std::getline(out, line)) break;

        // Format: "   0 C   :    0.123456   -0.234567    0.345678"
        // or just: "   0.123456   -0.234567    0.345678"
        std::istringstream iss(line);
        // Try to skip index and symbol
        int idx;
        std::string symbol;
        if (iss >> idx >> symbol) {
            // Skip colon if present
            if (symbol == ":") {
                // already skipped, read next token as first gradient component
                if (iss >> g(i, 0) >> g(i, 1) >> g(i, 2)) continue;
            } else {
                // symbol was element name, next should be colon
                char c;
                if (iss >> c && c == ':') {
                    if (iss >> g(i, 0) >> g(i, 1) >> g(i, 2)) continue;
                } else {
                    // No colon, try reading as numbers directly
                    iss.clear();
                    iss.str(line);
                    if (iss >> g(i, 0) >> g(i, 1) >> g(i, 2)) continue;
                }
            }
        }
        // If we can't parse, set to zero
        g.row(i).setZero();
    }

    return g;
}

Vector OrcaInterface::parseCharges(int natoms) const
{
    // Try JSON first
    if (!m_property_json.is_null() && m_property_json.contains("Calculation")) {
        try {
            auto calc = m_property_json["Calculation"];
            if (calc.contains("CalculationType") &&
                calc["CalculationType"].contains("Properties")) {
                auto props = calc["CalculationType"]["Properties"];
                if (props.contains("MullikenCharges")) {
                    auto ch = props["MullikenCharges"];
                    if (ch.is_array() && ch.size() == static_cast<size_t>(natoms)) {
                        Vector charges(natoms);
                        for (int i = 0; i < natoms; ++i)
                            charges(i) = ch[i].get<double>();
                        return charges;
                    }
                }
                if (props.contains("LoewdinCharges")) {
                    auto ch = props["LoewdinCharges"];
                    if (ch.is_array() && ch.size() == static_cast<size_t>(natoms)) {
                        Vector charges(natoms);
                        for (int i = 0; i < natoms; ++i)
                            charges(i) = ch[i].get<double>();
                        return charges;
                    }
                }
            }
        } catch (...) {
            // Fall through to text parsing
        }
    }

    return parseChargesFromText(natoms);
}

Vector OrcaInterface::parseChargesFromText(int natoms) const
{
    std::ifstream out(m_output_path);
    if (!out) return Vector::Zero(natoms);

    std::string line;
    bool found = false;

    while (std::getline(out, line)) {
        if (line.find("MULLIKEN ATOMIC CHARGES") != std::string::npos) {
            found = true;
            break;
        }
    }

    if (!found) {
        return Vector::Zero(natoms);
    }

    // Skip header lines (usually 1-2)
    for (int skip = 0; skip < 2; ++skip) {
        if (!std::getline(out, line)) break;
    }

    Vector charges(natoms);
    for (int i = 0; i < natoms; ++i) {
        if (!std::getline(out, line)) break;
        std::istringstream iss(line);
        int idx;
        std::string symbol;
        double q;
        if (iss >> idx >> symbol >> q) {
            charges(i) = q;
        } else {
            charges(i) = 0.0;
        }
    }

    return charges;
}

Position OrcaInterface::parseDipole() const
{
    if (!m_property_json.is_null() && m_property_json.contains("Calculation")) {
        try {
            auto calc = m_property_json["Calculation"];
            if (calc.contains("CalculationType") &&
                calc["CalculationType"].contains("Properties") &&
                calc["CalculationType"]["Properties"].contains("DipoleMoment")) {
                auto dip = calc["CalculationType"]["Properties"]["DipoleMoment"];
                if (dip.is_array() && dip.size() >= 3) {
                    return Position{dip[0].get<double>(),
                                    dip[1].get<double>(),
                                    dip[2].get<double>()};
                }
            }
        } catch (...) {
        }
    }
    return {0.0, 0.0, 0.0};
}

void OrcaInterface::citeOrca()
{
    CitationRegistry::cite("orca");
}

double OrcaInterface::calculate(bool gradient, int verbosity)
{
    clearError();
    m_property_json = json{};
    m_last_energy = 0.0;
    m_last_gradient = Matrix::Zero(0, 3);
    m_last_charges = Vector::Zero(0);
    m_last_dipole = {0.0, 0.0, 0.0};

    std::string orca_exe = findOrcaExecutable();
    if (orca_exe.empty()) {
        m_has_error = true;
        m_error_message = "ORCA executable not found. Set ORCA_PATH environment variable or use -orca_executable /full/path/to/orca.";
        return 0.0;
    }

    citeOrca();

    // Ensure gradient keyword is present if needed
    if (gradient) {
        // Read input, check for EnGrad
        std::ifstream in(m_input_path);
        if (in) {
            std::string content((std::istreambuf_iterator<char>(in)),
                                std::istreambuf_iterator<char>());
            in.close();
            if (content.find("EnGrad") == std::string::npos &&
                content.find("engrad") == std::string::npos) {
                // Append EnGrad to the ! line
                size_t bang = content.find("!");
                if (bang != std::string::npos) {
                    size_t eol = content.find("\n", bang);
                    content.insert(eol, " EnGrad");
                    std::ofstream out(m_input_path);
                    out << content;
                    out.close();
                }
            }
        }
    }

    // Run ORCA using the resolved executable path.
    // For parallel runs ORCA requires the absolute path to the binary.
    // Live output filter:
    //   verbosity >= 3: full ORCA output (debug)
    //   verbosity == 2: SCF iterations and final energy only
    std::stringstream cmd;
    cmd << "\"" << orca_exe << "\" \"" << m_input_path << "\" 2>&1";

#ifdef _WIN32
    FILE* pipe = _popen(cmd.str().c_str(), "r");
#else
    FILE* pipe = popen(cmd.str().c_str(), "r");
#endif
    if (!pipe) {
        m_has_error = true;
        m_error_message = "Failed to start ORCA process. Command: " + cmd.str();
        return 0.0;
    }

    std::ofstream out_file(m_output_path);
    char buffer[4096];
    bool in_scf_block = false;
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        std::string line(buffer);
        // Verbosity 3: full ORCA output.
        // Verbosity 2: SCF iterations and final results only.
        bool show = (verbosity >= 3);
        if (verbosity == 2) {
            // SCF block starts with "SCF Procedure:" or "ORCA LEAN-SCF"
            if (line.find("SCF Procedure:") != std::string::npos ||
                line.find("ORCA LEAN-SCF") != std::string::npos ||
                line.find("memory conserving SCF solver") != std::string::npos) {
                in_scf_block = true;
            }
            // SCF block ends after convergence/success messages
            if (line.find("SCF CONVERGED AFTER") != std::string::npos ||
                line.find("Gradient check signals convergence") != std::string::npos ||
                line.find("SUCCESS") != std::string::npos) {
                in_scf_block = false;
                show = true;  // show the convergence line itself
            } else if (in_scf_block) {
                show = true;
            }
            // Always show final results regardless of block state
            if (line.find("FINAL SINGLE POINT") != std::string::npos ||
                line.find("ORCA TERMINATED") != std::string::npos) {
                show = true;
            }
        }
        if (show) {
            std::cout << line;
        }
        if (out_file) {
            out_file << line;
        }
    }
    out_file.close();

#ifdef _WIN32
    int result = _pclose(pipe);
#else
    int status = pclose(pipe);
    int result = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
#endif
    if (result != 0) {
        m_has_error = true;
        m_error_message = "ORCA process exited with non-zero status (" + std::to_string(result) + ").";
        return 0.0;
    }

    // Try to get JSON output
    produceAndLoadJSON();

    // Parse results
    m_last_energy = parseEnergy();

    // Determine natoms from input (rough)
    // We'll read the input to count atoms in the xyz block
    int natoms = 0;
    {
        std::ifstream inp(m_input_path);
        if (inp) {
            std::string line;
            bool in_xyz = false;
            while (std::getline(inp, line)) {
                if (line.find("*xyz") != std::string::npos || line.find("* XYZ") != std::string::npos) {
                    in_xyz = true;
                    continue;
                }
                if (in_xyz) {
                    if (line.find("*") != std::string::npos && line.find("*xyz") == std::string::npos) {
                        break;
                    }
                    if (!line.empty() && line.find_first_not_of(" \t\r\n") != std::string::npos) {
                        natoms++;
                    }
                }
            }
        }
    }

    if (gradient && natoms > 0) {
        m_last_gradient = parseGradient(natoms);
    }
    if (natoms > 0) {
        m_last_charges = parseCharges(natoms);
    }
    m_last_dipole = parseDipole();

    return m_last_energy;
}
