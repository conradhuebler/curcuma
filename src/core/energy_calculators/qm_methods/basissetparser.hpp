#pragma once

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "dft_integrals.hpp"  // primitiveNorm, normalizeOrbitalSelfOverlap (WP1)
#include "GTOIntegrals.hpp"
#include "STOIntegrals.hpp"

// Forward declarations
namespace STO {
struct Orbital;
}
namespace GTO {
struct Orbital;
}

// Basis Set Parser namespace
namespace BasisSetParser {

// Shell type enumeration
enum ShellType {
    S_SHELL = 0,
    P_SHELL = 1,
    D_SHELL = 2,
    F_SHELL = 3,
    G_SHELL = 4
};

// Shell information structure
struct BasisShell {
    ShellType type; // Type of shell (S, P, D, etc.)
    std::vector<double> exponents; // Exponents
    std::vector<std::vector<double>> coefficients; // Coefficients for each contraction
    int numContractions; // Number of contractions in this shell
};

// Basis set for an element
struct ElementBasis {
    std::string symbol; // Element symbol
    int atomicNumber; // Atomic number
    std::vector<BasisShell> shells; // List of shells
    double vsip_s; // VSIP for s orbitals
    double vsip_p; // VSIP for p orbitals
    double vsip_d; // VSIP for d orbitals
};

// Map from element name/symbol to basis set
using BasisSetMap = std::map<std::string, ElementBasis>;

// Convert a shell type character to ShellType enum
inline ShellType charToShellType(char c)
{
    switch (std::toupper(c)) {
    case 'S':
        return S_SHELL;
    case 'P':
        return P_SHELL;
    case 'D':
        return D_SHELL;
    case 'F':
        return F_SHELL;
    case 'G':
        return G_SHELL;
    default:
        throw std::runtime_error("Unsupported shell type: " + std::string(1, c));
    }
}

// Convert ShellType enum to string
inline std::string shellTypeToString(ShellType type)
{
    switch (type) {
    case S_SHELL:
        return "S";
    case P_SHELL:
        return "P";
    case D_SHELL:
        return "D";
    case F_SHELL:
        return "F";
    case G_SHELL:
        return "G";
    default:
        return "Unknown";
    }
}

// Trim whitespace from a string
inline std::string trim(const std::string& str)
{
    size_t first = str.find_first_not_of(" \t\n\r\f\v");
    if (first == std::string::npos)
        return "";
    size_t last = str.find_last_not_of(" \t\n\r\f\v");
    return str.substr(first, (last - first + 1));
}

// Convert string to uppercase
inline std::string toUpper(std::string str)
{
    std::transform(str.begin(), str.end(), str.begin(),
        [](unsigned char c) { return std::toupper(c); });
    return str;
}

// Parse a TURBOMOLE/ORCA basis set file in the curcuma $DATA format.
// Claude Generated (WP1): rewritten as a single-pass line parser. The previous
// implementation misclassified indented shell lines (e.g. "S   3") as element
// headers and used a seekg-by-trimmed-length hack that broke on leading
// whitespace, so it never parsed the shipped def2-SV(P)/def2-SVP files. This
// version detects element headers (a single all-alpha token with no digits --
// full names like HYDROGEN or standard symbols like H/He) and keys the result
// map by the standard element symbol, so DFT can look up by symbol derived
// from the atomic number.
inline BasisSetMap parseBasisSetFile(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open basis set file: " + filename);
    }

    BasisSetMap basisSetMap;
    std::string line;

    // Skip header until the $DATA marker is found. Match the marker itself (a
    // line whose trimmed content starts with "$DATA"), not just any line that
    // happens to contain the substring in a comment.
    bool foundData = false;
    while (std::getline(file, line)) {
        std::string t = trim(line);
        if (t == "$DATA" || t.rfind("$DATA", 0) == 0) {
            foundData = true;
            break;
        }
    }
    if (!foundData) {
        throw std::runtime_error("Invalid basis set file format, $DATA not found.");
    }

    // Full element name / symbol (any case) -> standard symbol. Scope: H-Ne,
    // the elements present in the shipped def2-SVP file. Keys are uppercase;
    // looked up via toUpper(first).
    static const std::map<std::string, std::string> nameToSymbol = {
        { "HYDROGEN", "H" }, { "HELIUM", "He" }, { "LITHIUM", "Li" },
        { "BERYLLIUM", "Be" }, { "BORON", "B" }, { "CARBON", "C" },
        { "NITROGEN", "N" }, { "OXYGEN", "O" }, { "FLUORINE", "F" },
        { "NEON", "Ne" },
        { "H", "H" }, { "HE", "He" }, { "LI", "Li" }, { "BE", "Be" },
        { "B", "B" }, { "C", "C" }, { "N", "N" }, { "O", "O" },
        { "F", "F" }, { "NE", "Ne" }
    };
    // VSIP defaults per symbol (used only by the Hückel STO path; the DFT 1e
    // integrals ignore VSIP).
    static const std::map<std::string, std::tuple<int, double, double, double>> elemInfo = {
        { "H", { 1, -13.6, 0.0, 0.0 } }, { "He", { 2, -24.6, 0.0, 0.0 } },
        { "Li", { 3, -5.39, -3.67, 0.0 } }, { "Be", { 4, -9.32, -4.0, 0.0 } },
        { "B", { 5, -15.2, -8.3, 0.0 } }, { "C", { 6, -19.44, -10.67, 0.0 } },
        { "N", { 7, -26.0, -13.4, 0.0 } }, { "O", { 8, -32.3, -14.8, 0.0 } },
        { "F", { 9, -40.0, -18.0, 0.0 } }, { "Ne", { 10, -48.0, -23.0, 0.0 } }
    };

    auto isAllAlpha = [](const std::string& s) {
        if (s.empty()) return false;
        for (char c : s)
            if (!std::isalpha((unsigned char)c)) return false;
        return true;
    };

    ElementBasis current;
    bool haveCurrent = false;
    auto flush = [&]() {
        if (!haveCurrent) return;
        basisSetMap[current.symbol] = current;
        haveCurrent = false;
    };

    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty()) continue;
        if (line.find("$END") != std::string::npos) { flush(); break; }
        if (line[0] == '!') continue;  // comment

        std::istringstream iss(line);
        std::string first, second;
        iss >> first >> second;

        // Element header: a single all-alpha token (no digits, no count). A bare
        // shell letter (S/P/D/F/G) with no count is NOT a header in this format.
        bool isHeader = false;
        if (isAllAlpha(first)) {
            std::string up = toUpper(first);
            bool shellLetter = (up == "S" || up == "P" || up == "D" || up == "F" || up == "G");
            if (second.empty()) {
                isHeader = !shellLetter;
            } else {
                isHeader = !std::isdigit((unsigned char)second[0]);
            }
            if (isHeader && nameToSymbol.find(up) == nameToSymbol.end())
                isHeader = false;
        }

        if (isHeader) {
            flush();
            std::string sym = nameToSymbol.at(toUpper(first));
            current = ElementBasis();
            current.symbol = sym;
            auto ei = elemInfo.find(sym);
            if (ei != elemInfo.end()) {
                current.atomicNumber = std::get<0>(ei->second);
                current.vsip_s = std::get<1>(ei->second);
                current.vsip_p = std::get<2>(ei->second);
                current.vsip_d = std::get<3>(ei->second);
            } else {
                current.atomicNumber = 0;
                current.vsip_s = -10.0;
                current.vsip_p = -5.0;
                current.vsip_d = 0.0;
            }
            haveCurrent = true;
            continue;
        }

        // Shell line: "<S|P|D|F|G> <numPrim> [<numContractions>]"
        if (!haveCurrent)
            throw std::runtime_error("Shell line before any element header: " + line);
        if (first.empty())
            throw std::runtime_error("Empty shell type string");
        ShellType shellType = charToShellType(first[0]);
        int numPrimitives = 0;
        if (!(std::istringstream(second) >> numPrimitives))
            throw std::runtime_error("Failed to parse shell primitive count: " + line);

        BasisShell shell;
        shell.type = shellType;
        shell.numContractions = 1;
        {
            std::istringstream iss2(line);
            std::string tmp;
            int np2;
            iss2 >> tmp >> np2;
            if (iss2 >> shell.numContractions) { /* general contraction */ }
        }
        shell.coefficients.resize(shell.numContractions);

        for (int i = 0; i < numPrimitives; ++i) {
            if (!std::getline(file, line))
                throw std::runtime_error("Unexpected end of file while parsing basis set.");
            line = trim(line);
            std::istringstream pss(line);
            int idx;
            double exponent;
            if (!(pss >> idx >> exponent))
                throw std::runtime_error("Failed to parse primitive line: " + line);
            shell.exponents.push_back(exponent);
            for (int j = 0; j < shell.numContractions; ++j) {
                double coefficient;
                if (!(pss >> coefficient))
                    throw std::runtime_error("Missing coefficient for contraction "
                        + std::to_string(j + 1) + " in primitive " + std::to_string(i + 1));
                shell.coefficients[j].push_back(coefficient);
            }
        }
        current.shells.push_back(shell);
    }
    flush();
    return basisSetMap;
}

// Create a minimal STO-type basis from a GTO basis set
inline std::vector<STO::Orbital> createSTOFromGTOBasis(
    const ElementBasis& basis,
    double x, double y, double z,
    int atomIndex,
    bool includeAllOrbitals = false)
{
    std::vector<STO::Orbital> stoOrbitals;

    for (const BasisShell& shell : basis.shells) {
        // For STO approximation, we'll use only the first contraction of each shell
        // and extract an effective Slater exponent

        if (shell.exponents.empty() || shell.coefficients.empty() || shell.coefficients[0].empty())
            continue;

        // Calculate effective zeta for STO based on the first primitive GTO
        // This is a very rough approximation
        double effectiveZeta = std::sqrt(shell.exponents[0] / 2.0);

        // Based on shell type, create appropriate STO orbitals
        switch (shell.type) {
        case S_SHELL: {
            STO::Orbital orbital;
            orbital.x = x;
            orbital.y = y;
            orbital.z = z;
            orbital.type = STO::OrbitalType::S;
            orbital.zeta = effectiveZeta;
            orbital.VSIP = basis.vsip_s;
            orbital.atom = atomIndex;
            stoOrbitals.push_back(orbital);
        } break;
        case P_SHELL: {
            STO::Orbital orbital_px;
            orbital_px.x = x;
            orbital_px.y = y;
            orbital_px.z = z;
            orbital_px.type = STO::OrbitalType::PX;
            orbital_px.zeta = effectiveZeta;
            orbital_px.VSIP = basis.vsip_p;
            orbital_px.atom = atomIndex;
            stoOrbitals.push_back(orbital_px);

            STO::Orbital orbital_py;
            orbital_py.x = x;
            orbital_py.y = y;
            orbital_py.z = z;
            orbital_py.type = STO::OrbitalType::PY;
            orbital_py.zeta = effectiveZeta;
            orbital_py.VSIP = basis.vsip_p;
            orbital_py.atom = atomIndex;
            stoOrbitals.push_back(orbital_py);

            STO::Orbital orbital_pz;
            orbital_pz.x = x;
            orbital_pz.y = y;
            orbital_pz.z = z;
            orbital_pz.type = STO::OrbitalType::PZ;
            orbital_pz.zeta = effectiveZeta;
            orbital_pz.VSIP = basis.vsip_p;
            orbital_pz.atom = atomIndex;
            stoOrbitals.push_back(orbital_pz);
        } break;
        case D_SHELL:
            if (includeAllOrbitals) {
                // Add d orbitals only if explicitly requested
                STO::Orbital orbital_dxy;
                orbital_dxy.x = x;
                orbital_dxy.y = y;
                orbital_dxy.z = z;
                orbital_dxy.type = STO::OrbitalType::DXY;
                orbital_dxy.zeta = effectiveZeta;
                orbital_dxy.VSIP = basis.vsip_d;
                orbital_dxy.atom = atomIndex;
                stoOrbitals.push_back(orbital_dxy);

                STO::Orbital orbital_dyz;
                orbital_dyz.x = x;
                orbital_dyz.y = y;
                orbital_dyz.z = z;
                orbital_dyz.type = STO::OrbitalType::DYZ;
                orbital_dyz.zeta = effectiveZeta;
                orbital_dyz.VSIP = basis.vsip_d;
                orbital_dyz.atom = atomIndex;
                stoOrbitals.push_back(orbital_dyz);

                STO::Orbital orbital_dzx;
                orbital_dzx.x = x;
                orbital_dzx.y = y;
                orbital_dzx.z = z;
                orbital_dzx.type = STO::OrbitalType::DZX;
                orbital_dzx.zeta = effectiveZeta;
                orbital_dzx.VSIP = basis.vsip_d;
                orbital_dzx.atom = atomIndex;
                stoOrbitals.push_back(orbital_dzx);

                STO::Orbital orbital_dx2y2;
                orbital_dx2y2.x = x;
                orbital_dx2y2.y = y;
                orbital_dx2y2.z = z;
                orbital_dx2y2.type = STO::OrbitalType::DX2Y2;
                orbital_dx2y2.zeta = effectiveZeta;
                orbital_dx2y2.VSIP = basis.vsip_d;
                orbital_dx2y2.atom = atomIndex;
                stoOrbitals.push_back(orbital_dx2y2);

                STO::Orbital orbital_dz2;
                orbital_dz2.x = x;
                orbital_dz2.y = y;
                orbital_dz2.z = z;
                orbital_dz2.type = STO::OrbitalType::DZ2;
                orbital_dz2.zeta = effectiveZeta;
                orbital_dz2.VSIP = basis.vsip_d;
                orbital_dz2.atom = atomIndex;
                stoOrbitals.push_back(orbital_dz2);
            }
            break;
        // Higher angular momentum orbitals not commonly used in Extended Hückel
        default:
            break;
        }

        // For STO approximation, we typically only use the first shell of each type
        // unless includeAllOrbitals is true
        if (!includeAllOrbitals) {
            if (shell.type == S_SHELL || shell.type == P_SHELL || shell.type == D_SHELL) {
                // Skip other shells of the same type
                ShellType currentType = shell.type;
                while (shell.type == currentType) {
                    // Move to the next shell
                    if (&shell == &basis.shells.back())
                        break;
                    // shell = *(&shell + 1);
                }
            }
        }
    }

    return stoOrbitals;
}

// Create GTO orbitals from the parsed basis set.
// Claude Generated (WP1): iterates ALL contractions of each shell (forward-
// correct for generally-contracted bases such as cc-pVnZ; def2-SVP has
// numContractions==1 so numerically unchanged), emits the proper 6 cartesian d
// components in the order [DXX, DYY, DZZ, DXY, DXZ, DYZ] expected by the
// spherical transform, and pre-normalizes primitives (Turbomole basis files
// give coefficients for already-normalized primitives, so we multiply by the
// primitive norm to recover coefficients for the raw x^l exp(-a r^2) kernels)
// and then renormalizes each contracted AO so S_ii = 1. f/g shells throw.
inline std::vector<GTO::Orbital> createGTOFromBasis(
    const ElementBasis& basis,
    double x, double y, double z,
    int atomIndex)
{
    std::vector<GTO::Orbital> gtoOrbitals;

    // Build one pre-normalized contracted AO for a given cartesian component.
    auto makeOrbital = [&](GTO::OrbitalType type, int l, int m, int n,
                           const std::vector<double>& exps,
                           const std::vector<double>& rawCoeffs,
                           double vsip) {
        GTO::Orbital orbital;
        orbital.x = x;
        orbital.y = y;
        orbital.z = z;
        orbital.type = type;
        orbital.exponents = exps;
        orbital.coefficients.resize(exps.size());
        // Turbomole convention: rawCoeffs are for normalized primitives; multiply
        // by primitiveNorm to express the contraction over unnormalized primitives.
        for (size_t a = 0; a < exps.size(); ++a)
            orbital.coefficients[a] = rawCoeffs[a] * dft1e::primitiveNorm(exps[a], l, m, n);
        orbital.VSIP = vsip;
        orbital.atom = atomIndex;
        // Renormalize the whole contracted function so S_ii = 1 (matches ORCA).
        dft1e::normalizeOrbitalSelfOverlap(orbital);
        return orbital;
    };

    for (const BasisShell& shell : basis.shells) {
        switch (shell.type) {
        case S_SHELL: {
            for (int c = 0; c < shell.numContractions; ++c)
                gtoOrbitals.push_back(
                    makeOrbital(GTO::OrbitalType::S, 0, 0, 0,
                                shell.exponents, shell.coefficients[c], basis.vsip_s));
        } break;
        case P_SHELL: {
            for (int c = 0; c < shell.numContractions; ++c) {
                gtoOrbitals.push_back(makeOrbital(GTO::OrbitalType::PX, 1, 0, 0,
                                                  shell.exponents, shell.coefficients[c], basis.vsip_p));
                gtoOrbitals.push_back(makeOrbital(GTO::OrbitalType::PY, 0, 1, 0,
                                                  shell.exponents, shell.coefficients[c], basis.vsip_p));
                gtoOrbitals.push_back(makeOrbital(GTO::OrbitalType::PZ, 0, 0, 1,
                                                  shell.exponents, shell.coefficients[c], basis.vsip_p));
            }
        } break;
        case D_SHELL: {
            // Order MUST match dft_integrals.cpp buildSphericalTransform:
            //   [DXX, DYY, DZZ, DXY, DXZ, DYZ]
            for (int c = 0; c < shell.numContractions; ++c) {
                gtoOrbitals.push_back(makeOrbital(GTO::OrbitalType::DXX, 2, 0, 0,
                                                  shell.exponents, shell.coefficients[c], basis.vsip_d));
                gtoOrbitals.push_back(makeOrbital(GTO::OrbitalType::DYY, 0, 2, 0,
                                                  shell.exponents, shell.coefficients[c], basis.vsip_d));
                gtoOrbitals.push_back(makeOrbital(GTO::OrbitalType::DZZ, 0, 0, 2,
                                                  shell.exponents, shell.coefficients[c], basis.vsip_d));
                gtoOrbitals.push_back(makeOrbital(GTO::OrbitalType::DXY, 1, 1, 0,
                                                  shell.exponents, shell.coefficients[c], basis.vsip_d));
                gtoOrbitals.push_back(makeOrbital(GTO::OrbitalType::DXZ, 1, 0, 1,
                                                  shell.exponents, shell.coefficients[c], basis.vsip_d));
                gtoOrbitals.push_back(makeOrbital(GTO::OrbitalType::DYZ, 0, 1, 1,
                                                  shell.exponents, shell.coefficients[c], basis.vsip_d));
            }
        } break;
        case F_SHELL:
        case G_SHELL:
            throw std::runtime_error(
                "createGTOFromBasis: f/g shells not supported (def2-SVP H-Ne has none).");
        default:
            break;
        }
    }

    return gtoOrbitals;
}

// Set VSIP values for a basis set
inline void setVSIPValues(BasisSetMap& basisSet, const std::string& element,
    double vsip_s, double vsip_p, double vsip_d = 0.0)
{
    auto it = basisSet.find(element);
    if (it != basisSet.end()) {
        it->second.vsip_s = vsip_s;
        it->second.vsip_p = vsip_p;
        it->second.vsip_d = vsip_d;
    } else {
        std::cerr << "Warning: Element " << element << " not found in basis set." << std::endl;
    }
}

// Set VSIP values from a configuration file
inline void loadVSIPValues(BasisSetMap& basisSet, const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open VSIP configuration file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#')
            continue; // Skip empty lines and comments

        std::istringstream iss(line);
        std::string element;
        double vsip_s, vsip_p, vsip_d = 0.0;

        if (iss >> element >> vsip_s >> vsip_p) {
            // Optional d-VSIP
            iss >> vsip_d;
            setVSIPValues(basisSet, element, vsip_s, vsip_p, vsip_d);
        } else {
            std::cerr << "Warning: Invalid line in VSIP file: " << line << std::endl;
        }
    }
}

// Print basis set information
inline void printBasisSetInfo(const BasisSetMap& basisSet)
{
    std::cout << "Basis Set Information:" << std::endl;
    std::cout << "======================" << std::endl;

    for (const auto& [element, basis] : basisSet) {
        std::cout << "Element: " << element << " (Z=" << basis.atomicNumber << ")" << std::endl;
        std::cout << "  VSIP values: s=" << basis.vsip_s << ", p=" << basis.vsip_p
                  << ", d=" << basis.vsip_d << " eV" << std::endl;

        std::cout << "  Shells:" << std::endl;
        for (const auto& shell : basis.shells) {
            std::cout << "    " << shellTypeToString(shell.type)
                      << " shell with " << shell.exponents.size() << " primitive GTOs and "
                      << shell.numContractions << " contraction(s)" << std::endl;
        }
        std::cout << std::endl;
    }
}

// Example usage
inline void exampleUsage(const std::string& filename)
{
    std::cout << "Loading basis set from: " << filename << std::endl;

    try {
        // Parse the basis set file
        BasisSetMap basisSet = parseBasisSetFile(filename);

        // Set custom VSIP values for Extended Hückel
        setVSIPValues(basisSet, "H", -13.6, 0.0);
        setVSIPValues(basisSet, "C", -19.44, -10.67);

        // Print basis set info
        printBasisSetInfo(basisSet);

        // Create STO orbitals for carbon at origin
        if (basisSet.find("C") != basisSet.end()) {
            std::vector<STO::Orbital> carbonOrbitals = createSTOFromGTOBasis(basisSet.at("C"), 0.0, 0.0, 0.0, 0);

            std::cout << "Created " << carbonOrbitals.size()
                      << " STO orbitals for carbon." << std::endl;
        }

        // Create GTO orbitals for hydrogen at (1,1,1)
        if (basisSet.find("H") != basisSet.end()) {
            std::vector<GTO::Orbital> hydrogenOrbitals = createGTOFromBasis(basisSet.at("H"), 1.0, 1.0, 1.0, 1);

            std::cout << "Created " << hydrogenOrbitals.size()
                      << " GTO orbitals for hydrogen." << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

} // namespace BasisSetParser
