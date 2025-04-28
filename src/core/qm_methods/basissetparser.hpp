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
#include <vector>

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
ShellType charToShellType(char c)
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
std::string shellTypeToString(ShellType type)
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
std::string trim(const std::string& str)
{
    size_t first = str.find_first_not_of(" \t\n\r\f\v");
    if (first == std::string::npos)
        return "";
    size_t last = str.find_last_not_of(" \t\n\r\f\v");
    return str.substr(first, (last - first + 1));
}

// Convert string to uppercase
std::string toUpper(std::string str)
{
    std::transform(str.begin(), str.end(), str.begin(),
        [](unsigned char c) { return std::toupper(c); });
    return str;
}

// Parse a TURBOMOLE/ORCA basis set file
BasisSetMap parseBasisSetFile(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open basis set file: " + filename);
    }

    BasisSetMap basisSetMap;
    std::string line;

    // Skip header until $DATA is found
    bool foundData = false;
    while (std::getline(file, line)) {
        if (line.find("$DATA") != std::string::npos) {
            foundData = true;
            break;
        }
    }

    if (!foundData) {
        throw std::runtime_error("Invalid basis set file format, $DATA not found.");
    }

    // Parse basis sets for each element
    while (std::getline(file, line)) {
        std::cout << line << std::endl;

        line = trim(line);
        if (line.empty())
            continue;

        // Check for end of file
        if (line.find("$END") != std::string::npos) {
            break;
        }

        // Read element name
        std::string elementName = line;
        ElementBasis elementBasis;
        elementBasis.symbol = elementName;

        // Set default VSIP values (can be overridden later)
        // These are approximate values based on common Extended Hückel parameters
        if (elementName == "HYDROGEN" || elementName == "H") {
            elementBasis.atomicNumber = 1;
            elementBasis.vsip_s = -13.6; // eV
            elementBasis.vsip_p = 0.0; // Not applicable
            elementBasis.vsip_d = 0.0; // Not applicable
        } else if (elementName == "CARBON" || elementName == "C") {
            elementBasis.atomicNumber = 6;
            elementBasis.vsip_s = -19.44; // eV
            elementBasis.vsip_p = -10.67; // eV
            elementBasis.vsip_d = 0.0; // Not commonly used
        } else if (elementName == "NITROGEN" || elementName == "N") {
            elementBasis.atomicNumber = 7;
            elementBasis.vsip_s = -26.0; // eV
            elementBasis.vsip_p = -13.4; // eV
            elementBasis.vsip_d = 0.0; // Not commonly used
        } else if (elementName == "OXYGEN" || elementName == "O") {
            elementBasis.atomicNumber = 8;
            elementBasis.vsip_s = -32.3; // eV
            elementBasis.vsip_p = -14.8; // eV
            elementBasis.vsip_d = 0.0; // Not commonly used
        } else if (elementName == "SODIUM" || elementName == "Na") {
            elementBasis.atomicNumber = 11;
            elementBasis.vsip_s = -5.14; // eV
            elementBasis.vsip_p = -3.04; // eV
            elementBasis.vsip_d = 0.0; // Not commonly used
        } else {
            // For other elements, use approximate values
            elementBasis.atomicNumber = 0; // Unknown
            elementBasis.vsip_s = -10.0; // Placeholder
            elementBasis.vsip_p = -5.0; // Placeholder
            elementBasis.vsip_d = 0.0; // Placeholder
        }

        // Parse shells
        while (std::getline(file, line)) {
            line = trim(line);
            if (line.empty())
                continue;

            // Check if we've reached the next element
            if (line.find("$END") != std::string::npos || (std::isalpha(line[0]) && line.length() > 1 && !std::isdigit(line[1]))) {
                file.seekg(-static_cast<int>(line.length()) - 1, std::ios_base::cur);
                break;
            }

            // Parse shell type and number of primitive functions
            std::istringstream iss(line);
            std::string shellTypeStr;
            int numPrimitives;

            // Lese Shell-Typ und Anzahl der primitiven Funktionen
            if (!(iss >> shellTypeStr >> numPrimitives)) {
                throw std::runtime_error("Failed to parse shell type and number: " + line);
            }

            if (shellTypeStr.empty()) {
                throw std::runtime_error("Empty shell type string");
            }

            // Nur den ersten Buchstaben des Shell-Typs verwenden
            ShellType shellType = charToShellType(shellTypeStr[0]);

            BasisShell shell;
            shell.type = shellType;

            // Second number, if present, is the number of contractions
            shell.numContractions = 1; // Default: one contraction
            if (iss >> shell.numContractions) {
                // If present, it indicates the number of contractions
            }

            // Resize the coefficients vector based on number of contractions
            shell.coefficients.resize(shell.numContractions);

            // Read exponents and coefficients
            for (int i = 0; i < numPrimitives; ++i) {
                if (!std::getline(file, line)) {
                    throw std::runtime_error("Unexpected end of file while parsing basis set.");
                }

                line = trim(line);
                iss.clear();
                iss.str(line);

                int idx;
                double exponent;
                iss >> idx >> exponent;

                shell.exponents.push_back(exponent);

                // Read coefficients for each contraction
                for (int j = 0; j < shell.numContractions; ++j) {
                    double coefficient;
                    if (!(iss >> coefficient)) {
                        throw std::runtime_error("Missing coefficient for contraction " + std::to_string(j + 1) + " in primitive " + std::to_string(i + 1));
                    }
                    shell.coefficients[j].push_back(coefficient);
                }
            }

            elementBasis.shells.push_back(shell);
        }

        // Add the element basis to the map
        basisSetMap[elementName] = elementBasis;
        // Also add the element symbol as a key if it's different
        if (elementName.length() > 2) {
            // Extract standard element symbol (first or first two letters)
            std::string symbol;
            if (elementName.length() > 1 && std::islower(elementName[1])) {
                symbol = elementName.substr(0, 2);
            } else {
                symbol = elementName.substr(0, 1);
            }
            symbol[0] = std::toupper(symbol[0]);
            if (symbol.length() > 1) {
                symbol[1] = std::tolower(symbol[1]);
            }
            basisSetMap[symbol] = elementBasis;
        }
    }

    return basisSetMap;
}

// Create a minimal STO-type basis from a GTO basis set
std::vector<STO::Orbital> createSTOFromGTOBasis(
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

// Create GTO orbitals from the parsed basis set
std::vector<GTO::Orbital> createGTOFromBasis(
    const ElementBasis& basis,
    double x, double y, double z,
    int atomIndex)
{
    std::vector<GTO::Orbital> gtoOrbitals;

    for (const BasisShell& shell : basis.shells) {
        switch (shell.type) {
        case S_SHELL: {
            GTO::Orbital orbital;
            orbital.x = x;
            orbital.y = y;
            orbital.z = z;
            orbital.type = GTO::OrbitalType::S;
            orbital.exponents = shell.exponents;
            // Use the first contraction by default
            orbital.coefficients = shell.coefficients[0];
            orbital.VSIP = basis.vsip_s;
            orbital.atom = atomIndex;
            gtoOrbitals.push_back(orbital);
        } break;
        case P_SHELL: {
            // For each contraction, create three p-orbitals (px, py, pz)
            GTO::Orbital orbital_px;
            orbital_px.x = x;
            orbital_px.y = y;
            orbital_px.z = z;
            orbital_px.type = GTO::OrbitalType::PX;
            orbital_px.exponents = shell.exponents;
            orbital_px.coefficients = shell.coefficients[0];
            orbital_px.VSIP = basis.vsip_p;
            orbital_px.atom = atomIndex;
            gtoOrbitals.push_back(orbital_px);

            GTO::Orbital orbital_py;
            orbital_py.x = x;
            orbital_py.y = y;
            orbital_py.z = z;
            orbital_py.type = GTO::OrbitalType::PY;
            orbital_py.exponents = shell.exponents;
            orbital_py.coefficients = shell.coefficients[0];
            orbital_py.VSIP = basis.vsip_p;
            orbital_py.atom = atomIndex;
            gtoOrbitals.push_back(orbital_py);

            GTO::Orbital orbital_pz;
            orbital_pz.x = x;
            orbital_pz.y = y;
            orbital_pz.z = z;
            orbital_pz.type = GTO::OrbitalType::PZ;
            orbital_pz.exponents = shell.exponents;
            orbital_pz.coefficients = shell.coefficients[0];
            orbital_pz.VSIP = basis.vsip_p;
            orbital_pz.atom = atomIndex;
            gtoOrbitals.push_back(orbital_pz);
        } break;
        case D_SHELL: {
            // Create five d-orbitals
            GTO::Orbital orbital_dxy;
            orbital_dxy.x = x;
            orbital_dxy.y = y;
            orbital_dxy.z = z;
            orbital_dxy.type = GTO::OrbitalType::DXY;
            orbital_dxy.exponents = shell.exponents;
            orbital_dxy.coefficients = shell.coefficients[0];
            orbital_dxy.VSIP = basis.vsip_d;
            orbital_dxy.atom = atomIndex;
            gtoOrbitals.push_back(orbital_dxy);

            GTO::Orbital orbital_dyz;
            orbital_dyz.x = x;
            orbital_dyz.y = y;
            orbital_dyz.z = z;
            orbital_dyz.type = GTO::OrbitalType::DYZ;
            orbital_dyz.exponents = shell.exponents;
            orbital_dyz.coefficients = shell.coefficients[0];
            orbital_dyz.VSIP = basis.vsip_d;
            orbital_dyz.atom = atomIndex;
            gtoOrbitals.push_back(orbital_dyz);

            GTO::Orbital orbital_dzx;
            orbital_dzx.x = x;
            orbital_dzx.y = y;
            orbital_dzx.z = z;
            orbital_dzx.type = GTO::OrbitalType::DZX;
            orbital_dzx.exponents = shell.exponents;
            orbital_dzx.coefficients = shell.coefficients[0];
            orbital_dzx.VSIP = basis.vsip_d;
            orbital_dzx.atom = atomIndex;
            gtoOrbitals.push_back(orbital_dzx);

            GTO::Orbital orbital_dx2y2;
            orbital_dx2y2.x = x;
            orbital_dx2y2.y = y;
            orbital_dx2y2.z = z;
            orbital_dx2y2.type = GTO::OrbitalType::DX2Y2;
            orbital_dx2y2.exponents = shell.exponents;
            orbital_dx2y2.coefficients = shell.coefficients[0];
            orbital_dx2y2.VSIP = basis.vsip_d;
            orbital_dx2y2.atom = atomIndex;
            gtoOrbitals.push_back(orbital_dx2y2);

            GTO::Orbital orbital_dz2;
            orbital_dz2.x = x;
            orbital_dz2.y = y;
            orbital_dz2.z = z;
            orbital_dz2.type = GTO::OrbitalType::DZ2;
            orbital_dz2.exponents = shell.exponents;
            orbital_dz2.coefficients = shell.coefficients[0];
            orbital_dz2.VSIP = basis.vsip_d;
            orbital_dz2.atom = atomIndex;
            gtoOrbitals.push_back(orbital_dz2);
        } break;
        // Higher angular momentum orbitals can be added as needed
        default:
            // Skip unsupported shell types
            break;
        }
    }

    return gtoOrbitals;
}

// Set VSIP values for a basis set
void setVSIPValues(BasisSetMap& basisSet, const std::string& element,
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
void loadVSIPValues(BasisSetMap& basisSet, const std::string& filename)
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
void printBasisSetInfo(const BasisSetMap& basisSet)
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
void exampleUsage(const std::string& filename)
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
