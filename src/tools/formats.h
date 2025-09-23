/*
 * <Geometry tools for chemical structures.>
 * Copyright (C) 2020 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#pragma once

#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <fmt/color.h>
#include <fmt/core.h>

#include "src/core/elements.h"
#include "src/core/molecule.h"
// #include "src/core/fileiterator.h"

#include "src/tools/general.h"

namespace Files {

inline StringList SplitString(const std::string& string, const char* delim)
{
#pragma message("Has to be removed, as it is redudant")
    StringList elements;
    std::string element;
    for (const char& c : string) {
        if (*delim != c)
            element.push_back(c);
        else {
            if (element.size()) {
                elements.push_back(element);
            }
            element.clear();
        }
    }
    elements.push_back(element);
    return elements;
}

inline std::pair<int, std::array<double, 3>> Line2Atoms(const std::string& line)
{
    std::array<double, 3> vector;
    int element = 0;
    std::vector<std::string> elements = SplitString(line, " ");

    element = (Elements::String2Element(elements[0]));

    if (elements.size() >= 4) {
        double x = stod(elements[1]);
        double y = stod(elements[2]);
        double z = stod(elements[3]);
        vector[0] = x;
        vector[1] = y;
        vector[2] = z;
    }

    return std::pair<int, std::array<double, 3>>(element, vector);
}

inline Mol XYZString2Mol(const std::string& coord)
{
    Mol molecule;

    std::vector<std::string> lines = SplitString(coord, "\n");
    int atoms = 0;
    int index = 0;
    int i = 0;

    for (const auto& line : lines) {
        if (index == 0) {
            try {
                molecule.m_number_atoms = stoi(line);
                atoms = stoi(line);
            } catch (const std::invalid_argument& arg) {
                atoms = 0;
            }
        }

        if (i == 1)
            molecule.m_commentline = line;
        if (i > 1) {
            auto pair = Line2Atoms(line);
            molecule.m_atoms.push_back(pair.first);
            molecule.m_geometry.conservativeResize(molecule.m_geometry.rows() + 1, molecule.m_geometry.cols());
            molecule.m_geometry.row(molecule.m_geometry.rows() - 1) = Eigen::Vector3d(pair.second[0], pair.second[1], pair.second[2]);
            //            push_back(pair.second);
        }
        if (i - 1 == atoms) {
            break;
        }
        ++i;

        index++;
    }

    return molecule;
}

inline Mol XYZ2Mol(const std::string& filename)
{
    Mol molecule;

    std::vector<std::string> lines;
    int atoms = 0;
    int index = 0;
    int i = 0;

    auto file = new std::ifstream(filename);
    for (std::string line; getline(*file, line);) {
        if (index == 0) {
            try {
                molecule.m_number_atoms = stoi(line);
                atoms = stoi(line);
                molecule.m_geometry = Eigen::MatrixXd::Zero(atoms, 3);
            } catch (const std::invalid_argument& arg) {
                atoms = 0;
            }
        }

        if (i == 1)
            molecule.m_commentline = line;
        if (i > 1) {
            auto pair = Line2Atoms(line);
            molecule.m_atoms.push_back(pair.first);
            // molecule.m_geometry.push_back(pair.second);
            // molecule.m_geometry.conservativeResize(molecule.m_geometry.rows() +  1, molecule.m_geometry.cols());
            molecule.m_geometry.row(molecule.m_atoms.size() - 1) = Eigen::Vector3d(pair.second[0], pair.second[1], pair.second[2]);
        }
        if (i - 1 == atoms) {
            break;
        }
        ++i;

        index++;
    }

    return molecule;
}

inline Mol Coord2Mol(const std::string& filename)
{
    Mol molecule;
    molecule.m_geometry = Eigen::MatrixXd::Zero(0, 3);

    auto file = new std::ifstream(filename);
    for (std::string line; getline(*file, line);) {
        //if (readblock) {
        auto strings = SplitString(line, " ");
        if (strings.size() == 4) {
            try {
                molecule.m_atoms.push_back(Elements::String2Element(strings[3]));
                std::array<double, 3> vector;
                vector[0] = std::stod(strings[0]) * au;
                vector[1] = std::stod(strings[1]) * au;
                vector[2] = std::stod(strings[2]) * au;
                // molecule.m_geometry.push_back(vector);
                molecule.m_geometry.conservativeResize(molecule.m_atoms.size(), molecule.m_geometry.cols());
                molecule.m_geometry.row(molecule.m_geometry.rows() - 1) = Eigen::Vector3d(vector[0], vector[1], vector[2]);
            } catch (const std::invalid_argument& arg) {
            }
        }
    }

    return molecule;
}
inline Mol SDF2Mol(const std::string& filename)
{
    Mol molecule;
    int readblock = false;
    auto file = new std::ifstream(filename);
    for (std::string line; getline(*file, line);) {
        //if (readblock) {
        auto strings = SplitString(line, " ");
        if (strings.size() == 16) {
            try {
                molecule.m_atoms.push_back(Elements::String2Element(strings[3]));
                std::array<double, 3> vector;
                vector[0] = std::stod(strings[0]);
                vector[1] = std::stod(strings[1]);
                vector[2] = std::stod(strings[2]);
                // molecule.m_geometry.push_back(vector);
                molecule.m_geometry.conservativeResize(molecule.m_atoms.size(), molecule.m_geometry.cols());
                molecule.m_geometry.row(molecule.m_geometry.rows() - 1) = Eigen::Vector3d(vector[0], vector[1], vector[2]);
            } catch (const std::invalid_argument& arg) {
            }
        }
    }
    return molecule;
}

inline Mol Mol22Mol(const std::string& filename)
{
    Mol molecule;
    molecule.m_geometry = Eigen::MatrixXd::Zero(0, 3);

    int read_atom = false, read_bond = false;
    auto file = new std::ifstream(filename);
    for (std::string line; getline(*file, line);) {
        if (read_atom) {
            auto strings = SplitString(line, " ");
            if (strings.size() == 9) {
                try {
                    molecule.m_atoms.push_back(Elements::String2Element(strings[1]));
                    std::array<double, 3> vector;
                    vector[0] = std::stod(strings[2]);
                    vector[1] = std::stod(strings[3]);
                    vector[2] = std::stod(strings[4]);
                    // molecule.m_geometry.push_back(vector);
                    molecule.m_geometry.conservativeResize(molecule.m_atoms.size(), molecule.m_geometry.cols());
                    molecule.m_geometry.row(molecule.m_geometry.rows() - 1) = Eigen::Vector3d(vector[0], vector[1], vector[2]);
                } catch (const std::invalid_argument& arg) {
                }
            }
        }
        if (read_bond) {
            auto strings = SplitString(line, " ");
            if (strings.size() == 4) {
                try {
                    std::pair<int, int> vector(std::stoi(strings[1]), std::stoi(strings[2]));
                    molecule.m_bonds.push_back(vector);
                } catch (const std::invalid_argument& arg) {
                }
            }
        }
        if (line.compare("@<TRIPOS>ATOM") == 0) {
            read_atom = true;
            read_bond = false;
        }
        if (line.compare("@<TRIPOS>BOND") == 0) {
            read_atom = false;
            read_bond = true;
        }
    }
    return molecule;
}

// VTF (VMD Trajectory Format) reader for Coarse Graining - Claude Generated
inline Mol VTF2Mol(const std::string& filename)
{
    Mol molecule;
    molecule.m_geometry = Eigen::MatrixXd::Zero(0, 3);

    std::map<std::string, int> atom_types; // Map CG type names to element numbers
    std::vector<std::string> atom_names;   // Store atom names for each atom
    bool reading_coords = false;
    size_t coord_idx = 0;
    bool first_timestep_read = false;

    auto file = new std::ifstream(filename);
    for (std::string line; getline(*file, line);) {
        // Remove leading/trailing whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        if (line.empty() || line[0] == '#') continue; // Skip empty lines and comments

        auto tokens = SplitString(line, " ");
        if (tokens.empty()) continue;

        if (tokens[0] == "atom") {
            // Format: atom <index> radius <value> type <type_name> name <name>
            // Alternative: atom <index> name <type_name> [type <type_index>]
            std::string atom_type;

            if (tokens.size() >= 6 && tokens[4] == "type") {
                // Format: atom <index> radius <value> type <type_name> name <name>
                atom_type = tokens[5];
            } else if (tokens.size() >= 4 && tokens[2] == "name") {
                // Format: atom <index> name <type_name> [type <type_index>]
                atom_type = tokens[3];
            } else {
                continue; // Skip malformed atom lines
            }

            // Map all CG types to CG_ELEMENT (226)
            if (atom_types.find(atom_type) == atom_types.end()) {
                atom_types[atom_type] = CG_ELEMENT;
            }

            molecule.m_atoms.push_back(CG_ELEMENT);
            atom_names.push_back(atom_type);

            // Initialize geometry for this atom (will be set from coordinate data)
            molecule.m_geometry.conservativeResize(molecule.m_atoms.size(), 3);
            molecule.m_geometry.row(molecule.m_atoms.size() - 1) = Eigen::Vector3d(0, 0, 0);
        }
        else if (tokens[0] == "bond") {
            // bond <from>:<to>
            if (tokens.size() >= 2) {
                auto bond_spec = SplitString(tokens[1], ":");
                if (bond_spec.size() == 2) {
                    try {
                        int from = std::stoi(bond_spec[0]);
                        int to = std::stoi(bond_spec[1]);
                        molecule.m_bonds.push_back(std::make_pair(from, to));
                    } catch (const std::invalid_argument& arg) {
                        // Skip invalid bond specifications
                    }
                }
            }
        }
        else if (tokens[0] == "unitcell") {
            // unitcell <a> <b> <c> [alpha] [beta] [gamma]
            // Store unit cell info in comment line for now
            if (tokens.size() >= 4) {
                molecule.m_commentline = "unitcell " + tokens[1] + " " + tokens[2] + " " + tokens[3];
                if (tokens.size() >= 7) {
                    molecule.m_commentline += " " + tokens[4] + " " + tokens[5] + " " + tokens[6];
                }
            }
        }
        else if (tokens[0] == "timestep") {
            reading_coords = true;
            // Read only the first timestep for now
            if (first_timestep_read) break;
            first_timestep_read = true;
            coord_idx = 0; // Reset coordinate index for this timestep
            continue;
        }
        else if (reading_coords && tokens.size() >= 3) {
            // Coordinate line: <x> <y> <z>
            try {
                if (coord_idx < molecule.m_atoms.size()) {
                    double x = std::stod(tokens[0]);
                    double y = std::stod(tokens[1]);
                    double z = std::stod(tokens[2]);
                    molecule.m_geometry.row(coord_idx) = Eigen::Vector3d(x, y, z);
                    coord_idx++;
                }
            } catch (const std::invalid_argument& arg) {
                // Skip invalid coordinate lines
            }
        }
    }

    molecule.m_number_atoms = molecule.m_atoms.size();
    file->close();
    delete file;

    return molecule;
}

inline Molecule LoadFile(const std::string& filename)
{
    if (std::string(filename).find(".xyz") != std::string::npos || std::string(filename).find(".trj") != std::string::npos)
        return Molecule(XYZ2Mol(filename));
    else if (std::string(filename).find(".mol2") != std::string::npos)
        return Molecule(Mol22Mol(filename));
    else if (std::string(filename).find(".sdf") != std::string::npos)
        return Molecule(SDF2Mol(filename));
    else if (std::string(filename).find(".vtf") != std::string::npos)
        return Molecule(VTF2Mol(filename));
    else if (std::string(filename).find(".json") != std::string::npos) {
        Molecule molecule;
        molecule.ImportJson(filename);
        return molecule;
    } else if (std::string(filename).find("coord") != std::string::npos || std::string(filename).find("tmol") != std::string::npos)
        return Molecule(Coord2Mol(filename));
    else {
        fmt::print(fg(fmt::color::salmon) | fmt::emphasis::bold, "\nI dont understand the file type. Please use xyz (trj), sdf, mol2, vtf or turbomole coord files as input.\n");
        fmt::print(fg(fmt::color::salmon) | fmt::emphasis::bold, "\nTried to open " + filename + " and failed.");
        fmt::print("\n\n");
        return Molecule();
    }
}

inline Molecule LoadMol(const std::string& filename)
{
    if (std::string(filename).find(".xyz") != std::string::npos || std::string(filename).find(".trj") != std::string::npos)
        return (XYZ2Mol(filename));
    else if (std::string(filename).find(".mol2") != std::string::npos)
        return (Mol22Mol(filename));
    else if (std::string(filename).find(".sdf") != std::string::npos)
        return (SDF2Mol(filename));
    else if (std::string(filename).find(".vtf") != std::string::npos)
        return (VTF2Mol(filename));
    else if (std::string(filename).find("coord") != std::string::npos || std::string(filename).find("tmol") != std::string::npos)
        return (Coord2Mol(filename));
    else {
        fmt::print(fg(fmt::color::salmon) | fmt::emphasis::bold, "\nI dont understand the file type. Please use xyz (trj), sdf, mol2, vtf or turbomole coord files as input.\n");
        fmt::print(fg(fmt::color::salmon) | fmt::emphasis::bold, "\nTried to open " + filename + " and failed.");
        fmt::print("\n\n");
        return Molecule();
    }
}
/*
inline void xyz2allxyz(const std::string& xyzfile)
{
    std::string allxyz = xyzfile;
    allxyz.erase(allxyz.end() - 3, allxyz.end());

    std::ofstream input;
    input.open(allxyz + "allxyz", std::ios_base::app);
    FileIterator file(xyzfile);
    while (!file.AtEnd()) {
        Molecule mol = file.Next();
        input << mol.XYZString();
        if (!file.AtEnd())
            input << ">" << std::endl;
    }
    //input << "\n";
    input.close();
}
*/

}
