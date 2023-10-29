/*
 * <Geometry tools for chemical structures.>
 * Copyright (C) 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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
            molecule.m_geometry.push_back(pair.second);
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
            } catch (const std::invalid_argument& arg) {
                atoms = 0;
            }
        }

        if (i == 1)
            molecule.m_commentline = line;
        if (i > 1) {
            auto pair = Line2Atoms(line);
            molecule.m_atoms.push_back(pair.first);
            molecule.m_geometry.push_back(pair.second);
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
                molecule.m_geometry.push_back(vector);
            } catch (const std::invalid_argument& arg) {
            }
        }
        //}
        //if(line.compare("@<TRIPOS>ATOM") == 0)
        //    readblock = true;
        //if(line.compare("@<TRIPOS>BOND") == 0)
        //    readblock = false;
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
                molecule.m_geometry.push_back(vector);
            } catch (const std::invalid_argument& arg) {
            }
        }
        //}
        //if(line.compare("@<TRIPOS>ATOM") == 0)
        //    readblock = true;
        //if(line.compare("@<TRIPOS>BOND") == 0)
        //    readblock = false;
    }
    return molecule;
}

inline Mol Mol22Mol(const std::string& filename)
{
    Mol molecule;
    int readblock = false;
    auto file = new std::ifstream(filename);
    for (std::string line; getline(*file, line);) {
        if (readblock) {
            auto strings = SplitString(line, " ");
            if (strings.size() == 9) {
                try {
                    molecule.m_atoms.push_back(Elements::String2Element(strings[1]));
                    std::array<double, 3> vector;
                    vector[0] = std::stod(strings[2]);
                    vector[1] = std::stod(strings[3]);
                    vector[2] = std::stod(strings[4]);
                    molecule.m_geometry.push_back(vector);
                } catch (const std::invalid_argument& arg) {
                }
            }
        }
        if (line.compare("@<TRIPOS>ATOM") == 0)
            readblock = true;
        if (line.compare("@<TRIPOS>BOND") == 0)
            readblock = false;
    }
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
    else if (std::string(filename).find(".json") != std::string::npos) {
        Molecule molecule;
        molecule.ImportJson(filename);
        return molecule;
    } else if (std::string(filename).find("coord") != std::string::npos || std::string(filename).find("tmol") != std::string::npos)
        return Molecule(Coord2Mol(filename));
    else {
        fmt::print(fg(fmt::color::salmon) | fmt::emphasis::bold, "\nI dont understand the file type. Please use xyz (trj), sdf, mol2 or turbomole coord files as input.\n");
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
    else if (std::string(filename).find("coord") != std::string::npos || std::string(filename).find("tmol") != std::string::npos)
        return (Coord2Mol(filename));
    else {
        fmt::print(fg(fmt::color::salmon) | fmt::emphasis::bold, "\nI dont understand the file type. Please use xyz (trj), sdf, mol2 or turbomole coord files as input.\n");
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
