/*
 * <Load xyz files and iterate through them.>
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

#include "src/core/molecule.h"

#include "src/tools/formats.h"
#include "src/tools/general.h"

#include <fstream>
#include <iostream>
#include <string>

#include "fileiterator.h"

FileIterator::FileIterator(bool silent)
{
}

FileIterator::FileIterator(const std::string& filename, bool silent)
    : m_filename(filename)
{
    if (!silent)
        std::cerr << "Opening file " << m_filename << std::endl;
    m_basename = filename;
    m_basename.erase(m_basename.end() - 4, m_basename.end());
    m_file = new std::ifstream(m_filename);
    m_lines = CountLines();

    // Check if this is a VTF file - Claude Generated
    m_is_vtf_file = std::string(m_filename).find(".vtf") != std::string::npos;
    if (m_is_vtf_file) {
        m_mols = CountVTFTimesteps();
        ParseVTFHeader();
    }

    m_init = CheckNext();
}

FileIterator::FileIterator(char* filename, bool silent)
{
    m_filename = std::string(filename);
    if (!silent)
        std::cerr << "Opening file " << m_filename << std::endl;
    m_basename = std::string(filename);
    m_basename.erase(m_basename.end() - 4, m_basename.end());
    m_file = new std::ifstream(m_filename);
    m_lines = CountLines();

    // Check if this is a VTF file - Claude Generated
    m_is_vtf_file = std::string(m_filename).find(".vtf") != std::string::npos;
    if (m_is_vtf_file) {
        m_mols = CountVTFTimesteps();
        ParseVTFHeader();
    }

    m_init = CheckNext();
}

void FileIterator::setFile(const std::string& filename)
{
    m_filename = filename;
    m_basename = filename;
    m_basename.erase(m_basename.end() - 4, m_basename.end());
    m_file = new std::ifstream(m_filename);
    m_lines = CountLines();

    // Check if this is a VTF file - Claude Generated
    m_is_vtf_file = std::string(m_filename).find(".vtf") != std::string::npos;
    if (m_is_vtf_file) {
        m_mols = CountVTFTimesteps();
        ParseVTFHeader();
    }

    m_init = CheckNext();
}
FileIterator::~FileIterator()
{
}

Molecule FileIterator::Next()
{
    Molecule current = m_current;
    m_end = CheckNext();
    return current;
}

bool FileIterator::AtEnd()
{
    return m_end || m_init;
}
Molecule FileIterator::Current() const
{
    return m_current;
}

int FileIterator::MaxMolecules() const { return m_mols; }

int FileIterator::CurrentMolecule() const { return m_current_mol; }

std::string FileIterator::Basename() const { return m_basename; }

bool FileIterator::CheckNext()
{
    std::vector<std::string> lines;
    int atoms = 0;
    int index = 0;
    int i = 0;
    bool xyzfile = std::string(m_filename).find(".xyz") != std::string::npos || std::string(m_filename).find(".trj") != std::string::npos;

    // VTF file handling - Claude Generated
    if (m_is_vtf_file) {
        return ParseVTFTimestep();
    } else if (xyzfile) {
        Molecule mol(atoms, 0);
        for (std::string line; getline(*m_file, line);) {
            if (line.size() == 0 && i != 1)
                continue;
            if (index == 0 && xyzfile) {
                try {
                    atoms = stoi(line);
                } catch (const std::invalid_argument& arg) {
                    std::cerr << "FileIterator::CheckNext() Got some error at line " << line << "\n";
                    std::cerr << "Skipping molecules that follow after  " << m_current_mol << " molecule!" << std::endl;
                    return false;
                }
                m_mols = m_lines / (atoms + 2);
                mol = Molecule(atoms, 0);
                m_current_mol++;
            }
            if (xyzfile) {
                if (i == 1)
                    mol.setXYZComment(line);
                if (i > 1) {
                    try {
                        mol.setXYZ(line, i - 2);
                    } catch (const std::invalid_argument& arg) {
                        std::cerr << "FileIterator::CheckNext() Got some error at line " << line << "\n";
                        std::cerr << "Skipping molecules that follow after  " << m_current_mol << " molecule!" << std::endl;
                        return false;
                    }
                }
                if (i - 1 == atoms) {
                    mol.CalculateMass(); // Claude Generated 2025: Calculate mass after loading all atoms (fixes XYZ m_mass=0 bug)
                    m_current = mol;
                    index = 0;
                    return false;
                }
                ++i;
            } else {
                mol.setAtom(line, i);
            }
            index++;
        }
    } else {
        m_current = Files::LoadFile(m_filename);
        m_init = true;
        return false;
    }
    return true;
}

int FileIterator::CountLines() const
{
    std::ifstream inFile(m_filename);
    return std::count(std::istreambuf_iterator<char>(inFile),
        std::istreambuf_iterator<char>(), '\n');
}

// VTF-specific methods implementation - Claude Generated
int FileIterator::CountVTFTimesteps() const
{
    std::ifstream inFile(m_filename);
    int timestep_count = 0;
    for (std::string line; getline(inFile, line);) {
        if (line.find("timestep") == 0) {
            timestep_count++;
        }
    }
    return timestep_count;
}

bool FileIterator::ParseVTFHeader()
{
    // Reset file to beginning
    m_file->clear();
    m_file->seekg(0, std::ios::beg);

    // Use existing VTF2Mol to parse the structure (first timestep only)
    m_vtf_template = Files::VTF2Mol(m_filename);
    m_vtf_atom_count = m_vtf_template.AtomCount();

    // Reset file position after parsing
    m_file->clear();
    m_file->seekg(0, std::ios::beg);

    return m_vtf_atom_count > 0;
}

bool FileIterator::ParseVTFTimestep()
{
    std::string line;
    bool found_timestep = false;
    int coord_count = 0;

    // Create a copy of the template molecule
    Molecule current_mol = m_vtf_template;

    // Skip to next timestep
    while (getline(*m_file, line)) {
        if (line.find("timestep") == 0) {
            found_timestep = true;
            m_current_mol++;
            break;
        }
    }

    if (!found_timestep) {
        return true; // End of file reached
    }

    // Read coordinates for this timestep
    Geometry new_geometry = current_mol.getGeometry();
    for (int i = 0; i < m_vtf_atom_count && getline(*m_file, line); i++) {
        // Parse coordinate line: x y z
        std::vector<std::string> tokens = Files::SplitString(line, " ");
        if (tokens.size() >= 3) {
            try {
                double x = std::stod(tokens[0]);
                double y = std::stod(tokens[1]);
                double z = std::stod(tokens[2]);

                // Update geometry matrix
                new_geometry(i, 0) = x;
                new_geometry(i, 1) = y;
                new_geometry(i, 2) = z;
            } catch (const std::invalid_argument& arg) {
                std::cerr << "FileIterator::ParseVTFTimestep() Error parsing coordinate line: " << line << std::endl;
                return true; // Error - end iteration
            }
        }
    }

    // Set the updated geometry
    current_mol.setGeometry(new_geometry);

    m_current = current_mol;
    return false; // Successfully parsed timestep
}
