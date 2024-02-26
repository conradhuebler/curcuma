/*
 * <Load xyz files and iterate through them.>
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

#include "src/core/molecule.h"

#include "src/tools/formats.h"
#include "src/tools/general.h"

#include <fstream>
#include <iostream>
#include <string>

class FileIterator {
public:
    inline FileIterator(bool silent = false)

    {
    }

    inline FileIterator(const std::string& filename, bool silent = false)
        : m_filename(filename)
    {
        if (!silent)
            std::cerr << "Opening file " << m_filename << std::endl;
        m_basename = filename;
        m_basename.erase(m_basename.end() - 4, m_basename.end());
        m_file = new std::ifstream(m_filename);
        m_lines = CountLines();
        m_init = CheckNext();
    }

    inline FileIterator(char* filename, bool silent = false)
    {
        m_filename = std::string(filename);
        if (!silent)
            std::cerr << "Opening file " << m_filename << std::endl;
        m_basename = std::string(filename);
        ;
        m_basename.erase(m_basename.end() - 4, m_basename.end());
        m_file = new std::ifstream(m_filename);
        m_lines = CountLines();
        m_init = CheckNext();
    }

    inline void setFile(const std::string& filename)
    {
        m_filename = filename;
        m_basename = filename;
        m_basename.erase(m_basename.end() - 4, m_basename.end());
        m_file = new std::ifstream(m_filename);
        m_lines = CountLines();
        m_init = CheckNext();
    }

    inline Molecule Next()
    {
        Molecule current = m_current;
        m_end = CheckNext();
        return current;
    }

    inline bool AtEnd()
    {
        return m_end || m_init;
    }
    inline Molecule Current() const
    {
        return m_current;
    }

    inline int MaxMolecules() const { return m_mols; }

    inline int CurrentMolecule() const { return m_current_mol; }

    inline std::string Basename() const { return m_basename; }

private:
    bool CheckNext()
    {
        std::vector<std::string> lines;
        int atoms = 0;
        int index = 0;
        int i = 0;
        bool xyzfile = std::string(m_filename).find(".xyz") != std::string::npos || std::string(m_filename).find(".trj") != std::string::npos;
        if (xyzfile) {
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

    inline int CountLines() const
    {
        std::ifstream inFile(m_filename);
        return std::count(std::istreambuf_iterator<char>(inFile),
            std::istreambuf_iterator<char>(), '\n');
    }

    std::string m_filename, m_basename;
    std::ifstream* m_file;
    bool m_end = false, m_init = false;
    Molecule m_current;
    int m_lines = 0, m_current_mol = 0, m_mols = 0;
};
