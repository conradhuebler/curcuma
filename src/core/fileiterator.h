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

#include <fstream>
#include <iostream>
#include <string>

class FileIterator {
public:
    inline FileIterator(const std::string& filename)
        : m_filename(filename)
    {
        std::cerr << "Opening file " << m_filename << std::endl;
        m_file = new std::ifstream(m_filename);
        m_lines = CountLines();
        m_init = CheckNext();
    }

    inline FileIterator(char* filename)
    {
        m_filename = std::string(filename);
        std::cerr << "Opening file " << m_filename << std::endl;
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

private:
    bool CheckNext()
    {
        std::vector<std::string> lines;
        int atoms = 0;
        int index = 0;
        int i = 0;
        bool xyzfile = std::string(m_filename).find(".xyz") != std::string::npos || std::string(m_filename).find(".trj") != std::string::npos;
        Molecule mol(atoms, 0);
        for (std::string line; getline(*m_file, line);) {
            if (index == 0 && xyzfile) {
                atoms = stoi(line);
                m_mols = m_lines / (atoms + 2);
                mol = Molecule(atoms, 0);
                m_current_mol++;
            }
            if (xyzfile) {
                if (i == 1)
                    mol.setXYZComment(line);
                if (i > 1) {
                    mol.setXYZ(line, i - 2);
                }
                if (i - 1 == atoms) {
                    m_current = mol;
                    return false;
                }
                ++i;
            } else {
                mol.setAtom(line, i);
            }
            index++;
        }
        return true;
    }

    inline int CountLines() const
    {
        std::ifstream inFile(m_filename);
        return std::count(std::istreambuf_iterator<char>(inFile),
            std::istreambuf_iterator<char>(), '\n');
    }

    std::string m_filename;
    std::ifstream* m_file;
    bool m_end = false, m_init = false;
    Molecule m_current;
    int m_lines = 0, m_current_mol = 0, m_mols = 0;
};
