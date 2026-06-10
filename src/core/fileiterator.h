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

class FileIterator {
public:
    FileIterator(bool silent = false);
    FileIterator(const std::string& filename, bool silent = false);
    // const char* (not char*): a const char* such as std::string::c_str() cannot
    // bind to char* and would otherwise resolve to FileIterator(bool) (pointer ->
    // bool conversion), silently default-constructing an unusable iterator with an
    // empty filename. That bug caused an infinite LoadFile("") loop in ConfSearch.
    FileIterator(const char* filename, bool silent = false);
    ~FileIterator();

    void setFile(const std::string& filename);

    Molecule Next();

    bool AtEnd();
    Molecule Current() const;

    int MaxMolecules() const;

    int CurrentMolecule() const;

    std::string Basename() const;

private:
    bool CheckNext();

    int CountLines() const;

    // VTF-specific methods - Claude Generated
    int CountVTFTimesteps() const;
    bool ParseVTFHeader();
    bool ParseVTFTimestep();

    std::string m_filename, m_basename;
    std::ifstream* m_file = nullptr;
    bool m_end = false, m_init = false;
    bool m_nonxyz_loaded = false; // non-XYZ single-molecule file already delivered
    Molecule m_current;
    int m_lines = 0, m_current_mol = 0, m_mols = 0;

    // VTF-specific member variables - Claude Generated
    bool m_is_vtf_file = false;
    Molecule m_vtf_template;
    int m_vtf_atom_count = 0;
};
