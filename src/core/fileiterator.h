/*
 * <Load xyz files and iterate through them.>
 * Copyright (C) 2020 - 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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
    FileIterator(char* filename, bool silent = false);
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

    std::string m_filename, m_basename;
    std::ifstream* m_file;
    bool m_end = false, m_init = false;
    Molecule m_current;
    int m_lines = 0, m_current_mol = 0, m_mols = 0;
};
