/*
 * <Trajectory RMSD Analyse. >
 * Copyright (C) 2020 - 2022 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "src/core/molecule.h"

#include "curcumamethod.h"

class RMSDDriver;

#include "json.hpp"
using json = nlohmann::json;

const json RMSDTrajJson{
    { "writeUnique", false },
    { "writeAligned", false },
    { "rmsd", 1.5 },
    { "fragment", -1 },
    { "reference", "none" },
    { "second", "none" },
    { "heavy", false },
    { "pcafile", false },
    { "allxyz", false },
    { "RefFirst", false },
    { "noreorder", true },
    { "opt", false },
    { "filter", false },
    { "writeRMSD", true }
};

class RMSDTraj : public CurcumaMethod {
public:
    RMSDTraj(const json& controller = RMSDTrajJson, bool silent = true);
    virtual ~RMSDTraj();

    virtual bool Initialise();

    void setBaseName(const std::string& name) { m_filename = name; }
    void setFile(const std::string& filename) { m_filename = filename; }
    void setSecondFile(const std::string& filename)
    {
        m_second_file = filename;
        m_pairwise = true;
    }

    void setReferenceStructure(const std::string& reference) { m_reference = reference; }

    /*! \brief Set the index of the fragment that is used for rmsd calculation/atom reordering */
    inline void setFragment(int fragment) { m_fragment = fragment; }

    inline void setRMSDThreshold(double rmsd_threshold) { m_rmsd_threshold = rmsd_threshold; }

    inline void WriteUnique(bool write_unique) { m_writeUnique = write_unique; }
    inline void setHeavy(bool heavy) { m_heavy = heavy; }

    void start() override;

    bool CheckMolecule(Molecule* molecule);

    int StoredStructures() const { return m_stored_structures.size(); }
    void PostAnalyse();

    void Optimise();

    void Filter();

private:
    /* Read Controller has to be implemented for all */
    void LoadControlJson() override;

    /* Lets have this for all modules */
    nlohmann::json WriteRestartInformation() override { return json(); }

    /* Lets have this for all modules */
    bool LoadRestartInformation() override { return true; }

    StringList MethodName() const override { return { std::string("RMSDTraj") }; }

    /* Lets have all methods read the input/control file */
    void ReadControlFile() override {}

    std::string m_filename, m_reference, m_second_file, m_outfile;
    std::ofstream m_rmsd_file, m_pca_file, m_pairwise_file;
    std::vector<Molecule*> m_stored_structures;
    Molecule *m_initial, *m_previous;
    RMSDDriver* m_driver;
    std::vector<double> m_rmsd_vector, m_energy_vector;
    int m_fragment = -1;
    int m_currentIndex = 0;
    int m_atoms = -1;
    int m_max_lines = -1;
    bool m_writeUnique = false, m_pairwise = false, m_heavy = false, m_pcafile = false, m_writeAligned = false, m_ref_first = false, m_opt = false, m_filter = false, m_writeRMSD = true;
    bool m_allxyz = false;
    double m_rmsd_threshold = 1.0;
};
