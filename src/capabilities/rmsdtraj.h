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
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

#include "curcumamethod.h"

class RMSDDriver;

#include "json.hpp"
using json = nlohmann::json;

/* Claude Generated 2025: RMSDTraj Parameter Registry - replaces static RMSDTrajJson */
BEGIN_PARAMETER_DEFINITION(rmsdtraj)
    PARAM(heavy_only, Bool, false, "Use only heavy atoms.", "RMSD", {"heavy"})
    PARAM(rmsd_threshold, Double, 1.5, "RMSD threshold for clustering (Å).", "RMSD", {"rmsd"})
    PARAM(reference, String, "", "Reference structure file.", "Input", {})
    PARAM(second_trajectory, String, "", "Second trajectory file for comparison.", "Input", {"second"})
    PARAM(write_unique, Bool, false, "Write unique conformers to file.", "Output", {"writeUnique"})
    PARAM(write_aligned, Bool, false, "Write aligned trajectory.", "Output", {"writeAligned"})
    PARAM(write_rmsd, Bool, true, "Write RMSD values to a file.", "Output", {"writeRMSD"})
    PARAM(fragment, Int, -1, "Fragment for RMSD calculation (-1: all).", "RMSD", {})
    PARAM(pca_file, Bool, false, "Generate PCA output file.", "Analysis", {"pcafile"})
    PARAM(all_xyz, Bool, false, "Write all structures to separate files.", "Output", {"allxyz"})
    PARAM(ref_first, Bool, false, "Use first structure as reference.", "Algorithm", {"RefFirst"})
    PARAM(noreorder, Bool, true, "Disable atom reordering.", "Algorithm", {})
    PARAM(optimize, Bool, false, "Optimize structures before RMSD.", "Algorithm", {"opt"})
    PARAM(filter, Bool, false, "Apply structure filtering.", "Analysis", {})
    PARAM(offset, Int, 0, "Skip first N frames in trajectory.", "Input", {})
END_PARAMETER_DEFINITION

class RMSDTraj : public CurcumaMethod {
public:
    /**
     * @brief Constructor with JSON configuration (backward compatible)
     * Claude Generated: Phase 4 - ConfigManager Migration
     */
    RMSDTraj(const json& controller, bool silent = true);

    /**
     * @brief Constructor with ConfigManager configuration (new, preferred)
     * Claude Generated: Phase 4 - Native ConfigManager support
     */
    RMSDTraj(const ConfigManager& config, bool silent = true);
    virtual ~RMSDTraj();

    virtual bool Initialise() override;

    void setBaseName(const std::string& name) { setFile(name); }
    void setFile(const std::string& filename) override { CurcumaMethod::setFile(filename); }
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

    /**
     * @brief Prints detailed help information about the RMSDTraj module parameters
     */
    void printHelp() const
    {
        std::cout << "\n=== RMSDTraj: Trajectory RMSD Analysis Tool ===\n\n"
                  << "This module analyzes molecular dynamics trajectories by calculating RMSD values\n"
                  << "between structures, identifying unique conformers, and performing related analyses.\n\n"
                  << "Parameter           | Default     | Description\n"
                  << "-------------------|-------------|----------------------------------------------------\n"
                  << "writeUnique        | " << std::setw(11) << (m_defaults.at("writeUnique") ? "true" : "false") << " | Write only unique conformers to output file\n"
                  << "rmsd               | " << std::setw(11) << m_defaults.at("rmsd") << " | RMSD threshold for unique conformer detection (Å)\n"
                  << "writeAligned       | " << std::setw(11) << (m_defaults.at("writeAligned") ? "true" : "false") << " | Write aligned structures to output file\n"
                  << "reference          | " << std::setw(11) << m_defaults.at("reference") << " | Reference structure file for RMSD calculations\n"
                  << "second             | " << std::setw(11) << m_defaults.at("second") << " | Second trajectory file for pairwise comparison\n"
                  << "fragment           | " << std::setw(11) << m_defaults.at("fragment") << " | Fragment index for RMSD calculation (-1: use entire molecule)\n"
                  << "heavy              | " << std::setw(11) << (m_defaults.at("heavy") ? "true" : "false") << " | Only use heavy atoms (non-hydrogen) for RMSD\n"
                  << "pcafile            | " << std::setw(11) << (m_defaults.at("pcafile") ? "true" : "false") << " | Generate Principal Component Analysis output file\n"
                  << "allxyz             | " << std::setw(11) << (m_defaults.at("allxyz") ? "true" : "false") << " | Write all structures to separate XYZ files\n"
                  << "RefFirst           | " << std::setw(11) << (m_defaults.at("RefFirst") ? "true" : "false") << " | Use first structure in trajectory as reference\n"
                  << "noreorder          | " << std::setw(11) << (m_defaults.at("noreorder") ? "true" : "false") << " | Disable atom reordering during RMSD calculation\n"
                  << "opt                | " << std::setw(11) << (m_defaults.at("opt") ? "true" : "false") << " | Optimize structures before RMSD calculation\n"
                  << "filter             | " << std::setw(11) << (m_defaults.at("filter") ? "true" : "false") << " | Apply filtering to select structures\n"
                  << "writeRMSD          | " << std::setw(11) << (m_defaults.at("writeRMSD") ? "true" : "false") << " | Write RMSD values to output file\n"
                  << "offset             | " << std::setw(11) << m_defaults.at("offset") << " | Number of initial frames to skip in trajectory\n"
                  << "\n=== Output Files ===\n\n"
                  << "Several output files are generated with the base name derived from the input file:\n"
                  << "- [basename]_rmsd.dat     : RMSD values between structures\n"
                  << "- [basename]_unique.xyz   : Unique conformers (if writeUnique is true)\n"
                  << "- [basename]_aligned.xyz  : Aligned structures (if writeAligned is true)\n"
                  << "- [basename]_pca.dat      : PCA projection data (if pcafile is true)\n"
                  << "- [basename]_pairwise.dat : Pairwise RMSD comparison (if second file is provided)\n"
                  << "\n=== Common Usage Examples ===\n\n"
                  << "1. Extract unique conformers from a trajectory:\n"
                  << "   curcuma -rmsdtraj trajectory.xyz -writeUnique -rmsd 1.0\n\n"
                  << "2. Compare all frames to a reference structure:\n"
                  << "   curcuma -rmsdtraj trajectory.xyz -reference ref.xyz -writeRMSD\n\n"
                  << "3. Analyze only heavy atoms and generate PCA data:\n"
                  << "   curcuma -rmsdtraj trajectory.xyz -heavy -pcafile\n\n"
                  << "4. Compare two trajectories pairwise:\n"
                  << "   curcuma -rmsdtraj traj1.xyz -second traj2.xyz\n\n"
                  << "5. Focus analysis on a specific fragment:\n"
                  << "   curcuma -rmsdtraj trajectory.xyz -fragment 1\n\n"
                  << "Note: When comparing large trajectories, consider using the offset parameter\n"
                  << "to skip initial frames that may represent equilibration periods.\n"
                  << std::endl;
    }

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

    void ProcessSingleFile();
    void CompareTrajectories();

    std::string m_reference, m_second_file, m_outfile;
    std::ofstream m_rmsd_file, m_pca_file, m_pairwise_file;
    std::vector<Molecule*> m_stored_structures;
    Molecule *m_initial, *m_previous;
    RMSDDriver* m_driver;
    std::vector<double> m_rmsd_vector, m_energy_vector;
    int m_fragment = -1;
    int m_currentIndex = 0;
    int m_atoms = -1;
    int m_max_lines = -1;
    int m_offset = 0;
    bool m_writeUnique = false, m_pairwise = false, m_heavy = false, m_pcafile = false, m_writeAligned = false, m_ref_first = false, m_opt = false, m_filter = false, m_writeRMSD = true;
    bool m_allxyz = false;
    double m_rmsd_threshold = 1.0;
};
