/*
 * <Trajectory RMSD Analyse. >
 * Copyright (C) 2020 - 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "curcumamethod.h"
#include "src/global_config.h"

#include "src/capabilities/confscan.h"
#include "src/capabilities/curcumaopt.h"
#include "src/capabilities/rmsd.h"

#include "src/core/elements.h"
#include "src/core/fileiterator.h"
#include "src/core/molecule.h"

#include "src/tools/formats.h"
#include "src/tools/general.h"
#include "src/tools/geometry.h"

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "json.hpp"
using json = nlohmann::json;

#include "rmsdtraj.h"

RMSDTraj::RMSDTraj(const json& controller, bool silent)
    : CurcumaMethod(RMSDTrajJson, controller, silent)
{
    UpdateController(controller);
}
RMSDTraj::~RMSDTraj()
{
    for (int i = 0; i < m_stored_structures.size(); ++i)
        delete m_stored_structures[i];
    delete m_driver;
}

bool RMSDTraj::Initialise()
{
    // Synchronize CurcumaLogger verbosity with local verbosity
    CurcumaLogger::set_verbosity(m_verbosity);

    CurcumaLogger::header("Trajectory RMSD Analysis");

    // Enhanced parameter display with more context
    CurcumaLogger::param("Input trajectory", m_filename);
    CurcumaLogger::param("Write conformers", m_writeUnique ? "Yes" : "No");
    if (m_writeUnique) {
        CurcumaLogger::param("RMSD threshold", fmt::format("{:.3f} Å", m_rmsd_threshold));
    }
    CurcumaLogger::param("Write RMSD data", m_writeRMSD ? "Yes" : "No");
    CurcumaLogger::param("Write aligned structures", m_writeAligned ? "Yes" : "No");
    if (m_reference.compare("none") != 0) {
        CurcumaLogger::param("Reference structure", m_reference);
        m_stored_structures.push_back(new Molecule(Files::LoadFile(m_reference)));
        m_atoms = m_stored_structures[0]->AtomCount();
        CurcumaLogger::param("Reference atoms", std::to_string(m_atoms));
    }

    // Modern C++ way to remove file extension
    m_outfile = m_filename;
    size_t last_dot = m_outfile.find_last_of('.');
    if (last_dot != std::string::npos) {
        m_outfile = m_outfile.substr(0, last_dot);
    }

    // Initialize output files with user feedback
    if (m_writeRMSD) {
        m_rmsd_file.open(m_outfile + ".rmsd.dat");
        CurcumaLogger::success_fmt("RMSD data will be written to: {}.rmsd.dat", m_outfile);
    }

    if (m_pcafile) {
        m_pca_file.open(m_outfile + ".pca.dat");
        CurcumaLogger::success_fmt("PCA data will be written to: {}.pca.dat", m_outfile);
    }

    if (m_pairwise) {
        m_pairwise_file.open(m_outfile + ".pairwise.dat");
        CurcumaLogger::success_fmt("Pairwise RMSD will be written to: {}.pairwise.dat", m_outfile);
    }

    if (m_writeUnique) {
        CurcumaLogger::success_fmt("Unique conformers will be written to: {}.unique.xyz", m_outfile);
    }

    json RMSDJsonControl = {
        { "reorder", false },
        { "check", false },
        { "heavy", false },
        { "fragment", -1 },
        { "fragment_reference", -1 },
        { "fragment_target", -1 },
        { "init", -1 },
        { "pt", 0 },
        { "silent", true },
        { "storage", 1.0 },
        { "method", "incr" },
        //{ "noreorder", m_noreorder },
        { "threads", 1 }
    };
    m_driver = new RMSDDriver(RMSDJsonControl);
    m_driver->setProtons(!m_heavy);
    m_driver->setForceReorder(false);
    m_driver->setCheckConnections(false);
    m_driver->setFragment(m_fragment);
    std::ofstream export_file;
    if (m_writeUnique) {
        export_file.open(m_outfile + ".unique.xyz");
        export_file.close();
    }
    if (m_writeAligned) {
        export_file.open(m_outfile + ".aligned.xyz");
        export_file.close();
    }
    std::ifstream input(m_filename);
    std::vector<std::string> lines;
    //    int atoms = 0, atoms2 = 0;
    m_currentIndex = 0;
    //  int i = 0;
    //    int molecule = 0;
    std::ifstream inFile(m_filename);
    m_max_lines = std::count(std::istreambuf_iterator<char>(inFile),
        std::istreambuf_iterator<char>(), '\n');

    std::ifstream second;
    if (m_pairwise) {
        second.open(m_second_file);
    }
    return true;
}
void RMSDTraj::start()
{
    // Ensure CurcumaLogger verbosity is synchronized throughout execution
    CurcumaLogger::set_verbosity(m_verbosity);

    if (m_second_file.compare("none") == 0)
        ProcessSingleFile();
    else
        CompareTrajectories();
}

void RMSDTraj::ProcessSingleFile()
{
    auto start_time = std::chrono::high_resolution_clock::now();

    CurcumaLogger::success("Counting trajectory frames...");

    // Count actual molecules, not lines - this is more accurate for XYZ trajectories
    m_max_lines = 0;
    {
        FileIterator counter(m_filename);
        while (!counter.AtEnd()) {
            counter.Next(); // Skip the molecule data, just count
            m_max_lines++;
        }
    }

    CurcumaLogger::success_fmt("Found {} trajectory frames to analyze", m_max_lines);

    Molecule prev;
    FileIterator file(m_filename);
    // Improved progress tracking with better efficiency
    std::vector<bool> progress_reported(21, false); // 5% increments
    int structures_processed = 0;
    int structures_accepted = 0;
    int last_percent_reported = -1;

    CurcumaLogger::success("Starting trajectory analysis...");

    // Better debugging information
    if (m_verbosity >= 2) {
        CurcumaLogger::info_fmt("FileIterator initialized for: {}", m_filename);
        CurcumaLogger::info_fmt("Expected {} lines to process", m_max_lines);
    }

    while (!file.AtEnd()) {
        auto molecule = std::make_unique<Molecule>(file.Next());
        structures_processed++;

        // Debug: Always show first few structures at high verbosity
        if (structures_processed <= 5 && m_verbosity >= 3) {
            CurcumaLogger::info_fmt("Processing structure #{}", structures_processed);
        }

        bool check = CheckMolecule(molecule.get());
        if (check) {
            structures_accepted++;
            // Only report significant milestones, not every single structure
            if (structures_accepted % 10 == 1 || (structures_accepted < 10)) {
                CurcumaLogger::success_fmt("Structure #{} accepted (total: {})",
                    structures_processed, structures_accepted);
                std::cout.flush();
            }
        } else if (structures_processed <= 10 && m_verbosity >= 3) {
            CurcumaLogger::warn_fmt("Structure #{} rejected", structures_processed);
        }

        // More efficient progress reporting - 5% increments based on actual progress
        int current_percent = int((structures_processed / double(m_max_lines)) * 100);
        int progress_step = current_percent / 5; // 5% increments

        if (progress_step < progress_reported.size() && !progress_reported[progress_step] && current_percent != last_percent_reported && current_percent >= 5) {

            progress_reported[progress_step] = true;
            last_percent_reported = current_percent;

            CurcumaLogger::progress(structures_processed, m_max_lines,
                fmt::format("Processing trajectory - {}% done, {} structures accepted",
                    current_percent, structures_accepted));
            std::cout.flush();
        }

        if (CheckStop())
            break;
    }

    // Final timing and summary
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    CurcumaLogger::success_fmt("Trajectory analysis completed: {} structures processed, {} accepted",
        structures_processed, structures_accepted);
    CurcumaLogger::info_fmt("Analysis time: {:.3f} seconds", duration.count() / 1000.0);

    if (structures_processed > 0) {
        double acceptance_rate = (structures_accepted / double(structures_processed)) * 100.0;
        CurcumaLogger::param("Acceptance rate", fmt::format("{:.1f}%", acceptance_rate));

        if (duration.count() > 0) {
            double throughput = structures_processed / (duration.count() / 1000.0);
            CurcumaLogger::param("Throughput", fmt::format("{:.1f} structures/second", throughput));
        }
    }

    PostAnalyse();

    if (m_opt && m_writeUnique && !CheckStop()) {
        Optimise();
        if (m_filter)
            Filter();
    }
}

bool RMSDTraj::CheckMolecule(Molecule* molecule)
{
    bool result = false;
    // std::cout << molecule->Atom(0).second.transpose() << std::endl << std::endl;

    double energy = molecule->Energy();
    // Note: m_currentIndex removed - using structures_processed for progress instead

    if (m_stored_structures.size() == 0) {
        if (m_writeUnique) {
            molecule->Center();
            molecule->appendXYZFile(m_outfile + ".unique.xyz");
            // std::cout << "First structure added!" << std::endl;
            result = true;
        }
        m_stored_structures.push_back(new Molecule(molecule));
        m_initial = molecule;
        m_previous = molecule;
        return result;
    } else {
        m_previous = m_stored_structures[0];
        m_initial = m_previous;
        /*
        for (std::size_t i = 0; i < mol.GetFragments().size(); ++i)
            if (mol.getGeometryByFragment(i).rows() == atoms_target) {
                driver->setFragmentTarget(i);
                driver->setPartialRMSD(true);
            }
            */
    }
    m_driver->setScaling(1.3);
    if (m_pairwise == false) {
        std::vector<double> rmsd_results;

        if (!m_ref_first) // If we reference to the previouse (default) we have to pre-reorder and then calculate the rmsd with respect to the first structure and not the prevouise
        {
            m_driver->setReference(*m_previous);
            m_driver->setTarget(*molecule);
            m_driver->start();
            // m_previous = m_driver->TargetAligned();
            // molecule = m_driver->TargetAligned();
        }

        m_driver->setReference(*m_initial);
        m_driver->setTarget(*molecule);
        m_driver->start();

#ifdef CURCUMA_DEBUG
        if (m_driver->ReorderRules().size() && CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info_fmt("Reorder Rules: {}", Tools::Vector2String(m_driver->ReorderRules()));
        }
#endif

        m_rmsd_file << m_driver->RMSD() << "\t" << std::setprecision(10) << energy << std::endl;
        m_rmsd_vector.push_back(m_driver->RMSD());
        m_energy_vector.push_back(energy);
        if (m_writeAligned) {
            m_driver->TargetAligned().appendXYZFile(m_outfile + "_aligned.xyz");
        }

        if (m_pcafile) {
            Molecule mol2 = m_driver->TargetAligned();
            for (std::size_t j = 0; j < mol2.AtomCount(); ++j) {
                if (mol2.Atom(j).first != 1)
                    m_pca_file << mol2.Atom(j).second(0) << " " << mol2.Atom(j).second(1) << " " << mol2.Atom(j).second(2);
            }
            m_pca_file << std::endl;
        }

        double first_rmsd = m_driver->RMSD();
        if (m_writeUnique) {
            bool perform_rmsd = true;
            //  std::cout << std::endl;
            for (std::size_t mols = m_stored_structures.size() - 1; mols >= 0 && perform_rmsd && mols <= m_stored_structures.size(); --mols) {
                // m_driver->clear();
                m_driver->setReference(*m_stored_structures[mols]);
                m_driver->setTarget(*molecule);
                // std::cout << molecule->Atom(0).second.transpose() << " " << m_stored_structures[mols]->Atom(0).second.transpose() << std::endl;
                m_driver->start();

                rmsd_results.push_back(m_driver->RMSD());
                perform_rmsd = m_driver->RMSD() > m_rmsd_threshold;
                //  std::cout << m_driver->RMSD() << " ";
            }
            perform_rmsd = first_rmsd > m_rmsd_threshold && perform_rmsd;
            rmsd_results.push_back(first_rmsd);

            if (perform_rmsd) {
                molecule->LoadMolecule(m_driver->TargetAlignedReference());
                m_stored_structures.push_back(new Molecule(molecule));
                molecule->appendXYZFile(m_outfile + ".unique.xyz");
                //                std::cout << "New structure added ... ( " << m_stored_structures.size() << "). " << /*  int(m_currentIndex / double(m_max_lines) * 100) << " % done ...!" << */ std::endl;
                result = true;
            }
        }
    } else {
        CurcumaLogger::warn("Pairwise is disabled for now");
        /*
        m_driver->setReference(mol);
        m_driver->setTarget(mol_2);
        m_driver->start();
        m_pairwise_file << m_driver->RMSD() << std::endl;
        */
    }
    m_driver->clear();
    // i = -1;
    // mol = Molecule(atoms, 0);
    molecule++;
    // mol.setName(std::to_string(molecule));
    if (m_pairwise) {
        /*
        mol_2 = Molecule(atoms, 0);
        mol_2.setName(std::to_string(molecule));
        */
    }
    return result;
}

void RMSDTraj::CompareTrajectories()
{
    // Ensure CurcumaLogger verbosity is synchronized
    CurcumaLogger::set_verbosity(m_verbosity);

    CurcumaLogger::success("Comparing two trajectories...");

    FileIterator file1(m_filename);
    FileIterator file2(m_second_file);

    json RMSDJsonControl = {
        { "reorder", false },
        { "check", false },
        { "heavy", false },
        { "fragment", -1 },
        { "fragment_reference", -1 },
        { "fragment_target", -1 },
        { "init", -1 },
        { "pt", 0 },
        { "silent", true },
        { "storage", 1.0 },
        { "method", "inertia" },
        //{ "noreorder", m_noreorder },
        { "threads", 1 }
    };
    m_driver = new RMSDDriver(RMSDJsonControl);
    m_driver->setProtons(!m_heavy);
    m_driver->setForceReorder(false);
    m_driver->setCheckConnections(false);
    m_driver->setFragment(m_fragment);

    while (!file1.AtEnd() && !file2.AtEnd()) {
        Molecule* mol1 = new Molecule(file1.Next());
        Molecule* mol2 = new Molecule(file2.Next());

        //   std::cout << molecule->Atom(0).second.transpose() << std::endl;
        m_driver->setReference(*mol1);
        m_driver->setTarget(*mol2);
        m_driver->start();

        double rmsd = m_driver->RMSD();
        m_rmsd_vector.push_back(rmsd);
        CurcumaLogger::info_fmt("RMSD: {:.6f} Å", rmsd);
        delete mol1;
        delete mol2;
        if (CheckStop())
            break;
    }
    PostAnalyse();
}

void RMSDTraj::PostAnalyse()
{
    if (m_rmsd_vector.empty()) {
        CurcumaLogger::warn("No RMSD data available for analysis");
        return;
    }

    // Calculate comprehensive statistics
    double rmsd_mean = Tools::mean(m_rmsd_vector);
    double rmsd_median = Tools::median(m_rmsd_vector);
    double rmsd_std = Tools::stdev(m_rmsd_vector, rmsd_mean);
    auto rmsd_hist = Tools::Histogram(m_rmsd_vector, 100);
    double rmsd_shannon = Tools::ShannonEntropy(rmsd_hist);

    double energy_mean = 0.0, energy_median = 0.0, energy_std = 0.0, energy_shannon = 0.0;
    bool has_energy_data = !m_energy_vector.empty();

    if (has_energy_data) {
        energy_mean = Tools::mean(m_energy_vector);
        energy_median = Tools::median(m_energy_vector);
        energy_std = Tools::stdev(m_energy_vector, energy_mean);
        auto energy_hist = Tools::Histogram(m_energy_vector, 100);
        energy_shannon = Tools::ShannonEntropy(energy_hist);
    }

    // Display beautiful results to user
    CurcumaLogger::header("Trajectory Analysis Results");
    CurcumaLogger::info("");
    CurcumaLogger::success_fmt("Analyzed {} trajectory frames", m_rmsd_vector.size());
    CurcumaLogger::success_fmt("Unique conformers found: {}", m_stored_structures.size());

    CurcumaLogger::info("");
    CurcumaLogger::info("RMSD Statistics:");
    CurcumaLogger::param("Mean RMSD", fmt::format("{:.4f} Å", rmsd_mean));
    CurcumaLogger::param("Median RMSD", fmt::format("{:.4f} Å", rmsd_median));
    CurcumaLogger::param("Std. deviation", fmt::format("{:.4f} Å", rmsd_std));
    CurcumaLogger::param("Shannon entropy", fmt::format("{:.4f}", rmsd_shannon));
    CurcumaLogger::param("Min RMSD", fmt::format("{:.4f} Å", *std::min_element(m_rmsd_vector.begin(), m_rmsd_vector.end())));
    CurcumaLogger::param("Max RMSD", fmt::format("{:.4f} Å", *std::max_element(m_rmsd_vector.begin(), m_rmsd_vector.end())));

    if (has_energy_data) {
        CurcumaLogger::info("");
        CurcumaLogger::info("Energy Statistics:");
        CurcumaLogger::param("Mean energy", fmt::format("{:.6f} Eh", energy_mean));
        CurcumaLogger::param("Median energy", fmt::format("{:.6f} Eh", energy_median));
        CurcumaLogger::param("Std. deviation", fmt::format("{:.6f} Eh", energy_std));
        CurcumaLogger::param("Shannon entropy", fmt::format("{:.4f}", energy_shannon));
        CurcumaLogger::param("Min energy", fmt::format("{:.6f} Eh", *std::min_element(m_energy_vector.begin(), m_energy_vector.end())));
        CurcumaLogger::param("Max energy", fmt::format("{:.6f} Eh", *std::max_element(m_energy_vector.begin(), m_energy_vector.end())));
    }

    // Write to file for further analysis (preserve original functionality)
    if (m_rmsd_file.is_open()) {
        m_rmsd_file << "#Mean\t" << rmsd_mean << "\t" << energy_mean << std::endl;
        m_rmsd_file << "#Median\t" << rmsd_median << "\t" << energy_median << std::endl;
        m_rmsd_file << "#StdDev\t" << rmsd_std << "\t" << energy_std << std::endl;
        m_rmsd_file << "#Shannon\t" << rmsd_shannon << "\t" << energy_shannon << std::endl;
        CurcumaLogger::info_fmt("Statistical summary written to {}", m_outfile + "_rmsd.dat");
    }
}

void RMSDTraj::LoadControlJson()
{
    m_heavy = Json2KeyWord<bool>(m_defaults, "heavy");
    m_pcafile = Json2KeyWord<bool>(m_defaults, "pcafile");
    m_writeUnique = Json2KeyWord<bool>(m_defaults, "writeUnique");
    m_writeAligned = Json2KeyWord<bool>(m_defaults, "writeAligned");
    m_rmsd_threshold = Json2KeyWord<double>(m_defaults, "rmsd");
    m_fragment = Json2KeyWord<int>(m_defaults, "fragment");
    m_reference = Json2KeyWord<std::string>(m_defaults, "reference");
    m_second_file = Json2KeyWord<std::string>(m_defaults, "second");
    m_pairwise = (m_second_file.compare("none") != 0);
    m_allxyz = Json2KeyWord<bool>(m_defaults, "allxyz");
    m_ref_first = Json2KeyWord<bool>(m_defaults, "RefFirst");
    m_opt = Json2KeyWord<bool>(m_defaults, "opt");
    m_filter = Json2KeyWord<bool>(m_defaults, "filter");
    m_writeRMSD = Json2KeyWord<bool>(m_defaults, "writeRMSD");
    m_offset = m_defaults["offset"];
}

void RMSDTraj::Optimise()
{
    CurcumaOpt optimise(m_controller["rmsdtraj"], true);
    optimise.setFileName(m_outfile + ".unique.xyz");
    // optimise.setBaseName(m_outfile);
    optimise.start();
}

void RMSDTraj::Filter()
{
    ConfScan* scan = new ConfScan(m_controller["rmsdtraj"], false);
    scan->setFileName(m_outfile + ".opt.xyz");
    scan->start();
}
