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
    std::cout << "'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl;
    std::cout << "'    Scanning Trajectory file for RMSD and Conformers     '" << std::endl;
    std::cout << "'    Write Conformers ";
    if (m_writeUnique) {
        std::cout << "  Yes                              '" << std::endl;
        std::cout << "'    RMSD Threshold =  " << std::setprecision(3) << m_rmsd_threshold << "                                  '" << std::endl;
    } else
        std::cout << "  No                                '" << std::endl;
    std::cout << "'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl;
    if (m_reference.compare("none") != 0) {
        m_stored_structures.push_back(new Molecule(Files::LoadFile(m_reference)));
        m_atoms = m_stored_structures[0]->AtomCount();
    }

    m_outfile = Filename();
    for (int i = 0; i < 4; ++i)
        m_outfile.pop_back();

    if (m_writeRMSD)
        m_rmsd_file.open(m_outfile + "_rmsd.dat");

    if (m_pcafile)
        m_pca_file.open(m_outfile + "_pca.dat");

    if (m_pairwise) {
        m_pairwise_file.open(m_outfile + "_pairwise.dat");
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
        export_file.open(m_outfile + "_aligned.xyz");
        export_file.close();
    }
    std::ifstream input(Filename());
    std::vector<std::string> lines;
    //    int atoms = 0, atoms2 = 0;
    m_currentIndex = 0;
    //  int i = 0;
    //    int molecule = 0;
    std::ifstream inFile(Filename());
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
    if (m_second_file.compare("none") == 0)
        ProcessSingleFile();
    else
        CompareTrajectories();
}

void RMSDTraj::ProcessSingleFile()
{
    //  Molecule mol(m_atoms, 0);
    //  Molecule mol_2(m_atoms, 0);
    Molecule prev;
    FileIterator file(Filename());
    std::vector<int> progress(10, 0);
    while (!file.AtEnd()) {
        Molecule* molecule = new Molecule(file.Next());
        //   std::cout << molecule->Atom(0).second.transpose() << std::endl;
        bool check = CheckMolecule(molecule);
        if (check) {
            std::cout << "New structure added ... ( " << m_stored_structures.size() << "). " << /*  int(m_currentIndex / double(m_max_lines) * 100) << " % done ...!" << */ std::endl;
        } else {
        }
        if (progress[int((m_currentIndex / double(m_max_lines)) * 100) / 10] == 0) {
            progress[int((m_currentIndex / double(m_max_lines)) * 100) / 10] = 1;
            std::cout << int(m_currentIndex / double(m_max_lines) * 100) << " % done ...!" << std::endl;
        }
        delete molecule;
        if (CheckStop())
            break;
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
    m_currentIndex += molecule->AtomCount();

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

        if (m_driver->ReorderRules().size())
            std::cout << Tools::Vector2String(m_driver->ReorderRules()) << std::endl;

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
        std::cout << "Pairwise is disabled for now" << std::endl;
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
    FileIterator file1(Filename());
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
        { "method", "incr" },
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
        std::cout << rmsd << std::endl;
        delete mol1;
        delete mol2;
        if (CheckStop())
            break;
    }
    PostAnalyse();
}

void RMSDTraj::PostAnalyse()
{
    double rmsd_mean = Tools::mean(m_rmsd_vector);
    double rmsd_median = Tools::median(m_rmsd_vector);
    double rmsd_std = Tools::stdev(m_rmsd_vector, rmsd_mean);
    auto rmsd_hist = Tools::Histogram(m_rmsd_vector, 100);
    double rmsd_shannon = Tools::ShannonEntropy(rmsd_hist);

    double energy_mean = Tools::mean(m_energy_vector);
    double energy_median = Tools::median(m_energy_vector);
    double energy_std = Tools::stdev(m_energy_vector, energy_mean);
    auto energy_hist = Tools::Histogram(m_energy_vector, 100);
    double energy_shannon = Tools::ShannonEntropy(energy_hist);

    m_rmsd_file << "#" << rmsd_mean << "\t" << energy_mean << std::endl;
    m_rmsd_file << "#" << rmsd_median << "\t" << energy_median << std::endl;
    m_rmsd_file << "#" << rmsd_std << "\t" << energy_std << std::endl;
    m_rmsd_file << "#" << rmsd_shannon << "\t" << energy_shannon << std::endl;

    // delete driver;
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
