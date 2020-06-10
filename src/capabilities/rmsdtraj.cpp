/*
 * <Trajectory RMSD Analyse. >
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

#include "curcumamethod.h"

#include "src/capabilities/rmsd.h"

#include "src/core/elements.h"
#include "src/core/molecule.h"

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

void RMSDTraj::start()
{
    std::cout << "'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl;
    std::cout << "'    Scanning Trajectory file for RMSD and Conformers     '" << std::endl;
    std::cout << "'    Write Conformers ";
    if (m_write_unique) {
        std::cout << "  Yes                              '" << std::endl;
        std::cout << "'    RMSD Threshold =  " << std::setprecision(3) << m_rmsd_threshold << "                                  '" << std::endl;
    } else
        std::cout << "  No                                '" << std::endl;
    std::cout << "'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl;
    int atoms_target = -1;
    if (m_reference.compare("none") != 0) {
        m_stored_structures.push_back(Tools::LoadFile(m_reference));
        atoms_target = m_stored_structures[0].AtomCount();
    }

    std::string outfile = m_filename;
    for (int i = 0; i < 4; ++i)
        outfile.pop_back();

    m_rmsd_file.open(outfile + "_rmsd.dat");

    if (m_pcafile)
        m_pca_file.open(outfile + "_pca.dat");

    if (m_pairwise) {
        m_pairwise_file.open(outfile + "_pairwise.dat");
    }
    json rmsd = RMSDJson;
    RMSDDriver* driver = new RMSDDriver();
    driver->setSilent(true);
    driver->setProtons(!m_heavy);
    driver->setForceReorder(false);
    driver->setCheckConnections(false);
    driver->setFragment(m_fragment);
    std::ofstream export_file;
    if (m_write_unique) {
        export_file.open(outfile + "_unique.xyz");
        export_file.close();
    }
    std::ifstream input(m_filename);
    std::vector<std::string> lines;
    int atoms = 0, atoms2 = 0;
    int index = 0;
    int i = 0;
    int molecule = 0;
    std::ifstream inFile(m_filename);
    int max = std::count(std::istreambuf_iterator<char>(inFile),
        std::istreambuf_iterator<char>(), '\n');

    std::ifstream second;
    if (m_pairwise) {
        second.open(m_second_file);
    }

    bool xyzfile = std::string(m_filename).find(".xyz") != std::string::npos || std::string(m_filename).find(".trj") != std::string::npos;
    Molecule mol(atoms, 0);
    Molecule mol_2(atoms2, 0);
    for (std::string line; getline(input, line);) {
        std::string line2;
        if (m_pairwise) {
            getline(second, line2);
        }
        if (index == 0 && xyzfile) {
            atoms = stoi(line);
            mol = Molecule(atoms, 0);
            molecule++;
            mol.setName(std::to_string(molecule));

            if (m_pairwise) {
                atoms2 = stoi(line2);

                if (atoms2 != atoms) {
                    std::cout << "Different atom count in structures " << std::endl;
                    return;
                }
                mol_2 = Molecule(atoms2, 0);
                mol_2.setName(std::to_string(molecule));
            }
        }
        if (xyzfile) {
            if (i > 1) {
                mol.setXYZ(line, i - 2);
                if (m_pairwise)
                    mol_2.setXYZ(line2, i - 2);
            }
            if (i == 1) {
                StringList list = Tools::SplitString(line);
                if (list.size() == 4) {
                    if (list[0].compare("SCF") == 0 && list[1].compare("done") == 0) {
                        try {
                            mol.setEnergy(std::stod((list[2])));
                        } catch (const std::string& what_arg) {
                            mol.setEnergy(0);
                        }
                    } else {
                        if (list[3] == "")
                            mol.setEnergy(0);
                        else {
                            try {
                                mol.setEnergy(std::stod((list[3])));
                            } catch (const std::string& what_arg) {
                                mol.setEnergy(0);
                            }
                        }
                    }
                } else if (list.size() == 2) {
                    try {
                        mol.setEnergy(std::stod((list[0])));
                    } catch (const std::string& what_arg) {
                    }
                } else if (list.size() == 1) {
                    try {
                        mol.setEnergy(std::stod((list[0])));
                    } catch (const std::string& what_arg) {
                    }
                } else {
                    for (const std::string& s : list) {
                        double energy = 0;
                        if (Tools::isDouble(s)) {
                            energy = std::stod(s);
                            mol.setEnergy(energy);
                            break;
                        }
                    }
                }
            }
            if (i - 1 == atoms) {
                if (m_stored_structures.size() == 0) {
                    if (m_write_unique)
                        mol.appendXYZFile(outfile + "_unique.xyz");
                    m_stored_structures.push_back(mol);
                } else {
                    for (std::size_t i = 0; i < mol.GetFragments().size(); ++i)
                        if (mol.getGeometryByFragment(i).rows() == atoms_target) {
                            driver->setFragmentTarget(i);
                            driver->setPartialRMSD(true);
                        }
                }
                driver->setScaling(1.3);
                if (m_pairwise == false) {
                    std::vector<double> rmsd_results;
                    for (std::size_t mols = 0; mols < m_stored_structures.size(); ++mols) {
                        driver->setReference(m_stored_structures[mols]);
                        driver->setTarget(mol);
                        driver->start();
                        if (mols == 0) {
                            {
                                m_rmsd_file << driver->RMSD() << std::endl;
                                m_rmsd_vector.push_back(driver->RMSD());
                            }
                            Molecule mol2 = driver->TargetAligned();
                            if (m_pcafile) {
                                for (std::size_t j = 0; j < mol2.AtomCount(); ++j) {
                                    if (mol2.Atom(j).first != 1)
                                        m_pca_file << mol2.Atom(j).second(0) << " " << mol2.Atom(j).second(1) << " " << mol2.Atom(j).second(2);
                                }
                                m_pca_file << std::endl;
                            }
                        }
                        rmsd_results.push_back(driver->RMSD());
                    }
                    if (m_write_unique) {
                        int add = 0;
                        for (double rmsd : rmsd_results) {
                            add += rmsd > m_rmsd_threshold;
                        }
                        if (add == rmsd_results.size()) {
                            m_stored_structures.push_back(mol);
                            mol.appendXYZFile(outfile + "_unique.xyz");
                            std::cout << "New structure added ... ( " << m_stored_structures.size() << "). " << int(index / double(max) * 100) << " % done ...!" << std::endl;
                        }
                    }
                } else {
                    driver->setReference(mol);
                    driver->setTarget(mol_2);
                    driver->start();
                    m_pairwise_file << driver->RMSD() << std::endl;
                }
                i = -1;
                mol = Molecule(atoms, 0);
                molecule++;
                mol.setName(std::to_string(molecule));
                if (m_pairwise) {
                    mol_2 = Molecule(atoms, 0);
                    mol_2.setName(std::to_string(molecule));
                }
            }
            ++i;
        } else {
            mol.setAtom(line, i);
            if (m_pairwise)
                mol_2.setAtom(line2, i);
        }
        index++;
    }
    if (m_write_unique && m_allxyz)
        Tools::xyz2allxyz(outfile + "_unique.xyz");

    double mean = Tools::mean(m_rmsd_vector);
    double median = Tools::median(m_rmsd_vector);
    double std = Tools::stdev(m_rmsd_vector, mean);
    auto hist = Tools::Histogram(m_rmsd_vector, 100);
    double shannon = Tools::ShannonEntropy(hist);

    m_rmsd_file << "#" << mean << std::endl;
    m_rmsd_file << "#" << median << std::endl;
    m_rmsd_file << "#" << std << std::endl;
    m_rmsd_file << "#" << shannon << std::endl;

    delete driver;
}

void RMSDTraj::LoadControlJson()
{
    m_heavy = Json2KeyWord<bool>(m_defaults, "heavy");
    m_pcafile = Json2KeyWord<bool>(m_defaults, "pcafile");
    m_write_unique = Json2KeyWord<bool>(m_defaults, "write");
    m_rmsd_threshold = Json2KeyWord<double>(m_defaults, "rmsd");
    m_fragment = Json2KeyWord<int>(m_defaults, "fragment");
    m_reference = Json2KeyWord<std::string>(m_defaults, "reference");
    m_second_file = Json2KeyWord<std::string>(m_defaults, "second");
    m_pairwise = (m_second_file.compare("none") != 0);
    m_allxyz = Json2KeyWord<bool>(m_defaults, "allxyz");
}
