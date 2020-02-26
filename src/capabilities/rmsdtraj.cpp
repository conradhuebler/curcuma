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

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "src/capabilities/rmsd.h"

#include "src/core/elements.h"
#include "src/core/molecule.h"

#include "src/tools/general.h"
#include "src/tools/geometry.h"

#include "rmsdtraj.h"

RMSDTraj::RMSDTraj()
{
}

void RMSDTraj::AnalyseTrajectory()
{

    int atoms_target = -1;
    if (m_reference.size() != 0) {
        m_stored_structures.push_back(Tools::LoadFile(m_reference));
        atoms_target = m_stored_structures[0].AtomCount();
    }

    string outfile = m_filename;
    for (int i = 0; i < 4; ++i)
        outfile.pop_back();

    m_rmsd_file.open(outfile + "_rmsd.dat");
    m_pca_file.open(outfile + "_pca.dat");

    RMSDDriver* driver = new RMSDDriver;
    driver->setSilent(true);
    driver->setProtons(true);
    driver->setForceReorder(false);
    driver->setCheckConnections(false);
    driver->setFragment(m_fragment);
    std::ifstream input(m_filename);
    std::vector<std::string> lines;
    int atoms = 0;
    int index = 0;
    int i = 0;
    bool xyzfile = std::string(m_filename).find(".xyz") != std::string::npos || std::string(m_filename).find(".trj") != std::string::npos;
    Molecule mol(atoms, 0);
    for (std::string line; getline(input, line);) {
        if (index == 0 && xyzfile) {
            atoms = stoi(line);
            mol = Molecule(atoms, 0);
        }
        if (xyzfile) {
            if (i > 1) {
                mol.setXYZ(line, i - 2);
            }
            if (i - 1 == atoms) {

                if (m_stored_structures.size() == 0)
                    m_stored_structures.push_back(mol);
                else {
                    for (std::size_t i = 0; i < mol.GetFragments().size(); ++i)
                        if (mol.getGeometryByFragment(i).rows() == atoms_target) {
                            driver->setFragmentTarget(i);
                            driver->setPartialRMSD(true);
                        }
                }
                driver->setScaling(1.3);

                std::vector<double> rmsd_results;
                for (std::size_t mols = 0; mols < m_stored_structures.size(); ++mols) {
                    driver->setReference(m_stored_structures[mols]);
                    driver->setTarget(mol);
                    driver->AutoPilot();
                    if (mols == 0) {
                        m_rmsd_file << driver->RMSD() << std::endl;

                        Molecule mol2 = driver->TargetAligned();

                        for (std::size_t j = 0; j < mol2.AtomCount(); ++j) {
                            if (mol2.Atom(j).first != 1)
                                m_pca_file << mol2.Atom(j).second(0) << " " << mol2.Atom(j).second(1) << " " << mol2.Atom(j).second(2);
                        }
                        m_pca_file << std::endl;
                    }
                    rmsd_results.push_back(driver->RMSD());
                }
                int add = 0;
                for (double rmsd : rmsd_results) {
                    add += rmsd > 1;
                }
                if (add == rmsd_results.size()) {
                    m_stored_structures.push_back(mol);
                    mol.appendXYZFile(outfile + "_unique.xyz");
                    std::cout << "New structure added ... ( " << m_stored_structures.size() << ")." << std::endl;
                }
                i = -1;
                mol = Molecule(atoms, 0);
            }
            ++i;
        } else {
            mol.setAtom(line, i);
        }
        index++;
    }

    delete driver;
}
