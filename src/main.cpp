/*
 * <Internal Coordinate Handler for chemical structures.>
 * Copyright (C) 2019 - 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/core/fileiterator.h"
#include "src/core/molecule.h"
#include "src/core/xtbinterface.h"

#include "src/capabilities/confscan.h"
#include "src/capabilities/docking.h"
#include "src/capabilities/nebdocking.h"
#include "src/capabilities/pairmapper.h"
#include "src/capabilities/rmsd.h"
#include "src/capabilities/rmsdtraj.h"
#include "src/capabilities/simplemd.h"

#include "src/capabilities/optimiser/LBFGSInterface.h"

#include "src/tools/general.h"

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>


void Distance(const Molecule &mol, char **argv)
{
    int donor = stoi(std::string(argv[3]));
    int proton = stoi(std::string(argv[4]));
    int acceptor = stoi(std::string(argv[5]));
    std::cout << "Using atoms" << std::endl;
    std::cout << "Donor ";
    mol.printAtom(donor);
    std::cout << endl << "Proton: ";
    mol.printAtom(proton);
    std::cout << endl << "Acceptor: ";
    mol.printAtom(acceptor);
    std::cout << std::endl << "Hydrogen Bond Angle: "<<  mol.angle(donor, proton, acceptor) << std::endl;
    std::cout << "Hydrogen bond length " << mol.Distance(proton - 1, acceptor - 1) << std::endl;

}



int main(int argc, char **argv) {

    General::StartUp(argc, argv);

    RunTimer timer(true);

    if(argc < 2)
    {
        std::cerr << "No arguments given!" << std::endl;
        std::cerr << "Use:" << std::endl
                  << "-rmsd        * RMSD Calulator                *" << std::endl
                  << "-confscan    * Filter list of conformers     *" << std::endl
                  << "-dock        * Perform some docking          *" << std::endl;
        XTBInterface interface;
    }
    
    if(argc >= 2)
    {

        if(strcmp(argv[1], "-rmsd") == 0)
        {
            if (argc < 4) {
                std::cerr << "Please use curcuma for rmsd calcultion as follows\ncurcuma -rmsd A.xyz B.xyz" << std::endl;
                std::cerr << "Additonal arguments are:" << std::endl;
                std::cerr << "-reorder    **** Force reordering of structure! - It will be done automatically, if energies are close and rmsd is big." << std::endl;
                std::cerr << "-check      **** Check methyl group connectivity." << std::endl;
                std::cerr << "-fragment n **** Use n'th fragment. Bonds are determined from simple covalent radii for now!" << std::endl;

                exit(1);
            }
            Molecule mol1 = Tools::LoadFile(argv[2]);
            Molecule mol2 = Tools::LoadFile(argv[3]);

            bool reorder = false, check_connect = false, heavy = false;
            int fragment = -1, pt = 0;
            if (argc >= 5) {

                for (std::size_t i = 4; i < argc; ++i) {
                    // std::cout << argv[i] << " " << strcmp(argv[i], "-fragment") << std::endl;
                    if (strcmp(argv[i], "-fragment") == 0) {
                        if (i + 1 < argc) {
                            fragment = std::stoi(argv[i + 1]);
                            ++i;
                        }
                        continue;
                    }

                    if (strcmp(argv[i], "-pt") == 0) {
                        if (i + 1 < argc) {
                            pt = std::stoi(argv[i + 1]);
                            ++i;
                        }
                        continue;
                    }

                    if (strcmp(argv[i], "-reorder") == 0) {
                        reorder = true;
                        continue;
                    }

                    if (strcmp(argv[i], "-heavy") == 0) {
                        heavy = true;
                        continue;
                    }

                    if (strcmp(argv[i], "-check") == 0) {
                        check_connect = true;
                        continue;
                    }
                }
        }

        // mol1.print_geom();
        // mol2.print_geom();

        RMSDDriver *driver = new RMSDDriver(mol1, mol2);
        driver->setForceReorder(reorder);
        driver->setProtons(!heavy);
        driver->setFragment(fragment);
        driver->setCheckConnections(check_connect);
        driver->AutoPilot();
        std::cout << "RMSD for two molecules " << driver->RMSD() << std::endl;

        driver->ReferenceAligned().writeXYZFile("reference.xyz");
        driver->TargetAligned().writeXYZFile("target_align.xyz");
        driver->TargetReorderd().writeXYZFile("target_reorder.xyz");

        delete driver;

        exit(0);
        } else if (strcmp(argv[1], "-dock") == 0) {

            if (argc < 4) {
                std::cerr << "Please use curcuma for docking  as follows\ncurcuma -dock A.xyz B.xyz XXX YYY ZZZ" << std::endl;
                exit(1);
            }
            bool opt = false, filter = false;
            bool check = false;

            if (argc >= 4) {

                for (std::size_t i = 4; i < argc; ++i) {

                    if (strcmp(argv[i], "-opt") == 0) {
                        opt = true;
                        continue;
                    }

                    if (strcmp(argv[i], "-filter") == 0) {
                        filter = true;
                        continue;
                    }

                    if (strcmp(argv[i], "-check") == 0) {
                        check = true;
                        continue;
                    }
                }
            }

            Molecule mol1 = Tools::LoadFile(argv[2]);
            Molecule mol2 = Tools::LoadFile(argv[3]);

            Docking* docking = new Docking;
            docking->setHostStructure(mol1);
            docking->setGuestStructure(mol2);
            /*
            if (argc >= 7) {
                double XXX = stod(std::string(argv[4]));
                double YYY = stod(std::string(argv[5]));
                double ZZZ = stod(std::string(argv[6]));

                std::cout << "Docking Position:" << XXX << " " << YYY << " " << ZZZ << std::endl;
                docking->setAnchorPosition(Position{ XXX, YYY, ZZZ });
            }
            if (argc >= 7) {


            }*/
            docking->setCheck(check);
            docking->PerformDocking();

            StringList files = docking->Files();
            StringList optfile;

            if (opt) {
                for (auto file : files) {

                    FileIterator fileiter(file);

                    std::multimap<double, Molecule> results;
                    while (!fileiter.AtEnd()) {
                        Molecule mol = fileiter.Next();
                        Molecule mol2 = OptimiseGeometry(&mol, false, true, 1, 0.05);
                        results.insert(std::pair<double, Molecule>(mol2.Energy(), mol2));
                    }
                    for (int i = 0; i < 4; ++i)
                        file.pop_back();

                    file += "_opt.xyz";
                    optfile.push_back(file);
                    for (const auto& ref : results)
                        ref.second.appendXYZFile(file);
                }

                if (filter) {
                    for (auto file : optfile) {

                        ConfScan* scan = new ConfScan;
                        scan->setFileName(file);
                        scan->setHeavyRMSD(true);
                        //scan->setMaxRank(rank);
                        //scan->setWriteXYZ(writeXYZ);
                        //scan->setForceReorder(reorder);
                        //scan->setCheckConnections(check_connect);
                        //scan->setEnergyThreshold(energy);
                        scan->setNoName(true);
                        scan->scan();
                    }
                }
            }
        } else if (strcmp(argv[1], "-hbonds") == 0) {
            if(argc != 6)
            {
                std::cerr << "Please use curcuma for hydrogen bond analysis as follows\ncurcuma -hbonds A.xyz index_donor index_proton index_acceptor" << std::endl;
                return -1;
            }

            FileIterator file(argv[1]);
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                if (argc == 6) {
                    if (std::string(argv[1]).find("-hbonds") != std::string::npos) {
                        Distance(mol, argv);
                    }
                } else {
                    mol.print_geom();
                    std::cout << std::endl
                              << std::endl;
                    std::cout << mol.getGeometry() << std::endl;
                }
            }

            /*
            std::cerr << "Opening file " << argv[2] << std::endl;
            std::ifstream input( argv[2] );
            std::vector<std::string> lines;
            int atoms = 0;
            int index = 0;
            int i = 0;
            bool xyzfile = std::string(argv[2]).find(".xyz") != std::string::npos;
            Molecule mol(atoms, 0);
            for( std::string line; getline( input, line ); )
            {
                if(index == 0 && xyzfile)
                {
                    atoms = stoi(line);
                    mol = Molecule(atoms, 0);
                }
                if(xyzfile)
                {
                    if(i > 1)
                    {
                        mol.setXYZ(line, i-2);
                    }
                    if(i-1 == atoms)
                    {
                        if(argc == 6)
                        {
                            if(std::string(argv[1]).find("-hbonds") != std::string::npos)
                            {
                               Distance(mol, argv);
                            }
                        }else
                        {
                            mol.print_geom();
                            std::cout << std::endl << std::endl;
                            std::cout << mol.getGeometry() << std::endl;
                        }
                        i = -1;
                        mol = Molecule(atoms, 0);
                    }
                    ++i;
                }else
                {
                    mol.setAtom(line, i);
                }
                index++;
            }
          */
        } else if (strcmp(argv[1], "-confscan") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma for conformation scan and judge as follows\ncurcuma -confscan conffile.xyz" << std::endl;
                std::cerr << "Additonal arguments are:" << std::endl;
                std::cerr << "-writeXYZ  **** Write results to xyz files!" << std::endl;
                std::cerr << "-rank n    **** Write only the first n results!" << std::endl;
                std::cerr << "-reorder   **** Force reordering of structure! - It will be done automatically, if energies are close and rmsd is big." << std::endl;
                std::cerr << "-heavy     **** Use only heavy atoms for rmsd calculation." << std::endl;
                std::cerr << "-noname    **** Do not read possible name from xyz file." << std::endl;

                return -1;
            }
            bool writeXYZ = false;
            bool reorder = false, check_connect = false, heavy = false, noname = false;
            int rank = 1e10;
            double energy = 1.0;
            if (argc >= 4) {

                for (std::size_t i = 3; i < argc; ++i) {

                    if (strcmp(argv[i], "-rank") == 0) {
                        if (i + 1 < argc) {
                            rank = std::stoi(argv[i + 1]);
                            ++i;
                        }
                    }

                    if (strcmp(argv[i], "-energy") == 0) {
                        if (i + 1 < argc) {
                            energy = std::stod(argv[i + 1]);
                            ++i;
                        }
                    }

                    if (strcmp(argv[i], "-writeXYZ") == 0) {
                        writeXYZ = true;
                        continue;
                    }

                    if (strcmp(argv[i], "-reorder") == 0) {
                        reorder = true;
                        continue;
                    }

                    if (strcmp(argv[i], "-check") == 0) {
                        check_connect = true;
                        continue;
                    }

                    if (strcmp(argv[i], "-heavy") == 0) {
                        heavy = true;
                        continue;
                    }

                    if (strcmp(argv[i], "-noname") == 0) {
                        noname = true;
                        continue;
                    }
                }
            }

            std::cerr << "Opening file " << argv[2] << std::endl;

            ConfScan* scan = new ConfScan;
            scan->setFileName(argv[2]);
            scan->setHeavyRMSD(heavy);
            scan->setMaxRank(rank);
            scan->setWriteXYZ(writeXYZ);
            scan->setForceReorder(reorder);
            scan->setCheckConnections(check_connect);
            scan->setEnergyThreshold(energy);
            scan->setNoName(noname);
            scan->scan();

            return 0;
        } else if (strcmp(argv[1], "-led") == 0) {
            if (argc < 2) {
                std::cerr << "Please use curcuma for fragment assignment as follows:\ncurcuma -led input.xyz" << std::endl;
                return 0;
            }

            Molecule mol1 = Tools::LoadFile(argv[2]);
            mol1.printFragmente();

        } else if (strcmp(argv[1], "-hmap") == 0) {
            if (argc < 2) {
                std::cerr << "Please use curcuma for hydrogen bond mapping as follows:\ncurcuma -hmap trajectory.xyz" << std::endl;
                return 0;
            }

            std::vector<std::pair<int, int>> pairs, elements;

            if (argc >= 3) {

                for (std::size_t i = 3; i < argc; ++i) {

                    if (strcmp(argv[i], "-pair") == 0) {
                        if (i + 2 < argc) {
                            if (Tools::isInt(argv[i + 1]) && Tools::isInt(argv[i + 2])) {
                                int first = std::stoi(argv[i + 1]) - 1;
                                int second = std::stoi(argv[i + 2]) - 1;
                                ++i;
                                pairs.push_back(std::pair<int, int>(first, second));
                            } else {
                                int first = Elements::String2Element(argv[i + 1]);
                                int second = Elements::String2Element(argv[i + 2]);
                                ++i;
                                elements.push_back(std::pair<int, int>(first, second));
                            }
                        }
                    }

                    if (strcmp(argv[i], "-pairfile") == 0) {
                        if (i + 1 < argc) {

                            std::ifstream input(argv[i + 1]);

                            for (std::string line; getline(input, line);) {

                                std::vector<std::string> numbers = Tools::SplitString(line);
                                if (numbers.size() == 2) {
                                    if (Tools::isInt(numbers[0]) && Tools::isInt(numbers[1])) {
                                        int first = std::stoi(numbers[0]) - 1;
                                        int second = std::stoi(numbers[1]) - 1;
                                        pairs.push_back(std::pair<int, int>(first, second));
                                    }
                                }
                            }
                        }
                    }
                }
            }

            PairMapper mapper;
            mapper.setFile(argv[2]);
            for (const std::pair<int, int>& pair : pairs)
                mapper.addPair(pair);
            for (const std::pair<int, int>& pair : elements)
                mapper.addElementPair(pair);

            mapper.FindPairs();

        } else if (strcmp(argv[1], "-opt") == 0) {
            if (argc < 2) {
                std::cerr << "Please use curcuma for optimisation as follows:\ncurcuma -opt input.xyz" << std::endl;
                return 0;
            }

            string outfile = std::string(argv[2]);
            for (int i = 0; i < 4; ++i)
                outfile.pop_back();
            outfile += "_opt.xyz";
            FileIterator file(argv[2]);
            std::multimap<double, Molecule> results;
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                Molecule mol2 = OptimiseGeometry(&mol, false, true);
                //mol2.writeXYZFile(outfile);
                results.insert(std::pair<double, Molecule>(mol2.Energy(), mol2));
            }
            for (const auto& ref : results)
                ref.second.appendXYZFile(outfile);

            return 0;
        } else if (strcmp(argv[1], "-md") == 0) {
            if (argc < 2) {
                std::cerr << "Please use curcuma for test md assignment as follows:\ncurcuma -md input.xyz" << std::endl;
                return 0;
            }

            Molecule mol1 = Tools::LoadFile(argv[2]);

            SimpleMD md;
            md.setMolecule(mol1);
            md.Initialise();
            md.Dance();
        } else if (strcmp(argv[1], "-rmsdtraj") == 0) {
            if (argc <= 2) {
                std::cerr << "Please use curcuma for rmsd analysis of trajectories as follows:\ncurcuma -rmsdtraj input.xyz" << std::endl;
                std::cerr << "Additonal arguments are:" << std::endl;
                std::cerr << "-write        **** Write unique conformers!" << std::endl;
                std::cerr << "-rmsd d       **** Set rmsd threshold to d ( default = 1.0)!" << std::endl;
                std::cerr << "-fragment n   **** Set fragment to n." << std::endl;
                std::cerr << "-reference    **** Add different xyz structure as reference." << std::endl;
                std::cerr << "-second       **** Add second trajectory." << std::endl;
                std::cerr << "-heavy        **** Check only heavy atoms. Do not use with write." << std::endl;

                return 0;
            }
            int fragment = -1;
            string reference, second;
            double rmsd = 1;
            bool write_unique = false, heavy = false;
            bool isSecond = false;
            for (std::size_t i = 3; i < argc; ++i) {
                if (strcmp(argv[i], "-fragment") == 0) {
                    if (i + 1 < argc) {
                        fragment = std::stoi(argv[i + 1]);
                        ++i;
                    }
                    // continue;
                }
                if (strcmp(argv[i], "-reference") == 0) {
                    if (i + 1 < argc) {
                        reference = argv[i + 1];
                        ++i;
                    }
                    // continue;
                }

                if (strcmp(argv[i], "-second") == 0) {
                    if (i + 1 < argc) {
                        second = argv[i + 1];
                        isSecond = true;
                        ++i;
                    }
                    // continue;
                }

                if (strcmp(argv[i], "-rmsd") == 0) {
                    if (i + 1 < argc) {
                        rmsd = std::stod(argv[i + 1]);
                        ++i;
                    }
                    // continue;
                }

                if (strcmp(argv[i], "-write") == 0) {
                    write_unique = true;
                    continue;
                }

                if (strcmp(argv[i], "-heavy") == 0) {
                    heavy = true;
                    continue;
                }
            }

            RMSDTraj traj;
            traj.setReferenceStructure(reference);
            traj.WriteUnique(write_unique);
            traj.setRMSDThreshold(rmsd);
            traj.setFile(argv[2]);
            traj.setFragment(fragment);
            traj.setHeavy(heavy);
            if (isSecond)
                traj.setSecondFile(second);
            traj.AnalyseTrajectory();

        } else if (strcmp(argv[1], "-nebprep") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma for geometry preparation for nudge-elastic-band calculation follows:\ncurcuma -nebprep first.xyz second.xyz" << std::endl;
                return 0;
            }

            int pt = 0;
            for (std::size_t i = 3; i < argc; ++i) {
                if (strcmp(argv[i], "-pt") == 0) {
                    if (i + 1 < argc) {
                        pt = std::stoi(argv[i + 1]);
                        ++i;
                    }
                }
            }
            Molecule mol1 = Tools::LoadFile(argv[2]);
            Molecule mol2 = Tools::LoadFile(argv[3]);

            NEBDocking* nebdock = new NEBDocking;
            nebdock->setStructures(mol1, mol2);
            nebdock->setProtonTransfer(pt);
            nebdock->Prepare();
            delete nebdock;

        } else if (strcmp(argv[1], "-centroid") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma for centroid calculation of user definable fragments:\ncurcuma -centroid first.xyz" << std::endl;
                return 0;
            }

            int pt = 0, fragment = 0;
            std::vector<int> frag;
            for (std::size_t i = 3; i < argc; ++i) {
                if (strcmp(argv[i], "-fragment") == 0) {
                    if (i + 1 < argc) {
                        fragment = std::stoi(argv[i + 1]);
                        ++i;
                    }
                    // continue;
                }
                if (strcmp(argv[i], "-addfragment") == 0) {
                    bool loop = true;
                    while (i + 1 < argc && loop) {
                        StringList list = Tools::SplitString(argv[i + 1]);
                        if (list.size() > 1) {
                            for (auto element : list)
                                frag.push_back(std::stoi(element) - 1);
                        }
                        if (Tools::isInt(argv[i + 1]))
                            frag.push_back(std::stoi(argv[i + 1]) - 1);
                        else
                            loop = false;
                        ++i;
                    }
                    // continue;
                }
            }

            if (frag.size()) {
                std::cout << "Using fragment of atoms :";
                for (int atom : frag)
                    std::cout << atom + 1 << " ";
                std::cout << std::endl;
                std::cout << "to calculate centroid!" << std::endl;
            }

            std::ofstream result_file;

            FileIterator file(argv[2]);
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                if (frag.size()) {
                    result_file << GeometryTools::Centroid(mol.getGeometry(frag)).transpose() << std::endl;
                    std::cout << mol.getGeometry(frag) << std::endl;
                } else {
                    mol.GetFragments(1.2);
                    result_file << GeometryTools::Centroid(mol.getGeometryByFragment(fragment)).transpose() << std::endl;
                }
            }

            /*
            result_file.open("centroids.dat");
            std::ifstream input(argv[2]);
            std::vector<std::string> lines;
            int atoms = 0;
            int index = 0;
            int i = 0;
            bool xyzfile = std::string(argv[2]).find(".xyz") != std::string::npos || std::string(argv[1]).find(".trj") != std::string::npos;
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

                        if (frag.size()) {
                            result_file << GeometryTools::Centroid(mol.getGeometry(frag)).transpose() << std::endl;
                        } else {
                            mol.GetFragments(1.2);
                            result_file << GeometryTools::Centroid(mol.getGeometryByFragment(fragment)).transpose() << std::endl;
                        }
                        i = -1;
                        mol = Molecule(atoms, 0);
                    }
                    ++i;
                } else {
                    mol.setAtom(line, i);
                }
                index++;
            }*/

        } else {

            FileIterator file(argv[1]);
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                mol.setScaling(1.2);
                mol.CalculateRotationalConstants();
                mol.print_geom();
                mol.AnalyseIntermoleculeDistance();
                std::cout << std::endl
                          << std::endl;
            }
        }
    }
    return 0;
}
