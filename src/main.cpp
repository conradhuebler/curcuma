/*
 * <Curcuma main file.>
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
//#include "src/capabilities/optimiser/Proton.h"

#include "src/tools/general.h"

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>

#if __GNUC__
// Thanks to
// https://stackoverflow.com/questions/77005/how-to-automatically-generate-a-stacktrace-when-my-program-crashes
#include <execinfo.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void bt_handler(int sig)
{
    void* array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    fprintf(stderr, "Curcuma crashed. Although this is probably unintended, it happened anyway.\n Some kind of backtrace will be printed out!\n\n");
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    fprintf(stderr, "Good-By\n");
    exit(1);
}
#endif

#include "json.hpp"
using json = nlohmann::json;

void Distance(const Molecule &mol, char **argv)
{
    int donor = stoi(std::string(argv[3]));
    int proton = stoi(std::string(argv[4]));
    int acceptor = stoi(std::string(argv[5]));
    std::cout << "Using atoms " << donor << " " << proton << " " << acceptor << std::endl;
    std::cout << "Donor ";
    mol.printAtom(donor);
    std::cout << std::endl
              << "Proton: ";
    mol.printAtom(proton);
    std::cout << std::endl
              << "Acceptor: ";
    mol.printAtom(acceptor);
    std::cout << std::endl
              << "Hydrogen Bond Angle: " << mol.CalculateAngle(donor - 1, proton - 1, acceptor - 1) << std::endl;
    std::cout << "Hydrogen bond length " << mol.CalculateDistance(proton - 1, acceptor - 1) << std::endl;
}

int main(int argc, char **argv) {
#if __GNUC__
    signal(SIGSEGV, bt_handler);
    signal(SIGABRT, bt_handler);
#endif

    General::StartUp(argc, argv);
    RunTimer timer(true);
    if(argc < 2)
    {
        std::cerr << "No arguments given!" << std::endl;
        std::cerr << "Use:" << std::endl
                  << "-rmsd        * RMSD Calculator               *" << std::endl
                  << "-confscan    * Filter list of conformers     *" << std::endl
                  << "-dock        * Perform some docking          *" << std::endl
                  << "-opt         * LBFGS optimiser using xtb GFN *" << std::endl
                  << "-rmsdtraj    * Find unique structures        *" << std::endl;
        exit(1);
    }
    if(argc >= 2)
    {
        json controller = CLI2Json(argc, argv);

        if(strcmp(argv[1], "-rmsd") == 0)
        {
            if (argc < 4) {
                std::cerr << "Please use curcuma for rmsd calcultion as follows\ncurcuma -rmsd A.xyz B.xyz" << std::endl;
                std::cerr << "Additonal arguments are:" << std::endl;
                std::cerr << "-reorder    **** Force reordering of structure!" << std::endl;
                std::cerr << "-check      **** Check methyl group connectivity." << std::endl;
                std::cerr << "-heavy      **** Calculate RMSD for heavy atoms only. Affects Reordering." << std::endl;
                std::cerr << "-fragment n **** Use n'th fragment. Bonds are determined from simple covalent radii for now!" << std::endl;
                std::cerr << "-init n     **** Initialse Reordering with fixed fragement n" << std::endl;
                exit(1);
            }

            Molecule mol1 = Tools::LoadFile(argv[2]); // will only take first structure
            Molecule mol2 = Tools::LoadFile(argv[3]); // will only take first structure
            if (mol1.AtomCount() == 0 || mol2.AtomCount() == 0)
                exit(0);
            RMSDDriver* driver = new RMSDDriver(controller, false);
            driver->setReference(mol1);
            driver->setTarget(mol2);
            driver->start();
            std::cout << "RMSD for two molecules " << driver->RMSD() << std::endl;

            driver->ReferenceAligned().writeXYZFile("reference.xyz");
            driver->TargetAligned().writeXYZFile("target_align.xyz");
            driver->TargetReorderd().writeXYZFile("target_reorder.xyz");

            std::cout << Tools::Vector2String(driver->ReorderRules()) << std::endl;
            delete driver;
            exit(0);
        } else if (strcmp(argv[1], "-dock") == 0) {
            if (argc < 4) {
                std::cerr << "Please use curcuma for docking as follows\ncurcuma -dock -host A.xyz -guest B.xyz -Step_x 10 -Step_y 10 -Step_z 10" << std::endl;
                exit(1);
            }

            Docking* docking = new Docking(controller, false);
            if (docking->Initialise() == false) {
                docking->printError();
                return 0;
            }
            docking->start();
        } else if (strcmp(argv[1], "-hbonds") == 0) {
            if(argc != 6)
            {
                std::cerr << "Please use curcuma for hydrogen bond analysis as follows\ncurcuma -hbonds A.xyz index_donor index_proton index_acceptor" << std::endl;
                return -1;
            }
            FileIterator file(argv[2]);
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
        } else if (strcmp(argv[1], "-confscan") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma for conformation scan and judge as follows\ncurcuma -confscan conffile.xyz" << std::endl;
                std::cerr << "Additonal arguments are:" << std::endl;
                std::cerr << "-writeXYZ  **** Write results to xyz files!" << std::endl;
                std::cerr << "-rank n    **** Write only the first n results!" << std::endl;
                std::cerr << "-reorder   **** Force reordering of structure! - It will be done automatically, if energies are close and rmsd is big." << std::endl;
                std::cerr << "-heavy     **** Use only heavy atoms for rmsd calculation." << std::endl;
                std::cerr << "-noname    **** Do not read possible name from xyz file." << std::endl;
                std::cerr << "-noreorder **** Prevent reordering in any cases." << std::endl;
                std::cerr << "-norestart **** Prevent restarting in any cases." << std::endl;
                std::cerr << "-maxenergy **** Maximal energy difference between best and current conformer [kJ/mol] for a conformer to be analysed." << std::endl;
                std::cerr << "-energy    **** Energy threshold for identical structures [kJ/mol]." << std::endl;
                return -1;
            }
            ConfScan* scan = new ConfScan(controller);
            scan->setFileName(argv[2]);
            scan->start();
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
            std::string outfile = std::string(argv[2]);
            for (int i = 0; i < 4; ++i)
                outfile.pop_back();
            outfile += "_opt.xyz";

            json key = OptJson;
            key = MergeJson(key, controller["opt"]);

            FileIterator file(argv[2]);
            std::multimap<double, Molecule> results;
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                Molecule mol2 = OptimiseGeometry(&mol, key);
                //mol2.writeXYZFile(outfile);
                results.insert(std::pair<double, Molecule>(mol2.Energy(), mol2));
            }
            for (const auto& ref : results)
                ref.second.appendXYZFile(outfile);

            return 0;
        } /* else if (strcmp(argv[1], "-optp") == 0) {
            if (argc < 2) {
                std::cerr << "Please use curcuma for optimisation as follows:\ncurcuma -opt input.xyz" << std::endl;
                return 0;
            }
            std::string outfile = std::string(argv[2]);
            for (int i = 0; i < 4; ++i)
                outfile.pop_back();
            outfile += "_opt.xyz";

            json key = OptJson;
            key = MergeJson(key, controller["opt"]);

            FileIterator file(argv[2]);
            std::multimap<double, Molecule> results;
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                Molecule mol2 = OptimiseProtons(&mol, key);
                //mol2.writeXYZFile(outfile);
                results.insert(std::pair<double, Molecule>(mol2.Energy(), mol2));
            }
            for (const auto& ref : results)
                ref.second.appendXYZFile(outfile);

            return 0;
        } */
        else if (strcmp(argv[1], "-md") == 0) {
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
                std::cerr << "-heavy        **** Check only heavy atoms. Do not use with -write." << std::endl;
                return 0;
            }
            /*
            int fragment = -1;
            std::string reference, second;
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
    */
            RMSDTraj traj(controller, false);
            //traj.setReferenceStructure(reference);
            //traj.WriteUnique(write_unique);
            //traj.setRMSDThreshold(rmsd);
            traj.setFile(argv[2]);
            //traj.setFragment(fragment);
            //traj.setHeavy(heavy);
            //if (isSecond)
            //    traj.setSecondFile(second);
            traj.start();

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
            result_file.open("centroids.dat");
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
        } else if (strcmp(argv[1], "-split") == 0) {
            if (argc < 2) {
                std::cerr << "Please use curcuma to split supramolecular structures as follows:\ncurcuma -split molecule.xyz" << std::endl;
                return 0;
            }
            FileIterator file(argv[2]);
            int index = 1;
            std::string outfile = argv[2];
            for (int i = 0; i < 4; ++i)
                outfile.pop_back();
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                mol.setScaling(1.2);
                std::cout << file.MaxMolecules() << std::endl;
                if (file.MaxMolecules() <= 1)
                    mol.writeXYZFragments(outfile);
                else
                    mol.writeXYZFragments(outfile + "_M" + std::to_string(index));
                index++;
            }

        } else if (strcmp(argv[1], "-distance") == 0) {
            if (argc < 4) {
                std::cerr << "Please use curcuma to calculate distances as follows:\ncurcuma -distance molecule.xyz indexA indexB" << std::endl;
                return 0;
            }

            int indexA = 0, indexB = 0;
            try {
                indexA = std::stoi(argv[3]);
            } catch (const std::invalid_argument& arg) {
                std::cerr << "Please use curcuma to calculate distances as follows:\ncurcuma -distance molecule.xyz indexA indexB" << std::endl;
                return 0;
            }
            try {
                indexB = std::stoi(argv[4]);
            } catch (const std::invalid_argument& arg) {
                std::cerr << "Please use curcuma to calculate distances as follows:\ncurcuma -distance molecule.xyz indexA indexB" << std::endl;
                return 0;
            }
            FileIterator file(argv[2]);
            std::string outfile = argv[2];

            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                std::cout << ":: " << mol.CalculateDistance(indexA - 1, indexB - 1) << "::" << std::endl;
            }

        } else if (strcmp(argv[1], "-angle") == 0) {
            if (argc < 6) {
                std::cerr << "Please use curcuma to calculate angles as follows:\ncurcuma -distance molecule.xyz indexA indexB indexC" << std::endl;
                return 0;
            }

            int indexA = 0, indexB = 0, indexC = 0;
            try {
                indexA = std::stoi(argv[3]);
            } catch (const std::invalid_argument& arg) {
                std::cerr << "Please use curcuma to calculate angles as follows:\ncurcuma -distance molecule.xyz indexA indexB indexC" << std::endl;
                return 0;
            }
            try {
                indexB = std::stoi(argv[4]);
            } catch (const std::invalid_argument& arg) {
                std::cerr << "Please use curcuma to calculate angles as follows:\ncurcuma -distance molecule.xyz indexA indexB indexC" << std::endl;
                return 0;
            }
            try {
                indexC = std::stoi(argv[5]);
            } catch (const std::invalid_argument& arg) {
                std::cerr << "Please use curcuma to calculate angles as follows:\ncurcuma -distance molecule.xyz indexA indexB indexC" << std::endl;
                return 0;
            }
            FileIterator file(argv[2]);

            printf("\n  Angle\t\tr(%u,%u)\tr(%u,%u)\tr(%u,%u)\n", indexA - 1, indexB - 1, indexA - 1, indexC - 1, indexC - 1, indexB - 1);

            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                printf(":: %8.4f\t%8.4f\t%8.4f\t%8.4f ::\n", mol.CalculateAngle(indexA - 1, indexB - 1, indexC - 1), mol.CalculateDistance(indexA - 1, indexB - 1), mol.CalculateDistance(indexA - 1, indexC - 1), mol.CalculateDistance(indexC - 1, indexB - 1));
            }
            printf("\n\n");
        } else {
            bool centered = false;
            for (std::size_t i = 2; i < argc; ++i) {
                if (strcmp(argv[i], "-center") == 0) {
                    centered = true;
                }
            }

            FileIterator file(argv[1]);
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                mol.setScaling(1.2);
                mol.CalculateRotationalConstants();
                if (centered)
                    mol.setGeometry(GeometryTools::TranslateGeometry(mol.getGeometry(), GeometryTools::Centroid(mol.getGeometry()), Position{ 0, 0, 0 }));
                mol.print_geom();
                mol.AnalyseIntermoleculeDistance();
                std::cout << std::endl
                          << std::endl;
                // CompactTopo(mol.HydrogenBondMatrix(-1,-1));
                // std::cout << CompareTopoMatrix(mol.HydrogenBondMatrix(-1,-1),mol.HydrogenBondMatrix(-1,-1) ) << std::endl;
            }
        }
    }
    return 0;
}
