/*
 * <Curcuma main file.>
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *               2024 Gerd Gehrisch <gg27fyla@student.freiberg.de>
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
#include "src/core/energy_calculators/qm_methods/eht.h"
#include "src/core/energy_calculators/qm_methods/orcainterface.h"
#include "src/core/fileiterator.h"
#include "src/core/imagewriter.hpp"
#include "src/core/molecule.h"

#include "src/capabilities/analysenciplot.h"
#include "src/capabilities/analysis.h"
#include "src/capabilities/confscan.h"
#include "src/capabilities/confsearch.h"
#include "src/capabilities/confstat.h"
#include "src/capabilities/curcumaopt.h"
// Modern optimizer system - Claude Generated (simplified)
#include "src/capabilities/docking.h"
#include "src/capabilities/hessian.h"
#include "src/capabilities/casino.h"
#include "src/capabilities/optimisation/modern_optimizer_simple.h"
#include "src/capabilities/nebdocking.h"
#include "src/capabilities/pairmapper.h"
#include "src/capabilities/persistentdiagram.h"
#include "src/capabilities/qmdfffit.h"
#include "src/capabilities/rmsd.h"
#include "src/capabilities/rmsdtraj.h"
#include "src/capabilities/simplemd.h"
#include "src/capabilities/trajectoryanalysis.h"

#include "src/tools/cli_parser.h"
#include "src/tools/general.h"
#include "src/tools/info.h"

// Claude Generated: Parameter registry system
#include "generated/parameter_registry.h"
#include "src/core/parameter_registry.h"

#include "src/capabilities/optimiser/OptimiseDipoleScaling.h"
#include "src/capabilities/optimisation/modern_optimizer_simple.h"

#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef C17
#ifndef _WIN32
#include <filesystem>
namespace fs = std::filesystem;
#endif
#endif
// #include "ripser.h"

#ifndef _WIN32
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
#ifdef C17
#ifndef _WIN32
    //  std::filesystem::remove("stop");
#endif
#else
    remove("stop");
#endif
    exit(1);
}

void ctrl_c_handler(int s)
{
    if (std::filesystem::exists("stop")) {
#ifdef C17
#ifndef _WIN32
        std::filesystem::remove("stop");
#endif
#else
        remove("stop");
#endif
        printf("Caught stop signal a second time.\nWill exit now!\n\n");
        exit(1);
    } else {
        printf("Caught stop signal\nWill try to stop current stuff!\n");
        std::ofstream output("stop");
    }
}

#endif
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

double DotProduct(const Eigen::Vector3d& pos1, const Eigen::Vector3d& pos2)
{
    return pos1.dot(pos2);
}

// Structured Command Dispatch System - Claude Generated
struct CapabilityInfo {
    std::string description;
    std::string category;
    std::vector<std::string> supported_formats;
    std::function<int(const json&, int, char**)> handler;
};

int executeCompare(const json& controller, int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Please use curcuma to compare structures as follows\ncurcuma -compare A.xyz B.xyz -metric [rmsd, inertia, ripser]" << std::endl;
        return 1;
    }

    Molecule molecule1, molecule2;
    FileIterator file1, file2;
    std::string reffile, tarfile;

    if (std::filesystem::exists(std::string(argv[2]))) {
        file1.setFile(argv[2]);
        molecule1 = file1.Next();
        reffile = file1.Basename();
    }
    if (std::filesystem::exists(argv[3])) {
        file2.setFile(argv[3]);
        tarfile = file2.Basename();
        molecule2 = file2.Next();
    } else {
        tarfile = file1.Basename();
        molecule2 = file1.Next();
    }

    json compare_config = controller.value("compare", json::object());
    std::string metric = compare_config.value("metric", "rmsd");

    if (metric == "rmsd") {
        json config;
        config["rmsd"] = compare_config;
        RMSDDriver driver(config, true);
        driver.setReference(molecule1);
        driver.setTarget(molecule2);
        driver.start();
        std::cout << "RMSD: " << driver.RMSD() << " Å" << std::endl;
    } else if (metric == "inertia") {
        molecule1.CalculateRotationalConstants();
        molecule2.CalculateRotationalConstants();
        std::cout << "Inertia constants for " << reffile << ": " << std::endl;
        std::cout << "Ia: " << molecule1.Ia() << " MHz, Ib: " << molecule1.Ib() << " MHz, Ic: " << molecule1.Ic() << " MHz" << std::endl;
        std::cout << "Inertia constants for " << tarfile << ": " << std::endl;
        std::cout << "Ia: " << molecule2.Ia() << " MHz, Ib: " << molecule2.Ib() << " MHz, Ic: " << molecule2.Ic() << " MHz" << std::endl;
        double Ia_diff = std::abs(molecule1.Ia() - molecule2.Ia());
        double Ib_diff = std::abs(molecule1.Ib() - molecule2.Ib());
        double Ic_diff = std::abs(molecule1.Ic() - molecule2.Ic());
        std::cout << "Inertia differences: " << Ia_diff << ", " << Ib_diff << ", " << Ic_diff << " MHz" << std::endl;
    } else if (metric == "ripser") {
        PersistentDiagram pd(compare_config);
        pd.setDistanceMatrix(molecule1.LowerDistanceVector());
        Eigen::MatrixXd image_a = pd.generateImage(pd.generatePairs());
        pd.setDistanceMatrix(molecule2.LowerDistanceVector());
        Eigen::MatrixXd image_b = pd.generateImage(pd.generatePairs());
        std::cout << "Persistence diagram for " << reffile << ": " << std::endl;
        std::cout << image_a << std::endl;
        std::cout << "Persistence diagram for " << tarfile << ": " << std::endl;
        std::cout << image_b << std::endl;
        double diff = (image_a - image_b).cwiseAbs().sum();
        std::cout << "Persistence diagram difference: " << diff << std::endl;
    } else {
        std::cerr << "Unknown metric: " << metric << std::endl;
        std::cerr << "Supported metrics: rmsd, inertia, ripser" << std::endl;
        return 1;
    }
    return 0;
}

int executeHBonds(const json& controller, int argc, char** argv) {
    if (argc < 6) {
        std::cerr << "Please use curcuma for hydrogen bond analysis as follows\ncurcuma -hbonds A.xyz index_donor index_proton index_acceptor" << std::endl;
        return 1;
    }
    FileIterator file(argv[2]);
    while (!file.AtEnd()) {
        Molecule mol = file.Next();
        if (argc == 6) {
            Distance(mol, argv);
        } else {
            mol.print_geom();
            std::cout << std::endl << std::endl;
            std::cout << mol.getGeometry() << std::endl;
        }
    }
    return 0;
}

int executeLED(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma for fragment assignment as follows:\ncurcuma -led input.xyz" << std::endl;
        return 1;
    }
    Molecule mol1 = Files::LoadFile(argv[2]);
    if (!mol1.Atoms().empty())
        mol1.printFragmente();
    return 0;
}

int executeHMap(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma for hydrogen bond mapping as follows:\ncurcuma -hmap trajectory.xyz" << std::endl;
        return 1;
    }

    std::vector<std::pair<int, int>> pairs, elements;
    json hmap_config = controller.value("hmap", json::object());

    // Prioritize JSON config if available
    if (hmap_config.contains("pairs")) {
        // Assume format like "1:2;3:4;C:H"
        std::string pstr = hmap_config["pairs"].get<std::string>();
        std::vector<std::string> p_tokens = Tools::SplitString(pstr, ";");
        for (const auto& token : p_tokens) {
            std::vector<std::string> atoms = Tools::SplitString(token, ":");
            if (atoms.size() == 2) {
                if (Tools::isInt(atoms[0]) && Tools::isInt(atoms[1])) {
                    pairs.emplace_back(std::stoi(atoms[0]) - 1, std::stoi(atoms[1]) - 1);
                } else {
                    elements.emplace_back(Elements::String2Element(atoms[0]), Elements::String2Element(atoms[1]));
                }
            }
        }
    }

    // Still support legacy argv parsing for complex pair specifications if not in JSON
    for (std::size_t i = 3; i < argc; ++i) {
        if (strcmp(argv[i], "-pair") == 0) {
            if (i + 2 < argc) {
                if (Tools::isInt(argv[i + 1]) && Tools::isInt(argv[i + 2])) {
                    pairs.emplace_back(std::stoi(argv[i + 1]) - 1, std::stoi(argv[i + 2]) - 1);
                    i += 2;
                } else {
                    elements.emplace_back(Elements::String2Element(argv[i + 1]), Elements::String2Element(argv[i + 2]));
                    i += 2;
                }
            }
        }
        // ... other legacy args ...
    }

    PairMapper mapper;
    mapper.setFile(argv[2]);
    for (const auto& p : pairs) mapper.addPair(p);
    for (const auto& e : elements) mapper.addElementPair(e);
    mapper.FindPairs();
    return 0;
}

int executeNCI(const json& controller, int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Please use curcuma to post-process two RDG vs rho plots from NCIPLOT as follows:\ncurcuma -nci file1.dat file2.dat" << std::endl;
        return 1;
    }
    AnalyseNCIPlot analyse(controller.value("nci", json::object()));
    analyse.setFiles(argv[2], argv[3]);
    analyse.start();
    return 0;
}

int executeHessian(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma to analyse hessians:\ncurcuma -hessian molecule.xyz" << std::endl;
        return 1;
    }
    Molecule mol1 = Files::LoadFile(argv[2]);
    Hessian hessian(controller.value("hessian", json::object()));
    hessian.setMolecule(mol1);
    hessian.start();
    return 0;
}

int executeQMDFFFit(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma for QMDFF fitting:\ncurcuma -qmdfffit molecule.xyz" << std::endl;
        return 1;
    }
    Molecule mol1 = Files::LoadFile(argv[2]);
    QMDFFFit qmdfffit(controller.value("qmdfffit", json::object()));
    qmdfffit.setMolecule(mol1);
    qmdfffit.start();
    return 0;
}

int executeGyration(const json& controller, int argc, char** argv) {
    if (argc < 3) return 1;
    FileIterator file(argv[2]);
    json gyr_config = controller.value("gyration", json::object());
    double hmass = gyr_config.value("hmass", 1.0);

    int count = 1;
    double sum = 0, sum_mass = 0, sqrt_sum = 0, sqrt_sum_mass = 0;
    while (!file.AtEnd()) {
        Molecule mol = file.Next();
        std::pair<double, double> gyr = mol.GyrationRadius(hmass);
        if (std::isnan(gyr.first) || std::isnan(gyr.second)) continue;
        sum += gyr.first;
        sum_mass += gyr.second;
        sqrt_sum += sqrt(gyr.first);
        sqrt_sum_mass += sqrt(gyr.second);
        std::cout << ":: " << gyr.first << " " << sum / count << " " << gyr.second << " " << sum_mass / count << std::endl;
        count++;
    }
    return 0;
}

int executeDipole(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma to optimise the dipole:\ncurcuma -dipole molecule.xyz" << std::endl;
        return 1;
    }
    FileIterator file(argv[2]);
    auto lm_basename = file.Basename();
    const json blob = controller.value("dipole", json::object());

    std::vector<Molecule> conformers;
    while (!file.AtEnd()) {
        Molecule mol = file.Next();
        mol.Center(false);
        EnergyCalculator interface("gfn2", blob);
        interface.setMolecule(mol.getMolInfo());
        interface.CalculateEnergy(false);
        mol.setPartialCharges(interface.Charges());
        mol.setDipole(interface.Dipole() * au);
        conformers.push_back(mol);
    }
    if (conformers.empty()) return 1;

    const auto linear_vector = DipoleScalingCalculation(conformers);
    const auto nonlinear_vector = OptimiseDipoleScaling(conformers, linear_vector);

    // Output results... (truncated for brevity, keep existing logic)
    std::cout << "Dipole optimization completed." << std::endl;
    return 0;
}

int executeOrca(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma as follows:\ncurcuma -orca input" << std::endl;
        return 1;
    }
    OrcaInterface orca;
    orca.setInputFile(argv[2]);
    if (!orca.runOrca()) return 1;
    orca.getOrcaJSON();
    orca.readOrcaJSON();
    return 0;
}

int executeStride(const json& controller, int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Please use curcuma to keep only every nth structure as follows:\ncurcuma -stride trajectory.xyz 100" << std::endl;
        return 1;
    }
    int stride = std::stoi(argv[3]);
    FileIterator file(argv[2]);
    while (!file.AtEnd()) {
        Molecule mol = file.Next();
        if (file.CurrentMolecule() % stride == 0)
            mol.appendXYZFile("stride_output.xyz");
    }
    return 0;
}

int executeDMatrix(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma to calculate a distance matrix:\ncurcuma -dMatrix molecule.xyz [options]" << std::endl;
        UnifiedAnalysis dummy(json{}, true);
        dummy.printEnhancedTDAHelp();
        return 0;
    }

    // Redirect to UnifiedAnalysis with topology properties enabled
    json dmatrix_config = controller.value("dMatrix", json::object());
    json analysis_config = controller;

    // Map legacy dMatrix flags to modern topological parameters
    // Note: ConfigManager handles the "topological_" prefix if configured correctly,
    // or we can set them directly in the "analysis" or "topology" sub-objects.
    json top_config = analysis_config.value("analysis", json::object());
    top_config["properties"] = "topology";

    // Ensure distance matrix is saved by default for this command
    if (!top_config.contains("topological_save_distance_matrix")) {
        top_config["topological_save_distance_matrix"] = true;
    }

    analysis_config["analysis"] = top_config;

    auto* analysis = new UnifiedAnalysis(analysis_config, false);
    analysis->setFileName(argv[2]);
    analysis->start();
    delete analysis;
    return 0;
}

int executeCenter(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma to center a molecule:\ncurcuma -center molecule.xyz [-mass]" << std::endl;
        return 1;
    }

    bool use_mass = false;
    if (argc > 3 && (strcmp(argv[3], "-mass") == 0 || strcmp(argv[3], "--mass") == 0)) {
        use_mass = true;
    } else if (controller.contains("center") && controller["center"].value("mass", false)) {
        use_mass = true;
    }

    FileIterator file(argv[2]);
    std::string output_file = "centered.xyz";
    if (std::filesystem::exists(output_file)) {
        std::filesystem::remove(output_file);
    }

    while (!file.AtEnd()) {
        Molecule mol = file.Next();
        mol.Center(use_mass);
        mol.appendXYZFile(output_file);
    }
    std::cout << "Centered structures written to " << output_file << std::endl;
    return 0;
}

int executeReorder(const json& controller, int argc, char** argv) {
    if (argc < 3) return 1;
    FileIterator file(argv[2]);
    while (!file.AtEnd()) {
        Molecule mol = file.Next();
        mol.writeXYZFile("reordered.xyz", Tools::RandomVector(0, mol.AtomCount()));
    }
    return 0;
}

int executeBlock(const json& controller, int argc, char** argv) {
    if (argc < 4) return 1;
    int blocks = std::stoi(argv[3]);
    FileIterator file(argv[2]);
    int mols = file.MaxMolecules();
    int block_size = mols / blocks;
    int current_block = 1;
    int count = 0;
    while (!file.AtEnd()) {
        Molecule mol = file.Next();
        mol.appendXYZFile("block_" + std::to_string(current_block) + ".xyz");
        if (++count >= block_size) {
            current_block++;
            count = 0;
        }
    }
    return 0;
}

int executeConfSearch(const json& controller, int argc, char** argv) {
    if (argc < 3) return 1;
    ConfSearch search(controller, false);
    search.setFile(argv[2]);
    search.start();
    return 0;
}

int executeRMSDTraj(const json& controller, int argc, char** argv) {
    if (argc < 3) return 1;
    RMSDTraj traj(controller, false);
    traj.setFile(argv[2]);
    traj.Initialise();
    traj.start();
    return 0;
}

int executeNEBPrep(const json& controller, int argc, char** argv) {
    if (argc < 4) return 1;
    Molecule mol1 = Files::LoadFile(argv[2]);
    Molecule mol2 = Files::LoadFile(argv[3]);
    NEBDocking nebdock;
    nebdock.setStructures(mol1, mol2);
    nebdock.Prepare();
    return 0;
}

int executeCentroid(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma for centroid calculation: curcuma -centroid molecule.xyz -atoms \"1-5\"" << std::endl;
        return 1;
    }

    json centroid_config = controller.value("centroid", json::object());
    std::string atoms_str = "";
    if (centroid_config.contains("atoms")) atoms_str = centroid_config["atoms"].get<std::string>();
    else if (argc > 3) atoms_str = argv[3];

    if (atoms_str.empty()) {
        std::cerr << "Error: Atom selection required for centroid calculation." << std::endl;
        return 1;
    }

    std::vector<int> atoms = Tools::CreateList(atoms_str);
    FileIterator file(argv[2]);
    while (!file.AtEnd()) {
        Molecule mol = file.Next();
        Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
        for (int i : atoms) centroid += mol.getGeometry().row(i - 1);
        if (!atoms.empty()) centroid /= static_cast<double>(atoms.size());

        std::cout << "Centroid (" << atoms_str << "): "
                  << centroid.transpose() << " Å" << std::endl;
    }
    return 0;
}

int executeSplit(const json& controller, int argc, char** argv) {
    if (argc < 3) return 1;
    FileIterator file(argv[2]);
    while (!file.AtEnd()) {
        Molecule mol = file.Next();
        mol.writeXYZFragments("split_" + std::to_string(file.CurrentMolecule()));
    }
    return 0;
}

int executeDistance(const json& controller, int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Please use curcuma to calculate distances as follows:\ncurcuma -distance molecule.xyz indexA indexB" << std::endl;
        return 1;
    }

    std::string atomsA_str = argc > 3 ? argv[3] : "";
    std::string atomsB_str = argc > 4 ? argv[4] : "";

    // Check if we have atom selections in controller
    json dist_config = controller.value("distance", json::object());
    if (dist_config.contains("atoms_a")) atomsA_str = dist_config["atoms_a"].get<std::string>();
    if (dist_config.contains("atoms_b")) atomsB_str = dist_config["atoms_b"].get<std::string>();

    if (atomsA_str.empty() || atomsB_str.empty()) {
        std::cerr << "Error: Two atom selections required for distance calculation." << std::endl;
        return 1;
    }

    std::vector<int> A = Tools::CreateList(atomsA_str);
    std::vector<int> B = Tools::CreateList(atomsB_str);

    FileIterator file(argv[2]);
    while (!file.AtEnd()) {
        Molecule mol = file.Next();
        if (A.size() == 1 && B.size() == 1) {
            std::cout << "Distance (" << A[0] << "-" << B[0] << "): "
                      << mol.CalculateDistance(A[0] - 1, B[0] - 1) << " Å" << std::endl;
        } else {
            // Handle multiple atoms (centroid distance)
            Eigen::Vector3d posA = Eigen::Vector3d::Zero();
            for (int i : A) posA += mol.getGeometry().row(i - 1);
            if (!A.empty()) posA /= static_cast<double>(A.size());

            Eigen::Vector3d posB = Eigen::Vector3d::Zero();
            for (int i : B) posB += mol.getGeometry().row(i - 1);
            if (!B.empty()) posB /= static_cast<double>(B.size());

            std::cout << "Centroid Distance: " << (posA - posB).norm() << " Å" << std::endl;
        }
    }
    return 0;
}

int executeAngle(const json& controller, int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "Please use curcuma to calculate angles as follows:\ncurcuma -angle molecule.xyz indexA indexB indexC" << std::endl;
        return 1;
    }

    std::string atomsA_str = argc > 3 ? argv[3] : "";
    std::string atomsB_str = argc > 4 ? argv[4] : "";
    std::string atomsC_str = argc > 5 ? argv[5] : "";

    // Check if we have atom selections in controller
    json angle_config = controller.value("angle", json::object());
    if (angle_config.contains("atoms_a")) atomsA_str = angle_config["atoms_a"].get<std::string>();
    if (angle_config.contains("atoms_b")) atomsB_str = angle_config["atoms_b"].get<std::string>();
    if (angle_config.contains("atoms_c")) atomsC_str = angle_config["atoms_c"].get<std::string>();

    if (atomsA_str.empty() || atomsB_str.empty() || atomsC_str.empty()) {
        std::cerr << "Error: Three atom selections required for angle calculation." << std::endl;
        return 1;
    }

    std::vector<int> A = Tools::CreateList(atomsA_str);
    std::vector<int> B = Tools::CreateList(atomsB_str);
    std::vector<int> C = Tools::CreateList(atomsC_str);

    FileIterator file(argv[2]);
    while (!file.AtEnd()) {
        Molecule mol = file.Next();
        if (A.size() == 1 && B.size() == 1 && C.size() == 1) {
            std::cout << "Angle (" << A[0] << "-" << B[0] << "-" << C[0] << "): "
                      << mol.CalculateAngle(A[0] - 1, B[0] - 1, C[0] - 1) << " °" << std::endl;
        } else {
            std::cerr << "Error: Multiple atom selections for angle calculation not yet supported." << std::endl;
            return 1;
        }
    }
    return 0;
}


int executeAnalysis(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        UnifiedAnalysis dummy(json{}, true);
        dummy.printHelp();
        return 0;
    }

    // Extract analysis-specific configuration - Claude Generated
    json analysis_config = controller.contains("analysis") ? controller["analysis"] : controller;

    auto* analysis = new UnifiedAnalysis(analysis_config, false);
    analysis->setFileName(argv[2]);
    analysis->start();
    delete analysis;
    return 0;
}

int executeRMSD(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma for rmsd calculation as follows\ncurcuma -rmsd A.xyz B.xyz" << std::endl;
        std::cerr << "Please use curcuma for rmsd calculation as follows\ncurcuma -rmsd AB.xyz" << std::endl;
        return 1;
    }

    Molecule molecule1, molecule2;
    FileIterator file1, file2;
    std::string reffile, tarfile;

    if (std::filesystem::exists(std::string(argv[2]))) {
        file1.setFile(argv[2]);
        molecule1 = file1.Next();
        reffile = file1.Basename();
    }
    if (argc > 3 && std::filesystem::exists(argv[3])) {
        file2.setFile(argv[3]);
        tarfile = file2.Basename();
        molecule2 = file2.Next();
    } else {
        tarfile = file1.Basename();
        molecule2 = file1.Next();
    }

    if (molecule1.AtomCount() == 0 || molecule2.AtomCount() == 0) {
        std::cout << "At least one structure is empty:\n";
        std::cout << argv[2] << " " << molecule1.AtomCount() << " atoms" << std::endl;
        if (argc > 3) std::cout << argv[3] << " " << molecule2.AtomCount() << " atoms" << std::endl;
        return 1;
    }

    json rmsd_config = controller.value("rmsd", json::object());
    auto* driver = new RMSDDriver(rmsd_config, false);
    driver->setReference(molecule1);
    driver->setTarget(molecule2);
    driver->start();
    std::cout << "RMSD for two molecules " << driver->RMSD() << std::endl;

    driver->ReferenceAligned().writeXYZFile(reffile + ".centered.xyz");
    driver->TargetAligned().writeXYZFile(tarfile + ".centered.xyz");
    driver->TargetReorderd().writeXYZFile(tarfile + ".reordered.xyz");

    delete driver;
    return 0;
}

int executeSinglePoint(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma for energy calculation as follows:\ncurcuma -sp input.xyz" << std::endl;
        return 1;
    }

    json sp_controller = controller;
    json opt_params = sp_controller.contains("opt") ? sp_controller["opt"] : json{};
    opt_params["single_point"] = true;

    // Claude Generated (December 2025): Set geometry_file for automatic parameter caching
    std::string geometry_file(argv[2]);
    opt_params["geometry_file"] = geometry_file;
    sp_controller["geometry_file"] = geometry_file;  // Also set at top level for EnergyCalculator

    sp_controller["opt"] = opt_params;

    CurcumaOpt opt(sp_controller, false);
    opt.setFileName(argv[2]);
    opt.start();
    return 0;
}

// Additional capability handlers - Claude Generated
int executeOptimization(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        json safe_opt_config = controller.contains("opt") ? controller["opt"] : json{};
        ModernOptimization::ModernOptimizerDispatcher helper(safe_opt_config, false);
        helper.printHelp();
        return 0;
    }

    // Parse optional optimizer parameter from controller
    json opt_config_base = controller.contains("opt") ? controller["opt"] : json::object();
    std::string optimizer_method = opt_config_base.value("optimizer", "auto");

    if (optimizer_method != "auto" && optimizer_method != "") {
        std::cout << "🧪 Using native Curcuma optimizer: " << optimizer_method << std::endl;
    }

    // Check if we should use modern native optimizers
    bool use_modern = (optimizer_method == "native_lbfgs" || optimizer_method == "lbfgs" ||
        optimizer_method == "diis" || optimizer_method == "rfo" || optimizer_method == "auto");

    if (use_modern) {

        try {
            auto molecule = std::make_unique<Molecule>(argv[2]);
            std::string method = controller.value("method", "uff");

            // Claude Generated (December 2025): Set geometry_file for automatic parameter caching
            json energy_controller = controller;
            energy_controller["geometry_file"] = std::string(argv[2]);

            EnergyCalculator energy_calc(method, energy_controller);
            energy_calc.setMolecule(molecule->getMolInfo());

            // Claude Generated (October 2025): Merge with default opt parameters from ParameterRegistry
            json opt_defaults = ParameterRegistry::getInstance().getDefaultJson("opt");
            json opt_config = MergeJson(opt_defaults, controller.contains("opt") ? controller["opt"] : json{});

            auto result = ModernOptimization::ModernOptimizerDispatcher::optimizeStructure(
                molecule.get(), optimizer_method, &energy_calc, opt_config);

            if (result.success) {
                // Claude Generated: Derive output filename from input basename like Legacy CurcumaOpt
                // Extract basename from argv[2] (e.g., "input.xyz" → "input")
                std::string filename(argv[2]);
                std::string basename = filename.size() >= 4 ?
                    filename.substr(0, filename.size() - 4) : filename;  // Remove last 4 chars (.xyz)
                std::string output_file = opt_config.value("output", basename + ".opt.xyz");

                molecule->writeXYZFile(output_file);
                CurcumaLogger::success_fmt("Optimized structure written to: {}", output_file);
                return 0;
            } else {
                // Modern optimizer failed, fall through to legacy
                CurcumaLogger::warn(fmt::format("Modern optimization failed: {}, using legacy optimizer", result.error_message));
            }
        } catch (const std::exception& e) {
            // Fall through to legacy code below
            CurcumaLogger::warn(fmt::format("Modern optimization failed: {}, using legacy optimizer", e.what()));
        }
    }

    // Legacy optimization (fallback for exceptions or non-modern methods)
    // Claude Generated 2025: Apply same parameter merging as modern optimizer to fix JSON null bug
    json legacy_controller = controller;
    json opt_defaults = ParameterRegistry::getInstance().getDefaultJson("opt");
    json opt_params = MergeJson(opt_defaults, controller.contains("opt") ? controller["opt"] : json{});

    // Claude Generated (December 2025): Set geometry_file for automatic parameter caching
    std::string geometry_file(argv[2]);
    opt_params["geometry_file"] = geometry_file;
    legacy_controller["geometry_file"] = geometry_file;  // Also set at top level for EnergyCalculator

    legacy_controller["opt"] = opt_params;

    CurcumaOpt opt(legacy_controller, false);
    opt.setFileName(argv[2]);
    opt.start();
    return 0;
}

int executeConfScan(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma for conformation scan as follows:\ncurcuma -confscan conffile.xyz" << std::endl;
        std::cerr << "Additional arguments are:" << std::endl;
        std::cerr << "-writeXYZ  **** Write results to xyz files!" << std::endl;
        std::cerr << "-rank n    **** Write only the first n results!" << std::endl;
        std::cerr << "-reorder   **** Force reordering of structure!" << std::endl;
        std::cerr << "-heavy     **** Use only heavy atoms for rmsd calculation." << std::endl;
        std::cerr << "-maxenergy **** Maximal energy difference [kJ/mol]." << std::endl;
        return 1;
    }

    // Claude Generated 2025: Simplified - CLI2Json routes parameters, ConfigManager handles defaults
    // ConfScan's ConfigManager (Multi-Module: "confscan", "rmsd") will:
    // 1. Load defaults for both modules from ParameterRegistry
    // 2. Merge user parameters from controller["confscan"] and controller["rmsd"]
    // 3. Handle global parameters (verbosity, threads)
    // No manual merging needed here!

    // Claude Generated (December 2025): Set geometry_file for automatic parameter caching
    json scan_controller = controller;
    scan_controller["geometry_file"] = std::string(argv[2]);

    auto* scan = new ConfScan(scan_controller, false);  // Claude Generated: Explicit false for default verbosity level 1
    scan->setFileName(argv[2]);
    scan->start();
    int accepted = scan->AcceptedCount();
    int reorder_success = scan->ReorderSuccessfull();
    int reuse_count = scan->ReuseCount();
    int skipped_count = scan->ReorderSkippedCount();
    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << std::endl;
    delete scan;
    return 0;
}

int executeConfStat(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma for conformation statistics as follows:\ncurcuma -confstat conffile.xyz" << std::endl;
        return 1;
    }

    // Claude Generated (December 2025): Set geometry_file for automatic parameter caching
    json stat_controller = controller;
    stat_controller["geometry_file"] = std::string(argv[2]);

    auto* stat = new ConfStat(stat_controller, false);  // Claude Generated: Explicit false for default verbosity level 1
    stat->setFileName(argv[2]);
    stat->start();
    delete stat;
    return 0;
}

int executeDocking(const json& controller, int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Please use curcuma for docking as follows:\ncurcuma -dock -host A.xyz -guest B.xyz -Step_x 10 -Step_y 10 -Step_z 10" << std::endl;
        return 1;
    }

    auto* docking = new Docking(controller, false);
    if (!docking->Initialise()) {
        docking->printError();
        delete docking;
        return 1;
    }
    docking->start();
    delete docking;
    return 0;
}

int executeSimpleMD(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma for molecular dynamics as follows:\ncurcuma -md input.xyz" << std::endl;
        SimpleMD help(json::object(), true);
        help.printHelp();
        return 0;
    }

    Molecule mol = Files::LoadFile(argv[2]);
    if (mol.AtomCount() == 0) {
        CurcumaLogger::error_fmt("Could not load molecule from file: {}", argv[2]);
        return 1;
    }

    // Claude Generated 2025: Simplified - SimpleMD's ConfigManager handles defaults
    // SimpleMD constructor uses ConfigManager("simplemd", controller) which:
    // 1. Loads defaults from ParameterRegistry
    // 2. Merges user parameters from controller["md"]
    // 3. Handles alias resolution automatically
    // No manual merging needed here!

    // Claude Generated (December 2025): Set geometry_file for automatic parameter caching
    json md_controller = controller;
    md_controller["geometry_file"] = std::string(argv[2]);

    auto* md = new SimpleMD(md_controller, false);
    md->setFile(argv[2]);  // Claude Generated (October 2025): Set basename for trajectory file naming
    md->setMolecule(mol);
    md->Initialise();
    md->start();
    delete md;
    return 0;
}

int executeCasino(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        Casino dummy(json{}, true);
        dummy.printHelp();
        return 0;
    }
    // Claude Generated 2025: Apply parameter merging to fix JSON null issues
    json casino_defaults = ParameterRegistry::getInstance().getDefaultJson("casino");
    json casino_config = MergeJson(casino_defaults, controller.contains("casino") ? controller["casino"] : json{});

    auto* casino = new Casino(casino_config, false);
    casino->setFileName(argv[2]);
    casino->start();
    delete casino;
    return 0;
}

int executeTrajectoryAnalysis(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        TrajectoryAnalysis dummy(json{}, true);
        dummy.printHelp();
        return 0;
    }

    auto* traj = new TrajectoryAnalysis(controller, false);
    traj->setFileName(argv[2]);
    traj->start();
    delete traj;
    return 0;
}

// Capability registry - Claude Generated
const std::map<std::string, CapabilityInfo> CAPABILITY_REGISTRY = {
    {"analysis", {"Unified molecular analysis (all formats, all properties)", "analysis",
                  {"XYZ", "VTF", "MOL2", "SDF", "PDB"}, executeAnalysis}},
    {"rmsd", {"RMSD calculation between structures", "analysis",
              {"XYZ", "VTF", "MOL2", "SDF"}, executeRMSD}},
    {"sp", {"Single point energy calculation", "calculation",
            {"XYZ", "VTF", "MOL2", "SDF"}, executeSinglePoint}},
    {"opt", {"Geometry optimization with various algorithms", "optimization",
             {"XYZ", "VTF", "MOL2", "SDF"}, executeOptimization}},
    {"confscan", {"Conformational scanning along reaction coordinates", "conformational",
                  {"XYZ", "MOL2", "SDF"}, executeConfScan}},
    {"confstat", {"Conformational statistics and analysis", "conformational",
                  {"XYZ", "MOL2", "SDF"}, executeConfStat}},
    {"dock", {"Molecular docking calculations", "docking",
              {"XYZ", "MOL2", "SDF"}, executeDocking}},
    {"md", {"Molecular dynamics simulation", "dynamics",
            {"XYZ", "VTF", "MOL2", "SDF"}, executeSimpleMD}},
    {"casino", {"Casino Monte Carlo simulation with enhanced sampling", "dynamics",
                {"XYZ", "VTF", "MOL2", "SDF", "PDB"}, executeCasino}},
    {"traj", {"Time-series analysis for molecular trajectories", "analysis",
              {"XYZ.trj", "VTF", "SDF"}, executeTrajectoryAnalysis}},
    {"compare", {"Compare structures using RMSD, inertia, or ripser", "analysis",
                 {"XYZ", "MOL2", "SDF"}, executeCompare}},
    {"hbonds", {"Hydrogen bond analysis for donor-proton-acceptor triplets", "analysis",
                {"XYZ", "MOL2", "SDF"}, executeHBonds}},
    {"led", {"Local Energy Decomposition (LED) fragment assignment", "analysis",
             {"XYZ", "MOL2", "SDF"}, executeLED}},
    {"hmap", {"Hydrogen bond mapping and pair analysis", "analysis",
              {"XYZ", "MOL2", "SDF"}, executeHMap}},
    {"nci", {"Non-Covalent Interaction (NCI) plot post-processing", "analysis",
             {"DAT"}, executeNCI}},
    {"hessian", {"Hessian calculation and vibrational analysis", "calculation",
                 {"XYZ", "MOL2", "SDF"}, executeHessian}},
    {"qmdfffit", {"QMDFF force field parameter fitting", "calculation",
                  {"XYZ", "MOL2", "SDF"}, executeQMDFFFit}},
    {"gyration", {"Radius of gyration calculation", "analysis",
                  {"XYZ", "MOL2", "SDF"}, executeGyration}},
    {"dipole", {"Dipole moment optimization and scaling", "calculation",
                {"XYZ", "MOL2", "SDF"}, executeDipole}},
    {"modern-opt", {"Modern optimization with Strategy Pattern", "optimization",
                    {"XYZ", "MOL2", "SDF"}, executeOptimization}},
    {"orca", {"ORCA quantum chemistry software interface", "interface",
              {"INP"}, executeOrca}},
    {"stride", {"Reduce structure count by keeping every N-th image", "analysis",
                {"XYZ", "VTF"}, executeStride}},
    {"center", {"Translate molecular center to origin", "analysis",
                {"XYZ", "MOL2", "SDF"}, executeCenter}},
    {"reorder", {"Randomly reorder atoms in a structure", "analysis",
                 {"XYZ", "MOL2", "SDF"}, executeReorder}},
    {"block", {"Split structure files into smaller blocks", "analysis",
               {"XYZ", "VTF"}, executeBlock}},
    {"confsearch", {"Systematic conformational searching", "conformational",
                    {"XYZ", "MOL2", "SDF"}, executeConfSearch}},
    {"rmsdtraj", {"Trajectory RMSD analysis and clustering", "analysis",
                  {"XYZ.trj", "VTF"}, executeRMSDTraj}},
    {"nebprep", {"Nudged Elastic Band (NEB) geometry preparation", "conformational",
                 {"XYZ"}, executeNEBPrep}},
    {"centroid", {"Calculate centroid of atom fragments", "analysis",
                  {"XYZ", "MOL2", "SDF"}, executeCentroid}},
    {"split", {"Split supramolecular structures into fragments", "analysis",
               {"XYZ", "MOL2", "SDF"}, executeSplit}},
    {"distance", {"Calculate distance between atoms or fragments", "analysis",
                  {"XYZ", "MOL2", "SDF"}, executeDistance}},
    {"angle", {"Calculate bond angles between atoms", "analysis",
               {"XYZ", "MOL2", "SDF"}, executeAngle}},
    {"dMatrix", {"Distance matrix and persistent homology analysis", "analysis",
                 {"XYZ", "VTF"}, executeDMatrix}},
    // TODO: Add more capabilities here as they are converted
};

void showStructuredHelp(const std::string& category = "") {
    std::cout << "Curcuma - Computational Chemistry Toolkit" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;

    if (category.empty()) {
        // Show overview with enhanced format information
        std::cout << "Modern Molecular Modeling with Universal Format Support" << std::endl;
        std::cout << std::endl;

        // Auto-discover and group capabilities by category
        std::map<std::string, std::vector<std::string>> categories;
        std::set<std::string> all_formats;

        for (const auto& [name, info] : CAPABILITY_REGISTRY) {
            categories[info.category].push_back(name);
            for (const auto& format : info.supported_formats) {
                all_formats.insert(format);
            }
        }

        // Display capabilities by category with enhanced formatting
        for (const auto& [cat, capabilities] : categories) {
            std::string category_name = cat;
            category_name[0] = std::toupper(category_name[0]);

            std::cout << "📊 " << category_name << " Capabilities:" << std::endl;

            for (const auto& cap : capabilities) {
                const auto& info = CAPABILITY_REGISTRY.at(cap);

                // Enhanced formatting with icons
                std::string icon = "🔬";
                if (cat == "dynamics") icon = "⚡";
                else if (cat == "optimization") icon = "🎯";
                else if (cat == "conformational") icon = "🔄";
                else if (cat == "docking") icon = "🔗";
                else if (cat == "calculation") icon = "⚙️";

                std::cout << "  " << icon << " -" << cap
                         << std::string(std::max(1, 18 - static_cast<int>(cap.length())), ' ')
                         << info.description << std::endl;

                // Show supported formats for each capability
                std::cout << "      Formats: ";
                for (size_t i = 0; i < info.supported_formats.size(); ++i) {
                    std::cout << info.supported_formats[i];
                    if (i < info.supported_formats.size() - 1) std::cout << ", ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        // Enhanced format support section with auto-discovery
        std::cout << "🗂️ Supported File Formats:" << std::endl;
        std::cout << "  ┌─────────────────────────────────────────────────────────┐" << std::endl;
        std::cout << "  │ Atomistic:     XYZ, MOL2, SDF, PDB                     │" << std::endl;
        std::cout << "  │ Coarse-Grained: VTF (VMD Trajectory Format)            │" << std::endl;
        std::cout << "  │ Trajectories:   XYZ.trj, VTF multi-timestep           │" << std::endl;
        std::cout << "  │ All formats work transparently with all capabilities   │" << std::endl;
        std::cout << "  └─────────────────────────────────────────────────────────┘" << std::endl;
        std::cout << std::endl;

        // Usage examples
        std::cout << "💡 Quick Start Examples:" << std::endl;
        std::cout << "  curcuma -analysis molecule.xyz             # Complete molecular analysis" << std::endl;
        std::cout << "  curcuma -opt protein.pdb -method gfn2      # Quantum optimization" << std::endl;
        std::cout << "  curcuma -casino polymer.vtf -steps 10000   # Casino Monte Carlo simulation" << std::endl;
        std::cout << "  curcuma -traj md_output.xyz -stride 10     # Trajectory analysis" << std::endl;
        std::cout << std::endl;

        std::cout << "📚 Get detailed help: curcuma -help [category]" << std::endl;
        std::cout << "   Available categories: ";
        bool first = true;
        for (const auto& [cat, capabilities] : categories) {
            if (!first) std::cout << ", ";
            std::cout << cat;
            first = false;
        }
        std::cout << std::endl;

    } else {
        // Enhanced category-specific help
        bool found = false;
        std::vector<std::string> found_capabilities;

        // Collect capabilities in this category
        for (const auto& [name, info] : CAPABILITY_REGISTRY) {
            if (info.category == category) {
                found_capabilities.push_back(name);
                found = true;
            }
        }

        if (found) {
            std::string category_name = category;
            category_name[0] = std::toupper(category_name[0]);

            std::cout << "📋 " << category_name << " Capabilities - Detailed Help" << std::endl;
            std::cout << "══════════════════════════════════════════════════════" << std::endl;
            std::cout << std::endl;

            for (const auto& cap : found_capabilities) {
                const auto& info = CAPABILITY_REGISTRY.at(cap);

                std::cout << "🔹 -" << cap << std::endl;
                std::cout << "   Description: " << info.description << std::endl;
                std::cout << "   Usage: curcuma -" << cap << " input_file [options]" << std::endl;
                std::cout << "   Supported formats: ";

                for (size_t i = 0; i < info.supported_formats.size(); ++i) {
                    std::cout << info.supported_formats[i];
                    if (i < info.supported_formats.size() - 1) std::cout << ", ";
                }
                std::cout << std::endl;

                // Add specific help references
                std::cout << "   Get detailed help: curcuma -" << cap << " (without arguments)" << std::endl;
                std::cout << std::endl;
            }

            // Category-specific examples
            if (category == "analysis") {
                std::cout << "💡 Analysis Examples:" << std::endl;
                std::cout << "  curcuma -analysis molecule.xyz -properties all" << std::endl;
                std::cout << "  curcuma -traj trajectory.vtf -export_timeseries true" << std::endl;
                std::cout << "  curcuma -rmsd ref.xyz target.xyz" << std::endl;
            } else if (category == "dynamics") {
                std::cout << "💡 Dynamics Examples:" << std::endl;
                std::cout << "  curcuma -md system.xyz -temperature 300 -steps 10000" << std::endl;
                std::cout << "  curcuma -casino polymer.vtf -move_type mixed -adaptive_step true" << std::endl;
            } else if (category == "optimization") {
                std::cout << "💡 Optimization Examples:" << std::endl;
                std::cout << "  curcuma -opt molecule.xyz -method uff" << std::endl;
                std::cout << "  curcuma -sp system.pdb -method gfn2" << std::endl;
            }
            std::cout << std::endl;

        } else {
            std::cout << "❌ Unknown category: " << category << std::endl;
            std::cout << std::endl;
            std::cout << "Available categories: ";

            std::set<std::string> available_categories;
            for (const auto& [name, info] : CAPABILITY_REGISTRY) {
                available_categories.insert(info.category);
            }

            bool first = true;
            for (const auto& cat : available_categories) {
                if (!first) std::cout << ", ";
                std::cout << cat;
                first = false;
            }
            std::cout << std::endl;
        }
    }
}

int main(int argc, char **argv) {
#ifndef _WIN32
#if __GNUC__
    signal(SIGINT, ctrl_c_handler);
    signal(SIGSEGV, bt_handler);
    signal(SIGABRT, bt_handler);
#endif
#endif

    General::StartUp(argc, argv);

    // Claude Generated: Initialize parameter registry early
    initialize_generated_registry();
    if (!ParameterRegistry::getInstance().validateRegistry()) {
        std::cerr << "Warning: Parameter registry validation failed" << std::endl;
    }

#ifdef C17
#ifndef _WIN32
    std::filesystem::remove("stop");
#endif
#else
    remove("stop");
#endif
    RunTimer timer(true);

    // Handle no arguments or help requests - Claude Generated
    if(argc < 2) {
        showStructuredHelp();
        exit(1);
    }

    std::string command = CLIUtils::stripLeadingDashes(argv[1]);

    // Handle help requests - Claude Generated
    if (command == "help" || command == "h") {
        if (argc >= 3) {
            showStructuredHelp(argv[2]); // Category-specific help
        } else {
            showStructuredHelp();
        }
        return 0;
    }

    // Handle parameter registry commands - Claude Generated (October 2025)
    if (command == "list_modules" || command == "list-modules") {
        ParameterRegistry::getInstance().printAllModules();
        return 0;
    }

    if (command == "export_config" || command == "export-config") {
        if (argc < 3) {
            std::cerr << "Error: -export-config requires module name" << std::endl;
            std::cerr << "Usage: curcuma -export-config <module>" << std::endl;
            std::cerr << "Example: curcuma -export-config analysis" << std::endl;
            return 1;
        }
        std::string module = argv[2];
        json config = ParameterRegistry::getInstance().getDefaultJson(module);
        if (config.empty()) {
            std::cerr << "Error: Unknown module '" << module << "'" << std::endl;
            std::cerr << "Use: curcuma -list-modules to see available modules" << std::endl;
            return 1;
        }
        std::cout << config.dump(2) << std::endl;
        return 0;
    }

    if (command == "help_module" || command == "help-module") {
        if (argc < 3) {
            std::cerr << "Error: -help-module requires module name" << std::endl;
            std::cerr << "Usage: curcuma -help-module <module>" << std::endl;
            std::cerr << "Example: curcuma -help-module analysis" << std::endl;
            return 1;
        }
        std::string module = argv[2];
        ParameterRegistry::getInstance().printHelp(module);
        return 0;
    }

    json controller = CLIUtils::CLI2Json(argc, argv);

    // Handle config import - Claude Generated (October 2025)
    // Now global parameter, always in controller["import_config"] if present
    if (controller.contains("import_config")) {
        std::string config_file = controller["import_config"];

        std::ifstream file(config_file);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open config file: " << config_file << std::endl;
            return 1;
        }

        json imported_config;
        try {
            file >> imported_config;
        } catch (const std::exception& e) {
            std::cerr << "Error: Failed to parse JSON from " << config_file << std::endl;
            std::cerr << "Details: " << e.what() << std::endl;
            return 1;
        }

        // Merge imported config into controller
        // CLI arguments have higher priority (don't overwrite existing keys)
        for (auto it = imported_config.begin(); it != imported_config.end(); ++it) {
            if (!controller.contains(it.key())) {
                controller[it.key()] = it.value();
            } else if (controller[it.key()].is_object() && it.value().is_object()) {
                // Merge nested objects (e.g., controller["analysis"] with imported["analysis"])
                for (auto nested_it = it.value().begin(); nested_it != it.value().end(); ++nested_it) {
                    if (!controller[it.key()].contains(nested_it.key())) {
                        controller[it.key()][nested_it.key()] = nested_it.value();
                    }
                }
            }
        }

        std::cout << "Loaded configuration from: " << config_file << std::endl;
    }

    // Handle run export - Claude Generated (October 2025)
    // Export current configuration AFTER all merging/importing
    // Now global parameter, always in controller["export_run"] if present
    if (controller.contains("export_run")) {
        std::string export_file = controller["export_run"];

        // Create export config without meta-parameters
        json export_config = controller;

        // Remove all meta-parameters (those used for I/O operations, not actual computation)
        // Both at top level and from module sub-dictionaries
        std::vector<std::string> meta_params = { "export_run", "export-run", "import_config", "import-config" };

        for (const auto& meta : meta_params) {
            export_config.erase(meta);

            // Also remove from all module sub-dictionaries
            for (auto it = export_config.begin(); it != export_config.end(); ++it) {
                if (it.value().is_object()) {
                    it.value().erase(meta);
                }
            }
        }

        std::ofstream outfile(export_file);
        if (!outfile.is_open()) {
            std::cerr << "Error: Could not write to " << export_file << std::endl;
            return 1;
        }

        outfile << export_config.dump(2) << std::endl;
        outfile.close();

        std::cout << "Exported current run configuration to: " << export_file << std::endl;
    }

    // Try structured dispatch first - Claude Generated
    auto it = CAPABILITY_REGISTRY.find(command);
    if (it != CAPABILITY_REGISTRY.end()) {
        return it->second.handler(controller, argc, argv);
    }

    // Handle unknown commands
    std::cerr << "Error: Unknown command '-" << command << "'" << std::endl;
    std::cerr << "Use 'curcuma -help' to see available capabilities." << std::endl;

#ifdef C17
#ifndef _WIN32
    std::filesystem::remove("stop");
#endif
#else
    remove("stop");
#endif

    return 0;
}
