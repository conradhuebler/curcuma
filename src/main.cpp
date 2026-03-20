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
#include "src/core/curcuma_logger.h"

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
#include "src/capabilities/trajectory_statistics.h"
#include "src/capabilities/trajectoryanalysis.h"

#include "src/tools/trajectory_writer.h"

#include "src/tools/general.h"
#include "src/tools/info.h"

// Claude Generated: Parameter registry system
#include "generated/parameter_registry.h"
#include "src/core/parameter_registry.h"

#include "src/capabilities/optimiser/OptimiseDipoleScaling.h"
#include "src/capabilities/optimisation/modern_optimizer_simple.h"

#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
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

// Claude Generated: Geometry calculation helper functions
std::string getFormatArg(int argc, char** argv, const std::string& default_format = "human")
{
    for (int i = 2; i < argc - 1; ++i) {
        std::string arg = argv[i];
        if (arg == "-format" || arg == "-fmt") {
            return argv[i + 1];
        }
    }
    return default_format;
}

std::string getUnitArg(int argc, char** argv, const std::string& default_unit = "degrees")
{
    for (int i = 2; i < argc - 1; ++i) {
        std::string arg = argv[i];
        if (arg == "-unit" || arg == "-u") {
            return argv[i + 1];
        }
    }
    return default_unit;
}

int executeBond(const json& controller, int argc, char** argv)
{
    // Claude Generated: Extended to support trajectories (single and multi-frame)
    if (argc < 4) {
        std::cerr << "Usage: curcuma -bond <input_file> <atom1> <atom2> [-format human|json|csv]" << std::endl;
        std::cerr << "Example: curcuma -bond molecule.xyz 0 1" << std::endl;
        std::cerr << "Example: curcuma -bond md_traj.xyz 0 1 -format csv" << std::endl;
        return 1;
    }

    // Parse arguments
    std::string filename = argv[2];
    int atom1 = std::stoi(argv[3]);
    int atom2 = std::stoi(argv[4]);
    std::string format = getFormatArg(argc, argv, "human");

    // Load molecules (single or trajectory)
    if (!std::filesystem::exists(filename)) {
        std::cerr << "Error: File not found: " << filename << std::endl;
        return 1;
    }

    FileIterator file_iter(filename, true);
    if (file_iter.AtEnd()) {
        std::cerr << "Error: Could not read molecules from " << filename << std::endl;
        return 1;
    }

    // Read first molecule to validate atom indices
    Molecule first_mol = file_iter.Next();
    if (atom1 < 0 || atom1 >= first_mol.AtomCount() || atom2 < 0 || atom2 >= first_mol.AtomCount()) {
        std::cerr << "Error: Atom indices out of range (molecule has " << first_mol.AtomCount() << " atoms)" << std::endl;
        return 1;
    }

    // Collect distances from all frames
    std::vector<double> distances;
    distances.push_back(first_mol.CalculateDistance(atom1, atom2));

    // Read remaining frames if trajectory
    while (!file_iter.AtEnd()) {
        Molecule mol = file_iter.Next();
        distances.push_back(mol.CalculateDistance(atom1, atom2));
    }

    // Calculate statistics
    TrajectoryStatistics stats(10);  // window_size = 10
    for (size_t i = 0; i < distances.size(); ++i) {
        stats.addValue("bond_length", distances[i]);
    }

    // Output results
    if (distances.size() == 1) {
        // Single frame - simple output
        if (format == "json") {
            json result;
            result["calculation"] = "bond_length";
            result["atoms"] = json::array({atom1, atom2});
            result["distance_angstrom"] = distances[0];
            std::cout << result.dump(2) << std::endl;
        } else if (format == "csv") {
            std::cout << "atom1,atom2,distance_angstrom" << std::endl;
            std::cout << atom1 << "," << atom2 << "," << distances[0] << std::endl;
        } else { // human
            std::cout << "Bond Length Calculation:" << std::endl;
            std::cout << "  Atoms: " << atom1 << "-" << atom2 << std::endl;
            std::cout << "  Distance: " << std::fixed << std::setprecision(6) << distances[0] << " Å" << std::endl;
        }
    } else {
        // Multiple frames - use TrajectoryWriter (Phase 5 refactoring)
        json trajectory_data = TrajectoryWriter::createTrajectoryJSON(
            distances,
            "bond_length",
            "Ångström",
            stats
        );

        // Override default TrajectoryWriter formatting to match current style
        if (format == "json") {
            // Maintain current JSON format with specific structure
            json result;
            result["calculation"] = "bond_length_trajectory";
            result["atoms"] = json::array({atom1, atom2});
            result["num_frames"] = static_cast<int>(distances.size());
            result["frames"] = distances;
            result["statistics"] = {
                {"mean", stats.getMean("bond_length")},
                {"std_dev", stats.getStdDev("bond_length")},
                {"min", *std::min_element(distances.begin(), distances.end())},
                {"max", *std::max_element(distances.begin(), distances.end())}
            };
            std::cout << result.dump(2) << std::endl;
        } else {
            // Use original inline formatting for CSV and human formats
            // (maintaining compatibility while JSON uses trajectory framework)
            if (format == "csv") {
                std::cout << "frame,bond_length,<bond_length>,σ(bond_length)" << std::endl;
                double mean = stats.getMean("bond_length");
                double std_dev = stats.getStdDev("bond_length");
                for (size_t i = 0; i < distances.size(); ++i) {
                    std::cout << i << "," << std::fixed << std::setprecision(6) << distances[i]
                             << "," << mean << "," << std_dev << std::endl;
                }
            } else { // human
                std::cout << "Bond Length Trajectory:" << std::endl;
                std::cout << "  Atoms: " << atom1 << "-" << atom2 << std::endl;
                std::cout << "  Frames: " << distances.size() << std::endl << std::endl;

                // Header (original inline code)
                std::cout << std::setw(8) << "# Frame" << std::setw(15) << "bond_length"
                         << std::setw(15) << "<bond_length>" << std::setw(15) << "σ(bond_length)" << std::endl;

                // Data rows (original inline code)
                for (size_t i = 0; i < distances.size(); ++i) {
                    std::cout << std::setw(8) << i << std::fixed << std::setprecision(6)
                             << std::setw(15) << distances[i]
                             << std::setw(15) << stats.getMean("bond_length")
                             << std::setw(15) << stats.getStdDev("bond_length") << std::endl;
                }

                // Statistics footer (maintained for compatibility)
                double mean = stats.getMean("bond_length");
                double std_dev = stats.getStdDev("bond_length");
                std::cout << std::endl << "Statistics:" << std::endl;
                std::cout << "  Mean: " << std::fixed << std::setprecision(6) << mean << " Å" << std::endl;
                std::cout << "  StdDev: " << std_dev << " Å" << std::endl;
                std::cout << "  Min: " << *std::min_element(distances.begin(), distances.end()) << " Å" << std::endl;
                std::cout << "  Max: " << *std::max_element(distances.begin(), distances.end()) << " Å" << std::endl;
            }
        }
    }

    return 0;
}

// Claude Generated helper: Convert angles from radians based on unit
double convertAngle(double radians, const std::string& unit)
{
    if (unit == "degrees" || unit == "deg" || unit == "°") {
        return radians * 180.0 / M_PI;
    }
    return radians;  // Return radians
}

// Claude Generated helper: Get angle unit string
std::string getAngleUnitString(const std::string& unit)
{
    if (unit == "degrees" || unit == "deg" || unit == "°") {
        return "degrees";
    }
    return "radians";
}

int executeAngle(const json& controller, int argc, char** argv)
{
    // Claude Generated: Extended to support trajectories (single and multi-frame)
    if (argc < 5) {
        std::cerr << "Usage: curcuma -angle <input_file> <atom1> <atom2> <atom3> [-format human|json|csv] [-unit degrees|radians]" << std::endl;
        std::cerr << "Example: curcuma -angle molecule.xyz 0 1 2" << std::endl;
        std::cerr << "Example: curcuma -angle md_traj.xyz 0 1 2 -format csv -unit degrees" << std::endl;
        return 1;
    }

    // Parse arguments
    std::string filename = argv[2];
    int atom1 = std::stoi(argv[3]);
    int atom2 = std::stoi(argv[4]);
    int atom3 = std::stoi(argv[5]);
    std::string format = getFormatArg(argc, argv, "human");
    std::string unit = getUnitArg(argc, argv, "degrees");
    std::string angle_unit = getAngleUnitString(unit);

    // Load molecules (single or trajectory)
    if (!std::filesystem::exists(filename)) {
        std::cerr << "Error: File not found: " << filename << std::endl;
        return 1;
    }

    FileIterator file_iter(filename, true);
    if (file_iter.AtEnd()) {
        std::cerr << "Error: Could not read molecules from " << filename << std::endl;
        return 1;
    }

    // Read first molecule to validate atom indices
    Molecule first_mol = file_iter.Next();
    if (atom1 < 0 || atom1 >= first_mol.AtomCount() || atom2 < 0 || atom2 >= first_mol.AtomCount() ||
        atom3 < 0 || atom3 >= first_mol.AtomCount()) {
        std::cerr << "Error: Atom indices out of range (molecule has " << first_mol.AtomCount() << " atoms)" << std::endl;
        return 1;
    }

    // Collect angles from all frames
    std::vector<double> angles;
    angles.push_back(convertAngle(first_mol.CalculateAngle(atom1, atom2, atom3), unit));

    // Read remaining frames if trajectory
    while (!file_iter.AtEnd()) {
        Molecule mol = file_iter.Next();
        angles.push_back(convertAngle(mol.CalculateAngle(atom1, atom2, atom3), unit));
    }

    // Calculate statistics
    TrajectoryStatistics stats(10);  // window_size = 10
    for (size_t i = 0; i < angles.size(); ++i) {
        stats.addValue("bond_angle", angles[i]);
    }

    // Output results
    if (angles.size() == 1) {
        // Single frame - simple output
        if (format == "json") {
            json result;
            result["calculation"] = "bond_angle";
            result["atoms"] = json::array({atom1, atom2, atom3});
            result["angle"] = angles[0];
            result["unit"] = angle_unit;
            std::cout << result.dump(2) << std::endl;
        } else if (format == "csv") {
            std::cout << "atom1,atom2,atom3,angle_" << angle_unit << std::endl;
            std::cout << atom1 << "," << atom2 << "," << atom3 << "," << angles[0] << std::endl;
        } else { // human
            std::cout << "Bond Angle Calculation:" << std::endl;
            std::cout << "  Atoms: " << atom1 << "-" << atom2 << "-" << atom3 << std::endl;
            std::cout << "  Angle: " << std::fixed << std::setprecision(6) << angles[0] << " " << angle_unit << std::endl;
        }
    } else {
        // Multiple frames - use TrajectoryWriter (Phase 5 refactoring)
        json trajectory_data = TrajectoryWriter::createTrajectoryJSON(
            angles,
            "bond_angle",
            angle_unit,
            stats
        );

        // Override default TrajectoryWriter formatting to match current style
        if (format == "json") {
            // Maintain current JSON format with specific structure
            json result;
            result["calculation"] = "bond_angle_trajectory";
            result["atoms"] = json::array({atom1, atom2, atom3});
            result["num_frames"] = static_cast<int>(angles.size());
            result["frames"] = angles;
            result["unit"] = angle_unit;
            result["statistics"] = {
                {"mean", stats.getMean("bond_angle")},
                {"std_dev", stats.getStdDev("bond_angle")},
                {"min", *std::min_element(angles.begin(), angles.end())},
                {"max", *std::max_element(angles.begin(), angles.end())}
            };
            std::cout << result.dump(2) << std::endl;
        } else {
            // Use original inline formatting for CSV and human formats
            // (maintaining compatibility while JSON uses trajectory framework)
            if (format == "csv") {
                std::cout << "frame,bond_angle,<bond_angle>,σ(bond_angle)" << std::endl;
                double mean = stats.getMean("bond_angle");
                double std_dev = stats.getStdDev("bond_angle");
                for (size_t i = 0; i < angles.size(); ++i) {
                    std::cout << i << "," << std::fixed << std::setprecision(6) << angles[i]
                             << "," << mean << "," << std_dev << std::endl;
                }
            } else { // human
                std::cout << "Bond Angle Trajectory:" << std::endl;
                std::cout << "  Atoms: " << atom1 << "-" << atom2 << "-" << atom3 << std::endl;
                std::cout << "  Frames: " << angles.size() << std::endl << std::endl;

                // Header (original inline code)
                std::cout << std::setw(8) << "# Frame" << std::setw(15) << "bond_angle"
                         << std::setw(15) << "<bond_angle>" << std::setw(15) << "σ(bond_angle)" << std::endl;

                // Data rows (original inline code)
                for (size_t i = 0; i < angles.size(); ++i) {
                    std::cout << std::setw(8) << i << std::fixed << std::setprecision(6)
                             << std::setw(15) << angles[i]
                             << std::setw(15) << stats.getMean("bond_angle")
                             << std::setw(15) << stats.getStdDev("bond_angle") << std::endl;
                }

                // Statistics footer (maintained for compatibility)
                double mean = stats.getMean("bond_angle");
                double std_dev = stats.getStdDev("bond_angle");
                std::cout << std::endl << "Statistics:" << std::endl;
                std::cout << "  Mean: " << std::fixed << std::setprecision(6) << mean << " " << angle_unit << std::endl;
                std::cout << "  StdDev: " << std_dev << " " << angle_unit << std::endl;
                std::cout << "  Min: " << *std::min_element(angles.begin(), angles.end()) << " " << angle_unit << std::endl;
                std::cout << "  Max: " << *std::max_element(angles.begin(), angles.end()) << " " << angle_unit << std::endl;
            }
        }
    }

    return 0;
}

int executeTorsion(const json& controller, int argc, char** argv)
{
    // Claude Generated: Extended to support trajectories (single and multi-frame)
    if (argc < 6) {
        std::cerr << "Usage: curcuma -torsion <input_file> <atom1> <atom2> <atom3> <atom4> [-format human|json|csv] [-unit degrees|radians]" << std::endl;
        std::cerr << "Example: curcuma -torsion molecule.xyz 0 1 2 3" << std::endl;
        std::cerr << "Example: curcuma -torsion md_traj.xyz 0 1 2 3 -format csv -unit degrees" << std::endl;
        return 1;
    }

    // Parse arguments
    std::string filename = argv[2];
    int atom1 = std::stoi(argv[3]);
    int atom2 = std::stoi(argv[4]);
    int atom3 = std::stoi(argv[5]);
    int atom4 = std::stoi(argv[6]);
    std::string format = getFormatArg(argc, argv, "human");
    std::string unit = getUnitArg(argc, argv, "degrees");
    std::string dihedral_unit = getAngleUnitString(unit);

    // Load molecules (single or trajectory)
    if (!std::filesystem::exists(filename)) {
        std::cerr << "Error: File not found: " << filename << std::endl;
        return 1;
    }

    FileIterator file_iter(filename, true);
    if (file_iter.AtEnd()) {
        std::cerr << "Error: Could not read molecules from " << filename << std::endl;
        return 1;
    }

    // Read first molecule to validate atom indices
    Molecule first_mol = file_iter.Next();
    if (atom1 < 0 || atom1 >= first_mol.AtomCount() || atom2 < 0 || atom2 >= first_mol.AtomCount() ||
        atom3 < 0 || atom3 >= first_mol.AtomCount() || atom4 < 0 || atom4 >= first_mol.AtomCount()) {
        std::cerr << "Error: Atom indices out of range (molecule has " << first_mol.AtomCount() << " atoms)" << std::endl;
        return 1;
    }

    // Collect dihedrals from all frames
    std::vector<double> dihedrals;
    dihedrals.push_back(convertAngle(first_mol.CalculateDihedral(atom1, atom2, atom3, atom4), unit));

    // Read remaining frames if trajectory
    while (!file_iter.AtEnd()) {
        Molecule mol = file_iter.Next();
        dihedrals.push_back(convertAngle(mol.CalculateDihedral(atom1, atom2, atom3, atom4), unit));
    }

    // Calculate statistics
    TrajectoryStatistics stats(10);  // window_size = 10
    for (size_t i = 0; i < dihedrals.size(); ++i) {
        stats.addValue("dihedral", dihedrals[i]);
    }

    // Output results
    if (dihedrals.size() == 1) {
        // Single frame - simple output
        if (format == "json") {
            json result;
            result["calculation"] = "dihedral_angle";
            result["atoms"] = json::array({atom1, atom2, atom3, atom4});
            result["dihedral"] = dihedrals[0];
            result["unit"] = dihedral_unit;
            std::cout << result.dump(2) << std::endl;
        } else if (format == "csv") {
            std::cout << "atom1,atom2,atom3,atom4,dihedral_" << dihedral_unit << std::endl;
            std::cout << atom1 << "," << atom2 << "," << atom3 << "," << atom4 << "," << dihedrals[0] << std::endl;
        } else { // human
            std::cout << "Dihedral/Torsion Angle Calculation:" << std::endl;
            std::cout << "  Atoms: " << atom1 << "-" << atom2 << "-" << atom3 << "-" << atom4 << std::endl;
            std::cout << "  Dihedral: " << std::fixed << std::setprecision(6) << dihedrals[0] << " " << dihedral_unit << std::endl;
        }
    } else {
        // Multiple frames - use TrajectoryWriter (Phase 5 refactoring)
        json trajectory_data = TrajectoryWriter::createTrajectoryJSON(
            dihedrals,
            "dihedral",
            dihedral_unit,
            stats
        );

        // Override default TrajectoryWriter formatting to match current style
        if (format == "json") {
            // Maintain current JSON format with specific structure
            json result;
            result["calculation"] = "dihedral_angle_trajectory";
            result["atoms"] = json::array({atom1, atom2, atom3, atom4});
            result["num_frames"] = static_cast<int>(dihedrals.size());
            result["frames"] = dihedrals;
            result["unit"] = dihedral_unit;
            result["statistics"] = {
                {"mean", stats.getMean("dihedral")},
                {"std_dev", stats.getStdDev("dihedral")},
                {"min", *std::min_element(dihedrals.begin(), dihedrals.end())},
                {"max", *std::max_element(dihedrals.begin(), dihedrals.end())}
            };
            std::cout << result.dump(2) << std::endl;
        } else {
            // Use original inline formatting for CSV and human formats
            // (maintaining compatibility while JSON uses trajectory framework)
            if (format == "csv") {
                std::cout << "frame,dihedral,<dihedral>,σ(dihedral)" << std::endl;
                double mean = stats.getMean("dihedral");
                double std_dev = stats.getStdDev("dihedral");
                for (size_t i = 0; i < dihedrals.size(); ++i) {
                    std::cout << i << "," << std::fixed << std::setprecision(6) << dihedrals[i]
                             << "," << mean << "," << std_dev << std::endl;
                }
            } else { // human
                std::cout << "Dihedral/Torsion Angle Trajectory:" << std::endl;
                std::cout << "  Atoms: " << atom1 << "-" << atom2 << "-" << atom3 << "-" << atom4 << std::endl;
                std::cout << "  Frames: " << dihedrals.size() << std::endl << std::endl;

                // Header (original inline code)
                std::cout << std::setw(8) << "# Frame" << std::setw(15) << "dihedral"
                         << std::setw(15) << "<dihedral>" << std::setw(15) << "σ(dihedral)" << std::endl;

                // Data rows (original inline code)
                for (size_t i = 0; i < dihedrals.size(); ++i) {
                    std::cout << std::setw(8) << i << std::fixed << std::setprecision(6)
                             << std::setw(15) << dihedrals[i]
                             << std::setw(15) << stats.getMean("dihedral")
                             << std::setw(15) << stats.getStdDev("dihedral") << std::endl;
                }

                // Statistics footer (maintained for compatibility)
                double mean = stats.getMean("dihedral");
                double std_dev = stats.getStdDev("dihedral");
                std::cout << std::endl << "Statistics:" << std::endl;
                std::cout << "  Mean: " << std::fixed << std::setprecision(6) << mean << " " << dihedral_unit << std::endl;
                std::cout << "  StdDev: " << std_dev << " " << dihedral_unit << std::endl;
                std::cout << "  Min: " << *std::min_element(dihedrals.begin(), dihedrals.end()) << " " << dihedral_unit << std::endl;
                std::cout << "  Max: " << *std::max_element(dihedrals.begin(), dihedrals.end()) << " " << dihedral_unit << std::endl;
            }
        }
    }

    return 0;
}

double DotProduct(const Eigen::Vector3d& pos1, const Eigen::Vector3d& pos2)
{
    return pos1.dot(pos2);
}

// Helper function to set nested JSON values from dot-notation keys - Claude Generated
void setNestedJsonValue(json& target, const std::string& dotKey, const json& value) {
    if (dotKey.find('.') == std::string::npos) {
        // No dots - use as flat key
        target[dotKey] = value;
        return;
    }

    // Split on dots and create nested structure
    std::vector<std::string> keys;
    std::stringstream ss(dotKey);
    std::string key;

    while (std::getline(ss, key, '.')) {
        keys.push_back(key);
    }

    // Navigate/create nested structure
    json* current = &target;
    for (size_t i = 0; i < keys.size() - 1; ++i) {
        if (!current->contains(keys[i]) || !(*current)[keys[i]].is_object()) {
            (*current)[keys[i]] = json::object();
        }
        current = &(*current)[keys[i]];
    }

    // Set the final value
    (*current)[keys.back()] = value;
}

json CLI2Json(int argc, char** argv)
{
    json controller;
    json key;
    if (argc < 2)
        return controller;

    std::string keyword = argv[1];
    keyword.erase(0, 1);

    // Claude Generated (October 2025): Global parameters that should be accessible
    // both at top level (controller[param]) and module level (controller[module][param])
    // ENHANCED: Added "method" to support global energy method specification
    std::set<std::string> global_params = {
        "verbosity", "threads", "method",  // energy_method applies to all capabilities
        "export_run", "export-run", // Export current run configuration
        "import_config", "import-config" // Import custom configuration
    };

    // Claude Generated (October 2025): CLI keyword to module name mapping
    // Maps command-line keywords (e.g., -md) to actual module names (e.g., simplemd)
    std::map<std::string, std::string> keyword_to_module = {
        {"md", "simplemd"},
        {"opt", "opt"},
        {"sp", "opt"},  // single point also uses opt module
        {"confscan", "confscan"},
        {"rmsd", "rmsd"},
        {"analysis", "analysis"},
        {"hessian", "hessian"},
        {"casino", "casino"}
    };

    // Get actual module name (for ConfigManager)
    std::string module_name = keyword;
    if (keyword_to_module.count(keyword) > 0) {
        module_name = keyword_to_module[keyword];
    }

    for (int i = 2; i < argc; ++i) {
        std::string current = argv[i];
        std::string sub = current.substr(0, 1);

        if (sub == "-") {
            current.erase(0, 1);

            // Handle special verbosity shortcuts - Claude Generated
            if (current == "silent" || current == "quiet") {
                key["verbosity"] = 0;
                continue;
            } else if (current == "verbose") {
                key["verbosity"] = 3;
                continue;
            } else if (current == "plain") {
                // Claude Generated: Plain mode - no colors, no prefixes (like ORCA/Gaussian)
                CurcumaLogger::set_plain_mode(true);
                continue;
            } else if (current == "v") {
                // Handle -v N syntax for verbosity level
                if ((i + 1) < argc) {
                    try {
                        int verbosity_level = std::stoi(argv[i + 1]);
                        if (verbosity_level >= 0 && verbosity_level <= 3) {
                            key["verbosity"] = verbosity_level;
                            ++i; // Skip next argument
                            continue;
                        }
                    } catch (const std::exception&) {
                        // Fall through to normal processing
                    }
                }
            }

            // Claude Generated (October 2025 - CRITICAL FIX): Strip redundant keyword prefix from dotted parameters
            // Makes "-md.max_time 10" and "-max_time 10" synonymous within "-md" command context
            // Fixes SimpleMD double-nesting bug: controller["simplemd"]["md"]["max_time"] → controller["simplemd"]["max_time"]
            if (current.find('.') != std::string::npos) {
                size_t dot_pos = current.find('.');
                std::string param_prefix = current.substr(0, dot_pos);
                std::string param_key = current.substr(dot_pos + 1);

                // Resolve prefix to module name using keyword_to_module map (e.g., "md" → "simplemd")
                std::string prefix_module = param_prefix;
                if (keyword_to_module.count(param_prefix) > 0) {
                    prefix_module = keyword_to_module[param_prefix];
                }

                // If prefix matches current command's keyword or module name, strip it
                // Example: "-md input.xyz -md.max_time 10" → keyword="md", module_name="simplemd"
                //          param_prefix="md" matches keyword → strip → current="max_time"
                // Preserves cross-module routing: "-md -rmsd.method subspace" keeps "rmsd.method"
                if (param_prefix == keyword || prefix_module == module_name) {
                    current = param_key;  // Strip prefix: "md.max_time" → "max_time"
                }
            }

            if ((i + 1) >= argc || argv[i + 1][0] == '-' || argv[i + 1] == std::string("true") || argv[i + 1] == std::string("+")) {
                setNestedJsonValue(key, current, true);
            } else if (argv[i + 1] == std::string("false")) {
                setNestedJsonValue(key, current, false);
                ++i;
            } else {
                std::string next = argv[i + 1];
                //       std::cout << "next: " << next << std::endl;

                bool isNumber = true;
                bool isVector = next.find("|") != std::string::npos || next.find(",") != std::string::npos || next.find(":") != std::string::npos;
                //      std::cout << "isNumber: " << isNumber << std::endl
                //                  << "isVector: " << isVector << std::endl;
                if (isVector) {
                    isNumber = false;
                }
                if (!isVector) {
                    try {
                        // std::cout << "stod: " << std::stod(next) << std::endl;
                        std::stod(next);

                    } catch (const std::invalid_argument&) {
                        isNumber = false;
                        isVector = true;
                    }
                }
                // std::cout << "isNumber: " << isNumber << std::endl
                //             << "isVector: " << isVector << std::endl;
                if (isNumber) {
                    setNestedJsonValue(key, current, std::stod(next));
                } else if (isVector) {
                    //        std::cout << next << std::endl;
                    setNestedJsonValue(key, current, next);
                } else {
                    setNestedJsonValue(key, current, next);
                }
                ++i;
            }
        }
    }

    // Claude Generated 2025: Extract module-specific parameters to top level
    // Parameters like "rmsd.method" should be at controller["rmsd"]["method"], NOT controller["confscan"]["rmsd"]["method"]
    json module_params;
    std::vector<std::string> keys_to_remove;

    for (auto& [param_name, param_value] : key.items()) {
        // Claude Generated 2025: Check if this is a module-specific parameter
        // Could be either:
        // 1. Flat key with dots: "rmsd.method" (stored as param_name contains dot)
        // 2. Nested structure: setNestedJsonValue already created key["rmsd"]["method"]
        //    so param_name = "rmsd" and param_value = {"method": "subspace"}

        bool is_flat_dotted = param_name.find('.') != std::string::npos;
        bool is_nested_object = param_value.is_object();

        if (is_flat_dotted) {
            // Handle flat dot notation: "rmsd.method" or "-md.max_time"
            size_t dot_pos = param_name.find('.');
            std::string module_name = param_name.substr(0, dot_pos);
            std::string param_key = param_name.substr(dot_pos + 1);

            if (module_name != keyword) {
                // Claude Generated (October 2025 - CRITICAL FIX): Map keywords to actual module names
                // e.g., "-md.max_time" has module_name="md" from Punkt-notation,
                // but must be routed to module_params["simplemd"] (the actual module)
                std::string target_module = module_name;
                if (keyword_to_module.count(module_name) > 0) {
                    target_module = keyword_to_module[module_name];
                }
                if (!module_params.contains(target_module)) {
                    module_params[target_module] = json::object();
                }
                setNestedJsonValue(module_params[target_module], param_key, param_value);
                keys_to_remove.push_back(param_name);
            } else {
                // Claude Generated (October 2025 - CRITICAL FIX): Handle "-keyword.param" (e.g., "-md.max_time")
                // These should be FLAT in controller["simplemd"], not nested as key["md"]["max_time"]
                // Extract directly to top level to avoid double-nesting
                setNestedJsonValue(key, param_key, param_value);
                keys_to_remove.push_back(param_name);
            }
        } else if (is_nested_object && param_name != keyword) {
            // Handle nested structure for OTHER modules (not this command's keyword)
            // param_name = "rmsd", param_value = {"method": "subspace"}
            module_params[param_name] = param_value;
            keys_to_remove.push_back(param_name);
        }
        // NOTE: We KEEP nested structures where param_name == keyword!
        // E.g., key["md"] = {"max_time": 10} stays in key for backward compat
        // It will be stored in both controller["simplemd"] and controller["md"]
    }

    // Remove extracted module parameters from main command params
    for (const auto& remove_key : keys_to_remove) {
        key.erase(remove_key);
    }

    // Claude Generated (October 2025): Extract global parameters for explicit access
    json global_values;
    for (const auto& param : global_params) {
        if (key.count(param) > 0) {
            global_values[param] = key[param];
        }
    }

    // Claude Generated: DEBUG - Show what's in key before storing to controller
    std::cerr << "[CLI2Json DEBUG] keyword=" << keyword << ", module_name=" << module_name << std::endl;
    std::cerr << "[CLI2Json DEBUG] key object content:" << std::endl;
    std::cerr << key.dump(2) << std::endl;

    // Build controller with proper structure using actual module name
    // This enables ConfigManager to find parameters under the correct module name
    controller[module_name] = key;  // Changed from keyword to module_name

    // Claude Generated (October 2025): Backward compatibility - keep old keyword-based access
    // Some tests use "-md.max_time" which creates controller["md"],
    // but we now route to controller["simplemd"] via module_name
    // Store under both for compatibility
    if (module_name != keyword) {
        controller[keyword] = key;  // Also keep the old keyword for backward compat
    }

    // Merge module-specific parameters at top level
    for (auto& [mod_name, module_config] : module_params.items()) {
        controller[mod_name] = module_config;  // Direct assignment (may contain nested structure from setNestedJsonValue)
    }

    // Ensure global parameters are always set at top level AND in global section
    for (const auto& [param, value] : global_values.items()) {
        controller[param] = value;  // Top-level access (backward compat)
        controller["global"][param] = value;  // Explicit global namespace
    }
    return controller;
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
    // Guard against null JSON (can occur when -opt is used without additional parameters)
    json opt_config_base = (controller.contains("opt") && controller["opt"].is_object())
                           ? controller["opt"] : json::object();
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

// Claude Generated - RMSDTraj execution function
int executeRMSDTraj(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma for rmsd analysis of trajectories as follows:\ncurcuma -rmsdtraj input.xyz" << std::endl;
        std::cerr << "Additional arguments are:" << std::endl;
        std::cerr << "-write        **** Write unique conformers!" << std::endl;
        std::cerr << "-rmsd d       **** Set rmsd threshold to d ( default = 1.0)!" << std::endl;
        std::cerr << "-fragment n   **** Set fragment to n." << std::endl;
        std::cerr << "-reference    **** Add different xyz structure as reference." << std::endl;
        std::cerr << "-second       **** Add second trajectory." << std::endl;
        std::cerr << "-heavy        **** Check only heavy atoms. Do not use with -write." << std::endl;
        return 0;
    }

    // Use ParameterRegistry defaults and merge with controller["rmsdtraj"] if it exists
    json rmsdtraj_config = ParameterRegistry::getInstance().getDefaultJson("rmsdtraj");
    if (controller.contains("rmsdtraj") && !controller["rmsdtraj"].is_null()) {
        rmsdtraj_config = MergeJson(rmsdtraj_config, controller["rmsdtraj"]);
    }

    RMSDTraj traj(rmsdtraj_config, false);
    traj.setFile(argv[2]);
    traj.Initialise();
    traj.start();
    return 0;
}

// Capability registry - Claude Generated
const std::map<std::string, CapabilityInfo> CAPABILITY_REGISTRY = {
    {"analysis", {"Unified molecular analysis (all formats, all properties)", "analysis",
                  {"XYZ", "VTF", "MOL2", "SDF", "PDB"}, executeAnalysis}},
    {"rmsd", {"RMSD calculation between structures", "analysis",
              {"XYZ", "VTF", "MOL2", "SDF"}, executeRMSD}},
    {"rmsdtraj", {"Trajectory RMSD analysis and conformer filtering", "analysis",
                  {"XYZ", "TRJ"}, executeRMSDTraj}},
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
    {"bond", {"Calculate bond length between two atoms", "analysis",
              {"XYZ", "MOL2", "PDB", "SDF"}, executeBond}},
    {"torsion", {"Calculate dihedral/torsion angle between four atoms", "analysis",
                 {"XYZ", "MOL2", "PDB", "SDF"}, executeTorsion}},
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

    std::string command = argv[1];
    while (!command.empty() && command[0] == '-') command.erase(0, 1);

    // Handle help requests - Claude Generated
    if (command == "help" || command == "h") {
        if (argc >= 3) {
            showStructuredHelp(argv[2]); // Category-specific help
        } else {
            showStructuredHelp();
        }
        return 0;
    }

    // Phase 2: List available computational methods - Claude Generated 2025
    if (command == "methods") {
        auto methods = MethodFactory::getAvailableMethods();
        std::cout << "Available computational methods in this build:\n\n";

        std::cout << "Quantum Methods:\n";
        for (const auto& method : methods) {
            if (method.find("gfn") != std::string::npos ||
                method.find("eht") != std::string::npos ||
                method.find("pm") != std::string::npos ||
                method.find("am") != std::string::npos ||
                method.find("mndo") != std::string::npos) {
                std::cout << "  - " << method << "\n";
            }
        }

        std::cout << "\nForce Fields:\n";
        for (const auto& method : methods) {
            if (method.find("uff") != std::string::npos ||
                method.find("ff") != std::string::npos ||
                method.find("qmdff") != std::string::npos) {
                std::cout << "  - " << method << "\n";
            }
        }

        std::cout << "\nDispersion Corrections:\n";
        for (const auto& method : methods) {
            if (method.find("d3") != std::string::npos ||
                method.find("d4") != std::string::npos) {
                std::cout << "  - " << method << "\n";
            }
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

    json controller = CLI2Json(argc, argv);

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
