#include "polymerbuild.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energycalculator.h"
#include "curcumaopt.h"
#include "simplemd.h"
#include "src/tools/general.h"
#include "src/core/elements.h"

#include <regex>
#include <iostream>
#include <Eigen/Geometry>
#include <fmt/format.h>

PolymerBuild::PolymerBuild(const json& controller, bool silent)
    : m_config("polymerbuild", controller)
    , m_silent(silent)
{
    // Parse fragments map
    if (controller.contains("fragments") && controller["fragments"].is_object()) {
        for (auto& [name, path] : controller["fragments"].items()) {
            m_fragments[name] = path.get<std::string>();
        }
    } else {
        std::string fragments_str = m_config.get<std::string>("fragments", "{}");
        try {
            json frags = json::parse(fragments_str);
            for (auto& [name, path] : frags.items()) {
                m_fragments[name] = path.get<std::string>();
            }
        } catch (...) {
             if (fragments_str != "{}")
                CurcumaLogger::error("Failed to parse fragments JSON map: " + fragments_str);
        }
    }
}

PolymerBuild::~PolymerBuild() {}

void PolymerBuild::start()
{
    std::string sequence_str = m_config.get<std::string>("sequence", "");
    if (sequence_str.empty()) {
        CurcumaLogger::error("No polymer sequence provided. Use -sequence \"(A)10-B\"");
        return;
    }

    CurcumaLogger::info("Starting PolymerBuild");
    CurcumaLogger::param("sequence", sequence_str);

    std::vector<std::string> sequence = parseSequence(sequence_str);
    if (sequence.empty()) {
        CurcumaLogger::error("Parsed sequence is empty or invalid");
        return;
    }

    assemblePolymer(sequence);
}

void PolymerBuild::printHelp() const
{
    // Handled by ParameterRegistry
}

std::vector<std::string> PolymerBuild::parseSequence(const std::string& sequence)
{
    std::vector<std::string> result;
    // Simple parser for (A)n and A-B-C
    std::regex re("(\\(([A-Za-z0-9_]+)\\)([0-9]+))|([A-Za-z0-9_]+)");
    auto it = std::sregex_iterator(sequence.begin(), sequence.end(), re);
    auto end = std::sregex_iterator();

    for (; it != end; ++it) {
        std::smatch match = *it;
        if (match[1].matched) { // (NAME)NUMBER
            std::string name = match[2].str();
            int count = std::stoi(match[3].str());
            for (int i = 0; i < count; ++i) {
                result.push_back(name);
            }
        } else if (match[4].matched) { // NAME
            result.push_back(match[4].str());
        }
    }

    return result;
}

void PolymerBuild::assemblePolymer(const std::vector<std::string>& sequence)
{
    if (sequence.empty()) return;

    Molecule polymer;
    std::vector<std::pair<int, int>> interface_bonds;

    // Load first fragment
    std::string first_name = sequence[0];
    if (m_fragments.find(first_name) == m_fragments.end()) {
        CurcumaLogger::error("Fragment not found in map: " + first_name);
        return;
    }

    polymer.LoadMolecule(m_fragments[first_name]);
    polymer.setName("polymer_" + sequence[0]);

    for (size_t i = 1; i < sequence.size(); ++i) {
        std::string next_name = sequence[i];
        if (m_fragments.find(next_name) == m_fragments.end()) {
            CurcumaLogger::error("Fragment not found in map: " + next_name);
            continue;
        }

        CurcumaLogger::info(fmt::format("Adding fragment {}: {}", i, next_name));
        connectFragment(polymer, next_name, m_fragments[next_name], interface_bonds);

        // Update topology matrix for current polymer
        auto [dist, topo] = polymer.DistanceMatrix();
        // Claude Generated: interface_bonds are 0-indexed
        for (auto& bond : interface_bonds) {
            if (bond.first < topo.rows() && bond.second < topo.cols()) {
                topo(bond.first, bond.second) = 1;
                topo(bond.second, bond.first) = 1;
            }
        }
        polymer.setTopologyMatrix(topo);

        // Optional intermediate refinement
        if (m_config.get<bool>("optimize", true)) {
            CurcumaLogger::info("Running intermediate optimization...");
            json opt_ctrl;
            opt_ctrl["opt"]["method"] = m_config.get<std::string>("opt_method", "uff");
            opt_ctrl["opt"]["printOutput"] = false;
            opt_ctrl["opt"]["writeXYZ"] = false;

            CurcumaOpt opt(opt_ctrl, true);
            opt.addMolecule(polymer);
            opt.start();
            if (opt.Molecules()->size() > 0) {
                polymer = opt.Molecules()->back();
            }
        }

        if (m_config.get<bool>("dynamics", false)) {
            CurcumaLogger::info("Running intermediate dynamics...");
            json md_ctrl;
            md_ctrl["md"]["max_time"] = m_config.get<int>("md_steps", 1000);
            md_ctrl["md"]["method"] = m_config.get<std::string>("md_method", "uff");
            md_ctrl["md"]["verbosity"] = 0;

            SimpleMD md(md_ctrl, true);
            md.setMolecule(polymer);
            md.Initialise();
            md.start();
            auto unique = md.UniqueMolecules();
            if (!unique.empty()) {
                polymer = *unique.back();
            }
        }
    }

    applyCapping(polymer);

    // Final refresh of topology - Claude Generated: bonds are 0-indexed
    auto [dist, topo] = polymer.DistanceMatrix();
    for (auto& bond : interface_bonds) {
        if (bond.first < polymer.AtomCount() && bond.second < polymer.AtomCount()) {
            topo(bond.first, bond.second) = 1;
            topo(bond.second, bond.first) = 1;
        }
    }
    polymer.setTopologyMatrix(topo);

    polymer.writeXYZFile(polymer.Name() + "_final.xyz");
    CurcumaLogger::success("Polymer assembly completed: " + polymer.Name() + "_final.xyz");
}

void PolymerBuild::connectFragment(Molecule& current, const std::string& /*fragment_name*/,
                                   const std::string& fragment_file,
                                   std::vector<std::pair<int, int>>& interface_bonds)
{
    Molecule next;
    next.LoadMolecule(fragment_file);

    std::vector<int> current_xx;
    for (int i = 0; i < current.AtomCount(); ++i) {
        if (current.Atom(i).first == 0) current_xx.push_back(i);
    }

    std::vector<int> next_xx;
    for (int i = 0; i < next.AtomCount(); ++i) {
        if (next.Atom(i).first == 0) next_xx.push_back(i);
    }

    if (current_xx.empty()) {
        CurcumaLogger::error("Current structure has no Xx connection points");
        return;
    }
    if (next_xx.empty()) {
        CurcumaLogger::error("Next fragment has no Xx connection points");
        return;
    }

    int idx_out = current_xx.back(); // Use last Xx as outbound
    int idx_in = next_xx.front();   // Use first Xx as inbound

    // Claude Generated: Find closest non-H, non-Xx atom to each Xx
    // Xx atoms have covalent radius -1, so getConnectivtiy() doesn't work
    auto findClosestAtom = [](Molecule& mol, int xx_idx) -> int {
        Position xx_pos = mol.Atom(xx_idx).second;
        double min_dist = std::numeric_limits<double>::max();
        int closest = -1;

        for (int i = 0; i < mol.AtomCount(); ++i) {
            if (i == xx_idx) continue;
            int elem = mol.Atom(i).first;
            if (elem == 0 || elem == 1) continue; // Skip Xx and H

            double dist = (mol.Atom(i).second - xx_pos).norm();
            if (dist < min_dist) {
                min_dist = dist;
                closest = i;
            }
        }
        return closest;
    };

    int active_out = findClosestAtom(current, idx_out);
    int active_in = findClosestAtom(next, idx_in);

    if (active_out < 0) {
        CurcumaLogger::error("Could not find neighbor for Xx in current structure");
        return;
    }
    if (active_in < 0) {
        CurcumaLogger::error("Could not find neighbor for Xx in next fragment");
        return;
    }

    Position pos_out = current.Atom(idx_out).second;
    Position pos_active_out = current.Atom(active_out).second;
    Eigen::Vector3d v_out = (pos_out - pos_active_out).normalized();

    Position pos_in = next.Atom(idx_in).second;
    Position pos_active_in = next.Atom(active_in).second;
    Eigen::Vector3d v_in = (pos_in - pos_active_in).normalized();

    // Rotate next
    Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(v_in, -v_out);

    Geometry next_geom = next.Coords();
    for (int i = 0; i < next.AtomCount(); ++i) {
        Position p = next_geom.row(i);
        p -= pos_active_in;
        p = q * p;
        next_geom.row(i) = p;
    }

    double d = (Elements::CovalentRadius[current.Atom(active_out).first] +
                Elements::CovalentRadius[next.Atom(active_in).first]) * m_config.get<double>("bond_distance_scaling", 1.0);

    Position translation = pos_active_out + v_out * d;
    for (int i = 0; i < next.AtomCount(); ++i) {
        next_geom.row(i) += translation;
    }
    next.setGeometry(next_geom);

    int offset = current.AtomCount();
    for (int i = 0; i < next.AtomCount(); ++i) {
        current.addAtom(next.Atom(i));
    }

    // Record forced bond BEFORE removing atoms (0-indexed)
    interface_bonds.push_back({active_out, active_in + offset});

    // Remove Xx atoms
    std::vector<int> atoms_to_remove = {idx_out, idx_in + offset};
    current = current.AtomsRemoved(atoms_to_remove);

    // Claude Generated: Adjust interface_bonds for removed atom indices
    // Process in descending order for stable adjustment (larger indices first)
    int larger = std::max(idx_out, idx_in + offset);
    int smaller = std::min(idx_out, idx_in + offset);
    for (auto& bond : interface_bonds) {
        // Larger index removed first, smaller second
        if (bond.first > larger) bond.first--;
        if (bond.first > smaller) bond.first--;
        if (bond.second > larger) bond.second--;
        if (bond.second > smaller) bond.second--;
    }
}

void PolymerBuild::applyCapping(Molecule& mol)
{
    std::string cap_start = m_config.get<std::string>("cap_start", "H");
    std::string cap_end = m_config.get<std::string>("cap_end", "H");

    std::vector<int> xx_indices;
    for (int i = 0; i < mol.AtomCount(); ++i) {
        if (mol.Atom(i).first == 0) xx_indices.push_back(i);
    }

    for (size_t i = 0; i < xx_indices.size(); ++i) {
        int idx = xx_indices[i];
        std::string cap = (i == 0) ? cap_start : cap_end;
        if (cap == "none") continue;

        int element = Elements::String2Element(cap);
        if (element == 0) element = 1; // Default H

        Mol m = mol.getMolInfo();
        m.m_atoms[idx] = element;
        mol.LoadMolecule(m);
    }
}
