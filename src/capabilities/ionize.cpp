/*
 * Ionize — neutralise a solvated system with counter-ions.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 - see ionize.h.
 */

#include "ionize.h"

#include "src/core/citation_registry.h"
#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"
#include "src/core/parameter_registry.h"
#include "src/tools/formats.h"
#include "src/tools/pbc_utils.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <map>
#include <random>

namespace curcuma {

namespace {

/// Atomic number of an element symbol, 0 when it is none.
int atomicNumber(const std::string& symbol)
{
    for (std::size_t z = 1; z < Elements::ElementAbbr.size(); ++z)
        if (Elements::ElementAbbr[z] == symbol)
            return static_cast<int>(z);
    return 0;
}

/// "H2O" -> {1: 2, 8: 1}. Empty map when the formula is not made of element symbols and counts.
std::map<int, int> parseFormula(const std::string& formula)
{
    std::map<int, int> composition;
    std::size_t i = 0;
    while (i < formula.size()) {
        if (!std::isupper(static_cast<unsigned char>(formula[i])))
            return {};
        std::string symbol(1, formula[i++]);
        if (i < formula.size() && std::islower(static_cast<unsigned char>(formula[i])))
            symbol += formula[i++];
        int count = 0;
        const bool has_digits = i < formula.size() && std::isdigit(static_cast<unsigned char>(formula[i]));
        while (i < formula.size() && std::isdigit(static_cast<unsigned char>(formula[i])))
            count = count * 10 + (formula[i++] - '0');
        const int z = atomicNumber(symbol);
        if (z == 0)
            return {};
        composition[z] += has_digits ? count : 1;
    }
    return composition;
}

} // namespace

IonizeResult ionize(const Molecule& input, const IonizeOptions& options)
{
    IonizeResult result;

    const int delta = options.target_charge - options.charge;
    result.ion = delta > 0 ? options.cation : options.anion;
    result.final_charge = options.charge;
    if (delta == 0) {
        result.molecule = input;
        result.ok = true;
        return result;
    }
    const int ion_z = atomicNumber(result.ion);
    if (ion_z == 0) {
        result.error = "unknown ion element \"" + result.ion + "\"";
        return result;
    }
    const std::map<int, int> solvent = parseFormula(options.solvent);
    if (solvent.empty()) {
        result.error = "cannot read the solvent formula \"" + options.solvent + "\"";
        return result;
    }

    // Solvent molecules are connected fragments with exactly the solvent composition.
    Molecule mol = input;
    const std::vector<std::vector<int>> fragments = mol.GetFragments();
    const std::vector<int> atoms = mol.Atoms();
    std::vector<int> solvent_fragments;
    std::vector<bool> is_solute(atoms.size(), true);
    for (std::size_t f = 0; f < fragments.size(); ++f) {
        std::map<int, int> composition;
        for (int a : fragments[f])
            ++composition[atoms[a]];
        if (composition == solvent) {
            solvent_fragments.push_back(static_cast<int>(f));
            for (int a : fragments[f])
                is_solute[a] = false;
        }
    }
    result.solvent_molecules = static_cast<int>(solvent_fragments.size());
    if (solvent_fragments.empty()) {
        result.error = "no " + options.solvent + " molecule found in the structure";
        return result;
    }

    // Minimum image convention only with a cell; the inverse is computed once.
    const bool periodic = mol.hasPBC();
    const Eigen::Matrix3d cell = periodic ? mol.getUnitCell() : Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d cell_inv = cell.inverse();
    const auto distance = [&](const Position& a, const Position& b) {
        const Eigen::Vector3d r = a - b;
        return periodic ? PBCUtils::applyMinimumImage(r, cell, cell_inv).norm() : r.norm();
    };

    // Each solvent molecule is represented by its heaviest atom (O in water), the ion goes there.
    struct Candidate {
        int fragment;
        int centre;
    };
    std::vector<Candidate> candidates;
    for (int f : solvent_fragments) {
        int centre = fragments[f].front();
        for (int a : fragments[f])
            if (atoms[a] > atoms[centre])
                centre = a;
        const Position p = mol.Atom(centre).second;
        double nearest = 1e300;
        for (std::size_t a = 0; a < atoms.size(); ++a)
            if (is_solute[a])
                nearest = std::min(nearest, distance(p, mol.Atom(static_cast<int>(a)).second));
        if (nearest >= options.min_solute_distance)
            candidates.push_back({ f, centre });
    }
    result.candidates = static_cast<int>(candidates.size());

    std::mt19937 rng(options.seed);
    std::shuffle(candidates.begin(), candidates.end(), rng);

    const int wanted = std::abs(delta);
    std::vector<Candidate> chosen;
    for (const Candidate& c : candidates) {
        if (static_cast<int>(chosen.size()) == wanted)
            break;
        const Position p = mol.Atom(c.centre).second;
        bool free = true;
        for (const Candidate& other : chosen)
            if (distance(p, mol.Atom(other.centre).second) < options.min_ion_distance) {
                free = false;
                break;
            }
        if (free)
            chosen.push_back(c);
    }
    if (static_cast<int>(chosen.size()) < wanted) {
        result.error = "only " + std::to_string(chosen.size()) + " of " + std::to_string(wanted)
            + " ions fit (" + std::to_string(result.candidates) + " of "
            + std::to_string(result.solvent_molecules) + " solvent molecules are at least "
            + std::to_string(options.min_solute_distance) + " A from the solute); lower "
            "min_solute_distance or min_ion_distance, or use a larger box";
        return result;
    }

    // Rebuild: everything but the replaced molecules in the input order, ions last.
    std::vector<bool> removed(atoms.size(), false);
    for (const Candidate& c : chosen)
        for (int a : fragments[c.fragment]) {
            removed[a] = true;
            result.replaced_atoms.push_back(a);
        }
    std::sort(result.replaced_atoms.begin(), result.replaced_atoms.end());

    Molecule out;
    for (std::size_t a = 0; a < atoms.size(); ++a)
        if (!removed[a])
            out.addPair(mol.Atom(static_cast<int>(a)));
    for (const Candidate& c : chosen)
        out.addPair({ ion_z, mol.Atom(c.centre).second });
    if (periodic)
        out.setUnitCell(cell, true);
    result.final_charge = options.charge + (delta > 0 ? wanted : -wanted);
    out.setCharge(result.final_charge);

    result.molecule = out;
    result.ions_added = wanted;
    result.ok = true;
    return result;
}

Ionize::Ionize(const json& controller, bool silent)
    : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("ionize"), controller, silent)
{
    UpdateController(controller);
}

void Ionize::setMolecule(const Molecule& molecule)
{
    m_molecule = molecule;
    m_have_molecule = true;
}

void Ionize::LoadControlJson()
{
    m_options.charge = Json2KeyWord<int>(m_defaults, "charge");
    m_options.target_charge = Json2KeyWord<int>(m_defaults, "target_charge");
    m_options.cation = Json2KeyWord<std::string>(m_defaults, "cation");
    m_options.anion = Json2KeyWord<std::string>(m_defaults, "anion");
    m_options.solvent = Json2KeyWord<std::string>(m_defaults, "solvent");
    m_options.min_solute_distance = Json2KeyWord<double>(m_defaults, "min_solute_distance");
    m_options.min_ion_distance = Json2KeyWord<double>(m_defaults, "min_ion_distance");
    m_options.seed = static_cast<unsigned>(Json2KeyWord<int>(m_defaults, "seed"));
    m_output = Json2KeyWord<std::string>(m_defaults, "output");
}

void Ionize::start()
{
    m_results = json();
    m_error.clear();

    Molecule input;
    if (m_have_molecule) {
        input = m_molecule;
    } else {
        if (Filename().empty() || !std::filesystem::exists(Filename())) {
            m_error = "no such file: " + Filename();
            return;
        }
        input = Files::LoadFile(Filename());
    }

    m_result = ionize(input, m_options);
    if (!m_result.ok) {
        m_error = m_result.error;
        CurcumaLogger::error(m_error);
        return;
    }

    std::string name = m_output;
    if (name.empty())
        name = (Basename().empty() ? std::string("system") : Basename()) + "_ionized.xyz";
    m_result.molecule.writeXYZFile(outputPath(name));

    m_results["output"] = outputPath(name);
    m_results["ion"] = m_result.ion;
    m_results["ions_added"] = m_result.ions_added;
    m_results["solvent_molecules"] = m_result.solvent_molecules;
    m_results["candidates"] = m_result.candidates;
    m_results["final_charge"] = m_result.final_charge;
    m_results["atoms"] = static_cast<int>(m_result.molecule.AtomCount());
    m_results["replaced_atoms"] = m_result.replaced_atoms;

    // Replacing solvent molecules by ions follows the approach of GROMACS genion (reference in the citation
    // database; TODO: check it against the paper).
    CitationRegistry::cite("gromacs2015");
    CurcumaLogger::success(std::to_string(m_result.ions_added) + " " + m_result.ion
        + " ion(s) replace solvent molecules; total charge " + std::to_string(m_result.final_charge)
        + ", " + std::to_string(m_result.molecule.AtomCount()) + " atoms: " + outputPath(name));
}

} // namespace curcuma
