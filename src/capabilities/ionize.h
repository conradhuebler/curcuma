/*
 * Ionize — neutralise a solvated system with counter-ions.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 - A charged solute in explicit solvent (here: a protonated
 * polymer in water) needs counter-ions before a periodic simulation: the sum of
 * charges in a periodic cell must vanish. The ions replace solvent molecules, so the
 * number of atoms per volume stays what the solvent packing gave (the way GROMACS
 * `genion` works). Ions are monovalent: |charge - target_charge| of them are added,
 * anions for a positive system and cations for a negative one.
 *
 * Solvent molecules are the connected fragments whose composition equals the
 * `solvent` formula. Everything else is solute, including ions already present;
 * the given `charge` is the total charge of the input structure.
 */

#pragma once

#include "curcumamethod.h"

#include "src/core/molecule.h"
#include "src/core/parameter_macros.h"

#include <string>
#include <vector>

namespace curcuma {

/// What to add and where it may go. Defaults give a 0.15 nm-style spread suited to water boxes.
struct IonizeOptions {
    int charge = 0;                     ///< total charge of the input structure
    int target_charge = 0;              ///< total charge wanted afterwards
    std::string cation = "Na";          ///< element of the positive ion
    std::string anion = "Cl";           ///< element of the negative ion
    std::string solvent = "H2O";        ///< formula of the solvent molecule that is replaced
    double min_solute_distance = 6.0;   ///< Angstrom, nearest solute atom to a replaced solvent
    double min_ion_distance = 7.0;      ///< Angstrom, between two added ions
    unsigned seed = 1;                  ///< candidate order; the same seed gives the same ions
};

struct IonizeResult {
    bool ok = false;
    std::string error;
    Molecule molecule;                  ///< solvent kept, replaced molecules removed, ions appended last
    std::string ion;                    ///< element symbol of the added ions
    int ions_added = 0;
    int solvent_molecules = 0;          ///< solvent molecules found in the input
    int candidates = 0;                 ///< of those, far enough from the solute
    std::vector<int> replaced_atoms;    ///< input indices of the removed solvent atoms
    int final_charge = 0;
};

/**
 * @brief Replace solvent molecules by counter-ions.
 *
 * Distances use the minimum image convention when @p input has a unit cell. Fails
 * with an error text, not a partial result, when fewer ions fit than asked for.
 */
IonizeResult ionize(const Molecule& input, const IonizeOptions& options);

class Ionize : public CurcumaMethod {
public:
    explicit Ionize(const json& controller, bool silent = true);
    ~Ionize() = default;

    /// Structure in memory; used instead of the file given to setFile() when set.
    void setMolecule(const Molecule& molecule);

    /// Runs ionize(), writes `<output>` (default `<basename>_ionized.xyz`) and fills Results().
    void start() override;

    json Results() const override { return m_results; }
    const IonizeResult& result() const { return m_result; }
    std::string error() const { return m_error; }

private:
    nlohmann::json WriteRestartInformation() override { return json(); }
    bool LoadRestartInformation() override { return true; }
    StringList MethodName() const override { return { std::string("ionize") }; }
    void ReadControlFile() override {}
    void LoadControlJson() override;

    Molecule m_molecule;
    bool m_have_molecule = false;
    IonizeOptions m_options;
    std::string m_output;

    IonizeResult m_result;
    std::string m_error;
    json m_results;
};

} // namespace curcuma

namespace {
BEGIN_PARAMETER_DEFINITION(ionize)
MODULE_INFO("Neutralise a solvated system: replace solvent molecules by monovalent counter-ions.",
    "Assembly", { "ions", "neutralize", "genion" })

PARAM(charge, Int, 0, "Total charge of the input structure.", "Basic", {}, "tier=primary")
PARAM(target_charge, Int, 0, "Total charge after the ions are added.", "Basic", {})
PARAM(cation, String, "Na", "Element of the positive ion.", "Basic", {})
PARAM(anion, String, "Cl", "Element of the negative ion.", "Basic", {})
PARAM(solvent, String, "H2O", "Formula of the solvent molecule that is replaced.", "Basic", {})
PARAM(min_solute_distance, Double, 6.0, "Smallest distance from a replaced solvent molecule to any solute atom (Angstrom).", "Placement", {},
    "min=0")
PARAM(min_ion_distance, Double, 7.0, "Smallest distance between two added ions (Angstrom).", "Placement", {},
    "min=0")
PARAM(seed, Int, 1, "Seed for the order in which solvent molecules are tried.", "Placement", {},
    "min=0")
PARAM(output, String, "", "Output file; empty gives <input basename>_ionized.xyz.", "Output", {})
END_PARAMETER_DEFINITION
} // namespace
