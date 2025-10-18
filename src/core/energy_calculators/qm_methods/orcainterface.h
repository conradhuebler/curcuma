//
// Created by gerd on 11.12.24.
//

#pragma once
#ifndef ORCAINTERFACE_H
#define ORCAINTERFACE_H

#include "src/core/molecule.h"
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

// Claude Generated 2025: ORCA Parameter Registry - replaces static ORCASettings JSON
BEGIN_PARAMETER_DEFINITION(orca)
    // Tight-Binding (TB) / SCF Parameters
    PARAM(accuracy, Int, 1, "Accuracy level for ORCA calculations.", "SCF", {"tb_acc"})
    PARAM(max_iterations, Int, 250, "Maximum number of SCF iterations.", "SCF", {"tb_max_iter"})
    PARAM(damping, Double, 0.4, "SCF damping parameter for convergence.", "SCF", {"tb_damping"})
    PARAM(electronic_temperature, Double, 9.500e-4, "Electronic temperature in Hartree for Fermi smearing.", "SCF", {"tb_temp"})
    PARAM(initial_guess, String, "SAD", "Initial guess method (SAD=Superposition of Atomic Densities).", "SCF", {"tb_guess"})
    PARAM(verbosity, Int, 0, "ORCA verbosity level (0=minimal, 1=normal, 2=verbose).", "Output", {"tb_verbose"})

    // Solvation Models
    PARAM(cpcm_solvent, String, "none", "CPCM solvent name (none, water, acetone, etc.).", "Solvation", {"cpcm_solv"})
    PARAM(cpcm_epsilon, Double, -1.0, "CPCM dielectric constant. -1 = auto-detect from name.", "Solvation", {"cpcm_eps"})
    PARAM(alpb_solvent, String, "none", "ALPB solvent name (none, water, acetone, etc.).", "Solvation", {"alpb_solv"})
    PARAM(alpb_epsilon, Double, -1.0, "ALPB dielectric constant. -1 = auto-detect from name.", "Solvation", {"alpb_eps"})
END_PARAMETER_DEFINITION

class OrcaInterface {
public:
    explicit OrcaInterface(const ConfigManager& config);
    explicit OrcaInterface();
    ~OrcaInterface();

    // Setzt die ORCA Eingabedaten
    void setInputFile(const std::string& inputFile);

    // Führt ORCA aus und wartet auf die Beendigung
    bool runOrca();

    // Liest die Ergebnisse aus der ORCA-Ausgabedatei
    void readOrcaJSON();

    // Create Output.JSON
    bool getOrcaJSON();

private:
    std::string inputFilePath;
    std::string outputFilePath;
    json OrcaJSON;
    mutable ConfigManager m_config;

    // Hilfsfunktionen
    bool createInputFile(const std::string& content);
    bool executeOrcaProcess();
};

#endif // ORCAINTERFACE_H