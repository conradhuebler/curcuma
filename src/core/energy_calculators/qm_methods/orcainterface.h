//
// Created by gerd on 11.12.24.
//

#pragma once
#ifndef ORCAINTERFACE_H
#define ORCAINTERFACE_H

#include "src/core/molecule.h"

static json ORCASettings{
    { "tb_acc", 1 },
    { "tb_max_iter", 250 },
    { "tb_damping", 0.4 },
    { "tb_temp", 9.500e-4 },
    { "tb_verbose", 0 },
    { "tb_guess", "SAD" },
    { "cpcm_solv", "none" },
    { "alpb_solv", "none" },
    { "cpcm_eps", -1 },
    { "alpb_eps", -1 }
};

class OrcaInterface {
public:
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

    // Hilfsfunktionen
    bool createInputFile(const std::string& content);
    bool executeOrcaProcess();
};

#endif // ORCAINTERFACE_H