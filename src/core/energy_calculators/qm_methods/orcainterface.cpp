//
// Created by gerd on 11.12.24.
//

#include "orcainterface.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

OrcaInterface::OrcaInterface()
{
    inputFilePath = "orca.inp";
    outputFilePath = "orca.out";
}

OrcaInterface::~OrcaInterface()
{
    // Optionale Bereinigungsoperationen
}

void OrcaInterface::setInputFile(const std::string& inputFile)
{
    inputFilePath = inputFile;
}

bool OrcaInterface::createInputFile(const std::string& content)
{
    std::ofstream outFile(inputFilePath);
    if (!outFile) {
        std::cerr << "Fehler beim Erstellen der Eingabedatei!" << std::endl;
        return false;
    }
    outFile << content;
    outFile.close();
    return true;
}

bool OrcaInterface::executeOrcaProcess()
{
    // Hier rufen wir das ORCA-Programm über einen Systemaufruf auf
    std::stringstream command;
    command << "orca " << inputFilePath << " > " << outputFilePath;
    int result = std::system(command.str().c_str());

    // Überprüfen, ob der ORCA-Prozess erfolgreich ausgeführt wurde
    return (result == 0);
}

bool OrcaInterface::runOrca()
{
    // Starten Sie den ORCA-Prozess und warten Sie auf das Ergebnis
    std::cout << "Starte ORCA..." << std::endl;
    if (executeOrcaProcess()) {
        std::cout << "ORCA abgeschlossen!" << std::endl;
        return true;
    } else {
        std::cerr << "Fehler beim Ausführen von ORCA!" << std::endl;
        return false;
    }
}

void OrcaInterface::readOrcaJSON()
{
    // Liest die Ergebnisse aus der ORCA-Ausgabedatei
    std::ifstream property(inputFilePath + ".property.json");
    property >> OrcaJSON;
}

bool OrcaInterface::getOrcaJSON()
{
    // Hier rufen wir das ORCA_2JSON-Programm über einen Systemaufruf auf
    std::stringstream command;
    command << "orca_2json " << inputFilePath << " -property >> " << outputFilePath;
    const int result = std::system(command.str().c_str());

    // Überprüfen, ob der ORCA-Prozess erfolgreich ausgeführt wurde
    return (result == 0);
}