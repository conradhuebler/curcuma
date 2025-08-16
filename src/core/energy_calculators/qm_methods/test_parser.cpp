#include "GTOIntegrals.hpp"
#include "STOIntegrals.hpp"
#include "basissetparser.hpp"
#include <iostream>

// Beispielprogramm zur Demonstration des Basissatz-Parsers
int main(int argc, char** argv)
{
    // Überprüfen der Befehlszeilenargumente
    if (argc < 2) {
        std::cerr << "Verwendung: " << argv[0] << " <basissatz-datei>" << std::endl;
        return 1;
    }

    std::string basisSetFile = argv[1];

    try {
        // Basis-Set einlesen
        auto basisSetMap = BasisSetParser::parseBasisSetFile(basisSetFile);

        // VSIP-Werte für Extended Hückel-Berechnung setzen
        BasisSetParser::setVSIPValues(basisSetMap, "H", -13.6, 0.0);
        BasisSetParser::setVSIPValues(basisSetMap, "C", -19.44, -10.67);
        BasisSetParser::setVSIPValues(basisSetMap, "N", -26.0, -13.4);
        BasisSetParser::setVSIPValues(basisSetMap, "O", -32.3, -14.8);

        // Informationen zum Basissatz ausgeben
        BasisSetParser::printBasisSetInfo(basisSetMap);

        // Methan-Molekül als Beispiel erstellen
        std::cout << "\n=== Erstellung eines Methan-Moleküls (CH4) ===" << std::endl;

        // Tetraeder-Geometrie für Methan
        double factor = 1.09 / std::sqrt(3.0); // CH-Bindungslänge / sqrt(3)

        // Kohlenstoff im Ursprung
        std::vector<STO::Orbital> c_orbitals;
        if (basisSetMap.find("C") != basisSetMap.end()) {
            c_orbitals = BasisSetParser::createSTOFromGTOBasis(
                basisSetMap.at("C"), 0.0, 0.0, 0.0, 0);

            std::cout << "Erstellt: " << c_orbitals.size()
                      << " STO-Orbitale für Kohlenstoff." << std::endl;
        } else {
            std::cerr << "Element C nicht im Basissatz gefunden!" << std::endl;
            return 1;
        }

        // Vier Wasserstoff-Atome in Tetraeder-Anordnung
        std::vector<std::vector<STO::Orbital>> h_orbitals;
        if (basisSetMap.find("H") != basisSetMap.end()) {
            std::vector<std::vector<double>> h_positions = {
                { factor, factor, factor },
                { -factor, -factor, factor },
                { factor, -factor, -factor },
                { -factor, factor, -factor }
            };

            for (int i = 0; i < 4; ++i) {
                auto orbitals = BasisSetParser::createSTOFromGTOBasis(
                    basisSetMap.at("H"),
                    h_positions[i][0], h_positions[i][1], h_positions[i][2],
                    i + 1);

                h_orbitals.push_back(orbitals);
                std::cout << "Erstellt: " << orbitals.size()
                          << " STO-Orbitale für Wasserstoff " << i + 1 << "." << std::endl;
            }
        } else {
            std::cerr << "Element H nicht im Basissatz gefunden!" << std::endl;
            return 1;
        }

        // Alle Orbitale sammeln
        std::vector<STO::Orbital> all_orbitals = c_orbitals;
        for (const auto& h_orbital_set : h_orbitals) {
            all_orbitals.insert(all_orbitals.end(), h_orbital_set.begin(), h_orbital_set.end());
        }

        // Überlappungsmatrix berechnen
        std::cout << "\n=== Berechnung der Überlappungsmatrix ===" << std::endl;
        int n = all_orbitals.size();

        std::cout << "Anzahl der Orbitale: " << n << std::endl;
        std::cout << "Überlappungsmatrix:" << std::endl;

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                double overlap = STO::calculateOverlap(all_orbitals[i], all_orbitals[j]);
                std::cout << "S[" << i << "," << j << "] = " << overlap;
                if (j == i) {
                    std::cout << std::endl;
                } else {
                    std::cout << ", ";
                }
            }
        }

        // Beispiel für GTO-Basis
        std::cout << "\n=== GTO-Basis für Methan ===" << std::endl;

        std::vector<GTO::Orbital> c_gto_orbitals;
        if (basisSetMap.find("C") != basisSetMap.end()) {
            c_gto_orbitals = BasisSetParser::createGTOFromBasis(
                basisSetMap.at("C"), 0.0, 0.0, 0.0, 0);

            std::cout << "Erstellt: " << c_gto_orbitals.size()
                      << " GTO-Orbitale für Kohlenstoff." << std::endl;

            // Details zu GTO-Orbitalen anzeigen
            for (const auto& orbital : c_gto_orbitals) {
                std::string type;
                switch (orbital.type) {
                case GTO::OrbitalType::S:
                    type = "S";
                    break;
                case GTO::OrbitalType::PX:
                    type = "PX";
                    break;
                case GTO::OrbitalType::PY:
                    type = "PY";
                    break;
                case GTO::OrbitalType::PZ:
                    type = "PZ";
                    break;
                default:
                    type = "andere";
                    break;
                }

                std::cout << "  Orbital vom Typ " << type
                          << " mit " << orbital.exponents.size() << " Exponenten"
                          << std::endl;
            }
        }

        // GTO-Überlappung für C 2s-2s (Selbstüberlappung) berechnen
        if (!c_gto_orbitals.empty() && c_gto_orbitals[0].type == GTO::OrbitalType::S) {
            std::cout << "\nGTO-Selbstüberlappung für C 2s-2s: "
                      << GTO::calculateOverlap(c_gto_orbitals[0], c_gto_orbitals[0])
                      << std::endl;
        }

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Fehler: " << e.what() << std::endl;
        return 1;
    }
}
