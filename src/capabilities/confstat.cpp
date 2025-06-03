/*
 * <Statistics of conformers >
 * Copyright (C) 2020 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Contribution by AI Copilot Claude 3.5 Sonnet
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

#include "src/core/energycalculator.h"
#include "src/core/fileiterator.h"

#include "confstat.h"

ConfStat::ConfStat(const json& controller, bool silent)
    : CurcumaMethod(ConfStatJson, controller, silent)
{
    UpdateController(controller);
}

void ConfStat::LoadControlJson()
{
    m_cutoff = Json2KeyWord<double>(m_defaults, "Cutoff");
    m_temp = Json2KeyWord<double>(m_defaults, "Temp");
    m_method = Json2KeyWord<std::string>(m_defaults, "Method");
    m_parameter = Json2KeyWord<json>(m_defaults, "Parameter");
    m_print_threshold = Json2KeyWord<double>(m_defaults, "Threshold");
}

void ConfStat::start()
{
    if (!m_filename.empty()) {
        FileIterator file(m_filename);
        m_energies.clear();

        while (!file.AtEnd()) {
            auto mol = std::make_unique<Molecule>(file.Next());
            double energy = mol->Energy();

            if (m_method != "none") {
                EnergyCalculator calculator(m_method, m_controller);
                calculator.setParameter(m_parameter);
                calculator.setMolecule(mol->getMolInfo());
                energy = calculator.CalculateEnergy();
            }

            m_energies.push_back(energy);
        }
    }

    calculateStatistics();
    printStatistics();
}

void ConfStat::calculateStatistics()
{
    if (m_energies.empty())
        return;

    std::sort(m_energies.begin(), m_energies.end());

    const double RT = R * m_temp;
    const double hartree_to_kjmol = 2625.5;

    // Calculate partition function Q
    double Q = 0.0;
    for (const double energy : m_energies) {
        double deltaE = (energy - m_energies[0]) * hartree_to_kjmol;
        if (deltaE > m_cutoff)
            continue;
        Q += std::exp(-deltaE * 1000 / RT);
    }

    // Calculate populations and statistics
    double S = 0.0;
    double cumulative = 0.0;
    m_stats.lowest_energy = m_energies[0];
    m_stats.conformer_data.clear();

    for (const double energy : m_energies) {
        double deltaE = (energy - m_energies[0]) * hartree_to_kjmol;
        if (deltaE > m_cutoff)
            continue;

        // Boltzmann population
        double p_i = std::exp(-deltaE * 1000 / RT) / Q;

        // Entropy contribution
        if (p_i > 0) { // Avoid log(0)
            S -= p_i * std::log(p_i);
        }

        cumulative += p_i * 100.0;
        m_stats.conformer_data.emplace_back(
            energy - m_energies[0], // deltaE in Hartree
            p_i * 100.0, // population in %
            cumulative // cumulative population in %
        );
    }

    // Final thermodynamic quantities
    m_stats.free_energy_contribution = -RT * std::log(Q) / 1000.0; // kJ/mol
    m_stats.entropy_contribution = RT * S / 1000.0; // kJ/mol
}

void ConfStat::printStatistics() const
{
    fmt::print("\nConformer Statistics Analysis\n");
    fmt::print("===========================\n");
    fmt::print("Parameters:\n");
    fmt::print("  Temperature: {:.2f} K\n", m_temp);
    fmt::print("  Energy cutoff: {:.1f} kJ/mol\n", m_cutoff);
    fmt::print("  Method: {}\n\n", m_method);

    fmt::print("Lowest energy: {:.5f} Eh\n\n", m_stats.lowest_energy);

    fmt::print("{:^5} | {:^15} | {:^15} | {:^12}\n",
        "Conf", "ΔE (kJ/mol)", "Population (%)", "Cumulative (%)");
    fmt::print("----------------------------------------------------------\n");

    int conf_number = 1;
    for (const auto& [diff, pop, cum] : m_stats.conformer_data) {
        if (pop >= m_print_threshold) {
            fmt::print("{:5d} | {:15.2f} | {:15.2f} | {:12.2f}\n",
                conf_number, diff * 2625.5, pop, cum);
        }
        conf_number++;
    }

    fmt::print("\nThermodynamic Contributions\n");
    fmt::print("-------------------------\n");
    fmt::print("Free Energy: {:.2f} kJ/mol\n", m_stats.free_energy_contribution);
    fmt::print("Entropy (T*S): {:.2f} kJ/mol\n", m_stats.entropy_contribution);
}
