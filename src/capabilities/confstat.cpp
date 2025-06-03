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

    // If no degeneracies provided, assume all are 1
    if (m_degeneracies.empty()) {
        m_degeneracies.resize(m_energies.size(), 1);
    }

    // Sort energies and maintain degeneracy correlation
    std::vector<size_t> indices(m_energies.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
        [this](size_t a, size_t b) { return m_energies[a] < m_energies[b]; });

    std::vector<double> sorted_energies;
    std::vector<int> sorted_degeneracies;
    for (size_t idx : indices) {
        sorted_energies.push_back(m_energies[idx]);
        sorted_degeneracies.push_back(m_degeneracies[idx]);
    }
    m_energies = sorted_energies;
    m_degeneracies = sorted_degeneracies;

    const double RT = R * m_temp;
    const double hartree_to_kjmol = 2625.5;

    // Calculate partition function Q including degeneracy
    double Q = 0.0;
    for (size_t i = 0; i < m_energies.size(); ++i) {
        double deltaE = (m_energies[i] - m_energies[0]) * hartree_to_kjmol;
        if (deltaE > m_cutoff)
            continue;
        Q += m_degeneracies[i] * std::exp(-deltaE * 1000 / RT);
    }

    // Calculate populations and statistics
    double S = 0.0;
    double cumulative = 0.0;
    m_stats.lowest_energy = m_energies[0];
    m_stats.conformer_data.clear();

    for (size_t i = 0; i < m_energies.size(); ++i) {
        double deltaE = (m_energies[i] - m_energies[0]) * hartree_to_kjmol;
        if (deltaE > m_cutoff)
            continue;

        // Boltzmann population including degeneracy
        double p_i = m_degeneracies[i] * std::exp(-deltaE * 1000 / RT) / Q;

        // Entropy contribution
        if (p_i > 0) {
            S -= p_i * std::log(p_i / m_degeneracies[i]); // Modified for degeneracy
        }

        cumulative += p_i * 100.0;
        m_stats.conformer_data.emplace_back(
            m_energies[i] - m_energies[0], // deltaE in Hartree
            p_i * 100.0, // population in %
            cumulative, // cumulative population in %
            m_degeneracies[i] // degeneracy
        );
    }

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

    fmt::print("{:^5} | {:^15} | {:^15} | {:^12} | {:^10}\n",
        "Conf", "ΔE (kJ/mol)", "Population (%)", "Cumulative (%)", "Degeneracy");
    fmt::print("--------------------------------------------------------------\n");

    int conf_number = 1;
    for (const auto& [diff, pop, cum, deg] : m_stats.conformer_data) {
        if (pop >= m_print_threshold) {
            fmt::print("{:5d} | {:15.2f} | {:15.2f} | {:12.2f} | {:10d}\n",
                conf_number, diff * 2625.5, pop, cum, deg);
        }
        conf_number++;
    }

    fmt::print("\nThermodynamic Contributions\n");
    fmt::print("-------------------------\n");
    fmt::print("Free Energy: {:.2f} kJ/mol\n", m_stats.free_energy_contribution);
    fmt::print("Entropy (T*S): {:.2f} kJ/mol\n", m_stats.entropy_contribution);
}
