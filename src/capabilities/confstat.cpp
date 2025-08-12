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
#include "src/global_config.h"

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
    if (!Filename().empty()) {
        FileIterator file(Filename());
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
    m_stats.highest_energy = m_energies.back();
    m_stats.energy_span = (m_stats.highest_energy - m_stats.lowest_energy) * 2625.5;
    m_stats.total_states = std::accumulate(m_degeneracies.begin(),
        m_degeneracies.end(), 0);
    m_stats.unique_conformers = m_energies.size();
}

void ConfStat::printStatistics() const
{
    // Always show critical ConfStat results - these are primary outputs
    CurcumaLogger::success("═══ Conformer Statistics Analysis ═══");

    // Parameters - show at level 1 as they're essential for understanding results
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::param("Temperature", fmt::format("{:.2f} K", m_temp));
        CurcumaLogger::param("Energy cutoff", fmt::format("{:.1f} kJ/mol", m_cutoff));
        CurcumaLogger::param("Method", m_method);
    }

    CurcumaLogger::energy_abs(m_stats.lowest_energy, "Lowest energy");

    // Conformer table - critical data, always visible at level 1
    fmt::print("\n{:^5} | {:^15} | {:^15} | {:^12} | {:^10}\n",
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

    // Summary statistics - always visible as they are key results
    fmt::print("\nEnergy span: {:.2f} kJ/mol\n", m_stats.energy_span);
    fmt::print("Total states (including degeneracy): {}\n", m_stats.total_states);
    fmt::print("Free Energy: {:.2f} kJ/mol\n", m_stats.free_energy_contribution);
    fmt::print("Entropy (T*S): {:.2f} kJ/mol\n", m_stats.entropy_contribution);

    // Detailed distributions only at higher verbosity
    if (CurcumaLogger::get_verbosity() >= 2) {
        printEnergyDistribution();
        printPopulationDistribution();
    }
}

void ConfStat::setEnergiesWithDegeneracy(const std::vector<double>& energies,
    const std::vector<int>& degeneracies)
{
    if (energies.size() != degeneracies.size()) {
        throw std::invalid_argument(
            fmt::format("Mismatch in vector sizes: {} energies vs {} degeneracies",
                energies.size(), degeneracies.size()));
    }

    // Validate degeneracies
    if (std::any_of(degeneracies.begin(), degeneracies.end(),
            [](int deg) { return deg < 1; })) {
        throw std::invalid_argument("Degeneracy values must be positive integers");
    }

    m_energies = energies;
    m_degeneracies = degeneracies;
}

void ConfStat::printEnergyDistribution() const
{
    // Create energy histogram
    const int num_bins = 10;
    std::vector<int> histogram(num_bins, 0);
    double bin_width = m_stats.energy_span / num_bins;

    for (const auto& [diff, pop, cum, deg] : m_stats.conformer_data) {
        int bin = static_cast<int>(diff * 2625.5 / bin_width);
        if (bin >= 0 && bin < num_bins) {
            histogram[bin] += deg;
        }
    }

    CurcumaLogger::info("");
    CurcumaLogger::info("Energy Distribution:");
    for (int i = 0; i < num_bins; ++i) {
        double energy_start = i * bin_width;
        CurcumaLogger::info(fmt::format("{:5.1f} - {:5.1f} kJ/mol: {:3d} states | {:*<{}}",
            energy_start, energy_start + bin_width, histogram[i], "", histogram[i] / 2));
    }
}

void ConfStat::printPopulationDistribution() const
{
    CurcumaLogger::header("Population Distribution Analysis");

    // Calculate population statistics
    double total_pop = 0.0;
    int significant_conformers = 0; // Conformers with population > threshold

    // Population ranges for analysis (in %)
    const std::vector<double> ranges = { 50.0, 25.0, 10.0, 5.0, 1.0, 0.1 };
    std::vector<int> range_counts(ranges.size(), 0);

    for (const auto& [diff, pop, cum, deg] : m_stats.conformer_data) {
        total_pop += pop;

        if (pop >= m_print_threshold) {
            significant_conformers++;
        }

        // Count conformers in each population range
        for (size_t i = 0; i < ranges.size(); ++i) {
            if (pop >= ranges[i]) {
                range_counts[i]++;
            }
        }
    }

    // Print summary
    CurcumaLogger::param("Total conformers", std::to_string(m_stats.conformer_data.size()));
    CurcumaLogger::param("Significant conformers", fmt::format(">= {:.1f}%: {}", m_print_threshold, significant_conformers));

    // Print population ranges
    CurcumaLogger::info("");
    CurcumaLogger::info("Population ranges:");
    for (size_t i = 0; i < ranges.size(); ++i) {
        CurcumaLogger::info(fmt::format("Conformers with population ≥ {:5.1f}%: {:3d}",
            ranges[i], range_counts[i]));
    }

    // Print population distribution histogram
    CurcumaLogger::info("");
    CurcumaLogger::info("Population histogram:");
    const int num_bins = 10;
    std::vector<int> histogram(num_bins, 0);

    for (const auto& [diff, pop, cum, deg] : m_stats.conformer_data) {
        int bin = static_cast<int>(pop * num_bins / 100.0);
        if (bin >= num_bins)
            bin = num_bins - 1;
        if (bin < 0)
            bin = 0;
        histogram[bin]++;
    }

    // Print histogram
    for (int i = 0; i < num_bins; ++i) {
        double pop_start = i * 100.0 / num_bins;
        double pop_end = (i + 1) * 100.0 / num_bins;
        CurcumaLogger::info(fmt::format("{:4.1f}% - {:4.1f}%: {:3d} conformers | {:*<{}}",
            pop_start, pop_end, histogram[i], "", histogram[i]));
    }

    // Calculate and print statistical measures
    std::vector<double> populations;
    for (const auto& [diff, pop, cum, deg] : m_stats.conformer_data) {
        populations.push_back(pop);
    }

    double mean_pop = total_pop / m_stats.conformer_data.size();

    // Calculate standard deviation
    double variance = 0.0;
    for (double pop : populations) {
        variance += (pop - mean_pop) * (pop - mean_pop);
    }
    double std_dev = std::sqrt(variance / populations.size());

    CurcumaLogger::info("");
    CurcumaLogger::info("Statistical measures:");
    CurcumaLogger::param("Mean population", fmt::format("{:.2f}%", mean_pop));
    CurcumaLogger::param("Population std dev", fmt::format("{:.2f}%", std_dev));
}
