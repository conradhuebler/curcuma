/*
 * <Statistics of conformers >
 * Copyright (C) 2020 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *  Contribution by AI Copilot Claude 3.5 Sonnet
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

#pragma once

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "curcumamethod.h"
#include "src/core/energycalculator.h"
#include <fmt/format.h>

static json ConfStatJson{
    { "Cutoff", 10.0 }, // Energy cutoff in kJ/mol
    { "Temp", 298.15 }, // Temperature in Kelvin
    { "Threshold", 0.5 }, // Population threshold for printing
    { "Method", "none" }, // Energy calculation method
    { "Parameter", json::object() } // Additional parameters for energy calculation
};

struct Statistics {
    double lowest_energy{ 0.0 };
    double highest_energy{ 0.0 };
    double energy_span{ 0.0 };
    double free_energy_contribution{ 0.0 };
    double entropy_contribution{ 0.0 };
    int total_states{ 0 }; // Sum of all degeneracies
    int unique_conformers{ 0 }; // Number of unique energy levels
    std::vector<std::tuple<double, double, double, int>> conformer_data;
};

class ConfStat : public CurcumaMethod {
public:
    ConfStat(const json& controller = ConfStatJson, bool silent = true);

    // Set input filename
    void setFileName(const std::string& filename) { setFile(filename); }

    // Direct energy vector input
    void setEnergies(const std::vector<double>& energies) { m_energies = energies; }

    // Calculate population statistics
    void calculateStatistics();

    void start() override;

    void setEnergiesWithDegeneracy(const std::vector<double>& energies,
        const std::vector<int>& degeneracies);

    void printEnergyDistribution() const;
    void printPopulationDistribution() const;

    double getFirstPopulation() const
    {
        return !m_stats.conformer_data.empty() ? std::get<1>(m_stats.conformer_data[0]) : 0.0;
    }

    double getSecondPopulation() const
    {
        return m_stats.conformer_data.size() > 1 ? std::get<1>(m_stats.conformer_data[1]) : 0.0;
    }

    std::vector<double> getPopulations() const
    {
        std::vector<double> populations;
        for (const auto& [diff, pop, cum, deg] : m_stats.conformer_data) {
            populations.push_back(pop);
        }
        return populations;
    }

    double getTotalPopulation() const
    {
        return m_stats.conformer_data.empty() ? 0.0 : std::get<2>(m_stats.conformer_data.back());
    }

    int getSignificantConformerCount() const
    {
        return std::count_if(m_stats.conformer_data.begin(),
            m_stats.conformer_data.end(),
            [this](const auto& data) {
                return std::get<1>(data) >= m_print_threshold;
            });
    }

    double getEntropyContribution() const
    {
        return m_stats.entropy_contribution;
    }

    double getFreeEnergyContribution() const
    {
        return m_stats.free_energy_contribution;
    }

private:
    inline nlohmann::json WriteRestartInformation() override { return json(); }
    inline bool LoadRestartInformation() override { return true; }
    inline StringList MethodName() const override { return { std::string("ConfStat") }; }

    void ReadControlFile() override {};
    void LoadControlJson() override;
    void printStatistics() const;

    // Member variables with m_ prefix
    std::vector<double> m_energies;
    double m_temp{ 298.15 };
    double m_cutoff{ 10.0 }; // in kJ/mol
    double m_print_threshold{ 0.5 };
    std::string m_method{ "none" };
    json m_parameter;
    std::vector<int> m_degeneracies; // Degeneracy of each conformer

    Statistics m_stats;
};
