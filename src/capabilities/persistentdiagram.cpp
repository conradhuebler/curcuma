/*
 * <Generate Persisent Diagrams from Distance Matrices using ripser>
 * Copyright (C) 2021 - 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "ripser.h"

#include "json.hpp"

#include <iomanip>
#include <numeric>

#include "persistentdiagram.h"
#include "src/tools/general.h"
#include "src/core/parameter_registry.h"

using json = nlohmann::json;

PersistentDiagram::PersistentDiagram(const json& config, bool silent)
    : PersistentDiagram(ConfigManager("ripser", config), silent)
{
}

PersistentDiagram::PersistentDiagram(const ConfigManager& config, bool silent)
    : CurcumaMethod(json{}, config.exportConfig(), silent)
    , m_ratio(1.0)
    , m_xmax(4.0)
    , m_xmin(0.0)
    , m_ymax(4.0)
    , m_ymin(0.0)
    , m_bins(10)
    , m_scaling(0.1)
    , m_std_x(10.0)
    , m_std_y(10.0)
    , m_dimension(2)
    , m_epsilon(0.4)
{
    UpdateController(config.exportConfig());

    // Claude Generated 2025: Updated to use new canonical parameter names from ParameterRegistry
    m_ratio = m_defaults.value("ratio", 1.0);
    m_xmax = m_defaults.value("x_max", 4.0);
    m_xmin = m_defaults.value("x_min", 0.0);
    m_ymax = m_defaults.value("y_max", 4.0);
    m_ymin = m_defaults.value("y_min", 0.0);
    m_bins = m_defaults.value("bins", 10);
    m_scaling = m_defaults.value("scaling", 0.1);
    m_std_x = m_defaults.value("std_x", 10.0);
    m_std_y = m_defaults.value("std_y", 10.0);
    m_dimension = m_defaults.value("dimension", 2);
    m_epsilon = m_defaults.value("epsilon", 0.4);
    m_min = m_defaults.value("min", 0.0);
    m_threshold = std::numeric_limits<float>::max();
    checkHelp();
}

dpairs PersistentDiagram::generatePairs()
{
    // Claude Generated: Input validation and memory safety
    if (m_compressed_lower_distance_matrix.empty()) {
        return dpairs(); // Return empty pairs for invalid input
    }

    coefficient_t modulus = 2;
    std::vector<float> distance = m_compressed_lower_distance_matrix; // Keep copy for later use
    compressed_lower_distance_matrix dist(std::move(m_compressed_lower_distance_matrix)); // Use move for ripser

    value_t enclosing_radius = std::numeric_limits<value_t>::infinity();
    if (m_threshold == std::numeric_limits<value_t>::max()) {
        for (size_t i = 0; i < dist.size(); ++i) {
            value_t r_i = -std::numeric_limits<value_t>::infinity();
            for (size_t j = 0; j < dist.size(); ++j)
                r_i = std::max(r_i, dist(i, j));
            enclosing_radius = std::min(enclosing_radius, r_i);
        }
    }

    auto result = ripser<compressed_lower_distance_matrix>(std::move(dist), m_dimension, enclosing_radius,
        m_ratio, modulus)
                      .compute_barcodes();

    dpairs final;
    for (const auto& a : result)
        for (const auto& b : a.second) {
            final.push_back(b);
        }

    return final;
}

triples PersistentDiagram::generateTriples()
{
    // Claude Generated: Input validation and memory safety for generateTriples
    if (m_compressed_lower_distance_matrix.empty()) {
        return triples(); // Return empty triples for invalid input
    }

    if (m_en_scaling.empty()) {
        return triples(); // Return empty triples if EN scaling is missing
    }

    coefficient_t modulus = 2;
    std::vector<float> distance = m_compressed_lower_distance_matrix; // Keep copy for later use
    compressed_lower_distance_matrix dist(std::move(m_compressed_lower_distance_matrix)); // Use move for ripser

    value_t enclosing_radius = std::numeric_limits<value_t>::infinity();
    if (m_threshold == std::numeric_limits<value_t>::max()) {
        for (size_t i = 0; i < dist.size(); ++i) {
            value_t r_i = -std::numeric_limits<value_t>::infinity();
            for (size_t j = 0; j < dist.size(); ++j)
                r_i = std::max(r_i, dist(i, j));
            enclosing_radius = std::min(enclosing_radius, r_i);
        }
    }

    auto result = ripser<compressed_lower_distance_matrix>(std::move(dist), m_dimension, enclosing_radius,
        m_ratio, modulus)
                      .compute_barcodes();

    triples final;
    for (const auto& a : result)
        if (a.first != 0) {
            for (const auto& b : a.second) {
                triple t;
                t.start = b.first;
                t.end = b.second;
                final.push_back(t);
            }
        } else {
            int counter = 0;
            std::vector<int> indices;
            for (const auto& b : a.second)

            {
                bool found = false;
                for (int i = 0; i < distance.size(); ++i) {
                    // Claude Generated: Add bounds checking for m_en_scaling access
                    if (i >= static_cast<int>(m_en_scaling.size())) {
                        break; // Prevent out-of-bounds access to m_en_scaling
                    }

                    if (std::abs(distance[i] - b.second) < 1e-6 && std::find(indices.begin(), indices.end(), i) == indices.end()) {
                        //    std::cout << i << " " << distance[i] << " " << m_en_scaling[i] << " " << counter++ << std::endl;
                        indices.push_back(i);
                        triple t;
                        t.start = b.first;
                        t.end = b.second;
                        t.scaling = m_en_scaling[i] + m_epsilon;
                        final.push_back(t);

                        found = true;
                        break;
                    }
                }
                if (found == false) {
                    triple t;
                    t.start = b.first;
                    t.end = b.second;
                    final.push_back(t);
                }
            }
        }

    return final;
}

Eigen::MatrixXd PersistentDiagram::generateImage(const dpairs& pairs)
{
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(m_bins, m_bins);

    double i_bins = 1 / double(m_bins);
    for (int i = 0; i < matrix.cols(); ++i) // iteration over y values
    {
        double cur_y = m_ymax - (m_ymax - m_ymin) * i_bins * (i);
        for (int j = 0; j < matrix.rows(); ++j) // iteration over x values
        {
            double cur_x = (m_xmax - m_xmin) * i_bins * (j);
            for (const auto& pair : pairs) {
                double x = pair.first;
                double y = pair.second;
                if (x < m_min || y < m_min)
                    continue; // skip points outside the defined range
                matrix(i, j) += m_scaling * exp(-((cur_x - x) * (cur_x - x) * m_std_x)) * exp((-(cur_y - y) * (cur_y - y) * m_std_y));
            }
        }
    }
    return matrix;
}

Eigen::MatrixXd PersistentDiagram::generateImage(const triples& pairs)
{
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(m_bins, m_bins);

    double i_bins = 1 / double(m_bins);
    for (int i = 0; i < matrix.cols(); ++i) // iteration over y values
    {
        double cur_y = m_ymax - (m_ymax - m_ymin) * i_bins * (i);
        for (int j = 0; j < matrix.rows(); ++j) // iteration over x values
        {
            double cur_x = (m_xmax - m_xmin) * i_bins * (j);
            for (const auto& triple : pairs) {
                double x = triple.start;
                double y = triple.end;
                double std = triple.scaling;
                matrix(i, j) += m_scaling * exp(-((cur_x - x) * (cur_x - x) * m_std_x / std)) * exp((-(cur_y - y) * (cur_y - y) * m_std_y / std));
            }
        }
    }
    return matrix;
}

void PersistentDiagram::printHelp() const
{
    ParameterRegistry::getInstance().printHelp("persistentdiagram");
}