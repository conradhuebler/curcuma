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

#include <numeric>

#include "src/tools/general.h"

#include "persistentdiagram.h"

using json = nlohmann::json;

PersistentDiagram::PersistentDiagram(const json& config)
{
    json j = MergeJson(RipserJson, config);
    // std::cout << j << std::endl;
    m_ratio = j["ripser_ratio"];
    m_xmax = j["ripser_xmax"];
    m_xmin = j["ripser_xmin"];
    m_ymax = j["ripser_ymax"];
    m_ymin = j["ripser_ymin"];
    m_bins = j["ripser_bins"];
    m_scaling = j["ripser_scaling"];
    m_std_x = j["ripser_stdx"];
    m_std_y = j["ripser_stdy"];
    m_dimension = j["ripser_dimension"];
    m_epsilon = j["ripser_epsilon"];

    m_threshold = std::numeric_limits<float>::max();
}

dpairs PersistentDiagram::generatePairs()
{
    coefficient_t modulus = 2;
    std::vector<float> distance = m_compressed_lower_distance_matrix;
    compressed_lower_distance_matrix dist(std::move(m_compressed_lower_distance_matrix));

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
    coefficient_t modulus = 2;
    std::vector<float> distance = m_compressed_lower_distance_matrix;
    compressed_lower_distance_matrix dist(std::move(m_compressed_lower_distance_matrix));

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
