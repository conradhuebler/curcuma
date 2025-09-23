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

using json = nlohmann::json;

PersistentDiagram::PersistentDiagram(const json& config, bool silent)
    : CurcumaMethod(RipserJson, config, silent)
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
    UpdateController(config);
    std::cout << config << std::endl;

    m_ratio = m_defaults["ripser_ratio"];
    m_xmax = m_defaults["ripser_xmax"];
    m_xmin = m_defaults["ripser_xmin"];
    m_ymax = m_defaults["ripser_ymax"];
    m_ymin = m_defaults["ripser_ymin"];
    m_bins = m_defaults["ripser_bins"];
    m_scaling = m_defaults["ripser_scaling"];
    m_std_x = m_defaults["ripser_stdx"];
    m_std_y = m_defaults["ripser_stdy"];
    m_dimension = m_defaults["ripser_dimension"];
    m_epsilon = m_defaults["ripser_epsilon"];
    m_min = m_defaults["ripser_min"];
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
    std::cout << "\n=== Persistent Diagram Configuration Parameters ===\n\n"
              << "Parameter        | Default | Description\n"
              << "----------------|---------|----------------------------------------------------\n"
              << "ripser_xmax     | " << std::setw(7) << m_xmax << " | Maximum X value for the output image\n"
              << "ripser_xmin     | " << std::setw(7) << m_xmin << " | Minimum X value for the output image\n"
              << "ripser_ymax     | " << std::setw(7) << m_ymax << " | Maximum Y value for the output image\n"
              << "ripser_ymin     | " << std::setw(7) << m_ymin << " | Minimum Y value for the output image\n"
              << "ripser_bins     | " << std::setw(7) << m_bins << " | Number of bins in each dimension of the output image\n"
              << "ripser_scaling  | " << std::setw(7) << m_scaling << " | Scaling factor for the Gaussian values in the image\n"
              << "ripser_stdx     | " << std::setw(7) << m_std_x << " | Standard deviation (width) for X Gaussian kernel\n"
              << "ripser_stdy     | " << std::setw(7) << m_std_y << " | Standard deviation (width) for Y Gaussian kernel\n"
              << "ripser_ratio    | " << std::setw(7) << m_ratio << " | Ratio of persistence to be included\n"
              << "ripser_dimension| " << std::setw(7) << m_dimension << " | Homology dimension (0: components, 1: loops, 2: voids)\n"
              << "ripser_epsilon  | " << std::setw(7) << m_epsilon << " | Distance threshold for the persistent homology computation\n"
              << "\nExample configuration in JSON:\n"
              << "{\n"
              << "  \"ripser_bins\": 20,\n"
              << "  \"ripser_dimension\": 1,\n"
              << "  \"ripser_scaling\": 0.2,\n"
              << "  \"ripser_stdx\": 15,\n"
              << "  \"ripser_stdy\": 15\n"
              << "}\n\n"
              << "How it works:\n"
              << "1. The distance matrix is processed to extract topological features\n"
              << "2. Features are represented as birth-death pairs in persistence diagrams\n"
              << "3. These pairs are converted to a 2D image using Gaussian kernels\n"
              << "4. The image can be compared across molecules for similarity\n"
              << std::endl;
}