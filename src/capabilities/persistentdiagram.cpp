/*
 * <Generate Persisent Diagrams from Distance Matrices using ripser>
 * Copyright (C) 2021 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <numeric>

#include "persistentdiagram.h"

PersistentDiagram::PersistentDiagram()
{
    m_threshold = std::numeric_limits<float>::max();
}

dpairs PersistentDiagram::generatePairs()
{
    coefficient_t modulus = 2;

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

Eigen::MatrixXd PersistentDiagram::generateImage(const dpairs& pairs)
{
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(m_bins + 1, m_bins + 1);

    double scaling = 0.1;
    double i_bins = 1 / (m_bins);
    double i_c = 10;
    for (int i = 0; i < matrix.cols(); ++i) // iteration over y values
    {
        double cur_y = m_ymax - (m_ymax - m_ymin) * i_bins * (i);
        for (int j = 0; j < matrix.rows(); ++j) // iteration over x values
        {
            double cur_x = (m_xmax - m_xmin) * i_bins * (j);

            for (const auto& pair : pairs) {
                double x = pair.first;
                double y = pair.second;
                matrix(i, j) += scaling * exp(-((cur_x - x) * (cur_x - x) * i_c)) * exp((-(cur_y - y) * (cur_y - y) * i_c));
            }
        }
    }
    return matrix;
}
