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

#pragma once
//#include "ripser.h"

#include <vector>

#include <Eigen/Dense>

typedef std::pair<double, double> dpair;

typedef std::vector<dpair> dpairs;

class PersistentDiagram {
public:
    PersistentDiagram();

    inline void setDistanceMatrix(const std::vector<float>& vector)
    {
        m_compressed_lower_distance_matrix = vector;
    }

    inline void setDimension(int dim)
    {
        m_dimension = dim;
    }

    inline void setRatio(double ratio)
    {
        m_ratio = ratio;
    }

    inline void setThreshold(double threshold)
    {
        m_threshold = threshold;
    }

    dpairs generatePairs();

    inline void setXRange(double xmin, double xmax)
    {
        m_xmin = xmin;
        m_xmax = xmax;
    }

    inline void setYRange(double ymin, double ymax)
    {
        m_ymin = ymin;
        m_ymax = ymax;
    }

    Eigen::MatrixXd generateImage(const dpairs& pairs);

private:
    int m_dimension = 2;
    double m_threshold;
    double m_ratio = 1;
    double m_xmax, m_xmin, m_ymax, m_ymin;
    double m_bins = 10.0;
    std::vector<float> m_compressed_lower_distance_matrix;
};
