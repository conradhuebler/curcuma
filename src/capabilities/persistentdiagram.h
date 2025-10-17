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

#pragma once
// #include "ripser.h"

#include "src/capabilities/curcumamethod.h"
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

#include "json.hpp"

#include <vector>

#include <Eigen/Dense>

struct triple {
    double start = 0;
    double end = 0;
    double scaling = 1;
};

typedef std::pair<double, double> dpair;

typedef std::vector<dpair> dpairs;
typedef std::vector<triple> triples;

using json = nlohmann::json;

/* Claude Generated 2025: PersistentDiagram Parameter Registry - replaces static RipserJson */
BEGIN_PARAMETER_DEFINITION(ripser)
    PARAM(x_max, Double, 4.0, "Max x-value for persistence diagram.", "Diagram", {"ripser_xmax"})
    PARAM(x_min, Double, 0.0, "Min x-value for persistence diagram.", "Diagram", {"ripser_xmin"})
    PARAM(y_max, Double, 4.0, "Max y-value for persistence diagram.", "Diagram", {"ripser_ymax"})
    PARAM(y_min, Double, 0.0, "Min y-value for persistence diagram.", "Diagram", {"ripser_ymin"})
    PARAM(bins, Int, 10, "Number of bins for persistence image.", "Image", {"ripser_bins"})
    PARAM(scaling, Double, 0.1, "Scaling factor for persistence image.", "Image", {"ripser_scaling"})
    PARAM(std_x, Double, 10.0, "Standard deviation for x-axis in persistence image.", "Image", {"ripser_stdx"})
    PARAM(std_y, Double, 10.0, "Standard deviation for y-axis in persistence image.", "Image", {"ripser_stdy"})
    PARAM(ratio, Double, 1.0, "Ratio for ripser calculation.", "Algorithm", {"ripser_ratio"})
    PARAM(dimension, Int, 2, "Dimension for ripser calculation.", "Algorithm", {"ripser_dimension"})
    PARAM(epsilon, Double, 0.4, "Epsilon for ripser calculation.", "Algorithm", {"ripser_epsilon"})
    PARAM(min, Double, 0.0, "Minimum threshold.", "Algorithm", {"ripser_min"})
END_PARAMETER_DEFINITION

class PersistentDiagram : public CurcumaMethod {
    // This class generates persistent diagrams from distance matrices using ripser.
    // It can generate pairs and triples of points, and create images based on these pairs.
    // The parameters for the persistent diagram can be set through the constructor or setter methods.  {
public:
    /**
     * @brief Constructor with JSON configuration (backward compatible)
     * Claude Generated: Phase 4 - ConfigManager Migration
     */
    PersistentDiagram(const json& config, bool silent = true);

    /**
     * @brief Constructor with ConfigManager configuration (new, preferred)
     * Claude Generated: Phase 4 - Native ConfigManager support
     */
    PersistentDiagram(const ConfigManager& config, bool silent = true);

    void start() override
    {
        // This method is intentionally left empty as the main action is handled in generatePairs and generateTriples.
    }
    nlohmann::json WriteRestartInformation() override
    {
        return json(); // Return an empty JSON object as there is no restart information to write.
        // This method is intentionally left empty as there is no restart information to write.
    }
    bool LoadRestartInformation() override
    {
        // This method is intentionally left empty as there is no restart information to load.
        return true;
    }
    StringList MethodName() const override
    {
        return StringList({ "PD" });
    }

    void ReadControlFile() override
    {
        // This method is intentionally left empty as there is no control file to read.
    }

    void LoadControlJson() override
    {
        // This method is intentionally left empty as there is no control JSON to load.
    }

    inline void setDistanceMatrix(const std::vector<float>& vector)
    {
        m_compressed_lower_distance_matrix = vector;
    }

    inline void setENScaling(const std::vector<double>& en_scaling)
    {
        m_en_scaling = en_scaling;
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
    triples generateTriples();

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
    Eigen::MatrixXd generateImage(const triples& pairs);

    void setScaling(double scaling) { m_scaling = scaling; }
    void setBins(int bins) { m_bins = bins; }
    void setStdX(double std_x) { m_std_x = std_x; }
    void setStdY(double std_y) { m_std_x = std_y; }

    void printHelp() const override;

private:
    int m_dimension = 2;
    double m_threshold;
    double m_ratio = 1;
    double m_xmax = 4, m_xmin = 0, m_ymax = 4, m_ymin = 0;
    int m_bins = 10.0;
    double m_scaling = 0.1;
    double m_std_x = 10;
    double m_std_y = 10;
    double m_epsilon = 0.4;
    double m_min = 0;
    std::vector<float> m_compressed_lower_distance_matrix;
    std::vector<double> m_en_scaling;
};
