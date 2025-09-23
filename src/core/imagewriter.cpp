/*
 * <Image writer for Eigen matrices - Implementation>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 * AI Supported by Claude 4 Sonnet
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

#include "imagewriter.hpp"
#include <iostream>

// STB implementation - this should only be included ONCE across the entire project
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// ColorMaps implementations
namespace ColorMaps {
    RGB grayscale(double value) {
        unsigned char v = (unsigned char)(255 * std::max(0.0, std::min(1.0, value)));
        return {v, v, v};
    }

    RGB jet(double value) {
        value = std::max(0.0, std::min(1.0, value));

        if (value < 0.25) {
            return {0, (unsigned char)(255 * value * 4), 255};
        } else if (value < 0.5) {
            return {0, 255, (unsigned char)(255 * (2 - value * 4))};
        } else if (value < 0.75) {
            return {(unsigned char)(255 * (value * 4 - 2)), 255, 0};
        } else {
            return {255, (unsigned char)(255 * (4 - value * 4)), 0};
        }
    }

    RGB hot(double value) {
        value = std::max(0.0, std::min(1.0, value));

        if (value < 1.0/3.0) {
            return {(unsigned char)(255 * value * 3), 0, 0};
        } else if (value < 2.0/3.0) {
            return {255, (unsigned char)(255 * (value * 3 - 1)), 0};
        } else {
            return {255, 255, (unsigned char)(255 * (value * 3 - 2))};
        }
    }

    RGB viridis(double value) {
        value = std::max(0.0, std::min(1.0, value));

        double r = 0.267 * value + 0.004;
        double g = 0.865 * value + 0.004;
        double b = 0.565 * (-value * value + 1.5 * value) + 0.015;

        return {
            (unsigned char)(255 * std::max(0.0, std::min(1.0, r))),
            (unsigned char)(255 * std::max(0.0, std::min(1.0, g))),
            (unsigned char)(255 * std::max(0.0, std::min(1.0, b)))
        };
    }

    RGB coolwarm(double value) {
        value = std::max(0.0, std::min(1.0, value));

        if (value < 0.5) {
            // Cool (blau zu weiß)
            double t = value * 2;
            return {
                (unsigned char)(255 * (0.23 + 0.77 * t)),
                (unsigned char)(255 * (0.29 + 0.71 * t)),
                (unsigned char)(255 * (0.75 + 0.25 * t))
            };
        } else {
            // Warm (weiß zu rot)
            double t = (value - 0.5) * 2;
            return {
                (unsigned char)(255 * (1.0 - 0.13 * t)),
                (unsigned char)(255 * (1.0 - 0.69 * t)),
                (unsigned char)(255 * (1.0 - 0.84 * t))
            };
        }
    }
}

// Implementation of static member functions
Eigen::MatrixXd EigenImageWriter::adaptiveMatrixScaling(const Eigen::MatrixXd& matrix,
    double temperatureParam,
    double dampingStrength,
    bool preserveStructure)
{

    Eigen::MatrixXd result = matrix;

    double globalMax = result.maxCoeff();
    double globalMin = result.minCoeff();
    double threshold = globalMin + 0.8 * (globalMax - globalMin);

    Eigen::MatrixXd dampingMask = Eigen::MatrixXd::Ones(result.rows(), result.cols());

    for (int i = 0; i < result.rows(); ++i) {
        for (int j = 0; j < result.cols(); ++j) {
            if (result(i, j) > threshold) {
                double distanceFromPeak = std::abs(result(i, j) - globalMax) / (globalMax - globalMin);
                dampingMask(i, j) = 0.2 + 0.8 * exp(-dampingStrength * distanceFromPeak * distanceFromPeak);
            }
        }
    }

    result = result.cwiseProduct(dampingMask);

    if (temperatureParam > 0) {
        Eigen::MatrixXd enhanced = applySoftmaxEnhancement(result, temperatureParam);

        if (preserveStructure) {
            result = 0.6 * result + 0.4 * enhanced;
        } else {
            result = enhanced;
        }
    }

    return result;
}

Eigen::MatrixXd EigenImageWriter::applySoftmaxEnhancement(const Eigen::MatrixXd& matrix, double temperature)
{
    Eigen::MatrixXd result = matrix;

    for (int i = 0; i < matrix.rows(); ++i) {
        Eigen::VectorXd row = matrix.row(i);

        std::vector<double> rowValues(row.data(), row.data() + row.size());
        std::sort(rowValues.begin(), rowValues.end());
        double median = rowValues[rowValues.size() / 2];

        double sumExp = 0.0;
        Eigen::VectorXd expValues(row.size());

        for (int j = 0; j < row.size(); ++j) {
            expValues(j) = exp(temperature * (row(j) - median));
            sumExp += expValues(j);
        }

        double originalScale = row.mean();
        for (int j = 0; j < row.size(); ++j) {
            result(i, j) = originalScale * expValues(j) / sumExp * row.size();
        }
    }

    return result;
}

ColorMaps::RGB EigenImageWriter::polymerColorMap(double value)
{
    value = std::max(0.0, std::min(1.0, value));

    if (value < 0.3) {
        double t = value / 0.3;
        return { (unsigned char)(50 * t), (unsigned char)(100 * t), (unsigned char)(200 + 55 * t) };
    } else if (value < 0.7) {
        double t = (value - 0.3) / 0.4;
        return { (unsigned char)(255 * t), (unsigned char)(200 + 55 * t), (unsigned char)(100 * (1 - t)) };
    } else {
        double t = (value - 0.7) / 0.3;
        return { (unsigned char)(255), (unsigned char)(200 * (1 - t)), (unsigned char)(50 * (1 - t)) };
    }
}

bool EigenImageWriter::saveToFile(const std::string& filename, int width, int height,
    const std::vector<unsigned char>& pixels, int quality)
{
    std::string ext = filename.substr(filename.find_last_of(".") + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if (ext == "png") {
        return stbi_write_png(filename.c_str(), width, height, 3, pixels.data(), width * 3);
    } else if (ext == "jpg" || ext == "jpeg") {
        return stbi_write_jpg(filename.c_str(), width, height, 3, pixels.data(), quality);
    } else if (ext == "bmp") {
        return stbi_write_bmp(filename.c_str(), width, height, 3, pixels.data());
    } else if (ext == "tga") {
        return stbi_write_tga(filename.c_str(), width, height, 3, pixels.data());
    }

    std::cerr << "Unsupported format: " << ext << std::endl;
    return false;
}