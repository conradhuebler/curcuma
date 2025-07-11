/*
 * <Image writer fpr Eigen matrices>
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


#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

// Eigen einbinden
#include <Eigen/Dense>

// stb_image_write.h
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

namespace ColorMaps {
    struct RGB { unsigned char r, g, b; };

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

class EigenImageWriter {
public:
    enum ColorMap { GRAYSCALE,
        JET,
        HOT,
        VIRIDIS,
        COOLWARM,
        POLYMER };
    enum PostProcessing { NONE,
        ADAPTIVE,
        RING_FOCUSED };

private:
    static Eigen::MatrixXd applySoftmaxEnhancement(const Eigen::MatrixXd& matrix, double temperature);
    static Eigen::MatrixXd adaptiveMatrixScaling(const Eigen::MatrixXd& matrix,
        double temperatureParam = 2.0,
        double dampingStrength = 1.5,
        bool preserveStructure = true);
    static ColorMaps::RGB polymerColorMap(double value);
    static bool saveToFile(const std::string& filename, int width, int height,
        const std::vector<unsigned char>& pixels, int quality);

public:
    template <typename Derived>
    static bool saveMatrix(const Eigen::MatrixBase<Derived>& matrix,
        const std::string& filename,
        ColorMap colormap = GRAYSCALE,
        int quality = 90,
        bool symmetric_colormap = false,
        int width = -1,
        int height = -1,
        PostProcessing post_processing = NONE,
        double temperature = 2.0,
        double damping = 1.5,
        bool preserve_structure = true)
    {

        Eigen::MatrixXd processed_matrix = matrix.template cast<double>();

        switch (post_processing) {
        case ADAPTIVE:
            processed_matrix = adaptiveMatrixScaling(processed_matrix, temperature, damping, preserve_structure);
            break;
        case RING_FOCUSED:
            processed_matrix = adaptiveMatrixScaling(processed_matrix, 3.0, 2.0, false);
            break;
        case NONE:
        default:
            break;
        }

        int img_height = (height > 0) ? height : processed_matrix.rows();
        int img_width = (width > 0) ? width : processed_matrix.cols();

        std::vector<unsigned char> pixels(img_width * img_height * 3);

        double minVal = processed_matrix.minCoeff();
        double maxVal = processed_matrix.maxCoeff();

        if (symmetric_colormap) {
            double absMax = std::max(std::abs(minVal), std::abs(maxVal));
            minVal = -absMax;
            maxVal = absMax;
        }

        double range = maxVal - minVal;
        if (range == 0) range = 1;

        ColorMaps::RGB (*colorFunc)(double);
        switch (colormap) {
            case JET: colorFunc = ColorMaps::jet; break;
            case HOT: colorFunc = ColorMaps::hot; break;
            case VIRIDIS: colorFunc = ColorMaps::viridis; break;
            case COOLWARM: colorFunc = ColorMaps::coolwarm; break;
            case POLYMER:
                colorFunc = polymerColorMap;
                break;
            default: colorFunc = ColorMaps::grayscale; break;
        }

        for (int i = 0; i < img_height; ++i) {
            for (int j = 0; j < img_width; ++j) {
                double row_ratio = static_cast<double>(i) / img_height;
                double col_ratio = static_cast<double>(j) / img_width;
                int mat_row = static_cast<int>(row_ratio * processed_matrix.rows());
                int mat_col = static_cast<int>(col_ratio * processed_matrix.cols());

                double value = static_cast<double>(processed_matrix(mat_row, mat_col));
                double normalized = (value - minVal) / range;
                ColorMaps::RGB color = colorFunc(normalized);

                size_t idx = (i * img_width + j) * 3;
                pixels[idx] = color.r;
                pixels[idx + 1] = color.g;
                pixels[idx + 2] = color.b;
            }
        }

        return saveToFile(filename, img_width, img_height, pixels, quality);
    }

    template<typename Derived>
    static bool saveCorrelationMatrix(const Eigen::MatrixBase<Derived>& matrix,
                                    const std::string& filename,
                                    ColorMap colormap = COOLWARM) {
        return saveMatrix(matrix, filename, colormap, 90, true);
    }

    template<typename Derived>
    static bool saveMatrixSequence(const std::vector<Eigen::Matrix<typename Derived::Scalar,
                                                                   Derived::RowsAtCompileTime,
                                                                   Derived::ColsAtCompileTime>>& matrices,
                                  const std::string& baseFilename,
                                  ColorMap colormap = GRAYSCALE) {
        for (size_t i = 0; i < matrices.size(); ++i) {
            std::string filename = baseFilename + "_" +
                                 std::to_string(i) + ".png";
            if (!saveMatrix(matrices[i], filename, colormap)) {
                return false;
            }
        }
        return true;
    }

    template<typename Derived>
    static void printMatrixInfo(const Eigen::MatrixBase<Derived>& matrix,
                               const std::string& name = "Matrix") {
        std::cout << name << " Info:\n";
        std::cout << "  Größe: " << matrix.rows() << "x" << matrix.cols() << "\n";
        std::cout << "  Min: " << matrix.minCoeff() << "\n";
        std::cout << "  Max: " << matrix.maxCoeff() << "\n";
        std::cout << "  Mean: " << matrix.mean() << "\n";
        std::cout << "  Norm: " << matrix.norm() << "\n\n";
    }
};

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
