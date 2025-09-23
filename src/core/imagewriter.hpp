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

#pragma once 
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

// Eigen einbinden
#include <Eigen/Dense>

// Forward declaration - implementation in imagewriter.cpp
extern "C" {
int stbi_write_png(char const *filename, int w, int h, int comp, const void *data, int stride_in_bytes);
int stbi_write_jpg(char const *filename, int w, int h, int comp, const void *data, int quality);
int stbi_write_bmp(char const *filename, int w, int h, int comp, const void *data);
int stbi_write_tga(char const *filename, int w, int h, int comp, const void *data);
}

namespace ColorMaps {
    struct RGB { unsigned char r, g, b; };

    // Function declarations - implementations in imagewriter.cpp
    RGB grayscale(double value);
    RGB jet(double value);
    RGB hot(double value);
    RGB viridis(double value);
    RGB coolwarm(double value);
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

