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

// Namespace-Aliase für Lesbarkeit
using namespace Eigen;
using MatrixXd = Eigen::MatrixXd;
using MatrixXf = Eigen::MatrixXf;
using MatrixXi = Eigen::MatrixXi;

// Farbskalen (gleich wie vorher)
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

    // Coolwarm - gut für symmetrische Daten um 0
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

// Template-basierter Image Writer für alle Eigen Matrix-Typen
class EigenImageWriter {
public:
    enum ColorMap { GRAYSCALE, JET, HOT, VIRIDIS, COOLWARM };

    // Template-Funktion für alle Eigen Matrix-Typen
    template<typename Derived>
    static bool saveMatrix(const Eigen::MatrixBase<Derived>& matrix,
                          const std::string& filename,
                          ColorMap colormap = GRAYSCALE,
                          int quality = 90,
                          bool symmetric_colormap = false) {

        int height = matrix.rows();
        int width = matrix.cols();

        // RGB-Daten erstellen
        std::vector<unsigned char> pixels(width * height * 3);

        // Min/Max finden - Template-kompatibel
        double minVal = matrix.minCoeff();
        double maxVal = matrix.maxCoeff();

        // Für symmetrische Farbskalen (z.B. bei Korrelationsmatrizen)
        if (symmetric_colormap) {
            double absMax = std::max(std::abs(minVal), std::abs(maxVal));
            minVal = -absMax;
            maxVal = absMax;
        }

        double range = maxVal - minVal;
        if (range == 0) range = 1;

        // Farbfunktion auswählen
        ColorMaps::RGB (*colorFunc)(double);
        switch (colormap) {
            case JET: colorFunc = ColorMaps::jet; break;
            case HOT: colorFunc = ColorMaps::hot; break;
            case VIRIDIS: colorFunc = ColorMaps::viridis; break;
            case COOLWARM: colorFunc = ColorMaps::coolwarm; break;
            default: colorFunc = ColorMaps::grayscale; break;
        }

        // Pixel füllen
        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                double value = static_cast<double>(matrix(i, j));
                double normalized = (value - minVal) / range;
                ColorMaps::RGB color = colorFunc(normalized);

                size_t idx = (i * width + j) * 3;
                pixels[idx] = color.r;
                pixels[idx + 1] = color.g;
                pixels[idx + 2] = color.b;
            }
        }

        // Dateiformat erkennen und speichern
        return saveToFile(filename, width, height, pixels, quality);
    }

    // Spezielle Funktion für Korrelationsmatrizen mit symmetrischer Farbskala
    template<typename Derived>
    static bool saveCorrelationMatrix(const Eigen::MatrixBase<Derived>& matrix,
                                    const std::string& filename,
                                    ColorMap colormap = COOLWARM) {
        return saveMatrix(matrix, filename, colormap, 90, true);
    }

    // Funktion für Matrix-Sequenzen (z.B. bei iterativen Algorithmen)
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

    // Hilfsfunktion für Matrix-Info
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

private:
    static bool saveToFile(const std::string& filename, int width, int height,
                          const std::vector<unsigned char>& pixels, int quality) {
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
};

// Beispiel-Funktionen für verschiedene Matrix-Typen
namespace Examples {

    // Erstelle eine zufällige Matrix
    MatrixXd createRandomMatrix(int size) {
        return MatrixXd::Random(size, size);
    }

    // Erstelle Korrelationsmatrix
    MatrixXd createCorrelationMatrix(int size) {
        MatrixXd data = MatrixXd::Random(size, size * 2);  // Mehr Samples als Variablen
        MatrixXd centered = data.rowwise() - data.colwise().mean();
        MatrixXd cov = (centered.adjoint() * centered) / double(data.rows() - 1);

        // Normalisierung zu Korrelation
        VectorXd std_dev = cov.diagonal().array().sqrt();
        MatrixXd corr = cov.array() / (std_dev * std_dev.transpose()).array();

        return corr;
    }

    // Mandelbrot mit Eigen
    MatrixXd createMandelbrot(int size, double zoom = 1.0, double centerX = -0.5, double centerY = 0) {
        MatrixXd result(size, size);

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                double x = centerX + (j - size/2.0) * 4.0 / (size * zoom);
                double y = centerY + (i - size/2.0) * 4.0 / (size * zoom);

                std::complex<double> c(x, y);
                std::complex<double> z(0, 0);
                int iteration = 0;

                while (std::abs(z) < 2 && iteration < 100) {
                    z = z * z + c;
                    iteration++;
                }

                result(i, j) = iteration / 100.0;
            }
        }

        return result;
    }

    // Heatmap einer Funktion
    MatrixXd createHeatmap(int size, std::function<double(double, double)> func) {
        MatrixXd result(size, size);

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                double x = (j - size/2.0) / (size/8.0);
                double y = (i - size/2.0) / (size/8.0);
                result(i, j) = func(x, y);
            }
        }

        return result;
    }
}

int main() {
    std::cout << "Eigen Matrix zu Bild Konvertierung\n";
    std::cout << "=================================\n\n";

    // 1. Einfache zufällige Matrix
    std::cout << "1. Zufällige Matrix...\n";
    MatrixXd randomMatrix = Examples::createRandomMatrix(300);
    EigenImageWriter::printMatrixInfo(randomMatrix, "Random Matrix");
    EigenImageWriter::saveMatrix(randomMatrix, "eigen_random.png", EigenImageWriter::JET);

    // 2. Korrelationsmatrix mit symmetrischer Farbskala
    std::cout << "2. Korrelationsmatrix...\n";
    MatrixXd corrMatrix = Examples::createCorrelationMatrix(100);
    EigenImageWriter::printMatrixInfo(corrMatrix, "Correlation Matrix");
    EigenImageWriter::saveCorrelationMatrix(corrMatrix, "eigen_correlation.png");

    // 3. Mandelbrot-Set
    std::cout << "3. Mandelbrot-Set...\n";
    MatrixXd mandelbrot = Examples::createMandelbrot(400, 0.8);
    EigenImageWriter::saveMatrix(mandelbrot, "eigen_mandelbrot.png", EigenImageWriter::HOT);

    // 4. Mathematische Funktionen als Heatmaps
    std::cout << "4. Mathematische Funktionen...\n";
    auto gaussianFunc = [](double x, double y) {
        return std::exp(-(x*x + y*y)/2);
    };
    auto sineWaveFunc = [](double x, double y) {
        return std::sin(x) * std::cos(y);
    };

    MatrixXd gaussian = Examples::createHeatmap(200, gaussianFunc);
    MatrixXd sineWave = Examples::createHeatmap(200, sineWaveFunc);

    EigenImageWriter::saveMatrix(gaussian, "eigen_gaussian.png", EigenImageWriter::VIRIDIS);
    EigenImageWriter::saveMatrix(sineWave, "eigen_sinewave.png", EigenImageWriter::COOLWARM, 90, true);

    // 5. Verschiedene Matrix-Typen demonstrieren
    std::cout << "5. Verschiedene Matrix-Typen...\n";

    // Float Matrix
    MatrixXf floatMatrix = MatrixXf::Random(150, 150);
    EigenImageWriter::saveMatrix(floatMatrix, "eigen_float.png", EigenImageWriter::GRAYSCALE);

    // Integer Matrix (z.B. diskrete Daten)
    MatrixXi intMatrix = (100 * MatrixXd::Random(100, 100)).cast<int>();
    EigenImageWriter::saveMatrix(intMatrix, "eigen_int.png", EigenImageWriter::JET);

    // 6. Eigenwerte/Eigenvektoren visualisieren
    std::cout << "6. Eigenvektoren Matrix...\n";
    MatrixXd symMatrix = randomMatrix + randomMatrix.transpose(); // Symmetrisch machen
    Eigen::SelfAdjointEigenSolver<MatrixXd> solver(symMatrix);
    MatrixXd eigenvectors = solver.eigenvectors();
    VectorXd eigenvalues = solver.eigenvalues();

    EigenImageWriter::saveMatrix(eigenvectors, "eigen_eigenvectors.png", EigenImageWriter::COOLWARM, 90, true);

    // Eigenvalues als 1D "Matrix" visualisieren
    MatrixXd eigenvalueMatrix(eigenvalues.size(), 20);
    for (int i = 0; i < eigenvalues.size(); ++i) {
        eigenvalueMatrix.row(i).setConstant(eigenvalues[i]);
    }
    EigenImageWriter::saveMatrix(eigenvalueMatrix, "eigen_eigenvalues.png", EigenImageWriter::HOT);

    std::cout << "\nFertig! Erstellte Dateien:\n";
    std::cout << "- eigen_random.png (Zufällige Matrix, Jet)\n";
    std::cout << "- eigen_correlation.png (Korrelation, Coolwarm)\n";
    std::cout << "- eigen_mandelbrot.png (Mandelbrot, Hot)\n";
    std::cout << "- eigen_gaussian.png (Gauß-Funktion, Viridis)\n";
    std::cout << "- eigen_sinewave.png (Sinus-Wellen, Coolwarm)\n";
    std::cout << "- eigen_float.png (Float Matrix)\n";
    std::cout << "- eigen_int.png (Integer Matrix)\n";
    std::cout << "- eigen_eigenvectors.png (Eigenvektoren)\n";
    std::cout << "- eigen_eigenvalues.png (Eigenwerte)\n";

    return 0;
}
