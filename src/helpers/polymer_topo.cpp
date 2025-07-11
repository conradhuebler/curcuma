#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace Eigen;

// Struktur für Persistence Intervals
struct PersistenceInterval {
    double birth;
    double death;
    int dimension;

    double persistence() const { return death - birth; }
    double midpoint() const { return (birth + death) / 2.0; }
    double ratio() const { return death / birth; }
};

// Matrix-basierte Analyse-Klasse
class MatrixBasedTopologyAnalyzer {
private:
    MatrixXd persistenceMatrix;
    std::vector<PersistenceInterval> intervals;

public:
    MatrixBasedTopologyAnalyzer(const MatrixXd& matrix)
        : persistenceMatrix(matrix)
    {
        analyzeMatrix();
    }

    // Hauptanalyse der Matrix
    void analyzeMatrix()
    {
        detectPersistencePeaks();
        classifyIntervals();
    }

    // Erkenne charakteristische Peaks in der Persistence-Matrix
    void detectPersistencePeaks()
    {
        intervals.clear();

        // Annahme: Matrix hat Filtration-Werte als Zeilen/Spalten
        // Oder: Matrix repräsentiert Birth-Death Diagramm

        int rows = persistenceMatrix.rows();
        int cols = persistenceMatrix.cols();

        // Methode 1: Interpretiere als Birth-Death Diagramm
        if (rows == cols) {
            analyzeBirthDeathMatrix();
        } else {
            // Methode 2: Interpretiere als Filtration-Persistence Profil
            analyzeFiltrationProfile();
        }
    }

private:
    void analyzeBirthDeathMatrix()
    {
        int size = persistenceMatrix.rows();

        // Obere Dreiecksmatrix analysieren (Birth < Death)
        for (int i = 0; i < size; ++i) {
            for (int j = i + 1; j < size; ++j) {
                double value = persistenceMatrix(i, j);

                // Threshold für significance
                if (value > 1e-6) {
                    double birth = static_cast<double>(i) / size * 10.0; // Skalierung
                    double death = static_cast<double>(j) / size * 10.0;

                    // Dimension basierend auf Charakteristika schätzen
                    int dim = estimateDimensionFromPattern(birth, death, value);
                    intervals.push_back({ birth, death, dim });
                }
            }
        }
    }

    void analyzeFiltrationProfile()
    {
        int rows = persistenceMatrix.rows();
        int cols = persistenceMatrix.cols();

        // Interpretiere als [Filtration x Features] Matrix
        for (int col = 0; col < cols; ++col) {
            VectorXd profile = persistenceMatrix.col(col);

            // Finde lokale Maxima im Profil
            auto peaks = findLocalMaxima(profile);

            for (const auto& peak : peaks) {
                double birth = peak.first;
                double death = peak.first + peak.second; // width als proxy für death
                int dim = col % 3; // Rotation durch Dimensionen
                intervals.push_back({ birth, death, dim });
            }
        }
    }

    int estimateDimensionFromPattern(double birth, double death, double value)
    {
        double persistence = death - birth;
        double ratio = death / birth;

        // Heuristik basierend auf statistischen Mustern:
        // - Kurze, intensive Peaks → wahrscheinlich H0 (Connectivity)
        // - Mittlere Persistence mit moderater Intensität → H1 (Loops)
        // - Lange, schwache Persistence → H2 (Voids)

        if (persistence < 1.0 && value > 0.5) {
            return 0; // H0: Connected components
        } else if (persistence >= 1.0 && persistence < 3.0) {
            return 1; // H1: Loops (interessant für Polymere!)
        } else {
            return 2; // H2: Voids
        }
    }

    std::vector<std::pair<double, double>> findLocalMaxima(const VectorXd& profile)
    {
        std::vector<std::pair<double, double>> peaks;
        int size = profile.size();

        for (int i = 1; i < size - 1; ++i) {
            if (profile(i) > profile(i - 1) && profile(i) > profile(i + 1)) {
                double position = static_cast<double>(i) / size * 10.0;
                double width = estimatePeakWidth(profile, i);
                peaks.push_back({ position, width });
            }
        }

        return peaks;
    }

    double estimatePeakWidth(const VectorXd& profile, int peakIndex)
    {
        double peakValue = profile(peakIndex);
        double threshold = peakValue * 0.5; // Half-maximum

        int leftBound = peakIndex;
        int rightBound = peakIndex;

        // Finde linke Grenze
        while (leftBound > 0 && profile(leftBound) > threshold) {
            leftBound--;
        }

        // Finde rechte Grenze
        while (rightBound < profile.size() - 1 && profile(rightBound) > threshold) {
            rightBound++;
        }

        return static_cast<double>(rightBound - leftBound) / profile.size() * 10.0;
    }

    void classifyIntervals()
    {
        // Sortiere nach Dimension und Persistence
        std::sort(intervals.begin(), intervals.end(),
            [](const PersistenceInterval& a, const PersistenceInterval& b) {
                if (a.dimension != b.dimension)
                    return a.dimension < b.dimension;
                return a.persistence() > b.persistence();
            });
    }

public:
    // Matrix-basierte adaptive Skalierung
    MatrixXd adaptiveMatrixScaling(double temperatureParam = 2.0,
        double dampingStrength = 1.5,
        bool preserveStructure = true)
    {

        MatrixXd result = persistenceMatrix;

        // 1. Identifiziere dominante Strukturen (Peaks)
        double globalMax = result.maxCoeff();
        double globalMin = result.minCoeff();
        double threshold = globalMin + 0.8 * (globalMax - globalMin);

        // 2. Erstelle Dämpfungsmaske für dominante Features
        MatrixXd dampingMask = MatrixXd::Ones(result.rows(), result.cols());

        for (int i = 0; i < result.rows(); ++i) {
            for (int j = 0; j < result.cols(); ++j) {
                if (result(i, j) > threshold) {
                    // Gaußdämpfung um dominante Peaks
                    double distanceFromPeak = std::abs(result(i, j) - globalMax) / (globalMax - globalMin);
                    dampingMask(i, j) = 0.2 + 0.8 * exp(-dampingStrength * distanceFromPeak * distanceFromPeak);
                }
            }
        }

        // 3. Anwenden der Dämpfung
        result = result.cwiseProduct(dampingMask);

        // 4. Softmax-artige Enhancement für schwächere Signale
        if (temperatureParam > 0) {
            MatrixXd enhanced = applySoftmaxEnhancement(result, temperatureParam);

            if (preserveStructure) {
                // Gewichtete Kombination: Original-Struktur + Enhancement
                result = 0.6 * result + 0.4 * enhanced;
            } else {
                result = enhanced;
            }
        }

        return result;
    }

private:
    MatrixXd applySoftmaxEnhancement(const MatrixXd& matrix, double temperature)
    {
        MatrixXd result = matrix;

        // Berechne row-wise oder spaltenweise Softmax
        for (int i = 0; i < matrix.rows(); ++i) {
            VectorXd row = matrix.row(i);

            // Finde lokalen Kontext (Median als Referenz)
            std::vector<double> rowValues(row.data(), row.data() + row.size());
            std::sort(rowValues.begin(), rowValues.end());
            double median = rowValues[rowValues.size() / 2];

            // Softmax mit lokaler Referenz
            double sumExp = 0.0;
            VectorXd expValues(row.size());

            for (int j = 0; j < row.size(); ++j) {
                expValues(j) = exp(temperature * (row(j) - median));
                sumExp += expValues(j);
            }

            // Normalisierung mit Erhalt der Originalskala
            double originalScale = row.mean();
            for (int j = 0; j < row.size(); ++j) {
                result(i, j) = originalScale * expValues(j) / sumExp * row.size();
            }
        }

        return result;
    }

public:
    // Polymer-spezifische Analyse
    struct PolymerTopologyMetrics {
        double ringDensity; // H1 Feature-Dichte
        double chainConnectivity; // H0 Charakteristika
        double cavityStructure; // H2 Features
        double topologicalComplexity; // Gesamtkomplexität
        std::vector<double> ringDistribution; // Größenverteilung der Ringe
    };

    PolymerTopologyMetrics analyzePolymerTopology()
    {
        PolymerTopologyMetrics metrics;

        // Zähle Features nach Dimension
        std::map<int, int> dimensionCounts;
        std::vector<double> h1Persistences;

        for (const auto& interval : intervals) {
            dimensionCounts[interval.dimension]++;

            if (interval.dimension == 1) {
                h1Persistences.push_back(interval.persistence());
            }
        }

        // Ring-Dichte (H1 Features)
        metrics.ringDensity = dimensionCounts[1] / static_cast<double>(intervals.size());

        // Chain-Connectivity (H0 Features mit kurzer Persistence)
        int shortH0 = 0;
        for (const auto& interval : intervals) {
            if (interval.dimension == 0 && interval.persistence() < 2.0) {
                shortH0++;
            }
        }
        metrics.chainConnectivity = shortH0 / static_cast<double>(dimensionCounts[0] + 1);

        // Cavity-Struktur (H2)
        metrics.cavityStructure = dimensionCounts[2] / static_cast<double>(intervals.size() + 1);

        // Topological Complexity (Shannon-Entropie der Dimensionsverteilung)
        metrics.topologicalComplexity = calculateTopologicalComplexity(dimensionCounts);

        // Ring-Größenverteilung
        metrics.ringDistribution = h1Persistences;
        std::sort(metrics.ringDistribution.begin(), metrics.ringDistribution.end());

        return metrics;
    }

private:
    double calculateTopologicalComplexity(const std::map<int, int>& counts)
    {
        int total = 0;
        for (const auto& pair : counts) {
            total += pair.second;
        }

        if (total == 0)
            return 0.0;

        double entropy = 0.0;
        for (const auto& pair : counts) {
            double p = static_cast<double>(pair.second) / total;
            if (p > 0) {
                entropy -= p * log2(p);
            }
        }

        return entropy;
    }

public:
    // Visualisierung für Polymer-Analyse
    bool visualizePolymerTopology(const std::string& filenameBase,
        const PolymerTopologyMetrics& metrics)
    {

        // 1. Original Matrix
        visualizeMatrix(persistenceMatrix, filenameBase + "_original.png");

        // 2. Adaptiv skalierte Matrix
        MatrixXd adaptiveMatrix = adaptiveMatrixScaling(2.0, 1.5, true);
        visualizeMatrix(adaptiveMatrix, filenameBase + "_adaptive.png");

        // 3. Ring-fokussierte Ansicht
        MatrixXd ringFocused = adaptiveMatrixScaling(3.0, 2.0, false);
        visualizeMatrix(ringFocused, filenameBase + "_ring_focused.png");

        // 4. Barcode-Darstellung
        visualizeBarcodes(filenameBase + "_barcodes.png");

        // 5. Metriken-Dashboard
        createMetricsDashboard(metrics, filenameBase + "_metrics.png");

        return true;
    }

private:
    bool visualizeMatrix(const MatrixXd& matrix, const std::string& filename,
        bool useAdaptiveColormap = true)
    {

        int height = matrix.rows();
        int width = matrix.cols();
        std::vector<unsigned char> pixels(width * height * 3);

        double minVal = matrix.minCoeff();
        double maxVal = matrix.maxCoeff();
        double range = maxVal - minVal;
        if (range == 0)
            range = 1;

        // Adaptive Farbzuordnung für Polymer-Strukturen
        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                double normalized = (matrix(i, j) - minVal) / range;

                RGB color;
                if (useAdaptiveColormap) {
                    color = polymerColorMap(normalized);
                } else {
                    color = viridisColor(normalized);
                }

                size_t idx = (i * width + j) * 3;
                pixels[idx] = color.r;
                pixels[idx + 1] = color.g;
                pixels[idx + 2] = color.b;
            }
        }

        return stbi_write_png(filename.c_str(), width, height, 3,
            pixels.data(), width * 3);
    }

    struct RGB {
        unsigned char r, g, b;
    };

    RGB polymerColorMap(double value)
    {
        // Spezielle Farbkarte für Polymer-Topologie
        value = std::max(0.0, std::min(1.0, value));

        if (value < 0.3) {
            // Schwache Signale: Blautöne (potenzielle Ring-Regionen)
            double t = value / 0.3;
            return { (unsigned char)(50 * t), (unsigned char)(100 * t), (unsigned char)(200 + 55 * t) };
        } else if (value < 0.7) {
            // Mittlere Signale: Grün-Gelb (interessante Topologie)
            double t = (value - 0.3) / 0.4;
            return { (unsigned char)(255 * t), (unsigned char)(200 + 55 * t), (unsigned char)(100 * (1 - t)) };
        } else {
            // Starke Signale: Rot (dominante Strukturen)
            double t = (value - 0.7) / 0.3;
            return { (unsigned char)(255), (unsigned char)(200 * (1 - t)), (unsigned char)(50 * (1 - t)) };
        }
    }

    RGB viridisColor(double value)
    {
        // Standard Viridis für Vergleich
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

    bool visualizeBarcodes(const std::string& filename)
    {
        if (intervals.empty())
            return false;

        int width = 800;
        int height = 600;
        std::vector<unsigned char> pixels(width * height * 3, 255);

        // Sortiere nach Dimension und Persistence
        auto sortedIntervals = intervals;
        std::sort(sortedIntervals.begin(), sortedIntervals.end(),
            [](const PersistenceInterval& a, const PersistenceInterval& b) {
                if (a.dimension != b.dimension)
                    return a.dimension < b.dimension;
                return a.persistence() > b.persistence();
            });

        // Zeichne Barcodes
        double maxDeath = 10.0; // Anpassbar
        int barHeight = 4;
        int currentY = 20;

        for (const auto& interval : sortedIntervals) {
            int x1 = (int)(interval.birth / maxDeath * width);
            int x2 = (int)(interval.death / maxDeath * width);

            RGB color;
            switch (interval.dimension) {
            case 0:
                color = { 200, 100, 100 };
                break; // H0: Rot
            case 1:
                color = { 100, 200, 100 };
                break; // H1: Grün (Ringe!)
            case 2:
                color = { 100, 100, 200 };
                break; // H2: Blau
            default:
                color = { 150, 150, 150 };
                break;
            }

            // Zeichne Balken
            for (int y = currentY; y < currentY + barHeight && y < height; ++y) {
                for (int x = x1; x <= x2 && x < width; ++x) {
                    if (x >= 0 && y >= 0) {
                        int idx = (y * width + x) * 3;
                        pixels[idx] = color.r;
                        pixels[idx + 1] = color.g;
                        pixels[idx + 2] = color.b;
                    }
                }
            }

            currentY += barHeight + 2;
            if (currentY >= height - barHeight)
                break;
        }

        return stbi_write_png(filename.c_str(), width, height, 3,
            pixels.data(), width * 3);
    }

    bool createMetricsDashboard(const PolymerTopologyMetrics& metrics,
        const std::string& filename)
    {
        int width = 600;
        int height = 400;
        std::vector<unsigned char> pixels(width * height * 3, 255);

        // Einfaches Text-Dashboard (erweitert mit echter Text-Renderung)
        // Hier: Balkendiagramm für Metriken

        std::vector<double> values = {
            metrics.ringDensity,
            metrics.chainConnectivity,
            metrics.cavityStructure,
            metrics.topologicalComplexity / 3.0 // Normalisiert
        };

        std::vector<RGB> colors = {
            { 100, 200, 100 }, // Grün für Ring-Dichte
            { 200, 100, 100 }, // Rot für Chain-Connectivity
            { 100, 100, 200 }, // Blau für Cavity-Struktur
            { 150, 150, 100 } // Gelb für Komplexität
        };

        int barWidth = 80;
        int barSpacing = 120;
        int maxBarHeight = 300;

        for (size_t i = 0; i < values.size(); ++i) {
            int barHeight = (int)(values[i] * maxBarHeight);
            int startX = 50 + i * barSpacing;
            int startY = height - 50 - barHeight;

            for (int y = startY; y < height - 50; ++y) {
                for (int x = startX; x < startX + barWidth; ++x) {
                    if (x >= 0 && x < width && y >= 0 && y < height) {
                        int idx = (y * width + x) * 3;
                        pixels[idx] = colors[i].r;
                        pixels[idx + 1] = colors[i].g;
                        pixels[idx + 2] = colors[i].b;
                    }
                }
            }
        }

        return stbi_write_png(filename.c_str(), width, height, 3,
            pixels.data(), width * 3);
    }

public:
    // Getter für externe Analyse
    const std::vector<PersistenceInterval>& getIntervals() const { return intervals; }
    const MatrixXd& getMatrix() const { return persistenceMatrix; }
};

// Hilfsfunktion: Generiere Beispiel-Matrizen für verschiedene Polymer-Typen
MatrixXd generatePolymerExample(const std::string& type, int size = 200)
{
    MatrixXd matrix = MatrixXd::Zero(size, size);

    if (type == "linear_chain") {
        // Lineare Kette: Starke diagonale Struktur
        for (int i = 0; i < size - 1; ++i) {
            matrix(i, i + 1) = 1.0 - (double)i / size * 0.3; // Abnehmende Stärke
        }

    } else if (type == "ring_polymer") {
        // Ring-Polymer: Zyklische Struktur mit charakteristischen H1-Features
        for (int i = 0; i < size; ++i) {
            for (int j = i; j < size; ++j) {
                double distance = std::min(j - i, size - (j - i));
                if (distance < size / 10) {
                    matrix(i, j) = exp(-distance * distance / 20.0);
                    matrix(j, i) = matrix(i, j);
                }
            }
        }

        // Füge Ring-charakteristische Features hinzu
        int ringSize = size / 4;
        for (int ring = 0; ring < 3; ++ring) {
            int center = ring * ringSize + ringSize / 2;
            for (int i = -ringSize / 2; i < ringSize / 2; ++i) {
                for (int j = -ringSize / 2; j < ringSize / 2; ++j) {
                    int x = center + i;
                    int y = center + j;
                    if (x >= 0 && x < size && y >= 0 && y < size) {
                        double radius = sqrt(i * i + j * j);
                        if (radius > ringSize / 4 && radius < ringSize / 3) {
                            matrix(x, y) += 0.5 * exp(-pow(radius - ringSize / 3, 2) / 10);
                        }
                    }
                }
            }
        }

    } else if (type == "ladder_polymer") {
        // Leiter-Polymer: Zwei parallele Ketten mit Querverbindungen
        for (int i = 0; i < size / 2; ++i) {
            // Erste Kette
            if (i < size / 2 - 1) {
                matrix(i, i + 1) = 0.8;
                matrix(i + 1, i) = 0.8;
            }

            // Zweite Kette (parallel)
            int j = i + size / 2;
            if (j < size - 1) {
                matrix(j, j + 1) = 0.8;
                matrix(j + 1, j) = 0.8;
            }

            // Querverbindungen (Ladder-Rungs)
            if (i < size / 2 && j < size) {
                matrix(i, j) = 0.6;
                matrix(j, i) = 0.6;
            }
        }

    } else { // "statistical_network"
        // Statistisches Netzwerk: Zufällige Verbindungen mit clustering
        std::srand(42); // Reproduzierbare Zufälligkeit
        for (int i = 0; i < size; ++i) {
            for (int j = i; j < size; ++j) {
                double distance = abs(i - j);
                double prob = exp(-distance / 30.0) * (0.5 + 0.5 * (double)rand() / RAND_MAX);
                if (prob > 0.7) {
                    double strength = prob - 0.7;
                    matrix(i, j) = strength;
                    matrix(j, i) = strength;
                }
            }
        }
    }

    // Füge Rauschen hinzu
    MatrixXd noise = 0.05 * MatrixXd::Random(size, size);
    matrix += noise.cwiseAbs();

    return matrix;
}

int main()
{
    std::cout << "Matrix-basierte Topologie-Analyse für Polymere\n";
    std::cout << "==============================================\n\n";

    // Generiere verschiedene Polymer-Beispiele
    std::vector<std::string> polymerTypes = { "linear_chain", "ring_polymer",
        "ladder_polymer", "statistical_network" };

    for (const auto& type : polymerTypes) {
        std::cout << "Analysiere " << type << "...\n";

        // Generiere Beispiel-Matrix
        MatrixXd polymerMatrix = generatePolymerExample(type, 150);

        // Analysiere mit unserem Framework
        MatrixBasedTopologyAnalyzer analyzer(polymerMatrix);

        // Berechne Polymer-Metriken
        auto metrics = analyzer.analyzePolymerTopology();

        // Ausgabe der Metriken
        std::cout << "  Ring-Dichte: " << metrics.ringDensity << "\n";
        std::cout << "  Chain-Connectivity: " << metrics.chainConnectivity << "\n";
        std::cout << "  Cavity-Struktur: " << metrics.cavityStructure << "\n";
        std::cout << "  Topologische Komplexität: " << metrics.topologicalComplexity << "\n";
        std::cout << "  Anzahl H1-Features: " << metrics.ringDistribution.size() << "\n";

        // Visualisiere
        analyzer.visualizePolymerTopology(type, metrics);
        std::cout << "  → Visualisierungen erstellt: " << type << "_*.png\n\n";
    }

    std::cout << "Analyse-Empfehlungen:\n";
    std::cout << "- Ring-Polymere zeigen hohe H1-Dichte (Grüne Bereiche in Barcodes)\n";
    std::cout << "- Lineare Ketten haben niedrige topologische Komplexität\n";
    std::cout << "- Ladder-Polymere zeigen charakteristische H1-Muster\n";
    std::cout << "- Adaptive Skalierung verstärkt subtile Ring-Strukturen\n\n";

    std::cout << "Parameter-Tuning:\n";
    std::cout << "- Erhöhe Temperatur-Parameter für mehr Ring-Sensitivität\n";
    std::cout << "- Stärkere Dämpfung für bessere Auflösung schwacher Features\n";
    std::cout << "- Matrix-Größe anpassen je nach Polymer-Länge\n";

    return 0;
}
