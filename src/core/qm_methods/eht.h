/*
 * <Extendend Hückel Theory Implementation in Cucuma. >
 * Copyright (C) 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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
#include "src/core/global.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <set>
#include <vector>

#include <Eigen/Dense>

#include "GTOIntegrals.hpp"
#include "STOIntegrals.hpp"

#include "interface/abstract_interface.h"

#include "json.hpp"

struct ehtSTO_6GHs {
    std::vector<double> alpha = { 0.3552322122E+02, 0.6513143725E+01, 0.1822142904E+01, 0.6259552659E+00, 0.2430767471E+00, 0.1001124280E+00 };
    std::vector<double> coeff = { 0.9163596281E-02, 0.4936149294E-01, 0.1685383049E+00, 0.3705627997E+00, 0.4164915298E+00, 0.1303340841E+00 };
};

struct ehtSTO_6GCs {
    std::vector<double> alpha = { 0.3049723950E+02, 0.6036199601E+01, 0.1876046337E+01, 0.7217826470E+00, 0.3134706954E+00, 0.1436865550E+00 };
    std::vector<double> coeff = { -0.1325278809E-01, -0.4699171014E-01, -0.3378537151E-01, 0.2502417861E+00, 0.5951172526E+00, 0.2407061763E+00 };
};

struct ehtSTO_6GCp {
    std::vector<double> alpha = { 0.3049723950E+02, 0.6036199601E+01, 0.1876046337E+01, 0.7217826470E+00, 0.3134706954E+00, 0.1436865550E+00 };
    std::vector<double> coeff = { 0.3759696623E-02, 0.3767936984E-01, 0.1738967435E+00, 0.4180364347E+00, 0.4258595477E+00, 0.1017082955E+00 };
};

struct ehtSTO_6GNs {
    std::vector<double> alpha = { 0.3919880787E+02, 0.7758467071E+01, 0.2411325783E+01, 0.9277239437E+00, 0.4029111410E+00, 0.1846836552E+00 };
    std::vector<double> coeff = { -0.1325278809E-01, -0.4699171014E-01, -0.3378537151E-01, 0.2502417861E+00, 0.5951172526E+00, 0.2407061763E+00 };
};

struct ehtSTO_6GNp {
    std::vector<double> alpha = { 0.3919880787E+02, 0.7758467071E+01, 0.2411325783E+01, 0.9277239437E+00, 0.4029111410E+00, 0.1846836552E+00 };
    std::vector<double> coeff = { 0.3759696623E-02, 0.3767936984E-01, 0.1738967435E+00, 0.4180364347E+00, 0.4258595477E+00, 0.1017082955E+00 };
};

struct ehtSTO_6GOs {
    std::vector<double> alpha = { 0.5218776196E+02, 0.1032932006E+02, 0.3210344977E+01, 0.1235135428E+01, 0.5364201581E+00, 0.2458806060E+00 };
    std::vector<double> coeff = { -0.1325278809E-01, -0.4699171014E-01, -0.3378537151E-01, 0.2502417861E+00, 0.5951172526E+00, 0.2407061763E+00 };
};

struct ehtSTO_6GOp {
    std::vector<double> alpha = { 0.5218776196E+02, 0.1032932006E+02, 0.3210344977E+01, 0.1235135428E+01, 0.5364201581E+00, 0.2458806060E+00 };
    std::vector<double> coeff = { 0.3759696623E-02, 0.3767936984E-01, 0.1738967435E+00, 0.4180364347E+00, 0.4258595477E+00, 0.1017082955E+00 };
};

struct STO_6G {
    std::vector<double> alpha, coeff;
    int x = 0, y = 0, z = 0;
    int index = 0;
    std::string sym = "";
    double e = 0;
};

typedef std::vector<STO::Orbital> Basisset;

class EHT : public QMInterface {
public:
    EHT();
    virtual bool InitialiseMolecule(const Mol& molecule) override
    {
        m_mol = molecule;
        return InitialiseMolecule();
    }
    virtual bool InitialiseMolecule() override;
    virtual double Calculation(bool gradient = false, bool verbose = false);

    Matrix MolecularOrbitals() const { return m_mo; }
    Vector Energies() const { return m_energies; }
    int NumElectrons() const { return m_num_electrons; }

private:
    // std::vector<STO_6G> MakeBasis();
    Basisset MakeBasis();
    Mol m_mol;
    Matrix m_H, m_S;

    Matrix MakeOverlap(Basisset& basisset);
    Matrix MakeH(const Matrix& S, const Basisset& basisset);

    int m_num_electrons = 0;

    Matrix m_mo;
    Vector m_energies;

    bool m_verbose = false;
};

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <chrono>
#include <iostream>
#include <omp.h>

// Klasse für Leistungsmessung
class Timer {
private:
    std::chrono::high_resolution_clock::time_point start_time;
    std::string operation_name;

public:
    Timer(const std::string& name)
        : operation_name(name)
    {
        start_time = std::chrono::high_resolution_clock::now();
    }

    ~Timer()
    {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << operation_name << ": " << duration.count() << " ms" << std::endl;
    }
};

// Effizient parallelisierte Berechnung von S^(-1/2)
inline Matrix calculateS12Efficient(const Matrix& S)
{
    Timer t("S^(-1/2) Berechnung");

    // Maximale Thread-Anzahl ausgeben
    int max_threads = omp_get_max_threads();
    std::cout << "Verfügbare Threads: " << max_threads << std::endl;

    // 1. Eigenwertzerlegung (kann nicht direkt mit OpenMP parallelisiert werden)
    Eigen::SelfAdjointEigenSolver<Matrix> es(S);
    const Matrix& U = es.eigenvectors();
    const Vector& D = es.eigenvalues();

    // 2. Diagonalmatrix D^(-1/2) berechnen - explizit parallelisiert
    Vector D_sqrt_inv(D.size());

#pragma omp parallel for
    for (int i = 0; i < D.size(); ++i) {
        D_sqrt_inv(i) = 1.0 / std::sqrt(std::max(D(i), 1e-12));
    }

    // 3. Erste Matrix-Multiplikation: Temp = U * D^(-1/2)
    // Diese Multiplikation kann effizienter direkt implementiert werden
    Matrix Temp(U.rows(), U.cols());

#pragma omp parallel for
    for (int i = 0; i < U.rows(); ++i) {
        for (int j = 0; j < U.cols(); ++j) {
            Temp(i, j) = U(i, j) * D_sqrt_inv(j);
        }
    }

    // 4. Zweite Matrix-Multiplikation: S^(-1/2) = Temp * U^T
    Matrix S_1_2(S.rows(), S.cols());

#pragma omp parallel for
    for (int i = 0; i < S_1_2.rows(); ++i) {
        for (int j = 0; j < S_1_2.cols(); ++j) {
            double sum = 0.0;
            for (int k = 0; k < U.cols(); ++k) {
                sum += Temp(i, k) * U(j, k); // Beachte: U^T(k,j) = U(j,k)
            }
            S_1_2(i, j) = sum;
        }
    }

    return S_1_2;
}

// Effizient parallelisierte Diagonalisierung der Fock-Matrix
inline void solveFockMatrix(const Matrix& S, const Matrix& H,
    Vector& energies, Matrix& mo)
{
    Timer t("Gesamt-Diagonalisierung");

    // 1. S^(-1/2) berechnen
    Matrix S_1_2 = calculateS12Efficient(S);

    {
        Timer t2("Fock-Transformation & Diagonalisierung");

        // 2. Transformierte Fock-Matrix berechnen: F' = S^(-1/2)^T * H * S^(-1/2)
        // Dies können wir in zwei Schritten machen, um die Parallelisierung zu maximieren
        Matrix temp(H.rows(), S_1_2.cols());

#pragma omp parallel for
        for (int i = 0; i < H.rows(); ++i) {
            for (int j = 0; j < S_1_2.cols(); ++j) {
                double sum = 0.0;
                for (int k = 0; k < H.cols(); ++k) {
                    sum += H(i, k) * S_1_2(k, j);
                }
                temp(i, j) = sum;
            }
        }

        Matrix F_prime(S_1_2.rows(), temp.cols());

#pragma omp parallel for
        for (int i = 0; i < S_1_2.rows(); ++i) {
            for (int j = 0; j < temp.cols(); ++j) {
                double sum = 0.0;
                for (int k = 0; k < S_1_2.cols(); ++k) {
                    sum += S_1_2(k, i) * temp(k, j); // Beachten Sie: S_1_2^T(i,k) = S_1_2(k,i)
                }
                F_prime(i, j) = sum;
            }
        }

        // 3. Diagonalisieren von F'
        Eigen::SelfAdjointEigenSolver<Matrix> es_F(F_prime);
        energies = es_F.eigenvalues();

        // 4. Molekülorbitale berechnen: MO = S^(-1/2) * C'
        const Matrix& C_prime = es_F.eigenvectors();
        mo = Matrix(S_1_2.rows(), C_prime.cols());

#pragma omp parallel for
        for (int i = 0; i < mo.rows(); ++i) {
            for (int j = 0; j < mo.cols(); ++j) {
                double sum = 0.0;
                for (int k = 0; k < S_1_2.cols(); ++k) {
                    sum += S_1_2(i, k) * C_prime(k, j);
                }
                mo(i, j) = sum;
            }
        }
    }
}
