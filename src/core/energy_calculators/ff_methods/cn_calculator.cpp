/*
 * Coordination Number Calculator Implementation
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated - Consolidated from D3ParameterGenerator
 */

#include "cn_calculator.h"
#include "src/core/curcuma_logger.h"

#include <cmath>
#include <algorithm>

// Covalent radii from simple-dftd3 (Angstrom)
// Updated for elements H through Ar (18 elements)
// Data source: s-dftd3 atomic radii output
const std::vector<double> CNCalculator::COVALENT_RADII = {
    0.4267,  // 1 H
    0.6133,  // 2 He
    1.6000,  // 3 Li
    1.2533,  // 4 Be
    1.0267,  // 5 B
    1.0000,  // 6 C
    0.9467,  // 7 N
    0.8400,  // 8 O
    0.8533,  // 9 F
    0.8933,  // 10 Ne
    1.8667,  // 11 Na
    1.6667,  // 12 Mg
    1.5067,  // 13 Al
    1.3867,  // 14 Si
    1.4667,  // 15 P
    1.3600,  // 16 S
    1.3200,  // 17 Cl
    1.2800   // 18 Ar
};

std::vector<double> CNCalculator::calculateD3CN(
    const std::vector<int>& atoms,
    const Eigen::MatrixXd& geometry,
    double k1,
    double k2)
{
    // Calculate D3 coordination numbers using exponential counting function
    // Formula: CN_i = sum_{j≠i} 1 / (1 + exp(-k1 * (k2 * (R_cov_i + R_cov_j) / R_ij - 1)))
    // Reference: Grimme et al., J. Chem. Phys. 132, 154104 (2010)

    std::vector<double> cn_values(atoms.size(), 0.0);

    for (size_t i = 0; i < atoms.size(); ++i) {
        int elem_i = atoms[i] - 1;  // Convert to 0-based

        if (elem_i < 0 || elem_i >= static_cast<int>(COVALENT_RADII.size())) {
            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::warn("calculateCN: Element " + std::to_string(atoms[i]) +
                                   " out of range, using CN=0");
            }
            continue;
        }

        double rcov_i = COVALENT_RADII[elem_i];
        double cn = 0.0;

        for (size_t j = 0; j < atoms.size(); ++j) {
            if (i == j) continue;

            int elem_j = atoms[j] - 1;
            if (elem_j < 0 || elem_j >= static_cast<int>(COVALENT_RADII.size())) {
                continue;
            }

            double rcov_j = COVALENT_RADII[elem_j];

            // Calculate distance (geometry is in Angstrom)
            Eigen::Vector3d pos_i = geometry.row(i);
            Eigen::Vector3d pos_j = geometry.row(j);
            double r_ij = (pos_i - pos_j).norm();

            // D3 counting function
            double r_cov_sum = rcov_i + rcov_j;
            double arg = -k1 * (k2 * r_cov_sum / r_ij - 1.0);
            double count = 1.0 / (1.0 + std::exp(arg));

            cn += count;
        }

        cn_values[i] = cn;

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("CN[" + std::to_string(i) + "] (Z=" + std::to_string(atoms[i]) +
                               "): " + std::to_string(cn));
        }
    }

    return cn_values;
}

double CNCalculator::getCovalentRadius(int atomic_number)
{
    int idx = atomic_number - 1;  // Convert to 0-based

    if (idx < 0 || idx >= static_cast<int>(COVALENT_RADII.size())) {
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::warn("getCovalentRadius: Element " + std::to_string(atomic_number) +
                               " out of range");
        }
        return 1.0;  // Safe default
    }

    return COVALENT_RADII[idx];
}
