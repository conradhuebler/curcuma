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

// Covalent radii from Pyykkö & Atsumi (Chem. Eur. J. 15, 2009, 188-197) in Angstrom
// Reference: gfnff_param.f90:381-405 (covalentRadD3 base values before *aatoau*4/3 scaling)
// Used for D3 CN calculation (with 4/3 scaling applied in calculateGFNFFCN)
// Claude Generated (Feb 16, 2026): Extended from 18 to 86 elements to fix Z>18 CN=0 bug
const std::vector<double> CNCalculator::COVALENT_RADII = {
    0.32, 0.46,                                                 // H, He
    1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,           // Li-Ne
    1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96,           // Na-Ar
    1.76, 1.54,                                                 // K, Ca
    1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09, // Sc-Zn
    1.12, 1.09, 1.15, 1.10, 1.14, 1.17,                        // Ga-Kr
    1.89, 1.67,                                                 // Rb, Sr
    1.47, 1.39, 1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, // Y-Cd
    1.28, 1.26, 1.26, 1.23, 1.32, 1.31,                        // In-Xe
    2.09, 1.76,                                                 // Cs, Ba
    1.62, 1.47, 1.58, 1.57, 1.56, 1.55, 1.51,                  // La-Eu
    1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53,                  // Gd-Yb
    1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32, // Lu-Hg
    1.30, 1.30, 1.36, 1.31, 1.38, 1.42                         // Tl-Rn
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

std::vector<double> CNCalculator::calculateGFNFFCN(
    const std::vector<int>& atoms,
    const Eigen::MatrixXd& geometry_bohr,
    double threshold,
    double kn,
    double cnmax)
{
    // Calculate GFN-FF coordination numbers using error function with log transformation
    // Formula: erfCN = 0.5 * (1 + erf(kn * (r - rcov) / rcov))
    //          CN_final = log(1 + e^cnmax) - log(1 + e^(cnmax - CN_raw))
    // Reference: gfnff_cn.f90:66-126, Spicher & Grimme, J. Chem. Theory Comput. 2020
    // Claude Generated - December 21, 2025

    const double ANG2BOHR = 1.8897259886;  // 1 Angstrom = 1.8897259886 Bohr

    std::vector<double> cn_values(atoms.size(), 0.0);

    for (size_t i = 0; i < atoms.size(); ++i) {
        int elem_i = atoms[i] - 1;  // Convert to 0-based

        if (elem_i < 0 || elem_i >= static_cast<int>(COVALENT_RADII.size())) {
            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::warn("calculateGFNFFCN: Element " + std::to_string(atoms[i]) +
                                   " out of range, using CN=0");
            }
            continue;
        }

        double rcov_i_ang = COVALENT_RADII[elem_i];
        double rcov_i_bohr = rcov_i_ang * ANG2BOHR;
        double cn_raw = 0.0;

        // GFN-FF CN scaling factor (Reference: gfnff_param.f90:381-404)
        // Usually 4/3 * standard covalent radius
        const double k_scaled = 4.0 / 3.0;

        for (size_t j = 0; j < atoms.size(); ++j) {
            if (i == j) continue;

            int elem_j = atoms[j] - 1;
            if (elem_j < 0 || elem_j >= static_cast<int>(COVALENT_RADII.size())) {
                continue;
            }

            double rcov_j_ang = COVALENT_RADII[elem_j];
            double rcov_j_bohr = rcov_j_ang * ANG2BOHR;
            double rcov_ij = k_scaled * (rcov_i_bohr + rcov_j_bohr);

            // Calculate distance in Bohr (geometry is in Bohr)
            Eigen::Vector3d pos_i = geometry_bohr.row(i);
            Eigen::Vector3d pos_j = geometry_bohr.row(j);
            double distance_sq = (pos_i - pos_j).squaredNorm();

            if (distance_sq > threshold) continue;

            double distance = std::sqrt(distance_sq);

            // GFN-FF error function coordination number
            // dr = (r - rcov) / rcov
            double dr = (distance - rcov_ij) / rcov_ij;
            double erfCN = 0.5 * (1.0 + std::erf(kn * dr));

            cn_raw += erfCN;
        }

        // Apply logarithmic transformation for numerical stability
        // CN_final = log(1 + e^cnmax) - log(1 + e^(cnmax - CN_raw))
        // This keeps CN in range [0, cnmax]
        cn_values[i] = std::log(1.0 + std::exp(cnmax)) - std::log(1.0 + std::exp(cnmax - cn_raw));

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("GFNFF_CN[" + std::to_string(i) + "] (Z=" + std::to_string(atoms[i]) +
                               "): raw=" + std::to_string(cn_raw) + " final=" + std::to_string(cn_values[i]));
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
