/*
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This file is part of Curcuma - Native Solvation Module
 *
 * Extracted from Ulysses (Copyright (C) 2023- Filipe Menezes et al.)
 * Adapted for Curcuma by Claude (Anthropic AI Assistant)
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */

#ifndef SOLVATION_PARAMETERS_H
#define SOLVATION_PARAMETERS_H

#include <string>

namespace Curcuma {
namespace Solvation {

/**
 * @brief Get dielectric constant for a given solvent
 *
 * Provides dielectric constants (epsilon) for implicit solvation models.
 * Values taken from standard reference tables and optimized for GFN2/GFN1 methods.
 *
 * Reference: https://depts.washington.edu/eooptic/linkfiles/dielectric_chart[1].pdf
 *
 * @param solvent Solvent name (lowercase, e.g., "water", "acetone", "dmso")
 * @param method QM method name (e.g., "GFN2", "GFN1") - some solvents have method-specific values
 * @return Dielectric constant epsilon (dimensionless), returns 1.0 for unknown solvents
 *
 * @note Supported solvents:
 *   - Polar: water, methanol, ethanol, acetone, dmso, dmf, acetonitrile
 *   - Medium: thf, dichloromethane, chloroform, aniline
 *   - Non-polar: hexane, benzene, toluene, carbon disulfide
 *   - Special: "octanol wet" for partition coefficients
 */
inline double SolventDielectricConstant(const std::string& solvent, const std::string& method = "GFN2") {
    double epsilon = 1.0; // Vacuum default

    if (solvent == "acetone") {
        epsilon = 20.7;
    } else if (solvent == "acetonitrile") {
        epsilon = 37.5;
    } else if (solvent == "aniline") {
        epsilon = 6.8882;
    } else if (solvent == "benzaldehyde") {
        epsilon = 18.22;
    } else if (solvent == "benzene") {
        epsilon = (method == "GFN2") ? 7.0 : 2.27;
    } else if (solvent == "dichloromethane") {
        epsilon = 8.93;
    } else if (solvent == "chloroform") {
        epsilon = 4.81;
    } else if (solvent == "carbon disulfide") {
        epsilon = (method == "GFN2") ? 2.64 : 2.6;
    } else if (solvent == "dioxane") {
        epsilon = (method == "GFN2") ? 2.2099 : 2.25;
    } else if (solvent == "dmf") {
        epsilon = (method == "GFN2") ? 37.0 : 36.7;
    } else if (solvent == "dmso") {
        epsilon = 46.68;
    } else if (solvent == "ethanol") {
        epsilon = (method == "GFN2") ? 25.3 : 24.5;
    } else if (solvent == "diethyl ether") {
        epsilon = 4.33;
    } else if (solvent == "ethyl acetate") {
        epsilon = 6.02;
    } else if (solvent == "furane") {
        epsilon = 2.94;
    } else if (solvent == "hexadecane") {
        epsilon = (method == "GFN2") ? 1.88 : 2.09;
    } else if (solvent == "hexane") {
        epsilon = 1.88;
    } else if (solvent == "methanol") {
        epsilon = 32.7;
    } else if (solvent == "nitromethane") {
        epsilon = (method == "GFN2") ? 36.562 : 35.87;
    } else if (solvent == "octanol") {
        epsilon = 9.8629;
    } else if (solvent == "phenol") {
        epsilon = 12.4;
    } else if (solvent == "thf") {
        epsilon = (method == "GFN2") ? 7.6 : 7.58;
    } else if (solvent == "toluene") {
        epsilon = (method == "GFN2") ? 7.0 : 2.38;
    } else if (solvent == "water") {
        epsilon = (method == "GFN2") ? 80.2 : 80.1;
    } else if (solvent == "octanol wet") {
        epsilon = 9.8629;
    }

    return epsilon;
}

/**
 * @brief Get van der Waals radius for an atomic number
 *
 * Returns vdW radii in Angstrom for use in Born radii and SASA calculations.
 * Values optimized for GFN2/semi-empirical solvation models.
 *
 * @param atomic_number Atomic number (Z = 1 for H, 6 for C, etc.)
 * @return van der Waals radius in Angstrom
 *
 * @note Returns 2.0 Angstrom for unknown elements (safe default)
 */
inline double vdWaalsRadius(int atomic_number) {
    // Radii in Angstrom, optimized for GBSA/Born radii calculations
    // Source: Ulysses AtomicRadiipar.hpp
    switch (atomic_number) {
        case 1:  return 1.20;  // H
        case 2:  return 1.40;  // He
        case 3:  return 1.82;  // Li
        case 4:  return 1.53;  // Be
        case 5:  return 1.92;  // B
        case 6:  return 1.70;  // C
        case 7:  return 1.55;  // N
        case 8:  return 1.52;  // O
        case 9:  return 1.47;  // F
        case 10: return 1.54;  // Ne
        case 11: return 2.27;  // Na
        case 12: return 1.73;  // Mg
        case 13: return 1.84;  // Al
        case 14: return 2.10;  // Si
        case 15: return 1.80;  // P
        case 16: return 1.80;  // S
        case 17: return 1.75;  // Cl
        case 18: return 1.88;  // Ar
        case 19: return 2.75;  // K
        case 20: return 2.31;  // Ca
        case 25: return 2.05;  // Mn
        case 26: return 2.04;  // Fe
        case 29: return 1.40;  // Cu
        case 30: return 1.39;  // Zn
        case 35: return 1.85;  // Br
        case 53: return 1.98;  // I
        default: return 2.00;  // Safe default for unknown elements
    }
}

/**
 * @brief Get Born radius descreening parameter
 *
 * Descreening factors for different solvents and atomic types.
 * Used in GBOBC-II Born radii calculation.
 *
 * @param solvent Solvent name
 * @param atomic_number Atomic number
 * @return Descreening factor (dimensionless)
 *
 * @note Currently returns universal value 0.8 (safe default for all solvents/atoms)
 * @todo Extract full descreening parameter tables from Ulysses SolvationGFN2par.hpp
 */
inline double DeScreening(const std::string& solvent, int atomic_number) {
    // TODO: Extract full descreening parameter tables from Ulysses
    // For now: universal safe default that works reasonably for water
    return 0.8;
}

/**
 * @brief Get Born radius offset parameter
 *
 * Offset parameter used in GBOBC-II Born radii calculation.
 * Controls the relationship between vdW radius and effective Born radius.
 *
 * @return Born offset in Angstrom
 *
 * @note Value from Onufriev et al. 2004, optimized for GBOBC-II
 */
inline double BornOffset() {
    return 0.09; // Angstrom, from Ulysses implementation
}

/**
 * @brief Get Born radius scaling factor
 *
 * Solvent-specific scaling for final Born radii.
 *
 * @param solvent Solvent name
 * @return Scaling factor (dimensionless)
 *
 * @note Currently returns 1.0 (no scaling) for all solvents
 * @todo Extract solvent-specific scaling from Ulysses SolvationGFN2par.hpp
 */
inline double BornRadiusScale(const std::string& solvent) {
    // TODO: Extract solvent-specific scaling from Ulysses
    return 1.0; // No scaling for now
}

/**
 * @brief Get surface tension for SASA term
 *
 * Surface tension coefficient for solvent-accessible surface area energy.
 * Used in GBSA model: E_SASA = γ * SASA
 *
 * @param solvent Solvent name
 * @return Surface tension in kcal/(mol·Å²)
 *
 * @note Water default: 0.0072 kcal/(mol·Å²)
 * @todo Extract solvent-specific values from Ulysses
 */
inline double SurfaceTension(const std::string& solvent) {
    // TODO: Extract from Ulysses SolvationGFN2par.hpp
    if (solvent == "water") {
        return 0.0072; // kcal/(mol·Å²)
    }
    return 0.005; // Conservative default
}

} // namespace Solvation
} // namespace Curcuma

#endif // SOLVATION_PARAMETERS_H
