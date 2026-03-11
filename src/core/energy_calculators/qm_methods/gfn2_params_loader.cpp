/*
 * <GFN2 Parameter Loader from xtb reference>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Reference: Bannwarth et al. JCTC 2019, 15, 1652-1671
 * Parameters extracted from xtb-6.6.1 param_gfn2-xtb.txt
 *
 * This program is free software under GPL-3.0
 */

#include "gfn2_params_loader.h"
#include "gfn2_xtb_params.hpp"
#include "src/core/curcuma_logger.h"
#include "src/core/units.h"

#include <cmath>
#include <map>
#include <stdexcept>
#include <fmt/format.h>

namespace GFN2Params {

// =================================================================================
// GFN2-xTB Parameters from xtb reference implementation
// All parameters extracted from external/xtb/param_gfn2-xtb.txt
// =================================================================================

ParameterDatabase::ParameterDatabase() {}

bool ParameterDatabase::loadDefaultGFN2() {
    return loadCompleteGFN2();
}

bool ParameterDatabase::loadCompleteGFN2() {
    m_elements.clear();
    m_pairs.clear();

    // Add all 86 elements from xtb param_gfn2-xtb.txt
    for (int Z = 1; Z <= 86; ++Z) {
        ElementParams elem;
        elem.atomic_number = Z;

        // Determine number of shells based on element
        int n_shells = 1;
        if (Z > 2) n_shells = 2;    // Period 2+: s, p
        if (Z > 12) n_shells = 3;   // Period 3+: s, p, d

        // Load shell-resolved parameters from xtb arrays
        for (int shell = 0; shell < n_shells; ++shell) {
            ShellParams sp;

            // Self-energy (level): raw values are in eV, convert to Hartree
            // TBLite gfn2.f90 line 468: p_selfenergy * evtoau
            // Claude Generated (March 2026): Added missing eV→Hartree conversion
            sp.selfenergy = GFN2_SELF_ENERGY[Z-1][shell] * CurcumaUnit::Energy::EV_TO_HARTREE;

            // Slater exponent (zeta) - dimensionless, no conversion needed
            sp.zeta = GFN2_ZETA[Z-1][shell];

            // CN-shift coefficient (kcn): raw values are in eV, convert to Hartree
            // TBLite gfn2.f90 line 343: p_kcn * evtoau
            sp.kcn = GFN2_KCN[Z-1][shell] * CurcumaUnit::Energy::EV_TO_HARTREE;

            // shpoly coefficient (Hamiltonian scaling)
            // TBLite gfn2.f90 line 296: p_shpoly * 0.01
            // Claude Generated (March 2026): Added missing 0.01 factor
            sp.shpoly = GFN2_POLY[Z-1][shell] * 0.01;

            // Reference occupation
            sp.refocc = GFN2_REFOCC[Z-1][shell];

            elem.shells[shell] = sp;
        }

        // Element-specific parameters
        elem.gamma_ss = GFN2_HUBBARD[Z-1];
        elem.gamma_pp = GFN2_HUBBARD[Z-1];  // Simplified: same as s

        // Repulsion parameters
        elem.rep_alpha = GFN2_REP_ALPHA[Z-1];
        elem.rep_zeff = GFN2_REP_ZEFF[Z-1];

        // Multipole radius from element-specific p_rad array (TBLite gfn2.f90 line 529)
        // Claude Generated (March 2026): Use exact element-specific radii
        elem.rad = GFN2_RAD[Z-1];
        elem.multipole_rad = GFN2_RAD[Z-1];

        // Hubbard third derivative (from GAM3)
        // TBLite gfn2.f90 line 175: p_hubbard_derivs * 0.1
        // Claude Generated (March 2026): Added missing 0.1 factor
        elem.hubbard_deriv = GFN2_HUBBARD_DERIV3[Z-1] * 0.1;

        // AES2 parameters (dipole and quadrupole kernels)
        // TBLite applies 0.01 scaling factor: p_dkernel = 0.01 * [...] (gfn2.f90 line 473)
        // Claude Generated (March 2026): Apply missing 0.01 factor
        elem.dkernel = GFN2_DPOL[Z-1] * 0.01;
        elem.qkernel = GFN2_QPOL[Z-1] * 0.01;

        // Valence coordination number (estimated from period)
        if (Z <= 10) elem.valence_cn = 4.0;
        else if (Z <= 18) elem.valence_cn = 6.0;
        else elem.valence_cn = 8.0;

        // Dispersion parameters (stub for now)
        elem.c6_base = 1.0;
        elem.r4_over_r2 = 1.0;

        m_elements[Z] = elem;
    }

    // GFN2 uses kpair = 1.0 for ALL pairs (confirmed from TBLite gfn2.f90 line 727)
    // No custom pair parameters needed - default PairParams already has kpair=1.0
    // Claude Generated (March 2026): Removed incorrect custom kpair values

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("gfn2_params", fmt::format("Loaded {} elements from xtb reference (kpair=1.0 for all pairs)",
                            getNumElements()));
    }

    return true;
}

const ElementParams& ParameterDatabase::getElement(int Z) const {
    auto it = m_elements.find(Z);
    if (it == m_elements.end()) {
        throw std::runtime_error("No GFN2 params for Z=" + std::to_string(Z));
    }
    return it->second;
}

const PairParams& ParameterDatabase::getPair(int Z1, int Z2) const {
    static PairParams default_pair;

    // Try exact match first
    auto it = m_pairs.find({Z1, Z2});
    if (it != m_pairs.end()) {
        return it->second;
    }

    // Try symmetric match
    it = m_pairs.find({Z2, Z1});
    if (it != m_pairs.end()) {
        return it->second;
    }

    // Return default (kpair=1.0, kshell_*=1.0)
    return default_pair;
}

bool ParameterDatabase::hasElement(int Z) const {
    return m_elements.count(Z) > 0;
}

bool ParameterDatabase::hasPair(int /*Z1*/, int /*Z2*/) const {
    // We have defaults for all pairs
    return true;
}

double ParameterDatabase::getShellHubbard(int Z, int shell) const {
    if (Z < 1 || Z > 86 || shell < 0 || shell > 2) {
        return 0.0;
    }

    // Claude Generated (March 2026): Use exact shell-resolved Hubbard from TBLite
    // Formula: γ_shell = HUBBARD_PARAMETER[Z] * (1.0 + SHELL_HUBBARD_CORR[Z][shell])
    // Reference: TBLite gfn2.f90 lines 127-171, get_shell_hardness() in line 712
    double gamma_base = GFN2_HUBBARD[Z-1];

    // Get shell-specific correction and compute actual Hubbard parameter
    // The corrections are relative to 1.0, so we add to base value
    double shell_correction = GFN2_SHELL_HUBBARD_CORR[Z-1][shell];
    return gamma_base * (1.0 + shell_correction);
}

void ParameterDatabase::addPairParams(int Z1, int Z2, double kpair,
                                      double kss, double ksp, double kpp) {
    PairParams params;
    params.Z1 = Z1;
    params.Z2 = Z2;
    params.kpair = kpair;
    params.kshell_ss = kss;
    params.kshell_sp = ksp;
    params.kshell_pp = kpp;

    // Store both directions (symmetric)
    m_pairs[{Z1, Z2}] = params;
    m_pairs[{Z2, Z1}] = params;
}

bool ParameterDatabase::loadFromTOML(const std::string& /*toml_file*/) {
    // Not implemented - use loadCompleteGFN2() instead
    CurcumaLogger::warn("TOML loading not implemented, using hardcoded xtb parameters");
    return loadCompleteGFN2();
}

} // namespace GFN2Params
