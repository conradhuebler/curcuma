/*
 * <GFN1-xTB Parameter Database Implementation>
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Parameters extracted from TBLite (https://github.com/tblite/tblite)
 *   Source: src/tblite/xtb/gfn1.f90
 *   Copyright (C) 2019-2024 Sebastian Ehlert and contributors
 *   Licensed under LGPL-3.0-or-later
 *
 * Reference: S. Grimme, C. Bannwarth, P. Shushkov
 *   J. Chem. Theory Comput. 2017, 13, 1989-2009
 *   DOI: 10.1021/acs.jctc.7b00118
 *
 * This program is free software under GPL-3.0
 *
 * Claude Generated: Complete GFN1 parameter loader from TBLite arrays (March 2026)
 */

#include "gfn1_params_loader.h"
#include "gfn1_params.hpp"
#include "src/core/curcuma_logger.h"

#include <fmt/format.h>

namespace GFN1Params {

// =================================================================================
// Helper: Add symmetric pair parameters
// =================================================================================

void ParameterDatabase::addPairParams(int Z1, int Z2, double kpair,
                                      double kss, double ksp, double kpp)
{
    PairParams pair;
    pair.kpair = kpair;
    pair.kshell_ss = kss;
    pair.kshell_sp = ksp;
    pair.kshell_pp = kpp;

    // Store both orders for O(1) lookup
    m_pairs[{Z1, Z2}] = pair;
    m_pairs[{Z2, Z1}] = pair;
}

// =================================================================================
// TBLite get_dblock_row: returns d-block row (1,2,3) or 0 if not d-block
// TBLite gfn1.f90 line 756-770
// =================================================================================
static int getDBlockRow(int Z)
{
    if (Z > 20 && Z < 30) return 1;  // 3d: Sc-Cu
    if (Z > 38 && Z < 48) return 2;  // 4d: Y-Ag
    if (Z > 56 && Z < 80) return 3;  // 5d+4f: La-Au
    return 0;
}

// =================================================================================
// TBLite get_pair_param: kpair for element pair
// TBLite gfn1.f90 lines 723-754
// =================================================================================
static double getPairParam(int izp, int jzp)
{
    // kp array for d-block row pairs (TBLite gfn1.f90 line 727)
    static constexpr double kp[3] = {1.1, 1.2, 1.2};

    // Specific hardcoded pairs
    if (izp == 1 && jzp == 1)  return 0.96;
    if ((izp == 5 && jzp == 1) || (izp == 1 && jzp == 5))   return 0.95;
    if ((izp == 7 && jzp == 1) || (izp == 1 && jzp == 7))   return 1.04;
    if ((izp == 28 && jzp == 1) || (izp == 1 && jzp == 28)) return 0.90;
    if ((izp == 75 && jzp == 1) || (izp == 1 && jzp == 75)) return 0.80;
    if ((izp == 78 && jzp == 1) || (izp == 1 && jzp == 78)) return 0.80;
    if ((izp == 15 && jzp == 5) || (izp == 5 && jzp == 15)) return 0.97;
    if ((izp == 14 && jzp == 7) || (izp == 7 && jzp == 14)) return 1.01;

    // d-block row logic
    int itr = getDBlockRow(izp);
    int jtr = getDBlockRow(jzp);
    if (itr > 0 && jtr > 0) {
        return 0.5 * (kp[itr - 1] + kp[jtr - 1]);
    }

    return 1.0;
}

// =================================================================================
// TBLite kshell: shell-pair coupling factor
// TBLite gfn1.f90 lines 1037-1041
// kdiag = [1.85, 2.25, 2.0, 2.0, 2.0]  (indices 0..4)
// =================================================================================
static double getKShell(int k, int l)
{
    static constexpr double kdiag[5] = {1.85, 2.25, 2.0, 2.0, 2.0};
    // merge(2.08, (kdiag[l]+kdiag[k])/2, s-p pair)
    if ((k == 0 && l == 1) || (l == 0 && k == 1)) {
        return 2.08;
    }
    return (kdiag[l] + kdiag[k]) / 2.0;
}

// =================================================================================
// Load Default GFN1 Parameters — delegates to complete loader
// =================================================================================

bool ParameterDatabase::loadDefaultGFN1()
{
    return loadCompleteGFN1();
}

// =================================================================================
// Load Complete GFN1 Parameters (86 Elements from TBLite arrays)
// Claude Generated (March 2026): Replaces fabricated values with TBLite reference
// =================================================================================

bool ParameterDatabase::loadCompleteGFN1()
{
    m_elements.clear();
    m_pairs.clear();

    for (int Z = 1; Z <= 86; ++Z) {
        ElementParams elem;
        elem.atomic_number = Z;

        int nshell = GFN1_NSHELL[Z - 1];

        // Load shell-resolved parameters
        for (int shell_idx = 0; shell_idx < nshell; ++shell_idx) {
            int l = GFN1_ANG_SHELL[Z - 1][shell_idx];  // angular momentum of this shell

            ShellParams sp;

            // Self-energy: RAW eV → Hartree
            // TBLite gfn1.f90: p_selfenergy * evtoau
            // Indexed by [element][shell_position]
            sp.selfenergy = GFN1_SELF_ENERGY_RAW[Z - 1][shell_idx] * GFN1_EVTOAU;

            // CN-shift: RAW → Hartree (evtoau * 0.01)
            // TBLite gfn1.f90: p_kcn * evtoau * 0.01
            // Indexed by [element][shell_position]
            sp.kcn = GFN1_KCN_RAW[Z - 1][shell_idx] * GFN1_EVTOAU * 0.01;

            // Slater exponent (STO zeta) — dimensionless
            // Indexed by [element][shell_position]
            sp.gexp = GFN1_SLATER_EXPONENT[Z - 1][shell_idx];

            // Reference occupation — indexed by [element][angular_momentum]
            // For duplicate angular momenta (e.g. H's two s-shells), only the
            // first (valence) shell gets the occupation; polarization shells get 0.
            if (elem.shells.count(l) == 0) {
                sp.refocc = GFN1_REFERENCE_OCC[Z - 1][l];  // First shell with this l
            } else {
                sp.refocc = 0.0;  // Polarization shell (e.g. H's 2s)
            }

            // Shell polynomial — indexed by [element][angular_momentum]
            // TBLite gfn1.f90: p_shpoly * 0.01
            sp.shpoly = GFN1_SHPOLY_RAW[Z - 1][l] * 0.01;

            // Angular momentum, principal quantum number, and STO-NG primitives
            sp.angular_momentum = l;
            sp.principal_qn = GFN1_PRINCIPAL_QN[Z - 1][shell_idx];
            sp.num_primitives = GFN1_NUM_PRIMITIVES[Z - 1][shell_idx];

            // Shell list: ordered by shell index (matches TBLite)
            // Claude Generated (March 2026): Needed for H's two s-shells
            elem.shell_list.push_back(sp);

            // Angular momentum map: store first (valence) shell per l
            // For H (nshell=2, both s-shells): only store first s-shell
            if (elem.shells.count(l) == 0) {
                elem.shells[l] = sp;
            }
        }

        // Repulsion parameters (dimensionless)
        elem.rep_alpha = GFN1_REP_ALPHA[Z - 1];
        elem.rep_zeff = GFN1_REP_ZEFF[Z - 1];

        // Atomic radius: RAW Angstrom → Bohr
        // TBLite atomicrad.f90: atomic_rad * aatoau
        elem.rad = GFN1_ATOMIC_RAD_ANG[Z - 1] * GFN1_AATOAU;

        // Coulomb parameters: shell-resolved Hubbard
        // TBLite gfn1.f90 get_shell_hardness(): hubbard_parameter * shell_hubbard
        // shell_hubbard[l=0] is always 1.0, so gamma_ss = hubbard_parameter
        elem.gamma_ss = GFN1_HUBBARD[Z - 1] * GFN1_SHELL_HUBBARD[Z - 1][0];
        elem.gamma_sp = GFN1_HUBBARD[Z - 1] * GFN1_SHELL_HUBBARD[Z - 1][1];
        elem.gamma_pp = GFN1_HUBBARD[Z - 1] * GFN1_SHELL_HUBBARD[Z - 1][2];

        // Halogen bond parameters
        // TBLite gfn1.f90: halogen_bond * 0.1
        double xb_raw = GFN1_HALOGEN_BOND_RAW[Z - 1] * 0.1;
        if (xb_raw > 0.0) {
            // Store raw scaled value as strength; radius uses radscale * covalent radius
            elem.xb_strength = xb_raw;
            elem.xb_radius = GFN1_HALOGEN_RADSCALE * GFN1_ATOMIC_RAD_ANG[Z - 1];
        } else {
            elem.xb_strength = 0.0;
            elem.xb_radius = 0.0;
        }

        m_elements[Z] = elem;
    }

    // Load pair parameters from TBLite get_pair_param()
    // Only store pairs with kpair != 1.0 to save memory
    // The getPair() method returns default (1.0) for unstored pairs

    // Specific hardcoded pairs from TBLite gfn1.f90 lines 729-744
    addPairParams(1, 1, 0.96, getKShell(0, 0), getKShell(0, 1), getKShell(1, 1));
    addPairParams(1, 5, 0.95, getKShell(0, 0), getKShell(0, 1), getKShell(1, 1));
    addPairParams(1, 7, 1.04, getKShell(0, 0), getKShell(0, 1), getKShell(1, 1));
    addPairParams(1, 28, 0.90, getKShell(0, 0), getKShell(0, 2), getKShell(1, 2));
    addPairParams(1, 75, 0.80, getKShell(0, 0), getKShell(0, 2), getKShell(1, 2));
    addPairParams(1, 78, 0.80, getKShell(0, 0), getKShell(0, 2), getKShell(1, 2));
    addPairParams(5, 15, 0.97, getKShell(0, 0), getKShell(0, 1), getKShell(1, 1));
    addPairParams(7, 14, 1.01, getKShell(0, 0), getKShell(0, 1), getKShell(1, 1));

    // d-block pairs: all unique d-block element pairs
    // Row 1 (3d): Z=21-29, Row 2 (4d): Z=39-47, Row 3 (5d+4f): Z=57-79
    // kp = {1.1, 1.2, 1.2}
    // For same-row: kpair = kp[row]
    // For cross-row: kpair = 0.5*(kp[row_i] + kp[row_j])
    // We store these lazily — only generate for pairs that exist in the molecule
    // The getPair() fallback handles this via getPairParam()

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("gfn1_params", fmt::format("Loaded {} elements from TBLite reference",
                            getNumElements()));
    }

    return true;
}

// =================================================================================
// Pair access with TBLite fallback logic
// Claude Generated (March 2026)
// =================================================================================

bool ParameterDatabase::hasPair(int /*Z1*/, int /*Z2*/) const
{
    // Always return true — we have defaults for all pairs via getPairParam()
    return true;
}

PairParams ParameterDatabase::getPair(int Z1, int Z2) const
{
    // Check stored pairs first
    auto it = m_pairs.find({Z1, Z2});
    if (it != m_pairs.end()) {
        return it->second;
    }
    it = m_pairs.find({Z2, Z1});
    if (it != m_pairs.end()) {
        return it->second;
    }

    // Fallback: compute kpair from TBLite get_pair_param logic
    PairParams result;
    result.kpair = getPairParam(Z1, Z2);
    result.kshell_ss = getKShell(0, 0);
    result.kshell_sp = getKShell(0, 1);
    result.kshell_pp = getKShell(1, 1);
    return result;
}

} // namespace GFN1Params
