/*
 * <GFN1-xTB Parameter Database Loader>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Parameter extraction from TBLite (https://github.com/tblite/tblite)
 *   Copyright (C) 2019-2024 Sebastian Ehlert and contributors
 *   Licensed under LGPL-3.0-or-later
 *   Source: src/tblite/xtb/gfn1.f90
 *
 * Original GFN1-xTB method:
 *   S. Grimme, C. Bannwarth, P. Shushkov
 *   J. Chem. Theory Comput. 2017, 13, 1989-2009
 *   DOI: 10.1021/acs.jctc.7b00118
 *
 * This program is free software under GPL-3.0
 *
 * Claude Generated: GFN1 parameter loader infrastructure (November 2025)
 */

#pragma once

#include <map>
#include <string>
#include <vector>

/**
 * @brief GFN1-xTB Parameter Database
 *
 * This namespace contains structures and loaders for GFN1-xTB parameters
 * extracted from the TBLite reference implementation. Parameters are organized
 * by element and element pairs.
 *
 * Key Differences from GFN2:
 *   - Simpler shell structure (less shell-dependent parameters)
 *   - Halogen bond correction parameters
 *   - D3 dispersion instead of D4
 *   - Different parametrization philosophy (more empirical)
 *
 * Parameter Sources:
 *   - Element parameters: TBLite src/tblite/xtb/gfn1.f90
 *   - Halogen bond: Original GFN1 paper Table 4
 *   - Shell energies: Derived from TBLite TOML data
 */
namespace GFN1Params {

/**
 * @brief Shell-specific parameters for GFN1
 *
 * GFN1 uses a simpler shell structure than GFN2:
 *   - Self-energy with CN dependence (but simpler formula)
 *   - Gaussian exponent for STO basis
 *   - Reference occupation
 */
struct ShellParams {
    double selfenergy;   ///< Base self-energy E_ii (Hartree)
    double kcn;          ///< CN shift coefficient (Hartree/CN)
    double gexp;         ///< Gaussian exponent for STO
    double refocc;       ///< Reference occupation
};

/**
 * @brief Element-specific parameters for GFN1
 *
 * Stores all atomic parameters needed for GFN1-xTB calculations:
 *   - Shell parameters (s, p, optionally d)
 *   - Repulsion parameters
 *   - Coulomb parameters (gamma)
 *   - Halogen bond parameters (if applicable)
 */
struct ElementParams {
    int atomic_number;
    std::string symbol;

    // Shell parameters (key: 0=s, 1=p, 2=d)
    std::map<int, ShellParams> shells;

    // Repulsion parameters
    double rep_alpha;    ///< Repulsion exponent
    double rep_zeff;     ///< Effective nuclear charge for repulsion

    // Coulomb parameters
    double gamma_ss;     ///< Chemical hardness (s-s Coulomb integral)
    double gamma_sp;     ///< s-p Coulomb integral
    double gamma_pp;     ///< p-p Coulomb integral

    // Halogen bond parameters (for halogens only)
    double xb_radius;    ///< Halogen bond radius (Ångström)
    double xb_strength;  ///< Halogen bond strength (kcal/mol)
};

/**
 * @brief Element pair-specific parameters for GFN1
 *
 * Hamiltonian scaling factors for specific element pairs.
 * GFN1 uses simpler pair interactions than GFN2.
 */
struct PairParams {
    double kpair;        ///< Pairwise coupling strength
    double kshell_ss;    ///< s-s shell coupling
    double kshell_sp;    ///< s-p shell coupling
    double kshell_pp;    ///< p-p shell coupling
};

/**
 * @brief Complete GFN1 parameter database
 *
 * Container for all GFN1 parameters with lookup methods.
 * Maintains backward compatibility with legacy ArrayParameters.
 */
class ParameterDatabase {
public:
    ParameterDatabase() = default;

    // Load parameter sets
    bool loadDefaultGFN1();   ///< Load H, C, N, O (minimal set)
    bool loadCompleteGFN1();  ///< Load all 26 elements (Periods 1-5)

    // Element access
    bool hasElement(int Z) const { return m_elements.count(Z) > 0; }
    const ElementParams& getElement(int Z) const { return m_elements.at(Z); }

    // Pair access
    bool hasPair(int Z1, int Z2) const { return m_pairs.count({Z1, Z2}) > 0; }
    const PairParams& getPair(int Z1, int Z2) const { return m_pairs.at({Z1, Z2}); }

    // Statistics
    int getNumElements() const { return m_elements.size(); }
    int getNumPairs() const { return m_pairs.size() / 2; }  // Symmetric storage

private:
    // Helper: Add symmetric pair parameters
    void addPairParams(int Z1, int Z2, double kpair,
                      double kss, double ksp, double kpp);

    // Storage
    std::map<int, ElementParams> m_elements;  ///< Element parameters by Z
    std::map<std::pair<int,int>, PairParams> m_pairs;  ///< Pair parameters (symmetric)
};

} // namespace GFN1Params
