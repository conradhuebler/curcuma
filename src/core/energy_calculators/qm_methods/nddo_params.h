/*
 * <NDDO Parameter Database>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Centralized NDDO parameter database for MNDO, AM1, PM3 methods.
 * Parameters extracted from Ulysses (Menezes & Popowicz) MNDO.hpp.
 *
 * References:
 *   MNDO: M. J. S. Dewar, W. Thiel, JACS 99, 4899 (1977)
 *   AM1:  M. J. S. Dewar et al., JACS 107, 3902 (1985)
 *   PM3:  J. J. P. Stewart, J. Comput. Chem. 10, 209 (1989)
 *
 * Claude Generated: Centralized NDDO parameter database for Curcuma
 */

#pragma once

#include <map>
#include <string>
#include <vector>

namespace NDDOParams {

/**
 * @brief Complete NDDO parameter set for one element in one method.
 *
 * All energy parameters stored in eV, distances in Angstrom,
 * except Eisol which is in Hartree (atomic units).
 * Methods must convert to internal units as needed.
 */
struct ElementParams {
    // One-electron core integrals [eV]
    double U_ss = 0.0;
    double U_pp = 0.0;

    // Resonance integrals [eV]
    double beta_s = 0.0;
    double beta_p = 0.0;

    // Slater orbital exponents (method-specific, dimensionless)
    double zeta_s = 0.0;
    double zeta_p = 0.0;

    // Core repulsion parameter [1/Angstrom]
    double alpha = 0.0;

    // Multipole expansion parameters [Angstrom]
    double D1 = 0.0;
    double D2 = 0.0;

    // ERI exponents [Angstrom]
    double rho_s = 0.0;
    double rho_p = 0.0;

    // One-center two-electron integrals [eV]
    // These are the CRITICAL parameters that were previously missing
    double Gss = 0.0;   // (ss|ss)
    double Gpp = 0.0;   // (pp|pp)
    double Gsp = 0.0;   // (ss|pp) = (pp|ss)
    double Gp2 = 0.0;   // (pp|p'p') with different p directions
    double Hsp = 0.0;   // (sp|sp)

    // Isolated atom energy [Hartree] (for heat of formation)
    double Eisol = 0.0;

    // Gaussian correction parameters (AM1/PM3 core-core repulsion)
    // K in eV, L dimensionless, M in Angstrom
    std::vector<double> gauss_K;
    std::vector<double> gauss_L;
    std::vector<double> gauss_M;

    // Number of valence electrons
    int n_valence = 0;

    // Principal quantum number for valence shell
    int n_principal = 1;
};

/// Get complete MNDO parameter set (20 elements: H,Li,Be,B,C,N,O,F,Al,Si,P,S,Cl,Zn,Ge,Br,Sn,I,Hg,Pb)
const std::map<int, ElementParams>& getMNDOParams();

/// Get complete AM1 parameter set (23+ elements including Na,As,Se,Sb,Te)
const std::map<int, ElementParams>& getAM1Params();

/// Get complete PM3 parameter set (30 elements: H through Bi)
const std::map<int, ElementParams>& getPM3Params();

/// Check if element Z is supported for a given method
bool isElementSupported(const std::string& method, int Z);

/// Get the derived Hpp = 0.5*(Gpp - Gp2) for rotational invariance [eV]
inline double getHpp(const ElementParams& p) {
    double hpp = 0.5 * (p.Gpp - p.Gp2);
    return (hpp < 0.1) ? 0.1 : hpp;  // Enforce minimum as in Ulysses
}

/// Get number of core electrons for an element
int getCoreElectrons(int Z);

} // namespace NDDOParams
