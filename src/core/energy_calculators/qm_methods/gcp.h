/*
 * <Geometric counterpoise (gCP) + short-range basis (SRB) correction>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * The gCP model corrects the basis-set superposition error (BSSE) of a small
 * basis with a pairwise atomic term, and HF-3c adds a short-range basis
 * incompleteness term ("SRB", ORCA prints the sum as "gCP+bas"):
 *
 *   E_gCP  = sigma * sum_{A != B} emiss_A * xv_B * exp(-alpha R_AB^beta) / sqrt(S_AB)
 *   E_base = -qscal * sum_{A < B} (Z_A Z_B)^{3/2} * exp(-rscal (R0_AB)^{3/4} R_AB)
 *
 * S_AB is the overlap of two s-type Slater functions with exponents
 * eta * zeta_A, eta * zeta_B; xv_B = 1/sqrt(nbas_B - nel_B/2) is the number of
 * "virtual" basis functions of atom B; R0_AB is the D3 van-der-Waals pair radius.
 *
 * Literature:
 *   - H. Kruse, S. Grimme, J. Chem. Phys. 136, 154101 (2012)  (gCP model)
 *   - R. Sure, S. Grimme, J. Comput. Chem. 34, 1672 (2013)     (HF-3c, SRB term)
 * Ported from the reference implementation dftd3/simple-dftd3
 *   src/dftd3/gcp.f90 + src/dftd3/gcp/param.f90 (LGPL-3.0-or-later),
 * whose Python API (`dftd3.interface.GeometricCounterpoise(method="hf3c")`)
 * is the validation oracle (test_cases/qm_hf3c).
 *
 * Scope: H-Ne, the elements of the shipped MINIX basis. Only the parameter
 * set for HF-3c (hf/minix + base) is provided.
 *
 * Claude Generated: native gCP/SRB for HF-3c (Sep 2026)
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "src/core/global.h"

#include <string>
#include <vector>

namespace gcp {

/// Parameters of one gCP level (method/basis combination). Claude Generated.
struct Parameters {
    double sigma = 0.0;   ///< global scaling of the BSSE term
    double eta = 0.0;     ///< scaling of the Slater exponents
    double alpha = 0.0;   ///< exponential decay factor
    double beta = 0.0;    ///< power of the distance in the decay
    bool base = false;    ///< add the HF-3c short-range basis term
    double rscal = 0.0;   ///< SRB radius scaling
    double qscal = 0.0;   ///< SRB prefactor
    double cutoff = 60.0; ///< pair cutoff in Bohr (reference default)
};

/// HF-3c = HF/MINIX gCP + SRB (simple-dftd3 gcp/param.f90, case p_minix_bas).
Parameters hf3c();

/// True if every element of the molecule has gCP/SRB data (H-Ne).
bool supports(const std::vector<int>& atoms);

/**
 * @brief gCP (+ SRB if params.base) energy and, optionally, its gradient.
 *
 * @param atoms     atomic numbers
 * @param geom_bohr N x 3 Cartesian coordinates in Bohr
 * @param params    parameter set (e.g. gcp::hf3c())
 * @param gradient  if non-null, resized to N x 3 and filled with dE/dR in Eh/Bohr
 * @return energy in Hartree
 */
double energy(const std::vector<int>& atoms, const Matrix& geom_bohr,
              const Parameters& params, Matrix* gradient = nullptr);

/// The SRB/"base" part alone (Hartree) -- for the energy decomposition.
double baseEnergy(const std::vector<int>& atoms, const Matrix& geom_bohr,
                  const Parameters& params);

} // namespace gcp
