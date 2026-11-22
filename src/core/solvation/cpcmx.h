/*
 * <CPCM-X post-SCF solvation (extended conductor-like polarisable continuum model)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (Nov 2026): port of CPCM-X 1.1.0 (github.com/grimme-lab/CPCM-X, LGPL-3.0-or-later, M. Stahn,
 * S. Ehlert, S. Grimme) and of the SMD/surface part of numsa 0.2.0 (github.com/grimme-lab/numsa, LGPL-3.0-or-later)
 * as called by xtb 6.7.1 --cpcmx (xtb/src/solv/cpx.F90). LGPL-3 code may be used in this GPL-3 program; the solvent
 * database (CPCM-X DB/xtb/<solvent>.cosmo) is read at run time and is not part of the source tree.
 * AI-generated, machine-tested only; human production testing pending.
 *
 * Reference: M. Stahn, S. Ehlert, S. Grimme, J. Phys. Chem. A 127, 7036 (2023), doi:10.1021/acs.jpca.3c04382.
 */

#pragma once

#include <Eigen/Dense>

#include <string>
#include <vector>

namespace curcuma::solvation::cpcmx {

/**
 * CPCM-X in one paragraph: the solute is first converged in a perfect conductor (COSMO with epsilon = inf). Its
 * screening charge density sigma on the surface segments is the input of a COSMO-RS-like statistical
 * thermodynamics step (Klamt): the solvent, described by the sigma segments of its own conductor calculation (database),
 * is converged to self-consistent segment chemical potentials mu_S(sigma) (iterate_solvent with a Broyden mixer), and
 * the solute segments in that solvent give the chemical potential of the solute (compute_solute). The solvation free
 * energy is the sum of
 *   dG_is   ideal screening: E_COSMO(conductor) - E_gas - eta R T
 *   dG_av   averaging correction 0.4 sum_n A_n Phi_n (sigma'_n - sigma_n)
 *   dG_res  restoring free energy: mu_solute + omega N_ring
 *   dG_CDS  SMD cavity-dispersion-solvent-structure term (atomic surface tensions times SASA, numsa)
 *   dG_ss   state correction -R T ln(rho V_m / M) (1 bar gas to the pure solvent)
 *   dG_shift empirical shift (0 in the xtb parameter sets)
 * Units follow the Fortran code: segment data in Angstrom and e, energies in kcal/mol internally, results in Eh.
 */

/// One molecule as CPCM-X stores it (type.f90 molecule_data). Segment positions and atoms in Angstrom.
struct Molecule {
    std::vector<int> id;                 ///< 1-based atom of each segment
    Eigen::MatrixXd xyz;                 ///< nseg x 3 segment positions (Angstrom)
    Eigen::VectorXd area;                ///< segment area (Angstrom^2)
    Eigen::VectorXd su;                  ///< screening charge density sigma = q_n / A_n (e/Angstrom^2)
    Eigen::VectorXd pot;                 ///< solute potential on the segment (e/Bohr, read_cosmo: pot * BtoA)
    std::vector<int> z;                  ///< atomic numbers
    Eigen::MatrixXd atom_xyz;            ///< nat x 3 atom positions (Angstrom)
    double energy = 0.0;                 ///< COSMO (conductor) total energy, Eh
    double energy_gas = 0.0;             ///< gas-phase total energy, Eh
    // filled by the calculation
    Eigen::VectorXd sv, sv0, svt;        ///< averaged sigma (r_av), averaged sigma (2 r_av), orthogonalised sigma
};

/// COSMO-RS parameters of a CPCM-X parameter set (type.f90 parameter_type).
struct Parameters {
    double rav, aprime, fcorr, chb, shb, aeff, lambda, omega, eta, shift;
};

/// Contributions of the solvation free energy (cpxcalc.f90 calculation_type), all in Eh.
struct Result {
    double dG_is = 0.0, dG_cc = 0.0, dG_res = 0.0, dG_smd = 0.0, dG_ss = 0.0, dG_shift = 0.0;
    double chem_pot_sol = 0.0;           ///< solute chemical potential (kcal/mol)
    int near = 0;                        ///< effective number of ring atoms (det_rings)
    int solvent_iterations = 0;          ///< Broyden cycles of the solvent profile
    double dG() const { return dG_res + dG_cc + dG_is + dG_smd + dG_ss + dG_shift; }
};

/// xtb's conductor COSMO data of the solute, as cosmo.f90 writeCosmoFile would write it (Bohr, xtb autoaa).
struct ConductorData {
    std::vector<int> segment_atom;       ///< 0-based
    Eigen::MatrixXd segment_xyz;         ///< 3 x nseg (Bohr)
    Eigen::VectorXd charge;              ///< segment charge zeta (e)
    Eigen::VectorXd area;                ///< segment area (Angstrom^2)
    Eigen::VectorXd potential;           ///< potential (e/Bohr)
    std::vector<int> z;                  ///< atomic numbers
    Eigen::MatrixXd atom_xyz;            ///< nat x 3 (Bohr)
    double energy = 0.0;                 ///< total energy of the conductor SCF (Eh)
};

/**
 * @brief Solute input as read_cosmo (initialize.f90) gets it from xtb.cosmo.
 * @param reference_rounding true: round every number to the decimals xtb writes (segments f14.9, atom coordinates
 *        f19.14, energy f21.10) before the unit conversions of read_cosmo, so the result equals xtb --cpcmx, which
 *        passes the solute through that file; false: full precision (dG changes by <= 2e-10 Eh, measured on 4 molecules).
 */
Molecule soluteFromConductor(const ConductorData& data, bool reference_rounding);

/// Normalised Minnesota solvent name (data.f90 solvent_name), "" if unknown.
std::string solventName(const std::string& solvent);

/// Names accepted by solventName (the database file names), for error messages.
std::string knownSolvents();

/**
 * @brief Directory of the CPCM-X solvent database: the argument if not empty, else $CURCUMA_CPCMX_DB, else
 *        <executable dir>/cpcmx_db and <executable dir>/../share/curcuma/cpcmx_db. Returns "" if none exists.
 */
std::string findDatabase(const std::string& configured, std::vector<std::string>* tried = nullptr);

/**
 * @brief Read <db>/<name>.cosmo, the sigma segments of the solvent in xtb's conductor (internaldb.f90). The numbers
 *        are converted to single precision, as the Fortran include files that xtb compiles into its database are
 *        default-real literals.
 */
bool loadSolvent(const std::string& db_dir, const std::string& solvent, Molecule& out, std::string& error);

/// CPCM-X parameter set "xtb" for the solvent (internaldb.f90: xtb_water for water, xtb_other otherwise).
Parameters xtbParameters(const std::string& solvent_name);

/**
 * @brief The whole post-SCF step of xtb --cpcmx (cpx.F90 calculate_cpcmx): charge averaging, bonding and rings,
 *        restoring free energy, state correction and SMD CDS.
 * @param solute from soluteFromConductor(), energy_gas set
 * @param solvent from loadSolvent()
 * @param T temperature (xtb: 298.15 K), max_cycle / conv_crit of the solvent profile (xtb: 500, 1e-4 kcal/mol)
 * @param probe SMD probe radius (xtb: 0.4)
 */
/**
 * Oxygen radius rule of the SMD term for solvents with Abraham acidity alpha < 0.43 (SMD: R_O + 1.8 (0.43 - alpha)).
 * The amount 1.8 (0.43 - alpha) is added to a radius in Bohr, as in both Fortran versions.
 *  - Atom8:  the CPCM-X built into xtb 6.7.1 enlarges rad(8), the radius of the 8th atom whatever its element, and
 *            nothing for molecules with fewer than 8 atoms (measured: equals xtb 6.7.1 to 1e-10 Eh, see SOLVATION.md);
 *  - Oxygen: CPCM-X 1.1.0 enlarges every oxygen (where (to_number(symbols) == 8));
 *  - None:   no enlargement.
 */
enum class OxygenRadiusRule { Atom8, Oxygen, None };

/// Parse "atom8" | "oxygen" | "none" (case-insensitive); false if unknown.
bool parseOxygenRadiusRule(const std::string& text, OxygenRadiusRule& rule);

bool calculate(Molecule& solute, Molecule& solvent, const std::string& solvent_name, double T, int max_cycle,
               double conv_crit, double probe, OxygenRadiusRule rule, Result& result, std::string& error);

/// SMD CDS energy (Eh) of numsa 0.2.0 with CPCM-X's internal xtb parameters (smd.f90 calculate_cds_internal).
bool smdCds(const Molecule& solute, const std::string& solvent_name, double probe, OxygenRadiusRule rule, double& dG_cds,
            std::string& error);

} // namespace curcuma::solvation::cpcmx
