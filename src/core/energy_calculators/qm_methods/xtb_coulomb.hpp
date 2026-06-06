/*
 * <Isotropic shell-resolved Coulomb kernel for xTB>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Ports tblite's coulomb/charge/effective.f90 (Klopman–Ohno kernel with
 * averaged shell Hubbard parameters) for both GFN1 (harmonic average) and
 * GFN2 (arithmetic average). Exponent gexp=2 for both methods.
 *
 *   γ_ij(R) = 1 / ( R^gexp + γ̄_ij^(-gexp) )^(1/gexp)
 *
 * On-atom off-diagonal: γ_ij = γ̄_ij (pairwise average of shell Hubbards).
 * Diagonal:              γ_ii = γ_i   (raw shell Hubbard, not averaged with self).
 *
 * Per-shell Hubbard input follows tblite xtb/gfn{1,2}.f90:get_shell_hardness
 *   hardness(ish, isp) = hubbard_parameter(Z) * shell_hubbard(ang, Z)
 *
 * Claude Generated (Phase 3.3, Apr 2026). GPL-3.0.
 */

#pragma once
#ifndef CURCUMA_XTB_COULOMB_HPP_
#define CURCUMA_XTB_COULOMB_HPP_

#include "parameters/gfn1_params.hpp"
#include "parameters/gfn2_params.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <cstddef>
#include <vector>

namespace curcuma::xtb::coulomb {

enum class Method { GFN1, GFN2 };

/**
 * Per-shell Hubbard hardness after multiplication by the per-shell scaling.
 * Mirrors tblite xtb/gfn{1,2}.f90 :: get_shell_hardness.
 *
 * @param method GFN1 or GFN2 (selects the parameter tables).
 * @param z      atomic number (1..86).
 * @param ang    angular momentum (0=s, 1=p, 2=d) of the shell.
 */
inline double shell_hardness(Method method, int z, int ang)
{
    if (z < 1 || z > 86) return 0.0;
    if (ang < 0 || ang > 2) return 0.0;
    const double h = (method == Method::GFN1)
        ? curcuma::xtb::gfn1_params::hubbard_parameter[z - 1]
        : curcuma::xtb::gfn2_params::hubbard_parameter[z - 1];
    const double s = (method == Method::GFN1)
        ? curcuma::xtb::gfn1_params::shell_hubbard[z - 1][ang]
        : curcuma::xtb::gfn2_params::shell_hubbard[z - 1][ang];
    return h * s;
}

/** GFN1: harmonic average 2·gi·gj / (gi + gj). */
inline double harmonic_average(double gi, double gj)
{
    if (gi == 0.0 || gj == 0.0) return 0.0;
    return 2.0 / (1.0 / gi + 1.0 / gj);
}

/** GFN2: arithmetic average 0.5·(gi + gj). */
inline double arithmetic_average(double gi, double gj)
{
    return 0.5 * (gi + gj);
}

inline double average(Method method, double gi, double gj)
{
    return (method == Method::GFN1) ? harmonic_average(gi, gj)
                                    : arithmetic_average(gi, gj);
}

/**
 * Build the isotropic Coulomb matrix γ[nsh, nsh] per
 * tblite coulomb/charge/effective.f90 :: get_amat_0d.
 *
 *   Cross-atom (iat != jat), shell pair (ish, jsh):
 *       γ̄ = average( hardness(i,ish), hardness(j,jsh) )
 *       tmp = 1 / ( R^gexp + γ̄^(-gexp) )^(1/gexp)
 *       γ[ii+ish, jj+jsh] += tmp     (symmetric)
 *
 *   On-atom off-diagonal (iat == jat, ish != jsh):
 *       γ[ii+ish, ii+jsh] += average( hardness(i,ish), hardness(i,jsh) )
 *
 *   Diagonal:
 *       γ[ii+ish, ii+ish]  = hardness(i, ish)     (raw, no averaging)
 *
 * @param method        GFN1 or GFN2.
 * @param z_atoms       per-atom Z, length nat.
 * @param xyz_bohr      atom coords in bohr, length 3·nat, row-major.
 * @param shell_atom    shell → atom index, length nsh.
 * @param shell_ang     shell ang. mom., length nsh.
 */
inline Eigen::MatrixXd build_gamma_matrix(
    Method                   method,
    const std::vector<int>&  z_atoms,
    const std::vector<double>& xyz_bohr,
    const std::vector<int>&  shell_atom,
    const std::vector<int>&  shell_ang)
{
    const double gexp = (method == Method::GFN1)
        ? curcuma::xtb::gfn1_params::gexp
        : curcuma::xtb::gfn2_params::gexp;
    const double inv_gexp = 1.0 / gexp;

    const int nsh = static_cast<int>(shell_atom.size());
    Eigen::MatrixXd amat = Eigen::MatrixXd::Zero(nsh, nsh);

    // Precompute per-shell hardness.
    std::vector<double> g(nsh);
    for (int is = 0; is < nsh; ++is) {
        g[is] = shell_hardness(method, z_atoms[shell_atom[is]], shell_ang[is]);
    }

    for (int is = 0; is < nsh; ++is) {
        const int iat = shell_atom[is];
        for (int js = 0; js < nsh; ++js) {
            const int jat = shell_atom[js];
            if (is == js) {
                // Diagonal: raw hardness of this shell (not averaged with itself).
                amat(is, is) = g[is];
            } else if (iat == jat) {
                // On-atom cross-shell: pairwise average.
                amat(is, js) = average(method, g[is], g[js]);
            } else {
                // Cross-atom Klopman–Ohno kernel.
                const double dx = xyz_bohr[3 * iat + 0] - xyz_bohr[3 * jat + 0];
                const double dy = xyz_bohr[3 * iat + 1] - xyz_bohr[3 * jat + 1];
                const double dz = xyz_bohr[3 * iat + 2] - xyz_bohr[3 * jat + 2];
                const double r1 = std::sqrt(dx * dx + dy * dy + dz * dz);
                const double r1g = std::pow(r1, gexp);
                const double gam = average(method, g[is], g[js]);
                amat(is, js) = std::pow(r1g + std::pow(gam, -gexp), -inv_gexp);
            }
        }
    }
    return amat;
}

/**
 * Isotropic second-order electrostatic energy
 *   E_iso = 0.5 · qsh^T · γ · qsh
 * with qsh = reference shell population − actual shell population
 * (tblite sign convention, see scf/potential.f90).
 */
inline double energy_iso(const Eigen::MatrixXd& gamma,
                         const Eigen::VectorXd& q_sh)
{
    return 0.5 * q_sh.dot(gamma * q_sh);
}

/**
 * Potential v_sh[s] = Σ_s' γ[s, s'] · qsh[s']  (tblite add_to_vsh).
 */
inline Eigen::VectorXd potential_shell(const Eigen::MatrixXd& gamma,
                                       const Eigen::VectorXd& q_sh)
{
    return gamma * q_sh;
}

/**
 * Reference shell population n0_sh for a given element, per the parameter
 * table's reference_occ. Returns the first `nshell(Z)` entries.
 */
inline std::vector<double> reference_shell_populations(Method method, int z)
{
    std::vector<double> out;
    if (z < 1 || z > 86) return out;
    if (method == Method::GFN1) {
        const int ns = curcuma::xtb::gfn1_params::nshell[z - 1];
        for (int i = 0; i < ns; ++i) out.push_back(
            curcuma::xtb::gfn1_params::reference_occ[z - 1][i]);
    } else {
        const int ns = curcuma::xtb::gfn2_params::nshell[z - 1];
        for (int i = 0; i < ns; ++i) out.push_back(
            curcuma::xtb::gfn2_params::reference_occ[z - 1][i]);
    }
    return out;
}

/* ------------------------------------------------------------------------- *
 *  Third-order onsite correction — mirrors tblite coulomb/thirdorder.f90    *
 *                                                                           *
 *  GFN1 (atom-resolved):   E_3 = Σᵢ (1/3)·qᵢ³·Γ_i,   Γ_i = p_hubbard_derivs(Z_i) *
 *  GFN2 (shell-resolved):  E_3 = Σ_s (1/3)·qsh_s³·Γ_s,                     *
 *      Γ_s = p_hubbard_derivs(Z_{at(s)}) · shell_hubbard_derivs(ang_s)     *
 *                                                                           *
 *  Potential:                                                               *
 *      GFN1: v_at(i)  = qᵢ²·Γ_i                                            *
 *      GFN2: v_sh(s)  = qsh_s²·Γ_s                                         *
 * ------------------------------------------------------------------------- */

/** Per-shell Γ coefficient for third-order (shell-resolved, GFN2). */
inline double shell_hubbard_deriv_gfn2(int z, int ang)
{
    if (z < 1 || z > 86) return 0.0;
    if (ang < 0 || ang > 4) return 0.0;
    return curcuma::xtb::gfn2_params::p_hubbard_derivs[z - 1]
         * curcuma::xtb::gfn2_params::shell_hubbard_derivs[ang];
}

/** Per-atom Γ coefficient for third-order (atom-resolved, GFN1). */
inline double atom_hubbard_deriv_gfn1(int z)
{
    if (z < 1 || z > 86) return 0.0;
    return curcuma::xtb::gfn1_params::p_hubbard_derivs[z - 1];
}

/**
 * Third-order energy: GFN2 uses shell-resolved qsh, GFN1 uses atomic qat.
 *
 * @param method          GFN1 or GFN2.
 * @param z_atoms         per-atom Z.
 * @param shell_atom      shell → atom.
 * @param shell_ang       shell ang. mom.
 * @param q_at            atomic charges (used for GFN1).
 * @param q_sh            shell charges  (used for GFN2).
 */
inline double energy_third_order(
    Method                      method,
    const std::vector<int>&     z_atoms,
    const std::vector<int>&     shell_atom,
    const std::vector<int>&     shell_ang,
    const Eigen::VectorXd&      q_at,
    const Eigen::VectorXd&      q_sh)
{
    double e = 0.0;
    if (method == Method::GFN1) {
        for (std::size_t i = 0; i < z_atoms.size(); ++i) {
            const double q = q_at(static_cast<Eigen::Index>(i));
            e += q * q * q * atom_hubbard_deriv_gfn1(z_atoms[i]) / 3.0;
        }
    } else {
        for (std::size_t s = 0; s < shell_atom.size(); ++s) {
            const int    z  = z_atoms[shell_atom[s]];
            const double Gs = shell_hubbard_deriv_gfn2(z, shell_ang[s]);
            const double q  = q_sh(static_cast<Eigen::Index>(s));
            e += q * q * q * Gs / 3.0;
        }
    }
    return e;
}

/** Third-order atomic potential shift (GFN1 path). */
inline Eigen::VectorXd potential_third_order_atom(
    const std::vector<int>& z_atoms,
    const Eigen::VectorXd&  q_at)
{
    const Eigen::Index nat = static_cast<Eigen::Index>(z_atoms.size());
    Eigen::VectorXd v(nat);
    for (Eigen::Index i = 0; i < nat; ++i) {
        v(i) = q_at(i) * q_at(i) * atom_hubbard_deriv_gfn1(z_atoms[static_cast<std::size_t>(i)]);
    }
    return v;
}

/** Third-order shell potential shift (GFN2 path). */
inline Eigen::VectorXd potential_third_order_shell(
    const std::vector<int>& z_atoms,
    const std::vector<int>& shell_atom,
    const std::vector<int>& shell_ang,
    const Eigen::VectorXd&  q_sh)
{
    const Eigen::Index nsh = static_cast<Eigen::Index>(shell_atom.size());
    Eigen::VectorXd v(nsh);
    for (Eigen::Index s = 0; s < nsh; ++s) {
        const int z = z_atoms[shell_atom[static_cast<std::size_t>(s)]];
        const int l = shell_ang[static_cast<std::size_t>(s)];
        v(s) = q_sh(s) * q_sh(s) * shell_hubbard_deriv_gfn2(z, l);
    }
    return v;
}

} // namespace curcuma::xtb::coulomb

#endif // CURCUMA_XTB_COULOMB_HPP_
