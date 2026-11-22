/*
 * <CPCM-X: charge averaging, ring detection and the COSMO-RS restoring free energy>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Claude Generated (Nov 2026): port of CPCM-X 1.1.0 (LGPL-3.0-or-later)
 *   src/cpcmx/sigma_av.f90  average_one_charge, ortho_charge
 *   src/cpcmx/bonding.f90   det_bonds, det_rings, iter_ring (hb_grouping only feeds split_sigma, unused by 'crs')
 *   src/cpcmx/crs.f90       calcgas, E_dd, iterate_solvent, compute_solvent, compute_solute, state_correction
 *   src/cpcmx/broyden.f90   Broyden mixer of the solvent chemical potentials
 *   src/cpcmx/cpxcalc.f90   average_charge, init_bonding, calculate_solvation, dG
 * The Fortran evaluation order is kept where it changes rounding (the result is compared with xtb to 1e-10 Eh).
 * AI-generated, machine-tested only; human production testing pending.
 */

#include "cpcmx.h"
#include "cpcmx_internal.h"

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace curcuma::solvation::cpcmx {

namespace {

/// globals.f90 distance; the callers square it again (x**2.0_wp, folded to x*x by the compiler).
double distance(const Eigen::MatrixXd& a, int i, const Eigen::MatrixXd& b, int j)
{
    const double dx = a(i, 0) - b(j, 0), dy = a(i, 1) - b(j, 1), dz = a(i, 2) - b(j, 2);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * sigma_av.f90 average_one_charge: Gaussian-weighted average of the screening charge density over a radius r_av
 *   sigma'_i = sum_j sigma_j f_ij / sum_j f_ij,  f_ij = r_j^2 r_av^2 / (r_j^2 + r_av^2) exp(-d_ij^2 / (r_j^2 + r_av^2)),
 * with r_j^2 = A_j / pi the squared radius of a disc of the segment area (Klamt, J. Phys. Chem. 99, 2224 (1995)).
 */
Eigen::VectorXd averageOneCharge(double r_av, const Eigen::MatrixXd& xyz, const Eigen::VectorXd& charges,
                                 const Eigen::VectorXd& area)
{
    const int num = static_cast<int>(charges.size());
    Eigen::VectorXd av(num);
    const double r_av2 = r_av * r_av;
    for (int i = 0; i < num; ++i) {
        double counter = 0.0, denominator = 0.0;
        for (int j = 0; j < num; ++j) {
            const double r_u2 = area(j) / pi;
            const double d = distance(xyz, j, xyz, i);
            const double e = std::exp(-((d * d) / (r_u2 + r_av2)));
            counter = counter + (charges(j) * ((r_u2 * r_av2) / (r_u2 + r_av2)) * e);
            denominator = denominator + (((r_u2 * r_av2) / (r_u2 + r_av2)) * e);
        }
        av(i) = counter / denominator;
    }
    return av;
}

/// sigma_av.f90 ortho_charge: sigma_perp = sigma(2 r_av) - 0.816 sigma(r_av), the misfit correction descriptor.
Eigen::VectorXd orthoCharge(const Eigen::VectorXd& v, const Eigen::VectorXd& v0)
{
    Eigen::VectorXd vt(v.size());
    for (int i = 0; i < v.size(); ++i)
        vt(i) = v0(i) - 0.816 * v(i);
    return vt;
}

/// bonding.f90 det_bonds: atoms i < j are bonded if d_ij <= (r_cov,i + r_cov,j) * 1.15. Reproduced as in the source:
/// both radii are taken from atom i (radius=dict_get_key(cov_r,elements(i)) twice), 1.15 is a single-precision
/// literal, and only atoms up to the highest atom index that owns a segment are considered.
std::vector<std::vector<char>> detBonds(const Molecule& m)
{
    const int eles = *std::max_element(m.id.begin(), m.id.end());
    std::vector<std::vector<char>> bonded(eles, std::vector<char>(eles, 0));
    const double f115 = static_cast<double>(1.15f);
    for (int i = 0; i < eles; ++i) {
        const double cov_a = covalentRadius(m.z[i]);
        for (int j = i + 1; j < eles; ++j) {
            const double cov_b = covalentRadius(m.z[i]);
            if (distance(m.atom_xyz, i, m.atom_xyz, j) <= ((cov_a + cov_b) * f115)) {
                bonded[i][j] = 1;
                bonded[j][i] = 1;
            }
        }
    }
    return bonded;
}

/// bonding.f90 iter_ring: depth-first search for a path back to `check`; `visited` is not reset on backtracking
/// and the depth is limited to 11 bonds, as in the source.
void iterRing(int check, int i, int j, bool& is_ring, int eles, const std::vector<std::vector<char>>& bonded, int cut,
              std::vector<char>& visited)
{
    visited[j] = 1;
    for (int z = 0; z < eles; ++z) {
        if (bonded[j][z] && i != z && !visited[z]) {
            if (check == z) {
                is_ring = true;
                break;
            } else if (is_ring) {
                break;
            } else {
                if (cut >= 11)
                    break;
                iterRing(check, j, z, is_ring, eles, bonded, cut + 1, visited);
            }
        }
    }
}

/// bonding.f90 det_rings: ring atoms and N_ear, the effective number of ring atoms (omega N_ear enters dG_res):
/// a ring atom bonded to a non-ring atom counts 1, a ring atom without such a neighbour counts 2.
int detRings(const Molecule& m, const std::vector<std::vector<char>>& bonded)
{
    const int eles = static_cast<int>(bonded.size());
    std::vector<char> ring(eles, 0), visited(eles, 0);
    for (int i = 0; i < eles; ++i) {
        std::fill(visited.begin(), visited.end(), 0);
        bool is_ring = false;
        for (int j = 0; j < eles; ++j) {
            if (bonded[i][j])
                iterRing(i, i, j, is_ring, eles, bonded, 1, visited);
            if (is_ring)
                break;
        }
        ring[i] = is_ring;
    }
    int near = 0;
    for (int i = 0; i < eles; ++i) {
        if (!ring[i])
            continue;
        for (int j = 0; j < eles; ++j) {
            if (bonded[i][j] && !ring[j]) {
                near += 1;
                break;
            }
            if (j == eles - 1)
                near += 2;
        }
    }
    (void)m;
    return near;
}

/**
 * crs.f90 E_dd: interaction energy of two contacting segments (kcal/mol per effective contact area):
 *   E_hb     = c_hb max(0, sigma_acc - sigma_hb) min(0, sigma_don + sigma_hb)
 *   E_misfit = alpha'/2 (sigma_1 + sigma_2) ((sigma_1 + sigma_2) + f_corr (sigma_perp,1 + sigma_perp,2))
 */
double Edd(const Parameters& p, double sv1, double svt1, double sv2, double svt2)
{
    double svac, svdo;
    if (sv1 > sv2) {
        svac = sv1;
        svdo = sv2;
    } else {
        svac = sv2;
        svdo = sv1;
    }
    const double e_hb = p.chb * std::max(0.0, svac - p.shb) * std::min(0.0, svdo + p.shb);
    const double e_misfit = (p.aprime / 2.0) * (sv1 + sv2) * ((sv1 + sv2) + p.fcorr * (svt1 + svt2));
    return e_hb + e_misfit;
}

/// gfortran norm2 (scaled sum of squares), used by the Broyden mixer.
double norm2(const Eigen::VectorXd& x)
{
    double scale = 0.0, ssq = 1.0;
    for (int i = 0; i < x.size(); ++i) {
        if (x(i) != 0.0) {
            const double absx = std::abs(x(i));
            if (scale < absx) {
                const double r = scale / absx;
                ssq = 1.0 + ssq * r * r;
                scale = absx;
            } else {
                const double r = absx / scale;
                ssq += r * r;
            }
        }
    }
    return scale * std::sqrt(ssq);
}

/// broyden.f90 broyden_mixer (modified Broyden, D. D. Johnson, Phys. Rev. B 38, 12807 (1988)), 1d use only.
struct BroydenMixer {
    int ndim = 0, memory = 0, iter = 0;
    double damp = 0.4;
    Eigen::MatrixXd df, u, a;
    Eigen::VectorXd dq, dqlast, qlast_in, omega, q_in;

    BroydenMixer(int mem, int n, double d) : ndim(n), memory(mem), damp(d)
    {
        df.setZero(n, mem);
        u.setZero(n, mem);
        a.setZero(mem, mem);
        dq.setZero(n);
        dqlast.setZero(n);
        qlast_in.setZero(n);
        omega.setZero(mem);
        q_in.setZero(n);
    }
    void set(const Eigen::VectorXd& q) { q_in = q; }
    void diff(const Eigen::VectorXd& q) { dq = q - q_in; }
    void get(Eigen::VectorXd& q) const { q = q_in; }
    double error() const
    {
        double e = 0.0;
        const double n = static_cast<double>(dq.size());
        for (int i = 0; i < dq.size(); ++i)
            e = e + dq(i) * dq(i) / n;
        return std::sqrt(e);
    }
    void next()
    {
        ++iter;
        const int it1 = iter - 1;
        const double omega0 = 0.01, minw = 1.0, maxw = 100000.0, wfac = 0.01;
        if (iter == 1) {
            dqlast = dq;
            qlast_in = q_in;
            q_in = q_in + damp * dq;
            return;
        }
        const int k = it1 - 1;   // 0-based column of this step
        omega(k) = norm2(dq);
        if (omega(k) > (wfac / maxw))
            omega(k) = wfac / omega(k);
        else
            omega(k) = maxw;
        if (omega(k) < minw)
            omega(k) = minw;
        df.col(k) = dq - dqlast;
        double inv = std::max(norm2(df.col(k)), std::numeric_limits<double>::epsilon());
        inv = 1.0 / inv;
        df.col(k) = inv * df.col(k);
        Eigen::MatrixXd c(it1, 1);
        for (int i = 0; i < it1; ++i) {
            a(i, k) = df.col(i).dot(df.col(k));
            a(k, i) = a(i, k);
            c(i, 0) = omega(i) * df.col(i).dot(dq);
        }
        Eigen::MatrixXd beta(it1, it1);
        for (int i = 0; i < it1; ++i) {
            for (int r = 0; r < it1; ++r)
                beta(r, i) = omega(r) * omega(i) * a(r, i);
            beta(i, i) = beta(i, i) + omega0 * omega0;
        }
        // lineq: getrf + getrs(trans='t') solves beta^T x = c
        c = beta.transpose().partialPivLu().solve(c);
        u.col(k) = damp * df.col(k) + inv * (q_in - qlast_in);
        dqlast = dq;
        qlast_in = q_in;
        q_in = q_in + damp * dq;
        for (int i = 0; i < it1; ++i)
            q_in = q_in - omega(i) * c(i, 0) * u.col(i);
    }
};

/// crs.f90 iterate_solvent: one sweep of the self-consistent segment chemical potentials of the pure solvent
///   mu_j = -ln( 1/A sum_i A_i exp(-E_dd(i,j)/beta + mu_i) ),
/// updated in place (Gauss-Seidel order, a later j already sees the new mu_i of i < j), as in the source.
void iterateSolvent(Eigen::VectorXd& pot, const Eigen::VectorXd& area, const Eigen::MatrixXd& edd)
{
    double area_total = 0.0;
    for (int i = 0; i < area.size(); ++i)
        area_total += area(i);
    const int n = static_cast<int>(pot.size());
    for (int j = 0; j < n; ++j) {
        double temppot = 0.0;
        for (int i = 0; i < n; ++i)
            temppot = temppot + area(i) * std::exp(-edd(i, j) + pot(i));
        pot(j) = -std::log(temppot / area_total);
    }
}

/// crs.f90 compute_solvent: converge the solvent chemical potentials with the Broyden mixer (damping 0.4) until the
/// rms change is below conv_crit / autokcal.
bool computeSolvent(const Parameters& p, const Molecule& solvent, double T, int max_cycle, double conv_crit,
                    Eigen::VectorXd& pot, int& iterations, std::string& error)
{
    const int n = static_cast<int>(solvent.sv.size());
    const double beta = (R * jtokcal * T) / p.aeff;
    Eigen::MatrixXd edd(n, n);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            edd(i, j) = Edd(p, solvent.sv(j), solvent.svt(j), solvent.sv(i), solvent.svt(i)) / beta;
    pot.setZero(n);
    BroydenMixer mixer(max_cycle, n, 0.4);
    bool converged = false;
    int iter = 0;
    while (!converged) {
        if (iter > 0) {
            mixer.next();
            mixer.get(pot);
        }
        ++iter;
        mixer.set(pot);
        iterateSolvent(pot, solvent.area, edd);
        mixer.diff(pot);
        converged = mixer.error() < conv_crit / autokcal;
        if (iter >= max_cycle)
            break;
    }
    iterations = iter;
    if (!converged) {
        error = fmt::format("CPCM-X: solvent sigma profile did not converge after {} cycles", max_cycle);
        return false;
    }
    return true;
}

/// crs.f90 compute_solute: chemical potential of the solute segments in the converged solvent and of the solute
///   mu_X = beta sum_j A_j mu_j - lambda R T ln(A_solvent)   (kcal/mol).
double computeSolute(const Parameters& p, const Molecule& sol, const Molecule& solv, const Eigen::VectorXd& solv_pot,
                     double T)
{
    const double beta = (R * jtokcal * T) / p.aeff;
    double area_solv = 0.0;
    for (int i = 0; i < solv.area.size(); ++i)
        area_solv += solv.area(i);
    const int ns = static_cast<int>(sol.sv.size());
    Eigen::VectorXd sol_pot(ns);
    for (int j = 0; j < ns; ++j) {
        double temppot = 0.0;
        for (int i = 0; i < solv_pot.size(); ++i)
            temppot = temppot
                + (solv.area(i) * std::exp((-Edd(p, sol.sv(j), sol.svt(j), solv.sv(i), solv.svt(i)) / beta) + solv_pot(i)));
        sol_pot(j) = -std::log((1.0 / area_solv) * temppot);
    }
    double temppot = 0.0;
    for (int i = 0; i < ns; ++i)
        temppot = temppot + (sol.area(i) * sol_pot(i));
    return temppot * beta - (p.lambda * R * jtokcal * T * std::log(area_solv));
}

} // namespace

bool calculate(Molecule& solute, Molecule& solvent, const std::string& solvent_name, double T, int max_cycle,
               double conv_crit, double probe, OxygenRadiusRule rule, Result& res, std::string& error)
{
    res = Result{};
    const Parameters p = xtbParameters(solvent_name);
    if (solute.id.empty() || solvent.id.empty()) {
        error = "CPCM-X: no surface segments";
        return false;
    }

    // cpxcalc.f90 average_charge: sigma at r_av and 2 r_av, then the orthogonalised sigma
    if (p.rav == 0.0) {
        error = "CPCM-X: rav is zero";
        return false;
    }
    solvent.sv = averageOneCharge(p.rav, solvent.xyz, solvent.su, solvent.area);
    solute.sv = averageOneCharge(p.rav, solute.xyz, solute.su, solute.area);
    solvent.sv0 = averageOneCharge(p.rav * 2.0, solvent.xyz, solvent.su, solvent.area);
    solute.sv0 = averageOneCharge(p.rav * 2.0, solute.xyz, solute.su, solute.area);
    solvent.svt = orthoCharge(solvent.sv, solvent.sv0);
    solute.svt = orthoCharge(solute.sv, solute.sv0);

    // cpxcalc.f90 init_bonding: only det_rings of the solute reaches the energy (N_ear)
    res.near = detRings(solute, detBonds(solute));

    // crs.f90 calcgas: ideal screening energy and averaging correction
    {
        double ediel = 0.0, edielprime = 0.0;
        for (int i = 0; i < solute.sv.size(); ++i) {
            ediel = ediel + (solute.area(i) * solute.pot(i) * solute.su(i));
            edielprime = edielprime + (solute.area(i) * solute.pot(i) * solute.sv(i));
        }
        const double avcorr = (edielprime - ediel) / 2.0 * 0.8;
        const double dEreal = (solute.energy - solute.energy_gas) * autokcal;
        const double thermo = p.eta * R * jtokcal * T;
        res.dG_is = (dEreal - thermo) / autokcal;
        res.dG_cc = avcorr;
    }

    // crs.f90 compute_solvent / compute_solute, cpxcalc.f90 calculate_solvation
    Eigen::VectorXd solv_pot;
    if (!computeSolvent(p, solvent, T, max_cycle, conv_crit, solv_pot, res.solvent_iterations, error))
        return false;
    res.chem_pot_sol = computeSolute(p, solute, solvent, solv_pot, T);
    res.dG_shift = p.shift / autokcal;
    res.dG_res = (res.chem_pot_sol + p.omega * res.near) / autokcal;

    // crs.f90 state_correction: V_m = R T / 100 (L/mol at 1 bar), dG = -R T ln(rho V_m / M)
    const double density = solventDensity(solvent_name);
    double mass = 0.0;
    for (int z : solvent.z)
        mass = mass + atomicMass(z);
    if (density <= 0.0 || mass <= 0.0) {
        error = "CPCM-X: no density or mass for solvent " + solvent_name;
        return false;
    }
    const double V_m = (R * T) / 100.0;
    res.dG_ss = (-R * (jtokcal) * T * std::log((density * V_m) / mass)) / autokcal;

    // smd.f90 (numsa) CDS term
    return smdCds(solute, solvent_name, probe, rule, res.dG_smd, error);
}

} // namespace curcuma::solvation::cpcmx
