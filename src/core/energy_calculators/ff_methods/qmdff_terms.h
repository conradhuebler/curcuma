/*
 * <QMDFF-Terms for force field calculation>
 * Copyright (C) 2022 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * =============================================================================
 * ACKNOWLEDGMENT OF ORIGINAL WORK
 * =============================================================================
 *
 * QMDFF (Quantum Mechanically Derived Force Field) was developed by
 *   Prof. Stefan Grimme (University of Bonn)
 *   S. Grimme, J. Chem. Theory Comput. 2014, 10, 4497-4514
 *   DOI: 10.1021/ct500573f
 *
 * The energy expression implemented here is a line-by-line port of the Fortran
 * evaluator that used to ship with xtb as `src/qmdff.f90` (module xtb_qmdff).
 * It was deleted from xtb in commit b7dbd36 ("Refactoring of external drivers",
 * PR #568); the last version lives at `b7dbd36^:src/qmdff.f90` in the xtb git
 * history (a copy of the repository is vendored at external/xtb).
 *
 * Line references in the comments below ("qmdff.f90:NNN") point into that file.
 *
 * Claude Generated (August 2026): C++ port of ff_eg / ff_nonb / ff_hb and their
 * helper routines (abdamp, hbpara, eabhag, eabx, eabxag).
 *
 * =============================================================================
 * UNITS
 * =============================================================================
 * All kernels in this header work in ATOMIC UNITS: coordinates in Bohr,
 * energies in Hartree, charges in e. Callers holding Angstrom geometry (the
 * UFF/QMDFF convention inside FFWorkspace) must convert coordinates on the way
 * in and scale the returned gradient by 1.889726125 on the way out.
 *
 * The QMDFF-specific parameter structs below store lengths in Bohr. The generic
 * curcuma Bond/Angle structs keep their Angstrom convention, so the workspace
 * converts Bond::r0_ij at the call site.
 */

#pragma once

#include "qmdff_data.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>

namespace QMDFFTerms {

inline constexpr double kPi = 3.14159265358979323846264338327950;
inline constexpr double kTwoPi = 6.28318530717958623199592693708837;
inline constexpr double kSqrtPi = 1.77245385090551599275151910313925;

/// Number of Fourier terms a QMDFF torsion may carry (Fortran: ntterm = 4).
inline constexpr int kMaxTorsionTerms = 4;

// =============================================================================
// Parameter structs — one per Fortran term list
// =============================================================================

/**
 * @brief Bond stretch potential type (Fortran: bond(3,·)).
 *
 * QMDFF starts every 1,2 stretch as a Lennard-Jones-like term and promotes the
 * strong ones to Morse via `setmorse` (qmdff.f90:994-1019), so that weak bonds
 * keep the LJ shape while strong bonds dissociate correctly.
 */
enum class BondPotential : int {
    LennardJones = 0, ///< E = k (1 + (r0/r)^a - 2 (r0/r)^(a/2))
    Morse = 1         ///< E = k (1 - exp(-alpha (r - r0)))^2, alpha = a/(2 r0)
};

/**
 * @brief One QMDFF torsion or out-of-plane term (Fortran: tors/vtors).
 *
 * A regular torsion is a sum of up to four Fourier terms that is smoothly
 * switched between the +phi0 and the -phi0 branch by an error function
 * (qmdff.f90:432-451). `out_of_plane` selects the inversion branch
 * (tors(6,·) == 2, qmdff.f90:466-504), which uses a single term only.
 */
struct QMDFFTorsion {
    int i = 0, j = 0, k = 0, l = 0;
    int nterm = 1;                ///< tors(5,·): active Fourier terms (1..4)
    bool out_of_plane = false;    ///< tors(6,·) == 2
    double phi0 = 0.0;            ///< vtors(1,·): reference angle (radians)
    double scale = 1.0;           ///< vtors(2,·): overall barrier scaling
    std::array<double, kMaxTorsionTerms> rn{};    ///< vtors(3t+0,·): periodicity
    std::array<double, kMaxTorsionTerms> phase{}; ///< vtors(3t+1,·): phase (radians)
    std::array<double, kMaxTorsionTerms> v{};     ///< vtors(3t+2,·): amplitude (Hartree)
};

/**
 * @brief One non-covalent atom pair (Fortran: nci).
 *
 * `nk` is the topological class that selects the eps1/eps2 scaling factors,
 * see QMDFFData::eps1 / QMDFFData::eps2. The remaining members are the
 * element-pair constants that `ff_ini` (qmdff.f90:136-155) precomputes once, so
 * that the per-step loop stays a pure geometry evaluation.
 */
struct QMDFFNonCovalent {
    int i = 0, j = 0;
    int nk = 5;         ///< nci(3,·): 1 = 1,2 ... 5 = true NCI, 6 = no electrostatics
    double c6 = 0.0;    ///< c6ff(i,j) in atomic units (already scaled by s6)
    double r0_bj = 0.0; ///< r094ff: a1*sqrt(3 r2r4_i r2r4_j) + a2 (Bohr)
    double sr42 = 0.0;  ///< sr42ff: 3 * s8 * r2r4_i * r2r4_j
    double zab = 0.0;   ///< zabff: Z_i * Z_j (valence-electron product)
    double alpha = 0.0; ///< r0abff: repulsion range 16.5 / R0ab^1.5 (1/Bohr)
    double qq = 0.0;    ///< q_i * q_j (product of the QM-derived atomic charges)
};

/**
 * @brief One hydrogen- or halogen-bond triple (Fortran: hb/vhb).
 *
 * For a hydrogen bond `h` is the bridging hydrogen and both c1/c2 are used
 * (donor/acceptor charge scaling). For a halogen bond `h` is the halogen and
 * only c1 is used (qmdff.f90:244-254).
 */
struct QMDFFHBond {
    int a = 0, b = 0, h = 0;
    double c1 = 0.0, c2 = 0.0;
    bool halogen = false;
};

// =============================================================================
// Helper functions
// =============================================================================

/**
 * @brief Distance damping of bend and torsion terms (qmdff.f90:633-648).
 *
 * Bends and torsions must disappear when one of the participating bonds breaks,
 * otherwise dissociation curves are wrong. QMDFF multiplies every angle and
 * torsion by
 *      damp(r) = 1 / (1 + (r^2 / rcut)^2),
 *      rcut    = 3 * (R_A + R_B)^2      [Bohr^2, R = Mantina/Truhlar radii]
 * i.e. a soft cut-off at roughly twice the covalent bond length.
 *
 * @param zi        Atomic number of the first atom
 * @param zj        Atomic number of the second atom
 * @param r2        Squared distance in Bohr^2
 * @param damp      Out: damping factor in (0, 1]
 * @param ddamp     Out: 2 * d(damp)/d(r^2) — the extra factor 2 is the Fortran
 *                  convention so that d(damp)/dx_i = ddamp * (x_i - x_j)
 */
inline void abdamp(int zi, int zj, double r2, double& damp, double& ddamp)
{
    // Fortran: rcut = 3.0 * 3.5710642 * ((atomicRad(ati)+atomicRad(atj))*autoaa)**2
    // atomicRad is stored in Bohr there and converted back to Angstrom by autoaa;
    // our table is already in Angstrom, and 3.5710642 = (1/0.52917726)^2 turns
    // Angstrom^2 into Bohr^2. Keeping the literal constant reproduces the
    // reference arithmetic bit for bit.
    const double rad_sum = QMDFFData::atomicRadAngstrom(zi) + QMDFFData::atomicRadAngstrom(zj);
    const double rcut = 3.0 * 3.5710642 * rad_sum * rad_sum;
    const double ratio = r2 / rcut;
    const double rr = ratio * ratio;
    const double denom = 1.0 + rr;
    damp = 1.0 / denom;
    ddamp = -2.0 * 2.0 * rr / (r2 * denom * denom);
}

/**
 * @brief Charge-dependent HB/XB strength (qmdff.f90:753-760).
 *
 * c(q) = exp(-a q) / (exp(-a q) + b);  HB uses a = 10, b = 5,
 * XB uses a = -6.5, b = 1. A more negative partial charge on the donor or
 * acceptor therefore strengthens the hydrogen bond.
 */
inline double hbpara(double a, double b, double q)
{
    const double e = std::exp(-a * q);
    return e / (e + b);
}

/**
 * @brief Torsion angle in [0, pi] (qmdff.f90:1022-1082, function valijklff).
 *
 * QMDFF deliberately uses the acos branch (the `deter < 0` wrap to 2*pi - phi is
 * commented out in the reference) because only then is the analytic `dphidr`
 * gradient consistent. This matters here — unlike GFN-FF, the QMDFF torsion
 * energy is NOT even in phi (the erf switch distinguishes the two branches), so
 * the [0, pi] convention cannot be replaced by an atan2 angle.
 *
 * @param gradient Out: Matrix(4,3) rows = d phi / d r_{i,j,k,l}
 */
inline double torsionAngle(const Eigen::Vector3d& r_i,
                           const Eigen::Vector3d& r_j,
                           const Eigen::Vector3d& r_k,
                           const Eigen::Vector3d& r_l,
                           Eigen::Matrix<double, 4, 3>& gradient,
                           bool calculate_gradient)
{
    // Fortran ra/rb/rc
    const Eigen::Vector3d ra = r_j - r_i;
    const Eigen::Vector3d rb = r_k - r_j;
    const Eigen::Vector3d rc = r_l - r_k;

    Eigen::Vector3d na = ra.cross(rb);
    Eigen::Vector3d nb = rb.cross(rc);
    const double nan = na.norm();
    const double nbn = nb.norm();

    constexpr double eps = 1.0e-14;
    if (nan < eps || nbn < eps) {
        if (calculate_gradient)
            gradient.setZero();
        return 0.0;
    }

    const Eigen::Vector3d na_hat = na / nan;
    const Eigen::Vector3d nb_hat = nb / nbn;

    double snanb = na_hat.dot(nb_hat);
    if (std::abs(std::abs(snanb) - 1.0) < eps)
        snanb = (snanb < 0.0) ? -1.0 : 1.0;
    snanb = std::max(-1.0, std::min(1.0, snanb));

    const double phi = std::acos(snanb);

    if (!calculate_gradient)
        return phi;

    // -------------------------------------------------------------------------
    // dphidr (qmdff.f90 external, identical to gfnff_helpers.f90:514-583)
    // -------------------------------------------------------------------------
    gradient.setZero();

    const double cosphi = std::cos(phi);
    const double sinphi = std::sin(phi);
    const double nenner = nan * nbn * sinphi;
    if (std::abs(nenner) < eps)
        return phi;
    const double onenner = 1.0 / nenner;

    const Eigen::Vector3d rab = na.cross(rb);
    const Eigen::Vector3d rbb = nb.cross(rb);
    const Eigen::Vector3d rac = na.cross(rc);
    const Eigen::Vector3d rbc = nb.cross(rc);
    const Eigen::Vector3d rba = nb.cross(ra);
    const Eigen::Vector3d raa = na.cross(ra);

    const Eigen::Vector3d rapb = ra + rb;
    const Eigen::Vector3d rbpc = rb + rc;

    const Eigen::Vector3d rapba = rapb.cross(na);
    const Eigen::Vector3d rapbb = rapb.cross(nb);
    const Eigen::Vector3d rbpca = rbpc.cross(na);
    const Eigen::Vector3d rbpcb = rbpc.cross(nb);

    gradient.row(0) = onenner * (cosphi * nbn / nan * rab - rbb);
    gradient.row(1) = onenner * (cosphi * (nbn / nan * rapba + nan / nbn * rbc) - (rac + rapbb));
    gradient.row(2) = onenner * (cosphi * (nbn / nan * raa + nan / nbn * rbpcb) - (rba + rbpca));
    gradient.row(3) = onenner * (cosphi * nan / nbn * rbb - rab);

    return phi;
}

/**
 * @brief Out-of-plane (inversion) angle omega and its derivative.
 *
 * Port of `omega` / `domegadr` (xtb basic_geo.f90:366-394 / 401-...; byte
 * identical to the GFN-FF helper curcuma already ports in gfnff_geometry.h).
 * omega = asin(n_hat . v_hat), n = (r_i - r_j) x (r_k - r_j), v = r_l - r_i.
 *
 * @param gradient Out: Matrix(4,3) rows = d omega / d r_{i,j,k,l}
 */
inline double outOfPlaneAngle(const Eigen::Vector3d& r_i,
                              const Eigen::Vector3d& r_j,
                              const Eigen::Vector3d& r_k,
                              const Eigen::Vector3d& r_l,
                              Eigen::Matrix<double, 4, 3>& gradient,
                              bool calculate_gradient)
{
    const Eigen::Vector3d re = r_i - r_j;
    const Eigen::Vector3d rd = r_k - r_j;
    const Eigen::Vector3d rv = r_l - r_i;

    const Eigen::Vector3d rn = re.cross(rd);
    const double rnn = rn.norm();
    const double rvn = rv.norm();

    constexpr double eps = 1.0e-14;
    if (rnn < eps || rvn < eps) {
        if (calculate_gradient)
            gradient.setZero();
        return 0.0;
    }

    double sin_omega = (rn / rnn).dot(rv / rvn);
    sin_omega = std::max(-1.0, std::min(1.0, sin_omega));
    const double omega = std::asin(sin_omega);

    if (!calculate_gradient)
        return omega;

    gradient.setZero();
    const double nenner = rnn * rvn * std::cos(omega);
    if (std::abs(nenner) < eps)
        return omega;
    const double onenner = 1.0 / nenner;
    const double sinomega = std::sin(omega);

    const Eigen::Vector3d rdme = rd - re;
    const Eigen::Vector3d rve = rv.cross(re);
    const Eigen::Vector3d rne = rn.cross(re);
    const Eigen::Vector3d rdv = rd.cross(rv);
    const Eigen::Vector3d rdn = rd.cross(rn);
    const Eigen::Vector3d rvdme = rv.cross(rdme);
    const Eigen::Vector3d rndme = rn.cross(rdme);

    gradient.row(0) = onenner * (rdv - rn - sinomega * (rvn / rnn * rdn - rnn / rvn * rv));
    gradient.row(1) = onenner * (rvdme - sinomega * rvn / rnn * rndme);
    gradient.row(2) = onenner * (rve - sinomega * rvn / rnn * rne);
    gradient.row(3) = onenner * (rn - sinomega * rnn / rvn * rv);

    return omega;
}

// =============================================================================
// Energy terms — each returns the energy and (optionally) adds to a gradient
// =============================================================================

/**
 * @brief 1,2 bond stretch (qmdff.f90:316-354).
 *
 * LJ form   E = k [1 + (r0/r)^a - 2 (r0/r)^(a/2)]
 * Morse form E = k [1 - exp(-alpha (r - r0))]^2,  alpha = a / (2 r0)
 *
 * Both share the same three parameters (r0, k, a); `setmorse` only flips the
 * potential type, which is why a fitted parameter set can switch between them
 * without refitting.
 *
 * @param ri,rj      Positions in Bohr
 * @param r0         Equilibrium distance in Bohr (vbond(1,·))
 * @param k          Well depth / force constant in Hartree (vbond(2,·))
 * @param a          Shape exponent (vbond(3,·))
 * @param potential  LJ or Morse
 * @param gi,gj      Out: gradient contributions in Hartree/Bohr (accumulated)
 * @param do_gradient Compute gi/gj
 * @return Energy in Hartree
 */
inline double bondEnergy(const Eigen::Vector3d& ri, const Eigen::Vector3d& rj,
                         double r0, double k, double a, BondPotential potential,
                         Eigen::Vector3d& gi, Eigen::Vector3d& gj, bool do_gradient)
{
    const Eigen::Vector3d ra = ri - rj;
    const double r2 = ra.squaredNorm();
    const double r = std::sqrt(r2);
    if (r < 1.0e-12)
        return 0.0;

    double energy = 0.0;
    double fac = 0.0;

    if (potential == BondPotential::LennardJones) {
        const double u = r0 / r;
        const double ua = std::pow(u, a);
        const double ua2 = std::pow(u, 0.5 * a);
        energy = k * (1.0 + ua - 2.0 * ua2);
        fac = a * k * (-ua + ua2) / r2;
    } else {
        const double alpha = 0.5 * a / r0;
        const double e = std::exp(-alpha * (r - r0));
        energy = k * (1.0 - e) * (1.0 - e);
        fac = 2.0 * alpha * k * e * (1.0 - e) / r;
    }

    if (do_gradient) {
        gi += fac * ra;
        gj -= fac * ra;
    }
    return energy;
}

/**
 * @brief Angle bend with dissociation damping (qmdff.f90:357-402).
 *
 * Non-linear reference angle: E = k (cos(theta) - cos(theta0))^2
 * Linear reference angle (pi - theta0 < 1e-6): E = k (theta - theta0)^2
 * Both are multiplied by damp(r_ab) * damp(r_cb).
 *
 * @param ra,rb,rc   Positions of outer atom A, CENTRE B, outer atom C (Bohr).
 *                   Note the Fortran ordering: angl(1,·) is the centre.
 * @param za,zb,zc   Atomic numbers (for the damping radii)
 * @param theta0     Reference angle in radians (vangl(1,·))
 * @param k          Force constant in Hartree (vangl(2,·))
 * @param ga,gb,gc   Out: accumulated gradients in Hartree/Bohr
 */
inline double angleEnergy(const Eigen::Vector3d& ra, const Eigen::Vector3d& rb,
                          const Eigen::Vector3d& rc,
                          int za, int zb, int zc,
                          double theta0, double k,
                          Eigen::Vector3d& ga, Eigen::Vector3d& gb, Eigen::Vector3d& gc,
                          bool do_gradient)
{
    const Eigen::Vector3d vab = ra - rb;
    const Eigen::Vector3d vcb = rc - rb;
    const double rab2 = vab.squaredNorm();
    const double rcb2 = vcb.squaredNorm();
    if (rab2 < 1.0e-24 || rcb2 < 1.0e-24)
        return 0.0;

    const Eigen::Vector3d vp = vcb.cross(vab);
    const double rp = vp.norm() + 1.0e-14;

    double cosa = vab.dot(vcb) / std::sqrt(rab2 * rcb2);
    cosa = std::min(1.0, std::max(-1.0, cosa));
    const double theta = std::acos(cosa);

    double dampab = 0.0, ddampab = 0.0, dampcb = 0.0, ddampcb = 0.0;
    abdamp(za, zb, rab2, dampab, ddampab);
    abdamp(zc, zb, rcb2, dampcb, ddampcb);
    const double damp = dampab * dampcb;

    double ea = 0.0;
    double deddt = 0.0;
    if (kPi - theta0 < 1.0e-6) {
        // (near-)linear reference angle: harmonic in theta, cos would be flat here
        const double dt = theta - theta0;
        ea = k * dt * dt;
        deddt = 2.0 * k * dt;
    } else {
        const double dc = cosa - std::cos(theta0);
        ea = k * dc * dc;
        deddt = 2.0 * k * std::sin(theta) * (std::cos(theta0) - cosa);
    }

    if (do_gradient) {
        const Eigen::Vector3d deda = (-deddt / (rab2 * rp)) * vab.cross(vp);
        const Eigen::Vector3d dedc = (deddt / (rcb2 * rp)) * vcb.cross(vp);
        const Eigen::Vector3d dedb = deda + dedc;

        const Eigen::Vector3d term1 = ea * ddampab * dampcb * vab;
        const Eigen::Vector3d term2 = ea * ddampcb * dampab * vcb;

        ga += deda * damp + term1;
        gb -= dedb * damp + term1 + term2;
        gc += dedc * damp + term2;
    }

    return ea * damp;
}

/**
 * @brief Proper torsion, damped, with erf-switched +/-phi0 branches
 *        (qmdff.f90:413-464).
 *
 * E = s * sum_t v_t * [ (1-erf(phi-pi))/2 * (1 + cos(n_t (phi-phi0) + d_t))
 *                     + (1+erf(phi-pi))/2 * (1 + cos(n_t (phi+phi0-2pi) + d_t)) ]
 *
 * The error function switches smoothly between the two symmetry-equivalent
 * minima as phi passes pi, which keeps the potential differentiable although
 * phi itself is folded into [0, pi].
 *
 * @param ri..rl  Positions in Bohr, @param zi..zl atomic numbers
 * @param g       Out: Matrix(4,3) gradient contributions (accumulated)
 */
inline double torsionEnergy(const Eigen::Vector3d& ri, const Eigen::Vector3d& rj,
                            const Eigen::Vector3d& rk, const Eigen::Vector3d& rl,
                            int zi, int zj, int zk, int zl,
                            const QMDFFTorsion& t,
                            Eigen::Matrix<double, 4, 3>& g, bool do_gradient)
{
    const Eigen::Vector3d vab = ri - rj;
    const Eigen::Vector3d vcb = rj - rk;
    const Eigen::Vector3d vdc = rk - rl;

    double dampij = 0.0, damp2ij = 0.0;
    double dampjk = 0.0, damp2jk = 0.0;
    double dampkl = 0.0, damp2kl = 0.0;
    abdamp(zi, zj, vab.squaredNorm(), dampij, damp2ij);
    abdamp(zk, zj, vcb.squaredNorm(), dampjk, damp2jk);
    abdamp(zk, zl, vdc.squaredNorm(), dampkl, damp2kl);
    const double damp = dampij * dampjk * dampkl;

    Eigen::Matrix<double, 4, 3> dphi;
    const double phi = torsionAngle(ri, rj, rk, rl, dphi, do_gradient);

    double et = 0.0;
    double dij = 0.0;
    const double phipi = phi - kPi;
    const double ef = std::erf(phipi);
    const double expo = std::exp(-phipi * phipi) / kSqrtPi;

    const int nt = std::max(0, std::min(t.nterm, kMaxTorsionTerms));
    for (int it = 0; it < nt; ++it) {
        const double rn = t.rn[it];
        const double amp = t.v[it];
        const double dphi1 = phi - t.phi0;
        const double dphi2 = phi + t.phi0 - kTwoPi;
        const double c1 = rn * dphi1 + t.phase[it];
        const double c2 = rn * dphi2 + t.phase[it];
        const double e1 = amp * (1.0 + std::cos(c1));
        const double e2 = amp * (1.0 + std::cos(c2));
        et += 0.5 * (1.0 - ef) * e1 + (0.5 + 0.5 * ef) * e2;
        dij += -expo * e1 - 0.5 * (1.0 - ef) * amp * std::sin(c1) * rn
            + expo * e2 - (0.5 + 0.5 * ef) * amp * std::sin(c2) * rn;
    }

    et *= t.scale;
    dij *= t.scale * damp;

    if (do_gradient) {
        const Eigen::Vector3d term1 = et * damp2ij * dampjk * dampkl * vab;
        const Eigen::Vector3d term2 = et * damp2jk * dampij * dampkl * vcb;
        const Eigen::Vector3d term3 = et * damp2kl * dampij * dampjk * vdc;

        g.row(0) += dij * dphi.row(0) + term1.transpose();
        g.row(1) += dij * dphi.row(1) - term1.transpose() + term2.transpose();
        g.row(2) += dij * dphi.row(2) + term3.transpose() - term2.transpose();
        g.row(3) += dij * dphi.row(3) - term3.transpose();
    }

    return et * damp;
}

/**
 * @brief Out-of-plane / inversion term (qmdff.f90:466-504, tors(6,·) == 2).
 *
 * Two shapes, selected by the first periodicity entry:
 *   rn > 1e-6 : E = v (1 + cos(omega - omega0 + pi))   [single minimum]
 *   otherwise : E = v (cos(omega) - cos(omega0))^2     [double minimum]
 * Damping uses the three bonds that emanate from the central atom j.
 *
 * @param ri..rl  Positions in Bohr; j is the central atom
 * @param g       Out: Matrix(4,3) gradient contributions (accumulated)
 */
inline double outOfPlaneEnergy(const Eigen::Vector3d& ri, const Eigen::Vector3d& rj,
                               const Eigen::Vector3d& rk, const Eigen::Vector3d& rl,
                               int zi, int zj, int zk, int zl,
                               const QMDFFTorsion& t,
                               Eigen::Matrix<double, 4, 3>& g, bool do_gradient)
{
    const Eigen::Vector3d vab = rj - ri;
    const Eigen::Vector3d vcb = rj - rk;
    const Eigen::Vector3d vdc = rj - rl;

    double dampij = 0.0, damp2ij = 0.0;
    double dampjk = 0.0, damp2jk = 0.0;
    double dampjl = 0.0, damp2jl = 0.0;
    abdamp(zi, zj, vab.squaredNorm(), dampij, damp2ij);
    abdamp(zk, zj, vcb.squaredNorm(), dampjk, damp2jk);
    abdamp(zj, zl, vdc.squaredNorm(), dampjl, damp2jl);
    const double damp = dampij * dampjk * dampjl;

    Eigen::Matrix<double, 4, 3> domega;
    const double phi = outOfPlaneAngle(ri, rj, rk, rl, domega, do_gradient);

    double et = 0.0;
    double dij = 0.0;
    const double rn = t.rn[0];
    if (rn > 1.0e-6) {
        const double c1 = (phi - t.phi0) + kPi;
        et = (1.0 + std::cos(c1)) * t.scale;
        dij = -std::sin(c1) * t.scale * damp;
    } else {
        const double dc = std::cos(phi) - std::cos(t.phi0);
        et = t.scale * dc * dc;
        dij = 2.0 * t.scale * std::sin(phi) * (std::cos(t.phi0) - std::cos(phi)) * damp;
    }

    if (do_gradient) {
        const Eigen::Vector3d term1 = et * damp2ij * dampjk * dampjl * vab;
        const Eigen::Vector3d term2 = et * damp2jk * dampij * dampjl * vcb;
        const Eigen::Vector3d term3 = et * damp2jl * dampij * dampjk * vdc;

        g.row(0) += dij * domega.row(0) - term1.transpose();
        g.row(1) += dij * domega.row(1) + term1.transpose() + term2.transpose() + term3.transpose();
        g.row(2) += dij * domega.row(2) - term2.transpose();
        g.row(3) += dij * domega.row(3) - term3.transpose();
    }

    return et * damp;
}

/**
 * @brief Non-covalent pair: D3(BJ) dispersion + Coulomb + exponential repulsion
 *        (qmdff.f90:514-590, subroutine ff_nonb).
 *
 *   E_disp = -eps2 [ C6 / (r^6 + R0^6) + 3 s8 r2r4_i r2r4_j C6 / (r^8 + R0^8) ]
 *   E_es   = +eps1  q_i q_j / r
 *   E_rep  = +eps2  Z_i Z_j exp(-alpha r) / r        (only for r < 25 Bohr)
 *
 * The repulsion is what replaces the missing r^-12 wall of a classical FF: it
 * is a physically motivated Born-Mayer term with element-pair range parameters.
 *
 * @param r0_bj    Becke-Johnson radius R0 = a1 sqrt(3 r2r4_i r2r4_j) + a2 (Bohr)
 * @param sr42     3 * s8 * r2r4_i * r2r4_j
 * @param zab      Z_i * Z_j (valence-electron product)
 * @param alpha    Repulsion range 16.5 / R0ab^1.5 (1/Bohr)
 * @param qq       q_i * q_j
 * @param gi,gj    Out: accumulated gradients (Hartree/Bohr)
 * @param e_disp,e_es,e_rep Out: the three components (accumulated)
 */
inline void nonCovalentPair(const Eigen::Vector3d& ri, const Eigen::Vector3d& rj,
                            double c6, double r0_bj, double sr42,
                            double zab, double alpha, double qq,
                            int nk,
                            double& e_disp, double& e_es, double& e_rep,
                            Eigen::Vector3d& gi, Eigen::Vector3d& gj, bool do_gradient)
{
    const double eps1 = QMDFFData::eps1(nk);
    const double eps2 = QMDFFData::eps2(nk);

    const Eigen::Vector3d d = ri - rj;
    const double r2 = d.squaredNorm();
    const double r = std::sqrt(r2);
    if (r < 1.0e-12)
        return;
    const double oner = 1.0 / r;

    // --- London dispersion, D3 with Becke-Johnson damping ---
    if (eps2 != 0.0) {
        const double r4 = r2 * r2;
        const double r6 = r4 * r2;
        const double r06 = std::pow(r0_bj, 6);
        const double t6 = r6 + r06;
        const double t8 = r6 * r2 + r06 * r0_bj * r0_bj;
        const double c6t6 = c6 / t6;
        const double c6t8 = c6 / t8;
        const double t27 = sr42 * c6t8;

        e_disp -= (c6t6 + t27) * eps2;
        if (do_gradient) {
            const double drij = eps2 * (c6t6 * 6.0 * r4 / t6 + 8.0 * t27 * r6 / t8);
            gi += drij * d;
            gj -= drij * d;
        }
    }

    // --- Electrostatics from the QM-derived atomic charges ---
    if (eps1 != 0.0) {
        const double e0 = qq * oner * eps1;
        e_es += e0;
        if (do_gradient) {
            const double drij = e0 / r2;
            gi -= drij * d;
            gj += drij * d;
        }
    }

    // --- Born-Mayer repulsion (qmdff.f90:571-584) ---
    if (eps2 != 0.0 && r < 25.0) {
        const double t27 = zab * std::exp(-alpha * r);
        const double e0 = t27 * oner;
        e_rep += e0 * eps2;
        if (do_gradient) {
            const double drij = eps2 * t27 * (alpha * r + 1.0) * oner / r2;
            gi -= drij * d;
            gj += drij * d;
        }
    }
}

/**
 * @brief Hydrogen bond A...H-B with analytic gradient (qmdff.f90:766-881, eabhag).
 *
 * E = -DA(q) * damp_long(r_AB) * angle_term / r_AB^3
 *   damp_long = 1 / (1 + (r_AB / 8)^12)
 *   angle_term = [ (1 + cos(A-B,H)) / 2 ]^6
 *   DA        = (c_A r_AH^4 + c_B r_BH^4) / (r_AH^4 + r_BH^4)
 *
 * The r^-3 decay (rather than r^-4) and the switch of the reference angle to
 * whichever of A-H / B-H is shorter make the term symmetric in the two heavy
 * atoms, so a proton can be shared without a discontinuity.
 *
 * @param ra,rb,rh Positions of acceptor A, donor B and the bridging H (Bohr)
 * @param ga,gb,gh Out: accumulated gradients (Hartree/Bohr)
 * @return Energy in Hartree (0 if the interaction is negligible)
 */
inline double hydrogenBondEnergy(const Eigen::Vector3d& ra, const Eigen::Vector3d& rb,
                                 const Eigen::Vector3d& rh,
                                 double ca, double cb,
                                 Eigen::Vector3d& ga, Eigen::Vector3d& gb, Eigen::Vector3d& gh,
                                 bool do_gradient)
{
    constexpr double longcut = 8.0;
    constexpr double alp = 12.0;
    constexpr double alp3 = 6.0;

    const Eigen::Vector3d drah = ra - rh;
    const Eigen::Vector3d drbh = rb - rh;
    const Eigen::Vector3d drab = ra - rb;

    const double rab2 = drab.squaredNorm();
    const double rab = std::sqrt(rab2);
    const double rah2 = drah.squaredNorm();
    const double rah = std::sqrt(rah2);
    const double rbh2 = drbh.squaredNorm();
    const double rbh = std::sqrt(rbh2);
    if (rab < 1.0e-12 || rah < 1.0e-12 || rbh < 1.0e-12)
        return 0.0;

    // long-range damping, folded together with the 1/r_AB^3 decay
    const double ratio = std::pow(rab / longcut, alp);
    double rdampl = 1.0 / (1.0 + ratio);
    rdampl = rdampl / rab2 / rab;

    // reference angle uses the longer of A-H / B-H as the "donor" side
    double aprod = 0.0;
    double cosabh = 0.0;
    const bool ah_longer = (rah2 > rbh2);
    if (ah_longer) {
        aprod = 1.0 / rbh / rab;
        cosabh = -drbh.dot(drab) * aprod;
    } else {
        aprod = 1.0 / rah / rab;
        cosabh = drah.dot(drab) * aprod;
    }

    double aterm = 0.5 * (cosabh + 1.0);
    double apref = std::pow(aterm, alp3 - 1.0);
    aterm = aterm * apref;
    apref = alp3 * 0.5 * apref;

    const double rah4 = rah2 * rah2;
    const double rbh4 = rbh2 * rbh2;
    const double denom = 1.0 / (rah4 + rbh4);
    const double da = (ca * rah4 + cb * rbh4) * denom;

    const double eabh = -da * rdampl * aterm;
    if (eabh > -1.0e-8)
        return 0.0; // reference bails out on negligible interactions (qmdff.f90:828)

    if (do_gradient) {
        // donor-acceptor part
        double gi = 4.0 * (ca - cb) * rah2 * rbh4 * denom * denom;
        gi = -gi * rdampl * aterm;
        Eigen::Vector3d gav = gi * drah;

        gi = 4.0 * (cb - ca) * rbh2 * rah4 * denom * denom;
        gi = -gi * rdampl * aterm;
        Eigen::Vector3d gbv = gi * drbh;

        Eigen::Vector3d ghv = -gav - gbv;

        // long-range damping part
        gi = rdampl * rdampl * rab * (3.0 + (3.0 + alp) * ratio);
        gi = gi * da * aterm;
        const Eigen::Vector3d dg_long = gi * drab;
        gav += dg_long;
        gbv -= dg_long;

        // angle part
        Eigen::Vector3d dg, dg2;
        if (ah_longer) {
            dg = -aprod * drbh - cosabh * drab / rab2;
            dg2 = aprod * drab + cosabh * drbh / rbh2;
            const double gfac = -da * rdampl * apref;
            dg *= gfac;
            dg2 *= gfac;
            gav += dg;
            ghv += dg2;
            gbv -= dg + dg2;
        } else {
            dg = -aprod * drah + cosabh * drab / rab2;
            dg2 = -aprod * drab + cosabh * drah / rah2;
            const double gfac = -da * rdampl * apref;
            dg *= gfac;
            dg2 *= gfac;
            gbv += dg;
            ghv += dg2;
            gav -= dg + dg2;
        }

        ga += gav;
        gb += gbv;
        gh += ghv;
    }

    return eabh;
}

/**
 * @brief Halogen bond energy A...X-B (qmdff.f90:887-927, eabx).
 *
 * E = -c_A * damp_long(r_AB^2) * angle_term / r_XB^2, with a much longer
 * cut-off (longcut = 120 Bohr^2) and a softer alp = 6 than the HB term.
 * Only the energy is analytic in the reference; see halogenBondGradient().
 */
inline double halogenBondEnergy(const Eigen::Vector3d& ra, const Eigen::Vector3d& rb,
                                const Eigen::Vector3d& rx, double ca)
{
    constexpr double longcut = 120.0;
    constexpr double alp = 6.0;
    constexpr double alp3 = 6.0;

    const double rab2 = (ra - rb).squaredNorm();
    const double d2ik = (ra - rx).squaredNorm();
    const double d2jk = (rx - rb).squaredNorm();
    if (rab2 < 1.0e-24 || d2ik < 1.0e-24 || d2jk < 1.0e-24)
        return 0.0;

    const double dampl = 1.0 / (1.0 + std::pow(rab2 / longcut, alp));

    double term = 0.0;
    if (d2ik > d2jk) {
        const double xy = std::sqrt(rab2 * d2jk);
        term = 0.5 * (rab2 + d2jk - d2ik) / xy;
    } else {
        const double xy = std::sqrt(rab2 * d2ik);
        term = 0.5 * (rab2 + d2ik - d2jk) / xy;
    }
    const double aterm = std::pow(0.5 * (term + 1.0), alp3);

    return -ca * dampl * aterm / d2jk;
}

/**
 * @brief Halogen bond with central-difference gradient (qmdff.f90:933-986, eabxag).
 *
 * The reference differentiates eabx numerically with step 1e-6 Bohr and we keep
 * that, deliberately: the XB term contributes little to the total force and an
 * independent analytic derivation would be a silent source of disagreement with
 * the reference implementation.
 *
 * @param ga,gb,gx Out: accumulated gradients (Hartree/Bohr)
 */
inline double halogenBondEnergyGradient(const Eigen::Vector3d& ra, const Eigen::Vector3d& rb,
                                        const Eigen::Vector3d& rx, double ca,
                                        Eigen::Vector3d& ga, Eigen::Vector3d& gb,
                                        Eigen::Vector3d& gx, bool do_gradient)
{
    const double energy = halogenBondEnergy(ra, rb, rx, ca);
    if (!do_gradient)
        return energy;

    constexpr double step = 1.0e-6;
    const double inv = 1.0 / (2.0 * step);

    Eigen::Vector3d a = ra, b = rb, x = rx;
    for (int c = 0; c < 3; ++c) {
        a(c) += step;
        const double er = halogenBondEnergy(a, b, x, ca);
        a(c) -= 2.0 * step;
        const double el = halogenBondEnergy(a, b, x, ca);
        a(c) += step;
        ga(c) += (er - el) * inv;
    }
    for (int c = 0; c < 3; ++c) {
        b(c) += step;
        const double er = halogenBondEnergy(a, b, x, ca);
        b(c) -= 2.0 * step;
        const double el = halogenBondEnergy(a, b, x, ca);
        b(c) += step;
        gb(c) += (er - el) * inv;
    }
    for (int c = 0; c < 3; ++c) {
        x(c) += step;
        const double er = halogenBondEnergy(a, b, x, ca);
        x(c) -= 2.0 * step;
        const double el = halogenBondEnergy(a, b, x, ca);
        x(c) += step;
        gx(c) += (er - el) * inv;
    }
    return energy;
}

} // namespace QMDFFTerms
