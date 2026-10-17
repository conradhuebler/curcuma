/*
 * <Path collective variables (Branduardi s,z) - implementation.>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested pending - not TESTED/APPROVED).
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
 */

#include "src/capabilities/path_cv.h"

#include <cmath>
#include <limits>

#include "src/core/curcuma_logger.h"
#include "src/core/parameter_registry.h"

PathCV::PathCV(const std::vector<Geometry>& path, const Molecule& templ, double lambda)
    : m_path(path)
    , m_template(templ)
{
    if (m_path.size() < 2)
        CurcumaLogger::error("PathCV: need at least two reference frames.");
    m_lambda = (lambda > 0.0) ? lambda : computeLambdaFromPath(m_path);
}

double PathCV::computeLambdaFromPath(const std::vector<Geometry>& path)
{
    // lambda ~ 1/<d^2> with d the RMSD spacing of neighbouring frames, so that the
    // Gaussians of adjacent frames overlap appreciably. Using the plain (non-best-fit)
    // spacing here is sufficient for the scale estimate.
    if (path.size() < 2)
        return 1.0;
    const int nat = static_cast<int>(path[0].rows());
    double sum = 0.0;
    int cnt = 0;
    for (size_t i = 1; i < path.size(); ++i) {
        const double d2 = (path[i] - path[i - 1]).squaredNorm() / std::max(nat, 1);
        sum += d2;
        ++cnt;
    }
    const double mean_d2 = (cnt > 0) ? sum / cnt : 1.0;
    return (mean_d2 > 1e-12) ? 1.0 / mean_d2 : 1.0;
}

double PathCV::rmsdTo(int i, const Geometry& R, Geometry* grad) const
{
    // Claude Generated 2026: minimal RMSD of R onto frame i (Kabsch, exact) and its gradient
    // in the LABORATORY frame of R.
    //
    // Replaces the RMSDDriver call. That driver aligns the target onto the reference and
    // returns the gradient in the aligned (reference) frame, so the gradient was only right
    // for geometries already oriented like the reference. Measured: relative error of ds/dR
    // 4 % unrotated, 49 % after a 30 degree rotation, 99 % after 90 degrees (a
    // finite-difference check on a non-collinear molecule, see test_path_cv). Since NEB-MD
    // lets the images rotate freely (no COM/rotation removal), that error grew with time.
    //
    //   d^2 = 1/N sum_a |Rot (p_a - <p>) - (q_a - <q>)|^2 ,  Rot = V D U^T from SVD(P^T Q)
    //   d * dd/dp_a = Rot^T (Rot (p_a - <p>) - (q_a - <q>)) / N     (envelope theorem)
    //
    // The returned gradient keeps the sign convention of the old GradientTimesRMSD():
    // it is MINUS d*dd/dR, and the caller folds that sign into its chain rule.
    const Geometry& Q = m_path[i];
    const int N = static_cast<int>(R.rows());
    const Eigen::RowVector3d cp = R.colwise().mean();
    const Eigen::RowVector3d cq = Q.colwise().mean();
    const Eigen::MatrixXd P = R.rowwise() - cp;
    const Eigen::MatrixXd Qc = Q.rowwise() - cq;
    const Eigen::Matrix3d H = P.transpose() * Qc;
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
    if ((svd.matrixV() * svd.matrixU().transpose()).determinant() < 0)
        D(2, 2) = -1.0;
    const Eigen::Matrix3d Rot = svd.matrixV() * D * svd.matrixU().transpose();
    // rows: aligned_a = Rot * p_a  ->  P * Rot^T
    const Eigen::MatrixXd res = P * Rot.transpose() - Qc;
    const double d = std::sqrt(res.squaredNorm() / std::max(N, 1));
    if (grad)
        *grad = -(res * Rot) / static_cast<double>(std::max(N, 1)); // rows: Rot^T res_a
    return d;
}

void PathCV::evaluate(const Geometry& R, double* s_out, double* z_out, Geometry* ds_dR,
                      Geometry* dz_dR) const
{
    const int N = static_cast<int>(m_path.size());
    const int nat = static_cast<int>(R.rows());
    if (N < 2)
        return;

    std::vector<double> d(N, 0.0);
    std::vector<Geometry> gd(N);
    const bool need_grad = (ds_dR != nullptr) || (dz_dR != nullptr);
    for (int i = 0; i < N; ++i)
        d[i] = rmsdTo(i, R, need_grad ? &gd[i] : nullptr);

    // Softmax weights with the max subtracted for numerical stability.
    double emax = -std::numeric_limits<double>::max();
    for (int i = 0; i < N; ++i)
        emax = std::max(emax, -m_lambda * d[i] * d[i]);

    std::vector<double> w(N, 0.0);
    double wsum = 0.0;
    for (int i = 0; i < N; ++i) {
        w[i] = std::exp(-m_lambda * d[i] * d[i] - emax);
        wsum += w[i];
    }
    if (wsum < 1e-300)
        wsum = 1e-300;

    double s = 0.0;
    for (int i = 0; i < N; ++i)
        s += static_cast<double>(i) * w[i];
    s /= (wsum * (N - 1));

    if (s_out)
        *s_out = s;
    if (z_out)
        *z_out = -(1.0 / m_lambda) * (std::log(wsum) + emax);

    // dz/dR = sum_i (w_i/wsum) * 2 * d_i * dd_i/dR. The driver gradient is
    // -d(RMSD)/dR (see the sign note below), hence the leading minus here.
    if (dz_dR) {
        *dz_dR = Geometry::Zero(nat, 3);
        for (int i = 0; i < N; ++i) {
            if (gd[i].rows() != nat || gd[i].cols() != 3)
                continue;
            *dz_dR += (-(w[i] / wsum) * 2.0) * gd[i];
        }
    }

    if (!ds_dR)
        return;

    // ds/dR = sum_i [ i/(N-1) - s ] * (w_i/wsum) * (-2 lambda d_i) * dd_i/dR
    //
    // NOTE on the sign: RMSDDriver::Gradient() returns (ref - tar)/(RMSD*N), i.e.
    // MINUS d(RMSD)/dR_target (RMSD itself is sqrt(sum|delta|^2/N), see
    // RMSDFunctions::getRMSD). Folding that minus into the prefactor turns the
    // -2*lambda*d of the chain rule into +2*lambda*d. Getting this wrong shows up as
    // an exact factor -2 against finite differences, which is what the unit test
    // caught.
    *ds_dR = Geometry::Zero(nat, 3);
    for (int i = 0; i < N; ++i) {
        // The factor d_i that the chain rule produces is already contained in gd[i]
        // (= d_i * dd_i/dR, see rmsdTo), so it must NOT appear here a second time.
        // Keeping the two together is what removes the 1/RMSD singularity: a structure
        // sitting exactly on a reference frame is now an ordinary point.
        const double pref = (static_cast<double>(i) / (N - 1) - s) * (w[i] / wsum)
            * (2.0 * m_lambda);
        if (gd[i].rows() == nat && gd[i].cols() == 3)
            *ds_dR += pref * gd[i];
    }
}

double PathCV::zWallEnergy(const Geometry& R, double z_wall, double k, Geometry* grad,
                           double* z_out) const
{
    // Claude Generated 2026: flat-bottom wall on the distance-from-path z ("tube").
    //
    // Why: with s held by a constraint, the surface s = s0 is a (3N-7)-dimensional set that can
    // contain regions far from the reference path (retro-Diels-Alder: the dissociated ethylene
    // region at z ~ 1.4-3.8 A^2 next to the TS). Nothing in the window dynamics restrains z, so
    // the windows drift there. W confines them to z <= z_wall without touching the tube interior.
    //
    //   W(z)  = 0                          for z <= z_wall
    //   W(z)  = k/2 (z - z_wall)^2         for z >  z_wall
    //   dW/dR = W'(z) dz/dR = k (z - z_wall) dz/dR
    //
    // dz/dR comes from evaluate() in the laboratory frame (rmsdTo() returns lab-frame RMSD
    // gradients), so dW/dR can be added to dE/dR directly. W and dW/dR are continuous at
    // z = z_wall; the second derivative jumps there (harmless for velocity Verlet at the
    // stiffness chosen by zWallDefaultK).
    //
    // Note on the free energy: the sampled ensemble is exp(-(E+W)/kT), so the PMF obtained with
    // the wall is the PMF of the z-restricted system. It depends on z_wall.
    const int nat = static_cast<int>(R.rows());
    double s = 0.0, z = 0.0;
    Geometry ds, dz;
    evaluate(R, &s, &z, grad ? &ds : nullptr, grad ? &dz : nullptr);
    if (z_out)
        *z_out = z;
    const double excess = z - z_wall;
    if (excess <= 0.0 || k <= 0.0) {
        if (grad)
            *grad = Geometry::Zero(nat, 3);
        return 0.0;
    }
    if (grad)
        *grad = (k * excess) * dz;
    return 0.5 * k * excess * excess;
}

double PathCV::zWallDefaultK(double kT, double z_wall)
{
    // Claude Generated 2026: stiffness so that a thermal excursion is stopped close to the wall.
    // Require W = 3 kT at a penetration delta = 0.2 z_wall beyond the wall:
    //     k/2 delta^2 = 3 kT   ->   k = 6 kT / (0.2 z_wall)^2 = 150 kT / z_wall^2
    // The Boltzmann factor at that depth is exp(-3) = 0.05, so few samples reach z > 1.2 z_wall.
    // Example 300 K (kT = 9.50e-4 Eh), z_wall = 0.5 A^2: k = 0.57 Eh/A^4. Force at that depth:
    // |dW/dR| = k delta |dz/dR| with |dz/dR| ~ 2 d / sqrt(N) (whole 3N vector, d = RMSD ~
    // sqrt(z)); 16 atoms: 0.57 * 0.1 * 2 * 0.77 / 4 = 0.022 Eh/A, below typical bonded forces.
    if (z_wall <= 0.0)
        return 0.0;
    const double delta = 0.2 * z_wall;
    return 2.0 * 3.0 * kT / (delta * delta);
}

double PathCV::s(const Geometry& R) const
{
    double s_val = 0.0;
    evaluate(R, &s_val, nullptr, nullptr);
    return s_val;
}

double PathCV::z(const Geometry& R) const
{
    double z_val = 0.0;
    evaluate(R, nullptr, &z_val, nullptr);
    return z_val;
}
