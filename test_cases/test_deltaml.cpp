/*
 * < Unit tests for the Delta-ML descriptor and sparse kernel-ridge model >
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
 * Claude Generated (July 2026).
 *
 * These are the checks that cannot be done from the CLI with useful precision. The
 * derivative tests matter most: an energy and a gradient that disagree produce results
 * that look perfectly reasonable in a single point and destroy optimisation and dynamics.
 */

#include "src/core/ml/descriptor.h"
#include "src/core/ml/krr_model.h"

#include <cstdio>
#include <random>
#include <string>
#include <vector>

using namespace curcuma::ml;

namespace {

int g_run = 0;
int g_failed = 0;

void check(bool ok, const std::string& what, const std::string& detail = "")
{
    ++g_run;
    if (ok) {
        std::printf("  PASS  %s%s%s\n", what.c_str(), detail.empty() ? "" : "  ", detail.c_str());
    } else {
        ++g_failed;
        std::printf("  FAIL  %s%s%s\n", what.c_str(), detail.empty() ? "" : "  ", detail.c_str());
    }
}

/*! \brief A short-ranged, permutation-invariant analytic target a local model can represent. */
double analyticDelta(const std::vector<int>& z, const Geometry& g, Matrix* grad)
{
    const int n = static_cast<int>(z.size());
    if (grad)
        *grad = Matrix::Zero(n, 3);
    double e = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double amp = 0.001 * (z[i] + z[j]);
            const Position d(g(j, 0) - g(i, 0), g(j, 1) - g(i, 1), g(j, 2) - g(i, 2));
            const double r = d.norm();
            const double t = amp * std::exp(-0.5 * (r - 1.5) * (r - 1.5));
            e += t;
            if (grad) {
                const double dt = -t * (r - 1.5);
                const Position u = d / r;
                for (int c = 0; c < 3; ++c) {
                    (*grad)(j, c) += dt * u(c);
                    (*grad)(i, c) -= dt * u(c);
                }
            }
        }
    }
    return e;
}

MLStructure makeStructure(std::mt19937& rng, const std::vector<int>& elements)
{
    std::uniform_real_distribution<double> d(-2.2, 2.2);
    MLStructure s;
    s.elements = elements;
    const int n = static_cast<int>(elements.size());
    s.geometry = Geometry(n, 3);
    for (int i = 0; i < n; ++i)
        for (int g = 0; g < 3; ++g)
            s.geometry(i, g) = d(rng);

    Matrix dg;
    const double delta = analyticDelta(elements, s.geometry, &dg);
    // A deliberately large offset, to check the element-count baseline absorbs it.
    s.energy_base = -123.456;
    s.energy_ref = s.energy_base + delta;
    s.gradient_base = Matrix::Zero(n, 3);
    s.gradient_ref = dg;
    s.has_gradients = true;
    return s;
}

// =================================================================================
// Descriptor
// =================================================================================

void testDescriptor()
{
    std::printf("\n[descriptor]\n");

    const std::vector<int> elements = { 6, 1, 1, 7, 6, 8, 1, 6 };
    const int natoms = static_cast<int>(elements.size());

    std::mt19937 rng(20260726);
    std::uniform_real_distribution<double> dist(-2.6, 2.6);
    Geometry geom(natoms, 3);
    for (int i = 0; i < natoms; ++i)
        for (int g = 0; g < 3; ++g)
            geom(i, g) = dist(rng);

    ACSFConfig cfg = ACSFConfig::defaultFor({ 1, 6, 7, 8 });
    std::string err;
    check(cfg.validate(err), "default configuration is valid", err);

    const ACSFDescriptor desc(cfg);
    const int nfeat = cfg.featureCount();

    Matrix X;
    std::vector<AtomJacobian> jac;
    desc.compute(elements, geom, X, &jac);
    check(X.rows() == natoms && X.cols() == nfeat, "descriptor has the expected shape");
    check(X.allFinite(), "descriptor values are finite");

    // Translation invariance of the descriptor is an exact analytic identity: the
    // Jacobian summed over all partners must vanish. This is a sharper statement about
    // the derivative bookkeeping than any finite-difference comparison.
    double max_trans = 0.0;
    for (int i = 0; i < natoms; ++i)
        for (int g = 0; g < 3; ++g)
            for (int f = 0; f < nfeat; ++f) {
                double s = 0.0;
                for (std::size_t p = 0; p < jac[i].partners.size(); ++p)
                    s += jac[i].d(3 * static_cast<int>(p) + g, f);
                max_trans = std::max(max_trans, std::abs(s));
            }
    check(max_trans < 1e-12, "Jacobian is translation invariant",
        "max |sum| = " + std::to_string(max_trans));

    // Rotation invariance of the values.
    Eigen::Matrix3d R
        = Eigen::AngleAxisd(0.7391, Eigen::Vector3d(0.3, -0.5, 0.81).normalized()).toRotationMatrix();
    Geometry rot(natoms, 3);
    for (int i = 0; i < natoms; ++i) {
        const Eigen::Vector3d q = R * Eigen::Vector3d(geom(i, 0), geom(i, 1), geom(i, 2));
        rot(i, 0) = q(0);
        rot(i, 1) = q(1);
        rot(i, 2) = q(2);
    }
    Matrix Xr;
    desc.compute(elements, rot, Xr, nullptr);
    const double rot_dev = (X - Xr).cwiseAbs().maxCoeff();
    check(rot_dev < 1e-10, "descriptor is rotation invariant",
        "max dev = " + std::to_string(rot_dev));

    // Finite differences. The tolerance is set by the FD roundoff floor, not by the
    // analytic derivative, so it cannot be tightened much further.
    const double h = 1e-6;
    double max_abs = 0.0;
    for (int c = 0; c < natoms; ++c) {
        for (int g = 0; g < 3; ++g) {
            Geometry gp = geom, gm = geom;
            gp(c, g) += h;
            gm(c, g) -= h;
            Matrix Xp, Xm;
            desc.compute(elements, gp, Xp, nullptr);
            desc.compute(elements, gm, Xm, nullptr);
            for (int i = 0; i < natoms; ++i) {
                int slot = -1;
                for (std::size_t p = 0; p < jac[i].partners.size(); ++p)
                    if (jac[i].partners[p] == c)
                        slot = static_cast<int>(p);
                for (int f = 0; f < nfeat; ++f) {
                    const double numeric = (Xp(i, f) - Xm(i, f)) / (2.0 * h);
                    const double analytic = (slot >= 0) ? jac[i].d(3 * slot + g, f) : 0.0;
                    max_abs = std::max(max_abs, std::abs(numeric - analytic));
                }
            }
        }
    }
    check(max_abs < 1e-7, "Jacobian matches finite differences",
        "max dev = " + std::to_string(max_abs));

    // An element the descriptor was never configured for must contribute nothing rather
    // than silently borrowing another element's parameters.
    ACSFConfig narrow = ACSFConfig::defaultFor({ 1, 6 });
    const ACSFDescriptor nd(narrow);
    Matrix Xn;
    nd.compute(elements, geom, Xn, nullptr);
    bool unknown_zero = true;
    for (int i = 0; i < natoms; ++i)
        if (elements[i] != 1 && elements[i] != 6 && Xn.row(i).cwiseAbs().maxCoeff() != 0.0)
            unknown_zero = false;
    check(unknown_zero, "unknown elements produce a zero descriptor");
}

// =================================================================================
// Model
// =================================================================================

void testModel()
{
    std::printf("\n[model]\n");

    std::mt19937 rng(4711);
    const std::vector<int> elements = { 6, 1, 1, 7, 6, 8, 1, 6, 1, 6 };

    MLDataset train, test;
    for (int k = 0; k < 120; ++k)
        train.add(makeStructure(rng, elements));
    for (int k = 0; k < 40; ++k)
        test.add(makeStructure(rng, elements));

    DeltaMLTrainingOptions opt;
    opt.n_centers = 300;
    opt.lambda = 1e-10;
    opt.sigma_energy = 1e-5;
    opt.sigma_force = 1e-4;
    opt.use_forces = true;

    DeltaMLModel model;
    DeltaMLTrainingReport rep;
    const bool trained = model.train(train, opt, rep);
    check(trained, "training succeeds");
    if (!trained)
        return;
    check(rep.n_centers > 0 && rep.n_force_rows > 0, "energy and gradient rows were used",
        std::to_string(rep.n_energy_rows) + " energy, " + std::to_string(rep.n_force_rows)
            + " gradient");

    const DeltaMLErrors te = model.evaluate(test);
    check(te.force_rmse < te.base_force_rmse, "held-out gradient error beats the base method",
        std::to_string(te.force_rmse) + " vs " + std::to_string(te.base_force_rmse));
    check(te.spearman > te.base_spearman, "held-out ranking beats the base method",
        std::to_string(te.spearman) + " vs " + std::to_string(te.base_spearman));

    // The two prediction paths must agree exactly; if they drift apart, a capability that
    // asks for a gradient sees a different surface than one that does not.
    const MLStructure& probe = test[0];
    Matrix g_analytic;
    const double e_grad = model.predictEnergyAndGradient(probe.elements, probe.geometry, g_analytic);
    const double e_plain = model.predictEnergy(probe.elements, probe.geometry);
    check(std::abs(e_grad - e_plain) < 1e-12, "gradient path and plain path give the same energy",
        std::to_string(std::abs(e_grad - e_plain)));

    // Finite differences on the model itself.
    const double h = 1e-5;
    double max_dev = 0.0;
    for (int i = 0; i < probe.AtomCount(); ++i) {
        for (int c = 0; c < 3; ++c) {
            Geometry gp = probe.geometry, gm = probe.geometry;
            gp(i, c) += h;
            gm(i, c) -= h;
            const double ep = model.predictEnergy(probe.elements, gp);
            const double em = model.predictEnergy(probe.elements, gm);
            max_dev = std::max(max_dev, std::abs((ep - em) / (2.0 * h) - g_analytic(i, c)));
        }
    }
    check(max_dev < 1e-8, "model gradient matches finite differences",
        "max dev = " + std::to_string(max_dev) + " Eh/Angstrom");

    // Persistence must be lossless, otherwise a reloaded model silently describes a
    // different surface than the one that was validated.
    const std::string path = "test_deltaml_roundtrip.json";
    check(model.save(path), "model saves");
    DeltaMLModel reloaded;
    check(reloaded.load(path), "model loads");
    Matrix g2;
    const double e2 = reloaded.predictEnergyAndGradient(probe.elements, probe.geometry, g2);
    check(std::abs(e2 - e_grad) < 1e-12, "reloaded model reproduces the energy");
    check((g2 - g_analytic).cwiseAbs().maxCoeff() < 1e-12,
        "reloaded model reproduces the gradient");
    std::remove(path.c_str());

    // Energy-only training must interpolate its training energies; this separates a
    // regression defect from a mere energy/gradient weighting trade-off.
    DeltaMLTrainingOptions eonly = opt;
    eonly.use_forces = false;
    DeltaMLModel emodel;
    DeltaMLTrainingReport erep;
    if (emodel.train(train, eonly, erep)) {
        check(erep.train_errors.energy_mae_centered < 1e-5,
            "energy-only fit interpolates its training energies",
            std::to_string(erep.train_errors.energy_mae_centered * 627.5094740631) + " kcal/mol");
    } else {
        check(false, "energy-only training succeeds");
    }

    // An untrained model must be inert rather than guessing.
    DeltaMLModel empty;
    check(empty.predictEnergy(probe.elements, probe.geometry) == 0.0,
        "untrained model returns no correction");
    check(empty.uncertainty(probe.elements, probe.geometry) == 1.0,
        "untrained model reports maximum unfamiliarity");
}

} // namespace

int main()
{
    std::printf("========================================\n");
    std::printf("Delta-ML unit tests\n");
    std::printf("========================================\n");

    testDescriptor();
    testModel();

    std::printf("\n----------------------------------------\n");
    std::printf("%d checks, %d failed\n", g_run, g_failed);
    std::printf("%s\n", g_failed == 0 ? "ALL PASSED" : "FAILURES PRESENT");
    return g_failed == 0 ? 0 : 1;
}
