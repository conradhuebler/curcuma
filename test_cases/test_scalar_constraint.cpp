/*
 * <NVE energy conservation of SimpleMD's scalar (SHAKE/RATTLE) constraint.>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested).
 *
 * Constrains a plain interatomic distance (the simplest CV with a known, non-singular
 * gradient) on ethane/GFN-FF, NVE, and checks that (a) sigma stays at its target and
 * (b) the conserved energy Epot+Ekin does not drift more than the unconstrained run.
 * This separates the constraint machinery from the path CV, whose level sets are
 * singular where grad s -> 0.
 */
#include <cmath>
#include <iostream>

#include "generated/parameter_registry.h"
#include "src/capabilities/path_cv.h"
#include "src/capabilities/simplemd.h"
#include "src/core/molecule.h"

using curcuma::Molecule;

static struct RegistryInitializer {
    RegistryInitializer() { initialize_generated_registry(); }
} registry_initializer;

struct Result { double sigma_dev, e_span, e_drift; bool ok; };

static json makeCtrl(double dt, int rattle)
{
    json sd;
    sd["method"] = "gfnff"; sd["time_step"] = dt; sd["max_time"] = 1e6;
    sd["temperature"] = 300.0; sd["thermostat"] = "none"; sd["seed"] = 42;
    sd["no_center"] = true; sd["remove_com_mode"] = 0; sd["no_restart"] = true;
    sd["write_xyz"] = false; sd["write_initial_state"] = false;
    sd["write_restart_frequency"] = -1; sd["print_frequency"] = 1000000;
    sd["verbosity"] = 0; sd["unique"] = false; sd["wall_type"] = "none"; sd["threads"] = 1;
    sd["opt"] = false; sd["rmsd_mtd"] = false; sd["mtd"] = false;
    sd["rattle"] = rattle; sd["rattle_12"] = false; sd["rattle_13"] = false;
    json ctrl; ctrl["simplemd"] = sd;
    return ctrl;
}

/*! Potential energy of the molecule at geometry G (fresh SimpleMD, value from prepareRun()). */
static double epotAt(const std::string& xyz, const Geometry& G)
{
    Molecule mol(xyz);
    mol.setGeometry(G);
    SimpleMD md(makeCtrl(0.5, 1), true);
    md.setRestart(false);
    md.setMolecule(mol);
    md.Initialise();
    md.prepareRun();
    return md.potentialEnergy();
}

/*! trueGradient() must be dE/dR without the external force; gradient() must contain it. */
static bool checkTrueGradient(const std::string& xyz)
{
    Molecule mol(xyz);
    SimpleMD md(makeCtrl(0.5, 1), true);
    md.setRestart(false);
    md.setMolecule(mol);
    md.Initialise();
    md.prepareRun();
    Geometry F = Geometry::Zero(mol.AtomCount(), 3);
    for (int a = 0; a < F.rows(); ++a)
        F.row(a) << 0.01 * (a + 1), -0.02, 0.005 * a; // large compared with the true gradient
    for (int s = 0; s < 5; ++s) {
        md.applyExternalForces(F);
        md.step();
    }
    const Geometry R = md.positions();
    const Geometry gt = md.trueGradient();
    const Geometry gp = md.gradient();
    // the integrator's last gradient must still carry F (documents the contamination)
    const double contamination = (gp - gt - F).norm();
    double err = 0, nrm = 0;
    const double h = 1e-4;
    for (int a : { 0, 3 })
        for (int c = 0; c < 3; ++c) {
            Geometry Rp = R, Rm = R;
            Rp(a, c) += h;
            Rm(a, c) -= h;
            const double fd = (epotAt(xyz, Rp) - epotAt(xyz, Rm)) / (2 * h);
            err += (fd - gt(a, c)) * (fd - gt(a, c));
            nrm += fd * fd;
        }
    const double rel = std::sqrt(err / std::max(nrm, 1e-30));
    std::cout << "  trueGradient vs FD of Epot: rel err " << rel << " ; |gradient()-trueGradient()-F| = "
              << contamination << " ; |F| = " << F.norm() << std::endl;
    return rel < 1e-4 && contamination < 1e-8;
}

/*! \brief Claude Generated 2026: the z wall (PathCV::zWallEnergy via setExtraPotential) is part
 *  of the potential the integrator samples.
 *
 *  (1) trueGradient() == d(E + W)/dR by central finite difference of Epot + W, with the wall
 *      active, so NEB-MD's mean-force estimator <b . trueGradient()> averages b . grad(E + W);
 *  (2) NVE with the s-type distance constraint plus the wall: Epot + Ekin + W is conserved as
 *      well as Epot + Ekin without the wall (the wall is evaluated at the new positions inside
 *      Rattle(), not lagged by one step like applyExternalForces()). */
static bool checkExtraPotential(const std::string& xyz)
{
    Molecule mol(xyz);
    const Geometry R0 = mol.getGeometry();
    // two-frame reference path: R0 and R0 with one methyl H moved 0.3 A; z(R0) < 0
    Geometry R1 = R0;
    R1(2, 0) += 0.3;
    PathCV cv({ R0, R1 }, mol, 100.0);
    double z0 = 0.0;
    cv.evaluate(R0, nullptr, &z0, nullptr, nullptr);
    const double zw = z0 + 0.003, kw = 100.0; // narrow tube: thermal motion reaches the wall within fs
    auto wallfn = [&cv, zw, kw](const Geometry& R, Geometry& g) { return cv.zWallEnergy(R, zw, kw, &g); };
    auto dist = [](const Geometry& R, Geometry& g) {
        g = Geometry::Zero(R.rows(), 3);
        const Eigen::RowVector3d d = R.row(2) - R.row(5);
        const double r = d.norm();
        g.row(2) = d / r;
        g.row(5) = -d / r;
        return r;
    };

    // (1) gradient consistency
    bool ok = true;
    {
        SimpleMD md(makeCtrl(0.5, 1), true);
        md.setRestart(false);
        md.setMolecule(mol);
        md.Initialise();
        md.prepareRun();
        md.setExtraPotential(wallfn);
        int s = 0;
        for (; s < 400 && md.extraPotentialEnergy() <= 0.0; ++s)
            md.step();
        for (int t = 0; t < 3; ++t)
            md.step();
        const Geometry R = md.positions();
        Geometry gw;
        const double W = cv.zWallEnergy(R, zw, kw, &gw);
        const Geometry gt = md.trueGradient();
        double err = 0, nrm = 0;
        const double h = 1e-4;
        for (int a = 0; a < R.rows(); ++a)
            for (int c = 0; c < 3; ++c) {
                Geometry Rp = R, Rm = R;
                Rp(a, c) += h;
                Rm(a, c) -= h;
                const double fd = (epotAt(xyz, Rp) + cv.zWallEnergy(Rp, zw, kw, nullptr)
                                   - epotAt(xyz, Rm) - cv.zWallEnergy(Rm, zw, kw, nullptr)) / (2 * h);
                err += (fd - gt(a, c)) * (fd - gt(a, c));
                nrm += fd * fd;
            }
        const double rel = std::sqrt(err / std::max(nrm, 1e-30));
        std::cout << "  z wall: active after " << s << " steps, W = " << W << " Eh (md: "
                  << md.extraPotentialEnergy() << "), |dW/dR| = " << gw.norm()
                  << ", trueGradient vs FD of Epot+W: rel err " << rel << std::endl;
        ok = ok && W > 0.0 && std::abs(W - md.extraPotentialEnergy()) < 1e-14 && rel < 1e-4;
    }
    // (2) NVE conservation with constraint + wall
    auto nve = [&](bool with_wall, double& span, double& drift, double& wmax) {
        SimpleMD md(makeCtrl(0.5, 1), true);
        md.setRestart(false);
        md.setMolecule(mol);
        md.Initialise();
        md.prepareRun();
        Geometry g0;
        md.setScalarConstraint(dist, dist(md.positions(), g0), 1e-10);
        if (with_wall)
            md.setExtraPotential(wallfn);
        double e0 = 0, emin = 1e30, emax = -1e30, elast = 0;
        wmax = 0;
        for (int s = 0; s < 2000; ++s) {
            if (!md.step())
                return false;
            const double e = md.potentialEnergy() + md.kineticEnergy() + md.extraPotentialEnergy();
            wmax = std::max(wmax, md.extraPotentialEnergy());
            if (s == 10) e0 = e;
            if (s >= 10) { emin = std::min(emin, e); emax = std::max(emax, e); elast = e; }
        }
        span = emax - emin;
        drift = elast - e0;
        return true;
    };
    double sf = 0, df = 0, wf = 0, sw = 0, dw = 0, ww = 0;
    const bool r1 = nve(false, sf, df, wf), r2 = nve(true, sw, dw, ww);
    std::cout << "  NVE constrained, 1 ps: no wall span " << sf << " drift " << df << " Eh | wall span " << sw
              << " drift " << dw << " Eh, max W " << ww << " Eh" << std::endl;
    ok = ok && r1 && r2 && ww > 0.0 && std::abs(dw) < 10.0 * std::max(std::abs(df), 1e-5);
    return ok;
}

static Result run(const std::string& xyz, int a, int b, bool constrained, int nsteps, double dt, int rattle)
{
    json sd;
    sd["method"] = "gfnff"; sd["time_step"] = dt; sd["max_time"] = 1e6;
    sd["temperature"] = 300.0; sd["thermostat"] = "none"; sd["seed"] = 42;
    sd["no_center"] = true; sd["remove_com_mode"] = 0; sd["no_restart"] = true;
    sd["write_xyz"] = false; sd["write_initial_state"] = false;
    sd["write_restart_frequency"] = -1; sd["print_frequency"] = 1000000;
    sd["verbosity"] = 0; sd["unique"] = false; sd["wall_type"] = "none"; sd["threads"] = 1;
    sd["opt"] = false; sd["rmsd_mtd"] = false; sd["mtd"] = false;
    sd["rattle"] = rattle; sd["rattle_12"] = false; sd["rattle_13"] = false;
    json ctrl; ctrl["simplemd"] = sd;

    Molecule mol(xyz);
    SimpleMD md(ctrl, true);
    md.setRestart(false);
    md.setMolecule(mol);
    if (!md.Initialise())
        return { 0, 0, 0, false };
    md.prepareRun();

    auto dist = [&](const Geometry& R, Geometry& g) {
        g = Geometry::Zero(R.rows(), 3);
        const Eigen::RowVector3d d = R.row(a) - R.row(b);
        const double r = d.norm();
        g.row(a) = d / r;
        g.row(b) = -d / r;
        return r;
    };
    Geometry g0;
    const double target = dist(md.positions(), g0);
    if (constrained)
        md.setScalarConstraint(dist, target, 1e-10);

    double e0 = 0, emin = 1e30, emax = -1e30, dev = 0, elast = 0;
    for (int s = 0; s < nsteps; ++s) {
        if (!md.step())
            return { dev, 0, 0, false };
        const double e = md.potentialEnergy() + md.kineticEnergy();
        if (s == 10) e0 = e;
        if (s >= 10) { emin = std::min(emin, e); emax = std::max(emax, e); elast = e; }
        Geometry g;
        dev = std::max(dev, std::abs(dist(md.positions(), g) - target));
    }
    return { dev, emax - emin, elast - e0, true };
}

/*! \brief Claude Generated 2026: CSVR temperature of a (constrained) ethane run.
 *
 *  Returns <T> over `nsteps` (first 10 % discarded) with a 5-block standard error and the
 *  number of degrees of freedom SimpleMD uses. The constraint (H-H distance 2-5) removes one
 *  DOF; degreesOfFreedom() must reflect it (24 -> 23 with remove_com_mode 0). */
struct TResult { double meanT, seT; int dof; bool ok; };
static TResult csvrTemperature(const std::string& xyz, bool constrained, int nsteps)
{
    json ctrl = makeCtrl(0.5, constrained ? 1 : 0);
    ctrl["simplemd"]["thermostat"] = "csvr";
    Molecule mol(xyz);
    SimpleMD md(ctrl, true);
    md.setRestart(false);
    md.setMolecule(mol);
    if (!md.Initialise())
        return { 0, 0, 0, false };
    md.prepareRun();
    auto dist = [](const Geometry& R, Geometry& g) {
        g = Geometry::Zero(R.rows(), 3);
        const Eigen::RowVector3d d = R.row(2) - R.row(5);
        const double r = d.norm();
        g.row(2) = d / r;
        g.row(5) = -d / r;
        return r;
    };
    if (constrained) {
        Geometry g0;
        md.setScalarConstraint(dist, dist(md.positions(), g0), 1e-10);
    }
    const int skip = nsteps / 10;
    const int nb = 5;
    std::vector<double> bsum(nb, 0.0);
    std::vector<long> bn(nb, 0);
    for (int s = 0; s < nsteps; ++s) {
        if (!md.step())
            return { 0, 0, md.degreesOfFreedom(), false };
        if (s < skip)
            continue;
        const int b = std::min(nb - 1, (s - skip) * nb / (nsteps - skip));
        bsum[b] += md.currentTemperature();
        ++bn[b];
    }
    double m = 0, m2 = 0;
    for (int b = 0; b < nb; ++b) {
        const double v = bsum[b] / std::max<long>(bn[b], 1);
        m += v;
        m2 += v * v;
    }
    m /= nb;
    const double se = std::sqrt(std::max(m2 / nb - m * m, 0.0) / (nb - 1));
    return { m, se, md.degreesOfFreedom(), true };
}

int main()
{
    const std::string xyz = std::string(CURCUMA_SOURCE_DIR) + "/test_cases/cli/nebmd/01_gfnff_ethane/start.xyz";
    int failures = 0;
    auto check = [&](bool ok, const std::string& w) {
        std::cout << (ok ? "  PASS  " : "  FAIL  ") << w << std::endl;
        if (!ok) ++failures;
    };
    const int N = 2000;
    // atoms 2 and 5: H on different carbons (non-bonded distance, soft DOF, no bond stiffness)
    // 0-2: C-H distance (unequal masses, short distance = strongly curved constraint surface)
    // 0-1: C-C distance (equal masses, short distance)
    // The short-distance pairs caught the Oct 2026 SHAKE-direction bug (direction re-evaluated
    // at each iterate instead of grad(sigma)(R(t))): C-H drift was ~-6e-4 Eh, C-C ~+1e-4 Eh,
    // ~dt^2, while the long H...H distance looked fine. Unequal masses also make the SHAKE
    // displacement not parallel to grad(sigma), so the dv = dr/dt bookkeeping is exercised.
    for (double dt : { 0.5, 0.25 }) {
        for (const auto& ab : { std::pair<int, int>(2, 5), std::pair<int, int>(0, 2), std::pair<int, int>(0, 1) }) {
        std::cout << "  distance atoms " << ab.first << "-" << ab.second << "\n";
        const Result f = run(xyz, ab.first, ab.second, false, N, dt, 1);
        const Result c = run(xyz, ab.first, ab.second, true, N, dt, 1);
        std::cout << "  dt=" << dt << " fs  free: span " << f.e_span << " drift " << f.e_drift
                  << " Eh | constrained: span " << c.e_span << " drift " << c.e_drift
                  << " Eh, max |sigma-target| " << c.sigma_dev << " A\n";
        check(f.ok && c.ok, "runs complete (no instability)");
        check(c.sigma_dev < 1e-6, "distance held to 1e-6 A");
        const bool drift_ok = std::abs(c.e_drift) < 10.0 * std::max(std::abs(f.e_drift), 1e-5);
        check(drift_ok, "Epot+Ekin drift <= 10x the free run (floor 1e-5 Eh)");
        }
    }
    check(checkTrueGradient(xyz), "trueGradient() == dE/dR (FD), gradient() carries the external force");
    check(checkExtraPotential(xyz), "z wall: trueGradient() == d(E+W)/dR (FD); NVE Epot+Ekin+W drift <= 10x no-wall run");
    // Degrees of freedom of the constrained ensemble (Claude Generated 2026). The DOF
    // bookkeeping is checked exactly; <T> is reported (20 ps CSVR, 300 K, one seed) but not
    // asserted - see the AP0b note in docs: the CSVR noise term carries a separate bias.
    {
        const TResult tf = csvrTemperature(xyz, false, 40000);
        const TResult tc = csvrTemperature(xyz, true, 40000);
        std::cout << "  CSVR 300 K, 20 ps: free <T> = " << tf.meanT << " +- " << tf.seT << " K (dof "
                  << tf.dof << ") | constrained <T> = " << tc.meanT << " +- " << tc.seT << " K (dof "
                  << tc.dof << ")" << std::endl;
        check(tf.ok && tc.ok, "CSVR runs complete");
        check(tf.dof == 24 && tc.dof == 23, "scalar constraint removes exactly one DOF (24 -> 23)");
    }
    return failures == 0 ? 0 : 1;
}
