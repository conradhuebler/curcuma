/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 - the ionize capability.
 *
 * The system is built here: a methane "solute" in a lattice of 3.1 A-spaced water
 * molecules (10 x 10 x 10 sites, 31 A cube, a few sites left out around the solute).
 * What is checked are properties that follow from the construction (atom counts,
 * minimum distances, solvent molecules left whole), not numbers from an earlier run.
 */

#include "src/capabilities/ionize.h"
#include "src/core/molecule.h"

#include "generated/parameter_registry.h"

#include <cmath>
#include <cstdio>
#include <string>

struct RegistryInitializer {
    RegistryInitializer() { initialize_generated_registry(); }
} registry_initializer;

static int g_failed = 0;

static void check(bool ok, const std::string& what)
{
    std::printf("  %s  %s\n", ok ? "PASS" : "FAIL", what.c_str());
    if (!ok)
        ++g_failed;
}

static const double L = 31.0;

/// Methane at @p centre, water on the lattice except within 2.5 A of the methane carbon.
static Molecule buildBox(const Position& centre, int& waters, bool periodic)
{
    Molecule m;
    const double d = 0.63;   // C-H 1.09 A along the cube diagonals
    m.addPair({ 6, centre });
    const double s[4][3] = { { 1, 1, 1 }, { 1, -1, -1 }, { -1, 1, -1 }, { -1, -1, 1 } };
    for (const auto& h : s)
        m.addPair({ 1, centre + Position(h[0] * d, h[1] * d, h[2] * d) });
    waters = 0;
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            for (int k = 0; k < 10; ++k) {
                const Position o(3.1 * i, 3.1 * j, 3.1 * k);
                if ((o - centre).norm() < 2.5)
                    continue;
                m.addPair({ 8, o });
                m.addPair({ 1, o + Position(0.757, 0.586, 0.0) });
                m.addPair({ 1, o + Position(-0.757, 0.586, 0.0) });
                ++waters;
            }
    if (periodic)
        m.setUnitCell(Eigen::Matrix3d::Identity() * 3.1 * 10, true);
    return m;
}

static double minimumImage(double d, bool periodic)
{
    if (!periodic)
        return d;
    return d - L * std::round(d / L);
}

static double distance(const Position& a, const Position& b, bool periodic)
{
    const Position r = a - b;
    return Position(minimumImage(r(0), periodic), minimumImage(r(1), periodic),
        minimumImage(r(2), periodic)).norm();
}

int main()
{
    std::printf("ionize\n");
    int waters = 0;
    const Position centre(15.5, 15.5, 15.5);
    const Molecule box = buildBox(centre, waters, false);

    curcuma::IonizeOptions opt;
    opt.charge = 3;
    const curcuma::IonizeResult r = curcuma::ionize(box, opt);
    check(r.ok, "3 anions placed: " + r.error);
    check(r.solvent_molecules == waters,
        "solvent molecules found " + std::to_string(r.solvent_molecules) + " == built " + std::to_string(waters));
    check(r.ions_added == 3 && r.ion == "Cl", "three Cl ions");
    check(r.final_charge == 0 && r.molecule.Charge() == 0, "total charge 0 afterwards");
    check(static_cast<int>(r.molecule.AtomCount()) == static_cast<int>(box.AtomCount()) - 9 + 3,
        "atom count: three waters (9 atoms) replaced by three ions");
    check(r.replaced_atoms.size() == 9, "nine removed atoms reported");

    int cl = 0, o = 0, h = 0, c = 0;
    for (int z : r.molecule.Atoms()) {
        cl += z == 17;
        o += z == 8;
        h += z == 1;
        c += z == 6;
    }
    check(cl == 3 && o == waters - 3 && h == 4 + 2 * (waters - 3) && c == 1,
        "remaining waters whole, solute untouched");

    const int n = static_cast<int>(r.molecule.AtomCount());
    bool tail = true;
    for (int i = n - 3; i < n; ++i)
        tail = tail && r.molecule.Atoms()[i] == 17;
    check(tail, "ions are the last atoms");

    bool far = true, apart = true;
    for (int i = n - 3; i < n; ++i) {
        const Position p = r.molecule.Atom(i).second;
        for (int a = 0; a < 5; ++a)
            far = far && (p - r.molecule.Atom(a).second).norm() >= opt.min_solute_distance;
        for (int j = i + 1; j < n; ++j)
            apart = apart && (p - r.molecule.Atom(j).second).norm() >= opt.min_ion_distance;
    }
    check(far, "ions at least min_solute_distance from every solute atom");
    check(apart, "ions at least min_ion_distance from each other");

    const curcuma::IonizeResult again = curcuma::ionize(box, opt);
    bool same = true;
    for (int i = n - 3; i < n; ++i)
        same = same && (again.molecule.Atom(i).second - r.molecule.Atom(i).second).norm() < 1e-12;
    check(same, "same seed gives the same ions");
    opt.seed = 7;
    const curcuma::IonizeResult other = curcuma::ionize(box, opt);
    bool differs = false;
    for (int i = n - 3; i < n; ++i)
        differs = differs || (other.molecule.Atom(i).second - r.molecule.Atom(i).second).norm() > 1e-6;
    check(other.ok && differs, "another seed gives other ions");
    opt.seed = 1;

    // Negative system: cations.
    curcuma::IonizeOptions neg;
    neg.charge = -2;
    const curcuma::IonizeResult cat = curcuma::ionize(box, neg);
    check(cat.ok && cat.ion == "Na" && cat.ions_added == 2 && cat.final_charge == 0, "charge -2: two Na ions");

    // Already neutral: unchanged.
    curcuma::IonizeOptions zero;
    const curcuma::IonizeResult none = curcuma::ionize(box, zero);
    check(none.ok && none.ions_added == 0 && none.molecule.AtomCount() == box.AtomCount(), "neutral input unchanged");

    // More ions than fit.
    curcuma::IonizeOptions many;
    many.charge = 500;
    const curcuma::IonizeResult toomany = curcuma::ionize(box, many);
    check(!toomany.ok && toomany.error.find("fit") != std::string::npos, "too many ions: error, no partial result");

    // No solvent of the given formula.
    curcuma::IonizeOptions nosolv;
    nosolv.charge = 1;
    nosolv.solvent = "CH3OH";
    const curcuma::IonizeResult nos = curcuma::ionize(box, nosolv);
    check(!nos.ok, "no matching solvent: error");

    // Periodic cell: distances by minimum image. Solute in the corner, so waters across
    // the boundary are close to it although their coordinates are far.
    int pw = 0;
    const Molecule pbox = buildBox(Position(0.3, 0.3, 0.3), pw, true);
    curcuma::IonizeOptions popt;
    popt.charge = 4;
    popt.min_solute_distance = 6.0;
    const curcuma::IonizeResult pr = curcuma::ionize(pbox, popt);
    check(pr.ok, "periodic: 4 ions placed: " + pr.error);
    bool pfar = true, papart = true;
    const int pn = static_cast<int>(pr.molecule.AtomCount());
    for (int i = pn - 4; i < pn; ++i) {
        const Position p = pr.molecule.Atom(i).second;
        for (int a = 0; a < 5; ++a)
            pfar = pfar && distance(p, pr.molecule.Atom(a).second, true) >= popt.min_solute_distance - 1e-9;
        for (int j = i + 1; j < pn; ++j)
            papart = papart && distance(p, pr.molecule.Atom(j).second, true) >= popt.min_ion_distance - 1e-9;
    }
    check(pfar, "periodic: ions at least min_solute_distance away by minimum image");
    check(papart, "periodic: ions at least min_ion_distance apart by minimum image");
    check(pr.molecule.hasPBC(), "periodic: cell kept");

    std::printf("%s\n", g_failed ? "FAILED" : "all passed");
    return g_failed ? 1 : 0;
}
