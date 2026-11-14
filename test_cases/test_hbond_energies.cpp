/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 - per-triple hydrogen-bond energies of GFN-FF (CPU workspace).
 *
 * Checked against what must hold by construction, not against numbers from a run:
 *  - the triple energies sum to the hydrogen-bond term the decomposition reports,
 *  - the list has one entry per triple of getLastHBonds(), in that order,
 *  - nothing is recorded until asked for, and switching it on takes effect at the next evaluation,
 *  - the result does not depend on the thread count,
 *  - a hydrogen-bonded water dimer has a clearly negative triple energy that dies out with distance.
 */

#include "src/core/energycalculator.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/molecule.h"

#include "generated/parameter_registry.h"

#include <algorithm>
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

/// Two waters, O...O distance @p oo along x; the first donates a nearly linear hydrogen bond.
/// @p copies stacks further copies 3.2 A apart in z (more triples, more work for the threads).
static curcuma::Molecule waterDimers(double oo, int copies)
{
    curcuma::Molecule m;
    for (int c = 0; c < copies; ++c) {
        const double z = 3.2 * c;
        m.addPair({ 8, Position(0.0, 0.0, z) });
        m.addPair({ 1, Position(0.9572, 0.0, z) });
        m.addPair({ 1, Position(-0.2400, 0.9266, z) });
        m.addPair({ 8, Position(oo, 0.0, z) });
        m.addPair({ 1, Position(oo + 0.24, 0.93, z) });
        m.addPair({ 1, Position(oo + 0.24, -0.93, z) });
    }
    return m;
}

struct Run {
    std::vector<GFNFFHBondEnergy> energies;
    std::vector<GFNFFHydrogenBond> list;
    double hbondTerm = 0.0;
    bool ok = false;
};

static Run evaluate(const curcuma::Molecule& molecule, int threads, bool record)
{
    Run r;
    json config = { { "verbosity", 0 }, { "threads", threads }, { "gfnff", json::object() } };
    EnergyCalculator ec("gfnff", config);
    ec.setMolecule(molecule.getMolInfo());
    GFNFF* gfnff = ec.Interface() ? ec.Interface()->gfnffInstance() : nullptr;
    if (!gfnff)
        return r;
    if (record)
        gfnff->setRecordHBondEnergies(true);
    const double e = ec.CalculateEnergy(false);
    if (!std::isfinite(e))
        return r;
    r.energies = gfnff->hydrogenBondEnergies();
    r.list = gfnff->getLastHBonds();
    r.hbondTerm = ec.getEnergyDecomposition().value("HBond", 0.0);
    r.ok = true;
    return r;
}

static double sum(const std::vector<GFNFFHBondEnergy>& v)
{
    double s = 0.0;
    for (const auto& t : v)
        s += t.energy;
    return s;
}

int main()
{
    std::printf("hbond energies\n");

    const Run off = evaluate(waterDimers(2.9, 1), 1, false);
    check(off.ok, "GFN-FF runs on the water dimer");
    check(off.energies.empty(), "nothing is recorded unless asked for");

    const Run on = evaluate(waterDimers(2.9, 1), 1, true);
    check(on.ok && !on.energies.empty(), "recording on: one entry per triple");
    check(on.energies.size() == on.list.size(), "as many entries as getLastHBonds()");
    bool sameOrder = on.energies.size() == on.list.size();
    for (std::size_t i = 0; sameOrder && i < on.list.size(); ++i)
        sameOrder = on.energies[i].donor == on.list[i].i && on.energies[i].hydrogen == on.list[i].j
            && on.energies[i].acceptor == on.list[i].k && on.energies[i].case_type == on.list[i].case_type;
    check(sameOrder, "entries are in the order of the list, with the same atoms");
    check(std::fabs(sum(on.energies) - on.hbondTerm) < 1e-12,
        "triple energies sum to the HBond term (" + std::to_string(sum(on.energies)) + " vs " + std::to_string(on.hbondTerm) + " Eh)");
    check(std::fabs(off.hbondTerm - on.hbondTerm) < 1e-14, "recording does not change the energy");

    double lowest = 0.0;
    int lowestDonor = -1, lowestAcceptor = -1;
    for (const auto& t : on.energies)
        if (t.energy < lowest) {
            lowest = t.energy;
            lowestDonor = t.donor;
            lowestAcceptor = t.acceptor;
        }
    check(lowest < -1e-4, "the bonded dimer has a clearly attractive triple (" + std::to_string(lowest) + " Eh)");
    check(lowestDonor == 0 && lowestAcceptor == 3, "it is the donating water's O (atom 0) to the other water's O (atom 3)");

    const Run far = evaluate(waterDimers(6.0, 1), 1, true);
    double farLowest = 0.0;
    for (const auto& t : far.energies)
        farLowest = std::min(farLowest, t.energy);
    check(far.ok && std::fabs(farLowest) < 0.1 * std::fabs(lowest), "stretched to 6 A the triple energy is below a tenth");

    const Run many1 = evaluate(waterDimers(2.9, 6), 1, true);
    const Run many4 = evaluate(waterDimers(2.9, 6), 4, true);
    bool same = many1.ok && many4.ok && many1.energies.size() == many4.energies.size() && !many1.energies.empty();
    for (std::size_t i = 0; same && i < many1.energies.size(); ++i)
        same = std::fabs(many1.energies[i].energy - many4.energies[i].energy) < 1e-14;
    check(same, "six stacked dimers: the same entries with 1 and 4 threads");
    check(many1.ok && std::fabs(sum(many1.energies) - many1.hbondTerm) < 1e-11, "six dimers: the sum still equals the HBond term");

    std::printf("%s\n", g_failed ? "FAILED" : "all passed");
    return g_failed ? 1 : 0;
}
