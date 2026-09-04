/*
 * test_gfnff_react_filters.cpp — react-topology scan bookkeeping
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Sep 2026): pins the bookkeeping rules of the react-mode
 * hysteresis scan (docs/GFNFF_REACT_TOPOLOGY.md) that the MD tests cannot
 * isolate. Every case uses react_check_every=1 so each energy call scans once.
 *
 * A. Refractory period blocks exactly N scans after a break, and while it
 *    blocks, the FORCE FIELD topology has no bond either (the empty react set
 *    is authoritative; no geometric re-detection sneaks a bond back in).
 * B. Exchange resolution keeps the sigma bookkeeping consistent: breaking the
 *    bond between two over-valent atoms resolves BOTH; the partner must not
 *    lose a second bond in the same scan.
 * C. Event record + bond orders: consumeReactEvents() reports the pairs with a
 *    finite dE_jump, and reactiveBondOrders() reports chemical orders. The FT-Hueckel
 *    solver models ONE pi system, so its pi order is ~1.0 for an sp2-sp2 double bond
 *    AND for an sp-sp triple bond; the sp-sp case adds the second, degenerate pi
 *    system back. Pinned here for N2 (3), ethylene C=C (2) and the C-H bonds (1).
 */
#include <cmath>
#include <iomanip>
#include <iostream>

#include "src/core/energy_calculators/qm_methods/gfnff_method.h"
#include "src/core/energycalculator.h"
#include "src/core/molecule.h"

#include "json.hpp"
using json = nlohmann::json;

static GFNFF* gfnffOf(EnergyCalculator& calc)
{
    auto* m = dynamic_cast<GFNFFComputationalMethod*>(calc.Interface());
    return m ? m->getGFNFF() : nullptr;
}

static int g_failed = 0;
static void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << "\n";
    if (!ok)
        ++g_failed;
}

static json reactConfig(json extra)
{
    json gf = { { "topology_mode", "react" }, { "react_check_every", 1 }, { "react_check_disp_bohr", 0.0 } };
    for (auto& [k, v] : extra.items())
        gf[k] = v;
    return { { "verbosity", 0 }, { "threads", 1 }, { "gfnff", gf } };
}

// ---------------------------------------------------------------------------
// A. Refractory + authoritative empty bond set (H2: bond, break, wait, re-form)
// ---------------------------------------------------------------------------
static void caseRefractory()
{
    std::cout << "A. refractory period / empty react set\n";
    const int N = 3;
    curcuma::Molecule mol;
    mol.addPair({ 1, Position(0.00, 0.0, 0.0) });
    mol.addPair({ 1, Position(0.74, 0.0, 0.0) });

    EnergyCalculator calc("gfnff", reactConfig({ { "react_refractory_scans", N } }));
    calc.setMolecule(mol.getMolInfo());
    GFNFF* ff = gfnffOf(calc);
    check(ff != nullptr, "GFNFF instance reachable through the CPU wrapper");
    if (!ff)
        return;

    calc.CalculateEnergy(false); // call 1: seeded bond survives
    check(ff->reactiveBonds().size() == 1, "seeded H-H bond present");

    Matrix far = mol.getGeometry();
    far(1, 0) = 5.0;
    calc.updateGeometry(far);
    calc.CalculateEnergy(false); // call 2: break event, refractory armed
    check(ff->reactiveBonds().empty(), "bond broken at 5.0 A");

    Matrix close = mol.getGeometry();
    calc.updateGeometry(close);
    std::vector<size_t> nbonds;
    double e_blocked = 0.0;
    for (int s = 0; s < N; ++s) {
        e_blocked = calc.CalculateEnergy(false); // calls 3..N+2: blocked
        nbonds.push_back(ff->reactiveBonds().size());
    }
    bool all_blocked = true;
    for (size_t n : nbonds)
        all_blocked = all_blocked && (n == 0);
    check(all_blocked, "no re-formation during the " + std::to_string(N) + " refractory scans");

    // While the react set is empty the force-field topology must be bond-free too.
    const auto& topo = ff->getTopologyInfo();
    const bool topo_bond_free = topo.neighbor_lists.size() == 2
        && topo.neighbor_lists[0].empty() && topo.neighbor_lists[1].empty();
    check(topo_bond_free, "force-field topology carries no bond while the react set is empty");

    const double e_reformed = calc.CalculateEnergy(false); // call N+3: allowed again
    check(ff->reactiveBonds().size() == 1, "bond re-forms on the first scan after the refractory period");

    // The re-formed surface equals a fresh initialisation with that bond (bit-identical).
    Mol forced = mol.getMolInfo();
    forced.m_bonds = { { 0, 1 } };
    json cfg_auto = { { "verbosity", 0 }, { "threads", 1 }, { "gfnff", json::object() } };
    EnergyCalculator ref("gfnff", cfg_auto);
    ref.setMolecule(forced);
    const double e_ref = ref.CalculateEnergy(false);
    std::cout << std::scientific << std::setprecision(6)
              << "    E(blocked, no bond) = " << e_blocked << "  E(re-formed) = " << e_reformed
              << "  E(fresh init) = " << e_ref << "\n";
    check(std::abs(e_reformed - e_ref) < 1e-9, "re-formed surface == fresh init with the bond");
    check(e_blocked - e_reformed > 0.05, "bond-free surface lies far above the bonded one at 0.74 A");

    // Events: break, then formation with a finite dE_jump.
    auto events = ff->consumeReactEvents();
    bool ev_ok = events.size() == 2
        && events[0].formed.empty() && events[0].broken.size() == 1 && events[0].broken[0] == std::make_pair(0, 1)
        && events[1].broken.empty() && events[1].formed.size() == 1 && events[1].formed[0] == std::make_pair(0, 1)
        && std::isfinite(events[1].de_jump_eh) && events[1].de_jump_eh < -0.05;
    check(ev_ok, "event record: one break, one formation with finite negative dE_jump");
    if (events.size() == 2)
        std::cout << "    dE_jump(formation) = " << events[1].de_jump_eh << " Eh\n";
    check(ff->consumeReactEvents().empty(), "events are consumed (second call empty)");
}

// ---------------------------------------------------------------------------
// B. Exchange resolution sigma bookkeeping (H-H-H-H chain, both inner atoms over-valent)
// ---------------------------------------------------------------------------
static void caseExchangeBookkeeping()
{
    std::cout << "B. exchange resolution bookkeeping\n";
    // a-b 0.75, b-c 0.85 (weakest), c-d 0.75; seeded as a chain via forced bonds.
    curcuma::Molecule mol;
    mol.addPair({ 1, Position(0.00, 0.0, 0.0) });
    mol.addPair({ 1, Position(0.75, 0.0, 0.0) });
    mol.addPair({ 1, Position(1.60, 0.0, 0.0) });
    mol.addPair({ 1, Position(2.35, 0.0, 0.0) });
    Mol seeded = mol.getMolInfo();
    seeded.m_bonds = { { 0, 1 }, { 1, 2 }, { 2, 3 } };

    EnergyCalculator calc("gfnff", reactConfig({ { "react_exchange_scans", 1 }, { "react_refractory_scans", 0 } }));
    calc.setMolecule(seeded);
    GFNFF* ff = gfnffOf(calc);
    if (!ff) {
        check(false, "GFNFF instance reachable");
        return;
    }
    calc.CalculateEnergy(false); // scan 1: streak b=1, c=1 (not > 1)
    check(ff->reactiveBonds().size() == 3, "chain intact after the first scan");
    calc.CalculateEnergy(false); // scan 2: b resolves by breaking b-c; c must then be satisfied
    const auto& bonds = ff->reactiveBonds();
    const bool two_left = bonds.size() == 2
        && bonds[0] == std::make_pair(0, 1) && bonds[1] == std::make_pair(2, 3);
    check(two_left, "breaking the shared weakest bond resolves both atoms (a-b and c-d remain)");
    auto events = ff->consumeReactEvents();
    bool ev_ok = events.size() == 1 && events[0].formed.empty()
        && events[0].broken.size() == 1 && events[0].broken[0] == std::make_pair(1, 2);
    check(ev_ok, "exactly one break event (b-c), no pair both formed and broken");
}

// ---------------------------------------------------------------------------
// C. Bond orders from the Hueckel pi orders
// ---------------------------------------------------------------------------
static void caseBondOrders()
{
    std::cout << "C. bond orders\n";
    curcuma::Molecule n2;
    n2.addPair({ 7, Position(0.00, 0.0, 0.0) });
    n2.addPair({ 7, Position(1.10, 0.0, 0.0) });
    Mol seeded = n2.getMolInfo();
    seeded.m_bonds = { { 0, 1 } };
    EnergyCalculator calc("gfnff", reactConfig({}));
    calc.setMolecule(seeded);
    GFNFF* ff = gfnffOf(calc);
    if (!ff) {
        check(false, "GFNFF instance reachable");
        return;
    }
    calc.CalculateEnergy(false);
    const auto& orders = ff->reactiveBondOrders();
    check(orders.size() == ff->reactiveBonds().size(), "bond-order vector parallel to the bond list");
    check(orders.size() == 1 && orders[0] == 3, "N#N reported as order 3 (got "
        + (orders.empty() ? std::string("none") : std::to_string(orders[0])) + ")");
    check(ff->topologyMode() == "react", "topologyMode() reports react");

    // Ethylene: the sp2-sp2 double bond must stay 2 (same pi order as N2's triple
    // bond, so only the sp-sp rule may separate them), and every C-H stays 1.
    curcuma::Molecule eth;
    eth.addPair({ 6, Position(0.00, 0.00, 0.0) });
    eth.addPair({ 6, Position(1.33, 0.00, 0.0) });
    eth.addPair({ 1, Position(-0.55, 0.94, 0.0) });
    eth.addPair({ 1, Position(-0.55, -0.94, 0.0) });
    eth.addPair({ 1, Position(1.88, 0.94, 0.0) });
    eth.addPair({ 1, Position(1.88, -0.94, 0.0) });
    Mol eth_seed = eth.getMolInfo();
    eth_seed.m_bonds = { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 4 }, { 1, 5 } };
    EnergyCalculator calc_eth("gfnff", reactConfig({}));
    calc_eth.setMolecule(eth_seed);
    GFNFF* ff_eth = gfnffOf(calc_eth);
    if (!ff_eth) {
        check(false, "GFNFF instance reachable (ethylene)");
        return;
    }
    calc_eth.CalculateEnergy(false);
    const auto& eo = ff_eth->reactiveBondOrders();
    const auto& eb = ff_eth->reactiveBonds();
    int cc_order = -1, max_ch_order = -1;
    for (size_t k = 0; k < eb.size() && k < eo.size(); ++k) {
        const bool is_cc = eb[k] == std::make_pair(0, 1);
        if (is_cc)
            cc_order = eo[k];
        else
            max_ch_order = std::max(max_ch_order, eo[k]);
    }
    check(cc_order == 2, "ethylene C=C reported as order 2 (got " + std::to_string(cc_order) + ")");
    check(max_ch_order == 1, "ethylene C-H bonds stay order 1 (got " + std::to_string(max_ch_order) + ")");
}

int main()
{
    caseRefractory();
    caseExchangeBookkeeping();
    caseBondOrders();
    std::cout << (g_failed == 0 ? "PASS" : "FAIL") << " (" << g_failed << " failed)\n";
    return g_failed == 0 ? 0 : 1;
}
