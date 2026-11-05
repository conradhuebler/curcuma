/*
 * < Periodic MD: NPT ideal-gas test of the barostats (WP-PERIODIC-NONBONDED N2) >
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
 */

/*
 * Claude Generated (Nov 2026, WP-PERIODIC-NONBONDED N2).
 *
 * 64 neon atoms on a 4 x 4 x 4 lattice (spacing 50 Bohr) in a cubic cell of 200 Bohr, GFN-FF periodic (every atom
 * its own neutral fragment: EEQ charges 0, Coulomb 0; D4 at 50 Bohr ~4e-10 Eh per pair; collisions give a second
 * virial coefficient, measured as the virial part of the pressure: 1.1 % of P0 at L = 120 Bohr, ~0.2 % at 200 Bohr,
 * reported and required < 0.5 % of P0; the potential energy of the lattice is reported), CSVR 300 K (tau 100 fs), dt 1 fs, the
 * barostat every step with the compressibility 1/P0 (the volume relaxation time is then tau_p: 50 fs crescale,
 * 200 fs berendsen). Expectations with the full instantaneous pressure (2 K + tr Xi)/(3 V), K the kinetic energy of
 * the dof = 3N - 3 thermostatted momentum components (the cell removes the translation, D1; the COM momentum is 0):
 *   crescale:  the stationary density of the SDE at fixed scaled coordinates is V exp(-(H + P0 V)/kT) (Ito noise
 *              ~ 1/V); the momentum integral at fixed eps gives V^(dof/3) (each thermalised component scales with
 *              V^(1/3)), dV = V deps, hence rho(V) ~ V^(dof/3) exp(-P0 V/kT):
 *              <V> = N_eff kB T / P0 and <dV^2>/<V>^2 = 1/N_eff with N_eff = dof/3 + 1 (= N + 1 for dof = 3N, the
 *              value of Bernetti and Bussi; = N here).
 *   berendsen: <P> -> P0, i.e. <2K/(3V)> = P0: <V> = (dof/3) kB T/P0 (1 + <dV^2>/<V>^2) (Jensen), and the
 *              fluctuations are suppressed (relative variance well below 1/N_eff; what is left is the thermal noise
 *              of the kinetic pressure integrated over tau_p).
 * Statistics: block averages (blocks of 2 ps) after 5 ps of equilibration, criterion |<V> - V_exp| <= 3 SE (block
 * standard error), variance within +-25 % of 1/N_eff for crescale (n_eff ~ T/(2 tau_p) ~ 2000 at 200 ps: 3 sigma =
 * 10 %), below 0.5/N_eff for berendsen; crescale must also lie > 2 SE away from the Berendsen value (power check).
 * Usage: test_pbc_npt_idealgas <berendsen|crescale> [run_ps = 200]
 */

#include "src/capabilities/simplemd.h"
#include "src/core/curcuma_logger.h"
#include "src/core/molecule.h"
#include "src/core/units.h"
#include "json.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using json = nlohmann::json;

#include "generated/parameter_registry.h"
static struct RegistryInitializer {
    RegistryInitializer() { initialize_generated_registry(); }
} registry_initializer;

static int g_failures = 0;
static void check(bool ok, const std::string& what)
{
    std::printf("%s %s\n", ok ? "[PASS]" : "[FAIL]", what.c_str());
    if (!ok)
        ++g_failures;
}

static curcuma::Molecule neonLattice(int n_per_axis, double L_ang)
{
    curcuma::Molecule mol;
    const double d = L_ang / n_per_axis;
    for (int a = 0; a < n_per_axis; ++a)
        for (int b = 0; b < n_per_axis; ++b)
            for (int c = 0; c < n_per_axis; ++c)
                mol.addPair({ 10, Position{ (a + 0.5) * d, (b + 0.5) * d, (c + 0.5) * d } });
    mol.setUnitCell(Eigen::Matrix3d::Identity() * L_ang, true);
    return mol;
}

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    const std::string barostat = argc > 1 ? argv[1] : "crescale";
    const double run_ps = argc > 2 ? std::atof(argv[2]) : 200.0;
    const double tau_p = barostat == "berendsen" ? 200.0 : 50.0;  // berendsen: longer, so that its (biased) fluctuation is clearly below the ensemble one
    const int N = 64;
    const double T = 300.0;
    const double L0_bohr = 200.0;
    const double kT_eh = kb_Eh * T;  // Eh (global.h, the constant SimpleMD uses)
    const double V0 = L0_bohr * L0_bohr * L0_bohr;
    // P0 so that (N + 1) kB T / P0 = V0 (the start cell is the crescale mean)
    const double P0_bar = CurcumaUnit::Pressure::au_to_bar((N + 1) * kT_eh / V0);
    const double L_ang = L0_bohr * CurcumaUnit::Length::BOHR_TO_ANGSTROM;
    std::printf("64 Ne in L = %.2f Bohr, T %.0f K, kB T %.6e Eh, P0 %.4f bar, compressibility %.4f 1/bar\n", L0_bohr, T, kT_eh, P0_bar, 1.0 / P0_bar);

    json controller;
    json& md = controller["simplemd"];
    md["method"] = "gfnff";
    md["no_restart"] = true;
    md["thermostat"] = "csvr";
    md["coupling"] = 100.0;
    md["temperature"] = T;
    md["seed"] = 7;
    md["dt"] = 1.0;
    md["max_time"] = run_ps * 1000.0 + 5000.0;
    md["no_center"] = true;
    md["print_frequency"] = 1000000;
    md["dump_frequency"] = 1000000;
    md["write_restart_frequency"] = 1000000;
    md["write_xyz"] = false;
    md["write_unique"] = false;
    md["verbosity"] = 0;
    md["barostat"] = barostat;
    md["pressure"] = P0_bar;
    md["tau_p"] = tau_p;
    md["compressibility"] = 1.0 / P0_bar;
    md["pcoupl_interval"] = 1;
    // Ewald at fixed absolute parameters; a loose tolerance (the charges are exactly 0) keeps the k set small, the
    // matrix-free pcg replaces the dense N x N_k matrix build
    controller["gfnff"]["ewald_r_cut"] = 60.0;          // fixed radii far below any L/2 the volume fluctuations reach
    controller["gfnff"]["ewald_tolerance"] = 1e-4;      // (L = 200 +- 4 %: a crescale run at L = 120 once hit the 50 Bohr
    controller["gfnff"]["dispersion_cutoff_bohr"] = 40.0;  // D4 limit L > 100 after 75 ps, a legitimate 3.4 sigma volume excursion)
    controller["gfnff"]["eeq_ewald_dense_max_atoms"] = 0;
    controller["gfnff"]["frag_charge_model"] = "reference";
    controller["gfnff"]["cache_topology"] = false;

    SimpleMD sim(controller, true);
    sim.overrideBasename("pbc_npt_idealgas_" + barostat);
    sim.setMolecule(neonLattice(4, L_ang));
    check(sim.Initialise(), "set up (" + barostat + ")");
    if (g_failures)
        return 1;
    sim.prepareRun();
    const int dof = sim.degreesOfFreedom();
    const double epot0 = sim.potentialEnergy();
    std::printf("DOF %d (3N - 3 = %d), E_pot of the lattice %.3e Eh (ideal gas: 0)\n", dof, 3 * N - 3, epot0);
    check(std::abs(epot0) < 1e-5, "Ne-Ne interaction negligible (|E_pot| < 1e-5 Eh)");

    const double equil_fs = 5000.0, block_fs = 2000.0;
    std::vector<double> v_all, v_block;
    std::vector<double> blocks;
    double sum_virial_part = 0.0, sum_kin_part = 0.0;
    long nsamp = 0;
    int steps = 0;
    double t0 = 0.0;
    while (sim.step()) {
        ++steps;
        const double t = sim.currentTime();
        if (t < equil_fs)
            continue;
        if (v_block.empty())
            t0 = t;
        const double v = sim.cellVolumeAng3() / CurcumaUnit::Length::BOHR3_TO_ANGSTROM3;
        v_all.push_back(v);
        v_block.push_back(v);
        // split of the pressure: the virial part tr Xi/(3V) and the kinetic rest (K of the step before any crescale
        // velocity scaling: pressureBar() minus the virial part)
        const double p_vir = CurcumaUnit::Pressure::au_to_bar(sim.virialTraceEh() / (3.0 * v));
        sum_kin_part += sim.pressureBar() - p_vir;
        sum_virial_part += p_vir;
        ++nsamp;
        if (t - t0 >= block_fs) {
            double m = 0.0;
            for (double x : v_block)
                m += x;
            blocks.push_back(m / v_block.size());
            v_block.clear();
        }
    }
    check(sim.stopReason() == SimpleMD::StopReason::MaxTime, "run reached max_time");
    const long n = static_cast<long>(v_all.size());
    double mean = 0.0;
    for (double x : v_all)
        mean += x;
    mean /= n;
    double var = 0.0;
    for (double x : v_all)
        var += (x - mean) * (x - mean);
    var /= (n - 1);
    double bm = 0.0;
    for (double b : blocks)
        bm += b;
    bm /= blocks.size();
    double bvar = 0.0;
    for (double b : blocks)
        bvar += (b - bm) * (b - bm);
    bvar /= (blocks.size() - 1);
    const double se = std::sqrt(bvar / blocks.size());
    const double rel_var = var / (mean * mean);
    const double kT_over_P0 = kT_eh / (P0_bar / CurcumaUnit::Pressure::HARTREE_PER_BOHR3_TO_BAR);  // Bohr^3
    const double n_eff = dof / 3.0 + 1.0;
    const double v_cr = n_eff * kT_over_P0;                                      // crescale: N_eff kT/P0
    const double v_be = (dof / 3.0) * kT_over_P0 * (1.0 + rel_var);              // berendsen: <2K/3V> = P0 with Jensen
    const double v_exp = barostat == "berendsen" ? v_be : v_cr;
    std::printf("%s: %d steps, %ld samples after %.0f fs, %zu blocks of %.0f fs, tau_p %.0f fs\n", barostat.c_str(), steps, n, equil_fs, blocks.size(), block_fs, tau_p);
    std::printf("  <V> %.1f Bohr^3 = %.3f kT/P0 (block SE %.1f = %.3f %%), expected %.1f (%s), deviation %.2f SE; N_eff kT/P0 = %.1f, (dof/3) kT/P0 (1+var) = %.1f, (N+1) kT/P0 = %.1f\n",
        mean, mean / kT_over_P0, se, 100.0 * se / mean, v_exp, barostat == "berendsen" ? "(dof/3) kT/P0 (1 + var_rel)" : "N_eff kT/P0, N_eff = dof/3 + 1",
        (mean - v_exp) / se, v_cr, v_be, (N + 1) * kT_over_P0);
    std::printf("  relative variance %.5f, 1/N_eff = %.5f, ratio %.3f; <P_kin> %.4f bar, <P_virial> %.4e bar, P0 %.4f bar, <L> %.2f Bohr\n",
        rel_var, 1.0 / n_eff, rel_var * n_eff, sum_kin_part / nsamp, sum_virial_part / nsamp, P0_bar, std::cbrt(mean));
    check(std::abs(sum_virial_part / nsamp) < 5e-3 * P0_bar, "virial part of the pressure (collisions, second virial coefficient) < 0.5 % of P0: ideal gas");
    check(std::abs(mean - v_exp) <= 3.0 * se, "<V> within 3 block standard errors of the expectation");
    if (barostat == "berendsen") {
        check(rel_var < 0.5 / n_eff, "berendsen suppresses the volume fluctuation (relative variance < 0.5/N_eff)");
    } else {
        check(std::abs(rel_var * n_eff - 1.0) < 0.25, "crescale relative volume variance = 1/N_eff within 25 %");
        // the Berendsen expectation must be excluded, else the mean test has no power
        check(std::abs(mean - (dof / 3.0) * kT_over_P0) > 2.0 * se, "crescale <V> is not the Berendsen value (dof/3) kT/P0 (> 2 SE, power check)");
    }
    sim.finalizeRun();
    std::printf("%s: %d failure(s)\n", g_failures ? "FAILED" : "OK", g_failures);
    return g_failures ? 1 : 0;
}
