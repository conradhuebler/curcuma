/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Pressure terms of the periodic GFN-FF that belong to sharp cutoffs (docs/GPU_PERIODIC_GAPS.md A3) - Claude Generated
 * (Nov 2026). CPU path; the GPU hybrid gives the same virial (gfnff_gpu_pbc_virial).
 *
 * The N1 virial is -dE/d(eps) at a FIXED pair/triple set. A list with a sharp radius (HB/XB search radii and the
 * |E_HB| candidate filter, D4 evaluation radius, repulsion 20 Bohr) loses or gains members under a strain; when the
 * member's energy at the radius is not zero this adds an impulsive term to the pressure of the truncated model, which
 * a fixed-list derivative cannot see (V-LRC derives and adds it for the D4 tail). Here it is measured, per term:
 *   P_thermo,t = -(dE_t/d eps)/(3V) from a least-squares fit of E_t(eps) over isotropic strains with the lists
 *   re-derived at every strain (fresh instance), against P_vir,t = tr Xi_t/(3V) at eps = 0.
 * Also the radius dependence of the pressure: HB search radii scaled (hb_thr1_bohr2, hb_thr2_bohr2), candidate filter
 * off, repulsion lists unchanged (no PARAM).
 *
 * ctest mode (no arguments): 60 waters in L = 44 Bohr, checks only that the decomposition is consistent (sum of terms
 * = total, both for the energies and the virial; the fit residual is reported). Analysis mode:
 *   test_pbc_pressure_cutoffs --xyz FRAME.xyz [--threads N] [--strain] [--hb f1,f2,...]
 * Exit 0 when every check passes.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/fileiterator.h"
#include "src/core/molecule.h"
#include "src/core/units.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using json = nlohmann::json;

static int g_failures = 0;
static void check(bool ok, const std::string& what, double value, double limit)
{
    std::printf("  [%s] %-78s %.3e (limit %.1e)\n", ok ? "PASS" : "FAIL", what.c_str(), value, limit);
    if (!ok)
        ++g_failures;
}

static const double A2B = CurcumaUnit::Length::ANGSTROM_TO_BOHR;

enum Term { T_BOND, T_ANGLE, T_TORS, T_DISP, T_REP, T_COUL, T_HB, T_XB, T_BATM, T_TOTAL, NT };
static const char* kName[NT] = { "bond", "angle", "torsion+inv", "D4 pair+tail", "repulsion", "Coulomb", "hbond", "xbond", "batm", "TOTAL" };

struct Eval {
    bool ok = false;
    double e[NT] = {};
    double xi_tr[NT] = {};
    double lrc_tr = 0.0;
    int n_hb = 0;
};

static Eval evaluate(const Mol& m, const json& cfg, bool virial)
{
    Eval r;
    try {
        GFNFF ff(cfg);
        if (virial)
            ff.setWantVirial(true);
        if (!ff.InitialiseMolecule(m))
            return r;
        r.e[T_TOTAL] = ff.Calculation(virial);
        r.e[T_BOND] = ff.BondEnergy();
        r.e[T_ANGLE] = ff.AngleEnergy();
        r.e[T_TORS] = ff.DihedralEnergy() + ff.InversionEnergy() + ff.STorsEnergy();
        r.e[T_DISP] = ff.DispersionEnergy() + ff.DispersionTailEnergy();
        r.e[T_REP] = ff.BondedRepulsionEnergy() + ff.NonbondedRepulsionEnergy();
        r.e[T_COUL] = ff.CoulombEnergy();
        r.e[T_HB] = ff.HydrogenBondEnergy();
        r.e[T_XB] = ff.HalogenBondEnergy();
        r.e[T_BATM] = ff.BatmEnergy();
        r.n_hb = static_cast<int>(ff.getLastHBonds().size());
        if (virial) {
            const FFVirialComponents* v = ff.virialComponents();
            if (!v)
                return r;
            r.xi_tr[T_BOND] = v->bondTotal().trace();
            r.xi_tr[T_ANGLE] = v->angle.trace();
            r.xi_tr[T_TORS] = (v->dihedral + v->inversion + v->stors).trace();
            r.xi_tr[T_DISP] = v->dispersionTotal().trace();
            r.xi_tr[T_REP] = v->repulsion.trace();
            r.xi_tr[T_COUL] = v->coulombTotal().trace();
            r.xi_tr[T_HB] = v->hbond.trace();
            r.xi_tr[T_XB] = v->xbond.trace();
            r.xi_tr[T_BATM] = v->batm.trace();
            r.xi_tr[T_TOTAL] = v->total().trace();
            r.lrc_tr = v->tail_lrc.trace();
        }
        r.ok = true;
    } catch (const std::exception& e) {
        std::printf("  evaluate threw: %s\n", e.what());
    }
    return r;
}

/// The radii an NPT run uses: D4 evaluation radius frozen as GFNFF::freezeCutoffsForCellChanges(0.15) does.
static json config(const Mol& m, int threads)
{
    const double L = m.m_unit_cell(0, 0) * A2B;
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = threads;
    c["cache_topology"] = false;
    c["frag_charge_model"] = "reference";
    c["dispersion_cutoff_bohr"] = std::nextafter(0.5 * L * std::cbrt(1.0 / 1.15), 0.0);
    return c;
}

static Mol strained(const Mol& m, double eps)
{
    Mol s = m;
    s.m_geometry *= (1.0 + eps);
    s.m_unit_cell *= (1.0 + eps);
    return s;
}

static double bar(double eh, double v_bohr3) { return CurcumaUnit::Pressure::au_to_bar(eh / (3.0 * v_bohr3)); }

/// Least squares E = a + b x + c x^2; returns b and the rms residual
static void fitQuadratic(const std::vector<double>& x, const std::vector<double>& y, double& b, double& rms)
{
    Eigen::MatrixXd A(x.size(), 3);
    Eigen::VectorXd Y(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        A(i, 0) = 1.0;
        A(i, 1) = x[i];
        A(i, 2) = x[i] * x[i];
        Y(i) = y[i];
    }
    const Eigen::Vector3d c = A.colPivHouseholderQr().solve(Y);
    b = c(1);
    rms = std::sqrt((A * c - Y).squaredNorm() / x.size());
}

static void strainScan(const Mol& m, const json& cfg, const Eval& base, const std::vector<double>& eps_list)
{
    const double V = m.m_unit_cell.determinant() * std::pow(A2B, 3);
    std::vector<double> xs;
    std::vector<std::vector<double>> ys(NT);
    for (double eps : eps_list) {
        const Eval e = (eps == 0.0) ? base : evaluate(strained(m, eps), cfg, false);
        if (!e.ok) {
            check(false, "strained instance", 0, 1);
            return;
        }
        xs.push_back(eps);
        for (int t = 0; t < NT; ++t)
            ys[t].push_back(e.e[t]);
        std::printf("    eps %+.1e: E %.10f Eh, HB triples %d\n", eps, e.e[T_TOTAL], e.n_hb);
    }
    std::printf("  %-14s %12s %12s %12s %12s\n", "term", "P_vir (bar)", "P_thermo", "crossing", "fit rms Eh");
    for (int t = 0; t < NT; ++t) {
        double b, rms;
        fitQuadratic(xs, ys[t], b, rms);
        const double p_vir = bar(base.xi_tr[t], V);
        const double p_th = bar(-b, V);
        std::printf("  %-14s %12.3f %12.3f %12.3f %12.2e\n", kName[t], p_vir, p_th, p_th - p_vir, rms);
    }
    std::printf("  V-LRC of the D4 tail (in the reported pressure): %.3f bar\n", bar(base.lrc_tr, V));
}

static void hbScan(const Mol& m, json cfg, const Eval& base, const std::vector<double>& factors)
{
    const double V = m.m_unit_cell.determinant() * std::pow(A2B, 3);
    const double L = m.m_unit_cell(0, 0) * A2B;
    const double thr1 = 250.0, thr2 = 450.0;   // hb_accuracy 0.1
    std::printf("  HB radii scan (thr1, thr2 x f; P relative to f = 1, bar):\n  %6s %8s %9s %12s %12s %12s %10s\n", "f", "r2 Bohr", "triples",
        "dP_total", "dP_hbond", "dE_HB (Eh)", "dE/atom");
    for (double f : factors) {
        if (std::sqrt(thr2 * f) >= 0.5 * L - 0.2) {
            std::printf("  %6.2f  skipped (sqrt(thr2) %.2f >= L/2 - 0.2 = %.2f)\n", f, std::sqrt(thr2 * f), 0.5 * L - 0.2);
            continue;
        }
        json c = cfg;
        c["hb_thr1_bohr2"] = thr1 * f;
        c["hb_thr2_bohr2"] = thr2 * f;
        const Eval e = evaluate(m, c, true);
        if (!e.ok)
            continue;
        std::printf("  %6.2f %8.2f %9d %12.3f %12.3f %12.3e %10.2e\n", f, std::sqrt(thr2 * f), e.n_hb,
            bar(e.xi_tr[T_TOTAL] - base.xi_tr[T_TOTAL], V), bar(e.xi_tr[T_HB] - base.xi_tr[T_HB], V), e.e[T_HB] - base.e[T_HB],
            (e.e[T_TOTAL] - base.e[T_TOTAL]) / m.m_number_atoms);
    }
}

static void analyse(const Mol& m, int threads, bool strain, const std::vector<double>& hb, double filter)
{
    json cfg = config(m, threads);
    cfg["hb_min_pair_energy_eh"] = filter;
    const double V = m.m_unit_cell.determinant() * std::pow(A2B, 3);
    std::printf("%d atoms, L %.4f Bohr, V %.1f Bohr^3, D4 radius %.4f Bohr, hb_min_pair_energy_eh %.1e\n", m.m_number_atoms,
        m.m_unit_cell(0, 0) * A2B, V, cfg["dispersion_cutoff_bohr"].get<double>(), filter);
    const Eval base = evaluate(m, cfg, true);
    check(base.ok, "virial at the frame", base.ok ? 1 : 0, 1);
    if (!base.ok)
        return;
    std::printf("  virial pressure per term (bar):");
    for (int t = 0; t < NT; ++t)
        std::printf(" %s %.2f", kName[t], bar(base.xi_tr[t], V));
    std::printf("; V-LRC %.2f\n", bar(base.lrc_tr, V));
    if (strain)
        strainScan(m, cfg, base, { -8e-4, -4e-4, -2e-4, -1e-4, 0.0, 1e-4, 2e-4, 4e-4, 8e-4 });
    if (!hb.empty())
        hbScan(m, cfg, base, hb);
}

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    std::vector<std::string> files;
    int threads = 4;
    bool strain = false;
    double filter = 1e-9;
    std::vector<double> hb;
    for (int a = 1; a < argc; ++a) {
        const std::string s = argv[a];
        if (s == "--xyz") {
            while (a + 1 < argc && argv[a + 1][0] != '-')
                files.push_back(argv[++a]);
        } else if (s == "--threads" && a + 1 < argc)
            threads = std::atoi(argv[++a]);
        else if (s == "--strain")
            strain = true;
        else if (s == "--filter" && a + 1 < argc)
            filter = std::atof(argv[++a]);
        else if (s == "--hb" && a + 1 < argc) {
            std::stringstream ss(argv[++a]);
            std::string tok;
            while (std::getline(ss, tok, ','))
                hb.push_back(std::atof(tok.c_str()));
        }
    }
    if (!files.empty()) {
        for (const auto& f : files) {
            FileIterator it(f, true);
            Mol m = it.Next().getMolInfo();
            std::printf("== %s\n", f.c_str());
            if (!m.m_has_pbc) {
                std::printf("no unit cell\n");
                ++g_failures;
                continue;
            }
            analyse(m, threads, strain, hb, filter);
        }
        std::printf("%s: %d failure(s)\n", g_failures == 0 ? "PASS" : "FAIL", g_failures);
        return g_failures == 0 ? 0 : 1;
    }
    // ctest: consistency of the per-term split on a small random box
    const Mol dimer = TestMolecules::TestMoleculeRegistry::createMolecule("H2O_dimer", false).getMolInfo();
    std::mt19937 rng(20261123);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const double L = 44.0 / A2B;
    Mol m = dimer;
    m.m_atoms.clear();
    m.m_geometry = Geometry(0, 3);
    const Eigen::RowVector3d c0 = dimer.m_geometry.topRows(3).colwise().mean();
    for (int w = 0; w < 60; ++w) {
        for (int tries = 0; tries < 100000; ++tries) {
            Eigen::RowVector3d p(u(rng) * L, u(rng) * L, u(rng) * L);
            Geometry x(3, 3);
            for (int i = 0; i < 3; ++i)
                x.row(i) = dimer.m_geometry.row(i) - c0 + p;
            bool ok = true;
            for (int a = 0; a < 3 && ok; ++a)
                for (int b = 0; b < m.m_geometry.rows(); ++b) {
                    Eigen::RowVector3d d = x.row(a) - m.m_geometry.row(b);
                    for (int k = 0; k < 3; ++k)
                        d[k] -= L * std::rint(d[k] / L);
                    if (d.norm() < 2.6) {
                        ok = false;
                        break;
                    }
                }
            if (!ok)
                continue;
            const int n0 = m.m_geometry.rows();
            m.m_geometry.conservativeResize(n0 + 3, 3);
            m.m_geometry.bottomRows(3) = x;
            for (int i = 0; i < 3; ++i)
                m.m_atoms.push_back(dimer.m_atoms[i]);
            break;
        }
    }
    m.m_number_atoms = m.m_geometry.rows();
    m.m_has_pbc = true;
    m.m_unit_cell = Eigen::Matrix3d::Identity() * L;
    const Eval e = evaluate(m, config(m, 1), true);
    check(e.ok, "60 H2O: virial", e.ok ? 1 : 0, 1);
    if (e.ok) {
        double se = 0.0, sx = 0.0;
        for (int t = 0; t < T_TOTAL; ++t) {
            se += e.e[t];
            sx += e.xi_tr[t];
        }
        check(std::abs(se - e.e[T_TOTAL]) < 1e-9, "60 H2O: sum of term energies = total", std::abs(se - e.e[T_TOTAL]), 1e-9);
        check(std::abs(sx - e.xi_tr[T_TOTAL]) < 1e-9, "60 H2O: sum of term virials = total", std::abs(sx - e.xi_tr[T_TOTAL]), 1e-9);
    }
    std::printf("%s: %d failure(s)\n", g_failures == 0 ? "PASS" : "FAIL", g_failures);
    return g_failures == 0 ? 0 : 1;
}
