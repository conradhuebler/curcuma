/**
 * QMDFF analytical vs. numerical gradient test.
 *
 * Validates the C++ port of xtb's removed `src/qmdff.f90` (module xtb_qmdff,
 * last version at b7dbd36^:src/qmdff.f90). Every energy term of the QMDFF
 * expression is checked by comparing the analytic gradient produced by
 * ForceField/FFWorkspace against a central finite difference of the total
 * energy on the same parameter set.
 *
 * What this test DOES establish:
 *   - the analytic derivatives of bond / angle / torsion / out-of-plane /
 *     dispersion / electrostatics / repulsion / HB / XB are consistent with
 *     their energy expressions (internal consistency).
 *
 * What it does NOT establish:
 *   - agreement with published QMDFF numbers. Curcuma derives its own QMDFF
 *     parameters from the UFF-style topology; the reference parametrisation is
 *     a QM Hessian fit performed by a separate program. See docs/QMDFF.md.
 *
 * Usage:
 *   test_qmdff_gradients                 – run the default molecule set
 *   test_qmdff_gradients molecule.xyz    – test one XYZ file
 *   test_qmdff_gradients --no-charges    – drop electrostatics and HB/XB
 *   test_qmdff_gradients --no-inversions – zero the improper force constant
 *   test_qmdff_gradients --method=uff    – cross-check the shared UFF terms
 *   test_qmdff_gradients --verbose       – print the generated QMDFF term counts
 * The four diagnostic switches exist to localise a failure to a single term
 * group; the ctest invocation uses none of them.
 *
 * Exit code: 0 = all molecules PASS, 1 = at least one FAIL.
 *
 * Claude Generated (August 2026)
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */

#include "src/core/curcuma_logger.h"
#include "src/core/energycalculator.h"
#include "src/core/global.h"
#include "src/core/molecule.h"
#include "core/test_molecule_registry.h"
#include "json.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using json = nlohmann::json;

// Molecules exercising the different term groups.
//   H2O / CH4    : bond + angle only
//   CH3OH        : + torsion, polar (electrostatics)
//   CH3OCH3      : + several torsions, ether oxygen
//   C6H6         : ring torsions, out-of-plane inversions
//   H2O_dimer    : intermolecular NCI + hydrogen bond
static const std::vector<std::string> DEFAULT_MOLECULES = {
    "H2O", "CH4", "CH3OH", "CH3OCH3", "C6H6", "H2O_dimer",
};

namespace {

std::string g_method = "qmdff";
bool g_no_inversions = false;   // diagnostic: zero the UFF improper force constant
int g_verbosity = 0;            // diagnostic: show the generated QMDFF term counts

json makeConfig()
{
    json config = json::object();
    config["method"] = g_method;
    config["verbosity"] = g_verbosity;
    config["threads"] = 1;
    // Parameter caching would reuse a stale .param.json from a previous run and
    // silently change what is being differentiated.
    config["param_file"] = "none";
    config["writeparam"] = "none";
    if (g_no_inversions)
        config["inversion_force"] = 0.0;
    return config;
}

/// Assign crude partial charges so that the electrostatic and HB/XB terms are
/// actually populated. The magnitude is irrelevant for a derivative test — only
/// that the terms are non-zero.
void assignTestCharges(Mol& mol)
{
    const int n = static_cast<int>(mol.m_atoms.size());
    Vector q = Vector::Zero(n);
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        switch (mol.m_atoms[i]) {
        case 1:  q(i) =  0.15; break;  // H
        case 6:  q(i) = -0.10; break;  // C
        case 7:  q(i) = -0.45; break;  // N
        case 8:  q(i) = -0.40; break;  // O
        case 9:  q(i) = -0.25; break;  // F
        case 17: q(i) = -0.15; break;  // Cl
        default: q(i) =  0.00; break;
        }
        sum += q(i);
    }
    // Enforce neutrality so the electrostatics are not dominated by a net charge
    if (n > 0)
        q.array() -= sum / static_cast<double>(n);
    mol.m_partial_charges = q;
}

struct Result {
    std::string name;
    int n_atoms = 0;
    double energy = 0.0;
    double anal_norm = 0.0;
    double num_norm = 0.0;
    double max_diff = 0.0;
    double rms_diff = 0.0;
    int worst_atom = -1;
    int worst_comp = -1;
    bool pass = false;
    bool nan_detected = false;
};

Result runTest(const std::string& name, const curcuma::Molecule& cmol, double step, double tol,
               bool with_charges = true)
{
    Result res;
    res.name = name;
    res.n_atoms = cmol.AtomCount();

    Mol mol = cmol.getMolInfo();
    if (with_charges)
        assignTestCharges(mol);

    const json config = makeConfig();

    // --- evaluation geometry ------------------------------------------------
    // QMDFF takes r0 and theta0 from the structure it is parametrised on, so at
    // that structure every bonded term sits exactly in its minimum (E = 0, g = 0)
    // and the test would be vacuous. Parametrise at the reference geometry, then
    // evaluate away from it — the situation an optimiser or MD run produces.
    Matrix eval_geometry = mol.m_geometry;
    {
        // Deterministic pseudo-random displacement (LCG) so runs are reproducible
        unsigned int seed = 20260808u;
        auto next = [&seed]() {
            seed = seed * 1664525u + 1013904223u;
            return static_cast<double>(seed % 20001u) / 10000.0 - 1.0; // [-1, 1]
        };
        for (int a = 0; a < res.n_atoms; ++a)
            for (int c = 0; c < 3; ++c)
                eval_geometry(a, c) += 0.05 * next();
    }

    // --- analytic gradient -------------------------------------------------
    EnergyCalculator calc(g_method, config);
    calc.setMolecule(mol);
    calc.updateGeometry(eval_geometry);
    res.energy = calc.CalculateEnergy(true);
    const Matrix analytic = calc.Gradient();

    if (!std::isfinite(res.energy) || !analytic.allFinite()) {
        res.nan_detected = true;
        return res;
    }
    res.anal_norm = analytic.norm();

    // --- numerical gradient ------------------------------------------------
    // The parameter set is generated ONCE and only the geometry is displaced —
    // exactly what an optimiser or MD run sees, and exactly what the analytic
    // gradient differentiates. This matters for the dispersion term: QMDFF's
    // ff_ini computes the C6 coefficients once from the initial coordination
    // numbers (qmdff.f90:146-155) and ff_nonb then treats them as constants, so
    // neither the reference nor this port carries a dC6/dCN response term.
    // Regenerating the parameters at every displacement would silently add that
    // term to the numerical side and make a faithful port look wrong.
    Matrix numeric = Matrix::Zero(res.n_atoms, 3);
    for (int a = 0; a < res.n_atoms; ++a) {
        for (int c = 0; c < 3; ++c) {
            Matrix geom = eval_geometry;
            geom(a, c) += step;
            calc.updateGeometry(geom);
            const double ep = calc.CalculateEnergy(false);

            geom = eval_geometry;
            geom(a, c) -= step;
            calc.updateGeometry(geom);
            const double em = calc.CalculateEnergy(false);

            numeric(a, c) = (ep - em) / (2.0 * step);
        }
    }
    calc.updateGeometry(eval_geometry);
    res.num_norm = numeric.norm();

    double sum_sq = 0.0;
    for (int a = 0; a < res.n_atoms; ++a) {
        for (int c = 0; c < 3; ++c) {
            const double d = std::abs(analytic(a, c) - numeric(a, c));
            sum_sq += d * d;
            if (d > res.max_diff) {
                res.max_diff = d;
                res.worst_atom = a;
                res.worst_comp = c;
            }
        }
    }
    res.rms_diff = std::sqrt(sum_sq / static_cast<double>(3 * res.n_atoms));
    res.pass = (res.max_diff < tol);
    return res;
}

} // namespace

int main(int argc, char** argv)
{
    // parsed below; verbosity is applied once the flags are known

    // Step in Angstrom — the QMDFF geometry convention inside FFWorkspace.
    const double step = 1.0e-5;
    const double tol = 1.0e-6; // Eh/Angstrom

    std::vector<std::pair<std::string, curcuma::Molecule>> molecules;

    // --no-charges isolates the charge-independent terms (bonded + dispersion +
    // repulsion) from electrostatics and HB/XB — useful when localising a failure.
    bool with_charges = true;
    std::string file_arg;
    for (int a = 1; a < argc; ++a) {
        const std::string arg = argv[a];
        if (arg == "--no-charges")
            with_charges = false;
        else if (arg == "--verbose")
            g_verbosity = 2;
        else if (arg == "--no-inversions")
            g_no_inversions = true;
        else if (arg.rfind("--method=", 0) == 0)
            g_method = arg.substr(9); // cross-check the shared UFF terms

        else
            file_arg = arg;
    }

    if (!file_arg.empty()) {
        const std::string file = file_arg;
        curcuma::Molecule mol = curcuma::Molecule(file);
        if (mol.AtomCount() == 0) {
            std::cerr << "Could not read " << file << "\n";
            return 1;
        }
        molecules.emplace_back(file, mol);
    } else {
        for (const auto& name : DEFAULT_MOLECULES) {
            if (!TestMolecules::TestMoleculeRegistry::hasMolecule(name)) {
                std::cerr << "[SKIP] " << name << ": not in TestMoleculeRegistry\n";
                continue;
            }
            molecules.emplace_back(name, TestMolecules::TestMoleculeRegistry::createMolecule(name, false));
        }
    }

    CurcumaLogger::set_verbosity(g_verbosity);

    std::cout << "QMDFF analytic vs. numerical gradient (central difference, h = "
              << step << " A, tol = " << tol << " Eh/A)\n";
    std::cout << std::string(96, '-') << "\n";
    std::cout << std::left << std::setw(14) << "molecule"
              << std::right << std::setw(6) << "N"
              << std::setw(16) << "E [Eh]"
              << std::setw(14) << "|g_ana|"
              << std::setw(14) << "|g_num|"
              << std::setw(13) << "max diff"
              << std::setw(13) << "rms diff"
              << "  status\n";
    std::cout << std::string(96, '-') << "\n";

    int failed = 0;
    for (const auto& [name, mol] : molecules) {
        const Result r = runTest(name, mol, step, tol, with_charges);

        std::cout << std::left << std::setw(14) << r.name
                  << std::right << std::setw(6) << r.n_atoms
                  << std::setw(16) << std::fixed << std::setprecision(8) << r.energy
                  << std::setw(14) << std::setprecision(6) << r.anal_norm
                  << std::setw(14) << r.num_norm
                  << std::setw(13) << std::scientific << std::setprecision(3) << r.max_diff
                  << std::setw(13) << r.rms_diff;

        if (r.nan_detected) {
            std::cout << "  NAN\n";
            ++failed;
        } else if (r.pass) {
            std::cout << "  PASS\n";
        } else {
            std::cout << "  FAIL (atom " << r.worst_atom << " comp " << r.worst_comp << ")\n";
            ++failed;
        }
    }

    std::cout << std::string(96, '-') << "\n";
    if (failed == 0) {
        std::cout << "All " << molecules.size() << " molecules PASS\n";
        return 0;
    }
    std::cout << failed << " of " << molecules.size() << " molecules FAIL\n";
    return 1;
}
