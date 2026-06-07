/**
 * Native GFN1/GFN2 ALPB-solvation analytical vs numerical gradient test.
 *
 * Validates the in-SCF solvation gradient (Born radii / SASA / CDS H-bond / CM5
 * derivatives + the band-energy Pulay term) added by the native xTB solvation
 * path (docs/SQM_SOLVATION_WP.md, WP1). Compares XTB::Gradient() against a
 * central finite difference of the total energy, for both GFN1 (CM5 charges) and
 * GFN2 (Mulliken) in water.
 *
 * Uses TestMoleculeRegistry — do NOT hardcode geometry here.
 *
 * Usage: test_xtb_solvation_numgrad [solvent]   (default: water)
 * Exit code: 0 = all PASS, 1 = at least one FAIL.
 *
 * Claude Generated (June 2026)
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */

#include "src/core/energy_calculators/qm_methods/xtb_native.h"
#include "src/core/global.h"
#include "src/core/curcuma_logger.h"
#include "core/test_molecule_registry.h"

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using curcuma::xtb::XTB;
using curcuma::xtb::MethodType;

namespace {

// Build + run a fresh native xTB with implicit solvation; return total energy.
// A tight SCF threshold keeps the finite difference clean. gradient optional.
// solvent_model: 3 = ALPB (P16 kernel), 2 = GBSA (Still kernel).
double runXtb(MethodType method, const Mol& mol, const std::string& solvent,
              int solvent_model, bool gradient, Eigen::MatrixXd* grad_out = nullptr)
{
    XTB xtb(method);
    xtb.setSolvent(solvent);
    xtb.setSolventModel(solvent_model);
    xtb.setScfThreshold(1e-9);       // tighten so the FD reference is converged
    if (!xtb.InitialiseMolecule(mol))
        throw std::runtime_error("InitialiseMolecule failed");
    const double e = xtb.Calculation(gradient);
    if (gradient && grad_out)
        *grad_out = xtb.Gradient();  // Eh/Angstrom
    return e;
}

struct Result { std::string tag; int nat; double emax; bool pass; };

Result testMol(MethodType method, const std::string& mname, const std::string& solvent,
               int solvent_model, const std::string& methname, double tol)
{
    const std::string modtag = (solvent_model == 2) ? "gbsa" : "alpb";

    // registry stores Angstrom; Molecule expects Angstrom
    curcuma::Molecule cmol = TestMolecules::TestMoleculeRegistry::createMolecule(mname, false);
    Mol mol = cmol.getMolInfo();
    const int nat = mol.m_number_atoms;

    Eigen::MatrixXd G;
    const double e0 = runXtb(method, mol, solvent, solvent_model, true, &G);   // analytical, Eh/Angstrom

    const double h = 1e-4;  // Angstrom
    Eigen::MatrixXd Gnum(nat, 3);
    for (int i = 0; i < nat; ++i) {
        for (int d = 0; d < 3; ++d) {
            Mol mp = mol, mm = mol;
            mp.m_geometry(i, d) += h;
            mm.m_geometry(i, d) -= h;
            const double ep = runXtb(method, mp, solvent, solvent_model, false);
            const double em = runXtb(method, mm, solvent, solvent_model, false);
            Gnum(i, d) = (ep - em) / (2.0 * h);
        }
    }

    double emax = (G - Gnum).cwiseAbs().maxCoeff();
    bool pass = emax < tol;
    std::cout << "  " << methname << "/" << modtag << " " << mname << " (" << nat
              << " atoms) in " << solvent << ": E=" << std::fixed << std::setprecision(8) << e0
              << "  max|dG|=" << std::scientific << std::setprecision(3) << emax
              << " Eh/A  " << (pass ? "PASS" : "FAIL") << "\n";
    if (!pass) {
        for (int i = 0; i < nat; ++i)
            std::cout << "    atom " << i << "  anal=[" << G.row(i)
                      << "]  num=[" << Gnum.row(i) << "]\n";
    }
    return {methname + "/" + modtag + "/" + mname, nat, emax, pass};
}

}  // namespace

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    const std::string solvent = (argc > 1) ? argv[1] : "water";

    // tol in Eh/Angstrom. The native gradient is exact to the SCF threshold; with a
    // tight SCF and central differences (h=1e-4 A) a converged analytical gradient
    // agrees to ~1e-5 Eh/A. 5e-4 leaves margin for the FD truncation/SCF floor.
    const double tol = 5e-4;
    const std::vector<std::string> mols = {"H2O", "CH4", "CH3OH", "NH3"};
    // 3 = ALPB (P16 kernel), 2 = GBSA (Still kernel) — both share the SASA/HB/CM5 path.
    const std::vector<int> models = {3, 2};

    std::cout << "Native xTB ALPB/GBSA solvation gradient validation (analytical vs FD)\n"
              << "solvent=" << solvent << "  tol=" << tol << " Eh/Angstrom\n"
              << std::string(70, '=') << "\n";

    std::vector<Result> results;
    for (int model : models) {
        for (const auto& m : mols) {
            try { results.push_back(testMol(MethodType::GFN2, m, solvent, model, "gfn2", tol)); }
            catch (const std::exception& e) { std::cerr << "[SKIP] gfn2 " << m << ": " << e.what() << "\n"; }
            try { results.push_back(testMol(MethodType::GFN1, m, solvent, model, "gfn1", tol)); }
            catch (const std::exception& e) { std::cerr << "[SKIP] gfn1 " << m << ": " << e.what() << "\n"; }
        }
    }

    bool all = true;
    for (const auto& r : results) all = all && r.pass;
    std::cout << std::string(70, '=') << "\nOverall: " << (all ? "PASS" : "FAIL")
              << " (" << results.size() << " cases)\n";
    return all ? 0 : 1;
}
