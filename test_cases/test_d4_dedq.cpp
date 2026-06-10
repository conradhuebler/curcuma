/*
 * < D4 Charge-Derivative (dE_D4/dq) Finite-Difference Validation >
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Validates the first half of the D4 dispersion charge-response chain rule:
 * the per-atom analytical derivative dE_D4/dq_A returned by
 * D4Evaluator::computeEnergyAndGradient(with_dEdq=true).
 *
 * Method: inject a known charge vector via D4ParameterGenerator::
 * setTopologyCharges() (which getZetaCharges() then returns), perturb a single
 * atomic charge by +-delta, and compare the central finite difference of the
 * total D4 energy against the analytical dE_D4/dq_A.
 *
 * This isolates dE_D4/dq from the dq/dx response (validated separately by the
 * full-gradient FD test, test_xtb_gradient). GFN2 D4 parameters
 * (Caldeweyher 2019 BJ: s6=1, s8=2.7, a1=0.52, a2=5.0).
 *
 * Claude Generated (2026, AP ∂q/∂x — Phase 1). GPL-3.0.
 */

#include "src/core/config_manager.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/dispersion/d4_charge_model.h"
#include "src/core/energy_calculators/dispersion/d4_evaluator.h"
#include "src/core/energy_calculators/dispersion/d4param_generator.h"
#include "src/core/global.h"

#include "core/test_molecule_registry.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

static constexpr double AA_TO_AU = 1.0 / 0.529177210903;

// Finite-difference step on the atomic charge (electrons)
static constexpr double FD_DQ = 1.0e-5;

// Pass threshold: max |dEdq_anal - dEdq_num| per atom (Eh per electron)
static constexpr double DEDQ_TOL = 1.0e-6;

// Build a GFN2-parametrised D4 generator + evaluator pair (StandardBJ_D4).
static std::pair<std::unique_ptr<D4ParameterGenerator>,
                 std::unique_ptr<curcuma::dispersion::D4Evaluator>>
makeGfn2D4()
{
    nlohmann::json d4_config;
    d4_config["d4_s6"] = 1.0;
    d4_config["d4_s8"] = 2.7;
    d4_config["d4_a1"] = 0.52;
    d4_config["d4_a2"] = 5.0;
    d4_config["d4_alp"] = 16.0;
    ConfigManager cfg("d4param", d4_config);
    auto gen = std::make_unique<D4ParameterGenerator>(cfg);

    curcuma::dispersion::D4Params p;
    p.s6 = 1.0;
    p.s8 = 2.7;
    p.a1 = 0.52;
    p.a2 = 5.0;
    p.s9 = 5.0;
    p.alpha = 16.0;
    p.damping = curcuma::dispersion::DampingFormula::StandardBJ_D4;
    auto eval = std::make_unique<curcuma::dispersion::D4Evaluator>(gen.get(), p);
    return {std::move(gen), std::move(eval)};
}

// Energy-only evaluation with an injected charge vector.
static double energyWithCharges(D4ParameterGenerator& gen,
                                curcuma::dispersion::D4Evaluator& eval,
                                const std::vector<int>& atoms,
                                const Matrix& geom_bohr,
                                const Vector& charges)
{
    gen.setTopologyCharges(charges);
    Matrix grad;
    Vector dEdCN, dEdq;
    return eval.computeEnergyAndGradient(atoms, geom_bohr,
                                         /*with_gradient=*/false,
                                         grad, dEdCN, dEdq,
                                         /*with_dEdq=*/false);
}

struct MolResult {
    std::string name;
    int nat;
    double max_err;
    bool passed;
};

static MolResult validateMolecule(const std::string& name)
{
    curcuma::Molecule m = TestMolecules::TestMoleculeRegistry::createMolecule(name, false);
    Mol mol = m.getMolInfo();
    const int nat = mol.m_number_atoms;
    const std::vector<int> atoms = mol.m_atoms;
    const Matrix geom_bohr = mol.m_geometry * AA_TO_AU;

    auto [gen, eval] = makeGfn2D4();
    gen->GenerateParameters(atoms, geom_bohr);

    // Use the geometry-dependent EEQ charges as the reference charge vector.
    Vector q = gen->getEEQCharges();
    if (q.size() != nat) {
        q = Vector::Zero(nat);
    }

    // Analytical dEdq at the reference charges.
    gen->setTopologyCharges(q);
    Matrix grad;
    Vector dEdCN, dEdq_anal;
    eval->computeEnergyAndGradient(atoms, geom_bohr,
                                   /*with_gradient=*/true,
                                   grad, dEdCN, dEdq_anal,
                                   /*with_dEdq=*/true);

    // Central FD over each atomic charge.
    MolResult r{name, nat, 0.0, false};
    for (int a = 0; a < nat; ++a) {
        Vector qp = q;
        qp(a) += FD_DQ;
        const double Ep = energyWithCharges(*gen, *eval, atoms, geom_bohr, qp);

        Vector qm = q;
        qm(a) -= FD_DQ;
        const double Em = energyWithCharges(*gen, *eval, atoms, geom_bohr, qm);

        const double num = (Ep - Em) / (2.0 * FD_DQ);
        const double err = std::abs(num - dEdq_anal(a));
        if (err > r.max_err) r.max_err = err;
    }
    r.passed = r.max_err < DEDQ_TOL;
    return r;
}

// ── Phase 2: validate D4ChargeModel's analytical dq/dx (charge response) ──
// Compare the analytical Σ_A dEdq_A·∂q_A/∂R_m against a central FD of the
// single-shot EEQ charges. Uses a fixed pseudo-arbitrary dEdq vector so the
// test exercises every atom's response, not just the physical dispersion one.
static constexpr double FD_DX = 1.0e-5;             // Bohr
static constexpr double DQDX_TOL = 1.0e-6;          // Eh/Bohr

static MolResult validateChargeResponse(const std::string& name)
{
    curcuma::Molecule m = TestMolecules::TestMoleculeRegistry::createMolecule(name, false);
    Mol mol = m.getMolInfo();
    const int nat = mol.m_number_atoms;
    const std::vector<int> atoms = mol.m_atoms;
    const Matrix geom = mol.m_geometry * AA_TO_AU;

    // Fixed arbitrary "dEdq" so every atom contributes.
    Vector dEdq(nat);
    for (int i = 0; i < nat; ++i) dEdq(i) = 0.1 * std::sin(0.7 * (i + 1)) + 0.05;

    // Analytical response at the reference geometry.
    curcuma::dispersion::D4ChargeModel model;
    model.computeCharges(atoms, geom, 0.0);
    Matrix G_anal = Matrix::Zero(nat, 3);
    model.addChargeResponseGradient(dEdq, G_anal);

    // FD reference: dq/dR via central differences, contracted with dEdq.
    MolResult r{name, nat, 0.0, false};
    for (int mm = 0; mm < nat; ++mm) {
        for (int c = 0; c < 3; ++c) {
            Matrix gp = geom; gp(mm, c) += FD_DX;
            curcuma::dispersion::D4ChargeModel mp;
            Vector qp = mp.computeCharges(atoms, gp, 0.0);

            Matrix gm = geom; gm(mm, c) -= FD_DX;
            curcuma::dispersion::D4ChargeModel mmod;
            Vector qm = mmod.computeCharges(atoms, gm, 0.0);

            double gnum = 0.0;
            for (int a = 0; a < nat; ++a)
                gnum += dEdq(a) * (qp(a) - qm(a)) / (2.0 * FD_DX);

            const double err = std::abs(gnum - G_anal(mm, c));
            if (err > r.max_err) r.max_err = err;
        }
    }
    r.passed = r.max_err < DQDX_TOL;
    return r;
}

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);

    std::vector<std::string> mols;
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) mols.emplace_back(argv[i]);
    } else {
        mols = {"H2O", "CH4", "C6H6"};
    }

    std::cout << std::string(60, '=') << "\n"
              << "  D4 dE/dq Finite-Difference Validation (Phase 1)\n"
              << "  FD dq: " << FD_DQ << "   tolerance: " << DEDQ_TOL << " Eh/e\n"
              << std::string(60, '=') << "\n"
              << std::left << std::setw(12) << "Molecule"
              << std::setw(8) << "Atoms"
              << std::setw(16) << "MaxErr(Eh/e)"
              << "Result\n"
              << std::string(60, '-') << "\n";

    bool all_pass = true;
    for (const std::string& name : mols) {
        MolResult r = validateMolecule(name);
        std::cout << std::left << std::setw(12) << r.name
                  << std::setw(8) << r.nat
                  << std::setw(16) << std::scientific << std::setprecision(2) << r.max_err
                  << (r.passed ? "[PASS]" : "[FAIL]") << "\n";
        all_pass = all_pass && r.passed;
    }

    std::cout << std::string(60, '=') << "\n"
              << "  D4 single-shot EEQ dq/dx Validation (Phase 2)\n"
              << "  FD dx: " << FD_DX << " Bohr   tolerance: " << DQDX_TOL << " Eh/Bohr\n"
              << std::string(60, '=') << "\n"
              << std::left << std::setw(12) << "Molecule"
              << std::setw(8) << "Atoms"
              << std::setw(16) << "MaxErr(Eh/B)"
              << "Result\n"
              << std::string(60, '-') << "\n";
    for (const std::string& name : mols) {
        MolResult r = validateChargeResponse(name);
        std::cout << std::left << std::setw(12) << r.name
                  << std::setw(8) << r.nat
                  << std::setw(16) << std::scientific << std::setprecision(2) << r.max_err
                  << (r.passed ? "[PASS]" : "[FAIL]") << "\n";
        all_pass = all_pass && r.passed;
    }

    std::cout << std::string(60, '=') << "\n"
              << (all_pass ? "All tests PASSED.\n" : "Some tests FAILED.\n");
    return all_pass ? 0 : 1;
}
