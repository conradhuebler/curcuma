/*
 * <GFN2 CPSCF / Z-vector susceptibility unit test>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Gate 3b-2: validates applyOrbitalHessian and solveZVector.
 *
 * Test (A) — H2O only, ε = 1e-6:
 *   First-order PT δq_sh vs single-diag numerical δq_sh.
 *   Both paths ignore K, agreement to O(ε²/Δε²) expected.
 *   Tolerance 1e-9.
 *
 * Test (B) — H2O and CH4, ε = 1e-5:
 *   Z-vector CG self-consistency: |A·z_sol − rhs| / |rhs| < 1e-8.
 *   (Degeneracy in CH4 only affects test A, not the CG residual.)
 *
 * Claude Generated (Phase 3b-2, May 2026). GPL-3.0.
 */

#include "src/core/energy_calculators/qm_methods/xtb_native.h"
#include "src/core/curcuma_logger.h"
#include "core/test_molecule_registry.h"

#include <Eigen/Dense>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

using namespace curcuma::xtb;

static constexpr double TOL_SUSC  = 1.0e-9;  // test A: 1st-order vs single-diag
static constexpr double TOL_CGRES = 1.0e-8;  // test B: relative CG residual
static constexpr double EPS_A     = 1.0e-6;  // test A perturbation (keep 2nd-order < tol)
static constexpr double EPS_B     = 1.0e-5;  // test B perturbation

struct TestResult {
    std::string name;
    double      max_err;
    double      tol;
    bool        passed;
};

static void print(const TestResult& r)
{
    std::cout << "  [" << (r.passed ? "PASS" : "FAIL") << "] " << r.name
              << "  err=" << std::scientific << std::setprecision(3) << r.max_err
              << "  tol=" << r.tol << "\n";
}

// Build a zero-initialised Potential of the right dimensions.
static Potential makePotential(int nsh, int nat, int nao)
{
    Potential p;
    p.v_sh = Eigen::VectorXd::Zero(nsh);
    p.v_at = Eigen::VectorXd::Zero(nat);
    p.v_ao = Eigen::VectorXd::Zero(nao);
    p.v_dp = Eigen::MatrixXd::Zero(3, nat);
    p.v_qp = Eigen::MatrixXd::Zero(6, nat);
    return p;
}

// ----------------------------------------------------------------
// Run SCF and return (false, ...) on convergence failure.
// ----------------------------------------------------------------
static bool setupXTB(const Mol& mol, XTB& xtb)
{
    CurcumaLogger::set_verbosity(0);
    xtb.InitialiseMolecule(mol);
    xtb.Calculation(false);
    return xtb.isConverged();
}

// ================================================================
// Test A: first-order analytical z_0 vs single-diag numerical δq_sh
// Perturbation ε on shell 0. Applies to non-degenerate molecules only.
// ================================================================
static TestResult testA_firstOrder(const Mol& mol, const std::string& name)
{
    XTB xtb(MethodType::GFN2);
    if (!setupXTB(mol, xtb))
        return {"(A) " + name + " 1st-order vs single-diag", 1.0, TOL_SUSC, false};

    const int nao   = static_cast<int>(xtb.getMOCoefficients().rows());
    const int nocc  = xtb.getNumElectrons() / 2;
    const int nvirt = nao - nocc;
    const int nsh   = static_cast<int>(xtb.getShellCharges().size());
    const int nat   = mol.m_number_atoms;

    const Matrix& S  = xtb.getOverlap();
    const Matrix& F0 = xtb.getFock();

    // Re-diagonalise F0 to obtain a self-consistent (C0, eps0, P0) baseline.
    // The stored m_wfn comes from the DIIS-extrapolated Fock of the last SCF
    // iteration, which differs from the raw m_F by the residual DIIS error
    // (~1e-8). Diagonalising F0 here makes the analytical and numerical paths
    // share the SAME eigenbasis, so the only residual is O(ε²).
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> sol0(F0, S);
    if (sol0.info() != Eigen::Success)
        return {"(A) " + name + " 1st-order vs single-diag", 1.0, TOL_SUSC, false};
    const Matrix  C   = sol0.eigenvectors();
    const Vector  eps = sol0.eigenvalues();
    const Matrix  C_occ  = C.leftCols(nocc);
    const Matrix  C_virt = C.rightCols(nvirt);
    const Matrix  P0 = 2.0 * C_occ * C_occ.transpose();

    Potential dpot = makePotential(nsh, nat, nao);
    dpot.v_sh(0)   = EPS_A;

    const Matrix dF    = xtb.buildFockFromPotential(dpot);
    const Matrix dF_ov = C_occ.transpose() * dF * C_virt;

    // z_0 from first-order PT
    Matrix z_0(nocc, nvirt);
    for (int i = 0; i < nocc; ++i)
        for (int a = 0; a < nvirt; ++a) {
            const double de = eps(nocc + a) - eps(i);
            z_0(i, a) = -dF_ov(i, a) / (de > 1.0e-6 ? de : 1.0e-6);
        }

    const Matrix dP_anal = 2.0 * (C_occ * z_0 * C_virt.transpose()
                          + C_virt * z_0.transpose() * C_occ.transpose());
    Vector dq_sh_anal, dq_at_a; Eigen::MatrixXd ddp_a, dqp_a;
    xtb.mullikenFromDensity(dP_anal, dq_sh_anal, dq_at_a, ddp_a, dqp_a);

    // Numerical: single diag of F0 + δF
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> sol(F0 + dF, S);
    if (sol.info() != Eigen::Success)
        return {"(A) " + name + " 1st-order vs single-diag", 1.0, TOL_SUSC, false};

    const Matrix P1 = 2.0 * sol.eigenvectors().leftCols(nocc)
                    * sol.eigenvectors().leftCols(nocc).transpose();
    Vector dq_sh_num, dq_at_n; Eigen::MatrixXd ddp_n, dqp_n;
    xtb.mullikenFromDensity(P1 - P0, dq_sh_num, dq_at_n, ddp_n, dqp_n);

    const double err = (dq_sh_anal - dq_sh_num).cwiseAbs().maxCoeff();
    return {"(A) " + name + " 1st-order vs single-diag δq_sh",
            err, TOL_SUSC, err < TOL_SUSC};
}

// ================================================================
// Test B: Z-vector CG self-consistency |A·z_sol − rhs| / |rhs|
// ================================================================
static TestResult testB_zVector(const Mol& mol, const std::string& name)
{
    XTB xtb(MethodType::GFN2);
    if (!setupXTB(mol, xtb))
        return {"(B) " + name + " CG residual", 1.0, TOL_CGRES, false};

    const int nao  = static_cast<int>(xtb.getMOCoefficients().rows());
    const int nocc = xtb.getNumElectrons() / 2;
    const int nsh  = static_cast<int>(xtb.getShellCharges().size());
    const int nat  = mol.m_number_atoms;

    Potential dpot = makePotential(nsh, nat, nao);
    dpot.v_sh(0)   = EPS_B;

    const Matrix dF    = xtb.buildFockFromPotential(dpot);
    const Matrix C_occ = xtb.getMOCoefficients().leftCols(nocc);
    const Matrix C_virt = xtb.getMOCoefficients().rightCols(nao - nocc);
    const Matrix rhs   = -(C_occ.transpose() * dF * C_virt);

    const Matrix z_sol    = xtb.solveZVector(rhs);
    const Matrix residual = xtb.applyOrbitalHessian(z_sol) - rhs;
    const double err = residual.norm() / (rhs.norm() + 1.0e-30);
    return {"(B) " + name + " CG residual |A·z−rhs|/|rhs|",
            err, TOL_CGRES, err < TOL_CGRES};
}

// ================================================================
// Test C: FD validation of the full mulliken-source total gradient
// (Gate 3b-3). Central differences of the total GFN2 energy vs the
// analytical gradient with d4_charge_source="mulliken".
// ================================================================
// Analytical and central-difference total gradient for a given D4 charge source.
static void gradPair(const Mol& mol, const std::string& src,
                     Matrix& g_anal, Matrix& g_num)
{
    const int nat = mol.m_number_atoms;
    XTB xtb(MethodType::GFN2);
    xtb.setD4ChargeSource(src);
    xtb.InitialiseMolecule(mol);
    xtb.Calculation(true);
    g_anal = xtb.Gradient();

    const double h = 5.0e-5;  // Angstrom
    g_num.resize(nat, 3);
    for (int i = 0; i < nat; ++i)
        for (int k = 0; k < 3; ++k) {
            Mol mp = mol, mm = mol;
            mp.m_geometry(i, k) += h;
            mm.m_geometry(i, k) -= h;
            XTB xp(MethodType::GFN2), xm(MethodType::GFN2);
            xp.setD4ChargeSource(src); xm.setD4ChargeSource(src);
            xp.InitialiseMolecule(mp); const double Ep = xp.Calculation(false);
            xm.InitialiseMolecule(mm); const double Em = xm.Calculation(false);
            g_num(i, k) = (Ep - Em) / (2.0 * h);
        }
}

// Isolation test: the eeq path is FD-validated, so its (analytical−FD) error is
// the pre-existing GFN2-gradient baseline. The mulliken path shares that exact
// baseline; the ONLY difference is the q-response. Comparing
// (mull_anal − eeq_anal) against (mull_FD − eeq_FD) cancels the baseline and
// isolates the Mulliken charge-response gradient.
static TestResult testC_responseIsolation(const Mol& mol, const std::string& name, double tol)
{
    CurcumaLogger::set_verbosity(0);
    Matrix gm_a, gm_n, ge_a, ge_n;
    gradPair(mol, "mulliken", gm_a, gm_n);
    gradPair(mol, "eeq",      ge_a, ge_n);

    const Matrix d_anal = gm_a - ge_a;   // analytical q-response difference
    const Matrix d_num  = gm_n - ge_n;   // numerical q-response difference
    const double err = (d_anal - d_num).cwiseAbs().maxCoeff();
    if (std::getenv("CPSCF_DUMP")) {
        std::cout << "    --- " << name << " q-response vectors (anal | num | residual) ---\n";
        for (int i = 0; i < d_anal.rows(); ++i)
            for (int k = 0; k < 3; ++k)
                std::cout << "      atom " << i << " " << "xyz"[k] << ":  "
                          << std::setw(13) << d_anal(i,k) << "  "
                          << std::setw(13) << d_num(i,k) << "  "
                          << std::setw(13) << (d_anal(i,k)-d_num(i,k)) << "\n";
    }
    std::cout << "        [" << name << "] |Δanal|=" << std::scientific << std::setprecision(3)
              << d_anal.cwiseAbs().maxCoeff() << "  |Δnum|=" << d_num.cwiseAbs().maxCoeff()
              << "  ratio(a/n)=" << (d_num.cwiseAbs().maxCoeff() > 1e-12
                  ? d_anal.cwiseAbs().maxCoeff()/d_num.cwiseAbs().maxCoeff() : 0.0) << "\n";
    return {"(C) " + name + " response isolation (mull−eeq: anal vs FD)",
            err, tol, err < tol};
}

// ================================================================
// Test D: direct ∂q_A/∂x validation, fully decoupled from D4.
// Calling computeMullikenChargeResponse with dEdq = e_A yields
// ∂q_A/∂R (Eh-free, units 1/Bohr). FD-validate against the Mulliken
// atomic charge q_A. This isolates the CPSCF charge response itself.
// ================================================================
static constexpr double AU_PER_AA = 1.0 / 0.529177210903;  // AA_TO_AU
static TestResult testD_chargeResponse(const Mol& mol, const std::string& name, int targetA)
{
    CurcumaLogger::set_verbosity(0);
    const int nat = mol.m_number_atoms;

    // Tolerance 3e-2: the isotropic charge response (multipole terms gated off,
    // MP_RESPONSE_ENABLED=false) agrees with FD to ~12-25% on polar molecules;
    // this gate guards the validated RHS_SIGN=+1 (a sign regression gives ~150%).
    // Tightening to <1% is the subject of docs/PHASE3B5_MULTIPOLE_RESPONSE_WP.md.
    const double tolD = 3.0e-2;
    const double scf_thr = 1.0e-9;  // tight SCF → low-noise FD charge reference
    XTB xtb(MethodType::GFN2);
    xtb.setScfThreshold(scf_thr);
    xtb.InitialiseMolecule(mol);
    xtb.Calculation(false);
    if (!xtb.isConverged())
        return {"(D) " + name + " dq/dR", 1.0, tolD, false};

    Vector dEdq = Vector::Zero(nat);
    dEdq(targetA) = 1.0;
    Matrix g_anal = Matrix::Zero(nat, 3);
    xtb.computeMullikenChargeResponse(dEdq, g_anal);   // ∂q_A/∂R, 1/Bohr

    const double h = 2.0e-4;  // Angstrom (stable with tight SCF)
    Matrix g_num(nat, 3);
    for (int i = 0; i < nat; ++i)
        for (int k = 0; k < 3; ++k) {
            Mol mp = mol, mm = mol;
            mp.m_geometry(i, k) += h;
            mm.m_geometry(i, k) -= h;
            XTB xp(MethodType::GFN2), xm(MethodType::GFN2);
            xp.setScfThreshold(scf_thr); xm.setScfThreshold(scf_thr);
            xp.InitialiseMolecule(mp); xp.Calculation(false);
            xm.InitialiseMolecule(mm); xm.Calculation(false);
            const double qp = xp.getCharges()(targetA);
            const double qm = xm.getCharges()(targetA);
            g_num(i, k) = (qp - qm) / (2.0 * h) / AU_PER_AA;  // 1/Bohr
        }

    if (std::getenv("CPSCF_DUMP")) {
        std::cout << "    --- " << name << " dq_" << targetA
                  << "/dR (anal | num | res) ---\n";
        for (int i = 0; i < nat; ++i)
            for (int k = 0; k < 3; ++k)
                std::cout << "      " << i << " " << "xyz"[k] << ":  "
                          << std::setw(13) << g_anal(i,k) << "  "
                          << std::setw(13) << g_num(i,k) << "  "
                          << std::setw(13) << (g_anal(i,k)-g_num(i,k)) << "\n";
    }
    const double err = (g_anal - g_num).cwiseAbs().maxCoeff();
    std::cout << "        [" << name << " dq_" << targetA << "/dR] |anal|="
              << std::scientific << std::setprecision(3) << g_anal.cwiseAbs().maxCoeff()
              << "  |num|=" << g_num.cwiseAbs().maxCoeff() << "  err=" << err << "\n";
    return {"(D) " + name + " dq_" + std::to_string(targetA) + "/dR", err, tolD, err < tolD};
}

int main()
{
    std::cout << "=== test_xtb_cpscf: Gate 3b-2 — orbital Hessian + Z-vector ===\n";

    bool ok = true;

    // --- Test D: direct charge-response validation (gates RHS_SIGN) ---
    for (const char* nm : {"H2O", "HCN", "CH4"}) {
        curcuma::Molecule m =
            TestMolecules::TestMoleculeRegistry::createMolecule(nm, false);
        const TestResult r = testD_chargeResponse(m.getMolInfo(), nm, 0);
        print(r); ok &= r.passed;
    }

    // --- Test A: H2O only (no orbital degeneracy) ---
    {
        curcuma::Molecule m =
            TestMolecules::TestMoleculeRegistry::createMolecule("H2O", false);
        const TestResult r = testA_firstOrder(m.getMolInfo(), "H2O");
        print(r); ok &= r.passed;
    }

    // --- Test B: H2O and CH4 ---
    for (const char* name : {"H2O", "CH4"}) {
        curcuma::Molecule m =
            TestMolecules::TestMoleculeRegistry::createMolecule(name, false);
        const TestResult r = testB_zVector(m.getMolInfo(), name);
        print(r); ok &= r.passed;
    }

    // --- Test C (Gate 3b-3): Mulliken charge-response.
    //   GATED: H2 — homonuclear, q≡0 by symmetry → the response must contribute
    //   ≈0, i.e. the mulliken and eeq analytical gradients must coincide. This
    //   confirms the response routine does not corrupt the (already FD-validated)
    //   baseline. Multipole-free, so no deferred-term confound.
    {
        curcuma::Molecule m =
            TestMolecules::TestMoleculeRegistry::createMolecule("H2", false);
        Matrix gm_a, gm_n, ge_a, ge_n;
        gradPair(m.getMolInfo(), "mulliken", gm_a, gm_n);
        gradPair(m.getMolInfo(), "eeq",      ge_a, ge_n);
        const double err = (gm_a - ge_a).cwiseAbs().maxCoeff();
        const TestResult r{"(C) H2 response ≈ 0 (mulliken == eeq, symmetry)",
                           err, 1.0e-6, err < 1.0e-6};
        print(r); ok &= r.passed;
    }
    //   GATED (<5e-5): with RHS_SIGN=+1 the D4-mulliken charge-response gradient
    //   meets the target for H2O/HCN even with multipole terms off (the residual
    //   in the raw ∂q/∂x is immaterial here because ∂E_D4/∂q is small). The raw
    //   ∂q/∂x accuracy is gated separately by Test D.
    for (const char* name : {"H2O", "HCN"}) {
        curcuma::Molecule m =
            TestMolecules::TestMoleculeRegistry::createMolecule(name, false);
        const TestResult r = testC_responseIsolation(m.getMolInfo(), name, 5.0e-5);
        print(r); ok &= r.passed;
    }

    std::cout << "\n=== " << (ok ? "ALL PASSED" : "SOME FAILED") << " ===\n";
    return ok ? 0 : 1;
}
