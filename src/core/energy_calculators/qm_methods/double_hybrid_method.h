/*
 * <Native double-hybrid density functionals (WP12 stage 1: energies; stage 2: analytic gradient; CPU)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * A double hybrid adds a scaled second-order perturbation energy, evaluated with the Kohn-Sham orbitals and
 * orbital energies of a hybrid SCF, to that SCF energy (S. Grimme, J. Chem. Phys. 124, 034108 (2006); spin-component
 * scaled "DSD" form: S. Kozuch, D. Gruzman, J. M. L. Martin, J. Phys. Chem. C 114, 20801 (2010)):
 *
 *   E = E_KS[a_x HF exchange + sum_k w_k E_xc,k]  +  c_os E_os(KS)  +  c_ss E_ss(KS)  +  E_disp
 *
 *   step 1: closed-shell KS SCF of the hybrid (QMEngine, QMFunctional::DoubleHybrid, F = h + J + V_xc - (a_x/2) K),
 *   step 2: E_os, E_ss = the MP2 pair sums of qm_mp2.h with the KS orbitals and eigenvalues (non-canonical in the
 *           HF sense: the KS Fock operator is not the HF one, the single excitations are dropped, as in every
 *           double-hybrid definition): RI (QMEngine::computeRIMP2, default) or exact integrals (computeMP2),
 *   step 3: E_disp = D3(BJ) for the -d3bj rows (native D3ParameterGenerator, the HF-3c code; s9 = 0), D4 for the
 *           -d4 rows (native D4 of GFN2 with EEQ(2019) charges and ATM, curcuma::dispersion::D4DftDispersion, WP12
 *           stage 4).
 *
 * The functional table (name, a_x, kernels, c_os, c_ss, dispersion and the source of every number) is
 * kDoubleHybridTable in double_hybrid_method.cpp. Methods: b2plyp, b2plyp-d3bj, b2gp-plyp, b2gp-plyp-d3bj, dsd-blyp,
 * dsd-blyp-d3bj, dsd-pbep86, dsd-pbep86-d3bj, dsd-pbeb95-d3bj, revdsd-pbep86, revdsd-pbep86-d4, b2plyp-d4,
 * b2gp-plyp-d4.
 *
 * Gradient (stage 2, Oct 2026): QMEngine::computeDoubleHybridGradient (qm_dh_grad.cpp: KS gradient + PT2 part from
 * the Z-vector equations with the KS orbital Hessian incl. the XC kernel f_xc, RI or canonical PT2, frozen core) plus
 * the D3(BJ) gradient; CPU only. GPU: J/K and the RI-MP2 sum may run on the device through the existing paths; the
 * XC kernels P86/B95 are CPU-only (counted qm_xc fallback, QMEngine::ensureGpuXC), the PT2 gradient is a counted
 * qm_dh_grad fallback.
 *
 * Claude Generated (Oct 2026, WP12 stage 1, docs/DFT_ROADMAP/WP12_doppelhybride.md). AI-generated, machine-tested;
 * human production test pending.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "../computational_method.h"
#include "qm_engine.h"
#include "qm_mp2.h"
#include "src/core/molecule.h"
#include "src/core/parameter_macros.h"

#include <memory>
#include <string>
#include <vector>

class D3ParameterGenerator;
namespace curcuma::dispersion {
class D4DftDispersion;
struct D4DftParams;
}

// Claude Generated (Oct 2026, WP12 stage 1): double-hybrid settings of the qm module (PT2 options are shared with
// mp2/ri-mp2: -qm.frozen_core, -qm.ri_mp2_aux, -qm.gpu_mp2, -qm.mp2_batch_mb, -qm.mp2_max_memory_mb). Since Oct 31 2026
// -qm.frozen_core "auto" (its default) means true for the double hybrids (as ORCA), false for mp2 and the others.
BEGIN_PARAMETER_DEFINITION(qm)
    PARAM(dispersion, String, "auto",
          "Dispersion term of the native DFT methods. Double hybrids: auto = the one in the method name (none for b2plyp, dsd-pbep86, ...; D3(BJ) for the -d3bj names; D4 for the -d4 names) | none (drop it, e.g. to compare the KS + PT2 part with a reference without dispersion) | d3bj (only for a method whose table row has D3(BJ) parameters) | d4 (a -d4 name, or a name whose -d4 row has the same XC/PT2 numbers: b2plyp, b2gp-plyp, revdsd-pbep86). pbe/b3lyp/tpss/scan: auto = none unless the name ends in -d4, d4 adds D4, none drops it. D3(BJ) = the native one of hf-3c, two-body (s9 = 0); D4 = the native D4 of GFN2 with EEQ(2019) charges and ATM (dftd4 bj-eeq-atm, s9 = 1; dispersion/d4_dft.h), CPU.",
          "DoubleHybrid", {})
    PARAM(dh_pt2, String, "ri",
          "Double hybrids: PT2 integrals. ri = density-fitted (ia|jb) with -qm.ri_mp2_aux (default def2-universal-JKFIT = def2/JK, no def2-SVP/C ships with curcuma: a method variant, the label says so) | canonical = exact (ia|jb) from the n^4 AO tensor (as -method mp2).",
          "DoubleHybrid", {})
    // Default vwn5 since Oct 31 2026 (operator decision, WP12 doc; before: pz81). Behaviour change for dsd-pbep86*,
    // revdsd-pbep86*: energies move by 1e-3..6e-3 Eh.
    PARAM(p86_lda, String, "vwn5",
          "Double hybrids with P86 correlation (dsd-pbep86*, revdsd-pbep86*): local correlation under the P86 gradient term. vwn5 (default since Oct 31 2026; libxc GGA_C_P86VWN, what ORCA 6.1 uses, smooth) | pz81 = Perdew-Zunger 1981 (libxc GGA_C_P86 as in Gaussian, Q-Chem, Psi4, PySCF; discontinuous at r_s = 1, so energies jump by ~1e-9 Eh when a grid point crosses rho = 3/(4 pi) and -opt can stall at |g| ~ 3e-5). The revDSD-PBEP86 parameters were fitted with pz81 (Mehta et al. 2023), the default vwn5 uses them as ORCA's revDSD-PBEP86/2021 does. The two differ by 1e-3..6e-3 Eh for def2-SVP (table in the WP12 doc).",
          "DoubleHybrid", {})
END_PARAMETER_DEFINITION

/// Dispersion term of a double-hybrid table row.
enum class DHDispersion { None, D3BJ, D4 };

/**
 * @brief One double hybrid: everything that defines it, with the source of the numbers.
 *
 * E_xc(DFT) = sum terms[k].weight * E_kernel; exact exchange a_x; PT2 scaling c_os, c_ss; dispersion parameters
 * (D3(BJ): s6, a1, s8, a2 in Bohr). D4 rows leave s6..a2 at 0: their parameters live in ONE place,
 * curcuma::dispersion::d4DftParameters() (dispersion/d4_dft.cpp), under the row name without "-d4".
 * Claude Generated (Oct 2026, WP12).
 */
struct DoubleHybridDef {
    // Plain literals and a fixed term array (trivially destructible, like kMetaGGATable in qm_xc.cpp): no -Winline
    // diagnostics from implicit special members of a static table.
    const char* name;            ///< method name (lower case)
    const char* label;           ///< output label
    double a_x;                  ///< fraction of HF exchange
    qmxc::XCTerm terms[2];       ///< DFT exchange and correlation kernels with weights
    double c_os, c_ss;           ///< opposite- / same-spin PT2 scaling
    DHDispersion dispersion;
    double s6, a1, s8, a2;       ///< dispersion parameters (D3(BJ) / D4 rational damping)
    bool uses_p86;               ///< contains P86 correlation (-qm.p86_lda applies)
    const char* source;          ///< where every number comes from
};

/// The table (double_hybrid_method.cpp). Claude Generated (Oct 2026, WP12).
const std::vector<DoubleHybridDef>& doubleHybridTable();
/// Row by method name (case-insensitive); nullptr for an unknown name. Claude Generated (Oct 2026, WP12).
const DoubleHybridDef* findDoubleHybrid(const std::string& name);
/// The KS part of a row as an XCFunctional (P86 kernels mapped to the -qm.p86_lda variant: vwn5 -> P86VWN_C).
/// Claude Generated (Oct 2026, WP12).
qmxc::XCFunctional doubleHybridXC(const DoubleHybridDef& d, bool p86_vwn5);

class DoubleHybridMethod : public ComputationalMethod {
public:
    /// name: a method name of doubleHybridTable(); throws std::invalid_argument for an unknown name or option.
    DoubleHybridMethod(const std::string& name, const json& config = json{});
    ~DoubleHybridMethod() override;

    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    /// E_KS + c_os E_os + c_ss E_ss + E_disp (Eh); gradient = true also forms the analytic gradient (stage 2).
    double calculateEnergy(bool gradient = false) override;

    /// Eh/Angstrom (ComputationalMethod contract) of the last calculateEnergy(true); zero before.
    Matrix getGradient() const override;
    Vector getCharges() const override { return Vector::Zero(m_molecule.AtomCount()); }
    Vector getBondOrders() const override { return Vector::Zero(0); }
    Position getDipole() const override { return Position::Zero(); }

    std::string getMethodName() const override { return m_def->name; }
    bool isThreadSafe() const override { return true; }
    /// WP12 stage 2: analytic gradient for every table row (RI or canonical PT2, frozen core or not); D4 since
    /// stage 4 (Oct 2026).
    bool hasGradient() const override { return true; }
    void setThreadCount(int threads) override { (void)threads; }
    void setParameters(const json& params) override { (void)params; }
    json getParameters() const override { return json{}; }
    bool hasError() const override { return m_error; }
    std::string getErrorMessage() const override { return m_error_message; }

    /// ks, xc, os, ss (unscaled), c_os, c_ss, pt2 (scaled), disp, total, plus the settings used.
    json getEnergyDecomposition() const override;
    bool saveToFile(const std::string& filename) const override { (void)filename; return false; }

    const qmint::MP2Result& pt2Result() const { return m_pt2; }
    const DoubleHybridDef& definition() const { return *m_def; }
    /// Gradient (Eh/Bohr) of the last calculateEnergy(true) and its parts {1e, 2e, E_nn, XC basis, grid response,
    /// dispersion} (PT2 contributions summed into the first five). Claude Generated (WP12 stage 2).
    const Matrix& gradientBohr() const { return m_gradient_bohr; }
    const std::vector<Matrix>& gradientParts() const { return m_gradient_parts; }
    /// The KS engine (validation hooks only, e.g. QMEngine::freezeXCGrid in the finite-difference test).
    QMEngine& engine() { return *m_engine; }
    /// Z-vector iterations and timings of the last gradient.
    const qmint::MP2Result& pt2GradientResult() const { return m_pt2_grad; }

private:
    const DoubleHybridDef* m_def = nullptr;  ///< row of doubleHybridTable() (static storage)
    std::unique_ptr<QMEngine> m_engine;
    std::unique_ptr<D3ParameterGenerator> m_d3;
    std::unique_ptr<curcuma::dispersion::D4DftDispersion> m_d4;  ///< WP12 stage 4
    const curcuma::dispersion::D4DftParams* m_d4_params = nullptr;
    qmint::MP2Options m_options;
    qmint::MP2Result m_pt2;
    qmint::MP2Result m_pt2_grad;          ///< WP12 stage 2: result of computeDoubleHybridGradient
    Matrix m_gradient_bohr;               ///< WP12 stage 2: total gradient (Eh/Bohr)
    std::vector<Matrix> m_gradient_parts;
    DHDispersion m_dispersion = DHDispersion::None;  ///< after -qm.dispersion
    bool m_p86_vwn5 = false;
    Mol m_molecule;
    bool m_error = false;
    std::string m_error_message;
    double m_e_ks = 0.0, m_e_pt2 = 0.0, m_e_disp = 0.0, m_e_total = 0.0;
};
