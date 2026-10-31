/*
 * <Native double-hybrid density functionals (WP12 stage 1: energies, CPU)>
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
 *   step 3: E_disp = D3(BJ) for the -d3bj rows (native D3ParameterGenerator, the HF-3c code; s9 = 0).
 *
 * The functional table (name, a_x, kernels, c_os, c_ss, dispersion and the source of every number) is
 * kDoubleHybridTable in double_hybrid_method.cpp. Methods: b2plyp, b2plyp-d3bj, b2gp-plyp, b2gp-plyp-d3bj, dsd-blyp,
 * dsd-blyp-d3bj, dsd-pbep86, dsd-pbep86-d3bj, dsd-pbeb95-d3bj, revdsd-pbep86, revdsd-pbep86-d4.
 *
 * Energy only (stage 1): hasGradient() is false, a gradient request is an error (-opt/-md stop). GPU: J/K and the
 * RI-MP2 sum may run on the device through the existing paths; the XC kernels P86/B95 are CPU-only (counted qm_xc
 * fallback, QMEngine::ensureGpuXC).
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

// Claude Generated (Oct 2026, WP12 stage 1): double-hybrid settings of the qm module (PT2 options are shared with
// mp2/ri-mp2: -qm.frozen_core, -qm.ri_mp2_aux, -qm.gpu_mp2, -qm.mp2_batch_mb, -qm.mp2_max_memory_mb).
BEGIN_PARAMETER_DEFINITION(qm)
    PARAM(dispersion, String, "auto",
          "Double hybrids: dispersion term. auto = the one in the method name (none for b2plyp, dsd-pbep86, ...; D3(BJ) for the -d3bj names; D4 for revdsd-pbep86-d4) | none (drop it, e.g. to compare the KS + PT2 part with a reference without dispersion) | d3bj (only for a method whose table row has D3(BJ) parameters) | d4 (not available: native D4 exists only for GFN2; an error). The D3(BJ) term is the native one of hf-3c, two-body only (s9 = 0).",
          "DoubleHybrid", {})
    PARAM(dh_pt2, String, "ri",
          "Double hybrids: PT2 integrals. ri = density-fitted (ia|jb) with -qm.ri_mp2_aux (default def2-universal-JKFIT = def2/JK, no def2-SVP/C ships with curcuma: a method variant, the label says so) | canonical = exact (ia|jb) from the n^4 AO tensor (as -method mp2).",
          "DoubleHybrid", {})
    PARAM(p86_lda, String, "pz81",
          "Double hybrids with P86 correlation (dsd-pbep86*, revdsd-pbep86*): local correlation under the P86 gradient term. pz81 = Perdew-Zunger 1981 (libxc GGA_C_P86, Perdew 1986, Gaussian, Q-Chem, Psi4; the revDSD parameters were fitted with it) | vwn5 (libxc GGA_C_P86VWN; what ORCA 6.1 uses, it has no PZ81). The two differ by 2.0e-3 Eh (H2O) to 5.8e-3 Eh (benzene) for dsd-pbep86/def2-SVP (WP12 doc).",
          "DoubleHybrid", {})
END_PARAMETER_DEFINITION

/// Dispersion term of a double-hybrid table row.
enum class DHDispersion { None, D3BJ, D4 };

/**
 * @brief One double hybrid: everything that defines it, with the source of the numbers.
 *
 * E_xc(DFT) = sum terms[k].weight * E_kernel; exact exchange a_x; PT2 scaling c_os, c_ss; dispersion parameters
 * (D3(BJ): s6, a1, s8, a2 in Bohr; D4 parameters are listed for the record only, D4 is not implemented).
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
    /// E_KS + c_os E_os + c_ss E_ss + E_disp (Eh). gradient = true is an error (stage 1, no analytic gradient).
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override { return Matrix::Zero(m_molecule.AtomCount(), 3); }
    Vector getCharges() const override { return Vector::Zero(m_molecule.AtomCount()); }
    Vector getBondOrders() const override { return Vector::Zero(0); }
    Position getDipole() const override { return Position::Zero(); }

    std::string getMethodName() const override { return m_def->name; }
    bool isThreadSafe() const override { return true; }
    /// Stage 1: energies only (the double-hybrid gradient needs the Z-vector of the KS response: stage 2).
    bool hasGradient() const override { return false; }
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

private:
    const DoubleHybridDef* m_def = nullptr;  ///< row of doubleHybridTable() (static storage)
    std::unique_ptr<QMEngine> m_engine;
    std::unique_ptr<D3ParameterGenerator> m_d3;
    qmint::MP2Options m_options;
    qmint::MP2Result m_pt2;
    DHDispersion m_dispersion = DHDispersion::None;  ///< after -qm.dispersion
    bool m_p86_vwn5 = false;
    Mol m_molecule;
    bool m_error = false;
    std::string m_error_message;
    double m_e_ks = 0.0, m_e_pt2 = 0.0, m_e_disp = 0.0, m_e_total = 0.0;
};
