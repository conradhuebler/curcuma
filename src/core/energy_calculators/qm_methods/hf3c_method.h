/*
 * <Native HF-3c composite method>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * HF-3c (R. Sure, S. Grimme, J. Comput. Chem. 34, 1672 (2013)) is Hartree-Fock
 * in the minimal MINIX basis plus three atom-pairwise corrections:
 *
 *   E(HF-3c) = E(HF/MINIX) + E_D3(BJ) + E_gCP + E_SRB
 *
 *   - E_D3(BJ): D3 dispersion, Becke-Johnson damping, s6=1, s8=0.8777,
 *               a1=0.4171, a2=2.9149, two-body only   (D3ParameterGenerator)
 *   - E_gCP:    geometric counterpoise, basis-set superposition error of MINIX
 *   - E_SRB:    short-range basis incompleteness ("bas" in ORCA's "gCP+bas")
 *                                                        (gcp.h, both terms)
 *
 * ORCA's `! HF-3c` is `! HF MINIX D3BJ GCP(HF/MINIX) PATOM` and prints every
 * term on its own, which is how each piece here was validated
 * (docs/NATIVE_QM_IMPLEMENTATION.md, test_cases/qm_hf3c).
 *
 * The HF part is the native QMEngine with the basis forced to MINIX; the other
 * `-qm.*` settings (SCF threshold, guess, ...) are honoured.
 *
 * Scope: closed shell, H-Ne (MINIX file + gCP tables). Analytic gradient (Sep 2026):
 * HF (QMEngine, WP8) + D3 (with the CN chain rule) + gCP/SRB, returned in Eh/Angstrom.
 *
 * Claude Generated: native HF-3c composite (Sep 2026)
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "../computational_method.h"
#include "gcp.h"
#include "qm_engine.h"
#include "src/core/energy_calculators/ff_methods/d3param_generator.h"
#include "src/core/molecule.h"

#include <memory>
#include <string>

class HF3CMethod : public ComputationalMethod {
public:
    explicit HF3CMethod(const json& config = json{});
    ~HF3CMethod() = default;

    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override;
    Vector getCharges() const override { return Vector::Zero(m_molecule.AtomCount()); }
    Vector getBondOrders() const override { return Vector::Zero(0); }
    Position getDipole() const override { return Position::Zero(); }

    std::string getMethodName() const override { return "hf-3c"; }
    bool isThreadSafe() const override { return true; }

    bool hasGradient() const override { return true; }
    void setThreadCount(int threads) override { (void)threads; }
    void setParameters(const json& params) override { (void)params; }
    json getParameters() const override { return json{}; }
    bool hasError() const override { return m_error; }

    json getEnergyDecomposition() const override;
    bool saveToFile(const std::string& filename) const override { (void)filename; return false; }

private:
    std::unique_ptr<QMEngine> m_engine;           ///< HF/MINIX
    std::unique_ptr<D3ParameterGenerator> m_d3;   ///< HF-3c D3(BJ) preset
    gcp::Parameters m_gcp_params;                 ///< HF-3c gCP + SRB
    Mol m_molecule;
    bool m_calculation_done = false;
    bool m_error = false;

    Matrix m_gradient_bohr;                 ///< total dE/dR, Eh/Bohr
    std::vector<Matrix> m_gradient_parts;   ///< {HF, D3, gCP+SRB}, Eh/Bohr (diagnostics)

public:
    /// Gradient parts of the last calculateEnergy(true): {HF, D3, gCP+SRB}, Eh/Bohr.
    const std::vector<Matrix>& gradientParts() const { return m_gradient_parts; }

private:
    // Last energy and its parts (Hartree)
    double m_e_hf = 0.0, m_e_d3 = 0.0, m_e_gcp = 0.0, m_e_srb = 0.0, m_e_total = 0.0;
};
