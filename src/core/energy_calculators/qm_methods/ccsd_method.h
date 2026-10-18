/*
 * <Native canonical closed-shell CCSD method wrapper (WP13)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * ComputationalMethod for the method name
 *   ccsd     E = E(RHF) + E_corr(CCSD)
 * on top of the native QMEngine (RHF, -qm.basis, -qm.* SCF settings). The CCSD part is
 * QMEngine::computeCCSD (qm_ccsd.cpp); its settings are the qm-module parameters in qm_ccsd.h
 * (-qm.cc_max_iter, -qm.cc_threshold, -qm.cc_residual_threshold, -qm.cc_diis, -qm.cc_max_memory_mb)
 * plus -qm.frozen_core and -qm.mp2_max_memory_mb (AO tensor rebuild limit) from qm_mp2.h.
 *
 * Energy only: a gradient request is an error (no CCSD gradient). Closed shell, exact integrals
 * (-qm.ri j|jk is an error), H-Ne (the engine's basis files).
 *
 * Claude Generated (Oct 2026, WP13). AI-generated, machine-tested.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "../computational_method.h"
#include "qm_ccsd.h"
#include "qm_engine.h"
#include "src/core/molecule.h"

#include <memory>
#include <string>

class CCSDMethod : public ComputationalMethod {
public:
    explicit CCSDMethod(const json& config = json{});
    ~CCSDMethod() override = default;

    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override { return Matrix::Zero(m_molecule.AtomCount(), 3); }
    Vector getCharges() const override { return Vector::Zero(m_molecule.AtomCount()); }
    Vector getBondOrders() const override { return Vector::Zero(0); }
    Position getDipole() const override { return Position::Zero(); }

    std::string getMethodName() const override { return "ccsd"; }
    bool isThreadSafe() const override { return true; }

    bool hasGradient() const override { return false; }
    void setThreadCount(int threads) override { (void)threads; }
    void setParameters(const json& params) override { (void)params; }
    json getParameters() const override { return json{}; }
    bool hasError() const override { return m_error; }
    std::string getErrorMessage() const override { return m_error_message; }

    json getEnergyDecomposition() const override;
    bool saveToFile(const std::string& filename) const override { (void)filename; return false; }

    /// Last CCSD result (energies, iterations, T1 diagnostic, timings).
    const qmint::CCSDResult& ccsdResult() const { return m_result; }

private:
    std::unique_ptr<QMEngine> m_engine;
    qmint::CCSDOptions m_options;
    qmint::CCSDResult m_result;
    Mol m_molecule;
    bool m_error = false;
    std::string m_error_message;
    double m_e_hf = 0.0, m_e_total = 0.0;
};
