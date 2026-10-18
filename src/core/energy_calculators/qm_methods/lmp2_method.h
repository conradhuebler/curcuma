/*
 * <Native local PNO-RI-MP2 method wrapper (WP15)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * ComputationalMethod for the method names lmp2 (alias local-ri-mp2):
 *   E = E(RHF) + E_corr(local PNO-RI-MP2)
 * on top of the native QMEngine (RHF, -qm.basis, any SCF: exact or -qm.ri j|jk). The correlation part is
 * QMEngine::computeLMP2 (qm_lmp2.cpp); its settings are the -qm.lmp2_* parameters in qm_lmp2.h plus
 * -qm.frozen_core and -qm.ri_mp2_aux (qm_mp2.h). -qm.lmp2_check true also runs canonical RI-MP2 with the same
 * auxiliary basis and prints the difference.
 *
 * Energy only: a gradient request is an error. Closed shell, H-Ne (the engine's basis and aux files).
 *
 * Claude Generated (Oct 2026, WP15). AI-generated, machine-tested.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "../computational_method.h"
#include "qm_engine.h"
#include "qm_lmp2.h"
#include "qm_mp2.h"
#include "src/core/molecule.h"

#include <memory>
#include <string>

class LMP2Method : public ComputationalMethod {
public:
    explicit LMP2Method(const json& config = json{});
    ~LMP2Method() override = default;

    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override { return Matrix::Zero(m_molecule.AtomCount(), 3); }
    Vector getCharges() const override { return Vector::Zero(m_molecule.AtomCount()); }
    Vector getBondOrders() const override { return Vector::Zero(0); }
    Position getDipole() const override { return Position::Zero(); }

    std::string getMethodName() const override { return "lmp2"; }
    bool isThreadSafe() const override { return true; }

    bool hasGradient() const override { return false; }
    void setThreadCount(int threads) override { (void)threads; }
    void setParameters(const json& params) override { (void)params; }
    json getParameters() const override { return json{}; }
    bool hasError() const override { return m_error; }
    std::string getErrorMessage() const override { return m_error_message; }

    json getEnergyDecomposition() const override;
    bool saveToFile(const std::string& filename) const override { (void)filename; return false; }

    const qmint::LMP2Result& lmp2Result() const { return m_result; }

private:
    std::unique_ptr<QMEngine> m_engine;
    qmint::LMP2Options m_options;
    qmint::LMP2Result m_result;
    bool m_check = false;
    qmint::MP2Result m_check_result;
    Mol m_molecule;
    bool m_error = false;
    std::string m_error_message;
    double m_e_hf = 0.0, m_e_total = 0.0;
};
