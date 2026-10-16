/*
 * <Native canonical MP2 method wrapper (WP11)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * ComputationalMethod for the method names
 *   mp2      E = E(RHF) + E_OS + E_SS
 *   scs-mp2  E = E(RHF) + 1.2 E_OS + 1/3 E_SS   (Grimme 2003)
 *   sos-mp2  E = E(RHF) + 1.3 E_OS              (Jung et al. 2004)
 *   ri-mp2, scs-ri-mp2, sos-ri-mp2: the same with density-fitted (ia|jb) (QMEngine::computeRIMP2,
 *            qm_rimp2.cpp; auxiliary basis -qm.ri_mp2_aux, default def2-universal-JKFIT -- a JK fitting
 *            set as correlation fitting basis, a method variant; Claude Generated, Oct 2026)
 * on top of the native QMEngine (RHF, -qm.basis, -qm.* SCF settings). The MP2 part is
 * QMEngine::computeMP2 (qm_mp2.cpp); its settings are the qm-module parameters in qm_mp2.h
 * (-qm.frozen_core, -qm.mp2_max_memory_mb, -qm.mp2_batch_mb, -qm.gpu_mp2).
 *
 * Energy only: there is no MP2 gradient (needs the Z-vector / relaxed density). A gradient
 * request (-opt, -md, -hessian, -sp -gradient) is an error, never a number.
 * Closed shell, H-Ne (the engine's basis files).
 *
 * Claude Generated (Oct 2026, WP11). AI-generated, machine-tested.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "../computational_method.h"
#include "qm_engine.h"
#include "qm_mp2.h"
#include "src/core/molecule.h"

#include <memory>
#include <string>

class MP2Method : public ComputationalMethod {
public:
    /// variant: "mp2" | "scs-mp2" | "sos-mp2" | "ri-mp2" | "scs-ri-mp2" | "sos-ri-mp2"
    MP2Method(const std::string& variant, const json& config = json{});
    ~MP2Method() override = default;

    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override { return Matrix::Zero(m_molecule.AtomCount(), 3); }
    Vector getCharges() const override { return Vector::Zero(m_molecule.AtomCount()); }
    Vector getBondOrders() const override { return Vector::Zero(0); }
    Position getDipole() const override { return Position::Zero(); }

    std::string getMethodName() const override { return m_variant; }
    bool isThreadSafe() const override { return true; }

    bool hasGradient() const override { return false; }
    void setThreadCount(int threads) override { (void)threads; }
    void setParameters(const json& params) override { (void)params; }
    json getParameters() const override { return json{}; }
    bool hasError() const override { return m_error; }
    std::string getErrorMessage() const override { return m_error_message; }

    json getEnergyDecomposition() const override;
    bool saveToFile(const std::string& filename) const override { (void)filename; return false; }

    /// Last MP2 result (pair energies, timings).
    const qmint::MP2Result& mp2Result() const { return m_result; }

private:
    std::string m_variant;
    std::unique_ptr<QMEngine> m_engine;
    qmint::MP2Options m_options;
    qmint::MP2Result m_result;
    Mol m_molecule;
    bool m_error = false;
    std::string m_error_message;
    double m_e_hf = 0.0, m_e_total = 0.0;
};
