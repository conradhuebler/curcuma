/*
 * <Native ab-initio QM Method Wrapper (HF / KS-DFT)>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * ComputationalMethod wrapper for the native QM engine (QMEngine). The level
 * (hf/lda/pbe/b3lyp) is passed as a QMFunctional enum set by the method name in
 * MethodFactory, not as a parameter. Engine settings come from the `qm` scope
 * (`-qm.basis`, `-qm.scf_*`); the pre-Sep-2026 `dft` scope is still merged.
 * hasGradient() is false until WP8.
 *
 * Claude Generated: Native KS-DFT method wrapper for MethodFactory integration,
 *                   renamed DFTMethod -> QMMethod (Sep 2026)
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "../computational_method.h"
#include "qm_engine.h"
#include "src/core/molecule.h"
#include <memory>
#include <string>

class QMMethod : public ComputationalMethod {
public:
    explicit QMMethod(QMFunctional functional, const json& config = json{});
    ~QMMethod() = default;

    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override { return Vector::Zero(0); }
    Position getDipole() const override { return Position::Zero(); }

    std::string getMethodName() const override { return m_method_name; }
    bool isThreadSafe() const override { return true; }

    bool hasGradient() const override { return false; }  // WP8
    void setThreadCount(int threads) override { (void)threads; }
    void setParameters(const json& params) override { (void)params; }
    json getParameters() const override { return json{}; }
    bool hasError() const override { return false; }

    json getEnergyDecomposition() const override;
    bool saveToFile(const std::string& filename) const override { (void)filename; return false; }

    /**
     * @brief Engine configuration from a method config: registry defaults, then
     *        the top level, then the legacy `dft` scope, then the `qm` scope
     *        (later wins). Shared with HF3CMethod. Claude Generated (Sep 2026).
     */
    static json engineConfig(const json& config);

private:
    std::unique_ptr<QMEngine> m_engine;
    Mol m_molecule;
    bool m_calculation_done;
    double m_last_energy;
    std::string m_method_name;
};
