/*
 * <Native KS-DFT Method Wrapper>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * ComputationalMethod wrapper for the native KS-DFT engine. Delegates to DFT.
 * The functional is passed as a DFTFunctional enum (set by the method name in
 * MethodFactory), not as a parameter. hasGradient() is false until WP8.
 *
 * Claude Generated: Native KS-DFT method wrapper for MethodFactory integration
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "../computational_method.h"
#include "dft.h"
#include "src/core/molecule.h"
#include <memory>
#include <string>

class DFTMethod : public ComputationalMethod {
public:
    explicit DFTMethod(DFTFunctional functional, const json& config = json{});
    ~DFTMethod() = default;

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
    bool saveToFile(const std::string& filename) const override { return false; }

private:
    std::unique_ptr<DFT> m_dft;
    Mol m_molecule;
    bool m_calculation_done;
    double m_last_energy;
    std::string m_method_name;

    static json getDefaultConfig();
};