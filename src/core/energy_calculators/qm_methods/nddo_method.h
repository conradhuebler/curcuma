/*
 * <Unified NDDO Method Wrapper>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Single ComputationalMethod wrapper for all NDDO methods (MNDO, AM1, PM3, PM6).
 * Replaces pm3_method.h, am1_method.h, mndo_method.h, pm6_method.h.
 *
 * Claude Generated: Unified NDDO method wrapper for MethodFactory integration
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "../computational_method.h"
#include "nddo.h"
#include "src/core/molecule.h"
#include <memory>
#include <string>

class NDDOMethod : public ComputationalMethod {
public:
    explicit NDDOMethod(NDDOMethodType type, const json& config = json{});
    ~NDDOMethod() = default;

    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override { return Vector::Zero(0); }
    Position getDipole() const override { return Position::Zero(); }

    std::string getMethodName() const override { return m_method_name; }
    bool isThreadSafe() const override { return true; }

    bool hasGradient() const override { return true; }
    void setThreadCount(int threads) override { (void)threads; }
    void setParameters(const json& params) override { (void)params; }
    json getParameters() const override { return json{}; }
    bool hasError() const override { return false; }

    json getEnergyDecomposition() const;
    bool saveToFile(const std::string& filename) const override { return false; }

private:
    std::unique_ptr<NDDO> m_nddo;
    Mol m_molecule;
    bool m_calculation_done;
    double m_last_energy;
    std::string m_method_name;

    static json getDefaultConfig();
};
