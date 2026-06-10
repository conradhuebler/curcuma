/*
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 *
 * Claude Generated: PM6 method wrapper for MethodFactory integration
 */

#pragma once

#include "../computational_method.h"
#include "pm6.h"
#include "src/core/molecule.h"
#include <memory>

class PM6Method : public ComputationalMethod {
public:
    explicit PM6Method(const json& config = json{});
    ~PM6Method() = default;

    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override { return Vector::Zero(0); }
    Position getDipole() const override { return Position::Zero(); }

    std::string getMethodName() const override { return "pm6"; }
    bool isThreadSafe() const override { return true; }

    // Additional ComputationalMethod interface requirements
    bool hasGradient() const override { return true; }
    void setThreadCount(int threads) override { (void)threads; }  // PM6 uses global threads
    void setParameters(const json& params) override { (void)params; }  // Parameters set in constructor
    json getParameters() const override { return json{}; }
    bool hasError() const override { return false; }

    json getEnergyDecomposition() const;
    bool saveToFile(const std::string& filename) const override { return false; }

private:
    std::unique_ptr<PM6> m_pm6;
    Mol m_molecule;
    bool m_calculation_done;
    double m_last_energy;

    static json getDefaultConfig();
};
