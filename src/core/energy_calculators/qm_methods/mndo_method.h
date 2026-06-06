/*
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 *
 * Claude Generated: MNDO method wrapper for MethodFactory integration
 */

#pragma once

#include "../computational_method.h"
#include "mndo.h"
#include "src/core/molecule.h"
#include <memory>

class MNDOMethod : public ComputationalMethod {
public:
    explicit MNDOMethod(const json& config = json{});
    ~MNDOMethod() = default;

    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override { return Vector::Zero(0); }
    Position getDipole() const override { return Position::Zero(); }

    std::string getMethodName() const override { return "mndo"; }
    bool isThreadSafe() const override { return true; }

    // Additional ComputationalMethod interface requirements
    bool hasGradient() const override { return true; }
    void setThreadCount(int threads) override { (void)threads; }  // MNDO uses global threads
    void setParameters(const json& params) override { (void)params; }  // Parameters set in constructor
    json getParameters() const override { return json{}; }
    bool hasError() const override { return false; }

    json getEnergyDecomposition() const;
    bool saveToFile(const std::string& filename) const override { return false; }

private:
    std::unique_ptr<MNDO> m_mndo;
    Mol m_molecule;
    bool m_calculation_done;
    double m_last_energy;

    static json getDefaultConfig();
};
