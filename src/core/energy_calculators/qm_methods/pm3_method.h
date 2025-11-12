/*
 * <PM3 Method Wrapper>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 * This program is free software under GPL-3.0
 */

#pragma once

#include "../computational_method.h"
#include "pm3.h"
#include "src/core/molecule.h"
#include <memory>

class PM3Method : public ComputationalMethod {
public:
    explicit PM3Method(const json& config = json{});
    ~PM3Method() = default;

    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override { return Vector::Zero(0); }
    Vector getDipole() const override { return Vector::Zero(3); }

    std::string getMethodName() const override { return "pm3"; }
    bool isThreadSafe() const override { return true; }

    json getEnergyDecomposition() const;
    bool saveToFile(const std::string& filename) const override { return false; }

private:
    std::unique_ptr<PM3> m_pm3;
    Mol m_molecule;
    bool m_calculation_done;
    double m_last_energy;

    static json getDefaultConfig();
};
