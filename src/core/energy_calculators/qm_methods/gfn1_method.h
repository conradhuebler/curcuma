/*
 * <GFN1 Method Wrapper>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#pragma once

#include "../computational_method.h"
#include "gfn1.h"
#include "src/core/molecule.h"

#include <memory>

/**
 * @brief ComputationalMethod wrapper for native GFN1-xTB
 *
 * Wraps the native GFN1 implementation to provide the unified
 * ComputationalMethod interface for integration into Curcuma's
 * method factory system.
 *
 * Claude Generated: GFN1 method wrapper for Curcuma
 */
class GFN1Method : public ComputationalMethod {
public:
    explicit GFN1Method(const json& config = json{});
    ~GFN1Method() = default;

    // ComputationalMethod interface
    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override { return Vector::Zero(0); }
    Vector getDipole() const override { return Vector::Zero(3); }

    std::string getMethodName() const override { return "gfn1"; }
    bool isThreadSafe() const override { return true; }

    // GFN1-specific
    Vector getCoordinationNumbers() const;
    json getEnergyDecomposition() const;
    bool saveToFile(const std::string& filename) const override;

private:
    std::unique_ptr<GFN1> m_gfn1;
    Mol m_molecule;
    bool m_calculation_done;
    double m_last_energy;

    static json getDefaultConfig();
};
