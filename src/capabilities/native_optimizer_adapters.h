/*
 * <Native Optimizer Adapters - Wrap lbfgs.cpp algorithms as OptimizerDriver subclasses>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Apr 2026): Adapter classes for native L-BFGS, DIIS, RFO
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#pragma once

#include "optimizer_driver.h"
#include "src/capabilities/optimisation/lbfgs.h"

namespace Optimization {

/**
 * @brief Base adapter for native optimization algorithms from lbfgs.cpp
 * Claude Generated
 *
 * Wraps the LBFGS class (which provides LBFGSStep, DIISStep, RFOStep)
 * as an OptimizerDriver subclass. The scientific algorithm code in
 * lbfgs.cpp remains untouched — this adapter only bridges the interface.
 */
class NativeOptimizerAdapter : public OptimizerDriver {
protected:
    std::unique_ptr<LBFGS> m_lbfgs;
    LBFGS::Method m_native_method;

    // OptimizerDriver template method hooks
    bool InitializeOptimizerInternal() override;
    Vector CalculateOptimizationStep(const Vector& current_coordinates,
        const Vector& gradient) override;
    bool CheckMethodSpecificConvergence() const override;
    void UpdateOptimizerState(const Vector& new_coordinates,
        const Vector& new_gradient,
        double new_energy) override;
    void FinalizeOptimizationInternal() override;

    NativeOptimizerAdapter(LBFGS::Method method);

public:
    virtual ~NativeOptimizerAdapter() = default;
    bool supportsConstraints() const override { return true; }
    bool requiresHessian() const override { return false; }
};

/**
 * @brief Native L-BFGS optimizer — Nocedal & Wright two-loop recursion
 * Claude Generated
 */
class NativeLBFGSAdapter : public NativeOptimizerAdapter {
public:
    NativeLBFGSAdapter();
    std::string getName() const override { return "Native L-BFGS"; }
    OptimizerType getType() const override { return OptimizerType::NATIVE_LBFGS; }
    std::vector<std::string> getRequiredParameters() const override { return {}; }
};

/**
 * @brief Native DIIS optimizer — Pulay (1980) direct inversion of iterative subspace
 * Claude Generated
 */
class NativeDIISAdapter : public NativeOptimizerAdapter {
public:
    NativeDIISAdapter();
    std::string getName() const override { return "Native DIIS"; }
    OptimizerType getType() const override { return OptimizerType::NATIVE_DIIS; }
    std::vector<std::string> getRequiredParameters() const override { return {}; }
};

/**
 * @brief Native RFO optimizer — Banerjee et al. (1985) rational function optimization
 * Claude Generated
 */
class NativeRFOAdapter : public NativeOptimizerAdapter {
public:
    NativeRFOAdapter();
    std::string getName() const override { return "Native RFO"; }
    OptimizerType getType() const override { return OptimizerType::NATIVE_RFO; }
    std::vector<std::string> getRequiredParameters() const override { return {}; }
};

} // namespace Optimization
