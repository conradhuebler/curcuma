/*
 * <AncOpt Optimizer - Approximate Normal Coordinate Optimizer>
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This file contains a port of the AncOpt (Approximate Normal Coordinate Optimizer)
 * algorithm from XTB by Stefan Grimme and contributors.
 *
 * ORIGINAL IMPLEMENTATION:
 * File: external/xtb/src/optimizer.f90
 * Author: Stefan Grimme
 * Copyright (C) 2017-2020 Stefan Grimme
 * License: LGPL-3.0-or-later
 *
 * LITERATURE REFERENCE:
 * - Birkholz, A. B., & Schlegel, H. B. (2015)
 *   "Exploring the efficacy of the Bofill update in geometry optimization"
 *   Theoretical Chemistry Accounts, 135, 84
 * - AncOpt is unpublished work by Stefan Grimme
 *   Ported with respect to the original author
 *
 * This C++ port was created by Claude AI under instruction from Conrad Hübler.
 * All copyright remains with Conrad Hübler as project owner.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#pragma once

#include "optimizer_driver.h"
#include "src/core/molecule.h"
#include <Eigen/Dense>
#include <memory>
#include <vector>

using curcuma::Molecule;

namespace Optimization {

/**
 * @brief ANC (Approximate Normal Coordinate) data structure
 *
 * Contains the transformation between Cartesian and internal coordinates
 * based on approximate normal modes.
 */
struct ANCCoordinates {
    int n_atoms = 0; // Number of atoms
    int n3 = 0; // 3 * n_atoms (Cartesian dimension)
    int nvar = 0; // Number of internal coordinates (3N - 6 or 3N - 5)

    Matrix B; // Transformation matrix B (n3 x nvar)
    Vector coord; // Current internal coordinates
    Vector coord_reference; // Reference internal coordinates (Claude Nov 2025: Fix for ANC transform bug)
    Vector xyz_flat; // Current Cartesian coordinates (flattened)
    Vector xyz_reference; // Reference Cartesian coordinates (Claude Nov 2025: Fix for ANC transform bug)

    Matrix hess; // Hessian in internal coordinates (packed format)

    double hlow = 0.01; // Lower frequency cutoff
    double hmax = 5.0; // Upper frequency cutoff

    bool initialized = false;

    // Methods
    void allocate(int natoms, int num_vars, double h_low, double h_max);
    void deallocate();

    // Generate ANC from model Hessian
    bool generateANC(const Matrix& cartesian_hessian, const Vector& xyz, bool is_linear);

    // Transform between coordinate systems
    void getCartesian(Vector& xyz_out) const;
    void setCartesian(const Vector& xyz_in);
    Vector transformGradientToInternal(const Vector& cartesian_gradient) const;
};

/**
 * @brief Model Hessian parameters
 * Ported from XTB's model hessian settings
 */
struct ModelHessianParameters {
    enum Type {
        LINDH_1995, // Original Lindh (1995) model
        LINDH_2007, // Updated Lindh (2007) model
        LINDH_D2,   // Lindh with D2 dispersion corrections
        SWART,      // Swart model hessian
        READ_FILE   // Read from external file
    };

    Type model = LINDH_2007;
    double s6 = 30.0; // Scaling factor for model hessian

    static ModelHessianParameters fromJson(const json& config);
};

/**
 * @brief AncOpt Optimizer Implementation - Claude Generated
 *
 * Port of Stefan Grimme's AncOpt (Approximate Normal Coordinate Optimizer)
 * from XTB. Uses rational function optimization in approximate normal
 * coordinate space with BFGS Hessian updates.
 *
 * ALGORITHM OVERVIEW:
 * 1. Generate model Hessian (Lindh or similar)
 * 2. Transform to approximate normal coordinates (ANC)
 * 3. Optimize in ANC space using Rational Function (RF) method
 * 4. Update Hessian with BFGS or Powell
 * 5. Periodically regenerate ANC from updated Hessian
 *
 * SCIENTIFIC BACKGROUND:
 * - Uses eigenvectors of model Hessian as approximate normal modes
 * - Optimization in normal mode space reduces coupling between coordinates
 * - Rational Function method: Solves augmented eigenvalue problem
 *   | H  g |   | dx |       | dx |
 *   | g  0 | * | 1  | = λ * | 1  |
 * - BFGS updates maintain Hessian approximation between ANC regenerations
 */
class ANCOptimizer : public OptimizerDriver {
private:
    // ANC-specific parameters
    double m_maxdispl = 0.3; // Maximum displacement in ANC (Bohr)
    double m_hlow = 0.01; // Lower frequency cutoff
    double m_hmax = 5.0; // Upper frequency cutoff
    int m_maxmicro = 25; // Max micro-iterations before ANC regeneration
    int m_micro_current = 0; // Current micro-iteration

    ModelHessianParameters m_model_hess_params;

    // Hessian update method
    enum HessianUpdate {
        BFGS = 0,
        POWELL = 1
    };
    HessianUpdate m_hessian_update = BFGS;

    // ANC data structure
    std::unique_ptr<ANCCoordinates> m_anc;

    // Optimization state
    Vector m_displ; // Current displacement in internal coordinates
    Vector m_gint; // Gradient in internal coordinates
    Vector m_gint_old; // Previous gradient
    double m_gnorm = 0.0; // Gradient norm
    double m_gnorm_old = 0.0; // Previous gradient norm
    double m_depred = 0.0; // Predicted energy change

    // Convergence parameters (from XTB)
    double m_ethr = 5e-6; // Energy threshold (Eh)
    double m_gthr = 1e-3; // Gradient threshold (Eh/Bohr)
    double m_acc = 1.0; // SCC accuracy scaling

    // Micro-iteration control
    bool m_needs_anc_regeneration = true;

protected:
    // OptimizerDriver interface implementation
    bool InitializeOptimizerInternal() override;

    Vector CalculateOptimizationStep(const Vector& current_coordinates,
                                      const Vector& gradient) override;

    bool CheckMethodSpecificConvergence() const override;

    void UpdateOptimizerState(const Vector& new_coordinates,
                              const Vector& new_gradient,
                              double new_energy) override;

    void FinalizeOptimizationInternal() override;

    // AncOpt-specific methods

    /**
     * @brief Generate model Hessian in Cartesian coordinates
     * Ported from XTB's modhes() routine
     */
    Matrix generateModelHessian(const Molecule& mol);

    /**
     * @brief Generate ANC from Cartesian Hessian
     * Ported from XTB's tb_anc type
     */
    bool generateANCFromHessian(const Matrix& cart_hess, const Molecule& mol);

    /**
     * @brief Rational Function step calculation
     * Solves augmented eigenvalue problem for displacement
     * Ported from XTB optimizer.f90:676-729
     */
    Vector calculateRationalFunctionStep(const Vector& gradient_internal,
                                         const Matrix& hessian_internal);

    /**
     * @brief BFGS Hessian update
     * Ported from XTB bfgs.f90
     */
    void updateHessianBFGS(Matrix& hessian, const Vector& dx, const Vector& dg);

    /**
     * @brief Powell Hessian update
     * Ported from XTB broyden.f90
     */
    void updateHessianPowell(Matrix& hessian, const Vector& dx, const Vector& dg);

    /**
     * @brief Predict energy change from 2nd order model
     * ΔE = g·dx + 0.5·dx·H·dx
     */
    double predictEnergyChange(const Vector& gradient, const Vector& displacement,
                               const Matrix& hessian);

    /**
     * @brief Project out translations and rotations from Hessian
     * Ported from XTB trproj()
     */
    void projectTranslationsRotations(Matrix& hessian, const Molecule& mol);

    /**
     * @brief Get optimization thresholds based on level
     * Ported from XTB get_optthr()
     */
    void setOptimizationLevel(int level);

    // Helper methods
    void loadANCParameters(const json& config);
    bool checkLinearMolecule(const Molecule& mol) const;

public:
    ANCOptimizer();
    virtual ~ANCOptimizer() = default;

    // OptimizerInterface implementation
    std::string getName() const override {
        return "AncOpt (Approximate Normal Coordinate Optimizer)";
    }

    OptimizerType getType() const override {
        return OptimizerType::ANCOPT;
    }

    bool supportsConstraints() const override { return true; }
    bool requiresHessian() const override { return false; } // Generates model Hessian

    std::vector<std::string> getRequiredParameters() const override {
        return {"maxdispl", "hlow", "hmax", "maxmicro", "model_hessian"};
    }

    json GetDefaultConfiguration() const override;

    // Configuration setters
    void setMaxDisplacement(double maxd) { m_maxdispl = maxd; }
    void setFrequencyCutoffs(double hlow, double hmax) {
        m_hlow = hlow;
        m_hmax = hmax;
    }
    void setMicroIterations(int micro) { m_maxmicro = micro; }
    void setHessianUpdate(HessianUpdate method) { m_hessian_update = method; }
    void setModelHessianType(ModelHessianParameters::Type type) {
        m_model_hess_params.model = type;
    }
};

} // namespace Optimization
