/*
 * < Delta-ML corrected computational method >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Claude Generated (July 2026): method wrapper of the Delta-ML work package.
 */

#pragma once

#include "src/core/energy_calculators/computational_method.h"
#include "src/core/ml/krr_model.h"
#include "src/core/parameter_macros.h"

#include <memory>
#include <string>

BEGIN_PARAMETER_DEFINITION(deltaml)
PARAM(model, String, "",
    "Path to the delta-ML model JSON produced by '-activelearn'. Required. The shorthand "
    "'-method deltaml:<path>' sets this too and survives capabilities that rebuild their "
    "sub-configuration from the method name alone.",
    "Basic", {})
PARAM(base, String, "",
    "Base method the correction is added to. Empty uses the method stored in the model "
    "file, which is what you normally want: mixing a model with a different base than it "
    "was trained against produces silently wrong energies.",
    "Basic", {})
PARAM(enabled, Bool, true,
    "Evaluate the ML correction. Set false to run the bare base method through the same "
    "code path, for A/B comparisons.",
    "Basic", {})
PARAM(uncertainty_warn, Double, 0.9,
    "Warn once per run when the extrapolation indicator 1 - max_m k(x,z_m) exceeds this. "
    "Mind the scale before lowering it: the indicator is 1 - exp(-d^2/2 sigma^2) for a "
    "nearest centre at distance d, so 0.5 is reached at a perfectly ordinary d = 1.18 sigma "
    "and would fire on almost every structure. 0.9 corresponds to d = 2.1 sigma, where the "
    "correction is attenuated to about a tenth and the prediction really is extrapolation.",
    "Advanced", {})
PARAM(report_uncertainty, Bool, false,
    "Report the extrapolation indicator with every energy evaluation. The indicator itself "
    "is free (it falls out of the prediction), this only controls the logging.",
    "Advanced", {})
PARAM(uncertainty_record, String, "",
    "Append every structure whose indicator exceeds uncertainty_warn to this multi-XYZ file, "
    "with the indicator in the comment line. Empty disables it. This is what turns the "
    "warning into something actionable: the file is a ready-made input for another "
    "-activelearn run that extends the model exactly where it was out of its depth. "
    "'auto' derives the name from the model file.",
    "Advanced", {})
PARAM(uncertainty_abort, Double, 0.0,
    "Refuse to return an energy once the indicator exceeds this (0 disables). Intended for "
    "production optimisation, where a silently unsupported result is worse than no result. "
    "In dynamics it also serves as the trigger for stop-and-relabel: the run ends at the "
    "structure that left the trained region, which is precisely the one worth labelling.",
    "Advanced", {})
END_PARAMETER_DEFINITION

/**
 * @brief Base method plus a learned correction: E = E_base + Delta_ML.
 *
 * The whole point of implementing this as a ComputationalMethod is that every capability
 * in curcuma that takes '-method' -- single point, optimisation, dynamics, conformer
 * search -- then works on the corrected surface without any further change.
 *
 * Two details that are easy to get wrong and are handled here centrally:
 *
 * - **Topology anchor.** GFN-FF derives its bond list, hybridisation and pi system once,
 *   from whatever geometry it was initialised with. Different instances started from
 *   different conformers can therefore end up with different topologies, which injects
 *   noise into Delta that has nothing to do with the reference method. If the model
 *   carries an anchor geometry, the base method is always initialised from it and only
 *   then moved to the actual structure.
 *
 * - **Base method consistency.** A model trained against GFN-FF is meaningless on top of
 *   anything else. The base method is taken from the model file unless explicitly
 *   overridden, and a mismatch is reported.
 *
 * Claude Generated (July 2026)
 */
class DeltaMLMethod : public ComputationalMethod {
public:
    explicit DeltaMLMethod(const json& config);
    ~DeltaMLMethod() override = default;

    /*! \brief Extract the model path from a 'deltaml:<path>' method string, else empty. */
    static std::string modelPathFromMethodString(const std::string& method);

    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override { return m_gradient; }
    Vector getCharges() const override;
    Vector getBondOrders() const override;
    Position getDipole() const override;
    Vector getCN() const override;
    bool hasGradient() const override;

    std::string getMethodName() const override;
    bool isThreadSafe() const override;
    void setThreadCount(int threads) override;
    void setParameters(const json& params) override;
    json getParameters() const override { return m_parameters; }

    bool hasError() const override;
    std::string getErrorMessage() const override;
    void clearError() override;

    json getEnergyDecomposition() const override;

    void setWarmStart(bool on) override;
    void setIterativeMode(bool on) override;

    /*! \brief Learned correction of the last evaluation [Eh]. */
    double lastCorrection() const { return m_last_correction; }

    /*! \brief Extrapolation indicator of the last evaluation, or -1 if not computed. */
    double lastUncertainty() const { return m_last_uncertainty; }

    /*! \brief How many structures were written to the extrapolation record. */
    long long recordedStructures() const { return m_recorded; }

private:
    bool initialiseBase(const Mol& mol);

    /*! \brief Append the current geometry to the extrapolation record. */
    void recordExtrapolation(double indicator) const;

    /*! \brief Warn once per process when gradients are taken from an energy-only model.
     *  Claude Generated (October 2026). */
    void warnEnergyRoleOnce() const;
    int m_user_verbosity = 1; ///< logger level at construction, i.e. the caller's own

    std::unique_ptr<ComputationalMethod> m_base;
    curcuma::ml::DeltaMLModel m_model;

    json m_config; ///< full controller, forwarded to the base method

    std::string m_model_path;
    std::string m_base_method;
    bool m_enabled = true;
    bool m_report_uncertainty = false;
    double m_uncertainty_warn = 0.9;
    double m_uncertainty_abort = 0.0;
    std::string m_uncertainty_record;
    bool m_warned = false;
    mutable long long m_recorded = 0;

    std::vector<int> m_elements;
    Geometry m_geometry;

    Matrix m_gradient;
    double m_last_energy = 0.0;
    double m_last_base_energy = 0.0;
    double m_last_correction = 0.0;
    double m_last_uncertainty = -1.0;
};
