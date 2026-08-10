/*
 * <QMDFF Hessian fit for parametrisation. >
 * Copyright (C) 2023 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 */

#pragma once

#include "src/core/config_manager.h"
#include "src/core/energy_calculators/ff_methods/qmdff_par.h"
#include "src/core/energy_calculators/ff_methods/qmdff_parametrisation.h"
#include "src/core/molecule.h"

#include <vector>
#include "src/core/parameter_macros.h"

#include "curcumamethod.h"

/* Claude Generated 2025: QMDFFfit Parameter Registry - replaces static QMDFFFitJson */
BEGIN_PARAMETER_DEFINITION(qmdfffit)
    PARAM(method, String, "gfn2", "QM method supplying the reference Hessian and charges.", "General", {})
    PARAM(hessian_file, String, "hessian.json", "Input Hessian file.", "Input", {"hessian"})
    PARAM(charges_file, String, "scf.json", "Input charges file.", "Input", {"charges"})
    PARAM(potential, String, "qmdff", "Force field to parametrise - qmdff or gfnff. gfnff keeps GFN-FFs non-covalent block and fits only its bonded constants.", "General", {})
    PARAM(threads, Int, 1, "Number of threads.", "Performance", {})
    // Claude Generated (Aug 2026): controls of the linear Hessian fit, see docs/QMDFF.md
    PARAM(basins, Int, 1, "Number of fitting basins taken from a multi-frame input; 1 reproduces the single-structure fit.", "Fit", {})
    PARAM(basin_frames, String, "", "Explicit comma-separated frame indices to fit at; overrides basins.", "Fit", {})
    PARAM(weight_energy, Double, 1.0, "Weight of the relative-energy residual block, which is what constrains conformer ranking.", "Fit", {})
    PARAM(weight_gradient, Double, 0.0, "Weight of the gradient residual block; 0 by default because r0/theta0 belong to one basin and non-zero weights switch the force field off.", "Fit", {})
    PARAM(hessian_weight, String, "cartesian", "Inner product used by the fit - cartesian or mass.", "Fit", {})
    PARAM(lambda, Double, -1.0, "Tikhonov regularisation toward the Seminario guess; negative selects it automatically.", "Fit", {})
    PARAM(fd_step, Double, 1e-4, "Central-difference step in Bohr for the per-term unit Hessians.", "Fit", {})
    PARAM(lambda_torsion, Double, 1.0, "Separate Tikhonov prior keeping torsion and out-of-plane barriers near their rule-derived values; a Hessian does not determine barrier heights.", "Fit", {})
    PARAM(fit_torsions, Bool, true, "Fit the QMDFF torsion scale factors.", "Fit", {})
    PARAM(fit_oop, Bool, true, "Fit the QMDFF out-of-plane scale factors.", "Fit", {})
    PARAM(tie_torsions, Bool, true, "Tie all torsions about the same central bond to one shared scale.", "Fit", {})
    PARAM(nonnegative, Bool, true, "Constrain fitted force constants to their physical lower bounds.", "Fit", {})
    PARAM(verify_fd_step, Double, 2e-4, "Finite-difference step in Angstrom for the verification force-field Hessian.", "Fit", {})
    PARAM(verify, Bool, true, "Recompute one force-field Hessian to verify the fit and the parameter injection.", "Fit", {})
    PARAM(write_param, Bool, true, "Also write <input>.param.json so later runs load the fit automatically.", "Output", {})
END_PARAMETER_DEFINITION

class QMDFFFit : public CurcumaMethod {
public:
    /**
     * @brief Constructor with JSON configuration (backward compatible)
     * Claude Generated: Phase 4 - ConfigManager Migration
     */
    QMDFFFit(const json& controller, bool silent = true);

    /**
     * @brief Constructor with ConfigManager configuration (new, preferred)
     * Claude Generated: Phase 4 - Native ConfigManager support
     */
    QMDFFFit(const ConfigManager& config, bool silent = true);
    virtual ~QMDFFFit() {}
    void start() override;

    void setMolecule(const Molecule& molecule) { m_molecule = molecule; }

    /// Claude Generated (Aug 2026): all frames of a multi-structure input. When more than
    /// one basin is requested the fit is constrained by the Hessian, energy AND gradient
    /// at several of them at once — a single Hessian only fixes curvature, which cannot
    /// determine relative basin energies (docs/QMDFF_CONFORMER_USECASE.md).
    void setEnsemble(const std::vector<Molecule>& ensemble) { m_ensemble = ensemble; }

    /// Claude Generated (Aug 2026): input path, used for the <input>.param.json output.
    void setInputFile(const std::string& file) { m_input_file = file; }

private:
    /* Lets have this for all modules */
    nlohmann::json WriteRestartInformation() override { return json(); }

    /* Lets have this for all modules */
    bool LoadRestartInformation() override { return true; }

    StringList MethodName() const override { return { std::string("QMDFFFit") }; }

    /* Lets have all methods read the input/control file */
    void ReadControlFile() override {}

    /* Read Controller has to be implemented for all */
    void LoadControlJson() override;

    // --- Claude Generated (Aug 2026): steps of the rewritten fit --------------
    /// Obtain e0, QM partial charges and |g|max — from charges_file if it exists.
    bool prepareCharges(double& e0, Vector& charges, double& gmax);
    /// Obtain the RAW Cartesian reference Hessian (Hartree/Angstrom^2).
    bool prepareHessian(Matrix& hessian, Vector& qm_frequencies);
    /// Build the fit options from the registry parameters.
    curcuma::qmdff::FitOptions fitOptions() const;
    /// Which frames of the ensemble to use as fitting basins.
    std::vector<int> selectBasinFrames();
    /// One force-field Hessian on the fitted set; proves the injection and the fit.
    void verifyFit(const json& fitted, const curcuma::qmdff::FitReport& report,
        const Vector& qm_frequencies);
    /// Fit GFN-FF's bonded constants, keeping its non-covalent block untouched.
    void runGFNFFFit(double e0, const Vector& qm_frequencies);

    Molecule m_molecule;
    std::vector<Molecule> m_ensemble;
    std::string m_method = "gfn2", m_hessian_file, m_scf_file, m_potential = "qmdff";
    std::string m_input_file;
    Matrix m_hessian;

    // Fit controls
    std::string m_hessian_weight = "cartesian";
    double m_lambda = -1.0;
    double m_fd_step = 1e-4;
    double m_lambda_torsion = 1.0;
    bool m_fit_torsions = true, m_fit_oop = true, m_tie_torsions = true;
    bool m_nonnegative = true, m_verify = true, m_write_param = true;
    double m_verify_fd_step = 2e-4;
    int m_basins = 1;
    std::string m_basin_frames;
    double m_weight_energy = 1.0, m_weight_gradient = 1.0;

    int m_threads = 1;
};
