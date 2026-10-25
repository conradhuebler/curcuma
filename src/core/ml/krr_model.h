/*
 * < Sparse Kernel-Ridge model for Delta-Machine-Learning >
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
 * Claude Generated (July 2026): regression layer of the Delta-ML work package.
 */

#pragma once

#include "descriptor.h"
#include "ml_dataset.h"

#include "src/core/global.h"
#include "json.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace curcuma {
namespace ml {

/**
 * @brief Knobs for fitting the delta model.
 *
 * n_centers is the single lever that trades accuracy against evaluation cost. The model is
 * a sum of M kernel bumps, so a prediction costs O(M) regardless of how much training data
 * was collected. Small M gives a cheap model, large M approaches the dense kernel-ridge
 * limit at dense cost. Which end pays off depends entirely on how expensive the reference
 * method is relative to the base method -- that is an empirical question, so M is a
 * parameter and not a hard-coded choice.
 */
struct DeltaMLTrainingOptions {
    int n_centers = 1000;       ///< kernel centres M; <= 0 uses every training environment
    double kernel_sigma = 0.0;  ///< kernel width; <= 0 selects the median-distance heuristic
    /**
     * @brief RKHS regulariser strength, relative to the data term.
     *
     * The system solved is (Phi^T W Phi + lambda * s * K_MM) beta = Phi^T W y, where s is
     * the mean diagonal of the data block. Referring lambda to the data scale is what makes
     * it a usable dial: as an absolute weight it is swamped by ten orders of magnitude, and
     * rescaling W uniformly cancels out of the solution entirely.
     *
     * lambda = 1 regularises as strongly as the data constrains; small values approach
     * interpolation. Too small produces a rippled surface that is invisible in the energy
     * error and fails to conserve energy in dynamics -- raise it if an NVE run drifts.
     */
    double lambda = 1e-6;
    double jitter_relative = 1e-8; ///< ridge added to the normal equations, relative to their
                                   ///< diagonal scale; guards the factorisation at large M
    double sigma_energy = 1e-4; ///< assumed energy noise [Eh], sets the energy row weight
    double sigma_force = 1e-3;  ///< assumed force noise [Eh/Angstrom], sets the force weight
    bool use_forces = true;     ///< include gradient observations (3N per structure)
    unsigned seed = 42;         ///< seed for the farthest-point-sampling tie-breaks
    int verbosity = 1;
    /**
     * @brief GPU backend for the normal-equation accumulation ("cuda", "rocm", ...).
     *
     * Empty or "none" fits on the CPU. The fit is dominated by one symmetric rank-k update
     * per structure (M^2 * 3N), which is what gets offloaded; everything else stays where
     * it is. Falls back silently to the CPU when the plugin or the device is unavailable.
     */
    std::string gpu;
    /*! \brief Descriptor grid; default-constructed = built-in configuration (Claude Generated, Oct 2026). */
    ACSFSettings acsf;
};

/*! \brief Error metrics of a model on a structure set. */
struct DeltaMLErrors {
    int n_structures = 0;

    // Absolute agreement with the reference total energy.
    double energy_mae = 0.0;  ///< [Eh]
    double energy_rmse = 0.0; ///< [Eh]

    // Agreement after removing the mean signed error. This is the number that matters for
    // conformer work: a constant offset shifts every conformer alike and cancels in any
    // relative energy.
    double energy_mae_centered = 0.0; ///< [Eh]

    // Claude Generated (Oct 2026): signed errors (model minus reference). The minimum is the
    // deepest "energy hole" of the set, which a mean error averages away.
    double energy_mean_signed = 0.0; ///< [Eh]
    double energy_min_signed = 0.0;  ///< [Eh]

    double force_rmse = 0.0; ///< [Eh/Angstrom], 0 if the set carries no gradients
    double spearman = 0.0;   ///< rank correlation of predicted vs reference energies

    // Same metrics for the uncorrected base method, so the gain is visible rather than
    // asserted. Without this reference point an MAE is meaningless.
    double base_energy_mae_centered = 0.0;
    double base_force_rmse = 0.0;
    double base_spearman = 0.0;
};

/*! \brief What the fit did, for logging and for the model file. */
struct DeltaMLTrainingReport {
    int n_structures = 0;
    int n_environments = 0;
    int n_centers = 0;
    int n_features = 0;
    double kernel_sigma = 0.0;
    double lambda_used = 0.0;
    double jitter_used = 0.0;
    int n_energy_rows = 0;
    int n_force_rows = 0;
    double train_seconds = 0.0;
    bool solver_ok = false;
    DeltaMLErrors train_errors;
};

/**
 * @brief Sparse kernel-ridge model of the energy difference between two methods.
 *
 * Model form (subset-of-regressors, the same construction GAP uses):
 *
 *   Delta(A) = sum_Z n_Z(A) c_Z  +  sum_{i in A} sum_{m=1..M} beta_m k(x_i, z_m)
 *
 * The first term is an element-count baseline that absorbs the large constant offset
 * between two methods with different energy zeros. It has zero gradient and therefore
 * never touches forces. The second term is the learned part; the kernel is element-blocked,
 * k(x_i, z_m) = 0 unless the two environments belong to the same element.
 *
 * Training uses energies and gradients jointly. One reference calculation yields a single
 * energy but 3N gradient components, and that ratio is what makes a few hundred reference
 * points sufficient.
 */
class DeltaMLModel {
public:
    DeltaMLModel() = default;

    /**
     * @brief Fit the model.
     * @param data descriptor-bearing structures with base and reference labels
     * @param opt  training options
     * @param report filled with fit statistics
     * @return false if the data are unusable or the linear solve fails
     */
    bool train(const MLDataset& data, const DeltaMLTrainingOptions& opt, DeltaMLTrainingReport& report);

    bool isTrained() const { return m_trained; }

    /*! \brief Predicted energy correction [Eh]; optionally the extrapolation indicator too. */
    double predictEnergy(const std::vector<int>& elements, const Geometry& geometry,
        double* uncertainty = nullptr) const;

    /**
     * @brief Predicted correction and its gradient.
     * @param gradient output, natoms x 3 [Eh/Angstrom]; resized by the call
     * @param uncertainty optional output, the extrapolation indicator
     * @return energy correction [Eh]
     *
     * The indicator comes for free: the per-atom kernel row against the element block is
     * computed for the prediction anyway, and its maximum is exactly what the indicator
     * needs. Calling uncertainty() separately would repeat the whole descriptor pass, which
     * in dynamics means paying for it on every single step.
     */
    double predictEnergyAndGradient(const std::vector<int>& elements, const Geometry& geometry,
        Matrix& gradient, double* uncertainty = nullptr) const;

    /**
     * @brief Cheap extrapolation indicator, max over atoms of 1 - max_m k(x_i, z_m).
     *
     * Zero means every atomic environment sits on top of a kernel centre, one means the
     * structure is entirely outside the region the model has seen. This measures distance
     * in descriptor space, not error: it answers "have I been here before", never "am I
     * right here". Useful to steer active learning and to flag suspicious predictions,
     * useless as an accuracy guarantee.
     */
    double uncertainty(const std::vector<int>& elements, const Geometry& geometry) const;

    /*! \brief Error metrics on a labelled set, including the uncorrected base for contrast. */
    DeltaMLErrors evaluate(const MLDataset& data) const;

    bool save(const std::string& filename) const;
    bool load(const std::string& filename);
    json toJson() const;
    bool fromJson(const json& j);

    const ACSFConfig& descriptorConfig() const { return m_acsf; }
    int centerCount() const { return static_cast<int>(m_centers.rows()); }
    int featureCount() const { return static_cast<int>(m_centers.cols()); }
    double kernelSigma() const { return m_sigma; }

    void setBaseMethod(const std::string& m) { m_base_method = m; }
    const std::string& baseMethod() const { return m_base_method; }
    void setReferenceMethod(const std::string& m) { m_ref_method = m; }
    const std::string& referenceMethod() const { return m_ref_method; }

    /**
     * @brief What the model was fitted for: "joint" (energies and gradients) or "energy"
     *        (energies only, use_forces = false).
     *
     * An energy-only fit has a smaller relative energy error but its forces are no better than those of
     * the base method, so consumers that need gradients warn on it. Stored as "role" in the
     * model file; a file without the field (all models before October 2026) reads as "joint".
     * Claude Generated (October 2026).
     */
    void setRole(const std::string& r) { m_role = r; }
    const std::string& role() const { return m_role; }

    /**
     * @brief Anchor state: the bond graph every consumer must reuse.
     *
     * Only the graph is pinned. Everything geometry dependent -- coordination numbers, the
     * force-field parameters derived from them, charges -- has to be recomputed at the
     * structure actually being evaluated. Pinning those as well produces a potential whose
     * gradient no longer matches its own energy, measured at 4.7e-3 Eh/Angstrom against
     * 3.5e-5 for the unpinned method.
     *
     * The geometry is kept for provenance and as a fallback when no graph was captured.
     */
    void setAnchor(const std::vector<int>& elements, const Geometry& geometry,
        const std::vector<std::pair<int, int>>& bonds = {});

    /**
     * @brief Move prediction onto a GPU backend ("cuda", ...). Empty returns to the CPU.
     *
     * Uploads centres, weights and scaling once. Falls back silently when the plugin or the
     * device is unavailable, so callers need no capability check.
     */
    void setInferenceBackend(const std::string& backend);
    bool hasAnchor() const { return m_anchor_geometry.rows() > 0; }
    const std::vector<int>& anchorElements() const { return m_anchor_elements; }
    const Geometry& anchorGeometry() const { return m_anchor_geometry; }
    const std::vector<std::pair<int, int>>& anchorBonds() const { return m_anchor_bonds; }

private:
    /*! \brief Apply the stored feature standardisation in place. */
    void standardise(Matrix& X) const;

    /*! \brief Element-count baseline contribution for one structure. */
    double baselineEnergy(const std::vector<int>& elements) const;

    /**
     * @brief Kernel row of one environment against its element block.
     * @param start output, first centre index of the block
     * @param count output, number of centres in the block
     * @return kernel values, length count (empty if the element is unknown to the model)
     */
    Vector kernelRow(const Vector& x, int element, int& start, int& count) const;

    ACSFConfig m_acsf;

    Matrix m_centers;                   ///< M x F, standardised, grouped by element
    std::vector<int> m_center_elements; ///< element of each centre, sorted ascending
    std::map<int, std::pair<int, int>> m_element_block; ///< element -> (start, count)

    Vector m_beta; ///< M regression weights

    Vector m_feature_mean;  ///< F
    Vector m_feature_scale; ///< F, never zero

    std::map<int, double> m_element_energy; ///< element-count baseline c_Z [Eh]

    double m_sigma = 1.0;
    bool m_trained = false;

    std::string m_base_method;
    std::string m_ref_method;
    std::string m_role = "joint"; ///< "joint" or "energy", see role()

    std::vector<int> m_anchor_elements;
    Geometry m_anchor_geometry;
    std::vector<std::pair<int, int>> m_anchor_bonds;

    /// Prediction is const, so the device context is mutable state of an otherwise pure call.
    mutable std::shared_ptr<class GpuInference> m_gpu_infer;
};

} // namespace ml
} // namespace curcuma
