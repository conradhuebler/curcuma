/*
 * <Multi-CV Metadynamics Bias Engine>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (November 2025)
 *
 * MULTI-CV BIAS ENGINE
 * ====================
 *
 * MOTIVATION:
 * ----------
 * Current SimpleMD BiasThread only supports RMSD as single CV.
 * Modern metadynamics requires multiple collective variables for complex processes.
 *
 * DESIGN GOALS:
 * ------------
 * 1. Support arbitrary number of CVs (1 to ~3, limited by curse of dimensionality)
 * 2. Polymorphic CV interface (Distance, Angle, Dihedral, Gyration, Coordination)
 * 3. Efficient multi-dimensional Gaussian bias calculation
 * 4. Thread-safe for parallel bias evaluation
 * 5. Well-tempered and standard metadynamics support
 *
 * THEORETICAL BACKGROUND:
 * ----------------------
 * Multi-dimensional bias potential (see docs/MULTI_CV_METADYNAMICS.md):
 *
 *   V(s, t) = Σ_i w_i * exp(-Σ_α (s_α - s_α,i)² / (2σ_α²))
 *           = Σ_i w_i * ∏_α exp(-(s_α - s_α,i)² / (2σ_α²))
 *
 * where:
 *   s = (s₁, s₂, ..., s_d): Current CV values (d-dimensional)
 *   s_i: CV values when i-th Gaussian was deposited
 *   w_i: Gaussian height
 *   σ_α: Width along CV α
 *
 * BIAS FORCE CALCULATION (Multidimensional Chain Rule):
 *
 *   F_j^bias = -∇_r_j V(s, t)
 *           = -Σ_α (∂V/∂s_α) * (∂s_α/∂r_j)
 *
 * where:
 *   ∂V/∂s_α = Σ_i w_i * exp(...) * [-(s_α - s_α,i) / σ_α²]
 *   ∂s_α/∂r_j = CV gradient (from CV::CollectiveVariable::gradient())
 *
 * REFERENCES:
 * ----------
 * [1] Laio & Parrinello (2002). Escaping free-energy minima.
 *     Proc. Natl. Acad. Sci. USA 99, 12562-12566.
 *
 * [2] Barducci et al. (2008). Well-tempered metadynamics.
 *     Phys. Rev. Lett. 100, 020603.
 *
 * [3] Bussi et al. (2006). Equilibrium free energies from nonequilibrium metadynamics.
 *     Phys. Rev. Lett. 96, 090601.
 *
 * DESIGN PATTERN:
 * --------------
 * We use the **Strategy Pattern** to decouple bias calculation from MD engine.
 * BiasEngine can be used with SimpleMD, or any other MD code.
 */

#pragma once

#include "collective_variable.h"
#include "cv_factory.h"
#include "grid_accelerator.h"
#include "src/core/global.h"
#include "src/core/molecule.h"

#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include <chrono>

namespace CV {

/**
 * @struct GaussianHill
 * @brief Single Gaussian hill in multi-CV space
 *
 * STORAGE FORMAT:
 * Each Gaussian stores:
 * - Center coordinates in CV space: s = (s₁, s₂, ..., s_d)
 * - Widths along each CV: σ = (σ₁, σ₂, ..., σ_d)
 * - Height: w
 * - Time of deposition: t
 * - Well-tempered scaling factor
 *
 * MEMORY:
 * For d CVs: ~8*(2d + 4) bytes per Gaussian
 * Example: 2 CVs → ~64 bytes/Gaussian → 1000 Gaussians = 64 KB (negligible)
 */
struct GaussianHill {
    std::vector<double> cv_centers;  ///< s_i = (s₁,i, s₂,i, ..., s_d,i)
    std::vector<double> cv_sigmas;   ///< σ = (σ₁, σ₂, ..., σ_d)
    double height;                    ///< w_i (Gaussian amplitude)
    double time;                      ///< Deposition time (ps)
    double accumulated_bias;          ///< Σ V(s_i, t) for well-tempered MTD
    double well_tempered_factor;      ///< Scaling factor for WT-MTD
    int counter;                      ///< Number of times this region was visited

    /**
     * @brief Constructor
     *
     * @param centers CV values at deposition
     * @param sigmas Gaussian widths
     * @param h Initial height
     * @param t Deposition time
     */
    GaussianHill(const std::vector<double>& centers,
                 const std::vector<double>& sigmas,
                 double h,
                 double t = 0.0)
        : cv_centers(centers)
        , cv_sigmas(sigmas)
        , height(h)
        , time(t)
        , accumulated_bias(0.0)
        , well_tempered_factor(1.0)
        , counter(1)
    {
        if (centers.size() != sigmas.size()) {
            throw std::invalid_argument("GaussianHill: centers and sigmas must have same size");
        }
    }

    /**
     * @brief Get dimensionality (number of CVs)
     */
    int dimensionality() const { return cv_centers.size(); }
};

/**
 * @class BiasEngine
 * @brief Multi-CV metadynamics bias calculation engine
 *
 * THREAD-SAFETY:
 * -------------
 * This class is NOT thread-safe by default. For parallel MD:
 * - Create one BiasEngine per thread
 * - Synchronize Gaussian hills across threads periodically
 *
 * USAGE PATTERN:
 * -------------
 * ```cpp
 * // Setup
 * BiasEngine engine;
 * engine.addCV(CVFactory::create("distance", {0, 10}), 0.2);  // σ₁ = 0.2 Å
 * engine.addCV(CVFactory::create("angle", {5, 0, 10}), 15.0); // σ₂ = 15°
 * engine.setWellTempered(true, 3000.0);  // ΔT = 3000 K
 *
 * // During MD loop
 * auto [bias_energy, bias_gradient] = engine.computeBias(molecule);
 * // Apply bias_gradient to MD forces
 *
 * // Deposition
 * if (step % pace == 0) {
 *     engine.depositGaussian(molecule, height, current_time);
 * }
 * ```
 *
 * PERFORMANCE:
 * -----------
 * - Complexity: O(N_gaussians * N_cvs * N_atoms)
 * - Optimization: Early termination if Gaussian is >6σ away
 * - Future: Grid-based acceleration for >500 Gaussians
 *
 * @author Claude (Anthropic) & Conrad Hübler
 * @date November 2025
 */
class BiasEngine {
public:
    /**
     * @brief Constructor
     *
     * @param use_well_tempered Enable well-tempered metadynamics (recommended)
     * @param delta_T Bias temperature for WT-MTD (K)
     */
    BiasEngine(bool use_well_tempered = true, double delta_T = 3000.0)
        : m_use_well_tempered(use_well_tempered)
        , m_delta_T(delta_T)
        , m_initial_height(1.0)  // kJ/mol, will be set by user
        , m_total_bias(0.0)
        , m_num_gaussians(0)
        , m_grid(nullptr)
        , m_use_grid(false)
        , m_grid_auto(true)
        , m_grid_threshold(500)
        , m_grid_bins(100)
        , m_grid_initialized(false)
    {
    }

    /**
     * @brief Add a collective variable to the bias
     *
     * @param cv Unique pointer to CV object (moved into engine)
     * @param sigma Gaussian width for this CV
     *
     * NOTE: CVs must be added BEFORE depositing Gaussians
     *
     * EXAMPLE:
     * ```cpp
     * engine.addCV(CVFactory::create("distance", {0, 10}), 0.2);
     * ```
     */
    void addCV(CVPtr cv, double sigma) {
        if (sigma <= 0.0) {
            throw std::invalid_argument("BiasEngine: sigma must be > 0");
        }
        m_cvs.push_back(std::move(cv));
        m_cv_sigmas.push_back(sigma);
    }

    /**
     * @brief Set well-tempered metadynamics parameters
     *
     * @param use_wt Enable well-tempered (true) or standard MTD (false)
     * @param delta_T Bias temperature (K), typical: 5-10 × simulation temperature
     *
     * THEORY:
     * In WT-MTD, Gaussian height decreases as:
     *   w_i = w_0 * exp(-V(s_i, t_i) / (k_B * ΔT))
     *
     * This ensures convergence to:
     *   lim_{t→∞} V(s, t) = -(T + ΔT)/ΔT * F(s) + C
     */
    void setWellTempered(bool use_wt, double delta_T) {
        m_use_well_tempered = use_wt;
        m_delta_T = delta_T;
    }

    /**
     * @brief Set initial Gaussian height
     *
     * @param height Initial w_0 (kJ/mol)
     *
     * GUIDELINE:
     *   w_0 ≈ k_B * T / (20 to 100)
     *   For T = 300 K: w_0 ≈ 0.12 to 0.025 kJ/mol
     */
    void setInitialHeight(double height) {
        m_initial_height = height;
    }

    /**
     * @brief Compute bias potential and gradient
     *
     * @param mol Current molecular geometry
     * @return std::pair<double, Geometry> = {V_bias, ∇V_bias}
     *
     * ALGORITHM (Modified November 2025 - Claude Generated):
     * 1. Evaluate all CVs: s_α = CV_α.calculate(mol)
     * 2. Check if grid acceleration is active:
     *    - If yes: Use O(1) grid interpolation
     *    - If no: Use O(N_gaussians) direct summation
     * 3. Auto-enable grid if threshold reached
     *
     * COMPLEXITY:
     * - Direct: O(N_gaussians * N_cvs * N_atoms)
     * - Grid: O(N_cvs * N_atoms) = O(1) in N_gaussians
     */
    std::pair<double, Geometry> computeBias(const Molecule& mol) {
        if (m_cvs.empty()) {
            throw std::runtime_error("BiasEngine: No CVs added");
        }

        // Step 1: Evaluate all CVs
        std::vector<double> current_cv_values(m_cvs.size());
        for (size_t alpha = 0; alpha < m_cvs.size(); ++alpha) {
            current_cv_values[alpha] = m_cvs[alpha]->calculate(mol);
        }

        // Step 2: Check if we should enable grid (auto mode)
        if (m_grid_auto && !m_use_grid && m_num_gaussians >= m_grid_threshold) {
            m_use_grid = true;
        }

        // Step 3: Initialize grid if enabled but not yet initialized
        if (m_use_grid && !m_grid_initialized) {
            initializeGrid();
        }

        // Step 4: Compute bias (grid or direct)
        std::pair<double, Geometry> result;
        if (m_use_grid && m_grid_initialized) {
            result = computeBiasFromGrid(mol, current_cv_values);
        } else {
            result = computeBiasDirect(mol, current_cv_values);
        }

        m_total_bias = result.first;
        return result;
    }

    /**
     * @brief Deposit a new Gaussian at current position
     *
     * @param mol Current molecular geometry
     * @param height Gaussian height (if 0, use computed WT-MTD height)
     * @param time Current simulation time (ps)
     *
     * WELL-TEMPERED SCALING:
     * If use_well_tempered = true and height = 0:
     *   w_i = w_0 * exp(-V(s_i, t_i) / (k_B * ΔT))
     *
     * Modified November 2025 - Claude Generated:
     * Now also updates grid when grid acceleration is active.
     */
    void depositGaussian(const Molecule& mol, double height = 0.0, double time = 0.0) {
        if (m_cvs.empty()) {
            throw std::runtime_error("BiasEngine: No CVs added");
        }

        // Evaluate all CVs
        std::vector<double> cv_values(m_cvs.size());
        for (size_t alpha = 0; alpha < m_cvs.size(); ++alpha) {
            cv_values[alpha] = m_cvs[alpha]->calculate(mol);
        }

        // Compute height (well-tempered or fixed)
        double gaussian_height = height;
        if (m_use_well_tempered && height == 0.0) {
            // w_i = w_0 * exp(-V(s_i, t_i) / (k_B * ΔT))
            constexpr double kb_kJ = 0.0083144621;  // kJ/(mol·K)
            double scaling = std::exp(-m_total_bias / (kb_kJ * m_delta_T));
            gaussian_height = m_initial_height * scaling;
        } else if (height == 0.0) {
            gaussian_height = m_initial_height;
        }

        // Create and add Gaussian
        GaussianHill hill(cv_values, m_cv_sigmas, gaussian_height, time);
        m_gaussians.push_back(hill);
        m_num_gaussians++;

        // Update grid if active
        if (m_use_grid && m_grid_initialized) {
            m_grid->addGaussianToGrid(cv_values,
                                     gaussian_height * hill.well_tempered_factor,
                                     m_cv_sigmas);
        }
    }

    /**
     * @brief Get number of CVs
     */
    int getNumCVs() const { return m_cvs.size(); }

    /**
     * @brief Get number of deposited Gaussians
     */
    int getNumGaussians() const { return m_num_gaussians; }

    /**
     * @brief Get current total bias energy
     */
    double getTotalBias() const { return m_total_bias; }

    /**
     * @brief Get current CV values
     */
    std::vector<double> getCurrentCVValues(const Molecule& mol) const {
        std::vector<double> values(m_cvs.size());
        for (size_t i = 0; i < m_cvs.size(); ++i) {
            values[i] = m_cvs[i]->calculate(mol);
        }
        return values;
    }

    /**
     * @brief Get CV names/descriptions
     */
    std::vector<std::string> getCVDescriptions() const {
        std::vector<std::string> descs(m_cvs.size());
        for (size_t i = 0; i < m_cvs.size(); ++i) {
            descs[i] = m_cvs[i]->description();
        }
        return descs;
    }

    /**
     * @brief Write COLVAR file (time series of CV values and bias)
     *
     * @param filename Output filename
     * @param time Current time
     * @param mol Current molecule
     * @param append Append to existing file (true) or overwrite (false)
     *
     * FORMAT:
     * # time  cv1  cv2  ...  V_bias  N_gaussians
     */
    void writeCOLVAR(const std::string& filename, double time, const Molecule& mol, bool append = true) {
        std::ofstream file;
        if (append) {
            file.open(filename, std::ios::app);
        } else {
            file.open(filename);
            // Write header
            file << "# time";
            for (size_t i = 0; i < m_cvs.size(); ++i) {
                file << "  cv" << (i + 1);
            }
            file << "  V_bias  N_gaussians\n";
        }

        // Write data
        file << time;
        auto cv_values = getCurrentCVValues(mol);
        for (double val : cv_values) {
            file << "  " << val;
        }
        file << "  " << m_total_bias << "  " << m_num_gaussians << "\n";

        file.close();
    }

    /**
     * @brief Write HILLS file (Gaussian parameters)
     *
     * @param filename Output filename
     * @param append Append to existing file
     *
     * FORMAT (PLUMED-compatible):
     * # time  cv1_center  cv2_center  ...  sigma1  sigma2  ...  height  bias_factor
     */
    void writeHILLS(const std::string& filename, bool append = true) {
        std::ofstream file;
        if (append) {
            file.open(filename, std::ios::app);
        } else {
            file.open(filename);
            // Write header
            file << "# time";
            for (size_t i = 0; i < m_cvs.size(); ++i) {
                file << "  cv" << (i + 1) << "_center";
            }
            for (size_t i = 0; i < m_cvs.size(); ++i) {
                file << "  sigma" << (i + 1);
            }
            file << "  height  bias_factor\n";
        }

        // Write last deposited Gaussian
        if (!m_gaussians.empty()) {
            const auto& hill = m_gaussians.back();
            file << hill.time;
            for (double center : hill.cv_centers) {
                file << "  " << center;
            }
            for (double sigma : hill.cv_sigmas) {
                file << "  " << sigma;
            }
            file << "  " << hill.height << "  " << hill.well_tempered_factor << "\n";
        }

        file.close();
    }

    /**
     * @brief Enable grid-based acceleration
     *
     * @param bins_per_dim Number of grid bins per CV dimension (default: 100)
     * @param auto_enable Auto-enable when Gaussian count > threshold (default: true)
     * @param threshold Gaussian count to trigger auto-enable (default: 500)
     *
     * USAGE:
     * ```cpp
     * engine.enableGridAcceleration(100, true, 500);
     * ```
     *
     * PERFORMANCE:
     * - Without grid: O(N_gaussians) per MD step
     * - With grid: O(1) per MD step after setup
     * - Break-even: ~100-500 Gaussians (depends on CV dimensionality)
     *
     * MEMORY:
     * - 1D: ~1 KB (100 bins)
     * - 2D: ~80 KB (100×100)
     * - 3D: ~8 MB (100×100×100)
     */
    void enableGridAcceleration(int bins_per_dim = 100, bool auto_enable = true, int threshold = 500) {
        m_grid_bins = bins_per_dim;
        m_grid_auto = auto_enable;
        m_grid_threshold = threshold;
        m_use_grid = !auto_enable;  // If not auto, enable immediately
    }

    /**
     * @brief Check if grid acceleration is active
     */
    bool isGridActive() const {
        return m_use_grid && m_grid_initialized;
    }

    /**
     * @brief Get grid memory usage in bytes
     */
    size_t getGridMemoryUsage() const {
        if (m_grid) {
            return m_grid->getMemoryUsage();
        }
        return 0;
    }

private:
    std::vector<CVPtr> m_cvs;             ///< Collective variables
    std::vector<double> m_cv_sigmas;      ///< Gaussian widths (one per CV)
    std::vector<GaussianHill> m_gaussians; ///< Deposited Gaussians
    bool m_use_well_tempered;             ///< Use WT-MTD vs. standard MTD
    double m_delta_T;                     ///< Bias temperature for WT-MTD (K)
    double m_initial_height;              ///< Initial Gaussian height (kJ/mol)
    double m_total_bias;                  ///< Current bias energy (kJ/mol)
    int m_num_gaussians;                  ///< Number of deposited Gaussians

    // Grid acceleration
    std::unique_ptr<GridAccelerator> m_grid; ///< Grid-based bias accelerator
    bool m_use_grid;                      ///< Enable grid acceleration
    bool m_grid_auto;                     ///< Auto-enable grid at threshold
    int m_grid_threshold;                 ///< Gaussian count to enable grid
    int m_grid_bins;                      ///< Grid bins per CV dimension
    bool m_grid_initialized;              ///< Grid initialization flag

    /**
     * @brief Initialize grid with automatic bounds determination
     *
     * Claude Generated (November 2025)
     *
     * ALGORITHM:
     * 1. Determine CV ranges from deposited Gaussians
     * 2. Add ±3σ padding for extrapolation
     * 3. Create grid with m_grid_bins per dimension
     * 4. Populate grid with all existing Gaussians
     */
    void initializeGrid() {
        if (m_cvs.empty() || m_gaussians.empty()) {
            return;
        }

        // Create grid
        m_grid = std::make_unique<GridAccelerator>(m_cvs.size());

        // Determine bounds from Gaussian centers ±3σ
        for (size_t cv_idx = 0; cv_idx < m_cvs.size(); ++cv_idx) {
            double min_val = std::numeric_limits<double>::max();
            double max_val = std::numeric_limits<double>::lowest();

            for (const auto& hill : m_gaussians) {
                double center = hill.cv_centers[cv_idx];
                double sigma = hill.cv_sigmas[cv_idx];
                min_val = std::min(min_val, center - 3.0 * sigma);
                max_val = std::max(max_val, center + 3.0 * sigma);
            }

            // Add 10% padding for future Gaussians
            double range = max_val - min_val;
            min_val -= 0.1 * range;
            max_val += 0.1 * range;

            m_grid->setGridBounds(cv_idx, min_val, max_val, m_grid_bins);
        }

        // Initialize grid memory
        m_grid->initialize();

        // Populate grid with all existing Gaussians
        for (const auto& hill : m_gaussians) {
            m_grid->addGaussianToGrid(hill.cv_centers,
                                     hill.height * hill.well_tempered_factor,
                                     hill.cv_sigmas);
        }

        m_grid_initialized = true;
    }

    /**
     * @brief Compute bias using direct Gaussian summation (original method)
     *
     * Claude Generated (November 2025)
     */
    std::pair<double, Geometry> computeBiasDirect(const Molecule& mol,
                                                  const std::vector<double>& current_cv_values) {
        double bias_energy = 0.0;
        Geometry bias_gradient = Eigen::MatrixXd::Zero(mol.AtomCount(), 3);

        // Loop over all Gaussians
        for (const auto& hill : m_gaussians) {
            // Compute total exponential (product of Gaussians along each CV)
            double total_exp = 1.0;
            bool skip_gaussian = false;

            for (size_t alpha = 0; alpha < m_cvs.size(); ++alpha) {
                double ds = current_cv_values[alpha] - hill.cv_centers[alpha];

                // Early termination: if >6σ away, contribution is negligible
                if (std::abs(ds) > 6.0 * hill.cv_sigmas[alpha]) {
                    skip_gaussian = true;
                    break;
                }

                double sigma_sq = hill.cv_sigmas[alpha] * hill.cv_sigmas[alpha];
                total_exp *= std::exp(-ds * ds / (2.0 * sigma_sq));
            }

            if (skip_gaussian) continue;

            // Add to bias energy
            bias_energy += hill.height * total_exp * hill.well_tempered_factor;

            // Compute gradient contribution for each CV
            for (size_t alpha = 0; alpha < m_cvs.size(); ++alpha) {
                double ds = current_cv_values[alpha] - hill.cv_centers[alpha];
                double sigma_sq = hill.cv_sigmas[alpha] * hill.cv_sigmas[alpha];

                // ∂V/∂s_α = -w * exp(...) * (s_α - s_α,i) / σ_α²
                double dV_ds = -hill.height * total_exp * hill.well_tempered_factor * ds / sigma_sq;

                // Chain rule: ∇_r V += (∂V/∂s_α) * ∇_r s_α
                Geometry cv_gradient = m_cvs[alpha]->gradient(mol);
                bias_gradient += dV_ds * cv_gradient;
            }
        }

        return {bias_energy, bias_gradient};
    }

    /**
     * @brief Compute bias using grid interpolation (fast method)
     *
     * Claude Generated (November 2025)
     */
    std::pair<double, Geometry> computeBiasFromGrid(const Molecule& mol,
                                                    const std::vector<double>& current_cv_values) {
        // Interpolate bias energy from grid
        double bias_energy = m_grid->interpolateBias(current_cv_values);

        // Interpolate ∂V/∂s_α from grid
        Eigen::VectorXd dV_ds = m_grid->interpolateGradient(current_cv_values);

        // Apply chain rule: ∇_r V = Σ_α (∂V/∂s_α) * ∇_r s_α
        Geometry bias_gradient = Eigen::MatrixXd::Zero(mol.AtomCount(), 3);
        for (size_t alpha = 0; alpha < m_cvs.size(); ++alpha) {
            Geometry cv_gradient = m_cvs[alpha]->gradient(mol);
            bias_gradient += dV_ds(alpha) * cv_gradient;
        }

        return {bias_energy, bias_gradient};
    }
};

} // namespace CV
