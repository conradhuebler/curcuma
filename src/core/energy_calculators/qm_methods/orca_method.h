/**
 * @file orca_method.h
 * @brief ORCA method wrapper for ComputationalMethod interface
 *
 * Wraps the OrcaInterface as a ComputationalMethod so that ORCA composite
 * methods (HF-3c, B97-3c, r2SCAN-3c, PBEh-3c) and custom ORCA inputs are
 * available via MethodFactory::create().
 *
 * Supported method names:
 *   - hf-3c, b97-3c, r2scan-3c, pbeh-3c  (composite methods, auto-input)
 *   - orca                                (custom keywords via -orca_keywords)
 *
 * Claude Generated - June 2026
 */

#pragma once

#include "../computational_method.h"
#include "orcainterface.h"
#include "src/core/config_manager.h"

#include <memory>
#include <string>

/**
 * @brief ORCA method wrapper for ComputationalMethod interface.
 *
 * This wrapper adapts OrcaInterface to the unified ComputationalMethod
 * interface. ORCA is invoked as an external process; energy and gradients
 * are parsed from JSON or text output.
 */
class OrcaMethod : public ComputationalMethod {
public:
    /**
     * @param method_name One of: "hf-3c", "b97-3c", "r2scan-3c", "pbeh-3c", "orca"
     * @param config JSON configuration (orca-specific parameters)
     */
    OrcaMethod(const std::string& method_name, const json& config = json{});
    virtual ~OrcaMethod() = default;

    // =================================================================================
    // ComputationalMethod interface
    // =================================================================================

    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override;
       Position getDipole() const override;
    bool hasGradient() const override { return true; }

    std::string getMethodName() const override { return m_method_name; }
    bool isThreadSafe() const override { return false; }  // shell process is not thread-safe
    void setThreadCount(int threads) override;

    void setParameters(const json& params) override;
    json getParameters() const override;
    bool hasError() const override;
    void clearError() override;
    std::string getErrorMessage() const override;

    json getEnergyDecomposition() const override;

    // =================================================================================
    // ORCA-specific helpers
    // =================================================================================

    /**
     * @brief Set the electronic charge (overrides Molecule default).
     */
    void setCharge(int charge) { m_charge_override = charge; }

    /**
     * @brief Set the spin multiplicity (overrides Molecule default).
     */
    void setMultiplicity(int mult) { m_multiplicity_override = mult; }

    /**
     * @brief Set custom ORCA simple-input keywords (e.g. "PM3 PAL4").
     * Curcuma prepends '!'. Overrides orca_method and orca_extra_keywords.
     */
    void setCustomKeywords(const std::string& keywords);

    /**
     * @brief Get the last ORCA input file path written to disk.
     */
    std::string getInputPath() const;

    /**
     * @brief Get the last ORCA output file path.
     */
    std::string getOutputPath() const;

    /**
     * @brief Check if orca executable is available.
     */
    static bool isAvailable();

    /**
     * @brief Get list of supported method names.
     */
    static std::vector<std::string> getSupportedMethods();

private:
    std::unique_ptr<OrcaInterface> m_orca;
    std::string m_method_name;       ///< Canonical method name (e.g. "hf-3c")
    std::string m_orca_keyword;    ///< ORCA keyword for input (e.g. "HF-3c")
    Mol m_molecule;
    bool m_calculation_done = false;
    double m_last_energy = 0.0;

    int m_charge_override = 0;       ///< 0 = use molecule charge
    int m_multiplicity_override = 1; ///< 1 = singlet default

    /**
     * @brief Map canonical method name to ORCA keyword.
     */
    static std::string methodToOrcaKeyword(const std::string& method);

    /**
     * @brief Initialize OrcaInterface with current parameters.
     */
    bool initializeOrca();

    /**
     * @brief Build or update molecule from geometry matrix.
     */
    void updateMoleculeGeometry(const Matrix& geometry);
};
