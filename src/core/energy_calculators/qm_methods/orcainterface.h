//
// Created by gerd on 11.12.24.
// Extended by Claude for ComputationalMethod integration (June 2026)
//

#pragma once
#ifndef ORCAINTERFACE_H
#define ORCAINTERFACE_H

#include "src/core/molecule.h"
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

#include <string>
#include <vector>

// Claude Generated 2025: ORCA Parameter Registry - replaces static ORCASettings JSON
BEGIN_PARAMETER_DEFINITION(orca)
    // Method selection
    PARAM(orca_method, String, "hf-3c", "ORCA composite method keyword: HF-3c, B97-3c, r2SCAN-3c, PBEh-3c, or custom.", "Method", {})
    PARAM(orca_keywords, String, "", "Custom ORCA simple-input keywords (e.g. 'PM3 PAL4'). Curcuma prepends '!'. When set, overrides orca_method and orca_extra_keywords.", "Method", {})
    PARAM(orca_nprocs, Int, 1, "Number of CPU cores for ORCA (maps to %pal nprocs).", "Parallel", {})
    PARAM(orca_maxcore, Int, 500, "Max memory per core in MB (maps to %maxcore).", "Parallel", {})
    PARAM(orca_extra_keywords, String, "EnGrad TightSCF", "Additional ORCA simple-input keywords appended to the ! line (e.g. 'EnGrad TightSCF').", "Method", {})
    PARAM(orca_basename, String, "orca_calc", "Basename for ORCA input/output files written to disk.", "Output", {})
    PARAM(orca_executable, String, "", "Full path to the ORCA executable. If empty, the environment variable ORCA_PATH (or ORCA_ROOT) is checked, then 'orca' is searched via PATH.", "Method", {})
    PARAM(orca_allow_cg, Bool, false, "Allow ORCA inputs to include CG (coarse-grained) atoms. ORCA cannot parse element 226 (CG_ELEMENT). Leave false unless the molecule has been pre-converted to all-atom.", "Validation", {})

    // Solvation Models
    PARAM(cpcm_solvent, String, "none", "CPCM solvent name (none, water, acetone, etc.).", "Solvation", {"cpcm_solv"})
    PARAM(cpcm_epsilon, Double, -1.0, "CPCM dielectric constant. -1 = auto-detect from name.", "Solvation", {"cpcm_eps"})
    PARAM(alpb_solvent, String, "none", "ALPB solvent name (none, water, acetone, etc.).", "Solvation", {"alpb_solv"})
    PARAM(alpb_epsilon, Double, -1.0, "ALPB dielectric constant. -1 = auto-detect from name.", "Solvation", {"alpb_eps"})
END_PARAMETER_DEFINITION

/**
 * @brief ORCA external-process interface.
 *
 * Resolves the ORCA executable, writes input files, runs ORCA via popen,
 * and parses energy/gradient/charges from .property.json (preferred) or
 * the .out text (fallback).
 *
 * @note NOT thread-safe. Each thread needs its own instance with a unique
 *       orca_basename, because calculate() writes to m_input_path + ".inp" /
 *       ".out" / ".property.json" on disk. The wrapper OrcaMethod exposes
 *       isThreadSafe() == false; setThreadCount() re-creates the
 *       OrcaInterface, which does NOT guarantee unique basenames (caller's
 *       responsibility).
 */
class OrcaInterface {
public:
    explicit OrcaInterface(const ConfigManager& config);
    explicit OrcaInterface();
    /**
     * @brief Convenience constructor for raw JSON config.
     *        Equivalent to OrcaInterface(ConfigManager("orca", json_config)).
     */
    explicit OrcaInterface(const json& json_config);
    ~OrcaInterface();

    // =================================================================================
    // Original interface (kept for backward compatibility)
    // =================================================================================
    void setInputFile(const std::string& inputFile);
    bool runOrca();
    void readOrcaJSON();
    bool getOrcaJSON();

    // =================================================================================
    // ComputationalMethod support (Claude Generated - June 2026)
    // =================================================================================

    /**
     * @brief Run a user-supplied ORCA input file (legacy CLI path).
     *
     * For `curcuma -orca <input.inp>` style invocation. Skips generateInput()
     * and uses the file at m_input_path as-is.
     *
     * @param injectEnGrad If true, append "EnGrad" to the ! line if missing.
     * @param verbosity Curcuma verbosity (0=silent ... 3=full). Pass -1 for default.
     * @return true on success.
     *
     * @note This is the new public entry point for the legacy -orca CLI flow.
     *       Replaces the setInputFile+runOrca+getOrcaJSON sequence.
     */
    bool runExistingInput(bool injectEnGrad = false, int verbosity = -1);

    // =================================================================================
    // Helpers exposed for unit testing (@public-for-testing)
    // =================================================================================

    /**
     * @brief Build the *xyz block (Claude Generated - June 2026).
     * @public-for-testing Used by test_orca_interface.cpp.
     */
    std::string buildGeometryBlock(const Mol& mol) const;

    /**
     * @brief Parse gradient (natoms x 3) from .property.json or .out text.
     * @public-for-testing Used by test_orca_interface.cpp.
     */
    Matrix parseGradient(int natoms) const;

    /**
     * @brief Parse atomic charges from .property.json or .out text.
     * @public-for-testing Used by test_orca_interface.cpp.
     */
    Vector parseCharges(int natoms) const;

    /**
     * @brief Parse gradient from ORCA text output (fallback).
     * @public-for-testing Used by test_orca_interface.cpp.
     */
    Matrix parseGradientFromText(int natoms) const;

    /**
     * @brief Parse charges from ORCA text output (fallback).
     * @public-for-testing Used by test_orca_interface.cpp.
     */
    Vector parseChargesFromText(int natoms) const;

    /**
     * @brief Read .property.json into m_property_json.
     * @public-for-testing Tests inject a JSON mock via this method.
     */
    bool loadPropertyJSON();

    /**
     * @brief Parse energy from .property.json or .out text.
     * @public-for-testing Used by test_orca_interface.cpp.
     */
    double parseEnergy() const;

    /**
     * @brief Parse dipole moment from .property.json.
     * @public-for-testing Used by test_orca_interface.cpp.
     */
    Position parseDipole() const;

    /**
     * @brief Check if the 'orca' executable is available via system PATH.
     * @return true if orca can be found.
     */
    static bool checkOrcaExecutable();

    /**
     * @brief Resolve the full path to the ORCA executable.
     *
     * Resolution order:
     *   1. Config parameter orca_executable (if set)
     *   2. Environment variable ORCA_PATH (full path to binary)
     *   3. Environment variable ORCA_ROOT (directory, append /orca)
     *   4. "orca" via system PATH
     *
     * For parallel runs ORCA requires the absolute path to the binary.
     *
     * @return Absolute path to the ORCA executable, or empty string if not found.
     */
    std::string findOrcaExecutable() const;

    /**
     * @brief Generate and write an ORCA input file from a Molecule.
     *
     * Builds a complete input with:
     *   ! METHOD [extra_keywords]
     *   %pal nprocs N end
     *   %maxcore M
     *   *xyz charge multiplicity
     *     ...
     *   *
     *
     * If orca_keywords is set, uses ! orca_keywords instead of the default
     * ! method extra_keywords line.
     *
     * @param mol Molecule with geometry, charge, multiplicity.
     * @param method_keyword ORCA method keyword (e.g. "HF-3c").
     * @return true on success.
     */
    bool generateInput(const Mol& mol, const std::string& method_keyword = "");

    /**
     * @brief Overwrite geometry coordinates in an existing ORCA input file.
     *
     * Replaces the *xyz ... * block with current geometry.
     *
     * @note This function performs a non-atomic read-modify-write of the
     *       input file. The file write is the de-facto synchronization
     *       point; do not call from multiple threads on the same instance.
     *
     * @param mol New geometry.
     * @return true on success.
     */
    bool updateGeometryInInput(const Mol& mol);

    /**
     * @brief Run ORCA on the generated input and capture output.
     *
     * Writes energy/gradient/charges to member variables.
     *
     * @param gradient If true, ensure EnGrad keyword is present and parse gradients.
     * @param verbosity Curcuma verbosity level (0-3). Live ORCA output filtering:
     *        verbosity >= 3: full ORCA stdout; verbosity == 2: SCF iterations only.
     * @return Total energy in Hartree, or 0.0 on failure (check hasError()).
     */
    double calculate(bool gradient = false, int verbosity = 0);

    // Property accessors after calculate()
    double getEnergy() const { return m_last_energy; }
    Matrix getGradient() const { return m_last_gradient; }
    Vector getCharges() const { return m_last_charges; }
    Position getDipole() const { return m_last_dipole; }

    bool hasError() const { return m_has_error; }
    std::string getErrorMessage() const { return m_error_message; }
    void clearError() { m_has_error = false; m_error_message.clear(); }

    /**
     * @brief Run orca_2json to produce .property.json, then load it.
     * @public-for-testing Used by test_orca_interface.cpp.
     */
    bool produceAndLoadJSON();

    // File paths
    std::string getInputPath() const { return m_input_path; }
    std::string getOutputPath() const { return m_output_path; }

private:
    // Original members
    std::string inputFilePath;
    std::string outputFilePath;
    json OrcaJSON;
    mutable ConfigManager m_config;

    // Hilfsfunktionen (original)
    bool createInputFile(const std::string& content);
    bool executeOrcaProcess();

    // Citation is handled by the OrcaMethod wrapper (single source of truth).

    /**
     * @brief Run orca via popen and capture stdout/stderr into the .out file
     *        while applying the live-output verbosity filter.
     *
     * Extracted from calculate() so both calculate() and runExistingInput()
     * can share the same process-execution path.
     *
     * @param verbosity Curcuma verbosity (0=silent ... 3=full).
     * @return true on exit-status 0, false otherwise (sets m_has_error).
     */
    bool runProcessAndCapture(int verbosity);

    // State
    std::string m_input_path = "orca.inp";
    std::string m_output_path = "orca.out";
    json m_property_json;

    double m_last_energy = 0.0;
    Matrix m_last_gradient;
    Vector m_last_charges;
    Position m_last_dipole = {0.0, 0.0, 0.0};

    bool m_has_error = false;
    std::string m_error_message;
};

#endif // ORCAINTERFACE_H
