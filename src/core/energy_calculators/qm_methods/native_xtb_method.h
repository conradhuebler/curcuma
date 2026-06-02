/*
 * <Native xTB Method Wrapper — unified GFN1 / GFN2>
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated: a single ComputationalMethod adapter for the native
 * curcuma::xtb::XTB solver, parametrized by MethodType. It replaces the former
 * GFN1Method / GFN2Method wrappers, which were near-identical delegation to the
 * same XTB instance — GFN2 additionally carried try/catch error handling,
 * element validation, and orbital accessors. Those are kept here for both
 * methods (all backed by XTB, which supports GFN1 and GFN2 alike). The only
 * per-method differences are the default config (D3 vs D4) and the method name.
 */

#pragma once

#include "../computational_method.h"
#include "xtb_native.h"
#include "xtb_fragment_scf.h"
#include "src/core/molecule.h"

#include <memory>
#include <string>

/**
 * @brief Unified ComputationalMethod wrapper for native GFN1-xTB / GFN2-xTB.
 *
 * Construct with the MethodType; the factory passes GFN1 or GFN2. All
 * ComputationalMethod calls delegate to the wrapped curcuma::xtb::XTB.
 *
 * Claude Generated.
 */
class NativeXtbMethod : public ComputationalMethod {
public:
    explicit NativeXtbMethod(curcuma::xtb::MethodType method, const json& config = json{});
    ~NativeXtbMethod() override = default;

    // ---- ComputationalMethod interface ------------------------------------
    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override { return Vector{}; }     // not computed
    Position getDipole() const override { return Position{0.0, 0.0, 0.0}; }  // TODO
    bool hasGradient() const override { return true; }

    std::string getMethodName() const override;
    bool isThreadSafe() const override { return true; }
    // Forward the global -threads budget to the solver as its intra-molecule
    // thread count. The solver auto-gates it (serial when run under molecule-level
    // parallelism or for small systems), so this is safe to set unconditionally.
    void setThreadCount(int threads) override
    {
        m_thread_count = std::max(1, threads);
        if (m_xtb) m_xtb->setIntraThreads(m_thread_count);
        if (m_c1_driver) m_c1_driver->setIntraThreads(m_thread_count);
    }

    void setParameters(const json& params) override;
    json getParameters() const override { return m_parameters; }

    bool hasError() const override { return m_has_error; }
    void clearError() override { m_has_error = false; m_error_message.clear(); }
    std::string getErrorMessage() const override { return m_error_message; }

    // Optional advanced (QM) features — delegate to the native solver.
    Vector getOrbitalEnergies() const override;
    int getNumElectrons() const override;
    json getEnergyDecomposition() const override;
    bool saveToFile(const std::string& filename) const override;

    // ---- SCF warm-start / iterative-mode (Claude Generated) ---------------
    void setWarmStart(bool on) override;
    void setIterativeMode(bool on) override;

    // Access the underlying dense native solver (null when a large_system_mode
    // driver is active). Used by the GPU wrapper to install an external
    // eigensolver hook. Claude Generated (GPU port).
    curcuma::xtb::XTB* solver() { return m_xtb.get(); }

    // ---- Native xTB extras (not in the base interface) --------------------
    Matrix getMolecularOrbitals() const;
    double getHOMOLUMOGap() const;
    double getHOMOEnergy() const;
    double getLUMOEnergy() const;
    Vector getCoordinationNumbers() const;

private:
    curcuma::xtb::MethodType m_method;
    std::unique_ptr<curcuma::xtb::XTB> m_xtb;
    // C1 large-system driver (opt-in, c1_mode != none). When active, all
    // calculate/gradient/charge calls delegate here instead of m_xtb. Built in
    // setMolecule once the c1_mode is known. Claude Generated.
    std::unique_ptr<curcuma::xtb::FragmentScfDriver> m_c1_driver;
    Mol m_molecule;
    bool m_calculation_done = false;
    double m_last_energy = 0.0;

    // Read c1_mode from the controller ("xtb" scope, top-level fallback).
    std::string c1ModeString() const;

    // Push the controller settings into the native solver: the D4 charge-response
    // source (harmless for GFN1, which uses D3) and the SCF-convergence settings
    // (mode/guess/damping/DIIS/level-shift) via curcuma::xtb::applyXtbScfConfig.
    void applyConfig();
    void handleError(const std::string& operation);

    static json getDefaultConfig(curcuma::xtb::MethodType method);
    static bool isMoleculeSupported(const Mol& mol);  // native xTB params cover Z = 1..86
};
