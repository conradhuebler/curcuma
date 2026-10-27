/*
 * <Native canonical closed-shell CCSD method wrapper (WP13) -- implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Oct 2026, WP13). AI-generated, machine-tested.
 *
 * This program is free software under GPL-3.0
 */

#include "ccsd_method.h"
#include "qm_citations.h"
#include "qm_method.h"

#include "src/core/curcuma_logger.h"

#include <fmt/format.h>

#include <cmath>
#include <limits>

namespace {
bool readBool(const json& j, const char* key, bool fallback)
{
    if (!j.contains(key)) return fallback;
    const auto& v = j[key];
    if (v.is_boolean()) return v.get<bool>();
    if (v.is_number()) return v.get<double>() != 0.0;
    if (v.is_string()) {
        const std::string s = v.get<std::string>();
        return s == "true" || s == "1" || s == "on" || s == "yes";
    }
    return fallback;
}

double readNumber(const json& j, const char* key, double fallback)
{
    if (j.contains(key) && j[key].is_number()) return j[key].get<double>();
    if (j.contains(key) && j[key].is_string()) {
        try {
            return std::stod(j[key].get<std::string>());
        } catch (...) {
        }
    }
    return fallback;
}
}  // namespace

CCSDMethod::CCSDMethod(const json& config)
{
    const json ec = QMMethod::engineConfig(config);  // registry qm defaults < top level < dft < qm
    m_engine = std::make_unique<QMEngine>(QMFunctional::HF, ec);
    m_engine->setTightOrbitals(true);  // correlation energy is first order in the orbital error (Oct 25 2026)
    m_options.frozen_core = readBool(ec, "frozen_core", false);
    m_options.max_iter = (int)readNumber(ec, "cc_max_iter", m_options.max_iter);
    m_options.e_threshold = readNumber(ec, "cc_threshold", m_options.e_threshold);
    m_options.r_threshold = readNumber(ec, "cc_residual_threshold", m_options.r_threshold);
    m_options.diis = (int)readNumber(ec, "cc_diis", m_options.diis);
    m_options.max_memory_mb = readNumber(ec, "cc_max_memory_mb", m_options.max_memory_mb);
    m_options.ao_max_memory_mb = readNumber(ec, "mp2_max_memory_mb", m_options.ao_max_memory_mb);
    qmcite::ccsd();
}

bool CCSDMethod::setMolecule(const Mol& mol)
{
    m_molecule = mol;
    m_error = false;
    m_engine->resetForNewMolecule();  // a reused object must not keep the old molecule
    if (!m_engine->QMInterface::InitialiseMolecule(mol)) {
        m_error = true;
        m_error_message = "CCSD: QM engine initialisation failed";
        return false;
    }
    return true;
}

bool CCSDMethod::updateGeometry(const Matrix& geometry)
{
    m_molecule.m_geometry = geometry;
    return m_engine->UpdateMolecule(geometry);
}

double CCSDMethod::calculateEnergy(bool gradient)
{
    m_error = false;
    m_error_message.clear();
    if (gradient) {
        // No CCSD gradient. A zero or an HF gradient would silently optimise the wrong surface.
        m_error = true;
        m_error_message = "ccsd has no analytic gradient (energy only; use -sp, or mp2 / hf for -opt/-md/-hessian)";
        CurcumaLogger::error(m_error_message);
        return std::numeric_limits<double>::quiet_NaN();
    }
    m_e_hf = m_engine->Calculation(false);
    if (!m_engine->scfConverged()) {
        m_error = true;
        m_error_message = "CCSD: the RHF SCF did not converge";
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (!m_engine->computeCCSD(m_options, m_result)) {
        m_error = true;
        m_error_message = "CCSD: correlation energy failed (see the error above)";
        return std::numeric_limits<double>::quiet_NaN();
    }
    m_e_total = m_e_hf + m_result.e_corr;

    const int v = CurcumaLogger::get_verbosity();
    if (v >= 1) {
        CurcumaLogger::result(fmt::format("native ccsd -- RHF reference, canonical RCCSD ({} basis functions, {} frozen "
                                          "core orbital(s), {} active occupied, {} virtual)",
                                          m_result.n_basis, m_result.n_frozen, m_result.n_occ_active, m_result.n_virt));
        CurcumaLogger::result(fmt::format("  E(HF)          = {:.12f} Eh", m_e_hf));
        CurcumaLogger::result(fmt::format("  E(MP2 corr)    = {:.12f} Eh  (CCSD iteration 0)", m_result.e_mp2));
        CurcumaLogger::result(fmt::format("  E(CCSD corr)   = {:.12f} Eh", m_result.e_corr));
        CurcumaLogger::result(fmt::format("  CCSD iterations = {} (converged: last dE = {:.1e} Eh, max|R| = {:.1e})",
                                          m_result.iterations, m_result.last_de, m_result.max_residual));
        CurcumaLogger::result(fmt::format("  T1 diagnostic  = {:.8f}", m_result.t1_diagnostic));
        if (v >= 2)
            CurcumaLogger::info(fmt::format("CCSD timing: integrals {:.1f} ms, transformation {:.1f} ms, iterations {:.1f} ms ({})",
                                            m_result.t_integrals_ms, m_result.t_transform_ms, m_result.t_iterations_ms,
                                            m_result.info));
        CurcumaLogger::result(fmt::format("  E(ccsd)        = {:.12f} Eh", m_e_total));
        CurcumaLogger::energy_abs(m_e_total, "total ccsd energy");
    }
    return m_e_total;
}

json CCSDMethod::getEnergyDecomposition() const
{
    json d;
    d["hf"] = m_e_hf;
    d["mp2_corr"] = m_result.e_mp2;
    d["ccsd_corr"] = m_result.e_corr;
    d["t1_diagnostic"] = m_result.t1_diagnostic;
    d["ccsd_iterations"] = m_result.iterations;
    d["total"] = m_e_total;
    return d;
}
