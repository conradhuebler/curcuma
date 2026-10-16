/*
 * <Native canonical MP2 method wrapper (WP11) -- implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Oct 2026, WP11). AI-generated, machine-tested.
 *
 * This program is free software under GPL-3.0
 */

#include "mp2_method.h"
#include "qm_method.h"

#include "src/core/curcuma_logger.h"

#include <fmt/format.h>

#include <cmath>
#include <limits>

namespace {
std::string lowerString(const json& j, const char* key, const std::string& fallback)
{
    if (j.contains(key)) {
        if (j[key].is_string()) return j[key].get<std::string>();
        if (j[key].is_boolean()) return j[key].get<bool>() ? "on" : "off";
    }
    return fallback;
}

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
}  // namespace

MP2Method::MP2Method(const std::string& variant, const json& config)
    : m_variant(variant)
{
    const json ec = QMMethod::engineConfig(config);  // registry qm defaults < top level < dft < qm
    m_engine = std::make_unique<QMEngine>(QMFunctional::HF, ec);
    m_options.frozen_core = readBool(ec, "frozen_core", false);
    if (ec.contains("mp2_max_memory_mb") && ec["mp2_max_memory_mb"].is_number())
        m_options.max_memory_mb = ec["mp2_max_memory_mb"].get<double>();
    if (ec.contains("mp2_batch_mb") && ec["mp2_batch_mb"].is_number())
        m_options.batch_mb = ec["mp2_batch_mb"].get<double>();
    m_options.gpu_mp2 = lowerString(ec, "gpu_mp2", "auto");
}

bool MP2Method::setMolecule(const Mol& mol)
{
    m_molecule = mol;
    m_error = false;
    m_engine->resetForNewMolecule();  // a reused object must not keep the old molecule
    if (!m_engine->QMInterface::InitialiseMolecule(mol)) {
        m_error = true;
        m_error_message = "MP2: QM engine initialisation failed";
        return false;
    }
    return true;
}

bool MP2Method::updateGeometry(const Matrix& geometry)
{
    m_molecule.m_geometry = geometry;
    return m_engine->UpdateMolecule(geometry);
}

double MP2Method::calculateEnergy(bool gradient)
{
    m_error = false;
    m_error_message.clear();
    if (gradient) {
        // No MP2 gradient (Z-vector / relaxed density, own work package). A zero or an HF gradient
        // would silently optimise the wrong surface, so a gradient request is an error.
        m_error = true;
        m_error_message = m_variant + " has no analytic gradient (energy only; use -sp, or hf for -opt/-md/-hessian)";
        CurcumaLogger::error(m_error_message);
        return std::numeric_limits<double>::quiet_NaN();
    }

    m_e_hf = m_engine->Calculation(false);
    if (!m_engine->scfConverged()) {
        m_error = true;
        m_error_message = "MP2: the RHF SCF did not converge";
        return m_e_hf;
    }
    if (!m_engine->computeMP2(m_options, m_result)) {
        m_error = true;
        m_error_message = "MP2: correlation energy failed (see the error above)";
        return m_e_hf;
    }
    const double e_corr = (m_variant == "scs-mp2") ? m_result.e_scs()
        : (m_variant == "sos-mp2")                  ? m_result.e_sos()
                                                    : m_result.e_corr();
    m_e_total = m_e_hf + e_corr;

    const int v = CurcumaLogger::get_verbosity();
    if (v >= 1) {
        CurcumaLogger::result(fmt::format("native {} -- RHF reference, canonical MP2 ({} basis functions, {} frozen "
                                          "core orbital(s), {} active occupied, {} virtual)",
                                          m_variant, m_result.n_basis, m_result.n_frozen, m_result.n_occ_active,
                                          m_result.n_virt));
        CurcumaLogger::result(fmt::format("  E(HF)          = {:.12f} Eh", m_e_hf));
        CurcumaLogger::result(fmt::format("  E(corr, {:7s}) = {:.12f} Eh", m_variant, e_corr));
        if (v >= 2) {
            CurcumaLogger::result(fmt::format("  E(MP2 OS)      = {:.12f} Eh", m_result.e_os));
            CurcumaLogger::result(fmt::format("  E(MP2 SS)      = {:.12f} Eh", m_result.e_ss));
            CurcumaLogger::result(fmt::format("  E(MP2 corr)    = {:.12f} Eh", m_result.e_corr()));
            CurcumaLogger::result(fmt::format("  E(SCS-MP2 corr)= {:.12f} Eh  (1.2 OS + 1/3 SS)", m_result.e_scs()));
            CurcumaLogger::result(fmt::format("  E(SOS-MP2 corr)= {:.12f} Eh  (1.3 OS)", m_result.e_sos()));
            CurcumaLogger::info(fmt::format("MP2 timing: integrals {:.1f} ms, transformation {:.1f} ms, energy {:.1f} ms ({})",
                                            m_result.t_integrals_ms, m_result.t_transform_ms, m_result.t_energy_ms,
                                            m_result.info));
        }
        CurcumaLogger::result(fmt::format("  E({})  = {:.12f} Eh", m_variant, m_e_total));
        CurcumaLogger::energy_abs(m_e_total, "total " + m_variant + " energy");
    }
    return m_e_total;
}

json MP2Method::getEnergyDecomposition() const
{
    json d;
    d["hf"] = m_e_hf;
    d["mp2_os"] = m_result.e_os;
    d["mp2_ss"] = m_result.e_ss;
    d["mp2_corr"] = m_result.e_corr();
    d["scs_mp2_corr"] = m_result.e_scs();
    d["sos_mp2_corr"] = m_result.e_sos();
    d["total"] = m_e_total;
    return d;
}
