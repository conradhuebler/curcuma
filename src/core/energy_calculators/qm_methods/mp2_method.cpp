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
#include "src/core/units.h"

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
    // RI-MP2 variants (Claude Generated, Oct 2026): ri-mp2, scs-ri-mp2, sos-ri-mp2
    m_options.ri = m_variant.find("ri-mp2") != std::string::npos;
    m_options.ri_aux = lowerString(ec, "ri_mp2_aux", m_options.ri_aux);
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
    if (gradient && m_options.ri) {
        // No RI-MP2 gradient. A zero or an HF gradient would silently optimise the wrong surface, so a
        // gradient request is an error.
        m_error = true;
        m_error_message = m_variant + " has no analytic gradient (energy only; use -sp, or mp2 / hf for -opt/-md/-hessian)";
        CurcumaLogger::error(m_error_message);
        return std::numeric_limits<double>::quiet_NaN();
    }

    const bool scs = m_variant.rfind("scs-", 0) == 0, sos = m_variant.rfind("sos-", 0) == 0;
    // Spin-component scaling of the correlation energy whose gradient is formed (Claude Generated, Oct 2026)
    const double c_os = scs ? 1.2 : sos ? 1.3 : 1.0;
    const double c_ss = scs ? 1.0 / 3.0 : sos ? 0.0 : 1.0;
    m_e_hf = m_engine->Calculation(false);
    if (!m_engine->scfConverged()) {
        m_error = true;
        m_error_message = "MP2: the RHF SCF did not converge";
        return m_e_hf;
    }
    if (gradient) {
        // Canonical MP2 gradient (Claude Generated, Oct 2026): energy and gradient from one pass
        // (QMEngine::computeMP2Gradient, qm_mp2_grad.cpp). The engine works in Eh/Bohr.
        if (!m_engine->computeMP2Gradient(m_options, c_os, c_ss, m_result)) {
            m_error = true;
            m_error_message = m_engine->gradientError().empty() ? std::string("MP2 gradient failed") : m_engine->gradientError();
            m_gradient = Matrix::Zero(m_molecule.AtomCount(), 3);
            return m_e_hf;
        }
        m_gradient = m_engine->gradientBohr() * CurcumaUnit::Length::ANGSTROM_TO_BOHR;  // -> Eh/Angstrom
    } else if (!(m_options.ri ? m_engine->computeRIMP2(m_options, m_result) : m_engine->computeMP2(m_options, m_result))) {
        m_error = true;
        m_error_message = "MP2: correlation energy failed (see the error above)";
        return m_e_hf;
    }
    const double e_corr = scs ? m_result.e_scs() : sos ? m_result.e_sos() : m_result.e_corr();
    m_e_total = m_e_hf + e_corr;

    const int v = CurcumaLogger::get_verbosity();
    if (v >= 1) {
        if (m_options.ri) {
            // Label carries the auxiliary set: JK/J fitting sets used as correlation fitting basis (no def2-SVP/C)
            CurcumaLogger::result(fmt::format("native {} -- RHF reference ({}), RI-MP2 (aux {}; {} fitted functions; "
                                              "method variant, not RI-MP2/def2-SVP/C) ({} basis functions, {} frozen "
                                              "core orbital(s), {} active occupied, {} virtual)",
                                              m_variant, m_result.scf_ri.empty() ? "exact integrals" : m_result.scf_ri,
                                              m_result.aux_name, m_result.n_aux, m_result.n_basis, m_result.n_frozen,
                                              m_result.n_occ_active, m_result.n_virt));
        } else {
            CurcumaLogger::result(fmt::format("native {} -- RHF reference, canonical MP2 ({} basis functions, {} frozen "
                                              "core orbital(s), {} active occupied, {} virtual)",
                                              m_variant, m_result.n_basis, m_result.n_frozen, m_result.n_occ_active,
                                              m_result.n_virt));
        }
        CurcumaLogger::result(fmt::format("  E(HF)          = {:.12f} Eh{}", m_e_hf,
                                          m_result.scf_ri.empty() ? "" : "  (" + m_result.scf_ri + " SCF)"));
        CurcumaLogger::result(fmt::format("  E(corr, {:7s}) = {:.12f} Eh", m_variant, e_corr));
        if (v >= 2) {
            CurcumaLogger::result(fmt::format("  E(MP2 OS)      = {:.12f} Eh", m_result.e_os));
            CurcumaLogger::result(fmt::format("  E(MP2 SS)      = {:.12f} Eh", m_result.e_ss));
            CurcumaLogger::result(fmt::format("  E(MP2 corr)    = {:.12f} Eh", m_result.e_corr()));
            CurcumaLogger::result(fmt::format("  E(SCS-MP2 corr)= {:.12f} Eh  (1.2 OS + 1/3 SS)", m_result.e_scs()));
            CurcumaLogger::result(fmt::format("  E(SOS-MP2 corr)= {:.12f} Eh  (1.3 OS)", m_result.e_sos()));
            if (m_result.gradient)
                CurcumaLogger::info(fmt::format(
                    "MP2 gradient timing: amplitudes {:.1f} ms, back-transformation {:.1f} ms, Z-vector {:.1f} ms ({} iterations, "
                    "max residual {:.1e}), densities {:.1f} ms, derivative integrals 2e {:.1f} ms, 1e {:.1f} ms ({})",
                    m_result.t_grad_amplitudes_ms, m_result.t_grad_backtrans_ms, m_result.t_grad_zvector_ms,
                    m_result.zvector_iterations, m_result.zvector_residual, m_result.t_grad_density_ms,
                    m_result.t_grad_deriv2e_ms, m_result.t_grad_deriv1e_ms, m_result.info));
            else if (m_options.ri)
                CurcumaLogger::info(fmt::format("RI-MP2 timing: RI setup {:.1f} ms, B transformation {:.1f} ms, energy {:.1f} ms ({})",
                                                m_result.t_ri_setup_ms, m_result.t_transform_ms, m_result.t_energy_ms,
                                                m_result.info));
            else
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
