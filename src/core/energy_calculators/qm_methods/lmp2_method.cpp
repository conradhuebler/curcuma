/*
 * <Native local PNO-RI-MP2 method wrapper (WP15) -- implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Oct 2026, WP15). AI-generated, machine-tested.
 *
 * This program is free software under GPL-3.0
 */

#include "lmp2_method.h"
#include "qm_method.h"

#include "src/core/curcuma_logger.h"

#include <fmt/format.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdlib>
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

std::string readString(const json& j, const char* key, const std::string& fallback, bool lower)
{
    if (!j.contains(key) || !j[key].is_string()) return fallback;
    std::string s = j[key].get<std::string>();
    if (lower) std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return (char)std::tolower(c); });
    return s;
}
}  // namespace

LMP2Method::LMP2Method(const json& config)
{
    const json ec = QMMethod::engineConfig(config);  // registry qm defaults < top level < dft < qm
    m_engine = std::make_unique<QMEngine>(QMFunctional::HF, ec);
    qmint::LMP2Options& o = m_options;
    o.frozen_core = readBool(ec, "frozen_core", false);
    o.ri_aux = readString(ec, "ri_mp2_aux", o.ri_aux, false);
    o.localizer = readString(ec, "lmp2_localizer", o.localizer, true);
    o.localizer_algo = readString(ec, "lmp2_localizer_algo", o.localizer_algo, true);
    o.localizer_screen = readNumber(ec, "lmp2_localizer_screen", o.localizer_screen);
    o.bp_threshold = readNumber(ec, "lmp2_bp_threshold", o.bp_threshold);
    o.domain_radius = readNumber(ec, "lmp2_domain_radius", o.domain_radius);
    o.pno_cut = readNumber(ec, "lmp2_pno_cut", o.pno_cut);
    o.pno_correction = readBool(ec, "lmp2_pno_correction", o.pno_correction);
    o.pair_cut = readNumber(ec, "lmp2_pair_cut", o.pair_cut);
    o.distant_pairs = readString(ec, "lmp2_distant_pairs", o.distant_pairs, true);
    o.fock_cut = readNumber(ec, "lmp2_fock_cut", o.fock_cut);
    o.e_threshold = readNumber(ec, "lmp2_threshold", o.e_threshold);
    o.r_threshold = readNumber(ec, "lmp2_residual_threshold", o.r_threshold);
    o.max_iter = (int)readNumber(ec, "lmp2_max_iter", o.max_iter);
    o.diis = (int)readNumber(ec, "lmp2_diis", o.diis);
    o.max_memory_mb = readNumber(ec, "lmp2_max_memory_mb", o.max_memory_mb);
    o.fit = readString(ec, "lmp2_fit", o.fit, true);
    o.fit_domain = readString(ec, "lmp2_fit_domain", o.fit_domain, true);
    o.fit_radius = readNumber(ec, "lmp2_fit_radius", o.fit_radius);
    o.ao_cut = readNumber(ec, "lmp2_ao_cut", o.ao_cut);
    m_check = readBool(ec, "lmp2_check", false);
}

bool LMP2Method::setMolecule(const Mol& mol)
{
    m_molecule = mol;
    m_error = false;
    m_engine->resetForNewMolecule();  // a reused object must not keep the old molecule
    if (!m_engine->QMInterface::InitialiseMolecule(mol)) {
        m_error = true;
        m_error_message = "LMP2: QM engine initialisation failed";
        return false;
    }
    return true;
}

bool LMP2Method::updateGeometry(const Matrix& geometry)
{
    m_molecule.m_geometry = geometry;
    return m_engine->UpdateMolecule(geometry);
}

double LMP2Method::calculateEnergy(bool gradient)
{
    m_error = false;
    m_error_message.clear();
    if (gradient) {
        // No LMP2 gradient. A zero or an HF gradient would silently optimise the wrong surface.
        m_error = true;
        m_error_message = "lmp2 has no analytic gradient (energy only; use -sp, or mp2 / hf for -opt/-md/-hessian)";
        CurcumaLogger::error(m_error_message);
        return std::numeric_limits<double>::quiet_NaN();
    }
    const auto t0 = std::chrono::steady_clock::now();
    m_e_hf = m_engine->Calculation(false);
    const double t_scf = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    if (!m_engine->scfConverged()) {
        m_error = true;
        m_error_message = "LMP2: the RHF SCF did not converge";
        return std::numeric_limits<double>::quiet_NaN();
    }
    const auto t1 = std::chrono::steady_clock::now();
    if (!m_engine->computeLMP2(m_options, m_result)) {
        m_error = true;
        m_error_message = "LMP2: correlation energy failed (see the error above)";
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double t_lmp2 = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t1).count();
    const qmint::LMP2Result& r = m_result;
    m_e_total = m_e_hf + r.e_corr();

    // Diagnostic (env CURCUMA_LMP2_LOC_COMPARE): run LMP2 again on the same SCF with the other localisation algorithm
    // (fast <-> classic) and, if the value is a comma-separated list of numbers, with the fast route at each of these
    // lmp2_localizer_screen values in [0, 1); prints E_corr and the functional of every run. Claude Generated (Oct 2026).
    std::vector<std::string> loc_compare;
    if (const char* env = std::getenv("CURCUMA_LMP2_LOC_COMPARE")) {
        auto line = [&](const std::string& tag, const qmint::LMP2Result& r2) {
            return fmt::format("LMP2 localiser comparison: {:<22} E_corr {:.12f} F {:.12f} sweeps {:3d} t_loc {:9.1f} ms  dE_corr {:10.3e} "
                               "dF {:10.3e}",
                               tag, r2.e_corr(), r2.localization_functional, r2.localization_sweeps, r2.t_localize_ms,
                               r2.e_corr() - m_result.e_corr(), r2.localization_functional - m_result.localization_functional);
        };
        loc_compare.push_back(line(fmt::format("{} screen {:g} (run)", m_options.localizer_algo, m_options.localizer_screen), m_result));
        std::vector<qmint::LMP2Options> runs;
        qmint::LMP2Options o2 = m_options;
        o2.localizer_algo = m_options.localizer_algo == "fast" ? "classic" : "fast";
        runs.push_back(o2);
        std::string list(env);
        size_t pos = 0;
        while (pos < list.size()) {
            const size_t comma = list.find(',', pos);
            const std::string item = list.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos);
            pos = comma == std::string::npos ? list.size() : comma + 1;
            try {
                size_t used = 0;
                const double sc = std::stod(item, &used);
                if (used == item.size() && sc >= 0.0 && sc < 1.0) {  // "1" alone = only the other route
                    qmint::LMP2Options o3 = m_options;
                    o3.localizer_algo = "fast";
                    o3.localizer_screen = sc;
                    runs.push_back(o3);
                }
            } catch (...) {
            }
        }
        for (const auto& oo : runs) {
            qmint::LMP2Result r2;
            if (m_engine->computeLMP2(oo, r2))
                loc_compare.push_back(line(fmt::format("{} screen {:g}", oo.localizer_algo, oo.localizer_screen), r2));
            else
                loc_compare.push_back("LMP2 localiser comparison: run failed");
        }
    }

    double t_check = 0.0;
    if (m_check) {
        qmint::MP2Options mo;
        mo.frozen_core = m_options.frozen_core;
        mo.ri = true;
        mo.ri_aux = m_options.ri_aux;
        mo.gpu_mp2 = "off";
        const auto t2 = std::chrono::steady_clock::now();
        if (!m_engine->computeRIMP2(mo, m_check_result)) {
            m_error = true;
            m_error_message = "LMP2: the canonical RI-MP2 check failed (see the error above)";
            return std::numeric_limits<double>::quiet_NaN();
        }
        t_check = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t2).count();
    }

    const int v = CurcumaLogger::get_verbosity();
    if (v >= 1) {
        CurcumaLogger::result(fmt::format("native lmp2 -- RHF reference ({}), local PNO-RI-MP2 (aux {}; {} fitted functions; "
                                          "method variant) ({} atoms, {} basis functions, {} frozen core orbital(s), {} "
                                          "active occupied, {} virtual)",
                                          r.scf_ri.empty() ? "exact integrals" : r.scf_ri, r.aux_name, r.n_aux, r.n_atoms,
                                          r.n_basis, r.n_frozen, r.n_occ_active, r.n_virt));
        CurcumaLogger::result(fmt::format("  settings: localizer {}, BP {:g} + {:g} Bohr, T_CutPNO {:g}, pair cut {:g} Bohr ({}), "
                                          "Fock cut {:g}",
                                          m_options.localizer, m_options.bp_threshold, m_options.domain_radius, m_options.pno_cut,
                                          m_options.pair_cut, m_options.distant_pairs, m_options.fock_cut));
        if (m_options.fit == "local")
            CurcumaLogger::result(fmt::format("  fit: local, {} domains, radius {:g} Bohr, ao_cut {:g} (integral-direct 3-index integrals, "
                                              "no stored RI tensor); mean fit domain {:.0f} aux functions, {} metric factorisations, "
                                              "mean X rows {:.0f}",
                                              m_options.fit_domain, m_options.fit_radius, m_options.ao_cut, r.mean_fit_aux,
                                              r.n_fit_factors, r.mean_x_rows));
        else
            CurcumaLogger::result("  fit: global (stored RI tensor)");
        CurcumaLogger::result(fmt::format("  pairs kept {} of {} (distant {})", r.n_pairs_kept, r.n_pairs_total, r.n_pairs_distant));
        CurcumaLogger::result(fmt::format("  PNOs per kept pair: mean {:.2f}, max {} (PAO pair space mean {:.1f}; domain atoms mean {:.2f})",
                                          r.mean_pno, r.max_pno, r.mean_pao, r.mean_domain_atoms));
        CurcumaLogger::result(fmt::format("  amplitude fraction = {:.6e} (sum npno^2 / (o(o+1)/2 v^2)); {} pair couplings",
                                          r.amp_fraction, r.n_couplings));
        CurcumaLogger::result(fmt::format("  LMP2 iterations = {} (last dE = {:.1e} Eh, max|R| = {:.1e}); localisation sweeps {}",
                                          r.iterations, r.last_de, r.max_residual, r.localization_sweeps));
        CurcumaLogger::result(fmt::format("  E(HF)          = {:.12f} Eh{}", m_e_hf, r.scf_ri.empty() ? "" : "  (" + r.scf_ri + " SCF)"));
        CurcumaLogger::result(fmt::format("  E(LMP2 iter)   = {:.12f} Eh", r.e_iter));
        CurcumaLogger::result(fmt::format("  E(PNO corr)    = {:.12f} Eh", r.e_pno_corr));
        CurcumaLogger::result(fmt::format("  E(distant)     = {:.12f} Eh", r.e_distant));
        CurcumaLogger::result(fmt::format("  E(LMP2 corr)   = {:.12f} Eh", r.e_corr()));
        if (m_check) {
            const double ref = m_check_result.e_corr(), d = r.e_corr() - ref;
            CurcumaLogger::result(fmt::format("  E(RI-MP2 corr) = {:.12f} Eh  (canonical, same aux; check)", ref));
            CurcumaLogger::result(fmt::format("  LMP2 - RI-MP2  = {:.3e} Eh  ({:.3e} Eh/atom, recovery {:.4f} %)", d,
                                              d / std::max(1, r.n_atoms), 100.0 * r.e_corr() / ref));
        }
        if (v >= 2) {
            CurcumaLogger::info(fmt::format("LMP2 timing: SCF {:.1f} ms, LMP2 {:.1f} ms (RI setup {:.1f}, localisation {:.1f}, {} "
                                            "{:.1f}, pairs/PNOs {:.1f}{}, distant {:.1f}, overlaps {:.1f}, "
                                            "iterations {:.1f}){}; memory {:.0f} MB ({})",
                                            t_scf, t_lmp2, r.t_ri_setup_ms, r.t_localize_ms,
                                            m_options.fit == "local" ? "integral-direct (P|i r)" : "B transformation", r.t_transform_ms,
                                            r.t_pairs_ms,
                                            m_options.fit == "local" ? fmt::format(" (fit, thread-summed {:.1f})", r.t_fit_ms) : std::string(),
                                            r.t_distant_ms, r.t_overlap_ms, r.t_iter_ms,
                                            m_check ? fmt::format(", RI-MP2 check {:.1f} ms", t_check) : std::string(),
                                            r.memory_mb, r.info));
        }
        if (v >= 2)
            CurcumaLogger::info(fmt::format("LMP2 localisation ({}, {}): functional {:.12f}, {} sweeps; ms: {}", m_options.localizer,
                                            m_options.localizer_algo, r.localization_functional, r.localization_sweeps,
                                            r.localization_info));
        for (const auto& lc : loc_compare) CurcumaLogger::result(lc);
        CurcumaLogger::result(fmt::format("  E(lmp2)        = {:.12f} Eh", m_e_total));
        CurcumaLogger::energy_abs(m_e_total, "total lmp2 energy");
    }
    return m_e_total;
}

json LMP2Method::getEnergyDecomposition() const
{
    json d;
    d["hf"] = m_e_hf;
    d["lmp2_iter"] = m_result.e_iter;
    d["lmp2_pno_corr"] = m_result.e_pno_corr;
    d["lmp2_distant"] = m_result.e_distant;
    d["lmp2_corr"] = m_result.e_corr();
    if (m_check) d["rimp2_corr_check"] = m_check_result.e_corr();
    d["total"] = m_e_total;
    return d;
}
