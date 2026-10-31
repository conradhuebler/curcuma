/*
 * <Native double-hybrid density functionals (WP12 stage 1) -- table and implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Oct 2026, WP12 stage 1). AI-generated, machine-tested; human production test pending.
 *
 * This program is free software under GPL-3.0
 */

#include "double_hybrid_method.h"
#include "qm_citations.h"
#include "qm_method.h"

#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/d3param_generator.h"
#include "src/core/units.h"
#include "src/core/energy_calculators/ff_methods/cn_calculator.h"

#include <fmt/format.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>

using qmxc::XCKernel;

namespace {
std::string lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}

std::string readString(const json& j, const char* key, const std::string& fallback)
{
    if (j.contains(key) && j[key].is_string()) return lower(j[key].get<std::string>());
    return fallback;
}

bool readBool(const json& j, const char* key, bool fallback)
{
    if (!j.contains(key)) return fallback;
    const auto& v = j[key];
    if (v.is_boolean()) return v.get<bool>();
    if (v.is_number()) return v.get<double>() != 0.0;
    if (v.is_string()) {
        const std::string s = lower(v.get<std::string>());
        return s == "true" || s == "1" || s == "on" || s == "yes";
    }
    return fallback;
}

const char* dispersionName(DHDispersion d)
{
    return d == DHDispersion::D3BJ ? "D3(BJ)" : d == DHDispersion::D4 ? "D4" : "none";
}

// Source strings (WP12 doc, section "Parameter und Quellen"). "ORCA 6.1" = the values ORCA 6.1.0 prints for its
// keyword (ScalHFX, ScalDFX, ScalDFC, "Scaling for aa/bb pairs", "Scaling for ab pairs"), read from the output on
// Oct 31 2026; "Psi4" = psi4/driver/procrouting/dft/dh_functionals.py (master, Oct 2026); "s-dftd3" =
// simple-dftd3 assets/parameters.toml (main, Oct 2026). The papers of Kozuch/Martin 2013 and Santra et al. 2019
// were not read (paywalled); where a number rests on programs only, the row says so.
}  // namespace
#define DH_SRC_B2PLYP "Grimme, JCP 124, 034108 (2006), DOI 10.1063/1.2148954; same numbers in ORCA 6.1 and Psi4"
#define DH_SRC_B2GP "Karton et al., JPCA 112, 12868 (2008), DOI 10.1021/jp801805p; same numbers in ORCA 6.1 and Psi4"
#define DH_SRC_D3BJ2011 "D3(BJ): Goerigk, Grimme, PCCP 13, 6670 (2011), DOI 10.1039/c0cp02984j (s-dftd3 table)"
#define DH_SRC_KM2013 "Kozuch, Martin, JCC 34, 2327 (2013), DOI 10.1002/jcc.23391 (paper not read)"

// Out of line (the inline version in qm_citations.h exceeded GCC's inline limit, -Winline). Claude Generated (WP12).
namespace qmcite {
void doubleHybrid(const std::string& name, bool ri, bool d3bj)
{
    if (name.rfind("b2plyp", 0) == 0) CitationRegistry::cite("b2plyp");
    if (name.rfind("b2gp-plyp", 0) == 0) CitationRegistry::cite("b2gpplyp");
    if (name.rfind("dsd-", 0) == 0) {
        CitationRegistry::cite("dsd_blyp_2010");
        CitationRegistry::cite("dsd_2013");
    }
    if (name.rfind("dsd-pbep86", 0) == 0) CitationRegistry::cite("dsd_pbep86_2011");
    if (name.rfind("revdsd-", 0) == 0) {
        CitationRegistry::cite("revdsd_2019");
        CitationRegistry::cite("revdsd_pz81_2023");
    }
    mp2(ri ? "ri-mp2" : "mp2");
    if (d3bj) {
        CitationRegistry::cite("d3");
        CitationRegistry::cite("bj");
        if (name.rfind("b2", 0) == 0) CitationRegistry::cite("goerigk_grimme_2011");
    }
}
}  // namespace qmcite

// ---------------------------------------------------------------------------------------------------------------
// The table. Kernel weights are those of the DFT exchange (1 - a_x) and correlation (c_c); B88_X and PBE_X are the
// full GGA exchange (LDA part included, libxc GGA_X_B88 / GGA_X_PBE), LYP_C / P86_C / B95_C the full correlation.
// ---------------------------------------------------------------------------------------------------------------
const std::vector<DoubleHybridDef>& doubleHybridTable()
{
    static const std::vector<DoubleHybridDef> table = {
        { "b2plyp", "B2PLYP (0.53 HF + 0.47 B88, 0.73 LYP, 0.27 PT2)", 0.53,
          { { XCKernel::B88_X, 0.47 }, { XCKernel::LYP_C, 0.73 } }, 0.27, 0.27, DHDispersion::None, 0, 0, 0, 0, false,
          DH_SRC_B2PLYP },
        { "b2plyp-d3bj", "B2PLYP-D3(BJ)", 0.53, { { XCKernel::B88_X, 0.47 }, { XCKernel::LYP_C, 0.73 } }, 0.27, 0.27,
          DHDispersion::D3BJ, 0.64, 0.3065, 0.9147, 5.0570, false,
          DH_SRC_B2PLYP "; " DH_SRC_D3BJ2011 },
        { "b2gp-plyp", "B2GP-PLYP (0.65 HF + 0.35 B88, 0.64 LYP, 0.36 PT2)", 0.65,
          { { XCKernel::B88_X, 0.35 }, { XCKernel::LYP_C, 0.64 } }, 0.36, 0.36, DHDispersion::None, 0, 0, 0, 0, false,
          DH_SRC_B2GP },
        { "b2gp-plyp-d3bj", "B2GP-PLYP-D3(BJ)", 0.65, { { XCKernel::B88_X, 0.35 }, { XCKernel::LYP_C, 0.64 } }, 0.36,
          0.36, DHDispersion::D3BJ, 0.56, 0.0, 0.2597, 6.3332, false, DH_SRC_B2GP "; " DH_SRC_D3BJ2011 },
        // DSD without dispersion: the "not dispersion corrected" set of KM2013 (= ORCA's DSD-BLYP, DSD-BLYP/2013)
        { "dsd-blyp", "DSD-BLYP (0.75 HF + 0.25 B88, 0.53 LYP, 0.46 OS / 0.60 SS PT2)", 0.75,
          { { XCKernel::B88_X, 0.25 }, { XCKernel::LYP_C, 0.53 } }, 0.46, 0.60, DHDispersion::None, 0, 0, 0, 0, false,
          DH_SRC_KM2013 "; ORCA 6.1 DSD-BLYP and Psi4 DSD-BLYP agree" },
        { "dsd-blyp-d3bj", "DSD-BLYP-D3(BJ) (0.71 HF + 0.29 B88, 0.54 LYP, 0.47 OS / 0.40 SS PT2)", 0.71,
          { { XCKernel::B88_X, 0.29 }, { XCKernel::LYP_C, 0.54 } }, 0.47, 0.40, DHDispersion::D3BJ, 0.57, 0.0, 0.0, 5.4,
          false, DH_SRC_KM2013 "; XC, PT2 and D3(BJ) from Psi4 DSD-BLYP-D3BJ only (one source)" },
        { "dsd-pbep86", "DSD-PBEP86 (0.72 HF + 0.28 PBE, 0.44 P86, 0.51 OS / 0.36 SS PT2)", 0.72,
          { { XCKernel::PBE_X, 0.28 }, { XCKernel::P86_C, 0.44 } }, 0.51, 0.36, DHDispersion::None, 0, 0, 0, 0, true,
          DH_SRC_KM2013 "; ORCA 6.1 DSD-PBEP86 and Psi4 DSD-PBEP86 agree" },
        { "dsd-pbep86-d3bj", "DSD-PBEP86-D3(BJ) (0.69 HF + 0.31 PBE, 0.44 P86, 0.52 OS / 0.22 SS PT2)", 0.69,
          { { XCKernel::PBE_X, 0.31 }, { XCKernel::P86_C, 0.44 } }, 0.52, 0.22, DHDispersion::D3BJ, 0.48, 0.0, 0.0, 5.6,
          true, DH_SRC_KM2013 "; XC/PT2 from Psi4 DSD-PBEP86-D3BJ, D3(BJ) from Psi4 and s-dftd3 (agree)" },
        // dsd-pbeb95 (no dispersion) is NOT in the table: ORCA 6.1 prints 0.70 HF / 0.52 B95 and Psi4 has 0.69 / 0.54
        // (both 0.48 OS / 0.22 SS, both citing KM2013); without the paper the definition is not settled (WP12 doc).
        { "dsd-pbeb95-d3bj", "DSD-PBEB95-D3(BJ) (0.66 HF + 0.34 PBE, 0.55 B95, 0.46 OS / 0.09 SS PT2)", 0.66,
          { { XCKernel::PBE_X, 0.34 }, { XCKernel::B95_C, 0.55 } }, 0.46, 0.09, DHDispersion::D3BJ, 0.61, 0.0, 0.0, 6.2,
          false, DH_SRC_KM2013 "; XC/PT2 from Psi4 DSD-PBEB95-D3BJ, D3(BJ) from Psi4 and s-dftd3 (agree)" },
        // revDSD-PBEP86-D4 with PZ81-P86: Santra, Sylvetsky, Martin, JPCA 123, 5129 (2019), DOI 10.1021/acs.jpca.9b03157
        // (paper not read); numbers from Mehta, Santra, Martin, Can. J. Chem. 101, 656 (2023), DOI 10.1139/cjc-2022-0277
        // (arXiv 2211.03105, Table 1, PZ81 row: c_X 0.69, c_C 0.4224, c2ab 0.5935, c2ss 0.0566, D4 s6 0.5917, a1 0.371,
        // a2 4.201, s8 0); ORCA 6.1 revDSD-PBEP86/2021 prints the same four XC/PT2 numbers (but on VWN5-P86).
        // The dftd4 parameter file lists other D4 numbers for revdsdpbep86 (s6 0.5132, a1 0.44, a2 3.60).
        { "revdsd-pbep86", "revDSD-PBEP86 without D4 (0.69 HF + 0.31 PBE, 0.4224 P86, 0.5935 OS / 0.0566 SS PT2)", 0.69,
          { { XCKernel::PBE_X, 0.31 }, { XCKernel::P86_C, 0.4224 } }, 0.5935, 0.0566, DHDispersion::None, 0, 0, 0, 0,
          true, "Santra et al., JPCA 123, 5129 (2019), DOI 10.1021/acs.jpca.9b03157 (paper not read); numbers from Mehta "
                "et al., Can. J. Chem. 101, 656 (2023), DOI 10.1139/cjc-2022-0277, Table 1 (PZ81); = ORCA 6.1 "
                "revDSD-PBEP86/2021 XC/PT2" },
        { "revdsd-pbep86-d4", "revDSD-PBEP86-D4", 0.69, { { XCKernel::PBE_X, 0.31 }, { XCKernel::P86_C, 0.4224 } },
          0.5935, 0.0566, DHDispersion::D4, 0.5917, 0.371, 0.0, 4.201, true,
          "as revdsd-pbep86; D4 (s6 0.5917, a1 0.371, a2 4.201, s8 0, ATM 1) from Mehta et al. 2023 Table 1, NOT "
          "implemented" },
    };
    return table;
}

const DoubleHybridDef* findDoubleHybrid(const std::string& name)
{
    const std::string n = lower(name);
    for (const DoubleHybridDef& d : doubleHybridTable())
        if (n == d.name) return &d;
    return nullptr;
}

qmxc::XCFunctional doubleHybridXC(const DoubleHybridDef& d, bool p86_vwn5)
{
    qmxc::XCFunctional f;
    f.name = d.name;
    f.exact_exchange = d.a_x;
    f.terms.assign(std::begin(d.terms), std::end(d.terms));
    for (qmxc::XCTerm& t : f.terms)
        if (t.kernel == XCKernel::P86_C && p86_vwn5) t.kernel = XCKernel::P86VWN_C;
    std::string s;
    for (const qmxc::XCTerm& t : f.terms) s += fmt::format(" + {:g} {}", t.weight, qmxc::kernelName(t.kernel));
    f.description = fmt::format("{} KS part ({:g} HF{})", d.label, d.a_x, s);
    return f;
}

DoubleHybridMethod::DoubleHybridMethod(const std::string& name, const json& config)
{
    const DoubleHybridDef* d = findDoubleHybrid(name);
    if (!d) throw std::invalid_argument("double hybrid '" + name + "' is not in the table (double_hybrid_method.cpp)");
    m_def = d;
    const json ec = QMMethod::engineConfig(config);  // qm defaults < top level < dft < scf < qm

    const std::string p86 = readString(ec, "p86_lda", "pz81");
    if (p86 != "pz81" && p86 != "vwn5")
        throw std::invalid_argument("-qm.p86_lda '" + p86 + "' is unknown; use pz81 or vwn5");
    m_p86_vwn5 = (p86 == "vwn5");

    const std::string disp = readString(ec, "dispersion", "auto");
    if (disp == "auto")
        m_dispersion = m_def->dispersion;
    else if (disp == "none")
        m_dispersion = DHDispersion::None;
    else if (disp == "d3bj" || disp == "d3(bj)") {
        if (m_def->dispersion != DHDispersion::D3BJ)
            throw std::invalid_argument(fmt::format("-qm.dispersion d3bj: method {} has no D3(BJ) parameters in the "
                                                    "table (use the -d3bj method name, e.g. {}-d3bj, if it exists)",
                                                    m_def->name, m_def->name));
        m_dispersion = DHDispersion::D3BJ;
    } else if (disp == "d4")
        m_dispersion = DHDispersion::D4;
    else
        throw std::invalid_argument("-qm.dispersion '" + disp + "' is unknown; use auto, none, d3bj or d4");

    const std::string pt2 = readString(ec, "dh_pt2", "ri");
    if (pt2 != "ri" && pt2 != "canonical")
        throw std::invalid_argument("-qm.dh_pt2 '" + pt2 + "' is unknown; use ri or canonical");
    m_options.ri = (pt2 == "ri");
    m_options.frozen_core = readBool(ec, "frozen_core", false);
    if (ec.contains("mp2_max_memory_mb") && ec["mp2_max_memory_mb"].is_number())
        m_options.max_memory_mb = ec["mp2_max_memory_mb"].get<double>();
    if (ec.contains("mp2_batch_mb") && ec["mp2_batch_mb"].is_number())
        m_options.batch_mb = ec["mp2_batch_mb"].get<double>();
    m_options.gpu_mp2 = readString(ec, "gpu_mp2", "auto");
    if (ec.contains("ri_mp2_aux") && ec["ri_mp2_aux"].is_string()) m_options.ri_aux = ec["ri_mp2_aux"].get<std::string>();

    m_engine = std::make_unique<QMEngine>(QMFunctional::DoubleHybrid, ec);
    if (!m_engine->setXCFunctional(doubleHybridXC(*m_def, m_p86_vwn5)))
        throw std::runtime_error("double hybrid: the KS functional could not be set");
    m_engine->setTightOrbitals(true);  // the PT2 energy is first order in the orbital error (as mp2)

    if (m_dispersion == DHDispersion::D3BJ) {
        // Native D3(BJ) of hf-3c (D3ParameterGenerator), two-body only as the D3BJ fits of these functionals (Psi4
        // "d3bj2b"); rational damping R0 = a1 sqrt(C8/C6) + a2 (a1 = 0: a constant a2 in Bohr).
        json d3 = { { "d3_s6", m_def->s6 }, { "d3_s8", m_def->s8 }, { "d3_a1", m_def->a1 }, { "d3_a2", m_def->a2 },
                    { "d3_s9", 0.0 }, { "d3_alp", 14.0 } };
        m_d3 = std::make_unique<D3ParameterGenerator>(ConfigManager("d3param", d3));
    }
    qmcite::doubleHybrid(m_def->name, m_options.ri, m_dispersion == DHDispersion::D3BJ);
}

DoubleHybridMethod::~DoubleHybridMethod() = default;

bool DoubleHybridMethod::setMolecule(const Mol& mol)
{
    m_molecule = mol;
    m_error = false;
    m_engine->resetForNewMolecule();
    if (!m_engine->QMInterface::InitialiseMolecule(mol)) {
        m_error = true;
        m_error_message = std::string(m_def->name) + ": QM engine initialisation failed";
        return false;
    }
    return true;
}

bool DoubleHybridMethod::updateGeometry(const Matrix& geometry)
{
    m_molecule.m_geometry = geometry;
    return m_engine->UpdateMolecule(geometry);
}

double DoubleHybridMethod::calculateEnergy(bool gradient)
{
    m_error = false;
    m_error_message.clear();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (gradient) m_gradient_bohr = Matrix();  // no stale gradient after a failed request
    if (m_dispersion == DHDispersion::D4) {
        m_error = true;
        m_error_message = std::string(m_def->name) + ": D4 dispersion for double hybrids is not available in curcuma (native D4 is "
                                       "GFN2-only); use -qm.dispersion none for the KS + PT2 part (= " + m_def->name
            + " without D4) or a -d3bj double hybrid";
        CurcumaLogger::error(m_error_message);
        return nan;
    }

    // 1. hybrid KS SCF
    m_e_ks = m_engine->Calculation(false);
    if (!m_engine->scfConverged()) {
        m_error = true;
        m_error_message = std::string(m_def->name) + ": the KS SCF did not converge";
        return m_e_ks;
    }
    // 2. PT2 with the KS orbitals and eigenvalues
    if (!(m_options.ri ? m_engine->computeRIMP2(m_options, m_pt2) : m_engine->computeMP2(m_options, m_pt2))) {
        m_error = true;
        m_error_message = std::string(m_def->name) + ": PT2 part failed (see the error above)";
        return nan;
    }
    m_e_pt2 = m_def->c_os * m_pt2.e_os + m_def->c_ss * m_pt2.e_ss;
    // 3. dispersion (WP12 stage 2: with its gradient, pair part + CN chain as hf-3c, Eh/Bohr)
    m_e_disp = 0.0;
    const int nat = m_molecule.AtomCount();
    Matrix g_disp = Matrix::Zero(nat, 3);
    if (m_d3) {
        m_d3->prepareForEnergyGradient(m_molecule.m_atoms, m_molecule.m_geometry);
        Vector dEdcn = Vector::Zero(nat);
        m_e_disp = m_d3->getEnergyAndGradient(gradient, g_disp, dEdcn);
        if (gradient && dEdcn.size() == nat)
            CNCalculator::addD3CNGradient(m_molecule.m_atoms, m_molecule.m_geometry, dEdcn, g_disp, /*k1=*/16.0,
                                          /*k2=*/4.0 / 3.0,
                                          /*distance_unit_to_bohr=*/CurcumaUnit::Length::BOHR_TO_ANGSTROM);
    }
    m_e_total = m_e_ks + m_e_pt2 + m_e_disp;
    // 4. analytic gradient (WP12 stage 2): KS + PT2 (Z-vector with the KS orbital Hessian) + dispersion
    if (gradient) {
        if (!m_engine->computeDoubleHybridGradient(m_options, m_def->c_os, m_def->c_ss, m_pt2_grad)) {
            m_error = true;
            m_error_message = std::string(m_def->name) + ": " + m_engine->gradientError();
            m_gradient_bohr = Matrix::Zero(nat, 3);
            return nan;
        }
        m_gradient_bohr = m_engine->gradientBohr() + g_disp;
        m_gradient_parts = m_engine->gradientParts();
        m_gradient_parts.push_back(g_disp);
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info(fmt::format("{} gradient: |g| = {:.6e} Eh/Bohr; Z-vector {} iterations (max residual "
                                            "{:.1e}); PT2 amplitudes {:.1f} ms, back-transformation {:.1f} ms, Z-vector "
                                            "{:.1f} ms, 1e {:.1f} ms, 2e {:.1f} ms ({})",
                                            m_def->name, m_gradient_bohr.norm(), m_pt2_grad.zvector_iterations,
                                            m_pt2_grad.zvector_residual, m_pt2_grad.t_grad_amplitudes_ms,
                                            m_pt2_grad.t_grad_backtrans_ms, m_pt2_grad.t_grad_zvector_ms,
                                            m_pt2_grad.t_grad_deriv1e_ms, m_pt2_grad.t_grad_deriv2e_ms, m_pt2_grad.info));
    }

    if (CurcumaLogger::get_verbosity() >= 1) {
        const std::string pt2 = m_options.ri
            ? fmt::format("RI-PT2 (aux {}, {} functions; JK set as correlation fit, not def2-SVP/C)", m_pt2.aux_name,
                          m_pt2.n_aux)
            : std::string("canonical PT2 (exact integrals)");
        CurcumaLogger::result(fmt::format("native double hybrid {} -- {}{}; {}; {} frozen core orbital(s), {} active "
                                          "occupied, {} virtual; dispersion {}",
                                          m_def->name, m_def->label,
                                          m_def->uses_p86 ? (m_p86_vwn5 ? ", P86 on VWN5" : ", P86 on PZ81") : "", pt2,
                                          m_pt2.n_frozen, m_pt2.n_occ_active, m_pt2.n_virt, dispersionName(m_dispersion)));
        CurcumaLogger::result(fmt::format("  E(KS hybrid)   = {:.12f} Eh", m_e_ks));
        CurcumaLogger::result(fmt::format("  E(PT2)         = {:.12f} Eh  ({:g} x OS {:.12f} + {:g} x SS {:.12f})", m_e_pt2,
                                          m_def->c_os, m_pt2.e_os, m_def->c_ss, m_pt2.e_ss));
        if (m_dispersion != DHDispersion::None)
            CurcumaLogger::result(fmt::format("  E(disp {})  = {:.12f} Eh  (s6 {:g}, a1 {:g}, s8 {:g}, a2 {:g}, two-body)",
                                              dispersionName(m_dispersion), m_e_disp, m_def->s6, m_def->a1, m_def->s8,
                                              m_def->a2));
        CurcumaLogger::result(fmt::format("  E({})  = {:.12f} Eh", m_def->name, m_e_total));
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::param("E_xc", fmt::format("{:.12f} Eh", m_engine->xcEnergy()));
            CurcumaLogger::param("source", m_def->source);
            CurcumaLogger::info(fmt::format("PT2 timing: RI setup {:.1f} ms, transformation {:.1f} ms, energy {:.1f} ms ({})",
                                            m_pt2.t_ri_setup_ms, m_pt2.t_transform_ms, m_pt2.t_energy_ms, m_pt2.info));
        }
        CurcumaLogger::energy_abs(m_e_total, "total " + std::string(m_def->name) + " energy");
    }
    return m_e_total;
}

// Claude Generated (Oct 2026, WP12 stage 2).
Matrix DoubleHybridMethod::getGradient() const
{
    if (m_gradient_bohr.rows() != m_molecule.AtomCount()) return Matrix::Zero(m_molecule.AtomCount(), 3);
    return m_gradient_bohr * CurcumaUnit::Length::ANGSTROM_TO_BOHR;  // Eh/Bohr -> Eh/Angstrom
}

json DoubleHybridMethod::getEnergyDecomposition() const
{
    json d;
    d["ks"] = m_e_ks;
    d["xc"] = m_engine ? m_engine->xcEnergy() : 0.0;
    d["os"] = m_pt2.e_os;
    d["ss"] = m_pt2.e_ss;
    d["c_os"] = m_def->c_os;
    d["c_ss"] = m_def->c_ss;
    d["pt2"] = m_e_pt2;
    d["disp"] = m_e_disp;
    d["total"] = m_e_total;
    d["n_frozen"] = m_pt2.n_frozen;
    d["n_aux"] = m_pt2.n_aux;
    d["pt2_mode"] = m_options.ri ? "ri" : "canonical";
    d["dispersion"] = dispersionName(m_dispersion);
    d["p86_lda"] = m_p86_vwn5 ? "vwn5" : "pz81";
    d["scf_iterations"] = m_engine ? m_engine->scfIterations() : 0;
    return d;
}
