/*
 * <Native ab-initio QM Method Wrapper Implementation>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated: Native KS-DFT method wrapper for MethodFactory integration,
 *                   renamed DFTMethod -> QMMethod (Sep 2026)
 *
 * This program is free software under GPL-3.0
 */

#include "qm_method.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/method_factory.h"
#include "src/core/parameter_registry.h"
#include "src/core/units.h"
#include "src/core/energy_calculators/dispersion/d4_dft.h"
#include "qm_citations.h"

#include <fmt/format.h>

#include <algorithm>

json QMMethod::engineConfig(const json& config)
{
    json full_config = ParameterRegistry::getInstance().getDefaultJson("qm");
    // The CLI auto-routes every -qm.<param> (and its flat form) into the "qm"
    // module scope, i.e. controller["qm"], which is the ONLY place they land --
    // so the scope has to be merged explicitly. Reading the controller's top level
    // alone (as this wrapper once did) silently dropped every scoped flag.
    // Order: top level (flat/legacy fallback) < "dft" (the pre-Sep-2026 module
    // name, kept so old command lines and -import_config files still work) < "qm".
    //
    // Oct 27 2026 (master's standard for the shared SCF knobs, docs/SCF_MODES.md): the shared module "scf"
    // (scf_parameters.h; -scf_threshold, -scf.scf_threshold, -scf_max_iterations, ...) sits between "dft" and
    // "qm", as for gfn1/gfn2 (applyXtbScfConfig: method scope > scf > default). Keys of the qm and dft scopes
    // are mapped onto their canonical qm names, so the old spellings -qm.diis_start / -qm.diis_subspace keep
    // working; the top level is not alias-resolved (its diis_start / max_iterations are the optimiser's).
    auto& reg = ParameterRegistry::getInstance();
    auto merge_scope = [&](const json& scope) {
        for (const auto& [key, value] : scope.items()) {
            const std::string canon = reg.resolveAlias("qm", key);
            full_config[canon.empty() ? key : canon] = value;
        }
    };
    if (!config.empty())
        full_config.merge_patch(config);
    if (config.contains("dft") && config["dft"].is_object())
        merge_scope(config["dft"]);
    if (config.contains("scf") && config["scf"].is_object()) {
        json shared = config;  // aliases resolved and "method default" markers dropped (idempotent)
        MethodFactory::applySharedScfScope(shared, "hf");
        if (shared.contains("scf") && shared["scf"].is_object())
            for (const auto& [key, value] : shared["scf"].items())
                full_config[key] = value;  // the scf names are the qm names (scf_threshold, scf_mode, ...)
    }
    if (config.contains("qm") && config["qm"].is_object())
        merge_scope(config["qm"]);
    return full_config;
}

QMMethod::QMMethod(QMFunctional functional, const json& config, bool with_d4_name)
    : m_engine(nullptr)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    const json ec = engineConfig(config);
    m_engine = std::make_unique<QMEngine>(functional, ec);
    m_method_name = m_engine->getMethodNameStr();

    // Convert to lowercase for MethodFactory consistency
    std::transform(m_method_name.begin(), m_method_name.end(),
                   m_method_name.begin(), ::tolower);

    // WP12 stage 4 (Claude Generated, Oct 2026): native DFT-D4 for pbe/b3lyp/tpss/scan. -qm.dispersion auto = the
    // method name (-d4 or not), d4 = add it, none = drop it; d3bj has no parameters for these functionals here.
    std::string disp = "auto";
    if (ec.contains("dispersion") && ec["dispersion"].is_string()) disp = ec["dispersion"].get<std::string>();
    std::transform(disp.begin(), disp.end(), disp.begin(), ::tolower);
    bool want_d4 = with_d4_name;
    if (disp == "d4")
        want_d4 = true;
    else if (disp == "none")
        want_d4 = false;
    else if (disp != "auto")
        throw std::invalid_argument("-qm.dispersion '" + disp + "' is not available for " + m_method_name
                                    + " (auto, none or d4; D3(BJ) only for hf-3c and the -d3bj double hybrids)");
    if (want_d4) {
        const curcuma::dispersion::D4DftParams* p = curcuma::dispersion::findD4DftParameters(m_method_name);
        if (!p || functional == QMFunctional::HF || functional == QMFunctional::LDA)
            throw std::invalid_argument("-qm.dispersion d4: no D4 parameters for " + m_method_name
                                        + " (pbe, b3lyp, tpss, scan and the double hybrids b2plyp, b2gp-plyp, "
                                          "revdsd-pbep86)");
        m_d4 = std::make_unique<curcuma::dispersion::D4DftDispersion>(*p);
        m_method_name += "-d4";
        qmcite::d4Dft();
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("QMMethod initialized: " + m_engine->getMethodNameStr());
    }
}

QMMethod::~QMMethod() = default;

bool QMMethod::setMolecule(const Mol& mol)
{
    m_molecule = mol;
    m_calculation_done = false;

    m_engine->resetForNewMolecule();  // a reused object must not keep the old molecule
    return m_engine->QMInterface::InitialiseMolecule(mol);
}

bool QMMethod::updateGeometry(const Matrix& geometry)
{
    m_calculation_done = false;
    m_molecule.m_geometry = geometry;  // WP12 stage 4: the D4 term reads it
    // QMInterface::UpdateMolecule(Matrix) stores the geometry and calls the
    // engine's UpdateMolecule(), which invalidates the integral/SCF caches.
    return m_engine->UpdateMolecule(geometry);
}

double QMMethod::calculateEnergy(bool gradient)
{
    m_has_error = false;
    m_error_message.clear();
    m_last_energy = m_engine->Calculation(gradient);
    m_calculation_done = true;
    if (m_engine->spinStateRejected()) {
        // stage 5b (Nov 2026): electrons/-spin/-qm.reference contradict, or no open-shell variant exists
        m_has_error = true;
        m_error_message = m_engine->spinError();
    } else if (m_engine->riSetupFailed()) {
        // Claude Generated (Oct 2026, WP10): no energy without the fit the user asked for.
        m_has_error = true;
        m_error_message = "density fitting (-qm.ri " + m_engine->riMode() + ") setup failed, see the error above";
    } else if (m_engine->cosxFailed()) {
        // Claude Generated (Nov 2026): refused or failed -qm.exchange_method cosx, no exact-exchange energy under its label
        m_has_error = true;
        m_error_message = m_engine->cosxError();
    } else if (m_engine->xcSetupFailed()) {
        // Claude Generated (Oct 2026, WP5): no LDA/PBE/B3LYP energy without its grid
        m_has_error = true;
        m_error_message = "DFT grid (-qm.grid*) could not be built, see the error above";
    } else if ((m_engine->getFunctional() == QMFunctional::HF || m_engine->isKohnSham()) && !m_engine->scfConverged()
               && !m_engine->scfAllowUnconverged()) {
        // Claude Generated (Oct 27 2026): master's rule for every SCF method (docs/SCF_MODES.md "Failure handling"):
        // an unconverged SCF is an error, not a result. Before, HF printed an error and returned the energy (exit 0).
        m_has_error = true;
        m_error_message = fmt::format("native {} SCF did not converge within {} iterations; "
                                      "-scf_allow_unconverged true accepts the unconverged energy",
                                      m_engine->isKohnSham() ? "Kohn-Sham (" + m_engine->xcFunctional().name + ")"
                                                                : std::string("HF"),
                                      m_engine->scfMaxIter());
    } else if (gradient && m_engine->gradientFailed()) {
        // Claude Generated (Oct 2026): fail loud instead of handing back a zero gradient
        // (RI gradient not possible, e.g. memory limit for the 3-index densities).
        m_has_error = true;
        m_error_message = m_engine->gradientError();
    }
    if (m_d4 && !m_has_error) {
        // WP12 stage 4: native D4 (EEQ charges, ATM), energy Eh, gradient Eh/Bohr incl. dq/dR (Angstrom geometry
        // of m_molecule, kept current by updateGeometry).
        const Matrix geom_bohr = m_molecule.m_geometry * CurcumaUnit::Length::ANGSTROM_TO_BOHR;
        m_d4_gradient_bohr = Matrix();
        if (!m_d4->compute(m_molecule.m_atoms, geom_bohr, static_cast<double>(m_molecule.m_charge), gradient,
                           m_e_disp, m_d4_gradient_bohr)) {
            m_has_error = true;
            m_error_message = m_method_name + ": " + m_d4->error();
            return m_last_energy;
        }
        if (CurcumaLogger::get_verbosity() >= 1)
            CurcumaLogger::result(fmt::format("  E(disp D4)     = {:.12f} Eh  (two-body {:.12f}, ATM {:.12f}; EEQ "
                                              "charges; s6 {:g}, s8 {:g}, a1 {:g}, a2 {:g}, s9 {:g})",
                                              m_e_disp, m_d4->twoBodyEnergy(), m_d4->threeBodyEnergy(),
                                              m_d4->params().s6, m_d4->params().s8, m_d4->params().a1,
                                              m_d4->params().a2, m_d4->params().s9));
        m_last_energy += m_e_disp;
        if (CurcumaLogger::get_verbosity() >= 1)
            CurcumaLogger::energy_abs(m_last_energy, "total " + m_method_name + " energy");
    }
    return m_last_energy;
}

Matrix QMMethod::getGradient() const
{
    if (!m_calculation_done) {
        CurcumaLogger::warn("QMMethod: No calculation done yet");
        return Matrix::Zero(m_molecule.AtomCount(), 3);
    }
    if (!m_engine->hasGradient())
        return Matrix::Zero(m_molecule.AtomCount(), 3);  // no analytic gradient for this functional
    // The engine works in Eh/Bohr; ComputationalMethod::getGradient() is Eh/Angstrom
    // (dE/dx[A] = dE/dx[Bohr] * Bohr-per-Angstrom), see Known Issue #28.
    if (m_d4 && m_d4_gradient_bohr.rows() == m_molecule.AtomCount())
        return (m_engine->gradientBohr() + m_d4_gradient_bohr) * CurcumaUnit::Length::ANGSTROM_TO_BOHR;
    return m_engine->gradientBohr() * CurcumaUnit::Length::ANGSTROM_TO_BOHR;
}

Vector QMMethod::getCharges() const
{
    // No population analysis yet.
    return Vector::Zero(m_molecule.AtomCount());
}

json QMMethod::getEnergyDecomposition() const
{
    json decomp;
    decomp["total"] = m_last_energy;
    if (m_engine->getFunctional() == QMFunctional::HF && m_engine->scfConverged()) {
        // {ET, EV, EJ, Ex, E_elec}; E_nn is the remainder of the total.
        const std::vector<double> c = m_engine->energyComponents();
        decomp["kinetic"] = c[0];
        decomp["nuclear_attraction"] = c[1];
        decomp["coulomb"] = c[2];
        decomp["exchange"] = c[3];
        decomp["electronic"] = c[4];
        decomp["nuclear_repulsion"] = m_last_energy - c[4];
    } else if (m_engine->isKohnSham() && m_engine->scfConverged()) {
        // WP5 (Claude Generated, Oct 2026): {ET, EV, EJ, E_x(HF part, 0 for LDA/PBE), E_elec} + E_xc
        const std::vector<double> c = m_engine->energyComponents();
        decomp["kinetic"] = c[0];
        decomp["nuclear_attraction"] = c[1];
        decomp["coulomb"] = c[2];
        decomp["exchange_correlation"] = m_engine->xcEnergy();
        decomp["xc_exchange"] = m_engine->xcExchangeEnergy();
        decomp["xc_correlation"] = m_engine->xcCorrelationEnergy();
        if (m_engine->xcFunctional().exact_exchange != 0.0)
            decomp["exact_exchange"] = c[3];  // WP7 hybrid: -(a_x/4) Tr(PK)
        decomp["electronic"] = c[4];
        decomp["nuclear_repulsion"] = m_last_energy - c[4];
    } else {
        // unconverged or failed SCF (WP7: no functional returns E_nn alone any more): no components
    }
    if (m_engine->isOpenShell() && m_engine->scfConverged()) {
        // stage 5b: unrestricted reference, <S^2> and its exact value
        decomp["reference"] = m_engine->isKohnSham() ? "UKS" : "UHF";
        decomp["n_alpha"] = m_engine->nAlpha();
        decomp["n_beta"] = m_engine->nBeta();
        decomp["s2"] = m_engine->spinSquared();
        decomp["s2_exact"] = m_engine->spinSquaredExact();
    }
    if (!m_engine->stabilityReport().empty()) decomp["stability"] = m_engine->stabilityReport();  // -qm.scf_stability
    if (m_d4) {
        decomp["dispersion"] = m_e_disp;  // WP12 stage 4; nuclear_repulsion above then holds E_nn + E_disp
        if (decomp.contains("nuclear_repulsion")) decomp["nuclear_repulsion"] = decomp["nuclear_repulsion"].get<double>() - m_e_disp;
    }
    return decomp;
}
