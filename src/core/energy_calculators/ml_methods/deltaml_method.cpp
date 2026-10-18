/*
 * < Delta-ML corrected computational method >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Claude Generated (July 2026): method wrapper of the Delta-ML work package.
 */

#include "deltaml_method.h"

#include "src/core/citation_registry.h"
#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"
#include "src/core/energy_calculators/method_factory.h"
#include "src/core/units.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>

namespace {

/*! \brief Read a string setting from the deltaml block, a flat key, or a default. */
std::string readString(const json& config, const std::string& key, const std::string& fallback)
{
    if (config.contains("deltaml") && config["deltaml"].is_object()
        && config["deltaml"].contains(key) && config["deltaml"][key].is_string())
        return config["deltaml"][key].get<std::string>();
    const std::string flat = "deltaml_" + key;
    if (config.contains(flat) && config[flat].is_string())
        return config[flat].get<std::string>();
    return fallback;
}

bool readBool(const json& config, const std::string& key, bool fallback)
{
    if (config.contains("deltaml") && config["deltaml"].is_object() && config["deltaml"].contains(key))
        return config["deltaml"][key].get<bool>();
    const std::string flat = "deltaml_" + key;
    if (config.contains(flat))
        return config[flat].get<bool>();
    return fallback;
}

double readDouble(const json& config, const std::string& key, double fallback)
{
    if (config.contains("deltaml") && config["deltaml"].is_object() && config["deltaml"].contains(key))
        return config["deltaml"][key].get<double>();
    const std::string flat = "deltaml_" + key;
    if (config.contains(flat))
        return config[flat].get<double>();
    return fallback;
}

} // namespace

std::string DeltaMLMethod::modelPathFromMethodString(const std::string& method)
{
    const std::size_t colon = method.find(':');
    if (colon == std::string::npos)
        return std::string();
    return method.substr(colon + 1);
}

DeltaMLMethod::DeltaMLMethod(const json& config)
    : m_config(config)
{
    m_parameters = config;
    m_user_verbosity = CurcumaLogger::get_verbosity(); // see warnEnergyRoleOnce()

    // The method string may carry the model path. That form exists because several
    // capabilities rebuild their sub-configuration from a handful of keys and would
    // otherwise drop the deltaml block on the way down.
    const std::string method = config.value("method", std::string("deltaml"));
    m_model_path = modelPathFromMethodString(method);
    if (m_model_path.empty())
        m_model_path = readString(config, "model", "");

    m_enabled = readBool(config, "enabled", true);
    m_report_uncertainty = readBool(config, "report_uncertainty", false);
    m_uncertainty_warn = readDouble(config, "uncertainty_warn", 0.9);
    m_uncertainty_abort = readDouble(config, "uncertainty_abort", 0.0);
    m_uncertainty_record = readString(config, "uncertainty_record", "");

    if (m_model_path.empty()) {
        m_has_error = true;
        m_error_message = "deltaml requires a model file: pass -deltaml.model <file.json> "
                          "or use -method deltaml:<file.json>";
        CurcumaLogger::error(m_error_message);
        return;
    }

    if (m_uncertainty_record == "auto") {
        std::string base = m_model_path;
        const std::size_t dot = base.rfind(".json");
        if (dot != std::string::npos)
            base = base.substr(0, dot);
        m_uncertainty_record = base + ".extrapolation.xyz";
    }

    if (!m_model.load(m_model_path)) {
        m_has_error = true;
        m_error_message = "deltaml could not load model '" + m_model_path + "'";
        CurcumaLogger::error(m_error_message);
        return;
    }

    // A model is tied to the method it was trained against; silently pairing it with a
    // different base would produce plausible but wrong energies.
    const std::string override_base = readString(config, "base", "");
    m_base_method = override_base.empty() ? m_model.baseMethod() : override_base;
    if (m_base_method.empty())
        m_base_method = "gfnff";
    if (!override_base.empty() && !m_model.baseMethod().empty()
        && override_base != m_model.baseMethod()) {
        CurcumaLogger::warn_fmt(
            "deltaml: model was trained against '{}' but the base method is overridden to "
            "'{}'. The correction is only valid for the method it was fitted to.",
            m_model.baseMethod(), override_base);
    }
    if (m_base_method.rfind("deltaml", 0) == 0) {
        m_has_error = true;
        m_error_message = "deltaml cannot use another deltaml model as its base method";
        CurcumaLogger::error(m_error_message);
        return;
    }

    // The global -gpu setting drives the correction as well, not just the base method.
    const std::string gpu = config.value("gpu", std::string());
    if (!gpu.empty())
        m_model.setInferenceBackend(gpu);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("deltaml_model", m_model_path);
        CurcumaLogger::param("deltaml_base", m_base_method);
        CurcumaLogger::param("deltaml_reference", m_model.referenceMethod());
        CurcumaLogger::param("deltaml_centers", std::to_string(m_model.centerCount()));
        CurcumaLogger::param("deltaml_role", m_model.role());
    }
}

bool DeltaMLMethod::initialiseBase(const Mol& mol)
{
    json base_config = m_config;
    base_config["method"] = m_base_method;
    base_config.erase("deltaml");

    try {
        m_base = MethodFactory::create(m_base_method, base_config);
    } catch (const std::exception& e) {
        m_has_error = true;
        m_error_message = std::string("deltaml could not create base method '") + m_base_method
            + "': " + e.what();
        CurcumaLogger::error(m_error_message);
        return false;
    }
    if (!m_base) {
        m_has_error = true;
        m_error_message = "deltaml could not create base method '" + m_base_method + "'";
        CurcumaLogger::error(m_error_message);
        return false;
    }

    // Topology anchor: hand the base method the bond graph captured when the model was
    // fitted, and let it set itself up at the structure actually being evaluated. GFN-FF
    // takes Mol::m_bonds as an external graph and skips its own geometric detection, so
    // the connectivity is pinned while coordination numbers, the parameters derived from
    // them and the charges are all computed where they belong.
    //
    // The earlier construction -- initialise at the anchor geometry, then move -- pinned
    // far more than the graph: the force-field parameters stayed at the anchor while the
    // topology was rebuilt at the real geometry, and the mixture depended on how far the
    // structure had travelled. The result was a potential whose analytic gradient missed
    // its own finite differences by 4.7e-3 Eh/Angstrom, against 3.5e-5 for the same method
    // set up normally. Optimisation and dynamics on such a surface are meaningless.
    Mol anchored = mol;
    if (!m_model.anchorBonds().empty() && m_model.anchorElements() == mol.m_atoms)
        anchored.m_bonds = m_model.anchorBonds();
    else if (m_model.hasAnchor() && m_model.anchorElements() == mol.m_atoms)
        CurcumaLogger::warn("deltaml: model carries no bond graph (fitted before the anchor "
                            "was narrowed); topology is redetected per structure");

    if (!m_base->setMolecule(anchored)) {
        m_has_error = true;
        m_error_message = "deltaml: base method rejected the molecule";
        return false;
    }

    return true;
}

bool DeltaMLMethod::setMolecule(const Mol& mol)
{
    if (m_has_error && !m_base)
        return false;

    m_elements = mol.m_atoms;
    m_geometry = mol.m_geometry;
    m_gradient = Matrix::Zero(mol.m_number_atoms, 3);

    if (!initialiseBase(mol))
        return false;

    m_initialized = true;
    return true;
}

bool DeltaMLMethod::updateGeometry(const Matrix& geometry)
{
    m_geometry = geometry;
    if (!m_base)
        return false;
    return m_base->updateGeometry(geometry);
}

double DeltaMLMethod::calculateEnergy(bool gradient)
{
    // The correction rests on three distinct pieces of published work; the base method
    // registers its own citations separately.
    CitationRegistry::cite("deltaml");
    CitationRegistry::cite("acsf", "deltaml");
    CitationRegistry::cite("gap_sparse", "deltaml");

    if (!m_base) {
        m_has_error = true;
        m_error_message = "deltaml has no base method";
        return 0.0;
    }

    m_last_base_energy = m_base->calculateEnergy(gradient);
    m_last_correction = 0.0;
    m_last_uncertainty = -1.0;

    if (gradient) {
        m_gradient = m_base->getGradient();
        warnEnergyRoleOnce();
    } else {
        m_gradient = Matrix::Zero(static_cast<int>(m_elements.size()), 3);
    }

    if (m_enabled && m_model.isTrained()) {
        // The indicator falls out of the prediction, so it is always available at no cost.
        double indicator = -1.0;
        if (gradient) {
            Matrix dgrad;
            m_last_correction
                = m_model.predictEnergyAndGradient(m_elements, m_geometry, dgrad, &indicator);
            // Model and base method both report Eh/Angstrom: getGradient() has that contract
            // since Known Issue #28. An earlier version of this file assumed the base returned
            // Eh/Bohr and scaled the learned force by BOHR_TO_ANGSTROM; after the merge with
            // master that made the correction 0.529 times too small (found by comparing the
            // assembled gradient with central differences of the -sp energy). An error like
            // this is invisible in every energy, so keep the finite-difference gate.
            m_gradient += dgrad;
        } else {
            m_last_correction = m_model.predictEnergy(m_elements, m_geometry, &indicator);
        }
        m_last_uncertainty = indicator;

        if (indicator > m_uncertainty_warn) {
            if (!m_uncertainty_record.empty())
                recordExtrapolation(indicator);

            if (!m_warned) {
                m_warned = true;
                // Quantitative rather than alarming: state how far out the structure is and
                // by how much the correction is attenuated there, so the reader can judge
                // whether it matters instead of taking the word "outside" on trust.
                const double d_over_sigma = (indicator >= 1.0)
                    ? std::numeric_limits<double>::infinity()
                    : std::sqrt(std::max(0.0, -2.0 * std::log(1.0 - indicator)));
                CurcumaLogger::warn_fmt(
                    "deltaml: extrapolation indicator {:.3f} above {:.3f} -- nearest kernel "
                    "centre {:.1f} sigma away, correction attenuated to {:.0f}% there. The "
                    "correction decays smoothly to zero further out, so results degrade "
                    "towards the bare base method rather than diverging.{}",
                    indicator, m_uncertainty_warn, d_over_sigma, 100.0 * (1.0 - indicator),
                    m_uncertainty_record.empty()
                        ? " Set -deltaml.uncertainty_record auto to collect these structures."
                        : "");
            }
        }

        if (m_uncertainty_abort > 0.0 && indicator > m_uncertainty_abort) {
            m_has_error = true;
            m_error_message = fmt::format(
                "deltaml: extrapolation indicator {:.3f} exceeds the abort threshold {:.3f}; "
                "refusing to return an unsupported energy",
                indicator, m_uncertainty_abort);
            CurcumaLogger::error(m_error_message);
        }
    }

    m_last_energy = m_last_base_energy + m_last_correction;

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("deltaml_base_energy", fmt::format("{:.10f} Eh", m_last_base_energy));
        CurcumaLogger::param("deltaml_correction", fmt::format("{:.10f} Eh", m_last_correction));
    }

    return m_last_energy;
}

void DeltaMLMethod::warnEnergyRoleOnce() const
{
    // Claude Generated (October 2026). An energy-only model (role "energy", fitted with
    // use_forces = false) has a smaller relative energy error than the joint model but its
    // forces are measured no better than those of the base method (held-out gradient RMSE 18e-3 against
    // 18e-3 Eh/A on a 90-atom conformer set). Optimisation or dynamics on it would silently
    // give up the force gain of the joint model, so say so -- once per process, since
    // optimisers, MD and ConfSearch create many instances. No abort: the gradient is exact for
    // the surface the model defines, it is only not closer to the reference.
    //
    // Optimiser and MD loops silence the global logger while they evaluate (OptimizerDriver
    // drops it to 0 below verbosity 2), so a plain warn() there is swallowed. The level seen
    // when this method was constructed is the caller's own; if that one admits warnings, the
    // single line is printed through a thread-local override. A caller that was already
    // silent at construction (batch workers, -verbosity 0) stays silent.
    static std::atomic<bool> printed { false };
    if (m_model.role() != "energy" || printed.load())
        return;
    if (std::max(CurcumaLogger::get_verbosity(), m_user_verbosity) < 1)
        return;
    if (printed.exchange(true))
        return;
    const int saved_thread_level = CurcumaLogger::thread_verbosity();
    if (CurcumaLogger::get_verbosity() < 1)
        CurcumaLogger::set_thread_verbosity(1);
    CurcumaLogger::warn_fmt(
        "deltaml: model '{}' is an energy-only model (role 'energy'). Its forces are not "
        "more accurate than those of the base method '{}'. For -opt/-md use the joint model "
        "(<basename>.deltaml.json); this model is meant for relative energies (-sp, ConfScan "
        "energies).",
        m_model_path, m_base_method);
    CurcumaLogger::set_thread_verbosity(saved_thread_level);
}

void DeltaMLMethod::recordExtrapolation(double indicator) const
{
    // Several dynamics threads can share one record file, so serialise the append. The
    // file is opened per call rather than held open: a run that crashes mid-trajectory
    // should still leave the structures that led up to it on disk.
    static std::mutex record_mutex;
    const std::lock_guard<std::mutex> lock(record_mutex);

    std::ofstream out(m_uncertainty_record, std::ios::app);
    if (!out.is_open())
        return;

    const int natoms = static_cast<int>(m_elements.size());
    out << natoms << "\n";
    out << std::fixed << std::setprecision(6) << "uncertainty= " << indicator
        << " E_base= " << std::setprecision(8) << m_last_base_energy
        << " correction= " << m_last_correction << "\n";
    for (int i = 0; i < natoms; ++i) {
        const int z = m_elements[i];
        const std::string sym = (z > 0 && z < static_cast<int>(Elements::ElementAbbr.size()))
            ? Elements::ElementAbbr[z]
            : "X";
        out << std::left << std::setw(3) << sym << std::right << std::fixed << std::setprecision(8)
            << std::setw(16) << m_geometry(i, 0) << std::setw(16) << m_geometry(i, 1)
            << std::setw(16) << m_geometry(i, 2) << "\n";
    }
    ++m_recorded;
}

Vector DeltaMLMethod::getCharges() const
{
    return m_base ? m_base->getCharges() : Vector();
}

Vector DeltaMLMethod::getBondOrders() const
{
    return m_base ? m_base->getBondOrders() : Vector();
}

Position DeltaMLMethod::getDipole() const
{
    return m_base ? m_base->getDipole() : Position(0, 0, 0);
}

Vector DeltaMLMethod::getCN() const
{
    return m_base ? m_base->getCN() : Vector();
}

bool DeltaMLMethod::hasGradient() const
{
    // The learned part is always analytic; the base method decides.
    return m_base ? m_base->hasGradient() : false;
}

std::string DeltaMLMethod::getMethodName() const
{
    return "deltaml(" + m_base_method + ")";
}

bool DeltaMLMethod::isThreadSafe() const
{
    // Prediction only reads the model, so thread safety is entirely the base method's.
    return m_base ? m_base->isThreadSafe() : false;
}

void DeltaMLMethod::setThreadCount(int threads)
{
    m_thread_count = threads;
    if (m_base)
        m_base->setThreadCount(threads);
}

void DeltaMLMethod::setParameters(const json& params)
{
    m_parameters = params;
    if (m_base)
        m_base->setParameters(params);
}

bool DeltaMLMethod::hasError() const
{
    return m_has_error || (m_base && m_base->hasError());
}

std::string DeltaMLMethod::getErrorMessage() const
{
    if (m_has_error)
        return m_error_message;
    return m_base ? m_base->getErrorMessage() : std::string();
}

void DeltaMLMethod::clearError()
{
    m_has_error = false;
    m_error_message.clear();
    if (m_base)
        m_base->clearError();
}

json DeltaMLMethod::getEnergyDecomposition() const
{
    json j = m_base ? m_base->getEnergyDecomposition() : json::object();
    j["base_total"] = m_last_base_energy;
    j["ml_correction"] = m_last_correction;
    j["total"] = m_last_energy;
    if (m_last_uncertainty >= 0.0)
        j["ml_uncertainty"] = m_last_uncertainty;
    return j;
}

void DeltaMLMethod::setWarmStart(bool on)
{
    if (m_base)
        m_base->setWarmStart(on);
}

void DeltaMLMethod::setIterativeMode(bool on)
{
    if (m_base)
        m_base->setIterativeMode(on);
}
