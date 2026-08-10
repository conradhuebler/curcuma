/*
 * <QMDFF Hessian fit for parametrisation. >
 * Copyright (C) 2023 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * =============================================================================
 * Claude Generated (August 2026): rewritten around the linear Hessian fit.
 *
 * The QMDFF energy is exactly linear in every constant this determines (bond fc,
 * angle fc, torsion and out-of-plane scale), so the parametrisation is a single
 * least-squares solve — see qmdff_parametrisation.h. The previous implementation
 * ran Levenberg-Marquardt with a numerically differentiated FULL force-field
 * Hessian per trial vector, which was both unaffordable (~1e6 force-field
 * gradient calls per step at 50 atoms) and a guaranteed no-op, because the trial
 * constants never reached ForceField at all.
 * =============================================================================
 */

#include "src/capabilities/hessian.h"

#include "src/core/config_manager.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/forcefield.h"
#include "src/core/energy_calculators/ff_methods/forcefieldgenerator.h"
#include "src/core/energy_calculators/ff_methods/qmdff_par.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/energy_calculators/qm_methods/gfnff_method.h"
#include "src/core/energy_calculators/ff_methods/qmdff_parametrisation.h"
#include "src/core/energy_calculators/ff_methods/uff_par.h"
#include "src/core/energycalculator.h"
#include "src/core/parameter_registry.h"
#include "src/core/topology.h"

#include <algorithm>
#include <filesystem>
#include <sstream>
#include <fstream>

#include "qmdfffit.h"
#include "src/core/citation_registry.h"

using curcuma::qmdff::FitOptions;
using curcuma::qmdff::FitReport;
using curcuma::qmdff::QMDFFParametrisation;

QMDFFFit::QMDFFFit(const json& controller, bool silent)
    : QMDFFFit(ConfigManager("qmdfffit", controller), silent)
{
}

QMDFFFit::QMDFFFit(const ConfigManager& config, bool silent)
    : CurcumaMethod(json{}, config.exportConfig(), silent)
{
    UpdateController(config.exportConfig());
}

void QMDFFFit::LoadControlJson()
{
    m_method = Json2KeyWord<std::string>(m_defaults, "method");
    m_threads = Json2KeyWord<int>(m_defaults, "threads");
    // Claude Generated (Aug 2026): the registry stores the CANONICAL names
    // (hessian_file / charges_file); "hessian" and "charges" are only CLI aliases.
    // Json2KeyWord does not resolve aliases, so looking them up threw and made
    // every `-qmdfffit` run abort before it started.
    m_hessian_file = Json2KeyWord<std::string>(m_defaults, "hessian_file");
    m_scf_file = Json2KeyWord<std::string>(m_defaults, "charges_file");
    m_potential = Json2KeyWord<std::string>(m_defaults, "potential");

    m_hessian_weight = Json2KeyWord<std::string>(m_defaults, "hessian_weight");
    m_lambda = Json2KeyWord<double>(m_defaults, "lambda");
    m_fd_step = Json2KeyWord<double>(m_defaults, "fd_step");
    m_lambda_torsion = Json2KeyWord<double>(m_defaults, "lambda_torsion");
    m_fit_torsions = Json2KeyWord<bool>(m_defaults, "fit_torsions");
    m_fit_oop = Json2KeyWord<bool>(m_defaults, "fit_oop");
    m_tie_torsions = Json2KeyWord<bool>(m_defaults, "tie_torsions");
    m_nonnegative = Json2KeyWord<bool>(m_defaults, "nonnegative");
    m_verify = Json2KeyWord<bool>(m_defaults, "verify");
    m_verify_fd_step = Json2KeyWord<double>(m_defaults, "verify_fd_step");
    m_write_param = Json2KeyWord<bool>(m_defaults, "write_param");
    m_basins = Json2KeyWord<int>(m_defaults, "basins");
    m_basin_frames = Json2KeyWord<std::string>(m_defaults, "basin_frames");
    m_weight_energy = Json2KeyWord<double>(m_defaults, "weight_energy");
    m_weight_gradient = Json2KeyWord<double>(m_defaults, "weight_gradient");
}

FitOptions QMDFFFit::fitOptions() const
{
    FitOptions options;
    options.fd_step_bohr = m_fd_step;
    options.mass_weighted = (m_hessian_weight == "mass");
    options.lambda = m_lambda;
    options.lambda_torsion = m_lambda_torsion;
    options.fit_torsions = m_fit_torsions;
    options.fit_oop = m_fit_oop;
    options.tie_torsions_by_central_bond = m_tie_torsions;
    options.nonnegative = m_nonnegative;
    options.weight_energy = m_weight_energy;
    options.weight_gradient = m_weight_gradient;
    options.verbosity = CurcumaLogger::get_verbosity();
    return options;
}

// =============================================================================
// Step 1: reference charges
// =============================================================================

bool QMDFFFit::prepareCharges(double& e0, Vector& charges, double& gmax)
{
    e0 = 0.0;
    gmax = 0.0;

    if (std::filesystem::exists(m_scf_file)) {
        std::ifstream scffile(m_scf_file);
        json scfjson;
        scffile >> scfjson;
        if (scfjson.contains("e0"))
            e0 = scfjson["e0"];
        if (scfjson.contains("gradient_max"))
            gmax = scfjson["gradient_max"];
        if (scfjson.contains("charges"))
            charges = Tools::String2EigenVector(scfjson["charges"], "|");
        if (charges.size() == m_molecule.AtomCount()) {
            CurcumaLogger::info(fmt::format("Reusing charges from {}", m_scf_file));
            m_molecule.setPartialCharges(charges);
            return true;
        }
        CurcumaLogger::warn(fmt::format("{} does not contain usable charges — recomputing", m_scf_file));
    }

    EnergyCalculator energy(m_method, m_defaults);
    energy.setMolecule(m_molecule.getMolInfo());
    e0 = energy.CalculateEnergy(true);
    charges = energy.Charges();
    const Matrix gradient = energy.Gradient();
    if (gradient.size() > 0)
        gmax = gradient.cwiseAbs().maxCoeff();

    if (charges.size() != m_molecule.AtomCount()) {
        CurcumaLogger::error(fmt::format("{} returned {} charges for {} atoms — cannot continue",
            m_method, charges.size(), m_molecule.AtomCount()));
        return false;
    }

    // Claude Generated (Aug 2026): the previous code declared a local
    // `std::string charges = Tools::DoubleVector2String(...)` here, which SHADOWED the
    // outer Vector. On a first run the molecule therefore ended up with no partial
    // charges at all, silently disabling the QMDFF electrostatics and HB/XB terms.
    m_molecule.setPartialCharges(charges);

    json scfjson;
    scfjson["e0"] = e0;
    scfjson["charges"] = Tools::DoubleVector2String(charges);
    scfjson["gradient_max"] = gmax;
    scfjson["method"] = m_method;
    std::ofstream scffile(outputPath("scf.json"));
    scffile << scfjson;
    return true;
}

// =============================================================================
// Step 2: reference Hessian (RAW Cartesian, Hartree/Angstrom^2)
// =============================================================================

bool QMDFFFit::prepareHessian(Matrix& hessian, Vector& qm_frequencies)
{
    if (std::filesystem::exists(m_hessian_file)) {
        std::ifstream hess_file(m_hessian_file);
        json hessjson;
        hess_file >> hessjson;

        if (hessjson.value("mass_weighted", false)) {
            CurcumaLogger::error(fmt::format(
                "{} holds a MASS-WEIGHTED Hessian; the fit needs the raw Cartesian one. "
                "Delete the file to recompute it.", m_hessian_file));
            return false;
        }
        const int atoms = hessjson["atoms"];
        if (atoms != m_molecule.AtomCount()) {
            CurcumaLogger::error(fmt::format("{} is for {} atoms, the molecule has {}",
                m_hessian_file, atoms, m_molecule.AtomCount()));
            return false;
        }
        hessian = Matrix::Zero(3 * atoms, 3 * atoms);
        const std::vector<double> flat = Tools::String2DoubleVec(hessjson["hessian"], "|");
        int index = 0;
        for (int i = 0; i < 3 * atoms; ++i)
            for (int j = 0; j < 3 * atoms; ++j)
                hessian(i, j) = flat[index++];
        CurcumaLogger::info(fmt::format("Reusing Hessian from {}", m_hessian_file));
        return true;
    }

    Hessian hessian_calc(m_method, m_defaults, false);
    hessian_calc.setMolecule(m_molecule);
    hessian_calc.start();

    // getRawHessian(), NOT getHessian(): start() mass-weights m_hessian in place, and a
    // mass-weighted matrix would silently produce force constants scaled by sqrt(m_i m_j).
    hessian = hessian_calc.getRawHessian();
    qm_frequencies = hessian_calc.Frequencies();

    json hjson;
    hjson["atoms"] = m_molecule.AtomCount();
    hjson["units"] = "Hartree/Angstrom^2";
    hjson["mass_weighted"] = false;
    hjson["method"] = m_method;
    hjson["hessian"] = Tools::Matrix2String(hessian);
    std::ofstream hess_out(outputPath("hessian.json"));
    hess_out << hjson;
    return true;
}

// =============================================================================
// Step 5: verification
// =============================================================================

void QMDFFFit::verifyFit(const json& fitted, const FitReport& report, const Vector& qm_frequencies)
{
    json ff_controller = m_controller;
    ff_controller["method"] = "qmdff";
    ff_controller["threads"] = m_threads;
    // The auto-parameter cache must not substitute a stale <basename>.param.json for the
    // set we just fitted (ForceField::setParameter consults it first and returns on a hit).
    ff_controller["parameter_caching"] = false;
    // The default 0.005 Angstrom step leaves an O(step^2) truncation error of ~3e-5 relative
    // on the strongly anharmonic QMDFF bond term, which would swamp the deviation we are
    // trying to measure. A force-field gradient is cheap, so take a much tighter step.
    ff_controller["finite_diff_step"] = m_verify_fd_step;
    ff_controller["hessian"]["finite_diff_step"] = m_verify_fd_step;

    // Scoped: CurcumaMethod RAII-clamps the global CurcumaLogger verbosity for as long as
    // the Hessian object lives, so everything reported below must happen after it dies.
    Matrix raw_hessian;
    Vector ff_frequencies;
    {
        Hessian check("qmdff", ff_controller, true);
        check.setMolecule(m_molecule);
        check.setParameter(fitted);
        check.start();
        raw_hessian = check.getRawHessian();
        ff_frequencies = check.Frequencies();
    }

    const Matrix h_raw = QMDFFParametrisation::hessianAngstromToBohr(raw_hessian);
    // The finite-difference Hessian is only nearly symmetric (dg_j/dx_i and dg_i/dx_j differ
    // at O(step^2)), while the model is symmetric by construction. Compare like with like
    // and report the antisymmetric part separately as the finite-difference noise floor.
    const Matrix h_ff = 0.5 * (h_raw + h_raw.transpose());
    const double asym = (0.5 * (h_raw - h_raw.transpose())).norm();

    const double denom = std::max(1e-30, report.hessian_predicted.norm());
    const double rel = (h_ff - report.hessian_predicted).norm() / denom;

    CurcumaLogger::param("qmdff_fit_verification",
        fmt::format("||H_FF - H_predicted|| / ||H_predicted|| = {:.3e} "
                    "(finite-difference asymmetry {:.3e})",
            rel, asym / denom));

    if (rel > 1e-5) {
        CurcumaLogger::error(fmt::format(
            "Verification FAILED (relative deviation {:.3e}). Either the fitted parameters "
            "did not reach ForceField, or the unit-Hessian assembly is inconsistent with the "
            "force-field evaluation. The written parameter set is NOT trustworthy.", rel));
    } else {
        CurcumaLogger::success("Fit verified against an independently computed force-field Hessian");
    }

    if (qm_frequencies.size() == ff_frequencies.size() && qm_frequencies.size() > 6) {
        // Frequencies() returns eigenvalues of the mass-weighted Hessian, sorted ascending;
        // skip the six translation/rotation modes.
        double sum = 0.0, sum_low = 0.0;
        int n = 0, n_low = 0, n_imag = 0;
        for (int i = 6; i < qm_frequencies.size(); ++i) {
            sum += std::abs(qm_frequencies(i) - ff_frequencies(i));
            ++n;
            if (std::abs(qm_frequencies(i)) < 1.0) { // low-frequency part of the spectrum
                sum_low += std::abs(qm_frequencies(i) - ff_frequencies(i));
                ++n_low;
            }
            if (ff_frequencies(i) < 0.0)
                ++n_imag;
        }
        CurcumaLogger::param("qmdff_fit_freq_mad",
            fmt::format("all modes {:.5f}, low modes {:.5f} (mass-weighted eigenvalues)",
                n > 0 ? sum / n : 0.0, n_low > 0 ? sum_low / n_low : 0.0));
        if (n_imag > 0)
            CurcumaLogger::warn(fmt::format("{} imaginary force-field modes at the reference "
                                            "geometry — consider raising -lambda", n_imag));
    }
}



// =============================================================================
// GFN-FF parametrisation support
// =============================================================================
//
// GFN-FF hands native GFNFFParameterSet structs straight to ForceField and never emits
// the JSON term lists the fit engine reads, so we convert here. Only the keys the engine
// needs are written; everything else stays in the struct and is copied back afterwards.

namespace {

json gfnffParametersToJson(const GFNFFParameterSet& ps)
{
    json out;
    json bonds = json::array();
    for (const auto& b : ps.bonds) {
        json j;
        j["i"] = b.i; j["j"] = b.j;
        j["fc"] = b.fc;
        j["exponent"] = b.exponent;
        j["r0_ij"] = b.r0_ij;
        // dynamic r0 inputs — the fit must reproduce the r0 the force field actually uses
        j["z_i"] = b.z_i; j["z_j"] = b.z_j;
        j["r0_base_i"] = b.r0_base_i; j["r0_base_j"] = b.r0_base_j;
        j["cnfak_i"] = b.cnfak_i; j["cnfak_j"] = b.cnfak_j;
        j["rabshift"] = b.rabshift; j["ff"] = b.ff;
        j["nr_hb"] = b.nr_hb; j["hb_cn_H"] = b.hb_cn_H;
        bonds.push_back(j);
    }
    out["bonds"] = bonds;

    json angles = json::array();
    for (const auto& a : ps.angles) {
        json j;
        j["i"] = a.i; j["j"] = a.j; j["k"] = a.k;
        j["fc"] = a.fc;
        j["theta0_ijk"] = a.theta0_ijk;
        angles.push_back(j);
    }
    out["angles"] = angles;

    auto torsionList = [](const std::vector<Dihedral>& list) {
        json arr = json::array();
        for (const auto& d : list) {
            json j;
            j["i"] = d.i; j["j"] = d.j; j["k"] = d.k; j["l"] = d.l;
            j["V"] = d.V; j["n"] = d.n; j["phi0"] = d.phi0;
            j["is_nci"] = d.is_nci;
            arr.push_back(j);
        }
        return arr;
    };
    out["dihedrals"] = torsionList(ps.dihedrals);
    out["extra_dihedrals"] = torsionList(ps.extra_dihedrals);

    json inversions = json::array();
    for (const auto& inv : ps.inversions) {
        json j;
        j["i"] = inv.i; j["j"] = inv.j; j["k"] = inv.k; j["l"] = inv.l;
        j["fc"] = inv.fc;
        j["potential_type"] = inv.potential_type;
        j["omega0"] = inv.omega0;
        inversions.push_back(j);
    }
    out["inversions"] = inversions;
    out["method"] = "gfnff";
    return out;
}

/// Copy the fitted constants back onto the native structs; nothing else is touched.
void applyFittedToGFNFF(const json& fitted, GFNFFParameterSet& ps)
{
    const json bonds = fitted.value("bonds", json::array());
    for (std::size_t i = 0; i < ps.bonds.size() && i < bonds.size(); ++i)
        ps.bonds[i].fc = bonds[i].value("fc", ps.bonds[i].fc);

    const json angles = fitted.value("angles", json::array());
    for (std::size_t i = 0; i < ps.angles.size() && i < angles.size(); ++i)
        ps.angles[i].fc = angles[i].value("fc", ps.angles[i].fc);

    const json dih = fitted.value("dihedrals", json::array());
    for (std::size_t i = 0; i < ps.dihedrals.size() && i < dih.size(); ++i)
        ps.dihedrals[i].V = dih[i].value("V", ps.dihedrals[i].V);

    const json extra = fitted.value("extra_dihedrals", json::array());
    for (std::size_t i = 0; i < ps.extra_dihedrals.size() && i < extra.size(); ++i)
        ps.extra_dihedrals[i].V = extra[i].value("V", ps.extra_dihedrals[i].V);

    const json inv = fitted.value("inversions", json::array());
    for (std::size_t i = 0; i < ps.inversions.size() && i < inv.size(); ++i)
        ps.inversions[i].fc = inv[i].value("fc", ps.inversions[i].fc);
}

/// Raw Cartesian Hessian of GFN-FF evaluated through the PRODUCTION path with an
/// explicitly supplied parameter set.
///
/// A bare ForceField + setGFNFFParameters is NOT a GFN-FF: that function creates no
/// FFWorkspace and receives none of the state GFNFFComputationalMethod supplies, and was
/// measured 90% away from the real evaluation. Going through the wrapper instead means the
/// remainder H_nb and the verification both describe the force field the user actually runs.
///
/// Central differences of the analytic gradient. Geometry in Angstrom, gradient in
/// Eh/Angstrom (curcuma's convention), so the result is Hartree/Angstrom^2.
Matrix gfnffProductionHessian(const Mol& mol, const GFNFFParameterSet& ps,
    int threads, double step_angstrom, Vector* frozen_cn_out = nullptr,
    Matrix* frozen_dcn_out = nullptr)
{
    json controller;
    controller["threads"] = threads;
    controller["verbosity"] = 0;
    // Freeze the coordination numbers and EEQ charges. The fit models the force field with
    // its parameters held fixed, but GFN-FF normally lets CN (and hence the dynamic bond
    // r0) and the charges respond to every displacement. That response scales with the
    // bonded constants, so leaving it on makes the model and the force field disagree by a
    // k-dependent amount — measured as a 1.9e-1 verification residual. Frozen, the
    // production evaluation is exactly the quantity the fit represents.
    controller["static_cn"] = true;
    controller["static_charges"] = true;
    controller["gfnff"]["static_cn"] = true;
    controller["gfnff"]["static_charges"] = true;

    GFNFFComputationalMethod method("gfnff", controller);
    method.setExternalGFNFFParameters(ps);   // must precede setMolecule
    if (!method.setMolecule(mol))
        return Matrix();

    const int natoms = static_cast<int>(mol.m_atoms.size());
    const int n3 = 3 * natoms;
    Matrix hessian = Matrix::Zero(n3, n3);
    Matrix geometry = mol.m_geometry;

    // static_cn freezes the coordination numbers after the FIRST evaluation. Do that
    // evaluation at the REFERENCE geometry, so the frozen r0 is the one the unit Hessians
    // were built with — otherwise the state is captured at an already-displaced geometry.
    method.updateGeometry(geometry);
    method.calculateEnergy(true);
    if (frozen_cn_out && method.getGFNFF())
        *frozen_cn_out = method.getGFNFF()->getLastCN();
    if (frozen_dcn_out && method.getGFNFF()) {
        // The CN derivative is stored as a sparse apply-operator (diagonal + pair list).
        // Densify it column by column: applying it to the unit vector e_a yields the
        // (natoms x 3) field d cn_a / d x, which is one row of the matrix the fit wants.
        const CNDerivStore& store = method.getGFNFF()->getLastCNDerivatives();
        Matrix dcn = Matrix::Zero(natoms, 3 * natoms);
        if (!store.empty() && store.natoms == natoms) {
            Vector unit = Vector::Zero(natoms);
            GeoGradMatrix out(natoms, 3);
            for (int a = 0; a < natoms; ++a) {
                unit.setZero();
                unit(a) = 1.0;
                out.setZero();
                store.applyAdd(unit, out);
                for (int i = 0; i < natoms; ++i)
                    for (int c = 0; c < 3; ++c)
                        dcn(a, 3 * i + c) = out(i, c);
            }
        }
        *frozen_dcn_out = dcn;
    }

    for (int a = 0; a < natoms; ++a) {
        for (int c = 0; c < 3; ++c) {
            const double saved = geometry(a, c);

            // GFNFFComputationalMethod::getGradient() hands back GFNFF::Gradient()
            // unconverted, i.e. Hartree/BOHR — unlike the Hartree/Angstrom the rest of
            // curcuma passes around. Differencing it over an Angstrom step without this
            // factor leaves the Hessian a factor 1.8897 too small; measured directly as a
            // uniform ||dH_model|| / ||dH_real|| = 1.89 across bonds, angles and torsions.
            constexpr double grad_bohr_to_angstrom = 1.889726125;

            geometry(a, c) = saved + step_angstrom;
            method.updateGeometry(geometry);
            method.calculateEnergy(true);
            const Matrix gp = method.getGradient() * grad_bohr_to_angstrom;

            geometry(a, c) = saved - step_angstrom;
            method.updateGeometry(geometry);
            method.calculateEnergy(true);
            const Matrix gm = method.getGradient() * grad_bohr_to_angstrom;

            geometry(a, c) = saved;
            for (int b = 0; b < natoms; ++b)
                for (int d = 0; d < 3; ++d)
                    hessian(3 * a + c, 3 * b + d) = (gp(b, d) - gm(b, d)) / (2.0 * step_angstrom);
        }
    }
    return 0.5 * (hessian + hessian.transpose()).eval();
}

} // namespace

// =============================================================================
// Basin selection
// =============================================================================

std::vector<int> QMDFFFit::selectBasinFrames()
{
    std::vector<int> frames{ 0 };
    const int available = static_cast<int>(m_ensemble.size());
    if (available <= 1)
        return frames;

    if (!m_basin_frames.empty()) {
        frames.clear();
        std::stringstream stream(m_basin_frames);
        std::string item;
        while (std::getline(stream, item, ',')) {
            try {
                const int index = std::stoi(item);
                if (index >= 0 && index < available)
                    frames.push_back(index);
                else
                    CurcumaLogger::warn(fmt::format("basin frame {} out of range (0..{})",
                        index, available - 1));
            } catch (const std::exception&) {
                CurcumaLogger::warn(fmt::format("cannot parse basin frame '{}'", item));
            }
        }
        if (frames.empty())
            frames.push_back(0);
        return frames;
    }

    const int wanted = std::min(std::max(1, m_basins), available);
    if (wanted <= 1)
        return frames;

    // Spread the fitting points evenly over the input order. Without reference energies
    // for every frame there is no better criterion available here, and the sensitivity
    // study showed the NUMBER of fits matters more than their exact placement.
    frames.clear();
    for (int i = 0; i < wanted; ++i)
        frames.push_back(static_cast<int>(std::llround(
            static_cast<double>(i) * (available - 1) / std::max(1, wanted - 1))));
    frames.erase(std::unique(frames.begin(), frames.end()), frames.end());
    return frames;
}


// =============================================================================
// GFN-FF bonded-constant fit
// =============================================================================

void QMDFFFit::runGFNFFFit(double e0, const Vector& qm_frequencies)
{
    const Mol mol = m_molecule.getMolInfo();

    // 1) GFN-FF's own parameter generation — native structs, no JSON round trip
    json gfnff_controller = m_controller;
    gfnff_controller["method"] = "gfnff";
    gfnff_controller["threads"] = m_threads;

    GFNFFParameterSet parameters;
    Vector cn;
    {
        GFNFF gfnff(gfnff_controller);
        if (!gfnff.InitialiseMolecule(mol)) {
            CurcumaLogger::error("GFN-FF could not initialise the molecule");
            return;
        }
        gfnff.Calculation(true); // populates coordination numbers and EEQ charges
        parameters = gfnff.generateGFNFFParameterSet();
        cn = gfnff.getLastCN();
    }

    if (parameters.bonds.empty()) {
        CurcumaLogger::error("GFN-FF produced no bonded parameters");
        return;
    }

    const json jparams = gfnffParametersToJson(parameters);

    curcuma::qmdff::FitOptions options = fitOptions();
    options.source = curcuma::qmdff::TermSource::GFN_FF;

    QMDFFParametrisation parametrisation(m_molecule.Atoms(), m_molecule.getGeometry(),
        jparams, options);

    // 2) The non-fitted remainder, obtained BY DIFFERENCE. GFN-FF's non-covalent block has
    //    no JSON form, but linearity gives it for one force-field Hessian:
    //        H_nb = H_FF(k_current) - sum_p k_current_p U_p
    Vector cn_production;
    Matrix dcn_production;
    const Matrix h_ff_current = QMDFFParametrisation::hessianAngstromToBohr(
        gfnffProductionHessian(mol, parameters, m_threads, m_verify_fd_step,
            &cn_production, &dcn_production));

    // The coordination numbers and their derivative come from the run that produced
    // h_ff_current, not from the separate step-1 instance, so the fit is guaranteed to use
    // the same frozen state the force field did. Both REBUILD the term data, and both must
    // land before currentBondedHessian() below: the dynamic bond r0 depends on cn, and its
    // response dE/dcn * dcn/dx is part of the bond unit Hessian. Getting either wrong moves
    // the remainder instead, which the verification then reports.
    parametrisation.setCoordinationNumbers(cn_production.size() == cn.size() ? cn_production : cn);
    if (dcn_production.rows() == m_molecule.AtomCount())
        parametrisation.setCNDerivatives(dcn_production);
    if (h_ff_current.rows() != 3 * m_molecule.AtomCount()) {
        CurcumaLogger::error("Could not evaluate GFN-FF through the production path");
        return;
    }
    const Matrix h_bonded_current = parametrisation.currentBondedHessian();
    parametrisation.setExternalNonbonded(h_ff_current - h_bonded_current, 0.0,
        Matrix::Zero(m_molecule.AtomCount(), 3));

    // ---------------------------------------------------------------------------
    // Per-kind kernel check, independent of the remainder.
    //
    // Linearity says that scaling one term kind by (1+eps) must change the production
    // Hessian by exactly eps * sum_p k_p U_p over that kind. Comparing that against the
    // model isolates each kernel: no H_nb, no fitted constants, and any k-dependent
    // response the model omits (e.g. the dynamic r0 following the coordination numbers)
    // shows up as a deficit in exactly the kinds it affects.
    // ---------------------------------------------------------------------------
    {
        const double eps = 0.05;
        auto scaleKind = [&](curcuma::qmdff::TermKind kind) {
            GFNFFParameterSet sc = parameters;
            switch (kind) {
            case curcuma::qmdff::TermKind::Bond:
                // Skip the hydrogen-bridge bonds the fit also skips, otherwise the
                // reference change contains terms the model deliberately does not have.
                for (auto& b : sc.bonds)
                    if (b.nr_hb < 1)
                        b.fc *= (1.0 + eps);
                break;
            case curcuma::qmdff::TermKind::Angle:
                for (auto& a : sc.angles) a.fc *= (1.0 + eps);
                break;
            case curcuma::qmdff::TermKind::OutOfPlane:
                for (auto& i : sc.inversions) i.fc *= (1.0 + eps);
                break;
            default:
                for (auto& d : sc.dihedrals) d.V *= (1.0 + eps);
                for (auto& d : sc.extra_dihedrals) d.V *= (1.0 + eps);
                break;
            }
            return sc;
        };

        const Vector k_int = parametrisation.currentConstantsInternalPublic();
        for (auto kind : { curcuma::qmdff::TermKind::Bond, curcuma::qmdff::TermKind::Angle,
                 curcuma::qmdff::TermKind::Torsion, curcuma::qmdff::TermKind::OutOfPlane }) {
            Matrix d_model = Matrix::Zero(h_ff_current.rows(), h_ff_current.cols());
            int count = 0;
            for (const auto& t : parametrisation.termHessians()) {
                const bool match = (kind == curcuma::qmdff::TermKind::Torsion)
                    ? (t.kind == curcuma::qmdff::TermKind::Torsion
                        || t.kind == curcuma::qmdff::TermKind::ExtraTorsion)
                    : (t.kind == kind);
                if (!match)
                    continue;
                ++count;
                t.addTo(d_model, eps * k_int(t.group));
            }
            if (count == 0) {
                CurcumaLogger::param(fmt::format("gfnff_kernel_check_{}",
                                         curcuma::qmdff::termKindName(kind)),
                    "no terms of this kind were built");
                continue;
            }
            const Matrix h_scaled = QMDFFParametrisation::hessianAngstromToBohr(
                gfnffProductionHessian(mol, scaleKind(kind), m_threads, m_verify_fd_step));
            const Matrix d_real = h_scaled - h_ff_current;
            const double rel = (d_real - d_model).norm() / std::max(1e-30, d_real.norm());
            CurcumaLogger::param(fmt::format("gfnff_kernel_check_{}",
                                     curcuma::qmdff::termKindName(kind)),
                fmt::format("{} terms, rel = {:.3e}  (||dH_real|| = {:.3e}, ||dH_model|| = {:.3e})",
                    count, rel, d_real.norm(), d_model.norm()));

        }
    }

    // 3) fit
    FitReport report;
    const json fitted = parametrisation.fit(m_hessian, &report);
    {
        std::ofstream report_file(outputPath("gfnff_fit_report.json"));
        report_file << report.toJson();
    }
    CurcumaLogger::param("gfnff_fit_terms",
        fmt::format("{} bonds, {} angles, {} torsions, {} out-of-plane",
            report.n_bonds, report.n_angles, report.n_torsions, report.n_oop));

    // 4) verification — a force-field Hessian at the FITTED constants against the fit's own
    //    prediction. Because the remainder was taken at the OLD constants, any error in a
    //    unit Hessian survives as (k_new - k_old)(U - U_true) and shows up here.
    GFNFFParameterSet fitted_set = parameters;
    applyFittedToGFNFF(fitted, fitted_set);

    if (m_verify) {
        const Matrix h_ff_fitted = QMDFFParametrisation::hessianAngstromToBohr(
            gfnffProductionHessian(mol, fitted_set, m_threads, m_verify_fd_step));
        const double denom = std::max(1e-30, report.hessian_predicted.norm());
        const double rel = (h_ff_fitted - report.hessian_predicted).norm() / denom;
        CurcumaLogger::param("gfnff_fit_verification",
            fmt::format("||H_FF - H_predicted|| / ||H_predicted|| = {:.3e}", rel));
        if (rel > 1e-5) {
            CurcumaLogger::error(fmt::format(
                "Verification FAILED ({:.3e}). The model sum_p k_p U_p + H_nb does not "
                "reproduce the production GFN-FF Hessian at the fitted constants, so the "
                "constants are NOT trustworthy. On the reference molecules this check reaches "
                "~6e-09, so a failure means a real model/force-field mismatch. Read the "
                "per-kind gfnff_kernel_check_* lines above: they scale one term kind at a time "
                "and isolate which kernel disagrees. Known past causes, all fixed: the "
                "coordination numbers not reaching the term build (bonds fell back to the "
                "static r0_ij), the dynamic-r0 CN response missing from the bond unit Hessian, "
                "and getGradient() returning Hartree/Bohr where Hartree/Angstrom was assumed.",
                rel));
        } else {
            CurcumaLogger::success("GFN-FF unit terms verified against the force field");
        }
    }

    // 5) output — the native struct set, serialised for inspection
    {
        std::ofstream out(outputPath("gfnff_param.json"));
        out << gfnffParametersToJson(fitted_set);
    }
    CurcumaLogger::success(fmt::format("GFN-FF bonded constants fitted to the {} Hessian "
                                       "(R^2 = {:.4f}); non-covalent block untouched",
        m_method, report.r_squared));
}

// =============================================================================
// Driver
// =============================================================================

void QMDFFFit::start()
{
    CitationRegistry::cite("qmdff");

    if (m_potential != "qmdff" && m_potential != "quff" && m_potential != "gfnff") {
        CurcumaLogger::error(fmt::format(
            "-qmdfffit fits a force field to a QM Hessian; potential='{}' is not supported. "
            "Use 'qmdff' (QMDFF functional forms) or 'gfnff' (GFN-FF forms, non-covalent "
            "block kept as is).", m_potential));
        return;
    }

    CurcumaLogger::info(fmt::format("QMDFF parametrisation: reference method '{}', {} atoms",
        m_method, m_molecule.AtomCount()));

    // --- 1) charges ---------------------------------------------------------
    double e0 = 0.0, gmax = 0.0;
    Vector charges;
    if (!prepareCharges(e0, charges, gmax))
        return;

    if (gmax > 1.0e-3) {
        CurcumaLogger::warn(fmt::format(
            "|gradient|max = {:.3e} Eh/Angstrom — this geometry is not a stationary point of "
            "'{}'. QMDFF takes r0/theta0 from it, so the force field's minimum is defined "
            "somewhere the reference's is not and the fitted constants will be biased. "
            "Run -opt at the same level first.", gmax, m_method));
    }

    // --- 2) reference Hessian ----------------------------------------------
    Vector qm_frequencies;
    if (!prepareHessian(m_hessian, qm_frequencies))
        return;

    // --- 3) GFN-FF branch ----------------------------------------------------
    // Keep GFN-FF's non-covalent block (D4 dispersion, EEQ electrostatics, HB/XB) exactly
    // as it is — that is the part carrying ~80% of the conformational spread and the part
    // GFN-FF does better than QMDFF — and fit only the bonded constants to the Hessian.
    if (m_potential == "gfnff") {
        runGFNFFFit(e0, qm_frequencies);
        return;
    }

    // --- 3) initial QMDFF parameter set -------------------------------------
    m_controller["method"] = "qmdff";
    ConfigManager ff_config("forcefield", m_controller);
    ForceFieldGenerator generator(ff_config);
    generator.setMolecule(m_molecule.getMolInfo()); // carries the charges -> qq, HB/XB
    generator.Generate();

    json parameter = generator.getParameter();
    parameter["method"] = "qmdff"; // NOT "quff": the cache method check compares this
    parameter["e0"] = e0;
    {
        std::ofstream init_file(outputPath("qmdff_init_param.json"));
        init_file << parameter;
    }

    // --- 4) the fit — no force-field evaluation at all ------------------------
    QMDFFParametrisation parametrisation(m_molecule.Atoms(), m_molecule.getGeometry(),
        parameter, fitOptions());
    FitReport report;
    json fitted;

    std::vector<int> basin_frames = selectBasinFrames();
    if (basin_frames.size() <= 1) {
        fitted = parametrisation.fit(m_hessian, &report);
    } else {
        // Several basins at once: one parameter set constrained by the curvature AND the
        // relative energies AND the gradients of k structures. The relative-energy rows
        // are the ones that can fix conformer ranking; a Hessian alone cannot.
        std::vector<curcuma::qmdff::BasinData> basins;
        for (std::size_t i = 0; i < basin_frames.size(); ++i) {
            const int frame = basin_frames[i];
            curcuma::qmdff::BasinData basin;
            basin.geometry_angstrom = m_ensemble[frame].getGeometry();

            if (i == 0) {
                // basin 0 is the structure the topology came from — reuse its Hessian
                basin.hessian_angstrom = m_hessian;
                basin.energy = e0;
                basin.has_energy = true;
            } else {
                CurcumaLogger::info(fmt::format(
                    "basin {}/{}: {} Hessian at frame {} ({} gradient calls)",
                    i + 1, basin_frames.size(), m_method, frame, 6 * m_molecule.AtomCount()));
                Molecule mol = m_ensemble[frame];
                Matrix raw;
                {
                    Hessian hessian_calc(m_method, m_defaults, true);
                    hessian_calc.setMolecule(mol);
                    hessian_calc.start();
                    raw = hessian_calc.getRawHessian();
                }
                basin.hessian_angstrom = raw;

                EnergyCalculator energy(m_method, m_defaults);
                energy.setMolecule(mol.getMolInfo());
                basin.energy = energy.CalculateEnergy(true);
                basin.has_energy = true;
                basin.gradient_angstrom = energy.Gradient();
            }
            basins.push_back(basin);
        }
        fitted = parametrisation.fitMultiBasin(basins, &report);
    }

    {
        std::ofstream report_file(outputPath("qmdff_fit_report.json"));
        report_file << report.toJson();
    }

    CurcumaLogger::param("qmdff_fit_terms",
        fmt::format("{} bonds, {} angles, {} torsions, {} out-of-plane",
            report.n_bonds, report.n_angles, report.n_torsions, report.n_oop));

    // --- 5) verification ------------------------------------------------------
    if (m_verify)
        verifyFit(fitted, report, qm_frequencies);

    // --- 6) output ------------------------------------------------------------
    {
        std::ofstream parameterfile(outputPath("qmdff_param.json"));
        parameterfile << fitted;
    }

    if (m_write_param && !m_input_file.empty()) {
        const std::string auto_file = ForceField::generateParameterFileName(m_input_file);
        std::ofstream auto_out(auto_file);
        if (auto_out) {
            auto_out << fitted;
            CurcumaLogger::success(fmt::format(
                "Fitted parameters written to {} — 'curcuma -sp {} -method qmdff' will now "
                "load them automatically", auto_file, m_input_file));
        } else {
            CurcumaLogger::warn(fmt::format("Could not write {}", auto_file));
        }
    }
}
