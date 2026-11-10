/*
 * <NEB-like Molecular Dynamics - band of independently-thermostatted replicas
 *    coupled by spring forces along a reaction path, driven by MD.>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested pending - not TESTED/APPROVED).
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
 * See nebmd.h for the force model / timing caveat and docs/NEB_MD.md for usage.
 */

#include "src/capabilities/nebmd.h"
#include "src/core/units.h"
#include "src/capabilities/wham.h"
#include "src/capabilities/string_path.h"
#include "src/capabilities/z_filter.h"
#include "src/capabilities/fe_profile.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <utility>

#include <fmt/format.h>

#include "src/core/citation_registry.h"
#include "src/core/curcuma_logger.h"
#include "src/core/fileiterator.h"
#include "src/tools/formats.h"   // Files::LoadFile (transitively via general.h too)

// Claude Generated 2026: RAII guard that temporarily redirects std::cout to a null
// buffer so the per-image SimpleMD setup chatter (constructor -> LoadControlJson,
// Initialise, prepareRun) does not spam the console N times. std::cerr is left
// untouched so real errors stay visible. Exceptions are safe (dtor restores).
namespace {
class CoutSuppress {
public:
    CoutSuppress() : m_old(std::cout.rdbuf()) { std::cout.rdbuf(m_sink.rdbuf()); }
    ~CoutSuppress() { std::cout.rdbuf(m_old); }
private:
    std::streambuf* m_old;
    std::ostringstream m_sink;
};
} // namespace

NebMD::NebMD(const json& controller, bool silent)
    : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("nebmd"), controller, silent)
    , m_config("nebmd", controller)
{
    UpdateController(controller); // calls LoadControlJson() internally
}

NebMD::~NebMD()
{
    // Destroy in reverse construction order so the CurcumaMethod RAII verbosity
    // save/restore stack of the per-image SimpleMD instances stays balanced.
    for (auto it = m_images.rbegin(); it != m_images.rend(); ++it)
        if (*it) delete *it;
}

void NebMD::LoadControlJson()
{
    m_nimages = m_config.get<int>("nimages");
    if (m_nimages < 3) {
        CurcumaLogger::warn("NEB-MD: nimages < 3 makes no sense; raising to 3.");
        m_nimages = 3;
    }
    m_nudged = m_config.get<bool>("nudged");
    m_fixed_endpoints = m_config.get<bool>("fixed_endpoints");
    m_mass_weighted_spring = m_config.get<bool>("mass_weighted_spring");
    m_k_spring = m_config.get<double>("k_spring");
    m_dump = m_config.get<int>("dump_frequency");
    if (m_dump < 1) m_dump = 1;
    m_write_initial_path = m_config.get<bool>("write_initial_path");
    m_optimize = m_config.get<bool>("optimize");
    m_opt_iterations = m_config.get<int>("opt_iterations");
    m_opt_fmax = m_config.get<double>("opt_fmax");
    m_climbing = m_config.get<bool>("climbing");
    m_opt_k = m_config.get<double>("opt_k");
    m_path_file = m_config.get<std::string>("path_file", "");
    m_umbrella = m_config.get<bool>("umbrella");
    m_umbrella_k = m_config.get<double>("umbrella_k");
    m_umbrella_bins = m_config.get<int>("umbrella_bins");
    m_kz = m_config.get<double>("umbrella_kz");
    m_z_tolerance = m_config.get<double>("umbrella_z_tolerance");
    m_z_wall = std::max(m_config.get<double>("z_wall"), 0.0);       // Claude Generated 2026
    m_z_wall_k = std::max(m_config.get<double>("z_wall_k"), 0.0);
    m_window_trajectory = m_config.get<bool>("window_trajectory");
    m_umbrella_keep_springs = m_config.get<bool>("umbrella_keep_springs");
    m_reparametrize = m_config.get<int>("reparametrize");
    m_pmf = m_config.get<bool>("pmf");
    m_fixman = m_config.get<bool>("fixman");
    m_pmf_equilibration = m_config.get<int>("pmf_equilibration");

    // Decoupled scheme (Claude Generated 2026). "independent" is a thin, explicit layer
    // over the umbrella machinery: restrain s only, drop the neighbour springs, and switch
    // on the conditional mean-force estimator. It is named separately because it is a
    // different method, not a different parameter setting.
    m_scheme = m_config.get<std::string>("scheme");
    m_independent = (m_scheme == "independent");
    m_fe_bins = std::max(m_config.get<int>("fe_bins"), 4);
    m_div_probes = std::max(m_config.get<int>("div_probes"), 0);
    m_div_epsilon = m_config.get<double>("div_epsilon");
    m_ramp_steps = m_config.get<int>("ramp_steps");
    if (m_ramp_steps <= 0)
        m_ramp_steps = m_pmf_equilibration;
    m_bias_fmax = m_config.get<double>("bias_fmax");
    m_cv_trace = m_config.get<bool>("cv_trace");
    m_window_threads = std::max(m_config.get<int>("window_threads"), 1);
    m_max_temperature = m_config.get<double>("max_temperature");
    m_constrain = m_config.get<bool>("constrain");
    m_constrain_tol = m_config.get<double>("constrain_tolerance");
    m_adaptive_k = m_config.get<bool>("adaptive_k");
    m_adaptive_k_max = m_config.get<double>("adaptive_k_max");
    m_band_verbosity = CurcumaLogger::get_verbosity();
    // Finite-temperature string (Claude Generated 2026): restraint to the mean path only.
    m_string = (m_scheme == "string");
    m_string_k = m_config.get<double>("string_k");
    m_string_update = std::max(m_config.get<int>("string_update"), 1);
    m_string_alpha = m_config.get<double>("string_alpha");
    m_string_smooth = m_config.get<double>("string_smooth");
    m_string_tol = m_config.get<double>("string_tol");
    m_string_patience = std::max(m_config.get<int>("string_patience"), 1);
    m_string_max_time = m_config.get<double>("string_max_time");
    m_string_max_step = m_config.get<double>("string_max_step");          // Claude Generated 2026
    m_string_abort_on_dead = m_config.get<bool>("string_abort_on_dead"); // Claude Generated 2026
    m_string_prepare = m_config.get<bool>("string_prepare");             // Claude Generated 2026
    m_string_ftol = m_config.get<double>("string_ftol");
    m_string_blocks = std::max(m_config.get<int>("string_blocks"), 2);
    m_string_gain = m_config.get<double>("string_gain");
    m_string_gain_z = m_config.get<double>("string_gain_z");
    m_string_write_updates = m_config.get<bool>("string_write_updates");
    m_cv_histogram = m_config.get<std::string>("cv_histogram", "");
    m_cv_hist_zmax = m_config.get<double>("cv_hist_zmax");
    m_cv_frames = std::max(m_config.get<int>("cv_frames"), 0);   // Claude Generated 2026
    m_cv_extend = std::max(m_config.get<int>("cv_extend"), 0);   // Claude Generated 2026
    // Claude Generated 2026: extension as a fraction of the CV path length; when > 0 it takes
    // precedence over cv_extend (frames), the frame count is derived in buildCvPath.
    m_cv_extend_length = std::max(m_config.get<double>("cv_extend_length"), 0.0);
    m_cv_atoms_spec = m_config.get<std::string>("cv_atoms", "all"); // parsed in Initialise (needs the atoms)
    // Claude Generated 2026: distance-space path CV (validated in Initialise).
    m_cv_metric_name = m_config.get<std::string>("cv_metric", "rmsd");
    m_cv_pairs_spec = m_config.get<std::string>("cv_pairs", "auto");
    m_cv_pair_threshold = m_config.get<double>("cv_pair_threshold", 0.5);
    {
        const std::string t = m_config.get<std::string>("cv_pair_transform", "d");
        m_cv_transform = (t == "inv") ? PathCV::PairTransform::Inverse : PathCV::PairTransform::Distance;
    }
    if (m_cv_extend_length > 0.0 && m_cv_extend > 0)
        CurcumaLogger::warn(fmt::format("NEB-MD: both cv_extend ({}) and cv_extend_length ({}) are set - "
                                        "cv_extend_length takes precedence, cv_extend is ignored.",
                                        m_cv_extend, m_cv_extend_length));
    if (m_string) {
        if (m_umbrella || m_constrain)
            CurcumaLogger::warn("NEB-MD: scheme 'string' ignores umbrella/constrain (no path CV in stage 1).");
        m_umbrella = false;
        m_constrain = false;
        m_umbrella_keep_springs = false;
        if (m_reparametrize > 0)
            CurcumaLogger::warn("NEB-MD: -nebmd.reparametrize is ignored in scheme 'string' (only the "
                                "mean path is reparametrised, never the images).");
        m_reparametrize = 0;
    }
    if (m_independent) {
        m_umbrella = true;              // uses the frozen-path CV and its static restraint
        m_umbrella_keep_springs = false; // no neighbour coupling - that is the whole point
    }
    m_mtd_transverse = m_config.get<bool>("mtd_transverse");
    m_idpp = m_config.get<bool>("idpp");
    m_idpp_iterations = m_config.get<int>("idpp_iterations");
    m_repair_overlaps = m_config.get<bool>("repair_overlaps");
    m_pt = m_config.get<int>("proton_transfer");
    m_force_reorder = m_config.get<bool>("force_reorder");
    m_threads = m_config.get<int>("threads");
    // Claude Generated 2026: the energy method may arrive as a global (-method), as
    // -nebmd.method, or - most commonly - inside the MD sub-config (-md.method /
    // -simplemd.method), because that is where users put the MD settings. The band
    // optimisation and the endpoint single points must use the SAME method as the
    // images, otherwise a "-md.method gfn2" run would silently optimise with the
    // gfnff default (observed: GFN2 and GFN-FF returning a bit-identical barrier).
    m_method = m_config.get<std::string>("method", "gfnff");
    for (const char* scope : { "simplemd", "md" }) {
        if (m_controller.contains(scope) && m_controller[scope].is_object()
            && m_controller[scope].contains("method")
            && m_controller[scope]["method"].is_string()) {
            m_method = m_controller[scope]["method"].get<std::string>();
            break; // simplemd wins over md (it is the resolved module)
        }
    }

    // MD sub-config forwarded to every image's SimpleMD. Falls back to empty if
    // the user only passed flat globals (-method/-time_step/...); the per-image
    // ConfigManager("simplemd") picks those up via its flat-key backward-compat.
    m_md_config = m_controller.value("simplemd", json::object());
}

void NebMD::printHelp() const
{
    ParameterRegistry::getInstance().printHelp("nebmd");
    std::cout << "\nExample usage:\n"
              << "  curcuma -nebmd start.xyz end.xyz -nebmd.nimages 12 \\\n"
              << "      -md.method gfnff -md.time_step 0.5 -md.max_time 200 -md.temperature 300 \\\n"
              << "      -nebmd.k_spring 0.005 -nebmd.dump_frequency 20\n\n"
              << "MD parameters for the images: -simplemd.* (preferred), -md.*, or bare\n"
              << "-method/-time_step/-max_time/-thermostat all work (-temperature is multi-owner,\n"
              << "use -md.temperature/-simplemd.temperature). Energy method default: gfnff.\n"
              << "See docs/NEB_MD.md for the force model and (k, dt) guidance.\n"
              << std::endl;
}

void NebMD::setEndpoints(const Molecule& start, const Molecule& end)
{
    m_start = start;
    m_end = end;
    m_endpoints_set = true;
}

bool NebMD::Initialise()
{
    if (m_pmf) {
        CurcumaLogger::error("NEB-MD: -nebmd.pmf was removed (the spring-band mean force has no "
                             "well-defined coordinate). Use -nebmd.umbrella, optionally with "
                             "-nebmd.constrain.");
        return false;
    }
    if (!m_endpoints_set) {
        CurcumaLogger::error("NEB-MD: endpoints not set - call setEndpoints() first.");
        return false;
    }
    if (m_start.AtomCount() == 0 || m_end.AtomCount() == 0) {
        CurcumaLogger::error("NEB-MD: empty endpoint structure.");
        return false;
    }
    if (m_start.AtomCount() != m_end.AtomCount()) {
        CurcumaLogger::error("NEB-MD: endpoint atom counts differ ("
                           + std::to_string(m_start.AtomCount()) + " vs "
                           + std::to_string(m_end.AtomCount()) + ").");
        return false;
    }
    // Element multiset must match (same stoichiometry).
    std::vector<int> a = m_start.Atoms();
    std::vector<int> b = m_end.Atoms();
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    if (a != b) {
        CurcumaLogger::error("NEB-MD: endpoints have different element compositions.");
        return false;
    }

    // Claude Generated 2026: atom subset of the path CV (cv_atoms); empty = all atoms.
    {
        std::string err;
        if (!PathCV::parseAtomSelection(m_cv_atoms_spec, m_start, m_cv_atoms, err)) {
            CurcumaLogger::error("NEB-MD: " + err);
            return false;
        }
        const int nat = static_cast<int>(m_start.AtomCount());
        m_cv_w = PathCV::selectionWeights(m_cv_atoms, nat);
        m_cv_w_align = PathCV::selectionWeights(m_cv_atoms, nat, 1e-6); // see buildCvPath
        CurcumaLogger::result_fmt("NEB-MD: path CV atoms (cv_atoms) = {}: {} of {} atoms", m_cv_atoms_spec,
                                  m_cv_atoms.empty() ? nat : static_cast<int>(m_cv_atoms.size()), nat);

        // Claude Generated 2026: metric of the path CV (rmsd | distances) and its pair list.
        const std::string tr = m_config.get<std::string>("cv_pair_transform", "d");
        if (m_cv_metric_name != "rmsd" && m_cv_metric_name != "distances") {
            CurcumaLogger::error("NEB-MD: cv_metric '" + m_cv_metric_name + "' is unknown (rmsd | distances)");
            return false;
        }
        if (tr != "d" && tr != "inv") {
            CurcumaLogger::error("NEB-MD: cv_pair_transform '" + tr + "' is unknown (d | inv)");
            return false;
        }
        m_cv_metric = (m_cv_metric_name == "distances") ? PathCV::Metric::Distances : PathCV::Metric::Rmsd;
        m_cv_pairs.clear();
        m_cv_pairs_resolved = false;
        if (m_cv_metric == PathCV::Metric::Distances) {
            // cv_atoms only restricts the 'auto' candidates here: arclength and alignment of the
            // CV frames run over all atoms (the distance metric itself needs no fit).
            m_cv_w.clear();
            m_cv_w_align.clear();
            std::string spec = m_cv_pairs_spec;
            spec.erase(0, spec.find_first_not_of(" \t"));
            spec.erase(spec.find_last_not_of(" \t") + 1);
            if (spec != "auto") {
                if (!PathCV::parsePairList(spec, nat, m_cv_pairs, err)) {
                    CurcumaLogger::error("NEB-MD: " + err);
                    return false;
                }
                m_cv_pairs_resolved = true;
                if (!m_cv_atoms.empty())
                    CurcumaLogger::info("NEB-MD: explicit cv_pairs override cv_atoms.");
            }
            if (m_cv_pair_threshold < 0.0) {
                CurcumaLogger::error("NEB-MD: cv_pair_threshold must be >= 0");
                return false;
            }
            CurcumaLogger::result_fmt("NEB-MD: path CV metric = distances (x = {}, z in {}), cv_pairs = {}",
                                      m_cv_transform == PathCV::PairTransform::Inverse ? "1/d" : "d",
                                      m_cv_transform == PathCV::PairTransform::Inverse ? "A^-2" : "A^2", spec);
        } else {
            CurcumaLogger::result("NEB-MD: path CV metric = rmsd (z in A^2)");
        }
    }

    m_images.assign(m_nimages, nullptr);
    m_image_done.assign(m_nimages, false);
    m_image_geoms.assign(m_nimages, Geometry());
    resetPmfAccumulators(static_cast<int>(m_start.AtomCount()));
    return true;
}

void NebMD::resetPmfAccumulators(int /*nat*/)
{
    // Claude Generated 2026: single place where every PMF running sum is (re)initialised.
    // Previously duplicated in two call sites, which is how a new accumulator gets
    // forgotten in one of them.
    m_fe_count.assign(m_fe_bins, 0);
    m_fe_force.assign(m_fe_bins, 0.0);
    m_fe_div.assign(m_fe_bins, 0.0);
    m_fe_wall.assign(m_fe_bins, 0.0);
    m_fe_ssum.assign(m_fe_bins, 0.0);
    m_fe_wsum.assign(m_fe_bins, 0.0);
    m_win_s_sum.assign(m_nimages, 0.0);
    m_win_s2_sum.assign(m_nimages, 0.0);
    m_win_count.assign(m_nimages, 0);
    m_win_z_sum.assign(m_nimages, 0.0);
    m_win_z_max.assign(m_nimages, -1e300);
    m_win_zwall_n.assign(m_nimages, 0);
    m_win_wall_sum.assign(m_nimages, 0.0);
    m_win_w.assign(m_nimages, 0.0);
    m_win_wf.assign(m_nimages, 0.0);
    m_win_wd.assign(m_nimages, 0.0);
    m_win_wwall.assign(m_nimages, 0.0);
    m_prev_s.assign(m_nimages, -1.0);
    m_img_max_jump.assign(m_nimages, 0.0);
    m_img_max_jump_step.assign(m_nimages, -1);
    m_img_capped.assign(m_nimages, 0);
    m_window_died.assign(m_nimages, -1);
    m_window_died_T.assign(m_nimages, 0.0);
    m_calib_dg_sum.assign(m_nimages, 0.0);
    m_calib_n.assign(m_nimages, 0);
    m_fe_rejected_bin.assign(m_fe_bins, 0);
    m_win_rejected.assign(m_nimages, 0);
    m_probe_rng.clear();
    for (int i = 0; i < m_nimages; ++i)
        m_probe_rng.emplace_back(static_cast<std::mt19937::result_type>(20260802u + 7919u * (i + 1)));
}

Geometry NebMD::interpolate(const Geometry& a, const Geometry& b, double t) const
{
    return (1.0 - t) * a + t * b;
}

double NebMD::idppObjective(const Geometry& R, const Matrix& d_target, Geometry* grad) const
{
    // Claude Generated 2026: IDPP objective (Smidstrup et al., JCP 140, 214106 (2014), Eq. 7)
    //
    //   S(R) = sum_{a<b} w(d_ab) * (d_target_ab - d_ab)^2 ,  w(d) = 1/d^4
    //
    // The 1/d^4 weight is what makes this "repulsion aware": a pair that is much closer
    // than its target distance contributes enormously, so the optimiser pushes atoms
    // apart instead of letting them interpenetrate (which plain cartesian interpolation
    // happily does). Long-range pairs are damped, so the objective is dominated by the
    // local packing rather than by the overall molecular size.
    //
    // Gradient: dS/dR_a = sum_b [ -2*w*(dt-d) - 4*(dt-d)^2/d^5 ] * (R_a - R_b)/d
    const int nat = static_cast<int>(R.rows());
    double S = 0.0;
    if (grad)
        *grad = Geometry::Zero(nat, 3);

    for (int a = 0; a < nat; ++a) {
        for (int b = a + 1; b < nat; ++b) {
            Position rij = R.row(a).transpose() - R.row(b).transpose();
            double d = rij.norm();
            if (d < 1e-6)
                d = 1e-6;
            const double dt = d_target(a, b);
            const double diff = dt - d;
            const double d2 = d * d;
            const double d4 = d2 * d2;
            S += diff * diff / d4;
            if (grad) {
                // d/dd [ (dt-d)^2 / d^4 ] = -2(dt-d)/d^4 - 4(dt-d)^2/d^5
                const double dSdd = -2.0 * diff / d4 - 4.0 * diff * diff / (d4 * d);
                const Position g = (dSdd / d) * rij; // chain rule: dd/dR_a = rij/d
                grad->row(a) += g.transpose();
                grad->row(b) -= g.transpose();
            }
        }
    }
    return S;
}

bool NebMD::buildIDPPPath()
{
    // Claude Generated 2026: IDPP initial-path refinement.
    //
    // Linear cartesian interpolation ignores what the atoms do to each other: for any
    // non-trivial conformational change the intermediate images end up with atoms driven
    // into one another. IDPP instead interpolates the DISTANCE matrix linearly (which is
    // a far more sensible interpolation of a molecular structure) and then optimises each
    // image to reproduce its target distance matrix, weighted by 1/d^4 so that short
    // contacts dominate the objective.
    //
    // Neighbour coupling: this is run as a NEB on the IDPP surface (as in the paper), i.e.
    // each image additionally feels a spring to its neighbours, projected with the same
    // improved tangent used by the real band. So the images stay evenly spread and the
    // path stays connected while the repulsion is resolved.
    CitationRegistry::cite("idpp");

    const int nat = static_cast<int>(m_start_geom.rows());
    if (nat < 2)
        return true;

    // Endpoint distance matrices.
    Matrix d_start = Matrix::Zero(nat, nat);
    Matrix d_end = Matrix::Zero(nat, nat);
    for (int a = 0; a < nat; ++a) {
        for (int b = a + 1; b < nat; ++b) {
            const double ds = (m_start_geom.row(a) - m_start_geom.row(b)).norm();
            const double de = (m_end_geom.row(a) - m_end_geom.row(b)).norm();
            d_start(a, b) = d_start(b, a) = ds;
            d_end(a, b) = d_end(b, a) = de;
        }
    }

    // Per-image target distance matrix: linear interpolation in distance space.
    std::vector<Matrix> d_target(m_nimages);
    for (int i = 0; i < m_nimages; ++i) {
        const double t = (m_nimages == 1) ? 0.0
                                          : static_cast<double>(i) / static_cast<double>(m_nimages - 1);
        d_target[i] = (1.0 - t) * d_start + t * d_end;
    }

    double S_before = 0.0;
    for (int i = 1; i < m_nimages - 1; ++i)
        S_before += idppObjective(m_image_geoms[i], d_target[i], nullptr);

    // Steepest descent with a capped step on the IDPP + spring surface.
    const int max_iter = m_idpp_iterations;
    const double max_disp = 0.05;   // Angstrom per iteration
    const double k_idpp = 1.0;      // IDPP spring constant (objective units)

    for (int iter = 0; iter < max_iter; ++iter) {
        // Tangents from the CURRENT path (energy-free: use the plain bisecting tangent,
        // since the IDPP objective is not the physical energy).
        std::vector<Geometry> tau(m_nimages, Geometry::Zero(nat, 3));
        for (int i = 1; i < m_nimages - 1; ++i) {
            Geometry tp = m_image_geoms[i + 1] - m_image_geoms[i];
            Geometry tm = m_image_geoms[i] - m_image_geoms[i - 1];
            const double np = tp.norm();
            const double nm = tm.norm();
            if (np > 1e-12) tp /= np;
            if (nm > 1e-12) tm /= nm;
            Geometry t = tp + tm;
            const double n = t.norm();
            tau[i] = (n > 1e-12) ? Geometry(t / n) : Geometry(Geometry::Zero(nat, 3));
        }

        double gmax_all = 0.0;
        std::vector<Geometry> force(m_nimages, Geometry::Zero(nat, 3));
        for (int i = 1; i < m_nimages - 1; ++i) {
            Geometry g = Geometry::Zero(nat, 3);
            idppObjective(m_image_geoms[i], d_target[i], &g);
            Geometry F = -g; // steepest descent on the IDPP objective

            // Nudging: keep only the component perpendicular to the path ...
            double dot = 0.0;
            for (int a = 0; a < nat; ++a)
                dot += F.row(a).dot(tau[i].row(a));
            F -= dot * tau[i];

            // ... and add the parallel spring that keeps the images evenly spaced.
            const double dp = (m_image_geoms[i + 1] - m_image_geoms[i]).norm();
            const double dm = (m_image_geoms[i] - m_image_geoms[i - 1]).norm();
            F += k_idpp * (dp - dm) * tau[i];

            force[i] = F;
            for (int a = 0; a < nat; ++a)
                gmax_all = std::max(gmax_all, F.row(a).norm());
        }

        if (gmax_all < 1e-8)
            break;
        const double scale = std::min(max_disp / gmax_all, max_disp);
        for (int i = 1; i < m_nimages - 1; ++i)
            m_image_geoms[i] += scale * force[i];
    }

    double S_after = 0.0;
    for (int i = 1; i < m_nimages - 1; ++i)
        S_after += idppObjective(m_image_geoms[i], d_target[i], nullptr);

    if (getVerbosity() >= 1) {
        CurcumaLogger::result_fmt("NEB-MD: IDPP path refinement ({} iterations): objective {:.4e} -> {:.4e}",
                                  max_iter, S_before, S_after);
    }
    return S_after <= S_before;
}

int NebMD::repairPathOverlaps()
{
    // Claude Generated 2026: linear interpolation between two endpoints routinely drives
    // non-bonded atoms into each other at the intermediate images (the classic "broken
    // interpolated geometry"). With GFN-FF that shows up as extra bonds in the topology
    // and an image that relaxes BELOW the endpoints. Same remedy as
    // PolymerBuild::resolveOverlaps (polymerbuild.cpp): steepest descent on a soft
    // (r_cut/r)^4 clash potential over non-bonded pairs, with a capped displacement.
    //
    // The reference bond set is taken from image 0 (= the start structure), so genuinely
    // bonded neighbours are never pushed apart. Endpoints are left untouched.
    const int nat = static_cast<int>(m_start_geom.rows());
    const std::vector<int> atoms = m_start.Atoms();

    std::set<std::pair<int, int>> bonded;
    {
        Molecule ref = m_start;
        ref.setGeometry(m_image_geoms[0]);
        const auto [dist0, topo0] = ref.DistanceMatrix();
        for (int i = 0; i < topo0.rows(); ++i)
            for (int j = i + 1; j < topo0.cols(); ++j)
                if (topo0(i, j) > 0.5)
                    bonded.insert(std::make_pair(i, j));
    }

    const double max_displacement = 0.05; // Angstrom per step (conservative)
    const int max_steps = 200;
    int repaired = 0;

    for (int img = 1; img < m_nimages - 1; ++img) {
        Geometry coords = m_image_geoms[img];
        int clashes_initial = 0;

        for (int step = 0; step < max_steps; ++step) {
            Geometry grad = Geometry::Zero(nat, 3);
            int clashes = 0;

            for (int i = 0; i < nat; ++i) {
                for (int j = i + 1; j < nat; ++j) {
                    if (bonded.count(std::make_pair(i, j)))
                        continue; // bonded pairs are supposed to be close
                    Position diff = coords.row(i).transpose() - coords.row(j).transpose();
                    double d = diff.norm();
                    if (d < 1e-6)
                        d = 1e-6;
                    const double r_cov = Elements::CovalentRadius[atoms[i]] + Elements::CovalentRadius[atoms[j]];
                    const double r_cut = r_cov * 1.15; // matches the bond-detection cutoff
                    if (d >= r_cut)
                        continue;
                    ++clashes;
                    const double ratio = r_cut / d;
                    const double r4 = ratio * ratio * ratio * ratio;
                    const double gscale = 4.0 * r4 / (d * d);
                    const Position g = gscale * diff;
                    grad.row(i) += g.transpose();
                    grad.row(j) -= g.transpose();
                }
            }

            if (step == 0)
                clashes_initial = clashes;
            if (clashes == 0)
                break;

            double gmax = 0.0;
            for (int i = 0; i < nat; ++i)
                gmax = std::max(gmax, grad.row(i).norm());
            if (gmax < 1e-9)
                break;
            const double scale = std::min(max_displacement / gmax, max_displacement);
            coords += scale * grad;
        }

        if (clashes_initial > 0) {
            m_image_geoms[img] = coords;
            ++repaired;
            if (getVerbosity() >= 2)
                CurcumaLogger::info("NEB-MD: image " + std::to_string(img) + " had "
                    + std::to_string(clashes_initial) + " interpolation clash(es) - geometry relaxed.");
        }
    }

    if (repaired > 0)
        CurcumaLogger::warn("NEB-MD: repaired interpolation overlaps in " + std::to_string(repaired)
            + " image(s) (non-bonded atoms pushed apart). The initial path is an approximation - "
            + "check the final path, and prefer a pre-built path for demanding cases.");
    return repaired;
}

bool NebMD::checkPathTopology() const
{
    // Claude Generated 2026: GFN-FF derives bond topology from the geometry, so linear
    // interpolation between two endpoints can stretch a bond past the detection threshold
    // at an intermediate image (or bring two atoms into bonding range that are not bonded
    // at the endpoints). That image then gets a DIFFERENT topology -> discontinuous
    // energy/gradient along the band and a meaningless "minimum-energy path". This check
    // compares the bond-connectivity signature of every image against image 0 and warns.
    // For bond-breaking/forming reactions use a QM method (gfn2) or supply a better path.
    auto bondSignature = [](const Geometry& g, const Molecule& templ) -> std::set<std::pair<int, int>> {
        Molecule mol = templ;
        mol.setGeometry(g); // invalidates the distance/topology cache -> recomputed for g
        const auto [dist, topo] = mol.DistanceMatrix();
        std::set<std::pair<int, int>> bonds;
        for (int i = 0; i < topo.rows(); ++i)
            for (int j = i + 1; j < topo.cols(); ++j)
                if (topo(i, j) > 0.5)
                    bonds.insert(std::make_pair(i, j));
        return bonds;
    };

    const std::set<std::pair<int, int>> ref = bondSignature(m_image_geoms[0], m_start);
    for (int i = 1; i < m_nimages; ++i) {
        const auto s = bondSignature(m_image_geoms[i], m_start);
        if (s != ref) {
            CurcumaLogger::warn("NEB-MD: image " + std::to_string(i) + " has a different bond topology "
                + "than image 0 (" + std::to_string(s.size()) + " bonds vs " + std::to_string(ref.size())
                + "). Linear interpolation changed the connectivity, so the GFN-FF energy is "
                + "discontinuous along the band. For bond-breaking/forming reactions use a QM "
                + "method (e.g. gfn2) or supply a better initial path.");
            return false;
        }
    }
    return true;
}

bool NebMD::buildInitialPath()
{
    // Claude Generated 2026: an externally converged path (e.g. ORCA's neb_MEP_trj.xyz)
    // replaces the interpolation entirely. Sampling is only as good as the path it runs
    // on, and a linear/IDPP guess can converge to the wrong mechanism (documented for
    // helicene), so being able to hand in a real MEP is essential.
    if (!m_path_file.empty()) {
        std::vector<Geometry> frames;
        FileIterator it(m_path_file);
        while (!it.AtEnd()) {
            Molecule m = it.Next();
            if (m.AtomCount() == m_start.AtomCount())
                frames.push_back(m.getGeometry());
        }
        if (frames.size() < 3) {
            CurcumaLogger::error("NEB-MD: path_file '" + m_path_file + "' has fewer than 3 usable frames.");
            return false;
        }
        m_nimages = static_cast<int>(frames.size());
        m_image_geoms = frames;
        m_images.assign(m_nimages, nullptr);
        m_image_done.assign(m_nimages, false);
        resetPmfAccumulators(static_cast<int>(frames.front().rows()));
        m_start_geom = frames.front();
        m_end_geom = frames.back();
        std::vector<int> atoms = m_start.Atoms();
        m_masses.assign(atoms.size(), 0.0);
        for (size_t a = 0; a < atoms.size(); ++a)
            m_masses[a] = Elements::AtomicMass[atoms[a]];
        CurcumaLogger::result_fmt("NEB-MD: reference path read from {} ({} images)",
                                  m_path_file, m_nimages);
        return true;
    }

    // Align the end onto the start frame (Kabsch) so linear interpolation gives a
    // minimal-rotation path. Atom-wise interpolation REQUIRES both endpoints to
    // share the same atom order, so:
    //  - if the start/end atom sequences already match -> Kabsch align only
    //    (no_reorder, preserves the order);
    //  - if they differ and force_reorder is set -> attempt a Munkress reorder
    //    (method "inertia") and verify the result matches the start order;
    //  - otherwise -> error (the user must pre-align atom orders, e.g. -nebprep).
    // A dedicated rmsd config avoids the flat global "-method" (energy method)
    // leaking into the rmsd module's "method" parameter.
    std::vector<int> start_atoms = m_start.Atoms();
    std::vector<int> end_atoms = m_end.Atoms();
    const bool same_order = (start_atoms == end_atoms);

    json rmsd_mod = json::object({{"protons", true}, {"threads", 1}});
    if (same_order || !m_force_reorder) {
        rmsd_mod["no_reorder"] = true; // Kabsch align only, preserve atom order
    } else {
        rmsd_mod["method"] = "inertia"; // Munkress reorder + Kabsch
    }
    json rmsd_cfg = json::object();
    rmsd_cfg["rmsd"] = rmsd_mod;

    RMSDDriver driver(rmsd_cfg, true);
    driver.setReference(m_start);
    driver.setTarget(m_end);
    driver.setCheckConnections(true);
    driver.setProtonTransfer(m_pt);

    {
        CoutSuppress cs; // RMSD alignment is chatty
        driver.start();
    }

    Molecule end_aligned = (same_order || !m_force_reorder)
                                ? driver.TargetAligned()
                                : driver.TargetReorderd();
    if (end_aligned.AtomCount() == 0)
        end_aligned = driver.TargetAligned(); // reorder produced nothing -> fall back
    if (end_aligned.AtomCount() == 0) {
        CurcumaLogger::error("NEB-MD: RMSD alignment of the endpoints failed.");
        return false;
    }
    if (end_aligned.Atoms() != start_atoms) {
        CurcumaLogger::error("NEB-MD: end atom order does not match the start after "
                           "alignment. Pre-align the endpoints (e.g. 'curcuma -nebprep "
                           "start.xyz end.xyz') so the atom orders are consistent.");
        return false;
    }

    // RMSDDriver works in a centred frame; shift the aligned end back so that
    // image 0 equals the original start coordinates exactly.
    Geometry start_orig = m_start.getGeometry();
    Geometry start_centred = driver.ReferenceAligned().getGeometry();
    if (start_centred.rows() != start_orig.rows() || start_centred.rows() == 0) {
        CurcumaLogger::error("NEB-MD: could not obtain the aligned reference geometry.");
        return false;
    }
    Geometry shift = start_orig - start_centred; // constant per atom (the centroid)
    m_start_geom = start_orig;
    m_end_geom = end_aligned.getGeometry() + shift;

    // Per-atom masses (amu) for mass-weighted springs.
    std::vector<int> atoms = m_start.Atoms();
    m_masses.assign(atoms.size(), 0.0);
    for (size_t a = 0; a < atoms.size(); ++a)
        m_masses[a] = Elements::AtomicMass[atoms[a]];

    for (int i = 0; i < m_nimages; ++i) {
        double t = (m_nimages == 1) ? 0.0 : static_cast<double>(i) / static_cast<double>(m_nimages - 1);
        m_image_geoms[i] = interpolate(m_start_geom, m_end_geom, t);
    }
    return true;
}

void NebMD::checkExtensionEnergies(const std::vector<Geometry>& base,
                                   const std::vector<Geometry>& extended) const
{
    // Claude Generated 2026: an extrapolated frame is only meaningful as the continuation of a
    // basin if the energy rises beyond the endpoint. A LOWER energy means the path ends on a
    // slope (the endpoint is not a minimum along the continued direction) - warn, do not abort.
    const int m = m_cv_extend;
    const int N = static_cast<int>(base.size());
    if (m <= 0 || N < 2 || static_cast<int>(extended.size()) != N + 2 * m)
        return;
    const double kcal = 627.5094740631;
    auto sp = [this](const Geometry& g) -> double {
        EnergyCalculator ec(m_method, m_controller, Basename());
        Molecule mol = m_start;
        mol.setGeometry(g);
        ec.setMolecule(mol.getMolInfo());
        ec.updateGeometry(g);
        return ec.CalculateEnergy(false);
    };
    double e0 = 0.0, eN = 0.0;
    std::vector<double> eL(m), eR(m); // eL[j-1] = E(F_-j), eR[j-1] = E(F_N-1+j)
    {
        CoutSuppress cs;
        e0 = sp(base.front());
        eN = sp(base.back());
        for (int j = 1; j <= m; ++j) {
            eL[j - 1] = sp(extended[m - j]);
            eR[j - 1] = sp(extended[m + N - 1 + j]);
        }
    }
    CurcumaLogger::set_verbosity(m_band_verbosity);
    std::string l = "NEB-MD cv_extend: E(extension) - E(endpoint) [kcal/mol], start side:";
    std::string r = "NEB-MD cv_extend: E(extension) - E(endpoint) [kcal/mol], end side  :";
    bool lower = false;
    for (int j = 0; j < m; ++j) {
        l += fmt::format(" {:.3f}", (eL[j] - e0) * kcal);
        r += fmt::format(" {:.3f}", (eR[j] - eN) * kcal);
        if (eL[j] < e0 || eR[j] < eN)
            lower = true;
    }
    CurcumaLogger::result(l);
    CurcumaLogger::result(r);
    if (lower)
        CurcumaLogger::warn("NEB-MD cv_extend: an extrapolated CV frame lies LOWER in energy than the path "
                            "endpoint it continues. The endpoint is then no basin along the continued "
                            "direction; the extension is physically doubtful (check the endpoints are "
                            "minima, or use a smaller cv_extend).");
}

std::unique_ptr<PathCV> NebMD::makePathCV(const std::vector<Geometry>& frames) const
{
    // Claude Generated 2026: the one place that turns the CV settings into a PathCV.
    if (m_cv_metric == PathCV::Metric::Distances)
        return std::make_unique<PathCV>(frames, m_start, -1.0, std::vector<int>{}, PathCV::Metric::Distances,
                                        m_cv_pairs, m_cv_transform);
    return std::make_unique<PathCV>(frames, m_start, -1.0, m_cv_atoms);
}

double NebMD::cvSpacing(const Geometry& A, const Geometry& B) const
{
    // Claude Generated 2026: plain (not best-fit) spacing of two frames in the CV metric.
    if (m_cv_metric == PathCV::Metric::Distances)
        return std::sqrt(PathCV::pairSpacing2(A, B, m_cv_pairs, m_cv_transform));
    return curcuma::string_path::weightedRmsd(A, B, m_cv_w, false);
}

std::string NebMD::cvPairsText() const
{
    // Claude Generated 2026
    std::string t = fmt::format("x = {}, {} pairs:", m_cv_transform == PathCV::PairTransform::Inverse ? "1/d" : "d",
                                m_cv_pairs.size());
    for (const auto& p : m_cv_pairs)
        t += fmt::format(" {}-{}", p.first + 1, p.second + 1);
    return t;
}

bool NebMD::resolveCvPairs(const std::vector<Geometry>& frames)
{
    // Claude Generated 2026: (1) the pairs of the distance metric for cv_pairs 'auto', (2) the
    // pairs of the window diagnostics (J): the CV pairs for the distance metric, the heavy-atom
    // 'auto' pairs for the rmsd metric (no error if there are fewer than 2: J is then n/a).
    const int nat = static_cast<int>(m_start.AtomCount());
    std::vector<int> pool, heavy;
    if (m_cv_atoms.empty())
        for (int a = 0; a < nat; ++a)
            pool.push_back(a);
    else
        pool = m_cv_atoms;
    for (int a : pool)
        if (m_start.Atom(a).first > 1)
            heavy.push_back(a);
    auto listText = [](const std::vector<PathCV::AtomPair>& pairs, const std::vector<double>& spans) {
        std::string t;
        for (size_t p = 0; p < pairs.size(); ++p)
            t += fmt::format("{}{}-{} ({:.2f})", p == 0 ? "" : ", ", pairs[p].first + 1, pairs[p].second + 1,
                             p < spans.size() ? spans[p] : 0.0);
        return t;
    };
    std::vector<double> spans;
    if (m_cv_metric == PathCV::Metric::Distances) {
        if (!m_cv_pairs_resolved) {
            const std::vector<int>& cand = m_cv_atoms.empty() ? heavy : m_cv_atoms;
            m_cv_pairs = PathCV::selectPairsBySpan(frames, cand, m_cv_pair_threshold, &spans);
            if (m_cv_pairs.size() < 2) {
                CurcumaLogger::error(fmt::format(
                    "NEB-MD: cv_pairs auto found {} pair(s) whose distance spans more than cv_pair_threshold = {} A "
                    "over the CV frames ({} candidate atoms, {}); the distance-space CV needs at least 2. Lower "
                    "cv_pair_threshold, widen cv_atoms or give cv_pairs explicitly.",
                    m_cv_pairs.size(), m_cv_pair_threshold, cand.size(),
                    m_cv_atoms.empty() ? "heavy atoms" : "cv_atoms"));
                return false;
            }
            m_cv_pairs_resolved = true;
            CurcumaLogger::result_fmt("NEB-MD: cv_pairs auto ({} candidate atoms, span > {} A over {} CV frames): {} pairs, "
                                      "i-j (span A): {}", cand.size(), m_cv_pair_threshold, frames.size(),
                                      m_cv_pairs.size(), listText(m_cv_pairs, spans));
        } else {
            for (const auto& p : m_cv_pairs) {
                double dmin = 1e300, dmax = -1e300;
                for (const auto& g : frames) {
                    const double d = (g.row(p.first) - g.row(p.second)).norm();
                    dmin = std::min(dmin, d);
                    dmax = std::max(dmax, d);
                }
                spans.push_back(dmax - dmin);
            }
            CurcumaLogger::result_fmt("NEB-MD: cv_pairs explicit: {} pairs, i-j (span A over {} CV frames): {}",
                                      m_cv_pairs.size(), frames.size(), listText(m_cv_pairs, spans));
        }
        m_diag_pairs = m_cv_pairs;
        m_diag_transform = m_cv_transform;
    } else {
        m_diag_pairs = PathCV::selectPairsBySpan(frames, heavy, m_cv_pair_threshold, &spans);
        m_diag_transform = PathCV::PairTransform::Distance;
        if (m_diag_pairs.size() >= 2)
            CurcumaLogger::info(fmt::format("NEB-MD: window diagnostics: {} heavy-atom pairs (span > {} A): {}",
                                            m_diag_pairs.size(), m_cv_pair_threshold, listText(m_diag_pairs, spans)));
        else
            CurcumaLogger::info("NEB-MD: window diagnostics: fewer than 2 heavy-atom pairs change along the path; "
                                "the neighbour-jump index J is not available.");
    }
    return true;
}

void NebMD::accumulateWindowDiagnostics(int i, const Geometry& R)
{
    // Claude Generated 2026: called every kDiagInterval production steps from the parallel window
    // loop; touches only the entries of window i (pre-sized in setupUmbrella), so it is race-free.
    // No RNG, no change of positions: the dynamics is unaffected.
    if (!m_path_cv || i < 0 || i >= static_cast<int>(m_diag_n.size()))
        return;
    double z = 0.0;
    const int nn = m_path_cv->nearestFrame(R, &z);
    if (!std::isfinite(z))
        return;
    m_diag_z[i].push_back(z);
    m_diag_nn_sum[i] += nn;
    m_diag_nn2_sum[i] += static_cast<double>(nn) * nn;
    if (m_diag_pairs.size() >= 2)
        m_diag_x_sum[i] += PathCV::pairValues(R, m_diag_pairs, m_diag_transform);
    ++m_diag_n[i];
}

void NebMD::writeWindowDiagnostics() const
{
    // Claude Generated 2026: <base>.neb.windows.csv, one row per window.
    //   z_mean, z_p90, frac_beyond_wall: from z sampled every kDiagInterval steps (frac = nan without wall)
    //   nn_mean, nn_sd: index of the nearest CV frame (argmin_i d_i), 0-based on the CV frames
    //   J (row i, pair i,i+1): max_p |<x_p>_i - <x_p>_i+1| / max_k max_p |x_p,k+1 - x_p,k|
    // J > 3: the time-averaged pair coordinates of neighbouring windows differ by more than three
    // of the largest single CV-frame steps - the windows sample different structural families.
    CurcumaLogger::set_verbosity(m_band_verbosity);
    const int nw = static_cast<int>(m_diag_n.size());
    if (nw == 0)
        return;
    const bool have_j = m_diag_pairs.size() >= 2 && m_diag_x_step > 0.0;
    std::vector<double> J(nw, std::numeric_limits<double>::quiet_NaN());
    if (have_j)
        for (int i = 0; i + 1 < nw; ++i) {
            if (m_diag_n[i] < 1 || m_diag_n[i + 1] < 1)
                continue;
            const Vector a = m_diag_x_sum[i] / static_cast<double>(m_diag_n[i]);
            const Vector b = m_diag_x_sum[i + 1] / static_cast<double>(m_diag_n[i + 1]);
            J[i] = (a - b).cwiseAbs().maxCoeff() / m_diag_x_step;
        }
    const std::string path = outputPath(Basename() + ".neb.windows.csv");
    std::ofstream f(path);
    if (!f) {
        CurcumaLogger::warn("NEB-MD: cannot write " + path);
        return;
    }
    const std::string xunit = m_diag_transform == PathCV::PairTransform::Inverse ? "1/A" : "A";
    f << "# window diagnostics (every " << kDiagInterval << " production steps); cv_metric "
      << (m_cv_metric == PathCV::Metric::Distances ? "distances" : "rmsd") << ", z in "
      << (m_path_cv ? m_path_cv->zUnit() : std::string("A^2")) << ", nn = nearest CV frame (0-based, "
      << (m_path_cv ? m_path_cv->frames() : 0) << " frames)\n";
    if (have_j) {
        f << "# J pairs (" << (m_cv_metric == PathCV::Metric::Distances ? "CV pairs" : "heavy-atom auto pairs")
          << ", x = " << (m_diag_transform == PathCV::PairTransform::Inverse ? "1/d" : "d") << "):";
        for (const auto& p : m_diag_pairs)
            f << " " << p.first + 1 << "-" << p.second + 1;
        f << fmt::format("; reference step max_k max_p |x_p,k+1 - x_p,k| = {:.6f} {}\n", m_diag_x_step, xunit);
    } else {
        f << "# J not available (fewer than 2 diagnostic pairs)\n";
    }
    f << "window,s0,z_mean,z_p90,frac_beyond_wall,nn_mean,nn_sd,J\n";
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    double jmax = -1.0;
    int jmax_i = -1;
    for (int i = 0; i < nw; ++i) {
        const long n = m_diag_n[i];
        double zm = nan, zp = nan, fr = nan, nm = nan, nsd = nan;
        if (n > 0) {
            std::vector<double> zs = m_diag_z[i];
            std::sort(zs.begin(), zs.end());
            double zsum = 0.0;
            long nb = 0;
            for (double z : zs) {
                zsum += z;
                if (m_z_wall > 0.0 && z > m_z_wall)
                    ++nb;
            }
            zm = zsum / static_cast<double>(zs.size());
            const size_t r = static_cast<size_t>(std::ceil(0.9 * static_cast<double>(zs.size())));
            zp = zs[std::min(std::max<size_t>(r, 1), zs.size()) - 1]; // nearest-rank p90
            if (m_z_wall > 0.0)
                fr = static_cast<double>(nb) / static_cast<double>(zs.size());
            nm = m_diag_nn_sum[i] / static_cast<double>(n);
            nsd = std::sqrt(std::max(m_diag_nn2_sum[i] / static_cast<double>(n) - nm * nm, 0.0));
        }
        const double s0 = i < static_cast<int>(m_umbrella_s0.size()) ? m_umbrella_s0[i] : nan;
        f << fmt::format("{},{:.6f},{:.6f},{:.6f},{:.6f},{:.4f},{:.4f},{:.4f}\n", i, s0, zm, zp, fr, nm, nsd, J[i]);
        if (std::isfinite(J[i]) && J[i] > jmax) {
            jmax = J[i];
            jmax_i = i;
        }
        if (std::isfinite(J[i]) && J[i] > 3.0)
            CurcumaLogger::warn(fmt::format("NEB-MD: structural family changes between window {} and {} (J = {:.2f} > 3: "
                                            "<x_p> jumps by {:.3f} {} vs a CV-frame step of {:.3f} {}).",
                                            i, i + 1, J[i], J[i] * m_diag_x_step, xunit, m_diag_x_step, xunit));
    }
    if (jmax_i >= 0)
        CurcumaLogger::result_fmt("NEB-MD: window diagnostics written to {} (largest J = {:.2f} between windows {} and {})",
                                  path, jmax, jmax_i, jmax_i + 1);
    else
        CurcumaLogger::result_fmt("NEB-MD: window diagnostics written to {} (J not available)", path);
}

bool NebMD::buildCvPath(bool place_windows)
{
    // Claude Generated 2026: CV frames and window placement (see the block comment in nebmd.h).
    m_cv_path.clear();
    m_cv_orig_s.clear();
    m_cv_lo = 0;
    m_cv_hi = -1;
    // Claude Generated 2026: the CV pairs ('auto') and the diagnostic pairs are resolved on the
    // frames the CV is built on; only runs that build a path CV need them.
    const bool need_pairs = m_umbrella || !m_cv_histogram.empty();
    if (m_cv_frames <= 0 && m_cv_extend <= 0 && m_cv_extend_length <= 0.0)
        return need_pairs ? resolveCvPairs(m_image_geoms) : true; // default: CV on the image frames
    if (place_windows && !m_independent) {
        if (need_pairs && !resolveCvPairs(m_image_geoms))
            return false;
        if (m_string)
            CurcumaLogger::info("NEB-MD: cv_frames/cv_extend do not act on the string stage; they apply to "
                                "the path CV of scheme 'independent' (pipeline stage 2).");
        else
            CurcumaLogger::warn("NEB-MD: cv_frames/cv_extend are only used by scheme 'independent' and by "
                                "cv_histogram - ignored here.");
        return true;
    }
    if (m_image_geoms.size() < 2) {
        CurcumaLogger::error("NEB-MD: cv_frames/cv_extend need a path with at least two frames.");
        return false;
    }
    const int nat = static_cast<int>(m_image_geoms.front().rows());
    const std::vector<Geometry>& P = m_image_geoms;
    const int Np = static_cast<int>(P.size());

    // cv_atoms: arclength (re-sampling, L) in the selected-atom metric m_cv_w (0/1, empty = all).
    // The Kabsch fits (consecutive alignment, extension) use m_cv_w_align = 1 on the selected
    // atoms and 1e-6 on the others: the selected atoms decide the fit, the others only fix a
    // rotation the selection leaves free. Measured on the SN2 path (Cl-C-Cl, exactly collinear):
    // with pure 0/1 weights the H atoms jumped by up to 2.1 A between consecutive frames and the
    // extension frames lay up to 8165 kcal/mol above the endpoint; with the tie-break as with 'all'.
    std::vector<Geometry> base = curcuma::string_path::alignConsecutive(P, m_cv_w_align);
    if (m_cv_frames > 0) {
        // Part B: re-sample at equal arclength onto cv_frames frames (before the extension).
        int nf = m_cv_frames;
        if (nf < 3) {
            CurcumaLogger::warn(fmt::format("NEB-MD: cv_frames = {} is too few for a path CV; using 3.", nf));
            nf = 3;
        }
        const std::vector<Geometry> rs = curcuma::string_path::resamplePath(base, nf, m_cv_w);
        // Claude Generated 2026: 'auto' pairs from the RMSD-re-sampled frames (spec of cv_pairs),
        // then, for the distance metric, the final re-sampling at equal distance-space arclength.
        if (!resolveCvPairs(rs))
            return false;
        if (m_cv_metric == PathCV::Metric::Distances)
            base = PathCV::resampleByPairs(base, nf, m_cv_pairs, m_cv_transform);
        else
            base = rs;
    } else if (!resolveCvPairs(base)) {
        return false;
    }
    if (m_cv_extend_length > 0.0) {
        // Claude Generated 2026: extension given as a LENGTH (fraction of the arclength L of the
        // (re-sampled) CV path). m = ceil(f L / d) + 1 with d = L/(F-1); m_cv_extend then holds
        // the effective frame count for everything downstream (extension, windows, headers).
        const int F = static_cast<int>(base.size());
        double L = 0.0; // plain (already consecutively aligned) RMSD arclength, as PathCV sees it
        for (int k = 1; k < F; ++k)
            L += (m_cv_metric == PathCV::Metric::Distances)
                ? cvSpacing(base[k], base[k - 1]) // distance-space arclength (Claude Generated 2026)
                : curcuma::string_path::weightedRmsd(base[k], base[k - 1], m_cv_w, false);
        m_cv_extend = curcuma::string_path::extensionFrames(m_cv_extend_length, F);
        const double d = F > 1 ? L / static_cast<double>(F - 1) : 0.0;
        CurcumaLogger::result_fmt("NEB-MD cv_extend: {:.4g} of the path length = {} frames "
                                  "(path L = {:.4f} {} arclength, {} CV frames, spacing d = {:.4f}; "
                                  "windows reach {:.4f} = {:.4f} L beyond each end)",
                                  m_cv_extend_length, m_cv_extend, L,
                                  m_cv_metric == PathCV::Metric::Distances ? "distance-space" : "A RMSD", F, d,
                                  (m_cv_extend - 1) * d, L > 0.0 ? (m_cv_extend - 1) * d / L : 0.0);
    }
    std::vector<Geometry> cv = base;
    if (m_cv_extend > 0) {
        cv = curcuma::string_path::extendPath(base, m_cv_extend, m_cv_w_align);
        checkExtensionEnergies(base, cv);
    }
    m_cv_path = cv;
    const int M = static_cast<int>(cv.size());
    m_cv_lo = (m_cv_extend > 0) ? 1 : 0;
    m_cv_hi = M - 1 - m_cv_lo;

    const std::unique_ptr<PathCV> pcvp = makePathCV(m_cv_path); // Claude Generated 2026: active metric
    const PathCV& pcv = *pcvp;
    std::vector<double> sc(M, 0.0);
    for (int k = 0; k < M; ++k)
        sc[k] = pcv.s(cv[k]);
    m_cv_orig_s = { pcv.s(P[0]), pcv.s(P[std::min(1, Np - 1)]), pcv.s(P[std::max(Np - 2, 0)]), pcv.s(P[Np - 1]) };
    bool monotone = true;
    for (int k = m_cv_lo + 1; k <= m_cv_hi; ++k)
        if ((sc[k] - sc[k - 1]) * (sc[m_cv_hi] - sc[m_cv_lo]) <= 0.0)
            monotone = false;
    CurcumaLogger::result_fmt("NEB-MD: path CV on {} frames (path {} frames, cv_frames {}, cv_extend {}), "
                              "lambda = {:.3f}; s(CV frame 0) = {:.4f}, window range s = [{:.4f}, {:.4f}], "
                              "s(CV frame {}) = {:.4f}",
                              M, Np, m_cv_frames, m_cv_extend, pcv.lambda(), sc.front(), sc[m_cv_lo],
                              sc[m_cv_hi], M - 1, sc.back());
    if (m_cv_metric == PathCV::Metric::Distances) {
        // Claude Generated 2026: spacing and length of the CV frames in the distance metric.
        double L = 0.0, dmin = 1e300, dmax = 0.0;
        for (int k = m_cv_lo + 1; k <= m_cv_hi; ++k) {
            const double d = cvSpacing(cv[k], cv[k - 1]);
            L += d;
            dmin = std::min(dmin, d);
            dmax = std::max(dmax, d);
        }
        const double dm = L / std::max(m_cv_hi - m_cv_lo, 1);
        CurcumaLogger::result_fmt("NEB-MD: distance-space CV ({} pairs): spacing mean {:.4f} (min {:.4f}, max {:.4f}) {}, "
                                  "path length L = {:.4f} {} over the window range, lambda*d^2 = {:.3f}",
                                  m_cv_pairs.size(), dm, dmin, dmax,
                                  m_cv_transform == PathCV::PairTransform::Inverse ? "1/A" : "A", L,
                                  m_cv_transform == PathCV::PairTransform::Inverse ? "1/A" : "A",
                                  pcv.lambda() * dm * dm);
    }
    CurcumaLogger::result_fmt("NEB-MD: original path frames on this CV: s(first) = {:.4f}, s(second) = {:.4f}, "
                              "s(second-to-last) = {:.4f}, s(last) = {:.4f}",
                              m_cv_orig_s[0], m_cv_orig_s[1], m_cv_orig_s[2], m_cv_orig_s[3]);
    if (!monotone)
        CurcumaLogger::warn("NEB-MD: s is not monotonic over the window range of the CV frames - the path "
                            "folds back on itself in RMSD space; window placement is unreliable.");
    if (!place_windows)
        return true;

    std::vector<Geometry> win;
    if (m_cv_frames > 0) {
        // Part B: nwin = nimages windows, setpoints equidistant in s over [s(C_lo), s(C_hi)];
        // start geometry = point of the CV polyline at that s (bracketing frames k, k+1 and
        // t = (s* - s_k)/(s_k+1 - s_k)). setupUmbrella re-evaluates s0 on that geometry.
        const int nwin = std::max(m_nimages, 2);
        const double a = sc[m_cv_lo], b = sc[m_cv_hi];
        for (int i = 0; i < nwin; ++i) {
            if (i == 0) {
                win.push_back(cv[m_cv_lo]);
                continue;
            }
            if (i == nwin - 1) {
                win.push_back(cv[m_cv_hi]);
                continue;
            }
            const double target = a + (b - a) * static_cast<double>(i) / static_cast<double>(nwin - 1);
            int k = -1;
            for (int j = m_cv_lo; j < m_cv_hi; ++j)
                if ((sc[j] - target) * (sc[j + 1] - target) <= 0.0) {
                    k = j;
                    break;
                }
            if (k < 0) { // non-monotonic s: nearest CV frame
                k = m_cv_lo;
                for (int j = m_cv_lo; j <= m_cv_hi; ++j)
                    if (std::abs(sc[j] - target) < std::abs(sc[k] - target))
                        k = j;
                win.push_back(cv[k]);
                continue;
            }
            const double den = sc[k + 1] - sc[k];
            const double t = std::abs(den) > 1e-15 ? (target - sc[k]) / den : 0.0;
            win.push_back((1.0 - t) * cv[k] + t * cv[k + 1]);
        }
    } else {
        win.assign(cv.begin() + m_cv_lo, cv.begin() + m_cv_hi + 1);
    }
    m_nimages = static_cast<int>(win.size());
    m_image_geoms = win;
    m_images.assign(m_nimages, nullptr);
    m_image_done.assign(m_nimages, false);
    resetPmfAccumulators(nat);
    m_start_geom = win.front();
    m_end_geom = win.back();
    if (m_cv_extend > 0 && m_fixed_endpoints) {
        m_fixed_endpoints = false;
        CurcumaLogger::result("NEB-MD: cv_extend > 0 - every window is integrated (fixed_endpoints off); the "
                              "former path endpoints are ordinary windows now.");
    }
    CurcumaLogger::result_fmt("NEB-MD: {} windows on the CV, s from {:.4f} to {:.4f}", m_nimages,
                              sc[m_cv_lo], sc[m_cv_hi]);
    return true;
}

bool NebMD::prepareImages()
{
    int nat = static_cast<int>(m_start_geom.rows());
    if (nat == 0) {
        CurcumaLogger::error("NEB-MD: no path geometry - call buildInitialPath() first.");
        return false;
    }

    // Base simplemd sub-config shared by every image.
    json sd = m_md_config;
    sd["no_center"] = true;            // CRITICAL: do not recenter each image (would break inter-image geometry)
    sd["remove_com_mode"] = 0;         // no COM/rotation removal (keeps the shared frame)
    sd["no_restart"] = true;           // never load stale per-image restart snapshots (NEB-MD owns state)
    sd["write_xyz"] = m_window_trajectory; // per-image trajectories only on request (-nebmd.window_trajectory)
    sd["write_initial_state"] = false;
    sd["write_restart_frequency"] = -1; // disable per-image restart files (step() gates on `> -1` then `% freq`; 0 would SIGFPE)
    sd["opt"] = false;
    sd["rmsd_mtd"] = false;
    if (m_constrain) {
        // The scalar constraint lives inside SimpleMD::Rattle(), so that integrator has to
        // be selected. Bond constraints are switched off explicitly: they would change the
        // degrees of freedom and thereby the very ensemble whose free energy is measured.
        sd["rattle"] = 1;
        sd["rattle_12"] = false;
        sd["rattle_13"] = false;
    }
    // Claude Generated 2026: the string stage owns its run length (string_max_time); the
    // band loop ends earlier on convergence.
    if (m_string)
        sd["max_time"] = m_string_max_time;
    sd["mtd"] = false;
    // Claude Generated (Oct 2026): print_frequency below suppresses the per-image status lines,
    // so every image would only leave a header-only <image>.md.csv in the BMT dir.
    sd["write_csv"] = false;
    sd["wall_type"] = "none";
    sd["unique"] = false;
    sd["threads"] = m_threads;
    sd["method"] = m_method;
    if (!sd.contains("dump_frequency")) sd["dump_frequency"] = m_dump;
    // print_frequency must be > 0: step() does `int(step*dt) % print_frequency`
    // and %0 would SIGFPE. Use a large value so per-image status prints are suppressed.
    if (!sd.contains("print_frequency")) sd["print_frequency"] = 1000000;
    if (!sd.contains("verbosity")) sd["verbosity"] = 0;

    for (int i = 0; i < m_nimages; ++i) {
        const bool is_endpoint = (i == 0 || i == m_nimages - 1);
        if (m_fixed_endpoints && is_endpoint) {
            m_images[i] = nullptr; // clamped; positions stay in m_image_geoms[i]
            continue;
        }

        json img_ctrl = m_controller;
        // Each image needs its OWN random stream. Claude Generated 2026: handing every
        // image the same seed would give all windows identical thermostat noise, so the
        // replicas would be correlated by construction - the opposite of what an
        // independent-window scheme requires. Deriving the seed from the image index keeps
        // the run reproducible and makes the serial and the parallel result identical.
        {
            json sd_i = sd;
            const int base = sd.contains("seed") && sd["seed"].is_number_integer()
                ? sd["seed"].get<int>() : 0;
            sd_i["seed"] = (base > 0 ? base : 20260802) + 7919 * (i + 1);
            img_ctrl["simplemd"] = sd_i;
        }
        // Distinct seed per image so velocities are independent.
        int base_seed = sd.value("seed", -1);
        img_ctrl["simplemd"]["seed"] = (base_seed <= 0) ? (1000 + i * 101) : (base_seed + i);

        Molecule mol = m_start; // copies atoms / charge / name
        mol.setGeometry(m_image_geoms[i]);

        std::string imgbase = Basename() + ".image_"
                            + (i < 10 ? "0" : "") + std::to_string(i);

        SimpleMD* md = nullptr;
        bool ok = false;
        {
            CoutSuppress cs; // suppress per-image constructor/Initialise/prepareRun chatter
            md = new SimpleMD(img_ctrl, true);
            md->setRestart(false); // ensure InitVelocities runs and no restart is loaded
            md->overrideBasename(imgbase);
            if (!OutputDir().empty())
                md->setOutputDir(OutputDir()); // route image output into the NEB-MD BMT dir
            md->setMolecule(mol);
            ok = md->Initialise();
            if (ok) md->prepareRun();
        }
        if (!ok) {
            CurcumaLogger::error("NEB-MD: image " + std::to_string(i) + " initialisation failed.");
            delete md;
            m_images[i] = nullptr;
            return false;
        }
        m_images[i] = md;
    }

    computeEndpointEnergies();
    // Re-assert our verbosity: each image's EnergyCalculator setup clamps the global
    // CurcumaLogger verbosity to 0 (see simplemd.cpp / capabilities CLAUDE.md #1), which
    // would otherwise suppress the NEB-MD result/summary lines below.
    CurcumaLogger::set_verbosity(getVerbosity());
    return true;
}

void NebMD::computeEndpointEnergies()
{
    if (!m_fixed_endpoints)
        return; // endpoints are full images; their energy comes from potentialEnergy()

    auto sp = [this](const Geometry& g) -> double {
        EnergyCalculator ec(m_method, m_controller, Basename());
        Molecule mol = m_start;
        mol.setGeometry(g);
        ec.setMolecule(mol.getMolInfo());
        ec.updateGeometry(g);
        return ec.CalculateEnergy(false);
    };
    m_endpoint_energy_start = sp(m_start_geom);
    m_endpoint_energy_end = sp(m_end_geom);
}

std::vector<Geometry> NebMD::currentPath() const
{
    std::vector<Geometry> R(m_nimages);
    for (int i = 0; i < m_nimages; ++i)
        R[i] = m_images[i] ? m_images[i]->positions() : m_image_geoms[i];
    return R;
}

std::vector<double> NebMD::currentEnergies() const
{
    std::vector<double> E(m_nimages, 0.0);
    for (int i = 0; i < m_nimages; ++i) {
        if (m_images[i])
            E[i] = m_images[i]->potentialEnergy();
        else
            E[i] = (i == 0) ? m_endpoint_energy_start : m_endpoint_energy_end;
    }
    return E;
}

// Improved tangent (Henkelman-Jónsson 2000): energy-weighted mix of the forward
// and backward unit segments, returned as a unit vector in 3N space per image.
std::vector<Geometry> NebMD::computeTangents(const std::vector<Geometry>& R,
                                             const std::vector<double>& E) const
{
    int nat = static_cast<int>(R[0].rows());
    std::vector<Geometry> tau(m_nimages, Geometry::Zero(nat, 3));

    for (int i = 0; i < m_nimages; ++i) {
        Geometry t = Geometry::Zero(nat, 3);
        if (i == 0) {
            t = R[1] - R[0];
        } else if (i == m_nimages - 1) {
            t = R[m_nimages - 1] - R[m_nimages - 2];
        } else {
            Geometry tp = R[i + 1] - R[i];
            Geometry tm = R[i] - R[i - 1];
            double np = tp.norm();
            double nm = tm.norm();
            Geometry tp_n = Geometry::Zero(nat, 3);
            Geometry tm_n = Geometry::Zero(nat, 3);
            if (np > 1e-12) tp_n = tp / np;
            if (nm > 1e-12) tm_n = tm / nm;
            double dp = std::max(E[i + 1] - E[i], 0.0);
            double dm = std::max(E[i - 1] - E[i], 0.0);
            t = dp * tp_n + dm * tm_n;
        }
        double norm = t.norm();
        if (norm > 1e-12)
            t /= norm;
        else
            t.setZero();
        tau[i] = t;
    }
    return tau;
}

std::vector<Geometry> NebMD::computeSpringForces(const std::vector<Geometry>& R,
                                                 const std::vector<Geometry>& tangents) const
{
    int nat = static_cast<int>(R[0].rows());
    std::vector<Geometry> F(m_nimages, Geometry::Zero(nat, 3));

    for (int i = 1; i < m_nimages - 1; ++i) {
        Geometry fs = Geometry::Zero(nat, 3);
        for (int a = 0; a < nat; ++a) {
            double k_a = m_mass_weighted_spring ? m_k_spring * m_masses[a] : m_k_spring;
            fs.row(a) = k_a * (R[i + 1].row(a) + R[i - 1].row(a) - 2.0 * R[i].row(a));
        }
        if (m_nudged) {
            // Keep only the component parallel to the tangent (NEB nudging).
            double dot = 0.0;
            for (int a = 0; a < nat; ++a)
                dot += fs.row(a).dot(tangents[i].row(a));
            fs = dot * tangents[i];
        }
        F[i] = fs;
    }
    return F;
}

bool NebMD::stepBand()
{
    std::vector<Geometry> R = currentPath();
    std::vector<double> E = currentEnergies();
    std::vector<Geometry> tau = computeTangents(R, E);
    std::vector<Geometry> Fspring = computeSpringForces(R, tau);
    int nat = static_cast<int>(R[0].rows());

    // 1) Compute + inject the external (spring / nudging) force for every live image.
    //
    // Claude Generated 2026: the loop used to run over interior images only, so with
    // fixed_endpoints=false the endpoint images were created but never driven or advanced
    // - opt-in behaviour that silently did nothing. They are now included.
    //
    // This is NOT a change of the default: with fixed_endpoints=true the endpoints are
    // nullptr and the guard below skips them, so the default path is bit-identical. The
    // opt-in is the existing -nebmd.fixed_endpoints flag. It matters for short paths,
    // where a usable profile needs the endpoints sampled as well.
    for (int i = 0; i < m_nimages; ++i) {
        if (m_image_done[i])
            continue;
        SimpleMD* img = m_images[i];
        if (!img)
            continue;

        Geometry g = img->trueGradient(); // true gradient dE/dR at the current position (Eh/Å)
        Geometry Fext = Geometry::Zero(nat, 3);
        // applyExternalForces() adds its argument to m_eigen_gradient (dE/dR); the integrator
        // then applies force = -gradient. So to apply a desired FORCE F_desired we must pass
        // -F_desired (the gradient contribution). SimpleMD already applies F_true = -g.
        if (m_nudged) {
            // Desired total force = F_true_perp + F_spring_par = -g + (g·τ)τ + F_spring_par.
            // SimpleMD applies -g, so the extra force needed is (g·τ)τ + F_spring_par, i.e. pass
            // the negative of that as the gradient contribution.
            double gdot = 0.0;
            for (int a = 0; a < nat; ++a)
                gdot += g.row(a).dot(tau[i].row(a));
            Fext = -(gdot * tau[i] + Fspring[i]);
        } else {
            // Plain elastic band: extra force = F_spring, pass -F_spring.
            Fext = -Fspring[i];
        }
        // Umbrella mode: restrain s (progress along the path) and optionally z (distance
        // FROM the path). Whether the neighbour springs are kept is a choice:
        //   umbrella_keep_springs = false -> pure umbrella, windows are independent
        //   umbrella_keep_springs = true  -> images additionally feel their neighbours,
        //                                    which is what kept them from wandering off
        //                                    in the plain NEB-MD mode.
        // Keeping the springs biases the sampling (the bias is not removed by WHAM), so
        // it trades statistical rigour for staying on the path - report both.
        // Claude Generated 2026: scheme 'string' - no springs, no nudging; only the restraint
        // to the frame M_i of the (between updates frozen) mean path.
        if (m_string) {
            Fext.setZero();
            applyStringRestraint(i, img, Fext);
        } else if (m_umbrella) {
            if (!m_umbrella_keep_springs)
                Fext.setZero();
            applyUmbrellaBias(i, img, Fext);
        }

        img->applyExternalForces(Fext); // assignment; consumed by the next step()

        // Transverse metadynamics: hand the current tangent to the image so its RMSD-MTD
        // bias (if enabled via -simplemd.rmsd_mtd) only acts perpendicular to the path.
        // Without this the bias would drive the image ALONG the path and destroy the string.
        if (m_mtd_transverse)
            img->setMTDProjection(tau[i]);
    }

    // 2) Advance every live image by one MD step (independent thermostats/velocities).
    // Clamped endpoints are null (see the note above) and are skipped.
    bool any_running = false;
    // Claude Generated 2026: in the decoupled scheme the windows share nothing, so the MD
    // steps may run concurrently. Only step() is parallelised - it contains the energy and
    // gradient evaluation and therefore essentially all of the cost. Everything that
    // touches shared state (bias, accumulators, RNG, trace) stays on the serial path.
    std::vector<int> live;
    live.reserve(m_nimages);
    for (int i = 0; i < m_nimages; ++i)
        if (!m_image_done[i] && m_images[i])
            live.push_back(i);

    // Claude Generated 2026: scheme 'string' as well - between two updates every image is
    // coupled only to its frozen frame M_i (restraint applied above, serially).
    if ((m_independent || m_string) && m_window_threads > 1 && live.size() > 1) {
        std::vector<char> ok(live.size(), 1);
#pragma omp parallel for num_threads(m_window_threads) schedule(dynamic)
        for (int k = 0; k < static_cast<int>(live.size()); ++k)
            ok[k] = m_images[live[k]]->step() ? 1 : 0;
        for (size_t k = 0; k < live.size(); ++k) {
            if (!ok[k])
                m_image_done[live[k]] = true;
            else
                any_running = true;
        }
        killRunawayWindows(live);
    } else {
        for (int i : live) {
            const bool ok = m_images[i]->step(); // false at max_time / instability / stop-file
            if (!ok)
                m_image_done[i] = true;
            else
                any_running = true;
        }
        killRunawayWindows(live);
    }

    // Claude Generated 2026: string stage - a dead image is frozen (never integrated, its
    // mean-path frame never updated again), so every later update would measure a small
    // Delta_max on a partly unrelaxed path and report a convergence that does not hold.
    // With string_abort_on_dead the stage stops here; start() writes no string path and
    // failed() turns into a non-zero exit code (pipeline: stage 2 is not started).
    if (m_string && m_string_abort_on_dead) {
        std::string dead;
        for (int i : live)
            if (m_window_died[i] == m_step)
                dead += fmt::format("{}image {} (T = {:.4g} K)", dead.empty() ? "" : ", ", i, m_window_died_T[i]);
        if (!dead.empty()) {
            m_string_aborted = true;
            CurcumaLogger::set_verbosity(m_band_verbosity);
            CurcumaLogger::error(fmt::format(
                "NEB-MD string: ABORTED - dead image(s) at MD step {} after {} mean-path update(s): {}. "
                "A dead image is frozen and its part of the path would stay unrelaxed, so the string "
                "stage is stopped (string_abort_on_dead = true) and no string path is written. "
                "Hints: reduce string_alpha and/or string_max_step (now {:.3f} A) so the restraint "
                "centres move in smaller steps; check the endpoints and the start path (path_file) "
                "for clashes or a wrong atom order; lower string_k (now {:.4f}); see "
                "<base>.neb.string_history.csv for the update that preceded the failure.",
                m_step, m_string_updates, dead, m_string_max_step, m_string_k));
            return false;
        }
    }

    // (The constraint is no longer applied here - it lives inside SimpleMD::Rattle().)

    ++m_step;
    if (m_window_threads > 1)
        CurcumaLogger::set_verbosity(m_band_verbosity);

    // Claude Generated 2026: string reparametrisation - redistribute the images to equal
    // arclength so they stay well-defined umbrella windows (required for a meaningful PMF).
    if (m_reparametrize > 0 && m_step % m_reparametrize == 0)
        reparametrizeString();

    // Claude Generated 2026: scheme 'string' - accumulate the aligned image averages every
    // step and update the mean path every string_update steps; stop on convergence.
    if (m_string) {
        for (int i : live)
            if (m_images[i] && !m_image_done[i])
                accumulateStringMean(i);
        if (m_step % m_string_update == 0 && updateStringPath()) {
            if (m_step % m_dump == 0) {
                writePathSnapshot(m_step);
                writeEnergyTable(m_step);
            }
            return false; // converged
        }
    }

    // Umbrella histograms (after equilibration): record s(R) of every window.
    // Endpoints are included when they are integrated (fixed_endpoints=false).
    // Calibration phase: measure the local slope, then set the per-window stiffness once.
    if (m_adaptive_k && m_independent && m_path_cv) {
        if (m_step > m_pmf_equilibration / 2 && m_step < m_pmf_equilibration)
            accumulateCalibration(live);
        else if (m_step == m_pmf_equilibration)
            calibrateStiffness();
    }

    if (m_umbrella && m_path_cv && m_step > m_pmf_equilibration) {
        std::vector<int> acc;
        acc.reserve(m_nimages);
        for (int i = 0; i < m_nimages; ++i)
            if (m_images[i] && !m_image_done[i])
                acc.push_back(i);

        // Claude Generated 2026: the accumulation is the expensive half of the window
        // block - it evaluates the path CV with its gradient plus the finite-difference
        // probes for div b, which on a small molecule costs more than the force field
        // itself. It is therefore parallelised as well, with one buffer per image that is
        // merged afterwards; the histogram entries are already per image and race-free.
        std::vector<FeLocal> loc;
        if (m_independent) {
            loc.resize(acc.size());
            for (auto& L : loc) {
                L.count.assign(m_fe_bins, 0);
                L.force.assign(m_fe_bins, 0.0);
                L.div.assign(m_fe_bins, 0.0);
                L.wall.assign(m_fe_bins, 0.0);
                L.ssum.assign(m_fe_bins, 0.0);
                L.wsum.assign(m_fe_bins, 0.0);
                L.rejected_bin.assign(m_fe_bins, 0);
            }
        }

#pragma omp parallel for num_threads(m_window_threads) schedule(dynamic) \
    if (m_independent && m_window_threads > 1 && acc.size() > 1)
        for (int k = 0; k < static_cast<int>(acc.size()); ++k) {
            const int i = acc[k];
            const double s = m_path_cv->s(m_images[i]->positions());
            m_umbrella_ssum[i] += s;
            ++m_umbrella_n[i];
            const int b = static_cast<int>(s * m_umbrella_bins);
            if (b >= 0 && b < m_umbrella_bins)
                ++m_hist[i][b];
            if (m_independent)
                accumulateMeanForce(i, m_images[i], loc[k], m_probe_rng[i]);
            // Claude Generated 2026: window diagnostics (log/CSV only, no effect on the dynamics)
            if (m_independent && m_step % kDiagInterval == 0)
                accumulateWindowDiagnostics(i, m_images[i]->positions());
        }

        // Serial merge.
        for (size_t k = 0; k < loc.size(); ++k) {
            const int i = acc[k];
            for (int b = 0; b < m_fe_bins; ++b) {
                m_fe_count[b] += loc[k].count[b];
                m_fe_force[b] += loc[k].force[b];
                m_fe_div[b] += loc[k].div[b];
                m_fe_wall[b] += loc[k].wall[b];
                m_fe_ssum[b] += loc[k].ssum[b];
                m_fe_wsum[b] += loc[k].wsum[b];
            }
            m_win_count[i] += loc[k].win_n;
            m_win_s_sum[i] += loc[k].win_s;
            m_win_s2_sum[i] += loc[k].win_s2;
            m_win_z_sum[i] += loc[k].win_z;
            m_win_z_max[i] = std::max(m_win_z_max[i], loc[k].win_z_max);
            m_win_zwall_n[i] += loc[k].win_zwall_n;
            m_win_wall_sum[i] += loc[k].win_wall;
            m_win_w[i] += loc[k].win_w;
            m_win_wf[i] += loc[k].win_wf;
            m_win_wd[i] += loc[k].win_wd;
            m_win_wwall[i] += loc[k].win_wwall;
            m_fe_rejected += loc[k].rejected;
            m_win_rejected[i] += loc[k].rejected;
            for (int b = 0; b < m_fe_bins && b < static_cast<int>(loc[k].rejected_bin.size()); ++b)
                m_fe_rejected_bin[b] += loc[k].rejected_bin[b];
            m_fe_nonfinite += loc[k].nonfinite;
        }
    }

    if (m_step % m_dump == 0) {
        writePathSnapshot(m_step);
        writeEnergyTable(m_step);
    }
    return any_running;
}

Geometry NebMD::projectOutTangent(const Geometry& F, const Geometry& tau) const
{
    // Claude Generated 2026: remove the component of F along the (unit) tangent.
    double dot = 0.0;
    for (int a = 0; a < F.rows(); ++a)
        dot += F.row(a).dot(tau.row(a));
    return F - dot * tau;
}

double NebMD::bandEnergy(EnergyCalculator& ec, const Geometry& geom, Geometry* grad) const
{
    ec.updateGeometry(geom);
    const double e = ec.CalculateEnergy(grad != nullptr);
    if (grad)
        *grad = ec.Gradient();
    return e;
}

bool NebMD::optimizeBand()
{
    // Claude Generated 2026: classical NEB optimisation of the band (CI-NEB, FIRE).
    //
    // Why FIRE and not L-BFGS: the NEB force is
    //     F_i = -g_perp + (k(|R_{i+1}-R_i| - |R_i-R_{i-1}|)) tau
    // which is not the gradient of any scalar function, so a line search has no energy to
    // minimise along the search direction. FIRE (Bitzek 2006) is a damped-MD optimiser that
    // only ever consumes the force, which is why it is the standard NEB optimiser.
    CitationRegistry::cite("neb");
    CitationRegistry::cite("fire");
    if (m_climbing)
        CitationRegistry::cite("cineb");

    const int nat = static_cast<int>(m_start_geom.rows());
    EnergyCalculator ec(m_method, m_controller, Basename());
    Molecule probe = m_start;
    ec.setMolecule(probe.getMolInfo());

    // FIRE state
    std::vector<Geometry> vel(m_nimages, Geometry::Zero(nat, 3));
    double dt = 0.1, alpha = 0.1;
    int n_positive = 0;
    const double dt_max = 0.5, dt_min = 1e-5, f_inc = 1.1, f_dec = 0.5, alpha_start = 0.1, f_alpha = 0.99;
    const int n_min = 5;

    std::vector<double> E(m_nimages, 0.0);
    std::vector<Geometry> G(m_nimages, Geometry::Zero(nat, 3));
    bool converged = false;
    int iter = 0;
    double fmax_best = std::numeric_limits<double>::max();
    int stall_count = 0;
    int m_climb_locked = -1;                       // image index once the climber is latched
    const double m_climb_threshold = 10.0 * m_opt_fmax; // switch CI on once pre-relaxed

    for (; iter < m_opt_iterations; ++iter) {
        // Energies + gradients of all images (endpoints only once - they never move).
        for (int i = 0; i < m_nimages; ++i) {
            const bool endpoint = (i == 0 || i == m_nimages - 1);
            if (endpoint && iter > 0)
                continue;
            E[i] = bandEnergy(ec, m_image_geoms[i], &G[i]);
        }

        // Improved tangents (energy-weighted) from the current band.
        std::vector<Geometry> tau = computeTangents(m_image_geoms, E);

        // Climbing image: the highest interior image.
        int iclimb = -1;
        if (m_climbing) {
            // Claude Generated 2026: latch the climbing image instead of re-picking it
            // every iteration.
            //
            // Re-selecting each iteration is unstable whenever two images are nearly
            // degenerate in energy - which is exactly the case for a symmetric path
            // (helicene: the two central images are equal by symmetry). The identity of
            // the climber then flips back and forth, and with it the SIGN of its parallel
            // force, so the band oscillates between two states and can never converge.
            // Observed signature: fmax alternating between two fixed values (0.02502 /
            // 0.02738) with dt pinned at its floor.
            //
            // Standard practice (Henkelman, Uberuaga & Jonsson, JCP 113, 9901 (2000)):
            // relax the plain NEB first, then switch the climbing image on and keep it.
            if (m_climb_locked >= 0) {
                iclimb = m_climb_locked;
            } else if (fmax_best < m_climb_threshold) {
                iclimb = 1;
                for (int i = 2; i < m_nimages - 1; ++i)
                    if (E[i] > E[iclimb])
                        iclimb = i;
                m_climb_locked = iclimb;
                if (getVerbosity() >= 2)
                    CurcumaLogger::info(fmt::format(
                        "NEB: climbing image locked to {} at fmax = {:.5f}", iclimb, fmax_best));
            } else {
                iclimb = -1; // plain NEB until the band is pre-relaxed
            }
        }

        double fmax = 0.0;
        std::vector<Geometry> F(m_nimages, Geometry::Zero(nat, 3));
        for (int i = 1; i < m_nimages - 1; ++i) {
            const Geometry Ftrue = -G[i];
            double dot = 0.0;
            for (int a = 0; a < nat; ++a)
                dot += Ftrue.row(a).dot(tau[i].row(a));

            if (i == iclimb) {
                // CI-NEB: no spring, and the parallel component is INVERTED so the image
                // climbs uphill along the path onto the saddle point.
                F[i] = Ftrue - 2.0 * dot * tau[i];
            } else {
                const double dp = (m_image_geoms[i + 1] - m_image_geoms[i]).norm();
                const double dm = (m_image_geoms[i] - m_image_geoms[i - 1]).norm();
                // Claude Generated 2026: scale k with the image count. The spring force
                // is k*(dp-dm) and the spacing dp,dm shrink like 1/N, so with a fixed k
                // the spring weakens relative to the true force as N grows - the band
                // becomes effectively softer and the barrier drifts with nimages (the
                // artefact seen in the image-count scans). Scaling k by (N-1) keeps the
                // band stiffness invariant.
                // Reference: Henkelman, G.; Jonsson, H., J. Chem. Phys. 113, 9978 (2000),
                // Sec. II - the spring constant is defined per band so that the
                // distribution of images is independent of their number; see also
                // Henkelman, Uberuaga & Jonsson, J. Chem. Phys. 113, 9901 (2000).
                const double k_eff = m_opt_k * (m_nimages - 1);
                F[i] = (Ftrue - dot * tau[i]) + k_eff * (dp - dm) * tau[i];
            }
            for (int a = 0; a < nat; ++a)
                fmax = std::max(fmax, F[i].row(a).norm());
        }

        if (fmax < m_opt_fmax) {
            converged = true;
            break;
        }

        // Claude Generated 2026: stall detection. FIRE can settle into a state where the
        // force no longer decreases (typically the climbing image fighting the springs).
        // Continuing to opt_iterations then just wastes time, so stop and say so - an
        // honest "stalled at fmax=X" is more useful than a silent full-length run.
        if (fmax < fmax_best - 1e-6) {
            fmax_best = fmax;
            stall_count = 0;
        } else if (++stall_count >= 200) {
            CurcumaLogger::warn(fmt::format(
                "NEB optimisation stalled: fmax stuck at {:.5f} Eh/A (threshold {:.5f}) for 200 "
                "iterations - stopping. The band is relaxed but the criterion is not reachable "
                "with the current settings (try a larger opt_k, more images, or a looser opt_fmax).",
                fmax, m_opt_fmax));
            break;
        }

        // --- FIRE update (Bitzek et al., PRL 97, 170201 (2006), semi-implicit Euler) ---
        //
        // FIRE is damped MD:  v <- v + dt*F  (MD part), then the velocity is mixed
        // towards the force direction, and dt/alpha adapt to the sign of P = v.F.
        // The band is ONE system: P, |v| and |F| are global over all moving images,
        // and the trust radius must scale the whole band uniformly - scaling images
        // individually would distort the path.
        double P = 0.0;
        for (int i = 1; i < m_nimages - 1; ++i)
            for (int a = 0; a < nat; ++a)
                P += vel[i].row(a).dot(F[i].row(a));

        if (P > 0.0) {
            if (++n_positive > n_min) {
                dt = std::min(dt * f_inc, dt_max);
                alpha *= f_alpha;
            }
        } else {
            // Uphill: discard the velocity, shrink dt, reset the mixing.
            n_positive = 0;
            dt = std::max(dt * f_dec, dt_min);
            alpha = alpha_start;
            for (int i = 1; i < m_nimages - 1; ++i)
                vel[i].setZero();
        }

        // 1) MD half: v <- v + dt*F
        for (int i = 1; i < m_nimages - 1; ++i)
            vel[i] += dt * F[i];

        // 2) Mixing: v <- (1-alpha) v + alpha |v| Fhat, with |v| and |F| taken AFTER
        //    the MD step (they must refer to the same velocity that is being mixed -
        //    using the pre-update norms is inconsistent and was one reason the band
        //    kept overshooting).
        double vnorm = 0.0, fnorm = 0.0;
        for (int i = 1; i < m_nimages - 1; ++i) {
            vnorm += vel[i].squaredNorm();
            fnorm += F[i].squaredNorm();
        }
        vnorm = std::sqrt(vnorm);
        fnorm = std::sqrt(fnorm);
        if (fnorm > 1e-12) {
            for (int i = 1; i < m_nimages - 1; ++i)
                vel[i] = (1.0 - alpha) * vel[i] + (alpha * vnorm / fnorm) * F[i];
        }

        // 3) Position update with a GLOBAL trust radius: find the largest displacement
        //    anywhere in the band and scale every image by the same factor, so the
        //    relative geometry of the path is preserved.
        double dmax = 0.0;
        for (int i = 1; i < m_nimages - 1; ++i)
            dmax = std::max(dmax, (dt * vel[i]).rowwise().norm().maxCoeff());
        const double trust = 0.02; // Angstrom per iteration (conservative: large steps in a
                                   // stiff saddle region were driving the FIRE oscillation)
        const double scale = (dmax > trust) ? trust / dmax : 1.0;
        for (int i = 1; i < m_nimages - 1; ++i)
            m_image_geoms[i] += (dt * scale) * vel[i];

        if (getVerbosity() >= 2 && iter % 25 == 0)
            CurcumaLogger::info(fmt::format("NEB opt iter {:4d}: fmax = {:.5f} Eh/A, dt = {:.3f}",
                                            iter, fmax, dt));
    }

    // Report the optimised barrier - this is the classical dE(TS), independent of any MD.
    for (int i = 0; i < m_nimages; ++i)
        E[i] = bandEnergy(ec, m_image_geoms[i], nullptr);
    int imax = 0;
    for (int i = 1; i < m_nimages; ++i)
        if (E[i] > E[imax])
            imax = i;
    const double kcal = 627.5094740631;

    CurcumaLogger::result_fmt("NEB optimisation: {} after {} iterations",
                              converged ? "converged" : "NOT converged (raised opt_iterations?)", iter);
    // The barrier is E(highest image) - E(first image), i.e. a HEI estimate: no
    // interpolation between images, no saddle refinement. With the climbing image on,
    // that image is driven ONTO the saddle, so the estimate is meaningful only if the
    // band actually converged. Two failure modes to keep in mind:
    //   - band not converged  -> the HEI hangs above the MEP -> barrier TOO HIGH
    //   - climbing image off  -> the springs pull the top image off the saddle -> TOO LOW
    // Neither is detectable from the energy alone; a frequency calculation on the HEI is
    // the only real check (exactly one imaginary mode = saddle). Measured on helicene the
    // HEI had ZERO imaginary modes, i.e. the reported barrier was not a barrier at all.
    CurcumaLogger::result_fmt("NEB barrier dE(forward) = {:.2f} kcal/mol at image {} (climbing image: {})",
                              (E[imax] - E[0]) * kcal, imax, m_climbing ? "on" : "off");
    CurcumaLogger::result("  (highest-image estimate - verify with a frequency calculation: "
                          "exactly one imaginary mode means it really is a saddle point)");
    if (!converged)
        CurcumaLogger::warn("Band did not converge - the highest-image energy is an UPPER BOUND, "
                            "not the barrier.");
    CurcumaLogger::result_fmt("NEB barrier dE(reverse) = {:.2f} kcal/mol; reaction dE = {:.2f} kcal/mol",
                              (E[imax] - E[m_nimages - 1]) * kcal,
                              (E[m_nimages - 1] - E[0]) * kcal);

    // Persist the optimised path with its true energies (the MD images do not exist yet,
    // so writeBandFile() cannot supply them here).
    {
        const std::string path = outputPath(Basename() + ".neb.opt.xyz");
        for (int i = 0; i < m_nimages; ++i) {
            Molecule mol = m_start;
            mol.setGeometry(m_image_geoms[i]);
            mol.appendXYZFile(path, fmt::format("neb_opt image={} E={:.6f} dE={:.2f}kcal",
                                                i, E[i], (E[i] - E[0]) * kcal));
        }
        addBakFile(Basename() + ".neb.opt.xyz");
    }
    return converged;
}

void NebMD::setupUmbrella()
{
    // Claude Generated 2026: freeze the current (NEB-optimised) path as the reference
    // for the path collective variable and place one umbrella window per image.
    CitationRegistry::cite("pathcv");
    CitationRegistry::cite("wham");
    // Claude Generated 2026: the CV is built on cvFrames() - the image frames by default, the
    // re-sampled / extended frames when cv_frames / cv_extend are set (buildCvPath).
    const std::vector<Geometry>& cvp = cvFrames();
    m_path_cv = makePathCV(cvp); // Claude Generated 2026: rmsd or distance metric
    m_umbrella_s0.assign(m_nimages, 0.0);
    // Claude Generated 2026: the window centres must be the s values the reference frames
    // ACTUALLY have, not i/(N-1).
    //
    // s = sum_i (i-1) w_i / sum_i w_i / (N-1) equals i/(N-1) at frame i only in the limit
    // lambda -> infinity, where each frame's Gaussian excludes all others. At finite
    // lambda the neighbouring frames still carry weight and s is compressed towards the
    // middle: measured for a 3-frame ethane path with lambda*d^2 = 1, the first frame has
    // s = 0.148 and the last s = 0.855. Windows placed at 0 and 1 then target values that
    // no structure on the path can reach, so the restraint pushes permanently against an
    // unreachable setpoint and the image is driven off the path. Evaluating s at the
    // frames removes the problem exactly and costs one CV evaluation per image.
    for (int i = 0; i < m_nimages; ++i)
        m_umbrella_s0[i] = m_path_cv->s(m_image_geoms[i]);
    {
        const double span = m_umbrella_s0.back() - m_umbrella_s0.front();
        if (span < 0.5)
            CurcumaLogger::warn(fmt::format(
                "Path CV: the reference frames only span s = {:.3f} .. {:.3f}. With so few frames "
                "(or so small a lambda) the coordinate resolves the path poorly; consider more "
                "frames for the CV than windows for the sampling.",
                m_umbrella_s0.front(), m_umbrella_s0.back()));
    }
    // Claude Generated 2026: auto-select the force constant so that neighbouring windows
    // actually OVERLAP - without overlap WHAM has nothing to stitch and returns garbage.
    // The thermal width of a harmonic window is sigma = sqrt(kT/k); requiring
    // sigma ~ ds/2 with ds = 1/(N-1) the window spacing gives k = kT/(ds/2)^2.
    // (Measured failure of the previous fixed default k=200: sigma/ds = 0.02, nine gaps
    // in the histogram, WHAM barrier 351 kcal/mol - meaningless.)
    if (m_umbrella_k <= 0.0) {
        const double kT = thermalEnergy();
        const double ds = 1.0 / std::max(m_nimages - 1, 1);
        m_umbrella_k = kT / (0.5 * ds * 0.5 * ds);
        CurcumaLogger::info(fmt::format(
            "NEB-MD: umbrella k auto-set to {:.4f} (sigma_s ~ {:.3f}, window spacing {:.3f})",
            m_umbrella_k, std::sqrt(kT / m_umbrella_k), ds));
    }
    m_umbrella_kvec.assign(m_nimages, m_umbrella_k);
    // Hand the constraint to the integrator, where it belongs (Claude Generated 2026):
    // SimpleMD applies it inside Rattle(), i.e. before the kinetic energy is formed and
    // before the thermostat runs. Applying it from outside after step() left the injected
    // energy beyond the thermostat's reach and heated every constrained window to death.
    if (m_constrain) {
        for (int i = 0; i < m_nimages; ++i) {
            if (!m_images[i])
                continue;
            PathCV* cv = m_path_cv.get();
            m_images[i]->setScalarConstraint(
                [cv](const Geometry& R, Geometry& grad) {
                    double sv = 0.0, zv = 0.0;
                    cv->evaluate(R, &sv, &zv, &grad, nullptr);
                    return sv;
                },
                m_umbrella_s0[i], m_constrain_tol);
        }
        CurcumaLogger::result_fmt("NEB-MD: s held by a holonomic constraint (SHAKE/RATTLE inside "
                                  "the integrator), tolerance {:.1e}", m_constrain_tol);
    }
    // Claude Generated 2026: flat-bottom z wall (tube), W = k/2 (z - z_wall)^2 beyond z_wall.
    // Handed to the integrator as part of the potential (SimpleMD::setExtraPotential), so it is
    // evaluated at the new positions inside Verlet()/Rattle(), before the RATTLE velocity
    // projection, and its gradient is contained in trueGradient(): the mean-force estimator
    // then averages b.grad(E + W), the force of the ensemble exp(-(E+W)/kT) actually sampled.
    if (m_z_wall > 0.0) {
        if (m_z_wall_k <= 0.0)
            m_z_wall_k = PathCV::zWallDefaultK(thermalEnergy(), m_z_wall);
        for (int i = 0; i < m_nimages; ++i) {
            if (!m_images[i])
                continue;
            const PathCV* cv = m_path_cv.get();
            const double zw = m_z_wall, kw = m_z_wall_k;
            m_images[i]->setExtraPotential([cv, zw, kw](const Geometry& R, Geometry& grad) {
                return cv->zWallEnergy(R, zw, kw, &grad);
            });
        }
        const bool inv = m_path_cv->zUnit() == "A^-2"; // Claude Generated 2026: unit of z
        CurcumaLogger::result_fmt("NEB-MD: z wall at z = {:.4f} {}, k = {:.4f} {} (W = {:.2f} kT at "
                                  "z = 1.2 z_wall); the PMF is that of the z-restricted system",
                                  m_z_wall, m_path_cv->zUnit(), m_z_wall_k, inv ? "Eh A^4" : "Eh/A^4",
                                  0.5 * m_z_wall_k * 0.04 * m_z_wall * m_z_wall / thermalEnergy());
    }
    m_hist.assign(m_nimages, std::vector<int>(m_umbrella_bins, 0));
    m_umbrella_ssum.assign(m_nimages, 0.0);
    m_umbrella_n.assign(m_nimages, 0);
    // Claude Generated 2026: window diagnostics (windows.csv). Reference step of J = the largest
    // change of one diagnostic pair coordinate between consecutive CV frames.
    m_diag_z.assign(m_nimages, std::vector<double>());
    m_diag_nn_sum.assign(m_nimages, 0.0);
    m_diag_nn2_sum.assign(m_nimages, 0.0);
    m_diag_n.assign(m_nimages, 0);
    m_diag_x_sum.assign(m_nimages, Vector::Zero(static_cast<int>(m_diag_pairs.size())));
    m_diag_x_step = 0.0;
    if (m_diag_pairs.size() >= 2)
        for (size_t k = 1; k < cvp.size(); ++k) {
            const Vector dx = PathCV::pairValues(cvp[k], m_diag_pairs, m_diag_transform)
                - PathCV::pairValues(cvp[k - 1], m_diag_pairs, m_diag_transform);
            m_diag_x_step = std::max(m_diag_x_step, dx.cwiseAbs().maxCoeff());
        }
    CurcumaLogger::result_fmt("NEB-MD: umbrella sampling on the path CV (lambda = {:.3f}, k = {:.1f}, "
                              "{} windows)", m_path_cv->lambda(), m_umbrella_k, m_nimages);

    // Claude Generated 2026: conditioning of the path coordinate.
    //
    // The path CV can only tell two images apart if they are further apart than the
    // molecule wanders on its own. The reference spacing is known here; the thermal
    // amplitude is measured during the run and the ratio is reported at the end.
    {
        // Spacing and Z floor refer to the CV frames (identical to the image frames by default).
        // With cv_atoms the spacing runs over the selected atoms (the metric of the CV).
        const int nat = static_cast<int>(cvp[0].rows());
        const int ncv = static_cast<int>(cvp.size());
        double dsum = 0.0;
        for (int i = 1; i < ncv; ++i) {
            double d2 = 0.0;
            if (m_cv_metric == PathCV::Metric::Distances) { // Claude Generated 2026
                dsum += cvSpacing(cvp[i], cvp[i - 1]);
            } else if (m_cv_atoms.empty()) {
                for (int a = 0; a < nat; ++a)
                    d2 += (cvp[i].row(a) - cvp[i - 1].row(a)).squaredNorm();
                dsum += std::sqrt(d2 / nat);
            } else {
                for (int a : m_cv_atoms)
                    d2 += (cvp[i].row(a) - cvp[i - 1].row(a)).squaredNorm();
                dsum += std::sqrt(d2 / static_cast<double>(m_cv_atoms.size()));
            }
        }
        m_path_spacing = dsum / std::max(ncv - 1, 1);
        const char* unit = (m_cv_metric == PathCV::Metric::Distances && m_cv_transform == PathCV::PairTransform::Inverse)
            ? "1/A" : "A";
        CurcumaLogger::result_fmt("NEB-MD: reference path spacing {:.4f} {} per image "
                                  "(total {:.3f} {}, lambda*d^2 = {:.2f})",
                                  m_path_spacing, unit, dsum, unit, m_path_cv->lambda() * m_path_spacing * m_path_spacing);

        // Floor for |grad s|^2/m below which a sample is unusable. Evaluated at the
        // MIDPOINTS between reference frames: s has a plateau at each frame, so the
        // gradient is largest in between - that is the scale to normalise against.
        std::vector<double> Zmid;
        for (int i = 1; i < ncv; ++i) {
            const Geometry mid = 0.5 * (cvp[i] + cvp[i - 1]);
            double Z = 0.0;
            pathBiasField(mid, &Z);
            if (Z > 0.0)
                Zmid.push_back(Z);
        }
        if (!Zmid.empty()) {
            std::sort(Zmid.begin(), Zmid.end());
            m_fe_Zfloor = 1e-3 * Zmid[Zmid.size() / 2];
        }
    }
}


void NebMD::accumulateCalibration(const std::vector<int>& live)
{
    // Claude Generated 2026: read off the local slope of the profile during equilibration.
    // Legitimate because a restraint depending only on s leaves conditional averages
    // untouched, so <b.grad E> is the unbiased mean force whatever the current stiffness.
    if (!m_path_cv)
        return;
    for (int i : live) {
        if (!m_images[i] || m_image_done[i])
            continue;
        const Geometry R = m_images[i]->positions();
        const int nat = static_cast<int>(R.rows());
        double sv = 0.0, zv = 0.0;
        Geometry ds;
        m_path_cv->evaluate(R, &sv, &zv, &ds, nullptr);
        double Z = 0.0;
        for (int a = 0; a < nat; ++a)
            Z += ds.row(a).squaredNorm() / std::max(m_masses[a], 1e-12);
        if (Z < m_fe_Zfloor)
            continue;
        const Geometry g = m_images[i]->trueGradient();
        double bdotg = 0.0;
        for (int a = 0; a < nat; ++a)
            bdotg += ds.row(a).dot(g.row(a)) / (std::max(m_masses[a], 1e-12) * Z);
        if (!std::isfinite(bdotg))
            continue;
        m_calib_dg_sum[i] += bdotg;
        ++m_calib_n[i];
    }
}

void NebMD::calibrateStiffness()
{
    // k_i > |dG/ds|_i / (ds/2): stiff where the profile is steep, soft where it is flat.
    // A single global value cannot satisfy both ends at once, which is why the basin
    // windows held and the flank windows slid in every run of this work.
    double s0min = 1e30, s0max = -1e30;
    for (int i = 0; i < m_nimages; ++i) {
        s0min = std::min(s0min, m_umbrella_s0[i]);
        s0max = std::max(s0max, m_umbrella_s0[i]);
    }
    const double dsw = (s0max - s0min) / std::max(m_nimages - 1, 1);
    const double kmax = m_adaptive_k_max * m_umbrella_k;
    bool capped = false;
    CurcumaLogger::result("NEB-MD: calibrated per-window stiffness (from the measured local slope)");
    for (int i = 0; i < m_nimages; ++i) {
        if (m_calib_n[i] < 10)
            continue;
        const double dg = std::abs(m_calib_dg_sum[i] / static_cast<double>(m_calib_n[i]));
        double k_i = dg / std::max(0.5 * dsw, 1e-12);
        k_i = std::max(k_i, m_umbrella_k); // never softer than requested
        if (k_i > kmax) { k_i = kmax; capped = true; }
        m_umbrella_kvec[i] = k_i;
        CurcumaLogger::result(fmt::format("    window {:>2}: |dG/ds| = {:.4f} Eh -> k = {:.3f}",
                                          i, dg, k_i));
    }
    // Put every window back on its reference frame before production starts.
    //
    // Claude Generated 2026: without this the calibration arrives too late. During the
    // equilibration phase a window on a steep flank slides away, and once it has crossed a
    // maximum it cannot be recovered by a restraint sized from the LOCAL slope - measured
    // on helicene, the restoring force was already three times the local dG/ds while the
    // window kept sitting 1.45 spacings away, because returning would have meant climbing
    // back over the barrier. Resetting the geometry costs nothing statistically: the
    // equilibration samples are discarded anyway, and the reference frame is by definition
    // the place the window is meant to sample. Velocities are kept and only rescaled if the
    // discontinuous coordinate change has thrown the temperature off.
    int reset = 0;
    for (int i = 0; i < m_nimages; ++i) {
        if (!m_images[i] || m_image_done[i])
            continue;
        const double T_before = m_images[i]->currentTemperature();
        m_images[i]->setPositions(m_image_geoms[i]);
        const double T_target = m_images[i]->targetTemperature();
        if (T_before > 1e-6 && T_target > 1e-6 && T_before > 2.0 * T_target)
            m_images[i]->rescaleVelocities(std::sqrt(T_target / T_before));
        m_prev_s[i] = -1.0; // the jump detector must not see the reset as a CV jump
        ++reset;
    }
    if (reset > 0)
        CurcumaLogger::result_fmt("NEB-MD: {} windows reset onto their reference frames before "
                                  "production (calibration finished)", reset);

    if (capped)
        CurcumaLogger::warn("The calibrated stiffness hit its upper bound. The profile is steeper "
                            "than the restraint may be made without endangering the window; the "
                            "affected windows will not hold and the barrier is a lower bound.");
}


bool NebMD::constrainWindow(int i)
{
    // Claude Generated 2026: hold s(R) = s0 exactly (SHAKE), then remove the velocity
    // component along the constraint (RATTLE).
    //
    //   SHAKE:   R_a <- R_a - (g/Z) grad_a s / m_a   iterated until |g| = |s - s0| < tol
    //   RATTLE:  v_a <- v_a - (nu)  grad_a s / m_a ,  nu = (sum_a v_a.grad_a s) / Z
    //
    // with Z = sum_a |grad_a s|^2 / m_a. Both steps move along the same direction as the
    // constraint force, which is what makes them the discrete form of a Lagrange
    // multiplier rather than an arbitrary correction.
    //
    // This is applied after the integrator has taken its step, not inside its half-kicks,
    // so time reversibility is only approximate - the same caveat that the injected spring
    // force carries (see the note at the top of nebmd.h). It is nevertheless a genuine
    // constraint: s is exact at every point where anything is measured.
    SimpleMD* img = m_images[i];
    if (!img || !m_path_cv)
        return false;
    Geometry R = img->positions();
    const Geometry R_before = R;
    const int nat = static_cast<int>(R.rows());
    const double s_target = m_umbrella_s0[i];

    double sv = 0.0, zv = 0.0;
    Geometry ds;
    int it = 0;
    for (; it < 100; ++it) {
        m_path_cv->evaluate(R, &sv, &zv, &ds, nullptr);
        const double g = sv - s_target;
        if (std::abs(g) < m_constrain_tol)
            break;
        double Z = 0.0;
        for (int a = 0; a < nat; ++a)
            Z += ds.row(a).squaredNorm() / std::max(m_masses[a], 1e-12);
        if (Z < 1e-30)
            return false;
        const double mu = g / Z;
        for (int a = 0; a < nat; ++a)
            R.row(a) -= (mu / std::max(m_masses[a], 1e-12)) * ds.row(a);
    }
    if (it >= 100) {
        ++m_shake_fail;
        return false;
    }

    // The position correction must be carried into the velocities: dv = dr/dt, exactly as
    // SimpleMD::Rattle() does for its bond constraints (the m_dT_inverse factors there).
    //
    // Claude Generated 2026 - CORRECTED. The first version moved the positions and then
    // projected the velocities independently. That is not a constraint but two unrelated
    // corrections, and it pumps energy into the system at every step: measured, three of
    // four windows heated past 3000 K and were destroyed. Deriving the velocity change
    // from the position change is what makes the pair a discrete Lagrange multiplier.
    const Geometry dR = R - R_before;
    Geometry v = img->velocities();
    const double dt = std::max(img->timeStep(), 1e-12) * CurcumaUnit::Constants::FS_TO_MD_TIME; // velocities are in integrator units
    v += dR / dt;
    img->setPositions(R);

    // Second RATTLE stage: remove the remaining velocity component along the constraint,
    // so that ds/dt = 0 and the coordinate stays put over the following step.
    m_path_cv->evaluate(R, &sv, &zv, &ds, nullptr);
    double Z = 0.0, vdot = 0.0;
    for (int a = 0; a < nat; ++a) {
        const double m = std::max(m_masses[a], 1e-12);
        Z += ds.row(a).squaredNorm() / m;
        vdot += v.row(a).dot(ds.row(a));
    }
    if (Z > 1e-30) {
        const double nu = vdot / Z;
        for (int a = 0; a < nat; ++a)
            v.row(a) -= (nu / std::max(m_masses[a], 1e-12)) * ds.row(a);
    }
    img->setVelocities(v);
    return true;
}

void NebMD::killRunawayWindows(const std::vector<int>& live)
{
    // Claude Generated 2026: a window that has blown up must leave the averages.
    //
    // The characteristic failure in this work was not a crash but a transient: the
    // temperature jumped by many orders of magnitude for a short while, the thermostat
    // pulled it back, and the run continued with a damaged molecule. Nothing in the mean
    // values showed it - only the final potential energy did. Such a window is not noisy,
    // it is gone, and everything it contributes afterwards is meaningless.
    double limit = m_max_temperature;
    if (limit <= 0.0) {
        double T0 = 300.0;
        for (int i : live)
            if (m_images[i]) { T0 = std::max(m_images[i]->targetTemperature(), 1.0); break; }
        limit = 10.0 * T0;
    }
    for (int i : live) {
        if (!m_images[i] || m_image_done[i])
            continue;
        const double T = m_images[i]->currentTemperature();
        if (!std::isfinite(T) || T > limit) {
            m_image_done[i] = true;
            m_window_died[i] = m_step;
            if (i < static_cast<int>(m_window_died_T.size()))
                m_window_died_T[i] = T;
            CurcumaLogger::warn(fmt::format(
                "NEB-MD: window {} declared dead at step {} (T = {:.3g} K > {:.0f} K). It is "
                "removed from the integration and from all averages.", i, m_step, T, limit));
        }
    }
}

double NebMD::biasRampFactor() const
{
    // Claude Generated 2026: adiabatic switch-on of the restraints.
    //
    // Every image starts exactly ON the reference path and exactly AT its setpoint
    // (s0_i = s(gamma_i)), so at t = 0 no restraint force is needed at all. Switching the
    // full stiffness on immediately delivers a shock that can push the image off the path
    // - and once z has grown, all distances d_i grow together, the weights flatten, s
    // becomes ill-defined and grad s small and noisy. The restraint then pushes in a
    // meaningless direction with a magnitude that keeps growing with the deviation, which
    // is how helicene images were torn apart at k = 5.
    //
    // The ramp uses a raised cosine, so it starts and ends with zero slope: no kick at
    // switch-on and no discontinuity when production sampling begins.
    if (m_ramp_steps <= 0 || m_step >= m_ramp_steps)
        return 1.0;
    const double t = static_cast<double>(m_step) / static_cast<double>(m_ramp_steps);
    return 0.5 * (1.0 - std::cos(3.14159265358979323846 * t));
}

void NebMD::applyUmbrellaBias(int i, SimpleMD* img, Geometry& Fext) const
{
    // Progress restraint:  V_s = k_s/2 (s - s0_i)^2  ->  dV/dR = k_s (s - s0) ds/dR
    // applyExternalForces() adds to the GRADIENT, so dV/dR is added directly.
    if (!m_path_cv || !img)
        return;
    double s = 0.0, z = 0.0;
    Geometry ds, dz;
    m_path_cv->evaluate(img->positions(), &s, &z, &ds, m_kz > 0.0 ? &dz : nullptr);
    // Jump diagnostic (Claude Generated 2026): a best-fit RMSD can switch alignment
    // branch, which leaves s continuous but its gradient discontinuous. The restraint
    // then delivers an impulsive force - the suspected origin of the isolated
    // temperature excursions (one dump interval at 1e7 K, then recovery).
    const double gnorm = ds.norm();
    if (m_prev_s[i] >= 0.0) {
        const double jump = std::abs(s - m_prev_s[i]);
        if (jump > m_img_max_jump[i]) {
            m_img_max_jump[i] = jump;
            m_img_max_jump_step[i] = m_step;
        }
    }
    m_prev_s[i] = s;
    if (m_cv_trace) {
        if (!m_cv_trace_file.is_open()) {
            m_cv_trace_file.open(outputPath(Basename() + ".neb.cvtrace.csv"));
            m_cv_trace_file << "step,image,s,ds_step,grad_s_norm,z,T"
                            << (m_z_wall > 0.0 ? ",W_wall\n" : "\n"); // W_wall only with a z wall
        }
        m_cv_trace_file << m_step << "," << i << "," << fmt::format("{:.8f}", s) << ","
                        << fmt::format("{:.8f}", m_prev_s[i] < 0 ? 0.0 : std::abs(s - m_prev_s[i]))
                        << "," << fmt::format("{:.6f}", gnorm) << "," << fmt::format("{:.6f}", z)
                        << "," << fmt::format("{:.2f}", img->currentTemperature());
        if (m_z_wall > 0.0) { // Claude Generated 2026: wall energy in Eh from the same z
            const double ex = std::max(z - m_z_wall, 0.0);
            m_cv_trace_file << "," << fmt::format("{:.6e}", 0.5 * m_z_wall_k * ex * ex);
        }
        m_cv_trace_file << "\n";
    }

    if (m_constrain)
        return; // s is held exactly; an additional restraint would be double bookkeeping
    const double ramp = biasRampFactor();
    const double dev = s - m_umbrella_s0[i];
    const double k_i = windowK(i);
    Geometry bias = (ramp * k_i * dev) * ds;
    // Claude Generated 2026: cap the restraint force. Where the path CV saturates, |grad s|
    // can spike and the restraint then destroys the molecule instead of holding it. The
    // cap keeps the run alive; the counter makes sure the event is not silently ignored,
    // because frequent capping means the coordinate is badly conditioned.
    // The cap must scale the WHOLE vector by one factor, not each atom separately.
    // Per-atom clipping changes the DIRECTION of the force, so it is no longer parallel
    // to grad s - it then stops restraining s and instead deforms the molecule along a
    // direction with no physical meaning. Measured: with per-atom clipping and a stiff
    // restraint, helicene dissociated (Epot -3.2 .. -9.2 Eh instead of -9.2), while the
    // same run with a soft restraint stayed intact.
    const double kMaxBias = m_bias_fmax; // Eh/Angstrom, largest per-atom component
    double fmax = 0.0;
    for (int a = 0; a < bias.rows(); ++a)
        fmax = std::max(fmax, bias.row(a).norm());
    if (fmax > kMaxBias) {
        bias *= kMaxBias / fmax; // uniform scaling: direction preserved
        ++m_img_capped[i];
    }
    Fext += bias;

    // Tube restraint on the PERPENDICULAR coordinate: V_z = k_z/2 * z^2, applied only
    // beyond a tolerance so the image can fluctuate freely inside the tube.
    //
    // This is the piece that keeps the sampling ON the reference path. Without it the
    // images were measured to drift up to 1.11 A off a converged MEP, which is what
    // pushed the mean-force barrier up (sampling above the path) and emptied the
    // barrier bin for WHAM (sampling along it). z is exactly "distance from the path",
    // so restraining it is the direct remedy.
    if (m_kz > 0.0) {
        const double excess = z - m_z_tolerance;
        if (excess > 0.0)
            Fext += (ramp * m_kz * excess) * dz;
    }
}

void NebMD::writeWHAM() const
{
    // Claude Generated 2026: WHAM (Kumar et al. 1992) - self-consistently unbias the
    // window histograms into a single free-energy profile.
    //
    //   P(s) = sum_i H_i(s) / sum_i N_i exp(-(V_i(s) - f_i)/kT)
    //   exp(-f_i/kT) = sum_s P(s) exp(-V_i(s)/kT)
    //
    // iterated until the offsets f_i stop changing; G(s) = -kT ln P(s).
    if (!m_path_cv)
        return;
    const double kT = thermalEnergy(); // Eh at the thermostat setpoint
    const double kcal = 627.5094740631;
    const int B = m_umbrella_bins;

    std::vector<double> Ntot(m_nimages, 0.0), kvec(m_nimages, 0.0);
    for (int i = 0; i < m_nimages; ++i) {
        Ntot[i] = static_cast<double>(m_umbrella_n[i]);
        // Per-window stiffness: calibrateStiffness() sets m_umbrella_kvec[i], and that is
        // the k the dynamics used (applyUmbrellaBias). The global m_umbrella_k here
        // unbiased the histograms with the wrong potential whenever calibration was active.
        kvec[i] = windowK(i);
    }
    const std::vector<double> P = curcuma::whamSolve(m_hist, Ntot, m_umbrella_s0, kvec, kT);

    // G(s) = -kT ln P(s), shifted to zero at its minimum
    std::vector<double> G(B, std::numeric_limits<double>::quiet_NaN());
    double gmin = std::numeric_limits<double>::max();
    for (int b = 0; b < B; ++b) {
        if (P[b] > 0.0) {
            G[b] = -kT * std::log(P[b]) * kcal;
            gmin = std::min(gmin, G[b]);
        }
    }
    int nfilled = 0;
    double gmax = -std::numeric_limits<double>::max();
    int bmax = -1;
    for (int b = 0; b < B; ++b) {
        if (std::isfinite(G[b])) {
            G[b] -= gmin;
            ++nfilled;
            if (G[b] > gmax) { gmax = G[b]; bmax = b; }
        }
    }

    std::ofstream out(outputPath(Basename() + ".neb.wham.csv"));
    if (out) {
        out << "# WHAM free-energy profile on the path collective variable (Branduardi s)\n";
        out << "# umbrella_k=" << m_umbrella_k << " lambda=" << m_path_cv->lambda()
            << " windows=" << m_nimages << "\n";
        if (m_cv_metric == PathCV::Metric::Distances) // Claude Generated 2026 (default header unchanged)
            out << "# cv_metric distances, " << cvPairsText() << "\n";
        out << "s,G_kcal,counts\n";
        for (int b = 0; b < B; ++b) {
            long c = 0;
            for (int i = 0; i < m_nimages; ++i)
                c += m_hist[i][b];
            out << fmt::format("{:.4f}", (b + 0.5) / B) << ","
                << (std::isfinite(G[b]) ? fmt::format("{:.4f}", G[b]) : "nan") << ","
                << c << "\n";
        }
    }

    CurcumaLogger::result("");
    CurcumaLogger::result("=== WHAM free energy on the path CV ===");
    CurcumaLogger::result_fmt("  populated bins   : {} / {}", nfilled, B);
    for (int i = 0; i < m_nimages; ++i) {
        if (m_umbrella_n[i] > 0)
            CurcumaLogger::result_fmt("  window {:2d}: s0 = {:.3f}  <s> = {:.3f}  samples = {}",
                                      i, m_umbrella_s0[i],
                                      m_umbrella_ssum[i] / m_umbrella_n[i], m_umbrella_n[i]);
    }
    if (bmax >= 0)
        CurcumaLogger::result_fmt("  barrier dG       : {:.2f} kcal/mol at s = {:.3f}",
                                  gmax, (bmax + 0.5) / static_cast<double>(B));
    // Claude Generated 2026: quantify window overlap. WHAM stitches the windows together
    // through their shared region; with gaps in the histogram there is nothing to stitch
    // and the resulting profile is meaningless (not merely noisy).
    {
        int gaps = 0, last = -1;
        for (int b = 0; b < B; ++b) {
            long c = 0;
            for (int i = 0; i < m_nimages; ++i)
                c += m_hist[i][b];
            if (c > 0) {
                if (last >= 0 && b - last > 2)
                    ++gaps;
                last = b;
            }
        }
        // Also check whether the windows actually HELD their target: a restraint soft
        // enough to overlap may be too soft to keep an image on a steep barrier flank.
        double max_drift = 0.0;
        for (int i = 1; i < m_nimages - 1; ++i)
            if (m_umbrella_n[i] > 0)
                max_drift = std::max(max_drift,
                                     std::abs(m_umbrella_ssum[i] / m_umbrella_n[i] - m_umbrella_s0[i]));
        const double ds = 1.0 / std::max(m_nimages - 1, 1);

        if (gaps > 0) {
            CurcumaLogger::warn(fmt::format(
                "WHAM: {} gap(s) between populated bins - the umbrella windows do NOT overlap, "
                "so this profile is INVALID.", gaps));
        } else {
            CurcumaLogger::result("  window overlap   : OK (no gaps between populated bins)");
        }
        if (max_drift > ds) {
            CurcumaLogger::warn(fmt::format(
                "WHAM: windows drifted up to {:.3f} in s (spacing {:.3f}) - the restraint is too "
                "soft to hold them on the barrier flank, so parts of the path were never sampled.",
                max_drift, ds));
        }
        if (gaps > 0 || max_drift > ds) {
            // The two requirements pull in opposite directions and cannot both be met
            // with few windows; more windows relax the conflict (measured for a
            // ~60 kcal/mol barrier: k_hold/k_overlap = 11 at N=10, 5.3 at N=20,
            // 2.6 at N=40, 1.3 at N=80).
            CurcumaLogger::warn(
                "WHAM: overlap needs a SOFT restraint while holding a window on a steep flank "
                "needs a STIFF one. With few windows both cannot hold at once - use "
                "substantially more images (umbrella windows are not NEB images; 40-80 is "
                "typical) rather than tuning umbrella_k.");
        }
    }
}

double NebMD::reparametrizeString()
{
    // Claude Generated 2026: string reparametrisation (E, Ren & Vanden-Eijnden 2002;
    // Vanden-Eijnden & Venturoli 2009). Without it the images drift along the path during
    // the MD - they pile up in the minima and thin out over the barrier - so the arclength
    // coordinate of each "umbrella window" keeps changing and the mean-force integral is
    // integrated over a moving grid. Redistributing the images to equal arclength every
    // step is what turns the band into a set of well-defined windows.
    //
    // Implementation: piecewise-linear interpolation of the polyline through the current
    // images, re-sampled at equidistant arclength. (Linear rather than spline: with the
    // small per-step displacements of an MD the difference is negligible, and linear
    // interpolation cannot introduce the overshoot that a spline can - which would create
    // exactly the atom-atom clashes we fight elsewhere.)
    std::vector<Geometry> R = currentPath();

    // Claude Generated 2026: the polyline/arclength kernel now lives in
    // curcuma::string_path::reparametrizeEqualArclength (unweighted = the former inline code,
    // bit-identical, checked in test_string_path (d)); here only the application to the live
    // images remains.
    const std::vector<Geometry> targets = curcuma::string_path::reparametrizeEqualArclength(R, {});

    double disp2 = 0.0;
    int moved = 0;
    for (int i = 1; i < m_nimages - 1; ++i) {
        if (!m_images[i])
            continue; // clamped image: never moved

        const Geometry& target = targets[i];

        const double d = (target - R[i]).norm();
        if (d < 1e-10)
            continue;
        disp2 += d * d;
        ++moved;

        // Apply the new positions and restore the temperature: the coordinate change is
        // discontinuous, so the kinetic energy is preserved by rescaling rather than by
        // pretending the velocities still belong to the old geometry.
        const double T_before = m_images[i]->currentTemperature();
        m_images[i]->setPositions(target);
        if (T_before > 1e-6) {
            const double T_target = m_images[i]->targetTemperature();
            // Only correct gross drift; a mild mismatch is left to the thermostat.
            if (T_target > 1e-6 && T_before > 2.0 * T_target)
                m_images[i]->rescaleVelocities(std::sqrt(T_target / T_before));
        }
        m_image_geoms[i] = target;
    }

    return (moved > 0) ? std::sqrt(disp2 / moved) : 0.0;
}

/* ===================== Scheme 'string' (Claude Generated 2026) =====================
 * Finite-temperature string, stage 1 of the string -> path-CV chain. See nebmd.h for the
 * method; the pure geometry (Kabsch, averaging frames, reparametrisation) lives in
 * curcuma::string_path. */

namespace {
/*! Weighted squared norm sum_a w_a |d_a|^2 (uniform when w is empty). Claude Generated 2026. */
double stringWeightedSq(const Geometry& d, const std::vector<double>& w)
{
    if (w.empty())
        return d.squaredNorm();
    double s = 0.0;
    for (int a = 0; a < d.rows(); ++a)
        s += w[static_cast<size_t>(a)] * d.row(a).squaredNorm();
    return s;
}
} // namespace

void NebMD::prepareStringStartPath()
{
    // Claude Generated 2026: consecutive Kabsch alignment + equal-arclength reparametrisation of
    // the start path in the restraint metric, BEFORE the images are built, so every image starts
    // exactly on its prepared frame (energy/gradient consistent, no setPositions). Without it
    // update 1 measures and applies the reparametrisation of an unequally spaced / unaligned start
    // path (helicene ORCA MEP: 0.28 A, absorbed in 0.05 A trust-radius steps over 5 updates).
    namespace sp = curcuma::string_path;
    m_string_prep_shift = 0.0;
    if (m_image_geoms.size() < 3)
        return;
    const std::vector<double> w = m_mass_weighted_spring ? m_masses : std::vector<double>();
    std::vector<Geometry> P = sp::reparametrizeEqualArclength(sp::alignConsecutive(m_image_geoms, w), w);
    const double shift = sp::maxFrameDisplacement(P, m_image_geoms, w, true);
    m_string_prep_shift = shift;
    if (shift > 1e-6) {
        m_image_geoms = P;
        CurcumaLogger::result_fmt("NEB-MD string: start path prepared (aligned + equal arclength), "
                                  "largest best-fit shift {:.3e} A", shift);
    } else {
        CurcumaLogger::result_fmt("NEB-MD string: start path already aligned and equidistant "
                                  "(largest shift {:.2e} A <= 1e-6 A), left unchanged", shift);
    }
}

void NebMD::setupString()
{
    // Claude Generated 2026: freeze the starting mean path, choose k, open the history file.
    const int nat = static_cast<int>(m_image_geoms[0].rows());
    m_string_w = m_mass_weighted_spring ? m_masses : std::vector<double>();
    m_string_path = m_image_geoms; // the images start exactly on their frames (zero restraint)
    m_string_sum.assign(m_nimages, Geometry::Zero(nat, 3));
    m_string_n.assign(m_nimages, 0);
    m_string_capped.assign(m_nimages, 0);
    m_string_updates = 0;
    m_string_converged_count = 0;
    m_string_converged = false;
    m_string_last_delta = -1.0;
    m_string_aborted = false;
    m_string_last_scale = 1.0;
    m_string_block_sum.assign(m_nimages, std::vector<Geometry>(m_string_blocks, Geometry::Zero(nat, 3)));
    m_string_block_n.assign(m_nimages, std::vector<long>(m_string_blocks, 0));

    // Spacing d in the restraint metric (sqrt(amu)*A when mass-weighted) and the largest
    // parallel force along the starting path, from the step-0 energies.
    double dsum = 0.0;
    for (int i = 1; i < m_nimages; ++i)
        dsum += std::sqrt(stringWeightedSq(m_string_path[i] - m_string_path[i - 1], m_string_w));
    const double d = dsum / std::max(m_nimages - 1, 1);
    const std::vector<double> E = currentEnergies();
    double fpar = 0.0;
    for (int i = 1; i < m_nimages - 1; ++i)
        fpar = std::max(fpar, std::abs(E[i + 1] - E[i - 1]) / std::max(2.0 * d, 1e-12));
    if (m_string_k <= 0.0) {
        m_string_k = 4.0 * fpar / std::max(d, 1e-12);
        if (!(m_string_k > 0.0) || !std::isfinite(m_string_k)) {
            m_string_k = 0.05;
            CurcumaLogger::warn("NEB-MD string: automatic k undefined (flat start path); using 0.05.");
        } else {
            CurcumaLogger::result_fmt("NEB-MD string: k auto-set to {:.5f} from the start energies "
                                      "(F_par,max = {:.5f}, spacing d = {:.4f})", m_string_k, fpar, d);
        }
    }
    // Stability: restraint angular frequency sqrt(k) per MD time unit (1.9516 fs). The bound
    // k <= 0.2 at dt = 0.5 fs scales with 1/dt^2.
    double dt = 0.5;
    for (int i = 0; i < m_nimages; ++i)
        if (m_images[i]) {
            dt = m_images[i]->timeStep();
            break;
        }
    const double kmax = 0.2 * (0.5 / std::max(dt, 1e-6)) * (0.5 / std::max(dt, 1e-6));
    if (m_string_k > kmax)
        CurcumaLogger::warn(fmt::format("NEB-MD string: k = {:.4f} exceeds the stability guide {:.4f} for "
                                        "dt = {:.3f} fs (restraint period {:.1f} fs).",
                                        m_string_k, kmax, dt, 2.0 * M_PI * 1.9516 / std::sqrt(m_string_k)));
    CurcumaLogger::result_fmt("NEB-MD string: {} images, k = {:.4f} {} , update every {} steps, alpha = {:.2f}, "
                              "smooth = {:.2f}, tol = {:.3f} A x {}, max_time = {:.0f} fs",
                              m_nimages, m_string_k, m_mass_weighted_spring ? "Eh/(amu A^2)" : "Eh/A^2",
                              m_string_update, m_string_alpha, m_string_smooth, m_string_tol,
                              m_string_patience, m_string_max_time);
    CurcumaLogger::result_fmt("NEB-MD string: trust radius string_max_step = {:.3f} A{}, abort on dead image = {}",
                              m_string_max_step, m_string_max_step > 0.0 ? "" : " (off)",
                              m_string_abort_on_dead ? "yes" : "no");
    CurcumaLogger::result_fmt("NEB-MD string: holding estimate k*d/2 = {:.5f} vs F_par,max = {:.5f} "
                              "(start path, same metric)", 0.5 * m_string_k * d, fpar);
    if (m_string_ftol > 0.0)
        CurcumaLogger::result_fmt("NEB-MD string: force criterion ftol = {:.2e} Eh/A ({} blocks per update), "
                                  "gain = {:.1f} (z >= {:.1f}), start-path preparation {}",
                                  m_string_ftol, m_string_blocks, m_string_gain, m_string_gain_z,
                                  m_string_prepare ? fmt::format("shift {:.3e} A", m_string_prep_shift) : std::string("off"));
    else
        CurcumaLogger::warn("NEB-MD string: string_ftol = 0 - former criterion delta_max < string_tol "
                            "(false-positive prone, reproduction only)");

    std::ofstream h(outputPath(Basename() + ".neb.string_history.csv"));
    if (h)
        h << "update,step,delta_max,delta_rms,arclength,perp_max,below_tol,scale,n_dead,"
             "delta_perp_max,f_rms_max,f_se_max,f_sig_max,z_max,gain_max,prep_shift\n";
}

void NebMD::applyStringRestraint(int i, SimpleMD* img, Geometry& Fext) const
{
    // Claude Generated 2026: dV/dR_a = k w_a (R_a - M_aligned_a), M_i superimposed onto the
    // current R with the w-weighted Kabsch of the same sum (envelope theorem, see
    // test_string_path (b)). applyExternalForces() adds to the GRADIENT, so dV/dR is added
    // as is (Eh/A).
    if (!img || i < 0 || i >= static_cast<int>(m_string_path.size()))
        return;
    const Geometry& R = img->positions();
    const Geometry Ma = curcuma::string_path::applyAlignment(
        m_string_path[i], curcuma::string_path::kabsch(m_string_path[i], R, m_string_w));
    Geometry g(R.rows(), 3);
    double fmax = 0.0;
    for (int a = 0; a < R.rows(); ++a) {
        const double ka = m_string_k * (m_string_w.empty() ? 1.0 : m_string_w[static_cast<size_t>(a)]);
        g.row(a) = ka * (R.row(a) - Ma.row(a));
        fmax = std::max(fmax, g.row(a).norm());
    }
    // Same cap as the umbrella restraint (bias_fmax): uniform scaling keeps the direction.
    if (fmax > m_bias_fmax) {
        g *= m_bias_fmax / fmax;
        ++m_string_capped[i];
    }
    Fext += g;
}

void NebMD::accumulateStringMean(int i)
{
    // Claude Generated 2026: superimpose R_i onto M_i (same weights) and add it to the sum,
    // so the average lives in the frame of M_i and the path keeps one common frame.
    if (!m_images[i])
        return;
    const Geometry& R = m_images[i]->positions();
    const Geometry Ra = curcuma::string_path::applyAlignment(
        R, curcuma::string_path::kabsch(R, m_string_path[i], m_string_w));
    m_string_sum[i] += Ra;
    ++m_string_n[i];
    // Claude Generated 2026: block sums for the standard error of the mean force. Block index =
    // step within the update * blocks / string_update (stage-1 steps are counted from 1).
    if (i < static_cast<int>(m_string_block_sum.size()) && m_string_blocks > 0) {
        const long in_update = static_cast<long>((m_step - 1) % m_string_update);
        const int b = std::min(static_cast<int>(in_update * m_string_blocks / m_string_update), m_string_blocks - 1);
        m_string_block_sum[i][b] += Ra;
        ++m_string_block_n[i][b];
    }
}

bool NebMD::updateStringPath()
{
    // Claude Generated 2026: mean -> tau/perp -> block statistics (F, SE, z) -> gain -> mixing ->
    // smoothing -> reparametrisation -> delta/delta_perp -> trust radius, all in the pure
    // curcuma::string_path::updateMeanPath (unit-tested, test_string_path (i)); here only the
    // bookkeeping, the convergence criterion and the output. With string_gain <= 1 the new path
    // is bit-identical to the former inline update.
    namespace sp = curcuma::string_path;
    const int N = m_nimages;
    std::vector<Geometry> means(N);
    std::vector<std::vector<Geometry>> blocks(N);
    for (int i = 0; i < N; ++i) {
        if (m_string_n[i] <= 0)
            continue; // no samples (dead or clamped image): frame not updated
        means[i] = m_string_sum[i] / static_cast<double>(m_string_n[i]);
        if (i < static_cast<int>(m_string_block_sum.size()))
            for (int b = 0; b < m_string_blocks; ++b)
                if (m_string_block_n[i][b] > 0)
                    blocks[i].push_back(m_string_block_sum[i][b] / static_cast<double>(m_string_block_n[i][b]));
    }
    const sp::MeanPathUpdate U = sp::updateMeanPath(m_string_path, means, blocks, m_string_k, m_string_w,
                                                    m_string_alpha, m_string_gain, m_string_gain_z, m_string_ftol,
                                                    m_string_max_step, m_string_smooth, m_fixed_endpoints);
    const double delta_max = U.delta_max, scale = U.scale;
    m_string_last_scale = scale;
    int n_dead = 0;
    for (int i = 0; i < N; ++i)
        if (i < static_cast<int>(m_window_died.size()) && m_window_died[i] >= 0)
            ++n_dead;

    m_string_path = U.path;
    for (int i = 0; i < N; ++i) {
        m_string_sum[i].setZero();
        m_string_n[i] = 0;
        if (i < static_cast<int>(m_string_block_sum.size()))
            for (int b = 0; b < m_string_blocks; ++b) {
                m_string_block_sum[i][b].setZero();
                m_string_block_n[i][b] = 0;
            }
    }
    ++m_string_updates;
    if (m_string_write_updates) { // Claude Generated 2026: diagnostic history of M
        const std::vector<Geometry> frames = sp::alignConsecutive(m_string_path, {});
        const std::string upath = outputPath(Basename() + ".neb.string_updates.xyz");
        for (size_t i = 0; i < frames.size(); ++i) {
            Molecule mol = m_start;
            mol.setGeometry(frames[i]);
            mol.appendXYZFile(upath, fmt::format("string update={} image={}", m_string_updates, i));
        }
    }
    m_string_last_delta = delta_max;
    m_string_last_fsig = U.f_sig_max;
    m_string_last_se = U.f_se_max;
    m_string_last_z = U.z_max;
    // Claude Generated 2026: a scaled update (the path has not reached the proposed mean) or an
    // update with a frozen dead image (its frame is never updated) is no evidence of convergence.
    // Force criterion (string_ftol > 0): resolved (SE < ftol) and small (F_sig < ftol) perpendicular
    // mean force; delta_perp_max < string_tol is only a guard. string_ftol = 0: former criterion.
    const bool se_ok = std::isfinite(U.f_se_max) && U.f_se_max < m_string_ftol;
    bool below;
    if (m_string_ftol > 0.0)
        below = U.f_sig_max < m_string_ftol && se_ok && scale >= 1.0 && n_dead == 0
            && U.delta_perp_max < m_string_tol;
    else
        below = delta_max < m_string_tol && scale >= 1.0 && n_dead == 0;
    if (below)
        ++m_string_converged_count;
    else
        m_string_converged_count = 0;
    m_string_converged = (m_string_updates >= kStringMinUpdates && m_string_converged_count >= m_string_patience);

    auto fmtF = [](double v) { return std::isfinite(v) ? fmt::format("{:.6e}", v) : std::string("nan"); };
    std::ofstream h(outputPath(Basename() + ".neb.string_history.csv"), std::ios::app);
    if (h)
        h << m_string_updates << "," << m_step << "," << fmt::format("{:.6f}", delta_max) << ","
          << fmt::format("{:.6f}", U.delta_rms) << "," << fmt::format("{:.6f}", U.arclength) << ","
          << fmt::format("{:.6f}", U.perp_max) << "," << m_string_converged_count << ","
          << fmt::format("{:.6f}", scale) << "," << n_dead << ","
          << fmt::format("{:.6f}", U.delta_perp_max) << "," << fmtF(U.f_rms_max) << "," << fmtF(U.f_se_max) << ","
          << fmtF(U.f_sig_max) << "," << fmt::format("{:.3f}", U.z_max) << "," << fmt::format("{:.3f}", U.gain_max)
          << "," << fmt::format("{:.6f}", m_string_prep_shift) << "\n";
    CurcumaLogger::set_verbosity(m_band_verbosity);
    const std::string se_txt = std::isfinite(U.f_se_max) ? fmt::format("{:.2e}", U.f_se_max) : std::string("n/a");
    std::string noise;
    if (m_string_ftol > 0.0 && !se_ok)
        noise = std::isfinite(U.f_se_max)
            ? fmt::format(" [noise-limited: SE {:.2e} >= ftol {:.2e}; string_update x {:.1f} needed]", U.f_se_max,
                          m_string_ftol, (U.f_se_max / m_string_ftol) * (U.f_se_max / m_string_ftol))
            : std::string(" [noise-limited: SE n/a (< 2 blocks)]");
    CurcumaLogger::result_fmt("NEB-MD string update {:3d} (step {}): delta_max = {:.4f} A, delta_perp = {:.4f} A, "
                              "arclength = {:.4f} A, perp_max = {:.4f} A, F_sig = {:.2e} Eh/A, SE = {}, z = {:.1f}, "
                              "gain = {:.1f}, scale = {:.3f}, live images = {}/{}{}{}",
                              m_string_updates, m_step, delta_max, U.delta_perp_max, U.arclength, U.perp_max,
                              U.f_sig_max, se_txt, U.z_max, U.gain_max, scale, N - n_dead, N, noise,
                              m_string_converged ? "  -> converged" : "");
    return m_string_converged;
}

void NebMD::writeStringPath(const std::string& filename) const
{
    // Claude Generated 2026: consecutive UNWEIGHTED Kabsch alignment, so the plain frame
    // differences used by the stage-2 consumers (PathCV lambda, path spacing, Z floor) equal
    // the unweighted best-fit RMSD the path CV itself uses. Atom order = start structure.
    const std::vector<Geometry> frames = curcuma::string_path::alignConsecutive(m_string_path, {});
    const std::string path = outputPath(filename);
    std::remove(path.c_str());
    EnergyCalculator ec(m_method, m_controller, Basename());
    Molecule probe = m_start;
    ec.setMolecule(probe.getMolInfo());
    for (size_t i = 0; i < frames.size(); ++i) {
        const double E = bandEnergy(ec, frames[i], nullptr);
        Molecule mol = m_start;
        mol.setGeometry(frames[i]);
        mol.appendXYZFile(path, fmt::format("string update={} image={} E={:.8f}", m_string_updates, i, E));
    }
    CurcumaLogger::set_verbosity(m_band_verbosity);
}

void NebMD::histogramTrajectory(const std::string& traj) const
{
    // Claude Generated 2026: G_ref(s) = -kT ln P(s) of a free trajectory on the path CV. The
    // CV is the same object stage 2 builds (PathCV on m_image_geoms, lambda from the path),
    // so G_ref and the stage-2 mean-force G(s) refer to one coordinate.
    CitationRegistry::cite("pathcv");
    std::unique_ptr<PathCV> local;
    const PathCV* cv = m_path_cv.get();
    const std::vector<Geometry>& cvp = cvFrames(); // Claude Generated 2026: same CV as stage 2
    if (!cv) {
        local = makePathCV(cvp); // Claude Generated 2026: rmsd or distance metric
        cv = local.get();
    }
    const int nat = static_cast<int>(m_start.AtomCount());
    const double T = m_md_config.contains("temperature") && m_md_config["temperature"].is_number()
        ? m_md_config["temperature"].get<double>() : 300.0;
    const double kT = 3.166811563e-6 * T;
    const double kcal = 627.5094740631;

    // reference frames: interior range and the on-path z baseline z_path(s). z is offset
    // below zero ON the path, and the offset depends on the position (an end frame has one
    // neighbour, an interior frame two), so the filter works on the excess z - z_path(s) with
    // z_path linearly interpolated between the frames (sorted by s, clamped at the ends).
    // Interior range = window range of the CV frames [C_lo, C_hi] (= first/last frame by default).
    const int cv_lo = m_cv_path.empty() ? 0 : m_cv_lo;
    const int cv_hi = m_cv_path.empty() ? static_cast<int>(cvp.size()) - 1 : m_cv_hi;
    const double s_first = cv->s(cvp[cv_lo]), s_last = cv->s(cvp[cv_hi]);
    const double s_lo = std::min(s_first, s_last), s_hi = std::max(s_first, s_last);
    std::vector<std::pair<double, double>> zpath;
    for (const auto& g : cvp) {
        double s = 0.0, z = 0.0;
        cv->evaluate(g, &s, &z, nullptr, nullptr);
        zpath.emplace_back(s, z);
    }
    std::sort(zpath.begin(), zpath.end());
    auto zBase = [&zpath](double s) {
        if (s <= zpath.front().first)
            return zpath.front().second;
        if (s >= zpath.back().first)
            return zpath.back().second;
        size_t k = 1;
        while (k < zpath.size() - 1 && zpath[k].first < s)
            ++k;
        const double ds = zpath[k].first - zpath[k - 1].first;
        const double t = ds > 1e-15 ? (s - zpath[k - 1].first) / ds : 0.0;
        return (1.0 - t) * zpath[k - 1].second + t * zpath[k].second;
    };

    std::vector<double> sv, zv;
    std::vector<long> frame_index; // Claude Generated 2026: file position of each usable frame (TST pass)
    long wrong_nat = 0, iframe = -1;
    FileIterator it(traj);
    while (!it.AtEnd()) {
        Molecule m = it.Next();
        ++iframe;
        if (static_cast<int>(m.AtomCount()) != nat) {
            ++wrong_nat;
            continue;
        }
        double s = 0.0, z = 0.0;
        cv->evaluate(m.getGeometry(), &s, &z, nullptr, nullptr);
        if (!std::isfinite(s) || !std::isfinite(z)) {
            ++wrong_nat;
            continue;
        }
        sv.push_back(s);
        zv.push_back(z);
        frame_index.push_back(iframe);
    }
    const long nframes = static_cast<long>(sv.size());
    if (nframes == 0) {
        CurcumaLogger::error("NEB-MD cv_histogram: no usable frame in '" + traj + "'.");
        return;
    }
    // filter quantity: absolute z (user threshold) or the excess z - z_path(s) (automatic).
    // Claude Generated 2026: the automatic threshold separates the on-path population from a
    // distant one (frames in another basin, e.g. permuted labels) at the density valley of
    // log10(excess); unimodal -> only outliers > 10 x median excess (see z_filter.h).
    const bool auto_z = (m_cv_hist_zmax <= 0.0);
    std::vector<double> zf(nframes);
    double zmean = 0.0;
    for (long k = 0; k < nframes; ++k) {
        zf[k] = auto_z ? zv[k] - zBase(sv[k]) : zv[k];
        zmean += zf[k];
    }
    zmean /= static_cast<double>(nframes);
    curcuma::z_filter::ZThreshold zsel;
    if (auto_z)
        zsel = curcuma::z_filter::chooseThreshold(zf);
    else
        zsel.rule = "user cv_hist_zmax";
    const double zmax = auto_z ? zsel.threshold : m_cv_hist_zmax;
    long n_interior = 0, rej_interior = 0;
    for (long k = 0; k < nframes; ++k) {
        if (sv[k] < s_lo || sv[k] > s_hi)
            continue;
        ++n_interior;
        if (zf[k] > zmax)
            ++rej_interior;
    }

    const int B = m_fe_bins;
    const int nblocks = 5;
    std::vector<long> count(B, 0);
    std::vector<std::vector<long>> bcount(nblocks, std::vector<long>(B, 0));
    std::vector<long> btotal(nblocks, 0);
    long rejected = 0, used = 0;
    for (long k = 0; k < nframes; ++k) {
        if (zf[k] > zmax) {
            ++rejected;
            continue;
        }
        int b = static_cast<int>(std::floor(sv[k] * B));
        b = std::min(std::max(b, 0), B - 1);
        const int blk = static_cast<int>(std::min<long>(nblocks - 1, k * nblocks / nframes));
        ++count[b];
        ++bcount[blk][b];
        ++btotal[blk];
        ++used;
    }

    auto gOf = [&](long c, long tot) {
        return (c > 0 && tot > 0) ? -kT * std::log(static_cast<double>(c) / static_cast<double>(tot)) * kcal
                                  : std::numeric_limits<double>::quiet_NaN();
    };
    std::vector<double> G(B), Gse(B, std::numeric_limits<double>::quiet_NaN());
    double gmin = std::numeric_limits<double>::infinity();
    for (int b = 0; b < B; ++b) {
        G[b] = gOf(count[b], used);
        if (std::isfinite(G[b]))
            gmin = std::min(gmin, G[b]);
    }
    for (int b = 0; b < B; ++b) {
        if (std::isfinite(G[b]))
            G[b] -= gmin;
        // block SE: only where every block populated the bin (otherwise undefined)
        double m = 0.0, m2 = 0.0;
        bool all = true;
        for (int j = 0; j < nblocks; ++j) {
            const double g = gOf(bcount[j][b], btotal[j]);
            if (!std::isfinite(g)) {
                all = false;
                break;
            }
            m += g;
            m2 += g * g;
        }
        if (all) {
            m /= nblocks;
            Gse[b] = std::sqrt(std::max(m2 / nblocks - m * m, 0.0) / (nblocks - 1));
        }
    }

    // barrier / reaction dG over populated bins inside the interior range (AP4 convention)
    int b0 = -1, b1 = -1;
    double gmax_in = -std::numeric_limits<double>::infinity();
    for (int b = 0; b < B; ++b) {
        const double sc = (b + 0.5) / B;
        if (!std::isfinite(G[b]) || sc < s_lo || sc > s_hi)
            continue;
        if (b0 < 0) b0 = b;
        b1 = b;
        gmax_in = std::max(gmax_in, G[b]);
    }

    // Claude Generated 2026: basin-referenced barrier of G_ref, same function as the stage-2
    // profile (fe_profile.h basinBarrier) on the bin midpoints of the populated (count > 0, i.e.
    // after the z filter) bins inside the interior range; coverage from finite-difference slopes.
    std::vector<double> xb, gb;
    for (int b = 0; b < B; ++b) {
        const double sc = (b + 0.5) / B;
        if (count[b] <= 0 || !std::isfinite(G[b]) || sc < s_lo || sc > s_hi)
            continue;
        xb.push_back(sc);
        gb.push_back(G[b]);
    }
    curcuma::fe_profile::BasinBarrier basR;
    curcuma::fe_profile::BasinCoverage covR;
    if (xb.size() >= 3) {
        basR = curcuma::fe_profile::basinBarrier(xb, gb, kT * kcal);
        covR = curcuma::fe_profile::basinCoverage(xb, curcuma::fe_profile::finiteDifferenceSlope(xb, gb));
    } else {
        basR.warning = "fewer than 3 populated interior bins";
    }
    const bool covR_bad = xb.size() >= 3 && (!covR.left || !covR.right);

    // Claude Generated 2026: absolute TST rate of the free trajectory (fe_profile.h tstRate).
    // <Z^{1/2}> at s* from the frames themselves (free MD: weight 1, plain mean of Z^{1/2}), only
    // frames that pass the z filter, in the bins within s* +- 1.5 bin widths; if those hold fewer
    // than 50 frames, from the nearest populated bins on each side (zHalfNear fallback). Z needs
    // ds/dR, so the trajectory is read a second time and only the frames in the needed bins are
    // evaluated with the gradient. p(s*) from G_ref on the same bin midpoints as the basin barrier.
    curcuma::fe_profile::TstRate tstR;
    curcuma::fe_profile::ZHalfEstimate zhR;
    std::string tstR_note;
    if (!basR.valid) {
        tstR_note = "not computed - no basin barrier (" + basR.warning + ")";
    } else {
        const double hw = 1.5 / B;
        const long zmin_samples = 50;
        std::vector<char> need(B, 0);
        long nwin = 0;
        for (int b = 0; b < B; ++b)
            if (count[b] > 0 && std::abs((b + 0.5) / B - basR.s_star) <= hw) {
                need[b] = 1;
                nwin += count[b];
            }
        if (nwin < zmin_samples) { // the same bins zHalfNear's fallback will use (count = z-filtered frames)
            std::vector<double> sc(B);
            for (int b = 0; b < B; ++b)
                sc[b] = (b + 0.5) / B;
            for (size_t b : curcuma::fe_profile::zHalfFallbackBins(sc, count, basR.s_star, zmin_samples))
                need[b] = 1;
        }
        std::vector<double> zsum(B, 0.0);
        std::vector<long> zcnt(B, 0);
        // usable frame k (first pass) -> its bin, if it is needed and passes the z filter
        std::vector<int> frame_bin(nframes, -1);
        for (long k = 0; k < nframes; ++k) {
            const int b = std::min(std::max(static_cast<int>(std::floor(sv[k] * B)), 0), B - 1);
            if (need[b] && !(zf[k] > zmax))
                frame_bin[k] = b;
        }
        FileIterator it2(traj);
        long iframe2 = -1, k = 0;
        while (!it2.AtEnd() && k < nframes) {
            Molecule m = it2.Next();
            ++iframe2;
            if (iframe2 != frame_index[k])
                continue;
            const int b = frame_bin[k++];
            if (b < 0)
                continue;
            double s = 0.0, z = 0.0;
            Geometry ds;
            cv->evaluate(m.getGeometry(), &s, &z, &ds, nullptr);
            double Z = 0.0;
            for (int a = 0; a < nat; ++a)
                Z += ds.row(a).squaredNorm() / std::max(m_masses[a], 1e-12);
            if (std::isfinite(Z) && Z > 0.0) {
                zsum[b] += std::sqrt(Z);
                ++zcnt[b];
            }
        }
        std::vector<double> zs, zv;
        std::vector<long> zn;
        for (int b = 0; b < B; ++b)
            if (zcnt[b] > 0) {
                zs.push_back((b + 0.5) / B);
                zn.push_back(zcnt[b]);
                zv.push_back(zsum[b] / zcnt[b]);
            }
        zhR = curcuma::fe_profile::zHalfNear(zs, zn, zv, basR.s_star, hw, zmin_samples);
        if (!zhR.valid)
            tstR_note = "not computed - " + zhR.warning;
        else {
            std::vector<double> gb_eh(gb.size());
            for (size_t j = 0; j < gb.size(); ++j)
                gb_eh[j] = gb[j] / kcal;
            tstR = curcuma::fe_profile::tstRate(xb, gb_eh, kT, basR.s_star, zhR.value);
            if (!tstR.valid)
                tstR_note = "not computed - " + tstR.warning;
        }
    }

    const std::string out = outputPath(Basename() + ".neb.ref.csv");
    std::ofstream f(out);
    if (f) {
        f << "# Reference G_ref(s) = -kT ln P(s) of a free trajectory on the path CV (Claude Generated 2026)\n";
        f << "# trajectory " << traj << ", frames " << nframes << " (skipped wrong atom count / non-finite "
          << wrong_nat << "), T = " << T << " K, lambda = " << cv->lambda() << "\n";
        f << fmt::format("# z filter: {} <= {:.6g} A^2 ({}), mean of the filtered quantity {:.6g}, rejected {} = {:.4f}\n",
                         auto_z ? "z - z_path(s)" : "z", zmax, zsel.rule,
                         zmean, rejected, static_cast<double>(rejected) / static_cast<double>(nframes));
        f << fmt::format("# z filter interior: {} of {} interior frames rejected = {:.4f}\n", rej_interior, n_interior,
                         n_interior > 0 ? static_cast<double>(rej_interior) / static_cast<double>(n_interior) : 0.0);
        if (auto_z)
            f << fmt::format("# z filter detail: median excess {:.6g}, on-path mode {:.6g}, far peak {:.6g}, "
                             "peak score {:.3g}, bimodal {}\n",
                             zsel.median, zsel.mode, zsel.far_peak, zsel.score, zsel.bimodal ? 1 : 0);
        f << fmt::format("# interior s range [{:.6f}, {:.6f}] = s(first frame), s(last frame); block SE from {} contiguous blocks\n",
                         s_lo, s_hi, nblocks);
        if (basR.valid) {
            f << fmt::format("# G_ref barrier (basin, fwd/rev): {:.4f} / {:.4f} kcal/mol at s* = {:.4f} (G_X = -kT ln "
                             "int_X exp(-G/kT) ds over the populated interior bin midpoints, X = s <= s* / s >= s*; "
                             "ds in units of the CV coordinate s)\n", basR.barrier_fwd, basR.barrier_rev, basR.s_star);
            f << fmt::format("# G_ref reaction dG (basin): {:.4f} kcal/mol\n", basR.reaction);
        } else {
            f << "# G_ref barrier (basin, fwd/rev): not defined - " << basR.warning << "\n";
        }
        // Claude Generated 2026: absolute TST rate of the free trajectory
        if (tstR.valid) {
            f << fmt::format("# G_ref TST rate (fwd/rev): {:.6e} / {:.6e} 1/ps at s* = {:.4f} (k = 1/2 <|sdot|>(s*) p_X(s*), "
                             "p_X = exp(-G_ref(s*)/kT) / int_X exp(-G_ref/kT) ds over the populated interior bin midpoints; "
                             "independent of the scale of s)\n", tstR.k_fwd, tstR.k_rev, tstR.s_star);
            f << fmt::format("# G_ref TST dG (Eyring, fwd/rev): {:.4f} / {:.4f} kcal/mol (kT ln((kB T/h)/k), kB T/h = {:.6f} 1/ps)\n",
                             tstR.dG_fwd * kcal, tstR.dG_rev * kcal, tstR.eyring_prefactor);
            f << fmt::format("# G_ref TST detail: <|sdot|>(s*) {:.6e} s/ps = sqrt(2kT/pi) <Z^1/2>, <Z^1/2> {:.6e} amu^-1/2 A^-1 "
                             "(mean over z-filtered frames) from {} frames in {} bin(s){}; p(s*) fwd/rev {:.6e} / {:.6e} 1/s\n",
                             tstR.sdot_abs, zhR.value, zhR.samples, zhR.bins,
                             zhR.fallback ? " (FALLBACK: " + zhR.warning + ")" : "", tstR.p_fwd, tstR.p_rev);
        } else {
            f << "# G_ref TST rate: " << tstR_note << "\n";
        }
        if (covR_bad)
            f << fmt::format("# basin not covered: left {} right {} (dG/ds first {:.4e}, last {:.4e}, max|dG/ds| "
                             "{:.4e} kcal/mol)\n", covR.left ? "ok" : "NOT", covR.right ? "ok" : "NOT",
                             covR.slope_left, covR.slope_right, covR.max_abs);
        if (!m_cv_path.empty() && m_cv_extend_length > 0.0) // Claude Generated 2026 (only in length mode)
            f << fmt::format("# cv_extend_length {:.6g} of the path length = {} frames per end\n",
                             m_cv_extend_length, m_cv_extend);
        if (!m_cv_atoms.empty()) // Claude Generated 2026 (only with a selection: default header unchanged)
            f << fmt::format("# cv_atoms {} ({} atoms)\n", m_cv_atoms_spec, m_cv_atoms.size());
        if (m_cv_metric == PathCV::Metric::Distances) // Claude Generated 2026 (default header unchanged)
            f << "# cv_metric distances, " << cvPairsText() << "\n";
        if (!m_cv_path.empty() && m_cv_orig_s.size() == 4)
            f << fmt::format("# cv path: {} frames (cv_frames {}, cv_extend {}), lambda {:.6g}; original path frames "
                             "on this CV: s(first) {:.6f}, s(second) {:.6f}, s(second-to-last) {:.6f}, s(last) {:.6f}\n",
                             m_cv_path.size(), m_cv_frames, m_cv_extend, cv->lambda(),
                             m_cv_orig_s[0], m_cv_orig_s[1], m_cv_orig_s[2], m_cv_orig_s[3]);
        f << "s,count,P,G_kcal,G_se_kcal,interior\n";
        for (int b = 0; b < B; ++b) {
            const double sc = (b + 0.5) / B;
            f << fmt::format("{:.4f}", sc) << "," << count[b] << ","
              << fmt::format("{:.6e}", used > 0 ? static_cast<double>(count[b]) / used : 0.0) << ","
              << (std::isfinite(G[b]) ? fmt::format("{:.4f}", G[b]) : "nan") << ","
              << (std::isfinite(Gse[b]) ? fmt::format("{:.4f}", Gse[b]) : "nan") << ","
              << ((sc >= s_lo && sc <= s_hi) ? 1 : 0) << "\n";
        }
    }

    CurcumaLogger::set_verbosity(m_band_verbosity);
    CurcumaLogger::result("");
    CurcumaLogger::result("=== Reference histogram of the path CV (no MD) ===");
    CurcumaLogger::result_fmt("  trajectory         : {} ({} frames used of {}, {} skipped)", traj, used, nframes, wrong_nat);
    CurcumaLogger::result_fmt("  z filter           : {} <= {:.4g} A^2 ({}), rejected {} = {:.2f} %",
                              auto_z ? "z - z_path(s)" : "z", zmax, zsel.rule, rejected,
                              100.0 * rejected / std::max<long>(nframes, 1));
    CurcumaLogger::result_fmt("  z filter interior  : {} of {} interior frames rejected = {:.2f} %", rej_interior,
                              n_interior, 100.0 * rej_interior / std::max<long>(n_interior, 1));
    if (auto_z)
        CurcumaLogger::result_fmt("  z filter detail    : median excess {:.4g}, on-path mode {:.4g}, far peak {:.4g}, "
                                  "score {:.3g}", zsel.median, zsel.mode, zsel.far_peak, zsel.score);
    CurcumaLogger::result_fmt("  interior range     : s = [{:.4f}, {:.4f}]", s_lo, s_hi);
    if (b0 >= 0)
        CurcumaLogger::result_fmt("  G_ref barrier      : {:.2f} kcal/mol, reaction dG {:.2f} (first/last populated interior bin "
                                  "s = {:.3f} / {:.3f})",
                                  gmax_in - G[b0], G[b1] - G[b0], (b0 + 0.5) / B, (b1 + 0.5) / B);
    if (basR.valid) {
        CurcumaLogger::result_fmt("  G_ref barrier (basin, fwd/rev): {:.2f} / {:.2f} kcal/mol at s* = {:.4f}",
                                  basR.barrier_fwd, basR.barrier_rev, basR.s_star);
        CurcumaLogger::result_fmt("  G_ref reaction dG (basin): {:.2f} kcal/mol", basR.reaction);
    } else {
        CurcumaLogger::warn("NEB-MD cv_histogram: G_ref barrier (basin, fwd/rev) not defined - " + basR.warning);
    }
    if (tstR.valid) { // Claude Generated 2026
        CurcumaLogger::result_fmt("  G_ref TST rate (fwd/rev): {:.4e} / {:.4e} 1/ps at s* = {:.4f}",
                                  tstR.k_fwd, tstR.k_rev, tstR.s_star);
        CurcumaLogger::result_fmt("  G_ref TST dG (Eyring, fwd/rev): {:.2f} / {:.2f} kcal/mol", tstR.dG_fwd * kcal,
                                  tstR.dG_rev * kcal);
        CurcumaLogger::result_fmt("  G_ref TST detail   : <|sdot|>(s*) {:.4e} s/ps, p(s*) fwd/rev {:.4e} / {:.4e} 1/s, "
                                  "<Z^1/2> {:.4e} from {} frames in {} bin(s){}",
                                  tstR.sdot_abs, tstR.p_fwd, tstR.p_rev, zhR.value, zhR.samples, zhR.bins,
                                  zhR.fallback ? " (fallback)" : "");
        if (zhR.fallback)
            CurcumaLogger::warn("NEB-MD cv_histogram TST rate: " + zhR.warning);
    } else {
        CurcumaLogger::result("  G_ref TST rate     : " + tstR_note);
    }
    if (covR_bad)
        CurcumaLogger::warn(fmt::format("basin not covered ({}{}{}): the slope of G_ref at the outermost populated "
                                        "interior bin is still > 0.5 max|dG/ds| (left {:.3f}, right {:.3f}, max {:.3f} "
                                        "kcal/mol per unit s).", covR.left ? "" : "left",
                                        (!covR.left && !covR.right) ? ", " : "", covR.right ? "" : "right",
                                        covR.slope_left, covR.slope_right, covR.max_abs));
    CurcumaLogger::result_fmt("  written            : {}", out);
}


double NebMD::windowK(int i) const
{
    // Claude Generated 2026: restraint stiffness actually applied to window i.
    return (i >= 0 && i < static_cast<int>(m_umbrella_kvec.size())) ? m_umbrella_kvec[i] : m_umbrella_k;
}

double NebMD::thermalEnergy() const
{
    // kT in Hartree at the thermostat setpoint. Claude Generated 2026.
    // Taken from the images rather than hard-coded, so a run at a different temperature
    // does not silently keep using 300 K in the metric term.
    double T = 300.0;
    for (int i = 0; i < m_nimages; ++i)
        if (m_images[i]) {
            const double t = m_images[i]->targetTemperature();
            if (t > 1e-6)
                T = t;
            break;
        }
    return 3.166811563e-6 * T;
}

Geometry NebMD::pathBiasField(const Geometry& R, double* Z_out) const
{
    // Claude Generated 2026: the mass-weighted field
    //
    //     b_a = grad_a s / (m_a Z) ,      Z = sum_a |grad_a s|^2 / m_a
    //
    // which satisfies b.grad s = 1 by construction (the only requirement the master
    // equation makes) and whose Z is the inverse effective mass of the coordinate.
    Geometry b = Geometry::Zero(R.rows(), 3);
    if (!m_path_cv)
        return b;
    double s = 0.0, z = 0.0;
    Geometry ds;
    m_path_cv->evaluate(R, &s, &z, &ds, nullptr);
    const int nat = static_cast<int>(R.rows());
    double Z = 0.0;
    for (int a = 0; a < nat; ++a)
        Z += ds.row(a).squaredNorm() / std::max(m_masses[a], 1e-12);
    if (Z_out)
        *Z_out = Z;
    if (Z < 1e-30)
        return b;
    for (int a = 0; a < nat; ++a)
        b.row(a) = ds.row(a) / (std::max(m_masses[a], 1e-12) * Z);
    return b;
}

void NebMD::accumulateMeanForce(int /*i*/, SimpleMD* img, FeLocal& out, std::mt19937& rng) const
{
    // Claude Generated 2026: conditional mean force for the decoupled scheme.
    //
    //     dG/ds = <b.grad E>_s - kT <div b>_s
    //
    // Samples from every window go into the SAME global bins: a restraint that depends
    // only on s is constant on the surface s = const and therefore does not change the
    // conditional distribution there. No unbiasing, no window matching, no overlap
    // requirement - the windows only have to make sure every s is visited at all.
    if (!m_path_cv || !img)
        return;
    const Geometry R = img->positions();
    const int nat = static_cast<int>(R.rows());

    double s = 0.0, z = 0.0;
    Geometry ds;
    m_path_cv->evaluate(R, &s, &z, &ds, nullptr);
    double Z = 0.0;
    for (int a = 0; a < nat; ++a)
        Z += ds.row(a).squaredNorm() / std::max(m_masses[a], 1e-12);

    // The field b scales as 1/|grad s|, so it diverges wherever the coordinate stops
    // responding to the geometry. That is not a numerical artefact but the honest
    // statement that dG/ds is undefined there - and it would swamp the average. Such
    // samples are therefore discarded and counted; a large discarded fraction means the
    // coordinate is unusable for this system, which is a diagnosis, not a defect.
    if (Z < m_fe_Zfloor) {
        ++out.rejected;
        // Claude Generated 2026: the bin is known from s even for a rejected sample, so the
        // rejection fraction can be reported per support point (rejected_frac in fe.csv).
        if (std::isfinite(s) && !out.rejected_bin.empty()) {
            int rb = static_cast<int>(std::floor(s * m_fe_bins));
            rb = std::min(std::max(rb, 0), m_fe_bins - 1);
            ++out.rejected_bin[rb];
        }
        return;
    }

    // <b.grad E>: the true gradient, projected with the mass-weighted field.
    // With a z wall (z_wall > 0) trueGradient() is d(E + W)/dR: W is part of the sampled
    // potential (SimpleMD::setExtraPotential), so the estimator is that of the tube ensemble.
    const Geometry g = img->trueGradient(); // dE/dR in Eh/Angstrom
    double bdotg = 0.0;
    for (int a = 0; a < nat; ++a)
        bdotg += ds.row(a).dot(g.row(a)) / (std::max(m_masses[a], 1e-12) * Z);
    // Claude Generated 2026: the wall part b.grad W on its own (same integrator step as g), so
    // fe.csv can split the mean force into <b.grad E> + <b.grad W> - kT <div b>.
    const Geometry& gw = img->extraPotentialGradient();
    double bdotgw = 0.0;
    if (gw.rows() == nat)
        for (int a = 0; a < nat; ++a)
            bdotgw += ds.row(a).dot(gw.row(a)) / (std::max(m_masses[a], 1e-12) * Z);

    // <div b>: needs second derivatives of s. They are obtained from finite differences
    // of the ANALYTIC gradient along random directions (Hutchinson):
    //
    //     div b = tr(db/dR) = E_v [ v . (db/dR) v ] ,  v ~ N(0,1)
    //
    // so one probe costs two extra evaluations of ds/dR and NO energy evaluation. The
    // single-probe variance is large, but it averages out over the trajectory - which is
    // exactly what the estimator does anyway.
    double divb = 0.0;
    if (m_div_probes > 0) {
        std::normal_distribution<double> gauss(0.0, 1.0);
        const double eps = m_div_epsilon;
        for (int p = 0; p < m_div_probes; ++p) {
            Geometry v(nat, 3);
            for (int a = 0; a < nat; ++a)
                for (int c = 0; c < 3; ++c)
                    v(a, c) = gauss(rng);
            const Geometry bp = pathBiasField(R + eps * v);
            const Geometry bm = pathBiasField(R - eps * v);
            double acc = 0.0;
            for (int a = 0; a < nat; ++a)
                acc += v.row(a).dot(bp.row(a) - bm.row(a));
            divb += acc / (2.0 * eps);
        }
        divb /= static_cast<double>(m_div_probes);
    }

    // A window whose MD has gone unstable produces non-finite values. Without this guard
    // they propagate silently into the running sums and every reported number becomes NaN
    // - which is worse than a crash, because the run still finishes and prints a profile.
    // Measured: 4 of 8 helicene windows died this way and poisoned both estimators.
    if (!std::isfinite(s) || !std::isfinite(bdotg) || !std::isfinite(divb)) {
        ++out.nonfinite;
        return;
    }

    int bin = static_cast<int>(std::floor(s * m_fe_bins));
    bin = std::min(std::max(bin, 0), m_fe_bins - 1);
    // Ensemble weight. In the CONSTRAINED ensemble the sampled distribution differs from
    // the marginal one by Z^(-1/2), so every conditional average has to be reweighted:
    //     dG/ds = <Z^-1/2 (b.grad E - kT div b)> / <Z^-1/2>
    // For the restrained (unconstrained) scheme the weight is one and the expression
    // reduces to the plain average - both routes therefore share one code path.
    const double w = m_constrain ? 1.0 / std::sqrt(Z) : 1.0;
    out.count[bin] += 1;
    out.wsum[bin] += w;
    out.force[bin] += w * bdotg;
    out.div[bin] += w * divb;
    out.wall[bin] += w * bdotgw;
    out.ssum[bin] += s;

    {
        out.win_n += 1;
        out.win_s += s;
        out.win_s2 += s * s;
        // Thermal amplitude for the conditioning check. Claude Generated 2026.
        //
        // Measured through z, NOT through the Cartesian variance: the images translate and
        // rotate freely, so a variance of the raw coordinates is dominated by rigid-body
        // motion (measured 17 A for ethane, which is obviously not a vibration). z is a
        // best-fit RMSD distance from the path and therefore invariant under both, so
        // sqrt(<z>) is the transverse amplitude in the same metric as the path spacing.
        out.win_z += z;
        // Claude Generated 2026: z wall diagnostics (fraction outside the tube, wall energy).
        out.win_z_max = std::max(out.win_z_max, z);
        if (m_z_wall > 0.0 && z > m_z_wall) {
            ++out.win_zwall_n;
            out.win_wall += 0.5 * m_z_wall_k * (z - m_z_wall) * (z - m_z_wall);
        }
        out.win_w += w;
        out.win_wf += w * bdotg;
        out.win_wd += w * divb;
        out.win_wwall += w * bdotgw;
    }
}

void NebMD::writeMeanForceProfile() const
{
    CurcumaLogger::set_verbosity(m_band_verbosity); // see m_band_verbosity in nebmd.h
    // Claude Generated 2026: two independent estimators of the same profile.
    //
    //   (I)  binned conditional mean force   dG/ds = <b.grad E> - kT <div b>
    //   (II) mean restraint force per window dG/ds(<s>) = -k (<s> - s0)
    //
    // (II) carries the metric contribution implicitly, because the sampled distribution
    // of s IS the marginal density; (I) computes it explicitly through div b. Agreement
    // of the two is therefore the first internal check of the metric term available in
    // this code - previously it could only be judged by its order of magnitude.
    CitationRegistry::cite("denotter");
    const double kT = thermalEnergy();
    const double kcal = 627.5094740631;

    long total = 0;
    for (int b = 0; b < m_fe_bins; ++b)
        total += m_fe_count[b];
    if (total <= 0) {
        CurcumaLogger::warn("NEB-MD: no mean-force samples collected - run longer than "
                            "pmf_equilibration.");
        return;
    }

    // ---------- estimator (I): conditional mean force ----------
    std::vector<double> sI, dgI;
    std::vector<long> nI, rejI; // samples and Z-floor rejections per support point (Claude Generated 2026)
    std::vector<fe_profile::MeanForceSplit> splitI; // <b.grad E>, -kT<div b>, <b.grad W> (Claude Generated 2026)
    if (m_constrain) {
        // Constrained scheme: one support point per window, at the constrained value of s.
        // Binning would be meaningless here - s never moves, so a window fills exactly one
        // bin. This is the blue-moon form: <Z^-1/2 (b.grad E - kT div b)> / <Z^-1/2>.
        std::vector<std::pair<std::pair<double, double>, int>> pts;
        for (int i = 0; i < m_nimages; ++i) {
            if (m_win_count[i] < 10 || m_win_w[i] < 1e-300)
                continue;
            pts.push_back({ { m_umbrella_s0[i], (m_win_wf[i] - kT * m_win_wd[i]) / m_win_w[i] }, i });
        }
        std::sort(pts.begin(), pts.end());
        for (const auto& p : pts) {
            sI.push_back(p.first.first);
            dgI.push_back(p.first.second);
            nI.push_back(m_win_count[p.second]);
            rejI.push_back(p.second < static_cast<int>(m_win_rejected.size()) ? m_win_rejected[p.second] : 0);
            const int i = p.second;
            splitI.push_back(fe_profile::splitMeanForce(m_win_wf[i], m_win_wwall[i], m_win_wd[i], m_win_w[i], kT));
        }
        if (pts.size() < 3)
            CurcumaLogger::warn("Fewer than three constrained windows contributed; a profile "
                                "integrated over so few support points carries a large "
                                "discretisation error and the barrier should not be quoted.");
    } else {
        for (int b = 0; b < m_fe_bins; ++b) {
            if (m_fe_count[b] < 10) // an almost empty bin is noise, not information
                continue;
            const double n = static_cast<double>(m_fe_count[b]);
            const double wn = (m_fe_wsum[b] > 1e-300) ? m_fe_wsum[b] : n;
            sI.push_back(m_fe_ssum[b] / n);
            dgI.push_back(m_fe_force[b] / wn - kT * m_fe_div[b] / wn);
            nI.push_back(m_fe_count[b]);
            rejI.push_back(b < static_cast<int>(m_fe_rejected_bin.size()) ? m_fe_rejected_bin[b] : 0);
            splitI.push_back(fe_profile::splitMeanForce(m_fe_force[b], m_fe_wall[b], m_fe_div[b], wn, kT));
        }
    }
    auto integrate = [](const std::vector<double>& x, const std::vector<double>& y) {
        std::vector<double> G(x.size(), 0.0);
        for (size_t k = 1; k < x.size(); ++k)
            G[k] = G[k - 1] + 0.5 * (y[k] + y[k - 1]) * (x[k] - x[k - 1]);
        return G;
    };
    // Claude Generated 2026: barrier and reaction dG only from support points inside the
    // interior range [s(frame 0), s(frame N-1)] of the reference path. Outside it the path CV
    // saturates (no frame beyond the ends), so G there is not a free energy along the path.
    // With every support point inside (constrained scheme: s0_i = s(frame i)) the result is
    // the former max(G) - G.front() / G.back() - G.front() exactly.
    double s_lo = 0.0, s_hi = 1.0;
    if (!m_umbrella_s0.empty()) {
        s_lo = std::min(m_umbrella_s0.front(), m_umbrella_s0.back());
        s_hi = std::max(m_umbrella_s0.front(), m_umbrella_s0.back());
    }
    auto interiorRange = [&](const std::vector<double>& x, size_t& k0, size_t& k1) {
        bool found = false;
        for (size_t k = 0; k < x.size(); ++k) {
            if (x[k] < s_lo - 1e-12 || x[k] > s_hi + 1e-12)
                continue;
            if (!found) k0 = k;
            k1 = k;
            found = true;
        }
        return found;
    };
    auto barrier = [&](const std::vector<double>& x, const std::vector<double>& G) {
        size_t k0 = 0, k1 = 0;
        if (G.empty() || !interiorRange(x, k0, k1))
            return 0.0;
        double gmax = G[k0];
        for (size_t k = k0; k <= k1; ++k)
            if (x[k] >= s_lo - 1e-12 && x[k] <= s_hi + 1e-12)
                gmax = std::max(gmax, G[k]);
        return gmax - G[k0];
    };
    auto reaction = [&](const std::vector<double>& x, const std::vector<double>& G) {
        size_t k0 = 0, k1 = 0;
        if (G.empty() || !interiorRange(x, k0, k1))
            return 0.0;
        return G[k1] - G[k0];
    };
    const std::vector<double> GI = integrate(sI, dgI);

    // ---------- estimator (II): mean restraint force per window ----------
    // In a harmonic window the mean force of the restraint equals -dG/ds at <s>; the
    // window width additionally gives the curvature, G'' = kT/sigma^2 - k.
    std::vector<double> sII, dgII, curv;
    std::vector<int> winIdx;
    for (int i = 0; i < m_nimages && !m_constrain; ++i) {
        if (i >= static_cast<int>(m_win_count.size()) || m_win_count[i] < 10)
            continue;
        const double n = static_cast<double>(m_win_count[i]);
        const double sbar = m_win_s_sum[i] / n;
        const double var = std::max(m_win_s2_sum[i] / n - sbar * sbar, 1e-16);
        sII.push_back(sbar);
        dgII.push_back(-windowK(i) * (sbar - m_umbrella_s0[i]));
        curv.push_back(kT / var - windowK(i));
        winIdx.push_back(i);
    }
    // Windows may end up out of order if one slid past its neighbour; sort by <s>.
    std::vector<size_t> order(sII.size());
    for (size_t k = 0; k < order.size(); ++k)
        order[k] = k;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return sII[a] < sII[b]; });
    std::vector<double> sII_s, dgII_s;
    for (size_t k : order) {
        sII_s.push_back(sII[k]);
        dgII_s.push_back(dgII[k]);
    }
    const std::vector<double> GII = integrate(sII_s, dgII_s);

    // ---------- report ----------
    const double bI = barrier(sI, GI) * kcal, bII = barrier(sII_s, GII) * kcal;
    const double rI = reaction(sI, GI) * kcal;
    const double rII = reaction(sII_s, GII) * kcal;
    size_t ki0 = 0, ki1 = 0;
    const int n_inside = interiorRange(sI, ki0, ki1) ? static_cast<int>(ki1 - ki0 + 1) : 0;

    // Claude Generated 2026: grid-free barrier of estimator (I). The maximum over the discrete
    // support points underestimates the barrier whenever the true maximum falls between two
    // of them (coarse grids: 8 windows gave 1.94 vs 2.14 kcal/mol). A natural cubic spline of
    // dG/ds through the INTERIOR support points is integrated exactly and evaluated on a fine
    // grid (fe_profile.h). Same definition as the node barrier: max G - G(first interior
    // support point); reaction dG = G(last) - G(first). The G_kcal column stays trapezoidal.
    curcuma::fe_profile::SplineProfile splI;
    if (n_inside > 0) {
        const std::vector<double> s_in(sI.begin() + ki0, sI.begin() + ki1 + 1);
        const std::vector<double> d_in(dgI.begin() + ki0, dgI.begin() + ki1 + 1);
        splI = curcuma::fe_profile::integrateSpline(s_in, d_in);
        if (!splI.warning.empty())
            CurcumaLogger::warn("NEB-MD mean-force profile (I): " + splI.warning);
    }
    const double bI_spl = splI.barrier_forward * kcal, rI_spl = splI.reaction * kcal;

    // Claude Generated 2026: barrier referenced to the basin free energies (fe_profile.h,
    // basinBarrier) on the spline-integrated G(s): G_A = -kT ln int_{s<=s*} exp(-G/kT) ds,
    // G_B likewise for s >= s*; plus the coverage test on the outermost support points.
    curcuma::fe_profile::BasinBarrier basI;
    curcuma::fe_profile::BasinCoverage covI;
    if (n_inside >= 3) {
        basI = curcuma::fe_profile::basinBarrier(splI.s, splI.G, kT);
        const std::vector<double> s_in(sI.begin() + ki0, sI.begin() + ki1 + 1);
        const std::vector<double> d_in(dgI.begin() + ki0, dgI.begin() + ki1 + 1);
        covI = curcuma::fe_profile::basinCoverage(s_in, d_in);
    } else {
        basI.warning = "fewer than 3 interior support points";
    }

    // Claude Generated 2026: absolute TST rate across s* (fe_profile.h tstRate) - independent of
    // the scale of s, unlike the basin barrier. Only in the constrained ensemble, where the Z^{-1/2}
    // weight makes <Z^{1/2}>_s = count / wsum per bin; in the restrained ensemble the per-bin
    // samples are not drawn from the conditional ensemble at fixed s in the same sense, so no rate.
    curcuma::fe_profile::TstRate tstI;
    curcuma::fe_profile::ZHalfEstimate zhI;
    std::string tst_note;
    if (!m_constrain) {
        tst_note = "not computed - only in the constrained ensemble (-nebmd.constrain true), where <Z^1/2>(s*) = "
                   "count / sum Z^-1/2 per bin";
    } else if (!basI.valid) {
        tst_note = "not computed - no interior maximum s* (" + basI.warning + ")";
    } else {
        std::vector<double> zs, zv;
        std::vector<long> zn;
        for (int b = 0; b < m_fe_bins; ++b) {
            if (m_fe_count[b] < 10 || !(m_fe_wsum[b] > 1e-300))
                continue;
            const double n = static_cast<double>(m_fe_count[b]);
            zs.push_back(m_fe_ssum[b] / n);
            zn.push_back(m_fe_count[b]);
            zv.push_back(n / m_fe_wsum[b]);
        }
        zhI = curcuma::fe_profile::zHalfNear(zs, zn, zv, basI.s_star, 1.5 / m_fe_bins, 50);
        if (!zhI.valid)
            tst_note = "not computed - " + zhI.warning;
        else {
            tstI = curcuma::fe_profile::tstRate(splI.s, splI.G, kT, basI.s_star, zhI.value);
            if (!tstI.valid)
                tst_note = "not computed - " + tstI.warning;
        }
    }

    CurcumaLogger::result("");
    CurcumaLogger::result("=== Free energy from the decoupled scheme (independent windows) ===");
    CurcumaLogger::result_fmt("  samples (pooled)   : {} in {} of {} bins",
                              total, sI.size(), m_fe_bins);
    CurcumaLogger::result_fmt("  interior range     : s = [{:.4f}, {:.4f}] (s of the first / last path frame); "
                              "barrier and reaction dG use only support points inside ({} of {} for (I))",
                              s_lo, s_hi, n_inside, sI.size());
    CurcumaLogger::result_fmt("  (I)  mean force    : barrier (spline) {:.2f} kcal/mol at s = {:.4f}, reaction dG {:.2f}",
                              bI_spl, splI.s_max, rI_spl);
    CurcumaLogger::result_fmt("  (I)  mean force    : barrier (nodes) {:.2f} kcal/mol, reaction dG {:.2f} "
                              "(maximum over the support points)", bI, rI);
    if (basI.valid) {
        CurcumaLogger::result_fmt("  (I)  barrier (basin, fwd/rev): {:.2f} / {:.2f} kcal/mol at s* = {:.4f} "
                                  "(G(s*) - G_A / G(s*) - G_B, G_X = -kT ln int_X exp(-G/kT) ds)",
                                  basI.barrier_fwd * kcal, basI.barrier_rev * kcal, basI.s_star);
        CurcumaLogger::result_fmt("  (I)  reaction dG (basin): {:.2f} kcal/mol (G_B - G_A; scale depends on the CV "
                                  "parametrisation, compare only on the same CV)", basI.reaction * kcal);
    } else {
        CurcumaLogger::warn("NEB-MD mean-force profile (I): barrier (basin, fwd/rev) not defined - " + basI.warning);
    }
    if (n_inside >= 3 && (!covI.left || !covI.right))
        CurcumaLogger::warn(fmt::format(
            "basin not covered ({}{}{}): dG/ds at the outermost support point is still > 0.5 max|dG/ds| "
            "(left {:.3e}, right {:.3e}, max {:.3e} Eh), so the windows do not reach the basin bottom and "
            "the basin barrier is not referenced to the basin. Extend the CV (cv_extend, or a larger value) and "
            "check that the path ends in minima.",
            covI.left ? "" : "left", (!covI.left && !covI.right) ? ", " : "", covI.right ? "" : "right",
            covI.slope_left, covI.slope_right, covI.max_abs));
    // Claude Generated 2026: TST rate (absolute, independent of the scale of s)
    if (tstI.valid) {
        CurcumaLogger::result_fmt("  (I)  TST rate (fwd/rev): {:.4e} / {:.4e} 1/ps at s* = {:.4f}",
                                  tstI.k_fwd, tstI.k_rev, tstI.s_star);
        CurcumaLogger::result_fmt("  (I)  TST dG (Eyring, fwd/rev): {:.2f} / {:.2f} kcal/mol (kT ln((kB T/h)/k), "
                                  "kB T/h = {:.4f} 1/ps)", tstI.dG_fwd * kcal, tstI.dG_rev * kcal,
                                  tstI.eyring_prefactor);
        CurcumaLogger::result_fmt("  (I)  TST detail    : <|sdot|>(s*) {:.4e} s/ps, p(s*) fwd/rev {:.4e} / {:.4e} 1/s, "
                                  "<Z^1/2> {:.4e} from {} samples in {} bin(s){}",
                                  tstI.sdot_abs, tstI.p_fwd, tstI.p_rev, zhI.value, zhI.samples, zhI.bins,
                                  zhI.fallback ? " (fallback)" : "");
        if (zhI.fallback)
            CurcumaLogger::warn("NEB-MD TST rate: " + zhI.warning);
    } else {
        CurcumaLogger::result("  (I)  TST rate      : " + tst_note);
    }
    if (m_constrain)
        CurcumaLogger::result("  (II) restraint     : not defined - there is no restraint to read a "
                              "force from when s is held by a constraint");
    else
        CurcumaLogger::result_fmt("  (II) restraint     : barrier {:.2f} kcal/mol, reaction dG {:.2f}",
                                  bII, rII);
    const double spread = m_constrain ? 0.0 : std::abs(bI - bII);
    CurcumaLogger::result_fmt("  agreement          : {:.2f} kcal/mol between two independent "
                              "estimators of the same profile (node barriers of (I) and (II))", spread);
    if (spread > 0.1 * std::max(std::abs(bI), 1.0))
        CurcumaLogger::warn("The two estimators disagree by more than 10 %. They differ only in how "
                            "the metric contribution enters (explicitly via div b, implicitly via "
                            "the sampled density), so the discrepancy points at that term, at the "
                            "div-b probe settings, or at insufficient sampling in some bins.");

    // Metric contribution on its own - it is now a measured quantity, not an estimate.
    double dmin = 0.0, dmax = 0.0;
    bool first = true;
    for (int b = 0; b < m_fe_bins; ++b) {
        if (m_fe_count[b] < 10)
            continue;
        const double wn = (m_fe_wsum[b] > 1e-300) ? m_fe_wsum[b] : static_cast<double>(m_fe_count[b]);
        const double v = -kT * m_fe_div[b] / wn;
        if (first) { dmin = dmax = v; first = false; }
        dmin = std::min(dmin, v);
        dmax = std::max(dmax, v);
    }
    CurcumaLogger::result_fmt("  metric -kT<div b>  : {:.3e} .. {:.3e} Eh per unit s "
                              "({} probe(s), eps = {:.1e} A)", dmin, dmax, m_div_probes, m_div_epsilon);

    std::ofstream f(outputPath(Basename() + ".neb.fe.csv"));
    if (f) {
        f << "# Free energy from the decoupled scheme (Claude Generated 2026)\n";
        f << "# estimator I: binned conditional mean force  dG/ds = <b.grad E> - kT <div b>\n";
        f << "# a restraint depending only on s does not bias conditional averages, so all\n";
        f << "# windows are pooled into the same bins; no unbiasing and no overlap needed\n";
        // Claude Generated 2026: interior range + per-point sample count and Z-floor
        // rejection fraction rejected/(samples+rejected) (bin of a rejected sample from its s).
        f << fmt::format("# interior s range [{:.6f}, {:.6f}] = s(first frame), s(last frame); barrier and "
                         "reaction dG use only support points inside\n", s_lo, s_hi);
        f << fmt::format("# estimator I barrier (spline): {:.4f} kcal/mol at s = {:.6f}, reaction dG {:.4f} kcal/mol "
                         "(natural cubic spline of dG/ds through the interior support points, integrated exactly; "
                         "barrier = max G - G(first interior support point){})\n",
                         bI_spl, splI.s_max, rI_spl, splI.spline ? "" : "; FALLBACK trapezoid: " + splI.warning);
        f << fmt::format("# estimator I barrier (nodes): {:.4f} kcal/mol, reaction dG {:.4f} kcal/mol "
                         "(trapezoid, maximum over the support points; G_kcal column below)\n", bI, rI);
        // Claude Generated 2026: basin-referenced barrier (fe_profile.h basinBarrier) + CV info.
        if (basI.valid) {
            f << fmt::format("# estimator I barrier (basin, fwd/rev): {:.4f} / {:.4f} kcal/mol at s* = {:.6f} "
                             "(G_X = -kT ln int_X exp(-G/kT) ds over the spline, X = s <= s* / s >= s*; ds in units "
                             "of the CV coordinate s, so the absolute value depends on the CV parametrisation)\n",
                             basI.barrier_fwd * kcal, basI.barrier_rev * kcal, basI.s_star);
            f << fmt::format("# estimator I reaction dG (basin): {:.4f} kcal/mol\n", basI.reaction * kcal);
        } else {
            f << "# estimator I barrier (basin, fwd/rev): not defined - " << basI.warning << "\n";
        }
        // Claude Generated 2026: absolute TST rate + Eyring free energy (fe_profile.h tstRate)
        if (tstI.valid) {
            f << fmt::format("# TST rate (fwd/rev): {:.6e} / {:.6e} 1/ps at s* = {:.6f} (k = 1/2 <|sdot|>(s*) p_X(s*), "
                             "p_X = exp(-G(s*)/kT) / int_X exp(-G/kT) ds over the spline; independent of the scale of s)\n",
                             tstI.k_fwd, tstI.k_rev, tstI.s_star);
            f << fmt::format("# TST dG (Eyring, fwd/rev): {:.4f} / {:.4f} kcal/mol (kT ln((kB T/h)/k), kB T/h = {:.6f} 1/ps)\n",
                             tstI.dG_fwd * kcal, tstI.dG_rev * kcal, tstI.eyring_prefactor);
            f << fmt::format("# TST detail: <|sdot|>(s*) {:.6e} s/ps = sqrt(2kT/pi) <Z^1/2>, <Z^1/2> {:.6e} amu^-1/2 A^-1 "
                             "(count / sum Z^-1/2) from {} samples in {} bin(s){}; p(s*) fwd/rev {:.6e} / {:.6e} 1/s\n",
                             tstI.sdot_abs, zhI.value, zhI.samples, zhI.bins,
                             zhI.fallback ? " (FALLBACK: " + zhI.warning + ")" : "", tstI.p_fwd, tstI.p_rev);
        } else {
            f << "# TST rate: " << tst_note << "\n";
        }
        if (n_inside >= 3 && (!covI.left || !covI.right))
            f << fmt::format("# basin not covered: left {} right {} (dG/ds first {:.4e}, last {:.4e}, max|dG/ds| "
                             "{:.4e} Eh)\n", covI.left ? "ok" : "NOT", covI.right ? "ok" : "NOT",
                             covI.slope_left, covI.slope_right, covI.max_abs);
        if (!m_cv_path.empty() && m_cv_extend_length > 0.0) // Claude Generated 2026 (only in length mode)
            f << fmt::format("# cv_extend_length {:.6g} of the path length = {} frames per end\n",
                             m_cv_extend_length, m_cv_extend);
        if (!m_cv_atoms.empty()) // Claude Generated 2026 (only with a selection: default header unchanged)
            f << fmt::format("# cv_atoms {} ({} atoms)\n", m_cv_atoms_spec, m_cv_atoms.size());
        if (!m_cv_path.empty() && m_cv_orig_s.size() == 4)
            f << fmt::format("# cv path: {} frames (cv_frames {}, cv_extend {}), lambda {:.6g}; original path frames "
                             "on this CV: s(first) {:.6f}, s(second) {:.6f}, s(second-to-last) {:.6f}, s(last) {:.6f}\n",
                             m_cv_path.size(), m_cv_frames, m_cv_extend, m_path_cv ? m_path_cv->lambda() : 0.0,
                             m_cv_orig_s[0], m_cv_orig_s[1], m_cv_orig_s[2], m_cv_orig_s[3]);
        if (m_z_wall > 0.0) // Claude Generated 2026: only with a wall, so default output is unchanged
            f << fmt::format("# z wall {:.6g} {}, k {:.6g} {}: PMF of the z-restricted system, compare only "
                             "at equal z_wall\n", m_z_wall, m_path_cv ? m_path_cv->zUnit() : std::string("A^2"), m_z_wall_k,
                             (m_path_cv && m_path_cv->zUnit() == "A^-2") ? "Eh A^4" : "Eh/A^4");
        if (m_cv_metric == PathCV::Metric::Distances) // Claude Generated 2026 (only for distances: default header unchanged)
            f << "# cv_metric distances, " << cvPairsText() << "\n";
        // Claude Generated 2026: the last three columns split estimator I, dGds_Eh = gradE + metric + wall
        // (printed with 13 significant digits; dGds_Eh keeps its 8 decimals). Not defined for II: nan.
        f << "estimator,s,dGds_Eh,G_kcal,samples,rejected_frac,dGds_gradE_Eh,dGds_metric_Eh,dGds_wall_Eh\n";
        auto frac = [](long n, long r) { return (n + r) > 0 ? static_cast<double>(r) / static_cast<double>(n + r) : 0.0; };
        for (size_t k = 0; k < sI.size(); ++k)
            f << "I," << fmt::format("{:.6f}", sI[k]) << "," << fmt::format("{:.8f}", dgI[k]) << ","
              << fmt::format("{:.4f}", GI[k] * kcal) << "," << nI[k] << ","
              << fmt::format("{:.6f}", frac(nI[k], rejI[k])) << ","
              << fmt::format("{:.12e},{:.12e},{:.12e}", splitI[k].grad_e, splitI[k].metric, splitI[k].wall) << "\n";
        for (size_t k = 0; k < sII_s.size(); ++k) {
            const int i = winIdx[order[k]];
            const long r = i < static_cast<int>(m_win_rejected.size()) ? m_win_rejected[i] : 0;
            f << "II," << fmt::format("{:.6f}", sII_s[k]) << "," << fmt::format("{:.8f}", dgII_s[k])
              << "," << fmt::format("{:.4f}", GII[k] * kcal) << "," << m_win_count[i] << ","
              << fmt::format("{:.6f}", frac(m_win_count[i], r)) << ",nan,nan,nan\n";
        }
    }

    // Conditioning: can the coordinate resolve the images at all?
    if (m_fe_nonfinite > 0)
        CurcumaLogger::warn(fmt::format(
            "{} samples were dropped because the window's MD produced non-finite values. Those "
            "windows are dead - the profile below rests on the survivors only and must not be "
            "quoted without checking which windows contributed.", m_fe_nonfinite));
    if (m_fe_rejected > 0 || m_bias_capped > 0) {
        const double rej = 100.0 * m_fe_rejected / std::max<long>(total + m_fe_rejected, 1);
        CurcumaLogger::result_fmt("  discarded samples  : {:.1f} % (|grad s| collapsed); "
                                  "restraint force capped {} times", rej, m_bias_capped);
    }
    // Claude Generated 2026: per-window distance from the path (z) and the z wall statistics.
    for (int i = 0; i < m_nimages; ++i) {
        const double n = static_cast<double>(m_win_count[i]);
        if (n < 1.0)
            continue;
        CurcumaLogger::result_fmt("  window {:2d}: s0 {:.4f}  <z> {:.4f}  max z {:.4f} {}  z > z_wall {:.4f}  <W> {:.3e} Eh",
                                  i, m_umbrella_s0[i], m_win_z_sum[i] / n, m_win_z_max[i],
                                  m_path_cv ? m_path_cv->zUnit() : std::string("A^2"), m_win_zwall_n[i] / n, m_win_wall_sum[i] / n);
    }
    double amp = 0.0;
    int ampn = 0;
    for (size_t k = 0; k < winIdx.size(); ++k) {
        const int i = winIdx[k];
        const double n = static_cast<double>(m_win_count[i]);
        if (n < 1.0)
            continue;
        amp += std::sqrt(std::max(m_win_z_sum[i] / n, 0.0)); // z is an RMSD^2
        ++ampn;
    }
    if (ampn > 0 && m_path_spacing > 0.0) {
        amp /= ampn;
        const double ratio = m_path_spacing / std::max(amp, 1e-12);
        CurcumaLogger::result_fmt("  conditioning       : path spacing {:.4f} A vs thermal "
                                  "amplitude {:.4f} A -> ratio {:.2f}", m_path_spacing, amp, ratio);
        if (ratio < 2.0)
            CurcumaLogger::warn("The reference images lie closer together than the molecule moves on "
                                "its own. The path coordinate then cannot resolve them, s saturates "
                                "at 0 or 1, and no estimator can repair that. Use FEWER, more widely "
                                "spaced images - the mean-force route needs no window overlap, so "
                                "wide spacing costs nothing here (unlike WHAM, which demands the "
                                "opposite and is therefore not applicable to such a short path).");
    }

    // Window diagnostics: where the windows actually sat, and how wide they were.
    CurcumaLogger::result("  windows (target -> sampled, width, curvature):");
    for (size_t k = 0; k < winIdx.size(); ++k) {
        const int i = winIdx[k];
        const double n = static_cast<double>(m_win_count[i]);
        const double sbar = m_win_s_sum[i] / n;
        const double sd = std::sqrt(std::max(m_win_s2_sum[i] / n - sbar * sbar, 0.0));
        CurcumaLogger::result(fmt::format(
            "    {:>3}  s0 = {:.4f} -> <s> = {:.4f}   sigma = {:.4f}   G'' = {:>10.2f} kcal/mol",
            i, m_umbrella_s0[i], sbar, sd, curv[k] * kcal));
    }
    for (int i = 0; i < m_nimages; ++i) {
        m_bias_capped += m_img_capped[i];
        if (m_img_max_jump[i] > m_max_s_jump) {
            m_max_s_jump = m_img_max_jump[i];
            m_max_jump_step = static_cast<int>(m_img_max_jump_step[i]);
            m_max_jump_image = i;
        }
    }
    if (m_max_jump_step >= 0)
        CurcumaLogger::result_fmt("  largest jump of s  : {:.4f} in one step (image {}, step {}) "
                                  "- a large value indicates a discontinuous CV gradient",
                                  m_max_s_jump, m_max_jump_image, m_max_jump_step);
    // Holding criterion (Claude Generated 2026). Without it the other diagnostics are
    // degenerate: every one of them - estimator agreement, discard rate, conditioning -
    // improves monotonically as the restraint is weakened, and is perfect in the limit
    // k -> 0, where the windows have simply abandoned the barrier and the profile is flat.
    // A window that is not where it was put measures the basin it slid into, not the
    // barrier. Measured on helicene: at k = 0.5 all six windows had collapsed onto the two
    // path ends while every other indicator looked best of the whole screening.
    if (!winIdx.empty()) {
        double s0min = 1e30, s0max = -1e30;
        for (int i = 0; i < m_nimages; ++i) {
            s0min = std::min(s0min, m_umbrella_s0[i]);
            s0max = std::max(s0max, m_umbrella_s0[i]);
        }
        const double dsw = (s0max - s0min) / std::max(m_nimages - 1, 1);
        double worst_dev = 0.0;
        int worst_win = -1;
        for (size_t k = 0; k < winIdx.size(); ++k) {
            const int i = winIdx[k];
            const double sbar = m_win_s_sum[i] / static_cast<double>(m_win_count[i]);
            const double dev = std::abs(sbar - m_umbrella_s0[i]) / std::max(dsw, 1e-12);
            if (dev > worst_dev) { worst_dev = dev; worst_win = i; }
        }
        CurcumaLogger::result_fmt("  holding            : worst |<s>-s0| = {:.2f} window spacings "
                                  "(window {}); < 0.5 means the windows stay where they were put",
                                  worst_dev, worst_win);
        if (worst_dev > 0.5) {
            CurcumaLogger::warn(fmt::format(
                "The windows have left their setpoints by up to {:.1f} spacings. The restraint is "
                "too weak against the slope of the profile, so the barrier region is only crossed, "
                "not sampled, and the reported barrier is a lower bound at best. Holding requires "
                "roughly k > (dG/ds)/(ds/2); note that FEWER windows enlarge ds and therefore lower "
                "the stiffness needed - which also suits the path CV, whose reference frames must "
                "be spaced well above the thermal amplitude.", worst_dev));
        }
    }

    CurcumaLogger::warn("Decoupled scheme: AI-implemented and machine-tested only. Validate on a "
                        "system with a known barrier before quoting numbers.");
}




void NebMD::writeBandFile(const std::string& filename, int step, const std::string& tag) const
{
    std::string path = outputPath(filename);
    std::vector<Geometry> R = currentPath();
    std::vector<double> E = currentEnergies();
    for (int i = 0; i < m_nimages; ++i) {
        Molecule mol = m_start; // copies atom labels
        mol.setGeometry(R[i]);
        double T = m_images[i] ? m_images[i]->currentTemperature() : 0.0;
        std::string comment = tag + " step=" + std::to_string(step)
                            + " image=" + std::to_string(i)
                            + " E=" + fmt::format("{:.6f}", E[i])
                            + " T=" + fmt::format("{:.1f}", T);
        mol.appendXYZFile(path, comment);
    }
}

void NebMD::writeInitialPath() const
{
    writeBandFile(Basename() + ".neb.init.xyz", 0, "initial");
}

void NebMD::writePathSnapshot(int step) const
{
    writeBandFile(Basename() + ".neb.path.xyz", step, "path");
}

void NebMD::writeFinalPath() const
{
    writeBandFile(Basename() + ".neb.final.xyz", m_step, "final");
}

void NebMD::writeEnergyTable(int step) const
{
    std::string path = outputPath(Basename() + ".neb.energy.csv");
    std::vector<double> E = currentEnergies();
    std::ofstream f(path, std::ios::app);
    if (!f)
        return;
    for (int i = 0; i < m_nimages; ++i) {
        double ekin = m_images[i] ? m_images[i]->kineticEnergy() : 0.0;
        double T = m_images[i] ? m_images[i]->currentTemperature() : 0.0;
        f << step << "," << i << ","
          << fmt::format("{:.6f}", E[i]) << ","
          << fmt::format("{:.6f}", ekin) << ","
          << fmt::format("{:.2f}", T) << "\n";
    }
}

void NebMD::start()
{
    if (!Initialise())
        return;
    // Register band outputs for copy-back to CWD (BMT mode). Not in the analysis-only
    // cv_histogram mode, which writes none of them (Claude Generated 2026).
    if (m_cv_histogram.empty()) {
        addBakFile(Basename() + ".neb.path.xyz");
        addBakFile(Basename() + ".neb.energy.csv");
        addBakFile(Basename() + ".neb.final.xyz");
        if (m_write_initial_path)
            addBakFile(Basename() + ".neb.init.xyz");
    }
    if (!buildInitialPath())
        return;
    // Claude Generated 2026: analysis-only mode. The reference path is used exactly as read /
    // interpolated - IDPP, overlap repair and band optimisation would change it, and a
    // histogram of s is only comparable to stage 2 on the identical path.
    if (!m_cv_histogram.empty()) {
        if (!buildCvPath(false)) // same CV construction as stage 2 (Claude Generated 2026)
            return;
        histogramTrajectory(m_cv_histogram);
        addBakFile(Basename() + ".neb.ref.csv");
        return;
    }
    CitationRegistry::cite("neb"); // improved tangent (Henkelman & Jonsson 2000)
    if (m_idpp)
        buildIDPPPath();       // repulsion-aware, neighbour-coupled path refinement
    if (m_repair_overlaps)
        repairPathOverlaps();  // clean up any residual close contacts
    checkPathTopology();       // warn if the connectivity still differs from image 0
    if (m_optimize)
        optimizeBand();        // relax the band onto the MEP (classical CI-NEB) before the MD
    if (!buildCvPath(true))    // cv_frames / cv_extend: CV frames + windows (Claude Generated 2026)
        return;
    if (m_string && m_string_prepare)
        prepareStringStartPath(); // align + reparametrise BEFORE the images are built (Claude Generated 2026)
    if (!prepareImages())
        return;
    if (m_string) {
        setupString();     // freeze the starting mean path M, choose k (Claude Generated 2026)
        addBakFile(Basename() + ".neb.string.xyz");
        addBakFile(Basename() + ".neb.string_history.csv");
        if (m_string_write_updates) {
            std::remove(outputPath(Basename() + ".neb.string_updates.xyz").c_str());
            addBakFile(Basename() + ".neb.string_updates.xyz");
        }
    }
    if (m_umbrella)
        setupUmbrella();   // freeze the path, build the CV, place the windows
    if (m_write_initial_path)
        writeInitialPath();

    // CSV header (overwrite any stale file).
    {
        std::string path = outputPath(Basename() + ".neb.energy.csv");
        std::ofstream f(path);
        if (f) f << "step,image,Epot,Ekin,T\n";
    }
    writePathSnapshot(0);
    writeEnergyTable(0);

    CurcumaLogger::result_fmt("NEB-MD: {} images (nudged={}, fixed_endpoints={}, k={}, mass-weighted={})",
                              m_nimages, m_nudged, m_fixed_endpoints, m_k_spring, m_mass_weighted_spring);

    // Safety cap against a runaway loop; step() normally terminates at max_time.
    const long hard_cap = 100000000L;
    long guard = 0;
    while (stepBand()) {
        if (CheckStop())
            break;
        if (++guard > hard_cap) {
            CurcumaLogger::warn("NEB-MD: hard step cap reached, stopping band integration.");
            break;
        }
    }

    writeFinalPath();
    if (m_string && m_string_aborted) {
        // Claude Generated 2026: no string path from an aborted stage - remove a stale one (e.g.
        // from an earlier -no_bmt run in the same directory) so no consumer can pick it up.
        const std::string spath = outputPath(Basename() + ".neb.string.xyz");
        std::remove(spath.c_str());
        CurcumaLogger::error(fmt::format("NEB-MD string: stage aborted at step {} (dead image, string_abort_on_dead); "
                                         "no string path written. History: {}",
                                         m_step, outputPath(Basename() + ".neb.string_history.csv")));
        return;
    }
    if (m_string) {
        // Claude Generated 2026: stage-1 result = the mean path, consecutively aligned.
        writeStringPath(Basename() + ".neb.string.xyz");
        // Claude Generated 2026: with string_abort_on_dead = false the stage runs on with frozen
        // images; say loudly which frames of the written path were never relaxed.
        std::string frozen;
        int ndead = 0;
        for (int i = 0; i < m_nimages && i < static_cast<int>(m_window_died.size()); ++i)
            if (m_window_died[i] >= 0) {
                frozen += fmt::format("{}{} (dead at step {})", frozen.empty() ? "" : ", ", i, m_window_died[i]);
                ++ndead;
            }
        if (ndead > 0)
            CurcumaLogger::error(fmt::format("NEB-MD string: string stage finished with {} dead images: path frames {} "
                                             "are NOT relaxed (frozen since their death; the convergence "
                                             "criterion was never met with dead images).", ndead, frozen));
        long capped = 0;
        for (long c : m_string_capped)
            capped += c;
        // Claude Generated 2026: the end message names the criterion and its numbers.
        const bool se_ok = std::isfinite(m_string_last_se) && m_string_last_se < m_string_ftol;
        const std::string se_txt = std::isfinite(m_string_last_se) ? fmt::format("{:.2e}", m_string_last_se) : std::string("n/a");
        if (m_string_converged && m_string_ftol > 0.0)
            CurcumaLogger::result_fmt("NEB-MD string: converged (force criterion) after {} updates (step {}): "
                                      "F_sig_max = {:.2e} < ftol {:.2e} Eh/A, SE_max = {}, last delta_max = {:.4f} A",
                                      m_string_updates, m_step, m_string_last_fsig, m_string_ftol, se_txt,
                                      m_string_last_delta);
        else if (m_string_converged)
            CurcumaLogger::result_fmt("NEB-MD string: converged (former criterion delta_max < string_tol, "
                                      "string_ftol = 0) after {} updates (step {}), last delta_max = {:.4f} A",
                                      m_string_updates, m_step, m_string_last_delta);
        else if (m_string_ftol > 0.0 && !se_ok) {
            CurcumaLogger::warn(fmt::format("NEB-MD string: string stage ended at string_max_time: force criterion "
                                            "noise-limited (SE_max {} >= ftol {:.2e} Eh/A) after {} updates (step {}); "
                                            "resolving it needs string_update x {} or a longer average. Not an error; "
                                            "the written path is the last mean path.",
                                            se_txt, m_string_ftol, m_string_updates, m_step,
                                            std::isfinite(m_string_last_se) ? fmt::format("{:.1f}", (m_string_last_se / m_string_ftol) * (m_string_last_se / m_string_ftol)) : std::string("n/a")));
            if (m_string_last_z > 3.0 && m_string_last_fsig > m_string_ftol)
                CurcumaLogger::warn(fmt::format("NEB-MD string: path not relaxed - the perpendicular mean force is "
                                                "significant (z_max = {:.1f}, F_sig_max = {:.2e} > ftol {:.2e} Eh/A).",
                                                m_string_last_z, m_string_last_fsig, m_string_ftol));
        } else
            CurcumaLogger::warn(fmt::format("NEB-MD string: NOT converged after {} updates (step {}): F_sig_max = {:.2e} "
                                            "Eh/A (ftol {:.2e}), SE_max = {}, z_max = {:.1f}, last delta_max = {:.4f} A "
                                            "(tol {:.4f} A). The written path is the last mean path.",
                                            m_string_updates, m_step, m_string_last_fsig, m_string_ftol, se_txt,
                                            m_string_last_z, m_string_last_delta, m_string_tol));
        CurcumaLogger::result_fmt("NEB-MD string: restraint force capped {} times", capped);
        CurcumaLogger::result_fmt("NEB-MD string: mean path written to {} ({} frames)",
                                  outputPath(Basename() + ".neb.string.xyz"), m_nimages);
    }
    if (m_umbrella)
        writeWHAM();
    if (m_independent) {
        writeMeanForceProfile();
        writeWindowDiagnostics(); // Claude Generated 2026
    }

    // Final summary: per-image energy.
    {
        std::vector<double> E = currentEnergies();
        std::string line = "NEB-MD final Epot [Eh]:";
        for (int i = 0; i < m_nimages; ++i)
            line += " " + fmt::format("{:.6f}", E[i]);
        CurcumaLogger::result(line);
    }
    CurcumaLogger::result("NEB-MD: band integration finished after "
                        + std::to_string(m_step) + " steps. Outputs in "
                        + (OutputDir().empty() ? std::string("CWD") : OutputDir())
                        + " (*.neb.path.xyz, *.neb.energy.csv, *.neb.final.xyz"
                        + (m_write_initial_path ? ", *.neb.init.xyz" : "") + ").");
}

nlohmann::json NebMD::WriteRestartInformation()
{
    // Minimal band state (positions of all images); velocities are owned by the
    // per-image SimpleMD instances and are not serialised here in v1.
    json state;
    state["nimages"] = m_nimages;
    state["nudged"] = m_nudged;
    state["fixed_endpoints"] = m_fixed_endpoints;
    state["k_spring"] = m_k_spring;
    state["mass_weighted_spring"] = m_mass_weighted_spring;
    state["step"] = m_step;
    std::vector<Geometry> R = currentPath();
    json geoms = json::array();
    for (int i = 0; i < m_nimages; ++i) {
        json g = json::array();
        for (int a = 0; a < R[i].rows(); ++a)
            g.push_back({R[i](a, 0), R[i](a, 1), R[i](a, 2)});
        geoms.push_back(g);
    }
    state["geometries"] = geoms;
    return state;
}

bool NebMD::LoadRestartInformation()
{
    // v1 does not implement NEB-MD restart; the per-image MD restart machinery
    // is owned by SimpleMD. Returning true keeps the CurcumaMethod contract happy.
    return true;
}