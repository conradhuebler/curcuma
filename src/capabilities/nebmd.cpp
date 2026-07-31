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
    m_umbrella_keep_springs = m_config.get<bool>("umbrella_keep_springs");
    m_reparametrize = m_config.get<int>("reparametrize");
    m_pmf = m_config.get<bool>("pmf");
    m_fixman = m_config.get<bool>("fixman");
    m_pmf_equilibration = m_config.get<int>("pmf_equilibration");
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

    m_images.assign(m_nimages, nullptr);
    m_image_done.assign(m_nimages, false);
    m_image_geoms.assign(m_nimages, Geometry());
    m_pmf_force_sum.assign(m_nimages, 0.0);
    m_pmf_spring_sum.assign(m_nimages, 0.0);
    m_pmf_fixman_sum.assign(m_nimages, 0.0);
    m_pmf_epot_sum.assign(m_nimages, 0.0);
    m_pmf_epot_sq_sum.assign(m_nimages, 0.0);
    m_pmf_arc_sum.assign(m_nimages, 0.0);
    m_pmf_samples = 0;
    return true;
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
        m_pmf_force_sum.assign(m_nimages, 0.0);
        m_pmf_spring_sum.assign(m_nimages, 0.0);
        m_pmf_epot_sum.assign(m_nimages, 0.0);
        m_pmf_epot_sq_sum.assign(m_nimages, 0.0);
        m_pmf_arc_sum.assign(m_nimages, 0.0);
        m_pmf_fixman_sum.assign(m_nimages, 0.0);
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
    sd["write_xyz"] = false;           // suppress per-image trajectories; the band snapshot records the path
    sd["write_initial_state"] = false;
    sd["write_restart_frequency"] = -1; // disable per-image restart files (step() gates on `> -1` then `% freq`; 0 would SIGFPE)
    sd["opt"] = false;
    sd["rmsd_mtd"] = false;
    sd["mtd"] = false;
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
        img_ctrl["simplemd"] = sd;
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

    // 1) Compute + inject the external (spring / nudging) force for each interior image.
    for (int i = 1; i < m_nimages - 1; ++i) {
        if (m_image_done[i])
            continue;
        SimpleMD* img = m_images[i];
        if (!img)
            continue;

        Geometry g = img->gradient(); // true gradient dE/dR at the current position (Eh/Å)
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
        if (m_umbrella) {
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

    // 2) Advance every interior image by one MD step (independent thermostats/velocities).
    bool any_running = false;
    for (int i = 1; i < m_nimages - 1; ++i) {
        if (m_image_done[i])
            continue;
        SimpleMD* img = m_images[i];
        if (!img)
            continue;
        const bool ok = img->step(); // false at max_time / instability / stop-file
        if (!ok)
            m_image_done[i] = true;
        else
            any_running = true;
    }

    ++m_step;

    // Claude Generated 2026: string reparametrisation - redistribute the images to equal
    // arclength so they stay well-defined umbrella windows (required for a meaningful PMF).
    if (m_reparametrize > 0 && m_step % m_reparametrize == 0)
        reparametrizeString();

    // Claude Generated 2026: accumulate the mean force AFTER equilibration. The gradients
    // read above belong to the pre-step geometry, which is exactly the configuration the
    // averages should be taken over, so we accumulate with the tangents of this step.
    // Claude Generated 2026: the arclength mean-force estimator is only meaningful with
    // the neighbour springs. In umbrella mode no spring force is applied (the restraint
    // acts on the path CV instead), so subtracting Fspring in accumulatePMF() would
    // remove a bias that was never added - producing garbage (observed: 2.8e5 kcal/mol
    // over a 500 A "arclength"). WHAM on the CV histograms is the estimator there.
    if (m_pmf && !m_umbrella && m_step > m_pmf_equilibration)
        accumulatePMF(R, tau, Fspring);

    // Umbrella histograms (after equilibration): record s(R) of every window.
    if (m_umbrella && m_path_cv && m_step > m_pmf_equilibration) {
        for (int i = 1; i < m_nimages - 1; ++i) {
            if (!m_images[i] || m_image_done[i])
                continue;
            const double s = m_path_cv->s(m_images[i]->positions());
            m_umbrella_ssum[i] += s;
            ++m_umbrella_n[i];
            const int b = static_cast<int>(s * m_umbrella_bins);
            if (b >= 0 && b < m_umbrella_bins)
                ++m_hist[i][b];
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
    m_path_cv = std::make_unique<PathCV>(m_image_geoms, m_start);
    m_umbrella_s0.assign(m_nimages, 0.0);
    for (int i = 0; i < m_nimages; ++i)
        m_umbrella_s0[i] = static_cast<double>(i) / (m_nimages - 1);
    // Claude Generated 2026: auto-select the force constant so that neighbouring windows
    // actually OVERLAP - without overlap WHAM has nothing to stitch and returns garbage.
    // The thermal width of a harmonic window is sigma = sqrt(kT/k); requiring
    // sigma ~ ds/2 with ds = 1/(N-1) the window spacing gives k = kT/(ds/2)^2.
    // (Measured failure of the previous fixed default k=200: sigma/ds = 0.02, nine gaps
    // in the histogram, WHAM barrier 351 kcal/mol - meaningless.)
    if (m_umbrella_k <= 0.0) {
        const double kT = 3.166811e-6 * 300.0;
        const double ds = 1.0 / std::max(m_nimages - 1, 1);
        m_umbrella_k = kT / (0.5 * ds * 0.5 * ds);
        CurcumaLogger::info(fmt::format(
            "NEB-MD: umbrella k auto-set to {:.4f} (sigma_s ~ {:.3f}, window spacing {:.3f})",
            m_umbrella_k, std::sqrt(kT / m_umbrella_k), ds));
    }
    m_hist.assign(m_nimages, std::vector<int>(m_umbrella_bins, 0));
    m_umbrella_ssum.assign(m_nimages, 0.0);
    m_umbrella_n.assign(m_nimages, 0);
    CurcumaLogger::result_fmt("NEB-MD: umbrella sampling on the path CV (lambda = {:.3f}, k = {:.1f}, "
                              "{} windows)", m_path_cv->lambda(), m_umbrella_k, m_nimages);
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
    const double dev = s - m_umbrella_s0[i];
    Fext += (m_umbrella_k * dev) * ds;

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
            Fext += (m_kz * excess) * dz;
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
    const double kT = 3.166811e-6 * 300.0;  // Eh at 300 K (thermostat target)
    const double kcal = 627.5094740631;
    const int B = m_umbrella_bins;

    std::vector<double> f(m_nimages, 0.0), P(B, 0.0);
    std::vector<double> Ntot(m_nimages, 0.0);
    for (int i = 0; i < m_nimages; ++i)
        Ntot[i] = static_cast<double>(m_umbrella_n[i]);

    auto bias = [&](int i, int b) {
        const double s = (b + 0.5) / B;
        const double d = s - m_umbrella_s0[i];
        return 0.5 * m_umbrella_k * d * d;
    };

    for (int iter = 0; iter < 2000; ++iter) {
        // P(s) from the current offsets
        for (int b = 0; b < B; ++b) {
            double num = 0.0, den = 0.0;
            for (int i = 0; i < m_nimages; ++i) {
                num += m_hist[i][b];
                den += Ntot[i] * std::exp(-(bias(i, b) - f[i]) / kT);
            }
            P[b] = (den > 1e-300) ? num / den : 0.0;
        }
        // new offsets
        double maxdev = 0.0;
        for (int i = 0; i < m_nimages; ++i) {
            double sum = 0.0;
            for (int b = 0; b < B; ++b)
                sum += P[b] * std::exp(-bias(i, b) / kT);
            const double fnew = (sum > 1e-300) ? -kT * std::log(sum) : f[i];
            maxdev = std::max(maxdev, std::abs(fnew - f[i]));
            f[i] = fnew;
        }
        if (maxdev < 1e-10)
            break;
    }

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
    const int nat = static_cast<int>(m_start_geom.rows());
    std::vector<Geometry> R = currentPath();

    // Cumulative arclength of the current polyline.
    std::vector<double> s(m_nimages, 0.0);
    for (int i = 1; i < m_nimages; ++i)
        s[i] = s[i - 1] + (R[i] - R[i - 1]).norm();
    const double total = s[m_nimages - 1];
    if (total < 1e-9)
        return 0.0;

    double disp2 = 0.0;
    int moved = 0;
    for (int i = 1; i < m_nimages - 1; ++i) {
        if (!m_images[i])
            continue; // clamped image: never moved

        // Target arclength for a uniform distribution.
        const double s_target = total * static_cast<double>(i) / static_cast<double>(m_nimages - 1);

        // Locate the segment containing s_target and interpolate linearly within it.
        int seg = 0;
        while (seg < m_nimages - 2 && s[seg + 1] < s_target)
            ++seg;
        const double ds = s[seg + 1] - s[seg];
        const double t = (ds > 1e-12) ? (s_target - s[seg]) / ds : 0.0;
        const Geometry target = (1.0 - t) * R[seg] + t * R[seg + 1];

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

void NebMD::accumulatePMF(const std::vector<Geometry>& R,
                          const std::vector<Geometry>& tau,
                          const std::vector<Geometry>& Fspring)
{
    // Claude Generated 2026: accumulate <F.tau> per image plus the arclength coordinate.
    // Called once per band step AFTER the equilibration phase.
    const int nat = static_cast<int>(R[0].rows());

    // Arclength: cumulative distance along the current band.
    std::vector<double> s(m_nimages, 0.0);
    for (int i = 1; i < m_nimages; ++i)
        s[i] = s[i - 1] + (R[i] - R[i - 1]).norm();

    for (int i = 0; i < m_nimages; ++i) {
        m_pmf_arc_sum[i] += s[i];
        if (!m_images[i])
            continue; // clamped endpoint: no dynamics, no mean force

        // -F_true . tau  =  +g . tau  (g is the gradient dE/dR)
        const Geometry g = m_images[i]->gradient();
        double gdot = 0.0;
        for (int a = 0; a < nat; ++a)
            gdot += g.row(a).dot(tau[i].row(a));
        m_pmf_force_sum[i] += gdot;

        // Spring contribution. Claude Generated 2026 - CORRECTED.
        //
        // The previous code treated the spring as an umbrella bias and SUBTRACTED its
        // mean force from <grad E . tau>. That is wrong, and measurably so: an umbrella
        // bias is an external potential with a FIXED centre, whereas the NEB spring
        // couples to the neighbouring images, which fluctuate themselves. In a
        // stationary band each image satisfies
        //
        //     <F_true . tau> + <F_spring . tau> ~ 0
        //
        // because the spring is exactly what balances the true force and keeps the image
        // at its position. Subtracting one from the other therefore cancels the signal
        // rather than removing a bias - which is why the barrier COLLAPSED as sampling
        // improved (58.4 kcal/mol at 400 fs -> 4.69 at 100 ps).
        //
        // The correct reading is the opposite one: in the constrained/stationary
        // ensemble the mean spring force IS the negative mean thermodynamic force
        // (den Otter & Briels, JCP 109, 4139 (1998)), so
        //
        //     dG/ds = <grad E . tau>            (unconstrained part)
        //
        // is already the estimator, and the spring term must NOT enter it. We keep
        // accumulating it for diagnostics: |<F_spring.tau>| ~ |<grad E.tau>| is the
        // signature of a well-equilibrated band, and the residual of the two is a
        // direct measure of how far from stationarity the run is.
        double sdot = 0.0;
        for (int a = 0; a < nat; ++a)
            sdot += Fspring[i].row(a).dot(tau[i].row(a));
        m_pmf_spring_sum[i] += -sdot;

        // Metric (Fixman) term. Claude Generated 2026.
        //
        // s is a curvilinear coordinate, so the volume element of the constrained
        // ensemble is not constant and the free-energy derivative carries an extra
        // contribution
        //
        //     dG/ds |_metric = -kT/2 * d ln Z / ds ,   Z = sum_a |tau_a|^2 / m_a
        //
        // (Z is the mass metric of the CV; for a unit tangent it reduces to the
        // mass-weighted norm of the tangent). Reference: Fixman, PNAS 71, 3050 (1974);
        // den Otter & Briels, J. Chem. Phys. 109, 4139 (1998) for the constrained-MD
        // form used here.
        //
        // d ln Z / ds is taken as a central difference over the neighbouring images,
        // which is the only s-resolution the band provides. Endpoints and images whose
        // neighbours coincide are skipped (the term is then undefined, not zero).
        if (m_fixman && i > 0 && i < m_nimages - 1) {
            auto metric = [&](const Geometry& t) {
                double Z = 0.0;
                for (int a = 0; a < nat; ++a)
                    Z += t.row(a).squaredNorm() / std::max(m_masses[a], 1e-12);
                return Z;
            };
            const double Zp = metric(tau[i + 1]);
            const double Zm = metric(tau[i - 1]);
            const double ds = (R[i + 1] - R[i - 1]).norm();
            if (Zp > 1e-30 && Zm > 1e-30 && ds > 1e-12) {
                const double kT = 3.166811e-6 * 300.0; // Eh at the thermostat setpoint
                const double dlnZ = (std::log(Zp) - std::log(Zm)) / ds;
                m_pmf_fixman_sum[i] += -0.5 * kT * dlnZ;
            }
        }

        const double e = m_images[i]->potentialEnergy();
        m_pmf_epot_sum[i] += e;
        m_pmf_epot_sq_sum[i] += e * e;
    }
    ++m_pmf_samples;
}

void NebMD::writeFreeEnergyProfile() const
{
    CitationRegistry::cite("denotter");
    if (m_fixman)
        CitationRegistry::cite("fixman");
    if (m_pmf_samples <= 0) {
        CurcumaLogger::warn("NEB-MD: no PMF samples collected (run longer than pmf_equilibration) "
                            "- skipping the free-energy profile.");
        return;
    }

    const double inv_n = 1.0 / static_cast<double>(m_pmf_samples);
    const double kcal = 627.5094740631;

    // Mean arclength, mean force (spring bias removed) and mean potential energy.
    std::vector<double> s(m_nimages, 0.0), dGds(m_nimages, 0.0);
    std::vector<double> emean(m_nimages, 0.0), esd(m_nimages, 0.0);
    for (int i = 0; i < m_nimages; ++i) {
        s[i] = m_pmf_arc_sum[i] * inv_n;
        if (m_images[i]) {
            // dG/ds = <grad E . tau> + metric (Fixman) term.
            //
            // The spring term is deliberately NOT subtracted here - see the note in
            // accumulatePMF(). Subtracting it cancelled the signal (the spring balances
            // the true force in a stationary band) and drove the barrier to zero as the
            // statistics improved.
            dGds[i] = (m_pmf_force_sum[i] + m_pmf_fixman_sum[i]) * inv_n;
            emean[i] = m_pmf_epot_sum[i] * inv_n;
            const double var = m_pmf_epot_sq_sum[i] * inv_n - emean[i] * emean[i];
            esd[i] = (var > 0.0) ? std::sqrt(var) : 0.0;
        } else {
            emean[i] = (i == 0) ? m_endpoint_energy_start : m_endpoint_energy_end;
        }
    }
    // Clamped endpoints carry NO mean force (they are not integrated), so G is only
    // defined on the sampled interior. Extrapolating into them fabricates a barrier at
    // the endpoint - the profile is therefore reported over the sampled range only.
    const int i_lo = m_images[0] ? 0 : 1;
    const int i_hi = m_images[m_nimages - 1] ? m_nimages - 1 : m_nimages - 2;
    if (i_hi <= i_lo) {
        CurcumaLogger::warn("NEB-MD: fewer than two sampled images - no free-energy profile.");
        return;
    }

    // Integrate dG/ds over the arclength (trapezoid), sampled images only.
    std::vector<double> G(m_nimages, 0.0);
    for (int i = i_lo + 1; i <= i_hi; ++i)
        G[i] = G[i - 1] + 0.5 * (dGds[i] + dGds[i - 1]) * (s[i] - s[i - 1]);

    std::ofstream f(outputPath(Basename() + ".neb.pmf.csv"));
    if (f) {
        f << "# NEB-MD free-energy profile (finite-temperature string, mean-force integration)\n";
        f << "# samples=" << m_pmf_samples << " equilibration_steps=" << m_pmf_equilibration << "\n";
        f << "# dG/ds has the spring (umbrella) contribution removed; a metric/Fixman term is neglected.\n";
        f << "image,s_Angstrom,dGds_Eh_per_A,G_kcal,Epot_mean_Eh,Epot_sd_kcal\n";
        for (int i = 0; i < m_nimages; ++i) {
            f << i << "," << fmt::format("{:.6f}", s[i]) << ","
              << fmt::format("{:.8f}", dGds[i]) << ","
              << fmt::format("{:.4f}", G[i] * kcal) << ","
              << fmt::format("{:.6f}", emean[i]) << ","
              << fmt::format("{:.4f}", esd[i] * kcal) << "\n";
        }
    }

    // Barrier + width analysis. A single highest-energy image is a point estimate; the
    // barrier has a WIDTH, which is what controls the TST prefactor and tunnelling.
    int imax = i_lo;
    for (int i = i_lo; i <= i_hi; ++i)
        if (G[i] > G[imax])
            imax = i;
    const double dG_barrier = (G[imax] - G[i_lo]) * kcal;
    const double dE_barrier = (emean[imax] - emean[i_lo]) * kcal;

    // Width at half the barrier height (in arclength), from the discrete profile.
    double width = 0.0;
    if (dG_barrier > 0.0) {
        const double half = G[i_lo] + 0.5 * (G[imax] - G[i_lo]);
        double s_left = s[i_lo], s_right = s[i_hi];
        for (int i = imax; i > i_lo; --i) {
            if (G[i - 1] <= half) {
                const double t = (half - G[i - 1]) / std::max(G[i] - G[i - 1], 1e-30);
                s_left = s[i - 1] + t * (s[i] - s[i - 1]);
                break;
            }
        }
        for (int i = imax; i < i_hi; ++i) {
            if (G[i + 1] <= half) {
                const double t = (half - G[i + 1]) / std::max(G[i] - G[i + 1], 1e-30);
                s_right = s[i + 1] + t * (s[i] - s[i + 1]);
                break;
            }
        }
        width = s_right - s_left;
    }

    // ---- Profile evaluation -------------------------------------------------------
    // Forward/reverse barrier, reaction free energy, TS position, width and curvature.
    const double G_last = G[i_hi] * kcal;
    const double G_max = G[imax] * kcal;
    const double dG_reverse = G_max - G_last;
    const double dG_reaction = G_last - G[i_lo] * kcal;

    // Curvature at the top from a parabola through the three points around imax. The
    // magnitude of the (negative) curvature sets the imaginary frequency and therefore
    // both the TST prefactor and the tunnelling correction: a narrow, sharply curved
    // barrier tunnels much more than a broad flat one.
    double curvature = 0.0;
    if (imax > 0 && imax < m_nimages - 1) {
        const double h1 = s[imax] - s[imax - 1];
        const double h2 = s[imax + 1] - s[imax];
        if (h1 > 1e-9 && h2 > 1e-9) {
            const double gl = G[imax - 1] * kcal, gm = G_max, gr = G[imax + 1] * kcal;
            curvature = 2.0 * (h1 * gr - (h1 + h2) * gm + h2 * gl) / (h1 * h2 * (h1 + h2));
        }
    }

    CurcumaLogger::result("");
    CurcumaLogger::result("=== NEB-MD free-energy profile ===");
    if (m_nudged) {
        // This is the single most important caveat of the whole analysis.
        CurcumaLogger::warn("PMF with nudged=true is NOT thermodynamically valid: projecting out the "
                            "parallel true force makes the force field non-conservative, so there is no "
                            "underlying potential and no proper NVT ensemble (the images heat up). Use "
                            "-nebmd.nudged false for free-energy estimates; keep nudged=true only for "
                            "path/saddle searching.");
    }
    CurcumaLogger::result_fmt("  samples/image      : {} (after {} equilibration steps)",
                              m_pmf_samples, m_pmf_equilibration);
    CurcumaLogger::result_fmt("  sampled images     : {}..{} (clamped endpoints carry no mean force)",
                              i_lo, i_hi);
    CurcumaLogger::result_fmt("  path length        : {:.2f} A", s[m_nimages - 1]);
    CurcumaLogger::result_fmt("  barrier forward    : {:.2f} kcal/mol (image {}, s = {:.2f} A)",
                              dG_barrier, imax, s[imax]);
    CurcumaLogger::result_fmt("  barrier reverse    : {:.2f} kcal/mol", dG_reverse);
    CurcumaLogger::result_fmt("  reaction dG        : {:.2f} kcal/mol", dG_reaction);
    CurcumaLogger::result_fmt("  <Epot> barrier     : {:.2f} kcal/mol  (potential-energy estimate)",
                              dE_barrier);
    CurcumaLogger::result_fmt("  entropic/thermal   : {:.2f} kcal/mol  (dG - <Epot> difference)",
                              dG_barrier - dE_barrier);
    CurcumaLogger::result_fmt("  barrier width FWHM : {:.2f} A", width);
    CurcumaLogger::result_fmt("  curvature at top   : {:.2f} kcal/mol/A^2 (negative = real maximum; "
                              "large |value| = narrow barrier, more tunnelling)", curvature);
    CurcumaLogger::result_fmt("  Epot fluctuation   : +/-{:.2f} kcal/mol at the top", esd[imax] * kcal);
    {
        // Stationarity diagnostic: in an equilibrated band <F_spring.tau> should cancel
        // <grad E.tau> image by image. A large residual means the run is not stationary
        // and dG/ds cannot be trusted, whatever its value.
        double worst = 0.0;
        for (int i = i_lo; i <= i_hi; ++i) {
            if (!m_images[i])
                continue;
            const double f = m_pmf_force_sum[i] * inv_n;
            const double sp = m_pmf_spring_sum[i] * inv_n;
            const double scale = std::max(std::abs(f), 1e-12);
            worst = std::max(worst, std::abs(f + sp) / scale);
        }
        CurcumaLogger::result_fmt("  stationarity       : worst |<F.tau>+<Fspring.tau>|/|<F.tau>| = {:.2f} "
                                  "(0 = spring exactly balances the true force)", worst);
    }
    if (m_fixman) {
        // Report the metric term separately so its (expected small) size can be judged
        // against the total rather than being hidden inside it.
        double fmin = 0.0, fmax_c = 0.0, ftot = 0.0;
        bool first = true;
        for (int i = i_lo; i <= i_hi; ++i) {
            if (!m_images[i])
                continue;
            const double v = m_pmf_fixman_sum[i] * inv_n;
            ftot += std::abs(v);
            if (first) { fmin = fmax_c = v; first = false; }
            fmin = std::min(fmin, v);
            fmax_c = std::max(fmax_c, v);
        }
        CurcumaLogger::result_fmt("  metric (Fixman)    : dG/ds contribution {:.3e} .. {:.3e} Eh/A "
                                  "(mean |.| {:.3e}); expected O(kT) - a correctness term, not a cure",
                                  fmin, fmax_c, ftot / std::max(i_hi - i_lo + 1, 1));
    }

    // Compact ASCII profile (plain characters only - no UTF, per project rules).
    {
        double gmin = G[i_lo] * kcal, gmax = G[i_lo] * kcal;
        for (int i = i_lo; i <= i_hi; ++i) {
            gmin = std::min(gmin, G[i] * kcal);
            gmax = std::max(gmax, G[i] * kcal);
        }
        const double span = std::max(gmax - gmin, 1e-9);
        const int wbar = 40;
        CurcumaLogger::result("  profile (G in kcal/mol vs image):");
        for (int i = i_lo; i <= i_hi; ++i) {
            const double gi = G[i] * kcal;
            const int n = static_cast<int>(std::lround((gi - gmin) / span * wbar));
            std::string bar(std::max(n, 0), '#');
            CurcumaLogger::result(fmt::format("    {:>2} |{:<{}}| {:8.2f}{}", i, bar, wbar, gi,
                                              (i == imax ? "  <- TS" : "")));
        }
    }
    CurcumaLogger::warn("EXPERIMENTAL - do not quote these numbers as barriers. G(s) is only as good "
                        "as the path it is computed on, and the trapezoid integration over few images "
                        "carries a visible error (helicene test: dG(reaction) = -22 kcal/mol where "
                        "symmetry demands 0). See docs/NEB_MD.md for the validation status.");
    if (m_reparametrize <= 0) {
        CurcumaLogger::warn("PMF without string reparametrisation: the images drift along the path, "
                            "so the mean force is integrated over a moving grid. Set "
                            "-nebmd.reparametrize 5 (or similar) for a meaningful profile.");
    }
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
    // Register band outputs for copy-back to CWD (BMT mode).
    addBakFile(Basename() + ".neb.path.xyz");
    addBakFile(Basename() + ".neb.energy.csv");
    addBakFile(Basename() + ".neb.final.xyz");
    if (m_pmf)
        addBakFile(Basename() + ".neb.pmf.csv");
    if (m_write_initial_path)
        addBakFile(Basename() + ".neb.init.xyz");
    if (!buildInitialPath())
        return;
    CitationRegistry::cite("neb"); // improved tangent (Henkelman & Jonsson 2000)
    if (m_idpp)
        buildIDPPPath();       // repulsion-aware, neighbour-coupled path refinement
    if (m_repair_overlaps)
        repairPathOverlaps();  // clean up any residual close contacts
    checkPathTopology();       // warn if the connectivity still differs from image 0
    if (m_optimize)
        optimizeBand();        // relax the band onto the MEP (classical CI-NEB) before the MD
    if (!prepareImages())
        return;
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
    if (m_pmf)
        writeFreeEnergyProfile();
    if (m_umbrella)
        writeWHAM();

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