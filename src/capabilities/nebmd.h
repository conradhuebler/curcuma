/*
 * <NEB-like Molecular Dynamics: a chain of independently-thermostatted replicas
 *    coupled by spring forces along a reaction path. Each image runs its own MD
 *    (velocities + thermostat); the band is driven by MD, not by LBFGS.>
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
 */
/*
 * NEB-MD concept
 * --------------
 * N images (replicas) are placed along an initial path between a start and an
 * end structure (RMSD-reordered end -> linear interpolation). Each INTERIOR
 * image is driven by its own SimpleMD instance (independent velocities,
 * thermostat, kinetic energy). The images are coupled by spring forces to
 * their immediate neighbours so they move constrained along the path instead
 * of freely relaxing to the nearest minimum.
 *
 * Force model (per interior image i, atom a):
 *   g_i        = true gradient at R_i (Eh/Å, dE/dR)            [from image.gradient()]
 *   tau_i      = improved tangent (Henkelman-Jónsson 2000)     [unit-vector, per atom]
 *   F_spring   = k_a * (R_{i+1} + R_{i-1} - 2 R_i)             [full harmonic]
 *
 *   nudged   (default): external = (g_i·tau_i) tau_i + (F_spring·tau_i) tau_i
 *                       (SimpleMD already applies F_true = -g_i, so we inject the
 *                        negative parallel true-force component plus the parallel
 *                        spring => total F_true_perp + F_spring_parallel)
 *   non-nudged        : external = F_spring
 *                       (total = F_true + F_spring, full components)
 *
 * The external force is injected via SimpleMD::applyExternalForces() (one
 * assignment per step, units = gradient units = Eh/Å).
 *
 * Timing caveat: SimpleMD::Verlet() overwrites m_eigen_gradient with the
 * gradient at the NEW position mid-step (Energy() re-eval), so the injected
 * spring acts on the position update + first half-kick only, not on the
 * second half-kick. This breaks strict velocity-Verlet time symmetry by a
 * small, smooth perturbation - acceptable for this nonequilibrium
 * path-relaxation use with a thermostat. See docs/NEB_MD.md.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "curcumamethod.h"
#include "src/core/parameter_macros.h"     // PARAM block (build-time extracted)
#include "src/core/parameter_registry.h"    // getDefaultJson("nebmd")
#include "src/core/config_manager.h"        // type-safe parameter access
#include "src/core/molecule.h"
#include "src/core/elements.h"              // Elements::AtomicMass
#include "src/core/energycalculator.h"      // endpoint single-point energies
#include "src/capabilities/simplemd.h"      // one instance per image
#include "src/capabilities/rmsd.h"          // RMSDDriver: reorder end -> start
#include "src/capabilities/path_cv.h"        // Branduardi path collective variable

class NebMD : public CurcumaMethod {
public:
    NebMD(const json& controller = json(), bool silent = true);
    ~NebMD();

    /*! \brief Set the two endpoint structures (start, end). Atom counts / element
     *         multisets must match; the end is RMSD-reordered onto the start. */
    void setEndpoints(const Molecule& start, const Molecule& end);

    bool Initialise() override;
    void start() override;
    void printHelp() const override;

private:
    // --- path construction ---
    bool buildInitialPath();                 // RMSD-reorder end, linear-interpolate N images
    // Refine the linearly interpolated path with the Image-Dependent Pair Potential
    // (IDPP, Smidstrup et al., J. Chem. Phys. 140, 214106 (2014)): interpolate the
    // DISTANCE matrix between the endpoints and optimise every image against it with a
    // 1/d^4 weight, which penalises short contacts hard -> the path is repulsion-aware.
    // Run as a NEB on the IDPP surface, so each image also feels its neighbours
    // (spring force + tangent nudging). Returns false if the objective did not improve.
    bool buildIDPPPath();
    // IDPP objective S = sum_{a<b} (d_target - d)^2 / d^4 and its cartesian gradient.
    double idppObjective(const Geometry& R, const Matrix& d_target, Geometry* grad) const;
    // Check that every interpolated image has the same bond connectivity as image 0.
    // Linear interpolation can stretch a bond past the detection threshold at an
    // intermediate image -> that image gets a different GFN-FF topology (discontinuous
    // energy/gradient along the band). Advisory warn (no abort); returns false on mismatch.
    bool checkPathTopology() const;
    // Push apart non-bonded atom pairs that linear interpolation drove into each other
    // (steepest-descent on a soft (r_cut/r)^4 clash potential, same scheme as
    // PolymerBuild::resolveOverlaps). Only the interior images are touched; the endpoints
    // keep their input geometry. Returns the number of images that were repaired.
    int repairPathOverlaps();
    Geometry interpolate(const Geometry& a, const Geometry& b, double t) const;
    bool prepareImages();                    // construct N SimpleMD, setMolecule, Initialise, prepareRun
    void computeEndpointEnergies();          // single-point energies for fixed endpoints (improved tangent)

    // --- per-step coordinator ---
    std::vector<Geometry> currentPath() const;     // positions of all images (incl. fixed endpoints)
    std::vector<double> currentEnergies() const;   // Epot of all images (incl. fixed endpoints)
    // Improved tangent (Henkelman-Jónsson 2000), energy-weighted; one unit vector per image.
    std::vector<Geometry> computeTangents(const std::vector<Geometry>& R,
                                          const std::vector<double>& E) const;
    // Spring force per interior image (parallel component when nudged, full when not).
    std::vector<Geometry> computeSpringForces(const std::vector<Geometry>& R,
                                              const std::vector<Geometry>& tangents) const;
    bool stepBand();                         // one NEB-MD step; returns false when all interior images finished

    /* --- Free-energy profile G(s) (Claude Generated 2026) -------------------------
     * The band is a finite-temperature string: each image is an independently
     * thermostatted replica held near its slice of the path by the neighbour springs.
     * That is umbrella-LIKE (each image = a window) but NOT textbook umbrella sampling:
     * the restraint acts between neighbouring images in full 3N space and its centre
     * moves, so WHAM does not apply. Instead we accumulate the mean force along the
     * tangent and integrate it (umbrella-integration style):
     *
     *     dG/ds = -<F_true . tau> - <F_spring . tau>     (spring bias removed)
     *     G(s_i) = G(s_0) + integral_0^{s_i} dG/ds ds    (trapezoid over arclength)
     *
     * Caveats (documented, not hidden): the images must be equilibrated before the
     * averages mean anything (m_pmf_equilibration steps are discarded), the spring
     * removal assumes the restraint is harmonic and locally linear, and a metric
     * (Fixman) term is neglected - small for cartesian coupling but not exactly zero.
     */
    /*! \brief Optimise the band to the minimum energy path (classical NEB) before any MD.
     *
     *  Uses FIRE (Bitzek et al., PRL 97, 170201 (2006)) rather than L-BFGS: the nudged
     *  force is NOT the gradient of any potential (the parallel true force is projected
     *  out), so a line search - which needs an energy to compare against - is not
     *  applicable. FIRE only ever uses the force, which makes it the standard choice for
     *  NEB. Optionally promotes the highest image to a climbing image (CI-NEB,
     *  Henkelman, Uberuaga & Jonsson, JCP 113, 9901 (2000)) so it converges onto the
     *  saddle point instead of being held back by the springs.
     *  \return true when the force criterion was met before max_iterations. */
    bool optimizeBand();
    /*! \brief Energy + gradient of one geometry via a shared EnergyCalculator. */
    double bandEnergy(EnergyCalculator& ec, const Geometry& geom, Geometry* grad) const;

    /*! \brief Redistribute the images equidistantly along the current path (string
     *  reparametrisation, Vanden-Eijnden & Venturoli 2009 / E, Ren & Vanden-Eijnden).
     *
     *  This is the piece that makes the band a well-defined set of umbrella windows: the
     *  MD lets the images drift along the path (they bunch up in the minima and thin out
     *  over the barrier), so the arclength coordinate keeps changing under the integral.
     *  Reparametrisation restores a uniform spacing by cubic-Hermite interpolation of the
     *  polyline through the images, re-sampled at equal arclength.
     *
     *  Endpoints are never moved. Velocities are rescaled, not copied, so the image
     *  temperature is preserved across the (discontinuous) coordinate change.
     *  \return the RMS displacement applied [Angstrom] (0 if the path was already uniform). */
    double reparametrizeString();

    /*! \brief Umbrella sampling on the path collective variable (Branduardi s).
     *
     *  Replaces the neighbour springs by a harmonic restraint on a FROZEN reference
     *  path: image i is held at s_0,i = i/(N-1) by V = k/2 (s(R) - s_0,i)^2. Because
     *  the reference frames do not move, every window is stationary - which is exactly
     *  what the arclength coupling fails to provide - so the histograms can be combined
     *  with WHAM into a proper G(s). */
    void setupUmbrella();                      ///< freeze the path, build the PathCV
    void applyUmbrellaBias(int i, SimpleMD* img, Geometry& Fext) const; ///< add dV/dR
    void writeWHAM() const;                    ///< unbias the histograms -> G(s)

    void accumulatePMF(const std::vector<Geometry>& R,
                       const std::vector<Geometry>& tau,
                       const std::vector<Geometry>& Fspring);
    void writeFreeEnergyProfile() const;     // -> <basename>.neb.pmf.csv + barrier summary

    /* --- Transverse RMSD metadynamics (Claude Generated 2026) ---------------------
     * A plain RMSD-MTD bias would push each image out of its own basin ALONG the path,
     * destroying the string. Projecting the bias force onto the plane perpendicular to
     * the tangent lets an image explore its transverse degrees of freedom (the entropic
     * part of the barrier) while the path direction stays governed by the springs.
     */
    Geometry projectOutTangent(const Geometry& F, const Geometry& tau) const;

    // --- output (routed through outputPath() / BMT) ---
    void writeInitialPath() const;           // -> <basename>.neb.init.xyz
    void writePathSnapshot(int step) const;  // append all images -> <basename>.neb.path.xyz
    void writeEnergyTable(int step) const;   // append rows to <basename>.neb.energy.csv
    void writeFinalPath() const;             // -> <basename>.neb.final.xyz
    // Shared writer: append all N image geometries as a multi-XYZ block to `filename`,
    // each frame commented with "tag step=.. image=.. E=.. T=..".
    void writeBandFile(const std::string& filename, int step, const std::string& tag) const;

    // --- CurcumaMethod contract ---
    nlohmann::json WriteRestartInformation() override;
    bool LoadRestartInformation() override;
    void LoadControlJson() override;
    void ReadControlFile() override {}
    StringList MethodName() const override { return { "NEBMD" }; }

    // --- members ---
    Molecule m_start, m_end;
    Geometry m_start_geom, m_end_geom;       // aligned / reordered, in start frame (Å)
    std::vector<double> m_masses;            // per-atom mass (amu), from m_start.Atoms()
    std::vector<SimpleMD*> m_images;         // size N; nullptr for fixed endpoints
    std::vector<bool> m_image_done;          // interior images whose step() returned false
    std::vector<Geometry> m_image_geoms;     // current path (fixed-endpoint positions + restart)

    // PMF accumulators (index = image). Running sums over the production steps.
    std::vector<double> m_pmf_force_sum;     // sum of -(F_true . tau)  [Eh/Angstrom]
    std::vector<double> m_pmf_spring_sum;    // sum of -(F_spring . tau) (bias, subtracted)
    std::vector<double> m_pmf_fixman_sum;    // sum of the metric (Fixman) contribution
    std::vector<double> m_pmf_epot_sum;      // sum of E_pot [Eh]
    std::vector<double> m_pmf_epot_sq_sum;   // sum of E_pot^2 (for the fluctuation width)
    std::vector<double> m_pmf_arc_sum;       // sum of the arclength coordinate s_i [Angstrom]
    long m_pmf_samples = 0;                  // number of accumulated production steps

    double m_endpoint_energy_start = 0.0;    // SP energy of fixed start endpoint
    double m_endpoint_energy_end = 0.0;      // SP energy of fixed end   endpoint
    bool m_endpoints_set = false;
    int m_step = 0;

    // --- parameters (loaded in LoadControlJson) ---
    int m_nimages = 8;
    bool m_nudged = true;
    bool m_fixed_endpoints = true;
    bool m_mass_weighted_spring = true;
    double m_k_spring = 0.005;               // Eh/Å² (or Eh/(amu·Å²) when mass-weighted)
    int m_dump = 50;                         // path-snapshot / energy-table cadence (band steps)
    bool m_write_initial_path = true;
    bool m_optimize = false;                 // run a classical NEB optimisation before the MD
    int m_opt_iterations = 300;              // FIRE iterations for the band optimisation
    double m_opt_fmax = 0.02;                // convergence: max force per atom [Eh/Angstrom]
    bool m_climbing = true;                  // promote the highest image (CI-NEB)
    double m_opt_k = 0.1;                    // NEB spring constant during optimisation [Eh/A^2]
    // Umbrella sampling on the path CV (Claude Generated 2026)
    std::string m_path_file;                 // external reference path (multi-XYZ)
    bool m_umbrella = false;                 // use path-CV umbrella instead of neighbour springs
    double m_umbrella_k = 0.0;               // restraint strength (0 = auto, see setupUmbrella)
    double m_kz = 0.0;                       // tube restraint on the perpendicular CV z (0 = off)
    double m_z_tolerance = 0.05;             // free tube radius before the z restraint bites
    bool m_umbrella_keep_springs = false;    // also keep the neighbour springs in umbrella mode
    int m_umbrella_bins = 100;               // histogram bins over s in [0,1]
    std::unique_ptr<PathCV> m_path_cv;       // frozen reference path
    std::vector<double> m_umbrella_s0;       // window centres, one per image
    std::vector<std::vector<int>> m_hist;    // per-image histogram over s
    std::vector<double> m_umbrella_ssum;     // running <s> per image (diagnostics)
    std::vector<long> m_umbrella_n;          // samples per image

    int m_reparametrize = 0;                 // reparametrise the string every N steps (0 = off)
    bool m_pmf = false;                      // experimental mean-force profile (opt-in, not quantitative)
    bool m_fixman = true;                    // include the metric (Fixman) term in dG/ds
    int m_pmf_equilibration = 200;           // band steps discarded before averaging
    bool m_mtd_transverse = true;            // project the RMSD-MTD bias perpendicular to tau
    bool m_idpp = true;                      // IDPP refinement of the interpolated path
    int m_idpp_iterations = 200;             // steepest-descent iterations for IDPP
    bool m_repair_overlaps = true;           // push apart clashing atoms in interpolated images
    int m_pt = 0;                            // proton-transfer limit for RMSD reorder
    bool m_force_reorder = true;
    int m_threads = 1;
    std::string m_method = "uff";            // energy method (forwarded to each image)
    json m_md_config;                        // controller["simplemd"] sub-object for the images

    ConfigManager m_config;                  // "nebmd" parameter access

    // vvvvvvvvvvvv PARAMETER DEFINITION BLOCK vvvvvvvvvvvv
    // Claude Generated 2026 - Parameter Registry Integration
    BEGIN_PARAMETER_DEFINITION(nebmd)

    // --- Band ---
    PARAM(nimages, Int, 8, "Number of replicas along the path (including endpoints).", "Band", {"N"})
    PARAM(nudged, Bool, true, "Nudged NEB projection: true force perpendicular + spring parallel. false = full spring + full true force (plain elastic band).", "Band", {})
    PARAM(fixed_endpoints, Bool, true, "Clamp endpoint images (do not integrate them). false = integrate endpoints as full images.", "Band", {})
    PARAM(k_spring, Double, 0.005, "Spring constant in Eh/Å² (mass-weighted: Eh/(amu·Å²)). Raise for a stiffer band, lower for looser coupling.", "Band", {"k"})
    PARAM(mass_weighted_spring, Bool, true, "Scale k by per-atom mass so the spring frequency is mass-independent (MD-friendly; keeps H atoms stable).", "Band", {})
    // --- Band optimisation (classical NEB) ---
    PARAM(optimize, Bool, false, "Optimise the band to the minimum energy path with a classical NEB (FIRE optimiser) before starting the MD. Gives the potential barrier dE and a relaxed path for the subsequent sampling.", "Optimisation", {"neb"})
    PARAM(opt_iterations, Int, 300, "Maximum FIRE iterations for the NEB band optimisation.", "Optimisation", {})
    PARAM(opt_fmax, Double, 0.02, "Convergence threshold for the NEB optimisation: maximum force per atom in Eh/Angstrom.", "Optimisation", {})
    PARAM(climbing, Bool, true, "Climbing-image NEB: the highest image inverts its parallel force so it converges onto the saddle point instead of being pulled back by the springs.", "Optimisation", {})
    PARAM(opt_k, Double, 0.1, "Spring constant used during the NEB optimisation (Eh/Angstrom^2). Independent of the MD spring k_spring.", "Optimisation", {})

    // --- Free energy / sampling ---
    PARAM(path_file, String, "", "Read the reference path from a multi-XYZ file (e.g. an ORCA neb_MEP_trj.xyz) instead of interpolating between the endpoints. The frame count then sets nimages. Use this to sample on an externally converged MEP.", "Band", {})
    PARAM(umbrella, Bool, false, "Umbrella sampling on the path collective variable (Branduardi s) instead of neighbour springs. The reference path is frozen, so the windows are stationary and WHAM applies - this is the statistically correct route to G(s).", "FreeEnergy", {})
    PARAM(umbrella_k, Double, 0.0, "Force constant of the harmonic restraint on the path CV s (Eh per unit s^2). 0 = choose automatically so neighbouring windows overlap (sigma ~ half the window spacing), which is what WHAM requires.", "FreeEnergy", {})
    PARAM(umbrella_kz, Double, 0.0, "Force constant of a tube restraint on the PERPENDICULAR path CV z (distance from the reference path). 0 = off. Keeps the sampling on the path; without it images were measured to drift up to 1.1 A off a converged MEP.", "FreeEnergy", {})
    PARAM(umbrella_z_tolerance, Double, 0.05, "Free tube radius: the z restraint only acts beyond this distance from the path, so the image can fluctuate inside the tube.", "FreeEnergy", {})
    PARAM(umbrella_keep_springs, Bool, false, "Keep the neighbour springs in umbrella mode so an image still feels its neighbours. Helps it stay on the path, but the spring bias is NOT removed by WHAM - diagnostics, not production free energies.", "FreeEnergy", {})
    PARAM(umbrella_bins, Int, 100, "Number of histogram bins over s in [0,1] for the WHAM analysis.", "FreeEnergy", {})
    PARAM(reparametrize, Int, 0, "Redistribute the images to equal arclength every N band steps (string reparametrisation). 0 = off. Required for a meaningful free-energy profile; a value of 1-10 keeps the images as well-defined umbrella windows.", "FreeEnergy", {})
    PARAM(pmf, Bool, false, "EXPERIMENTAL and currently NOT quantitative: accumulate the mean force along the tangent and integrate it to G(s). Requires string reparametrisation (not implemented), so the absolute numbers are unreliable - see docs/NEB_MD.md. Written to <base>.neb.pmf.csv.", "FreeEnergy", {})
    PARAM(fixman, Bool, true, "Include the metric (Fixman) term -kT/2 dlnZ/ds in dG/ds. s is a curvilinear coordinate, so the constrained-ensemble volume element is not constant. Expected magnitude is of order kT (~0.6 kcal/mol at 300 K) - it is a correctness term, not a fix for large errors.", "FreeEnergy", {})
    PARAM(pmf_equilibration, Int, 200, "Band steps discarded before the mean-force averaging starts (equilibration).", "FreeEnergy", {})
    PARAM(mtd_transverse, Bool, true, "When RMSD metadynamics is enabled (-simplemd.rmsd_mtd true), project its bias force perpendicular to the path tangent so it samples transverse degrees of freedom without destroying the string.", "FreeEnergy", {})

    PARAM(idpp, Bool, true, "Refine the interpolated path with the Image-Dependent Pair Potential (Smidstrup 2014): interpolate the distance matrix and optimise each image with a 1/d^4 weight, so short contacts are resolved and neighbours are accounted for.", "Band", {})
    PARAM(idpp_iterations, Int, 200, "Steepest-descent iterations for the IDPP path refinement.", "Band", {})
    PARAM(repair_overlaps, Bool, true, "Push apart non-bonded atoms that remain too close after interpolation (soft clash relaxation, as in polymerbuild). Runs after the IDPP refinement.", "Band", {})

    // --- Output ---
    PARAM(dump_frequency, Int, 50, "Write a path snapshot and append to the energy table every N band steps.", "Output", {"dump"})
    PARAM(write_initial_path, Bool, true, "Write the initial interpolated path as a multi-XYZ before integration.", "Output", {})

    // --- RMSD endpoint reorder ---
    PARAM(proton_transfer, Int, 0, "Allowed proton-transfer count for the RMSD endpoint reorder.", "RMSD", {"pt"})
    PARAM(force_reorder, Bool, true, "Force atom reordering of the end structure onto the start atom order.", "RMSD", {})

    // --- System ---
    PARAM(threads, Int, 1, "Threads per image (keep 1: gfnff CitationRegistry race at >1).", "System", {})

    END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^
};