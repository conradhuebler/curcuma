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

#include <fstream>
#include <memory>
#include <random>
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

    /*! \brief True when the run ended in an error state that must not be consumed further
     *         (currently: string stage aborted on a dead image, string_abort_on_dead).
     *         The caller turns it into a non-zero exit code and skips pipeline stage 2.
     *         Claude Generated 2026. */
    bool failed() const { return m_string_aborted; }

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

    /* Scheme 'string': finite-temperature string, stage 1 of the string -> path-CV chain.
     * Claude Generated 2026 (AI-generated, machine-tested).
     *
     * Each live image i is restrained to the CURRENT mean-path frame M_i only - no
     * neighbour springs (Maragliano & Vanden-Eijnden, CPL 426, 168 (2006)):
     *
     *   V_i(R) = k/2 sum_a w_a |R_a - (Rot* M_i,a + t*)|^2 ,  w_a = m_a (mass_weighted_spring)
     *
     * with (Rot*, t*) the w-weighted Kabsch minimisers of the SAME sum, so by the envelope
     * theorem dV/dR_a = k w_a (R_a - M_aligned_a) exactly; V is invariant under rigid motion
     * of R (images rotate freely, remove_com_mode = 0). Between two updates M is FIXED, so
     * every image samples a canonical ensemble with a known static bias. Every string_update
     * steps: M_i <- (1-alpha) M_i + alpha <R_i aligned onto M_i>, optional Laplacian
     * smoothing, equal-(mass-weighted)-arclength reparametrisation of M only (the MD is never
     * touched: no setPositions, no velocity rescale). Fixed endpoints are never updated.
     * Convergence (Claude Generated 2026, force criterion, string_ftol > 0): per update the
     * perpendicular mean force F_ia = k w_a perp_ia (perp = part of <R_i> - M_i normal to the
     * path) and its standard error from string_blocks block means (string_path::updateMeanPath);
     * converged when in string_patience consecutive updates max_i F_sig < string_ftol AND
     * max_i SE < string_ftol AND scale == 1 AND no dead image AND delta_perp_max < string_tol
     * (at least kStringMinUpdates updates). SE >= ftol = 'noise-limited' (not an error: the
     * stage ends at string_max_time). string_ftol = 0 restores the former criterion
     * Delta_max < string_tol (false-positive prone: Delta_max ~ alpha F/(k m) is small for
     * any stiff restraint). Optional gain string_gain on the perpendicular part.
     * Start path preparation (string_prepare): before the images are built, the path is
     * consecutively Kabsch-aligned and reparametrised to equal arclength, so update 1 does not
     * have to absorb the reparametrisation.
     * Robustness (Claude Generated 2026): an update larger than string_max_step is scaled by
     * one common factor (string_path::limitStep, trust radius); a scaled update or one with a
     * dead image does not count towards convergence; a dead image aborts the stage by default
     * (string_abort_on_dead).
     *
     * The result is a k-dependent REFERENCE path for the path CV, not an exact MFEP. */
    static constexpr int kStringMinUpdates = 5;  // updates before convergence may be declared
    bool m_string = false;                       // m_scheme == "string"
    double m_string_k = 0.05;                    // Eh/(amu A^2) (mass-weighted) or Eh/A^2
    int m_string_update = 2000;                  // MD steps between mean-path updates
    double m_string_alpha = 0.5;                 // mixing of the new mean into M
    double m_string_smooth = 0.0;                // Laplacian smoothing kappa per update
    double m_string_tol = 0.03;                  // convergence threshold on Delta_max [A]
    int m_string_patience = 3;                   // consecutive updates below string_tol
    double m_string_max_time = 20000.0;          // MD length of stage 1 [fs]
    std::vector<Geometry> m_string_path;         // mean path M_i (common frame)
    std::vector<Geometry> m_string_sum;          // running sum of R_i aligned onto M_i
    std::vector<long> m_string_n;                // samples in the running sum
    std::vector<double> m_string_w;              // per-atom metric weights (masses or empty)
    int m_string_updates = 0;                    // updates performed
    int m_string_converged_count = 0;            // consecutive updates below string_tol
    bool m_string_converged = false;
    double m_string_last_delta = -1.0;           // Delta_max of the last update [A]
    mutable std::vector<long> m_string_capped;   // restraint-force cap hits per image
    // Robustness (Claude Generated 2026): trust radius on an update and the dead-image policy.
    // An update counts as "below tol" only if it was NOT scaled (scale == 1) and no image is dead.
    double m_string_max_step = 0.05;             // trust radius of one update [A] (0 = off)
    bool m_string_abort_on_dead = true;          // abort stage 1 when an image is declared dead
    bool m_string_aborted = false;               // stage 1 aborted because of a dead image
    double m_string_last_scale = 1.0;            // trust-radius factor s of the last update
    // Force criterion, gain, start-path preparation (Claude Generated 2026)
    bool m_string_prepare = true;                // align + reparametrise the start path
    double m_string_ftol = 5e-4;                 // force threshold [Eh/A RMS per atom]; 0 = old criterion
    int m_string_blocks = 5;                     // blocks per update for the standard error
    double m_string_gain = 1.0;                  // perpendicular gain (1 = off)
    double m_string_gain_z = 5.0;                // gain only for z >= gain_z (and F_sig >= ftol)
    double m_string_prep_shift = 0.0;            // max best-fit shift of the start-path preparation [A]
    bool m_string_write_updates = false;         // append M after every update to <base>.neb.string_updates.xyz
    std::vector<std::vector<Geometry>> m_string_block_sum; // [image][block] aligned sums
    std::vector<std::vector<long>> m_string_block_n;        // [image][block] samples
    // diagnostics of the last update (for the end message)
    double m_string_last_fsig = 0.0, m_string_last_se = 0.0, m_string_last_z = 0.0;

    void prepareStringStartPath();                                    // Claude Generated 2026

    void setupString();                                               // Claude Generated 2026
    void applyStringRestraint(int i, SimpleMD* img, Geometry& Fext) const; // Claude Generated 2026
    void accumulateStringMean(int i);                                 // Claude Generated 2026
    bool updateStringPath();                                          // Claude Generated 2026
    void writeStringPath(const std::string& filename) const;          // Claude Generated 2026

    /*! \brief Reference free energy G_ref(s) from a free trajectory (Claude Generated 2026).
     *
     *  Reads a multi-XYZ trajectory, evaluates s,z of every frame on the reference path,
     *  discards frames with z above the threshold (cv_hist_zmax, see PARAM), histograms s
     *  into fe_bins bins and writes G_ref = -kT ln P (kcal/mol, minimum 0), a standard error
     *  from 5 contiguous blocks, and the interior range [s(frame 0), s(frame N-1)] to
     *  <base>.neb.ref.csv. No MD is run. */
    void histogramTrajectory(const std::string& traj) const;
    std::string m_cv_histogram;                  // trajectory for histogramTrajectory ("" = off)
    double m_cv_hist_zmax = 0.0;                 // z filter (0 = automatic)

    /* --- Path CV decoupled from the windows (cv_frames / cv_extend). Claude Generated 2026 ---
     *
     * Default (cv_frames = 0, cv_extend = 0): the path CV is built on the image frames, one
     * window per frame, endpoints fixed - exactly the former behaviour (m_cv_path stays empty).
     *
     * Otherwise (scheme 'independent' and cv_histogram only; ignored elsewhere):
     *  1. P = path frames (path_file / interpolation / NEB result), consecutively Kabsch-aligned
     *     (string_path::alignConsecutive, weights m_cv_w: uniform, or 0/1 on the cv_atoms selection).
     *  2. cv_frames > 0: P is re-sampled at equal arclength onto cv_frames frames
     *     (string_path::resamplePath; interpolates or decimates, endpoints kept).
     *  3. cv_extend = m > 0: m frames are added beyond each end by straight continuation of the
     *     end segment (string_path::extendPath): F_{-j} = R_0 + j (R_0 - R_1'), F_{N-1+j} =
     *     R_{N-1} + j (R_{N-1} - R_{N-2}'), primes = Kabsch-aligned onto the neighbour. The SP
     *     energy of every extension frame is compared with its endpoint (warning if lower).
     *     cv_extend_length = f > 0 (takes precedence over cv_extend): m is derived from a LENGTH,
     *     m = string_path::extensionFrames(f, F) = ceil(f (F-1)) + 1 (F = frames after step 2), so
     *     the windows reach >= f L beyond each end independent of cv_frames (Claude Generated 2026).
     *     Steps 1-3 measure RMSD and arclength over the cv_atoms selection, the metric of the CV.
     *  4. The PathCV (lambda from the plain spacings) is built on the resulting frames C_0..C_{M-1}.
     *     Window range: CV frames lo..hi with lo = 1 if m > 0 else 0, hi = M-1-lo, i.e. only the
     *     OUTERMOST extension frame per side is excluded (it is the plateau of s).
     *  5. Windows: cv_frames = 0 -> one window per CV frame lo..hi (start geometry = that frame).
     *     cv_frames > 0 -> the window count stays nimages (or the frame count of path_file); the
     *     targets s*_i are equidistant in s between s(C_lo) and s(C_hi), the start geometry is the
     *     point of the CV polyline at s*_i (bracketing CV frames, linear in s between them), and
     *     setupUmbrella then uses s0_i = PathCV.s(start geometry) (not s*_i).
     *  6. With m > 0 every window is integrated (fixed_endpoints is switched off): the former
     *     endpoints (minima) become ordinary windows whose mean force is integrated.
     * histogramTrajectory uses the same CV (steps 1-4); its interior range is [s(C_lo), s(C_hi)]. */
    bool buildCvPath(bool place_windows);                   // Claude Generated 2026
    void checkExtensionEnergies(const std::vector<Geometry>& base,
                                const std::vector<Geometry>& extended) const; // Claude Generated 2026
    /*! \brief Frames the path CV is built on: m_cv_path if set, else the image frames. */
    const std::vector<Geometry>& cvFrames() const { return m_cv_path.empty() ? m_image_geoms : m_cv_path; }
    int m_cv_frames = 0;                         // CV frames after re-sampling (0 = path as given)
    int m_cv_extend = 0;                         // extra CV frames per path end (0 = off); effective value
                                                 // derived from cv_extend_length in buildCvPath when that is > 0
    double m_cv_extend_length = 0.0;             // extension per end as fraction of the CV path length (0 = off)
    std::string m_cv_atoms_spec = "all";         // cv_atoms as given (Claude Generated 2026)
    std::vector<int> m_cv_atoms;                 // atoms of the path CV RMSD, 0-based (empty = all)
    std::vector<double> m_cv_w;                  // the same as 0/1 weights for string_path (empty = all)
    std::vector<double> m_cv_w_align;            // alignment weights: 1 selected, 1e-6 others (empty = all)
    std::vector<Geometry> m_cv_path;             // CV frames (empty = image frames, default)
    int m_cv_lo = 0, m_cv_hi = -1;               // window range on the CV frames (hi < 0: last frame)
    std::vector<double> m_cv_orig_s;             // s of original path frames 0, 1, N-2, N-1 on the CV

    /* --- Distance-space path CV and window diagnostics. Claude Generated 2026 ---
     * cv_metric distances: d_i is the RMS difference of atom-pair coordinates (PathCV::Metric).
     * The pairs are resolved once per run (resolveCvPairs, in buildCvPath): an explicit list in
     * Initialise, 'auto' on the CV frames after the RMSD re-sampling and before the extension.
     * The CV frames are then re-sampled at equal distance-space arclength (PathCV::resampleByPairs).
     * Diagnostics (both metrics, no effect on the dynamics): every kDiagInterval steps of the
     * production phase each window records z, the nearest CV frame and the pair coordinates
     * (CV pairs; rmsd metric: heavy-atom 'auto' pairs m_diag_pairs). writeWindowDiagnostics()
     * writes <base>.neb.windows.csv and warns where the neighbour-jump index J > 3. */
    static constexpr int kDiagInterval = 10;
    std::string m_cv_metric_name = "rmsd";       // cv_metric as given
    std::string m_cv_pairs_spec = "auto";        // cv_pairs as given
    double m_cv_pair_threshold = 0.5;            // span threshold of 'auto' pairs [A]
    PathCV::Metric m_cv_metric = PathCV::Metric::Rmsd;
    PathCV::PairTransform m_cv_transform = PathCV::PairTransform::Distance;
    std::vector<PathCV::AtomPair> m_cv_pairs;    // pairs of the distance metric (0-based)
    std::vector<PathCV::AtomPair> m_diag_pairs;  // pairs of the J diagnostic (= m_cv_pairs for distances)
    PathCV::PairTransform m_diag_transform = PathCV::PairTransform::Distance;
    bool m_cv_pairs_resolved = false;
    std::vector<std::vector<double>> m_diag_z;   // per window: z every kDiagInterval steps
    std::vector<double> m_diag_nn_sum, m_diag_nn2_sum; // per window: nearest-frame index sums
    std::vector<Vector> m_diag_x_sum;            // per window: sum of the diagnostic pair coordinates
    std::vector<long> m_diag_n;                  // per window: diagnostic samples
    double m_diag_x_step = 0.0;                  // max_k max_p |x_p,k+1 - x_p,k| over the CV frames

    bool resolveCvPairs(const std::vector<Geometry>& frames);          // Claude Generated 2026
    std::unique_ptr<PathCV> makePathCV(const std::vector<Geometry>& frames) const; // Claude Generated 2026
    double cvSpacing(const Geometry& A, const Geometry& B) const;      // plain spacing in the CV metric
    std::string cvPairsText() const;                                   // "x = d, 3 pairs: 1-4 1-9 4-9"
    void accumulateWindowDiagnostics(int i, const Geometry& R);         // Claude Generated 2026
    void writeWindowDiagnostics() const;                               // Claude Generated 2026


    void resetPmfAccumulators(int nat);      // Claude Generated 2026

    /* Decoupled scheme: conditional mean force. Claude Generated 2026.
     *
     * Master equation (den Otter & Briels): for any smooth field b with b.grad(s) = 1,
     *
     *     dG/ds = <b.grad E>_s  -  kT <div b>_s
     *
     * The mass-weighted choice b_a = grad_a s / (m_a Z), Z = sum_a |grad_a s|^2 / m_a,
     * satisfies the condition by construction and makes Z the effective inverse mass of
     * the coordinate, which is also what the TST prefactor needs.
     *
     * Crucially, a bias that depends only on s is CONSTANT on the surface s = const and
     * therefore leaves the conditional distribution there untouched. Conditional averages
     * are thus unbiased regardless of the restraint, no unbiasing is required, and samples
     * from every window may be pooled into one global set of bins. */
    std::vector<long> m_fe_count;            // samples per bin
    std::vector<double> m_fe_force;          // sum of b.grad E            [Eh]
    std::vector<double> m_fe_div;            // sum of div b               [dimensionless]
    std::vector<double> m_fe_wall;           // sum of b.grad W (z wall part of m_fe_force) [Eh]
    std::vector<double> m_fe_ssum;           // sum of s (refined bin centre)
    std::vector<double> m_fe_wsum;           // sum of the ensemble weight (Z^-1/2 when constrained)
    std::vector<double> m_win_s_sum;         // per window: sum of s
    std::vector<double> m_win_s2_sum;        // per window: sum of s^2
    std::vector<long> m_win_count;           // per window: samples
    std::vector<double> m_win_z_sum;         // per window: sum of z (transverse amplitude)
    std::vector<double> m_win_z_max;         // per window: largest z (Claude Generated 2026)
    std::vector<long> m_win_zwall_n;         // per window: samples with z > z_wall
    std::vector<double> m_win_wall_sum;      // per window: sum of the wall energy W [Eh]
    /* Weighted per-window sums for the constrained scheme. Claude Generated 2026.
     * With s held exactly, every window contributes to a single value of s, so the binned
     * profile of the restrained scheme degenerates to one occupied bin per window. The
     * profile is therefore built from the windows themselves, which is also the textbook
     * form of a blue-moon calculation: one mean force per constrained value of the
     * coordinate, integrated over the grid of those values. */
    std::vector<double> m_win_w, m_win_wf, m_win_wd;
    std::vector<double> m_win_wwall;         // Z^-1/2 weighted sum of b.grad W (Claude Generated 2026)
    std::mt19937 m_rng{ 20260801u };         // fixed seed: reproducible probes
    /* Per-image state for the parallel window loop. Claude Generated 2026.
     *
     * Everything the per-image block writes is kept per image, so no worker ever touches a
     * counter another worker owns and no critical section is needed. The values are merged
     * only when they are reported. The probe generators are seeded per image for the same
     * reason the MD generators are: one shared stream would both race and correlate the
     * replicas. */
    std::vector<std::mt19937> m_probe_rng;
    mutable std::vector<double> m_img_max_jump;      // largest |s(t)-s(t-1)| per image
    mutable std::vector<long> m_img_max_jump_step;   // where it happened
    mutable std::vector<long> m_img_capped;          // restraint-force cap hits per image

    /*! \brief Per-image accumulation buffer, merged serially after the parallel loop. */
    struct FeLocal {
        std::vector<long> count;
        std::vector<double> force, div, ssum, wsum;
        std::vector<double> wall;                    // b.grad W part of force (Claude Generated 2026)
        double win_s = 0.0, win_s2 = 0.0, win_z = 0.0;
        double win_z_max = -1e300, win_wall = 0.0;  // z wall diagnostics (Claude Generated 2026)
        long win_zwall_n = 0;
        double win_w = 0.0, win_wf = 0.0, win_wd = 0.0; // Z^-1/2 weighted sums
        double win_wwall = 0.0;                      // ... of b.grad W (Claude Generated 2026)
        long win_n = 0, rejected = 0, nonfinite = 0;
        std::vector<long> rejected_bin;              // Z-floor rejections per s bin (Claude Generated 2026)
    };
    double m_fe_Zfloor = 0.0;                // reject samples whose |grad s| has collapsed
    long m_fe_rejected = 0;                  // how often that happened
    std::vector<long> m_fe_rejected_bin;     // ... per s bin (bin from s, also when rejected; Claude Generated 2026)
    std::vector<long> m_win_rejected;        // ... per window (Claude Generated 2026)
    long m_fe_nonfinite = 0;                 // samples dropped because the window went unstable
    double m_path_spacing = 0.0;             // mean adjacent-frame RMSD of the reference path
    mutable long m_bias_capped = 0;          // restraint-force cap hits
    bool m_cv_trace = false;                 // per-step trace of s and grad s
    mutable std::vector<double> m_prev_s;    // previous s per image (jump detection)
    mutable std::ofstream m_cv_trace_file;
    mutable double m_max_s_jump = 0.0;       // largest |s(t) - s(t-1)| seen
    mutable int m_max_jump_step = -1, m_max_jump_image = -1;

    Geometry pathBiasField(const Geometry& R, double* Z_out = nullptr) const; // Claude Generated 2026
    double windowK(int i) const;             // restraint k actually applied to window i
    double thermalEnergy() const;            // kT in Eh at the thermostat setpoint
    double biasRampFactor() const;
    void killRunawayWindows(const std::vector<int>& live); // Claude Generated 2026
    void accumulateCalibration(const std::vector<int>& live);  // Claude Generated 2026
    void calibrateStiffness();                                 // Claude Generated 2026
    bool constrainWindow(int i);                               // Claude Generated 2026           // 0 -> 1 over m_ramp_steps (Claude Generated 2026)
    void accumulateMeanForce(int i, SimpleMD* img, FeLocal& out, std::mt19937& rng) const;
    void writeMeanForceProfile() const;                 // Claude Generated 2026

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
    double m_z_wall = 0.0;                   // flat-bottom wall on z [A^2], 0 = off (Claude Generated 2026)
    double m_z_wall_k = 0.0;                 // its stiffness [Eh/A^4] (0 in the config = automatic)
    bool m_window_trajectory = false;        // per-image .trj.xyz (Claude Generated 2026)
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
    std::string m_scheme = "band";           // "band" (spring-coupled) or "independent"
    bool m_independent = false;              // m_scheme == "independent"
    int m_fe_bins = 100;                     // bins over s for the conditional mean force
    int m_div_probes = 1;                    // Hutchinson probes per step for div b
    double m_div_epsilon = 1e-3;             // FD displacement for div b [Angstrom]
    int m_ramp_steps = 0;                    // adiabatic switch-on of the restraints
    int m_window_threads = 1;                // windows advanced concurrently (independent scheme)
    double m_max_temperature = 0.0;          // kill threshold for a runaway window [K]
    bool m_adaptive_k = false;               // per-window stiffness from the measured slope
    double m_adaptive_k_max = 20.0;          // cap, as a multiple of umbrella_k
    std::vector<double> m_umbrella_kvec;     // per-window stiffness (all = umbrella_k unless adaptive)
    bool m_constrain = false;                // hard constraint on s instead of a restraint
    double m_constrain_tol = 1e-8;
    long m_shake_fail = 0;                   // SHAKE iterations that did not converge
    std::vector<double> m_calib_dg_sum;      // calibration: sum of b.grad E per window
    std::vector<long> m_calib_n;
    std::vector<int> m_window_died;          // step at which a window was declared dead (-1 = alive)
    std::vector<double> m_window_died_T;     // temperature [K] at which it was declared dead (Claude Generated 2026)
    /* Verbosity owned by the band. Claude Generated 2026.
     *
     * CurcumaLogger's level is a shared static, and the RAII save/restore of CurcumaMethod
     * cannot span concurrent workers: with several SimpleMD instances stepping at once, the
     * restores interleave and the level is left at the innermost (silent) value. The band
     * then finishes correctly but prints nothing - the free-energy profile was computed and
     * silently discarded. The band therefore remembers its own level and re-asserts it
     * before it reports. */
    int m_band_verbosity = 1;
    double m_bias_fmax = 0.5;                // cap on the largest per-atom restraint force
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
    PARAM(umbrella_kz, Double, 0.0, "Force constant of a tube restraint on the PERPENDICULAR path CV z. 0 = off. CAUTION - z is a squared distance (Angstrom^2), so this constant has units Eh/Angstrom^4 and innocent-looking values are enormous: the resulting force is kz*(z-z0)*|grad z|, and with |grad z| ~ 2d ~ 1.4 A a value of 2.0 gives ~2.8 Eh/Angstrom, roughly thirty times a typical bonded force. Measured: kz = 2.0 destroyed every window regardless of the time step, while the same runs without the tube were stable. Sane values are of order 0.05.", "FreeEnergy", {})
    PARAM(umbrella_z_tolerance, Double, 0.05, "Free tube radius in the units of z, i.e. Angstrom^2: the restraint acts only beyond it, so the image fluctuates freely inside. Must exceed the thermal transverse amplitude, otherwise real physics is frozen out instead of drift being suppressed - measured ~0.24 Angstrom^2 for helicene at 300 K, so the previously used 0.05 was far too tight and 1.0 is a reasonable tube.", "FreeEnergy", {})
    PARAM(umbrella_keep_springs, Bool, false, "Keep the neighbour springs in umbrella mode so an image still feels its neighbours. Helps it stay on the path, but the spring bias is NOT removed by WHAM - diagnostics, not production free energies.", "FreeEnergy", {})
    PARAM(umbrella_bins, Int, 100, "Number of histogram bins over s in [0,1] for the WHAM analysis.", "FreeEnergy", {})
    PARAM(reparametrize, Int, 0, "Redistribute the images to equal arclength every N band steps (string reparametrisation). 0 = off. Required for a meaningful free-energy profile; a value of 1-10 keeps the images as well-defined umbrella windows.", "FreeEnergy", {})
    PARAM(pmf, Bool, false, "REMOVED (Oct 2026): the neighbour-spring mean-force estimator had no well-defined coordinate (arclength of a fluctuating band) and read a force-contaminated gradient. Setting it is an error; use -nebmd.umbrella (with -nebmd.constrain for the conditional mean force).", "FreeEnergy", {})
    PARAM(fixman, Bool, true, "Include the metric (Fixman) term -kT/2 dlnZ/ds in dG/ds. s is a curvilinear coordinate, so the constrained-ensemble volume element is not constant. Expected magnitude is of order kT (~0.6 kcal/mol at 300 K) - it is a correctness term, not a fix for large errors.", "FreeEnergy", {})
    PARAM(pmf_equilibration, Int, 200, "Band steps discarded before the mean-force averaging starts (equilibration).", "FreeEnergy", {})
    // --- Decoupled scheme (reformulation, Claude Generated 2026) ---
    PARAM(scheme, String, "band", "Sampling scheme. 'band' = neighbour-spring coupled images (the original NEB-MD; kept for reproducibility). 'independent' = each image is restrained ONLY on the path CV s and feels no neighbours at all, so every window is a proper canonical ensemble and the windows are fully independent. 'independent' additionally switches on the conditional mean-force estimator. 'string' = finite-temperature string (stage 1 of the string -> path-CV chain): every image is restrained in full Cartesian space (rotation-invariant, Kabsch) to its frame of a mean path M, M is updated from the image averages and reparametrised to equal arclength every string_update steps; output <base>.neb.string.xyz, a k-dependent reference path for the path CV (not an exact MFEP).", "FreeEnergy", {})
    // --- Finite-temperature string (scheme 'string', Claude Generated 2026) ---
    PARAM(string_k, Double, 0.05, "Stiffness of the restraint that holds each image at its mean-path frame, V = k/2 sum_a m_a |R_a - M_aligned_a|^2, in Eh/(amu*Angstrom^2) (Eh/Angstrom^2 with mass_weighted_spring false). Holding requires roughly k*(d/2) > F_par,max (d = image spacing in the same metric); measured on ethane with 12 images 0.005 let the images slide. 0 = automatic from the start energies, k = 4 F_par,max/d with F_par,max = max_i |E(i+1)-E(i-1)|/(2d). The restraint frequency is sqrt(k) per 1.95 fs; values above ~0.2 at dt = 0.5 fs are warned about.", "String", {})
    PARAM(string_update, Int, 2000, "MD steps between two updates of the mean path M. Between updates M is fixed, so every image samples a canonical ensemble with a known static bias; the update replaces M_i by the mixture with the aligned average of image i.", "String", {})
    PARAM(string_alpha, Double, 0.5, "Mixing factor of an update: M_i <- (1-alpha) M_i + alpha <R_i aligned onto M_i>. 1 = replace by the average, smaller values damp the noise of short averaging windows.", "String", {})
    PARAM(string_smooth, Double, 0.0, "Laplacian smoothing per update, M_i <- (1-kappa) M_i + kappa/2 (M_i-1 + M_i+1) on interior frames before the reparametrisation. 0 = off.", "String", {})
    PARAM(string_tol, Double, 0.03, "Convergence threshold (Angstrom) on Delta_max = max_i mass-weighted best-fit RMSD between the mean path before and after an update (after reparametrisation). This is a noise estimate, not a measured value; check the history file.", "String", {})
    PARAM(string_patience, Int, 3, "Number of consecutive updates with Delta_max < string_tol required to declare convergence (at least 5 updates are always performed).", "String", {})
    PARAM(string_max_step, Double, 0.05, "Trust radius of one mean-path update in Angstrom (0 = off). The proposed update (after mixing, smoothing and equal-arclength reparametrisation) is measured frame by frame in the metric of the convergence test (mass-weighted best-fit RMSD between the old and the new frame); if its largest value Delta_max exceeds string_max_step, the step of EVERY frame is scaled by one common factor s = string_max_step/Delta_max (shape of the step kept, endpoints fixed), so no restraint centre jumps by more than string_max_step. Why: an abrupt jump of the restraint centres by Delta pumps roughly k/2 m Delta^2 per atom into the image at once; on helicene (GFN2, k = 0.05, Delta = 0.28 A after the first update of an ORCA MEP) this heated three images above the kill temperature. A scaled update never counts towards convergence; the factor s is logged per update and written as column 'scale' of the history file. Typical update displacements of a relaxed string are far below 0.05 A (ethane: <= 0.004 A), so the default only acts on large corrections of the start path.", "String", {})
    PARAM(string_abort_on_dead, Bool, true, "Abort the string stage when an image is declared dead (temperature above max_temperature). A dead image is frozen: it is no longer integrated and its mean-path frame is never updated again, so the remaining updates are small and would report a convergence that does not hold for that part of the path. true = stop immediately with an error (image numbers, step, temperature), write no <base>.neb.string.xyz, exit code != 0, and in the pipeline stage 2 is not started. false = continue as before, but such updates never count towards convergence and the end of the stage reports which frames are NOT relaxed.", "String", {})
    PARAM(string_prepare, Bool, true, "Prepare the start path of the string stage before the images are built: consecutive mass-weighted Kabsch alignment and equal-(mass-weighted)-arclength reparametrisation, so the images start exactly on the prepared frames. Without it the first updates only absorb the reparametrisation of an unequally spaced start path (helicene ORCA MEP: 0.28 A in 0.05 A trust-radius steps over 5 updates). The largest best-fit shift (Angstrom) is logged and written as column prep_shift; shifts <= 1e-6 A leave the path bit-identical.", "String", {})
    PARAM(string_ftol, Double, 5e-4, "Force convergence threshold of the string stage in Eh/Angstrom (RMS per atom of the perpendicular mean force F_a = k m_a perp_a, perp = part of <R_i> - M_i normal to the path). Converged when in string_patience consecutive updates max_i F_sig < string_ftol (F_sig = noise-corrected force) AND max_i SE < string_ftol (standard error from string_blocks block means: the force is resolved) AND the update was not scaled AND no image is dead AND the perpendicular step < string_tol. If SE >= string_ftol the criterion is noise-limited: the stage ends at string_max_time with a warning that names the factor (SE/ftol)^2 by which string_update would have to grow; this is not an error. 0 = former criterion (Delta_max < string_tol), kept only for reproduction - it reports convergence for any stiff restraint because Delta_max ~ alpha F/(k m).", "String", {})
    PARAM(string_blocks, Int, 5, "Number of blocks per mean-path update for the standard error of the mean force (block index = step within the update * blocks / string_update). At least 2; with fewer populated blocks SE is undefined (n/a) and the criterion counts as noise-limited.", "String", {})
    PARAM(string_gain, Double, 1.0, "Gain on the perpendicular part of a mean-path update: M_i <- M_i + alpha dev_i + alpha (G_i - 1) perp_i with G_i = min(string_gain, string_max_step/(alpha rmsd(perp_i))) (>= 1), applied only where the force is significant (z_i >= string_gain_z) AND above the convergence level (F_sig >= string_ftol > 0); elsewhere the plain mixing (1-alpha) M + alpha <R>. The fixed point is unchanged (perp = 0 there). 1 = off. Model (unit test): gain 20 reaches the valley in 0.33-0.52 of the plain updates; ethane: no benefit. WARNING: UNSTABLE on helicene (GFN2, gain 20): frames next to the endpoints rose to 73-83 kcal/mol within 4 updates - the gain assumes perpendicular curvatures far below k*m, which fails for bond stretches. Keep 1 on molecules.", "String", {})
    PARAM(string_gain_z, Double, 5.0, "Significance threshold z = F_sig/SE for the string gain. With 5 blocks z is t-distributed with 4 dof: z >= 3 alone triggered in 3.7 % of zero-signal frame-updates in the unit test; with z >= 5 and F_sig >= string_ftol 0/1100.", "String", {})
    PARAM(string_write_updates, Bool, false, "Append the mean path after every update (consecutively aligned, unweighted Kabsch, no energies) to <base>.neb.string_updates.xyz, comment 'string update=U image=i'. Diagnostic for the relaxation history (e.g. distance to a reference path per update).", "String", {})
    PARAM(string_max_time, Double, 20000.0, "MD length of the string stage in fs. Overrides simplemd.max_time for the images in scheme 'string'; the run ends earlier when the string has converged.", "String", {})
    PARAM(cv_histogram, String, "", "Analysis only, no MD: evaluate the path CV s,z of every frame of this multi-XYZ trajectory (e.g. an unbiased -md run) on the reference path (path_file, else the interpolated path; optimize/idpp/repair_overlaps are skipped so the path is used as written), histogram s into fe_bins bins and write G_ref(s) = -kT ln P(s) with a 5-block standard error to <base>.neb.ref.csv. kT from simplemd.temperature (default 300 K). Frames with z > cv_hist_zmax are discarded and counted.", "FreeEnergy", {})
    PARAM(cv_hist_zmax, Double, 0.0, "z filter for cv_histogram in Angstrom^2: frames with z above this value are discarded. 0 = automatic: the filter acts on the excess z - z_path(s), z_path linearly interpolated between the z values of the reference frames (z is offset below zero ON the path and the offset depends on s). The log10(excess) histogram is searched for a second, separated population (frames far from the path, e.g. another basin with permuted atom labels); if one is found the threshold is placed in the density valley between the on-path mode and that population, otherwise only outliers with excess > 10 x median excess are discarded (no filter when the median excess <= 0). Assumes the on-path frames are the majority. Threshold, rejected fraction and rejected fraction of the interior frames are reported in the console and the ref.csv header.", "FreeEnergy", {})
    PARAM(cv_frames, Int, 0, "Number of frames of the path collective variable, decoupled from the number of windows. 0 = the CV is built on the path frames as given (one window per frame, the former behaviour). > 0: the path (path_file, the stage-1 string, or the interpolated/NEB path) is consecutively Kabsch-aligned and re-sampled at equal arclength onto cv_frames frames (linear interpolation; more frames than the path interpolates, fewer decimates; endpoints kept). lambda of the CV follows from the new spacing, so the CV no longer changes with nimages. The window count stays nimages (or the frame count of path_file): the window setpoints are equidistant in s between the first and the last window-capable CV frame, the start geometry of a window is the point of the CV polyline at that s, and s0 is then re-evaluated on that geometry. Applies to scheme 'independent' and to cv_histogram (same construction there); ignored by the other schemes.", "PathCV", {})
    PARAM(cv_extend, Int, 0, "Number of additional CV frames beyond EACH end of the path (0 = off). The end segment is continued straight: F_-j = R_0 + j (R_0 - R_1'), F_N-1+j = R_N-1 + j (R_N-1 - R_N-2') for j = 1..cv_extend, with R_1' (R_N-2') Kabsch-aligned onto its neighbour, i.e. every extension step has the length of the last path segment. Motivation: with fixed endpoints the windows only span s(frame 1)..s(frame N-2), and the endpoints are plateaus of s (grad s -> 0), so the basin bottoms are never sampled. With cv_extend > 0 only the outermost extension frame per side is kept as plateau; windows may lie on all other CV frames, and the former endpoints (minima) become ordinary, integrated windows (fixed_endpoints is switched off). Without cv_frames there is one window per CV frame except the outermost per side. The single-point energy of every extension frame is compared with the endpoint; a lower energy means the path does not end in a basin there and is warned about (no abort). The basin free energies of the profile, G_A = -kT ln int_A exp(-G/kT) ds, are integrals over s, so their absolute value depends on how s is parametrised (cv_frames, cv_extend, lambda): compare basin barriers only between runs with the same CV definition. Applies to scheme 'independent' and to cv_histogram; ignored by the other schemes.", "PathCV", {})
    PARAM(cv_atoms, String, "all", "Atoms the path CV s,z is defined on (stage 2 / scheme 'independent', z wall, cv_trace, cv_histogram). 'all' (default) = every atom, the original CV. 'heavy' = all atoms with Z > 1. Otherwise a 1-based list with ranges, e.g. '1-6,9'. RMSD d_i, centroid and Kabsch rotation run over the selected atoms only (mean square over them); the CV gradient of the other atoms is exactly zero, but all atoms stay dynamic. The CV frames (cv_frames re-sampling, cv_extend/cv_extend_length) are aligned and spaced in the same selected-atom RMSD, and lambda follows from those spacings. At least 3 atoms are required (fewer cannot fix the rotation about their axis). Stage 1 (the string with its Kabsch restraint) is not affected. Why: with all atoms counted equally, thermal motion of light peripheral atoms (methyl H) inflates z, and the displacement of one heavy atom is diluted by the average over all atoms.", "PathCV", {})
    PARAM(cv_metric, String, "rmsd", "Metric of the path CV distances d_i (Claude Generated 2026). 'rmsd' (default) = best-fit RMSD to each CV frame (over cv_atoms), the original CV. 'distances' = root-mean-square difference of selected atom-pair coordinates x_p (cv_pairs, cv_pair_transform): d_i^2 = 1/P sum_p (x_p(R) - x_p(C_i))^2. No Kabsch fit; rotation-, translation- and mirror-invariant; atoms in no pair (e.g. methyl H) do not affect s and z. With 'distances' the CV frames are re-sampled (cv_frames) at equal distance-space arclength and lambda = 1/<d_k,k+1^2> in that metric; z is then in A^2 (d) or A^-2 (inv). Use it when the reactive atoms carry little of the RMSD path metric (rule of thumb: below ~90 %). Enantiomerisation paths are invisible in distance space.", "PathCV", {})
    PARAM(cv_pairs, String, "auto", "Atom pairs of the distance metric (cv_metric distances; Claude Generated 2026). 'auto' = every pair of the cv_atoms selection (cv_atoms 'all': every heavy atom, Z > 1) whose distance spans more than cv_pair_threshold Angstrom (max - min) over the CV frames (after the cv_frames re-sampling by RMSD arclength, before the extension); fewer than 2 such pairs is an error. Otherwise a 1-based list of pairs i-j, e.g. '1-4,1-9,4-9' ('1-4' is the PAIR (1,4), not atoms 1..4); i != j, no duplicates, at least 2 pairs; an explicit list overrides cv_atoms. The selected pairs and their spans are logged.", "PathCV", {})
    PARAM(cv_pair_threshold, Double, 0.5, "Minimum span (Angstrom, max - min of the distance over the CV frames) of a pair selected by cv_pairs 'auto' (Claude Generated 2026). The same rule on the heavy atoms gives the pairs of the window diagnostics (<base>.neb.windows.csv, neighbour-jump index J) in the rmsd metric.", "PathCV", {})
    PARAM(cv_pair_transform, String, "d", "Pair coordinate of the distance metric (Claude Generated 2026): 'd' = x = |R_a - R_b| (z in A^2), 'inv' = x = 1/|R_a - R_b| (z in A^-2; short distances weigh more, long ones fade out).", "PathCV", {})
    PARAM(cv_extend_length, Double, 0.0, "Extension of the path CV beyond EACH end given as a LENGTH: the fraction of the total arclength L of the CV path (after the cv_frames re-sampling) that the windows reach beyond the original endpoint. 0 = off (then cv_extend, in frames, applies unchanged). > 0 takes precedence over cv_extend (a warning is printed when both are set). Why: cv_extend counts FRAMES, and a frame covers a length that depends on cv_frames (more CV frames = shorter spacing), so a fixed cv_extend shrinks the sampled basin when cv_frames grows. Rule: with F CV frames and mean spacing d = L/(F-1) (RMSD arclength of the consecutively aligned path, the metric of the extension and of the path CV), the number of extension frames per end is m = ceil(cv_extend_length * L / d) + 1 = ceil(cv_extend_length * (F-1)) + 1; the +1 is the outermost extension frame, which remains the plateau of s and gets no window. The windows then cover at least cv_extend_length * L beyond each endpoint. The resulting m is logged ('cv_extend: <fraction> of the path length = m frames') and used exactly like cv_extend = m (straight continuation of the end segment, energy check of the extension frames, fixed_endpoints switched off). Without cv_frames the end segment may differ from the mean spacing d, so the covered length is only approximately cv_extend_length * L. Typical: 0.05-0.1. Applies to scheme 'independent' and to cv_histogram (same CV construction); ignored by the other schemes.", "PathCV", {})
    PARAM(pipeline, Bool, false, "Two-stage chain in one call: stage 1 = scheme 'string' (length string_max_time) produces <base>.neb.string.xyz; stage 2 reads it as path_file and runs scheme 'independent' with umbrella + constrain (conditional mean force along the frozen path CV, length simplemd.max_time). Both stages write into the same BMT directory; stage-2 files carry the basename <base>.stage2. In stage 2 the path is used as written: optimize, idpp and repair_overlaps are switched off there. -nebmd.constrain false keeps a restraint instead of the constraint in stage 2.", "Pipeline", {})
    PARAM(fe_bins, Int, 100, "Number of bins over s in [0,1] for the conditional mean-force profile. A bias that depends only on s is constant on the surface s=const, so conditional averages are unbiased and samples from ALL windows may be pooled into the same bins.", "FreeEnergy", {})
    PARAM(div_probes, Int, 1, "Random probes per step for the divergence term -kT<div b> (Hutchinson estimator). The term needs second derivatives of s, which are obtained here from finite differences of the analytic gradient - no extra energy evaluations. One probe per step suffices because the estimate is averaged over the whole trajectory.", "FreeEnergy", {})
    PARAM(div_epsilon, Double, 1e-3, "Displacement (Angstrom) for the finite-difference probe of div b. Must be small against the curvature scale of the path but large enough to stay above round-off.", "FreeEnergy", {})
    PARAM(ramp_steps, Int, 0, "Adiabatic switch-on: band steps over which the restraint constants (on s and on z) are raised smoothly from zero to their final values. 0 = use pmf_equilibration, so the switch-on finishes exactly when production sampling begins. Each image starts ON the path and AT its setpoint, where no restraint force is needed; applying full stiffness immediately delivers a shock that can push the image off the path, and once it is off, s is ill-defined and the restraint pushes in a meaningless direction with growing magnitude.", "FreeEnergy", {})
    PARAM(constrain, Bool, false, "Hold s by a genuine holonomic constraint (SHAKE on the positions, RATTLE on the velocities) instead of by a harmonic restraint. This is the only scheme in which holding is exact by construction: there is no stiffness that could be too small to hold the window or too large for the molecule to survive. Measured motivation: on a barrier flank the two requirements were k > 3 and k < 1 respectively, so no restraint exists. The price is the evaluation, which needs the Z^(-1/2) weighting of the constrained ensemble.", "FreeEnergy", {})
    PARAM(constrain_tolerance, Double, 1e-8, "Convergence threshold of the SHAKE iteration on |s - s0|.", "FreeEnergy", {})
    PARAM(z_wall, Double, 0.0, "Flat-bottom wall ('tube') on the distance-from-path z in the stage-2 windows (same units as z: Angstrom^2 for cv_metric rmsd and for distances with cv_pair_transform d, Angstrom^-2 for inv): W = 0 for z <= z_wall, W = z_wall_k/2 (z - z_wall)^2 beyond. 0 = off. Applied inside the integrator as part of the potential, constrained and restraint mode alike; its gradient enters the mean force. The resulting PMF is that of the z-restricted system and depends on z_wall: compare only runs with equal z_wall. Choose it above the on-path z of the windows (e.g. 3x the largest mean z of windows that stay on the path).", "FreeEnergy", {})
    PARAM(z_wall_k, Double, 0.0, "Stiffness of the z wall (Eh/Angstrom^4). 0 = automatic: W = 3 kT at a penetration of 0.2 z_wall, i.e. k = 150 kT / z_wall^2 (300 K, z_wall 0.5: 0.57).", "FreeEnergy", {})
    PARAM(window_trajectory, Bool, false, "Write a trajectory per window/image (<base>.image_NN.trj.xyz) every simplemd dump_frequency steps. Off by default: the files grow with images x steps.", "Output", {})
    PARAM(adaptive_k, Bool, false, "Per-window restraint stiffness, calibrated from the locally measured slope of the profile. A single global k cannot hold every window: the requirement is k > (dG/ds)/(ds/2), and dG/ds is near zero in the basins but large on the barrier flanks - measured on helicene, the two flank windows slid by 0.9 and 1.4 window spacings while the two basin windows held at the same k. The calibration is legitimate because the conditional mean force is unbiased regardless of the restraint, so it can be read off during the equilibration phase and used for production.", "FreeEnergy", {})
    PARAM(adaptive_k_max, Double, 20.0, "Upper bound for the calibrated stiffness, as a multiple of umbrella_k. Prevents the calibration from prescribing a value at which the window no longer survives; if the bound binds, the profile is steeper than the sampling can hold and that is reported.", "FreeEnergy", {})
    PARAM(max_temperature, Double, 0.0, "A window whose instantaneous temperature exceeds this value (K) is declared dead: it is removed from the integration and from all averages, and reported by name. 0 = ten times the thermostat setpoint. Without this a blown-up window keeps contributing - the failure is invisible in the averages and only the final energy betrays it, which is how several runs in this work produced plausible-looking profiles from damaged windows.", "FreeEnergy", {})
    PARAM(window_threads, Int, 1, "Number of windows advanced concurrently. Meaningful with scheme=independent and scheme=string (between two mean-path updates a string image is coupled only to its frozen frame M_i, so the images are independent there): no tangents, no springs, no neighbour positions. Parallelised is exclusively the MD step (energy and gradient), which is where the time goes; bias evaluation and accumulation stay serial, so no shared accumulator is touched from a worker. Because the windows are genuinely independent, the parallel result must be bit-identical to the serial one - any difference proves hidden shared state.", "FreeEnergy", {})
    PARAM(cv_trace, Bool, false, "Diagnostic: write one line per image and step with s, its change since the previous step and |grad s|, so that a temperature excursion can be correlated with a discontinuity of the collective variable. Best-fit RMSD can switch alignment branch on a symmetric molecule, which leaves s continuous but its GRADIENT jumping - a restraint built on it then delivers an impulsive force. Produces a large file; for diagnosis only.", "FreeEnergy", {})
    PARAM(bias_fmax, Double, 0.5, "Upper bound (Eh/Angstrom) for the largest per-atom component of the restraint force. The whole force vector is scaled uniformly when the bound is hit, so its direction - and thereby its meaning as a restraint on s - is preserved. Note that a SUSTAINED force of this size is far above typical bonded forces (~0.1); frequent capping is a warning, not a fix.", "FreeEnergy", {})
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