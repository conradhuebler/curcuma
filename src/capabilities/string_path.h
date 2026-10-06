/*
 * <Pure path/string helpers for the finite-temperature string (NEB-MD scheme 'string'):
 *  weighted Kabsch, weighted RMSD, equal-arclength reparametrisation, consecutive
 *  alignment and smoothing of a discretised path.>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested - not TESTED/APPROVED).
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
 * Why a separate module: the string method (E, Ren & Vanden-Eijnden, PRB 66, 052301 (2002);
 * finite-temperature string, Maragliano & Vanden-Eijnden, CPL 426, 168 (2006)) consists of
 * two parts - sampling (MD, lives in NebMD) and pure geometry on a list of frames (averaging,
 * reparametrisation, alignment). The geometric part has no state and no MD, so it is kept here
 * as free functions that can be unit-tested on analytic model paths (test_string_path).
 *
 * Conventions
 *  - A Geometry is (natoms x 3) in Angstrom (row = atom).
 *  - Weights w (one per atom) enter every metric as  |dR|_w^2 = sum_a w_a |dR_a|^2.
 *    An EMPTY weight vector means uniform weights 1 (plain Euclidean / Frobenius norm).
 *  - A rotation `rot` acts on column vectors; for row-stored atoms the aligned geometry is
 *        x' = (x - c_src) * rot^T + c_dst .
 */

#pragma once

#include <vector>

#include <Eigen/Dense>

#include "src/core/global.h"

namespace curcuma {
namespace string_path {

/*! \brief Rigid-body superposition src -> dst (Claude Generated 2026). */
struct Alignment {
    Eigen::Matrix3d rot = Eigen::Matrix3d::Identity(); ///< proper rotation (det +1)
    Eigen::RowVector3d c_src = Eigen::RowVector3d::Zero(); ///< weighted centroid of src
    Eigen::RowVector3d c_dst = Eigen::RowVector3d::Zero(); ///< weighted centroid of dst
};

/*! \brief Weighted Kabsch: the proper rotation + translation minimising
 *         sum_a w_a | rot (src_a - c_src) + c_dst - dst_a |^2   (Claude Generated 2026).
 *
 *  Centroids are w-weighted; the covariance is H = sum_a w_a (src_a-c_src)^T (dst_a-c_dst),
 *  rot = V D U^T from SVD(H) = U S V^T with D fixing det(rot) = +1 (no reflections).
 *  Template: PathCV::rmsdTo (unweighted). */
Alignment kabsch(const Geometry& src, const Geometry& dst, const std::vector<double>& w);

/*! \brief Apply an Alignment to every row of g (Claude Generated 2026). */
Geometry applyAlignment(const Geometry& g, const Alignment& A);

/*! \brief Weighted RMSD  sqrt( sum_a w_a |a_a - b_a|^2 / sum_a w_a ) (Claude Generated 2026).
 *  best_fit = true superimposes a onto b (weighted Kabsch) first. */
double weightedRmsd(const Geometry& a, const Geometry& b, const std::vector<double>& w, bool best_fit);

/*! \brief Re-sample a polyline at equal (weighted) arclength (Claude Generated 2026).
 *
 *  Segment length |P_{i+1} - P_i|_w (plain differences, NO alignment - the frames must share
 *  one coordinate frame), cumulative arclength s_i, targets s_i* = L i/(N-1), linear
 *  interpolation inside the segment containing s_i*. Endpoints are returned bit-identical.
 *  With empty w this reproduces the former NebMD::reparametrizeString() kernel exactly. */
std::vector<Geometry> reparametrizeEqualArclength(const std::vector<Geometry>& path,
                                                  const std::vector<double>& w);

/*! \brief Superimpose every frame onto its (already aligned) predecessor (Claude Generated 2026).
 *
 *  Frame 0 is returned unchanged. Afterwards the PLAIN difference of consecutive frames equals
 *  their best-fit difference in the metric of w - which is what consumers that use plain
 *  differences (PathCV::computeLambdaFromPath, NebMD path spacing) silently assume. */
std::vector<Geometry> alignConsecutive(const std::vector<Geometry>& path, const std::vector<double>& w);

/*! \brief max_i weightedRmsd(a_i, b_i) over two equally long frame lists (Claude Generated 2026). */
double maxFrameDisplacement(const std::vector<Geometry>& a, const std::vector<Geometry>& b,
                            const std::vector<double>& w, bool best_fit);

/*! \brief One Laplacian smoothing pass  P_i <- (1-kappa) P_i + kappa/2 (P_{i-1} + P_{i+1})
 *         for the interior frames; endpoints unchanged; kappa <= 0 returns the input
 *         (Claude Generated 2026). */
std::vector<Geometry> smoothPath(const std::vector<Geometry>& path, double kappa);

/*! \brief Re-sample a polyline onto nout frames at equal (weighted) arclength (Claude Generated 2026).
 *
 *  Generalisation of reparametrizeEqualArclength to a different frame count: segment lengths
 *  |P_{i+1} - P_i|_w from PLAIN differences (the frames must share one coordinate frame, e.g.
 *  after alignConsecutive), targets s_k* = L k/(nout-1), linear interpolation inside the
 *  segment containing s_k*. nout > N interpolates, nout < N decimates. The first and last
 *  frame are returned bit-identical. For nout == N the arithmetic is that of
 *  reparametrizeEqualArclength (bit-identical result). nout < 2 or N < 2 returns the input. */
std::vector<Geometry> resamplePath(const std::vector<Geometry>& path, int nout,
                                   const std::vector<double>& w);

/*! \brief Extend a path by m frames beyond each end by linear extrapolation (Claude Generated 2026).
 *
 *  With R_1' = R_1 superimposed onto R_0 and R_{N-2}' = R_{N-2} superimposed onto R_{N-1}
 *  (weighted Kabsch, weights w):
 *
 *      F_{-j}    = R_0     + j (R_0     - R_1'    )      j = 1..m
 *      F_{N-1+j} = R_{N-1} + j (R_{N-1} - R_{N-2}')
 *
 *  i.e. the last segment is continued straight with its own (best-fit) length. Returned order:
 *  F_{-m}, ..., F_{-1}, R_0, ..., R_{N-1}, F_N, ..., F_{N-1+m}; the original frames are copied
 *  bit-identically, atom order is unchanged. Because R_1' is optimally superimposed onto R_0,
 *  the plain difference F_{-1} - R_0 satisfies the Kabsch stationarity condition, so if the input
 *  is consecutively aligned the output is as well (plain RMSD = best-fit RMSD of neighbours).
 *  No physical check is made here: an extrapolated frame can be unphysical (compressed bonds)
 *  if the end segment is long; callers should check its energy. m <= 0 or N < 2 returns input. */
std::vector<Geometry> extendPath(const std::vector<Geometry>& path, int m, const std::vector<double>& w);

/*! \brief Number of extension frames per path end for an extension given as a LENGTH
 *         (Claude Generated 2026).
 *
 *  For a path of nframes frames at (equal-arclength) spacing d = L/(nframes-1), an extension by
 *  length_fraction * L beyond each end needs
 *
 *      m = ceil(length_fraction * L / d) + 1 = ceil(length_fraction * (nframes - 1)) + 1
 *
 *  frames: extendPath continues the end segment with spacing d, and the outermost extension
 *  frame stays a plateau of s (no window), so the windows reach (m-1) d >= length_fraction * L
 *  beyond the original endpoint. The ceil is taken with a 1e-9 tolerance so that exact
 *  multiples of d (e.g. 0.1 * 10 = 1.0000000000000002) are not rounded up by round-off.
 *  length_fraction <= 0 or nframes < 2 returns 0 (extension off). Pure, unit-tested
 *  (test_string_path). */
int extensionFrames(double length_fraction, int nframes);

/*! \brief Result of limitStep (Claude Generated 2026). */
struct StepLimit {
    std::vector<Geometry> path; ///< path after the trust-radius limitation
    double scale = 1.0;         ///< common scaling factor s in (0, 1]; 1 = proposed path unchanged
    double delta_max = 0.0;     ///< max_i weightedRmsd(proposed_i, old_i, w, best_fit) BEFORE scaling
};

/*! \brief Trust radius for a mean-path update of the string (Claude Generated 2026).
 *
 *  Measures the proposed step frame by frame in the metric of the string convergence test,
 *  r_i = weightedRmsd(proposed_i, old_i, w, best_fit = true), delta_max = max_i r_i (endpoints
 *  excluded when fix_endpoints). If max_step <= 0 or delta_max <= max_step, the proposed path is
 *  returned UNCHANGED (bit-identical, scale = 1). Otherwise ONE common factor
 *
 *      s = max_step / delta_max
 *
 *  scales the step of every frame:  out_i = old_i + s (P_i - old_i),  P_i = proposed_i
 *  superimposed onto old_i (weighted Kabsch, weights w). The shape of the step is kept (each
 *  frame moves parallel to its best-fit step); only the rigid-body part of the step, which the
 *  rotation-invariant restraint does not see, is dropped. Because P_i is the optimal
 *  superposition, the cross-covariance of (old_i, P_i) is symmetric positive semidefinite and so is
 *  the one of (old_i, out_i): out_i is optimally superimposed on old_i as well, i.e. its best-fit
 *  displacement is exactly s r_i <= max_step. With fix_endpoints the first and last frame are
 *  copied bit-identically from old. Pure, unit-tested (test_string_path (h)). */
StepLimit limitStep(const std::vector<Geometry>& old_path, const std::vector<Geometry>& proposed,
                    double max_step, const std::vector<double>& w, bool fix_endpoints);

/*! \brief Result + diagnostics of one finite-temperature-string update (Claude Generated 2026).
 *
 *  Per-frame vectors have one entry per frame; frames that were not updated (no samples, or a
 *  fixed endpoint) carry 0 there and are excluded from the maxima. Forces in the units of
 *  k * w * length (NebMD: Eh/A with k in Eh/(amu A^2) and w = masses). */
struct MeanPathUpdate {
    std::vector<Geometry> path;     ///< new mean path (after gain, mixing, smoothing, reparametrisation, limitStep)
    std::vector<double> perp_rmsd;  ///< rmsd_w of the perpendicular part of <R_i> - M_i [length]
    std::vector<double> f_rms;      ///< RMS per atom of the perpendicular mean force k w_a perp_a
    std::vector<double> f_se;       ///< standard error of f_rms from the block means (+inf if < 2 blocks)
    std::vector<double> f_sig;      ///< sqrt(max(f_rms^2 - f_se^2, 0)), noise-corrected force
    std::vector<double> z;          ///< f_sig / f_se (0 if f_se is not finite or 0)
    std::vector<double> gain;       ///< applied gain G_i (1 = plain mixing)
    double perp_max = 0.0;          ///< max_i perp_rmsd (the former 'perp_max' column)
    double f_rms_max = 0.0, f_se_max = 0.0, f_sig_max = 0.0, z_max = 0.0, gain_max = 1.0;
    double delta_max = 0.0;         ///< max_i best-fit rmsd_w(M_new_i, M_old_i) BEFORE the trust radius
    double delta_rms = 0.0;         ///< RMS of the same over the interior frames
    double delta_perp_max = 0.0;    ///< max_i rmsd_w of the perpendicular part of the best-fit step (before the trust radius)
    double scale = 1.0;             ///< trust-radius factor of limitStep (1 = not limited)
    double arclength = 0.0;         ///< sum of plain rmsd_w between consecutive frames of the new path
    int n_updated = 0;              ///< frames that received a mean
    MeanPathUpdate() = default;
    MeanPathUpdate(const MeanPathUpdate&) = default;
    MeanPathUpdate(MeanPathUpdate&&) = default;
    MeanPathUpdate& operator=(const MeanPathUpdate&) = default;
    MeanPathUpdate& operator=(MeanPathUpdate&&) = default;
    ~MeanPathUpdate(); ///< out of line (avoids -Winline on the vector members)
};

/*! \brief One mean-path update of the finite-temperature string, with force/noise statistics
 *         and an optional perpendicular gain (Claude Generated 2026).
 *
 *  For every frame i with a mean (means[i].rows() > 0; endpoints skipped when fix_endpoints):
 *    dev_i  = mean_i - M_i                     (mean already superimposed onto M_i)
 *    tau_i  = (M_{i+1} - M_{i-1}) / |.|_w      (plain difference; the path must be aligned)
 *    perp_i = dev_i - (dev_i . tau_i)_w tau_i
 *    F_ia   = k w_a perp_ia,  f_rms_i = sqrt(sum_{a,c} F_iac^2 / N_at)
 *  Stiff harmonic restraint: <R> - M = F_mean/(k m), so F is the perpendicular mean force
 *  (k-independent). Its noise: the same F from each block mean, s2_ac = sample variance over the
 *  n_b blocks, f_se_i = sqrt(sum_{a,c} s2_ac / n_b / N_at); f_sig_i = sqrt(max(f_rms^2 - f_se^2, 0))
 *  (E[f_rms^2] = F_true^2 + f_se^2), z_i = f_sig_i / f_se_i.
 *  Gain (the fixed point is unchanged because only the perpendicular part is amplified, and it is
 *  zero at the fixed point): G_i = 1 if gain <= 1, ftol <= 0, z_i < gain_z or f_sig_i < ftol, else
 *  G_i = max(1, min(gain, max_step / (alpha rmsd_w(perp_i)))).
 *    G_i == 1:  M_i' = (1-alpha) M_i + alpha mean_i               (bit-identical to the former update)
 *    G_i  > 1:  M_i' = M_i + alpha dev_i + alpha (G_i - 1) perp_i
 *  Then smoothPath(smooth) (if > 0), reparametrizeEqualArclength(w), endpoints copied exactly
 *  (fix_endpoints), delta/delta_perp measured, limitStep(max_step) if delta_max > max_step > 0.
 *  With gain <= 1 the returned path is bit-identical to the former NebMD::updateStringPath. */
MeanPathUpdate updateMeanPath(const std::vector<Geometry>& old_path, const std::vector<Geometry>& means,
                              const std::vector<std::vector<Geometry>>& block_means, double k,
                              const std::vector<double>& w, double alpha, double gain, double gain_z,
                              double ftol, double max_step, double smooth, bool fix_endpoints);

} // namespace string_path
} // namespace curcuma
