/*
 * <Automatic z threshold for the path-CV reference histogram (population separation).>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested).
 * Header-only so that NebMD::histogramTrajectory and the unit test share one implementation.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace curcuma {
namespace z_filter {

/*! \brief Result of the automatic z threshold choice (Claude Generated 2026). */
struct ZThreshold {
    double threshold = std::numeric_limits<double>::infinity(); ///< keep frames with excess <= threshold
    bool bimodal = false;      ///< true: a separated far population was found, threshold sits in the valley
    double median = 0.0;       ///< median of the excess over all frames
    double mode = 0.0;         ///< excess at the on-path mode (centre of the densest low-excess log bin)
    double far_peak = 0.0;     ///< excess at the far-population peak (bimodal only)
    double score = 0.0;        ///< significance of the far peak over the valley (Poisson z-like score)
    std::string rule;          ///< ASCII description of the rule that set the threshold
};

/*! \brief Choose a z threshold that separates on-path frames from a distant population.
 *
 * Input is the excess e = z - z_path(s) of every frame (A^2). Frames of a free MD that sit
 * near the reference path have a chi^2-like excess distribution; frames in a different
 * basin (e.g. a permuted methyl group of ethane, i.e. swapped atom labels) sit far away
 * and form a second, separated population. A threshold relative to the mean (the former
 * 3 x mean rule) cuts INTO the on-path tail when the far population inflates the mean,
 * and lets far frames through otherwise. This rule looks for the gap instead:
 *
 *  1. floor = 1e-3 * median(e); L = log10(max(e, floor)); histogram L in 0.1-decade bins,
 *     smoothed by a 3-bin running SUM (keeps Poisson statistics: var(h) ~ h).
 *  2. On-path mode: the densest smoothed bin at or below the median bin. This assumes the
 *     on-path frames are the majority (> 50 %); a narrow far population can be denser per
 *     log bin than the broad on-path one, so the global maximum is NOT used.
 *  3. Scanning right of the mode with the running minimum m of h, the far peak is the bin j
 *     with the largest score (h_j - m) / sqrt(h_j + m + 1). Bimodal if score >= 5 and
 *     m <= 0.5 h_j. The threshold is the centre of the longest run of bins at the valley
 *     minimum between mode and far peak (the middle of an empty gap, or the density minimum).
 *  4. Otherwise (unimodal): only true outliers are discarded, excess > 10 x median.
 *     A median <= 0, or fewer than 20 frames, disables the filter.
 *
 * \param excess  z - z_path(s) per frame, A^2
 * \return        threshold (frames with excess > threshold are discarded) and diagnostics
 */
inline ZThreshold chooseThreshold(const std::vector<double>& excess)
{
    ZThreshold r;
    const size_t n = excess.size();
    if (n < 20) {
        r.rule = "auto: fewer than 20 frames, no filter";
        return r;
    }
    std::vector<double> sorted(excess);
    std::sort(sorted.begin(), sorted.end());
    r.median = (n % 2) ? sorted[n / 2] : 0.5 * (sorted[n / 2 - 1] + sorted[n / 2]);
    if (!(r.median > 0.0)) {
        r.rule = "auto: median excess <= 0, no filter";
        return r;
    }
    const double floor_v = 1e-3 * r.median;
    const double lo = std::log10(floor_v);
    const double hi = std::log10(std::max(sorted.back(), floor_v));
    double w = 0.1; // bin width in decades
    int nb = static_cast<int>(std::ceil((hi - lo) / w)) + 1;
    if (nb > 400) { // pathological range: widen the bins instead of allocating more
        w = (hi - lo) / 399.0;
        nb = 400;
    }
    auto binOf = [&](double e) {
        const double L = std::log10(std::max(e, floor_v));
        return std::min(nb - 1, std::max(0, static_cast<int>(std::floor((L - lo) / w))));
    };
    auto centre = [&](double b) { return std::pow(10.0, lo + (b + 0.5) * w); };
    std::vector<double> c(nb, 0.0), h(nb, 0.0);
    for (double e : excess)
        c[binOf(e)] += 1.0;
    for (int b = 0; b < nb; ++b)
        h[b] = c[b] + (b > 0 ? c[b - 1] : 0.0) + (b + 1 < nb ? c[b + 1] : 0.0);

    const int bm = binOf(r.median);
    int mode = 0;
    for (int b = 1; b <= bm; ++b)
        if (h[b] > h[mode])
            mode = b;
    r.mode = centre(mode);

    // far peak: largest rise above the running minimum right of the mode
    double mv = h[mode], best = -1.0;
    int jbest = -1;
    double mbest = 0.0;
    for (int j = mode + 1; j < nb; ++j) {
        mv = std::min(mv, h[j]);
        const double sc = (h[j] - mv) / std::sqrt(h[j] + mv + 1.0);
        if (sc > best) {
            best = sc;
            jbest = j;
            mbest = mv;
        }
    }
    r.score = std::max(best, 0.0);
    if (jbest > mode && best >= 5.0 && mbest <= 0.5 * h[jbest]) {
        // valley: longest contiguous run of bins at the minimum between mode and far peak
        int run_a = -1, run_len = 0, cur_a = -1, cur_len = 0;
        for (int b = mode + 1; b < jbest; ++b) {
            if (h[b] == mbest) {
                if (cur_len == 0)
                    cur_a = b;
                ++cur_len;
                if (cur_len > run_len) {
                    run_len = cur_len;
                    run_a = cur_a;
                }
            } else
                cur_len = 0;
        }
        if (run_a >= 0) {
            r.bimodal = true;
            r.far_peak = centre(jbest);
            r.threshold = std::pow(10.0, lo + (run_a + 0.5 * run_len) * w);
            r.rule = "auto: valley between on-path and far population";
            return r;
        }
    }
    r.threshold = 10.0 * r.median;
    r.rule = "auto: unimodal, outliers > 10 x median excess";
    return r;
}

} // namespace z_filter
} // namespace curcuma
