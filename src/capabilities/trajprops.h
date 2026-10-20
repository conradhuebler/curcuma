/*
 * TrajectoryProperties — streaming properties of a long trajectory.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 - A solvated polymer run is tens of thousands of atoms and
 * thousands of frames. Measurement keeps every frame in memory and UnifiedAnalysis
 * keeps every selected frame, so neither is usable there. This capability reads one
 * frame at a time (FileIterator), resolves the atom selection once on the first
 * frame of the range, and keeps only the per-frame numbers.
 *
 * Quantity so far: `gyration`, R_g^2 = sum m_i |r_i - R|^2 / sum m_i with R the
 * centre of mass (unit weights with mass_weighted=false), of the selected atoms.
 * The coordinates are used as stored: a chain that was wrapped around a periodic
 * boundary gives a wrong R_g, and an XYZ trajectory carries no cell to unwrap with.
 *
 * The statistics report the plain standard deviation over frames and a block
 * standard error of the mean. Consecutive frames of an MD run are correlated, so
 * the standard deviation divided by sqrt(frames) understates the uncertainty of
 * the mean; the block error (consecutive blocks, spread of the block means) does
 * not assume independent frames but needs blocks longer than the correlation time.
 */

#pragma once

#include "curcumamethod.h"

#include "src/core/molecule.h"
#include "src/core/parameter_macros.h"

#include <string>
#include <vector>

namespace curcuma {

/// R_g of @p atoms (all atoms when empty). Mass weighted about the centre of mass, or geometric.
double gyrationRadius(const Molecule& molecule, const std::vector<int>& atoms, bool massWeighted);

/// Which of @p total frames are used: the last `lastFraction` of them when that is in (0, 1],
/// otherwise first..last (last < 0 means the end), then every `stride`-th. Empty when nothing qualifies.
std::vector<int> resolveFrames(int total, int first, int last, double lastFraction, int stride);

struct SeriesStatistics {
    int n = 0;
    double mean = 0.0;
    double stddev = 0.0;        ///< sample standard deviation over the frames
    double min = 0.0;
    double max = 0.0;
    double median = 0.0;
    int blocks = 0;             ///< blocks actually used
    double blockError = 0.0;    ///< standard error of the mean from the block means; 0 with fewer than 2 blocks
};

/// Statistics of @p values; the block error uses @p blocks consecutive blocks (fewer when there are fewer values).
SeriesStatistics seriesStatistics(const std::vector<double>& values, int blocks);

class TrajectoryProperties : public CurcumaMethod {
public:
    explicit TrajectoryProperties(const json& controller, bool silent = true);
    ~TrajectoryProperties() = default;

    /// Read the trajectory given to setFile(), write `<basename>_<quantity>.csv` and fill Results().
    void start() override;

    json Results() const override { return m_results; }
    std::string error() const { return m_error; }

private:
    nlohmann::json WriteRestartInformation() override { return json(); }
    bool LoadRestartInformation() override { return true; }
    StringList MethodName() const override { return { std::string("trajprops") }; }
    void ReadControlFile() override {}
    void LoadControlJson() override {}

    json m_results;
    std::string m_error;
};

} // namespace curcuma

namespace {
BEGIN_PARAMETER_DEFINITION(trajprops)
MODULE_INFO("Streaming properties of a long trajectory, one frame in memory at a time.",
    "Analysis", { "trajectory", "gyration" })

PARAM(quantity, String, "gyration", "What to compute per frame.", "Basic", {},
    "tier=primary; enum=gyration")
PARAM(selection, Selection, "", "Atoms the quantity is computed for (selection grammar, resolved once on the "
                                "first frame used); empty means all atoms.", "Basic", {}, "tier=primary")
PARAM(mass_weighted, Bool, true, "Weight by atomic mass (false: geometric).", "Basic", {})
PARAM(first_frame, Int, 0, "First frame (0-based).", "Frames", {}, "min=0")
PARAM(last_frame, Int, -1, "Last frame; -1 runs to the end.", "Frames", {})
PARAM(stride, Int, 1, "Use every N-th frame of the range.", "Frames", {}, "min=1")
PARAM(last_fraction, Double, 0.0, "If in (0,1]: use only the last fraction of the trajectory (0.3333 for the last "
                                  "10 ns of 30 ns); overrides first_frame and last_frame.", "Frames", {}, "min=0")
PARAM(blocks, Int, 5, "Number of blocks for the block standard error.", "Statistics", {}, "min=1")
PARAM(output, String, "", "Per-frame CSV; empty gives <basename>_<quantity>.csv.", "Output", {})
END_PARAMETER_DEFINITION
} // namespace
