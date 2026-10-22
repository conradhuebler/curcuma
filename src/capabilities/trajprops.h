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
 * Quantities so far:
 *  - `gyration`: R_g^2 = sum m_i |r_i - R|^2 / sum m_i with R the centre of mass
 *    (unit weights with mass_weighted=false), of the selected atoms.
 *  - `hbonds`: geometric hydrogen bonds D-H...A. A hydrogen belongs to the nearest
 *    donor atom (element in hb_elements) within hb_bond_length; the bond is counted
 *    when the D...A distance is at most hb_distance and the angle D-H...A at least
 *    hb_angle. Between two groups (`selection`, `selection_b`) both directions are
 *    counted separately; with one group, bonds inside it.
 * The coordinates are used as stored: a chain that was wrapped around a periodic
 * boundary gives a wrong R_g, no hydrogen bond is found across a boundary, and an
 * XYZ trajectory carries no cell to unwrap with.
 *
 * The statistics report the plain standard deviation over frames and a block
 * standard error of the mean. Consecutive frames of an MD run are correlated, so
 * the standard deviation divided by sqrt(frames) understates the uncertainty of
 * the mean; the block error (consecutive blocks, spread of the block means) does
 * not assume independent frames but needs blocks longer than the correlation time.
 */

#pragma once

#include "curcumamethod.h"

#include "src/core/energy_calculators/ff_methods/gfnff_parameters.h"
#include "src/core/molecule.h"
#include "src/core/parameter_macros.h"

#include <string>
#include <vector>

namespace curcuma {

/// R_g of @p atoms (all atoms when empty). Mass weighted about the centre of mass, or geometric.
double gyrationRadius(const Molecule& molecule, const std::vector<int>& atoms, bool massWeighted);

/// Geometric hydrogen-bond criteria. Defaults are the common choice (D...A <= 3.5 A, D-H...A >= 150 deg).
struct HBondCriteria {
    double distance = 3.5;       ///< Angstrom, donor to acceptor
    double angle = 150.0;        ///< degrees, D-H...A at the hydrogen, at least
    double bondLength = 1.3;     ///< Angstrom, a hydrogen within this of a donor atom belongs to it
    std::vector<int> elements { 7, 8 };   ///< atomic numbers of donor and acceptor atoms (N, O)
};

struct HBondCount {
    int aDonates = 0;            ///< donor in group A, acceptor in group B (all bonds when B is empty)
    int bDonates = 0;            ///< donor in group B, acceptor in group A
    int total() const { return aDonates + bDonates; }
};

/// Hydrogen bonds between @p groupA and @p groupB, or inside @p groupA when @p groupB is empty.
/// An empty @p groupA means all atoms. Each (hydrogen, acceptor) pair is counted once.
HBondCount countHydrogenBonds(const Molecule& molecule, const std::vector<int>& groupA,
    const std::vector<int>& groupB, const HBondCriteria& criteria);

/// Energy criterion: a hydrogen-bond triple of GFN-FF counts when its energy is below a limit.
struct HBondEnergyOptions {
    double limitKJ = 0.0;   ///< kJ/mol; triples with energy < limit count (0: every attractive triple)
    int charge = 0;         ///< total charge of the system (GFN-FF needs it for the EEQ charges)
    int threads = 1;
};

struct HBondEnergyCount {
    bool ok = false;
    std::string error;
    HBondCount count;                     ///< triples below the limit, per direction (as countHydrogenBonds)
    double energySumKJ = 0.0;             ///< sum of the energies of all triples in the groups, kJ/mol
    std::vector<double> energiesKJ;       ///< energy of every triple in the groups that is not zero, kJ/mol
};

/// The counting step of the energy criterion alone, for triples that are already known (a trajectory frame
/// evaluated here, or a running MD). @p atomCount is the number of atoms the indices refer to; groups as in
/// countHydrogenBonds(); a triple counts when its energy is below @p limitKJ (kJ/mol).
HBondEnergyCount classifyHBondEnergies(const std::vector<GFNFFHBondEnergy>& triples, int atomCount,
    const std::vector<int>& groupA, const std::vector<int>& groupB, double limitKJ);

/// Hydrogen bonds by the GFN-FF triple energies (needs one GFN-FF evaluation of the frame, which sets up
/// topology and charges from scratch: far more expensive than the geometric count). Groups as in
/// countHydrogenBonds(). The triple list of GFN-FF is a candidate list; only triples whose energy is below
/// the limit count.
HBondEnergyCount countHydrogenBondsByEnergy(const Molecule& molecule, const std::vector<int>& groupA,
    const std::vector<int>& groupB, const HBondEnergyOptions& options);

/// Atoms of @p selection resolved against @p molecule. A leading '!' means every atom NOT in the
/// rest of the selection ("!F1": all but the first fragment). Empty selection: all atoms.
std::vector<int> resolveSelection(Molecule molecule, const std::string& selection);

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
    "tier=primary; enum=gyration|hbonds")
PARAM(selection, Selection, "", "Atoms the quantity is computed for (selection grammar, resolved once on the "
                                "first frame used; a leading ! means all atoms except the rest, e.g. !F1); "
                                "empty means all atoms.", "Basic", {}, "tier=primary")
PARAM(selection_b, Selection, "", "hbonds: second group. Bonds between the two groups are counted, per "
                                  "direction; empty counts bonds inside `selection`.", "Basic", {})
PARAM(mass_weighted, Bool, true, "Weight by atomic mass (false: geometric).", "Basic", {})
PARAM(first_frame, Int, 0, "First frame (0-based).", "Frames", {}, "min=0")
PARAM(last_frame, Int, -1, "Last frame; -1 runs to the end.", "Frames", {})
PARAM(stride, Int, 1, "Use every N-th frame of the range.", "Frames", {}, "min=1")
PARAM(last_fraction, Double, 0.0, "If in (0,1]: use only the last fraction of the trajectory (0.3333 for the last "
                                  "10 ns of 30 ns); overrides first_frame and last_frame.", "Frames", {}, "min=0")
PARAM(hb_distance, Double, 3.5, "hbonds: largest donor...acceptor distance (Angstrom).", "Hydrogen bonds", {}, "min=0")
PARAM(hb_angle, Double, 150.0, "hbonds: smallest D-H...A angle (degrees).", "Hydrogen bonds", {}, "min=0")
PARAM(hb_bond_length, Double, 1.3, "hbonds: a hydrogen within this distance of a donor atom belongs to it (Angstrom).",
    "Hydrogen bonds", {}, "min=0")
PARAM(hb_elements, String, "N,O", "hbonds: element symbols of donor and acceptor atoms, comma separated.",
    "Hydrogen bonds", {})
PARAM(hb_criterion, String, "geometric", "hbonds: geometric (distance and angle) or energy (GFN-FF triple energy below "
    "hb_energy_limit; one GFN-FF evaluation per frame, use stride).", "Hydrogen bonds", {}, "enum=geometric|energy")
PARAM(hb_energy_limit, Double, 0.0, "hbonds, criterion energy: a triple counts when its energy is below this (kJ/mol, "
    "negative = attractive). 0 counts every attractive triple; calibrate with the histogram.", "Hydrogen bonds", {})
PARAM(hb_hist_bins, Int, 40, "hbonds, criterion energy: bins of the triple-energy histogram over [-hb_hist_max, 0].",
    "Hydrogen bonds", {}, "min=1")
PARAM(hb_hist_max, Double, 20.0, "hbonds, criterion energy: lower edge of the histogram (kJ/mol, positive number).",
    "Hydrogen bonds", {}, "min=0")
PARAM(charge, Int, 0, "hbonds, criterion energy: total charge of the system for GFN-FF.", "Hydrogen bonds", {})
PARAM(threads, Int, 1, "hbonds, criterion energy: threads of the GFN-FF evaluation.", "Hydrogen bonds", {}, "min=1")
PARAM(blocks, Int, 5, "Number of blocks for the block standard error.", "Statistics", {}, "min=1")
PARAM(output, String, "", "Per-frame CSV; empty gives <basename>_<quantity>.csv.", "Output", {})
END_PARAMETER_DEFINITION
} // namespace
