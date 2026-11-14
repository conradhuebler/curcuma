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
 *  - `rdf`: radial distribution function g(r) between two atom groups with running
 *    coordination number n(r), optional element filters, minimum image convention in an
 *    orthorhombic cell given by the user (an XYZ file carries none).
 *  - `sasa`: solvent-accessible surface area of an atom group by the Shrake-Rupley
 *    test-point method with a probe sphere (default 1.4 A, a water molecule).
 * The coordinates are used as stored: a chain that was wrapped around a periodic
 * boundary gives a wrong R_g, no hydrogen bond is found across a boundary, and an
 * XYZ trajectory carries no cell to unwrap with.
 *
 * LITERATURE (TODO: written by the AI from memory, not yet checked against the papers by the
 * operator; also in docs/TRAJPROPS.md, the citation database and TODO.md):
 *  - block error: Flyvbjerg, Petersen, J. Chem. Phys. 91, 461 (1989)
 *  - geometric hydrogen bond: cf. Luzar, Chandler, Nature 379, 55 (1996); Kumar, Schmidt,
 *    Skinner, J. Chem. Phys. 126, 204107 (2007)
 *  - g(r), n(r): Allen, Tildesley, Computer Simulation of Liquids (2nd ed., 2017)
 *  - Shrake-Rupley: J. Mol. Biol. 79, 351 (1973); probe sphere: Lee, Richards, J. Mol. Biol.
 *    55, 379 (1971); radii: Bondi, J. Phys. Chem. 68, 441 (1964); test points on a sphere:
 *    Gonzalez, Math. Geosci. 42, 49 (2010)
 *  - the setup this serves: Min et al., Soft Matter (2015), doi 10.1039/c4sm02242d
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

#include <functional>
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

/// Element symbols "N,O" -> atomic numbers. False with @p error set for an unknown symbol.
bool parseElementList(const std::string& list, std::vector<int>& numbers, std::string* error);

/// Keep the atoms of @p atoms whose element is in @p elements (all atoms of the molecule when @p atoms is empty;
/// no filtering when @p elements is empty).
std::vector<int> filterByElement(const Molecule& molecule, const std::vector<int>& atoms, const std::vector<int>& elements);

// --- radial distribution function -------------------------------------------------------------------
struct RdfOptions {
    double rMax = 10.0;                 ///< Angstrom
    double binWidth = 0.05;             ///< Angstrom
    bool periodic = false;
    double cell[3] = { 0.0, 0.0, 0.0 }; ///< orthorhombic edges (Angstrom), minimum image when periodic
};

/// Add the ordered pairs (i in A, j in B, i != j) of one frame to @p histogram (bins of binWidth up to rMax).
/// Groups must be identical or disjoint (checked by the caller).
void accumulateRdf(const Molecule& molecule, const std::vector<int>& groupA, const std::vector<int>& groupB,
    const RdfOptions& options, std::vector<double>& histogram);

struct RdfResult {
    std::vector<double> r;              ///< bin centres
    std::vector<double> g;              ///< g(r)
    std::vector<double> gError;         ///< standard error from the blocks, 0 with fewer than 2 blocks
    std::vector<double> coordination;   ///< n(r): neighbours of B per atom of A inside the bin's upper edge
};

/// g(r) and n(r) from per-block histograms of ordered pairs. @p framesPerBlock[b] frames went into block b;
/// @p identicalGroups: A and B are the same atoms (the density uses N_B - 1); @p volume in A^3.
RdfResult normalizeRdf(const std::vector<std::vector<double>>& blockHistograms, const std::vector<int>& framesPerBlock,
    int atomsA, int atomsB, bool identicalGroups, double volume, double binWidth);

// --- solvent-accessible surface area ----------------------------------------------------------------
/// van der Waals radius in Angstrom (Bondi 1964) for the elements H, C, N, O, F, P, S, Cl, Br, I; 0 for any other.
double bondiRadius(int atomicNumber);

struct SasaOptions {
    double probe = 1.4;                 ///< probe sphere radius, Angstrom (1.4: a water molecule)
    int points = 960;                   ///< test points per atom
    std::vector<std::pair<int, double>> radii;   ///< per-element radius overrides (atomic number, Angstrom)
};

/// Solvent-accessible surface area (A^2) of @p surface, with @p occluders covering it (the atoms of @p surface
/// themselves are always occluders). Shrake-Rupley: an atom's area is 4 pi (r+probe)^2 times the fraction of
/// its test points outside every other atom's (r+probe) sphere. Fails (negative result, @p error set) for an
/// element without a radius.
double solventAccessibleSurface(const Molecule& molecule, const std::vector<int>& surface,
    const std::vector<int>& occluders, const SasaOptions& options, std::string* error);

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

    /// Quantities with their own loop; start() dispatches to them after the common checks.
    void startRdf();
    void startSasa();
    /// Read the trajectory once, frame by frame; @p visit gets the molecule, its frame number and the
    /// 0-based count of frames used, and returns false to stop (after setting m_error). False if stopped.
    bool streamFrames(const std::vector<int>& wanted, const std::function<bool(Molecule&, int, int)>& visit);
    /// Frames to use from the parameters, or empty with m_error set.
    std::vector<int> framesFromParameters(int* totalInFile);

    json m_results;
    std::string m_error;
};

} // namespace curcuma

namespace {
BEGIN_PARAMETER_DEFINITION(trajprops)
MODULE_INFO("Streaming properties of a long trajectory, one frame in memory at a time.",
    "Analysis", { "trajectory", "gyration" })

PARAM(quantity, String, "gyration", "What to compute per frame.", "Basic", {},
    "tier=primary; enum=gyration|hbonds|rdf|sasa")
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
PARAM(rdf_elements_a, String, "", "rdf: keep only these elements of the first group (symbols, comma separated); empty = all.", "RDF", {})
PARAM(rdf_elements_b, String, "", "rdf: keep only these elements of the second group; empty = all.", "RDF", {})
PARAM(rdf_r_max, Double, 10.0, "rdf: largest distance (Angstrom); with a cell at most half its shortest edge.", "RDF", {}, "min=0")
PARAM(rdf_bin_width, Double, 0.05, "rdf: bin width (Angstrom).", "RDF", {}, "min=0")
PARAM(cell, String, "", "rdf: orthorhombic cell \"Lx,Ly,Lz\" in Angstrom (minimum image convention and the volume for g(r)). An XYZ file has none.", "RDF", {})
PARAM(rdf_volume, Double, 0.0, "rdf without a cell: volume (A^3) for the density of the second group; the distances are then not wrapped.", "RDF", {}, "min=0")
PARAM(sasa_probe, Double, 1.4, "sasa: probe sphere radius (Angstrom); 1.4 is a water molecule.", "SASA", {}, "min=0")
PARAM(sasa_points, Int, 960, "sasa: test points per atom.", "SASA", {}, "min=10")
PARAM(sasa_radii, String, "", "sasa: radius overrides, e.g. \"Cl:1.75,O:1.6\" (Angstrom); the default set is Bondi (H, C, N, O, F, P, S, Cl, Br, I).", "SASA", {})
PARAM(sasa_occluders, String, "group", "sasa: what covers the surface atoms: group (the selected atoms only, the solute in a probe's view) or all (every atom, e.g. the solvent too).", "SASA", {}, "enum=group|all")
PARAM(blocks, Int, 5, "Number of blocks for the block standard error.", "Statistics", {}, "min=1")
PARAM(output, String, "", "Per-frame CSV; empty gives <basename>_<quantity>.csv.", "Output", {})
END_PARAMETER_DEFINITION
} // namespace
