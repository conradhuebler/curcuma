/*
 * Measurement — one capability for the simple geometric quantities.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 - The pieces were all here and the bracket was missing.
 * `GeometryTools` has the primitives, `TrajectoryStatistics` has mean, deviation,
 * min, max and median, and `FileIterator` walks a trajectory — but six separate
 * functions in main.cpp each repeated the same argv parsing, the same existence
 * check, the same frame loop and three output formats by hand, with no PARAM
 * block, no schema and no Results() between them.
 *
 * This is that bracket. One Kind enum, one parameter block, one frame loop. The
 * six commands become thin dispatches onto it, and the same measurement is then
 * reachable from the CLI, from an embedded caller and from a tool layer without
 * three descriptions of it drifting apart.
 */

#pragma once

#include "curcumamethod.h"

#include "src/core/parameter_macros.h"

#include "src/capabilities/trajectory_statistics.h"
#include "src/core/molecule.h"

#include <string>
#include <vector>

namespace curcuma {

/// What is being measured. Which atoms a kind needs is a property of the kind,
/// not of the caller: requiredAtoms() is the single place that knows.
enum class MeasurementKind {
    Distance,          ///< between two atoms
    Angle,             ///< over three
    Dihedral,          ///< over four, signed
    Gyration,          ///< radius of gyration of an atom set (all atoms when empty)
    Centroid,          ///< the set's centroid; reported as a position, not a scalar
    RmsdToReference    ///< against the first frame, or against a given reference
};

/// How many atoms @p kind needs. 0 means "any number, including all of them".
int requiredAtoms(MeasurementKind kind);

/// Name of @p kind as it appears in the controller and in Results().
std::string measurementKindName(MeasurementKind kind);
/// Parse a name back; returns false when it is not one of the six.
bool parseMeasurementKind(const std::string& name, MeasurementKind& out);

class Measurement : public CurcumaMethod {
public:
    explicit Measurement(const json& controller, bool silent = true);
    ~Measurement() = default;

    /// Measure over a structure file or a trajectory. Every frame is measured.
    void setFile(const std::string& path) { m_file = path; }
    /// Measure a structure already in memory. Ignored when a file is set.
    void setMolecule(const Molecule& molecule);
    /// Reference for RmsdToReference. Without one the run's first frame is used.
    void setReference(const Molecule& reference);

    void start() override;

    /**
     * @brief The value per frame and the statistics over them.
     *
     * Pure: no file is written. `values` carries every frame, `statistics` the
     * mean, standard deviation, min, max and median, and `unit` says what the
     * numbers are. A Centroid measurement reports `positions` instead of
     * `values`, because a centroid is not a scalar and rounding it into one is
     * how a wrong number gets quoted later.
     */
    json Results() const override { return m_results; }

    /// Frames actually measured. 0 after a failure; @a error then says why.
    int frameCount() const { return m_frames; }
    std::string error() const { return m_error; }

private:
    nlohmann::json WriteRestartInformation() override { return json(); }
    bool LoadRestartInformation() override { return true; }
    StringList MethodName() const override { return { std::string("measurement") }; }
    void ReadControlFile() override {}
    void LoadControlJson() override;

    /// One frame. Returns false and fills m_error on an index out of range.
    bool measureFrame(const Molecule& molecule, double& value, Position& position);

    std::string m_file;
    Molecule m_molecule;
    bool m_have_molecule = false;
    Molecule m_reference;
    bool m_have_reference = false;

    MeasurementKind m_kind = MeasurementKind::Distance;
    std::vector<int> m_atoms;
    std::string m_unit;
    int m_window = 10;
    int m_first_frame = 0;
    int m_last_frame = -1;   ///< -1: to the end

    int m_frames = 0;
    std::string m_error;
    json m_results;
};

} // namespace curcuma

namespace {
BEGIN_PARAMETER_DEFINITION(measurement)
MODULE_INFO("Simple geometric measurements over a structure or a trajectory: distance, angle, "
            "dihedral, radius of gyration, centroid and RMSD, with statistics over the frames.",
    "Analysis", { "distance", "angle", "torsion", "bond", "gyration", "centroid" })

PARAM(kind, String, "distance", "What to measure.", "Basic", {},
    "tier=primary; enum=distance|angle|dihedral|gyration|centroid|rmsd")
PARAM(atoms, Selection, "", "Atoms to measure, in the selection grammar. Two for a distance, "
                            "three for an angle, four for a dihedral; any number for gyration "
                            "and centroid, empty meaning all of them.", "Basic", {},
    "tier=primary")
PARAM(unit, String, "degrees", "Unit for the angular kinds: degrees or radians. Lengths are "
                               "always Angstrom.", "Basic", {},
    "enum=degrees|radians")
PARAM(window, Int, 10, "Moving-average window for the per-frame statistics.", "Statistics", {},
    "min=1")
PARAM(first_frame, Int, 0, "First frame to measure (0-based).", "Frames", {},
    "min=0")
PARAM(last_frame, Int, -1, "Last frame to measure; -1 runs to the end.", "Frames", {})
END_PARAMETER_DEFINITION
} // namespace
