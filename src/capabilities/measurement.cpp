/*
 * Measurement — one capability for the simple geometric quantities.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated 2026
 */

#include "measurement.h"

#include "src/capabilities/rmsd/rmsd_functions.h"
#include "src/core/curcuma_logger.h"
#include "src/core/fileiterator.h"
#include "src/core/parameter_registry.h"
#include "src/tools/geometry.h"

#include <algorithm>
#include <cmath>
#include <filesystem>

namespace curcuma {

int requiredAtoms(MeasurementKind kind)
{
    switch (kind) {
    case MeasurementKind::Distance: return 2;
    case MeasurementKind::Angle:    return 3;
    case MeasurementKind::Dihedral: return 4;
    case MeasurementKind::Gyration:
    case MeasurementKind::Centroid:
    case MeasurementKind::CentroidDistance:
    case MeasurementKind::RmsdToReference: return 0;
    }
    return 0;
}

std::string measurementKindName(MeasurementKind kind)
{
    switch (kind) {
    case MeasurementKind::Distance: return "distance";
    case MeasurementKind::Angle:    return "angle";
    case MeasurementKind::Dihedral: return "dihedral";
    case MeasurementKind::Gyration: return "gyration";
    case MeasurementKind::Centroid: return "centroid";
    case MeasurementKind::CentroidDistance: return "centroid_distance";
    case MeasurementKind::RmsdToReference: return "rmsd";
    }
    return "distance";
}

bool parseMeasurementKind(const std::string& name, MeasurementKind& out)
{
    if (name == "distance") { out = MeasurementKind::Distance; return true; }
    if (name == "angle")    { out = MeasurementKind::Angle; return true; }
    if (name == "dihedral" || name == "torsion") { out = MeasurementKind::Dihedral; return true; }
    if (name == "gyration") { out = MeasurementKind::Gyration; return true; }
    if (name == "centroid") { out = MeasurementKind::Centroid; return true; }
    if (name == "centroid_distance") { out = MeasurementKind::CentroidDistance; return true; }
    if (name == "rmsd")     { out = MeasurementKind::RmsdToReference; return true; }
    return false;
}

Measurement::Measurement(const json& controller, bool silent)
    : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("measurement"), controller,
        silent)
{
    UpdateController(controller);
}

void Measurement::setMolecule(const Molecule& molecule)
{
    m_molecule = molecule;
    m_have_molecule = true;
}

void Measurement::setReference(const Molecule& reference)
{
    m_reference = reference;
    m_have_reference = true;
}

void Measurement::LoadControlJson()
{
    const std::string kind = Json2KeyWord<std::string>(m_defaults, "kind");
    if (!parseMeasurementKind(kind, m_kind)) {
        m_error = "unknown kind \"" + kind + "\"";
        return;
    }
    m_unit = Json2KeyWord<std::string>(m_defaults, "unit");
    m_dihedral_positive = Json2KeyWord<std::string>(m_defaults, "dihedral_range") == "positive";
    m_window = Json2KeyWord<int>(m_defaults, "window");
    m_first_frame = Json2KeyWord<int>(m_defaults, "first_frame");
    m_last_frame = Json2KeyWord<int>(m_defaults, "last_frame");

    // The selection is resolved per frame, not here: the atom count is a property
    // of the structure and a file has not been opened yet.
    m_atoms.clear();
}

bool Measurement::measureFrame(const Molecule& molecule, double& value, Position& position)
{
    // The selection is resolved against this frame. GetFragments() first, always:
    // FragString2Indicies reads a cache that is otherwise empty and every "Fn"
    // matches nothing without a word.
    std::vector<int> atoms = m_atoms;
    const std::string selection = Json2KeyWord<std::string>(m_defaults, "atoms");
    if (!selection.empty()) {
        Molecule copy = molecule;
        copy.GetFragments();
        atoms = copy.FragString2Indicies(selection);
        if (atoms.empty()) {
            m_error = "selection \"" + selection + "\" matched no atoms";
            return false;
        }
    }

    const int needed = requiredAtoms(m_kind);
    if (needed > 0 && static_cast<int>(atoms.size()) != needed) {
        m_error = measurementKindName(m_kind) + " needs exactly " + std::to_string(needed)
            + " atoms, the selection gave " + std::to_string(atoms.size());
        return false;
    }
    for (int index : atoms) {
        if (index < 0 || index >= molecule.AtomCount()) {
            m_error = "atom index " + std::to_string(index) + " is outside 0.."
                + std::to_string(molecule.AtomCount() - 1);
            return false;
        }
    }

    const bool radians = m_unit == "radians";
    const auto toUnit = [radians](double degrees) {
        return radians ? degrees * pi / 180.0 : degrees;
    };

    switch (m_kind) {
    case MeasurementKind::Distance:
        value = GeometryTools::Distance(molecule.Atom(atoms[0]).second,
            molecule.Atom(atoms[1]).second);
        return true;
    case MeasurementKind::Angle:
        value = toUnit(GeometryTools::Angle(molecule.Atom(atoms[0]).second,
            molecule.Atom(atoms[1]).second, molecule.Atom(atoms[2]).second));
        return true;
    case MeasurementKind::Dihedral: {
        double degrees = GeometryTools::Dihedral(molecule.Atom(atoms[0]).second,
            molecule.Atom(atoms[1]).second, molecule.Atom(atoms[2]).second,
            molecule.Atom(atoms[3]).second);
        // The CLI's -torsion has always reported [0, 360). Signed is the IUPAC
        // convention and the better default, but changing what an existing command
        // prints is not a side effect anyone should get for free.
        if (m_dihedral_positive && degrees < 0.0)
            degrees += 360.0;
        value = toUnit(degrees);
        return true;
    }
    case MeasurementKind::Gyration:
    case MeasurementKind::Centroid: {
        Geometry subset(atoms.empty() ? molecule.AtomCount() : int(atoms.size()), 3);
        if (atoms.empty()) {
            subset = molecule.getGeometry();
        } else {
            for (size_t i = 0; i < atoms.size(); ++i) {
                const Position p = molecule.Atom(atoms[i]).second;
                subset(int(i), 0) = p(0);
                subset(int(i), 1) = p(1);
                subset(int(i), 2) = p(2);
            }
        }
        if (m_kind == MeasurementKind::Gyration)
            value = GeometryTools::GyrationRadius(subset);
        else
            position = GeometryTools::Centroid(subset);
        return true;
    }
    case MeasurementKind::CentroidDistance: {
        const std::string other = Json2KeyWord<std::string>(m_defaults, "atoms_b");
        if (other.empty()) {
            m_error = "centroid_distance needs a second selection in atoms_b";
            return false;
        }
        Molecule copy = molecule;
        copy.GetFragments();
        const std::vector<int> setB = copy.FragString2Indicies(other);
        if (setB.empty()) {
            m_error = "selection \"" + other + "\" matched no atoms";
            return false;
        }
        const auto centreOf = [&molecule](const std::vector<int>& set) {
            Position c { 0, 0, 0 };
            for (int index : set)
                c += molecule.Atom(index).second;
            return set.empty() ? c : Position(c / double(set.size()));
        };
        for (int index : setB) {
            if (index < 0 || index >= molecule.AtomCount()) {
                m_error = "atom index " + std::to_string(index) + " is outside the structure";
                return false;
            }
        }
        value = (centreOf(atoms) - centreOf(setB)).norm();
        return true;
    }
    case MeasurementKind::RmsdToReference: {
        if (!m_have_reference) {
            m_error = "no reference structure for an rmsd measurement";
            return false;
        }
        if (m_reference.AtomCount() != molecule.AtomCount()) {
            m_error = "the reference has " + std::to_string(m_reference.AtomCount())
                + " atoms and this frame " + std::to_string(molecule.AtomCount());
            return false;
        }
        // Centroids removed, then the best-fit rotation: drift and tumbling are not
        // a change of structure.
        Geometry a = m_reference.getGeometry();
        Geometry b = molecule.getGeometry();
        a = GeometryTools::TranslateGeometry(a, GeometryTools::Centroid(a), Position { 0, 0, 0 });
        b = GeometryTools::TranslateGeometry(b, GeometryTools::Centroid(b), Position { 0, 0, 0 });
        value = RMSDFunctions::getRMSD(a, RMSDFunctions::applyRotation(b,
            RMSDFunctions::BestFitRotation(a, b)));
        return true;
    }
    }
    return false;
}

void Measurement::start()
{
    m_results = json();
    m_frames = 0;
    if (!m_error.empty())
        return;

    std::vector<Molecule> frames;
    if (!m_file.empty()) {
        if (!std::filesystem::exists(m_file)) {
            m_error = "no such file: " + m_file;
            return;
        }
        FileIterator iterator(m_file, true);
        int index = 0;
        while (!iterator.AtEnd()) {
            Molecule molecule = iterator.Next();
            const bool wanted = index >= m_first_frame
                && (m_last_frame < 0 || index <= m_last_frame);
            if (wanted)
                frames.push_back(molecule);
            ++index;
        }
    } else if (m_have_molecule) {
        frames.push_back(m_molecule);
    }

    if (frames.empty()) {
        if (m_error.empty())
            m_error = "nothing to measure: give a file or a structure";
        return;
    }
    // Without an explicit reference, an RMSD measurement is against where the run
    // started -- which is what "how far has it moved" means on a trajectory.
    if (m_kind == MeasurementKind::RmsdToReference && !m_have_reference)
        setReference(frames.front());

    const std::string name = measurementKindName(m_kind);
    TrajectoryStatistics statistics(std::max(1, m_window));
    json values = json::array();
    json positions = json::array();

    for (const Molecule& molecule : frames) {
        double value = 0.0;
        Position position { 0, 0, 0 };
        if (!measureFrame(molecule, value, position))
            return;
        if (m_kind == MeasurementKind::Centroid) {
            positions.push_back(json::array({ position(0), position(1), position(2) }));
        } else {
            values.push_back(value);
            statistics.addValue(name, value);
        }
        ++m_frames;
    }

    m_results["kind"] = name;
    m_results["frames"] = m_frames;
    m_results["unit"] = (m_kind == MeasurementKind::Angle || m_kind == MeasurementKind::Dihedral)
        ? m_unit
        : std::string("Angstrom");
    if (m_kind == MeasurementKind::Centroid) {
        // A centroid is a position, not a scalar. Reducing it to one and running
        // statistics over that is how a meaningless number gets quoted later.
        m_results["positions"] = positions;
    } else {
        m_results["values"] = values;
        if (m_frames > 1)
            m_results["statistics"] = statistics.exportStatistics(name);
    }
}

} // namespace curcuma
