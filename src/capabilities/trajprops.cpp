/*
 * TrajectoryProperties — see trajprops.h.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026
 */

#include "trajprops.h"

#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"
#include "src/core/fileiterator.h"
#include "src/core/parameter_registry.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>

namespace curcuma {

double gyrationRadius(const Molecule& molecule, const std::vector<int>& atoms, bool massWeighted)
{
    const std::vector<int> numbers = molecule.Atoms();
    std::vector<int> subset = atoms;
    if (subset.empty()) {
        subset.resize(numbers.size());
        for (std::size_t i = 0; i < numbers.size(); ++i)
            subset[i] = static_cast<int>(i);
    }
    if (subset.empty())
        return 0.0;

    const auto weight = [&](int i) { return massWeighted ? Elements::AtomicMass[numbers[i]] : 1.0; };
    double total = 0.0;
    Eigen::Vector3d centre = Eigen::Vector3d::Zero();
    for (int i : subset) {
        const Eigen::Vector3d p = molecule.Atom(i).second;
        centre += weight(i) * p;
        total += weight(i);
    }
    centre /= total;
    double second = 0.0;
    for (int i : subset) {
        const Eigen::Vector3d p = molecule.Atom(i).second;
        second += weight(i) * (p - centre).squaredNorm();
    }
    return std::sqrt(second / total);
}

std::vector<int> resolveFrames(int total, int first, int last, double lastFraction, int stride)
{
    std::vector<int> frames;
    if (total <= 0)
        return frames;
    int begin = std::max(0, first);
    int end = last < 0 ? total - 1 : std::min(last, total - 1);
    if (lastFraction > 0.0 && lastFraction <= 1.0) {
        const int count = std::max(1, static_cast<int>(std::ceil(lastFraction * total)));
        begin = total - count;
        end = total - 1;
    }
    for (int i = begin; i <= end; i += std::max(1, stride))
        frames.push_back(i);
    return frames;
}

SeriesStatistics seriesStatistics(const std::vector<double>& values, int blocks)
{
    SeriesStatistics s;
    s.n = static_cast<int>(values.size());
    if (values.empty())
        return s;
    double sum = 0.0;
    for (double v : values)
        sum += v;
    s.mean = sum / s.n;
    double squares = 0.0;
    for (double v : values)
        squares += (v - s.mean) * (v - s.mean);
    s.stddev = s.n > 1 ? std::sqrt(squares / (s.n - 1)) : 0.0;
    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    s.min = sorted.front();
    s.max = sorted.back();
    s.median = s.n % 2 ? sorted[s.n / 2] : 0.5 * (sorted[s.n / 2 - 1] + sorted[s.n / 2]);

    // Block averaging: consecutive blocks of equal length (a remainder at the end is dropped).
    const int used = std::max(1, std::min(blocks, s.n));
    const int length = s.n / used;
    s.blocks = used;
    if (used >= 2 && length >= 1) {
        std::vector<double> means;
        for (int b = 0; b < used; ++b) {
            double block = 0.0;
            for (int i = 0; i < length; ++i)
                block += values[b * length + i];
            means.push_back(block / length);
        }
        double m = 0.0;
        for (double v : means)
            m += v;
        m /= used;
        double var = 0.0;
        for (double v : means)
            var += (v - m) * (v - m);
        s.blockError = std::sqrt(var / (used - 1) / used);
    }
    return s;
}

TrajectoryProperties::TrajectoryProperties(const json& controller, bool silent)
    : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("trajprops"), controller, silent)
{
    UpdateController(controller);
}

void TrajectoryProperties::start()
{
    m_results = json();
    m_error.clear();

    const std::string quantity = Json2KeyWord<std::string>(m_defaults, "quantity");
    if (quantity != "gyration") {
        m_error = "unknown quantity \"" + quantity + "\" (gyration)";
        CurcumaLogger::error(m_error);
        return;
    }
    if (Filename().empty() || !std::filesystem::exists(Filename())) {
        m_error = "no such file: " + Filename();
        CurcumaLogger::error(m_error);
        return;
    }

    const std::string selection = Json2KeyWord<std::string>(m_defaults, "selection");
    const bool massWeighted = Json2KeyWord<bool>(m_defaults, "mass_weighted");

    FileIterator iterator(Filename(), true);
    const int total = iterator.MaxMolecules();
    const std::vector<int> wanted = resolveFrames(total, Json2KeyWord<int>(m_defaults, "first_frame"),
        Json2KeyWord<int>(m_defaults, "last_frame"), Json2KeyWord<double>(m_defaults, "last_fraction"),
        Json2KeyWord<int>(m_defaults, "stride"));
    if (wanted.empty()) {
        m_error = "no frame in the requested range (the file has " + std::to_string(total) + ")";
        CurcumaLogger::error(m_error);
        return;
    }

    std::vector<int> atoms;
    std::vector<double> values;
    std::vector<int> used;
    std::size_t next = 0;
    int frame = 0;
    while (!iterator.AtEnd() && next < wanted.size()) {
        if (frame != wanted[next]) {
            iterator.Next();   // outside the range: read and dropped
            ++frame;
            continue;
        }
        Molecule molecule = iterator.Next();
        if (next == 0 && !selection.empty()) {
            Molecule copy = molecule;
            copy.GetFragments();
            atoms = copy.FragString2Indicies(selection);
            if (atoms.empty()) {
                m_error = "selection \"" + selection + "\" matched no atoms";
                CurcumaLogger::error(m_error);
                return;
            }
        }
        for (int i : atoms) {
            if (i < 0 || i >= static_cast<int>(molecule.AtomCount())) {
                m_error = "frame " + std::to_string(frame) + " has fewer atoms than the selection needs";
                CurcumaLogger::error(m_error);
                return;
            }
        }
        values.push_back(gyrationRadius(molecule, atoms, massWeighted));
        used.push_back(frame);
        ++next;
        ++frame;
    }

    const SeriesStatistics s = seriesStatistics(values, Json2KeyWord<int>(m_defaults, "blocks"));
    std::string name = Json2KeyWord<std::string>(m_defaults, "output");
    if (name.empty())
        name = (Basename().empty() ? std::string("trajectory") : Basename()) + "_" + quantity + ".csv";
    {
        std::ofstream csv(outputPath(name));
        csv << "frame," << quantity << "_angstrom\n";
        for (std::size_t i = 0; i < values.size(); ++i)
            csv << used[i] << "," << values[i] << "\n";
    }

    m_results["quantity"] = quantity;
    m_results["unit"] = "Angstrom";
    m_results["mass_weighted"] = massWeighted;
    m_results["atoms"] = atoms.empty() ? -1 : static_cast<int>(atoms.size());   // -1: all atoms
    m_results["frames_in_file"] = total;
    m_results["frames_used"] = s.n;
    m_results["first_frame_used"] = used.front();
    m_results["last_frame_used"] = used.back();
    m_results["mean"] = s.mean;
    m_results["stddev"] = s.stddev;
    m_results["min"] = s.min;
    m_results["max"] = s.max;
    m_results["median"] = s.median;
    m_results["blocks"] = s.blocks;
    m_results["block_error"] = s.blockError;
    m_results["series_file"] = outputPath(name);

    CurcumaLogger::success(quantity + ": mean " + std::to_string(s.mean) + " A over " + std::to_string(s.n)
        + " frames (block error " + std::to_string(s.blockError) + ", " + std::to_string(s.blocks) + " blocks): "
        + outputPath(name));
}

} // namespace curcuma
