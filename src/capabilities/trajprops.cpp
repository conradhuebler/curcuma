/*
 * TrajectoryProperties — see trajprops.h.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026
 */

#include "trajprops.h"

#include "src/core/citation_registry.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energycalculator.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/units.h"
#include "src/core/elements.h"
#include "src/core/fileiterator.h"
#include "src/core/parameter_registry.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <unordered_map>

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

std::vector<int> resolveSelection(Molecule molecule, const std::string& selection)
{
    std::vector<int> atoms;
    if (selection.empty())
        return atoms;
    const bool complement = selection.front() == '!';
    molecule.GetFragments();   // FragString2Indicies reads this cache; without it every "Fn" matches nothing
    const std::vector<int> picked = molecule.FragString2Indicies(complement ? selection.substr(1) : selection);
    if (!complement)
        return picked;
    std::vector<bool> in(molecule.AtomCount(), false);
    for (int i : picked)
        if (i >= 0 && i < static_cast<int>(in.size()))
            in[i] = true;
    for (int i = 0; i < static_cast<int>(in.size()); ++i)
        if (!in[i])
            atoms.push_back(i);
    return atoms;
}

HBondCount countHydrogenBonds(const Molecule& molecule, const std::vector<int>& groupA,
    const std::vector<int>& groupB, const HBondCriteria& criteria)
{
    HBondCount count;
    const std::vector<int> numbers = molecule.Atoms();
    const int n = static_cast<int>(numbers.size());
    const Geometry geometry = molecule.getGeometry();
    const auto position = [&](int i) { return Eigen::Vector3d(geometry(i, 0), geometry(i, 1), geometry(i, 2)); };

    // 0: in neither group, bit 1: in A, bit 2: in B. With one group, A is all atoms when empty.
    std::vector<unsigned char> group(n, 0);
    if (groupA.empty()) {
        std::fill(group.begin(), group.end(), 1);
    } else {
        for (int i : groupA)
            if (i >= 0 && i < n)
                group[i] |= 1;
    }
    const bool twoGroups = !groupB.empty();
    for (int i : groupB)
        if (i >= 0 && i < n)
            group[i] |= 2;

    const auto capable = [&](int z) {
        return std::find(criteria.elements.begin(), criteria.elements.end(), z) != criteria.elements.end();
    };

    // Cell list of the donor/acceptor-capable heavy atoms. One grid serves both searches:
    // the cell edge is the larger of the two cut-offs.
    const double cell = std::max(criteria.distance, criteria.bondLength);
    // Exact cell key (21 bits per axis): a hash with collisions would scan one cell twice and double-count.
    const auto key = [&](const Eigen::Vector3d& p, int dx, int dy, int dz) {
        const auto axis = [&](double v, int d) {
            return static_cast<unsigned long long>(static_cast<long long>(std::floor(v / cell)) + d + (1LL << 20)) & 0x1FFFFFULL;
        };
        return static_cast<long long>((axis(p.x(), dx) << 42) | (axis(p.y(), dy) << 21) | axis(p.z(), dz));
    };
    std::unordered_map<long long, std::vector<int>> grid;
    for (int i = 0; i < n; ++i)
        if (capable(numbers[i]))
            grid[key(position(i), 0, 0, 0)].push_back(i);

    const double cosLimit = std::cos(criteria.angle * M_PI / 180.0);
    for (int h = 0; h < n; ++h) {
        if (numbers[h] != 1)
            continue;
        const Eigen::Vector3d ph = position(h);

        // The donor: nearest capable atom within bondLength.
        int donor = -1;
        double nearest = criteria.bondLength;
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    const auto it = grid.find(key(ph, dx, dy, dz));
                    if (it == grid.end())
                        continue;
                    for (int d : it->second) {
                        const double r = (position(d) - ph).norm();
                        if (r <= nearest) {
                            nearest = r;
                            donor = d;
                        }
                    }
                }
        if (donor < 0 || group[donor] == 0)
            continue;
        const Eigen::Vector3d pd = position(donor);

        // Acceptors around the donor.
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    const auto it = grid.find(key(pd, dx, dy, dz));
                    if (it == grid.end())
                        continue;
                    for (int a : it->second) {
                        if (a == donor)
                            continue;
                        const Eigen::Vector3d pa = position(a);
                        if ((pa - pd).norm() > criteria.distance)
                            continue;
                        const Eigen::Vector3d hd = pd - ph, ha = pa - ph;
                        const double cosine = hd.dot(ha) / (hd.norm() * ha.norm());
                        if (cosine > cosLimit)   // angle smaller than the limit
                            continue;
                        // Direction bookkeeping.
                        if (!twoGroups) {
                            if ((group[donor] & 1) && (group[a] & 1))
                                ++count.aDonates;
                        } else {
                            if ((group[donor] & 1) && (group[a] & 2))
                                ++count.aDonates;
                            else if ((group[donor] & 2) && (group[a] & 1))
                                ++count.bDonates;
                        }
                    }
                }
    }
    return count;
}

bool parseElementList(const std::string& list, std::vector<int>& numbers, std::string* error)
{
    numbers.clear();
    std::stringstream stream(list);
    std::string symbol;
    while (std::getline(stream, symbol, ',')) {
        symbol.erase(std::remove(symbol.begin(), symbol.end(), ' '), symbol.end());
        if (symbol.empty())
            continue;
        int z = 0;
        for (std::size_t e = 1; e < Elements::ElementAbbr.size(); ++e)
            if (Elements::ElementAbbr[e] == symbol)
                z = static_cast<int>(e);
        if (z == 0) {
            if (error)
                *error = "unknown element \"" + symbol + "\"";
            return false;
        }
        numbers.push_back(z);
    }
    return true;
}

std::vector<int> filterByElement(const Molecule& molecule, const std::vector<int>& atoms, const std::vector<int>& elements)
{
    const std::vector<int> numbers = molecule.Atoms();
    std::vector<int> pool = atoms;
    if (pool.empty()) {
        pool.resize(numbers.size());
        for (std::size_t i = 0; i < numbers.size(); ++i)
            pool[i] = static_cast<int>(i);
    }
    if (elements.empty())
        return pool;
    std::vector<int> kept;
    for (int i : pool)
        if (i >= 0 && i < static_cast<int>(numbers.size())
            && std::find(elements.begin(), elements.end(), numbers[i]) != elements.end())
            kept.push_back(i);
    return kept;
}

void accumulateRdf(const Molecule& molecule, const std::vector<int>& groupA, const std::vector<int>& groupB,
    const RdfOptions& options, std::vector<double>& histogram)
{
    const Geometry geometry = molecule.getGeometry();
    const double rMax2 = options.rMax * options.rMax;
    const int bins = static_cast<int>(histogram.size());
    for (int i : groupA) {
        for (int j : groupB) {
            if (i == j)
                continue;
            double d[3] = { geometry(j, 0) - geometry(i, 0), geometry(j, 1) - geometry(i, 1), geometry(j, 2) - geometry(i, 2) };
            if (options.periodic)
                for (int k = 0; k < 3; ++k)
                    d[k] -= options.cell[k] * std::round(d[k] / options.cell[k]);
            const double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            if (r2 >= rMax2)
                continue;
            const int bin = static_cast<int>(std::sqrt(r2) / options.binWidth);
            if (bin >= 0 && bin < bins)
                histogram[bin] += 1.0;
        }
    }
}

RdfResult normalizeRdf(const std::vector<std::vector<double>>& blockHistograms, const std::vector<int>& framesPerBlock,
    int atomsA, int atomsB, bool identicalGroups, double volume, double binWidth)
{
    RdfResult result;
    if (blockHistograms.empty())
        return result;
    const std::size_t bins = blockHistograms.front().size();
    const double density = (identicalGroups ? atomsB - 1 : atomsB) / volume;
    const auto shell = [&](std::size_t k) {
        const double inner = k * binWidth, outer = (k + 1) * binWidth;
        return 4.0 / 3.0 * M_PI * (outer * outer * outer - inner * inner * inner);
    };
    int frames = 0;
    for (int f : framesPerBlock)
        frames += f;

    std::vector<double> total(bins, 0.0);
    for (const auto& block : blockHistograms)
        for (std::size_t k = 0; k < bins; ++k)
            total[k] += block[k];

    result.r.resize(bins);
    result.g.resize(bins);
    result.gError.assign(bins, 0.0);
    result.coordination.resize(bins);
    double cumulative = 0.0;
    for (std::size_t k = 0; k < bins; ++k) {
        result.r[k] = (k + 0.5) * binWidth;
        result.g[k] = total[k] / (static_cast<double>(frames) * atomsA * density * shell(k));
        cumulative += total[k];
        result.coordination[k] = cumulative / (static_cast<double>(frames) * atomsA);
    }
    // Standard error of g(r) from the spread of the block values.
    std::vector<std::size_t> usable;
    for (std::size_t b = 0; b < blockHistograms.size(); ++b)
        if (framesPerBlock[b] > 0)
            usable.push_back(b);
    if (usable.size() >= 2) {
        const double n = static_cast<double>(usable.size());
        for (std::size_t k = 0; k < bins; ++k) {
            double mean = 0.0;
            std::vector<double> values;
            for (std::size_t b : usable) {
                values.push_back(blockHistograms[b][k] / (static_cast<double>(framesPerBlock[b]) * atomsA * density * shell(k)));
                mean += values.back();
            }
            mean /= n;
            double var = 0.0;
            for (double v : values)
                var += (v - mean) * (v - mean);
            result.gError[k] = std::sqrt(var / (n - 1.0) / n);
        }
    }
    return result;
}

double bondiRadius(int z)
{
    switch (z) {
    case 1: return 1.20;
    case 6: return 1.70;
    case 7: return 1.55;
    case 8: return 1.52;
    case 9: return 1.47;
    case 15: return 1.80;
    case 16: return 1.80;
    case 17: return 1.75;
    case 35: return 1.85;
    case 53: return 1.98;
    default: return 0.0;
    }
}

double solventAccessibleSurface(const Molecule& molecule, const std::vector<int>& surface,
    const std::vector<int>& occluders, const SasaOptions& options, std::string* error)
{
    const std::vector<int> numbers = molecule.Atoms();
    const Geometry geometry = molecule.getGeometry();
    const int n = static_cast<int>(numbers.size());

    // Radius per atom: override, else Bondi.
    std::vector<double> radius(n, 0.0);
    for (int i = 0; i < n; ++i) {
        radius[i] = bondiRadius(numbers[i]);
        for (const auto& o : options.radii)
            if (o.first == numbers[i])
                radius[i] = o.second;
    }
    // Everything that can cover: the occluders and the surface atoms themselves.
    std::vector<int> cover = occluders;
    cover.insert(cover.end(), surface.begin(), surface.end());
    std::sort(cover.begin(), cover.end());
    cover.erase(std::unique(cover.begin(), cover.end()), cover.end());
    for (int i : cover) {
        if (i < 0 || i >= n || radius[i] <= 0.0) {
            if (error)
                *error = "no van der Waals radius for element " + (i >= 0 && i < n ? Elements::ElementAbbr[numbers[i]] : std::string("?"))
                    + " (atom " + std::to_string(i + 1) + "); give it in sasa_radii";
            return -1.0;
        }
    }

    // Test points on a unit sphere: Fibonacci lattice, z = 1 - (2k+1)/N, azimuth k * golden angle.
    const int points = std::max(10, options.points);
    std::vector<Eigen::Vector3d> unit(points);
    const double golden = M_PI * (3.0 - std::sqrt(5.0));
    for (int k = 0; k < points; ++k) {
        const double z = 1.0 - (2.0 * k + 1.0) / points;
        const double rho = std::sqrt(std::max(0.0, 1.0 - z * z));
        unit[k] = Eigen::Vector3d(rho * std::cos(golden * k), rho * std::sin(golden * k), z);
    }

    double total = 0.0;
    std::vector<int> near;
    for (int i : surface) {
        const Eigen::Vector3d ci(geometry(i, 0), geometry(i, 1), geometry(i, 2));
        const double Ri = radius[i] + options.probe;
        near.clear();
        for (int j : cover) {
            if (j == i)
                continue;
            const Eigen::Vector3d cj(geometry(j, 0), geometry(j, 1), geometry(j, 2));
            const double reach = Ri + radius[j] + options.probe;
            if ((cj - ci).squaredNorm() < reach * reach)
                near.push_back(j);
        }
        int accessible = 0;
        for (int k = 0; k < points; ++k) {
            const Eigen::Vector3d p = ci + Ri * unit[k];
            bool free = true;
            for (int j : near) {
                const Eigen::Vector3d cj(geometry(j, 0), geometry(j, 1), geometry(j, 2));
                const double Rj = radius[j] + options.probe;
                if ((p - cj).squaredNorm() < Rj * Rj) {
                    free = false;
                    break;
                }
            }
            if (free)
                ++accessible;
        }
        total += 4.0 * M_PI * Ri * Ri * accessible / points;
    }
    return total;
}

HBondEnergyCount classifyHBondEnergies(const std::vector<GFNFFHBondEnergy>& triples, int atomCount,
    const std::vector<int>& groupA, const std::vector<int>& groupB, double limitKJ)
{
    HBondEnergyCount result;
    const int n = atomCount;
    std::vector<unsigned char> group(n, 0);
    if (groupA.empty()) {
        std::fill(group.begin(), group.end(), 1);
    } else {
        for (int i : groupA)
            if (i >= 0 && i < n)
                group[i] |= 1;
    }
    const bool twoGroups = !groupB.empty();
    for (int i : groupB)
        if (i >= 0 && i < n)
            group[i] |= 2;

    for (const GFNFFHBondEnergy& t : triples) {
        if (t.donor < 0 || t.donor >= n || t.acceptor < 0 || t.acceptor >= n)
            continue;
        const bool donorA = group[t.donor] & 1, donorB = group[t.donor] & 2;
        const bool acceptorA = group[t.acceptor] & 1, acceptorB = group[t.acceptor] & 2;
        bool aDonates = false, bDonates = false;
        if (!twoGroups) {
            aDonates = donorA && acceptorA;
        } else {
            aDonates = donorA && acceptorB;
            bDonates = !aDonates && donorB && acceptorA;
        }
        if (!aDonates && !bDonates)
            continue;
        const double kj = CurcumaUnit::Energy::hartree_to_kjmol(t.energy);
        result.energySumKJ += kj;
        if (kj != 0.0)
            result.energiesKJ.push_back(kj);
        if (kj < limitKJ) {
            if (aDonates)
                ++result.count.aDonates;
            else
                ++result.count.bDonates;
        }
    }
    result.ok = true;
    return result;
}

HBondEnergyCount countHydrogenBondsByEnergy(const Molecule& molecule, const std::vector<int>& groupA,
    const std::vector<int>& groupB, const HBondEnergyOptions& options)
{
    HBondEnergyCount result;
    Molecule copy = molecule;
    copy.setCharge(options.charge);
    json config = { { "verbosity", 0 }, { "threads", options.threads }, { "gfnff", json::object() } };
    EnergyCalculator calculator("gfnff", config);
    calculator.setMolecule(copy.getMolInfo());
    GFNFF* gfnff = calculator.Interface() ? calculator.Interface()->gfnffInstance() : nullptr;
    if (!gfnff) {
        result.error = "GFN-FF is not available for the energy criterion";
        return result;
    }
    gfnff->setRecordHBondEnergies(true);
    if (!std::isfinite(calculator.CalculateEnergy(false))) {
        result.error = "GFN-FF evaluation failed";
        return result;
    }
    const std::vector<GFNFFHBondEnergy> triples = gfnff->hydrogenBondEnergies();
    return classifyHBondEnergies(triples, static_cast<int>(molecule.AtomCount()), groupA, groupB, options.limitKJ);
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

    const auto fail = [this](const std::string& message) {
        m_error = message;
        CurcumaLogger::error(m_error);
    };

    const std::string quantity = Json2KeyWord<std::string>(m_defaults, "quantity");
    if (quantity != "gyration" && quantity != "hbonds" && quantity != "rdf" && quantity != "sasa") {
        fail("unknown quantity \"" + quantity + "\" (gyration, hbonds, rdf, sasa)");
        return;
    }
    if (Filename().empty() || !std::filesystem::exists(Filename())) {
        fail("no such file: " + Filename());
        return;
    }
    if (quantity == "rdf") {
        startRdf();
        return;
    }
    if (quantity == "sasa") {
        startSasa();
        return;
    }

    const std::string selection = Json2KeyWord<std::string>(m_defaults, "selection");
    const std::string selectionB = Json2KeyWord<std::string>(m_defaults, "selection_b");
    const bool massWeighted = Json2KeyWord<bool>(m_defaults, "mass_weighted");

    HBondCriteria criteria;
    if (quantity == "hbonds") {
        criteria.distance = Json2KeyWord<double>(m_defaults, "hb_distance");
        criteria.angle = Json2KeyWord<double>(m_defaults, "hb_angle");
        criteria.bondLength = Json2KeyWord<double>(m_defaults, "hb_bond_length");
        criteria.elements.clear();
        std::stringstream list(Json2KeyWord<std::string>(m_defaults, "hb_elements"));
        std::string symbol;
        while (std::getline(list, symbol, ',')) {
            int z = 0;
            for (std::size_t e = 1; e < Elements::ElementAbbr.size(); ++e)
                if (Elements::ElementAbbr[e] == symbol)
                    z = static_cast<int>(e);
            if (z == 0) {
                fail("hb_elements: unknown element \"" + symbol + "\"");
                return;
            }
            criteria.elements.push_back(z);
        }
        if (criteria.elements.empty()) {
            fail("hb_elements is empty");
            return;
        }
    }

    const std::string criterion = Json2KeyWord<std::string>(m_defaults, "hb_criterion");
    const bool byEnergy = quantity == "hbonds" && criterion == "energy";
    if (quantity == "hbonds" && criterion != "geometric" && criterion != "energy") {
        fail("hb_criterion must be geometric or energy, not \"" + criterion + "\"");
        return;
    }
    HBondEnergyOptions energyOptions;
    energyOptions.limitKJ = Json2KeyWord<double>(m_defaults, "hb_energy_limit");
    energyOptions.charge = Json2KeyWord<int>(m_defaults, "charge");
    energyOptions.threads = Json2KeyWord<int>(m_defaults, "threads");
    const int histBins = std::max(1, Json2KeyWord<int>(m_defaults, "hb_hist_bins"));
    const double histMax = Json2KeyWord<double>(m_defaults, "hb_hist_max");
    std::vector<long long> histogram(histBins, 0);
    long long belowRange = 0, aboveZero = 0;
    std::vector<double> energySums;

    FileIterator iterator(Filename(), true);
    const int total = iterator.MaxMolecules();
    const std::vector<int> wanted = resolveFrames(total, Json2KeyWord<int>(m_defaults, "first_frame"),
        Json2KeyWord<int>(m_defaults, "last_frame"), Json2KeyWord<double>(m_defaults, "last_fraction"),
        Json2KeyWord<int>(m_defaults, "stride"));
    if (wanted.empty()) {
        fail("no frame in the requested range (the file has " + std::to_string(total) + ")");
        return;
    }

    std::vector<int> atoms, atomsB;
    std::vector<double> values;              // gyration, or hbonds total
    std::vector<double> donatedA, donatedB;  // hbonds only
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
        if (next == 0) {
            atoms = resolveSelection(molecule, selection);
            if (!selection.empty() && atoms.empty()) {
                fail("selection \"" + selection + "\" matched no atoms");
                return;
            }
            atomsB = resolveSelection(molecule, selectionB);
            if (!selectionB.empty() && atomsB.empty()) {
                fail("selection_b \"" + selectionB + "\" matched no atoms");
                return;
            }
        }
        for (const std::vector<int>* group : { &atoms, &atomsB }) {
            for (int i : *group) {
                if (i < 0 || i >= static_cast<int>(molecule.AtomCount())) {
                    fail("frame " + std::to_string(frame) + " has fewer atoms than the selection needs");
                    return;
                }
            }
        }
        if (quantity == "gyration") {
            values.push_back(gyrationRadius(molecule, atoms, massWeighted));
        } else if (byEnergy) {
            const HBondEnergyCount c = countHydrogenBondsByEnergy(molecule, atoms, atomsB, energyOptions);
            if (!c.ok) {
                fail("frame " + std::to_string(frame) + ": " + c.error);
                return;
            }
            values.push_back(c.count.total());
            donatedA.push_back(c.count.aDonates);
            donatedB.push_back(c.count.bDonates);
            energySums.push_back(c.energySumKJ);
            for (double kj : c.energiesKJ) {
                if (kj > 0.0) {
                    ++aboveZero;
                } else if (kj < -histMax) {
                    ++belowRange;
                } else {
                    const int bin = std::min(histBins - 1, static_cast<int>((kj + histMax) / histMax * histBins));
                    ++histogram[bin];
                }
            }
        } else {
            const HBondCount c = countHydrogenBonds(molecule, atoms, atomsB, criteria);
            values.push_back(c.total());
            donatedA.push_back(c.aDonates);
            donatedB.push_back(c.bDonates);
        }
        used.push_back(frame);
        ++next;
        ++frame;
    }

    const int blocks = Json2KeyWord<int>(m_defaults, "blocks");
    const SeriesStatistics s = seriesStatistics(values, blocks);
    std::string name = Json2KeyWord<std::string>(m_defaults, "output");
    if (name.empty())
        name = (Basename().empty() ? std::string("trajectory") : Basename()) + "_" + quantity + ".csv";
    {
        std::ofstream csv(outputPath(name));
        if (quantity == "gyration")
            csv << "frame,gyration_angstrom\n";
        else
            csv << "frame,hbonds_total,a_donates,b_donates" << (byEnergy ? ",energy_sum_kj_mol" : "") << "\n";
        for (std::size_t i = 0; i < values.size(); ++i) {
            csv << used[i] << "," << values[i];
            if (quantity == "hbonds")
                csv << "," << donatedA[i] << "," << donatedB[i];
            if (byEnergy)
                csv << "," << energySums[i];
            csv << "\n";
        }
    }
    std::string histogramName;
    if (byEnergy) {
        // Triples per frame in each energy bin (counts summed over the frames, divided by the frame count).
        histogramName = (Basename().empty() ? std::string("trajectory") : Basename()) + "_hbonds_energy_histogram.csv";
        std::ofstream csv(outputPath(histogramName));
        csv << "bin_low_kj_mol,bin_high_kj_mol,triples_per_frame\n";
        for (int b = 0; b < histBins; ++b)
            csv << (-histMax + histMax * b / histBins) << "," << (-histMax + histMax * (b + 1) / histBins) << ","
                << static_cast<double>(histogram[b]) / static_cast<double>(values.size()) << "\n";
    }

    m_results["quantity"] = quantity;
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
    if (quantity == "gyration") {
        m_results["unit"] = "Angstrom";
        m_results["mass_weighted"] = massWeighted;
    } else {
        m_results["unit"] = "hydrogen bonds per frame";
        m_results["atoms_b"] = atomsB.empty() ? -1 : static_cast<int>(atomsB.size());
        m_results["criterion"] = criterion;
        if (byEnergy) {
            m_results["energy_limit_kj_mol"] = energyOptions.limitKJ;
            m_results["charge"] = energyOptions.charge;
            const SeriesStatistics es = seriesStatistics(energySums, blocks);
            m_results["energy_sum_kj_mol"] = { { "mean", es.mean }, { "stddev", es.stddev }, { "block_error", es.blockError } };
            m_results["histogram_file"] = outputPath(histogramName);
            m_results["triples_above_zero_per_frame"] = static_cast<double>(aboveZero) / static_cast<double>(values.size());
            m_results["triples_below_histogram_range_per_frame"] = static_cast<double>(belowRange) / static_cast<double>(values.size());
        } else {
            m_results["criteria"] = { { "distance", criteria.distance }, { "angle", criteria.angle },
                { "bond_length", criteria.bondLength } };
        }
        const SeriesStatistics a = seriesStatistics(donatedA, blocks);
        const SeriesStatistics b = seriesStatistics(donatedB, blocks);
        m_results["a_donates"] = { { "mean", a.mean }, { "stddev", a.stddev }, { "block_error", a.blockError } };
        m_results["b_donates"] = { { "mean", b.mean }, { "stddev", b.stddev }, { "block_error", b.blockError } };
    }

    // References of the methods used (kept in the citation database; TODO: check them against the papers).
    if (s.n > 1)
        CitationRegistry::cite("flyvbjerg1989");
    if (quantity == "hbonds") {
        if (byEnergy) {
            CitationRegistry::cite("gfnff");
        } else {
            CitationRegistry::cite("luzar1996");
            CitationRegistry::cite("kumar2007hbond");
        }
    }

    CurcumaLogger::success(quantity + ": mean " + std::to_string(s.mean) + " over " + std::to_string(s.n)
        + " frames (block error " + std::to_string(s.blockError) + ", " + std::to_string(s.blocks) + " blocks): "
        + outputPath(name));
}


std::vector<int> TrajectoryProperties::framesFromParameters(int* totalInFile)
{
    FileIterator iterator(Filename(), true);
    const int total = iterator.MaxMolecules();
    if (totalInFile)
        *totalInFile = total;
    const std::vector<int> wanted = resolveFrames(total, Json2KeyWord<int>(m_defaults, "first_frame"),
        Json2KeyWord<int>(m_defaults, "last_frame"), Json2KeyWord<double>(m_defaults, "last_fraction"),
        Json2KeyWord<int>(m_defaults, "stride"));
    if (wanted.empty()) {
        m_error = "no frame in the requested range (the file has " + std::to_string(total) + ")";
        CurcumaLogger::error(m_error);
    }
    return wanted;
}

bool TrajectoryProperties::streamFrames(const std::vector<int>& wanted, const std::function<bool(Molecule&, int, int)>& visit)
{
    FileIterator iterator(Filename(), true);
    std::size_t next = 0;
    int frame = 0;
    while (!iterator.AtEnd() && next < wanted.size()) {
        if (frame != wanted[next]) {
            iterator.Next();   // outside the range: read and dropped
            ++frame;
            continue;
        }
        Molecule molecule = iterator.Next();
        if (!visit(molecule, frame, static_cast<int>(next)))
            return false;
        ++next;
        ++frame;
    }
    return true;
}

void TrajectoryProperties::startRdf()
{
    const auto fail = [this](const std::string& message) {
        m_error = message;
        CurcumaLogger::error(m_error);
    };

    RdfOptions options;
    options.rMax = Json2KeyWord<double>(m_defaults, "rdf_r_max");
    options.binWidth = Json2KeyWord<double>(m_defaults, "rdf_bin_width");
    if (options.rMax <= 0.0 || options.binWidth <= 0.0 || options.binWidth > options.rMax) {
        fail("rdf_r_max and rdf_bin_width must be positive, the bin not wider than r_max");
        return;
    }
    double volume = Json2KeyWord<double>(m_defaults, "rdf_volume");
    const std::string cell = Json2KeyWord<std::string>(m_defaults, "cell");
    if (!cell.empty()) {
        std::stringstream stream(cell);
        std::string part;
        int k = 0;
        while (std::getline(stream, part, ',') && k < 3) {
            try {
                options.cell[k++] = std::stod(part);
            } catch (...) {
                k = 99;
            }
        }
        if (k != 3 || options.cell[0] <= 0.0 || options.cell[1] <= 0.0 || options.cell[2] <= 0.0) {
            fail("cell must be three positive edges \"Lx,Ly,Lz\" in Angstrom, got \"" + cell + "\"");
            return;
        }
        options.periodic = true;
        volume = options.cell[0] * options.cell[1] * options.cell[2];
        if (options.rMax > 0.5 * std::min({ options.cell[0], options.cell[1], options.cell[2] })) {
            fail("rdf_r_max must not exceed half the shortest cell edge (minimum image convention)");
            return;
        }
    }
    if (volume <= 0.0) {
        fail("rdf needs the volume for the density: give cell \"Lx,Ly,Lz\" or rdf_volume (an XYZ file has no cell)");
        return;
    }

    std::vector<int> elementsA, elementsB;
    std::string parseError;
    if (!parseElementList(Json2KeyWord<std::string>(m_defaults, "rdf_elements_a"), elementsA, &parseError)
        || !parseElementList(Json2KeyWord<std::string>(m_defaults, "rdf_elements_b"), elementsB, &parseError)) {
        fail("rdf_elements: " + parseError);
        return;
    }
    const std::string selection = Json2KeyWord<std::string>(m_defaults, "selection");
    const std::string selectionB = Json2KeyWord<std::string>(m_defaults, "selection_b");

    int total = 0;
    const std::vector<int> wanted = framesFromParameters(&total);
    if (wanted.empty())
        return;
    const int nBlocks = std::max(1, std::min(Json2KeyWord<int>(m_defaults, "blocks"), static_cast<int>(wanted.size())));
    const int bins = static_cast<int>(std::floor(options.rMax / options.binWidth));
    std::vector<std::vector<double>> blockHistograms(nBlocks, std::vector<double>(bins, 0.0));
    std::vector<int> framesPerBlock(nBlocks, 0);
    std::vector<int> groupA, groupB;
    bool identical = false;
    std::vector<int> used;

    const bool finished = streamFrames(wanted, [&](Molecule& molecule, int frame, int k) {
        if (k == 0) {
            const std::vector<int> baseA = resolveSelection(molecule, selection);
            if (!selection.empty() && baseA.empty()) {
                fail("selection \"" + selection + "\" matched no atoms");
                return false;
            }
            groupA = filterByElement(molecule, baseA, elementsA);
            if (selectionB.empty() && elementsB.empty()) {
                groupB = groupA;
            } else {
                const std::vector<int> baseB = selectionB.empty() ? baseA : resolveSelection(molecule, selectionB);
                if (!selectionB.empty() && baseB.empty()) {
                    fail("selection_b \"" + selectionB + "\" matched no atoms");
                    return false;
                }
                groupB = filterByElement(molecule, baseB, elementsB);
            }
            if (groupA.empty() || groupB.empty()) {
                fail("a group is empty after the selection and the element filter");
                return false;
            }
            identical = groupA == groupB;
            if (!identical) {
                std::vector<int> common;
                std::set_intersection(groupA.begin(), groupA.end(), groupB.begin(), groupB.end(), std::back_inserter(common));
                if (!common.empty()) {
                    fail("the two groups overlap without being identical (" + std::to_string(common.size())
                        + " shared atoms); the density of the second group would be ambiguous");
                    return false;
                }
            }
        }
        const int block = std::min(nBlocks - 1, k * nBlocks / static_cast<int>(wanted.size()));
        accumulateRdf(molecule, groupA, groupB, options, blockHistograms[block]);
        ++framesPerBlock[block];
        used.push_back(frame);
        return true;
    });
    if (!finished)
        return;

    const RdfResult rdf = normalizeRdf(blockHistograms, framesPerBlock, static_cast<int>(groupA.size()),
        static_cast<int>(groupB.size()), identical, volume, options.binWidth);
    std::string name = Json2KeyWord<std::string>(m_defaults, "output");
    if (name.empty())
        name = (Basename().empty() ? std::string("trajectory") : Basename()) + "_rdf.csv";
    {
        std::ofstream csv(outputPath(name));
        csv << "r_angstrom,g,g_error,coordination\n";
        for (std::size_t k = 0; k < rdf.r.size(); ++k)
            csv << rdf.r[k] << "," << rdf.g[k] << "," << rdf.gError[k] << "," << rdf.coordination[k] << "\n";
    }
    m_results["quantity"] = "rdf";
    m_results["atoms_a"] = static_cast<int>(groupA.size());
    m_results["atoms_b"] = static_cast<int>(groupB.size());
    m_results["identical_groups"] = identical;
    m_results["frames_in_file"] = total;
    m_results["frames_used"] = static_cast<int>(used.size());
    m_results["first_frame_used"] = used.front();
    m_results["last_frame_used"] = used.back();
    m_results["r_max"] = options.rMax;
    m_results["bin_width"] = options.binWidth;
    m_results["periodic"] = options.periodic;
    m_results["volume_angstrom3"] = volume;
    m_results["blocks"] = nBlocks;
    m_results["coordination_at_r_max"] = rdf.coordination.empty() ? 0.0 : rdf.coordination.back();
    m_results["series_file"] = outputPath(name);

    CitationRegistry::cite("allen_tildesley");
    if (used.size() > 1)
        CitationRegistry::cite("flyvbjerg1989");
    CurcumaLogger::success("rdf: " + std::to_string(groupA.size()) + " x " + std::to_string(groupB.size()) + " atoms over "
        + std::to_string(used.size()) + " frames, n(r_max) = " + std::to_string(rdf.coordination.back()) + ": " + outputPath(name));
}

void TrajectoryProperties::startSasa()
{
    const auto fail = [this](const std::string& message) {
        m_error = message;
        CurcumaLogger::error(m_error);
    };

    SasaOptions options;
    options.probe = Json2KeyWord<double>(m_defaults, "sasa_probe");
    options.points = Json2KeyWord<int>(m_defaults, "sasa_points");
    std::stringstream overrides(Json2KeyWord<std::string>(m_defaults, "sasa_radii"));
    std::string item;
    while (std::getline(overrides, item, ',')) {
        const std::size_t colon = item.find(':');
        std::vector<int> z;
        std::string error;
        double value = 0.0;
        try {
            value = colon == std::string::npos ? 0.0 : std::stod(item.substr(colon + 1));
        } catch (...) {
            value = 0.0;
        }
        if (colon == std::string::npos || value <= 0.0 || !parseElementList(item.substr(0, colon), z, &error) || z.size() != 1) {
            fail("sasa_radii: cannot read \"" + item + "\" (expected Element:radius, e.g. Cl:1.75)");
            return;
        }
        options.radii.push_back({ z.front(), value });
    }
    const std::string occluderMode = Json2KeyWord<std::string>(m_defaults, "sasa_occluders");
    if (occluderMode != "group" && occluderMode != "all") {
        fail("sasa_occluders must be group or all, not \"" + occluderMode + "\"");
        return;
    }
    const std::string selection = Json2KeyWord<std::string>(m_defaults, "selection");

    int total = 0;
    const std::vector<int> wanted = framesFromParameters(&total);
    if (wanted.empty())
        return;
    std::vector<int> surface, occluders, used;
    std::vector<double> values;
    const bool finished = streamFrames(wanted, [&](Molecule& molecule, int frame, int k) {
        if (k == 0) {
            surface = resolveSelection(molecule, selection);
            if (!selection.empty() && surface.empty()) {
                fail("selection \"" + selection + "\" matched no atoms");
                return false;
            }
            if (surface.empty())
                surface = filterByElement(molecule, {}, {});
            if (occluderMode == "all")
                occluders = filterByElement(molecule, {}, {});
        }
        std::string error;
        const double area = solventAccessibleSurface(molecule, surface, occluders, options, &error);
        if (area < 0.0) {
            fail(error);
            return false;
        }
        values.push_back(area);
        used.push_back(frame);
        return true;
    });
    if (!finished)
        return;

    const int blocks = Json2KeyWord<int>(m_defaults, "blocks");
    const SeriesStatistics s = seriesStatistics(values, blocks);
    std::string name = Json2KeyWord<std::string>(m_defaults, "output");
    if (name.empty())
        name = (Basename().empty() ? std::string("trajectory") : Basename()) + "_sasa.csv";
    {
        std::ofstream csv(outputPath(name));
        csv << "frame,sasa_angstrom2\n";
        for (std::size_t i = 0; i < values.size(); ++i)
            csv << used[i] << "," << values[i] << "\n";
    }
    m_results["quantity"] = "sasa";
    m_results["unit"] = "Angstrom^2";
    m_results["atoms"] = static_cast<int>(surface.size());
    m_results["occluders"] = occluderMode;
    m_results["probe"] = options.probe;
    m_results["points_per_atom"] = options.points;
    m_results["radii"] = "Bondi 1964 (H C N O F P S Cl Br I)" + std::string(options.radii.empty() ? "" : ", with overrides from sasa_radii");
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

    CitationRegistry::cite("shrake1973");
    CitationRegistry::cite("lee1971");
    CitationRegistry::cite("bondi1964");
    CitationRegistry::cite("gonzalez2010fibonacci");
    if (s.n > 1)
        CitationRegistry::cite("flyvbjerg1989");
    CurcumaLogger::success("sasa: mean " + std::to_string(s.mean) + " A^2 over " + std::to_string(s.n) + " frames (block error "
        + std::to_string(s.blockError) + "): " + outputPath(name));
}

} // namespace curcuma
