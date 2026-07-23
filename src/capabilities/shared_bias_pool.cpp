/*
 * <Shared Bias Pool for Parallel RMSD-MTD Conformational Search>
 * Copyright (C) 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Apr 2026)
 */

#include "src/capabilities/shared_bias_pool.h"

#include "src/core/fileiterator.h"
#include "src/core/molecule.h"

#include <algorithm>
#include <sstream>

std::vector<BiasStructure> SharedBiasPool::snapshot() const
{
    std::shared_lock<std::shared_mutex> lock(m_mutex);
    return m_structures;
}

int SharedBiasPool::biasStructureCount() const
{
    return m_global_count.load(std::memory_order_relaxed);
}

int SharedBiasPool::depositBiasStructure(const BiasStructure& structure)
{
    std::unique_lock<std::shared_mutex> lock(m_mutex);
    int index = static_cast<int>(m_structures.size());
    BiasStructure bs = structure;
    bs.index = index;
    m_structures.push_back(std::move(bs));
    m_global_count.store(static_cast<int>(m_structures.size()), std::memory_order_release);
    return index;
}

int SharedBiasPool::depositBatch(const std::vector<BiasStructure>& structures)
{
    std::unique_lock<std::shared_mutex> lock(m_mutex);
    int first_index = static_cast<int>(m_structures.size());
    for (int i = 0; i < static_cast<int>(structures.size()); ++i) {
        BiasStructure bs = structures[i];
        bs.index = first_index + i;
        m_structures.push_back(std::move(bs));
    }
    m_global_count.store(static_cast<int>(m_structures.size()), std::memory_order_release);
    return first_index;
}

void SharedBiasPool::registerVisits(const std::vector<std::pair<int, double>>& updates)
{
    if (updates.empty())
        return;
    std::unique_lock<std::shared_mutex> lock(m_mutex);
    for (const auto& [idx, wt_weight] : updates) {
        if (idx >= 0 && idx < static_cast<int>(m_structures.size())) {
            m_structures[idx].counter++;          // exploration: hill height W = k*counter
            m_structures[idx].factor += wt_weight; // opt-in well-tempered output weight
        }
    }
}

void SharedBiasPool::pruneByCounter(int min_counter)
{
    std::unique_lock<std::shared_mutex> lock(m_mutex);
    auto it = std::remove_if(m_structures.begin(), m_structures.end(),
        [min_counter](const BiasStructure& bs) {
            // Persistent structures (fed-back optimised minima) are never pruned -- they
            // represent real basins the search should keep biasing against. Claude Generated (Jun 2026).
            return bs.counter < min_counter && !bs.persistent;
        });
    m_structures.erase(it, m_structures.end());
    // Re-index after pruning
    for (int i = 0; i < static_cast<int>(m_structures.size()); ++i) {
        m_structures[i].index = i;
    }
    m_global_count.store(static_cast<int>(m_structures.size()), std::memory_order_release);
}

void SharedBiasPool::pruneNonPersistent()
{
    std::unique_lock<std::shared_mutex> lock(m_mutex);
    auto it = std::remove_if(m_structures.begin(), m_structures.end(),
        [](const BiasStructure& bs) { return !bs.persistent; });
    m_structures.erase(it, m_structures.end());
    for (int i = 0; i < static_cast<int>(m_structures.size()); ++i)
        m_structures[i].index = i;
    m_global_count.store(static_cast<int>(m_structures.size()), std::memory_order_release);
}

int SharedBiasPool::capToSize(int max_size)
{
    // Claude Generated (Jul 2026): enforce rmsd_mtd_max_gaussians. Keep every persistent (fed-back
    // optimised) minimum plus the highest-counter non-persistent snapshots up to max_size total.
    if (max_size <= 0)
        return 0;
    std::unique_lock<std::shared_mutex> lock(m_mutex);
    const int n = static_cast<int>(m_structures.size());
    if (n <= max_size)
        return 0;

    std::vector<int> nonpersistent;
    int persistent_count = 0;
    for (int i = 0; i < n; ++i) {
        if (m_structures[i].persistent)
            ++persistent_count;
        else
            nonpersistent.push_back(i);
    }
    int keep_np = max_size - persistent_count; // how many non-persistent snapshots we may keep
    if (keep_np < 0)
        keep_np = 0;                            // persistent minima alone already fill (or exceed) the cap

    // Highest-counter (most-visited) non-persistent snapshots first; drop the tail.
    std::sort(nonpersistent.begin(), nonpersistent.end(),
        [this](int a, int b) { return m_structures[a].counter > m_structures[b].counter; });
    std::vector<char> drop(n, 0);
    for (int k = keep_np; k < static_cast<int>(nonpersistent.size()); ++k)
        drop[nonpersistent[k]] = 1;

    std::vector<BiasStructure> kept;
    kept.reserve(n);
    for (int i = 0; i < n; ++i)
        if (!drop[i])
            kept.push_back(m_structures[i]);
    const int removed = n - static_cast<int>(kept.size());
    m_structures = std::move(kept);
    for (int i = 0; i < static_cast<int>(m_structures.size()); ++i)
        m_structures[i].index = i; // re-index (safe between cycles; each MD run re-Initialises)
    m_global_count.store(static_cast<int>(m_structures.size()), std::memory_order_release);
    return removed;
}

void SharedBiasPool::setPermutations(const std::vector<std::vector<int>>& permutations)
{
    std::unique_lock<std::shared_mutex> lock(m_mutex);
    m_permutations = permutations;
}

std::vector<std::vector<int>> SharedBiasPool::permutations() const
{
    std::shared_lock<std::shared_mutex> lock(m_mutex);
    return m_permutations;
}

void SharedBiasPool::setWeights(const std::vector<double>& weights)
{
    std::unique_lock<std::shared_mutex> lock(m_mutex);
    m_weights = weights;
}

std::vector<double> SharedBiasPool::weights() const
{
    std::shared_lock<std::shared_mutex> lock(m_mutex);
    return m_weights;
}

void SharedBiasPool::clear()
{
    std::unique_lock<std::shared_mutex> lock(m_mutex);
    m_structures.clear();
    m_global_count.store(0, std::memory_order_release);
}

nlohmann::json SharedBiasPool::serializeMetadata() const
{
    std::shared_lock<std::shared_mutex> lock(m_mutex);
    nlohmann::json bias_array = nlohmann::json::array();
    for (const auto& bs : m_structures) {
        nlohmann::json entry;
        entry["time"] = bs.time;
        entry["rmsd_reference"] = bs.rmsd_reference;
        entry["energy"] = bs.energy;
        entry["factor"] = bs.factor;
        entry["index"] = bs.index;
        entry["counter"] = bs.counter;
        entry["temperature"] = bs.temperature;
        entry["persistent"] = bs.persistent;
        bias_array.push_back(entry);
    }
    return bias_array;
}

std::string SharedBiasPool::serializeGeometry() const
{
    std::shared_lock<std::shared_mutex> lock(m_mutex);
    std::ostringstream oss;
    for (const auto& bs : m_structures) {
        // Write each bias structure geometry as an XYZ frame
        Molecule mol;
        // Reconstruct molecule from geometry -- we need the reference atom types
        // This is handled by storing the full Molecule separately
        // For now, we store as raw coordinates
        int natoms = bs.geometry.rows();
        oss << natoms << "\n";
        oss << "Bias structure " << bs.index << " time=" << bs.time << "\n";
        for (int i = 0; i < natoms; ++i) {
            // Placeholder: coordinates only. Atom types must come from reference.
            oss << "X " << bs.geometry(i, 0) << " " << bs.geometry(i, 1) << " " << bs.geometry(i, 2) << "\n";
        }
    }
    return oss.str();
}

void SharedBiasPool::deserializeMetadata(const nlohmann::json& metadata)
{
    std::unique_lock<std::shared_mutex> lock(m_mutex);
    m_structures.clear();
    for (const auto& entry : metadata) {
        BiasStructure bs;
        bs.time = entry.value("time", 0.0);
        bs.rmsd_reference = entry.value("rmsd_reference", 0.0);
        bs.energy = entry.value("energy", 0.0);
        bs.factor = entry.value("factor", 1.0);
        bs.index = entry.value("index", 0);
        bs.counter = entry.value("counter", 0);
        bs.temperature = entry.value("temperature", 0.0);
        bs.persistent = entry.value("persistent", false);
        m_structures.push_back(std::move(bs));
    }
    m_global_count.store(static_cast<int>(m_structures.size()), std::memory_order_release);
}

void SharedBiasPool::deserializeGeometry(const std::string& xyz_data)
{
    // Geometries are restored by the caller who has the reference molecule
    // and can properly reconstruct atom types. This method is a placeholder
    // for the full serialization path.
    // For cross-temperature propagation, geometries are passed via BiasStructure
    // objects directly, not through string serialization.
}

void SharedBiasPool::restoreStructures(const std::vector<BiasStructure>& structures)
{
    // Claude Generated (Jun 2026): one-shot full-state restore for ConfSearch restart.
    // The caller has already rebuilt complete BiasStructure objects (geometry + counter +
    // energy + index + persistent flag) from the checkpoint, so we just take them verbatim.
    std::unique_lock<std::shared_mutex> lock(m_mutex);
    m_structures = structures;
    m_global_count.store(static_cast<int>(m_structures.size()), std::memory_order_release);
}