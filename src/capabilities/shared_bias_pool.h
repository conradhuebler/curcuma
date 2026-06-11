/*
 * <Shared Bias Pool for Parallel RMSD-MTD Conformational Search>
 * Copyright (C) 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Thread-safe container for sharing bias structures across parallel
 * MD workers during ConfSearch. Uses std::shared_mutex for read-heavy
 * access: every MD step reads all structures (shared lock), but new
 * structures are deposited infrequently (unique lock).
 *
 * Claude Generated (Apr 2026)
 */

#pragma once

#include <shared_mutex>
#include <atomic>
#include <vector>

#include "src/capabilities/simplemd.h"  // BiasStructure
#include "json.hpp"

class SharedBiasPool {
public:
    explicit SharedBiasPool() = default;
    ~SharedBiasPool() = default;

    // Non-copyable, non-movable (contains mutex)
    SharedBiasPool(const SharedBiasPool&) = delete;
    SharedBiasPool& operator=(const SharedBiasPool&) = delete;
    SharedBiasPool(SharedBiasPool&&) = delete;
    SharedBiasPool& operator=(SharedBiasPool&&) = delete;

    /** Snapshot all current bias structures (shared_lock, O(N) copy).
     *  Callers use this for bias evaluation; releases the lock immediately.
     *  The returned vector is a deep copy -- safe to use without locking. */
    std::vector<BiasStructure> snapshot() const;

    /** Atomic read of the current global structure count.
     *  Uses memory_order_relaxed -- sufficient for count checks
     *  in the deposition criterion. */
    int biasStructureCount() const;

    /** Deposit a new bias structure (unique_lock, exclusive access).
     *  Called when a SimpleMD instance decides to add a new Gaussian.
     *  Returns the assigned global index. */
    int depositBiasStructure(const BiasStructure& structure);

    /** Batch deposit multiple structures at once (unique_lock).
     *  Used for initial seeding and cross-temperature propagation.
     *  Returns the index of the first deposited structure. */
    int depositBatch(const std::vector<BiasStructure>& structures);

    /** Atomically register a visit to the given structures: counter++ (raises the
     *  exploration hill height W = k*counter) and factor += weight (the opt-in
     *  well-tempered OUTPUT weight; pass 0 when well-tempering is off, so the force
     *  and deposition are never affected). One unique_lock per batch -- called at most
     *  once per MD step, so it does not defeat the read-heavy snapshot() design.
     *  updates: (structure index, well-tempered weight increment) pairs. */
    void registerVisits(const std::vector<std::pair<int, double>>& updates);

    /** Prune structures whose counter is below threshold.
     *  Called between temperature cycles (no concurrent access).
     *  Removes rarely-visited regions to keep pool size manageable. */
    void pruneByCounter(int min_counter);

    /** Remove all non-persistent (raw MD snapshot) structures, keeping only
     *  the optimised minima (persistent=true). Called after feeding back
     *  optimised minima so the raw snapshots are replaced by their converged
     *  geometries and the pool stays clean for the next MD cycle. */
    void pruneNonPersistent();

    /** Claude Generated (Jun 2026): set the symmetry/atom-permutation set (full-atom reorder
     *  rules discovered by ConfScan). Set between temperature cycles; read every MD step.
     *  Empty (default) -> the RMSD-MTD bias uses the identity only (unpermuted, bit-identical). */
    void setPermutations(const std::vector<std::vector<int>>& permutations);

    /** Snapshot of the current permutation set (shared_lock, deep copy). The MTD sums a
     *  Gaussian over each image (smooth), so a discontinuous hard min is never formed. */
    std::vector<std::vector<int>> permutations() const;

    /** Claude Generated (Jun 2026): per-atom flexibility weights (1/sigma^2) for the RMSF-weighted
     *  RMSD-MTD bias (Phase C "weighted"). Set between cycles; read every MD step. Empty (default)
     *  -> uniform weights -> the standard unweighted best-fit RMSD (bit-identical). */
    void setWeights(const std::vector<double>& weights);

    /** Snapshot of the current per-atom weights (shared_lock, deep copy). */
    std::vector<double> weights() const;

    /** Remove all structures. Used when starting fresh
     *  or when resetting between temperature cycles. */
    void clear();

    /** Claude Generated (Jun 2026): full-state restore for ConfSearch restart.
     *  Replaces the pool contents with the given structures (metadata AND geometry),
     *  preserving each BiasStructure::index and counter exactly. Used when resuming a
     *  ConfSearch run from a checkpoint. Unlike deserializeMetadata()/deserializeGeometry()
     *  (which split the two and leave geometry as a placeholder), this restores ready-to-use
     *  bias structures in one call. */
    void restoreStructures(const std::vector<BiasStructure>& structures);

    /** Serialize metadata (counter, energy, factor, index, temperature)
     *  to JSON. Matches existing restart format for backward compatibility. */
    nlohmann::json serializeMetadata() const;

    /** Serialize geometries as XYZ string (multi-frame). */
    std::string serializeGeometry() const;

    /** Restore from JSON metadata. Geometries must be loaded separately
     *  via deserializeGeometry(). */
    void deserializeMetadata(const nlohmann::json& metadata);

    /** Restore geometries from XYZ string. Must match the order
     *  from serializeGeometry(). */
    void deserializeGeometry(const std::string& xyz_data);

private:
    mutable std::shared_mutex m_mutex;
    std::vector<BiasStructure> m_structures;
    std::atomic<int> m_global_count{0};
    std::vector<std::vector<int>> m_permutations;  // Claude Generated (Jun 2026): symmetry reorder rules (full-atom)
    std::vector<double> m_weights;                 // Claude Generated (Jun 2026): per-atom RMSF weights (empty = uniform)
};