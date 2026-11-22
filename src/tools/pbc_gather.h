/*
 * <Thread-count independent gather form of symmetric pair loops>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Claude Generated (Nov 2026, CPU threading pass of WP-PERIODIC-XTB, stage 2).
 *
 * The serial periodic pair loops of the native xTB and the D3/D4 helpers have the form
 *
 *     for i in 0..N-1:  for j in 0..i (or i-1):  for every image T of (i, j):
 *         out[i] += c(i, j, T);   out[j] += c'(i, j, T)          (c' = c or -c, j != i only)
 *
 * Atom a receives its contributions in this serial order: first the pairs (a, j) of its own row, j ascending,
 * then the pairs (i, a) of the later rows i > a, i ascending (rows i < a never touch a). gatherSymmetricPairs()
 * visits, for every atom a, exactly these pairs in exactly this order, so a per-atom sum built from them has the
 * bits of the serial loop for every thread count; the atoms are spread over the threads in contiguous ranges (each
 * written by one thread only). The price is that every pair is evaluated twice (once per member).
 *
 * Bitwise trap (measured): with -march=native GCC fuses a product into its addition (FMA) only when both are in
 * one basic block. A visit that adds "first ? +g : -g" through a branch changed the last bit (1e-19 in the D3
 * gradient). Keep the two update statements of the serial loop verbatim and point the end that is not atom a at a
 * scratch value: double* gi = first ? &out[a] : &scr; double* gj = first ? &scr : &out[a]; gi += g; gj -= g.
 *
 * Scalars (an energy) and the 3x3 virial are one running sum over all pairs in the serial loop. Callers record
 * the factors of every term during the visit with first == true (pair (a, j) of row a, in serial order) in a
 * per-row buffer (RowRecords) and replay the buffers in row order on one thread with the statement of the serial
 * loop: the same additions in the same order, so the same bits for every thread count. The factors are stored,
 * not their product, because the compiler may fuse the last multiplication with the addition (FMA).
 */

#pragma once

#include <array>
#include <functional>
#include <vector>

namespace PBCGather {

/// Runs worker(tid, nthreads) on nthreads threads (the caller's pool); nullptr = serial worker(0, 1).
using StripeRunner = std::function<void(const std::function<void(int, int)>&)>;

inline void runStripes(const StripeRunner* run, const std::function<void(int, int)>& worker)
{
    if (run && *run)
        (*run)(worker);
    else
        worker(0, 1);
}

/**
 * @brief Visit, per atom a, every pair of the serial symmetric loop that touches a, in serial order.
 * @param n          number of atoms
 * @param with_self  the serial row runs j = 0..i (true) or j = 0..i-1 (false)
 * @param visit      visit(a, i, j, first): pair (i, j) of row i; first == true means a == i (own row,
 *                   includes j == i when with_self), first == false means a == j < i
 */
template <class Visit>
void gatherSymmetricPairs(const StripeRunner* run, int n, bool with_self, Visit&& visit)
{
    // contiguous atom ranges per thread (the per-atom cost is uniform, every atom visits all n pairs): neighbouring
    // output rows stay in one thread, no false sharing on the per-atom accumulators
    runStripes(run, [&](int tid, int nth) {
        const int a0 = static_cast<int>((static_cast<long long>(n) * tid) / nth);
        const int a1 = static_cast<int>((static_cast<long long>(n) * (tid + 1)) / nth);
        for (int a = a0; a < a1; ++a) {
            const int jmax = with_self ? a : a - 1;
            for (int j = 0; j <= jmax; ++j)
                visit(a, a, j, true);
            for (int i = a + 1; i < n; ++i)
                visit(a, i, a, false);
        }
    });
}

/// Per-row record buffers for the deferred serial reduction of energies and virials (see above).
template <int K>
struct RowRecords {
    std::vector<std::vector<std::array<double, K>>> rows;
    explicit RowRecords(int n, bool active = true) : rows(active ? n : 0) {}
    bool active() const { return !rows.empty(); }
    void push(int row, const std::array<double, K>& r) { rows[row].push_back(r); }
    /// replay(r) for every record, rows in order, records in insertion order (one thread)
    template <class F>
    void replay(F&& f) const
    {
        for (const auto& row : rows)
            for (const auto& r : row) f(r);
    }
};

} // namespace PBCGather
