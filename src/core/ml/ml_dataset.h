/*
 * < Labelled structure container for Delta-Machine-Learning >
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
 * Claude Generated (July 2026): dataset layer of the Delta-ML work package.
 */

#pragma once

#include "src/core/global.h"
#include "json.hpp"

#include <string>
#include <vector>

using json = nlohmann::json;

namespace curcuma {
namespace ml {

/**
 * @brief One labelled structure: geometry plus base and reference results.
 *
 * Unit convention follows curcuma throughout: geometries in Angstrom, energies in Hartree,
 * hence gradients in Eh/Angstrom (verified against EnergyCalculator::NumGrad, which steps
 * the Angstrom geometry directly).
 *
 * The learning target is the difference, never the absolute values:
 *   delta        = E_ref - E_base
 *   deltaGradient = G_ref - G_base
 */
struct MLStructure {
    std::vector<int> elements; ///< atomic numbers
    Geometry geometry;         ///< natoms x 3 [Angstrom]
    int charge = 0;
    double spin = 0.0;

    double energy_base = 0.0; ///< [Eh]
    double energy_ref = 0.0;  ///< [Eh]
    Matrix gradient_base;     ///< natoms x 3 [Eh/Angstrom]
    Matrix gradient_ref;      ///< natoms x 3 [Eh/Angstrom]
    bool has_gradients = false;

    int cycle = -1; ///< active-learning cycle that produced this point, -1 for offline data

    int AtomCount() const { return static_cast<int>(elements.size()); }
    double delta() const { return energy_ref - energy_base; }
    Matrix deltaGradient() const { return gradient_ref - gradient_base; }

    json toJson() const;
    static MLStructure fromJson(const json& j);
};

/**
 * @brief Collection of labelled structures with JSON persistence.
 *
 * Kept deliberately simple: the dataset is small by ML standards (hundreds of points),
 * because that is the entire point of delta-learning.
 */
class MLDataset {
public:
    void add(const MLStructure& s) { m_structures.push_back(s); }
    std::size_t size() const { return m_structures.size(); }
    bool empty() const { return m_structures.empty(); }
    void clear() { m_structures.clear(); }

    const MLStructure& operator[](std::size_t i) const { return m_structures[i]; }
    MLStructure& operator[](std::size_t i) { return m_structures[i]; }

    const std::vector<MLStructure>& structures() const { return m_structures; }

    /*! \brief Sorted, unique atomic numbers occurring anywhere in the set. */
    std::vector<int> species() const;

    /*! \brief True if every structure carries both gradients. */
    bool hasGradients() const;

    /**
     * @brief Deterministic train/test split.
     * @param test_fraction share of structures routed into @p test
     * @param seed          RNG seed, fixed so results are reproducible
     *
     * Random rather than by index: consecutive structures in a trajectory or conformer
     * file are correlated, so a contiguous tail split would flatter the model.
     */
    void split(double test_fraction, unsigned seed, MLDataset& train, MLDataset& test) const;

    bool save(const std::string& filename) const;
    bool load(const std::string& filename);

    /*! \brief Convenience dump of the geometries as multi-XYZ for visual inspection. */
    bool writeXYZ(const std::string& filename) const;

    json toJson() const;
    bool fromJson(const json& j);

private:
    std::vector<MLStructure> m_structures;
};

} // namespace ml
} // namespace curcuma
