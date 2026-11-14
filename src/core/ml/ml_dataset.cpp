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

#include "ml_dataset.h"

#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>
#include <set>

namespace curcuma {
namespace ml {

// =================================================================================
// MLStructure
// =================================================================================

namespace {

    /*! \brief Flatten a natoms x 3 matrix into a row-major vector for JSON storage. */
    std::vector<double> matrixToVector(const Matrix& m)
    {
        std::vector<double> out;
        out.reserve(static_cast<std::size_t>(m.rows() * m.cols()));
        for (int i = 0; i < m.rows(); ++i)
            for (int j = 0; j < m.cols(); ++j)
                out.push_back(m(i, j));
        return out;
    }

    Matrix vectorToMatrix(const std::vector<double>& v, int rows, int cols)
    {
        Matrix m = Matrix::Zero(rows, cols);
        if (static_cast<int>(v.size()) != rows * cols)
            return m;
        int k = 0;
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                m(i, j) = v[k++];
        return m;
    }

} // namespace

json MLStructure::toJson() const
{
    json j;
    j["elements"] = elements;
    j["geometry"] = matrixToVector(geometry);
    j["charge"] = charge;
    j["spin"] = spin;
    j["energy_base"] = energy_base;
    j["energy_ref"] = energy_ref;
    j["has_gradients"] = has_gradients;
    j["cycle"] = cycle;
    if (has_gradients) {
        j["gradient_base"] = matrixToVector(gradient_base);
        j["gradient_ref"] = matrixToVector(gradient_ref);
    }
    return j;
}

MLStructure MLStructure::fromJson(const json& j)
{
    MLStructure s;
    s.elements = j.value("elements", std::vector<int>{});
    const int n = static_cast<int>(s.elements.size());
    s.geometry = vectorToMatrix(j.value("geometry", std::vector<double>{}), n, 3);
    s.charge = j.value("charge", 0);
    s.spin = j.value("spin", 0.0);
    s.energy_base = j.value("energy_base", 0.0);
    s.energy_ref = j.value("energy_ref", 0.0);
    s.has_gradients = j.value("has_gradients", false);
    s.cycle = j.value("cycle", -1);
    if (s.has_gradients) {
        s.gradient_base = vectorToMatrix(j.value("gradient_base", std::vector<double>{}), n, 3);
        s.gradient_ref = vectorToMatrix(j.value("gradient_ref", std::vector<double>{}), n, 3);
    }
    return s;
}

// =================================================================================
// MLDataset
// =================================================================================

std::vector<int> MLDataset::species() const
{
    std::set<int> s;
    for (const auto& st : m_structures)
        for (int z : st.elements)
            s.insert(z);
    return std::vector<int>(s.begin(), s.end());
}

bool MLDataset::hasGradients() const
{
    if (m_structures.empty())
        return false;
    for (const auto& s : m_structures)
        if (!s.has_gradients)
            return false;
    return true;
}

void MLDataset::split(double test_fraction, unsigned seed, MLDataset& train, MLDataset& test) const
{
    train.clear();
    test.clear();

    std::vector<std::size_t> order(m_structures.size());
    std::iota(order.begin(), order.end(), 0);
    std::mt19937 rng(seed);
    std::shuffle(order.begin(), order.end(), rng);

    const std::size_t n_test
        = static_cast<std::size_t>(std::round(test_fraction * static_cast<double>(m_structures.size())));
    for (std::size_t k = 0; k < order.size(); ++k) {
        if (k < n_test)
            test.add(m_structures[order[k]]);
        else
            train.add(m_structures[order[k]]);
    }
}

json MLDataset::toJson() const
{
    json j;
    j["format"] = "curcuma-deltaml-dataset";
    j["version"] = 1;
    j["n_structures"] = m_structures.size();
    json arr = json::array();
    for (const auto& s : m_structures)
        arr.push_back(s.toJson());
    j["structures"] = arr;
    return j;
}

bool MLDataset::fromJson(const json& j)
{
    m_structures.clear();
    if (!j.contains("structures") || !j["structures"].is_array()) {
        CurcumaLogger::error("MLDataset: JSON has no 'structures' array");
        return false;
    }
    for (const auto& e : j["structures"])
        m_structures.push_back(MLStructure::fromJson(e));
    return true;
}

bool MLDataset::save(const std::string& filename) const
{
    std::ofstream out(filename);
    if (!out.is_open()) {
        CurcumaLogger::error_fmt("MLDataset: cannot write '{}'", filename);
        return false;
    }
    out << toJson().dump(1);
    return true;
}

bool MLDataset::load(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in.is_open()) {
        CurcumaLogger::error_fmt("MLDataset: cannot read '{}'", filename);
        return false;
    }
    json j;
    try {
        in >> j;
    } catch (const std::exception& e) {
        CurcumaLogger::error_fmt("MLDataset: malformed JSON in '{}': {}", filename, e.what());
        return false;
    }
    return fromJson(j);
}

bool MLDataset::writeXYZ(const std::string& filename) const
{
    std::ofstream out(filename);
    if (!out.is_open()) {
        CurcumaLogger::error_fmt("MLDataset: cannot write '{}'", filename);
        return false;
    }
    for (const auto& s : m_structures) {
        out << s.AtomCount() << "\n";
        out << std::fixed << std::setprecision(8) << "E_base= " << s.energy_base
            << " E_ref= " << s.energy_ref << " delta= " << s.delta() << " cycle= " << s.cycle << "\n";
        for (int i = 0; i < s.AtomCount(); ++i) {
            const int z = s.elements[i];
            const std::string sym
                = (z > 0 && z < static_cast<int>(Elements::ElementAbbr.size()))
                ? Elements::ElementAbbr[z]
                : "X";
            out << std::left << std::setw(3) << sym << std::right << std::fixed << std::setprecision(8)
                << std::setw(16) << s.geometry(i, 0) << std::setw(16) << s.geometry(i, 1)
                << std::setw(16) << s.geometry(i, 2) << "\n";
        }
    }
    return true;
}

} // namespace ml
} // namespace curcuma
