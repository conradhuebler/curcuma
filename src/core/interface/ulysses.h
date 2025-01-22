/*
 * < C++ Ulysses Interface >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 */

#include "Eigen/Dense"

#include <vector>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Geometry;

class BSet;
class QCbasis;

class UlyssesObject {
public:
    UlyssesObject();
    ~UlyssesObject();

    void Calculate(bool gradient, bool verbose);
    void setMethod(const std::string& method);
    void setMolecule(const Geometry& geom, const std::vector<int>& atm, int chrge, int multpl, std::string pg);
    void UpdateGeometry(const Geometry& geom);
    double Energy() const { return m_energy; }
    Geometry Gradient() const { return m_gradient; }
    void setTele(double Tele) { m_Tele = Tele; }
    void setMaxIter(int maxiter) { m_SCFmaxiter = maxiter; }

private:
    BSet* m_bset;
    QCbasis* m_electron;
    std::string m_method;
    double m_Tele;
    int m_SCFmaxiter;
    double m_energy;
    Geometry m_gradient;
};
