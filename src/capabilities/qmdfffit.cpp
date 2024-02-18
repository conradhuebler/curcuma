/*
 * <QMDFF Hessian fit for parametrisation. >
 * Copyright (C) 2023 - 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/capabilities/optimiser/LevMarQMDFFFit.h"

#include "src/capabilities/hessian.h"

#include "src/core/energycalculator.h"
#include "src/core/qmdff_par.h"
#include "src/core/topology.h"
#include "src/core/uff_par.h"

#include "qmdfffit.h"

QMDFFFit::QMDFFFit(const json& controller, bool silent)
    : CurcumaMethod(QMDFFFitJson, controller, silent)
{
    UpdateController(controller);
}

void QMDFFFit::LoadControlJson()
{
    m_method = Json2KeyWord<std::string>(m_defaults, "method");
    m_threads = Json2KeyWord<int>(m_defaults, "threads");
}

void QMDFFFit::start()
{

    std::cout << "Parametrising QMDFF (see S. Grimmme, J. Chem. Theory Comput. 2014, 10, 10, 4497–4514 [10.1021/ct500573f]) for the original publication!" << std::endl;

    std::cout << "Starting with the hessian ..." << std::endl;
    std::cout << m_defaults << m_controller << std::endl;
    Hessian hessian(m_defaults);
    hessian.setMolecule(m_molecule);
    m_atom_types = m_molecule.Atoms();
    m_geometry = m_molecule.getGeometry();
    hessian.start();
    m_hessian = hessian.getHessian();
    // std::cout << m_hessian << std::endl;
    Initialise();
    json parameter;
    parameter["bonds"] = Bonds();
    parameter["angles"] = Angles();
    std::ofstream parameterfile_init("qmdff_init_param.json");
    parameterfile_init << parameter;
    json qmdff_init = MergeJson(m_controller, QMDFFFitJson);
    qmdff_init["threads"] = m_threads;
    qmdff_init["variable"] = false;
    qmdff_init["method"] = "qmdff";
    EnergyCalculator calculator("qmdff", qmdff_init);
    calculator.setMolecule(m_molecule);
    calculator.setParameter(parameter);
    calculator.CalculateEnergy(false, false);
    Hessian const_hessian("qmdff", qmdff_init);
    /*
    parameter["bonds"] = bonds;
    parameter["angles"] = angles;
    */
    const_hessian.setMolecule(m_molecule);
    const_hessian.setParameter(parameter);
    const_hessian.start();
    Matrix const_hessian_matrix = const_hessian.getHessian();
    for (int i = 0; i < const_hessian_matrix.cols(); ++i)
        for (int j = i; j < const_hessian_matrix.cols(); ++j) {
            if (std::isnan(const_hessian_matrix(i, j))) {
                const_hessian_matrix(i, j) = 0;
                const_hessian_matrix(j, i) = 0;
            }
        }
    // std::cout << const_hessian_matrix << std::endl;
    int counter = 1;
    for (int start = 0; start < 10 && counter; ++start) {
        parameter["bonds"] = Bonds();
        parameter["angles"] = Angles();

        m_fc_parameter.resize(m_qmdffbonds.size() + m_qmdffangle.size());
        int i = 0;
        for (const auto& b : m_qmdffbonds) {
            m_fc_parameter(i) = std::abs(b.kAB);
            ++i;
        }
        for (const auto& b : m_qmdffangle) {
            m_fc_parameter(i) = std::abs(b.kabc);
            ++i;
        }
        qmdff_init["variable"] = true;
        qmdff_init["const"] = false;

        Vector vec = OptimiseFC(m_molecule, m_hessian, const_hessian_matrix, m_fc_parameter, parameter, qmdff_init);
        std::cout << vec << std::endl;
        json bonds = parameter["bonds"];
        json angles = parameter["angles"];
        int index = 0;
        for (int i = 0; i < bonds.size(); ++i) {
            bonds[i]["kAB"] = vec(index);
            m_qmdffbonds[i].kAB = vec(index);
            index++;
        }
        for (int i = 0; i < angles.size(); ++i) {
            angles[i]["kabc"] = vec(index);
            m_qmdffangle[i].kabc = vec(index);
            index++;
        }
        json hc = HessianJson;
        hc["method"] = "qmdff";
        hc["threads"] = m_threads;

        auto cache = m_qmdffbonds;
        m_qmdffbonds.clear();
        counter = 0;
        for (auto c : cache) {
            if ((c.kAB > 0 && c.distance == 1) || c.distance == 0)
                m_qmdffbonds.push_back(c);
            else
                counter++;
        }
        auto cache2 = m_qmdffangle;
        m_qmdffangle.clear();
        for (auto c : cache2) {
            if (c.kabc > 0)
                m_qmdffangle.push_back(c);
            else
                counter++;
        }
        std::cout << "Skipping " << counter << " force constants " << std::endl;
        // std::cout << he2.getHessian() << std::endl;
    }

    std::cout << "Eigenvectors" << std::endl;
    std::cout << hessian.Frequencies().transpose() << std::endl;
    qmdff_init["variable"] = true;
    qmdff_init["const"] = true;

    Hessian he2(qmdff_init, false);

    parameter["bonds"] = Bonds();
    parameter["angles"] = Angles();
    he2.setMolecule(m_molecule);
    he2.setParameter(parameter);
    he2.start();
    std::cout << he2.Frequencies().transpose() << std::endl;

    std::cout << parameter << std::endl;
    std::ofstream parameterfile("qmdff_param.json");
    parameterfile << parameter;
}

bool QMDFFFit::Initialise()
{
    // std::cout << "Initialising QMDFF (see S. Grimmme, J. Chem. Theory Comput. 2014, 10, 10, 4497–4514 [10.1021/ct500573f]) for the original publication!" << std::endl;

    // m_uff_atom_types = std::vector<int>(m_atom_types.size(), 0);
    m_coordination = std::vector<int>(m_atom_types.size(), 0);
    std::vector<std::set<int>> ignored_vdw;
    // m_topo = Eigen::MatrixXd::Zero(m_atom_types.size(), m_atom_types.size());
    TContainer bonds, nonbonds, angles, dihedrals, inversions;
    m_scaling = 1.4;
    // m_gradient = Eigen::MatrixXd::Zero(m_atom_types.size(), 3);
    for (int i = 0; i < m_atom_types.size(); ++i) {
        m_stored_bonds.push_back(std::vector<int>());
        ignored_vdw.push_back(std::set<int>({ i }));
        for (int j = 0; j < m_atom_types.size() && m_stored_bonds[i].size() < CoordinationNumber[m_atom_types[i]]; ++j) {
            if (i == j)
                continue;
            double x_i = m_geometry(i, 0) * m_au;
            double x_j = m_geometry(j, 0) * m_au;

            double y_i = m_geometry(i, 1) * m_au;
            double y_j = m_geometry(j, 1) * m_au;

            double z_i = m_geometry(i, 2) * m_au;
            double z_j = m_geometry(j, 2) * m_au;

            double r_ij = sqrt((((x_i - x_j) * (x_i - x_j)) + ((y_i - y_j) * (y_i - y_j)) + ((z_i - z_j) * (z_i - z_j))));

            if (r_ij <= (Elements::CovalentRadius[m_atom_types[i]] + Elements::CovalentRadius[m_atom_types[j]]) * m_scaling * m_au) {
                if (bonds.insert({ std::min(i, j), std::max(i, j) })) {
                    m_coordination[i]++;
                    m_stored_bonds[i].push_back(j);
                    ignored_vdw[i].insert(j);
                }
                //   m_topo(i, j) = 1;
                //   m_topo(j, i) = 1;
            }
        }
    }

    if (m_rings)
        m_identified_rings = Topology::FindRings(m_stored_bonds, m_atom_types.size());

    bonds.clean();
    setBonds(bonds, ignored_vdw, angles, dihedrals, inversions);

    angles.clean();
    setAngles(angles, ignored_vdw);
    /*
    dihedrals.clean();
    setDihedrals(dihedrals);

    inversions.clean();
    setInversions(inversions);

    nonbonds.clean();
    setvdWs(ignored_vdw);

    */
    return true;
}

void QMDFFFit::setBonds(const TContainer& bonds, std::vector<std::set<int>>& ignored_vdw, TContainer& angels, TContainer& dihedrals, TContainer& inversions)
{
    for (const auto& bond : bonds.Storage()) {
        QMDFFBond b;

        b.a = bond[0];
        b.b = bond[1];
        double xa = (m_geometry)(b.a, 0) * m_au;
        double xb = (m_geometry)(b.b, 0) * m_au;

        double ya = (m_geometry)(b.a, 1) * m_au;
        double yb = (m_geometry)(b.b, 1) * m_au;

        double za = (m_geometry)(b.a, 2) * m_au;
        double zb = (m_geometry)(b.b, 2) * m_au;

        b.reAB = Topology::Distance(xa, xb, ya, yb, za, zb);
        b.kAB = 1;
        double kaA = ka(m_atom_types[b.a]);
        double kaB = ka(m_atom_types[b.b]);
        double dEN = Elements::PaulingEN[m_atom_types[b.a]] - Elements::PaulingEN[m_atom_types[b.b]];
        b.exponA = kaA * kaB + kEN * dEN * dEN;
        b.distance = 0;
        m_qmdffbonds.push_back(b);

        int i = bond[0];
        int j = bond[1];

        std::vector<int> k_bodies;
        for (auto t : m_stored_bonds[i]) {
            k_bodies.push_back(t);

            if (t == j)
                continue;
            angels.insert({ i, std::min(t, j), std::max(j, t) });
            ignored_vdw[i].insert(t);
        }

        std::vector<int> l_bodies;
        for (auto t : m_stored_bonds[j]) {
            l_bodies.push_back(t);

            if (t == i)
                continue;
            angels.insert({ std::min(i, t), j, std::max(t, i) });
            ignored_vdw[j].insert(t);
        }

        for (int k : k_bodies) {
            for (int l : l_bodies) {
                if (k == i || k == j || k == l || i == j || i == l || j == l)
                    continue;
                dihedrals.insert({ k, i, j, l });

                ignored_vdw[i].insert(k);
                ignored_vdw[i].insert(l);
                ignored_vdw[j].insert(k);
                ignored_vdw[j].insert(l);
                ignored_vdw[k].insert(l);
                ignored_vdw[l].insert(k);
            }
        }
        if (m_stored_bonds[i].size() == 3) {
            inversions.insert({ i, m_stored_bonds[i][0], m_stored_bonds[i][1], m_stored_bonds[i][2] });
        }
        if (m_stored_bonds[j].size() == 3) {
            inversions.insert({ j, m_stored_bonds[j][0], m_stored_bonds[j][1], m_stored_bonds[j][2] });
        }
    }
}

void QMDFFFit::setAngles(const TContainer& angles, const std::vector<std::set<int>>& ignored_vdw)
{
    for (const auto& angle : angles.Storage()) {

        if (angle[0] == angle[1] || angle[0] == angle[2] || angle[1] == angle[2])
            continue;

        double xa = (m_geometry)(angle[0], 0) * m_au;
        double xb = (m_geometry)(angle[1], 0) * m_au;
        double xc = (m_geometry)(angle[2], 0) * m_au;

        double ya = (m_geometry)(angle[0], 1) * m_au;
        double yb = (m_geometry)(angle[1], 1) * m_au;
        double yc = (m_geometry)(angle[2], 1) * m_au;

        double za = (m_geometry)(angle[0], 2) * m_au;
        double zb = (m_geometry)(angle[1], 2) * m_au;
        double zc = (m_geometry)(angle[2], 2) * m_au;

        double rab = Topology::Distance(xa, xb, ya, yb, za, zb);
        double rac = Topology::Distance(xa, xc, ya, yc, za, zc);
        double rbc = Topology::Distance(xb, xc, yb, yc, zb, zc);

        QMDFFAngle a;

        QMDFFBond b;
        if (rab < rac && rbc < rac) {
            a.a = angle[1];
            a.b = angle[0];
            a.c = angle[2];

            b.a = angle[1];
            b.b = angle[2];
        } else if (rac < rab && rbc < rab) {
            a.a = angle[1];
            a.b = angle[0];
            a.c = angle[2];

            b.a = angle[1];
            b.b = angle[2];
        } else if (rac < rbc && rab < rbc) {
            a.a = angle[0];
            a.b = angle[1];
            a.c = angle[2];

            b.a = angle[0];
            b.b = angle[2];
        }
        b.distance = 1;
#pragma message("which distance should be used")

        b.reAB = Topology::Distance(xa, xb, ya, yb, za, zb) + Topology::Distance(xa, xc, ya, yc, za, zc);
        b.kAB = 0.10;
        double kaA = ka(m_atom_types[b.a]);
        double kaB = ka(m_atom_types[b.b]);
        double dEN = Elements::PaulingEN[m_atom_types[b.a]] - Elements::PaulingEN[m_atom_types[b.b]];
        b.exponA = ka13 + kb13 * kaA * kaB;
        b.distance = 1;
        // m_qmdffbonds.push_back(b);

        a.reAB = Topology::Distance(xa, xb, ya, yb, za, zb);
        a.reAC = Topology::Distance(xa, xc, ya, yc, za, zc);
        auto atom_0 = m_geometry.row(a.a);
        auto atom_1 = m_geometry.row(a.b);
        auto atom_2 = m_geometry.row(a.c);

        auto vec_1 = atom_0 - atom_1;
        auto vec_2 = atom_0 - atom_2;

        a.thetae = acos((vec_1.dot(vec_2) / (sqrt(vec_1.dot(vec_1)) * sqrt(vec_2.dot(vec_2))))) * 180 / pi;
        a.kabc = 0.1;

        m_qmdffangle.push_back(a);
    }
}

json QMDFFFit::Bonds() const
{
    json bonds;
    for (int i = 0; i < m_qmdffbonds.size(); ++i) {
        json bond;
        bond["a"] = m_qmdffbonds[i].a;
        bond["b"] = m_qmdffbonds[i].b;
        bond["reAB"] = m_qmdffbonds[i].reAB;
        bond["kAB"] = m_qmdffbonds[i].kAB;
        bond["exponA"] = m_qmdffbonds[i].exponA;
        bond["distance"] = m_qmdffbonds[i].distance;
        bonds[i] = bond;
    }
    return bonds;
}

json QMDFFFit::Angles() const
{
    json angles;

    for (int i = 0; i < m_qmdffangle.size(); ++i) {
        json angle;
        angle["a"] = m_qmdffangle[i].a;
        angle["b"] = m_qmdffangle[i].b;
        angle["c"] = m_qmdffangle[i].c;

        angle["kabc"] = m_qmdffangle[i].kabc;
        angle["thetae"] = m_qmdffangle[i].thetae;
        angle["reAB"] = m_qmdffangle[i].reAB;
        angle["reAC"] = m_qmdffangle[i].reAC;

        angles[i] = angle;
    }
    return angles;
}
