/*
 * < Generic force field class for curcuma . >
 * Copyright (C) 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#pragma once

#include "src/core/global.h"

#include "hbonds.h"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#ifdef USE_D3
#include "src/core/dftd3interface.h"
#endif

#ifdef USE_D4
#include "src/core/dftd4interface.h"
#endif

#include "src/core/qmdff_par.h"
#include "src/core/uff_par.h"

#include <functional>
#include <set>
#include <vector>

#include <Eigen/Dense>

#include "json.hpp"
using json = nlohmann::json;

struct Bond {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0, k = 0, distance = 0;
    double fc = 0, exponent = 0, r0_ij = 0, r0_ik = 0;
};

struct Angle {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0, k = 0;
    double fc = 0, r0_ij = 0, r0_ik = 0, theta0_ijk = 0;
    double C0 = 0, C1 = 0, C2 = 0;
};

struct Dihedral {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0, k = 0, l = 0;
    double V = 0, n = 0, phi0 = 0;
};

struct Inversion {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0, k = 0, l = 0;
    double fc = 0, C0 = 0, C1 = 0, C2 = 0;
};

struct vdW {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0;
    double C_ij = 0, r0_ij = 0;
};

struct EQ {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0;
    double C_ij = 0, r0_ij = 0;
};
class ForceFieldThread : public CxxThread {

public:
    ForceFieldThread(int thread, int threads);
    virtual int execute() override;
    virtual int Type() const { return 1; }
    void addBond(const Bond& bonds);
    void addAngle(const Angle& angles);
    void addDihedral(const Dihedral& dihedrals);
    void addInversion(const Inversion& inversions);
    void addvdW(const vdW& vdWs);
    void addEQ(const EQ& EQs);

    inline void UpdateGeometry(const Matrix& geometry, bool gradient)
    {
        m_geometry = geometry;
        m_calculate_gradient = gradient;
        m_gradient = Eigen::MatrixXd::Zero(m_geometry.rows(), 3);
    }

    inline void setGeometry(const Matrix& geometry, bool gradient)
    {
        m_geometry = geometry;
        m_calculate_gradient = gradient;
        m_gradient = Eigen::MatrixXd::Zero(m_geometry.rows(), 3);
    }

    inline void setMethod(int method)
    {
        m_method = method;
    }
    double BondEnergy() { return m_bond_energy; }
    double AngleEnergy() { return m_angle_energy; }
    double DihedralEnergy() { return m_dihedral_energy; }
    double InversionEnergy() { return m_inversion_energy; }
    double VdWEnergy() { return m_vdw_energy; }
    double RepEnergy() { return m_rep_energy; }
    double EQEnergy() { return m_eq_energy; }

    Matrix Gradient() const { return m_gradient; }

private:
    void CalculateUFFBondContribution();
    void CalculateUFFAngleContribution();
    void CalculateUFFDihedralContribution();
    void CalculateUFFInversionContribution();
    void CalculateUFFvdWContribution();

    void CalculateQMDFFBondContribution();
    void CalculateQMDFFAngleContribution();
    void CalculateQMDFFDihedralContribution();

    // double HarmonicBondStretching();

    // double LJBondStretching();

    std::function<double()> CalculateBondContribution;
    std::function<double()> CalculateAngleContribution;
    std::function<double()> CalculateTorsionContribution;
    std::function<double()> CalculateInversionContribution;
    std::function<double()> CalculateVdWContribution;
    std::function<double()> CalculateEQContribution;
    std::function<double()> CalculateHBondContribution;

    std::vector<Bond> m_uff_bonds, m_qmdff_bonds;
    std::vector<Angle> m_uff_angles, m_qmdff_angles;
    std::vector<Dihedral> m_uff_dihedrals, m_qmdff_dihedrals;
    std::vector<Inversion> m_uff_inversions, m_qmdff_inversions;
    std::vector<vdW> m_uff_vdWs;
    std::vector<EQ> m_qmdff_EQs;

protected:
    Matrix m_geometry, m_gradient;
    double m_energy = 0, m_bond_energy = 0.0, m_angle_energy = 0.0, m_dihedral_energy = 0.0, m_inversion_energy = 0.0, m_vdw_energy = 0.0, m_rep_energy = 0.0, m_eq_energy = 0.0;

    double m_final_factor = 1;
    double m_bond_scaling = 1, m_angle_scaling = 1, m_dihedral_scaling = 1, m_inversion_scaling = 1, m_vdw_scaling = 1, m_rep_scaling = 1;
    double m_au = 1;
    double m_d = 1e-3;
    int m_calc_gradient = 1;
    int m_thread = 0, m_threads = 0, m_method = 1;
    bool m_calculate_gradient = true;
};

class D3Thread : public ForceFieldThread {

public:
    D3Thread(int thread, int threads);
    ~D3Thread();
    virtual int Type() const { return 2; }

    void setParamater(const json& parameter)
    {
#ifdef USE_D3
        m_d3->UpdateParametersD3(parameter);
#endif
    }

    void Initialise(const std::vector<int>& atom_types)
    {
        m_atom_types = atom_types;
#ifdef USE_D3
        m_d3->InitialiseMolecule(m_atom_types);
#endif
    }
    virtual int execute() override;

private:
#ifdef USE_D3
    DFTD3Interface* m_d3;
#endif
    std::vector<int> m_atom_types;
};

class H4Thread : public ForceFieldThread {

public:
    H4Thread(int thread, int threads);
    ~H4Thread();
    virtual int Type() const { return 3; }

    void setParamater(const json& parameter)
    {
        m_h4correction.set_OH_O(parameter["h4_oh_o"].get<double>());
        m_h4correction.set_OH_N(parameter["h4_oh_n"].get<double>());
        m_h4correction.set_NH_O(parameter["h4_nh_o"].get<double>());
        m_h4correction.set_NH_N(parameter["h4_nh_n"].get<double>());

        m_h4correction.set_WH_O(parameter["h4_wh_o"].get<double>());
        m_h4correction.set_NH4(parameter["h4_nh4"].get<double>());
        m_h4correction.set_COO(parameter["h4_coo"].get<double>());
        m_h4correction.set_HH_Rep_K(parameter["hh_rep_k"].get<double>());
        m_h4correction.set_HH_Rep_E(parameter["hh_rep_e"].get<double>());
        m_h4correction.set_HH_Rep_R0(parameter["hh_rep_r0"].get<double>());
    }

    void Initialise(const std::vector<int>& atom_types)
    {
        m_atom_types = atom_types;
        m_h4correction.allocate(m_atom_types.size());
    }
    virtual int execute() override;

private:
    hbonds4::H4Correction m_h4correction;
    std::vector<int> m_atom_types;
};
