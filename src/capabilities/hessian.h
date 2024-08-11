/*
 * < General Calculator for the Hessian>
 * Copyright (C) 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/core/global.h"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include "src/core/energycalculator.h"

#include "src/capabilities/curcumamethod.h"

#include <functional>

#pragma once

/* some of the stuff might be important later, but are not used now */

static const json HessianJson = {
    { "hess_read_file", "hessian" },
    { "hess_read_xyz", "none" },
    { "hess_calc", true },
    { "hess_read", false },
    { "hess_write_file", "none" },
    { "freq_scale", 1 },
    { "thermo", 298.15 },
    { "freq_cutoff", 50 },
    { "hess", 1 },
    { "method", "uff" },
    { "threads", 1 }
};

class HessianThread : public CxxThread {
public:
    HessianThread(const json& controller, int i, int j, int xi, int xj, bool fullnumerical = true);
    ~HessianThread();

    void setMolecule(const Molecule& molecule);
    void setParameter(const json& parameter) { m_parameter = parameter; }
    int execute() override;

    int I() const { return m_i; }
    int J() const { return m_j; }
    int XI() const { return m_xi; }
    int XJ() const { return m_xj; }
    double DD() const { return m_dd; }
    Matrix Gradient() const { return m_gradient; }

private:
    void Numerical();
    void Seminumerical();
    std::function<void(void)> m_schema;

    // EnergyCalculator* m_calculator;
    std::string m_method;
    json m_controller, m_parameter;
    Molecule m_molecule;
    Matrix m_gradient;
    Geometry m_geom_ip_jp, m_geom_im_jp, m_geom_ip_jm, m_geom_im_jm;
    int m_i, m_j, m_xi, m_xj;
    bool m_fullnumerical = true;
    double m_dd = 0;
    double m_d = 5e-3;
};

class Hessian : public CurcumaMethod {
public:
    Hessian(const std::string& method, const json& controller, bool silent = true);
    Hessian(const json& controller, bool silent = true);

    void setMolecule(const Molecule& molecule);

    void CalculateHessian(int type = 1)
    {
        m_hess_calc = type;
        start();
    }
    void PrintVibrations() { PrintVibrationsPlain(m_frequencies); }

    void PrintVibrations(Vector& eigenvalues, const Vector& projected);
    void PrintVibrationsPlain(const Vector& eigenvalues);

    void start() override;

    Matrix getHessian() const { return m_hessian; }
    void setParameter(const json& parameter) { m_parameter = parameter; }
    Vector Frequencies() { return m_frequencies; }

private:
    /* Lets have this for all modules */
    nlohmann::json WriteRestartInformation() override { return json(); }

    /* Lets have this for all modules */
    bool LoadRestartInformation() override { return true; }

    StringList MethodName() const override { return { std::string("Hessian") }; }

    /* Lets have all methods read the input/control file */
    void ReadControlFile() override {}

    /* Read Controller has to be implemented for all */
    void LoadControlJson() override;

    void CalculateHessianNumerical();
    void CalculateHessianSemiNumerical();
    void FiniteDiffHess();
    std::function<double(double)> m_scale_functions;

    void LoadMolecule(const std::string& file);
    void LoadHessian(const std::string& file);

    Matrix ProjectHessian(const Matrix& hessian);

    Vector ConvertHessian(Matrix& hessian);

    Matrix m_eigen_geometry, m_eigen_gradient, m_hessian;
    Vector m_frequencies;
    Molecule m_molecule;
    std::string m_method;
    json m_controller, m_parameter;
    std::vector<int> m_atoms_i, m_atoms_j;
    int m_threads = 1;
    int m_atom_count = 0;
    double m_freq_scale = 1, m_thermo = 298.5, m_freq_cutoff = 50;
    bool m_hess_calc = true, m_hess_write = false, m_hess_read = false;
    int m_hess = 1;
    std::string m_read_file = "none", m_write_file = "none", m_read_xyz = "none";
};
