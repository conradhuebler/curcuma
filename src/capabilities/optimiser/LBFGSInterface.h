/*
 * <Geometrie optimisation using external LBFGS and xtb. >
 * Copyright (C) 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/NonLinearOptimization>

#include <iomanip>
#include <iostream>

#include "src/capabilities/rmsd.h"

#include "src/core/elements.h"
#include "src/core/global.h"
#include "src/core/molecule.h"
#include "src/core/xtbinterface.h"

#include "src/tools/geometry.h"

#include "json.hpp"
#include <LBFGS.h>
#include <LBFGSB.h>

using json = nlohmann::json;

const json OptJson{
    { "writeXYZ", false },
    { "printOutput", true },
    { "dE", 0.75 },
    { "dRMSD", 0.01 },
    { "GFN", 2 },
    { "InnerLoop", 20 },
    { "OuterLoop", 100 },
    { "LBFGS_eps", 1e-5 }
};

using Eigen::VectorXd;
using namespace LBFGSpp;

class LBFGSInterface {
public:
    LBFGSInterface(int n_)
        : n(n_)
    {
    }
    double operator()(const VectorXd& x, VectorXd& grad)
    {
        double fx = 0.0;
        int attyp[m_atoms];
        double charge = 0;
        std::vector<int> atoms = m_molecule->Atoms();
        double coord[3 * m_atoms];
        double gradient[3 * m_atoms];
        double dist_gradient[3 * m_atoms];

        for (int i = 0; i < m_atoms; ++i) {
            attyp[i] = m_molecule->Atoms()[i];
            coord[3 * i + 0] = x(3 * i + 0) / au;
            coord[3 * i + 1] = x(3 * i + 1) / au;
            coord[3 * i + 2] = x(3 * i + 2) / au;
        }
        m_interface->UpdateMolecule(coord);
        //try {
        fx = m_interface->GFNCalculation(m_method, gradient);
        //} catch (const std::logic_error& error) {
        // }

        for (int i = 0; i < m_atoms; ++i) {
            grad[3 * i + 0] = gradient[3 * i + 0];
            grad[3 * i + 1] = gradient[3 * i + 1];
            grad[3 * i + 2] = gradient[3 * i + 2];
        }
        m_energy = fx;
        m_parameter = x;
        return fx;
    }

    double LastEnergy() const { return m_energy; }

    double m_energy = 0, m_last_change = 0, m_last_rmsd = 0;
    Vector Parameter() const { return m_parameter; }
    void setMolecule(const Molecule* molecule)
    {
        m_molecule = molecule;
        m_atoms = m_molecule->AtomCount();
    }
    void setInterface(XTBInterface* interface) { m_interface = interface; }
    void setMethod(int method) { m_method = method; }

private:
    int m_iter = 0;
    int m_atoms = 0;
    int n;
    int m_method = 2;
    XTBInterface* m_interface;
    Vector m_parameter;
    const Molecule* m_molecule;
};

inline Molecule OptimiseGeometry(const Molecule* host, const json& controller)
{
    //PrintController(controller);
    bool writeXYZ = Json2KeyWord<bool>(controller, "writeXYZ");
    bool printOutput = Json2KeyWord<bool>(controller, "printOutput");
    double dE = Json2KeyWord<double>(controller, "dE");
    double dRMSD = Json2KeyWord<double>(controller, "dRMSD");
    double LBFGS_eps = Json2KeyWord<double>(controller, "LBFGS_eps");
    int method = Json2KeyWord<int>(controller, "GFN");
    int InnerLoop = Json2KeyWord<int>(controller, "InnerLoop");
    int OuterLoop = Json2KeyWord<int>(controller, "OuterLoop");

    Geometry geometry = host->getGeometry();
    Molecule tmp(host);
    Molecule h(host);
    Vector parameter(3 * host->AtomCount());

    for (int i = 0; i < host->AtomCount(); ++i) {
        parameter(3 * i) = geometry(i, 0);
        parameter(3 * i + 1) = geometry(i, 1);
        parameter(3 * i + 2) = geometry(i, 2);
    }

    XTBInterface interface;
    interface.InitialiseMolecule(host);

    double final_energy = interface.GFNCalculation(method);

    LBFGSParam<double> param;
    param.epsilon = LBFGS_eps;
    param.max_iterations = InnerLoop;

    LBFGSSolver<double> solver(param);
    LBFGSInterface fun(3 * host->AtomCount());
    fun.setMolecule(host);
    fun.setInterface(&interface);
    fun.setMethod(method);
    double fx;

    json RMSDJsonControl = {
        { "reorder", false },
        { "check", false },
        { "heavy", false },
        { "fragment", -1 },
        { "fragment_reference", -1 },
        { "fragment_target", -1 },
        { "init", -1 },
        { "pt", 0 },
        { "silent", true },
        { "storage", 1.0 },
        { "method", "incr" },
        { "noreorder", true },
        { "threads", 1 }
    };

    RMSDDriver* driver = new RMSDDriver(RMSDJsonControl);

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now(), end;
    if (printOutput) {
        std::cout << "Step\tCurrent Energy [Eh]\tEnergy Change\tRMSD Change\tt [s]" << std::endl;
    }
    int atoms_count = host->AtomCount();
    for (int outer = 0; outer < OuterLoop; ++outer) {
        int niter = 0;
        try {
            niter = solver.minimize(fun, parameter, fx);
        } catch (const std::logic_error& error) {
        }
        parameter = fun.Parameter();

        for (int i = 0; i < atoms_count; ++i) {
            geometry(i, 0) = parameter(3 * i);
            geometry(i, 1) = parameter(3 * i + 1);
            geometry(i, 2) = parameter(3 * i + 2);
        }
        h.setGeometry(geometry);

        driver->setReference(tmp);
        driver->setTarget(h);
        driver->start();
        if (printOutput) {
            end = std::chrono::system_clock::now();
            std::cout << outer << "\t" << std::setprecision(9) << fun.m_energy << "\t\t" << std::setprecision(5) << (fun.m_energy - final_energy) * 2625.5 << "\t\t" << std::setprecision(6) << driver->RMSD() << "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0 << std::endl;
            start = std::chrono::system_clock::now();
        }
        final_energy = fun.m_energy;
        tmp = h;
        if (writeXYZ) {
            h.setEnergy(final_energy);
            h.appendXYZFile("curcuma_optim.xyz");
        }
        if (((fun.m_energy - final_energy) * 2625.5 < dE && driver->RMSD() < dRMSD))
            break;
    }

    for (int i = 0; i < host->AtomCount(); ++i) {
        geometry(i, 0) = parameter(3 * i);
        geometry(i, 1) = parameter(3 * i + 1);
        geometry(i, 2) = parameter(3 * i + 2);
    }
    h.setEnergy(final_energy);
    h.setGeometry(geometry);
    return h;
}

inline void OptimiseGeometryThreaded(const Molecule* host, std::string* result_string, Molecule* result_molecule, const json& controller)
{
    //PrintController(controller);
    bool writeXYZ = Json2KeyWord<bool>(controller, "writeXYZ");
    bool printOutput = Json2KeyWord<bool>(controller, "printOutput");
    double dE = Json2KeyWord<double>(controller, "dE");
    double dRMSD = Json2KeyWord<double>(controller, "dRMSD");
    int method = Json2KeyWord<int>(controller, "GFN");
    int InnerLoop = Json2KeyWord<int>(controller, "InnerLoop");
    int OuterLoop = Json2KeyWord<int>(controller, "OuterLoop");
    double LBFGS_eps = Json2KeyWord<double>(controller, "LBFGS_eps");

    std::stringstream ss;

    Geometry geometry = host->getGeometry();
    Molecule tmp(host);
    Molecule h(host);
    Vector parameter(3 * host->AtomCount());

    for (int i = 0; i < host->AtomCount(); ++i) {
        parameter(3 * i) = geometry(i, 0);
        parameter(3 * i + 1) = geometry(i, 1);
        parameter(3 * i + 2) = geometry(i, 2);
    }
    std::chrono::time_point<std::chrono::system_clock> init = std::chrono::system_clock::now(), start = std::chrono::system_clock::now(), end;

    ss << "Step\tCurrent Energy [Eh]\tEnergy Change\tRMSD Change\tt [s]" << std::endl;
    XTBInterface interface;
    interface.InitialiseMolecule(host);

    double final_energy = interface.GFNCalculation(method);

    LBFGSParam<double> param;
    param.epsilon = LBFGS_eps;
    param.max_iterations = InnerLoop;

    LBFGSSolver<double> solver(param);
    LBFGSInterface fun(3 * host->AtomCount());
    fun.setMolecule(host);
    fun.setInterface(&interface);
    fun.setMethod(method);

    double fx;
    json RMSDJsonControl = {
        { "reorder", false },
        { "check", false },
        { "heavy", false },
        { "fragment", -1 },
        { "fragment_reference", -1 },
        { "fragment_target", -1 },
        { "init", -1 },
        { "pt", 0 },
        { "silent", true },
        { "storage", 1.0 },
        { "method", "incr" },
        { "noreorder", true },
        { "threads", 1 }
    };
    RMSDDriver* driver = new RMSDDriver(RMSDJsonControl);

    int atoms_count = host->AtomCount();
    for (int outer = 0; outer < OuterLoop; ++outer) {
        int niter = solver.minimize(fun, parameter, fx);
        parameter = fun.Parameter();

        for (int i = 0; i < atoms_count; ++i) {
            geometry(i, 0) = parameter(3 * i);
            geometry(i, 1) = parameter(3 * i + 1);
            geometry(i, 2) = parameter(3 * i + 2);
        }
        h.setGeometry(geometry);
        //h.setEnergy(final_energy);
        //h.appendXYZFile("curcuma_optim.xyz");
        driver->setReference(tmp);
        driver->setTarget(h);
        driver->start();
        end = std::chrono::system_clock::now();
        ss << outer << "\t" << std::setprecision(9) << fun.m_energy << "\t\t" << std::setprecision(5) << (fun.m_energy - final_energy) * 2625.5 << "\t\t" << std::setprecision(6) << driver->RMSD() << "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0 << std::endl;
        start = std::chrono::system_clock::now();

        final_energy = fun.m_energy;
        tmp = h;
        if (((fun.m_energy - final_energy) * 2625.5 < dE && driver->RMSD() < dRMSD))
            break;
    }

    for (int i = 0; i < host->AtomCount(); ++i) {
        geometry(i, 0) = parameter(3 * i);
        geometry(i, 1) = parameter(3 * i + 1);
        geometry(i, 2) = parameter(3 * i + 2);
    }
    h.setEnergy(final_energy);
    h.setGeometry(geometry);
    ss << "Final"
       << "\t" << std::setprecision(9) << final_energy << "\t\t" << std::setprecision(5) << (fun.m_energy - final_energy) * 2625.5 << "\t\t" << std::setprecision(6) << driver->RMSD() << "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(end - init).count() / 1000.0 << std::endl;

    result_string->append(ss.str());
    result_molecule->LoadMolecule(h);
}

class LBFGSBInterface {
public:
    LBFGSBInterface(int n_)
        : n(n_)
    {
    }
    double operator()(const VectorXd& x, VectorXd& grad)
    {
        double fx = 0.0;
        int attyp[m_atoms];
        double charge = 0;
        std::vector<int> atoms = m_molecule->Atoms();
        double coord[3 * m_atoms];
        double gradient[3 * m_atoms];
        double dist_gradient[3 * m_atoms];

        for (int i = 0; i < m_atoms; ++i) {
            attyp[i] = m_molecule->Atoms()[i];
            coord[3 * i + 0] = x(3 * i + 0) / au;
            coord[3 * i + 1] = x(3 * i + 1) / au;
            coord[3 * i + 2] = x(3 * i + 2) / au;
            std::cout << x(3 * i + 0) << "\t" << x(3 * i + 1) << "\t" << x(3 * i + 2) << std::endl;
        }
        m_interface->UpdateMolecule(coord);
        fx = m_interface->GFNCalculation(m_method, gradient);

        for (int i = 0; i < m_atoms; ++i) {
            grad[3 * i + 0] = gradient[3 * i + 0];
            grad[3 * i + 1] = gradient[3 * i + 1];
            grad[3 * i + 2] = gradient[3 * i + 2];
        }
        m_energy = fx;
        m_parameter = x;
        return fx;
    }

    double LastEnergy() const { return m_energy; }

    double m_energy = 0, m_last_change = 0, m_last_rmsd = 0;
    Vector Parameter() const { return m_parameter; }
    void setMolecule(const Molecule* molecule)
    {
        m_molecule = molecule;
        m_atoms = m_molecule->AtomCount();
    }
    void setInterface(XTBInterface* interface) { m_interface = interface; }
    void setMethod(int method) { m_method = method; }

private:
    int m_iter = 0;
    int m_atoms = 0;
    int n;
    int m_method = 2;
    XTBInterface* m_interface;
    Vector m_parameter;
    const Molecule* m_molecule;
};

inline Molecule OptimiseProtons(const Molecule* host, const json& controller)
{
    //PrintController(controller);
    bool writeXYZ = Json2KeyWord<bool>(controller, "writeXYZ");
    bool printOutput = Json2KeyWord<bool>(controller, "printOutput");
    double dE = Json2KeyWord<double>(controller, "dE");
    double dRMSD = Json2KeyWord<double>(controller, "dRMSD");
    double LBFGS_eps = Json2KeyWord<double>(controller, "LBFGS_eps");
    int method = Json2KeyWord<int>(controller, "GFN");
    int InnerLoop = Json2KeyWord<int>(controller, "InnerLoop");
    int OuterLoop = Json2KeyWord<int>(controller, "OuterLoop");

    Geometry geometry = host->getGeometry();
    Molecule tmp(host);
    Molecule h(host);
    Vector parameter(3 * host->AtomCount()), lb(3 * host->AtomCount()), ub(3 * host->AtomCount());
    double delta = 100;
    for (int i = 0; i < host->AtomCount(); ++i) {
        parameter(3 * i) = geometry(i, 0);
        parameter(3 * i + 1) = geometry(i, 1);
        parameter(3 * i + 2) = geometry(i, 2);

        if (host->Atom(i).first != 1) {
            lb(3 * i) = geometry(i, 0) - delta;
            lb(3 * i + 1) = geometry(i, 1) - delta;
            lb(3 * i + 2) = geometry(i, 2) - delta;

            ub(3 * i) = geometry(i, 0) + delta;
            ub(3 * i + 1) = geometry(i, 1) + delta;
            ub(3 * i + 2) = geometry(i, 2) + delta;
        } else {
            lb(3 * i) = std::numeric_limits<double>::infinity();
            lb(3 * i + 1) = std::numeric_limits<double>::infinity();
            lb(3 * i + 2) = std::numeric_limits<double>::infinity();

            ub(3 * i) = std::numeric_limits<double>::infinity();
            ub(3 * i + 1) = std::numeric_limits<double>::infinity();
            ub(3 * i + 2) = std::numeric_limits<double>::infinity();
        }
    }

    XTBInterface interface;
    interface.InitialiseMolecule(host);

    double final_energy = interface.GFNCalculation(method);

    LBFGSBParam<double> param;
    param.epsilon = LBFGS_eps;
    param.max_iterations = InnerLoop;

    LBFGSBSolver<double> solver(param);
    LBFGSBInterface fun(3 * host->AtomCount());
    fun.setMolecule(host);
    fun.setInterface(&interface);
    fun.setMethod(method);
    double fx;

    json RMSDJsonControl = {
        { "reorder", false },
        { "check", false },
        { "heavy", false },
        { "fragment", -1 },
        { "fragment_reference", -1 },
        { "fragment_target", -1 },
        { "init", -1 },
        { "pt", 0 },
        { "silent", true },
        { "storage", 1.0 },
        { "method", "incr" },
        { "noreorder", true },
        { "threads", 1 }
    };

    RMSDDriver* driver = new RMSDDriver(RMSDJsonControl);

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now(), end;
    if (printOutput) {
        std::cout << "Step\tCurrent Energy [Eh]\tEnergy Change\tRMSD Change\tt [s]" << std::endl;
    }
    int atoms_count = host->AtomCount();
    for (int outer = 0; outer < OuterLoop; ++outer) {
        int niter = solver.minimize(fun, parameter, fx, lb, ub);
        parameter = fun.Parameter();

        for (int i = 0; i < atoms_count; ++i) {
            geometry(i, 0) = parameter(3 * i);
            geometry(i, 1) = parameter(3 * i + 1);
            geometry(i, 2) = parameter(3 * i + 2);
        }
        h.setGeometry(geometry);

        driver->setReference(tmp);
        driver->setTarget(h);
        driver->start();
        if (printOutput) {
            end = std::chrono::system_clock::now();
            std::cout << outer << "\t" << std::setprecision(9) << fun.m_energy << "\t\t" << std::setprecision(5) << (fun.m_energy - final_energy) * 2625.5 << "\t\t" << std::setprecision(6) << driver->RMSD() << "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0 << std::endl;
            start = std::chrono::system_clock::now();
        }
        final_energy = fun.m_energy;
        tmp = h;
        if (writeXYZ) {
            h.setEnergy(final_energy);
            h.appendXYZFile("curcuma_optim.xyz");
        }
        if (((fun.m_energy - final_energy) * 2625.5 < dE && driver->RMSD() < dRMSD))
            break;
    }

    for (int i = 0; i < host->AtomCount(); ++i) {
        geometry(i, 0) = parameter(3 * i);
        geometry(i, 1) = parameter(3 * i + 1);
        geometry(i, 2) = parameter(3 * i + 2);
    }
    h.setEnergy(final_energy);
    h.setGeometry(geometry);
    return h;
}
