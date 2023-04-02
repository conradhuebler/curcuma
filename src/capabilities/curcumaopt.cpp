/*
 * <Handling optimisation of structures. >
 * Copyright (C) 2020 - 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/capabilities/rmsd.h"

#include "src/core/elements.h"
#include "src/core/energycalculator.h"
#include "src/core/fileiterator.h"
#include "src/core/global.h"
#include "src/core/hessian.h"
#include "src/core/molecule.h"

#include <LBFGS.h>
#include <LBFGSB.h>


#include <iomanip>
#include <iostream>

#include <fmt/color.h>
#include <fmt/core.h>

#ifdef test
#include "optimiser/CppNumericalSolversInterface.h"
#endif

#include "optimiser/LBFGSppInterface.h"

#include "curcumaopt.h"

using Eigen::VectorXd;
using namespace LBFGSpp;

int OptThread::execute()
{
    m_final = CurcumaOpt::LBFGSOptimise(&m_molecule, m_controller, m_result, &m_intermediate);
    return 0;
}

int SPThread::execute()
{
    auto start = std::chrono::system_clock::now();

    double energy = CurcumaOpt::SinglePoint(&m_molecule, m_controller, m_result);
    m_final = m_molecule;
    m_final.setEnergy(energy);
    auto end = std::chrono::system_clock::now();
    m_result = fmt::format("Single Point Energy = {0} Eh ({1} secs)\n", energy, std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
    return 0;
}

CurcumaOpt::CurcumaOpt(const json& controller, bool silent)
    : CurcumaMethod(CurcumaOptJson, controller, silent)
{
    UpdateController(controller);
}

void CurcumaOpt::LoadControlJson()
{
    m_threads = Json2KeyWord<int>(m_defaults, "threads");
    m_writeXYZ = Json2KeyWord<bool>(m_defaults, "writeXYZ");
    m_printoutput = Json2KeyWord<bool>(m_defaults, "printOutput");
    m_dE = Json2KeyWord<double>(m_defaults, "dE");
    m_dRMSD = Json2KeyWord<double>(m_defaults, "dRMSD");
    m_method = Json2KeyWord<std::string>(m_defaults, "method");
    m_charge = Json2KeyWord<double>(m_defaults, "Charge");
    m_spin = Json2KeyWord<double>(m_defaults, "Spin");
    m_singlepoint = Json2KeyWord<bool>(m_defaults, "SinglePoint");
    m_serial = Json2KeyWord<bool>(m_defaults, "serial");
    m_hessian = Json2KeyWord<bool>(m_defaults, "hessian");
    if (m_method == "GFNFF")
        m_threads = 1;
}

void CurcumaOpt::start()
{
    if (m_file_set) {
        FileIterator file(m_filename);
        std::multimap<double, Molecule> results;
        while (!file.AtEnd()) {
            Molecule mol = file.Next();
            mol.setCharge(m_charge);
            mol.setSpin(m_spin);
            m_molecules.push_back(mol);
        }
    }
    if (!m_serial)
        ProcessMolecules(m_molecules);
    else {
        ProcessMoleculesSerial(m_molecules);
    }
}

void CurcumaOpt::ProcessMoleculesSerial(const std::vector<Molecule>& molecules)
{
    EnergyCalculator interface(Json2KeyWord<std::string>(m_defaults, "method"), m_controller["sp"]);

    auto iter = molecules.begin();
    interface.setMolecule(*iter);

    while (iter != molecules.end()) {
        if (iter == molecules.end())
            continue;
        if (iter->AtomCount() == 0)
            continue;
        auto start = std::chrono::system_clock::now();
        interface.updateGeometry(iter->Coords());
        double energy = interface.CalculateEnergy(true, true);
        if (m_hessian) {
            Hessian hess(m_method, m_defaults, m_threads);
            hess.setMolecule(*iter);
            hess.CalculateHessian(true);
        }
        auto end = std::chrono::system_clock::now();
        std::cout << fmt::format("Single Point Energy = {0} Eh ({1} secs)\n", energy, std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
        ++iter;
    }
}

void CurcumaOpt::ProcessMolecules(const std::vector<Molecule>& molecules)
{
    int threads = m_threads;

    CxxThreadPool* pool = new CxxThreadPool;
    pool->setProgressBar(CxxThreadPool::ProgressBarType::Continously);
    pool->setActiveThreadCount(threads);
    pool->StaticPool();
    std::vector<SPThread*> thread_block;
    auto iter = molecules.begin();
    while (iter != molecules.end()) {
        for (int i = 0; i < threads; ++i) {
            if (iter == molecules.end())
                continue;
            if (iter->AtomCount() == 0)
                continue;

            SPThread* th;
            if (!m_singlepoint)
                th = new OptThread;
            else
                th = new SPThread;

            th->setMolecule(*iter);
            th->setController(m_defaults);
            thread_block.push_back(th);
            pool->addThread(th);

            ++iter;
        }
    }
    pool->StartAndWait();
    m_molecules.clear();
    for (auto t : pool->OrderedList()) {
        const SPThread* thread = static_cast<const SPThread*>(t.second);
        if (!thread->Finished()) {
            std::cout << " not finished " << thread->getMolecule().Energy() << std::endl;
            continue;
        }
        // if (m_threads > 1)
        std::cout << thread->Output();

        Molecule* mol2 = new Molecule(thread->getMolecule());
        if (m_hessian) {
            Hessian hess(m_method, m_defaults, m_threads);
            hess.setMolecule(mol2);
            hess.CalculateHessian(true);
        }
        if (!m_singlepoint)
            mol2->appendXYZFile(Optfile());
        m_molecules.push_back(Molecule(mol2));
        if (m_writeXYZ) {
            for (const auto& m : *(thread->Intermediates()))
                m.appendXYZFile(Trjfile());
        }
    }
    delete pool;
}

void CurcumaOpt::clear()
{
    m_molecules.clear();
}

double CurcumaOpt::SinglePoint(const Molecule* initial, const json& controller, std::string& output)
{
    std::string method = Json2KeyWord<std::string>(controller, "method");

    Geometry geometry = initial->getGeometry();
    Molecule tmp(initial);
    Molecule h(initial);
    Vector parameter(3 * initial->AtomCount());

    for (int i = 0; i < initial->AtomCount(); ++i) {
        parameter(3 * i) = geometry(i, 0);
        parameter(3 * i + 1) = geometry(i, 1);
        parameter(3 * i + 2) = geometry(i, 2);
    }

    EnergyCalculator interface(method, controller);
    interface.setMolecule(*initial);

    return interface.CalculateEnergy(true, true);
}

Molecule CurcumaOpt::LBFGSOptimise(const Molecule* initial, const json& controller, std::string& output, std::vector<Molecule>* intermediate)
{
    bool printOutput = Json2KeyWord<bool>(controller, "printOutput");
    double dE = Json2KeyWord<double>(controller, "dE");
    double dRMSD = Json2KeyWord<double>(controller, "dRMSD");
    double LBFGS_eps = Json2KeyWord<double>(controller, "LBFGS_eps");
    double GradNorm = Json2KeyWord<double>(controller, "GradNorm");

    std::string method = Json2KeyWord<std::string>(controller, "method");
    int MaxIter = Json2KeyWord<int>(controller, "MaxIter");
    int StoreIntermediate = Json2KeyWord<int>(controller, "StoreIntermediate");
    int ConvCount = Json2KeyWord<int>(controller, "ConvCount");
    int SingleStep = Json2KeyWord<int>(controller, "SingleStep");

    if(Json2KeyWord<int>(controller, "threads") > 1 )
    {
        printOutput = false;
    }
    bool optH = Json2KeyWord<bool>(controller, "optH");
    std::vector<int> constrain;
    Geometry geometry = initial->getGeometry();
    intermediate->push_back(initial);
    Molecule previous(initial);
    Molecule next(initial);
    Vector parameter(3 * initial->AtomCount()), old_parameter(3 * initial->AtomCount());

    for (int i = 0; i < initial->AtomCount(); ++i) {
        parameter(3 * i) = geometry(i, 0);
        parameter(3 * i + 1) = geometry(i, 1);
        parameter(3 * i + 2) = geometry(i, 2);
        constrain.push_back(initial->Atom(i).first == 1);
    }

    EnergyCalculator interface(method, controller);
    interface.setMolecule(*initial);

    double final_energy = interface.CalculateEnergy(true);

    LBFGSParam<double> param;
    param.epsilon = LBFGS_eps;
    param.m = StoreIntermediate;
    /*
        param.linesearch = 3;
        param.ftol = 1e-6;
        param.wolfe = 0.1;
    */

    LBFGSSolver<double> solver(param);
    LBFGSInterface fun(3 * initial->AtomCount());
    fun.setMolecule(initial);
    fun.setInterface(&interface);
    if (optH)
        fun.setConstrains(constrain);

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
    output += fmt::format("\nCharge {} Spin {}\n\n", initial->Charge(), initial->Spin());
    output += fmt::format("{2: ^{1}} {3: ^{1}} {4: ^{1}} {5: ^{1}} {6: ^{1}} {7: ^{1}}\n", "", 15, "Step", "Current Energy", "Energy Change", "RMSD Change", "Gradient Norm", "time");
    output += fmt::format("{2: ^{1}} {3: ^{1}} {4: ^{1}} {5: ^{1}} {6: ^{1}} {7: ^{1}}\n", "", 15, " ", "[Eh]", "[kJ/mol]", "[A]", "[A]", "[s]");

    if (printOutput) {
        std::cout << output;
        output.clear();
    }
    bool perform_optimisation = true;
    bool error = false;

    int atoms_count = initial->AtomCount();
    int converged = solver.InitializeSingleSteps(fun, parameter, fx);
    int iteration = 0;
    if (converged)
        perform_optimisation = false;

    for (iteration = 0; iteration <= MaxIter && perform_optimisation; ++iteration) {
        old_parameter = parameter;
        try {
            solver.SingleStep(fun, parameter, fx);
            if (fun.isError()) {
                perform_optimisation = false;
            }
        } catch (const std::logic_error& error_result) {
            if (solver.Step() < 1e-8) {
                perform_optimisation = false;
            } else {
                perform_optimisation = false;
                error = true;
                output += fmt::format("LBFGS interface signalled some logic error!\n");
                output += fmt::format(" -- {0: ^75} --\n", error_result.what());
                output += fmt::format("{0: ^75}\n\n", "*** Geometry Optimisation Not Really converged ***");

                if (printOutput) {
                    std::cout << output;
                    output.clear();
                }
                for (int i = 0; i < atoms_count; ++i) {
                    geometry(i, 0) = old_parameter(3 * i);
                    geometry(i, 1) = old_parameter(3 * i + 1);
                    geometry(i, 2) = old_parameter(3 * i + 2);
                }
                next.setGeometry(geometry);
            }
        } catch (const std::runtime_error& error_result) {
            output += fmt::format("LBFGS interface signalled some runtime error!\n");
            output += fmt::format(" -- {0: ^75} --\n", error_result.what());
            output += fmt::format("{0: ^75}\n\n", "*** Geometry Optimisation Not Really converged ***");

            if (printOutput) {
                std::cout << output;
                output.clear();
            }
            error = true;
            perform_optimisation = false;
            for (int i = 0; i < atoms_count; ++i) {
                geometry(i, 0) = old_parameter(3 * i);
                geometry(i, 1) = old_parameter(3 * i + 1);
                geometry(i, 2) = old_parameter(3 * i + 2);
            }
            next.setGeometry(geometry);
        }
        if ((iteration % SingleStep == 0 && perform_optimisation) || fun.isError()) {
            parameter = fun.Parameter();

            for (int i = 0; i < atoms_count; ++i) {
                geometry(i, 0) = parameter(3 * i);
                geometry(i, 1) = parameter(3 * i + 1);
                geometry(i, 2) = parameter(3 * i + 2);
            }
            next.setGeometry(geometry);

            driver->setReference(previous);
            driver->setTarget(next);
            driver->start();
            end = std::chrono::system_clock::now();

#ifdef GCC
            output += fmt::format("{1: ^{0}} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f} {6: ^{0}f}\n", 15, iteration, fun.m_energy, (fun.m_energy - final_energy) * 2625.5, driver->RMSD(), solver.final_grad_norm(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
#else
            output += fmt::format("{1} {2} {3} {4} {5}\n", 15, iteration, fun.m_energy, (fun.m_energy - final_energy) * 2625.5, driver->RMSD(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);

#endif
        start = std::chrono::system_clock::now();
        if (printOutput) {
            std::cout << output;
            output.clear();
        }
        /*
         * Energy = 1
         * RMSD = 2
         * LBFGS Conv = 4
         * Gradient Norm = 8
         * */
        converged = 1 * (abs(fun.m_energy - final_energy) * 2625.5 < dE)
            + 2 * (driver->RMSD() < dRMSD)
            + 4 * (solver.isConverged())
            + 8 * (solver.final_grad_norm() < GradNorm);
        perform_optimisation = ((converged & ConvCount) != ConvCount) && (fun.isError() == 0);
        /*
        std::cout << (abs(fun.m_energy - final_energy) * 2625.5 < 0.05)
                  << " " << (int(driver->RMSD() < 0.01))
                  << " " << solver.isConverged()
                  << " " << (solver.final_grad_norm()  < 0.0002)
                  << std::endl;

        std::cout << converged << " " << minConverged << " " << perform_optimisation << std::endl;
        */
        /*
        if ((abs(fun.m_energy - final_energy) * 2625.5) < 0.05 && driver->RMSD() < 0.01) {
           perform_optimisation = false;
            break;
        }
        */
        final_energy = fun.m_energy;
        if (next.Check() == 0) {
            previous = next;
            previous.setEnergy(final_energy);
            intermediate->push_back(next);
        } else {
            perform_optimisation = false;
            error = true;
        }
        }
    }
    end = std::chrono::system_clock::now();
    if (iteration >= MaxIter) {
        output += fmt::format("{0: ^75}\n\n", "*** Maximum number of iterations reached! ***");
        error = true;
    }
    if (error == false) {
        output += fmt::format("{1: ^{0}} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f} {6: ^{0}f}\n", 15, iteration, fun.m_energy, (fun.m_energy - final_energy) * 2625.5, driver->RMSD(), solver.final_grad_norm(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
        output += fmt::format("{0: ^75}\n\n", "*** Geometry Optimisation converged ***");
        if (printOutput) {
            std::cout << output;
            output.clear();
        }
    } else {
        output += fmt::format("{1: ^{0}} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f} {6: ^{0}f}\n", 15, iteration, fun.m_energy, (fun.m_energy - final_energy) * 2625.5, driver->RMSD(), solver.final_grad_norm(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
        output += fmt::format("{0: ^75}\n\n", "*** Geometry Optimisation Not Really converged ***");
        if (printOutput) {
            std::cout << output;
            output.clear();
        }
    }

    if (next.Check() == 0) {
        for (int i = 0; i < initial->AtomCount(); ++i) {
            geometry(i, 0) = parameter(3 * i);
            geometry(i, 1) = parameter(3 * i + 1);
            geometry(i, 2) = parameter(3 * i + 2);
        }
        previous.setEnergy(final_energy);
        previous.setGeometry(geometry);
    }
    return previous;
}
