/*
 * <Handling optimisation of structures. >
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

#include "src/capabilities/rmsd.h"

#include "src/core/elements.h"
#include "src/core/fileiterator.h"
#include "src/core/global.h"
#include "src/core/molecule.h"
#include "src/core/xtbinterface.h"

#include <LBFGS.h>
#include <LBFGSB.h>

#ifdef test
#include "cppoptlib/function.h"
#include "cppoptlib/solver/lbfgs.h"
#endif

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
    // OptimiseGeometryThreaded(&m_molecule, &m_result, &m_final, m_controller);
    m_final = CurcumaOpt::LBFGSOptimise(&m_molecule, m_controller, m_result, &m_intermediate);
    std::cout << m_result << std::endl
              << std::endl;
    return 0;
}

int SPThread::execute()
{
    // OptimiseGeometryThreaded(&m_molecule, &m_result, &m_final, m_controller);
    auto start = std::chrono::system_clock::now();

    double energy = CurcumaOpt::SinglePoint(&m_molecule, m_controller, m_result);
    m_final = m_molecule;
    m_final.setEnergy(energy);
    auto end = std::chrono::system_clock::now();

    m_result = fmt::format("Single Point Energy = {0} Eh ({1} secs)\n", energy, std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
    std::cout << m_result;
    return 0;
}

CurcumaOpt::CurcumaOpt(const json& controller, bool silent)
    : CurcumaMethod(CurcumaOptJson, controller, silent)
{
    std::cout << controller << std::endl;
    UpdateController(controller);
}

void CurcumaOpt::LoadControlJson()
{
    m_threads = Json2KeyWord<int>(m_defaults, "threads");
    m_writeXYZ = Json2KeyWord<bool>(m_defaults, "writeXYZ");
    m_printoutput = Json2KeyWord<bool>(m_defaults, "printOutput");
    m_dE = Json2KeyWord<double>(m_defaults, "dE");
    m_dRMSD = Json2KeyWord<double>(m_defaults, "dRMSD");
    m_GFNmethod = Json2KeyWord<int>(m_defaults, "GFN");
    m_charge = Json2KeyWord<double>(m_defaults, "Charge");
    m_spin = Json2KeyWord<double>(m_defaults, "Spin");
    m_singlepoint = Json2KeyWord<bool>(m_defaults, "SinglePoint");
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

    ProcessMolecules(m_molecules);
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
        if (!thread->Finished())
            continue;
        if (m_threads > 1)
            std::cout << thread->Output();

        Molecule* mol2 = new Molecule(thread->getMolecule());
        mol2->appendXYZFile(Optfile());
        m_molecules.push_back(Molecule(mol2));
        if (m_writeXYZ) {
            for (const auto& m : *(thread->Intermediates()))
                m.appendXYZFile(Trjfile());
        }
    }
    delete pool;
}

#ifdef test
Molecule CurcumaOpt::CppNumSolvOptimise(const Molecule* host, const json& controller)
{
    using Solver = cppoptlib::solver::Lbfgs<CppNumSolvInterface>;

    PrintController(controller);
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

    CppNumSolvInterface xtb_function;
    xtb_function.setMolecule(host);
    xtb_function.setInterface(&interface);
    xtb_function.setMethod(method); //cppoptlib::solver::State
    Solver solver;
    double fx;
    // Evaluate
    auto state = xtb_function.Eval(parameter);
    std::cout << xtb_function(parameter) << " = " << state.value << std::endl;
    std::cout << state.x << std::endl;
    std::cout << state.gradient << std::endl;
    std::cout << state.hessian << std::endl;

    auto [solution, solver_state] = solver.Minimize(xtb_function, parameter);
    /*
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
        printf("Step\tCurrent Energy\t\tEnergy Change\t\tRMSD Change\t   t\n");
        printf("\t[Eh]\t\t\t[kJ/mol]\t\t [A]\t\t   [s]\n");
    }
    int atoms_count = host->AtomCount();
    for (int outer = 0; outer < OuterLoop; ++outer) {
        int niter = 0;
        try {
            niter = solver.Minimize(xtb_function, parameter);
        } catch (const std::logic_error& error) {
            break;
        } catch (const std::runtime_error& error) {
            break;
        }

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
            printf("%i\t%8.5f\t\t%8.5f\t\t%8.5f\t%8.5f\n", outer, fun.m_energy, (fun.m_energy - final_energy) * 2625.5 ,  driver->RMSD(),std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0  );
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
    */
    parameter = solution.x;
    for (int i = 0; i < host->AtomCount(); ++i) {
        geometry(i, 0) = parameter(3 * i);
        geometry(i, 1) = parameter(3 * i + 1);
        geometry(i, 2) = parameter(3 * i + 2);
    }
    h.setEnergy(final_energy);
    h.setGeometry(geometry);
    return h;
}
#endif

double CurcumaOpt::SinglePoint(const Molecule* initial, const json& controller, std::string& output)
{
    bool printOutput = Json2KeyWord<bool>(controller, "printOutput");
    double dE = Json2KeyWord<double>(controller, "dE");
    double dRMSD = Json2KeyWord<double>(controller, "dRMSD");
    double LBFGS_eps = Json2KeyWord<double>(controller, "LBFGS_eps");
    int method = Json2KeyWord<int>(controller, "GFN");
    int InnerLoop = Json2KeyWord<int>(controller, "InnerLoop");
    int OuterLoop = Json2KeyWord<int>(controller, "OuterLoop");
    if (Json2KeyWord<int>(controller, "threads") > 1) {
        printOutput = false;
    }

    Geometry geometry = initial->getGeometry();
    Molecule tmp(initial);
    Molecule h(initial);
    Vector parameter(3 * initial->AtomCount());

    for (int i = 0; i < initial->AtomCount(); ++i) {
        parameter(3 * i) = geometry(i, 0);
        parameter(3 * i + 1) = geometry(i, 1);
        parameter(3 * i + 2) = geometry(i, 2);
    }

    XTBInterface interface;
    interface.InitialiseMolecule(initial);

    return interface.GFNCalculation(method);
}

Molecule CurcumaOpt::LBFGSOptimise(const Molecule* initial, const json& controller, std::string& output, std::vector<Molecule>* intermediate)
{
    bool printOutput = Json2KeyWord<bool>(controller, "printOutput");
    double dE = Json2KeyWord<double>(controller, "dE");
    double dRMSD = Json2KeyWord<double>(controller, "dRMSD");
    double LBFGS_eps = Json2KeyWord<double>(controller, "LBFGS_eps");
    int method = Json2KeyWord<int>(controller, "GFN");
    int InnerLoop = Json2KeyWord<int>(controller, "InnerLoop");
    int OuterLoop = Json2KeyWord<int>(controller, "OuterLoop");
    if(Json2KeyWord<int>(controller, "threads") > 1 )
    {
        printOutput = false;
    }

    Geometry geometry = initial->getGeometry();
    intermediate->push_back(initial);
    Molecule tmp(initial);
    Molecule h(initial);
    Vector parameter(3 * initial->AtomCount());

    for (int i = 0; i < initial->AtomCount(); ++i) {
        parameter(3 * i) = geometry(i, 0);
        parameter(3 * i + 1) = geometry(i, 1);
        parameter(3 * i + 2) = geometry(i, 2);
    }

    XTBInterface interface;
    interface.InitialiseMolecule(initial);

    double final_energy = interface.GFNCalculation(method);

    LBFGSParam<double> param;
    param.epsilon = LBFGS_eps;
    param.max_iterations = InnerLoop;

    LBFGSSolver<double> solver(param);
    LBFGSInterface fun(3 * initial->AtomCount());
    fun.setMolecule(initial);
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
    output += fmt::format("Charge {} Spin {}\n\n", initial->Charge(), initial->Spin());
    output += fmt::format("{2: ^{1}} {3: ^{1}} {4: ^{1}} {5: ^{1}} {6: ^{1}}\n", "", 15, "Step", "Current Energy", "Energy Change", "RMSD Change", "time");
    output += fmt::format("{2: ^{1}} {3: ^{1}} {4: ^{1}} {5: ^{1}} {6: ^{1}}\n", "", 15, " ", "[Eh]", "[kJ/mol]", "[A]", "[s]");

    if (printOutput) {
        std::cout << output;
        output.clear();
    }

    int atoms_count = initial->AtomCount();
    for (int outer = 0; outer < OuterLoop; ++outer) {
        int niter = 0;
        try {
            niter = solver.minimize(fun, parameter, fx);
        } catch (const std::logic_error& error) {
            output += fmt::format("LBFGS interface signalled some logic error!\n");
            output += fmt::format("{0: ^75}\n\n", "*** Geometry Optimisation Not Really converged ***");

            if (printOutput) {
                std::cout << output;
                output.clear();
            }
            break;
        } catch (const std::runtime_error& error) {
            output += fmt::format("LBFGS interface signalled some runtime error!\n");
            output += fmt::format("{0: ^75}\n\n", "*** Geometry Optimisation Not Really converged ***");

            if (printOutput) {
                std::cout << output;
                output.clear();
            }
            break;
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
        end = std::chrono::system_clock::now();
        output += fmt::format("{1: ^{0}} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f}\n", 15, outer, fun.m_energy, (fun.m_energy - final_energy) * 2625.5, driver->RMSD(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
        start = std::chrono::system_clock::now();
        if (printOutput) {
            std::cout << output;
            output.clear();
        }
        final_energy = fun.m_energy;
        if (tmp.Check() == 0) {
            tmp = h;
            h.setEnergy(final_energy);
            intermediate->push_back(h);
        } else {
            output += fmt::format("{0: ^75}\n\n", "*** Geometry Optimisation  Not Really converged ***");
            if (printOutput) {
                std::cout << output;
                output.clear();
            }
            break;
        }
        if (((fun.m_energy - final_energy) * 2625.5 < dE && driver->RMSD() < dRMSD)) {
            output += fmt::format("{0: ^75}\n\n", "*** Geometry Optimisation converged ***");
            if (printOutput) {
                std::cout << output;
                output.clear();
            }
            break;
        }
    }

    if (tmp.Check() == 0) {
        for (int i = 0; i < initial->AtomCount(); ++i) {
            geometry(i, 0) = parameter(3 * i);
            geometry(i, 1) = parameter(3 * i + 1);
            geometry(i, 2) = parameter(3 * i + 2);
        }
        h.setEnergy(final_energy);
        h.setGeometry(geometry);
    }
    return h;
}
