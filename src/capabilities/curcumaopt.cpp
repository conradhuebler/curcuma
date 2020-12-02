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
#include "src/core/global.h"
#include "src/core/molecule.h"
#include "src/core/xtbinterface.h"

#include <LBFGS.h>
#include <LBFGSB.h>

#include "function.h"
#include "solver/lbfgs.h"

#include <iomanip>
#include <iostream>

#include "optimiser/CppNumericalSolversInterface.h"
#include "optimiser/LBFGSppInterface.h"

#include "curcumaopt.h"

using Eigen::VectorXd;
using namespace LBFGSpp;

int OptThread::execute()
{
    // OptimiseGeometryThreaded(&m_molecule, &m_result, &m_final, m_controller);
    m_final = CurcumaOpt::LBFGSOptimise(&m_molecule, m_controller);
    return 0;
}

CurcumaOpt::CurcumaOpt(const json& controller, bool silent)
    : CurcumaMethod(CurcumaOptJson, controller, silent)
{
    UpdateController(controller);
}

void CurcumaOpt::LoadControlJson()
{
    // m_cutoff = Json2KeyWord<double>(m_defaults, "Cutoff");
    // m_temp = Json2KeyWord<double>(m_defaults, "Temp");
}

void CurcumaOpt::start()
{
    int threads = 12; //MaxThreads();
    if (m_file_set) {
        std::string outfile = std::string(m_filename);
        for (int i = 0; i < 4; ++i)
            outfile.pop_back();
        outfile += "_opt.xyz";

        json key = CurcumaOptJson;
        key = MergeJson(key, m_controller);
        std::vector<Molecule> mols;
        FileIterator file(m_filename);
        std::multimap<double, Molecule> results;
        while (!file.AtEnd()) {
            mols.push_back(file.Next());
            //Molecule mol = file.Next();
            //Molecule mol2 = LBFGSOptimise(&mol, key);
            // CppNumSolvOptimise(&mol, key);
            //mol2.writeXYZFile(outfile);
            //results.insert(std::pair<double, Molecule>(mol2.Energy(), mol2));
        }

        CxxThreadPool* pool = new CxxThreadPool;
        //pool->RedirectOutput(&m_curcuma_progress);
        pool->setProgressBar(CxxThreadPool::ProgressBarType::Continously);
        pool->setActiveThreadCount(threads);
        std::vector<OptThread*> thread_block;
        auto iter = mols.begin();
        while (iter != mols.end()) {
            for (int i = 0; i < threads; ++i) {
                if (iter == mols.end())
                    continue;

                auto pair = *iter;

                OptThread* th = new OptThread;
                th->setMolecule(pair);
                th->setController(OptJsonPrivate);
                thread_block.push_back(th);
                pool->addThread(th);

                ++iter;
            }
        }

        std::cout << "Batch calculation started! " << std::endl;
        pool->StartAndWait();
        std::cout << "Batch evaluation ... " << std::endl;
        for (auto* t : pool->Finished()) {
            const OptThread* thread = static_cast<const OptThread*>(t);
            Molecule* mol2 = new Molecule(thread->getMolecule());
            mol2->appendXYZFile(outfile);
            //delete thread;
        }
        delete pool;
    }
}

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

Molecule CurcumaOpt::LBFGSOptimise(const Molecule* host, const json& controller)
{
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
        printf("Step\tCurrent Energy\t\tEnergy Change\t\tRMSD Change\t   t\n");
        printf("\t[Eh]\t\t\t[kJ/mol]\t\t [A]\t\t   [s]\n");
    }
    int atoms_count = host->AtomCount();
    for (int outer = 0; outer < OuterLoop; ++outer) {
        int niter = 0;
        try {
            niter = solver.minimize(fun, parameter, fx);
        } catch (const std::logic_error& error) {
            break;
        } catch (const std::runtime_error& error) {
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
        if (printOutput) {
            end = std::chrono::system_clock::now();
            printf("%i\t%8.5f\t\t%8.5f\t\t%8.5f\t%8.5f\n", outer, fun.m_energy, (fun.m_energy - final_energy) * 2625.5, driver->RMSD(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
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
