/*
 * <Handling optimisation of structures. >
 * Copyright (C) 2020 - 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/capabilities/hessian.h"
#include "src/capabilities/rmsd.h"

#include "src/core/elements.h"
#include "src/core/energycalculator.h"
#include "src/core/fileiterator.h"
#include "src/core/global.h"
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

int SPThread::execute()
{
    auto start = std::chrono::system_clock::now();
    std::vector<double> charges;

    double energy = CurcumaOpt::SinglePoint(&m_molecule, m_controller, m_result, m_param, charges);
    m_final = m_molecule;
    m_final.setEnergy(energy);
    m_scf["e0"] = energy;
    if (charges.size())
        m_scf["charges"] = Tools::DoubleVector2String(charges);
    auto end = std::chrono::system_clock::now();
    m_result = fmt::format("Single Point Energy = {0} Eh ({1} secs)\n", energy, std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
    return 0;
}

int OptThread::execute()
{
    std::vector<double> charges;
    m_final = CurcumaOpt::LBFGSOptimise(&m_molecule, m_controller, m_result, &m_intermediate, m_param, charges, ThreadId(), Basename() + ".opt.trj");
    m_scf["e0"] = m_final.Energy();
    if (charges.size())
        m_scf["charges"] = Tools::DoubleVector2String(charges);
    return 0;
}

int OptMThread::execute()
{
    std::vector<double> charges;

    for (const Molecule& molecule : m_molecules) {
        // auto m = CurcumaOpt::LBFGSOptimise(&molecule, m_controller, m_result, &m_intermediate, m_param, charges, ThreadId(), Basename() + ".opt.trj");
        // m.appendXYZFile(Basename() + ".t" + std::to_string(ThreadId() + ".xyz");
        // m_finals.push_back(m);
    }

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
    m_hessian = Json2KeyWord<int>(m_defaults, "hessian");
    if (m_method.compare("GFNFF") == 0)
        m_threads = 1;
}

void CurcumaOpt::start()
{
    if (m_file_set) {
        getBasename(m_filename);
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
    std::string method = Json2KeyWord<std::string>(m_defaults, "method");

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
#ifdef USE_TBLITE
        if (method.compare("gfn2") == 0) {
            auto dipole = interface.Dipole();
            std::cout << std::endl
                      << std::endl
                      << "Dipole momement (GFN2)" << dipole[0] << " " << dipole[1] << " " << dipole[2] << " : " << sqrt(dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2]) << std::endl;
        }
#endif

        if (m_hessian) {
            Hessian hess(m_method, m_defaults, false);
            hess.setMolecule(*iter);
            hess.setParameter(interface.Parameter());
            hess.CalculateHessian(m_hessian);

            auto hessian = hess.getHessian();
            std::string hessian_string = Tools::Matrix2String(hessian);

            json hjson;
            hjson["atoms"] = hessian.cols() / 3;
            hjson["hessian"] = hessian_string;
            std::ofstream hess_file("hessian.json");
            hess_file << hjson;

            json scfjson;
            scfjson["e0"] = energy;
            if (interface.Charges().size()) {
                std::string charges = Tools::DoubleVector2String(interface.Charges());
                scfjson["charges"] = charges;
            }
            std::ofstream scffile("scf.json");
            scffile << scfjson;
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

            th->setBaseName(Basename());

            th->setMolecule(*iter);
            th->setController(m_defaults);
            th->setThreadId(i);
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
            std::cout << m_defaults << std::endl;
            Hessian hess(m_method, m_defaults, false);
            hess.setParameter(thread->Parameter());
            hess.setMolecule(mol2);
            hess.CalculateHessian(m_hessian);
            auto hessian = hess.getHessian();
            std::string hessian_string = Tools::Matrix2String(hessian);

            json hjson;
            hjson["atoms"] = hessian.cols() / 3;
            hjson["hessian"] = hessian_string;
            std::ofstream hess_file("hessian.json");
            hess_file << hjson;

            std::ofstream scffile("scf.json");
            scffile << thread->SCF();
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

double CurcumaOpt::SinglePoint(const Molecule* initial, const json& controller, std::string& output, json& param, std::vector<double>& charges)
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
    param = interface.Parameter();

    double energy = interface.CalculateEnergy(true, true);
    double store = 0;
#ifdef USE_TBLITE
    if (method.compare("gfn2") == 0) {
        auto dipole = interface.Dipole();
        std::cout << std::endl
                  << std::endl
                  << "Dipole momement (GNF2)" << dipole[0] << " " << dipole[1] << " " << dipole[2] << " : " << sqrt(dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2]) * 2.5418 << std::endl;
        store = sqrt(dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2]);
    }
#endif
    charges = interface.Charges();
    return energy;
}

Molecule CurcumaOpt::LBFGSOptimise(Molecule* initial, const json& controller, std::string& output, std::vector<Molecule>* intermediate, json& parameters, std::vector<double>& charges, int thread, const std::string& basename)
{
    bool printOutput = Json2KeyWord<bool>(controller, "printOutput");
    bool fusion = Json2KeyWord<bool>(controller, "fusion");
    double dE = Json2KeyWord<double>(controller, "dE");
    double dRMSD = Json2KeyWord<double>(controller, "dRMSD");
    double GradNorm = Json2KeyWord<double>(controller, "GradNorm");
    double maxrise = Json2KeyWord<double>(controller, "maxrise");
    std::string method = Json2KeyWord<std::string>(controller, "method");
    int MaxIter = Json2KeyWord<int>(controller, "MaxIter");
    int ConvCount = Json2KeyWord<int>(controller, "ConvCount");
    int SingleStep = Json2KeyWord<int>(controller, "SingleStep");
    /*
        if(Json2KeyWord<int>(controller, "threads") > 1 )
        {
            printOutput = false;
        }
        */
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
    parameters = interface.Parameter();
    double final_energy = interface.CalculateEnergy(true);
    initial->setEnergy(final_energy);
    initial->writeXYZFile(basename + ".t" + std::to_string(thread) + ".xyz");
    std::cout << "Initial energy " << final_energy << "Eh" << std::endl;
    LBFGSParam<double> param;
    param.m = Json2KeyWord<int>(controller, "LBFGS_m");

    param.epsilon = Json2KeyWord<double>(controller, "LBFGS_eps_abs");
    param.epsilon_rel = Json2KeyWord<double>(controller, "LBFGS_eps_rel");
    param.past = Json2KeyWord<int>(controller, "LBFGS_past");
    param.delta = Json2KeyWord<double>(controller, "LBFGS_delta");
    param.linesearch = Json2KeyWord<int>(controller, "LBFGS_LST");
    param.max_linesearch = Json2KeyWord<double>(controller, "LBFGS_ls_iter");
    param.min_step = Json2KeyWord<double>(controller, "LBFGS_min_step");
    param.ftol = Json2KeyWord<double>(controller, "LBFGS_ftol");
    param.wolfe = Json2KeyWord<double>(controller, "LBFGS_wolfe");

    // param.linesearch = Json2KeyWord<int>(controller, "LBFGS_LS");

    LBFGSSolver<double, LineSearchBacktracking> solver(param);
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

    for (iteration = 1; iteration <= MaxIter && perform_optimisation; ++iteration) {
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
                output += fmt::format("{1: ^25} {2: ^{0}f}\n", 2, "FINAL SINGLE POINT ENERGY", final_energy);

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
            output += fmt::format("{1: ^25} {2: ^{0}f}\n", 2, "FINAL SINGLE POINT ENERGY", final_energy);

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
        if ((fun.m_energy - final_energy) * 2625.5 > maxrise && iteration > 10) {
            if (printOutput) {
                output += fmt::format("Energy rises too much!\n");
                output += fmt::format("{0: ^75}\n\n", "*** Geometry Optimisation sufficiantly converged ***");
                output += fmt::format("{1: ^25} {2: ^{0}f}\n", 2, "FINAL SINGLE POINT ENERGY", final_energy);

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
        std::ifstream test_file("stop");
        bool result = test_file.is_open();
        test_file.close();
        if (result) {
            perform_optimisation = false;
            error = true;
        }
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
        if (next.Check() == 0 || (next.Check() == 1 && fusion)) {
            previous = next;
            next.setEnergy(final_energy);
            intermediate->push_back(next);
            next.appendXYZFile(basename + ".t" + std::to_string(thread) + ".xyz");

        } else {
            output += fmt::format("{0: ^75}\n\n", "*** Check next failed! ***");

            perform_optimisation = false;
            error = true;
        }
        }
    }
    end = std::chrono::system_clock::now();
    charges = interface.Charges();

    if (iteration >= MaxIter) {
        output += fmt::format("{0: ^75}\n\n", "*** Maximum number of iterations reached! ***");
        output += fmt::format("{1: ^25} {2: ^{0}f}\n", 2, "FINAL SINGLE POINT ENERGY", final_energy);

        error = true;
    }
    if (error == false) {
        output += fmt::format("{1: ^{0}} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f} {6: ^{0}f}\n", 15, iteration, fun.m_energy, (fun.m_energy - final_energy) * 2625.5, driver->RMSD(), solver.final_grad_norm(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
        output += fmt::format("{0: ^75}\n", "*** Geometry Optimisation converged ***");
        output += fmt::format("{1: ^25} {2: ^{0}f}\n", 2, "FINAL SINGLE POINT ENERGY", final_energy);

        if (printOutput) {
            std::cout << output;
            output.clear();
        }
    } else {
        output += fmt::format("{1: ^{0}} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f} {6: ^{0}f}\n", 15, iteration, fun.m_energy, (fun.m_energy - final_energy) * 2625.5, driver->RMSD(), solver.final_grad_norm(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
        output += fmt::format("{0: ^75}\n\n", "*** Geometry Optimisation Not Really converged ***");
        output += fmt::format("{1: ^25} {2: ^{0}f}\n", 2, "FINAL SINGLE POINT ENERGY", final_energy);

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
