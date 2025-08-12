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

#include "src/capabilities/optimiser/lbfgs.h"

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

#include "curcumaopt.h"

using Eigen::VectorXd;
using namespace LBFGSpp;

double LBFGSInterface::operator()(const VectorXd& x, VectorXd& grad)
{
    double fx = 0.0;
    double charge = 0;
    m_interface->updateGeometry(x);
    if (m_interface->HasNan()) {
        m_error = true;
        return 0;
    }
    fx = m_interface->CalculateEnergy(true);
    Geometry gradient = m_interface->Gradient();
    m_error = std::isnan(fx);

    for (int i = 0; i < m_atoms; ++i) {
        grad[3 * i + 0] = gradient.data()[3 * i + 0] * (m_constrains[i]);
        grad[3 * i + 1] = gradient.data()[3 * i + 1] * (m_constrains[i]);
        grad[3 * i + 2] = gradient.data()[3 * i + 2] * (m_constrains[i]);
    }
    m_energy = fx;
    m_parameter = x;
    return fx;
}

void LBFGSInterface::setMolecule(const Molecule* molecule)
{
    m_molecule = molecule;
    m_atoms = m_molecule->AtomCount();
    for (int i = 0; i < m_atoms; ++i)
        m_constrains.push_back(1);
}

int SPThread::execute()
{
    auto start = std::chrono::system_clock::now();
    Vector charges;
    double energy = m_curcumaOpt->SinglePoint(&m_molecule, m_result, charges);

    m_final = m_molecule;
    m_final.setEnergy(energy);
    m_scf["e0"] = energy;
    if (charges.rows())
        m_scf["charges"] = Tools::DoubleVector2String(charges);
    auto end = std::chrono::system_clock::now();
    m_result = fmt::format("Single Point Energy = {0} Eh ({1} secs)\n", energy, std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
    return 0;
}

int OptThread::execute()
{
    Vector charges;
    if (m_optimethod == 0)
        m_final = m_curcumaOpt->LBFGSOptimise(&m_molecule, m_result, &m_intermediate, charges, getThreadId(), Basename() + ".opt.trj");
    else
        m_final = m_curcumaOpt->GPTLBFGS(&m_molecule, m_result, &m_intermediate, charges, getThreadId(), Basename() + ".opt.trj");

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
    m_GradNorm = Json2KeyWord<double>(m_defaults, "GradNorm");
    m_ConvCount = Json2KeyWord<int>(m_defaults, "ConvCount");

    m_method = Json2KeyWord<std::string>(m_defaults, "method");
    m_charge = Json2KeyWord<double>(m_defaults, "Charge");
    m_spin = Json2KeyWord<double>(m_defaults, "Spin");
    m_singlepoint = Json2KeyWord<bool>(m_defaults, "SinglePoint");
    m_serial = Json2KeyWord<bool>(m_defaults, "serial");
    m_hessian = Json2KeyWord<int>(m_defaults, "hessian");
    m_optH = Json2KeyWord<bool>(m_defaults, "optH");
    m_maxiter = Json2KeyWord<int>(m_defaults, "maxiter");
    m_maxrise = Json2KeyWord<int>(m_defaults, "maxrise");
    m_fusion = Json2KeyWord<bool>(m_defaults, "fusion");
    m_optimethod = Json2KeyWord<int>(m_defaults, "optimethod");
    m_inithess = Json2KeyWord<bool>(m_defaults, "inithess");
    m_lambda = Json2KeyWord<double>(m_defaults, "lambda");
    m_diis_hist = Json2KeyWord<int>(m_defaults, "diis_hist");
    m_diis_start = Json2KeyWord<int>(m_defaults, "diis_start");
    m_mo_scheme = Json2KeyWord<bool>(m_defaults, "mo_scheme");
    m_mo_scale = Json2KeyWord<double>(m_defaults, "mo_scale");
    m_mo_homo = Json2KeyWord<int>(m_defaults, "mo_homo");
    m_mo_lumo = Json2KeyWord<int>(m_defaults, "mo_lumo");

    if (m_optimethod == 0) {
        std::cout << "Using external lBFGS module" << std::endl;
    } else {
        std::cout << "Using gpt coded optimisation module" << std::endl;
    }

    if (m_method.compare("GFNFF") == 0)
        m_threads = 1;
}

void CurcumaOpt::start()
{
    if (m_file_set) {
        FileIterator file(Filename());
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
    // Claude Generated: Use new constructor with basename for parameter caching
    EnergyCalculator interface(Json2KeyWord<std::string>(m_defaults, "method"), m_controller["sp"], Basename());
    std::string method = Json2KeyWord<std::string>(m_defaults, "method");

    auto iter = molecules.begin();
    interface.setMolecule(iter->getMolInfo());

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
            if (!m_singlepoint) {
                th = new OptThread(this);
                th->setOptiMethod(m_optimethod);
            } else
                th = new SPThread(this);

            th->setBaseName(Basename());

            th->setMolecule(*iter);
            th->setThreadId(i);
            thread_block.push_back(th);
            pool->addThread(th);

            ++iter;
        }
    }
    pool->StartAndWait();
    m_molecules.clear();
#pragma message "TODO: Check the following code block"
    for (auto t : pool->getFinishedThreads()) {
        const SPThread* thread = static_cast<const SPThread*>(t);
        if (!thread->isFinished()) {
            std::cout << " not finished " << thread->getMolecule().Energy() << std::endl;
            continue;
        }
        // if (m_threads > 1)
        std::cout << thread->Output();

        Molecule* mol2 = new Molecule(thread->getMolecule());
        if (m_hessian) {
            std::cout << m_defaults << std::endl;
            Hessian hess(m_method, m_defaults, false);
            hess.setParameter(m_parameters);
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

double CurcumaOpt::SinglePoint(const Molecule* initial, std::string& output, Vector& charges)
{
    std::string method = m_method; // Json2KeyWord<std::string>(m_controller, "method");

    Geometry geometry = initial->getGeometry();
    Molecule tmp(initial);
    Molecule h(initial);
    Vector parameter(3 * initial->AtomCount());

    for (int i = 0; i < initial->AtomCount(); ++i) {
        parameter(3 * i) = geometry(i, 0);
        parameter(3 * i + 1) = geometry(i, 1);
        parameter(3 * i + 2) = geometry(i, 2);
    }

    // Claude Generated: Use new constructor with basename for parameter caching
    EnergyCalculator interface(method, m_controller["sp"], Basename());
    interface.setMolecule(initial->getMolInfo());
    json param = interface.Parameter();
    double energy = interface.CalculateEnergy(false, true);
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
    if (m_mo_scheme) {
        m_orbital_energies = interface.Energies();
        m_num_electrons = interface.NumElectrons();
        WriteMO(m_mo_homo, m_mo_lumo);
        WriteMOAscii();
    }
    return energy;
}

void CurcumaOpt::WriteMO(int n, int m)
{
    std::ofstream file(Basename() + ".inc");
    std::ostream& out = std::cout;
    double spacing = 0.1;
    auto write_tikz = [&](std::ostream& os) {
        os << "\\begin{tikzpicture}\n";
        int highest_occupied_orbital = m_num_electrons / 2 - 1;
        int lowest_unoccupied_orbital = highest_occupied_orbital + 1;
        int start_orbital = std::max(0, highest_occupied_orbital - m);
        int end_orbital = std::min(static_cast<int>(m_orbital_energies.size()), highest_occupied_orbital + n + 1);

        for (int i = start_orbital; i < end_orbital; ++i) {
            double energy = m_orbital_energies[i] * m_mo_scale;
            os << std::fixed << std::setprecision(4);

            // Check if the orbital is HOMO or LUMO
            if (i == highest_occupied_orbital) {
                os << "\t\\draw[thick, color=red] (1," << energy << ") -- (0.0," << energy << ");\n";
                os << "\t\\node[label=below:{\\tiny " << m_orbital_energies[i] << "}] at (0.5," << energy << ") {};\n";
            } else if (i == lowest_unoccupied_orbital) {
                os << "\t\\draw[thick, color=blue] (1," << energy << ") -- (0.0," << energy << ");\n";
                os << "\t\\node[label=below:{\\tiny " << m_orbital_energies[i] << "}] at (0.5," << energy << ") {};\n";
            } else {
                os << "\t\\draw[thick] (1," << energy << ") -- (0.0," << energy << ");\n";
                os << "\t\\node[label=below:{\\tiny " << m_orbital_energies[i] << "}] at (0.5," << energy << ") {};\n";
            }

            if (i <= highest_occupied_orbital) {
                os << "\t\\draw[thick,<-, color=orange] (0.25," << energy + spacing << ") -- (0.25," << energy - spacing << ");\n";
                os << "\t\\draw[thick,->, color=orange] (0.75," << energy + spacing << ") -- (0.75," << energy - spacing << ");\n";
            }
        }
        os << "\t\\node at (0.5,-1) {{" + Basename() + "}};\n";
        os << "\t\\node at (0.5,-1.25) {\\hphantom{\\tiny (" + Basename() + ")}};\n";
        os << "\\end{tikzpicture}\n";
    };

    write_tikz(out);
    write_tikz(file);

    file.close();
}

void CurcumaOpt::WriteMOAscii()
{
    std::vector<std::string> levels;
    for (size_t i = 0; i < m_orbital_energies.size(); ++i) {
        double energy = m_orbital_energies[i] * m_mo_scale;
        std::string level = std::to_string(m_orbital_energies[i]);
        level.resize(10, ' '); // Adjust the width for alignment
        if (i < m_num_electrons / 2) {
            levels.push_back("|" + std::string(8, '-') + "|" + level + "|<--->|");
        } else {
            levels.push_back("|" + std::string(8, '-') + "|" + level + "|     |");
        }
    }

    std::reverse(levels.begin(), levels.end()); // Reverse to print from top to bottom
    for (const auto& level : levels) {
        std::cout << level << std::endl;
    }
}

Molecule CurcumaOpt::LBFGSOptimise(Molecule* initial, std::string& output, std::vector<Molecule>* intermediate, Vector& charges, int thread, const std::string& basename)
{
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

    // Claude Generated: Use new constructor with basename for parameter caching
    EnergyCalculator interface(m_method, m_controller["opt"], Basename());

    interface.setMolecule(initial->getMolInfo());
    m_parameters = interface.Parameter();
    double final_energy = interface.CalculateEnergy(true);
    initial->setEnergy(final_energy);
    initial->writeXYZFile(basename + ".t" + std::to_string(thread) + ".xyz");
    std::cout << "Initial energy " << final_energy << "Eh" << std::endl;
    LBFGSParam<double> param;
    param.m = Json2KeyWord<int>(m_defaults, "LBFGS_m");

    param.epsilon = Json2KeyWord<double>(m_defaults, "LBFGS_eps_abs");
    param.epsilon_rel = Json2KeyWord<double>(m_defaults, "LBFGS_eps_rel");
    param.past = Json2KeyWord<int>(m_defaults, "LBFGS_past");
    param.delta = Json2KeyWord<double>(m_defaults, "LBFGS_delta");
    param.linesearch = Json2KeyWord<int>(m_defaults, "LBFGS_LST");
    param.max_linesearch = Json2KeyWord<double>(m_defaults, "LBFGS_ls_iter");
    param.min_step = Json2KeyWord<double>(m_defaults, "LBFGS_min_step");
    param.ftol = Json2KeyWord<double>(m_defaults, "LBFGS_ftol");
    param.wolfe = Json2KeyWord<double>(m_defaults, "LBFGS_wolfe");

    LBFGSSolver<double, LineSearchBacktracking> solver(param);
    LBFGSInterface fun(3 * initial->AtomCount());
    fun.setMolecule(initial);
    fun.setInterface(&interface);
    if (m_optH)
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

    if (m_printoutput) {
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

    for (iteration = 1; iteration <= m_maxiter && perform_optimisation; ++iteration) {
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

                if (m_printoutput) {
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

            if (m_printoutput) {
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
        if ((fun.m_energy - final_energy) * 2625.5 > m_maxrise && iteration > 10) {
            if (m_printoutput) {
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

        if ((perform_optimisation) || fun.isError()) {
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
        if (m_printoutput) {
            std::cout << output;
            output.clear();
        }
        /*
         * Energy = 1
         * RMSD = 2
         * LBFGS Conv = 4
         * Gradient Norm = 8
         * */
        converged = 1 * (abs(fun.m_energy - final_energy) * 2625.5 < m_dE)
            + 2 * (driver->RMSD() < m_dRMSD)
            + 4 * (solver.isConverged())
            + 8 * (solver.final_grad_norm() < m_GradNorm);
        perform_optimisation = ((converged & m_ConvCount) != m_ConvCount) && (fun.isError() == 0);
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
        if (next.Check() == 0 || (next.Check() == 1 && m_fusion)) {
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

    if (iteration >= m_maxiter) {
        output += fmt::format("{0: ^75}\n\n", "*** Maximum number of iterations reached! ***");
        output += fmt::format("{1: ^25} {2: ^{0}f}\n", 2, "FINAL SINGLE POINT ENERGY", final_energy);

        error = true;
    }
    if (error == false) {
        output += fmt::format("{1: ^{0}} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f} {6: ^{0}f}\n", 15, iteration, fun.m_energy, (fun.m_energy - final_energy) * 2625.5, driver->RMSD(), solver.final_grad_norm(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
        output += fmt::format("{0: ^75}\n", "*** Geometry Optimisation converged ***");
        output += fmt::format("{1: ^25} {2: ^{0}f}\n", 2, "FINAL SINGLE POINT ENERGY", final_energy);

        if (m_printoutput) {
            std::cout << output;
            output.clear();
        }
    } else {
        output += fmt::format("{1: ^{0}} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f} {6: ^{0}f}\n", 15, iteration, fun.m_energy, (fun.m_energy - final_energy) * 2625.5, driver->RMSD(), solver.final_grad_norm(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
        output += fmt::format("{0: ^75}\n\n", "*** Geometry Optimisation Not Really converged ***");
        output += fmt::format("{1: ^25} {2: ^{0}f}\n", 2, "FINAL SINGLE POINT ENERGY", final_energy);

        if (m_printoutput) {
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

Molecule CurcumaOpt::GPTLBFGS(Molecule* initial, std::string& output, std::vector<Molecule>* intermediate, Vector& charges, int thread, const std::string& basename)
{
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

    // Claude Generated: Use new constructor with basename for parameter caching
    EnergyCalculator interface(m_method, m_controller["opt"], Basename());

    interface.setMolecule(initial->getMolInfo());
    m_parameters = interface.Parameter();
    double final_energy = interface.CalculateEnergy(true);
    initial->setEnergy(final_energy);
    initial->writeXYZFile(basename + ".t" + std::to_string(thread) + ".xyz");
    std::cout << "Initial energy " << final_energy << "Eh" << std::endl;
    LBFGSParam<double> param;
    param.m = Json2KeyWord<int>(m_defaults, "LBFGS_m");

    param.epsilon = Json2KeyWord<double>(m_defaults, "LBFGS_eps_abs");
    param.epsilon_rel = Json2KeyWord<double>(m_defaults, "LBFGS_eps_rel");
    param.past = Json2KeyWord<int>(m_defaults, "LBFGS_past");
    param.delta = Json2KeyWord<double>(m_defaults, "LBFGS_delta");
    param.linesearch = Json2KeyWord<int>(m_defaults, "LBFGS_LST");
    param.max_linesearch = Json2KeyWord<double>(m_defaults, "LBFGS_ls_iter");
    param.min_step = Json2KeyWord<double>(m_defaults, "LBFGS_min_step");
    param.ftol = Json2KeyWord<double>(m_defaults, "LBFGS_ftol");
    param.wolfe = Json2KeyWord<double>(m_defaults, "LBFGS_wolfe");

    // LBFGSSolver<double, LineSearchBacktracking> solver(param);
    // LBFGSInterface fun(3 * initial->AtomCount());
    // fun.setMolecule(initial);
    // fun.setInterface(&interface);
    // if (m_optH)
    // gptfgs.setConstrains(constrain);

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

    if (m_printoutput) {
        std::cout << output;
        output.clear();
    }
    bool perform_optimisation = true;
    bool error = false;

    int atoms_count = initial->AtomCount();
    int converged = 0;
    int iteration = 0;
    if (converged)
        perform_optimisation = false;
    std::vector<double> mass = std::vector<double>(3 * initial->AtomCount(), 0);
    for (int i = 0; i < initial->AtomCount(); ++i) {
        mass[3 * i + 0] = Elements::AtomicMass[initial->Atom(i).first];
        mass[3 * i + 1] = Elements::AtomicMass[initial->Atom(i).first];
        mass[3 * i + 2] = Elements::AtomicMass[initial->Atom(i).first];
    }

    LBFGS gptfgs(Json2KeyWord<int>(m_defaults, "LBFGS_m"));
    gptfgs.initialize(atoms_count, parameter);
    gptfgs.setEnergyCalculator(&interface);
    gptfgs.setOptimMethod(m_optimethod);
    gptfgs.setLambda(m_lambda);
    gptfgs.setMasses(mass);
    gptfgs.setDIIS(m_diis_hist, m_diis_start);
    std::cout << m_lambda << std::endl
              << std::endl;
    if (m_inithess || m_optimethod == 3) {
        Hessian hess(m_method, m_defaults, false);
        hess.setMolecule(*initial);
        hess.setParameter(interface.Parameter());
        hess.CalculateHessian(1);

        auto hessian = hess.getHessian();
        gptfgs.setHessian(hessian);
    }

    for (iteration = 1; iteration <= m_maxiter && perform_optimisation; ++iteration) {
        old_parameter = parameter;

        parameter = gptfgs.step();
        if (gptfgs.isError()) {
            perform_optimisation = false;
        }

        if ((gptfgs.Energy() - final_energy) * 2625.5 > m_maxrise && iteration > 10) {
            if (m_printoutput) {
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

        // if ((iteration % m_singlepoint == 0 && perform_optimisation) || gptfgs.isError()) {
        // parameter = fun.Parameter();

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
        output += fmt::format("{1: ^{0}} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f} {6: ^{0}f}\n", 15, iteration, gptfgs.Energy(), (gptfgs.Energy() - final_energy) * 2625.5, driver->RMSD(), gptfgs.getCurrentGradient().norm(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
#else
        output += fmt::format("{1} {2} {3} {4} {5}\n", 15, iteration, fun.m_energy, (fun.m_energy - final_energy) * 2625.5, driver->RMSD(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);

#endif
        start = std::chrono::system_clock::now();
        if (m_printoutput) {
            std::cout << output;
            output.clear();
        }
        /*
         * Energy = 1
         * RMSD = 2
         * LBFGS Conv = 4
         * Gradient Norm = 8
         * */
        converged = 1 * (abs(gptfgs.Energy() - final_energy) * 2625.5 < m_dE)
            + 2 * (driver->RMSD() < m_dRMSD)
            + 4 * (gptfgs.isConverged())
            + 8 * (gptfgs.getCurrentGradient().norm() < m_GradNorm);
        perform_optimisation = ((converged & m_ConvCount) != m_ConvCount) && (gptfgs.isError() == 0);
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

        // if ((abs(gptfgs.Energy() - final_energy) * 2625.5) < 0.05 && driver->RMSD() < 0.01) {
        //    perform_optimisation = false;
        //     break;
        // }

        final_energy = gptfgs.Energy();
        if (next.Check() == 0 || (next.Check() == 1 && m_fusion)) {
            previous = next;
            next.setEnergy(final_energy);
            intermediate->push_back(next);
            next.appendXYZFile(basename + ".t" + std::to_string(thread) + ".xyz");

        } else {
            output += fmt::format("{0: ^75}\n\n", "*** Check next failed! ***");

            perform_optimisation = false;
            error = true;
        }

        if (m_optimethod == 3 && iteration % 20 == 0) {
            Hessian hess(m_method, m_defaults, false);
            hess.setMolecule(next);
            hess.setParameter(interface.Parameter());
            hess.CalculateHessian(1);

            auto hessian = hess.getHessian();
            gptfgs.setHessian(hessian);
        }
        //}
    }
    end = std::chrono::system_clock::now();
    charges = interface.Charges();

    if (iteration >= m_maxiter) {
        output += fmt::format("{0: ^75}\n\n", "*** Maximum number of iterations reached! ***");
        output += fmt::format("{1: ^25} {2: ^{0}f}\n", 2, "FINAL SINGLE POINT ENERGY", final_energy);

        error = true;
    }
    if (error == false) {
        output += fmt::format("{1: ^{0}} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f} {6: ^{0}f}\n", 15, iteration, gptfgs.Energy(), (gptfgs.Energy() - final_energy) * 2625.5, driver->RMSD(), gptfgs.getCurrentGradient().norm(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
        output += fmt::format("{0: ^75}\n", "*** Geometry Optimisation converged ***");
        output += fmt::format("{1: ^25} {2: ^{0}f}\n", 2, "FINAL SINGLE POINT ENERGY", final_energy);

        if (m_printoutput) {
            std::cout << output;
            output.clear();
        }
    } else {
        output += fmt::format("{1: ^{0}} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f} {6: ^{0}f}\n", 15, iteration, gptfgs.Energy(), (gptfgs.Energy() - final_energy) * 2625.5, driver->RMSD(), gptfgs.getCurrentGradient().norm(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
        output += fmt::format("{0: ^75}\n\n", "*** Geometry Optimisation Not Really converged ***");
        output += fmt::format("{1: ^25} {2: ^{0}f}\n", 2, "FINAL SINGLE POINT ENERGY", final_energy);

        if (m_printoutput) {
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
