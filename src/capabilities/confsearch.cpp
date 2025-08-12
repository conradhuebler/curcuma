/*
 * <Conformational Search based on Molecular Dynamics>
 * Copyright (C) 2022 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/global_config.h"

#include "src/capabilities/confscan.h"
#include "src/capabilities/curcumaopt.h"
#include "src/capabilities/simplemd.h"

#include "src/core/fileiterator.h"
#include "src/core/molecule.h"

#include "src/tools/general.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <fstream>
#include <iostream>
#include <stdio.h>

#include "confsearch.h"
using curcuma::Molecule;

ConfSearch::ConfSearch(const json& controller, bool silent)
    : CurcumaMethod(ConfSearchJson, controller, silent)
{
    UpdateController(controller);
}

ConfSearch::~ConfSearch()
{
}

void ConfSearch::setFile(const std::string& filename)
{
    CurcumaMethod::setFile(filename);

    FileIterator file(Filename());
    while (!file.AtEnd()) {
        Molecule* mol = new Molecule(file.Next());
        m_in_stack.push_back(mol);
        m_topo_matrix = mol->DistanceMatrix().second;
    }
}

bool ConfSearch::Initialise()
{
    return true;
}

void ConfSearch::start()
{
    nlohmann::json md = m_defaults;
    md["unique"] = true;
    for (m_currentT = m_startT; m_currentT >= m_endT; m_currentT -= m_deltaT) {
        std::cout << std::endl
                  << std::endl
                  << " *** Performing MD simulation at " << m_currentT << " K with " << m_repeat << " repetitions and " << m_in_stack.size() << " structures (" << m_repeat * m_in_stack.size() << ")" << std::endl
                  << std::endl;
        md["T"] = m_currentT;
        md["impuls"] = m_currentT;
        std::cout << md << std::endl;
        // std::vector<Molecule*> uniques;
        PerformMolecularDynamics(m_in_stack, md);

        nlohmann::json opt = CurcumaOptJson;
        opt["method"] = m_method;
        opt["threads"] = m_threads;
        PerformOptimisation("ff", opt);

        nlohmann::json scan = ConfSearchJson;
        scan["rmsdmethod"] = "hybrid";
        scan["fewerFile"] = true;
        scan["threads"] = m_threads;
        scan["method"] = m_method;

        PerformFilter("ff", scan);

        for (int i = 0; i < m_in_stack.size(); ++i)
            delete m_in_stack[i];
        m_in_stack.clear();
        double energy = 0;
        FileIterator file("confsearch.unique.opt.accepted.xyz");
        while (!file.AtEnd()) {
            Molecule* mol = new Molecule(file.Next());
            if ((m_topo_matrix - mol->DistanceMatrix().second).cwiseAbs().sum() != 0) {
                delete mol;
                continue;
            }
            if (energy < 0) {
                if ((mol->Energy() < energy) * 2625.5 < m_energy_window)
                    m_in_stack.push_back(mol);
                else
                    delete mol;
            } else {
                m_in_stack.push_back(mol);
                energy = mol->Energy();
            }
        }
    }
}

std::string ConfSearch::PerformMolecularDynamics(const std::vector<Molecule*>& molecules, const nlohmann::json& parameter)
{
    CxxThreadPool* pool = new CxxThreadPool;
    int index = 0;
    std::cout << "Filling pool" << std::endl;
    for (int repeat = 0; repeat < m_repeat; ++repeat) {
        for (int i = 0; i < molecules.size(); ++i) {
            MDThread* thread = new MDThread(parameter);
            thread->setThreadId(index);
            thread->setBasename("confsearch");
            thread->setMolecule(molecules[i * repeat]);
            index++;
            pool->addThread(thread);
        }
    }
    if (m_method.compare("gfnff") == 0)
        pool->setActiveThreadCount(1);
    else
        pool->setActiveThreadCount(m_threads);
    pool->StartAndWait();

    std::string file = "confsearch.unique.xyz";
    std::ofstream result_file;
    result_file.open(file);
    result_file.close();

    for (const auto& thread : pool->getFinishedThreads()) {
        auto structures = static_cast<MDThread*>(thread)->MDDriver()->UniqueMolecules();
        int index = 0;
        for (const auto* molecule : structures) {
            if (index != 0 || molecules.size() == 0)
                molecule->appendXYZFile(file);
            index++;
        }
    }
    delete pool;
    return file;
}

std::string ConfSearch::PerformOptimisation(const std::string& f, const nlohmann::json& parameter)
{
    std::string file = "confsearch.unique";
    std::ofstream result_file;
    result_file.open(file);
    result_file.close();

    CurcumaOpt optimise(parameter, false);
    optimise.setFileName("confsearch.unique.xyz");
    optimise.overrideBasename(file);
    optimise.start();

    return file;
}

std::string ConfSearch::PerformFilter(const std::string& f, const nlohmann::json& parameter)
{
    ConfScan* scan = new ConfScan(parameter, false);
    scan->setFileName("confsearch.unique.opt.xyz");
    scan->start();
    return std::string("fff");
}

nlohmann::json ConfSearch::WriteRestartInformation()
{
    nlohmann::json restart;
    return restart;
}

bool ConfSearch::LoadRestartInformation()
{
    return true;
}

void ConfSearch::ReadControlFile()
{
}

void ConfSearch::LoadControlJson()
{
    m_method = Json2KeyWord<std::string>(m_defaults, "method");
    m_thermostat = Json2KeyWord<std::string>(m_defaults, "thermostat");
    m_rattle = Json2KeyWord<bool>(m_defaults, "rattle");
    m_spin = Json2KeyWord<int>(m_defaults, "spin");
    m_charge = Json2KeyWord<int>(m_defaults, "charge");
    //    m_single_step = Json2KeyWord<double>(m_defaults, "dT"); // * fs2amu;
    m_time = Json2KeyWord<double>(m_defaults, "time");
    m_startT = Json2KeyWord<double>(m_defaults, "startT");
    m_endT = Json2KeyWord<double>(m_defaults, "endT");
    m_deltaT = Json2KeyWord<double>(m_defaults, "deltaT");
    m_repeat = Json2KeyWord<int>(m_defaults, "repeat");
    m_rmsd = Json2KeyWord<double>(m_defaults, "rmsd");
    m_threads = Json2KeyWord<int>(m_defaults, "threads");
    m_energy_window = Json2KeyWord<double>(m_defaults, "energy_window");
    m_dT = Json2KeyWord<double>(m_defaults, "dT");
}
