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

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include <fstream>
#include <iostream>
#include <stdio.h>

#include "confsearch.h"

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
    m_filename = filename;
    std::string name = std::string(filename);
    for (int i = 0; i < 4; ++i)
        name.pop_back();
    m_basename = name;

    FileIterator file(m_filename);
    while (!file.AtEnd()) {
        Molecule* mol = new Molecule(file.Next());
        m_in_stack.push_back(mol);
    }
}

bool ConfSearch::Initialise()
{
    return true;
}

void ConfSearch::start()
{
    nlohmann::json md = CurcumaMDJson;
    md["gfn"] = m_gfn;
    md["impuls_scaling"] = 0.75;
    md["dt"] = 0.5;
    md["velo"] = 4;
    md["thermostat_steps"] = 400;
    md["hmass"] = 1;
    md["berendson"] = 200;
    md["maxtime"] = m_time;
    md["rmsd"] = m_rmsd;
    md["resuce"] = true;
    md["print"] = m_time / 3.0;
    for (m_currentT = m_startT; m_currentT >= m_endT; m_currentT -= m_deltaT) {
        std::cout << std::endl
                  << std::endl
                  << " *** Performing MD simulation at " << m_currentT << " K with " << m_repeat << " repetitions and " << m_in_stack.size() << " structures (" << m_repeat * m_in_stack.size() << ")" << std::endl
                  << std::endl;
        md["T"] = m_currentT;
        md["impuls"] = m_currentT;
        std::cout << md << std::endl;
        std::vector<Molecule*> uniques;
        PerformMolecularDynamics(m_in_stack, md);
        //  for (Molecule* r : result)
        //      uniques.push_back(r);

        nlohmann::json opt = CurcumaOptJson;
        opt["gfn"] = m_gfn;
        opt["threads"] = m_threads;
        PerformOptimisation(std::string("fff"), opt);

        nlohmann::json scan = ConfSearchJson;
        scan["rmsdmethod"] = "hybrid";
        scan["fewerFile"] = true;
        scan["threads"] = m_threads;
        scan["gfn"] = m_gfn;

        PerformFilter(std::string("fff"), scan);
        /*
        for (int i = 0; i < m_in_stack.size(); ++i)
            delete m_in_stack[i];
        m_in_stack.clear();
        */
        std::string name = m_basename;
        setFile(m_basename + "." + std::to_string(m_currentT) + ".opt.accepted.xyz");
        m_basename = name;

        /*
        for (int i = 0; i < uniques.size(); ++i)
            delete uniques[i];

        for (Molecule* r : optimised)
        {
            std::cout << r->Energy() << std::endl;
            m_in_stack.push_back(new Molecule(r));
            delete r;
        }
        */
    }
}

std::string ConfSearch::PerformMolecularDynamics(const std::vector<Molecule*>& molecules, const nlohmann::json& parameter)
{
    CxxThreadPool* pool = new CxxThreadPool;
    int index = 0;
    for (int repeat = 0; repeat < m_repeat; ++repeat) {
        for (int i = 0; i < molecules.size(); ++i) {
            MDThread* thread = new MDThread(index, parameter);
            thread->setMolecule(molecules[i * repeat]);
            index++;
            pool->addThread(thread);
        }
    }
    if (m_gfn == 66)
        pool->setActiveThreadCount(1);
    else
        pool->setActiveThreadCount(m_threads);
    pool->StartAndWait();
    std::string file = m_basename + "." + std::to_string(m_currentT) + ".unqiues.xyz";

    for (const auto& thread : pool->Finished()) {
        auto structures = static_cast<MDThread*>(thread)->MDDriver()->UniqueMolecules();
        for (const auto* molecule : structures) {
            molecule->appendXYZFile(file);
        }
    }
    delete pool;
    return file;
}

std::string ConfSearch::PerformOptimisation(const std::string& f, const nlohmann::json& parameter)
{
    std::string file = m_basename + "." + std::to_string(m_currentT) + ".unqiues.xyz";

    CurcumaOpt optimise(parameter, false);
    optimise.setFileName(file);
    optimise.setBaseName(m_basename + "." + std::to_string(m_currentT));
    optimise.start();

    return file;
    /*
    CurcumaOpt opt(parameter, false);
    for (Molecule* molecule : molecules) {
        opt.addMolecule(Molecule(molecule));
        delete molecule;
    }
    opt.start();
    */
    /*
    std::vector<Molecule*> optimised;
    auto mols = opt.Molecules();
    for (const Molecule& mol : *mols) {
        optimised.push_back(new Molecule(mol));
        mol.appendXYZFile(m_basename + std::to_string(m_currentT) + ".optimised.xyz");
    }
    return optimised;
    */
}

std::string ConfSearch::PerformFilter(const std::string& f, const nlohmann::json& parameter)
{
    ConfScan* scan = new ConfScan(parameter, false);
    scan->setFileName(m_basename + "." + std::to_string(m_currentT) + ".opt.xyz");
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
    m_gfn = Json2KeyWord<int>(m_defaults, "GFN");
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
}
