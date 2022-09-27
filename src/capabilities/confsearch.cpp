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
    for (double T = m_startT; T >= m_endT; T -= m_deltaT) {
        md["T"] = T;
        md["impuls"] = T;
        std::vector<Molecule*> uniques;
        for (int i = 0; i < m_repeat; ++i) {
            auto result = PerformMolecularDynamics(m_in_stack, md);
            for (Molecule* r : result)
                uniques.push_back(r);
        }

        nlohmann::json opt = CurcumaOptJson;
        opt["gfn"] = m_gfn;
        std::vector<Molecule*> optimised = PerformOptimisation(uniques, opt);
        /*
        nlohmann::json scan;
        PerformFilter(optimised, scan);
        */
        for (int i = 0; i < m_in_stack.size(); ++i)
            delete m_in_stack[i];

        for (int i = 0; i < uniques.size(); ++i)
            delete uniques[i];

        for (Molecule* r : optimised)
            m_in_stack.push_back(r);
    }
}

std::vector<Molecule*> ConfSearch::PerformMolecularDynamics(const std::vector<Molecule*>& molecules, const nlohmann::json& parameter)
{

    CxxThreadPool* pool = new CxxThreadPool;
    for (int i = 0; i < molecules.size(); ++i) {
        MDThread* thread = new MDThread(i, parameter);
        thread->setMolecule(molecules[i]);
        pool->addThread(thread);
    }
    pool->setActiveThreadCount(1);
    pool->StartAndWait();

    std::vector<Molecule*> unique;
    for (const auto& thread : pool->Finished()) {
        auto structures = static_cast<MDThread*>(thread)->MDDriver()->UniqueMolecules();
        for (const auto* molecule : structures) {
            unique.push_back(new Molecule(molecule));
            molecule->appendXYZFile("unqiues.xyz");
        }
    }
    return unique;
}

std::vector<Molecule*> ConfSearch::PerformOptimisation(const std::vector<Molecule*>& molecules, const nlohmann::json& parameter)
{
    CurcumaOpt opt(parameter, false);
    for (Molecule* molecule : molecules) {
        opt.addMolecule(Molecule(molecule));
        delete molecule;
    }
    opt.start();
    std::vector<Molecule*> optimised;
    auto mols = opt.Molecules();
    for (const Molecule& mol : *mols) {
        optimised.push_back(new Molecule(mol));
        mol.appendXYZFile("optimised.xyz");
    }
    return optimised;
}

void ConfSearch::PerformFilter(const std::vector<Molecule*>& molecules, const nlohmann::json& parameter)
{
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
    m_time = Json2KeyWord<double>(m_defaults, "time") * fs2amu;
    m_startT = Json2KeyWord<double>(m_defaults, "startT");
    m_endT = Json2KeyWord<double>(m_defaults, "endT");
    m_deltaT = Json2KeyWord<double>(m_defaults, "deltaT");
    m_repeat = Json2KeyWord<int>(m_defaults, "repeat");
}
