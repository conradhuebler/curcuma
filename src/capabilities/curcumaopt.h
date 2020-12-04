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

#pragma once

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include "curcumamethod.h"

static json CurcumaOptJson{
    { "Solver", 1 },
    { "writeXYZ", false },
    { "printOutput", true },
    { "dE", 0.1 },
    { "dRMSD", 0.01 },
    { "GFN", 2 },
    { "InnerLoop", 20 },
    { "OuterLoop", 100 },
    { "LBFGS_eps", 1e-5 },
    { "Threads", 1 }
};

const json OptJsonPrivate{
    { "writeXYZ", false },
    { "printOutput", true },
    { "dE", 0.1 },
    { "dRMSD", 0.01 },
    { "GFN", 2 },
    { "InnerLoop", 20 },
    { "OuterLoop", 100 },
    { "LBFGS_eps", 1e-5 }
};

class OptThread : public CxxThread {
public:
    OptThread() = default;
    ~OptThread() = default;

    inline void setMolecule(const Molecule& molecule) { m_molecule = molecule; }
    inline Molecule getMolecule() const { return m_final; }
    int execute();

    inline void setController(const json& controller) { m_controller = controller; }
    std::string Output() const { return m_result; }

private:
    std::string m_result;
    Molecule m_molecule, m_final;
    json m_controller = OptJsonPrivate;
};

class CurcumaOpt : public CurcumaMethod {
public:
    CurcumaOpt(const json& controller, bool silent);

    void setFileName(const std::string& filename)
    {
        m_filename = filename;
        m_file_set = true;
        m_mol_set = false;
        m_mols_set = false;
    }

    void start() override; // TODO make pure virtual and move all main action here

    static Molecule LBFGSOptimise(const Molecule* host, const json& controller);
#ifdef test
    static Molecule CppNumSolvOptimise(const Molecule* host, const json& controller);
#endif
private:
    /* Lets have this for all modules */
    inline nlohmann::json WriteRestartInformation() override { return json(); }

    /* Lets have this for all modules */
    inline bool LoadRestartInformation() override { return true; }

    inline std::string MethodName() const override { return std::string("opt"); }

    /* Lets have all methods read the input/control file */
    void ReadControlFile() override{};

    /* Read Controller has to be implemented for all */
    void LoadControlJson() override;

    std::string m_filename;
    Molecule m_molecule;
    std::vector<Molecule> m_molecules;
    bool m_file_set = false, m_mol_set = false, m_mols_set = false, m_writeXYZ = true, m_printoutput = true;
    int m_threads = 1, m_GFNmethod = 2;
    double m_dE = 0.1, m_dRMSD = 0.01;
};
