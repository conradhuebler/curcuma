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
    { "writeXYZ", true },
    { "printOutput", true },
    { "dE", 0.1 },
    { "dRMSD", 0.01 },
    { "GFN", 2 },
    { "InnerLoop", 20 },
    { "OuterLoop", 100 },
    { "LBFGS_eps", 1e-5 },
    { "Threads", 1 },
    { "Charge", 0 },
    { "Spin", 0 }
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
    const std::vector<Molecule>* Intermediates() const { return &m_intermediate; }

private:
    std::string m_result;
    Molecule m_molecule, m_final;
    json m_controller = OptJsonPrivate;
    std::vector<Molecule> m_intermediate;
};

class CurcumaOpt : public CurcumaMethod {
public:
    CurcumaOpt(const json& controller, bool silent);

    void setFileName(const std::string& filename)
    {
        m_filename = filename;
        m_basename = std::string(m_filename);
        for (int i = 0; i < 4; ++i)
            m_basename.pop_back();

        m_file_set = true;
        m_mol_set = false;
        m_mols_set = false;
    }

    inline void addMolecule(const Molecule& molecule)
    {
        m_molecules.push_back(molecule);
        m_file_set = false;
        m_mol_set = false;
        m_mols_set = true;
    }

    inline void setBaseName(const std::string& basename)
    {
        m_basename = basename;

        std::ofstream tfile1;
        tfile1.open(Optfile());
        tfile1.close();

        std::ofstream tfile2;
        tfile2.open(Trjfile());
        tfile2.close();
    }

    inline std::string Optfile() const { return std::string(m_basename + ".opt.xyz"); }
    inline std::string Trjfile() const { return std::string(m_basename + ".trj.xyz"); }

    void start() override; // TODO make pure virtual and move all main action here

    inline const std::vector<Molecule>* Molecules() const { return &m_molecules; }

    static Molecule LBFGSOptimise(const Molecule* host, const json& controller, std::string& output, std::vector<Molecule>* intermediate);
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

    void ProcessMolecules(const std::vector<Molecule>& molecule);

    std::string m_filename, m_basename = "curcuma_job";
    Molecule m_molecule;
    std::vector<Molecule> m_molecules;
    bool m_file_set = false, m_mol_set = false, m_mols_set = false, m_writeXYZ = true, m_printoutput = true;
    int m_threads = 1, m_GFNmethod = 2;
    double m_dE = 0.1, m_dRMSD = 0.01;
    int m_charge = 0, m_spin = 0;
};
