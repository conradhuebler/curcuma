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

#pragma once

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include "optimiser/LBFGSppInterface.h"

#include "curcumamethod.h"

static json CurcumaOptJson{
    { "writeXYZ", true },
    { "printOutput", true },
    { "dE", 0.1 },
    { "dRMSD", 0.01 },
    { "method", "uff" },
    { "MaxIter", 5000 },
 //   { "LBFGS_LS", 3},
    { "LBFGS_m", 2000 },
    { "LBFGS_past", 0},
    { "LBFGS_eps_abs", 1e-5 },
    { "LBFGS_eps_rel", 1e-5 },
    { "LBFGS_delta", 0 },
    { "LBFGS_LST", 3},
    { "LBFGS_ls_iter", 2 },
    { "LBFGS_min_step", 1e-4 },
    { "LBFGS_ftol", 1e-4 },
    { "LBFGS_wolfe", 0.9 },
    { "SingleStep", 1 },
    { "ConvCount", 11 },
    { "GradNorm", 1e-4 },
    { "Threads", 1 },
    { "Charge", 0 },
    { "Spin", 0 },
    { "SinglePoint", false },
    { "optH", false },
    { "serial", false },
    { "hessian", 0 }
};

const json OptJsonPrivate{
    { "writeXYZ", false },
    { "printOutput", true },
    { "dE", 0.1 },
    { "dRMSD", 0.01 },
    { "GFN", 2 },
    { "MaxIter", 2000 },
    { "LBFGS_eps", 1e-5 }
};

class SPThread : public CxxThread {
public:
    SPThread() = default;
    ~SPThread() = default;

    inline void setMolecule(const Molecule& molecule) { m_molecule = molecule; }
    inline Molecule getMolecule() const { return m_final; }
    virtual int execute() override;

    inline void setController(const json& controller) { m_controller = controller; }
    std::string Output() const { return m_result; }
    const std::vector<Molecule>* Intermediates() const { return &m_intermediate; }
    void setBaseName(const std::string& basename) { m_basename = basename; }
    std::string Basename() const { return m_basename; }

protected:
    std::string m_result;
    Molecule m_molecule, m_final;
    json m_controller = OptJsonPrivate;
    std::vector<Molecule> m_intermediate;
    std::string m_basename;
};

class OptThread : public SPThread {
public:
    OptThread() = default;
    ~OptThread() = default;

    int execute() override;
};

class CurcumaOpt : public CurcumaMethod {
public:
    CurcumaOpt(const json& controller, bool silent);

    void setFileName(const std::string& filename)
    {
        m_filename = filename;

        getBasename(filename);
        std::ofstream tfile1;
        tfile1.open(Optfile());
        tfile1.close();

        std::ofstream tfile2;
        tfile2.open(Trjfile());
        tfile2.close();

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
    /*
    inline void setBaseName()
    {
        //m_basename = basename;

        std::ofstream tfile1;
        tfile1.open(Optfile());
        tfile1.close();

        std::ofstream tfile2;
        tfile2.open(Trjfile());
        tfile2.close();
    }
*/

    inline std::string Optfile() const { return std::string(Basename() + ".opt.xyz"); }
    inline std::string Trjfile() const { return std::string(Basename() + ".trj.xyz"); }
    void start() override; // TODO make pure virtual and move all main action here

    void setSinglePoint(bool sp) { m_singlepoint = sp; }
    inline const std::vector<Molecule>* Molecules() const { return &m_molecules; }

    static Molecule LBFGSOptimise(Molecule* host, const json& controller, std::string& output, std::vector<Molecule>* intermediate, int thread = -1, const std::string& basename = "base");
    static double SinglePoint(const Molecule* initial, const json& controller, std::string& output);

    void clear();

private:
    /* Lets have this for all modules */
    inline nlohmann::json WriteRestartInformation() override { return json(); }

    /* Lets have this for all modules */
    inline bool LoadRestartInformation() override { return true; }

    inline StringList MethodName() const override { return { std::string("opt"), std::string("sp") }; }

    /* Lets have all methods read the input/control file */
    void ReadControlFile() override{};

    /* Read Controller has to be implemented for all */
    void LoadControlJson() override;

    void ProcessMolecules(const std::vector<Molecule>& molecule);
    void ProcessMoleculesSerial(const std::vector<Molecule>& molecule);

    std::string m_filename;
    std::string m_method = "UFF";
    Molecule m_molecule;
    std::vector<Molecule> m_molecules;
    bool m_file_set = false, m_mol_set = false, m_mols_set = false, m_writeXYZ = true, m_printoutput = true, m_singlepoint = false;
    int m_hessian = 0;
    int m_threads = 1;
    double m_dE = 0.1, m_dRMSD = 0.01;
    int m_charge = 0, m_spin = 0;
    int m_serial = false;
};
