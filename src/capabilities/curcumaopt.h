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

#pragma once

#include "src/core/energycalculator.h"

#include "src/capabilities/optimiser/lbfgs.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <LBFGS.h>
#include <LBFGSB.h>

#include "curcumamethod.h"

class CurcumaOpt;

static json CurcumaOptJson{
    { "writeXYZ", true },
    { "printOutput", true },
    { "dE", 0.1 },
    { "dRMSD", 0.01 },
    { "method", "uff" },
    { "MaxIter", 5000 },
    //   { "LBFGS_LS", 3},
    { "LBFGS_m", 2000 },
    { "LBFGS_past", 0 },
    { "LBFGS_eps_abs", 1e-5 },
    { "LBFGS_eps_rel", 1e-5 },
    { "LBFGS_delta", 0 },
    { "LBFGS_LST", 3 },
    { "LBFGS_ls_iter", 2 },
    { "LBFGS_min_step", 1e-4 },
    { "LBFGS_ftol", 1e-4 },
    { "LBFGS_wolfe", 0.9 },
    { "SingleStep", 1 },
    { "ConvCount", 11 },
    { "GradNorm", 5e-4 },
    { "threads", 1 },
    { "Charge", 0 },
    { "Spin", 0 },
    { "SinglePoint", false },
    { "optH", false },
    { "serial", false },
    { "hessian", 0 },
    { "fusion", false },
    { "maxrise", 100 },
    { "optimethod", 0 },
    { "inithess", false },
    { "lambda", 0.1 },
    { "diis_hist", 5 },
    { "diis_start", 5 },
    { "mo_scheme", false },
    { "mo_scale", 1.0 },
    { "mo_homo", -1 },
    { "mo_lumo", -1 }

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

using Eigen::VectorXd;
using namespace LBFGSpp;

class LBFGSInterface {
public:
    LBFGSInterface(int n_)
    //: n(n_)
    {
    }
    double operator()(const VectorXd& x, VectorXd& grad);

    double LastEnergy() const { return m_energy; }

    double m_energy = 0, m_last_change = 0, m_last_rmsd = 0;
    Vector Parameter() const { return m_parameter; }
    void setMolecule(const Molecule* molecule);

    void setConstrains(const std::vector<int> constrains) { m_constrains = constrains; }
    void setInterface(EnergyCalculator* interface) { m_interface = interface; }
    void setMethod(int method) { m_method = method; }
    bool isError() const { return m_error; }

private:
    //  int m_iter = 0;
    int m_atoms = 0;
    //  int n;
    int m_method = 2;
    std::vector<int> m_constrains;
    EnergyCalculator* m_interface;
    Vector m_parameter;
    const Molecule* m_molecule;
    bool m_error = false;
};

class SPThread : public CxxThread {
public:
    SPThread(CurcumaOpt* curcumaOpt)
        : m_curcumaOpt(curcumaOpt)
    {
    }
    ~SPThread() = default;

    inline void setMolecule(const Molecule& molecule) { m_molecule = molecule; }
    inline Molecule getMolecule() const { return m_final; }
    virtual int execute() override;

    //  inline void setController(const json& controller) { m_controller = controller; }
    std::string Output() const { return m_result; }
    const std::vector<Molecule>* Intermediates() const { return &m_intermediate; }
    void setBaseName(const std::string& basename) { m_basename = basename; }
    std::string Basename() const { return m_basename; }
    // inline json Parameter() const { return m_param; }
    inline json SCF() const { return m_scf; }
    void setOptiMethod(int method) { m_optimethod = method; }

protected:
    std::string m_result;
    Molecule m_molecule, m_final;
    json m_scf;
    std::vector<Molecule> m_intermediate;
    std::string m_basename;
    CurcumaOpt* m_curcumaOpt;
    int m_optimethod = 0;
};

class OptThread : public SPThread {
public:
    OptThread(CurcumaOpt* curcumaOpt)
        : SPThread(curcumaOpt)
    {
    }
    ~OptThread() = default;

    int execute() override;

private:
};

class OptMThread : public SPThread {
public:
    OptMThread(CurcumaOpt* curcumaOpt)
        : SPThread(curcumaOpt)
    {
    }
    ~OptMThread() = default;

    int execute() override;
    inline void setMolecules(const std::vector<Molecule>& molecules) { m_molecules = molecules; }
    std::vector<Molecule> Molecules() const { return m_finals; }

private:
    std::vector<Molecule> m_molecules, m_finals;
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

    Molecule LBFGSOptimise(Molecule* host, std::string& output, std::vector<Molecule>* intermediate, Vector& charges, int thread = -1, const std::string& basename = "base");
    Molecule GPTLBFGS(Molecule* host, std::string& output, std::vector<Molecule>* intermediate, Vector& charges, int thread = -1, const std::string& basename = "base");

    double SinglePoint(const Molecule* initial, std::string& output, Vector& charges);

    void clear();
    void WriteMO(int n, int m);
    void WriteMOAscii();

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
    json m_parameters;
    std::vector<Molecule> m_molecules;
    Vector m_orbital_energies;
    Matrix m_molecular_orbitals;

    bool m_file_set = false, m_mol_set = false, m_mols_set = false, m_writeXYZ = true, m_printoutput = true, m_singlepoint = false, m_fusion = false, m_optH = false, m_inithess = false, m_mo_scheme = false;
    int m_hessian = 0, m_num_electrons = 0, m_mo_homo = -1, m_mo_lumo = -1;
    int m_threads = 1;
    int m_ConvCount = 1;
    double m_dE = 0.1, m_dRMSD = 0.01, m_maxenergy = 100, m_GradNorm = 1e-5, m_lambda = 0.1, m_mo_scale = 1.0;
    int m_charge = 0, m_spin = 0;
    int m_serial = false;
    int m_maxiter = 100, m_maxrise = 10, m_optimethod = 1, m_diis_hist = 10, m_diis_start = 10;
};
