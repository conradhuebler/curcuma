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
#include "src/core/parameter_macros.h"
#include "src/core/parameter_registry.h"

#include "src/capabilities/optimisation/lbfgs.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <LBFGS.h>
#include <LBFGSB.h>

#include "curcumamethod.h"

class CurcumaOpt;

// Claude Generated (October 2025): Replace static JSON with ParameterRegistry
// This ensures consistency with the parameter definitions above and enables
// automatic help generation, validation, and JSON export/import capabilities
inline json getCurcumaOptJson() {
    return ParameterRegistry::getInstance().getDefaultJson("opt");
}

// Legacy static reference for backward compatibility - Claude Generated
// NOTE: This creates the JSON on first access, so it's safe for initialization order
#define CurcumaOptJson getCurcumaOptJson()

const json OptJsonPrivate{
    { "writeXYZ", false },
    { "printOutput", true },
    { "dE", 0.1 },
    { "dRMSD", 0.01 },
    { "GFN", 2 },
    { "MaxIter", 2000 },
    { "LBFGS_eps", 1e-5 }
};

// ============================================================================
// Parameter Registry Integration - Claude Generated (October 2025)
// ============================================================================
namespace {
    BEGIN_PARAMETER_DEFINITION(opt)

    // Basic optimization options
    PARAM(write_xyz, Bool, true, "Write XYZ trajectory files", "Basic", { "writeXYZ" })
    PARAM(print_output, Bool, true, "Print optimization progress", "Basic", { "printOutput" })
    PARAM(method, String, "uff", "Computational method for energy/gradient", "Basic", {})
    PARAM(threads, Int, 1, "Number of parallel threads", "Basic", {})
    PARAM(single_point, Bool, false, "Single point energy calculation only", "Basic", { "SinglePoint" })
    PARAM(serial, Bool, false, "Force serial execution", "Basic", {})

    // Convergence criteria
    PARAM(max_iter, Int, 5000, "Maximum optimization iterations", "Convergence", { "MaxIter" })
    PARAM(d_e, Double, 0.1, "Energy change threshold [kJ/mol]", "Convergence", { "dE" })
    PARAM(d_rmsd, Double, 0.01, "RMSD change threshold [Angstrom]", "Convergence", { "dRMSD" })
    PARAM(grad_norm, Double, 5e-4, "Gradient norm convergence threshold", "Convergence", { "GradNorm" })
    PARAM(conv_count, Int, 11, "Number of consecutive converged steps", "Convergence", { "ConvCount" })
    PARAM(single_step, Int, 1, "Single step mode", "Convergence", { "SingleStep" })

    // LBFGS optimizer parameters
    PARAM(lbfgs_m, Int, 2000, "LBFGS memory size", "LBFGS", { "LBFGS_m" })
    PARAM(lbfgs_past, Int, 0, "LBFGS past iterations for delta convergence", "LBFGS", { "LBFGS_past" })
    PARAM(lbfgs_eps_abs, Double, 1e-5, "LBFGS absolute epsilon", "LBFGS", { "LBFGS_eps_abs" })
    PARAM(lbfgs_eps_rel, Double, 1e-5, "LBFGS relative epsilon", "LBFGS", { "LBFGS_eps_rel" })
    PARAM(lbfgs_delta, Double, 0.0, "LBFGS delta for convergence", "LBFGS", { "LBFGS_delta" })
    PARAM(lbfgs_lst, Int, 3, "LBFGS line search type", "LBFGS", { "LBFGS_LST" })
    PARAM(lbfgs_ls_iter, Int, 2, "LBFGS line search max iterations", "LBFGS", { "LBFGS_ls_iter" })
    PARAM(lbfgs_min_step, Double, 1e-4, "LBFGS minimum step size", "LBFGS", { "LBFGS_min_step" })
    PARAM(lbfgs_ftol, Double, 1e-4, "LBFGS function tolerance", "LBFGS", { "LBFGS_ftol" })
    PARAM(lbfgs_wolfe, Double, 0.9, "LBFGS Wolfe condition parameter", "LBFGS", { "LBFGS_wolfe" })

    // Molecular properties
    PARAM(charge, Int, 0, "Total molecular charge", "Molecular", { "Charge" })
    PARAM(spin, Int, 0, "Spin multiplicity", "Molecular", { "Spin" })

    // Advanced optimization options
    PARAM(opt_h, Bool, false, "Optimize only hydrogen atoms", "Advanced", { "optH" })
    PARAM(hessian, Int, 0, "Hessian calculation: 0=none, 1=numerical, 2=analytical", "Advanced", {})
    PARAM(fusion, Bool, false, "Enable fusion optimization", "Advanced", {})
    PARAM(max_rise, Int, 100, "Maximum energy rise allowed", "Advanced", { "maxrise" })
    PARAM(opti_method, Int, 0, "Optimization method: 0=LBFGS external, 1=internal", "Advanced", { "optimethod" })
    PARAM(init_hess, Bool, false, "Calculate initial Hessian", "Advanced", { "inithess" })
    PARAM(lambda, Double, 0.1, "Lambda parameter for optimization", "Advanced", {})
    PARAM(diis_hist, Int, 5, "DIIS history size", "Advanced", {})
    PARAM(diis_start, Int, 5, "DIIS start iteration", "Advanced", {})

    // Molecular orbital output
    PARAM(mo_scheme, Bool, false, "Generate MO occupation scheme", "MO Output", {})
    PARAM(mo_scale, Double, 1.0, "MO visualization scale factor", "MO Output", {})
    PARAM(mo_homo, Int, -1, "HOMO orbital index for output", "MO Output", {})
    PARAM(mo_lumo, Int, -1, "LUMO orbital index for output", "MO Output", {})

    END_PARAMETER_DEFINITION
}
// ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^

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
        setFile(filename);

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
