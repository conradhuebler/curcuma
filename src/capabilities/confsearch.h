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

#pragma once

#include <limits>
#include <string>
#include <vector>

#include "src/tools/general.h"

#include "src/core/molecule.h"
#include "src/core/intra_parallel_context.h"

#include "json.hpp"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include "src/capabilities/curcumamethod.h"
#include "src/capabilities/shared_bias_pool.h"
#include "src/capabilities/optimizer_factory.h"
#include "src/core/energycalculator.h"
using namespace curcuma;
static const nlohmann::json ConfSearchJson{
    { "method", "gfnff" },
    { "startT", 600 },
    { "endT", 300 },
    { "deltaT", 50 },
    { "repeat", 10 },
    { "time", 1e4 }, // 2 ps
    { "rmsd", 1.25 },
    { "threads", 1 },
    { "energy_window", 100 },
    { "writeXYZ", true },
    { "printOutput", true },
    { "MaxTime", 5000 },
    { "T", 298.15 },
    { "dt", 1 }, // single step in fs
    { "rm_COM", 100 }, // remove translation and rotation every x fs
    { "charge", 0 },
    { "Spin", 0 },
    { "rmrottrans", 0 },
    { "nocenter", false },
    { "dump", 50 },
    { "print", 1000 },
    { "unique", false },
    { "rmsd", 1.5 },
    { "opt", false },
    { "hmass", 1 },
    { "velo", 1 },
    { "rescue", false },
    { "coupling", 10 },
    { "MaxTopoDiff", 15 },
    { "impuls_scaling", 0.75 },
    { "writeinit", false },
    { "initfile", "none" },
    { "norestart", false },
    { "writerestart", 1000 },
    { "rattle", false },
    { "rattle_tolerance", 1e-2 },
    { "rattle_maxiter", 10 },
    { "thermostat", "csvr" },
    { "respa", 1 },
    { "dipole", false },
    { "seed", -1 },
    { "max_bias_export", 1000 },
    // Claude Generated (Jun 2026): ConfSearch efficiency/robustness controls
    { "rattle_threshold_temp", 400 }, // K; at/above this cycle temperature RATTLE is auto-enabled
    { "rattle_hot_mode", 2 }, // RATTLE mode used above the threshold (2 = constrain X-H only)
    { "topo_check", false }, // opt-in: abort an MD run when the molecule fragments (GetFragments grows)
    { "topo_check_interval", 0 }, // steps between topology checks (0 -> use the MD dump frequency)
    { "seed_energy_window", 50 }, // kJ/mol vs. the running global minimum; selects next-cycle MD seeds
    { "seed_window_schedule", "static" }, // "static" or "exp" (funnel: window shrinks each cycle)
    { "seed_window_decay", 0.5 }, // per-cycle multiplier for the exp funnel schedule
    { "epot_abort", false }, // opt-in: abort an MD run when the running-mean potential climbs too high
    { "epot_abort_window", 250 }, // kJ/mol above the run's starting energy (must exceed thermal baseline)
    { "opt_feedback_bias", true }, // deposit optimised minima back into the shared bias pool
    { "opt_feedback_height", 5 }, // hill counter (height = k*counter) assigned to fed-back minima
    { "mtd_permutation", true }, // feed ConfScan's symmetry reorder rules into the RMSD-MTD bias (smooth, sum-over-images)
    { "cleanenergy", false },
    { "wall", "none" }, // can be spheric or rect
    { "wall_type", "logfermi" }, // can be logfermi or harmonic
    { "wall_spheric_radius", 0 },
    { "wall_xl", 0 },
    { "wall_yl", 0 },
    { "wall_zl", 0 },
    { "wall_x_min", 0 },
    { "wall_x_max", 0 },
    { "wall_y_min", 0 },
    { "wall_y_max", 0 },
    { "wall_z_min", 0 },
    { "wall_z_max", 0 },
    { "wall_temp", 298.15 },
    { "wall_beta", 6 }
};

// Claude Generated (Apr 2026): Thread class for parallel geometry optimization
class OptThread : public CxxThread {
public:
    OptThread(const Molecule& molecule, const json& parameter)
        : m_molecule(molecule), m_parameter(parameter)
    {
        setAutoDelete(false);
    }

    virtual int execute() override
    {
        // This OptThread runs as one task among many under a molecule-level pool.
        // Suppress intra-molecule fan-out so methods that honor the flag (native
        // gfn1/gfn2) stay serial and the cores are not oversubscribed.
        curcuma::SuppressIntraParallel intra_guard;

        std::string method = m_parameter.value("method", std::string("gfnff"));
        EnergyCalculator energy_calc(method, m_parameter);
        m_result = Optimization::OptimizationDispatcher::optimizeStructure(
            &m_molecule, Optimization::OptimizerType::LBFGSPP, &energy_calc, m_parameter);
        return m_result.success ? 0 : 1;
    }

    const Optimization::OptimizationResult& result() const { return m_result; }

private:
    Molecule m_molecule;
    json m_parameter;
    Optimization::OptimizationResult m_result;
};

// Claude Generated (May 2026, ICX-build): forward decl must be inside the namespace.
// ICX rejects `class curcuma::Molecule;` as a nested-name forward decl; GCC accepts it.
namespace curcuma { class Molecule; }

class ConfSearch : public CurcumaMethod {
public:
    ConfSearch(const json& controller, bool silent);
    ~ConfSearch();

    void setFile(const std::string& file);
    virtual bool Initialise() override;

    virtual void start() override;

    virtual void printHelp() const override
    {
        std::cout << "Usage: curcuma -confsearch input.xyz [parameters]\n\n";
        std::cout << "Conformational Search Parameters (from ConfSearchJson):\n";
        for (auto& el : ConfSearchJson.items()) {
            std::cout << "  -" << el.key() << "  (default: " << el.value() << ")\n";
        }
        std::cout << "\nRMSD-MTD is enabled by default. Use -rmsd_mtd false to disable.\n";
    }

private:
    void PerformMolecularDynamics(const std::vector<Molecule*>& molecules, const nlohmann::json& parameter);

    std::string PerformOptimisation(const std::string& filename, const nlohmann::json& parameter);

    std::string PerformFilter(const std::string& filename, const nlohmann::json& parameter);

    /* Lets have this for all modules */
    virtual nlohmann::json WriteRestartInformation() override;

    /* Lets have this for all modules */
    virtual bool LoadRestartInformation() override;

    virtual StringList MethodName() const override
    {
        return { "ConfSearch" };
    }

    /* Lets have all methods read the input/control file */
    virtual void ReadControlFile() override;

    /* Read Controller has to be implemented for all */
    virtual void LoadControlJson() override;

    StringList m_error_list;
    std::string m_method, m_thermostat;
    bool m_silent = true, m_rattle = true;
    double m_dT = 4;
    std::vector<Molecule*> m_in_stack, m_final_stack;
    int m_spin = 0, m_charge = 0, m_repeat = 5, m_threads = 1, m_max_bias_export = 1000;
    double m_time = 1e4, m_startT = 500, m_endT = 300, m_deltaT = 50, m_currentT = 0, m_rmsd = 1.25, m_energy_window = 100;
    // Claude Generated (Jun 2026): efficiency/robustness controls (see ConfSearchJson for meaning)
    double m_rattle_threshold_temp = 400, m_seed_energy_window = 50, m_seed_window_decay = 0.5, m_epot_abort_window = 250;
    int m_rattle_hot_mode = 2, m_topo_check_interval = 0, m_opt_feedback_height = 5;
    bool m_topo_check = false, m_epot_abort = false, m_opt_feedback_bias = true, m_mtd_permutation = true;
    std::string m_seed_window_schedule = "static";
    double m_global_min = std::numeric_limits<double>::infinity(); // running lowest energy across all cycles
    std::vector<std::vector<int>> m_permutation_cache; // Claude Generated (Jun 2026): symmetry reorder rules from ConfScan, fed into MTD
    Matrix m_topo_matrix;
    SharedBiasPool* m_bias_pool = nullptr;  // Claude Generated (Apr 2026): shared bias pool for parallel ConfSearch
};
