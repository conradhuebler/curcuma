/*
 * <QMDFF Hessian fit for parametrisation. >
 * Copyright (C) 2023 - 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/core/energy_calculators/ff_methods/qmdff_par.h"
#include "src/core/molecule.h"
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

#include "curcumamethod.h"

/* Claude Generated 2025: QMDFFfit Parameter Registry - replaces static QMDFFFitJson */
BEGIN_PARAMETER_DEFINITION(qmdfffit)
    PARAM(method, String, "gfn2", "QM method for initial force field.", "General", {})
    PARAM(hessian_file, String, "hessian.json", "Input Hessian file.", "Input", {"hessian"})
    PARAM(charges_file, String, "scf.json", "Input charges file.", "Input", {"charges"})
    PARAM(potential, String, "uff", "Potential for initial FF (uff or qmdff).", "General", {})
    PARAM(threads, Int, 1, "Number of threads.", "Performance", {})
END_PARAMETER_DEFINITION

class QMDFFFit : public CurcumaMethod {
public:
    /**
     * @brief Constructor with JSON configuration (backward compatible)
     * Claude Generated: Phase 4 - ConfigManager Migration
     */
    QMDFFFit(const json& controller, bool silent = true);

    /**
     * @brief Constructor with ConfigManager configuration (new, preferred)
     * Claude Generated: Phase 4 - Native ConfigManager support
     */
    QMDFFFit(const ConfigManager& config, bool silent = true);
    virtual ~QMDFFFit() {}
    void start() override;

    void setMolecule(const Molecule& molecule) { m_molecule = molecule; }

private:
    // void FitDihedrals(json &dihedrals);

    /* Lets have this for all modules */
    nlohmann::json WriteRestartInformation() override { return json(); }

    /* Lets have this for all modules */
    bool LoadRestartInformation() override { return true; }

    StringList MethodName() const override { return { std::string("QMDFFFit") }; }

    /* Lets have all methods read the input/control file */
    void ReadControlFile() override {}

    /* Read Controller has to be implemented for all */
    void LoadControlJson() override;
    Molecule m_molecule;
    std::string m_method = "gfn2", m_hessian_file, m_scf_file;
    Matrix m_hessian, m_geometry;
    std::vector<int> m_atom_types, m_coordination, m_topo;
    std::vector<std::vector<int>> m_stored_bonds;
    std::vector<std::vector<int>> m_identified_rings;
    Vector m_fc_parameter;
    json m_bonds, m_angles;

    double m_scaling;
    bool m_rings = true;
    double m_au = 1;
    int m_threads = 1;
};
