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

#include "src/capabilities/optimiser/LevMarQMDFFFit.h"

#include "src/capabilities/hessian.h"

#include "src/core/energycalculator.h"
#include "src/core/forcefieldgenerator.h"
#include "src/core/qmdff_par.h"
#include "src/core/topology.h"
#include "src/core/uff_par.h"

#include "qmdfffit.h"

QMDFFFit::QMDFFFit(const json& controller, bool silent)
    : CurcumaMethod(QMDFFFitJson, controller, silent)
{
    UpdateController(controller);
}

void QMDFFFit::LoadControlJson()
{
    m_method = Json2KeyWord<std::string>(m_defaults, "method");
    m_threads = Json2KeyWord<int>(m_defaults, "threads");
}

void QMDFFFit::start()
{
    std::cout << "Parametrising QMDFF (see S. Grimmme, J. Chem. Theory Comput. 2014, 10, 10, 4497–4514 [10.1021/ct500573f]) for the original publication!" << std::endl;
    std::cout << "Starting with the hessian ..." << std::endl;
    std::cout << m_defaults << m_controller << std::endl;
    Hessian hessian(m_defaults, false);
    hessian.setMolecule(m_molecule);
    m_atom_types = m_molecule.Atoms();
    m_geometry = m_molecule.getGeometry();
    hessian.start();
    m_hessian = hessian.getHessian();
    // hessian.PrintVibrations();
    //  std::cout << m_hessian << std::endl;
    // Initialise();
    ForceFieldGenerator ff(m_controller);
    ff.setMolecule(m_molecule);
    ff.Generate();
    json parameter = ff.getParameter();
    // parameter["bonds"] = Bonds();
    // parameter["angles"] = Angles();
    std::ofstream parameterfile_init("qmdff_init_param.json");
    parameterfile_init << parameter;
    // std::cout << parameter << std::endl;
    parameterfile_init.close();
    json qmdff_init = MergeJson(m_controller, QMDFFFitJson);
    qmdff_init["threads"] = m_threads;
    qmdff_init["variable"] = false;
    qmdff_init["method"] = "qmdff";
    EnergyCalculator calculator("qmdff", qmdff_init);
    calculator.setMolecule(m_molecule);
    calculator.setParameter(parameter);
    calculator.CalculateEnergy(false, false);
    Hessian const_hessian("qmdff", qmdff_init);
    json bonds = parameter["bonds"];
    json angles = parameter["angles"];

    /*
    parameter["bonds"] = bonds;
    parameter["angles"] = angles;
    */
    const_hessian.setMolecule(m_molecule);
    const_hessian.setParameter(parameter);
    const_hessian.start();
    Matrix const_hessian_matrix = const_hessian.getHessian();
    for (int i = 0; i < const_hessian_matrix.cols(); ++i)
        for (int j = i; j < const_hessian_matrix.cols(); ++j) {
            // if (std::isnan(const_hessian_matrix(i, j)))
            {
                const_hessian_matrix(i, j) = 0;
                const_hessian_matrix(j, i) = 0;
            }
        }
    // std::cout << const_hessian_matrix << std::endl;
    int counter = 1;
    for (int start = 0; start < 10 && counter; ++start) {
        std::cout << bonds.size() << " " << angles.size() << std::endl;
        std::cout << bonds << std::endl;
        parameter["bonds"] = bonds;
        parameter["angles"] = angles;
        std::ofstream parameterfile_init("qmdff_" + std::to_string(start) + "_param.json");
        parameterfile_init << parameter;
        // std::cout << parameter << std::endl;
        parameterfile_init.close();
        m_fc_parameter.resize(bonds.size() + angles.size());
        int i = 0;
        for (const auto& b : bonds) {
            m_fc_parameter(i) = b["fc"];
            ++i;
        }
        for (const auto& b : angles) {
            m_fc_parameter(i) = b["fc"];
            ++i;
        }
        qmdff_init["variable"] = true;
        qmdff_init["const"] = false;
        std::cout << m_fc_parameter << std::endl;

        Vector vec = OptimiseFC(m_molecule, m_hessian, const_hessian_matrix, m_fc_parameter, parameter, qmdff_init);
        std::cout << vec << std::endl;
        bonds = parameter["bonds"];
        angles = parameter["angles"];
        int index = 0;
        for (int i = 0; i < bonds.size(); ++i) {
            bonds[i]["fc"] = vec(index);
            // m_qmdffbonds[i].kAB = vec(index);
            index++;
        }
        for (int i = 0; i < angles.size(); ++i) {
            angles[i]["fc"] = vec(index);
            // m_qmdffangle[i].kabc = vec(index);
            index++;
        }
        json hc = HessianJson;
        hc["method"] = "qmdff";
        hc["threads"] = m_threads;

        auto cache = bonds;

        bonds.clear();
        counter = 0;
        for (auto c : cache) {
            if ((c["fc"] > 0 && c["distance"] == 1) || c["distance"] == 0)
                bonds.push_back(c);
            else {
                // c["fc"] = 0;
                // bonds.push_back(c);
                counter++;
            }
        }
        auto cache2 = angles;
        angles.clear();
        for (auto c : cache2) {
            if (c["fc"] > 0)
                angles.push_back(c);
            else {
                // c["fc"] = 0;
                // angles.push_back(c);
                counter++;
            }
        }
        std::cout << "Skipping " << counter << " force constants " << std::endl;
        std::cout << bonds.size() << " " << angles.size() << std::endl;

        // std::cout << he2.getHessian() << std::endl;
    }
    /*
        std::cout << "Eigenvectors" << std::endl;
        std::cout << hessian.Frequencies().transpose() << std::endl;
        qmdff_init["variable"] = true;
        qmdff_init["const"] = true;
    */
    Hessian he2(qmdff_init, false);

    parameter["bonds"] = bonds;
    parameter["angles"] = angles;
    he2.setMolecule(m_molecule);
    he2.setParameter(parameter);
    he2.start();
    std::cout << he2.Frequencies().transpose() << std::endl;

    // std::cout << parameter << std::endl;
    std::ofstream parameterfile("qmdff_param.json");
    parameterfile << parameter;
}
