/*
 * <Simple MD Module for Cucuma. >
 * Copyright (C) 2020 - 2022 Conrad Hübler <Conrad.Huebler@gmx.net>
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
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "src/capabilities/curcumaopt.h"
#include "src/capabilities/rmsdtraj.h"

#include "src/core/elements.h"
#include "src/core/global.h"
#include "src/core/molecule.h"
#include "src/core/tbliteinterface.h"

#include "src/tools/geometry.h"

#include "simplemd.h"

SimpleMD::SimpleMD(const json& controller, bool silent)
    : CurcumaMethod(CurcumaMDJson, controller, silent)
{
    UpdateController(controller);
    m_interface = new TBLiteInterface;
}

SimpleMD::~SimpleMD()
{
}

void SimpleMD::LoadControlJson()
{
    m_gfn = Json2KeyWord<int>(m_defaults, "GFN");
    m_spin = Json2KeyWord<int>(m_defaults, "spin");
    m_charge = Json2KeyWord<int>(m_defaults, "charge");
    m_single_step = Json2KeyWord<int>(m_defaults, "dT"); // * fs2amu;
    m_maxsteps = Json2KeyWord<int>(m_defaults, "MaxSteps"); // 1e4 * fs2amu;
    m_temperatur = Json2KeyWord<double>(m_defaults, "T");
    m_centered = Json2KeyWord<bool>(m_defaults, "centered");
    m_dumb = Json2KeyWord<int>(m_defaults, "dump");
    m_print = Json2KeyWord<int>(m_defaults, "print");
    m_timestep = m_single_step * fs2amu / 100.0;
    m_berendson = m_timestep; // Automatically choose Velocity Rescaling

    m_rmsd = Json2KeyWord<double>(m_defaults, "rmsd");
    m_opt = Json2KeyWord<bool>(m_defaults, "opt");
    m_scale_velo = Json2KeyWord<double>(m_defaults, "velo");
    m_rescue = Json2KeyWord<bool>(m_defaults, "rescue");
    m_thermostat_steps = Json2KeyWord<int>(m_defaults, "thermostat_steps");
    m_berendson = Json2KeyWord<double>(m_defaults, "berendson") * m_timestep;
    if (m_berendson < m_timestep)
        m_berendson = m_timestep;
    m_writeXYZ = Json2KeyWord<bool>(m_defaults, "writeXYZ");
}

bool SimpleMD::Initialise()
{
    LoadRestartInformation();
    if (m_molecule.AtomCount() == 0)
        return false;

    if (!m_restart) {
        std::ofstream result_file;
        result_file.open(m_basename + ".trj.xyz");
        result_file.close();
    }
    m_natoms = m_molecule.AtomCount();
    m_molecule.setCharge(0);

    m_mass = std::vector<double>(3 * m_natoms, 0);
    m_atomtype = std::vector<int>(m_natoms, 0);

    srand(time(NULL));
    if (!m_restart) {
        m_current_geometry = std::vector<double>(3 * m_natoms, 0);
        m_velocities = std::vector<double>(3 * m_natoms, 0);
        m_currentStep = 0;
    } else {
        //   std::cout << Tools::DoubleVector2String(m_current_geometry)  << std::endl;
        //   std::cout << Tools::DoubleVector2String(m_velocities)  << std::endl;
    }
    /*
    if(m_opt)
    {
        json js = CurcumaOptJson;
        js["writeXYZ"] = false;
        js["gfn"] = m_gfn;
        CurcumaOpt optimise(js, true);
        optimise.addMolecule(&m_molecule);
        optimise.start();
        auto mol = optimise.Molecules();
        std::cout << mol->size() << std::endl;
        auto molecule = ((*mol)[0]);
        for(int i = 0; i < molecule.AtomCount(); ++i)
            m_current_geometry[i] = molecule.Atom(i).second;
    }
    */
    for (int i = 0; i < m_natoms; ++i) {
        m_atomtype[i] = m_molecule.Atom(i).first;

        if (!m_restart) {
            Position pos = m_molecule.Atom(i).second;
            m_current_geometry[3 * i + 0] = pos(0) / au;
            m_current_geometry[3 * i + 1] = pos(1) / au;
            m_current_geometry[3 * i + 2] = pos(2) / au;
            /*
                        double factor;

                        if (m_atomtype[i] == 1)
                            factor = sqrt(kb * m_temperatur / (Elements::AtomicMass[m_atomtype[i]] * m_hmass * amu2au)) * 2;
                        else
                            factor = sqrt(kb * m_temperatur / (Elements::AtomicMass[m_atomtype[i]] * amu2au)) * 2;

                        m_velocities[3 * i + 0] = ((rand() / double(RAND_MAX) * 2) - 1) * factor;
                        m_velocities[3 * i + 1] = ((rand() / double(RAND_MAX) * 2) - 1) * factor;
                        m_velocities[3 * i + 2] = ((rand() / double(RAND_MAX) * 2) - 1) * factor;*/
        }
        if (m_atomtype[i] == 1) {
            m_mass[3 * i + 0] = Elements::AtomicMass[m_atomtype[i]] * amu2au * m_hmass;
            m_mass[3 * i + 1] = Elements::AtomicMass[m_atomtype[i]] * amu2au * m_hmass;
            m_mass[3 * i + 2] = Elements::AtomicMass[m_atomtype[i]] * amu2au * m_hmass;
        } else {
            m_mass[3 * i + 0] = Elements::AtomicMass[m_atomtype[i]] * amu2au;
            m_mass[3 * i + 1] = Elements::AtomicMass[m_atomtype[i]] * amu2au;
            m_mass[3 * i + 2] = Elements::AtomicMass[m_atomtype[i]] * amu2au;
        }
    }
    if (!m_restart)
        InitVelocities(m_scale_velo);

    m_molecule.setCharge(m_charge);
    m_molecule.setSpin(m_spin);
    m_interface->InitialiseMolecule(m_molecule);

    if (m_writeUnique) {
        json rmsdtraj = RMSDTrajJson;
        rmsdtraj["writeUnique"] = true;
        rmsdtraj["rmsd"] = m_rmsd;
        rmsdtraj["writeRMSD"] = false;
        m_unqiue = new RMSDTraj(rmsdtraj, true);
        m_unqiue->setBaseName(m_basename + ".xyz");
        m_unqiue->Initialise();
    }

    m_initialised = true;
    return true;
}

void SimpleMD::InitVelocities(double scaling)
{
    for (int i = 0; i < m_natoms; ++i) {
        double factor;
        if (m_atomtype[i] == 1)
            factor = sqrt(kb * m_temperatur / (Elements::AtomicMass[m_atomtype[i]] * m_hmass * amu2au)) * scaling;
        else
            factor = sqrt(kb * m_temperatur / (Elements::AtomicMass[m_atomtype[i]] * amu2au)) * scaling;

        m_velocities[3 * i + 0] = ((rand() / double(RAND_MAX) * 2) - 1) * factor;
        m_velocities[3 * i + 1] = ((rand() / double(RAND_MAX) * 2) - 1) * factor;
        m_velocities[3 * i + 2] = ((rand() / double(RAND_MAX) * 2) - 1) * factor;
    }
}

nlohmann::json SimpleMD::WriteRestartInformation()
{
    nlohmann::json restart;
    restart["GFN"] = m_gfn;
    restart["dT"] = m_timestep;
    restart["MaxSteps"] = m_maxsteps;
    restart["T"] = m_temperatur;
    restart["currentStep"] = m_currentStep;
    restart["velocities"] = Tools::DoubleVector2String(m_velocities);
    restart["geometry"] = Tools::DoubleVector2String(m_current_geometry);
    restart["centered"] = m_centered;
    restart["average_T"] = m_aver_Temp;
    restart["average_Epot"] = m_aver_Epot;
    restart["average_Ekin"] = m_aver_Ekin;
    restart["average_Etot"] = m_aver_Etot;
    restart["berendson"] = m_berendson;
    restart["thermostat_steps"] = m_thermostat_steps;
    return restart;
};

bool SimpleMD::LoadRestartInformation()
{
    if (!Restart())
        return false;
    StringList files = RestartFiles();
    int error = 0;
    for (const auto& f : files) {

        std::ifstream file(f);
        json restart;
        try {
            file >> restart;
        } catch (json::type_error& e) {
            error++;
            continue;
        } catch (json::parse_error& e) {
            error++;
            continue;
        }

        json md;
        try {
            md = restart[MethodName()[0]];
        } catch (json::type_error& e) {
            error++;
            continue;
        }
        return LoadRestartInformation(md);
        /*
        std::string geometry, velocities;

        try {
            m_gfn = md["GFN"];
        } catch (json::type_error& e) {
        }
        try {
            m_timestep = md["dT"];
        } catch (json::type_error& e) {
        }
        try {
            m_maxsteps = md["MaxSteps"];
        } catch (json::type_error& e) {
        }
        try {
            m_centered = md["centered"];
        } catch (json::type_error& e) {
        }
        try {
            m_temperatur = md["T"];
        } catch (json::type_error& e) {
        }
        try {
            m_currentStep = md["currentStep"];
        } catch (json::type_error& e) {
        }

        try {
            m_aver_Epot = md["average_Epot"];
        } catch (json::type_error& e) {
        }

        try {
            m_aver_Ekin = md["average_Ekin"];
        } catch (json::type_error& e) {
        }

        try {
            m_aver_Etot = md["average_Etot"];
        } catch (json::type_error& e) {
        }
        try {
            m_aver_Temp = md["average_T"];
        } catch (json::type_error& e) {
        }

        try {
            geometry = md["geometry"];
        } catch (json::type_error& e) {
        }
        try {
            velocities = md["velocities"];
        } catch (json::type_error& e) {
        }
        if (geometry.size()) {
            m_current_geometry = Tools::String2DoubleVec(geometry);
        }
        if (velocities.size()) {
            m_velocities = Tools::String2DoubleVec(velocities);
        }
        m_restart = geometry.size() && velocities.size();
        */
    }
    return true;
};

bool SimpleMD::LoadRestartInformation(const json& state)
{
    std::string geometry, velocities;

    try {
        m_gfn = state["GFN"];
    } catch (json::type_error& e) {
    }
    try {
        m_timestep = state["dT"];
    } catch (json::type_error& e) {
    }
    try {
        m_maxsteps = state["MaxSteps"];
    } catch (json::type_error& e) {
    }
    try {
        m_centered = state["centered"];
    } catch (json::type_error& e) {
    }
    try {
        m_temperatur = state["T"];
    } catch (json::type_error& e) {
    }
    try {
        m_currentStep = state["currentStep"];
    } catch (json::type_error& e) {
    }

    try {
        m_aver_Epot = state["average_Epot"];
    } catch (json::type_error& e) {
    }

    try {
        m_aver_Ekin = state["average_Ekin"];
    } catch (json::type_error& e) {
    }

    try {
        m_aver_Etot = state["average_Etot"];
    } catch (json::type_error& e) {
    }
    try {
        m_aver_Temp = state["average_T"];
    } catch (json::type_error& e) {
    }

    try {
        m_berendson = state["berendson"];
    } catch (json::type_error& e) {
    }

    try {
        m_thermostat_steps = state["thermostat_steps"];
    } catch (json::type_error& e) {
    }

    try {
        geometry = state["geometry"];
    } catch (json::type_error& e) {
    }
    try {
        velocities = state["velocities"];
    } catch (json::type_error& e) {
    }
    if (geometry.size()) {
        m_current_geometry = Tools::String2DoubleVec(geometry);
    }
    if (velocities.size()) {
        m_velocities = Tools::String2DoubleVec(velocities);
    }
    m_restart = geometry.size() && velocities.size();

    return true;
}

void SimpleMD::start()
{
    if (m_initialised == false)
        return;
    auto unix_timestamp = std::chrono::seconds(std::time(NULL));
    m_unix_started = std::chrono::milliseconds(unix_timestamp).count();
    double coord[3 * m_natoms];
    double gradient_prev[3 * m_natoms], gradient_current[3 * m_natoms];
    std::vector<json> states;
    for (int i = 0; i < 3 * m_natoms; ++i) {
        coord[i] = m_current_geometry[i];
        gradient_prev[i] = 0;
        gradient_current[i] = 0;
    }

#ifdef GCC
    //         std::cout << fmt::format("{0: ^{0}} {1: ^{1}} {2: ^{2}} {3: ^{3}} {4: ^{4}}\n", "Step", "Epot", "Ekin", "Etot", "T");
    // std::cout << fmt::format("{1: ^{0}} {1: ^{1}} {1: ^{2}} {1: ^{3}} {1: ^{4}}\n", "", "Eh", "Eh", "Eh", "K");
#else
    std::cout << "Step"
              << "\t"
              << "Epot"
              << "\t"
              << "Ekin"
              << "\t"
              << "Etot"
              << "\t"
              << "T" << std::endl;
    std::cout << "  "
              << "\t"
              << "Eh"
              << "\t"
              << "Eh"
              << "\t"
              << "Eh"
              << "\t"
              << "T" << std::endl;
#endif
    Gradient(coord, gradient_prev);
    for (; m_currentStep < m_maxsteps; ++m_currentStep) {
        if (CheckStop() == true) {
            TriggerWriteRestart();
            return;
        }
        // Gradient(coord, gradient_prev);

        if (m_currentStep % m_dumb == 0) {
            bool write = WriteGeometry();
            //   std::cout << write << std::endl;
            if (write) {
                states.push_back(WriteRestartInformation());
                m_current_rescue = 0;
            } else if (!write && m_rescue && states.size() > (1 - m_current_rescue)) {
                std::cout << "Molecule exploded, resetting to previous state ..." << std::endl;
                LoadRestartInformation(states[states.size() - 1 - m_current_rescue]);
                Geometry geometry = m_molecule.getGeometry();
                for (int i = 0; i < m_natoms; ++i) {
                    geometry(i, 0) = m_current_geometry[3 * i + 0] * au;
                    geometry(i, 1) = m_current_geometry[3 * i + 1] * au;
                    geometry(i, 2) = m_current_geometry[3 * i + 2] * au;
                }
                m_molecule.setGeometry(geometry);
                m_molecule.GetFragments();
                InitVelocities(-1);
                for (int i = 0; i < 3 * m_natoms; ++i) {
                    coord[i] = m_current_geometry[i];
                }
                Gradient(coord, gradient_prev);
                m_Ekin = EKin();
                m_Etot = m_Epot + m_Ekin;
                m_current_rescue++;
                PrintStatus();
            }
        }
        SimpleIntegrator(coord, gradient_prev, gradient_current);
        // VelocityVerlet(coord, gradient_prev, gradient_current);

        m_Ekin = EKin();

        if (m_currentStep % m_print == 0) {
            m_Etot = m_Epot + m_Ekin;
            PrintStatus();
        }
        if (m_currentStep % m_thermostat_steps == 0) {
            Berendson();
        }
        if (m_current_rescue >= m_max_rescue) {
            std::cout << "Nothing really helps" << std::endl;
            break;
        }
    }
}
void SimpleMD::SimpleIntegrator(double* coord, double* grad_prev, double* gradient_current)
{
    UpdatePosition(grad_prev, coord);
    m_Epot = Gradient(coord, gradient_current);

    UpdateVelocities(grad_prev, gradient_current);
}

void SimpleMD::VelocityVerlet(double* coord, double* grad_prev, double* gradient_current)
{
    for (int i = 0; i < 3 * m_natoms; ++i) {
        m_velocities[i] += 0.5 * m_timestep * (grad_prev[i] * grad_prev[i]) / m_mass[i];
    }

    for (int i = 0; i < 3 * m_natoms; ++i) {
        m_current_geometry[i] += m_timestep * m_velocities[i];
        coord[i] = m_current_geometry[i];
    }

    m_Epot = Gradient(coord, gradient_current);

    for (int i = 0; i < 3 * m_natoms; ++i) {
        m_velocities[i] += 0.5 * m_timestep * (gradient_current[i] * gradient_current[i]) / m_mass[i];
        grad_prev[i] = gradient_current[i];
    }
}

void SimpleMD::PrintStatus() const
{
    auto unix_timestamp = std::chrono::seconds(std::time(NULL));

    int current = std::chrono::milliseconds(unix_timestamp).count();
    double duration = (current - m_unix_started) / double(m_currentStep);
    double remaining = (m_maxsteps - m_currentStep) / duration;
    //  double remaining = 0;
#ifdef GCC
    std::cout << fmt::format("{1: ^{0}} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f} {6: ^{0}f} {7: ^{0}f} {8: ^{0}f} {9: ^{0}f} {10: ^{0}f}\n", 15, m_currentStep, m_Epot, m_aver_Epot, m_Ekin, m_aver_Ekin, m_Etot, m_aver_Etot, m_T, m_aver_Temp, remaining);
#else
    std::cout << m_currentStep << " " << m_Epot << " " << m_Ekin << " " << m_Epot + m_Ekin << m_T << std::endl;

#endif
}

void SimpleMD::PrintMatrix(const double* matrix)
{
    std::cout << "Print Matrix" << std::endl;
    for (int i = 0; i < m_natoms; ++i) {
        std::cout << matrix[3 * i] << " " << matrix[3 * i + 1] << " " << matrix[3 * i + 2] << std::endl;
    }
    std::cout << std::endl;
}

void SimpleMD::UpdatePosition(const double* gradient, double* coord)
{
    for (int i = 0; i < 3 * m_natoms; ++i) {
        m_current_geometry[i] += m_timestep * m_velocities[i] + gradient[i] / (2 * m_mass[i]);
        //     m_current_geometry[i] += m_timestep * m_velocities[i]; // + gradient[i] / (2 * m_mass[i]);
        coord[i] = m_current_geometry[i];
    }
}

void SimpleMD::UpdateVelocities(double* gradient_prev, const double* gradient_curr)
{
    for (int i = 0; i < 3 * m_natoms; ++i) {
        m_velocities[i] -= m_timestep * ((gradient_prev[i] + gradient_curr[i]) / (4 * m_mass[i]));
        gradient_prev[i] = gradient_curr[i];
    }
}

double SimpleMD::Gradient(const double* coord, double* grad)
{
    m_interface->UpdateMolecule(coord);

    double Energy = m_interface->GFNCalculation(m_gfn, grad);
    /* for(int i = 0; i < 3 * m_natoms; ++i)
     {
         grad[i] /= au;

     }*/
    return Energy;
}

double SimpleMD::EKin()
{
    double ekin = 0;
    for (int i = 0; i < 3 * m_natoms; ++i) {
        ekin += m_mass[i] * m_velocities[i] * m_velocities[i];
    }
    ekin *= 0.5;
    m_T = 2.0 * ekin / (kb * 3 * m_natoms);
    m_aver_Temp = (m_T + (m_currentStep)*m_aver_Temp) / (m_currentStep + 1);
    m_aver_Epot = (m_Epot + (m_currentStep)*m_aver_Epot) / (m_currentStep + 1);
    m_aver_Ekin = (m_Ekin + (m_currentStep)*m_aver_Ekin) / (m_currentStep + 1);
    m_aver_Etot = (m_Etot + (m_currentStep)*m_aver_Etot) / (m_currentStep + 1);

    return ekin;
}

bool SimpleMD::WriteGeometry()
{
    Geometry geometry = m_molecule.getGeometry();
    for (int i = 0; i < m_natoms; ++i) {
        geometry(i, 0) = m_current_geometry[3 * i + 0] * au;
        geometry(i, 1) = m_current_geometry[3 * i + 1] * au;
        geometry(i, 2) = m_current_geometry[3 * i + 2] * au;
    }
    int f1 = m_molecule.GetFragments().size();
    m_molecule.setGeometry(geometry);

    int f2 = m_molecule.GetFragments().size();
    //  std::cout << f1 << " ... " << f2 << std::endl;
    m_prev_index = std::abs(f2 - f1);
    if (f1 != f2) {
        return false;
    }
    if (m_writeXYZ) {
        m_molecule.setEnergy(m_Epot);
        m_molecule.setName(std::to_string(m_currentStep));
        if (m_centered)
            m_molecule.Center();
        m_molecule.appendXYZFile(m_basename + ".trj.xyz");
    }
    if (m_writeUnique) {
        m_unqiue->CheckMolecule(new Molecule(m_molecule));
    }
    return true;
}

void SimpleMD::Berendson()
{
    double lambda = sqrt(1 + m_timestep / m_berendson * (m_temperatur / m_T - 1));
    for (int i = 0; i < 3 * m_natoms; ++i) {
        m_velocities[i] *= lambda;
    }
}
