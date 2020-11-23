/*
 * <Statistics of conformers >
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

#include "confstat.h"

ConfStat::ConfStat(const json& controller, bool silent)
    : CurcumaMethod(ConfStatJson, controller, silent)
{
    UpdateController(controller);
}

void ConfStat::LoadControlJson()
{
    // m_cutoff = Json2KeyWord<double>(m_defaults, "Cutoff");
    // m_temp = Json2KeyWord<double>(m_defaults, "Temp");
}

void ConfStat::start()
{
    FileIterator file(m_filename);
    while (!file.AtEnd()) {
        Molecule* mol = new Molecule(file.Next());
        double energy = mol->Energy();
        /*if (std::abs(energy) < 1e-5 || m_gfn != -1) {
            XTBInterface interface; // As long as xtb leaks, we have to put it heare
                // I might not leak really, but was unable to clear everything
            if (m_gfn == -1)
                m_gfn = 2;
            interface.InitialiseMolecule(mol);
            energy = interface.GFNCalculation(m_gfn, 0);
            //interface.clear();
        }*/
        m_energies.push_back(energy);
        /*
        m_ordered_list.insert(std::pair<double, int>(energy, molecule));
        molecule++;
        if (m_noname)
            mol->setName(NamePattern(molecule));
        mol->CalculateRotationalConstants();
        std::pair<std::string, Molecule*> pair(mol->Name(), mol);
        m_molecules.push_back(pair);
        */
    }

    if (m_energies.size() == 0)
        return;

    std::sort(m_energies.begin(), m_energies.end());
    double factor = -1 * 1000 * 2625.5 / (R * m_temp);
    double sum = 0;
    double max = 0;
    double CR = 0;
    double consistency = 0.0;
    for (int i = 0; i < m_energies.size(); ++i) {
        double diff = m_energies[i] - m_energies[0];
        double pi = exp((diff)*factor);
        max += pi;
    }
    printf("\n\nLowest energy is %8.5f Eh!\n\n", m_energies[0]);
    printf("   Diff E (in Eh and kJ/mol) |Intermediates    | Population\n\n");
    for (int i = 0; i < m_energies.size(); ++i) {
        double diff = m_energies[i] - m_energies[0];
        double tmp = exp(diff * factor);
        double pi = exp(-1 * (diff)*2625.5 * 1000 / R / m_temp) / max;
        CR -= pi * log(pi);
        sum += tmp;
        if (pi * 100 >= m_print_threshold || i < 10)
            printf("%8.5f %8.2f | %8.5f %8.5f %8.5f | %8.5f\n", diff, (diff)*2625.5, tmp, sum, -R * m_temp * log(1 + sum) / 1000.0, pi * 100);
        // consistency +=  pi/max * 100;
        // std::cout << pi << std::endl;
        //std::cout << diff << " " << (diff)*2625.5 << " " << tmp << " " << sum << " " << -R * m_temp * log(1 + sum) / 1000.0 << " " << exp(-1 * (diff)*2625.5 * 1000 / R / m_temp) / max * 100 << " %" << std::endl;
    }
    std::cout << "Finally add " << -R * m_temp * log(sum) / 1000 << " kJ/mol to the free enthalpy / energy !" << std::endl;
    std::cout << "Scr " << R * CR << " J/mol*K to the entropy = " << -1 * CR * R * m_temp / 1000 << " kJ/mol !" << std::endl;
    // std::cout << consistency << std::endl;
}
