/*
 * <RMSD Test application within curcuma.>
 * Copyright (C) 2019 - 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/core/molecule.h"

#include "src/capabilities/rmsd.h"

#include "src/tools/general.h"

#include <iostream>
#include <string>

#include "json.hpp"
using json = nlohmann::json;

int main(int argc, char** argv)
{
    int threads = MaxThreads();

    Molecule m1("input_aa.xyz");
    Molecule m2("input_ab.xyz");

    json controller = RMSDJson;
    controller["threads"] = threads;
    RMSDDriver* driver = new RMSDDriver(controller, false);
    driver->setReference(m1);
    driver->setTarget(m2);
    driver->start();
    if (abs(driver->RMSD() - 2.05127) < 1e-5) {
        std::cout << "RMSD calculation with reordering passed (" << driver->RMSD() << ")." << std::endl;
        return 0;
    } else {
        std::cout << "RMSD calculation with reordering failed (" << driver->RMSD() << ")." << std::endl;
        return -1;
    }
}
