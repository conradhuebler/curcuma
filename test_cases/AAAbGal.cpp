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

int AAAbGal_dtemplate() // template
{
    int threads = MaxThreads();

    Molecule m1("A.xyz");
    Molecule m2("B.xyz");

    json controller = RMSDJson;
    controller["threads"] = threads;
    controller["reorder"] = true;
    controller["method"] = "dtemplate";
    RMSDDriver* driver = new RMSDDriver(controller, false);
    driver->setReference(m1);
    driver->setTarget(m2);
    driver->start();
    std::cout << driver->RMSD() << std::endl;
    if (abs(driver->RMSD() - 0.457061) < 1e-5) {
        std::cout << "RMSD calculation with reordering passed (" << driver->RMSD() << ")." << std::endl;
        return 0;
    } else {
        std::cout << "RMSD calculation with reordering failed (" << driver->RMSD() << ")." << std::endl;
        return -1;
    }
}

int AAAbGal_free() // hybrid
{
    int threads = MaxThreads();

    Molecule m1("A.xyz");
    Molecule m2("B.xyz");

    json controller = RMSDJson;
    controller["threads"] = threads;
    controller["reorder"] = true;
    controller["method"] = "free";
    RMSDDriver* driver = new RMSDDriver(controller, false);
    driver->setReference(m1);
    driver->setTarget(m2);
    driver->start();
    if (abs(driver->RMSD() - 0.457061) < 1e-5) {
        std::cout << "RMSD calculation with reordering passed (" << driver->RMSD() << ")." << std::endl;
        return EXIT_SUCCESS;
    } else {
        std::cout << "RMSD calculation with reordering failed (" << driver->RMSD() << ")." << std::endl;
        return EXIT_FAILURE;
    }
}

int AAAbGal_template() // template
{
    int threads = MaxThreads();

    Molecule m1("A.xyz");
    Molecule m2("B.xyz");

    json controller = RMSDJson;
    controller["threads"] = threads;
    controller["reorder"] = true;
    controller["method"] = "template";
    controller["nomunkres"] = true;
    RMSDDriver* driver = new RMSDDriver(controller, false);
    driver->setReference(m1);
    driver->setTarget(m2);
    driver->start();
    if (abs(driver->RMSD() - 0.457061) < 1e-5) {
        std::cout << "RMSD calculation with reordering passed (" << driver->RMSD() << ")." << std::endl;
        return 0;
    } else {
        std::cout << "RMSD calculation with reordering failed (" << driver->RMSD() << ")." << std::endl;
        return -1;
    }
}

int AAAbGal_subspace() // subspace
{
    int threads = MaxThreads();

    Molecule m1("A.xyz");
    Molecule m2("B.xyz");

    json controller = RMSDJson;
    controller["threads"] = threads;
    controller["reorder"] = true;
    controller["method"] = "subspace";
    controller["nomunkres"] = true;

    RMSDDriver* driver = new RMSDDriver(controller, false);
    driver->setReference(m1);
    driver->setTarget(m2);
    driver->start();
    if (abs(driver->RMSD() - 0.457061) < 1e-5) {
        std::cout << "RMSD calculation with reordering passed (" << driver->RMSD() << ")." << std::endl;
        return EXIT_SUCCESS;
    } else {
        std::cout << "RMSD calculation with reordering failed (" << driver->RMSD() << ")." << std::endl;
        return EXIT_FAILURE;
    }
}

int AAAbGal_incr() // incremental
{
    int threads = MaxThreads();

    Molecule m1("A.xyz");
    Molecule m2("B.xyz");

    json controller = RMSDJson;
    controller["threads"] = threads;
    controller["reorder"] = true;
    controller["method"] = "incr";
    RMSDDriver* driver = new RMSDDriver(controller, false);
    driver->setReference(m1);
    driver->setTarget(m2);
    driver->start();
    if (abs(driver->RMSD() - 0.457061) < 1e-5) {
        std::cout << "RMSD calculation with reordering passed (" << driver->RMSD() << ")." << std::endl;
        return EXIT_SUCCESS;
    } else {
        std::cout << "RMSD calculation with reordering failed (" << driver->RMSD() << ")." << std::endl;
        return EXIT_FAILURE;
    }
}


int main(int argc, char** argv)
{
    if(argc == 1)
        return EXIT_FAILURE;
    if(std::string(argv[1]).compare("template") == 0)
        return AAAbGal_template();
    else if (std::string(argv[1]).compare("subspace") == 0)
        return AAAbGal_subspace();
    else if(std::string(argv[1]).compare("incr") == 0)
        return AAAbGal_incr();
    else if (std::string(argv[1]).compare("free") == 0)
        return AAAbGal_free();
    else if (std::string(argv[1]).compare("dtemplate") == 0)
        return AAAbGal_dtemplate();
}
