/*
 * <RMSD Test application within curcuma.>
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/capabilities/confscan.h"
#include "src/capabilities/rmsd.h"
#include "src/core/parameter_registry.h"  // Claude Generated - For ParameterRegistry access
#include "src/tools/general.h"

#include <iostream>
#include <string>

#include "json.hpp"
using json = nlohmann::json;

// Claude Generated - Initialize ParameterRegistry from auto-generated definitions
#include "generated/parameter_registry.h"
static struct RegistryInitializer {
    RegistryInitializer() { initialize_generated_registry(); }
} registry_initializer;

// Claude Generated (Oct 2026): the expected counters are golden values; on a mismatch print actual and expected,
// otherwise the failure carries no information.
static int checkCounters(const char* name, int accepted, int reorder, int reuse, int skipped,
    int e_accepted, int e_reorder, int e_reuse, int e_skipped)
{
    if (accepted == e_accepted && reorder == e_reorder && reuse == e_reuse && skipped == e_skipped)
        return EXIT_SUCCESS;
    std::cerr << name << ": counters differ from the golden values (actual / expected): accepted " << accepted << " / " << e_accepted
              << ", reorder successful " << reorder << " / " << e_reorder << ", reuse " << reuse << " / " << e_reuse
              << ", reorder skipped " << skipped << " / " << e_skipped << std::endl;
    return EXIT_FAILURE;
}

int free()
{
    int threads = MaxThreads();

    // Claude Generated 2025: Proper Multi-Module JSON for ConfigManager
    // ConfScan constructor expects { "confscan": {...}, "rmsd": {...} }
    json controller;
    controller["confscan"] = ParameterRegistry::getInstance().getDefaultJson("confscan");
    controller["rmsd"] = ParameterRegistry::getInstance().getDefaultJson("rmsd");
    controller["confscan"]["threads"] = threads;
    controller["confscan"]["silent"] = true;
    controller["confscan"]["restart"] = false;
    controller["rmsd"]["method"] = "free";

    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    return checkCounters(__func__, accepted, reorder_success, reuse_count, skipped_count, 14, 5, 1, 237);
}

int subspace()
{
    int threads = MaxThreads();

    // Claude Generated 2025: Correct Multi-Module JSON with overrides only
    json controller;
    controller["confscan"] = ParameterRegistry::getInstance().getDefaultJson("confscan");
    controller["rmsd"] = ParameterRegistry::getInstance().getDefaultJson("rmsd");
    controller["confscan"]["threads"] = threads;
    controller["confscan"]["silent"] = true;
    controller["confscan"]["restart"] = false;
    controller["rmsd"]["method"] = "subspace";

    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << " " << std::endl;

    return checkCounters(__func__, accepted, reorder_success, reuse_count, skipped_count, 14, 5, 1, 237);
}

int template_method()
{
    int threads = MaxThreads();

    // Claude Generated 2025: Correct Multi-Module JSON with overrides only
    json controller;
    controller["confscan"] = ParameterRegistry::getInstance().getDefaultJson("confscan");
    controller["rmsd"] = ParameterRegistry::getInstance().getDefaultJson("rmsd");
    controller["confscan"]["threads"] = threads;
    controller["confscan"]["silent"] = true;
    controller["confscan"]["restart"] = false;
    controller["rmsd"]["method"] = "template";

    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();

    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << " " << std::endl;

    return checkCounters(__func__, accepted, reorder_success, reuse_count, skipped_count, 15, 4, 1, 246);
}

int dtemplate()
{
    int threads = MaxThreads();

    // Claude Generated 2025: Correct Multi-Module JSON with overrides only
    json controller;
    controller["confscan"] = ParameterRegistry::getInstance().getDefaultJson("confscan");
    controller["rmsd"] = ParameterRegistry::getInstance().getDefaultJson("rmsd");
    controller["confscan"]["threads"] = threads;
    controller["confscan"]["silent"] = true;
    controller["confscan"]["restart"] = false;
    controller["rmsd"]["method"] = "dtemplate";

    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    return checkCounters(__func__, accepted, reorder_success, reuse_count, skipped_count, 14, 5, 1, 237);
}

int molalign()
{
    int threads = 1;

    // Claude Generated 2025: Correct Multi-Module JSON with overrides only
    json controller;
    controller["confscan"] = ParameterRegistry::getInstance().getDefaultJson("confscan");
    controller["rmsd"] = ParameterRegistry::getInstance().getDefaultJson("rmsd");
    controller["confscan"]["threads"] = threads;
    controller["confscan"]["silent"] = true;
    controller["confscan"]["restart"] = false;
    controller["rmsd"]["method"] = "molalign";

    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    return checkCounters(__func__, accepted, reorder_success, reuse_count, skipped_count, 20, 0, 0, 349);
}

int sLX1()
{
    int threads = MaxThreads();

    // Claude Generated 2025: Correct Multi-Module JSON with overrides only
    json controller;
    controller["confscan"] = ParameterRegistry::getInstance().getDefaultJson("confscan");
    controller["rmsd"] = ParameterRegistry::getInstance().getDefaultJson("rmsd");
    controller["confscan"]["threads"] = threads;
    controller["confscan"]["silent"] = true;
    controller["confscan"]["restart"] = false;
    controller["confscan"]["slx"] = "1.0";
    controller["rmsd"]["method"] = "subspace";

    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << " " << std::endl;
    return checkCounters(__func__, accepted, reorder_success, reuse_count, skipped_count, 16, 3, 1, 138);
}

int sLX2()
{
    int threads = MaxThreads();

    // Claude Generated 2025: Correct Multi-Module JSON with overrides only
    json controller;
    controller["confscan"] = ParameterRegistry::getInstance().getDefaultJson("confscan");
    controller["rmsd"] = ParameterRegistry::getInstance().getDefaultJson("rmsd");
    controller["confscan"]["threads"] = threads;
    controller["confscan"]["silent"] = true;
    controller["confscan"]["restart"] = false;
    controller["confscan"]["slx"] = "2.0";
    controller["rmsd"]["method"] = "subspace";

    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << " " << std::endl;
    return checkCounters(__func__, accepted, reorder_success, reuse_count, skipped_count, 14, 5, 1, 101);
}

int sLX2Reset()
{
    int threads = MaxThreads();

    // Claude Generated 2025: Correct Multi-Module JSON with overrides only
    json controller;
    controller["confscan"] = ParameterRegistry::getInstance().getDefaultJson("confscan");
    controller["rmsd"] = ParameterRegistry::getInstance().getDefaultJson("rmsd");
    controller["confscan"]["threads"] = threads;
    controller["confscan"]["silent"] = true;
    controller["confscan"]["restart"] = false;
    controller["confscan"]["slx"] = "2.0";
    controller["confscan"]["reset"] = true;
    controller["rmsd"]["method"] = "subspace";
    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << " " << std::endl;
    return checkCounters(__func__, accepted, reorder_success, reuse_count, skipped_count, 16, 5, 10, 101);
}

int sLX20()
{
    int threads = MaxThreads();

    // Claude Generated 2025: Correct Multi-Module JSON with overrides only
    json controller;
    controller["confscan"] = ParameterRegistry::getInstance().getDefaultJson("confscan");
    controller["rmsd"] = ParameterRegistry::getInstance().getDefaultJson("rmsd");
    controller["confscan"]["threads"] = threads;
    controller["confscan"]["silent"] = true;
    controller["confscan"]["restart"] = false;
    controller["confscan"]["slx"] = "2.0";
    controller["rmsd"]["method"] = "subspace";
    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    return checkCounters(__func__, accepted, reorder_success, reuse_count, skipped_count, 14, 5, 1, 101);
}

int sLX20Reset()
{
    int threads = MaxThreads();

    // Claude Generated 2025: Correct Multi-Module JSON with overrides only
    json controller;
    controller["confscan"] = ParameterRegistry::getInstance().getDefaultJson("confscan");
    controller["rmsd"] = ParameterRegistry::getInstance().getDefaultJson("rmsd");
    controller["confscan"]["threads"] = threads;
    controller["confscan"]["silent"] = true;
    controller["confscan"]["restart"] = false;
    controller["confscan"]["slx"] = "2.0";
    controller["confscan"]["reset"] = true;
    controller["rmsd"]["method"] = "subspace";
    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    return checkCounters(__func__, accepted, reorder_success, reuse_count, skipped_count, 16, 5, 10, 101);
}

int main(int argc, char** argv)
{
    if (argc == 1)
        return EXIT_FAILURE;
    if (std::string(argv[1]).compare("free") == 0)
        return free();
    else if (std::string(argv[1]).compare("subspace") == 0)
        return subspace();
    else if (std::string(argv[1]).compare("dtemplate") == 0)
        return dtemplate();
    else if (std::string(argv[1]).compare("template") == 0)
        return template_method();
    else if (std::string(argv[1]).compare("molalign") == 0)
        return molalign();
    else if (std::string(argv[1]).compare("sLX1") == 0)
        return sLX1();
    else if (std::string(argv[1]).compare("sLX2") == 0)
        return sLX2();
    else if (std::string(argv[1]).compare("sLX2Reset") == 0)
        return sLX2Reset();
    else if (std::string(argv[1]).compare("sLX20") == 0)
        return sLX20();
    else if (std::string(argv[1]).compare("sLX20Reset") == 0)
        return sLX20Reset();
}
