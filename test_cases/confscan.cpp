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

#include "src/capabilities/confscan.h"
#include "src/capabilities/rmsd.h"
#include "src/tools/general.h"

#include <iostream>
#include <string>

#include "json.hpp"
using json = nlohmann::json;

int free()
{
    int threads = MaxThreads();

    json controller = ConfScanJson;
    controller["method"] = "free";
    controller["threads"] = threads;
    controller["silent"] = true;
    controller["restart"] = false;

    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    if (accepted == 14 && reorder_success == 5 && reuse_count == 1 && skipped_count == 237)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}

int subspace()
{
    int threads = MaxThreads();

    json controller = ConfScanJson;
    controller["method"] = "subspace";
    controller["threads"] = threads;
    controller["silent"] = true;
    controller["restart"] = false;

    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << " " << std::endl;

    if (accepted == 14 && reorder_success == 5 && reuse_count == 1 && skipped_count == 237)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}

int template_method()
{
    int threads = MaxThreads();

    json controller = ConfScanJson;
    controller["method"] = "template";
    controller["threads"] = threads;
    controller["silent"] = true;
    controller["restart"] = false;

    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();

    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << " " << std::endl;

    if (accepted == 15 && reorder_success == 4 && reuse_count == 1 && skipped_count == 246)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}

int dtemplate()
{
    int threads = MaxThreads();

    json controller = ConfScanJson;
    controller["method"] = "dtemplate";
    controller["threads"] = threads;
    controller["silent"] = true;
    controller["restart"] = false;

    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << " " << std::endl;
    if (accepted == 15 && reorder_success == 4 && reuse_count == 1 && skipped_count == 248)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}

int molalign()
{
    int threads = 1;

    json controller = ConfScanJson;
    controller["method"] = "molalign";
    controller["threads"] = threads;
    controller["silent"] = true;
    controller["restart"] = false;

    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << " " << std::endl;
    if (accepted == 15 && reorder_success == 5 && reuse_count == 0 && skipped_count == 237)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}

int sLX1()
{
    int threads = MaxThreads();

    json controller = ConfScanJson;
    controller["method"] = "subspace";
    controller["threads"] = threads;
    controller["silent"] = true;
    controller["restart"] = false;
    controller["sLX"] = "1.0";

    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << " " << std::endl;
    if (accepted == 16 && reorder_success == 3 && reuse_count == 1 && skipped_count == 138)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}

int sLX2()
{
    int threads = MaxThreads();

    json controller = ConfScanJson;
    controller["method"] = "subspace";
    controller["threads"] = threads;
    controller["silent"] = true;
    controller["restart"] = false;
    controller["sLX"] = "2.0";

    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << " " << std::endl;
    if (accepted == 14 && reorder_success == 5 && reuse_count == 1 && skipped_count == 101)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}

int sLX2Reset()
{
    int threads = MaxThreads();

    json controller = ConfScanJson;
    controller["method"] = "subspace";
    controller["threads"] = threads;
    controller["silent"] = true;
    controller["restart"] = false;
    controller["sLX"] = "2.0";
    controller["reset"] = true;
    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << " " << std::endl;
    if (accepted == 16 && reorder_success == 5 && reuse_count == 10 && skipped_count == 101)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}

int sLX20()
{
    int threads = MaxThreads();

    json controller = ConfScanJson;
    controller["method"] = "subspace";
    controller["threads"] = threads;
    controller["silent"] = true;
    controller["restart"] = false;
    controller["sLX"] = "2.0";
    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << " " << std::endl;
    if (accepted == 15 && reorder_success == 4 && reuse_count == 0 && skipped_count == 0)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}

int sLX20Reset()
{
    int threads = MaxThreads();

    json controller = ConfScanJson;
    controller["method"] = "subspace";
    controller["threads"] = threads;
    controller["silent"] = true;
    controller["restart"] = false;
    controller["sLX"] = "2.0";
    controller["reset"] = true;
    ConfScan* confscan = new ConfScan(controller);
    confscan->setFileName("input.xyz");
    confscan->start();
    int accepted = confscan->AcceptedCount();
    int reorder_success = confscan->ReorderSuccessfull();
    int reuse_count = confscan->ReuseCount();
    int skipped_count = confscan->ReorderSkippedCount();
    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << " " << std::endl;
    if (accepted == 16 && reorder_success == 4 && reuse_count == 6 && skipped_count == 0)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
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
