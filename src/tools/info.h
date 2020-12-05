/*
 * <Geometry tools for chemical structures.>
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

#pragma once

#include <fmt/color.h>
#include <fmt/core.h>

#include "src/version.h"

namespace General {

inline void StartUp(int argc, char** argv)
{
    std::string version_string = "Curcuma - Version:  " + git_tag;
    std::string hastag = "Git Commit Hash: " + git_commit_hash;
    fmt::print(
        "*{0:*^{1}}*\n"
        "*{2: ^{1}}*\n"
        "*{0: ^{1}}*\n"
        "*{3: ^{1}}*\n"
        "*{0: ^{1}}*\n"
        "*{12: ^{1}}*\n"
        "*{0: ^{1}}*\n"
        "*{4: ^{1}}*\n"
        "*{5: ^{1}}*\n"
        "*{0: ^{1}}*\n"
        "*{6: ^{1}}*\n"
        "*{7: ^{1}}*\n"
        "*{0: ^{1}}*\n"
        "*{8: ^{1}}*\n"
        "*{9: ^{1}}*\n"
        "*{10: ^{1}}*\n"
        "*{0: ^{1}}*\n"
        "*{11: ^{1}}*\n"
        "*{0: ^{1}}*\n"
        "*{0:*^{1}}*\n",
        "", 60,
        "Curcuma - Simple Molecular Modelling tool!",
        version_string,
        "Visit the website for initial usage",
        "https://github.com/conradhuebler/curcuma",
        "This program comes without any warranty",
        "It might even be total useless",
        "If you obtain some usefull results with",
        "curcuma, please cite",
        "http://doi.org/10.5281/zenodo.4302722",
        hastag,
        "Written by Conrad Hübler TU Freiberg");

    fmt::print("\nCurcuma was started as:\n");
    for (int index = 0; index < argc; ++index)
        fmt::print(fg(fmt::color::light_green) | fmt::emphasis::bold, "{} ", argv[index]);

    fmt::print("\n\n");
}

}
