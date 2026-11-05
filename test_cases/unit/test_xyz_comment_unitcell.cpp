/*
 * Unit test: unit cell from the XYZ comment line
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Claude Generated (Nov 2026): "unitcell a b c" without three numeric angles must give an
 * orthorhombic cell (90 90 90) instead of no cell; "unitcell abc" gives no cell and a warning
 * (the warning text is checked by the ctest PASS_REGULAR_EXPRESSION).
 * Usage: test_xyz_comment_unitcell <scratch dir>
 */

#include "src/core/curcuma_logger.h"
#include "src/core/molecule.h"
#include "src/tools/formats.h"

#include "test_molecule_registry.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

namespace {

Molecule loadWithComment(const std::string& dir, const std::string& name, const std::string& comment)
{
    // geometry from the structure library (registry H2O, Angstrom); only the comment line varies
    const curcuma::Molecule water = TestMolecules::TestMoleculeRegistry::createMolecule("H2O", false);
    const auto geom = water.getGeometry();
    const auto atoms = water.Atoms();
    const std::string path = dir + "/" + name + ".xyz";
    std::ofstream out(path);
    out << atoms.size() << "\n"
        << comment << "\n";
    for (std::size_t i = 0; i < atoms.size(); ++i)
        out << Elements::ElementAbbr[atoms[i]] << " " << fmt::format("{:.8f} {:.8f} {:.8f}", geom(i, 0), geom(i, 1), geom(i, 2)) << "\n";
    out.close();
    return Files::LoadFile(path);
}

bool expectCell(const std::string& dir, const std::string& name, const std::string& comment)
{
    Molecule mol = loadWithComment(dir, name, comment);
    if (!mol.hasPBC()) {
        std::printf("FAIL [%s]: no cell for '%s'\n", name.c_str(), comment.c_str());
        return false;
    }
    const Eigen::Matrix3d cell = mol.getUnitCell();
    const Eigen::Matrix3d ref = 10.0 * Eigen::Matrix3d::Identity();
    const double dev = (cell - ref).cwiseAbs().maxCoeff();
    if (dev > 1e-10) {
        std::printf("FAIL [%s]: cell deviates from 10x10x10 orthorhombic by %.3e for '%s'\n",
            name.c_str(), dev, comment.c_str());
        return false;
    }
    std::printf("ok   [%s]: 10x10x10, angles 90 for '%s'\n", name.c_str(), comment.c_str());
    return true;
}

bool expectNoCell(const std::string& dir, const std::string& name, const std::string& comment)
{
    Molecule mol = loadWithComment(dir, name, comment);
    if (mol.hasPBC()) {
        std::printf("FAIL [%s]: unexpected cell for '%s'\n", name.c_str(), comment.c_str());
        return false;
    }
    std::printf("ok   [%s]: no cell for '%s'\n", name.c_str(), comment.c_str());
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    const std::string dir = argc > 1 ? argv[1] : ".";
    CurcumaLogger::set_verbosity(1);

    int failed = 0;
    failed += !expectCell(dir, "uc_abc", "unitcell 10 10 10");
    failed += !expectCell(dir, "uc_angles", "unitcell 10 10 10 90 90 90");
    failed += !expectCell(dir, "uc_text", "unitcell 10 10 10 water box");
    failed += !expectCell(dir, "uc_text_long", "unitcell 10 10 10 water box 432");
    failed += !expectCell(dir, "lattice_abc", "lattice 10 10 10");
    failed += !expectNoCell(dir, "uc_bad", "unitcell abc");

    if (failed) {
        std::printf("%d CHECKS FAILED\n", failed);
        return 1;
    }
    std::printf("ALL 6 CHECKS PASSED\n");
    return 0;
}
