/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 — the CIF reader and the supercell builder.
 *
 * The case that matters is symmetry expansion. A CIF usually holds only the
 * asymmetric unit, and a reader that ignores the operations returns a perfectly
 * well-formed molecule with the wrong number of atoms — an error that survives
 * into whatever is computed from it. So the file here is a face-centred cell with
 * one site, whose four positions can be worked out on paper, and that is what the
 * result is compared against.
 */

#include "src/tools/cif.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

static int g_failed = 0;

static void check(bool ok, const std::string& what)
{
    std::printf("  %s  %s\n", ok ? "PASS" : "FAIL", what.c_str());
    if (!ok)
        ++g_failed;
}

static bool near(double a, double b, double tol = 1e-4)
{
    return std::fabs(a - b) < tol;
}

static std::string write(const std::string& name, const std::string& text)
{
    std::ofstream file(name);
    file << text;
    file.close();
    return name;
}

int main()
{
    std::printf("cif reader\n");

    // A cubic cell of 4 A with one site at the origin and the four face-centring
    // operations. On paper that is atoms at (0,0,0), (2,2,0), (2,0,2), (0,2,2).
    const std::string faceCentred = write("test_fcc.cif", R"CIF(
data_test
_cell_length_a    4.0
_cell_length_b    4.0
_cell_length_c    4.0
_cell_angle_alpha 90.0
_cell_angle_beta  90.0
_cell_angle_gamma 90.0
loop_
_symmetry_equiv_pos_as_xyz
  'x, y, z'
  'x+1/2, y+1/2, z'
  'x+1/2, y, z+1/2'
  'x, y+1/2, z+1/2'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
  Na1  Na  0.0  0.0  0.0
)CIF");

    {
        const curcuma::CifResult r = curcuma::ReadCif(faceCentred);
        check(r.ok(), "the file parses" + (r.ok() ? std::string() : " [" + r.error + "]"));
        check(near(r.cell.a, 4.0) && near(r.cell.gamma, 90.0), "the cell is read");
        check(r.asymmetric_atoms == 1, "one site is written in the file");
        check(r.symmetry_operations == 4, "and four symmetry operations are found");
        check(r.molecule.AtomCount() == 4,
            "so the cell holds four atoms -- ignoring the operations would have given one");

        // The four face-centring positions, in Angstrom.
        bool origin = false, ab = false, ac = false, bc = false;
        for (int i = 0; i < r.molecule.AtomCount(); ++i) {
            const Position p = r.molecule.Atom(i).second;
            if (near(p(0), 0) && near(p(1), 0) && near(p(2), 0)) origin = true;
            if (near(p(0), 2) && near(p(1), 2) && near(p(2), 0)) ab = true;
            if (near(p(0), 2) && near(p(1), 0) && near(p(2), 2)) ac = true;
            if (near(p(0), 0) && near(p(1), 2) && near(p(2), 2)) bc = true;
        }
        check(origin && ab && ac && bc,
            "and they sit exactly where the face-centring puts them");
    }

    // A site at a special position is generated more than once by the operations;
    // the duplicates have to go, or the cell comes out with too many atoms.
    {
        const std::string special = write("test_special.cif", R"CIF(
data_test
_cell_length_a    4.0
_cell_length_b    4.0
_cell_length_c    4.0
_cell_angle_alpha 90.0
_cell_angle_beta  90.0
_cell_angle_gamma 90.0
loop_
_symmetry_equiv_pos_as_xyz
  'x, y, z'
  '-x, -y, -z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
  C1  C  0.0  0.0  0.0
)CIF");
        const curcuma::CifResult r = curcuma::ReadCif(special);
        check(r.molecule.AtomCount() == 1,
            "an atom on an inversion centre is generated twice and counted once");
    }

    // No operations in the file: taken as P1, and said so rather than assumed.
    {
        const std::string p1 = write("test_p1.cif", R"CIF(
data_test
_cell_length_a    3.0
_cell_length_b    3.0
_cell_length_c    3.0
_cell_angle_alpha 90.0
_cell_angle_beta  90.0
_cell_angle_gamma 90.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
  O1  O  0.25  0.25  0.25
)CIF");
        const curcuma::CifResult r = curcuma::ReadCif(p1);
        check(r.molecule.AtomCount() == 1 && r.symmetry_operations == 1,
            "a file without symmetry operations is read as P1");
        bool said = false;
        for (const std::string& note : r.notes) {
            if (note.find("P1") != std::string::npos)
                said = true;
        }
        check(said, "and the note says so, rather than leaving it to be assumed");
        const Position p = r.molecule.Atom(0).second;
        check(near(p(0), 0.75) && near(p(1), 0.75) && near(p(2), 0.75),
            "fractional coordinates are turned into Angstrom with the cell");
    }

    // Replication.
    {
        const curcuma::CifResult r = curcuma::ReadCif(faceCentred);
        std::string error;
        const curcuma::Molecule big = curcuma::Supercell(r.molecule, 2, 1, 1, &error);
        check(error.empty() && big.AtomCount() == 8,
            "a 2x1x1 supercell of four atoms holds eight");
        check(near(big.getUnitCell()(0, 0), 8.0),
            "and carries the enlarged cell, so it can be replicated again");

        const curcuma::Molecule cube = curcuma::Supercell(r.molecule, 2, 2, 2, &error);
        check(cube.AtomCount() == 32, "2x2x2 holds thirty-two");

        curcuma::Molecule cellless;
        cellless.addPair({ 6, Position(0, 0, 0) });
        curcuma::Supercell(cellless, 2, 2, 2, &error);
        check(!error.empty() && error.find("no unit cell") != std::string::npos,
            "a structure without a cell is refused by name, not silently returned unchanged");
    }

    std::printf("%s (%d failed)\n", g_failed ? "FAIL" : "PASS", g_failed);
    return g_failed ? 1 : 0;
}
