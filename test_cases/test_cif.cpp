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

#include <clocale>
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

    // Disorder: two alternatives of one site, 0.03 A apart, in groups 1 and 2.
    // The old reader compared every generated atom with every earlier one of the
    // same element and swallowed the second alternative; a structure refined in two
    // conformations lost atoms that way. Claude Generated (Sep 2026).
    {
        const std::string file = write("test_disorder.cif", R"CIF(
data_test
_cell_length_a    10.0
_cell_length_b    10.0
_cell_length_c    10.0
_cell_angle_alpha 90.0
_cell_angle_beta  90.0
_cell_angle_gamma 90.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_disorder_assembly
_atom_site_disorder_group
  C1   C  0.100  0.100  0.100  1          .  .
  C2A  C  0.300  0.300  0.300  0.634(18)  A  1
  C2B  C  0.303  0.300  0.300  0.366(18)  A  2
)CIF");
        const curcuma::CifData data = curcuma::ReadCifData(file);
        check(data.ok() && data.sites.size() == 3, "three sites are read, occupancy columns and all");
        check(data.disorder_groups.size() == 2 && data.majorDisorderGroup() == 1,
            "two disorder groups, and the one with occupancy 0.634 is the major");
        check(curcuma::ReadCif(file).molecule.AtomCount() == 3,
            "with every alternative, both stay -- 0.03 A apart is not one atom");

        curcuma::CifBuildOptions major;
        major.select_disorder_group = true;
        major.disorder_group = 1;
        const curcuma::Molecule one = curcuma::BuildCif(data, major);
        curcuma::CifBuildOptions minor = major;
        minor.disorder_group = 2;
        const curcuma::Molecule two = curcuma::BuildCif(data, minor);
        check(one.AtomCount() == 2 && two.AtomCount() == 2,
            "choosing a group keeps the ordered site and that group's alternative");
        check(near(one.Atom(1).second(0), 3.0) && near(two.Atom(1).second(0), 3.03),
            "and each conformation brings its own position");
    }

    // A site 0.00001 from a cell face under an inversion: -x gives 0.99999, the same
    // position one cell over. The old check compared Cartesian distances and kept
    // both, one on each face.
    {
        const std::string file = write("test_face.cif", R"CIF(
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
  O1  O  0.00001  0.0  0.5
)CIF");
        check(curcuma::ReadCif(file).molecule.AtomCount() == 1,
            "a site on an inversion centre at a cell face is one atom, not one per face");
    }

    // The asymmetric unit is shown as written: a coordinate outside the cell stays
    // there, so a molecule written across a face is not cut apart. The unit cell
    // wraps it.
    {
        const std::string file = write("test_asym.cif", R"CIF(
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
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
  C1  C  1.25  0.0  0.0
)CIF");
        const curcuma::CifData data = curcuma::ReadCifData(file);
        curcuma::CifBuildOptions asym;
        asym.content = curcuma::CifContent::AsymmetricUnit;
        const curcuma::Molecule unit = curcuma::BuildCif(data, asym);
        check(unit.AtomCount() == 1 && near(unit.Atom(0).second(0), 5.0),
            "the asymmetric unit keeps the file's coordinate, 1.25 a = 5 A");
        const curcuma::Molecule cell = curcuma::BuildCif(data, curcuma::CifBuildOptions());
        check(cell.AtomCount() == 2 && near(cell.Atom(0).second(0), 1.0),
            "the unit cell holds both images, wrapped: 0.25 a = 1 A");
    }

    // Two data blocks: one structure is read, not the union of both.
    {
        const std::string file = write("test_blocks.cif", R"CIF(
data_first
_cell_length_a    3.0
_cell_length_b    3.0
_cell_length_c    3.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
  Na1  Na  0.0  0.0  0.0
data_second
_cell_length_a    9.0
_cell_length_b    9.0
_cell_length_c    9.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
  Cl1  Cl  0.5  0.5  0.5
)CIF");
        const curcuma::CifResult r = curcuma::ReadCif(file);
        bool said = false;
        for (const std::string& note : r.notes) {
            if (note.find("more than one data block") != std::string::npos)
                said = true;
        }
        check(r.molecule.AtomCount() == 1 && near(r.cell.a, 3.0),
            "of two data blocks only the first is read, with its own cell");
        check(said, "and a note says the rest was left out");
    }

    // Displacement parameters. A 4-fold rotation (y, -x, z) in a tetragonal cell
    // must swap U11 and U22 -- the ellipsoid turns with its atom. B is converted
    // to U. Claude Generated (Sep 2026).
    {
        const std::string file = write("test_adp.cif", R"CIF(
data_test
_cell_length_a    4.0
_cell_length_b    4.0
_cell_length_c    6.0
_cell_angle_alpha 90.0
_cell_angle_beta  90.0
_cell_angle_gamma 90.0
loop_
_symmetry_equiv_pos_as_xyz
  'x, y, z'
  'y, -x, z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_B_iso_or_equiv
  C1  C  0.1  0.2  0.3  1.5791367
  H1  H  0.3  0.1  0.3  3.9478418
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_23
  C1  0.010  0.030  0.020  0  0  0
)CIF");
        const curcuma::CifData data = curcuma::ReadCifData(file);
        check(data.sites.size() == 2 && data.sites[0].has_aniso && !data.sites[1].has_aniso,
            "the aniso loop is matched to its site by label");
        check(near(data.sites[1].u_iso, 0.05), "B_iso 3.948 is U_iso 0.05 (B = 8 pi^2 U)");
        std::vector<Eigen::Matrix3d> u;
        const curcuma::Molecule cell = curcuma::BuildCif(data, curcuma::CifBuildOptions(), &u);
        check(u.size() == size_t(cell.AtomCount()) && u.size() == 4,
            "one displacement tensor per atom, in atom order");
        check(near(u[0](0, 0), 0.010) && near(u[0](1, 1), 0.030) && near(u[0](2, 2), 0.020),
            "the original image carries U as written (orthogonal cell: U_cart = U)");
        check(near(u[1](0, 0), 0.030) && near(u[1](1, 1), 0.010) && near(u[1](2, 2), 0.020),
            "the 4-fold image has U11 and U22 swapped");
        check(near(u[2](0, 0), 0.05) && near(u[2](1, 1), 0.05) && near(u[2](0, 1), 0.0),
            "an isotropic atom gets U_iso times the identity");
    }

    // Oblique axes: in a monoclinic cell trace(U_cart)/3 must be the textbook U_eq,
    // U_eq = 1/3 [U11 a^2 a*^2 + U22 b^2 b*^2 + U33 c^2 c*^2 + 2 U13 a c a* c* cos(beta)].
    {
        const std::string file = write("test_adp_mono.cif", R"CIF(
data_test
_cell_length_a    7.0
_cell_length_b    9.0
_cell_length_c    11.0
_cell_angle_alpha 90.0
_cell_angle_beta  100.0
_cell_angle_gamma 90.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
  S1  S  0.25  0.25  0.25
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_23
  S1  0.020  0.030  0.040  0  0.005  0
)CIF");
        std::vector<Eigen::Matrix3d> u;
        curcuma::BuildCif(curcuma::ReadCifData(file), curcuma::CifBuildOptions(), &u);
        const double beta = 100.0 * M_PI / 180.0;
        const double sb = std::sin(beta);
        // Monoclinic reciprocal lengths: a* = 1/(a sin b), b* = 1/b, c* = 1/(c sin b).
        const double a = 7.0, b = 9.0, c = 11.0;
        const double as = 1.0 / (a * sb), bs = 1.0 / b, cs = 1.0 / (c * sb);
        const double ueq = (0.020 * a * a * as * as + 0.030 * b * b * bs * bs
                               + 0.040 * c * c * cs * cs + 2.0 * 0.005 * a * c * as * cs * std::cos(beta))
            / 3.0;
        check(u.size() == 1 && near(u[0].trace() / 3.0, ueq, 1e-9),
            "monoclinic: trace(U_cart)/3 equals the textbook U_eq, cos(beta) term and all");
        check(u.size() == 1 && near((u[0] - u[0].transpose()).norm(), 0.0, 1e-12),
            "and U_cart is symmetric");
    }

    // Complete molecules: a C=O pair across the x face. Wrapped image by image the
    // two atoms land at opposite faces; completed, they are bonded again and the
    // pair's centroid sits in the cell. The header values come along.
    // Claude Generated (Sep 2026).
    {
        const std::string file = write("test_complete.cif", R"CIF(
data_test
_space_group_name_H-M_alt 'P 1'
_space_group_IT_number 1
_cell_formula_units_Z 1
_chemical_formula_sum 'C O'
_chemical_formula_moiety
;
C O
;
_cell_length_a    10.0
_cell_length_b    10.0
_cell_length_c    10.0
_cell_angle_alpha 90.0
_cell_angle_beta  90.0
_cell_angle_gamma 90.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
  C1  C  0.93  0.5  0.5
  O1  O  0.05  0.5  0.5
)CIF");
        const curcuma::CifData data = curcuma::ReadCifData(file);
        check(data.space_group == "P 1" && data.space_group_number == 1 && data.formula_units_z == 1
                && data.formula_sum == "C O" && data.formula_moiety == "C O",
            "space group, IT number, Z and the formula are read as the file states them");
        const curcuma::Molecule cut = curcuma::BuildCif(data, curcuma::CifBuildOptions());
        check(near((cut.Atom(0).second - cut.Atom(1).second).norm(), 8.8),
            "wrapped image by image, the two atoms sit 8.8 A apart at opposite faces");
        curcuma::CifBuildOptions complete;
        complete.complete_molecules = true;
        const curcuma::Molecule whole = curcuma::BuildCif(data, complete);
        const Position centre = 0.5 * (whole.Atom(0).second + whole.Atom(1).second);
        check(near((whole.Atom(0).second - whole.Atom(1).second).norm(), 1.2),
            "completed, they are bonded again at 1.2 A");
        check(centre(0) >= 0.0 && centre(0) < 10.0, "and the molecule's centre lies in the cell");
    }

    // The numbers are read the same under a locale whose decimal separator is a
    // comma. A GUI host (qurcuma, via QApplication) adopts the system locale, and
    // std::stod then read "0.25" as 0 and "3.5" as 3. Skipped, not failed, where
    // no such locale is installed. Claude Generated (Sep 2026).
    {
        const char* german = std::setlocale(LC_NUMERIC, "de_DE.UTF-8");
        if (!german)
            german = std::setlocale(LC_NUMERIC, "de_DE.utf8");
        if (!german) {
            std::printf("  SKIP  no de_DE locale installed; comma-decimal parsing not tested\n");
        } else {
            const std::string file = write("test_locale.cif", R"CIF(
data_test
_cell_length_a    3.5
_cell_length_b    3.5
_cell_length_c    3.5
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
            const curcuma::CifResult r = curcuma::ReadCif(file);
            std::setlocale(LC_NUMERIC, "C");
            const Position p = r.ok() && r.molecule.AtomCount() == 1
                ? r.molecule.Atom(0).second : Position(0, 0, 0);
            check(near(r.cell.a, 3.5) && near(p(0), 0.875) && near(p(2), 0.875),
                "under a comma-decimal locale the cell and the coordinates still read as written");
        }
    }

    std::printf("%s (%d failed)\n", g_failed ? "FAIL" : "PASS", g_failed);
    return g_failed ? 1 : 0;
}
