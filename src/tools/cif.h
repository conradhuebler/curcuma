/*
 * CIF reader and supercell builder.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 - Coordinates, the unit cell, and replication.
 *
 * The one thing a CIF reader must not do quietly: most CIF files contain only the
 * ASYMMETRIC UNIT, and the rest of the cell follows from the symmetry operations
 * listed alongside it. A reader that ignores them returns a fraction of the
 * structure -- a perfectly well-formed molecule with the wrong number of atoms,
 * which is the kind of error that survives all the way into a published number.
 * So the operations are applied, and what was applied is reported.
 *
 * Not a general CIF parser. It reads what is needed to get atoms into a cell:
 * cell lengths and angles, the atom_site loop (fractional or Cartesian), and the
 * symmetry operation loop. Anything else in the file is ignored, and the notes
 * say so rather than leaving it to be assumed.
 */

#pragma once

#include "src/core/molecule.h"

#include <string>
#include <vector>

namespace curcuma {

/// The unit cell as the file states it, plus the lattice it implies.
struct CifCell {
    double a = 0.0, b = 0.0, c = 0.0;              ///< Angstrom
    double alpha = 90.0, beta = 90.0, gamma = 90.0; ///< degrees
    /// Columns are the a, b and c vectors in Angstrom. Identity when invalid.
    Eigen::Matrix3d lattice = Eigen::Matrix3d::Identity();
    bool valid = false;
};

/// One atom site as the file writes it. Claude Generated (Sep 2026).
struct CifSite {
    std::string label;
    std::string element;
    Eigen::Vector3d coordinate = Eigen::Vector3d::Zero(); ///< fractional, or Angstrom when Cartesian
    double occupancy = 1.0;          ///< _atom_site_occupancy
    std::string disorder_assembly;   ///< _atom_site_disorder_assembly, empty for '.'
    int disorder_group = 0;          ///< _atom_site_disorder_group (SHELX PART); 0 = ordered
    /// Displacement parameters (Angstrom^2, B converted to U). Anisotropic: U11
    /// U22 U33 U12 U13 U23 as the aniso loop writes them, i.e. on the reciprocal-
    /// scaled crystal axes. Isotropic: U_iso or U_equiv from the site loop.
    bool has_aniso = false;
    double u_aniso[6] {};
    bool has_iso = false;
    double u_iso = 0.0;
};

/// A disorder group over the whole file: every site of SHELX PART n (sign
/// dropped -- PART -1 only changes how SHELX generates bonds at special
/// positions). In a structure refined in two conformations these are the two.
struct CifDisorderGroup {
    int group = 0;
    int sites = 0;                   ///< sites carrying this group
    double mean_occupancy = 0.0;     ///< mean _atom_site_occupancy of those sites
};

/// A symmetry operation: fractional x' = matrix * x + translation.
struct CifSymmetryOperation {
    double matrix[3][3] {};
    double translation[3] {};
};

/// Everything read from the file, before any operation is applied.
struct CifData {
    CifCell cell;
    std::vector<CifSite> sites;
    std::vector<CifSymmetryOperation> operations;  ///< identity alone when the file lists none
    bool fractional = true;          ///< false when the file gave Cartesian coordinates
    std::vector<CifDisorderGroup> disorder_groups; ///< sorted by group; empty when ordered
    /// As the file states them, empty when absent: _space_group_name_H-M_alt (or
    /// _symmetry_space_group_name_H-M), _space_group_IT_number, _cell_formula_units_Z,
    /// _chemical_formula_sum.
    std::string space_group;
    int space_group_number = 0;
    int formula_units_z = 0;
    std::string formula_sum;
    /// _chemical_formula_moiety, e.g. "C16 H41 N7 Si2 4+, 4(Cl -)" -- the one place a
    /// CIF says which particles carry which charge.
    std::string formula_moiety;
    std::string error;               ///< empty on success
    std::vector<std::string> notes;  ///< what was ignored, in plain words

    bool ok() const { return error.empty(); }
    /// The group with the highest mean occupancy (the major conformation); 0 when
    /// the structure is ordered.
    int majorDisorderGroup() const;
};

enum class CifContent {
    AsymmetricUnit, ///< the sites as written, positions not wrapped -- what the file shows
    UnitCell        ///< every symmetry image, wrapped into the cell -- what the crystal holds
};

struct CifBuildOptions {
    CifContent content = CifContent::UnitCell;
    /// Which sites: all of them (every disorder alternative at once) when false,
    /// or the ordered sites plus those of @a disorder_group when true.
    bool select_disorder_group = false;
    int disorder_group = 0;
    /// Unit cell only: put every molecule back together across the cell faces
    /// (bonded by covalent radii under periodic boundaries) and move it so its
    /// centroid lies in the cell -- Mercury's "complete molecules". Without it each
    /// image is wrapped on its own and molecules at a face come out in pieces.
    bool complete_molecules = false;
};

/// Read @p filename without building anything. Only the first data block that
/// holds atom sites is read; a note says so when there are more.
CifData ReadCifData(const std::string& filename);

/**
 * @brief Build a molecule from @p data.
 *
 * @p displacements, when given, receives one Cartesian displacement tensor per
 * atom, in atom order (Angstrom^2; the zero matrix where the file gives none).
 * Anisotropic U is taken to Cartesian axes as U_cart = M N U N M^T with M the
 * lattice and N = diag(a*, b*, c*), and a symmetry image turns it with the
 * operation's rotation, U* -> R U* R^T -- a thermal ellipsoid follows its atom.
 * Isotropic U is U_iso times the identity.
 *
 * For the unit cell every site is mapped through every operation and wrapped into
 * the cell; images of the SAME site closer than 0.1 Angstrom (periodically, so 0.0
 * and 0.99999 are one position) are one atom -- that is what a special position
 * is. Different sites are never merged: two disorder alternatives 0.05 Angstrom
 * apart are two atoms, and choosing between them is @a options' job.
 * The molecule carries the cell, so Supercell() can replicate a unit cell.
 */
Molecule BuildCif(const CifData& data, const CifBuildOptions& options,
    std::vector<Eigen::Matrix3d>* displacements = nullptr);

struct CifResult {
    Molecule molecule;
    CifCell cell;
    int asymmetric_atoms = 0;   ///< atoms as written in the file
    int symmetry_operations = 0;///< operations found and applied (1 = P1, i.e. none)
    bool fractional = true;     ///< false when the file gave Cartesian coordinates
    std::vector<CifDisorderGroup> disorder_groups; ///< empty when ordered
    std::string error;          ///< empty on success
    std::vector<std::string> notes;  ///< what was ignored, in plain words

    bool ok() const { return error.empty(); }
};

/**
 * @brief Read @p filename as CIF: the full unit cell with every site.
 *
 * ReadCifData() + BuildCif() with the unit cell and all disorder alternatives --
 * nothing is dropped silently; a disordered structure gets a note saying that its
 * alternatives overlap. @a symmetry_operations says how many were used, so "did
 * this expand" has an answer rather than an assumption.
 */
CifResult ReadCif(const std::string& filename);

/**
 * @brief Replicate @p molecule @p na x @p nb x @p nc along its cell vectors.
 *
 * Requires a unit cell on the molecule; without one the input is returned
 * unchanged and @p error, when given, says so. The result carries the enlarged
 * cell, so a supercell can be replicated again.
 */
Molecule Supercell(const Molecule& molecule, int na, int nb, int nc,
    std::string* error = nullptr);

} // namespace curcuma
