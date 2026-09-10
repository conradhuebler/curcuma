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

struct CifResult {
    Molecule molecule;
    CifCell cell;
    int asymmetric_atoms = 0;   ///< atoms as written in the file
    int symmetry_operations = 0;///< operations found and applied (1 = P1, i.e. none)
    bool fractional = true;     ///< false when the file gave Cartesian coordinates
    std::string error;          ///< empty on success
    std::vector<std::string> notes;  ///< what was ignored, in plain words

    bool ok() const { return error.empty(); }
};

/**
 * @brief Read @p filename as CIF.
 *
 * Symmetry operations are applied to the asymmetric unit and the results are
 * deduplicated by position (0.1 Angstrom), which is what turns the contents of the
 * file into the contents of the cell. @a symmetry_operations says how many were
 * used, so "did this expand" has an answer rather than an assumption.
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
