/*
 * <Implicit Solvation Model Interface>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (June 2026): Common interface for implicit solvation models
 * (ALPB, GBSA, CPCM) shared by the native GFN1/GFN2 SCF, the GPU SCF and GFN-FF.
 */

#pragma once

#include "src/core/global.h"

#include <string>
#include <vector>

/**
 * @brief Abstract interface for implicit (continuum) solvation models.
 *
 * The interface exposes the three hooks a self-consistent SCF needs:
 *
 *  1. @ref addPotential — the in-SCF Fock contribution v_i = dE_solv/dq_i.
 *     For all Born/GB/ALPB models the charge-dependent solvation energy is
 *     E = 1/2 * q^T B q (B = symmetric Born interaction matrix, already scaled
 *     by the dielectric factor keps), so the potential is exactly v = B q.
 *     Charge-independent parts (non-polar surface / CDS, free-energy shift)
 *     do NOT contribute to the potential — they only enter @ref energy and
 *     @ref addGradient.
 *
 *  2. @ref energy — the total solvation free energy for a given charge vector
 *     (in the 1/2 form, so it can be added to the SCF total energy WITHOUT
 *     double-counting the band-structure term, mirroring the Coulomb ES2 term).
 *
 *  3. @ref addGradient — the nuclear gradient of the solvation free energy,
 *     evaluated at the converged charges.
 *
 * All quantities are in atomic units (Bohr, Hartree). Atomic charges are the
 * per-atom partial charges (Mulliken for GFN1/GFN2, EEQ for GFN-FF).
 *
 * Claude Generated (June 2026).
 */
class ImplicitSolvationModel {
public:
    virtual ~ImplicitSolvationModel() = default;

    /**
     * @brief Initialise the model for a molecule, solvent and host method.
     * @param atomic_numbers Atomic numbers Z per atom.
     * @param solvent        Solvent name (e.g. "water", "dmso").
     * @param method         Host method selecting the parameter set
     *                       ("gfn1", "gfn2", "gfnff").
     * @return true on success; false if solvent/method has no parameters.
     */
    virtual bool init(const std::vector<int>& atomic_numbers,
                      const std::string& solvent,
                      const std::string& method) = 0;

    /**
     * @brief Recompute geometry-dependent state (Born radii, SASA, Born matrix).
     *        Must be called once per geometry before the SCF.
     * @param atomic_numbers Atomic numbers Z per atom.
     * @param xyz_bohr       Coordinates in Bohr (N x 3).
     */
    virtual void update(const std::vector<int>& atomic_numbers,
                        const Matrix& xyz_bohr) = 0;

    /**
     * @brief Add the in-SCF solvation potential v_i = dE_solv/dq_i to v_at.
     *        Called every SCF iteration with the current atomic charges.
     * @param q_at  Current atomic partial charges (length nat).
     * @param v_at  [in/out] Per-atom potential shift (length nat), accumulated.
     */
    virtual void addPotential(const Vector& q_at, Vector& v_at) const = 0;

    /**
     * @brief Total solvation free energy (1/2 form) for the given charges.
     * @param q_at Atomic partial charges.
     * @return Solvation free energy in Hartree.
     */
    virtual double energy(const Vector& q_at) const = 0;

    /**
     * @brief Add the solvation nuclear gradient to an existing gradient.
     * @param atomic_numbers Atomic numbers Z per atom.
     * @param xyz_bohr       Coordinates in Bohr (N x 3).
     * @param q_at           Converged atomic partial charges.
     * @param gradient       [in/out] Molecular gradient (Eh/Bohr), accumulated.
     */
    virtual void addGradient(const std::vector<int>& atomic_numbers,
                             const Matrix& xyz_bohr,
                             const Vector& q_at,
                             Matrix& gradient) = 0;

    /**
     * @brief Device (GPU) hook: the symmetric Born matrix B for v_at += B·q_at.
     *
     * Returns a pointer to the current nat×nat Born interaction matrix (the same B
     * that @ref addPotential contracts), valid after @ref update, so the GPU
     * potential build can add the reaction field on-device. Returns nullptr when the
     * device path is not applicable — in particular when the solute charge is not the
     * plain atomic charge (e.g. GFN1 CM5, where q_solute = q_at + cm5), so callers
     * must fall back to the host-driven loop. Default nullptr. Claude Generated (WP4b).
     */
    virtual const double* deviceBornMatrix() const { return nullptr; }
};
