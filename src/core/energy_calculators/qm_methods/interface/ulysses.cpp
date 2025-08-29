/*
 * < C++ Ulysses Interface >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/global_config.h"

#ifdef USE_BLAS
#define EIGEN_USE_BLAS
#endif

#ifdef USE_MKL
#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACK
#endif

#include "Eigen/Dense"
#include <vector>

#include "QC.hpp"

#include "BSet.hpp"
#include "GFN.hpp"
#include "MNDOd.hpp"
#include "Molecule.hpp"

#include <algorithm>
#include <iostream>
#include <string>

#include "src/core/curcuma_logger.h"
#include "ulysses.h"

static matrixE Geom2Matrix(const Geometry& geometry)
{
    matrixE matrix(geometry.rows(), 3);
    for (int i = 0; i < geometry.rows(); ++i) {
        matrix.thematrix(i, 0) = geometry(i, 0);
        matrix.thematrix(i, 1) = geometry(i, 1);
        matrix.thematrix(i, 2) = geometry(i, 2);
    }
    return matrix;
}

static Geometry Matrix2Geom(const matrixE& matrix)
{
    Geometry geometry(matrix.rows() / 3, 3);
    for (int i = 0; i < matrix.rows(); ++i) {
        geometry.data()[i] = matrix.thematrix(i);
    }
    return geometry;
}

UlyssesObject::UlyssesObject()
{
}

UlyssesObject::~UlyssesObject()
{
    delete m_bset;
    delete m_electron;
}

void UlyssesObject::Calculate(bool gradient, bool verbose)
{
    if (m_solvent != "none")
        m_electron->setSolvent(m_solvent);
    if (m_method == "gfn2")
        m_electron->setElectronTemp(m_Tele);
    // Claude Generated: Universal calculation for all methods
    // All Ulysses methods use the same Calculate interface
    m_electron->Calculate(int(verbose), m_SCFmaxiter);
    m_energy = m_electron->getEnergy();

    if (gradient) {
        matrixE grad;
        m_electron->AnalyticalGrad(grad);
        m_gradient = Matrix2Geom(grad);
    }
}

void UlyssesObject::setMethod(std::string& method)
{
    // Claude Generated: Simplified - no string parsing, just set the base method
    if (method == "ugfn2")
        m_method = "gfn2";
    else
        m_method = method;

    // Correction is now set separately in setMolecule
}

void UlyssesObject::setMolecule(const Geometry& geom, const std::vector<int>& atm, int chrge, int multpl, std::string pg, std::string correction)
{
    // Claude Generated: Set correction from parameter instead of string parsing
    m_correction = correction;

    Molecule mol;
    auto matrix = Geom2Matrix(geom);
    std::vector<size_t> atom_numbers;
    for (int i = 0; i < atm.size(); i++) {
        atom_numbers.push_back(atm[i]);
    }
    mol.set2System(matrix, atom_numbers, chrge, multpl, pg);

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Creating Ulysses BSet and Electron objects");
        CurcumaLogger::param("method", m_method);
        CurcumaLogger::param("correction", m_correction);
        CurcumaLogger::param("atoms", (int)atm.size());
        CurcumaLogger::param("charge", chrge);
        CurcumaLogger::param("multiplicity", multpl);
    }

    m_bset = new BSet(mol, m_method);

    if (m_method == "gfn2") {
        m_electron = new GFN2(*m_bset, mol);
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Created GFN2 electron object (no correction support)");
        }
    } else if (m_method == "pm6") {
        m_electron = new PM6(*m_bset, mol, "0", m_correction);
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Created PM6 electron object");
            CurcumaLogger::param("correction_applied", m_correction);
        }
    } else if (m_method == "am1") {
        m_electron = new AM1(*m_bset, mol, "0", m_correction);
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Created AM1 electron object");
            CurcumaLogger::param("correction_applied", m_correction);
        }
    } else if (m_method == "pm3") {
        m_electron = new PM3(*m_bset, mol, "0", m_correction);
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Created PM3 electron object");
            CurcumaLogger::param("correction_applied", m_correction);
        }
    } else if (m_method == "mndo") {
        m_electron = new MNDO(*m_bset, mol, "0", m_correction);
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Created MNDO electron object");
            CurcumaLogger::param("correction_applied", m_correction);
        }
    } else if (m_method == "mndod") {
        m_electron = new MNDOd(*m_bset, mol, "0", m_correction);
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Created MNDOd electron object");
            CurcumaLogger::param("correction_applied", m_correction);
        }
    } else if (m_method == "rm1") {
        m_electron = new RM1(*m_bset, mol, "0", m_correction);
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Created RM1 electron object");
            CurcumaLogger::param("correction_applied", m_correction);
        }
    } else if (m_method == "pm3pddg") {
        m_electron = new PM3PDDG(*m_bset, mol, "0", m_correction);
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Created PM3PDDG electron object");
            CurcumaLogger::param("correction_applied", m_correction);
        }
    } else if (m_method == "mndopddg") {
        m_electron = new MNDOPDDG(*m_bset, mol, "0", m_correction);
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Created MNDOPDDG electron object");
            CurcumaLogger::param("correction_applied", m_correction);
        }
    } else if (m_method == "pm3bp") {
        m_electron = new PM3BP(*m_bset, mol, "0", m_correction);
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Created PM3BP electron object");
            CurcumaLogger::param("correction_applied", m_correction);
        }
    } else {
        // Fallback to PM6 for unknown methods
        CurcumaLogger::warn("Unknown Ulysses method '" + m_method + "', falling back to PM6");
        m_electron = new PM6(*m_bset, mol, "0", m_correction);
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Created PM6 electron object (fallback)");
            CurcumaLogger::param("correction_applied", m_correction);
        }
    }
}

void UlyssesObject::UpdateGeometry(const Geometry& geom)
{
    auto matrix = Geom2Matrix(geom);
    // All methods use the same setGeometry interface
    m_electron->setGeometry(matrix);
}

Vector UlyssesObject::Charges() const
{
    std::vector<double> AtmCharge = m_electron->getQAtoms();
    Vector charges(AtmCharge.size());

    for (int i = 0; i < AtmCharge.size(); i++) {
        charges(i) = AtmCharge[i];
    }

    return Vector();
}

Vector UlyssesObject::OrbitalEnergies() const
{
    std::vector<double> OrbE = m_electron->EnMOs();
    Vector OrbE_vector(OrbE.size());
    for (int i = 0; i < OrbE.size(); i++) {
        OrbE_vector(i) = OrbE[i];
    }
    return OrbE_vector;
}

Vector UlyssesObject::OrbitalOccupations() const
{
    std::vector<double> OrbOcc = m_electron->occ();
    Vector OrbOcc_vector(OrbOcc.size());
    for (int i = 0; i < OrbOcc.size(); i++) {
        OrbOcc_vector(i) = OrbOcc[i];
    }
    return OrbOcc_vector;
}