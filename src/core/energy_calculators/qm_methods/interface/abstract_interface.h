/*
 * < C++ Abstract Class for QM Interface >
 * Copyright (C) 2020 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/core/global.h"

class QMInterface {
public:
    QMInterface() = default;

    virtual ~QMInterface() = default;

    virtual bool InitialiseMolecule(const Mol& molecule)
    {
        m_atomcount = molecule.m_number_atoms;
        m_geometry = molecule.m_geometry;
        m_spin = molecule.m_spin;
        m_charge = molecule.m_charge;
        m_atoms = molecule.m_atoms;
        m_gradient = Matrix::Zero(m_atomcount, 3);

        return InitialiseMolecule();
    }
    virtual bool InitialiseMolecule(const Mol* molecule)
    {
        m_atomcount = molecule->m_number_atoms;
        m_geometry = molecule->m_geometry;
        m_spin = molecule->m_spin;
        m_charge = molecule->m_charge;
        m_atoms = molecule->m_atoms;
        m_gradient = Matrix::Zero(m_atomcount, 3);
        return InitialiseMolecule();
    }
    virtual bool InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin)
    {
        m_atomcount = natoms;
        m_charge = charge;
        m_spin = spin;
        m_atoms = std::vector<int>(attyp, attyp + natoms);
        m_geometry = Matrix::Zero(natoms, 3);
        for (int i = 0; i < natoms; ++i) {
            m_geometry(i, 0) = coord[i * 3];
            m_geometry(i, 1) = coord[i * 3 + 1];
            m_geometry(i, 2) = coord[i * 3 + 2];
        }
        m_gradient = Matrix::Zero(m_atomcount, 3);

        return InitialiseMolecule();
    }

    virtual bool InitialiseMolecule() { return true; }

    virtual bool UpdateMolecule(const Mol& molecule)
    {
        m_atomcount = molecule.m_number_atoms;
        m_geometry = molecule.m_geometry;
        m_spin = molecule.m_spin;
        m_charge = molecule.m_charge;
        m_atoms = molecule.m_atoms;

        return UpdateMolecule();
    }
    virtual bool UpdateMolecule(const Mol* molecule)
    {
        m_atomcount = molecule->m_number_atoms;
        m_geometry = molecule->m_geometry;
        m_spin = molecule->m_spin;
        m_charge = molecule->m_charge;
        m_atoms = molecule->m_atoms;

        return UpdateMolecule();
    }
    virtual bool UpdateMolecule(const double* coord)
    {
        for (int i = 0; i < m_atomcount; ++i) {
            m_geometry(i, 0) = coord[3 * i];
            m_geometry(i, 1) = coord[3 * i + 1];
            m_geometry(i, 2) = coord[3 * i + 2];
        }

        return UpdateMolecule();
    }
    virtual bool UpdateMolecule(const Matrix& geometry)
    {
        m_geometry = geometry;

        return UpdateMolecule();
    }
    virtual bool UpdateMolecule(const Vector& geometry)
    {
        for (int i = 0; i < m_atomcount; ++i) {
            m_geometry(i, 0) = geometry(3 * i);
            m_geometry(i, 1) = geometry(3 * i + 1);
            m_geometry(i, 2) = geometry(3 * i + 2);
        }

        return UpdateMolecule();
    }

    virtual bool UpdateMolecule()
    {
        // std::cout << m_geometry << std::endl;
        return true;
    }
    virtual bool Error() { return false; };
    virtual double Calculation(bool gradient = false) = 0;

    virtual void clear() {}

    virtual Vector Charges() const { return Vector{}; }
    virtual Vector Dipole() const { return Vector{}; }

    virtual Vector BondOrders() const { return Vector{}; }
    virtual Vector OrbitalEnergies() const { return Vector{}; }
    virtual Vector OrbitalOccupations() const { return Vector{}; }

    virtual void setMethod(const std::string& method) { m_method = method; }
    virtual Geometry Gradient() const { return m_gradient; }
    void setMult(int multi)
    {
        m_spin = multi - 1;
        m_muli = multi;
    }

    virtual bool hasGradient() const { return false; }

protected:
    bool m_initialised = false;
    Matrix m_geometry, m_gradient;
    double m_charge, m_spin = 0, m_muli = 1;
    std::vector<int> m_atoms;

    std::string m_method = "none";
    int m_method_switch = 0;
    int m_atomcount = 0;
};
