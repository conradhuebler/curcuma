/*
 * < C++ XTB and tblite Interface >
 * Copyright (C) 2020 - 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "interface/abstract_interface.h"

#ifndef tblite_delete
#include "tblite.h"
#endif

static json TBLiteSettings{
    { "tb_acc", 1 },
    { "SCFmaxiter", 100 },
    { "tb_damping", 0.4 },
    { "Tele", 300 },
    { "verbose", 0 },
    { "tb_guess", "SAD" },
    { "solvent_model", 0 }, // 0 - none; 1 - CPCM; 2 - GB, 3 - ALPB
    { "solvent_eps", -1 },
    { "solvent", "none" },
    { "solvent_gb_version", 0 }, // 0 - GBSA, 1 - ALPB
    { "solvent_gb_kernel", 1 },
    { "solvent_alpb_version", 12 },
    { "solvent_alpb_reference", 1 },
    { "spin", 0 }
};

class TBLiteInterface : public QMInterface {
public:
    TBLiteInterface(const json& tblitesettings = TBLiteSettings);
    ~TBLiteInterface();

    bool InitialiseMolecule(const Mol& mol) override;
    bool InitialiseMolecule(const Mol* mol) override;
    bool InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin) override;
    bool InitialiseMolecule();

    bool UpdateMolecule(const Mol& mol) override;
    bool UpdateMolecule(const Mol* mol) override;
    bool UpdateMolecule(const Matrix& geometry) override;
    bool UpdateMolecule(const double* coord) override;
    bool UpdateMolecule();

    bool Error() override { return m_error_count >= 10; }
    double Calculation(bool gradient = 0, bool verbose = false);

    void clear() override;

    Vector Charges() const override;
    Vector Dipole() const override;

    Vector BondOrders() const;
    Vector OrbitalEnergies() const override;
    Vector OrbitalOccupations() const override;
    void setMethod(const std::string& method) override;

private:
    void ApplySolvation();

    void tbliteError();
    void tbliteContextError();
    double* m_coord;
    int* m_attyp;

    double m_thr = 1.0e-10;
    int m_acc = 2;
    int m_SCFmaxiter = 100;
    int m_verbose = 0;
    int m_guess = 0;
    int m_error_count = 0;
    double m_damping = 0.5;
    double m_Tele = 298;
    double m_solvent_eps = -1;
    int m_solvent_model = 0;
    char* m_solvent = NULL;
    int m_solvent_gb_version = 0;
    int m_solvent_gb_kernel = 1;
    int m_solvent_alpb_version = 12;
    int m_solvent_alpb_reference = 1;

    tblite_error m_error = NULL;
    tblite_structure m_tblite_mol = NULL;
    tblite_result m_tblite_res = NULL;
    tblite_context m_ctx = NULL;
    tblite_calculator m_tblite_calc = NULL;
    tblite_container m_tb_cont = NULL;

    bool m_initialised = false, m_calculator = false;
    json m_tblitesettings;
};
