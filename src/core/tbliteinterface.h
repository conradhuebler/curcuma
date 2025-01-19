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

#include "src/core/global.h"

#ifndef tblite_delete
#include "tblite.h"
#endif

#include "src/core/interface/abstract_interface.h"

static json TBLiteSettings{
    { "tb_acc", 1 },
    { "tb_max_iter", 250 },
    { "tb_damping", 0.4 },
    { "tb_temp", 9.500e-4 },
    { "tb_verbose", 0 },
    { "tb_guess", "SAD" },
    { "solv", "none" },
    { "solv_eps", -1 }
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
    double Calculation(double* gradient = 0, bool verbose = false);

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
    int m_maxiter = 100;
    int m_verbose = 0;
    int m_guess = 0;
    int m_error_count = 0;
    double m_damping = 0.5;
    double m_temp = 1000;
    double m_solv_eps = -1;
    char* m_alpb_solv = NULL;
    int m_gb_type = 0, m_born_kernel = 1;
    int m_solv_param = 12;
    int m_solv_ref = 1;
    bool m_cpcm = false, m_alpb = false;
    tblite_error m_error = NULL;
    tblite_structure m_tblite_mol = NULL;
    tblite_result m_tblite_res = NULL;
    tblite_context m_ctx = NULL;
    tblite_calculator m_tblite_calc = NULL;
    tblite_container m_tb_cont = NULL;

    bool m_initialised = false, m_calculator = false;
    json m_tblitesettings;
};
