/*
 * <CPCM-X: constants and internal helpers shared by cpcmx_*.cpp>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0 (see cpcmx.h).
 *
 * Claude Generated (Nov 2026). AI-generated, machine-tested only; human production testing pending.
 */

#pragma once

#include <string>
#include <vector>

namespace curcuma::solvation::cpcmx {

// CPCM-X globals.f90. Several constants are default-real (single precision) literals in the Fortran source and are
// reproduced with that precision: BtoA = 0.529177 and jtokcal = 0.000239006 are floats widened to double.
inline const double btoa = static_cast<double>(0.529177f);
inline const double jtokcal = static_cast<double>(0.000239006f);
constexpr double autokcal = 627.509;   ///< kcal/mol per Eh (627.509_wp)
constexpr double R = 8.314;            ///< J/(mol K) (8.314_wp)
constexpr double pi = 3.14159265358979323846;

// mctc-lib io/convert.f90 (used by numsa): autoaa = hbar / (m_e c alpha) * 1e10 from CODATA 2018, evaluated as the
// Fortran parameter expression; the value printed by mctc-lib with 17 digits is 0.529177210904492434.
constexpr double mctc_autoaa = 0.529177210904492434;
constexpr double mctc_aatoau = 1.0 / mctc_autoaa;

// cpcmx_data.cpp
double solventDensity(const std::string& name);   ///< kg/m^3, -1 if unknown (data.f90 density)
double atomicMass(int z);                         ///< data.f90 atomicMassNist
double covalentRadius(int z);                     ///< initialize.f90 setup_cov (Angstrom, 0 if not listed)
bool smdSolventProperties(const std::string& name, double& n, double& alpha, double& beta, double& msurft,
                          double& arom, double& fclbr);
const std::vector<std::string>& smdTextWater();
const std::vector<std::string>& smdTextOther();

} // namespace curcuma::solvation::cpcmx
