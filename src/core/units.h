/*
 * <Centralized Unit System and Physical Constants for Curcuma>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated: Educational unit system designed under Conrad's instruction
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * EDUCATIONAL PURPOSE: All physical constants and unit conversions centralized
 * for transparency and consistency across the Curcuma codebase
 */

#pragma once

namespace CurcumaUnit {

// ===== FUNDAMENTAL PHYSICAL CONSTANTS =====
// Reference: CODATA-2018 internationally recommended values
// Units: Atomic units (Hartree, Bohr, electron mass, elementary charge)

namespace Constants {
    // Mathematical constants
    constexpr double PI = 3.141592653589793238462643383279502884;
    constexpr double E_EULER = 2.718281828459045235360287471352662498;

    // Speed of light (atomic units)
    constexpr double SPEED_OF_LIGHT_AU = 137.035999084;

    // Boltzmann constant in different units
    constexpr double BOLTZMANN_HARTREE_K = 3.1668115634556e-6; // Eh/K
    constexpr double BOLTZMANN_SI = 1.380649e-23; // J/K
    constexpr double BOLTZMANN_EV_K = 8.617333262e-5; // eV/K

    // Avogadro's number
    constexpr double AVOGADRO = 6.02214076e23; // mol^-1

    // Gas constant
    constexpr double GAS_CONSTANT = 8.31446261815324; // J/(mol·K)

    // Atomic mass unit
    constexpr double ATOMIC_MASS_UNIT = 1.66053906660e-27; // kg
    constexpr double AMU_TO_AU = 1822.8884850; // m_e units

    // Time conversion factor
    constexpr double FS_TO_AMU = 41.34137314; // fs·amu conversion
    constexpr double ATOMIC_TIME_TO_FS = 24.188843265857; // aut to fs
}

// ===== UNIT CONVERSION FUNCTIONS =====
// Educational focus: Clear function names and comprehensive coverage

namespace Energy {
    // Base conversion factors (CODATA-2018)
    constexpr double HARTREE_TO_KJMOL = 2625.4996394798;
    constexpr double HARTREE_TO_KCALMOL = 627.5094740631;
    constexpr double HARTREE_TO_EV = 27.211386245988;
    constexpr double HARTREE_TO_WAVENUMBER = 219474.6313632; // cm^-1
    constexpr double EV_TO_HARTREE = 1.0 / HARTREE_TO_EV;

    // Forward conversions (from Hartree)
    inline constexpr double hartree_to_kjmol(double eh) { return eh * HARTREE_TO_KJMOL; }
    inline constexpr double hartree_to_kcalmol(double eh) { return eh * HARTREE_TO_KCALMOL; }
    inline constexpr double hartree_to_ev(double eh) { return eh * HARTREE_TO_EV; }
    inline constexpr double hartree_to_wavenumber(double eh) { return eh * HARTREE_TO_WAVENUMBER; }

    // Reverse conversions (to Hartree)
    inline constexpr double kjmol_to_hartree(double kjmol) { return kjmol / HARTREE_TO_KJMOL; }
    inline constexpr double kcalmol_to_hartree(double kcalmol) { return kcalmol / HARTREE_TO_KCALMOL; }
    inline constexpr double ev_to_hartree(double ev) { return ev * EV_TO_HARTREE; }
    inline constexpr double wavenumber_to_hartree(double cm1) { return cm1 / HARTREE_TO_WAVENUMBER; }

    // Inter-unit conversions (bypassing Hartree for efficiency)
    inline constexpr double kjmol_to_kcalmol(double kjmol) { return kjmol / 4.184; }
    inline constexpr double kcalmol_to_kjmol(double kcalmol) { return kcalmol * 4.184; }
}

namespace Length {
    // Base conversion factor (CODATA-2018)
    constexpr double BOHR_TO_ANGSTROM = 0.529177210903;
    constexpr double ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM;
    constexpr double BOHR_TO_METER = 5.29177210903e-11;
    constexpr double ANGSTROM_TO_METER = 1.0e-10;

    // Length conversions
    inline constexpr double bohr_to_angstrom(double bohr) { return bohr * BOHR_TO_ANGSTROM; }
    inline constexpr double angstrom_to_bohr(double ang) { return ang * ANGSTROM_TO_BOHR; }
    inline constexpr double bohr_to_meter(double bohr) { return bohr * BOHR_TO_METER; }
    inline constexpr double angstrom_to_meter(double ang) { return ang * ANGSTROM_TO_METER; }
    inline constexpr double bohr_to_nm(double bohr) { return bohr * BOHR_TO_ANGSTROM * 0.1; }
    inline constexpr double angstrom_to_nm(double ang) { return ang * 0.1; }
}

namespace Time {
    // Base conversion factors
    constexpr double ATOMIC_TIME_TO_FS = 24.188843265857; // aut to fs
    constexpr double FS_TO_ATOMIC_TIME = 1.0 / ATOMIC_TIME_TO_FS;

    // Time conversions
    inline constexpr double atomic_time_to_fs(double aut) { return aut * ATOMIC_TIME_TO_FS; }
    inline constexpr double fs_to_atomic_time(double fs) { return fs * FS_TO_ATOMIC_TIME; }
    inline constexpr double atomic_time_to_ps(double aut) { return aut * ATOMIC_TIME_TO_FS * 1e-3; }
    inline constexpr double ps_to_atomic_time(double ps) { return ps * 1e3 / ATOMIC_TIME_TO_FS; }
}

namespace ElectricDipole {
    // Dipole moment conversions
    // Reference: 1 Debye = 3.33564e-30 C·m = 0.208194 e·Å
    constexpr double DEBYE_TO_E_ANGSTROM = 0.208194;
    constexpr double E_ANGSTROM_TO_DEBYE = 1.0 / DEBYE_TO_E_ANGSTROM;

    inline constexpr double debye_to_e_angstrom(double debye) { return debye * DEBYE_TO_E_ANGSTROM; }
    inline constexpr double e_angstrom_to_debye(double e_ang) { return e_ang * E_ANGSTROM_TO_DEBYE; }
}

namespace Temperature {
    // Temperature conversions
    inline constexpr double celsius_to_kelvin(double celsius) { return celsius + 273.15; }
    inline constexpr double kelvin_to_celsius(double kelvin) { return kelvin - 273.15; }
    inline constexpr double fahrenheit_to_kelvin(double fahrenheit)
    {
        return (fahrenheit - 32.0) * 5.0 / 9.0 + 273.15;
    }
    inline constexpr double kelvin_to_fahrenheit(double kelvin)
    {
        return (kelvin - 273.15) * 9.0 / 5.0 + 32.0;
    }
}

// ===== FORCE FIELD UNIT CONVERSIONS =====
// Special conversions needed for different force field implementations

namespace ForceField {
    // Common force field energy unit conversions
    constexpr double HARTREE_TO_KCALMOL_FACTOR = 1.0 / 627.5094740631; // For FF parameters
    constexpr double KCALMOL_TO_HARTREE_FACTOR = 627.5094740631;

    // UFF/QMDFF specific conversion (historical factor used in codebase)
    constexpr double LEGACY_CONVERSION_FACTOR = 1.0 / 2625.15 * 4.19; // Used in legacy FF code

    inline constexpr double kcalmol_force_to_hartree(double kcalmol)
    {
        return kcalmol * HARTREE_TO_KCALMOL_FACTOR;
    }
    inline constexpr double hartree_force_to_kcalmol(double hartree)
    {
        return hartree * KCALMOL_TO_HARTREE_FACTOR;
    }
}

// ===== UTILITY FUNCTIONS =====
// Smart unit formatting for educational output

namespace Format {
    // Automatic unit selection based on value ranges
    inline std::string format_energy_relative(double hartree_diff)
    {
        double kjmol = Energy::hartree_to_kjmol(std::abs(hartree_diff));
        if (kjmol < 1.0) {
            return std::to_string(Energy::hartree_to_ev(hartree_diff) * 1000.0) + " meV";
        } else if (kjmol < 100.0) {
            return std::to_string(kjmol) + " kJ/mol";
        } else {
            return std::to_string(Energy::hartree_to_kcalmol(hartree_diff)) + " kcal/mol";
        }
    }

    inline std::string format_energy_absolute(double hartree)
    {
        return std::to_string(hartree) + " Eh";
    }

    inline std::string format_length(double bohr)
    {
        double angstrom = Length::bohr_to_angstrom(bohr);
        if (angstrom < 100.0) {
            return std::to_string(angstrom) + " Å";
        } else {
            return std::to_string(Length::bohr_to_nm(bohr)) + " nm";
        }
    }

    inline std::string format_dipole(double e_angstrom)
    {
        return std::to_string(ElectricDipole::e_angstrom_to_debye(e_angstrom)) + " D";
    }
}

// ===== COMPATIBILITY ALIASES =====
// For backward compatibility with existing code

namespace Legacy {
    // Aliases for existing scattered constants
    constexpr double au = Length::BOHR_TO_ANGSTROM; // From global.h
    constexpr double eV2Eh = Energy::EV_TO_HARTREE; // From global.h
    constexpr double fs2amu = Constants::FS_TO_AMU; // From global.h
    constexpr double au2eV = Energy::HARTREE_TO_EV; // From gfn2-xtb_param.hpp
}

} // namespace CurcumaUnit