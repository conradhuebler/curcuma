/*
 * <Extended Hückel Theory Parameters Module>
 * Copyright (C) 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 */

#pragma once

#include "STOIntegrals.hpp"
#include <map>
#include <string>
#include <vector>

namespace EHTParams {

/**
 * @brief Atomic orbital parameters for Extended Hückel Theory
 *
 * Contains the essential parameters for each atomic orbital type:
 * - zeta: Slater orbital exponent
 * - VSIP: Valence State Ionization Potential (in eV)
 * - n_electrons: Number of electrons contributed by this orbital
 */
struct OrbitalParameters {
    double zeta; ///< Slater orbital exponent
    double VSIP; ///< Valence State Ionization Potential in eV
    int n_electrons; ///< Number of electrons contributed by this orbital

    OrbitalParameters(double z = 0.0, double vsip = 0.0, int n_elec = 0)
        : zeta(z)
        , VSIP(vsip)
        , n_electrons(n_elec)
    {
    }
};

/**
 * @brief Complete set of orbital parameters for a given element
 *
 * Maps orbital types (S, PX, PY, PZ, etc.) to their respective parameters
 */
using ElementParameters = std::map<STO::OrbitalType, OrbitalParameters>;

/**
 * @brief EHT parameter database for different elements
 *
 * This class provides a centralized repository of Extended Hückel Theory
 * parameters for various elements. Parameters are based on literature values
 * from Hoffmann and other standard sources.
 */
class ParameterDatabase {
public:
    /**
     * @brief Get orbital parameters for a specific element and orbital type
     * @param atomic_number Atomic number of the element
     * @param orbital_type Type of orbital (S, PX, PY, PZ, etc.)
     * @return OrbitalParameters containing zeta, VSIP, and electron count
     * @throws std::runtime_error if element or orbital type not found
     */
    static OrbitalParameters getOrbitalParameters(int atomic_number, STO::OrbitalType orbital_type);

    /**
     * @brief Get all orbital parameters for a specific element
     * @param atomic_number Atomic number of the element
     * @return ElementParameters map containing all orbital types for the element
     * @throws std::runtime_error if element not found
     */
    static ElementParameters getElementParameters(int atomic_number);

    /**
     * @brief Check if parameters exist for a given element
     * @param atomic_number Atomic number to check
     * @return true if parameters exist, false otherwise
     */
    static bool hasElement(int atomic_number);

    /**
     * @brief Get the total number of valence electrons for an element
     * @param atomic_number Atomic number of the element
     * @return Total number of valence electrons
     */
    static int getValenceElectrons(int atomic_number);

    /**
     * @brief Get the Wolfsberg-Helmholz constant K
     * @return K value (typically 1.75)
     */
    static double getWolfsbergHelmholzConstant() { return 1.75; }

    /**
     * @brief Get element name from atomic number
     * @param atomic_number Atomic number
     * @return Element symbol as string
     */
    static std::string getElementSymbol(int atomic_number);

private:
    /**
     * @brief Initialize the parameter database
     *
     * Sets up the internal parameter maps with literature values for:
     * - H (Z=1): 1s orbital
     * - C (Z=6): 2s, 2px, 2py, 2pz orbitals
     * - N (Z=7): 2s, 2px, 2py, 2pz orbitals
     * - O (Z=8): 2s, 2px, 2py, 2pz orbitals
     * - F (Z=9): 2s, 2px, 2py, 2pz orbitals
     */
    static void initializeDatabase();

    /**
     * @brief Internal parameter storage
     * Maps atomic number -> orbital type -> parameters
     */
    static std::map<int, ElementParameters> s_parameters;
    static bool s_initialized;
};

/**
 * @brief STO-6G basis set coefficients for different elements and orbitals
 *
 * Contains the Gaussian primitive coefficients and exponents needed to
 * represent Slater-Type Orbitals using 6 Gaussian functions (STO-6G basis).
 */
struct STO6GParameters {
    std::vector<double> alpha; ///< Gaussian exponents
    std::vector<double> coeff; ///< Contraction coefficients

    STO6GParameters() = default;
    STO6GParameters(const std::vector<double>& a, const std::vector<double>& c)
        : alpha(a)
        , coeff(c)
    {
    }
};

/**
 * @brief STO-6G basis set database
 *
 * Provides STO-6G basis set parameters for different elements and orbital types.
 * These are used for more accurate integral evaluation when needed.
 */
class STO6GDatabase {
public:
    /**
     * @brief Get STO-6G parameters for a specific element and orbital type
     * @param atomic_number Atomic number of the element
     * @param orbital_type Type of orbital (S, PX, PY, PZ)
     * @return STO6GParameters containing alpha and coefficient vectors
     */
    static STO6GParameters getParameters(int atomic_number, STO::OrbitalType orbital_type);

    /**
     * @brief Check if STO-6G parameters exist for given element and orbital
     * @param atomic_number Atomic number
     * @param orbital_type Orbital type
     * @return true if parameters exist
     */
    static bool hasParameters(int atomic_number, STO::OrbitalType orbital_type);

private:
    static void initializeDatabase();
    static std::map<std::pair<int, STO::OrbitalType>, STO6GParameters> s_sto6g_parameters;
    static bool s_initialized;
};

} // namespace EHTParams