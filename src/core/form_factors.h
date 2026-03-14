/*
 * <Atomic and CG form factors for scattering calculations - Claude Generated 2025>
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

#pragma once

#include <cmath>
#include <array>

/*! \brief Form factor database and calculation functions - Claude Generated 2025
 *
 * Provides atomic form factors (Cromer-Mann parameterization) and CG bead form factors
 * for scattering calculations (P(q), S(q)).
 *
 * Educational Focus:
 * - Transparent mathematical implementation (Cromer-Mann 9-parameter formula)
 * - Direct access to tabulated values (IUCr standard)
 * - Minimal abstraction - just data and calculation functions
 *
 * Data Source: International Tables for Crystallography, Vol. C, Table 6.1.1.4
 * Coverage: Elements H (Z=1) to Xe (Z=54) - sufficient for organic/polymer chemistry
 */

namespace FormFactors {

// CG element constant (defined in CG system)
constexpr int CG_ELEMENT = 226;

/*! \brief Cromer-Mann 9-parameter atomic form factor data
 *
 * Formula: f(q) = c + Σᵢ₌₁⁴ aᵢ exp(-bᵢ(q/4π)²)
 * where q is the scattering vector magnitude in Å⁻¹
 */
struct CromerMannParams {
    double a1, b1;  // First Gaussian amplitude and decay
    double a2, b2;  // Second Gaussian amplitude and decay
    double a3, b3;  // Third Gaussian amplitude and decay
    double a4, b4;  // Fourth Gaussian amplitude and decay
    double c;       // Constant term
};

/*! \brief Cromer-Mann parameters for elements Z=1 to Z=54 - Claude Generated 2025
 *
 * Data from International Tables for Crystallography, Vol. C (2006)
 * Index by atomic number: CROMER_MANN_DATA[Z] for element with atomic number Z
 *
 * NOTE: Index 0 is unused (no element with Z=0), array starts at index 1
 */
constexpr std::array<CromerMannParams, 55> CROMER_MANN_DATA = {{
    // Index 0 - placeholder (unused)
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},

    // Z=1: H (Hydrogen)
    {0.493002, 10.5109, 0.322912, 26.1257, 0.140191, 3.14236, 0.040810, 57.7997, 0.003038},

    // Z=2: He (Helium)
    {0.873400, 9.10370, 0.630900, 3.35680, 0.311200, 22.9276, 0.178000, 0.98210, 0.006400},

    // Z=3: Li (Lithium)
    {1.128200, 3.95460, 0.750800, 1.05240, 0.617500, 85.3905, 0.465300, 168.261, 0.037700},

    // Z=4: Be (Beryllium)
    {1.591900, 43.6427, 1.127800, 1.86230, 0.539100, 103.483, 0.702900, 0.54200, 0.038500},

    // Z=5: B (Boron)
    {2.054500, 23.2185, 1.332600, 1.02100, 1.097900, 60.3498, 0.706800, 0.14030, -0.19320},

    // Z=6: C (Carbon)
    {2.310000, 20.8439, 1.020000, 10.2075, 1.588600, 0.56870, 0.865000, 51.6512, 0.215600},

    // Z=7: N (Nitrogen)
    {12.2126, 0.00570, 3.132200, 9.89330, 2.012500, 28.9975, 1.166300, 0.58260, -11.529},

    // Z=8: O (Oxygen)
    {3.048500, 13.2771, 2.286800, 5.70110, 1.546300, 0.32390, 0.867000, 32.9089, 0.250800},

    // Z=9: F (Fluorine)
    {3.539200, 10.2825, 2.641200, 4.29440, 1.517000, 0.26150, 1.024300, 26.1476, 0.277600},

    // Z=10: Ne (Neon)
    {3.955300, 8.40420, 3.112500, 3.42620, 1.454600, 0.23060, 1.125100, 21.7184, 0.351500},

    // Z=11: Na (Sodium)
    {4.762600, 3.28500, 3.173600, 8.84220, 1.267400, 0.31360, 1.112800, 129.424, 0.676000},

    // Z=12: Mg (Magnesium)
    {5.420400, 2.82750, 2.173500, 79.2611, 1.226900, 0.38080, 2.307300, 7.19370, 0.858400},

    // Z=13: Al (Aluminum)
    {6.420200, 3.03870, 1.900200, 0.74260, 1.593600, 31.5472, 1.964600, 85.0886, 1.115100},

    // Z=14: Si (Silicon)
    {6.291500, 2.43860, 3.035300, 32.3337, 1.989100, 0.67850, 1.541000, 81.6937, 1.140700},

    // Z=15: P (Phosphorus)
    {6.434500, 1.90670, 4.179100, 27.1570, 1.780000, 0.52600, 1.490800, 68.1645, 1.114900},

    // Z=16: S (Sulfur)
    {6.905300, 1.46790, 5.203400, 22.2151, 1.437900, 0.25360, 1.586300, 56.172, 0.866900},

    // Z=17: Cl (Chlorine)
    {11.4604, 0.01040, 7.196400, 1.16620, 6.255600, 18.5194, 1.645500, 47.7784, -9.5574},

    // Z=18: Ar (Argon)
    {7.484500, 0.90720, 6.772300, 14.8407, 0.653900, 43.8983, 1.644200, 33.3929, 1.444500},

    // Z=19: K (Potassium)
    {8.218600, 12.7949, 7.439800, 0.77480, 1.051900, 213.187, 0.865900, 41.6841, 1.422800},

    // Z=20: Ca (Calcium)
    {8.626600, 10.4421, 7.387300, 0.65990, 1.589900, 85.7484, 1.021100, 178.437, 1.375100},

    // Z=21: Sc (Scandium)
    {9.189000, 9.02130, 7.367900, 0.57290, 1.640900, 136.108, 1.468000, 51.3531, 1.332900},

    // Z=22: Ti (Titanium)
    {9.759500, 7.85080, 7.355800, 0.50000, 1.699100, 35.6338, 1.902100, 116.105, 1.280700},

    // Z=23: V (Vanadium)
    {10.2971, 6.86570, 7.351100, 0.43850, 2.070300, 26.8938, 2.057100, 102.478, 1.219900},

    // Z=24: Cr (Chromium)
    {10.6406, 6.10380, 7.353700, 0.39200, 3.324000, 20.2626, 1.492200, 98.7399, 1.183200},

    // Z=25: Mn (Manganese)
    {11.2819, 5.34090, 7.357300, 0.34320, 3.019300, 17.8674, 2.244100, 83.7543, 1.089600},

    // Z=26: Fe (Iron)
    {11.7695, 4.76110, 7.357300, 0.30720, 3.522200, 15.3535, 2.304500, 76.8805, 1.036900},

    // Z=27: Co (Cobalt)
    {12.2841, 4.27910, 7.340900, 0.27840, 4.003400, 13.5359, 2.348800, 71.1692, 1.011800},

    // Z=28: Ni (Nickel)
    {12.8376, 3.87850, 7.292000, 0.25650, 4.443800, 12.1763, 2.380000, 66.3421, 1.034100},

    // Z=29: Cu (Copper)
    {13.3380, 3.58280, 7.167600, 0.24700, 5.615800, 11.3966, 1.673500, 64.8126, 1.191000},

    // Z=30: Zn (Zinc)
    {14.0743, 3.26550, 7.031800, 0.23330, 5.162500, 10.3163, 2.410000, 58.7097, 1.304100},

    // Z=31: Ga (Gallium)
    {15.2354, 3.06690, 6.700600, 0.24120, 4.359100, 10.7805, 2.962300, 61.4135, 1.718900},

    // Z=32: Ge (Germanium)
    {16.0816, 2.85090, 6.374700, 0.25160, 3.706900, 11.4468, 3.683000, 54.7625, 2.131300},

    // Z=33: As (Arsenic)
    {16.6723, 2.63450, 6.070100, 0.26470, 3.431300, 12.9479, 4.277900, 47.7972, 2.531000},

    // Z=34: Se (Selenium)
    {17.0006, 2.40980, 5.819600, 0.27260, 3.973100, 15.2372, 4.354300, 43.8163, 2.840900},

    // Z=35: Br (Bromine)
    {17.1789, 2.17230, 5.235800, 16.5796, 5.637700, 0.26090, 3.985100, 41.4328, 2.955700},

    // Z=36: Kr (Krypton)
    {17.3555, 1.93840, 6.728600, 16.5623, 5.549300, 0.22610, 3.537500, 39.3972, 2.825000},

    // Z=37: Rb (Rubidium)
    {17.1784, 1.78880, 9.643500, 17.3151, 5.139900, 0.27480, 1.529200, 164.934, 3.487300},

    // Z=38: Sr (Strontium)
    {17.5663, 1.55640, 9.818400, 14.0988, 5.422000, 0.16640, 2.669400, 132.376, 2.506400},

    // Z=39: Y (Yttrium)
    {17.7760, 1.40290, 10.2946, 12.8006, 5.726290, 0.12560, 3.265880, 104.354, 1.912130},

    // Z=40: Zr (Zirconium)
    {17.8765, 1.27618, 10.9480, 11.9160, 5.417320, 0.11762, 3.657210, 87.6627, 2.069290},

    // Z=41: Nb (Niobium)
    {17.6142, 1.18865, 12.0144, 11.7660, 4.041830, 0.20478, 3.533460, 69.7957, 3.755910},

    // Z=42: Mo (Molybdenum)
    {3.70250, 0.27720, 17.2356, 1.09580, 12.8876, 11.0040, 3.742900, 61.6584, 4.387500},

    // Z=43: Tc (Technetium)
    {19.1301, 0.86413, 11.0948, 8.14487, 4.649010, 21.5707, 2.712630, 86.8472, 5.404280},

    // Z=44: Ru (Ruthenium)
    {19.2674, 0.80852, 12.9182, 8.43467, 4.863370, 24.7997, 1.567560, 94.2928, 5.378740},

    // Z=45: Rh (Rhodium)
    {19.2957, 0.75155, 14.3501, 8.21758, 4.734250, 25.8749, 1.289180, 98.6062, 5.328000},

    // Z=46: Pd (Palladium)
    {19.3319, 0.69866, 15.5017, 7.98929, 5.295370, 25.2052, 0.605844, 76.8986, 5.265930},

    // Z=47: Ag (Silver)
    {19.2808, 0.64460, 16.6885, 7.47260, 4.804500, 24.6605, 1.046300, 99.8156, 5.179000},

    // Z=48: Cd (Cadmium)
    {19.2214, 0.59460, 17.6444, 6.90890, 4.461000, 24.7008, 1.602900, 87.4825, 5.069400},

    // Z=49: In (Indium)
    {19.1624, 0.54760, 18.5596, 6.37760, 4.294800, 25.8499, 2.039600, 92.8029, 4.939100},

    // Z=50: Sn (Tin)
    {19.1889, 5.83030, 19.1005, 0.50310, 4.458500, 26.8909, 2.466300, 83.9571, 4.782100},

    // Z=51: Sb (Antimony)
    {19.6418, 5.30340, 19.0455, 0.46070, 5.037100, 27.9074, 2.682700, 75.2825, 4.590900},

    // Z=52: Te (Tellurium)
    {19.9644, 4.81742, 19.0138, 0.42080, 6.144870, 28.5284, 2.523900, 70.8403, 4.352000},

    // Z=53: I (Iodine)
    {20.1472, 4.34700, 18.9949, 0.38140, 7.513800, 27.7660, 2.273500, 66.8776, 4.071200},

    // Z=54: Xe (Xenon)
    {20.2933, 3.92820, 19.0298, 0.34400, 8.976700, 26.4659, 1.990000, 64.2658, 3.711800}
}};

/*! \brief Get atomic form factor using Cromer-Mann parameterization - Claude Generated 2025
 *
 * Calculates the atomic scattering form factor f(q) for X-ray scattering.
 *
 * Formula: f(q) = c + Σᵢ₌₁⁴ aᵢ exp(-bᵢ(q/4π)²)
 *
 * \param Z Atomic number (1-54 supported, H to Xe)
 * \param q Scattering vector magnitude in Å⁻¹ (typically 0.01 - 10.0)
 * \return Form factor f(q) in electron units
 *
 * Reference: International Tables for Crystallography, Vol. C, Table 6.1.1.4
 * Accuracy: ~1% for q < 2 Å⁻¹, degrades at higher q
 *
 * Educational Note:
 * - Form factor represents how an atom scatters X-rays
 * - Decreases with q as interference effects reduce scattering
 * - At q=0, f(0) ≈ Z (number of electrons)
 */
inline double getAtomicFormFactor(int Z, double q) {
    // Bounds checking
    if (Z < 1 || Z > 54) {
        return 0.0;  // Out of range, return zero scattering
    }

    const auto& p = CROMER_MANN_DATA[Z];

    // Calculate s = q/(4π) for Cromer-Mann formula
    const double s = q / (4.0 * M_PI);
    const double s2 = s * s;

    // Cromer-Mann 9-parameter formula
    // f(q) = c + Σᵢ₌₁⁴ aᵢ exp(-bᵢs²)
    return p.c +
           p.a1 * std::exp(-p.b1 * s2) +
           p.a2 * std::exp(-p.b2 * s2) +
           p.a3 * std::exp(-p.b3 * s2) +
           p.a4 * std::exp(-p.b4 * s2);
}

/*! \brief Get CG bead form factor for uniform sphere - Claude Generated 2025
 *
 * Calculates the form factor for a uniform density sphere (typical CG bead model).
 *
 * Formula: P(q) = [3(sin(qR) - qR·cos(qR))/(qR)³]²
 * Simplified: f(q) = 3(sin(qR) - qR·cos(qR))/(qR)³
 *
 * \param radius Sphere radius in Ångström (typically 2-5 Å for CG beads)
 * \param q Scattering vector magnitude in Å⁻¹
 * \return Form factor f(q) normalized to f(0)=1
 *
 * Educational Note:
 * - Models CG bead as homogeneous sphere
 * - First minimum at qR ≈ 4.49 (characteristic of sphere)
 * - Oscillates at higher q (form factor oscillations)
 * - Small q limit: f(q) ≈ 1 - (qR)²/10
 *
 * Physical Interpretation:
 * - Radius determines q-range of flat scattering
 * - Larger beads → scattering drops at smaller q
 * - Guinier region: ln(P(q)) ≈ -(qRg)²/3 where Rg = √(3/5)R
 */
inline double getCGSphereFormFactor(double radius, double q) {
    const double qR = q * radius;

    // Small qR limit: use Taylor expansion to avoid 0/0
    // f(q) ≈ 1 - (qR)²/10 + (qR)⁴/280 - ...
    if (qR < 1e-4) {
        return 1.0;
    }

    // Standard formula: f(q) = 3(sin(qR) - qR·cos(qR))/(qR)³
    const double sin_qR = std::sin(qR);
    const double cos_qR = std::cos(qR);
    const double qR3 = qR * qR * qR;

    return 3.0 * (sin_qR - qR * cos_qR) / qR3;
}

/*! \brief Get form factor amplitude for either atomic or CG element - Claude Generated 2025
 *
 * Convenience function that dispatches to correct form factor calculation.
 * Auto-detects CG element (226) vs atomic elements.
 *
 * \param element Element number (1-54 for atomic, 226 for CG)
 * \param q Scattering vector magnitude in Å⁻¹
 * \param cg_radius Radius for CG sphere (used only if element==226)
 * \return Form factor f(q)
 *
 * Educational Note:
 * - Simplifies code by hiding CG vs atomic distinction
 * - Caller doesn't need to check element type
 */
inline double getFormFactor(int element, double q, double cg_radius = 3.0) {
    if (element == CG_ELEMENT) {
        return getCGSphereFormFactor(cg_radius, q);
    } else {
        return getAtomicFormFactor(element, q);
    }
}

} // namespace FormFactors
