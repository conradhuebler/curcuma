#pragma once

#include "LofthusOverlap.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <array>

/**
 * @brief Namespace containing functions for calculating Slater-Type Orbital integrals
 *
 * This namespace provides functions to calculate overlap integrals between
 * Slater-Type Orbitals (STOs) for quantum chemistry calculations.
 * Implements s, p, and d orbital overlaps with correct mathematical formulas.
 *
 * References:
 * - Mulliken, R. S., Rieke, C. A., Orloff, D., & Orloff, H. (1949).
 *   J. Chem. Phys., 17, 1248.
 * - Roothaan, C. C. J. (1951). J. Chem. Phys., 19, 1445.
 */
namespace STO {

/**
 * @brief Calculates factorial of an integer
 * @param n Input integer (must be non-negative)
 * @return n!
 */
static inline double factorial(int n)
{
    if (n < 0) {
        throw std::invalid_argument("Factorial requires non-negative integer");
    }
    if (n <= 1) {
        return 1.0;
    }
    return n * factorial(n - 1);
}

/**
 * @brief Calculates binomial coefficient (n choose k)
 * @param n Total number of items
 * @param k Number of items to choose
 * @return Binomial coefficient
 */
static inline double binomial(int n, int k)
{
    if (k < 0 || k > n) {
        return 0.0;
    }

    // Optimize calculation for large values
    if (k > n - k) {
        k = n - k;
    }

    double result = 1.0;
    for (int i = 0; i < k; ++i) {
        result *= (n - i);
        result /= (i + 1);
    }

    return result;
}

/**
 * @brief Direction cosines represent normalized direction vector components
 */
struct DirectionCosines {
    double l, m, n;

    std::string toString() const
    {
        std::stringstream ss;
        ss << "(" << std::fixed << std::setprecision(6)
           << l << ", " << m << ", " << n << ")";
        return ss.str();
    }
};

/**
 * @brief Enumerate orbital types
 */
enum OrbitalType {
    S = 0, // s orbital
    PX = 1, // px orbital
    PY = 2, // py orbital
    PZ = 3, // pz orbital
    DXY = 4, // dxy orbital
    DYZ = 5, // dyz orbital
    DZX = 6, // dzx orbital
    DX2Y2 = 7, // dx²-y² orbital
    DZ2 = 8 // dz² orbital
};

/**
 * @brief Convert orbital type to string representation
 * @param type Orbital type
 * @return String name of the orbital
 */
static inline std::string orbitalTypeToString(OrbitalType type)
{
    switch (type) {
    case S:
        return "S";
    case PX:
        return "PX";
    case PY:
        return "PY";
    case PZ:
        return "PZ";
    case DXY:
        return "DXY";
    case DYZ:
        return "DYZ";
    case DZX:
        return "DZX";
    case DX2Y2:
        return "DX2Y2";
    case DZ2:
        return "DZ2";
    default:
        return "Unknown";
    }
}

/**
 * @brief Represents an atomic orbital
 *
 * Claude Generated (January 2025): Added shell field for GFN2 shell-resolved calculations
 * Shell indices: 0=s, 1=p, 2=d (matches TBLite convention)
 */
struct Orbital {
    double x = 0, y = 0, z = 0; // Position
    double zeta; // Slater exponent
    double VSIP; // Valence State Ionization Potential
    OrbitalType type; // Orbital type
    int atom; // Atom index
    int shell = 0; // Shell index (0=s, 1=p, 2=d) for shell-resolved calculations
    int principal_n = 1; // Principal quantum number (Claude Generated March 2026)

    std::string toString() const
    {
        std::stringstream ss;
        ss << orbitalTypeToString(type) << " orbital at ("
           << std::fixed << std::setprecision(6)
           << x << ", " << y << ", " << z
           << ") with zeta=" << zeta << ", VSIP=" << VSIP
           << ", shell=" << shell << ", n=" << principal_n;
        return ss.str();
    }
};

/**
 * @brief Calculate direction cosines between two orbitals
 * @param orb1 First orbital
 * @param orb2 Second orbital
 * @return Direction cosines (normalized direction vector)
 */
static inline DirectionCosines getDirectionCosines(const Orbital& orb1, const Orbital& orb2)
{
    double dx = orb2.x - orb1.x;
    double dy = orb2.y - orb1.y;
    double dz = orb2.z - orb1.z;
    double R = std::sqrt(dx * dx + dy * dy + dz * dz);

    // Handle the case where atoms are very close to each other
    if (R < 1e-10) {
        // Default direction for degenerate case (arbitrary choice)
        return { 0.0, 0.0, 1.0 };
    }

    // Return the normalized direction vector
    return { dx / R, dy / R, dz / R };
}

/**
 * @brief Calculate overlap between two s orbitals
 * @param zeta1 Slater exponent of first orbital
 * @param zeta2 Slater exponent of second orbital
 * @param R Distance between orbital centers
 * @param debug Whether to print debug information
 * @return Overlap integral value
 *
 * Formula: S_{ss} = (2√(ζ₁ζ₂))/(ζ₁+ζ₂)^(3/2) · e^(-ρ) · (1 + ρ + ρ²/3)
 * where ρ = R(ζ₁+ζ₂)/2
 */
static inline double calculateSSOverlap(double zeta1, double zeta2, double R_angstrom, bool debug = false)
{
    // Claude Generated: 1s-1s STO overlap integral
    // S = N * exp(-p*R) * (1 + p*R + (p*R)²/3)
    // where p = (ζ₁+ζ₂)/2, R in Bohr, N = [4ζ₁ζ₂/(ζ₁+ζ₂)²]^(3/2)
    // Reference: Mulliken et al., J. Chem. Phys. 17, 1248 (1949)

    constexpr double bohr = 0.52917721092;  // Angstrom per Bohr
    double R = R_angstrom / bohr;  // Convert to Bohr

    // Average exponent p = (ζ₁+ζ₂)/2 — NOT reduced mass ζ₁ζ₂/(ζ₁+ζ₂)
    double p = (zeta1 + zeta2) / 2.0;
    double pR = p * R;

    // Normalization factor
    double N = pow(4.0 * zeta1 * zeta2 / pow(zeta1 + zeta2, 2.0), 1.5);

    double exp_term = exp(-pR);
    double poly_term = 1.0 + pR + pR * pR / 3.0;
    double result = N * exp_term * poly_term;

    if (debug) {
        std::cout << "S-S overlap: zeta1=" << zeta1 << ", zeta2=" << zeta2
                  << ", R_ang=" << R_angstrom << ", R_bohr=" << R
                  << ", p=" << p << ", N=" << N << ", exp=" << exp_term
                  << ", poly=" << poly_term
                  << ", result=" << result << std::endl;
    }

    return result;
}

/**
 * @brief Calculate SS overlap with principal quantum number correction
 * Claude Generated (March 2026): Corrects for different n values (e.g. 1s-2s).
 *
 * The base formula computes the 1s-1s overlap (n1=n2=1). For different n,
 * we scale by the ratio of exact same-center overlaps:
 *   S(n1, n2, R) ≈ S_1s1s(R) * S_exact(n1, n2, R=0) / S_1s1s(R=0)
 *
 * Same-center exact formula:
 *   S(n1, n2, ζ1, ζ2, R=0) = sqrt((2ζ1)^(2n1+1)/(2n1)!) * sqrt((2ζ2)^(2n2+1)/(2n2)!)
 *                              * (n1+n2)! / (ζ1+ζ2)^(n1+n2+1)
 */
static inline double calculateSSOverlapN(double zeta1, double zeta2, double R,
                                          int n1, int n2, bool debug = false)
{
    // 1s-1s overlap (base)
    double S_1s1s = calculateSSOverlap(zeta1, zeta2, R, debug);

    if (n1 == 1 && n2 == 1)
        return S_1s1s;

    // Exact same-center overlap for (n1, n2):
    // N_n = sqrt((2ζ)^(2n+1) / (2n)!)
    // S(R=0) = N1 * N2 * (n1+n2)! / (ζ1+ζ2)^(n1+n2+1)
    double N1 = std::sqrt(std::pow(2.0 * zeta1, 2*n1 + 1) / factorial(2*n1));
    double N2 = std::sqrt(std::pow(2.0 * zeta2, 2*n2 + 1) / factorial(2*n2));
    double S_exact_R0 = N1 * N2 * factorial(n1 + n2) / std::pow(zeta1 + zeta2, n1 + n2 + 1);

    // Same-center overlap for (1, 1):
    double N1_1s = std::sqrt(std::pow(2.0 * zeta1, 3) / 2.0);
    double N2_1s = std::sqrt(std::pow(2.0 * zeta2, 3) / 2.0);
    double S_1s1s_R0 = N1_1s * N2_1s * 2.0 / std::pow(zeta1 + zeta2, 3);

    // Correction factor
    double correction = (S_1s1s_R0 > 1e-15) ? S_exact_R0 / S_1s1s_R0 : 1.0;

    return S_1s1s * correction;
}

/**
 * @brief Calculate overlap between s and p orbitals
 * @param zeta_s Slater exponent of s orbital
 * @param zeta_p Slater exponent of p orbital
 * @param R Distance between orbital centers
 * @param l Direction cosine component for p orbital
 * @param debug Whether to print debug information
 * @return Overlap integral value
 *
 * Formula: S_{sp} = (2√(ζₛζₚ))/(ζₛ+ζₚ)^(3/2) · λ · e^(-ρ) · ρ(1 + ρ/3)
 * where ρ = R(ζₛ+ζₚ)/2 and λ is the direction cosine
 */
static inline double calculateSPOverlap(double zeta1, double zeta2, double R_angstrom,
    double direction, bool debug = false)
{
    // Claude Generated: 1s-2p_σ STO overlap integral
    // S_sp = direction * N * pR * exp(-pR) * (1 + pR/3)
    // where p = (ζ_s+ζ_p)/2, R in Bohr
    // Reference: Mulliken et al., J. Chem. Phys. 17, 1248 (1949)

    if (R_angstrom < 1e-10)
        return 0.0;

    constexpr double bohr = 0.52917721092;
    double R = R_angstrom / bohr;

    double p = (zeta1 + zeta2) / 2.0;
    double pR = p * R;

    // Normalization: N = √(ζ_s ζ_p) * [2/(ζ_s+ζ_p)]^(3/2)
    double N = sqrt(zeta1 * zeta2) * pow(2.0 / (zeta1 + zeta2), 1.5);

    double exp_term = exp(-pR);
    double S_sp = N * pR * exp_term * (1.0 + pR / 3.0);
    double result = direction * S_sp;

    if (debug) {
        std::cout << "S-P overlap: zeta1=" << zeta1 << ", zeta2=" << zeta2
                  << ", R_ang=" << R_angstrom << ", R_bohr=" << R
                  << ", direction=" << direction
                  << ", p=" << p << ", N=" << N
                  << ", exp=" << exp_term
                  << ", S_sp=" << S_sp
                  << ", result=" << result << std::endl;
    }

    return result;
}

/**
 * @brief Calculate overlap between two p orbitals
 * @param zeta1 Slater exponent of first orbital
 * @param zeta2 Slater exponent of second orbital
 * @param R Distance between orbital centers
 * @param cos_angle Cosine of angle between p orbital directions
 * @param same_axis Whether p orbitals are along the same axis
 * @param debug Whether to print debug information
 * @return Overlap integral value
 */
static inline double calculatePPOverlap(double zeta1, double zeta2, double R_angstrom,
    double cos_angle, bool same_axis,
    bool debug = false)
{
    // Claude Generated: 2p-2p STO overlap integrals (σ and π components)
    // S_pp_σ = N * exp(-pR) * (1 - pR + (pR)²/3)
    // S_pp_π = N * exp(-pR) * (1 + pR + (pR)²/3)
    // where p = (ζ₁+ζ₂)/2, R in Bohr
    // Reference: Mulliken et al., J. Chem. Phys. 17, 1248 (1949)

    if (R_angstrom < 1e-10)
        return same_axis ? 1.0 : 0.0;

    constexpr double bohr = 0.52917721092;
    double R = R_angstrom / bohr;

    double p = (zeta1 + zeta2) / 2.0;
    double pR = p * R;

    // Normalization: N = ζ₁ζ₂ * [2/(ζ₁+ζ₂)]³
    double N = zeta1 * zeta2 * pow(2.0 / (zeta1 + zeta2), 3.0);

    double exp_term = exp(-pR);
    double result;

    if (same_axis) {
        double poly = 1.0 - pR + pR * pR / 3.0;
        result = N * exp_term * poly;
    } else if (fabs(cos_angle) < 1e-10) {
        double poly = 1.0 + pR + pR * pR / 3.0;
        result = N * exp_term * poly;
    } else {
        double poly_sigma = 1.0 - pR + pR * pR / 3.0;
        double poly_pi = 1.0 + pR + pR * pR / 3.0;
        double cos2 = cos_angle * cos_angle;
        result = N * exp_term * (cos2 * poly_sigma + (1.0 - cos2) * poly_pi);
    }

    // Debug-Ausgabe wenn angefordert
    if (debug) {
        std::cout << "P-P overlap: zeta1=" << zeta1 << ", zeta2=" << zeta2
                  << ", R=" << R << ", cos_angle=" << cos_angle
                  << ", same_axis=" << (same_axis ? "true" : "false")
                  << ", p=" << p
                  << ", N=" << N << ", exp=" << exp_term
                  << ", result=" << result << std::endl;
    }

    return result;
}

/**
 * @brief Calculate overlap between s and d orbitals
 * @param zeta_s Slater exponent of s orbital
 * @param zeta_d Slater exponent of d orbital
 * @param R Distance between orbital centers
 * @param dc Direction cosines
 * @param d_type Type of d orbital
 * @param debug Whether to print debug information
 * @return Overlap integral value
 */
static inline double calculateSDOverlap(double zeta_s, double zeta_d, double R,
    const DirectionCosines& dc, OrbitalType d_type,
    bool debug = false)
{
    // For identical centers (R=0)
    if (R < 1e-10)
        return 0.0;

    // Calculate rho parameter
    double rho = R * (zeta_s + zeta_d) / 2.0;

    // Calculate normalization factor
    double numerator = 2.0 * std::sqrt(zeta_s * zeta_d);
    double denominator = std::pow(zeta_s + zeta_d, 1.5);
    double N = numerator / denominator;

    // Calculate exponential and polynomial terms
    double exp_term = std::exp(-rho);
    double poly_term = rho * rho * (1.0 + rho / 5.0);

    // Calculate angular term based on d orbital type
    double angular_term = 0.0;

    switch (d_type) {
    case DXY:
        angular_term = std::sqrt(3.0) * dc.l * dc.m;
        break;
    case DYZ:
        angular_term = std::sqrt(3.0) * dc.m * dc.n;
        break;
    case DZX:
        angular_term = std::sqrt(3.0) * dc.n * dc.l;
        break;
    case DX2Y2:
        angular_term = (std::sqrt(3.0) / 2.0) * (dc.l * dc.l - dc.m * dc.m);
        break;
    case DZ2:
        angular_term = (3.0 * dc.n * dc.n - 1.0) / 2.0;
        break;
    default:
        throw std::invalid_argument("Invalid d orbital type for S-D overlap");
    }

    // Final result
    double result = N * exp_term * poly_term * angular_term;

    // Debug output if requested
    if (debug) {
        std::cout << "S-D overlap: zeta_s=" << zeta_s << ", zeta_d=" << zeta_d
                  << ", R=" << R << ", d_type=" << orbitalTypeToString(d_type)
                  << ", dc=" << dc.toString()
                  << ", rho=" << rho << ", N=" << N
                  << ", exp=" << exp_term
                  << ", poly=" << poly_term
                  << ", angular=" << angular_term
                  << ", result=" << result << std::endl;
    }

    return result;
}

/**
 * @brief Calculate overlap between p and d orbitals
 * @param zeta_p Slater exponent of p orbital
 * @param zeta_d Slater exponent of d orbital
 * @param R Distance between orbital centers
 * @param dc Direction cosines
 * @param p_type Type of p orbital
 * @param d_type Type of d orbital
 * @param debug Whether to print debug information
 * @return Overlap integral value
 */
static inline double calculatePDOverlap(double zeta_p, double zeta_d, double R,
    const DirectionCosines& dc,
    OrbitalType p_type, OrbitalType d_type,
    bool debug = false)
{
    // For identical centers (R=0)
    if (R < 1e-10)
        return 0.0;

    // Validate p orbital type
    if (p_type < PX || p_type > PZ) {
        throw std::invalid_argument("Invalid p orbital type for P-D overlap");
    }

    // Calculate rho parameter
    double rho = R * (zeta_p + zeta_d) / 2.0;

    // Calculate normalization factor
    double numerator = 2.0 * std::sqrt(zeta_p * zeta_d);
    double denominator = std::pow(zeta_p + zeta_d, 1.5);
    double N = numerator / denominator;

    // Calculate exponential and polynomial terms
    double exp_term = std::exp(-rho);
    double poly_term = rho * (1.0 + rho / 3.0 + (rho * rho) / 15.0);

    // Determine angular term based on orbital types
    double angular_term = 0.0;

    switch (d_type) {
    case DXY:
        if (p_type == PX)
            angular_term = std::sqrt(3.0) * dc.m;
        else if (p_type == PY)
            angular_term = std::sqrt(3.0) * dc.l;
        else
            angular_term = 0.0; // PZ with DXY is zero
        break;

    case DYZ:
        if (p_type == PX)
            angular_term = 0.0; // PX with DYZ is zero
        else if (p_type == PY)
            angular_term = std::sqrt(3.0) * dc.n;
        else if (p_type == PZ)
            angular_term = std::sqrt(3.0) * dc.m;
        break;

    case DZX:
        if (p_type == PX)
            angular_term = std::sqrt(3.0) * dc.n;
        else if (p_type == PY)
            angular_term = 0.0; // PY with DZX is zero
        else if (p_type == PZ)
            angular_term = std::sqrt(3.0) * dc.l;
        break;

    case DX2Y2:
        if (p_type == PX)
            angular_term = std::sqrt(3.0) * dc.l;
        else if (p_type == PY)
            angular_term = -std::sqrt(3.0) * dc.m;
        else
            angular_term = 0.0; // PZ with DX2Y2 is zero
        break;

    case DZ2:
        if (p_type == PX)
            angular_term = -(std::sqrt(3.0) / 2.0) * dc.l;
        else if (p_type == PY)
            angular_term = -(std::sqrt(3.0) / 2.0) * dc.m;
        else if (p_type == PZ)
            angular_term = std::sqrt(3.0) * dc.n;
        break;

    default:
        throw std::invalid_argument("Invalid d orbital type for P-D overlap");
    }

    // Final result
    double result = N * exp_term * poly_term * angular_term;

    // Debug output if requested
    if (debug) {
        std::cout << "P-D overlap: p_type=" << orbitalTypeToString(p_type)
                  << ", d_type=" << orbitalTypeToString(d_type)
                  << ", zeta_p=" << zeta_p << ", zeta_d=" << zeta_d
                  << ", R=" << R << ", dc=" << dc.toString()
                  << ", rho=" << rho << ", N=" << N
                  << ", exp=" << exp_term
                  << ", poly=" << poly_term
                  << ", angular=" << angular_term
                  << ", result=" << result << std::endl;
    }

    return result;
}

/**
 * @brief Calculate overlap between two d orbitals
 * @param zeta1 Slater exponent of first orbital
 * @param zeta2 Slater exponent of second orbital
 * @param R Distance between orbital centers
 * @param dc Direction cosines
 * @param d_type1 Type of first d orbital
 * @param d_type2 Type of second d orbital
 * @param debug Whether to print debug information
 * @return Overlap integral value
 */
static inline double calculateDDOverlap(double zeta1, double zeta2, double R,
    const DirectionCosines& dc,
    OrbitalType d_type1, OrbitalType d_type2,
    bool debug = false)
{
    // For identical centers (R=0)
    if (R < 1e-10)
        return (d_type1 == d_type2) ? 1.0 : 0.0;

    // Validate d orbital types
    if (d_type1 < DXY || d_type1 > DZ2 || d_type2 < DXY || d_type2 > DZ2) {
        throw std::invalid_argument("Invalid d orbital type for D-D overlap");
    }

    // Calculate rho parameter
    double rho = R * (zeta1 + zeta2) / 2.0;

    // Calculate normalization factor
    double numerator = 2.0 * std::sqrt(zeta1 * zeta2);
    double denominator = std::pow(zeta1 + zeta2, 1.5);
    double N = numerator / denominator;

    // Calculate exponential and polynomial terms
    double exp_term = std::exp(-rho);
    double poly_term = (1.0 + rho + (rho * rho) / 5.0 + (rho * rho * rho) / 35.0);

    // Calculate angular term based on d orbital types
    double angular_term = 0.0;

    // Same d orbital types
    if (d_type1 == d_type2) {
        switch (d_type1) {
        case DXY:
            angular_term = (3.0 * dc.l * dc.l * dc.m * dc.m - (rho * rho) / 35.0);
            break;
        case DYZ:
            angular_term = (3.0 * dc.m * dc.m * dc.n * dc.n - (rho * rho) / 35.0);
            break;
        case DZX:
            angular_term = (3.0 * dc.n * dc.n * dc.l * dc.l - (rho * rho) / 35.0);
            break;
        case DX2Y2:
            angular_term = (3.0 * std::pow(dc.l * dc.l - dc.m * dc.m, 2) / 4.0 - (rho * rho) / 35.0);
            break;
        case DZ2:
            angular_term = (std::pow(3.0 * dc.n * dc.n - 1.0, 2) / 4.0 - (rho * rho) / 35.0);
            break;
        default:
            break; // Should never reach here due to earlier validation
        }
    }
    // Different d orbital types - only implementing the most common ones
    else if ((d_type1 == DXY && d_type2 == DX2Y2) || (d_type1 == DX2Y2 && d_type2 == DXY)) {
        angular_term = (3.0 * std::sqrt(3.0) / 4.0) * dc.l * dc.m * (dc.l * dc.l - dc.m * dc.m);
    } else if ((d_type1 == DXY && d_type2 == DZ2) || (d_type1 == DZ2 && d_type2 == DXY)) {
        angular_term = (3.0 * std::sqrt(3.0) / 4.0) * dc.l * dc.m * (3.0 * dc.n * dc.n - 1.0);
    } else if ((d_type1 == DYZ && d_type2 == DZ2) || (d_type1 == DZ2 && d_type2 == DYZ)) {
        angular_term = (3.0 * std::sqrt(3.0) / 4.0) * dc.m * dc.n * (3.0 * dc.n * dc.n - 1.0);
    } else if ((d_type1 == DZX && d_type2 == DZ2) || (d_type1 == DZ2 && d_type2 == DZX)) {
        angular_term = (3.0 * std::sqrt(3.0) / 4.0) * dc.l * dc.n * (3.0 * dc.n * dc.n - 1.0);
    } else if ((d_type1 == DX2Y2 && d_type2 == DZ2) || (d_type1 == DZ2 && d_type2 == DX2Y2)) {
        angular_term = (3.0 * std::sqrt(3.0) / 4.0) * (dc.l * dc.l - dc.m * dc.m) * (3.0 * dc.n * dc.n - 1.0);
    }

    // Final result
    double result = N * exp_term * poly_term * angular_term;

    // Debug output if requested
    if (debug) {
        std::cout << "D-D overlap: d_type1=" << orbitalTypeToString(d_type1)
                  << ", d_type2=" << orbitalTypeToString(d_type2)
                  << ", zeta1=" << zeta1 << ", zeta2=" << zeta2
                  << ", R=" << R << ", dc=" << dc.toString()
                  << ", rho=" << rho << ", N=" << N
                  << ", exp=" << exp_term
                  << ", poly=" << poly_term
                  << ", angular=" << angular_term
                  << ", result=" << result << std::endl;
    }

    return result;
}

/**
 * @brief Calculate overlap between two orbitals of any supported type
 * @param orb1 First orbital
 * @param orb2 Second orbital
 * @param debug Whether to print debug information
 * @return Overlap integral value
 */
static inline double calculateOverlap(const Orbital& orb1, const Orbital& orb2, bool debug = false)
{
    // Calculate distance between orbitals
    double dx = orb2.x - orb1.x;
    double dy = orb2.y - orb1.y;
    double dz = orb2.z - orb1.z;
    double R = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (debug) {
        std::cout << "Calculating overlap between orbitals: " << std::endl
                  << "  Orbital 1: " << orb1.toString() << std::endl
                  << "  Orbital 2: " << orb2.toString() << std::endl
                  << "  Distance R = " << R << std::endl;
    }

    // Get direction cosines
    DirectionCosines dc = getDirectionCosines(orb1, orb2);

    if (debug) {
        std::cout << "Direction cosines: " << dc.toString() << std::endl;
    }

    // Special case for identical orbitals
    if (R < 1e-10 && orb1.type == orb2.type && std::abs(orb1.zeta - orb2.zeta) < 1e-10) {
        if (debug) {
            std::cout << "Identical orbitals detected - returning 1.0" << std::endl;
        }
        return 1.0;
    }

    OrbitalType type1 = orb1.type;
    OrbitalType type2 = orb2.type;
    double zeta1 = orb1.zeta;
    double zeta2 = orb2.zeta;
    int n1 = orb1.principal_n;
    int n2 = orb2.principal_n;

    double result = 0.0;

    // Claude Generated (April 2026): Exact Lofthus/Pople STO overlaps for s/p orbitals
    // Uses confocal elliptical coordinates with Ak/Bk auxiliary integrals.
    // IMPORTANT: argB = 0.5*R*(ζ_bra - ζ_ket) is asymmetric, so orbital ordering
    // must match Lofthus convention: bra=orb1 (left), ket=orb2 (right).
    // Ref: E. Lofthus, Mol. Phys. 5, 105 (1962); Pople & Beveridge, "Approx. MO Theory" (1970)
    constexpr double bohr = 0.52917721092;
    double R_bohr = R / bohr;

    // Determine orbital category pair (without swapping — order matters for Lofthus)
    bool is_s1 = (type1 == S), is_s2 = (type2 == S);
    bool is_p1 = (type1 >= PX && type1 <= PZ), is_p2 = (type2 >= PX && type2 <= PZ);

    if (is_s1 && is_s2) {
        // SS overlap: symmetric in orbital order
        result = Lofthus::exactSSOverlap(n1, zeta1, n2, zeta2, R_bohr);
    }
    else if (is_s1 && is_p2) {
        // SP case: S on orb1 (bra), P on orb2 (ket)
        // Following Ulysses: result = -SProt(k,1) * sigma_SP
        // SProt(k,1) = direction cosine for component k
        double direction = 0.0;
        if (type2 == PX) direction = dc.l;
        else if (type2 == PY) direction = dc.m;
        else direction = dc.n;

        // exactSPSigmaOverlap(n_s, zeta_s, n_p, zeta_p, R)
        double sp_sigma = Lofthus::exactSPSigmaOverlap(n1, zeta1, n2, zeta2, R_bohr);
        result = -sp_sigma * direction;  // negative sign per Ulysses convention
    }
    else if (is_p1 && is_s2) {
        // PS case: P on orb1 (bra), S on orb2 (ket)
        // Following Ulysses: result = +SProt(k,1) * sigma_PS
        // Uses separate PSOvIntIndex polynomial with [1,+1] factor
        double direction = 0.0;
        if (type1 == PX) direction = dc.l;
        else if (type1 == PY) direction = dc.m;
        else direction = dc.n;

        // exactPSSigmaOverlap(n_p, zeta_p, n_s, zeta_s, R)
        double ps_sigma = Lofthus::exactPSSigmaOverlap(n1, zeta1, n2, zeta2, R_bohr);
        result = ps_sigma * direction;  // positive sign per Ulysses convention
    }
    else if (is_p1 && is_p2) {
        // PP overlap: combine sigma and pi components with rotation
        double pp_sigma = Lofthus::exactPPSigmaOverlap(n1, zeta1, n2, zeta2, R_bohr);
        double pp_pi = Lofthus::exactPPPiOverlap(n1, zeta1, n2, zeta2, R_bohr);

        // Rotation from diatomic to lab frame using SProt matrix
        // S(i,j) = -sigma*SProt(i,0)*SProt(j,0) + pi*(SProt(i,1)*SProt(j,1) + SProt(i,2)*SProt(j,2))
        double cost = dc.n;
        double sint = std::sqrt(dc.l * dc.l + dc.m * dc.m);
        double cosp = 1.0, sinp = 0.0;
        if (sint > 1e-8) {
            cosp = dc.l / sint;
            sinp = dc.m / sint;
        }

        double SP[3][3] = {
            {sint * cosp, cost * cosp, -sinp},
            {sint * sinp, cost * sinp,  cosp},
            {cost,        -sint,         0.0 }
        };

        int i1 = type1 - PX;
        int i2 = type2 - PX;

        result = -pp_sigma * SP[i1][0] * SP[i2][0]
                + pp_pi * (SP[i1][1] * SP[i2][1] + SP[i1][2] * SP[i2][2]);
    }
    else if ((is_s1 || is_p1) && !(is_s2 || is_p2)) {
        // SD or PD overlap: use old approximate formulas
        if (is_s1) {
            result = calculateSDOverlap(zeta1, zeta2, R, dc, type2, debug);
        } else {
            result = calculatePDOverlap(zeta1, zeta2, R, dc, type1, type2, debug);
        }
    }
    else if (!(is_s1 || is_p1) && (is_s2 || is_p2)) {
        // DS or DP overlap: swap and use old formulas
        if (is_s2) {
            result = calculateSDOverlap(zeta2, zeta1, R, dc, type1, debug);
        } else {
            result = calculatePDOverlap(zeta2, zeta1, R, dc, type2, type1, debug);
        }
    }
    else {
        // DD overlap
        result = calculateDDOverlap(zeta1, zeta2, R, dc, type1, type2, debug);
    }

    if (debug) {
        std::cout << "Final overlap result: " << result << std::endl
                  << std::endl;
    }

    return result;
}

/**
 * @brief Calculate radial derivative of overlap with respect to distance R (dS/dR)
 * @param orb1 First orbital
 * @param orb2 Second orbital
 * @param R Distance between orbitals
 * @return dS/dR
 *
 * Claude Generated (February 2026): For analytical gradients
 */
static inline double calculateOverlapDerivative(const Orbital& orb1, const Orbital& orb2, double R)
{
    // Claude Generated: Overlap derivatives with correct p = (ζ₁+ζ₂)/2 and R in Bohr
    // dS/dR_angstrom = dS/dR_bohr * dR_bohr/dR_ang = dS/dR_bohr / bohr
    if (R < 1e-10) return 0.0;

    constexpr double bohr = 0.52917721092;
    double R_bohr = R / bohr;

    double zeta1 = orb1.zeta;
    double zeta2 = orb2.zeta;
    OrbitalType type1 = orb1.type;
    OrbitalType type2 = orb2.type;

    if (type1 > type2) {
        std::swap(type1, type2);
        std::swap(zeta1, zeta2);
    }

    double p = (zeta1 + zeta2) / 2.0;
    double pR = p * R_bohr;
    double exp_term = std::exp(-pR);

    // All derivatives are dS/dR_angstrom, so divide by bohr at the end
    if (type1 == S && type2 == S) {
        double N = std::pow(4.0 * zeta1 * zeta2 / std::pow(zeta1 + zeta2, 2.0), 1.5);
        // d/dR[exp(-pR)(1+pR+p²R²/3)] = exp(-pR)*(-p²R²/3 * p - p + p + p²R - p³R²/3)
        // = exp(-pR) * p * (-pR/3 - p²R²/3) = -p*exp(-pR)*pR*(1/3 + pR/3)
        return N * p * exp_term * (-pR / 3.0 - pR * pR / 3.0) / bohr;
    }
    else if (type1 == S && type2 <= PZ) {
        double N = std::sqrt(zeta1 * zeta2) * std::pow(2.0 / (zeta1 + zeta2), 1.5);
        // d/dR[pR*exp(-pR)*(1+pR/3)] = p*exp(-pR)*(1+pR/3) + pR*exp(-pR)*(-p*(1+pR/3) + p/3)
        // = p*exp(-pR)*(1 - pR/3 - p²R²/3)
        return N * p * exp_term * (1.0 - pR / 3.0 - pR * pR / 3.0) / bohr;
    }
    else if (type1 <= PZ && type2 <= PZ) {
        double N = zeta1 * zeta2 * std::pow(2.0 / (zeta1 + zeta2), 3.0);
        if (type1 == type2) { // Sigma: d/dR[exp(-pR)(1-pR+p²R²/3)]
            return N * exp_term * (-2.0*p + 5.0/3.0*p*pR - p*pR*pR/3.0) / bohr;
        } else { // Pi: d/dR[exp(-pR)(1+pR+p²R²/3)]
            return N * p * exp_term * (-pR/3.0 - pR*pR/3.0) / bohr;
        }
    }

    return 0.0;
}

/**
 * @brief Calculate full gradient of overlap with respect to position of orb1
 * @param orb1 First orbital
 * @param orb2 Second orbital
 * @return 3D gradient vector
 */
static inline std::array<double, 3> calculateOverlapFullGradient(const Orbital& orb1, const Orbital& orb2)
{
    double dx = orb2.x - orb1.x;
    double dy = orb2.y - orb1.y;
    double dz = orb2.z - orb1.z;
    double R = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (R < 1e-10) return {0.0, 0.0, 0.0};

    double dSdR = calculateOverlapDerivative(orb1, orb2, R);

    // dS/dx_A = dS/dR * dR/dx_A = dS/dR * (x_A - x_B) / R
    return { dSdR * (orb1.x - orb2.x) / R,
             dSdR * (orb1.y - orb2.y) / R,
             dSdR * (orb1.z - orb2.z) / R };
}

} // namespace STO
