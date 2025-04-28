#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

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
 */
struct Orbital {
    double x = 0, y = 0, z = 0; // Position
    double zeta; // Slater exponent
    double VSIP; // Valence State Ionization Potential
    OrbitalType type; // Orbital type
    int atom; // Atom index

    std::string toString() const
    {
        std::stringstream ss;
        ss << orbitalTypeToString(type) << " orbital at ("
           << std::fixed << std::setprecision(6)
           << x << ", " << y << ", " << z
           << ") with zeta=" << zeta << ", VSIP=" << VSIP;
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
static inline double calculateSSOverlap(double zeta1, double zeta2, double R, bool debug = false)
{
    // For identical centers (R=0)
    if (R < 1e-10)
        return 1.0;

    // Calculate rho parameter
    double rho = R * (zeta1 + zeta2) / 2.0;

    // Calculate normalization factor
    double numerator = 2.0 * std::sqrt(zeta1 * zeta2);
    double denominator = std::pow(zeta1 + zeta2, 1.5);
    double N = numerator / denominator;

    // Calculate polynomial and exponential terms
    double poly_term = 1.0 + rho + rho * rho / 3.0;
    double exp_term = std::exp(-rho);

    // Final result
    double result = N * exp_term * poly_term;

    // Debug output if requested
    if (debug) {
        std::cout << "S-S overlap: zeta1=" << zeta1 << ", zeta2=" << zeta2 << ", R=" << R
                  << ", rho=" << rho
                  << ", num=" << numerator
                  << ", denom=" << denominator
                  << ", N=" << N
                  << ", exp=" << exp_term
                  << ", poly=" << poly_term
                  << ", result=" << result << std::endl;
    }

    return result;
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
static inline double calculateSPOverlap(double zeta_s, double zeta_p, double R, double l, bool debug = false)
{
    // For identical centers (R=0)
    if (R < 1e-10)
        return 0.0;

    // Calculate rho parameter
    double rho = R * (zeta_s + zeta_p) / 2.0;

    // Calculate normalization factor
    double numerator = 2.0 * std::sqrt(zeta_s * zeta_p);
    double denominator = std::pow(zeta_s + zeta_p, 1.5);
    double N = numerator / denominator;

    // Calculate polynomial and exponential terms
    double poly_term = rho * (1.0 + rho / 3.0);
    double exp_term = std::exp(-rho);

    // Final result
    double result = N * l * exp_term * poly_term;

    // Debug output if requested
    if (debug) {
        std::cout << "S-P overlap: zeta_s=" << zeta_s << ", zeta_p=" << zeta_p
                  << ", R=" << R << ", l=" << l << ", rho=" << rho
                  << ", N=" << N
                  << ", exp=" << exp_term
                  << ", poly=" << poly_term
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
static inline double calculatePPOverlap(double zeta1, double zeta2, double R,
    double cos_angle, bool same_axis,
    bool debug = false)
{
    // For identical centers (R=0)
    if (R < 1e-10)
        return same_axis ? 1.0 : 0.0;

    // Calculate rho parameter
    double rho = R * (zeta1 + zeta2) / 2.0;

    // Calculate normalization factor
    double numerator = 2.0 * std::sqrt(zeta1 * zeta2);
    double denominator = std::pow(zeta1 + zeta2, 1.5);
    double N = numerator / denominator;

    // Calculate exponential term
    double exp_term = std::exp(-rho);
    double poly_term;
    double result;

    if (same_axis) {
        // P-P overlap, same axis
        poly_term = (cos_angle * cos_angle) * (1.0 + rho + (2.0 * rho * rho) / 5.0) - (rho * rho) / 5.0;
    } else {
        // P-P overlap, different axes
        poly_term = cos_angle * cos_angle * (1.0 + rho + (2.0 * rho * rho) / 5.0);
    }

    // Final result
    result = N * exp_term * poly_term;

    // Debug output if requested
    if (debug) {
        std::cout << "P-P overlap: zeta1=" << zeta1 << ", zeta2=" << zeta2
                  << ", R=" << R << ", cos_angle=" << cos_angle
                  << ", same_axis=" << (same_axis ? "true" : "false")
                  << ", rho=" << rho << ", N=" << N
                  << ", exp=" << exp_term
                  << ", poly=" << poly_term
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

    // Sort orbital types for consistent calculation
    OrbitalType type1 = orb1.type;
    OrbitalType type2 = orb2.type;
    double zeta1 = orb1.zeta;
    double zeta2 = orb2.zeta;

    // Always use the lower type first for symmetry
    if (type1 > type2) {
        std::swap(type1, type2);
        std::swap(zeta1, zeta2);

        if (debug) {
            std::cout << "Swapped orbital order for calculation" << std::endl;
        }
    }

    double result = 0.0;

    // Calculate appropriate overlap based on orbital types
    if (type1 == S) {
        if (type2 == S) {
            result = calculateSSOverlap(zeta1, zeta2, R, debug);
        } else if (type2 <= PZ) {
            double direction = 0.0;
            if (type2 == PX)
                direction = dc.l;
            else if (type2 == PY)
                direction = dc.m;
            else
                direction = dc.n;

            result = calculateSPOverlap(zeta1, zeta2, R, direction, debug);
        } else {
            result = calculateSDOverlap(zeta1, zeta2, R, dc, type2, debug);
        }
    } else if (type1 <= PZ) {
        if (type2 <= PZ) {
            bool same_axis = (type1 == type2);

            double l1 = 0.0, l2 = 0.0;
            if (type1 == PX)
                l1 = dc.l;
            else if (type1 == PY)
                l1 = dc.m;
            else
                l1 = dc.n;

            if (type2 == PX)
                l2 = dc.l;
            else if (type2 == PY)
                l2 = dc.m;
            else
                l2 = dc.n;

            result = calculatePPOverlap(zeta1, zeta2, R, l1 * l2, same_axis, debug);
        } else {
            result = calculatePDOverlap(zeta1, zeta2, R, dc, type1, type2, debug);
        }
    } else {
        result = calculateDDOverlap(zeta1, zeta2, R, dc, type1, type2, debug);
    }

    if (debug) {
        std::cout << "Final overlap result: " << result << std::endl
                  << std::endl;
    }

    return result;
}

} // namespace STO
