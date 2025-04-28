#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

/* Integrals for Gaussian Type Orbitals */
namespace GTO {

static inline double factorial(int n)
{
    if (n <= 1)
        return 1.0;
    return n * factorial(n - 1);
}

static inline double binomial(int n, int k)
{
    return factorial(n) / (factorial(k) * factorial(n - k));
}

// Direction cosines
struct DirectionCosines {
    double l, m, n;
};

enum OrbitalType {
    S = 0,
    PX = 1,
    PY = 2,
    PZ = 3,
    DXY = 4,
    DYZ = 5,
    DZX = 6,
    DX2Y2 = 7,
    DZ2 = 8
};

// Structure for a contracted Gaussian function
struct Orbital {
    double x = 0, y = 0, z = 0; // Position of orbital center
    double VSIP; // Valence State Ionization Potential
    OrbitalType type; // Type of orbital (s, px, py, etc.)
    std::vector<double> exponents; // Gaussian exponents
    std::vector<double> coefficients; // Contraction coefficients
    int atom; // Index of atom to which orbital belongs
};

// Get direction cosines between two orbital centers
static inline DirectionCosines getDirectionCosines(const Orbital& orb1, const Orbital& orb2)
{
    double dx = orb2.x - orb1.x;
    double dy = orb2.y - orb1.y;
    double dz = orb2.z - orb1.z;
    double R = std::sqrt(dx * dx + dy * dy + dz * dz);

    // Handle the case where the atoms are very close to each other
    if (R < 1e-10) {
        // Default direction for degenerate case (arbitrary choice)
        return { 0.0, 0.0, 1.0 };
    }

    // Return the correct direction cosines
    return { dx / R, dy / R, dz / R };
}

// Calculate overlap between two primitive Gaussian functions
static inline double calculatePrimitiveOverlap(double alpha1, double alpha2, double R2,
    int l1, int m1, int n1, int l2, int m2, int n2,
    double Ax, double Ay, double Az,
    double Bx, double By, double Bz)
{
    // Gaussian product theorem
    double gamma = alpha1 + alpha2;
    double prefactor = std::pow(M_PI / gamma, 1.5) * std::exp(-alpha1 * alpha2 * R2 / gamma);

    // Compute center of product Gaussian
    double Px = (alpha1 * Ax + alpha2 * Bx) / gamma;
    double Py = (alpha1 * Ay + alpha2 * By) / gamma;
    double Pz = (alpha1 * Az + alpha2 * Bz) / gamma;

    // Calculate overlap integrals along each axis
    auto overlap_1D = [](int l1, int l2, double PA, double PB, double gamma) {
        double sum = 0.0;

        for (int i = 0; i <= l1; i++) {
            for (int j = 0; j <= l2; j++) {
                if ((i + j) % 2 == 0) { // Only even terms contribute
                    double binom_l1_i = binomial(l1, i);
                    double binom_l2_j = binomial(l2, j);
                    double prefactor = binom_l1_i * binom_l2_j * factorial(i + j - 1) / std::pow(2 * gamma, (i + j) / 2.0);

                    sum += prefactor * std::pow(PA, l1 - i) * std::pow(PB, l2 - j);
                }
            }
        }

        return sum;
    };

    // Calculate overlap along each coordinate
    double overlap_x = overlap_1D(l1, l2, Px - Ax, Px - Bx, gamma);
    double overlap_y = overlap_1D(m1, m2, Py - Ay, Py - By, gamma);
    double overlap_z = overlap_1D(n1, n2, Pz - Az, Pz - Bz, gamma);

    return prefactor * overlap_x * overlap_y * overlap_z;
}

// Convert orbital type to angular momentum components
static inline void orbitalTypeToComponents(OrbitalType type, int& l, int& m, int& n)
{
    switch (type) {
    case S:
        l = 0;
        m = 0;
        n = 0;
        break;
    case PX:
        l = 1;
        m = 0;
        n = 0;
        break;
    case PY:
        l = 0;
        m = 1;
        n = 0;
        break;
    case PZ:
        l = 0;
        m = 0;
        n = 1;
        break;
    case DXY:
        l = 1;
        m = 1;
        n = 0;
        break;
    case DYZ:
        l = 0;
        m = 1;
        n = 1;
        break;
    case DZX:
        l = 1;
        m = 0;
        n = 1;
        break;
    case DX2Y2:
        l = 2;
        m = 0;
        n = 0;
        break; // This is an approximation - we need proper linear combinations
    case DZ2:
        l = 0;
        m = 0;
        n = 2;
        break; // This is an approximation - we need proper linear combinations
    default:
        l = 0;
        m = 0;
        n = 0;
        break;
    }
}

// Calculate overlap between two contracted Gaussian orbitals
static inline double calculateOverlap(const Orbital& orb1, const Orbital& orb2)
{
    // Calculate distance between centers
    double dx = orb2.x - orb1.x;
    double dy = orb2.y - orb1.y;
    double dz = orb2.z - orb1.z;
    double R2 = dx * dx + dy * dy + dz * dz;

    std::cout << "Calculating GTO overlap between orbitals at positions: "
              << "(" << orb1.x << "," << orb1.y << "," << orb1.z << ") and "
              << "(" << orb2.x << "," << orb2.y << "," << orb2.z << ")"
              << " with R² = " << R2 << std::endl;

    // Special case for identical orbitals at the same center
    if (R2 < 1e-10 && orb1.type == orb2.type && orb1.exponents.size() == orb2.exponents.size()) {
        bool identical = true;
        for (size_t i = 0; i < orb1.exponents.size(); i++) {
            if (std::abs(orb1.exponents[i] - orb2.exponents[i]) > 1e-10 || std::abs(orb1.coefficients[i] - orb2.coefficients[i]) > 1e-10) {
                identical = false;
                break;
            }
        }
        if (identical) {
            std::cout << "Identical orbitals detected - returning 1.0" << std::endl;
            return 1.0;
        }
    }

    // Convert orbital types to angular momentum components
    int l1, m1, n1, l2, m2, n2;
    orbitalTypeToComponents(orb1.type, l1, m1, n1);
    orbitalTypeToComponents(orb2.type, l2, m2, n2);

    // Special handling for d orbitals (linear combinations)
    // For simplicity we've approximated them above, but proper treatment involves
    // linear combinations of cartesian GTOs

    // Calculate overlap for each pair of primitive Gaussians
    double overlap = 0.0;

    for (size_t i = 0; i < orb1.exponents.size(); i++) {
        for (size_t j = 0; j < orb2.exponents.size(); j++) {
            double alpha1 = orb1.exponents[i];
            double alpha2 = orb2.exponents[j];
            double c1 = orb1.coefficients[i];
            double c2 = orb2.coefficients[j];

            // Normalize primitive Gaussians
            double norm1 = std::pow(2 * alpha1 / M_PI, 0.75) * std::sqrt(std::pow(4 * alpha1, l1 + m1 + n1) / (factorial(2 * l1 - 1) * factorial(2 * m1 - 1) * factorial(2 * n1 - 1)));

            double norm2 = std::pow(2 * alpha2 / M_PI, 0.75) * std::sqrt(std::pow(4 * alpha2, l2 + m2 + n2) / (factorial(2 * l2 - 1) * factorial(2 * m2 - 1) * factorial(2 * n2 - 1)));

            // Calculate primitive overlap
            double primitive_overlap = calculatePrimitiveOverlap(
                alpha1, alpha2, R2,
                l1, m1, n1, l2, m2, n2,
                orb1.x, orb1.y, orb1.z,
                orb2.x, orb2.y, orb2.z);

            // Add contribution from this pair of primitives
            overlap += c1 * c2 * norm1 * norm2 * primitive_overlap;

            std::cout << "  Primitive overlap: alpha1=" << alpha1 << ", alpha2=" << alpha2
                      << ", c1=" << c1 << ", c2=" << c2
                      << ", norm1=" << norm1 << ", norm2=" << norm2
                      << ", primitive_overlap=" << primitive_overlap
                      << ", contribution=" << c1 * c2 * norm1 * norm2 * primitive_overlap
                      << std::endl;
        }
    }

    std::cout << "Final GTO overlap result: " << overlap << std::endl
              << std::endl;
    return overlap;
}

// Create a contracted Gaussian orbital from standard basis set specifications
static inline Orbital createOrbital(double x, double y, double z, OrbitalType type,
    const std::vector<double>& exponents,
    const std::vector<double>& coefficients,
    double vsip, int atom_index)
{
    Orbital orb;
    orb.x = x;
    orb.y = y;
    orb.z = z;
    orb.type = type;
    orb.exponents = exponents;
    orb.coefficients = coefficients;
    orb.VSIP = vsip;
    orb.atom = atom_index;
    return orb;
}

// Create a standard STO-3G s-type orbital
static inline Orbital createSTO3G_S(double x, double y, double z, double zeta, double vsip, int atom_index)
{
    // STO-3G contraction coefficients and exponents for s-type orbitals
    std::vector<double> exponents(3);
    std::vector<double> coefficients(3);

    // These parameters are for STO-3G with zeta = 1.0
    // So we need to scale the exponents by zeta²
    const double alpha1 = 0.16885540;
    const double alpha2 = 0.62391373;
    const double alpha3 = 3.42525091;
    const double c1 = 0.44463454;
    const double c2 = 0.53532814;
    const double c3 = 0.15432897;

    exponents[0] = alpha1 * zeta * zeta;
    exponents[1] = alpha2 * zeta * zeta;
    exponents[2] = alpha3 * zeta * zeta;

    coefficients[0] = c1;
    coefficients[1] = c2;
    coefficients[2] = c3;

    return createOrbital(x, y, z, OrbitalType::S, exponents, coefficients, vsip, atom_index);
}

// Create a standard STO-3G p-type orbital
static inline Orbital createSTO3G_P(double x, double y, double z, OrbitalType p_type,
    double zeta, double vsip, int atom_index)
{
    // STO-3G contraction coefficients and exponents for p-type orbitals
    std::vector<double> exponents(3);
    std::vector<double> coefficients(3);

    // These parameters are for STO-3G with zeta = 1.0
    const double alpha1 = 0.0941539;
    const double alpha2 = 0.4057711;
    const double alpha3 = 2.2271139;
    const double c1 = 0.44463454;
    const double c2 = 0.53532814;
    const double c3 = 0.15432897;

    exponents[0] = alpha1 * zeta * zeta;
    exponents[1] = alpha2 * zeta * zeta;
    exponents[2] = alpha3 * zeta * zeta;

    coefficients[0] = c1;
    coefficients[1] = c2;
    coefficients[2] = c3;

    return createOrbital(x, y, z, p_type, exponents, coefficients, vsip, atom_index);
}

// Create standard basis sets for common elements
namespace STO3G {
    // Hydrogen 1s orbital
    static inline Orbital createH_1s(double x, double y, double z, int atom_index)
    {
        return createSTO3G_S(x, y, z, 1.24, -13.6, atom_index); // zeta value for Hydrogen
    }

    // Carbon 2s orbital
    static inline Orbital createC_2s(double x, double y, double z, int atom_index)
    {
        return createSTO3G_S(x, y, z, 1.6083, -19.44, atom_index); // zeta value for Carbon 2s
    }

    // Carbon 2p orbitals
    static inline Orbital createC_2px(double x, double y, double z, int atom_index)
    {
        return createSTO3G_P(x, y, z, OrbitalType::PX, 1.5679, -10.67, atom_index);
    }

    static inline Orbital createC_2py(double x, double y, double z, int atom_index)
    {
        return createSTO3G_P(x, y, z, OrbitalType::PY, 1.5679, -10.67, atom_index);
    }

    static inline Orbital createC_2pz(double x, double y, double z, int atom_index)
    {
        return createSTO3G_P(x, y, z, OrbitalType::PZ, 1.5679, -10.67, atom_index);
    }

    // Oxygen 2s orbital
    static inline Orbital createO_2s(double x, double y, double z, int atom_index)
    {
        return createSTO3G_S(x, y, z, 2.2458, -32.3, atom_index);
    }

    // Oxygen 2p orbitals
    static inline Orbital createO_2px(double x, double y, double z, int atom_index)
    {
        return createSTO3G_P(x, y, z, OrbitalType::PX, 2.2266, -14.8, atom_index);
    }

    static inline Orbital createO_2py(double x, double y, double z, int atom_index)
    {
        return createSTO3G_P(x, y, z, OrbitalType::PY, 2.2266, -14.8, atom_index);
    }

    static inline Orbital createO_2pz(double x, double y, double z, int atom_index)
    {
        return createSTO3G_P(x, y, z, OrbitalType::PZ, 2.2266, -14.8, atom_index);
    }

    // Nitrogen 2s orbital
    static inline Orbital createN_2s(double x, double y, double z, int atom_index)
    {
        return createSTO3G_S(x, y, z, 1.9170, -26.0, atom_index);
    }

    // Nitrogen 2p orbitals
    static inline Orbital createN_2px(double x, double y, double z, int atom_index)
    {
        return createSTO3G_P(x, y, z, OrbitalType::PX, 1.9170, -13.4, atom_index);
    }

    static inline Orbital createN_2py(double x, double y, double z, int atom_index)
    {
        return createSTO3G_P(x, y, z, OrbitalType::PY, 1.9170, -13.4, atom_index);
    }

    static inline Orbital createN_2pz(double x, double y, double z, int atom_index)
    {
        return createSTO3G_P(x, y, z, OrbitalType::PZ, 1.9170, -13.4, atom_index);
    }
} // namespace STO3G

} // namespace GTO
