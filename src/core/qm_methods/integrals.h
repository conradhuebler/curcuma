#pragma once

#include <cmath>
#include <iostream>

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

/* Integrals for Slater Type Orbitals */
namespace STO {

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

struct Orbital {
    double x = 0, y = 0, z = 0; // Position
    double zeta; // Slater exponent
    double VSIP; // Valence State Ionization Potential
    OrbitalType type; // Orbital type
    int atom;
};

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

static inline double calculateSSOverlap(double zeta1, double zeta2, double R)
{
    // For identical centers (R=0)
    if (R < 1e-10)
        return 1.0;

    // Korrigierte Formel für S-S Überlappung
    double rho = R * (zeta1 + zeta2) / 2.0;

    // Explizit berechneter Normierungsfaktor
    double numerator = 2.0 * std::sqrt(zeta1 * zeta2);
    double denominator = std::pow(zeta1 + zeta2, 1.5);
    double N = numerator / denominator;

    double poly_term = 1.0 + rho + rho * rho / 3.0;
    double exp_term = std::exp(-rho);

    double result = N * exp_term * poly_term;

    // Debug-Ausgabe
    std::cout << "S-S overlap: zeta1=" << zeta1 << ", zeta2=" << zeta2 << ", R=" << R
              << ", rho=" << rho
              << ", numerator=" << numerator
              << ", denominator=" << denominator
              << ", N=" << N
              << ", exp_term=" << exp_term
              << ", poly_term=" << poly_term
              << ", result=" << result << std::endl;

    return result;
}

static inline double calculateSPOverlap(double zeta_s, double zeta_p, double R, double l)
{
    // For identical centers (R=0)
    if (R < 1e-10)
        return 0.0;

    // Korrigierte Formel für S-P Überlappung
    double rho = R * (zeta_s + zeta_p) / 2.0;

    // Explizit berechneter Normierungsfaktor
    double numerator = 2.0 * std::sqrt(zeta_s * zeta_p);
    double denominator = std::pow(zeta_s + zeta_p, 1.5);
    double N = numerator / denominator;

    double poly_term = rho * (1.0 + rho / 3.0);
    double exp_term = std::exp(-rho);

    double result = N * l * exp_term * poly_term;

    // Debug-Ausgabe
    std::cout << "S-P overlap: zeta_s=" << zeta_s << ", zeta_p=" << zeta_p
              << ", R=" << R << ", l=" << l << ", rho=" << rho
              << ", N=" << N
              << ", exp_term=" << exp_term
              << ", poly_term=" << poly_term
              << ", result=" << result << std::endl;

    return result;
}

static inline double calculatePPOverlap(double zeta1, double zeta2, double R, double cos_angle, bool same_axis)
{
    // For identical centers (R=0)
    if (R < 1e-10)
        return same_axis ? 1.0 : 0.0;

    // Korrigierte Formel für P-P Überlappung
    double rho = R * (zeta1 + zeta2) / 2.0;

    // Explizit berechneter Normierungsfaktor
    double numerator = 2.0 * std::sqrt(zeta1 * zeta2);
    double denominator = std::pow(zeta1 + zeta2, 1.5);
    double N = numerator / denominator;

    double exp_term = std::exp(-rho);
    double result = 0.0;

    if (same_axis) {
        double poly_term = (cos_angle * cos_angle) * (1.0 + rho + (2.0 * rho * rho) / 5.0) - (rho * rho) / 5.0;
        result = N * exp_term * poly_term;
    } else {
        double poly_term = cos_angle * cos_angle * (1.0 + rho + (2.0 * rho * rho) / 5.0);
        result = N * exp_term * poly_term;
    }

    // Debug-Ausgabe
    std::cout << "P-P overlap: zeta1=" << zeta1 << ", zeta2=" << zeta2
              << ", R=" << R << ", cos_angle=" << cos_angle
              << ", same_axis=" << (same_axis ? "true" : "false")
              << ", rho=" << rho << ", N=" << N
              << ", exp_term=" << exp_term
              << ", result=" << result << std::endl;

    return result;
}

static inline double calculateSDOverlap(double zeta_s, double zeta_d, double R,
    const DirectionCosines& dc, OrbitalType d_type)
{
    // For identical centers (R=0)
    if (R < 1e-10)
        return 0.0;

    // Korrigierte Formel für S-D Überlappung
    double rho = R * (zeta_s + zeta_d) / 2.0;

    // Explizit berechneter Normierungsfaktor
    double numerator = 2.0 * std::sqrt(zeta_s * zeta_d);
    double denominator = std::pow(zeta_s + zeta_d, 1.5);
    double N = numerator / denominator;

    double exp_term = std::exp(-rho);
    double poly_term = rho * rho * (1.0 + rho / 5.0);

    double angular_term = 0.0;
    std::string d_type_str;

    switch (d_type) {
    case DXY:
        angular_term = std::sqrt(3.0) * dc.l * dc.m;
        d_type_str = "DXY";
        break;
    case DYZ:
        angular_term = std::sqrt(3.0) * dc.m * dc.n;
        d_type_str = "DYZ";
        break;
    case DZX:
        angular_term = std::sqrt(3.0) * dc.n * dc.l;
        d_type_str = "DZX";
        break;
    case DX2Y2:
        angular_term = (std::sqrt(3.0) / 2.0) * (dc.l * dc.l - dc.m * dc.m);
        d_type_str = "DX2Y2";
        break;
    case DZ2:
        angular_term = (3.0 * dc.n * dc.n - 1.0) / 2.0;
        d_type_str = "DZ2";
        break;
    default:
        d_type_str = "unknown";
    }

    double result = N * exp_term * poly_term * angular_term;

    // Debug-Ausgabe
    std::cout << "S-D overlap: zeta_s=" << zeta_s << ", zeta_d=" << zeta_d
              << ", R=" << R << ", d_type=" << d_type_str
              << ", dc=(" << dc.l << "," << dc.m << "," << dc.n << ")"
              << ", rho=" << rho << ", N=" << N
              << ", exp_term=" << exp_term
              << ", poly_term=" << poly_term
              << ", angular_term=" << angular_term
              << ", result=" << result << std::endl;

    return result;
}

static inline double calculatePDOverlap(double zeta_p, double zeta_d, double R,
    const DirectionCosines& dc,
    OrbitalType p_type, OrbitalType d_type)
{
    // For identical centers (R=0)
    if (R < 1e-10)
        return 0.0;

    // Korrigierte Formel für P-D Überlappung
    double rho = R * (zeta_p + zeta_d) / 2.0;

    // Explizit berechneter Normierungsfaktor
    double numerator = 2.0 * std::sqrt(zeta_p * zeta_d);
    double denominator = std::pow(zeta_p + zeta_d, 1.5);
    double N = numerator / denominator;

    double exp_term = std::exp(-rho);
    double poly_term = rho * (1.0 + rho / 3.0 + (rho * rho) / 15.0);

    // Mapping von P-Orbital-Richtungen
    double p_dir = 0.0;
    std::string p_type_str, d_type_str;

    if (p_type == PX) {
        p_dir = dc.l;
        p_type_str = "PX";
    } else if (p_type == PY) {
        p_dir = dc.m;
        p_type_str = "PY";
    } else if (p_type == PZ) {
        p_dir = dc.n;
        p_type_str = "PZ";
    }

    double angular_term = 0.0;

    switch (d_type) {
    case DXY:
        d_type_str = "DXY";
        if (p_type == PX)
            angular_term = std::sqrt(3.0) * dc.m;
        else if (p_type == PY)
            angular_term = std::sqrt(3.0) * dc.l;
        else
            angular_term = 0.0; // PZ with DXY is zero
        break;

    case DYZ:
        d_type_str = "DYZ";
        if (p_type == PX)
            angular_term = 0.0; // PX with DYZ is zero
        else if (p_type == PY)
            angular_term = std::sqrt(3.0) * dc.n;
        else if (p_type == PZ)
            angular_term = std::sqrt(3.0) * dc.m;
        break;

    case DZX:
        d_type_str = "DZX";
        if (p_type == PX)
            angular_term = std::sqrt(3.0) * dc.n;
        else if (p_type == PY)
            angular_term = 0.0; // PY with DZX is zero
        else if (p_type == PZ)
            angular_term = std::sqrt(3.0) * dc.l;
        break;

    case DX2Y2:
        d_type_str = "DX2Y2";
        if (p_type == PX)
            angular_term = std::sqrt(3.0) * dc.l;
        else if (p_type == PY)
            angular_term = -std::sqrt(3.0) * dc.m;
        else
            angular_term = 0.0; // PZ with DX2Y2 is zero
        break;

    case DZ2:
        d_type_str = "DZ2";
        if (p_type == PX)
            angular_term = -(std::sqrt(3.0) / 2.0) * dc.l;
        else if (p_type == PY)
            angular_term = -(std::sqrt(3.0) / 2.0) * dc.m;
        else if (p_type == PZ)
            angular_term = std::sqrt(3.0) * dc.n;
        break;

    default:
        d_type_str = "unknown";
    }

    double result = N * exp_term * poly_term * angular_term;

    // Debug-Ausgabe
    std::cout << "P-D overlap: p_type=" << p_type_str << ", d_type=" << d_type_str
              << ", zeta_p=" << zeta_p << ", zeta_d=" << zeta_d
              << ", R=" << R << ", dc=(" << dc.l << "," << dc.m << "," << dc.n << ")"
              << ", rho=" << rho << ", N=" << N
              << ", exp_term=" << exp_term
              << ", poly_term=" << poly_term
              << ", angular_term=" << angular_term
              << ", result=" << result << std::endl;

    return result;
}

static inline double calculateDDOverlap(double zeta1, double zeta2, double R,
    const DirectionCosines& dc,
    OrbitalType d_type1, OrbitalType d_type2)
{
    // For identical centers (R=0)
    if (R < 1e-10)
        return (d_type1 == d_type2) ? 1.0 : 0.0;

    // Korrigierte Formel für D-D Überlappung
    double rho = R * (zeta1 + zeta2) / 2.0;

    // Explizit berechneter Normierungsfaktor
    double numerator = 2.0 * std::sqrt(zeta1 * zeta2);
    double denominator = std::pow(zeta1 + zeta2, 1.5);
    double N = numerator / denominator;

    double exp_term = std::exp(-rho);
    double poly_term = (1.0 + rho + (rho * rho) / 5.0 + (rho * rho * rho) / 35.0);

    double angular_term = 0.0;
    std::string d_type1_str, d_type2_str;

    // Konvertiere die Enum-Werte in Strings für Debug-Ausgabe
    switch (d_type1) {
    case DXY:
        d_type1_str = "DXY";
        break;
    case DYZ:
        d_type1_str = "DYZ";
        break;
    case DZX:
        d_type1_str = "DZX";
        break;
    case DX2Y2:
        d_type1_str = "DX2Y2";
        break;
    case DZ2:
        d_type1_str = "DZ2";
        break;
    default:
        d_type1_str = "unknown";
    }

    switch (d_type2) {
    case DXY:
        d_type2_str = "DXY";
        break;
    case DYZ:
        d_type2_str = "DYZ";
        break;
    case DZX:
        d_type2_str = "DZX";
        break;
    case DX2Y2:
        d_type2_str = "DX2Y2";
        break;
    case DZ2:
        d_type2_str = "DZ2";
        break;
    default:
        d_type2_str = "unknown";
    }

    // Die D-D Überlappungsformeln sind komplex und fallspezifisch
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
            angular_term = (3.0 * pow(dc.l * dc.l - dc.m * dc.m, 2) / 4.0 - (rho * rho) / 35.0);
            break;
        case DZ2:
            angular_term = (pow(3.0 * dc.n * dc.n - 1.0, 2) / 4.0 - (rho * rho) / 35.0);
            break;
        }
    } else if ((d_type1 == DXY && d_type2 == DX2Y2) || (d_type1 == DX2Y2 && d_type2 == DXY)) {
        angular_term = (3.0 * std::sqrt(3.0) / 4.0) * dc.l * dc.m * (dc.l * dc.l - dc.m * dc.m);
    } else if ((d_type1 == DXY && d_type2 == DZ2) || (d_type1 == DZ2 && d_type2 == DXY)) {
        angular_term = (3.0 * std::sqrt(3.0) / 4.0) * dc.l * dc.m * (3.0 * dc.n * dc.n - 1.0);
    }
    // Weitere D-D Überlappungsfälle würden hier implementiert

    double result = N * exp_term * poly_term * angular_term;

    // Debug-Ausgabe
    std::cout << "D-D overlap: d_type1=" << d_type1_str << ", d_type2=" << d_type2_str
              << ", zeta1=" << zeta1 << ", zeta2=" << zeta2
              << ", R=" << R << ", dc=(" << dc.l << "," << dc.m << "," << dc.n << ")"
              << ", rho=" << rho << ", N=" << N
              << ", exp_term=" << exp_term
              << ", poly_term=" << poly_term
              << ", angular_term=" << angular_term
              << ", result=" << result << std::endl;

    return result;
}

static inline double calculateOverlap(const Orbital& orb1, const Orbital& orb2)
{
    double dx = orb2.x - orb1.x;
    double dy = orb2.y - orb1.y;
    double dz = orb2.z - orb1.z;
    double R = std::sqrt(dx * dx + dy * dy + dz * dz);

    std::cout << "Calculating overlap between orbitals at positions: "
              << "(" << orb1.x << "," << orb1.y << "," << orb1.z << ") and "
              << "(" << orb2.x << "," << orb2.y << "," << orb2.z << ")"
              << " with distance R=" << R << std::endl;

    // Get direction cosines
    DirectionCosines dc = getDirectionCosines(orb1, orb2);

    std::cout << "Direction cosines: (" << dc.l << ", " << dc.m << ", " << dc.n << ")" << std::endl;

    // Special case for identical orbitals
    if (R < 1e-10 && orb1.type == orb2.type && std::abs(orb1.zeta - orb2.zeta) < 1e-10) {
        std::cout << "Identical orbitals detected - returning 1.0" << std::endl;
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

        std::cout << "Swapped orbital order for calculation" << std::endl;
    }

    double result = 0.0;

    // Calculate appropriate overlap based on orbital types
    if (type1 == S) {
        if (type2 == S) {
            result = calculateSSOverlap(zeta1, zeta2, R);
        } else if (type2 <= PZ) {
            double direction = 0.0;
            if (type2 == PX)
                direction = dc.l;
            else if (type2 == PY)
                direction = dc.m;
            else
                direction = dc.n;

            result = calculateSPOverlap(zeta1, zeta2, R, direction);
        } else {
            result = calculateSDOverlap(zeta1, zeta2, R, dc, type2);
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

            result = calculatePPOverlap(zeta1, zeta2, R, l1 * l2, same_axis);
        } else {
            result = calculatePDOverlap(zeta1, zeta2, R, dc, type1, type2);
        }
    } else {
        result = calculateDDOverlap(zeta1, zeta2, R, dc, type1, type2);
    }

    std::cout << "Final overlap result: " << result << std::endl
              << std::endl;
    return result;
}

} // namespace STO
