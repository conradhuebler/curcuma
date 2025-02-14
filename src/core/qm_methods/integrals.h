

#pragma once

#include <cmath>
static inline double factorial(int n)
{
    if (n <= 1)
        return 1;
    return n * factorial(n - 1);
}

static inline double binomial(int n, int k)
{
    return factorial(n) / (factorial(k) * factorial(n - k));
}

// Richtungskosinusse
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

    if (std::abs(dx) < 1e-8)
        return { 1.0, 0.0, 0.0 };
    else if (std::abs(dy) < 1e-8)
        return { 0.0, 1.0, 0.0 };
    else if (std::abs(dz) < 1e-8)
        return { 0.0, 0.0, 1.0 };
    else
        return { dx / R, dy / R, dz / R };
}

// Slater-Knotenfunktionen
static inline double SlaterRadial(int n, double zeta, double r)
{
    double norm = std::pow(2 * zeta, n + 0.5) / std::sqrt(factorial(2 * n));
    return norm * std::pow(r, n - 1) * std::exp(-zeta * r);
}

static inline double calculateSSOverlap(double zeta1, double zeta2, double R)
{
    std::cout << "zeta1: " << zeta1 << " zeta2: " << zeta2 << " R: " << R << std::endl;
    /*
    double zeta_sum = zeta1 + zeta2;
    double prefactor = 8.0 * pow(zeta1 * zeta2, 1.5) / pow(zeta_sum, 3);
    double exponent = exp(-zeta_sum * R / 2.0);
    double polynomial = 1.0 + (zeta_sum * R)/2.0 + (pow(zeta_sum, 2) * pow(R, 2))/12.0;
    return prefactor * exponent * polynomial;
    */
    double N = std::pow(zeta1 * zeta2, 0.75);
    return N * std::exp(-0.5 * (zeta1 + zeta2) * R);
}

static inline double calculateSPOverlap(double zeta_s, double zeta_p, double R, double l)
{
    /*
        double zeta_sum = zeta_s + zeta_p;
        double prefactor = sqrt(pow(zeta_s, 3) * pow(zeta_p, 5)) * 8.0 / pow(zeta_sum, 4);
        double exponent = exp(-zeta_sum * R / 2.0);
        double polynomial = (zeta_sum * R)/2.0 * (1.0 + (zeta_sum * R)/6.0);
        return l*prefactor * exponent * polynomial;
        */

    double N = std::pow(zeta_s * zeta_p, 0.75);
    double rho = (zeta_s + zeta_p) * R / 2.0;
    return N * l * R * std::exp(-rho) * (1.0 + rho);
}

static inline double calculatePPOverlap(double zeta1, double zeta2, double R, double l, bool same_axis)
{
    double N = std::pow(zeta1 * zeta2, 0.75);
    double rho = (zeta1 + zeta2) * R / 2.0;
    if (same_axis) {
        if (R < 1e-10)
            return N;
        std::cout << N << "* " << std::exp(-rho) << "* (" << l * l << "* )" << 1.0 << "+ " << rho << " +" << rho * rho / 3.0 << "  - "
                  << (rho * rho / 3.0) << "=" << N * std::exp(-rho) * (l * l * (1.0 + rho + rho * rho / 3.0) - (rho * rho / 3.0)) << std::endl;
        return N * std::exp(-rho) * (l * l * (1.0 + rho + rho * rho / 3.0) - (rho * rho / 3.0));
    } else {
        return N * l * l * std::exp(-rho) * (1.0 + rho + rho * rho / 3.0);
    }
}

static inline double calculateSDOverlap(double zeta1, double zeta2, double R,
    const DirectionCosines& dc, OrbitalType d_type)
{
    double N = std::pow(zeta1 * zeta2, 0.75);
    double rho = (zeta1 + zeta2) * R / 2.0;

    switch (d_type) {
    case DXY:
        return N * std::sqrt(3.0) * dc.l * dc.m * R * R * std::exp(-rho) * (1.0 + rho);
    case DYZ:
        return N * std::sqrt(3.0) * dc.m * dc.n * R * R * std::exp(-rho) * (1.0 + rho);
    case DZX:
        return N * std::sqrt(3.0) * dc.n * dc.l * R * R * std::exp(-rho) * (1.0 + rho);
    case DX2Y2:
        return N * 0.5 * std::sqrt(3.0) * (dc.l * dc.l - dc.m * dc.m) * R * R * std::exp(-rho) * (1.0 + rho);
    case DZ2:
        return N * (3.0 * dc.n * dc.n - 1.0) * R * R * std::exp(-rho) * (1.0 + rho);
    default:
        return 0.0;
    }
}

static inline double calculatePDOverlap(double zeta1, double zeta2, double R,
    const DirectionCosines& dc,
    OrbitalType p_type, OrbitalType d_type)
{
    double N = std::pow(zeta1 * zeta2, 0.75);
    double rho = (zeta1 + zeta2) * R / 2.0;
    double exp_term = std::exp(-rho);

    // Beispiel für px mit verschiedenen d-Orbitalen
    if (p_type == PX) {
        switch (d_type) {
        case DXY:
            return N * std::sqrt(3.0) * dc.m * R * exp_term * (1.0 + rho + rho * rho / 3.0);
        case DYZ:
            return 0.0; // Symmetriebedingt
        case DZX:
            return N * std::sqrt(3.0) * dc.n * R * exp_term * (1.0 + rho + rho * rho / 3.0);
        case DX2Y2:
            return N * std::sqrt(3.0) * dc.l * R * exp_term * (1.0 + rho + rho * rho / 3.0);
        case DZ2:
            return N * std::sqrt(3.0) * dc.l * R * exp_term * (1.0 + rho + rho * rho / 3.0);
        }
    } else if (p_type == PY) {
        switch (d_type) {
        case DXY:
            return N * std::sqrt(3.0) * dc.l * R * exp_term * (1.0 + rho + rho * rho / 3.0);
        case DYZ:
            return N * std::sqrt(3.0) * dc.n * R * exp_term * (1.0 + rho + rho * rho / 3.0);
        case DZX:
            return 0.0; // Symmetriebedingt
        case DX2Y2:
            return -N * std::sqrt(3.0) * dc.m * R * exp_term * (1.0 + rho + rho * rho / 3.0);
        case DZ2:
            return -N * std::sqrt(3.0) * dc.m * R * exp_term * (1.0 + rho + rho * rho / 3.0);
        }
    } else if (p_type == PZ) {
        switch (d_type) {
        case DXY:
            return 0.0; // Symmetriebedingt
        case DYZ:
            return N * std::sqrt(3.0) * dc.m * R * exp_term * (1.0 + rho + rho * rho / 3.0);
        case DZX:
            return N * std::sqrt(3.0) * dc.l * R * exp_term * (1.0 + rho + rho * rho / 3.0);
        case DX2Y2:
            return 0.0; // Symmetriebedingt
        case DZ2:
            return N * 2.0 * dc.n * R * exp_term * (1.0 + rho + rho * rho / 3.0);
        }
    }
    return 0.0;
}

static inline double calculateDDOverlap(double zeta1, double zeta2, double R,
    const DirectionCosines& dc,
    OrbitalType d_type1, OrbitalType d_type2)
{
    double N = std::pow(zeta1 * zeta2, 0.75);
    double rho = (zeta1 + zeta2) * R / 2.0;
    double exp_term = std::exp(-rho);

    // Matrix-Elemente für d-d Überlappung
    if (d_type1 == d_type2) {
        switch (d_type1) {
        case DXY:
        case DYZ:
        case DZX:
            return N * exp_term * (3.0 * std::pow(dc.l * dc.m, 2) * (1.0 + rho + rho * rho / 3.0 + rho * rho * rho / 15.0) - (rho * rho / 3.0 + rho * rho * rho / 15.0));
        case DX2Y2:
            return N * exp_term * (0.75 * std::pow(dc.l * dc.l - dc.m * dc.m, 2) * (1.0 + rho + rho * rho / 3.0 + rho * rho * rho / 15.0));
        case DZ2:
            return N * exp_term * (std::pow(3.0 * dc.n * dc.n - 1.0, 2) * (1.0 + rho + rho * rho / 3.0 + rho * rho * rho / 15.0));
        }
    }

    // D-D Kreuzterme
    // DXY mit anderen
    if (d_type1 == DXY && d_type2 == DYZ)
        return N * 3.0 * dc.l * dc.m * dc.m * dc.n * exp_term * (1.0 + rho + rho * rho / 3.0 + rho * rho * rho / 15.0);

    if (d_type1 == DXY && d_type2 == DZX)
        return N * 3.0 * dc.l * dc.l * dc.m * dc.n * exp_term * (1.0 + rho + rho * rho / 3.0 + rho * rho * rho / 15.0);

    if (d_type1 == DXY && d_type2 == DX2Y2)
        return N * 1.5 * dc.l * dc.m * (dc.l * dc.l - dc.m * dc.m) * exp_term * (1.0 + rho + rho * rho / 3.0 + rho * rho * rho / 15.0);

    if (d_type1 == DXY && d_type2 == DZ2)
        return N * std::sqrt(3.0) * dc.l * dc.m * (3.0 * dc.n * dc.n - 1.0) * exp_term * (1.0 + rho + rho * rho / 3.0 + rho * rho * rho / 15.0);

    // DYZ mit anderen
    if (d_type1 == DYZ && d_type2 == DZX)
        return N * 3.0 * dc.l * dc.m * dc.n * dc.n * exp_term * (1.0 + rho + rho * rho / 3.0 + rho * rho * rho / 15.0);

    if (d_type1 == DYZ && d_type2 == DX2Y2)
        return N * 1.5 * dc.m * dc.n * (dc.l * dc.l - dc.m * dc.m) * exp_term * (1.0 + rho + rho * rho / 3.0 + rho * rho * rho / 15.0);

    if (d_type1 == DYZ && d_type2 == DZ2)
        return N * std::sqrt(3.0) * dc.m * dc.n * (3.0 * dc.n * dc.n - 1.0) * exp_term * (1.0 + rho + rho * rho / 3.0 + rho * rho * rho / 15.0);

    // DZX mit anderen
    if (d_type1 == DZX && d_type2 == DX2Y2)
        return N * 1.5 * dc.n * dc.l * (dc.l * dc.l - dc.m * dc.m) * exp_term * (1.0 + rho + rho * rho / 3.0 + rho * rho * rho / 15.0);

    if (d_type1 == DZX && d_type2 == DZ2)
        return N * std::sqrt(3.0) * dc.n * dc.l * (3.0 * dc.n * dc.n - 1.0) * exp_term * (1.0 + rho + rho * rho / 3.0 + rho * rho * rho / 15.0);

    // DX2Y2 mit DZ2
    if (d_type1 == DX2Y2 && d_type2 == DZ2)
        return N * 0.5 * std::sqrt(3.0) * (dc.l * dc.l - dc.m * dc.m) * (3.0 * dc.n * dc.n - 1.0) * exp_term * (1.0 + rho + rho * rho / 3.0 + rho * rho * rho / 15.0);

    // Symmetrische Matrix: wenn d_type1 > d_type2, rekursiver Aufruf mit vertauschten Argumenten
    if (d_type1 > d_type2)
        return calculateDDOverlap(zeta2, zeta1, R, dc, d_type2, d_type1);

    return 0.0;
}

static inline double calculateOverlap(const Orbital& orb1, const Orbital& orb2)
{
    // Wenn gleiche Orbitale am gleichen Ort -> Überlappung = 1
    /* if (orb1.type == orb2.type &&
        std::abs(orb1.x - orb2.x) < 1e-10 &&
        std::abs(orb1.y - orb2.y) < 1e-10 &&
        std::abs(orb1.z - orb2.z) < 1e-10) {
        return 1.0;
    }*/

    DirectionCosines dc = getDirectionCosines(orb1, orb2);
    double R = std::sqrt(std::pow(orb2.x - orb1.x, 2) + std::pow(orb2.y - orb2.y, 2) + std::pow(orb2.z - orb2.z, 2));
    std::cout << "R: " << R << " ";
    // if(R < 1e-10)
    //     return 1.0;
    //  Sortierung der Orbitale nach Typ für eindeutige Behandlung
    OrbitalType type1 = orb1.type;
    OrbitalType type2 = orb2.type;
    double zeta1 = orb1.zeta;
    double zeta2 = orb2.zeta;

    // Symmetrische Behandlung: Immer niedrigeren Typ zuerst
    if (type1 > type2) {
        std::swap(type1, type2);
        std::swap(zeta1, zeta2);
    }

    // Auswahl des richtigen Überlappungsintegrals basierend auf Orbital-Typen
    if (type1 == S) {
        std::cout << " S - ";
        if (type2 == S) {
            std::cout << "S ";
            // s-s Überlappung
            return calculateSSOverlap(zeta1, zeta2, R);
        } else if (type2 <= PZ) {
            std::cout << "P ";
            // s-p Überlappung
            double direction = 0.0;
            if (type2 == PX)
                direction = dc.l;
            else if (type2 == PY)
                direction = dc.m;
            else
                direction = dc.n;
            return calculateSPOverlap(zeta1, zeta2, R, direction);
        } else {
            std::cout << "D ";
            // s-d Überlappung
            return calculateSDOverlap(zeta1, zeta2, R, dc, type2);
        }
    } else if (type1 <= PZ) {
        std::cout << " P - ";
        if (type2 <= PZ) {
            std::cout << "P ";
            // p-p Überlappung
            bool same_axis = (type1 == type2);
            double l1, l2;
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
            return calculatePPOverlap(zeta1, zeta2, R, l1 * l2, same_axis);
        } else {
            std::cout << "D ";
            // p-d Überlappung
            return calculatePDOverlap(zeta1, zeta2, R, dc, type1, type2);
        }
    } else {
        std::cout << " D - D ";
        // d-d Überlappung
        return calculateDDOverlap(zeta1, zeta2, R, dc, type1, type2);
    }
}
} // namespace STO