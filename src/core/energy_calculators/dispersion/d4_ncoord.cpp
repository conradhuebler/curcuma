/*
 * D4 electronegativity-weighted covalent coordination number — implementation
 *
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 — exact port of dftd4's EN-weighted covalent CN.
 * The rad[] / en[] tables and the kn/k4/k5/k6 constants are copied verbatim
 * from external/cpp-d4/src/dftd_ncoord.cpp so the values are bit-identical to
 * the dftd4 reference that tblite's GFN2 D4 uses (no re-derivation, no reuse of
 * curcuma's Pyykkö COVALENT_RADII table).
 */

#include "d4_ncoord.h"

#include <cmath>

namespace curcuma::dispersion {

namespace {

// D3 covalent radii used to construct the coordination number, already scaled
// as rad = covalent_rad_d3 · 4/3 · Å→Bohr (Bohr), indexed by atomic number Z
// (index 0 is a dummy). Verbatim from external/cpp-d4/src/dftd_ncoord.cpp:73.
constexpr double rad[119] = {
    0.0,
    0.8062831465047213,  1.1590320231005369,  3.0235617993927044,  2.3684567428576182,
    1.9401188212769855,  1.8897261246204402,  1.7889407313073502,  1.5873699446811698,
    1.6125662930094427,  1.6881553379942602,  3.5274887659581555,  3.1495435410340673,
    2.8471873610947966,  2.6204202261403440,  2.7715983161099791,  2.5700275294837986,
    2.4944384844989811,  2.4188494395141635,  4.4345573057759662,  3.8802376425539711,
    3.3511143276602473,  3.0739544960492493,  3.0487581477209771,  2.7715983161099791,
    2.6960092711251620,  2.6204202261403440,  2.5196348328272538,  2.4944384844989811,
    2.5448311811555264,  2.7464019677817069,  2.8219910127665244,  2.7464019677817069,
    2.8975800577513415,  2.7715983161099791,  2.8723837094230693,  2.9479727544078864,
    4.7621098340435095,  4.2077901708215135,  3.7038632042560633,  3.5022924176298824,
    3.3259179793319751,  3.1243471927057946,  2.8975800577513415,  2.8471873610947966,
    2.8471873610947966,  2.7212056194534342,  2.8975800577513415,  3.0991508443775224,
    3.2251325860188853,  3.1747398893623395,  3.1747398893623395,  3.0991508443775224,
    3.3259179793319751,  3.3007216310037024,  5.2660368006089602,  4.4345573057759662,
    4.0818084291801515,  3.7038632042560633,  3.9810230358670613,  3.9558266875387886,
    3.9306303392105164,  3.9054339908822433,  3.8046485975691535,  3.8298449458974257,
    3.8046485975691535,  3.7794522492408804,  3.7542559009126082,  3.7542559009126082,
    3.7290595525843355,  3.8550412942256984,  3.6786668559277906,  3.4518997209733380,
    3.3007216310037024,  3.0991508443775224,  2.9731691027361595,  2.9227764060796142,
    2.7967946644382522,  2.8219910127665244,  2.8471873610947966,  3.3259179793319751,
    3.2755252826754302,  3.2755252826754302,  3.4267033726450653,  3.3007216310037024,
    3.4770960693016097,  3.5778814626147004,  5.0644660139827797,  4.5605390474173291,
    4.2077901708215135,  3.9810230358670613,  3.8298449458974257,  3.8550412942256984,
    3.8802376425539711,  3.9054339908822433,  3.7542559009126082,  3.7542559009126082,
    3.8046485975691535,  3.8046485975691535,  3.7290595525843355,  3.7794522492408804,
    3.9306303392105164,  3.9810230358670613,  3.6534705075995175,  3.5526851142864277,
    3.3763106759885204,  3.2503289343471575,  3.1999362376906122,  3.0487581477209771,
    2.9227764060796142,  2.8975800577513415,  2.7464019677817069,  3.0739544960492493,
    3.4267033726450653,  3.6030778109429726,  3.6786668559277906,  3.9810230358670613,
    3.7290595525843355,  3.9558266875387886};

// Pauling electronegativities, indexed by atomic number Z (index 0 dummy).
// Verbatim from external/cpp-d4/src/dftd_ncoord.cpp:195 (Fr–Og are dummies=1.5).
constexpr double en[119] = {
    0.0,
    2.20, 3.00,                                                 //  H,He
    0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 4.50,             //  Li-Ne
    0.93, 1.31, 1.61, 1.90, 2.19, 2.58, 3.16, 3.50,             //  Na-Ar
    0.82, 1.00,                                                 //  K,Ca
    1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.90, 1.65, //  Sc-Zn
    1.81, 2.01, 2.18, 2.55, 2.96, 3.00,                         //  Ga-Kr
    0.82, 0.95,                                                 //  Rb,Sr
    1.22, 1.33, 1.60, 2.16, 1.90, 2.20, 2.28, 2.20, 1.93, 1.69, //  Y-Cd
    1.78, 1.96, 2.05, 2.10, 2.66, 2.60,                         //  In-Xe
    0.79, 0.89,                                                 //  Cs,Ba
    1.10, 1.12, 1.13, 1.14, 1.15, 1.17, 1.18,                   //  La-Eu
    1.20, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26,                   //  Gd-Yb
    1.27, 1.30, 1.50, 2.36, 1.90, 2.20, 2.20, 2.28, 2.54, 2.00, //  Lu-Hg
    1.62, 2.33, 2.02, 2.00, 2.20, 2.20,                         //  Tl-Rn
    1.50, 1.50,                                                 //  Fr,Ra (dummy)
    1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,                   //  Ac-Am (dummy)
    1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,                   //  Cm-No (dummy)
    1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, //  Rf-Cn (dummy)
    1.50, 1.50, 1.50, 1.50, 1.50, 1.50};                        //  Nh-Og (dummy)

constexpr double kn = 7.5;
constexpr double k4 = 4.10451;
constexpr double k5 = 19.08857;
const double k6 = 2.0 * std::pow(11.28174, 2);
const double hlfosqrtpi = 1.0 / 1.7724538509055159;  // 1/sqrt(pi)

constexpr int MAX_Z = 118;

// erf counting function: 0.5·(1 + erf(−k·(rr−1)))
inline double erf_count(double k, double rr)
{
    return 0.5 * (1.0 + std::erf(-k * (rr - 1.0)));
}

// d(erf_count)/d(rr) = −k/sqrt(π)·exp(−(k·(rr−1))²)
inline double derf_count(double k, double rr)
{
    return -k * hlfosqrtpi * std::exp(-std::pow(k * (rr - 1.0), 2));
}

}  // namespace

std::vector<double> computeD4CovalentCN(const std::vector<int>& atoms,
                                        const Matrix& geom_bohr,
                                        double cutoff_bohr)
{
    const int N = static_cast<int>(atoms.size());
    std::vector<double> cn(N, 0.0);

    for (int i = 0; i < N; ++i) {
        const int zi = atoms[i];
        if (zi < 1 || zi > MAX_Z) continue;
        for (int j = 0; j < i; ++j) {
            const int zj = atoms[j];
            if (zj < 1 || zj > MAX_Z) continue;

            const double r = (geom_bohr.row(i) - geom_bohr.row(j)).norm();
            if (r > cutoff_bohr) continue;

            const double rcovij = rad[zi] + rad[zj];
            const double rr = r / rcovij;
            const double den = k4 * std::exp(-std::pow(std::fabs(en[zi] - en[zj]) + k5, 2) / k6);
            const double countf = den * erf_count(kn, rr);

            cn[i] += countf;
            cn[j] += countf;
        }
    }
    return cn;
}

void addD4CovalentCNGradient(const std::vector<int>& atoms,
                             const Matrix& geom_bohr,
                             const Vector& dEdcn,
                             Matrix& grad_out,
                             double cutoff_bohr)
{
    const int N = static_cast<int>(atoms.size());
    if (dEdcn.size() != N) return;

    for (int i = 0; i < N; ++i) {
        const int zi = atoms[i];
        if (zi < 1 || zi > MAX_Z) continue;
        for (int j = 0; j < i; ++j) {
            const int zj = atoms[j];
            if (zj < 1 || zj > MAX_Z) continue;

            const Eigen::Vector3d rij = geom_bohr.row(i) - geom_bohr.row(j);
            const double r = rij.norm();
            if (r > cutoff_bohr || r < 1e-12) continue;

            const double rcovij = rad[zi] + rad[zj];
            const double rr = r / rcovij;
            const double den = k4 * std::exp(-std::pow(std::fabs(en[zi] - en[zj]) + k5, 2) / k6);
            // d(countf_ij)/dr  (rr = r/rcovij ⇒ d(rr)/dr = 1/rcovij)
            const double dcountf = den * derf_count(kn, rr) / rcovij;

            const double factor = (dEdcn(i) + dEdcn(j)) * dcountf / r;
            grad_out.row(i) += factor * rij.transpose();
            grad_out.row(j) -= factor * rij.transpose();
        }
    }
}

}  // namespace curcuma::dispersion
