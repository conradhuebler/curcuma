/*
 * STO-to-CGTO Conversion and CGTO Overlap Integrals
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Implements Stewart's STO-NG expansion for converting Slater-type orbitals
 * to contracted Gaussian-type orbitals, and computes CGTO overlap integrals.
 *
 * Reference: R. F. Stewart, J. Chem. Phys. 52, 431-438 (1970)
 * TBLite implementation: src/tblite/basis/slater.f90
 *
 * Claude Generated (March 2026): Direct port from TBLite for GFN1/GFN2 overlap accuracy
 */

#pragma once

#include <array>
#include <cmath>
#include <vector>

namespace CGTO {

// Double factorial: (2l-1)!! for angular momentum normalization
constexpr double dfactorial[] = {1.0, 1.0, 3.0, 15.0, 105.0, 945.0, 10395.0, 135135.0};
constexpr double two_over_pi = 2.0 / M_PI;

// A single contracted GTO shell
struct Shell {
    int ang;                          // Angular momentum (0=s, 1=p, 2=d)
    int nprim;                        // Number of primitives
    std::vector<double> alpha;        // GTO exponents (already scaled by zeta^2)
    std::vector<double> coeff;        // Contraction coefficients (normalized)
};

// ============================================================================
// Stewart STO-NG Coefficients (from TBLite basis/slater.f90)
// Tables indexed by [primitive][ityp], where ityp = nlm_to_ityp(n, l)
// ============================================================================

// ityp mapping: 1s=0, 2s=1, 3s=2, 4s=3, 5s=4, 2p=5, 3p=6, 4p=7, 5p=8,
//               3d=9, 4d=10, 5d=11, 4f=12, 5f=13, 5g=14
inline int nlm_to_ityp(int n, int l)
{
    switch (l) {
    case 0: return n - 1;      // s: 1s=0, 2s=1, 3s=2, ...
    case 1: return 3 + n;      // p: 2p=5, 3p=6, ...
    case 2: return 6 + n;      // d: 3d=9, 4d=10, ...
    case 3: return 8 + n;      // f: 4f=12, 5f=13
    case 4: return 9 + n;      // g: 5g=14
    default: return -1;
    }
}

// --- STO-3G tables ---
constexpr double pAlpha3[3][15] = {
    {2.227660584e+0, 2.581578398e+0, 5.641487709e-1, 2.267938753e-1, 1.080198458e-1, 9.192379002e-1, 2.692880368e+0, 4.859692220e-1, 2.127482317e-1, 5.229112225e-1, 1.777717219e-1, 4.913352950e-1, 3.483826963e-1, 1.649233885e-1, 2.545432122e-1},
    {4.057711562e-1, 1.567622104e-1, 6.924421391e-2, 4.448178019e-2, 4.408119382e-2, 2.359194503e-1, 1.489359592e-1, 7.430216918e-2, 4.729648620e-2, 1.639595876e-1, 8.040647350e-2, 7.329090601e-2, 1.249380537e-1, 7.487066646e-2, 1.006544376e-1},
    {1.098175104e-1, 6.018332272e-2, 3.269529097e-2, 2.195294664e-2, 2.610811810e-2, 8.009805746e-2, 5.739585040e-2, 3.653340923e-2, 2.604865324e-2, 6.386630021e-2, 3.949855551e-2, 3.594209290e-2, 5.349995725e-2, 3.735787219e-2, 4.624463922e-2}
};
constexpr double pCoeff3[3][15] = {
    {1.543289673e-1, -5.994474934e-2, -1.782577972e-1, -3.349048323e-1, -6.617401158e-1, 1.623948553e-1, -1.061945788e-2, -6.147823411e-2, -1.389529695e-1, 1.686596060e-1, 2.308552718e-1, -2.010175008e-2, 1.737856685e-1, 1.909729355e-1, 1.780980905e-1},
    {5.353281423e-1, 5.960385398e-1, 8.612761663e-1, 1.056744667e+0, 7.467595004e-1, 5.661708862e-1, 5.218564264e-1, 6.604172234e-1, 8.076691064e-1, 5.847984817e-1, 6.042409177e-1, 5.899370608e-1, 5.973380628e-1, 6.146060459e-1, 6.063757846e-1},
    {4.446345422e-1, 4.581786291e-1, 2.261841969e-1, 1.256661680e-1, 7.146490945e-1, 4.223071752e-1, 5.450015143e-1, 3.932639495e-1, 2.726029342e-1, 4.056779523e-1, 2.595768926e-1, 4.658445960e-1, 3.929395614e-1, 3.059611271e-1, 3.828552923e-1}
};

// --- STO-4G tables ---
constexpr double pAlpha4[4][15] = {
    {5.216844534e+0, 1.161525551e+1, 1.513265591e+0, 3.242212833e-1, 8.602284252e-1, 1.798260992e+0, 1.853180239e+0, 1.492607880e+0, 3.962838833e-1, 9.185846715e-1, 1.995825422e+0, 4.230617826e-1, 5.691670217e-1, 2.017831152e-1, 3.945205573e-1},
    {9.546182760e-1, 2.000243111e+0, 4.262497508e-1, 1.663217177e-1, 1.189050200e-1, 4.662622228e-1, 1.915075719e-1, 4.327619272e-1, 1.838858552e-1, 2.920461109e-1, 1.823461280e-1, 8.293863702e-2, 2.074585819e-1, 1.001952178e-1, 1.588100623e-1},
    {2.652034102e-1, 1.607280687e-1, 7.643320863e-2, 5.081097451e-2, 3.446076176e-2, 1.643718620e-1, 8.655487938e-2, 7.553156064e-2, 4.958248334e-2, 1.187568890e-1, 8.197240896e-2, 4.590326388e-2, 9.298346885e-2, 5.447006630e-2, 7.646521729e-2},
    {8.801862774e-2, 6.125744532e-2, 3.760545063e-2, 2.829066600e-2, 1.974798796e-2, 6.543927065e-2, 4.184253862e-2, 3.706272183e-2, 2.750222273e-2, 5.286755896e-2, 4.000634951e-2, 2.628744797e-2, 4.473508853e-2, 3.037569283e-2, 3.898703611e-2}
};
constexpr double pCoeff4[4][15] = {
    {5.675242080e-2, -1.198411747e-2, -3.295496352e-2, -1.120682822e-1, 1.103657561e-2, 5.713170255e-2, -1.434249391e-2, -6.035216774e-3, -1.801459207e-2, 5.799057705e-2, -2.816702620e-3, -2.421626009e-2, 5.902730589e-2, 9.174268830e-2, 6.010484250e-2},
    {2.601413550e-1, -5.472052539e-2, -1.724516959e-1, -2.845426863e-1, -5.606519023e-1, 2.857455515e-1, 2.755177589e-1, -6.013310874e-2, -1.360777372e-1, 3.045581349e-1, 2.177095871e-1, 3.937644956e-1, 3.191828952e-1, 4.023496947e-1, 3.309738329e-1},
    {5.328461143e-1, 5.805587176e-1, 7.518511194e-1, 8.909873788e-1, 1.179429987e+0, 5.517873105e-1, 5.846750879e-1, 6.451518200e-1, 7.533973719e-1, 5.601358038e-1, 6.058047348e-1, 5.489520286e-1, 5.639423893e-1, 4.937432100e-1, 5.655207585e-1},
    {2.916254405e-1, 4.770079976e-1, 3.589627317e-1, 3.517811205e-1, 1.734974376e-1, 2.632314924e-1, 2.144986514e-1, 4.117923820e-1, 3.409304859e-1, 2.432423313e-1, 2.717811257e-1, 1.190436963e-1, 2.284796537e-1, 1.254001522e-1, 2.171122608e-1}
};

// --- STO-6G tables ---
constexpr double pAlpha6[6][15] = {
    {2.310303149e+1, 2.768496241e+1, 3.273031938e+0, 3.232838646e+0, 1.410128298e+0, 5.868285913e+0, 5.077973607e+0, 2.389722618e+0, 3.778623374e+0, 2.488296923e+0, 4.634239420e+0, 8.820520428e-1, 1.357718039e+0, 1.334096840e+0, 8.574668996e-1},
    {4.235915534e+0, 5.077140627e+0, 9.200611311e-1, 3.605788802e-1, 5.077878915e-1, 1.530329631e+0, 1.340786940e+0, 7.960947826e-1, 3.499121109e-1, 7.981487853e-1, 1.341648295e+0, 3.410838409e-1, 5.004907278e-1, 2.372312347e-1, 3.497184772e-1},
    {1.185056519e+0, 1.426786050e+0, 3.593349765e-1, 1.717905487e-1, 1.847926858e-1, 5.475665231e-1, 2.248434849e-1, 3.415541380e-1, 1.683175469e-1, 3.311327490e-1, 2.209593028e-1, 9.204308840e-2, 2.296565064e-1, 1.269485744e-1, 1.727917060e-1},
    {4.070988982e-1, 2.040335729e-1, 8.636686991e-2, 5.277666487e-2, 1.061070594e-1, 2.288932733e-1, 1.131741848e-1, 8.847434525e-2, 5.404070736e-2, 1.559114463e-1, 1.101467943e-1, 5.472831774e-2, 1.173146814e-1, 7.290318381e-2, 9.373643151e-2},
    {1.580884151e-1, 9.260298399e-2, 4.797373812e-2, 3.163400284e-2, 3.669584901e-2, 1.046655969e-1, 6.076408893e-2, 4.958248334e-2, 3.328911801e-2, 7.877734732e-2, 5.904190370e-2, 3.391202830e-2, 6.350097171e-2, 4.351355997e-2, 5.340032759e-2},
    {6.510953954e-2, 4.416183978e-2, 2.724741144e-2, 1.874093091e-2, 2.213558430e-2, 4.948220127e-2, 3.315424265e-2, 2.816929784e-2, 2.063815019e-2, 4.058484363e-2, 3.232628887e-2, 2.108227374e-2, 3.474556673e-2, 2.598071843e-2, 3.057364464e-2}
};
constexpr double pCoeff6[6][15] = {
    {9.163596280e-3, -4.151277819e-3, -6.775596947e-3, 1.374817488e-3, 2.695439582e-3, 7.924233646e-3, -3.329929840e-3, -1.665913575e-3, 1.163246387e-4, 7.283828112e-3, -4.749842876e-4, -4.097377019e-3, 6.930234381e-3, -9.486751531e-4, 6.729778096e-3},
    {4.936149294e-2, -2.067024148e-2, -5.639325779e-2, -8.666390043e-2, 1.850157487e-2, 5.144104825e-2, -1.419488340e-2, -1.657464971e-2, -2.920771322e-2, 5.386799363e-2, -3.566777891e-3, -2.508271857e-2, 5.634263745e-2, 4.624275998e-2, 5.874145170e-2},
    {1.685383049e-1, -5.150303337e-2, -1.587856086e-1, -3.130627309e-1, -9.588628125e-2, 1.898400060e-1, 1.639395770e-1, -5.958513378e-2, -1.381051233e-1, 2.072139149e-1, 1.108670481e-1, 2.648458555e-1, 2.217065797e-1, 2.373699784e-1, 2.339955227e-1},
    {3.705627997e-1, 3.346271174e-1, 5.534527651e-1, 7.812787397e-1, -5.200673560e-1, 4.049863191e-1, 4.485358256e-1, 4.053115554e-1, 5.706134877e-1, 4.266269092e-1, 4.159646930e-1, 5.097437054e-1, 4.411388883e-1, 4.589112231e-1, 4.512983737e-1},
    {4.164915298e-1, 5.621061301e-1, 5.015351020e-1, 4.389247988e-1, 1.087619490e+0, 4.012362861e-1, 3.908813050e-1, 5.433958189e-1, 4.768808140e-1, 3.843100204e-1, 4.621672517e-1, 2.654483467e-1, 3.688112625e-1, 3.205010548e-1, 3.552053926e-1},
    {1.303340841e-1, 1.712994697e-1, 7.223633674e-2, 2.487178756e-2, 3.103964343e-1, 1.051855189e-1, 7.411456232e-2, 1.204970491e-1, 6.021665516e-2, 8.902827546e-2, 1.081250196e-1, 2.623132212e-2, 7.787514504e-2, 5.077063693e-2, 6.974153145e-2}
};

// --- STO-5G tables ---
constexpr double pAlpha5[5][15] = {
    {1.130563696e+1, 8.984956862e+0, 4.275877914e+0, 2.980263783e+0, 7.403763257e-1, 3.320386533e+0, 6.466803859e+0, 1.091977298e+0, 3.422168934e-1, 1.539033958e+0, 1.522122079e+0, 9.702946470e-1, 8.925960415e-1, 1.670735676e+0, 5.895429375e-1},
    {2.071728178e+0, 1.673710636e+0, 1.132409433e+0, 3.792228833e-1, 1.367990863e-1, 8.643257633e-1, 1.555914802e+0, 3.719985051e-1, 1.665099900e-1, 4.922090297e-1, 2.173041823e-1, 3.603270196e-1, 3.277589120e-1, 2.072477219e-1, 2.393343780e-1},
    {5.786484833e-1, 1.944726668e-1, 4.016256968e-1, 1.789717224e-1, 9.135301779e-2, 3.079819284e-1, 1.955925255e-1, 8.590019352e-2, 5.443732013e-2, 2.029756928e-1, 1.084876577e-1, 8.668717752e-2, 1.492869962e-1, 1.024709357e-1, 1.172646904e-1},
    {1.975724573e-1, 8.806345634e-2, 7.732370620e-2, 5.002110360e-2, 3.726907315e-2, 1.273309895e-1, 8.809647701e-2, 4.786503860e-2, 3.367775277e-2, 9.424112917e-2, 5.836797641e-2, 4.833708379e-2, 7.506099109e-2, 5.537913898e-2, 6.254074479e-2},
    {7.445271746e-2, 4.249068522e-2, 3.800708627e-2, 2.789361681e-2, 2.241490836e-2, 5.606243164e-2, 4.234835707e-2, 2.730479990e-2, 2.091949042e-2, 4.569058269e-2, 3.206682246e-2, 2.751899341e-2, 3.892475795e-2, 3.072866652e-2, 3.411243214e-2}
};
constexpr double pCoeff5[5][15] = {
    {2.214055312e-2, -1.596349096e-2, -3.920358850e-3, 1.513948997e-3, 1.375523371e-2, 2.079051117e-2, -2.329023747e-3, -1.143929558e-2, -3.113958289e-2, 2.020869128e-2, -3.673711876e-3, -3.231527611e-3, 1.999839052e-2, 1.909729355e-1, 1.780980905e-1},
    {1.135411520e-1, -5.685884883e-2, -4.168430506e-2, -7.316801518e-2, -3.097344179e-1, 1.235472099e-1, -1.357395221e-2, -6.322651538e-2, -1.374007017e-1, 1.321157923e-1, 1.167122499e-1, -2.434931372e-2, 1.395427440e-1, 6.146060459e-1, 6.063757846e-1},
    {3.318161484e-1, 3.698265599e-1, -1.637440990e-1, -3.143703799e-1, -3.199192259e-1, 3.667738886e-1, 2.632185383e-1, 4.398907721e-1, 5.573881018e-1, 3.911240346e-1, 4.216476416e-1, 3.440817054e-1, 4.091508237e-1, 3.059611271e-1, 3.828552923e-1},
    {4.825700713e-1, 5.480512593e-1, 7.419373723e-1, 9.032615169e-1, 1.084547038e+0, 4.834930290e-1, 5.880427024e-1, 5.245859166e-1, 4.855428100e-1, 4.779609701e-1, 4.547673415e-1, 5.693674376e-1, 4.779609701e-1, 1.254001522e-1, 2.171122608e-1},
    {1.935721966e-1, 1.472634893e-1, 3.724364929e-1, 3.294210848e-1, 3.345288361e-1, 1.653444074e-1, 2.242794445e-1, 1.017072253e-1, 6.605423564e-2, 1.463662294e-1, 1.037803318e-1, 1.511340183e-1, 1.463662294e-1, 5.077063693e-2, 6.974153145e-2}
};

/**
 * @brief Convert an STO basis function to a contracted GTO (STO-NG expansion)
 *
 * Uses Stewart's STO-NG coefficients (J. Chem. Phys. 52, 431-438, 1970).
 * TBLite applies normalization to coefficients so the overlap between
 * normalized primitives gives the correct CGTO overlap directly.
 *
 * @param n Principal quantum number (1-5)
 * @param l Angular momentum (0=s, 1=p, 2=d)
 * @param zeta Slater exponent
 * @param ng Number of Gaussians (3, 4, 5, or 6)
 * @return CGTO Shell with normalized coefficients
 */
inline Shell slater_to_gauss(int n, int l, double zeta, int ng)
{
    Shell s;
    s.ang = l;
    s.nprim = ng;
    s.alpha.resize(ng);
    s.coeff.resize(ng);

    int ityp = nlm_to_ityp(n, l);
    if (ityp < 0 || ityp >= 15) {
        // Fallback: single Gaussian with zeta^2 exponent
        s.nprim = 1;
        s.alpha = {zeta * zeta};
        s.coeff = {1.0};
        return s;
    }

    // Get raw exponents and coefficients from Stewart tables
    const double (*alpha_table)[15] = nullptr;
    const double (*coeff_table)[15] = nullptr;

    switch (ng) {
    case 3: alpha_table = pAlpha3; coeff_table = pCoeff3; break;
    case 4: alpha_table = pAlpha4; coeff_table = pCoeff4; break;
    case 5: alpha_table = pAlpha5; coeff_table = pCoeff5; break;
    case 6: alpha_table = pAlpha6; coeff_table = pCoeff6; break;
    default:
        // Fallback to STO-6G
        alpha_table = pAlpha6; coeff_table = pCoeff6; ng = 6;
        s.nprim = 6; s.alpha.resize(6); s.coeff.resize(6);
        break;
    }

    double zeta2 = zeta * zeta;
    for (int i = 0; i < ng; ++i) {
        s.alpha[i] = alpha_table[i][ityp] * zeta2;
        s.coeff[i] = coeff_table[i][ityp];
    }

    // Apply normalization (TBLite basis/slater.f90 line 504-505, norm=.true.)
    // N = (2α/π)^(3/4) * (4α)^(l/2) / sqrt((2l-1)!!)
    // Note: dfactorial index is l (C++ 0-based), NOT l+1
    // Fortran: dfactorial(l+1) is 1-indexed, so dfactorial(1+1)=dfactorial(2)=1.0 for l=1
    // C++:     dfactorial[l]   is 0-indexed, so dfactorial[1]=1.0 for l=1
    for (int i = 0; i < ng; ++i) {
        double a = s.alpha[i];
        s.coeff[i] *= std::pow(two_over_pi * a, 0.75)
                    * std::pow(4.0 * a, 0.5 * l)
                    / std::sqrt(dfactorial[l]);
    }

    return s;
}

/**
 * @brief Compute overlap between two normalized primitive Gaussians (s-type)
 *
 * S = (π/(α+β))^(3/2) * exp(-αβ/(α+β) * R²)
 * But since coefficients already contain normalization:
 * The full primitive overlap is just the un-normalized Gaussian overlap.
 */
inline double primitive_ss_overlap(double alpha_a, double alpha_b, double R2)
{
    double gamma = alpha_a + alpha_b;
    return std::pow(M_PI / gamma, 1.5) * std::exp(-alpha_a * alpha_b / gamma * R2);
}

/**
 * @brief Compute overlap between s and p normalized primitive Gaussians
 *
 * For px: S = (Px-Bx) * factor, where P is the Gaussian product center
 * For py, pz: similar with y, z components
 *
 * @param alpha_s Exponent of s-function
 * @param alpha_p Exponent of p-function
 * @param R2 Squared distance
 * @param component Direction (Px-Bx, Py-By, or Pz-Bz) from product center to p-center
 */
inline double primitive_sp_overlap(double alpha_s, double alpha_p, double R2,
                                    double PB_component)
{
    double gamma = alpha_s + alpha_p;
    double S00 = std::pow(M_PI / gamma, 1.5) * std::exp(-alpha_s * alpha_p / gamma * R2);
    return S00 * PB_component;
}

/**
 * @brief Compute overlap between two p-type normalized primitive Gaussians
 *
 * For same axis (e.g., px-px): S = (PA_x * PB_x + 1/(2γ)) * S_00
 * For different axes (e.g., px-py): S = PA_x * PB_y * S_00
 */
inline double primitive_pp_overlap(double alpha_a, double alpha_b, double R2,
                                    double PA_comp, double PB_comp, bool same_axis)
{
    double gamma = alpha_a + alpha_b;
    double S00 = std::pow(M_PI / gamma, 1.5) * std::exp(-alpha_a * alpha_b / gamma * R2);
    if (same_axis)
        return S00 * (PA_comp * PB_comp + 0.5 / gamma);
    else
        return S00 * PA_comp * PB_comp;
}

/**
 * @brief Compute contracted GTO overlap between two shells at given positions
 *
 * S = Σ_i Σ_j c_i * c_j * S_prim(α_i, α_j, R²)
 *
 * @param shell_a First shell
 * @param shell_b Second shell
 * @param xa,ya,za Position of shell a (Bohr)
 * @param xb,yb,zb Position of shell b (Bohr)
 * @param type_a Orbital type (0=s, 1=px, 2=py, 3=pz)
 * @param type_b Orbital type
 * @return Overlap integral
 */
inline double cgto_overlap(const Shell& shell_a, const Shell& shell_b,
                            double xa, double ya, double za,
                            double xb, double yb, double zb,
                            int type_a, int type_b)
{
    double dx = xb - xa, dy = yb - ya, dz = zb - za;
    double R2 = dx * dx + dy * dy + dz * dz;

    double S = 0.0;

    for (int i = 0; i < shell_a.nprim; ++i) {
        double ai = shell_a.alpha[i];
        double ci = shell_a.coeff[i];
        for (int j = 0; j < shell_b.nprim; ++j) {
            double aj = shell_b.alpha[j];
            double cj = shell_b.coeff[j];
            double gamma = ai + aj;

            // Gaussian product center
            double Px = (ai * xa + aj * xb) / gamma;
            double Py = (ai * ya + aj * yb) / gamma;
            double Pz = (ai * za + aj * zb) / gamma;

            if (type_a == 0 && type_b == 0) {
                // s-s
                S += ci * cj * primitive_ss_overlap(ai, aj, R2);
            } else if (type_a == 0 && type_b >= 1 && type_b <= 3) {
                // s-p
                double PB;
                if (type_b == 1) PB = Px - xb;       // px
                else if (type_b == 2) PB = Py - yb;   // py
                else PB = Pz - zb;                      // pz
                S += ci * cj * primitive_sp_overlap(ai, aj, R2, PB);
            } else if (type_a >= 1 && type_a <= 3 && type_b == 0) {
                // p-s (swap and negate direction)
                double PA;
                if (type_a == 1) PA = Px - xa;
                else if (type_a == 2) PA = Py - ya;
                else PA = Pz - za;
                S += ci * cj * primitive_sp_overlap(aj, ai, R2, PA);
            } else if (type_a >= 1 && type_a <= 3 && type_b >= 1 && type_b <= 3) {
                // p-p
                double PA, PB;
                int axis_a, axis_b;
                if (type_a == 1) { PA = Px - xa; axis_a = 0; }
                else if (type_a == 2) { PA = Py - ya; axis_a = 1; }
                else { PA = Pz - za; axis_a = 2; }
                if (type_b == 1) { PB = Px - xb; axis_b = 0; }
                else if (type_b == 2) { PB = Py - yb; axis_b = 1; }
                else { PB = Pz - zb; axis_b = 2; }
                S += ci * cj * primitive_pp_overlap(ai, aj, R2, PA, PB, axis_a == axis_b);
            }
        }
    }

    return S;
}

/**
 * @brief Gram-Schmidt orthogonalize shell_b to shell_a (same angular momentum)
 *
 * TBLite basis/ortho.f90: When an atom has two shells with the same angular momentum
 * (e.g., H has 1s and 2s), the second shell must be orthogonalized to the first.
 * Without this, the overlap matrix is near-singular (S(1s,2s) ≈ 0.98).
 *
 * The orthogonalized shell_b gets shell_a's primitives appended with -overlap * coeff,
 * effectively encoding: |2s'⟩ = (|2s⟩ - S₁₂|1s⟩) / √(1 - S₁₂²)
 *
 * @param shell_a First shell (unchanged)
 * @param shell_b Second shell (modified in place — gets more primitives)
 */
inline void orthogonalize(const Shell& shell_a, Shell& shell_b)
{
    if (shell_a.ang != shell_b.ang) return;

    // Step 1: Compute overlap ⟨a|b⟩ between the two contracted shells (same center)
    // For same-center s-functions: S = Σ_ij c_i * c_j * (π/(α_i+α_j))^(3/2)
    double overlap = 0.0;
    for (int i = 0; i < shell_a.nprim; ++i) {
        for (int j = 0; j < shell_b.nprim; ++j) {
            double eab = shell_a.alpha[i] + shell_b.alpha[j];
            double kab = std::pow(M_PI / eab, 1.5);
            overlap += shell_a.coeff[i] * shell_b.coeff[j] * kab;
        }
    }

    // Step 2: Append shell_a's primitives to shell_b with -overlap * coefficient
    // |b'⟩ = |b⟩ - overlap * |a⟩
    int old_nprim = shell_b.nprim;
    int new_nprim = old_nprim + shell_a.nprim;
    shell_b.alpha.resize(new_nprim);
    shell_b.coeff.resize(new_nprim);
    for (int i = 0; i < shell_a.nprim; ++i) {
        shell_b.alpha[old_nprim + i] = shell_a.alpha[i];
        shell_b.coeff[old_nprim + i] = -overlap * shell_a.coeff[i];
    }
    shell_b.nprim = new_nprim;

    // Step 3: Renormalize — compute ⟨b'|b'⟩ and divide coefficients by √(⟨b'|b'⟩)
    double norm = 0.0;
    for (int i = 0; i < shell_b.nprim; ++i) {
        for (int j = 0; j < shell_b.nprim; ++j) {
            double eab = shell_b.alpha[i] + shell_b.alpha[j];
            double kab = std::pow(M_PI / eab, 1.5);
            norm += shell_b.coeff[i] * shell_b.coeff[j] * kab;
        }
    }
    double inv_sqrt_norm = 1.0 / std::sqrt(norm);
    for (int i = 0; i < shell_b.nprim; ++i) {
        shell_b.coeff[i] *= inv_sqrt_norm;
    }
}

} // namespace CGTO
