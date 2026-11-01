/*
 * <Geometric counterpoise (gCP) + short-range basis (SRB) correction>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Ported from dftd3/simple-dftd3 src/dftd3/gcp.f90 (gcp_energy, gcp_deriv,
 * srb_energy, srb_deriv, dsovl, aaux_all, baux_all, bint_all) and
 * src/dftd3/gcp/param.f90 (emiss_hf_minix, nbas_minix, slater_s/p/d, the
 * p_minix_bas branch, effective_atomic_number), LGPL-3.0-or-later. D3 pair radii from
 * src/dftd3/data/vdwrad.f90. See gcp.h for the model.
 *
 * Claude Generated: native gCP/SRB for HF-3c (Sep 2026)
 *
 * This program is free software under GPL-3.0
 */

#include "gcp.h"

#include "src/core/units.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace gcp {
namespace {

// Claude Generated (Nov 2026, stage 5a): tables extended from H-Ne to H-Kr, the range of the reference
// tables (param.f90 emiss_hf_minix, nbas_minix: 36 entries); heavier elements use an effective atomic number
// (effectiveZ, param.f90 effective_atomic_number), exactly as the reference does.
constexpr int kMaxZ = 36;      // H-Kr, the range of the MINIX gCP tables
constexpr int kMaxAtomZ = 86;  // H-Rn: elements accepted through effectiveZ (the engine's range)

// BSSE per atom for HF/MINIX: param.f90 emiss_hf_minix =
// [emiss_hf_minis(1:2), 0.177871, 0.171596, emiss_hf_minis(5:10), 1.114110, 1.271150,
//  emiss_hf_minisd(13:18), emiss_hf_sv(19:30), emiss_hf_svp(31:36)].
// NOTE: the digits matter -- a table rounded to 5 places (the earlier Python
// witness) is off by up to 2.6e-7 Eh on NH3.
// Ga-Kr (emiss_hf_svp): the reference writes these literals WITHOUT the `_wp` kind (default REAL),
// so they are single precision widened to double; kEmissSvpF keeps that (see kSlaterS below).
constexpr float kEmissSvpF[6] = { 0.350021f, 0.345779f, 0.349532f, 0.367305f, 0.382008f, 0.399709f };
constexpr double kEmissMinix[kMaxZ + 1] = {
    0.0,
    0.042400, 0.028324,
    0.177871, 0.171596, 0.224237, 0.279950, 0.357906, 0.479012, 0.638518, 0.832349,
    1.114110, 1.271150, 1.446950, 1.610980, 1.766610, 1.988230, 2.228450, 2.487960,
    0.374252, 0.460972,
    0.444886, 0.404993, 0.378406, 0.373439, 0.361245, 0.360014, 0.362928, 0.243801, 0.405299, 0.396510,
    kEmissSvpF[0], kEmissSvpF[1], kEmissSvpF[2], kEmissSvpF[3], kEmissSvpF[4], kEmissSvpF[5]
};

// Number of basis functions per atom in MINIX (param.f90 nbas_minix): the reference's numbers, not
// the size of the actual basis (K-Zn take the SV count, Ga-Kr the SVP count).
constexpr double kNbasMinix[kMaxZ + 1] = {
    0, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5,
    9, 9, 14, 14, 14, 14, 14, 14,
    11, 11, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
    32, 32, 32, 32, 32, 32
};

// Slater exponents (param.f90 slater_s/p/d = gcp setzet ZS/ZP/ZD); slaterExponent() forms the average.
// FLOAT on purpose: the reference declares `[real(wp) :: 1.2000, ...]` with no
// `_wp` kind suffix, so Fortran reads each literal as default (single-precision)
// REAL and only then widens it -- e.g. 2.5644 becomes 2.5643999576568604. Using
// the exact decimals instead shifts E_gCP by ~2e-8 relative (6e-10 Eh on F2).
// The averages are formed in double, like `(slater_s + slater_p)/2` in wp.
constexpr float kSlaterS[kMaxZ + 1] = {
    0.0f, 1.2000f, 1.6469f, 0.6534f, 1.0365f, 1.3990f, 1.7210f, 2.0348f, 2.2399f, 2.5644f, 2.8812f,
    0.8675f, 1.1935f, 1.5143f, 1.7580f, 1.9860f, 2.1362f, 2.3617f, 2.5796f, 0.9362f, 1.2112f,
    1.2870f, 1.3416f, 1.3570f, 1.3804f, 1.4761f, 1.5465f, 1.5650f, 1.5532f, 1.5781f, 1.7778f,
    2.0675f, 2.2702f, 2.4546f, 2.5680f, 2.7523f, 2.9299f
};
constexpr float kSlaterP[kMaxZ + 1] = {
    0.0f, 0.0000f, 0.0000f, 0.5305f, 0.8994f, 1.2685f, 1.6105f, 1.9398f, 2.0477f, 2.4022f, 2.7421f,
    0.6148f, 0.8809f, 1.1660f, 1.4337f, 1.6755f, 1.7721f, 2.0176f, 2.2501f, 0.6914f, 0.9329f,
    0.9828f, 1.0104f, 0.9947f, 0.9784f, 1.0641f, 1.1114f, 1.1001f, 1.0594f, 1.0527f, 1.2448f,
    1.5073f, 1.7680f, 1.9819f, 2.0548f, 2.2652f, 2.4617f
};
// slater_d: only Sc-Zn are used (both references have unused entries for Ga-Kr).
constexpr float kSlaterD[kMaxZ + 1] = {
    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
    2.4341f, 2.6439f, 2.7809f, 2.9775f, 3.2208f, 3.4537f, 3.6023f, 3.7017f, 3.8962f, 2.0477f,
    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
};

// D3 van-der-Waals pair radii R0_AB in Angstrom, packed lower triangle
// (vdwrad.f90: index num1 + num2*(num2-1)/2 with num1 <= num2, 1-based), H-Kr.
constexpr double kR0abAngstrom[kMaxZ * (kMaxZ + 1) / 2] = {
    2.1823, 1.8547, 1.7347, 2.9086, 2.5732, 3.4956, 2.3550, 2.5095, 2.9802, 3.0982,
    2.5141, 2.3917, 2.9977, 2.9484, 3.2160, 2.4492, 2.2527, 3.1933, 3.0214, 2.9531,
    2.9103, 2.3667, 2.1328, 2.8784, 2.7660, 2.7776, 2.7063, 2.6225, 2.1768, 2.0625,
    2.6395, 2.6648, 2.6482, 2.5697, 2.4846, 2.4817, 2.0646, 1.9891, 2.5086, 2.6908,
    2.6233, 2.4770, 2.3885, 2.3511, 2.2996, 1.9892, 1.9251, 2.4190, 2.5473, 2.4994,
    2.4091, 2.3176, 2.2571, 2.1946, 2.1374, 2.9898, 2.6397, 3.6031, 3.1219, 3.7620,
    3.2485, 2.9357, 2.7093, 2.5781, 2.4839, 3.7082, 2.5129, 2.7321, 3.1052, 3.2962,
    3.1331, 3.2000, 2.9586, 3.0822, 2.8582, 2.7120, 3.2570, 3.4839, 2.8766, 2.7427,
    3.2776, 3.2363, 3.5929, 3.2826, 3.0911, 2.9369, 2.9030, 2.7789, 3.3921, 3.3970,
    4.0106, 2.8884, 2.6605, 3.7513, 3.1613, 3.3605, 3.3325, 3.0991, 2.9297, 2.8674,
    2.7571, 3.8129, 3.3266, 3.7105, 3.7917, 2.8304, 2.5538, 3.3932, 3.1193, 3.1866,
    3.1245, 3.0465, 2.8727, 2.7664, 2.6926, 3.4608, 3.2984, 3.5142, 3.5418, 3.5017,
    2.6190, 2.4797, 3.1331, 3.0540, 3.0651, 2.9879, 2.9054, 2.8805, 2.7330, 2.6331,
    3.2096, 3.5668, 3.3684, 3.3686, 3.3180, 3.3107, 2.4757, 2.4019, 2.9789, 3.1468,
    2.9768, 2.8848, 2.7952, 2.7457, 2.6881, 2.5728, 3.0574, 3.3264, 3.3562, 3.2529,
    3.1916, 3.1523, 3.1046, 2.3725, 2.3289, 2.8760, 2.9804, 2.9093, 2.8040, 2.7071,
    2.6386, 2.5720, 2.5139, 2.9517, 3.1606, 3.2085, 3.1692, 3.0982, 3.0352, 2.9730,
    2.9148, 3.2147, 2.8315, 3.8724, 3.4621, 3.8823, 3.3760, 3.0746, 2.8817, 2.7552,
    2.6605, 3.9740, 3.6192, 3.6569, 3.9586, 3.6188, 3.3917, 3.2479, 3.1434, 4.2411,
    2.7597, 3.0588, 3.3474, 3.6214, 3.4353, 3.4729, 3.2487, 3.3200, 3.0914, 2.9403,
    3.4972, 3.7993, 3.6773, 3.8678, 3.5808, 3.8243, 3.5826, 3.4156, 3.8765, 4.1035,
    2.7361, 2.9765, 3.2475, 3.5004, 3.4185, 3.4378, 3.2084, 3.2787, 3.0604, 2.9187,
    3.4037, 3.6759, 3.6586, 3.8327, 3.5372, 3.7665, 3.5310, 3.3700, 3.7788, 3.9804,
    3.8903, 2.6832, 2.9060, 3.2613, 3.4359, 3.3538, 3.3860, 3.1550, 3.2300, 3.0133,
    2.8736, 3.4024, 3.6142, 3.5979, 3.5295, 3.4834, 3.7140, 3.4782, 3.3170, 3.7434,
    3.9623, 3.8181, 3.7642, 2.6379, 2.8494, 3.1840, 3.4225, 3.2771, 3.3401, 3.1072,
    3.1885, 2.9714, 2.8319, 3.3315, 3.5979, 3.5256, 3.4980, 3.4376, 3.6714, 3.4346,
    3.2723, 3.6859, 3.8985, 3.7918, 3.7372, 3.7211, 2.9230, 2.6223, 3.4161, 2.8999,
    3.0557, 3.3308, 3.0555, 2.8508, 2.7385, 2.6640, 3.5263, 3.0277, 3.2990, 3.7721,
    3.5017, 3.2751, 3.1368, 3.0435, 3.7873, 3.2858, 3.2140, 3.1727, 3.2178, 3.4414,
    2.5490, 2.7623, 3.0991, 3.3252, 3.1836, 3.2428, 3.0259, 3.1225, 2.9032, 2.7621,
    3.2490, 3.5110, 3.4429, 3.3845, 3.3574, 3.6045, 3.3658, 3.2013, 3.6110, 3.8241,
    3.7090, 3.6496, 3.6333, 3.0896, 3.5462, 2.4926, 2.7136, 3.0693, 3.2699, 3.1272,
    3.1893, 2.9658, 3.0972, 2.8778, 2.7358, 3.2206, 3.4566, 3.3896, 3.3257, 3.2946,
    3.5693, 3.3312, 3.1670, 3.5805, 3.7711, 3.6536, 3.5927, 3.5775, 3.0411, 3.4885,
    3.4421, 2.4667, 2.6709, 3.0575, 3.2357, 3.0908, 3.1537, 2.9235, 3.0669, 2.8476,
    2.7054, 3.2064, 3.4519, 3.3593, 3.2921, 3.2577, 3.2161, 3.2982, 3.1339, 3.5606,
    3.7582, 3.6432, 3.5833, 3.5691, 3.0161, 3.4812, 3.4339, 3.4327, 2.4515, 2.6338,
    3.0511, 3.2229, 3.0630, 3.1265, 2.8909, 3.0253, 2.8184, 2.6764, 3.1968, 3.4114,
    3.3492, 3.2691, 3.2320, 3.1786, 3.2680, 3.1036, 3.5453, 3.7259, 3.6090, 3.5473,
    3.5327, 3.0018, 3.4413, 3.3907, 3.3593, 3.3462, 2.4413, 2.6006, 3.0540, 3.1987,
    3.0490, 3.1058, 2.8643, 2.9948, 2.7908, 2.6491, 3.1950, 3.3922, 3.3316, 3.2585,
    3.2136, 3.1516, 3.2364, 3.0752, 3.5368, 3.7117, 3.5941, 3.5313, 3.5164, 2.9962,
    3.4225, 3.3699, 3.3370, 3.3234, 3.3008, 2.4318, 2.5729, 3.0416, 3.1639, 3.0196,
    3.0843, 2.8413, 2.7436, 2.7608, 2.6271, 3.1811, 3.3591, 3.3045, 3.2349, 3.1942,
    3.1291, 3.2111, 3.0534, 3.5189, 3.6809, 3.5635, 3.5001, 3.4854, 2.9857, 3.3897,
    3.3363, 3.3027, 3.2890, 3.2655, 3.2309, 2.8502, 2.6934, 3.2467, 3.1921, 3.5663,
    3.2541, 3.0571, 2.9048, 2.8657, 2.7438, 3.3547, 3.3510, 3.9837, 3.6871, 3.4862,
    3.3389, 3.2413, 3.1708, 3.6096, 3.6280, 3.6860, 3.5568, 3.4836, 3.2868, 3.3994,
    3.3476, 3.3170, 3.2950, 3.2874, 3.2606, 3.9579, 2.9226, 2.6838, 3.7867, 3.1732,
    3.3872, 3.3643, 3.1267, 2.9541, 2.8505, 2.7781, 3.8475, 3.3336, 3.7359, 3.8266,
    3.5733, 3.3959, 3.2775, 3.1915, 3.9878, 3.8816, 3.5810, 3.5364, 3.5060, 3.8097,
    3.3925, 3.3348, 3.3019, 3.2796, 3.2662, 3.2464, 3.7136, 3.8619, 2.9140, 2.6271,
    3.4771, 3.1774, 3.2560, 3.1970, 3.1207, 2.9406, 2.8322, 2.7571, 3.5455, 3.3514,
    3.5837, 3.6177, 3.5816, 3.3902, 3.2604, 3.1652, 3.7037, 3.6283, 3.5858, 3.5330,
    3.4884, 3.5789, 3.4094, 3.3473, 3.3118, 3.2876, 3.2707, 3.2521, 3.5570, 3.6496,
    3.6625, 2.7300, 2.5870, 3.2471, 3.1487, 3.1667, 3.0914, 3.0107, 2.9812, 2.8300,
    2.7284, 3.3259, 3.3182, 3.4707, 3.4748, 3.4279, 3.4182, 3.2547, 3.1353, 3.5116,
    3.9432, 3.8828, 3.8303, 3.7880, 3.3760, 3.7218, 3.3408, 3.3059, 3.2698, 3.2446,
    3.2229, 3.4422, 3.5023, 3.5009, 3.5268, 2.6026, 2.5355, 3.1129, 3.2863, 3.1029,
    3.0108, 2.9227, 2.8694, 2.8109, 2.6929, 3.1958, 3.4670, 3.4018, 3.3805, 3.3218,
    3.2815, 3.2346, 3.0994, 3.3937, 3.7266, 3.6697, 3.6164, 3.5730, 3.2522, 3.5051,
    3.4686, 3.4355, 3.4084, 3.3748, 3.3496, 3.3692, 3.4052, 3.3910, 3.3849, 3.3662,
    2.5087, 2.4814, 3.0239, 3.1312, 3.0535, 2.9457, 2.8496, 2.7780, 2.7828, 2.6532,
    3.1063, 3.3143, 3.3549, 3.3120, 3.2421, 3.1787, 3.1176, 3.0613, 3.3082, 3.5755,
    3.5222, 3.4678, 3.4231, 3.1684, 3.3528, 3.3162, 3.2827, 3.2527, 3.2308, 3.2029,
    3.3173, 3.3343, 3.3092, 3.2795, 3.2452, 3.2096
};

// Effective atomic number of the reference (param.f90 effective_atomic_number): elements beyond Kr take
// the parameters of the lighter element with the same valence shell (Rb-Xe -> K-Kr, Cs-La -> K-Sc,
// Ce-Lu -> Sc, Hf-Rn -> Ti-Kr). Used for emiss, nbas, the electron count, the Slater exponent, the
// shell, and (for HF-3c's base term) the charge and radius of the SRB pair. Identity for H-Kr.
int effectiveZ(int Z)
{
    if (Z >= 37 && Z <= 54) return Z - 18;
    if (Z >= 55 && Z <= 57) return Z - 36;
    if (Z >= 58 && Z <= 71) return 21;
    if (Z >= 72 && Z <= 89) return Z - 50;
    return Z;
}

// Z is the effective atomic number here (and in the helpers below).
// Claude Generated (Nov 2026): s for H/He, (s+p)/2 for Li-Ca AND Ga-Kr, (s+p+d)/3 for Sc-Zn -- the rule of the
// gCP program ORCA 6.1 runs (gcp setzet: `case(3:20,31:) (ZS+ZP)/2d0`, `case(21:30) (ZS+ZP+ZD)/3d0`, data from
// default-REAL literals as below). simple-dftd3 (param.f90 slater_exp) uses (s+2d)/3 and s instead: E_gCP of Sc-Kr
// differs by up to 1e-3 Eh between the two references (reference_gcp.json); ORCA's rule wins. H-Ca are the same.
double slaterExponent(int Z)
{
    if (Z <= 2) return static_cast<double>(kSlaterS[Z]);
    if (Z >= 21 && Z <= 30)
        return (static_cast<double>(kSlaterS[Z]) + static_cast<double>(kSlaterP[Z]) + static_cast<double>(kSlaterD[Z])) / 3.0;
    return 0.5 * (static_cast<double>(kSlaterS[Z]) + static_cast<double>(kSlaterP[Z]));
}

double r0abBohr(int za, int zb)
{
    const int lo = std::min(za, zb), hi = std::max(za, zb);
    return CurcumaUnit::Length::angstrom_to_bohr(kR0abAngstrom[lo + hi * (hi - 1) / 2 - 1]);
}

// Principal quantum number of the valence s shell used in the overlap (1 for
// H/He, 2 for Li-Ne, 3 from Na on: "4s, 5s will be treated as 3s"): the reference `shell` data.
int shellOf(int Z) { return (Z <= 2) ? 1 : (Z <= 10 ? 2 : 3); }

// A_k(x) = int_1^inf t^k e^{-xt} dt, upward recursion (aaux_all).
void aauxAll(double x, double* a, int kmax)
{
    const double ex = std::exp(-x), rx = 1.0 / x;
    a[0] = ex * rx;
    for (int k = 1; k <= kmax; ++k)
        a[k] = (k * a[k - 1] + ex) * rx;
}

// B_k(x) = int_-1^1 t^k e^{-xt} dt, closed form (baux_all).
void bauxAll(double x, double* b, int kmax)
{
    const double ep = std::exp(x), em = std::exp(-x), rx = 1.0 / x;
    for (int k = 0; k <= kmax; ++k) {
        double term = rx;
        double sgn = (k % 2 == 0) ? 1.0 : -1.0;
        double sp = sgn * term, sm = term;
        for (int j = 1; j <= k; ++j) {
            term *= (k - j + 1) * rx;
            sgn = -sgn;
            sp += sgn * term;
            sm += term;
        }
        b[k] = ep * sp - em * sm;
    }
}

// B_k(x) as the truncated 12-term Taylor series (bint_all). The reference uses it
// whenever the two Slater exponents differ by less than 0.1, where the closed
// form cancels badly -- this includes like-element pairs (x == 0) and, in H-Ne,
// the He-C pair. Ported with the same truncation so the energy matches exactly.
void bintAll(double x, double* b, int kmax)
{
    if (std::abs(x) < 1.0e-6) {
        for (int k = 0; k <= kmax; ++k)
            b[k] = (k % 2 == 0) ? 2.0 / (k + 1.0) : 0.0;
        return;
    }
    constexpr int nterm = 12;
    double pw[nterm + 1];
    pw[0] = 1.0;
    for (int i = 1; i <= nterm; ++i)
        pw[i] = pw[i - 1] * (-x) / i;
    for (int k = 0; k <= kmax; ++k) {
        double acc = 0.0;
        for (int i = k % 2; i <= nterm; i += 2)
            acc += pw[i] / (k + i + 1.0);
        b[k] = 2.0 * acc;
    }
}

// Overlap S(R) of two normalized ns Slater functions and dS/dR (dsovl), in the
// elliptic-coordinate form S = cnorm R^m sum_t w_t A_{p_t}(ax) B_{q_t}(bx).
// Za/Zb are the effective element numbers (they select 1s, 2s or 3s), za/zb the exponents.
void slaterOverlap(double r, int Za, int Zb, double za, double zb, double& s0, double& s1)
{
    const bool lsame = std::abs(za - zb) < 0.1;
    const int ii = shellOf(Za) * shellOf(Zb);

    int m = 0, nterm = 0;
    double wt[6] = {}, cnorm = 0.0;
    int pa[6] = {}, qb[6] = {};
    switch (ii) {
    case 1:  // <1s|1s>
        m = 3; nterm = 2;
        wt[0] = 1.0; wt[1] = -1.0;
        pa[0] = 2; pa[1] = 0;
        qb[0] = 0; qb[1] = 2;
        cnorm = 0.25 * std::sqrt(std::pow(za * zb, 3));
        break;
    case 2:  // <1s|2s>: the 1s function is always "a"
        if (shellOf(Za) >= shellOf(Zb)) std::swap(za, zb);
        m = 4; nterm = 4;
        wt[0] = 1.0; wt[1] = -1.0; wt[2] = 1.0; wt[3] = -1.0;
        pa[0] = 3; pa[1] = 0; pa[2] = 2; pa[3] = 1;
        qb[0] = 0; qb[1] = 3; qb[2] = 1; qb[3] = 2;
        cnorm = std::sqrt(1.0 / 3.0) * std::sqrt(std::pow(za, 3) * std::pow(zb, 5)) * 0.125;
        break;
    case 3:  // <1s|3s>: the 1s function is always "a" (Claude Generated, Nov 2026, dsovl case(3))
        if (shellOf(Za) >= shellOf(Zb)) std::swap(za, zb);
        m = 5; nterm = 4;
        wt[0] = 1.0; wt[1] = -1.0; wt[2] = 2.0; wt[3] = -2.0;
        pa[0] = 4; pa[1] = 0; pa[2] = 3; pa[3] = 1;
        qb[0] = 0; qb[1] = 4; qb[2] = 1; qb[3] = 3;
        cnorm = std::sqrt(std::pow(za, 3) * std::pow(zb, 7) / 7.5) * 0.0625 / std::sqrt(3.0);
        break;
    case 4:  // <2s|2s>
        m = 5; nterm = 3;
        wt[0] = 1.0; wt[1] = 1.0; wt[2] = -2.0;
        pa[0] = 4; pa[1] = 0; pa[2] = 2;
        qb[0] = 0; qb[1] = 4; qb[2] = 2;
        cnorm = std::sqrt(std::pow(za * zb, 5)) * 0.0625 / 3.0;
        break;
    case 6:  // <2s|3s>: the 2s function is always "a" (dsovl case(6))
        if (shellOf(Za) >= shellOf(Zb)) std::swap(za, zb);
        m = 6; nterm = 6;
        wt[0] = 1.0; wt[1] = 1.0; wt[2] = -2.0; wt[3] = -2.0; wt[4] = 1.0; wt[5] = 1.0;
        pa[0] = 5; pa[1] = 4; pa[2] = 3; pa[3] = 2; pa[4] = 1; pa[5] = 0;
        qb[0] = 0; qb[1] = 1; qb[2] = 2; qb[3] = 3; qb[4] = 4; qb[5] = 5;
        cnorm = std::sqrt(std::pow(za, 5) * std::pow(zb, 7) / 7.5) * 0.03125 / 3.0;
        break;
    case 9:  // <3s|3s> (dsovl case(9))
        m = 7; nterm = 4;
        wt[0] = 1.0; wt[1] = -3.0; wt[2] = 3.0; wt[3] = -1.0;
        pa[0] = 6; pa[1] = 4; pa[2] = 2; pa[3] = 0;
        qb[0] = 0; qb[1] = 2; qb[2] = 4; qb[3] = 6;
        cnorm = std::sqrt(std::pow(za * zb, 7)) / 1440.0;
        break;
    default:
        throw std::runtime_error("gCP: invalid shell combination");
    }

    const double ha = 0.5 * (za + zb), hb = 0.5 * (zb - za);
    double av[9], bv[9];
    aauxAll(ha * r, av, 8);
    if (lsame)
        bintAll(hb * r, bv, 8);
    else
        bauxAll(hb * r, bv, 8);

    double f0 = 0.0, f1 = 0.0;
    for (int t = 0; t < nterm; ++t) {
        const int p = pa[t], q = qb[t];
        f0 += wt[t] * av[p] * bv[q];
        f1 -= wt[t] * (ha * av[p + 1] * bv[q] + hb * av[p] * bv[q + 1]);
    }
    const double rm = std::pow(r, m);
    s0 = cnorm * rm * f0;
    s1 = cnorm * (m * std::pow(r, m - 1) * f0 + rm * f1);
}

// Number of virtual MINIX functions of atom Z, as 1/sqrt(nbas - nel/2) (0 if
// the atom has none, e.g. He with a single 1s function).
double inverseSqrtVirtuals(int Z)
{
    const double xv = kNbasMinix[Z] - 0.5 * Z;
    return (xv >= 0.5) ? 1.0 / std::sqrt(xv) : 0.0;
}

// HF-3c base term of one pair, with the EFFECTIVE atomic numbers for the charge product and the D3 radius.
// Claude Generated (Nov 2026): pairs with an element beyond Ar contribute nothing. This is the rule of the gCP
// program ORCA 6.1 runs for `! HF-3c` (otool_gcp 1.06; gcp basegrad: `if(iz(i).gt.18) cycle`, the same in the
// grimme-lab/gcp source of 2016), checked element by element (test_cases/qm_hf3c/reference_gcp.json).
// simple-dftd3 (gcp.f90 srb_energy, the H-Ne oracle of Sep 2026) has no such cut: a reference split, ORCA wins.
double srbPair(const std::vector<int>& atoms, const Matrix& xyz, const Parameters& p,
               int i, int j, Matrix* gradient)
{
    const int zi = effectiveZ(atoms[i]), zj = effectiveZ(atoms[j]);
    if (zi > 18 || zj > 18) return 0.0;
    const double dx = xyz(i, 0) - xyz(j, 0), dy = xyz(i, 1) - xyz(j, 1), dz = xyz(i, 2) - xyz(j, 2);
    const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (r > p.cutoff || r < 1.0e-14) return 0.0;
    const double r0 = p.rscal * std::pow(r0abBohr(zi, zj), 0.75);
    const double ff = -std::pow(static_cast<double>(zi) * zj, 1.5);
    const double expt = std::exp(-r0 * r);
    if (gradient) {
        const double g = -p.qscal * ff * r0 * expt / r;  // (dE/dR) / R
        (*gradient)(i, 0) += g * dx; (*gradient)(j, 0) -= g * dx;
        (*gradient)(i, 1) += g * dy; (*gradient)(j, 1) -= g * dy;
        (*gradient)(i, 2) += g * dz; (*gradient)(j, 2) -= g * dz;
    }
    return p.qscal * ff * expt;
}

} // namespace

Parameters hf3c()
{
    // gcp/param.f90, case(p_minix_bas) with method hf3c: HF/MINIX gCP parameters
    // and base = .true. with rscal = 0.7, qscal = 0.03 (Sure & Grimme 2013).
    Parameters p;
    p.sigma = 0.1290;
    p.eta = 1.1526;
    p.alpha = 1.1549;
    p.beta = 1.1763;
    p.base = true;
    p.rscal = 0.7;
    p.qscal = 0.03;
    return p;
}

bool supports(const std::vector<int>& atoms)
{
    for (int Z : atoms)
        if (Z < 1 || Z > kMaxAtomZ) return false;
    return true;
}

double baseEnergy(const std::vector<int>& atoms, const Matrix& xyz, const Parameters& p)
{
    if (!p.base) return 0.0;
    double e = 0.0;
    const int n = static_cast<int>(atoms.size());
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < i; ++j)
            e += srbPair(atoms, xyz, p, i, j, nullptr);
    return e;
}

double energy(const std::vector<int>& atoms, const Matrix& xyz, const Parameters& p,
              Matrix* gradient)
{
    if (!supports(atoms))
        throw std::runtime_error("gCP: element outside H-Rn (MINIX parameters)");

    const int n = static_cast<int>(atoms.size());
    if (gradient) *gradient = Matrix::Zero(n, 3);

    double e = 0.0;
    for (int i = 0; i < n; ++i) {
        const int Zi = effectiveZ(atoms[i]);
        const double xvi = inverseSqrtVirtuals(Zi);
        for (int j = 0; j < i; ++j) {
            const int Zj = effectiveZ(atoms[j]);
            const double xvj = inverseSqrtVirtuals(Zj);
            // BSSE that atom i suffers from j's virtuals, and vice versa.
            const double emij = p.sigma * (kEmissMinix[Zi] * xvj + kEmissMinix[Zj] * xvi);

            const double dx = xyz(i, 0) - xyz(j, 0), dy = xyz(i, 1) - xyz(j, 1), dz = xyz(i, 2) - xyz(j, 2);
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (r > p.cutoff || r < 1.0e-14) continue;

            double s = 0.0, ds = 0.0;
            slaterOverlap(r, Zi, Zj, p.eta * slaterExponent(Zi), p.eta * slaterExponent(Zj), s, ds);
            const double expv = std::exp(-p.alpha * std::pow(r, p.beta));
            const double bsse = expv / std::sqrt(s);
            e += emij * bsse;

            if (gradient) {
                // d/dR [exp(-a R^b) S^{-1/2}] = bsse * (-a b R^{b-1} - S'/(2S)); times vec/R.
                const double dEdR = emij * bsse * (-p.alpha * p.beta * std::pow(r, p.beta - 1.0) - 0.5 * ds / s);
                const double g = dEdR / r;
                (*gradient)(i, 0) += g * dx; (*gradient)(j, 0) -= g * dx;
                (*gradient)(i, 1) += g * dy; (*gradient)(j, 1) -= g * dy;
                (*gradient)(i, 2) += g * dz; (*gradient)(j, 2) -= g * dz;
            }
        }
    }

    if (p.base)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < i; ++j)
                e += srbPair(atoms, xyz, p, i, j, gradient);
    return e;
}

} // namespace gcp
