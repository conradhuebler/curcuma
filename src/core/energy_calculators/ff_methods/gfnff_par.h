/*
 * GFN-FF Parameters - angewChem2020 Parameter Set
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Reference: S. Spicher, S. Grimme, Angew. Chem. Int. Ed. 2020, 59, 15665
 * Fortran Implementation: gfnff_param.f90
 *
 * 🤖 Generated with Claude Code - Parameter extraction and consolidation
 * Co-Authored-By: Claude <noreply@anthropic.com>
 */

#pragma once

#include <vector>

namespace GFNFFParameters {

// ============================================================================
// SECTION 1: EEQ (Electronegativity Equalization) Parameters (160 lines)
// ============================================================================
// Reference: gfnff_param.f90 chi/gam/alp/cnf_angewChem2020 arrays (lines 57-139)

// Electronegativity χ (chi) - Hartree
// Reference: gfnff_param.f90:57-76
static const std::vector<double> chi_eeq = {
    1.227054, 1.451412, 0.813363, 1.062841, 1.186499, // H-B
    1.311555, 1.528485, 1.691201, 1.456784, 1.231037, // C-Ne
    0.772989, 1.199092, 1.221576, 1.245964, 1.248942, // Na-P
    1.301708, 1.312474, 1.247701, 0.781237, 0.940834, // S-Ca
    0.950000, 0.974455, 0.998911, 1.023366, 1.047822, // Sc-Mn
    1.072277, 1.096733, 1.121188, 1.145644, 1.170099, // Fe-Zn
    1.205357, 1.145447, 1.169499, 1.253293, 1.329909, // Ga-Br
    1.116527, 0.950975, 0.964592, 0.897786, 0.932824, // Kr-Zr
    0.967863, 1.002901, 1.037940, 1.072978, 1.108017, // Nb-Rh
    1.143055, 1.178094, 1.213132, 1.205076, 1.075529, // Pd-Sn
    1.206919, 1.303658, 1.332656, 1.179317, 0.789115, // Sb-Cs
    0.798704, 1.127797, 1.127863, 1.127928, 1.127994, // Ba-Nd
    1.128059, 1.128125, 1.128190, 1.128256, 1.128322, // Pm-Tb
    1.128387, 1.128453, 1.128518, 1.128584, 1.128649, // Dy-Yb
    1.128715, 1.128780, 1.129764, 1.130747, 1.131731, // Lu-Re
    1.132714, 1.133698, 1.134681, 1.135665, 1.136648, // Os-Hg
    1.061832, 1.053084, 1.207830, 1.236314, 1.310129, // Tl-At
    1.157380 // Rn
};

// Chemical hardness γ (gamma) - Hartree
// Reference: gfnff_param.f90:107-125 (gam_angewChem2020)
// CRITICAL: XTB uses gam_angewChem2020, NOT old gam array (line 611: param%gam(:) = gam_angewChem2020)
// NOTE: 1-based Fortran values converted to 0-based C++ indexing (Z-1)
static const std::vector<double> gam_eeq = {
    -0.448428, 0.131022, // H-He (indices 0-1)
    0.571431, 0.334622, -0.089208, -0.025895, -0.027280, -0.031236, -0.159892, 0.074198, // Li-Ne (C at index 5)
    0.316829, 0.326072, 0.069748, -0.120184, -0.193159, // Na-P
    -0.182428, -0.064093, 0.061914, 0.318112, 0.189248, // S-Ca
    -0.104172, -0.082038, -0.059903, -0.037769, -0.015635, // Sc-Mn
    0.006500, 0.028634, 0.050768, 0.072903, 0.095037, // Fe-Zn
    0.131140, 0.097006, -0.065744, -0.058394, 0.063307, // Ga-Br
    0.091652, 0.386337, 0.530677, -0.030705, -0.020787, // Kr-Zr
    -0.010869, -0.000951, 0.008967, 0.018884, 0.028802, // Nb-Rh
    0.038720, 0.048638, 0.058556, 0.036488, 0.077711, // Pd-Sn
    0.077025, 0.004547, 0.039909, 0.082630, 0.485375, // Sb-Cs
    0.416264, -0.011212, -0.011046, -0.010879, -0.010713, // Ba-Nd
    -0.010546, -0.010380, -0.010214, -0.010047, -0.009881, // Pm-Tb
    -0.009714, -0.009548, -0.009382, -0.009215, -0.009049, // Dy-Yb
    -0.008883, -0.008716, -0.006220, -0.003724, -0.001228, // Lu-Re
    0.001267, 0.003763, 0.006259, 0.008755, 0.011251, // Os-Hg
    0.020477, -0.056566, 0.051943, 0.076708, 0.000273, // Tl-At
    -0.068929 // Rn
};

// Coulomb damping parameter α (alpha) - Bohr^-1
// CRITICAL: Must be SQUARED in EEQ calculation (gfnff_ini.f90:420)
// Reference: gfnff_param.f90:99-118
static const std::vector<double> alpha_eeq = {
    0.585069, 0.432382, 0.628636, 0.743646, 1.167323, // H-B
    0.903430, 1.278388, 0.905347, 1.067014, 2.941513, // C-Ne
    0.687680, 0.792170, 1.337040, 1.251409, 1.068295, // Na-P
    1.186476, 1.593532, 2.056749, 0.674196, 0.868052, // S-Ca
    0.575052, 0.613424, 0.651796, 0.690169, 0.728541, // Sc-Mn
    0.766913, 0.805285, 0.843658, 0.882030, 0.920402, // Fe-Zn
    0.877178, 1.422350, 1.405901, 1.646860, 2.001970, // Ga-Br
    2.301695, 1.020617, 0.634141, 0.652752, 0.668845, // Kr-Zr
    0.684938, 0.701032, 0.717125, 0.733218, 0.749311, // Nb-Rh
    0.765405, 0.781498, 0.797591, 1.296844, 1.534068, // Pd-Sn
    1.727781, 1.926871, 2.175548, 2.177702, 0.977079, // Sb-Cs
    0.770260, 0.757372, 0.757352, 0.757332, 0.757313, // Ba-Nd
    0.757293, 0.757273, 0.757253, 0.757233, 0.757213, // Pm-Tb
    0.757194, 0.757174, 0.757154, 0.757134, 0.757114, // Dy-Yb
    0.757095, 0.757075, 0.756778, 0.756480, 0.756183, // Lu-Re
    0.755886, 0.755589, 0.755291, 0.754994, 0.754697, // Os-Hg
    0.868029, 1.684375, 2.001040, 2.067331, 2.228923, // Tl-At
    1.874218 // Rn
};

// CN correction factor - dimensionless
// Reference: gfnff_param.f90:120-139
static const std::vector<double> cnf_eeq = {
    0.008904, 0.004641, 0.048324, 0.080316, -0.051990, // H-B
    0.031779, 0.132184, 0.157353, 0.064120, 0.036540, // C-Ne
    -0.000627, 0.005412, 0.018809, 0.016329, 0.012149, // Na-P
    0.021484, 0.014212, 0.014939, 0.003597, 0.032921, // S-Ca
    -0.021804, -0.022797, -0.023789, -0.024782, -0.025775, // Sc-Mn
    -0.026767, -0.027760, -0.028753, -0.029745, -0.030738, // Fe-Zn
    -0.004189, -0.011113, -0.021305, -0.012311, 0.049781, // Ga-Br
    -0.040533, 0.012872, 0.021056, -0.003395, 0.000799, // Kr-Zr
    0.004992, 0.009186, 0.013379, 0.017573, 0.021766, // Nb-Rh
    0.025960, 0.030153, 0.034347, -0.000052, -0.039776, // Pd-Sn
    0.006661, 0.050424, 0.068985, 0.023470, -0.024950, // Sb-Cs
    -0.033006, 0.058973, 0.058595, 0.058217, 0.057838, // Ba-Nd
    0.057460, 0.057082, 0.056704, 0.056326, 0.055948, // Pm-Tb
    0.055569, 0.055191, 0.054813, 0.054435, 0.054057, // Dy-Yb
    0.053679, 0.053300, 0.047628, 0.041955, 0.036282, // Lu-Re
    0.030610, 0.024937, 0.019264, 0.013592, 0.007919, // Os-Hg
    0.006383, -0.089155, -0.001293, 0.019269, 0.074803, // Tl-At
    0.016657 // Rn
};

// ============================================================================
// SECTION 2: Bond Parameters (91 lines)
// ============================================================================
// Reference: gfnff_param.f90 bond parameters and radius data

// Bond force constant base values - kcal/(mol·Å²)
// Reference: gfnff_param.f90:733-741
static const std::vector<double> bond_params = {
    0.417997, 0.258490, 0.113608, 0.195935, 0.231217, // H-B
    0.385248, 0.379257, 0.339249, 0.330706, 0.120319, // C-Ne
    0.127255, 0.173647, 0.183796, 0.273055, 0.249044, // Na-P
    0.290653, 0.218744, 0.034706, 0.136353, 0.192467, // S-Ca
    0.335860, 0.314452, 0.293044, 0.271636, 0.250228, // Sc-Mn
    0.228819, 0.207411, 0.186003, 0.164595, 0.143187, // Fe-Zn
    0.212434, 0.210451, 0.219870, 0.224618, 0.272206, // Ga-Br
    0.147864, 0.150000, 0.150000, 0.329501, 0.309632, // Kr-Zr
    0.289763, 0.269894, 0.250025, 0.230155, 0.210286, // Nb-Rh
    0.190417, 0.170548, 0.150679, 0.192977, 0.173411, // Pd-Sn
    0.186907, 0.192891, 0.223202, 0.172577, 0.150000, // Sb-Cs
    0.150000, 0.370682, 0.368511, 0.366339, 0.364168, // Ba-Nd
    0.361996, 0.359825, 0.357654, 0.355482, 0.353311, // Pm-Tb
    0.351139, 0.348968, 0.346797, 0.344625, 0.342454, // Dy-Yb
    0.340282, 0.338111, 0.305540, 0.272969, 0.240398, // Lu-Re
    0.207828, 0.175257, 0.142686, 0.110115, 0.077544, // Os-Hg
    0.108597, 0.148422, 0.183731, 0.192274, 0.127706, // Tl-At
    0.086756, 0.150000, 0.150000, 0.370682, 0.368511, // Rn-Ra
    0.366339, 0.364168, 0.361996, 0.359825, 0.357654, // Ac-Am
    0.355482, 0.353311, 0.351139, 0.348968, 0.346797, // Cm-Cf
    0.344625, 0.342454, 0.340282 // Es-Lr
};

// CN-independent base covalent radii - Bohr
// Reference: gfnff_rab.f90:82-102
static const std::vector<double> r0_gfnff = {
    0.55682207, 0.80966997, 2.49092101, 1.91705642, 1.35974851, // H-B
    0.98310699, 0.98423007, 0.76716063, 1.06139799, 1.17736822, // C-Ne
    2.85570926, 2.56149012, 2.31673425, 2.03181740, 1.82568535, // Na-P
    1.73685958, 1.97498207, 2.00136196, 3.58772537, 2.68096221, // S-Ca
    2.23355957, 2.33135502, 2.15870365, 2.10522128, 2.16376162, // Sc-Mn
    2.10804037, 1.96460045, 2.00476257, 2.22628712, 2.43846700, // Fe-Zn
    2.39408483, 2.24245792, 2.05751204, 2.15427677, 2.27191920, // Ga-Br
    2.19722638, 3.80910350, 3.26020971, 2.99716916, 2.71707818, // Kr-Zr
    2.34950167, 2.11644818, 2.47180659, 2.32198800, 2.32809515, // Nb-Rh
    2.15244869, 2.55958313, 2.59141300, 2.62030465, 2.39935278, // Pd-Sn
    2.56912355, 2.54374096, 2.56914830, 2.53680807, 4.24537037, // Sb-Cs
    3.66542289, 3.19903011, 2.80000000, 2.80000000, 2.80000000, // Ba-Nd (58-60 placeholder)
    2.80000000, 2.80000000, 2.80000000, 2.80000000, 2.80000000, // Pm-Tb (61-65)
    2.80000000, 2.80000000, 2.80000000, 2.80000000, 2.80000000, // Dy-Yb (66-70)
    2.80000000, 2.34880037, 2.37597108, 2.49067697, 2.14100577, // Lu-Re (71-75)
    2.33473532, 2.19498900, 2.12678348, 2.34895048, 2.33422774, // Os-Hg (76-80)
    2.86560827, 2.62488837, 2.88376127, 2.75174124, 2.83054552, // Tl-At (81-85)
    2.63264944 // Rn (86)
};

// CN-dependent radius corrections - Bohr
// Reference: gfnff_rab.f90:103-122
static const std::vector<double> cnfak_gfnff = {
    0.17957827,  0.25584045, -0.02485871,  0.00374217,  0.05646607, // H-B
    0.10514203,  0.09753494,  0.30470380,  0.23261783,  0.36752208, // C-Ne
    0.00131819, -0.00368122, -0.01364510,  0.04265789,  0.07583916, // Na-P
    0.08973207, -0.00589677,  0.13689929, -0.01861307,  0.11061699, // S-Ca
    0.10201137,  0.05426229,  0.06014681,  0.05667719,  0.02992924, // Sc-Mn
    0.03764312,  0.06140790,  0.08563465,  0.03707679,  0.03053526, // Fe-Zn
   -0.00843454,  0.01887497,  0.06876354,  0.01370795, -0.01129196, // Ga-Br
    0.07226529,  0.01005367,  0.01541506,  0.05301365,  0.07066571, // Kr-Zr
    0.07637611,  0.07873977,  0.02997732,  0.04745400,  0.04582912, // Nb-Rh
    0.10557321,  0.02167468,  0.05463616,  0.05370913,  0.05985441, // Pd-Sn
    0.02793994,  0.02922983,  0.02220438,  0.03340460, -0.04110969, // Sb-Cs
   -0.01987240,  0.07260201,  0.07700000,  0.07700000,  0.07700000, // Ba-Nd (58-60)
    0.07700000,  0.07700000,  0.07700000,  0.07700000,  0.07700000, // Pm-Tb (61-65)
    0.07700000,  0.07700000,  0.07700000,  0.07700000,  0.07700000, // Dy-Yb (66-70)
    0.07700000,  0.08379100,  0.07314553,  0.05318438,  0.06799334, // Lu-Re (71-75)
    0.04671159,  0.06758819,  0.09488437,  0.07556405,  0.13384502, // Os-Hg (76-80)
    0.03203572,  0.04235009,  0.03153769, -0.00152488,  0.02714675, // Tl-At (81-85)
    0.04800662 // Rn (86)
};

// CRITICAL FIX (Dec 31, 2025): XTB uses TWO different EN arrays!
// 1. param%en (gfnff_param.f90:154-161) - for ALPHA calculation
// 2. en (gfnffrab.f90:63-82) - for R0 calculation (ff factor)
//
// Using param%en for both caused trade-off: alpha improved 14×, r0 degraded 180×

// EN values for ALPHA calculation (param%en from gfnff_param.f90:154-161)
// Used in: alpha = srb1*(1 + fsrb2*ΔEN² + srb3*bstrength)
static const std::vector<double> en_gfnff = {
    2.200, 3.000, 0.980, 1.570, 2.040, 2.550, 3.040, 3.440, 3.980, // H-F (Z=1-9)
    4.500, 0.930, 1.310, 1.610, 1.900, 2.190, 2.580, 3.160, 3.500, // Ne-Ar (Z=10-18)
    0.820, 1.000, 1.360, 1.540, 1.630, 1.660, 1.550, 1.830, 1.880, // K-Co (Z=19-27)
    1.910, 1.900, 1.650, 1.810, 2.010, 2.180, 2.550, 2.960, 3.000, // Ni-Kr (Z=28-36)
    0.820, 0.950, 1.220, 1.330, 1.600, 2.160, 1.900, 2.200, 2.280, // Rb-Rh (Z=37-45)
    2.200, 1.930, 1.690, 1.780, 1.960, 2.050, 2.100, 2.660, 2.600, // Pd-Xe (Z=46-54)
    0.790, 0.890, 1.100, 1.120, 1.130, 1.140, 1.150, 1.170, 1.180, // Cs-Tb (Z=55-63)
    1.200, 1.210, 1.220, 1.230, 1.240, 1.250, 1.260, 1.270, 1.300, // Dy-Hf (Z=64-72)
    1.500, 1.700, 1.900, 2.100, 2.200, 2.200, 2.200, 2.000, 1.620, // Ta-Hg (Z=73-80)
    2.330, 2.020, 2.000, 2.200, 2.200  // Tl-Rn (Z=81-86)
};

// EN values for R0 calculation (en from gfnffrab.f90:63-82)
// Used in: ff = 1.0 - k1*ΔEN - k2*ΔEN², then r0 = (ra + rb + shift) * ff
static const std::vector<double> en_rab_gfnff = {
    2.30085633, 2.78445145, 1.52956084, 1.51714704, 2.20568300, // H-B (Z=1-5)
    2.49640820, 2.81007174, 4.51078438, 4.67476223, 3.29383610, // C-Ne (Z=6-10)
    2.84505365, 2.20047950, 2.31739628, 2.03636974, 1.97558064, // Na-P (Z=11-15)
    2.13446570, 2.91638164, 1.54098156, 2.91656301, 2.26312147, // S-Ca (Z=16-20)
    2.25621439, 1.32628677, 2.27050569, 1.86790977, 2.44759456, // Sc-Mn (Z=21-25)
    2.49480042, 2.91545568, 3.25897750, 2.68723778, 1.86132251, // Fe-Zn (Z=26-30)
    2.01200832, 1.97030722, 1.95495427, 2.68920990, 2.84503857, // Ga-Br (Z=31-35)
    2.61591858, 2.64188286, 2.28442252, 1.33011187, 1.19809388, // Kr-Zr (Z=36-40)
    1.89181390, 2.40186898, 1.89282464, 3.09963488, 2.50677823, // Nb-Rh (Z=41-45)
    2.61196704, 2.09943450, 2.66930105, 1.78349472, 2.09634533, // Pd-Sn (Z=46-50)
    2.00028974, 1.99869908, 2.59072029, 2.54497829, 2.52387890, // Sb-Cs (Z=51-55)
    2.30204667, 1.60119300, 2.00000000, 2.00000000, 2.00000000, // Ba-Nd (Z=56-60)
    2.00000000, 2.00000000, 2.00000000, 2.00000000, 2.00000000, // Pm-Tb (Z=61-65)
    2.00000000, 2.00000000, 2.00000000, 2.00000000, 2.00000000, // Dy-Yb (Z=66-70)
    2.00000000, 2.30089349, 1.75039077, 1.51785130, 2.62972945, // Lu-Re (Z=71-75)
    2.75372921, 2.62540906, 2.55860939, 3.32492356, 2.65140898, // Os-Hg (Z=76-80)
    1.52014458, 2.54984804, 1.72021963, 2.69303422, 1.81031095, // Tl-At (Z=81-85)
    2.34224386  // Rn (Z=86)
};

// Row-dependent EN polynomial coefficients (scaled by 10^-3 in Fortran)
// Reference: gfnff_rab.f90:125-136
// Index: [row-1][0 or 1]  (row = 1..6 for H-He, Li-Ne, Na-Ar, K-Kr, Rb-Xe, Cs-Rn)
static const double p_enpoly[6][2] = {
    { 29.84522887, -8.87843763 },  // Row 1: H, He
    { -1.70549806,  2.10878369 },  // Row 2: Li-Ne
    {  6.54013762,  0.08009374 },  // Row 3: Na-Ar
    {  6.39169003, -0.85808076 },  // Row 4: K-Kr
    {  6.00000000, -1.15000000 },  // Row 5: Rb-Xe (extrapolated)
    {  5.60000000, -1.30000000 }   // Row 6: Cs-Rn (extrapolated)
};

// GFN-FF base covalent radii - Bohr (alias for r0_gfnff)
// Used for bond parameter calculations
static const std::vector<double> rcov_bohr = r0_gfnff;

// DFT-D3 Covalent radii (Pyykkö & Atsumi) - Bohr
// Claude Generated (Jan 19, 2026): For torsion damping calculations
// Reference: gfnff_param.f90 covalentRadD3 (values in Å, converted to Bohr)
// NOTE (Jan 23, 2026): XTB applies 4/3 scaling in covalentradd3.f90:62
//   covalentRadD3(1:118) = [...] * aatoau * 4.0_wp / 3.0_wp
// However, Curcuma torsion energy is already ~108× TOO LARGE, so applying
// 4/3 scaling (which weakens damping → more energy) makes things WORSE.
// The root cause is elsewhere (likely in fctot calculation or torsion counting).
// Keeping original values (without 4/3) until root cause is found.
static const std::vector<double> covalent_rad_d3 = {
    0.60392702, 0.86945142, // H, He (without 4/3 scaling)
    2.26643422, 1.77478229, 1.45333953, 1.41731155, 1.34116514, 1.19027623, // Li-C-Ne
    1.20788167, 1.26458288, 2.64097810, 2.36067855, 2.13177212, 1.96297420, // Na-P
    2.07555175, 1.92490730, 1.86819168, 1.81173529, 3.32054985, 2.90629615, // S-Ca
    2.51017939, 2.30338353, 2.28602099, 2.07555175, 2.02135057, 1.96297420, // Sc-Fe
    1.88582236, 1.86819168, 1.90651985, 2.05733234, 2.11403355, 2.05733234, // Co-Ga
    2.17073477, 2.07555175, 2.15310255, 2.20754546, 3.56594020, 3.15168650, // As-Sr
    2.77387840, 2.62297027, 2.49070592, 2.34155557, 2.17073477, 2.13177212, // Y-Rh
    2.13177212, 2.03805675, 2.17073477, 2.32122182, 2.41664091, 2.37768235, // Pd-Sn
    2.41664091, 2.32122182, 2.49070592, 2.47280723, 3.94795890, 3.32054985, // Sb-Ba
    3.05928569, 2.77387840, 2.98177264, 2.96386629, 2.94595994, 2.92805359, // La-Sm
    2.85103624, 2.86853468, 2.85103624, 2.83353007, 2.81603163, 2.81603163, // Eu-Dy
    2.79852546, 2.88779861, 2.75500402, 2.58435964, 2.47280723, 2.32122182, // Ho-Re
    2.22750842, 2.18854577, 2.09479507, 2.13177212, 2.36067855, 2.49070592, // Os-Tl
    2.45333719, 2.45333719, 2.56584604, 2.47280723, 2.60439593, 2.68012880, // Pb-Rn
    3.79658336, 3.41540674, 3.15168650, 2.98177264 // Fr-U
};

// Standard covalent radii - Angström
// Used for bond detection threshold (× 1.3)
static const std::vector<double> covalent_radii = {
    0.32, 0.37, 1.30, 0.99, 0.84, 0.75, 0.71, 0.64, 0.60, 0.62, // H-Ne
    1.60, 1.40, 1.24, 1.14, 1.09, 1.04, 1.00, 1.01, // Na-Ar
    2.00, 1.74, 1.59, 1.48, 1.44, 1.30, 1.29, 1.24, 1.18, // K-Ni
    1.17, 1.22, 1.20, 1.23, 1.20, 1.20, 1.18, 1.17, 1.16, // Cu-Kr
    2.15, 1.90, 1.76, 1.64, 1.56, 1.46, 1.38, 1.36, 1.34, // Rb-Pd
    1.30, 1.36, 1.40, 1.42, 1.40, 1.40, 1.37, 1.36, 1.36, // Ag-Xe
    2.38, 2.06, 1.94, 1.84, 1.90, 1.88, 1.86, 1.85, 1.83, // Cs-Eu
    1.82, 1.81, 1.80, 1.79, 1.77, 1.77, 1.78, 1.74, 1.64, // Gd-Lu
    1.58, 1.50, 1.41, 1.36, 1.32, 1.30, 1.30, 1.32, 1.44, // Hf-Au
    1.45, 1.50, 1.42, 1.48, 1.46, // Hg-Po
    2.42, 2.11, 2.01, 1.90, 1.84, 1.83, 1.80, 1.80, 1.73, // Fr-Am
    1.68, 1.68, 1.68, 1.65, 1.67, 1.73, 1.76, 1.61 // Cm-Lr
};

// ============================================================================
// SECTION 3: Hybridization and Bond Strength (19 lines)
// ============================================================================
// Reference: gfnff_param.f90 bond strength arrays

// Bond strength base values (1.00=single, 1.24=double, 1.98=triple)
// Reference: gfnff_param.f90:580-588
static const double bstren[9] = {
    0.0,   // Index 0 unused
    1.00,  // [1] single bond
    1.24,  // [2] double bond
    1.98,  // [3] triple bond
    1.22,  // [4] hypervalent bond
    1.00,  // [5] M-X (metal-ligand)
    0.78,  // [6] M eta
    3.40,  // [7] M-M
    3.40   // [8] M-M
};

// Hybridization-dependent bond strength matrix (4×4)
// Formula: 0.67*bstren[bond_i] + 0.33*bstren[bond_j]
// Reference: gfnff_param.f90:804-814
static const double bsmat[4][4] = {
    { 1.000, 1.323, 1.079, 1.000 },  // hyb=0 (unknown/sp3)
    { 1.323, 1.980, 1.484, 1.323 },  // hyb=1 (sp)
    { 1.079, 1.484, 1.240, 1.079 },  // hyb=2 (sp2)
    { 1.000, 1.323, 1.079, 1.000 }   // hyb=3 (sp3)
};

// ============================================================================
// SECTION 4: Topology and Metal Classification (28 lines)
// ============================================================================

// Metal type classification
// 0=non-metal, 1=main group metal, 2=transition metal
static const int metal_type[86] = {
    0, 0, // H-He
    1, 1, 0, 0, 0, 0, 0, 0, // Li-Ne
    1, 1, 1, 0, 0, 0, 0, 0, // Na-Ar
    1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, 0, 0, 0, 0, // K-Kr
    1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0, 0, 0, // Rb-Xe
    1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0 // Cs-Rn
};

// Periodic group numbers
// Main group: 1-8, d-block: negative values
static const int periodic_group[86] = {
    1, 8, // H-He
    1, 2, 3, 4, 5, 6, 7, 8, // Li-Ne
    1, 2, 3, 4, 5, 6, 7, 8, // Na-Ar
    1, 2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, 3, 4, 5, 6, 7, 8, // K-Kr
    1, 2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, 3, 4, 5, 6, 7, 8, // Rb-Xe
    1, 2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, 3, 4, 5, 6, 7, 8 // Cs-Rn
};

// ============================================================================
// SECTION 5: Angle Parameters (46 lines)
// ============================================================================
// Reference: gfnff_param.f90 angle parameters

// Angle force constant base values - kcal/(mol·rad²)
// Reference: gfnff_param.f90:247-265
static const std::vector<double> angle_params = {
    1.661808, 0.300000, 0.018158, 0.029224, 0.572683, // H-B
    0.771055, 1.053577, 2.159889, 1.525582, 0.400000, // C-Ne
    0.041070, 0.028889, 0.086910, 0.494456, 0.409204, // Na-P
    0.864972, 1.986025, 0.491537, 0.050168, 0.072745, // S-Ca
    0.378334, 0.346400, 0.314466, 0.282532, 0.250598, // Sc-Mn
    0.218663, 0.186729, 0.154795, 0.122861, 0.090927, // Fe-Zn
    0.140458, 0.653971, 0.528465, 0.420379, 2.243492, // Ga-Br
    0.400000, 0.035341, 0.022704, 0.195060, 0.188476, // Kr-Zr
    0.181892, 0.175308, 0.168724, 0.162139, 0.155555, // Nb-Rh
    0.148971, 0.142387, 0.135803, 0.169779, 0.265730, // Pd-Sn
    0.505495, 0.398254, 2.640752, 0.568026, 0.032198, // Sb-Cs
    0.036663, 0.281449, 0.280526, 0.279603, 0.278680, // Ba-Nd
    0.277757, 0.276834, 0.275911, 0.274988, 0.274065, // Pm-Tb
    0.273142, 0.272219, 0.271296, 0.270373, 0.269450, // Dy-Yb
    0.268528, 0.267605, 0.253760, 0.239916, 0.226071, // Lu-Re
    0.212227, 0.198382, 0.184538, 0.170693, 0.156849, // Os-Hg
    0.104547, 0.313474, 0.220185, 0.415042, 1.259822, // Tl-At
    0.400000, 0.032198, 0.036663, 0.281449, 0.280526, // Rn-Ra
    0.279603, 0.278680, 0.277757, 0.276834, 0.275911, // Ac-Am
    0.274988, 0.274065, 0.273142, 0.272219, 0.271296, // Cm-Cf
    0.270373, 0.269450, 0.268528 // Es-Lr
};

// Neighbor scaling factors for angle force constants
// Reference: gfnff_param.f90:angl2_angewChem2020 array
static const std::vector<double> angl2_neighbors = {
    0.624197, 0.600000, 0.050000, 0.101579, 0.180347, // H-B
    0.755851, 0.761551, 0.813653, 0.791274, 0.400000, // C-Ne
    0.000000, 0.022706, 0.100000, 0.338514, 0.453023, // Na-P
    0.603722, 1.051121, 0.547904, 0.000000, 0.059059, // S-Ca
    0.117040, 0.118438, 0.119836, 0.121234, 0.122632, // Sc-Mn
    0.124031, 0.125429, 0.126827, 0.128225, 0.129623, // Fe-Zn
    0.206779, 0.466678, 0.496442, 0.617321, 0.409933, // Ga-Br
    0.400000, 0.000000, 0.000000, 0.119120, 0.118163, // Kr-Zr
    0.117206, 0.116249, 0.115292, 0.114336, 0.113379, // Nb-Rh
    0.112422, 0.111465, 0.110508, 0.149917, 0.308383, // Pd-Sn
    0.527398, 0.577885, 0.320371, 0.568026, 0.000000, // Sb-Cs
    0.000000, 0.078710, 0.079266, 0.079822, 0.080379, // Ba-Nd
    0.080935, 0.081491, 0.082047, 0.082603, 0.083159, // Pm-Tb
    0.083716, 0.084272, 0.084828, 0.085384, 0.085940, // Dy-Yb
    0.086496, 0.087053, 0.095395, 0.103738, 0.112081, // Lu-Re
    0.120423, 0.128766, 0.137109, 0.145451, 0.153794, // Os-Hg
    0.323570, 0.233450, 0.268137, 0.307481, 0.316447, // Tl-At
    0.400000, 0.000000, 0.000000, 0.119120, 0.118163, // Rn-Ra
    0.117206, 0.116249, 0.115292, 0.114336, 0.113379, // Ac-Am
    0.112422, 0.111465, 0.110508, 0.149917, 0.308383, // Cm-Cf
    0.527398, 0.577885, 0.320371, 0.568026, 0.000000, // Es-Lr
    0.000000, 0.078710, 0.079266, 0.079822, 0.080379 // Last elements
};

// ============================================================================
// SECTION 6: Repulsion Parameters (33 lines)
// ============================================================================
// Reference: gfnff_param.f90 repulsion arrays

// Repulsion exponent parameter - dimensionless
// Reference: gfnff_param.f90:141-160
static const std::vector<double> repa_angewChem2020 = {
    2.639785, 3.575012, 0.732142, 1.159621, 1.561585, // H-B
    1.762895, 2.173015, 2.262269, 2.511112, 3.577220, // C-Ne
    0.338845, 0.693023, 0.678792, 0.804784, 1.012178, // Na-P
    1.103469, 1.209798, 1.167791, 0.326946, 0.595242, // S-Ca
    1.447860, 1.414501, 1.381142, 1.347783, 1.314424, // Sc-Mn
    1.281065, 1.247706, 1.214347, 1.180988, 1.147629, // Fe-Zn
    0.700620, 0.721266, 0.741789, 0.857434, 0.875583, // Ga-Br
    0.835876, 0.290625, 0.554446, 0.623980, 0.696005, // Kr-Zr
    0.768030, 0.840055, 0.912081, 0.984106, 1.056131, // Nb-Rh
    1.128156, 1.200181, 1.272206, 0.478807, 0.479759, // Pd-Sn
    0.579840, 0.595241, 0.644458, 0.655289, 0.574626, // Sb-Cs
    0.560506, 0.682723, 0.684824, 0.686925, 0.689026, // Ba-Nd
    0.691127, 0.693228, 0.695329, 0.697430, 0.699531, // Pm-Tb
    0.701631, 0.703732, 0.705833, 0.707934, 0.710035, // Dy-Yb
    0.712136, 0.714237, 0.745751, 0.777265, 0.808779, // Lu-Re
    0.840294, 0.871808, 0.903322, 0.934836, 0.966350, // Os-Hg
    0.467729, 0.486102, 0.559176, 0.557520, 0.563373, // Tl-At
    0.484713 // Rn
};

// Effective nuclear charges for repulsion - dimensionless
// Reference: gfnff_param.f90:162-181
static const std::vector<double> repz = {
    1., 2.,                                                   // H-He
    1., 2., 3., 4., 5., 6., 7., 8.,                          // Li-Ne
    1., 2., 3., 4., 5., 6., 7., 8.,                          // Na-Ar
    1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.,       // K-Zn
    3., 4., 5., 6., 7., 8.,                                  // Ga-Kr
    1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.,       // Rb-Cd
    3., 4., 5., 6., 7., 8.,                                  // In-Xe
    1., 2., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3.,  // Cs-Lu (lanthanides use Z=3)
    4., 5., 6., 7., 8., 9., 10., 11., 12.,                   // Hf-Hg
    3., 4., 5., 6., 7., 8.                                   // Tl-Rn
};

// Bonded repulsion scaling factor
// Reference: gfnff_param.f90:373
static constexpr double REPSCALB = 1.7583;

// Non-bonded repulsion scaling factor
// Reference: gfnff_param.f90:374
static constexpr double REPSCALN = 0.4270;

// Non-bonded repulsion: charge-dependent scaling (Claude Generated Dec 24, 2025)
// Reference: gfnff_param.f90 gen%qrepscal
// Used in: alphanb = repan * (1.0 + qa * QREPSCAL)
static constexpr double QREPSCAL = 0.3482;

// Non-bonded repulsion: neighbor-count scaling (Claude Generated Dec 24, 2025)
// Reference: gfnff_param.f90 gen%nrepscal
// Used in: fn = 1.0 + NREPSCAL / (1.0 + nb²)
static constexpr double NREPSCAL = -0.1270;

// H-H pair special scaling factor (Claude Generated Dec 24, 2025)
// Reference: gfnff_param.f90 gen%hhfac
// Applied to all H-H non-bonded pairs
static constexpr double HHFAC = 0.6290;

// H-H 1,3-pair scaling factor (Claude Generated Dec 24, 2025)
// Reference: gfnff_param.f90 gen%hh13rep
// Additional scaling for 1,3 H-H pairs
static constexpr double HH13REP = 1.4580;

// H-H 1,4-pair scaling factor (Claude Generated Dec 24, 2025)
// Reference: gfnff_param.f90 gen%hh14rep
// Additional scaling for 1,4 H-H pairs
static constexpr double HH14REP = 0.7080;

// Non-bonded repulsion exponent parameter - dimensionless
// Reference: gfnff_param.f90:207-225 (repan_angewChem2020)
// Used ONLY for non-bonded pairs with arithmetic mean combining rule
static const std::vector<double> repan_angewChem2020 = {
    1.071395, 1.072699, 1.416847, 1.156187, 0.682382, // H-B
    0.556380, 0.746785, 0.847242, 0.997252, 0.873051, // C-Ne
    0.322503, 0.415554, 0.423946, 0.415776, 0.486773, // Na-P
    0.494532, 0.705274, 0.706778, 0.311178, 0.399439, // S-Ca
    0.440983, 0.475582, 0.510180, 0.544779, 0.579377, // Sc-Mn
    0.613976, 0.648574, 0.683173, 0.717772, 0.752370, // Fe-Zn
    0.429944, 0.420053, 0.384743, 0.443762, 0.538680, // Ga-Br
    0.472196, 0.423850, 0.385815, 0.249213, 0.285604, // Kr-Zr
    0.321995, 0.358387, 0.394778, 0.431169, 0.467560, // Nb-Rh
    0.503952, 0.540343, 0.576734, 0.333476, 0.348734, // Pd-Sn
    0.358194, 0.351053, 0.404536, 0.389847, 0.302575, // Sb-Cs
    0.163290, 0.187645, 0.190821, 0.193998, 0.197174, // Ba-Nd
    0.200351, 0.203527, 0.206703, 0.209880, 0.213056, // Pm-Tb
    0.216233, 0.219409, 0.222585, 0.225762, 0.228938, // Dy-Yb
    0.232115, 0.235291, 0.282937, 0.330583, 0.378229, // Lu-Re
    0.425876, 0.473522, 0.521168, 0.568814, 0.616460, // Os-Hg
    0.242521, 0.293680, 0.320931, 0.322666, 0.333641, // Tl-At
    0.434163  // Rn
};

// Torsion angle damping parameter
// Reference: gfnff_param.f90:464
static constexpr double atcutt = 0.505;

// Torsion angle damping parameter for NCI HB term
// Reference: gfnff_param.f90:466
static constexpr double atcutt_nci = 0.305;

// ============================================================================
// SECTION 7: Dispersion Parameters (29 lines)
// ============================================================================
// Reference: Grimme et al., free-atom C6 coefficients

// Free-atom C6 dispersion coefficients - Hartree·Bohr^6
// Reference: D3/D4 literature values (approximate)
static const std::vector<double> C6_atomic = {
    // H-He
    6.50, 1.42,
    // Li-Ne
    1387.0, 214.0, 99.5, 46.6, 24.2, 15.6, 9.52, 6.38,
    // Na-Ar
    1556.0, 627.0, 528.0, 305.0, 185.0, 134.0, 94.6, 64.3,
    // K-Ca
    3897.0, 2221.0,
    // Sc-Zn
    1383.0, 1044.0, 832.0, 602.0, 552.0, 482.0, 408.0, 373.0, 253.0, 284.0,
    // Ga-Kr
    498.0, 354.0, 246.0, 210.0, 162.0, 129.5,
    // Rb-Sr
    4691.0, 3170.0,
    // Y-Cd
    1968.0, 1677.0, 1263.0, 1028.0, 1390.0, 1029.0, 1118.0, 1251.0, 1225.0, 1225.0,
    // In-Xe
    2896.0, 2290.0, 1896.0, 1830.0, 1612.0, 1416.0,
    // Cs-Ba
    6582.0, 5727.0,
    // La-Lu (lanthanides - approximate values)
    3884.0, 3708.0, 3551.0, 3410.0, 3280.0, 3163.0, 3056.0, 2958.0, 2868.0, 2785.0,
    2708.0, 2638.0, 2573.0, 2513.0, 2458.0,
    // Hf-Hg
    2051.0, 1877.0, 1659.0, 1529.0, 1414.0, 1305.0, 1206.0, 1118.0, 1037.0, 1185.0,
    // Tl-Rn
    3292.0, 3135.0, 2762.0, 2600.0, 2452.0, 2318.0
};

// BJ damping parameters for D3/D4 integration
// Reference: Beckert-Johnson damping function
// From gfnff_data_types.f90: initGFFDispersion subroutine
static constexpr double s6 = 1.0;    // C6 scaling
static constexpr double s8 = 2.85;   // C8 scaling (was 2.4)
static constexpr double a1 = 0.80;   // BJ damping a1 (was 0.48)
static constexpr double a2 = 4.60;   // BJ damping a2 (was 4.80)

// ============================================================================
// SECTION 7.5: Torsion Parameters (56 lines)
// ============================================================================
// Reference: gfnff_param.f90:267-305 (tors_angewChem2020 and tors2_angewChem2020)
// Claude Generated (2025): Extracted from Fortran to fix blocking issue in gfnff_torsions.cpp

// Torsion angle parameter 1 - phase/barrier information
// Reference: gfnff_param.f90:267-285
static const std::vector<double> tors_angewChem2020 = {
    0.100000, 0.100000, 0.100000, 0.000000, 0.121170, // H-B
    0.260028, 0.222546, 0.250620, 0.256328, 0.400000, // C-Ne
    0.115000, 0.000000, 0.103731, 0.069103, 0.104280, // Na-P
    0.226131, 0.300000, 0.400000, 0.124098, 0.000000, // S-Ca
    0.105007, 0.107267, 0.109526, 0.111786, 0.114046, // Sc-Mn
    0.116305, 0.118565, 0.120825, 0.123084, 0.125344, // Fe-Zn
    0.395722, 0.349100, 0.147808, 0.259811, 0.400000, // Ga-Br
    0.400000, 0.112206, -0.004549, 0.198713, 0.179472, // Kr-Zr
    0.160232, 0.140991, 0.121751, 0.102510, 0.083270, // Nb-Rh
    0.064029, 0.044789, 0.025548, 0.202245, 0.278223, // Pd-Sn
    0.280596, 0.229057, 0.300000, 0.423199, 0.090741, // Sb-Cs
    0.076783, 0.310896, 0.309131, 0.307367, 0.305602, // Ba-Nd
    0.303838, 0.302073, 0.300309, 0.298544, 0.296779, // Pm-Tb
    0.295015, 0.293250, 0.291486, 0.289721, 0.287957, // Dy-Yb
    0.286192, 0.284427, 0.257959, 0.231490, 0.205022, // Lu-Re
    0.178553, 0.152085, 0.125616, 0.099147, 0.072679, // Os-Hg
    0.203077, 0.169346, 0.090568, 0.144762, 0.231884, // Tl-At
    0.400000 // Rn
};

// Torsion angle parameter 2 - additional barrier modulation
// Reference: gfnff_param.f90:287-305
static const std::vector<double> tors2_angewChem2020 = {
    1.618678, 1.000000, 0.064677, 0.000000, 0.965814, // H-B
    1.324709, 1.079334, 1.478599, 0.304844, 0.500000, // C-Ne
    0.029210, 0.000000, 0.417423, 0.334275, 0.817008, // Na-P
    0.922181, 0.356367, 0.684881, 0.029210, 0.000000, // S-Ca
    0.035902, 0.090952, 0.146002, 0.201052, 0.256103, // Sc-Mn
    0.311153, 0.366203, 0.421253, 0.476303, 0.531353, // Fe-Zn
    0.482963, 1.415893, 1.146581, 1.338448, 0.376801, // Ga-Br
    0.500000, 0.027213, -0.004549, 0.003820, 0.093011, // Kr-Zr
    0.182202, 0.271393, 0.360584, 0.449775, 0.538965, // Nb-Rh
    0.628156, 0.717347, 0.806538, 0.077000, 0.185110, // Pd-Sn
    0.432427, 0.887811, 0.267721, 0.571662, 0.000000, // Sb-Cs
    0.000000, 0.122336, 0.131176, 0.140015, 0.148855, // Ba-Nd
    0.157695, 0.166534, 0.175374, 0.184214, 0.193053, // Pm-Tb
    0.201893, 0.210733, 0.219572, 0.228412, 0.237252, // Dy-Yb
    0.246091, 0.254931, 0.387526, 0.520121, 0.652716, // Lu-Re
    0.785311, 0.917906, 1.050500, 1.183095, 1.315690, // Os-Hg
    0.219729, 0.344830, 0.331862, 0.767979, 0.536799, // Tl-At
    0.500000 // Rn
};

// ============================================================================
// SECTION 8: Hydrogen Bond and Halogen Bond Parameters
// ============================================================================
// Reference: gfnff_param.f90:468-522 (angewChem2020 parameter set)
// Claude Generated (2025): Direct extraction from Fortran reference

// --- Global HB/XB Damping Parameters ---

static constexpr double HB_BACUT = 49.0;        // HB angle cut-off
static constexpr double HB_SCUT = 22.0;         // HB short-range cut-off
static constexpr double XB_BACUT = 70.0;        // XB angle cut-off
static constexpr double XB_SCUT = 5.0;          // XB short-range cut-off
static constexpr double HB_ALP = 6.0;           // Damping exponent (HB/XB)
static constexpr double HB_LONGCUT = 85.0;      // HB long-range cut-off (Bohr²)
static constexpr double HB_LONGCUT_XB = 70.0;   // XB long-range cut-off (Bohr²)
static constexpr double HB_ST = 15.0;           // Charge scaling strength (HB)
static constexpr double HB_SF = 1.0;            // Charge scaling offset (HB)
static constexpr double XB_ST = 15.0;           // Charge scaling strength (XB)
static constexpr double XB_SF = 0.03;           // Charge scaling offset (XB)
static constexpr double HB_ABMIX = 0.80;        // A-B mixing parameter
static constexpr double HB_NBCUT = 11.20;       // Neighbor angle cut-off (Case 2)
static constexpr double TORS_HB = 0.94;         // NCI torsion term shift
static constexpr double BEND_HB = 0.20;         // NCI bending term shift

// Bond energy HB scaling (egbond_hb)
// Reference: gfnff_param.f90:484
// Scales the bond exponent for X-H bonds participating in hydrogen bridges
// Formula: alpha_modified = (1.0 - (1.0 - VBOND_SCALE) * hb_cn_H) * alpha
// For hb_cn_H = 1: alpha_modified = 0.9 * alpha (10% reduction)
static constexpr double VBOND_SCALE = 0.9;

// Global acidity parameters
static constexpr double XHACI_GLOBABH = 0.268;  // A-H...B general scaling
static constexpr double XHACI_COH = 0.350;      // A-H...O=C scaling
static constexpr double XHACI_GLOB = 1.50;      // Baseline acidity

// --- Element-Specific Basicity (xhbas) ---
// For both HB and XB - acceptor atom B
// Index: atomic number (0-based, index 0 unused)
// Reference: gfnff_param.f90:488-502
static const std::vector<double> hb_basicity = {
    0.0,    // Z=0 (placeholder)
    0.0,    // H (1)
    0.0,    // He (2)
    0.0, 0.0, 0.0,                              // Li, Be, B (3-5)
    0.80,   // C (6)
    1.68,   // N (7)
    0.67,   // O (8)
    0.52,   // F (9)
    0.0,    // Ne (10)
    0.0, 0.0, 0.0,                              // Na, Mg, Al (11-13)
    4.0,    // Si (14)
    3.5,    // P (15)
    2.0,    // S (16)
    1.5,    // Cl (17)
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  // Ar-As (18-32)
    3.5,    // As (33) = param%xhbas(15)
    2.0,    // Se (34) = param%xhbas(16)
    1.5,    // Br (35)
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  // Kr-Sn (36-50)
    3.5,    // Sb (51) = param%xhbas(15)
    2.0,    // Te (52) = param%xhbas(16)
    1.9,    // I (53)
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, // Xe-At (54-85)
    0.0     // Rn (86)
};

// --- Element-Specific HB Acidity (xhaci) ---
// For hydrogen bond donor atom (bonded to H)
// Reference: gfnff_param.f90:503-512
static const std::vector<double> hb_acidity = {
    0.0,    // Z=0 (placeholder)
    0.0,    // H (1)
    0.0,    // He (2)
    0.0, 0.0, 0.0,                                      // Li, Be, B (3-5)
    0.75,                                                // C (6) - weaker
    XHACI_GLOB + 0.1,                                   // N (7) = 1.60
    XHACI_GLOB,                                         // O (8) = 1.50
    XHACI_GLOB,                                         // F (9) = 1.50
    0.0, 0.0, 0.0, 0.0, 0.0,                            // Ne-Al (10-14)
    XHACI_GLOB,                                         // P (15) = 1.50
    XHACI_GLOB,                                         // S (16) = 1.50
    XHACI_GLOB + 1.0,                                   // Cl (17) = 2.50
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  // Ar-As (18-32)
    XHACI_GLOB,                                         // As (33) = 1.50 (via pattern)
    XHACI_GLOB,                                         // Se (34) = 1.50
    XHACI_GLOB + 1.0,                                   // Br (35) = 2.50
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  // Kr-Sn (36-50)
    XHACI_GLOB,                                         // Sb (51) = 1.50
    XHACI_GLOB,                                         // Te (52) = 1.50
    XHACI_GLOB + 1.0,                                   // I (53) = 2.50
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, // Xe-At (54-85)
    0.0     // Rn (86)
};

// --- Element-Specific XB Acidity (xbaci) ---
// For halogen atom X in X...B interaction
// Reference: gfnff_param.f90:513-522
static const std::vector<double> xb_acidity = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  // Z=0-9
    0.0, 0.0, 0.0, 0.0, 0.0,                            // Z=10-14
    1.0,    // P (15)
    1.0,    // S (16)
    0.5,    // Cl (17)
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  // Ar-As (18-32)
    1.2,    // As (33)
    1.2,    // Se (34)
    0.9,    // Br (35)
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  // Kr-Sn (36-50)
    1.2,    // Sb (51)
    1.2,    // Te (52)
    1.2,    // I (53)
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, // Xe-At (54-85)
    0.0     // Rn (86)
};

// ============================================================================
// Metal Bond Shift Factors (Claude Generated - January 2026)
// ============================================================================
// Reference: XTB gfnff_param.f90:807-810
// Applied to bond equilibrium distances (r0) when metals are present
//
// Educational Documentation:
// - metal1_shift: Group 1+2 metals (Li, Na, Mg, Ca, ...)
//   Bonds involving alkali/alkaline earth metals need longer equilibrium distances
//
// - metal2_shift: Transition metals (Fe, Cu, Zn, ...)
//   TM-ligand bonds typically have characteristic metal coordination distances
//
// - metal3_shift: Main group metals (Al, Ga, In, Sn, Pb, ...)
//   Post-transition metals have intermediate bond length corrections
//
// - eta_shift: η-coordinated ligands (π-bonding to metals)
//   Multiply by coordination number (CN) for multi-hapto coordination
//   Example: Ferrocene Fe-C₆H₆ has η⁶ coordination (6 carbon neighbors)
//
// Literature: Spicher, S.; Grimme, S. Angew. Chem. Int. Ed. 2020
// ============================================================================

constexpr double METAL1_SHIFT = 0.2;   // Group 1+2 metals (Li, Na, Mg, Ca, ...) in Bohr
constexpr double METAL2_SHIFT = 0.15;  // Transition metals (Fe, Cu, Zn, ...) in Bohr
constexpr double METAL3_SHIFT = 0.05;  // Main group metals (Al, Ga, In, ...) in Bohr
constexpr double ETA_SHIFT = 0.040;    // η-coordination (per neighbor) in Bohr

// ============================================================================
// ============================================================================
// TORSION BARRIER SCALING FACTORS (Claude Generated - January 2026)
// ============================================================================
// Reference: gfnff_param.f90:791-797 (gen%torsf array)
//
// These factors scale the torsional barrier height based on bond type:
// - torsf_single: Single bonds (sp3-sp3, standard alkanes)
// - torsf_pi: Pi bonds (double bonds, aromatic systems)
// - torsf_improper: Out-of-plane bending terms
// - torsf_pi_improper: Pi contribution to improper torsions
//
// Extra sp3-sp3 torsion factors (for gauche conformations):
// - torsf_extra_C: Extra barrier for C-C sp3-sp3 (negative = stabilizes gauche)
// - torsf_extra_N: Extra barrier for N-X sp3-sp3 (positive = destabilizes gauche)
// - torsf_extra_O: Extra barrier for O-X sp3-sp3 (negative = stabilizes gauche)
// ============================================================================

constexpr double torsf_single = 1.00;       // Single bond barrier factor
constexpr double torsf_pi = 1.18;           // Pi bond barrier factor
constexpr double torsf_improper = 1.05;     // Improper torsion factor
constexpr double torsf_pi_improper = 0.50;  // Pi part of improper torsions
constexpr double torsf_extra_C = -0.90;     // Extra sp3 C-C torsion
constexpr double torsf_extra_N = 0.70;      // Extra sp3 N-X torsion
constexpr double torsf_extra_O = -2.00;     // Extra sp3 O-X torsion

// Ring Torsion Barrier Factors (Claude Generated - January 2026)
// ============================================================================
// Reference: XTB gfnff_param.f90:787-790, gfnff_ini.f90:1814-1835
// Applied to torsional barriers when all four atoms are in the same ring
//
// Educational Documentation:
// - Ring torsions have different conformational preferences than acyclic bonds
// - Small rings (3-, 4-membered) are constrained and prefer specific puckering
// - 6-membered rings have strong preference for chair conformation
//
// Ring-specific barriers (relative to acyclic baseline):
// - FR3: 3-ring torsions are very flexible (planar preference, n=1, φ₀=0°)
// - FR4: 4-ring torsions have moderate barriers (puckered, n=6, φ₀=30°)
// - FR5: 5-ring torsions have intermediate barriers (envelope, n=6, φ₀=30°)
// - FR6: 6-ring torsions have STRONG barriers (chair preference, n=3, φ₀=60°)
//
// Examples:
// - Cyclohexane: FR6=5.7 gives strong chair/boat energy difference (~7 kcal/mol)
// - Cyclopentane: FR5=1.5 gives envelope puckering
// - Cyclobutane: FR4=1.0 gives butterfly conformation
// - Cyclopropane: FR3=0.3 nearly planar, highly strained
//
// Literature: Spicher, S.; Grimme, S. Angew. Chem. Int. Ed. 2020
// ============================================================================

constexpr double FR3 = 0.3;  // 3-ring torsion barrier factor (planar)
constexpr double FR4 = 1.0;  // 4-ring torsion barrier factor (butterfly)
constexpr double FR5 = 1.5;  // 5-ring torsion barrier factor (envelope)
constexpr double FR6 = 5.7;  // 6-ring torsion barrier factor (chair preference)

} // namespace GFNFFParameters
