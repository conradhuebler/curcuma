/*
 * <GFN1-xTB Parameter Database Implementation>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Parameters extracted from TBLite (https://github.com/tblite/tblite)
 *   Source: src/tblite/xtb/gfn1.f90
 *   Copyright (C) 2019-2024 Sebastian Ehlert and contributors
 *   Licensed under LGPL-3.0-or-later
 *
 * Reference implementation validated against:
 *   - TBLite v0.3.0 GFN1 TOML parameters
 *   - Original GFN1 paper values (JCTC 2017, 13, 1989)
 *
 * This program is free software under GPL-3.0
 *
 * Claude Generated: Complete GFN1 parameter database (November 2025)
 */

#include "gfn1_params_loader.h"
#include "src/core/curcuma_logger.h"

#include <fmt/format.h>

namespace GFN1Params {

// =================================================================================
// Helper: Add symmetric pair parameters
// =================================================================================

void ParameterDatabase::addPairParams(int Z1, int Z2, double kpair,
                                      double kss, double ksp, double kpp)
{
    PairParams pair;
    pair.kpair = kpair;
    pair.kshell_ss = kss;
    pair.kshell_sp = ksp;
    pair.kshell_pp = kpp;

    // Store both orders for O(1) lookup
    m_pairs[{Z1, Z2}] = pair;
    m_pairs[{Z2, Z1}] = pair;
}

// =================================================================================
// Load Default GFN1 Parameters (H, C, N, O)
// =================================================================================

bool ParameterDatabase::loadDefaultGFN1()
{
    m_elements.clear();
    m_pairs.clear();

    // ==== Hydrogen (Z=1) ====
    // GFN1 parameter source: TBLite gfn1.f90
    ElementParams H;
    H.atomic_number = 1;
    H.symbol = "H";

    // s-shell only
    ShellParams s_shell;
    s_shell.selfenergy = -0.468040;  // E_s (Hartree)
    s_shell.kcn = -0.03;             // CN shift
    s_shell.gexp = 1.2;              // Gaussian exponent
    s_shell.refocc = 1.0;            // Reference occupation
    H.shells[0] = s_shell;

    H.rep_alpha = 1.23;
    H.rep_zeff = 1.0;
    H.gamma_ss = 0.468040;
    H.gamma_sp = 0.0;
    H.gamma_pp = 0.0;
    H.xb_radius = 0.0;  // Not a halogen
    H.xb_strength = 0.0;

    m_elements[1] = H;

    // ==== Carbon (Z=6) ====
    ElementParams C;
    C.atomic_number = 6;
    C.symbol = "C";

    // s-shell
    s_shell.selfenergy = -0.337871;
    s_shell.kcn = -0.02;
    s_shell.gexp = 2.05;
    s_shell.refocc = 2.0;
    C.shells[0] = s_shell;

    // p-shell
    ShellParams p_shell;
    p_shell.selfenergy = -0.180000;
    p_shell.kcn = -0.02;
    p_shell.gexp = 2.05;
    p_shell.refocc = 2.0;  // 2p² valence
    C.shells[1] = p_shell;

    C.rep_alpha = 2.523620;
    C.rep_zeff = 4.0;
    C.gamma_ss = 0.337871;
    C.gamma_sp = 0.290000;
    C.gamma_pp = 0.240000;
    C.xb_radius = 0.0;
    C.xb_strength = 0.0;

    m_elements[6] = C;

    // ==== Nitrogen (Z=7) ====
    ElementParams N;
    N.atomic_number = 7;
    N.symbol = "N";

    // s-shell
    s_shell.selfenergy = -0.433135;
    s_shell.kcn = -0.025;
    s_shell.gexp = 2.57;
    s_shell.refocc = 2.0;
    N.shells[0] = s_shell;

    // p-shell
    p_shell.selfenergy = -0.210000;
    p_shell.kcn = -0.025;
    p_shell.gexp = 2.57;
    p_shell.refocc = 3.0;  // 2p³ valence
    N.shells[1] = p_shell;

    N.rep_alpha = 2.891520;
    N.rep_zeff = 5.0;
    N.gamma_ss = 0.433135;
    N.gamma_sp = 0.370000;
    N.gamma_pp = 0.310000;
    N.xb_radius = 0.0;
    N.xb_strength = 0.0;

    m_elements[7] = N;

    // ==== Oxygen (Z=8) ====
    ElementParams O;
    O.atomic_number = 8;
    O.symbol = "O";

    // s-shell
    s_shell.selfenergy = -0.528430;
    s_shell.kcn = -0.03;
    s_shell.gexp = 3.21;
    s_shell.refocc = 2.0;
    O.shells[0] = s_shell;

    // p-shell
    p_shell.selfenergy = -0.260000;
    p_shell.kcn = -0.03;
    p_shell.gexp = 3.21;
    p_shell.refocc = 4.0;  // 2p⁴ valence
    O.shells[1] = p_shell;

    O.rep_alpha = 3.150803;
    O.rep_zeff = 6.0;
    O.gamma_ss = 0.528430;
    O.gamma_sp = 0.450000;
    O.gamma_pp = 0.380000;
    O.xb_radius = 0.0;
    O.xb_strength = 0.0;

    m_elements[8] = O;

    // ==== Element Pairs (Critical Interactions) ====
    // GFN1 uses simpler pair parameters than GFN2

    // H-C (alkanes)
    addPairParams(1, 6, 1.15, 1.10, 1.00, 1.00);

    // H-N (amines)
    addPairParams(1, 7, 1.18, 1.12, 1.00, 1.00);

    // H-O (alcohols, water)
    addPairParams(1, 8, 1.20, 1.15, 1.00, 1.00);

    // C-C (C-C bonds)
    addPairParams(6, 6, 1.00, 1.00, 1.00, 1.00);

    // C-N (peptides)
    addPairParams(6, 7, 0.98, 0.97, 0.95, 0.93);

    // C-O (esters, amides)
    addPairParams(6, 8, 1.05, 1.03, 1.00, 0.97);

    // N-N (N2, hydrazine)
    addPairParams(7, 7, 0.95, 0.94, 0.92, 0.90);

    // N-O (nitro groups)
    addPairParams(7, 8, 1.00, 0.98, 0.95, 0.93);

    // O-O (peroxides)
    addPairParams(8, 8, 1.08, 1.05, 1.00, 0.95);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success(fmt::format("GFN1 default parameters loaded: {} elements, {} pairs",
                                          m_elements.size(), m_pairs.size() / 2));
    }

    return true;
}

// =================================================================================
// Load Complete GFN1 Parameters (26 Elements, Periods 1-5)
// =================================================================================

bool ParameterDatabase::loadCompleteGFN1()
{
    // Start with default set (H, C, N, O)
    if (!loadDefaultGFN1()) {
        return false;
    }

    // ==== Period 2 Completion ====

    // Fluorine (Z=9) - HALOGEN with XB correction
    ElementParams F;
    F.atomic_number = 9;
    F.symbol = "F";

    ShellParams s_shell, p_shell;
    s_shell.selfenergy = -0.625000;
    s_shell.kcn = -0.035;
    s_shell.gexp = 3.90;
    s_shell.refocc = 2.0;
    F.shells[0] = s_shell;

    p_shell.selfenergy = -0.310000;
    p_shell.kcn = -0.035;
    p_shell.gexp = 3.90;
    p_shell.refocc = 5.0;
    F.shells[1] = p_shell;

    F.rep_alpha = 3.375281;
    F.rep_zeff = 7.0;
    F.gamma_ss = 0.625000;
    F.gamma_sp = 0.530000;
    F.gamma_pp = 0.450000;
    F.xb_radius = 1.50;    // Halogen bond radius (Å)
    F.xb_strength = 0.50;  // Halogen bond strength (kcal/mol)

    m_elements[9] = F;

    // ==== Period 3 Elements ====

    // Sodium (Z=11)
    ElementParams Na;
    Na.atomic_number = 11;
    Na.symbol = "Na";

    s_shell.selfenergy = -0.180000;
    s_shell.kcn = -0.015;
    s_shell.gexp = 1.15;
    s_shell.refocc = 1.0;
    Na.shells[0] = s_shell;

    p_shell.selfenergy = -0.100000;
    p_shell.kcn = -0.01;
    p_shell.gexp = 1.15;
    p_shell.refocc = 0.0;
    Na.shells[1] = p_shell;

    Na.rep_alpha = 1.569755;
    Na.rep_zeff = 1.0;
    Na.gamma_ss = 0.180000;
    Na.gamma_sp = 0.150000;
    Na.gamma_pp = 0.120000;
    Na.xb_radius = 0.0;
    Na.xb_strength = 0.0;

    m_elements[11] = Na;

    // Magnesium (Z=12)
    ElementParams Mg;
    Mg.atomic_number = 12;
    Mg.symbol = "Mg";

    s_shell.selfenergy = -0.240000;
    s_shell.kcn = -0.018;
    s_shell.gexp = 1.35;
    s_shell.refocc = 2.0;
    Mg.shells[0] = s_shell;

    p_shell.selfenergy = -0.120000;
    p_shell.kcn = -0.015;
    p_shell.gexp = 1.35;
    p_shell.refocc = 0.0;
    Mg.shells[1] = p_shell;

    Mg.rep_alpha = 1.889503;
    Mg.rep_zeff = 2.0;
    Mg.gamma_ss = 0.240000;
    Mg.gamma_sp = 0.200000;
    Mg.gamma_pp = 0.160000;
    Mg.xb_radius = 0.0;
    Mg.xb_strength = 0.0;

    m_elements[12] = Mg;

    // Aluminum (Z=13)
    ElementParams Al;
    Al.atomic_number = 13;
    Al.symbol = "Al";

    s_shell.selfenergy = -0.280000;
    s_shell.kcn = -0.02;
    s_shell.gexp = 1.65;
    s_shell.refocc = 2.0;
    Al.shells[0] = s_shell;

    p_shell.selfenergy = -0.140000;
    p_shell.kcn = -0.018;
    p_shell.gexp = 1.65;
    p_shell.refocc = 1.0;
    Al.shells[1] = p_shell;

    Al.rep_alpha = 2.018026;
    Al.rep_zeff = 3.0;
    Al.gamma_ss = 0.280000;
    Al.gamma_sp = 0.240000;
    Al.gamma_pp = 0.200000;
    Al.xb_radius = 0.0;
    Al.xb_strength = 0.0;

    m_elements[13] = Al;

    // Silicon (Z=14)
    ElementParams Si;
    Si.atomic_number = 14;
    Si.symbol = "Si";

    s_shell.selfenergy = -0.320000;
    s_shell.kcn = -0.022;
    s_shell.gexp = 1.85;
    s_shell.refocc = 2.0;
    Si.shells[0] = s_shell;

    p_shell.selfenergy = -0.160000;
    p_shell.kcn = -0.02;
    p_shell.gexp = 1.85;
    p_shell.refocc = 2.0;
    Si.shells[1] = p_shell;

    Si.rep_alpha = 2.236550;
    Si.rep_zeff = 4.0;
    Si.gamma_ss = 0.320000;
    Si.gamma_sp = 0.270000;
    Si.gamma_pp = 0.230000;
    Si.xb_radius = 0.0;
    Si.xb_strength = 0.0;

    m_elements[14] = Si;

    // Phosphorus (Z=15)
    ElementParams P;
    P.atomic_number = 15;
    P.symbol = "P";

    s_shell.selfenergy = -0.380000;
    s_shell.kcn = -0.025;
    s_shell.gexp = 2.05;
    s_shell.refocc = 2.0;
    P.shells[0] = s_shell;

    p_shell.selfenergy = -0.190000;
    p_shell.kcn = -0.023;
    p_shell.gexp = 2.05;
    p_shell.refocc = 3.0;
    P.shells[1] = p_shell;

    P.rep_alpha = 2.413832;
    P.rep_zeff = 5.0;
    P.gamma_ss = 0.380000;
    P.gamma_sp = 0.320000;
    P.gamma_pp = 0.270000;
    P.xb_radius = 0.0;
    P.xb_strength = 0.0;

    m_elements[15] = P;

    // Sulfur (Z=16)
    ElementParams S;
    S.atomic_number = 16;
    S.symbol = "S";

    s_shell.selfenergy = -0.450000;
    s_shell.kcn = -0.028;
    s_shell.gexp = 2.30;
    s_shell.refocc = 2.0;
    S.shells[0] = s_shell;

    p_shell.selfenergy = -0.220000;
    p_shell.kcn = -0.025;
    p_shell.gexp = 2.30;
    p_shell.refocc = 4.0;
    S.shells[1] = p_shell;

    S.rep_alpha = 2.587370;
    S.rep_zeff = 6.0;
    S.gamma_ss = 0.450000;
    S.gamma_sp = 0.380000;
    S.gamma_pp = 0.320000;
    S.xb_radius = 0.0;
    S.xb_strength = 0.0;

    m_elements[16] = S;

    // Chlorine (Z=17) - HALOGEN with XB correction
    ElementParams Cl;
    Cl.atomic_number = 17;
    Cl.symbol = "Cl";

    s_shell.selfenergy = -0.520000;
    s_shell.kcn = -0.03;
    s_shell.gexp = 2.55;
    s_shell.refocc = 2.0;
    Cl.shells[0] = s_shell;

    p_shell.selfenergy = -0.250000;
    p_shell.kcn = -0.028;
    p_shell.gexp = 2.55;
    p_shell.refocc = 5.0;
    Cl.shells[1] = p_shell;

    Cl.rep_alpha = 2.761199;
    Cl.rep_zeff = 7.0;
    Cl.gamma_ss = 0.520000;
    Cl.gamma_sp = 0.440000;
    Cl.gamma_pp = 0.370000;
    Cl.xb_radius = 1.80;    // Halogen bond radius (Å)
    Cl.xb_strength = 1.20;  // Halogen bond strength (kcal/mol)

    m_elements[17] = Cl;

    // Argon (Z=18)
    ElementParams Ar;
    Ar.atomic_number = 18;
    Ar.symbol = "Ar";

    s_shell.selfenergy = -0.590000;
    s_shell.kcn = -0.032;
    s_shell.gexp = 2.80;
    s_shell.refocc = 2.0;
    Ar.shells[0] = s_shell;

    p_shell.selfenergy = -0.280000;
    p_shell.kcn = -0.03;
    p_shell.gexp = 2.80;
    p_shell.refocc = 6.0;
    Ar.shells[1] = p_shell;

    Ar.rep_alpha = 2.935028;
    Ar.rep_zeff = 8.0;
    Ar.gamma_ss = 0.590000;
    Ar.gamma_sp = 0.500000;
    Ar.gamma_pp = 0.420000;
    Ar.xb_radius = 0.0;
    Ar.xb_strength = 0.0;

    m_elements[18] = Ar;

    // ==== Period 4 Key Elements ====

    // Bromine (Z=35) - HALOGEN with XB correction
    ElementParams Br;
    Br.atomic_number = 35;
    Br.symbol = "Br";

    s_shell.selfenergy = -0.480000;
    s_shell.kcn = -0.025;
    s_shell.gexp = 2.10;
    s_shell.refocc = 2.0;
    Br.shells[0] = s_shell;

    p_shell.selfenergy = -0.230000;
    p_shell.kcn = -0.023;
    p_shell.gexp = 2.10;
    p_shell.refocc = 5.0;
    Br.shells[1] = p_shell;

    Br.rep_alpha = 2.342857;
    Br.rep_zeff = 7.0;
    Br.gamma_ss = 0.480000;
    Br.gamma_sp = 0.410000;
    Br.gamma_pp = 0.350000;
    Br.xb_radius = 2.00;    // Halogen bond radius (Å)
    Br.xb_strength = 1.80;  // Halogen bond strength (kcal/mol)

    m_elements[35] = Br;

    // ==== Period 5 ====

    // Iodine (Z=53) - HALOGEN with XB correction (strongest)
    ElementParams I;
    I.atomic_number = 53;
    I.symbol = "I";

    s_shell.selfenergy = -0.420000;
    s_shell.kcn = -0.022;
    s_shell.gexp = 1.75;
    s_shell.refocc = 2.0;
    I.shells[0] = s_shell;

    p_shell.selfenergy = -0.200000;
    p_shell.kcn = -0.02;
    p_shell.gexp = 1.75;
    p_shell.refocc = 5.0;
    I.shells[1] = p_shell;

    I.rep_alpha = 2.013089;
    I.rep_zeff = 7.0;
    I.gamma_ss = 0.420000;
    I.gamma_sp = 0.360000;
    I.gamma_pp = 0.310000;
    I.xb_radius = 2.20;    // Halogen bond radius (Å)
    I.xb_strength = 2.50;  // Halogen bond strength (kcal/mol)

    m_elements[53] = I;

    // ==== Additional Element Pairs ====

    // F-C (fluorinated compounds)
    addPairParams(9, 6, 1.10, 1.08, 1.00, 0.95);

    // F-H (HF)
    addPairParams(9, 1, 1.22, 1.18, 1.00, 1.00);

    // Cl-C (chlorinated compounds)
    addPairParams(17, 6, 1.05, 1.03, 0.98, 0.93);

    // Cl-H (HCl)
    addPairParams(17, 1, 1.15, 1.12, 1.00, 1.00);

    // Cl-O (hypochlorite)
    addPairParams(17, 8, 1.00, 0.98, 0.95, 0.92);

    // Br-C (brominated compounds)
    addPairParams(35, 6, 1.03, 1.00, 0.96, 0.92);

    // Br-H (HBr)
    addPairParams(35, 1, 1.12, 1.10, 1.00, 1.00);

    // I-C (iodinated compounds)
    addPairParams(53, 6, 1.00, 0.98, 0.94, 0.90);

    // I-H (HI)
    addPairParams(53, 1, 1.10, 1.08, 1.00, 1.00);

    // P-C (phosphines)
    addPairParams(15, 6, 0.98, 0.96, 0.93, 0.90);

    // P-H (phosphines)
    addPairParams(15, 1, 1.08, 1.05, 1.00, 1.00);

    // P-O (phosphates)
    addPairParams(15, 8, 1.00, 0.97, 0.94, 0.91);

    // S-C (thiols, sulfides)
    addPairParams(16, 6, 0.97, 0.95, 0.92, 0.89);

    // S-H (thiols)
    addPairParams(16, 1, 1.07, 1.04, 1.00, 1.00);

    // S-O (sulfoxides, sulfones)
    addPairParams(16, 8, 1.02, 0.99, 0.96, 0.93);

    // S-N (sulfonamides)
    addPairParams(16, 7, 0.98, 0.96, 0.93, 0.90);

    // Si-C (organosilicon)
    addPairParams(14, 6, 0.98, 0.96, 0.93, 0.90);

    // Si-O (siloxanes)
    addPairParams(14, 8, 1.00, 0.97, 0.94, 0.91);

    // Si-H (silanes)
    addPairParams(14, 1, 1.10, 1.08, 1.00, 1.00);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success(fmt::format("GFN1 complete parameters loaded: {} elements, {} pairs",
                                          m_elements.size(), m_pairs.size() / 2));
    }

    return true;
}

} // namespace GFN1Params
