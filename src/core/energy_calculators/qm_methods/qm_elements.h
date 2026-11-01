/*
 * <Native QM: element symbols and basis-file element names>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Nov 2026). AI-generated, machine-tested.
 *
 * One table for the native QM engine (QMEngine, BasisSetParser, qmri aux parser, ECP parser): standard symbol and
 * the upper-case element name the Basis Set Exchange writes as block header in the '$DATA' (GAMESS-US) export
 * (ALUMINIUM, CAESIUM, SULFUR spelled as BSE spells them). Z = 1..103.
 */

#pragma once

#include <string>

namespace qmelem {

inline constexpr int kMaxZ = 103;

/// Standard symbol of element Z (1..103), nullptr outside.
inline const char* symbol(int Z)
{
    static const char* const s[kMaxZ + 1] = { nullptr,
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
        "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce",
        "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
        "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
        "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr" };
    return (Z >= 1 && Z <= kMaxZ) ? s[Z] : nullptr;
}

/// Upper-case element name as in the BSE '$DATA' export (HYDROGEN ... LAWRENCIUM), nullptr outside 1..103.
inline const char* bseName(int Z)
{
    static const char* const n[kMaxZ + 1] = { nullptr,
        "HYDROGEN", "HELIUM", "LITHIUM", "BERYLLIUM", "BORON", "CARBON", "NITROGEN", "OXYGEN", "FLUORINE", "NEON",
        "SODIUM", "MAGNESIUM", "ALUMINIUM", "SILICON", "PHOSPHORUS", "SULFUR", "CHLORINE", "ARGON", "POTASSIUM",
        "CALCIUM", "SCANDIUM", "TITANIUM", "VANADIUM", "CHROMIUM", "MANGANESE", "IRON", "COBALT", "NICKEL", "COPPER",
        "ZINC", "GALLIUM", "GERMANIUM", "ARSENIC", "SELENIUM", "BROMINE", "KRYPTON", "RUBIDIUM", "STRONTIUM",
        "YTTRIUM", "ZIRCONIUM", "NIOBIUM", "MOLYBDENUM", "TECHNETIUM", "RUTHENIUM", "RHODIUM", "PALLADIUM", "SILVER",
        "CADMIUM", "INDIUM", "TIN", "ANTIMONY", "TELLURIUM", "IODINE", "XENON", "CAESIUM", "BARIUM", "LANTHANUM",
        "CERIUM", "PRASEODYMIUM", "NEODYMIUM", "PROMETHIUM", "SAMARIUM", "EUROPIUM", "GADOLINIUM", "TERBIUM",
        "DYSPROSIUM", "HOLMIUM", "ERBIUM", "THULIUM", "YTTERBIUM", "LUTETIUM", "HAFNIUM", "TANTALUM", "TUNGSTEN",
        "RHENIUM", "OSMIUM", "IRIDIUM", "PLATINUM", "GOLD", "MERCURY", "THALLIUM", "LEAD", "BISMUTH", "POLONIUM",
        "ASTATINE", "RADON", "FRANCIUM", "RADIUM", "ACTINIUM", "THORIUM", "PROTACTINIUM", "URANIUM", "NEPTUNIUM",
        "PLUTONIUM", "AMERICIUM", "CURIUM", "BERKELIUM", "CALIFORNIUM", "EINSTEINIUM", "FERMIUM", "MENDELEVIUM",
        "NOBELIUM", "LAWRENCIUM" };
    return (Z >= 1 && Z <= kMaxZ) ? n[Z] : nullptr;
}

/// Atomic number of a symbol or BSE name (any letter case; ALUMINUM/CESIUM/SULPHUR accepted too), 0 if unknown.
inline int atomicNumber(const std::string& token)
{
    std::string up;
    for (char c : token) up += static_cast<char>((c >= 'a' && c <= 'z') ? c - 32 : c);
    if (up == "ALUMINUM") return 13;
    if (up == "SULPHUR") return 16;
    if (up == "CESIUM") return 55;
    for (int Z = 1; Z <= kMaxZ; ++Z) {
        if (up == bseName(Z)) return Z;
        std::string s = symbol(Z);
        for (char& c : s) c = static_cast<char>((c >= 'a' && c <= 'z') ? c - 32 : c);
        if (up == s) return Z;
    }
    return 0;
}

}  // namespace qmelem
