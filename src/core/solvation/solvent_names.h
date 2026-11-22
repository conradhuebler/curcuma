/*
 * <Solvent names and solution states for the native ALPB/GBSA solvation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Claude Generated (Nov 2026). AI-generated, machine-tested; human production testing pending.
 */

#pragma once

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>

#include "src/core/units.h"

namespace curcuma {
namespace solvation {

/**
 * @brief Canonical solvent key from any accepted spelling.
 *
 * The input is lower-cased (xtb does the same: `Water`, `CHCl3` work there) and mapped with the alias list of the xtb source tree
 * (`normalizeSolventName`, xtb/src/solv/model.f90; the 6.7.1 binary knows a shorter list) plus `dichlormethane` (accepted by the xtb
 * parameter loader). curcuma
 * additionally accepts the spellings that docs/SOLVATION.md has always listed: `diethyl ether`, `carbon disulfide`,
 * `ethyl acetate`, `octanol wet`, `wet octanol`, `wetoctanol`. Spaces and underscores are removed before the lookup.
 * The result is the file-name key of the parameter tables (`ch2cl2`, `chcl3`, `cs2`, `ether`, `hexane`, `woctanol`, ...);
 * an unknown name is returned unchanged (lower-cased), so the caller's table lookup fails.
 */
inline std::string normalizeSolventName(const std::string& input)
{
    std::string s;
    s.reserve(input.size());
    for (unsigned char c : input) {
        if (c == ' ' || c == '_' || c == '\t')
            continue;
        s += static_cast<char>(std::tolower(c));
    }
    struct Alias { const char* canon; const char* names; };  // names separated by '|'
    static const Alias kAliases[] = {
        { "acetonitrile", "methylcyanide" },
        { "ch2cl2", "methylenechloride|dichloromethane|dcm|dichlormethane" },
        { "chcl3", "chloroform|trichloromethane|tcm" },
        { "cs2", "carbondisulfide" },
        { "ether", "diethylether" },
        { "dmf", "dimethylformamide" },
        { "dmso", "dimethylsulfoxide" },
        { "furane", "furan" },
        { "hexane", "nhexan|n-hexan|nhexane|n-hexane" },
        { "thf", "tetrahydrofuran" },
        { "water", "h2o" },
        { "ethylacetate", "ethyl-acetate" },
        { "woctanol", "octanolwet|wetoctanol|octanol(wet)" },
    };
    for (const Alias& a : kAliases) {
        const std::string names = std::string("|") + a.names + "|";
        if (names.find("|" + s + "|") != std::string::npos)
            return a.canon;
    }
    return s;
}

/// True for "no solvent": none, vacuum, gas, gasphase or an empty name (case-insensitive). docs/SOLVATION.md lists `vacuum`
/// as a synonym of `none`.
inline bool isGasPhase(const std::string& name)
{
    std::string s;
    for (unsigned char c : name)
        if (c != ' ' && c != '_')
            s += static_cast<char>(std::tolower(c));
    return s.empty() || s == "none" || s == "vacuum" || s == "gas" || s == "gasphase";
}

/// Reference state of the solvation free energy (xtb `--alpb water [reference|bar1M|gsolv]`, tblite `solution_state`).
enum class SolutionState {
    GSolv,      ///< gas phase and solution at the same concentration (xtb and tblite default)
    Bar1M,      ///< 1 bar ideal gas to 1 mol/L solution
    Reference   ///< 1 bar ideal gas to the solvent at its own concentration (xtb "reference")
};

/// Parse a state name (case-insensitive): gsolv, bar1M, reference. Returns false for anything else.
inline bool parseSolutionState(const std::string& name, SolutionState& out)
{
    std::string s;
    for (unsigned char c : name)
        s += static_cast<char>(std::tolower(c));
    if (s.empty() || s == "gsolv") { out = SolutionState::GSolv; return true; }
    if (s == "bar1m" || s == "bar1mol" || s == "mol1bar") { out = SolutionState::Bar1M; return true; }
    if (s == "reference") { out = SolutionState::Reference; return true; }
    return false;
}

inline const char* solutionStateName(SolutionState s)
{
    return s == SolutionState::Bar1M ? "bar1M" : (s == SolutionState::Reference ? "reference" : "gsolv");
}

/// Ideal-gas molar volume in L/mol at 298.15 K and 1 bar, the rounded value of xtb 6.7.1 and pprcht/gfnff.
constexpr double kIdealGasMolVolume = 24.79;

/**
 * @brief Shift of the solvation free energy for the reference state, in Hartree.
 *
 * xtb/state.f90 and tblite/solvation/shift.f90 (get_stateshift):
 *   gsolv:     0
 *   bar1M:     T kB ln(Vig T / T0)
 *   reference: T kB [ ln(Vig T / T0) + ln( (density/refDensity) * (refMolMass/molMass) ) ]
 * with T0 = 298.15 K, refDensity = 1e-3 and refMolMass = 1 as written in the reference sources (the "au" comments there
 * are not conversions), and `density` the solvent number density in the reference's atomic units (rhos/smass * 8.92388e-2).
 * One implementation for every method (GFN1, GFN2, GFN-FF): the ideal-gas molar volume is xtb's 24.79 L/mol
 * (kIdealGasMolVolume; tblite derives 24.78960 from R, a 1.6e-5 relative difference, 1.4e-8 Eh in the bar1M shift, not
 * ported separately). With CODATA-2018 kB the bar1M shift is 3.0312492e-3 Eh, xtb 6.7.1 gives 3.0312463e-3 Eh.
 *
 * @param density     solvent density as stored by the loaders (rhos/smass * molcm3toau)
 * @param molar_mass  solvent molar mass in g/mol
 */
inline double solutionStateShift(SolutionState state, double temperature, double density, double molar_mass)
{
    constexpr double vig = kIdealGasMolVolume;
    constexpr double kB_au = CurcumaUnit::Constants::BOLTZMANN_HARTREE_K;   // CODATA 2018 kB / Hartree, Eh/K
    constexpr double T0 = 298.15;
    constexpr double refDensity = 1.0e-3, refMolMass = 1.0;
    const double bar = std::log(vig * temperature / T0);
    switch (state) {
    case SolutionState::Bar1M:
        return temperature * kB_au * bar;
    case SolutionState::Reference:
        return temperature * kB_au * (bar + std::log(density / refDensity * refMolMass / molar_mass));
    default:
        return 0.0;
    }
}


} // namespace solvation
} // namespace curcuma
