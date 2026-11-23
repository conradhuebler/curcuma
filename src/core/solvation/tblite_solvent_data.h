/*
 * <Dielectric constants of tblite's solvent table, used by the native CPCM>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Claude Generated (Nov 2026). AI-generated, machine-tested; human production testing pending.
 *
 * Source: tblite 0.4.0 (commit 9ca469a), src/tblite/solvation/data.f90, get_solvent_data: the first column
 * (dielectric constant, Minnesota solvent descriptor database) of every case label, aliases included; the labels
 * "n-hexan"/"n-hexane" are left out because tblite's compact() removes the '-' before the lookup, so they never match.
 * tblite's `--cpcm <name>` takes epsilon from this table, not from the ALPB parameter sets (water: 78.36 here,
 * 80.2 in the ALPB set that ALPB and xtb --cosmo use).
 */

#pragma once

#include <cctype>
#include <limits>
#include <string>

#include "src/core/solvation/solvent_names.h"

namespace curcuma {
namespace solvation {

namespace detail {
struct TbliteSolventEps {
    const char* name;
    double eps;
};
inline constexpr TbliteSolventEps kTbliteSolventEps[] = {
    { "inf", std::numeric_limits<double>::max() },   // tblite: huge(1.0_wp), conductor limit
    { "2methylpyridine", 9.9533 },
    { "4methyl2pentanone", 12.8871 },
    { "aceticacid", 6.2528 },
    { "acetone", 20.493 },
    { "acetonitrile", 35.6881 },
    { "mecn", 35.6881 },
    { "acetophenone", 17.44 },
    { "aniline", 6.8882 },
    { "anisole", 4.2247 },
    { "benzaldehyde", 17.85 },
    { "benzene", 2.2706 },
    { "benzonitrile", 25.592 },
    { "benzylalcohol", 12.4569 },
    { "bromobenzene", 5.3954 },
    { "bromoethane", 9.01 },
    { "bromoform", 4.2488 },
    { "bromooctane", 5.0244 },
    { "butanol", 17.3323 },
    { "butanone", 18.2457 },
    { "butylacetate", 4.9941 },
    { "butylbenzene", 2.36 },
    { "carbondisulfide", 2.6105 },
    { "cs2", 2.6105 },
    { "carbontet", 2.228 },
    { "chlorobenzene", 5.6968 },
    { "chloroform", 4.7113 },
    { "chcl3", 4.7113 },
    { "chlorohexane", 5.9491 },
    { "cyclohexane", 2.0165 },
    { "cyclohexanone", 15.6186 },
    { "decalin", 2.196 },
    { "decane", 1.9846 },
    { "decanol", 7.5305 },
    { "dibromoethane", 4.9313 },
    { "dibutylether", 3.0473 },
    { "dichloroethane", 10.125 },
    { "diethylether", 4.24 },
    { "ether", 4.24 },
    { "diisopropylether", 3.38 },
    { "dimethylacetamide", 37.7807 },
    { "dimethylformamide", 37.219 },
    { "dmf", 37.219 },
    { "dimethylpyridine", 7.1735 },
    { "dimethylsulfoxide", 46.826 },
    { "dmso", 46.826 },
    { "dioxane", 2.2099 },
    { "dodecane", 2.006 },
    { "ethanol", 24.852 },
    { "ethoxybenzene", 4.1797 },
    { "ethylacetate", 5.9867 },
    { "ethylbenzene", 2.4339 },
    { "fluorobenzene", 5.42 },
    { "fluoroctane", 3.89 },
    { "furan", 2.953 },
    { "furane", 2.953 },
    { "heptane", 1.9113 },
    { "heptanol", 11.321 },
    { "hexadecane", 2.0402 },
    { "hexadecyliodide", 3.5338 },
    { "hexane", 1.8819 },
    { "nhexan", 1.8819 },
    { "nhexane", 1.8819 },
    { "hexanol", 12.5102 },
    { "iodobenzene", 4.547 },
    { "isobutanol", 16.7766 },
    { "isooctane", 1.9358 },
    { "isopropanol", 19.2645 },
    { "isopropylbenzene", 2.3712 },
    { "mcresol", 12.44 },
    { "mesitylene", 2.265 },
    { "methanol", 32.613 },
    { "methoxyethanol", 17.2 },
    { "methylenechloride", 8.93 },
    { "dichloromethane", 8.93 },
    { "dichlormethane", 8.93 },
    { "ch2cl2", 8.93 },
    { "methylformamide", 181.5619 },
    { "nitrobenzene", 34.8091 },
    { "nitroethane", 28.2896 },
    { "nitromethane", 36.5623 },
    { "nonane", 1.9605 },
    { "nonanol", 8.5991 },
    { "octane", 1.9406 },
    { "octanol", 9.8629 },
    { "odichlorobenzene", 9.9949 },
    { "onitrotoluene", 25.6692 },
    { "pentadecane", 2.0333 },
    { "pentane", 1.8371 },
    { "pentanol", 15.13 },
    { "perfluorobenzene", 2.029 },
    { "phenol", 12.4 },
    { "phenylether", 3.73 },
    { "propanol", 20.5237 },
    { "pyridine", 12.9776 },
    { "secbutanol", 15.9436 },
    { "secbutylbenzene", 2.3446 },
    { "tbutylbenzene", 2.3447 },
    { "tetrachloroethene", 2.268 },
    { "tetrahydrofuran", 7.4257 },
    { "thf", 7.4257 },
    { "tetrahydrothiophenedioxide", 43.9622 },
    { "tetralin", 2.771 },
    { "toluene", 2.3741 },
    { "tributylphosphate", 8.1781 },
    { "triethylamine", 2.3832 },
    { "trimethylbenzene", 2.3653 },
    { "undecane", 1.991 },
    { "water", 78.36 },
    { "h2o", 78.36 },
    { "woctanol", 9.8629 },
    { "xylene", 2.3879 },
};
} // namespace detail

/**
 * @brief Dielectric constant of a named solvent as tblite's get_solvent_data returns it; -1 if unknown.
 *
 * The name is compacted as tblite does (lower case, '-', ',' and blanks removed; curcuma also drops '_'). A name that
 * is not in the table is retried with curcuma's canonical key (normalizeSolventName), so the xtb aliases accepted by
 * ALPB and COSMO (`dcm`, `tcm`, `wet octanol`, ...) work here as well. Claude Generated (Nov 2026).
 */
inline double tbliteSolventEpsilon(const std::string& input)
{
    std::string s;
    for (unsigned char c : input) {
        if (c == '-' || c == ',' || c == ' ' || c == '_' || c == '\t')
            continue;
        s += static_cast<char>(std::tolower(c));
    }
    for (const auto& e : detail::kTbliteSolventEps)
        if (s == e.name)
            return e.eps;
    const std::string key = normalizeSolventName(input);
    for (const auto& e : detail::kTbliteSolventEps)
        if (key == e.name)
            return e.eps;
    return -1.0;
}

} // namespace solvation
} // namespace curcuma
