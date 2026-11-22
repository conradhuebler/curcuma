/*
 * <CPCM-X: parameters, solvent tables, solvent database and solute input>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Claude Generated (Nov 2026): data of CPCM-X 1.1.0 (LGPL-3.0-or-later):
 *   src/cpcmx/internaldb.f90  parameter sets xtb_water / xtb_other, internal SMD parameter text, internalcosmo
 *   src/cpcmx/data.f90        solvent_name, density, atomicMassNist
 *   src/cpcmx/initialize.f90  read_cosmo (solute file), load_solvent, setup_cov (covalent radii)
 * and of numsa 0.2.0 src/smd/init.f90 (SMD solvent descriptors). The tables were extracted from the Fortran sources by
 * script; the solvent database files themselves are read at run time (loadSolvent).
 * AI-generated, machine-tested only; human production testing pending.
 */

#include "cpcmx.h"
#include "cpcmx_internal.h"

#include "src/core/elements.h"
#include "src/core/energy_calculators/ff_methods/alpb_parameters.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

#if defined(__linux__)
#include <climits>
#include <unistd.h>
#endif

namespace curcuma::solvation::cpcmx {

namespace {

// data.f90 solvent_name: alias -> Minnesota name (generated from CPCM-X 1.1.0)
const std::pair<const char*, const char*> kSolventNames[] = {
    {"2methylpyridine", "2methylpyridine"}, {"4methyl2pentanone", "4methyl2pentanone"}, {"aceticacid", "aceticacid"},
    {"ethanoicacid", "aceticacid"}, {"acetonitrile", "acetonitrile"}, {"methylcyanide", "acetonitrile"},
    {"acetophenone", "acetophenone"}, {"aniline", "aniline"}, {"anisole", "anisole"}, {"methoxybenzene", "anisole"},
    {"benzene", "benzene"}, {"benzonitrile", "benzonitrile"}, {"cyanobenzene", "benzonitrile"},
    {"benzylalcohol", "benzylalcohol"}, {"bromobenzene", "bromobenzene"}, {"bromoethane", "bromoethane"},
    {"bromoform", "bromoform"}, {"tribromomethane", "bromoform"}, {"bromooctane", "bromooctane"},
    {"butanol", "butanol"}, {"butanone", "butanone"}, {"butylacetate", "butylacetate"},
    {"butylbenzene", "butylbenzene"}, {"carbondisulfide", "carbondisulfide"}, {"cs2", "carbondisulfide"},
    {"carbontet", "carbontet"}, {"carbontetrachloride", "carbontet"}, {"ccl4", "carbontet"},
    {"tetrachloromethane", "carbontet"}, {"chlorobenzene", "chlorobenzene"}, {"chloroform", "chloroform"},
    {"chcl3", "chloroform"}, {"trichloromethane", "chloroform"}, {"tcm", "chloroform"},
    {"chlorohexane", "chlorohexane"}, {"cyclohexane", "cyclohexane"}, {"cyclohexanone", "cyclohexanone"},
    {"decalin", "decalin"}, {"decane", "decane"}, {"decanol", "decanol"}, {"dibromoethane", "dibromoethane"},
    {"dibutylether", "dibutylether"}, {"dichloroethane", "dichloroethane"}, {"diethylether", "diethylether"},
    {"ether", "diethylether"}, {"diisopropylether", "diisopropylether"}, {"dimethylacetamide", "dimethylacetamide"},
    {"dimethylformamide", "dimethylformamide"}, {"dmf", "dimethylformamide"},
    {"dimethylpyridine", "dimethylpyridine"}, {"dimethylsulfoxide", "dimethylsulfoxide"},
    {"dmso", "dimethylsulfoxide"}, {"dodecane", "dodecane"}, {"ethanol", "ethanol"},
    {"ethoxybenzene", "ethoxybenzene"}, {"ethylacetate", "ethylacetate"}, {"ethylbenzene", "ethylbenzene"},
    {"fluorobenzene", "fluorobenzene"}, {"fluoroctane", "fluoroctane"}, {"heptane", "heptane"},
    {"heptanol", "heptanol"}, {"hexadecane", "hexadecane"}, {"hexadecyliodide", "hexadecyliodide"},
    {"iodohexadecane", "hexadecyliodide"}, {"1-iodohexadecane", "hexadecyliodide"}, {"hexane", "hexane"},
    {"hexanol", "hexanol"}, {"iodobenzene", "iodobenzene"}, {"isobutanol", "isobutanol"}, {"isooctane", "isooctane"},
    {"2,2,4-trimethylpentane", "isooctane"}, {"isopropanol", "isopropanol"}, {"isopropylbenzene", "isopropylbenzene"},
    {"isopropyltoluene", "isopropyltoluene"}, {"mcresol", "mcresol"}, {"3-methylphenol", "mcresol"},
    {"mesitylene", "mesitylene"}, {"1,3,5-trimethylbenzene", "mesitylene"}, {"methoxyethanol", "methoxyethanol"},
    {"methylenechloride", "methylenechloride"}, {"dichloromethane", "methylenechloride"},
    {"dcm", "methylenechloride"}, {"ch2cl2", "methylenechloride"}, {"methylformamide", "methylformamide"},
    {"nitrobenzene", "nitrobenzene"}, {"nitroethane", "nitroethane"}, {"nitromethane", "nitromethane"},
    {"nonane", "nonane"}, {"nonanol", "nonanol"}, {"octane", "octane"}, {"octanol", "octanol"},
    {"odichlorobenzene", "odichlorobenzene"}, {"odcb", "odichlorobenzene"}, {"onitrotoluene", "onitrotoluene"},
    {"pentadecane", "pentadecane"}, {"pentane", "pentane"}, {"pentanol", "pentanol"},
    {"perfluorobenzene", "perfluorobenzene"}, {"hexaflurobenzene", "perfluorobenzene"},
    {"phenylether", "phenylether"}, {"propanol", "propanol"}, {"pyridine", "pyridine"}, {"secbutanol", "secbutanol"},
    {"secbutylbenzene", "secbutylbenzene"}, {"tbutylbenzene", "tbutylbenzene"},
    {"tetrachloroethene", "tetrachloroethene"}, {"c2cl4", "tetrachloroethene"},
    {"tetrahydrofuran", "tetrahydrofuran"}, {"thf", "tetrahydrofuran"},
    {"tetrahydrothiophenedioxide", "tetrahydrothiophenedioxide"}, {"sulfolan", "tetrahydrothiophenedioxide"},
    {"tetralin", "tetralin"}, {"1,2,3,4-tetrahydronaphthalene", "tetralin"}, {"toluene", "toluene"},
    {"tributylphosphate", "tributylphosphate"}, {"triethylamine", "triethylamine"},
    {"trimethylbenzene", "trimethylbenzene"}, {"1,2,3-trimethylbenzene", "trimethylbenzene"},
    {"undecane", "undecane"}, {"water", "water"}, {"h2o", "water"}, {"xylene", "xylene"}, {"methanol", "methanol"},
};
// data.f90 density (kg/m^3); "f" marks a default-real literal (single precision in the Fortran source)
const std::pair<const char*, double> kDensity[] = {
    {"2methylpyridine", 944.0}, {"4methyl2pentanone", 802.0}, {"aceticacid", 1051.0}, {"acetonitrile", 787.0},
    {"acetophenone", 1028.0}, {"aniline", 1022.0}, {"anisole", 995.6}, {"benzene", 879.0}, {"benzonitrile", 1009.0},
    {"benzylalcohol", 1041.0}, {"bromobenzene", 1490.0}, {"bromoethane", 1460.0}, {"bromoform", 2891.2},
    {"bromooctane", 1118.0}, {"butanol", 810.0}, {"butanone", 806.0}, {"butylacetate", 921.2},
    {"butylbenzene", 860.1}, {"carbondisulfide", 1263.2}, {"carbontet", 1590.0}, {"chlorobenzene", 1105.8},
    {"chloroform", 1483.2}, {"chlorohexane", 879.0}, {"cyclohexane", 779.0}, {"cyclohexanone", 945.0},
    {"decalin", 890.0}, {"decane", 730.0}, {"decanol", 840.0}, {"dibromoethane", 2055.5}, {"dibutylether", 767.0},
    {"dichloroethane", 1253.0}, {"diethylether", 714.0}, {"diisopropylether", 724.0}, {"dimethylacetamide", 943.0},
    {"dimethylformamide", 950.0}, {"dimethylpyridine", 945.0}, {"dimethylsulfoxide", 1101.0}, {"dodecane", 748.7},
    {"ethanol", 790.0}, {"ethoxybenzene", 967.0}, {"ethylacetate", 902.0}, {"ethylbenzene", 867.0},
    {"fluorobenzene", 1022.5}, {"fluoroctane", 814.0}, {"heptane", 683.8}, {"heptanol", 822.0}, {"hexadecane", 773.2},
    {"hexadecyliodide", 1121.0}, {"hexane", 659.0}, {"hexanol", 850.0}, {"iodobenzene", 1808.0},
    {"isobutanol", 802.0}, {"isooctane", 691.9}, {"isopropanol", 785.0}, {"isopropylbenzene", 866.0},
    {"isopropyltoluene", 857.0}, {"mcresol", 1033.6}, {"mesitylene", 860.0}, {"methoxyethanol", 966.0},
    {"methylenechloride", 1322.0}, {"methylformamide", 1011.0}, {"nitrobenzene", 1203.7}, {"nitroethane", 1050.0},
    {"nitromethane", 1139.0}, {"nonane", 718.0}, {"nonanol", 827.0}, {"octane", 703.0}, {"octanol", 829.0},
    {"odichlorobenzene", 1306.0}, {"onitrotoluene", 1162.2}, {"pentadecane", 768.5}, {"pentane", 626.0},
    {"pentanol", 818.0}, {"perfluorobenzene", 1610.0}, {"phenylether", 1070.0}, {"propanol", 803.0},
    {"pyridine", 983.0}, {"secbutanol", 810.0}, {"secbutylbenzene", 858.0}, {"tbutylbenzene", 866.9},
    {"tetrachloroethene", 1630.0}, {"tetrahydrofuran", 888.0}, {"tetrahydrothiophenedioxide", 1260.0},
    {"tetralin", 974.0}, {"toluene", 867.0}, {"tributylphosphate", 982.0}, {"triethylamine", 729.0},
    {"trimethylbenzene", 894.4}, {"undecane", 740.1}, {"water", 998.0}, {"xylene", 880.0}, {"methanol", 792.0},
};
// numsa 0.2.0 smd/init.f90 init_smd_def: n, alpha, beta, macroscopic surface tension, aromatic fraction, F/Cl/Br fraction
struct SmdSolventProps { const char* name; double n, alpha, beta, msurft, arom, fclbr; };
const SmdSolventProps kSmdProps[] = {
    {"methanol", 1.3288, 0.43, 0.47, 31.77, 0.0, 0.0}, {"2methylpyridine", 1.4957, 0.0, 0.58, 47.5, 0.509796, 0.0},
    {"4methyl2pentanone", 1.3962, 0.0, 0.51, 33.83, 0.0, 0.0}, {"aceticacid", 1.372, 0.61, 0.44, 39.01, 0.0, 0.0},
    {"acetonitrile", 1.3442, 0.07, 0.32, 41.25, 0.0, 0.0}, {"acetophenone", 1.5372, 0.0, 0.48, 56.19, 0.444889, 0.0},
    {"aniline", 1.5863, 0.26, 0.41, 60.62, 0.734449, 0.0}, {"anisole", 1.5174, 0.0, 0.29, 50.52, 0.5625, 0.0},
    {"benzene", 1.5011, 0.0, 0.14, 40.62, 1.0, 0.0}, {"benzylalcohol", 1.5396, 0.33, 0.56, 52.96, 0.5625, 0.0},
    {"bromobenzene", 1.5597, 0.0, 0.09, 50.72, 0.734449, 0.020449},
    {"bromoethane", 1.4239, 0.0, 0.12, 34.0, 0.0, 0.110889}, {"bromoform", 1.6005, 0.15, 0.06, 64.58, 0.0, 0.5625},
    {"bromooctane", 1.4524, 0.0, 0.12, 41.28, 0.0, 0.0121}, {"butanol", 1.3993, 0.37, 0.48, 35.88, 0.0, 0.0},
    {"butanone", 1.3788, 0.0, 0.51, 34.5, 0.0, 0.0}, {"butylacetate", 1.3941, 0.0, 0.45, 35.81, 0.0, 0.0},
    {"butylbenzene", 1.4898, 0.0, 0.15, 41.33, 0.36, 0.0}, {"carbondisulfide", 1.6319, 0.0, 0.07, 45.45, 0.0, 0.0},
    {"carbontet", 1.4601, 0.0, 0.0, 38.04, 0.0, 0.64},
    {"chlorobenzene", 1.5241, 0.0, 0.07, 47.48, 0.734449, 0.020449},
    {"chloroform", 1.4459, 0.15, 0.02, 38.39, 0.0, 0.5625}, {"chlorohexane", 1.4199, 0.0, 0.1, 37.03, 0.0, 0.020449},
    {"cyclohexane", 1.4266, 0.0, 0.0, 35.48, 0.0, 0.0}, {"cyclohexanone", 1.4507, 0.0, 0.56, 49.76, 0.0, 0.0},
    {"decalin", 1.4528, 0.0, 0.0, 43.82, 0.0, 0.0}, {"decane", 1.4102, 0.0, 0.0, 33.64, 0.0, 0.0},
    {"decanol", 1.4372, 0.37, 0.48, 41.04, 0.0, 0.0}, {"dibromoethane", 1.5387, 0.1, 0.17, 56.93, 0.0, 0.25},
    {"dibutylether", 1.3992, 0.0, 0.45, 35.98, 0.0, 0.0}, {"dichloroethane", 1.4448, 0.1, 0.11, 45.86, 0.0, 0.25},
    {"diethylether", 1.3526, 0.0, 0.41, 23.96, 0.0, 0.0}, {"diisopropylether", 1.3679, 0.0, 0.41, 24.86, 0.0, 0.0},
    {"dimethylacetamide", 1.438, 0.0, 0.78, 47.62, 0.0, 0.0},
    {"dimethylformamide", 1.4305, 0.0, 0.74, 49.56, 0.0, 0.0},
    {"dimethylpyridine", 1.4953, 0.0, 0.63, 44.64, 0.390625, 0.0},
    {"dimethylsulfoxide", 1.417, 0.0, 0.88, 61.78, 0.0, 0.0}, {"dmso", 1.417, 0.0, 0.88, 61.78, 0.0, 0.0},
    {"DMSO", 1.417, 0.0, 0.88, 61.78, 0.0, 0.0}, {"dodecane", 1.4216, 0.0, 0.0, 35.85, 0.0, 0.0},
    {"ethanol", 1.3611, 0.37, 0.48, 31.62, 0.0, 0.0}, {"ethoxybenzene", 1.5076, 0.0, 0.32, 46.65, 0.444889, 0.0},
    {"ethylacetate", 1.3723, 0.0, 0.45, 33.67, 0.0, 0.0}, {"ethylbenzene", 1.4959, 0.0, 0.15, 41.38, 0.5625, 0.0},
    {"fluorobenzene", 1.4684, 0.0, 0.1, 38.37, 0.734449, 0.020449},
    {"fluoroctane", 1.3935, 0.0, 0.1, 33.92, 0.0, 0.012321}, {"heptane", 1.3878, 0.0, 0.0, 28.28, 0.0, 0.0},
    {"hexadecane", 1.4345, 0.0, 0.0, 38.93, 0.0, 0.0}, {"hexadecyliodide", 1.4806, 0.0, 0.15, 46.48, 0.0, 0.0},
    {"hexane", 1.3749, 0.0, 0.0, 25.7495, 0.0, 0.0}, {"hexanol", 1.4178, 0.37, 0.48, 37.15, 0.0, 0.0},
    {"iodobenzene", 1.62, 0.0, 0.12, 55.72, 0.734449, 0.0}, {"isobutanol", 1.3955, 0.37, 0.48, 32.38, 0.0, 0.0},
    {"isooctane", 1.3915, 0.0, 0.0, 26.38, 0.0, 0.0}, {"isopropanol", 1.3776, 0.33, 0.56, 30.13, 0.0, 0.0},
    {"isopropylbenzene", 1.4915, 0.0, 0.16, 39.85, 0.444889, 0.0},
    {"isopropyltoluene", 1.4909, 0.0, 0.19, 38.34, 0.36, 0.0}, {"mcresol", 1.5438, 0.57, 0.34, 51.37, 0.5625, 0.0},
    {"mesitylene", 1.4994, 0.0, 0.19, 39.65, 0.444889, 0.0}, {"methoxyethanol", 1.4024, 0.3, 0.84, 44.39, 0.0, 0.0},
    {"methylenechloride", 1.4242, 0.1, 0.05, 39.15, 0.0, 0.444889},
    {"methylformamide", 1.4319, 0.4, 0.55, 55.4372, 0.0, 0.0},
    {"nitrobenzene", 1.5562, 0.0, 0.28, 57.54, 0.444889, 0.0}, {"nitroethane", 1.3917, 0.02, 0.33, 46.25, 0.0, 0.0},
    {"nitromethane", 1.3817, 0.06, 0.31, 52.58, 0.0, 0.0}, {"nonane", 1.4054, 0.0, 0.0, 32.21, 0.0, 0.0},
    {"nonanol", 1.4333, 0.37, 0.48, 40.14, 0.0, 0.0}, {"octane", 1.3974, 0.0, 0.0, 30.4273, 0.0, 0.0},
    {"octanol", 1.4295, 0.37, 0.48, 39.01, 0.0, 0.0}, {"odichlorobenzene", 1.5515, 0.0, 0.04, 52.72, 0.5625, 0.0625},
    {"onitrotoluene", 1.545, 0.0, 0.27, 59.12, 0.36, 0.0}, {"pentadecane", 1.4315, 0.0, 0.0, 38.34, 0.0, 0.0},
    {"pentane", 1.3575, 0.0, 0.0, 22.3, 0.0, 0.0}, {"pentanol", 1.4101, 0.37, 0.48, 36.5, 0.0, 0.0},
    {"perfluorobenzene", 1.3777, 0.0, 0.0, 31.74, 0.25, 0.25}, {"phenylether", 1.5787, 0.0, 0.2, 38.5, 0.851929, 0.0},
    {"propanol", 1.385, 0.37, 0.48, 33.57, 0.0, 0.0}, {"pyridine", 1.5095, 0.0, 0.52, 52.62, 0.693889, 0.0},
    {"secbutanol", 1.3978, 0.33, 0.56, 32.44, 0.0, 0.0}, {"secbutylbenzene", 1.4895, 0.0, 0.16, 40.35, 0.36, 0.0},
    {"tbutylbenzene", 1.4927, 0.0, 0.16, 39.78, 0.36, 0.0},
    {"tetrachloroethene", 1.5053, 0.0, 0.0, 45.19, 0.0, 0.444889},
    {"tetrahydrofuran", 1.405, 0.0, 0.48, 39.44, 0.0, 0.0},
    {"tetrahydrothiophenedioxide", 1.4833, 0.0, 0.88, 87.49, 0.0, 0.0},
    {"tetralin", 1.5413, 0.0, 0.19, 47.74, 0.36, 0.0}, {"toluene", 1.4961, 0.0, 0.14, 40.2, 0.734449, 0.0},
    {"tributylphosphate", 1.4224, 0.0, 1.21, 27.55, 0.0, 0.0}, {"triethylamine", 1.401, 0.0, 0.79, 29.1, 0.0, 0.0},
    {"trimethylbenzene", 1.5048, 0.0, 0.19, 42.03, 0.444889, 0.0}, {"undecane", 1.4398, 0.0, 0.0, 34.85, 0.0, 0.0},
    {"xylene", 1.4995, 0.0, 0.16, 41.38, 0.5625, 0.0},
};


/// data.f90 atomicMassNist (Z = 1..118).
const double kAtomicMass[118] = {
    1.00794075, 4.00260193, 6.94003660, 9.01218307, 10.81102805, 12.01073590, 14.00670321, 15.99940492,
    18.99840316, 20.18004638, 22.98976928, 24.30505162, 26.98153853, 28.08549871, 30.97376200, 32.06478741,
    35.45293758, 39.94779856, 39.09830091, 40.07802251, 44.95590828, 47.86674496, 50.94146504, 51.99613176,
    54.93804391, 55.84514443, 58.93319429, 58.69334711, 63.54603995, 65.37778253, 69.72306607, 72.62755016,
    74.92159457, 78.95938856, 79.90352778, 83.79800000, 85.46766360, 87.61664447, 88.90584030, 91.22364160,
    92.90637300, 95.95978854, 97.90721240, 101.06494014, 102.90549800, 106.41532751, 107.86814963, 112.41155782,
    114.81808663, 118.71011259, 121.75978367, 127.60312648, 126.90447190, 131.29276145, 132.90545196, 137.32689163,
    138.90546887, 140.11573074, 140.90765760, 144.24159603, 144.91275590, 150.36635571, 151.96437813, 157.25213065,
    158.92535470, 162.49947282, 164.93032880, 167.25908265, 168.93421790, 173.05415017, 174.96681496, 178.48497872,
    180.94787564, 183.84177755, 186.20670455, 190.22485963, 192.21605165, 195.08445686, 196.96656879, 200.59916703,
    204.38341284, 207.21690806, 208.98039910, 208.98243080, 209.98714790, 222.01757820, 223.01973600, 226.02541030,
    227.02775230, 232.03805580, 231.03588420, 238.02891046, 237.04817360, 244.06420530, 243.06138130, 247.07035410,
    247.07030730, 251.07958860, 252.08298000, 257.09510610, 258.09843150, 259.10103000, 262.10961000, 267.12179000,
    269.12791000, 271.13393000, 270.13336000, 276.14846000, 276.15159000, 280.16131000, 282.16912000, 284.17416000,
    284.17873000, 289.19042000, 288.19274000, 293.20449000, 292.20746000, 294.21392000,
};

std::string toLower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

/// Fortran default-real (single precision) literal: decimal -> float (correctly rounded), then widened to double.
double singleLiteral(const std::string& token) { return static_cast<double>(std::strtof(token.c_str(), nullptr)); }

/// Round as a Fortran Fw.d edit descriptor writes (nearest, on the exact binary value) and read the text back.
double roundFixed(double x, int decimals)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", decimals, x);
    return std::strtod(buf, nullptr);
}

} // namespace

// ---------------------------------------------------------------------------------------------------------------
// internal data used by the other cpcmx_*.cpp files (cpcmx_internal.h)
// ---------------------------------------------------------------------------------------------------------------

double solventDensity(const std::string& name)
{
    for (const auto& d : kDensity)
        if (name == d.first)
            return d.second;
    return -1.0;
}

double atomicMass(int z) { return (z > 0 && z <= 118) ? kAtomicMass[z - 1] : -1.0; }

/// initialize.f90 setup_cov: hard-coded covalent radii (Angstrom) for H..Ti; every other element gets 0
/// (dict_get_key returns DICT_NULL = 0 for a missing key).
double covalentRadius(int z)
{
    static const double r[22] = { 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66,
                                  1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76, 1.70, 1.60 };
    return (z >= 1 && z <= 22) ? r[z - 1] : 0.0;
}

bool smdSolventProperties(const std::string& name, double& n, double& alpha, double& beta, double& msurft,
                          double& arom, double& fclbr)
{
    for (const auto& p : kSmdProps)
        if (name == p.name) {
            n = p.n;
            alpha = p.alpha;
            beta = p.beta;
            msurft = p.msurft;
            arom = p.arom;
            fclbr = p.fclbr;
            return true;
        }
    return false;
}

/// internaldb.f90 xtb_water_smd: the SMD parameter text CPCM-X writes to ./smd_h2o and numsa reads back.
const std::vector<std::string>& smdTextWater()
{
    static const std::vector<std::string> t = {
        "#Zk", "H \t-77.70436564", "C \t-44.32921125", "N \t-25.97069034", "O \t-50.70232515", "F \t-6.43754644",
        "Cl\t-36.05536752", "Br\t-20.79864972", "S \t-18.13699434",
        "#Zkk Matrix",
        "H\t00.00\t76.79778940\t00.00\t00.00\t00.00",
        "C\t00.00\t-1.94928804\t00.00\t00.00\t00.00",
        "O\t00.00\t12.35744639\t00.00\t165.85557924\t12.23268205",
        "N\t00.00\t-24.36906946\t00.00\t00.00\t00.00",
        "P\t00.00\t00.00\t00.00\t00.00\t00.00",
        "#rzkk Matrix",
        "H  00.00  1.84000000  00.00  1.55000000  00.00  00.00  00.00  00.00  00.00  00.00",
        "C  1.55000000  1.84000000  1.84000000  1.84000000  1.84000000  2.20000000  2.20000000  2.10000000  2.30000000  2.60000000",
        "N  00.00  1.84000000  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "O  00.00  1.33000000  1.50000000  1.80000000  00.00  2.10000000  00.00  00.00  00.00  00.00",
        "F  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "P  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "S  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "Cl  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "Br  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "I  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "#drzkk Matrix",
        "H  00.00  0.30000000  00.00  0.30000000  00.00  00.00  00.00  00.00  00.00  00.00",
        "C  0.30000000  0.30000000  0.30000000  0.30000000  0.30000000  0.30000000  0.30000000  0.30000000  0.30000000  0.30000000",
        "N  00.00  0.30000000  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "O  00.00  0.10000000  0.30000000  0.30000000  0.30000000  00.00  00.00  00.00  00.00  00.00",
        "F  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "P  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "S  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "Cl  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "Br  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "I  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "#NC3", "33.62391370", "1.22500000", "0.06500000"
    };
    return t;
}

/// internaldb.f90 xtb_other_smd: the SMD parameter text CPCM-X writes to ./smd_ot for all other solvents.
const std::vector<std::string>& smdTextOther()
{
    static const std::vector<std::string> t = {
        "#Zk", "C\t-3.70891582", "O\t38.16986455", "N\t11.67657766", "Cl\t-26.23511892", "Br\t-19.44296237",
        "S\t-15.85282828", "H\t-5.61002654", "P\t-31.18737493",
        "#Zk2", "C\t8.80543594", "O\t-2.68778654",
        "#Zk3", "C\t-5.80252796", "O\t-87.53912863",
        "#Zkk Matrix",
        "H\t00.00\t-6.10646558\t-72.27657574\t00.00",
        "C\t00.00\t-15.05405921\t00.00\t-16.54455408",
        "O\t00.00\t-20.55673110\t00.00\t00.00",
        "N\t00.00\t00.00\t00.00\t00.00",
        "#Zkk2 Matrix",
        "O\t00.00\t65.41437683\t00.00",
        "C\t00.00\t00.00\t-2.81967860",
        "N\t00.00\t-32.31978386\t00.00",
        "#Zkk3 Matrix",
        "O\t-40.24403232\t67.09645304",
        "N\t00.00\t00.00",
        "#rzkk Matrix",
        "H  00.00  1.55100494  00.00  1.53632796  00.00  00.00  00.00  00.00  00.00  00.00",
        "C  1.54582282  1.83993169  1.83959886  2.07175304  1.84000000  2.24400000  2.20000000  2.10000000  2.30000000  2.60000000",
        "N  00.00  1.84000000  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "O  00.00  1.33181762  1.50000000  1.80000000  00.00  2.10000000  00.00  00.00  00.00  00.00",
        "F  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "P  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "S  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "Cl  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "Br  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "I  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "#drzkk Matrix",
        "H  00.00  0.30023238  00.00  0.30023238  00.00  00.00  00.00  00.00  00.00  00.00",
        "C  0.30023238  0.30023238  0.30023238  0.30023238  0.30023238  0.30023238  0.30023238  0.30023238  0.30023238  0.30023238",
        "N  00.00  0.30023238  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "O  00.00  0.11286083  0.30023238  0.30023238  0.30023238  00.00  00.00  00.00  00.00  00.00",
        "F  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "P  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "S  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "Cl  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "Br  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "I  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00  00.00",
        "#NC3", "00.00", "1.22500000", "0.06500000",
        "#Solvent", "0.08563085", "-7.21289853", "-14.52469035", "-5.13736584"
    };
    return t;
}

// ---------------------------------------------------------------------------------------------------------------
// public API
// ---------------------------------------------------------------------------------------------------------------

std::string solventName(const std::string& solvent)
{
    const std::string s = toLower(solvent);
    for (const auto& p : kSolventNames)
        if (s == p.first)
            return p.second;
    return "";
}

std::string knownSolvents()
{
    std::vector<std::string> names;
    for (const auto& p : kSolventNames)
        if (std::find(names.begin(), names.end(), p.second) == names.end())
            names.push_back(p.second);
    std::string out;
    for (const auto& n : names)
        out += (out.empty() ? "" : ", ") + n;
    return out;
}

/// internaldb.f90 load_param("xtb", solvent): water/h2o -> xtb_water, every other solvent -> xtb_other.
Parameters xtbParameters(const std::string& solvent_name)
{
    if (solvent_name == "water" || solvent_name == "h2o")
        return { 0.30880726, 7540.06830745, -1.26119686, 15813.08917624, 0.00845974, 16.58382447, -0.23040291,
                 -0.13909842, -10.60951107, 0.0 };
    return { 0.53725564, 1977.88495090, 1.61721497, 2664.61735599, 0.00806820, 57.11402448, -0.12356319, -0.08426734,
             -8.41974122, 0.0 };
}

std::string findDatabase(const std::string& configured, std::vector<std::string>* tried)
{
    namespace fs = std::filesystem;
    std::vector<std::string> dirs;
    if (!configured.empty() && configured != "auto") {
        dirs.push_back(configured);
    } else {
        if (const char* e = std::getenv("CURCUMA_CPCMX_DB"); e && *e)
            dirs.emplace_back(e);
#if defined(__linux__)
        char buf[PATH_MAX];
        const ssize_t n = ::readlink("/proc/self/exe", buf, sizeof(buf) - 1);
        if (n > 0) {
            buf[n] = '\0';
            const fs::path exe_dir = fs::path(buf).parent_path();
            dirs.push_back((exe_dir / "cpcmx_db").string());
            dirs.push_back((exe_dir / ".." / "share" / "curcuma" / "cpcmx_db").string());
        }
#endif
    }
    for (const auto& d : dirs) {
        if (tried)
            tried->push_back(d);
        std::error_code ec;
        if (fs::is_directory(d, ec))
            return d;
    }
    return "";
}

/// The solvent side of internaldb.f90 internalcosmo. The arrays of DB/xtb/<name>.fh are the columns of
/// DB/xtb/<name>.cosmo written as default-real literals (checked for all 88 files), so every number is single
/// precision: sid = atom, sxyz = btoa * x(Bohr), sarea = area, su = charge/area, spot = potential (no conversion),
/// axyz = btoa * atom coordinate (Bohr), aelement = element.
bool loadSolvent(const std::string& db_dir, const std::string& solvent, Molecule& out, std::string& error)
{
    const std::string name = solventName(solvent);
    if (name.empty()) {
        error = solvent + " is not a valid CPCM-X solvent (known: " + knownSolvents() + ")";
        return false;
    }
    const std::string path = (std::filesystem::path(db_dir) / (name + ".cosmo")).string();
    std::ifstream in(path);
    if (!in) {
        error = "Solvent " + name + " was not found in the CPCM-X database (missing file " + path + ")";
        return false;
    }
    std::vector<std::string> lines;
    for (std::string l; std::getline(in, l);)
        lines.push_back(l);
    out = Molecule{};
    size_t k = 0;
    while (k < lines.size() && lines[k].rfind("#atom", 0) != 0)
        ++k;
    std::vector<std::array<double, 3>> axyz;
    for (++k; k < lines.size() && lines[k].find('$') == std::string::npos; ++k) {
        std::istringstream ls(lines[k]);
        std::string idx, x, y, z, el;
        if (!(ls >> idx >> x >> y >> z >> el))
            continue;
        axyz.push_back({ btoa * singleLiteral(x), btoa * singleLiteral(y), btoa * singleLiteral(z) });
        out.z.push_back(Elements::String2Element(el));
    }
    while (k < lines.size() && lines[k] != "$segment_information")
        ++k;
    std::vector<std::array<double, 6>> seg;
    for (++k; k < lines.size(); ++k) {
        if (lines[k].empty() || lines[k][0] == '#')
            continue;
        std::istringstream ls(lines[k]);
        std::string n, at, x, y, z, q, a, s, p;
        if (!(ls >> n >> at >> x >> y >> z >> q >> a >> s >> p))
            continue;
        out.id.push_back(std::stoi(at));
        seg.push_back({ btoa * singleLiteral(x), btoa * singleLiteral(y), btoa * singleLiteral(z), singleLiteral(a),
                        singleLiteral(s), singleLiteral(p) });
    }
    if (axyz.empty() || seg.empty()) {
        error = "CPCM-X database file " + path + " has no atoms or no segments";
        return false;
    }
    const int nat = static_cast<int>(axyz.size());
    out.atom_xyz.resize(nat, 3);
    for (int i = 0; i < nat; ++i)
        for (int c = 0; c < 3; ++c)
            out.atom_xyz(i, c) = axyz[i][c];
    const int ns = static_cast<int>(seg.size());
    out.xyz.resize(ns, 3);
    out.area.resize(ns);
    out.su.resize(ns);
    out.pot.resize(ns);
    for (int i = 0; i < ns; ++i) {
        for (int c = 0; c < 3; ++c)
            out.xyz(i, c) = seg[i][c];
        out.area(i) = seg[i][3];
        out.su(i) = seg[i][4];
        out.pot(i) = seg[i][5];
    }
    return true;
}

/// read_cosmo (initialize.f90) applied to the xtb.cosmo file that writeCosmoFile (xtb cosmo.f90) would write:
///   segment line '(2i5, 7(1x, f14.9))': n, atom, x, y, z (Bohr), charge, area, charge/area, phi/autoaa;
///   read_cosmo keeps area and charge/area and multiplies position and potential by BtoA (single-precision 0.529177);
///   atom line '(i4, 3(1x, f19.14), ...)' in Bohr, multiplied by btoa; energy line f21.10.
Molecule soluteFromConductor(const ConductorData& d, bool reference_rounding)
{
    auto r = [reference_rounding](double x, int dec) { return reference_rounding ? roundFixed(x, dec) : x; };
    const double autoaa_xtb = ALPBParameters::autoaa;   // xtb mctc convert.f90: 0.52917726
    Molecule m;
    const int ns = static_cast<int>(d.segment_atom.size());
    m.id.resize(ns);
    m.xyz.resize(ns, 3);
    m.area.resize(ns);
    m.su.resize(ns);
    m.pot.resize(ns);
    for (int i = 0; i < ns; ++i) {
        m.id[i] = d.segment_atom[i] + 1;
        for (int c = 0; c < 3; ++c)
            m.xyz(i, c) = r(d.segment_xyz(c, i), 9) * btoa;
        m.area(i) = r(d.area(i), 9);
        m.su(i) = r(d.charge(i) / d.area(i), 9);
        m.pot(i) = r(d.potential(i) / autoaa_xtb, 9) * btoa;
    }
    const int nat = static_cast<int>(d.z.size());
    m.z = d.z;
    m.atom_xyz.resize(nat, 3);
    for (int i = 0; i < nat; ++i)
        for (int c = 0; c < 3; ++c)
            m.atom_xyz(i, c) = r(d.atom_xyz(i, c), 14) * btoa;
    m.energy = r(d.energy, 10);
    return m;
}

} // namespace curcuma::solvation::cpcmx
