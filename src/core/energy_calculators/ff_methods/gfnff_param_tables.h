/*
 * < Runtime-overridable GFN-FF parameter tables >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Sep 2026, rev-gfnff WP1a)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

/**
 * @file gfnff_param_tables.h
 * @brief Per-instance copies of the GFN-FF parameter tables that a fit may change.
 *
 * gfnff_par.h holds the published GFN-FF parameters as compile-time constants and stays
 * the source of truth. A parameter fit (rev-gfnff) needs to evaluate the force field with
 * modified numbers without recompiling, so GFNFFTables carries a RUNTIME copy of the tables
 * that stages 1-3 of the roadmap touch: the "gen%" scalars that used to be function-local
 * literals in gfnff_method.cpp, and the element-indexed EEQ, bond, angle, repulsion and
 * torsion tables. Every GFNFF instance owns a shared pointer to one; the pristine defaults
 * are built once from gfnff_par.h and shared, an override is a deep-merged copy.
 *
 * Override JSON (sparse, flat):
 * @code
 * {"gen":    {"rabshift": -0.11, "srb1": 0.3731},
 *  "tables": {"bond_params": {"6": 0.385, "7": 0.379},      // keys are atomic numbers Z
 *             "repa": [ ... full array ... ]},                // or the whole table
 *  "rev":    {"bo_center": 2.0, "p_over": {"1": 0.05}}}      // free-form, read by the rev stages
 * @endcode
 * Unknown keys are an error (fail loud), so a typo cannot silently fit nothing.
 *
 * Element convention: every table here is indexed [Z-1] exactly like its gfnff_par.h
 * original (bond_params[z-1], chi_eeq[z-1], repa[z-1], ...).
 */

#include "json.hpp"

#include <array>
#include <memory>
#include <string>
#include <vector>

using json = nlohmann::json;

/// Global ("gen%") scalars of GFN-FF that were literals in the generator code.
struct GFNFFGen {
    // bond r0 shifts (gfnff_ini.f90:1143-1149, 1268-1276)
    double rabshift = -0.110;
    double rabshifth = -0.050;
    double hyper_shift = 0.030;
    double hshift3 = -0.110;
    double hshift4 = -0.110;
    double hshift5 = -0.060;
    // Hueckel pi corrections of the bond term (gfnff_ini.f90:1174, 1285) and the Hueckel diagonal
    double hueckelp = 0.340;
    double bzref = 0.370;
    double hueckelp2 = 1.00;
    double bzref2 = 0.315;
    double hueckelp3 = -0.24;
    // bond charge factor
    double qfacbm0 = 0.047;
    std::array<double, 5> qfacbm { 1.0, -0.2, -0.2, 0.70, 0.50 };
    // bond exponent alpha = srb1 (1 + fsrb2 dEN^2 + srb3 bstrength)
    double srb1 = 0.3731;
    double srb2 = 0.3171;
    double srb3 = 0.2538;
    // bond strength table (gen%bstren) and the hybridisation matrix built from it
    std::array<double, 9> bstren { 0.0, 1.00, 1.24, 1.98, 1.22, 1.00, 0.78, 3.40, 3.40 };
    double bsmat[4][4] = {
        { 1.0000, 1.3234, 1.0792, 1.0000 },
        { 1.3234, 1.9800, 1.4842, 1.3234 },
        { 1.0792, 1.4842, 1.2400, 1.0792 },
        { 1.0000, 1.3234, 1.0792, 1.0000 }
    };
    // angle / torsion force-constant threshold (gen%fcthr)
    double fcthr = 0.001;
    // 3-/4-body damping cutoffs (gen%atcuta, gen%atcutt, gen%atcutt_nci)
    double atcuta = 0.595;
    double atcutt = 0.505;
    double atcutt_nci = 0.305;
    // hydrogen-bond mixing (gen%hbabmix): p_bh = 1 + hbabmix, p_ab = -hbabmix
    double hbabmix = 0.8;
    // repulsion scalars (gfnff_par.h REPSCALB ...)
    double repscalb = 1.7583;
    double repscaln = 0.4270;
    double qrepscal = 0.3480;
    double nrepscal = -0.1270;
    double hhfac = 0.6290;
    double hh13rep = 1.4580;
    double hh14rep = 0.7080;
    // torsion barrier factors
    double torsf_single = 1.00;
    double torsf_pi = 1.18;
    double torsf_improper = 1.05;
    double torsf_pi_improper = 0.50;
    double torsf_extra_C = -0.90;
    double torsf_extra_N = 0.70;
    double torsf_extra_O = -2.00;
    double fr3 = 0.3;
    double fr4 = 1.0;
    double fr5 = 1.5;
    double fr6 = 5.7;
};

/// Runtime copy of the fit-relevant GFN-FF tables plus the gen scalars.
struct GFNFFTables {
    GFNFFGen gen;
    // EEQ (indexed [Z-1])
    std::vector<double> chi_eeq, gam_eeq, alpha_eeq, cnf_eeq;
    // bond term (indexed [Z-1])
    std::vector<double> bond_params, r0_gfnff, cnfak_gfnff, en_gfnff, en_rab_gfnff;
    // angle term
    std::vector<double> angle_params, angl2_neighbors;
    // repulsion (repa/repan = bonded/non-bonded exponents, repz = effective charges)
    std::vector<double> repa, repan, repz;
    // torsion
    std::vector<double> tors, tors2;
    /// Free-form parameters of the rev-gfnff stages (bond-order switch, over-coordination, charge model).
    json rev = json::object();
    /// Empty for the pristine defaults; otherwise a hash of the override document (goes into the topology fingerprint).
    std::string hash;
    /// The override document this table set was built from (empty object for the defaults).
    json overrides = json::object();

    /// The published GFN-FF values, copied once from gfnff_par.h and shared.
    static std::shared_ptr<const GFNFFTables> defaults();
    /// Defaults deep-merged with a sparse override document; throws std::runtime_error on unknown keys.
    static std::shared_ptr<const GFNFFTables> fromOverrides(const json& overrides);
    /// Every table and scalar as JSON (full arrays, Z-keyed objects are not used on output).
    json toJSON() const;
    /// Names of all overridable "gen" scalars and "tables" arrays (for help output and tests).
    static std::vector<std::string> genNames();
    static std::vector<std::string> tableNames();
};
