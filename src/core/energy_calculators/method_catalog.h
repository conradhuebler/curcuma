/*
 * <Method catalog: which parameters each energy method reads, per-method help>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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

/**
 * @file method_catalog.h
 * @brief One table that says, per energy method, what it is and which parameters it reads.
 *
 * Claude Generated (Oct 31, 2026), WP "Methodenhilfe vereinheitlichen".
 *
 * @par Design: help and parameters are organised per method, not per implementation family
 * Operator decision of 2026-10-31. A user thinks in methods (gfn2, gfnff, b2plyp, hf-3c, uff),
 * not in the code family that implements them (there is no "qm" vs "sqm" help category). Rules:
 *  - ONE source: kMethodCatalog (method_catalog.cpp) declares every MethodFactory name with its
 *    implementation, status, gradient, GPU, restrictions and the parameter groups it reads.
 *    `curcuma -help <method>`, `curcuma -help methods`, the warning about parameters the selected
 *    method does not read, and the flat-flag routing (MethodFactory::parameterScopeForMethod,
 *    readsSharedScfScope) are all derived from it. Nothing else may keep a second list.
 *  - Warning, never an error: a parameter that is unknown or not read by the selected method is
 *    reported at run time (CurcumaLogger::warn); existing command lines must not break.
 *  - Every new method and every new parameter of a method module must be declared here, so that
 *    it is reachable through `-help <method>` and its flat name. `curcuma -check_method_catalog`
 *    (ctest `method_catalog_consistency`) fails when a factory name or a registered parameter of
 *    a method module is missing.
 *  - The shared SCF knobs (module scf: -scf_threshold, -scf_max_iterations, ...) mean the same for
 *    every method that reads them; where a method interprets one differently the group note and
 *    docs/METHODS_OVERVIEW.md (table knob x method) say so. Method-scoped spellings (-qm.scf_*,
 *    -xtb.*) stay valid as aliases with higher precedence.
 * See docs/METHODS_OVERVIEW.md.
 */

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "external/json.hpp"

/// A set of registry parameters that a method reads, shown under one topic in `-help <method>`.
/// Claude Generated (Oct 2026).
struct MethodParamGroup {
    std::string id;                    ///< key used in MethodCatalogEntry::groups, e.g. "qm_grid"
    std::string topic;                 ///< heading in the help ("SCF", "Grid / DFT", "PT2" ...); empty = registry category
    std::string module;                ///< registry module the names belong to
    std::vector<std::string> names;    ///< explicit parameter names; empty = every parameter of the module
    std::vector<std::string> except;   ///< removed from a whole-module group
    std::string note;                  ///< one line printed under the topic (when it applies, how it differs)
    ~MethodParamGroup();               ///< out of line (method_catalog.cpp): avoids -Winline noise
};

/// One row of the catalog: a method (and its variants/aliases) as the user sees it.
/// Claude Generated (Oct 2026).
struct MethodCatalogEntry {
    std::vector<std::string> names;    ///< MethodFactory names covered by this row (lower case, primary first)
    std::string family;                ///< e.g. "native xTB", "native ab-initio QM", "external: Ulysses"
    std::string implementation;        ///< where the code lives: "native (curcuma)", "TBLite", "xtb binary", "Ulysses", "ORCA" ...
    std::string status;                ///< "native, AI-generated, machine-tested" or "external library, curcuma wrapper"
    std::string gradient;              ///< "analytic", "none", ...
    std::string gpu;                   ///< "cuda, rocm, vulkan" / "cuda" / "no"
    std::string summary;               ///< one sentence
    std::string restrictions;          ///< closed shell, element range, ...
    std::string validation;            ///< documentation of the validation (path relative to the repo)
    std::string route_scope;           ///< registry module that wins for a flat flag several modules own ("" = none)
    std::vector<std::string> groups;   ///< ids of MethodParamGroup the method reads
    std::vector<std::pair<std::string, std::string>> variants;  ///< name -> what differs from the primary name
    ~MethodCatalogEntry();             ///< out of line (method_catalog.cpp): avoids -Winline noise
};

/// Shared SCF knob x method: how one knob acts on one method family (docs table and help).
struct ScfKnobRow {
    std::string knob;                   ///< flat flag name, e.g. "scf_threshold"
    std::vector<std::string> effect;    ///< one cell per ScfKnobColumns() entry
};

/// The method catalog. All functions are static; the data lives in method_catalog.cpp.
/// Claude Generated (Oct 2026).
class MethodCatalog {
public:
    static const std::vector<MethodCatalogEntry>& entries();
    static const std::vector<MethodParamGroup>& groups();
    /// Row for a method name (case-insensitive), nullptr if unknown.
    static const MethodCatalogEntry* find(const std::string& method_name);
    static const MethodParamGroup* findGroup(const std::string& id);

    /// Expanded list (module, canonical name) of every parameter @p entry reads, in group order.
    static std::vector<std::pair<std::string, std::string>> parametersOf(const MethodCatalogEntry& entry);
    /// True when the method of @p entry reads parameter @p name (canonical or alias) of @p module.
    /// The legacy scope "dft" counts as "qm", merged modules via ParameterRegistry::canonicalModule().
    static bool reads(const MethodCatalogEntry& entry, const std::string& module, const std::string& name);

    /// Module that wins for a flat flag registered in several modules ("" = none).
    static std::string routeScope(const std::string& method_name);
    /// True when the method reads the shared registry module "scf".
    static bool readsSharedScf(const std::string& method_name);

    /// `curcuma -help <method>`. Returns false when @p name is not a method.
    static bool printMethodHelp(const std::string& name);
    /// `curcuma -help methods`: every factory name with family, status, gradient, GPU, one sentence.
    static void printMethodsOverview();
    /// The shared SCF knob table (columns = method families).
    static const std::vector<std::string>& scfKnobColumns();
    static const std::vector<ScfKnobRow>& scfKnobTable();

    /// Warnings for parameters in @p controller that the selected method(s) do not read or that are
    /// not registered at all. @p methods are the method names the run uses (-method, -md_method, ...),
    /// @p command_module the module of the CLI command. Never throws; returns the messages.
    static std::vector<std::string> unreadParameterWarnings(const nlohmann::json& controller,
                                                            const std::vector<std::string>& methods,
                                                            const std::string& command_module);

    /// Consistency check for ctest: every MethodFactory name has exactly one row, every row name is a
    /// factory name, every group parameter is registered, every registered parameter of a method module
    /// is read by some row or listed with a reason in the unclaimed list. Empty = consistent.
    static std::vector<std::string> consistencyProblems();
};
