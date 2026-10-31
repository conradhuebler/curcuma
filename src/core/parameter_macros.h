/*
 * <Parameter definition macros for automated help and validation>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 *
 */

#pragma once

// Claude Generated: Parameter definition macros for curcuma_param_parser
//
// These macros are used to define parameters directly in capability headers.
// The curcuma_param_parser tool extracts these definitions at build time
// and generates a unified parameter registry (generated/parameter_registry.h).
//
// Usage example - see src/capabilities/analysis.h for real implementation
//
// PARAM macro arguments:
//   1. name          - Canonical parameter name (snake_case recommended)
//   2. type          - ParamType: String, Int, Double, Bool
//   3. default_value - Default value (appropriate for type)
//   4. help_text     - Human-readable description
//   5. category      - Grouping for help display (e.g., "Basic", "Advanced", "Output")
//   6. aliases       - Array-like list of alternative names: {} or {"alias1", "alias2"}
//
// Note: These macros expand to nothing during normal compilation.
//       They only serve as markers for the curcuma_param_parser tool.

#define BEGIN_PARAMETER_DEFINITION(module)
#define PARAM(name, type, default_value, help_text, category, ...) // Variadic to support multiple aliases
#define END_PARAMETER_DEFINITION

// Claude Generated 2026 - Optional trailing annotation on PARAM, and module
// metadata. Both are additive: every existing PARAM line keeps working untouched,
// and a module without MODULE_INFO simply has no description.
//
//   MODULE_INFO("what this module configures", "Category", {"cli", "verbs"})
//
// The CLI verb a module serves was recorded nowhere. Only ten of the thirty-six
// commands share a name with their module -- `md` is `simplemd`, `dock` is
// `docking` -- so anything generating an interface had to guess.
//
// PARAM gains an optional seventh argument: one string of "key=value" pairs
// separated by ';'. Recognised keys:
//
//   tier=primary|advanced|expert   how prominently to offer it (default advanced)
//   enum=a|b|c                     the permitted values
//   unit=K                         the unit, which used to live in the help prose
//   min=0  max=1                   numeric bounds
//   requires=wall_type!=none       when the parameter matters at all
//   deprecated                     flag, no value
//   replaced_by=new_name
//
//   PARAM(thermostat, String, "csvr", "Thermostat type.", "Thermostat", {},
//         "tier=primary; enum=csvr|berendsen|andersen|nosehover|none")
//
// New ParamTypes: StringList and Json (their default is written as JSON text),
// Selection (FragString grammar) and Path.
#define MODULE_INFO(description, category, ...)
