/*
 * <Accuracy Profile implementation>
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

#include "accuracy_profile.h"

#include <algorithm>
#include <cctype>
#include <unordered_map>

// ============================================================================
// Helper: case-insensitive string comparison
// ============================================================================
static std::string toLower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c) { return std::tolower(c); });
    return s;
}

// ============================================================================
// Parse accuracy level
// ============================================================================
AccuracyLevel AccuracyProfile::parseLevel(const std::string& level)
{
    std::string l = toLower(level);
    if (l == "loose") return AccuracyLevel::Loose;
    if (l == "normal") return AccuracyLevel::Normal;
    if (l == "medium") return AccuracyLevel::Medium;
    if (l == "high") return AccuracyLevel::High;
    return AccuracyLevel::Unknown;
}

// ============================================================================
// Set JSON value only if key is absent (user values take precedence)
// ============================================================================
void AccuracyProfile::setIfAbsent(json& config, const std::string& key, const json& value)
{
    if (!config.contains(key)) {
        config[key] = value;
    }
}

// ============================================================================
// Module-specific parameter mappings
// ============================================================================

static void apply_gfnff(json& config, AccuracyLevel level)
{
    switch (level) {
    case AccuracyLevel::Loose:
        AccuracyProfile::setIfAbsent(config, "eeq_max_iterations", 100);
        AccuracyProfile::setIfAbsent(config, "eeq_tolerance", 1e-6);
        AccuracyProfile::setIfAbsent(config, "eeq_accuracy", 1e-4);
        AccuracyProfile::setIfAbsent(config, "eeq_distance_cutoff", 20.0);
        AccuracyProfile::setIfAbsent(config, "cn_cutoff_bohr", 4.0);
        AccuracyProfile::setIfAbsent(config, "cn_accuracy", 0.5);
        AccuracyProfile::setIfAbsent(config, "allow_unconverged_charges", true);
        break;
    case AccuracyLevel::Normal:
        AccuracyProfile::setIfAbsent(config, "eeq_max_iterations", 200);
        AccuracyProfile::setIfAbsent(config, "eeq_tolerance", 1e-8);
        AccuracyProfile::setIfAbsent(config, "eeq_accuracy", 1e-6);
        AccuracyProfile::setIfAbsent(config, "eeq_distance_cutoff", 30.0);
        AccuracyProfile::setIfAbsent(config, "cn_cutoff_bohr", 6.0);
        AccuracyProfile::setIfAbsent(config, "cn_accuracy", 1.0);
        AccuracyProfile::setIfAbsent(config, "allow_unconverged_charges", false);
        break;
    case AccuracyLevel::Medium:
        AccuracyProfile::setIfAbsent(config, "eeq_max_iterations", 500);
        AccuracyProfile::setIfAbsent(config, "eeq_tolerance", 1e-10);
        AccuracyProfile::setIfAbsent(config, "eeq_accuracy", 1e-8);
        AccuracyProfile::setIfAbsent(config, "eeq_distance_cutoff", 40.0);
        AccuracyProfile::setIfAbsent(config, "cn_cutoff_bohr", 8.0);
        AccuracyProfile::setIfAbsent(config, "cn_accuracy", 1.5);
        AccuracyProfile::setIfAbsent(config, "allow_unconverged_charges", false);
        break;
    case AccuracyLevel::High:
        AccuracyProfile::setIfAbsent(config, "eeq_max_iterations", 1000);
        AccuracyProfile::setIfAbsent(config, "eeq_tolerance", 1e-12);
        AccuracyProfile::setIfAbsent(config, "eeq_accuracy", 1e-10);
        AccuracyProfile::setIfAbsent(config, "eeq_distance_cutoff", 50.0);
        AccuracyProfile::setIfAbsent(config, "cn_cutoff_bohr", 10.0);
        AccuracyProfile::setIfAbsent(config, "cn_accuracy", 2.0);
        AccuracyProfile::setIfAbsent(config, "allow_unconverged_charges", false);
        break;
    default:
        break;
    }
}

static void apply_eeq_solver(json& config, AccuracyLevel level)
{
    switch (level) {
    case AccuracyLevel::Loose:
        AccuracyProfile::setIfAbsent(config, "max_pcg_iterations", 100);
        AccuracyProfile::setIfAbsent(config, "pcg_tolerance", 1e-6);
        AccuracyProfile::setIfAbsent(config, "pcg_large_system_tol_factor", 1e-4);
        AccuracyProfile::setIfAbsent(config, "pcg_large_threshold", 100);
        AccuracyProfile::setIfAbsent(config, "pcg_large_system_iterations", 2000);
        AccuracyProfile::setIfAbsent(config, "convergence_threshold", 1e-4);
        AccuracyProfile::setIfAbsent(config, "max_iterations", 20);
        AccuracyProfile::setIfAbsent(config, "eeq_distance_cutoff", 20.0);
        AccuracyProfile::setIfAbsent(config, "allow_unconverged_charges", true);
        break;
    case AccuracyLevel::Normal:
        AccuracyProfile::setIfAbsent(config, "max_pcg_iterations", 100);
        AccuracyProfile::setIfAbsent(config, "pcg_tolerance", 1e-10);
        AccuracyProfile::setIfAbsent(config, "pcg_large_system_tol_factor", 1e-6);
        AccuracyProfile::setIfAbsent(config, "pcg_large_threshold", 500);
        AccuracyProfile::setIfAbsent(config, "pcg_large_system_iterations", 100);
        AccuracyProfile::setIfAbsent(config, "convergence_threshold", 1e-6);
        AccuracyProfile::setIfAbsent(config, "max_iterations", 50);
        AccuracyProfile::setIfAbsent(config, "eeq_distance_cutoff", 30.0);
        AccuracyProfile::setIfAbsent(config, "allow_unconverged_charges", false);
        break;
    case AccuracyLevel::Medium:
        AccuracyProfile::setIfAbsent(config, "max_pcg_iterations", 500);
        AccuracyProfile::setIfAbsent(config, "pcg_tolerance", 1e-12);
        AccuracyProfile::setIfAbsent(config, "pcg_large_system_tol_factor", 1e-8);
        AccuracyProfile::setIfAbsent(config, "pcg_large_threshold", 1000);
        AccuracyProfile::setIfAbsent(config, "pcg_large_system_iterations", 10000);
        AccuracyProfile::setIfAbsent(config, "convergence_threshold", 1e-8);
        AccuracyProfile::setIfAbsent(config, "max_iterations", 100);
        AccuracyProfile::setIfAbsent(config, "eeq_distance_cutoff", 40.0);
        AccuracyProfile::setIfAbsent(config, "allow_unconverged_charges", false);
        break;
    case AccuracyLevel::High:
        AccuracyProfile::setIfAbsent(config, "max_pcg_iterations", 1000);
        AccuracyProfile::setIfAbsent(config, "pcg_tolerance", 1e-14);
        AccuracyProfile::setIfAbsent(config, "pcg_large_system_tol_factor", 1e-10);
        AccuracyProfile::setIfAbsent(config, "pcg_large_threshold", 2000);
        AccuracyProfile::setIfAbsent(config, "pcg_large_system_iterations", 20000);
        AccuracyProfile::setIfAbsent(config, "convergence_threshold", 1e-10);
        AccuracyProfile::setIfAbsent(config, "max_iterations", 200);
        AccuracyProfile::setIfAbsent(config, "eeq_distance_cutoff", 50.0);
        AccuracyProfile::setIfAbsent(config, "allow_unconverged_charges", false);
        break;
    default:
        break;
    }
}

static void apply_opt(json& config, AccuracyLevel level)
{
    switch (level) {
    case AccuracyLevel::Loose:
        AccuracyProfile::setIfAbsent(config, "energy_threshold", 1.0);
        AccuracyProfile::setIfAbsent(config, "gradient_threshold", 1e-3);
        AccuracyProfile::setIfAbsent(config, "rmsd_threshold", 0.1);
        AccuracyProfile::setIfAbsent(config, "max_iterations", 2000);
        AccuracyProfile::setIfAbsent(config, "lbfgs_eps_abs", 1e-4);
        AccuracyProfile::setIfAbsent(config, "lbfgs_ftol", 1e-3);
        AccuracyProfile::setIfAbsent(config, "opt_level", -1);
        AccuracyProfile::setIfAbsent(config, "maxmicro", 10);
        break;
    case AccuracyLevel::Normal:
        AccuracyProfile::setIfAbsent(config, "energy_threshold", 0.1);
        AccuracyProfile::setIfAbsent(config, "gradient_threshold", 5e-4);
        AccuracyProfile::setIfAbsent(config, "rmsd_threshold", 0.01);
        AccuracyProfile::setIfAbsent(config, "max_iterations", 5000);
        AccuracyProfile::setIfAbsent(config, "lbfgs_eps_abs", 1e-5);
        AccuracyProfile::setIfAbsent(config, "lbfgs_ftol", 1e-4);
        AccuracyProfile::setIfAbsent(config, "opt_level", 0);
        AccuracyProfile::setIfAbsent(config, "maxmicro", 20);
        break;
    case AccuracyLevel::Medium:
        AccuracyProfile::setIfAbsent(config, "energy_threshold", 0.01);
        AccuracyProfile::setIfAbsent(config, "gradient_threshold", 1e-4);
        AccuracyProfile::setIfAbsent(config, "rmsd_threshold", 0.001);
        AccuracyProfile::setIfAbsent(config, "max_iterations", 10000);
        AccuracyProfile::setIfAbsent(config, "lbfgs_eps_abs", 1e-6);
        AccuracyProfile::setIfAbsent(config, "lbfgs_ftol", 1e-5);
        AccuracyProfile::setIfAbsent(config, "opt_level", 1);
        AccuracyProfile::setIfAbsent(config, "maxmicro", 40);
        break;
    case AccuracyLevel::High:
        AccuracyProfile::setIfAbsent(config, "energy_threshold", 0.001);
        AccuracyProfile::setIfAbsent(config, "gradient_threshold", 1e-5);
        AccuracyProfile::setIfAbsent(config, "rmsd_threshold", 0.0001);
        AccuracyProfile::setIfAbsent(config, "max_iterations", 20000);
        AccuracyProfile::setIfAbsent(config, "lbfgs_eps_abs", 1e-7);
        AccuracyProfile::setIfAbsent(config, "lbfgs_ftol", 1e-6);
        AccuracyProfile::setIfAbsent(config, "opt_level", 2);
        AccuracyProfile::setIfAbsent(config, "maxmicro", 80);
        break;
    default:
        break;
    }
}

static void apply_simplemd(json& config, AccuracyLevel level)
{
    switch (level) {
    case AccuracyLevel::Loose:
        AccuracyProfile::setIfAbsent(config, "rattle_tol_12", 1e-3);
        AccuracyProfile::setIfAbsent(config, "rattle_tol_13", 1e-2);
        AccuracyProfile::setIfAbsent(config, "rattle_max_iterations", 50);
        AccuracyProfile::setIfAbsent(config, "time_step", 2.0);
        break;
    case AccuracyLevel::Normal:
        AccuracyProfile::setIfAbsent(config, "rattle_tol_12", 1e-4);
        AccuracyProfile::setIfAbsent(config, "rattle_tol_13", 1e-3);
        AccuracyProfile::setIfAbsent(config, "rattle_max_iterations", 100);
        AccuracyProfile::setIfAbsent(config, "time_step", 1.0);
        break;
    case AccuracyLevel::Medium:
        AccuracyProfile::setIfAbsent(config, "rattle_tol_12", 1e-5);
        AccuracyProfile::setIfAbsent(config, "rattle_tol_13", 1e-4);
        AccuracyProfile::setIfAbsent(config, "rattle_max_iterations", 200);
        AccuracyProfile::setIfAbsent(config, "time_step", 0.5);
        break;
    case AccuracyLevel::High:
        AccuracyProfile::setIfAbsent(config, "rattle_tol_12", 1e-6);
        AccuracyProfile::setIfAbsent(config, "rattle_tol_13", 1e-5);
        AccuracyProfile::setIfAbsent(config, "rattle_max_iterations", 500);
        AccuracyProfile::setIfAbsent(config, "time_step", 0.2);
        break;
    default:
        break;
    }
}

static void apply_ancopt(json& config, AccuracyLevel level)
{
    switch (level) {
    case AccuracyLevel::Loose:
        AccuracyProfile::setIfAbsent(config, "energy_threshold", 5e-4);
        AccuracyProfile::setIfAbsent(config, "gradient_threshold", 1e-2);
        AccuracyProfile::setIfAbsent(config, "opt_level", -1);
        AccuracyProfile::setIfAbsent(config, "maxmicro", 10);
        break;
    case AccuracyLevel::Normal:
        AccuracyProfile::setIfAbsent(config, "energy_threshold", 5e-6);
        AccuracyProfile::setIfAbsent(config, "gradient_threshold", 1.890e-3);
        AccuracyProfile::setIfAbsent(config, "opt_level", 0);
        AccuracyProfile::setIfAbsent(config, "maxmicro", 20);
        break;
    case AccuracyLevel::Medium:
        AccuracyProfile::setIfAbsent(config, "energy_threshold", 1e-6);
        AccuracyProfile::setIfAbsent(config, "gradient_threshold", 8e-4);
        AccuracyProfile::setIfAbsent(config, "opt_level", 1);
        AccuracyProfile::setIfAbsent(config, "maxmicro", 40);
        break;
    case AccuracyLevel::High:
        AccuracyProfile::setIfAbsent(config, "energy_threshold", 1e-7);
        AccuracyProfile::setIfAbsent(config, "gradient_threshold", 2e-4);
        AccuracyProfile::setIfAbsent(config, "opt_level", 2);
        AccuracyProfile::setIfAbsent(config, "maxmicro", 80);
        break;
    default:
        break;
    }
}

static void apply_modern_optimizer(json& config, AccuracyLevel level)
{
    switch (level) {
    case AccuracyLevel::Loose:
        AccuracyProfile::setIfAbsent(config, "d_e", 1e-4);
        AccuracyProfile::setIfAbsent(config, "grad_norm", 1e-3);
        AccuracyProfile::setIfAbsent(config, "trust_radius", 0.1);
        AccuracyProfile::setIfAbsent(config, "max_iterations", 2000);
        break;
    case AccuracyLevel::Normal:
        AccuracyProfile::setIfAbsent(config, "d_e", 1e-6);
        AccuracyProfile::setIfAbsent(config, "grad_norm", 5e-4);
        AccuracyProfile::setIfAbsent(config, "trust_radius", 0.05);
        AccuracyProfile::setIfAbsent(config, "max_iter", 5000);
        break;
    case AccuracyLevel::Medium:
        AccuracyProfile::setIfAbsent(config, "d_e", 1e-7);
        AccuracyProfile::setIfAbsent(config, "grad_norm", 1e-4);
        AccuracyProfile::setIfAbsent(config, "trust_radius", 0.02);
        AccuracyProfile::setIfAbsent(config, "max_iter", 10000);
        break;
    case AccuracyLevel::High:
        AccuracyProfile::setIfAbsent(config, "d_e", 1e-8);
        AccuracyProfile::setIfAbsent(config, "grad_norm", 1e-5);
        AccuracyProfile::setIfAbsent(config, "trust_radius", 0.01);
        AccuracyProfile::setIfAbsent(config, "max_iter", 20000);
        break;
    default:
        break;
    }
}

// ============================================================================
// Apply module mapping
// ============================================================================
void AccuracyProfile::applyModule(json& config, const std::string& module, AccuracyLevel level)
{
    std::string m = toLower(module);
    if (m == "gfnff") {
        apply_gfnff(config, level);
    } else if (m == "eeq_solver") {
        apply_eeq_solver(config, level);
    } else if (m == "opt") {
        apply_opt(config, level);
    } else if (m == "simplemd" || m == "md") {
        apply_simplemd(config, level);
    } else if (m == "ancopt") {
        apply_ancopt(config, level);
    } else if (m == "modern_optimizer") {
        apply_modern_optimizer(config, level);
    }
}

// ============================================================================
// Public API
// ============================================================================
void AccuracyProfile::apply(json& config, const std::string& module, const std::string& level)
{
    AccuracyLevel l = parseLevel(level);
    if (l == AccuracyLevel::Unknown) {
        return; // Invalid level — leave defaults untouched
    }
    applyModule(config, module, l);
}

void AccuracyProfile::applyToAll(json& controller)
{
    if (!controller.contains("accuracy")) {
        return;
    }

    std::string level = controller.value("accuracy", "normal");
    AccuracyLevel l = parseLevel(level);
    if (l == AccuracyLevel::Unknown) {
        return;
    }

    // Apply to all known sub-controllers
    const std::vector<std::string> modules = {
        "opt", "simplemd", "gfnff", "eeq_solver", "ancopt", "modern_optimizer"
    };

    for (const auto& mod : modules) {
        if (controller.contains(mod) && controller[mod].is_object()) {
            applyModule(controller[mod], mod, l);
        } else {
            // Ensure sub-controller exists so parameters are available
            controller[mod] = controller.value(mod, json::object());
            applyModule(controller[mod], mod, l);
        }
    }

    // Claude Generated (April 2026): Cross-module parameter forwarding.
    // Parameters like skip_phase2 and allow_unconverged_charges may be specified
    // at the capability level (opt, simplemd) or top-level but must reach gfnff/eeq_solver.
    const std::vector<std::string> forward_keys = {"skip_phase2", "allow_unconverged_charges"};
    const std::vector<std::string> dst_modules = {"gfnff", "eeq_solver"};
    for (const auto& key : forward_keys) {
        // Collect the value from the first source that has it (top-level or any module)
        json value;
        bool has_value = false;
        // Check top-level first
        if (controller.contains(key) && !controller[key].is_object()) {
            value = controller[key];
            has_value = true;
        }
        // Then check modules
        for (const auto& src_mod : modules) {
            if (!controller.contains(src_mod) || !controller[src_mod].is_object())
                continue;
            if (controller[src_mod].contains(key)) {
                value = controller[src_mod][key];
                has_value = true;
                break; // First match wins
            }
        }
        if (!has_value)
            continue;
        // Forward to destination modules
        for (const auto& dst_mod : dst_modules) {
            if (!controller.contains(dst_mod) || !controller[dst_mod].is_object())
                continue;
            if (!controller[dst_mod].contains(key)) {
                controller[dst_mod][key] = value;
            }
        }
    }
}

json AccuracyProfile::getMappedValue(const std::string& module, const std::string& param, AccuracyLevel level)
{
    json config = json::object();
    applyModule(config, module, level);
    if (config.contains(param)) {
        return config[param];
    }
    return json(nullptr);
}
