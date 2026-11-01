/*
 * <Convergence presets for geometry optimisation (single source of truth)>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (Aug 2026)
 */

#pragma once

#include <string>

namespace Optimization {

/**
 * @brief Thresholds behind -convergence_preset loose|normal|tight|verytight.
 *
 * One decade per step in every criterion, so the name says what accuracy to expect:
 *
 *  | preset    | dE [kJ/mol] | dRMSD [A] | |g| [Eh/Bohr] | max steps          |
 *  |-----------|-------------|-----------|---------------|--------------------|
 *  | loose     |  1          | 0.1       | 5e-3          | max(50,     1 x N) |
 *  | normal    |  0.1        | 0.01      | 5e-4          | max(500,   10 x N) |
 *  | tight     |  0.01       | 1e-3      | 5e-5          | max(1000,  25 x N) |
 *  | verytight |  0.001      | 1e-4      | 5e-6          | max(5000, 100 x N) |
 *
 * N = number of atoms (see scaledMaxIterations). An explicit -max_iterations switches the size
 * scaling off and is used verbatim.
 *
 * The energy column is the threshold on the energy CHANGE between two steps, and it is also
 * roughly the accuracy the preset delivers: measured on seven MD snapshots of a 107-atom peptide
 * (GFN-FF, identical starting geometries), the old `loose` landed 0.22 kJ/mol (max 0.87) and the
 * old `normal` 0.10 kJ/mol (max 0.28) above the verytight result.
 *
 * The gradient criterion is the NORM of the full 3N vector, so it is system-size dependent:
 * 5e-4 means 7.7e-5 per component for butane but 2.8e-5 for a 107-atom molecule. A large molecule
 * is therefore converged more strictly than a small one at the same preset -- unlike the
 * size-independent max/RMS-component criteria of ORCA or Gaussian. Kept as is for now because
 * every stored convergence behaviour (and every golden test value) rests on it.
 *
 * The step caps are deliberately low: measured on the same snapshots, the thresholds themselves
 * cost almost nothing (656 steps at `loose` vs 750 at `verytight`, i.e. three decades of gradient
 * for 13 % more work), because a cartesian L-BFGS spends its steps getting NEAR the minimum, not
 * on the last digits. What actually makes a crude pre-optimisation cheap is the step cap, so that
 * is where the presets differ by more than a decade.
 *
 * MEASURED across five systems (GFN-FF, five perturbed starting geometries each, deviation from
 * the verytight result; "cap" = how many of the five stopped on the step cap):
 *
 *  | system              | loose        | normal        | tight        | verytight |
 *  |---------------------|--------------|---------------|--------------|-----------|
 *  | butane, 14          | 1.0  (5 cap) | 0.013         | 0.011        | 134 steps |
 *  | caffeine, 24        | 2.1  (4 cap) | 0.93          | 0.001        | 232 steps |
 *  | water cluster, 24   | 252  (3 cap) | 208   (1 cap) | 5.2  (2 cap) | 1056 steps|
 *  | helicene, 42        | 2.3  (5 cap) | 0.003         | 0.000        | 149 steps |
 *  | peptide, 107 (MD)   | 184  (7 cap) | 16.9  (6 cap) | 0.004        | 748 steps |
 *
 * So the caps fit compact molecules up to ~50 atoms and are TOO SMALL for two cases: a weakly
 * bound cluster (the water octamer needs ~1000 steps for its last kJ/mol -- `tight` truncates it)
 * and anything much beyond 100 atoms (`normal` truncates a 107-atom peptide at 17 kJ/mol). Both
 * are reported ("not converged after N"), never silent. A size- or stiffness-aware cap is the
 * obvious next step; the numbers above are the data for choosing it.
 */
struct ConvergencePreset {
    double energy_threshold;   ///< |dE| between steps, kJ/mol
    double rmsd_threshold;     ///< geometry change between steps, Angstrom
    double gradient_threshold; ///< norm of the full 3N gradient, Eh/Bohr
    int max_iterations;        ///< step cap for a small molecule (floor of the size-scaled cap)
    int iterations_per_atom;   ///< size term; effective cap = max(max_iterations, this * atoms)
};

/**
 * @brief Effective step cap for a molecule of @p atoms atoms.
 *
 * A fixed cap cannot fit both ends: the step count needed to reach a minimum grows with the number
 * of degrees of freedom, so 500 steps that fully converge butane truncate a 107-atom peptide at
 * 17 kJ/mol. Measured steps to full convergence (GFN-FF, verytight): butane/14 134, caffeine/24
 * 232, helicene/42 149, peptide/107 748 -- roughly 10 per atom, which is where the `normal` factor
 * comes from; the other presets scale around it by the same decade logic as the thresholds.
 *
 * This fixes the SIZE case, not the STIFFNESS case: a water octamer (24 atoms) needs ~1100 steps
 * because it is floppy and hydrogen-bonded, and no function of the atom count predicts that. Such a
 * run still stops at the cap -- and says so ("not converged after N").
 */
inline int scaledMaxIterations(const ConvergencePreset& preset, int atoms)
{
    if (preset.iterations_per_atom <= 0 || atoms <= 0)
        return preset.max_iterations;
    const long scaled = static_cast<long>(preset.iterations_per_atom) * atoms;
    return static_cast<int>(scaled > preset.max_iterations ? scaled : preset.max_iterations);
}

/**
 * @brief Look up a preset by name.
 * @return false for an unknown name, leaving @p out untouched.
 *
 * Claude Generated (Aug 2026): the four presets used to be written out three times
 * (OptimizationContext::fromJson, buildPresetConfig, CurcumaOpt::LoadControlJson) -- identical by
 * accident, not by construction. This is the single definition all three now read.
 */
inline bool lookupConvergencePreset(const std::string& name, ConvergencePreset& out)
{
    if (name == "loose")
        out = { 1.0, 0.1, 5e-3, 50, 1 };
    else if (name == "normal")
        out = { 0.1, 0.01, 5e-4, 500, 10 };
    else if (name == "tight")
        out = { 0.01, 1e-3, 5e-5, 1000, 25 };
    else if (name == "verytight")
        out = { 0.001, 1e-4, 5e-6, 5000, 100 };
    else
        return false;
    return true;
}

} // namespace Optimization
