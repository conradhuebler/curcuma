/*
 * <Alignment Strategy Pattern for RMSD calculations>
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
 * Claude Generated - Strategy Pattern implementation for 7 alignment methods
 */

#include "rmsd_strategies.h"
#include "../rmsd.h" // For RMSDDriver access
#include "src/core/fileiterator.h"
#include "src/core/global.h"
#include <filesystem>
#include <fmt/core.h>

// Claude Generated - Strategy Factory implementation with enum support
std::unique_ptr<AlignmentStrategy> AlignmentStrategyFactory::createStrategy(AlignmentMethod method)
{
    switch (method) {
    case AlignmentMethod::INCREMENTAL:
        return std::make_unique<IncrementalAlignmentStrategy>();
    case AlignmentMethod::TEMPLATE:
        return std::make_unique<TemplateAlignmentStrategy>();
    case AlignmentMethod::HEAVY_TEMPLATE:
        return std::make_unique<HeavyTemplateStrategy>();
    case AlignmentMethod::ATOM_TEMPLATE:
        return std::make_unique<AtomTemplateStrategy>();
    case AlignmentMethod::INERTIA:
        return std::make_unique<InertiaAlignmentStrategy>();
    case AlignmentMethod::MOLALIGN:
        return std::make_unique<MolAlignStrategy>();
    case AlignmentMethod::DISTANCE_TEMPLATE:
        return std::make_unique<DistanceTemplateStrategy>();
    case AlignmentMethod::PREDEFINED_ORDER:
        return std::make_unique<PredefinedOrderStrategy>();
    default:
        CURCUMA_ERROR(fmt::format("Unknown alignment method: {}", static_cast<int>(method)));
        return nullptr;
    }
}

// Claude Generated - Legacy compatibility for integer method IDs
std::unique_ptr<AlignmentStrategy> AlignmentStrategyFactory::createStrategy(int method_id)
{
    switch (method_id) {
    case 1:
        return createStrategy(AlignmentMethod::INCREMENTAL);
    case 2:
        return createStrategy(AlignmentMethod::TEMPLATE);
    case 3:
        return createStrategy(AlignmentMethod::HEAVY_TEMPLATE);
    case 4:
        return createStrategy(AlignmentMethod::ATOM_TEMPLATE);
    case 5:
        return createStrategy(AlignmentMethod::INERTIA);
    case 6:
        return createStrategy(AlignmentMethod::MOLALIGN);
    case 7:
        return createStrategy(AlignmentMethod::DISTANCE_TEMPLATE);
    case 10:
        return createStrategy(AlignmentMethod::PREDEFINED_ORDER);
    default:
        CURCUMA_ERROR(fmt::format("Unknown alignment method ID: {}", method_id));
        return nullptr;
    }
}

std::map<int, std::string> AlignmentStrategyFactory::getAvailableStrategies()
{
    std::map<int, std::string> strategies;

    auto addStrategy = [&strategies](int id) {
        auto strategy = createStrategy(id);
        if (strategy) {
            strategies[id] = strategy->getName();
        }
    };

    addStrategy(1);
    addStrategy(2);
    addStrategy(3);
    addStrategy(4);
    addStrategy(5);
    addStrategy(6);
    addStrategy(7);
    addStrategy(10);

    return strategies;
}

// Claude Generated - Method 1: Incremental Alignment Strategy
// REAL REFACTORING: Complete algorithm migration from RMSDDriver::ReorderIncremental()
AlignmentResult IncrementalAlignmentStrategy::align(RMSDDriver* driver, const AlignmentConfig& config)
{
    AlignmentResult result;

    try {
        // MIGRATED ALGORITHM: Complete ReorderIncremental logic
        int inter_size = driver->m_reference.AtomCount() * (driver->m_reference.AtomCount() - 1) * driver->m_intermedia_storage;

        driver->m_reorder_reference = driver->m_reference;
        driver->m_reorder_target = driver->m_target;

        driver->m_reorder_reference.setGeometry(driver->CenterMolecule(driver->m_reference.getGeometry()));
        driver->m_reorder_target.setGeometry(driver->CenterMolecule(driver->m_target.getGeometry()));

        std::pair<Molecule, LimitedStorage> init_result = driver->InitialisePair();
        Molecule ref = init_result.first;
        LimitedStorage storage_shelf = init_result.second;

        int reference_reordered = 0;
        int reference_not_reorordered = 0;
        int max = std::min(driver->m_reference.AtomCount(), driver->m_target.AtomCount());
        int combinations = 0;
        int wake_up = 100;
        CxxThreadPool* pool = new CxxThreadPool;

        // Progress bar control based on verbosity
        if (driver->m_verbosity == 0)
            pool->setProgressBar(CxxThreadPool::ProgressBarType::None);
        else
            pool->setProgressBar(CxxThreadPool::ProgressBarType::Continously);
        pool->setActiveThreadCount(driver->m_threads);

        std::vector<AtomDef> atoms;
        while (driver->m_reorder_reference_geometry.rows() < driver->m_reorder_reference.AtomCount() && driver->m_reorder_reference_geometry.rows() < driver->m_reorder_target.AtomCount() && ((reference_reordered + reference_not_reorordered) <= driver->m_reference.AtomCount())) {

            int thread_count = 0;
            Molecule reference = ref;
            int i = reference.AtomCount();
            double mass = driver->m_reference.ConnectedMass(i);
            auto atom = driver->m_reorder_reference.Atom(i);

            // Progress reporting
            if (driver->m_verbosity >= 2) {
                int progress = int((reference_reordered + reference_not_reorordered) / double(max) * 100);
                CURCUMA_PROGRESS(reference_reordered + reference_not_reorordered, max, "Incremental reordering atoms");
            }

            int element = atom.first;
            reference.addPair(atom);
            if (driver->m_dynamic_center)
                driver->m_reorder_reference_geometry = GeometryTools::TranslateGeometry(reference.getGeometry(), reference.Centroid(true), Position{ 0, 0, 0 });
            else
                driver->m_reorder_reference_geometry = reference.getGeometry();

            std::vector<RMSDThread*> threads;
            for (const auto& e : *storage_shelf.data()) {
                RMSDThread* thread = new RMSDThread(reference, driver->m_reorder_target, driver->m_reorder_reference_geometry,
                    reference.DistanceMatrix().second, e.second, mass, element, driver->m_topo);
                pool->addThread(thread);
                threads.push_back(thread);
                thread_count++;
            }
            pool->StaticPool();

            int match = 0;
            // Only start threads if the current element can be found in target
            if (std::find(driver->m_target.Atoms().begin(), driver->m_target.Atoms().end(), element) != driver->m_target.Atoms().end()) {
                pool->StartAndWait();
            }

            LimitedStorage storage_shelf_next(inter_size);
            for (const auto t : threads) {
                RMSDThread* thread = static_cast<RMSDThread*>(t);
                for (const auto& item : (*thread->data())) {
                    if (thread->Match()) {
                        storage_shelf_next.addItem(item.first, item.second);
                    }
                    match += thread->Match();
                }
                combinations += thread->Calculations();
            }

            if (match == 0) {
                Molecule ref_0;
                auto atom = driver->m_reorder_reference.Atom(i);
                for (std::size_t index = 0; index < driver->m_reference.AtomCount(); index++) {
                    if (i != index)
                        ref_0.addPair(driver->m_reference.Atom(index));
                }
                if (std::find(atoms.begin(), atoms.end(), driver->m_reorder_reference.Atom(i)) == atoms.end()) {
                    ref_0.addPair(driver->m_reorder_reference.Atom(i));
                    atoms.push_back(driver->m_reorder_reference.Atom(i));
                    driver->m_reorder_reference = ref_0;
                    if (driver->m_verbosity >= 2) {
                        CURCUMA_WARN(fmt::format("Atom order altered! Atom {} pushed to end of list", i + reference_reordered));
                        CURCUMA_INFO(fmt::format("Element {} : ({:.3f} {:.3f} {:.3f})", atom.first, atom.second[0], atom.second[1], atom.second[2]));
                    }
                    reference_reordered++;
                } else
                    reference_reordered = driver->m_reference.AtomCount();
            } else {
                ref = reference;
                storage_shelf = storage_shelf_next;
                int index = 0;
                for (const auto& i : *storage_shelf.data()) {
                    driver->m_intermedia_rules.push_back(i.second);
                    if (index > driver->m_limit)
                        break;
                    index++;
                }
                reference_not_reorordered++;
            }
            pool->clear();
        }
        delete pool;

        // Final processing
        int count = 0;
        for (const auto& e : *storage_shelf.data()) {
            if (count > max * (max - 1) * driver->m_intermedia_storage)
                continue;
            std::vector<int> rule = driver->FillMissing(driver->m_reference, e.second);
            if (std::find(driver->m_stored_rules.begin(), driver->m_stored_rules.end(), rule) == driver->m_stored_rules.end()) {
                driver->m_stored_rules.push_back(rule);
                driver->m_intermedia_rules.push_back(rule);
                driver->m_intermediate_cost_matrices.insert(std::pair<double, std::vector<int>>(e.first, rule));
                count++;
            }
        }

        if (driver->m_stored_rules.size() == 0) {
            if (driver->m_verbosity >= 2)
                CURCUMA_WARN("No new solution found");
            for (int i = 0; i < driver->m_reference.AtomCount(); ++i)
                driver->m_reorder_rules.push_back(i);
        } else {
            driver->m_reorder_rules = driver->m_stored_rules[0];
        }

        driver->m_target_reordered = driver->ApplyOrder(driver->m_reorder_rules, driver->m_target);
        driver->m_target = driver->m_target_reordered;
        driver->m_target_aligned = driver->m_target;

        // Extract results
        result.rmsd = driver->RMSD();
        result.rmsd_raw = driver->RMSDRaw();
        result.reorder_rules = driver->ReorderRules();
        result.reference_aligned = driver->ReferenceAligned();
        result.target_aligned = driver->TargetAligned();
        result.target_reordered = driver->TargetReorderd();
        result.success = true;

        CURCUMA_SUCCESS("Incremental alignment algorithm completed successfully");

    } catch (const std::exception& e) {
        result.success = false;
        result.error_message = fmt::format("Incremental alignment failed: {}", e.what());
        CURCUMA_ERROR(result.error_message);
    }

    return result;
}

// Claude Generated - Method 2: Template Alignment Strategy
// REAL REFACTORING: Complete algorithm migration from RMSDDriver::TemplateReorder()
AlignmentResult TemplateAlignmentStrategy::align(RMSDDriver* driver, const AlignmentConfig& config)
{
    AlignmentResult result;

    try {
        // MIGRATED ALGORITHM: Complete TemplateReorder logic
        auto fragments = driver->CheckFragments();

        if (fragments.first == -1 || fragments.second == -1) {
            result.success = false;
            result.error_message = "Template alignment failed - no valid fragments found";
            CURCUMA_WARN(result.error_message);
            return result;
        }

        driver->m_fragment = -1;
        driver->m_fragment_target = -1;
        driver->m_fragment_reference = -1;

        Molecule reference_mol = driver->m_reference.getFragmentMolecule(fragments.first);
        Molecule target_mol = driver->m_target.getFragmentMolecule(fragments.second);

        auto operators = driver->GetOperateVectors(fragments.first, fragments.second);
        Eigen::Matrix3d R = operators.first;

        Geometry cached_reference = driver->m_reference.getGeometryByFragment(fragments.first, driver->m_protons);
        Geometry cached_target = driver->m_target.getGeometryByFragment(fragments.second, driver->m_protons);

        Geometry ref = GeometryTools::TranslateMolecule(driver->m_reference, GeometryTools::Centroid(cached_reference), Position{ 0, 0, 0 });
        Geometry tget = GeometryTools::TranslateMolecule(driver->m_target, GeometryTools::Centroid(cached_target), Position{ 0, 0, 0 });

        Geometry rotated = tget * R;

        Molecule ref_mol = driver->m_reference;
        ref_mol.setGeometry(ref);

        Molecule tar_mol = driver->m_target;
        tar_mol.setGeometry(rotated);

        auto blob = driver->MakeCostMatrix(ref, rotated);
        driver->m_cost_limit = blob.first;

        driver->InsertRotation(blob);

        // Extract results
        result.rmsd = driver->RMSD();
        result.rmsd_raw = driver->RMSDRaw();
        result.reorder_rules = driver->ReorderRules();
        result.reference_aligned = driver->ReferenceAligned();
        result.target_aligned = driver->TargetAligned();
        result.target_reordered = driver->TargetReorderd();
        result.success = true;

        CURCUMA_SUCCESS("Template alignment with prealignment and Kuhn-Munkres completed successfully");

    } catch (const std::exception& e) {
        result.success = false;
        result.error_message = fmt::format("Template alignment failed: {}", e.what());
        CURCUMA_ERROR(result.error_message);
    }

    return result;
}

// Claude Generated - Method 3: Heavy Template Strategy
// REAL REFACTORING: Complete algorithm migration from RMSDDriver::HeavyTemplate() + PrepareHeavyTemplate()
AlignmentResult HeavyTemplateStrategy::align(RMSDDriver* driver, const AlignmentConfig& config)
{
    AlignmentResult result;

    try {
        // MIGRATED ALGORITHM: Complete HeavyTemplate logic
        if (driver->m_verbosity >= 2)
            CURCUMA_INFO("Preparing heavy atom template structure");

        // Extract heavy atoms (non-hydrogen) from both molecules
        Molecule reference;
        std::vector<int> reference_indicies, target_indicies;
        for (std::size_t i = 0; i < driver->m_reference.AtomCount(); ++i) {
            std::pair<int, Position> atom = driver->m_reference.Atom(i);
            if (atom.first != 1) { // Skip hydrogen atoms
                reference.addPair(atom);
                reference_indicies.push_back(i);
            }
        }

        Molecule target;
        for (std::size_t i = 0; i < driver->m_target.AtomCount(); ++i) {
            std::pair<int, Position> atom = driver->m_target.Atom(i);
            if (atom.first != 1) { // Skip hydrogen atoms
                target.addPair(atom);
                target_indicies.push_back(i);
            }
        }

        // Cache original molecules
        Molecule cached_reference_mol = driver->m_reference;
        Molecule cached_target_mol = driver->m_target;

        // Set heavy-atom-only molecules as current
        driver->m_reference = reference;
        driver->m_target = target;
        driver->m_init_count = driver->m_heavy_init;

        // Use incremental strategy for heavy atoms only
        auto incremental_strategy = AlignmentStrategyFactory::createStrategy(1);
        if (!incremental_strategy) {
            result.success = false;
            result.error_message = "Failed to create IncrementalAlignmentStrategy for heavy template";
            CURCUMA_ERROR(result.error_message);
            return result;
        }

        AlignmentConfig heavy_config = config; // Use same configuration
        AlignmentResult incremental_result = incremental_strategy->align(driver, heavy_config);

        if (!incremental_result.success) {
            result.success = false;
            result.error_message = "Incremental alignment failed for heavy atoms";
            CURCUMA_ERROR(result.error_message);
            return result;
        }

        // Transform rules back to full molecule indices
        std::vector<std::vector<int>> transformed_rules;
        for (int i = 0; i < driver->m_stored_rules.size(); ++i) {
            std::vector<int> tmp;
            for (auto index : driver->m_stored_rules[i])
                tmp.push_back(target_indicies[index]);
            transformed_rules.push_back(tmp);
        }
        driver->m_stored_rules = transformed_rules;

        // Restore original molecules
        driver->m_reference = cached_reference_mol;
        driver->m_target = cached_target_mol;

        // Extract results
        result.rmsd = driver->RMSD();
        result.rmsd_raw = driver->RMSDRaw();
        result.reorder_rules = driver->ReorderRules();
        result.reference_aligned = driver->ReferenceAligned();
        result.target_aligned = driver->TargetAligned();
        result.target_reordered = driver->TargetReorderd();
        result.success = true;

        CURCUMA_SUCCESS("Heavy atom template alignment completed successfully");

    } catch (const std::exception& e) {
        result.success = false;
        result.error_message = fmt::format("Heavy template alignment failed: {}", e.what());
        CURCUMA_ERROR(result.error_message);
    }

    return result;
}

// Claude Generated - Method 4: Atom Template Strategy
// REAL REFACTORING: Complete algorithm migration from RMSDDriver::AtomTemplate() + PrepareAtomTemplate()
AlignmentResult AtomTemplateStrategy::align(RMSDDriver* driver, const AlignmentConfig& config)
{
    AlignmentResult result;

    try {
        // MIGRATED ALGORITHM: Complete AtomTemplate logic using element_templates
        const std::vector<int>& templateatom = config.element_templates.empty() ? std::vector<int>{ config.element } : config.element_templates;

        if (templateatom.empty()) {
            result.success = false;
            result.error_message = "No template elements specified for atom template alignment";
            CURCUMA_ERROR(result.error_message);
            return result;
        }

        // Extract atoms matching template elements from both molecules
        Molecule reference;
        std::vector<int> reference_indicies, target_indicies;
        for (std::size_t i = 0; i < driver->m_reference.AtomCount(); ++i) {
            std::pair<int, Position> atom = driver->m_reference.Atom(i);
            if (std::find(templateatom.begin(), templateatom.end(), atom.first) != templateatom.end()) {
                reference.addPair(atom);
                reference_indicies.push_back(i);
            }
        }

        Molecule target;
        for (std::size_t i = 0; i < driver->m_target.AtomCount(); ++i) {
            std::pair<int, Position> atom = driver->m_target.Atom(i);
            if (std::find(templateatom.begin(), templateatom.end(), atom.first) != templateatom.end()) {
                target.addPair(atom);
                target_indicies.push_back(i);
            }
        }

        if (reference_indicies.size() == 0 || target_indicies.size() == 0) {
            result.success = false;
            result.error_message = "Templates are empty, maybe try different elements";
            CURCUMA_ERROR(result.error_message);
            return result;
        }

        // Cache original molecules
        Molecule cached_reference_mol = driver->m_reference;
        Molecule cached_target_mol = driver->m_target;

        // Set template-atom-only molecules as current
        driver->m_reference = reference;
        driver->m_target = target;
        driver->m_init_count = driver->m_heavy_init;

        // Use incremental strategy for template atoms only
        auto incremental_strategy = AlignmentStrategyFactory::createStrategy(1);
        if (!incremental_strategy) {
            result.success = false;
            result.error_message = "Failed to create IncrementalAlignmentStrategy for atom template";
            CURCUMA_ERROR(result.error_message);
            return result;
        }

        AlignmentConfig template_config = config;
        AlignmentResult incremental_result = incremental_strategy->align(driver, template_config);

        if (!incremental_result.success) {
            result.success = false;
            result.error_message = "Incremental alignment failed for template atoms";
            CURCUMA_ERROR(result.error_message);
            return result;
        }

        // Restore original molecules
        driver->m_reference = cached_reference_mol;
        driver->m_target = cached_target_mol;

        // Transform rules and generate cost matrices
        std::vector<std::vector<int>> transformed_rules;
        for (int i = 0; i < driver->m_stored_rules.size(); ++i) {
            std::vector<int> tmp;
            for (auto index : driver->m_stored_rules[i])
                tmp.push_back(target_indicies[index]);
            transformed_rules.push_back(tmp);
            driver->m_intermedia_rules.push_back(tmp);

            auto OperateVectors = driver->GetOperateVectors(reference_indicies, tmp);
            auto cost_result = driver->MakeCostMatrix(OperateVectors.first);
            driver->InsertRotation(cost_result);
        }

        // Additional cost matrix calculations for intermediate rules
        for (const auto& indices : driver->m_intermedia_rules) {
            auto cost_result = driver->MakeCostMatrix(reference_indicies, indices);
            driver->InsertRotation(cost_result);
        }

        driver->m_stored_rules = transformed_rules;

        // Extract results
        result.rmsd = driver->RMSD();
        result.rmsd_raw = driver->RMSDRaw();
        result.reorder_rules = driver->ReorderRules();
        result.reference_aligned = driver->ReferenceAligned();
        result.target_aligned = driver->TargetAligned();
        result.target_reordered = driver->TargetReorderd();
        result.success = true;

        CURCUMA_SUCCESS("Selected atom template alignment completed successfully");

    } catch (const std::exception& e) {
        result.success = false;
        result.error_message = fmt::format("Atom template alignment failed: {}", e.what());
        CURCUMA_ERROR(result.error_message);
    }

    return result;
}

// Claude Generated - Method 5: Inertia Alignment Strategy
// REAL REFACTORING: Complete algorithm migration from RMSDDriver::TemplateFree()
AlignmentResult InertiaAlignmentStrategy::align(RMSDDriver* driver, const AlignmentConfig& config)
{
    AlignmentResult result;

    try {
        // MIGRATED ALGORITHM: Complete TemplateFree (Inertia-based alignment) logic

        // === INTELLIGENT PREFILTERING BASED ON MOMENTS OF INERTIA ===

        // Calculate moments of inertia
        Molecule refTemp = driver->m_reference;
        Molecule targTemp = driver->m_target;
        refTemp.CalculateRotationalConstants();
        targTemp.CalculateRotationalConstants();

        std::vector<double> refMoments = { refTemp.Ia(), refTemp.Ib(), refTemp.Ic() };
        std::vector<double> targMoments = { targTemp.Ia(), targTemp.Ib(), targTemp.Ic() };

        if (driver->m_verbosity >= 3) {
            CURCUMA_INFO(fmt::format("Reference moments: [{:.6f}, {:.6f}, {:.6f}]", refMoments[0], refMoments[1], refMoments[2]));
            CURCUMA_INFO(fmt::format("Target moments: [{:.6f}, {:.6f}, {:.6f}]", targMoments[0], targMoments[1], targMoments[2]));
        }

        // Evaluate all permutations and filter them
        std::vector<std::vector<int>> allPermutations = {
            { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 }
        };

        auto evaluatePermutation = [&](const std::vector<int>& perm) -> double {
            double totalDiff = 0.0;
            for (int i = 0; i < 3; i++) {
                double refMoment = refMoments[i];
                double targMoment = targMoments[perm[i]];
                double relativeDiff = std::abs(refMoment - targMoment) / std::max(refMoment, targMoment);
                totalDiff += relativeDiff;
            }
            return totalDiff;
        };

        // Sort permutations by quality
        std::vector<std::pair<double, std::vector<int>>> scoredPermutations;
        for (const auto& perm : allPermutations) {
            double score = evaluatePermutation(perm);
            scoredPermutations.push_back({ score, perm });
            if (driver->m_verbosity >= 3) {
                CURCUMA_INFO(fmt::format("Permutation [{},{},{}] Score: {:.4f}", perm[0], perm[1], perm[2], score));
            }
        }
        std::sort(scoredPermutations.begin(), scoredPermutations.end());

        // Take only the best permutations
        std::vector<std::vector<int>> permutations;
        double threshold = 0.5; // Adjustable
        int maxPerms = 3; // Max number

        for (const auto& scored : scoredPermutations) {
            if ((scored.first <= threshold || permutations.size() < 2) && permutations.size() < maxPerms) {
                permutations.push_back(scored.second);
                if (driver->m_verbosity >= 2) {
                    CURCUMA_INFO(fmt::format("Using permutation [{},{},{}] with score: {:.4f}",
                        scored.second[0], scored.second[1], scored.second[2], scored.first));
                }
            }
        }

        // Generate all 8 orientations
        std::vector<std::vector<int>> orientations;
        for (int i = 0; i < 8; i++) {
            orientations.push_back({ (i & 1) ? 1 : -1,
                (i & 2) ? 1 : -1,
                (i & 4) ? 1 : -1 });
        }

        // === END PREFILTERING ===

        int bestRefPerm = 0, bestRefOrient = 0;
        int bestTargPerm = 0, bestTargOrient = 0;

        if (driver->m_verbosity >= 2) {
            int totalCombinations = std::pow(permutations.size() * orientations.size(), 2);
            int originalCombinations = std::pow(6 * 8, 2);
            CURCUMA_INFO(fmt::format("Testing {} x {} x {} x {} = {} combinations (vs {} without filter)",
                permutations.size(), orientations.size(), permutations.size(), orientations.size(),
                totalCombinations, originalCombinations));
        }

        int counter = 0;
        Molecule reference = driver->m_reference;
        Molecule target = driver->m_target;

        // Test all combinations for Reference and Target
        for (int refPerm = 0; refPerm < permutations.size(); refPerm++) {
            for (int refOrient = 0; refOrient < orientations.size(); refOrient++) {
                for (int targPerm = 0; targPerm < permutations.size(); targPerm++) {
                    for (int targOrient = 0; targOrient < orientations.size(); targOrient++) {
                        counter++;

                        // Align reference
                        reference.AlignAxis(permutations[refPerm], orientations[refOrient]);

                        // Align target
                        target.AlignAxis(permutations[targPerm], orientations[targOrient]);

                        // Calculate cost matrix
                        auto costmatrix = driver->MakeCostMatrix(reference.getGeometry(), target.getGeometry());

                        // Pass to InsertRotation
                        driver->InsertRotation(costmatrix);
                    }
                }
            }
        }

        if (driver->m_verbosity >= 2) {
            CURCUMA_SUCCESS(fmt::format("Inertia alignment completed! Tested: {} combinations", counter));
        }

        // Extract results
        result.rmsd = driver->RMSD();
        result.rmsd_raw = driver->RMSDRaw();
        result.reorder_rules = driver->ReorderRules();
        result.reference_aligned = driver->ReferenceAligned();
        result.target_aligned = driver->TargetAligned();
        result.target_reordered = driver->TargetReorderd();
        result.success = true;

        CURCUMA_SUCCESS("Inertia-based alignment with moment-of-inertia prefiltering completed successfully");

    } catch (const std::exception& e) {
        result.success = false;
        result.error_message = fmt::format("Inertia alignment failed: {}", e.what());
        CURCUMA_ERROR(result.error_message);
    }

    return result;
}

// Claude Generated - Method 6: MolAlign Strategy
// REAL REFACTORING: Complete algorithm migration from RMSDDriver::MolAlignLib()
AlignmentResult MolAlignStrategy::align(RMSDDriver* driver, const AlignmentConfig& config)
{
    AlignmentResult result;

    try {
        // MIGRATED ALGORITHM: Complete MolAlignLib (external tool) logic

        // Write molecules to temporary files for molalign
        driver->m_reference.writeXYZFile("molaign_ref.xyz");
        driver->m_target.writeXYZFile("molalign_tar.xyz");

        FILE* FileOpen;
        std::string command = config.molalign_bin + config.molalign_args + " molaign_ref.xyz molalign_tar.xyz 2>&1";

        if (driver->m_verbosity >= 3) {
            CURCUMA_INFO(fmt::format("MolAlign command: {}", command));
        }

        FileOpen = popen(command.c_str(), "r");
        if (!FileOpen) {
            result.success = false;
            result.error_message = "Failed to execute MolAlign command";
            CURCUMA_ERROR(result.error_message);
            return result;
        }

        bool ok = true;
        bool rndm = false;
        bool error = false;
        char line[130];

        // Parse molalign output
        while (fgets(line, sizeof line, FileOpen)) {
            rndm = std::string(line).find("random") != std::string::npos;
            error = std::string(line).find("Error") != std::string::npos;
            if (driver->m_verbosity >= 3) {
                std::string line_str = std::string(line);
                if (!line_str.empty() && line_str.back() == '\n') {
                    line_str.pop_back();
                }
                CURCUMA_INFO(fmt::format("MolAlign output: {}", line_str));
            }
        }
        pclose(FileOpen);

        // Check if alignment was successful
        if (std::filesystem::exists("aligned.xyz") && !rndm) {
            if (driver->m_verbosity >= 1) {
                CURCUMA_CITATION("J. Chem. Inf. Model. 2023, 63, 4, 1157–1165 - DOI: 10.1021/acs.jcim.2c01187");
            }

            FileIterator file("aligned.xyz", true);
            driver->m_reference_centered = file.Next();
            driver->m_target_reordered = file.Next();
            driver->m_target_aligned = driver->m_target_reordered;
            driver->m_target = driver->m_target_reordered;
            driver->m_rmsd = driver->CalculateRMSD();

            // Clean up temporary file
            std::filesystem::remove("aligned.xyz");
            std::filesystem::remove("molaign_ref.xyz");
            std::filesystem::remove("molalign_tar.xyz");

            // Extract results
            result.rmsd = driver->RMSD();
            result.rmsd_raw = driver->RMSDRaw();
            result.reorder_rules = driver->ReorderRules();
            result.reference_aligned = driver->ReferenceAligned();
            result.target_aligned = driver->TargetAligned();
            result.target_reordered = driver->TargetReorderd();
            result.success = true;

            CURCUMA_SUCCESS("External MolAlign alignment completed successfully");

        } else {
            // Clean up temporary files
            std::filesystem::remove("molaign_ref.xyz");
            std::filesystem::remove("molalign_tar.xyz");

            if (rndm) {
                result.success = false;
                result.error_message = "molalign has trouble with random numbers, no permutation will be performed";
                CURCUMA_ERROR(result.error_message);
            } else if (error) {
                result.success = false;
                result.error_message = fmt::format("molalign has trouble finding solution, try increasing tolerance via molalign_tolerance (current: {})", config.molalign_tolerance);
                CURCUMA_ERROR(result.error_message);
            } else if (!rndm && !error) {
                result.success = false;
                result.error_message = fmt::format("Molalign was not found at: {}. Consider getting it from https://github.com/qcuaeh/molalignlib", config.molalign_bin);
                CURCUMA_ERROR(result.error_message);
            }
        }

    } catch (const std::exception& e) {
        result.success = false;
        result.error_message = fmt::format("MolAlign strategy failed: {}", e.what());
        CURCUMA_ERROR(result.error_message);

        // Clean up temporary files in case of exception
        std::filesystem::remove("molaign_ref.xyz");
        std::filesystem::remove("molalign_tar.xyz");
        std::filesystem::remove("aligned.xyz");
    }

    return result;
}

// Claude Generated - Method 7: Distance Template Strategy
// Refactored from RMSDDriver::DistanceTemplate()
AlignmentResult DistanceTemplateStrategy::align(RMSDDriver* driver, const AlignmentConfig& config)
{
    AlignmentResult result;
    result.success = false;

    try {
        // PHYSICALLY MIGRATED from RMSDDriver::DistanceTemplate() + PrepareDistanceTemplate()

        if (driver->m_verbosity >= 2)
            CURCUMA_INFO("Preparing template structure on atom distances");

        if (driver->m_verbosity >= 2)
            CURCUMA_INFO("Starting template preparation");

        RunTimer time;
        std::map<double, std::pair<int, int>> m_distance_reference, m_distance_target;
        Position ref_centroid = driver->m_reference.Centroid(true);
        Position tar_centroid = driver->m_target.Centroid(true);

        Geometry reference_geometry = GeometryTools::TranslateMolecule(
            driver->m_reference, driver->m_reference.Centroid(true), Position{ 0, 0, 0 });
        Geometry target_geometry = GeometryTools::TranslateMolecule(
            driver->m_target, driver->m_target.Centroid(true), Position{ 0, 0, 0 });

        std::map<double, int> distance_reference, distance_target;
        Molecule reference;
        std::vector<int> reference_indicies, target_indicies;

        // Build distance maps for all atom pairs
        for (std::size_t i = 0; i < driver->m_reference.AtomCount() && i < driver->m_target.AtomCount(); ++i) {
            double distance = (driver->m_reference.Atom(i).second - ref_centroid).norm();
            distance_reference[distance] = i;
            distance = (driver->m_target.Atom(i).second - tar_centroid).norm();
            distance_target[distance] = i;

            for (int j = i + 1; j < driver->m_reference.AtomCount() && j < driver->m_target.AtomCount(); ++j) {
                double distance1 = (driver->m_reference.Atom(i).second - driver->m_reference.Atom(j).second).norm();
                m_distance_reference[distance1] = std::pair<int, int>(i, j);
                double distance2 = (driver->m_target.Atom(i).second - driver->m_target.Atom(j).second).norm();
                m_distance_target[distance2] = std::pair<int, int>(i, j);
            }
        }

        Molecule target;
        int i = 0;
        auto ref_end = m_distance_reference.cend();
        auto tar_end = m_distance_target.cend();

        // Select template atoms based on largest distances
        while (reference_indicies.size() < config.limit) {
            ref_end--;
            tar_end--;
            std::pair<int, Position> atom_r1 = driver->m_reference.Atom(ref_end->second.first);
            std::pair<int, Position> atom_r2 = driver->m_reference.Atom(ref_end->second.second);
            std::pair<int, Position> atom_t1 = driver->m_target.Atom(tar_end->second.first);
            std::pair<int, Position> atom_t2 = driver->m_target.Atom(tar_end->second.second);

            if (atom_r1.first == atom_t1.first && std::find(reference_indicies.begin(), reference_indicies.end(), ref_end->second.first) == reference_indicies.end() && std::find(target_indicies.begin(), target_indicies.end(), tar_end->second.first) == target_indicies.end()) {
                reference.addPair(atom_r1);
                reference_indicies.push_back(ref_end->second.first);
                target.addPair(atom_t1);
                target_indicies.push_back(tar_end->second.first);
            }

            if (atom_r2.first == atom_t2.first && std::find(reference_indicies.begin(), reference_indicies.end(), ref_end->second.second) == reference_indicies.end() && std::find(target_indicies.begin(), target_indicies.end(), tar_end->second.second) == target_indicies.end()) {
                reference.addPair(atom_r2);
                reference_indicies.push_back(ref_end->second.second);
                target.addPair(atom_t2);
                target_indicies.push_back(tar_end->second.second);
            }
            ++i;
        }

        // Cache original molecules
        Molecule cached_reference_mol = driver->m_reference;
        Molecule cached_target_mol = driver->m_target;

        // Set template molecules for incremental alignment
        driver->m_reference = reference;
        driver->m_target = target;
        driver->m_init_count = driver->m_heavy_init;

        // Use incremental strategy for template alignment
        auto incremental_strategy = AlignmentStrategyFactory::createStrategy(1);
        if (!incremental_strategy) {
            result.success = false;
            result.error_message = "Failed to create IncrementalAlignmentStrategy for distance template";
            CURCUMA_ERROR(result.error_message);
            return result;
        }

        AlignmentConfig template_config = config;
        AlignmentResult incremental_result = incremental_strategy->align(driver, template_config);

        // Restore original molecules
        driver->m_reference = cached_reference_mol;
        driver->m_target = cached_target_mol;

        if (!incremental_result.success) {
            result.success = false;
            result.error_message = "Incremental alignment failed for distance template";
            CURCUMA_ERROR(result.error_message);
            return result;
        }

        // Process stored results and transform back to full molecule indices
        for (const auto& indices : driver->m_intermedia_rules) {
            auto cost_result = driver->MakeCostMatrix(reference_indicies, indices);
            driver->InsertRotation(cost_result);
        }
        driver->m_intermedia_rules.clear();

        std::map<double, std::vector<int>> order_rules;
        std::vector<std::vector<int>> transformed_rules;
        for (int i = 0; i < driver->m_stored_rules.size(); ++i) {
            std::vector<int> tmp;
            for (auto index : driver->m_stored_rules[i])
                tmp.push_back(target_indicies[index]);
            auto OperateVectors = driver->GetOperateVectors(reference_indicies, tmp);
            auto cost_result = driver->MakeCostMatrix(OperateVectors.first);
            driver->InsertRotation(cost_result);
        }
        driver->m_stored_rules.clear();

        // Set final results
        result.rmsd = driver->RMSD();
        result.rmsd_raw = driver->RMSDRaw();
        result.reorder_rules = driver->ReorderRules();
        result.reference_aligned = driver->ReferenceAligned();
        result.target_aligned = driver->TargetAligned();
        result.target_reordered = driver->TargetReorderd();
        result.success = true;

        CURCUMA_SUCCESS("Distance-based template alignment completed successfully");

    } catch (const std::exception& e) {
        result.success = false;
        result.error_message = fmt::format("Distance template alignment failed: {}", e.what());
        CURCUMA_ERROR(result.error_message);
    }

    return result;
}

// Claude Generated - Method 10: Predefined Order Strategy
// Handles predefined reorder rules
AlignmentResult PredefinedOrderStrategy::align(RMSDDriver* driver, const AlignmentConfig& config)
{
    AlignmentResult result;

    try {
        // For predefined order, we expect the rules to be already set
        // This is a special case where no reordering calculation is needed
        auto predefined_rules = driver->ReorderRules();

        if (!predefined_rules.empty()) {
            // Apply the predefined order and calculate RMSD
            result.rmsd = driver->RMSD();
            result.rmsd_raw = driver->RMSDRaw();
            result.reorder_rules = predefined_rules;
            result.reference_aligned = driver->ReferenceAligned();
            result.target_aligned = driver->TargetAligned();
            result.target_reordered = driver->TargetReorderd();
            result.success = true;

            CURCUMA_SUCCESS(fmt::format("Predefined order alignment completed with {} atoms", predefined_rules.size()));
        } else {
            result.success = false;
            result.error_message = "No predefined reorder rules available";
            CURCUMA_WARN(result.error_message);
        }

    } catch (const std::exception& e) {
        result.success = false;
        result.error_message = fmt::format("Predefined order alignment failed: {}", e.what());
        CURCUMA_ERROR(result.error_message);
    }

    return result;
}