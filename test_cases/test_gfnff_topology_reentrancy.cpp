/**
 * @file test_gfnff_topology_reentrancy.cpp
 * @brief Regression guard: calculateTopologyInfoOnce() must be re-entrant.
 *
 * Calling the topology build twice on identical input must produce identical
 * Phase-1 topology charges and identical Phase-2 EEQ charges. It did not:
 * EEQSolver::m_pending_geometry / m_pending_cn are written only by Phase 2
 * (calculateFinalCharges) and were never cleared, and solveWithSchurCholesky
 * uses "are those buffers non-empty?" as its implicit test for "am I Phase 2?".
 * From the second call onward Phase 1 passed that test and consumed (or, after
 * an invalidateCholeskyCache(), poisoned) the Phase-2 Cholesky factor — solving
 * the TOPOLOGICAL-distance system with the factor of the GEOMETRIC matrix.
 *
 * Wrong Phase-1 qa propagates through dgam / alpeeq / gam_corrected into the
 * Phase-2 charges and hence into the GFN-FF Coulomb energy.
 *
 * Two scenarios are covered, because they failed differently:
 *   A) plain repeat call            -> Phase 1 consumed the Phase-2 factor
 *   B) invalidateCholeskyCache() in -> Phase 1 persisted its own factor, which
 *      between the calls (what          Phase 2 of that same call then consumed
 *      gfnff_method.cpp:920 does        (a different wrong value)
 *      on every MD/opt rebuild)
 *
 * Claude Generated (Jul 2026).
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/molecule.h"

using curcuma::Molecule;

namespace {

constexpr double kTol = 1e-12;

/// Max abs deviation between two vectors; -1.0 if the sizes disagree.
double maxDiff(const Eigen::VectorXd& a, const Eigen::VectorXd& b)
{
    if (a.size() != b.size() || a.size() == 0)
        return -1.0;
    return (a - b).cwiseAbs().maxCoeff();
}

void report(const std::string& label, double d, bool& ok)
{
    const bool pass = (d >= 0.0 && d <= kTol);
    std::cout << "  " << std::left << std::setw(44) << label
              << (pass ? "PASS" : "FAIL") << "   max|diff| = "
              << std::scientific << std::setprecision(3) << d << std::endl;
    if (!pass)
        ok = false;
}

/// Build topology twice and compare. `invalidate_between` mirrors what
/// GFNFF::getCachedTopology() does on a full MD/opt topology rebuild.
bool checkMolecule(const std::string& xyz, bool invalidate_between)
{
    Molecule mol(xyz);
    if (mol.AtomCount() == 0) {
        std::cout << "  SKIP (could not load " << xyz << ")" << std::endl;
        return true;
    }

    GFNFF gfnff;
    if (!gfnff.InitialiseMolecule(mol.getMolInfo())) {
        std::cout << "  FAIL (InitialiseMolecule)" << std::endl;
        return false;
    }

    const GFNFF::TopologyInfo first = gfnff.calculateTopologyInfoOnce();
    if (invalidate_between)
        gfnff.invalidateEEQCachesForTest();
    const GFNFF::TopologyInfo second = gfnff.calculateTopologyInfoOnce();

    bool ok = true;
    report("Phase-1 topology charges (qa)",
        maxDiff(first.topology_charges, second.topology_charges), ok);
    report("Phase-2 EEQ charges", maxDiff(first.eeq_charges, second.eeq_charges), ok);
    report("dxi corrections", maxDiff(first.dxi, second.dxi), ok);
    report("dgam corrections", maxDiff(first.dgam, second.dgam), ok);
    return ok;
}

} // namespace

int main(int argc, char** argv)
{
    const std::string root = (argc > 1) ? argv[1] : std::string("test_cases");
    const std::vector<std::string> molecules = {
        root + "/sqm_reference/molecules/H2O.xyz",       // the reported case
        root + "/sqm_reference/molecules/CH4.xyz",       // apolar control
        root + "/sqm_reference/molecules/acetic_acid_dimer.xyz", // nfrag=2, polar
    };

    std::cout << "\n=== GFN-FF topology re-entrancy ===\n" << std::endl;

    bool all_ok = true;
    for (const auto& xyz : molecules) {
        for (int mode = 0; mode < 2; ++mode) {
            const bool invalidate = (mode == 1);
            std::cout << xyz << (invalidate ? "  [cache invalidated between calls]"
                                            : "  [plain repeat call]")
                      << std::endl;
            all_ok &= checkMolecule(xyz, invalidate);
            std::cout << std::endl;
        }
    }

    std::cout << (all_ok ? "ALL PASS" : "FAILURES PRESENT") << "\n" << std::endl;
    return all_ok ? 0 : 1;
}
