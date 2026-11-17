/*
 * EEQSolver::gpuSupports / gpuMethodNames: which solver names a -gpu backend accepts (name logic only,
 * no device needed). Claude Generated (Oct 2026). AI-generated, machine-tested.
 *   cuda   : every name
 *   rocm   : auto, cholesky, ppcg (Oct 2026: device projected PCG), sytrf (Nov 2026: device dsytrf +
 *            host dsytrs); no pcg/batched/lu/ldlt
 *   vulkan : every name (GFN-FF runs on the CPU there)
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>. GPL-3.0.
 */
#include "src/core/energy_calculators/ff_methods/eeq_solver.h"

#include <cstdio>
#include <string>

static int fails = 0;
static void check(bool ok, const std::string& what)
{
    std::printf("%s  %s\n", ok ? "PASS" : "FAIL", what.c_str());
    if (!ok) ++fails;
}

int main()
{
    const char* all[] = { "auto", "cholesky", "ppcg", "pcg", "lu", "ldlt", "sytrf", "batched" };
    for (const char* n : all) {
        const EEQSolveMethod m = EEQSolver::parseSolveMethod(n);
        check(std::string(EEQSolver::methodName(m)) == n, std::string("name round trip ") + n);
        const bool low = std::string(n) == "auto" || std::string(n) == "cholesky" || std::string(n) == "ppcg"
            || std::string(n) == "sytrf";
        check(EEQSolver::gpuSupports(m, "cuda"), std::string("cuda accepts ") + n);
        check(EEQSolver::gpuSupports(m, "vulkan"), std::string("vulkan accepts ") + n);
        check(EEQSolver::gpuSupports(m, "rocm") == low, std::string("rocm ") + (low ? "accepts " : "rejects ") + n);
    }
    check(EEQSolver::gpuMethodNames("rocm") == "auto, cholesky, ppcg, sytrf", "rocm names: " + EEQSolver::gpuMethodNames("rocm"));
    // Nov 2026: auto never resolves to sytrf, below or above the ppcg thresholds; an explicit sytrf is kept.
    check(EEQSolver::resolveMethod(EEQSolveMethod::Auto, 24, 1, 500, 1) == EEQSolveMethod::SchurCholesky, "auto, 24 atoms -> cholesky");
    check(EEQSolver::resolveMethod(EEQSolveMethod::Auto, 5000, 3, 500, 1) == EEQSolveMethod::ProjectedPCG, "auto, 5000 atoms -> ppcg");
    check(EEQSolver::resolveMethod(EEQSolveMethod::Sytrf, 5000, 3, 500, 1) == EEQSolveMethod::Sytrf, "explicit sytrf kept above the thresholds");
    std::printf("cuda names: %s\nvulkan names: %s\n", EEQSolver::gpuMethodNames("cuda").c_str(),
                EEQSolver::gpuMethodNames("vulkan").c_str());
    return fails ? 1 : 0;
}
