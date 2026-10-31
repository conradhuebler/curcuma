/*
 * EEQSolver::gpuSupports / gpuMethodNames: which solver names a -gpu backend accepts (name logic only,
 * no device needed). Claude Generated (Oct 2026). AI-generated, machine-tested.
 *   cuda   : every name
 *   rocm   : auto, cholesky, ppcg (Oct 2026: device projected PCG; no pcg/batched/lu/ldlt)
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
    const char* all[] = { "auto", "cholesky", "ppcg", "pcg", "lu", "ldlt", "batched" };
    for (const char* n : all) {
        const EEQSolveMethod m = EEQSolver::parseSolveMethod(n);
        const bool low = std::string(n) == "auto" || std::string(n) == "cholesky" || std::string(n) == "ppcg";
        check(EEQSolver::gpuSupports(m, "cuda"), std::string("cuda accepts ") + n);
        check(EEQSolver::gpuSupports(m, "vulkan"), std::string("vulkan accepts ") + n);
        check(EEQSolver::gpuSupports(m, "rocm") == low, std::string("rocm ") + (low ? "accepts " : "rejects ") + n);
    }
    check(EEQSolver::gpuMethodNames("rocm") == "auto, cholesky, ppcg", "rocm names: " + EEQSolver::gpuMethodNames("rocm"));
    std::printf("cuda names: %s\nvulkan names: %s\n", EEQSolver::gpuMethodNames("cuda").c_str(),
                EEQSolver::gpuMethodNames("vulkan").c_str());
    return fails ? 1 : 0;
}
