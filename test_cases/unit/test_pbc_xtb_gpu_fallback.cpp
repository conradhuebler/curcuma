/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Counted host fallback of the periodic CUDA device stages of native GFN1/GFN2 - Claude Generated (Nov 23, 2026,
 * docs/GPU_PERIODIC_GAPS.md A4/A7).
 *
 * For every periodic device stage token the test injects a failure (env CURCUMA_PBC_GPU_FAIL=<token>: the stage skips
 * its device call and takes its failure path) and checks
 *   (a) the fallback is counted under "xTB periodic device stage <token>" (curcuma::gpuFallbackCounts(), the registry
 *       that -gpu_strict acts on) for the stages the water-dimer run reaches,
 *   (b) energy, gradient and virial equal the host variant of the same stage (CURCUMA_PBC_GPU_HOST=<token>; for the
 *       failed resident step "resstep" the host-driven loop, token resloop): |dE| <= 1e-9 Eh, |dG| <= 1e-8 Eh/A,
 *       |dXi| <= 1e-8 Eh,
 *   (c) a normal run counts no stage fallback,
 *   (d) -gpu_strict ends a child process at the first injected stage failure with exit code 3.
 * Water dimer (registry) in an 8 A cube, scf_threshold 1e-9, gfn1 and gfn2. Without the CUDA plugin or a device the
 * test exits 77 (ctest SKIP_RETURN_CODE), never 0.
 * Argument: the backend, "cuda" (default) or "rocm" (ctest pbc_xtb_gpu_fallback_rocm, ROCm parity A3, Nov 24, 2026: the
 * tokens of the stages that have a HIP port; the GFN1 resident tokens run on the host there; ints since G1B and
 * GFN2 devpot, resloop, resstep, potdl since R25 have a HIP port, Nov 24, 2026).
 */
#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "src/core/gpu_fallback.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <sys/wait.h>
#include <unistd.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <string>

using json = nlohmann::json;

static int g_failures = 0;
static void check(bool ok, const std::string& what)
{
    if (!ok) {
        ++g_failures;
        std::printf("FAIL: %s\n", what.c_str());
    }
}

static const std::string kStagePrefix = "xTB periodic device stage ";

static int stageCount(const std::string& token)
{
    for (const auto& kv : curcuma::gpuFallbackCounts())
        if (kv.first == kStagePrefix + token + ": host loop") return kv.second;
    return 0;
}
static int stageTotal()
{
    int n = 0;
    for (const auto& kv : curcuma::gpuFallbackCounts())
        if (kv.first.rfind(kStagePrefix, 0) == 0) n += kv.second;
    return n;
}
static int otherTotal()
{
    int n = 0;
    for (const auto& kv : curcuma::gpuFallbackCounts())
        if (kv.first.rfind(kStagePrefix, 0) != 0) n += kv.second;
    return n;
}

static Mol dimerInCell(double L)
{
    curcuma::Molecule m = TestMolecules::TestMoleculeRegistry::createMolecule("H2O_dimer", false);
    Mol mol = m.getMolInfo();
    mol.m_unit_cell = Eigen::Matrix3d::Identity() * L;
    mol.m_has_pbc = true;
    return mol;
}

struct Result {
    double e = 0.0;
    Matrix g;
    Eigen::Matrix3d xi = Eigen::Matrix3d::Zero();
    bool ok = false;
};

static std::string g_backend = "cuda";

static Result run(const char* method, const Mol& mol)
{
    json cfg = { { "verbosity", 0 }, { "threads", 1 }, { "gpu", g_backend }, { "xtb", { { "scf_threshold", 1.0e-9 } } } };
    EnergyCalculator calc(method, cfg);
    calc.setMolecule(mol);
    calc.setWantVirial(true);
    Result r;
    r.e = calc.CalculateEnergy(true);
    r.ok = !calc.Error() && calc.virialValid();
    r.g = calc.Gradient();
    r.xi = calc.virial();
    return r;
}

static void setEnv(const char* name, const std::string& v)
{
    if (v.empty()) unsetenv(name);
    else setenv(name, v.c_str(), 1);
}

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    if (argc > 1) g_backend = argv[1];
    const bool rocm = g_backend == "rocm";
    const Mol mol = dimerInCell(8.0);
    int strict_code = -1;
    // (d) -gpu_strict ends the run at the first injected stage failure (child process, exit code 3); forked before
    //     the parent touches the device (a CUDA context does not survive fork)
    {
        std::fflush(stdout);
        const pid_t pid = fork();
        if (pid == 0) {
            setEnv("CURCUMA_PBC_GPU_FAIL", "ewgrad");
            curcuma::setGpuStrict(true);
            run("gfn2", mol);
            std::_Exit(0);
        }
        int status = 0;
        waitpid(pid, &status, 0);
        const int code = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
        std::printf("-gpu_strict with an injected ewgrad failure: child exit code %d (expected %d)\n", code,
                    curcuma::kGpuStrictExitCode);
        strict_code = code;
    }


    // (c) normal run: the GPU path must run (else SKIP) and count no stage fallback
    {
        const Result r = run("gfn2", mol);
        if (otherTotal() > 0) {
            std::printf("test_pbc_xtb_gpu_fallback: SKIP (the -gpu %s request fell back to the CPU: no plugin or device)\n", g_backend.c_str());
            return 77;
        }
        check(r.ok, "gfn2 normal " + g_backend + " run");
    }
    {
        const Result r = run("gfn1", mol);
        check(r.ok, "gfn1 normal " + g_backend + " run");
    }
    check(strict_code == curcuma::kGpuStrictExitCode, "-gpu_strict stops at the first periodic stage fallback");
    std::printf("normal runs: %d stage fallbacks counted\n", stageTotal());
    check(stageTotal() == 0, "a normal periodic " + g_backend + " run counts no stage fallback");

    // Tokens per method; reached = the water-dimer run in an 8 A cube passes through the stage (measured, Nov 23, 2026).
    // extra = tokens kept on the host in both runs (wmat is reached only when h0grad, which builds W itself, is on the host)
    struct Stage { const char* token; const char* host_token; const char* extra = ""; };
    const Stage gfn2_stages[] = { { "ints", "ints" }, { "devpot", "devpot" }, { "resloop", "resloop" },
        { "resstep", "resloop" }, { "potdl", "potdl" }, { "d4dedq", "d4dedq" }, { "gamma", "gamma" },
        { "mpmat", "mpmat" }, { "cn", "cn" }, { "d4cn", "d4cn" }, { "d4pair", "d4pair" }, { "atm", "atm" },
        { "cnfold", "cnfold" }, { "rep", "rep" }, { "mpgrad", "mpgrad" }, { "h0grad", "h0grad" }, { "wmat", "wmat", "h0grad" },
        { "ewgrad", "ewgrad" }, { "kograd", "kograd" }, { "cnchain", "cnchain" } };
    const Stage gfn1_stages[] = { { "ints", "ints" }, { "resloop", "resloop" }, { "resstep", "resloop" },
        { "gamma", "gamma" }, { "cn", "cn" }, { "d3", "d3" }, { "rep", "rep" }, { "h0grad", "h0grad" },
        { "wmat", "wmat", "h0grad" }, { "ewgrad", "ewgrad" }, { "kograd", "kograd" }, { "cnchain", "cnchain" } };
    // ROCm: the stages with a HIP port (rocm/xtb_hip_periodic_host.hiph, G1B, R25)
    const Stage gfn2_rocm[] = { { "ints", "ints" }, { "devpot", "devpot" }, { "resloop", "resloop" },
        { "resstep", "resloop" }, { "potdl", "potdl" }, { "d4dedq", "d4dedq" }, { "gamma", "gamma" },
        { "mpmat", "mpmat" }, { "cn", "cn" },
        { "d4cn", "d4cn" }, { "d4pair", "d4pair" }, { "atm", "atm" }, { "cnfold", "cnfold" }, { "rep", "rep" },
        { "mpgrad", "mpgrad" }, { "h0grad", "h0grad" }, { "wmat", "wmat", "h0grad" }, { "ewgrad", "ewgrad" },
        { "kograd", "kograd" }, { "cnchain", "cnchain" } };
    const Stage gfn1_rocm[] = { { "ints", "ints" }, { "resloop", "resloop" }, { "resstep", "resloop" }, { "gamma", "gamma" }, { "cn", "cn" }, { "d3", "d3" }, { "rep", "rep" },
        { "h0grad", "h0grad" }, { "wmat", "wmat", "h0grad" },
        { "ewgrad", "ewgrad" }, { "kograd", "kograd" }, { "cnchain", "cnchain" } };
    int reached = 0;
    for (const char* method : { "gfn2", "gfn1" }) {
        const bool g2 = std::string(method) == "gfn2";
        const Stage* begin = rocm ? (g2 ? std::begin(gfn2_rocm) : std::begin(gfn1_rocm))
                                  : (g2 ? std::begin(gfn2_stages) : std::begin(gfn1_stages));
        const Stage* end = rocm ? (g2 ? std::end(gfn2_rocm) : std::end(gfn1_rocm))
                                : (g2 ? std::end(gfn2_stages) : std::end(gfn1_stages));
        for (const Stage* s = begin; s != end; ++s) {
            const std::string extra(s->extra);
            setEnv("CURCUMA_PBC_GPU_HOST", extra.empty() ? std::string(s->host_token) : extra + "," + s->host_token);
            const Result host = run(method, mol);
            setEnv("CURCUMA_PBC_GPU_HOST", extra);
            const int c0 = stageCount(s->token);
            setEnv("CURCUMA_PBC_GPU_FAIL", s->token);
            const Result inj = run(method, mol);
            setEnv("CURCUMA_PBC_GPU_FAIL", "");
            setEnv("CURCUMA_PBC_GPU_HOST", "");
            const int counted = stageCount(s->token) - c0;
            const double de = std::abs(inj.e - host.e);
            const double dg = (inj.g.rows() == host.g.rows() && inj.g.size() > 0) ? (inj.g - host.g).cwiseAbs().maxCoeff() : 1.0;
            const double dx = (inj.xi - host.xi).cwiseAbs().maxCoeff();
            std::printf("%s %-8s counted %3d  |dE| %.2e Eh  max|dG| %.2e Eh/A  max|dXi| %.2e Eh (vs host token %s)\n",
                        method, s->token, counted, de, dg, dx, s->host_token);
            const std::string tag = std::string(method) + " stage " + s->token;
            check(host.ok && inj.ok, tag + ": both runs succeed");
            check(counted >= 1, tag + ": the injected failure is counted");
            check(de <= 1.0e-9, tag + ": energy after fallback = host variant (1e-9 Eh)");
            check(dg <= 1.0e-8, tag + ": gradient after fallback = host variant (1e-8 Eh/A)");
            check(dx <= 1.0e-8, tag + ": virial after fallback = host variant (1e-8 Eh)");
            if (counted >= 1) ++reached;
        }
    }

    std::printf("stages with a counted injected fallback: %d\n", reached);
    if (g_failures == 0)
        std::printf("test_pbc_xtb_gpu_fallback: all checks passed\n");
    else
        std::printf("test_pbc_xtb_gpu_fallback: %d failures\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
