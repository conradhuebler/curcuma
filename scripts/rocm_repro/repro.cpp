// Agent N reproducer: repeated calculations in one process, optional preceding CPU calculation.
// usage: repro <molecule> <method> <nrep> <pre> <threads> [json-extra] [L]
//   pre: none | cpu | gpu (one calculation before the series, not counted)
// prints one line per GPU run: energy %.17g and FNV hash of energy+gradient bytes.
#include "src/core/energycalculator.h"
#include "src/core/curcuma_logger.h"
#include "test_molecule_registry.h"
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cfenv>
#include <xmmintrin.h>
#include <pmmintrin.h>

using json = nlohmann::json;

static uint64_t fnv(const void* p, size_t n, uint64_t h)
{
    const unsigned char* c = static_cast<const unsigned char*>(p);
    for (size_t i = 0; i < n; ++i) { h ^= c[i]; h *= 1099511628211ULL; }
    return h;
}

static void fpstate(const char* tag)
{
    unsigned csr = _mm_getcsr();
    std::fenv_t env; std::fegetenv(&env);
    unsigned short cw; asm volatile("fnstcw %0" : "=m"(cw));
    std::fprintf(stderr, "FPSTATE %s mxcsr=0x%04x (FTZ=%d DAZ=%d RC=%d) x87cw=0x%04x\n", tag, csr & 0xffc0,
                 (csr >> 15) & 1, (csr >> 6) & 1, (csr >> 13) & 3, cw);
}

static double runOnce(const Mol& mol, const char* method, bool gpu, int threads, const json& extra, uint64_t& h)
{
    json cfg = { { "verbosity", 0 }, { "threads", threads }, { "xtb", { { "scf_threshold", 1.0e-9 } } } };
    cfg.merge_patch(extra);
    if (gpu) cfg["gpu"] = "rocm";
    EnergyCalculator calc(method, cfg);
    calc.setMolecule(mol);
    const bool vir = std::getenv("REPRO_VIRIAL") != nullptr;
    if (vir) calc.setWantVirial(true);
    double e = calc.CalculateEnergy(true);
    Matrix g = calc.Gradient();
    h = fnv(&e, sizeof(double), 1469598103934665603ULL);
    h = fnv(g.data(), sizeof(double) * g.size(), h);
    if (vir) { Eigen::Matrix3d xi = calc.virial(); h = fnv(xi.data(), sizeof(double) * 9, h); }
    return e;
}

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    if (argc < 6) { std::fprintf(stderr, "usage\n"); return 2; }
    const char* name = argv[1]; const char* method = argv[2];
    int nrep = std::atoi(argv[3]); std::string pre = argv[4]; int threads = std::atoi(argv[5]);
    json extra = argc > 6 ? json::parse(argv[6]) : json::object();
    double L = argc > 7 ? std::atof(argv[7]) : 0.0;
    curcuma::Molecule m = TestMolecules::TestMoleculeRegistry::createMolecule(name, false);
    Mol mol = m.getMolInfo();
    if (L > 0) { mol.m_unit_cell = Eigen::Matrix3d::Identity() * L; mol.m_has_pbc = true; }
    fpstate("start");
    uint64_t h;
    if (pre == "cpu" || pre == "gpu") {
        const int pt = std::getenv("REPRO_PRE_THREADS") ? std::atoi(std::getenv("REPRO_PRE_THREADS")) : threads;
        double e = runOnce(mol, method, pre == "gpu", pt, extra, h);
        std::printf("PRE %s %.17g %016llx\n", pre.c_str(), e, (unsigned long long)h);
        fpstate("after_pre");
    }
    const bool gpu = !(std::getenv("REPRO_CPU_SERIES"));
    for (int i = 0; i < nrep; ++i) {
        double e = runOnce(mol, method, gpu, threads, extra, h);
        std::printf("RUN %d %.17g %016llx\n", i, e, (unsigned long long)h);
        std::fflush(stdout);
    }
    fpstate("end");
    return 0;
}
