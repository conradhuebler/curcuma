/*
 * <Native QM engine -- density fitting (RI-J / RI-JK), WP10>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * See qm_ri.h for the equations and literature. Contents:
 *   1. Turbomole-format auxiliary basis parser and file resolution
 *   2. Shells, Hermite expansion, Coulomb auxiliary R, pair-pair kernel
 *   3. Harmonic (pure) transform for any L
 *   4. RIJK: 3-index tensor, metric, factorisation, J/K builds
 *   5. QMEngine glue: -qm.ri parameters, buildJK dispatch (exact path untouched)
 *
 * Claude Generated (Oct 2026, WP10). AI-generated, machine-tested.
 *
 * This program is free software under GPL-3.0
 */

#include "qm_ri.h"
#include "qm_engine.h"
#include "qm_citations.h"
#include "qm_basis_embedded.h"  // built-in fitting sets (Nov 2026)
#include "qm_elements.h"
#include "qm_gpu_ri_backend.h"

#include "src/core/blas_threads.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/gpu_plugin.h"
#include "src/core/gpu_fallback.h"
#include "src/core/units.h"

#include <Eigen/Dense>
#include <fmt/format.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace qmri {

namespace {
constexpr double PI = 3.14159265358979323846264338327950288;

double msSince(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

std::string lower(std::string s)
{
    for (char& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

bool fileExists(const std::string& p)
{
    std::ifstream f(p);
    return f.good();
}
}  // namespace

// =================================================================================
// 1. Auxiliary basis file
// =================================================================================

// Turbomole format (Basis Set Exchange export "turbomole"):
//   $jbas            (or $jkbas, $basis ...)
//   *
//   h def2-universal-JFIT
//   *
//       3   s
//        15.6752927              0.0186886
//        ...
//   *
//   he ...
//   $end
// Lines starting with '#' are comments. Coefficients refer to normalised primitives.
// Claude Generated (Nov 2026): `path` may be the pseudo path "[embedded]/<name>.turbomole" of a built-in set
// (resolveAuxBasisFile); its text goes through the same stream parser as a file.
AuxBasisMap parseTurbomoleBasis(const std::string& path, int max_z)
{
    if (qmbasis::isEmbeddedPath(path)) {
        std::string text;
        if (!qmbasis::embeddedText(path, text, max_z)) throw std::runtime_error("no built-in text for " + path);
        std::istringstream in(text);
        return parseTurbomoleStream(in, path);
    }
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open auxiliary basis file " + path);
    AuxBasisMap map = parseTurbomoleStream(in, path);
    // Claude Generated (Nov 2026): elements the file lacks (the H-Ne .turbomole files of the source directory) come
    // from the built-in set of the same name, if there is one; the file wins for every element it has.
    std::string stem = path.substr(path.find_last_of('/') == std::string::npos ? 0 : path.find_last_of('/') + 1);
    if (stem.size() > 10 && stem.compare(stem.size() - 10, 10, ".turbomole") == 0) stem.resize(stem.size() - 10);
    if (const auto* e = qmbasis::findEmbedded(stem, "turbomole"); e && max_z > 10) {
        std::istringstream bin(qmbasis::embeddedSetText(*e, max_z));
        for (auto& [el, shells] : parseTurbomoleStream(bin, path + " (built-in fill-in)"))
            if (!map.count(el)) map[el] = shells;
    }
    return map;
}

AuxBasisMap parseTurbomoleStream(std::istream& in, const std::string& path)
{
    AuxBasisMap map;
    std::string line, element;
    bool expect_element = false;
    const std::string letters = "spdfghi";
    int lineno = 0;
    while (std::getline(in, line)) {
        ++lineno;
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#') continue;
        std::string t = line.substr(first);
        if (t[0] == '$') {
            if (t.rfind("$end", 0) == 0) break;
            element.clear();
            expect_element = false;
            continue;
        }
        if (t[0] == '*') {
            // "*" closes the previous element (if any); the next header may follow
            element.clear();
            expect_element = true;
            continue;
        }
        std::istringstream ss(t);
        if (std::isalpha((unsigned char)t[0])) {
            if (!expect_element) throw std::runtime_error(fmt::format("{}:{}: element header without leading '*'", path, lineno));
            ss >> element;
            element = lower(element);
            map[element];  // create
            expect_element = false;
            // swallow the "*" that follows the header
            std::string next;
            while (std::getline(in, next)) {
                ++lineno;
                const auto f2 = next.find_first_not_of(" \t\r");
                if (f2 == std::string::npos || next[f2] == '#') continue;
                if (next[f2] != '*') throw std::runtime_error(fmt::format("{}:{}: '*' expected after element header", path, lineno));
                break;
            }
            continue;
        }
        if (element.empty()) throw std::runtime_error(fmt::format("{}:{}: shell outside an element block", path, lineno));
        expect_element = false;
        int nprim = 0;
        std::string lch;
        if (!(ss >> nprim >> lch) || lch.size() != 1 || letters.find(lch[0]) == std::string::npos || nprim <= 0)
            throw std::runtime_error(fmt::format("{}:{}: shell header 'n l' expected, got '{}'", path, lineno, t));
        AuxShellData sh;
        sh.L = (int)letters.find(lch[0]);
        for (int k = 0; k < nprim; ++k) {
            if (!std::getline(in, line)) throw std::runtime_error(path + ": file ends inside a shell");
            ++lineno;
            for (char& c : line)
                if (c == 'D' || c == 'd') c = 'E';  // Fortran exponents
            std::istringstream ps(line);
            double e = 0, c = 0;
            if (!(ps >> e >> c)) throw std::runtime_error(fmt::format("{}:{}: 'exponent coefficient' expected", path, lineno));
            sh.exps.push_back(e);
            sh.coefs.push_back(c);
        }
        map[element].push_back(sh);
    }
    if (map.empty()) throw std::runtime_error("no basis data in " + path);
    return map;
}

std::string resolveAuxBasisFile(const std::string& name_in, const std::string& mode)
{
    std::string name = name_in;
    const std::string ln = lower(name);
    if (ln.empty() || ln == "auto") name = (mode == "j") ? "def2-universal-JFIT" : "def2-universal-JKFIT";
    else if (ln == "def2/j" || ln == "def2-universal-jfit" || ln == "weigend" || ln == "def2-jfit") name = "def2-universal-JFIT";
    else if (ln == "def2/jk" || ln == "def2-universal-jkfit" || ln == "def2-jkfit") name = "def2-universal-JKFIT";
    else {
        // Claude Generated (Nov 2026): ORCA names of the RI-MP2 fitting sets, def2-<X>/C = def2-<X>-RIFIT (the same
        // numbers for H-Ne, checked by scripts/qm_basis_export_tables.py --compare-orca), def2-mTZVPP/J =
        // def2-mTZVPP-RIJ, and any built-in fitting set name in other letter case.
        static const char* const kOrcaC[][2] = {
            { "def2-svp/c", "def2-SVP-RIFIT" }, { "def2-svpd/c", "def2-SVPD-RIFIT" },
            { "def2-tzvp/c", "def2-TZVP-RIFIT" }, { "def2-tzvpp/c", "def2-TZVPP-RIFIT" },
            { "def2-tzvpd/c", "def2-TZVPD-RIFIT" }, { "def2-tzvppd/c", "def2-TZVPPD-RIFIT" },
            { "def2-qzvpp/c", "def2-QZVPP-RIFIT" }, { "def2-qzvppd/c", "def2-QZVPPD-RIFIT" },
            { "def2-mtzvpp/j", "def2-mTZVPP-RIJ" } };
        for (const auto& a : kOrcaC)
            if (ln == a[0]) name = a[1];
        if (name == name_in && name.find('/') == std::string::npos)
            if (const auto* e = qmbasis::findEmbeddedNoCase(name, "turbomole")) name = e->name;
    }
    // Claude Generated (Nov 2026): the built-in copy (qm_basis_embedded.h) comes LAST, a file always wins;
    // CURCUMA_QM_BASIS_EMBEDDED_ONLY=1 skips the file search (test hook).
    const bool verbose = CurcumaLogger::get_verbosity() >= 2;
    const std::string file = name + ".turbomole";
    std::vector<std::string> cand;
    if (!qmbasis::embeddedOnly()) {
        if (name.find('/') != std::string::npos && fileExists(name)) {
            qmcite::auxBasis(name);
            if (verbose) CurcumaLogger::info("QM: auxiliary basis " + name + ": file " + name);
            return name;
        }
        for (const char* var : { "CURCUMA_QM_BASIS", "CURCUMA_DFT_BASIS", "CURCUMA_DATA" })
            if (const char* env = std::getenv(var); env && env[0] != '\0') cand.push_back(std::string(env) + "/" + file);
#ifdef CURCUMA_SOURCE_DIR
        cand.push_back(std::string(CURCUMA_SOURCE_DIR) + "/src/core/energy_calculators/qm_methods/" + file);
#endif
        cand.push_back(file);
        for (const auto& c : cand)
            if (fileExists(c)) {
                qmcite::auxBasis(name);   // literature of the fitting basis actually loaded
                if (verbose) CurcumaLogger::info("QM: auxiliary basis " + name + ": file " + c);
                return c;
            }
    }
    if (const auto* e = qmbasis::findEmbedded(name, "turbomole")) {
        qmcite::auxBasis(name);
        if (verbose) CurcumaLogger::info("QM: auxiliary basis " + name + ": embedded copy (built into the binary)");
        return qmbasis::embeddedPath(*e, ".turbomole");
    }
    std::string looked;
    for (const auto& c : cand) looked += (looked.empty() ? "" : ", ") + c;
    throw std::runtime_error(fmt::format(
        "auxiliary basis '{}' not found: looked for {} (search path: explicit path, $CURCUMA_QM_BASIS, $CURCUMA_DATA, "
        "the source directory, the current directory){}; built in are {}",
        name_in, looked.empty() ? file : looked,
        qmbasis::embeddedOnly() ? " -- file search disabled by CURCUMA_QM_BASIS_EMBEDDED_ONLY" : "",
        qmbasis::embeddedNames("turbomole")));
}

// =================================================================================
// 2. Shells and the McMurchie-Davidson pair-pair kernel
// =================================================================================

Shell::Shell() = default;
Shell::Shell(const Shell&) = default;
Shell::Shell(Shell&&) noexcept = default;
Shell& Shell::operator=(const Shell&) = default;
Shell& Shell::operator=(Shell&&) noexcept = default;
Shell::~Shell() = default;

std::vector<Shell> orbitalShells(const std::vector<GTO::Orbital>& basis)
{
    // Same grouping rule as qm_integrals.cpp groupShells: consecutive AOs on one centre with
    // the same exponents and L form a shell (at most (L+1)(L+2)/2 components).
    std::vector<Shell> shells;
    for (size_t i = 0; i < basis.size(); ++i) {
        int l, m, k;
        GTO::orbitalTypeToComponents(basis[i].type, l, m, k);
        const int L = l + m + k;
        bool join = false;
        if (!shells.empty()) {
            const Shell& s = shells.back();
            const GTO::Orbital& prev = basis[i - 1];
            join = s.L == L && prev.atom == basis[i].atom && prev.exponents == basis[i].exponents && prev.x == basis[i].x
                && prev.y == basis[i].y && prev.z == basis[i].z && (int)s.lmn.size() < (L + 1) * (L + 2) / 2;
        }
        if (!join) {
            Shell s;
            s.first = (int)i;
            s.L = L;
            s.atom = basis[i].atom;
            s.x = basis[i].x;
            s.y = basis[i].y;
            s.z = basis[i].z;
            s.exps = basis[i].exponents;
            shells.push_back(std::move(s));
        }
        shells.back().lmn.push_back({ l, m, k });
        shells.back().coef.push_back(basis[i].coefficients);
    }
    return shells;
}

namespace {
// lower-case symbol, the key of the Turbomole aux parser (qm_elements.h; Claude Generated, Nov 2026: Z = 1..103)
std::string auxKey(int Z)
{
    std::string s = qmelem::symbol(Z);
    for (char& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

// Cartesian components of one shell in PySCF order: x-power descending, then y-power descending.
std::vector<std::array<int, 3>> cartesianComponents(int L)
{
    std::vector<std::array<int, 3>> c;
    for (int a = L; a >= 0; --a)
        for (int b = L - a; b >= 0; --b) c.push_back({ a, b, L - a - b });
    return c;
}

// Primitive normalisation of x^L exp(-alpha r^2) (the radial factor shared by all components).
double radialNorm(double alpha, int L)
{
    double df = 1.0;  // (2L-1)!!
    for (int k = 2 * L - 1; k > 1; k -= 2) df *= k;
    return std::pow(2.0 * alpha / PI, 0.75) * std::pow(4.0 * alpha, 0.5 * L) / std::sqrt(df);
}
}  // namespace

std::vector<Shell> auxiliaryShells(const AuxBasisMap& aux, const std::vector<int>& atomZ, const Matrix& pos)
{
    std::vector<Shell> shells;
    int first = 0;
    for (size_t a = 0; a < atomZ.size(); ++a) {
        const int Z = atomZ[a];
        if (qmelem::symbol(Z) == nullptr) throw std::runtime_error(fmt::format("RI: element Z={} outside 1..103", Z));
        auto it = aux.find(auxKey(Z));
        if (it == aux.end()) throw std::runtime_error(fmt::format("RI: no auxiliary functions for element {}", qmelem::symbol(Z)));
        for (const AuxShellData& d : it->second) {
            Shell s;
            s.x = pos(a, 0);
            s.y = pos(a, 1);
            s.z = pos(a, 2);
            s.L = d.L;
            s.atom = (int)a;
            s.first = first;
            s.exps = d.exps;
            s.lmn = cartesianComponents(d.L);
            std::vector<double> c(d.exps.size());
            for (size_t k = 0; k < d.exps.size(); ++k) c[k] = d.coefs[k] * radialNorm(d.exps[k], d.L);
            s.coef.assign(s.lmn.size(), c);
            first += (int)s.lmn.size();
            shells.push_back(std::move(s));
        }
    }
    return shells;
}

namespace {

// Hermite expansion coefficients E^{ij}_t along one axis (Helgaker 9.5.5/9.5.6), flat
// [t][i][j] with i <= iA, j <= iB, t <= iA+iB; E^{00}_0 = K (Gaussian product prefactor).
std::vector<double> hermite(int iA, int iB, double PA, double PB, double p, double K)
{
    const int tmax = iA + iB, si = iB + 1, st = (iA + 1) * (iB + 1);
    std::vector<double> E((size_t)(tmax + 1) * st, 0.0);
    auto at = [&](int t, int i, int j) -> double {
        if (t < 0 || t > tmax || i < 0 || j < 0) return 0.0;
        return E[(size_t)t * st + i * si + j];
    };
    const double g2 = 0.5 / p;
    E[0] = K;
    for (int i = 1; i <= iA; ++i)
        for (int t = 0; t <= tmax; ++t)
            E[(size_t)t * st + i * si] = g2 * at(t - 1, i - 1, 0) + PA * at(t, i - 1, 0) + (t + 1) * at(t + 1, i - 1, 0);
    for (int j = 1; j <= iB; ++j)
        for (int i = 0; i <= iA; ++i)
            for (int t = 0; t <= tmax; ++t)
                E[(size_t)t * st + i * si + j] = g2 * at(t - 1, i, j - 1) + PB * at(t, i, j - 1) + (t + 1) * at(t + 1, i, j - 1);
    return E;
}

// R^0_{tuv} on the simplex t+u+v <= L (Helgaker 9.9.13 with the two-electron rho):
// R^n_000 = (-2 rho)^n F_n(T); R^n_{t+1,u,v} = t R^{n+1}_{t-1,u,v} + W_x R^{n+1}_{t,u,v}, W = P - Q.
// Layout R0[(t*(L+1) + u)*(L+1) + v]. Same recursion as qm_integrals.cpp buildR0Simplex.
void coulombR0(int L, double Wx, double Wy, double Wz, double rho, const std::vector<double>& F, std::vector<double>& R0,
               std::vector<double>& lvl)
{
    const int s = L + 1;
    const size_t sz = (size_t)s * s * s;
    R0.resize(sz);
    lvl.resize(sz);
    auto at = [s](int t, int u, int v) { return ((size_t)t * s + u) * s + v; };
    double bpow[32];
    bpow[0] = 1.0;
    for (int n = 1; n <= L; ++n) bpow[n] = bpow[n - 1] * (-2.0 * rho);
    std::vector<double>* cur = &R0;
    std::vector<double>* nxt = &lvl;
    (*cur)[0] = bpow[L] * F[L];
    for (int n = L - 1; n >= 0; --n) {
        std::vector<double>& c = *cur;
        std::vector<double>& x = *nxt;
        x[0] = bpow[n] * F[n];
        const int top = L - n;
        for (int t = 0; t <= top; ++t)
            for (int u = 0; u <= top - t; ++u)
                for (int v = 0; v <= top - t - u; ++v) {
                    if (t + u + v == 0) continue;
                    double r;
                    if (v > 0) r = (v >= 2 ? (v - 1) * c[at(t, u, v - 2)] : 0.0) + Wz * c[at(t, u, v - 1)];
                    else if (u > 0) r = (u >= 2 ? (u - 1) * c[at(t, u - 2, v)] : 0.0) + Wy * c[at(t, u - 1, v)];
                    else r = (t >= 2 ? (t - 1) * c[at(t - 2, u, v)] : 0.0) + Wx * c[at(t - 1, u, v)];
                    x[at(t, u, v)] = r;
                }
        std::swap(cur, nxt);
    }
    if (cur != &R0) R0.swap(*cur);
}

struct PrimPair {
    int ia = 0, ib = 0;
    double p = 0, Px = 0, Py = 0, Pz = 0;
    std::vector<double> Ex, Ey, Ez;  ///< flat Hermite tables, K folded into Ex
    // Special members out of line: the implicit inline ones trip -Winline (RI gradient, Oct 2026).
    PrimPair();
    PrimPair(PrimPair&&) noexcept;
    ~PrimPair();
};
PrimPair::PrimPair() = default;
PrimPair::PrimPair(PrimPair&&) noexcept = default;
PrimPair::~PrimPair() = default;

/// A shell pair (A, B). An auxiliary function is the pair (C, unit s function on C).
struct Pair {
    const Shell* A = nullptr;
    const Shell* B = nullptr;
    std::vector<PrimPair> prims;
    double schwarz = 0.0;  ///< sqrt(max_comp |(ab|ab)|)
};

Pair makePair(const Shell& A, const Shell& B)
{
    Pair pr;
    pr.A = &A;
    pr.B = &B;
    for (size_t ia = 0; ia < A.exps.size(); ++ia)
        for (size_t ib = 0; ib < B.exps.size(); ++ib) {
            const double a = A.exps[ia], b = B.exps[ib];
            PrimPair pp;
            pp.ia = (int)ia;
            pp.ib = (int)ib;
            pp.p = a + b;
            pp.Px = (a * A.x + b * B.x) / pp.p;
            pp.Py = (a * A.y + b * B.y) / pp.p;
            pp.Pz = (a * A.z + b * B.z) / pp.p;
            const double ab2 = (A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) + (A.z - B.z) * (A.z - B.z);
            const double K = std::exp(-a * b / pp.p * ab2);
            if (K == 0.0) continue;
            pp.Ex = hermite(A.L, B.L, pp.Px - A.x, pp.Px - B.x, pp.p, K);
            pp.Ey = hermite(A.L, B.L, pp.Py - A.y, pp.Py - B.y, pp.p, 1.0);
            pp.Ez = hermite(A.L, B.L, pp.Pz - A.z, pp.Pz - B.z, pp.p, 1.0);
            pr.prims.push_back(std::move(pp));
        }
    return pr;
}

Shell unitShell(const Shell& C)
{
    Shell u;
    u.x = C.x;
    u.y = C.y;
    u.z = C.z;
    u.L = 0;
    u.atom = C.atom;
    u.exps = { 0.0 };
    u.lmn = { { 0, 0, 0 } };
    u.coef = { { 1.0 } };
    return u;
}

/// Contracted (ab|cd) for every component quartet of the two pairs, out[(ca*nB + cb)*nCD + cc*nD + cd]:
///   (ab|cd) = sum_prim c_a c_b c_c c_d 2 pi^(5/2) / (p q sqrt(p+q))
///             sum_{tuv} E^ab_t E^ab_u E^ab_v sum_{tau nu phi} (-1)^(tau+nu+phi) E^cd_tau E^cd_nu E^cd_phi R_{t+tau,u+nu,v+phi}
/// (Helgaker 9.9.33). The ket sum is formed once per ket component pair into G_cd(t,u,v) and
/// accumulated over the ket primitives before the bra contraction.
void pairBlock(const Pair& bra, const Pair& ket, std::vector<double>& out)
{
    const Shell& A = *bra.A;
    const Shell& B = *bra.B;
    const Shell& C = *ket.A;
    const Shell& D = *ket.B;
    const int nA = (int)A.lmn.size(), nB = (int)B.lmn.size(), nC = (int)C.lmn.size(), nD = (int)D.lmn.size();
    const int nCD = nC * nD;
    out.assign((size_t)nA * nB * nCD, 0.0);
    const int LAB = A.L + B.L, LCD = C.L + D.L, Ltot = LAB + LCD;
    const int sAB_t = (A.L + 1) * (B.L + 1), sAB_i = B.L + 1;
    const int sCD_t = (C.L + 1) * (D.L + 1), sCD_i = D.L + 1;
    const int s = Ltot + 1;
    auto at = [s](int t, int u, int v) { return ((size_t)t * s + u) * s + v; };
    const int g = LAB + 1;
    const size_t gsz = (size_t)g * g * g;
    thread_local std::vector<double> F, R, lvl, G;
    G.resize((size_t)nCD * gsz);

    for (const PrimPair& pb : bra.prims) {
        std::fill(G.begin(), G.end(), 0.0);
        for (const PrimPair& pk : ket.prims) {
            const double p = pb.p, q = pk.p, rho = p * q / (p + q);
            const double Wx = pb.Px - pk.Px, Wy = pb.Py - pk.Py, Wz = pb.Pz - pk.Pz;
            qmint::boysFunction(Ltot, rho * (Wx * Wx + Wy * Wy + Wz * Wz), F);
            coulombR0(Ltot, Wx, Wy, Wz, rho, F, R, lvl);
            const double pref = 2.0 * std::pow(PI, 2.5) / (p * q * std::sqrt(p + q));
            for (int cc = 0; cc < nC; ++cc) {
                const auto& lc = C.lmn[cc];
                const double wc = pref * C.coef[cc][pk.ia];
                for (int cd = 0; cd < nD; ++cd) {
                    const double w = wc * D.coef[cd][pk.ib];
                    if (w == 0.0) continue;
                    const auto& ld = D.lmn[cd];
                    const int tx = lc[0] + ld[0], ty = lc[1] + ld[1], tz = lc[2] + ld[2];
                    double ex[16], ey[16], ez[16];
                    for (int k = 0; k <= tx; ++k) ex[k] = ((k & 1) ? -1.0 : 1.0) * pk.Ex[k * sCD_t + lc[0] * sCD_i + ld[0]];
                    for (int k = 0; k <= ty; ++k) ey[k] = ((k & 1) ? -1.0 : 1.0) * pk.Ey[k * sCD_t + lc[1] * sCD_i + ld[1]];
                    for (int k = 0; k <= tz; ++k) ez[k] = ((k & 1) ? -1.0 : 1.0) * pk.Ez[k * sCD_t + lc[2] * sCD_i + ld[2]];
                    double* Gcd = &G[(size_t)(cc * nD + cd) * gsz];
                    for (int t = 0; t <= LAB; ++t)
                        for (int u = 0; u <= LAB - t; ++u)
                            for (int v = 0; v <= LAB - t - u; ++v) {
                                double acc = 0.0;
                                for (int k1 = 0; k1 <= tx; ++k1) {
                                    if (ex[k1] == 0.0) continue;
                                    double a2 = 0.0;
                                    for (int k2 = 0; k2 <= ty; ++k2) {
                                        if (ey[k2] == 0.0) continue;
                                        const double* Rr = &R[at(t + k1, u + k2, v)];
                                        double a3 = 0.0;
                                        for (int k3 = 0; k3 <= tz; ++k3) a3 += ez[k3] * Rr[k3];
                                        a2 += ey[k2] * a3;
                                    }
                                    acc += ex[k1] * a2;
                                }
                                Gcd[((size_t)t * g + u) * g + v] += w * acc;
                            }
                }
            }
        }
        for (int ca = 0; ca < nA; ++ca) {
            const auto& la = A.lmn[ca];
            const double c_a = A.coef[ca][pb.ia];
            for (int cb = 0; cb < nB; ++cb) {
                const double c_ab = c_a * B.coef[cb][pb.ib];
                if (c_ab == 0.0) continue;
                const auto& lb = B.lmn[cb];
                const int tx = la[0] + lb[0], ty = la[1] + lb[1], tz = la[2] + lb[2];
                double ex[16], ey[16], ez[16];
                for (int k = 0; k <= tx; ++k) ex[k] = pb.Ex[k * sAB_t + la[0] * sAB_i + lb[0]];
                for (int k = 0; k <= ty; ++k) ey[k] = pb.Ey[k * sAB_t + la[1] * sAB_i + lb[1]];
                for (int k = 0; k <= tz; ++k) ez[k] = pb.Ez[k * sAB_t + la[2] * sAB_i + lb[2]];
                double* o = &out[(size_t)(ca * nB + cb) * nCD];
                for (int k = 0; k < nCD; ++k) {
                    const double* Gk = &G[(size_t)k * gsz];
                    double val = 0.0;
                    for (int t = 0; t <= tx; ++t) {
                        if (ex[t] == 0.0) continue;
                        for (int u = 0; u <= ty; ++u) {
                            const double etu = ex[t] * ey[u];
                            if (etu == 0.0) continue;
                            const double* row = Gk + ((size_t)t * g + u) * g;
                            double sv = 0.0;
                            for (int v = 0; v <= tz; ++v) sv += ez[v] * row[v];
                            val += etu * sv;
                        }
                    }
                    o[k] += c_ab * val;
                }
            }
        }
    }
}

double schwarzOf(const Pair& pr)
{
    std::vector<double> buf;
    pairBlock(pr, pr, buf);
    const int nA = (int)pr.A->lmn.size(), nB = (int)pr.B->lmn.size(), nAB = nA * nB;
    double mx = 0.0;
    for (int k = 0; k < nAB; ++k) mx = std::max(mx, std::abs(buf[(size_t)k * nAB + k]));
    return std::sqrt(mx);
}

/// Contract axis k (dims d[0..2]) of a 3-index block with T (d[k] x m), in place.
void transformAxis3(std::vector<double>& blk, std::vector<double>& tmp, int d[3], int k, const Eigen::MatrixXd& T)
{
    size_t outer = 1, inner = 1;
    for (int m = 0; m < k; ++m) outer *= d[m];
    for (int m = k + 1; m < 3; ++m) inner *= d[m];
    const int nin = d[k], nout = (int)T.cols();
    tmp.assign(outer * nout * inner, 0.0);
    for (size_t o = 0; o < outer; ++o)
        for (int a = 0; a < nin; ++a)
            for (int j = 0; j < nout; ++j) {
                const double t = T(a, j);
                if (t == 0.0) continue;
                const double* x = &blk[(o * nin + a) * inner];
                double* y = &tmp[(o * nout + j) * inner];
                for (size_t r = 0; r < inner; ++r) y[r] += t * x[r];
            }
    blk.swap(tmp);
    d[k] = nout;
}

// Dense kernels on the row-major tensor (N rows of length ncol = n^2). Each thread owns a range of
// columns for all rows and runs ONE level-3 call on it (Eigen; with -DUSE_BLAS Eigen hands it to the
// linked BLAS, which is single-threaded MKL in the release build, so the threading is ours). The
// column-major view of a column chunk is X = T(:, c0:c0+w)^T (w x N, leading dimension ncol), so a left
// operation on T is a right operation on X. Up to Oct 2026 these were scalar loops (factor/B was ~half
// of the GPU RI-JK wall time); the blocked form changes only the summation order (last bits).
// Claude Generated (Oct 2026).
constexpr Eigen::Index kRowChunk = 256;
using StridedMap = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>;

/// In place T := L^-1 T (L lower triangular N x N)  <=>  X := X L^-T per column chunk (dtrsm).
void lowerSolveRows(const Eigen::MatrixXd& L, Matrix& T, int threads)
{
    const Eigen::Index N = T.rows(), ncol = T.cols(), nchunks = (ncol + kRowChunk - 1) / kRowChunk;
    double* base = T.data();
#pragma omp parallel for schedule(dynamic) num_threads(std::max(1, threads))
    for (Eigen::Index b = 0; b < nchunks; ++b) {
        const Eigen::Index c0 = b * kRowChunk, w = std::min(kRowChunk, ncol - c0);
        StridedMap X(base + c0, w, N, Eigen::OuterStride<>(ncol));
        L.transpose().triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheRight>(X);
    }
}

/// out := F T (F: r x N, T: N x ncol row-major)  <=>  Out_chunk = X F^T (dgemm).
void leftMultiplyRows(const Eigen::MatrixXd& F, const Matrix& T, Matrix& out, int threads)
{
    const Eigen::Index N = T.rows(), r = F.rows(), ncol = T.cols(), nchunks = (ncol + kRowChunk - 1) / kRowChunk;
    out.resize(r, ncol);
    const Eigen::MatrixXd Ft = F.transpose();
#pragma omp parallel for schedule(dynamic) num_threads(std::max(1, threads))
    for (Eigen::Index b = 0; b < nchunks; ++b) {
        const Eigen::Index c0 = b * kRowChunk, w = std::min(kRowChunk, ncol - c0);
        Eigen::Map<const Eigen::MatrixXd, 0, Eigen::OuterStride<>> X(T.data() + c0, w, N, Eigen::OuterStride<>(ncol));
        StridedMap O(out.data() + c0, w, r, Eigen::OuterStride<>(ncol));
        O.noalias() = X * Ft;
    }
}

/// Where each shell lands in the output basis: first index, width, transform (identity -> empty).
struct OutMap {
    int n = 0;
    std::vector<int> first, width;
    std::vector<Eigen::MatrixXd> T;  ///< empty = identity
    OutMap();
    OutMap(const OutMap&);
    OutMap& operator=(const OutMap&);
    ~OutMap();  // out of line: the implicit inline one trips -Winline
};
OutMap::OutMap() = default;
OutMap::OutMap(const OutMap&) = default;
OutMap& OutMap::operator=(const OutMap&) = default;
OutMap::~OutMap() = default;

OutMap orbitalOutMap(const std::vector<Shell>& sh, int ncart, const Matrix* Q)
{
    OutMap o;
    const bool sph = Q != nullptr && Q->size() != 0;
    o.n = sph ? (int)Q->cols() : ncart;
    for (const Shell& s : sh) {
        const int nc = (int)s.lmn.size();
        if (!sph) {
            o.first.push_back(s.first);
            o.width.push_back(nc);
            o.T.emplace_back();
            continue;
        }
        int lo = o.n, hi = -1;
        for (int j = 0; j < o.n; ++j)
            for (int a = 0; a < nc; ++a)
                if ((*Q)(s.first + a, j) != 0.0) {
                    lo = std::min(lo, j);
                    hi = std::max(hi, j);
                }
        Eigen::MatrixXd T = Eigen::MatrixXd::Zero(nc, hi - lo + 1);
        for (int a = 0; a < nc; ++a)
            for (int j = 0; j < T.cols(); ++j) T(a, j) = (*Q)(s.first + a, lo + j);
        o.first.push_back(lo);
        o.width.push_back((int)T.cols());
        if (T.rows() == T.cols() && T.isIdentity(0.0)) T.resize(0, 0);
        o.T.push_back(T);
    }
    return o;
}

OutMap auxOutMap(const std::vector<Shell>& sh, bool spherical)
{
    OutMap o;
    for (const Shell& s : sh) {
        o.first.push_back(o.n);
        if (spherical && s.L >= 2) {
            o.T.push_back(harmonicTransform(s.L));
            o.width.push_back(2 * s.L + 1);
        } else {
            o.T.emplace_back();
            o.width.push_back((int)s.lmn.size());
        }
        o.n += o.width.back();
    }
    return o;
}

}  // namespace

// =================================================================================
// 3. Harmonic transform
// =================================================================================

// The pure functions of degree L are r^L Y_lm exp(-alpha r^2) = h(x,y,z) exp(-alpha r^2) with h a
// homogeneous HARMONIC polynomial (Laplacian h = 0). With all cartesian components sharing one
// radial factor, a coefficient vector over the components is a coefficient vector of h, so the
// pure subspace is the null space of the Laplacian map from degree L to degree L-2 monomials:
//   Lap(x^a y^b z^c) = a(a-1) x^(a-2) y^b z^c + b(b-1) x^a y^(b-2) z^c + c(c-1) x^a y^b z^(c-2).
// Dimension (L+1)(L+2)/2 - L(L-1)/2 = 2L+1. The null-space basis from the SVD is orthonormal.
Eigen::MatrixXd harmonicTransform(int L)
{
    const auto comps = cartesianComponents(L);
    const int nc = (int)comps.size();
    if (L <= 1) return Eigen::MatrixXd::Identity(nc, nc);
    const auto low = cartesianComponents(L - 2);
    auto index = [&](int a, int b, int c) {
        for (size_t k = 0; k < low.size(); ++k)
            if (low[k][0] == a && low[k][1] == b && low[k][2] == c) return (int)k;
        return -1;
    };
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero((int)low.size(), nc);
    for (int k = 0; k < nc; ++k) {
        const int a = comps[k][0], b = comps[k][1], c = comps[k][2];
        if (a >= 2) M(index(a - 2, b, c), k) += a * (a - 1);
        if (b >= 2) M(index(a, b - 2, c), k) += b * (b - 1);
        if (c >= 2) M(index(a, b, c - 2), k) += c * (c - 1);
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullV);
    // M has full row rank L(L-1)/2; the last 2L+1 right singular vectors span its null space.
    return svd.matrixV().rightCols(2 * L + 1);
}

Matrix coulombMetric(const std::vector<Shell>& aux, bool spherical, int threads)
{
    const OutMap o = auxOutMap(aux, spherical);
    std::vector<Shell> units;
    units.reserve(aux.size());
    for (const Shell& s : aux) units.push_back(unitShell(s));
    std::vector<Pair> pairs;
    pairs.reserve(aux.size());
    for (size_t c = 0; c < aux.size(); ++c) pairs.push_back(makePair(aux[c], units[c]));
    Matrix J = Matrix::Zero(o.n, o.n);
    const int ns = (int)aux.size();
#pragma omp parallel for schedule(dynamic) num_threads(std::max(1, threads))
    for (int c = 0; c < ns; ++c) {
        std::vector<double> buf, tmp;
        for (int d = c; d < ns; ++d) {
            pairBlock(pairs[c], pairs[d], buf);
            int dims[3] = { (int)aux[c].lmn.size(), 1, (int)aux[d].lmn.size() };
            if (o.T[c].size()) transformAxis3(buf, tmp, dims, 0, o.T[c]);
            if (o.T[d].size()) transformAxis3(buf, tmp, dims, 2, o.T[d]);
            for (int i = 0; i < dims[0]; ++i)
                for (int j = 0; j < dims[2]; ++j) {
                    const double v = buf[(size_t)i * dims[2] + j];
                    J(o.first[c] + i, o.first[d] + j) = v;
                    J(o.first[d] + j, o.first[c] + i) = v;
                }
        }
    }
    return J;
}

// =================================================================================
// 4. RIJK
// =================================================================================

// Packed lower triangle -> full symmetric block (Claude Generated, Oct 2026): the lower rows are contiguous
// in the packed vector (memcpy), the upper triangle is mirrored in 32 x 32 tiles (cache-friendly transpose).
void unpackSymmetric(const double* p, int n, double* f)
{
    for (int m = 0; m < n; ++m) std::memcpy(f + (size_t)m * n, p + (size_t)m * (m + 1) / 2, sizeof(double) * (size_t)(m + 1));
    constexpr int bs = 32;
    for (int i0 = 0; i0 < n; i0 += bs)
        for (int j0 = 0; j0 <= i0; j0 += bs) {
            const int i1 = std::min(n, i0 + bs), j1 = std::min(n, j0 + bs);
            for (int i = i0; i < i1; ++i)
                for (int j = j0, je = std::min(j1, i); j < je; ++j) f[(size_t)j * n + i] = f[(size_t)i * n + j];
        }
}

void packContractionWeights(const Matrix& P, std::vector<double>& w)
{
    const int n = (int)P.rows();
    w.resize((size_t)n * (n + 1) / 2);
    for (int m = 0; m < n; ++m) {
        double* wm = w.data() + (size_t)m * (m + 1) / 2;
        for (int k = 0; k < m; ++k) wm[k] = P(m, k) + P(k, m);
        wm[m] = P(m, m);
    }
}

// Shells, pairs and out-maps of one geometry, kept by RIJK when rows have to be recomputed (batched storage, FP32
// setup in column chunks, the FP64 phase after mixed precision). The Pair objects point into osh / units, so the
// vectors are filled once (reserved) and never reallocated. Claude Generated (Oct 2026).
struct RIJK::Direct {
    std::vector<Shell> osh, units;
    OutMap oo, ao;
    std::vector<Pair> opairs;                 ///< canonical orbital shell pairs A <= B not screened for every aux shell
    std::vector<std::pair<int, int>> oidx;    ///< (A, B) of opairs
    std::vector<Pair> apairs;                 ///< (aux shell, unit s) pairs
};

void RIJK::fillRows(const std::vector<int>* pairs, int c0, int c1, double* out, Eigen::Index ld, Eigen::Index col0, long& blocks,
                    long& skipped) const
{
    const Direct& d = *m_dd;
    const int threads = std::max(1, m_opt.threads);
    const int nk = pairs ? (int)pairs->size() : (int)d.opairs.size();
    const int row0 = d.ao.first[c0];
    const double scr = m_opt.screening;
    long bl = 0, sk = 0;
#pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+ : bl, sk)
    for (int kk = 0; kk < nk; ++kk) {
        const int k = pairs ? (*pairs)[kk] : kk;
        const int A = d.oidx[k].first, B = d.oidx[k].second;
        thread_local std::vector<double> buf, tmp;
        for (int c = c0; c < c1; ++c) {
            if (d.opairs[k].schwarz * d.apairs[c].schwarz < scr) {
                ++sk;
                continue;
            }
            ++bl;
            pairBlock(d.opairs[k], d.apairs[c], buf);
            int dm[3] = { (int)d.osh[A].lmn.size(), (int)d.osh[B].lmn.size(), (int)m_aux[c].lmn.size() };
            if (d.oo.T[A].size()) transformAxis3(buf, tmp, dm, 0, d.oo.T[A]);
            if (d.oo.T[B].size()) transformAxis3(buf, tmp, dm, 1, d.oo.T[B]);
            if (d.ao.T[c].size()) transformAxis3(buf, tmp, dm, 2, d.ao.T[c]);
            for (int a = 0; a < dm[0]; ++a)
                for (int b = 0; b < dm[1]; ++b) {
                    const int m = d.oo.first[A] + a, nn = d.oo.first[B] + b;
                    if (A == B && m < nn) continue;
                    const double* v = &buf[((size_t)a * dm[1] + b) * dm[2]];
                    const Eigen::Index col = pairIndex(m, nn) - col0;
                    double* o = out + (size_t)(d.ao.first[c] - row0) * ld + col;
                    for (int p = 0; p < dm[2]; ++p) o[(size_t)p * ld] = v[p];
                }
        }
    }
    blocks += bl;
    skipped += sk;
}

Eigen::VectorXd RIJK::metricSolve(const Eigen::VectorXd& g) const
{
    if (m_gL.size()) {
        Eigen::VectorXd y = m_gL.triangularView<Eigen::Lower>().solve(g);
        return m_gL.transpose().triangularView<Eigen::Upper>().solve(y);
    }
    if (m_gU.size()) return m_gU * (m_gw.cwiseInverse().asDiagonal() * (m_gU.transpose() * g));
    throw std::logic_error("RIJK::metricSolve: no metric factor");
}

std::string RIJK::storageUsed() const
{
    if (m_batched) return "batched";
    return m_T32.empty() ? "stored" : "stored FP32";
}

bool RIJK::setup(const std::vector<GTO::Orbital>& basis, const Matrix* Q, const std::vector<Shell>& aux,
                 const RIOptions& opt, std::string& why)
{
    m_opt = opt;
    m_aux = aux;
    m_gL.resize(0, 0);
    m_gU.resize(0, 0);
    m_gw.resize(0);
    m_T.resize(0, 0);
    std::vector<float>().swap(m_T32);
    m_batched = false;
    m_dd.reset();
    m_aux_batches.clear();
    curcuma::ScopedBlasThreads blas_guard(std::max(1, opt.threads));
    const auto t0 = std::chrono::steady_clock::now();
    auto dd = std::make_shared<Direct>();
    dd->osh = orbitalShells(basis);
    dd->oo = orbitalOutMap(dd->osh, (int)basis.size(), Q);
    dd->ao = auxOutMap(aux, opt.spherical_aux);
    const std::vector<Shell>& osh = dd->osh;
    const OutMap& oo = dd->oo;
    const OutMap& ao = dd->ao;
    m_n = oo.n;
    m_naux = ao.n;
    const int n = m_n, N = m_naux;
    const Eigen::Index np = npair();  // packed columns n(n+1)/2 (Oct 2026; before n^2)
    const int threads = std::max(1, opt.threads);
    const int naux_sh = (int)aux.size();

    // Storage (Claude Generated, Oct 2026): stored FP64 as before; FP32 (mode jk, ri_mixed_precision); batched (no
    // tensor). auto takes stored when it fits the limit, else batched.
    const bool jk = opt.mode == "jk";
    const bool eig = opt.metric == "eigen";
    const std::string storage = lower(opt.storage);
    if (storage != "stored" && storage != "batched" && storage != "auto") {
        why = "unknown -qm.ri_storage '" + opt.storage + "' (use stored|batched|auto)";
        return false;
    }
    const double need_mb = (double)N * np * 8.0 / 1.0e6 * (eig && jk ? 2.0 : 1.0);
    // FP32 setup in column chunks: FP32 tensor + one FP64 chunk (two with eigen)
    const Eigen::Index chunk_cols = std::max<Eigen::Index>(1, (Eigen::Index)(256.0e6 / (8.0 * std::max(1, N))));
    const double need32_mb = (double)N * np * 4.0 / 1.0e6 + (double)N * std::min(np, chunk_cols) * 8.0 / 1.0e6 * (eig ? 2.0 : 1.0);
    bool batched = storage == "batched" || (opt.fit_radius > 0.0 && jk), fp32 = false;  // local fit: always batched
    if (!batched) {
        if (opt.fp32 && jk) {
            if (need32_mb <= opt.max_memory_mb)
                fp32 = true;
            else if (storage == "auto")
                batched = true;
            else {
                why = fmt::format("FP32 3-index tensor needs {:.0f} MB (naux {} x n(n+1)/2 {}, packed), above -qm.ri_max_memory_mb {:g}",
                                  need32_mb, N, (long)np, opt.max_memory_mb);
                return false;
            }
        } else if (need_mb > opt.max_memory_mb) {
            if (storage == "auto")
                batched = true;
            else {
                why = fmt::format("3-index tensor needs {:.0f} MB (naux {} x n(n+1)/2 {}, packed), above -qm.ri_max_memory_mb {:g}", need_mb,
                                  N, (long)np, opt.max_memory_mb);
                return false;
            }
        }
    }

    // Pairs and Schwarz factors
    std::vector<Pair> opairs;
    std::vector<std::pair<int, int>> oidx;
    for (int A = 0; A < (int)osh.size(); ++A)
        for (int B = A; B < (int)osh.size(); ++B) {
            opairs.push_back(makePair(osh[A], osh[B]));
            oidx.push_back({ A, B });
        }
    dd->units.reserve(aux.size());
    for (const Shell& s : aux) dd->units.push_back(unitShell(s));
    std::vector<Pair>& apairs = dd->apairs;
    apairs.reserve(aux.size());
    for (size_t c = 0; c < aux.size(); ++c) apairs.push_back(makePair(m_aux[c], dd->units[c]));  // m_aux: kept with the object
    const int npairs = (int)opairs.size();
#pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int k = 0; k < npairs; ++k) opairs[k].schwarz = schwarzOf(opairs[k]);
    for (Pair& p : apairs) p.schwarz = schwarzOf(p);
    m_osch.resize(npairs);
    for (int k = 0; k < npairs; ++k) m_osch[k] = opairs[k].schwarz;
    m_asch.resize(naux_sh);
    for (int c = 0; c < naux_sh; ++c) m_asch[c] = apairs[c].schwarz;
    // Keep only the pairs that survive the screen for at least one aux shell (the others are skipped for every c,
    // so the tensor does not change); their count goes into the skipped blocks of the statistics.
    double amax = 0.0;
    for (const Pair& p : apairs) amax = std::max(amax, p.schwarz);
    long pre_skipped = 0;
    {
        int kept = 0;
        for (int k = 0; k < npairs; ++k)
            if (opairs[k].schwarz * amax >= opt.screening) ++kept;
        dd->opairs.reserve(kept);
        dd->oidx.reserve(kept);
        for (int k = 0; k < npairs; ++k) {
            if (opairs[k].schwarz * amax >= opt.screening) {
                dd->opairs.push_back(std::move(opairs[k]));
                dd->oidx.push_back(oidx[k]);
            } else
                pre_skipped += naux_sh;
        }
        std::vector<Pair>().swap(opairs);
    }
    m_dd = dd;
    const double t_schwarz = msSince(t0);

    // Metric (raw) and the Coulomb normalisation: every aux function scaled to (P|P) = 1. The fit is invariant to
    // it; it only puts the metric eigenvalues on a common scale for the eigen cutoff and the report.
    const auto t2 = std::chrono::steady_clock::now();
    Matrix Jraw = coulombMetric(aux, opt.spherical_aux, threads);
    const double t_2c = msSince(t2);
    Eigen::VectorXd dn(N);
    for (int P = 0; P < N; ++P) {
        if (!(Jraw(P, P) > 0.0)) {
            why = fmt::format("auxiliary function {} has (P|P) = {:.3e}", P, Jraw(P, P));
            return false;
        }
        dn(P) = 1.0 / std::sqrt(Jraw(P, P));
    }
    m_dn = dn;

    // 3-index tensor T(P, pairIndex(m, nn)) = (P|m nn) in the active orbital basis, lower triangle only. In a
    // diagonal shell block (A == B) both (a,b) and (b,a) are computed; the one with m >= nn is kept, which is
    // the value the full n^2 layout kept for both entries (last write in loop order), so the packed tensor
    // holds the same numbers. Stored storage only; batched / FP32 compute the rows later (fillRows).
    const auto t1 = std::chrono::steady_clock::now();
    Matrix T3;
    long blocks = 0, skipped = pre_skipped;
    if (!batched && !fp32) {
        T3 = Matrix::Zero(N, np);
        fillRows(nullptr, 0, naux_sh, T3.data(), np, 0, blocks, skipped);
    }
    double t_3c = msSince(t1);
    if (!batched && !fp32) {
        if (const char* dump = std::getenv("CURCUMA_QM_RI_DUMP"); dump && dump[0] != '\0') {
            m_raw_metric = Jraw;
            m_raw_T3.resize(N, (Eigen::Index)n * n);  // full n^2 columns for the dump (diagnostic)
            for (int P = 0; P < N; ++P) unpackSymmetric(T3.data() + (size_t)P * np, n, m_raw_T3.data() + (size_t)P * n * n);
        }
        for (int P = 0; P < N; ++P) T3.row(P) *= dn(P);
    }
    Eigen::MatrixXd Jn = dn.asDiagonal() * Eigen::MatrixXd(Jraw) * dn.asDiagonal();
    Matrix().swap(Jraw);
    m_Jn.resize(0, 0);
    m_dom_cache.clear();
    if (opt.fit_radius > 0.0 && jk) {
        // local fit: metric blocks of the domains, function -> atom maps, atom positions (from the shell centres)
        m_Jn = Jn;
        int natoms = 0;
        for (const Shell& s : aux) natoms = std::max(natoms, s.atom + 1);
        for (const Shell& s : osh) natoms = std::max(natoms, s.atom + 1);
        m_atom_xyz = Eigen::MatrixXd::Zero(natoms, 3);
        for (const Shell& s : aux) m_atom_xyz.row(s.atom) << s.x, s.y, s.z;
        for (const Shell& s : osh) m_atom_xyz.row(s.atom) << s.x, s.y, s.z;
        m_ao_atom.assign(n, 0);
        for (size_t s = 0; s < osh.size(); ++s)
            for (int j = 0; j < oo.width[s]; ++j) m_ao_atom[oo.first[s] + j] = osh[s].atom;
        m_auxfn_atom.assign(N, 0);
        for (int c = 0; c < naux_sh; ++c)
            for (int j = 0; j < ao.width[c]; ++j) m_auxfn_atom[ao.first[c] + j] = aux[c].atom;
    }

    const auto t3 = std::chrono::steady_clock::now();
    std::string metric = opt.metric;
    if (metric != "cholesky" && metric != "eigen") {
        why = "unknown -qm.ri_metric '" + metric + "' (use cholesky|eigen)";
        return false;
    }
    double wmin = std::numeric_limits<double>::quiet_NaN(), wmax = wmin;
    double t_chol = 0.0;  // Cholesky of the metric alone (diagnostic split of factor/B)
    std::string solve_where = "host";
    int dropped = 0;
    m_metric_used = metric;
    Eigen::MatrixXd Ffit;  // eigen + (batched or FP32): w^-1/2 U^T
    if (metric == "cholesky") {
        Eigen::LLT<Eigen::MatrixXd> llt(Jn);
        if (llt.info() != Eigen::Success) {
            metric = "eigen";
            m_metric_used = "eigen (Cholesky failed)";
        } else {
            m_L = llt.matrixL();
            m_gL = m_L;  // gradient: F^T = L^-T in both modes
            t_chol = msSince(t3);
            if (opt.mode == "jk" && !batched && !fp32) {
                // B = L^-1 T (triangular solve, no explicit inverse): on the device if the caller asks for it
                // and it works, else on the host. Both are the same solve; only the summation order differs.
                int dev = 0;
                if (opt.device_lower_solve) {
                    std::string dwhy;
                    dev = opt.device_lower_solve(m_L, T3, dwhy);
                    if (dev < 0) {
                        why = "device metric solve failed after writing part of the tensor: " + dwhy;
                        return false;
                    }
                    solve_where = dev == 1 ? "device" : "host (device failed: " + dwhy + ")";
                }
                if (dev != 1) lowerSolveRows(m_L, T3, threads);
            }
            if (opt.mode == "jk") m_L.resize(0, 0);
            m_T.swap(T3);
            m_nfit = N;
        }
    }
    if (metric == "eigen") {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Jn);
        if (es.info() != Eigen::Success) {
            why = "eigendecomposition of the RI metric failed";
            return false;
        }
        const Eigen::VectorXd& w = es.eigenvalues();
        wmin = w(0);
        wmax = w(N - 1);
        std::vector<int> keep;
        for (int k = 0; k < N; ++k)
            if (w(k) > opt.metric_cutoff) keep.push_back(k);
        dropped = N - (int)keep.size();
        Eigen::MatrixXd U((Eigen::Index)N, (Eigen::Index)keep.size());
        Eigen::VectorXd wk((Eigen::Index)keep.size());
        for (size_t k = 0; k < keep.size(); ++k) {
            U.col(k) = es.eigenvectors().col(keep[k]);
            wk(k) = w(keep[k]);
        }
        m_nfit = (int)keep.size();
        m_gU = U;  // gradient: F^T = U w^-1/2 in both modes
        m_gw = wk;
        if (opt.mode == "jk") {
            // B = w^-1/2 U^T T  (J^-1 ~ U w^-1 U^T = (w^-1/2 U^T)^T (w^-1/2 U^T))
            const Eigen::MatrixXd F = wk.cwiseSqrt().cwiseInverse().asDiagonal() * U.transpose();
            if (!batched && !fp32) {
                leftMultiplyRows(F, T3, m_T, threads);
                Matrix().swap(T3);
            } else
                Ffit = F;
        } else {
            m_U = U;
            m_w = wk;
            m_T.swap(T3);
        }
    }
    double t_fac = msSince(t3);

    // Batched rows (also prepared for FP32, whose final SCF iterations run batched): aux shell batches of at most
    // batch_mb, and the memory check of the batched builds (Claude Generated, Oct 2026).
    std::string storage_info = "stored";
    if (batched || fp32) {
        const double cap_mb = opt.batch_mb > 0.0 ? opt.batch_mb : std::min(256.0, opt.max_memory_mb / 8.0);
        const double row_mb = (double)np * 8.0 / 1.0e6;
        long rows = 0, maxrows = 0;
        int c0 = 0;
        for (int c = 0; c < naux_sh; ++c) {
            const int w = ao.width[c];
            if (c > c0 && (rows + w) * row_mb > cap_mb) {
                m_aux_batches.push_back({ c0, c });
                maxrows = std::max(maxrows, rows);
                c0 = c;
                rows = 0;
            }
            rows += w;
        }
        m_aux_batches.push_back({ c0, naux_sh });
        m_batch_rows = std::max(maxrows, rows);
        // metric factor (N^2) + row batch + per-thread K (n^2) + per-thread packed J for FP32 builds
        m_fixed_mb = ((double)N * N + (double)m_batch_rows * np + (double)threads * n * n) * 8.0 / 1.0e6;
        const double ycol_mb = (double)N * n * 8.0 / 1.0e6 * (Ffit.size() ? 2.0 : 1.0);  // Y per orbital column (eigen: + X)
        double held_mb = fp32 ? need32_mb : 0.0;  // the FP32 tensor stays next to the batched buffers
        if (fp32 && storage == "auto" && m_fixed_mb + ycol_mb + held_mb > opt.max_memory_mb) {
            fp32 = false;  // FP32 tensor plus the buffers of its exact builds do not fit: batched FP64 instead
            batched = true;
            held_mb = 0.0;
        }
        if (jk && m_fixed_mb + ycol_mb + held_mb > opt.max_memory_mb) {
            why = fmt::format("batched RI-JK needs at least {:.0f} MB (metric factor, a {}-row batch of 3-index integrals, per-thread K, "
                              "one orbital column of Y = {:.0f} MB{}), above -qm.ri_max_memory_mb {:g}",
                              m_fixed_mb + ycol_mb + held_mb, m_batch_rows, ycol_mb,
                              fp32 ? fmt::format(", FP32 tensor {:.0f} MB", need32_mb) : std::string(), opt.max_memory_mb);
            return false;
        }
        storage_info = fmt::format("batched ({} batches of <= {} rows, {:.0f} MB each{})", m_aux_batches.size(), m_batch_rows,
                                   m_batch_rows * row_mb,
                                   opt.fit_radius > 0.0 && jk ? fmt::format("; local exchange fit R = {:g} Bohr", opt.fit_radius) : std::string());
        m_batched = batched;
    }

    // FP32 tensor (mode jk): column chunks of whole orbital rows m (shells B), each computed in FP64, normalised,
    // fitted (L^-1 or w^-1/2 U^T, columns are independent) and rounded to FP32. Peak = FP32 tensor + one chunk.
    if (fp32) {
        const auto tf = std::chrono::steady_clock::now();
        const int ns = (int)osh.size();
        for (int s = 0; s + 1 < ns; ++s)
            if (oo.first[s + 1] != oo.first[s] + oo.width[s]) {
                why = "FP32 RI setup needs contiguous orbital shells";
                return false;
            }
        std::vector<std::vector<int>> byB(ns);
        for (int k = 0; k < (int)dd->oidx.size(); ++k) byB[dd->oidx[k].second].push_back(k);
        m_T32.assign((size_t)m_nfit * np, 0.0f);
        Matrix C, Cf;
        double t_int = 0.0;
        for (int s0 = 0; s0 < ns;) {
            int s1 = s0;
            Eigen::Index f0 = oo.first[s0], cols = 0;
            std::vector<int> plist;
            do {
                plist.insert(plist.end(), byB[s1].begin(), byB[s1].end());
                ++s1;
                const Eigen::Index f1 = oo.first[s1 - 1] + oo.width[s1 - 1];
                cols = f1 * (f1 + 1) / 2 - f0 * (f0 + 1) / 2;
            } while (s1 < ns && cols < chunk_cols);
            const Eigen::Index col0 = f0 * (f0 + 1) / 2;
            C = Matrix::Zero(N, cols);
            const auto ti = std::chrono::steady_clock::now();
            fillRows(&plist, 0, naux_sh, C.data(), cols, col0, blocks, skipped);
            t_int += msSince(ti);
            for (int P = 0; P < N; ++P) C.row(P) *= dn(P);
            const Matrix* src = &C;
            if (m_gL.size()) {
                lowerSolveRows(m_gL, C, threads);
            } else {
                leftMultiplyRows(Ffit, C, Cf, threads);
                src = &Cf;
            }
#pragma omp parallel for schedule(static) num_threads(threads)
            for (int Qr = 0; Qr < m_nfit; ++Qr) {
                const double* a = src->data() + (size_t)Qr * cols;
                float* b = m_T32.data() + (size_t)Qr * np + col0;
                for (Eigen::Index j = 0; j < cols; ++j) b[j] = (float)a[j];
            }
            s0 = s1;
        }
        t_3c = t_int;
        t_fac += msSince(tf) - t_int;
        storage_info = fmt::format("stored FP32 ({:.1f} MB; exact FP64 corrections {})", (double)m_T32.size() * 4.0 / 1.0e6, storage_info);
    }
    if (!batched && !fp32) m_dd.reset();  // stored: the pair tables are not needed after the setup
    if (CurcumaLogger::get_verbosity() >= 3 && std::isnan(wmin)) {
        // report the conditioning (diagnostic only; the numerical path does not depend on it)
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Jn, Eigen::EigenvaluesOnly);
        wmin = es.eigenvalues()(0);
        wmax = es.eigenvalues()(N - 1);
    }
    const std::string cond = std::isnan(wmin)
        ? std::string("metric eigenvalues at verbosity 3")
        : fmt::format("eigenvalues of the Coulomb-normalised metric {:.3e} .. {:.3e}, cond {:.3e}{}", wmin, wmax, wmax / wmin,
                      dropped ? fmt::format(", {} dropped below {:.1e}", dropped, opt.metric_cutoff) : std::string());
    m_stats = fmt::format(
        "naux {} (fit rows {}), n {}, tensor {:.1f} MB; metric {} ({}); setup: Schwarz {:.1f} ms, 3-index {:.1f} ms ({} of {} shell "
        "blocks), metric {:.1f} ms, factor/B {:.1f} ms{}; storage {}",
        N, m_nfit, n, tensorMB(), m_metric_used, cond, t_schwarz, t_3c, blocks, blocks + skipped, t_2c, t_fac,
        t_chol > 0.0 ? fmt::format(" (Cholesky {:.1f} ms, solve {})", t_chol, solve_where) : std::string(), storage_info);
    return true;
}

Eigen::VectorXd RIJK::fitCoefficients(const Eigen::VectorXd& g) const
{
    if (m_L.size()) {
        Eigen::VectorXd y = m_L.triangularView<Eigen::Lower>().solve(g);
        return m_L.transpose().triangularView<Eigen::Upper>().solve(y);
    }
    if (m_U.size()) return m_U * (m_w.cwiseInverse().asDiagonal() * (m_U.transpose() * g));
    return g;  // mode jk: the tensor already carries the metric (B), c = B . D
}

// Claude Generated (Oct 2026, RI-MP2): the same metric factor as the jk setup applies to the
// n^2 columns, here applied to the (usually much shorter) transformed rows.
void RIJK::fitRows(Matrix& X) const
{
    if (!storedFP64()) throw std::logic_error("RIJK::fitRows: needs the stored FP64 tensor (ri_storage " + storageUsed() + ")");
    if (m_opt.mode == "jk") return;  // tensor already fitted
    const int threads = std::max(1, m_opt.threads);
    if (m_L.size()) {
        lowerSolveRows(m_L, X, threads);
        return;
    }
    if (m_U.size()) {
        const Eigen::MatrixXd F = m_w.cwiseSqrt().cwiseInverse().asDiagonal() * m_U.transpose();
        Matrix out;
        leftMultiplyRows(F, X, out, threads);
        X.swap(out);
    }
}

Matrix RIJK::coulomb(const Matrix& P) const
{
    if (m_batched || !m_T32.empty()) {  // Claude Generated (Oct 2026): no FP64 tensor in memory
        Matrix J;
        if (m_batched)
            buildBatched(P, J, nullptr, nullptr);
        else
            buildFP32(P, J, nullptr, nullptr);
        return J;
    }
    const int n = m_n, threads = std::max(1, m_opt.threads);
    std::vector<double> Pw;  // packed contraction weights (P_mm, P_mk + P_km)
    packContractionWeights(P, Pw);
    const Eigen::Index ncol = npair();
    // g_Q = sum_mn T_Q,mn P_mn (one dot product per row). Rows of the tensor, not m_nfit: in mode j
    // with the eigen metric the tensor keeps all naux rows while m_nfit counts the kept eigenvalues
    // (fixed Oct 2026; before, a dropped eigenvalue in mode j read past g).
    const int rows = (int)m_T.rows();
    Eigen::VectorXd g(rows);
#pragma omp parallel for schedule(static) num_threads(threads)
    for (int Q = 0; Q < rows; ++Q) {
        const double* t = m_T.data() + (size_t)Q * ncol;
        double s = 0.0;
        for (Eigen::Index k = 0; k < ncol; ++k) s += t[k] * Pw[k];
        g(Q) = s;
    }
    const Eigen::VectorXd c = fitCoefficients(g);
    // J_p = sum_Q c_Q T_Q,p over the packed pairs (column chunks in parallel, Q summed in order), then unpacked
    std::vector<double> Jp((size_t)ncol);
    const Eigen::Index chunk = 1024, nchunks = (ncol + chunk - 1) / chunk;
#pragma omp parallel for schedule(static) num_threads(threads)
    for (Eigen::Index b = 0; b < nchunks; ++b) {
        const Eigen::Index c0 = b * chunk, w = std::min(chunk, ncol - c0);
        double* out = Jp.data() + c0;
        for (Eigen::Index k = 0; k < w; ++k) out[k] = 0.0;
        for (int Q = 0; Q < rows; ++Q) {
            const double cq = c(Q);
            const double* t = m_T.data() + (size_t)Q * ncol + c0;
            for (Eigen::Index k = 0; k < w; ++k) out[k] += cq * t[k];
        }
    }
    Matrix J(n, n);
    unpackSymmetric(Jp.data(), n, J.data());
    return J;
}

void RIJK::factorDensity(const Matrix& P, Eigen::MatrixXd& V, int& npos)
{
    const Eigen::MatrixXd Pc(P);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Pc);
    const Eigen::VectorXd& w = es.eigenvalues();
    const double wmax = w.cwiseAbs().maxCoeff();
    std::vector<int> pos, neg;
    for (int k = 0; k < w.size(); ++k) {
        if (std::abs(w(k)) <= 1.0e-13 * wmax) continue;
        (w(k) > 0.0 ? pos : neg).push_back(k);
    }
    V.resize(Pc.rows(), (Eigen::Index)(pos.size() + neg.size()));
    int col = 0;
    for (int k : pos) V.col(col++) = es.eigenvectors().col(k) * std::sqrt(w(k));
    for (int k : neg) V.col(col++) = es.eigenvectors().col(k) * std::sqrt(-w(k));
    npos = (int)pos.size();
}

void RIJK::build(const Matrix& P, Matrix& J, Matrix* K, const Eigen::MatrixXd* Vocc) const
{
    if (!K) {
        J = coulomb(P);
        return;
    }
    if (m_opt.mode != "jk") throw std::logic_error("RIJK::build: exchange needs mode jk");
    if (m_batched) {  // Claude Generated (Oct 2026): recomputed integrals, K from Y = T V
        buildBatched(P, J, K, Vocc);
        return;
    }
    if (!m_T32.empty()) {  // Claude Generated (Oct 2026): FP32 tensor, FP32 GEMM/SYRK, FP64 accumulation
        buildFP32(P, J, K, Vocc);
        return;
    }
    const int n = m_n, threads = std::max(1, m_opt.threads);
    // Density factor P = V S V^T: the occupied orbitals (V = sqrt(2) C_occ, all signs +, given by the caller
    // when P = 2 C_occ C_occ^T exactly) or the eigen-factorisation of P (n^3 on the host).
    Eigen::MatrixXd Vown;
    const Eigen::MatrixXd* Vp = Vocc;
    int npos = 0;
    if (Vp && Vp->rows() == n && Vp->cols() > 0) {
        npos = (int)Vp->cols();
    } else {
        const auto tf = std::chrono::steady_clock::now();
        factorDensity(P, Vown, npos);
        m_factor_ms += msSince(tf);
        Vp = &Vown;
    }
    const Eigen::MatrixXd& V = *Vp;
    const int r = (int)V.cols(), nneg = r - npos;
    // One pass over B for J and K (Claude Generated, Oct 2026; before: J = coulomb(P), two more passes over the
    // naux x n^2 tensor, then K). Per batch of nq fit rows, per row q:
    //   X^q = B^q V   (n x r; B^q is stored as its packed lower triangle, copied into a full-storage buffer of
    //                  which only the triangle is valid, and multiplied as a symmetric matrix, dsymm)
    //   g_q = Tr(B^q P) = sum_k s_k v_k . X^q_{:,k}      (the Coulomb fit coefficient; mode jk: c = g)
    //   J_p += g_q B^q_p                                 (packed, while B^q is still in cache)
    // X^q is row block q of Y (nq n x r, column-major), then
    //   K  += Y+ Y+^T - Y- Y-^T                          (rank updates, dsyrk; the first npos columns of Y
    //                                                     read as one column-major n x (nq npos) matrix are
    //                                                     [.. (X^q)_{:,k} ..])
    // The same SYRK as the device (cuda/qm_gpu_ri.cu). Single-threaded BLAS calls, our own OpenMP over the
    // batches; per thread private packed J and lower-triangle K, batches statically distributed (deterministic
    // for a fixed thread count), summed in thread order. J equals coulomb(P) up to the summation order and the
    // eigenvalues factorDensity drops (|lambda| <= 1e-13 max|lambda|; none with the orbital factor).
    const Eigen::Index np = npair();
    const int nq = std::clamp(512 / std::max(1, r) + 1, 1, 64);  // SYRK depth nq * r >= 512
    const int nbatch = (m_nfit + nq - 1) / nq;
    std::vector<Eigen::MatrixXd> Kt(threads);
    std::vector<std::vector<double>> Jt(threads);
#pragma omp parallel num_threads(threads)
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        Eigen::MatrixXd& Kl = Kt[tid];
        std::vector<double>& Jl = Jt[tid];
        Kl = Eigen::MatrixXd::Zero(n, n);
        Jl.assign((size_t)np, 0.0);
        Eigen::MatrixXd Y, F = Eigen::MatrixXd::Zero(n, n);  // only the upper triangle is ever written/read
#pragma omp for schedule(static)
        for (int b = 0; b < nbatch; ++b) {
            const int Q0 = b * nq, cnt = std::min(nq, m_nfit - Q0);
            Y.resize((Eigen::Index)cnt * n, r);
            for (int q = 0; q < cnt; ++q) {
                const double* tq = m_T.data() + (size_t)(Q0 + q) * np;
                // column-major F(k, m) = B^q_mk for k <= m: the upper triangle of F holds the packed rows
                for (int m = 0; m < n; ++m) std::memcpy(F.data() + (size_t)m * n, tq + (size_t)m * (m + 1) / 2, sizeof(double) * (size_t)(m + 1));
                Y.middleRows((Eigen::Index)q * n, n).noalias() = F.selfadjointView<Eigen::Upper>() * V;
                double g = 0.0;
                for (int k = 0; k < r; ++k) {
                    const double d = V.col(k).dot(Y.col(k).segment((Eigen::Index)q * n, n));
                    g += k < npos ? d : -d;
                }
                double* jl = Jl.data();
                for (Eigen::Index p = 0; p < np; ++p) jl[p] += g * tq[p];
            }
            if (npos > 0) {
                Eigen::Map<const Eigen::MatrixXd> Yp(Y.data(), n, (Eigen::Index)cnt * npos);
                Kl.selfadjointView<Eigen::Lower>().rankUpdate(Yp, 1.0);
            }
            if (nneg > 0) {
                Eigen::Map<const Eigen::MatrixXd> Yn(Y.data() + (size_t)cnt * n * npos, n, (Eigen::Index)cnt * nneg);
                Kl.selfadjointView<Eigen::Lower>().rankUpdate(Yn, -1.0);
            }
        }
    }
    Eigen::MatrixXd Ks = Eigen::MatrixXd::Zero(n, n);
    std::vector<double> Js((size_t)np, 0.0);
    for (int t = 0; t < threads; ++t) {
        if (Kt[t].size())
            for (int q = 0; q < n; ++q)
                for (int m = q; m < n; ++m) Ks(m, q) += Kt[t](m, q);
        if (!Jt[t].empty())
            for (Eigen::Index p = 0; p < np; ++p) Js[p] += Jt[t][p];
    }
    Matrix Kf(n, n), Jf(n, n);
    for (int q = 0; q < n; ++q)
        for (int m = q; m < n; ++m) {
            Kf(m, q) = Ks(m, q);
            Kf(q, m) = Ks(m, q);
        }
    unpackSymmetric(Js.data(), n, Jf.data());
    J = Jf;
    *K = Kf;
}

// ---------------------------------------------------------------------------------
// 4a. Resource-saving J/K builds (Claude Generated, Oct 2026; docs WP10 "RI-JK ressourcensparend")
// ---------------------------------------------------------------------------------
//
// BATCHED (no stored tensor). With T_n the Coulomb-normalised raw rows (P|mn), J_n = L L^T and P = V S V^T:
//   g_P = T_n,P . w(P),   c = J_n^-1 g,   J_mn = sum_P c_P T_n,P,mn                 (= sum_Q B^Q (B^Q . P))
//   Y^P = T_n^P V (n x r),   X = L^-1 Y (auxiliary index),   K = sum_Q sum_k s_k X^Q_k X^Q_k^T  (= sum_Q B^Q P B^Q)
// The rows are recomputed in aux-shell batches (fillRows, bitwise the stored setup's values). Pass 0 forms g and
// Y for orbital batch 0, pass 1 forms J (c is known after pass 0) and Y for orbital batch 1, and so on: max(2,
// #orbital batches) passes. The orbital batch width rb follows -qm.ri_max_memory_mb (Y is naux n rb doubles), each
// batch carries one sign. B is never formed: a forward substitution over row batches would need all earlier rows
// of B (L is a full triangle), so batching B itself saves nothing; the occupied half-transform does.
// FP32 (mixed precision): the stored B in FP32, per row X^q = B^q V in FP32 (ssymm), J and g accumulated in FP64,
// per batch K_b = sum X X^T in FP32 (ssyrk) added into an FP64 K, so the FP32 rounding does not grow with naux.

std::string RIJK::takeBuildStats() const
{
    std::string s;
    if (m_bt_builds > 0)
        s = fmt::format("batched builds {}: passes {}, orbital batches {}; 3-index {:.0f} ms, Y {:.0f} ms, L^-1 Y {:.0f} ms, SYRK {:.0f} ms, "
                        "J {:.0f} ms{}",
                        m_bt_builds, m_bt_passes, m_bt_ybatches, m_bt_int_ms, m_bt_y_ms, m_bt_solve_ms, m_bt_syrk_ms, m_bt_j_ms,
                        m_bt_dev_ms > 0 ? fmt::format(", device solve+SYRK {:.0f} ms", m_bt_dev_ms) : std::string());
    if (m_fp_builds > 0)
        s += fmt::format("{}FP32 builds {} in {:.0f} ms", s.empty() ? "" : "; ", m_fp_builds, m_fp_ms);
    if (m_loc_builds > 0)
        s += fmt::format("; local fit R = {:g} Bohr: {} builds ({} with the global fit), mean domain {:.0f} of {} aux functions, {} new "
                         "domain factors; pivoted Cholesky + domains {:.0f} ms, Y {:.0f} ms, local solves + SYRK {:.0f} ms",
                         m_opt.fit_radius, m_loc_builds, m_loc_global, m_loc_orbitals ? (double)m_loc_aux_sum / m_loc_orbitals : 0.0, m_naux,
                         m_loc_new_domains, m_loc_factor_ms, m_loc_y_ms, m_loc_solve_ms);
    m_loc_builds = m_loc_global = m_loc_orbitals = m_loc_aux_sum = m_loc_new_domains = 0;
    m_loc_factor_ms = m_loc_y_ms = m_loc_solve_ms = 0.0;
    m_bt_builds = m_bt_passes = m_bt_ybatches = 0;
    m_bt_int_ms = m_bt_y_ms = m_bt_solve_ms = m_bt_syrk_ms = m_bt_j_ms = m_bt_dev_ms = 0.0;
    m_fp_builds = 0;
    m_fp_ms = 0.0;
    return s;
}

void RIJK::buildBatched(const Matrix& P, Matrix& J, Matrix* K, const Eigen::MatrixXd* Vocc) const
{
    using clk = std::chrono::steady_clock;
    const Direct& d = *m_dd;
    const int n = m_n, N = m_naux, threads = std::max(1, m_opt.threads);
    const int naux_sh = (int)m_aux.size();
    const Eigen::Index np = npair();
    if (K && m_opt.fit_radius > 0.0 && buildBatchedLocal(P, J, *K, Vocc)) {
        ++m_bt_builds;
        return;
    }
    if (K && m_opt.fit_radius > 0.0) ++m_loc_global;  // indefinite P or a domain metric not positive definite: global fit
    curcuma::ScopedBlasThreads blas_guard(1);  // our own OpenMP over single-threaded BLAS calls (as build())
    std::vector<double> Pw;
    packContractionWeights(P, Pw);

    // Density factor (as build()): sqrt(2) C_occ from the caller or the eigen-factorisation of P
    Eigen::MatrixXd Vown;
    const Eigen::MatrixXd* Vp = nullptr;
    int npos = 0;
    if (K) {
        if (Vocc && Vocc->rows() == n && Vocc->cols() > 0) {
            Vp = Vocc;
            npos = (int)Vocc->cols();
        } else {
            const auto tf = clk::now();
            factorDensity(P, Vown, npos);
            m_factor_ms += msSince(tf);
            Vp = &Vown;
        }
    }
    const int r = Vp ? (int)Vp->cols() : 0;
    const bool chol = m_gL.size() != 0;
    Eigen::MatrixXd F;  // eigen metric: X = F Y, F = w^-1/2 U^T
    if (K && !chol) F = m_gw.cwiseSqrt().cwiseInverse().asDiagonal() * m_gU.transpose();
    // Orbital batches with one sign each; width from the memory limit (Y: naux n rb, eigen also X: nfit n rb)
    struct OrbBatch {
        int k0, cnt;
        double sign;
    };
    std::vector<OrbBatch> obs;
    if (K && r > 0) {
        const double ycol_mb = (double)N * n * 8.0 / 1.0e6 * (chol ? 1.0 : 2.0);
        // the FP32 tensor (mixed precision) stays in memory next to Y
        const double held_mb = (double)m_T32.size() * 4.0 / 1.0e6;
        const int rb = std::clamp((int)std::floor((m_opt.max_memory_mb - m_fixed_mb - held_mb) / ycol_mb), 1, r);
        for (int k0 = 0; k0 < npos; k0 += rb) obs.push_back({ k0, std::min(rb, npos - k0), 1.0 });
        for (int k0 = npos; k0 < r; k0 += rb) obs.push_back({ k0, std::min(rb, r - k0), -1.0 });
    }
    const int nob = (int)obs.size();
    const int npass = std::max(2, nob);

    Eigen::VectorXd g = Eigen::VectorXd::Zero(N), c;
    std::vector<double> Jp((size_t)np, 0.0);
    Matrix Tb((Eigen::Index)m_batch_rows, np);
    Matrix Y, X;
    std::vector<Eigen::MatrixXd> Kt(K ? threads : 0);
    for (auto& k : Kt) k = Eigen::MatrixXd::Zero(n, n);
    Eigen::MatrixXd Kdev;  // device-accumulated part (full)
    long blocks = 0, skipped = 0;
    for (int pass = 0; pass < npass; ++pass) {
        const OrbBatch* ob = pass < nob ? &obs[pass] : nullptr;
        const bool doG = pass == 0, doJ = pass == 1;
        if (!doG && !doJ && !ob) continue;
        ++m_bt_passes;
        Eigen::MatrixXd Vb;
        if (ob) {
            Vb = Vp->middleCols(ob->k0, ob->cnt);
            Y.resize(N, (Eigen::Index)n * ob->cnt);
            ++m_bt_ybatches;
        }
        for (const auto& cb : m_aux_batches) {
            const int row0 = d.ao.first[cb.first];
            const int rows = (cb.second < naux_sh ? d.ao.first[cb.second] : N) - row0;
            const auto ti = clk::now();
            Tb.topRows(rows).setZero();
            fillRows(nullptr, cb.first, cb.second, Tb.data(), np, 0, blocks, skipped);
            for (int q = 0; q < rows; ++q) Tb.row(q) *= m_dn(row0 + q);
            m_bt_int_ms += msSince(ti);
            if (doG) {
#pragma omp parallel for schedule(static) num_threads(threads)
                for (int q = 0; q < rows; ++q) {
                    const double* t = Tb.data() + (size_t)q * np;
                    double s = 0.0;
                    for (Eigen::Index k = 0; k < np; ++k) s += t[k] * Pw[k];
                    g(row0 + q) = s;
                }
            }
            if (doJ) {
                const auto tj = clk::now();
                const Eigen::Index chunk = 1024, nchunks = (np + chunk - 1) / chunk;
#pragma omp parallel for schedule(static) num_threads(threads)
                for (Eigen::Index b = 0; b < nchunks; ++b) {
                    const Eigen::Index c0 = b * chunk, w = std::min(chunk, np - c0);
                    double* out = Jp.data() + c0;
                    for (int q = 0; q < rows; ++q) {
                        const double cq = c(row0 + q);
                        const double* t = Tb.data() + (size_t)q * np + c0;
                        for (Eigen::Index k = 0; k < w; ++k) out[k] += cq * t[k];
                    }
                }
                m_bt_j_ms += msSince(tj);
            }
            if (ob) {
                const auto ty = clk::now();
#pragma omp parallel num_threads(threads)
                {
                    Eigen::MatrixXd Fu = Eigen::MatrixXd::Zero(n, n);  // only the upper triangle is written/read
#pragma omp for schedule(static)
                    for (int q = 0; q < rows; ++q) {
                        const double* tq = Tb.data() + (size_t)q * np;
                        for (int m = 0; m < n; ++m)
                            std::memcpy(Fu.data() + (size_t)m * n, tq + (size_t)m * (m + 1) / 2, sizeof(double) * (size_t)(m + 1));
                        Eigen::Map<Eigen::MatrixXd> Yq(Y.data() + (size_t)(row0 + q) * n * ob->cnt, n, ob->cnt);
                        Yq.noalias() = Fu.selfadjointView<Eigen::Upper>() * Vb;
                    }
                }
                m_bt_y_ms += msSince(ty);
            }
        }
        if (doG) c = metricSolve(g);
        if (!ob) continue;
        // X = L^-1 Y (or F Y) and the rank updates; the device hook does both for the Cholesky metric
        if (chol && m_dev_exchange) {
            const auto td = clk::now();
            if (Kdev.size() == 0) Kdev = Eigen::MatrixXd::Zero(n, n);
            if (m_dev_exchange(Y, ob->cnt, ob->sign, Kdev)) {
                m_bt_dev_ms += msSince(td);
                continue;
            }
            m_dev_exchange = nullptr;  // device error: this batch and the rest on the CPU (the caller counted it)
            Kdev.setZero();            // a failed call may have added part of the batch
            // earlier batches done on the device are lost: redo the build on the CPU
            Kdev.resize(0, 0);
            buildBatched(P, J, K, Vocc);
            return;
        }
        const auto ts = clk::now();
        const Matrix* Xs = &Y;
        if (chol)
            lowerSolveRows(m_gL, Y, threads);
        else {
            leftMultiplyRows(F, Y, X, threads);
            Xs = &X;
        }
        m_bt_solve_ms += msSince(ts);
        const auto tk = clk::now();
        const int cnt = ob->cnt, rowsX = (int)Xs->rows();
        const int nq = std::clamp(512 / std::max(1, cnt) + 1, 1, 64);
        const int nbatch = (rowsX + nq - 1) / nq;
#pragma omp parallel num_threads(threads)
        {
#ifdef _OPENMP
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
#pragma omp for schedule(static)
            for (int b = 0; b < nbatch; ++b) {
                const int Q0 = b * nq, qn = std::min(nq, rowsX - Q0);
                Eigen::Map<const Eigen::MatrixXd> Xb(Xs->data() + (size_t)Q0 * n * cnt, n, (Eigen::Index)qn * cnt);
                Kt[tid].selfadjointView<Eigen::Lower>().rankUpdate(Xb, ob->sign);
            }
        }
        m_bt_syrk_ms += msSince(tk);
    }
    J.resize(n, n);
    unpackSymmetric(Jp.data(), n, J.data());
    if (K) {
        Eigen::MatrixXd Ks = Eigen::MatrixXd::Zero(n, n);
        for (int t = 0; t < (int)Kt.size(); ++t)
            for (int q = 0; q < n; ++q)
                for (int m = q; m < n; ++m) Ks(m, q) += Kt[t](m, q);
        if (Kdev.size())
            for (int q = 0; q < n; ++q)
                for (int m = q; m < n; ++m) Ks(m, q) += Kdev(m, q);
        K->resize(n, n);
        for (int q = 0; q < n; ++q)
            for (int m = q; m < n; ++m) {
                (*K)(m, q) = Ks(m, q);
                (*K)(q, m) = Ks(m, q);
            }
    }
    ++m_bt_builds;
}

void RIJK::buildFP32(const Matrix& P, Matrix& J, Matrix* K, const Eigen::MatrixXd* Vocc) const
{
    const auto t0 = std::chrono::steady_clock::now();
    const int n = m_n, threads = std::max(1, m_opt.threads);
    const Eigen::Index np = npair();
    curcuma::ScopedBlasThreads blas_guard(1);
    Eigen::MatrixXd Vown;
    const Eigen::MatrixXd* Vp = Vocc;
    int npos = 0;
    if (!K) {
        // Coulomb only: g_Q = B_Q . w, J = sum_Q g_Q B_Q (FP64 accumulation of FP32 rows)
        std::vector<double> Pw;
        packContractionWeights(P, Pw);
        Eigen::VectorXd g(m_nfit);
#pragma omp parallel for schedule(static) num_threads(threads)
        for (int Q = 0; Q < m_nfit; ++Q) {
            const float* t = m_T32.data() + (size_t)Q * np;
            double s = 0.0;
            for (Eigen::Index k = 0; k < np; ++k) s += (double)t[k] * Pw[k];
            g(Q) = s;
        }
        std::vector<double> Jp((size_t)np, 0.0);
        const Eigen::Index chunk = 1024, nchunks = (np + chunk - 1) / chunk;
#pragma omp parallel for schedule(static) num_threads(threads)
        for (Eigen::Index b = 0; b < nchunks; ++b) {
            const Eigen::Index c0 = b * chunk, w = std::min(chunk, np - c0);
            double* out = Jp.data() + c0;
            for (int Q = 0; Q < m_nfit; ++Q) {
                const double gq = g(Q);
                const float* t = m_T32.data() + (size_t)Q * np + c0;
                for (Eigen::Index k = 0; k < w; ++k) out[k] += gq * (double)t[k];
            }
        }
        J.resize(n, n);
        unpackSymmetric(Jp.data(), n, J.data());
        ++m_fp_builds;
        m_fp_ms += msSince(t0);
        return;
    }
    if (Vp && Vp->rows() == n && Vp->cols() > 0) {
        npos = (int)Vp->cols();
    } else {
        const auto tf = std::chrono::steady_clock::now();
        factorDensity(P, Vown, npos);
        m_factor_ms += msSince(tf);
        Vp = &Vown;
    }
    const Eigen::MatrixXf Vf = Vp->cast<float>();
    const Eigen::MatrixXd& Vd = *Vp;
    const int r = (int)Vf.cols(), nneg = r - npos;
    const int nq = std::clamp(512 / std::max(1, r) + 1, 1, 64);
    const int nbatch = (m_nfit + nq - 1) / nq;
    std::vector<Eigen::MatrixXd> Kt(threads);
    std::vector<std::vector<double>> Jt(threads);
#pragma omp parallel num_threads(threads)
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        Eigen::MatrixXd& Kl = Kt[tid];
        std::vector<double>& Jl = Jt[tid];
        Kl = Eigen::MatrixXd::Zero(n, n);
        Jl.assign((size_t)np, 0.0);
        Eigen::MatrixXf Y, Ff, Kb;
        Eigen::MatrixXd Yd, Fd;
        if (m_fp32_arith) {
            Ff = Eigen::MatrixXf::Zero(n, n);
            Kb.resize(n, n);
        } else
            Fd = Eigen::MatrixXd::Zero(n, n);
#pragma omp for schedule(static)
        for (int b = 0; b < nbatch; ++b) {
            if (!m_fp32_arith) {
                // FP32 storage, FP64 arithmetic (second phase): the rows widened to FP64, then exactly the stored build
                const int Q0 = b * nq, cnt = std::min(nq, m_nfit - Q0);
                Yd.resize((Eigen::Index)cnt * n, r);
                for (int q = 0; q < cnt; ++q) {
                    const float* tq = m_T32.data() + (size_t)(Q0 + q) * np;
                    for (int m = 0; m < n; ++m) {
                        const float* src = tq + (size_t)m * (m + 1) / 2;
                        double* dst = Fd.data() + (size_t)m * n;
                        for (int k = 0; k <= m; ++k) dst[k] = (double)src[k];
                    }
                    Yd.middleRows((Eigen::Index)q * n, n).noalias() = Fd.selfadjointView<Eigen::Upper>() * Vd;
                    double gq = 0.0;
                    for (int k = 0; k < r; ++k) {
                        const double d = Vd.col(k).dot(Yd.col(k).segment((Eigen::Index)q * n, n));
                        gq += k < npos ? d : -d;
                    }
                    double* jl = Jl.data();
                    for (Eigen::Index p = 0; p < np; ++p) jl[p] += gq * (double)tq[p];
                }
                if (npos > 0) {
                    Eigen::Map<const Eigen::MatrixXd> Yp(Yd.data(), n, (Eigen::Index)cnt * npos);
                    Kl.selfadjointView<Eigen::Lower>().rankUpdate(Yp, 1.0);
                }
                if (nneg > 0) {
                    Eigen::Map<const Eigen::MatrixXd> Yn(Yd.data() + (size_t)cnt * n * npos, n, (Eigen::Index)cnt * nneg);
                    Kl.selfadjointView<Eigen::Lower>().rankUpdate(Yn, -1.0);
                }
                continue;
            }
            const int Q0 = b * nq, cnt = std::min(nq, m_nfit - Q0);
            Y.resize((Eigen::Index)cnt * n, r);
            for (int q = 0; q < cnt; ++q) {
                const float* tq = m_T32.data() + (size_t)(Q0 + q) * np;
                for (int m = 0; m < n; ++m) std::memcpy(Ff.data() + (size_t)m * n, tq + (size_t)m * (m + 1) / 2, sizeof(float) * (size_t)(m + 1));
                Y.middleRows((Eigen::Index)q * n, n).noalias() = Ff.selfadjointView<Eigen::Upper>() * Vf;
                // g_q = Tr(B^q P) = sum_k s_k v_k . X^q_k, accumulated in FP64
                double gq = 0.0;
                for (int k = 0; k < r; ++k) {
                    const float* y = Y.data() + (size_t)k * Y.rows() + (size_t)q * n;
                    const double* v = Vd.data() + (size_t)k * n;
                    double dsum = 0.0;
                    for (int m = 0; m < n; ++m) dsum += v[m] * (double)y[m];
                    gq += k < npos ? dsum : -dsum;
                }
                double* jl = Jl.data();
                for (Eigen::Index p = 0; p < np; ++p) jl[p] += gq * (double)tq[p];
            }
            // K_b in FP32 (ssyrk), added into the FP64 K of this thread
            Kb.triangularView<Eigen::Lower>().setZero();
            if (npos > 0) {
                Eigen::Map<const Eigen::MatrixXf> Yp(Y.data(), n, (Eigen::Index)cnt * npos);
                Kb.selfadjointView<Eigen::Lower>().rankUpdate(Yp, 1.0f);
            }
            if (nneg > 0) {
                Eigen::Map<const Eigen::MatrixXf> Yn(Y.data() + (size_t)cnt * n * npos, n, (Eigen::Index)cnt * nneg);
                Kb.selfadjointView<Eigen::Lower>().rankUpdate(Yn, -1.0f);
            }
            for (int q = 0; q < n; ++q)
                for (int m = q; m < n; ++m) Kl(m, q) += (double)Kb(m, q);
        }
    }
    Eigen::MatrixXd Ks = Eigen::MatrixXd::Zero(n, n);
    std::vector<double> Js((size_t)np, 0.0);
    for (int t = 0; t < threads; ++t) {
        if (Kt[t].size())
            for (int q = 0; q < n; ++q)
                for (int m = q; m < n; ++m) Ks(m, q) += Kt[t](m, q);
        if (!Jt[t].empty())
            for (Eigen::Index p = 0; p < np; ++p) Js[p] += Jt[t][p];
    }
    J.resize(n, n);
    unpackSymmetric(Js.data(), n, J.data());
    K->resize(n, n);
    for (int q = 0; q < n; ++q)
        for (int m = q; m < n; ++m) {
            (*K)(m, q) = Ks(m, q);
            (*K)(q, m) = Ks(m, q);
        }
    ++m_fp_builds;
    m_fp_ms += msSince(t0);
}

// ---------------------------------------------------------------------------------
// 4a'. Local exchange fit (Claude Generated, Oct 2026; -qm.ri_fit_radius > 0, an APPROXIMATION)
// ---------------------------------------------------------------------------------
//
// P = L L^T with L the pivoted Cholesky factor of P (columns l_k localised around their pivot function for an
// insulator; exact factorisation, no orbital localisation iterations). K = sum_k (B l_k)(B l_k)^T is the global fit;
// the local fit replaces, per column k, the global auxiliary space by the domain F_k (auxiliary functions on atoms
// within R of the core atoms of l_k: the pivot atom plus the atoms with the largest sum_{m on A} l_mk^2 until they hold
// 99.9 % of sum_m l_mk^2):
//   Y_k(:, P) = T_n^P l_k (P in F_k),   K_k = Y_k J_F^-1 Y_k^T = X_k X_k^T,  X_k = Y_k L_F^-T,   J_F = L_F L_F^T.
// For one density the robust (Dunlap) and the symmetric local fit coincide (both sides use F_k), so the error per
// column is the Coulomb norm of the fitting residual of the densities chi_m l_k outside F_k; R -> infinity gives the
// global fit (tested). Domains follow the pivots, so the Fock matrix is only piecewise smooth in P (a pivot change
// changes a domain); J keeps the global fit (passes as in buildBatched). The metric factor of each distinct domain is
// cached for the geometry.

bool RIJK::buildBatchedLocal(const Matrix& P, Matrix& J, Matrix& K, const Eigen::MatrixXd* Vocc) const
{
    constexpr double kCoreFraction = 0.999;  // share of sum_m l_mk^2 on the core atoms of column k
    (void)Vocc;
    using clk = std::chrono::steady_clock;
    const Direct& d = *m_dd;
    const int n = m_n, N = m_naux, threads = std::max(1, m_opt.threads);
    const int naux_sh = (int)m_aux.size();
    const Eigen::Index np = npair();
    const int natoms = (int)m_atom_xyz.rows();
    curcuma::ScopedBlasThreads blas_guard(1);

    // 1. pivoted Cholesky of P (positive semidefinite for every density runSCF builds; an indefinite P falls back)
    const auto tf = clk::now();
    Eigen::MatrixXd Pd(P);
    Eigen::VectorXd dg = Pd.diagonal();
    const double dmax = dg.maxCoeff();
    if (!(dmax > 0.0)) return false;
    const double tol = 1.0e-12 * dmax;
    Eigen::MatrixXd Lc(n, n);
    std::vector<int> piv;
    int r = 0;
    const bool frozen = m_loc_frozen && (int)m_loc_piv.size() <= n;
    for (; r < n; ++r) {
        Eigen::Index p = 0;
        double dp = 0.0;
        if (frozen) {
            // fixed pivot order (Cholesky of P in a fixed permutation: smooth in P); stop where the old rank ended
            if (r >= (int)m_loc_piv.size()) break;
            p = m_loc_piv[r];
            dp = dg(p);
            if (dp <= tol) {  // rank changed: re-pivot from scratch
                m_loc_frozen = false;
                m_loc_unchanged = 0;
                return buildBatchedLocal(P, J, K, Vocc);
            }
        } else {
            dp = dg.maxCoeff(&p);
            if (dp <= tol) break;
        }
        Eigen::VectorXd col = Pd.col(p);
        if (r > 0) col.noalias() -= Lc.leftCols(r) * Lc.row(p).leftCols(r).transpose();
        col /= std::sqrt(dp);
        Lc.col(r) = col;
        dg -= col.cwiseAbs2();
        dg(p) = 0.0;
        piv.push_back((int)p);
    }
    if (dg.minCoeff() < -1.0e-8 * dmax) return false;  // indefinite
    if (frozen && dg.maxCoeff() > 1.0e-10 * dmax) {  // the frozen rank misses part of P: re-pivot
        m_loc_frozen = false;
        m_loc_unchanged = 0;
        return buildBatchedLocal(P, J, K, Vocc);
    }
    // 2. core atoms and domains (atom sets), grouped. Domains follow the pivots until two consecutive builds have the
    // same pivots and domains (or after 10 builds); from then on both stay fixed for this geometry, so the Fock matrix
    // is a smooth function of P (with moving domains the SCF did not converge at R = 4 Bohr for benzene).
    std::vector<std::vector<int>> dom(r);
    for (int k = 0; k < r && frozen; ++k) dom[k] = m_loc_dom[k];
    for (int k = 0; k < r && !frozen; ++k) {
        std::vector<double> w(natoms, 0.0);
        double tot = 0.0;
        for (int m = 0; m < n; ++m) {
            const double v = Lc(m, k) * Lc(m, k);
            w[m_ao_atom[m]] += v;
            tot += v;
        }
        // core: the pivot atom plus the atoms with the largest weights until they hold kCoreFraction of sum_m l_mk^2
        std::vector<int> order(natoms);
        for (int A = 0; A < natoms; ++A) order[A] = A;
        std::sort(order.begin(), order.end(), [&](int a, int b) { return w[a] > w[b] || (w[a] == w[b] && a < b); });
        std::vector<int> core{ m_ao_atom[piv[k]] };
        double acc = w[core[0]];
        for (int A : order) {
            if (acc >= kCoreFraction * tot) break;
            if (A == core[0]) continue;
            core.push_back(A);
            acc += w[A];
        }
        for (int B = 0; B < natoms; ++B)
            for (int A : core)
                if ((m_atom_xyz.row(A) - m_atom_xyz.row(B)).norm() <= m_opt.fit_radius) {
                    dom[k].push_back(B);
                    break;
                }
    }
    if (!frozen) {
        m_loc_unchanged = (piv == m_loc_piv && dom == m_loc_dom) ? m_loc_unchanged + 1 : 0;
        m_loc_piv = piv;
        m_loc_dom = dom;
        if (m_loc_unchanged >= 1 || m_loc_builds + 1 >= 10) m_loc_frozen = true;
    }
    // aux functions per atom (contiguous), domain function lists and offsets
    std::vector<int> afirst(natoms, -1), acount(natoms, 0);
    for (int Pi = 0; Pi < N; ++Pi) {
        const int A = m_auxfn_atom[Pi];
        if (afirst[A] < 0) afirst[A] = Pi;
        ++acount[A];
    }
    std::vector<int> fk(r, 0);
    std::vector<std::vector<int>> off(r);  // off[k][atom] = first local column of that atom in F_k, -1 if absent
    std::vector<std::shared_ptr<Eigen::MatrixXd>> LF(r);
    std::vector<std::vector<std::pair<int, int>>> byAtom(natoms);  // atom -> (orbital k, local offset)
    long new_dom = 0;
    for (int k = 0; k < r; ++k) {
        off[k].assign(natoms, -1);
        for (int B : dom[k]) {
            if (acount[B] == 0) continue;
            off[k][B] = fk[k];
            byAtom[B].push_back({ k, fk[k] });
            fk[k] += acount[B];
        }
        auto it = m_dom_cache.find(dom[k]);
        if (it == m_dom_cache.end()) {
            std::vector<int> fn;
            for (int B : dom[k])
                for (int j = 0; j < acount[B]; ++j) fn.push_back(afirst[B] + j);
            Eigen::MatrixXd JF((Eigen::Index)fn.size(), (Eigen::Index)fn.size());
            for (size_t a = 0; a < fn.size(); ++a)
                for (size_t b = 0; b < fn.size(); ++b) JF(a, b) = m_Jn(fn[a], fn[b]);
            Eigen::LLT<Eigen::MatrixXd> llt(JF);
            if (llt.info() != Eigen::Success) return false;
            it = m_dom_cache.emplace(dom[k], std::make_shared<Eigen::MatrixXd>(llt.matrixL())).first;
            ++new_dom;
        }
        LF[k] = it->second;
    }
    m_loc_factor_ms += msSince(tf);
    m_loc_new_domains += new_dom;
    // 3. passes: 0 = g and Y_k, 1 = J
    std::vector<Eigen::MatrixXd> Y(r);
    for (int k = 0; k < r; ++k) Y[k].resize(n, fk[k]);
    Eigen::VectorXd g = Eigen::VectorXd::Zero(N), c;
    std::vector<double> Jp((size_t)np, 0.0);
    Matrix Tb((Eigen::Index)m_batch_rows, np);
    std::vector<double> Pw;
    packContractionWeights(P, Pw);
    long blocks = 0, skipped = 0;
    for (int pass = 0; pass < 2; ++pass) {
        ++m_bt_passes;
        for (const auto& cb : m_aux_batches) {
            const int row0 = d.ao.first[cb.first];
            const int rows = (cb.second < naux_sh ? d.ao.first[cb.second] : N) - row0;
            const auto ti = clk::now();
            Tb.topRows(rows).setZero();
            fillRows(nullptr, cb.first, cb.second, Tb.data(), np, 0, blocks, skipped);
            for (int q = 0; q < rows; ++q) Tb.row(q) *= m_dn(row0 + q);
            m_bt_int_ms += msSince(ti);
            if (pass == 0) {
                const auto ty = clk::now();
#pragma omp parallel num_threads(threads)
                {
                    Eigen::MatrixXd Fu = Eigen::MatrixXd::Zero(n, n), Vb, Yb;
#pragma omp for schedule(dynamic)
                    for (int q = 0; q < rows; ++q) {
                        const int Pi = row0 + q, A = m_auxfn_atom[Pi];
                        const double* tq = Tb.data() + (size_t)q * np;
                        double s = 0.0;
                        for (Eigen::Index k = 0; k < np; ++k) s += tq[k] * Pw[k];
                        g(Pi) = s;
                        const auto& ks = byAtom[A];
                        if (ks.empty()) continue;
                        for (int m = 0; m < n; ++m)
                            std::memcpy(Fu.data() + (size_t)m * n, tq + (size_t)m * (m + 1) / 2, sizeof(double) * (size_t)(m + 1));
                        Vb.resize(n, (Eigen::Index)ks.size());
                        for (size_t j = 0; j < ks.size(); ++j) Vb.col(j) = Lc.col(ks[j].first);
                        Yb.noalias() = Fu.selfadjointView<Eigen::Upper>() * Vb;
                        const int jP = Pi - afirst[A];
                        for (size_t j = 0; j < ks.size(); ++j) Y[ks[j].first].col(ks[j].second + jP) = Yb.col(j);
                    }
                }
                m_loc_y_ms += msSince(ty);
            } else {
                const auto tj = clk::now();
                const Eigen::Index chunk = 1024, nchunks = (np + chunk - 1) / chunk;
#pragma omp parallel for schedule(static) num_threads(threads)
                for (Eigen::Index b = 0; b < nchunks; ++b) {
                    const Eigen::Index c0 = b * chunk, w = std::min(chunk, np - c0);
                    double* out = Jp.data() + c0;
                    for (int q = 0; q < rows; ++q) {
                        const double cq = c(row0 + q);
                        const double* t = Tb.data() + (size_t)q * np + c0;
                        for (Eigen::Index kk = 0; kk < w; ++kk) out[kk] += cq * t[kk];
                    }
                }
                m_bt_j_ms += msSince(tj);
            }
        }
        if (pass == 0) {
            c = metricSolve(g);
            // per column k: X_k = Y_k L_F^-T, K += X_k X_k^T
            const auto ts = clk::now();
            std::vector<Eigen::MatrixXd> Kt(threads, Eigen::MatrixXd());
#pragma omp parallel num_threads(threads)
            {
#ifdef _OPENMP
                const int tid = omp_get_thread_num();
#else
                const int tid = 0;
#endif
                Kt[tid] = Eigen::MatrixXd::Zero(n, n);
#pragma omp for schedule(dynamic)
                for (int k = 0; k < r; ++k) {
                    if (fk[k] == 0) continue;
                    LF[k]->transpose().triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheRight>(Y[k]);
                    Kt[tid].selfadjointView<Eigen::Lower>().rankUpdate(Y[k], 1.0);
                    Eigen::MatrixXd().swap(Y[k]);
                }
            }
            Eigen::MatrixXd Ks = Eigen::MatrixXd::Zero(n, n);
            for (int t = 0; t < threads; ++t)
                if (Kt[t].size())
                    for (int q = 0; q < n; ++q)
                        for (int m = q; m < n; ++m) Ks(m, q) += Kt[t](m, q);
            K.resize(n, n);
            for (int q = 0; q < n; ++q)
                for (int m = q; m < n; ++m) {
                    K(m, q) = Ks(m, q);
                    K(q, m) = Ks(m, q);
                }
            m_loc_solve_ms += msSince(ts);
        }
    }
    J.resize(n, n);
    unpackSymmetric(Jp.data(), n, J.data());
    ++m_loc_builds;
    m_loc_orbitals += r;
    for (int k = 0; k < r; ++k) m_loc_aux_sum += fk[k];
    return true;
}

// =================================================================================
// 4b. Integral-direct 3-index blocks (Claude Generated, Oct 2026, WP15 local fitting)
// =================================================================================
// The same helpers as RIJK::setup (makePair, unitShell, schwarzOf, pairBlock, transformAxis3, the OutMaps),
// called block by block; RIJK itself is not touched. Shell pairs are kept per first shell A in a sorted list of
// partners B >= A, so a lookup is a binary search.

struct ThreeIndexDirect::Impl {
    std::vector<Shell> osh, aux, units;
    OutMap oo, ao;
    std::vector<std::vector<int>> partner;  ///< partner[A] = sorted B >= A with a built pair
    std::vector<std::vector<int>> neigh;    ///< neigh[A] = sorted B (either order) with a built pair
    std::vector<std::vector<Pair>> pairs;   ///< pairs[A][k] = pair (A, partner[A][k])
    std::vector<Pair> apairs;
    double screening = 1e-12;
    long npairs = 0;
    const Pair* find(int A, int B) const
    {
        const std::vector<int>& p = partner[A];
        auto it = std::lower_bound(p.begin(), p.end(), B);
        if (it == p.end() || *it != B) return nullptr;
        return &pairs[A][it - p.begin()];
    }
};

ThreeIndexDirect::ThreeIndexDirect() = default;
ThreeIndexDirect::~ThreeIndexDirect() = default;

bool ThreeIndexDirect::setup(const std::vector<GTO::Orbital>& basis, const Matrix* Q, const std::vector<Shell>& aux,
                             bool spherical_aux, double screening, int threads, std::string& why)
{
    auto d = std::make_shared<Impl>();
    d->osh = orbitalShells(basis);
    d->oo = orbitalOutMap(d->osh, (int)basis.size(), Q);
    d->aux = aux;
    d->ao = auxOutMap(d->aux, spherical_aux);
    d->screening = screening;
    const int ns = (int)d->osh.size(), na = (int)d->aux.size();
    if (ns == 0 || na == 0) {
        why = "empty orbital or auxiliary basis";
        return false;
    }
    d->units.reserve(na);  // apairs point into units: no reallocation after this
    for (const Shell& s : d->aux) d->units.push_back(unitShell(s));
    d->apairs.reserve(na);
    for (int c = 0; c < na; ++c) d->apairs.push_back(makePair(d->aux[c], d->units[c]));
    for (Pair& p : d->apairs) p.schwarz = schwarzOf(p);
    d->partner.assign(ns, {});
    d->pairs.resize(ns);
    constexpr double kPrefactorCut = 1.0e-22;  // exp(-ab/(a+b) R^2) below this for every primitive pair: not built
#pragma omp parallel for schedule(dynamic) num_threads(std::max(1, threads))
    for (int A = 0; A < ns; ++A) {
        const Shell& sA = d->osh[A];
        for (int B = A; B < ns; ++B) {
            const Shell& sB = d->osh[B];
            const double R2 = (sA.x - sB.x) * (sA.x - sB.x) + (sA.y - sB.y) * (sA.y - sB.y) + (sA.z - sB.z) * (sA.z - sB.z);
            double kmax = 0.0;
            for (double a : sA.exps)
                for (double b : sB.exps) kmax = std::max(kmax, std::exp(-a * b / (a + b) * R2));
            if (kmax < kPrefactorCut) continue;
            d->partner[A].push_back(B);
            d->pairs[A].push_back(makePair(sA, sB));
            d->pairs[A].back().schwarz = schwarzOf(d->pairs[A].back());
        }
    }
    d->neigh.assign(ns, {});
    for (int A = 0; A < ns; ++A) {
        d->npairs += (long)d->partner[A].size();
        for (int B : d->partner[A]) {
            d->neigh[A].push_back(B);
            if (B != A) d->neigh[B].push_back(A);
        }
    }
    for (auto& v : d->neigh) std::sort(v.begin(), v.end());
    m_impl = d;
    return true;
}

int ThreeIndexDirect::n() const { return m_impl->oo.n; }
int ThreeIndexDirect::naux() const { return m_impl->ao.n; }
int ThreeIndexDirect::nShells() const { return (int)m_impl->osh.size(); }
int ThreeIndexDirect::shellFirst(int A) const { return m_impl->oo.first[A]; }
int ThreeIndexDirect::shellWidth(int A) const { return m_impl->oo.width[A]; }
int ThreeIndexDirect::shellAtom(int A) const { return m_impl->osh[A].atom; }
int ThreeIndexDirect::nAuxShells() const { return (int)m_impl->aux.size(); }
int ThreeIndexDirect::auxFirst(int c) const { return m_impl->ao.first[c]; }
int ThreeIndexDirect::auxWidth(int c) const { return m_impl->ao.width[c]; }
int ThreeIndexDirect::auxAtom(int c) const { return m_impl->aux[c].atom; }
double ThreeIndexDirect::auxSchwarz(int c) const { return m_impl->apairs[c].schwarz; }
long ThreeIndexDirect::nPairs() const { return m_impl->npairs; }
const std::vector<int>& ThreeIndexDirect::neighbours(int A) const { return m_impl->neigh[A]; }

double ThreeIndexDirect::pairSchwarz(int A, int B) const
{
    const Pair* p = m_impl->find(std::min(A, B), std::max(A, B));
    return p ? p->schwarz : 0.0;
}

bool ThreeIndexDirect::block(int c, int A, int B, std::vector<double>& out) const
{
    const Impl& d = *m_impl;
    const int A0 = std::min(A, B), B0 = std::max(A, B);
    const Pair* pr = d.find(A0, B0);
    if (!pr || pr->schwarz * d.apairs[c].schwarz < d.screening) return false;
    thread_local std::vector<double> buf, tmp;
    pairBlock(*pr, d.apairs[c], buf);
    int dims[3] = { (int)d.osh[A0].lmn.size(), (int)d.osh[B0].lmn.size(), (int)d.aux[c].lmn.size() };
    if (d.oo.T[A0].size()) transformAxis3(buf, tmp, dims, 0, d.oo.T[A0]);
    if (d.oo.T[B0].size()) transformAxis3(buf, tmp, dims, 1, d.oo.T[B0]);
    if (d.ao.T[c].size()) transformAxis3(buf, tmp, dims, 2, d.ao.T[c]);
    if (A == A0) {
        out.swap(buf);
        return true;
    }
    // requested (A, B) with A > B: out[(a * wB + b) * wc + p] = buf[(b * wA + a) * wc + p] with wA = dims[1]
    const int wB = dims[0], wA = dims[1], wc = dims[2];
    out.resize((size_t)wA * wB * wc);
    for (int a = 0; a < wA; ++a)
        for (int b = 0; b < wB; ++b)
            std::memcpy(&out[((size_t)a * wB + b) * wc], &buf[((size_t)b * wA + a) * wc], sizeof(double) * wc);
    return true;
}

// =================================================================================
// 6. Analytic nuclear gradient of the RI energy (Claude Generated, Oct 2026)
// =================================================================================
//
// Closed-shell RHF, spin-summed density P = 2 C_o C_o^T. Raw auxiliary basis {P}, Coulomb metric
// J_PQ = (P|Q), 3-index integrals (mn|P), g_P = sum_mn (P|mn) P_mn.
//
//   E_J = 1/2 g^T J^-1 g                                   (RI-J and RI-JK)
//   E_K = -1/4 sum_mnls P_mn P_ls sum_PQ (ml|P) [J^-1]_PQ (Q|ns)   (RI-JK only)
//
// The RI-SCF is variational for its own energy (F = H + J_RI - 1/2 K_RI with J_RI = dE_J/dP,
// K_RI = -2 dE_K/dP), so only the explicit dependence on the integrals is differentiated; the
// orbital (orthonormality) response is the -W dS term of the one-electron part, with W built from
// the RI-SCF orbital energies. With d(J^-1) = -J^-1 (dJ) J^-1:
//
//   dE_J/dx = sum_P c_P sum_mn (mn|P)^x P_mn - 1/2 sum_PQ c_P c_Q (P|Q)^x,       c = J^-1 g
//   dE_K/dx = -2 sum_P sum_mn (mn|P)^x Gam^P_mn + sum_PQ M_PQ (P|Q)^x
//     Gam^P = 1/4 sum_Q [J^-1]_PQ (P (Q|..) P),   M = 1/4 J^-1 A J^-1,  A_RS = Tr((R|..) P (S|..) P)
//
// (P = 2 C_o C_o^T turns Gam^P_mn = sum_ij C_mi C_nj d^P_ij, d^P_ij = sum_Q [J^-1]_PQ (Q|ij) and
// M_PQ = sum_ij d^P_ij d^Q_ij into AO quantities.) Weigend, Haeser, Patzelt, Ahlrichs, Chem. Phys.
// Lett. 294, 143 (1998) (RI-J gradient); Weigend, Phys. Chem. Chem. Phys. 4, 4285 (2002) (RI-JK);
// the same expressions as PySCF df/grad/rhf.py (auxbasis_response). Together:
//
//   dE/dx = sum_P sum_mn G^P_mn (mn|P)^x - sum_PQ M2_PQ (P|Q)^x,
//   G^P = c_P P - 2 Gam^P,   M2 = 1/2 c c^T - M      (mode j: G^P = c_P P, M2 = 1/2 c c^T).
//
// With the stored fitted tensor B = F T_n (T_n = D T_raw, D_PP = 1/sqrt((P|P)), F = L^-1 for the
// Cholesky metric J_n = L L^T, F = w^-1/2 U^T for the eigen metric), J_raw^-1 = D F^T F D:
//   c = D F^T (B . P),  Gam = 1/4 D F^T (P B P),  M = 1/4 D F^T Z F D,  Z_RS = Tr(B^R P B^S P),
// formed with the density factor P = V S V^T (Y^Q = V (S W^Q S) V^T, W^Q = V^T B^Q V, Z = W_s W^T).
// With the eigen metric and dropped eigenvalues J^-1 is a pseudo-inverse and the formula is not the
// exact derivative (not used by default; nothing is dropped at the default cutoff in the tests).
//
// Derivative integrals: d/dA_k of a cartesian Gaussian = 2 alpha (l+1_k) - l_k (l-1_k) (same
// normalisation), from Hermite tables raised by one on both orbital shells (as gradientTwoElectron in
// qm_integrals.cpp). (mn|P): derivatives on the centres of m and n, the one on the auxiliary centre
// from translational invariance, d/dC = -(d/dA + d/dB). (P|Q): derivative on P's centre, Q's centre
// gets the negative. Only first derivatives, so the raised orbital tables need L+1 (d -> f on the
// derivative side) and the raised auxiliary tables L+1 (g -> h) only in the 2-index part.

namespace {

/// As makePair, Hermite tables raised by one on BOTH shells (centre derivatives of either shell).
Pair makePairRaised(const Shell& A, const Shell& B)
{
    Pair pr;
    pr.A = &A;
    pr.B = &B;
    for (size_t ia = 0; ia < A.exps.size(); ++ia)
        for (size_t ib = 0; ib < B.exps.size(); ++ib) {
            const double a = A.exps[ia], b = B.exps[ib];
            PrimPair pp;
            pp.ia = (int)ia;
            pp.ib = (int)ib;
            pp.p = a + b;
            pp.Px = (a * A.x + b * B.x) / pp.p;
            pp.Py = (a * A.y + b * B.y) / pp.p;
            pp.Pz = (a * A.z + b * B.z) / pp.p;
            const double ab2 = (A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) + (A.z - B.z) * (A.z - B.z);
            const double K = std::exp(-a * b / pp.p * ab2);
            if (K == 0.0) continue;
            pp.Ex = hermite(A.L + 1, B.L + 1, pp.Px - A.x, pp.Px - B.x, pp.p, K);
            pp.Ey = hermite(A.L + 1, B.L + 1, pp.Py - A.y, pp.Py - B.y, pp.p, 1.0);
            pp.Ez = hermite(A.L + 1, B.L + 1, pp.Pz - A.z, pp.Pz - B.z, pp.p, 1.0);
            pr.prims.push_back(std::move(pp));
        }
    return pr;
}

/// Derivative integrals of one bra pair against one ket pair, contracted with cartesian weights:
///   gA_k = sum_{a,b,cd} W[(a nB + b) nCD + cd] d/dA_k (ab|cd),  gB_k the same for the centre of B.
/// bra: makePairRaised; ket: makePair (not differentiated here). The ket Hermite sum G_cd(t,u,v)
/// is the one of pairBlock (bra Hermite range raised by one); per bra component pair it is first
/// contracted with the weights, H = sum_cd W_cd G_cd, so each derivative term is one E E E H sum.
void derivBlock(const Pair& bra, const Pair& ket, const double* W, bool wantB, double gA[3], double gB[3])
{
    const Shell& A = *bra.A;
    const Shell& B = *bra.B;
    const Shell& C = *ket.A;
    const Shell& D = *ket.B;
    const int nA = (int)A.lmn.size(), nB = (int)B.lmn.size(), nC = (int)C.lmn.size(), nD = (int)D.lmn.size();
    const int nCD = nC * nD;
    for (int k = 0; k < 3; ++k) gA[k] = gB[k] = 0.0;
    const int LAB = A.L + B.L + 1, LCD = C.L + D.L, Ltot = LAB + LCD;
    const int sAB_t = (A.L + 2) * (B.L + 2), sAB_i = B.L + 2;  // raised tables
    const int sCD_t = (C.L + 1) * (D.L + 1), sCD_i = D.L + 1;
    const int s = Ltot + 1;
    auto at = [s](int t, int u, int v) { return ((size_t)t * s + u) * s + v; };
    const int g = LAB + 1;
    const size_t gsz = (size_t)g * g * g;
    thread_local std::vector<double> F, R, lvl, G, H;
    G.resize((size_t)nCD * gsz);
    H.resize(gsz);

    for (const PrimPair& pb : bra.prims) {
        const double alpha = A.exps[pb.ia], beta = B.exps[pb.ib];
        std::fill(G.begin(), G.end(), 0.0);
        for (const PrimPair& pk : ket.prims) {
            const double p = pb.p, q = pk.p, rho = p * q / (p + q);
            const double Wx = pb.Px - pk.Px, Wy = pb.Py - pk.Py, Wz = pb.Pz - pk.Pz;
            qmint::boysFunction(Ltot, rho * (Wx * Wx + Wy * Wy + Wz * Wz), F);
            coulombR0(Ltot, Wx, Wy, Wz, rho, F, R, lvl);
            const double pref = 2.0 * std::pow(PI, 2.5) / (p * q * std::sqrt(p + q));
            for (int cc = 0; cc < nC; ++cc) {
                const auto& lc = C.lmn[cc];
                const double wc = pref * C.coef[cc][pk.ia];
                for (int cd = 0; cd < nD; ++cd) {
                    const double w = wc * D.coef[cd][pk.ib];
                    if (w == 0.0) continue;
                    const auto& ld = D.lmn[cd];
                    const int tx = lc[0] + ld[0], ty = lc[1] + ld[1], tz = lc[2] + ld[2];
                    double ex[16], ey[16], ez[16];
                    for (int k = 0; k <= tx; ++k) ex[k] = ((k & 1) ? -1.0 : 1.0) * pk.Ex[k * sCD_t + lc[0] * sCD_i + ld[0]];
                    for (int k = 0; k <= ty; ++k) ey[k] = ((k & 1) ? -1.0 : 1.0) * pk.Ey[k * sCD_t + lc[1] * sCD_i + ld[1]];
                    for (int k = 0; k <= tz; ++k) ez[k] = ((k & 1) ? -1.0 : 1.0) * pk.Ez[k * sCD_t + lc[2] * sCD_i + ld[2]];
                    double* Gcd = &G[(size_t)(cc * nD + cd) * gsz];
                    for (int t = 0; t <= LAB; ++t)
                        for (int u = 0; u <= LAB - t; ++u)
                            for (int v = 0; v <= LAB - t - u; ++v) {
                                double acc = 0.0;
                                for (int k1 = 0; k1 <= tx; ++k1) {
                                    if (ex[k1] == 0.0) continue;
                                    double a2 = 0.0;
                                    for (int k2 = 0; k2 <= ty; ++k2) {
                                        if (ey[k2] == 0.0) continue;
                                        const double* Rr = &R[at(t + k1, u + k2, v)];
                                        double a3 = 0.0;
                                        for (int k3 = 0; k3 <= tz; ++k3) a3 += ez[k3] * Rr[k3];
                                        a2 += ey[k2] * a3;
                                    }
                                    acc += ex[k1] * a2;
                                }
                                Gcd[((size_t)t * g + u) * g + v] += w * acc;
                            }
                }
            }
        }
        // (a'b'|H) for bra powers a, b (either one shifted by +-1 on one axis)
        auto val = [&](const int a[3], const int b[3]) {
            double sum = 0.0;
            for (int t = 0; t <= a[0] + b[0]; ++t) {
                const double et = pb.Ex[t * sAB_t + a[0] * sAB_i + b[0]];
                if (et == 0.0) continue;
                for (int u = 0; u <= a[1] + b[1]; ++u) {
                    const double etu = et * pb.Ey[u * sAB_t + a[1] * sAB_i + b[1]];
                    if (etu == 0.0) continue;
                    const double* hrow = &H[((size_t)t * g + u) * g];
                    double sv = 0.0;
                    for (int v = 0; v <= a[2] + b[2]; ++v) sv += pb.Ez[v * sAB_t + a[2] * sAB_i + b[2]] * hrow[v];
                    sum += etu * sv;
                }
            }
            return sum;
        };
        for (int ca = 0; ca < nA; ++ca) {
            const auto& la = A.lmn[ca];
            for (int cb = 0; cb < nB; ++cb) {
                const double c_ab = A.coef[ca][pb.ia] * B.coef[cb][pb.ib];
                if (c_ab == 0.0) continue;
                const double* w = W + (size_t)(ca * nB + cb) * nCD;
                bool any = false;
                for (int cd = 0; cd < nCD; ++cd) {
                    if (w[cd] == 0.0) continue;
                    const double* Gk = &G[(size_t)cd * gsz];
                    if (!any) {
                        for (size_t x = 0; x < gsz; ++x) H[x] = w[cd] * Gk[x];
                        any = true;
                    } else {
                        for (size_t x = 0; x < gsz; ++x) H[x] += w[cd] * Gk[x];
                    }
                }
                if (!any) continue;
                const auto& lb = B.lmn[cb];
                int a[3] = { la[0], la[1], la[2] }, b[3] = { lb[0], lb[1], lb[2] };
                for (int k = 0; k < 3; ++k) {
                    a[k] += 1;
                    double dA = 2.0 * alpha * val(a, b);
                    a[k] -= 2;
                    if (la[k] > 0) dA -= la[k] * val(a, b);
                    a[k] += 1;
                    gA[k] += c_ab * dA;
                    if (!wantB) continue;
                    b[k] += 1;
                    double dB = 2.0 * beta * val(a, b);
                    b[k] -= 2;
                    if (lb[k] > 0) dB -= lb[k] * val(a, b);
                    b[k] += 1;
                    gB[k] += c_ab * dB;
                }
            }
        }
    }
}

/// In place T := L^-T T (L lower triangular N x N)  <=>  X := X L^-1 per column chunk (dtrsm, see lowerSolveRows).
void upperSolveRowsLT(const Eigen::MatrixXd& L, Matrix& T, int threads)
{
    const Eigen::Index N = T.rows(), ncol = T.cols(), nchunks = (ncol + kRowChunk - 1) / kRowChunk;
    double* base = T.data();
#pragma omp parallel for schedule(dynamic) num_threads(std::max(1, threads))
    for (Eigen::Index b = 0; b < nchunks; ++b) {
        const Eigen::Index c0 = b * kRowChunk, w = std::min(kRowChunk, ncol - c0);
        StridedMap X(base + c0, w, N, Eigen::OuterStride<>(ncol));
        L.triangularView<Eigen::Lower>().solveInPlace<Eigen::OnTheRight>(X);
    }
}

}  // namespace

void RIJK::applyFitTranspose(Matrix& X) const
{
    const int threads = std::max(1, m_opt.threads);
    if (m_gL.size()) {
        upperSolveRowsLT(m_gL, X, threads);
        return;
    }
    if (m_gU.size()) {
        const Eigen::MatrixXd F2 = m_gU * m_gw.cwiseSqrt().cwiseInverse().asDiagonal();  // N x nfit
        Matrix out;
        leftMultiplyRows(F2, X, out, threads);
        X.swap(out);
        return;
    }
    throw std::logic_error("RIJK::applyFitTranspose: no metric factor (setup not run)");
}

bool RIJK::gradientDensities(const Matrix& P, double kx, Eigen::VectorXd& c, Matrix& G, Eigen::MatrixXd& M2, std::string& why) const
{
    const int n = m_n, N = m_naux, threads = std::max(1, m_opt.threads);
    const int rows = (int)m_T.rows();  // mode jk: nfit (B), mode j: naux (T_n)
    const Matrix Pc = P;               // contiguous row-major
    const Eigen::Index ncol = npair(); // packed tensor and G columns
    std::vector<double> Pw;
    packContractionWeights(Pc, Pw);
    // g_Q = sum_mn T_Q,mn P_mn = sum_p T_Q,p w_p
    Eigen::VectorXd gq(rows);
#pragma omp parallel for schedule(static) num_threads(threads)
    for (int Q = 0; Q < rows; ++Q) {
        const double* t = m_T.data() + (size_t)Q * ncol;
        double s = 0.0;
        for (Eigen::Index k = 0; k < ncol; ++k) s += t[k] * Pw[k];
        gq(Q) = s;
    }
    Eigen::VectorXd cn;
    if (m_opt.mode == "j") {
        cn = fitCoefficients(gq);  // J_n^-1 g_n
    } else {
        Matrix x(rows, 1);
        for (int Q = 0; Q < rows; ++Q) x(Q, 0) = gq(Q);
        applyFitTranspose(x);  // F^T (B . P)
        cn = x.col(0);
    }
    c = m_dn.cwiseProduct(cn);  // raw auxiliary basis
    M2 = 0.5 * c * c.transpose();
    G.resize(0, 0);
    if (m_opt.mode == "j" || kx == 0.0) return true;  // Coulomb only (RI-J, or a pure density functional)

    const double need_mb = (double)N * ncol * 8.0 / 1.0e6;
    if (need_mb + tensorMB() > m_opt.max_memory_mb) {
        why = fmt::format("RI-JK gradient needs {:.1f} MB for the 3-index densities in addition to the {:.1f} MB tensor, above "
                          "-qm.ri_max_memory_mb {:g}",
                          need_mb, tensorMB(), m_opt.max_memory_mb);
        return false;
    }
    Eigen::MatrixXd V;
    int npos = 0;
    factorDensity(P, V, npos);
    const int r = (int)V.cols(), r2 = r * r;
    std::vector<double> Vr((size_t)n * r);  // row-major n x r
    for (int l = 0; l < n; ++l)
        for (int k = 0; k < r; ++k) Vr[(size_t)l * r + k] = V(l, k);
    // Per fit row Q: X = B^Q V (n x r), W^Q = V^T X (r x r), Ws^Q = (s s^T) o W^Q, so that
    // P B^Q P = V Ws^Q V^T. F^T acts on the auxiliary index only, hence F^T (P B P) = V (F^T Ws) V^T:
    // the metric solve runs on the r^2 columns of Ws instead of the n^2 columns of P B P.
    std::vector<double> Wr((size_t)rows * r2);
    Matrix Ws(rows, r2);
#pragma omp parallel num_threads(threads)
    {
        std::vector<double> X((size_t)n * r), Bfull((size_t)n * n);
#pragma omp for schedule(static)
        for (int Q = 0; Q < rows; ++Q) {
            unpackRow(Q, Bfull.data());
            const double* B = Bfull.data();
            std::fill(X.begin(), X.end(), 0.0);
            for (int m = 0; m < n; ++m) {
                double* xm = &X[(size_t)m * r];
                const double* bm = B + (size_t)m * n;
                for (int l = 0; l < n; ++l) {
                    const double b = bm[l];
                    if (b == 0.0) continue;
                    const double* vl = &Vr[(size_t)l * r];
                    for (int k = 0; k < r; ++k) xm[k] += b * vl[k];
                }
            }
            double* W = &Wr[(size_t)Q * r2];
            std::fill(W, W + r2, 0.0);
            for (int m = 0; m < n; ++m) {
                const double* vm = &Vr[(size_t)m * r];
                const double* xm = &X[(size_t)m * r];
                for (int k = 0; k < r; ++k) {
                    const double v = vm[k];
                    if (v == 0.0) continue;
                    double* wk = W + (size_t)k * r;
                    for (int l = 0; l < r; ++l) wk[l] += v * xm[l];
                }
            }
            double* ws = Ws.data() + (size_t)Q * r2;
            for (int k = 0; k < r; ++k)
                for (int l = 0; l < r; ++l) ws[(size_t)k * r + l] = ((k < npos) == (l < npos)) ? W[(size_t)k * r + l] : -W[(size_t)k * r + l];
        }
    }
    // Z_RS = Tr(B^R P B^S P) = sum_kl Ws^R_kl W^S_kl
    Matrix Z(rows, rows);
#pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int R = 0; R < rows; ++R) {
        const double* a = Ws.data() + (size_t)R * r2;
        for (int S = 0; S <= R; ++S) {
            const double* b = &Wr[(size_t)S * r2];
            double sum = 0.0;
            for (int k = 0; k < r2; ++k) sum += a[k] * b[k];
            Z(R, S) = sum;
            Z(S, R) = sum;
        }
    }
    std::vector<double>().swap(Wr);
    // Gam = 1/4 D F^T (P B P) = 1/4 D V (F^T Ws) V^T  ->  G^P = c_P P - 1/2 D_PP V (F^T Ws)^P V^T
    applyFitTranspose(Ws);  // N rows
    G.resize(N, ncol);
#pragma omp parallel num_threads(threads)
    {
        std::vector<double> U((size_t)n * r);
#pragma omp for schedule(static)
        for (int Pi = 0; Pi < N; ++Pi) {
            const double* ws = Ws.data() + (size_t)Pi * r2;
            // U = V Ws (n x r)
            std::fill(U.begin(), U.end(), 0.0);
            for (int m = 0; m < n; ++m) {
                const double* vm = &Vr[(size_t)m * r];
                double* um = &U[(size_t)m * r];
                for (int k = 0; k < r; ++k) {
                    const double v = vm[k];
                    if (v == 0.0) continue;
                    const double* wk = ws + (size_t)k * r;
                    for (int l = 0; l < r; ++l) um[l] += v * wk[l];
                }
            }
            // G^P = c_P P - 1/2 D_PP U V^T (symmetric: packed lower triangle; P symmetrised, which changes
            // nothing for an exactly symmetric density)
            double* y = G.data() + (size_t)Pi * ncol;
            const double cp = c(Pi), dp = -0.5 * kx * m_dn(Pi);  // kx: WP8 hybrid (1 for RHF)
            for (int m = 0; m < n; ++m) {
                const double* um = &U[(size_t)m * r];
                double* ym = y + (size_t)m * (m + 1) / 2;
                for (int l = 0; l <= m; ++l) {
                    const double* vl = &Vr[(size_t)l * r];
                    double sum = 0.0;
                    for (int k = 0; k < r; ++k) sum += um[k] * vl[k];
                    ym[l] = cp * (l == m ? Pc(m, m) : 0.5 * (Pc(m, l) + Pc(l, m))) + dp * sum;
                }
            }
        }
    }
    Matrix().swap(Ws);
    // M = 1/4 D F^T Z F D  (F^T Z, transpose, F^T again; Z symmetric)
    applyFitTranspose(Z);  // N x rows
    Matrix Zt = Z.transpose();
    applyFitTranspose(Zt);  // N x N
    for (int Pi = 0; Pi < N; ++Pi)
        for (int Qi = 0; Qi < N; ++Qi) M2(Pi, Qi) -= kx * 0.25 * m_dn(Pi) * m_dn(Qi) * 0.5 * (Zt(Pi, Qi) + Zt(Qi, Pi));
    return true;
}

namespace {

/// Shell data shared by the host and the device 3-index gradient (built once per gradient).
struct GradSetup {
    std::vector<Shell> osh;
    OutMap oo, ao;
    std::vector<Eigen::MatrixXd> oT, aT;  ///< transposed transforms (active/spherical -> cartesian)
    std::vector<Shell> units;
    std::vector<Pair> apairs;              ///< (C, unit s) per auxiliary shell
    std::vector<std::pair<int, int>> oidx; ///< canonical orbital shell pairs, order of RIJK::setup
    GradSetup();
    ~GradSetup();
};
GradSetup::GradSetup() = default;
GradSetup::~GradSetup() = default;

void buildGradSetup(const std::vector<GTO::Orbital>& basis, const Matrix* Q, const std::vector<Shell>& aux, bool spherical,
                    GradSetup& gs)
{
    gs.osh = orbitalShells(basis);
    gs.oo = orbitalOutMap(gs.osh, (int)basis.size(), Q);
    gs.ao = auxOutMap(aux, spherical);
    const int nsh = (int)gs.osh.size(), naux_sh = (int)aux.size();
    gs.oT.assign(nsh, Eigen::MatrixXd());
    gs.aT.assign(naux_sh, Eigen::MatrixXd());
    for (int s = 0; s < nsh; ++s)
        if (gs.oo.T[s].size()) gs.oT[s] = gs.oo.T[s].transpose();
    for (int s = 0; s < naux_sh; ++s)
        if (gs.ao.T[s].size()) gs.aT[s] = gs.ao.T[s].transpose();
    gs.units.clear();
    gs.units.reserve(aux.size());
    for (const Shell& s : aux) gs.units.push_back(unitShell(s));
    gs.apairs.clear();
    gs.apairs.reserve(aux.size());
    for (int s = 0; s < naux_sh; ++s) gs.apairs.push_back(makePair(aux[s], gs.units[s]));
    gs.oidx.clear();
    for (int A = 0; A < nsh; ++A)
        for (int B = A; B < nsh; ++B) gs.oidx.push_back({ A, B });
}

/// Cartesian weight block W[a][b][c] of the 3-index density for orbital pair k and auxiliary shell cs
/// (G^P_mn, or c_P P_mn when G is empty), back-transformed from the active/spherical functions. False
/// when the block is skipped: one-centre block, or below the density-weighted Schwarz screen.
bool weightBlock(const GradSetup& gs, const std::vector<Shell>& aux, int k, int cs, int n, const Matrix& G, const Eigen::VectorXd& c,
                 const Matrix& Pc, double schwarz, double thr, std::vector<double>& blk, std::vector<double>& tmp)
{
    const int A = gs.oidx[k].first, B = gs.oidx[k].second;
    if (gs.osh[A].atom == gs.osh[B].atom && gs.osh[B].atom == aux[cs].atom) return false;  // no gradient
    const int wA = gs.oo.width[A], wB = gs.oo.width[B], wC = gs.ao.width[cs];
    const bool fullG = G.size() != 0;
    blk.resize((size_t)wA * wB * wC);
    double gmax = 0.0;
    for (int i = 0; i < wA; ++i)
        for (int j = 0; j < wB; ++j) {
            const int m = gs.oo.first[A] + i, l = gs.oo.first[B] + j;
            const Eigen::Index col = fullG ? pairIndex(m, l) : (Eigen::Index)m * n + l;  // G packed, P full
            double* o = &blk[((size_t)i * wB + j) * wC];
            for (int s = 0; s < wC; ++s) {
                const int Pi = gs.ao.first[cs] + s;
                const double v = fullG ? G(Pi, col) : c(Pi) * Pc.data()[col];
                o[s] = v;
                gmax = std::max(gmax, std::abs(v));
            }
        }
    if (schwarz >= 0.0 && schwarz * gmax < thr) return false;
    int d[3] = { wA, wB, wC };
    if (gs.oT[A].size()) transformAxis3(blk, tmp, d, 0, gs.oT[A]);
    if (gs.oT[B].size()) transformAxis3(blk, tmp, d, 1, gs.oT[B]);
    if (gs.aT[cs].size()) transformAxis3(blk, tmp, d, 2, gs.aT[cs]);
    return true;
}

}  // namespace

Matrix RIJK::threeIndexGradient(const std::vector<GTO::Orbital>& basis, const Matrix* Q, const Eigen::VectorXd& c, const Matrix& P,
                                const Matrix& G, int natoms) const
{
    const int threads = std::max(1, m_opt.threads);
    GradSetup gs;
    buildGradSetup(basis, Q, m_aux, m_opt.spherical_aux, gs);
    const int n = m_n, naux_sh = (int)m_aux.size(), npairs = (int)gs.oidx.size();
    // The derivative integrals are not bounded by the Schwarz factors of the undifferentiated ones
    // (factors 2 alpha), so the density-weighted screen runs two orders below the RI screening.
    const double thr = m_opt.screening * 1.0e-2;
    const bool haveSchwarz = (int)m_osch.size() == npairs && (int)m_asch.size() == naux_sh;
    const Matrix Pc = P;
    // Oct 2026 (WP8): fixed round-robin chunks with own accumulators, added in order (thread-count independent)
    const int nchunk = std::max(1, std::min(256, npairs));
    std::vector<Matrix> gt(nchunk, Matrix::Zero(natoms, 3));
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (int ch = 0; ch < nchunk; ++ch) {
        Matrix& g = gt[ch];
        std::vector<double> blk, tmp;
        for (int k = ch; k < npairs; k += nchunk) {
            const int A = gs.oidx[k].first, B = gs.oidx[k].second;
            const int aA = gs.osh[A].atom, aB = gs.osh[B].atom;
            Pair braR;
            bool built = false;
            for (int cs = 0; cs < naux_sh; ++cs) {
                const double sch = haveSchwarz ? m_osch[k] * m_asch[cs] : -1.0;
                if (!weightBlock(gs, m_aux, k, cs, n, G, c, Pc, sch, thr, blk, tmp)) continue;
                if (!built) {
                    braR = makePairRaised(gs.osh[A], gs.osh[B]);
                    built = true;
                }
                double gA[3], gB[3];
                derivBlock(braR, gs.apairs[cs], blk.data(), true, gA, gB);
                const int aC = m_aux[cs].atom;
                const double w = (A != B) ? 2.0 : 1.0;  // (A,B) and (B,A) blocks, symmetric G
                for (int x = 0; x < 3; ++x) {
                    g(aA, x) += w * gA[x];
                    g(aB, x) += w * gB[x];
                    g(aC, x) -= w * (gA[x] + gB[x]);
                }
            }
        }
    }
    Matrix grad = Matrix::Zero(natoms, 3);
    for (int t = 0; t < nchunk; ++t) grad += gt[t];
    return grad;
}

// Claude Generated (Oct 2026): the same quantity as threeIndexGradient with the integrals on the device
// (cuda/qm_gpu_ri_grad.cu). The host builds the shells, the Hermite tables of every primitive pair
// (bra raised, kets plain), the screened task list and the cartesian weight blocks (weightBlock).
bool RIJK::threeIndexGradientDevice(const std::vector<GTO::Orbital>& basis, const Matrix* Q, const Eigen::VectorXd& c, const Matrix& P,
                                    const Matrix& G, int natoms, const std::string& backend, int device, Matrix& grad,
                                    std::string& info, std::string& why) const
{
    const auto t0 = std::chrono::steady_clock::now();
    const int threads = std::max(1, m_opt.threads);
    GradSetup gs;
    buildGradSetup(basis, Q, m_aux, m_opt.spherical_aux, gs);
    const int n = m_n, nsh = (int)gs.osh.size(), naux_sh = (int)m_aux.size(), npairs = (int)gs.oidx.size();
    const double thr = m_opt.screening * 1.0e-2;
    const bool haveSchwarz = (int)m_osch.size() == npairs && (int)m_asch.size() == naux_sh;
    const Matrix Pc = P;

    // shells: orbital then auxiliary
    std::vector<int> sh_L, sh_atom, sh_ncomp, sh_lmn_off, lmn, sh_nprim, sh_coef_off, sh_exp_off;
    std::vector<double> coef, exps;
    auto addShell = [&](const Shell& s) {
        sh_L.push_back(s.L);
        sh_atom.push_back(s.atom);
        sh_ncomp.push_back((int)s.lmn.size());
        sh_lmn_off.push_back((int)lmn.size());
        for (const auto& t : s.lmn) lmn.insert(lmn.end(), t.begin(), t.end());
        sh_nprim.push_back((int)s.exps.size());
        sh_coef_off.push_back((int)coef.size());
        for (const auto& cc : s.coef) coef.insert(coef.end(), cc.begin(), cc.end());
        sh_exp_off.push_back((int)exps.size());
        exps.insert(exps.end(), s.exps.begin(), s.exps.end());
    };
    for (const Shell& s : gs.osh) addShell(s);
    for (const Shell& s : m_aux) addShell(s);
    // primitive pairs and Hermite tables
    std::vector<int> pp_ia, pp_ib, pair_pp_off(npairs), pair_npp(npairs), aux_pp_off(naux_sh), aux_npp(naux_sh), pair_A(npairs),
        pair_B(npairs);
    std::vector<double> pp_geom, E;
    std::vector<long long> pp_E_off;
    auto addPrims = [&](const Pair& pr) {
        for (const PrimPair& pp : pr.prims) {
            pp_ia.push_back(pp.ia);
            pp_ib.push_back(pp.ib);
            pp_geom.insert(pp_geom.end(), { pp.p, pp.Px, pp.Py, pp.Pz });
            pp_E_off.push_back((long long)E.size());
            E.insert(E.end(), pp.Ex.begin(), pp.Ex.end());
            E.insert(E.end(), pp.Ey.begin(), pp.Ey.end());
            E.insert(E.end(), pp.Ez.begin(), pp.Ez.end());
        }
    };
    for (int k = 0; k < npairs; ++k) {
        pair_A[k] = gs.oidx[k].first;
        pair_B[k] = gs.oidx[k].second;
        const Pair pr = makePairRaised(gs.osh[pair_A[k]], gs.osh[pair_B[k]]);
        pair_pp_off[k] = (int)pp_ia.size();
        pair_npp[k] = (int)pr.prims.size();
        addPrims(pr);
    }
    for (int cs = 0; cs < naux_sh; ++cs) {
        aux_pp_off[cs] = (int)pp_ia.size();
        aux_npp[cs] = (int)gs.apairs[cs].prims.size();
        addPrims(gs.apairs[cs]);
    }
    const double tables_ms = msSince(t0);
    // tasks and weight blocks (per pair in parallel, concatenated in pair order)
    std::vector<std::vector<int>> taux(npairs);
    std::vector<std::vector<double>> tw(npairs);
#pragma omp parallel num_threads(threads)
    {
        std::vector<double> blk, tmp;
#pragma omp for schedule(dynamic)
        for (int k = 0; k < npairs; ++k)
            for (int cs = 0; cs < naux_sh; ++cs) {
                const double sch = haveSchwarz ? m_osch[k] * m_asch[cs] : -1.0;
                if (!weightBlock(gs, m_aux, k, cs, n, G, c, Pc, sch, thr, blk, tmp)) continue;
                taux[k].push_back(cs);
                tw[k].insert(tw[k].end(), blk.begin(), blk.end());
            }
    }
    std::vector<long long> task_start(npairs + 1, 0), task_woff;
    std::vector<int> task_aux;
    std::vector<long long> pair_woff(npairs + 1, 0);
    for (int k = 0; k < npairs; ++k) {
        const int A = pair_A[k], B = pair_B[k];
        const size_t bsz = gs.osh[A].lmn.size() * gs.osh[B].lmn.size();
        size_t off = 0;
        for (int cs : taux[k]) {
            task_aux.push_back(cs);
            task_woff.push_back(pair_woff[k] + (long long)off);
            off += bsz * m_aux[cs].lmn.size();
        }
        pair_woff[k + 1] = pair_woff[k] + (long long)tw[k].size();
        task_start[k + 1] = (long long)task_aux.size();
    }
    std::vector<double> W((size_t)pair_woff[npairs]);
#pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int k = 0; k < npairs; ++k) {
        std::copy(tw[k].begin(), tw[k].end(), W.begin() + pair_woff[k]);
        std::vector<double>().swap(tw[k]);
    }
    const double prep_ms = msSince(t0);

    qmint::QMGpuRIGradInput in;
    in.natoms = natoms;
    in.nshell = nsh + naux_sh;
    in.aux_base = nsh;
    in.sh_L = sh_L.data();
    in.sh_atom = sh_atom.data();
    in.sh_ncomp = sh_ncomp.data();
    in.sh_lmn_off = sh_lmn_off.data();
    in.lmn = lmn.data();
    in.nlmn = (int)lmn.size();
    in.sh_nprim = sh_nprim.data();
    in.sh_coef_off = sh_coef_off.data();
    in.coef = coef.data();
    in.ncoef = (int)coef.size();
    in.sh_exp_off = sh_exp_off.data();
    in.exps = exps.data();
    in.nexp = (int)exps.size();
    in.npairs = npairs;
    in.pair_A = pair_A.data();
    in.pair_B = pair_B.data();
    in.pair_pp_off = pair_pp_off.data();
    in.pair_npp = pair_npp.data();
    in.aux_pp_off = aux_pp_off.data();
    in.aux_npp = aux_npp.data();
    in.npp = (int)pp_ia.size();
    in.pp_ia = pp_ia.data();
    in.pp_ib = pp_ib.data();
    in.pp_geom = pp_geom.data();
    in.pp_E_off = pp_E_off.data();
    in.nE = (long long)E.size();
    in.E = E.data();
    in.task_start = task_start.data();
    in.ntask = (long long)task_aux.size();
    in.task_aux = task_aux.data();
    in.task_woff = task_woff.data();
    in.nW = (long long)W.size();
    in.W = W.data();
    Matrix g = Matrix::Zero(natoms, 3);  // row-major natoms x 3, as the device writes it
    qmint::QMGpuRIGradOutput out;
    if (!gpu_plugin::qmRiGrad(backend, in, g.data(), out, device, why)) return false;
    grad = g;
    info = fmt::format("device: host prep {:.1f} ms (tables {:.1f} ms; {} tasks, {:.1f} MB weights), upload {:.1f} ms, kernel {:.1f} ms",
                       prep_ms, tables_ms,
                       (long long)task_aux.size(), W.size() * 8.0 / 1.0e6, out.ms_upload, out.ms_kernel);
    return true;
}

Matrix RIJK::metricGradient(const Eigen::MatrixXd& M2, int natoms) const
{
    const int threads = std::max(1, m_opt.threads);
    const OutMap ao = auxOutMap(m_aux, m_opt.spherical_aux);
    const int ns = (int)m_aux.size();
    std::vector<Eigen::MatrixXd> aT(ns);
    for (int s = 0; s < ns; ++s)
        if (ao.T[s].size()) aT[s] = ao.T[s].transpose();
    std::vector<Shell> units;
    units.reserve(m_aux.size());
    for (const Shell& s : m_aux) units.push_back(unitShell(s));
    std::vector<Pair> apairs, apairsR;
    apairs.reserve(ns);
    apairsR.reserve(ns);
    for (int s = 0; s < ns; ++s) {
        apairs.push_back(makePair(m_aux[s], units[s]));
        apairsR.push_back(makePairRaised(m_aux[s], units[s]));
    }
    const int nchunk = std::max(1, std::min(256, ns));  // Oct 2026 (WP8): thread-count independent sum
    std::vector<Matrix> gt(nchunk, Matrix::Zero(natoms, 3));
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (int ch = 0; ch < nchunk; ++ch) {
        Matrix& g = gt[ch];
        std::vector<double> blk, tmp;
        for (int cs = ch; cs < ns; cs += nchunk) {
            const int ac = m_aux[cs].atom;
            for (int ds = cs + 1; ds < ns; ++ds) {
                const int ad = m_aux[ds].atom;
                if (ac == ad) continue;  // same centre: no gradient
                const int wc = ao.width[cs], wd = ao.width[ds];
                blk.resize((size_t)wc * wd);
                for (int i = 0; i < wc; ++i)
                    for (int j = 0; j < wd; ++j) blk[(size_t)i * wd + j] = M2(ao.first[cs] + i, ao.first[ds] + j);
                int d[3] = { wc, 1, wd };
                if (aT[cs].size()) transformAxis3(blk, tmp, d, 0, aT[cs]);
                if (aT[ds].size()) transformAxis3(blk, tmp, d, 2, aT[ds]);
                double gA[3], gB[3];
                derivBlock(apairsR[cs], apairs[ds], blk.data(), false, gA, gB);
                // (C,D) and (D,C) blocks equal (symmetric M2); d/dD = -d/dC
                for (int x = 0; x < 3; ++x) {
                    g(ac, x) += 2.0 * gA[x];
                    g(ad, x) -= 2.0 * gA[x];
                }
            }
        }
    }
    Matrix grad = Matrix::Zero(natoms, 3);
    for (int t = 0; t < nchunk; ++t) grad += gt[t];
    return grad;
}

}  // namespace qmri

// =================================================================================
// 5. QMEngine glue (Claude Generated, Oct 2026, WP10)
// =================================================================================

namespace {
std::string optString(const json& config, const char* key, const std::string& fallback)
{
    if (!config.contains(key)) return fallback;
    const json& v = config[key];
    if (v.is_string()) return v.get<std::string>();
    if (v.is_boolean()) return v.get<bool>() ? "on" : "off";
    return fallback;
}
std::string lowerStr(std::string s)
{
    for (char& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}
double msGlue(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}
}  // namespace

void QMEngine::readRIConfig(const json& config)
{
    m_ri_mode = lowerStr(optString(config, "ri", m_ri_mode));
    if (m_ri_mode == "on" || m_ri_mode == "true") m_ri_mode = "jk";
    if (m_ri_mode == "false" || m_ri_mode == "none") m_ri_mode = "off";
    if (m_ri_mode != "off" && m_ri_mode != "j" && m_ri_mode != "jk") {
        CurcumaLogger::warn(fmt::format("QM: unknown -qm.ri '{}' (use off|j|jk); using off", m_ri_mode));
        m_ri_mode = "off";
    }
    m_ri_aux = optString(config, "ri_aux", m_ri_aux);
    m_ri_metric = lowerStr(optString(config, "ri_metric", m_ri_metric));
    m_gpu_ri = lowerStr(optString(config, "gpu_ri", m_gpu_ri));
    m_gpu_ri_grad = lowerStr(optString(config, "gpu_ri_grad", m_gpu_ri_grad));
    m_gpu_ri_fit = lowerStr(optString(config, "gpu_ri_fit", m_gpu_ri_fit));
    if (config.contains("ri_metric_cutoff") && config["ri_metric_cutoff"].is_number())
        m_ri_metric_cutoff = config["ri_metric_cutoff"].get<double>();
    if (config.contains("ri_max_memory_mb") && config["ri_max_memory_mb"].is_number())
        m_ri_max_memory_mb = config["ri_max_memory_mb"].get<double>();
    // Resource-saving RI-JK (Claude Generated, Oct 2026)
    m_ri_storage = lowerStr(optString(config, "ri_storage", m_ri_storage));
    m_gpu_ri_stream = lowerStr(optString(config, "gpu_ri_stream", m_gpu_ri_stream));
    m_ri_mixed_precision = lowerStr(optString(config, "ri_mixed_precision", m_ri_mixed_precision));
    if (m_ri_mixed_precision == "on") m_ri_mixed_precision = "true";
    if (m_ri_mixed_precision == "off") m_ri_mixed_precision = "false";
    if (m_ri_mixed_precision != "true" && m_ri_mixed_precision != "false" && m_ri_mixed_precision != "auto") {
        CurcumaLogger::warn(fmt::format("QM: unknown -qm.ri_mixed_precision '{}' (use true|false|auto); using false", m_ri_mixed_precision));
        m_ri_mixed_precision = "false";
    }
    auto num = [&](const char* key, double& v) {
        if (config.contains(key) && config[key].is_number()) v = config[key].get<double>();
    };
    num("ri_batch_mb", m_ri_batch_mb);
    num("gpu_ri_chunk_mb", m_gpu_ri_chunk_mb);
    num("ri_fp32_switch", m_ri_fp32_switch);
    num("ri_fit_radius", m_ri_fit_radius);
}

namespace {
/// MemAvailable of /proc/meminfo in MB, or -1 when it cannot be read (non-Linux). Claude Generated (Oct 25 2026).
double hostMemAvailableMB()
{
    std::ifstream in("/proc/meminfo");
    std::string key;
    double kb = 0.0;
    std::string unit;
    while (in >> key >> kb >> unit)
        if (key == "MemAvailable:") return kb / 1024.0;
    return -1.0;
}
}  // namespace

// Claude Generated (Oct 25 2026): -qm.ri_max_memory_mb 0 (default) = what the machine can hold now, read when the RI
// setup runs: MemAvailable minus max(8000 MB, 15 %) as margin for the rest of the process and the system. Replaces the
// fixed 16000 MB that sent C40H82 RI-JK into batched recomputation and stopped RI-MP2 (docs/QM_TIMING_ORCA_PYSCF_2026-10.md).
// The auto value is read ONCE per RI setup (cleared by resetRI) and then kept: MemAvailable shrinks by the 3-index
// tensor as soon as it is allocated, and a later check of the peak (tensor + (ia|P) factors in RI-MP2) against a
// re-read value counted the tensor twice (RI-MP2 C80H162: "needs 160572 MB, limit 100731 MB" although the setup
// limit was 209778 MB). The peak includes the tensor, so it belongs to the value measured before the tensor existed.
double QMEngine::riMemoryLimitMB() const
{
    if (m_ri_max_memory_mb > 0.0) return m_ri_max_memory_mb;
    if (m_ri_auto_limit_mb > 0.0) return m_ri_auto_limit_mb;
    const double avail = hostMemAvailableMB();
    if (avail <= 0.0) return 16000.0;
    m_ri_auto_limit_mb = std::max(1000.0, avail - std::max(8000.0, 0.15 * avail));
    m_ri_auto_avail_mb = avail;
    return m_ri_auto_limit_mb;
}

std::string QMEngine::riMemoryLimitInfo() const
{
    if (m_ri_max_memory_mb > 0.0) return fmt::format("-qm.ri_max_memory_mb {:g}", m_ri_max_memory_mb);
    const double limit = riMemoryLimitMB();   // fixes the value on the first call
    if (m_ri_auto_avail_mb <= 0.0) return "-qm.ri_max_memory_mb auto = 16000 MB (/proc/meminfo not readable)";
    return fmt::format("-qm.ri_max_memory_mb auto = {:.0f} MB (MemAvailable {:.0f} MB at the start of the RI setup minus the margin max(8000 MB, 15 %))",
                       limit, m_ri_auto_avail_mb);
}

void QMEngine::resetRI()
{
    m_ri_auto_limit_mb = -1.0;   // re-read MemAvailable for the next setup
    m_ri_auto_avail_mb = -1.0;
    m_ri.reset();
    m_ri_failed = false;
    m_gpu_ri_dev.reset();
    m_gpu_ri_checked = false;
    m_ri_occ_factor.resize(0, 0);
    m_ri_occ_density.resize(0, 0);
    m_ri_mp_stage = 0;
    m_ri_mp_refresh = m_ri_mp_last_exact = false;
    m_ri_mp_dJ.resize(0, 0);
    m_ri_mp_dK.resize(0, 0);
}

// The orbital factor stored by runSCF if P is exactly the density it belongs to (bitwise; a damped or
// extrapolated density is not), else nullptr (the J/K build then factorises P itself). Claude Generated (Oct 2026).
const Eigen::MatrixXd* QMEngine::riOccFactor(const Matrix& P) const
{
    if (m_ri_occ_factor.size() == 0 || m_ri_occ_density.rows() != P.rows() || m_ri_occ_density.cols() != P.cols()) return nullptr;
    if (!(m_ri_occ_density.array() == P.array()).all()) return nullptr;
    return &m_ri_occ_factor;
}

namespace {
// Auxiliary shells of the molecule for an auxiliary basis file (geometry in Angstrom). Claude Generated (Oct 2026).
std::vector<qmri::Shell> riAuxShells(const std::string& file, const std::vector<int>& atoms, const Matrix& geomAngstrom)
{
    const qmri::AuxBasisMap aux = qmri::parseTurbomoleBasis(file, qmbasis::parseMaxZ(atoms));
    Matrix pos = Matrix::Zero((Eigen::Index)atoms.size(), 3);
    for (int i = 0; i < (int)atoms.size(); ++i)
        for (int k = 0; k < 3; ++k) pos(i, k) = CurcumaUnit::Length::angstrom_to_bohr(geomAngstrom(i, k));
    return qmri::auxiliaryShells(aux, atoms, pos);
}

// Number of active auxiliary functions (spherical: 2L+1 for L >= 2).
long riAuxCount(const std::vector<qmri::Shell>& sh, bool spherical)
{
    long N = 0;
    for (const auto& s : sh) N += (spherical && s.L >= 2) ? 2 * s.L + 1 : (long)s.lmn.size();
    return N;
}

// Batched RI-JK exchange on the device (-gpu, batched storage, Cholesky metric): the hook creates the solver on its
// first call; any failure is a counted fallback "qm_ri" and the CPU takes over. Claude Generated (Oct 2026).
void installDeviceExchange(const std::shared_ptr<qmri::RIJK>& ri, const std::string& backend, int device, double max_dev_mb, int threads)
{
    if (!ri || !(ri->batched() || ri->fp32()) || ri->metricCholesky().size() == 0 || ri->mode() != "jk") return;
    auto solver = std::make_shared<std::shared_ptr<qmint::QMGpuRIYSolver>>();
    auto dead = std::make_shared<bool>(false);
    qmri::RIJK* raw = ri.get();
    ri->setDeviceExchange([raw, solver, dead, backend, device, max_dev_mb, threads](const Matrix& Y, int rb, double sign, Eigen::MatrixXd& K) {
        if (*dead) return false;
        if (!*solver) {
            std::string why;
            const Eigen::MatrixXd& L = raw->metricCholesky();
            auto s = gpu_plugin::createQmRIYSolver(backend, L.data(), (int)L.rows(), raw->n(), device, max_dev_mb, threads, why);
            if (!s) {
                *dead = true;
                curcuma::reportGpuFallback("qm_ri", "device batched RI exchange not possible: " + why);
                return false;
            }
            *solver = std::shared_ptr<qmint::QMGpuRIYSolver>(std::move(s));
            if (CurcumaLogger::get_verbosity() >= 2) CurcumaLogger::info("QM: batched RI exchange (L^-1 Y, SYRK) on " + (*solver)->name());
        }
        std::vector<double> Kr((size_t)K.rows() * K.cols(), 0.0);
        if (!(*solver)->exchange(Y.data(), rb, sign, Kr.data())) {
            *dead = true;
            curcuma::reportGpuFallback("qm_ri", "device batched RI exchange failed");
            return false;
        }
        K += Eigen::Map<const Eigen::MatrixXd>(Kr.data(), K.rows(), K.cols());  // symmetric: row/column major alike
        return true;
    });
}
}  // namespace

// Claude Generated (Oct 2026): two-electron part of the RI gradient (qmri::RIJK gradientDensities, threeIndexGradient[Device], metricGradient; derivation
// in section 6 above); mode j adds the exact exchange gradient (same canonical-quartet kernel as the
// exact path, exchange density only). Pc: cartesian density for the exact part.
// On failure `info` holds the reason. The gradient needs the stored FP64 tensor: after a batched or FP32 SCF a
// separate stored setup is made for it (same aux basis and options; error if it does not fit, Oct 2026).
bool QMEngine::riGradientTwoElectron(const Matrix& Pc, Matrix& g2, std::string& info) const
{
    // WP8 (Oct 2026): Kohn-Sham weights the exchange part with a_x (0 for LDA/GGA, 0.2 for B3LYP)
    const double kx = isKohnSham() ? m_xc.exact_exchange : 1.0;
    Matrix g;
    if (!riGradientForDensity(m_density, kx, g, info)) return false;
    if (m_ri_mode == "j" && kx != 0.0) {
        const auto tk = std::chrono::steady_clock::now();
        if (kx == 1.0)
            g += qmint::gradientTwoElectronExchange(m_gto_basis, Pc, m_atomcount, m_threads, m_eri_screening,
                                                    m_Q.size() != 0 ? &m_Q : nullptr);
        else
            g += kx * qmint::gradientTwoElectronExchange(m_gto_basis, Pc, m_atomcount, m_threads, m_eri_screening,
                                                         m_Q.size() != 0 ? &m_Q : nullptr);
        info += fmt::format(", exact exchange {:.1f} ms", msGlue(tk));
    }
    g2 = g;
    return true;
}

// Claude Generated (Nov 2026, stage 5b): unrestricted RI gradient. The RI energy is quadratic in each density, so
//   E = 1/2 Tr(Pt J[Pt]) - (kx/2) sum_s Tr(Ps K[Ps]),
// and D(P, kx) = riGradientForDensity(P, kx) is the gradient of 1/2 Tr(P J[P]) - (kx/4) Tr(P K[P]); since
// D(2 Ps, kx) - D(2 Ps, 0) = grad(-kx Tr(Ps K[Ps])),
//   grad E = D(Pt, 0) + 1/2 sum_s [ D(2 Ps, kx) - D(2 Ps, 0) ],
// five RI gradient passes on positive semidefinite densities. RI-J (mode j): Coulomb fitted (D(Pt, 0)), exchange exact
// (gradientTwoElectronUHF with the Coulomb part removed). Pac, Pbc: cartesian spin densities.
bool QMEngine::riGradientTwoElectronUHF(const Matrix& Pac, const Matrix& Pbc, Matrix& g2, std::string& info) const
{
    const double kx = isKohnSham() ? m_xc.exact_exchange : 1.0;
    const Matrix Pt = m_density_a + m_density_b;
    Matrix g;
    std::string i0;
    if (!riGradientForDensity(Pt, 0.0, g, i0)) {
        info = i0;
        return false;
    }
    info = "J: " + i0;
    if (kx == 0.0) {
        g2 = g;
        return true;
    }
    if (m_ri_mode == "j") {
        // exact exchange part only: UHF weight minus its Coulomb part
        const Matrix* Qs = m_Q.size() != 0 ? &m_Q : nullptr;
        const Matrix gx = qmint::gradientTwoElectronUHF(m_gto_basis, Pac, Pbc, kx, m_atomcount, m_threads, m_eri_screening, Qs)
            - qmint::gradientTwoElectronUHF(m_gto_basis, Pac, Pbc, 0.0, m_atomcount, m_threads, m_eri_screening, Qs);
        g2 = g + gx;
        info += ", exact exchange per spin";
        return true;
    }
    for (const Matrix* Ps : { &m_density_a, &m_density_b }) {
        Matrix gk, g0;
        std::string ik;
        const Matrix P2 = 2.0 * (*Ps);
        if (!riGradientForDensity(P2, kx, gk, ik) || !riGradientForDensity(P2, 0.0, g0, ik)) {
            info = ik;
            return false;
        }
        g += 0.5 * (gk - g0);
    }
    info += ", exchange per spin (4 more RI passes)";
    g2 = g;
    return true;
}

// The RI part of the gradient of 1/2 J[P].P - (kx/4) K[P].P for the active-basis density P (mode j: Coulomb only, the
// caller adds the exact exchange). On failure `info` holds the reason. The gradient needs the stored FP64 tensor: after
// a batched or FP32 SCF a separate stored setup is made for it (same aux basis and options; error if it does not fit,
// Oct 2026). Body of the former riGradientTwoElectron (Oct 2026), split out unchanged in Nov 2026 (stage 5b).
bool QMEngine::riGradientForDensity(const Matrix& P, double kx, Matrix& g2, std::string& info) const
{
    if (!ensureRI()) {
        info = "RI gradient: density fitting setup failed (-qm.ri " + m_ri_mode + ")";
        return false;
    }
    if (m_ri_fit_radius > 0.0 && m_ri_mode == "jk") {
        info = "RI gradient not implemented for the local fit (-qm.ri_fit_radius > 0): the SCF energy is not the global-fit RI energy";
        return false;
    }
    std::shared_ptr<qmri::RIJK> ri = m_ri;
    std::string rebuilt;
    if (!m_ri->storedFP64()) {
        const auto tr = std::chrono::steady_clock::now();
        std::string why;
        try {
            qmri::RIOptions opt;
            opt.mode = m_ri_mode;
            opt.metric = m_ri_metric;
            opt.metric_cutoff = m_ri_metric_cutoff;
            opt.screening = m_eri_screening;
            opt.max_memory_mb = riMemoryLimitMB();
            opt.spherical_aux = !m_cartesian_d;
            opt.threads = m_threads;
            opt.storage = "stored";
            auto r = std::make_shared<qmri::RIJK>();
            if (r->setup(m_gto_basis, m_Q.size() != 0 ? &m_Q : nullptr, riAuxShells(m_ri_file, m_atoms, m_geometry), opt, why)) ri = r;
        } catch (const std::exception& e) {
            why = e.what();
        }
        if (ri == m_ri) {
            info = "RI gradient not possible: the SCF ran with ri_storage " + m_ri->storageUsed() +
                   " and the gradient needs the stored FP64 tensor: " + why;
            return false;
        }
        rebuilt = fmt::format(", stored FP64 tensor rebuilt for the gradient in {:.1f} ms", msGlue(tr));
    }
    const auto t0 = std::chrono::steady_clock::now();
    std::string why;
    Eigen::VectorXd c;
    Matrix G;
    Eigen::MatrixXd M2;
    const Matrix* Q = m_Q.size() != 0 ? &m_Q : nullptr;
    if (!ri->gradientDensities(P, kx, c, G, M2, why)) {
        info = "RI gradient not possible: " + why;
        return false;
    }
    const auto t1 = std::chrono::steady_clock::now();
    // 3-index part: device (-qm.gpu_ri_grad) or CPU; any device failure is a counted fallback.
    Matrix g3;
    std::string dev_info;
    bool on_device = false;
    const std::string be = gpuGradBackend(m_gpu_ri_grad, "gpu_ri_grad", "qm_ri_grad", m_gpu_ri_grad_failed);
    if (!be.empty()) {
        std::string dwhy;
        on_device = ri->threeIndexGradientDevice(m_gto_basis, Q, c, P, G, m_atomcount, be, m_gpu_device, g3, dev_info, dwhy);
        if (!on_device) {
            m_gpu_ri_grad_failed = true;
            curcuma::reportGpuFallback("qm_ri_grad", "device RI gradient not possible: " + dwhy);
        }
    }
    if (!on_device) g3 = ri->threeIndexGradient(m_gto_basis, Q, c, P, G, m_atomcount);
    const double t_3c = msGlue(t1);
    if (on_device) {
        if (const char* chk = std::getenv("CURCUMA_QM_GPU_RI_GRAD_CHECK"); chk && chk[0] == '1') {
            const auto c0 = std::chrono::steady_clock::now();
            const Matrix gc = ri->threeIndexGradient(m_gto_basis, Q, c, P, G, m_atomcount);
            CurcumaLogger::result(fmt::format("QM GPU RI grad check: natoms={} max|dev-cpu|={:.3e} max|cpu|={:.3e} dev_ms={:.1f} cpu_ms={:.1f}",
                                              m_atomcount, (g3 - gc).cwiseAbs().maxCoeff(), gc.cwiseAbs().maxCoeff(), t_3c, msGlue(c0)));
        }
    }
    Matrix().swap(G);
    const auto t2 = std::chrono::steady_clock::now();
    Matrix g = g3 - ri->metricGradient(M2, m_atomcount);
    const double t_2c = msGlue(t2), t_dens = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double t_ri = msGlue(t0);
    g2 = g;
    info = fmt::format("RI-{}: densities {:.1f} ms, 3-index {:.1f} ms {}, metric {:.1f} ms (total {:.1f} ms){}",
                       m_ri_mode == "j" ? "J" : "JK", t_dens, t_3c, on_device ? "(" + dev_info + ")" : std::string("(CPU)"), t_2c, t_ri,
                       rebuilt);
    return true;
}

// One line per SCF at verbosity 2: number and time of the RI J/K builds (diagnostic).
void QMEngine::riReport() const
{
    if (m_ri && CurcumaLogger::get_verbosity() >= 2 && m_ri_builds > 0) {
        std::string dev = "CPU";
        if (m_gpu_ri_dev) {
            dev = "contractions on " + m_gpu_ri_dev->name();
            if (auto* lm = dynamic_cast<qmint::QMGpuRILowMem*>(m_gpu_ri_dev.get())) dev += " (" + lm->stats() + ")";
        }
        CurcumaLogger::info(fmt::format("QM: RI-{} J/K builds: {} in {:.1f} ms ({:.2f} ms each; density factorisation {:.1f} ms; {}; storage {})",
                                        m_ri_mode == "j" ? "J" : "JK", m_ri_builds, m_ri_build_ms, m_ri_build_ms / m_ri_builds,
                                        m_ri_factor_ms + m_ri->takeFactorMs(), dev, m_ri->storageUsed()));
        if (m_ri_mp_exact_builds > 0)
            CurcumaLogger::info(fmt::format("QM: RI mixed precision: {} iterations on the FP32 tensor, {} exact FP64 builds (corrections)", m_ri_fp32_builds,
                                            m_ri_mp_exact_builds));
        const std::string bs = m_ri->takeBuildStats();
        if (!bs.empty()) CurcumaLogger::info("QM: RI " + bs);
    }
    m_ri_build_ms = 0.0;
    m_ri_factor_ms = 0.0;
    m_ri_builds = 0;
}

bool QMEngine::riActive() const { return m_ri_mode == "j" || m_ri_mode == "jk"; }

// Mixed precision (Claude Generated, Oct 2026). Stage 1: FP32 tensor, FP32 GEMM/SYRK (fast, noisy: ||dP|| stalls near
// 1e-5 for benzene). Stage 2: the same FP32 tensor with FP64 arithmetic (a smooth map) plus the correction
// dJ = J_exact(P_ref) - J_32(P_ref), dK likewise, from one exact batched FP64 build at P_ref. The correction depends on
// P only through the rounding of B (relative 6e-8), so F_32(P) + dF(P_ref) differs from F_exact(P) by ~6e-8 |F(P - P_ref)|.
// Convergence is accepted only on an exact build: when stage 2 reaches ||dP|| < scf_threshold the next build is exact
// (refreshing the correction), and the final build of finishConvergedSCF is exact as well. runSCF never accepts an
// iteration whose Fock matrix did not come from an exact FP64 build.
bool QMEngine::riLowPrecision() const { return m_ri && m_ri_mp_stage > 0 && !m_ri_mp_last_exact; }

void QMEngine::riRequestExactBuild() const
{
    if (m_ri_mp_stage > 0) m_ri_mp_refresh = true;
}

bool QMEngine::riEndLowPrecision(double delta_P) const
{
    if (!m_ri || m_ri_mp_stage == 0) return false;
    ++m_ri_fp32_builds;
    if (m_ri_fp32_best <= 0.0 || delta_P < 0.5 * m_ri_fp32_best) {
        m_ri_fp32_best = delta_P;
        m_ri_fp32_stall = 0;
    } else
        ++m_ri_fp32_stall;
    const bool stalled = m_ri_fp32_stall >= 3;
    if (m_ri_mp_stage == 1) {
        if (!(delta_P < m_ri_fp32_switch || delta_P < m_scf_threshold || stalled)) return false;
        m_ri_mp_stage = 2;
        m_ri_mp_refresh = true;
        m_ri->setFP32Arithmetic(false);
        if (auto* lm = dynamic_cast<qmint::QMGpuRILowMem*>(m_gpu_ri_dev.get())) lm->setFP64Arithmetic(true);
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info(fmt::format("QM: RI mixed precision: ||dP|| = {:.3e} after {} FP32 iterations{}; FP64 arithmetic on the FP32 tensor "
                                            "plus exact (batched) corrections from now on",
                                            delta_P, m_ri_fp32_builds, stalled ? " (stalled)" : ""));
        m_ri_fp32_best = 0.0;
        m_ri_fp32_stall = 0;
        return true;
    }
    // stage 2: verify a converged corrected map with an exact build (and refresh after a stall)
    if ((delta_P < m_scf_threshold && !m_ri_mp_last_exact) || stalled) {
        m_ri_mp_refresh = true;
        m_ri_fp32_stall = 0;
    }
    return false;
}

bool QMEngine::ensureRI() const
{
    if (m_ri) return true;
    if (m_ri_failed) return false;
    const auto t0 = std::chrono::steady_clock::now();
    std::string why;
    try {
        const std::string file = qmri::resolveAuxBasisFile(m_ri_aux, m_ri_mode);
        const std::vector<qmri::Shell> sh = riAuxShells(file, m_atoms, m_geometry);
        qmri::RIOptions opt;
        opt.mode = m_ri_mode;
        opt.metric = m_ri_metric;
        opt.metric_cutoff = m_ri_metric_cutoff;
        opt.screening = m_eri_screening;
        const double ri_limit = riMemoryLimitMB();  // -qm.ri_max_memory_mb, 0 = auto from MemAvailable (Oct 25 2026)
        opt.max_memory_mb = ri_limit;
        if (CurcumaLogger::get_verbosity() >= 2) CurcumaLogger::info("QM: RI host memory limit " + riMemoryLimitInfo());
        opt.spherical_aux = !m_cartesian_d;
        opt.threads = m_threads;
        opt.storage = m_ri_storage;
        opt.batch_mb = m_ri_batch_mb;
        opt.fit_radius = m_ri_fit_radius;
        if (m_ri_fit_radius > 0.0 && m_ri_mode == "jk" && m_ri_mixed_precision == "true")
            CurcumaLogger::warn("QM: -qm.ri_fit_radius > 0 runs batched FP64 builds; -qm.ri_mixed_precision is ignored");
        // -qm.ri_mixed_precision (Claude Generated, Oct 2026): true = FP32 tensor; auto = FP32 only when the FP64
        // tensor does not fit (host limit, or -qm.gpu_max_memory_mb with -gpu) but the FP32 one does.
        bool fp32 = m_ri_mode == "jk" && m_ri_mixed_precision == "true" && !(m_ri_fit_radius > 0.0);
        if (m_ri_mode == "jk" && m_ri_mixed_precision == "auto" && !(m_ri_fit_radius > 0.0)) {
            const double np = (double)m_nbf * (m_nbf + 1) / 2.0, N = (double)riAuxCount(sh, !m_cartesian_d);
            const double mb64 = N * np * 8.0 / 1.0e6, mb32 = mb64 / 2.0;
            const bool dev_cap = gpuRequested() && m_gpu_max_memory_mb > 0.0;
            fp32 = (mb64 > ri_limit && mb32 <= ri_limit) ||
                   (dev_cap && mb64 > m_gpu_max_memory_mb && mb32 <= m_gpu_max_memory_mb && mb32 <= ri_limit);
        }
        opt.fp32 = fp32;
        // -qm.gpu_ri_fit (Claude Generated, Oct 2026): the B = L^-1 T solve on the device when the RI
        // contractions will run there too. Failure before any write: host solve; after a partial write:
        // setup repeated without the device. Both counted as qm_ri_fit.
        const bool off = [](const std::string& v) { return v == "off" || v == "false"; }(m_gpu_ri_fit);
        const bool ri_off = m_gpu_ri == "off" || m_gpu_ri == "false";
        std::string fit_backend;
        if (m_ri_mode == "jk" && !off && !ri_off && !m_gpu_ri_fit_failed && gpuRequested()) {
            fit_backend = m_gpu;
            if (fit_backend == "none" || !gpu_plugin::available(fit_backend)) {
                curcuma::reportGpuFallback("qm_ri_fit", "no GPU plugin for -gpu " + m_gpu);
                fit_backend.clear();
            }
        }
        if (!fit_backend.empty()) {
            const int device = m_gpu_device, host_threads = std::max(1, m_threads);
            opt.device_lower_solve = [fit_backend, device, host_threads](const Eigen::MatrixXd& L, Matrix& T, std::string& w) {
                bool touched = false;
                if (gpu_plugin::qmRiLowerSolve(fit_backend, T.data(), (int)T.rows(), (long)T.cols(), L.data(), device, host_threads,
                                               touched, w))
                    return 1;
                return touched ? -1 : 0;
            };
        }
        auto ri = std::make_shared<qmri::RIJK>();
        bool ok = ri->setup(m_gto_basis, m_Q.size() != 0 ? &m_Q : nullptr, sh, opt, why);
        if (opt.device_lower_solve && (!ok || ri->stats().find("solve host (device failed") != std::string::npos)) {
            m_gpu_ri_fit_failed = true;
            curcuma::reportGpuFallback("qm_ri_fit", ok ? "device metric solve failed, host solve used" : why);
            if (!ok && why.rfind("device metric solve", 0) == 0) {
                opt.device_lower_solve = nullptr;
                ri = std::make_shared<qmri::RIJK>();
                ok = ri->setup(m_gto_basis, m_Q.size() != 0 ? &m_Q : nullptr, sh, opt, why);
            }
        }
        if (ok) {
            m_ri = ri;
            m_ri_file = file;
        }
    } catch (const std::exception& e) {
        why = e.what();
    }
    if (!m_ri) {
        m_ri_failed = true;
        CurcumaLogger::error("QM: density fitting (-qm.ri " + m_ri_mode + ") not possible: " + why + " [limit: " + riMemoryLimitInfo() + "]");
        return false;
    }
    m_ri_mp_stage = m_ri->fp32() ? 1 : 0;
    m_ri_mp_refresh = m_ri_mp_last_exact = false;
    m_ri_mp_exact_builds = 0;
    m_ri_fp32_best = 0.0;
    m_ri_fp32_stall = 0;
    m_ri_fp32_builds = 0;
    if ((m_ri->batched() || m_ri->fp32()) && gpuRequested() && m_gpu_ri != "off" && m_gpu_ri != "false" && !m_gpu_ri_failed) {
        const std::string backend = m_gpu;
        if (backend == "none" || !gpu_plugin::available(backend))
            curcuma::reportGpuFallback("qm_ri", "no GPU plugin for -gpu " + m_gpu);
        else
            installDeviceExchange(m_ri, backend, m_gpu_device, m_gpu_max_memory_mb, m_threads);
    }
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("QM: RI-{} with {} ({:.1f} ms total)", m_ri_mode == "j" ? "J" : "JK", m_ri_file, msGlue(t0)));
        CurcumaLogger::info("QM: RI " + m_ri->stats());
    }
    if (const char* dump = std::getenv("CURCUMA_QM_RI_DUMP"); dump && dump[0] != '\0') riDump(dump);
    return true;
}

// CURCUMA_QM_RI_DUMP=<prefix>: writes <prefix>_metric.txt (raw (P|Q)), <prefix>_3c.txt (raw (P|mn), active
// orbital basis, row P = n*n values), <prefix>_ovl.txt (active overlap). Diagnostic for the PySCF comparison.
void QMEngine::riDump(const std::string& prefix) const
{
    auto write = [](const std::string& path, const Matrix& M) {
        std::ofstream f(path);
        f << M.rows() << " " << M.cols() << "\n";
        f.precision(17);
        for (Eigen::Index i = 0; i < M.rows(); ++i) {
            for (Eigen::Index j = 0; j < M.cols(); ++j) f << M(i, j) << (j + 1 < M.cols() ? " " : "");
            f << "\n";
        }
    };
    write(prefix + "_metric.txt", m_ri->rawMetric());
    write(prefix + "_3c.txt", m_ri->rawThreeIndex());
    write(prefix + "_ovl.txt", m_S);
}

// J(P), K(P): exact path (stored tensor / integral-direct / device, qm_scf.cpp buildJKExact) or the RI fit.
// RI-J: Coulomb from the fit, exchange from the exact path (the exact path forms J alongside, so RI-J saves
// no time for HF; it is the Coulomb part a pure functional would need). RI-JK: both from the fit.
void QMEngine::buildJK(const Matrix& P, Matrix& J, Matrix& K) const
{
    if (!riActive()) {
        buildJKExact(P, J, K);
        return;
    }
    if (!ensureRI()) {
        // runSCF checks ensureRI() first; never return exact numbers under the RI label
        const double nan = std::numeric_limits<double>::quiet_NaN();
        J = Matrix::Constant(m_nbf, m_nbf, nan);
        K = Matrix::Constant(m_nbf, m_nbf, nan);
        return;
    }
    const auto t0 = std::chrono::steady_clock::now();
    if (m_ri_mode == "j") {
        Matrix Jx;
        buildJKExact(P, Jx, K);
        if (!gpuRIBuild(P, J, nullptr)) J = m_ri->coulomb(P);
    } else {
        const Eigen::MatrixXd* Vocc = riOccFactor(P);
        m_ri_mp_last_exact = false;
        if (m_ri_mp_stage == 2 && m_ri_mp_refresh) {
            // mixed precision: exact FP64 build (batched) and the new correction exact - FP32-tensor build at this P
            Matrix J32, K32;
            if (!gpuRIBuild(P, J32, &K32, Vocc)) m_ri->build(P, J32, &K32, Vocc);
            m_ri->buildExact(P, J, K, Vocc);
            m_ri_mp_dJ = J - J32;
            m_ri_mp_dK = K - K32;
            m_ri_mp_refresh = false;
            m_ri_mp_last_exact = true;
            ++m_ri_mp_exact_builds;
            if (CurcumaLogger::get_verbosity() >= 3)
                CurcumaLogger::info(fmt::format("QM: RI mixed precision: exact build {}, max|dJ| {:.3e}, max|dK| {:.3e}", m_ri_mp_exact_builds,
                                                m_ri_mp_dJ.cwiseAbs().maxCoeff(), m_ri_mp_dK.cwiseAbs().maxCoeff()));
        } else {
            if (!gpuRIBuild(P, J, &K, Vocc)) m_ri->build(P, J, &K, Vocc);
            if (m_ri_mp_stage == 2 && m_ri_mp_dJ.size()) {
                J += m_ri_mp_dJ;
                K += m_ri_mp_dK;
            }
        }
    }
    m_ri_build_ms += msGlue(t0);
    ++m_ri_builds;
}

// -qm.gpu_ri auto|on|off: the RI contractions (J from the tensor, K as batched GEMM + SYRK) on the device.
// Any device error: counted fallback "qm_ri", CPU from then on. Since Oct 2026 the tensor may also be FP32
// (ri_mixed_precision) or streamed from host memory (-qm.gpu_ri_stream, qm_gpu_ri_lowmem.cu); batched storage
// has no tensor, there the device only takes the exchange solve (installDeviceExchange).
bool QMEngine::gpuRIBuild(const Matrix& P, Matrix& J, Matrix* K, const Eigen::MatrixXd* Vocc) const
{
    if (m_gpu_ri_failed || m_gpu_ri == "off" || m_gpu_ri == "false" || !gpuRequested()) return false;
    if (m_ri->batched()) return false;
    if (m_gpu_ri != "auto" && m_gpu_ri != "on" && m_gpu_ri != "true") {
        static bool warned = false;
        if (!warned) {
            warned = true;
            CurcumaLogger::warn(fmt::format("QM: unknown -qm.gpu_ri '{}' (use auto|on|off); using auto", m_gpu_ri));
        }
    }
    auto fail = [&](const std::string& what) {
        m_gpu_ri_failed = true;
        m_gpu_ri_dev.reset();
        curcuma::reportGpuFallback("qm_ri", what);
        return false;
    };
    if (!m_gpu_ri_dev) {
        const std::string backend = m_gpu;
        if (backend == "none" || !gpu_plugin::available(backend)) return fail("no GPU plugin for -gpu " + m_gpu);
        std::string why;
        const auto t0 = std::chrono::steady_clock::now();
        const int stream = (m_gpu_ri_stream == "on" || m_gpu_ri_stream == "true") ? 1
            : (m_gpu_ri_stream == "off" || m_gpu_ri_stream == "false")            ? 2
                                                                                    : 0;
        const bool fp32 = m_ri->fp32();
        const bool capped = m_gpu_max_memory_mb > 0.0 && m_ri->tensorMB() > m_gpu_max_memory_mb;
        const int rows = fp32 ? m_ri->nfit() : (int)m_ri->tensor().rows();
        const void* T = fp32 ? (const void*)m_ri->tensor32().data() : (const void*)m_ri->tensor().data();
        if (fp32 || stream == 1 || (capped && stream == 0)) {
            m_gpu_ri_dev = gpu_plugin::createQmRILowMem(backend, T, fp32, rows, m_ri->n(), m_gpu_device, stream, m_gpu_max_memory_mb,
                                                        m_gpu_ri_chunk_mb, m_threads, why);
        } else if (capped && stream == 2) {
            why = fmt::format("RI tensor {:.1f} MB above -qm.gpu_max_memory_mb {:g} and -qm.gpu_ri_stream off", m_ri->tensorMB(),
                              m_gpu_max_memory_mb);
        } else {
            m_gpu_ri_dev = gpu_plugin::createQmRI(backend, m_ri->tensor().data(), rows, m_ri->n(), m_gpu_device, why);
            if (!m_gpu_ri_dev && stream == 0 && why.find("does not fit") != std::string::npos) {
                std::string why2;
                m_gpu_ri_dev = gpu_plugin::createQmRILowMem(backend, T, false, rows, m_ri->n(), m_gpu_device, 1, m_gpu_max_memory_mb,
                                                            m_gpu_ri_chunk_mb, m_threads, why2);
                if (!m_gpu_ri_dev) why += "; streaming not possible: " + why2;
            }
        }
        if (!m_gpu_ri_dev) return fail("device RI setup not possible: " + why);
        if (CurcumaLogger::get_verbosity() >= 2) {
            auto* lm = dynamic_cast<qmint::QMGpuRILowMem*>(m_gpu_ri_dev.get());
            CurcumaLogger::info(fmt::format("QM: RI contractions on {} ({:.1f} MB {}tensor {} in {:.1f} ms)", m_gpu_ri_dev->name(),
                                            m_ri->tensorMB(), fp32 ? "FP32 " : "",
                                            lm && lm->streamed() ? "kept in pinned host memory, streamed per build" : "uploaded",
                                            msGlue(t0)));
        }
    }
    const int n = m_ri->n();
    auto* lowmem = dynamic_cast<qmint::QMGpuRILowMem*>(m_gpu_ri_dev.get());
    if (lowmem && K) {
        // one pass over the (streamed or FP32) tensor for g, J and K (mode jk: c = g)
        Eigen::MatrixXd Vown;
        const Eigen::MatrixXd* Vp = Vocc;
        int npos = Vocc ? (int)Vocc->cols() : 0;
        if (!Vp) {
            const auto tf = std::chrono::steady_clock::now();
            qmri::RIJK::factorDensity(P, Vown, npos);
            m_ri_factor_ms += msGlue(tf);
            Vp = &Vown;
        }
        std::vector<double> Pw, Jp((size_t)m_ri->npair()), g((size_t)m_ri->nfit());
        qmri::packContractionWeights(P, Pw);
        Matrix Kd = Matrix::Zero(n, n);
        if (!lowmem->buildJK(Pw.data(), Vp->data(), (int)Vp->cols(), npos, g.data(), Jp.data(), Kd.data()))
            return fail("device streamed/FP32 RI J/K build failed");
        J = Matrix::Zero(n, n);
        qmri::unpackSymmetric(Jp.data(), n, J.data());
        *K = Kd;
        if (!m_gpu_ri_checked) {
            m_gpu_ri_checked = true;
            if (const char* chk = std::getenv("CURCUMA_QM_GPU_RI_CHECK"); chk && chk[0] == '1') {
                Matrix Jc, Kc;
                m_ri->build(P, Jc, &Kc);
                CurcumaLogger::result(fmt::format("QM GPU RI check: n={} naux={} max|dJ|={:.3e} max|dK|={:.3e} max|J|={:.3e} max|K|={:.3e}", n,
                                                  m_ri->nfit(), (J - Jc).cwiseAbs().maxCoeff(), (Kd - Kc).cwiseAbs().maxCoeff(),
                                                  Jc.cwiseAbs().maxCoeff(), Kc.cwiseAbs().maxCoeff()));
            }
        }
        return true;
    }
    J = Matrix::Zero(n, n);
    // Coulomb on the packed tensor: g = T . w (w = packed contraction weights of P) on the device; mode j needs
    // c = J^-1 g on the host (N^2), then J_p = T^T c, unpacked on the host
    std::vector<double> Pw, Jp((size_t)m_ri->npair());
    qmri::packContractionWeights(P, Pw);
    std::vector<double> g((size_t)(m_ri->fp32() ? m_ri->nfit() : m_ri->tensor().rows()));  // mode j: naux rows (see RIJK::coulomb)
    if (!m_gpu_ri_dev->contract(Pw.data(), g.data())) return fail("device g = T w(P) failed");
    Eigen::VectorXd c = m_ri->fitCoefficients(Eigen::Map<const Eigen::VectorXd>(g.data(), (Eigen::Index)g.size()));
    if (!m_gpu_ri_dev->coulomb(c.data(), Jp.data())) return fail("device J = T^T c failed");
    qmri::unpackSymmetric(Jp.data(), n, J.data());
    if (!K) return true;
    Eigen::MatrixXd Vown;
    const Eigen::MatrixXd* Vp = Vocc;
    int npos = Vocc ? (int)Vocc->cols() : 0;
    if (!Vp) {
        const auto tf = std::chrono::steady_clock::now();
        qmri::RIJK::factorDensity(P, Vown, npos);
        m_ri_factor_ms += msGlue(tf);
        Vp = &Vown;
    }
    Matrix Kd = Matrix::Zero(n, n);
    if (!m_gpu_ri_dev->exchange(Vp->data(), (int)Vp->cols(), npos, Kd.data())) return fail("device exchange build failed");
    *K = Kd;
    if (!m_gpu_ri_checked) {
        m_gpu_ri_checked = true;
        if (const char* chk = std::getenv("CURCUMA_QM_GPU_RI_CHECK"); chk && chk[0] == '1') {
            Matrix Jc, Kc;
            m_ri->build(P, Jc, &Kc);
            CurcumaLogger::result(fmt::format("QM GPU RI check: n={} naux={} max|dJ|={:.3e} max|dK|={:.3e} max|J|={:.3e} max|K|={:.3e}", n,
                                              m_ri->nfit(), (J - Jc).cwiseAbs().maxCoeff(), (Kd - Kc).cwiseAbs().maxCoeff(),
                                              Jc.cwiseAbs().maxCoeff(), Kc.cwiseAbs().maxCoeff()));
        }
    }
    return true;
}
