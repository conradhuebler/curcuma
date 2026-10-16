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
AuxBasisMap parseTurbomoleBasis(const std::string& path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open auxiliary basis file " + path);
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
    if (name.find('/') != std::string::npos && fileExists(name)) return name;
    const std::string file = name + ".turbomole";
    std::vector<std::string> cand;
    for (const char* var : { "CURCUMA_QM_BASIS", "CURCUMA_DFT_BASIS", "CURCUMA_DATA" })
        if (const char* env = std::getenv(var); env && env[0] != '\0') cand.push_back(std::string(env) + "/" + file);
#ifdef CURCUMA_SOURCE_DIR
    cand.push_back(std::string(CURCUMA_SOURCE_DIR) + "/src/core/energy_calculators/qm_methods/" + file);
#endif
    cand.push_back(file);
    for (const auto& c : cand)
        if (fileExists(c)) return c;
    throw std::runtime_error("auxiliary basis '" + name_in + "' not found (looked for " + file + ")");
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
const char* kSymbols[] = { nullptr, "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne" };

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
        if (Z < 1 || Z > 10) throw std::runtime_error(fmt::format("RI: element Z={} outside H-Ne", Z));
        auto it = aux.find(kSymbols[Z]);
        if (it == aux.end()) throw std::runtime_error(fmt::format("RI: no auxiliary functions for element {}", kSymbols[Z]));
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
};

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

// Dense kernels on the row-major tensor (N rows of length ncol = n^2). Written as plain loops over
// column chunks (each thread owns a column range for all rows), because the build links the system
// BLAS, which may be the single-threaded reference implementation (it is on the development machine).

/// In place T := L^-1 T (L lower triangular N x N): forward substitution row by row.
void lowerSolveRows(const Eigen::MatrixXd& L, Matrix& T, int threads)
{
    const int N = (int)T.rows();
    const Eigen::Index ncol = T.cols(), chunk = 128, nchunks = (ncol + chunk - 1) / chunk;
    double* base = T.data();
#pragma omp parallel for schedule(dynamic) num_threads(std::max(1, threads))
    for (Eigen::Index b = 0; b < nchunks; ++b) {
        const Eigen::Index c0 = b * chunk, w = std::min(chunk, ncol - c0);
        for (int P = 0; P < N; ++P) {
            double* rp = base + (size_t)P * ncol + c0;
            for (int Q = 0; Q < P; ++Q) {
                const double l = L(P, Q);
                if (l == 0.0) continue;
                const double* rq = base + (size_t)Q * ncol + c0;
                for (Eigen::Index k = 0; k < w; ++k) rp[k] -= l * rq[k];
            }
            const double inv = 1.0 / L(P, P);
            for (Eigen::Index k = 0; k < w; ++k) rp[k] *= inv;
        }
    }
}

/// out := F T (F: r x N, T: N x ncol row-major).
void leftMultiplyRows(const Eigen::MatrixXd& F, const Matrix& T, Matrix& out, int threads)
{
    const int N = (int)T.rows(), r = (int)F.rows();
    const Eigen::Index ncol = T.cols(), chunk = 128, nchunks = (ncol + chunk - 1) / chunk;
    out = Matrix::Zero(r, ncol);
#pragma omp parallel for schedule(dynamic) num_threads(std::max(1, threads))
    for (Eigen::Index b = 0; b < nchunks; ++b) {
        const Eigen::Index c0 = b * chunk, w = std::min(chunk, ncol - c0);
        for (int P = 0; P < N; ++P) {
            const double* tp = T.data() + (size_t)P * ncol + c0;
            for (int i = 0; i < r; ++i) {
                const double f = F(i, P);
                double* o = out.data() + (size_t)i * ncol + c0;
                for (Eigen::Index k = 0; k < w; ++k) o[k] += f * tp[k];
            }
        }
    }
}

/// Where each shell lands in the output basis: first index, width, transform (identity -> empty).
struct OutMap {
    int n = 0;
    std::vector<int> first, width;
    std::vector<Eigen::MatrixXd> T;  ///< empty = identity
    OutMap();
    OutMap(const OutMap&);
    ~OutMap();  // out of line: the implicit inline one trips -Winline
};
OutMap::OutMap() = default;
OutMap::OutMap(const OutMap&) = default;
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

bool RIJK::setup(const std::vector<GTO::Orbital>& basis, const Matrix* Q, const std::vector<Shell>& aux,
                 const RIOptions& opt, std::string& why)
{
    m_opt = opt;
    curcuma::ScopedBlasThreads blas_guard(std::max(1, opt.threads));
    const auto t0 = std::chrono::steady_clock::now();
    const std::vector<Shell> osh = orbitalShells(basis);
    const OutMap oo = orbitalOutMap(osh, (int)basis.size(), Q);
    const OutMap ao = auxOutMap(aux, opt.spherical_aux);
    m_n = oo.n;
    m_naux = ao.n;
    const int n = m_n, N = m_naux;
    const double need_mb = (double)N * n * n * 8.0 / 1.0e6 * (opt.metric == "eigen" && opt.mode == "jk" ? 2.0 : 1.0);
    if (need_mb > opt.max_memory_mb) {
        why = fmt::format("3-index tensor needs {:.0f} MB (naux {} x n^2 {}), above -qm.ri_max_memory_mb {:g}", need_mb, N,
                          (long)n * n, opt.max_memory_mb);
        return false;
    }

    // Pairs and Schwarz factors
    std::vector<Pair> opairs;
    std::vector<std::pair<int, int>> oidx;
    for (int A = 0; A < (int)osh.size(); ++A)
        for (int B = A; B < (int)osh.size(); ++B) {
            opairs.push_back(makePair(osh[A], osh[B]));
            oidx.push_back({ A, B });
        }
    std::vector<Shell> units;
    units.reserve(aux.size());
    for (const Shell& s : aux) units.push_back(unitShell(s));
    std::vector<Pair> apairs;
    apairs.reserve(aux.size());
    for (size_t c = 0; c < aux.size(); ++c) apairs.push_back(makePair(aux[c], units[c]));
    const int threads = std::max(1, opt.threads);
    const int npairs = (int)opairs.size(), naux_sh = (int)aux.size();
#pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int k = 0; k < npairs; ++k) opairs[k].schwarz = schwarzOf(opairs[k]);
    for (Pair& p : apairs) p.schwarz = schwarzOf(p);
    const double t_schwarz = msSince(t0);

    // 3-index tensor T(P, m*n + nn) = (P|m nn) in the active orbital basis
    const auto t1 = std::chrono::steady_clock::now();
    Matrix T3 = Matrix::Zero(N, (Eigen::Index)n * n);
    long blocks = 0, skipped = 0;
#pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+ : blocks, skipped)
    for (int k = 0; k < npairs; ++k) {
        const int A = oidx[k].first, B = oidx[k].second;
        std::vector<double> buf, tmp;
        for (int c = 0; c < naux_sh; ++c) {
            if (opairs[k].schwarz * apairs[c].schwarz < opt.screening) {
                ++skipped;
                continue;
            }
            ++blocks;
            pairBlock(opairs[k], apairs[c], buf);
            int d[3] = { (int)osh[A].lmn.size(), (int)osh[B].lmn.size(), (int)aux[c].lmn.size() };
            if (oo.T[A].size()) transformAxis3(buf, tmp, d, 0, oo.T[A]);
            if (oo.T[B].size()) transformAxis3(buf, tmp, d, 1, oo.T[B]);
            if (ao.T[c].size()) transformAxis3(buf, tmp, d, 2, ao.T[c]);
            for (int a = 0; a < d[0]; ++a)
                for (int b = 0; b < d[1]; ++b) {
                    const int m = oo.first[A] + a, nn = oo.first[B] + b;
                    const double* v = &buf[((size_t)a * d[1] + b) * d[2]];
                    for (int p = 0; p < d[2]; ++p) {
                        T3(ao.first[c] + p, (Eigen::Index)m * n + nn) = v[p];
                        T3(ao.first[c] + p, (Eigen::Index)nn * n + m) = v[p];
                    }
                }
        }
    }
    const double t_3c = msSince(t1);

    // Metric
    const auto t2 = std::chrono::steady_clock::now();
    Matrix Jraw = coulombMetric(aux, opt.spherical_aux, threads);
    const double t_2c = msSince(t2);
    if (const char* dump = std::getenv("CURCUMA_QM_RI_DUMP"); dump && dump[0] != '\0') {
        m_raw_metric = Jraw;
        m_raw_T3 = T3;
    }

    // Coulomb normalisation: every aux function scaled to (P|P) = 1. The fit is invariant to it;
    // it only puts the metric eigenvalues on a common scale for the eigen cutoff and the report.
    Eigen::VectorXd dn(N);
    for (int P = 0; P < N; ++P) {
        if (!(Jraw(P, P) > 0.0)) {
            why = fmt::format("auxiliary function {} has (P|P) = {:.3e}", P, Jraw(P, P));
            return false;
        }
        dn(P) = 1.0 / std::sqrt(Jraw(P, P));
    }
    Eigen::MatrixXd Jn = dn.asDiagonal() * Eigen::MatrixXd(Jraw) * dn.asDiagonal();
    for (int P = 0; P < N; ++P) T3.row(P) *= dn(P);

    const auto t3 = std::chrono::steady_clock::now();
    std::string metric = opt.metric;
    if (metric != "cholesky" && metric != "eigen") {
        why = "unknown -qm.ri_metric '" + metric + "' (use cholesky|eigen)";
        return false;
    }
    double wmin = std::numeric_limits<double>::quiet_NaN(), wmax = wmin;
    int dropped = 0;
    m_metric_used = metric;
    if (metric == "cholesky") {
        Eigen::LLT<Eigen::MatrixXd> llt(Jn);
        if (llt.info() != Eigen::Success) {
            metric = "eigen";
            m_metric_used = "eigen (Cholesky failed)";
        } else {
            m_L = llt.matrixL();
            if (opt.mode == "jk") {
                lowerSolveRows(m_L, T3, threads);  // B = L^-1 T (forward substitution, no explicit inverse)
                m_L.resize(0, 0);
            }
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
        if (opt.mode == "jk") {
            // B = w^-1/2 U^T T  (J^-1 ~ U w^-1 U^T = (w^-1/2 U^T)^T (w^-1/2 U^T))
            const Eigen::MatrixXd F = wk.cwiseSqrt().cwiseInverse().asDiagonal() * U.transpose();
            leftMultiplyRows(F, T3, m_T, threads);
            Matrix().swap(T3);
        } else {
            m_U = U;
            m_w = wk;
            m_T.swap(T3);
        }
    }
    const double t_fac = msSince(t3);
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
        "blocks), metric {:.1f} ms, factor/B {:.1f} ms",
        N, m_nfit, n, tensorMB(), m_metric_used, cond, t_schwarz, t_3c, blocks, blocks + skipped, t_2c, t_fac);
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
    const int n = m_n, threads = std::max(1, m_opt.threads);
    const Matrix Pc = P;  // contiguous row-major
    const Eigen::Index ncol = (Eigen::Index)n * n;
    // g_Q = sum_mn T_Q,mn P_mn (one dot product per row)
    Eigen::VectorXd g(m_nfit);
#pragma omp parallel for schedule(static) num_threads(threads)
    for (int Q = 0; Q < m_nfit; ++Q) {
        const double* t = m_T.data() + (size_t)Q * ncol;
        double s = 0.0;
        for (Eigen::Index k = 0; k < ncol; ++k) s += t[k] * Pc.data()[k];
        g(Q) = s;
    }
    const Eigen::VectorXd c = fitCoefficients(g);
    // J_mn = sum_Q c_Q T_Q,mn (column chunks in parallel, Q summed in order)
    Matrix J(n, n);
    const Eigen::Index chunk = 1024, nchunks = (ncol + chunk - 1) / chunk;
#pragma omp parallel for schedule(static) num_threads(threads)
    for (Eigen::Index b = 0; b < nchunks; ++b) {
        const Eigen::Index c0 = b * chunk, w = std::min(chunk, ncol - c0);
        double* out = J.data() + c0;
        for (Eigen::Index k = 0; k < w; ++k) out[k] = 0.0;
        for (int Q = 0; Q < m_nfit; ++Q) {
            const double cq = c(Q);
            const double* t = m_T.data() + (size_t)Q * ncol + c0;
            for (Eigen::Index k = 0; k < w; ++k) out[k] += cq * t[k];
        }
    }
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

void RIJK::build(const Matrix& P, Matrix& J, Matrix* K) const
{
    J = coulomb(P);
    if (!K) return;
    if (m_opt.mode != "jk") throw std::logic_error("RIJK::build: exchange needs mode jk");
    const int n = m_n, threads = std::max(1, m_opt.threads);
    Eigen::MatrixXd V;
    int npos = 0;
    const auto tf = std::chrono::steady_clock::now();
    factorDensity(P, V, npos);
    m_factor_ms += msSince(tf);
    const int r = (int)V.cols();
    // V row-major (n x r) for contiguous rows; sign s_k = +1 (k < npos) or -1
    std::vector<double> Vr((size_t)n * r);
    for (int l = 0; l < n; ++l)
        for (int k = 0; k < r; ++k) Vr[(size_t)l * r + k] = V(l, k);
    // K_mn = sum_Q sum_k s_k X^Q_mk X^Q_nk with X^Q = B^Q V. Plain loops on purpose: the build links
    // whatever BLAS the system has (the reference netlib BLAS on the development machine, single-threaded
    // and ~10x slower than these loops); per thread a private lower triangle of K, Q statically
    // distributed (deterministic for a fixed thread count), summed afterwards.
    std::vector<std::vector<double>> Kt(threads);
#pragma omp parallel num_threads(threads)
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        std::vector<double>& Kl = Kt[tid];
        Kl.assign((size_t)n * n, 0.0);
        std::vector<double> X((size_t)n * r), Xs((size_t)n * r);
#pragma omp for schedule(static)
        for (int Q = 0; Q < m_nfit; ++Q) {
            const double* B = m_T.data() + (size_t)Q * n * n;
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
                double* sm = &Xs[(size_t)m * r];
                for (int k = 0; k < r; ++k) sm[k] = k < npos ? xm[k] : -xm[k];
            }
            for (int m = 0; m < n; ++m) {
                const double* sm = &Xs[(size_t)m * r];
                double* km = &Kl[(size_t)m * n];
                for (int q = 0; q <= m; ++q) {
                    const double* xq = &X[(size_t)q * r];
                    double s = 0.0;
                    for (int k = 0; k < r; ++k) s += sm[k] * xq[k];
                    km[q] += s;
                }
            }
        }
    }
    Matrix Kf = Matrix::Zero(n, n);
    for (int t = 0; t < threads; ++t)
        if (!Kt[t].empty())
            for (int m = 0; m < n; ++m)
            for (int q = 0; q <= m; ++q) Kf(m, q) += Kt[t][(size_t)m * n + q];
    for (int m = 0; m < n; ++m)
        for (int q = 0; q < m; ++q) Kf(q, m) = Kf(m, q);
    *K = Kf;
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
    if (config.contains("ri_metric_cutoff") && config["ri_metric_cutoff"].is_number())
        m_ri_metric_cutoff = config["ri_metric_cutoff"].get<double>();
    if (config.contains("ri_max_memory_mb") && config["ri_max_memory_mb"].is_number())
        m_ri_max_memory_mb = config["ri_max_memory_mb"].get<double>();
}

void QMEngine::resetRI()
{
    m_ri.reset();
    m_ri_failed = false;
    m_gpu_ri_dev.reset();
    m_gpu_ri_checked = false;
}

// One line per SCF at verbosity 2: number and time of the RI J/K builds (diagnostic).
void QMEngine::riReport() const
{
    if (m_ri && CurcumaLogger::get_verbosity() >= 2 && m_ri_builds > 0)
        CurcumaLogger::info(fmt::format("QM: RI-{} J/K builds: {} in {:.1f} ms ({:.2f} ms each; density factorisation {:.1f} ms; {})",
                                        m_ri_mode == "j" ? "J" : "JK", m_ri_builds, m_ri_build_ms, m_ri_build_ms / m_ri_builds,
                                        m_ri_factor_ms + m_ri->takeFactorMs(),
                                        m_gpu_ri_dev ? "contractions on " + m_gpu_ri_dev->name() : std::string("CPU")));
    m_ri_build_ms = 0.0;
    m_ri_factor_ms = 0.0;
    m_ri_builds = 0;
}

bool QMEngine::riActive() const { return m_ri_mode == "j" || m_ri_mode == "jk"; }

bool QMEngine::ensureRI() const
{
    if (m_ri) return true;
    if (m_ri_failed) return false;
    const auto t0 = std::chrono::steady_clock::now();
    std::string why;
    try {
        const std::string file = qmri::resolveAuxBasisFile(m_ri_aux, m_ri_mode);
        const qmri::AuxBasisMap aux = qmri::parseTurbomoleBasis(file);
        Matrix pos = Matrix::Zero(m_atomcount, 3);
        for (int i = 0; i < m_atomcount; ++i)
            for (int k = 0; k < 3; ++k) pos(i, k) = CurcumaUnit::Length::angstrom_to_bohr(m_geometry(i, k));
        const std::vector<qmri::Shell> sh = qmri::auxiliaryShells(aux, m_atoms, pos);
        qmri::RIOptions opt;
        opt.mode = m_ri_mode;
        opt.metric = m_ri_metric;
        opt.metric_cutoff = m_ri_metric_cutoff;
        opt.screening = m_eri_screening;
        opt.max_memory_mb = m_ri_max_memory_mb;
        opt.spherical_aux = !m_cartesian_d;
        opt.threads = m_threads;
        auto ri = std::make_shared<qmri::RIJK>();
        if (ri->setup(m_gto_basis, m_Q.size() != 0 ? &m_Q : nullptr, sh, opt, why)) {
            m_ri = ri;
            m_ri_file = file;
        }
    } catch (const std::exception& e) {
        why = e.what();
    }
    if (!m_ri) {
        m_ri_failed = true;
        CurcumaLogger::error("QM: density fitting (-qm.ri " + m_ri_mode + ") not possible: " + why);
        return false;
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
    } else if (!gpuRIBuild(P, J, &K)) {
        m_ri->build(P, J, &K);
    }
    m_ri_build_ms += msGlue(t0);
    ++m_ri_builds;
}

// -qm.gpu_ri auto|on|off: the RI contractions (J from the tensor, K as batched GEMM + SYRK) on the device.
// Any device error: counted fallback "qm_ri", CPU from then on.
bool QMEngine::gpuRIBuild(const Matrix& P, Matrix& J, Matrix* K) const
{
    if (m_gpu_ri_failed || m_gpu_ri == "off" || m_gpu_ri == "false" || !gpuRequested()) return false;
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
        const std::string backend = (m_gpu == "auto") ? gpu_plugin::firstAvailable() : m_gpu;
        if (backend == "none" || !gpu_plugin::available(backend)) return fail("no GPU plugin for -gpu " + m_gpu);
        std::string why;
        const auto t0 = std::chrono::steady_clock::now();
        const Matrix& T = m_ri->tensor();
        m_gpu_ri_dev = gpu_plugin::createQmRI(backend, T.data(), (int)T.rows(), m_ri->n(), m_gpu_device, why);
        if (!m_gpu_ri_dev) return fail("device RI setup not possible: " + why);
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info(fmt::format("QM: RI contractions on {} ({:.1f} MB tensor uploaded in {:.1f} ms)", m_gpu_ri_dev->name(),
                                            m_ri->tensorMB(), msGlue(t0)));
    }
    const int n = m_ri->n();
    const Matrix Pc = P;
    J = Matrix::Zero(n, n);
    // Coulomb: g = T . vec(P) on the device; mode j needs c = J^-1 g on the host (N^2), then J = T^T c
    std::vector<double> g((size_t)m_ri->nfit());
    if (!m_gpu_ri_dev->contract(Pc.data(), g.data())) return fail("device g = T vec(P) failed");
    Eigen::VectorXd c = m_ri->fitCoefficients(Eigen::Map<const Eigen::VectorXd>(g.data(), (Eigen::Index)g.size()));
    if (!m_gpu_ri_dev->coulomb(c.data(), J.data())) return fail("device J = T^T c failed");
    if (!K) return true;
    Eigen::MatrixXd V;
    int npos = 0;
    const auto tf = std::chrono::steady_clock::now();
    qmri::RIJK::factorDensity(P, V, npos);
    m_ri_factor_ms += msGlue(tf);
    Matrix Kd = Matrix::Zero(n, n);
    if (!m_gpu_ri_dev->exchange(V.data(), (int)V.cols(), npos, Kd.data())) return fail("device exchange build failed");
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
