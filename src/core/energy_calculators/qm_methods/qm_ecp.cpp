/*
 * <Native QM: semi-local effective core potentials (def2-ECP)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Nov 2026). AI-generated, machine-tested only. Formulas and conventions: qm_ecp.h.
 *
 * Layout of this file:
 *   1. parser of the Turbomole '$ecp' block
 *   2. building blocks: scaled modified spherical Bessel functions, Legendre polynomials, real spherical harmonics,
 *      Lebedev angular tables T_ijk,lambda,lm, radial Gauss-Chebyshev grid
 *   3. "functions" = lists of cartesian primitives on one centre (a basis function or its centre derivative)
 *   4. local part (three-centre overlaps, analytic)
 *   5. semi-local part (angular projections on the radial grid, quadrature-order doubling)
 *   6. public entry points buildECP / gradientECP
 */

#include "qm_ecp.h"

#include "qm_grid.h"

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <map>
#include <sstream>
#include <stdexcept>

namespace qmecp {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kExpCut = 46.0;  // e^-46 ~ 1e-20: terms below are dropped

std::string lower(std::string s)
{
    for (char& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

int letterL(char c)
{
    static const std::string letters = "spdfghik";
    const auto p = letters.find((char)std::tolower((unsigned char)c));
    return p == std::string::npos ? -1 : (int)p;
}
}  // namespace

// =====================================================================================================
// 1. Parser
// =====================================================================================================

bool parseTurbomoleECP(const std::string& block, ElementECP& out, std::string& why)
{
    out = ElementECP();
    std::istringstream in(block);
    std::string line;
    bool have_header = false;
    int channel = -2;  // -1 = local, l >= 0 = semi-local, -2 = none yet
    std::map<int, std::vector<EcpTerm>> semi;
    while (std::getline(in, line)) {
        const auto f = line.find_first_not_of(" \t\r");
        if (f == std::string::npos || line[f] == '#' || line[f] == '*' || line[f] == '$') continue;
        std::string t = lower(line.substr(f));
        if (t.rfind("ncore", 0) == 0) {
            // "ncore = 28 lmax = 3"
            std::string s = t;
            std::replace(s.begin(), s.end(), '=', ' ');
            std::istringstream ss(s);
            std::string k1, k2;
            if (!(ss >> k1 >> out.ncore >> k2 >> out.lmax) || k2 != "lmax" || out.ncore < 0 || out.lmax < 0 ||
                out.lmax > 6) {
                why = "ECP: cannot read '" + line + "'";
                return false;
            }
            have_header = true;
            continue;
        }
        if (!have_header) continue;  // element name line ("rb def2-ecp")
        std::istringstream ss(t);
        std::vector<std::string> tok;
        for (std::string x; ss >> x;) tok.push_back(x);
        if (tok.size() == 1 && std::isalpha((unsigned char)tok[0][0])) {
            const std::string& h = tok[0];
            if (h.size() == 1) {
                if (letterL(h[0]) != out.lmax) {
                    why = "ECP: local channel '" + h + "' is not lmax = " + std::to_string(out.lmax);
                    return false;
                }
                channel = -1;
            } else if (h.size() == 3 && h[1] == '-' && letterL(h[2]) == out.lmax && letterL(h[0]) >= 0 &&
                       letterL(h[0]) < out.lmax) {
                channel = letterL(h[0]);
            } else {
                why = "ECP: unknown channel header '" + h + "'";
                return false;
            }
            continue;
        }
        if (tok.size() != 3 || channel == -2) {
            why = "ECP: cannot read term line '" + line + "'";
            return false;
        }
        EcpTerm term;
        try {
            term.coef = std::stod(tok[0]);
            term.n = std::stoi(tok[1]);
            term.zeta = std::stod(tok[2]);
        } catch (...) {
            why = "ECP: cannot read term line '" + line + "'";
            return false;
        }
        if (term.n != 2) {
            why = "ECP: r^(n-2) terms with n != 2 are not implemented (line '" + line + "')";
            return false;
        }
        if (!(term.zeta > 0.0)) {
            why = "ECP: non-positive exponent in '" + line + "'";
            return false;
        }
        if (term.coef == 0.0) continue;  // e.g. the empty local h channel of the lanthanides
        if (channel == -1)
            out.local.push_back(term);
        else
            semi[channel].push_back(term);
    }
    if (!have_header) {
        why = "ECP: no 'ncore = .. lmax = ..' line";
        return false;
    }
    out.semilocal.assign(out.lmax, {});
    for (auto& [l, v] : semi) out.semilocal[l] = v;
    return true;
}

// =====================================================================================================
// 2. Building blocks
// =====================================================================================================

namespace {

// e^(-x) i_l(x), l = 0..lmax (modified spherical Bessel function of the first kind, scaled). Series for x <= 50
// (positive terms, no cancellation), the exact finite form above:
//   i_l(x) = 1/(2x) [ e^x sum_k (-1)^k a_k (2x)^-k + (-1)^(l+1) e^-x sum_k a_k (2x)^-k ],  a_k = (l+k)!/(k!(l-k)!)
void scaledBessel(double x, int lmax, double* out)
{
    if (x < 1.0e-14) {
        out[0] = 1.0;
        for (int l = 1; l <= lmax; ++l) out[l] = 0.0;
        return;
    }
    if (x <= 50.0) {
        const double ex = std::exp(-x), h = 0.5 * x * x;
        double xl = 1.0, df = 1.0;  // x^l, (2l+1)!!
        for (int l = 0; l <= lmax; ++l) {
            if (l > 0) {
                xl *= x;
                df *= (2 * l + 1);
            }
            double t = xl / df, s = t;
            for (int k = 0; k < 2000; ++k) {
                t *= h / ((k + 1.0) * (2.0 * l + 2.0 * k + 3.0));
                s += t;
                if (t < 1.0e-17 * s) break;
            }
            out[l] = ex * s;
        }
        return;
    }
    const double e2 = std::exp(-2.0 * x), inv = 1.0 / (2.0 * x);
    for (int l = 0; l <= lmax; ++l) {
        double s1 = 0.0, s2 = 0.0, a = 1.0, p = 1.0;  // a_k, (2x)^-k
        for (int k = 0; k <= l; ++k) {
            if (k > 0) {
                a *= (double)(l + k) * (l - k + 1) / k;
                p *= inv;
            }
            s1 += ((k & 1) ? -a : a) * p;
            s2 += a * p;
        }
        out[l] = inv * (s1 + ((l & 1) ? 1.0 : -1.0) * e2 * s2);
    }
}

void legendre(double t, int lmax, double* P)
{
    P[0] = 1.0;
    if (lmax >= 1) P[1] = t;
    for (int l = 2; l <= lmax; ++l) P[l] = ((2 * l - 1) * t * P[l - 1] - (l - 1) * P[l - 2]) / l;
}

// Real orthonormal spherical harmonics Y_lm on the unit sphere, l = 0..lmax, stored at offset l*l + (m + l).
// Any orthonormal real basis of each degree-l space gives the same projector sum_m |lm><lm|.
void realSphericalHarmonics(double x, double y, double z, int lmax, double* Y)
{
    const double st2 = std::max(0.0, 1.0 - z * z), phi = std::atan2(y, x);
    // associated Legendre P_l^m(z) without Condon-Shortley phase
    std::vector<double> P((lmax + 1) * (lmax + 1), 0.0);
    auto at = [&](int l, int m) -> double& { return P[l * (lmax + 1) + m]; };
    double pmm = 1.0;
    for (int m = 0; m <= lmax; ++m) {
        if (m > 0) pmm *= (2 * m - 1) * std::sqrt(st2);
        at(m, m) = pmm;
        if (m + 1 <= lmax) at(m + 1, m) = z * (2 * m + 1) * pmm;
        for (int l = m + 2; l <= lmax; ++l) at(l, m) = ((2 * l - 1) * z * at(l - 1, m) - (l + m - 1) * at(l - 2, m)) / (l - m);
    }
    for (int l = 0; l <= lmax; ++l) {
        const double base = (2 * l + 1) / (4.0 * kPi);
        Y[l * l + l] = std::sqrt(base) * at(l, 0);
        for (int m = 1; m <= l; ++m) {
            double fr = 1.0;  // (l-m)!/(l+m)!
            for (int k = l - m + 1; k <= l + m; ++k) fr /= k;
            const double nrm = std::sqrt(2.0 * base * fr) * at(l, m);
            Y[l * l + l + m] = nrm * std::cos(m * phi);
            Y[l * l + l - m] = nrm * std::sin(m * phi);
        }
    }
}

// Monomial index of x^i y^j z^k, i + j + k <= nmax.
struct MonoIndex {
    int nmax = 0;
    std::vector<int> idx;  // [(i*(nmax+1)+j)*(nmax+1)+k]
    int count = 0;
    explicit MonoIndex(int n = 0) : nmax(n), idx((n + 1) * (n + 1) * (n + 1), -1)
    {
        for (int t = 0; t <= n; ++t)
            for (int i = t; i >= 0; --i)
                for (int j = t - i; j >= 0; --j) idx[(i * (n + 1) + j) * (n + 1) + (t - i - j)] = count++;
    }
    int operator()(int i, int j, int k) const { return idx[(i * (nmax + 1) + j) * (nmax + 1) + k]; }
};

// Smallest built-in Lebedev grid that integrates polynomials of degree `deg` exactly.
int lebedevPointsForDegree(int deg)
{
    for (int n : { 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730,
                   2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 })
        if (qmgrid::lebedevDegree(n) >= deg) return n;
    return 0;
}

// Angular tables of one (function centre A, ECP centre C) pair:
//   T[mono][lambda * nlm + (l*l + m + l)] = (2 lambda + 1) int dOmega x^i y^j z^k P_lambda(Omega . u) Y_lm(Omega)
// (u = D/|D|; on-centre only lambda = 0 with P_0 = 1), lambda <= n + l, zero unless n + lambda + l is even.
struct AngularTable {
    int nmax = 0, lproj = 0, lam_max = 0, nlm = 0;
    bool on_centre = false;
    MonoIndex mono;
    std::vector<double> T;  // mono.count x (lam_max+1) x nlm
    double& at(int m, int lam, int lm) { return T[((size_t)m * (lam_max + 1) + lam) * nlm + lm]; }
    double get(int m, int lam, int lm) const { return T[((size_t)m * (lam_max + 1) + lam) * nlm + lm]; }
};

AngularTable angularTable(int nmax, int lproj, const double D[3])
{
    AngularTable a;
    a.nmax = nmax;
    a.lproj = lproj;
    a.lam_max = nmax + lproj;
    a.nlm = lproj * lproj;  // l = 0..lproj-1
    a.mono = MonoIndex(nmax);
    const double d = std::sqrt(D[0] * D[0] + D[1] * D[1] + D[2] * D[2]);
    a.on_centre = d < 1.0e-12;
    a.T.assign((size_t)a.mono.count * (a.lam_max + 1) * std::max(1, a.nlm), 0.0);
    if (lproj == 0) return a;
    const int npts = lebedevPointsForDegree(2 * (nmax + lproj));
    if (npts == 0) throw std::runtime_error("ECP: no Lebedev grid of degree " + std::to_string(2 * (nmax + lproj)));
    Matrix xyz;
    Vector w;
    qmgrid::lebedevGrid(npts, xyz, w);
    const double u[3] = { a.on_centre ? 0.0 : D[0] / d, a.on_centre ? 0.0 : D[1] / d, a.on_centre ? 0.0 : D[2] / d };
    std::vector<double> Y(a.nlm > 0 ? (size_t)lproj * lproj : 1), P(a.lam_max + 1);
    std::vector<double> px(nmax + 1), py(nmax + 1), pz(nmax + 1);
    for (int q = 0; q < npts; ++q) {
        const double x = xyz(q, 0), y = xyz(q, 1), z = xyz(q, 2), wq = 4.0 * kPi * w(q);
        realSphericalHarmonics(x, y, z, lproj - 1, Y.data());
        legendre(x * u[0] + y * u[1] + z * u[2], a.lam_max, P.data());
        px[0] = py[0] = pz[0] = 1.0;
        for (int k = 1; k <= nmax; ++k) {
            px[k] = px[k - 1] * x;
            py[k] = py[k - 1] * y;
            pz[k] = pz[k - 1] * z;
        }
        for (int t = 0; t <= nmax; ++t)
            for (int i = t; i >= 0; --i)
                for (int j = t - i; j >= 0; --j) {
                    const int k = t - i - j, mi = a.mono(i, j, k);
                    const double mon = wq * px[i] * py[j] * pz[k];
                    for (int l = 0; l < lproj; ++l)
                        for (int lam = 0; lam <= t + l; ++lam) {
                            if ((t + lam + l) & 1) continue;
                            if (a.on_centre && lam > 0) break;
                            const double f = (2 * lam + 1) * mon * (a.on_centre ? 1.0 : P[lam]);
                            for (int m = 0; m < 2 * l + 1; ++m) a.at(mi, lam, l * l + m) += f * Y[l * l + m];
                        }
                }
    }
    return a;
}

// Gauss-Chebyshev (2nd kind) radial grid with n points on r = ln(2/(1-x))/ln 2; w includes dr/dx (not r^2).
void radialGrid(int n, std::vector<double>& r, std::vector<double>& w)
{
    r.resize(n);
    w.resize(n);
    const double ln2 = std::log(2.0);
    for (int i = 1; i <= n; ++i) {
        const double th = i * kPi / (n + 1), x = std::cos(th);
        r[i - 1] = std::log(2.0 / (1.0 - x)) / ln2;
        w[i - 1] = kPi / (n + 1) * std::sin(th) / ((1.0 - x) * ln2);
    }
}

// =====================================================================================================
// 3. Functions: sum of cartesian primitives c (x-Ax)^a (y-Ay)^b (z-Az)^c exp(-alpha |r-A|^2) on one centre
// =====================================================================================================

struct Term {
    double c, alpha;
    int p[3];
};
struct Fn {
    int atom;
    double A[3];
    std::vector<Term> terms;
    int maxPower() const
    {
        int m = 0;
        for (const auto& t : terms) m = std::max(m, t.p[0] + t.p[1] + t.p[2]);
        return m;
    }
};

Fn basisFunction(const GTO::Orbital& o)
{
    Fn f;
    f.atom = o.atom;
    f.A[0] = o.x;
    f.A[1] = o.y;
    f.A[2] = o.z;
    int l, m, n;
    GTO::orbitalTypeToComponents(o.type, l, m, n);
    for (size_t k = 0; k < o.exponents.size(); ++k) f.terms.push_back({ o.coefficients[k], o.exponents[k], { l, m, n } });
    return f;
}

// d phi / d A_d = sum_k c_k [2 alpha_k g(p + e_d) - p_d g(p - e_d)]
Fn derivativeFunction(const Fn& f, int d)
{
    Fn g = f;
    g.terms.clear();
    for (const auto& t : f.terms) {
        Term up = t;
        up.c = 2.0 * t.alpha * t.c;
        up.p[d] += 1;
        g.terms.push_back(up);
        if (t.p[d] > 0) {
            Term dn = t;
            dn.c = -t.p[d] * t.c;
            dn.p[d] -= 1;
            g.terms.push_back(dn);
        }
    }
    return g;
}

double binom(int n, int k)
{
    double b = 1.0;
    for (int i = 1; i <= k; ++i) b = b * (n - k + i) / i;
    return b;
}

// =====================================================================================================
// 4. Local part: <a| sum_k c_k exp(-zeta_k r_C^2) |b>, analytic three-centre overlap
// =====================================================================================================

double gaussMoment(int k, double g)  // int t^k exp(-g t^2) dt over R
{
    if (k & 1) return 0.0;
    double df = 1.0;
    for (int i = k - 1; i > 1; i -= 2) df *= i;
    return df / std::pow(2.0 * g, k / 2) * std::sqrt(kPi / g);
}

double threeCentreOverlap(const Term& a, const double* A, const Term& b, const double* B, double zeta, const double* C)
{
    const double g = a.alpha + b.alpha + zeta;
    double ab2 = 0, ac2 = 0, bc2 = 0, Q[3];
    for (int d = 0; d < 3; ++d) {
        ab2 += (A[d] - B[d]) * (A[d] - B[d]);
        ac2 += (A[d] - C[d]) * (A[d] - C[d]);
        bc2 += (B[d] - C[d]) * (B[d] - C[d]);
        Q[d] = (a.alpha * A[d] + b.alpha * B[d] + zeta * C[d]) / g;
    }
    const double ex = (a.alpha * b.alpha * ab2 + a.alpha * zeta * ac2 + b.alpha * zeta * bc2) / g;
    if (ex > kExpCut) return 0.0;
    double v = std::exp(-ex);
    for (int d = 0; d < 3; ++d) {
        double s = 0.0;
        const double qa = Q[d] - A[d], qb = Q[d] - B[d];
        for (int i = 0; i <= a.p[d]; ++i)
            for (int j = 0; j <= b.p[d]; ++j) {
                if ((i + j) & 1) continue;
                s += binom(a.p[d], i) * binom(b.p[d], j) * std::pow(qa, a.p[d] - i) * std::pow(qb, b.p[d] - j) *
                     gaussMoment(i + j, g);
            }
        v *= s;
    }
    return v;
}

double localElement(const Fn& a, const Fn& b, const EcpCentre& c)
{
    const double C[3] = { c.x, c.y, c.z };
    double v = 0.0;
    for (const auto& ta : a.terms)
        for (const auto& tb : b.terms)
            for (const auto& k : c.ecp->local) v += ta.c * tb.c * k.coef * threeCentreOverlap(ta, a.A, tb, b.A, k.zeta, C);
    return v;
}

// =====================================================================================================
// 5. Semi-local part
// =====================================================================================================

// Radius beyond which every semi-local channel of the centre is below 1e-20 (|c| e^(-zeta r^2)).
double semilocalCutoff(const ElementECP& e)
{
    double rc = 0.0;
    for (const auto& ch : e.semilocal)
        for (const auto& t : ch) rc = std::max(rc, std::sqrt(std::max(0.0, std::log(std::abs(t.coef)) + kExpCut) / t.zeta));
    return rc;
}

// Projections F_f,lm(r_k) of the functions `fns` (all on atom A) onto the centre C for l < lproj: out[f] is a
// (lproj^2) x npts block, row l*l + m + l. Grid points beyond `rcut` stay zero.
void projections(const std::vector<const Fn*>& fns, const EcpCentre& c, int lproj, const std::vector<double>& r,
                 double rcut, std::vector<Matrix>& out)
{
    out.assign(fns.size(), Matrix());
    if (fns.empty()) return;
    const int npts = (int)r.size(), nlm = lproj * lproj;
    const double* A = fns[0]->A;
    const double D[3] = { A[0] - c.x, A[1] - c.y, A[2] - c.z };
    const double d = std::sqrt(D[0] * D[0] + D[1] * D[1] + D[2] * D[2]);
    int nmax = 0;
    for (const Fn* f : fns) nmax = std::max(nmax, f->maxPower());
    const AngularTable ang = angularTable(nmax, lproj, D);
    const int lam_max = ang.lam_max;

    // radial factors per exponent: R[(n*(lam_max+1)+lam)*npts + k] = r^n e^(-alpha (r-D)^2) e^(-x) i_lam(x)
    std::map<double, std::vector<double>> radial;
    auto radialFor = [&](double alpha) -> const std::vector<double>& {
        auto it = radial.find(alpha);
        if (it != radial.end()) return it->second;
        std::vector<double> R((size_t)(nmax + 1) * (lam_max + 1) * npts, 0.0);
        std::vector<double> bes(lam_max + 1);
        for (int k = 0; k < npts; ++k) {
            const double rk = r[k];
            if (rk > rcut) continue;
            const double ex = alpha * (rk - d) * (rk - d);
            if (ex > kExpCut) continue;
            scaledBessel(2.0 * alpha * d * rk, ang.on_centre ? 0 : lam_max, bes.data());
            const double e = std::exp(-ex);
            double rn = 1.0;
            for (int n = 0; n <= nmax; ++n) {
                for (int lam = 0; lam <= (ang.on_centre ? 0 : lam_max); ++lam)
                    R[((size_t)n * (lam_max + 1) + lam) * npts + k] = rn * e * bes[lam];
                rn *= rk;
            }
        }
        return radial.emplace(alpha, std::move(R)).first->second;
    };

    for (size_t fi = 0; fi < fns.size(); ++fi) {
        Matrix F = Matrix::Zero(nlm, npts);
        for (const Term& t : fns[fi]->terms) {
            // whole term negligible inside rcut?
            if (d > rcut && t.alpha * (d - rcut) * (d - rcut) > kExpCut) continue;
            const std::vector<double>& R = radialFor(t.alpha);
            // angular coefficients of this term: Acoef[(n*(lam_max+1)+lam)*nlm + lm]
            std::vector<double> Ac((size_t)(nmax + 1) * (lam_max + 1) * nlm, 0.0);
            for (int i = 0; i <= t.p[0]; ++i)
                for (int j = 0; j <= t.p[1]; ++j)
                    for (int k = 0; k <= t.p[2]; ++k) {
                        const double B = binom(t.p[0], i) * std::pow(-D[0], t.p[0] - i) * binom(t.p[1], j) *
                                         std::pow(-D[1], t.p[1] - j) * binom(t.p[2], k) * std::pow(-D[2], t.p[2] - k);
                        if (B == 0.0) continue;
                        const int n = i + j + k, mi = ang.mono(i, j, k);
                        for (int lam = 0; lam <= lam_max; ++lam)
                            for (int lm = 0; lm < nlm; ++lm) Ac[((size_t)n * (lam_max + 1) + lam) * nlm + lm] += B * ang.get(mi, lam, lm);
                    }
            for (int n = 0; n <= nmax; ++n)
                for (int lam = 0; lam <= lam_max; ++lam) {
                    const double* Rr = &R[((size_t)n * (lam_max + 1) + lam) * npts];
                    for (int lm = 0; lm < nlm; ++lm) {
                        const double a = t.c * Ac[((size_t)n * (lam_max + 1) + lam) * nlm + lm];
                        if (a == 0.0) continue;
                        double* row = &F(lm, 0);
                        for (int k = 0; k < npts; ++k) row[k] += a * Rr[k];
                    }
                }
        }
        out[fi] = std::move(F);
    }
}

// Semi-local contribution of one centre at quadrature order n: V (ncart x ncart) and, if dfn != nullptr, the
// gradient pieces G[d](a) = sum_b P_ab <d_d a| P dU P |b> (d = 0..2).
struct SemiLocalResult {
    Matrix V;
    Matrix G;  // ncart x 3
};

SemiLocalResult semilocal(const std::vector<Fn>& fns, const std::vector<std::array<Fn, 3>>* dfns, const Matrix* P,
                          const EcpCentre& c, int npts, int threads)
{
    const ElementECP& e = *c.ecp;
    const int nf = (int)fns.size(), lproj = e.lmax;
    SemiLocalResult res;
    res.V = Matrix::Zero(nf, nf);
    if (dfns) res.G = Matrix::Zero(nf, 3);
    bool any = false;
    for (const auto& ch : e.semilocal) any = any || !ch.empty();
    if (!any) return res;
    std::vector<double> r, w;
    radialGrid(npts, r, w);
    const double rcut = semilocalCutoff(e);

    // group functions by atom (one angular table and radial cache per atom)
    std::map<int, std::vector<int>> by_atom;
    for (int i = 0; i < nf; ++i) by_atom[fns[i].atom].push_back(i);
    std::vector<int> atoms;
    for (auto& [a, v] : by_atom) atoms.push_back(a);

    std::vector<Matrix> F(nf);                         // nlm x npts per function
    std::vector<std::array<Matrix, 3>> dF(dfns ? nf : 0);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
    for (int ia = 0; ia < (int)atoms.size(); ++ia) {
        const auto& idx = by_atom.at(atoms[ia]);
        std::vector<const Fn*> list;
        for (int i : idx) list.push_back(&fns[i]);
        if (dfns)
            for (int i : idx)
                for (int d = 0; d < 3; ++d) list.push_back(&(*dfns)[i][d]);
        std::vector<Matrix> out;
        projections(list, c, lproj, r, rcut, out);
        for (size_t k = 0; k < idx.size(); ++k) F[idx[k]] = std::move(out[k]);
        if (dfns)
            for (size_t k = 0; k < idx.size(); ++k)
                for (int d = 0; d < 3; ++d) dF[idx[k]][d] = std::move(out[idx.size() + 3 * k + d]);
    }

    for (int l = 0; l < lproj; ++l) {
        if (e.semilocal[l].empty()) continue;
        const int ml = 2 * l + 1, M = ml * npts;
        std::vector<double> om(npts);
        for (int k = 0; k < npts; ++k) {
            double u = 0.0;
            for (const auto& t : e.semilocal[l]) u += t.coef * std::exp(-t.zeta * r[k] * r[k]);
            om[k] = w[k] * r[k] * r[k] * u;
        }
        Eigen::MatrixXd Pi(nf, M), PiW(nf, M);
        for (int i = 0; i < nf; ++i)
            for (int m = 0; m < ml; ++m)
                for (int k = 0; k < npts; ++k) {
                    const double v = F[i](l * l + m, k);
                    Pi(i, m * npts + k) = v;
                    PiW(i, m * npts + k) = v * om[k];
                }
        const Eigen::MatrixXd Vl = PiW * Pi.transpose();
        res.V += Vl;
        if (dfns) {
            // G[d](a) = sum_b P_ab sum_col dPi_d(a,col) om(col) Pi(b,col) = sum_col dPi_d(a,col) (P PiW)(a,col)
            const Eigen::MatrixXd X = Eigen::MatrixXd(*P) * PiW;
            for (int i = 0; i < nf; ++i)
                for (int d = 0; d < 3; ++d) {
                    double s = 0.0;
                    for (int m = 0; m < ml; ++m)
                        for (int k = 0; k < npts; ++k) s += dF[i][d](l * l + m, k) * X(i, m * npts + k);
                    res.G(i, d) += s;
                }
        }
    }
    return res;
}

}  // namespace

// =====================================================================================================
// 6. Public entry points
// =====================================================================================================

bool buildECP(const std::vector<GTO::Orbital>& basis, const std::vector<EcpCentre>& centres, const EcpOptions& opt,
              Matrix& V, EcpLevels& levels, std::string& why)
{
    const int n = (int)basis.size();
    V = Matrix::Zero(n, n);
    levels.level.assign(centres.size(), 0);
    std::vector<Fn> fns(n);
    for (int i = 0; i < n; ++i) fns[i] = basisFunction(basis[i]);
    for (size_t ci = 0; ci < centres.size(); ++ci) {
        const EcpCentre& c = centres[ci];
        // local part
        if (!c.ecp->local.empty()) {
#pragma omp parallel for schedule(dynamic, 1) num_threads(opt.threads > 0 ? opt.threads : 1)
            for (int i = 0; i < n; ++i)
                for (int j = 0; j <= i; ++j) {
                    const double v = localElement(fns[i], fns[j], c);
                    V(i, j) += v;
                    if (j != i) V(j, i) += v;
                }
        }
        // semi-local part, quadrature order doubled until converged
        Matrix prev;
        bool ok = false;
        for (int lev = opt.min_level; lev <= opt.max_level; ++lev) {
            const int npts = (1 << lev) - 1;
            Matrix cur = semilocal(fns, nullptr, nullptr, c, npts, opt.threads).V;
            if (std::getenv("CURCUMA_QM_ECP_DEBUG") && prev.size() != 0)
                fmt::print("ECP centre {} level {} npts {} max|dV| {:.3e} max|V| {:.3e}\n", c.atom + 1, lev, npts,
                           (cur - prev).cwiseAbs().maxCoeff(), cur.cwiseAbs().maxCoeff());
            // accept when the change is below radial_tol or the double-precision floor of the block (2e-15 max|V|;
            // CeF4 def2-SVP: max|V| 528 Eh, the change stalls at 2e-13 from 4095 points on)
            const double floor_tol = std::max(opt.radial_tol, 2.0e-15 * cur.cwiseAbs().maxCoeff());
            if (prev.size() != 0 && (cur - prev).cwiseAbs().maxCoeff() < floor_tol) {
                V += cur;
                levels.level[ci] = lev;
                ok = true;
                break;
            }
            prev = std::move(cur);
        }
        if (!ok) {
            why = fmt::format("ECP: radial quadrature of atom {} not converged to {:.1e} with {} points", c.atom + 1,
                              opt.radial_tol, (1 << opt.max_level) - 1);
            return false;
        }
    }
    return true;
}

Matrix gradientECP(const std::vector<GTO::Orbital>& basis, const std::vector<EcpCentre>& centres, const Matrix& P,
                   int natoms, const EcpOptions& opt, const EcpLevels& levels)
{
    const int n = (int)basis.size();
    Matrix g = Matrix::Zero(natoms, 3);
    if (centres.empty()) return g;
    std::vector<Fn> fns(n);
    std::vector<std::array<Fn, 3>> dfns(n);
    for (int i = 0; i < n; ++i) {
        fns[i] = basisFunction(basis[i]);
        for (int d = 0; d < 3; ++d) dfns[i][d] = derivativeFunction(fns[i], d);
    }
    for (size_t ci = 0; ci < centres.size(); ++ci) {
        const EcpCentre& c = centres[ci];
        Matrix G = Matrix::Zero(n, 3);  // G(a,d) = sum_b P_ab <d_d a|U_C|b>
        if (!c.ecp->local.empty()) {
#pragma omp parallel for schedule(dynamic, 1) num_threads(opt.threads > 0 ? opt.threads : 1)
            for (int i = 0; i < n; ++i)
                for (int d = 0; d < 3; ++d) {
                    double s = 0.0;
                    for (int j = 0; j < n; ++j)
                        if (P(i, j) != 0.0) s += P(i, j) * localElement(dfns[i][d], fns[j], c);
                    G(i, d) += s;
                }
        }
        int lev = ci < levels.level.size() ? levels.level[ci] : 0;
        if (lev <= 0) {
            EcpLevels tmp;
            Matrix Vtmp;
            std::string why;
            buildECP(basis, { c }, opt, Vtmp, tmp, why);
            lev = tmp.level.empty() || tmp.level[0] <= 0 ? opt.max_level : tmp.level[0];
        }
        G += semilocal(fns, &dfns, &P, c, (1 << lev) - 1, opt.threads).G;
        // dE/dA(a) = 2 G(a); the centre C takes -2 G(a) (translational invariance of each integral)
        for (int i = 0; i < n; ++i)
            for (int d = 0; d < 3; ++d) {
                g(fns[i].atom, d) += 2.0 * G(i, d);
                g(c.atom, d) -= 2.0 * G(i, d);
            }
    }
    return g;
}

}  // namespace qmecp
