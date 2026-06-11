/*
 * <Domain-decomposition ddCOSMO core for the native CPCM solvation model>
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
 * Claude Generated (June 2026): C++ port of tblite cpcm_dd.f90 (ddCOSMO core).
 */

#include "src/core/solvation/cpcm_dd.h"

#include <cmath>
#include <stdexcept>

namespace Curcuma {
namespace Solvation {

namespace {

constexpr double pi = 3.14159265358979323846;
constexpr double fourpi = 4.0 * pi;
const double sq2 = std::sqrt(2.0);

// Root-mean-square and max norm of a flat vector (tblite rmsvec).
void rmsvec(int n, const double* v, double& vrms, double& vmax)
{
    vrms = 0.0;
    vmax = 0.0;
    for (int i = 0; i < n; ++i) {
        vmax = std::max(vmax, std::abs(v[i]));
        vrms += v[i] * v[i];
    }
    vrms = std::sqrt(vrms / static_cast<double>(n));
}

// Gauss-Jordan elimination with full pivoting (tblite gjinv). Inverts the
// action of `a` (n x n) onto the rhs columns of `b` (n x nrhs) in place.
bool gjinv(int n, int nrhs, Eigen::MatrixXd& a, Eigen::MatrixXd& b)
{
    std::vector<int> indxc(n), indxr(n), piv(n, 0);
    int irow = 0, icol = 0;
    for (int i = 0; i < n; ++i) {
        double big = 0.0;
        for (int j = 0; j < n; ++j)
            if (piv[j] != 1)
                for (int k = 0; k < n; ++k)
                    if (piv[k] == 0)
                        if (std::abs(a(j, k)) > big) {
                            big = std::abs(a(j, k));
                            irow = j;
                            icol = k;
                        }
        piv[icol] += 1;
        if (piv[icol] > 1)
            return false;
        if (irow != icol) {
            for (int c = 0; c < n; ++c)
                std::swap(a(irow, c), a(icol, c));
            for (int c = 0; c < nrhs; ++c)
                std::swap(b(irow, c), b(icol, c));
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a(icol, icol) == 0.0)
            return false;
        const double pinv = 1.0 / a(icol, icol);
        a(icol, icol) = 1.0;
        for (int c = 0; c < n; ++c)
            a(icol, c) *= pinv;
        for (int c = 0; c < nrhs; ++c)
            b(icol, c) *= pinv;
        for (int j = 0; j < n; ++j)
            if (j != icol) {
                const double dum = a(j, icol);
                a(j, icol) = 0.0;
                for (int c = 0; c < n; ++c)
                    a(j, c) -= a(icol, c) * dum;
                for (int c = 0; c < nrhs; ++c)
                    b(j, c) -= b(icol, c) * dum;
            }
    }
    for (int j = n - 1; j >= 0; --j)
        if (indxr[j] != indxc[j])
            for (int r = 0; r < n; ++r)
                std::swap(a(r, indxr[j]), a(r, indxc[j]));
    return true;
}

// Build the DIIS B matrix (tblite makeb). nmat is 1-based.
void makeb(int /*n*/, int nmat, int /*ndiis*/, const Eigen::MatrixXd& e, Eigen::MatrixXd& b)
{
    if (nmat == 1) {
        b(0, 0) = 0.0;
        b(0, 1) = 1.0;
        b(1, 0) = 1.0;
        b(1, 1) = e.col(0).dot(e.col(0));
    } else {
        b(nmat, 0) = 1.0;
        b(0, nmat) = 1.0;
        for (int i = 1; i <= nmat - 1; ++i) {
            const double bij = e.col(i - 1).dot(e.col(nmat - 1));
            b(nmat, i) = bij;
            b(i, nmat) = bij;
        }
        b(nmat, nmat) = e.col(nmat - 1).dot(e.col(nmat - 1));
    }
}

// DIIS extrapolation step (tblite diis). nmat is 1-based and updated in place.
void diis(int n, int& nmat, int ndiis, Eigen::MatrixXd& x, Eigen::MatrixXd& e,
          Eigen::MatrixXd& b, double* xnew)
{
    if (nmat >= ndiis) {
        for (int j = 2; j <= nmat - 10; ++j)
            for (int k = 2; k <= nmat - 10; ++k)
                b(j - 1, k - 1) = b(j + 10 - 1, k + 10 - 1);
        for (int j = 1; j <= nmat - 10; ++j) {
            x.col(j - 1) = x.col(j + 10 - 1);
            e.col(j - 1) = e.col(j + 10 - 1);
        }
        nmat = nmat - 10;
    }
    const int nmat1 = nmat + 1;
    Eigen::MatrixXd bloc(nmat1, nmat1), cex(nmat1, 1);
    makeb(n, nmat, ndiis, e, b);
    bloc = b.topLeftCorner(nmat1, nmat1);
    cex.setZero();
    cex(0, 0) = 1.0;
    if (!gjinv(nmat1, 1, bloc, cex)) {
        nmat = 1;
        return;
    }
    Eigen::Map<Eigen::VectorXd> xn(xnew, n);
    xn.setZero();
    for (int i = 1; i <= nmat; ++i)
        xn += cex(i, 0) * x.col(i - 1);
    nmat = nmat + 1;
}

} // namespace

double DomainDecomposition::fsw(double t, double s, double eta)
{
    const double f6 = 6.0, f10 = 10.0, f12 = 12.0, f15 = 15.0;
    const double x = t - (s + 1.0) * eta / 2.0;
    const double flow = 1.0 - eta;
    if (x >= 1.0)
        return 0.0;
    if (x <= flow)
        return 1.0;
    const double a = f15 * eta - f12;
    const double b = f10 * eta * eta - f15 * eta + f6;
    return ((x - 1.0) * (x - 1.0) * (1.0 - x) * (f6 * x * x + a * x + b)) / std::pow(eta, 5);
}

double DomainDecomposition::dfsw(double t, double s, double eta)
{
    const double f30 = 30.0;
    const double x = t - (s + 1.0) * eta / 2.0;
    const double flow = 1.0 - eta;
    if (x >= 1.0)
        return 0.0;
    if (x <= flow)
        return 0.0;
    return f30 * (1.0 - x) * (x - 1.0) * (x - 1.0 + eta) * (x - 1.0 + eta) / std::pow(eta, 5);
}

void DomainDecomposition::trgev(double x, double y, double* cx, double* sx) const
{
    cx[0] = 1.0;
    sx[0] = 0.0;
    if (lmax >= 1) {
        cx[1] = x;
        sx[1] = y;
    }
    for (int m = 2; m <= lmax; ++m) {
        cx[m] = 2.0 * x * cx[m - 1] - cx[m - 2];
        sx[m] = 2.0 * x * sx[m - 1] - sx[m - 2];
    }
}

void DomainDecomposition::polleg(double x, double y, double* plm) const
{
    double fact = 1.0;
    double pmm = 1.0;
    const double somx2 = y;
    for (int m = 0; m <= lmax; ++m) {
        plm[idx(m, m)] = pmm;
        if (m == lmax)
            return;
        const double fm = static_cast<double>(m);
        double pmm1 = x * (2.0 * fm + 1.0) * pmm;
        plm[idx(m + 1, m)] = pmm1;
        const double pmmo = pmm;
        double pmm_loc = pmm;
        for (int l = m + 2; l <= lmax; ++l) {
            const double fl = static_cast<double>(l);
            const double pll = (x * (2.0 * fl - 1.0) * pmm1 - (fl + fm - 1.0) * pmm_loc) / (fl - fm);
            plm[idx(l, m)] = pll;
            pmm_loc = pmm1;
            pmm1 = pll;
        }
        pmm = -pmmo * fact * somx2;
        fact = fact + 2.0;
    }
}

void DomainDecomposition::ylmbas(const double x[3], double* basloc, double* vplm,
                                 double* vcos, double* vsin) const
{
    const double cthe = x[2];
    const double sthe = std::sqrt(1.0 - cthe * cthe);
    double cphi, sphi;
    if (sthe != 0.0) {
        cphi = x[0] / sthe;
        sphi = x[1] / sthe;
    } else {
        cphi = 0.0;
        sphi = 0.0;
    }
    if (sthe != 0.0) {
        trgev(cphi, sphi, vcos, vsin);
    } else {
        for (int m = 0; m <= lmax; ++m) {
            vcos[m] = 1.0;
            vsin[m] = 0.0;
        }
    }
    polleg(cthe, sthe, vplm);
    for (int l = 0; l <= lmax; ++l) {
        const int ind = idx(l, 0);
        basloc[ind] = m_facs[ind] * vplm[ind];
        for (int m = 1; m <= l; ++m) {
            const double plm = vplm[ind + m];
            basloc[ind + m] = m_facs[ind + m] * plm * vcos[m];
            basloc[ind - m] = m_facs[ind - m] * plm * vsin[m];
        }
    }
}

void DomainDecomposition::dbasis(const double x[3], double* basloc, double* dbsloc,
                                 double* vplm, double* vcos, double* vsin) const
{
    const double cthe = x[2];
    const double sthe = std::sqrt(1.0 - cthe * cthe);
    double cphi, sphi;
    if (sthe != 0.0) {
        cphi = x[0] / sthe;
        sphi = x[1] / sthe;
    } else {
        cphi = 1.0;
        sphi = 0.0;
    }
    double et[3] = {cthe * cphi, cthe * sphi, -sthe};
    double ep[3];
    if (sthe != 0.0) {
        ep[0] = -sphi / sthe;
        ep[1] = cphi / sthe;
        ep[2] = 0.0;
    } else {
        ep[0] = 0.0;
        ep[1] = 1.0;
        ep[2] = 0.0;
    }
    if (sthe != 0.0) {
        trgev(cphi, sphi, vcos, vsin);
    } else {
        for (int m = 0; m <= lmax; ++m) {
            vcos[m] = 1.0;
            vsin[m] = 0.0;
        }
    }
    const double VC = 0.0;
    const double VS = cthe;
    polleg(cthe, sthe, vplm);

    for (int k = 0; k < nylm; ++k) {
        basloc[k] = 0.0;
        dbsloc[0 + 3 * k] = 0.0;
        dbsloc[1 + 3 * k] = 0.0;
        dbsloc[2 + 3 * k] = 0.0;
    }
    for (int l = 0; l <= lmax; ++l) {
        const int ind = idx(l, 0);
        double fln = m_facs[ind];
        basloc[ind] = fln * vplm[ind];
        if (l > 0) {
            for (int c = 0; c < 3; ++c)
                dbsloc[c + 3 * ind] = fln * vplm[ind + 1] * et[c];
        }
        for (int m = 1; m <= l; ++m) {
            fln = m_facs[ind + m];
            const double plm = fln * vplm[ind + m];
            double pp1 = 0.0;
            if (m < l)
                pp1 = -0.5 * vplm[ind + m + 1];
            const double pm1 = 0.5 * (static_cast<double>(l + m) * static_cast<double>(l - m + 1) * vplm[ind + m - 1]);
            const double pp = pp1 + pm1;

            basloc[ind + m] = plm * vcos[m];
            if (sthe != 0.0) {
                for (int c = 0; c < 3; ++c)
                    dbsloc[c + 3 * (ind + m)] = -fln * pp * vcos[m] * et[c] - static_cast<double>(m) * plm * vsin[m] * ep[c];
            } else {
                for (int c = 0; c < 3; ++c)
                    dbsloc[c + 3 * (ind + m)] = -fln * pp * vcos[m] * et[c] - fln * pp * ep[c] * VC;
            }

            basloc[ind - m] = plm * vsin[m];
            if (sthe != 0.0) {
                for (int c = 0; c < 3; ++c)
                    dbsloc[c + 3 * (ind - m)] = -fln * pp * vsin[m] * et[c] + static_cast<double>(m) * plm * vcos[m] * ep[c];
            } else {
                for (int c = 0; c < 3; ++c)
                    dbsloc[c + 3 * (ind - m)] = -fln * pp * vsin[m] * et[c] - fln * pp * ep[c] * VS;
            }
        }
    }
}

double DomainDecomposition::intmlp(double t, const double* sigma_col, const double* basloc) const
{
    double tt = 1.0;
    double ss = 0.0;
    for (int l = 0; l <= lmax; ++l) {
        const int ind = idx(l, 0);
        const double fac = tt / m_facl[ind];
        double dot = 0.0;
        for (int k = l * l; k <= l * l + 2 * l; ++k)
            dot += basloc[k] * sigma_col[k];
        ss += fac * dot;
        tt *= t;
    }
    return ss;
}

void DomainDecomposition::init(const DomainDecompositionInput& input,
                               const std::vector<double>& rvdw,
                               const std::vector<LebedevPoint>& grid)
{
    nat = static_cast<int>(rvdw.size());
    ngrid = static_cast<int>(grid.size());
    lmax = input.lmax;
    conv = input.conv;
    eta = input.eta;
    m_grad = true;
    nylm = (lmax + 1) * (lmax + 1);

    m_rvdw = rvdw;
    m_w.assign(ngrid, 0.0);
    m_grid.resize(3, ngrid);
    m_basis.resize(nylm, ngrid);
    m_fact.assign(std::max(2 * lmax + 1, 2), 0.0);
    m_facl.assign(nylm, 0.0);
    m_facs.assign(nylm, 0.0);
    m_inl.assign(nat + 1, 0);
    m_nl.assign(static_cast<size_t>(nat) * nngmax, 0);
    m_fi.resize(ngrid, nat);
    m_ui.resize(ngrid, nat);
    m_zi.assign(static_cast<size_t>(3) * ngrid * nat, 0.0);
    m_xyz.resize(3, nat);

    // factorials k! : fact[0]=1, fact[1]=1, fact[k]=k*fact[k-1]
    m_fact[0] = 1.0;
    m_fact[1] = 1.0;
    for (int i = 2; i < static_cast<int>(m_fact.size()); ++i)
        m_fact[i] = static_cast<double>(i) * m_fact[i - 1];

    for (int l = 0; l <= lmax; ++l) {
        const int ind = idx(l, 0);
        const double fl = (2.0 * static_cast<double>(l) + 1.0) / fourpi;
        const double ffl = std::sqrt(fl);
        for (int k = l * l; k <= l * l + 2 * l; ++k)
            m_facl[k] = fl;
        m_facs[ind] = ffl;
        for (int m = 1; m <= l; ++m) {
            double fnorm = sq2 * ffl * std::sqrt(m_fact[l - m] / m_fact[l + m]);
            if (m % 2 == 1)
                fnorm = -fnorm;
            m_facs[ind + m] = fnorm;
            m_facs[ind - m] = fnorm;
        }
    }

    for (int i = 0; i < ngrid; ++i) {
        m_grid(0, i) = grid[i].x;
        m_grid(1, i) = grid[i].y;
        m_grid(2, i) = grid[i].z;
        m_w[i] = grid[i].weight * fourpi;
    }

    std::vector<double> vplm(nylm), vcos(lmax + 1), vsin(lmax + 1);
    for (int i = 0; i < ngrid; ++i) {
        const double xg[3] = {m_grid(0, i), m_grid(1, i), m_grid(2, i)};
        ylmbas(xg, m_basis.col(i).data(), vplm.data(), vcos.data(), vsin.data());
    }
}

void DomainDecomposition::update(const Mat& xyz)
{
    m_xyz = xyz;
    mknnl();
    mkfiui();
    mkccav();
}

void DomainDecomposition::mknnl()
{
    int ii = 0;
    for (int iat = 0; iat < nat; ++iat) {
        m_inl[iat] = ii;
        for (int jat = 0; jat < nat; ++jat) {
            if (iat != jat) {
                const double dx = m_xyz(0, iat) - m_xyz(0, jat);
                const double dy = m_xyz(1, iat) - m_xyz(1, jat);
                const double dz = m_xyz(2, iat) - m_xyz(2, jat);
                const double d2 = dx * dx + dy * dy + dz * dz;
                const double r = m_rvdw[iat] + m_rvdw[jat];
                if (d2 <= r * r) {
                    m_nl[ii] = jat;
                    ++ii;
                }
            }
        }
    }
    m_inl[nat] = ii;
}

void DomainDecomposition::mkfiui()
{
    m_fi.setZero();
    m_ui.setZero();
    std::fill(m_zi.begin(), m_zi.end(), 0.0);
    const double swthr = 1.0 + (se + 1.0) * eta / 2.0;

    for (int iat = 0; iat < nat; ++iat) {
        for (int i = 0; i < ngrid; ++i) {
            for (int ii = m_inl[iat]; ii < m_inl[iat + 1]; ++ii) {
                const int jat = m_nl[ii];
                const double v0 = m_xyz(0, iat) + m_rvdw[iat] * m_grid(0, i) - m_xyz(0, jat);
                const double v1 = m_xyz(1, iat) + m_rvdw[iat] * m_grid(1, i) - m_xyz(1, jat);
                const double v2 = m_xyz(2, iat) + m_rvdw[iat] * m_grid(2, i) - m_xyz(2, jat);
                const double vv = std::sqrt(v0 * v0 + v1 * v1 + v2 * v2);
                const double t = vv / m_rvdw[jat];
                const double xt = fsw(t, se, eta);
                if (m_grad && t < swthr && t > swthr - eta) {
                    const double fac = dfsw(t, se, eta) / m_rvdw[jat];
                    const size_t base = static_cast<size_t>(3) * (i + static_cast<size_t>(ngrid) * iat);
                    m_zi[base + 0] += fac * v0 / vv;
                    m_zi[base + 1] += fac * v1 / vv;
                    m_zi[base + 2] += fac * v2 / vv;
                }
                m_fi(i, iat) += xt;
            }
            if (m_fi(i, iat) <= 1.0)
                m_ui(i, iat) = 1.0 - m_fi(i, iat);
        }
    }
}

void DomainDecomposition::mkccav()
{
    ncav = 0;
    for (int iat = 0; iat < nat; ++iat)
        for (int i = 0; i < ngrid; ++i)
            if (m_ui(i, iat) > 0.0)
                ++ncav;

    m_ccav.resize(3, ncav);
    m_cavAtom.assign(ncav, 0);
    int ii = 0;
    for (int iat = 0; iat < nat; ++iat)
        for (int i = 0; i < ngrid; ++i)
            if (m_ui(i, iat) > 0.0) {
                m_ccav(0, ii) = m_xyz(0, iat) + m_rvdw[iat] * m_grid(0, i);
                m_ccav(1, ii) = m_xyz(1, iat) + m_rvdw[iat] * m_grid(1, i);
                m_ccav(2, ii) = m_xyz(2, iat) + m_rvdw[iat] * m_grid(2, i);
                m_cavAtom[ii] = iat;
                ++ii;
            }
}

void DomainDecomposition::intrhs(int /*iat*/, const double* x_grid, double* xlm) const
{
    for (int k = 0; k < nylm; ++k)
        xlm[k] = 0.0;
    for (int ig = 0; ig < ngrid; ++ig) {
        const double wx = m_w[ig] * x_grid[ig];
        const double* b = m_basis.col(ig).data();
        for (int k = 0; k < nylm; ++k)
            xlm[k] += b[k] * wx;
    }
}

void DomainDecomposition::calcv(bool first, int iat, double* pot, const Mat& sigma,
                                double* basloc, double* vplm, double* vcos, double* vsin) const
{
    for (int ig = 0; ig < ngrid; ++ig)
        pot[ig] = 0.0;
    if (first)
        return;

    for (int its = 0; its < ngrid; ++its) {
        if (m_ui(its, iat) < 1.0) {
            for (int ij = m_inl[iat]; ij < m_inl[iat + 1]; ++ij) {
                const int jat = m_nl[ij];
                const double v0 = m_xyz(0, iat) + m_rvdw[iat] * m_grid(0, its) - m_xyz(0, jat);
                const double v1 = m_xyz(1, iat) + m_rvdw[iat] * m_grid(1, its) - m_xyz(1, jat);
                const double v2 = m_xyz(2, iat) + m_rvdw[iat] * m_grid(2, its) - m_xyz(2, jat);
                const double vv = std::sqrt(v0 * v0 + v1 * v1 + v2 * v2);
                const double tij = vv / m_rvdw[jat];
                if (tij < 1.0) {
                    const double sij[3] = {v0 / vv, v1 / vv, v2 / vv};
                    const double xij = fsw(tij, se, eta);
                    const double oij = (m_fi(its, iat) > 1.0) ? xij / m_fi(its, iat) : xij;
                    ylmbas(sij, basloc, vplm, vcos, vsin);
                    pot[its] += oij * intmlp(tij, sigma.col(jat).data(), basloc);
                }
            }
        }
    }
}

void DomainDecomposition::adjrhs(int iat, const Mat& xi, double* vlm, double* basloc,
                                 double* vplm, double* vcos, double* vsin) const
{
    for (int ij = m_inl[iat]; ij < m_inl[iat + 1]; ++ij) {
        const int jat = m_nl[ij];
        for (int ig = 0; ig < ngrid; ++ig) {
            const double v0 = m_xyz(0, jat) + m_rvdw[jat] * m_grid(0, ig) - m_xyz(0, iat);
            const double v1 = m_xyz(1, jat) + m_rvdw[jat] * m_grid(1, ig) - m_xyz(1, iat);
            const double v2 = m_xyz(2, jat) + m_rvdw[jat] * m_grid(2, ig) - m_xyz(2, iat);
            const double vvji = std::sqrt(v0 * v0 + v1 * v1 + v2 * v2);
            const double tji = vvji / m_rvdw[iat];
            if (tji < (1.0 + (se + 1.0) / 2.0 * eta)) {
                const double sji[3] = {v0 / vvji, v1 / vvji, v2 / vvji};
                const double xji = fsw(tji, se, eta);
                const double oji = (m_fi(ig, jat) > 1.0) ? xji / m_fi(ig, jat) : xji;
                ylmbas(sji, basloc, vplm, vcos, vsin);
                double t = 1.0;
                const double fac = m_w[ig] * xi(ig, jat) * oji;
                for (int l = 0; l <= lmax; ++l) {
                    const int ind = idx(l, 0);
                    const double ffac = fac * t / m_facl[ind];
                    for (int m = -l; m <= l; ++m)
                        vlm[ind + m] += ffac * basloc[ind + m];
                    t *= tji;
                }
            }
        }
    }
}

void DomainDecomposition::lx(const Mat& x, Mat& y) const
{
    y.setZero();
    std::vector<double> pot(ngrid), vplm(nylm), basloc(nylm), vcos(lmax + 1), vsin(lmax + 1);
    for (int iat = 0; iat < nat; ++iat) {
        calcv(false, iat, pot.data(), x, basloc.data(), vplm.data(), vcos.data(), vsin.data());
        intrhs(iat, pot.data(), y.col(iat).data());
        y.col(iat) *= -1.0;
    }
}

void DomainDecomposition::lstarx(const Mat& x, Mat& y) const
{
    y.setZero();
    Mat xi(ngrid, nat);
    for (int iat = 0; iat < nat; ++iat)
        for (int ig = 0; ig < ngrid; ++ig)
            xi(ig, iat) = x.col(iat).dot(m_basis.col(ig));

    std::vector<double> vplm(nylm), basloc(nylm), vcos(lmax + 1), vsin(lmax + 1);
    for (int iat = 0; iat < nat; ++iat) {
        adjrhs(iat, xi, y.col(iat).data(), basloc.data(), vplm.data(), vcos.data(), vsin.data());
        y.col(iat) *= -1.0;
    }
}

void DomainDecomposition::ldm1x(const Mat& x, Mat& y) const
{
    for (int iat = 0; iat < nat; ++iat)
        for (int k = 0; k < nylm; ++k)
            y(k, iat) = m_facl[k] * x(k, iat);
}

void DomainDecomposition::hsnorm(const double* u, double& unorm) const
{
    unorm = 0.0;
    for (int l = 0; l <= lmax; ++l) {
        const int ind = idx(l, 0);
        const double fac = 1.0 / (1.0 + static_cast<double>(l));
        for (int m = -l; m <= l; ++m)
            unorm += fac * u[ind + m] * u[ind + m];
    }
    unorm = std::sqrt(unorm);
}

double DomainDecomposition::hnorm(const Mat& x) const
{
    std::vector<double> u(nat);
    for (int iat = 0; iat < nat; ++iat)
        hsnorm(x.col(iat).data(), u[iat]);
    double vrms, vmax;
    rmsvec(nat, u.data(), vrms, vmax);
    return vrms;
}

void DomainDecomposition::jacobiDiis(int n, double tol, const Mat& rhs, Mat& x, int& n_iter,
                                     bool& ok, const std::function<void(const Mat&, Mat&)>& matvec) const
{
    const int diis_max = ndiis;
    const bool dodiis = (diis_max != 0);
    const int maxiter = n_iter;
    Eigen::MatrixXd x_diis, e_diis, bmat;
    int nmat = 1;
    if (dodiis) {
        const int lenb = diis_max + 1;
        x_diis = Eigen::MatrixXd::Zero(n, diis_max);
        e_diis = Eigen::MatrixXd::Zero(n, diis_max);
        bmat = Eigen::MatrixXd::Zero(lenb, lenb);
    }
    Mat y(rhs.rows(), rhs.cols());
    Mat x_new(rhs.rows(), rhs.cols());
    ok = false;
    int it = 0;
    for (it = 1; it <= maxiter; ++it) {
        matvec(x, y);
        y = rhs - y;
        ldm1x(y, x_new);
        if (dodiis) {
            x_diis.col(nmat - 1) = Eigen::Map<const Eigen::VectorXd>(x_new.data(), n);
            e_diis.col(nmat - 1) = Eigen::Map<const Eigen::VectorXd>(x_new.data(), n) - Eigen::Map<const Eigen::VectorXd>(x.data(), n);
            diis(n, nmat, diis_max, x_diis, e_diis, bmat, x_new.data());
        }
        x = x_new - x;
        const double rms_norm = hnorm(x);
        ok = (rms_norm < tol);
        x = x_new;
        if (ok)
            break;
    }
    n_iter = (it > maxiter) ? maxiter : it;
}

void DomainDecomposition::wghpot(const std::vector<double>& phi, Mat& g) const
{
    int ic = 0;
    g.setZero();
    for (int iat = 0; iat < nat; ++iat)
        for (int ig = 0; ig < ngrid; ++ig)
            if (m_ui(ig, iat) != 0.0) {
                g(ig, iat) = -m_ui(ig, iat) * phi[ic];
                ++ic;
            }
}

void DomainDecomposition::solveCosmoDirect(bool cart, const std::vector<double>& phi,
                                           const Mat& glm, Mat& sigma, bool restart) const
{
    const double tol = conv;
    int n_iter = 200;
    bool ok = false;
    Mat rhs(nylm, nat);
    if (cart) {
        Mat g(ngrid, nat);
        wghpot(phi, g);
        std::vector<double> col(ngrid);
        for (int iat = 0; iat < nat; ++iat) {
            for (int ig = 0; ig < ngrid; ++ig)
                col[ig] = g(ig, iat);
            intrhs(iat, col.data(), rhs.col(iat).data());
        }
    } else {
        rhs = glm;
    }
    if (!restart)
        for (int iat = 0; iat < nat; ++iat)
            for (int k = 0; k < nylm; ++k)
                sigma(k, iat) = m_facl[k] * rhs(k, iat);

    jacobiDiis(nat * nylm, tol, rhs, sigma, n_iter, ok,
               [this](const Mat& a, Mat& b) { lx(a, b); });
    if (!ok)
        throw std::runtime_error("direct ddCOSMO did not converge");
}

void DomainDecomposition::solveCosmoAdjoint(const Mat& psi, Mat& sigma, bool restart,
                                            double accuracy) const
{
    const double tol = (accuracy > 0.0) ? accuracy : conv;
    int n_iter = 200;
    bool ok = false;
    if (!restart)
        for (int iat = 0; iat < nat; ++iat)
            for (int k = 0; k < nylm; ++k)
                sigma(k, iat) = m_facl[k] * psi(k, iat);

    jacobiDiis(nat * nylm, tol, psi, sigma, n_iter, ok,
               [this](const Mat& a, Mat& b) { lstarx(a, b); });
    if (!ok)
        throw std::runtime_error("adjoint ddCOSMO did not converge");
}

void DomainDecomposition::fdoka(int iat, const Mat& sigma, const double* xi_col, double* basloc,
                                double* dbsloc, double* vplm, double* vcos, double* vsin, double fx[3]) const
{
    const double tlow = 1.0 - 0.5 * (1.0 - se) * eta;
    const double thigh = 1.0 + 0.5 * (1.0 + se) * eta;
    for (int ig = 0; ig < ngrid; ++ig) {
        double va[3] = {0.0, 0.0, 0.0};
        for (int ij = m_inl[iat]; ij < m_inl[iat + 1]; ++ij) {
            const int jat = m_nl[ij];
            const double vij[3] = {
                m_xyz(0, iat) + m_rvdw[iat] * m_grid(0, ig) - m_xyz(0, jat),
                m_xyz(1, iat) + m_rvdw[iat] * m_grid(1, ig) - m_xyz(1, jat),
                m_xyz(2, iat) + m_rvdw[iat] * m_grid(2, ig) - m_xyz(2, jat)};
            const double vvij = std::sqrt(vij[0] * vij[0] + vij[1] * vij[1] + vij[2] * vij[2]);
            const double tij = vvij / m_rvdw[jat];
            if (tij >= thigh)
                continue;
            const double sij[3] = {vij[0] / vvij, vij[1] / vvij, vij[2] / vvij};
            dbasis(sij, basloc, dbsloc, vplm, vcos, vsin);
            double alp[3] = {0.0, 0.0, 0.0};
            double t = 1.0;
            const double* sig = sigma.col(jat).data();
            for (int l = 1; l <= lmax; ++l) {
                const int ind = idx(l, 0);
                const double fl = static_cast<double>(l);
                const double fac = t / m_facl[ind];
                for (int m = -l; m <= l; ++m) {
                    const double f2 = fac * sig[ind + m];
                    const double f1 = f2 * fl * basloc[ind + m];
                    for (int c = 0; c < 3; ++c)
                        alp[c] += f1 * sij[c] + f2 * dbsloc[c + 3 * (ind + m)];
                }
                t *= tij;
            }
            const double beta = intmlp(tij, sig, basloc);
            const double xij = fsw(tij, se, eta);
            double oij, f2;
            if (m_fi(ig, iat) > 1.0) {
                oij = xij / m_fi(ig, iat);
                f2 = -oij / m_fi(ig, iat);
            } else {
                oij = xij;
                f2 = 0.0;
            }
            const double f1 = oij / m_rvdw[jat];
            const size_t zbase = static_cast<size_t>(3) * (ig + static_cast<size_t>(ngrid) * iat);
            for (int c = 0; c < 3; ++c)
                va[c] += f1 * alp[c] + beta * f2 * m_zi[zbase + c];
            if (tij > tlow) {
                double f3 = beta * dfsw(tij, se, eta) / m_rvdw[jat];
                if (m_fi(ig, iat) > 1.0)
                    f3 = f3 / m_fi(ig, iat);
                for (int c = 0; c < 3; ++c)
                    va[c] += f3 * sij[c];
            }
        }
        for (int c = 0; c < 3; ++c)
            fx[c] += m_w[ig] * xi_col[ig] * va[c];
    }
}

void DomainDecomposition::fdokb(int iat, const Mat& sigma, const Mat& xi, double* basloc,
                                double* dbsloc, double* vplm, double* vcos, double* vsin, double fx[3]) const
{
    const double tlow = 1.0 - 0.5 * (1.0 - se) * eta;
    const double thigh = 1.0 + 0.5 * (1.0 + se) * eta;
    for (int ig = 0; ig < ngrid; ++ig) {
        double vb[3] = {0.0, 0.0, 0.0};
        double vc[3] = {0.0, 0.0, 0.0};
        for (int ji = m_inl[iat]; ji < m_inl[iat + 1]; ++ji) {
            const int jat = m_nl[ji];
            const double vji[3] = {
                m_xyz(0, jat) + m_rvdw[jat] * m_grid(0, ig) - m_xyz(0, iat),
                m_xyz(1, jat) + m_rvdw[jat] * m_grid(1, ig) - m_xyz(1, iat),
                m_xyz(2, jat) + m_rvdw[jat] * m_grid(2, ig) - m_xyz(2, iat)};
            const double vvji = std::sqrt(vji[0] * vji[0] + vji[1] * vji[1] + vji[2] * vji[2]);
            const double tji = vvji / m_rvdw[iat];
            if (tji > thigh)
                continue;
            const double sji[3] = {vji[0] / vvji, vji[1] / vvji, vji[2] / vvji};
            dbasis(sji, basloc, dbsloc, vplm, vcos, vsin);
            double alp[3] = {0.0, 0.0, 0.0};
            double t = 1.0;
            const double* sig_i = sigma.col(iat).data();
            for (int l = 1; l <= lmax; ++l) {
                const int ind = idx(l, 0);
                const double fl = static_cast<double>(l);
                const double fac = t / m_facl[ind];
                for (int m = -l; m <= l; ++m) {
                    const double f2 = fac * sig_i[ind + m];
                    const double f1 = f2 * fl * basloc[ind + m];
                    for (int c = 0; c < 3; ++c)
                        alp[c] += f1 * sji[c] + f2 * dbsloc[c + 3 * (ind + m)];
                }
                t *= tji;
            }
            const double xji = fsw(tji, se, eta);
            const double oji = (m_fi(ig, jat) > 1.0) ? xji / m_fi(ig, jat) : xji;
            const double f1 = oji / m_rvdw[iat];
            for (int c = 0; c < 3; ++c)
                vb[c] += f1 * alp[c] * xi(ig, jat);
            if (tji > tlow) {
                const double beta = intmlp(tji, sig_i, basloc);
                double di, fac;
                if (m_fi(ig, jat) > 1.0) {
                    di = 1.0 / m_fi(ig, jat);
                    fac = di * xji;
                    bool proc = false;
                    double b = 0.0;
                    for (int jk = m_inl[jat]; jk < m_inl[jat + 1]; ++jk) {
                        const int kat = m_nl[jk];
                        const double vjk[3] = {
                            m_xyz(0, jat) + m_rvdw[jat] * m_grid(0, ig) - m_xyz(0, kat),
                            m_xyz(1, jat) + m_rvdw[jat] * m_grid(1, ig) - m_xyz(1, kat),
                            m_xyz(2, jat) + m_rvdw[jat] * m_grid(2, ig) - m_xyz(2, kat)};
                        const double vvjk = std::sqrt(vjk[0] * vjk[0] + vjk[1] * vjk[1] + vjk[2] * vjk[2]);
                        const double tjk = vvjk / m_rvdw[kat];
                        if (kat != iat && tjk <= thigh) {
                            proc = true;
                            const double sjk[3] = {vjk[0] / vvjk, vjk[1] / vvjk, vjk[2] / vvjk};
                            ylmbas(sjk, basloc, vplm, vcos, vsin);
                            const double g1 = intmlp(tjk, sigma.col(kat).data(), basloc);
                            const double xjk = fsw(tjk, se, eta);
                            b += g1 * xjk;
                        }
                    }
                    if (proc) {
                        const double g1 = di * di * dfsw(tji, se, eta) / m_rvdw[iat];
                        const double g2 = g1 * xi(ig, jat) * b;
                        for (int c = 0; c < 3; ++c)
                            vc[c] += g2 * sji[c];
                    }
                } else {
                    di = 1.0;
                    fac = 0.0;
                }
                const double f2 = (1.0 - fac) * di * dfsw(tji, se, eta) / m_rvdw[iat];
                for (int c = 0; c < 3; ++c)
                    vb[c] += f2 * xi(ig, jat) * beta * sji[c];
            }
        }
        for (int c = 0; c < 3; ++c)
            fx[c] -= m_w[ig] * (vb[c] - vc[c]);
    }
}

void DomainDecomposition::fdoga(int iat, const Mat& xi, const Mat& phi, double fx[3]) const
{
    const double swthr = 1.0 + (se + 1.0) * eta / 2.0;
    for (int ig = 0; ig < ngrid; ++ig) {
        double alp[3] = {0.0, 0.0, 0.0};
        if (m_ui(ig, iat) > 0.0 && m_ui(ig, iat) < 1.0) {
            const size_t zbase = static_cast<size_t>(3) * (ig + static_cast<size_t>(ngrid) * iat);
            for (int c = 0; c < 3; ++c)
                alp[c] += phi(ig, iat) * xi(ig, iat) * m_zi[zbase + c];
        }
        for (int ji = m_inl[iat]; ji < m_inl[iat + 1]; ++ji) {
            const int jat = m_nl[ji];
            const double vji[3] = {
                m_xyz(0, jat) + m_rvdw[jat] * m_grid(0, ig) - m_xyz(0, iat),
                m_xyz(1, jat) + m_rvdw[jat] * m_grid(1, ig) - m_xyz(1, iat),
                m_xyz(2, jat) + m_rvdw[jat] * m_grid(2, ig) - m_xyz(2, iat)};
            const double vvji = std::sqrt(vji[0] * vji[0] + vji[1] * vji[1] + vji[2] * vji[2]);
            const double tji = vvji / m_rvdw[iat];
            if (tji < swthr && tji > swthr - eta && m_ui(ig, jat) > 0.0) {
                const double sji[3] = {vji[0] / vvji, vji[1] / vvji, vji[2] / vvji};
                const double fac = -dfsw(tji, se, eta) / m_rvdw[iat];
                for (int c = 0; c < 3; ++c)
                    alp[c] += fac * phi(ig, jat) * xi(ig, jat) * sji[c];
            }
        }
        for (int c = 0; c < 3; ++c)
            fx[c] += m_w[ig] * alp[c];
    }
}

void DomainDecomposition::getDeriv(double keps, const std::vector<double>& phi,
                                   const Mat& sigma, const Mat& s, Mat& gradient) const
{
    Mat xi(ngrid, nat);
    for (int iat = 0; iat < nat; ++iat)
        for (int ig = 0; ig < ngrid; ++ig)
            xi(ig, iat) = s.col(iat).dot(m_basis.col(ig));

    Mat phiexp = Mat::Zero(ngrid, nat);
    int ii = 0;
    for (int iat = 0; iat < nat; ++iat)
        for (int ig = 0; ig < ngrid; ++ig)
            if (m_ui(ig, iat) > 0.0) {
                phiexp(ig, iat) = phi[ii];
                ++ii;
            }

    gradient.setZero();
    std::vector<double> basloc(nylm), dbsloc(3 * nylm), vplm(nylm), vcos(lmax + 1), vsin(lmax + 1);
    std::vector<double> xi_col(ngrid);
    for (int iat = 0; iat < nat; ++iat) {
        double fx[3] = {0.0, 0.0, 0.0};
        for (int ig = 0; ig < ngrid; ++ig)
            xi_col[ig] = xi(ig, iat);
        fdoka(iat, sigma, xi_col.data(), basloc.data(), dbsloc.data(), vplm.data(), vcos.data(), vsin.data(), fx);
        fdokb(iat, sigma, xi, basloc.data(), dbsloc.data(), vplm.data(), vcos.data(), vsin.data(), fx);
        fdoga(iat, xi, phiexp, fx);
        gradient(0, iat) = fx[0];
        gradient(1, iat) = fx[1];
        gradient(2, iat) = fx[2];
    }
    gradient *= keps;
}

void DomainDecomposition::getZeta(double keps, const Mat& s, std::vector<double>& zeta) const
{
    zeta.assign(ncav, 0.0);
    int ii = 0;
    for (int iat = 0; iat < nat; ++iat)
        for (int its = 0; its < ngrid; ++its)
            if (m_ui(its, iat) > 0.0) {
                zeta[ii] = keps * m_w[its] * m_ui(its, iat) * m_basis.col(its).dot(s.col(iat));
                ++ii;
            }
}

} // namespace Solvation
} // namespace Curcuma
