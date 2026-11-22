/*
 * <ddCOSMO implicit solvation for the native GFN1/GFN2 methods>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Claude Generated (Nov 2026): line-by-line port of
 *   xtb/src/solv/cosmo.f90          (TCosmo: radii, SASA term, addShift, getEnergy, addGradient)
 *   xtb/src/solv/ddcosmo/core.f90   (cavity, switching function, spherical harmonics, force kernels)
 *   xtb/src/solv/ddcosmo/solver.f90 (Jacobi/DIIS solver, L and L* actions, H^-1/2 norm)
 *   xtb/src/solv/sasa.f90           (numerical solvent accessible surface and its derivative)
 * Fortran indices are 1-based, here 0-based: (l, m) lives at l*l + l + m.
 * Claude Generated (Nov 2026): Model::TbliteCpcm, the CPCM of tblite 0.4.0 (src/tblite/solvation/cpcm.f90 on the
 * ddCOSMO core cpcm_dd.f90, which is xtb's core.f90 + solver.f90 with gradient instead of force sign convention):
 * only init() differs (epsilon table, keps, grid, eta, no surface term).
 * Claude Generated (Nov 2026): Model::XtbCosmoInf, xtb's conductor set gfn_cosmo_inf (include/param_cosmo_inf.fh, xtb
 * --cosmo inf), the COSMO step of xtb --cpcmx; cosmoSegments() returns the segment data xtb writes to xtb.cosmo
 * (cosmo.f90 writeCosmoFile) for the CPCM-X post-SCF step (cpcmx.h).
 * AI-generated, machine-tested only; human production testing pending.
 */

#include "ddcosmo.h"

#include "lebedev_xtb_grids.h"
#include "solvent_names.h"
#include "tblite_solvation_params.h"
#include "tblite_solvent_data.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/alpb_parameters.h"

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

using ALPBParameters::aatoau;
using ALPBParameters::fourpi;

constexpr double kPi = ALPBParameters::pi;
constexpr double kSe = -1.0;   // ddCOSMO switch shift: interior (core.f90 'se')
constexpr int kDiisMax = 25;   // cosmo.f90 ndiis
constexpr int kMaxIter = 200;  // cosmo.f90 n_iter

// Surface smoothing of the SASA integration (cosmo.f90 / sasa.f90: w = 0.3 Angstrom)
constexpr double kW = ALPBParameters::w_smooth;
constexpr double kAh0 = ALPBParameters::ah0;
constexpr double kAh1 = ALPBParameters::ah1;
constexpr double kAh3 = ALPBParameters::ah3;
constexpr double kTolSesp = 1.0e-6;          // sasa.f90 tolsesp
constexpr double kSurfaceTension = 1.0e-5;   // model.f90 surfaceTension
constexpr double kProbeRadAA = 1.3;          // param_cosmo.fh rprobe (Angstrom)
constexpr double kRadScale = 1.0;            // param_cosmo.fh c1 (bornScale)
constexpr double kROffsetAA = 2.0;           // model.f90 srcut, passed as rOffset

// tblite CPCM (cpcm.f90): keps is divided by 1 + alpha_alpb, alpha_alpb = 0.571412 (the ALPB shape constant);
// cpcm_input defaults nang = grid_size(6) = 74, eta = 2.0, lmax 6, conv 1e-8, rscale 1.
constexpr double kAlphaAlpb = 0.571412;
constexpr int kCpcmGrid = 74;
constexpr double kCpcmEta = 2.0;
constexpr int kCosmoGrid = 230;
constexpr double kCosmoEta = 0.2;

// xtb conductor set gfn_cosmo_inf (include/param_cosmo_inf.fh), used by xtb --cosmo inf and by the COSMO step of
// --cpcmx: epsilon = +inf (keps = 1/2), probe radius 2.94903643 Angstrom, element radius scales sx (vdwRad =
// vanDerWaalsRadCosmo * sx, model.f90 paramToModel), surface tensions gamscale. Elements not listed: gamscale 0, sx 1.
constexpr double kCosmoInfProbeAA = 2.94903643;
double cosmoInfGamscale(int z)
{
    switch (z) {
    case 1: return -0.35147423;
    case 6: return -0.26053405;
    case 7: return -0.85100961;
    case 8: return -0.38455056;
    case 9: return 0.02328270;
    case 15: return -0.03236934;
    case 16: return -0.25944025;
    case 17: return 0.01711509;
    default: return 0.0;
    }
}
double cosmoInfRadiusScale(int z)
{
    switch (z) {
    case 1: return 0.93069759;
    case 6: return 1.11953986;
    case 7: return 1.21082767;
    case 8: return 1.06487129;
    case 9: return 0.95536183;
    case 15: return 1.00338080;
    case 16: return 1.17965821;
    case 17: return 1.16340197;
    default: return 1.0;
    }
}

/// COSMO-optimised van-der-Waals radii in Angstrom (xtb/src/param/vdwradcosmo.f90), Z = 1..94.
constexpr double kStub = 2.223;
constexpr double kCosmoRadAA[94] = {
    1.3000, 1.6380, 1.5700, 1.0530, 2.0480, 2.0000, 1.8300, 1.7200, //  H-O
    1.7200, 1.8018, 1.8000, 1.6380, 2.1530, 2.2000, 2.1060, 2.1600, //  F-S
    2.0500, 2.2000, 2.2230, kStub,  kStub,  2.2930, kStub,  kStub,  // Cl-Cr
    kStub,  kStub,  kStub,  kStub,  kStub,  1.6260, kStub,  2.7000, // Mn-Ge
    2.3500, 2.2000, 2.1600, 2.3630, kStub,  kStub,  kStub,  kStub,  // As-Zr
    kStub,  kStub,  kStub,  kStub,  kStub,  kStub,  kStub,  kStub,  // Nb-Cd
    2.2580, 2.5500, 2.4100, 2.4100, 2.3200, 2.5270, kStub,  kStub,  // In-Ba
    kStub,  kStub,  kStub,  kStub,  kStub,  kStub,  kStub,  kStub,  // La-Gd
    kStub,  kStub,  kStub,  kStub,  kStub,  kStub,  kStub,  kStub,  // Tb-Hf
    kStub,  kStub,  kStub,  kStub,  kStub,  kStub,  kStub,  kStub,  // Ta-Hg
    kStub,  2.3600, 2.4220, 2.3050, 2.3630, 2.5740, kStub,  kStub,  // Tl-Ra
    kStub,  kStub,  kStub,  kStub,  kStub,  kStub                   // Ac-Pu
};

/// Dielectric constant from a plain number ("78.4"); false if the string is not a number.
bool parseNumber(const std::string& s, double& value)
{
    try {
        size_t pos = 0;
        value = std::stod(s, &pos);
        return pos == s.size();
    } catch (...) {
        return false;
    }
}

} // namespace

// ═════════════════════════════════════════════════════════════════════════════
// Parameters
// ═════════════════════════════════════════════════════════════════════════════

DdCosmoSolvation::~DdCosmoSolvation() = default;

void DdCosmoSolvation::setModel(Model model)
{
    m_model = model;
    m_nang_request = (model == Model::TbliteCpcm) ? kCpcmGrid : kCosmoGrid;
    m_eta = (model == Model::TbliteCpcm) ? kCpcmEta : kCosmoEta;
}

bool DdCosmoSolvation::solventAvailable(const std::string& method, const std::string& solvent, Model model)
{
    double eps = 0.0;
    if (model == Model::XtbCosmoInf)
        return true;   // the conductor needs no solvent parameters
    if (parseNumber(solvent, eps))
        return eps > 0.0;
    if (model == Model::TbliteCpcm)
        return curcuma::solvation::tbliteSolventEpsilon(solvent) > 0.0;
    return curcuma::solvation::getTbliteSolvationParam(method, true, curcuma::solvation::normalizeSolventName(solvent))
        != nullptr;
}

bool DdCosmoSolvation::init(const std::vector<int>& atomic_numbers, const std::string& solvent,
                            const std::string& method)
{
    m_initialized = false;
    m_nat = static_cast<int>(atomic_numbers.size());
    m_solvent = solvent;
    m_method = method;

    // Dielectric constant and surface tensions as xtb model.f90 loadInternalParam with cosmo=.true.: a named solvent
    // takes epsv and gamscale of its ALPB set of the method level (GFN1 or GFN2), everything else from gfn_cosmo
    // (probe 1.3 A, radius scale 1, COSMO radii, no H-bond term, no shift). A number is the dielectric constant
    // and keeps the gamscale of gfn_cosmo, which equals the gfn1/ALPB/water set (checked against param_cosmo.fh).
    const bool cpcm = (m_model == Model::TbliteCpcm);
    const bool conductor = (m_model == Model::XtbCosmoInf);
    const std::string key = curcuma::solvation::normalizeSolventName(solvent);
    const auto* named = curcuma::solvation::getTbliteSolvationParam(method, true, key);
    const auto* deflt = curcuma::solvation::getTbliteSolvationParam("gfn1", true, "water");
    double eps_number = 0.0;
    const bool is_number = parseNumber(solvent, eps_number);
    const double* gamscale = nullptr;
    const double cpcm_eps = cpcm ? curcuma::solvation::tbliteSolventEpsilon(solvent) : -1.0;
    if (conductor) {
        // xtb model.f90 with solvent 'inf': param = gfn_cosmo_inf, epsv = +inf; the solvent name is not used
        m_eps = std::numeric_limits<double>::infinity();
        m_param_source = "xtb conductor set gfn_cosmo_inf (--cosmo inf)";
    } else if (cpcm) {
        // tblite app/cli.f90 --cpcm: a solvent name from get_solvent_data, else the argument as a number
        if (m_eps_override > 0.0) {
            m_eps = m_eps_override;
            m_param_source = "explicit epsilon";
        } else if (cpcm_eps > 0.0) {
            m_eps = cpcm_eps;
            m_param_source = "tblite solvent table/" + solvent;
        } else if (is_number && eps_number > 0.0) {
            m_eps = eps_number;
            m_param_source = "epsilon from the solvent argument";
        } else {
            CurcumaLogger::error("CPCM: solvent '" + solvent + "' is neither in tblite's solvent table nor a dielectric constant");
            return false;
        }
    } else if (m_eps_override > 0.0) {
        m_eps = m_eps_override;
        gamscale = named ? named->gamscale : deflt->gamscale;
        m_param_source = named ? method + "/ALPB/" + key + " surface tensions, explicit epsilon"
                               : std::string("xtb default COSMO surface tensions, explicit epsilon");
    } else if (named) {
        m_eps = named->epsv;
        gamscale = named->gamscale;
        m_param_source = method + "/ALPB/" + key;
    } else if (is_number && eps_number > 0.0) {
        m_eps = eps_number;
        gamscale = deflt->gamscale;
        m_param_source = "xtb default COSMO set, epsilon from the solvent argument";
    } else {
        CurcumaLogger::error("COSMO: solvent '" + solvent + "' is neither a solvent with parameters nor a dielectric constant");
        return false;
    }
    if (!(m_eps > 0.0)) {
        CurcumaLogger::error("COSMO: dielectric constant must be positive");
        return false;
    }
    m_keps = 0.5 * (1.0 - 1.0 / m_eps);   // cosmo.f90 addShift (1/inf = 0: keps = 1/2 for the conductor)
    if (cpcm)
        m_keps /= (1.0 + kAlphaAlpb);     // tblite cpcm.f90 new_cpcm: -0.5 (1/eps - 1) / (1 + alpha_alpb)

    if (m_lmax < 1 || m_lmax > 20) {
        CurcumaLogger::error("COSMO: lmax must be between 1 and 20");
        return false;
    }

    // radii (cosmo.f90 initCosmo)
    m_rvdw.resize(m_nat);
    m_vdwsa.resize(m_nat);
    m_wrp.resize(m_nat);
    m_trj2.resize(2, m_nat);
    m_gamsasa.resize(m_nat);
    const double probe = (conductor ? kCosmoInfProbeAA : kProbeRadAA) * aatoau;
    for (int i = 0; i < m_nat; ++i) {
        const int z = atomic_numbers[i];
        if (z < 1 || z > 94) {
            CurcumaLogger::error(fmt::format("COSMO: no radius for element Z = {}", z));
            return false;
        }
        // model.f90 paramToModel: vdwRad = vanDerWaalsRadCosmo * sx (sx = 1 in gfn_cosmo, element scales in gfn_cosmo_inf)
        const double vdw = conductor ? (kCosmoRadAA[z - 1] * aatoau) * cosmoInfRadiusScale(z) : kCosmoRadAA[z - 1] * aatoau;
        m_rvdw(i) = kRadScale * vdw;
        m_vdwsa(i) = vdw + probe;
        m_trj2(0, i) = (m_vdwsa(i) - kW) * (m_vdwsa(i) - kW);
        m_trj2(1, i) = (m_vdwsa(i) + kW) * (m_vdwsa(i) + kW);
        double r = m_vdwsa(i) + kW;
        m_wrp(i) = (0.25 / kW + 3.0 * kAh3 * (0.2 * r * r - 0.5 * r * m_vdwsa(i) + m_vdwsa(i) * m_vdwsa(i) / 3.0)) * r * r * r;
        r = m_vdwsa(i) - kW;
        m_wrp(i) -= (0.25 / kW + 3.0 * kAh3 * (0.2 * r * r - 0.5 * r * m_vdwsa(i) + m_vdwsa(i) * m_vdwsa(i) / 3.0)) * r * r * r;
        const double gam = conductor ? cosmoInfGamscale(z) : (cpcm ? 0.0 : gamscale[z - 1]);
        m_gamsasa(i) = cpcm ? 0.0 : gam * fourpi * kSurfaceTension;   // CPCM: no surface term
    }
    m_srcut = 2.0 * (kW + m_vdwsa.maxCoeff()) + kROffsetAA * aatoau;

    // grids: COSMO cavity grid nearest to nAng/2, surface grid nearest to nAng (bisectSearch as xtb);
    // CPCM cavity grid nearest to nAng (tblite list_bisection), no surface grid
    namespace leb = curcuma::solvation::lebedev;
    const int icav = leb::bisectGridIndex(cpcm ? m_nang_request : m_nang_request / 2);
    const int isas = cpcm ? icav : leb::bisectGridIndex(m_nang_request);
    std::vector<double> gx, gw;
    if (icav < 0 || isas < 0 || !leb::getAngGrid(icav, gx, gw)) {
        CurcumaLogger::error(fmt::format("COSMO: no Lebedev grid for {} points", m_nang_request));
        return false;
    }
    m_ngrid = static_cast<int>(gw.size());
    m_grid.resize(3, m_ngrid);
    m_w.resize(m_ngrid);
    for (int i = 0; i < m_ngrid; ++i) {
        for (int c = 0; c < 3; ++c)
            m_grid(c, i) = gx[3 * i + c];
        m_w(i) = gw[i] * 4.0 * kPi;   // core.f90: weights are normalised
    }
    if (cpcm) {
        gx.clear();
        gw.clear();
    } else {
        leb::getAngGrid(isas, gx, gw);
    }
    m_sgrid.resize(3, static_cast<int>(gw.size()));
    m_sweight.resize(static_cast<int>(gw.size()));
    for (int i = 0; i < static_cast<int>(gw.size()); ++i) {
        for (int c = 0; c < 3; ++c)
            m_sgrid(c, i) = gx[3 * i + c];
        m_sweight(i) = gw[i];
    }

    // spherical-harmonics normalisation (core.f90 initDomainDecomposition)
    m_nylm = (m_lmax + 1) * (m_lmax + 1);
    m_fact.assign(2 * m_lmax + 1, 1.0);   // m_fact[k] = k!
    for (int k = 1; k < 2 * m_lmax + 1; ++k)
        m_fact[k] = k * m_fact[k - 1];
    m_facl.assign(m_nylm, 0.0);
    m_facs.assign(m_nylm, 0.0);
    for (int l = 0; l <= m_lmax; ++l) {
        const int i0 = ind(l);
        const double fl = (2.0 * l + 1.0) / (4.0 * kPi);
        const double ffl = std::sqrt(fl);
        for (int m = -l; m <= l; ++m)
            m_facl[i0 + m] = fl;
        m_facs[i0] = ffl;
        for (int m = 1; m <= l; ++m) {
            double fnorm = std::sqrt(2.0) * ffl * std::sqrt(m_fact[l - m] / m_fact[l + m]);
            if (m % 2 == 1)
                fnorm = -fnorm;
            m_facs[i0 + m] = fnorm;
            m_facs[i0 - m] = fnorm;
        }
    }
    m_basis.resize(m_nylm, m_ngrid);
    std::vector<double> vplm(m_nylm), vcos(m_lmax + 1), vsin(m_lmax + 1);
    for (int i = 0; i < m_ngrid; ++i)
        ylmbas(m_grid.col(i).data(), m_basis.col(i).data(), vplm.data(), vcos.data(), vsin.data());

    m_sigma = Eigen::MatrixXd::Zero(m_nylm, m_nat);
    m_s = Eigen::MatrixXd::Zero(m_nylm, m_nat);
    m_psi = Eigen::MatrixXd::Zero(m_nylm, m_nat);
    m_have_sigma = m_have_s = false;
    m_initialized = true;
    if (CurcumaLogger::get_verbosity() >= 2)
        printInfo();
    return true;
}

void DdCosmoSolvation::printInfo() const
{
    CurcumaLogger::info(m_model == Model::TbliteCpcm   ? "Native CPCM (ddCOSMO solver, port of tblite 0.4.0)"
                        : m_model == Model::XtbCosmoInf ? "Native COSMO conductor (ddCOSMO, xtb 6.7.1 --cosmo inf, CPCM-X step)"
                                                        : "Native COSMO (ddCOSMO, port of xtb 6.7.1)");
    CurcumaLogger::param("solvent", m_solvent);
    CurcumaLogger::param("parameters", m_param_source);
    CurcumaLogger::param("dielectric constant", fmt::format("{:.4f}", m_eps));
    CurcumaLogger::param("lmax / eta / conv", fmt::format("{} / {} / {:.1e}", m_lmax, m_eta, m_conv));
    CurcumaLogger::param("grid points (cavity / surface)", fmt::format("{} / {}", m_ngrid, m_sweight.size()));
}

// ═════════════════════════════════════════════════════════════════════════════
// Geometry update: SASA, neighbour list, characteristic functions, cavity, J
// ═════════════════════════════════════════════════════════════════════════════

void DdCosmoSolvation::update(const std::vector<int>& atomic_numbers, const Matrix& xyz_bohr)
{
    (void)atomic_numbers;
    if (!m_initialized)
        return;
    m_xyz.resize(3, m_nat);
    for (int i = 0; i < m_nat; ++i)
        for (int c = 0; c < 3; ++c)
            m_xyz(c, i) = xyz_bohr(i, c);

    computeSurface(xyz_bohr);

    // cosmo.f90 update deallocates sigma and s: no restart across geometries
    m_have_sigma = m_have_s = false;

    buildNeighbourList();
    buildCharacteristicFunctions();
    buildCavity();

    // getCoulombMatrix: J(ic, j) = 1/|ccav_ic - R_j|
    m_jmat.resize(m_ncav, m_nat);
    for (int j = 0; j < m_nat; ++j)
        for (int ic = 0; ic < m_ncav; ++ic)
            m_jmat(ic, j) = 1.0 / (m_ccav.col(ic) - m_xyz.col(j)).norm();
}

/// Numerical solvent accessible surface (sasa.f90 compute_numsa, compute_w_sp) and its contracted
/// derivative dsdr(:, j) = sum_i gamma_i d A_i / d R_j. Neighbour list as cosmo.f90 update_nnlist_sasa.
void DdCosmoSolvation::computeSurface(const Matrix& xyz)
{
    std::vector<std::vector<int>> nn(m_nat);
    const double srcut2 = m_srcut * m_srcut;
    for (int i1 = 0; i1 < m_nat; ++i1)
        for (int i2 = 0; i2 < i1; ++i2) {
            const double dx = xyz(i1, 0) - xyz(i2, 0), dy = xyz(i1, 1) - xyz(i2, 1), dz = xyz(i1, 2) - xyz(i2, 2);
            if (dx * dx + dy * dy + dz * dz < srcut2) {
                nn[i1].push_back(i2);
                nn[i2].push_back(i1);
            }
        }

    m_sasa = Eigen::VectorXd::Zero(m_nat);
    m_dsdr = Eigen::MatrixXd::Zero(3, m_nat);
    const int nang = static_cast<int>(m_sweight.size());
    Eigen::MatrixXd grads(3, m_nat);
    std::vector<double> grds;
    std::vector<int> grdi;
    for (int iat = 0; iat < m_nat; ++iat) {
        const double rsas = m_vdwsa(iat);
        grads.setZero();
        double sasai = 0.0;
        const int nno = static_cast<int>(nn[iat].size());
        grds.assign(3 * std::max(nno, 1), 0.0);
        grdi.assign(std::max(nno, 1), 0);
        for (int ip = 0; ip < nang; ++ip) {
            const double px = xyz(iat, 0) + rsas * m_sgrid(0, ip);
            const double py = xyz(iat, 1) + rsas * m_sgrid(1, ip);
            const double pz = xyz(iat, 2) + rsas * m_sgrid(2, ip);
            // compute_w_sp
            int nni = 0;
            double sasap = 1.0;
            for (int k = 0; k < nno; ++k) {
                const int ia = nn[iat][k];
                const double tx = px - xyz(ia, 0), ty = py - xyz(ia, 1), tz = pz - xyz(ia, 2);
                const double tj2 = tx * tx + ty * ty + tz * tz;
                if (tj2 < m_trj2(1, ia)) {
                    if (tj2 <= m_trj2(0, ia)) {
                        sasap = 0.0;
                        break;
                    }
                    const double sqtj = std::sqrt(tj2);
                    const double uj = sqtj - m_vdwsa(ia);
                    const double ah3uj2 = kAh3 * uj * uj;
                    double dsasaij = kAh1 + 3.0 * ah3uj2;
                    double sasaij = kAh0 + (kAh1 + ah3uj2) * uj;
                    sasaij = std::max(0.0, std::min(sasaij, 1.0));
                    if (sasaij < std::numeric_limits<double>::min()) {
                        sasap = 0.0;
                        break;
                    }
                    sasap *= sasaij;
                    dsasaij = dsasaij / (sasaij * sqtj);
                    grdi[nni] = ia;
                    grds[3 * nni + 0] = dsasaij * tx;
                    grds[3 * nni + 1] = dsasaij * ty;
                    grds[3 * nni + 2] = dsasaij * tz;
                    ++nni;
                }
            }
            if (sasap > kTolSesp) {
                const double wsa = m_sweight(ip) * m_wrp(iat) * sasap;
                sasai += wsa;
                for (int jj = 0; jj < nni; ++jj)
                    for (int c = 0; c < 3; ++c) {
                        const double d = wsa * grds[3 * jj + c];
                        grads(c, iat) += d;
                        grads(c, grdi[jj]) -= d;
                    }
            }
        }
        m_sasa(iat) = sasai;
        m_dsdr += grads * m_gamsasa(iat);   // mctc_gemv(dsdrt, gamsasa, dsdr)
    }
    m_gsasa = m_sasa.dot(m_gamsasa);
}

/// core.f90 mknnl: spheres i and j are neighbours if they overlap, |R_i - R_j| <= r_i + r_j.
void DdCosmoSolvation::buildNeighbourList()
{
    m_inl.assign(m_nat + 1, 0);
    m_nl.clear();
    for (int iat = 0; iat < m_nat; ++iat) {
        m_inl[iat] = static_cast<int>(m_nl.size());
        for (int jat = 0; jat < m_nat; ++jat) {
            if (iat == jat)
                continue;
            const double d2 = (m_xyz.col(iat) - m_xyz.col(jat)).squaredNorm();
            const double r2 = (m_rvdw(iat) + m_rvdw(jat)) * (m_rvdw(iat) + m_rvdw(jat));
            if (d2 <= r2)
                m_nl.push_back(jat);
        }
    }
    m_inl[m_nat] = static_cast<int>(m_nl.size());
}

/// core.f90 mkfiui: fi(n,i) = sum_j chi(t_n^ij), U_i(n) = 1 - fi for fi <= 1, and the derivative zi.
void DdCosmoSolvation::buildCharacteristicFunctions()
{
    m_fi = Eigen::MatrixXd::Zero(m_ngrid, m_nat);
    m_ui = Eigen::MatrixXd::Zero(m_ngrid, m_nat);
    m_zi.assign(3 * m_ngrid * m_nat, 0.0);
    const double swthr = 1.0 + (kSe + 1.0) * m_eta / 2.0;
    for (int iat = 0; iat < m_nat; ++iat)
        for (int i = 0; i < m_ngrid; ++i) {
            for (int ii = m_inl[iat]; ii < m_inl[iat + 1]; ++ii) {
                const int jat = m_nl[ii];
                const Eigen::Vector3d v = m_xyz.col(iat) + m_rvdw(iat) * m_grid.col(i) - m_xyz.col(jat);
                const double vv = v.norm();
                const double t = vv / m_rvdw(jat);
                const double xt = fsw(t, m_eta);
                if (t < swthr && t > swthr - m_eta) {
                    const double fac = dfsw(t, m_eta) / m_rvdw(jat);
                    for (int c = 0; c < 3; ++c)
                        m_zi[c + 3 * (i + m_ngrid * iat)] += fac * v(c) / vv;
                }
                m_fi(i, iat) += xt;
            }
            if (m_fi(i, iat) <= 1.0)
                m_ui(i, iat) = 1.0 - m_fi(i, iat);
        }
}

/// core.f90 mkccav: exposed grid points (U > 0) in sphere-major order.
void DdCosmoSolvation::buildCavity()
{
    m_ncav = 0;
    for (int iat = 0; iat < m_nat; ++iat)
        for (int i = 0; i < m_ngrid; ++i)
            if (m_ui(i, iat) > 0.0)
                ++m_ncav;
    m_ccav.resize(3, m_ncav);
    int ii = 0;
    for (int iat = 0; iat < m_nat; ++iat)
        for (int i = 0; i < m_ngrid; ++i)
            if (m_ui(i, iat) > 0.0)
                m_ccav.col(ii++) = m_xyz.col(iat) + m_rvdw(iat) * m_grid.col(i);
}

// ═════════════════════════════════════════════════════════════════════════════
// Switching function and spherical harmonics (core.f90)
// ═════════════════════════════════════════════════════════════════════════════

/// 5th-degree switch chi(t): 1 inside (t <= 1-eta), 0 outside (t >= 1), interior shift se = -1.
double DdCosmoSolvation::fsw(double t, double eta)
{
    const double x = t - (kSe + 1.0) * eta / 2.0;
    const double flow = 1.0 - eta;
    if (x >= 1.0)
        return 0.0;
    if (x <= flow)
        return 1.0;
    const double a = 15.0 * eta - 12.0;
    const double b = 10.0 * eta * eta - 15.0 * eta + 6.0;
    return ((x - 1.0) * (x - 1.0) * (1.0 - x) * (6.0 * x * x + a * x + b)) / std::pow(eta, 5);
}

double DdCosmoSolvation::dfsw(double t, double eta)
{
    const double x = t - (kSe + 1.0) * eta / 2.0;
    const double flow = 1.0 - eta;
    if (x >= 1.0 || x <= flow)
        return 0.0;
    return 30.0 * (1.0 - x) * (x - 1.0) * (x - 1.0 + eta) * (x - 1.0 + eta) / std::pow(eta, 5);
}

/// cos(m phi), sin(m phi) by the Chebyshev recurrence (core.f90 trgev).
void DdCosmoSolvation::trgev(double x, double y, double* cx, double* sx) const
{
    cx[0] = 1.0;
    sx[0] = 0.0;
    cx[1] = x;
    sx[1] = y;
    for (int m = 2; m <= m_lmax; ++m) {
        cx[m] = 2.0 * x * cx[m - 1] - cx[m - 2];
        sx[m] = 2.0 * x * sx[m - 1] - sx[m - 2];
    }
}

/// Associated Legendre functions P_l^m(cos theta) (core.f90 polleg), stored at (l, m >= 0).
void DdCosmoSolvation::polleg(double x, double y, double* plm) const
{
    double fact = 1.0, pmm = 1.0;
    const double somx2 = y;
    for (int m = 0; m <= m_lmax; ++m) {
        const int i0 = (m + 1) * (m + 1) - 1;   // (l = m, m)
        plm[i0] = pmm;
        if (m == m_lmax)
            return;
        const double fm = m;
        double pmm1 = x * (2.0 * fm + 1.0) * pmm;
        plm[i0 + 2 * m + 2] = pmm1;             // (l = m+1, m)
        const double pmmo = pmm;
        for (int l = m + 2; l <= m_lmax; ++l) {
            const double fl = l;
            const double pll = (x * (2.0 * fl - 1.0) * pmm1 - (fl + fm - 1.0) * pmm) / (fl - fm);
            plm[l * l + l + m] = pll;
            pmm = pmm1;
            pmm1 = pll;
        }
        pmm = -pmmo * fact * somx2;
        fact += 2.0;
    }
}

/// Real spherical harmonics Y_lm at the unit vector x (core.f90 ylmbas).
void DdCosmoSolvation::ylmbas(const double* x, double* basloc, double* vplm, double* vcos, double* vsin) const
{
    const double cthe = x[2];
    const double sthe = std::sqrt(std::max(0.0, 1.0 - cthe * cthe));
    double cphi = 0.0, sphi = 0.0;
    if (sthe != 0.0) {
        cphi = x[0] / sthe;
        sphi = x[1] / sthe;
        trgev(cphi, sphi, vcos, vsin);
    } else {
        for (int m = 0; m <= m_lmax; ++m) {
            vcos[m] = 1.0;
            vsin[m] = 0.0;
        }
    }
    polleg(cthe, sthe, vplm);
    for (int l = 0; l <= m_lmax; ++l) {
        const int i0 = ind(l);
        basloc[i0] = m_facs[i0] * vplm[i0];
        for (int m = 1; m <= l; ++m) {
            const double plm = vplm[i0 + m];
            basloc[i0 + m] = m_facs[i0 + m] * plm * vcos[m];
            basloc[i0 - m] = m_facs[i0 - m] * plm * vsin[m];
        }
    }
}

/// Y_lm and their gradients on the unit sphere, dbsloc 3 x nylm (core.f90 dbasis).
void DdCosmoSolvation::dbasis(const double* x, double* basloc, double* dbsloc, double* vplm, double* vcos,
                              double* vsin) const
{
    const double cthe = x[2];
    const double sthe = std::sqrt(std::max(0.0, 1.0 - cthe * cthe));
    double cphi, sphi;
    if (sthe != 0.0) {
        cphi = x[0] / sthe;
        sphi = x[1] / sthe;
    } else {
        cphi = 1.0;
        sphi = 0.0;
    }
    const double et[3] = { cthe * cphi, cthe * sphi, -sthe };
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
        for (int m = 0; m <= m_lmax; ++m) {
            vcos[m] = 1.0;
            vsin[m] = 0.0;
        }
    }
    // core.f90 declares VC, VS as integers: VC = 0, VS = int(cthe) (only used at the poles, where cthe = +-1)
    const double VC = 0.0;
    const double VS = static_cast<double>(static_cast<int>(cthe));
    polleg(cthe, sthe, vplm);

    std::fill(basloc, basloc + m_nylm, 0.0);
    std::fill(dbsloc, dbsloc + 3 * m_nylm, 0.0);
    for (int l = 0; l <= m_lmax; ++l) {
        const int i0 = ind(l);
        double fln = m_facs[i0];
        basloc[i0] = fln * vplm[i0];
        if (l > 0)
            for (int c = 0; c < 3; ++c)
                dbsloc[3 * i0 + c] = fln * vplm[i0 + 1] * et[c];
        for (int m = 1; m <= l; ++m) {
            fln = m_facs[i0 + m];
            const double plm = fln * vplm[i0 + m];
            const double pp1 = (m < l) ? -0.5 * vplm[i0 + m + 1] : 0.0;
            const double pm1 = 0.5 * (static_cast<double>(l + m) * static_cast<double>(l - m + 1) * vplm[i0 + m - 1]);
            const double pp = pp1 + pm1;
            basloc[i0 + m] = plm * vcos[m];
            for (int c = 0; c < 3; ++c) {
                if (sthe != 0.0)
                    dbsloc[3 * (i0 + m) + c] = -fln * pp * vcos[m] * et[c] - m * plm * vsin[m] * ep[c];
                else
                    dbsloc[3 * (i0 + m) + c] = -fln * pp * vcos[m] * et[c] - fln * pp * ep[c] * VC;
            }
            basloc[i0 - m] = plm * vsin[m];
            for (int c = 0; c < 3; ++c) {
                if (sthe != 0.0)
                    dbsloc[3 * (i0 - m) + c] = -fln * pp * vsin[m] * et[c] + m * plm * vcos[m] * ep[c];
                else
                    dbsloc[3 * (i0 - m) + c] = -fln * pp * vsin[m] * et[c] - fln * pp * ep[c] * VS;
            }
        }
    }
}

/// sum_lm 4 pi/(2l+1) t^l Y_lm(s) sigma_lm (core.f90 intmlp).
double DdCosmoSolvation::intmlp(double t, const double* sigma, const double* basloc) const
{
    double tt = 1.0, ss = 0.0;
    for (int l = 0; l <= m_lmax; ++l) {
        const int i0 = ind(l);
        const double fac = tt / m_facl[i0];
        double d = 0.0;
        for (int m = -l; m <= l; ++m)
            d += basloc[i0 + m] * sigma[i0 + m];
        ss += fac * d;
        tt *= t;
    }
    return ss;
}

/// Projection of a grid function on Y_lm: xlm = sum_n Y_lm(s_n) w_n x_n (core.f90 intrhs).
void DdCosmoSolvation::intrhs(const double* x, double* xlm) const
{
    std::fill(xlm, xlm + m_nylm, 0.0);
    for (int ig = 0; ig < m_ngrid; ++ig) {
        const double f = m_w(ig) * x[ig];
        for (int k = 0; k < m_nylm; ++k)
            xlm[k] += m_basis(k, ig) * f;
    }
}

/// Potential of the neighbours' expansions at the buried points of sphere iat (core.f90 calcv).
void DdCosmoSolvation::calcv(int iat, const Eigen::MatrixXd& sigma, double* pot) const
{
    std::vector<double> basloc(m_nylm), vplm(m_nylm), vcos(m_lmax + 1), vsin(m_lmax + 1);
    std::fill(pot, pot + m_ngrid, 0.0);
    for (int its = 0; its < m_ngrid; ++its) {
        if (!(m_ui(its, iat) < 1.0))
            continue;
        for (int ij = m_inl[iat]; ij < m_inl[iat + 1]; ++ij) {
            const int jat = m_nl[ij];
            const Eigen::Vector3d vij = m_xyz.col(iat) + m_rvdw(iat) * m_grid.col(its) - m_xyz.col(jat);
            const double vvij = vij.norm();
            const double tij = vvij / m_rvdw(jat);
            if (tij < 1.0) {
                const Eigen::Vector3d sij = vij / vvij;
                const double xij = fsw(tij, m_eta);
                const double oij = (m_fi(its, iat) > 1.0) ? xij / m_fi(its, iat) : xij;
                ylmbas(sij.data(), basloc.data(), vplm.data(), vcos.data(), vsin.data());
                pot[its] += oij * intmlp(tij, sigma.col(jat).data(), basloc.data());
            }
        }
    }
}

/// Negative action of the off-diagonal blocks of L^* on sphere iat (core.f90 adjrhs).
void DdCosmoSolvation::adjrhs(int iat, const Eigen::MatrixXd& xi, double* vlm) const
{
    std::vector<double> basloc(m_nylm), vplm(m_nylm), vcos(m_lmax + 1), vsin(m_lmax + 1);
    for (int ij = m_inl[iat]; ij < m_inl[iat + 1]; ++ij) {
        const int jat = m_nl[ij];
        for (int ig = 0; ig < m_ngrid; ++ig) {
            const Eigen::Vector3d vji = m_xyz.col(jat) + m_rvdw(jat) * m_grid.col(ig) - m_xyz.col(iat);
            const double vvji = vji.norm();
            const double tji = vvji / m_rvdw(iat);
            if (tji < (1.0 + (kSe + 1.0) / 2.0 * m_eta)) {
                const Eigen::Vector3d sji = vji / vvji;
                const double xji = fsw(tji, m_eta);
                const double oji = (m_fi(ig, jat) > 1.0) ? xji / m_fi(ig, jat) : xji;
                ylmbas(sji.data(), basloc.data(), vplm.data(), vcos.data(), vsin.data());
                double t = 1.0;
                const double fac = m_w(ig) * xi(ig, jat) * oji;
                for (int l = 0; l <= m_lmax; ++l) {
                    const int i0 = ind(l);
                    const double ffac = fac * t / m_facl[i0];
                    for (int m = -l; m <= l; ++m)
                        vlm[i0 + m] += ffac * basloc[i0 + m];
                    t *= tji;
                }
            }
        }
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// Linear solver (solver.f90)
// ═════════════════════════════════════════════════════════════════════════════

/// y = L x, off-diagonal blocks only (solver.f90 lx).
void DdCosmoSolvation::lx(const Eigen::MatrixXd& x, Eigen::MatrixXd& y) const
{
    y = Eigen::MatrixXd::Zero(m_nylm, m_nat);
    std::vector<double> pot(m_ngrid);
    for (int iat = 0; iat < m_nat; ++iat) {
        calcv(iat, x, pot.data());
        intrhs(pot.data(), y.col(iat).data());
        y.col(iat) = -y.col(iat);
    }
}

/// y = L^* x, off-diagonal blocks only (solver.f90 lstarx).
void DdCosmoSolvation::lstarx(const Eigen::MatrixXd& x, Eigen::MatrixXd& y) const
{
    y = Eigen::MatrixXd::Zero(m_nylm, m_nat);
    const Eigen::MatrixXd xi = m_basis.transpose() * x;   // xi(ig, iat) = sum_lm x(lm, iat) Y_lm(s_ig)
    for (int iat = 0; iat < m_nat; ++iat) {
        adjrhs(iat, xi, y.col(iat).data());
        y.col(iat) = -y.col(iat);
    }
}

/// rms over spheres of the H^-1/2 norm sqrt(sum_l 1/(l+1) sum_m u_lm^2) (solver.f90 hnorm, core.f90 hsnorm).
double DdCosmoSolvation::hnorm(const Eigen::MatrixXd& x) const
{
    double vrms = 0.0;
    for (int iat = 0; iat < m_nat; ++iat) {
        double u = 0.0;
        for (int l = 0; l <= m_lmax; ++l) {
            const double fac = 1.0 / (1.0 + l);
            for (int m = -l; m <= l; ++m)
                u += fac * x(ind(l) + m, iat) * x(ind(l) + m, iat);
        }
        u = std::sqrt(u);
        vrms += u * u;
    }
    return std::sqrt(vrms / static_cast<double>(m_nat));
}

/**
 * @brief Jacobi iterations with DIIS extrapolation (solver.f90 jacobi_diis, diis, makeb, gjinv).
 *
 * L = D + O with the diagonal blocks D = (4 pi / (2l+1)) I, so x_new = D^-1 (rhs - O x); the increment
 * is measured in the H^-1/2 norm (norm = 4 in xtb). DIIS keeps up to 25 vectors and drops the oldest 10
 * when full. Returns false if not converged in 200 iterations.
 */
bool DdCosmoSolvation::jacobiDiis(bool adjoint, double tol, const Eigen::MatrixXd& rhs, Eigen::MatrixXd& x) const
{
    const int n = m_nylm * m_nat;
    Eigen::MatrixXd xd(n, kDiisMax), ed(n, kDiisMax), b = Eigen::MatrixXd::Zero(kDiisMax + 1, kDiisMax + 1);
    int nmat = 1;   // 1-based count as in the Fortran
    Eigen::MatrixXd y, xnew(m_nylm, m_nat);
    auto B = [&](int i, int j) -> double& { return b(i - 1, j - 1); };

    for (int it = 1; it <= kMaxIter; ++it) {
        if (adjoint)
            lstarx(x, y);
        else
            lx(x, y);
        y = rhs - y;
        for (int iat = 0; iat < m_nat; ++iat)
            for (int k = 0; k < m_nylm; ++k)
                xnew(k, iat) = m_facl[k] * y(k, iat);   // ldm1x

        // DIIS
        xd.col(nmat - 1) = Eigen::Map<const Eigen::VectorXd>(xnew.data(), n);
        ed.col(nmat - 1) = Eigen::Map<const Eigen::VectorXd>(xnew.data(), n) - Eigen::Map<const Eigen::VectorXd>(x.data(), n);
        if (nmat >= kDiisMax) {
            for (int j = 2; j <= nmat - 10; ++j)
                for (int k = 2; k <= nmat - 10; ++k)
                    B(j, k) = B(j + 10, k + 10);
            for (int j = 1; j <= nmat - 10; ++j) {
                xd.col(j - 1) = xd.col(j + 9);
                ed.col(j - 1) = ed.col(j + 9);
            }
            nmat -= 10;
        }
        const int nmat1 = nmat + 1;
        // makeb
        if (nmat == 1) {
            B(1, 1) = 0.0;
            B(1, 2) = 1.0;
            B(2, 1) = 1.0;
            B(2, 2) = ed.col(0).squaredNorm();
        } else {
            B(nmat + 1, 1) = 1.0;
            B(1, nmat + 1) = 1.0;
            for (int i = 1; i <= nmat - 1; ++i) {
                const double bij = ed.col(i - 1).dot(ed.col(nmat - 1));
                B(nmat + 1, i + 1) = bij;
                B(i + 1, nmat + 1) = bij;
            }
            B(nmat + 1, nmat + 1) = ed.col(nmat - 1).squaredNorm();
        }
        // gjinv: Gauss-Jordan elimination with full pivoting on bloc, one right-hand side cex = e_1
        Eigen::MatrixXd a = b.topLeftCorner(nmat1, nmat1);
        Eigen::VectorXd cex = Eigen::VectorXd::Zero(nmat1);
        cex(0) = 1.0;
        {
            std::vector<int> indxc(nmat1), indxr(nmat1), piv(nmat1, 0);
            int irow = 0, icol = 0;
            bool ok = true;
            for (int i = 0; i < nmat1 && ok; ++i) {
                double big = 0.0;
                for (int j = 0; j < nmat1; ++j)
                    if (piv[j] != 1)
                        for (int k = 0; k < nmat1; ++k)
                            if (piv[k] == 0 && std::abs(a(j, k)) > big) {
                                big = std::abs(a(j, k));
                                irow = j;
                                icol = k;
                            }
                piv[icol] += 1;
                if (piv[icol] > 1) { ok = false; break; }
                if (irow != icol) {
                    a.row(irow).swap(a.row(icol));
                    std::swap(cex(irow), cex(icol));
                }
                indxr[i] = irow;
                indxc[i] = icol;
                if (a(icol, icol) == 0.0) { ok = false; break; }
                const double pinv = 1.0 / a(icol, icol);
                a(icol, icol) = 1.0;
                a.row(icol) *= pinv;
                cex(icol) *= pinv;
                for (int j = 0; j < nmat1; ++j)
                    if (j != icol) {
                        const double dum = a(j, icol);
                        a(j, icol) = 0.0;
                        a.row(j) -= a.row(icol) * dum;
                        cex(j) -= cex(icol) * dum;
                    }
            }
            if (!ok)
                throw std::runtime_error("ddCOSMO: singular DIIS matrix");
        }
        Eigen::VectorXd xe = Eigen::VectorXd::Zero(n);
        for (int i = 1; i <= nmat; ++i)
            xe += cex(i) * xd.col(i - 1);
        nmat += 1;
        Eigen::Map<Eigen::VectorXd>(xnew.data(), n) = xe;

        const double rms_norm = hnorm(xnew - x);
        x = xnew;
        if (rms_norm < tol)
            return true;
    }
    return false;
}

// ═════════════════════════════════════════════════════════════════════════════
// COSMO model (cosmo.f90)
// ═════════════════════════════════════════════════════════════════════════════

/// Direct problem L X = g(Phi) (cosmo.f90 solveCosmoDirect with cart = .true.).
void DdCosmoSolvation::solveDirect(const Eigen::VectorXd& phi, bool restart) const
{
    // wghpot: g(n, i) = -U_i(n) Phi(n) on the exposed points, then intrhs per sphere
    Eigen::MatrixXd rhs(m_nylm, m_nat);
    std::vector<double> g(m_ngrid);
    int ic = 0;
    for (int iat = 0; iat < m_nat; ++iat) {
        for (int ig = 0; ig < m_ngrid; ++ig) {
            g[ig] = 0.0;
            if (m_ui(ig, iat) != 0.0) {
                g[ig] = -m_ui(ig, iat) * phi(ic);
                ++ic;
            }
        }
        intrhs(g.data(), rhs.col(iat).data());
    }
    if (!restart)
        for (int iat = 0; iat < m_nat; ++iat)
            for (int k = 0; k < m_nylm; ++k)
                m_sigma(k, iat) = m_facl[k] * rhs(k, iat);
    if (!jacobiDiis(false, m_conv, rhs, m_sigma))
        throw std::runtime_error("direct ddCOSMO did not converge");
    m_have_sigma = true;
}

/// Adjoint problem L^* S = Psi (cosmo.f90 solveCosmoAdjoint).
void DdCosmoSolvation::solveAdjoint(bool restart, double tol) const
{
    if (!restart)
        for (int iat = 0; iat < m_nat; ++iat)
            for (int k = 0; k < m_nylm; ++k)
                m_s(k, iat) = m_facl[k] * m_psi(k, iat);
    if (!jacobiDiis(true, tol, m_psi, m_s))
        throw std::runtime_error("adjoint ddCOSMO did not converge");
    m_have_s = true;
}

/// Direct solution X and Psi for the charges q (restart from the previous X, as addShift does).
void DdCosmoSolvation::ensureSolution(const Vector& q) const
{
    if (m_have_sigma && m_q_sigma.size() == q.size() && m_q_sigma == q)
        return;
    const Eigen::VectorXd phi = m_jmat * q;   // getPhi
    solveDirect(phi, m_have_sigma);
    m_psi.setZero();                          // getPsi
    for (int iat = 0; iat < m_nat; ++iat)
        m_psi(0, iat) = std::sqrt(fourpi) * q(iat);
    m_q_sigma = q;
}

/// zeta(n) = keps w_n U_i(n) sum_lm Y_lm(s_n) S_i,lm on the exposed points (cosmo.f90 getZeta).
void DdCosmoSolvation::getZeta(Eigen::VectorXd& zeta) const
{
    zeta.resize(m_ncav);
    int ii = 0;
    for (int iat = 0; iat < m_nat; ++iat)
        for (int its = 0; its < m_ngrid; ++its)
            if (m_ui(its, iat) > 0.0) {
                zeta(ii) = m_keps * m_w(its) * m_ui(its, iat) * m_basis.col(its).dot(m_s.col(iat));
                ++ii;
            }
}

/// Segment data as xtb cosmo.f90 writeCosmoFile: zeta(ii) = w(ig) ui(ig,iat) <basis(:,ig), s(:,iat)> (no keps),
/// area(ii) = w(ig) autoaa^2 rvdw(iat)^2, potential = getPhi(q). The adjoint is re-solved first with the tighter
/// threshold of addGradient (conv * 1e-3); Psi is that of q_at. Claude Generated (Nov 2026).
bool DdCosmoSolvation::segments(const Vector& q_at, Segments& out) const
{
    if (!m_initialized || q_at.size() != m_nat)
        return false;
    ensureSolution(q_at);
    solveAdjoint(m_have_s, m_conv * 1.0e-3);
    m_q_s = q_at;
    const double autoaa = ALPBParameters::autoaa;
    const Eigen::VectorXd phi = m_jmat * q_at;
    out.atom.resize(m_ncav);
    out.position = m_ccav;
    out.charge.resize(m_ncav);
    out.area.resize(m_ncav);
    out.potential = phi;
    out.radius.resize(m_nat);
    out.atom_position = m_xyz;
    int ii = 0;
    for (int iat = 0; iat < m_nat; ++iat) {
        out.radius(iat) = m_rvdw(iat) * autoaa;
        for (int ig = 0; ig < m_ngrid; ++ig)
            if (m_ui(ig, iat) > 0.0) {
                out.atom[ii] = iat;
                out.charge(ii) = m_w(ig) * m_ui(ig, iat) * m_basis.col(ig).dot(m_s.col(iat));
                out.area(ii) = m_w(ig) * autoaa * autoaa * m_rvdw(iat) * m_rvdw(iat);
                ++ii;
            }
    }
    return ii == m_ncav;
}

void DdCosmoSolvation::addPotential(const Vector& q_at, Vector& v_at) const
{
    if (!m_initialized)
        return;
    // cosmo.f90 addShift: direct solve, adjoint solve, then v = -J^T zeta + keps sqrt(4 pi) X_00
    ensureSolution(q_at);
    if (!(m_have_s && m_q_s.size() == q_at.size() && m_q_s == q_at)) {
        solveAdjoint(m_have_s, m_conv);
        m_q_s = q_at;
    }
    Eigen::VectorXd zeta;
    getZeta(zeta);
    v_at.noalias() -= m_jmat.transpose() * zeta;
    for (int iat = 0; iat < m_nat; ++iat)
        v_at(iat) += m_keps * std::sqrt(fourpi) * m_sigma(0, iat);
}

double DdCosmoSolvation::electrostaticEnergy(const Vector& q_at) const
{
    if (!m_initialized)
        return 0.0;
    ensureSolution(q_at);
    return m_keps * (m_sigma.array() * m_psi.array()).sum();
}

/// cosmo.f90 getEnergy: E = 1/2 (1 - 1/eps) <X, Psi> + E_sasa.
double DdCosmoSolvation::energy(const Vector& q_at) const
{
    if (!m_initialized)
        return 0.0;
    return electrostaticEnergy(q_at) + m_gsasa;
}

/// Force kernels of <S, L^(x) X> and <S, g^(x)> (cosmo.f90 forces).
void DdCosmoSolvation::forces(const Eigen::VectorXd& phi, Eigen::MatrixXd& fx) const
{
    const Eigen::MatrixXd xi = m_basis.transpose() * m_s;   // xi(ig, iat)
    Eigen::MatrixXd phiexp = Eigen::MatrixXd::Zero(m_ngrid, m_nat);
    int ii = 0;
    for (int iat = 0; iat < m_nat; ++iat)
        for (int ig = 0; ig < m_ngrid; ++ig)
            if (m_ui(ig, iat) > 0.0)
                phiexp(ig, iat) = phi(ii++);
    fx = Eigen::MatrixXd::Zero(3, m_nat);
    for (int iat = 0; iat < m_nat; ++iat) {
        fdoka(iat, xi, fx.col(iat).data());
        fdokb(iat, xi, fx.col(iat).data());
        fdoga(iat, xi, phiexp, fx.col(iat).data());
    }
    fx *= m_keps;
}

/// First part of <S, L^(x) X> (core.f90 fdoka).
void DdCosmoSolvation::fdoka(int iat, const Eigen::MatrixXd& xi, double* fx) const
{
    std::vector<double> basloc(m_nylm), dbsloc(3 * m_nylm), vplm(m_nylm), vcos(m_lmax + 1), vsin(m_lmax + 1);
    const double tlow = 1.0 - 0.5 * (1.0 - kSe) * m_eta;
    const double thigh = 1.0 + 0.5 * (1.0 + kSe) * m_eta;
    for (int ig = 0; ig < m_ngrid; ++ig) {
        double va[3] = { 0.0, 0.0, 0.0 };
        for (int ij = m_inl[iat]; ij < m_inl[iat + 1]; ++ij) {
            const int jat = m_nl[ij];
            const Eigen::Vector3d vij = m_xyz.col(iat) + m_rvdw(iat) * m_grid.col(ig) - m_xyz.col(jat);
            const double vvij = vij.norm();
            const double tij = vvij / m_rvdw(jat);
            if (tij >= thigh)
                continue;
            const Eigen::Vector3d sij = vij / vvij;
            dbasis(sij.data(), basloc.data(), dbsloc.data(), vplm.data(), vcos.data(), vsin.data());
            double alp[3] = { 0.0, 0.0, 0.0 };
            double t = 1.0;
            for (int l = 1; l <= m_lmax; ++l) {
                const int i0 = ind(l);
                const double fl = l;
                const double fac = t / m_facl[i0];
                for (int m = -l; m <= l; ++m) {
                    const double f2 = fac * m_sigma(i0 + m, jat);
                    const double f1 = f2 * fl * basloc[i0 + m];
                    for (int c = 0; c < 3; ++c)
                        alp[c] += f1 * sij(c) + f2 * dbsloc[3 * (i0 + m) + c];
                }
                t *= tij;
            }
            const double beta = intmlp(tij, m_sigma.col(jat).data(), basloc.data());
            const double xij = fsw(tij, m_eta);
            double oij, f2;
            if (m_fi(ig, iat) > 1.0) {
                oij = xij / m_fi(ig, iat);
                f2 = -oij / m_fi(ig, iat);
            } else {
                oij = xij;
                f2 = 0.0;
            }
            const double f1 = oij / m_rvdw(jat);
            for (int c = 0; c < 3; ++c)
                va[c] += f1 * alp[c] + beta * f2 * zi(c, ig, iat);
            if (tij > tlow) {
                double f3 = beta * dfsw(tij, m_eta) / m_rvdw(jat);
                if (m_fi(ig, iat) > 1.0)
                    f3 /= m_fi(ig, iat);
                for (int c = 0; c < 3; ++c)
                    va[c] += f3 * sij(c);
            }
        }
        for (int c = 0; c < 3; ++c)
            fx[c] -= m_w(ig) * xi(ig, iat) * va[c];
    }
}

/// Second part of <S, L^(x) X> (core.f90 fdokb).
void DdCosmoSolvation::fdokb(int iat, const Eigen::MatrixXd& xi, double* fx) const
{
    std::vector<double> basloc(m_nylm), dbsloc(3 * m_nylm), vplm(m_nylm), vcos(m_lmax + 1), vsin(m_lmax + 1);
    const double tlow = 1.0 - 0.5 * (1.0 - kSe) * m_eta;
    const double thigh = 1.0 + 0.5 * (1.0 + kSe) * m_eta;
    for (int ig = 0; ig < m_ngrid; ++ig) {
        double vb[3] = { 0.0, 0.0, 0.0 }, vc[3] = { 0.0, 0.0, 0.0 };
        for (int ji = m_inl[iat]; ji < m_inl[iat + 1]; ++ji) {
            const int jat = m_nl[ji];
            const Eigen::Vector3d vji = m_xyz.col(jat) + m_rvdw(jat) * m_grid.col(ig) - m_xyz.col(iat);
            const double vvji = vji.norm();
            const double tji = vvji / m_rvdw(iat);
            if (tji > thigh)
                continue;
            const Eigen::Vector3d sji = vji / vvji;
            dbasis(sji.data(), basloc.data(), dbsloc.data(), vplm.data(), vcos.data(), vsin.data());
            double alp[3] = { 0.0, 0.0, 0.0 };
            double t = 1.0;
            for (int l = 1; l <= m_lmax; ++l) {
                const int i0 = ind(l);
                const double fl = l;
                const double fac = t / m_facl[i0];
                for (int m = -l; m <= l; ++m) {
                    const double f2 = fac * m_sigma(i0 + m, iat);
                    const double f1 = f2 * fl * basloc[i0 + m];
                    for (int c = 0; c < 3; ++c)
                        alp[c] += f1 * sji(c) + f2 * dbsloc[3 * (i0 + m) + c];
                }
                t *= tji;
            }
            const double xji = fsw(tji, m_eta);
            const double oji = (m_fi(ig, jat) > 1.0) ? xji / m_fi(ig, jat) : xji;
            const double f1 = oji / m_rvdw(iat);
            for (int c = 0; c < 3; ++c)
                vb[c] += f1 * alp[c] * xi(ig, jat);
            if (tji > tlow) {
                const double beta = intmlp(tji, m_sigma.col(iat).data(), basloc.data());
                double di, fac;
                if (m_fi(ig, jat) > 1.0) {
                    di = 1.0 / m_fi(ig, jat);
                    fac = di * xji;
                    bool proc = false;
                    double b = 0.0;
                    for (int jk = m_inl[jat]; jk < m_inl[jat + 1]; ++jk) {
                        const int kat = m_nl[jk];
                        const Eigen::Vector3d vjk = m_xyz.col(jat) + m_rvdw(jat) * m_grid.col(ig) - m_xyz.col(kat);
                        const double vvjk = vjk.norm();
                        const double tjk = vvjk / m_rvdw(kat);
                        if (kat != iat && tjk <= thigh) {
                            proc = true;
                            const Eigen::Vector3d sjk = vjk / vvjk;
                            ylmbas(sjk.data(), basloc.data(), vplm.data(), vcos.data(), vsin.data());
                            const double g1 = intmlp(tjk, m_sigma.col(kat).data(), basloc.data());
                            const double xjk = fsw(tjk, m_eta);
                            b += g1 * xjk;
                        }
                    }
                    if (proc) {
                        const double g1 = di * di * dfsw(tji, m_eta) / m_rvdw(iat);
                        const double g2 = g1 * xi(ig, jat) * b;
                        for (int c = 0; c < 3; ++c)
                            vc[c] += g2 * sji(c);
                    }
                } else {
                    di = 1.0;
                    fac = 0.0;
                }
                const double f2 = (1.0 - fac) * di * dfsw(tji, m_eta) / m_rvdw(iat);
                for (int c = 0; c < 3; ++c)
                    vb[c] += f2 * xi(ig, jat) * beta * sji(c);
            }
        }
        for (int c = 0; c < 3; ++c)
            fx[c] += m_w(ig) * (vb[c] - vc[c]);
    }
}

/// U^(x) Phi contribution to <S, g^(x)> (core.f90 fdoga).
void DdCosmoSolvation::fdoga(int iat, const Eigen::MatrixXd& xi, const Eigen::MatrixXd& phi, double* fx) const
{
    const double swthr = 1.0 + (kSe + 1.0) * m_eta / 2.0;
    for (int ig = 0; ig < m_ngrid; ++ig) {
        double alp[3] = { 0.0, 0.0, 0.0 };
        if (m_ui(ig, iat) > 0.0 && m_ui(ig, iat) < 1.0)
            for (int c = 0; c < 3; ++c)
                alp[c] += phi(ig, iat) * xi(ig, iat) * zi(c, ig, iat);
        for (int ji = m_inl[iat]; ji < m_inl[iat + 1]; ++ji) {
            const int jat = m_nl[ji];
            const Eigen::Vector3d vji = m_xyz.col(jat) + m_rvdw(jat) * m_grid.col(ig) - m_xyz.col(iat);
            const double vvji = vji.norm();
            const double tji = vvji / m_rvdw(iat);
            if (tji < swthr && tji > swthr - m_eta && m_ui(ig, jat) > 0.0) {
                const Eigen::Vector3d sji = vji / vvji;
                const double fac = -dfsw(tji, m_eta) / m_rvdw(iat);
                for (int c = 0; c < 3; ++c)
                    alp[c] += fac * phi(ig, jat) * xi(ig, jat) * sji(c);
            }
        }
        for (int c = 0; c < 3; ++c)
            fx[c] -= m_w(ig) * alp[c];
    }
}

/**
 * @brief Nuclear gradient of E_el + E_sasa at the converged charges (cosmo.f90 addGradient).
 *
 * dE/dR = dE_sasa/dR - fx, with fx = keps [<S, L^(x) X> + <S, U^(x) Phi>] plus the two
 * field terms: zeta times the solute field at the cavity points and the field of zeta at the nuclei
 * times q. The adjoint S is re-converged to conv * 1e-3 first, as in xtb.
 */
void DdCosmoSolvation::addGradient(const std::vector<int>& atomic_numbers, const Matrix& xyz_bohr,
                                   const Vector& q_at, Matrix& gradient)
{
    (void)atomic_numbers;
    (void)xyz_bohr;
    if (!m_initialized)
        return;
    for (int j = 0; j < m_nat; ++j)
        for (int c = 0; c < 3; ++c)
            gradient(j, c) += m_dsdr(c, j);

    ensureSolution(q_at);
    solveAdjoint(m_have_s && m_q_s.size() == q_at.size() && m_q_s == q_at, m_conv * 1.0e-3);
    m_q_s = q_at;

    const Eigen::VectorXd phi = m_jmat * q_at;
    Eigen::MatrixXd fx;
    forces(phi, fx);
    Eigen::VectorXd zeta;
    getZeta(zeta);

    // 1. solute field at the cavity points times zeta
    int ii = 0;
    for (int iat = 0; iat < m_nat; ++iat)
        for (int ig = 0; ig < m_ngrid; ++ig)
            if (m_ui(ig, iat) > 0.0) {
                Eigen::Vector3d ef = Eigen::Vector3d::Zero();
                for (int k = 0; k < m_nat; ++k) {
                    const Eigen::Vector3d v = m_ccav.col(ii) - m_xyz.col(k);
                    const double r2 = v.squaredNorm();
                    ef += q_at(k) / (r2 * std::sqrt(r2)) * v;
                }
                fx.col(iat) -= zeta(ii) * ef;
                ++ii;
            }
    // 2. field of zeta at the nuclei times the charges
    for (int iat = 0; iat < m_nat; ++iat) {
        Eigen::Vector3d ef = Eigen::Vector3d::Zero();
        for (int ic = 0; ic < m_ncav; ++ic) {
            const Eigen::Vector3d v = m_xyz.col(iat) - m_ccav.col(ic);
            const double r2 = v.squaredNorm();
            ef += zeta(ic) / (r2 * std::sqrt(r2)) * v;
        }
        fx.col(iat) -= ef * q_at(iat);
    }
    for (int j = 0; j < m_nat; ++j)
        for (int c = 0; c < 3; ++c)
            gradient(j, c) -= fx(c, j);
}
