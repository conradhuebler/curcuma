/*
 * <CPCM-X: SMD cavity-dispersion-solvent-structure (CDS) term>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Claude Generated (Nov 2026): port of the CDS path that CPCM-X 1.1.0 (src/cpcmx/smd.f90 calculate_cds_internal)
 * runs through numsa 0.2.0 (LGPL-3.0-or-later):
 *   src/smd/io.f90       read_parameters_h2o / read_parameters_ot (the parameter text, kept in memory here; CPCM-X
 *                        writes it to ./smd_h2o or ./smd_ot and numsa reads that file back, curcuma writes no file)
 *   src/smd/init.f90     init_smd_h2o / init_smd_ot (solvent-descriptor combination of the parameters)
 *   src/smd/sigma.f90    calc_surft, T (atomic surface tensions with the SMD switching function)
 *   src/smd/cds.f90      calc_cds
 *   src/numsa/surface.f90 new_surface_integrator, update_nnlist, compute_numsa, compute_w_sp (smoothed SASA)
 *   src/numsa/data.f90   vdw_rad_smd (Bondi radii of SMD)
 * Model: G_CDS = sum_k sigma_k A_k + sigma^[M] sum_k A_k (Marenich, Cramer, Truhlar, J. Phys. Chem. B 113, 6378
 * (2009)), A_k the solvent accessible surface of atom k, sigma_k = Z_k + sum_k' Z_kk' T(R_kk') from the geometry.
 * AI-generated, machine-tested only; human production testing pending.
 */

#include "cpcmx.h"
#include "cpcmx_internal.h"

#include "lebedev_xtb_grids.h"
#include "src/core/elements.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace curcuma::solvation::cpcmx {

namespace {

constexpr int kMaxElem = 94;

/// numsa smd_init smd_param (1-based element indices as in Fortran; index 0 unused).
struct SmdParam {
    std::vector<double> zk = std::vector<double>(kMaxElem + 1, 0.0);
    std::vector<std::vector<double>> zkk = std::vector<std::vector<double>>(kMaxElem + 1, std::vector<double>(kMaxElem + 1, 0.0));
    std::vector<std::vector<double>> rzkk = zkk, drzkk = zkk;
    double nc3 = 0.0, rnc3 = 0.0, drnc3 = 0.0, s_m = 0.0, alpha = 0.0;
};

int toNumber(const std::string& symbol) { return Elements::String2Element(symbol); }

std::vector<std::string> tokens(const std::string& line)
{
    std::istringstream ss(line);
    std::vector<std::string> t;
    for (std::string w; ss >> w;)
        t.push_back(w);
    return t;
}

/// io.f90: read "Symbol value" lines into zk(row, Z) until `stop` (the first line of the text is skipped by the caller).
size_t readZk(const std::vector<std::string>& t, size_t k, const std::string& stop, std::vector<double>& zk)
{
    for (; k < t.size() && t[k] != stop; ++k) {
        const auto w = tokens(t[k]);
        zk[toNumber(w[0])] = std::stod(w[1]);
    }
    return k + 1;   // skip the stop line
}

/// io.f90: a square "Symbol v1 v2 ..." block until `stop`: m(sym_i, sym_j) = value(i, j).
size_t readMatrix(const std::vector<std::string>& t, size_t k, const std::string& stop, std::vector<std::vector<double>>& m)
{
    std::vector<std::vector<std::string>> rows;
    for (; k < t.size() && t[k] != stop; ++k)
        rows.push_back(tokens(t[k]));
    const size_t n = rows.size();
    std::vector<int> sym(n);
    for (size_t i = 0; i < n; ++i)
        sym[i] = toNumber(rows[i][0]);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            m[sym[i]][sym[j]] = std::stod(rows[i][1 + j]);
    return k + 1;
}

/// init_smd_h2o with read_parameters_h2o of the internal water text (alpha = 0.82, s_m = 0).
SmdParam smdWater()
{
    SmdParam p;
    const auto& t = smdTextWater();
    size_t k = 1;
    k = readZk(t, k, "#Zkk Matrix", p.zk);
    k = readMatrix(t, k, "#rzkk Matrix", p.zkk);
    k = readMatrix(t, k, "#drzkk Matrix", p.rzkk);
    k = readMatrix(t, k, "#NC3", p.drzkk);
    p.nc3 = std::stod(t[k]);
    p.rnc3 = std::stod(t[k + 1]);
    p.drnc3 = std::stod(t[k + 2]);
    p.alpha = 0.82;
    p.s_m = 0.0;
    return p;
}

/// init_smd_ot with read_parameters_ot of the internal text: zk = zk1 n + zk2 alpha + zk3 beta (same for zkk) and
/// s_m = s_g gamma + s_r2 arom^2 + s_p2 fclbr^2 + s_b2 beta^2 from the solvent descriptors of init_smd_def.
bool smdOther(const std::string& solvent, SmdParam& p)
{
    double n, alpha, beta, msurft, arom, fclbr;
    if (!smdSolventProperties(solvent, n, alpha, beta, msurft, arom, fclbr))
        return false;
    std::vector<std::vector<double>> zk(3, std::vector<double>(kMaxElem + 1, 0.0));
    std::vector<std::vector<std::vector<double>>> zkk(3, p.zkk);
    const auto& t = smdTextOther();
    size_t k = 1;
    k = readZk(t, k, "#Zk2", zk[0]);
    k = readZk(t, k, "#Zk3", zk[1]);
    k = readZk(t, k, "#Zkk Matrix", zk[2]);
    k = readMatrix(t, k, "#Zkk2 Matrix", zkk[0]);
    k = readMatrix(t, k, "#Zkk3 Matrix", zkk[1]);
    k = readMatrix(t, k, "#rzkk Matrix", zkk[2]);
    k = readMatrix(t, k, "#drzkk Matrix", p.rzkk);
    k = readMatrix(t, k, "#NC3", p.drzkk);
    p.nc3 = std::stod(t[k]);
    p.rnc3 = std::stod(t[k + 1]);
    p.drnc3 = std::stod(t[k + 2]);
    const double sg = std::stod(t[k + 4]), sr = std::stod(t[k + 5]), sp = std::stod(t[k + 6]), sb = std::stod(t[k + 7]);
    for (int z = 1; z <= kMaxElem; ++z) {
        p.zk[z] = zk[0][z] * n + zk[1][z] * alpha + zk[2][z] * beta;
        for (int z2 = 1; z2 <= kMaxElem; ++z2)
            p.zkk[z][z2] = zkk[0][z][z2] * n + zkk[1][z][z2] * alpha + zkk[2][z][z2] * beta;
    }
    p.s_m = sg * msurft + sr * (arom * arom) + sp * (fclbr * fclbr) + sb * (beta * beta);
    p.alpha = alpha;
    return true;
}

/// numsa data.f90 vdw_rad_smd (Angstrom) for Z = 1..88, converted with mctc aatoau; -1 (missing) beyond.
double smdRadius(int z)
{
    static const double s = 2.0;   // smdstub
    static const double r[88] = {
        1.20, 1.40, 1.81, 1.53, 1.92, 1.85, 1.89, 1.52, 1.73, 1.54, 2.27, 1.73, 1.84, 2.47, 2.12, 2.49,
        2.38, 1.88, 2.75, 2.31, s,    s,    s,    s,    s,    s,    s,    s,    s,    s,    1.87, 2.11,
        1.85, 1.90, 3.06, 2.02, 3.03, 2.49, s,    s,    s,    s,    s,    s,    s,    s,    s,    s,
        1.93, 2.17, 2.06, 2.06, 1.98, 2.16, 3.43, 2.68, s,    s,    s,    s,    s,    s,    s,    s,
        s,    s,    s,    s,    s,    s,    s,    s,    s,    s,    s,    s,    s,    s,    s,    s,
        1.96, 2.02, 2.07, 1.97, 2.02, 2.20, 3.48, 2.83 };
    return (z >= 1 && z <= 88) ? mctc_aatoau * r[z - 1] : -1.0;
}

/// sigma.f90 T: SMD switching function exp(dr / (R - dr - r)) for R < r + dr (R in Angstrom), else 0.
double switchT(const Eigen::MatrixXd& xyz, int a, int b, double rzkk, double drzkk)
{
    const double dx = xyz(0, a) - xyz(0, b), dy = xyz(1, a) - xyz(1, b), dz = xyz(2, a) - xyz(2, b);
    const double R = std::sqrt(dx * dx + dy * dy + dz * dz) * mctc_autoaa;
    if (R < (rzkk + drzkk)) {
        const double denom = R - drzkk - rzkk;
        return std::exp(drzkk / denom);
    }
    return 0.0;
}

/// sigma.f90 calc_surft: atomic surface tensions sigma_k (cal/mol/A^2) for H, C, N, O; Z_k for every other element.
std::vector<double> calcSurft(const Eigen::MatrixXd& xyz, const std::vector<int>& Z, const SmdParam& p)
{
    const int nat = static_cast<int>(Z.size());
    std::vector<double> sk(nat, 0.0);
    for (int i = 0; i < nat; ++i) {
        double s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0;
        switch (Z[i]) {
        case 1:
            for (int j = 0; j < nat; ++j) {
                if (Z[j] == 6)
                    s1 = s1 + switchT(xyz, i, j, p.rzkk[1][6], p.drzkk[1][6]);
                else if (Z[j] == 8)
                    s2 = s2 + switchT(xyz, i, j, p.rzkk[1][8], p.drzkk[1][8]);
            }
            sk[i] = p.zk[1] + p.zkk[1][6] * s1 + p.zkk[1][8] * s2;
            break;
        case 6:
            for (int j = 0; j < nat; ++j) {
                if (Z[j] == 6) {
                    if (i != j)
                        s1 = s1 + switchT(xyz, i, j, p.rzkk[6][6], p.drzkk[6][6]);
                } else if (Z[j] == 7) {
                    s2 = s2 + switchT(xyz, i, j, p.rzkk[6][7], p.drzkk[6][7]);
                }
            }
            sk[i] = p.zk[6] + p.zkk[6][6] * s1 + p.zkk[6][7] * (s2 * s2);
            break;
        case 7:
            for (int j = 0; j < nat; ++j) {
                if (Z[j] != 6)
                    continue;
                double add = switchT(xyz, i, j, p.rzkk[6][6], p.drzkk[6][6]);
                double nc = 0.0;
                for (int k = 0; k < nat; ++k)   // all interactions of that C with the other atoms
                    if (k != j && k != i)
                        nc = nc + switchT(xyz, j, k, p.rzkk[6][Z[k]], p.drzkk[6][Z[k]]);
                add = add * (nc * nc);
                s1 = s1 + add;
                s2 = s2 + switchT(xyz, i, j, p.rnc3, p.drnc3);   // additional NC3 correction
            }
            sk[i] = p.zk[7] + p.zkk[7][6] * std::pow(s1, 1.3) + p.nc3 * s2;
            break;
        case 8:
            for (int j = 0; j < nat; ++j) {
                if (Z[j] == 6)
                    s1 = s1 + switchT(xyz, i, j, p.rzkk[8][6], p.drzkk[8][6]);
                else if (Z[j] == 7)
                    s2 = s2 + switchT(xyz, i, j, p.rzkk[8][7], p.drzkk[8][7]);
                else if (Z[j] == 8) {
                    if (i != j)
                        s3 = s3 + switchT(xyz, i, j, p.rzkk[8][8], p.drzkk[8][8]);
                } else if (Z[j] == 15)
                    s4 = s4 + switchT(xyz, i, j, p.rzkk[8][15], p.drzkk[8][15]);
            }
            sk[i] = p.zk[8] + p.zkk[8][6] * s1 + p.zkk[8][7] * s2 + p.zkk[8][8] * s3 + p.zkk[8][15] * s4;
            break;
        default:
            sk[i] = (Z[i] >= 1 && Z[i] <= kMaxElem) ? p.zk[Z[i]] : 0.0;
        }
    }
    return sk;
}

/// numsa surface.f90: smoothed numerical SASA per atom (Bohr^2) on a Lebedev grid of 110 points.
std::vector<double> numsaSurface(const Eigen::MatrixXd& xyz, const std::vector<double>& rad_of_atom, double probe)
{
    const int nat = static_cast<int>(rad_of_atom.size());
    // new_surface_integrator: the local w (default_smoothing = 0.3 * aatoau, 0.3 a single-precision literal)
    const double w = static_cast<double>(0.3f) * mctc_aatoau;
    const double w3 = w * (w * w);
    const double ah0 = 0.5, ah1 = 3.0 / (4.0 * w), ah3 = -1.0 / (4.0 * w3);
    std::vector<double> vdwsa(nat), wrp(nat), trj2a(nat), trj2b(nat);
    double vmax = -1e300;
    for (int i = 0; i < nat; ++i) {
        vdwsa[i] = rad_of_atom[i] + probe;
        trj2a[i] = (vdwsa[i] - w) * (vdwsa[i] - w);
        trj2b[i] = (vdwsa[i] + w) * (vdwsa[i] + w);
        double r = vdwsa[i] + w;
        wrp[i] = (0.25 / w + 3.0 * ah3 * (0.2 * r * r - 0.5 * r * vdwsa[i] + vdwsa[i] * vdwsa[i] / 3.0)) * r * r * r;
        r = vdwsa[i] - w;
        wrp[i] = wrp[i] - (0.25 / w + 3.0 * ah3 * (0.2 * r * r - 0.5 * r * vdwsa[i] + vdwsa[i] * vdwsa[i] / 3.0)) * r * r * r;
        vmax = std::max(vmax, vdwsa[i]);
    }
    const double srcut = 2 * (w + vmax) + 2.0 * mctc_aatoau;   // default_offset
    // grid_size(8) = 110 points; numsa multiplies the weights (sum 1) by 4*pi as (w * 4) * pi
    std::vector<double> gx, gw;
    curcuma::solvation::lebedev::getAngGrid(7, gx, gw);
    const int ngrid = static_cast<int>(gw.size());
    for (double& v : gw)
        v = v * 4 * pi;

    // update_nnlist: pairs (i, j < i) in this order; neighbours in the order the pairs are met
    std::vector<std::vector<int>> nn(nat);
    const double srcut2 = srcut * srcut;
    for (int i1 = 0; i1 < nat; ++i1)
        for (int i2 = 0; i2 < i1; ++i2) {
            const double x = xyz(0, i1) - xyz(0, i2), y = xyz(1, i1) - xyz(1, i2), z = xyz(2, i1) - xyz(2, i2);
            const double dr2 = x * x + y * y + z * z;
            if (dr2 < srcut2) {
                nn[i1].push_back(i2);
                nn[i2].push_back(i1);
            }
        }

    // compute_numsa / compute_w_sp
    std::vector<double> surface(nat, 0.0);
    for (int iat = 0; iat < nat; ++iat) {
        const double rsas = vdwsa[iat];
        double sasai = 0.0;
        for (int ip = 0; ip < ngrid; ++ip) {
            const double px = xyz(0, iat) + rsas * gx[3 * ip], py = xyz(1, iat) + rsas * gx[3 * ip + 1],
                         pz = xyz(2, iat) + rsas * gx[3 * ip + 2];
            double sasap = 1.0;
            for (int ia : nn[iat]) {
                const double tx = px - xyz(0, ia), ty = py - xyz(1, ia), tz = pz - xyz(2, ia);
                const double tj2 = tx * tx + ty * ty + tz * tz;
                if (tj2 < trj2b[ia]) {
                    if (tj2 <= trj2a[ia]) {
                        sasap = 0.0;
                        break;
                    }
                    const double sqtj = std::sqrt(tj2);
                    const double uj = sqtj - vdwsa[ia];
                    const double ah3uj2 = ah3 * uj * uj;
                    const double sasaij = ah0 + (ah1 + ah3uj2) * uj;
                    sasap = sasap * sasaij;
                }
            }
            if (sasap > 1.0e-6) {   // tolsesp
                const double wsa = gw[ip] * wrp[iat] * sasap;
                sasai = sasai + wsa;
            }
        }
        surface[iat] = sasai;
    }
    return surface;
}

} // namespace

/// smd.f90 calculate_cds_internal as CPCM-X calls it (cpxcalc.f90 cds_default): species = unique(segment atoms),
/// coordinates coord(j, :) * (1/BtoA) for j = 1..nspecies (the first nspecies atoms, as in the source), SMD radii
/// per atom with the oxygen rule (OxygenRadiusRule, cpcmx.h), probe radius added to the radii.
bool parseOxygenRadiusRule(const std::string& text, OxygenRadiusRule& rule)
{
    std::string t = text;
    std::transform(t.begin(), t.end(), t.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (t == "atom8" || t == "xtb")
        rule = OxygenRadiusRule::Atom8;
    else if (t == "oxygen" || t == "cpcmx")
        rule = OxygenRadiusRule::Oxygen;
    else if (t == "none")
        rule = OxygenRadiusRule::None;
    else
        return false;
    return true;
}

bool smdCds(const Molecule& solute, const std::string& solvent_name, double probe, OxygenRadiusRule rule, double& dG_cds,
            std::string& error)
{
    std::vector<int> species;   // sort.f90 unique_int: order of first appearance
    for (int id : solute.id)
        if (std::find(species.begin(), species.end(), id) == species.end())
            species.push_back(id);
    const int nsp = static_cast<int>(species.size());

    SmdParam param;
    if (solvent_name == "water" || solvent_name == "h2o") {
        param = smdWater();
    } else if (!smdOther(solvent_name, param)) {
        error = "CPCM-X: no SMD solvent properties for " + solvent_name + " (numsa init_smd stops here as well)";
        return false;
    }

    // rad = get_vdw_rad_smd(symbols), per atom; then the oxygen rule (OxygenRadiusRule); 0.43 and 1.8 are
    // single-precision literals in the source
    const int nat = static_cast<int>(solute.z.size());
    std::vector<double> rad(nat);
    for (int k = 0; k < nat; ++k)
        rad[k] = smdRadius(solute.z[k]);
    if (param.alpha < static_cast<double>(0.43f)) {
        const double add = static_cast<double>(1.8f) * (static_cast<double>(0.43f) - param.alpha);
        if (rule == OxygenRadiusRule::Oxygen) {
            for (int k = 0; k < nat; ++k)
                if (solute.z[k] == 8)
                    rad[k] = rad[k] + add;
        } else if (rule == OxygenRadiusRule::Atom8 && nat >= 8) {
            rad[7] = rad[7] + add;   // xtb 6.7.1: rad(8), the 8th atom
        }
    }

    Eigen::MatrixXd coord(3, nsp);
    std::vector<double> rad_sp(nsp);
    std::vector<int> Z(nsp);
    for (int j = 0; j < nsp; ++j) {
        for (int c = 0; c < 3; ++c)
            coord(c, j) = solute.atom_xyz(j, c) * (1 / btoa);
        rad_sp[j] = rad[species[j] - 1];   // new_surface_integrator: vdwsa(iat) = rad(num(iat)) + probe
        Z[j] = solute.z[species[j] - 1];   // init_info: Z = to_number(ident(species(i)))
    }

    const std::vector<double> surface = numsaSurface(coord, rad_sp, probe);
    const std::vector<double> sk = calcSurft(coord, Z, param);

    // cds.f90 calc_cds
    double sum_surface = 0.0, sum_cds = 0.0;
    for (int i = 0; i < nsp; ++i)
        sum_surface += surface[i];
    const double cds_sm = param.s_m * sum_surface * (mctc_autoaa * mctc_autoaa);
    for (int i = 0; i < nsp; ++i)
        sum_cds += sk[i] * surface[i] * (mctc_autoaa * mctc_autoaa);
    const double dG_disp = (sum_cds + cds_sm) / 1000;
    dG_cds = dG_disp / autokcal;
    return true;
}

} // namespace curcuma::solvation::cpcmx
