/*
 * <ALPB Solvation Model for GFN-FF>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): Port of ALPB solvation from Fortran gfnff/gbsa
 *
 * References:
 * - Fortran: external/gfnff/src/gbsa/gbsa.f90 (TBorn type, update, getEnergy, addGradient)
 * - Fortran: external/gfnff/src/gbsa/model.f90 (paramToModel, newBornModel)
 */

#include "alpb_solvation.h"
#include "alpb_parameters.h"
#include "src/core/curcuma_logger.h"
#include "src/core/units.h"

#include <algorithm>
#include <cmath>
#include <numeric>

using namespace ALPBParameters;

// ═══════════════════════════════════════════════════════════════════
// Lebedev 230-pt grid generator (ported from lebedev.f90)
// ═══════════════════════════════════════════════════════════════════

namespace {

/// genOh1: 6 points from (±1,0,0) permutations
void genOh1(int& n, double* x, double* w, double v) {
    double a = 1.0;
    // x is 3×nang layout: x[0..2] = first point xyz
    x[0]= a; x[1]= 0; x[2]= 0; w[0]=v;
    x[3]=-a; x[4]= 0; x[5]= 0; w[1]=v;
    x[6]= 0; x[7]= a; x[8]= 0; w[2]=v;
    x[9]= 0; x[10]=-a; x[11]=0; w[3]=v;
    x[12]=0; x[13]=0; x[14]= a; w[4]=v;
    x[15]=0; x[16]=0; x[17]=-a; w[5]=v;
    n += 6;
}

/// genOh3: 8 points from (±a,±a,±a), a=1/sqrt(3)
void genOh3(int& n, double* x, double* w, double v) {
    double a = std::sqrt(1.0/3.0);
    int k = 0;
    for (int s1 : {1, -1})
        for (int s2 : {1, -1})
            for (int s3 : {1, -1}) {
                x[k*3+0] = s1*a; x[k*3+1] = s2*a; x[k*3+2] = s3*a;
                w[k] = v;
                k++;
            }
    n += 8;
}

/// genOh4: 24 points from permutations of (±a,±a,±b), b=sqrt(1-2a²)
void genOh4(int& n, double* x, double* w, double a, double v) {
    double b = std::sqrt(1.0 - 2.0*a*a);
    double vals[3] = {a, a, b};
    int k = 0;
    // All permutations of (a,a,b) with sign changes
    for (int perm = 0; perm < 3; perm++) {
        double v0 = vals[perm]; double v1 = vals[(perm+1)%3]; double v2 = vals[(perm+2)%3];
        for (int s0 : {1,-1}) for (int s1 : {1,-1}) for (int s2 : {1,-1}) {
            x[k*3+0] = s0*v0; x[k*3+1] = s1*v1; x[k*3+2] = s2*v2;
            w[k] = v;
            k++;
        }
    }
    n += 24;
}

/// genOh5: 24 points from permutations of (±a,±b,0), b=sqrt(1-a²)
void genOh5(int& n, double* x, double* w, double a, double v) {
    double b = std::sqrt(1.0 - a*a);
    double vals[3][2] = {{a,b},{b,a},{a,b}};  // Will be placed in different plane pairs
    int k = 0;
    // (a,b,0) permutations
    for (int s0 : {1,-1}) for (int s1 : {1,-1}) {
        x[k*3+0]=s0*a; x[k*3+1]=s1*b; x[k*3+2]=0; w[k]=v; k++;
    }
    for (int s0 : {1,-1}) for (int s1 : {1,-1}) {
        x[k*3+0]=s0*b; x[k*3+1]=s1*a; x[k*3+2]=0; w[k]=v; k++;
    }
    for (int s0 : {1,-1}) for (int s1 : {1,-1}) {
        x[k*3+0]=s0*a; x[k*3+1]=0; x[k*3+2]=s1*b; w[k]=v; k++;
    }
    for (int s0 : {1,-1}) for (int s1 : {1,-1}) {
        x[k*3+0]=s0*b; x[k*3+1]=0; x[k*3+2]=s1*a; w[k]=v; k++;
    }
    for (int s0 : {1,-1}) for (int s1 : {1,-1}) {
        x[k*3+0]=0; x[k*3+1]=s0*a; x[k*3+2]=s1*b; w[k]=v; k++;
    }
    for (int s0 : {1,-1}) for (int s1 : {1,-1}) {
        x[k*3+0]=0; x[k*3+1]=s0*b; x[k*3+2]=s1*a; w[k]=v; k++;
    }
    n += 24;
}

/// genOh6: 48 points from permutations of (±a,±b,±c), c=sqrt(1-a²-b²)
void genOh6(int& n, double* x, double* w, double a, double b, double v) {
    double c = std::sqrt(1.0 - a*a - b*b);
    double vals[3] = {a, b, c};
    int k = 0;
    // All 6 permutations × 8 sign changes = 48
    int perms[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    for (auto& p : perms) {
        for (int s0 : {1,-1}) for (int s1 : {1,-1}) for (int s2 : {1,-1}) {
            x[k*3+0] = s0*vals[p[0]]; x[k*3+1] = s1*vals[p[1]]; x[k*3+2] = s2*vals[p[2]];
            w[k] = v;
            k++;
        }
    }
    n += 48;
}

/// Generate 230-point Lebedev grid (ld0230 from lebedev.f90)
void generateLebedev230(Eigen::MatrixXd& grid, Eigen::VectorXd& weights) {
    // Temporary flat storage: grid points as 3*230 doubles
    std::vector<double> x(3 * 230, 0.0);
    std::vector<double> w(230, 0.0);
    int n = 0;

    double a, b, v;

    // Exact sequence from ld0230
    v = -0.5522639919727325e-1;
    genOh1(n, &x[n*3], &w[n-0], v);  // n was 0, now 6

    // Actually, need to track offset differently:
    // Reset n to track current insertion point
    n = 0;
    v = -0.5522639919727325e-1;
    genOh1(n, x.data(), w.data(), v);

    v = 0.4450274607445226e-2;
    genOh3(n, x.data() + n*3, w.data() + n, v);

    a = 0.4492044687397611e+0;
    v = 0.4496841067921404e-2;
    genOh4(n, x.data() + n*3, w.data() + n, a, v);

    a = 0.2520419490210201e+0;
    v = 0.5049153450478750e-2;
    genOh4(n, x.data() + n*3, w.data() + n, a, v);

    a = 0.6981906658447242e+0;
    v = 0.3976408018051883e-2;
    genOh4(n, x.data() + n*3, w.data() + n, a, v);

    a = 0.6587405243460960e+0;
    v = 0.4401400650381014e-2;
    genOh4(n, x.data() + n*3, w.data() + n, a, v);

    a = 0.4038544050097660e-1;
    v = 0.1724544350544401e-1;
    genOh4(n, x.data() + n*3, w.data() + n, a, v);

    a = 0.5823842309715585e+0;
    v = 0.4231083095357343e-2;
    genOh5(n, x.data() + n*3, w.data() + n, a, v);

    a = 0.3545877390518688e+0;
    v = 0.5198069864064399e-2;
    genOh5(n, x.data() + n*3, w.data() + n, a, v);

    a = 0.2272181808998187e+0;
    b = 0.4864661535886647e+0;
    v = 0.4695720972568883e-2;
    genOh6(n, x.data() + n*3, w.data() + n, a, b, v);

    // Copy to Eigen
    grid.resize(3, n);
    weights.resize(n);
    for (int i = 0; i < n; ++i) {
        grid(0, i) = x[i*3+0];
        grid(1, i) = x[i*3+1];
        grid(2, i) = x[i*3+2];
        weights(i) = w[i];
    }
}

} // anonymous namespace


// ═══════════════════════════════════════════════════════════════════
// ALPBSolvation implementation
// ═══════════════════════════════════════════════════════════════════

ALPBSolvation::ALPBSolvation() = default;
ALPBSolvation::~ALPBSolvation() = default;

bool ALPBSolvation::init(const std::vector<int>& atomic_numbers, const std::string& solvent)
{
    m_solvent = solvent;
    m_nat = static_cast<int>(atomic_numbers.size());
    m_ntpair = m_nat * (m_nat - 1) / 2;

    // Load solvent parameters
    if (!loadSolventParameters(solvent)) {
        CurcumaLogger::error("ALPB: Unknown solvent '" + solvent + "'");
        return false;
    }

    // Setup per-atom radii and SASA parameters
    setupRadii(atomic_numbers);

    // Generate Lebedev angular grid
    generateLebedevGrid();

    // Setup pair indices
    m_ppind.resize(m_ntpair);
    int ij = 0;
    for (int i = 0; i < m_nat; ++i) {
        for (int j = 0; j < i; ++j) {
            m_ppind[ij] = {i, j};
            ij++;
        }
    }

    // Allocate arrays
    m_ddpair = Eigen::MatrixXd::Zero(4, m_ntpair);
    m_nnlistr.reserve(m_ntpair);
    m_nnsas.resize(m_nat, 0);
    m_nnlists.resize(m_nat);

    m_brad = Eigen::VectorXd::Zero(m_nat);
    m_brdr = Eigen::MatrixXd::Zero(3 * m_nat, m_nat);

    m_sasa = Eigen::VectorXd::Zero(m_nat);
    m_dsdr = Eigen::MatrixXd::Zero(3, m_nat);
    m_dsdrt.resize(m_nat);
    for (int i = 0; i < m_nat; ++i) {
        m_dsdrt[i] = Eigen::MatrixXd::Zero(3, m_nat);
    }

    // HB setup
    m_hbmag = Eigen::VectorXd::Zero(m_nat);
    for (int i = 0; i < m_nat; ++i) {
        int z = atomic_numbers[i];
        if (z >= 1 && z <= MAX_ELEM) {
            m_hbmag(i) = m_hb_elem[z - 1];
        }
    }
    m_lhb = (m_hbmag.array() < 0.0).any();
    if (m_lhb) {
        m_hbw = Eigen::VectorXd::Zero(m_nat);
        m_dhbdw = Eigen::VectorXd::Zero(m_nat);
    }

    m_born_mat = Eigen::MatrixXd::Zero(m_nat, m_nat);

    m_initialized = true;

    if (CurcumaLogger::get_verbosity() >= 1) {
        printInfo();
    }

    return true;
}

bool ALPBSolvation::loadSolventParameters(const std::string& solvent)
{
    const SolventParameters* p = getSolventParameters(solvent);
    if (!p) return false;

    m_dielectric_const = p->epsv;
    m_molar_mass = p->smass;
    m_density = p->rhos / p->smass * molcm3toau;
    m_gshift = p->gshift * kcaltoau;
    m_born_scale = p->c1;
    m_born_offset = p->soset * 0.1 * aatoau;
    m_probe_rad = p->rprobe * aatoau;

    // ALPB: alpbet = 0.571412 / eps (Fortran gbsa.f90:324-325)
    // For GFN-FF, param%alpha is always 0, but ALPB is always enabled (P16 kernel)
    m_alpbet = alpb_factor / m_dielectric_const;
    m_keps = (1.0 / m_dielectric_const - 1.0) / (1.0 + m_alpbet);

    // Born radii cutoff
    m_lrcut = lrcut;

    // Element-specific parameters
    m_gamscale_elem.resize(MAX_ELEM);
    m_sx_elem.resize(MAX_ELEM);
    m_hb_elem.resize(MAX_ELEM);

    for (int i = 0; i < MAX_ELEM; ++i) {
        // Surface tension: gamscale * 4π * 1e-5 (Fortran model.f90:520)
        m_gamscale_elem[i] = p->gamscale[i] * fourpi * surfaceTension;
        // Descreening
        m_sx_elem[i] = p->sx[i];
        // HB strength: -(tmp²) * kcaltoau (Fortran model.f90:524-525)
        if (std::abs(p->tmp[i]) > 1.0e-3) {
            m_hb_elem[i] = -(p->tmp[i] * p->tmp[i]) * kcaltoau;
        } else {
            m_hb_elem[i] = 0.0;
        }
    }

    return true;
}

void ALPBSolvation::setupRadii(const std::vector<int>& atomic_numbers)
{
    m_vdwr.resize(m_nat);
    m_rho.resize(m_nat);
    m_svdw.resize(m_nat);
    m_vdwsa.resize(m_nat);
    m_trj2.resize(2, m_nat);
    m_wrp.resize(m_nat);
    m_gamsasa.resize(m_nat);

    for (int i = 0; i < m_nat; ++i) {
        int z = atomic_numbers[i];
        m_vdwr(i) = getVdwRadD3(z);
        m_rho(i) = m_vdwr(i) * (z >= 1 && z <= MAX_ELEM ? m_sx_elem[z-1] : 0.8);
        m_svdw(i) = m_vdwr(i) - m_born_offset;

        m_vdwsa(i) = m_vdwr(i) + m_probe_rad;
        m_trj2(0, i) = (m_vdwsa(i) - w_smooth) * (m_vdwsa(i) - w_smooth);
        m_trj2(1, i) = (m_vdwsa(i) + w_smooth) * (m_vdwsa(i) + w_smooth);

        // Radial weight wrp (Fortran gbsa.f90:393-399)
        double r = m_vdwsa(i) + w_smooth;
        double wrp1 = (0.25 / w_smooth + 3.0 * ah3 * (0.2*r*r - 0.5*r*m_vdwsa(i) +
                       m_vdwsa(i)*m_vdwsa(i)/3.0)) * r*r*r;
        r = m_vdwsa(i) - w_smooth;
        double wrp2 = (0.25 / w_smooth + 3.0 * ah3 * (0.2*r*r - 0.5*r*m_vdwsa(i) +
                       m_vdwsa(i)*m_vdwsa(i)/3.0)) * r*r*r;
        m_wrp(i) = wrp1 - wrp2;

        m_gamsasa(i) = (z >= 1 && z <= MAX_ELEM) ? m_gamscale_elem[z-1] : 0.0;
    }

    // SASA cutoff (Fortran gbsa.f90:402)
    m_srcut = 2.0 * (w_smooth + m_vdwsa.maxCoeff()) + rOffset;
}

void ALPBSolvation::generateLebedevGrid()
{
    generateLebedev230(m_ang_grid, m_ang_weight);
    m_nang = static_cast<int>(m_ang_weight.size());
}


// ═══════════════════════════════════════════════════════════════════
// Update: Neighbor lists, Born radii, SASA, Born matrix
// Reference: Fortran gbsa.f90:426-473 (update subroutine)
// ═══════════════════════════════════════════════════════════════════

void ALPBSolvation::update(const std::vector<int>& atomic_numbers, const Matrix& xyz_bohr)
{
    // 1. Update neighbor lists and pair distances
    updateNeighborList(xyz_bohr);

    // 2. Compute Born radii
    computeBornRadii();

    // 3. Compute SASA
    computeSASA(xyz_bohr);

    // 4. Contract SASA gradient: dsdr = dsdrt * gamsasa
    m_dsdr = Eigen::MatrixXd::Zero(3, m_nat);
    for (int i = 0; i < m_nat; ++i) {
        for (int j = 0; j < m_nat; ++j) {
            m_dsdr.col(j) += m_dsdrt[i].col(j) * m_gamsasa(i);
        }
    }
    m_gsasa_total = m_sasa.dot(m_gamsasa);

    // 5. HB correction
    if (m_lhb) {
        computeHBCorrection();
    }

    // 6. ALPB shape descriptor
    if (m_alpbet > 0.0) {
        computeADet(xyz_bohr);
    }

    // 7. Build Born interaction matrix
    m_born_mat.setZero();
    buildBornMatrix();
}


// ═══════════════════════════════════════════════════════════════════
// Neighbor list update (Fortran gbsa.f90:945-997)
// ═══════════════════════════════════════════════════════════════════

void ALPBSolvation::updateNeighborList(const Matrix& xyz_bohr)
{
    double lrcut2 = m_lrcut * m_lrcut;
    double srcut2 = m_srcut * m_srcut;

    m_nnrad = 0;
    m_nnlistr.clear();
    for (int i = 0; i < m_nat; ++i) {
        m_nnsas[i] = 0;
        m_nnlists[i].clear();
    }

    // Compute pair distances and build neighbor lists
    for (int kk = 0; kk < m_ntpair; ++kk) {
        int i = m_ppind[kk].first;
        int j = m_ppind[kk].second;

        // xyz_bohr is N×3 row-major
        double dx = xyz_bohr(i, 0) - xyz_bohr(j, 0);
        double dy = xyz_bohr(i, 1) - xyz_bohr(j, 1);
        double dz = xyz_bohr(i, 2) - xyz_bohr(j, 2);
        double dr2 = dx*dx + dy*dy + dz*dz;
        double r = std::sqrt(dr2);

        m_ddpair(0, kk) = r;
        m_ddpair(1, kk) = dx;
        m_ddpair(2, kk) = dy;
        m_ddpair(3, kk) = dz;

        if (dr2 < lrcut2) {
            m_nnlistr.push_back({i, j, kk});
            m_nnrad++;

            if (dr2 < srcut2) {
                m_nnsas[i]++;
                m_nnsas[j]++;
                m_nnlists[i].push_back(j);
                m_nnlists[j].push_back(i);
            }
        }
    }
}


// ═══════════════════════════════════════════════════════════════════
// Born radii (born.f90: compute_psi + compute_bornr)
// ═══════════════════════════════════════════════════════════════════

void ALPBSolvation::computePsi()
{
    // Psi accumulation and derivatives (Fortran born.f90:115-384)
    // m_brad is used to store psi temporarily, then transformed to Born radii

    m_brad.setZero();
    m_brdr.setZero();

    // Diagonal derivative accumulator (dpsitr in Fortran)
    Eigen::MatrixXd dpsitr = Eigen::MatrixXd::Zero(3, m_nat);

    for (int kk = 0; kk < m_nnrad; ++kk) {
        int ii = m_nnlistr[kk][0];
        int jj = m_nnlistr[kk][1];
        int nn = m_nnlistr[kk][2];

        double r = m_ddpair(0, nn);
        double dx = m_ddpair(1, nn);
        double dy = m_ddpair(2, nn);
        double dz = m_ddpair(3, nn);

        double rhoi = m_rho(ii);
        double rhoj = m_rho(jj);
        double rvdwi = m_vdwr(ii);
        double rvdwj = m_vdwr(jj);

        int ovij = (r < rvdwi + rhoj) ? 1 : 0;
        int ovji = (r < rhoi + rvdwj) ? 1 : 0;
        int ov = ovij + 10 * ovji;

        auto addPsiNoOverlap = [&](double rho_use, double rvdw_use, int atom_self, int atom_other, int sign) {
            // Non-overlapping contribution
            double r1 = 1.0 / r;
            double ap = r + rho_use;
            double am = r - rho_use;
            double ab = ap * am;
            double rhab = rho_use / ab;
            double lnab = 0.5 * std::log(am / ap) * r1;
            double gi = rhab + lnab;
            double dgi = -2.0 * rhab / ab + (rhab - lnab) * r1 * r1;

            m_brad(atom_self) += gi;

            double drjj_x = dgi * dx * sign;
            double drjj_y = dgi * dy * sign;
            double drjj_z = dgi * dz * sign;

            dpsitr(0, atom_self) += drjj_x;
            dpsitr(1, atom_self) += drjj_y;
            dpsitr(2, atom_self) += drjj_z;

            // Off-diagonal: brdr(3*nat × nat), stored as brdr(3*j+{0,1,2}, i) for ∂psi_i/∂xyz_j
            m_brdr(3*atom_other+0, atom_self) -= drjj_x;
            m_brdr(3*atom_other+1, atom_self) -= drjj_y;
            m_brdr(3*atom_other+2, atom_self) -= drjj_z;
        };

        auto addPsiOverlap = [&](double rho_use, double rvdw_check, int atom_self, int atom_other, int sign) {
            // Overlapping contribution (born.f90:248-273 or 277-302)
            if (r + rho_use <= rvdw_check) return;  // Fully inside

            double r1 = 1.0 / r;
            double r12 = 0.5 * r1;
            double r24 = r12 * r12;
            double ap = r + rho_use;
            double am = r - rho_use;
            double rh1 = 1.0 / rvdw_check;
            double rhr1 = 1.0 / ap;
            double aprh1 = ap * rh1;
            double lnab = std::log(aprh1);

            double gi = rh1 - rhr1 + r12 * (0.5 * am * (rhr1 - rh1 * aprh1) - lnab);

            double dgi = rhr1*rhr1 * (1.0 - 0.25*am*r1*(1.0 + aprh1*aprh1)) +
                         rho_use * r24 * (rhr1 - rh1*aprh1) +
                         r12 * (r1*lnab - rhr1);
            dgi *= r1;

            m_brad(atom_self) += gi;

            double drjj_x = dgi * dx * sign;
            double drjj_y = dgi * dy * sign;
            double drjj_z = dgi * dz * sign;

            dpsitr(0, atom_self) += drjj_x;
            dpsitr(1, atom_self) += drjj_y;
            dpsitr(2, atom_self) += drjj_z;

            m_brdr(3*atom_other+0, atom_self) -= drjj_x;
            m_brdr(3*atom_other+1, atom_self) -= drjj_y;
            m_brdr(3*atom_other+2, atom_self) -= drjj_z;
        };

        switch (ov) {
        case 0:  // No overlap
            if (std::abs(rhoi - rhoj) < 1e-8) {
                // Equal radii: symmetric contribution
                double r1 = 1.0 / r;
                double ap = r + rhoj;
                double am = r - rhoj;
                double ab = ap * am;
                double rhab = rhoj / ab;
                double lnab = 0.5 * std::log(am / ap) * r1;
                double gi = rhab + lnab;
                double dgi = -2.0 * rhab / ab + (rhab - lnab) * r1 * r1;

                m_brad(ii) += gi;
                m_brad(jj) += gi;

                double drx = dgi * dx;
                double dry = dgi * dy;
                double drz = dgi * dz;

                dpsitr(0,ii) += drx; dpsitr(1,ii) += dry; dpsitr(2,ii) += drz;
                m_brdr(3*jj+0,ii) -= drx; m_brdr(3*jj+1,ii) -= dry; m_brdr(3*jj+2,ii) -= drz;
                dpsitr(0,jj) -= drx; dpsitr(1,jj) -= dry; dpsitr(2,jj) -= drz;
                m_brdr(3*ii+0,jj) += drx; m_brdr(3*ii+1,jj) += dry; m_brdr(3*ii+2,jj) += drz;
            } else {
                // Unequal radii
                addPsiNoOverlap(rhoj, rvdwi, ii, jj, +1);
                addPsiNoOverlap(rhoi, rvdwj, jj, ii, -1);
            }
            break;

        case 10:  // ij no overlap, ji overlap
            addPsiNoOverlap(rhoj, rvdwi, ii, jj, +1);
            addPsiOverlap(rhoi, rvdwj, jj, ii, -1);
            break;

        case 1:   // ij overlap, ji no overlap
            addPsiOverlap(rhoj, rvdwi, ii, jj, +1);
            addPsiNoOverlap(rhoi, rvdwj, jj, ii, -1);
            break;

        case 11:  // Both overlap
            addPsiOverlap(rhoj, rvdwi, ii, jj, +1);
            addPsiOverlap(rhoi, rvdwj, jj, ii, -1);
            break;
        }
    }

    // Save diagonal terms
    for (int i = 0; i < m_nat; ++i) {
        m_brdr(3*i+0, i) = dpsitr(0, i);
        m_brdr(3*i+1, i) = dpsitr(1, i);
        m_brdr(3*i+2, i) = dpsitr(2, i);
    }
}

void ALPBSolvation::computeBornRadii()
{
    // Step 1: Compute psi values and derivatives
    computePsi();

    // Step 2: Transform psi to Born radii (Fortran born.f90:79-112)
    for (int i = 0; i < m_nat; ++i) {
        double br = m_brad(i);
        double svdwi = m_svdw(i);
        double vdwri = m_vdwr(i);
        double s1 = 1.0 / svdwi;
        double v1 = 1.0 / vdwri;
        double s2 = 0.5 * svdwi;

        br = br * s2;

        double arg2 = br * (obc_gamma * br - obc_beta);
        double arg = br * (obc_alpha + arg2);
        arg2 = 2.0 * arg2 + obc_alpha + obc_gamma * br * br;

        double th = std::tanh(arg);
        double ch = std::cosh(arg);

        double born_radius = m_born_scale / (s1 - v1 * th);

        double dpsi = ch * (s1 - v1 * th);
        dpsi = s2 * v1 * arg2 / (dpsi * dpsi);
        dpsi = m_born_scale * dpsi;

        m_brad(i) = born_radius;

        // Scale brdr derivatives by dpsi
        for (int j = 0; j < m_nat; ++j) {
            m_brdr(3*j+0, i) *= dpsi;
            m_brdr(3*j+1, i) *= dpsi;
            m_brdr(3*j+2, i) *= dpsi;
        }
    }
}


// ═══════════════════════════════════════════════════════════════════
// SASA (sasa.f90: compute_numsa + compute_w_sp)
// ═══════════════════════════════════════════════════════════════════

void ALPBSolvation::computeWSP(int nat, const int* nnlists_i, int nno,
                                const Eigen::Vector3d& xyzp,
                                const Eigen::MatrixXd& xyza,
                                double& sasap, Eigen::MatrixXd& grds,
                                int& nni, std::vector<int>& grdi) const
{
    // Fortran sasa.f90:145-199
    nni = 0;
    sasap = 1.0;

    for (int i = 0; i < nno; ++i) {
        int ia = nnlists_i[i];

        double tj0 = xyzp(0) - xyza(ia, 0);
        double tj1 = xyzp(1) - xyza(ia, 1);
        double tj2 = xyzp(2) - xyza(ia, 2);
        double tj2_sq = tj0*tj0 + tj1*tj1 + tj2*tj2;

        if (tj2_sq < m_trj2(1, ia)) {
            if (tj2_sq <= m_trj2(0, ia)) {
                sasap = 0.0;
                return;
            }

            double sqtj = std::sqrt(tj2_sq);
            double uj = sqtj - m_vdwsa(ia);
            double ah3uj2 = ah3 * uj * uj;
            double dsasaij = ah1 + 3.0 * ah3uj2;
            double sasaij = ah0 + (ah1 + ah3uj2) * uj;

            sasap *= sasaij;
            dsasaij = dsasaij / (sasaij * sqtj);

            grdi.push_back(ia);
            grds(0, nni) = dsasaij * tj0;
            grds(1, nni) = dsasaij * tj1;
            grds(2, nni) = dsasaij * tj2;
            nni++;
        }
    }
}

void ALPBSolvation::computeSASA(const Matrix& xyz_bohr)
{
    // Fortran sasa.f90:47-142
    m_sasa.setZero();
    for (int i = 0; i < m_nat; ++i) {
        m_dsdrt[i].setZero();
    }

    // Need xyz in column-major for SASA: xyz_col(atom, coord)
    // Our xyz_bohr is row-major N×3

    int max_nn = *std::max_element(m_nnsas.begin(), m_nnsas.end());
    if (max_nn == 0) max_nn = 1;

    for (int iat = 0; iat < m_nat; ++iat) {
        double rsas = m_vdwsa(iat);
        int nno = m_nnsas[iat];

        Eigen::MatrixXd grads = Eigen::MatrixXd::Zero(3, m_nat);
        double sasai = 0.0;

        Eigen::Vector3d xyza_i(xyz_bohr(iat, 0), xyz_bohr(iat, 1), xyz_bohr(iat, 2));
        double wr = m_wrp(iat);

        Eigen::MatrixXd grds = Eigen::MatrixXd::Zero(3, std::max(nno, 1));
        std::vector<int> grdi;
        grdi.reserve(nno);

        for (int ip = 0; ip < m_nang; ++ip) {
            Eigen::Vector3d xyzp = xyza_i + rsas * m_ang_grid.col(ip);

            double sasap;
            int nni;
            grdi.clear();

            // Prepare neighbor list as array
            const int* nn_data = m_nnlists[iat].empty() ? nullptr : m_nnlists[iat].data();

            computeWSP(m_nat, nn_data, nno, xyzp, xyz_bohr, sasap, grds, nni, grdi);

            if (sasap > tolsesp) {
                double wsa = m_ang_weight(ip) * wr * sasap;
                sasai += wsa;

                for (int jj = 0; jj < nni; ++jj) {
                    int nnj = grdi[jj];
                    double drx = wsa * grds(0, jj);
                    double dry = wsa * grds(1, jj);
                    double drz = wsa * grds(2, jj);
                    grads(0, iat) += drx;
                    grads(1, iat) += dry;
                    grads(2, iat) += drz;
                    grads(0, nnj) -= drx;
                    grads(1, nnj) -= dry;
                    grads(2, nnj) -= drz;
                }
            }
        }

        m_sasa(iat) = sasai;
        m_dsdrt[iat] = grads;
    }
}


// ═══════════════════════════════════════════════════════════════════
// Born interaction matrix (P16 kernel) — kernel.f90:713-767
// ═══════════════════════════════════════════════════════════════════

void ALPBSolvation::buildBornMatrix()
{
    // P16 kernel (kernel.f90: addBornMatP16)
    for (int kk = 0; kk < m_ntpair; ++kk) {
        double r1 = m_ddpair(0, kk);
        int iat = m_ppind[kk].first;
        int jat = m_ppind[kk].second;

        double ab = std::sqrt(m_brad(iat) * m_brad(jat));
        double arg = ab / (ab + zetaP16o16 * r1);
        arg = arg * arg;  // ^2
        arg = arg * arg;  // ^4
        arg = arg * arg;  // ^8
        arg = arg * arg;  // ^16
        double fgb = r1 + ab * arg;
        double dfgb = 1.0 / fgb;

        m_born_mat(iat, jat) += m_keps * dfgb;
        m_born_mat(jat, iat) += m_keps * dfgb;
    }

    // Self-energy
    for (int i = 0; i < m_nat; ++i) {
        m_born_mat(i, i) += m_keps / m_brad(i);
    }

    // HB correction on diagonal (Fortran gbsa.f90:508-510)
    if (m_lhb) {
        for (int i = 0; i < m_nat; ++i) {
            m_born_mat(i, i) += 2.0 * m_hbw(i);
        }
    }

    // ALPB shape correction (Fortran gbsa.f90:514-517)
    if (m_alpbet > 0.0) {
        double alpb_corr = m_keps * m_alpbet / m_adet;
        m_born_mat.array() += alpb_corr;
    }
}


// ═══════════════════════════════════════════════════════════════════
// HB correction (Fortran gbsa.f90:1003-1042)
// ═══════════════════════════════════════════════════════════════════

void ALPBSolvation::computeHBCorrection()
{
    for (int i = 0; i < m_nat; ++i) {
        double smaxd = 1.0 / (m_vdwsa(i) * m_vdwsa(i));
        double sasad = m_sasa(i) * smaxd;
        m_hbw(i) = m_hbmag(i) * sasad;
        m_dhbdw(i) = m_hbmag(i) * smaxd;
    }
}

void ALPBSolvation::addGradientHBond(const Vector& charges, double& ghb, Matrix& gradient)
{
    // Fortran gbsa.f90:1046-1091
    ghb = 0.0;
    for (int i = 0; i < m_nat; ++i) {
        ghb += m_hbw(i) * charges(i) * charges(i);
    }

    for (int i = 0; i < m_nat; ++i) {
        double dhbed = m_dhbdw(i);
        if (std::abs(dhbed) <= 0.0) continue;
        dhbed *= charges(i) * charges(i);
        for (int j = 0; j < m_nat; ++j) {
            gradient(j, 0) += m_dsdrt[i](0, j) * dhbed;
            gradient(j, 1) += m_dsdrt[i](1, j) * dhbed;
            gradient(j, 2) += m_dsdrt[i](2, j) * dhbed;
        }
    }
}


// ═══════════════════════════════════════════════════════════════════
// ALPB shape descriptor (Fortran gbsa.f90:677-719)
// ═══════════════════════════════════════════════════════════════════

void ALPBSolvation::computeADet(const Matrix& xyz_bohr)
{
    constexpr double tof = 2.0 / 5.0;

    double totRad3 = 0.0;
    Eigen::Vector3d center = Eigen::Vector3d::Zero();

    for (int i = 0; i < m_nat; ++i) {
        double rad2 = m_vdwr(i) * m_vdwr(i);
        double rad3 = rad2 * m_vdwr(i);
        totRad3 += rad3;
        center(0) += xyz_bohr(i, 0) * rad3;
        center(1) += xyz_bohr(i, 1) * rad3;
        center(2) += xyz_bohr(i, 2) * rad3;
    }
    center /= totRad3;

    Eigen::Matrix3d inertia = Eigen::Matrix3d::Zero();
    for (int i = 0; i < m_nat; ++i) {
        double rad2 = m_vdwr(i) * m_vdwr(i);
        double rad3 = rad2 * m_vdwr(i);
        Eigen::Vector3d vec(xyz_bohr(i,0) - center(0),
                           xyz_bohr(i,1) - center(1),
                           xyz_bohr(i,2) - center(2));
        double r2 = vec.squaredNorm();
        inertia += rad3 * ((r2 + tof * rad2) * Eigen::Matrix3d::Identity() - vec * vec.transpose());
    }

    double det = inertia.determinant();
    m_adet = std::sqrt(std::pow(det, 1.0/3.0) / (tof * totRad3));
}

void ALPBSolvation::addADetDeriv(const Matrix& xyz_bohr, double kEps_alpbet,
                                  const Vector& charges, Matrix& gradient)
{
    // Fortran gbsa.f90:722-789
    constexpr double tof = 2.0 / 5.0;

    double qtotal = charges.sum();
    double totRad3 = 0.0;
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (int i = 0; i < m_nat; ++i) {
        double rad2 = m_vdwr(i) * m_vdwr(i);
        double rad3 = rad2 * m_vdwr(i);
        totRad3 += rad3;
        center(0) += xyz_bohr(i, 0) * rad3;
        center(1) += xyz_bohr(i, 1) * rad3;
        center(2) += xyz_bohr(i, 2) * rad3;
    }
    center /= totRad3;

    Eigen::Matrix3d inertia = Eigen::Matrix3d::Zero();
    for (int i = 0; i < m_nat; ++i) {
        double rad2 = m_vdwr(i) * m_vdwr(i);
        double rad3 = rad2 * m_vdwr(i);
        Eigen::Vector3d vec(xyz_bohr(i,0)-center(0), xyz_bohr(i,1)-center(1), xyz_bohr(i,2)-center(2));
        double r2 = vec.squaredNorm();
        inertia += rad3 * ((r2 + tof*rad2) * Eigen::Matrix3d::Identity() - vec * vec.transpose());
    }

    double det = inertia.determinant();
    double aDet_local = std::sqrt(std::pow(det, 1.0/3.0) / (tof * totRad3));

    // Cofactor matrix (adjugate / transpose of cofactors)
    // Fortran gbsa.f90:769-779
    Eigen::Matrix3d aDeriv;
    aDeriv(0,0) = inertia(0,0)*(inertia(1,1)+inertia(2,2)) - inertia(0,1)*inertia(0,1) - inertia(0,2)*inertia(0,2);
    aDeriv(0,1) = inertia(0,1)*inertia(2,2) - inertia(0,2)*inertia(1,2);
    aDeriv(0,2) = inertia(0,2)*inertia(1,1) - inertia(0,1)*inertia(2,1);
    aDeriv(1,0) = aDeriv(0,1);
    aDeriv(1,1) = inertia(1,1)*(inertia(0,0)+inertia(2,2)) - inertia(0,1)*inertia(0,1) - inertia(1,2)*inertia(1,2);
    aDeriv(1,2) = inertia(0,0)*inertia(1,2) - inertia(0,1)*inertia(0,2);
    aDeriv(2,0) = aDeriv(0,2);
    aDeriv(2,1) = aDeriv(1,2);
    aDeriv(2,2) = inertia(2,2)*(inertia(0,0)+inertia(1,1)) - inertia(0,2)*inertia(0,2) - inertia(1,2)*inertia(1,2);

    // Scale factor (Fortran gbsa.f90:779-780)
    double scale = (250.0 / (48.0 * totRad3*totRad3*totRad3 * std::pow(aDet_local, 5.0)))
                 * (-0.5 * kEps_alpbet * qtotal * qtotal / (aDet_local * aDet_local));
    aDeriv *= scale;

    for (int i = 0; i < m_nat; ++i) {
        double rad2 = m_vdwr(i) * m_vdwr(i);
        double rad3 = rad2 * m_vdwr(i);
        Eigen::Vector3d vec(xyz_bohr(i,0)-center(0), xyz_bohr(i,1)-center(1), xyz_bohr(i,2)-center(2));
        Eigen::Vector3d dE = rad3 * (aDeriv * vec);
        gradient(i, 0) += dE(0);
        gradient(i, 1) += dE(1);
        gradient(i, 2) += dE(2);
    }
}


// ═══════════════════════════════════════════════════════════════════
// Energy computation (Fortran gbsa.f90:549-615)
// ═══════════════════════════════════════════════════════════════════

double ALPBSolvation::getEnergy(const Vector& charges) const
{
    ALPBEnergyParts parts = getEnergyParts(charges);
    return parts.total();
}

ALPBEnergyParts ALPBSolvation::getEnergyParts(const Vector& charges) const
{
    ALPBEnergyParts parts;

    // Born energy: E = 0.5 * q^T * B * q
    // (including HB and ALPB correction which are already in bornMat)
    Vector shift = 0.5 * m_born_mat * charges;
    double total_born = charges.dot(shift);

    // HB energy (separate for decomposition)
    if (m_lhb) {
        for (int i = 0; i < m_nat; ++i) {
            parts.ghb += m_hbw(i) * charges(i) * charges(i);
        }
    }

    // gborn = total_born - ghb (HB is counted in bornMat diagonal)
    parts.gborn = total_born - parts.ghb;
    parts.gsasa = m_gsasa_total;
    parts.gshift = m_gshift;

    return parts;
}


// ═══════════════════════════════════════════════════════════════════
// Gradient (Fortran gbsa.f90:618-674)
// ═══════════════════════════════════════════════════════════════════

void ALPBSolvation::addGradientP16(const Vector& charges, double& energy, Matrix& gradient)
{
    // P16 gradient (kernel.f90: addGradientP16, lines 283-383)
    Eigen::VectorXd dEdbr = Eigen::VectorXd::Zero(m_nat);
    double egb = 0.0;

    // Pair contributions
    for (int kk = 0; kk < m_ntpair; ++kk) {
        double r1 = m_ddpair(0, kk);
        double r2 = r1 * r1;
        int iat = m_ppind[kk].first;
        int jat = m_ppind[kk].second;
        double qq = charges(iat) * charges(jat);

        double ab = std::sqrt(m_brad(iat) * m_brad(jat));
        double arg1 = ab / (ab + zetaP16o16 * r1);
        double arg16 = arg1 * arg1;
        arg16 = arg16 * arg16;
        arg16 = arg16 * arg16;
        arg16 = arg16 * arg16;

        double fgb = r1 + ab * arg16;
        double dfgb = 1.0 / fgb;
        double dfgb2 = dfgb * dfgb;

        egb += qq * m_keps * dfgb;

        // Direct gradient
        double ap = (1.0 - zetaP16 * arg1 * arg16) * dfgb2;
        double dGx = ap * m_ddpair(1, kk) * m_keps / r1 * qq;
        double dGy = ap * m_ddpair(2, kk) * m_keps / r1 * qq;
        double dGz = ap * m_ddpair(3, kk) * m_keps / r1 * qq;
        gradient(iat, 0) -= dGx; gradient(iat, 1) -= dGy; gradient(iat, 2) -= dGz;
        gradient(jat, 0) += dGx; gradient(jat, 1) += dGy; gradient(jat, 2) += dGz;

        // Born radii derivatives
        double bp = -0.5 * (r1 * zetaP16 / ab * arg1 + 1.0) / ab * arg16 * dfgb2;
        dEdbr(iat) += m_brad(jat) * bp * m_keps * qq;
        dEdbr(jat) += m_brad(iat) * bp * m_keps * qq;
    }

    // Self-energy
    for (int i = 0; i < m_nat; ++i) {
        double bp = 1.0 / m_brad(i);
        double qq = charges(i) * bp;
        egb += 0.5 * charges(i) * qq * m_keps;
        double dEdbri = -0.5 * m_keps * qq * bp;
        dEdbr(i) += dEdbri * charges(i);
    }

    // Contract Born radii derivatives: gradient += brdr * dEdbr
    // brdr is (3*nat × nat), dEdbr is (nat)
    // For each atom j: gradient(j,:) += sum_i brdr(3*j+{0,1,2}, i) * dEdbr(i)
    for (int i = 0; i < m_nat; ++i) {
        if (std::abs(dEdbr(i)) < 1e-30) continue;
        for (int j = 0; j < m_nat; ++j) {
            gradient(j, 0) += m_brdr(3*j+0, i) * dEdbr(i);
            gradient(j, 1) += m_brdr(3*j+1, i) * dEdbr(i);
            gradient(j, 2) += m_brdr(3*j+2, i) * dEdbr(i);
        }
    }

    energy = egb;
}

void ALPBSolvation::addGradient(const std::vector<int>& atomic_numbers,
                                 const Matrix& xyz_bohr,
                                 const Vector& charges,
                                 Matrix& gradient)
{
    double gborn = 0.0;

    // P16 Born gradient
    addGradientP16(charges, gborn, gradient);

    // SASA gradient: dsdr already contracted
    for (int j = 0; j < m_nat; ++j) {
        gradient(j, 0) += m_dsdr(0, j);
        gradient(j, 1) += m_dsdr(1, j);
        gradient(j, 2) += m_dsdr(2, j);
    }

    // HB gradient
    double ghb = 0.0;
    if (m_lhb) {
        addGradientHBond(charges, ghb, gradient);
    }

    // ALPB shape gradient
    if (m_alpbet > 0.0) {
        gborn += charges.sum() * charges.sum() * m_alpbet / m_adet * m_keps;
        addADetDeriv(xyz_bohr, m_keps * m_alpbet, charges, gradient);
    }
}


// ═══════════════════════════════════════════════════════════════════
// Info output
// ═══════════════════════════════════════════════════════════════════

void ALPBSolvation::printInfo() const
{
    CurcumaLogger::info("ALPB Solvation Model");
    CurcumaLogger::param("Solvent", m_solvent);
    CurcumaLogger::param("Dielectric constant", std::to_string(m_dielectric_const));
    CurcumaLogger::param("Born scaling (c1)", std::to_string(m_born_scale));
    CurcumaLogger::param("keps", std::to_string(m_keps));
    CurcumaLogger::param("alpbet", std::to_string(m_alpbet));
    CurcumaLogger::param("Free energy shift", std::to_string(m_gshift) + " Eh");
    CurcumaLogger::param("Grid points", std::to_string(m_nang) + " per atom");
    CurcumaLogger::param("HB correction", m_lhb ? "true" : "false");
    CurcumaLogger::param("Kernel", "P16");
}
