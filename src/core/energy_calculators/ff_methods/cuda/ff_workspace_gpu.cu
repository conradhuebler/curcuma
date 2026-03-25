/*
 * <FFWorkspaceGPU Implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): GPU implementation of FFWorkspaceGPU.
 *
 * This file defines:
 *   - FFWorkspaceGPUImpl   (CUDA-internal: all SoA buffers + stream)
 *   - SoA upload methods   (DispersionSoA::upload, etc.)
 *   - FFWorkspaceGPU class (constructor, calculate, postProcessCPU, getters)
 *
 * Kernel launch layout: blockDim = 256, gridDim = ceil(n/256).
 * All gradient accumulation via atomicAdd on double (requires compute ≥ 6.0).
 *
 * Reference: Spicher/Grimme J. Chem. Theory Comput. 2020 (GFN-FF)
 */

#include "ff_workspace_gpu.h"
#include "gfnff_soa.h"
#include "gfnff_kernels.cuh"

#include "../gfnff_parameters.h"
#include "../forcefieldthread.h"   // Bond, Angle, Dihedral, Inversion structs
#include "src/core/energy_calculators/ff_methods/gfnff_par.h"  // covalent_rad_d3

#include "src/core/curcuma_logger.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstring>
#include <stdexcept>
#include <vector>

// Matrix, Vector, SpMatrix are defined in global.h (via ff_workspace_gpu.h -> ff_workspace.h)
// global.h: Matrix = Eigen::Matrix<double, Dynamic, Dynamic, RowMajor>
// global.h: Vector = Eigen::VectorXd
// global.h: SpMatrix = Eigen::SparseMatrix<double>

// ============================================================================
// D3 Covalent radii for angle/dihedral distance damping (Bohr, WITHOUT 4/3 factor)
// Copied from gfnff_par.h since nvcc has trouble linking static const std::vector
// Reference: Pyykkö & Atsumi, Chem. Eur. J. 15, 2009, 188-197
// Values × aatoau where aatoau = 1/0.52917726 (Bohr conversion)
// ============================================================================
static const double s_rcov_d3_bohr[87] = {
    0.32/0.52917726, 0.46/0.52917726,                                         // H, He
    1.20/0.52917726, 0.94/0.52917726, 0.77/0.52917726, 0.75/0.52917726,       // Li-C
    0.71/0.52917726, 0.63/0.52917726, 0.64/0.52917726, 0.67/0.52917726,       // N-Ne
    1.40/0.52917726, 1.25/0.52917726, 1.13/0.52917726, 1.04/0.52917726,       // Na-Si
    1.10/0.52917726, 1.02/0.52917726, 0.99/0.52917726, 0.96/0.52917726,       // P-Ar
    1.76/0.52917726, 1.54/0.52917726,                                         // K, Ca
    1.33/0.52917726, 1.22/0.52917726, 1.21/0.52917726, 1.10/0.52917726,       // Sc-Cr
    1.07/0.52917726, 1.04/0.52917726, 1.00/0.52917726, 0.99/0.52917726,       // Mn-Ni
    1.01/0.52917726, 1.09/0.52917726,                                         // Cu, Zn
    1.12/0.52917726, 1.09/0.52917726, 1.15/0.52917726, 1.10/0.52917726,       // Ga-Se
    1.14/0.52917726, 1.17/0.52917726,                                         // Br, Kr
    1.64/0.52917726, 1.46/0.52917726, 1.31/0.52917726, 1.26/0.52917726,       // Rb-Sr
    1.23/0.52917726, 1.22/0.52917726, 1.22/0.52917726, 1.20/0.52917726,       // Y-Mo
    1.19/0.52917726, 1.19/0.52917726, 1.18/0.52917726, 1.17/0.52917726,       // Tc-Pd
    1.18/0.52917726, 1.20/0.52917726, 1.21/0.52917726, 1.23/0.52917726,       // Ag-Cd
    1.28/0.52917726, 1.28/0.52917726, 1.28/0.52917726, 1.27/0.52917726,       // In-Te
    1.27/0.52917726, 1.34/0.52917726,                                         // I, Xe
    1.94/0.52917726, 1.71/0.52917726, 1.58/0.52917726, 1.51/0.52917726,       // Cs-Nd
    1.44/0.52917726, 1.44/0.52917726, 1.44/0.52917726, 1.43/0.52917726,       // Pm-Eu
    1.43/0.52917726, 1.43/0.52917726, 1.43/0.52917726, 1.40/0.52917726,       // Gd-Dy
    1.39/0.52917726, 1.39/0.52917726, 1.40/0.52917726, 1.41/0.52917726,       // Ho-Hg
    1.39/0.52917726, 1.39/0.52917726, 1.38/0.52917726, 1.38/0.52917726,       // Tl-Po
    1.38/0.52917726, 1.42/0.52917726,                                         // At, Rn
    2.01/0.52917726, 1.81/0.52917726,                                         // Fr, Ra
    1.67/0.52917726, 1.58/0.52917726                                          // Ac, Th
};

// ============================================================================
// GPU thread launch helper
// ============================================================================

static inline int gridFor(int n, int block = 256) {
    return (n + block - 1) / block;
}

// ============================================================================
// CUDA error check (host-side, after sync)
// ============================================================================

static void checkCuda(cudaError_t err, const char* msg) {
    if (err != cudaSuccess)
        throw std::runtime_error(std::string(msg) + ": " + cudaGetErrorString(err));
}

// ============================================================================
// SoA upload implementations (declared in gfnff_soa.h)
// ============================================================================

// ---------------------------------------------------------------------------
void DispersionSoA::upload(const std::vector<GFNFFDispersion>& v, cudaStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    std::vector<int>    h_i(n), h_j(n);
    std::vector<double> h_C6(n), h_r4r2(n), h_r0sq(n), h_zeta(n), h_rcut(n);

    for (int k = 0; k < n; ++k) {
        h_i[k]    = v[k].i;
        h_j[k]    = v[k].j;
        h_C6[k]   = v[k].C6;
        h_r4r2[k] = v[k].r4r2ij;
        h_r0sq[k] = v[k].r0_squared;
        h_zeta[k] = v[k].zetac6;
        h_rcut[k] = v[k].r_cut;
    }

    idx_i.upload(h_i, stream);
    idx_j.upload(h_j, stream);
    C6.upload(h_C6, stream);
    r4r2ij.upload(h_r4r2, stream);
    r0_sq.upload(h_r0sq, stream);
    zetac6.upload(h_zeta, stream);
    r_cut.upload(h_rcut, stream);
}

// ---------------------------------------------------------------------------
void RepulsionSoA::upload(const std::vector<GFNFFRepulsion>& v, cudaStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    std::vector<int>    h_i(n), h_j(n);
    std::vector<double> h_alpha(n), h_repab(n), h_rcut(n);

    for (int k = 0; k < n; ++k) {
        h_i[k]     = v[k].i;
        h_j[k]     = v[k].j;
        h_alpha[k] = v[k].alpha;
        h_repab[k] = v[k].repab;
        h_rcut[k]  = v[k].r_cut;
    }

    idx_i.upload(h_i, stream);
    idx_j.upload(h_j, stream);
    alpha.upload(h_alpha, stream);
    repab.upload(h_repab, stream);
    r_cut.upload(h_rcut, stream);
}

// ---------------------------------------------------------------------------
void CoulombSoA::upload(const std::vector<GFNFFCoulomb>& v, cudaStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    std::vector<int>    h_i(n), h_j(n);
    std::vector<double> h_gamma(n), h_rcut(n);

    for (int k = 0; k < n; ++k) {
        h_i[k]     = v[k].i;
        h_j[k]     = v[k].j;
        h_gamma[k] = v[k].gamma_ij;
        h_rcut[k]  = v[k].r_cut;
    }

    idx_i.upload(h_i, stream);
    idx_j.upload(h_j, stream);
    gamma_ij.upload(h_gamma, stream);
    r_cut.upload(h_rcut, stream);
}

// ---------------------------------------------------------------------------
void BondSoA::upload(const std::vector<Bond>& v, const std::vector<int>& atom_types,
                     cudaStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    std::vector<int>    h_i(n), h_j(n), h_nr_hb(n), h_hb_H(n);
    std::vector<double> h_r0(n), h_fc(n), h_alpha(n), h_rabshift(n), h_ff(n);
    std::vector<double> h_cnfak_i(n), h_cnfak_j(n), h_rb_i(n), h_rb_j(n);
    std::vector<double> h_hb_cn_H(n);

    for (int k = 0; k < n; ++k) {
        h_i[k]       = v[k].i;
        h_j[k]       = v[k].j;
        h_r0[k]      = v[k].r0_ij;
        h_fc[k]      = v[k].fc;
        h_alpha[k]   = v[k].exponent;
        h_rabshift[k]= v[k].rabshift;
        h_ff[k]      = v[k].ff;
        h_cnfak_i[k] = v[k].cnfak_i;
        h_cnfak_j[k] = v[k].cnfak_j;
        h_rb_i[k]    = v[k].r0_base_i;
        h_rb_j[k]    = v[k].r0_base_j;
        h_nr_hb[k]   = v[k].nr_hb;
        h_hb_cn_H[k] = v[k].hb_cn_H;
        // Determine which atom is H for HB alpha gradient
        if (v[k].nr_hb >= 1) {
            if (v[k].i < static_cast<int>(atom_types.size()) && atom_types[v[k].i] == 1)
                h_hb_H[k] = v[k].i;
            else if (v[k].j < static_cast<int>(atom_types.size()) && atom_types[v[k].j] == 1)
                h_hb_H[k] = v[k].j;
            else
                h_hb_H[k] = -1;
        } else {
            h_hb_H[k] = -1;
        }
    }

    idx_i.upload(h_i, stream);
    idx_j.upload(h_j, stream);
    r0.upload(h_r0, stream);
    fc.upload(h_fc, stream);
    alpha.upload(h_alpha, stream);
    rabshift.upload(h_rabshift, stream);
    ff.upload(h_ff, stream);
    cnfak_i.upload(h_cnfak_i, stream);
    cnfak_j.upload(h_cnfak_j, stream);
    r0_base_i.upload(h_rb_i, stream);
    r0_base_j.upload(h_rb_j, stream);
    nr_hb.upload(h_nr_hb, stream);
    hb_cn_H.upload(h_hb_cn_H, stream);
    hb_H_atom.upload(h_hb_H, stream);
}

// ---------------------------------------------------------------------------
void HBAlphaSoA::upload(const std::vector<int>& h_idx, const std::vector<int>& b_idx,
                        const std::vector<double>& rcov, cudaStream_t stream)
{
    n = static_cast<int>(h_idx.size());
    if (n == 0) return;
    idx_H.upload(h_idx, stream);
    idx_B.upload(b_idx, stream);
    rcov_sum.upload(rcov, stream);
}

// ---------------------------------------------------------------------------
void AngleSoA::upload(const std::vector<Angle>& v,
                      const std::vector<int>& atom_types,
                      cudaStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    const int nat = static_cast<int>(atom_types.size());
    std::vector<int>    h_i(n), h_j(n), h_k(n);
    std::vector<int>    h_ati(n), h_atj(n), h_atk(n);
    std::vector<double> h_fc(n), h_th(n);

    for (int m = 0; m < n; ++m) {
        h_i[m]   = v[m].i;
        h_j[m]   = v[m].j;
        h_k[m]   = v[m].k;
        h_ati[m] = (v[m].i < nat) ? atom_types[v[m].i] : 0;
        h_atj[m] = (v[m].j < nat) ? atom_types[v[m].j] : 0;
        h_atk[m] = (v[m].k < nat) ? atom_types[v[m].k] : 0;
        h_fc[m]  = v[m].fc;
        h_th[m]  = v[m].theta0_ijk;
    }

    idx_i.upload(h_i, stream); idx_j.upload(h_j, stream); idx_k.upload(h_k, stream);
    ati.upload(h_ati, stream); atj.upload(h_atj, stream); atk.upload(h_atk, stream);
    fc.upload(h_fc, stream);
    theta0.upload(h_th, stream);
}

// ---------------------------------------------------------------------------
void DihedralSoA::upload(const std::vector<Dihedral>& standard,
                         const std::vector<Dihedral>& extra,
                         const std::vector<int>& atom_types,
                         cudaStream_t stream)
{
    const int ns  = static_cast<int>(standard.size());
    const int ne  = static_cast<int>(extra.size());
    n = ns + ne;
    if (n == 0) return;

    const int nat = static_cast<int>(atom_types.size());
    std::vector<int>    h_i(n), h_j(n), h_k(n), h_l(n);
    std::vector<int>    h_ati(n), h_atj(n), h_atk(n), h_atl(n);
    std::vector<double> h_V(n), h_phi0(n), h_nper(n);
    std::vector<int>    h_isnci(n);

    auto fill = [&](int offset, const std::vector<Dihedral>& v) {
        for (int m = 0; m < static_cast<int>(v.size()); ++m) {
            const int idx = offset + m;
            h_i[idx]     = v[m].i;  h_j[idx] = v[m].j;
            h_k[idx]     = v[m].k;  h_l[idx] = v[m].l;
            h_ati[idx]   = (v[m].i < nat) ? atom_types[v[m].i] : 0;
            h_atj[idx]   = (v[m].j < nat) ? atom_types[v[m].j] : 0;
            h_atk[idx]   = (v[m].k < nat) ? atom_types[v[m].k] : 0;
            h_atl[idx]   = (v[m].l < nat) ? atom_types[v[m].l] : 0;
            h_V[idx]     = v[m].V;
            h_phi0[idx]  = v[m].phi0;
            h_nper[idx]  = v[m].n;
            h_isnci[idx] = v[m].is_nci ? 1 : 0;
        }
    };

    fill(0,  standard);
    fill(ns, extra);

    idx_i.upload(h_i, stream); idx_j.upload(h_j, stream);
    idx_k.upload(h_k, stream); idx_l.upload(h_l, stream);
    ati.upload(h_ati, stream); atj.upload(h_atj, stream);
    atk.upload(h_atk, stream); atl.upload(h_atl, stream);
    V.upload(h_V, stream);
    phi0.upload(h_phi0, stream);
    n_period.upload(h_nper, stream);
    is_nci.upload(h_isnci, stream);
}

// ---------------------------------------------------------------------------
void InversionSoA::upload(const std::vector<Inversion>& v,
                           const std::vector<int>& atom_types,
                           cudaStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    const int nat = static_cast<int>(atom_types.size());
    std::vector<int>    h_i(n), h_j(n), h_k(n), h_l(n), h_pt(n);
    std::vector<int>    h_ati(n), h_atj(n), h_atk(n), h_atl(n);
    std::vector<double> h_fc(n), h_om(n), h_C0(n), h_C1(n), h_C2(n);

    for (int m = 0; m < n; ++m) {
        h_i[m]  = v[m].i;  h_j[m] = v[m].j;
        h_k[m]  = v[m].k;  h_l[m] = v[m].l;
        h_ati[m] = (v[m].i < nat) ? atom_types[v[m].i] : 0;
        h_atj[m] = (v[m].j < nat) ? atom_types[v[m].j] : 0;
        h_atk[m] = (v[m].k < nat) ? atom_types[v[m].k] : 0;
        h_atl[m] = (v[m].l < nat) ? atom_types[v[m].l] : 0;
        h_fc[m] = v[m].fc;
        h_om[m] = v[m].omega0;
        h_C0[m] = v[m].C0;
        h_C1[m] = v[m].C1;
        h_C2[m] = v[m].C2;
        h_pt[m] = v[m].potential_type;
    }

    idx_i.upload(h_i, stream); idx_j.upload(h_j, stream);
    idx_k.upload(h_k, stream); idx_l.upload(h_l, stream);
    ati.upload(h_ati, stream); atj.upload(h_atj, stream);
    atk.upload(h_atk, stream); atl.upload(h_atl, stream);
    fc.upload(h_fc, stream);
    omega0.upload(h_om, stream);
    C0.upload(h_C0, stream);
    C1.upload(h_C1, stream);
    C2.upload(h_C2, stream);
    potential_type.upload(h_pt, stream);
}

// ---------------------------------------------------------------------------
void STorsionSoA::upload(const std::vector<GFNFFSTorsion>& v, cudaStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    std::vector<int>    h_i(n), h_j(n), h_k(n), h_l(n);
    std::vector<double> h_eref(n);

    for (int m = 0; m < n; ++m) {
        h_i[m]    = v[m].i;
        h_j[m]    = v[m].j;
        h_k[m]    = v[m].k;
        h_l[m]    = v[m].l;
        h_eref[m] = v[m].erefhalf;
    }

    idx_i.upload(h_i, stream);  idx_j.upload(h_j, stream);
    idx_k.upload(h_k, stream);  idx_l.upload(h_l, stream);
    erefhalf.upload(h_eref, stream);
}

// ---------------------------------------------------------------------------
void BATMSoA::upload(const std::vector<GFNFFBatmTriple>& v, cudaStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    std::vector<int>    h_i(n), h_j(n), h_k(n);
    std::vector<double> h_zi(n), h_zj(n), h_zk(n);

    for (int m = 0; m < n; ++m) {
        h_i[m]  = v[m].i;
        h_j[m]  = v[m].j;
        h_k[m]  = v[m].k;
        h_zi[m] = v[m].zb3atm_i;
        h_zj[m] = v[m].zb3atm_j;
        h_zk[m] = v[m].zb3atm_k;
    }

    idx_i.upload(h_i, stream);  idx_j.upload(h_j, stream);  idx_k.upload(h_k, stream);
    zb3atm_i.upload(h_zi, stream);
    zb3atm_j.upload(h_zj, stream);
    zb3atm_k.upload(h_zk, stream);
}

// ---------------------------------------------------------------------------
void ATMSoA::upload(const std::vector<ATMTriple>& v,
                     const std::vector<int>& atom_types,
                     cudaStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    const int nat = static_cast<int>(atom_types.size());
    std::vector<int>    h_i(n), h_j(n), h_k(n);
    std::vector<int>    h_ati(n), h_atj(n), h_atk(n);
    std::vector<double> h_c6ij(n), h_c6ik(n), h_c6jk(n);
    std::vector<double> h_s9(n), h_a1(n), h_a2(n), h_alp(n), h_tscale(n);

    for (int m = 0; m < n; ++m) {
        h_i[m]      = v[m].i;
        h_j[m]      = v[m].j;
        h_k[m]      = v[m].k;
        h_ati[m]    = (v[m].i < nat) ? atom_types[v[m].i] : 0;
        h_atj[m]    = (v[m].j < nat) ? atom_types[v[m].j] : 0;
        h_atk[m]    = (v[m].k < nat) ? atom_types[v[m].k] : 0;
        h_c6ij[m]   = v[m].C6_ij;
        h_c6ik[m]   = v[m].C6_ik;
        h_c6jk[m]   = v[m].C6_jk;
        h_s9[m]     = v[m].s9;
        h_a1[m]     = v[m].a1;
        h_a2[m]     = v[m].a2;
        h_alp[m]    = v[m].alp;
        h_tscale[m] = v[m].triple_scale;
    }

    idx_i.upload(h_i, stream);  idx_j.upload(h_j, stream);  idx_k.upload(h_k, stream);
    ati.upload(h_ati, stream);  atj.upload(h_atj, stream);  atk.upload(h_atk, stream);
    C6_ij.upload(h_c6ij, stream);
    C6_ik.upload(h_c6ik, stream);
    C6_jk.upload(h_c6jk, stream);
    s9.upload(h_s9, stream);
    a1.upload(h_a1, stream);
    a2.upload(h_a2, stream);
    alp.upload(h_alp, stream);
    triple_scale.upload(h_tscale, stream);
}

// ---------------------------------------------------------------------------
void XBondSoA::upload(const std::vector<GFNFFHalogenBond>& v,
                       const std::vector<int>& atom_types,
                       cudaStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    const int nat = static_cast<int>(atom_types.size());
    std::vector<int>    h_i(n), h_j(n), h_k(n);
    std::vector<int>    h_eA(n), h_eB(n);
    std::vector<double> h_qX(n), h_qB(n), h_aX(n), h_rcut(n);

    for (int m = 0; m < n; ++m) {
        h_i[m]    = v[m].i;   // donor A
        h_j[m]    = v[m].j;   // halogen X
        h_k[m]    = v[m].k;   // acceptor B
        h_eA[m]   = (v[m].i < nat) ? atom_types[v[m].i] : 0;
        h_eB[m]   = (v[m].k < nat) ? atom_types[v[m].k] : 0;
        h_qX[m]   = v[m].q_X;
        h_qB[m]   = v[m].q_B;
        h_aX[m]   = v[m].acidity_X;
        h_rcut[m] = v[m].r_cut;
    }

    idx_i.upload(h_i, stream);  idx_j.upload(h_j, stream);  idx_k.upload(h_k, stream);
    elem_A.upload(h_eA, stream);  elem_B.upload(h_eB, stream);
    q_X.upload(h_qX, stream);
    q_B.upload(h_qB, stream);
    acidity_X.upload(h_aX, stream);
    r_cut.upload(h_rcut, stream);
}

// ---------------------------------------------------------------------------
void HBondSoA::upload(const std::vector<GFNFFHydrogenBond>& v,
                       const std::vector<int>& atom_types,
                       cudaStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    const int nat = static_cast<int>(atom_types.size());
    std::vector<int>    h_i(n), h_j(n), h_k(n);
    std::vector<int>    h_eA(n), h_eB(n), h_case(n);
    std::vector<double> h_qH(n), h_qA(n), h_qB(n);
    std::vector<double> h_basA(n), h_basB(n), h_aciA(n), h_aciB(n);
    std::vector<double> h_rcut(n);
    std::vector<int>    h_nbB_offset(n), h_nbB_count(n);
    std::vector<int>    h_acc_parent(n);
    std::vector<int>    h_nbC_offset(n), h_nbC_count(n);
    std::vector<double> h_repz(n);

    // First pass: count total neighbors
    int total_B = 0, total_C = 0;
    for (int m = 0; m < n; ++m) {
        total_B += static_cast<int>(v[m].neighbors_B.size());
        total_C += static_cast<int>(v[m].neighbors_C.size());
    }

    std::vector<int> h_nbB_flat(total_B);
    std::vector<int> h_nbC_flat(total_C);
    int offB = 0, offC = 0;

    for (int m = 0; m < n; ++m) {
        h_i[m]     = v[m].i;   // donor A
        h_j[m]     = v[m].j;   // hydrogen H
        h_k[m]     = v[m].k;   // acceptor B
        h_eA[m]    = (v[m].i < nat) ? atom_types[v[m].i] : 0;
        h_eB[m]    = (v[m].k < nat) ? atom_types[v[m].k] : 0;
        h_qH[m]    = v[m].q_H;
        h_qA[m]    = v[m].q_A;
        h_qB[m]    = v[m].q_B;
        h_basA[m]  = v[m].basicity_A;
        h_basB[m]  = v[m].basicity_B;
        h_aciA[m]  = v[m].acidity_A;
        h_aciB[m]  = v[m].acidity_B;
        h_rcut[m]  = v[m].r_cut;
        h_case[m]  = v[m].case_type;
        h_acc_parent[m] = v[m].acceptor_parent_index;

        // repz_B for case 4
        int zB = (v[m].k < nat) ? atom_types[v[m].k] : 0;
        h_repz[m] = (zB >= 1 && zB <= 86) ? GFNFFParameters::repz[zB - 1] : 1.0;

        // Pack neighbors_B
        h_nbB_offset[m] = offB;
        h_nbB_count[m]  = static_cast<int>(v[m].neighbors_B.size());
        for (int nb : v[m].neighbors_B) h_nbB_flat[offB++] = nb;

        // Pack neighbors_C
        h_nbC_offset[m] = offC;
        h_nbC_count[m]  = static_cast<int>(v[m].neighbors_C.size());
        for (int nb : v[m].neighbors_C) h_nbC_flat[offC++] = nb;
    }

    total_nb_B = total_B;
    total_nb_C = total_C;

    idx_i.upload(h_i, stream);  idx_j.upload(h_j, stream);  idx_k.upload(h_k, stream);
    elem_A.upload(h_eA, stream);  elem_B.upload(h_eB, stream);
    q_H.upload(h_qH, stream);  q_A.upload(h_qA, stream);  q_B.upload(h_qB, stream);
    basicity_A.upload(h_basA, stream);  basicity_B.upload(h_basB, stream);
    acidity_A.upload(h_aciA, stream);   acidity_B.upload(h_aciB, stream);
    r_cut.upload(h_rcut, stream);
    case_type.upload(h_case, stream);
    nb_B_offset.upload(h_nbB_offset, stream);
    nb_B_count.upload(h_nbB_count, stream);
    if (total_B > 0) nb_B_flat.upload(h_nbB_flat, stream);
    acceptor_parent.upload(h_acc_parent, stream);
    nb_C_offset.upload(h_nbC_offset, stream);
    nb_C_count.upload(h_nbC_count, stream);
    if (total_C > 0) nb_C_flat.upload(h_nbC_flat, stream);
    repz_B.upload(h_repz, stream);
}

// ============================================================================
// FFWorkspaceGPUImpl — CUDA-internal state (not visible outside this .cu)
// ============================================================================

struct FFWorkspaceGPUImpl {
    // Static SoA buffers (uploaded once at construction)
    DispersionSoA  disp;            ///< D4 dispersion pairs
    RepulsionSoA   bonded_rep;      ///< Bonded repulsion pairs
    RepulsionSoA   nonbonded_rep;   ///< Non-bonded repulsion pairs
    CoulombSoA     coulomb;         ///< Coulomb pairs (gamma + cutoff only; charges dynamic)
    BondSoA        bonds;           ///< Bond stretching
    AngleSoA       angles;          ///< Angle bending
    DihedralSoA    dihedrals;       ///< Standard + extra torsions (combined)
    InversionSoA   inversions;      ///< Out-of-plane inversions

    // Phase 2 SoA buffers (March 2026): formerly CPU residual, now all GPU
    STorsionSoA    storsions;       ///< Triple bond torsions
    BATMSoA        batm;            ///< Bonded ATM 3-body
    ATMSoA         atm;             ///< Axilrod-Teller-Muto 3-body dispersion
    XBondSoA       xbonds;          ///< Halogen bonds
    HBondSoA       hbonds;          ///< Hydrogen bonds

    // Dynamic buffers (reset + upload each calculation step)
    CudaBuffer<double> d_coords;    ///< [N*3] row-major (x,y,z per atom)
    CudaBuffer<double> d_charges;   ///< [N] EEQ charges
    CudaBuffer<double> d_cn;        ///< [N] D3 coordination numbers
    CudaBuffer<double> d_topo_charges; ///< [N] topology charges (for BATM)
    CudaBuffer<double> d_grad;      ///< [N*3] gradient output (atomicAdd)
    CudaBuffer<double> d_dEdcn;     ///< [N] bond dE/dCN (atomicAdd)
    CudaBuffer<double> d_dlogdcn;  ///< [N] logistic squashing factor for CN chain-rule
    CudaBuffer<double> d_energy;    ///< [1] total GPU energy (atomicAdd)

    // CN chain-rule pair list (uploaded once, used every step)
    CudaBuffer<int>    d_cn_idx_i;     ///< [n_cn_pairs]
    CudaBuffer<int>    d_cn_idx_j;     ///< [n_cn_pairs]
    CudaBuffer<double> d_cn_rcov_sum;  ///< [n_cn_pairs]
    int                n_cn_pairs = 0;

    // HB alpha chain-rule: per-atom zz accumulator + HB neighbor pair list
    CudaBuffer<double> d_zz_hb;        ///< [N] per-atom zz accumulator (zeroed each step)
    HBAlphaSoA         hb_alpha;       ///< (H, B) pairs for alpha gradient

    // Coulomb self-energy parameters (uploaded once)
    CudaBuffer<double> d_coul_chi_base; ///< [N]
    CudaBuffer<double> d_coul_cnf;      ///< [N]
    CudaBuffer<double> d_coul_gam;      ///< [N]
    CudaBuffer<double> d_coul_alp;      ///< [N]
    bool               coul_self_on_gpu = false;

    // Per-term energy accumulators — consolidated into single contiguous buffer
    // Claude Generated (March 2026): Reduces 14 individual cudaMemcpy to 1
    enum EnergySlot {
        E_DISP = 0, E_BREP = 1, E_NBREP = 2, E_COUL = 3, E_BOND = 4,
        E_ANGLE = 5, E_DIHED = 6, E_INV = 7, E_STORS = 8, E_BATM = 9,
        E_ATM = 10, E_XBOND = 11, E_HBOND = 12, E_COUL_SELF = 13,
        N_ENERGY_SLOTS = 14
    };
    CudaBuffer<double> d_energies;      ///< [N_ENERGY_SLOTS] all per-term energies

    // Diagnostic snapshot buffers (GPU-side, no pipeline stall)
    // Claude Generated (March 2026): Device-to-device copy instead of sync+download
    CudaBuffer<double> d_dEdcn_snapshot;   ///< [N] copy of dEdcn before qtmp
    CudaBuffer<double> d_grad_snapshot;    ///< [3*N] copy of gradient before CN chain-rule

    int N = 0;
    cudaStream_t stream = nullptr;
};

// ============================================================================
// Helper: extract Eigen N×3 RowMajor matrix to flat [x0,y0,z0,x1,...] array
// Note: Matrix (global.h) is already RowMajor — element-by-element copy is safe
// ============================================================================

static std::vector<double> toRowMajor(const Matrix& m) {
    const int N = static_cast<int>(m.rows());
    std::vector<double> v(3 * N);
    for (int i = 0; i < N; ++i) {
        v[3*i+0] = m(i, 0);
        v[3*i+1] = m(i, 1);
        v[3*i+2] = m(i, 2);
    }
    return v;
}

// ============================================================================
// FFWorkspaceGPU Constructor
// ============================================================================

FFWorkspaceGPU::FFWorkspaceGPU(const GFNFFParameterSet& params,
                                int natoms,
                                const std::vector<int>& atom_types)
    : m_natoms(natoms)
{
    // Verify CUDA device
    int dev_count = 0;
    checkCuda(cudaGetDeviceCount(&dev_count), "cudaGetDeviceCount");
    if (dev_count == 0)
        throw std::runtime_error("FFWorkspaceGPU: no CUDA device found");

    // Pre-allocate CPU staging buffers BEFORE any CUDA operations.
    // CUDA corrupts C++ heap metadata; subsequent heap allocs crash.
    // These buffers are reused every step without new allocations.
    m_h_coords.resize(3 * natoms);
    m_h_grad.resize(3 * natoms);

    m_impl = std::make_unique<FFWorkspaceGPUImpl>();
    m_impl->N = natoms;

    // Create dedicated CUDA stream
    checkCuda(cudaStreamCreate(&m_impl->stream), "cudaStreamCreate");
    cudaStream_t stream = m_impl->stream;

    // Upload covalent radii to constant memory (used by angle/dihedral distance damping)
    upload_rcov_d3(s_rcov_d3_bohr, 87);

    // Upload covalent radii for HB/XB vdW radii lookup (Å units)
    upload_covalent_radii(GFNFFParameters::covalent_radii.data(),
                          static_cast<int>(GFNFFParameters::covalent_radii.size()));

    // Store dispersion pair indices on host for per-step dc6dcn extraction
    {
        const int nd = static_cast<int>(params.dispersions.size());
        m_disp_idx_i_host.resize(nd);
        m_disp_idx_j_host.resize(nd);
        m_h_dc6dcn_ij.resize(nd, 0.0);
        m_h_dc6dcn_ji.resize(nd, 0.0);
        for (int k = 0; k < nd; ++k) {
            m_disp_idx_i_host[k] = params.dispersions[k].i;
            m_disp_idx_j_host[k] = params.dispersions[k].j;
        }
    }

    // Store topology charges for BATM kernel (uploaded to GPU each calculate() call)
    m_topology_charges = params.topology_charges;

    // --- Upload static SoA interaction lists ---
    m_impl->disp.upload(params.dispersions, stream);
    m_impl->bonded_rep.upload(params.bonded_repulsions, stream);
    m_impl->nonbonded_rep.upload(params.nonbonded_repulsions, stream);
    m_impl->coulomb.upload(params.coulombs, stream);
    m_impl->bonds.upload(params.bonds, atom_types, stream);
    m_impl->angles.upload(params.angles, atom_types, stream);
    m_impl->dihedrals.upload(params.dihedrals, params.extra_dihedrals, atom_types, stream);
    m_impl->inversions.upload(params.inversions, atom_types, stream);

    // Phase 2 SoA uploads (March 2026): HB/XB/ATM/BATM/sTors — all on GPU
    m_impl->storsions.upload(params.storsions, stream);
    m_impl->batm.upload(params.batm_triples, stream);
    m_impl->atm.upload(params.atm_triples, atom_types, stream);
    m_impl->xbonds.upload(params.xbonds, atom_types, stream);
    m_impl->hbonds.upload(params.hbonds, atom_types, stream);

    // --- Allocate dynamic per-step buffers ---
    const int N3 = 3 * natoms;
    m_impl->d_coords.alloc(N3);
    m_impl->d_charges.alloc(natoms);
    m_impl->d_cn.alloc(natoms);
    m_impl->d_topo_charges.alloc(natoms);
    m_impl->d_grad.alloc(N3);
    m_impl->d_dEdcn.alloc(natoms);
    m_impl->d_energy.alloc(1);

    // Per-term energy accumulators — single contiguous buffer
    m_impl->d_energies.alloc(FFWorkspaceGPUImpl::N_ENERGY_SLOTS);

    // Diagnostic snapshot buffers (device-to-device copy, no sync stall)
    m_impl->d_dEdcn_snapshot.alloc(natoms);
    m_impl->d_grad_snapshot.alloc(N3);

    // --- HB alpha chain-rule: allocate zz buffer + build (H,B) pair list ---
    m_impl->d_zz_hb.alloc(natoms);
    if (!params.bond_hb_data.empty()) {
        const auto& rcov_d3 = GFNFFParameters::covalent_rad_d3;  // already in Bohr
        constexpr double rcov_43 = 4.0 / 3.0;
        constexpr double rcov_scal = 1.78;

        std::vector<int> hb_h_idx, hb_b_idx;
        std::vector<double> hb_rcov;
        for (const auto& entry : params.bond_hb_data) {
            int H = entry.H;
            int ati = (H < natoms) ? atom_types[H] : 0;
            for (int B : entry.B_atoms) {
                int atj = (B < natoms) ? atom_types[B] : 0;
                if (ati < 1 || atj < 1) continue;
                double rcovij = rcov_scal * rcov_43 * (rcov_d3[ati - 1] + rcov_d3[atj - 1]);  // already in Bohr
                hb_h_idx.push_back(H);
                hb_b_idx.push_back(B);
                hb_rcov.push_back(rcovij);
            }
        }
        m_impl->hb_alpha.upload(hb_h_idx, hb_b_idx, hb_rcov, stream);
    }

    // --- Pre-allocate CPU result buffers ---
    m_result_gradient.resize(natoms, 3);
    m_result_gradient.setZero();
    m_dEdcn_total.resize(natoms);
    m_dEdcn_total.setZero();

    // --- Extract per-atom Coulomb self-energy params (for CPU postProcess TERM 2+3) ---
    // Mirrors FFWorkspace::setInteractionLists() lines 90-118
    if (!params.coulombs.empty()) {
        m_coul_chi_base  = Vector::Zero(natoms);
        m_coul_gam       = Vector::Zero(natoms);
        m_coul_alp       = Vector::Zero(natoms);
        m_coul_cnf       = Vector::Zero(natoms);
        m_coul_chi_static= Vector::Zero(natoms);

        std::vector<bool> seen(natoms, false);
        for (const auto& c : params.coulombs) {
            if (c.i < natoms && !seen[c.i]) {
                m_coul_chi_base(c.i)   = c.chi_base_i;
                m_coul_gam(c.i)        = c.gam_i;
                m_coul_alp(c.i)        = c.alp_i;
                m_coul_cnf(c.i)        = c.cnf_i;
                m_coul_chi_static(c.i) = c.chi_i;
                seen[c.i] = true;
            }
            if (c.j < natoms && !seen[c.j]) {
                m_coul_chi_base(c.j)   = c.chi_base_j;
                m_coul_gam(c.j)        = c.gam_j;
                m_coul_alp(c.j)        = c.alp_j;
                m_coul_cnf(c.j)        = c.cnf_j;
                m_coul_chi_static(c.j) = c.chi_j;
                seen[c.j] = true;
            }
        }
        // Upload Coulomb self-energy parameters to GPU for k_coulomb_self kernel
        setCoulombSelfEnergyParams(m_coul_chi_base, m_coul_gam, m_coul_alp,
                                   m_coul_cnf, m_coul_chi_static);
    }

    // Store initial charges and E0 from parameter set
    if (params.eeq_charges.size() == natoms)
        m_eeq_charges = params.eeq_charges;
    m_e0 = params.e0;

    // Synchronise: wait for all uploads to finish
    checkCuda(cudaStreamSynchronize(stream), "init sync");

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format(
            "FFWorkspaceGPU: {} atoms | disp={} brep={} nbrep={} coul={} bond={} ang={} dih={} inv={}"
            " stors={} batm={} atm={} xb={} hb={}",
            natoms,
            m_impl->disp.n, m_impl->bonded_rep.n, m_impl->nonbonded_rep.n,
            m_impl->coulomb.n, m_impl->bonds.n, m_impl->angles.n,
            m_impl->dihedrals.n, m_impl->inversions.n,
            m_impl->storsions.n, m_impl->batm.n, m_impl->atm.n,
            m_impl->xbonds.n, m_impl->hbonds.n));
    }
}

// ============================================================================
// Destructor
// ============================================================================

FFWorkspaceGPU::~FFWorkspaceGPU()
{
    if (m_impl && m_impl->stream) {
        // Synchronize before destroying: pending async ops must finish
        // before CudaBuffer destructors free GPU memory
        cudaStreamSynchronize(m_impl->stream);
        cudaStreamDestroy(m_impl->stream);
    }
}

// ============================================================================
// Per-step state setters
// ============================================================================

void FFWorkspaceGPU::setEEQCharges(const Vector& q)
{
    m_eeq_charges = q;
}

void FFWorkspaceGPU::setTopologyCharges(const Vector& q)
{
    m_topology_charges = q;
}

void FFWorkspaceGPU::setD3CN(const Vector& cn)
{
    m_cn = cn;
}

void FFWorkspaceGPU::setCNDerivatives(const Vector& cn, const Vector& cnf,
                                       const std::vector<SpMatrix>& /* dcn */)
{
    // Only cnf is needed (for k_subtract_qtmp). CN pair list replaces sparse dcn.
    m_cnf = cnf;
}

void FFWorkspaceGPU::setE0(double e0)
{
    m_e0 = e0;
}

void FFWorkspaceGPU::updateBondHBCN(const std::vector<double>& hb_cn_values)
{
    const int nb = m_impl->bonds.n;
    if (nb == 0 || static_cast<int>(hb_cn_values.size()) != nb) return;
    m_impl->bonds.hb_cn_H.upload(hb_cn_values.data(), nb);
}

void FFWorkspaceGPU::setCoulombSelfEnergyParams(const Vector& chi_base, const Vector& gam,
                                                  const Vector& alp,     const Vector& cnf,
                                                  const Vector& chi_static)
{
    m_coul_chi_base   = chi_base;
    m_coul_gam        = gam;
    m_coul_alp        = alp;
    m_coul_cnf        = cnf;
    m_coul_chi_static = chi_static;

    // Upload Coulomb self-energy parameters to GPU
    // Claude Generated (March 2026): GPU-side Coulomb TERM 2+3 + qtmp
    const int N = m_natoms;
    if (chi_base.size() == N && gam.size() == N && alp.size() == N && cnf.size() == N) {
        auto& impl = *m_impl;
        impl.d_coul_chi_base.alloc(N);
        impl.d_coul_cnf.alloc(N);
        impl.d_coul_gam.alloc(N);
        impl.d_coul_alp.alloc(N);
        impl.d_coul_chi_base.upload(chi_base.data(), N);
        impl.d_coul_cnf.upload(cnf.data(), N);
        impl.d_coul_gam.upload(gam.data(), N);
        impl.d_coul_alp.upload(alp.data(), N);
        impl.coul_self_on_gpu = true;
    }
}

void FFWorkspaceGPU::updateHBonds(const std::vector<GFNFFHydrogenBond>& hbonds,
                                    const std::vector<int>& atom_types)
{
    // Re-upload HBond SoA after dynamic re-detection.
    // Claude Generated (March 2026): Needed because updateHBXBIfNeeded() may find
    // additional HBonds on re-detection that weren't in the initial parameter set.
    m_impl->hbonds.upload(hbonds, atom_types, m_impl->stream);
}

void FFWorkspaceGPU::updateXBonds(const std::vector<GFNFFHalogenBond>& xbonds,
                                    const std::vector<int>& atom_types)
{
    m_impl->xbonds.upload(xbonds, atom_types, m_impl->stream);
}

void FFWorkspaceGPU::setDispersionEnabled(bool v)  { m_dispersion_enabled = v; }
void FFWorkspaceGPU::setHBondEnabled(bool v)        { m_hbond_enabled = v; }
void FFWorkspaceGPU::setRepulsionEnabled(bool v)    { m_repulsion_enabled = v; }
void FFWorkspaceGPU::setCoulombEnabled(bool v)      { m_coulomb_enabled = v; }

// ============================================================================
// GPU-only CN chain-rule setters (Claude Generated March 2026)
// ============================================================================

void FFWorkspaceGPU::setCNPairList(const std::vector<int>& idx_i,
                                    const std::vector<int>& idx_j,
                                    const std::vector<double>& rcov_sum)
{
    const int np = static_cast<int>(idx_i.size());
    if (np == 0 || idx_j.size() != idx_i.size() || rcov_sum.size() != idx_i.size()) return;

    m_cn_pair_i    = idx_i;
    m_cn_pair_j    = idx_j;
    m_cn_pair_rcov = rcov_sum;
    m_cn_n_pairs   = np;

    // Upload to GPU
    auto& impl = *m_impl;
    impl.d_cn_idx_i.alloc(np);
    impl.d_cn_idx_j.alloc(np);
    impl.d_cn_rcov_sum.alloc(np);
    impl.d_cn_idx_i.upload(idx_i.data(), np);
    impl.d_cn_idx_j.upload(idx_j.data(), np);
    impl.d_cn_rcov_sum.upload(rcov_sum.data(), np);
    impl.n_cn_pairs = np;
    m_cn_pairs_on_gpu = true;

    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::info(fmt::format("  GPU CN chain-rule: {} pairs uploaded", np));
}

void FFWorkspaceGPU::setDlogDCN(const Vector& dlogdcn)
{
    m_dlogdcn = dlogdcn;
    // Upload to GPU
    auto& impl = *m_impl;
    const int N = m_natoms;
    if (dlogdcn.size() != N) return;
    if (!impl.d_dlogdcn.ptr) impl.d_dlogdcn.alloc(N);
    impl.d_dlogdcn.upload(dlogdcn.data(), N);
}

void FFWorkspaceGPU::updateDispersionDC6DCN(const Matrix& dc6dcn)
{
    // Claude Generated (March 2026): Extract per-pair dc6dcn values and upload to GPU.
    // Eliminates the CPU round-trip (download dEdcn, add dispersion, re-upload).
    const int nd = static_cast<int>(m_disp_idx_i_host.size());
    if (nd == 0 || dc6dcn.size() == 0) return;

    for (int k = 0; k < nd; ++k) {
        int i = m_disp_idx_i_host[k];
        int j = m_disp_idx_j_host[k];
        if (i < dc6dcn.rows() && j < dc6dcn.cols()) {
            m_h_dc6dcn_ij[k] = dc6dcn(i, j);
            m_h_dc6dcn_ji[k] = dc6dcn(j, i);
        } else {
            m_h_dc6dcn_ij[k] = 0.0;
            m_h_dc6dcn_ji[k] = 0.0;
        }
    }

    m_impl->disp.dc6dcn_ij.upload(m_h_dc6dcn_ij.data(), nd);
    m_impl->disp.dc6dcn_ji.upload(m_h_dc6dcn_ji.data(), nd);
}

// ============================================================================
// calculate() — main entry point
// ============================================================================

double FFWorkspaceGPU::calculate(bool gradient)
{
    auto& impl = *m_impl;
    const int N = m_natoms;
    cudaStream_t stream = impl.stream;
    constexpr int BLOCK = 256;

    // =========================================================================
    // 1. Zero GPU accumulators (total + per-term)
    // =========================================================================
    cudaMemsetAsync(impl.d_energy.ptr, 0, sizeof(double),         stream);
    cudaMemsetAsync(impl.d_grad.ptr,   0, 3*N*sizeof(double),     stream);
    cudaMemsetAsync(impl.d_dEdcn.ptr,  0, N*sizeof(double),       stream);
    // Zero all per-term energy accumulators in single call
    cudaMemsetAsync(impl.d_energies.ptr, 0,
                    FFWorkspaceGPUImpl::N_ENERGY_SLOTS * sizeof(double), stream);

    // =========================================================================
    // 2. Upload per-step charges and CN
    //    Geometry was already uploaded by setGeometry() before this call.
    //    d_coords is valid on the GPU at this point.
    // =========================================================================

    // Upload EEQ charges
    if (m_eeq_charges.size() == N) {
        impl.d_charges.upload(m_eeq_charges.data(), N, stream);
    }

    // Upload D3 CN
    if (m_cn.size() == N) {
        impl.d_cn.upload(m_cn.data(), N, stream);
    }

    // Upload topology charges (used by BATM kernel)
    if (m_topology_charges.size() == N) {
        impl.d_topo_charges.upload(m_topology_charges.data(), N, stream);
    }

    // =========================================================================
    // 3. Launch CUDA kernels (only for non-empty lists and enabled terms)
    // =========================================================================

    // --- Dispersion (D4 modified BJ formula + dEdcn chain-rule) ---
    if (m_dispersion_enabled && impl.disp.n > 0) {
        k_dispersion<<<gridFor(impl.disp.n, BLOCK), BLOCK, 0, stream>>>(
            impl.disp.n,
            impl.disp.idx_i.ptr,  impl.disp.idx_j.ptr,
            impl.disp.C6.ptr,     impl.disp.r4r2ij.ptr,
            impl.disp.r0_sq.ptr,  impl.disp.zetac6.ptr,
            impl.disp.r_cut.ptr,
            gradient ? impl.disp.dc6dcn_ij.ptr : nullptr,
            gradient ? impl.disp.dc6dcn_ji.ptr : nullptr,
            impl.d_coords.ptr,
            impl.d_dEdcn.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_DISP]);
    }

    // --- Bonded repulsion ---
    if (m_repulsion_enabled && impl.bonded_rep.n > 0) {
        k_repulsion<<<gridFor(impl.bonded_rep.n, BLOCK), BLOCK, 0, stream>>>(
            impl.bonded_rep.n,
            impl.bonded_rep.idx_i.ptr,  impl.bonded_rep.idx_j.ptr,
            impl.bonded_rep.alpha.ptr,  impl.bonded_rep.repab.ptr,
            impl.bonded_rep.r_cut.ptr,
            impl.d_coords.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_BREP]);
    }

    // --- Non-bonded repulsion ---
    if (m_repulsion_enabled && impl.nonbonded_rep.n > 0) {
        k_repulsion<<<gridFor(impl.nonbonded_rep.n, BLOCK), BLOCK, 0, stream>>>(
            impl.nonbonded_rep.n,
            impl.nonbonded_rep.idx_i.ptr,  impl.nonbonded_rep.idx_j.ptr,
            impl.nonbonded_rep.alpha.ptr,  impl.nonbonded_rep.repab.ptr,
            impl.nonbonded_rep.r_cut.ptr,
            impl.d_coords.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_NBREP]);
    }

    // --- Coulomb TERM 1 (pairwise erf-damped) ---
    if (m_coulomb_enabled && impl.coulomb.n > 0) {
        k_coulomb<<<gridFor(impl.coulomb.n, BLOCK), BLOCK, 0, stream>>>(
            impl.coulomb.n,
            impl.coulomb.idx_i.ptr,  impl.coulomb.idx_j.ptr,
            impl.coulomb.gamma_ij.ptr,
            impl.coulomb.r_cut.ptr,
            impl.d_coords.ptr,
            impl.d_charges.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_COUL]);
    }

    // --- Bond stretching (with CN-dependent r0, dEdcn, and HB alpha zz) ---
    if (impl.bonds.n > 0) {
        // Zero HB alpha zz buffer if HB pairs exist
        if (gradient && impl.hb_alpha.n > 0) {
            impl.d_zz_hb.zero(N, stream);
        }
        k_bonds<<<gridFor(impl.bonds.n, BLOCK), BLOCK, 0, stream>>>(
            impl.bonds.n,
            impl.bonds.idx_i.ptr,     impl.bonds.idx_j.ptr,
            impl.bonds.r0.ptr,
            impl.bonds.r0_base_i.ptr, impl.bonds.r0_base_j.ptr,
            impl.bonds.cnfak_i.ptr,   impl.bonds.cnfak_j.ptr,
            impl.bonds.rabshift.ptr,
            impl.bonds.ff.ptr,
            impl.bonds.fc.ptr,
            impl.bonds.alpha.ptr,
            impl.bonds.nr_hb.ptr,
            impl.bonds.hb_cn_H.ptr,
            impl.bonds.hb_H_atom.ptr,
            impl.d_coords.ptr,
            impl.d_cn.ptr,
            impl.d_grad.ptr,
            impl.d_dEdcn.ptr,
            (gradient && impl.hb_alpha.n > 0) ? impl.d_zz_hb.ptr : nullptr,
            &impl.d_energies.ptr[impl.E_BOND]);
    }

    // --- Angle bending (cosine + distance damping) ---
    if (impl.angles.n > 0) {
        k_angles<<<gridFor(impl.angles.n, BLOCK), BLOCK, 0, stream>>>(
            impl.angles.n,
            impl.angles.idx_i.ptr,  impl.angles.idx_j.ptr,  impl.angles.idx_k.ptr,
            impl.angles.ati.ptr,    impl.angles.atj.ptr,    impl.angles.atk.ptr,
            impl.angles.fc.ptr,     impl.angles.theta0.ptr,
            impl.d_coords.ptr,
            (double*)nullptr,       // rcov_d3 unused (constant memory)
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_ANGLE]);
    }

    // --- Dihedrals (standard + extra, Fourier + distance damping) ---
    if (impl.dihedrals.n > 0) {
        k_dihedrals<<<gridFor(impl.dihedrals.n, BLOCK), BLOCK, 0, stream>>>(
            impl.dihedrals.n,
            impl.dihedrals.idx_i.ptr, impl.dihedrals.idx_j.ptr,
            impl.dihedrals.idx_k.ptr, impl.dihedrals.idx_l.ptr,
            impl.dihedrals.ati.ptr,   impl.dihedrals.atj.ptr,
            impl.dihedrals.atk.ptr,   impl.dihedrals.atl.ptr,
            impl.dihedrals.V.ptr,
            impl.dihedrals.phi0.ptr,
            impl.dihedrals.n_period.ptr,
            impl.dihedrals.is_nci.ptr,
            impl.d_coords.ptr,
            (double*)nullptr,         // rcov_d3 unused (constant memory)
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_DIHED]);
    }

    // --- Out-of-plane inversions ---
    if (impl.inversions.n > 0) {
        k_inversions<<<gridFor(impl.inversions.n, BLOCK), BLOCK, 0, stream>>>(
            impl.inversions.n,
            impl.inversions.idx_i.ptr, impl.inversions.idx_j.ptr,
            impl.inversions.idx_k.ptr, impl.inversions.idx_l.ptr,
            impl.inversions.ati.ptr,   impl.inversions.atj.ptr,
            impl.inversions.atk.ptr,   impl.inversions.atl.ptr,
            impl.inversions.fc.ptr,    impl.inversions.omega0.ptr,
            impl.inversions.C0.ptr,    impl.inversions.C1.ptr,    impl.inversions.C2.ptr,
            impl.inversions.potential_type.ptr,
            impl.d_coords.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_INV]);
    }

    // --- Triple bond torsions ---
    if (impl.storsions.n > 0) {
        k_storsions<<<gridFor(impl.storsions.n, BLOCK), BLOCK, 0, stream>>>(
            impl.storsions.n,
            impl.storsions.idx_i.ptr, impl.storsions.idx_j.ptr,
            impl.storsions.idx_k.ptr, impl.storsions.idx_l.ptr,
            impl.storsions.erefhalf.ptr,
            impl.d_coords.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_STORS]);
    }

    // --- Bonded ATM (3-body charge-scaled) ---
    if (impl.batm.n > 0) {
        k_batm<<<gridFor(impl.batm.n, BLOCK), BLOCK, 0, stream>>>(
            impl.batm.n,
            impl.batm.idx_i.ptr, impl.batm.idx_j.ptr, impl.batm.idx_k.ptr,
            impl.batm.zb3atm_i.ptr, impl.batm.zb3atm_j.ptr, impl.batm.zb3atm_k.ptr,
            impl.d_coords.ptr,
            impl.d_topo_charges.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_BATM]);
    }

    // --- ATM 3-body dispersion ---
    if (impl.atm.n > 0) {
        k_atm<<<gridFor(impl.atm.n, BLOCK), BLOCK, 0, stream>>>(
            impl.atm.n,
            impl.atm.idx_i.ptr, impl.atm.idx_j.ptr, impl.atm.idx_k.ptr,
            impl.atm.ati.ptr,   impl.atm.atj.ptr,   impl.atm.atk.ptr,
            impl.atm.C6_ij.ptr, impl.atm.C6_ik.ptr, impl.atm.C6_jk.ptr,
            impl.atm.s9.ptr, impl.atm.a1.ptr, impl.atm.a2.ptr, impl.atm.alp.ptr,
            impl.atm.triple_scale.ptr,
            impl.d_coords.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_ATM]);
    }

    // --- Halogen bonds (3-body A-X...B) ---
    if (impl.xbonds.n > 0) {
        k_xbonds<<<gridFor(impl.xbonds.n, BLOCK), BLOCK, 0, stream>>>(
            impl.xbonds.n,
            impl.xbonds.idx_i.ptr, impl.xbonds.idx_j.ptr, impl.xbonds.idx_k.ptr,
            impl.xbonds.elem_A.ptr, impl.xbonds.elem_B.ptr,
            impl.xbonds.q_X.ptr, impl.xbonds.q_B.ptr, impl.xbonds.acidity_X.ptr,
            impl.xbonds.r_cut.ptr,
            impl.d_coords.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_XBOND]);
    }

    // --- Hydrogen bonds (3-body A-H...B, all 4 cases) ---
    if (m_hbond_enabled && impl.hbonds.n > 0) {
        k_hbonds<<<gridFor(impl.hbonds.n, BLOCK), BLOCK, 0, stream>>>(
            impl.hbonds.n,
            impl.hbonds.idx_i.ptr, impl.hbonds.idx_j.ptr, impl.hbonds.idx_k.ptr,
            impl.hbonds.elem_A.ptr, impl.hbonds.elem_B.ptr,
            impl.hbonds.q_H.ptr, impl.hbonds.q_A.ptr, impl.hbonds.q_B.ptr,
            impl.hbonds.basicity_A.ptr, impl.hbonds.basicity_B.ptr,
            impl.hbonds.acidity_A.ptr,  impl.hbonds.acidity_B.ptr,
            impl.hbonds.r_cut.ptr,
            impl.hbonds.case_type.ptr,
            impl.hbonds.nb_B_offset.ptr, impl.hbonds.nb_B_count.ptr,
            impl.hbonds.nb_B_flat.ptr,
            impl.hbonds.acceptor_parent.ptr,
            impl.hbonds.nb_C_offset.ptr, impl.hbonds.nb_C_count.ptr,
            impl.hbonds.nb_C_flat.ptr,
            impl.hbonds.repz_B.ptr,
            impl.d_coords.ptr,
            impl.d_grad.ptr,
            &impl.d_energies.ptr[impl.E_HBOND]);
    }

    // --- HB alpha-modulation chain-rule gradient ---
    // Claude Generated (March 2026): Missing gradient for d(alpha)/d(hb_cn_H) chain rule
    // zz_hb accumulated by k_bonds, now apply d(hb_cn_H)/dx via error-function CN derivative
    if (gradient && impl.hb_alpha.n > 0 && impl.d_zz_hb.ptr) {
        constexpr double hb_kn = 27.5;  // HB CN decay constant (different from D3 kn=-7.5)
        k_hb_alpha_chainrule<<<gridFor(impl.hb_alpha.n, BLOCK), BLOCK, 0, stream>>>(
            impl.hb_alpha.n,
            impl.hb_alpha.idx_H.ptr,
            impl.hb_alpha.idx_B.ptr,
            impl.hb_alpha.rcov_sum.ptr,
            impl.d_coords.ptr,
            impl.d_zz_hb.ptr,
            impl.d_grad.ptr,
            hb_kn);
    }

    // =========================================================================
    // 4. GPU postprocess: Coulomb self-energy, dispersion dEdcn, qtmp, CN chain-rule
    // Claude Generated (March 2026): Replaces postProcessCPU for full GPU consistency
    // =========================================================================

    // Coulomb TERM 2+3 self-energy on GPU (energy only, no gradient)
    if (m_coulomb_enabled && impl.coul_self_on_gpu) {
        // E_COUL_SELF slot already zeroed with d_energies above
        k_coulomb_self<<<(N+255)/256, 256, 0, impl.stream>>>(
            N,
            impl.d_charges.ptr,        // eeq_charges
            impl.d_coul_chi_base.ptr,
            impl.d_coul_cnf.ptr,
            impl.d_cn.ptr,
            impl.d_coul_gam.ptr,
            impl.d_coul_alp.ptr,
            &impl.d_energies.ptr[impl.E_COUL_SELF]);
    }

    // Dispersion dEdcn is now computed directly in k_dispersion kernel on GPU
    // (via dc6dcn_ij/dc6dcn_ji SoA uploaded by updateDispersionDC6DCN).
    // No CPU round-trip needed.

    // Snapshot dEdcn BEFORE qtmp modification (async device-to-device, no pipeline stall)
    // Claude Generated (March 2026): Replaced sync+download with D2D copy
    if (gradient) {
        cudaMemcpyAsync(impl.d_dEdcn_snapshot.ptr, impl.d_dEdcn.ptr,
                        N * sizeof(double), cudaMemcpyDeviceToDevice, stream);
    }

    // Subtract qtmp from dEdcn on GPU (Coulomb TERM 1b correction)
    if (gradient && m_coulomb_enabled && m_cnf.size() == N && m_eeq_charges.size() == N) {
        k_subtract_qtmp<<<(N+255)/256, 256, 0, impl.stream>>>(
            N,
            impl.d_charges.ptr,     // eeq_charges
            impl.d_coul_cnf.ptr,    // cnf
            impl.d_cn.ptr,          // cn
            impl.d_dEdcn.ptr);      // modified in-place
    }

    // CN chain-rule gradient on GPU (replaces dcn * dEdcn_combined)
    if (gradient && m_cn_pairs_on_gpu && impl.n_cn_pairs > 0 && impl.d_dlogdcn.ptr) {
        // Snapshot gradient BEFORE CN chain-rule (async device-to-device, no pipeline stall)
        // Claude Generated (March 2026): Replaced sync+download with D2D copy
        cudaMemcpyAsync(impl.d_grad_snapshot.ptr, impl.d_grad.ptr,
                        3*N*sizeof(double), cudaMemcpyDeviceToDevice, stream);

        k_cn_chainrule<<<(impl.n_cn_pairs + 255)/256, 256, 0, impl.stream>>>(
            impl.n_cn_pairs,
            impl.d_cn_idx_i.ptr,
            impl.d_cn_idx_j.ptr,
            impl.d_cn_rcov_sum.ptr,
            impl.d_coords.ptr,
            impl.d_dEdcn.ptr,       // combined dEdcn (bond + disp - qtmp)
            impl.d_dlogdcn.ptr,
            impl.d_grad.ptr,
            m_kn);
    }

    // =========================================================================
    // 5. Synchronise and download final results
    // =========================================================================
    checkCuda(cudaGetLastError(), "postprocess kernel launch check");
    checkCuda(cudaDeviceSynchronize(), "postprocess device sync");

    // Download all per-term GPU energies in single transfer
    // Claude Generated (March 2026): Replaces 14 individual cudaMemcpy with 1
    double h_energies[FFWorkspaceGPUImpl::N_ENERGY_SLOTS] = {};
    checkCuda(cudaMemcpy(h_energies, impl.d_energies.ptr,
                         FFWorkspaceGPUImpl::N_ENERGY_SLOTS * sizeof(double),
                         cudaMemcpyDeviceToHost), "per-term energy download");
    const double e_disp      = h_energies[impl.E_DISP];
    const double e_brep      = h_energies[impl.E_BREP];
    const double e_nbrep     = h_energies[impl.E_NBREP];
    const double e_coul1     = h_energies[impl.E_COUL];
    const double e_bond      = h_energies[impl.E_BOND];
    const double e_angle     = h_energies[impl.E_ANGLE];
    const double e_dihedral  = h_energies[impl.E_DIHED];
    const double e_inversion = h_energies[impl.E_INV];
    const double e_stors     = h_energies[impl.E_STORS];
    const double e_batm      = h_energies[impl.E_BATM];
    const double e_atm       = h_energies[impl.E_ATM];
    const double e_xbond     = h_energies[impl.E_XBOND];
    const double e_hbond     = h_energies[impl.E_HBOND];
    const double e_coul_self = h_energies[impl.E_COUL_SELF];

    // Download gradient (using pre-allocated staging buffer — no heap allocs)
    if (gradient) {
        checkCuda(cudaMemcpy(m_h_grad.data(), impl.d_grad.ptr,
                             3*N*sizeof(double), cudaMemcpyDeviceToHost),
                  "gradient download");

        for (int i = 0; i < N; ++i) {
            m_result_gradient(i, 0) = m_h_grad[3*i+0];
            m_result_gradient(i, 1) = m_h_grad[3*i+1];
            m_result_gradient(i, 2) = m_h_grad[3*i+2];
        }

        // Download diagnostic snapshots from GPU (no extra sync needed — already synced)
        // Claude Generated (March 2026): Snapshots taken via D2D copy, downloaded here
        m_dEdcn_total.resize(N);
        std::vector<double> h_dEdcn(N);
        checkCuda(cudaMemcpy(h_dEdcn.data(), impl.d_dEdcn_snapshot.ptr,
                             N * sizeof(double), cudaMemcpyDeviceToHost), "dEdcn snapshot download");
        for (int i = 0; i < N; ++i)
            m_dEdcn_total[i] = h_dEdcn[i];

        // Download pre-CN gradient snapshot
        std::vector<double> h_grad_snap(3*N);
        checkCuda(cudaMemcpy(h_grad_snap.data(), impl.d_grad_snapshot.ptr,
                             3*N*sizeof(double), cudaMemcpyDeviceToHost), "pre-cn grad snapshot download");
        m_grad_before_cn.resize(N, 3);
        for (int i = 0; i < N; ++i) {
            m_grad_before_cn(i, 0) = h_grad_snap[3*i+0];
            m_grad_before_cn(i, 1) = h_grad_snap[3*i+1];
            m_grad_before_cn(i, 2) = h_grad_snap[3*i+2];
        }
    }

    // =========================================================================
    // 6. Populate energy components
    // =========================================================================
    m_result_energy.reset();
    m_result_energy.dispersion    = e_disp;
    m_result_energy.bonded_rep    = e_brep;
    m_result_energy.nonbonded_rep = e_nbrep;
    m_result_energy.coulomb       = e_coul1 + e_coul_self;  // TERM 1 + TERM 2+3
    m_result_energy.bond          = e_bond;
    m_result_energy.angle         = e_angle;
    m_result_energy.dihedral      = e_dihedral;
    m_result_energy.inversion     = e_inversion;
    m_result_energy.stors         = e_stors;
    m_result_energy.batm          = e_batm;
    m_result_energy.atm           = e_atm;
    m_result_energy.xbond         = e_xbond;
    m_result_energy.hbond         = e_hbond;

    // =========================================================================
    // 7. Add E0 (baseline energy from parameter set) and return
    // =========================================================================
    return m_result_energy.total() + m_e0;
}

// ============================================================================
// setGeometry — extract N×3 RowMajor matrix to flat array and upload to GPU
// ============================================================================

void FFWorkspaceGPU::setGeometry(const Matrix& geom)
{
    auto& impl = *m_impl;
    const int N = m_natoms;
    if (geom.rows() != N) return;

    // Fill pre-allocated staging buffer in-place (replaces toRowMajor heap alloc)
    for (int i = 0; i < N; ++i) {
        m_h_coords[3*i+0] = geom(i, 0);
        m_h_coords[3*i+1] = geom(i, 1);
        m_h_coords[3*i+2] = geom(i, 2);
    }
    impl.d_coords.upload(m_h_coords.data(), 3*N, impl.stream);
}

// ============================================================================
// Getters
// ============================================================================

const Matrix& FFWorkspaceGPU::gradient() const
{
    // Debug: Check if result is properly sized
    if (m_result_gradient.rows() != m_natoms || m_result_gradient.cols() != 3) {
        std::cout << "[DEBUG FFWorkspaceGPU::gradient] WARNING: gradient size mismatch! "
                  << "rows=" << m_result_gradient.rows() << " cols=" << m_result_gradient.cols()
                  << " expected (" << m_natoms << ", 3)\n" << std::flush;
    }
    return m_result_gradient;
}

const FFEnergyComponents& FFWorkspaceGPU::energyComponents() const
{
    return m_result_energy;
}

const Vector& FFWorkspaceGPU::dEdcnTotal() const
{
    return m_dEdcn_total;
}

int FFWorkspaceGPU::dispersionCount() const
{
    return m_impl ? m_impl->disp.n : 0;
}

int FFWorkspaceGPU::bondCount() const
{
    return m_impl ? m_impl->bonds.n : 0;
}
