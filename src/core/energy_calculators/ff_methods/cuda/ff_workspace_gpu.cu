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
#include <Eigen/Sparse>

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
void BondSoA::upload(const std::vector<Bond>& v, cudaStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    std::vector<int>    h_i(n), h_j(n);
    std::vector<double> h_r0(n), h_fc(n), h_alpha(n), h_rabshift(n), h_ff(n);
    std::vector<double> h_cnfak_i(n), h_cnfak_j(n), h_rb_i(n), h_rb_j(n);

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
void InversionSoA::upload(const std::vector<Inversion>& v, cudaStream_t stream)
{
    n = static_cast<int>(v.size());
    if (n == 0) return;

    std::vector<int>    h_i(n), h_j(n), h_k(n), h_l(n), h_pt(n);
    std::vector<double> h_fc(n), h_om(n), h_C0(n), h_C1(n), h_C2(n);

    for (int m = 0; m < n; ++m) {
        h_i[m]  = v[m].i;  h_j[m] = v[m].j;
        h_k[m]  = v[m].k;  h_l[m] = v[m].l;
        h_fc[m] = v[m].fc;
        h_om[m] = v[m].omega0;
        h_C0[m] = v[m].C0;
        h_C1[m] = v[m].C1;
        h_C2[m] = v[m].C2;
        h_pt[m] = v[m].potential_type;
    }

    idx_i.upload(h_i, stream); idx_j.upload(h_j, stream);
    idx_k.upload(h_k, stream); idx_l.upload(h_l, stream);
    fc.upload(h_fc, stream);
    omega0.upload(h_om, stream);
    C0.upload(h_C0, stream);
    C1.upload(h_C1, stream);
    C2.upload(h_C2, stream);
    potential_type.upload(h_pt, stream);
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

    // Dynamic buffers (reset + upload each calculation step)
    CudaBuffer<double> d_coords;    ///< [N*3] row-major (x,y,z per atom)
    CudaBuffer<double> d_charges;   ///< [N] EEQ charges
    CudaBuffer<double> d_cn;        ///< [N] D3 coordination numbers
    CudaBuffer<double> d_grad;      ///< [N*3] gradient output (atomicAdd)
    CudaBuffer<double> d_dEdcn;     ///< [N] bond dE/dCN (atomicAdd)
    CudaBuffer<double> d_energy;    ///< [1] total GPU energy (atomicAdd)

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

    m_impl = std::make_unique<FFWorkspaceGPUImpl>();
    m_impl->N = natoms;

    // Create dedicated CUDA stream
    checkCuda(cudaStreamCreate(&m_impl->stream), "cudaStreamCreate");
    cudaStream_t stream = m_impl->stream;

    // Upload covalent radii to constant memory (used by angle/dihedral distance damping)
    upload_rcov_d3(s_rcov_d3_bohr, 87);

    // --- Upload static SoA interaction lists ---
    m_impl->disp.upload(params.dispersions, stream);
    m_impl->bonded_rep.upload(params.bonded_repulsions, stream);
    m_impl->nonbonded_rep.upload(params.nonbonded_repulsions, stream);
    m_impl->coulomb.upload(params.coulombs, stream);
    m_impl->bonds.upload(params.bonds, stream);
    m_impl->angles.upload(params.angles, atom_types, stream);
    m_impl->dihedrals.upload(params.dihedrals, params.extra_dihedrals, atom_types, stream);
    m_impl->inversions.upload(params.inversions, stream);

    // --- Allocate dynamic per-step buffers ---
    const int N3 = 3 * natoms;
    m_impl->d_coords.alloc(N3);
    m_impl->d_charges.alloc(natoms);
    m_impl->d_cn.alloc(natoms);
    m_impl->d_grad.alloc(N3);
    m_impl->d_dEdcn.alloc(natoms);
    m_impl->d_energy.alloc(1);

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
    }

    // Store initial charges and E0 from parameter set
    if (params.eeq_charges.size() == natoms)
        m_eeq_charges = params.eeq_charges;
    m_e0 = params.e0;

    // Synchronise: wait for all uploads to finish
    checkCuda(cudaStreamSynchronize(stream), "init sync");

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format(
            "FFWorkspaceGPU: {} atoms | disp={} brep={} nbrep={} coul={} bond={} ang={} dih={} inv={}",
            natoms,
            m_impl->disp.n, m_impl->bonded_rep.n, m_impl->nonbonded_rep.n,
            m_impl->coulomb.n, m_impl->bonds.n, m_impl->angles.n,
            m_impl->dihedrals.n, m_impl->inversions.n));
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
    if (m_cpu_residual) m_cpu_residual->setEEQCharges(q);
}

void FFWorkspaceGPU::setTopologyCharges(const Vector& q)
{
    m_topology_charges = q;
    if (m_cpu_residual) m_cpu_residual->setTopologyCharges(q);
}

void FFWorkspaceGPU::setD3CN(const Vector& cn)
{
    // Store D3 CN for bond r0 calculation on GPU
    m_cn = cn;
    if (m_cpu_residual) m_cpu_residual->setD3CN(cn);
}

void FFWorkspaceGPU::setCNDerivatives(const Vector& cn, const Vector& cnf,
                                       const std::vector<SpMatrix>& dcn)
{
    m_cn  = cn;
    m_cnf = cnf;
    m_dcn = dcn;
    if (m_cpu_residual) m_cpu_residual->setCNDerivatives(cn, cnf, dcn);
}

void FFWorkspaceGPU::setDC6DCNPtr(const Matrix* ptr)
{
    m_dc6dcn_ptr = ptr;
    if (m_cpu_residual) m_cpu_residual->setDC6DCNPtr(ptr);
}

void FFWorkspaceGPU::setE0(double e0)
{
    m_e0 = e0;
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
}

void FFWorkspaceGPU::setDispersionEnabled(bool v)  { m_dispersion_enabled = v; }
void FFWorkspaceGPU::setHBondEnabled(bool v)        { m_hbond_enabled = v; }
void FFWorkspaceGPU::setRepulsionEnabled(bool v)    { m_repulsion_enabled = v; }
void FFWorkspaceGPU::setCoulombEnabled(bool v)      { m_coulomb_enabled = v; }

void FFWorkspaceGPU::setCPUResidualWorkspace(FFWorkspace* cpu_ws)
{
    m_cpu_residual = cpu_ws;
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
    // 1. Zero GPU accumulators
    // =========================================================================
    cudaMemsetAsync(impl.d_energy.ptr, 0, sizeof(double),         stream);
    cudaMemsetAsync(impl.d_grad.ptr,   0, 3*N*sizeof(double),     stream);
    cudaMemsetAsync(impl.d_dEdcn.ptr,  0, N*sizeof(double),       stream);

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

    // =========================================================================
    // 3. Launch CUDA kernels (only for non-empty lists and enabled terms)
    // =========================================================================

    // --- Dispersion (D4 modified BJ formula) ---
    if (m_dispersion_enabled && impl.disp.n > 0) {
        k_dispersion<<<gridFor(impl.disp.n, BLOCK), BLOCK, 0, stream>>>(
            impl.disp.n,
            impl.disp.idx_i,  impl.disp.idx_j,
            impl.disp.C6,     impl.disp.r4r2ij,
            impl.disp.r0_sq,  impl.disp.zetac6,
            impl.disp.r_cut,
            impl.d_coords,
            impl.d_grad,
            impl.d_energy);
    }

    // --- Bonded repulsion ---
    if (m_repulsion_enabled && impl.bonded_rep.n > 0) {
        k_repulsion<<<gridFor(impl.bonded_rep.n, BLOCK), BLOCK, 0, stream>>>(
            impl.bonded_rep.n,
            impl.bonded_rep.idx_i,  impl.bonded_rep.idx_j,
            impl.bonded_rep.alpha,  impl.bonded_rep.repab,
            impl.bonded_rep.r_cut,
            impl.d_coords,
            impl.d_grad,
            impl.d_energy);
    }

    // --- Non-bonded repulsion ---
    if (m_repulsion_enabled && impl.nonbonded_rep.n > 0) {
        k_repulsion<<<gridFor(impl.nonbonded_rep.n, BLOCK), BLOCK, 0, stream>>>(
            impl.nonbonded_rep.n,
            impl.nonbonded_rep.idx_i,  impl.nonbonded_rep.idx_j,
            impl.nonbonded_rep.alpha,  impl.nonbonded_rep.repab,
            impl.nonbonded_rep.r_cut,
            impl.d_coords,
            impl.d_grad,
            impl.d_energy);
    }

    // --- Coulomb TERM 1 (pairwise erf-damped) ---
    if (m_coulomb_enabled && impl.coulomb.n > 0) {
        k_coulomb<<<gridFor(impl.coulomb.n, BLOCK), BLOCK, 0, stream>>>(
            impl.coulomb.n,
            impl.coulomb.idx_i,  impl.coulomb.idx_j,
            impl.coulomb.gamma_ij,
            impl.coulomb.r_cut,
            impl.d_coords,
            impl.d_charges,
            impl.d_grad,
            impl.d_energy);
    }

    // --- Bond stretching (with CN-dependent r0 and dEdcn accumulation) ---
    if (impl.bonds.n > 0) {
        k_bonds<<<gridFor(impl.bonds.n, BLOCK), BLOCK, 0, stream>>>(
            impl.bonds.n,
            impl.bonds.idx_i,     impl.bonds.idx_j,
            impl.bonds.r0,
            impl.bonds.r0_base_i, impl.bonds.r0_base_j,
            impl.bonds.cnfak_i,   impl.bonds.cnfak_j,
            impl.bonds.rabshift,
            impl.bonds.ff,
            impl.bonds.fc,
            impl.bonds.alpha,
            impl.d_coords,
            impl.d_cn,
            impl.d_grad,
            impl.d_dEdcn,
            impl.d_energy);
    }

    // --- Angle bending (cosine + distance damping) ---
    if (impl.angles.n > 0) {
        k_angles<<<gridFor(impl.angles.n, BLOCK), BLOCK, 0, stream>>>(
            impl.angles.n,
            impl.angles.idx_i,  impl.angles.idx_j,  impl.angles.idx_k,
            impl.angles.ati,    impl.angles.atj,    impl.angles.atk,
            impl.angles.fc,     impl.angles.theta0,
            impl.d_coords,
            nullptr,            // rcov_d3 unused (constant memory)
            impl.d_grad,
            impl.d_energy);
    }

    // --- Dihedrals (standard + extra, Fourier + distance damping) ---
    if (impl.dihedrals.n > 0) {
        k_dihedrals<<<gridFor(impl.dihedrals.n, BLOCK), BLOCK, 0, stream>>>(
            impl.dihedrals.n,
            impl.dihedrals.idx_i, impl.dihedrals.idx_j,
            impl.dihedrals.idx_k, impl.dihedrals.idx_l,
            impl.dihedrals.ati,   impl.dihedrals.atj,
            impl.dihedrals.atk,   impl.dihedrals.atl,
            impl.dihedrals.V,
            impl.dihedrals.phi0,
            impl.dihedrals.n_period,
            impl.dihedrals.is_nci,
            impl.d_coords,
            nullptr,              // rcov_d3 unused (constant memory)
            impl.d_grad,
            impl.d_energy);
    }

    // --- Out-of-plane inversions ---
    if (impl.inversions.n > 0) {
        k_inversions<<<gridFor(impl.inversions.n, BLOCK), BLOCK, 0, stream>>>(
            impl.inversions.n,
            impl.inversions.idx_i, impl.inversions.idx_j,
            impl.inversions.idx_k, impl.inversions.idx_l,
            impl.inversions.fc,    impl.inversions.omega0,
            impl.inversions.C0,    impl.inversions.C1,    impl.inversions.C2,
            impl.inversions.potential_type,
            impl.d_coords,
            impl.d_grad,
            impl.d_energy);
    }

    // =========================================================================
    // 4. Synchronise and download results
    // =========================================================================
    // Synchronize device-wide: ensures all async kernel operations complete
    // and any kernel errors (illegal memory access, etc.) propagate as exceptions.
    checkCuda(cudaGetLastError(), "kernel launch check");
    checkCuda(cudaDeviceSynchronize(), "calculate device sync");

    // Download total GPU energy
    double e_gpu = 0.0;
    checkCuda(cudaMemcpy(&e_gpu, impl.d_energy.ptr, sizeof(double),
                         cudaMemcpyDeviceToHost), "energy download");

    // Download gradient
    if (gradient) {
        std::vector<double> h_grad(3 * N);
        checkCuda(cudaMemcpy(h_grad.data(), impl.d_grad.ptr,
                             3*N*sizeof(double), cudaMemcpyDeviceToHost),
                  "gradient download");

        m_result_gradient.resize(N, 3);
        for (int i = 0; i < N; ++i) {
            m_result_gradient(i, 0) = h_grad[3*i+0];
            m_result_gradient(i, 1) = h_grad[3*i+1];
            m_result_gradient(i, 2) = h_grad[3*i+2];
        }

        // Download dEdcn (for CN chain-rule in postProcessCPU)
        m_dEdcn_total.resize(N);
        checkCuda(cudaMemcpy(m_dEdcn_total.data(), impl.d_dEdcn.ptr,
                             N*sizeof(double), cudaMemcpyDeviceToHost),
                  "dEdcn download");
    }

    // =========================================================================
    // 5. Populate energy components and run CPU post-processing
    // =========================================================================
    m_result_energy.reset();
    // Store raw GPU total in dispersion field temporarily (all terms combined)
    // Per-term decomposition would require separate GPU accumulators.
    m_result_energy.dispersion = e_gpu;  // placeholder: GPU total in one bucket

    postProcessCPU(gradient);

    // =========================================================================
    // 6. CPU residual: H-bonds, X-bonds, ATM, BATM (Phase 1 stays on CPU)
    // =========================================================================
    if (m_cpu_residual) {
        double e_cpu = m_cpu_residual->calculate(gradient);
        const auto& cpu_comp = m_cpu_residual->energyComponents();
        m_result_energy.hbond += cpu_comp.hbond;
        m_result_energy.xbond += cpu_comp.xbond;
        m_result_energy.atm   += cpu_comp.atm;
        m_result_energy.batm  += cpu_comp.batm;
        if (gradient) {
            m_result_gradient += m_cpu_residual->gradient();
        }
    }

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

    // Extract RowMajor matrix to flat [x0,y0,z0,x1,...] array for GPU kernels
    std::vector<double> h_coords = toRowMajor(geom);
    impl.d_coords.upload(h_coords.data(), 3*N, impl.stream);

    // Forward geometry to CPU residual workspace (HB/XB/ATM/BATM)
    if (m_cpu_residual) m_cpu_residual->setGeometry(geom);
}

// ============================================================================
// postProcessCPU — mirrors FFWorkspace::postProcess()
// ============================================================================

void FFWorkspaceGPU::postProcessCPU(bool gradient)
{
    const int N = m_natoms;

    // =========================================================================
    // Coulomb TERM 2+3: Self-energy (O(N), sequential)
    // Reference: Fortran gfnff_engrad.F90:1678-1679
    // Mirrors: ff_workspace.cpp::postProcess() lines 354-383
    // =========================================================================
    if (m_coul_gam.size() == N && m_eeq_charges.size() == N && m_coulomb_enabled) {
        const double sqrt_2_over_pi = 0.797884560802865;
        const bool has_cn = (m_cn.size() == N);
        double E_en = 0.0, E_self = 0.0;

        for (int i = 0; i < N; ++i) {
            if (m_coul_alp(i) <= 0.0) continue;
            const double q = m_eeq_charges(i);
            if (std::isnan(q)) continue;

            double chi;
            if (m_coul_cnf(i) != 0.0 && has_cn) {
                chi = m_coul_chi_base(i) + m_coul_cnf(i) * std::sqrt(std::max(m_cn(i), 0.0));
            } else {
                chi = m_coul_chi_static(i);
            }
            E_en   -= q * chi;
            E_self += 0.5 * q * q * (m_coul_gam(i) + sqrt_2_over_pi / std::sqrt(m_coul_alp(i)));
        }
        m_result_energy.coulomb += E_en + E_self;

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format(
                "  Coulomb self-energy (GPU workspace): EN={:+.12f}, Self={:+.12f} Eh",
                E_en, E_self));
        }
    }

    // =========================================================================
    // dEdcn chain-rule gradient + Coulomb TERM 1b
    // Reference: Fortran gfnff_engrad.F90:418-422 (bond/disp), 449-454 (coulomb)
    // Mirrors: ff_workspace.cpp::postProcess() lines 385-422
    // =========================================================================
    if (gradient && !m_dcn.empty() && m_dcn.size() == 3) {
        // Compute TERM 1b qtmp: Coulomb gradient via CN chain-rule
        Vector qtmp = Vector::Zero(N);
        const bool has_term1b = (m_eeq_charges.size() == N &&
                                 m_cnf.size() == N &&
                                 m_cn.size() == N);
        if (has_term1b && m_coulomb_enabled) {
            for (int i = 0; i < N; ++i) {
                const double cn_i = std::max(m_cn(i), 0.0);
                qtmp(i) = m_eeq_charges(i) * m_cnf(i) / (2.0 * std::sqrt(cn_i) + 1e-16);
            }
        }

        // Combined: gradient += dcn * (dEdcn_total - qtmp)
        Vector dEdcn_combined = has_term1b
            ? (m_dEdcn_total - qtmp).eval()
            : m_dEdcn_total;

        for (int dim = 0; dim < 3; ++dim) {
            if (m_dcn[dim].rows() == N && m_dcn[dim].cols() == N) {
                m_result_gradient.col(dim) += m_dcn[dim] * dEdcn_combined;
            }
        }
    }
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
