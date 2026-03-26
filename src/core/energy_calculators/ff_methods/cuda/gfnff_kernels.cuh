/*
 * <GFN-FF CUDA Kernel Declarations>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): GPU kernels for GFN-FF energy and gradient.
 * Minimum compute capability: 6.0 (Pascal) — native double atomicAdd.
 *
 * Each kernel computes one energy term, accumulates energy via atomicAdd
 * on a single double, and accumulates gradients via atomicAdd on grad[N*3].
 *
 * Thread layout: blockDim.x = 256, gridDim.x = ceil(n/256)
 */

#pragma once

#ifdef USE_CUDA

#include <cuda_runtime.h>
#include "gfnff_soa.h"

// ============================================================================
// Pairwise kernels: 1 thread = 1 pair
// ============================================================================

/// GFN-FF D4 dispersion (BJ-damped, GFN-FF modified formula)
/// E = -C6 * zetac6 * (t6 + 2*r4r2ij*t8)  where t6=1/(r6+R06), t8=1/(r8+R08)
/// Reference: Fortran gfnff_gdisp0.f90:365-377
__global__ void k_dispersion(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ C6,
    const double* __restrict__ r4r2ij,
    const double* __restrict__ r0_sq,
    const double* __restrict__ zetac6,
    const double* __restrict__ r_cut,
    const double* __restrict__ dc6dcn_ij, ///< [n] dC6/dCN(i) per pair (nullptr if no dEdcn)
    const double* __restrict__ dc6dcn_ji, ///< [n] dC6/dCN(j) per pair (nullptr if no dEdcn)
    const double* __restrict__ coords,   ///< [N*3] row-major (x,y,z per atom)
    double*                    dEdcn,    ///< [N] CN chain-rule accumulator (atomicAdd, nullptr ok)
    double*                    grad,     ///< [N*3] atomicAdd
    double*                    energy    ///< [1]   atomicAdd
);

/// GFN-FF bonded or non-bonded repulsion (same formula, different params)
/// E = repab * exp(-alpha * r^1.5) / r
/// Reference: Fortran gfnff_engrad.F90:467-495 (bonded), 255-276 (nonbonded)
__global__ void k_repulsion(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ alpha,
    const double* __restrict__ repab,
    const double* __restrict__ r_cut,
    const double* __restrict__ coords,
    double*                    grad,
    double*                    energy
);

/// GFN-FF Coulomb TERM 1 (pairwise, dynamic EEQ charges)
/// E = qi * qj * erf(gamma_ij * r) / r
/// Reference: Fortran gfnff_engrad.F90:1378-1389
__global__ void k_coulomb(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ gamma_ij,
    const double* __restrict__ r_cut,
    const double* __restrict__ coords,
    const double* __restrict__ charges,  ///< [N] EEQ charges per atom
    double*                    grad,
    double*                    energy
);

// ============================================================================
// Bonded kernels: 1 thread = 1 interaction
// ============================================================================

/// GFN-FF bond stretching (exponential potential, CN-dependent r0)
/// E = fc * exp(-alpha * (r - r0)^2)
/// Reference: Fortran gfnff_engrad.F90:675-721
__global__ void k_bonds(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ r0,
    const double* __restrict__ r0_base_i,
    const double* __restrict__ r0_base_j,
    const double* __restrict__ cnfak_i,
    const double* __restrict__ cnfak_j,
    const double* __restrict__ rabshift,
    const double* __restrict__ ff,
    const double* __restrict__ fc,
    const double* __restrict__ alpha,
    const int*    __restrict__ nr_hb,      ///< HB count per bond
    const double* __restrict__ hb_cn_H,    ///< HB coordination number of H
    const int*    __restrict__ hb_H_atom,  ///< Index of H atom per bond (-1 if no HB)
    const double* __restrict__ coords,
    const double* __restrict__ cn,     ///< [N] D3 coordination numbers
    double*                    grad,
    double*                    dEdcn,  ///< [N] CN chain-rule accumulator (atomicAdd)
    double*                    zz_hb,  ///< [N] HB alpha chain-rule zz accumulator (atomicAdd, nullable)
    double*                    energy
);

/// GFN-FF angle bending (cosine + distance damping)
/// E = fc * (cos(theta) - cos(theta0))^2 * damp_ij * damp_jk
/// Reference: Fortran gfnff_engrad.F90:857-916
__global__ void k_angles(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ ati,
    const int*    __restrict__ atj,
    const int*    __restrict__ atk,
    const double* __restrict__ fc,
    const double* __restrict__ theta0,
    const double* __restrict__ coords,
    const double* __restrict__ rcov_d3, ///< [MAX_ELEM] covalent radii table (Bohr)
    double*                    grad,
    double*                    energy
);

/// GFN-FF dihedral torsions (Fourier + distance damping, standard + extra)
/// E = V * (1 + cos(n*(phi - phi0) + pi)) * damp
/// Reference: Fortran gfnff_engrad.F90:1041-1122
__global__ void k_dihedrals(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ idx_l,
    const int*    __restrict__ ati,
    const int*    __restrict__ atj,
    const int*    __restrict__ atk,
    const int*    __restrict__ atl,
    const double* __restrict__ V,
    const double* __restrict__ phi0,
    const double* __restrict__ n_period,
    const int*    __restrict__ is_nci,
    const double* __restrict__ coords,
    const double* __restrict__ rcov_d3,
    double*                    grad,
    double*                    energy
);

/// GFN-FF out-of-plane inversions
/// E = fc * (1 - cos(omega)) * damp  or  fc * (cos(omega) - cos(omega0))^2 * damp
/// Reference: Fortran gfnff_ini.f90 inversion potential
__global__ void k_inversions(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ idx_l,
    const int*    __restrict__ ati,
    const int*    __restrict__ atj,
    const int*    __restrict__ atk,
    const int*    __restrict__ atl,
    const double* __restrict__ fc,
    const double* __restrict__ omega0,
    const double* __restrict__ C0,
    const double* __restrict__ C1,
    const double* __restrict__ C2,
    const int*    __restrict__ potential_type,
    const double* __restrict__ coords,
    double*                    grad,
    double*                    energy
);

// ============================================================================
// Phase 2 kernels: CPU residual terms ported to GPU (March 2026)
// ============================================================================

/// Triple bond torsions: E = -erefhalf * cos(2φ) + erefhalf
/// Reference: ff_workspace_gfnff.cpp:calcSTorsions, Fortran gfnff_engrad.F90:3454
__global__ void k_storsions(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ idx_l,
    const double* __restrict__ erefhalf,
    const double* __restrict__ coords,
    double*                    grad,
    double*                    energy
);

/// Bonded ATM (BATM): 3-body charge-scaled angular term for 1,4-pairs
/// Reference: ff_workspace_gfnff.cpp:calcBATM, Fortran gfnff_engrad.F90:3267-3334
__global__ void k_batm(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const double* __restrict__ zb3atm_i,
    const double* __restrict__ zb3atm_j,
    const double* __restrict__ zb3atm_k,
    const double* __restrict__ coords,
    const double* __restrict__ topo_charges, ///< [N] Phase-1 topology charges
    double*                    grad,
    double*                    energy
);

/// ATM (Axilrod-Teller-Muto): 3-body dispersion with BJ damping (energy + gradient)
/// Reference: ff_workspace_gfnff.cpp:calcATM+calcATMGradient
__global__ void k_atm(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ ati,
    const int*    __restrict__ atj,
    const int*    __restrict__ atk,
    const double* __restrict__ C6_ij,
    const double* __restrict__ C6_ik,
    const double* __restrict__ C6_jk,
    const double* __restrict__ s9,
    const double* __restrict__ a1,
    const double* __restrict__ a2,
    const double* __restrict__ alp,
    const double* __restrict__ triple_scale,
    const double* __restrict__ coords,
    double*                    grad,
    double*                    energy
);

/// Halogen bonds (3-body A-X...B): distance damped electrostatic
/// Reference: ff_workspace_gfnff.cpp:calcHalogenBonds, Fortran rbxgfnff_eg
__global__ void k_xbonds(
    int n,
    const int*    __restrict__ idx_i,     ///< donor A
    const int*    __restrict__ idx_j,     ///< halogen X
    const int*    __restrict__ idx_k,     ///< acceptor B
    const int*    __restrict__ elem_A,
    const int*    __restrict__ elem_B,
    const double* __restrict__ q_X,
    const double* __restrict__ q_B,
    const double* __restrict__ acidity_X,
    const double* __restrict__ r_cut,
    const double* __restrict__ coords,
    double*                    grad,
    double*                    energy
);

/// Hydrogen bonds (3-body A-H...B): multi-case with neighbor damping
/// Reference: ff_workspace_gfnff.cpp:calcHydrogenBonds, Fortran abhgfnff_eg*
__global__ void k_hbonds(
    int n,
    const int*    __restrict__ idx_i,     ///< donor A
    const int*    __restrict__ idx_j,     ///< hydrogen H
    const int*    __restrict__ idx_k,     ///< acceptor B
    const int*    __restrict__ elem_A,
    const int*    __restrict__ elem_B,
    const double* __restrict__ q_H,
    const double* __restrict__ q_A,
    const double* __restrict__ q_B,
    const double* __restrict__ basicity_A,
    const double* __restrict__ basicity_B,
    const double* __restrict__ acidity_A,
    const double* __restrict__ acidity_B,
    const double* __restrict__ r_cut,
    const int*    __restrict__ case_type,
    const int*    __restrict__ nb_B_offset,
    const int*    __restrict__ nb_B_count,
    const int*    __restrict__ nb_B_flat,
    const int*    __restrict__ acceptor_parent,
    const int*    __restrict__ nb_C_offset,
    const int*    __restrict__ nb_C_count,
    const int*    __restrict__ nb_C_flat,
    const double* __restrict__ repz_B,
    const double* __restrict__ coords,
    double*                    grad,
    double*                    energy
);

// ============================================================================
// GPU-only postprocess kernels (replace postProcessCPU)
// Claude Generated (March 2026): Full GPU gradient consistency
// ============================================================================

/// Coulomb TERM 2+3 self-energy (O(N), energy only — no gradient contribution)
/// E_en = -Σ qi * chi_eff_i,  E_self = 0.5 * Σ qi² * (gam_i + sqrt(2/pi)/sqrt(alp_i))
/// Reference: ff_workspace.cpp::postProcess() lines 354-383
__global__ void k_coulomb_self(
    int N,
    const double* __restrict__ eeq_charges,   ///< [N] dynamic EEQ charges
    const double* __restrict__ chi_base,       ///< [N] base electronegativity
    const double* __restrict__ cnf,            ///< [N] CN-dependent chi correction
    const double* __restrict__ cn,             ///< [N] coordination numbers
    const double* __restrict__ gam,            ///< [N] chemical hardness
    const double* __restrict__ alp,            ///< [N] Coulomb exponent
    double*                    energy           ///< [1] atomicAdd
);

/// Subtract qtmp from dEdcn in-place: dEdcn[i] -= q[i]*cnf[i]/(2*sqrt(cn[i])+eps)
/// This implements Coulomb TERM 1b chain-rule correction.
/// Reference: ff_workspace.cpp::postProcess() lines 391-399
__global__ void k_subtract_qtmp(
    int N,
    const double* __restrict__ eeq_charges,
    const double* __restrict__ cnf,
    const double* __restrict__ cn,
    double*                    dEdcn            ///< [N] modified in-place
);

/// CN chain-rule gradient: pairwise kernel over CN-relevant atom pairs.
/// For each pair (i,j): grad_i += fac*(ri-rj), grad_j -= fac*(ri-rj)
/// where fac = dS/dr / rij * (dEdcn[i]*dlogdcn[i] + dEdcn[j]*dlogdcn[j])
/// Replaces: dcn[dim] * dEdcn_combined sparse matrix-vector multiply
/// Reference: gfnff_method.cpp:calculateCoordinationNumberDerivatives (CN formula)
__global__ void k_cn_chainrule(
    int n_pairs,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ rcov_sum,       ///< [n_pairs] scaled cov. radius sum
    const double* __restrict__ coords,
    const double* __restrict__ dEdcn,          ///< [N] combined dE/dCN (bond+disp-qtmp)
    const double* __restrict__ dlogdcn,        ///< [N] logistic squashing factor
    double*                    grad,
    double        kn                            ///< CN decay constant (-7.5)
);

/// HB alpha-modulation chain-rule gradient: pairwise kernel over (H, B) pairs.
/// For each pair: grad_H += zz_H * dS/dr * (rH-rB)/r, grad_B -= same
/// where dS/dr = (-kn / (rcov*sqrt(pi))) * exp(-(kn*dr)^2), dr=(r-rcov)/rcov
/// Reference: Fortran gfnff_engrad.F90:1054-1063
__global__ void k_hb_alpha_chainrule(
    int n_pairs,
    const int*    __restrict__ idx_H,
    const int*    __restrict__ idx_B,
    const double* __restrict__ rcov_sum,       ///< [n_pairs] scaled cov. radius sum
    const double* __restrict__ coords,
    const double* __restrict__ zz_hb,          ///< [N] per-H zz accumulator from k_bonds
    double*                    grad,
    double        kn                            ///< HB CN decay constant (27.5)
);

// ============================================================================
// CN Compute kernel - GFN-FF coordination number calculation
// Claude Generated (March 2026): GPU implementation for Phase 1 GPU migration
// ============================================================================

/// Compute GFN-FF coordination numbers on GPU.
/// CN_raw[i] = sum_{j≠i} 0.5 * (1 + erf(kn * (r_ij - rcov_ij) / rcov_ij))
/// CN_final[i] = log(1 + e^cnmax) - log(1 + e^(cnmax - CN_raw[i]))
/// Thread layout: 1 thread per atom, each loops over all other atoms.
/// Uses constant memory d_rcov_d3 for covalent radii (uploaded via upload_rcov_d3)
/// Reference: gfnff_cn.f90:66-126, Spicher & Grimme J. Chem. Theory Comput. 2020
__global__ void k_cn_compute(
    int natoms,
    const double* __restrict__ coords,      ///< [3*N] in Bohr
    const int*    __restrict__ atom_types,   ///< [N] 1-based atomic numbers
    double*       __restrict__ cn_raw,      ///< [N] output: raw CN (erf sum)
    double*       __restrict__ cn_final,    ///< [N] output: log-transformed CN
    double        kn,                        ///< decay constant (-7.5)
    double        cnmax,                     ///< squashing limit (4.4)
    double        threshold_sq);             ///< distance cutoff squared (Bohr²)

// ============================================================================
// DC6DCN per-pair kernel (Phase 2 GPU optimization)
// Claude Generated (March 2026): Compute dc6/dcn directly per dispersion pair,
// eliminating O(N²) CPU matrix construction.
// ============================================================================

/// Compute dC6(i,j)/dCN(i) and dC6(i,j)/dCN(j) for each dispersion pair.
/// Formula: dc6dcn(i,j) = Σ_{ri,rj} dgw(i,ri) * gw(j,rj) * C6_ref(Zi,Zj,ri,rj)
/// Thread layout: 1 thread per dispersion pair.
/// Reference: d4param_generator.cpp:computeDC6DCN(), Fortran gfnff_gdisp0.f90:262-305
__global__ void k_dc6dcn_per_pair(
    int n_pairs,
    const int*    __restrict__ idx_i,          ///< [n] pair atom i indices
    const int*    __restrict__ idx_j,          ///< [n] pair atom j indices
    const int*    __restrict__ atom_types,     ///< [N] atomic numbers (1-based)
    const double* __restrict__ gw,             ///< [N * MAX_REF] Gaussian weights (padded)
    const double* __restrict__ dgw,            ///< [N * MAX_REF] weight derivatives
    const double* __restrict__ c6_flat,        ///< [MAX_ELEM² * MAX_REF²] C6 reference table
    const int*    __restrict__ refn,           ///< [MAX_ELEM] nref per element
    double*       __restrict__ dc6dcn_ij,      ///< [n] output: dC6(i,j)/dCN(i)
    double*       __restrict__ dc6dcn_ji       ///< [n] output: dC6(i,j)/dCN(j)
);

// ============================================================================
// Utility: zero device array
// ============================================================================
__global__ void k_zero_double(double* arr, int n);

// ============================================================================
// Upload covalent radii table to GPU constant memory
// ============================================================================
void upload_rcov_d3(const double* rcov, int n);

// ============================================================================
// Upload covalent radii table for HB/XB vdW radii to GPU constant memory
// ============================================================================
void upload_covalent_radii(const double* radii, int n);

#endif // USE_CUDA
