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
 * Thread layout: dynamic block size (32-512), grid = ceil(n/blockSize)
 * Phase 6 (March 2026): Warp-level reduction + adaptive block sizing.
 */

#pragma once

#ifdef USE_CUDA

#include <cuda_runtime.h>
#include "gfnff_soa.h"

// ============================================================================
// Launch bounds for optimal GPU occupancy (Phase 6: March 2026)
//
// GFNFF_KERNEL_BOUNDS: 512 threads, min 2 blocks/SM for good occupancy
// - Targets compute capability 6.0+ (Pascal and newer)
// - Balance between register pressure and occupancy
// - Allows dynamic block sizing via getLaunchConfig() for better throughput
// - 512 threads gives better occupancy on modern GPUs (RTX 3080, A100, etc.)
// ============================================================================
#define GFNFF_KERNEL_BOUNDS __launch_bounds__(512, 2)
// Light kernel bound for simple element-wise kernels (k_build_eeq_rhs, k_zero, etc.)
#define GFNFF_KERNEL_BOUNDS_LIGHT __launch_bounds__(256, 4)

// D4 constants (matching d4param_generator.h)
#define D4_MAX_ELEM 118
#define D4_MAX_REF  7

// ============================================================================
// Pairwise kernels: 1 thread = 1 pair
// ============================================================================

/// GFN-FF D4 dispersion (BJ-damped, GFN-FF modified formula)
/// E = -C6 * zetac6 * (t6 + 2*r4r2ij*t8)  where t6=1/(r6+R06), t8=1/(r8+R08)
/// Reference: Fortran gfnff_gdisp0.f90:365-377
__global__ GFNFF_KERNEL_BOUNDS void k_dispersion(
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
    const double* __restrict__ cx,       ///< [N] x-coordinates (Bohr, SoA)
    const double* __restrict__ cy,       ///< [N] y-coordinates (Bohr, SoA)
    const double* __restrict__ cz,       ///< [N] z-coordinates (Bohr, SoA)
    double*                    dEdcn,    ///< [N] CN chain-rule accumulator (atomicAdd, nullptr ok)
    double*                    grad,     ///< [N*3] atomicAdd
    double*                    energy    ///< [1]   atomicAdd
);

/// GFN-FF bonded or non-bonded repulsion (same formula, different params)
/// E = repab * exp(-alpha * r^1.5) / r
/// Reference: Fortran gfnff_engrad.F90:467-495 (bonded), 255-276 (nonbonded)
__global__ GFNFF_KERNEL_BOUNDS void k_repulsion(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ alpha,
    const double* __restrict__ repab,
    const double* __restrict__ r_cut,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy
);

/// Mixed-precision variant: FP32 intermediates, FP64 accumulation
__global__ GFNFF_KERNEL_BOUNDS void k_repulsion_mixed(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ alpha,
    const double* __restrict__ repab,
    const double* __restrict__ r_cut,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy
);

/// GFN-FF Coulomb TERM 1 (pairwise, dynamic EEQ charges)
/// E = qi * qj * erf(gamma_ij * r) / r
/// Reference: Fortran gfnff_engrad.F90:1378-1389
__global__ GFNFF_KERNEL_BOUNDS void k_coulomb(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ gamma_ij,
    const double* __restrict__ r_cut,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
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
__global__ GFNFF_KERNEL_BOUNDS void k_bonds(
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
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ cn,     ///< [N] D3 coordination numbers
    double*                    grad,
    double*                    dEdcn,  ///< [N] CN chain-rule accumulator (atomicAdd)
    double*                    zz_hb,  ///< [N] HB alpha chain-rule zz accumulator (atomicAdd, nullable)
    double*                    energy
);

/// GFN-FF angle bending (cosine + distance damping)
/// E = fc * (cos(theta) - cos(theta0))^2 * damp_ij * damp_jk
/// Reference: Fortran gfnff_engrad.F90:857-916
__global__ GFNFF_KERNEL_BOUNDS void k_angles(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ ati,
    const int*    __restrict__ atj,
    const int*    __restrict__ atk,
    const double* __restrict__ fc,
    const double* __restrict__ theta0,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ rcov_d3, ///< [MAX_ELEM] covalent radii table (Bohr)
    double*                    grad,
    double*                    energy
);

/// GFN-FF dihedral torsions (Fourier + distance damping, standard + extra)
/// E = V * (1 + cos(n*(phi - phi0) + pi)) * damp
/// Reference: Fortran gfnff_engrad.F90:1041-1122
__global__ GFNFF_KERNEL_BOUNDS void k_dihedrals(
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
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ rcov_d3,
    double*                    grad,
    double*                    energy
);

/// GFN-FF out-of-plane inversions
/// E = fc * (1 - cos(omega)) * damp  or  fc * (cos(omega) - cos(omega0))^2 * damp
/// Reference: Fortran gfnff_ini.f90 inversion potential
__global__ GFNFF_KERNEL_BOUNDS void k_inversions(
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
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy
);

// ============================================================================
// Phase 2 kernels: CPU residual terms ported to GPU (March 2026)
// ============================================================================

/// Triple bond torsions: E = -erefhalf * cos(2φ) + erefhalf
/// Reference: ff_workspace_gfnff.cpp:calcSTorsions, Fortran gfnff_engrad.F90:3454
__global__ GFNFF_KERNEL_BOUNDS void k_storsions(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ idx_l,
    const double* __restrict__ erefhalf,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy
);

/// Bonded ATM (BATM): 3-body charge-scaled angular term for 1,4-pairs
/// Reference: ff_workspace_gfnff.cpp:calcBATM, Fortran gfnff_engrad.F90:3267-3334
__global__ GFNFF_KERNEL_BOUNDS void k_batm(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const double* __restrict__ zb3atm_i,
    const double* __restrict__ zb3atm_j,
    const double* __restrict__ zb3atm_k,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ topo_charges, ///< [N] Phase-1 topology charges
    double*                    grad,
    double*                    energy
);

/// Mixed-precision variant of k_batm: FP32 intermediates, FP64 accumulation
__global__ GFNFF_KERNEL_BOUNDS void k_batm_mixed(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const double* __restrict__ zb3atm_i,
    const double* __restrict__ zb3atm_j,
    const double* __restrict__ zb3atm_k,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ topo_charges,
    double*                    grad,
    double*                    energy
);

/// Mixed-precision variant of k_xbonds: FP32 intermediates, FP64 accumulation
__global__ GFNFF_KERNEL_BOUNDS void k_xbonds_mixed(
    int n,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const int*    __restrict__ idx_k,
    const int*    __restrict__ elem_A,
    const int*    __restrict__ elem_B,
    const double* __restrict__ q_X,
    const double* __restrict__ q_B,
    const double* __restrict__ acidity_X,
    const double* __restrict__ r_cut,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy
);

/// ATM (Axilrod-Teller-Muto): 3-body dispersion with BJ damping (energy + gradient)
/// Reference: ff_workspace_gfnff.cpp:calcATM+calcATMGradient
__global__ GFNFF_KERNEL_BOUNDS void k_atm(
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
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy
);

/// Halogen bonds (3-body A-X...B): distance damped electrostatic
/// Reference: ff_workspace_gfnff.cpp:calcHalogenBonds, Fortran rbxgfnff_eg
__global__ GFNFF_KERNEL_BOUNDS void k_xbonds(
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
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*                    grad,
    double*                    energy
);

/// Claude Generated (June 2026): device-resident HB charge refresh. Gathers the live
/// per-atom EEQ charges (q_atom) into the HB SoA arrays by donor/H/acceptor index, so the
/// next step's k_hbonds uses current charges instead of values frozen at topology build.
__global__ void k_gather_hb_charges(
    int n,
    const int* __restrict__ idx_i, const int* __restrict__ idx_j, const int* __restrict__ idx_k,
    const double* __restrict__ q_atom,
    double* __restrict__ q_A, double* __restrict__ q_H, double* __restrict__ q_B);

/// Hydrogen bonds (3-body A-H...B): multi-case with neighbor damping
/// Reference: ff_workspace_gfnff.cpp:calcHydrogenBonds, Fortran abhgfnff_eg*
__global__ GFNFF_KERNEL_BOUNDS void k_hbonds(
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
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
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
__global__ GFNFF_KERNEL_BOUNDS void k_coulomb_self(
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
__global__ GFNFF_KERNEL_BOUNDS void k_subtract_qtmp(
    int N,
    const double* __restrict__ eeq_charges,
    const double* __restrict__ cnf,
    const double* __restrict__ cn,
    double*                    dEdcn            ///< [N] modified in-place
);

/// Fused Coulomb self-energy + qtmp subtraction (single O(N) pass)
/// Replaces separate k_coulomb_self + k_subtract_qtmp launches.
/// Uses blockReduceAddEnergy for energy (fixes pattern inconsistency).
__global__ GFNFF_KERNEL_BOUNDS void k_coulomb_postprocess(
    int N,
    const double* __restrict__ eeq_charges,
    const double* __restrict__ chi_base,
    const double* __restrict__ cnf,
    const double* __restrict__ cn,
    const double* __restrict__ gam,
    const double* __restrict__ alp,
    double*                    dEdcn,
    double*                    energy,
    bool                       do_subtract
);

/// CN chain-rule gradient: pairwise kernel over CN-relevant atom pairs.
/// For each pair (i,j): grad_i += fac*(ri-rj), grad_j -= fac*(ri-rj)
/// where fac = dS/dr / rij * (dEdcn[i]*dlogdcn[i] + dEdcn[j]*dlogdcn[j])
/// Replaces: dcn[dim] * dEdcn_combined sparse matrix-vector multiply
/// Reference: gfnff_method.cpp:calculateCoordinationNumberDerivatives (CN formula)
__global__ GFNFF_KERNEL_BOUNDS void k_cn_chainrule(
    int n_pairs,
    const int*    __restrict__ idx_i,
    const int*    __restrict__ idx_j,
    const double* __restrict__ rcov_sum,       ///< [n_pairs] scaled cov. radius sum
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const double* __restrict__ dEdcn,          ///< [N] combined dE/dCN (bond+disp-qtmp)
    const double* __restrict__ dlogdcn,        ///< [N] logistic squashing factor
    double*                    grad,
    double        kn                            ///< CN decay constant (-7.5)
);

/// HB alpha-modulation chain-rule gradient: pairwise kernel over (H, B) pairs.
/// For each pair: grad_H += zz_H * dS/dr * (rH-rB)/r, grad_B -= same
/// where dS/dr = (-kn / (rcov*sqrt(pi))) * exp(-(kn*dr)^2), dr=(r-rcov)/rcov
/// Reference: Fortran gfnff_engrad.F90:1054-1063
__global__ GFNFF_KERNEL_BOUNDS void k_hb_alpha_chainrule(
    int n_pairs,
    const int*    __restrict__ idx_H,
    const int*    __restrict__ idx_B,
    const double* __restrict__ rcov_sum,       ///< [n_pairs] scaled cov. radius sum
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
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
__global__ GFNFF_KERNEL_BOUNDS void k_cn_compute(
    int natoms,
    const double* __restrict__ cx,           ///< [N] x-coordinates (Bohr, SoA)
    const double* __restrict__ cy,           ///< [N] y-coordinates (Bohr, SoA)
    const double* __restrict__ cz,           ///< [N] z-coordinates (Bohr, SoA)
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
__global__ GFNFF_KERNEL_BOUNDS void k_dc6dcn_per_pair(
    int n_pairs,
    const int*    __restrict__ idx_i,          ///< [n] pair atom i indices
    const int*    __restrict__ idx_j,          ///< [n] pair atom j indices
    const int*    __restrict__ atom_types,     ///< [N] atomic numbers (1-based)
    const double* __restrict__ gw,             ///< [N * MAX_REF] Gaussian weights (padded)
    const double* __restrict__ dgw,            ///< [N * MAX_REF] weight derivatives
    const double* __restrict__ c6_flat,        ///< [MAX_ELEM² * MAX_REF²] C6 reference table
    double*       __restrict__ dc6dcn_ij,      ///< [n] output: dC6(i,j)/dCN(i)
    double*       __restrict__ dc6dcn_ji       ///< [n] output: dC6(i,j)/dCN(j)
    // refn read from d_refn_const (constant memory, uploaded via upload_refn_const)
);

// ============================================================================
// GPU dlogdcn computation (Apr 2026)
// Compute logistic squashing derivative directly on GPU after k_cn_compute.
// dlogdcn[i] = exp(cnmax) / (exp(cnmax) + exp(cn_raw[i]))
// Eliminates CPU loop + H2D upload in gfnff_gpu_method.cpp.
// ============================================================================

__global__ GFNFF_KERNEL_BOUNDS void k_dlogdcn(
    int natoms,
    const double* __restrict__ cn_raw,  ///< [N] raw CN values (erf sum)
    double*       __restrict__ dlogdcn, ///< [N] output: logistic derivative
    double        exp_cnmax             ///< pre-computed exp(cnmax)
);

// ============================================================================
// GPU Gaussian weight computation (Phase 6: March 2026)
// Claude Generated: Compute gw and dgw/dCN on GPU, eliminating CPU computation
// + flatten + H2D upload.
// ============================================================================

/// Compute normalized Gaussian weights and their CN-derivatives per atom.
/// gw(ref) = exp(-wf*(CN-CN_ref)^2) / norm,  dgw/dCN via quotient rule.
/// Thread layout: 1 thread per atom, each loops over MAX_REF references.
/// Reference: d4param_generator.cpp:precomputeGaussianWeights() + computeGaussianWeightDerivatives()
__global__ GFNFF_KERNEL_BOUNDS void k_gaussian_weights(
    int natoms,
    const double* __restrict__ cn,           ///< [N] coordination numbers
    const int*    __restrict__ atom_types,   ///< [N] atomic numbers (1-based)
    double*       __restrict__ gw,           ///< [N * MAX_REF] output: normalized weights
    double*       __restrict__ dgw           ///< [N * MAX_REF] output: weight derivatives
    // refcn and refn read from d_refcn_const/d_refn_const (constant memory)
);

// ============================================================================
// Utility: zero device array
// ============================================================================
__global__ void k_zero_double(double* arr, int n);

// ============================================================================
// GPU topology displacement check (flag-based, March 2026)
// ============================================================================
__global__ void k_check_displacement(
    int N,
    const double* __restrict__ cx,           ///< [N] current x-coordinates (SoA)
    const double* __restrict__ cy,           ///< [N] current y-coordinates (SoA)
    const double* __restrict__ cz,           ///< [N] current z-coordinates (SoA)
    const double* __restrict__ rx,           ///< [N] reference x-coordinates (SoA)
    const double* __restrict__ ry,           ///< [N] reference y-coordinates (SoA)
    const double* __restrict__ rz,           ///< [N] reference z-coordinates (SoA)
    double        threshold_sq,
    int*          exceeded_flag);

// ============================================================================
// Upload covalent radii table to GPU constant memory
// ============================================================================
void upload_rcov_d3(const double* rcov, int n);

// ============================================================================
// Upload covalent radii table for HB/XB vdW radii to GPU constant memory
// ============================================================================
void upload_covalent_radii(const double* radii, int n);

// ============================================================================
// Upload D4 reference CN and nref tables to GPU constant memory (once at init)
// Claude Generated (March 2026): Phase 8 — constant memory broadcast for refcn/refn
// ============================================================================
void upload_refcn_const(const double* data, int n);
void upload_refn_const(const int* data, int n);

// ============================================================================
// GPU CN pair list generation (Apr 2026)
// Two-pass kernels: count valid pairs, then write (i, j, rcov_sum).
// Replaces CPU generateCNPairList() O(N^2) loop.
// ============================================================================

__global__ GFNFF_KERNEL_BOUNDS void k_generate_cn_pairs_count(
    int N,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const int*    __restrict__ atom_types,
    double        cutoff_factor,
    int*          d_count);

__global__ GFNFF_KERNEL_BOUNDS void k_generate_cn_pairs_write(
    int N,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const int*    __restrict__ atom_types,
    double        cutoff_factor,
    int*          d_counter,
    int*          idx_i,
    int*          idx_j,
    double*       rcov_sum);

// ============================================================================
// WP-A (Jun 2026): on-device D4 dispersion pair-list build.
// Two-pass enumeration (count + build) mirroring the CN pair list, plus the
// per-pair static data (C6 contraction + r4r2ij + R0² + zetac6) so the host
// O(N²) GenerateDispersionPairsNative loop + the per-build H2D upload are not
// needed when -gfnff.gpu_disp_pairs_on_device is on. Reference (bit-identical
// target): d4param_generator.cpp:GenerateDispersionPairsNative (line 1789).
// ============================================================================

/// Pass 1: count i<j pairs with r² ≤ cutoff_sq and valid elements (1..118).
__global__ GFNFF_KERNEL_BOUNDS void k_disp_pairs_count(
    int N,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const int*    __restrict__ atom_types,
    double        cutoff_sq,
    int*          d_count);

/// Pass 2: write idx_i/idx_j + per-pair static data for each surviving pair.
/// C6 = Σ_{ri,rj} gw_i·gw_j·c6ref (CN-only weighting, == getChargeWeightedC6);
/// r4r2ij = 3·sqrtZr4r2_i·sqrtZr4r2_j; r0_sq = (a1·√r4r2ij + a2)²;
/// zetac6 = zetaChargeScale(Zi,q_i)·zetaChargeScale(Zj,q_j); r_cut = 50.
__global__ GFNFF_KERNEL_BOUNDS void k_disp_pairs_build(
    int N,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    const int*    __restrict__ atom_types,
    double        cutoff_sq,
    const double* __restrict__ gw,          ///< [N*MAX_REF] device Gaussian weights
    const double* __restrict__ c6_flat,     ///< C6 reference table
    const double* __restrict__ sqrtzr4r2,   ///< [118] sqrt(Z·r4/r2) per element
    const double* __restrict__ topo_q,      ///< [N] Phase-1 topology charges
    const double* __restrict__ zeta_zeff,   ///< [86] zeta zeff per element
    const double* __restrict__ zeta_c,      ///< [86] zeta c per element
    double        a1,
    double        a2,
    int*          d_counter,
    int*          out_i,
    int*          out_j,
    double*       out_C6,
    double*       out_r4r2,
    double*       out_r0sq,
    double*       out_zeta,
    double*       out_rcut
    // refn read from d_refn_const (constant memory)
);

// ============================================================================
// GPU HB CN per-bond computation (Apr 2026)
// Replaces CPU HB CN loop in gfnff_gpu_method.cpp.
// ============================================================================

/// 1 thread = 1 (H,B) pair. Atomically accumulates into per-H buffer.
__global__ GFNFF_KERNEL_BOUNDS void k_hb_cn_per_atom(
    int n_pairs,
    const int*    __restrict__ idx_H,
    const int*    __restrict__ idx_B,
    const double* __restrict__ rcov_sum,
    const double* __restrict__ cx,
    const double* __restrict__ cy,
    const double* __restrict__ cz,
    double*       __restrict__ hb_cn_per_atom, // [N] zeroed before launch
    double        kn,                          // 27.5
    double        thr_sq);                     // 900.0

/// 1 thread = 1 bond. Reads hb_H_atom from BondSoA and scatters per-H CN into hb_cn_H.
__global__ GFNFF_KERNEL_BOUNDS void k_hb_cn_map_to_bonds(
    int nb,
    const int*    __restrict__ nr_hb,
    const int*    __restrict__ hb_H_atom,
    const double* __restrict__ hb_cn_per_atom,
    double*       __restrict__ hb_cn_H);       // [nb] output, BondSoA buffer

// ============================================================================
// GPU EEQ Schur complement (Apr 2026)
// Replaces CPU Schur loop after Cholesky solve in gfnff_gpu_method.cpp.
// For nfrag == 1 only (common case). nfrag > 1 falls back to CPU path.
// ============================================================================

/// Block-reduce S = sum(Z2_col0) and Cz1 = sum(z1) on GPU.
/// Launch with 1 block. Shared mem: 2 * blockSize * sizeof(double).
/// Writes 2 scalars to d_sums[0] = Cz1, d_sums[1] = S.
__global__ void k_eeq_reduce_sums(
    int N,
    const double* __restrict__ d_rhs,  ///< [N * nrhs] column-major, nrhs >= 2
    double*       __restrict__ d_sums); ///< [2] output: [Cz1, S]

/// Element-wise Schur complement for nfrag == 1.
/// charges[i] = z1[i] - Z2[i] * lambda
__global__ void k_eeq_schur_nfrag1(
    int N,
    const double* __restrict__ d_rhs,   ///< [N * 2] column-major: col0=z1, col1=Z2
    double        lambda,
    double*       __restrict__ charges); ///< [N] output

// ============================================================================
// WP7-A: General Schur complement for nfrag > 1 (May 2026)
// Replaces CPU Schur loop in gfnff_gpu_method.cpp for multi-fragment systems.
// k_eeq_reduce_fragment_sums: Cz1[f] = Σ_{i∈frag_f} z1[i],
//                             S[f,g] = Σ_{i∈frag_f} Z2[i,g]
// k_eeq_schur_general:        q[i] = z1[i] - Σ_g Z2[i,g] * λ[g]
// ============================================================================

/// Block-reduce per-fragment sums of z1 (column 0) and Z2[:,g] (columns 1..nfrag).
/// Each block reduces a tile of atoms in shared memory, then atomicAdd's
/// (nfrag + 1) values into d_Cz1[atom_frag[i]] and d_S[atom_frag[i]*nfrag + g].
/// Caller must zero d_Cz1 and d_S before launch (cudaMemsetAsync).
__global__ void k_eeq_reduce_fragment_sums(
    int N,
    int nfrag,
    const double* __restrict__ d_rhs,        ///< [N*(nfrag+1)] column-major
    const int*    __restrict__ d_atom_frag,  ///< [N] 0-indexed fragment id per atom
    double*       __restrict__ d_Cz1,        ///< [nfrag] output: per-fragment Σ z1
    double*       __restrict__ d_S);         ///< [nfrag*nfrag] row-major: S[f*nfrag+g]

/// Element-wise Schur apply for nfrag > 1.
/// charges[i] = z1[i] - Σ_{g=0..nfrag-1} Z2[i,g] * lambda[g]
__global__ void k_eeq_schur_general(
    int N,
    int nfrag,
    const double* __restrict__ d_rhs,        ///< [N*(nfrag+1)] column-major
    const double* __restrict__ d_lambda,     ///< [nfrag] Lagrange multipliers
    double*       __restrict__ charges);     ///< [N] output (in-place on d_rhs[0..N-1] OK)

// ============================================================================
// WP7-C: GPU PCG kernels (May 2026)
// 1-thread-per-atom helpers used by EEQSolverGPU::solveSinglePCG.
// Matvec / dot / axpy go via cuBLAS (cublasDsymv / Ddot / Daxpy).
// ============================================================================

/// Extract A's diagonal and invert it: M_inv[i] = 1/A[i,i].
/// A is column-major N×N (k_eeq_build_matrix output, both triangles filled).
__global__ void k_pcg_extract_diag_inv(
    int N,
    const double* __restrict__ d_A,
    double*       __restrict__ d_M_inv);

/// Compute initial residual r = b − A·x.
/// d_Ax must hold A·x (precomputed via cublasDsymv).
__global__ void k_pcg_init_residual(
    int N,
    const double* __restrict__ d_b,
    const double* __restrict__ d_Ax,
    double*       __restrict__ d_r);

/// Apply Jacobi preconditioner: z = M_inv ⊙ r (elementwise multiply).
__global__ void k_pcg_apply_precond(
    int N,
    const double* __restrict__ d_M_inv,
    const double* __restrict__ d_r,
    double*       __restrict__ d_z);

/// PCG direction update: p_out = z + β·p_in. Single fused kernel.
__global__ void k_pcg_dir_update(
    int N,
    const double* __restrict__ d_z,
    double        beta,
    const double* __restrict__ d_p_in,
    double*       __restrict__ d_p_out);

// ============================================================================
// WP7-D: GPU block-Jacobi preconditioner for the PCG (Jun 2026)
// Port of the CPU EEQSolver::BlockJacobiPC (per-fragment exact inverse). Replaces
// the diagonal Jacobi z = M_inv⊙r with z = blockdiag(A_ff^-1)·r for nfrag>=2,
// cutting PCG iterations from ~30-100 to ~2-5 on many-fragment (solvent) systems.
// Reference: eeq_solver.cpp BlockJacobiPC::apply / buildBlockJacobi.
// ============================================================================

/// Symmetrize each per-fragment block in place (mirror lower → upper triangle).
/// cusolverDnDpotri(LOWER) leaves the upper triangle untouched; the block-Jacobi
/// apply does a full GEMV, so the upper triangle must be filled. One block / fragment.
__global__ void k_eeq_symmetrize_blocks(
    int           nfrag,
    double*       __restrict__ d_A_blocks,    ///< [sum N_f²] explicit inverses, column-major per block
    const int*    __restrict__ frag_sizes,    ///< [nfrag] N_f
    const int*    __restrict__ frag_offsets_A);///< [nfrag] start of fragment f's N_f² block

/// Apply the block-Jacobi preconditioner: z = blockdiag(A_ff^-1)·r.
/// One thread-block per fragment; r/z stay in global atom order, gather/scatter via
/// frag_atom_map. d_Ainv_blocks holds the (symmetric) per-fragment inverse blocks.
/// Dynamic shared memory: m_max_frag_N doubles (staged r_f).
__global__ void k_eeq_block_jacobi_apply(
    int           nfrag,
    const double* __restrict__ d_Ainv_blocks, ///< [sum N_f²] symmetric inverses, column-major per block
    const int*    __restrict__ frag_sizes,    ///< [nfrag] N_f
    const int*    __restrict__ frag_offsets_A,///< [nfrag] start of fragment f's N_f² block
    const int*    __restrict__ frag_atom_offsets,///< [nfrag+1] sorted-position start per fragment
    const int*    __restrict__ frag_atom_map, ///< [N] sorted-position → global atom index
    const double* __restrict__ d_r,           ///< [N] residual (global order)
    double*       __restrict__ d_z);          ///< [N] preconditioned residual (global order)

// ============================================================================
// WP2: GPU-side EEQ RHS construction
// ============================================================================

/// Build EEQ RHS vector on GPU: rhs[i] = chi_corr[i] + cnf[i]*sqrt(max(cn[i],0))
/// Launch after k_cn_compute on same stream — eliminates finalizeCNForCPU() sync
/// for the EEQ RHS. chi_corr and cnf are topology-constant (upload once per topo).
/// Reference: gfnff_method.cpp prepareEEQParametersForGPU lines 1029-1037
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_build_eeq_rhs(
    int N,
    const double* __restrict__ d_cn,       ///< [N] log-transformed CN (d_cn_final)
    const double* __restrict__ d_chi_corr, ///< [N] -chi+dxi+amide_corr (topology-const)
    const double* __restrict__ d_cnf,      ///< [N] cnf_eeq per atom (topology-const)
    double*       __restrict__ d_rhs       ///< [N] output RHS
);

// ============================================================================
// WP5-C: GPU-side D4 dc6dcn skip check
// ============================================================================

/// Stage 1: per-block max reduction of |cn_cur - cn_ref| → d_block_max[]
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_check_dc6dcn_skip(
    int N,
    const double* __restrict__ cn_cur,
    const double* __restrict__ cn_ref,
    double* __restrict__ d_block_max);

/// Stage 2: single-block reduction of d_block_max[] → skip_flag (0 or 1)
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_check_dc6dcn_skip_final(
    int n_blocks,
    const double* __restrict__ d_block_max,
    double abs_threshold,
    int* __restrict__ skip_flag);

// ============================================================================
// WP3: Pair-list-based CN computation — O(N²) → O(n_pairs)
// ============================================================================

/// Compute raw GFN-FF CN from pre-built pair list: 1 thread per (i,j) pair.
/// atomicAdd into d_cn_raw[i] and d_cn_raw[j] — caller must zero d_cn_raw first.
/// Uses same erf formula as k_cn_compute: contrib = 0.5*(1+erf(kn*(r/rcov-1))).
/// Reuses CN pair list from generateCNPairListOnGPU() (built once per topology).
/// Reference: gfnff_cn.f90:66-126
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_cn_compute_pairs(
    int           n_pairs,
    const int*    __restrict__ idx_i,    ///< [n_pairs] atom i index
    const int*    __restrict__ idx_j,    ///< [n_pairs] atom j index
    const double* __restrict__ rcov_sum, ///< [n_pairs] sum of 4/3-scaled covalent radii
    const double* __restrict__ cx,       ///< [N] x-coordinates (Bohr, SoA)
    const double* __restrict__ cy,       ///< [N] y-coordinates (Bohr, SoA)
    const double* __restrict__ cz,       ///< [N] z-coordinates (Bohr, SoA)
    double*       __restrict__ cn_raw,   ///< [N] output: erf sum (atomicAdd, zero before)
    double        kn                     ///< CN decay constant (-7.5)
);

/// Apply log squashing transform from cn_raw to cn_final: 1 thread per atom.
/// cn_final[i] = log(1+exp(cnmax)) - log(1+exp(cnmax-cn_raw[i]))
/// Called after k_cn_compute_pairs; replaces the embedded transform in k_cn_compute.
/// Reference: gfnff_cn.f90:93-96
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_logcn(
    int natoms,
    const double* __restrict__ cn_raw,   ///< [N] raw erf-sum CN
    double*       __restrict__ cn_final, ///< [N] output: log-squashed CN
    double        cnmax                  ///< squashing limit (4.4)
);

// ============================================================================
// WP6: Batched per-fragment EEQ kernels (nfrag > 1)
// Claude Generated (May 2026): Independent N_f×N_f Coulomb blocks per fragment.
// Cross-fragment Coulomb set to zero (valid for well-separated fragments in MD).
// ============================================================================

/**
 * @brief Build independent per-fragment N_f×N_f Coulomb blocks for batched EEQ.
 *
 * 1 thread per lower-triangle element across ALL fragments (packed).
 * total_pairs = sum_f N_f*(N_f+1)/2 total threads.
 * Thread maps to (fragment_f, local_i, local_j) via binary search on frag_offsets_pair.
 * Writes column-major block into d_A_blocks[frag_offsets_A[f] + local_j*N_f + local_i].
 *
 * Reference: same erf(gamma*r)/r formula as k_eeq_build_matrix, intra-fragment only.
 */
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_eeq_build_fragment_matrices(
    int           total_pairs,               ///< sum_f N_f*(N_f+1)/2
    const double* __restrict__ cx,           ///< [N] x-coordinates (Bohr, SoA, global order)
    const double* __restrict__ cy,           ///< [N] y-coordinates
    const double* __restrict__ cz,           ///< [N] z-coordinates
    const double* __restrict__ alpha,        ///< [N] alpha² per atom (global order)
    const double* __restrict__ gam,          ///< [N] gam_corrected per atom (global order)
    const int*    __restrict__ frag_sizes,       ///< [nfrag] N_f per fragment
    const int*    __restrict__ frag_offsets_A,   ///< [nfrag] start in d_A_blocks
    const int*    __restrict__ frag_offsets_pair,///< [nfrag+1] prefix sum N_f*(N_f+1)/2
    const int*    __restrict__ frag_atom_offsets,///< [nfrag+1] prefix sum N_f (atom start)
    const int*    __restrict__ frag_atom_map,    ///< [N] sorted-position → global atom index
    double*       __restrict__ d_A_blocks,       ///< [sum N_f²] output, column-major per block
    int           nfrag,
    double        cutoff_sq                      ///< distance cutoff² (0 = no cutoff)
);

/**
 * @brief Gather per-atom RHS (b_atoms) from global order into per-fragment RHS blocks.
 *
 * 1 thread per sorted atom position k (0..N-1).
 * Reads d_rhs_global[frag_atom_map[k]] and writes to d_rhs_blocks col0.
 * Col1 (constraint = all-ones per fragment) is pre-filled at uploadFragmentTopology.
 */
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_eeq_gather_rhs_fragments(
    int           N,
    const double* __restrict__ d_rhs_global,      ///< [N] from k_build_eeq_rhs (global order)
    const int*    __restrict__ frag_atom_map,      ///< [N] sorted-position → global index
    const int*    __restrict__ frag_atom_offsets,  ///< [nfrag+1] prefix sum N_f
    const int*    __restrict__ frag_offsets_rhs,   ///< [nfrag] start in d_rhs_blocks
    const int*    __restrict__ frag_sizes,         ///< [nfrag] N_f per fragment
    double*       __restrict__ d_rhs_blocks,       ///< [sum N_f*2] col0 written here
    int           nfrag
);

#endif // USE_CUDA
