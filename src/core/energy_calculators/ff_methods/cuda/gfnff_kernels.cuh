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
    const double* __restrict__ coords,   ///< [N*3] row-major (x,y,z per atom)
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
    const double* __restrict__ coords,
    const double* __restrict__ cn,     ///< [N] D3 coordination numbers
    double*                    grad,
    double*                    dEdcn,  ///< [N] CN chain-rule accumulator (atomicAdd)
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
// Utility: zero device array
// ============================================================================
__global__ void k_zero_double(double* arr, int n);

// ============================================================================
// Upload covalent radii table to GPU constant memory
// ============================================================================
void upload_rcov_d3(const double* rcov, int n);

#endif // USE_CUDA
