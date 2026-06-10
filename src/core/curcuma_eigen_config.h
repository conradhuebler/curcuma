/**
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * curcuma_eigen_config.h — Centralized Eigen BLAS/LAPACK activation
 *
 * This file MUST be included BEFORE any Eigen headers. It is pulled in by
 * pch_base.h (precompiled header) and global.h, so Eigen optimizations are
 * active for the entire build without per-file boilerplate.
 *
 * CMake sets USE_BLAS / USE_MKL as target compile definitions; this header
 * translates them into the Eigen macros.
 *
 * Note: EIGEN_USE_LAPACKE is deliberately NOT defined here. It is applied
 * per-file in selected TUs (rf_solver.cpp, lbfgs.cpp) to avoid Eigen 3.4
 * complex-type conflicts that break other translation units.
 */

#ifndef CURCUMA_EIGEN_CONFIG_H
#define CURCUMA_EIGEN_CONFIG_H

#ifdef USE_BLAS
#    define EIGEN_USE_BLAS
#endif

#ifdef USE_MKL
// EIGEN_USE_MKL_ALL disabled — incompatible with Eigen 3.4 (mkl_direct_call.h
// void* → double* cast errors). MKL provides standard BLAS, so
// EIGEN_USE_BLAS suffices for Eigen to call MKL's BLAS routines.
#    define EIGEN_USE_BLAS
#endif

#endif // CURCUMA_EIGEN_CONFIG_H