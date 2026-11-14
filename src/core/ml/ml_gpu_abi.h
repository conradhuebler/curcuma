/*
 * < Plain-C interface between the Delta-ML host code and its GPU plugin >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (October 2026): device-resident ACSF descriptor, kernel and force
 * contraction. Shared by ml_gpu.cpp (host, compiled into curcuma_core) and the CUDA plugin;
 * plain old data only, so both compilers lay it out identically.
 */

#pragma once

/* Upper bounds of the descriptor grid the device kernels accept. A model outside them is
 * refused by the host wrapper with a message; it never runs a truncated grid. */
#define CURCUMA_ML_MAX_RADIAL 32
#define CURCUMA_ML_MAX_ZETA 8
#define CURCUMA_ML_MAX_LAMBDA 8
/* Neighbours per atom inside max(rcut_radial, rcut_angular). Exceeding it is an error
 * (CURCUMA_ML_ERR_NEIGHBOURS), never a silently shortened list. */
#define CURCUMA_ML_MAX_NEIGHBOURS 512

/* Return codes of the plugin entry points. */
#define CURCUMA_ML_OK 0
#define CURCUMA_ML_ERR_ARGS 1
#define CURCUMA_ML_ERR_ALLOC 2
#define CURCUMA_ML_ERR_CUDA 3
#define CURCUMA_ML_ERR_BLAS 4
#define CURCUMA_ML_ERR_NOT_READY 5
#define CURCUMA_ML_ERR_NEIGHBOURS 10

/*! \brief The ACSF grid of descriptor.h (ACSFConfig) in a flat, fixed-size form. */
struct CurcumaMLDeviceConfig {
    int n_species;
    int n_radial;
    int n_zeta;
    int n_lambda;
    double rcut_radial;
    double rcut_angular;
    double radial_eta;
    double angular_eta;
    double radial_rs[CURCUMA_ML_MAX_RADIAL];
    double zeta[CURCUMA_ML_MAX_ZETA];
    double zeta_pref[CURCUMA_ML_MAX_ZETA]; /* 2^(1-zeta), computed on the host as in descriptor.cpp */
    int zeta_int[CURCUMA_ML_MAX_ZETA];     /* zeta if it is an integer in [1, 64], else 0: powers by multiplication */
    double lambda[CURCUMA_ML_MAX_LAMBDA];
};
