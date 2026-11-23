/*
 * <GFN-FF device virial slot layout>
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
 * Claude Generated (Nov 2026, docs/GPU_PERIODIC_GAPS.md A1): layout of the device virial buffer, shared by the
 * CUDA kernels (gfnff_kernels.cu) and the host wrapper (gfnff_gpu_method_impl.h). No runtime-API types.
 */

#pragma once

/// Nine row-major values Xi_ab (Eh; a = position, b = gradient component) per term. VIR_CN is the whole CN chain
/// rule (bond dr0/dCN, dC6/dCN incl. the D4 tail, TERM 1b): the device keeps one dE/dCN vector.
enum GfnffVirialSlot { VIR_BOND = 0, VIR_ANGLE, VIR_DIHED, VIR_INV, VIR_STORS, VIR_DISP, VIR_REP, VIR_HBOND, VIR_XBOND,
                       VIR_BATM, VIR_CN, VIR_N_TENSORS };
/// One energy value after the tensors: the C8 part of the homogeneous D4 tail (pressure correction E6/E8)
constexpr int GFNFF_VIR_TAIL_E8 = 9 * VIR_N_TENSORS;
/// Doubles in the buffer
constexpr int GFNFF_VIR_BUFFER = 9 * VIR_N_TENSORS + 1;
