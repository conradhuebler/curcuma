/*
 * <Global command-line parameters>
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
 */

#pragma once

#include "src/core/parameter_macros.h"

// Claude Generated (Nov 2026): the names CLI2Json treats as global. Every name in this module
// is copied to controller[name] and controller["global"][name] and is never auto-routed into a
// capability module, so this list is the single place that decides what "global" means.
//
// The defaults here document the value used when the flag is absent. They are NOT merged into
// -export_run: several GPU flags are read with contains() ("given or not"), so a default written
// into an exported file would turn into an explicit user value on replay.
// verbosity, threads, method, charge and spin are also parameters of capability modules, with
// the defaults of that capability; the value given on the command line applies to all of them.
BEGIN_PARAMETER_DEFINITION(global)
MODULE_INFO("Parameters accepted by every command.", "General", {})
    PARAM(verbosity, Int, 1, "Output level: 0=silent, 1=table, 2=detailed, 3=debug.", "General", {})
    PARAM(threads, Int, 1, "Number of CPU threads.", "General", {})
    PARAM(method, String, "gfnff", "Energy method for every capability that computes energies (gfnff, gfn2, uff, ...). In rmsd and confscan, method names the alignment method instead.", "General", {})
    PARAM(charge, Int, 0, "Total charge of the molecule.", "General", {})
    PARAM(spin, Int, 0, "Spin of the molecule, the integer stored by Molecule::setSpin.", "General", {})
    PARAM(export_run, String, "", "Write the resolved configuration of this run to this JSON file; curcuma -import_config file.json replays it.", "General", {})
    PARAM(import_config, String, "", "Read a configuration JSON written by -export_run; command-line values win.", "General", {})
    PARAM(bak, String, "", "Output files to copy from the BMT directory back to the working directory.", "Output", {})
    PARAM(no_bmt, Bool, false, "Write output files to the working directory instead of a Basename.Method.Timestamp directory.", "Output", {})
    PARAM(noprogress, Bool, false, "Disable live progress bars (also off when stdout is not a terminal).", "Output", {})
    PARAM(gpu, String, "none", "GPU backend: none (CPU), cuda, rocm or vulkan. A backend that is not compiled in warns and runs on the CPU.", "GPU", {})
    PARAM(gpu_device, Int, 0, "Device index of the calculation for -gpu.", "GPU", {})
    PARAM(gpu_devices, String, "", "Device pool for batch workers and the device set of one calculation, for example 0,2. Empty: all visible devices.", "GPU", {})
    PARAM(gpu_workers_per_device, Int, 1, "Batch workers sharing one device.", "GPU", {})
    PARAM(gpu_strict, Bool, false, "Stop at the first GPU fallback to the CPU (exit code 3) instead of counting it.", "GPU", {})
    PARAM(gpu_memory_check, Bool, true, "Refuse a basis whose estimated device memory exceeds the free memory (native GFN1/GFN2). false skips the estimate.", "GPU", {})
    PARAM(gpu_sparse_integrals, String, "auto", "Screened pair storage of S, H0 and the multipole integrals: auto, on or off (native GFN1/GFN2).", "GPU", {})
    PARAM(gpu_eigensolver_devices, String, "", "Devices for the distributed eigensolve: all, a list such as 0,1,2,3, or none. Empty: all usable devices when several are usable.", "GPU", {})
    PARAM(gpu_eigensolver_backend, String, "auto", "Distributed eigensolve backend: auto, mp (cuSOLVERMp over NCCL) or mg (cusolverMg).", "GPU", {})
    PARAM(gpu_eigensolver_block, Int, 128, "Column block size of the distributed eigensolve.", "GPU", {})
    PARAM(gpu_eigensolver_min_nao, Int, 4000, "Below this many basis functions the single-GPU eigensolver is used.", "GPU", {})
    PARAM(gpu_eigensolver_fp32, Bool, true, "Also distribute the FP32 iterations of the mixed-precision SCF (backend mp only).", "GPU", {})
    PARAM(gpu_eigensolver_verify, Bool, true, "Check the distributed eigensolve against the single-GPU result.", "GPU", {})
    PARAM(gpu_density_devices, String, "", "Devices over which the pattern density is split. Empty: not split.", "GPU", {})
    PARAM(gpu_density_min_nao, Int, 4000, "Below this many basis functions the density is not split.", "GPU", {})
END_PARAMETER_DEFINITION
