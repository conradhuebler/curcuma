/*
 * D3 Reference Data Declarations
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * External declarations for D3 reference data from s-dftd3
 * Data is split across multiple .cpp files for faster compilation:
 * - d3_reference_cn.cpp: 721 coordination number reference values
 * - d3_reference_c6.cpp: 262,444 C6 dispersion coefficients
 *
 * Claude Generated - December 20, 2025
 */

#pragma once

#include <array>

// External declarations for D3 reference data
// These are defined in d3_reference_cn.cpp and d3_reference_c6.cpp
//
// Claude Generated (Jul 2026): these are `const std::array`, NOT `std::vector`. A global
// std::vector runs a constructor at static-init time (heap-allocates + copies the whole
// table on EVERY process start, ~2 MB for the C6 table — pure startup tax even when D3 is
// never used). A `const std::array` of literals is *constant-initialized* into .rodata:
// zero startup cost, mapped read-only. Only .size()/operator[] are used, so the switch is
// transparent to callers.

/// Coordination number reference data (824 values, >= 721 used)
/// Structure: 103 elements × 7 reference states
/// Data source: s-dftd3 reference.f90
extern const std::array<double, 824> reference_cn_data_complete;

/// C6 dispersion coefficient data (262,444 values)
/// Structure: 5,356 element pairs × 7×7 reference combinations
/// Data source: s-dftd3 reference.f90
extern const std::array<double, 262444> reference_c6_data_complete;
