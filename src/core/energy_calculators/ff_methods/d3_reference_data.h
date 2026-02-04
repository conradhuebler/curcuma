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

#include <vector>

// External declarations for D3 reference data
// These are defined in d3_reference_cn.cpp and d3_reference_c6.cpp

/// Coordination number reference data (721 values)
/// Structure: 103 elements × 7 reference states
/// Data source: s-dftd3 reference.f90
extern std::vector<double> reference_cn_data_complete;

/// C6 dispersion coefficient data (262,444 values)
/// Structure: 5,356 element pairs × 7×7 reference combinations
/// Data source: s-dftd3 reference.f90
extern std::vector<double> reference_c6_data_complete;
