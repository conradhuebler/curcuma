/*
 * <BMT Output Directory Utility>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Standalone utility for creating and managing BMT (Basename.Method.Timestamp)
 * output directories. Works with both CurcumaMethod-based and standalone
 * command handlers.
 *
 * Claude Generated 2026
 */

#pragma once

#include <string>
#include <vector>

#include "src/core/global.h" // for nlohmann::json

namespace BMTUtils {

// Create BMT directory: basename.keyword.YYYYMMDD_HHMMSS
// Returns the created directory path
std::string createBMTDir(const std::string& basename, const std::string& keyword);

// Write metadata.txt inside BMT directory with calculation info
void writeMetadata(const std::string& bmt_dir,
                   const std::string& basename,
                   const std::string& method,
                   const std::string& input_file);

// Copy registered files from BMT directory back to working directory
void processBakFiles(const std::string& bmt_dir,
                     const std::vector<std::string>& bak_files);

// Build full output path: bmt_dir/filename (or just filename if bmt_dir empty)
std::string outputPath(const std::string& bmt_dir, const std::string& filename);

// Strip file extension from filename (handles .xyz, .mol2, .sdf, .pdb, etc.)
// Returns the basename without extension
std::string stripExtension(const std::string& filename);

// Collect -bak file names from JSON controller
// Handles both string and array forms of the "bak" parameter
std::vector<std::string> collectBakFiles(const nlohmann::json& controller);

} // namespace BMTUtils