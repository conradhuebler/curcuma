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

#include "src/tools/bmt_utils.h"
#include "src/core/curcuma_logger.h"

#include <chrono>
#include <ctime>
#include <fstream>
#include <string>
#include <vector>

#ifdef C17
#ifndef _WIN32
#include <filesystem>
namespace fs = std::filesystem;
#endif
#endif

namespace BMTUtils {

std::string createBMTDir(const std::string& basename, const std::string& keyword)
{
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    std::tm* local_tm = std::localtime(&now_time);
    char timestamp[32];
    std::strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", local_tm);

    std::string bmt_dir = basename + "." + keyword + "." + timestamp;

#ifdef C17
#ifndef _WIN32
    std::filesystem::create_directories(bmt_dir);
#endif
#endif

    CurcumaLogger::info_fmt("Output directory: {}", bmt_dir);
    return bmt_dir;
}

void writeMetadata(const std::string& bmt_dir,
                   const std::string& basename,
                   const std::string& method,
                   const std::string& input_file)
{
    std::string meta_path = bmt_dir + "/metadata.txt";
    std::ofstream meta(meta_path);
    if (!meta.is_open()) {
        CurcumaLogger::warn_fmt("Could not create metadata file: {}", meta_path);
        return;
    }

    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    std::tm* local_tm = std::localtime(&now_time);
    char timestamp[32];
    std::strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", local_tm);

    meta << "# Curcuma Calculation Metadata\n";
    meta << "basename: " << basename << "\n";
    meta << "method: " << method << "\n";
    meta << "timestamp: " << timestamp << "\n";
    meta << "input_file: " << input_file << "\n";
    meta.close();
}

void processBakFiles(const std::string& bmt_dir,
                     const std::vector<std::string>& bak_files)
{
    if (bmt_dir.empty() || bak_files.empty())
        return;

#ifdef C17
#ifndef _WIN32
    for (const auto& bak_file : bak_files) {
        std::string src = bmt_dir + "/" + bak_file;
        if (std::filesystem::exists(src)) {
            std::filesystem::copy_file(src, bak_file,
                std::filesystem::copy_options::overwrite_existing);
            CurcumaLogger::info_fmt("Copied {} to working directory", bak_file);
        } else {
            CurcumaLogger::warn_fmt("Backup file {} not found in output directory", bak_file);
        }
    }
#endif
#endif
}

std::string outputPath(const std::string& bmt_dir, const std::string& filename)
{
    if (bmt_dir.empty())
        return filename;
    return bmt_dir + "/" + filename;
}

} // namespace BMTUtils