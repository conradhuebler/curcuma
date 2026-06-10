/*
 * test_bmt_utils.cpp - Unit tests for BMT Output Directory System
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Tests: createBMTDir, outputPath, writeMetadata, processBakFiles
 *        (existing source, missing source, empty list)
 *
 * Claude Generated 2026
 */

#include "src/tools/bmt_utils.h"
#include "src/core/curcuma_logger.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define TEST_ASSERT(condition, msg) \
    do { \
        tests_run++; \
        if (condition) { \
            std::cout << "  [PASS] " << msg << std::endl; \
            tests_passed++; \
        } else { \
            std::cout << "  [FAIL] " << msg << std::endl; \
            tests_failed++; \
        } \
    } while (0)

// Helper: remove directory recursively
static void cleanup_dir(const std::string& path)
{
#ifdef C17
#ifndef _WIN32
    std::filesystem::remove_all(path);
#endif
#endif
}

void test_stripExtension()
{
    std::cout << "\n=== test_stripExtension ===" << std::endl;

    // Basic .xyz stripping
    TEST_ASSERT(BMTUtils::stripExtension("input.xyz") == "input",
                "stripExtension('input.xyz') == 'input'");

    // .mol2 extension
    TEST_ASSERT(BMTUtils::stripExtension("molecule.mol2") == "molecule",
                "stripExtension('molecule.mol2') == 'molecule'");

    // .sdf extension
    TEST_ASSERT(BMTUtils::stripExtension("complex.sdf") == "complex",
                "stripExtension('complex.sdf') == 'complex'");

    // .pdb extension
    TEST_ASSERT(BMTUtils::stripExtension("protein.pdb") == "protein",
                "stripExtension('protein.pdb') == 'protein'");

    // No extension
    TEST_ASSERT(BMTUtils::stripExtension("noext") == "noext",
                "stripExtension('noext') == 'noext'");

    // Multi-dot filename (only last extension stripped)
    // Note: std::filesystem::path::stem() strips only the last extension
    // "input.opt.xyz" -> "input.opt" (stem strips .xyz, not all extensions)
    std::string multi_dot_result = BMTUtils::stripExtension("input.opt.xyz");
    TEST_ASSERT(multi_dot_result == "input.opt",
                "stripExtension('input.opt.xyz') == 'input.opt'");
}

void test_outputPath()
{
    std::cout << "\n=== test_outputPath ===" << std::endl;

    // Empty BMT dir: passthrough
    TEST_ASSERT(BMTUtils::outputPath("", "file.xyz") == "file.xyz",
                "outputPath('', 'file.xyz') == 'file.xyz'");

    // With BMT dir: concatenation
    TEST_ASSERT(BMTUtils::outputPath("input.opt.20260609", "file.xyz") == "input.opt.20260609/file.xyz",
                "outputPath('input.opt.20260609', 'file.xyz') == 'input.opt.20260609/file.xyz'");
}

void test_createBMTDir()
{
    std::cout << "\n=== test_createBMTDir ===" << std::endl;

    // Create a BMT directory
    std::string bmt_dir = BMTUtils::createBMTDir("test_mol", "opt");

    // Verify directory name pattern: basename.keyword.YYYYMMDD_HHMMSS
    TEST_ASSERT(bmt_dir.find("test_mol.opt.") == 0,
                "BMT dir starts with 'test_mol.opt.'");

    // Verify directory was created
#ifdef C17
#ifndef _WIN32
    TEST_ASSERT(std::filesystem::exists(bmt_dir),
                "BMT directory exists on disk");
#endif
#endif

    // Cleanup
    cleanup_dir(bmt_dir);
}

void test_writeMetadata()
{
    std::cout << "\n=== test_writeMetadata ===" << std::endl;

    std::string bmt_dir = BMTUtils::createBMTDir("test_meta", "md");
    BMTUtils::writeMetadata(bmt_dir, "test_meta", "md", "test_meta.xyz");

    // Read metadata file and verify contents
    std::string meta_path = bmt_dir + "/metadata.txt";
    std::ifstream meta(meta_path);
    TEST_ASSERT(meta.is_open(), "metadata.txt file exists and is readable");

    if (meta.is_open()) {
        std::string line;
        bool found_basename = false, found_method = false, found_input = false;
        while (std::getline(meta, line)) {
            if (line.find("basename: test_meta") != std::string::npos) found_basename = true;
            if (line.find("method: md") != std::string::npos) found_method = true;
            if (line.find("input_file: test_meta.xyz") != std::string::npos) found_input = true;
        }
        TEST_ASSERT(found_basename, "metadata.txt contains basename");
        TEST_ASSERT(found_method, "metadata.txt contains method");
        TEST_ASSERT(found_input, "metadata.txt contains input_file");
    }

    // Cleanup
    meta.close();
    cleanup_dir(bmt_dir);
}

void test_processBakFiles_existing_source()
{
    std::cout << "\n=== test_processBakFiles (existing source) ===" << std::endl;

    std::string bmt_dir = BMTUtils::createBMTDir("test_bak", "analysis");

    // Create a file inside BMT dir
    std::string src_path = BMTUtils::outputPath(bmt_dir, "result.dat");
    {
        std::ofstream out(src_path);
        out << "test data for bak" << std::endl;
    }

    // Process bak files: should copy result.dat from BMT dir to CWD
    std::vector<std::string> bak_files = {"result.dat"};
    BMTUtils::processBakFiles(bmt_dir, bak_files);

    // Verify file was copied to CWD
    TEST_ASSERT(std::ifstream("result.dat").is_open(),
                "processBakFiles copies existing file to CWD");

    // Cleanup
    std::remove("result.dat");
    cleanup_dir(bmt_dir);
}

void test_processBakFiles_missing_source()
{
    std::cout << "\n=== test_processBakFiles (missing source) ===" << std::endl;

    std::string bmt_dir = BMTUtils::createBMTDir("test_bak_missing", "opt");

    // Request a file that doesn't exist in BMT dir
    std::vector<std::string> bak_files = {"nonexistent_file.xyz"};
    // Should not crash, just warn
    BMTUtils::processBakFiles(bmt_dir, bak_files);

    TEST_ASSERT(true, "processBakFiles with missing source does not crash");

    // Cleanup
    cleanup_dir(bmt_dir);
}

void test_processBakFiles_empty_list()
{
    std::cout << "\n=== test_processBakFiles (empty list) ===" << std::endl;

    std::string bmt_dir = BMTUtils::createBMTDir("test_bak_empty", "md");

    // Empty bak_files list: should be a no-op
    std::vector<std::string> bak_files;
    BMTUtils::processBakFiles(bmt_dir, bak_files);

    TEST_ASSERT(true, "processBakFiles with empty list does not crash");

    // Cleanup
    cleanup_dir(bmt_dir);
}

void test_collectBakFiles()
{
    std::cout << "\n=== test_collectBakFiles ===" << std::endl;

    // Test with no "bak" key
    nlohmann::json ctrl1 = nlohmann::json::object();
    auto result1 = BMTUtils::collectBakFiles(ctrl1);
    TEST_ASSERT(result1.empty(), "collectBakFiles with no 'bak' key returns empty");

    // Test with string "bak"
    nlohmann::json ctrl2 = nlohmann::json::object();
    ctrl2["bak"] = "result.xyz";
    auto result2 = BMTUtils::collectBakFiles(ctrl2);
    TEST_ASSERT(result2.size() == 1 && result2[0] == "result.xyz",
                "collectBakFiles with string 'bak' returns single-element vector");

    // Test with array "bak"
    nlohmann::json ctrl3 = nlohmann::json::object();
    ctrl3["bak"] = nlohmann::json::array({"file1.xyz", "file2.dat"});
    auto result3 = BMTUtils::collectBakFiles(ctrl3);
    TEST_ASSERT(result3.size() == 2 && result3[0] == "file1.xyz" && result3[1] == "file2.dat",
                "collectBakFiles with array 'bak' returns correct vector");
}

int main()
{
    std::cout << "=== BMT Utils Unit Tests ===" << std::endl;

    test_stripExtension();
    test_outputPath();
    test_createBMTDir();
    test_writeMetadata();
    test_processBakFiles_existing_source();
    test_processBakFiles_missing_source();
    test_processBakFiles_empty_list();
    test_collectBakFiles();

    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Tests run:    " << tests_run << std::endl;
    std::cout << "Tests passed:  " << tests_passed << std::endl;
    std::cout << "Tests failed:  " << tests_failed << std::endl;

    return tests_failed > 0 ? 1 : 0;
}