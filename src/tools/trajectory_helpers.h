/*
 * trajectory_helpers.h - Claude Generated 2025
 *
 * JSON schema conversion helper functions for TrajectoryWriter
 * Used by geometry commands (bond, angle, torsion) to standardize output format
 */

#pragma once

#include <vector>
#include <string>
#include "json.hpp"
#include "trajectory_statistics.h"

using json = nlohmann::json;

// ---------- Geometry Command JSON Schema Conversion ----------

/**
 * @brief Convert a vector of values and statistics to TrajectoryWriter-compatible JSON format
 * @param values Vector of measured values across trajectory
 * @param metric_name Name of the metric (e.g., "bond_length", "bond_angle", "dihedral")
 *param stats TrajectoryStatistics object with calculated statistics
 *param additional_data Optional additional metadata to include
 *return JSON object ready for TrajectoryWriter
 */
inline json convertToTrajectoryJSON(
    const std::vector<double>& values,
    const std::string& metric_name,
    const TrajectoryStatistics& stats,
    const json& additional_data = {});

/**
 * @brief Create complete trajectory JSON schema for simple geometric analysis
 * @param values Vector of geometric measurements (bond lengths, angles, dihedrals)
 * @param metric_name Name of the measurement
 *param unit Unit string (e.g., "Å", "°")
 *param stats TrajectoryStatistics with statistical analysis
 *param atom_info Optional atom information for each frame
 *param frame_numbers Frame numbers associated with each value
 *return Complete JSON object with trajectory data and statistics
 */
inline json createTrajectoryJSON(
    const std::vector<double>& values,
    const std::string& metric_name,
    const std::string& unit,
    const TrajectoryStatistics& stats,
    const std::vector<json>& atom_info = {},
    const std::vector<int>& frame_numbers = {});