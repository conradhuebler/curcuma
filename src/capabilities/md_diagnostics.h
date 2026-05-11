/*
 * < MD Diagnostics JSONL Writer (WP-S2, May 2026) >
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#pragma once

#include "src/core/global.h"
#include "json.hpp"

#include <fstream>
#include <string>

using json = nlohmann::json;

/**
 * @brief Streaming writer for per-step MD diagnostics in JSONL format.
 *
 * Claude Generated (WP-S2, May 2026): writes one JSON object per line to
 * `<basename>.diag.jsonl`. Used by SimpleMD to record energy decomposition,
 * charges, CN, gradient norms, HB/XB counts at every dump-frequency step.
 *
 * File is opened in append mode and flushed after each snapshot for
 * crash-robust streaming.
 */
class MDDiagnosticsWriter {
public:
    explicit MDDiagnosticsWriter(const std::string& path);
    ~MDDiagnosticsWriter() = default;

    bool isOpen() const { return m_out.is_open(); }

    /**
     * @brief Append one snapshot record to the JSONL file.
     * Per-atom Vectors with size 0 produce empty JSON arrays (not an error).
     */
    void writeSnapshot(int step, double time_fs,
                       const json& energy_decomp,
                       const Vector& charges,
                       const Vector& cn,
                       const Matrix& gradient,
                       int hb_count, int xb_count);

private:
    std::ofstream m_out;
};
