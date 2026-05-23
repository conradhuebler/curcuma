/*
 * < MD Diagnostics JSONL Writer (WP-S2, May 2026) >
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#include "md_diagnostics.h"

#include <vector>

MDDiagnosticsWriter::MDDiagnosticsWriter(const std::string& path)
{
    m_out.open(path, std::ios::out | std::ios::app);
}

void MDDiagnosticsWriter::writeSnapshot(int step, double time_fs,
                                        const json& energy_decomp,
                                        const Vector& charges,
                                        const Vector& cn,
                                        const Matrix& gradient,
                                        int hb_count, int xb_count,
                                        const json& timing)
{
    if (!m_out.is_open()) {
        return;
    }

    json rec;
    rec["step"] = step;
    rec["time_fs"] = time_fs;
    rec["energy"] = energy_decomp;

    rec["charges"] = std::vector<double>(
        charges.data(), charges.data() + charges.size());
    rec["cn"] = std::vector<double>(
        cn.data(), cn.data() + cn.size());

    std::vector<double> gnorm(gradient.rows());
    for (Eigen::Index i = 0; i < gradient.rows(); ++i) {
        gnorm[i] = gradient.row(i).norm();
    }
    rec["gradient_norm"] = gnorm;

    rec["hb_count"] = hb_count;
    rec["xb_count"] = xb_count;

    // WP-P1 (May 2026): optional per-phase wall-clock breakdown
    if (!timing.empty()) {
        rec["timing_ms"] = timing;
    }

    m_out << rec.dump() << "\n";
    m_out.flush();
}
