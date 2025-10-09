/*
 * <Online statistics for trajectory analysis using Welford's algorithm.>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 */

#include "trajectory_statistics.h"

#include <algorithm>
#include <numeric>

TrajectoryStatistics::TrajectoryStatistics(int moving_window_size)
    : m_window_size(moving_window_size)
{
}

void TrajectoryStatistics::addValue(const std::string& name, double value)
{
    // Claude Generated 2025: Welford's online algorithm for numerically stable mean/variance
    // Reference: Donald Knuth, TAOCP Vol 2, section 4.2.2
    //
    // Algorithm: For each new value x:
    //   count = count + 1
    //   delta = x - mean
    //   mean = mean + delta / count
    //   delta2 = x - mean  (using updated mean!)
    //   M2 = M2 + delta * delta2
    //
    // Then: variance = M2 / count, std_dev = sqrt(variance)

    auto& stats = m_stats[name]; // Creates entry if doesn't exist

    stats.count++;
    double delta = value - stats.mean;
    stats.mean += delta / stats.count;
    double delta2 = value - stats.mean; // Use updated mean
    stats.M2 += delta * delta2;

    // Update moving average window
    stats.window.push_back(value);
    if (static_cast<int>(stats.window.size()) > m_window_size) {
        stats.window.pop_front(); // Remove oldest value
    }
}

double TrajectoryStatistics::getMean(const std::string& name) const
{
    auto it = m_stats.find(name);
    if (it == m_stats.end() || it->second.count == 0) {
        return 0.0;
    }
    return it->second.mean;
}

double TrajectoryStatistics::getStdDev(const std::string& name) const
{
    auto it = m_stats.find(name);
    if (it == m_stats.end() || it->second.count < 2) {
        return 0.0; // Need at least 2 values for std dev
    }

    // Sample standard deviation: sqrt(M2 / (n-1))
    // Using n-1 (Bessel's correction) for unbiased estimate
    double variance = it->second.M2 / (it->second.count - 1);
    return std::sqrt(variance);
}

double TrajectoryStatistics::getMovingAverage(const std::string& name) const
{
    auto it = m_stats.find(name);
    if (it == m_stats.end() || it->second.window.empty()) {
        return 0.0;
    }

    // Simple arithmetic mean of values in window
    const auto& window = it->second.window;
    double sum = std::accumulate(window.begin(), window.end(), 0.0);
    return sum / window.size();
}

int TrajectoryStatistics::getCount(const std::string& name) const
{
    auto it = m_stats.find(name);
    if (it == m_stats.end()) {
        return 0;
    }
    return it->second.count;
}

void TrajectoryStatistics::reset(const std::string& name)
{
    auto it = m_stats.find(name);
    if (it != m_stats.end()) {
        it->second.count = 0;
        it->second.mean = 0.0;
        it->second.M2 = 0.0;
        it->second.window.clear();
    }
}

void TrajectoryStatistics::resetAll()
{
    m_stats.clear();
}

void TrajectoryStatistics::setWindowSize(int size)
{
    m_window_size = size;

    // Trim existing windows to new size
    for (auto& [name, stats] : m_stats) {
        while (static_cast<int>(stats.window.size()) > m_window_size) {
            stats.window.pop_front();
        }
    }
}
