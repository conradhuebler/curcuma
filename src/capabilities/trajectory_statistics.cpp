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
#include <stdexcept>

// -----------------
// Existing methods (unchanged)
// -----------------

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

    // Update min/max tracking
    stats.min_val = std::min(stats.min_val, value);
    stats.max_val = std::max(stats.max_val, value);

    // Update moving average window
    stats.window.push_back(value);
    if (static_cast<int>(stats.window.size()) > m_window_size) {
        stats.window.pop_front(); // Remove oldest value
    }

    // Optionally store full series data
    if (m_store_full_series) {
        // Check if this metric should be stored (if filter is active)
        bool should_store = m_full_series_metrics.empty() ||
                           std::find(m_full_series_metrics.begin(), m_full_series_metrics.end(), name) != m_full_series_metrics.end();

        if (should_store) {
            stats.full_series.push_back(value);
        }
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

// -----------------
// Extended methods for TrajectoryWriter integration
// -----------------

double TrajectoryStatistics::getMin(const std::string& name) const
{
    auto it = m_stats.find(name);
    if (it == m_stats.end() || it->second.count == 0) {
        return 0.0;
    }
    return it->second.min_val;
}

double TrajectoryStatistics::getMax(const std::string& name) const
{
    auto it = m_stats.find(name);
    if (it == m_stats.end() || it->second.count == 0) {
        return 0.0;
    }
    return it->second.max_val;
}

double TrajectoryStatistics::getMedian(const std::string& name) const
{
    auto it = m_stats.find(name);
    if (it == m_stats.end() || it->second.full_series.empty()) {
        return 0.0;
    }

    return calculateMedian(it->second.full_series);
}

double TrajectoryStatistics::getRange(const std::string& name) const
{
    return getMax(name) - getMin(name);
}

double TrajectoryStatistics::getVariance(const std::string& name) const
{
    auto it = m_stats.find(name);
    if (it == m_stats.end() || it->second.count == 0) {
        return 0.0;
    }

    // Population variance (not sample variance)
    return it->second.M2 / it->second.count;
}

json TrajectoryStatistics::exportStatistics(const std::string& name) const
{
    auto it = m_stats.find(name);
    if (it == m_stats.end()) {
        return json{};
    }

    const auto& stats = it->second;
    json result;

    result["count"] = stats.count;
    result["mean"] = stats.mean;
    result["std"] = getStdDev(name);
    result["variance"] = getVariance(name);
    result["min"] = stats.min_val;
    result["max"] = stats.max_val;
    result["range"] = stats.max_val - stats.min_val;

    if (!stats.full_series.empty()) {
        result["median"] = calculateMedian(stats.full_series);
    } else {
        result["median"] = 0.0;  // Not available without full series
    }

    result["moving_avg"] = getMovingAverage(name);

    return result;
}

json TrajectoryStatistics::exportAllStatistics() const
{
    json result;
    for (const auto& [name, stats] : m_stats) {
        result[name] = exportStatistics(name);
    }
    return result;
}

std::vector<double> TrajectoryStatistics::getSeries(const std::string& name) const
{
    auto it = m_stats.find(name);
    if (it == m_stats.end()) {
        return {};
    }
    return it->second.full_series;
}

void TrajectoryStatistics::setStoreFullSeries(bool store, const std::vector<std::string>& metrics)
{
    m_store_full_series = store;
    m_full_series_metrics = metrics;

    if (!store) {
        // Clear any existing full series data
        for (auto& [name, stats] : m_stats) {
            stats.full_series.clear();
            stats.full_series.shrink_to_fit();
        }
    }
}

bool TrajectoryStatistics::isStoringFullSeries() const
{
    return m_store_full_series;
}

size_t TrajectoryStatistics::getSeriesLength(const std::string& name) const
{
    auto it = m_stats.find(name);
    if (it == m_stats.end()) {
        return 0;
    }
    return it->second.full_series.size();
}

double TrajectoryStatistics::calculateMedian(const std::vector<double>& data)
{
    if (data.empty()) {
        return 0.0;
    }

    std::vector<double> sorted_data = data;  // Copy to avoid modifying original
    std::sort(sorted_data.begin(), sorted_data.end());

    size_t n = sorted_data.size();
    if (n % 2 == 1) {
        return sorted_data[n / 2];  // Middle element for odd length
    } else {
        // Average of two middle elements for even length
        return (sorted_data[n / 2 - 1] + sorted_data[n / 2]) / 2.0;
    }
}
