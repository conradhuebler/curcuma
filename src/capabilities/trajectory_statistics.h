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

#pragma once

#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <string>
#include <vector>
#include "json.hpp"

using json = nlohmann::json;

/*! \brief Online statistics calculator for trajectory time-series data - Claude Generated 2025
 *
 * Implements Welford's numerically stable online algorithm for computing
 * cumulative mean and variance, plus moving average with configurable window.
 *
 * Key features:
 * - O(1) memory per metric (except moving average window)
 * - Numerically stable for large datasets
 * - Multiple independent metrics tracked simultaneously
 * - Thread-safe for single-writer scenarios
 *
 * Reference: Donald Knuth, Art of Computer Programming Vol 2, section 4.2.2
 */
class TrajectoryStatistics {
public:
    TrajectoryStatistics(int moving_window_size = 10);
    ~TrajectoryStatistics() = default;

    /*! \brief Add new value for a named metric
     * \param name Metric identifier (e.g., "gyration_unweighted", "rout")
     * \param value New data point to add to statistics
     *
     * Updates cumulative mean, variance, and moving average window
     */
    void addValue(const std::string& name, double value);

    /*! \brief Get cumulative mean (from all timesteps)
     * \param name Metric identifier
     * \return Mean of all values seen so far
     */
    double getMean(const std::string& name) const;

    /*! \brief Get cumulative standard deviation
     * \param name Metric identifier
     * \return Sample standard deviation (sqrt(variance / (n-1)))
     */
    double getStdDev(const std::string& name) const;

    /*! \brief Get moving average over last N timesteps
     * \param name Metric identifier
     * \return Mean of values in current window, or cumulative mean if window not full
     */
    double getMovingAverage(const std::string& name) const;

    /*! \brief Get number of data points seen for a metric
     * \param name Metric identifier
     * \return Count of addValue() calls for this metric
     */
    int getCount(const std::string& name) const;

    /*! \brief Reset all statistics for a metric
     * \param name Metric identifier
     */
    void reset(const std::string& name);

    /*! \brief Reset all statistics for all metrics */
    void resetAll();

    /*! \brief Set moving window size (affects future calculations)
     * \param size Window size for moving average
     */
    void setWindowSize(int size);

    // -----------------
    // Extended methods for TrajectoryWriter integration
    // -----------------

    /*! \brief Get minimum value seen for a metric
     * \param name Metric identifier
     * \return Minimum value observed so far
     */
    double getMin(const std::string& name) const;

    /*! \brief Get maximum value seen for a metric
     * \param name Metric identifier
     * \return Maximum value observed so far
     */
    double getMax(const std::string& name) const;

    /*! \brief Get median value (requires full series storage)
     * \param name Metric identifier
     * \return Median of all values seen so far
     */
    double getMedian(const std::string& name) const;

    /*! \brief Get range of values for a metric
     * \param name Metric identifier
     * \return Difference between maximum and minimum
     */
    double getRange(const std::string& name) const;

    /*! \brief Get variance for a metric
     * \param name Metric identifier
     * \return Population variance (M2 / count)
     */
    double getVariance(const std::string& name) const;

    /*! \brief Export complete statistics for a metric as JSON
     * \param name Metric identifier
     * \return JSON with all statistics (mean, std, min, max, median, etc.)
     */
    json exportStatistics(const std::string& name) const;

    /*! \brief Export all statistics for all metrics as JSON
     * \return JSON with complete statistics for all metrics
     */
    json exportAllStatistics() const;

    /*! \brief Get full time series data for a metric
     * \param name Metric identifier
     * \return Vector of all values seen (requires full series storage)
     */
    std::vector<double> getSeries(const std::string& name) const;

    /*! \brief Set whether to store full time series data
     * \param store Enable/disable full series storage
     * \param metrics Optional list of metrics to store (empty = all)
     */
    void setStoreFullSeries(bool store, const std::vector<std::string>& metrics = {});

    /*! \brief Check if full series storage is enabled
     * \return True if storing complete time series
     */
    bool isStoringFullSeries() const;

    /*! \brief Get length of stored time series for a metric
     * \param name Metric identifier
     * \return Number of values stored
     */
    size_t getSeriesLength(const std::string& name) const;

private:
    /*! \brief Statistics for a single metric - Welford's algorithm state */
    struct MetricStats {
        int count = 0; ///< Number of values seen
        double mean = 0.0; ///< Cumulative mean
        double M2 = 0.0; ///< Sum of squared differences from mean
        double min_val = std::numeric_limits<double>::max(); ///< Minimum value observed
        double max_val = -std::numeric_limits<double>::max(); ///< Maximum value observed
        std::deque<double> window; ///< Circular buffer for moving average
        std::vector<double> full_series; ///< Complete time series (optional, memory-intensive)
    };

    std::map<std::string, MetricStats> m_stats; ///< Per-metric statistics
    int m_window_size; ///< Moving average window size
    bool m_store_full_series = false; ///< Enable full series storage
    std::vector<std::string> m_full_series_metrics; ///< Which metrics to store fully

    // Helper methods
    static double calculateMedian(const std::vector<double>& data);
};
