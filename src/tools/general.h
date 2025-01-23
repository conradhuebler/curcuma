/*
 * <Geometry tools for chemical structures.>
 * Copyright (C) 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <Eigen/Dense>

#include <chrono>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <fmt/color.h>
#include <fmt/core.h>

#include "src/core/elements.h"
#include "src/core/global.h"
#include "src/core/molecule.h"
using namespace curcuma;

#include "src/tools/formats.h"

class RunTimer {
public:
    RunTimer(bool print = false)
        : m_print(print)
    {
        m_start = std::chrono::system_clock::now();
        if (m_print) {
            std::time_t end_time = std::chrono::system_clock::to_time_t(m_start);
            std::cout << "Started computation at " << std::ctime(&end_time) << std::endl;
        }
    }

    ~RunTimer()
    {
        if (m_print) {
            std::cout << "\n\nFinished after " << Elapsed() / 1000 << " seconds!" << std::endl;
            std::time_t end_time = std::chrono::system_clock::to_time_t(m_end);
            std::cout << "Finished computation at " << std::ctime(&end_time) << std::endl;
        }
    }

    inline int Elapsed()
    {
        m_end = std::chrono::system_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(m_end - m_start).count();
    }

    inline void Reset()
    {
        m_start = std::chrono::system_clock::now();
    }

private:
    std::chrono::time_point<std::chrono::system_clock> m_start, m_end;
    bool m_print;
};



namespace Tools {

inline StringList SplitString(const std::string& string)
{
    StringList elements;
    std::string element;
    const char* delim = " ";
    for (const char& c : string) {
        if (*delim != c)
            element.push_back(c);
        else {
            if (element.size()) {
                elements.push_back(element);
            }
            element.clear();
        }
    }
    elements.push_back(element);
    return elements;
}

inline StringList SplitString(const std::string& string, const char* delim)
{
    StringList elements;
    std::string element;
    for (const char& c : string) {
        if (*delim != c)
            element.push_back(c);
        else {
            if (element.size()) {
                elements.push_back(element);
            }
            element.clear();
        }
    }
    elements.push_back(element);
    return elements;
}

inline bool isInt(const std::string& input)
{
    return std::all_of(input.begin(), input.end(), ::isdigit);
}

inline bool isDouble(const std::string& input)
{
    bool isD = true;
    try {
        std::stod(input);
    } catch (const std::string& what_arg) {
        isD = false;
    } catch (const std::invalid_argument& argument) {
        isD = false;
    }
    return isD;
    /*
    const char* delim = ".";
    StringList list = SplitString(input, delim);
    if (list.size() != 2)
        return false;

    bool left = isInt(list[0]);
    bool right = isInt(list[1]);
    return left && right;
    */
    // return std::all_of(input.begin(), input.end(), ::isdigit);
}

inline bool StringContains(const std::string& string, const std::string& search)
{
    return string.find(search) != std::string::npos;
}

inline int VectorDifference(const std::vector<int>& tmp_a, const std::vector<int>& tmp_b)
{
    int difference = 0;
    std::vector<int> a, b;
    if (tmp_a.size() < tmp_b.size()) {
        a = tmp_a;
        b = tmp_b;
    } else {
        b = tmp_a;
        a = tmp_b;
    }

    difference += b.size() - a.size();
    for (int i : a) {
        std::vector<int>::iterator it = std::find(b.begin(), b.end(), i);
        difference += it == b.end();
    }
    return difference;
}

inline double median(const std::vector<double>& vector)
{
    if (vector.size() == 0)
        return 0;
    double median = 0;

    if (vector.size() % 2)
        median = vector[vector.size() / 2];
    else {
        double left = vector[vector.size() / 2 - 1];
        double right = vector[vector.size() / 2];
        median = (left + right) * 0.5;
    }
    return median;
}

inline double mean(const std::vector<double>& vector)
{
    double mean = 0;
    for (double val : vector)
        mean += val;
    return mean / double(vector.size());
}

inline double stdev(const std::vector<double>& vector, double mean)
{
    double stdev = 0;
    for (double val : vector)
        stdev += (val - mean) * (val - mean);
    return sqrt(stdev / double(vector.size() - 1));
}

inline std::vector<std::pair<double, double>> Histogram(const std::vector<double>& vector, int bins)
{
    std::vector<std::pair<double, double>> histogram;
    if (vector.size() == 0)
        return histogram;

    double min = vector[0];
    double max = vector[0];

    for (double val : vector) {
        min = std::min(min, val);
        max = std::max(max, val);
    }

    std::vector<std::pair<double, int>> hist;
    double h = (max - min) / double(bins);

    std::vector<double> x(bins, 0);

    for (int j = 0; j < bins; j++) {
        x[j] = min + h / 2. + j * h;
        std::pair<double, int> bin;
        bin.second = 0;
        bin.first = min + h / 2. + j * h;
        hist.push_back(bin);
    }

    for (double val : vector) {
        int jStar = std::floor((val - min) / h); // use the floor function to get j* for the nearest point to x_j* to phi
        if (jStar >= bins || jStar < 0)
            continue; // if you are outside the grid continue
        hist[jStar].second++;
    }
    int max_val = 0;
    for (auto element : hist)
        max_val = std::max(max_val, element.second);

    for (auto element : hist)
        histogram.push_back(std::pair<double, double>(element.first, element.second / double(max_val)));

    return histogram;
}

inline double ShannonEntropy(const std::vector<std::pair<double, double>>& histogram)
{
    if (histogram.size() == 0)
        return 0;
    double entropy = 0;
    double sum = 0.0;

    for (int i = 0; i < histogram.size(); ++i) {
        sum += histogram[i].second;
    }

    double d = (histogram[histogram.size() - 1].first - histogram[0].first) / double(histogram.size());
    double lower = 1 / double(histogram.size());

    for (int i = 0; i < histogram.size(); ++i) {
        if (histogram[i].second < lower)
            continue;
        const double value = histogram[i].second / sum * d;
        entropy += value * log2(value);
    }

    return -1 * entropy;
}

inline std::string Vector2String(const std::vector<int>& vector, const std::string& delim = "|")
{
    std::string result = " ";

    for (auto i : vector)
        result += std::to_string(i) + delim;
    result.pop_back();

    result += " ";

    return result;
}

inline std::string DoubleVector2String(const std::vector<double>& vector, const std::string& delim = "|")
{
    std::string result = "";

    for (auto i : vector)
        result += std::to_string(i) + delim;
    result.pop_back();

    result += " ";

    return result;
}

inline std::string DoubleVector2String(const Vector& vector, const std::string& delim = "|")
{
    std::string result = "";

    for (int i = 0; i < vector.rows(); ++i)
        result += std::to_string(vector[i]) + delim;
    result.pop_back();

    result += " ";

    return result;
}

inline std::string Matrix2String(const Matrix& matrix)
{
    std::string result = "";

    for (int i = 0; i < matrix.cols(); ++i) {
        for (int j = 0; j < matrix.rows(); ++j)
            result += std::to_string(matrix(i, j)) + "|";
        result += "|";
    }
    result.pop_back();
    result.pop_back();

    result += "";

    return result;
}

inline std::vector<int> String2Vector(const std::string& string)
{
    std::vector<int> vector;
    StringList tmp = SplitString(string, "|");
    for (auto i : tmp) {
        try {
            vector.push_back(std::stoi(i));
        } catch (const std::invalid_argument& argument) {
            continue;
        }
    }
    return vector;
}

inline Vector String2EigenVector(const std::string& string, const char* delim = "|")
{
    std::vector<double> tmpvector;

    StringList tmp = SplitString(string, delim);
    for (auto i : tmp) {
        try {
            tmpvector.push_back(std::stod(i));
        } catch (const std::invalid_argument& argument) {
            continue;
        }
    }
    return Eigen::Map<Vector>(tmpvector.data(), tmpvector.size());
    ;
}

inline std::string VectorVector2String(const std::vector<std::vector<int>>& vector)
{
    std::string result;
    if (vector.size() == 0)
        return result;

    for (auto vec : vector) {
        result += Vector2String(vec) + ";";
    }
    result.pop_back();
    return result;
}

inline std::vector<std::vector<int>> String2VectorVector(const std::string& string)
{
    std::vector<std::vector<int>> result;

    StringList vectors = SplitString(string, ";");
    for (const auto& element : vectors)
        result.push_back(String2Vector(element));
    return result;
}

inline std::vector<int> RandomVector(int start, int end)
{
    /* fast taken from
     * https://www.geeksforgeeks.org/generating-random-number-range-c/ */
    std::vector<int> sequence;
    while(sequence.size() < (end - start))
    {
        int num = (rand() % (end - start)) + start;
        if(std::find(sequence.begin(), sequence.end(), num) == sequence.end())
            sequence.push_back(num);
    }
    return sequence;
}

inline std::string DoublePtr2String(const double* array, int size)
{
    std::string result;
    if (size == 0)
        return result;

    for (int i = 0; i < size; ++i) {
        result += std::to_string(array[i]) + ";";
    }
    result.pop_back();
    return result;
}

inline std::vector<double> String2DoubleVec(const std::string& str, const char* delim = ";")
{
    std::vector<double> result;

    StringList vectors = SplitString(str, delim);
    for (const auto& element : vectors)
        result.push_back(std::stod(element));

    return result;
}

inline std::vector<int> CreateList(const std::string& list)
{
    std::vector<int> result;

    StringList tmp_list = Tools::SplitString(list, ",");
    for (const std::string& single : tmp_list) {
        if (Tools::StringContains(single, ":")) {
            auto sub = Tools::SplitString(single, ":");
            if (sub.size() != 2)
                continue;
            if (!Tools::isInt(sub[0]) || !Tools::isInt(sub[1]))
                continue;
            int start = std::stoi(sub[0]);
            int ende = std::stoi(sub[1]);

            while (start <= ende) {
                result.push_back(start);
                start++;
            }
        } else {
            if (Tools::isInt(single)) {
                int i = std::stoi(single);
                result.push_back(i);
            }
        }
    }
    return result;
}
}
