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

#include "src/core/fileiterator.h"
#include "src/core/global.h"
#include "src/core/molecule.h"

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

private:
    std::chrono::time_point<std::chrono::system_clock> m_start, m_end;
    bool m_print;
};

namespace General {

inline void StartUp(int argc, char** argv)
{
    std::cout << "*************************************************" << std::endl
              << "*   Curcuma - Simple Molecular Modelling tool   *" << std::endl
              << "*                                               *" << std::endl
              << "*    Curcuma - Version " << git_tag << "                   *" << std::endl
              << "*                                               *" << std::endl
              << "*    Written by Conrad Hübler TU Freiberg       *" << std::endl
              << "*                                               *" << std::endl
              << "*    Visit the website for initial usage        *" << std::endl
              << "*    https://github.com/conradhuebler/curcuma   *" << std::endl
              << "*                                               *" << std::endl
              << "*    This program comes without any warranty    *" << std::endl
              << "*    It might even be total useless             *" << std::endl
              << "*                                               *" << std::endl
              << "*    Nothing to cite yet ...                    *" << std::endl
              << "*    Git Commit Hash: " << git_commit_hash << "                   *" << std::endl
              //<< "*    Git Commit Date: " << git_date << "                   *" << std::endl
              << "*                                               *" << std::endl
              << "*************************************************" << std::endl;

    for (int index = 0; index < argc; ++index)
        std::cout << argv[index] << " ";
    std::cout << std::endl;
}

}

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

inline Molecule LoadFile(const std::string& filename)
{
    bool xyzfile = std::string(filename).find(".xyz") != std::string::npos || std::string(filename).find(".trj") != std::string::npos;

    if (xyzfile == false)
        throw 1;

    std::vector<std::string> lines;
    std::ifstream input(filename);

    int atoms = 0;
    int index = 0;
    int i = 0;

    Molecule mol(atoms, 0);
    for (std::string line; getline(input, line);) {
        if (index == 0 && xyzfile) {
            atoms = stoi(line);
            if (atoms < 1)
                throw 2;
            mol = Molecule(atoms, 0);
        }
        if (i > 1 && mol.AtomCount() < atoms) { // Load only the first file, TODO make it more flexible
            mol.setXYZ(line, i - 2);
        }
        index++;
        ++i;
    }
    return mol;
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

inline std::string Vector2String(const std::vector<int>& vector)
{
    std::string result = " ";

    for (auto i : vector)
        result += std::to_string(i) + "|";
    result.pop_back();

    result += " ";

    return result;
}

inline std::vector<int> String2Vector(const std::string& string)
{
    std::vector<int> vector;
    StringList tmp = SplitString(string, "|");
    for (auto i : tmp)
        vector.push_back(std::stoi(i));
    return vector;
}

inline std::string VectorVector2String(const std::vector<std::vector<int>>& vector)
{
    std::string result;
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

inline void xyz2allxyz(const std::string& xyzfile)
{
    std::string allxyz = xyzfile;
    allxyz.erase(allxyz.end() - 3, allxyz.end());

    std::ofstream input;
    input.open(allxyz + "allxyz", std::ios_base::app);
    FileIterator file(xyzfile);
    while (!file.AtEnd()) {
        Molecule mol = file.Next();
        input << mol.XYZString();
        if (!file.AtEnd())
            input << ">" << std::endl;
    }
    //input << "\n";
    input.close();
}
}
