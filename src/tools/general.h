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

#include "src/core/global.h"
#include "src/core/molecule.h"

typedef std::vector<std::string> StringList;

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
            std::cout << "Finished after " << Elapsed() / 1000 << " seconds!" << std::endl;
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
    const char* delim = ".";
    StringList list = SplitString(input, delim);
    if (list.size() != 2)
        return false;

    bool left = isInt(list[0]);
    bool right = isInt(list[1]);
    return left && right;
    // return std::all_of(input.begin(), input.end(), ::isdigit);
}

inline Molecule LoadFile(const string& filename)
{

    bool xyzfile = std::string(filename).find(".xyz") != std::string::npos;

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
        if (i > 1) {
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

}
