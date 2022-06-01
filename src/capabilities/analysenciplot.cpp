/*
 * <Post-process NCI RDG vs sign(lambda2)rho plots>
 * Copyright (C) 2021 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <map>

#include "src/tools/general.h"

#include "analysenciplot.h"

AnalyseNCIPlot::AnalyseNCIPlot(const json& controller)
    : CurcumaMethod(AnalyseNCIPlotJson, controller, false)

{
    UpdateController(controller);

    std::cout << "NCIPlot Analysis" << std::endl;
}

void AnalyseNCIPlot::start()
{
    int scaling = m_bins;
    std::cout << "NCIPlot Analysis - Loading files" << std::endl;

    m_NCI1 = LoadFile(m_file1);
    m_NCI2 = LoadFile(m_file2);

    std::string file1 = m_file1;
    std::string file2 = m_file2;

    file1.erase(file1.end() - 4, file1.end());
    file2.erase(file2.end() - 4, file2.end());

    file1 = file1 + ".depleted.dat";
    file2 = file2 + ".depleted.dat";

    std::vector<double> distances_1, distances_2;
    std::map<int, std::map<double, double>::iterator> index_mapper_1, index_mapper_2;
    for (auto iter = m_NCI1.begin(); iter != m_NCI1.end(); ++iter) {

        if (index_mapper_1.find(int(iter->first * scaling)) == index_mapper_1.end())
            index_mapper_1.insert(std::pair<int, std::map<double, double>::iterator>(int(iter->first * scaling), iter));

        auto next = iter;
        next++;
        if (next == m_NCI1.end())
            break;

        distances_1.push_back(sqrt((iter->first - next->first) * (iter->first - next->first) + (iter->second - next->second) * (iter->second - next->second)));
    }
    index_mapper_1.insert(std::pair<int, std::map<double, double>::iterator>(index_mapper_1.size(), m_NCI1.end()));

    for (auto iter = m_NCI2.begin(); iter != m_NCI2.end(); ++iter) {

        if (index_mapper_2.find(int(iter->first * scaling)) == index_mapper_2.end())
            index_mapper_2.insert(std::pair<int, std::map<double, double>::iterator>(int(iter->first * scaling), iter));

        auto next = iter;
        next++;
        if (next == m_NCI2.end())
            break;

        distances_2.push_back(sqrt((iter->first - next->first) * (iter->first - next->first) + (iter->second - next->second) * (iter->second - next->second)));
    }
    index_mapper_2.insert(std::pair<int, std::map<double, double>::iterator>(index_mapper_2.size(), m_NCI2.end()));

    double mean1 = Tools::mean(distances_1);
    double mean2 = Tools::mean(distances_2);

    std::map<double, double> final, left, right;

    int count = 0;

    for (auto iter_i = m_NCI1.begin(); iter_i != m_NCI1.end(); ++iter_i) {
        bool match = false;

        int index_start = int(iter_i->first * scaling);
        while (index_mapper_2.find(index_start) == index_mapper_2.end())
            index_start--;
        int index_end = int(iter_i->first * scaling) + 1;

        if (index_end <= m_NCI2.size())
            break;

        while (index_mapper_2.find(index_end) == index_mapper_2.end())
            index_end++;
        if (m_local_distance) {
            mean1 = getDistance(m_NCI1, index_mapper_2[index_start], index_mapper_2[index_end]);
            std::cout << mean1 << std::endl;
        }
        for (auto iter_j = index_mapper_2[index_start]; iter_j != index_mapper_2[index_end]; ++iter_j) {
            double d = (sqrt((iter_i->first - iter_j->first) * (iter_i->first - iter_j->first) + (iter_i->second - iter_j->second) * (iter_i->second - iter_j->second)));
            if (d < mean1 / m_scale_d1) {
                match = true;
                count++;
                break;
            } else if (d > 1) {
                match = true;
                break;
            }
        }
        if (!match) {
            final.insert(std::pair<double, double>(iter_i->first, iter_i->second));
            left.insert(std::pair<double, double>(iter_i->first, iter_i->second));
        }
    }

    for (auto iter_i = m_NCI2.begin(); iter_i != m_NCI2.end(); ++iter_i) {
        bool match = false;

        int index_start = int(iter_i->first * scaling);
        while (index_mapper_1.find(index_start) == index_mapper_1.end())
            index_start--;
        int index_end = int(iter_i->first * scaling) + 1;

        if (index_end <= m_NCI1.size())
            break;

        while (index_mapper_1.find(index_end) == index_mapper_1.end())
            index_end++;

        if (m_local_distance) {
            mean2 = getDistance(m_NCI2, index_mapper_1[index_start], index_mapper_1[index_end]);
        }

        for (auto iter_j = index_mapper_1[index_start]; iter_j != index_mapper_1[index_end]; ++iter_j) {
            double d = (sqrt((iter_i->first - iter_j->first) * (iter_i->first - iter_j->first) + (iter_i->second - iter_j->second) * (iter_i->second - iter_j->second)));
            if (d < mean2 / m_scale_d2) {
                match = true;
                count++;
                break;
            } else if (d > 1) {
                match = true;
                break;
            }
        }
        if (!match) {
            final.insert(std::pair<double, double>(iter_i->first, iter_i->second));
            right.insert(std::pair<double, double>(iter_i->first, iter_i->second));
        }
    }

    std::ofstream input;
    input.open("combined.dat", std::ios::out);

    for (auto i : final)
        input << i.first << " " << i.second << std::endl;
    input.close();

    input.open(file1, std::ios::out);

    for (auto i : left)
        input << i.first << " " << i.second << std::endl;
    input.close();

    input.open(file2, std::ios::out);

    for (auto i : right)
        input << i.first << " " << i.second << std::endl;
    input.close();
    std::cout << "NCIPlot Analysis - Loading files done" << std::endl;
}
std::map<double, double> AnalyseNCIPlot::LoadFile(const std::string& file)
{
    std::map<double, double> map;
    std::ifstream input(file);
    std::vector<int> ignore_index;
    int index = 0;
    for (std::string line; getline(input, line);) {
        index++;
        StringList list = Tools::SplitString(line, " ");
        if (list.size() == 0)
            continue;
        if (std::stod(list[0]) >= -0.07 && std::stod(list[0]) <= 0.07)
            map.insert(std::pair<double, double>(std::stod(list[0]), std::stod(list[1])));
    }
    return map;
}

void AnalyseNCIPlot::LoadControlJson()
{
    m_bins = Json2KeyWord<double>(m_defaults, "bins");
    m_scale_d1 = Json2KeyWord<double>(m_defaults, "scale_d1");
    m_scale_d2 = Json2KeyWord<double>(m_defaults, "scale_d2");

    m_local_distance = Json2KeyWord<bool>(m_defaults, "local_distance");
}
double AnalyseNCIPlot::getDistance(const std::map<double, double>& map, const std::map<double, double>::iterator& begin, const std::map<double, double>::iterator& end)
{
    std::vector<double> distances;
    for (auto iter = begin; iter != end; ++iter) {
        auto next = iter;
        next++;
        if (next == map.end())
            break;

        distances.push_back(sqrt((iter->first - next->first) * (iter->first - next->first) + (iter->second - next->second) * (iter->second - next->second)));
    }
    return Tools::mean(distances);
}
