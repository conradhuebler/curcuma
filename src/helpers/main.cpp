/*
 * <Internal Coordinate Handler for chemical structures.>
 * Copyright (C) 2019 - 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/tools/general.h"
#include "src/tools/info.h"

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
    General::StartUp(argc, argv);

    RunTimer timer(true);

    if (argc < 2) {
        std::cerr << "No arguments given!" << std::endl;
        std::cerr << "Use:" << std::endl
                  << "-statistic   * Some statistical calulation on rows of numbers *" << std::endl;
    }

    if (argc >= 2) {
        if (strcmp(argv[1], "-statistic") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma for statistical tools follows\ncurcuma_tools -statistic file.dat" << std::endl;
                std::cerr << "Additonal arguments are:" << std::endl;
                std::cerr << "-bins       **** Number of bins to use (default = 100)." << std::endl;
                std::cerr << "-force      **** Force append, even if last lines exist." << std::endl;

                exit(1);
            }
            int bins = 100;
            bool force = false;
            if (argc >= 3) {
                for (std::size_t i = 3; i < argc; ++i) {
                    if (strcmp(argv[i], "-bins") == 0) {
                        if (i + 1 < argc) {
                            bins = std::stoi(argv[i + 1]);
                            ++i;
                        }
                        continue;
                    }
                    if (strcmp(argv[i], "-force") == 0) {
                        force = true;
                        continue;
                    }
                }
            }
            std::vector<std::vector<double>> storage;
            std::vector<std::string> lines;
            std::ifstream input(argv[2]);
            std::vector<int> ignore_index;
            int index = 0;
            int columns = 0;
            for (std::string line; getline(input, line);) {
                index++;
                StringList list = Tools::SplitString(line, "\t");
                if (list.size() == 0)
                    continue;
                if (list[0].find("#") != std::string::npos) {
                    ignore_index.push_back(index);
                    continue;
                }
                if (storage.size() == 0) {
                    storage = std::vector<std::vector<double>>(list.size());
                }
                if (storage.size() == list.size()) {
                    int sub = 0;
                    std::vector<double> row;
                    for (auto element : list) {
                        double val = 0.0;
                        bool ok = true;
                        try {
                            val = std::stod(element);
                        } catch (const std::string& what_arg) {
                            ok = false;
                        } catch (const std::invalid_argument& arg) {
                            ok = false;
                        }
                        if (ok) {
                            storage[sub].push_back(val);
                            sub++;
                        }
                    }
                }
            }

            if (ignore_index.size() == 0 || force) {
                std::ofstream input;
                input.open(argv[2], std::ios_base::app);

                std::vector<double> mean, median, stdev, entropy;

                for (std::vector<double> vector : storage) {
                    std::sort(vector.begin(), vector.end());
                    double m = Tools::mean(vector);
                    double stde = Tools::stdev(vector, m);
                    auto hist = Tools::Histogram(vector, bins);
                    double shannon = Tools::ShannonEntropy(hist);
                    double med = Tools::median(vector);

                    mean.push_back(m);
                    median.push_back(med);
                    stdev.push_back(stde);
                    entropy.push_back(shannon);
                }

                input << "#";

                for (auto val : mean) {
                    input << std::setprecision(6) << val << "   ";
                }
                input << std::endl
                      << "#";

                for (auto val : median) {
                    input << std::setprecision(6) << val << "   ";
                }
                input << std::endl
                      << "#";

                for (auto val : stdev) {
                    input << std::setprecision(6) << val << "   ";
                }
                input << std::endl
                      << "#";

                for (auto val : entropy) {
                    input << std::setprecision(6) << val << "   ";
                }
                input << std::endl;
            }
        } else if (strcmp(argv[1], "-allxyz") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma for converting xyz files into allxyz files (orca format) as follows\ncurcuma_tools -allxyz input.xyz" << std::endl;
                exit(1);
            }
            Tools::xyz2allxyz(argv[2]);
        }
    }
    return 0;
}
