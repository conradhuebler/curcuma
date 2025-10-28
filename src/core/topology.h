/*
 * <Topology Stuff. >
 * Copyright (C) 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <algorithm>
#include <vector>

namespace Topology {

inline double Distance(double x1, double x2, double y1, double y2, double z1, double z2)
{
    return sqrt((((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)) + ((z1 - z2) * (z1 - z2))));
}

inline std::vector<std::vector<int>> FindRings(const std::vector<std::vector<int>>& stored_bonds, int atoms, int maxsize = 30, int storage = 10)
{
    std::vector<std::vector<int>> identified_rings;
    std::vector<int> done;

    for (int i = 0; i < atoms; ++i) {
        //  if (std::find(done.begin(), done.end(), i) != done.end())
        //      continue;
        if (stored_bonds[i].size() == 1) {
            done.push_back(i);
            continue;
        }
        bool loop = true;
        std::vector<int> bonded = stored_bonds[i];
        std::vector<int> knots, outer;
        // std::vector< std::pair<int, std::vector<int> > > mknots;
        /*
        for (int atom : bonded) {
            if (stored_bonds[atom].size() == 1)
                done.push_back(atom);
            knots.push_back(atom);
            outer.push_back(atom);
        }*/
        std::vector<std::vector<int>> stash;
        stash.push_back(std::vector<int>{ i });
        // done.push_back(std::vector<int>());
        int index = -1;
        while (stash.size() && done.size() < storage * atoms) { // this is just a hack ...
            for (auto tmp : knots) {
                auto it = std::find(done.begin(), done.end(), tmp);
                if (it != done.end())
                    done.erase(it);
            }
            for (int s = 0; s < stash.size(); ++s) {
                int outeratom = stash[s][stash[s].size() - 1];
                {
                    std::vector<int> bonded = stored_bonds[outeratom];
                    std::vector<int> vacant;
                    bool close_ring = false;
                    for (int atom : bonded) {
                        // std::cout << atom << " " << stash[s][stash[s].size() - 2] << std::endl;
                        // Claude Fix 2025-10-19: Prevent out-of-bounds access when stash[s].size() < 2
                        // stash[s] starts with size 1 (line 56), so size()-2 would be -1 (18446744073709551615 unsigned)
                        if (stash[s].size() >= 2 && stash[s][stash[s].size() - 2] == atom)
                            continue;

                        if (stash[s][0] == atom) {
                            vacant.push_back(atom);
                            close_ring = true;
                            break;
                        }
                        if (stored_bonds[atom].size() == 1) // ignore atoms with only one bond
                        {
                            done.push_back(atom);
                            continue;
                        }
                        int counter = 0;
                        for (auto a : stash[s]) // check if atom is already in the list
                            counter += a == atom;
                        bool isKnot = std::find(knots.begin(), knots.end(), atom) != knots.end(); // check if atom is a knot
                        if (counter) // if atom is in list
                        {
                            if (isKnot) // take care of knot-atoms
                            {
                                if (counter >= stored_bonds[atom].size()) // if fewer counts than binding partners, except
                                    continue;
                            } else
                                continue;
                        }

                        if (std::find(done.begin(), done.end(), atom) != done.end()) {
                            continue;
                        }

                        vacant.push_back(atom);
                        if (stash[s].size() == 1)
                            break;
                    }
                    if (vacant.size() == 0 || close_ring) {
                        loop = false;
                        if (stash[s].size() < 3) {
                            stash.erase(stash.begin() + s);
                            break;
                        }
                        int first = stash[s][0];
                        int last = stash[s][stash[s].size() - 1];
                        bool connected = false;
                        for (int a : stored_bonds[first]) {
                            if (a == last) {
                                connected = true;
                                break;
                            }
                        }
                        if (connected) {
                            index = s;
                            auto tmp = stash[s];
                            std::sort(tmp.begin(), tmp.end());
                            identified_rings.push_back(tmp);
                            for (int a : stash[s]) {
                                if (stash[s].size() < maxsize)
                                    done.push_back(a);
                                // std::cout << a << " ";
                            }
                            // std::cout << std::endl;
                            stash.erase(stash.begin() + s);
                        } else {
                            //    for (int a : stash[s])
                            //        done.push_back(a);
                            stash.erase(stash.begin() + s);
                        }
                    } else if (vacant.size() == 1)
                        stash[s].push_back(vacant[0]);
                    else {
                        auto currentstash = stash[s];
                        stash.erase(stash.begin() + s);
                        if (currentstash.size() > maxsize)
                            break;
                        // auto currdone = done[s];
                        // done.erase(done.begin() + s);
                        knots.push_back(outeratom);
                        for (int atom : vacant) {
                            stash.push_back(currentstash);
                            stash[stash.size() - 1].push_back(atom);
                            //  done.push_back(currdone);
                        }
                    }
                }
            }
        }
    }
    for (auto a : identified_rings) {
        // if (a.size() < 10) {
        // for (auto i : a)
        //    std::cout << i << " ";
        // std::cout << std::endl;
        //}
    }
    return identified_rings;
}

}
