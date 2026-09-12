/*
 * < Runtime-overridable GFN-FF parameter tables >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Sep 2026, rev-gfnff WP1a)
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
 */

#include "gfnff_param_tables.h"
#include "gfnff_par.h"

#include <functional>
#include <map>
#include <stdexcept>

namespace {

// name -> pointer-to-member for the gen scalars (one table drives merge, dump and help)
using GenField = double GFNFFGen::*;
const std::map<std::string, GenField>& genFields()
{
    static const std::map<std::string, GenField> f = {
        { "rabshift", &GFNFFGen::rabshift }, { "rabshifth", &GFNFFGen::rabshifth }, { "hyper_shift", &GFNFFGen::hyper_shift },
        { "hshift3", &GFNFFGen::hshift3 }, { "hshift4", &GFNFFGen::hshift4 }, { "hshift5", &GFNFFGen::hshift5 },
        { "hueckelp", &GFNFFGen::hueckelp }, { "bzref", &GFNFFGen::bzref }, { "hueckelp2", &GFNFFGen::hueckelp2 },
        { "bzref2", &GFNFFGen::bzref2 }, { "hueckelp3", &GFNFFGen::hueckelp3 }, { "qfacbm0", &GFNFFGen::qfacbm0 },
        { "srb1", &GFNFFGen::srb1 }, { "srb2", &GFNFFGen::srb2 }, { "srb3", &GFNFFGen::srb3 }, { "fcthr", &GFNFFGen::fcthr },
        { "atcuta", &GFNFFGen::atcuta }, { "atcutt", &GFNFFGen::atcutt }, { "atcutt_nci", &GFNFFGen::atcutt_nci },
        { "hbabmix", &GFNFFGen::hbabmix }, { "repscalb", &GFNFFGen::repscalb }, { "repscaln", &GFNFFGen::repscaln },
        { "qrepscal", &GFNFFGen::qrepscal }, { "nrepscal", &GFNFFGen::nrepscal }, { "hhfac", &GFNFFGen::hhfac },
        { "hh13rep", &GFNFFGen::hh13rep }, { "hh14rep", &GFNFFGen::hh14rep }, { "torsf_single", &GFNFFGen::torsf_single },
        { "torsf_pi", &GFNFFGen::torsf_pi }, { "torsf_improper", &GFNFFGen::torsf_improper },
        { "torsf_pi_improper", &GFNFFGen::torsf_pi_improper }, { "torsf_extra_C", &GFNFFGen::torsf_extra_C },
        { "torsf_extra_N", &GFNFFGen::torsf_extra_N }, { "torsf_extra_O", &GFNFFGen::torsf_extra_O },
        { "fr3", &GFNFFGen::fr3 }, { "fr4", &GFNFFGen::fr4 }, { "fr5", &GFNFFGen::fr5 }, { "fr6", &GFNFFGen::fr6 },
    };
    return f;
}

using TableField = std::vector<double> GFNFFTables::*;
const std::map<std::string, TableField>& tableFields()
{
    static const std::map<std::string, TableField> f = {
        { "chi_eeq", &GFNFFTables::chi_eeq }, { "gam_eeq", &GFNFFTables::gam_eeq }, { "alpha_eeq", &GFNFFTables::alpha_eeq },
        { "cnf_eeq", &GFNFFTables::cnf_eeq }, { "bond_params", &GFNFFTables::bond_params }, { "r0_gfnff", &GFNFFTables::r0_gfnff },
        { "cnfak_gfnff", &GFNFFTables::cnfak_gfnff }, { "en_gfnff", &GFNFFTables::en_gfnff }, { "en_rab_gfnff", &GFNFFTables::en_rab_gfnff },
        { "angle_params", &GFNFFTables::angle_params }, { "angl2_neighbors", &GFNFFTables::angl2_neighbors },
        { "repa", &GFNFFTables::repa }, { "repan", &GFNFFTables::repan }, { "repz", &GFNFFTables::repz },
        { "tors", &GFNFFTables::tors }, { "tors2", &GFNFFTables::tors2 },
    };
    return f;
}

void applyTable(std::vector<double>& table, const std::string& name, const json& value)
{
    if (value.is_array()) {
        if (value.size() != table.size())
            throw std::runtime_error("GFN-FF parameter override: table '" + name + "' expects " + std::to_string(table.size())
                + " entries, got " + std::to_string(value.size()));
        for (size_t i = 0; i < table.size(); ++i)
            table[i] = value[i].get<double>();
        return;
    }
    if (!value.is_object())
        throw std::runtime_error("GFN-FF parameter override: table '" + name + "' must be an array or a {\"Z\": value} object");
    for (const auto& [key, v] : value.items()) {
        int z = 0;
        try {
            z = std::stoi(key);
        } catch (...) {
            throw std::runtime_error("GFN-FF parameter override: table '" + name + "' key '" + key + "' is not an atomic number");
        }
        if (z < 1 || z > static_cast<int>(table.size()))
            throw std::runtime_error("GFN-FF parameter override: table '" + name + "' has no entry for Z=" + key);
        table[z - 1] = v.get<double>();
    }
}

} // namespace

std::shared_ptr<const GFNFFTables> GFNFFTables::defaults()
{
    static const std::shared_ptr<const GFNFFTables> instance = [] {
        namespace P = GFNFFParameters; // qualified: the unqualified names would resolve to the members
        auto t = std::make_shared<GFNFFTables>();
        t->chi_eeq = P::chi_eeq;
        t->gam_eeq = P::gam_eeq;
        t->alpha_eeq = P::alpha_eeq;
        t->cnf_eeq = P::cnf_eeq;
        t->bond_params = P::bond_params;
        t->r0_gfnff = P::r0_gfnff;
        t->cnfak_gfnff = P::cnfak_gfnff;
        t->en_gfnff = P::en_gfnff;
        t->en_rab_gfnff = P::en_rab_gfnff;
        t->angle_params = P::angle_params;
        t->angl2_neighbors = P::angl2_neighbors;
        t->repa = P::repa_angewChem2020;
        t->repan = P::repan_angewChem2020;
        t->repz = P::repz;
        t->tors = P::tors_angewChem2020;
        t->tors2 = P::tors2_angewChem2020;
        // gen scalars that gfnff_par.h defines (the rest are the literals of the generator code)
        t->gen.repscalb = P::REPSCALB;
        t->gen.repscaln = P::REPSCALN;
        t->gen.qrepscal = P::QREPSCAL;
        t->gen.nrepscal = P::NREPSCAL;
        t->gen.hhfac = P::HHFAC;
        t->gen.hh13rep = P::HH13REP;
        t->gen.hh14rep = P::HH14REP;
        t->gen.atcutt = P::atcutt;
        t->gen.atcutt_nci = P::atcutt_nci;
        t->gen.torsf_single = P::torsf_single;
        t->gen.torsf_pi = P::torsf_pi;
        t->gen.torsf_improper = P::torsf_improper;
        t->gen.torsf_pi_improper = P::torsf_pi_improper;
        t->gen.torsf_extra_C = P::torsf_extra_C;
        t->gen.torsf_extra_N = P::torsf_extra_N;
        t->gen.torsf_extra_O = P::torsf_extra_O;
        t->gen.fr3 = P::FR3;
        t->gen.fr4 = P::FR4;
        t->gen.fr5 = P::FR5;
        t->gen.fr6 = P::FR6;
        for (int i = 0; i < 9; ++i)
            t->gen.bstren[i] = P::bstren[i];
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                t->gen.bsmat[i][j] = P::bsmat[i][j];
        return t;
    }();
    return instance;
}

std::shared_ptr<const GFNFFTables> GFNFFTables::fromOverrides(const json& overrides)
{
    if (!overrides.is_object() || overrides.empty())
        return defaults();
    auto t = std::make_shared<GFNFFTables>(*defaults());
    for (const auto& [section, body] : overrides.items()) {
        if (section == "gen") {
            for (const auto& [name, v] : body.items()) {
                auto it = genFields().find(name);
                if (it == genFields().end())
                    throw std::runtime_error("GFN-FF parameter override: unknown gen scalar '" + name + "'");
                t->gen.*(it->second) = v.get<double>();
            }
        } else if (section == "tables") {
            for (const auto& [name, v] : body.items()) {
                auto it = tableFields().find(name);
                if (it == tableFields().end())
                    throw std::runtime_error("GFN-FF parameter override: unknown table '" + name + "'");
                applyTable(t.get()->*(it->second), name, v);
            }
        } else if (section == "rev") {
            if (!body.is_object())
                throw std::runtime_error("GFN-FF parameter override: 'rev' must be an object");
            for (const auto& [name, v] : body.items())
                t->rev[name] = v;
        } else {
            throw std::runtime_error("GFN-FF parameter override: unknown section '" + section + "' (expected gen, tables, rev)");
        }
    }
    t->overrides = overrides;
    t->hash = std::to_string(std::hash<std::string> {}(overrides.dump()));
    return t;
}

json GFNFFTables::toJSON() const
{
    json j;
    for (const auto& [name, field] : genFields())
        j["gen"][name] = gen.*field;
    j["gen"]["bstren"] = std::vector<double>(gen.bstren.begin(), gen.bstren.end());
    for (const auto& [name, field] : tableFields())
        j["tables"][name] = this->*field;
    j["rev"] = rev;
    j["hash"] = hash;
    return j;
}

std::vector<std::string> GFNFFTables::genNames()
{
    std::vector<std::string> out;
    for (const auto& [name, _] : genFields())
        out.push_back(name);
    return out;
}

std::vector<std::string> GFNFFTables::tableNames()
{
    std::vector<std::string> out;
    for (const auto& [name, _] : tableFields())
        out.push_back(name);
    return out;
}
