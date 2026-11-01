/*
 * Built-in basis sets through the engine's own parsers (ctest qm_basis_*, Nov 2026).
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 *   dump_embedded_basis --list                       name, format, element range, number of elements, characters
 *   dump_embedded_basis --text NAME                  the assembled text of set NAME (all elements) on stdout
 *   dump_embedded_basis --write NAME DIR             the same text into DIR/NAME.dat ('$DATA') or .turbomole
 *   dump_embedded_basis --values NAME                parse the full text with the parser the engine uses for that
 *                                                    format (BasisSetParser for '$DATA', qmri for Turbomole) and
 *                                                    print "el L iprim exponent coefficient" with 17 digits
 *   dump_embedded_basis --time NAME REPEAT           mean parse time of the full text and of the H-Ne part, ms
 *   dump_embedded_basis --raw NAME Z                 the raw element block (e.g. def2-ECP, Z = 53)
 *
 * check_qm_basis_embed.py compares the --values output with the numbers it reads from --text (token -> double),
 * so the parsers are shown to read every embedded number exactly; the generator checks the text against BSE.
 *
 * Claude Generated (Nov 2026). GPL-3.0.
 */
#include "src/core/energy_calculators/qm_methods/basissetparser.hpp"
#include "src/core/energy_calculators/qm_methods/qm_basis_embedded.h"
#include "src/core/energy_calculators/qm_methods/qm_ri.h"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace {

const qmbasis::tables::EmbeddedBasisSet* findAny(const std::string& name)
{
    for (const char* f : { "data", "turbomole", "ecp" })
        if (const auto* s = qmbasis::findEmbedded(name, f)) return s;
    return nullptr;
}

const char* letterOf(BasisSetParser::ShellType t)
{
    static std::string s;
    s = BasisSetParser::shellTypeToString(t);
    return s.c_str();
}

double parseMs(const qmbasis::tables::EmbeddedBasisSet& s, int max_z, int repeat)
{
    const std::string text = qmbasis::embeddedSetText(s, max_z);
    size_t sink = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (int r = 0; r < repeat; ++r) {
        if (std::strcmp(s.format, "data") == 0) {
            sink += BasisSetParser::parseBasisSetString(text).size();
        } else {
            std::istringstream in(text);
            sink += qmri::parseTurbomoleStream(in, s.name).size();
        }
    }
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    if (sink == 0) std::fprintf(stderr, "nothing parsed\n");
    return ms / repeat;
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "usage: dump_embedded_basis --list | --text NAME | --write NAME DIR | --values NAME | --time NAME N | --raw NAME Z\n";
        return 1;
    }
    const std::string mode = argv[1];
    if (mode == "--list") {
        for (int i = 0; i < qmbasis::tables::kEmbeddedBasisSetCount; ++i) {
            const auto& s = qmbasis::tables::kEmbeddedBasisSets[i];
            std::printf("%s\t%s\t%s\t%zu\t%zu\t%s\n", s.name, s.format, s.elements, qmbasis::embeddedElements(s).size(),
                        qmbasis::embeddedSetText(s).size(), s.ecp);
        }
        return 0;
    }
    if (argc < 3) return 1;
    const auto* s = findAny(argv[2]);
    if (!s) {
        std::cerr << "no built-in set " << argv[2] << "\n";
        return 2;
    }
    if (mode == "--text") {
        std::fputs(qmbasis::embeddedSetText(*s).c_str(), stdout);
        return 0;
    }
    if (mode == "--write" && argc >= 4) {
        const std::string ext = std::strcmp(s->format, "data") == 0 ? ".dat" : ".turbomole";
        std::ofstream f(std::string(argv[3]) + "/" + s->name + ext);
        f << qmbasis::embeddedSetText(*s);
        return f.good() ? 0 : 3;
    }
    if (mode == "--values") {
        if (std::strcmp(s->format, "data") == 0) {
            const auto map = BasisSetParser::parseBasisSetString(qmbasis::embeddedSetText(*s));
            for (const auto& [sym, eb] : map)
                for (const auto& sh : eb.shells)
                    for (size_t k = 0; k < sh.exponents.size(); ++k)
                        std::printf("%s %s %zu %.17g %.17g\n", sym.c_str(), letterOf(sh.type), k, sh.exponents[k],
                                    sh.coefficients[0][k]);
        } else if (std::strcmp(s->format, "turbomole") == 0) {
            std::istringstream in(qmbasis::embeddedSetText(*s));
            const auto map = qmri::parseTurbomoleStream(in, s->name);
            for (const auto& [sym, shells] : map)
                for (const auto& sh : shells)
                    for (size_t k = 0; k < sh.exps.size(); ++k)
                        std::printf("%s %c %zu %.17g %.17g\n", sym.c_str(), "SPDFGHI"[sh.L], k, sh.exps[k], sh.coefs[k]);
        } else {
            std::cerr << "no parser for format " << s->format << " (raw access only: --raw)\n";
            return 4;
        }
        return 0;
    }
    if (mode == "--time" && argc >= 4) {
        const int repeat = std::max(1, std::atoi(argv[3]));
        const double full = parseMs(*s, 0, repeat), hne = parseMs(*s, 10, repeat);
        std::printf("%s full %.3f ms H-Ne %.3f ms (%zu / %zu characters)\n", s->name, full, hne,
                    qmbasis::embeddedSetText(*s).size(), qmbasis::embeddedSetText(*s, 10).size());
        return 0;
    }
    if (mode == "--raw" && argc >= 4) {
        const char* b = qmbasis::embeddedElementBlock(s->name, std::atoi(argv[3]));
        if (!b) return 5;
        std::fputs(b, stdout);
        return 0;
    }
    return 1;
}
