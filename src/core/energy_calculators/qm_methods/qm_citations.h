/*
 * <Citations of the native ab-initio QM basis sets>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Oct 2026): registers the literature of the basis sets a native QM run actually loads, so that
 * CitationRegistry prints it at program end and writes it to curcuma_citations.bib. The entries live in
 * src/core/citation_database.cpp (DOIs checked against Crossref). Attribution is by basis set NAME; a user-supplied
 * file with another name is not cited (nothing is guessed). The shipped files are Basis Set Exchange downloads
 * (def2-*, JFIT/JKFIT) or the ORCA export of MINIX (HF-3c), see their headers.
 */

#pragma once

#include "src/core/citation_registry.h"

#include <algorithm>
#include <cctype>
#include <string>

namespace qmcite {

inline std::string lowerStem(const std::string& name)
{
    std::string s = name.substr(name.find_last_of('/') == std::string::npos ? 0 : name.find_last_of('/') + 1);
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    for (const char* ext : { ".dat", ".turbomole" })
        if (s.size() > std::string(ext).size() && s.compare(s.size() - std::string(ext).size(), std::string(ext).size(), ext) == 0)
            s.erase(s.size() - std::string(ext).size());
    return s;
}

inline void bse()
{
    CitationRegistry::cite("bse_2019");
    CitationRegistry::cite("bse_2007");
    CitationRegistry::cite("bse_1996");
}

/// Orbital basis set by name or file (def2-SV(P), def2-SVP, ...; MINIX).
inline void orbitalBasis(const std::string& name_or_file)
{
    const std::string s = lowerStem(name_or_file);
    if (s.rfind("def2-", 0) == 0) {
        CitationRegistry::cite("def2_basis");
        if (s.rfind("def2-sv", 0) == 0)
            CitationRegistry::cite("def2_sv_1992");
        bse();
    } else if (s == "minix") {
        CitationRegistry::cite("hf3c");   // MINIX is defined in the HF-3c paper
    }
}

/// Auxiliary (density fitting) basis set by name or file (def2-universal-JFIT / -JKFIT).
inline void auxBasis(const std::string& name_or_file)
{
    const std::string s = lowerStem(name_or_file);
    CitationRegistry::cite("ri_hf");
    if (s.find("jkfit") != std::string::npos) {
        CitationRegistry::cite("def2_jkfit");
        bse();
    } else if (s.find("jfit") != std::string::npos) {
        CitationRegistry::cite("def2_jfit");
        bse();
    }
}

/// Closed-shell HF/SCF machinery of the native engine: Roothaan-Hall equations, McMurchie-Davidson integrals, Pulay DIIS.
inline void hf()
{
    CitationRegistry::cite("roothaan");
    CitationRegistry::cite("md_integrals");
    CitationRegistry::cite("diis");
}

/// Integral-direct SCF with Schwarz / density-weighted screening.
inline void directScf()
{
    CitationRegistry::cite("direct_scf");
    CitationRegistry::cite("schwarz_screening");
}

/// MP2 variants by method name: mp2, scs-mp2, sos-mp2, ri-mp2, scs-ri-mp2, sos-ri-mp2.
inline void mp2(const std::string& variant, bool gradient = false)
{
    CitationRegistry::cite("mp2");
    if (variant.find("scs") != std::string::npos) CitationRegistry::cite("scs_mp2");
    if (variant.find("sos") != std::string::npos) CitationRegistry::cite("sos_mp2");
    if (variant.find("ri-mp2") != std::string::npos) {
        CitationRegistry::cite("ri_mp2_feyereisen");
        CitationRegistry::cite("ri_vahtras");
        CitationRegistry::cite("ri_mp2_weigend");
    }
    if (gradient) CitationRegistry::cite("zvector");
}

/// Canonical closed-shell CCSD.
inline void ccsd()
{
    CitationRegistry::cite("ccsd_purvis");
    CitationRegistry::cite("ccsd_scuseria");
    CitationRegistry::cite("ccsd_stanton");
}

/// Local PNO-RI-MP2: localiser "pm" or "boys", Boughton-Pulay domains, pair natural orbitals.
inline void lmp2(const std::string& localizer)
{
    mp2("ri-mp2");
    CitationRegistry::cite("lmp2_pulay");
    CitationRegistry::cite("lmp2_saebo");
    CitationRegistry::cite("boughton_pulay");
    CitationRegistry::cite(localizer == "boys" ? "foster_boys" : "pipek_mezey");
    CitationRegistry::cite("pno_edmiston");
    CitationRegistry::cite("pno_meyer");
    CitationRegistry::cite("pno_neese");
    CitationRegistry::cite("dlpno_mp2");
}

/// HF-3c composite (HF/MINIX, D3(BJ), gCP, SRB).
inline void hf3c()
{
    CitationRegistry::cite("hf3c");
    CitationRegistry::cite("gcp");
    CitationRegistry::cite("d3");
    CitationRegistry::cite("bj");
}

} // namespace qmcite
