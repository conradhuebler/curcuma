/*
 * <Native QM: lookup of the built-in basis sets>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Nov 2026). AI-generated, machine-tested.
 *
 * Built in (GENERATED headers qm_basis_data_*.h and qm_basis_tables.h, written by scripts/qm_basis_export_tables.py
 * from Basis Set Exchange data at a pinned commit, see those headers and docs/THIRD_PARTY_NOTICES.md):
 *   orbital ('$DATA' text, BasisSetParser): def2-SVP, -SV(P), -SVPD, -TZVP, -TZVPP, -TZVPD, -TZVPPD, -QZVP, -QZVPP,
 *     -QZVPD, -QZVPPD, -mTZVP, -mTZVPP (every element BSE has), MINIX (H-Ne)
 *   fitting (Turbomole text, qmri::parseTurbomoleBasis): def2-universal-JFIT, -JKFIT, def2-SV(P)-JKFIT,
 *     def2-mTZVPP-RIJ, def2-{SVP,SV(P),SVPD,TZVP,TZVPP,TZVPD,TZVPPD,QZVP,QZVPP,QZVPPD}-RIFIT
 *   ECP (Turbomole '$ecp' text): def2-ECP, Rb-Rn; read by qm_ecp.cpp (embeddedElementBlock("def2-ECP", Z)).
 * The engine takes H-Rn (kEngineMaxZ) with s, p, d, f, g shells (createGTOFromBasis; f/g since Nov 2026, Z > 10
 * since Nov 2026); MINIX stays H-Ne.
 *
 * Lookup order of QMEngine (orbital, `-qm.basis`) and qmri::resolveAuxBasisFile (`-qm.ri_aux`, `-qm.ri_mp2_aux`):
 *   1. an explicit path (contains '/') that exists
 *   2. $CURCUMA_QM_BASIS (old name $CURCUMA_DFT_BASIS), $CURCUMA_DATA
 *   3. <CURCUMA_SOURCE_DIR>/src/core/energy_calculators/qm_methods/ (compile-time path)
 *   4. the current directory
 *   5. the built-in copy, by exact basis-set name (LAST: a file always wins)
 * CURCUMA_QM_BASIS_EMBEDDED_ONLY=1 skips steps 1-4 (test hook; ctest qm_basis_*).
 *
 * A resolved built-in set is passed around as the pseudo path "[embedded]/<name>.dat" or ".turbomole";
 * embeddedText() assembles its text (header + element blocks + footer; the engine asks for the blocks up to parseMaxZ(atoms))
 * and the parsers read it through the same stream parser as a file; for H-Ne the parse result equals that of a
 * file written from the full text (ctest qm_basis_sets_*). The data are stored per element and live in
 * one translation unit (qm_basis_embedded.cpp), so only that file compiles the several MB of text.
 */

#pragma once

#include <string>
#include <vector>

namespace qmbasis {

namespace tables {
/// One built-in set (filled by the generated qm_basis_tables.h).
struct EmbeddedBasisSet {
    const char* name;           ///< exact name as used by -qm.basis / -qm.ri_aux (def2-TZVP, def2-TZVP-RIFIT, ...)
    const char* format;         ///< "data" ('$DATA', BasisSetParser), "turbomole" (qmri), "ecp" (Turbomole $ecp)
    const char* header;         ///< text before the first element block
    const char* footer;         ///< text after the last element block
    const char* const* blocks;  ///< blocks[Z-1]: text of element Z, nullptr if the set does not have it
    int max_z;                  ///< length of `blocks`
    const char* elements;       ///< element range, e.g. "H-La, Hf-Rn"
    const char* ecp;            ///< ECP the BSE combines with the set for heavy elements ("def2-ECP") or ""
};
extern const EmbeddedBasisSet kEmbeddedBasisSets[];
extern const int kEmbeddedBasisSetCount;
}  // namespace tables

inline constexpr const char kEmbeddedPrefix[] = "[embedded]/";

/// CURCUMA_QM_BASIS_EMBEDDED_ONLY set to anything but "" or "0": no file lookup, built-in sets only.
bool embeddedOnly();

/// Built-in set with exactly this name and format ("data" orbital, "turbomole" fitting, "ecp"), else nullptr.
const tables::EmbeddedBasisSet* findEmbedded(const std::string& name, const char* format);

/// As findEmbedded, but the name compared without letter case (aliases of the fitting sets, qm_ri.cpp).
const tables::EmbeddedBasisSet* findEmbeddedNoCase(const std::string& name, const char* format);

/// Names of the built-in sets of one format, comma separated (for error messages and help).
std::string embeddedNames(const char* format);

bool isEmbeddedPath(const std::string& path);

/// Pseudo path of a built-in set: "[embedded]/<name><extension>".
std::string embeddedPath(const tables::EmbeddedBasisSet& s, const char* extension);

/// Full text of a set: header, the blocks of the elements 1..max_z (0 = all), footer.
std::string embeddedSetText(const tables::EmbeddedBasisSet& s, int max_z = 0);

/// Highest Z the engine handles (H-Rn, the def2 range; Nov 2026, H-Ne before). The parsers get only the blocks up
/// to the highest Z of the molecule (parseMaxZ): parsing the full text costs up to 2.1 ms (def2-universal-JKFIT,
/// 86 elements) against 0.1 ms for H-Ne, per RI setup.
inline constexpr int kEngineMaxZ = 86;

/// Number of element blocks a parser needs for these atoms: max(10, highest Z), capped at kEngineMaxZ. Never below
/// 10, so an H-Ne molecule reads exactly the text it read before Nov 2026. Claude Generated (Nov 2026).
inline int parseMaxZ(const std::vector<int>& atoms)
{
    int z = 10;
    for (int a : atoms) z = a > z ? a : z;
    return z > kEngineMaxZ ? kEngineMaxZ : z;
}

/// Text of the built-in set behind a pseudo path into `text` (elements 1..max_z, 0 = all); false if the path is
/// not one. QMEngine and qmri pass parseMaxZ(atoms).
bool embeddedText(const std::string& path, std::string& text, int max_z = 0);

/// Raw, parser-agnostic text of element Z in the built-in set `name` (any format), nullptr if absent.
/// For "def2-ECP" this is the Turbomole '$ecp' block ("ncore = .. lmax = ..", potentials), read by QMEngine::setupECP.
const char* embeddedElementBlock(const std::string& name, int Z);

/// Atomic numbers present in a built-in set (ascending).
std::vector<int> embeddedElements(const tables::EmbeddedBasisSet& s);

}  // namespace qmbasis
