/*
 * <Native QM: lookup of the built-in basis sets>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Nov 2026). AI-generated, machine-tested.
 *
 * The only translation unit that includes the generated data (qm_basis_tables.h -> qm_basis_data_*.h); see
 * qm_basis_embedded.h for the lookup order and the list of sets.
 */

#include "qm_basis_embedded.h"

#include "qm_basis_tables.h"

#include <cctype>
#include <cstdlib>
#include <cstring>

namespace qmbasis {

bool embeddedOnly()
{
    const char* v = std::getenv("CURCUMA_QM_BASIS_EMBEDDED_ONLY");
    return v != nullptr && v[0] != '\0' && std::strcmp(v, "0") != 0;
}

const tables::EmbeddedBasisSet* findEmbedded(const std::string& name, const char* format)
{
    for (int i = 0; i < tables::kEmbeddedBasisSetCount; ++i) {
        const auto& s = tables::kEmbeddedBasisSets[i];
        if (name == s.name && std::strcmp(format, s.format) == 0) return &s;
    }
    return nullptr;
}

const tables::EmbeddedBasisSet* findEmbeddedNoCase(const std::string& name, const char* format)
{
    auto low = [](std::string t) {
        for (char& c : t) c = (char)std::tolower((unsigned char)c);
        return t;
    };
    const std::string ln = low(name);
    for (int i = 0; i < tables::kEmbeddedBasisSetCount; ++i) {
        const auto& s = tables::kEmbeddedBasisSets[i];
        if (ln == low(s.name) && std::strcmp(format, s.format) == 0) return &s;
    }
    return nullptr;
}

std::string embeddedNames(const char* format)
{
    std::string out;
    for (int i = 0; i < tables::kEmbeddedBasisSetCount; ++i) {
        const auto& s = tables::kEmbeddedBasisSets[i];
        if (std::strcmp(format, s.format) == 0) out += (out.empty() ? "" : ", ") + std::string(s.name);
    }
    return out;
}

bool isEmbeddedPath(const std::string& path) { return path.rfind(kEmbeddedPrefix, 0) == 0; }

std::string embeddedPath(const tables::EmbeddedBasisSet& s, const char* extension)
{
    return std::string(kEmbeddedPrefix) + s.name + extension;
}

std::string embeddedSetText(const tables::EmbeddedBasisSet& s, int max_z)
{
    const int n = (max_z <= 0 || max_z > s.max_z) ? s.max_z : max_z;
    size_t len = std::strlen(s.header) + std::strlen(s.footer);
    for (int z = 0; z < n; ++z)
        if (s.blocks[z]) len += std::strlen(s.blocks[z]);
    std::string text;
    text.reserve(len);
    text += s.header;
    for (int z = 0; z < n; ++z)
        if (s.blocks[z]) text += s.blocks[z];
    text += s.footer;
    return text;
}

bool embeddedText(const std::string& path, std::string& text, int max_z)
{
    if (!isEmbeddedPath(path)) return false;
    const std::string name = path.substr(std::strlen(kEmbeddedPrefix));
    for (const char* ext : { ".dat", ".turbomole" }) {
        const size_t n = std::strlen(ext);
        if (name.size() > n && name.compare(name.size() - n, n, ext) == 0) {
            const tables::EmbeddedBasisSet* s = findEmbedded(name.substr(0, name.size() - n), ext[1] == 'd' ? "data" : "turbomole");
            if (!s) return false;
            text = embeddedSetText(*s, max_z);
            return true;
        }
    }
    return false;
}

const char* embeddedElementBlock(const std::string& name, int Z)
{
    for (int i = 0; i < tables::kEmbeddedBasisSetCount; ++i) {
        const auto& s = tables::kEmbeddedBasisSets[i];
        if (name == s.name) return (Z >= 1 && Z <= s.max_z) ? s.blocks[Z - 1] : nullptr;
    }
    return nullptr;
}

std::vector<int> embeddedElements(const tables::EmbeddedBasisSet& s)
{
    std::vector<int> z;
    for (int i = 0; i < s.max_z; ++i)
        if (s.blocks[i]) z.push_back(i + 1);
    return z;
}

}  // namespace qmbasis
