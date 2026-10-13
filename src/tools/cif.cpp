/*
 * CIF reader and supercell builder.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated 2026
 */

#include "cif.h"

#include "src/core/elements.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cctype>
#include <cmath>
#include <fstream>
#include <locale>
#include <map>
#include <sstream>

namespace curcuma {

namespace {

    std::string trim(const std::string& text)
    {
        const auto first = text.find_first_not_of(" \t\r\n");
        if (first == std::string::npos)
            return std::string();
        const auto last = text.find_last_not_of(" \t\r\n");
        return text.substr(first, last - first + 1);
    }

    /// A number as CIF writes it, always with a '.' decimal point. std::from_chars,
    /// not std::stod: stod follows the process locale, and inside a GUI that has
    /// adopted de_DE it read "0.5957" as 0 -- every coordinate truncated, Cl on
    /// top of Na. from_chars ignores the locale and does not throw. A leading '+'
    /// is accepted (from_chars would not); trailing text ends the number, as with
    /// stod. Claude Generated (Sep 2026).
    ///
    /// Floating-point from_chars is missing from older libc++ (Apple Clang), which
    /// does not define __cpp_lib_to_chars then; there a stream fixed to the classic
    /// "C" locale does the same job, also independent of the process locale.
    bool parseNumber(const std::string& text, double& value)
    {
        const char* begin = text.data();
        const char* end = begin + text.size();
        while (begin < end && std::isspace(static_cast<unsigned char>(*begin)))
            ++begin;
        if (begin < end && *begin == '+')
            ++begin;
#if defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L && !defined(CURCUMA_CIF_STREAM_NUMBERS)
        const auto [ptr, ec] = std::from_chars(begin, end, value);
        return ec == std::errc() && ptr != begin;
#else
        std::istringstream stream(std::string(begin, end));
        stream.imbue(std::locale::classic());
        double parsed = 0.0;
        if (!(stream >> parsed))
            return false;
        value = parsed;
        return true;
#endif
    }

    /// CIF numbers carry their uncertainty in brackets: "5.4307(2)". Drop it.
    double toNumber(const std::string& token, bool* ok = nullptr)
    {
        std::string clean;
        for (char c : token) {
            if (c == '(')
                break;
            clean += c;
        }
        double value = 0.0;
        const bool parsed = parseNumber(clean, value);
        if (ok)
            *ok = parsed;
        return parsed ? value : 0.0;
    }

    /// Split a CIF data line into tokens, honouring single and double quotes.
    std::vector<std::string> tokenise(const std::string& line)
    {
        std::vector<std::string> tokens;
        std::string current;
        char quote = 0;
        for (char c : line) {
            if (quote) {
                if (c == quote) {
                    quote = 0;
                    tokens.push_back(current);
                    current.clear();
                } else {
                    current += c;
                }
                continue;
            }
            if (c == '\'' || c == '"') {
                quote = c;
                continue;
            }
            if (std::isspace(static_cast<unsigned char>(c))) {
                if (!current.empty()) {
                    tokens.push_back(current);
                    current.clear();
                }
                continue;
            }
            current += c;
        }
        if (!current.empty())
            tokens.push_back(current);
        return tokens;
    }

    /// "Fe3+" -> "Fe", "O1" -> "O", "C12A" -> "C". Leading letters only, and the
    /// second one only when it is lower case, which is how element symbols are
    /// spelled.
    std::string elementOf(const std::string& token)
    {
        std::string symbol;
        for (char c : token) {
            if (std::isalpha(static_cast<unsigned char>(c))) {
                if (symbol.empty())
                    symbol += char(std::toupper(static_cast<unsigned char>(c)));
                else if (symbol.size() == 1 && std::islower(static_cast<unsigned char>(c)))
                    symbol += c;
                else
                    break;
            } else {
                break;
            }
        }
        return symbol;
    }

    /// Lattice vectors from the six cell parameters, in the standard setting:
    /// a along x, b in the xy plane, c completing it.
    Eigen::Matrix3d latticeFrom(double a, double b, double c,
        double alpha, double beta, double gamma)
    {
        const double d2r = pi / 180.0;
        const double ca = std::cos(alpha * d2r);
        const double cb = std::cos(beta * d2r);
        const double cg = std::cos(gamma * d2r);
        const double sg = std::sin(gamma * d2r);

        Eigen::Matrix3d lattice = Eigen::Matrix3d::Zero();
        lattice(0, 0) = a;
        lattice(0, 1) = b * cg;
        lattice(1, 1) = b * sg;
        lattice(0, 2) = c * cb;
        lattice(1, 2) = c * (ca - cb * cg) / sg;
        const double zz = c * c - lattice(0, 2) * lattice(0, 2) - lattice(1, 2) * lattice(1, 2);
        lattice(2, 2) = zz > 0.0 ? std::sqrt(zz) : 0.0;
        return lattice;
    }

    /// One component of a symmetry operation, e.g. "-x+1/2" or "y".
    /// Returns the three coefficients and the translation.
    bool parseComponent(const std::string& text, double coefficient[3], double& translation)
    {
        coefficient[0] = coefficient[1] = coefficient[2] = 0.0;
        translation = 0.0;

        std::string s;
        for (char c : text) {
            if (!std::isspace(static_cast<unsigned char>(c)))
                s += char(std::tolower(static_cast<unsigned char>(c)));
        }
        if (s.empty())
            return false;

        size_t i = 0;
        while (i < s.size()) {
            int sign = 1;
            if (s[i] == '+') {
                ++i;
            } else if (s[i] == '-') {
                sign = -1;
                ++i;
            }
            if (i >= s.size())
                return false;

            // A bare axis, or a number that may be a fraction.
            if (s[i] == 'x' || s[i] == 'y' || s[i] == 'z') {
                coefficient[s[i] - 'x'] += sign;
                ++i;
                continue;
            }
            std::string number;
            while (i < s.size() && (std::isdigit(static_cast<unsigned char>(s[i]))
                       || s[i] == '.' || s[i] == '/')) {
                number += s[i];
                ++i;
            }
            if (number.empty())
                return false;
            double value = 0.0;
            const auto slash = number.find('/');
            // A malformed token ("." or "/2") is an unreadable operation, which the
            // caller counts and reports; stod threw here, uncaught.
            if (slash != std::string::npos) {
                double numerator = 0.0;
                double denominator = 0.0;
                if (!parseNumber(number.substr(0, slash), numerator)
                    || !parseNumber(number.substr(slash + 1), denominator)
                    || denominator == 0.0)
                    return false;
                value = numerator / denominator;
            } else if (!parseNumber(number, value)) {
                return false;
            }
            // "1/2x" is a coefficient, "1/2" on its own is a translation.
            if (i < s.size() && (s[i] == 'x' || s[i] == 'y' || s[i] == 'z')) {
                coefficient[s[i] - 'x'] += sign * value;
                ++i;
            } else {
                translation += sign * value;
            }
        }
        return true;
    }

    using SymmetryOperation = CifSymmetryOperation;

    bool parseOperation(const std::string& text, SymmetryOperation& out)
    {
        std::vector<std::string> parts;
        std::string current;
        for (char c : text) {
            if (c == ',') {
                parts.push_back(current);
                current.clear();
            } else {
                current += c;
            }
        }
        parts.push_back(current);
        if (parts.size() != 3)
            return false;
        for (int row = 0; row < 3; ++row) {
            double coefficient[3];
            double translation = 0.0;
            if (!parseComponent(parts[size_t(row)], coefficient, translation))
                return false;
            for (int col = 0; col < 3; ++col)
                out.matrix[row][col] = coefficient[col];
            out.translation[row] = translation;
        }
        return true;
    }

    double wrap(double value)
    {
        double v = std::fmod(value, 1.0);
        if (v < 0.0)
            v += 1.0;
        // 0.9999999 is 0, not a separate site.
        if (v > 1.0 - 1e-6)
            v = 0.0;
        return v;
    }

    /// Claude Generated (Sep 2026) - Reassemble molecules cut by the cell faces.
    /// Two atoms are bonded when their minimum-image distance is below 1.15 x
    /// the sum of their covalent radii (curcuma's table, as GetFragments uses).
    /// A breadth-first walk then gives each atom the image that sits next to the
    /// atom it was reached from, and every finished molecule is shifted by a
    /// lattice vector so its centroid lies in [0,1)^3. O(N^2) in the cell atoms.
    template <typename Atoms>
    void completeMolecules(Atoms& atoms, const Eigen::Matrix3d& lattice)
    {
        const size_t n = atoms.size();
        std::vector<double> radius(n, 1.5);
        for (size_t i = 0; i < n; ++i) {
            const int z = Elements::String2Element(atoms[i].element);
            if (z > 0 && z < int(Elements::CovalentRadius.size()))
                radius[i] = Elements::CovalentRadius[size_t(z)];
        }
        const auto minimumImage = [](Eigen::Vector3d d) {
            for (int k = 0; k < 3; ++k)
                d(k) -= std::round(d(k));
            return d;
        };
        std::vector<std::vector<size_t>> neighbours(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                const Eigen::Vector3d d = minimumImage(atoms[j].fractional - atoms[i].fractional);
                if ((lattice * d).norm() < 1.15 * (radius[i] + radius[j])) {
                    neighbours[i].push_back(j);
                    neighbours[j].push_back(i);
                }
            }
        }
        std::vector<bool> placed(n, false);
        for (size_t start = 0; start < n; ++start) {
            if (placed[start])
                continue;
            std::vector<size_t> molecule { start };
            placed[start] = true;
            for (size_t head = 0; head < molecule.size(); ++head) {
                const size_t current = molecule[head];
                for (size_t next : neighbours[current]) {
                    if (placed[next])
                        continue;
                    atoms[next].fractional = atoms[current].fractional
                        + minimumImage(atoms[next].fractional - atoms[current].fractional);
                    placed[next] = true;
                    molecule.push_back(next);
                }
            }
            Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
            for (size_t index : molecule)
                centroid += atoms[index].fractional;
            centroid /= double(molecule.size());
            Eigen::Vector3d shift;
            for (int k = 0; k < 3; ++k)
                shift(k) = -std::floor(centroid(k));
            for (size_t index : molecule)
                atoms[index].fractional += shift;
        }
    }

} // namespace

int CifData::majorDisorderGroup() const
{
    int best = 0;
    double bestOccupancy = -1.0;
    for (const CifDisorderGroup& group : disorder_groups) {
        if (group.mean_occupancy > bestOccupancy) {
            bestOccupancy = group.mean_occupancy;
            best = group.group;
        }
    }
    return best;
}

CifData ReadCifData(const std::string& filename)
{
    CifData result;
    std::ifstream file(filename);
    if (!file.is_open()) {
        result.error = "cannot open " + filename;
        return result;
    }

    std::vector<CifSite>& sites = result.sites;
    std::map<std::string, std::array<double, 6>> aniso;   // label -> U11..U23
    std::vector<SymmetryOperation>& operations = result.operations;
    bool cartesian = false;
    int unreadableOperations = 0;
    int dataBlocks = 0;
    bool skippedBlocks = false;

    std::string line;
    while (std::getline(file, line)) {
        const std::string stripped = trim(line);
        if (stripped.empty() || stripped[0] == '#')
            continue;

        // Claude Generated (Sep 2026) - One structure per read. Several data
        // blocks used to be read as one: every block's atoms and operations were
        // collected, under the cell of the last. The first block that has atom
        // sites is the structure; a block before it without sites (a "global"
        // block of a multi-structure file) must not leave its cell or operations
        // behind.
        if (stripped.rfind("data_", 0) == 0) {
            ++dataBlocks;
            if (!sites.empty()) {
                skippedBlocks = true;
                break;
            }
            result.cell = CifCell();
            operations.clear();
            unreadableOperations = 0;
            continue;
        }

        if (stripped.rfind("_cell_length_a", 0) == 0
            || stripped.rfind("_cell_length_b", 0) == 0
            || stripped.rfind("_cell_length_c", 0) == 0
            || stripped.rfind("_cell_angle_", 0) == 0) {
            const std::vector<std::string> tokens = tokenise(stripped);
            if (tokens.size() < 2)
                continue;
            bool ok = false;
            const double value = toNumber(tokens[1], &ok);
            if (!ok)
                continue;
            if (tokens[0] == "_cell_length_a") result.cell.a = value;
            else if (tokens[0] == "_cell_length_b") result.cell.b = value;
            else if (tokens[0] == "_cell_length_c") result.cell.c = value;
            else if (tokens[0] == "_cell_angle_alpha") result.cell.alpha = value;
            else if (tokens[0] == "_cell_angle_beta") result.cell.beta = value;
            else if (tokens[0] == "_cell_angle_gamma") result.cell.gamma = value;
            continue;
        }

        // Claude Generated (Sep 2026) - What the file says the structure is, so a
        // reader can check the atoms against it (Z x formula).
        {
            std::vector<std::string> tokens = tokenise(stripped);
            // A text field may sit on the following lines between ';' lines.
            if (tokens.size() == 1 && (tokens[0] == "_chemical_formula_moiety"
                                          || tokens[0] == "_chemical_formula_sum")) {
                std::streampos back = file.tellg();
                std::string next;
                while (std::getline(file, next) && trim(next).empty())
                    back = file.tellg();
                if (!trim(next).empty() && trim(next)[0] == ';') {
                    std::string text = trim(trim(next).substr(1));
                    while (std::getline(file, next) && trim(next).rfind(";", 0) != 0)
                        text += (text.empty() ? "" : " ") + trim(next);
                    tokens.push_back(text);
                } else {
                    file.clear();
                    file.seekg(back);
                }
            }
            if (tokens.size() >= 2) {
                std::string rest;
                for (size_t k = 1; k < tokens.size(); ++k)
                    rest += (k > 1 ? " " : "") + tokens[k];
                if (tokens[0] == "_space_group_name_H-M_alt"
                    || (tokens[0] == "_symmetry_space_group_name_H-M" && result.space_group.empty())) {
                    result.space_group = rest;
                    continue;
                }
                if (tokens[0] == "_space_group_IT_number" || tokens[0] == "_symmetry_Int_Tables_number") {
                    result.space_group_number = int(std::lround(toNumber(tokens[1])));
                    continue;
                }
                if (tokens[0] == "_cell_formula_units_Z") {
                    result.formula_units_z = int(std::lround(toNumber(tokens[1])));
                    continue;
                }
                if (tokens[0] == "_chemical_formula_sum") {
                    result.formula_sum = rest;
                    continue;
                }
                if (tokens[0] == "_chemical_formula_moiety") {
                    result.formula_moiety = rest;
                    continue;
                }
            }
        }

        if (stripped != "loop_")
            continue;

        // Collect the loop's tags, then its rows.
        std::vector<std::string> tags;
        std::streampos mark = file.tellg();
        while (std::getline(file, line)) {
            const std::string tag = trim(line);
            if (tag.empty() || tag[0] == '#')
                continue;
            if (tag[0] != '_') {
                file.seekg(mark);
                break;
            }
            tags.push_back(tokenise(tag)[0]);
            mark = file.tellg();
        }
        const auto indexOf = [&tags](const std::string& name) {
            const auto it = std::find(tags.begin(), tags.end(), name);
            return it == tags.end() ? -1 : int(it - tags.begin());
        };

        const int symmetryColumn = indexOf("_symmetry_equiv_pos_as_xyz") >= 0
            ? indexOf("_symmetry_equiv_pos_as_xyz")
            : indexOf("_space_group_symop_operation_xyz");
        const int typeColumn = indexOf("_atom_site_type_symbol");
        const int labelColumn = indexOf("_atom_site_label");
        const int occupancyColumn = indexOf("_atom_site_occupancy");
        const int assemblyColumn = indexOf("_atom_site_disorder_assembly");
        const int groupColumn = indexOf("_atom_site_disorder_group");
        // Displacement parameters: U or B (B = 8 pi^2 U), isotropic in the site
        // loop, anisotropic in a loop of their own keyed by label.
        const int uIsoColumn = indexOf("_atom_site_U_iso_or_equiv");
        const int bIsoColumn = indexOf("_atom_site_B_iso_or_equiv");
        const int anisoLabelColumn = indexOf("_atom_site_aniso_label");
        const char* const uTags[6] = { "_atom_site_aniso_U_11", "_atom_site_aniso_U_22",
            "_atom_site_aniso_U_33", "_atom_site_aniso_U_12", "_atom_site_aniso_U_13",
            "_atom_site_aniso_U_23" };
        const char* const bTags[6] = { "_atom_site_aniso_B_11", "_atom_site_aniso_B_22",
            "_atom_site_aniso_B_33", "_atom_site_aniso_B_12", "_atom_site_aniso_B_13",
            "_atom_site_aniso_B_23" };
        int anisoColumns[6];
        bool anisoIsB = false;
        for (int k = 0; k < 6; ++k)
            anisoColumns[k] = indexOf(uTags[k]);
        if (anisoColumns[0] < 0) {
            for (int k = 0; k < 6; ++k)
                anisoColumns[k] = indexOf(bTags[k]);
            anisoIsB = anisoColumns[0] >= 0;
        }
        int xColumn = indexOf("_atom_site_fract_x");
        int yColumn = indexOf("_atom_site_fract_y");
        int zColumn = indexOf("_atom_site_fract_z");
        if (xColumn < 0) {
            xColumn = indexOf("_atom_site_Cartn_x");
            yColumn = indexOf("_atom_site_Cartn_y");
            zColumn = indexOf("_atom_site_Cartn_z");
            if (xColumn >= 0)
                cartesian = true;
        }

        while (std::getline(file, line)) {
            const std::string row = trim(line);
            if (row.empty() || row[0] == '#')
                continue;
            if (row[0] == '_' || row == "loop_" || row.rfind("data_", 0) == 0) {
                // Back up: this line belongs to whatever comes next.
                file.seekg(mark);
                break;
            }
            mark = file.tellg();
            const std::vector<std::string> tokens = tokenise(row);
            const auto column = [&tokens](int index) -> std::string {
                return index >= 0 && index < int(tokens.size()) ? tokens[size_t(index)] : std::string();
            };

            if (anisoLabelColumn >= 0 && anisoColumns[0] >= 0 && xColumn < 0) {
                std::array<double, 6> u {};
                bool complete = true;
                for (int k = 0; k < 6; ++k) {
                    bool ok = false;
                    u[size_t(k)] = toNumber(column(anisoColumns[k]), &ok);
                    complete = complete && ok;
                }
                if (complete) {
                    if (anisoIsB) {
                        for (double& value : u)
                            value /= 8.0 * pi * pi;
                    }
                    aniso[column(anisoLabelColumn)] = u;
                }
                continue;
            }
            if (symmetryColumn >= 0 && int(tokens.size()) > symmetryColumn) {
                SymmetryOperation operation;
                if (parseOperation(tokens[size_t(symmetryColumn)], operation))
                    operations.push_back(operation);
                else
                    ++unreadableOperations;
                continue;
            }
            if (xColumn >= 0 && int(tokens.size()) > std::max({ xColumn, yColumn, zColumn })) {
                CifSite site;
                site.label = column(labelColumn);
                site.element = elementOf(typeColumn >= 0 ? column(typeColumn) : site.label);
                if (site.element.empty())
                    continue;
                site.coordinate = Eigen::Vector3d(toNumber(column(xColumn)),
                    toNumber(column(yColumn)), toNumber(column(zColumn)));
                // '.' and '?' mean "not given": fully occupied, not disordered.
                bool ok = false;
                const double occupancy = toNumber(column(occupancyColumn), &ok);
                if (ok)
                    site.occupancy = occupancy;
                const std::string assembly = column(assemblyColumn);
                if (assembly != "." && assembly != "?")
                    site.disorder_assembly = assembly;
                const double group = toNumber(column(groupColumn), &ok);
                if (ok)
                    site.disorder_group = int(std::lround(group));
                const double uIso = toNumber(column(uIsoColumn), &ok);
                if (ok) {
                    site.has_iso = true;
                    site.u_iso = uIso;
                } else {
                    const double bIso = toNumber(column(bIsoColumn), &ok);
                    if (ok) {
                        site.has_iso = true;
                        site.u_iso = bIso / (8.0 * pi * pi);
                    }
                }
                sites.push_back(site);
            }
        }
    }

    if (sites.empty()) {
        result.error = "no atom sites found in " + filename;
        return result;
    }
    result.fractional = !cartesian;
    // The aniso loop may come before or after the sites; match by label.
    for (CifSite& site : sites) {
        const auto it = aniso.find(site.label);
        if (it == aniso.end())
            continue;
        site.has_aniso = true;
        for (int k = 0; k < 6; ++k)
            site.u_aniso[k] = it->second[size_t(k)];
    }
    if (skippedBlocks || dataBlocks > 1) {
        result.notes.push_back("the file holds more than one data block; only the first with "
                               "atom sites was read");
    }

    result.cell.valid = result.cell.a > 0.0 && result.cell.b > 0.0 && result.cell.c > 0.0;
    if (result.cell.valid) {
        result.cell.lattice = latticeFrom(result.cell.a, result.cell.b, result.cell.c,
            result.cell.alpha, result.cell.beta, result.cell.gamma);
    } else if (!cartesian) {
        result.error = "fractional coordinates without a usable cell";
        return result;
    } else {
        result.notes.push_back("no unit cell in the file; the coordinates are taken as they are "
                               "and no supercell can be built from them");
    }

    if (unreadableOperations > 0) {
        result.notes.push_back(std::to_string(unreadableOperations)
            + " symmetry operation(s) could not be parsed and were skipped -- the structure may "
              "be incomplete");
    }
    if (operations.empty()) {
        SymmetryOperation identity;
        for (int i = 0; i < 3; ++i)
            identity.matrix[i][i] = 1.0;
        operations.push_back(identity);
        result.notes.push_back("no symmetry operations in the file; the atom sites are taken as "
                               "the whole cell (P1)");
    }

    // Disorder groups over the whole file, for the caller to choose between.
    std::map<int, std::pair<int, double>> groups;   // group -> (sites, occupancy sum)
    for (const CifSite& site : sites) {
        if (site.disorder_group == 0)
            continue;
        auto& entry = groups[std::abs(site.disorder_group)];
        ++entry.first;
        entry.second += site.occupancy;
    }
    int disorderedSites = 0;
    for (const auto& [group, entry] : groups) {
        CifDisorderGroup summary;
        summary.group = group;
        summary.sites = entry.first;
        summary.mean_occupancy = entry.second / entry.first;
        result.disorder_groups.push_back(summary);
        disorderedSites += entry.first;
    }
    if (!groups.empty()) {
        result.notes.push_back(std::to_string(disorderedSites) + " of "
            + std::to_string(sites.size()) + " sites are disordered, in "
            + std::to_string(groups.size()) + " group(s)");
    }
    return result;
}

Molecule BuildCif(const CifData& data, const CifBuildOptions& options,
    std::vector<Eigen::Matrix3d>* displacements)
{
    Molecule molecule;
    if (displacements)
        displacements->clear();
    if (!data.ok())
        return molecule;

    // U* = N U N, on the fractional axes; U_cart = M U* M^T. N holds the
    // reciprocal lengths a*, b*, c* -- the row norms of M^-1.
    const Eigen::Matrix3d& cellMatrix = data.cell.lattice;
    Eigen::Matrix3d reciprocal = Eigen::Matrix3d::Identity();
    if (data.cell.valid) {
        const Eigen::Matrix3d inverse = cellMatrix.inverse();
        for (int i = 0; i < 3; ++i)
            reciprocal(i, i) = inverse.row(i).norm();
    }
    const auto fractionalU = [&reciprocal](const CifSite& site) {
        const double* u = site.u_aniso;
        Eigen::Matrix3d U;
        U << u[0], u[3], u[4],
             u[3], u[1], u[5],
             u[4], u[5], u[2];
        return Eigen::Matrix3d(reciprocal * U * reciprocal);
    };
    // The tensor of one atom: the site's, turned by @p rotation (fractional).
    const auto displacement = [&](const CifSite& site, const Eigen::Matrix3d& rotation) {
        if (site.has_aniso && data.cell.valid && data.fractional) {
            const Eigen::Matrix3d turned = rotation * fractionalU(site) * rotation.transpose();
            return Eigen::Matrix3d(cellMatrix * turned * cellMatrix.transpose());
        }
        if (site.has_iso)
            return Eigen::Matrix3d(site.u_iso * Eigen::Matrix3d::Identity());
        return Eigen::Matrix3d(Eigen::Matrix3d::Zero());
    };

    const Eigen::Matrix3d& lattice = data.cell.lattice;
    const bool cartesian = !data.fractional;
    const double tolerance = 0.1;   // Angstrom
    const auto selected = [&options](const CifSite& site) {
        return !options.select_disorder_group || site.disorder_group == 0
            || std::abs(site.disorder_group) == std::abs(options.disorder_group);
    };
    const auto add = [&molecule](const std::string& element, const Eigen::Vector3d& position) {
        molecule.addPair({ Elements::String2Element(element),
            Position(position(0), position(1), position(2)) });
    };

    // Unit-cell atoms are collected first: completing molecules moves them after
    // all images exist. Atom order is kept either way.
    struct CellAtom {
        std::string element;
        Eigen::Vector3d fractional;
        Eigen::Matrix3d u;
    };
    std::vector<CellAtom> cellAtoms;

    const Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
    for (const CifSite& site : data.sites) {
        if (!selected(site))
            continue;
        // As written: no cell, or the asymmetric unit, which keeps the file's
        // coordinates unwrapped so a molecule stays whole.
        if (cartesian) {
            add(site.element, site.coordinate);
            if (displacements)
                displacements->push_back(displacement(site, identity));
            continue;
        }
        if (options.content == CifContent::AsymmetricUnit) {
            add(site.element, lattice * site.coordinate);
            if (displacements)
                displacements->push_back(displacement(site, identity));
            continue;
        }

        // Every image of this site; images of the same site that coincide
        // (periodically) are one atom. Only this site's images are compared:
        // merging across sites would swallow a disorder alternative that happens
        // to lie close to its partner.
        std::vector<Eigen::Vector3d> images;
        for (const SymmetryOperation& operation : data.operations) {
            Eigen::Vector3d fractional;
            for (int row = 0; row < 3; ++row) {
                fractional(row) = operation.matrix[row][0] * site.coordinate(0)
                    + operation.matrix[row][1] * site.coordinate(1)
                    + operation.matrix[row][2] * site.coordinate(2) + operation.translation[row];
            }
            for (int i = 0; i < 3; ++i)
                fractional(i) = wrap(fractional(i));
            bool duplicate = false;
            for (const Eigen::Vector3d& image : images) {
                Eigen::Vector3d delta = fractional - image;
                for (int i = 0; i < 3; ++i)
                    delta(i) -= std::round(delta(i));   // minimum image: 0.0 and 0.99999 meet
                if ((lattice * delta).norm() < tolerance) {
                    duplicate = true;
                    break;
                }
            }
            if (duplicate)
                continue;
            images.push_back(fractional);
            Eigen::Matrix3d rotation;
            for (int row = 0; row < 3; ++row)
                for (int col = 0; col < 3; ++col)
                    rotation(row, col) = operation.matrix[row][col];
            cellAtoms.push_back({ site.element, fractional, displacement(site, rotation) });
        }
    }

    if (options.complete_molecules && !cellAtoms.empty())
        completeMolecules(cellAtoms, lattice);
    for (const CellAtom& atom : cellAtoms) {
        add(atom.element, lattice * atom.fractional);
        if (displacements)
            displacements->push_back(atom.u);
    }
    if (data.cell.valid)
        molecule.setUnitCell(lattice, true);
    return molecule;
}

CifResult ReadCif(const std::string& filename)
{
    const CifData data = ReadCifData(filename);
    CifResult result;
    result.cell = data.cell;
    result.asymmetric_atoms = int(data.sites.size());
    result.symmetry_operations = int(data.operations.size());
    result.fractional = data.fractional;
    result.disorder_groups = data.disorder_groups;
    result.error = data.error;
    result.notes = data.notes;
    if (!data.ok())
        return result;
    result.molecule = BuildCif(data, CifBuildOptions());
    if (!data.disorder_groups.empty()) {
        result.notes.push_back("all disorder alternatives are included, so they overlap; "
                               "choose one group with ReadCifData() and BuildCif()");
    }
    return result;
}

Molecule Supercell(const Molecule& molecule, int na, int nb, int nc, std::string* error)
{
    const auto fail = [error, &molecule](const std::string& message) {
        if (error)
            *error = message;
        return molecule;
    };
    if (na < 1 || nb < 1 || nc < 1)
        return fail("a supercell needs at least 1 x 1 x 1");

    // hasPBC(), not the determinant: an unset cell is the IDENTITY matrix, whose
    // determinant is 1, so a determinant test would have replicated a molecule
    // along three imaginary 1 Angstrom axes and produced a pile of overlapping
    // atoms instead of an error.
    const Eigen::Matrix3d lattice = molecule.getUnitCell();
    if (!molecule.hasPBC() || lattice.determinant() <= 0.0)
        return fail("this structure carries no unit cell, so there is nothing to replicate");
    if (na == 1 && nb == 1 && nc == 1) {
        if (error)
            error->clear();
        return molecule;
    }

    Molecule out;
    for (int ia = 0; ia < na; ++ia) {
        for (int ib = 0; ib < nb; ++ib) {
            for (int ic = 0; ic < nc; ++ic) {
                const Eigen::Vector3d shift = lattice * Eigen::Vector3d(ia, ib, ic);
                for (int i = 0; i < molecule.AtomCount(); ++i) {
                    const std::pair<int, Position> atom = molecule.Atom(i);
                    out.addPair({ atom.first,
                        Position(atom.second(0) + shift(0), atom.second(1) + shift(1),
                            atom.second(2) + shift(2)) });
                }
            }
        }
    }
    Eigen::Matrix3d enlarged = lattice;
    enlarged.col(0) *= na;
    enlarged.col(1) *= nb;
    enlarged.col(2) *= nc;
    out.setUnitCell(enlarged, true);
    if (error)
        error->clear();
    return out;
}

} // namespace curcuma
