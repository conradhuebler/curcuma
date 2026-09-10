/*
 * CIF reader and supercell builder.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated 2026
 */

#include "cif.h"

#include "src/core/elements.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
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

    /// CIF numbers carry their uncertainty in brackets: "5.4307(2)". Drop it.
    double toNumber(const std::string& token, bool* ok = nullptr)
    {
        std::string clean;
        for (char c : token) {
            if (c == '(')
                break;
            clean += c;
        }
        try {
            const double value = std::stod(clean);
            if (ok)
                *ok = true;
            return value;
        } catch (...) {
            if (ok)
                *ok = false;
            return 0.0;
        }
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
            if (slash != std::string::npos) {
                const double numerator = std::stod(number.substr(0, slash));
                const double denominator = std::stod(number.substr(slash + 1));
                if (denominator == 0.0)
                    return false;
                value = numerator / denominator;
            } else {
                value = std::stod(number);
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

    struct SymmetryOperation {
        double matrix[3][3] {};
        double translation[3] {};
    };

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

} // namespace

CifResult ReadCif(const std::string& filename)
{
    CifResult result;
    std::ifstream file(filename);
    if (!file.is_open()) {
        result.error = "cannot open " + filename;
        return result;
    }

    struct Site {
        std::string element;
        double x = 0.0, y = 0.0, z = 0.0;
    };
    std::vector<Site> sites;
    std::vector<SymmetryOperation> operations;
    bool cartesian = false;
    int unreadableOperations = 0;

    std::string line;
    while (std::getline(file, line)) {
        const std::string stripped = trim(line);
        if (stripped.empty() || stripped[0] == '#')
            continue;

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

            if (symmetryColumn >= 0 && int(tokens.size()) > symmetryColumn) {
                SymmetryOperation operation;
                if (parseOperation(tokens[size_t(symmetryColumn)], operation))
                    operations.push_back(operation);
                else
                    ++unreadableOperations;
                continue;
            }
            if (xColumn >= 0 && int(tokens.size()) > std::max({ xColumn, yColumn, zColumn })) {
                Site site;
                if (typeColumn >= 0 && int(tokens.size()) > typeColumn)
                    site.element = elementOf(tokens[size_t(typeColumn)]);
                else if (labelColumn >= 0 && int(tokens.size()) > labelColumn)
                    site.element = elementOf(tokens[size_t(labelColumn)]);
                if (site.element.empty())
                    continue;
                site.x = toNumber(tokens[size_t(xColumn)]);
                site.y = toNumber(tokens[size_t(yColumn)]);
                site.z = toNumber(tokens[size_t(zColumn)]);
                sites.push_back(site);
            }
        }
    }

    if (sites.empty()) {
        result.error = "no atom sites found in " + filename;
        return result;
    }
    result.asymmetric_atoms = int(sites.size());
    result.fractional = !cartesian;

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
    result.symmetry_operations = int(operations.size());

    // Apply the operations and drop the duplicates they produce at special
    // positions. Without this the file's asymmetric unit would be mistaken for the
    // whole structure.
    std::vector<std::pair<std::string, Eigen::Vector3d>> atoms;
    const double tolerance = 0.1;   // Angstrom
    for (const Site& site : sites) {
        for (const SymmetryOperation& operation : operations) {
            Eigen::Vector3d fractional;
            for (int row = 0; row < 3; ++row) {
                fractional(row) = operation.matrix[row][0] * site.x
                    + operation.matrix[row][1] * site.y
                    + operation.matrix[row][2] * site.z + operation.translation[row];
            }
            Eigen::Vector3d position;
            if (cartesian) {
                position = Eigen::Vector3d(site.x, site.y, site.z);
            } else {
                for (int i = 0; i < 3; ++i)
                    fractional(i) = wrap(fractional(i));
                position = result.cell.lattice * fractional;
            }
            bool duplicate = false;
            for (const auto& existing : atoms) {
                if (existing.first == site.element
                    && (existing.second - position).norm() < tolerance) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate)
                atoms.emplace_back(site.element, position);
            if (cartesian)
                break;   // no cell, so the operations have nothing to act on
        }
    }

    for (const auto& atom : atoms) {
        const int z = Elements::String2Element(atom.first);
        result.molecule.addPair({ z, Position(atom.second(0), atom.second(1), atom.second(2)) });
    }
    if (result.cell.valid)
        result.molecule.setUnitCell(result.cell.lattice, true);
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
