/*
 * External potentials for SimpleMD — a configured, persistent bias.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated 2026
 */

#include "external_potentials.h"

#include <cmath>

namespace curcuma {

namespace {

    /// Resolve a selection expression against @p molecule.
    ///
    /// GetFragments() first, always: FragString2Indicies reads the fragment cache
    /// directly and every "Fn" matches nothing while that cache is empty. It fails
    /// silently, which is how it caught out the first caller that used it.
    std::vector<int> resolve(const Molecule& molecule, const std::string& expression)
    {
        Molecule copy = molecule;
        copy.GetFragments();
        return copy.FragString2Indicies(expression);
    }

    std::array<double, 3> centroid(const Geometry& geometry, const std::vector<int>& atoms)
    {
        std::array<double, 3> c { { 0.0, 0.0, 0.0 } };
        if (atoms.empty())
            return c;
        for (int index : atoms) {
            c[0] += geometry(index, 0);
            c[1] += geometry(index, 1);
            c[2] += geometry(index, 2);
        }
        const double n = static_cast<double>(atoms.size());
        return { { c[0] / n, c[1] / n, c[2] / n } };
    }

    bool readVector(const json& entry, const char* key, std::array<double, 3>& out)
    {
        if (!entry.contains(key) || !entry[key].is_array() || entry[key].size() != 3)
            return false;
        for (int i = 0; i < 3; ++i)
            out[i] = entry[key][i].get<double>();
        return true;
    }

} // namespace

std::vector<ExternalPotential> parseExternalPotentials(const json& list,
    const Molecule& molecule, std::string* error)
{
    const auto fail = [error](const std::string& message) {
        if (error)
            *error = message;
        return std::vector<ExternalPotential> {};
    };

    std::vector<ExternalPotential> potentials;
    if (list.is_null())
        return potentials;
    if (!list.is_array())
        return fail("external_potentials must be a list of objects");

    int number = 0;
    for (const json& entry : list) {
        ++number;
        const std::string where = "external_potentials[" + std::to_string(number - 1) + "]";
        if (!entry.is_object())
            return fail(where + " is not an object");

        ExternalPotential potential;
        potential.label = entry.value("label", "potential " + std::to_string(number));

        const std::string kind = entry.value("kind", std::string());
        if (kind == "constant_force")
            potential.kind = ExternalPotential::Kind::ConstantForce;
        else if (kind == "centroid_harmonic")
            potential.kind = ExternalPotential::Kind::CentroidHarmonic;
        else if (kind == "distance_harmonic")
            potential.kind = ExternalPotential::Kind::DistanceHarmonic;
        else {
            return fail(where + ": kind must be constant_force, centroid_harmonic or "
                                "distance_harmonic, not \"" + kind + "\"");
        }

        potential.selection = entry.value("atoms", std::string());
        if (potential.selection.empty())
            return fail(where + ": needs \"atoms\", a selection like \"F2\" or \"1:20\"");
        potential.atoms = resolve(molecule, potential.selection);
        if (potential.atoms.empty()) {
            return fail(where + ": selection \"" + potential.selection
                + "\" matched no atoms -- a bias on nothing is a configuration error");
        }

        switch (potential.kind) {
        case ExternalPotential::Kind::ConstantForce: {
            if (!readVector(entry, "direction", potential.direction))
                return fail(where + ": constant_force needs \"direction\" as [x, y, z]");
            const double norm = std::sqrt(potential.direction[0] * potential.direction[0]
                + potential.direction[1] * potential.direction[1]
                + potential.direction[2] * potential.direction[2]);
            if (norm < 1e-12)
                return fail(where + ": \"direction\" has zero length");
            for (int i = 0; i < 3; ++i)
                potential.direction[i] /= norm;
            potential.magnitude = entry.value("magnitude", 0.0);
            if (potential.magnitude == 0.0)
                return fail(where + ": constant_force needs a non-zero \"magnitude\" in Eh/Angstrom");
            break;
        }
        case ExternalPotential::Kind::CentroidHarmonic: {
            if (!readVector(entry, "target", potential.target))
                return fail(where + ": centroid_harmonic needs \"target\" as [x, y, z]");
            potential.k = entry.value("k", 0.0);
            if (potential.k <= 0.0)
                return fail(where + ": centroid_harmonic needs a positive \"k\" in Eh/Angstrom^2");
            break;
        }
        case ExternalPotential::Kind::DistanceHarmonic: {
            potential.selection_b = entry.value("atoms_b", std::string());
            if (potential.selection_b.empty())
                return fail(where + ": distance_harmonic needs a second selection \"atoms_b\"");
            potential.atoms_b = resolve(molecule, potential.selection_b);
            if (potential.atoms_b.empty()) {
                return fail(where + ": selection \"" + potential.selection_b
                    + "\" matched no atoms");
            }
            potential.k = entry.value("k", 0.0);
            if (potential.k <= 0.0)
                return fail(where + ": distance_harmonic needs a positive \"k\" in Eh/Angstrom^2");
            if (!entry.contains("r0"))
                return fail(where + ": distance_harmonic needs \"r0\" in Angstrom");
            potential.r0 = entry.value("r0", 0.0);
            break;
        }
        }

        potentials.push_back(potential);
    }

    if (error)
        error->clear();
    return potentials;
}

double applyExternalPotentials(std::vector<ExternalPotential>& potentials,
    const Geometry& geometry, const Geometry& previous, Geometry& gradient)
{
    if (potentials.empty())
        return 0.0;

    const bool trackWork = previous.rows() == geometry.rows()
        && previous.cols() == geometry.cols() && geometry.rows() > 0;
    double total = 0.0;

    for (ExternalPotential& potential : potentials) {
        // The force on each atom of this potential, so the work it does over the
        // step can be summed without evaluating anything twice.
        std::vector<std::array<double, 3>> forces(potential.atoms.size(), { { 0, 0, 0 } });
        std::vector<std::array<double, 3>> forcesB(potential.atoms_b.size(), { { 0, 0, 0 } });

        switch (potential.kind) {
        case ExternalPotential::Kind::ConstantForce: {
            // E = -sum_i F . r_i, so dE/dr_i = -F and the force is +F.
            double dot = 0.0;
            for (size_t a = 0; a < potential.atoms.size(); ++a) {
                const int i = potential.atoms[a];
                for (int c = 0; c < 3; ++c) {
                    const double f = potential.magnitude * potential.direction[c];
                    forces[a][c] = f;
                    gradient(i, c) -= f;
                    dot += f * geometry(i, c);
                }
            }
            // Origin-dependent by nature; the accumulated work below is the
            // quantity that means something.
            potential.energy = -dot;
            potential.value = 0.0;
            break;
        }
        case ExternalPotential::Kind::CentroidHarmonic: {
            const std::array<double, 3> c = centroid(geometry, potential.atoms);
            double d2 = 0.0;
            std::array<double, 3> delta { { 0, 0, 0 } };
            for (int j = 0; j < 3; ++j) {
                delta[j] = c[j] - potential.target[j];
                d2 += delta[j] * delta[j];
            }
            potential.value = std::sqrt(d2);
            potential.energy = 0.5 * potential.k * d2;
            const double n = static_cast<double>(potential.atoms.size());
            for (size_t a = 0; a < potential.atoms.size(); ++a) {
                const int i = potential.atoms[a];
                for (int j = 0; j < 3; ++j) {
                    const double g = potential.k * delta[j] / n;   // dE/dr_i
                    gradient(i, j) += g;
                    forces[a][j] = -g;
                }
            }
            break;
        }
        case ExternalPotential::Kind::DistanceHarmonic: {
            const std::array<double, 3> ca = centroid(geometry, potential.atoms);
            const std::array<double, 3> cb = centroid(geometry, potential.atoms_b);
            std::array<double, 3> delta { { 0, 0, 0 } };
            double d2 = 0.0;
            for (int j = 0; j < 3; ++j) {
                delta[j] = ca[j] - cb[j];
                d2 += delta[j] * delta[j];
            }
            const double d = std::sqrt(d2);
            potential.value = d;
            if (d < 1e-9) {
                // Coincident centroids: the direction is undefined, and a restraint
                // that cannot say which way to push has to do nothing rather than
                // divide by zero.
                potential.energy = 0.5 * potential.k * potential.r0 * potential.r0;
                break;
            }
            const double diff = d - potential.r0;
            potential.energy = 0.5 * potential.k * diff * diff;
            const double na = static_cast<double>(potential.atoms.size());
            const double nb = static_cast<double>(potential.atoms_b.size());
            for (size_t a = 0; a < potential.atoms.size(); ++a) {
                const int i = potential.atoms[a];
                for (int j = 0; j < 3; ++j) {
                    const double g = potential.k * diff * delta[j] / d / na;
                    gradient(i, j) += g;
                    forces[a][j] = -g;
                }
            }
            for (size_t b = 0; b < potential.atoms_b.size(); ++b) {
                const int i = potential.atoms_b[b];
                for (int j = 0; j < 3; ++j) {
                    const double g = -potential.k * diff * delta[j] / d / nb;
                    gradient(i, j) += g;
                    forcesB[b][j] = -g;
                }
            }
            break;
        }
        }

        if (trackWork) {
            double dw = 0.0;
            for (size_t a = 0; a < potential.atoms.size(); ++a) {
                const int i = potential.atoms[a];
                for (int j = 0; j < 3; ++j)
                    dw += forces[a][j] * (geometry(i, j) - previous(i, j));
            }
            for (size_t b = 0; b < potential.atoms_b.size(); ++b) {
                const int i = potential.atoms_b[b];
                for (int j = 0; j < 3; ++j)
                    dw += forcesB[b][j] * (geometry(i, j) - previous(i, j));
            }
            potential.work += dw;
        }

        total += potential.energy;
    }
    return total;
}

json describeExternalPotentials(const std::vector<ExternalPotential>& potentials)
{
    json list = json::array();
    for (const ExternalPotential& potential : potentials) {
        json entry;
        entry["label"] = potential.label;
        entry["atoms"] = potential.selection;
        entry["atom_count"] = static_cast<int>(potential.atoms.size());
        switch (potential.kind) {
        case ExternalPotential::Kind::ConstantForce:
            entry["kind"] = "constant_force";
            entry["direction"] = potential.direction;
            entry["magnitude"] = potential.magnitude;
            break;
        case ExternalPotential::Kind::CentroidHarmonic:
            entry["kind"] = "centroid_harmonic";
            entry["target"] = potential.target;
            entry["k"] = potential.k;
            entry["displacement"] = potential.value;
            break;
        case ExternalPotential::Kind::DistanceHarmonic:
            entry["kind"] = "distance_harmonic";
            entry["atoms_b"] = potential.selection_b;
            entry["atom_count_b"] = static_cast<int>(potential.atoms_b.size());
            entry["k"] = potential.k;
            entry["r0"] = potential.r0;
            entry["distance"] = potential.value;
            break;
        }
        entry["energy"] = potential.energy;
        entry["work"] = potential.work;
        list.push_back(entry);
    }
    return list;
}

} // namespace curcuma
