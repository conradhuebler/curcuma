/*
 * <Native xTB large-system fragmentation driver (C1) — header>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (C1, June 2026): a shared fragmentation scaffold that scales
 * the native GFN1/GFN2 SCF to large systems by exploiting locality. The dense
 * O(N^3) eigensolve in curcuma::xtb::XTB is the wall above ~1000 atoms; this
 * driver runs the *existing, validated* dense XTB on sub-systems and assembles
 * the result, trading a controlled accuracy loss for tractable large-N runs.
 *
 * Three opt-in modes (default `none` = the normal dense path, untouched):
 *   - Fragments (C1c): disconnected-fragment SCF. Partition by bond
 *     connectivity, one dense SCF per fragment, sum energies + block-diagonal
 *     gradient. Exact for non-interacting clusters; neglects all inter-fragment
 *     coupling. A connected molecule yields one fragment -> falls back to dense.
 *   - DC (C1a): divide-and-conquer overlapping fragments (general large-system
 *     method, connected systems). [implemented in a later step]
 *   - Sparse (C1b): sparse build + non-orthogonal density purification, the
 *     O(N) path for gapped systems at 0 K. [implemented in a later step]
 *
 * The per-fragment engine is curcuma::xtb::XTB reused verbatim — all the
 * validated SCF / gradient / dispersion code comes for free.
 */

#pragma once

#include "src/core/global.h"      // Mol, Matrix, Vector, json
#include "src/core/molecule.h"    // Molecule::GetFragments (connectivity)
#include "xtb_native.h"           // curcuma::xtb::XTB, MethodType

#include <string>
#include <vector>

namespace curcuma::xtb {

/* ------------------------------------------------------------------------- *
 *  C1 fragmentation mode selector
 * ------------------------------------------------------------------------- */
enum class C1Mode {
    None,       // dense reference (default) — driver not engaged
    Fragments,  // C1c — disconnected-fragment SCF
    DC,         // C1a — divide-and-conquer overlapping fragments
    Sparse,     // C1b — sparse build + non-orthogonal purification
};

/// Parse a user string to a C1Mode. Unknown / empty -> None (dense).
/// Aliases: frag->fragments, divide-conquer/divide_and_conquer->dc.
inline C1Mode parseC1Mode(const std::string& s)
{
    if (s == "fragments" || s == "frag" || s == "fragment") return C1Mode::Fragments;
    if (s == "dc" || s == "divide-conquer" || s == "divide_and_conquer") return C1Mode::DC;
    if (s == "sparse") return C1Mode::Sparse;
    return C1Mode::None;
}

inline const char* c1ModeName(C1Mode m)
{
    switch (m) {
    case C1Mode::None:      return "none";
    case C1Mode::Fragments: return "fragments";
    case C1Mode::DC:        return "dc";
    case C1Mode::Sparse:    return "sparse";
    }
    return "none";
}

/* ------------------------------------------------------------------------- *
 *  FragmentScfDriver
 *
 *  Orchestrates the dense curcuma::xtb::XTB over sub-systems. Owns the full
 *  molecule, the chosen C1 mode and its accuracy knobs, and after calculate()
 *  exposes the assembled energy, gradient, charges and decomposition.
 * ------------------------------------------------------------------------- */
class FragmentScfDriver {
public:
    FragmentScfDriver(MethodType method, const json& config);

    // Lifecycle — mirrors the ComputationalMethod calls the wrapper forwards.
    bool   setMolecule(const Mol& mol);
    bool   updateGeometry(const Matrix& geometry);
    double calculate(bool gradient);

    // Assembled results (valid after calculate()).
    double         getEnergy() const               { return m_energy; }
    Matrix         getGradient() const             { return m_gradient; }
    Vector         getCharges() const              { return m_charges; }
    nlohmann::json getEnergyDecomposition() const  { return m_decomp; }
    bool           converged() const               { return m_converged; }
    int            numFragments() const            { return m_num_fragments; }

    void setIntraThreads(int n) { m_intra_threads = (n < 1) ? 1 : n; }

    bool hasError() const { return m_has_error; }
    const std::string& errorMessage() const { return m_error_message; }

private:
    // Per-mode drivers.
    double runDense(bool gradient);          // single dense XTB on the whole system
    double runFragments(bool gradient);      // C1c
    double runDivideConquer(bool gradient);  // C1a  (later step)
    double runSparse(bool gradient);         // C1b  (later step)

    // Helpers.
    void configureXtb(XTB& xtb) const;       // apply D4 source + SCF config to a fragment XTB
    Mol  buildSubMol(const std::vector<int>& atoms, int charge) const;  // slice the global Mol

    MethodType m_method;
    json       m_config;
    C1Mode     m_mode = C1Mode::None;

    Mol      m_mol;          // full system (atoms/geometry/charge)
    Molecule m_molecule;     // same system as a Molecule, for GetFragments()

    // C1 knobs (from the xtb config scope; Bohr units for radii).
    double m_buffer_bohr     = 10.0;  // C1a DC buffer radius
    double m_cell_bohr       = 12.0;  // C1a DC cell size
    double m_sparse_threshold = 1.0e-6; // C1b drop tolerance
    int    m_intra_threads   = 1;     // forwarded -threads budget

    // Assembled outputs.
    double         m_energy = 0.0;
    Matrix         m_gradient;
    Vector         m_charges;
    nlohmann::json m_decomp;
    bool           m_converged = false;
    int            m_num_fragments = 1;

    bool        m_has_error = false;
    std::string m_error_message;
};

} // namespace curcuma::xtb
