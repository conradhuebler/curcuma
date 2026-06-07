#!/usr/bin/env python3
# Claude Generated (June 2026): extract tblite ALPB/GBSA solvation parameters.
#
# Parses the Fortran include files under
#   external/tblite/src/tblite/solvation/data/{alpb,cds,shift}/param_*.fh
# and emits a single C++ header
#   src/core/solvation/tblite_solvation_params.h
# containing the per-(method, model, solvent) parameter set used by the native
# GFN1/GFN2 ALPB solvation path (see docs/SQM_SOLVATION_WP.md, WP0.2).
#
# The three directories define identically-named parameters (e.g. gfn2_alpb_water
# appears as an alpb_parameter, a cds_parameter and a shift_parameter), so the
# three records are joined on the Fortran variable name.
#
#   alpb_parameter  = { epsv, c1, soset, sx(94) }
#   cds_parameter   = { rprobe, gamscale(94), sqrtGhbond(94) }
#   shift_parameter = { smass, rhos, gshift }
#
# Name -> (method, model, solvent):  gfn{1,2}[_alpb]_<solvent>
#   "alpb" present  -> model = alpb   (ALPB, analytical linearized PB)
#   "alpb" absent   -> model = gbsa   (plain generalized Born + SASA)

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "external", "tblite", "src", "tblite", "solvation", "data")
OUT = os.path.join(ROOT, "src", "core", "solvation", "tblite_solvation_params.h")

NELEM = 94

# A Fortran real literal: optional sign, digits/decimal, optional d/D/e/E exponent,
# trailing _wp kind suffix.
NUM_RE = re.compile(r"[-+]?\d*\.\d+(?:[dDeE][-+]?\d+)?_wp|[-+]?\d+\.(?:[dDeE][-+]?\d+)?_wp")
# A parameter definition header: "type(X_parameter),parameter :: NAME = X_parameter ("
DEF_RE = re.compile(r"parameter\s*::\s*([A-Za-z0-9_]+)\s*=\s*[A-Za-z0-9_]+_parameter\s*\(",
                    re.IGNORECASE)


def parse_numbers(body):
    """Return all Fortran reals in `body` in order, as Python floats."""
    out = []
    for m in NUM_RE.finditer(body):
        tok = m.group(0).replace("_wp", "").replace("d", "e").replace("D", "e")
        out.append(float(tok))
    return out


def split_definitions(text):
    """Yield (name, body) for each `... :: NAME = X_parameter ( body )` block.

    The body is everything up to the matching closing paren of the constructor.
    """
    for m in DEF_RE.finditer(text):
        name = m.group(1)
        # Balance parentheses starting at the '(' just matched.
        i = m.end() - 1  # index of the opening '('
        depth = 0
        j = i
        while j < len(text):
            c = text[j]
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        body = text[i + 1:j]
        yield name, body


def load_dir(subdir, expected_counts):
    """Parse every param_*.fh in data/<subdir>; return {name: [floats]}.

    expected_counts: set of acceptable number-counts for a sanity check.
    """
    result = {}
    d = os.path.join(DATA, subdir)
    for fn in sorted(os.listdir(d)):
        if not fn.endswith(".fh"):
            continue
        with open(os.path.join(d, fn)) as f:
            text = f.read()
        for name, body in split_definitions(text):
            nums = parse_numbers(body)
            if expected_counts and len(nums) not in expected_counts:
                sys.exit(f"ERROR {subdir}/{fn}: {name} has {len(nums)} numbers "
                         f"(expected one of {sorted(expected_counts)})")
            result[name] = nums
    return result


def parse_name(name):
    """gfn{1,2}[_alpb]_<solvent> -> (method, model, solvent)."""
    m = re.match(r"(gfn[12])_(alpb_)?(.+)$", name)
    if not m:
        sys.exit(f"ERROR: cannot parse parameter name '{name}'")
    method = m.group(1)
    model = "alpb" if m.group(2) else "gbsa"
    solvent = m.group(3)
    return method, model, solvent


def fmt_array(vals):
    """Format a 94-element double array as C++ initializer rows of 5."""
    assert len(vals) == NELEM
    rows = []
    for r in range(0, NELEM, 5):
        chunk = vals[r:r + 5]
        rows.append("        " + ", ".join(f"{v: .8e}" for v in chunk) + ",")
    s = "\n".join(rows)
    return "{\n" + s + "\n    }"


def main():
    alpb = load_dir("alpb", {3 + NELEM})            # epsv,c1,soset + sx[94]
    cds = load_dir("cds", {1 + 2 * NELEM})          # rprobe + gamscale[94] + sqrtGhbond[94]
    shift = load_dir("shift", {3})                  # smass,rhos,gshift

    names = sorted(set(alpb) & set(cds) & set(shift))
    missing = (set(alpb) ^ set(cds)) | (set(alpb) ^ set(shift))
    if missing:
        sys.exit(f"ERROR: parameter names not present in all three dirs: {sorted(missing)}")

    entries = []
    for name in names:
        method, model, solvent = parse_name(name)
        a = alpb[name]
        c = cds[name]
        s = shift[name]
        entries.append(dict(
            key=name, method=method, model=model, solvent=solvent,
            epsv=a[0], c1=a[1], soset=a[2], sx=a[3:3 + NELEM],
            rprobe=c[0], gamscale=c[1:1 + NELEM], sqrtGhbond=c[1 + NELEM:1 + 2 * NELEM],
            smass=s[0], rhos=s[1], gshift=s[2],
        ))

    # Spot-check (docs/SQM_SOLVATION_WP.md WP0.2): gfn2_alpb_water.
    chk = next(e for e in entries if e["key"] == "gfn2_alpb_water")
    assert abs(chk["epsv"] - 80.2) < 1e-9, chk["epsv"]
    assert abs(chk["c1"] - 1.47438678) < 1e-8, chk["c1"]
    assert abs(chk["sx"][0] - 0.18678116) < 1e-8, chk["sx"][0]

    with open(OUT, "w") as f:
        f.write(HEADER_TOP)
        for e in entries:
            f.write(
                f"inline const TbliteSolvationParam tbsolv_{e['key']} = {{\n"
                f"    \"{e['key']}\", \"{e['method']}\", \"{e['model']}\", \"{e['solvent']}\",\n"
                f"    {e['epsv']: .8e}, {e['c1']: .8e}, {e['soset']: .8e},  // epsv, c1, soset\n"
                f"    {e['rprobe']: .8e},  // rprobe (Angstrom)\n"
                f"    {e['smass']: .8e}, {e['rhos']: .8e}, {e['gshift']: .8e},  // smass, rhos, gshift(kcal/mol)\n"
                f"    /* sx */ {fmt_array(e['sx'])},\n"
                f"    /* gamscale */ {fmt_array(e['gamscale'])},\n"
                f"    /* sqrtGhbond */ {fmt_array(e['sqrtGhbond'])}\n"
                f"}};\n\n")

        # Registry array of pointers.
        f.write("inline const TbliteSolvationParam* const kTbliteSolvationParams[] = {\n")
        for e in entries:
            f.write(f"    &tbsolv_{e['key']},\n")
        f.write("};\n")
        f.write(f"inline constexpr int kTbliteSolvationParamCount = {len(entries)};\n\n")
        f.write(LOOKUP_FN)

    print(f"Wrote {OUT}: {len(entries)} entries "
          f"({sum(1 for e in entries if e['model']=='alpb')} alpb, "
          f"{sum(1 for e in entries if e['model']=='gbsa')} gbsa).")
    print("Spot-check gfn2_alpb_water OK: "
          f"epsv={chk['epsv']}, c1={chk['c1']}, sx[0]={chk['sx'][0]}")


HEADER_TOP = """/*
 * <tblite ALPB/GBSA solvation parameters>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * AUTO-GENERATED by scripts/extract_tblite_solvation_params.py — DO NOT EDIT.
 * Source: external/tblite/src/tblite/solvation/data/{alpb,cds,shift}/param_*.fh
 *
 * Per-(method, model, solvent) ALPB/GBSA parameters for the native GFN1/GFN2
 * solvation path. Units match the tblite Fortran reference (see the loaders
 * data/{alpb,cds,shift}.f90): epsv dimensionless; c1/soset/sx feed
 *   born_offset = soset*0.1*aatoau, born_scale = c1, descreening = sx;
 * rprobe in Angstrom (probe = rprobe*aatoau); tension = gamscale*1e-5;
 * hbond_elem = -kcaltoau*sqrtGhbond^2; gshift in kcal/mol (*kcaltoau).
 */

#pragma once

#include <string>

namespace curcuma {
namespace solvation {

/// Combined ALPB(Born) + CDS(non-polar) + shift parameters for one
/// (method, model, solvent). Claude Generated (June 2026).
struct TbliteSolvationParam {
    const char* key;       ///< e.g. "gfn2_alpb_water"
    const char* method;    ///< "gfn1" | "gfn2"
    const char* model;     ///< "alpb" | "gbsa"
    const char* solvent;   ///< canonical solvent name
    double epsv;           ///< dielectric constant
    double c1;             ///< Born radius scaling (born_scale)
    double soset;          ///< Born offset (0.1 Angstrom units)
    double rprobe;         ///< CDS probe radius (Angstrom)
    double smass;          ///< solvent molar mass (g/mol)
    double rhos;           ///< solvent density (g/cm^3)
    double gshift;         ///< empirical free-energy shift (kcal/mol)
    double sx[94];         ///< per-element dielectric descreening
    double gamscale[94];   ///< per-element CDS surface tension scale
    double sqrtGhbond[94]; ///< per-element sqrt(H-bond strength)
};

"""

LOOKUP_FN = """/**
 * @brief Look up the tblite solvation parameter set.
 * @param method  "gfn1" or "gfn2".
 * @param alpb    true for ALPB, false for GBSA.
 * @param solvent Solvent name (aliases resolved, case-insensitive on caller side).
 * @return pointer to the parameter set, or nullptr if unknown.
 *
 * Solvent aliases mirror tblite get_alpb_param (data/alpb.f90).
 * Claude Generated (June 2026).
 */
inline const TbliteSolvationParam* getTbliteSolvationParam(
    const std::string& method, bool alpb, const std::string& solvent)
{
    // Resolve solvent aliases to the canonical file name used above.
    std::string s = solvent;
    auto alias = [&](const char* canon, std::initializer_list<const char*> names) {
        for (const char* n : names)
            if (s == n) { s = canon; return; }
    };
    alias("ch2cl2", {"dichlormethane", "methylenechloride", "dcm"});
    alias("chcl3", {"chloroform"});
    alias("cs2", {"carbondisulfide"});
    alias("dmf", {"dimethylformamide"});
    alias("dmso", {"dimethylsulfoxide"});
    alias("ether", {"diethylether"});
    alias("thf", {"tetrahydrofuran"});
    alias("hexane", {"nhexan", "n-hexan", "nhexane", "n-hexane"});
    alias("water", {"h2o"});

    const char* model = alpb ? "alpb" : "gbsa";
    for (int i = 0; i < kTbliteSolvationParamCount; ++i) {
        const TbliteSolvationParam* p = kTbliteSolvationParams[i];
        if (method == p->method && std::string(model) == p->model && s == p->solvent)
            return p;
    }
    return nullptr;
}

} // namespace solvation
} // namespace curcuma
"""


if __name__ == "__main__":
    main()
