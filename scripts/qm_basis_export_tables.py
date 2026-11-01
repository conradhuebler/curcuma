#!/usr/bin/env python3
"""Export the built-in basis sets of the native QM engine (H-Ne) into a C++ header.

Writes src/core/energy_calculators/qm_methods/qm_basis_tables.h with the text of five basis sets in the formats
the engine's parsers read (BasisSetParser::parseBasisSetStream for '$DATA', qmri::parseTurbomoleBasisStream for
Turbomole):

    def2-SVP, def2-SV(P)                      orbital, '$DATA' blocks   <- BSE tests/sources/*.tm.bz2
    def2-universal-JFIT, def2-universal-JKFIT auxiliary, Turbomole     <- BSE tests/sources/*.tm.bz2
    MINIX (HF-3c)                             orbital, '$DATA' blocks   <- BSE data/huzinaga/MINI-scaled.0.json
                                                                           + the Li/Be p shells of MINIX (below)

Source: the Basis Set Exchange repository https://github.com/MolSSI-BSE/basis_set_exchange (BSD-3-Clause,
Copyright (c) 2020 The Molecular Sciences Software Institute, Virginia Tech), pinned to commit BSE_COMMIT. The
license text is reproduced in the generated header and in docs/THIRD_PARTY_NOTICES.md.

What is copied how:
  * def2 orbital and fitting sets: every exponent and coefficient token of H-Ne is copied unchanged from the BSE
    file (only a Fortran 'D' exponent marker would be turned into 'E'; there is none for H-Ne). Shell order as in
    the BSE file.
  * MINIX: Sure & Grimme, J. Comput. Chem. 34, 1672-1685 (2013), Table 1 (composition of the MINIX basis):
        H-He, B-Ne : MINIS
        Li-Be      : MINIS + 1(p)
        Na-Mg      : MINIS + 1(p)
        Al-Ar      : MINIS + 1(d)
        K-Zn       : SV
        Ga-Kr      : SVP
        Rb-Xe      : def2-SV(P) with ECP
    Only H-Ne is exported (the range of the native engine; Rb-Xe would need ECPs, which the engine lacks).
    MINIS = BSE 'Scaled MINI' (MINI-scaled.0.json, original BSE contributor David Feller). Exponents: the BSE
    values; coefficients: the BSE values with each contraction renormalised to unit self-overlap of the
    contracted function (normalised primitives), printed with 10 decimals. This is the normalisation ORCA 6.1
    applies on export (orca_exportbasis -b minix); the result equals the ORCA export and the shipped
    qm_methods/MINIX.dat token for token (checked below). The added p shell of Li (exponent 0.1) and Be
    (exponent 0.25), coefficient 1, is not in BSE: it is taken from ORCA 6.1 (orca_exportbasis -b minix, the
    origin of MINIX.dat).

Checks done on every run (also with --check):
  * the parsed numbers of the embedded text equal, as IEEE doubles, the numbers of the shipped files in
    src/core/energy_calculators/qm_methods/ for H-Ne (same shells, same order) -- so the C++ parse result of the
    embedded copy and of the file is the same;
  * --compare-orca DIR: compares with ORCA exports (orca_exportbasis -b minix|def2-svp|def2-sv(p)|def2/j|def2/jk,
    files minix.bas, def2-svp.bas, def2-sv(p).bas, def2_j.bas, def2_jk.bas in DIR) and prints the maximum relative
    differences of exponents and of the contraction coefficients after normalising both sides (ORCA renormalises
    contractions; for the fitting sets it also sets single-primitive coefficients to 1).

Usage:
    python3 scripts/qm_basis_export_tables.py [--cache DIR] [--compare-orca DIR]   # write the header
    python3 scripts/qm_basis_export_tables.py --check [--cache DIR]                # fail if the header differs
Without --cache the files are downloaded (curl) from raw.githubusercontent.com at BSE_COMMIT; with --cache they
are read from DIR (and downloaded into it when missing).

Claude Generated (Nov 2026). Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>. GPL-3.0.
"""
import argparse
import bz2
import json
import math
import os
import re
import subprocess
import sys
import urllib.parse

BSE_COMMIT = "4adaf1372c7101620ca1a9f3130be9ae97fb8f30"  # master on 2026-08-07, downloaded 2026-11-01
BSE_RAW = "https://raw.githubusercontent.com/MolSSI-BSE/basis_set_exchange/" + BSE_COMMIT + "/"

HERE = os.path.dirname(os.path.abspath(__file__))
QMDIR = os.path.join(HERE, "..", "src", "core", "energy_calculators", "qm_methods")
OUT = os.path.join(QMDIR, "qm_basis_tables.h")

ELEMENTS = ["h", "he", "li", "be", "b", "c", "n", "o", "f", "ne"]
NAMES = ["HYDROGEN", "HELIUM", "LITHIUM", "BERYLLIUM", "BORON", "CARBON", "NITROGEN", "OXYGEN", "FLUORINE", "NEON"]
LETTERS = "spdfghi"

# (embedded name, kind, BSE path, repository file)
SETS = [
    ("def2-SVP", "orbital", "basis_set_exchange/tests/sources/def2-svp.1.tm.bz2", "def2-SVP.dat"),
    ("def2-SV(P)", "orbital", "basis_set_exchange/tests/sources/def2-sv(p).1.tm.bz2", "def2-SV(P).dat"),
    ("def2-universal-JFIT", "jbas", "basis_set_exchange/tests/sources/def2-universal-jfit.1.tm.bz2",
     "def2-universal-JFIT.turbomole"),
    ("def2-universal-JKFIT", "jkbas", "basis_set_exchange/tests/sources/def2-universal-jkfit.1.tm.bz2",
     "def2-universal-JKFIT.turbomole"),
    ("MINIX", "orbital", "basis_set_exchange/data/huzinaga/MINI-scaled.0.json", "MINIX.dat"),
]
LICENSE_PATH = "LICENSE"

# MINIX = MINIS + 1(p) for Li and Be (Sure & Grimme 2013, Table 1). Numbers from ORCA 6.1 orca_exportbasis -b minix.
MINIX_EXTRA_P = {"li": ("0.1000000000", "1.0000000000"), "be": ("0.2500000000", "1.0000000000")}


# ------------------------------------------------------------------------------------------------- input
def fetch(path, cache):
    """Bytes of one file of the BSE repository at BSE_COMMIT (from the cache directory if present)."""
    local = os.path.join(cache, os.path.basename(path)) if cache else None
    if local and os.path.exists(local):
        with open(local, "rb") as f:
            return f.read()
    url = BSE_RAW + urllib.parse.quote(path)
    data = subprocess.run(["curl", "-sSfL", url], check=True, capture_output=True).stdout
    if local:
        os.makedirs(cache, exist_ok=True)
        with open(local, "wb") as f:
            f.write(data)
    return data


def parse_tm(text):
    """Turbomole basis text -> {element: [(L, [(exp_token, coef_token), ...]), ...]} (tokens as written)."""
    out, elem, expect = {}, None, False
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        t = lines[i].strip()
        i += 1
        if not t or t.startswith("#"):
            continue
        if t.startswith("$"):
            if t.startswith("$end") or t.startswith("$ecp"):  # the ECP section of the BSE files is not read
                break
            elem, expect = None, False
            continue
        if t.startswith("*"):
            elem, expect = None, True
            continue
        if t[0].isalpha():
            assert expect, "element header without '*': " + t
            elem = t.split()[0].lower()
            out[elem] = []
            while lines[i].strip() == "" or lines[i].strip().startswith("#"):
                i += 1
            assert lines[i].strip().startswith("*"), lines[i]
            i += 1
            expect = False
            continue
        n, l = t.split()[:2]
        prims = []
        for _ in range(int(n)):
            tok = lines[i].split()
            i += 1
            prims.append((tok[0], tok[1]))
        out[elem].append((LETTERS.index(l.lower()), prims))
    return out


def parse_dat(text):
    """curcuma '$DATA' text -> {element: [(L, [(exp_token, coef_token), ...]), ...]} (H-Ne, full names)."""
    out, elem = {}, None
    lines = text.split("$DATA", 1)[1].splitlines()
    i = 0
    while i < len(lines):
        t = lines[i].split()
        i += 1
        if not t or t[0].startswith("!"):
            continue
        if t[0].startswith("$END"):
            break
        if t[0].upper() in NAMES:
            elem = ELEMENTS[NAMES.index(t[0].upper())]
            out[elem] = []
            continue
        if len(t) == 1 or t[0].isalpha() and not t[1][0].isdigit():
            elem = None  # element outside H-Ne (skipped like the C++ parser's map)
            continue
        if elem is None:
            i += int(t[1])
            continue
        prims = []
        for _ in range(int(t[1])):
            tok = lines[i].split()
            i += 1
            prims.append((tok[1], tok[2]))
        out[elem].append(("SPDFG".index(t[0].upper()), prims))
    return out


def num(tok):
    return float(tok.replace("D", "E").replace("d", "e"))


# ------------------------------------------------------------------------------------------------- MINIX
def self_overlap(l, exps, coefs):
    s = 0.0
    for a, ca in zip(exps, coefs):
        for b, cb in zip(exps, coefs):
            s += ca * cb * (2.0 * math.sqrt(a * b) / (a + b)) ** (l + 1.5)
    return s


def normalised(l, exps, coefs):
    n = math.sqrt(self_overlap(l, exps, coefs))
    return [c / n for c in coefs]


def minix_from_bse(js):
    """MINIX H-Ne: BSE MINI-scaled (renormalised, 10 decimals) + the Li/Be p shell."""
    out = {}
    for z, el in enumerate(ELEMENTS, start=1):
        shells = []
        for sh in js["elements"][str(z)]["electron_shells"]:
            (l,) = sh["angular_momentum"]
            assert len(sh["coefficients"]) == 1, "general contraction in MINI-scaled"
            exps = [float(x) for x in sh["exponents"]]
            coefs = normalised(l, exps, [float(x) for x in sh["coefficients"][0]])
            shells.append((l, [("%.10f" % e, "%.10f" % c) for e, c in zip(exps, coefs)]))
            for tok, e in zip(shells[-1][1], exps):
                assert float(tok[0]) == e, "exponent not exact at 10 decimals"
        if el in MINIX_EXTRA_P:
            shells.append((1, [MINIX_EXTRA_P[el]]))
        out[el] = shells
    return out


# ------------------------------------------------------------------------------------------------- text
def dat_block(el, shells):
    lines = ["", NAMES[ELEMENTS.index(el)]]
    for l, prims in shells:
        lines.append("    %s   %d" % ("SPDFG"[l], len(prims)))
        for k, (e, c) in enumerate(prims, start=1):
            lines.append("    %d  %20s  %20s" % (k, e.replace("D", "E"), c.replace("D", "E")))
    return "\n".join(lines) + "\n"


def tm_block(el, name, shells):
    lines = ["*", "%s %s" % (el, name), "*"]
    for l, prims in shells:
        lines.append("    %d   %s" % (len(prims), LETTERS[l]))
        for e, c in prims:
            lines.append("    %20s  %20s" % (e, c))
    return "\n".join(lines) + "\n"


def compare(name, emb, ref, what):
    """Assert identical doubles and shell structure for H-Ne; returns the number of compared values."""
    count = 0
    for el in ELEMENTS:
        a, b = emb[el], ref[el]
        assert [(l, len(p)) for l, p in a] == [(l, len(p)) for l, p in b], (name, el, what, "shell structure differs")
        for (l, pa), (_, pb) in zip(a, b):
            for (ea, ca), (eb, cb) in zip(pa, pb):
                assert num(ea) == num(eb) and num(ca) == num(cb), (name, el, what, ea, eb, ca, cb)
                count += 2
    return count


def compare_orca(name, emb, orca_file):
    """Max relative differences against an ORCA NewGTO export (exponents; normalised coefficient vectors)."""
    data, el, cur = {}, None, None
    for line in open(orca_file):
        t = line.split()
        if not t:
            continue
        if t[0].startswith("New") and len(t) >= 2:
            el = t[1].lower()
            data[el] = []
            continue
        if el is None:
            continue
        if t[0].startswith("end"):
            el = None
            continue
        if t[0].upper() in "SPDFGHI" and len(t) == 2:
            cur = (("SPDFGHI".index(t[0].upper())), [])
            data[el].append(cur)
            continue
        if len(t) == 3 and t[0].isdigit():
            cur[1].append((t[1], t[2]))
    de = dc_raw = dc_norm = ae = ac = 0.0
    for e in ELEMENTS:
        a, b = emb[e], data[e]
        assert [(l, len(p)) for l, p in a] == [(l, len(p)) for l, p in b], (name, e, "ORCA shell structure differs")
        for (l, pa), (_, pb) in zip(a, b):
            ea = [num(x) for x, _ in pa]
            eb = [num(x) for x, _ in pb]
            ca = [num(y) for _, y in pa]
            cb = [num(y) for _, y in pb]
            de = max(de, max(abs(x - y) / abs(y) for x, y in zip(ea, eb)))
            ae = max(ae, max(abs(x - y) for x, y in zip(ea, eb)))
            dc_raw = max(dc_raw, max(abs(x - y) / abs(y) for x, y in zip(ca, cb)))
            na, nb = normalised(l, ea, ca), normalised(l, eb, cb)
            dc_norm = max(dc_norm, max(abs(x - y) / abs(y) for x, y in zip(na, nb)))
            # ORCA prints 10 decimals: compare in ORCA's scale (its coefficients rescaled to the embedded norm)
            scale = math.sqrt(self_overlap(l, eb, cb))
            ac = max(ac, max(abs(x * scale - y) for x, y in zip(na, cb)))
    print("ORCA %-22s max rel diff: exponents %.1e, coefficients raw %.1e, normalised %.1e;"
          " max abs diff (ORCA prints 10 decimals): exponents %.1e, coefficients %.1e"
          % (name, de, dc_raw, dc_norm, ae, ac))
    return de, dc_raw, dc_norm


def cxx_literal(text):
    """Adjacent raw string literals, one per element block (each well below the 16 KB MSVC limit)."""
    out = []
    for chunk in text.split("\x00"):
        if chunk:
            assert ")BASIS\"" not in chunk
            out.append('    R"BASIS(' + chunk + ')BASIS"')
    return "\n".join(out)


def citations():
    """BibTeX entries of the embedded sets, copied from src/core/citation_database.cpp."""
    src = open(os.path.join(HERE, "..", "src", "core", "citation_database.cpp"), encoding="utf-8").read()
    out = []
    lines = src.splitlines()
    for key in ("def2_basis", "def2_sv_1992", "def2_jfit", "def2_jkfit", "hf3c", "bse_2019"):
        start = [i for i, l in enumerate(lines) if l.strip().startswith('{ "%s", {' % key)]
        assert len(start) == 1, key
        i = start[0]
        while not lines[i].strip().startswith('"@'):
            i += 1
        entry = []
        while True:
            s = lines[i].strip()
            assert s.startswith('"') and s.endswith('"'), (key, s)
            entry.append(re.sub(r"\\(.)", lambda m: {"n": "", '"': '"', "\\": "\\"}[m.group(1)], s[1:-1]))
            if s == '"}"':
                break
            i += 1
        out.append("\n".join(entry))
    return out


def build(cache):
    license_text = fetch(LICENSE_PATH, cache).decode()
    sets = []
    nchecked = 0
    for name, kind, path, repo in SETS:
        raw = fetch(path, cache)
        if path.endswith(".bz2"):
            data = parse_tm(bz2.decompress(raw).decode())
        else:
            data = minix_from_bse(json.loads(raw.decode()))
        data = {el: data[el] for el in ELEMENTS}
        repo_text = open(os.path.join(QMDIR, repo), encoding="utf-8").read()
        ref = parse_dat(repo_text) if repo.endswith(".dat") else parse_tm(repo_text)
        nchecked += compare(name, data, ref, repo)
        if kind == "orbital":
            text = "!----------------------------------------------------------------------\n" \
                   "! %s (H-Ne), built into curcuma; see qm_basis_tables.h for the source\n" \
                   "!----------------------------------------------------------------------\n\n$DATA\n\x00" % name
            text += "\x00".join(dat_block(el, data[el]) for el in ELEMENTS) + "\x00\n$END\n"
        else:
            text = "$%s\n\x00" % kind + "\x00".join(tm_block(el, name, data[el]) for el in ELEMENTS) + "\x00*\n$end\n"
        sets.append((name, "data" if kind == "orbital" else "turbomole", text, data))
    return license_text, sets, nchecked


def render(license_text, sets):
    w = []
    w.append("/*")
    w.append(" * <Native QM: built-in basis sets def2-SVP, def2-SV(P), def2-universal-JFIT/-JKFIT, MINIX (H-Ne)>")
    w.append(" * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>")
    w.append(" *")
    w.append(" * GENERATED by scripts/qm_basis_export_tables.py - do not edit.")
    w.append(" * Claude Generated (Nov 2026). AI-generated, machine-tested.")
    w.append(" *")
    w.append(" * Source: Basis Set Exchange repository, https://github.com/MolSSI-BSE/basis_set_exchange")
    w.append(" *   commit " + BSE_COMMIT + " (master of 2026-08-07, downloaded 2026-11-01)")
    for name, kind, path, repo in SETS:
        w.append(" *   %-22s <- %s" % (name, path))
    w.append(" * def2 sets: unchanged from the BSE download for H-Ne (every exponent/coefficient token copied as is; the")
    w.append(" *   BSE file headers state 'Data from Turbomole 7.3'; original BSE contributor Florian Weigend, PI F. Weigend")
    w.append(" *   and R. Ahlrichs).")
    w.append(" * MINIX section: Sure & Grimme, J. Comput. Chem. 34, 1672-1685 (2013), Table 1 (MINIX composition):")
    w.append(" *     H-He, B-Ne: MINIS; Li-Be: MINIS+1(p); Na-Mg: MINIS+1(p); Al-Ar: MINIS+1(d); K-Zn: SV;")
    w.append(" *     Ga-Kr: SVP; Rb-Xe: def2-SV(P) with ECP.  Only H-Ne is built in.")
    w.append(" *   MINIS = BSE 'Scaled MINI' (data/huzinaga/MINI-scaled.0.json, original BSE contributor David Feller):")
    w.append(" *   exponents as in BSE; coefficients = BSE coefficients renormalised to a unit-norm contracted function,")
    w.append(" *   10 decimals (= ORCA 6.1 orca_exportbasis -b minix = qm_methods/MINIX.dat). Li/Be extra p shell")
    w.append(" *   (exponent 0.1 / 0.25, coefficient 1): from ORCA 6.1 orca_exportbasis -b minix (not in BSE).")
    w.append(" * The parsed numbers equal those of the files in qm_methods/ (checked by the generator); a file found in the")
    w.append(" * basis search path takes precedence over these copies (qm_basis_embedded.h).")
    w.append(" *")
    w.append(" * References (BibTeX from src/core/citation_database.cpp):")
    for entry in citations():
        for line in entry.splitlines():
            w.append(" *   " + line)
    w.append(" *")
    w.append(" * Basis Set Exchange license (BSD-3-Clause), reproduced verbatim:")
    w.append(" *")
    for line in license_text.rstrip("\n").splitlines():
        w.append((" * " + line).rstrip())
    w.append(" *")
    w.append(" * This program is free software under GPL-3.0 (the derived header); the BSE data keep their license above.")
    w.append(" */")
    w.append("")
    w.append("#pragma once")
    w.append("")
    w.append("namespace qmbasis {")
    w.append("namespace tables {")
    w.append("")
    w.append("/// One built-in basis set: name as used by -qm.basis / -qm.ri_aux, format of `text`")
    w.append("/// (\"data\": '$DATA' blocks for BasisSetParser, \"turbomole\": Turbomole for qmri), and the file text.")
    w.append("struct EmbeddedBasisSet {")
    w.append("    const char* name;")
    w.append("    const char* format;")
    w.append("    const char* text;")
    w.append("};")
    w.append("")
    for i, (name, fmt_, text, _) in enumerate(sets):
        w.append("// %s (H-Ne)" % name)
        w.append("inline constexpr const char kText%d[] =" % i)
        w.append(cxx_literal(text) + ";")
        w.append("")
    w.append("inline constexpr EmbeddedBasisSet kEmbeddedBasisSets[] = {")
    for i, (name, fmt_, _, _) in enumerate(sets):
        w.append('    { "%s", "%s", kText%d },' % (name, fmt_, i))
    w.append("};")
    w.append("")
    w.append("}  // namespace tables")
    w.append("}  // namespace qmbasis")
    return "\n".join(w) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cache", help="directory with (or for) the BSE files")
    ap.add_argument("--check", action="store_true", help="fail if the checked-in header differs")
    ap.add_argument("--compare-orca", metavar="DIR", help="directory with orca_exportbasis outputs")
    a = ap.parse_args()
    license_text, sets, nchecked = build(a.cache)
    print("embedded numbers equal to the files in qm_methods/: %d values (H-Ne, 5 sets)" % nchecked)
    if a.compare_orca:
        files = {"def2-SVP": "def2-svp.bas", "def2-SV(P)": "def2-sv(p).bas", "def2-universal-JFIT": "def2_j.bas",
                 "def2-universal-JKFIT": "def2_jk.bas", "MINIX": "minix.bas"}
        for name, _, _, data in sets:
            compare_orca(name, data, os.path.join(a.compare_orca, files[name]))
    text = render(license_text, sets)
    if a.check:
        old = open(OUT, encoding="utf-8").read() if os.path.exists(OUT) else ""
        if old != text:
            print("qm_basis_tables.h differs from the generator output; rerun scripts/qm_basis_export_tables.py")
            return 1
        print("qm_basis_tables.h is up to date")
        return 0
    with open(OUT, "w", encoding="utf-8") as f:
        f.write(text)
    print("wrote", os.path.relpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
