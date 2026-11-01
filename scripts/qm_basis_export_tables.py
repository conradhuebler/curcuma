#!/usr/bin/env python3
"""Export the built-in basis sets of the native QM engine into generated C++ headers.

Stage 1 (Nov 2026) built in def2-SVP, def2-SV(P), def2-universal-JFIT/-JKFIT and MINIX for H-Ne. Stage 2 (Nov 2026)
built in every def2 orbital and auxiliary set of the Basis Set Exchange and the def2 ECP, for every element the BSE
has; MINIX stays H-Ne. Output (src/core/energy_calculators/qm_methods/):

    qm_basis_data_def2_sv.h     def2-SVP, def2-SV(P), def2-SVPD                         orbital, '$DATA' blocks
    qm_basis_data_def2_tzv.h    def2-TZVP, def2-TZVPP, def2-TZVPD, def2-TZVPPD           orbital, '$DATA' blocks
    qm_basis_data_def2_qzv.h    def2-QZVP, def2-QZVPP, def2-QZVPD, def2-QZVPPD           orbital, '$DATA' blocks
    qm_basis_data_def2_mtzv.h   def2-mTZVP, def2-mTZVPP                                 orbital, '$DATA' blocks
    qm_basis_data_def2_jk.h     def2-universal-JFIT, -JKFIT, def2-SV(P)-JKFIT, def2-mTZVPP-RIJ   Turbomole
    qm_basis_data_def2_rifit.h  def2-{SVP,SV(P),SVPD,TZVP,TZVPP,TZVPD,TZVPPD,QZVP,QZVPP,QZVPPD}-RIFIT Turbomole
    qm_basis_data_def2_ecp.h    def2-ECP                                                Turbomole '$ecp'
    qm_basis_data_minix.h       MINIX (H-Ne)                                            orbital, '$DATA' blocks
    qm_basis_tables.h           registry kEmbeddedBasisSets[] (included only by qm_basis_embedded.cpp)

Every set is stored per element (blocks[Z-1], nullptr for an element the set does not have), so that the C++
lookup (qm_basis_embedded.cpp) can hand the parsers the elements they need and give raw access per element.

Source: the Basis Set Exchange repository https://github.com/MolSSI-BSE/basis_set_exchange (BSD-3-Clause,
Copyright (c) 2020 The Molecular Sciences Software Institute, Virginia Tech), pinned to commit BSE_COMMIT. The
license text is reproduced in every generated header and in docs/THIRD_PARTY_NOTICES.md.

What is copied how:
  * def2 orbital, fitting and ECP sets: the Turbomole files basis_set_exchange/tests/sources/<name>.1.tm.bz2 (the
    files the BSE curated its def2 data from; header 'Data from Turbomole 7.3'). Every exponent/coefficient token is
    copied unchanged except that a Fortran 'D' exponent marker becomes 'E' (the '$DATA' parser reads only 'E').
    Shell order as in the .tm file. The ECP lines are copied verbatim (whitespace normalised).
  * MINIX: Sure & Grimme, J. Comput. Chem. 34, 1672-1685 (2013), Table 1 (composition of the MINIX basis):
        H-He, B-Ne : MINIS;  Li-Be, Na-Mg : MINIS + 1(p);  Al-Ar : MINIS + 1(d);  K-Zn : SV;  Ga-Kr : SVP;
        Rb-Xe      : def2-SV(P) with ECP
    Only H-Ne is exported (Rb-Xe would need def2-SV(P) plus the def2 ECP in the SCF, which the engine lacks).
    MINIS = BSE 'Scaled MINI' (MINI-scaled.0.json, original BSE contributor David Feller); coefficients renormalised
    to a unit-norm contraction, 10 decimals (= ORCA 6.1 orca_exportbasis -b minix = qm_methods/MINIX.dat); the Li/Be
    p shell (exponent 0.1 / 0.25, coefficient 1) is not in BSE and is taken from ORCA 6.1.

Checks done on every run (also with --check):
  (1) every set against the BSE json data composed the way BSE composes it (data/<name>.1.table.json ->
      element.json -> component files, basis_set_exchange/compose.py): same elements, and per element the same
      shells with the same numbers as IEEE doubles (shells matched as a multiset; a different shell ORDER and a
      different spelling of the same number are counted and printed, not an error); ECP: ncore, potentials,
      r powers, exponents and coefficients;
  (2) the ECP sections inside the orbital .tm files equal def2-ecp.1.tm.bz2;
  (3) the five stage-1 sets equal the files in qm_methods/ for H-Ne (IEEE doubles, same order), and the elements
      beyond Ne where the shipped def2-SV(P).dat / def2-SVP.dat differ from BSE are listed;
  (4) --compare-orca DIR: H-Ne against ORCA 6.1 exports (orca_exportbasis -b <name> > <name with / -> _>.bas), max
      relative differences of exponents and of the coefficients after normalising both sides.

Usage:
    python3 scripts/qm_basis_export_tables.py [--cache DIR] [--compare-orca DIR]   # write the headers
    python3 scripts/qm_basis_export_tables.py --check [--cache DIR]                # fail if a header differs
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
BSE_DATA = "basis_set_exchange/data/"
BSE_TM = "basis_set_exchange/tests/sources/"

HERE = os.path.dirname(os.path.abspath(__file__))
QMDIR = os.path.join(HERE, "..", "src", "core", "energy_calculators", "qm_methods")

SYMBOLS = ("h he li be b c n o f ne na mg al si p s cl ar k ca sc ti v cr mn fe co ni cu zn ga ge as se br kr "
           "rb sr y zr nb mo tc ru rh pd ag cd in sn sb te i xe cs ba la ce pr nd pm sm eu gd tb dy ho er tm yb "
           "lu hf ta w re os ir pt au hg tl pb bi po at rn fr ra ac th pa u np pu am cm bk cf es fm md no lr").split()
# element names as in the BSE '$DATA' (GAMESS) export
NAMES = ("HYDROGEN HELIUM LITHIUM BERYLLIUM BORON CARBON NITROGEN OXYGEN FLUORINE NEON SODIUM MAGNESIUM ALUMINIUM "
         "SILICON PHOSPHORUS SULFUR CHLORINE ARGON POTASSIUM CALCIUM SCANDIUM TITANIUM VANADIUM CHROMIUM MANGANESE "
         "IRON COBALT NICKEL COPPER ZINC GALLIUM GERMANIUM ARSENIC SELENIUM BROMINE KRYPTON RUBIDIUM STRONTIUM "
         "YTTRIUM ZIRCONIUM NIOBIUM MOLYBDENUM TECHNETIUM RUTHENIUM RHODIUM PALLADIUM SILVER CADMIUM INDIUM TIN "
         "ANTIMONY TELLURIUM IODINE XENON CAESIUM BARIUM LANTHANUM CERIUM PRASEODYMIUM NEODYMIUM PROMETHIUM "
         "SAMARIUM EUROPIUM GADOLINIUM TERBIUM DYSPROSIUM HOLMIUM ERBIUM THULIUM YTTERBIUM LUTETIUM HAFNIUM "
         "TANTALUM TUNGSTEN RHENIUM OSMIUM IRIDIUM PLATINUM GOLD MERCURY THALLIUM LEAD BISMUTH POLONIUM ASTATINE "
         "RADON FRANCIUM RADIUM ACTINIUM THORIUM PROTACTINIUM URANIUM NEPTUNIUM PLUTONIUM AMERICIUM CURIUM "
         "BERKELIUM CALIFORNIUM EINSTEINIUM FERMIUM MENDELEVIUM NOBELIUM LAWRENCIUM").split()
assert len(SYMBOLS) == len(NAMES) == 103
HNE = SYMBOLS[:10]
LETTERS = "spdfghi"

# (embedded name, kind, BSE .tm file (or json for MINIX), header family)
#   kind: orbital ('$DATA' text), jbas / jkbas / cbas (Turbomole fitting text), ecp (Turbomole '$ecp' text)
SETS = [
    ("def2-SVP", "orbital", "def2-svp.1.tm.bz2", "def2_sv"),
    ("def2-SV(P)", "orbital", "def2-sv(p).1.tm.bz2", "def2_sv"),
    ("def2-SVPD", "orbital", "def2-svpd.1.tm.bz2", "def2_sv"),
    ("def2-TZVP", "orbital", "def2-tzvp.1.tm.bz2", "def2_tzv"),
    ("def2-TZVPP", "orbital", "def2-tzvpp.1.tm.bz2", "def2_tzv"),
    ("def2-TZVPD", "orbital", "def2-tzvpd.1.tm.bz2", "def2_tzv"),
    ("def2-TZVPPD", "orbital", "def2-tzvppd.1.tm.bz2", "def2_tzv"),
    ("def2-QZVP", "orbital", "def2-qzvp.1.tm.bz2", "def2_qzv"),
    ("def2-QZVPP", "orbital", "def2-qzvpp.1.tm.bz2", "def2_qzv"),
    ("def2-QZVPD", "orbital", "def2-qzvpd.1.tm.bz2", "def2_qzv"),
    ("def2-QZVPPD", "orbital", "def2-qzvppd.1.tm.bz2", "def2_qzv"),
    ("def2-mTZVP", "orbital", "def2-mtzvp.1.tm.bz2", "def2_mtzv"),
    ("def2-mTZVPP", "orbital", "def2-mtzvpp.1.tm.bz2", "def2_mtzv"),
    ("def2-universal-JFIT", "jbas", "def2-universal-jfit.1.tm.bz2", "def2_jk"),
    ("def2-universal-JKFIT", "jkbas", "def2-universal-jkfit.1.tm.bz2", "def2_jk"),
    ("def2-SV(P)-JKFIT", "jkbas", "def2-sv(p)-jkfit.1.tm.bz2", "def2_jk"),
    ("def2-mTZVPP-RIJ", "jbas", "def2-mtzvpp-rij.1.tm.bz2", "def2_jk"),
    ("def2-SVP-RIFIT", "cbas", "def2-svp-rifit.1.tm.bz2", "def2_rifit"),
    ("def2-SV(P)-RIFIT", "cbas", "def2-sv(p)-rifit.1.tm.bz2", "def2_rifit"),
    ("def2-SVPD-RIFIT", "cbas", "def2-svpd-rifit.1.tm.bz2", "def2_rifit"),
    ("def2-TZVP-RIFIT", "cbas", "def2-tzvp-rifit.1.tm.bz2", "def2_rifit"),
    ("def2-TZVPP-RIFIT", "cbas", "def2-tzvpp-rifit.1.tm.bz2", "def2_rifit"),
    ("def2-TZVPD-RIFIT", "cbas", "def2-tzvpd-rifit.1.tm.bz2", "def2_rifit"),
    ("def2-TZVPPD-RIFIT", "cbas", "def2-tzvppd-rifit.1.tm.bz2", "def2_rifit"),
    ("def2-QZVP-RIFIT", "cbas", "def2-qzvp-rifit.1.tm.bz2", "def2_rifit"),
    ("def2-QZVPP-RIFIT", "cbas", "def2-qzvpp-rifit.1.tm.bz2", "def2_rifit"),
    ("def2-QZVPPD-RIFIT", "cbas", "def2-qzvppd-rifit.1.tm.bz2", "def2_rifit"),
    ("def2-ECP", "ecp", "def2-ecp.1.tm.bz2", "def2_ecp"),
    ("MINIX", "orbital", "huzinaga/MINI-scaled.0.json", "minix"),
]
FAMILIES = [
    ("def2_sv", "def2-SVP, def2-SV(P), def2-SVPD"),
    ("def2_tzv", "def2-TZVP, def2-TZVPP, def2-TZVPD, def2-TZVPPD"),
    ("def2_qzv", "def2-QZVP, def2-QZVPP, def2-QZVPD, def2-QZVPPD"),
    ("def2_mtzv", "def2-mTZVP, def2-mTZVPP"),
    ("def2_jk", "def2-universal-JFIT, def2-universal-JKFIT, def2-SV(P)-JKFIT, def2-mTZVPP-RIJ"),
    ("def2_rifit", "def2-*-RIFIT (SVP, SV(P), SVPD, TZVP, TZVPP, TZVPD, TZVPPD, QZVP, QZVPP, QZVPPD)"),
    ("def2_ecp", "def2-ECP (effective core potentials; NOT used by the engine yet)"),
    ("minix", "MINIX (H-Ne)"),
]
# stage-1 sets: embedded H-Ne must equal the shipped files (same shells, same order)
REPO_FILES = {"def2-SVP": "def2-SVP.dat", "def2-SV(P)": "def2-SV(P).dat", "MINIX": "MINIX.dat",
              "def2-universal-JFIT": "def2-universal-JFIT.turbomole",
              "def2-universal-JKFIT": "def2-universal-JKFIT.turbomole"}
# ORCA 6.1 names (orca_exportbasis -b) of the sets ORCA has; compared for H-Ne with --compare-orca
ORCA_NAMES = {"def2-SVP": "def2-svp", "def2-SV(P)": "def2-sv(p)", "def2-SVPD": "def2-svpd", "def2-TZVP": "def2-tzvp",
              "def2-TZVPP": "def2-tzvpp", "def2-TZVPD": "def2-tzvpd", "def2-TZVPPD": "def2-tzvppd",
              "def2-QZVP": "def2-qzvp", "def2-QZVPP": "def2-qzvpp", "def2-QZVPD": "def2-qzvpd",
              "def2-QZVPPD": "def2-qzvppd", "def2-mTZVP": "def2-mtzvp", "def2-mTZVPP": "def2-mtzvpp",
              "def2-universal-JFIT": "def2/j", "def2-universal-JKFIT": "def2/jk", "def2-mTZVPP-RIJ": "def2-mtzvpp/j",
              "def2-SVP-RIFIT": "def2-svp/c", "def2-SVPD-RIFIT": "def2-svpd/c", "def2-TZVP-RIFIT": "def2-tzvp/c",
              "def2-TZVPP-RIFIT": "def2-tzvpp/c", "def2-TZVPD-RIFIT": "def2-tzvpd/c",
              "def2-TZVPPD-RIFIT": "def2-tzvppd/c", "def2-QZVPP-RIFIT": "def2-qzvpp/c",
              "def2-QZVPPD-RIFIT": "def2-qzvppd/c", "MINIX": "minix"}
# BSE table file stem where it is not the set name
TABLE_NAMES = {"def2-mTZVP": "def2-mtzvp", "def2-mTZVPP": "def2-mtzvpp", "def2-mTZVPP-RIJ": "def2-mtzvpp-rij"}
LICENSE_PATH = "LICENSE"
MINIX_EXTRA_P = {"li": ("0.1000000000", "1.0000000000"), "be": ("0.2500000000", "1.0000000000")}
# literature per family (keys of src/core/citation_database.cpp; BibTeX copied into the header)
FAMILY_CITES = {
    "def2_sv": ["def2_basis", "def2_sv_1992", "def2_diffuse", "bse_2019"],
    "def2_tzv": ["def2_basis", "def2_diffuse", "bse_2019"],
    "def2_qzv": ["def2_basis", "def2_qzvp_2003", "def2_diffuse", "bse_2019"],
    "def2_mtzv": ["def2_mtzvp", "def2_mtzvpp", "bse_2019"],
    "def2_jk": ["def2_jfit", "def2_jkfit", "def2_mtzvpp", "bse_2019"],
    "def2_rifit": ["def2_rifit_1998", "def2_rifit_2005", "def2_rifit_2007", "def2_rifit_diffuse", "bse_2019"],
    "def2_ecp": ["bse_2019"],
    "minix": ["hf3c", "bse_2019"],
}


# ------------------------------------------------------------------------------------------------- input
def fetch(path, cache):
    """Bytes of one file of the BSE repository at BSE_COMMIT (from the cache directory if present)."""
    local = os.path.join(cache, path.replace("/", "_")) if cache else None
    old_local = os.path.join(cache, os.path.basename(path)) if cache else None  # stage-1 cache layout
    for f in (local, old_local):
        if f and os.path.exists(f):
            with open(f, "rb") as h:
                return h.read()
    url = BSE_RAW + urllib.parse.quote(path)
    data = subprocess.run(["curl", "-sSfL", url], check=True, capture_output=True).stdout
    if local:
        os.makedirs(cache, exist_ok=True)
        with open(local, "wb") as f:
            f.write(data)
    return data


def fetch_json(path, cache):
    return json.loads(fetch(path, cache).decode())


def parse_tm(text):
    """Turbomole basis text -> ({element: [(L, [(exp_token, coef_token), ...]), ...]}, {element: [ecp lines]}).

    Tokens as written. Primitive lines with three columns (index, exponent, coefficient: def2-mTZVPP-RIJ) lose the
    index. The '$ecp' section is returned line by line (whitespace normalised)."""
    out, ecp, elem, expect, sect = {}, {}, None, False, None
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        t = lines[i].strip()
        i += 1
        if not t or t.startswith("#"):
            continue
        if t.startswith("$"):
            sect = t.split()[0].lower()
            if sect == "$end":
                break
            elem, expect = None, False
            continue
        if t.startswith("*"):
            elem, expect = None, True
            continue
        if expect and t[0].isalpha():
            elem = t.split()[0].lower()
            while lines[i].strip() == "" or lines[i].strip().startswith("#"):
                i += 1
            assert lines[i].strip().startswith("*"), lines[i]
            i += 1
            expect = False
            if sect == "$ecp":
                assert elem not in ecp, elem
                ecp[elem] = []
            else:
                assert elem not in out, elem
                out[elem] = []
            continue
        assert elem is not None, "data outside an element block: " + t
        if sect == "$ecp":
            ecp[elem].append(" ".join(t.split()))
            continue
        n, l = t.split()[:2]
        prims = []
        for _ in range(int(n)):
            tok = lines[i].split()
            i += 1
            if len(tok) == 3:
                tok = tok[1:]
            assert len(tok) == 2, lines[i - 1]
            prims.append((tok[0], tok[1]))
        out[elem].append((LETTERS.index(l.lower()), prims))
    return out, ecp


def parse_dat(text):
    """curcuma '$DATA' text -> {element: [(L, [(exp_token, coef_token), ...]), ...]} (all elements in NAMES)."""
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
        if len(t) == 1 and t[0].upper() in NAMES:
            elem = SYMBOLS[NAMES.index(t[0].upper())]
            out[elem] = []
            continue
        assert elem is not None, t
        prims = []
        for _ in range(int(t[1])):
            tok = lines[i].split()
            i += 1
            prims.append((tok[1], tok[2]))
        out[elem].append(("SPDFGHI".index(t[0].upper()), prims))
    return out


def num(tok):
    return float(tok.replace("D", "E").replace("d", "e"))


def etok(tok):
    """Token as embedded: Fortran 'D' exponent marker -> 'E' (the only change)."""
    return tok.replace("D", "E").replace("d", "e")


# ------------------------------------------------------------------------------------------------- BSE json
def compose_bse(table_name, cache):
    """BSE composition (basis_set_exchange/compose.py, no manipulation): {element: (shells, ecp)}.

    shells: [(L, [(exp, coef), ...])] split per general-contraction column, in component order;
    ecp: None or (ncore, [(L, [(coef, r, exp), ...])])."""
    table = fetch_json(BSE_DATA + table_name, cache)
    out = {}
    comp_cache = {}
    for z, efile in table["elements"].items():
        el = fetch_json(BSE_DATA + efile, cache)["elements"][z]
        shells, ecp = [], None
        for comp in el["components"]:
            if comp not in comp_cache:
                comp_cache[comp] = fetch_json(BSE_DATA + comp, cache)["elements"]
            cd = comp_cache[comp][z]
            for sh in cd.get("electron_shells", []):
                assert len(sh["angular_momentum"]) == 1, (table_name, z, "sp shell")
                for col in sh["coefficients"]:
                    shells.append((sh["angular_momentum"][0], list(zip(sh["exponents"], col))))
            if "ecp_potentials" in cd:
                assert ecp is None
                pots = []
                for p in cd["ecp_potentials"]:
                    assert len(p["coefficients"]) == 1 and len(p["angular_momentum"]) == 1
                    pots.append((p["angular_momentum"][0],
                                 list(zip(p["coefficients"][0], p["r_exponents"], p["gaussian_exponents"]))))
                ecp = (cd["ecp_electrons"], pots)
        out[SYMBOLS[int(z) - 1]] = (shells, ecp)
    return out


def tm_ecp_to_values(lines):
    """Turbomole ECP lines of one element -> (ncore, [(L, [(coef, r, exp)])]) with floats."""
    m = re.match(r"ncore\s*=\s*(\d+)\s+lmax\s*=\s*(\d+)$", lines[0])
    assert m, lines[0]
    ncore, lmax = int(m.group(1)), int(m.group(2))
    pots, cur = [], None
    for ln in lines[1:]:
        t = ln.split()
        if len(t) == 1 and re.fullmatch(r"[spdfghi](-[spdfghi])?", t[0]):
            L = LETTERS.index(t[0][0]) if "-" in t[0] else lmax
            if "-" not in t[0]:
                assert LETTERS.index(t[0]) == lmax
            cur = (L, [])
            pots.append(cur)
            continue
        assert len(t) == 3 and cur is not None, ln
        cur[1].append((float(t[0]), int(t[1]), float(t[2])))
    return ncore, pots


def shell_key(l, prims):
    return (l, tuple((num(e), num(c)) for e, c in prims))


def compare_json(name, kind, tm, tm_ecp, bse):
    """Check (1). Returns a dict of findings; raises on a value difference."""
    res = {"elements": len(bse), "only_tm": [], "only_bse": [], "order": [], "spelling": 0, "values": 0, "ecp": 0}
    if kind == "ecp":
        res["only_tm"] = sorted(set(tm_ecp) - set(bse), key=SYMBOLS.index)
        res["only_bse"] = sorted(set(bse) - set(tm_ecp), key=SYMBOLS.index)
        for el in tm_ecp:
            ncore, pots = tm_ecp_to_values(tm_ecp[el])
            jcore, jpots = bse[el][1]
            jv = [(L, [(float(c), int(r), float(e)) for c, r, e in p]) for L, p in jpots]
            assert ncore == jcore and pots == jv, (name, el, "ECP differs from BSE json")
            res["ecp"] += 1
            res["values"] += sum(3 * len(p) for _, p in pots)
        return res
    res["only_tm"] = sorted(set(tm) - set(bse), key=SYMBOLS.index)
    res["only_bse"] = sorted(set(bse) - set(tm), key=SYMBOLS.index)
    for el in tm:
        if el not in bse:
            continue
        a = [shell_key(l, p) for l, p in tm[el]]
        b = [(l, tuple((float(e), float(c)) for e, c in p)) for l, p in bse[el][0]]
        assert sorted(a) == sorted(b), (name, el, "shells differ from BSE json")
        if a != b:
            res["order"].append(el)
        for (l, pa), (_, pb) in zip(sorted(tm[el], key=lambda s: shell_key(*s)),
                                    sorted(bse[el][0], key=lambda s: (s[0], tuple((float(e), float(c)) for e, c in s[1])))):
            for (ea, ca), (eb, cb) in zip(pa, pb):
                res["values"] += 2
                res["spelling"] += (etok(ea) != eb) + (etok(ca) != cb)
    return res


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
    for z, el in enumerate(HNE, start=1):
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
    lines = ["", NAMES[SYMBOLS.index(el)]]
    for l, prims in shells:
        lines.append("    %s   %d" % ("SPDFGHI"[l], len(prims)))
        for k, (e, c) in enumerate(prims, start=1):
            lines.append("    %d  %20s  %20s" % (k, etok(e), etok(c)))
    return "\n".join(lines) + "\n"


def tm_block(el, name, shells):
    lines = ["*", "%s %s" % (el, name), "*"]
    for l, prims in shells:
        lines.append("    %d   %s" % (len(prims), LETTERS[l]))
        for e, c in prims:
            lines.append("    %20s  %20s" % (etok(e), etok(c)))
    return "\n".join(lines) + "\n"


def ecp_block(el, name, lines):
    return "\n".join(["*", "%s %s" % (el, name.lower()), "*"] + lines) + "\n"


def compare_repo(name, emb, ref, elements):
    """Identical doubles and shell structure (same order); returns the number of compared values."""
    count = 0
    for el in elements:
        a, b = emb[el], ref[el]
        assert [(l, len(p)) for l, p in a] == [(l, len(p)) for l, p in b], (name, el, "shell structure differs")
        for (l, pa), (_, pb) in zip(a, b):
            for (ea, ca), (eb, cb) in zip(pa, pb):
                assert num(ea) == num(eb) and num(ca) == num(cb), (name, el, ea, eb, ca, cb)
                count += 2
    return count


def differing_elements(emb, ref):
    """(elements whose shells differ as a multiset of (L, numbers), elements with the same shells in another order)."""
    values, order = [], []
    for el in SYMBOLS:
        if el in emb and el in ref:
            a = [shell_key(l, p) for l, p in emb[el]]
            b = [shell_key(l, p) for l, p in ref[el]]
            if sorted(a) != sorted(b):
                values.append(el)
            elif a != b:
                order.append(el)
    return values, order


def read_orca(orca_file):
    data, el = {}, None
    for line in open(orca_file):
        t = line.split()
        if not t:
            continue
        if t[0].startswith("New") and t[0].endswith("GTO") and len(t) >= 2:  # NewGTO, NewAuxJGTO, NewAuxCGTO ...
            el = t[1].lower()
            data[el] = []
            continue
        if el is None:
            continue
        if t[0] == "end;" or t[0] == "end":
            el = None
            continue
        if t[0].upper() in "SPDFGHI" and len(t) == 2 and t[1].isdigit():
            data[el].append(("SPDFGHI".index(t[0].upper()), []))
            continue
        if len(t) == 3 and t[0].isdigit():
            data[el][-1][1].append((t[1], t[2]))
    return data


def compare_orca(name, emb, orca_file, elements):
    """Max relative differences against an ORCA NewGTO export (shells matched by L, length and exponents)."""
    data = read_orca(orca_file)
    de = dc_norm = ac = 0.0
    nsh = 0
    missing = []
    for e in elements:
        if e not in emb:
            continue
        if e not in data:
            missing.append(e)
            continue
        key = lambda s: (s[0], len(s[1]), tuple(-num(x) for x, _ in s[1]))
        a, b = sorted(emb[e], key=key), sorted(data[e], key=key)
        if [(l, len(p)) for l, p in a] != [(l, len(p)) for l, p in b]:
            return ("structure differs for %s: %s vs ORCA %s" % (e, [(l, len(p)) for l, p in a],
                                                               [(l, len(p)) for l, p in b]))
        for (l, pa), (_, pb) in zip(a, b):
            ea, eb = [num(x) for x, _ in pa], [num(x) for x, _ in pb]
            ca, cb = [num(y) for _, y in pa], [num(y) for _, y in pb]
            de = max(de, max(abs(x - y) / abs(y) for x, y in zip(ea, eb)))
            na, nb = normalised(l, ea, ca), normalised(l, eb, cb)
            dc_norm = max(dc_norm, max(abs(x - y) / max(abs(y), 1e-300) for x, y in zip(na, nb)))
            # ORCA prints 10 decimals: absolute difference in ORCA's scale (embedded rescaled to ORCA's norm)
            scale = math.sqrt(self_overlap(l, eb, cb))
            ac = max(ac, max(abs(x * scale - y) for x, y in zip(na, cb)))
            nsh += 1
    msg = ("%d shells, max rel diff exponents %.1e, normalised coefficients %.1e (max abs %.1e in ORCA's 10-decimal"
           " scale)" % (nsh, de, dc_norm, ac))
    if missing:
        msg += "; not in ORCA: " + " ".join(missing)
    return msg


# ------------------------------------------------------------------------------------------------- C++
def cxx_pieces(text, limit=15000):
    """One raw string literal per piece, pieces split at line ends below `limit` characters (MSVC 16 KB)."""
    assert ")BASIS\"" not in text
    pieces, cur = [], ""
    for line in text.splitlines(keepends=True):
        if cur and len(cur) + len(line) > limit:
            pieces.append(cur)
            cur = ""
        cur += line
    if cur:
        pieces.append(cur)
    assert sum(len(p) for p in pieces) < 60000, "element block too large"
    return "\n".join('        R"BASIS(' + p + ')BASIS"' for p in pieces)


def citations(keys):
    """BibTeX entries copied from src/core/citation_database.cpp."""
    src = open(os.path.join(HERE, "..", "src", "core", "citation_database.cpp"), encoding="utf-8").read()
    out = []
    lines = src.splitlines()
    for key in keys:
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


def ident(name):
    return re.sub(r"[^A-Za-z0-9]", "_", name)


def notice(w, title, license_text, family_sets, cite_keys):
    w.append("/*")
    w.append(" * <Native QM: built-in basis sets, %s>" % title)
    w.append(" * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>")
    w.append(" *")
    w.append(" * GENERATED by scripts/qm_basis_export_tables.py - do not edit.")
    w.append(" * Claude Generated (Nov 2026). AI-generated, machine-tested.")
    w.append(" *")
    w.append(" * Source: Basis Set Exchange repository, https://github.com/MolSSI-BSE/basis_set_exchange")
    w.append(" *   commit " + BSE_COMMIT + " (master of 2026-08-07, downloaded 2026-11-01)")
    for s in family_sets:
        path = (BSE_DATA if s["kind"] == "orbital" and s["name"] == "MINIX" else BSE_TM) + s["src"]
        w.append(" *   %-22s <- %s  (%s)" % (s["name"], path, s["range"]))
    if any(s["name"] != "MINIX" for s in family_sets):
        w.append(" * def2 data: the BSE Turbomole sources, every exponent/coefficient token unchanged except Fortran 'D' -> 'E'")
        w.append(" *   (shell order as in the source). The BSE files state 'Data from Turbomole 7.3'; the BSE BSD-3 notice")
        w.append(" *   below is the licence the data are used under; Turbomole's own terms were NOT verified. Checked by the")
        w.append(" *   generator: the numbers equal the BSE json data (data/<name>.1.table.json composed as BSE composes it).")
    if any(s["name"] == "MINIX" for s in family_sets):
        w.append(" * MINIX: Sure & Grimme, J. Comput. Chem. 2013, 34, 1672-1685, Table 1 (MINIX composition):")
        w.append(" *     H-He, B-Ne: MINIS; Li-Be: MINIS+1(p); Na-Mg: MINIS+1(p); Al-Ar: MINIS+1(d); K-Zn: SV;")
        w.append(" *     Ga-Kr: SVP; Rb-Xe: def2-SV(P) with ECP.  Only H-Ne is built in (Rb-Xe needs the def2 ECP in the SCF).")
        w.append(" *   MINIS = BSE 'Scaled MINI' (data/huzinaga/MINI-scaled.0.json, original BSE contributor David Feller):")
        w.append(" *   exponents as in BSE; coefficients = BSE coefficients renormalised to a unit-norm contracted function,")
        w.append(" *   10 decimals (= ORCA 6.1 orca_exportbasis -b minix = qm_methods/MINIX.dat). Li/Be extra p shell")
        w.append(" *   (exponent 0.1 / 0.25, coefficient 1): from ORCA 6.1 orca_exportbasis -b minix (not in BSE).")
    w.append(" * A file found in the basis search path takes precedence over these copies (qm_basis_embedded.h).")
    w.append(" *")
    w.append(" * References (BibTeX from src/core/citation_database.cpp):")
    for entry in citations(cite_keys):
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


def render_family(fam, title, license_text, sets):
    fs = [s for s in sets if s["family"] == fam]
    w = []
    notice(w, title, license_text, fs, FAMILY_CITES[fam])
    w.append("")
    w.append("#pragma once")
    w.append("")
    w.append("// Data only; included by qm_basis_tables.h (and through it only by qm_basis_embedded.cpp).")
    w.append("namespace qmbasis {")
    w.append("namespace tables {")
    w.append("")
    for s in fs:
        idn = ident(s["name"])
        w.append("// %s: %s, %d elements (%s), %d characters" % (s["name"], s["format"], len(s["blocks"]), s["range"],
                                                               s["chars"]))
        w.append("static const char kHeader_%s[] = R\"BASIS(%s)BASIS\";" % (idn, s["header"]))
        w.append("static const char kFooter_%s[] = R\"BASIS(%s)BASIS\";" % (idn, s["footer"]))
        w.append("static const char* const kBlocks_%s[%d] = {" % (idn, s["max_z"]))
        for z in range(1, s["max_z"] + 1):
            el = SYMBOLS[z - 1]
            if el in s["blocks"]:
                w.append("    // %d %s" % (z, el))
                w.append(cxx_pieces(s["blocks"][el]) + ",")
            else:
                w.append("    nullptr,  // %d %s: not in %s" % (z, el, s["name"]))
        w.append("};")
        w.append("")
    w.append("}  // namespace tables")
    w.append("}  // namespace qmbasis")
    return "\n".join(w) + "\n"


def render_registry(license_text, sets):
    w = []
    notice(w, "registry of all built-in sets", license_text, sets,
           sorted({k for v in FAMILY_CITES.values() for k in v}, key=lambda k: list(
               dict.fromkeys(k2 for v in FAMILY_CITES.values() for k2 in v)).index(k)))
    w.append("")
    w.append("#pragma once")
    w.append("")
    w.append("// Include ONLY from qm_basis_embedded.cpp: the data headers define the text (several MB) with internal")
    w.append("// linkage, kEmbeddedBasisSets/kEmbeddedBasisSetCount are defined here (declared in qm_basis_embedded.h).")
    w.append("#include \"qm_basis_embedded.h\"")
    w.append("")
    for fam, _ in FAMILIES:
        w.append("#include \"qm_basis_data_%s.h\"" % fam)
    w.append("")
    w.append("namespace qmbasis {")
    w.append("namespace tables {")
    w.append("")
    w.append("const EmbeddedBasisSet kEmbeddedBasisSets[] = {")
    for s in sets:
        idn = ident(s["name"])
        w.append('    { "%s", "%s", kHeader_%s, kFooter_%s, kBlocks_%s, %d, "%s", "%s" },'
                 % (s["name"], s["format"], idn, idn, idn, s["max_z"], s["range"], s["ecp"]))
    w.append("};")
    w.append("const int kEmbeddedBasisSetCount = %d;" % len(sets))
    w.append("")
    w.append("}  // namespace tables")
    w.append("}  // namespace qmbasis")
    return "\n".join(w) + "\n"


def element_range(els):
    zs = sorted(SYMBOLS.index(e) + 1 for e in els)
    # compact ranges with element symbols, e.g. "H-La, Hf-Rn"
    parts, start, prev = [], zs[0], zs[0]
    for z in zs[1:] + [None]:
        if z is not None and z == prev + 1:
            prev = z
            continue
        a, b = SYMBOLS[start - 1].capitalize(), SYMBOLS[prev - 1].capitalize()
        parts.append(a if start == prev else a + "-" + b)
        if z is not None:
            start = prev = z
    return ", ".join(parts)


# ------------------------------------------------------------------------------------------------- build
def build(cache, report):
    license_text = fetch(LICENSE_PATH, cache).decode()
    ecp_ref = None
    sets = []
    for name, kind, src, fam in SETS:
        s = {"name": name, "kind": kind, "src": src, "family": fam, "ecp": ""}
        if name == "MINIX":
            data, ecp = minix_from_bse(fetch_json(BSE_DATA + src, cache)), {}
            s["bse"] = None
        else:
            data, ecp = parse_tm(bz2.decompress(fetch(BSE_TM + src, cache)).decode())
            bse = compose_bse(TABLE_NAMES.get(name, name) + ".1.table.json", cache)
            s["check"] = compare_json(name, kind, data, ecp, bse)
            s["bse_ecp"] = sorted((el for el in bse if bse[el][1] is not None), key=SYMBOLS.index)
            s["bse"] = bse
        s["data"], s["ecp_lines"] = data, ecp
        if kind == "ecp":
            ecp_ref = ecp
            s["format"] = "ecp"
            s["header"], s["footer"] = "$ecp\n", "*\n$end\n"
            s["blocks"] = {el: ecp_block(el, name, ecp[el]) for el in ecp}
        elif kind == "orbital":
            s["format"] = "data"
            s["header"] = ("!----------------------------------------------------------------------\n"
                           "! %s, built into curcuma; see qm_basis_data_%s.h for the source\n"
                           "!----------------------------------------------------------------------\n\n$DATA\n"
                           % (name, fam))
            s["footer"] = "\n$END\n"
            s["blocks"] = {el: dat_block(el, data[el]) for el in data}
        else:
            s["format"] = "turbomole"
            s["header"], s["footer"] = "$%s\n" % kind, "*\n$end\n"
            s["blocks"] = {el: tm_block(el, name, data[el]) for el in data}
        els = s["blocks"].keys()
        s["max_z"] = max(SYMBOLS.index(e) + 1 for e in els)
        s["range"] = element_range(els)
        s["chars"] = len(s["header"]) + len(s["footer"]) + sum(len(b) for b in s["blocks"].values())
        sets.append(s)

    # (2) ECP sections of the orbital .tm files == def2-ECP; json ECP components
    for s in sets:
        if s["kind"] == "orbital" and s["name"] != "MINIX":
            if s["ecp_lines"]:
                for el, lines in s["ecp_lines"].items():
                    assert lines == ecp_ref[el], (s["name"], el, "ECP section differs from def2-ECP")
                s["ecp"] = "def2-ECP"
            s["ecp_note"] = "ECP in .tm: %d elements; ECP in BSE json: %d elements" % (len(s["ecp_lines"]), len(s["bse_ecp"]))
            assert sorted(s["ecp_lines"], key=SYMBOLS.index) == s["bse_ecp"], (s["name"], "ECP element lists differ")

    # (3) stage-1 files
    nrepo = 0
    repo_diff = {}
    for s in sets:
        if s["name"] not in REPO_FILES:
            continue
        repo_text = open(os.path.join(QMDIR, REPO_FILES[s["name"]]), encoding="utf-8").read()
        ref = parse_dat(repo_text) if REPO_FILES[s["name"]].endswith(".dat") else parse_tm(repo_text)[0]
        nrepo += compare_repo(s["name"], s["data"], ref, HNE)
        repo_diff[s["name"]] = (sorted(ref, key=SYMBOLS.index), differing_elements(s["data"], ref))

    if report:
        print("BSE commit %s" % BSE_COMMIT)
        print("%-21s %-9s %-24s %6s %9s  %s" % ("set", "format", "elements", "count", "chars", "check vs BSE json"))
        for s in sets:
            c = s.get("check")
            if c is None:
                chk = "MINIX: BSE MINI-scaled renormalised + ORCA Li/Be p (equals qm_methods/MINIX.dat)"
            elif s["kind"] == "ecp":
                chk = "ECP of %d elements equal (%d values)" % (c["ecp"], c["values"])
            else:
                chk = "%d values equal; spelling differs %d; shell order differs for %d elements" % (
                    c["values"], c["spelling"], len(c["order"]))
                if c["only_tm"] or c["only_bse"]:
                    chk += "; ONLY in .tm: %s, ONLY in json: %s" % (c["only_tm"], c["only_bse"])
            print("%-21s %-9s %-24s %6d %9d  %s" % (s["name"], s["format"], s["range"], len(s["blocks"]), s["chars"], chk))
        for s in sets:
            if s.get("ecp_note"):
                print("  %-19s %s%s" % (s["name"], s["ecp_note"], "" if s["ecp"] else " (no ECP for any element)"))
        for s in sets:
            c = s.get("check")
            if c and c["order"]:
                hne = [e for e in c["order"] if e in HNE]
                print("  %-19s shell order differs from the BSE json composition (base + add files) for: %s%s" % (
                    s["name"], " ".join(c["order"]), ("  [H-Ne: %s]" % " ".join(hne)) if hne else ""))
        print("stage-1 sets equal to the files in qm_methods/ for H-Ne: %d values" % nrepo)
        for name, (els, (dval, dord)) in repo_diff.items():
            print("  %-19s file %s has %d elements (%s); numbers differ from BSE for %d%s; same numbers, other shell"
                  " order for %d%s" % (
                      name, REPO_FILES[name], len(els), element_range(els), len(dval),
                      (": " + " ".join(e.capitalize() for e in dval)) if dval else "", len(dord),
                      (": " + " ".join(e.capitalize() for e in dord)) if dord else ""))
    return license_text, sets


def outputs(license_text, sets):
    files = {}
    for fam, title in FAMILIES:
        files["qm_basis_data_%s.h" % fam] = render_family(fam, title, license_text, sets)
    files["qm_basis_tables.h"] = render_registry(license_text, sets)
    return files


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cache", help="directory with (or for) the BSE files")
    ap.add_argument("--check", action="store_true", help="fail if a checked-in header differs")
    ap.add_argument("--compare-orca", metavar="DIR", help="directory with orca_exportbasis outputs")
    a = ap.parse_args()
    license_text, sets = build(a.cache, True)
    if a.compare_orca:
        for s in sets:
            if s["name"] not in ORCA_NAMES:
                print("ORCA %-21s no ORCA 6.1 set of that name" % s["name"])
                continue
            f = os.path.join(a.compare_orca, ORCA_NAMES[s["name"]].replace("/", "_") + ".bas")
            if not os.path.exists(f):
                print("ORCA %-21s missing %s" % (s["name"], f))
                continue
            print("ORCA %-21s (%s) H-Ne: %s" % (s["name"], ORCA_NAMES[s["name"]], compare_orca(s["name"], s["data"], f, HNE)))
    files = outputs(license_text, sets)
    total = 0
    bad = []
    for fn, text in files.items():
        path = os.path.join(QMDIR, fn)
        total += len(text.encode())
        if a.check:
            old = open(path, encoding="utf-8").read() if os.path.exists(path) else ""
            if old != text:
                bad.append(fn)
        else:
            with open(path, "w", encoding="utf-8") as f:
                f.write(text)
    print("%d generated files, %d bytes" % (len(files), total))
    if a.check:
        if bad:
            print("differ from the generator output (rerun scripts/qm_basis_export_tables.py):", " ".join(bad))
            return 1
        print("all generated headers are up to date")
        return 0
    print("wrote", " ".join(sorted(files)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
