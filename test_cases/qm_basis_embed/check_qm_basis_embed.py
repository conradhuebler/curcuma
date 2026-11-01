#!/usr/bin/env python3
"""Built-in (embedded) basis sets of the native QM engine vs the basis files (ctest label qm_basis).

Cases (--case):
  bitwise   run dump_qm_energy twice: (1) the basis files copied into a temporary directory and found through
            $CURCUMA_QM_BASIS, (2) CURCUMA_QM_BASIS_EMBEDDED_ONLY=1 (no file lookup at all). The 17-digit energies
            must be identical, and the verbosity-2 line must name the temporary file in (1) and 'embedded' in (2).
  empty     `curcuma -sp` from an empty directory with every basis variable unset and embedded-only mode on; the
            printed energy must equal the run with the files (all printed digits).
  priority  a copy of def2-SVP.dat with one exponent changed in $CURCUMA_QM_BASIS must win over the embedded copy
            (energy changes, the log names the file).
  missing   an unknown orbital basis (def2-XZVP) or fitting basis must fail and name the search path.
  (stage 2, Nov 2026)
  parse     every built-in set: the numbers the engine's parser reads from the full text (dump_embedded_basis
            --values; '$DATA' parser keeps H-Ne, Turbomole parser all elements) equal the tokens of the text read as
            doubles, every number of the text (H-Ne for '$DATA') is read, the parse of the full text succeeds; prints
            the parse time of the full text and of its H-Ne part.
  sets      for one molecule: each case of STAGE2_CASES run with the basis files written from the built-in text into
            a temporary $CURCUMA_QM_BASIS and embedded-only; 17-digit energies must be identical.
  orca      for one molecule: each case of reference_orca.json embedded-only against ORCA 6.1, |dE| <= --tol (RI-MP2:
            2e-8, ORCA's RI-MP2 is ~1e-8 off) and, for RI-MP2, against PySCF DFMP2 at 1e-9.
  ecp       def2-ECP is built in for Rb-Rn as raw text (ncore/lmax line), has no parser and no user in the engine.

Claude Generated (Nov 2026). GPL-3.0.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, "..", ".."))
QMDIR = os.path.join(SRC, "src", "core", "energy_calculators", "qm_methods")
XYZ = os.path.join(SRC, "test_cases", "qm_1e", "H2O.xyz")
MOLS = {"H2O": XYZ, "NH3": os.path.join(SRC, "test_cases", "qm_1e", "NH3.xyz"),
        "BeH2": os.path.join(SRC, "test_cases", "qm_1e", "BeH2.xyz")}
# stage-2 bitwise cases: (method, qm settings); every set the engine can run (s, p, d orbital sets; any fitting set)
STAGE2_CASES = [("hf", ["basis=def2-SVPD"]), ("hf", ["basis=def2-mTZVP"]), ("hf", ["basis=def2-mTZVPP"]), ("hf", ["basis=def2-SV(P)", "ri=jk", "ri_aux=def2-SV(P)-JKFIT"]),
                ("hf", ["ri=j", "ri_aux=def2-mTZVPP-RIJ"]), ("hf", ["ri=jk", "ri_aux=def2-TZVP-RIFIT"])] + [
    ("ri-mp2", ["basis=%s" % ("def2-SVPD" if "VPD" in x or "VPPD" in x else "def2-SVP"), "ri_mp2_aux=def2-%s-RIFIT" % x])
    for x in ("SVP", "SV(P)", "SVPD", "TZVP", "TZVPP", "TZVPD", "TZVPPD", "QZVP", "QZVPP", "QZVPPD")]
# BeH2: the TZVP family has at most d shells for H-Be, so it runs there
BEH2_CASES = [("hf", ["basis=" + b]) for b in ("def2-TZVP", "def2-TZVPP", "def2-TZVPD", "def2-TZVPPD")] + [
    ("ri-mp2", ["basis=def2-TZVP", "ri_mp2_aux=def2-TZVP-RIFIT"])]
BASIS_VARS = ("CURCUMA_QM_BASIS", "CURCUMA_DFT_BASIS", "CURCUMA_DATA", "CURCUMA_QM_BASIS_EMBEDDED_ONLY")
FILES = ["def2-SVP.dat", "def2-SV(P).dat", "MINIX.dat", "def2-universal-JFIT.turbomole",
         "def2-universal-JKFIT.turbomole"]


def env_with(**kw):
    env = {k: v for k, v in os.environ.items() if k not in BASIS_VARS}
    env.update(kw)
    return env


def run(cmd, env, cwd=None):
    p = subprocess.run(cmd, env=env, cwd=cwd, capture_output=True, text=True, timeout=1200)
    return p.returncode, p.stdout + p.stderr


def fail(msg, out=""):
    print("FAIL:", msg)
    if out:
        print(out[-4000:])
    sys.exit(1)


def energy17(out):
    m = re.search(r"^energy (\S+)$", out, re.M)
    return m.group(1) if m else None


def qm_args(pairs):
    args = []
    for p in pairs:
        k, v = p.split("=", 1)
        args += ["--qm", k, v]
    return args


def case_bitwise(a):
    tmp = tempfile.mkdtemp(prefix="qm_basis_files_")
    try:
        for f in FILES:
            shutil.copy(os.path.join(QMDIR, f), tmp)
        cmd = [a.dump, XYZ, a.method] + qm_args(a.qm) + ["--verbosity", "2"]
        rc1, out1 = run(cmd, env_with(CURCUMA_QM_BASIS=tmp))
        rc2, out2 = run(cmd, env_with(CURCUMA_QM_BASIS_EMBEDDED_ONLY="1"))
        e1, e2 = energy17(out1), energy17(out2)
        if rc1 or e1 is None:
            fail("file run failed", out1)
        if rc2 or e2 is None:
            fail("embedded run failed", out2)
        if not re.search(r"basis \S+: file " + re.escape(tmp), out1):
            fail("file run did not report the temporary basis file", out1)
        if re.search(r"basis \S+: embedded", out1):
            fail("file run used an embedded copy", out1)
        if not re.search(r"basis \S+: embedded copy", out2):
            fail("embedded run did not report the embedded copy", out2)
        if re.search(r"basis \S+: file ", out2):
            fail("embedded-only run read a file", out2)
        print("file     ", e1)
        print("embedded ", e2)
        if e1 != e2:
            fail("energies differ")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def case_empty(a):
    tmp = tempfile.mkdtemp(prefix="qm_basis_empty_")
    files = tempfile.mkdtemp(prefix="qm_basis_files_")
    try:
        xyz = os.path.join(files, "H2O.xyz")
        shutil.copy(XYZ, xyz)
        for f in FILES:
            shutil.copy(os.path.join(QMDIR, f), files)
        cmd = [a.curcuma, "-sp", xyz, "-method", a.method, "-no_bmt", "-verbosity", "2"] + a.extra.split()
        rc1, out1 = run(cmd, env_with(CURCUMA_QM_BASIS=files), cwd=files)
        rc2, out2 = run(cmd, env_with(CURCUMA_QM_BASIS_EMBEDDED_ONLY="1"), cwd=tmp)
        pat = r"Single Point Energy = (-?\d+\.\d+) Eh"
        m1, m2 = re.search(pat, out1), re.search(pat, out2)
        if rc1 or not m1:
            fail("file run failed", out1)
        if rc2 or not m2:
            fail("embedded run from an empty directory failed", out2)
        if not re.search(r"basis \S+: embedded copy", out2) or re.search(r"basis \S+: file ", out2):
            fail("empty-directory run did not use the embedded copy only", out2)
        left = [f for f in os.listdir(tmp) if f != "curcuma_citations.bib"]
        print("file     ", m1.group(1))
        print("embedded ", m2.group(1), "(run in", tmp + ", files left:", left, ")")
        if m1.group(1) != m2.group(1):
            fail("energies differ")
        if a.cite and not re.search(a.cite, out2):
            fail("citation '%s' missing in the embedded run" % a.cite, out2)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
        shutil.rmtree(files, ignore_errors=True)


def case_priority(a):
    tmp = tempfile.mkdtemp(prefix="qm_basis_mod_")
    try:
        text = open(os.path.join(QMDIR, "def2-SVP.dat"), encoding="utf-8").read()
        mod = text.replace("13.0107010", "14.0107010", 1)  # first s exponent of H
        assert mod != text
        with open(os.path.join(tmp, "def2-SVP.dat"), "w", encoding="utf-8") as f:
            f.write(mod)
        cmd = [a.dump, XYZ, "hf", "--verbosity", "2"]
        rc1, out1 = run(cmd, env_with(CURCUMA_QM_BASIS=tmp))
        rc2, out2 = run(cmd, env_with(CURCUMA_QM_BASIS_EMBEDDED_ONLY="1"))
        e1, e2 = energy17(out1), energy17(out2)
        if rc1 or rc2 or e1 is None or e2 is None:
            fail("runs failed", out1 + out2)
        if not re.search(r"basis def2-SVP: file " + re.escape(tmp), out1):
            fail("the modified file was not used", out1)
        print("modified file", e1)
        print("embedded     ", e2)
        if abs(float(e1) - float(e2)) < 1e-6:
            fail("the modified file did not change the energy: the embedded copy won")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def case_missing(a):
    tmp = tempfile.mkdtemp(prefix="qm_basis_none_")
    try:
        for qm, what in (([ "basis=def2-XZVP"], "def2-XZVP"), (["ri=jk", "ri_aux=def2-universal-XFIT"], "def2-universal-XFIT")):
            rc, out = run([a.dump, XYZ, "hf"] + qm_args(qm), env_with(CURCUMA_QM_BASIS=tmp))
            if rc == 0:
                fail("unknown basis %s did not fail" % what, out)
            if what not in out or "CURCUMA_QM_BASIS" not in out or "built in" not in out:
                fail("error message for %s does not name the search path and the built-in sets" % what, out)
            print(what, "->", [l for l in out.splitlines() if what in l][0].strip())
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def dump(a, *args):
    p = subprocess.run([a.dump_basis] + list(args), capture_output=True, text=True, timeout=600)
    return p.returncode, p.stdout, p.stderr


def element_tables():
    """(BSE upper-case names, lower-case symbols) for Z = 1..103, read from qm_elements.h (the table the engine uses)."""
    here = os.path.dirname(os.path.abspath(__file__))
    hdr = os.path.join(here, "..", "..", "src", "core", "energy_calculators", "qm_methods", "qm_elements.h")
    txt = open(hdr).read()
    blocks = re.findall(r"\{ nullptr,(.*?)\};", txt, re.S)
    syms = [x.lower() for x in re.findall(r'"([A-Za-z]+)"', blocks[0])]
    names = re.findall(r'"([A-Z]+)"', blocks[1])
    return names, syms


def text_numbers(text, fmt):
    """{(element symbol lower, L letter, iprim within shell): [(exp, coef)]} read from the text as Python doubles,
    in the order of the text (shell index implied by order)."""
    names = {}
    out = []
    lines = text.splitlines()
    if fmt == "data":
        hne, sym = element_tables()
        dropped = set()
        i = lines.index("$DATA") + 1
        el = None
        while i < len(lines):
            t = lines[i].split()
            i += 1
            if not t or t[0].startswith("!"):
                continue
            if t[0] == "$END":
                break
            if len(t) == 1:
                el = sym[hne.index(t[0])] if t[0] in hne else None
                continue
            if t[0].upper() not in "SPDFG":
                # shell above g (h: def2-QZVP family, Ce-Yb): the parser drops the whole element
                out[:] = [x for x in out if x[0] != el]
                dropped.add(el)
                el = None
            for k in range(int(t[1])):
                tok = lines[i].split()
                i += 1
                if el:
                    out.append((el, t[0], k, float(tok[1]), float(tok[2])))
    else:
        el, expect, after_header = None, False, False
        i = 0
        while i < len(lines):
            t = lines[i].strip()
            i += 1
            if not t or t.startswith("$"):
                continue
            if t == "*":
                expect, after_header = not after_header, False  # "*" before a header opens, the one after closes
                continue
            if expect:
                el, expect, after_header = t.split()[0].lower(), False, True
                continue
            n, l = t.split()
            for k in range(int(n)):
                tok = lines[i].split()
                i += 1
                out.append((el, l.upper(), k, float(tok[0]), float(tok[1])))
    return out


def case_parse(a):
    rc, listing, err = dump(a, "--list")
    if rc:
        fail("--list failed", err)
    nsets = 0
    for row in listing.strip().splitlines():
        name, fmt, elements, nel, nchar, ecp = (row.split("\t") + [""])[:6]
        rc, text, err = dump(a, "--text", name)
        if rc:
            fail("--text %s failed" % name, err)
        if fmt == "ecp":
            print("%-21s ecp       %-26s %3s elements, %7s characters: raw text only" % (name, elements, nel, nchar))
            continue
        rc, vals, err = dump(a, "--values", name)
        if rc:
            fail("parse of the full text of %s failed" % name, err)
        ref = text_numbers(text, fmt)
        got = []
        for ln in vals.splitlines():
            el, l, k, e, c = ln.split()
            got.append((el.lower(), l, int(k), float(e), float(c)))
        # the parsers return a map keyed by element (sorted); compare per element in text order
        for coll in (ref, got):
            coll.sort(key=lambda x: x[0])
        def per_el(xs):
            d = {}
            for x in xs:
                d.setdefault(x[0], []).append(x[1:])
            return d
        r, g = per_el(ref), per_el(got)
        if r.keys() != g.keys():
            fail("%s: elements differ: text %s, parser %s" % (name, sorted(r), sorted(g)))
        for el in r:
            if r[el] != g[el]:
                bad = [(x, y) for x, y in zip(r[el], g[el]) if x != y][:3]
                fail("%s %s: parsed numbers differ from the text: %s (counts %d/%d)" % (name, el, bad, len(r[el]), len(g[el])))
        rc, tm, err = dump(a, "--time", name, "20")
        print("%-21s %-9s %-26s %3s elements, %7s characters: %5d primitives of %2d elements read exactly; %s" % (
            name, fmt, elements, nel, nchar, len(ref), len(r), tm.strip().split(" ", 1)[1] if not rc else "time n/a"))
        nsets += 1
    if nsets < 25:
        fail("only %d parsable sets listed" % nsets)


def write_all(a, d):
    rc, listing, err = dump(a, "--list")
    for row in listing.strip().splitlines():
        name, fmt = row.split("\t")[:2]
        if fmt in ("data", "turbomole"):
            rc, _, err = dump(a, "--write", name, d)
            if rc:
                fail("--write %s failed" % name, err)


def case_sets(a):
    xyz = MOLS[a.molecule]
    tmp = tempfile.mkdtemp(prefix="qm_basis_written_")
    try:
        write_all(a, tmp)
        for method, qm in (BEH2_CASES if a.molecule == "BeH2" else STAGE2_CASES):
            cmd = [a.dump, xyz, method] + qm_args(qm) + ["--verbosity", "2"]
            rc1, out1 = run(cmd, env_with(CURCUMA_QM_BASIS=tmp))
            rc2, out2 = run(cmd, env_with(CURCUMA_QM_BASIS_EMBEDDED_ONLY="1"))
            e1, e2 = energy17(out1), energy17(out2)
            if rc1 or e1 is None:
                fail("file run failed: %s %s" % (method, qm), out1)
            if rc2 or e2 is None:
                fail("embedded run failed: %s %s" % (method, qm), out2)
            if re.search(r"basis \S+: embedded", out1) or not re.search(r"basis \S+: file " + re.escape(tmp), out1):
                fail("file run did not read the written files only: %s %s" % (method, qm), out1)
            if re.search(r"basis \S+: file ", out2):
                fail("embedded-only run read a file", out2)
            print("%-4s %-7s %-55s file %s embedded %s" % (a.molecule, method, " ".join(qm), e1, e2))
            if e1 != e2:
                fail("energies differ")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def case_orca(a):
    import json
    ref = json.load(open(os.path.join(HERE, "reference_orca.json")))[a.molecule]
    xyz = MOLS[a.molecule]
    worst = 0.0
    for case, r in ref.items():
        if case == "xyz":
            continue
        qm = ["%s=%s" % kv for kv in r["curcuma"].items()] + ["scf_threshold=1e-12"]
        rc, out = run([a.dump, xyz, r["method"]] + qm_args(qm), env_with(CURCUMA_QM_BASIS_EMBEDDED_ONLY="1"))
        e = energy17(out)
        if rc or e is None:
            fail("run failed: " + case, out)
        d = float(e) - r["orca"]
        worst = max(worst, abs(d))
        # ORCA's RI-MP2 is ~1e-8 Eh off curcuma and PySCF, which agree to ~1e-11 (reference_orca.json _comment)
        tol = a.tol if r["method"] != "ri-mp2" else max(a.tol, 2e-8)
        extra = ""
        if "pyscf" in r:
            dp = float(e) - r["pyscf"]
            extra = "  PySCF diff %+.2e" % dp
            if abs(dp) > 1e-9:
                fail("%s: |dE| vs PySCF = %.2e > 1e-9" % (case, abs(dp)))
        print("%-4s %-40s curcuma %.12f ORCA %.12f diff %+.2e%s" % (a.molecule, case, float(e), r["orca"], d, extra))
        if abs(d) > tol:
            fail("%s: |dE| = %.2e > %.1e" % (case, abs(d), tol))
    print("max |dE| %.2e Eh" % worst)


def case_ecp(a):
    rc, listing, err = dump(a, "--list")
    row = [r for r in listing.splitlines() if r.startswith("def2-ECP\t")]
    if not row or row[0].split("\t")[1] != "ecp" or row[0].split("\t")[3] != "50":
        fail("def2-ECP not listed with 50 elements", listing)
    for z, ncore in ((37, 28), (53, 28), (58, 28), (86, 60)):
        rc, raw, err = dump(a, "--raw", "def2-ECP", str(z))
        if rc or not re.search(r"ncore\s*=\s*%d\s+lmax" % ncore, raw):
            fail("def2-ECP Z=%d: raw block missing or ncore != %d" % (z, ncore), raw + err)
        print("def2-ECP Z=%d:" % z, raw.splitlines()[1], "|", raw.splitlines()[3])
    for z in (1, 10, 36):
        rc, raw, err = dump(a, "--raw", "def2-ECP", str(z))
        if rc == 0:
            fail("def2-ECP has a block for Z=%d" % z)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True,
                    choices=["bitwise", "empty", "priority", "missing", "parse", "sets", "orca", "ecp"])
    ap.add_argument("--dump-basis", help="dump_embedded_basis binary")
    ap.add_argument("--molecule", default="H2O", choices=sorted(MOLS))
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--dump", help="dump_qm_energy binary")
    ap.add_argument("--curcuma", help="curcuma binary")
    ap.add_argument("--method", default="hf")
    ap.add_argument("--qm", action="append", default=[], help="KEY=VALUE for the qm module (dump_qm_energy)")
    ap.add_argument("--extra", default="", help="extra curcuma CLI arguments (case empty)")
    ap.add_argument("--cite", default="", help="regex that must appear in the citation block (case empty)")
    a = ap.parse_args()
    {"bitwise": case_bitwise, "empty": case_empty, "priority": case_priority, "missing": case_missing,
     "parse": case_parse, "sets": case_sets, "orca": case_orca, "ecp": case_ecp}[a.case](a)
    print("PASS")


if __name__ == "__main__":
    main()
