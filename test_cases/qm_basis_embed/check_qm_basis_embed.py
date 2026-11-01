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
  missing   an unknown orbital basis (def2-TZVP) or fitting basis must fail and name the search path.

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
        for qm, what in (([ "basis=def2-TZVP"], "def2-TZVP"), (["ri=jk", "ri_aux=def2-universal-XFIT"], "def2-universal-XFIT")):
            rc, out = run([a.dump, XYZ, "hf"] + qm_args(qm), env_with(CURCUMA_QM_BASIS=tmp))
            if rc == 0:
                fail("unknown basis %s did not fail" % what, out)
            if what not in out or "CURCUMA_QM_BASIS" not in out or "built in" not in out:
                fail("error message for %s does not name the search path and the built-in sets" % what, out)
            print(what, "->", [l for l in out.splitlines() if what in l][0].strip())
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True, choices=["bitwise", "empty", "priority", "missing"])
    ap.add_argument("--dump", help="dump_qm_energy binary")
    ap.add_argument("--curcuma", help="curcuma binary")
    ap.add_argument("--method", default="hf")
    ap.add_argument("--qm", action="append", default=[], help="KEY=VALUE for the qm module (dump_qm_energy)")
    ap.add_argument("--extra", default="", help="extra curcuma CLI arguments (case empty)")
    ap.add_argument("--cite", default="", help="regex that must appear in the citation block (case empty)")
    a = ap.parse_args()
    {"bitwise": case_bitwise, "empty": case_empty, "priority": case_priority, "missing": case_missing}[a.case](a)
    print("PASS")


if __name__ == "__main__":
    main()
