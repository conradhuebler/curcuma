#!/usr/bin/env python3
# Claude Generated (Oct 2026). AI-generated, machine-tested only.
"""Generate the CTest entries of the structure property tests from the structure library.

Reads test_cases/structures/manifest.json and test_cases/struct_test_policy.json, selects (structure, method) pairs per family and writes a CMake file
with one add_test per pair. Called at configure time by test_cases/CMakeLists.txt; run it by hand with --plan to see what
is selected and why a structure is left out.

    gen_struct_tests.py --out FILE.cmake [--cuda ON|OFF] [--rocm ON|OFF]   write the CMake include file
    gen_struct_tests.py --plan                                 print the selection and the reasons for every skipped structure

Charge and spin come from the manifest. When the manifest has none, the structure is used only if the policy lists it under
`assumed_neutral_closed_shell` and its electron count is even (an AI proposal that the operator has to confirm); otherwise it
is skipped with the reason. Rules of the library: test_cases/structures/README.md.
"""
import argparse
import json
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIB = os.path.join(REPO, "test_cases", "structures")
Z = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13,
     'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,
     'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
     'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49,
     'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54}


def electrons(path):
    with open(path) as f:
        lines = f.read().split("\n")
    n = int(lines[0].split()[0])
    total = 0
    for l in lines[2:2 + n]:
        e = l.split()[0]
        if e not in Z:
            return None
        total += Z[e]
    return total


def select(manifest, policy, family):
    """Yield (structure, charge, spin, note) and (structure, None, None, reason) for skipped ones."""
    fam = policy["families"][family]
    assumed = set(policy.get("assumed_neutral_closed_shell", {}).get("ids", []))
    for s in manifest["structures"]:
        sid = s["id"]
        if s["class"] not in fam["classes"]:
            continue
        if s.get("role") == "invalid-input":
            yield s, None, None, "role invalid-input"
            continue
        if not (2 <= s["natoms"] <= fam["max_atoms"]):
            continue
        ne = electrons(os.path.join(LIB, s["file"]))
        q, sp = s.get("charge"), s.get("spin")
        if q is not None and sp is not None:
            yield s, q, sp, ""
        elif sid in assumed:
            if ne is None or ne % 2 != 0:
                yield s, None, None, "listed as assumed neutral closed shell, but the electron count is odd or unknown"
            else:
                yield s, 0, 0, "assumed neutral closed shell (operator to confirm)"
        else:
            yield s, None, None, "no charge/spin in the manifest and not in the assumed list"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    ap.add_argument("--cuda", default="OFF")
    ap.add_argument("--rocm", default="OFF")
    ap.add_argument("--plan", action="store_true")
    args = ap.parse_args()
    manifest = json.load(open(os.path.join(LIB, "manifest.json")))
    policy = json.load(open(os.path.join(REPO, "test_cases", "struct_test_policy.json")))
    out = []
    n_tests = 0
    skipped = []
    for family, fam in policy["families"].items():
        for s, q, sp, note in select(manifest, policy, family):
            if q is None:
                skipped.append((s["id"], note))
                continue
            for m in fam["methods"]:
                prefix = fam.get("name_prefix", "structinv")
                name = "%s_%s_%s" % (prefix, s["id"].replace(".", "_").replace("-", "_"), m)
                cmd = ('add_test(NAME %s COMMAND ${CURCUMA_STRUCT_PYTHON} ${CURCUMA_STRUCT_REPO}/scripts/struct_invariants.py '
                       '--curcuma $<TARGET_FILE:curcuma> --structure %s --method %s --charge %d --spin %d --tests %s)'
                       % (name, s["id"], m, q, sp, fam["tests"]))
                out.append("# %s%s" % (s["id"], (": " + note) if note else ""))
                out.append(cmd)
                out.append('set_tests_properties(%s PROPERTIES TIMEOUT 300 LABELS "structinv;%s")' % (name, m))
                n_tests += 1
                for backend, flag in (("cuda", args.cuda), ("rocm", args.rocm)):
                    if flag.upper() in ("ON", "1", "TRUE") and fam.get("gpu", True):
                        gname = name.replace(prefix + "_", prefix + ("_gpu_" if backend == "cuda" else "_rocm_"), 1)
                        out.append(cmd.replace(name, gname, 1).replace("--tests %s" % fam["tests"],
                                                                      "--tests %s --gpu %s" % (fam.get("gpu_tests", "gpu"), backend)))
                        out.append('set_tests_properties(%s PROPERTIES TIMEOUT 300 RESOURCE_LOCK "gpu" LABELS "structinv;gpu;%s;%s")'
                                   % (gname, backend, m))
                        n_tests += 1
    if args.plan:
        print("%d tests selected, %d structures skipped" % (n_tests, len(skipped)))
        for sid, why in skipped:
            print("  skipped %-40s %s" % (sid, why))
        return 0
    header = ("# Claude Generated by scripts/gen_struct_tests.py at configure time; do not edit.\n"
              "# Source: test_cases/structures/manifest.json and test_cases/struct_test_policy.json. %d tests.\n" % n_tests)
    with open(args.out, "w") as f:
        f.write(header + "\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
