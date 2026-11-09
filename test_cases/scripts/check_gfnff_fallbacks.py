#!/usr/bin/env python3
"""Fallback literals of gfnff_method.cpp against the gfnff PARAM defaults.

GFNFF::GFNFF(const json&) reads many keys with value("key", literal). A literal that differs
from the PARAM default in gfnff.h is a second, silent default (docs/debt_2026-11/PARAMS.md H6).
Exit 1 lists every registered key whose literal fallback differs from its PARAM default.
Claude Generated (Nov 2026).
"""
import re
import sys
from pathlib import Path

root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[2]
ff = root / "src/core/energy_calculators/ff_methods"

params = {}
in_block = False
for line in (ff / "gfnff.h").read_text(errors="ignore").splitlines():
    if "BEGIN_PARAMETER_DEFINITION(gfnff)" in line:
        in_block = True
    elif "END_PARAMETER_DEFINITION" in line:
        in_block = False
    m = re.match(r'\s*PARAM\(\s*(\w+)\s*,\s*(\w+)\s*,\s*(.+?)\s*,\s*"', line)
    if in_block and m:
        params[m.group(1)] = (m.group(2), m.group(3))


def literal(text, ptype):
    """Python value of a C++ literal argument, None when it is not a plain literal."""
    t = text.strip()
    t = re.sub(r'^std::string\(\s*("(?:[^"\\]|\\.)*")\s*\)$', r"\1", t)
    if t in ("true", "false"):
        return t == "true"
    if re.fullmatch(r'"(?:[^"\\]|\\.)*"', t):
        return t[1:-1]
    try:
        return float(t.rstrip("fF"))
    except ValueError:
        return None


src = (ff / "gfnff_method.cpp").read_text(errors="ignore")
call = re.compile(r'\.value\(\s*"(\w+)"\s*,\s*((?:[^()]|\((?:[^()]|\([^()]*\))*\))*?)\)', re.S)
bad, checked = [], 0
for m in call.finditer(src):
    key, arg = m.group(1), m.group(2)
    if key not in params:
        continue
    ptype, default = params[key]
    a, b = literal(arg, ptype), literal(default, ptype)
    if a is None or b is None:
        continue
    checked += 1
    if isinstance(a, float) and isinstance(b, float):
        same = a == b
    else:
        same = a == b
    if not same:
        line = src.count("\n", 0, m.start()) + 1
        bad.append(f"gfnff_method.cpp:{line}: {key}: fallback {arg.strip()} but PARAM default {default}")

print(f"checked {checked} literal fallbacks of registered keys")
for b in bad:
    print("MISMATCH", b)
sys.exit(1 if bad else 0)
