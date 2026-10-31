#!/usr/bin/env python3
"""
gpu_parity_check.py - compare what the CUDA and ROCm backends contain, from the source text only.

Claude Generated (Oct 2026). Nothing is compiled or run. It finds paths that exist on one backend
and not on the other; it says nothing about speed or correctness.

Sections
  seam     GpuScfBackend `supports*` overrides (GFN1/GFN2), CUDA: xtb_gpu_method.cpp,
           ROCm: xtb_hip_method.cpp plus the shared template in xtb_gpu_adapter.h.
           "switch" means the ROCm override is not a plain `return true;` (a runtime switch or a
           fallback decides), so presence is not the same as "on by default".
  traits   GFN-FF backend traits (static constexpr in gfnff_gpu_method.h / gfnff_hip_method.h).
  kernels  `__global__` function names (CUDA .cu/.cuh, ROCm .hip/.hiph/.h), names normalised by
           removing `_hip`/`_cuda` tokens; same function under a different name shows as a pair
           of one-sided entries, so read the list, do not count it.
  env      "CURCUMA_*" string literals in the backend files, prefixes ROCM_/CUDA_/GPU_ removed.
  params   PARAM names `rocm_*` and `cuda_*` (prefix removed); `gpu_*` PARAMs are shared.

Only seam and traits are counted: they are typed declarations, so a difference there is a real
difference. kernels, env and params are printed for reading and never counted: the two backends
name and group their kernels and switches differently (ROCm has one rocm_* switch per feature,
CUDA mostly has the feature always on), so a one-sided name is a hint, not a finding.

Known differences are kept in scripts/gpu_parity_known.json ({section: {key: reason}}); a one-sided
entry with a reason is printed as known, without one as UNEXPLAINED. {"kernel_aliases": {name: name}}
maps a kernel name of either backend onto a common name so renamed pairs are not listed.

Usage
  python3 scripts/gpu_parity_check.py [--section seam,traits,kernels,env,params]
                                      [--known FILE] [--fail-on-unexplained]
Exit status 1 with --fail-on-unexplained when an unexplained difference exists.
"""
import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EC = ROOT / "src/core/energy_calculators"

CUDA_DIRS = [EC / "ff_methods/cuda", EC / "qm_methods/cuda"]
ROCM_DIRS = [EC / "ff_methods/rocm", EC / "qm_methods/rocm"]
CUDA_EXTRA = [EC / "qm_methods/xtb_gpu_method.cpp", EC / "qm_methods/gfnff_gpu_method.cpp",
              EC / "qm_methods/gfnff_gpu_method.h"]
ROCM_EXTRA = [EC / "qm_methods/xtb_hip_method.cpp", EC / "qm_methods/gfnff_hip_method.cpp",
              EC / "qm_methods/gfnff_hip_method.h"]
SRC_EXT = {".cu", ".cuh", ".cpp", ".h", ".hip", ".hiph", ".hpp"}


def read(p):
    try:
        return Path(p).read_text(errors="replace")
    except OSError:
        return ""


def files(dirs, extra):
    out = []
    for d in dirs:
        if d.is_dir():
            out += [p for p in sorted(d.rglob("*")) if p.suffix in SRC_EXT]
    out += [p for p in extra if p.exists()]
    return out


def brace_body(text, start):
    """Text of the {...} block that opens at or after `start`, without the outer braces."""
    i = text.find("{", start)
    if i < 0:
        return ""
    depth = 0
    for j in range(i, len(text)):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                return text[i + 1:j]
    return ""


# ---------------------------------------------------------------- seam
def seam():
    base = sorted(set(re.findall(r"virtual bool (supports\w+)\(\) const", read(EC / "qm_methods/xtb_native.h"))))

    def overrides(paths):
        found = {}
        for p in paths:
            t = read(p)
            for m in re.finditer(r"bool (supports\w+)\(\) const override", t):
                body = brace_body(t, m.end()).strip()
                kind = "yes" if re.fullmatch(r"return\s+true\s*;", body) else "switch"
                # the first definition wins; a plain `true` beats a conditional one
                if m.group(1) not in found or kind == "yes":
                    found[m.group(1)] = kind
        return found

    cuda = overrides([EC / "qm_methods/xtb_gpu_method.cpp"])
    rocm = overrides([EC / "qm_methods/xtb_hip_method.cpp", EC / "qm_methods/xtb_gpu_adapter.h"])
    rows = []
    for name in base:
        c, r = cuda.get(name), rocm.get(name)
        rows.append((name, c or "-", r or "-"))
    return rows


# ---------------------------------------------------------------- traits
def traits():
    def parse(p):
        return dict(re.findall(r"static constexpr (?:bool|int) (\w+)\s*=\s*([^;]+);", read(p)))

    c = parse(EC / "qm_methods/gfnff_gpu_method.h")
    r = parse(EC / "qm_methods/gfnff_hip_method.h")
    return [(k, c.get(k, "-"), r.get(k, "-")) for k in sorted(set(c) | set(r))]


# ---------------------------------------------------------------- kernels
KERNEL = re.compile(r"__global__\s+(?:[A-Z_0-9]+\s+|__launch_bounds__\([^)]*\)\s+)*(?:static\s+|inline\s+)*void\s+(\w+)")


def norm_kernel(n, aliases):
    n = aliases.get(n, n)
    n = re.sub(r"(?:^|_)(hip|cuda)(?=_|$)", "", n)
    n = re.sub(r"^k_res_", "k_", n)   # ROCm resident-loop copies of CUDA kernels (k_res_broyden_* ~ k_broyden_*)
    return n


def kernels(aliases):
    def collect(paths):
        s = {}
        for p in paths:
            for m in KERNEL.finditer(read(p)):
                s.setdefault(norm_kernel(m.group(1), aliases), set()).add(m.group(1))
        return s

    return collect(files(CUDA_DIRS, CUDA_EXTRA)), collect(files(ROCM_DIRS, ROCM_EXTRA))


# ---------------------------------------------------------------- env, params
def norm_env(n):
    return re.sub(r"^CURCUMA_(?:ROCM_|CUDA_|GPU_)", "CURCUMA_", n)


def env():
    pat = re.compile(r'"(CURCUMA_[A-Z0-9_]+)"')

    def collect(paths):
        return {norm_env(m): m for p in paths for m in pat.findall(read(p))}

    return collect(files(CUDA_DIRS, CUDA_EXTRA)), collect(files(ROCM_DIRS, ROCM_EXTRA))


def params():
    c, r, shared = {}, {}, set()
    for p in sorted((ROOT / "src").rglob("*.h")):
        t = read(p)
        for m in re.finditer(r"PARAM\((cuda|rocm|gpu)_(\w+?)\s*,", t):
            kind, rest = m.groups()
            if kind == "cuda":
                c[rest] = f"cuda_{rest}"
            elif kind == "rocm":
                r[rest] = f"rocm_{rest}"
            else:
                shared.add(rest)
    return c, r, len(shared)


# ---------------------------------------------------------------- report
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--section", default="seam,traits,kernels,env,params")
    ap.add_argument("--known", default=str(ROOT / "scripts/gpu_parity_known.json"))
    ap.add_argument("--fail-on-unexplained", action="store_true")
    a = ap.parse_args()
    secs = set(a.section.split(","))
    known = {}
    kp = Path(a.known)
    if kp.exists():
        known = json.loads(kp.read_text())
    unexplained = 0

    counted = {"seam", "traits"}

    def tag(section, key):
        nonlocal unexplained
        why = known.get(section, {}).get(key)
        if why:
            return f"known: {why}"
        if section not in counted:
            return "info"
        unexplained += 1
        return "UNEXPLAINED"

    if "seam" in secs:
        print("## seam: GpuScfBackend supports* (GFN1/GFN2)\n")
        print("| path | CUDA | ROCm | note |\n|---|---|---|---|")
        for name, c, r in seam():
            note = ""
            if (c == "-") != (r == "-"):
                note = tag("seam", name)
            elif "switch" in (c, r) and c != r:
                note = "ROCm decided at runtime (switch or fallback)"
            print(f"| {name} | {c} | {r} | {note} |")
        print()

    if "traits" in secs:
        print("## traits: GFN-FF backend traits\n")
        print("| trait | CUDA | ROCm | note |\n|---|---|---|---|")
        for k, c, r in traits():
            note = tag("traits", k) if c != r else ""
            print(f"| {k} | {c} | {r} | {note} |")
        print()

    if "kernels" in secs:
        c, r = kernels(known.get('kernel_aliases', {}))
        oc, orr = sorted(set(c) - set(r)), sorted(set(r) - set(c))
        print(f"## kernels: `__global__` names (CUDA {len(c)}, ROCm {len(r)}, common {len(set(c) & set(r))})\n")
        print(f"CUDA only ({len(oc)}):")
        for k in oc:
            print(f"- {sorted(c[k])[0]}  [{tag('kernels', k)}]")
        print(f"\nROCm only ({len(orr)}):")
        for k in orr:
            print(f"- {sorted(r[k])[0]}  [{tag('kernels', k)}]")
        print()

    if "env" in secs:
        c, r = env()
        oc, orr = sorted(set(c) - set(r)), sorted(set(r) - set(c))
        print(f"## env: CURCUMA_* literals (CUDA {len(c)}, ROCm {len(r)}, common {len(set(c) & set(r))})\n")
        print(f"CUDA only ({len(oc)}):")
        for k in oc:
            print(f"- {c[k]}  [{tag('env', k)}]")
        print(f"\nROCm only ({len(orr)}):")
        for k in orr:
            print(f"- {r[k]}  [{tag('env', k)}]")
        print()

    if "params" in secs:
        c, r, shared = params()
        oc, orr = sorted(set(c) - set(r)), sorted(set(r) - set(c))
        print(f"## params: PARAM switches (cuda_* {len(c)}, rocm_* {len(r)}, shared gpu_* {shared})\n")
        print(f"cuda_* without a rocm_* counterpart ({len(oc)}):")
        for k in oc:
            print(f"- {c[k]}  [{tag('params', k)}]")
        print(f"\nrocm_* without a cuda_* counterpart ({len(orr)}):")
        for k in orr:
            print(f"- {r[k]}  [{tag('params', k)}]")
        print()

    print(f"UNEXPLAINED differences: {unexplained}")
    if a.fail_on_unexplained and unexplained:
        sys.exit(1)


if __name__ == "__main__":
    main()
