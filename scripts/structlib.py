#!/usr/bin/env python3
"""Structure library for the curcuma test suite (test_cases/structures).

Claude Generated (Oct 2026). Rules: test_cases/structures/README.md.

  structlib.py check            validate manifest and files (exit 1 on errors; used by the pre-commit hook)
  structlib.py report           counts per class and provenance, share of unknown provenance
  structlib.py usage [ID ...]   which tests use which structure (derived from legacy paths and CMake/script text)
  structlib.py migrate TESTDIR  replace the local structure copies of a test directory by structures.txt (byte-identical)
  structlib.py retire [PATH ...] remove byte-identical legacy copies (all that remain without PATH)
  structlib.py path ID ...      absolute path of the library file(s)
  structlib.py stage DEST N=ID  copy structures into DEST under the names N (for scripts)
  structlib.py add FILE ...     register a new structure (see --help); rewrites the comment line with id/charge/spin/level/source

Standard library only.
"""
import argparse
import collections
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys

ROOT = subprocess.check_output(["git", "-C", os.path.dirname(os.path.abspath(__file__)), "rev-parse", "--show-toplevel"], text=True).strip()
LIB = os.path.join(ROOT, "test_cases", "structures")
MANIFEST = os.path.join(LIB, "manifest.json")

CLASSES = ["atoms", "small", "medium", "large", "clusters", "metals", "bulk", "ensembles", "cg"]
KINDS = ["optimized", "program_output", "literature", "database", "experimental", "constructed", "unknown"]
ROLES = ["equilibrium", "non-equilibrium", "transition-state", "ensemble", "stress", "invalid-input", "unspecified"]
ID_RE = re.compile(r"^[a-z0-9][a-z0-9_.+-]*$")
MAX_FILE = 1_500_000      # bytes per structure file
FORBIDDEN_EXT = {".out", ".log", ".bib", ".csv", ".png", ".txt"}
ELEMENTS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
            "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
            "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr").split()
ELEMSET = set(ELEMENTS)
TRANSITION = set("Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os Ir Pt Au Hg La".split())

# required provenance fields per kind (new entries); legacy entries may be partial but must name the kind
REQUIRED = {
    "optimized": ["program", "method"],
    "program_output": ["program"],
    "literature": ["reference"],
    "database": ["reference"],
    "experimental": ["reference"],
    "constructed": ["description"],
    "unknown": [],
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def load_manifest():
    with open(MANIFEST, encoding="utf-8") as f:
        return json.load(f)


def save_manifest(m):
    m["structures"].sort(key=lambda e: e["id"])
    with open(MANIFEST, "w", encoding="utf-8") as f:
        json.dump(m, f, indent=1, ensure_ascii=False)
        f.write("\n")


def parse_xyz(path):
    """Return (frames, natoms, symbols_of_first_frame, coords_first, comments) or raise ValueError."""
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.read().split("\n")
    i, frames, first, comments = 0, 0, None, []
    natoms0 = None
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            n = int(lines[i].split()[0])
        except ValueError:
            raise ValueError(f"line {i + 1}: atom count expected")
        if i + 1 + n >= len(lines) + 0 and i + 2 + n > len(lines):
            raise ValueError(f"line {i + 1}: truncated frame")
        comments.append(lines[i + 1] if i + 1 < len(lines) else "")
        sym, xyz = [], []
        for k in range(n):
            p = lines[i + 2 + k].split()
            if len(p) < 4:
                raise ValueError(f"line {i + 3 + k}: atom line expected")
            sym.append(p[0])
            xyz.append((float(p[1]), float(p[2]), float(p[3])))
        if natoms0 is None:
            natoms0, first = n, (sym, xyz)
        elif n != natoms0:
            raise ValueError(f"frame {frames + 1} has {n} atoms, first frame {natoms0}")
        frames += 1
        i += 2 + n
    if frames == 0:
        raise ValueError("no frame")
    return frames, natoms0, first[0], first[1], comments


def hill_formula(symbols):
    c = collections.Counter(s if not s.isdigit() else f"#{s}" for s in symbols)
    order = []
    if "C" in c:
        order = ["C"] + (["H"] if "H" in c else [])
        order += sorted(k for k in c if k not in ("C", "H"))
    else:
        order = sorted(c)
    return "".join(k + (str(c[k]) if c[k] > 1 else "") for k in order)


def min_neighbour(coords):
    best = []
    for i, a in enumerate(coords):
        d = min((math.dist(a, b) for j, b in enumerate(coords) if j != i), default=None)
        if d is not None:
            best.append(d)
    return (min(best), max(best)) if best else (None, None)


def check(args):
    errors, warnings = [], []
    if not os.path.exists(MANIFEST):
        print("ERROR   no manifest at", os.path.relpath(MANIFEST, ROOT))
        return 1
    try:
        m = load_manifest()
    except json.JSONDecodeError as e:
        print("ERROR   manifest.json is not valid JSON:", e)
        return 1
    if m.get("version") != 1:
        errors.append("manifest version must be 1")
    entries = m.get("structures", [])
    seen_ids, seen_sha, listed = {}, {}, set()
    for e in entries:
        i = e.get("id", "<no id>")
        where = f"[{i}]"
        if not ID_RE.match(i or ""):
            errors.append(f"{where} id must match {ID_RE.pattern}")
        if i in seen_ids:
            errors.append(f"{where} duplicate id")
        seen_ids[i] = e
        cls = e.get("class")
        if cls not in CLASSES:
            errors.append(f"{where} class {cls!r} not in {CLASSES}")
        rel = e.get("file", "")
        path = os.path.join(LIB, rel)
        listed.add(os.path.normpath(rel))
        if os.path.dirname(rel) != cls:
            errors.append(f"{where} file {rel} is not in the directory of its class {cls}")
        if os.path.splitext(rel)[1] in FORBIDDEN_EXT:
            errors.append(f"{where} forbidden file type {rel}")
        if not os.path.exists(path):
            errors.append(f"{where} file missing: {rel}")
            continue
        size = os.path.getsize(path)
        if size > MAX_FILE:
            (warnings if e.get("legacy") else errors).append(f"{where} {size} bytes exceeds {MAX_FILE}; large structures belong to a fetched set")
        h = sha256(path)
        if h != e.get("sha256"):
            errors.append(f"{where} sha256 differs from the manifest (a changed geometry needs a new id)")
        if h in seen_sha:
            errors.append(f"{where} same content as {seen_sha[h]}")
        seen_sha[h] = i
        fmt = os.path.splitext(rel)[1].lstrip(".")
        if fmt == "xyz":
            try:
                frames, n, sym, xyz, comments = parse_xyz(path)
            except ValueError as ex:
                errors.append(f"{where} xyz parse error: {ex}")
                continue
            if n != e.get("natoms"):
                errors.append(f"{where} natoms {n} in file, {e.get('natoms')} in manifest")
            if frames != e.get("frames"):
                errors.append(f"{where} frames {frames} in file, {e.get('frames')} in manifest")
            if frames > 1 and cls != "ensembles":
                errors.append(f"{where} multi-frame file must be in class ensembles")
            bad = [s for s in set(sym) if s not in ELEMSET and not s.isdigit()]
            if bad and e.get("role") != "invalid-input":
                errors.append(f"{where} unknown element symbols {bad}")
            if hill_formula(sym) != e.get("formula"):
                errors.append(f"{where} formula {hill_formula(sym)} in file, {e.get('formula')} in manifest")
            if any(not all(math.isfinite(v) for v in p) for p in xyz):
                errors.append(f"{where} non-finite coordinate")
            lo, hi = min_neighbour(xyz)
            if lo is not None and (lo < 0.5 or hi > 6.0):
                warnings.append(f"{where} nearest-neighbour distances {lo:.2f}..{hi:.2f} A: units or fragment far away?")
            if not e.get("legacy"):
                c1 = comments[0]
                if f"id={i}" not in c1:
                    errors.append(f"{where} comment line must carry id={i} (new structures)")
        p = e.get("provenance") or {}
        k = p.get("kind")
        if k not in KINDS:
            errors.append(f"{where} provenance.kind {k!r} not in {KINDS}")
        else:
            if k == "unknown" and not e.get("legacy"):
                errors.append(f"{where} provenance unknown is only allowed for legacy entries")
            if not e.get("legacy"):
                for fld in REQUIRED[k]:
                    if not p.get(fld):
                        errors.append(f"{where} provenance.{fld} is required for kind {k}")
        if not e.get("legacy"):
            for fld in ("charge", "spin"):
                if not isinstance(e.get(fld), int):
                    errors.append(f"{where} {fld} must be an integer for new structures")
            if e.get("role") not in ROLES:
                errors.append(f"{where} role {e.get('role')!r} not in {ROLES}")
        if e.get("needs_name"):
            warnings.append(f"{where} generic legacy name, give it a meaningful id")
        for lp in e.get("legacy_paths", []):
            lpath = os.path.join(ROOT, lp)
            if not os.path.exists(lpath):
                errors.append(f"{where} legacy path no longer exists: {lp} (remove it from legacy_paths)")
            elif sha256(lpath) != h:
                errors.append(f"{where} legacy path {lp} differs from the library file (migration must be byte-identical)")
    # structures.txt lists in the tests
    ids = {e["id"] for e in entries}
    for dp, dn, fn in os.walk(os.path.join(ROOT, "test_cases")):
        if "structures.txt" in fn and "GMTKN55-testset" not in dp and os.sep + "release" not in dp:
            lst = os.path.join(dp, "structures.txt")
            rel = os.path.relpath(lst, ROOT)
            names = set()
            for n, line in enumerate(open(lst, encoding="utf-8"), 1):
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                mm = re.match(r"^(\S+)(?:\s+as\s+(\S+))?$", s)
                if not mm:
                    errors.append(f"{rel}:{n}: expected '<id>' or '<id> as <file name>'")
                    continue
                sid, name = mm.group(1), mm.group(2) or (mm.group(1) + ".xyz")
                if sid not in ids:
                    errors.append(f"{rel}:{n}: structure {sid} is not in the manifest")
                if name in names:
                    errors.append(f"{rel}:{n}: file name {name} listed twice")
                names.add(name)
                if os.path.exists(os.path.join(dp, name)):
                    errors.append(f"{rel}:{n}: {name} also exists as a local file in the test directory")
    # stray files
    for dp, dn, fn in os.walk(LIB):
        for f in fn:
            rel = os.path.normpath(os.path.relpath(os.path.join(dp, f), LIB))
            if rel in ("manifest.json", "README.md", "structures.cmake") or rel in listed:
                continue
            if re.search(r"\.(topo|param)\.json$", rel):
                warnings.append(f"cache file next to a library structure (delete it before re-measuring): {rel}")
                continue
            errors.append(f"file not in manifest: {rel}")
    for msg in errors:
        print("ERROR  ", msg)
    for msg in warnings:
        print("WARNING", msg)
    print(f"structlib check: {len(entries)} structures, {len(errors)} error(s), {len(warnings)} warning(s)")
    return 1 if errors else 0


def report(args):
    m = load_manifest()
    E = m["structures"]
    print(f"structures: {len(E)} (legacy {sum(1 for e in E if e.get('legacy'))}), "
          f"fetched sets: {len(m.get('sets', []))}")
    print("by class:", dict(collections.Counter(e["class"] for e in E)))
    kinds = collections.Counter(e["provenance"]["kind"] for e in E)
    print("by provenance kind:", dict(kinds))
    lvl = sum(1 for e in E if e["provenance"].get("method"))
    print(f"method recorded: {lvl} of {len(E)} ({100 * lvl // max(len(E), 1)} %)")
    unk = kinds.get("unknown", 0)
    print(f"provenance unknown: {unk} of {len(E)} ({100 * unk // max(len(E), 1)} %)")
    print("with reference calculations on record:", sum(1 for e in E if e.get("reference_calculations")))
    print("needs a meaningful name:", sum(1 for e in E if e.get("needs_name")))
    if args.list_unknown:
        for e in E:
            if e["provenance"]["kind"] == "unknown":
                print("  unknown:", e["id"], "| comment:", repr(e["provenance"].get("comment", ""))[:60])
    return 0


def test_id(path):
    p = path.replace("\\", "/")
    m = re.match(r"test_cases/cli/([^/]+)/([^/]+)/", p)
    return f"cli_{m.group(1)}_{m.group(2)}" if m else None


def usage(args):
    m = load_manifest()
    corpus = []
    for dp, dn, fn in os.walk(os.path.join(ROOT, "test_cases")):
        if "structures" in dp.split(os.sep) or "GMTKN55-testset" in dp or os.sep + "release" in dp:
            continue
        for f in fn:
            if f.endswith((".txt", ".cmake", ".sh", ".py", ".cpp")) or f == "CMakeLists.txt":
                try:
                    corpus.append((os.path.relpath(os.path.join(dp, f), ROOT),
                                   open(os.path.join(dp, f), encoding="utf-8", errors="replace").read()))
                except OSError:
                    pass
    corpus.append(("CMakeLists.txt", open(os.path.join(ROOT, "CMakeLists.txt"), encoding="utf-8", errors="replace").read()))
    lists = collections.defaultdict(list)
    for dp, dn, fn in os.walk(os.path.join(ROOT, "test_cases")):
        if "structures.txt" in fn and "GMTKN55-testset" not in dp and os.sep + "release" not in dp:
            who = test_id(os.path.relpath(dp, ROOT).replace(os.sep, "/") + "/") or os.path.relpath(dp, ROOT)
            for line in open(os.path.join(dp, "structures.txt"), encoding="utf-8"):
                s = line.split("#")[0].strip()
                if s:
                    lists[s.split()[0]].append(who)
    want = set(args.ids)
    for e in m["structures"]:
        if want and e["id"] not in want:
            continue
        users = set(lists.get(e["id"], []))
        for lp in e.get("legacy_paths", []):
            t = test_id(lp)
            if t:
                users.add(t)
            for name, text in corpus:
                if lp in text or lp.replace("test_cases/", "") in text:
                    users.add(name)
        print(f"{e['id']}: {', '.join(sorted(users)) or '-'}")
    return 0


def add(args):
    src = os.path.abspath(args.file)
    if not re.match(ID_RE, args.id):
        sys.exit(f"id must match {ID_RE.pattern}")
    if args.cls not in CLASSES:
        sys.exit(f"class must be one of {CLASSES}")
    if args.kind == "unknown":
        sys.exit("kind unknown is not allowed for new structures")
    m = load_manifest()
    if any(e["id"] == args.id for e in m["structures"]):
        sys.exit(f"id {args.id} exists; a changed geometry needs a new id (use --supersedes)")
    frames, n, sym, xyz, comments = parse_xyz(src)
    if frames != 1:
        sys.exit("only single-frame xyz can be added with this command; ensembles need the class ensembles and --keep-comments")
    level = "/".join(x for x in [(args.program or "") + (("-" + args.version) if args.version else ""), args.method or ""] if x) or "n/a"
    comment = f"id={args.id} charge={args.charge} spin={args.spin} level={level} source={args.kind}" + (f":{args.reference}" if args.reference else "")
    with open(src, encoding="utf-8") as f:
        lines = f.read().split("\n")
    lines[1] = comment
    dest_dir = os.path.join(LIB, args.cls)
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, args.id + ".xyz")
    with open(dest, "w", encoding="utf-8") as f:
        f.write("\n".join(lines).rstrip("\n") + "\n")
    prov = {"kind": args.kind}
    for k in ("program", "version", "method", "basis", "solvent", "convergence", "reference", "license", "description", "evidence"):
        v = getattr(args, k)
        if v:
            prov[k] = v
    if args.energy is not None:
        prov["energy_eh"] = args.energy
    entry = {"id": args.id, "class": args.cls, "file": f"{args.cls}/{args.id}.xyz", "format": "xyz",
             "formula": hill_formula(sym), "natoms": n, "frames": 1, "charge": args.charge, "spin": args.spin,
             "role": args.role, "size_bytes": os.path.getsize(dest), "sha256": sha256(dest), "provenance": prov}
    if args.derived_from:
        entry["derived_from"] = args.derived_from
    if args.supersedes:
        entry["supersedes"] = args.supersedes
    if args.variant_of:
        entry["variant_of"] = args.variant_of
    if args.notes:
        entry["notes"] = args.notes
    m["structures"].append(entry)
    save_manifest(m)
    print("added", entry["file"])
    return check(argparse.Namespace())


def migrate(args):
    """Replace the local copies of library structures in a test directory by a structures.txt list."""
    tdir = os.path.relpath(os.path.abspath(args.testdir), ROOT).replace(os.sep, "/")
    m = load_manifest()
    by_path = {}
    for e in m["structures"]:
        for lp in e.get("legacy_paths", []):
            by_path[lp] = e
    tracked = subprocess.check_output(["git", "ls-files", tdir], text=True, cwd=ROOT).split("\n")
    moved, lines = [], []
    for p in sorted(x for x in tracked if x):
        e = by_path.get(p)
        if e is None:
            continue
        name = os.path.basename(p)
        lines.append(f"{e['id']}" + ("" if name == e["id"] + os.path.splitext(name)[1] else f" as {name}"))
        moved.append((p, e))
    if not moved:
        print("nothing to migrate in", tdir)
        return 0
    for p, e in moved:
        if sha256(os.path.join(ROOT, p)) != e["sha256"]:
            sys.exit(f"{p} differs from the library file of {e['id']}; migration must be byte-identical")
    lst = os.path.join(ROOT, tdir, "structures.txt")
    existing = []
    if os.path.exists(lst):
        existing = [l.rstrip("\n") for l in open(lst, encoding="utf-8") if l.strip() and not l.startswith("#")]
    out = sorted(set(existing + lines))
    print(f"{tdir}: {len(moved)} file(s) -> structures.txt ({len(out)} lines)")
    if args.dry_run:
        for l in out:
            print("  ", l)
        return 0
    with open(lst, "w", encoding="utf-8") as f:
        f.write("# Structures from test_cases/structures (id, optionally 'as <file name>'); see structures/README.md\n")
        f.write("\n".join(out) + "\n")
    for p, e in moved:
        subprocess.check_call(["git", "rm", "-q", p], cwd=ROOT)
        e["legacy_paths"].remove(p)
        if not e["legacy_paths"]:
            del e["legacy_paths"]
    subprocess.check_call(["git", "add", os.path.join(tdir, "structures.txt")], cwd=ROOT)
    save_manifest(m)
    return 0


def retire(args):
    """git rm legacy copies of library structures (byte-identical only) and drop them from legacy_paths."""
    m = load_manifest()
    by_path = {lp: e for e in m["structures"] for lp in e.get("legacy_paths", [])}
    paths = args.paths or sorted(by_path)
    for p in paths:
        e = by_path.get(p)
        if e is None:
            sys.exit(f"{p} is not a legacy path of any structure")
        if sha256(os.path.join(ROOT, p)) != e["sha256"]:
            sys.exit(f"{p} differs from the library file of {e['id']}; not retiring")
    for p in paths:
        e = by_path[p]
        subprocess.check_call(["git", "rm", "-q", p], cwd=ROOT)
        e["legacy_paths"].remove(p)
        if not e["legacy_paths"]:
            del e["legacy_paths"]
    save_manifest(m)
    print(f"retired {len(paths)} legacy file(s)")
    return 0


def library_path(sid):
    """Absolute path of the library file of a structure id (for scripts)."""
    for e in load_manifest()["structures"]:
        if e["id"] == sid:
            return os.path.join(LIB, e["file"])
    raise KeyError(f"structure {sid} is not in the library")


def path_cmd(args):
    for sid in args.ids:
        print(library_path(sid))
    return 0


def stage_cmd(args):
    """Copy library structures into DEST under the given names (like curcuma_stage_structures in CMake)."""
    os.makedirs(args.dest, exist_ok=True)
    for item in args.items:
        name, _, sid = item.partition("=")
        src = library_path(sid)
        shutil.copyfile(src, os.path.join(args.dest, name + os.path.splitext(src)[1]))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check").set_defaults(fn=check)
    r = sub.add_parser("report"); r.add_argument("--list-unknown", action="store_true"); r.set_defaults(fn=report)
    u = sub.add_parser("usage"); u.add_argument("ids", nargs="*"); u.set_defaults(fn=usage)
    rt = sub.add_parser("retire"); rt.add_argument("paths", nargs="*"); rt.set_defaults(fn=retire)
    pa = sub.add_parser("path"); pa.add_argument("ids", nargs="+"); pa.set_defaults(fn=path_cmd)
    st = sub.add_parser("stage"); st.add_argument("dest"); st.add_argument("items", nargs="+", help="NAME=ID"); st.set_defaults(fn=stage_cmd)
    mg = sub.add_parser("migrate"); mg.add_argument("testdir"); mg.add_argument("--dry-run", action="store_true"); mg.set_defaults(fn=migrate)
    a = sub.add_parser("add"); a.set_defaults(fn=add)
    a.add_argument("file"); a.add_argument("--id", required=True); a.add_argument("--class", dest="cls", required=True)
    a.add_argument("--charge", type=int, required=True); a.add_argument("--spin", type=int, required=True, help="number of unpaired electrons")
    a.add_argument("--role", default="unspecified", choices=ROLES)
    a.add_argument("--kind", required=True, choices=[k for k in KINDS if k != "unknown"])
    for k in ("program", "version", "method", "basis", "solvent", "convergence", "reference", "license", "description", "evidence", "derived-from", "supersedes", "variant-of", "notes"):
        a.add_argument("--" + k)
    a.add_argument("--energy", type=float, help="energy in Eh at the recorded level")
    args = ap.parse_args()
    sys.exit(args.fn(args))


if __name__ == "__main__":
    main()
