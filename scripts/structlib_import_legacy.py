#!/usr/bin/env python3
"""One-time import of the structure files that were scattered over test_cases/ into test_cases/structures/.

Claude Generated (Oct 2026). Phase 1 of the structure library: files are copied byte for byte (tests are not
changed), content duplicates are stored once, and every entry lists its `legacy_paths`. Provenance is only DERIVED
from what the repository records (comment lines, xtb .out files, reference_data/*.json); everything else stays
`unknown`. Nothing is guessed. Re-running regenerates the library from the tracked legacy files.

Usage: python3 scripts/structlib_import_legacy.py [--dry-run]
"""
import collections
import glob
import json
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import structlib as sl  # noqa: E402
import pyranoside_stereo as ps  # noqa: E402

ROOT = sl.ROOT
DRY = "--dry-run" in sys.argv

SKIP_RE = re.compile(
    r"^\.|\.(centered|reordered|accepted|rejected|initial|thresh|reuse)\.xyz$|\.reorder\.\d+\.xyz$|\.trj\.xyz$"
    r"|^input\.opt\.xyz$|^optimized\.xyz$|^xtbopt\.xyz$|^aligned_|^reorder_|\.opt\.trj\.")
GENERIC = {"input", "ref", "target", "conformers", "struc", "mol", "molecule", "test", "start", "a", "b", "c"}
RADII = {"H": .31, "B": .84, "C": .76, "N": .71, "O": .66, "F": .57, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02,
         "Br": 1.2, "I": 1.39, "Li": 1.28, "Na": 1.66, "Mg": 1.41, "Al": 1.21}


def tracked_structures():
    out = subprocess.check_output(["git", "ls-files", "test_cases"], text=True, cwd=ROOT).split("\n")
    res, skipped = [], []
    for p in out:
        if not p.lower().endswith((".xyz", ".vtf")) or p.startswith("test_cases/structures/"):
            continue
        (skipped if SKIP_RE.search(os.path.basename(p)) else res).append(p)
    return res, skipped


def fragments(sym, xyz):
    n = len(sym)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i in range(n):
        for j in range(i + 1, n):
            ri, rj = RADII.get(sym[i], 1.4), RADII.get(sym[j], 1.4)
            if sl.math.dist(xyz[i], xyz[j]) < 1.25 * (ri + rj) + 0.1:
                parent[find(i)] = find(j)
    return len({find(i) for i in range(n)})


def classify(path, frames, sym, xyz):
    if path.endswith(".vtf"):
        return "cg"
    if sym and all(s.isdigit() for s in sym):
        return "cg"
    if frames > 1:
        return "ensembles"
    if set(sym) & sl.TRANSITION:
        return "metals"
    n = len(sym)
    if "polymer" in os.path.basename(path).lower() or n > 500:
        return "bulk"
    if n == 1:
        return "atoms"
    if fragments(sym, xyz) > 1:
        return "clusters"
    return "small" if n <= 12 else ("medium" if n <= 60 else "large")


def provenance_from_comment(comment):
    """Derive what the comment line states. Returns (provenance dict, charge, spin)."""
    c = comment.strip()
    prov = {"kind": "unknown"}
    charge = spin = None
    if c:
        prov["comment"] = c[:200]
    m = re.search(r"energy:\s*(-?\d+\.\d+)\s+gnorm:\s*(\d+\.\d+)\s+xtb:\s*([0-9.]+)", c)
    if m:
        prov.update(kind="optimized", program="xtb", version=m.group(3), energy_eh=float(m.group(1)), gnorm=float(m.group(2)),
                    level_note="xtbopt.xyz comment (energy, gnorm, xtb version); Hamiltonian and convergence level not recorded")
        return prov, charge, spin
    m = re.search(r"\*\* Energy =\s*(-?\d+\.\d+) Eh \*\* Charge = (-?\d+) \*\* Spin = (-?\d+) \*\* Curcuma [0-9.]+ \(([^)]*)\)", c)
    if m:
        prov.update(kind="program_output", program="curcuma", version=m.group(4), energy_eh=float(m.group(1)),
                    level_note="curcuma output comment; the method used is not recorded")
        return prov, int(m.group(2)), int(m.group(3))
    m = re.match(r"^-?\d+\.\d+\s*$", c)
    if m:
        prov.update(level_note="a bare number, no unit or program recorded")
        return prov, charge, spin
    if re.search(r"r\([A-Za-z]+\)\s*=|ang\([A-Za-z]+\)\s*=", c):
        prov.update(kind="constructed", description=c[:200],
                    level_note="idealised geometry parameters stated in the comment; no calculation level")
        return prov, charge, spin
    m = re.search(r"\b(gfn2|gfn1|gfn-ff|gfnff|uff|pm3|am1)\b.*\b(relax|minimum|optimi[sz]ed|opt)\b|\b(relax|minimum|optimi[sz]ed)\b.*\b(gfn2|gfn1|gfn-ff|gfnff)\b", c, re.I)
    if m:
        meth = re.search(r"gfn2|gfn1|gfn-ff|gfnff|uff|pm3|am1", c, re.I).group(0).lower().replace("gfn-ff", "gfnff")
        prov.update(kind="optimized", method=meth, level_note="method named in the comment line; program and convergence not recorded")
    return prov, charge, spin


def out_file_info(path):
    """xtb .out next to the structure (same directory and stem): a reference calculation, not an optimisation."""
    base = os.path.splitext(path)[0] + ".out"
    full = os.path.join(ROOT, base)
    if not os.path.exists(full):
        return None
    txt = open(full, encoding="utf-8", errors="replace").read()
    v = re.search(r"xtb version ([0-9.]+)", txt)
    call = re.search(r"program call\s*:\s*\S*xtb\s+(.*)", txt)
    if not (v and call):
        return None
    flags = " ".join(call.group(1).split()[1:])
    return {"kind": "xtb-output", "program": "xtb", "version": v.group(1), "flags": flags, "evidence": base,
            "note": "single-point reference at this geometry (no --opt in the call)" if "--opt" not in flags else "includes --opt"}


def json_index():
    idx = collections.defaultdict(list)
    files = glob.glob(os.path.join(ROOT, "test_cases/reference_data/**/*.json"), recursive=True) + \
        glob.glob(os.path.join(ROOT, "test_cases/sqm_reference/reference_data/*.json"))
    for f in sorted(files):
        try:
            d = json.load(open(f, encoding="utf-8"))
        except Exception:
            continue
        if not isinstance(d, dict):
            continue
        mol = d.get("molecule")
        g = d.get("geometry_file") or (mol.get("geometry_file") if isinstance(mol, dict) else None)
        if not g:
            continue
        rel = os.path.relpath(f, ROOT)
        if "xtb_version" in d:
            info = {"kind": "xtb-sqm", "program": "xtb", "version": d.get("xtb_version"), "method": d.get("method"), "file": rel}
        elif isinstance(mol, dict):
            md = d.get("metadata") or {}
            info = {"kind": "gfnff-terms", "generator": md.get("generator"), "date": md.get("date"), "file": rel}
        else:
            continue
        for k in ("charge", "spin_unpaired", "solvent"):
            if k in d:
                info[k] = d[k]
        idx[os.path.normpath(g)].append(info)
    return idx


# Derived by reading the README next to the structures; applied only to the files named here.
README_EVIDENCE = {
    "test_cases/optimisation/makrocyclic/AnGrad.xyz": {
        "kind": "optimized", "program": "curcuma", "method": "uff",
        "version": "commit 4e19c68f8c252babdcd7c87e95287c8f34177e9a (LBFGSpp commit 7e3848617795ddd0e25f4b772e679adfee583229)",
        "evidence": "test_cases/optimisation/makrocyclic/README.md",
        "level_note": "optimised with analytical gradients; UFF is inferred from the README text (universal force field), convergence not recorded"},
    "test_cases/optimisation/makrocyclic/NumGrad.xyz": {
        "kind": "optimized", "program": "curcuma", "method": "uff",
        "version": "commit 4e19c68f8c252babdcd7c87e95287c8f34177e9a (LBFGSpp commit 7e3848617795ddd0e25f4b772e679adfee583229)",
        "evidence": "test_cases/optimisation/makrocyclic/README.md",
        "level_note": "optimised with numerical gradients; UFF is inferred from the README text (universal force field), convergence not recorded"},
}


# Meaningful ids for the legacy entries. Key: id generated by this script, value: (new id, description or None).
# Each name rests on evidence in the repository: composition/fragment analysis (formula), the comment line, the file or
# directory name, or the docs. Where identity is inferred, the description says how. Names ending in .alt<k> mark
# variants of a substance whose difference to the main variant is not recorded anywhere (needs_name stays set).
RENAME = {
    "polymer_2x": ("peo201x2-water1500", "two chains of composition C402H806O202 (= H-(OCH2CH2)201-OH, PEO-like) and 1500 H2O; from fragment analysis, chain identity inferred from the composition"),
    "polymer_2x_gfn2_opt": ("peo201x2-water1500.gfn2-relaxed", "same system as peo201x2-water1500; comment: GFN2 relaxed (400 LBFGS steps from the GFN-FF minimum)"),
    "polymer_2x_gfnff_opt": ("peo201x2-water1500.gfnff-min", "same system as peo201x2-water1500; comment: GFN-FF minimum"),
    "unnamed.c402h806o202.9ae88f": ("peo201-chain", "one chain C402H806O202 (= H-(OCH2CH2)201-OH, PEO-like); legacy name polymer; chain identity inferred from the composition"),
    "mixture2": ("urea400-water1000", "400 CH4N2O and 1000 H2O (fragment analysis)"),
    "complex": ("macrocycle-bgal", "host C76H108N12O8 plus guest C7H14O6; the guest is methyl beta-D-galactopyranoside by configuration analysis (see analysis); the operator remembered bGlc (2026-10-01), which the coordinates do not support"),
    "angrad": ("macrocycle-bgal.uff-analytic-grad", "optimised with analytical gradients, see optimisation/makrocyclic/README.md"),
    "numgrad": ("macrocycle-bgal.uff-numeric-grad", "optimised with numerical gradients, see optimisation/makrocyclic/README.md"),
    "makrocyclic.input": ("macrocycle-bgal.conf2", "comment line: input_2, a conformer from a curcuma run"),
    "aaa-bgal.a": ("aaa-bglc.conf76", "comment line: input_76, a conformer from a curcuma run; host C36H48N6 plus guest C7H14O6; the guest is methyl beta-D-glucopyranoside by configuration analysis (see analysis), although the directory is named AAA-bGal"),
    "aaa-bgal.b": ("aaa-bglc.conf96", "comment line: input_96, a conformer from a curcuma run; host C36H48N6 plus guest C7H14O6; the guest is methyl beta-D-glucopyranoside by configuration analysis (see analysis), although the directory is named AAA-bGal"),
    "jl22-bgal": ("jl22-bgal.conf6", "comment line: input_6, a conformer from a curcuma run; host C58H78N10 plus guest C7H14O6; the guest is methyl beta-D-galactopyranoside by configuration analysis (see analysis)"),
    "gfn-2": ("jl22-bgal.xtb-gfn2-opt", "file name GFN-2.xyz and xtb 6.6.0 optimisation comment; method taken from the file name"),
    "gfn-ff": ("jl22-bgal.xtb-gfnff-opt", "file name GFN-FF.xyz and xtb 6.6.0 optimisation comment; method taken from the file name"),
    "conf": ("aaa-host-frames353", "353 frames of C36H48N6, the formula of the host fragment of aaa-bgal; identity inferred from the formula (the operator is not sure that the host is AAA)"),
    "unnamed.c36h48n6.614365": ("aaa-host.rmsd-b", "C36H48N6 as used by the cli/rmsd tests (second of the two variants)"),
    "unnamed.c36h48n6.f841c5": ("aaa-host.rmsd-a", "C36H48N6 as used by the cli/rmsd tests (first of the two variants)"),
    "unnamed.c40h59n9o6.a1c30b": ("c33h45n9-agal-frames44", "44 frames of a complex of host C33H45N9 and guest C7H14O6; the guest is methyl alpha-D-galactopyranoside in all 44 frames (see analysis); curcuma conformer output"),
    "helicen": ("helicene-frames17", "17 frames of C26H16 (formula of hexahelicene); what the frames are is not recorded"),
    "triose": ("trisaccharide-c18h32o16", "formula C18H32O16 equals three hexose units minus two H2O; legacy name triose"),
    "triose.input": ("trisaccharide-c18h32o16.opt-start", "input of the optimisation test (directory optimisation/triose); the geometry level is not recorded"),
    "acetic_acid_dimer": ("acetic-acid-dimer", None),
    "unnamed.c4h8o4.42eaf7": ("acetic-acid-dimer.cyclic-hb", "comment line: cyclic hydrogen-bonded acetic acid dimer"),
    "unnamed.c4h10.327f4a": ("butane.anti", "comment line: anti conformation for torsion testing"),
    "unnamed.c2h6.caf8f7": ("ethane", None),
    "unnamed.ch4o.f00834": ("methanol.alt1", None),
    "unnamed.h2o.c08b18": ("water.fast-cli", "comment line: water molecule for fast CLI tests"),
    "unnamed.h3n.7c9e2b": ("ammonia.alt1", None),
    "unnamed.h4.66abcc": ("h-atoms4", "four H atoms (cli/simplemd react recombination test)"),
    "unnamed.cg.fc737e": ("cg-spheres2", "comment line: CG system with 2 spheres (element 226 beads)"),
    "unnamed.x.406d22": ("cg-vtf.analysis-parallel", "VTF input of cli/analysis/01_parallel_equivalence"),
    "unnamed.x.b9d192": ("cg-vtf.single-point", "VTF input of cli/cg/01_single_point"),
    "water": ("water", None),
    "h2o": ("water.ideal-c2v", "comment line: water C2v, r(OH)=0.9572 A"),
    "h2o.mol-larger": ("water.alt1", None),
    "ch4": ("methane", None),
    "ch4.sqmref": ("methane.ideal-td", "comment line: methane Td, r(CH)=1.091 A"),
    "c6h6": ("benzene", None),
    "c6h6.sqmref": ("benzene.ideal-d6h", "comment line: benzene D6h, r(CC)=1.397 A"),
    "benzene": ("benzene.validation", "from the validation/ directory"),
    "c6h5cooh": ("benzoic-acid", None),
    "ch3oh": ("methanol", None),
    "ch3och3": ("dimethyl-ether", None),
    "hcn": ("hydrogen-cyanide", None),
    "hcn.cli-gpu_gradient": ("hydrogen-cyanide.gpu-qresponse", "same coordinates as hydrogen-cyanide; only the comment line differs (cli/gpu_gradient test)"),
    "nh3": ("ammonia.ideal-c3v", "comment line: ammonia C3v, r(NH)=1.012 A, ang(HNH)=106.7 deg"),
    "hcl": ("hydrogen-chloride", None),
    "hcl.mol-larger": ("hydrogen-chloride.alt1", None),
    "hcl.sqmref": ("hydrogen-chloride.d-shell-validation", "comment line: HCl d-shell validation (X-I1)"),
    "h2": ("dihydrogen.r0.741", "comment line: H2 equilibrium bond length 0.741 A"),
    "hh": ("dihydrogen.r0.47", "bond length 0.47 A measured from the coordinates"),
    "lih": ("lithium-hydride.r1.595", "comment line: LiH equilibrium bond length 1.595 A"),
    "he2": ("helium-dimer.r3.0", "comment line: He2 non-bonded test, 3.0 A separation"),
    "h2s": ("hydrogen-sulfide.d-shell-validation", "comment line: H2S d-shell validation (X-I1)"),
    "ph3": ("phosphine.d-shell-validation", "comment line: PH3 d-shell validation (X-I1)"),
    "sih4": ("silane.d-shell-validation", "comment line: SiH4 d-shell validation (X-I1)"),
    "o3": ("ozone", None),
    "oh": ("hydroxyl", None),
    "nacl": ("na2cl2-isolated-ions", "comment line: 4 isolated Na+/Cl- ions"),
    "ethene": ("ethene", None),
    "fi_pyr": ("fi-pyridine", "formula C5H5FIN; GMTKN55 HAL59 FI_pyr (a halogen-bond complex, see docs/KNOWN_ISSUES_ARCHIVE.md #26); identity from the name and the formula"),
    "cl2m_ea25": ("cl2-anion", "GMTKN55 G21EA/EA_25, the dichlorine radical anion (docs/KNOWN_ISSUES_ARCHIVE.md #17); charge -1 per that entry, not in the file"),
    "h3op_h2o2": ("hydronium-water2", "GMTKN55 WATER27 H3O+(H2O)2 (docs/KNOWN_ISSUES_ARCHIVE.md #31); charge +1 per that entry, not in the file"),
    "fragment_2xx": ("polymerbuild-fragment-2xx", "comment line: test fragment with 2 Xx connection atoms; deliberately invalid element Xx"),
    "caffeine": ("caffeine", None),
    "caffeine.opt": ("caffeine.opt", None),
    "water8_cluster": ("water8-cluster", None),
}


def tag_for(path):
    p = path.split("/")
    if p[1] == "cli":
        return "cli-" + p[2]
    if p[1] == "sqm_reference":
        return "sqmref"
    if p[1] == "molecules":
        return "mol-" + (p[2] if len(p) > 3 else "x")
    return p[1].lower()


def main():
    files, skipped = tracked_structures()
    groups = collections.defaultdict(list)
    for p in files:
        groups[sl.sha256(os.path.join(ROOT, p))].append(p)
    jidx = json_index()
    # base ids
    base = {}
    for h, paths in groups.items():
        stems = collections.Counter(os.path.splitext(os.path.basename(p))[0].lower().replace(" ", "_") for p in paths)
        stem = stems.most_common(1)[0][0]
        base[h] = stem
    entries = []
    by_base = collections.defaultdict(list)
    for h, paths in groups.items():
        by_base[base[h]].append(h)
    plan = {}
    for b, hs in by_base.items():
        hs.sort(key=lambda h: (-len(groups[h]), sorted(groups[h])[0]))
        for rank, h in enumerate(hs):
            paths = sorted(groups[h])
            ext = os.path.splitext(paths[0])[1]
            generic = b in GENERIC or len(b) <= 1
            dirs = {os.path.dirname(p) for p in paths}
            if generic:
                if len(dirs) == 1 and "/cli/" not in paths[0]:
                    cid = f"{os.path.basename(next(iter(dirs))).lower()}.{b}"
                    needs = False
                else:
                    cid, needs = None, True
            else:
                cid = b if rank == 0 else f"{b}.{tag_for(paths[0])}"
                needs = False
            plan[h] = (cid, needs, paths, ext)
    used = set()
    unparsable = []
    for h, (cid, needs, paths, ext) in sorted(plan.items(), key=lambda kv: kv[1][2][0]):
        full = os.path.join(ROOT, paths[0])
        if ext == ".xyz":
            try:
                frames, n, sym, xyz, comments = sl.parse_xyz(full)
            except ValueError as ex:
                unparsable.append(f"{paths[0]}: {ex}")
                continue
        else:
            frames, n, sym, xyz, comments = 1, 0, [], [], [""]
        cls = classify(paths[0], frames, sym, xyz)
        formula = sl.hill_formula(sym) if sym else None
        if cid is None:
            cid = f"unnamed.{('cg' if (formula or '').startswith('#') else (formula or 'x')).lower()}.{h[:6]}"
        if cid in used:
            cid = f"{cid}.{h[:6]}"
        cid = re.sub(r"[^a-z0-9_.+-]", "-", cid.lower())
        used.add(cid)
        prov, charge, spin = provenance_from_comment(comments[0])
        bad = [s for s in set(sym) if s not in sl.ELEMSET and not s.isdigit()]
        role = "invalid-input" if bad else ("ensemble" if frames > 1 else "unspecified")
        e = {"id": cid, "class": cls, "file": f"{cls}/{cid}{ext}", "format": ext.lstrip("."), "formula": formula,
             "natoms": n, "frames": frames, "charge": charge, "spin": spin, "role": role,
             "size_bytes": os.path.getsize(full), "sha256": h, "provenance": prov, "legacy": True,
             "legacy_paths": paths}
        stems = sorted({os.path.splitext(os.path.basename(p))[0] for p in paths})
        if len(stems) > 1:
            e["aliases"] = stems
        if needs:
            e["needs_name"] = True
        for p in paths:
            if p in README_EVIDENCE:
                old = e["provenance"]
                e["provenance"] = {**README_EVIDENCE[p], **({"comment": old["comment"]} if "comment" in old else {}),
                                   **({"energy_eh": old["energy_eh"]} if "energy_eh" in old else {})}
        if ext == ".xyz" and n >= 30:
            frames_res = ps.analyse_file(full)
            guests = [r for fr in frames_res for r in fr if "series" in r]
            if guests:
                cnt = collections.Counter((r["series"], r["sugar"], r["anomer"], r["axial_eq"], r["faces"]) for r in guests)
                e["analysis"] = {
                    "method": "scripts/pyranoside_stereo.py (ring faces relative to CH2OH, Haworth orientation); validated on Open Babel built "
                              "reference glycosides and cross-checked with Open Babel's canonical isomeric SMILES of the extracted guest",
                    "guest_formula": "C7H14O6", "frames": len(frames_res),
                    "guest_configuration": [{"series": k[0], "sugar": k[1], "anomer": k[2], "axial_equatorial_C1_to_C5": k[3],
                                             "faces_C1_to_C5": k[4], "frames": v} for k, v in cnt.most_common()]}
        refs = []
        for p in paths:
            info = out_file_info(p)
            if info and info not in refs:
                refs.append(info)
            for info in jidx.get(os.path.normpath(p), []):
                if info not in refs:
                    refs.append(info)
        if refs:
            e["reference_calculations"] = refs
        entries.append(e)
    by_old = {e["id"]: e for e in entries}
    for old, (new, desc) in RENAME.items():
        e = by_old.get(old)
        if e is None:
            continue
        aliases = set(e.get("aliases", [])) | {old}
        e["id"], e["file"] = new, f"{e['class']}/{new}.{e['format']}"
        e["aliases"] = sorted(a for a in aliases if a != new)
        if desc:
            e["description"] = desc
        e.pop("needs_name", None)
        if new.endswith(".alt1"):
            e["needs_name"] = True
    ids = [e["id"] for e in entries]
    dup = [i for i, n in collections.Counter(ids).items() if n > 1]
    if dup:
        sys.exit(f"duplicate ids after renaming: {dup}")
    for e in entries:
        if e["id"].startswith("unnamed."):
            print("still unnamed:", e["id"])
    sets = [
        {"id": "gmtkn55", "kind": "database", "fetch": "python scripts/fetch_testset.py fetch gmtkn55",
         "directory": "test_cases/GMTKN55-testset", "tracked": False, "doc": "docs/GMTKN55_VALIDATION.md",
         "note": "Grimme-group benchmark set; ignored by git, retrieved on demand"},
        {"id": "mor41", "kind": "database", "fetch": "python scripts/fetch_testset.py fetch mor41",
         "directory": "test_cases/MOR41-testset", "tracked": False, "doc": "docs/MOR41_VALIDATION.md",
         "note": "organometallic reaction set; retrieved on demand"},
        {"id": "s30l", "kind": "database", "fetch": "manual (paywalled Supporting Information), see docs/TESTSET_RETRIEVAL.md",
         "directory": "test_cases/s30l_test_set", "tracked": False, "doc": "docs/S30L_GFNNF_VALIDATION.md"},
        {"id": "s30lci", "kind": "database", "fetch": "supplied manually",
         "directory": "test_cases/s30lci_test_set", "tracked": "README and reference file only",
         "doc": "test_cases/s30lci_test_set/README", "note": "counterion variant of S30L"},
    ]
    manifest = {"version": 1,
                "note": "Generated for legacy entries by scripts/structlib_import_legacy.py; maintained with scripts/structlib.py. Rules: test_cases/structures/README.md",
                "structures": entries, "sets": sets,
                "not_imported": {"program_output": sorted(skipped), "unparsable": unparsable}}
    print(f"files considered: {len(files)}, skipped as program output: {len(skipped)}, unique contents: {len(entries)}")
    if DRY:
        print("dry run, nothing written")
        return
    if os.path.isdir(sl.LIB):
        for cls in sl.CLASSES + ["polymers"]:
            shutil.rmtree(os.path.join(sl.LIB, cls), ignore_errors=True)
    os.makedirs(sl.LIB, exist_ok=True)
    for e in entries:
        dest = os.path.join(sl.LIB, e["file"])
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copyfile(os.path.join(ROOT, e["legacy_paths"][0]), dest)
    sl.save_manifest(manifest)


if __name__ == "__main__":
    main()
