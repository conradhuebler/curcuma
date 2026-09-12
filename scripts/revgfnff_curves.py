#!/usr/bin/env python3
# Claude Generated (Sep 2026) - rev-gfnff WP2 diagnostics
"""GFN-FF against the r2SCAN-3c reference curves (classes A, C, E of scripts/revgfnff_ref.py).

For every reference system the GFN-FF energies are computed twice with the batch single
point: (1) with the topology of the reference (equilibrium) geometry kept along the whole
curve ("kept": what a reactive run that has the bond sees), and (2) with a fresh topology
perception at every point ("fresh": what a plain single point sees, i.e. where the bond
list flips). Both are compared with the reference relative to its own minimum.

Per curve the report gives: r_eq of reference and GFN-FF (the minimum of the sampled
points), the well depth D_e = E(largest r) - E(min) of the reference (UKS series) and of
GFN-FF, the RMS error of the relative energies in the bonded region (r <= 1.5 r_eq) and
in the dissociation tail, and the point where the fresh perception drops the bond.

Usage:
    python scripts/revgfnff_curves.py                # all classes, writes the report
    python scripts/revgfnff_curves.py --classes A --only h2 ch3oh --method gfnff --param ov.json
"""
import argparse
import json
import math
import subprocess
import sys
import tempfile
from collections import OrderedDict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from revgfnff_data import REF, iter_systems, read_points_xyz  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
CURCUMA = REPO / "release" / "curcuma"
AU2KCAL = 627.509474
OUT = REPO / "test_cases" / "revgfnff" / "ref" / "_results"


def write_frames(path, frames, comments=None):
    with path.open("w") as f:
        for k, fr in enumerate(frames):
            f.write(f"{len(fr)}\n{(comments or [''] * len(frames))[k]}\n")
            f.write("".join(f"{s} {x:.8f} {y:.8f} {z:.8f}\n" for s, x, y, z in fr))


def batch_energies(frames, charge, spin, method, param_file=None, reuse_topology=False, extra=None):
    """GFN-FF (or any method) energies of the frames via the batch single point, kcal/mol list."""
    work = Path(tempfile.mkdtemp(prefix="revgfnff_curve_"))
    xyz = work / "frames.xyz"
    write_frames(xyz, frames)
    cmd = [str(CURCUMA), "-sp", "frames.xyz", "-method", method, "-charge", str(charge), "-spin", str(spin),
           "-batch", "true", "-batch_out", "frames.jsonl", "-no_bmt", "-verbosity", "0", "-threads", "1"]
    if reuse_topology:
        cmd += ["-batch_reuse_topology", "true"]
    if param_file:
        cmd += ["-gfnff.param_file", str(param_file)]
    if extra:
        cmd += extra
    subprocess.run(cmd, capture_output=True, text=True, cwd=work)
    out = []
    jl = work / "frames.jsonl"
    if jl.exists():
        for line in jl.read_text().splitlines():
            rec = json.loads(line)
            out.append(rec.get("energy_eh"))
    for f in work.glob("*"):
        f.unlink()
    work.rmdir()
    return out


def r_of_label(label):
    if label.startswith("r=") or label.startswith("d="):
        return float(label[2:])
    return None


def rms(v):
    v = [x for x in v if x is not None]
    return math.sqrt(sum(x * x for x in v) / len(v)) if v else float("nan")


def analyse_system(cls, sdirs, method, param_file, ref_geom_cache):
    """sdirs: the directories of one curve (RKS and/or UKS series); the reference at each r is their minimum."""
    if not isinstance(sdirs, list):
        sdirs = [sdirs]
    info = json.loads((sdirs[0] / "energies.json").read_text())
    by_r = OrderedDict()
    for sdir in sdirs:
        inf = json.loads((sdir / "energies.json").read_text())
        frs = read_points_xyz(sdir / "points.xyz")
        for p, fr in zip(inf["points"], frs):
            r = r_of_label(p["label"])
            key = round(r, 4) if r is not None else p["point"]
            e = p["energy_eh"]
            if key not in by_r:
                by_r[key] = [p["point"], r, e, fr]
            elif e is not None and (by_r[key][2] is None or e < by_r[key][2]):
                by_r[key][2] = e
    pts = [(v[0], v[1], v[2]) for v in by_r.values()]
    frames = [v[3] for v in by_r.values()]
    # reference geometry of the molecule (for the kept topology): first token of the system name
    mol = info["system"].split("_")[0]
    ref_xyz = REF / "_geom" / mol / "ref.xyz"
    if not ref_xyz.exists():
        return None
    ref_atoms = read_points_xyz(ref_xyz)[0]
    spin = info["mult"] - 1
    charge = info["charge"]
    # GFN-FF: kept topology (reference geometry prepended so frame 0 defines the topology; class C
    # has an extra atom, there the first point (largest d? no: d=1.0) is used as it is)
    if cls == "A" or cls == "E":
        kept = batch_energies([ref_atoms] + frames, charge, spin, method, param_file, reuse_topology=True)[1:]
    else:
        far = frames[-1]  # d = 3.0 A: attacker not bonded; keeps the host topology + free attacker
        kept = batch_energies([far] + frames, charge, spin, method, param_file, reuse_topology=True)[1:]
    fresh = batch_energies(frames, charge, spin, method, param_file, reuse_topology=False)
    rows = []
    for (k, r, e_ref), fr, ek, ef in zip(pts, frames, kept, fresh):
        rows.append({"r": r, "ref": e_ref, "kept": ek, "fresh": ef})
    rows = [row for row in rows if row["r"] is not None]
    rows.sort(key=lambda row: row["r"])
    ok = [row for row in rows if row["ref"] is not None and row["kept"] is not None]
    if len(ok) < 4:
        return {"system": info["system"], "class": cls, "n": len(ok), "note": "too few points"}
    ref_min = min(ok, key=lambda row: row["ref"])
    kept_min = min(ok, key=lambda row: row["kept"])
    r_eq_ref, r_eq_ff = ref_min["r"], kept_min["r"]
    de_ref = (ok[-1]["ref"] - ref_min["ref"]) * AU2KCAL
    de_ff = (ok[-1]["kept"] - kept_min["kept"]) * AU2KCAL
    # relative energies to each curve's own minimum
    err_bond, err_tail, err_fresh = [], [], []
    flip = None
    for row in ok:
        rel_ref = (row["ref"] - ref_min["ref"]) * AU2KCAL
        rel_ff = (row["kept"] - kept_min["kept"]) * AU2KCAL
        (err_bond if row["r"] <= 1.5 * r_eq_ref else err_tail).append(rel_ff - rel_ref)
        if row["fresh"] is not None:
            err_fresh.append((row["fresh"] - row["kept"]) * AU2KCAL)
            if flip is None and abs(row["fresh"] - row["kept"]) * AU2KCAL > 1.0 and row["r"] > r_eq_ref:
                flip = row["r"]
    return {"system": info["system"], "class": cls, "tag": info.get("tag", ""), "n": len(ok),
            "r_eq_ref": r_eq_ref, "r_eq_ff": r_eq_ff, "de_ref": de_ref, "de_ff": de_ff,
            "rms_bond": rms(err_bond), "rms_tail": rms(err_tail), "max_tail": max(err_tail, key=abs) if err_tail else float("nan"),
            "fresh_flip_r": flip, "rows": [{k: (v * AU2KCAL if k in ("ref", "kept", "fresh") and v is not None else v)
                                            for k, v in row.items()} for row in ok],
            "e_ref_min": ref_min["ref"], "e_ff_min": kept_min["kept"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--classes", nargs="+", default=["A", "C", "E"])
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--method", default="gfnff")
    ap.add_argument("--param", help="GFN-FF parameter override file (-gfnff.param_file)")
    ap.add_argument("--tag", default="", help="suffix for the output files")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    results = []
    systems = [(cls, sdir) for cls, sdir in iter_systems(args.classes)
               if not args.only or any(sdir.name.startswith(o) for o in args.only)]
    # class A: the RKS and UKS series of one bond form one variational curve, E_ref(r) = min(RKS, UKS)
    merged, seen = [], set()
    for cls, sdir in systems:
        base = sdir.name[:-4] if sdir.name.endswith(("_rks", "_uks")) else sdir.name
        if (cls, base) in seen:
            continue
        seen.add((cls, base))
        partners = [d for c, d in systems if c == cls and (d.name == base or d.name[:-4] == base)]
        merged.append((cls, base, partners))
    for cls, base, partners in merged:
        res = analyse_system(cls, partners, args.method, args.param, {})
        if res:
            res["system"] = base
        if res:
            results.append(res)
            if "note" not in res:
                print(f"{cls} {res['system']:28s} r_eq {res['r_eq_ref']:.3f}/{res['r_eq_ff']:.3f}  De {res['de_ref']:7.1f}/{res['de_ff']:7.1f}"
                      f"  rms bond {res['rms_bond']:6.2f} tail {res['rms_tail']:6.2f}  flip {res['fresh_flip_r']}")
    name = f"curves_{args.method}{('_' + args.tag) if args.tag else ''}"
    (OUT / f"{name}.json").write_text(json.dumps(results, indent=1))
    lines = [f"# {args.method} vs r2SCAN-3c reference curves" + (f" (overrides {args.param})" if args.param else ""), "",
             "AI-generated (scripts/revgfnff_curves.py), machine-evaluated. kcal/mol, Angstrom. 'kept' = topology of the",
             "equilibrium structure kept along the curve; 'fresh' = topology re-perceived at every point.",
             "D_e = E(largest r) - E(min) on the sampled grid (reference: UKS series where available).", "",
             "| class | system | n | r_eq ref / FF | D_e ref / FF | RMS bonded (r<=1.5 r_eq) | RMS tail | max tail | fresh flips at r |",
             "|---|---|---:|---|---|---:|---:|---:|---:|"]
    for r in results:
        if "note" in r:
            lines.append(f"| {r['class']} | {r['system']} | {r['n']} | {r['note']} | | | | | |")
            continue
        lines.append(f"| {r['class']} | {r['system']} | {r['n']} | {r['r_eq_ref']:.3f} / {r['r_eq_ff']:.3f} | {r['de_ref']:.1f} / {r['de_ff']:.1f} | "
                     f"{r['rms_bond']:.2f} | {r['rms_tail']:.2f} | {r['max_tail']:+.1f} | {r['fresh_flip_r'] if r['fresh_flip_r'] is not None else '-'} |")
    (OUT / f"{name}.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT / name}.md/.json")


if __name__ == "__main__":
    main()
