#!/usr/bin/env python3
# Claude Generated (Sep 2026) - rev-gfnff WP2/WP1d
"""Loader for the r2SCAN-3c reference data written by scripts/revgfnff_ref.py.

iter_points(classes=None) yields one dict per reference point that has an energy and a
gradient: {class, system, tag, label, charge, mult, atoms [(sym,x,y,z) Angstrom],
energy_eh, gradient_eh_ang (3N list), s2, hirshfeld}. Points still being computed are
simply absent, so a fitter can be run on a partial campaign.

    python scripts/revgfnff_data.py          # summary of what is on disk
"""
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
REF = REPO / "test_cases" / "revgfnff" / "ref"
BOHR = 0.529177210903


def read_points_xyz(path):
    lines = path.read_text().splitlines()
    frames, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        frames.append([(t[0], float(t[1]), float(t[2]), float(t[3])) for t in (ln.split() for ln in lines[i + 2:i + 2 + n])])
        i += n + 2
    return frames


def iter_systems(classes=None):
    for cls_dir in sorted(REF.glob("[A-Z]")):
        if classes and cls_dir.name not in classes:
            continue
        for sdir in sorted(cls_dir.iterdir()):
            e = sdir / "energies.json"
            if e.exists():
                yield cls_dir.name, sdir


def iter_points(classes=None):
    for cls, sdir in iter_systems(classes):
        info = json.loads((sdir / "energies.json").read_text())
        grads = json.loads((sdir / "gradients.json").read_text())["gradients"] if (sdir / "gradients.json").exists() else []
        frames = read_points_xyz(sdir / "points.xyz")
        for p in info["points"]:
            k = p["point"]
            if p["energy_eh"] is None or k >= len(grads) or grads[k] is None or k >= len(frames):
                continue
            yield {"class": cls, "system": info["system"], "tag": info.get("tag", ""), "label": p["label"],
                   "charge": info["charge"], "mult": info["mult"], "atoms": frames[k], "energy_eh": p["energy_eh"],
                   "gradient_eh_ang": [g / BOHR for row in grads[k] for g in row], "s2": p.get("s2"),
                   "hirshfeld": p.get("hirshfeld")}


def summary():
    by_cls = {}
    for cls, sdir in iter_systems():
        info = json.loads((sdir / "energies.json").read_text())
        by_cls.setdefault(cls, []).append((info["system"], info["n_ok"], info["n_points"], info.get("wall_s", 0)))
    total = 0
    for cls in sorted(by_cls):
        rows = by_cls[cls]
        ok = sum(r[1] for r in rows)
        total += ok
        print(f"class {cls}: {len(rows)} systems, {ok}/{sum(r[2] for r in rows)} points, {sum(r[3] for r in rows) / 3600:.1f} h ORCA wall")
        for r in rows:
            flag = "" if r[1] == r[2] else "  INCOMPLETE"
            print(f"    {r[0]:32s} {r[1]:3d}/{r[2]:3d}  {r[3]:7.0f} s{flag}")
    print(f"total usable points: {total}")


if __name__ == "__main__":
    summary()
