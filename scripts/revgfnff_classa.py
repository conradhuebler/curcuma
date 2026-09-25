#!/usr/bin/env python3
# Claude Generated (Sep 2026) - rev-gfnff class-A bond-well harness
"""rev-gfnff class-A bond-well harness: the model's bond well against the r2SCAN-3c reference.

One command for the measurement that four earlier jobs each rebuilt in their own scratchpad
(_log/FABLE_ROADMAP_REVIEW.md section 2.1/2.2, _log/CLASSA_FROZENCN.md, _log/R0_FIX_STATUS.md,
_log/OUTLIER_STATUS.md).  Stage 3a (iii) - the well form - needs it on every iteration.

For every class-A bond type under test_cases/revgfnff/ref/A/ the model is run along the reference
grid and, per bond type, the harness reports (all kcal/mol, Angstrom, and relative to each curve's
own minimum):

    rms      RMS of (model - reference) over the common grid points
    r_eq     radius of the sampled minimum
    D_e      E(largest r) - E(min) on the sampled grid
    r50/r90  interpolated radius where the curve's own rise reaches 50 / 90 % of its own D_e,
             expressed in units of its own r_eq
    k        second derivative at the sampled minimum (3-point Lagrange), kcal/mol/A^2

model and reference side by side plus the SIGNED deviation of each, so "too shallow" and "too
steep" stay distinguishable.

    python scripts/revgfnff_classa.py --all --method revgfnff --json /tmp/classa.json
    python scripts/revgfnff_classa.py --systems ch4_C-H h2o_O-H --mode fresh

Definitions were fixed by reproducing the recorded numbers, not by choosing them: r_eq/D_e/k and
the "own D_e, own r_eq" rise radii reproduce _log/CLASSA_FROZENCN.md table 2 (ch4_C-H 115.2 /
1.584 / 2.224 / 776.4; h2o_O-H 121.5 / 1.541 / 2.082 / 1220.1) and _log/OUTLIER_STATUS.md
section C (the same quantities in Angstrom: 2.4270 = 2.224 * 1.0914) to the printed digit.

TWO PROTOCOL TRAPS, both of which have silently corrupted measurements in this tree

  1. kept vs fresh.  GFN-FF re-perceives the bond graph per geometry by default, and above
     roughly 1.3 x the covalent sum a stretched pair is not perceived as a bond at all
     (getnb).  The energy difference between "the pair is a bond" and "the pair is not a bond" is
     18-117 kcal/mol on these curves (measured: _log/OUTLIER_STATUS.md section F) - one to two
     orders of magnitude above the effects the class-A data are read for.  --mode kept (the
     DEFAULT here) prepends the frame at the reference r_eq to the batch and then passes
     -batch_reuse_topology true AND -gfnff.reuse_topology_check false, so frame 0's bond graph
     is kept for the whole scan: the run still sees the bond it is stretching.  --mode fresh
     re-perceives every frame and is offered only because it answers the OTHER question ("where
     does a plain single point stop seeing the bond"); it is the wrong protocol for a
     bond-stretch scan, and a well-shape number taken with it measures the perception threshold,
     not the well.
     Why the second flag is not optional (commit 76e7f83a, "Fix batch calculator reuse running
     every frame on the first frame's topology", cherry-picked on top of ceb2d160): that commit
     made -batch_reuse_topology true RE-PERCEIVE the bond graph per frame by default, because a
     reused calculator's frames are not necessarily one trajectory (it also serves conformer
     series and multi-molecule batches).  batch reuse therefore no longer provides the
     kept-topology semantics this mode is named after: without the opt-out, --mode kept silently
     measures the fresh protocol while still being recorded as "kept" - the exact confusion
     _log/OUTLIER_STATUS.md section F is about, and the reason the recorded class-A numbers use
     the kept protocol.  -gfnff.reuse_topology_check false restores "trust frame 0's graph", i.e.
     the protocol every recorded class-A number was measured with.  See _log/CLASSA_HARNESS.md
     ("Protocol repair" section) for the before/after verification.
     A caller can still override: if --extra names reuse_topology_check, this harness adds no
     flag of its own and the caller's value decides.
     A third protocol statement sits between the two: -method revgfnff with
     -gfnff.topology_mode react re-detects the bond graph with hysteresis inside the kept run.
     That is what _log/CLASSA_FROZENCN.md mode 1 and _log/OUTLIER_STATUS.md "react" are, and it
     is reached with --extra "-gfnff.topology_mode react"; a plain revgfnff run differs from it
     (1.9 kcal/mol on ch4_C-H).

  2. <basename>.topo.json.  The on-disk topology cache's fingerprint is the element list plus the
     bond graph - with NO geometry - so a stale cache applies one geometry's perceived topology
     (and its Phase-1 EEQ charges) to another.  It has corrupted measurements here twice
     (_log/R0_FIX_STATUS.md section 2, Known Issue #11/#21(c)).  This harness therefore gives
     every (bond type, mode, method) its own fresh directory AND passes
     -gfnff.cache_topology false, and it asserts that no *.topo.json was left behind.

Reference side (unchanged from the four earlier jobs, so the numbers stay comparable):
pointwise min(RKS, UKS) at every radius - the same convention as
scripts/revgfnff_curves.py::analyse_system and scripts/revgfnff_data.py.  Do not "improve" it.
On top of it the exclusions of test_cases/revgfnff/ref/QUALITY.md are applied (see the two
tables below).

AI-generated, machine-evaluated.  No src/ changes, no build.
"""
import argparse
import datetime
import hashlib
import json
import math
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from revgfnff_data import read_points_xyz  # noqa: E402  (shared frame reader)

REPO = Path(__file__).resolve().parents[1]
REF = REPO / "test_cases" / "revgfnff" / "ref"
AU2KCAL = 627.509474

# ---------------------------------------------------------------------------------------------
# QUALITY.md exclusions.  Audit these against test_cases/revgfnff/ref/QUALITY.md before editing.
# ---------------------------------------------------------------------------------------------

# QUALITY.md section 3: five class-A UKS series where two runs of the same input land on
# different broken-symmetry solutions (listed there with dE and dS^2 per radius).  12 radii.
# The whole radius is dropped from the merged curve: min(RKS, UKS) has no trustworthy value
# there - the UKS branch is unstable and the RKS branch is the substitution QUALITY.md section 5
# warns against.  Consequence: class H cannot stand in for class A on these five series.
QUALITY_UNSTABLE = {
    "c2h2_CTC_uks": [2.4008],
    "f2_F-F_uks": [1.5400, 1.8200, 1.9600, 2.8000],
    "hcn_CTN_uks": [3.4520],
    "n2_NTN_uks": [3.2822, 1.6411],
    "n2h2_NDN_uks": [1.7335, 1.8573, 1.9811, 2.4764],
}

# QUALITY.md section 4 / section 5 point 2: the class-A UKS trees of these two series were
# overwritten by the SlowConv retry and have no far point (16/20).  Filling the gap from the RKS
# branch would make the reference well too deep by about +53.1 (of2) / +51.5 (clf) kcal/mol, so
# every radius without a UKS number is dropped rather than substituted.
# cl2_Cl-Cl added Sep 22, 2026 (CL2_WELLFIT_P1_STATUS.md): its class-A UKS tree originally
# converged on only 4/20 points (all far); a --uks-inside-out --slowconv recompute raised that
# to 15/20 (r = 1.52-4.06 A, missing only the far tail r >= 4.57 A). Without this entry the
# fallback silently used RKS on the whole dissociating side, which cannot dissociate (rises to
# +34.5 kcal/mol at 4.57 A relative to the minimum) and fitted the Cl-Cl mg well ~15 kcal/mol too
# deep (D_e 70.1 vs the corrected ~55).
QUALITY_REQUIRE_UKS = {"of2_O-F", "clf_F-Cl", "cl2_Cl-Cl"}

# Fixed r/r_eq ratios at which the per-point residual is reported (the earlier logs index their
# tables this way).  The reference grid is r_eq * 1.0668^k, so these land on grid points.
REPORT_RATIOS = (1.0, 1.2, 1.3, 1.4, 1.6, 2.0, 2.5)


def r_of_label(label):
    """Radius from a reference label ('r=1.0691', 'd=3.0'); None for anything else."""
    if label.startswith("r=") or label.startswith("d="):
        return float(label[2:])
    return None


def bond_types():
    """The class-A bond types on disk, derived from the series directory names."""
    bases = set()
    for d in (REF / "A").iterdir():
        if (d / "energies.json").exists():
            bases.add(d.name.rsplit("_", 1)[0])
    return sorted(bases)


def series_of(base):
    return sorted(d for d in (REF / "A").iterdir() if d.name.startswith(base + "_"))


def load_reference(base):
    """Merged reference curve of one bond type.

    Pointwise min(RKS, UKS) keyed by the r label, exactly as scripts/revgfnff_curves.py
    ::analyse_system does it.  QUALITY.md exclusions are applied on top; per-series coverage and
    exclusion counts are returned so that "fewer points" is never silent.
    """
    entries = {}          # rounded r -> dict(r, energy, frame, series, has_uks)
    coverage = []         # (series, n_points, n_ok, n_missing, n_excluded_quality)
    excluded_radii = set()
    n_quality = 0
    for d in series_of(base):
        info = json.loads((d / "energies.json").read_text())
        frames = read_points_xyz(d / "points.xyz")
        excluded = {round(x, 4) for x in QUALITY_UNSTABLE.get(d.name, [])}
        n_ok = n_missing = n_excl = 0
        for k, p in enumerate(info["points"]):
            r = r_of_label(p["label"])
            if r is None or k >= len(frames):
                continue
            key = round(r, 4)
            if key in excluded:
                # the radius goes, not just this branch: QUALITY.md section 3 says a consumer of
                # min(RKS, UKS) must drop the radius.  Dropping only the UKS number would silently
                # substitute RKS, i.e. exactly the substitution QUALITY.md section 5 warns about.
                excluded_radii.add(key)
                n_excl += 1
                continue
            if p["energy_eh"] is None:
                n_missing += 1
                continue
            n_ok += 1
            cur = entries.get(key)
            if cur is None:
                entries[key] = {"r": r, "energy": p["energy_eh"], "frame": frames[k],
                                "series": d.name, "has_uks": d.name.endswith("_uks")}
            else:
                cur["has_uks"] = cur["has_uks"] or d.name.endswith("_uks")
                if p["energy_eh"] < cur["energy"]:
                    cur.update({"energy": p["energy_eh"], "frame": frames[k], "series": d.name})
        n_quality += n_excl
        coverage.append((d.name, len(info["points"]), n_ok, n_missing, n_excl))
    rows = sorted((e for k, e in entries.items() if k not in excluded_radii), key=lambda e: e["r"])
    if base in QUALITY_REQUIRE_UKS:
        # a radius with no UKS number at all is dropped, not filled in from RKS
        before = len(rows)
        rows = [e for e in rows if e["has_uks"]]
        n_quality += before - len(rows)
    info0 = json.loads((series_of(base)[0] / "energies.json").read_text())
    return {"base": base, "rows": rows, "coverage": coverage,
            "charge": info0["charge"], "mult": info0["mult"], "tag": info0.get("tag", ""),
            "n_quality_excluded": n_quality}


def second_derivative(rs, es, i):
    """3-point Lagrange second derivative at grid index i, kcal/mol/A^2 (nan at the end points)."""
    if i <= 0 or i >= len(rs) - 1:
        return float("nan")
    h1, h2 = rs[i] - rs[i - 1], rs[i + 1] - rs[i]
    d2 = 2.0 * (es[i - 1] * h2 - es[i] * (h1 + h2) + es[i + 1] * h1) / (h1 * h2 * (h1 + h2))
    return d2 * AU2KCAL


def rise_radius(rs, rel, r_eq, de, frac):
    """Radius (linear interpolation on the ascending branch) where the curve's rise is frac * D_e."""
    target = frac * de
    for i in range(len(rs) - 1):
        if rs[i] < r_eq - 1e-9:
            continue
        if rel[i] <= target <= rel[i + 1]:
            return rs[i] + (rs[i + 1] - rs[i]) * (target - rel[i]) / (rel[i + 1] - rel[i])
    return float("nan")   # the curve never reaches that level on the sampled grid


def rise_radius_ref_de(rs, rel, r_eq, de_ref):
    """(r50, r90) in Angstrom against the REFERENCE D_e as the level.

    The second of the two conventions in the earlier logs: _log/OUTLIER_STATUS.md section C
    measures the rise against the reference well depth and prints Angstrom, while
    _log/CLASSA_FROZENCN.md (and _log/FABLE_ROADMAP_REVIEW.md section 2.2, whose Morse r50/r90
    come from the curve's own D_e and k_e) measure it against the curve's OWN D_e and print
    r/r_eq.  Both are reported so either can be reproduced; the table uses the CLASSA one.
    """
    return tuple(rise_radius(rs, rel, r_eq, de_ref, f) for f in (0.5, 0.9))


def curve_metrics(rs, eh):
    """r_eq, D_e, r50, r90 (in r/r_eq), k of one curve sampled on the shared grid."""
    if len(rs) < 4:
        return None
    i = min(range(len(rs)), key=lambda j: eh[j])
    r_eq = rs[i]
    de = (eh[-1] - eh[i]) * AU2KCAL
    rel = [(e - eh[i]) * AU2KCAL for e in eh]
    return {"r_eq": r_eq, "d_e": de,
            "r50": rise_radius(rs, rel, r_eq, de, 0.5) / r_eq,
            "r90": rise_radius(rs, rel, r_eq, de, 0.9) / r_eq,
            "k": second_derivative(rs, eh, i), "rel": rel, "i_min": i}


def write_frames(path, frames, comments):
    with path.open("w") as f:
        for k, fr in enumerate(frames):
            f.write(f"{len(fr)}\n{comments[k]}\n")
            f.write("".join(f"{s} {x:.8f} {y:.8f} {z:.8f}\n" for s, x, y, z in fr))


def run_model(binary, base, ref, method, mode, extra, workdir, verbosity):
    """Model energies (Eh) on [r_eq frame] + ascending reference grid, one per frame."""
    rows = ref["rows"]
    imin = min(range(len(rows)), key=lambda j: rows[j]["energy"])
    frames = [rows[imin]["frame"]] + [e["frame"] for e in rows]
    comments = [f"prepended r_eq frame ({base}, r={rows[imin]['r']})"]
    comments += [f"reference grid point {k} (r={e['r']})" for k, e in enumerate(rows)]
    if workdir:
        d = Path(workdir) / base / f"{method}_{mode}"
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
        cleanup = False
    else:
        d = Path(tempfile.mkdtemp(prefix="revgfnff_classa_"))
        cleanup = True
    try:
        xyz = d / "scan.xyz"
        write_frames(xyz, frames, comments)
        out = d / "out.jsonl"
        cmd = [str(binary), "-sp", "scan.xyz", "-method", method,
               "-batch", "true", "-batch_out", "out.jsonl",
               "-batch_reuse_topology", "true" if mode == "kept" else "false",
               "-gfnff.cache_topology", "false",
               "-charge", str(ref["charge"]), "-spin", str(ref["mult"] - 1),
               "-no_bmt", "-threads", "1", "-verbosity", str(verbosity)]
        # Claude Generated (Sep 14, 2026): kept topology is no longer what -batch_reuse_topology
        # true does by itself.  Commit 76e7f83a made calculator reuse re-validate the carried-over
        # bond graph against every frame (PARAM reuse_topology_check) and rebuild when the two
        # differ, so a batch that is not one homogeneous trajectory - which a bond-stretch scan
        # never is - re-perceives each frame.  "kept" would then silently be "fresh".  The
        # opt-out is the documented way back to frame-0 semantics (main.cpp scans the raw argv
        # for it, so it must not be reordered away or folded into the controller).
        if mode == "kept" and not any("reuse_topology_check" in chunk for chunk in extra):
            cmd += ["-gfnff.reuse_topology_check", "false"]
        for chunk in extra:
            cmd += shlex.split(chunk)
        proc = subprocess.run(cmd, capture_output=True, text=True, cwd=d)
        stale = sorted(p.name for p in d.glob("*.topo.json"))
        if stale:
            raise RuntimeError(f"topology cache written into {d}: {stale} - the measurement is unsafe")
        energies = []
        if out.exists():
            for line in out.read_text().splitlines():
                try:
                    energies.append(json.loads(line).get("energy_eh"))
                except json.JSONDecodeError:
                    energies.append(None)
        if proc.returncode != 0 and not energies:
            sys.stderr.write(f"[{base}] curcuma exit {proc.returncode}: "
                             f"{proc.stderr.strip()[-300:]}\n")
        n = len(frames)
        energies = (energies + [None] * n)[:n]
        return energies[1:], energies[0]     # grid, prepended frame
    finally:
        if cleanup:
            shutil.rmtree(d, ignore_errors=True)


def measure(binary, base, method, mode, extra, workdir, verbosity):
    ref = load_reference(base)
    e_grid, e_pre = run_model(binary, base, ref, method, mode, extra, workdir, verbosity)
    rows = ref["rows"]
    rs = [e["r"] for e in rows]
    e_ref = [e["energy"] for e in rows]
    keep = [k for k in range(len(rows)) if e_grid[k] is not None and e_ref[k] is not None]
    rec = {"bond": base, "tag": ref["tag"], "n_grid": len(rows), "n_ok": len(keep),
           "n_quality_excluded": ref["n_quality_excluded"],
           "coverage": [{"series": s, "n_points": n, "n_ok": ok, "n_missing": miss,
                         "n_quality_excluded": q} for s, n, ok, miss, q in ref["coverage"]],
           "prepended_r": rows[min(range(len(rows)), key=lambda j: rows[j]["energy"])]["r"]}
    if len(keep) < 4:
        rec["note"] = "too few usable points"
        return rec
    rs = [rs[k] for k in keep]
    e_ref = [e_ref[k] for k in keep]
    e_mod = [e_grid[k] for k in keep]
    if e_pre is None:
        rec["note"] = "prepended r_eq frame failed"
        return rec
    m_ref = curve_metrics(rs, e_ref)
    m_mod = curve_metrics(rs, e_mod)
    rms = math.sqrt(sum((a - b) ** 2 for a, b in zip(m_mod["rel"], m_ref["rel"])) / len(rs))
    rec.update({"rms": rms,
                "r_eq": {"ref": m_ref["r_eq"], "model": m_mod["r_eq"]},
                "d_e": {"ref": m_ref["d_e"], "model": m_mod["d_e"]},
                "r50": {"ref": m_ref["r50"], "model": m_mod["r50"]},
                "r90": {"ref": m_ref["r90"], "model": m_mod["r90"]},
                "k": {"ref": m_ref["k"], "model": m_mod["k"]},
                "rise_ref_de": {
                    "ref": rise_radius_ref_de(rs, m_ref["rel"], m_ref["r_eq"], m_ref["d_e"]),
                    "model": rise_radius_ref_de(rs, m_mod["rel"], m_mod["r_eq"], m_ref["d_e"]),
                    "unit": "Angstrom, level = 50/90 % of the REFERENCE D_e "
                            "(_log/OUTLIER_STATUS.md section C convention)"},
                "residual_at": {}})
    for name in ("r_eq", "d_e", "r50", "r90", "k"):
        rec[name]["dev"] = rec[name]["model"] - rec[name]["ref"]
    # signed residual at the fixed r/r_eq ratios, on the reference r_eq (as the earlier logs do)
    for ratio in REPORT_RATIOS:
        tgt = m_ref["r_eq"] * ratio
        k = min(range(len(rs)), key=lambda j: abs(rs[j] - tgt))
        if abs(rs[k] - tgt) > 0.02 * m_ref["r_eq"]:
            rec["residual_at"][f"{ratio:g}"] = None
            continue
        rec["residual_at"][f"{ratio:g}"] = m_mod["rel"][k] - m_ref["rel"][k]
    rec["residual_at"]["last"] = m_mod["rel"][-1] - m_ref["rel"][-1]   # largest sampled r
    return rec


def fmt(x, nd=2):
    if x is None:
        return "-"
    if isinstance(x, float) and math.isnan(x):
        return "nan"
    return f"{x:.{nd}f}"


def table(records, method, mode, extra):
    hdr = (f"# class-A bond-well harness: {method}"
           f"{' + ' + ' '.join(extra) if extra else ''}, mode {mode}")
    lines = [hdr, "",
             "AI-generated (scripts/revgfnff_classa.py), machine-evaluated. kcal/mol, Angstrom.",
             "Reference: pointwise min(RKS, UKS) per radius, QUALITY.md exclusions applied.",
             "Every curve is relative to its own minimum; r50/r90 are in units of that curve's own",
             "r_eq (interpolated on the ascending branch). k = 3-point Lagrange at the sampled minimum.",
             "dev = model - reference (signed).", "",
             "| bond | tag | n_ok | excl | rms | r_eq ref/model | D_e ref/model | dev D_e |"
             " r50 ref/model | dev r50 | r90 ref/model | dev r90 | k ref/model | dev k | dev@1.4 | dev@1.6 |",
             "|---|---|---:|---:|---:|---|---|---:|---|---:|---|---:|---|---:|---:|---:|"]
    for r in records:
        if "note" in r:
            lines.append(f"| {r['bond']} | | {r['n_ok']}/{r['n_grid']} | {r['n_quality_excluded']} |"
                         f" | | | | | | | | | | | {r['note']} |")
            continue
        lines.append(
            f"| {r['bond']} | {r['tag']} | {r['n_ok']}/{r['n_grid']} | {r['n_quality_excluded']} |"
            f" {r['rms']:.1f} |"
            f" {r['r_eq']['ref']:.4f}/{r['r_eq']['model']:.4f} |"
            f" {r['d_e']['ref']:.1f}/{r['d_e']['model']:.1f} | {r['d_e']['dev']:+.1f} |"
            f" {fmt(r['r50']['ref'],3)}/{fmt(r['r50']['model'],3)} | {fmt(r['r50']['dev'],3)} |"
            f" {fmt(r['r90']['ref'],3)}/{fmt(r['r90']['model'],3)} | {fmt(r['r90']['dev'],3)} |"
            f" {fmt(r['k']['ref'],0)}/{fmt(r['k']['model'],0)} | {fmt(r['k']['dev'],0)} |"
            f" {fmt(r['residual_at'].get('1.4'),2)} | {fmt(r['residual_at'].get('1.6'),2)} |")
    return "\n".join(lines)


def aggregate(records):
    ok = [r for r in records if "note" not in r]
    lines = ["", "## Aggregate (median / max of the signed deviations)"]
    for name, unit in (("d_e", "kcal/mol"), ("r50", "r/r_eq"), ("r90", "r/r_eq"), ("k", "kcal/mol/A^2")):
        vals = [r[name]["dev"] for r in ok if not math.isnan(r[name]["dev"])]
        if not vals:
            continue
        vals_sorted = sorted(vals)
        med = vals_sorted[len(vals) // 2]
        worst = max(vals, key=abs)
        rb = max(ok, key=lambda r: abs(r[name]["dev"]))
        lines.append(f"- dev {name} [{unit}]: median {med:+.3f}, max |{worst:+.3f}| ({rb['bond']}),"
                     f" nan {len(ok) - len(vals)}/{len(ok)}")
    rms = sorted(r["rms"] for r in ok)
    if rms:
        lines.append(f"- rms of the curve: median {rms[len(rms)//2]:.2f}, max {rms[-1]:.2f} "
                     f"({max(ok, key=lambda r: r['rms'])['bond']}) kcal/mol")
    lines.append(f"- {len(ok)}/{len(records)} bond types measured")
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", default=str(REPO / "release" / "curcuma"),
                    help="curcuma binary (default release/curcuma)")
    ap.add_argument("--method", default="gfnff", help="model method (default gfnff)")
    ap.add_argument("--systems", nargs="+", help="class-A bond types, e.g. ch4_C-H h2o_O-H")
    ap.add_argument("--all", action="store_true", help="every class-A bond type on disk")
    ap.add_argument("--mode", choices=["kept", "fresh"], default="kept",
                    help="kept = bond graph of the r_eq frame kept along the scan (default); "
                         "fresh = re-perceived at every point")
    ap.add_argument("--workdir", help="keep the batch files here (per bond/mode); default temp")
    ap.add_argument("--json", help="write the full per-bond records to this file")
    ap.add_argument("--extra", action="append", default=[],
                    help="extra curcuma arguments, shell-quoted, repeatable "
                         "(e.g. --extra \"-gfnff.topology_mode react\")")
    ap.add_argument("--verbosity", type=int, default=0)
    args = ap.parse_args()

    if not args.systems and not args.all:
        ap.error("give --systems NAME... or --all")
    binary = Path(args.binary).resolve()
    if not binary.exists():
        ap.error(f"binary not found: {binary}")
    want = bond_types() if args.all else args.systems
    unknown = [b for b in want if not series_of(b)]
    if unknown:
        ap.error(f"not a class-A bond type: {' '.join(unknown)}")

    st = binary.stat()
    digest = hashlib.md5(binary.read_bytes()).hexdigest()
    mtime = datetime.datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds")
    print(f"# binary {binary}  md5 {digest}  mtime {mtime}")
    print(f"# method {args.method}  mode {args.mode}  extra {' '.join(args.extra) or '-'}"
          f"  ({len(want)} bond types)")

    records = []
    for base in want:
        rec = measure(binary, base, args.method, args.mode, args.extra, args.workdir, args.verbosity)
        records.append(rec)
        if "note" in rec:
            print(f"  {base:16s} {rec['n_ok']}/{rec['n_grid']}  {rec['note']}")
        else:
            print(f"  {base:16s} n_ok {rec['n_ok']:2d}/{rec['n_grid']}  excl {rec['n_quality_excluded']}"
                  f"  rms {rec['rms']:6.2f}  D_e dev {rec['d_e']['dev']:+7.2f}"
                  f"  r90 dev {fmt(rec['r90']['dev'],3):>7s}  k dev {fmt(rec['k']['dev'],0):>7s}")
        for c in rec["coverage"]:
            if c["n_quality_excluded"] or c["n_missing"]:
                print(f"      {c['series']:24s} n_ok {c['n_ok']:2d}/{c['n_points']}"
                      f"  not converged {c['n_missing']}  QUALITY-excluded {c['n_quality_excluded']}")

    excl_total = sum(r["n_quality_excluded"] for r in records)
    print(f"# QUALITY.md exclusions applied: {excl_total} points over {len(records)} bond types")
    text = table(records, args.method, args.mode, args.extra)
    text += "\n" + "\n".join(aggregate(records)) + "\n"
    if args.json:
        payload = {"binary": str(binary), "md5": digest, "mtime": mtime,
                   "method": args.method, "mode": args.mode, "extra": args.extra,
                   "protocol": "kept" if args.mode == "kept" else "fresh",
                   "reference_branch": "pointwise min(RKS, UKS)",
                   "quality_exclusions": {"unstable_radii": QUALITY_UNSTABLE,
                                          "require_uks": sorted(QUALITY_REQUIRE_UKS),
                                          "n_points_excluded": excl_total},
                   "records": records}
        Path(args.json).write_text(json.dumps(payload, indent=1))
        print(f"# wrote {args.json}")
    else:
        print(text)


if __name__ == "__main__":
    main()
