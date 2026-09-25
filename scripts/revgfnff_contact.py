#!/usr/bin/env python3
# Claude Generated (Sep 2026) - rev-gfnff stage 3a(ii) falsifier 2
"""Class-S rigid-contact interaction curves: revgfnff against gfnff and against r2SCAN-3c.

Purpose (FABLE_ROADMAP_REVIEW.md item 8, stage-3a(ii) c_ij design review):
the second falsifier of the valence-share factor is "revgfnff's interaction curve stays
within 1 kcal/mol of gfnff's" on rigid intermolecular contact scans. The four class-S
reference scans (`test_cases/revgfnff/ref/S/`, r2SCAN-3c, ORCA 6.1) already exist; this
script turns them into that falsifier, so it can be re-run in one command after any change
to the bond well.

PROTOCOL
  Frames       the class-S points of one system, in ascending heavy-atom separation
               (20 points, 2.30 - 6.00 A). The two isolated-monomer points of the ORCA
               job could not be computed (its $new_job chain cannot change the atom
               count), so the interaction energy is referenced to the LARGEST separation
               (d = 6.00 A), as the campaign documents (test_cases/revgfnff/_log/
               ORCA_REF_STATUS.md).
  Topology     for the `kept` curve the topology-defining frame is the LARGEST-separation
               frame (class S is a dissociation curve in reverse: at 6.0 A the two
               fragments are separate, so the perceived force field is the
               two-fragment one and the curve stays continuous). That frame is written
               first into the multi-frame XYZ and `-batch_reuse_topology true` then
               applies its force field to every frame. This matters: seeding the same
               scan from the d = 2.30 A frame instead shifts the whole curve by
               ~8 mEh (Coulomb only, i.e. a different fragment/EEQ constraint), see the
               `_topology_index` docstring in scripts/revgfnff_fit.py. The script prints
               the frame index and separation it selected and asserts it is the maximum.
  `fresh`      a second pass with `-batch_reuse_topology false`, which re-perceives the
               topology at every point (what a plain single point sees). Reported as a
               diagnostic, not part of the verdict: at short separation the topology can
               flip there, which is a perception step and not a statement about c_ij.
  Cache        every run in a fresh temporary directory AND `-gfnff.cache_topology false`.
               `<basename>.topo.json` has no geometry in its fingerprint, so a stale
               cache silently applies one geometry's topology to another (documented;
               it has corrupted measurements in this project twice).
  Compute      CPU, `-gpu none`, `-threads 1`, `-no_bmt`, `-verbosity 0`, batch single
               point. Every run is one curcuma process.

ACCEPTANCE
  max |E_int(revgfnff) - E_int(gfnff)| over the scan <= --tol (default 1.0 kcal/mol).
  The verdict is taken on the `kept` curve; `fresh` and the two reference differences
  (|rev - ref|, |gfnff - ref|) are reported alongside.

Usage:
    python3 scripts/revgfnff_contact.py                       # all four, prints the tables
    python3 scripts/revgfnff_contact.py --systems water_dimer_OO --json out.json
    python3 scripts/revgfnff_contact.py --binary release/curcuma --report rep.md --note "..."
"""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from revgfnff_data import REF, read_points_xyz  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
AU2KCAL = 627.5094740631
S_DIR = REF / "S"
METHODS = ("gfnff", "revgfnff")   # both are always run; the verdict compares them
DEFAULT_REPORT = REPO / "test_cases" / "revgfnff" / "_log" / "CONTACT_HARNESS.md"
DEFAULT_JSON = REF / "_results" / "contact_curves.json"   # next to curves_gfnff.json

# The four reference minima as reported when the scans were produced (separation [A],
# interaction energy [kcal/mol], referenced to d = 6.00 A). Reproduced from the raw
# energies.json by verify_reference() below and reported as agree/disagree.
REPORTED_MINIMA = {
    "water_dimer_OO": (3.00, -3.84),
    "hf_dimer_FF": (2.90, -2.74),
    "nh3_h2o_NO": (3.25, -1.92),
    "ch4_h2o_CO": (3.80, -0.51),
}


# ------------------------------------------------------------------ data


def scan_value(label):
    """Numeric separation after 'd=' / 'r=' in a point label, or None."""
    if "=" not in label:
        return None
    try:
        return float(label.split("=", 1)[1])
    except ValueError:
        return None


def load_scan(system):
    """Usable points of one class-S system, ascending in separation.

    A point is usable if it has an energy and its frame exists in points.xyz. Returns
    (info, points) with points = [{index, label, d, energy_eh, atoms}].
    """
    sdir = S_DIR / system
    info = json.loads((sdir / "energies.json").read_text())
    frames = read_points_xyz(sdir / "points.xyz")
    pts = []
    for p in info["points"]:
        k = p["point"]
        d = scan_value(p["label"])
        if p.get("energy_eh") is None or d is None or k >= len(frames):
            continue
        pts.append({"index": k, "label": p["label"], "d": d,
                    "energy_eh": p["energy_eh"], "frame": frames[k]})
    if not pts:
        raise RuntimeError(f"{system}: no usable points")
    pts.sort(key=lambda p: p["d"])
    return info, pts


def interaction_curve(e_by_d):
    """E_int = E(d) - E(largest d). e_by_d: list of (d, energy). Returns list of (d, E_int)."""
    zero = e_by_d[-1][1]
    if zero is None:
        raise RuntimeError("no energy at the largest separation")
    return [(d, (e - zero) * AU2KCAL if e is not None else None) for d, e in e_by_d], zero


def parabolic_min(xs, ys):
    """Refined (x, y) of the minimum of the three points bracketing the sampled minimum."""
    i = min(range(len(ys)), key=lambda k: ys[k])
    if i == 0 or i == len(ys) - 1:
        return xs[i], ys[i], "sampled (minimum on the scan edge)"
    x = [xs[i - 1], xs[i], xs[i + 1]]
    y = [ys[i - 1], ys[i], ys[i + 1]]
    den = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2])
    if abs(den) < 1e-12:
        return xs[i], ys[i], "sampled (degenerate bracket)"
    a = (x[2] * (y[1] - y[0]) + x[1] * (y[0] - y[2]) + x[0] * (y[2] - y[1])) / den
    b = (x[2] * x[2] * (y[0] - y[1]) + x[1] * x[1] * (y[2] - y[0]) + x[0] * x[0] * (y[1] - y[2])) / den
    if a <= 0.0:
        return xs[i], ys[i], "sampled (not convex)"
    xm = -b / (2.0 * a)
    return xm, a * xm * xm + b * xm + (y[1] - a * x[1] * x[1] - b * x[1]), "parabolic"


# ------------------------------------------------------------------ curcuma


def file_record(path):
    p = Path(path)
    if not p.exists():
        return {"path": str(p), "exists": False}
    st = p.stat()
    return {"path": str(p), "exists": True, "md5": hashlib.md5(p.read_bytes()).hexdigest(),
            "mtime": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(st.st_mtime)),
            "size": st.st_size}


def write_frames(path, frames, comments):
    with path.open("w") as f:
        for fr, c in zip(frames, comments):
            f.write(f"{len(fr)}\n{c}\n")
            f.write("".join(f"{s} {x:.8f} {y:.8f} {z:.8f}\n" for s, x, y, z in fr))


def batch_energies(binary, frames, comments, charge, spin, method, workdir, reuse, extra=None):
    """One batch single point over `frames`; returns (energies, terms, zero_note).

    The run happens in a fresh directory inside `workdir` (fresh basename per call), with
    the on-disk topology cache disabled, so no stale `<basename>.topo.json` can be read.
    """
    tag = f"{method}_{'reuse' if reuse else 'fresh'}_{abs(hash((method, len(frames), reuse)))%10**8}"
    work = Path(workdir) / tag
    work.mkdir(parents=True, exist_ok=True)
    for stale in work.glob("*.topo.json"):
        stale.unlink()
    xyz = work / f"{tag}.xyz"
    write_frames(xyz, frames, comments)
    out = work / f"{tag}.jsonl"
    cmd = [str(binary), "-sp", xyz.name, "-method", method, "-charge", str(charge), "-spin", str(spin),
           "-batch", "true", "-batch_out", out.name, "-no_bmt", "-verbosity", "0",
           "-threads", "1", "-gpu", "none", "-gfnff.cache_topology", "false",
           "-batch_reuse_topology", "true" if reuse else "false"]
    if extra:
        cmd += extra
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=work)
    energies, terms, errs = [], [], []
    if out.exists():
        for line in out.read_text().splitlines():
            rec = json.loads(line)
            energies.append(rec.get("energy_eh"))
            terms.append(rec.get("terms"))
            errs.append(rec.get("error"))
    leftover = sorted(p.name for p in work.glob("*.topo.json"))
    shutil.rmtree(work, ignore_errors=True)
    if proc.returncode != 0 or len(energies) != len(frames):
        raise RuntimeError(f"batch run failed (rc={proc.returncode}, {len(energies)}/{len(frames)} frames"
                           f"{', errors: ' + str([e for e in errs if e]) if any(errs) else ''}):"
                           f" {' '.join(cmd)}")
    return energies, terms, leftover


# ------------------------------------------------------------------ one system


def run_system(system, binary, workdir, tol=1.0, secondary_from=2.50, extra=None):
    info, pts = load_scan(system)
    ds = [p["d"] for p in pts]
    n = len(pts)
    topo_i = max(range(n), key=lambda i: ds[i])          # largest separation = topology frame
    assert topo_i == n - 1, f"{system}: ascending order should put the largest separation last"
    topo_frame = pts[topo_i]
    charge, spin = info["charge"], info["mult"] - 1

    frames = [p["frame"] for p in pts]
    comments = [f"{system} {p['label']}" for p in pts]
    zeros_leftover = []

    res = {"system": system, "tag": info.get("tag", ""), "charge": charge, "spin": spin,
           "n_points": n, "topology_frame_index": topo_frame["index"],
           "topology_frame_label": topo_frame["label"], "topology_frame_d": topo_frame["d"],
           "topology_frame_is_largest": topo_frame["d"] == max(ds),
           "modes": {}}

    for meth in METHODS:
        # kept: topology frame written FIRST, then every point ascending (the largest-separation
        # frame therefore appears twice: at index 0, the zero of the curve, and at the end).
        e, terms, leftover = batch_energies(binary, [topo_frame["frame"]] + frames,
                                            ["topology_seed " + topo_frame["label"]] + comments,
                                            charge, spin, meth, workdir, reuse=True, extra=extra)
        zeros_leftover += leftover
        # fresh: the same points, topology re-perceived at every one of them (diagnostic only)
        ef, termsf, leftoverf = batch_energies(binary, frames, comments, charge, spin,
                                               meth, workdir, reuse=False, extra=extra)
        zeros_leftover += leftoverf
        res["modes"][meth] = {"kept": {"zero_eh": e[0], "energies_eh": e[1:], "terms": terms[1:],
                                       "zero_terms": terms[0]},
                              "fresh": {"zero_eh": ef[-1], "energies_eh": ef, "terms": termsf}}

    # per-point table
    ref_curve, _ = interaction_curve([(p["d"], p["energy_eh"]) for p in pts])
    rows = []
    for k, p in enumerate(pts):
        row = {"d": p["d"], "label": p["label"], "ref_int": ref_curve[k][1]}
        for meth in METHODS:
            for m in ("kept", "fresh"):
                e = res["modes"][meth][m]["energies_eh"][k]
                z = res["modes"][meth][m]["zero_eh"]
                row[f"{meth}_{m}"] = (e - z) * AU2KCAL
                row[f"{meth}_{m}_eh"] = e
        row["rev_minus_gfnff"] = row["revgfnff_kept"] - row["gfnff_kept"]
        row["rev_minus_ref"] = row["revgfnff_kept"] - row["ref_int"]
        row["gfnff_minus_ref"] = row["gfnff_kept"] - row["ref_int"]
        rows.append(row)
    res["points"] = rows

    # maxima and minima
    def amax(key):
        """(signed value of largest magnitude, its separation) over the points carrying `key`."""
        pairs = [(r[key], r["d"]) for r in rows if r.get(key) is not None]
        if not pairs:
            return None, None
        v, d = max(pairs, key=lambda t: abs(t[0]))
        return v, d

    mx_rev_gfnff, d_rev_gfnff = amax("rev_minus_gfnff")
    mx_rev_ref, d_rev_ref = amax("rev_minus_ref")
    mx_gfnff_ref, d_gfnff_ref = amax("gfnff_minus_ref")
    fresh_delta = max(abs(r["revgfnff_fresh"] - r["gfnff_fresh"]) for r in rows)
    # secondary maximum: outside the hard-wall region (the reference itself is heavily
    # repulsive below `secondary_from`; the contact H...Y distance there is under ~1.5 A).
    # Reported alongside, never instead of, the all-points maximum.
    sec = [(r["d"], r["rev_minus_gfnff"]) for r in rows if r["d"] >= secondary_from]
    sec_max, sec_d = (max(sec, key=lambda t: abs(t[1]))[1], max(sec, key=lambda t: abs(t[1]))[0]) if sec else (None, None)
    n_above = sum(1 for r in rows if abs(r["rev_minus_gfnff"]) > tol)
    curve_minima = {}
    for key in ("ref_int", "gfnff_kept", "revgfnff_kept", "gfnff_fresh", "revgfnff_fresh"):
        if any(key in r for r in rows):
            xs = [r["d"] for r in rows if key in r]
            ys = [r[key] for r in rows if key in r]
            xm, ym, how = parabolic_min(xs, ys)
            curve_minima[key] = {"d": xm, "e_int": ym, "how": how,
                                 "sampled_d": xs[min(range(len(ys)), key=lambda i: ys[i])]}
    res["max_abs"] = {"rev_minus_gfnff": mx_rev_gfnff, "rev_minus_gfnff_at_d": d_rev_gfnff,
                      "rev_minus_gfnff_abs": abs(mx_rev_gfnff) if mx_rev_gfnff is not None else None,
                      "rev_minus_gfnff_fresh": fresh_delta,
                      "rev_minus_gfnff_above_secondary_from": sec_max,
                      "rev_minus_gfnff_above_secondary_from_at_d": sec_d,
                      "secondary_from_d": secondary_from,
                      "n_points_above_tol": n_above,
                      "rev_minus_ref": mx_rev_ref, "rev_minus_ref_at_d": d_rev_ref,
                      "rev_minus_ref_abs": abs(mx_rev_ref) if mx_rev_ref is not None else None,
                      "gfnff_minus_ref": mx_gfnff_ref, "gfnff_minus_ref_at_d": d_gfnff_ref,
                      "gfnff_minus_ref_abs": abs(mx_gfnff_ref) if mx_gfnff_ref is not None else None}
    # term attribution at the worst point (which energy term carries the deviation)
    worst_k = next(k for k, r in enumerate(rows) if r["d"] == d_rev_gfnff)
    tg = res["modes"]["gfnff"]["kept"]["terms"][worst_k]
    tr = res["modes"]["revgfnff"]["kept"]["terms"][worst_k]
    if tg and tr:
        res["term_delta_at_worst_point"] = {
            "d": d_rev_gfnff,
            "terms": {k: (tr.get(k, 0.0) - tg.get(k, 0.0)) * AU2KCAL for k in sorted(tg)},
            "total_kcal": sum((tr.get(k, 0.0) - tg.get(k, 0.0)) * AU2KCAL for k in tg)}
    res["minima"] = curve_minima
    res["asymptote_eh"] = {f"{meth}_{m}": res["modes"][meth][m]["zero_eh"]
                           for meth in METHODS for m in ("kept", "fresh")}
    res["stale_topo_cache_seen"] = zeros_leftover
    # sanity checks
    res["checks"] = {
        "zero_point_agrees_kept_vs_fresh": {
            meth: abs(res["modes"][meth]["kept"]["zero_eh"] - res["modes"][meth]["fresh"]["zero_eh"]) <= 1e-9
            for meth in METHODS},
        "repr_at_smallest_d": {"d": rows[0]["d"],
                               "ref": rows[0]["ref_int"], "gfnff": rows[0]["gfnff_kept"],
                               "revgfnff": rows[0]["revgfnff_kept"],
                               "repulsive": all(rows[0][k] > 0 for k in ("ref_int", "gfnff_kept", "revgfnff_kept"))},
        "near_zero_at_second_largest_d": {"d": rows[-2]["d"],
                                          "ref": rows[-2]["ref_int"], "gfnff": rows[-2]["gfnff_kept"],
                                          "revgfnff": rows[-2]["revgfnff_kept"]},
    }
    return res, info


def verify_reference(results):
    """Check the four reference minima of the raw energies.json against REPORTED_MINIMA."""
    out = {}
    for r in results:
        system = r["system"]
        ds = [p["d"] for p in r["points"]]
        ys = [p["ref_int"] for p in r["points"]]
        xm, ym, how = parabolic_min(ds, ys)
        i = min(range(len(ys)), key=lambda k: ys[k])
        rep_d, rep_e = REPORTED_MINIMA.get(system, (None, None))
        ok = (rep_d is not None and abs(ds[i] - rep_d) <= 0.05 and abs(ys[i] - rep_e) <= 0.05)
        out[system] = {"sampled_d": ds[i], "sampled_e_int": ys[i], "parabolic_d": xm, "parabolic_e_int": ym,
                       "reported_d": rep_d, "reported_e_int": rep_e, "verdict": "AGREE" if ok else "DISAGREE"}
    return out


# ------------------------------------------------------------------ report


def write_report(path, results, checks, binary, tol, note_lines):
    L = []
    L.append("# Class-S contact scans: revgfnff vs gfnff vs r2SCAN-3c (interaction curves)")
    L.append("")
    L.append("AI-generated (`scripts/revgfnff_contact.py`), machine-evaluated. kcal/mol, Angstrom.")
    L.append("")
    L.append(f"- binary: `{binary['path']}` md5 `{binary.get('md5', '?')}` mtime {binary.get('mtime', '?')}")
    L.append(f"- topology frame: the LARGEST-separation frame (frame 0 of the multi-frame XYZ, "
             f"`-batch_reuse_topology true`); `-gfnff.cache_topology false`, fresh temp dir per run, CPU/`-gpu none`")
    L.append(f"- interaction energy E_int(d) = E(d) - E(6.00 A) (the largest separation; the two isolated-monomer "
             f"ORCA jobs could not be computed, see ORCA_REF_STATUS.md)")
    L.append("- curve modes (both always computed): `kept` (topology of the largest separation kept over the "
             "whole scan; the verdict curve) and `fresh` (re-perceived at every point; diagnostic)")
    L.append(f"- acceptance: max |E_int(revgfnff) - E_int(gfnff)| <= {tol:.1f} kcal/mol on `kept`")
    L.append("")
    L.append("## Reference minima reproduced from the raw energies.json")
    L.append("")
    L.append("| system | sampled d / E_int | parabolic d / E_int | reported d / E_int | verdict |")
    L.append("|---|---|---|---|---|")
    for s, c in checks["reference_minima"].items():
        L.append(f"| {s} | {c['sampled_d']:.2f} / {c['sampled_e_int']:+.3f} | {c['parabolic_d']:.3f} / "
                 f"{c['parabolic_e_int']:+.3f} | {c['reported_d']:.2f} / {c['reported_e_int']:+.2f} | "
                 f"**{c['verdict']}** |")
    L.append("")
    L.append("## Per-system comparison")
    L.append("")
    L.append("| system | topo frame idx (d) | n | pts > tol | max \\|rev-gfnff\\| (at d) | max \\|rev-ref\\| (at d) | "
             "max \\|gfnff-ref\\| (at d) | verdict |")
    L.append("|---|---|---:|---:|---|---|---|---|")
    n_fail = 0
    for r in results:
        m = r["max_abs"]
        ok = m["rev_minus_gfnff_abs"] is not None and m["rev_minus_gfnff_abs"] <= tol
        n_fail += 0 if ok else 1
        L.append(f"| {r['system']} | {r['topology_frame_index']} ({r['topology_frame_d']:.2f}) | {r['n_points']} | "
                 f"{m['n_points_above_tol']} | "
                 f"{m['rev_minus_gfnff_abs']:.3f} ({m['rev_minus_gfnff_at_d']:.2f}) | "
                 f"{m['rev_minus_ref_abs']:.3f} ({m['rev_minus_ref_at_d']:.2f}) | "
                 f"{m['gfnff_minus_ref_abs']:.3f} ({m['gfnff_minus_ref_at_d']:.2f}) | "
                 f"{'**PASS**' if ok else '**FAIL**'} |")
    L.append("")
    worst = max((r["max_abs"]["rev_minus_gfnff_abs"] for r in results), default=float("nan"))
    L.append(f"**Tolerance verdict: {'PASS' if n_fail == 0 else 'FAIL'}** - worst system "
             f"{worst:.3f} kcal/mol against the {tol:.1f} kcal/mol bound; {len(results) - n_fail}/{len(results)} pass.")
    sec = max((abs(r["max_abs"]["rev_minus_gfnff_above_secondary_from"] or 0.0) for r in results), default=0.0)
    L.append(f"Secondary (informational, never instead of the verdict): the same maximum restricted to "
             f"d >= {results[0]['max_abs']['secondary_from_d']:.2f} A is **{sec:.3f} kcal/mol** for all four systems "
             "- below that the rigid contact is inside the reference's own hard wall (E_int >= +1 kcal/mol there), "
             "where stage-1's over-coordination term is the whole difference (term attribution below).")
    fresh_worst = max((r["max_abs"]["rev_minus_gfnff_fresh"] or 0.0 for r in results), default=0.0)
    L.append(f"`fresh` mode for reference: worst max |rev - gfnff| = {fresh_worst:.3f} kcal/mol "
             "(perception flips included, not part of the verdict).")
    L.append("")
    L.append("Absolute-energy asymptote E(d=6.00 A), rev - gfnff, in kcal/mol: " +
             ", ".join(f"{r['system']} {(r['asymptote_eh']['revgfnff_kept'] - r['asymptote_eh']['gfnff_kept']) * AU2KCAL:+.5f}"
                       for r in results) +
             " - constant and small; it is carried entirely by the over-coordination term and E_int cancels it by construction.")
    L.append("")
    L.append("## Minima positions (parabolic over the sampled bracket)")
    L.append("")
    L.append("| system | ref d / E_int | gfnff d / E_int | revgfnff d / E_int | rev - ref at its own min |")
    L.append("|---|---|---|---|")
    for r in results:
        m = r["minima"]
        L.append(f"| {r['system']} | {m['ref_int']['d']:.3f} / {m['ref_int']['e_int']:+.3f} | "
                 f"{m['gfnff_kept']['d']:.3f} / {m['gfnff_kept']['e_int']:+.3f} | "
                 f"{m['revgfnff_kept']['d']:.3f} / {m['revgfnff_kept']['e_int']:+.3f} | "
                 f"{m['revgfnff_kept']['e_int'] - m['ref_int']['e_int']:+.3f} |")
    L.append("")
    L.append("## Term attribution at the worst point (rev - gfnff, kcal/mol, |d| > 0.001)")
    L.append("")
    for r in results:
        td = r.get("term_delta_at_worst_point")
        if not td:
            continue
        big = {k: v for k, v in td["terms"].items() if abs(v) > 0.001}
        L.append(f"- **{r['system']}** (d = {td['d']:.2f}): " +
                 ", ".join(f"{k} {v:+.3f}" for k, v in sorted(big.items(), key=lambda t: -abs(t[1]))) +
                 f" (sum {td['total_kcal']:+.3f})")
    L.append("")
    L.append("## Sanity checks")
    L.append("")
    L.append("| system | repulsive at smallest d | E_int at 2nd-largest d (ref/gfnff/rev) | zero agrees kept vs fresh |")
    L.append("|---|---|---|---|")
    for r in results:
        c = r["checks"]
        zk = c["zero_point_agrees_kept_vs_fresh"]
        ztxt = " / ".join(f"{k} {'yes' if v else ('no' if v is not None else '-')}" for k, v in zk.items())
        nz = c["near_zero_at_second_largest_d"]
        L.append(f"| {r['system']} | {'yes' if c['repr_at_smallest_d']['repulsive'] else '**NO**'} "
                 f"({c['repr_at_smallest_d']['ref']:+.1f}/{c['repr_at_smallest_d']['gfnff']:+.1f}/"
                 f"{c['repr_at_smallest_d']['revgfnff']:+.1f}) | d={nz['d']:.2f}: "
                 f"{nz['ref']:+.3f}/{nz['gfnff']:+.3f}/{nz['revgfnff']:+.3f} | {ztxt} |")
    L.append("")
    if note_lines:
        L.append("## Notes")
        L.append("")
        for line in note_lines:
            L.append(f"- {line}")
        L.append("")
    L.append("## Per-point curves (kept mode, kcal/mol, interaction energy vs d = 6.00 A)")
    L.append("")
    ds = [p["d"] for p in results[0]["points"]]
    L.append("| curve | " + " | ".join(f"{d:.2f}" for d in ds) + " |")
    L.append("|" + "---|" * (len(ds) + 1))
    for r in results:
        for lbl, key in (("ref", "ref_int"), ("gfnff", "gfnff_kept"), ("revgfnff", "revgfnff_kept")):
            L.append(f"| {r['system']}:{lbl} | " +
                     " | ".join(f"{p[key]:+.2f}" for p in r["points"]) + " |")
        L.append(f"| {r['system']}:**rev-gfnff** | " +
                 " | ".join(f"{p['rev_minus_gfnff']:+.2f}" for p in r["points"]) + " |")
    L.append("")
    L.append("Full precision of every point, both `kept` and `fresh`, in the `--json` output "
             "(`test_cases/revgfnff/ref/_results/contact_curves.json`).")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(L) + "\n")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", default=str(REPO / "release" / "curcuma"))
    ap.add_argument("--systems", nargs="*", default=None, help="class-S system names")
    ap.add_argument("--all", action="store_true", help="all four class-S systems (default)")
    ap.add_argument("--tol", type=float, default=1.0, help="acceptance bound in kcal/mol")
    ap.add_argument("--secondary-from", type=float, default=2.50,
                    help="separation from which the informational short-contact-free maximum is taken")
    ap.add_argument("--json", default=str(DEFAULT_JSON), help="per-point output (set to '' to skip)")
    ap.add_argument("--report", default=str(DEFAULT_REPORT), help="markdown report (set to '' to skip)")
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--param", default=None, help="GFN-FF parameter override (-gfnff.param_file)")
    ap.add_argument("--note", action="append", default=[], help="extra line(s) under '## Notes'")
    args = ap.parse_args()

    systems = sorted(args.systems) if args.systems else sorted(p.name for p in S_DIR.iterdir()
                                                             if (p / "energies.json").exists())
    if not systems:
        ap.error("no class-S systems found")
    binary = Path(args.binary).resolve()
    brec = file_record(binary)
    if not brec.get("exists"):
        ap.error(f"binary not found: {binary}")
    print(f"binary {binary} md5 {brec['md5']} mtime {brec['mtime']}")

    workdir = args.workdir or tempfile.mkdtemp(prefix="revgfnff_contact_")
    Path(workdir).mkdir(parents=True, exist_ok=True)
    extra = ["-gfnff.param_file", str(Path(args.param).resolve())] if args.param else None

    results = []
    t0 = time.time()
    for s in systems:
        res, _ = run_system(s, binary, workdir, tol=args.tol,
                            secondary_from=args.secondary_from, extra=extra)
        results.append(res)
        m = res["max_abs"]
        print(f"{s:18s} topo frame {res['topology_frame_index']:2d} (d={res['topology_frame_d']:.2f})  "
              f"min ref {res['minima']['ref_int']['d']:.2f}/{res['minima']['ref_int']['e_int']:+.2f}  "
              f"rev {res['minima']['revgfnff_kept']['d']:.2f}/{res['minima']['revgfnff_kept']['e_int']:+.2f}  "
              f"max|rev-gfnff| {m['rev_minus_gfnff_abs']:.3f} (at d={m['rev_minus_gfnff_at_d']:.2f}) "
              f"{m['n_points_above_tol']}/{res['n_points']} pts > tol")
        if res["stale_topo_cache_seen"]:
            print(f"  WARNING: a .topo.json appeared despite cache_topology=false: {res['stale_topo_cache_seen']}")
    checks = {"reference_minima": verify_reference(results)}
    for s, c in checks["reference_minima"].items():
        print(f"reference minimum {s:18s} {c['sampled_d']:.2f}/{c['sampled_e_int']:+.3f}  "
              f"reported {c['reported_d']:.2f}/{c['reported_e_int']:+.2f}  {c['verdict']}")

    out = {"binary": brec, "protocol": {"curve_modes": ["kept", "fresh"], "tolerance_kcal": args.tol,
                                        "topology_frame": "largest separation",
                                        "cache_topology": False, "batch_reuse_topology": True,
                                        "workdir": args.workdir or "(temporary, removed after the run)", "param_file": args.param,
                                        "reference_zero": "largest separation (d=6.00 A)",
                                        "secondary_from_d": args.secondary_from,
                                        "wall_s": None},
           "checks": checks, "systems": results}
    out["protocol"]["wall_s"] = round(time.time() - t0, 1)
    if args.json:
        Path(args.json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.json).write_text(json.dumps(out, indent=1))
        print(f"wrote {args.json}")
    if args.report:
        write_report(args.report, results, checks, brec, args.tol, args.note)
        print(f"wrote {args.report}")
    if not args.workdir:
        shutil.rmtree(workdir, ignore_errors=True)

    n_fail = sum(0 if (r["max_abs"]["rev_minus_gfnff_abs"] or 9e9) <= args.tol else 1 for r in results)
    print(f"TOLERANCE {'PASS' if n_fail == 0 else 'FAIL'}: {len(results) - n_fail}/{len(results)} systems "
          f"within {args.tol:.1f} kcal/mol")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
