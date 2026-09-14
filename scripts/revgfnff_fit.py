#!/usr/bin/env python3
# Claude Generated (Sep 2026) - rev-gfnff WP1d
"""Parameter fitter for rev-gfnff, over the r2SCAN-3c reference data (scripts/revgfnff_data.py).

numpy only (no scipy, not installed here). Evaluates a parameter vector by writing a sparse
GFN-FF override JSON (see gfnff_param_tables.h) and running one `release/curcuma -sp ... -batch
true -batch_reuse_topology true -method revgfnff -gfnff.param_file ...` call PER REFERENCE
SYSTEM, in parallel (each curcuma process pinned to -threads 1). Every system is written once
as a multi-frame XYZ whose first frame is the topology-defining geometry; only the override
JSON changes between evaluations.

Residuals are RELATIVE energies within a system (dE = E - E_topology), because the reference
is absolute r2SCAN-3c and not comparable to force-field absolute energies. Gradients are
compared absolutely, in kcal/mol/Angstrom. See docs/REV_GFNFF_ROADMAP.md WP1d.

Usage:
    python scripts/revgfnff_fit.py --write-example-config config.json
    python scripts/revgfnff_fit.py --config config.json --evaluate-only
    python scripts/revgfnff_fit.py --config config.json --dry-run
    python scripts/revgfnff_fit.py --config config.json --out fit_out --max-iter 8 --jobs 8
"""
import argparse
import json
import math
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import revgfnff_data as rd
import gmtkn55_reactions as gr  # Claude Generated (Sep 2026) - WP3 barrier dataset

REPO = Path(__file__).resolve().parents[1]
CURCUMA = REPO / "release" / "curcuma"
HARTREE_TO_KCAL = 627.5094740631

# ------------------------------------------------------------------ reference-data assembly


@dataclass
class Point:
    label: str
    atoms: list
    energy_eh: float
    gradient_eh_ang: list  # flat 3N list, or None


@dataclass
class System:
    cls: str
    name: str
    charge: int
    mult: int
    points: list  # Point, ordered [topology_frame, ...rest]
    topology_mode: str = "static"  # "static" (no -gfnff.topology_mode flag) or "react"


def _scan_value(label):
    """Numeric value after '=' in a point label ('r=1.23' -> 1.23), or None if unparseable."""
    if "=" not in label:
        return None
    try:
        return float(label.split("=", 1)[1])
    except ValueError:
        return None


def _topology_index(cls, points):
    """Index of the topology-defining frame within one system's point list.

    A/E/L (dissociation-type curves): smallest scan value (bonded geometry).
    C (hyper-coordination) and S (rigid intermolecular contact): largest scan value
    (fragments separated) -- a class-S contact scan is a dissociation curve in reverse.
    D (MD snapshots): the first point, unconditionally.
    Falls back to index 0 if the labels don't carry a parseable scan value.
    """
    if cls == "D":
        return 0
    values = [_scan_value(p["label"]) for p in points]
    if all(v is not None for v in values):
        if cls in ("C", "S"):
            return max(range(len(points)), key=lambda i: values[i])
        return min(range(len(points)), key=lambda i: values[i])
    return 0


def effective_topology_mode(cls, requested):
    """Class D is always static (MD snapshots have no reaction coordinate to order by);
    every other class follows the CLI --topology choice."""
    return "static" if cls == "D" else requested


def _order_points(cls, mode, pts, idx):
    """Order the points of one system into the multi-frame XYZ, topology frame first.

    react mode needs a monotonic reaction-coordinate trajectory so the react topology
    machinery (bond formation/breaking, WP-react) sees a continuous scan instead of jumps:
    class A/E/L go shortest -> longest scan value, class C longest -> shortest (fragments
    approaching), class S likewise (a contact scan run in the association direction, so that
    frame 0 is the separated geometry the _topology_index helper names as the topology frame).
    static mode (and any class whose labels don't carry a scan value) keeps
    the topology frame first and the rest in discovery order, as order is immaterial there.
    """
    if mode == "react" and cls != "D":
        values = [_scan_value(p["label"]) for p in pts]
        if all(v is not None for v in values):
            order = sorted(range(len(pts)), key=lambda i: values[i], reverse=(cls in ("C", "S")))
            return [pts[i] for i in order]
    return [pts[idx]] + [p for i, p in enumerate(pts) if i != idx]


def _merge_class_a(raw_a):
    """Merge the RKS/UKS series of each class-A bond into one curve.

    raw_a: dict system_name ("h2_H-H_rks", ...) -> list of point dicts from iter_points.
    Point-wise, at matching label, the lower-energy (ground-state) series wins; a label
    present in only one series is kept as-is. Returns dict base_name -> merged point list.
    """
    families = {}
    for sysname, pts in raw_a.items():
        base = sysname
        for suf in ("_rks", "_uks"):
            if base.endswith(suf):
                base = base[: -len(suf)]
                break
        families.setdefault(base, {})[sysname] = pts
    merged = {}
    for base, variants in families.items():
        by_label = {}
        for pts in variants.values():
            for p in pts:
                lbl = p["label"]
                if lbl not in by_label or p["energy_eh"] < by_label[lbl]["energy_eh"]:
                    by_label[lbl] = p
        merged[base] = list(by_label.values())
    return merged


def _build_system(cls, name, pts, topology_mode):
    idx = _topology_index(cls, pts)
    mode = effective_topology_mode(cls, topology_mode)
    ordered = _order_points(cls, mode, pts, idx)
    points = [Point(p["label"], p["atoms"], p["energy_eh"], p["gradient_eh_ang"]) for p in ordered]
    return System(cls=cls, name=name, charge=ordered[0]["charge"], mult=ordered[0]["mult"],
                  points=points, topology_mode=mode)


def load_systems(classes, topology_mode="react"):
    """Load every reference system of the requested classes, class-A RKS/UKS merged.

    topology_mode: "static" or "react" (see effective_topology_mode -- class D is always static).
    """
    raw = {}
    for p in rd.iter_points(classes):
        raw.setdefault((p["class"], p["system"]), []).append(p)

    systems = []
    raw_a = {sysname: pts for (cls, sysname), pts in raw.items() if cls == "A"}
    for base, pts in sorted(_merge_class_a(raw_a).items()):
        systems.append(_build_system("A", base, pts, topology_mode))
    for (cls, sysname), pts in sorted(raw.items()):
        if cls == "A":
            continue
        systems.append(_build_system(cls, sysname, pts, topology_mode))
    return systems


# ------------------------------------------------------------------ GMTKN55 barrier dataset (Claude Generated, Sep 2026, WP3)
#
# docs/REV_GFNFF_ROADMAP.md WP3 "Barrier acceptance MEASURED and NOT MET": stage-1's
# over-coordination term was fitted only against class-C rigid approach curves, and the
# fitted N/O penalties wrecked WCPT18/PX13 (proton-transfer TSs, where a real TS atom is
# transiently over-coordinated). This dataset lets the fitter see actual barrier heights,
# not just the artificial hyper-coordination curves.
#
# A barrier "point" is a reaction energy (TS - reactant, or a general tmer2++ stoichiometry),
# built from single-point energies of the GMTKN55 structures involved. Unlike classes A/C/D/E/L
# (one system = one molecule at several geometries, topology shared across frames), each
# GMTKN55 structure is an independent molecule and must get its OWN topology
# (`-batch_reuse_topology false`); structures are batched per (charge, spin) bucket since one
# curcuma process applies a single -charge/-spin to every frame it is given (see main.cpp's
# batch-mode comment).


@dataclass
class BarrierReaction:
    subset: str
    label: str
    keys: list    # ["{structure_dir}/{name}", ...] aligned with coeffs, structure_dir(subset) from gmtkn55_reactions
    coeffs: list
    ref: float    # kcal/mol, published reference barrier/reaction energy


@dataclass
class BarrierGroup:
    charge: int
    mult: int
    keys: list    # ordered unique structure keys sharing this (charge, mult)

    @property
    def name(self):
        return f"barrier_c{self.charge}_m{self.mult}"


def read_single_xyz(path):
    """One-frame XYZ -> [(sym, x, y, z), ...] (Angstrom, as written by the GMTKN55 testset)."""
    lines = path.read_text().splitlines()
    n = int(lines[0].split()[0])
    return [(t[0], float(t[1]), float(t[2]), float(t[3])) for t in (ln.split() for ln in lines[2:2 + n])]


def load_barrier_data(subsets):
    """GMTKN55 barrier reactions + charge/spin batching groups for the given subsets.

    Reuses scripts/gmtkn55_reactions.py for stoichiometry/reference parsing (.res + upstream
    CSV cross-check) and structure_meta() for charge/UHF; aborts on the same inconsistencies
    gmtkn55_reactions.py itself would abort on (mismatched .res/CSV reference, missing rows).
    Returns (reactions: list[BarrierReaction], groups: list[BarrierGroup]).
    """
    rxs, problems = gr.load_reactions(subsets)
    if problems:
        sys.exit("barrier reaction list inconsistent:\n  " + "\n  ".join(problems))
    keys_meta = {}  # "{structure_dir}/{name}" -> (charge, mult)
    reactions = []
    for rx in rxs:
        keys = [f"{gr.structure_dir(rx.subset)}/{s}" for s in rx.species]
        for k, s in zip(keys, rx.species):
            if k not in keys_meta:
                charge, uhf = gr.structure_meta(rx.subset, s)
                keys_meta[k] = (charge, uhf + 1)
        reactions.append(BarrierReaction(rx.subset, rx.label, keys, list(rx.coeffs), rx.ref))
    buckets = {}
    for k, cm in keys_meta.items():
        buckets.setdefault(cm, []).append(k)
    groups = [BarrierGroup(charge=cm[0], mult=cm[1], keys=sorted(ks)) for cm, ks in sorted(buckets.items())]
    return reactions, groups


def write_barrier_group_xyz(group, path):
    lines = []
    for k in group.keys:
        subset_dir, name = k.split("/", 1)
        atoms = read_single_xyz(gr.TESTSET / subset_dir / name / "struc.xyz")
        lines.append(str(len(atoms)))
        lines.append(k)
        for sym, x, y, z in atoms:
            lines.append(f"{sym} {x:.10f} {y:.10f} {z:.10f}")
    path.write_text("\n".join(lines) + "\n")


def run_barrier_group(curcuma, workdir, group, override_path):
    """One curcuma batch single point over a (charge, mult) bucket's structures.

    One topology per structure (-batch_reuse_topology false): unlike run_system()'s classes
    A/C/D/E/L (one molecule, several geometries), every frame here is an unrelated molecule.

    -gfnff.cache_topology false is REQUIRED here, not cosmetic: the on-disk .topo.json cache
    is keyed by a fingerprint of atom count/Z-list/bond graph (gfnff_method.cpp's
    computeTopologyFingerprint()), not geometry or charge. All frames of one batch process
    share one cache-file basename (group.name), and a reactant-complex/TS/product triple
    along a reaction path routinely PERCEIVES THE SAME bond graph at different geometries --
    e.g. GMTKN55 BH76 fch3fcomp and fch3fts. Without this flag, the second frame gets a
    false cache hit and silently reuses the first frame's Phase-1 EEQ charges, which is
    exactly wrong for the barrier dataset (found by comparing a batched fch3fts energy,
    -1.4586 Eh, against a fresh single-structure run, -1.32475 Eh -- the latter matches the
    per-structure reference cache in test_cases/revgfnff/fit_work/cache/terms_revgfnff_default.json).
    """
    xyz_path = workdir / f"{group.name}.xyz"
    out_path = workdir / f"{group.name}.jsonl"
    cmd = [
        str(curcuma), "-sp", str(xyz_path), "-batch", "true", "-batch_out", str(out_path),
        "-batch_reuse_topology", "false", "-method", "revgfnff",
        "-gfnff.param_file", str(override_path), "-gfnff.cache_topology", "false",
        "-charge", str(group.charge), "-spin", str(group.mult - 1),
        "-gradient", "false", "-threads", "1", "-verbosity", "0", "-no_bmt",
    ]
    try:
        subprocess.run(cmd, cwd=str(workdir), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                        timeout=600, check=False)
    except subprocess.TimeoutExpired:
        return []
    if not out_path.exists():
        return []
    frames = []
    for line in out_path.read_text().splitlines():
        line = line.strip()
        if line:
            try:
                frames.append(json.loads(line))
            except json.JSONDecodeError:
                frames.append({"error": "unparseable JSONL line"})
    return frames


def barrier_pairs_from_energies(reactions, energy_by_key):
    """{subset: [(model_kcal, ref_kcal), ...]} from per-structure energies (Eh); reactions
    with any missing structure energy are skipped (not counted as a silent zero)."""
    out = {}
    for rx in reactions:
        if not all(k in energy_by_key for k in rx.keys):
            continue
        model = sum(c * energy_by_key[k] for k, c in zip(rx.keys, rx.coeffs)) * HARTREE_TO_KCAL
        out.setdefault(rx.subset, []).append((model, rx.ref))
    return out


def barrier_subset_stats(pairs_by_subset):
    stats = {}
    for subset, pairs in pairs_by_subset.items():
        n = len(pairs)
        mad = sum(abs(m - r) for m, r in pairs) / n
        rms = math.sqrt(sum((m - r) ** 2 for m, r in pairs) / n)
        stats[subset] = {"n": n, "MAD": mad, "RMS": rms}
    return stats


# ------------------------------------------------------------------ parameter vector <-> override JSON


def build_override(param_defs, x):
    """Turn a flat parameter vector into the sparse {"rev":{...},"gen":{...},"tables":{...}} doc.

    A dotted name like "rev.p_over.1" becomes override["rev"]["p_over"]["1"] = value.
    """
    override = {}
    for pdef, val in zip(param_defs, x):
        parts = pdef["name"].split(".")
        node = override
        for part in parts[:-1]:
            node = node.setdefault(part, {})
        node[parts[-1]] = float(val)
    return override


def p0_vector(param_defs):
    return np.array([p["p0"] for p in param_defs], dtype=float)


def bounds(param_defs):
    lo = np.array([p["lo"] for p in param_defs], dtype=float)
    hi = np.array([p["hi"] for p in param_defs], dtype=float)
    return lo, hi


def scales(param_defs):
    return np.array([p["scale"] for p in param_defs], dtype=float)


# ------------------------------------------------------------------ running curcuma


def write_system_xyz(system, path):
    lines = []
    for p in system.points:
        lines.append(str(len(p.atoms)))
        lines.append(f"{system.name} {p.label}")
        for sym, x, y, z in p.atoms:
            lines.append(f"{sym} {x:.10f} {y:.10f} {z:.10f}")
    path.write_text("\n".join(lines) + "\n")


def run_system(curcuma, workdir, system, override_path):
    """Run one curcuma batch single point over a system's multi-frame XYZ; return the frame list."""
    xyz_path = workdir / f"{system.name}.xyz"
    out_path = workdir / f"{system.name}.jsonl"
    cmd = [
        str(curcuma), "-sp", str(xyz_path), "-batch", "true", "-batch_out", str(out_path),
        "-batch_reuse_topology", "true", "-method", "revgfnff",
        "-gfnff.param_file", str(override_path),
        "-charge", str(system.charge), "-spin", str(system.mult - 1),
        "-gradient", "true", "-threads", "1", "-verbosity", "0", "-no_bmt",
    ]
    if system.topology_mode == "react":
        cmd += ["-gfnff.topology_mode", "react"]
    try:
        subprocess.run(cmd, cwd=str(workdir), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                        timeout=300, check=False)
    except subprocess.TimeoutExpired:
        return []
    if not out_path.exists():
        return []
    frames = []
    for line in out_path.read_text().splitlines():
        line = line.strip()
        if line:
            try:
                frames.append(json.loads(line))
            except json.JSONDecodeError:
                frames.append({"error": "unparseable JSONL line"})
    return frames


# ------------------------------------------------------------------ evaluation / loss


@dataclass
class Evaluation:
    loss: float
    residuals: np.ndarray
    class_stats: dict          # cls -> {"n_E", "rms_E", "n_G", "rms_G"}
    guard_d: dict              # class-D dE-RMS guard, drives a loss penalty: {"baseline","current","limit","ok"} or {}
    guard_grad: dict           # old class-A gradient-RMS guard, REPORTED ONLY (no penalty): same shape or {}
    n_failed_systems: int
    n_failed_frames: int
    barrier_stats: dict = field(default_factory=dict)   # subset -> {"n", "MAD", "RMS"} (Claude Generated, Sep 2026, WP3)


class FitContext:
    def __init__(self, param_defs, dataset_defs, lam, guard_factor, systems, workdir, curcuma, jobs,
                 guard_d_systems=None, grad_guard_factor=1.5, barrier_reactions=None, barrier_groups=None):
        self.param_defs = param_defs
        self.dataset_defs = dataset_defs
        self.lam = lam
        self.guard_factor = guard_factor            # class-D dE-RMS guard threshold (drives the penalty)
        self.grad_guard_factor = grad_guard_factor   # old class-A gradient guard threshold (reported only)
        self.systems = systems
        self.guard_d_systems = guard_d_systems or []
        # Claude Generated (Sep 2026, WP3): GMTKN55 barrier dataset -- reactions (stoichiometry +
        # reference) and the (charge, spin) batching groups that supply their single points.
        self.barrier_reactions = barrier_reactions or []
        self.barrier_groups = barrier_groups or []
        self.workdir = workdir
        self.curcuma = curcuma
        self.jobs = jobs
        self.cache = {}
        self.n_evals = 0

        # systems actually run every evaluation: the loss systems plus whatever class-D
        # guard systems aren't already part of them (avoids double-running class D if a
        # config ever puts "D" in its own datasets).
        self._run_systems = list(self.systems)
        seen = {s.name for s in self._run_systems}
        for s in self.guard_d_systems:
            if s.name not in seen:
                self._run_systems.append(s)
                seen.add(s.name)

        self.grad_guard_baseline = None
        if any(s.cls == "A" for s in systems):
            self.grad_guard_baseline = self._compute_grad_guard_baseline()
        self.d_guard_baseline = None
        if self.guard_d_systems:
            self.d_guard_baseline = self._compute_d_guard_baseline()

    def _run_all(self, systems, override_path):
        results = {}
        with ThreadPoolExecutor(max_workers=self.jobs) as ex:
            futs = {ex.submit(run_system, self.curcuma, self.workdir, s, override_path): s for s in systems}
            for fut in as_completed(futs):
                s = futs[fut]
                results[s.name] = fut.result()
        return results

    def _run_barrier_groups(self, override_path):
        """Claude Generated (Sep 2026, WP3): run every (charge, spin) barrier bucket, return
        {structure_key: energy_eh} for every structure that computed cleanly."""
        results = {}
        with ThreadPoolExecutor(max_workers=self.jobs) as ex:
            futs = {ex.submit(run_barrier_group, self.curcuma, self.workdir, g, override_path): g
                    for g in self.barrier_groups}
            for fut in as_completed(futs):
                g = futs[fut]
                results[g.name] = fut.result()
        energy_by_key = {}
        n_failed = 0
        for g in self.barrier_groups:
            frames = results.get(g.name)
            if not frames or len(frames) != len(g.keys):
                n_failed += len(g.keys)
                continue
            for k, fr in zip(g.keys, frames):
                if "error" in fr or "energy_eh" not in fr:
                    n_failed += 1
                    continue
                energy_by_key[k] = fr["energy_eh"]
        return energy_by_key, n_failed

    def _compute_grad_guard_baseline(self):
        """Gradient RMS (Eh/Angstrom) at every class-A topology frame, DEFAULT parameters.

        Reported only (see Evaluation.guard_grad) -- kept as a second number, no longer
        the guard that feeds the loss penalty (that is now the class-D dE-RMS guard below).
        """
        empty = self.workdir / "override_baseline_grad.json"
        empty.write_text("{}")
        a_systems = [s for s in self.systems if s.cls == "A"]
        results = self._run_all(a_systems, empty)
        sq, n = 0.0, 0
        for s in a_systems:
            frames = results.get(s.name)
            if frames and "gradient_eh_ang" in frames[0]:
                g = np.array(frames[0]["gradient_eh_ang"], dtype=float).reshape(-1)
                sq += float(np.dot(g, g))
                n += g.size
        return math.sqrt(sq / n) if n else None

    def _compute_d_guard_baseline(self):
        """rms(dE_model) over every class-D point, DEFAULT parameters.

        Self-consistent (model vs its own topology frame, no reference needed): this guards
        against a trial parameter set blowing up the off-equilibrium (MD-snapshot) part of
        the PES relative to the current, already-tested defaults.
        """
        empty = self.workdir / "override_baseline_d.json"
        empty.write_text("{}")
        results = self._run_all(self.guard_d_systems, empty)
        rms = self._d_rms_from_frames(results)
        return rms

    def _d_rms_from_frames(self, frames_by_system):
        sq, n = 0.0, 0
        for s in self.guard_d_systems:
            frames = frames_by_system.get(s.name)
            if not frames or len(frames) != len(s.points) or "error" in frames[0] or "energy_eh" not in frames[0]:
                continue
            e0 = frames[0]["energy_eh"]
            for fr in frames:
                if "error" in fr or "energy_eh" not in fr:
                    continue
                dE = (fr["energy_eh"] - e0) * HARTREE_TO_KCAL
                sq += dE * dE
                n += 1
        return math.sqrt(sq / n) if n else None

    def evaluate(self, x):
        key = tuple(round(float(v), 10) for v in x)
        if key in self.cache:
            return self.cache[key]
        self.n_evals += 1
        override_path = self.workdir / "override_iter.json"
        override_path.write_text(json.dumps(build_override(self.param_defs, x)))
        frames_by_system = self._run_all(self._run_systems, override_path)

        class_pts = {}   # cls -> list of (dE_model, dE_ref)
        class_grad = {}  # cls -> list of (diff_array,) -- kcal/mol/A
        a_topo_sq, a_topo_n = 0.0, 0
        n_failed_systems = 0
        n_failed_frames = 0

        for s in self.systems:
            frames = frames_by_system.get(s.name)
            if not frames or len(frames) != len(s.points) or "error" in frames[0] or "energy_eh" not in frames[0]:
                n_failed_systems += 1
                continue
            e0_model = frames[0]["energy_eh"]
            e0_ref = s.points[0].energy_eh
            if s.cls == "A" and "gradient_eh_ang" in frames[0]:
                g = np.array(frames[0]["gradient_eh_ang"], dtype=float).reshape(-1)
                a_topo_sq += float(np.dot(g, g))
                a_topo_n += g.size
            for pt, fr in zip(s.points, frames):
                if "error" in fr or "energy_eh" not in fr:
                    n_failed_frames += 1
                    continue
                dE_model = (fr["energy_eh"] - e0_model) * HARTREE_TO_KCAL
                dE_ref = (pt.energy_eh - e0_ref) * HARTREE_TO_KCAL
                class_pts.setdefault(s.cls, []).append((dE_model, dE_ref))
                if "gradient_eh_ang" in fr and pt.gradient_eh_ang is not None:
                    gm = np.array(fr["gradient_eh_ang"], dtype=float).reshape(-1) * HARTREE_TO_KCAL
                    g_ref = np.array(pt.gradient_eh_ang, dtype=float) * HARTREE_TO_KCAL
                    class_grad.setdefault(s.cls, []).append(gm - g_ref)

        class_stats = {}
        for cls in sorted(set(class_pts) | set(class_grad)):
            pairs = class_pts.get(cls, [])
            n_e = len(pairs)
            rms_e = math.sqrt(sum((m - r) ** 2 for m, r in pairs) / n_e) if n_e else None
            diffs = class_grad.get(cls, [])
            n_g = len(diffs)
            rms_g = math.sqrt(sum(float(np.dot(d, d)) / d.size for d in diffs) / n_g) if n_g else None
            class_stats[cls] = {"n_E": n_e, "rms_E": rms_e, "n_G": n_g, "rms_G": rms_g}

        # Claude Generated (Sep 2026, WP3): GMTKN55 barrier reactions -- one curcuma batch call
        # per (charge, spin) bucket, then the reaction/barrier sums from those single points.
        n_failed_barrier = 0
        barrier_pairs_by_subset = {}
        if self.barrier_groups:
            energy_by_key, n_failed_barrier = self._run_barrier_groups(override_path)
            barrier_pairs_by_subset = barrier_pairs_from_energies(self.barrier_reactions, energy_by_key)
        barrier_stats = barrier_subset_stats(barrier_pairs_by_subset)

        residual_chunks = []
        for ds in self.dataset_defs:
            if "barriers" in ds:
                w_r = ds.get("weight_R", 1.0)
                pairs = [pr for subset in ds["barriers"] for pr in barrier_pairs_by_subset.get(subset, [])]
                if pairs and w_r:
                    n_r = len(pairs)
                    residual_chunks.append(np.array([(m - r) for m, r in pairs]) * math.sqrt(w_r / n_r))
                continue
            w_e, w_g = ds.get("weight_E", 1.0), ds.get("weight_G", 0.0)
            pairs = [pr for c in ds["classes"] for pr in class_pts.get(c, [])]
            if pairs and w_e:
                n_e = len(pairs)
                residual_chunks.append(np.array([(m - r) for m, r in pairs]) * math.sqrt(w_e / n_e))
            diffs = [d for c in ds["classes"] for d in class_grad.get(c, [])]
            if diffs and w_g:
                n_g = len(diffs)
                for d in diffs:
                    residual_chunks.append(d * math.sqrt(w_g / (n_g * d.size)))

        if self.lam:
            p0 = p0_vector(self.param_defs)
            sc = scales(self.param_defs)
            residual_chunks.append(math.sqrt(self.lam) * (np.asarray(x) - p0) / sc)

        # class-D dE-RMS guard (drives the penalty)
        guard_d = {}
        if self.d_guard_baseline is not None:
            current_d = self._d_rms_from_frames(frames_by_system)
            if current_d is not None:
                limit = self.guard_factor * self.d_guard_baseline
                excess = max(0.0, current_d - limit)
                guard_d = {"baseline": self.d_guard_baseline, "current": current_d, "limit": limit,
                           "ok": excess == 0.0}
                residual_chunks.append(np.array([math.sqrt(1e3) * excess]))

        # old class-A gradient guard (reported only -- no residual appended)
        guard_grad = {}
        if self.grad_guard_baseline is not None:
            current_g = math.sqrt(a_topo_sq / a_topo_n) if a_topo_n else None
            if current_g is not None:
                limit = self.grad_guard_factor * self.grad_guard_baseline
                guard_grad = {"baseline": self.grad_guard_baseline, "current": current_g, "limit": limit,
                              "ok": current_g <= limit}

        residuals = np.concatenate(residual_chunks) if residual_chunks else np.zeros(1)
        loss = float(np.dot(residuals, residuals))
        result = Evaluation(loss, residuals, class_stats, guard_d, guard_grad,
                             n_failed_systems, n_failed_frames + n_failed_barrier, barrier_stats)
        self.cache[key] = result
        return result


# ------------------------------------------------------------------ optimizers (numpy only)


def levenberg_marquardt(ctx, x0, lo, hi, sc, max_iter, log):
    x = np.clip(x0.copy(), lo, hi)
    ev = ctx.evaluate(x)
    loss = ev.loss
    history = [loss]
    lam = 1e-2
    n = len(x)
    for it in range(max_iter):
        # forward-difference Jacobian of the residual vector
        r0 = ev.residuals
        cols = []
        for j in range(n):
            h = 1e-3 * sc[j]
            xp = x.copy()
            xp[j] = min(xp[j] + h, hi[j])
            if xp[j] == x[j]:  # clipped to no-op at the upper bound; step down instead
                xp[j] = max(x[j] - h, lo[j])
                h = xp[j] - x[j]
            evp = ctx.evaluate(xp)
            rp = evp.residuals
            m = min(len(rp), len(r0))
            cols.append((rp[:m] - r0[:m]) / h)
        m = min(len(c) for c in cols)
        J = np.column_stack([c[:m] for c in cols])
        r0m = r0[:m]
        JTJ = J.T @ J
        JTr = J.T @ r0m

        accepted = False
        for _ in range(12):
            try:
                delta = np.linalg.solve(JTJ + lam * np.eye(n), -JTr)
            except np.linalg.LinAlgError:
                lam *= 10
                continue
            x_new = np.clip(x + delta, lo, hi)
            ev_new = ctx.evaluate(x_new)
            if ev_new.loss < loss:
                x, ev, loss = x_new, ev_new, ev_new.loss
                lam = max(lam / 3, 1e-12)
                accepted = True
                break
            lam *= 10
            if lam > 1e10:
                break
        history.append(loss)
        log(f"  LM iter {it + 1}/{max_iter}: loss={loss:.6g} lambda={lam:.2e} accepted={accepted} evals={ctx.n_evals}")
        if not accepted:
            break
        if len(history) >= 2 and history[-2] > 0 and abs(history[-2] - history[-1]) / history[-2] < 1e-4:
            log("  LM: relative loss change below 1e-4, stopping")
            break
    return x, ev, history


def nelder_mead(ctx, x0, lo, hi, sc, max_iter, log):
    n = len(x0)

    def f(x):
        return ctx.evaluate(np.clip(x, lo, hi)).loss

    simplex = [x0.copy()]
    for i in range(n):
        xi = x0.copy()
        xi[i] = np.clip(xi[i] + sc[i], lo[i], hi[i])
        simplex.append(xi)
    fvals = [f(p) for p in simplex]
    history = [min(fvals)]
    alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, 0.5
    for it in range(max_iter):
        order = np.argsort(fvals)
        simplex = [simplex[i] for i in order]
        fvals = [fvals[i] for i in order]
        centroid = np.mean(simplex[:-1], axis=0)
        xr = np.clip(centroid + alpha * (centroid - simplex[-1]), lo, hi)
        fr = f(xr)
        if fvals[0] <= fr < fvals[-2]:
            simplex[-1], fvals[-1] = xr, fr
        elif fr < fvals[0]:
            xe = np.clip(centroid + gamma * (xr - centroid), lo, hi)
            fe = f(xe)
            simplex[-1], fvals[-1] = (xe, fe) if fe < fr else (xr, fr)
        else:
            xc = np.clip(centroid + rho * (simplex[-1] - centroid), lo, hi)
            fc = f(xc)
            if fc < fvals[-1]:
                simplex[-1], fvals[-1] = xc, fc
            else:
                best = simplex[0]
                for i in range(1, len(simplex)):
                    simplex[i] = np.clip(best + sigma * (simplex[i] - best), lo, hi)
                    fvals[i] = f(simplex[i])
        history.append(min(fvals))
        log(f"  NM iter {it + 1}/{max_iter}: loss={history[-1]:.6g} evals={ctx.n_evals}")
        if len(history) >= 2 and history[-2] > 0 and abs(history[-2] - history[-1]) / history[-2] < 1e-4:
            log("  NM: relative loss change below 1e-4, stopping")
            break
    best_i = int(np.argmin(fvals))
    x_best = np.clip(simplex[best_i], lo, hi)
    return x_best, ctx.evaluate(x_best), history


# ------------------------------------------------------------------ reporting


def format_class_stats(class_stats):
    lines = []
    for cls, s in sorted(class_stats.items()):
        e = f"{s['rms_E']:.4f}" if s["rms_E"] is not None else "n/a"
        g = f"{s['rms_G']:.4f}" if s["rms_G"] is not None else "n/a"
        lines.append(f"    class {cls}: n_E={s['n_E']:4d} rms_dE={e:>9s} kcal/mol   "
                      f"n_G={s['n_G']:4d} rms_grad={g:>9s} kcal/mol/A")
    return "\n".join(lines)


def format_barrier_stats(barrier_stats):
    """Claude Generated (Sep 2026, WP3): per-subset MAD/RMS of GMTKN55 barrier reactions."""
    lines = []
    for subset, s in sorted(barrier_stats.items()):
        lines.append(f"    barrier {subset:10s}: n={s['n']:4d} MAD={s['MAD']:8.2f} kcal/mol   RMS={s['RMS']:8.2f} kcal/mol")
    return "\n".join(lines)


def format_guards(ev):
    lines = []
    if ev.guard_d:
        g = ev.guard_d
        lines.append(f"    guard (class D dE-RMS, penalized): current {g['current']:.4f} kcal/mol, "
                      f"baseline {g['baseline']:.4f}, limit {g['limit']:.4f}, ok={g['ok']}")
    if ev.guard_grad:
        g = ev.guard_grad
        lines.append(f"    guard (class A grad-RMS, reported only): current {g['current']:.6f} Eh/A, "
                      f"baseline {g['baseline']:.6f}, limit {g['limit']:.6f}, ok={g['ok']}")
    return "\n".join(lines)


# rev.bo2_* bounds are the original wide ones: the repulsion blend now runs on its own
# switch (rev_bo4_center/width, fixed at 1.5/-10 so revgfnff == gfnff at equilibrium),
# so rev.bo2_* only shapes the over-coordination term E_over and no longer needs the
# tight, near-fixed bounds it had while it also drove the repulsion blend.
EXAMPLE_CONFIG = {
    "parameters": [
        {"name": "rev.bo_center", "p0": 2.0, "lo": 1.5, "hi": 3.0, "scale": 0.2},
        {"name": "rev.bo_width", "p0": -7.5, "lo": -15.0, "hi": -3.0, "scale": 1.0},
        {"name": "rev.bo2_center", "p0": 1.4, "lo": 1.2, "hi": 1.8, "scale": 0.1},
        {"name": "rev.bo2_width", "p0": -8.0, "lo": -16.0, "hi": -4.0, "scale": 1.0},
        {"name": "rev.p_over.1", "p0": 0.3, "lo": 0.02, "hi": 2.0, "scale": 0.1},
        {"name": "rev.p_over.6", "p0": 0.3, "lo": 0.02, "hi": 2.0, "scale": 0.1},
        {"name": "rev.p_over.7", "p0": 0.3, "lo": 0.02, "hi": 2.0, "scale": 0.1},
        {"name": "rev.p_over.8", "p0": 0.3, "lo": 0.02, "hi": 2.0, "scale": 0.1},
        {"name": "rev.over_shift", "p0": 0.5, "lo": 0.0, "hi": 1.0, "scale": 0.1},
    ],
    "datasets": [
        {"classes": ["A"], "weight_E": 1.0, "weight_G": 0.1},
        {"classes": ["C"], "weight_E": 1.0, "weight_G": 0.1},
    ],
    "lambda": 0.01,
    "guard_factor": 1.10,       # class-D dE-RMS guard (drives the penalty)
    "grad_guard_factor": 1.5,   # old class-A gradient guard (reported only)
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, help="parameter/dataset config JSON")
    ap.add_argument("--write-example-config", action="store_true",
                     help="write the stage-1 example config to --config (default revgfnff_fit_example.json) and exit")
    ap.add_argument("--workdir", type=Path, default=REPO / "test_cases" / "revgfnff" / "fit_work",
                     help="scratch directory for per-system XYZ/JSONL/override files")
    ap.add_argument("--out", type=Path, help="output directory for fit_result.json / override_fitted.json")
    ap.add_argument("--jobs", type=int, default=8, help="parallel curcuma processes (each -threads 1)")
    ap.add_argument("--max-iter", type=int, default=20)
    ap.add_argument("--method", choices=["lm", "nm"], default="lm")
    ap.add_argument("--evaluate-only", action="store_true", help="print per-class RMS at p0 and exit")
    ap.add_argument("--dry-run", action="store_true", help="smoke test: evaluate p0 on 2 systems only")
    ap.add_argument("--curcuma", type=Path, default=CURCUMA)
    ap.add_argument("--topology", choices=["static", "react"], default="react",
                     help="topology mode for classes A/C/E/L (class D is always static): "
                          "react adds -gfnff.topology_mode react and orders frames as a "
                          "monotonic reaction coordinate so bond formation/breaking can fire")
    args = ap.parse_args()

    if args.write_example_config:
        out = args.config or Path("revgfnff_fit_example.json")
        out.write_text(json.dumps(EXAMPLE_CONFIG, indent=2) + "\n")
        print(f"wrote example config to {out}")
        return

    if not args.config:
        ap.error("--config is required (or --write-example-config)")
    config = json.loads(args.config.read_text())
    param_defs = config["parameters"]
    dataset_defs = config["datasets"]
    lam = config.get("lambda", 0.0)
    guard_factor = config.get("guard_factor", 1.10)
    grad_guard_factor = config.get("grad_guard_factor", 1.5)

    classes_needed = sorted({c for d in dataset_defs for c in d.get("classes", [])})
    # Claude Generated (Sep 2026, WP3): GMTKN55 barrier subsets requested by any dataset entry
    # of the form {"barriers": [...], "weight_R": ...} -- see load_barrier_data().
    barrier_subsets_needed = sorted({s for d in dataset_defs for s in d.get("barriers", [])})

    systems = []
    if classes_needed:
        print(f"loading reference data for classes {classes_needed} (topology={args.topology}) ...")
        systems = load_systems(classes_needed, topology_mode=args.topology)
        print(f"loaded {len(systems)} reference systems "
              f"({sum(len(s.points) for s in systems)} points)")

    barrier_reactions, barrier_groups = [], []
    if barrier_subsets_needed and not args.dry_run:
        print(f"loading GMTKN55 barrier reactions for {barrier_subsets_needed} ...")
        barrier_reactions, barrier_groups = load_barrier_data(barrier_subsets_needed)
        print(f"loaded {len(barrier_reactions)} barrier reactions "
              f"({len(barrier_groups)} charge/spin batches, "
              f"{sum(len(g.keys) for g in barrier_groups)} structures)")

    if not systems and not barrier_reactions:
        print("no reference systems or barrier reactions found for the requested datasets -- aborting")
        sys.exit(1)

    # class-D guard systems: always loaded (regardless of the requested datasets) unless
    # this is just the --dry-run mechanics smoke test, since the guard needs the whole
    # class-D set to establish its default-parameter baseline. D is always static (part 1).
    guard_d_systems = []
    if not args.dry_run:
        guard_d_systems = load_systems(["D"], topology_mode="static")
        if guard_d_systems:
            print(f"loaded {len(guard_d_systems)} class-D guard systems "
                  f"({sum(len(s.points) for s in guard_d_systems)} points, static topology)")

    if args.dry_run:
        picked, seen_cls = [], set()
        for s in systems:
            if s.cls not in seen_cls:
                picked.append(s)
                seen_cls.add(s.cls)
            if len(picked) == 2:
                break
        if len(picked) < 2:
            picked = systems[:2]
        systems = picked
        print(f"--dry-run: restricted to {[s.name for s in systems]}")

    args.workdir.mkdir(parents=True, exist_ok=True)
    # A stale GFN-FF topology cache next to a reused batch XYZ can silently replay an
    # earlier evaluation's parameters (the fingerprint is supposed to include the table
    # hash, but a workdir surviving across curcuma rebuilds/parameter schema changes is
    # not worth the risk) -- always start a run from a clean topology cache.
    stale_topo = list(args.workdir.glob("*.topo.json"))
    for f in stale_topo:
        f.unlink()
    if stale_topo:
        print(f"removed {len(stale_topo)} stale .topo.json cache file(s) from {args.workdir}")
    for s in systems + guard_d_systems:
        write_system_xyz(s, args.workdir / f"{s.name}.xyz")
    for g in barrier_groups:
        write_barrier_group_xyz(g, args.workdir / f"{g.name}.xyz")

    ctx = FitContext(param_defs, dataset_defs, lam, guard_factor, systems, args.workdir, args.curcuma, args.jobs,
                      guard_d_systems=guard_d_systems, grad_guard_factor=grad_guard_factor,
                      barrier_reactions=barrier_reactions, barrier_groups=barrier_groups)
    x0 = p0_vector(param_defs)

    t0 = time.time()
    ev0 = ctx.evaluate(x0)
    print(f"\np0 evaluation: loss={ev0.loss:.6g}  "
          f"failed systems={ev0.n_failed_systems}/{len(systems)}  failed frames={ev0.n_failed_frames}")
    print(format_class_stats(ev0.class_stats))
    if ev0.barrier_stats:
        print(format_barrier_stats(ev0.barrier_stats))
    print(format_guards(ev0))

    if args.evaluate_only or args.dry_run:
        print(f"\nevals={ctx.n_evals}  wall={time.time() - t0:.1f}s")
        return

    lo, hi = bounds(param_defs)
    sc = scales(param_defs)

    def log(msg):
        print(msg)

    print(f"\nfitting with {args.method} (max_iter={args.max_iter}, jobs={args.jobs}) ...")
    if args.method == "lm":
        x_final, ev_final, history = levenberg_marquardt(ctx, x0, lo, hi, sc, args.max_iter, log)
    else:
        x_final, ev_final, history = nelder_mead(ctx, x0, lo, hi, sc, args.max_iter, log)
    wall = time.time() - t0

    print(f"\nfinal evaluation: loss={ev_final.loss:.6g}  "
          f"failed systems={ev_final.n_failed_systems}/{len(systems)}  failed frames={ev_final.n_failed_frames}")
    print(format_class_stats(ev_final.class_stats))
    if ev_final.barrier_stats:
        print(format_barrier_stats(ev_final.barrier_stats))
    print(format_guards(ev_final))

    print("\nparameters (name: p0 -> final):")
    for pdef, v0, vf in zip(param_defs, x0, x_final):
        print(f"    {pdef['name']:20s} {v0: .6g} -> {vf: .6g}")

    print(f"\nevaluations: {ctx.n_evals}   wall time: {wall:.1f} s")

    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)
        result = {
            "method": args.method,
            "max_iter": args.max_iter,
            "n_evaluations": ctx.n_evals,
            "wall_s": wall,
            "loss_history": history,
            "parameters": [
                {"name": p["name"], "p0": float(v0), "final": float(vf), "lo": p["lo"], "hi": p["hi"], "scale": p["scale"]}
                for p, v0, vf in zip(param_defs, x0, x_final)
            ],
            "class_stats_before": ev0.class_stats,
            "class_stats_after": ev_final.class_stats,
            "barrier_stats_before": ev0.barrier_stats,
            "barrier_stats_after": ev_final.barrier_stats,
            "guard_d_before": ev0.guard_d,
            "guard_d_after": ev_final.guard_d,
            "guard_grad_before": ev0.guard_grad,
            "guard_grad_after": ev_final.guard_grad,
            "loss_before": ev0.loss,
            "loss_after": ev_final.loss,
            "topology": args.topology,
        }
        (args.out / "fit_result.json").write_text(json.dumps(result, indent=2))
        (args.out / "override_fitted.json").write_text(json.dumps(build_override(param_defs, x_final), indent=2))
        print(f"\nwrote {args.out / 'fit_result.json'} and {args.out / 'override_fitted.json'}")


if __name__ == "__main__":
    main()
