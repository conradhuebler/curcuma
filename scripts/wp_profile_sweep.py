#!/usr/bin/env python3
"""
WP-P2: GFN-FF Cross-System Benchmark Sweep

Drives curcuma over (system x path x mode) combinations, collects the
per-step timing_ms blocks from <basename>.diag.jsonl files produced by
WP-P1's md_diagnostics_timing PARAM, and aggregates them into:

  <out_dir>/wp_profile_results.csv     -- one row per JSONL snapshot
  <out_dir>/wp_profile_aggregate.csv   -- one row per config, mean over records
  <out_dir>/wp_profile_summary.md      -- human-readable summary tables

Usage:
  python3 scripts/wp_profile_sweep.py --systems caffeine,polymer,mixture \\
       --paths cpu,gpu --modes baseline,static_all,cutoff_auto,both \\
       --max-time 200 --output /tmp/wp_p2_sweep

Stdlib-only (json, subprocess, argparse, statistics, csv, pathlib, time).
"""

import argparse
import csv
import json
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent

SYSTEMS = {
    "caffeine": REPO_ROOT / "test_cases/molecules/larger/caffeine.xyz",
    "polymer":  REPO_ROOT / "test_cases/molecules/larger/polymer.xyz",
    "mixture":  REPO_ROOT / "test_cases/molecules/larger/mixture.xyz",
}

MODES = {
    "baseline":    {},
    "static_all":  {"static_all": True},
    "cutoff_auto": {"eeq_distance_cutoff_auto": True},
    "both":        {"static_all": True, "eeq_distance_cutoff_auto": True},
}

TIMING_FIELDS = [
    "cn", "eeq_topo", "cnf", "dcn", "d4_gw", "eeq_solve",
    "charge_dist", "total", "ff_total", "integrator",
    "hbxb_update", "step_total",
]
GPU_FIELDS = [
    "dispersion", "bonds", "angles", "dihedrals", "inversions",
    "coulomb", "hbond", "xbond", "atm", "batm",
    "bonded_rep", "nonbonded_rep", "stors",
]


def find_curcuma(explicit: str | None) -> Path:
    """Pick the curcuma binary: explicit path, ./curcuma, then release/curcuma."""
    candidates = []
    if explicit:
        candidates.append(Path(explicit))
    candidates.extend([
        Path.cwd() / "curcuma",
        REPO_ROOT / "release" / "curcuma",
        REPO_ROOT / "release-icx" / "curcuma",
    ])
    for c in candidates:
        if c.is_file() and (c.stat().st_mode & 0o111):
            return c.resolve()
    raise FileNotFoundError(f"No curcuma binary found in: {candidates}")


def detect_cuda(curcuma: Path) -> bool:
    """Best-effort detection: curcuma --help typically lists 'gpu' option when USE_CUDA=ON."""
    try:
        out = subprocess.run([str(curcuma), "--help"], capture_output=True, text=True, timeout=15)
        text = (out.stdout + out.stderr).lower()
        return "cuda" in text or "use_cuda" in text or "-gpu" in text
    except Exception:
        return False


def run_one(curcuma, system_path, label, path, mode_flags, out_dir, max_time, dt, threads, timeout):
    """Run a single MD config; returns (wall_seconds, exit_code, jsonl_path)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = {
        "gfnff": dict(mode_flags),
        "simplemd": {
            "method": "gfnff", "max_time": max_time, "dt": dt,
            "thermostat": "csvr", "dump_frequency": 10,
            "md_diagnostics": True, "md_diagnostics_timing": True,
            "seed": 42, "threads": threads,
        },
        "global": {
            "method": "gfnff", "threads": threads,
        },
    }
    if path == "gpu":
        cfg["global"]["gpu"] = "cuda"

    cfg_path = out_dir / "cfg.json"
    cfg_path.write_text(json.dumps(cfg, indent=2))

    # Stage the input XYZ next to cfg so curcuma writes outputs into the dir
    local_xyz = out_dir / system_path.name
    if not local_xyz.exists() or local_xyz.read_bytes() != system_path.read_bytes():
        shutil.copy(system_path, local_xyz)

    t0 = time.monotonic()
    try:
        cp = subprocess.run(
            [str(curcuma), "-md", local_xyz.name, "-import_config", cfg_path.name],
            cwd=str(out_dir), capture_output=True, text=True, timeout=timeout)
        wall = time.monotonic() - t0
        (out_dir / "stdout.log").write_text(cp.stdout)
        (out_dir / "stderr.log").write_text(cp.stderr)
        rc = cp.returncode
    except subprocess.TimeoutExpired as e:
        wall = time.monotonic() - t0
        (out_dir / "stdout.log").write_text(e.stdout or "")
        (out_dir / "stderr.log").write_text((e.stderr or "") + f"\n[TIMEOUT after {timeout}s]\n")
        rc = -1

    diag_jsonl = out_dir / f"{system_path.stem}.diag.jsonl"
    return wall, rc, diag_jsonl


def parse_jsonl(path: Path):
    if not path.is_file():
        return []
    out = []
    for line in path.open():
        line = line.strip()
        if line:
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                pass
    return out


def aggregate_records(records, skip_first=1):
    """Returns dict {field: mean} for top-level timing_ms + nested gpu fields."""
    if len(records) <= skip_first:
        return {}
    usable = records[skip_first:]
    agg = {}
    # Top-level timing fields
    for f in TIMING_FIELDS:
        vals = [r.get("timing_ms", {}).get(f) for r in usable]
        vals = [v for v in vals if isinstance(v, (int, float))]
        if vals:
            agg[f] = statistics.mean(vals)
    # GPU sub-block
    for f in GPU_FIELDS:
        vals = [r.get("timing_ms", {}).get("gpu", {}).get(f) for r in usable]
        vals = [v for v in vals if isinstance(v, (int, float))]
        if vals:
            agg[f"gpu_{f}"] = statistics.mean(vals)
    return agg


def write_results_csv(rows, path: Path):
    if not rows:
        return
    fieldnames = sorted({k for r in rows for k in r.keys()})
    head_cols = ["system", "path", "mode", "threads", "step", "time_fs"]
    head_cols = [c for c in head_cols if c in fieldnames]
    rest = [c for c in fieldnames if c not in head_cols]
    cols = head_cols + rest
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_summary_md(aggregates, summary_path, max_time, dt, threads):
    """Aggregates is list of dicts: {system, path, mode, threads, wall, exit_code, <fields>}."""
    summary_path = Path(summary_path)
    by_system_path = {}
    for a in aggregates:
        by_system_path.setdefault((a["system"], a["path"]), []).append(a)

    lines = [
        "# WP-P2 Cross-System Sweep Summary",
        "",
        f"- max_time = {max_time} fs  |  dt = {dt} fs  |  threads = {threads}  |  thermostat = csvr",
        f"- Mean over records 1.. (Record 0 = init step is skipped)",
        "",
    ]

    for (sys_name, path), entries in sorted(by_system_path.items()):
        lines.append(f"## {sys_name} ({path})")
        lines.append("")
        # Pick which columns to display
        show_cols = ["step_total", "dcn", "eeq_solve", "cn", "d4_gw", "ff_total", "integrator"]
        # Find baseline for speedup
        baseline = next((e for e in entries if e["mode"] == "baseline"), None)
        base_step = baseline.get("step_total") if baseline else None

        header = "| Modus      | " + " | ".join(c for c in show_cols) + " | wall (s) | exit | speedup |"
        sep    = "|------------|" + "|".join("-" * (len(c) + 2) for c in show_cols) + "|----------|------|---------|"
        lines.append(header)
        lines.append(sep)
        for e in entries:
            cells = []
            for c in show_cols:
                v = e.get(c)
                cells.append(f"{v:.2f}" if isinstance(v, (int, float)) else "-")
            wall_str = f"{e.get('wall', 0):.1f}" if isinstance(e.get("wall"), (int, float)) else "-"
            exit_str = str(e.get("exit_code", "-"))
            sp_val = e.get("step_total")
            if base_step and isinstance(sp_val, (int, float)) and sp_val > 0:
                speedup = base_step / sp_val
                sp_str = f"{speedup:.2f}x"
            else:
                sp_str = "-"
            lines.append(f"| {e['mode']:<10} | " + " | ".join(cells) + f" | {wall_str} | {exit_str} | {sp_str} |")
        lines.append("")

        # GPU subblock if any
        has_gpu = any("gpu_dispersion" in e or "gpu_coulomb" in e for e in entries)
        if has_gpu:
            gpu_cols = ["gpu_dispersion", "gpu_bonds", "gpu_coulomb", "gpu_hbond", "gpu_atm"]
            lines.append("### GPU per-term timings (ms)")
            lines.append("")
            lines.append("| Modus      | " + " | ".join(c for c in gpu_cols) + " |")
            lines.append("|------------|" + "|".join("-" * (len(c) + 2) for c in gpu_cols) + "|")
            for e in entries:
                cells = [f"{e.get(c, 0):.2f}" if isinstance(e.get(c), (int, float)) else "-" for c in gpu_cols]
                lines.append(f"| {e['mode']:<10} | " + " | ".join(cells) + " |")
            lines.append("")

    summary_path.write_text("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--curcuma", help="Path to curcuma binary (default: discover)")
    ap.add_argument("--systems", default="caffeine,polymer,mixture",
                    help="Comma-separated system names")
    ap.add_argument("--paths", default="cpu,gpu",
                    help="Comma-separated path labels (cpu, gpu)")
    ap.add_argument("--modes", default=",".join(MODES.keys()),
                    help="Comma-separated mode names")
    ap.add_argument("--max-time", type=float, default=200.0, help="MD max_time in fs")
    ap.add_argument("--dt", type=float, default=1.0, help="MD timestep in fs")
    ap.add_argument("--threads", type=int, default=4, help="CPU threads per run")
    ap.add_argument("--timeout", type=int, default=600, help="Per-run timeout in seconds")
    ap.add_argument("--output", default="/tmp/wp_p2_sweep", help="Output directory")
    ap.add_argument("--thread-sweep", help="Optional comma-separated thread counts; runs each baseline only")
    args = ap.parse_args()

    curcuma = find_curcuma(args.curcuma)
    print(f"[wp_p2] curcuma = {curcuma}", flush=True)

    have_cuda = detect_cuda(curcuma)
    print(f"[wp_p2] USE_CUDA detected: {have_cuda}", flush=True)

    out_root = Path(args.output)
    out_root.mkdir(parents=True, exist_ok=True)

    systems = [s.strip() for s in args.systems.split(",") if s.strip()]
    paths = [p.strip() for p in args.paths.split(",") if p.strip()]
    modes = [m.strip() for m in args.modes.split(",") if m.strip()]

    if not have_cuda and "gpu" in paths:
        print("[wp_p2] No CUDA detected — dropping 'gpu' from paths", flush=True)
        paths = [p for p in paths if p != "gpu"]

    # Validate systems
    valid_systems = []
    for s in systems:
        if s not in SYSTEMS:
            print(f"[wp_p2] WARNING: unknown system {s}, skipping", flush=True)
            continue
        if not SYSTEMS[s].is_file():
            print(f"[wp_p2] WARNING: {SYSTEMS[s]} not found, skipping", flush=True)
            continue
        valid_systems.append(s)
    if not valid_systems:
        print("[wp_p2] No usable systems — aborting", flush=True)
        sys.exit(1)

    thread_sweep = None
    if args.thread_sweep:
        thread_sweep = [int(t) for t in args.thread_sweep.split(",") if t.strip()]

    all_records = []
    aggregates = []

    def configs():
        if thread_sweep:
            # Threading audit: one system, baseline-only, multiple thread counts
            for s in valid_systems:
                for p in paths:
                    for t in thread_sweep:
                        yield s, p, "baseline", t
        else:
            for s in valid_systems:
                for p in paths:
                    for m in modes:
                        yield s, p, m, args.threads

    for sys_name, path, mode, threads in configs():
        system_path = SYSTEMS[sys_name]
        label = f"{sys_name}_{path}_{mode}_t{threads}"
        sub = out_root / label
        print(f"[wp_p2] run {label} ... ", end="", flush=True)
        wall, rc, jsonl = run_one(
            curcuma, system_path, label, path, MODES[mode], sub,
            args.max_time, args.dt, threads, args.timeout)
        print(f"wall={wall:.1f}s rc={rc} jsonl={'yes' if jsonl.is_file() else 'no'}", flush=True)
        records = parse_jsonl(jsonl)
        agg = aggregate_records(records)
        # Append rows
        for idx, rec in enumerate(records):
            row = {
                "system": sys_name, "path": path, "mode": mode,
                "threads": threads, "step": rec.get("step"),
                "time_fs": rec.get("time_fs"),
            }
            for f, v in (rec.get("timing_ms") or {}).items():
                if isinstance(v, (int, float)):
                    row[f] = v
            gpu = (rec.get("timing_ms") or {}).get("gpu") or {}
            for f, v in gpu.items():
                if isinstance(v, (int, float)):
                    row[f"gpu_{f}"] = v
            all_records.append(row)
        # Aggregate
        a = {"system": sys_name, "path": path, "mode": mode, "threads": threads,
             "wall": wall, "exit_code": rc, "n_records": len(records)}
        a.update(agg)
        aggregates.append(a)

    # Write CSVs
    write_results_csv(all_records, out_root / "wp_profile_results.csv")
    write_results_csv(aggregates,   out_root / "wp_profile_aggregate.csv")
    write_summary_md(aggregates, out_root / "wp_profile_summary.md",
                     args.max_time, args.dt, args.threads)

    print(f"[wp_p2] outputs in {out_root}", flush=True)
    print(f"        wp_profile_results.csv    ({len(all_records)} rows)", flush=True)
    print(f"        wp_profile_aggregate.csv  ({len(aggregates)} rows)", flush=True)
    print(f"        wp_profile_summary.md", flush=True)


if __name__ == "__main__":
    main()
