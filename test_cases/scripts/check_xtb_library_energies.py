#!/usr/bin/env python3
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Oct 2026)
"""Regression gate: curcuma (CPU) single points against stored xtb 6.7.1 energies.

Reference: test_cases/reference_data/xtb_library_energies.json (no xtb run, no new reference).
For every structure with status "ok" the first frame is written to a fresh temporary
directory as m.xyz and curcuma runs `-sp m.xyz -method gfnff|gfn2` (4 threads, -no_bmt).

Rules (tolerances are checked against the stored xtb energies):
  GFN-FF: |new - xtb| <= 1e-5 Eh, except the frozen KNOWN_DEVIATIONS table: there
          |(new - xtb) - d_frozen| <= 1e-5 Eh.
          For charge != 0 the run uses -gfnff.frag_charge_model reference.
  GFN2:   natoms <= 500 -> |new - xtb| <= 1e-4 Eh.
          natoms > 500 (peo201-chain) only with --include-large, then as known deviation.

Exit code: 0 = no FAIL, 1 = at least one FAIL, 2 = setup error (binary/files missing).

Usage: check_xtb_library_energies.py --curcuma <binary> [--reference <json>]
           [--structures <dir>] [--only <id>] [--include-large]
"""
import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

GFNFF_TOL = 1e-5
GFNFF_KNOWN_LIMIT = 1e-5
GFN2_TOL = 1e-4
GFN2_LARGE_NATOMS = 500
KNOWN_TOL = 1e-5

# Known deviations curcuma - xtb 6.7.1 (Eh), frozen from the measurement of 2026-10-11 (curcuma at
# 5d34e7ed..573db3ea). They are NOT derived from the reference file, so an edited xtb value cannot
# turn a regression into a "known deviation". The GFN-FF glycoside systems (173/231 atoms) are not
# explained (TODO.md); peo201-chain GFN2 is the D4 ATM cutoff split (+2.6e-5 Eh with
# -xtb.d4_atm_cutoff 40.0). The test guards against a change of the deviation, not against its size.
KNOWN_DEVIATIONS = {
    ("helicene-frames17", "gfnff"): +1.0725e-05,
    ("jl22-bgal.conf6", "gfnff"): +1.281529e-03,
    ("jl22-bgal.xtb-gfn2-opt", "gfnff"): +1.232441e-03,
    ("jl22-bgal.xtb-gfnff-opt", "gfnff"): +8.04045e-04,
    ("macrocycle-bgal", "gfnff"): +3.97418e-04,
    ("macrocycle-bgal.conf2", "gfnff"): +4.08097e-04,
    ("macrocycle-bgal.uff-analytic-grad", "gfnff"): +1.205250e-03,
    ("macrocycle-bgal.uff-numeric-grad", "gfnff"): +5.59316e-04,
    ("peo201-chain", "gfn2"): +1.602879e-03,
}
RUN_TIMEOUT_S = 600

ANSI_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
ENERGY_RE = re.compile(r"Single Point Energy\s*=\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*Eh")


def first_frame_xyz(path):
    """Return the text of the first xyz frame (natoms line, comment line, atom lines)."""
    lines = path.read_text().splitlines(keepends=True)
    natoms = int(lines[0].split()[0])
    if len(lines) < natoms + 2:
        raise ValueError(f"{path.name}: fewer lines than the first frame needs")
    return "".join(lines[: natoms + 2])


def run_curcuma(binary, xyz_text, method, charge, spin):
    """Run one curcuma single point in a fresh directory; return the energy in Eh."""
    workdir = Path(tempfile.mkdtemp(prefix="xtblib_"))
    try:
        (workdir / "m.xyz").write_text(xyz_text)
        cmd = [str(binary), "-sp", "m.xyz", "-method", method, "-threads", "4",
               "-charge", str(charge), "-spin", str(spin), "-no_bmt", "-verbosity", "1"]
        if method == "gfnff" and charge != 0:
            cmd += ["-gfnff.frag_charge_model", "reference"]
        try:
            proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True,
                                  timeout=RUN_TIMEOUT_S)
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"timeout after {RUN_TIMEOUT_S} s")
        text = ANSI_RE.sub("", proc.stdout + proc.stderr)
        matches = ENERGY_RE.findall(text)
        if not matches:
            raise RuntimeError(f"no 'Single Point Energy' line (exit {proc.returncode})")
        return float(matches[-1])
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def fmt_line(sid, method, natoms, new, xtb, de, tol, verdict, mark):
    return (f"{sid:<26} {method:<6} {natoms:>5}  curcuma={new:+.10f}  xtb={xtb:+.10f}  "
            f"dE={de:+.3e}  tol={tol:.1e}  {verdict}{mark}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    here = Path(__file__).resolve().parent
    ap.add_argument("--curcuma", required=True, help="curcuma binary")
    ap.add_argument("--reference", default=str(here.parent / "reference_data" / "xtb_library_energies.json"))
    ap.add_argument("--structures", default=str(here.parent / "structures"))
    ap.add_argument("--only", default=None, help="check only this structure id")
    ap.add_argument("--include-large", action="store_true",
                    help="also run GFN2 on structures with natoms > 500 (known-deviation mode)")
    args = ap.parse_args()

    binary = Path(args.curcuma)
    ref_path = Path(args.reference)
    struct_dir = Path(args.structures)
    if not binary.is_file():
        print(f"ERROR: curcuma binary not found: {binary}", file=sys.stderr)
        return 2
    if not ref_path.is_file():
        print(f"ERROR: reference file not found: {ref_path}", file=sys.stderr)
        return 2
    manifest_path = struct_dir / "manifest.json"
    if not manifest_path.is_file():
        print(f"ERROR: manifest not found: {manifest_path}", file=sys.stderr)
        return 2

    ref = json.loads(ref_path.read_text())["structures"]
    manifest = {s["id"]: s for s in json.loads(manifest_path.read_text())["structures"]}
    if args.only is not None:
        if args.only not in ref:
            print(f"ERROR: --only id not in reference: {args.only}", file=sys.stderr)
            return 2
        ref = {args.only: ref[args.only]}

    checked = passed = known = skipped = failed = 0
    known_list = []
    fail_list = []
    print("xtb 6.7.1 library energies: curcuma (CPU, 4 threads) vs stored xtb")
    for sid in sorted(ref):
        entry = ref[sid]
        natoms = entry["natoms"]
        charge = entry["charge"] or 0
        spin = entry["spin"] or 0
        if entry["status"] != "ok":
            print(f"{sid:<26} -      {natoms:>5}  SKIP (status {entry['status']}: {entry['note'][:60]})")
            skipped += 2
            continue
        meta = manifest.get(sid)
        if meta is None:
            print(f"{sid:<26} -      {natoms:>5}  SKIP (not in manifest)")
            skipped += 2
            continue
        if meta["format"] != "xyz":
            print(f"{sid:<26} -      {natoms:>5}  SKIP (format {meta['format']} not handled)")
            skipped += 2
            continue
        src = struct_dir / meta["file"]
        if not src.is_file():
            print(f"ERROR: structure file missing: {src}", file=sys.stderr)
            return 2
        xyz_text = first_frame_xyz(src)

        # Each method: (name, stored xtb, stored curcuma baseline, tolerance mode)
        methods = []
        if entry["gfnff_xtb"] is None:
            print(f"{sid:<26} gfnff  {natoms:>5}  SKIP (no stored gfnff_xtb)")
            skipped += 1
        else:
            base = entry["gfnff_curcuma"]
            if charge != 0 and entry.get("gfnff_curcuma_fragref") is not None:
                base = entry["gfnff_curcuma_fragref"]
            methods.append(("gfnff", entry["gfnff_xtb"], base, "gfnff"))
        if entry["gfn2_xtb"] is None:
            print(f"{sid:<26} gfn2   {natoms:>5}  SKIP (no stored gfn2_xtb)")
            skipped += 1
        elif natoms > GFN2_LARGE_NATOMS and not args.include_large:
            print(f"{sid:<26} gfn2   {natoms:>5}  SKIP (natoms > {GFN2_LARGE_NATOMS}, use --include-large)")
            skipped += 1
        else:
            methods.append(("gfn2", entry["gfn2_xtb"], entry["gfn2_curcuma"], "gfn2"))

        for name, xtb, stored_cur, rule in methods:
            checked += 1
            try:
                new = run_curcuma(binary, xyz_text, name, charge, spin)
            except RuntimeError as exc:
                failed += 1
                fail_list.append(f"{sid} {name}")
                print(f"{sid:<26} {name:<6} {natoms:>5}  FAIL (curcuma run: {exc})")
                continue
            de = new - xtb
            d_known = KNOWN_DEVIATIONS.get((sid, name))
            if d_known is not None:
                tol = KNOWN_TOL
                ok = abs(de - d_known) <= KNOWN_TOL
                mark = f"  known deviation (frozen d={d_known:+.4e})"
                known += 1
                if ok:
                    known_list.append(f"{sid} {name} d={d_known:+.4e}")
            else:
                tol = GFNFF_TOL if rule == "gfnff" else GFN2_TOL
                ok = abs(de) <= tol
                mark = ""
            print(fmt_line(sid, name, natoms, new, xtb, de, tol, "PASS" if ok else "FAIL", mark))
            if ok:
                passed += 1
            else:
                failed += 1
                fail_list.append(f"{sid} {name}")

    print("")
    print(f"SUMMARY: geprueft={checked} bestanden={passed} (davon bekannt={known}) "
          f"FAIL={failed} uebersprungen={skipped}")
    if fail_list:
        print("FAILED: " + ", ".join(fail_list))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
