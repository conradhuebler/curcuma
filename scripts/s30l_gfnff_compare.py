#!/usr/bin/env python3
# Claude Generated (Jul 2026)
"""S30L host-guest benchmark: compare Curcuma native GFN-FF vs xtb 6.6.1 GFN-FF.

Runs single-point GFN-FF on host (A), guest (B), complex (AB) of all 30 S30L
complexes with both engines on identical xyz geometries (converted once from
the Turbomole coord files, Bohr -> Angstrom), computes association energies
DE = E(AB) - E(A) - E(B) [kcal/mol], and compares:
  - Primary:   curcuma vs xtb            (does native reproduce the reference impl?)
  - Secondary: both vs reference_s30l    (context; GFN-FF does not match lit. values)

Charged systems (23-30) carry .CHRG files; charge is passed to curcuma via
-charge and to xtb via a .CHRG file in the xtb working directory.

Read-only w.r.t. the test set; writes only under _run/.
"""
import csv
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TESTSET = REPO / "test_cases" / "s30l_test_set"
RUNDIR = TESTSET / "_run"
CURCUMA = REPO / "release" / "curcuma"
XTB = Path.home() / "Downloads" / "xtb-6.6.1" / "bin" / "xtb"

BOHR = 0.52917721067          # Bohr -> Angstrom
AU2KCAL = 627.509474           # Hartree -> kcal/mol

# Only --gfnff is needed for a single point (xtb defaults to SP when no --opt).
XTB_ARGS = ["--gfnff", "--sp"]


def read_coord(path):
    """Parse Turbomole $coord (Bohr) -> list of (element, x, y, z) in Angstrom.

    Lines starting with '$' are directives and skipped. Each data line is
    'x y z element' (element lowercase, e.g. 'c','h','n')."""
    atoms = []
    in_coord = False
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("$"):
            in_coord = s.startswith("$coord")
            continue
        if not in_coord:
            continue
        parts = s.split()
        if len(parts) < 4:
            continue
        el = parts[3]
        x, y, z = float(parts[0]) * BOHR, float(parts[1]) * BOHR, float(parts[2]) * BOHR
        atoms.append((el.capitalize(), x, y, z))
    return atoms


def read_charge(folder):
    """Read integer charge from folder/.CHRG if present, else 0."""
    chrg = folder / ".CHRG"
    if chrg.exists():
        # int() tolerates leading '+'/'-'; .CHRG holds e.g. '+1',' -2 '
        return int(chrg.read_text().strip())
    return 0


def write_xyz(path, atoms, comment=""):
    with path.open("w") as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for el, x, y, z in atoms:
            f.write(f"{el} {x:.10f} {y:.10f} {z:.10f}\n")


def parse_curcuma_energy(stdout):
    m = re.search(r"Single Point Energy\s*=\s*(-?\d+\.\d+)\s*Eh", stdout)
    if not m:
        return None
    return float(m.group(1))


def parse_xtb_energy(stdout):
    # xtb 6.6.1 prints "| TOTAL ENERGY  ...  Eh |" in a bordered table, plus
    # a plain "total energy   :  -123.456789 Eh" summary line. Try both.
    m = re.search(r"TOTAL ENERGY[^\n]*?(-?\d+\.\d+)\s*Eh", stdout)
    if not m:
        m = re.search(r"total energy\s*:\s*(-?\d+\.\d+)\s*Eh", stdout)
    if not m:
        # last-ditch: any "Eh" float on a line mentioning energy
        for line in stdout.splitlines():
            if "energy" in line.lower():
                mm = re.search(r"(-?\d+\.\d+)\s*Eh", line)
                if mm:
                    return float(mm.group(1))
    return float(m.group(1)) if m else None


def run_curcuma(xyz_path, charge):
    cmd = [str(CURCUMA), "-sp", str(xyz_path), "-method", "gfnff",
           "-charge", str(charge), "-verbosity", "0", "-no_bmt"]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    return parse_curcuma_energy(proc.stdout), proc.stdout + proc.stderr


def run_xtb(xyz_path, charge, workdir):
    """Copy xyz + .CHRG into a clean workdir and run xtb --gfnff --sp there."""
    workdir.mkdir(parents=True, exist_ok=True)
    local_xyz = workdir / xyz_path.name
    shutil.copy(xyz_path, local_xyz)
    (workdir / ".CHRG").write_text(f"{charge:+d}\n")
    cmd = [str(XTB), str(local_xyz)] + XTB_ARGS
    env = dict(os.environ)
    env["XTBPATH"] = ""  # avoid picking up stray param files
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600,
                          cwd=str(workdir), env=env)
    return parse_xtb_energy(proc.stdout), proc.stdout + proc.stderr


def main():
    only = set(int(x) for x in sys.argv[1:]) if len(sys.argv) > 1 else None
    RUNDIR.mkdir(parents=True, exist_ok=True)
    refs = [float(x) for x in (TESTSET / "reference_s30l").read_text().split()]
    assert len(refs) == 30, f"expected 30 reference values, got {len(refs)}"

    rows = []
    for i in range(1, 31):
        if only and i not in only:
            continue
        syscomp = TESTSET / str(i)
        if not syscomp.is_dir():
            print(f"[{i:2d}] missing folder, skip", flush=True)
            continue
        energies = {}   # (engine, part) -> Eh or None
        logdir = RUNDIR / str(i)
        logdir.mkdir(parents=True, exist_ok=True)
        charges = {}
        for part in ("A", "B", "AB"):
            coord = syscomp / part / "coord"
            if not coord.exists():
                print(f"[{i:2d}] {part}: no coord, skip part", flush=True)
                continue
            atoms = read_coord(coord)
            if not atoms:
                print(f"[{i:2d}] {part}: parsed 0 atoms, skip", flush=True)
                continue
            q = read_charge(syscomp / part)
            charges[part] = q
            xyz = logdir / f"{part}.xyz"
            write_xyz(xyz, atoms, comment=f"S30L-{i}-{part} q={q}")
            # curcuma
            e_cur, out_cur = run_curcuma(xyz, q)
            energies[("cur", part)] = e_cur
            (logdir / f"{part}.cur.log").write_text(out_cur)
            # xtb
            xwdir = logdir / f"{part}_xtb"
            e_xtb, out_xtb = run_xtb(xyz, q, xwdir)
            energies[("xtb", part)] = e_xtb
            (logdir / f"{part}.xtb.log").write_text(out_xtb)
            # clean xtb scratch files
            for scratch in ("xtbrestart", "charges", "wbo", "gfnff_topo",
                            "gfnff_charges", "gfnff_adjacency"):
                sf = xwdir / scratch
                if sf.exists():
                    try:
                        sf.unlink()
                    except OSError:
                        pass
            print(f"[{i:2d}] {part} q={q:+d}  cur={e_cur}  xtb={e_xtb}", flush=True)

        def assoc(engine):
            eab = energies.get((engine, "AB"))
            ea = energies.get((engine, "A"))
            eb = energies.get((engine, "B"))
            if eab is None or ea is None or eb is None:
                return None
            return (eab - ea - eb) * AU2KCAL

        de_cur = assoc("cur")
        de_xtb = assoc("xtb")
        de_ref = refs[i - 1]
        d_cur_xtb = (de_cur - de_xtb) if (de_cur is not None and de_xtb is not None) else None
        d_cur_ref = (de_cur - de_ref) if de_cur is not None else None
        d_xtb_ref = (de_xtb - de_ref) if de_xtb is not None else None
        expected = [("cur", p) for p in ("A", "B", "AB")] + \
                  [("xtb", p) for p in ("A", "B", "AB")]
        status = "ok" if all(energies.get(k) is not None for k in expected) else "FAIL"
        rows.append({
            "i": i, "qAB": charges.get("AB", 0),
            "de_cur": de_cur, "de_xtb": de_xtb, "de_ref": de_ref,
            "d_cur_xtb": d_cur_xtb, "d_cur_ref": d_cur_ref, "d_xtb_ref": d_xtb_ref,
            "status": status,
        })

    # write CSV
    csv_path = RUNDIR / "results.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["i", "qAB", "de_cur", "de_xtb", "de_ref",
                    "d_cur_xtb", "d_cur_ref", "d_xtb_ref", "status"])
        for r in rows:
            w.writerow([r["i"], r["qAB"],
                        f"{r['de_cur']:.4f}" if r["de_cur"] is not None else "",
                        f"{r['de_xtb']:.4f}" if r["de_xtb"] is not None else "",
                        f"{r['de_ref']:.4f}",
                        f"{r['d_cur_xtb']:.4f}" if r["d_cur_xtb"] is not None else "",
                        f"{r['d_cur_ref']:.4f}" if r["d_cur_ref"] is not None else "",
                        f"{r['d_xtb_ref']:.4f}" if r["d_xtb_ref"] is not None else "",
                        r["status"]])

    # summary stats
    def stats(key):
        vals = [abs(r[key]) for r in rows if r[key] is not None]
        if not vals:
            return None
        n = len(vals)
        mad = sum(vals) / n
        mx = max(vals)
        rms = math.sqrt(sum(v * v for v in vals) / n)
        return n, mad, mx, rms

    md_path = RUNDIR / "results.md"
    with md_path.open("w") as f:
        f.write("# S30L: Curcuma GFN-FF vs xtb 6.6.1 GFN-FF\n\n")
        f.write("All association energies in kcal/mol. DE = E(AB)-E(A)-E(B).\n")
        f.write("d_cur_xtb = curcuma - xtb (primary); d_*_ref = engine - reference_s30l.\n\n")
        f.write("| i | qAB | DE_cur | DE_xtb | DE_ref | d_cur_xtb | d_cur_ref | d_xtb_ref | status |\n")
        f.write("|---|----|-------:|-------:|-------:|---------:|---------:|---------:|:------|\n")
        for r in rows:
            f.write("| {i} | {q:+d} | {c} | {x} | {ref:.2f} | {cx} | {cr} | {xr} | {s} |\n".format(
                i=r["i"], q=r["qAB"],
                c=f"{r['de_cur']:.2f}" if r["de_cur"] is not None else "-",
                x=f"{r['de_xtb']:.2f}" if r["de_xtb"] is not None else "-",
                ref=r["de_ref"],
                cx=f"{r['d_cur_xtb']:.2f}" if r["d_cur_xtb"] is not None else "-",
                cr=f"{r['d_cur_ref']:.2f}" if r["d_cur_ref"] is not None else "-",
                xr=f"{r['d_xtb_ref']:.2f}" if r["d_xtb_ref"] is not None else "-",
                s=r["status"]))
        f.write("\n## Summary (|deviation|, kcal/mol)\n\n")
        for label, key in [("curcuma vs xtb", "d_cur_xtb"),
                           ("curcuma vs ref", "d_cur_ref"),
                           ("xtb vs ref", "d_xtb_ref")]:
            s = stats(key)
            if s is None:
                f.write(f"- {label}: no data\n")
            else:
                n, mad, mx, rms = s
                f.write(f"- {label}: n={n}  MAD={mad:.3f}  max={mx:.3f}  RMSD={rms:.3f}\n")
        # flag outliers
        outliers = [r for r in rows if r["d_cur_xtb"] is not None and abs(r["d_cur_xtb"]) > 1.0]
        f.write(f"\n## Outliers (|d_cur_xtb| > 1.0 kcal/mol): {len(outliers)}\n\n")
        for r in outliers:
            f.write(f"- system {r['i']} (qAB={r['qAB']:+d}): d_cur_xtb={r['d_cur_xtb']:.3f} "
                    f"(DE_cur={r['de_cur']:.2f}, DE_xtb={r['de_xtb']:.2f})\n")

    print(f"\nWrote {csv_path}\nWrote {md_path}")
    # also print summary to stdout
    s = stats("d_cur_xtb")
    if s:
        n, mad, mx, rms = s
        print(f"\n=== curcuma vs xtb (n={n}): MAD={mad:.3f}  max={mx:.3f}  RMSD={rms:.3f} kcal/mol ===")


if __name__ == "__main__":
    main()