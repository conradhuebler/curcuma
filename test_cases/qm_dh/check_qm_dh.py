#!/usr/bin/env python3
"""ctest driver for the native double hybrids (WP12 stage 1), see CMakeLists.txt in this directory.

Gates (definitions as in check_qm_mgga.py; measured values in docs/DFT_ROADMAP/WP12_doppelhybride.md):
  kernel     scaled error |d| / max(|ref|, scale) <= 1e-10 for eps, v_rho, v_sigma, v_tau against libxc, where scale is
             the size of the cancelling terms (|eps_LDA(rho)| for eps and v_rho, rho |eps_LDA| / sigma for v_sigma,
             rho |eps_LDA| / tau for v_tau); FD of rho eps <= 1e-6 relative to max(|v|, scale) for rho >= 1e-8,
             1e-5 below (noise-limited stencil)
  energy     |E_KS - E_KS(PySCF)| <= 1e-8 Eh, |E_xc| <= 1e-7 Eh (first order in the density), |E_os|, |E_ss| <= 1e-9 Eh
  fc / canonical  the same PT2 gates for frozen core / exact-integral PT2
  d3         |E_D3(BJ) - ORCA 6.1| <= 1e-9 Eh (ORCA prints 12 digits)
  threads    1 vs 4 threads: total energy bitwise; update: UpdateMolecule A -> B = fresh B to 1e-10 Eh
  gradient / cli_opt  since stage 2 (Oct 2026): the gradient request and -opt must succeed (values in ../qm_dh_grad)
  cli_defaults  Oct 31 2026 defaults: DH frozen core (1 orbital for H2O) and VWN5-P86, explicit false wins, mp2 no fc
  d4         (WP12 stage 4) native D4 term alone vs dftd4 4.3.1 (reference_d4.json): |E|, |E_2body|, |E_ATM| <= 1e-12 Eh,
             charges and CN <= 1e-10, gradient <= 1e-11 Eh/Bohr; vs ORCA 6.1 (reference_orca_d4.json, 9 printed
             decimals) |E|, |g| <= 1e-9; analytic gradient vs 4-point FD of the D4 energy <= 1e-10 Eh/Bohr
  dhd4 / ksd4  method with D4 - method without = the D4 term (energy 1e-12 Eh, gradient 1e-12 Eh/Bohr); dhd4 with
             --fd: whole-method gradient vs 4-point FD of its energy <= 5e-10 Eh/Bohr (moving grid, VWN5-P86)
Claude Generated (Oct 2026, WP12 stage 1).
"""
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.dirname(HERE)
CODATA_LDA = (3.0 / math.pi) ** (1.0 / 3.0)


def run(cmd, timeout=900):
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return p.returncode, p.stdout, p.stderr


def driver(a, args):
    rc, out, err = run([a.exe] + args)
    if rc != 0:
        print(err)
        sys.exit(f"driver failed ({rc}): {' '.join(args)}")
    # CurcumaLogger::error prints even at verbosity 0 (e.g. the gradient error): the JSON is the last "{" block
    # (ANSI colour reset in front of it), so start at the first "{" that ends a line
    start = out.find("{\n")
    return json.loads(out[start:])


def split_case(case):
    """'H2O_dsd-pbep86@vwn5_sg1' -> molecule, method, p86 variant, grid. Without '@' the P86 variant is pz81: the
    PySCF/libxc references (GGA_C_P86 = PZ81) need it, so it is passed explicitly (curcuma's default is vwn5 since
    Oct 31 2026); frozen core is passed explicitly as well (--fc or not; curcuma's DH default is fc since then)."""
    parts = case.split("_")
    mol, grid = parts[0], parts[-1] if len(parts) > 2 else "sg1"
    fn = "_".join(parts[1:-1]) if len(parts) > 2 else parts[1]
    p86 = "pz81"
    if "@" in fn:
        fn, p86 = fn.split("@")
    return mol, fn, p86, grid


def lda_scale(rho):
    return 0.75 * CODATA_LDA * rho ** (1.0 / 3.0)


def check_kernel(a, ref):
    out = driver(a, ["kernel", a.ref])
    worst_all = 0.0
    fails = 0
    for block, keys in (("gga", ("p86", "p86vwn")), ("mgga", ("b95",))):
        g = ref["kernel"][block]
        rho, sigma = g["rho"], g["sigma"]
        tau = g.get("tau")
        for key in keys:
            r, c = g[key], out[key]
            comps = ["eps", "vrho", "vsigma"] + (["vtau"] if block == "mgga" else [])
            for comp in comps:
                worst = wfd = 0.0
                for i in range(len(rho)):
                    if block == "mgga" and g["kind"][i] != "admissible" and comp != "eps":
                        pass  # the clamped points are compared as well (both lower sigma to 8 rho tau)
                    el = lda_scale(rho[i])
                    s = max(sigma[i], g[key]["sigma_floor"])
                    scale = el if comp in ("eps", "vrho") else (rho[i] * el / s if comp == "vsigma" else rho[i] * el / max(tau[i], 1e-20))
                    d = abs(c[comp][i] - r[comp][i]) / max(abs(r[comp][i]), scale)
                    worst = max(worst, d)
                    fd = c.get(comp + "_fd")
                    if fd is not None and fd[i] is not None:
                        dfd = abs(fd[i] - c[comp][i]) / max(abs(c[comp][i]), scale)
                        wfd = max(wfd, dfd)
                        # rho < 1e-8: the 4-point stencil is noise-limited (cancellation of the PW92 parts of B95;
                        # the analytic v agrees with libxc to 5e-15 there), gate 1e-5 instead of 1e-6
                        if dfd > (1e-6 if rho[i] >= 1e-8 else 1e-5):
                            fails += 1
                            print(f"FD {key} {comp} point {i} rho {rho[i]:.1e} sigma {sigma[i]:.2e}: {dfd:.2e}")
                print(f"{key:7s} {comp:7s} max scaled |d| vs libxc {worst:.2e}, vs FD {wfd:.2e}")
                worst_all = max(worst_all, worst)
    if worst_all > 1e-10 or fails:
        sys.exit(f"kernel: max {worst_all:.2e} (gate 1e-10), {fails} FD failures")
    print("OK")


def energy_args(mol, fn, p86, grid, extra=()):
    xyz = os.path.join(TC, REF_XYZ[mol])
    args = ["energy", xyz, "--method", fn, "--grid", grid, "--p86", p86, "--dispersion", "none"] + list(extra)
    return args


def compare(tag, got, want, tol):
    d = got - want
    print(f"  {tag:5s} curcuma {got:.12f} ref {want:.12f} d {d:+.2e} (gate {tol:.0e})")
    return abs(d) <= tol


def check_energy(a, ref, case, mode):
    mol, fn, p86, grid = split_case(case)
    r = ref["energies"][case]
    extra = []
    key = "ri_full"
    if mode == "fc":
        extra, key = ["--fc"], "ri_fc"
    elif mode == "canonical":
        extra, key = ["--pt2", "canonical"], "canonical_full"
    out = driver(a, energy_args(mol, fn, p86, grid, extra))
    if out.get("error") or "exception" in out:
        sys.exit(f"error: {out}")
    ok = compare("E_KS", out["ks"], r["e_ks"], 1e-8)
    ok &= compare("E_xc", out["xc"], r["exc"], 1e-7)
    ok &= compare("E_os", out["os"], r[key]["e_os"], 1e-9)
    ok &= compare("E_ss", out["ss"], r[key]["e_ss"], 1e-9)
    if mode == "fc" and out["n_frozen"] != r["n_frozen_core"]:
        ok = False
    print(f"  iterations {out['scf_iterations']}, n_aux {out['n_aux']}, frozen {out['n_frozen']}")
    if not ok:
        sys.exit("FAIL")
    print("OK")


def check_d3(a, orca, case):
    mol, fn = case.split("_", 1)
    want = orca["d3"][case]
    out = driver(a, ["energy", os.path.join(TC, REF_XYZ[mol]), "--method", fn, "--dispersion", "auto"])
    if out.get("error"):
        sys.exit(f"error: {out}")
    if not compare("E_D3", out["disp"], want, 1e-9):
        sys.exit("FAIL")
    print("OK")


def check_simple(a, case, chk):
    mol, fn = case.split("_", 1)
    xyz = os.path.join(TC, REF_XYZ[mol])
    if chk == "gradient":
        # WP12 stage 2 (Oct 2026): the gradient exists now; a request must succeed (values: test_cases/qm_dh_grad)
        out = driver(a, ["energy", xyz, "--method", fn, "--gradient"])
        print(out.get("message"))
        if not out.get("has_gradient") or out.get("error"):
            sys.exit("a gradient request must succeed and hasGradient() be true (WP12 stage 2)")
    elif chk == "threads":
        e1 = driver(a, ["energy", xyz, "--method", fn, "--threads", "1"])["energy"]
        e4 = driver(a, ["energy", xyz, "--method", fn, "--threads", "4"])["energy"]
        print(f"  1 thread {e1!r}, 4 threads {e4!r}")
        if e1 != e4:
            sys.exit("not bitwise")
    elif chk == "update":
        xyzB = os.path.join(TC, "qm_grad", "H2O_distorted.xyz")
        if not os.path.exists(xyzB):
            xyzB = os.path.join(TC, "qm_1e", "H2O.xyz")
        out = driver(a, ["update", xyz, xyzB, "--method", fn])
        d = abs(out["moving"] - out["fresh"])
        print(f"  moving {out['moving']:.12f} fresh {out['fresh']:.12f} d {d:.1e}")
        if d > 1e-10:
            sys.exit("FAIL")
    print("OK")


def curcuma_run(a, args, workdir):
    rc, out, err = run([a.curcuma] + args, timeout=900)
    return rc, out + err


def check_cli(a, chk):
    xyz = os.path.join(TC, "qm_1e", "H2O.xyz")
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(xyz, os.path.join(tmp, "H2O.xyz"))
        cwd = os.getcwd()
        os.chdir(tmp)
        try:
            if chk == "cli_sp":
                rc, txt = curcuma_run(a, ["-sp", "H2O.xyz", "-method", "b2plyp", "-no_bmt"], tmp)
                m = re.search(r"E\(b2plyp\)\s+=\s+(-?\d+\.\d+)", txt)
                print(txt[-1500:])
                if rc != 0 or not m:
                    sys.exit(f"-sp failed (exit {rc})")
                print("  E(b2plyp) =", m.group(1))
            elif chk == "cli_opt":
                # WP12 stage 2 (Oct 2026): -opt runs on the analytic gradient (geometry checks: test_cases/qm_dh_grad)
                rc, txt = curcuma_run(a, ["-opt", "H2O.xyz", "-method", "b2plyp", "-no_bmt"], tmp)
                print(txt[-1500:])
                if rc != 0 or not os.path.exists(os.path.join(tmp, "H2O.opt.xyz")) or "no analytic gradient" in txt:
                    sys.exit(f"-opt with b2plyp must run on the analytic gradient (exit {rc})")
            elif chk == "cli_defaults":
                # Oct 31 2026 defaults: double hybrids frozen core + VWN5-P86, explicit false still wins, mp2 keeps no fc
                def frozen(args):
                    rc, txt = curcuma_run(a, ["-sp", "H2O.xyz", "-no_bmt", "-verbosity", "1"] + args, tmp)
                    txt = re.sub(r"\x1b\[[0-9;]*m", "", txt)
                    m = re.search(r"(\d+) frozen", txt)
                    return rc, (int(m.group(1)) if m else None), txt
                rc1, n1, t1 = frozen(["-method", "dsd-pbep86"])
                rc2, n2, _ = frozen(["-method", "dsd-pbep86", "-qm.frozen_core", "false"])
                rc3, n3, _ = frozen(["-method", "mp2"])
                rc4, n4, _ = frozen(["-method", "mp2", "-qm.frozen_core", "true"])
                vwn = "P86 on VWN5" in t1
                print(f"  dsd-pbep86 default: {n1} frozen, VWN5 {vwn}; -qm.frozen_core false: {n2}; mp2 default: {n3}; "
                      f"mp2 -qm.frozen_core true: {n4}")
                if (rc1, rc2, rc3, rc4) != (0, 0, 0, 0) or (n1, n2, n3, n4) != (1, 0, 0, 1) or not vwn:
                    print(t1[-1500:])
                    sys.exit("defaults: expected 1/0/0/1 frozen orbitals and P86 on VWN5")
            elif chk == "cli_d4":
                # WP12 stage 4 (Oct 2026): revdsd-pbep86-d4 runs with native D4; a DH without D4 parameters is an error
                rc, txt = curcuma_run(a, ["-sp", "H2O.xyz", "-method", "revdsd-pbep86-d4", "-no_bmt"], tmp)
                print(txt[-1500:])
                m = re.search(r"E\(disp D4\)\s+=\s+(-?\d+\.\d+)", txt)
                if rc != 0 or not m or not re.search(r"E\(revdsd-pbep86-d4\)\s+=", txt):
                    sys.exit(f"revdsd-pbep86-d4 must run and print E(disp D4) (exit {rc})")
                rc2, txt2 = curcuma_run(a, ["-sp", "H2O.xyz", "-method", "dsd-blyp", "-qm.dispersion", "d4", "-no_bmt"],
                                        tmp)
                print(txt2[-800:])
                if "has no D4 parameters" not in txt2:
                    sys.exit(f"-qm.dispersion d4 on dsd-blyp must be an error naming the D4 methods (exit {rc2})")
            elif chk == "cli_ks_d4":
                # WP12 stage 4: pbe-d4 = pbe + D4 on the command line, and -qm.dispersion d4 on pbe is the same
                vals = {}
                for tag, args in (("pbe", ["-method", "pbe"]), ("pbe-d4", ["-method", "pbe-d4"]),
                                  ("flag", ["-method", "pbe", "-qm.dispersion", "d4"])):
                    rc, txt = curcuma_run(a, ["-sp", "H2O.xyz", "-no_bmt"] + args, tmp)
                    m = re.findall(r"Single Point Energy = (-?\d+\.\d+) Eh", txt)
                    d = re.findall(r"E\(disp D4\)\s+=\s+(-?\d+\.\d+)", txt)
                    if rc != 0 or not m:
                        print(txt[-1500:])
                        sys.exit(f"-sp {args} failed (exit {rc})")
                    vals[tag] = (float(m[-1]), float(d[-1]) if d else 0.0)
                print(f"  pbe {vals['pbe']}, pbe-d4 {vals['pbe-d4']}, pbe -qm.dispersion d4 {vals['flag']}")
                e0, (e4, d4), (ef, df) = vals["pbe"][0], vals["pbe-d4"], vals["flag"]
                if d4 == 0.0 or abs(e4 - e0 - d4) > 2e-8 or abs(ef - e4) > 1e-12:
                    sys.exit("pbe-d4 must equal pbe + E(disp D4) (8 printed decimals) and pbe -qm.dispersion d4")
            elif chk == "cli_opt_d4":
                # WP12 stage 4: -opt with revdsd-pbep86-d4 (analytic gradient incl. D4) from a distorted H2O
                shutil.copy(os.path.join(TC, "qm_grad", "H2O_distorted.xyz"), os.path.join(tmp, "H2Od.xyz"))
                rc, txt = curcuma_run(a, ["-opt", "H2Od.xyz", "-method", "revdsd-pbep86-d4", "-no_bmt",
                                          "-opt.convergence_preset", "tight"], tmp)
                print(txt[-1500:])
                if rc != 0 or not os.path.exists(os.path.join(tmp, "H2Od.opt.xyz")):
                    sys.exit(f"-opt with revdsd-pbep86-d4 failed (exit {rc})")
            elif chk == "cli_unknown_option":
                rc, txt = curcuma_run(a, ["-sp", "H2O.xyz", "-method", "b2plyp", "-qm.dispersion", "d3bj", "-no_bmt"],
                                      tmp)
                print(txt[-1200:])
                if "has no D3(BJ) parameters" not in txt:
                    sys.exit(f"-qm.dispersion d3bj on b2plyp must be an error naming b2plyp-d3bj (exit {rc})")
        finally:
            os.chdir(cwd)
    print("OK")


def check_d4(a, case):
    """WP12 stage 4: the D4 term alone against dftd4 4.3.1 and ORCA 6.1, analytic gradient against FD."""
    ref = json.load(open(os.path.join(HERE, "reference_d4.json")))["cases"][case]
    orca = json.load(open(os.path.join(HERE, "reference_orca_d4.json")))["disp"].get(case)
    out = driver(a, ["d4", os.path.join(TC, ref["xyz"]), "--method", ref["functional"], "--fd"])
    if out.get("error"):
        sys.exit(f"error: {out}")
    md = lambda x, y: max(abs(p - q) for p, q in zip(x, y))
    ok = compare("E_D4", out["energy"], ref["energy"], 1e-12)
    ok &= compare("E2", out["two_body"], ref["two_body"], 1e-12)
    ok &= compare("E3", out["three_body"], ref["three_body"], 1e-12)
    rows = [("q", md(out["charges"], ref["charges"]), 1e-10), ("CN", md(out["cn"], ref["cn"]), 1e-10),
            ("g-dftd4", md(out["gradient"], ref["gradient"]), 1e-11),
            ("g-FD", md(out["gradient"], out["fd_gradient"]), 1e-10)]
    if orca:
        ok &= compare("E_ORCA", out["energy"], orca["energy"], 1e-9)
        rows.append(("g-ORCA", md(out["gradient"], orca["gradient"]), 1e-9))
    for tag, d, tol in rows:
        print(f"  {tag:8s} max |d| {d:.2e} (gate {tol:.0e})")
        ok &= d <= tol
    if not ok:
        sys.exit("FAIL")
    print("OK")


def check_with_without(a, case, mode):
    """WP12 stage 4: E(method with D4) - E(without) and the gradients = the D4 term alone; dhd4 also vs FD."""
    mol, fn = case.split("_", 1)
    xyz = os.path.join(TC, REF_XYZ.get(mol, {"H2CO": "qm_grad/formaldehyde_hf3c_min_reference.xyz"}.get(mol, "")))
    args = [mode, xyz, "--method", fn]
    if mode == "dhd4":
        args += ["--p86", "vwn5", "--fc"]
        if mol == "H2O":
            args += ["--fd"]
    out = driver(a, args)
    w, n, d4 = out["with"], out["none"], out["d4"]
    if w["error"] or n["error"]:
        sys.exit(f"error: {w['message']} / {n['message']}")
    md = lambda x, y: max(abs(p - q) for p, q in zip(x, y))
    de = w["energy"] - n["energy"] - d4["energy"]
    dg = md([p - q for p, q in zip(w["gradient"], n["gradient"])], d4["gradient"])
    print(f"  {w['name']}: E {w['energy']:.12f}, E(none) {n['energy']:.12f}, E_D4 {d4['energy']:.12f}; "
          f"E - E(none) - E_D4 {de:+.1e}; max |g - g(none) - g_D4| {dg:.1e}")
    ok = abs(de) <= 1e-12 and dg <= 1e-12
    if "fd_gradient" in out:
        f = md(w["gradient"], out["fd_gradient"])
        print(f"  max |g - FD(E)| {f:.2e} Eh/Bohr (gate 5e-10)")
        ok &= f <= 5e-10
    if not ok:
        sys.exit("FAIL")
    print("OK")


def check_gpu(a, case):
    mol, fn = case.split("_", 1)
    xyz = os.path.join(TC, REF_XYZ[mol])
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(xyz, os.path.join(tmp, "in.xyz"))
        cwd = os.getcwd()
        os.chdir(tmp)
        try:
            base = ["-sp", "in.xyz", "-method", fn, "-no_bmt", "-qm.scf_threshold", "1e-10"]
            rc0, cpu = curcuma_run(a, base, tmp)
            rc1, gpu = curcuma_run(a, base + ["-gpu", "cuda"], tmp)
            if "no GPU plugin" in gpu or "plugin" in gpu and "not available" in gpu:
                print("no CUDA plugin/device: skipped")
                sys.exit(77)
            rx = re.compile(r"E\(" + re.escape(fn) + r"\)\s+=\s+(-?\d+\.\d+)")
            mc, mg = rx.search(cpu), rx.search(gpu)
            if not mc or not mg:
                print(gpu[-2000:])
                sys.exit("no energy line")
            d = abs(float(mc.group(1)) - float(mg.group(1)))
            print(f"  CPU {mc.group(1)} GPU {mg.group(1)} d {d:.1e}")
            # WP12 stage 3 (Oct 31 2026): P86/B95 run on the device, no qm_xc fallback; -gpu_strict true runs through
            if "qm_xc" in gpu or "GPU fallbacks during this run" in gpu:
                print(gpu[-2000:])
                sys.exit("expected the device XC path (no counted qm_xc fallback)")
            rc2, strict = curcuma_run(a, base + ["-gpu", "cuda", "-gpu_strict", "true"], tmp)
            ms = rx.search(strict)
            print(f"  -gpu_strict true exit {rc2}")
            if d > 1e-8 or rc2 != 0 or not ms or abs(float(ms.group(1)) - float(mc.group(1))) > 1e-8:
                print(strict[-2000:])
                sys.exit("FAIL")
        finally:
            os.chdir(cwd)
    print("OK")


REF_XYZ = {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--orca", required=True)
    ap.add_argument("--curcuma", required=True)
    ap.add_argument("--check", required=True)
    ap.add_argument("--case", default="")
    a = ap.parse_args()
    ref = json.load(open(a.ref))
    for e in ref["energies"].values():
        REF_XYZ[e["molecule"]] = e["xyz"]
    if a.check == "kernel":
        check_kernel(a, ref)
    elif a.check in ("energy", "fc", "canonical"):
        check_energy(a, ref, a.case, a.check)
    elif a.check == "d3":
        check_d3(a, json.load(open(a.orca)), a.case)
    elif a.check in ("gradient", "threads", "update"):
        check_simple(a, a.case, a.check)
    elif a.check == "d4":
        check_d4(a, a.case)
    elif a.check in ("dhd4", "ksd4"):
        check_with_without(a, a.case, a.check)
    elif a.check.startswith("cli_"):
        check_cli(a, a.check)
    elif a.check == "gpu":
        check_gpu(a, a.case)
    else:
        sys.exit("unknown check " + a.check)


if __name__ == "__main__":
    main()
