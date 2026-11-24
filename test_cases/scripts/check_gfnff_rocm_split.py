#!/usr/bin/env python3
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026, docs/ROCM_PARITY_A6.md)
"""GFN-FF multi-GPU split (EEQ, Coulomb tiles, D4 dispersion) on ROCm, checked on ONE device.

The ROCm port of the projected-PCG EEQ split (gpu_split_devices, opt-in -gfnff.rocm_multi_gpu_split true)
cannot be run over several GPUs on a machine with one GPU. The test hook CURCUMA_ROCM_SPLIT_SAME_DEVICE=1 keeps
helper entries that name the calculation's own device, so "-gpu_split_devices 0,0" runs the split code with two
participants (two column blocks, own streams, events and peer copies) on device 0, and "0,0,0" with three.

Checks (short NVE MD of polymer_2x.xyz, 7320 atoms, 1502 fragments, projected PCG):
  1. split (2 participants) == single device to the CUDA test limits (charges 1e-9 e, energy terms 1e-9 Eh,
     gradient norm 1e-8): the split changes only the summation order of the matvec;
  2. split with 2 and with 3 participants are bitwise identical (the column dot products of
     k_split_coldot_*_hip do not depend on the column partition);
  3. without the opt-in an explicit device list is refused: warning, counted GPU fallback, -gpu_strict exits 3.
Coulomb-tile and D4 split alone (Nov 2026, docs/ROCM_PARITY_SPL.md; EEQ split kept off with a huge
gpu_split_eeq_min_atoms), 2 and 4 participants:
  4. with the fixed-point accumulation (-gfnff.rocm_fixed_point, default) bitwise equal to the single device: every
     tile row block and gather row is evaluated as in the unsplit launch and integer sums do not depend on the order;
  5. with FP64 atomics (-gfnff.rocm_fixed_point false) within 1e-8 Eh per energy term and 1e-8 in the gradient norm;
  6. -sp -dump_gradient of polymer.xyz (1410 atoms), 3 participants: gradient component-wise within 1e-8 Eh/A
     (bitwise with fixed point), and the status lines name the Coulomb tiles and the dispersion on devices [0,0,0].

Not covered: any run with more than one physical device (peer access, P2P copies between cards).

Exit 0 = pass, 1 = fail, 77 = skipped (no ROCm plugin or no HIP device).
Usage: check_gfnff_rocm_split.py <curcuma-binary> <polymer_2x.xyz> [<polymer.xyz>]
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

TOL_Q = 1e-9      # e (as check_gfnff_gpu_split.py)
TOL_E = 1e-9      # Eh, per energy term
TOL_G = 1e-8      # gradient norm


def flat(x):
    if isinstance(x, list):
        o = []
        for y in x:
            o.extend(flat(y))
        return o
    return [x]


def run(binary, xyz, workdir, extra, env_extra=None, expect_rc=0):
    os.makedirs(workdir, exist_ok=True)
    shutil.copy(xyz, workdir)
    base = os.path.splitext(os.path.basename(xyz))[0]
    cmd = [binary, "-md", os.path.basename(xyz), "-method", "gfnff", "-gpu", "rocm",
           "-gpu_strict", "true", "-no_bmt", "-verbosity", "2", "-maxtime", "4", "-dt", "1",
           "-thermostat", "none", "-seed", "7", "-md_diagnostics", "true", "-dump", "1",
           "-gpu_split_min_atoms", "0"] + extra
    env = dict(os.environ)
    env.update(env_extra or {})
    r = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=900, env=env)
    if r.returncode != expect_rc:
        print(f"run {' '.join(extra)}: rc={r.returncode}, expected {expect_rc}\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
        return None, r
    if expect_rc != 0:
        return [], r
    with open(os.path.join(workdir, base + ".diag.jsonl")) as f:
        return [json.loads(l) for l in f], r


def run_sp(binary, xyz, workdir, extra, env_extra=None):
    """-sp -dump_gradient: (energy, gradient rows, stdout) or None."""
    os.makedirs(workdir, exist_ok=True)
    shutil.copy(xyz, workdir)
    cmd = [binary, "-sp", os.path.basename(xyz), "-method", "gfnff", "-gpu", "rocm", "-gpu_strict", "true",
           "-no_bmt", "-verbosity", "2", "-dump_gradient", "grad.txt", "-gpu_split_min_atoms", "0"] + extra
    env = dict(os.environ)
    env.update(env_extra or {})
    r = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=900, env=env)
    if r.returncode != 0:
        print(f"sp {' '.join(extra)}: rc={r.returncode}\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
        return None
    energy, rows = None, []
    with open(os.path.join(workdir, "grad.txt")) as f:
        for line in f:
            if line.startswith("# energy"):
                energy = float(line.split()[2])
            elif not line.startswith("#") and line.strip():
                rows.append([float(x) for x in line.split()])
    return energy, rows, r.stdout


def identical(a_list, b_list):
    return len(a_list) == len(b_list) and all(a == b for a, b in zip(a_list, b_list))


def worst(a_list, b_list):
    w = {"q": 0.0, "E": 0.0, "g": 0.0}
    for a, b in zip(a_list, b_list):
        w["q"] = max(w["q"], max(abs(x - y) for x, y in zip(flat(a["charges"]), flat(b["charges"]))))
        w["E"] = max(w["E"], max(abs(a["energy"][k] - b["energy"][k]) for k in a["energy"]))
        w["g"] = max(w["g"], max(abs(x - y) for x, y in zip(flat(a["gradient_norm"]), flat(b["gradient_norm"]))))
    return w


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    binary, xyz = sys.argv[1], sys.argv[2]
    plugin = os.path.join(os.path.dirname(os.path.abspath(binary)), "libcurcuma_rocm.so")
    if not os.path.exists(plugin):
        print("SKIP: no ROCm plugin next to the binary")
        return 77
    tmp = tempfile.mkdtemp(prefix="gfnff_rocm_split_")
    hook = {"CURCUMA_ROCM_SPLIT_SAME_DEVICE": "1"}
    optin = ["-gfnff.rocm_multi_gpu_split", "true"]
    try:
        single, r = run(binary, xyz, os.path.join(tmp, "single"), ["-gpu_split_devices", "none"])
        if single is None:
            if "no HIP device" in (r.stdout + r.stderr) or "No HIP" in (r.stdout + r.stderr):
                print("SKIP: no HIP device")
                return 77
            return 1
        two, r2 = run(binary, xyz, os.path.join(tmp, "two"), ["-gpu_split_devices", "0,0"] + optin, hook)
        three, r3 = run(binary, xyz, os.path.join(tmp, "three"), ["-gpu_split_devices", "0,0,0"] + optin, hook)
        if two is None or three is None or not single or len(two) != len(single) or len(three) != len(single):
            print("FAIL: missing diagnostics")
            return 1
        # The status line must say that the split ran (verbosity 2), not that it fell back.
        for name, rr in (("two", r2), ("three", r3)):
            if "EEQ matrix + PCG matvec on devices" not in rr.stdout:
                print(f"FAIL: {name}: no split status line in the output")
                return 1
        ok = True
        w = worst(two, single)
        print(f"split(2) vs single, {len(two)} MD records: max|dq| {w['q']:.2e} e, max|dE term| {w['E']:.2e} Eh, "
              f"max|dgrad| {w['g']:.2e}")
        ok &= w["q"] <= TOL_Q and w["E"] <= TOL_E and w["g"] <= TOL_G
        w23 = worst(two, three)
        print(f"split(2) vs split(3): max|dq| {w23['q']:.2e} e, max|dE term| {w23['E']:.2e} Eh, "
              f"max|dgrad| {w23['g']:.2e} (expected 0)")
        ok &= w23["q"] == 0.0 and w23["E"] == 0.0 and w23["g"] == 0.0
        # Refusal without the opt-in: an explicit list is warned about and counted (-gpu_strict -> exit 3).
        _, rr = run(binary, xyz, os.path.join(tmp, "refused"), ["-gpu_split_devices", "0,1"], None, expect_rc=3)
        if rr.returncode != 3 or "rocm_multi_gpu_split" not in (rr.stdout + rr.stderr):
            print(f"FAIL: refusal without opt-in: rc={rr.returncode}")
            ok = False
        else:
            print("refusal without -gfnff.rocm_multi_gpu_split: exit 3 with message")
        # Coulomb-tile + D4 split alone (EEQ split off), fixed point (default) and FP64 atomics.
        no_eeq = ["-gfnff.gpu_split_eeq_min_atoms", "1000000000"]
        fixed = True   # -gfnff.rocm_fixed_point defaults to true (gfnff.h, RocmGfnffOptions)
        for fp_name, fp_args in (("fixed point", []),
                                 ("FP64 atomics", ["-gfnff.rocm_fixed_point", "false"])):
            ref, _ = run(binary, xyz, os.path.join(tmp, "cd_single" + str(len(fp_args))),
                         ["-gpu_split_devices", "none"] + fp_args)
            for parts in ("0,0", "0,0,0,0"):
                got, rg = run(binary, xyz, os.path.join(tmp, "cd" + parts.replace(",", "") + str(len(fp_args))),
                              ["-gpu_split_devices", parts] + optin + no_eeq + fp_args, hook)
                if ref is None or got is None or len(got) != len(ref):
                    print(f"FAIL: Coulomb/D4 split {parts} ({fp_name}): missing diagnostics")
                    ok = False
                    continue
                for tag in ("Coulomb tiles on devices", "dispersion on devices"):
                    if tag not in rg.stdout:
                        print(f"FAIL: Coulomb/D4 split {parts} ({fp_name}): no status line '{tag}'")
                        ok = False
                w = worst(got, ref)
                same = identical(got, ref)
                print(f"Coulomb/D4 split {parts} vs single ({fp_name}): max|dq| {w['q']:.2e} e, "
                      f"max|dE term| {w['E']:.2e} Eh, max|dgrad| {w['g']:.2e}, records identical: {same}")
                if not fp_args:
                    ok &= same
                else:
                    ok &= w["q"] <= TOL_Q and w["E"] <= 1e-8 and w["g"] <= 1e-8
        # Same-device helper on its own stream and event (hook value "streams", the path of a real helper device),
        # 2 participants: bitwise equal to the sequential partitions with fixed point.
        seq, _ = run(binary, xyz, os.path.join(tmp, "cd_seq2"), ["-gpu_split_devices", "0,0"] + optin + no_eeq, hook)
        strm, rs = run(binary, xyz, os.path.join(tmp, "cd_strm2"), ["-gpu_split_devices", "0,0"] + optin + no_eeq,
                       {"CURCUMA_ROCM_SPLIT_SAME_DEVICE": "streams"})
        if seq is None or strm is None or "Coulomb tiles on devices [0,0]" not in rs.stdout:
            print("FAIL: Coulomb/D4 split with helper streams")
            ok = False
        else:
            w = worst(strm, seq)
            same = identical(strm, seq)
            print(f"Coulomb/D4 split 0,0 helper stream vs sequential: max|dE term| {w['E']:.2e} Eh, "
                  f"max|dgrad| {w['g']:.2e}, records identical: {same}")
            ok &= same if fixed else (w["E"] <= 1e-8 and w["g"] <= 1e-8)
        # Component-wise gradient of one single point (polymer.xyz), 3 participants.
        if len(sys.argv) > 3:
            xyz1 = sys.argv[3]
            a = run_sp(binary, xyz1, os.path.join(tmp, "sp_single"), ["-gpu_split_devices", "none"])
            b = run_sp(binary, xyz1, os.path.join(tmp, "sp_three"), ["-gpu_split_devices", "0,0,0"] + optin + no_eeq,
                       hook)
            if a is None or b is None or len(a[1]) != len(b[1]):
                print("FAIL: -sp runs")
                ok = False
            else:
                dg = max(abs(x - y) for ra, rb in zip(a[1], b[1]) for x, y in zip(ra, rb))
                de = abs(a[0] - b[0])
                st = "Coulomb tiles on devices [0,0,0]" in b[2] and "dispersion on devices [0,0,0]" in b[2]
                print(f"-sp polymer split(3) vs single: |dE| {de:.2e} Eh, max|dG| {dg:.2e} Eh/A, status lines: {st}")
                ok &= st and dg <= 1e-8 and de <= 1e-8
                if fixed:
                    ok &= dg == 0.0 and de == 0.0
        print("PASS" if ok else "FAIL")
        return 0 if ok else 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
