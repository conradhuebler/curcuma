#!/usr/bin/env python3
"""Claude Generated (Oct 2026): -gfnff.charge_groups and -gfnff.cn_zero_guard on the CUDA GFN-FF path (no CPU fallback).

Uses only the curcuma binary (the same binary runs the CPU reference without -gpu); every structure is built at run time.
Every GPU run uses -gpu_strict true, so a GPU fallback (also the CPU method fallback) fails the test.
  1. cn_zero_guard: Cl- ... H2O (H-Cl 3.2 A, Cl first so frag_charge_model reference gives it -1, CN(Cl) = 0 exactly).
     GPU gradient against the CPU gradient within 1e-10 Eh/A with the switch on and off; the switch changes the GPU
     gradient by > 1e-4 Eh/A; with it the GPU gradient of the Cl and of one H agrees with a central finite difference of
     the GPU energy within 2e-5 Eh/A (6 components). The single point reuses the CPU topology charges, so this part
     checks the device gradient kernels, not the device EEQ.
  2. charge groups: 5-fs MD (10 steps, CSVR, same seed) GPU against CPU - here the device EEQ solves every step after
     the first. (AB)2 chain + polymerbuild group file (charge 2, frag_charge_model reference), the same chain with
     charge_groups auto, NH4+ Cl- with charge_groups ionpairs (charge 0), the chain + the group file + cn_zero_guard.
     Per-frame charges within 1e-9 e, every group sum within 1e-8 of its charge on the GPU, mean Epot within 1e-8 Eh,
     the device EEQ path reported ("[GPU CUDA").
Exit 0 = pass, 1 = fail, 77 = skipped (no CUDA device or no CUDA plugin next to the binary).
Usage: check_gpu_charge_groups.py CURCUMA_BINARY SOURCE_DIR
"""
import json, os, re, subprocess, sys, tempfile
import numpy as np

curcuma, src = sys.argv[1], sys.argv[2]
failed = 0
E_RE = re.compile(r"Single Point Energy = (-?\d+\.\d+)")
GPU = ["-gpu", "cuda", "-gpu_strict", "true"]
REF = ["-gfnff.frag_charge_model", "reference"]


def check(ok, what):
    global failed
    print(("  PASS  " if ok else "  FAIL  ") + what)
    failed += 0 if ok else 1


def visible_gpus():
    try:
        out = subprocess.run(["nvidia-smi", "-L"], capture_output=True, text=True, timeout=30)
        return sum(1 for l in out.stdout.splitlines() if l.startswith("GPU "))
    except Exception:
        return 0


if visible_gpus() < 1:
    print("SKIP: no CUDA device")
    sys.exit(77)
if not os.path.exists(os.path.join(os.path.dirname(os.path.abspath(curcuma)), "libcurcuma_cuda.so")):
    print("SKIP: no CUDA plugin next to the binary")
    sys.exit(77)

work = tempfile.mkdtemp(prefix="curcuma_gpu_charge_groups_")
os.chdir(work)
count = [0]


def write(path, el, X):
    with open(path, "w") as f:
        f.write(f"{len(el)}\nx\n" + "".join("%s %.10f %.10f %.10f\n" % (e, *x) for e, x in zip(el, X)))


def read_xyz(path):
    L = open(path).read().split("\n")
    m = int(L[0])
    return [l.split()[0] for l in L[2:2 + m]], np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + m]])


def lib(id_):
    return subprocess.run([sys.executable, src + "/scripts/structlib.py", "path", id_], capture_output=True,
                          text=True).stdout.strip()


def sp(el, X, charge, extra=(), grad=False):
    """single point in a fresh directory; returns (energy, gradient Eh/A or None, return code)"""
    count[0] += 1
    d = os.path.join(work, f"sp{count[0]}")
    os.mkdir(d)
    write(os.path.join(d, "in.xyz"), el, X)
    cmd = [curcuma, "-sp", "in.xyz", "-method", "gfnff", "-charge", str(charge), "-verbosity", "0", "-no_bmt",
           "-threads", "1", "-gfnff.cache_topology", "false"] + list(extra) + (["-dump_gradient", "g.txt"] if grad else [])
    p = subprocess.run(cmd, cwd=d, capture_output=True, text=True)
    m = E_RE.findall(p.stdout + p.stderr)
    G = np.loadtxt(os.path.join(d, "g.txt"), comments="#").reshape(len(el), 3) if grad and p.returncode == 0 else None
    return (float(m[-1]) if m else None), G, p.returncode


def md(xyz, charge, extra):
    """10 MD steps with per-step diagnostics; returns (rc, charges per frame, mean Epot, log text)"""
    count[0] += 1
    d = os.path.join(work, f"md{count[0]}")
    os.mkdir(d)
    subprocess.run(["/bin/cp", "-f", xyz, os.path.join(d, "in.xyz")], check=True)
    cmd = [curcuma, "-md", "in.xyz", "-method", "gfnff", "-charge", str(charge), "-no_bmt", "-verbosity", "1",
           "-threads", "2", "-gfnff.cache_topology", "false", "-md.max_time", "5", "-md.time_step", "0.5",
           "-md.thermostat", "csvr", "-md.temperature", "303", "-md.seed", "303", "-md.dump_frequency", "1",
           "-md.md_diagnostics", "true", "-md.write_restart_frequency", "-1", "-md.no_restart", "true"] + list(extra)
    p = subprocess.run(cmd, cwd=d, capture_output=True, text=True)
    txt = re.sub(r"\x1b\[[0-9;]*m", "", p.stdout + p.stderr)
    diag, fj = os.path.join(d, "in.diag.jsonl"), os.path.join(d, "in.snapshots", "in.final.json")
    Q = [np.array(json.loads(l)["charges"]) for l in open(diag)] if os.path.exists(diag) else []
    avg = json.load(open(fj))["MD"]["average_Epot"] if os.path.exists(fj) else None
    return p.returncode, Q, avg, txt


# ---- 1. cn_zero_guard: device TERM 1b ----------------------------------------------------------------------------
print("1. cn_zero_guard on the device (Cl- ... H2O, CN(Cl) = 0)")
el_w, X_w = read_xyz(lib("water"))
o = el_w.index("O")
order = [o] + [i for i in range(3) if i != o]
el_w, X_w = [el_w[i] for i in order], X_w[order] - X_w[o]
v = X_w[1] - X_w[0]
X4 = np.vstack([X_w, X_w[1] + 3.2 * v / np.linalg.norm(v)]) + np.random.RandomState(11).normal(0, 0.06, (4, 3))
el_cl, X_cl = ["Cl"] + el_w, np.vstack([X4[3:], X4[:3]])          # Cl first: the reference rule gives it the -1
ON = REF + ["-gfnff.cn_zero_guard", "true"]
res = {}
for tag, extra in (("cpu off", REF), ("gpu off", REF + GPU), ("cpu on", ON), ("gpu on", ON + GPU)):
    res[tag] = sp(el_cl, X_cl, -1, extra, grad=True)
check(all(r[2] == 0 and r[1] is not None for r in res.values()), "four single points ran, GPU without fallback (rc %s)"
      % [r[2] for r in res.values()])
if all(r[1] is not None for r in res.values()):
    d_off = np.abs(res["gpu off"][1] - res["cpu off"][1]).max()
    d_on = np.abs(res["gpu on"][1] - res["cpu on"][1]).max()
    d_sw = np.abs(res["gpu on"][1] - res["gpu off"][1]).max()
    check(d_off < 1e-10 and d_on < 1e-10, "GPU against CPU gradient: %.1e (off), %.1e (on) Eh/A (< 1e-10)" % (d_off, d_on))
    check(d_sw > 1e-4, "the switch changes the GPU gradient by %.2e Eh/A (> 1e-4)" % d_sw)
    check(res["gpu on"][0] == res["gpu off"][0] == res["cpu on"][0], "energy unchanged by the switch and equal to the CPU")
    err, h = 0.0, 0.002
    G = res["gpu on"][1]
    for i, k in ((0, 0), (0, 1), (0, 2), (2, 0), (2, 1), (2, 2)):
        Y = X_cl.copy(); Y[i, k] += h; ep = sp(el_cl, Y, -1, ON + GPU)[0]
        Y = X_cl.copy(); Y[i, k] -= h; em = sp(el_cl, Y, -1, ON + GPU)[0]
        err = max(err, abs((ep - em) / (2 * h) - G[i, k]))
    check(err <= 2e-5, "GPU, switch on: max |FD - analytic| = %.2e Eh/A on Cl and H (<= 2e-5)" % err)

# ---- 2. charge groups: device EEQ during MD -----------------------------------------------------------------------
print("2. charge groups on the device EEQ (10 MD steps, GPU against CPU)")
frags = json.dumps({"A": src + "/src/molecules/pdmaema.xyz", "B": src + "/src/molecules/pdmaema_h.xyz"})
subprocess.run([curcuma, "-polymerbuild", "-fragments", frags, "-sequence", "(AB)2", "-verbosity", "0", "-no_bmt"],
               capture_output=True, text=True)
chain, chain_groups = os.path.join(work, "polymer_A_final.xyz"), os.path.join(work, "polymer_A_final.groups.json")
el_n, X_n = read_xyz(lib("ammonia.ideal-c3v"))
write("nh4cl.xyz", el_n + ["H", "Cl"], np.vstack([X_n, [[0.0, 0.0, 1.012], [0.0, 0.0, -3.3]]]))
nh4cl = os.path.join(work, "nh4cl.xyz")
cases = []
if os.path.exists(chain_groups):
    G1 = ["-gfnff.charge_groups", chain_groups]
    cases += [("chain, group file", chain, 2, REF + G1, json.load(open(chain_groups))["groups"]),
              ("chain, auto", chain, 2, REF + ["-gfnff.charge_groups", "auto"], None),
              ("chain, group file + cn_zero_guard", chain, 2, REF + G1 + ["-gfnff.cn_zero_guard", "true"],
               json.load(open(chain_groups))["groups"])]
else:
    check(False, "polymerbuild wrote the group file")
cases.append(("NH4+ Cl-, ionpairs", nh4cl, 0, ["-gfnff.charge_groups", "ionpairs"], [{"atoms": [5], "charge": -1.0}]))
for name, xyz, q, extra, groups in cases:
    rg, Qg, ag, tg = md(xyz, q, extra + GPU)
    rc, Qc, ac, _ = md(xyz, q, extra)
    ok = rg == 0 and rc == 0 and len(Qg) == len(Qc) and len(Qg) > 2 and ag is not None and ac is not None
    check(ok, "%s: both runs finished, GPU without fallback (rc %d / %d, %d frames)" % (name, rg, rc, len(Qg)))
    if not ok:
        continue
    dq = max(float(np.abs(a - b).max()) for a, b in zip(Qg, Qc))
    check(dq < 1e-9, "%s: charges GPU against CPU %.1e e (< 1e-9)" % (name, dq))
    check(abs(ag - ac) < 1e-8, "%s: mean Epot GPU against CPU %.1e Eh (< 1e-8)" % (name, abs(ag - ac)))
    check(max(abs(a.sum() - q) for a in Qg) < 1e-8, "%s: GPU total charge %d" % (name, q))
    if groups:
        dev = max(abs(a[g["atoms"]].sum() - g["charge"]) for a in Qg for g in groups)
        check(dev < 1e-8, "%s: GPU group sums, max |sum - Q| %.1e (< 1e-8)" % (name, dev))
    check("[GPU CUDA" in tg, "%s: device EEQ path reported (%s)" % (name, next((l.strip()[:110] for l in tg.split("\n")
          if "EEQ solver (final" in l), "no solver line")))

print(f"{count[0]} curcuma runs, {failed} failure(s)")
sys.exit(1 if failed else 0)
