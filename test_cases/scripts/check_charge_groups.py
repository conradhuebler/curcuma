#!/usr/bin/env python3
"""Claude Generated (Oct 2026): GFN-FF EEQ charge groups (-gfnff.charge_groups) on a protonated PDMAEMA (AB)2 chain.

Uses only the curcuma binary. Builds the chain (polymerbuild also writes the group file), then checks
  1. every group sums to its prescribed charge, the total to the molecule charge (the ungrouped EEQ does not);
  2. analytic gradient against a central finite difference of the own energy (perturbed geometry);
  3. broken input (an atom in two groups) is rejected: error message, energy equal to the ungrouped one.
Usage: check_charge_groups.py CURCUMA_BINARY SOURCE_DIR
"""
import json, os, re, subprocess, sys, tempfile
import numpy as np

curcuma, src = sys.argv[1], sys.argv[2]
failed = 0

def check(ok, what):
    global failed
    print(("  PASS  " if ok else "  FAIL  ") + what)
    failed += 0 if ok else 1

work = tempfile.mkdtemp(prefix="curcuma_charge_groups_")
os.chdir(work)
frags = json.dumps({"A": src + "/src/molecules/pdmaema.xyz", "B": src + "/src/molecules/pdmaema_h.xyz"})
subprocess.run([curcuma, "-polymerbuild", "-fragments", frags, "-sequence", "(AB)2", "-verbosity", "0", "-no_bmt"],
               capture_output=True, text=True, check=True)
check(os.path.exists("polymer_A_final.groups.json"), "polymerbuild wrote polymer_A_final.groups.json")
groups = json.load(open("polymer_A_final.groups.json"))["groups"]
check(len(groups) == 2 and all(g["charge"] == 1.0 for g in groups), "two groups of charge +1 (one per protonated unit)")

lines = open("polymer_A_final.xyz").read().split("\n")
n = int(lines[0])
el = [l.split()[0] for l in lines[2:2 + n]]
X0 = np.array([[float(v) for v in l.split()[1:4]] for l in lines[2:2 + n]])
np.random.seed(1)
X = X0 + np.random.normal(0, 0.06, X0.shape)

def write(path, X):
    with open(path, "w") as f:
        f.write(f"{n}\nx\n")
        for i in range(n):
            f.write("%s %.10f %.10f %.10f\n" % (el[i], *X[i]))

def sp(path, groups_file=None, dump=None, grad=None):
    cmd = [curcuma, "-sp", path, "-method", "gfnff", "-charge", "2", "-verbosity", "0", "-no_bmt"]
    if groups_file:
        cmd += ["-gfnff.charge_groups", groups_file]
    if dump:
        cmd += ["-dump_atomic_charges", dump]
    if grad:
        cmd += ["-dump_gradient", grad]
    out = subprocess.run(cmd, capture_output=True, text=True)
    m = re.search(r"Single Point Energy = (-?\d+\.\d+)", out.stdout)
    return (float(m.group(1)) if m else None), out.stdout + out.stderr

write("base.xyz", X0)
e_free, _ = sp("base.xyz", None, "q_free.txt")
e_grp, _ = sp("base.xyz", "polymer_A_final.groups.json", "q_grp.txt")
qf, qg = np.loadtxt("q_free.txt"), np.loadtxt("q_grp.txt")
sums_free = [qf[g["atoms"]].sum() for g in groups]
sums_grp = [qg[g["atoms"]].sum() for g in groups]
check(all(abs(s - 1.0) < 1e-6 for s in sums_grp), "grouped EEQ: group sums %s" % ", ".join("%.6f" % s for s in sums_grp))
check(abs(qg.sum() - 2.0) < 1e-6, "grouped EEQ: total charge %.6f" % qg.sum())
check(all(abs(s - 1.0) > 0.1 for s in sums_free), "ungrouped EEQ spreads the charge (group sums %s)" % ", ".join("%.3f" % s for s in sums_free))
check(e_free is not None and e_grp is not None and abs(e_free - e_grp) > 1e-4, "groups change the energy (%.6f against %.6f Eh)" % (e_grp, e_free))

# gradient against a central finite difference of the own energy, perturbed geometry
write("pert.xyz", X)
e0, _ = sp("pert.xyz", "polymer_A_final.groups.json", None, "pert.grad")
G = np.loadtxt("pert.grad", comments="#")
inside = groups[0]["atoms"]
outside = [i for i in range(n) if i not in set(sum([g["atoms"] for g in groups], []))]
worst, scale = 0.0, 0.0
h = 0.002
for i in (inside[0], inside[10], outside[0]):
    for k in range(3):
        Xp, Xm = X.copy(), X.copy()
        Xp[i, k] += h
        Xm[i, k] -= h
        write("p.xyz", Xp)
        write("m.xyz", Xm)
        fd = (sp("p.xyz", "polymer_A_final.groups.json")[0] - sp("m.xyz", "polymer_A_final.groups.json")[0]) / (2 * h)
        worst = max(worst, abs(fd - G[i, k]))
        scale = max(scale, abs(G[i, k]))
check(worst < 5e-5, "gradient against finite difference: max |diff| %.2e Eh/A (largest gradient component %.3f)" % (worst, scale))

# broken input
json.dump({"groups": [{"atoms": [0, 1, 2], "charge": 1.0}, {"atoms": [2, 3], "charge": 1.0}]}, open("bad.json", "w"))
e_bad, text = sp("base.xyz", "bad.json")
check("ignored" in text and e_bad is not None and abs(e_bad - e_free) < 1e-9, "atom in two groups: rejected, energy equals the ungrouped one")

print("FAILED" if failed else "all passed")
sys.exit(1 if failed else 0)
