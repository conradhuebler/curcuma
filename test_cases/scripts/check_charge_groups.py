#!/usr/bin/env python3
"""Claude Generated (Oct 2026): GFN-FF EEQ charge groups (-gfnff.charge_groups) on a protonated PDMAEMA (AB)2 chain.

Uses only the curcuma binary. Builds the chain (polymerbuild also writes the group file), then checks
  1. every group sums to its prescribed charge, the total to the molecule charge (the ungrouped EEQ does not);
  2. analytic gradient against a central finite difference of the own energy (perturbed geometry);
  3. broken input (an atom in two groups) is rejected: error message, energy equal to the ungrouped one;
  4. -gfnff.charge_groups auto (onium perception): two groups of 14 atoms (N, N-H, three alpha-C, their eight H,
     beta-C), each inside its protonated unit, sums exactly 1; the same groups and charges for a permuted atom order;
     neutral molecule / charge 1 with two centres: no groups, energy equal to the ungrouped one; gradient against a
     central finite difference.
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

def sp(path, groups_file=None, dump=None, grad=None, charge=2, verbosity=0):
    # the topology cache <basename>.topo.json does not carry the molecule charge: a file written at charge 2 is
    # reused at charge 0 (measured: identical energy). Every run here therefore builds its topology afresh.
    cmd = [curcuma, "-sp", path, "-method", "gfnff", "-charge", str(charge), "-verbosity", str(verbosity), "-no_bmt",
           "-gfnff.cache_topology", "false"]
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

# ---- auto mode (onium perception) ----
def auto_groups(text):
    """group atom lists printed by the auto mode at verbosity >= 2"""
    clean = re.sub(r"\x1b\[[0-9;]*m", "", text)
    return [sorted(int(a) for a in m.group(1).split()) for m in re.finditer(r"charge_groups auto group \d+: centre \d+ atoms([ \d]+)", clean)]

e_auto, text = sp("base.xyz", "auto", "q_auto.txt", verbosity=2)
qa = np.loadtxt("q_auto.txt")
ag = auto_groups(text)
check("2 onium centres" in text and len(ag) == 2 and all(len(g) == 14 for g in ag),
      "auto: two groups of 14 atoms (sizes %s)" % [len(g) for g in ag])
elem_ok = all(sorted(el[a] for a in g) == sorted(["N"] + ["C"] * 4 + ["H"] * 9) for g in ag)
check(elem_ok, "auto: each group is N + 4 C + 9 H (N-H, three alpha-C, eight alpha-H, beta-C)")
check(all(any(set(g) <= set(u["atoms"]) for u in groups) for g in ag), "auto: each group lies inside one protonated unit of the polymerbuild group file")
sums_auto = [qa[g].sum() for g in ag]
check(all(abs(s - 1.0) < 1e-6 for s in sums_auto) and abs(qa.sum() - 2.0) < 1e-6,
      "auto: group sums %s, total %.6f" % (", ".join("%.8f" % s for s in sums_auto), qa.sum()))

# (ii) atom order: permute the xyz, the perceived groups and the charges must map onto the original ones
perm = np.random.RandomState(7).permutation(n)      # new atom k is old atom perm[k]
el_orig = el
el = [el_orig[p] for p in perm]
write("perm.xyz", X0[perm])
e_auto_p, text_p = sp("perm.xyz", "auto", "q_auto_p.txt", verbosity=2)
e_free_p, _ = sp("perm.xyz", None, "q_free_p.txt")
el = el_orig
qap = np.loadtxt("q_auto_p.txt")
agp = sorted(sorted(int(perm[a]) for a in g) for g in auto_groups(text_p))
check(agp == sorted(ag), "auto, permuted atom order: the same groups (mapped back)")
qfp = np.loadtxt("q_free_p.txt")
sums_auto_p = [qap[g].sum() for g in auto_groups(text_p)]
check(all(abs(s - 1.0) < 1e-6 for s in sums_auto_p), "auto, permuted atom order: group sums %s" % ", ".join("%.8f" % s for s in sums_auto_p))
# GFN-FF itself is not invariant under renumbering (Known Issue #34): on this chain the UNGROUPED energy moves by
# ~6.5e-5 Eh and the charges by ~7e-4 e (measured Oct 2026, one permutation). The groups are identical (checked
# above), so the test asks that the groups add no order dependence of their own: the permutation-induced change with
# groups stays within 1e-5 Eh / 2e-4 e of the change without (measured 2.7e-6 Eh, 8e-5 e).
dq_auto = np.abs(qap - qa[perm]).max()
dq_free = np.abs(qfp - qf[perm]).max()
check(abs(dq_auto - dq_free) < 2e-4, "auto, permuted atom order: max |dq| %.1e e with groups, %.1e e without (#34)" % (dq_auto, dq_free))
dE_auto, dE_free = e_auto_p - e_auto, e_free_p - e_free
check(abs(dE_auto - dE_free) < 1e-5,
      "auto, permuted atom order: dE %.2e Eh with groups, %.2e Eh without (#34)" % (dE_auto, dE_free))

# (iii) neutral molecule: condition (d) finds no centre; charge 1 with two centres: rejected
e_auto0, text0 = sp("base.xyz", "auto", charge=0, verbosity=1)
e_free0, _ = sp("base.xyz", None, charge=0)
check("0 onium centres" in text0 and e_auto0 is not None and abs(e_auto0 - e_free0) < 1e-9,
      "auto, charge 0: no centres, energy equals the ungrouped one")
# charge 1 with two centres, closed shell: remove one H from a carbon outside both groups (carbanion, the second,
# invisible charge carrier; same electron count as the +2 chain); the base chain at charge 1 is open shell and
# condition (a) would block it before the count rule is reached
in_ag = set(sum(ag, []))
h_out = next(i for i in range(n) if el[i] == "H" and i not in in_ag)
keep = [i for i in range(n) if i != h_out]
with open("minus_h.xyz", "w") as f:
    f.write(f"{n - 1}\nx\n")
    for i in keep:
        f.write("%s %.10f %.10f %.10f\n" % (el[i], *X0[i]))
e_auto1, text1 = sp("minus_h.xyz", "auto", charge=1, verbosity=1)
e_free1, _ = sp("minus_h.xyz", None, charge=1)
check("auto ignored" in text1 and e_auto1 is not None and abs(e_auto1 - e_free1) < 1e-9,
      "auto, charge 1 with two centres: rejected, energy equals the ungrouped one")
# charge 1, open shell (spin 0 requested, odd electron count): condition (a), no centres
e_auto1o, text1o = sp("base.xyz", "auto", charge=1, verbosity=1)
check("0 onium centres" in text1o, "auto, odd electron count: no centres (condition a)")

# several fragments: library structure hydronium-water2 (H3O+(H2O)2, three GFN-FF fragments) -> rejected
h3o = subprocess.run([sys.executable, src + "/scripts/structlib.py", "path", "hydronium-water2"],
                     capture_output=True, text=True).stdout.strip()
if h3o and os.path.exists(h3o):
    e_h3o_auto, text_h3o = sp(h3o, "auto", charge=1, verbosity=1)
    e_h3o_free, _ = sp(h3o, None, charge=1)
    check("auto ignored" in text_h3o and "fragments" in text_h3o and e_h3o_auto is not None
          and abs(e_h3o_auto - e_h3o_free) < 1e-9, "auto, H3O+(H2O)2 (3 fragments): rejected, energy equals the ungrouped one")
else:
    check(False, "library structure hydronium-water2 not found")

# (iv) gradient against a central finite difference with auto, perturbed geometry
e0, _ = sp("pert.xyz", "auto", None, "pert_auto.grad")
G = np.loadtxt("pert_auto.grad", comments="#")
in_auto = set(sum(ag, []))
outside_auto = [i for i in range(n) if i not in in_auto]
worst, scale = 0.0, 0.0
for i in (ag[0][0], ag[1][-1], outside_auto[0]):
    for k in range(3):
        Xp, Xm = X.copy(), X.copy()
        Xp[i, k] += h
        Xm[i, k] -= h
        write("p.xyz", Xp)
        write("m.xyz", Xm)
        fd = (sp("p.xyz", "auto")[0] - sp("m.xyz", "auto")[0]) / (2 * h)
        worst = max(worst, abs(fd - G[i, k]))
        scale = max(scale, abs(G[i, k]))
check(worst < 5e-5, "auto: gradient against finite difference: max |diff| %.2e Eh/A (largest gradient component %.3f)" % (worst, scale))

print("FAILED" if failed else "all passed")
sys.exit(1 if failed else 0)
