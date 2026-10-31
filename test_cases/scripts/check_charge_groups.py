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
  5. EXPERIMENTAL charge_groups ionpairs (parity rule) on NH4+ Cl-, NH4+ CH3O- and pyridinium Cl- built at run time
     from library structures: fragment sums +1/-1 (and +-0.9), the energy equals that of the same groups given as a
     file, permuted atom order, gradient against a central finite difference (12 components, 6 per fragment), no effect
     for one fragment / a charged system / an ambiguous anion count, parity keeps a neutral NH3 neutral,
     ionpair_aromatic only with the switch and not for neutral pyridine.
  6. file key "rest": "per_fragment" on the (AB)2 chain in 24 waters with 2 Cl- (built at run time, -ionize,
     scripts/charge_groups_for_solvated.py): unit sums +1, Cl- -1, every water 0; no groups: units not +1; gradient
     against a central finite difference (12 components); rejections (partly covered fragment, wrong sum, atom in two
     groups, unknown key value); permuted atom order; "rest": "single" equals the file without the key.
  7. stored library structure pdmaema-ab2-h2-water24-cl2.grad-defect with its group file (rest per_fragment): the Cl- whose
     nearest H is 3.18 A away sits at GFN-FF CN = 0 exactly and carries q = -1; the analytic gradient on it and on that H
     against a central finite difference (tolerance 5e-5 Eh/A; 0.93 / 1.05 Eh/A before Oct 31, 2026: Coulomb TERM 1b
     q*cnf/(2*sqrt(0) + 1e-16), Known Issue #43), and -gfnff.cn_zero_guard true changes nothing with groups (always on).
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

def sp(path, groups_file=None, dump=None, grad=None, charge=2, verbosity=0, extra=()):
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
    cmd += list(extra)
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

# ---- (v) EXPERIMENTAL -gfnff.charge_groups ionpairs: formal charges on the whole fragments of a neutral ion pair ----
# The ion pairs are built here at run time from library structures (no new structure file): NH4+ Cl- from the ideal
# C3v ammonia (fourth H on the C3 axis at 1.012 A, Cl 3.3 A below N, nearest H...Cl 3.07 A), NH4+ CH3O- from the same
# NH4 and the methanol without its hydroxyl H, and pyridinium chloride from the pyridine of fi-pyridine (F and I
# removed, H on N at 1.01 A, Cl 3.0 A beyond that H). Not optimised: test geometries, not physical ones.
# The rule: cation = fragment with an onium centre (+q), anions = fragments without a centre with an odd neutral
# electron count (-q each, exactly as many as centres), even-count fragments neutral; one EEQ group per fragment.
def lib(id_):
    return subprocess.run([sys.executable, src + "/scripts/structlib.py", "path", id_], capture_output=True, text=True).stdout.strip()

def read_xyz(path):
    L = open(path).read().split("\n")
    m = int(L[0])
    return [l.split()[0] for l in L[2:2 + m]], np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + m]])

def write_mol(path, els, XX):
    with open(path, "w") as f:
        f.write(f"{len(els)}\nx\n")
        for e_, x_ in zip(els, XX):
            f.write("%s %.10f %.10f %.10f\n" % (e_, *x_))

def ip_line(text):
    m = re.search(r"GFNFF charge_groups ionpairs: ([^\n]*)", re.sub(r"\x1b\[[0-9;]*m", "", text))
    return m.group(1) if m else ""

IP = ["-gfnff.charge_groups", "ionpairs"]
el_nh3, X_nh3 = read_xyz(lib("ammonia.ideal-c3v"))
el_ip = el_nh3 + ["H", "Cl"]
X_ip = np.vstack([X_nh3, [[0.0, 0.0, 1.012], [0.0, 0.0, -3.3]]])
write_mol("nh4cl.xyz", el_ip, X_ip)
e_ip, t_ip = sp("nh4cl.xyz", None, "q_ip.txt", charge=0, verbosity=1, extra=IP)
e_ip0, _ = sp("nh4cl.xyz", None, "q_ip0.txt", charge=0)
q_ip, q_ip0 = np.loadtxt("q_ip.txt"), np.loadtxt("q_ip0.txt")
check("1:5/+1.00" in ip_line(t_ip) and "2:1/-1.00" in ip_line(t_ip), "ionpairs, NH4+ Cl-: two fragments, +1 on NH4, -1 on Cl (%s)" % ip_line(t_ip)[:160])
check("NOT for reactions with a shared proton" in t_ip, "ionpairs: the run states that the rule is not meant for shared-proton reactions")
check(abs(q_ip[:5].sum() - 1.0) < 1e-6 and abs(q_ip[5] + 1.0) < 1e-6 and abs(q_ip0[:5].sum()) < 1e-6,
      "ionpairs, NH4+ Cl-: fragment sums %.8f / %.8f (default %.2e / %.2e)" % (q_ip[:5].sum(), q_ip[5], q_ip0[:5].sum(), q_ip0[5]))
check(e_ip is not None and e_ip0 is not None and abs(e_ip - e_ip0) > 1e-3, "ionpairs changes the energy (%.6f against %.6f Eh)" % (e_ip, e_ip0))
# the same groups written as a file (one group per fragment) give the identical energy: one mechanism
json.dump({"groups": [{"atoms": [0, 1, 2, 3, 4], "charge": 1.0}, {"atoms": [5], "charge": -1.0}]}, open("nh4cl.groups.json", "w"))
e_ipf, _ = sp("nh4cl.xyz", "nh4cl.groups.json", charge=0)
check(e_ipf == e_ip, "ionpairs == the same fragment groups from a file (%.8f / %.8f Eh)" % (e_ip, e_ipf))
e_ip9, t_ip9 = sp("nh4cl.xyz", None, "q_ip9.txt", charge=0, verbosity=1, extra=IP + ["-gfnff.ionpair_charge", "0.9"])
q_ip9 = np.loadtxt("q_ip9.txt")
check(abs(q_ip9[:5].sum() - 0.9) < 1e-6 and abs(q_ip9[5] + 0.9) < 1e-6, "ionpair_charge 0.9: fragment sums %.8f / %.8f" % (q_ip9[:5].sum(), q_ip9[5]))
# permutation: same roles (mapped back) and no order dependence of its own (Known Issue #34)
perm_ip = np.array([5, 2, 0, 4, 1, 3])  # new atom k is old atom perm_ip[k]
write_mol("nh4cl_p.xyz", [el_ip[p] for p in perm_ip], X_ip[perm_ip])
e_ipp, t_ipp = sp("nh4cl_p.xyz", None, "q_ipp.txt", charge=0, verbosity=1, extra=IP)
e_ipp0, _ = sp("nh4cl_p.xyz", None, None, charge=0)
q_ipp = np.loadtxt("q_ipp.txt")
check("1:1/-1.00" in ip_line(t_ipp) and "2:5/+1.00" in ip_line(t_ipp) and np.abs(q_ipp - q_ip[perm_ip]).max() < 2e-4,
      "ionpairs, permuted atom order: Cl first, the same fragment charges (max |dq| %.1e e)" % np.abs(q_ipp - q_ip[perm_ip]).max())
check(abs((e_ipp - e_ip) - (e_ipp0 - e_ip0)) < 1e-5, "ionpairs, permuted: dE %.2e Eh with the rule, %.2e Eh without" % (e_ipp - e_ip, e_ipp0 - e_ip0))
# gradient against a central finite difference: NH4+ CH3O- (methoxide: 17 electrons when neutral, odd), perturbed;
# N and one H of the cation, O and C of the anion: 6 components per fragment
el_me, X_me = read_xyz(lib("methanol"))
iO = el_me.index("O")
iHO = min((i for i, e_ in enumerate(el_me) if e_ == "H"), key=lambda i: np.linalg.norm(X_me[i] - X_me[iO]))
keep_me = [i for i in range(len(el_me)) if i != iHO]
el_mo, X_mo = [el_me[i] for i in keep_me], X_me[keep_me]
X_mo = X_mo - X_mo[el_mo.index("O")] + np.array([0.0, 0.0, 1.012 + 3.0])  # O 3.0 A beyond the axial H of NH4
el_am = el_ip[:5] + el_mo
X_am = np.vstack([X_ip[:5], X_mo]) + np.random.RandomState(3).normal(0, 0.06, (5 + len(el_mo), 3))
write_mol("nh4_meo.xyz", el_am, X_am)
e_am, t_am = sp("nh4_meo.xyz", None, "q_am.txt", "am.grad", charge=0, verbosity=1, extra=IP)
q_am = np.loadtxt("q_am.txt")
check("1:5/+1.00" in ip_line(t_am) and "2:5/-1.00" in ip_line(t_am) and abs(q_am[:5].sum() - 1.0) < 1e-6,
      "ionpairs, NH4+ CH3O-: +1 / -1 (%s)" % ip_line(t_am)[:120])
G = np.loadtxt("am.grad", comments="#")
iC = 5 + el_mo.index("C")
iOm = 5 + el_mo.index("O")
worst, scale, ncomp = 0.0, 0.0, 0
for i in (0, 4, iOm, iC):
    for k in range(3):
        Xp_, Xm_ = X_am.copy(), X_am.copy()
        Xp_[i, k] += h
        Xm_[i, k] -= h
        write_mol("p.xyz", el_am, Xp_)
        write_mol("m.xyz", el_am, Xm_)
        fd = (sp("p.xyz", None, charge=0, extra=IP)[0] - sp("m.xyz", None, charge=0, extra=IP)[0]) / (2 * h)
        worst = max(worst, abs(fd - G[i, k]))
        scale = max(scale, abs(G[i, k]))
        ncomp += 1
check(worst < 5e-5, "ionpairs, NH4+ CH3O-: gradient against finite difference, %d components: max |diff| %.2e Eh/A (largest component %.3f)" % (ncomp, worst, scale))
# no effect: one fragment (NH4+ alone), a charged system (NH4+ + NH3 at +1), ambiguous (NH4+ + 3 Cl at charge 0: 3 odd
# candidates for 1 centre); parity: NH4+ Cl- + NH3 fires with the NH3 neutral (even electron count)
write_mol("nh4.xyz", el_ip[:5], X_ip[:5])
e1, t1 = sp("nh4.xyz", None, charge=1, verbosity=1, extra=IP)
e1d, _ = sp("nh4.xyz", None, charge=1)
check("fragment(s)" in ip_line(t1) and "no effect" in ip_line(t1) and e1 == e1d, "ionpairs, NH4+ alone (one fragment): no effect, energy identical")
X_far = X_nh3 + np.array([8.0, 0.0, 0.0])
write_mol("nh4_nh3.xyz", el_ip[:5] + el_nh3, np.vstack([X_ip[:5], X_far]))
e2, t2 = sp("nh4_nh3.xyz", None, charge=1, verbosity=1, extra=IP)
e2d, _ = sp("nh4_nh3.xyz", None, charge=1)
check("charge 1" in ip_line(t2) and "no effect" in ip_line(t2) and e2 == e2d, "ionpairs, NH4+ + NH3 at +1 (charged): no effect, energy identical")
write_mol("nh4cl3.xyz", el_ip + ["Cl", "Cl"], np.vstack([X_ip, [[8.0, 0.0, 0.0], [-8.0, 0.0, 0.0]]]))
e3, t3 = sp("nh4cl3.xyz", None, charge=0, verbosity=1, extra=IP)
e3d, _ = sp("nh4cl3.xyz", None, charge=0)
check("ambiguous" in ip_line(t3) and e3 == e3d, "ionpairs, NH4+ + 3 Cl at charge 0 (3 candidates for 1 centre): ambiguous, energy identical")
write_mol("nh4cl_nh3.xyz", el_ip + el_nh3, np.vstack([X_ip, X_far]))
e4, t4 = sp("nh4cl_nh3.xyz", None, "q_3.txt", charge=0, verbosity=1, extra=IP)
q_3 = np.loadtxt("q_3.txt")
check("3:4/+0.00" in ip_line(t4) and abs(q_3[:5].sum() - 1.0) < 1e-6 and abs(q_3[5] + 1.0) < 1e-6 and abs(q_3[6:].sum()) < 1e-6,
      "ionpairs, NH4+ Cl- + NH3: parity keeps the NH3 neutral, sums %.6f / %.6f / %.6f" % (q_3[:5].sum(), q_3[5], q_3[6:].sum()))
# ionpair_aromatic: pyridinium chloride fires only with the switch; neutral pyridine + HCl never
el_fp, X_fp = read_xyz(lib("fi-pyridine"))
keep_fp = [i for i, e_ in enumerate(el_fp) if e_ not in ("F", "I")]
el_py, X_py = [el_fp[i] for i in keep_fp], X_fp[keep_fp]
iN = el_py.index("N")
write_mol("pyhcl.xyz", el_py + ["H", "Cl"], np.vstack([X_py, X_py[iN] + [0, 0, 1.01], X_py[iN] + [0, 0, 4.01]]))
eA_, tA_ = sp("pyhcl.xyz", None, charge=0, verbosity=1, extra=IP)
eB_, tB_ = sp("pyhcl.xyz", None, "q_py.txt", charge=0, verbosity=1, extra=IP + ["-gfnff.ionpair_aromatic", "true"])
e0_, _ = sp("pyhcl.xyz", None, charge=0)
q_py = np.loadtxt("q_py.txt")
check("no cation fragment" in ip_line(tA_) and eA_ == e0_, "ionpairs without ionpair_aromatic, pyridinium Cl-: no centre, energy identical")
check("aromatic rings 1" in ip_line(tB_) and abs(q_py[:-1].sum() - 1.0) < 1e-6 and abs(q_py[-1] + 1.0) < 1e-6,
      "ionpair_aromatic, pyridinium Cl-: one ring, fragment sums %.8f / %.8f" % (q_py[:-1].sum(), q_py[-1]))
el_hcl, X_hcl = read_xyz(lib("hydrogen-chloride"))
write_mol("py_hcl.xyz", el_py + el_hcl, np.vstack([X_py, X_hcl + X_py[iN] + [0, 0, 4.0]]))
eC_, tC_ = sp("py_hcl.xyz", None, charge=0, verbosity=1, extra=IP + ["-gfnff.ionpair_aromatic", "true"])
eC0, _ = sp("py_hcl.xyz", None, charge=0)
check("aromatic rings 0" in ip_line(tC_) and eC_ == eC0, "ionpair_aromatic, neutral pyridine + HCl (6 pi electrons): no centre, energy identical")

# ---- (vi) solvated chain: file key "rest": "per_fragment" (Claude Generated, Oct 2026) ----
# Built at run time (no structure file): the (AB)2 chain of section 1 (+2), 26 copies of the library water (random
# rotations, lattice points 2.9 A or more from the chain), then -ionize replaces 2 waters by Cl- (total charge 0), then
# scripts/charge_groups_for_solvated.py writes the groups: the two unit groups, the rest of the chain (charge 0), each Cl-
# (-1). GFN-FF makes every water its own group with charge 0.
el_w, X_w = read_xyz(lib("water"))
X_w = X_w - X_w[el_w.index("O")]
rng = np.random.RandomState(11)
def rotation():
    qv = rng.normal(size=4)
    qv /= np.linalg.norm(qv)
    w_, x_, y_, z_ = qv
    return np.array([[1 - 2 * (y_ * y_ + z_ * z_), 2 * (x_ * y_ - z_ * w_), 2 * (x_ * z_ + y_ * w_)],
                     [2 * (x_ * y_ + z_ * w_), 1 - 2 * (x_ * x_ + z_ * z_), 2 * (y_ * z_ - x_ * w_)],
                     [2 * (x_ * z_ - y_ * w_), 2 * (y_ * z_ + x_ * w_), 1 - 2 * (x_ * x_ + y_ * y_)]])
el_s, X_s = list(el), X0.copy()
lo, hi = X0.min(0) - 5.0, X0.max(0) + 5.0
grid = np.mgrid[lo[0]:hi[0]:3.1, lo[1]:hi[1]:3.1, lo[2]:hi[2]:3.1].reshape(3, -1).T
dg = np.linalg.norm(grid[:, None] - X0[None], axis=2).min(1)
grid = grid[dg >= 2.9]
grid = grid[np.argsort(np.linalg.norm(grid[:, None] - X0[None], axis=2).min(1), kind="stable")]
for c in grid:
    if len(el_s) == n + 3 * 26:
        break
    W3 = X_w @ rotation().T + c
    if np.linalg.norm(W3[:, None] - X_s[None], axis=2).min() >= 2.2:
        el_s += el_w
        X_s = np.vstack([X_s, W3])
write_mol("solv.xyz", el_s, X_s)
subprocess.run([curcuma, "-ionize", "solv.xyz", "-charge", "2", "-min_solute_distance", "2.5", "-min_ion_distance", "5",
                "-output", "solv_ion.xyz", "-no_bmt", "-verbosity", "0"], capture_output=True, text=True)
check(os.path.exists("solv_ion.xyz"), "solvated chain: -ionize wrote solv_ion.xyz")
el_i, X_i = read_xyz("solv_ion.xyz")
n_i = len(el_i)
i_cl = [i for i in range(n_i) if el_i[i] == "Cl"]
n_wat = sum(1 for e_ in el_i if e_ == "O") - sum(1 for e_ in el if e_ == "O")
check(len(i_cl) == 2 and n_wat >= 20 and n_i == n + 3 * n_wat + 2,
      "solvated chain: %d atoms = chain %d + %d water + %d Cl-" % (n_i, n, n_wat, len(i_cl)))
cg = subprocess.run([sys.executable, src + "/scripts/charge_groups_for_solvated.py", "solv_ion.xyz", "--polymer-groups",
                     "polymer_A_final.groups.json", "--polymer-offset", "0", "--anion", "Cl", "--charge", "0", "--out", "solv.groups.json"],
                    capture_output=True, text=True)
check(cg.returncode == 0 and os.path.exists("solv.groups.json"), "charge_groups_for_solvated.py: %s" % cg.stdout.strip().split("\n")[-1])
sg = json.load(open("solv.groups.json"))
check(sg.get("rest") == "per_fragment" and len(sg["groups"]) == 5, "group file: rest per_fragment, 5 groups (2 units, chain rest, 2 Cl-)")
cg_bad = subprocess.run([sys.executable, src + "/scripts/charge_groups_for_solvated.py", "solv_ion.xyz", "--polymer-groups",
                         "polymer_A_final.groups.json", "--anion", "Cl", "--charge", "1", "--out", "x.json"], capture_output=True, text=True)
check(cg_bad.returncode != 0 and not os.path.exists("x.json"), "charge_groups_for_solvated.py: wrong --charge is an error, nothing written")
# (i) group sums: units +1 exactly, Cl- -1, every water 0
e_s, t_s = sp("solv_ion.xyz", "solv.groups.json", "q_s.txt", charge=0, verbosity=1)
q_s = np.loadtxt("q_s.txt")
unit_sums = [q_s[g_["atoms"]].sum() for g_ in sg["groups"][:2]]
wat_sums = [q_s[i:i + 3].sum() for i in range(n, n_i) if el_i[i] == "O"]
check("ignored" not in t_s and all(abs(s - 1.0) < 1e-6 for s in unit_sums) and all(abs(q_s[i] + 1.0) < 1e-6 for i in i_cl)
      and max(abs(s) for s in wat_sums) < 1e-6 and abs(q_s.sum()) < 1e-6,
      "rest per_fragment: unit sums %s, Cl %s, %d waters max |sum| %.1e, total %.1e" % (
          ", ".join("%.8f" % s for s in unit_sums), ", ".join("%.6f" % q_s[i] for i in i_cl), len(wat_sums), max(abs(s) for s in wat_sums), q_s.sum()))
# (ii) without groups the units do not carry +1 (every fragment of a neutral system is 0 in GFN-FF, Cl- included)
e_s0, _ = sp("solv_ion.xyz", None, "q_s0.txt", charge=0)
q_s0 = np.loadtxt("q_s0.txt")
check(e_s0 is not None and all(abs(q_s0[g_["atoms"]].sum() - 1.0) > 0.5 for g_ in sg["groups"][:2]),
      "no groups: unit sums %s, Cl %s (E %.6f, with groups %.6f Eh)" % (", ".join("%.3f" % q_s0[g_["atoms"]].sum() for g_ in sg["groups"][:2]),
                                                                        ", ".join("%.3f" % q_s0[i] for i in i_cl), e_s0, e_s))
# (iii) gradient against a central finite difference, perturbed (sigma 0.06 A): 12 components, 4 in a unit group, 2 on a
# Cl-, 4 on a water (O and H), 2 on the chain rest
X_ip_s = X_i + np.random.RandomState(5).normal(0, 0.06, X_i.shape)
write_mol("solv_p.xyz", el_i, X_ip_s)
sp("solv_p.xyz", "solv.groups.json", None, "solv_p.grad", charge=0)
G = np.loadtxt("solv_p.grad", comments="#")
u0 = sg["groups"][0]["atoms"]
o_w = next(i for i in range(n, n_i) if el_i[i] == "O")
comps = [(u0[0], 0), (u0[0], 1), (u0[5], 0), (u0[5], 1), (i_cl[0], 0), (i_cl[0], 1), (o_w, 0), (o_w, 1), (o_w + 1, 0), (o_w + 1, 2),
         (sg["groups"][2]["atoms"][0], 0), (sg["groups"][2]["atoms"][0], 1)]
worst, scale = 0.0, 0.0
for i, k in comps:
    Xp_, Xm_ = X_ip_s.copy(), X_ip_s.copy()
    Xp_[i, k] += h
    Xm_[i, k] -= h
    write_mol("p.xyz", el_i, Xp_)
    write_mol("m.xyz", el_i, Xm_)
    fd = (sp("p.xyz", "solv.groups.json", charge=0)[0] - sp("m.xyz", "solv.groups.json", charge=0)[0]) / (2 * h)
    worst = max(worst, abs(fd - G[i, k]))
    scale = max(scale, abs(G[i, k]))
check(worst < 5e-5, "rest per_fragment: gradient against finite difference, %d components: max |diff| %.2e Eh/A (largest %.3f)" % (len(comps), worst, scale))
# (iv) rejections: run continues without groups, energy equal to the ungrouped one
def rejected(groups_json, tag):
    json.dump(groups_json, open("rej.json", "w"))
    e_r, t_r = sp("solv_ion.xyz", "rej.json", charge=0)
    check("ignored" in t_r and e_r == e_s0, "rest per_fragment rejects %s, energy equals the ungrouped one" % tag)
    return t_r
t_part = rejected({"rest": "per_fragment", "groups": sg["groups"][:2] + sg["groups"][3:]}, "a partly covered chain (rest group left out)")
check("covered completely or not at all" in t_part, "partly covered: message names the rule")
rejected({"rest": "per_fragment", "groups": sg["groups"] + [{"atoms": [o_w], "charge": 0.0}]}, "a partly covered water (only its O)")
t_sum = rejected({"rest": "per_fragment", "groups": sg["groups"][:4]}, "a wrong charge sum (one Cl- left out: +1)")
check("sum to the molecule charge" in t_sum, "wrong sum: message names the rule")
rejected({"rest": "per_fragment", "groups": sg["groups"] + [{"atoms": [i_cl[0]], "charge": 0.0}]}, "an atom in two groups")
rejected({"rest": "every_atom", "groups": sg["groups"]}, "an unknown value of rest")
# (v) permutation of the atom order: same group sums, energy change as without groups (Known Issue #34)
perm_s = np.random.RandomState(9).permutation(n_i)   # new atom k is old atom perm_s[k]
inv_s = np.argsort(perm_s)
write_mol("solv_perm.xyz", [el_i[p] for p in perm_s], X_i[perm_s])
json.dump({"rest": "per_fragment", "groups": [{"atoms": sorted(int(inv_s[a]) for a in g_["atoms"]), "charge": g_["charge"]} for g_ in sg["groups"]]},
          open("solv_perm.groups.json", "w"))
e_sp, _ = sp("solv_perm.xyz", "solv_perm.groups.json", "q_sp.txt", charge=0)
e_sp0, _ = sp("solv_perm.xyz", None, None, charge=0)
q_sp = np.loadtxt("q_sp.txt")[inv_s]   # back to the original numbering
sums_p = [q_sp[g_["atoms"]].sum() for g_ in sg["groups"]] + [q_sp[i:i + 3].sum() for i in range(n, n_i) if el_i[i] == "O"]
sums_o = [q_s[g_["atoms"]].sum() for g_ in sg["groups"]] + wat_sums
check(max(abs(a - b) for a, b in zip(sums_p, sums_o)) < 1e-6,
      "rest per_fragment, permuted atom order: group and water sums equal (max diff %.1e e)" % max(abs(a - b) for a, b in zip(sums_p, sums_o)))
check(abs((e_sp - e_s) - (e_sp0 - e_s0)) < 1e-5,
      "rest per_fragment, permuted: dE %.2e Eh with groups, %.2e Eh without (#34), max |dq| %.1e e" % (e_sp - e_s, e_sp0 - e_s0, np.abs(q_sp - q_s).max()))
# (vi) "rest": "single" is the default: the same file with the key spelled out gives the identical energy (one-fragment chain)
json.dump({"rest": "single", "groups": groups}, open("single.json", "w"))
e_single, _ = sp("base.xyz", "single.json")
check(e_single == e_grp, "rest single (explicit) == groups without the key (%.8f / %.8f Eh)" % (e_single, e_grp))

# (vii) Claude Generated (Oct 31, 2026): the stored defect structure (library id pdmaema-ab2-h2-water24-cl2.grad-defect,
# group file test_cases/reference_data/charge_groups/) - Cl- 180 at CN = 0 exactly (nearest H 98 at 3.18 A) with q = -1
# from the groups. Before the fix the analytic gradient on it was (0.928, -0.049, -1.025) against FD (-0.009, -0.009, 0.021)
# Eh/A: Coulomb TERM 1b q*cnf/(2*sqrt(CN) + 1e-16) at CN = 0 (Known Issue #43); the fix applies the cn_zero_guard
# whenever a fragment-charge override is active. Components: Cl- 180 (x, y, z) and its nearest H 98 (x, z).
el_d, X_d = read_xyz(lib("pdmaema-ab2-h2-water24-cl2.grad-defect"))
groups_d = src + "/test_cases/reference_data/charge_groups/pdmaema-ab2-h2-water24-cl2.groups.json"
write_mol("defect.xyz", el_d, X_d)
e_d, _ = sp("defect.xyz", groups_d, "q_d.txt", "defect.grad", charge=0)
q_d = np.loadtxt("q_d.txt")
G = np.loadtxt("defect.grad", comments="#")
check(el_d[180] == "Cl" and abs(q_d[180] + 1.0) < 1e-6, "defect structure: atom 180 is a Cl- with q = %.6f from the groups" % q_d[180])
worst, scale = 0.0, 0.0
for i, k in ((180, 0), (180, 1), (180, 2), (98, 0), (98, 2)):
    Xp_, Xm_ = X_d.copy(), X_d.copy()
    Xp_[i, k] += h
    Xm_[i, k] -= h
    write_mol("p.xyz", el_d, Xp_)
    write_mol("m.xyz", el_d, Xm_)
    fd = (sp("p.xyz", groups_d, charge=0)[0] - sp("m.xyz", groups_d, charge=0)[0]) / (2 * h)
    worst = max(worst, abs(fd - G[i, k]))
    scale = max(scale, abs(G[i, k]))
check(worst < 5e-5, "defect structure, rest per_fragment, Cl- at CN = 0: gradient against finite difference, 5 components: max |diff| %.2e Eh/A (largest %.3f)" % (worst, scale))
e_dg, _ = sp("defect.xyz", groups_d, None, "defect_guard.grad", charge=0, extra=["-gfnff.cn_zero_guard", "true"])
G_g = np.loadtxt("defect_guard.grad", comments="#")
check(e_dg == e_d and np.array_equal(G, G_g), "defect structure: -gfnff.cn_zero_guard true changes neither energy nor gradient with groups (guard always on with an override)")

print("FAILED" if failed else "all passed")
sys.exit(1 if failed else 0)
