#!/usr/bin/env python3
"""rev_excess_bond_extend adversarial checks (X2_COMPRESSED_SURVEY_STATUS.md section 8.4): does the
2c-3e candidate rule fire where it must not (SN2 path, anion...HX, water across the contact, q = -2)
and still fire where it should (Cl2- bare and with a water at one end)? Usage: this.py BIN
Prints fired count and dE(recX - rec) in kcal/mol; exit status = number of unexpected cases.
Claude Generated (Sep 27, 2026). AI-generated, machine-tested only."""
import json, os, subprocess, sys, tempfile, math
from concurrent.futures import ThreadPoolExecutor
B = sys.argv[1]
REC = ["-gfnff.rev_charge_model", "sqe", "-gfnff.rev_sqe_phase1", "true", "-gfnff.rev_excess_electron", "true",
       "-gfnff.rev_excess_mode", "harris", "-gfnff.frag_charge_model", "ensemble", "-gfnff.frag_charge_s_max", "1.2",
       "-gfnff.rev_sqe_virtual_pairs", "true", "-gfnff.rev_sqe_group_pairs_only", "true",
       "-gfnff.frag_charge_atomic_ea", "true", "-gfnff.rev_pi_excess_electron", "true"]
X = REC + ["-gfnff.rev_excess_bond_extend", "1.8"]
T = os.environ.get("X2S_SCRATCH", "/var/tmp/x2s") + "/adv"; os.makedirs(T, exist_ok=True)


def run(atoms, q, flags, verb=2):
    d = tempfile.mkdtemp(dir=T)
    open(d + "/m.xyz", "w").write(f"{len(atoms)}\n\n" + "\n".join(f"{a} {x:.6f} {y:.6f} {z:.6f}" for a, x, y, z in atoms) + "\n")
    p = subprocess.run([B, "-sp", "m.xyz", "-method", "revgfnff", "-batch", "true", "-batch_out", "o.jsonl",
                        "-gfnff.cache_topology", "false", "-charge", str(q), "-threads", "1", "-verbosity", str(verb),
                        "-no_bmt"] + flags, cwd=d, capture_output=True, text=True)
    e = json.loads(open(d + "/o.jsonl").readline())["energy_eh"]
    return e, p.stdout.count("2c-3e candidate bond")


def water(ox, oy, oz, toward):  # water with one H pointing along -x (toward the anion) if toward
    s = -1.0 if toward else 1.0
    return [("O", ox, oy, oz), ("H", ox + s * 0.96, oy, oz), ("H", ox - 0.24 * s, oy + 0.93, oz)]


cases = {}
# must NOT fire
cases["Cl- ... Cl-, q=-2, 4.0 A"] = ([("Cl", 0, 0, 0), ("Cl", 0, 0, 4.0)], -2, False)
cases["Cl- ... Cl, q=-1, water ON the axis between (O at mid)"] = ([("Cl", 0, 0, 0), ("Cl", 0, 0, 4.4),
    ("O", 0, 0, 2.2), ("H", 0.76, 0, 2.2 + 0.59), ("H", -0.76, 0, 2.2 + 0.59)], -1, False)
cases["F- ... H-F (F bonded)"] = ([("F", 0, 0, 0), ("H", 0, 0, 1.6), ("F", 0, 0, 2.5)], -1, False)
cases["Cl- ... CH3Cl entrance complex"] = ([("Cl", 0, 0, -3.2), ("C", 0, 0, 0), ("H", 1.03, 0, -0.35), ("H", -0.515, 0.892, -0.35),
    ("H", -0.515, -0.892, -0.35), ("Cl", 0, 0, 1.80)], -1, False)
cases["Cl- ... O=O (not a same-row pair)"] = ([("Cl", 0, 0, 0), ("O", 0, 0, 3.3), ("O", 0, 0, 4.5)], -1, False)
# SN2 path Cl...CH3...Cl, asymmetric scan (the TS of BH76/clch3clts is the symmetric point)
for d1 in (2.0, 2.3, 2.6, 3.0, 3.5):
    d2 = 4.6 - d1 if d1 < 2.3 else 4.6 - 2.3 + (d1 - 2.3) * 0.2
    cases[f"SN2 [Cl..CH3..Cl]- d1 {d1:.1f} d2 {d2:.2f}"] = ([("Cl", 0, 0, -d1), ("C", 0, 0, 0), ("H", 1.08, 0, 0), ("H", -0.54, 0.935, 0),
        ("H", -0.54, -0.935, 0), ("Cl", 0, 0, d2)], -1, False)
# MUST fire
for r in (2.7, 3.3, 4.0):
    cases[f"Cl2- {r} A bare"] = ([("Cl", 0, 0, 0), ("Cl", 0, 0, r)], -1, True)
    cases[f"Cl2- {r} A + water H-bonded at one end (X..H 2.2 A on axis)"] = ([("Cl", 0, 0, 0), ("Cl", 0, 0, r),
        ("H", 0, 0, r + 2.2), ("O", 0, 0, r + 3.16), ("H", 0.93, 0, r + 3.40)], -1, True)
    cases[f"Cl2- {r} A + water beside one end (perpendicular, X..H 2.2 A)"] = ([("Cl", 0, 0, 0), ("Cl", 0, 0, r),
        ("H", 2.2, 0, r), ("O", 3.16, 0, r), ("H", 3.40, 0.93, r)], -1, True)

with ThreadPoolExecutor(20) as ex:
    J = {k: (ex.submit(run, a, q, REC, 0), ex.submit(run, a, q, X, 2), want) for k, (a, q, want) in cases.items()}
bad = 0
for k, (a, b, want) in J.items():
    e0, _ = a.result(); e1, n = b.result()
    ok = (n > 0) == want
    bad += not ok
    print(f"{'ok ' if ok else 'BAD'} {k:62s} fired {n}  dE {(e1 - e0) * 627.5094740631:+9.3f}")
print("unexpected:", bad)
sys.exit(bad)
