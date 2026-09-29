#!/usr/bin/env python3
"""
Q5 (FABLE_BOND_STATE_2.md section 7) offline verification -- item 1 and item 5 of the
orchestrator's task: a structural count, NOT an energy calculation.

Reuses collect_refset()/run_curcuma()/BOND_RE from scripts/revgfnff_bondgate_sweep.py (no
modification to that file) to get the perceived GFN-FF/rev-gfnff topology (BOND i(Zi)-j(Zj)
lines, CURCUMA_BONDDUMP=1) for every one of the 2647 reference-set structures (GMTKN55 2462 +
MOR41 95 + S30L-CI 90), and asks two purely structural questions per structure:

  (A) does any atom with Z==1 (hydrogen) have exactly two bonded partners in the perceived
      topology ("a two-coordinate H")?
  (B) does any 3-ring (triangle i-j-k, all three pairs bonded) contain at least one Z==1 atom
      ("an H-bridged triangle")?

This is exactly Fable's Section 7.3 claim: "80 structures have a hydrogen with two listed
partners ... 38 of them with an H-bridged triangle. Everything else is untouched." No energy,
no rule application -- just the count, verified against the actual binary's topology
perception (release/curcuma, md5 53b32b06 at time of writing).
"""
import concurrent.futures as cf
import itertools
import sys
import tempfile
import shutil
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import revgfnff_bondgate_sweep as gate  # noqa: E402


def structure_h_scope(label, xyz, charge, spin):
    with tempfile.TemporaryDirectory() as td:
        local = Path(td) / "struc.xyz"
        shutil.copy(xyz, local)
        out = gate.run_curcuma(local, td, charge=charge, spin=spin)
    if out is None:
        return label, dict(ok=False, note="TIMEOUT")
    lines = gate.clean_lines(out)
    pairs = set()
    Z = {}
    for ln in lines:
        m = gate.BOND_RE.match(ln)
        if m:
            i, zi, j, zj = int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4))
            pairs.add((min(i, j), max(i, j)))
            Z[i] = zi
            Z[j] = zj
    if not pairs:
        return label, dict(ok=False, note="NO BOND line found (parse failure or crash)")

    # adjacency
    adj = {}
    for (i, j) in pairs:
        adj.setdefault(i, set()).add(j)
        adj.setdefault(j, set()).add(i)

    two_coord_h = [i for i in adj if Z.get(i) == 1 and len(adj[i]) == 2]

    # 3-ring detection: triangle i<j<k all pairwise bonded, containing >=1 H
    h_bridged_triangle = False
    triangle_atoms = set()
    atoms = sorted(adj.keys())
    for i in atoms:
        neighbours = sorted(n for n in adj[i] if n > i)
        for a, b in itertools.combinations(neighbours, 2):
            if b in adj.get(a, ()):
                tri = (i, a, b)
                if any(Z.get(x) == 1 for x in tri):
                    h_bridged_triangle = True
                    triangle_atoms.add(tri)

    return label, dict(ok=True, two_coord_h=two_coord_h, n_two_coord_h=len(two_coord_h),
                        h_bridged_triangle=h_bridged_triangle, triangles=sorted(triangle_atoms),
                        n_atoms=len(atoms), n_bonds=len(pairs))


def main():
    jobs = gate.collect_refset(limit=0)
    print("H-scope structural sweep: %d structures" % len(jobs), flush=True)
    results = []
    done = 0
    with cf.ThreadPoolExecutor(max_workers=24) as ex:
        futs = {ex.submit(structure_h_scope, label, xyz, c, s): label
                for (label, xyz, c, s) in jobs}
        for fut in cf.as_completed(futs):
            label = futs[fut]
            try:
                lbl, r = fut.result()
            except Exception as e:
                lbl, r = label, dict(ok=False, note="EXC: %s" % e)
            results.append((lbl, r))
            done += 1
            if done % 500 == 0:
                print("  ... %d/%d" % (done, len(jobs)), flush=True)

    results.sort(key=lambda x: x[0])
    n_fail = sum(1 for _, r in results if not r["ok"])
    with_h2 = [(lbl, r) for lbl, r in results if r["ok"] and r["n_two_coord_h"] > 0]
    with_tri = [(lbl, r) for lbl, r in with_h2 if r["h_bridged_triangle"]]

    print("\n=== SUMMARY ===")
    print("total structures        : %d" % len(results))
    print("parse/run failures       : %d" % n_fail)
    print("structures w/ 2-coord H  : %d" % len(with_h2))
    print("  of those, w/ H-triangle: %d" % len(with_tri))

    print("\n=== per-structure list (2-coord H) ===")
    for lbl, r in with_h2:
        print("%-40s n_two_coord_h=%d atoms=%s triangle=%s"
              % (lbl, r["n_two_coord_h"], r["two_coord_h"], r["h_bridged_triangle"]))

    if n_fail:
        print("\n=== FAILURES ===")
        for lbl, r in results:
            if not r["ok"]:
                print("%-40s %s" % (lbl, r.get("note", "?")))

    out = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if out:
        import json
        out.write_text(json.dumps(
            {lbl: {k: v for k, v in r.items() if k != "triangles" or True} for lbl, r in results},
            default=list, indent=0))
        print("\nwrote %s" % out)


if __name__ == "__main__":
    main()
