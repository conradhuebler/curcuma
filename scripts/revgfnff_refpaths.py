#!/usr/bin/env python3
# Claude Generated (Sep 2026) - rev-gfnff WP2 item 9: relaxed reference paths
"""Relaxed r2SCAN-3c minimum-energy paths (ORCA NEB) for the BH76 RKT hydrogen-transfer
reactions.

Why this exists
---------------
The operator fixed the acceptance criterion for radical hydrogen transfer as the **path
shape** against the r2SCAN-3c reference (rms along the path + barrier position), not the
absolute barrier height (docs/REV_GFNFF_ROADMAP.md decision 6).  Only relaxed paths support
that criterion.  Before this script the campaign held relaxed NEB paths for six class-B
systems, of which two are RKT (`rkt06_h_h2`, `rkt14_h_oh_h2_o`); the other RKT entries in
`ref/B/` are either `*_neb_failed` or hand-interpolated single-point paths (`--b-mode
ts-points`), and the latter are not minimum-energy paths.

What it does (per RKT transition state, one reaction)
----------------------------------------------------
1.  Reads the reaction from `test_cases/GMTKN55-testset/BH76/.res`: the first entry names the
    two reactant fragments, the second the two product fragments (both are benchmark-supplied
    equilibrium structures in the same directory).
2.  Maps the reactant-fragment pair and the product-fragment pair onto the benchmark TS
    geometry (`map_onto_ts`: same-element permutation search filtered by the fragment's own
    bond graph, so a spectator-H swap cannot win on RMSD alone).
3.  Separates the two fragments of each side rigidly along their closest-contact vector to
    `SEPARATION` Angstrom, giving relaxed endpoints in the TS atom order (the fragments keep
    their own equilibrium geometries, so no endpoint relaxation run is needed).
4.  Runs ORCA `NEB-TS` with the benchmark TS as `NEB_TS_XYZFILE`, then an `EnGrad` multi-job
    pass over band + TS, in the same energies.json/gradients.json/points.xyz/meta.json layout
    as classes A/C/D/E/B (plus the path descriptors below).
5.  Records points, converged points, barrier height and barrier position (image, path
    fraction, value of the transfer coordinate xi = d(A-H*) - d(B-H*)), reaction energy,
    charge/multiplicity/UKS, <S**2> per point and the wall time.

Band convergence counts as success.  ORCA refuses to run the NEB-TS step when the band has no
interior maximum ("No barrier was found", nonzero return code) and writes the converged band
anyway - that is a property of the r2SCAN-3c surface for the very small BH76 H-transfer
barriers, not a driver failure, and the band is still a relaxed path.

Failure modes of the earlier class-B campaign (WP2_STATUS.md / WP2_LOST_SCANS.md), checked
before spending the budget, and what they mean here:
  * "no-barrier band" (rkt01, rkt10 and every other small-barrier step): the band converges
    correctly, r2SCAN-3c puts the maximum at the reactant endpoint.  Reproduced, kept, flagged.
  * "idpp failed to converge" (n2_h2_n2h2): happens only without a TS guess.  Every RKT
    reaction has a benchmark TS, so the seeded IDPP path converges in 0 iterations.
  * "TS optimisation not converged" (PX13, N2Hx chain): needs the NEB-TS step, which ORCA
    skips when there is no barrier.  Not reachable here.
  * `--b-mode ts-points` is NOT used: on rkt10 its chained `$new_job` MO carry-over landed on a
    different SCF solution for the TS point (28 kcal/mol high, <S**2> 0.7522 vs 0.7576), so
    those barrier numbers are single-point artefacts.  The same carry-over corrupts the band
    points here (rkt02: two images 19 and 46 kcal/mol high), so `run` flags them and
    `rescf` re-evaluates EVERY point as its own ORCA job, which is the canonical profile; the
    band's own and the chained values stay in `band_vs_fresh` for the branch comparison.

Usage:
    python scripts/revgfnff_refpaths.py plan
    python scripts/revgfnff_refpaths.py run --jobs 4 --nprocs 4
    python scripts/revgfnff_refpaths.py run --only RKT02 RKT03
    python scripts/revgfnff_refpaths.py status
    python scripts/revgfnff_refpaths.py rescf     # fresh single point per geometry, canonical
"""
import argparse
import gzip
import itertools
import json
import math
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import revgfnff_ref as R  # noqa: E402  (driver helpers: orca_input/parse_orca/run_orca/...)

REF = R.REF
P_DIR = REF / "P"
GMTKN = R.GMTKN
K = 627.5094740631  # Eh -> kcal/mol

SEPARATION = 4.5      # Angstrom between the two fragments of an endpoint
NIMAGES = 10          # images of the first NEB attempt
RETRY_IMAGES = 14     # second attempt
JUMP = 3.0            # kcal/mol: flag a band point this far off the line through its neighbours
KEEP = ["job.inp", "job.out.gz", "job_retry.inp", "job_retry.out.gz", "job.NEB.log",
        "job_retry.NEB.log", "job_MEP_trj.xyz", "job_retry_MEP_trj.xyz",
        "job_NEB-TS_converged.xyz", "job_retry_NEB-TS_converged.xyz"]
COV = {"H": 0.32, "C": 0.75, "N": 0.71, "O": 0.63, "F": 0.64, "P": 1.11, "S": 1.05, "Cl": 0.99}


# ------------------------------------------------------------------ .res / reaction list


def parse_bh76_rkt():
    """BH76/.res -> {TSNAME: [(frags..., ref_barrier_kcal), ...]} in file order.

    Line shape: `$tmer  h/$f  hcl/$f  RKT01/$f  x  -1  -1  1  $w  6.1` -- the last `$f` token
    is the TS, the ones before it are the fragments of that reaction entry, the trailing float
    is the benchmark barrier for that direction."""
    rxns = {}
    for line in (GMTKN / "BH76" / ".res").read_text().splitlines():
        if not line.startswith("$tmer"):
            continue
        toks = line.split()
        names = [t[:-3] for t in toks if t.endswith("/$f")]
        if len(names) < 2:
            continue
        ts = names[-1].upper()
        if not ts.startswith("RKT"):
            continue
        try:
            ref = float(toks[-1])
        except ValueError:
            ref = None
        rxns.setdefault(ts, []).append((tuple(names[:-1]), ref))
    return rxns


def _band_energies(path):
    """Per-frame energies from an ORCA NEB band trajectory (they sit in the comment line)."""
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    out, i = [], 0
    while i < len(lines):
        n = int(lines[i].split()[0])
        e = None
        for t in lines[i + 1].split():
            try:
                e = float(t)
            except ValueError:
                pass
        out.append(e)
        i += n + 2
    return out


def load_frag(name):
    d = GMTKN / "BH76" / name
    atoms = R.read_xyz(d / "struc.xyz")
    charge = int((d / ".CHRG").read_text().strip()) if (d / ".CHRG").exists() else 0
    return atoms, charge


# ------------------------------------------------------------------ fragment -> TS mapping


def _neighbours(atoms, scale=1.3):
    n = len(atoms)
    adj = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if math.dist(atoms[i][1:], atoms[j][1:]) < scale * (COV[atoms[i][0]] + COV[atoms[j][0]]):
                adj[i].add(j)
                adj[j].add(i)
    return adj


def map_onto_ts(parts, ts_atoms):
    """Map the fragments of one side onto `ts_atoms`; returns (ranked, n_valid) where `ranked`
    is the list of surviving mappings sorted by score, each (mean per-fragment Kabsch rmsd,
    inv, pieces).  `combine()` turns an entry into one geometry in TS slot order.

    Same-element permutations are enumerated and filtered by the fragments' own bond graphs
    (heavy-atom skeleton strictly, H-involving bonds generously, so a *forming* H-H bond at
    1.26 A survives), then scored by the MEAN PER-FRAGMENT Kabsch rmsd with each fragment
    placed independently on its TS atoms.  Per-fragment placement is what makes the score
    meaningful: the concatenated fragment list is not a geometry (every fragment sits at its
    own origin), so a single global Kabsch on it would score the arbitrary relative
    orientation and reject correct mappings.  The winning placement is the endpoint geometry,
    so the TS-like relative orientation comes out of the same step.
    """
    frag_atoms = [x for p in parts for x in p]
    if sorted(a[0] for a in frag_atoms) != sorted(a[0] for a in ts_atoms):
        raise ValueError("element multiset differs between fragments and TS")
    frag_bonds, off = [], 0
    for p in parts:
        adj = _neighbours(p)
        frag_bonds += [(i + off, j + off) for i in range(len(p)) for j in adj[i] if j > i]
        off += len(p)
    ts_pts = [a[1:] for a in ts_atoms]
    groups = []
    for e in sorted({a[0] for a in ts_atoms}):
        groups.append(([i for i, a in enumerate(ts_atoms) if a[0] == e],
                       [i for i, a in enumerate(frag_atoms) if a[0] == e]))
    n, nfrag = len(ts_atoms), len(frag_atoms)
    ranked = []
    for combo in itertools.product(*[list(itertools.permutations(f)) for _, f in groups]):
        inv = [0] * nfrag           # inv[frag index] = TS slot
        for (ts_idx, _), perm in zip(groups, combo):
            for slot, fidx in zip(ts_idx, perm):
                inv[fidx] = slot
        ok = True
        for i, j in frag_bonds:
            d = math.dist(ts_pts[inv[i]], ts_pts[inv[j]])
            lim = (COV[frag_atoms[i][0]] + COV[frag_atoms[j][0]]) + (1.2 if "H" in
                   (frag_atoms[i][0], frag_atoms[j][0]) else 0.35 * (COV[frag_atoms[i][0]] + COV[frag_atoms[j][0]]))
            if d > lim:
                ok = False
                break
        if not ok:
            continue
        pieces, score, off2 = [], 0.0, 0
        for p in parts:
            Qk = np.array([ts_pts[inv[off2 + i]] for i in range(len(p))], dtype=float)
            Pk = np.array([p[i][1:] for i in range(len(p))], dtype=float)
            al, r = R._kabsch(Pk, Qk)
            pieces.append([(p[i][0], float(al[i][0]), float(al[i][1]), float(al[i][2]))
                           for i in range(len(p))])
            score += r
            off2 += len(p)
        ranked.append((score / len(parts), list(inv), pieces))
    if not ranked:
        raise ValueError("no connectivity-preserving mapping found")
    ranked.sort(key=lambda t: t[0])
    return ranked, len(ranked)


def combine(inv, pieces, n):
    """Assemble per-fragment placements into one geometry in TS slot order (`inv` maps the
    fragment atom index -> TS slot).  This is the order every endpoint and the TS guess share,
    which the NEB needs."""
    combined = [None] * n
    off = 0
    for piece in pieces:
        for i in range(len(piece)):
            combined[inv[off + i]] = piece[i]
        off += len(piece)
    return combined


def separate(frag_lo, frag_hi, sep=SEPARATION):
    """Rigidly translate the second fragment along its closest-contact vector to the first so
    the minimum inter-fragment distance becomes `sep` Angstrom; returns (atoms, original dmin).
    Both arguments are fragment placements (not a TS-slot-ordered mix), which is why the
    callers keep the per-fragment pieces around instead of slicing the combined geometry."""
    dmin, ai, bi = None, None, None
    for i in range(len(frag_lo)):
        for j in range(len(frag_hi)):
            d = math.dist(frag_lo[i][1:], frag_hi[j][1:])
            if dmin is None or d < dmin:
                dmin, ai, bi = d, i, j
    v = [frag_hi[bi][c + 1] - frag_lo[ai][c + 1] for c in range(3)]
    nrm = math.sqrt(sum(c * c for c in v)) or 1.0
    shift = [(sep - dmin) * c / nrm for c in v]
    moved = [(s, x + shift[0], y + shift[1], z + shift[2]) for s, x, y, z in frag_hi]
    return frag_lo + moved, dmin


def _place(inv, pieces, n):
    """Separate the two fragments of one side and assemble the result in TS slot order."""
    lo, hi = pieces[0], pieces[1]
    merged, dmin = separate(lo, hi)
    n0 = len(lo)
    return combine(inv, [merged[:n0], merged[n0:]], n), dmin


def transfer_coordinate(reactant, product):
    """(h*, A, B): the transferring H, the atom it leaves (A) and the atom it goes to (B).
    Maximises  [d_react(B,h) - d_react(A,h)] + [d_prod(A,h) - d_prod(B,h)]  over H atoms and
    ordered pairs A != B -- i.e. the bond that breaks minus the bond that forms, in the same
    atom-index frame (both geometries are in TS atom order).  A != B is required: the earlier
    "nearest neighbour on each side" rule degenerates to xi == 0 whenever the H's nearest
    neighbour is the same atom at both ends."""
    best = None
    n = len(reactant)
    for i, a in enumerate(reactant):
        if a[0] != "H":
            continue
        for A in range(n):
            if A == i:
                continue
            for B in range(n):
                if B == i or B == A:
                    continue
                s = (math.dist(reactant[B][1:], reactant[i][1:]) - math.dist(reactant[A][1:], reactant[i][1:])
                     + math.dist(product[A][1:], product[i][1:]) - math.dist(product[B][1:], product[i][1:]))
                if best is None or s > best[0]:
                    best = (s, i, A, B)
    return best[1], best[2], best[3]


def xi(image, h, A, B):
    return math.dist(image[A][1:], image[h][1:]) - math.dist(image[B][1:], image[h][1:])


def _rmsd(a, b):
    return math.sqrt(sum(math.dist(a[k][1:], b[k][1:]) ** 2 for k in range(len(a))) / len(a))


# ------------------------------------------------------------------ job


class PathJob:
    """One relaxed NEB path for one RKT hydrogen-transfer reaction."""

    def __init__(self, ts_name, fwd, rev, nimages=NIMAGES):
        self.ts_name = ts_name
        self.fwd_frags, self.ref_fwd = fwd
        self.rev_frags, self.ref_rev = rev
        self.nimages = nimages
        self.system = ts_name.lower()
        self.uks = False
        self.broken_sym = False
        self.tag = (f"relaxed NEB path {ts_name}: {' + '.join(self.fwd_frags)} -> "
                    f"{' + '.join(self.rev_frags)} (fragments separated to {SEPARATION} A; "
                    f"BH76 barriers {self.ref_fwd}/{self.ref_rev} kcal/mol)")

    @property
    def dir(self):
        return P_DIR / self.system

    def done(self):
        e = self.dir / "energies.json"
        if not e.exists():
            return False
        try:
            return json.loads(e.read_text()).get("status") in ("ok", "failed", "blocked")
        except Exception:
            return False

    # -------------------------------------------------------------- build
    def build(self):
        ts_atoms, charge, mult = R.read_gmtkn("BH76", self.ts_name)
        n = len(ts_atoms)
        out = {}
        for side, frags in (("reactant", self.fwd_frags), ("product", self.rev_frags)):
            parts, qsum = [], 0
            for name in frags:
                a, q = load_frag(name)
                parts.append(a)
                qsum += q
            if qsum != charge:
                raise ValueError(f"{side} fragment charges sum to {qsum}, TS charge is {charge}")
            ranked, nvalid = map_onto_ts(parts, ts_atoms)
            score, inv, pieces = ranked[0]
            geom, dmin = _place(inv, pieces, n)
            out[side] = (geom, score, nvalid, dmin, ranked)
        # H + H2 -> H + H2 (RKT06): both .res entries name the same fragments, so the best
        # mapping would give reactant == product and a trivial band.  Take the best mapping
        # whose placement differs from the reactant's -- that is the H-exchange the benchmark
        # means.
        if [f.lower() for f in self.fwd_frags] == [f.lower() for f in self.rev_frags]:
            r_geom = out["reactant"][0]
            for score, inv, pieces in out["product"][4][1:]:
                geom, dmin = _place(inv, pieces, n)
                if _rmsd(r_geom, geom) > 1.0:
                    out["product"] = (geom, score, out["product"][2], dmin, out["product"][4])
                    break
        react, prod = out["reactant"], out["product"]
        if sorted(a[0] for a in react[0]) != sorted(a[0] for a in prod[0]):
            raise ValueError("reactant/product element multisets differ")
        return dict(ts=ts_atoms, reactant=react[0], product=prod[0], charge=charge, mult=mult,
                    r=(react[1], prod[1]), nvalid=(react[2], prod[2]),
                    dmin=(react[3], prod[3]),
                    rrmsd=round(_rmsd(react[0], prod[0]), 3))

    # -------------------------------------------------------------- rescf
    def rescf(self, nprocs, log):
        """Re-evaluate every path point as its OWN ORCA job (no `$new_job` MO carry-over) and
        make those energies canonical.

        Why: the multi-job `EnGrad` pass carries MOs from point to point, and for these radical
        H-transfer surfaces that can land a point on a different SCF solution (measured on
        rkt02: band images 5-7 came out 19 and 46 kcal/mol high with <S**2> 0.7535/0.7538
        instead of 0.7579/0.7601).  A fresh single point is deterministic and reproducible, so it
        is the defensible profile; the ORCA-NEB band's own per-image energies and the chained
        values are both kept next to it, which is exactly the branch comparison a consumer of
        radical UKS data needs.
        """
        d = self.dir
        f = d / "energies.json"
        if not f.exists():
            return "missing"
        dd = json.loads(f.read_text())
        if dd.get("status") != "ok":
            return dd.get("status")
        old = {p["point"]: p for p in dd["points"]}
        bandf = d / "job_MEP_trj.xyz" if (d / "job_MEP_trj.xyz").exists() else d / "job_retry_MEP_trj.xyz"
        bandE = _band_energies(bandf)
        frames = R.read_xyz_frames(d / "points.xyz")
        charge, mult = dd["charge"], dd["mult"]
        fresh = []
        for k, geom in enumerate(frames):
            d2 = d / "fresh" / f"pt{k}"
            d2.mkdir(parents=True, exist_ok=True)
            (d2 / "job.inp").write_text(R.orca_input([geom], charge, mult, nprocs, uks=mult > 1))
            R.run_orca(d2, "job.inp")
            out = (d2 / "job.out").read_text(errors="replace")
            (d2 / "job.out").unlink()
            rec = R.parse_orca(out)
            rec = rec[0] if rec else {"energy_eh": None, "gradient_eh_bohr": None, "s2": None}
            fresh.append(rec)
        nb = dd["n_band_frames"]
        cmp_rows = []
        for k in range(len(frames)):
            b = bandE[k] if k < len(bandE) else None
            c = old.get(k, {}).get("energy_eh")
            e = fresh[k]["energy_eh"]
            cmp_rows.append({"point": k, "band_eh": b, "chain_eh": c, "fresh_eh": e,
                             "band_minus_fresh_kcal": None if (b is None or e is None) else round((e - b) * K, 3),
                             "chain_minus_fresh_kcal": None if (c is None or e is None) else round((e - c) * K, 3),
                             "fresh_s2": fresh[k]["s2"], "old_s2": old.get(k, {}).get("s2")})
        E = [r["fresh_eh"] for r in cmp_rows]
        finite = [k for k in range(nb) if E[k] is not None]
        kmax = max(finite, key=lambda k: E[k]) if finite else None
        interior = kmax is not None and 0 < kmax < nb - 1
        g = R.read_xyz_frames(d / "points.xyz")
        h, A, B = dd["transfer_atoms"]["h"], dd["transfer_atoms"]["A"], dd["transfer_atoms"]["B"]
        pts = []
        for k, row in enumerate(cmp_rows):
            p = dict(old.get(k, {}))
            e = row["fresh_eh"]
            p.update({"energy_eh": e, "s2": row["fresh_s2"], "xi_ang": round(xi(g[k], h, A, B), 4)})
            p.pop("gradient_eh_bohr", None)
            pts.append(p)
        grads = []
        for k in range(len(frames)):
            grads.append(None if E[k] is None else fresh[k]["gradient_eh_bohr"])
        (d / "gradients.json").write_text(json.dumps({"unit": "Eh/Bohr", "gradients": grads}))
        with (d / "points.xyz").open("w") as fx:
            for k, atoms in enumerate(frames):
                fx.write(f"{len(atoms)}\nE={E[k] if E[k] is not None else 'nan'} charge={charge} "
                         f"mult={mult} point={k} role={pts[k].get('role')} xi={pts[k]['xi_ang']}\n")
                fx.write("".join(f"{s} {x:.8f} {y:.8f} {z:.8f}\n" for s, x, y, z in atoms))
        new = dict(dd)
        new.update({
            "energy_source": "fresh independent single point per geometry (no MO carry-over)",
            "n_converged": sum(1 for e in E if e is not None),
            "barrier_kcal": round((E[kmax] - E[0]) * K, 3) if kmax is not None and E[0] is not None else None,
            "barrier_image": kmax,
            "barrier_xi_ang": None if kmax is None else round(xi(g[kmax], h, A, B), 4),
            "no_interior_barrier": (not interior),
            "xi_reactant": round(xi(g[0], h, A, B), 4),
            "xi_product": round(xi(g[nb - 1], h, A, B), 4),
            "reaction_energy_kcal": round((E[nb - 1] - E[0]) * K, 3)
                if E[0] is not None and E[nb - 1] is not None else None,
            "ts_benchmark_energy_kcal": round((E[nb] - E[0]) * K, 3) if nb < len(E) and E[nb] is not None else None,
            "neb_ts_energy_kcal": round((E[nb + 1] - E[0]) * K, 3) if nb + 1 < len(E) and E[nb + 1] is not None else None,
            "band_vs_fresh": cmp_rows,
            "max_band_vs_fresh_kcal": max((abs(r["band_minus_fresh_kcal"]) for r in cmp_rows
                                           if r["band_minus_fresh_kcal"] is not None), default=None),
            "max_chain_vs_fresh_kcal": max((abs(r["chain_minus_fresh_kcal"]) for r in cmp_rows
                                            if r["chain_minus_fresh_kcal"] is not None), default=None),
            "points": pts,
        })
        (d / "energies.json").write_text(json.dumps(new, indent=1))
        log(f"  P {self.system}: rescf done, barrier {new['barrier_kcal']} at image {kmax}/{nb - 1}, "
            f"max|band-fresh| {new['max_band_vs_fresh_kcal']} kcal/mol")
        return "ok"

    # -------------------------------------------------------------- run
    def _attempt(self, d, tag, reactants, charge, mult, nprocs, nimages, maxiter):
        inp = R.orca_neb_input(reactants, charge, mult, nprocs, nimages, self.uks,
                               self.broken_sym, True, maxiter)
        (d / f"{tag}.inp").write_text(inp)
        rc, wall = R.run_orca(d, f"{tag}.inp", extra_keep=set(KEEP))
        out = (d / f"{tag}.out").read_text(errors="replace")
        with gzip.open(d / f"{tag}.out.gz", "wt") as gz:
            gz.write(out)
        (d / f"{tag}.out").unlink()
        if tag != "job":
            for src, dst in (("job_MEP_trj.xyz", f"{tag}_MEP_trj.xyz"),
                             ("job_NEB-TS_converged.xyz", f"{tag}_NEB-TS_converged.xyz"),
                             ("job.NEB.log", f"{tag}.NEB.log")):
                f = d / src
                if f.exists():
                    f.rename(d / dst)
        band = d / ("job_MEP_trj.xyz" if tag == "job" else f"{tag}_MEP_trj.xyz")
        return band, wall, rc, out, ("No barrier was found" in out), nimages

    def run(self, nprocs, log):
        d = self.dir
        d.mkdir(parents=True, exist_ok=True)
        try:
            b = self.build()
        except Exception as exc:
            (d / "energies.json").write_text(json.dumps(
                {"class": "P", "system": self.system, "tag": self.tag, "status": "blocked",
                 "n_points": 0, "n_ok": 0, "wall_s": 0.0, "error": f"build failed: {exc}"}, indent=1))
            log(f"  P {self.system}: BLOCKED (build) {exc}")
            return "blocked", 0, 0.0
        charge, mult = b["charge"], b["mult"]
        self.charge, self.mult = charge, mult
        self.uks = mult > 1
        ts, reactant, product = b["ts"], b["reactant"], b["product"]
        R.write_xyz(d / "reactant.xyz", reactant, f"{self.system} reactant ({SEPARATION} A)")
        R.write_xyz(d / "product.xyz", product, f"{self.system} product ({SEPARATION} A)")
        R.write_xyz(d / "ts_guess.xyz", ts, f"{self.system} benchmark TS (NEB_TS_XYZFILE)")
        build_note = (f"endpoints: {' + '.join(self.fwd_frags)} / {' + '.join(self.rev_frags)} "
                      f"benchmark fragment files mapped onto the TS (Kabsch rmsd "
                      f"{b['r'][0]:.4f}/{b['r'][1]:.4f} A, {b['nvalid'][0]}/{b['nvalid'][1]} valid "
                      f"permutations), fragments separated to {SEPARATION} A (closest contact "
                      f"{b['dmin'][0]:.3f}/{b['dmin'][1]:.3f} A in the TS)")
        t0 = time.time()
        band, wall, rc, out, nob, nim = self._attempt(d, "job", reactant, charge, mult, nprocs,
                                                      self.nimages, 300)
        attempts, no_barrier = 1, nob
        if not band.exists():
            log(f"  P {self.system}: attempt 1 produced no band (rc={rc}), retry with {RETRY_IMAGES} images")
            band, wall2, rc, out, nob, nim = self._attempt(d, "job_retry", reactant, charge, mult,
                                                           nprocs, RETRY_IMAGES, 400)
            wall += wall2
            attempts, no_barrier = 2, nob
        if not band.exists():
            (d / "energies.json").write_text(json.dumps(
                {"class": "P", "system": self.system, "tag": self.tag, "charge": charge, "mult": mult,
                 "uks": self.uks, "ref_kind": "UKS" if self.uks else "RKS", "status": "failed",
                 "n_points": 0, "n_ok": 0, "attempts": attempts, "wall_s": round(wall, 1),
                 "note": build_note, "orca_rc": rc}, indent=1))
            (d / "meta.json").write_text(json.dumps(R.meta(
                {"status": "failed", "attempts": attempts, "note": build_note, "orca_rc": rc,
                 "keywords": R.KEYWORDS + " NEB-TS", "nprocs": nprocs, "charge": charge, "mult": mult,
                 "uks": self.uks, "ref_kind": "UKS" if self.uks else "RKS"}), indent=1))
            log(f"  P {self.system}: FAILED (no band) after {attempts} attempt(s), {wall:.0f} s")
            return "failed", 0, wall

        images = R.read_xyz_frames(band)
        nb = len(images)
        refined = d / ("job_NEB-TS_converged.xyz" if attempts == 1 else "job_retry_NEB-TS_converged.xyz")
        points = list(images) + [ts] + ([R.read_xyz(refined)] if refined.exists() else [])
        roles = ["reactant"] + ["image"] * (nb - 2) + ["product", "ts_benchmark"] \
            + (["neb_ts"] if refined.exists() else [])

        eg_inp = R.orca_input(points, charge, mult, nprocs, uks=self.uks)
        (d / "job_engrad.inp").write_text(eg_inp)
        rc2, wall2 = R.run_orca(d, "job_engrad.inp", extra_keep=set(KEEP) | {"job_engrad.inp"})
        eg_out = (d / "job_engrad.out").read_text(errors="replace")
        with gzip.open(d / "job_engrad.out.gz", "wt") as gz:
            gz.write(eg_out)
        (d / "job_engrad.out").unlink()
        recs = R.parse_orca(eg_out)
        wall += wall2

        # Band points whose chained-$new_job MO carry-over may have landed on another SCF
        # solution: re-run them alone and keep both branches in the record.
        isolated, branches = [], {}
        for k in range(1, nb - 1):
            ek = recs[k - 1]["energy_eh"] if k - 1 < len(recs) else None
            ep = recs[k]["energy_eh"] if k < len(recs) else None
            en = recs[k + 1]["energy_eh"] if k + 1 < len(recs) else None
            if None in (ek, ep, en) or abs(ep - 0.5 * (ek + en)) * K > JUMP:
                isolated.append(k)
        for k in isolated:
            d2 = d / "isolated" / f"pt{k}"
            d2.mkdir(parents=True, exist_ok=True)
            (d2 / "job.inp").write_text(R.orca_input([points[k]], charge, mult, nprocs, uks=self.uks))
            R.run_orca(d2, "job.inp")
            t2 = (d2 / "job.out").read_text(errors="replace")
            (d2 / "job.out").unlink()
            r2 = R.parse_orca(t2)
            if r2 and r2[0]["energy_eh"] is not None:
                branches[str(k)] = {"chain_eh": recs[k]["energy_eh"], "fresh_eh": r2[0]["energy_eh"],
                                    "chain_s2": recs[k]["s2"], "fresh_s2": r2[0]["s2"],
                                    "d_kcal": round((r2[0]["energy_eh"] - recs[k]["energy_eh"]) * K, 3)}
                recs[k] = r2[0]

        h, A, B = transfer_coordinate(reactant, product)
        energies, ok_n, E = [], 0, []
        with (d / "points.xyz").open("w") as fx:
            for k, atoms in enumerate(points):
                rec = recs[k] if k < len(recs) else {"energy_eh": None, "gradient_eh_bohr": None, "s2": None}
                e = rec["energy_eh"] if rec.get("gradient_eh_bohr") is not None else None
                if e is not None:
                    ok_n += 1
                E.append(e)
                energies.append({"point": k, "label": f"role={roles[k]}_{k}", "role": roles[k],
                                 "energy_eh": e, "s2": rec.get("s2"),
                                 "xi_ang": round(xi(atoms, h, A, B), 4),
                                 "gradient_eh_bohr": rec.get("gradient_eh_bohr") if e is not None else None})
                fx.write(f"{len(atoms)}\nE={e if e is not None else 'nan'} charge={charge} mult={mult} "
                         f"point={k} role={roles[k]} xi={energies[-1]['xi_ang']}\n")
                fx.write("".join(f"{s} {x:.8f} {y:.8f} {z:.8f}\n" for s, x, y, z in atoms))
        (d / "gradients.json").write_text(json.dumps(
            {"unit": "Eh/Bohr", "gradients": [e["gradient_eh_bohr"] for e in energies]}))

        finite = [k for k in range(nb) if E[k] is not None]
        kmax = max(finite, key=lambda k: E[k]) if finite else None
        arc = [0.0]
        for k in range(1, nb):
            arc.append(arc[-1] + _rmsd(images[k - 1], images[k]))
        total_arc = arc[-1] if arc[-1] > 0 else 1.0
        interior = kmax is not None and 0 < kmax < nb - 1
        desc = {
            "b_mode": "neb-relaxed", "n_band_frames": nb, "n_points": len(points),
            "n_converged": ok_n, "attempts": attempts, "nimages": nim, "orca_rc": rc,
            "no_interior_barrier": (not interior), "orca_no_barrier_message": no_barrier,
            "barrier_kcal": round((E[kmax] - E[0]) * K, 3) if kmax is not None and E[0] is not None else None,
            "barrier_image": kmax,
            "barrier_path_fraction": round(arc[kmax] / total_arc, 4) if kmax is not None else None,
            "barrier_xi_ang": round(xi(images[kmax], h, A, B), 4) if kmax is not None else None,
            "xi_reactant": round(xi(images[0], h, A, B), 4),
            "xi_product": round(xi(images[nb - 1], h, A, B), 4),
            "reaction_energy_kcal": round((E[nb - 1] - E[0]) * K, 3)
                if E[0] is not None and E[nb - 1] is not None else None,
            "ts_benchmark_energy_kcal": round((E[nb] - E[0]) * K, 3) if nb < len(E) and E[nb] is not None else None,
            "neb_ts_energy_kcal": round((E[nb + 1] - E[0]) * K, 3) if nb + 1 < len(E) and E[nb + 1] is not None else None,
            "ts_to_path_rmsd_ang": round(min(_rmsd(ts, images[k]) for k in range(nb)), 4),
            "path_arc_ang": round(total_arc, 4),
            "transfer_atoms": {"h": h, "A": A, "B": B},
            "s2": [e["s2"] for e in energies],
            "isolated_reruns": isolated, "branch_differences": branches,
            "charge": charge, "mult": mult, "uks": self.uks,
            "ref_kind": "UKS" if self.uks else "RKS",
            "bh76_ref_barrier_kcal": {"fwd": self.ref_fwd, "rev": self.ref_rev},
            "wall_s": round(wall, 1), "note": build_note,
        }
        for e in energies:
            e.pop("gradient_eh_bohr")
        (d / "energies.json").write_text(json.dumps(
            {"class": "P", "system": self.system, "tag": self.tag, "status": "ok",
             **desc, "points": energies}, indent=1))
        m = R.meta({"status": "ok", "note": build_note, "keywords": R.KEYWORDS + " NEB-TS then EnGrad",
                    "nprocs": nprocs, "charge": charge, "mult": mult, "uks": self.uks,
                    "ref_kind": "UKS" if self.uks else "RKS", "nimages": nim, "attempts": attempts,
                    "no_interior_barrier": (not interior)})
        m["script"] = "scripts/revgfnff_refpaths.py"
        (d / "meta.json").write_text(json.dumps(m, indent=1))
        log(f"  P {self.system}: band {nb} frames, {ok_n}/{len(points)} points ok, "
            f"barrier {desc['barrier_kcal']} kcal/mol at image {kmax}/{nb - 1}, {wall:.0f} s")
        return "ok", ok_n, wall


# ------------------------------------------------------------------ driver


def build_jobs(only, log):
    rxns = parse_bh76_rkt()
    jobs = []
    for ts, entries in sorted(rxns.items()):
        if ts == "RKT22":  # C5H8 isomerisation, not a hydrogen transfer
            continue
        if len(entries) < 2:
            log(f"  {ts}: {len(entries)} .res entry, skipped")
            continue
        if only and ts not in only:
            continue
        jobs.append(PathJob(ts, entries[0], entries[1]))
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["plan", "run", "status", "rescf"])
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--nprocs", type=int, default=4)
    args = ap.parse_args()

    def log(msg):
        print(msg, flush=True)

    jobs = build_jobs(args.only, log)
    if args.cmd == "plan":
        for j in jobs:
            print(f"  {'done' if j.done() else 'todo'} P {j.system:26s} {j.tag}")
        return
    if args.cmd == "status":
        for j in jobs:
            f = j.dir / "energies.json"
            if not f.exists():
                print(f"  {j.system:26s} not started")
                continue
            dd = json.loads(f.read_text())
            print(f"  {j.system:26s} {dd.get('status')} barrier={dd.get('barrier_kcal')} "
                  f"img={dd.get('barrier_image')}/{dd.get('n_band_frames', 0) - 1} "
                  f"mult={dd.get('mult')} ok={dd.get('n_ok')}/{dd.get('n_points')}")
        return

    if args.cmd == "rescf":
        todo = [j for j in jobs if (j.dir / "energies.json").exists()
                and json.loads((j.dir / "energies.json").read_text()).get("status") == "ok"
                and not (j.dir / "energies.json").read_text().find('"energy_source"') >= 0]
        log(f"rescf: {len(todo)} reactions")
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            list(ex.map(lambda j: j.rescf(args.nprocs, log), todo))
        log("rescf finished")
        return

    todo = [j for j in jobs if not j.done()]
    log(f"{len(jobs)} reactions planned, {len(todo)} to run, {args.jobs} x {args.nprocs} cores")

    def work(j):
        try:
            j.run(args.nprocs, log)
        except Exception as exc:
            log(f"  P {j.system}: FAILED {exc!r}")

    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        list(ex.map(work, todo))
    log("campaign finished")


if __name__ == "__main__":
    main()
