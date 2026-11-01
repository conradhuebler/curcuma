#!/usr/bin/env python3
"""Where did the good conformers come from? -- provenance analysis of a ConfSearch run.

Since Aug 2026 every structure carries where it was born:

    cycle02_T550K_r3.gfn2#001 <- cycle02_T550K_r3_w3_i5027_t120
                                 stage        rep walker hill  time [fs]

This script reads that back and answers the questions a search cannot answer while it
runs: which walker was worth its cores, and at what point of a trajectory the useful
structures actually appear.

MEASURED ON THE FIRST RUN THAT CARRIED PROVENANCE (107-atom peptide, gfn2/gfn2,
1747 structures with provenance):

  deposit time of all structures      median 880 fs
  deposit time of the 200 deepest     median 220 fs
  deposit time of the  10 deepest     median  55 fs

  i.e. the deeper a structure, the earlier it was deposited. Every MD starts at an
  already optimised seed and the RMSD bias drives it away from there, so the first
  steps still sit near a good minimum while later ones explore breadth at higher
  energy. The long trajectory buys coverage, not depth.

  yield per walker: the best structure came from walker 3 (+1.0 kJ/mol against the
  reference), walker 9 never got below +42.7 -- at comparable structure counts
  (121 to 232). A factor of 40 in quality between seeds of the same phase.

WHAT IT DOES NOT SHOW: which seed will be productive. Walker 3 was not the
lowest-energy seed of its phase but the seventh. The analysis is diagnostic, not
predictive -- treat a single run as one sample.

Usage:
    confsearch_provenance.py <run-directory-or-cumulative.xyz> [reference.xyz]

With a reference structure the energies are reported relative to it; without, relative
to the deepest structure of the run.
"""

import os
import re
import sys
from collections import defaultdict

HARTREE_TO_KJ = 2625.4996
PROVENANCE = re.compile(r"<- (cycle\d+)_T(\d+)K_r(\d+)_w(-?\d+)_i(\d+)_t(\d+)")
ENERGY = re.compile(r"Energy\s*=\s*(-?\d+\.\d+)")


def read_frames(path):
    """(energy, comment) per frame -- geometry is not needed here."""
    out = []
    with open(path) as handle:
        lines = handle.read().split("\n")
    i = 0
    while i < len(lines) and lines[i].strip().isdigit():
        n = int(lines[i])
        comment = lines[i + 1]
        match = ENERGY.search(comment)
        if match:
            out.append((float(match.group(1)), comment))
        i += n + 2
    return out


def median(values):
    if not values:
        return float("nan")
    ordered = sorted(values)
    mid = len(ordered) // 2
    return ordered[mid] if len(ordered) % 2 else 0.5 * (ordered[mid - 1] + ordered[mid])


def quantile(values, q):
    if not values:
        return float("nan")
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(q * len(ordered)))]


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    target = sys.argv[1]
    path = target
    if os.path.isdir(target):
        candidates = [f for f in os.listdir(target) if f.endswith("cumulative.opt.xyz")]
        if not candidates:
            print(f"no cumulative.opt.xyz in {target}")
            return 1
        path = os.path.join(target, candidates[0])

    frames = read_frames(path)
    if not frames:
        print(f"no structures in {path}")
        return 1

    reference = None
    if len(sys.argv) > 2:
        ref_frames = read_frames(sys.argv[2])
        if ref_frames:
            reference = ref_frames[0][0]

    rows = []          # (energy, stage, temperature, repetition, walker, hill, time)
    without = 0
    for energy, comment in frames:
        match = PROVENANCE.search(comment)
        if not match:
            without += 1
            continue
        rows.append((energy, match.group(1), int(match.group(2)), int(match.group(3)),
                     int(match.group(4)), int(match.group(5)), int(match.group(6))))

    print(f"{os.path.basename(path)}: {len(rows)} structure(s) with provenance, "
          f"{without} without (recombination products, or written before Aug 2026)")
    if not rows:
        return 0

    anchor = reference if reference is not None else min(r[0] for r in rows)
    label = "reference" if reference is not None else "deepest structure"
    rel = [(r[0] - anchor) * HARTREE_TO_KJ for r in rows]
    times = [r[6] for r in rows]

    print(f"\nDeposit time [fs] -- when along its trajectory a structure was laid down")
    print(f"  all {len(rows)}: median {median(times):.0f}, quartiles "
          f"{quantile(times, 0.25):.0f} / {quantile(times, 0.75):.0f}, max {max(times)}")
    order = sorted(range(len(rows)), key=lambda i: rel[i])
    for k in (10, 50, 200):
        if k > len(order):
            continue
        picked = [times[i] for i in order[:k]]
        print(f"  {k:>4} deepest: median {median(picked):.0f}")
    print("  (a falling median means the useful structures appear EARLY -- the long "
          "trajectory then buys coverage, not depth)")

    print(f"\nYield per walker (one MD run per seed), energies relative to the {label}")
    per_walker = defaultdict(list)
    for row, r in zip(rows, rel):
        per_walker[row[4]].append(r)
    for walker in sorted(per_walker):
        values = per_walker[walker]
        print(f"  w{walker:>3}: {len(values):>5} structure(s), deepest {min(values):+8.1f} kJ/mol, "
              f"median {median(values):+8.1f}")

    print(f"\nYield per stage and repetition")
    per_stage = defaultdict(list)
    for row, r in zip(rows, rel):
        per_stage[(row[1], row[2], row[3])].append(r)
    for key in sorted(per_stage):
        values = per_stage[key]
        print(f"  {key[0]}_T{key[1]}K r{key[2]}: {len(values):>5} structure(s), "
              f"deepest {min(values):+8.1f} kJ/mol")
    return 0


if __name__ == "__main__":
    sys.exit(main())
