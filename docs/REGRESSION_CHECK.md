# Regression check against the reference sets (old binary vs new binary)

> 🤖 AI-generated, machine-tested. Numbers below are measurements from Sep 17/18, 2026.

`scripts/gmtkn55_compare.py` and `scripts/mor41_validation.py` compare curcuma
against **xtb**, which answers "how good is the method". This document is about
the other question, the one every change has to answer first: **did my change
move any energy at all, and where?** For that, the reference is not xtb but
curcuma itself, built from the commit the work started at. That also works on a
machine without xtb installed.

## How

```bash
python scripts/fetch_testset.py fetch mor41        # once per machine
python scripts/fetch_testset.py fetch gmtkn55

git worktree add /tmp/ref <commit-before-the-work>
cd /tmp/ref && mkdir release && cd release
cmake .. -DCMAKE_BUILD_TYPE=Release && make -j16 curcuma

python scripts/refset_regression.py --set mor41   --method gfnff \
    --old /tmp/ref/release/curcuma --new release/curcuma
python scripts/refset_regression.py --set gmtkn55 --method gfn2 ...
```

`refset_regression.py` runs every structure through both binaries single-threaded,
in its own scratch directory (every GMTKN55 structure file is called `struc.xyz`
and GFN-FF caches its topology next to that name - Known Issue #28), applies the
`.CHRG`/`.UHF` of each directory, and reports how many structures differ at all,
the MAD and the worst offenders.

## What the outcome means

- **0 differences** is the expected result for a change that is meant to be
  numerically neutral - a threading change with disjoint writes, a new code path
  that computes the same terms, a pure performance rewrite.
- **Last-digit differences** (1e-6 kcal/mol and below) are reassociated sums or a
  different BLAS routine, and are fine if intended.
- **Tens of kcal/mol on a few structures** needs a look at `scan_convergence.py`
  before drawing conclusions: a non-converged SCF wanders on any perturbation and
  says nothing about the change. GMTKN55/gfn2 has exactly three such structures,
  the single-atom cations `G21IP/b+`, `be+` and `c+`, whose SCF oscillates between
  about -0.79 and -0.93 Eh for all 151 iterations.

## Measured, Sep 2026

Reference `f36d6c11`, checked against `e4ba43ff` - i.e. across the FP64/FP32
reduction split, `-threads` driving the eigensolve, the three threaded
memory-bound SCF loops, the implicit GFN-FF Coulomb path, and the FP32
stagnation guard:

| set | method | structures | differ | MAD | max |
|---|---|---:|---:|---:|---:|
| MOR41 | gfnff | 95 | 0 | 0 | 0 |
| MOR41 | gfn1 | 95 | 0 | 0 | 0 |
| MOR41 | gfn2 | 95 | 0 | 0 | 0 |
| GMTKN55 | gfnff | 2462 | 0 | 0 | 0 |
| GMTKN55 | gfn1 | 2462 | 1 | 2.5e-09 | 6.3e-06 |
| GMTKN55 | gfn2 | 2462 | 2 | 2.2e-02 | 52.1 |

The gfn1 outlier is `CARBHB12/3CL_B` in its last printed digit. Both gfn2
outliers are non-converged single-atom cations (above), verified to be
non-converged in the OLD binary as well. The GFN-FF rows are the strongest
statement in the table, because those runs take a genuinely different code path
(implicit Coulomb pairs) and still reproduce every energy bit for bit.

## Traps this caught (do not repeat them)

- **Rebuilding while a comparison runs.** The first GMTKN55 pass of this session
  overlapped with a `make`, so part of it used one binary and part another. Any
  number from such a run is worthless; the comparison was repeated with a stable
  binary.
- **Reading a subset from a sample.** Characterising a subset from three
  structures once reported `MB16-43` an order of magnitude better than it was
  (Known Issue #23). Run the whole set.
- **Believing a cached energy.** `gmtkn55_compare.py` caches in
  `_run/energies.json`; `refset_regression.py` deliberately has no cache.
