# Benchmark Test-Set Retrieval (on demand)

Grimme-group benchmark sets used to validate native GFN1/GFN2/GFN-FF are large,
carry their own citation terms, and are not committed to git. `scripts/fetch_testset.py`
downloads them into the exact layout the existing validation harnesses already
expect, and `scripts/testset_perf.py` benchmarks curcuma's CPU/threading/GPU
performance on whatever set is fetched.

## Fetching a test set

```bash
python scripts/fetch_testset.py list                 # registered sets + local status
python scripts/fetch_testset.py fetch mor41           # download + extract
python scripts/fetch_testset.py fetch gmtkn55 s30l    # multiple at once
python scripts/fetch_testset.py fetch mor41 --force   # re-fetch even if present
```

| name | source | destination | automatable |
|---|---|---|---|
| `mor41` | [chemie.uni-bonn.de/grimme](https://www.chemie.uni-bonn.de/grimme/de/software/mor41) tar.gz | `test_cases/MOR41-testset/<name>/mol.xyz` | yes |
| `gmtkn55` | [github.com/grimme-lab/GMTKN55](https://github.com/grimme-lab/GMTKN55) (shallow clone) | `test_cases/GMTKN55-testset/` | yes |
| `s30l` | ACS Supporting Information (10.1021/acs.jctc.5b00296), blocks automated fetches | `test_cases/s30l_test_set/<1..30>/{A,B,AB}/coord` + `.CHRG` + `reference_s30l` | no - manual, instructions printed by the script |

Fetched structures are gitignored (only the small hand-curated driver files -
`reactions.dat`, `reference_s30l` - are tracked); each destination gets a
`PROVENANCE.txt` with source URL, fetch timestamp and citation. Existing
harnesses need no changes: `scripts/mor41_validation.py` and `scripts/s30l_*.py`
already read from these exact paths.

```bash
python scripts/fetch_testset.py fetch mor41
python scripts/mor41_validation.py --only 5   # smoke test
```

`gmtkn55` is cloned as published (54 subsets on disk, 2462 structures,
geometries + reference data + the upstream `eval.py`).

## GMTKN55 vs xtb

`scripts/gmtkn55_compare.py` runs the same curcuma-vs-xtb reproduction check as
`mor41_validation.py`/`s30l_gfnff_compare.py`, but per-structure single points
(not weighted reaction energies - the `.res` files are `tmer2++` shell scripts,
a separate and much larger undertaking, not implemented here) across all 54
subsets:

```bash
python scripts/gmtkn55_compare.py --subset ACONF S22   # smoke test, fast
python scripts/gmtkn55_compare.py --method gfnff        # one method, all subsets
python scripts/gmtkn55_compare.py                       # full sweep (~30-60 min)
```

Both engines' binaries are located automatically (`$XTB_BIN` env var, else
`which xtb`, else a few common install paths - override with `--xtb PATH`).
Output: `test_cases/GMTKN55-testset/_run/gmtkn55_{results,summary}_<method>.csv|md`
(per-subset MAD/RMSD/max + global stats + outlier list), cached incrementally
in `_run/energies.json` so `--subset`/`--limit` reruns and a resumed full sweep
never recompute.

**Open-shell limitation (verified Sep 2026, see CLAUDE.md Known Issues #9)**:
native `-sp` has no working path to request UHF occupation for gfn1/gfn2 -
`-spin N` sets inert `Molecule` metadata that no `energy_calculators/` code
reads. `gmtkn55_compare.py` therefore skips every structure with a nonzero
`.UHF` file for gfn1/gfn2 (about a third of GMTKN55: RC21, RSE43, G21EA/G21IP,
parts of BH76, ...) rather than silently report a wrong closed-shell number;
gfnff is unaffected (no explicit open-shell term) and runs everything.

**Results (full sweep, Sep 2026)**: gfn2 MAD 0.000 / gfn1 MAD 0.047 / gfnff
MAD 19.98 kcal/mol vs xtb. The GFN-FF number is dominated by a found bug -
native GFN-FF returns exactly 0.0 Eh for any single free atom (missing EEQ
self-energy, see CLAUDE.md Known Issues #8); excluding those 107 structures
GFN-FF MAD drops to 3.54 kcal/mol. Full analysis, per-subset breakdown and
outlier categories: [docs/GMTKN55_VALIDATION.md](GMTKN55_VALIDATION.md).

## Performance benchmark

`scripts/testset_perf.py` times `curcuma -sp`/`-opt` on the largest structures
of a fetched set across a `-threads N` grid and every GPU backend whose plugin
(`release/libcurcuma_<backend>.so`) is actually present on the machine, so a
run never silently reports CPU numbers as GPU numbers.

```bash
python scripts/testset_perf.py                                    # MOR41, cpu + detected GPU plugins
python scripts/testset_perf.py --threads 1,2,4,8,16 --n-structures 5
python scripts/testset_perf.py --gpu cpu,cuda --opt --repeats 5    # heavier -opt workload
```

Output goes to `<root>/_run/perf_<timestamp>/{results.csv,summary.md}`
(gitignored): per-structure median wall time, threading speedup/parallel
efficiency relative to `threads=1`, and GPU-vs-best-CPU speedup per backend.
Numbers are machine-local and not committed anywhere; re-run on the target
machine for the paper's performance section rather than reusing a quoted
figure.

## Adding a new registry entry

Edit the `REGISTRY` dict in `scripts/fetch_testset.py` - `kind` is `tar`
(download + extract archive), `git` (shallow clone), or `manual` (print
instructions, e.g. for paywalled Supporting Information). `check` is a
callable that confirms the expected layout landed correctly.
