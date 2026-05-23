# MD Diagnostics — Per-Step JSONL Dump

**Status**: ✅ Implemented (WP-S2, May 2026) — ⚙️ machine-tested, human production
testing pending.

`simplemd` can optionally append a streaming-friendly JSONL file
`<basename>.diag.jsonl` parallel to the XYZ trajectory. Each line is one
JSON object describing the state at one snapshot — useful for static-mode
validation, charge/CN drift forensics, and energy-decomposition post-mortems.

## Activation

Single bool PARAM `md_diagnostics` (default `false`). Frequency follows
`dump_frequency` (one record whenever `WriteGeometry()` writes a frame).

```bash
# JSON-config snippet
{
  "simplemd": {
    "method": "gfnff",
    "max_time": 1000.0,
    "dt": 1.0,
    "dump_frequency": 5,
    "md_diagnostics": true
  }
}

# Run
./curcuma -md mol.xyz -import_config config.json
# → produces mol.trj.xyz   (existing)
# → produces mol.diag.jsonl (new)
```

## Format

One JSON object per line, append-mode + flush after each snapshot
(crash-robust):

```json
{
  "step": 5,
  "time_fs": 5.0,
  "energy": {
    "Bond": -0.123, "Angle": -0.045, "Torsion": -0.001,
    "Inversion": 0.0, "Dispersion": -0.0034, "Coulomb": -4.76,
    "HBond": -0.012, "XBond": 0.0, "ATM": -1e-7, "BATM": 5e-6
  },
  "charges":       [-0.231, 0.114, 0.117, ...],
  "cn":            [3.98, 1.02, 1.01, ...],
  "gradient_norm": [2.3e-4, 1.1e-4, ...],
  "hb_count": 3,
  "xb_count": 0
}
```

All `charges` / `cn` / `gradient_norm` arrays have length `atom_count`.
Energy components in Hartree. `time_fs` in femtoseconds. Empty arrays
(non-GFN-FF methods) are valid.

## Size

- Polymer (N≈470): ~22 KB/snapshot full → 1000 snapshots ≈ 22 MB.
- Caffeine (N=24): ~1 KB/snapshot → 1000 snapshots ≈ 1 MB.

For long production runs, raise `dump_frequency` accordingly.

## Reader

`scripts/analyse_diag.py` is a minimal one-line-at-a-time loader showing
the pattern:

```bash
python3 scripts/analyse_diag.py mol.diag.jsonl
#     step    time_fs        E_total   max|q|  max(CN)  hb  xb
#        0       0.00      -4.972771   0.4697    3.499 358   0
#        5       5.00      -4.947431   0.4768    3.506 358   0
#       ...
```

Custom analysis: iterate `open(path)` line-by-line, `json.loads(line)` —
no need to load the entire file into memory.

## Method-Coverage

| Method | charges | cn | hb/xb | energy decomposition |
|--------|---------|----|-------|-----------------------|
| `gfnff` (native) | ✅ | ✅ | ✅ | full (10 terms) |
| `uff` / `uff-d3` | ✅ | empty | 0 | basic FF terms |
| `gfn2` / xtb / tblite | ✅ | empty | 0 | zeros (QM, no native decomp) |

Hooks (`getCN`, `getHBCount`, `getXBCount`) live in
`ComputationalMethod` with empty defaults; only
`GFNFFComputationalMethod` overrides them.

## Out of Scope (planned follow-ups)

- HDF5-backend for fast post-processing.
- `md_diagnostics_freq` independent of `dump_frequency`.
- `md_diagnostics_fields` selective field filter.
- Live-plotting hook.
- Trajectory-replay through JSONL.

## Implementation Pointers

- PARAM: `src/capabilities/simplemd.h` (within `BEGIN_PARAMETER_DEFINITION(simplemd)`).
- Writer: `src/capabilities/md_diagnostics.{h,cpp}`.
- SimpleMD hook: `src/capabilities/simplemd.cpp::step()`, inside the
  existing `m_step % m_dump == 0` block (alongside `WriteGeometry()`).
- ComputationalMethod hooks: `src/core/energy_calculators/computational_method.h`.
- GFN-FF override: `src/core/energy_calculators/qm_methods/gfnff_method.{h,cpp}`.
- EnergyCalculator forwarders: `src/core/energycalculator.{h,cpp}` —
  `CN()`, `HBCount()`, `XBCount()`.
