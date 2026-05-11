# WP-S2: MD-Diagnostics — JSONL-Streaming pro Step

**Kategorie**: Klein
**Aufwand**: ~1 Tag (~100 LoC)
**Wirkung**: post-mortem-Analyse von Charges/CN/Energien-Trajektorie; Pflicht für Static-Mode-Validation
**Abhängigkeiten**: WP-S1 (gemeinsam genutzt für Drift-Validation), unabhängig nutzbar
**Status**: 🤖 Geplant

---

## Motivation

Für Static-Mode-Validation (WP-S1) und allgemeine GFN-FF MD-Forensik fehlt eine
streaming-fähige Pro-Step-Aufzeichnung interner Größen. `simplemd` schreibt
heute nur `.xyz`-Trajektorie + globale Energie — Charge-Drift, CN-Drift,
Energie-Dekomposition pro Term sind unsichtbar.

Lösung: parallel zur Trajektorie eine JSONL-Datei (`<basename>.diag.jsonl`),
eine Zeile pro Snapshot, streaming-append-fähig.

---

## Output-Format

`<run>.diag.jsonl` — eine JSON-Object-Zeile pro Snapshot:

```json
{
  "step": 1000,
  "time_ps": 1.0,
  "energy": {
    "total": -45.123456,
    "bonds": -12.34,
    "angles": -5.67,
    "torsions": -2.10,
    "inversions": -0.05,
    "dispersion": -1.23,
    "repulsion_bonded": 0.45,
    "repulsion_nonbonded": 0.87,
    "coulomb": -22.31,
    "hbond": -2.05,
    "xbond": 0.0,
    "storsion": 0.0,
    "batm": -0.001,
    "atm": -0.0001
  },
  "charges":      [-0.231, 0.114, 0.117, ...],
  "cn":           [3.98, 1.02, 1.01, ...],
  "gradient_norm":[2.3e-4, 1.1e-4, ...],
  "hb_count": 3,
  "xb_count": 0
}
```

Größe-Schätzung:
- Polymer (N=472), ~22 KB / Zeile bei `charges+cn+grad_norm+energies` (full)
- 1000 Steps → ~22 MB
- 10000 Steps → 220 MB (vertretbar; optional HDF5-Backend)

---

## Implementierungspunkte

### 1. Neue PARAMs (`simplemd.h` oder `gfnff.h`)

```cpp
PARAM(md_diagnostics_freq, Int, 0,
      "Frequency (in MD steps) to write per-step diagnostics to <basename>.diag.jsonl. "
      "0 = disabled. Records full energy decomposition, charges, CN, gradient norms, "
      "HB/XB counts per snapshot. Useful for static-mode validation and forensic analysis.",
      "Output", {})

PARAM(md_diagnostics_fields, String, "energies,charges,cn,grad_norm,hb_xb",
      "Comma-separated list of fields to dump. Options: energies, charges, cn, "
      "grad_norm, hb_xb, fragments. Use 'energies' alone for minimal overhead.",
      "Output", {})

PARAM(md_diagnostics_format, String, "jsonl",
      "Output format: 'jsonl' (one JSON object per line, streaming-friendly) or "
      "'hdf5' (compressed binary, faster post-processing — requires HDF5 build).",
      "Output", {})
```

### 2. Diagnostics-Writer (`src/capabilities/md_diagnostics.h/cpp` neu)

Minimaler Wrapper um `std::ofstream` mit append-Mode:

```cpp
class MDDiagnosticsWriter {
public:
    MDDiagnosticsWriter(const std::string& basename,
                        const std::set<std::string>& fields);
    void writeSnapshot(int step, double time_ps,
                       const json& energy_decomp,
                       const Vector& charges,
                       const Vector& cn,
                       const Matrix& gradient,
                       int hb_count, int xb_count);
private:
    std::ofstream m_out;
    std::set<std::string> m_fields;
};
```

### 3. Hook in `simplemd.cpp`

Nach jeder MD-Step-Output-Phase (parallel zu `WriteGeometry()`):

```cpp
if (m_diag_writer && step % m_diag_freq == 0) {
    // Energy decomposition vom letzten EnergyCalculator-Call abrufen
    json energy_decomp = m_method->getEnergyDecomposition();   // neuer ComputationalMethod-Hook
    m_diag_writer->writeSnapshot(step, time_ps, energy_decomp,
                                 m_method->getCharges(),
                                 m_method->getCN(),         // neuer Hook
                                 m_method->getGradient(),
                                 m_method->getHBCount(),
                                 m_method->getXBCount());
}
```

### 4. ComputationalMethod-Hooks (neu)

```cpp
// computational_method.h
virtual json getEnergyDecomposition() const { return {}; }
virtual Vector getCN() const                { return {}; }
virtual int getHBCount() const              { return 0; }
virtual int getXBCount() const              { return 0; }
```

Standardimplementierung leer; nur `GFNFFMethod` überschreibt.

### 5. Post-Processing-Helper (`scripts/analyse_diag.py` neu)

Optional, aber empfohlen:

```python
# scripts/analyse_diag.py — quick line-by-line stream loader
import json, sys
for line in open(sys.argv[1]):
    rec = json.loads(line)
    print(rec["step"], rec["energy"]["total"], max(rec["charges"]))
```

---

## Tests

`test_cases/cli/simplemd/11_diagnostics_dump`:
- 50-step NVE polymer
- `md_diagnostics_freq=5` → erwarte 10 JSONL-Zeilen
- Validation: jede Zeile JSON-parsbar, `charges.length == atom_count`, `energy.total`
  monoton (NVE)
- Reproducibility: zwei Läufe → byte-identische JSONL-Dateien

---

## Out of Scope

- Live-Plotting / Streaming-Visualisierung: separates Folge-WP
- HDF5-Backend: optional, kann später als Drop-in über `md_diagnostics_format=hdf5`
- Trajectory-Replay über JSONL: separates WP, evtl. via `curcuma -replay`

---

## Acceptance-Checkliste

- [ ] 3 PARAMs registriert, in `make GenerateParams` ohne Warning
- [ ] `MDDiagnosticsWriter` + Append-Logik in `simplemd.cpp`
- [ ] 4 neue Hooks in `ComputationalMethod` + Override in `GFNFFMethod`
- [ ] CLI-Test `11_diagnostics_dump` grün
- [ ] `scripts/analyse_diag.py` minimaler Reader committed
- [ ] Dokumentation in `docs/MD_DIAGNOSTICS.md` (neue Datei) + Link aus
      `docs/GFNFF_STATUS.md`
