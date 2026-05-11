# WP-P3: GFN-FF Bond-Term Microbench + Threading-Audit

**Kategorie**: Mittel
**Aufwand**: ~1.5 Tage (~150 LoC C++ Microbench + ~50 LoC Threading-Audit + Optimierung falls Hypothese bestätigt)
**Wirkung**: Klärt das 75-µs/Bond-Mai-2026-Mysterium quantitativ; identifiziert reale CPU-Hot-Spots in CN/EEQ-serial-Block; nutzt Erkenntnisse aus WP-P1/P2 als Eingangsdaten
**Abhängigkeiten**: WP-P1 (Timing-Hooks zum quantitativen Vergleich); WP-P2 (Sweep-Daten für Cross-System-Verifikation)
**Status**: 🤖 Geplant
**Priorität**: zuletzt — die Optimierungs-Entscheidungen ergeben sich aus WP-P2-Daten

---

## Motivation

Die Performance-Roadmap (`docs/GFNFF-PERFORMANCE-ROADMAP.md:198`) listet
zwei "ungeklärte" CPU-Hotspots:

1. **P4c — Bond-Term**: `359 ms / 4800 Bonds = 75 µs/Bond` bei mixture.xyz
   ist "~100× zu langsam für eine simple `exp((r-r0)·α)`-Auswertung".
   Hypothese: Cache-Misses durch Pointer-Indirection oder Eigen-Matrix-Ops
   pro Bond statt Batch. **Nie isoliert profiled.**
2. **P4b — Threading-Audit**: "CN + EEQ (serial CPU)" macht 53 % der
   mixture.xyz-Wall-Clock aus (1766 ms / 3341 ms). Unsere Static-Mode-Tests
   zeigten 0 % Speedup für `static_cn` allein — Indiz dass CN nicht der
   Bottleneck ist wie angenommen, oder dass die Skalierung kaputt ist.

Beide Lücken sind essentielle Eingangsdaten für jede ernsthafte
CPU-Optimierungs-Entscheidung. WP-P3 klärt sie mit isolierten Microbenchmarks.

---

## Implementierungspunkte

### Teil A — Bond-Term Microbench

Neuer Unit-Test `test_cases/test_gfnff_timing.cpp`:

```cpp
#include "src/core/energy_calculators/ff_methods/forcefieldthread.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "core/test_molecule_registry.h"
#include <chrono>
#include <iostream>

static void benchmarkBondOnly(const std::string& molecule_name, int repeats) {
    using namespace TestMolecules;
    auto mol = TestMoleculeRegistry::createMolecule(molecule_name, false);

    GFNFF gfnff(default_gfnff_params());
    gfnff.InitialiseMolecule(mol);

    // Get bond list from prepared parameters
    const auto& bond_list = gfnff.getForceField()->getGFNFFBonds();
    int n_bonds = bond_list.size();

    // Repeatedly call only the bond contribution kernel
    ForceFieldThread thread(0, gfnff.getForceField());
    thread.setMethod(3);  // GFN-FF
    thread.setGFNFFBonds(bond_list);

    auto t0 = std::chrono::high_resolution_clock::now();
    for (int r = 0; r < repeats; ++r) {
        thread.resetEnergy();
        thread.CalculateGFNFFBondContribution();
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double us_per_bond = (total_ms * 1000.0) / (repeats * n_bonds);

    std::cout << "Bond microbench " << molecule_name
              << ": n_bonds=" << n_bonds
              << " repeats=" << repeats
              << " total=" << total_ms << " ms"
              << " us/bond=" << us_per_bond
              << "\n";

    // Roadmap hypothesis: 75 us/bond. Realistic for simple exp: < 5 us/bond.
    if (us_per_bond > 10.0) {
        std::cout << "  [HOT] >10 us/bond — confirms P4c hypothesis\n";
    }
}

int main() {
    benchmarkBondOnly("CH4",        100000);   // 4 bonds, very cold
    benchmarkBondOnly("CH3OCH3",    100000);   // 8 bonds, small
    benchmarkBondOnly("triose",     10000);    // ~70 bonds, medium
    benchmarkBondOnly("monosaccharide", 10000); // ~30 bonds
    return 0;
}
```

CTest-Eintrag:
```cmake
add_executable(test_gfnff_timing test_gfnff_timing.cpp)
target_link_libraries(test_gfnff_timing curcuma_core)
add_test(NAME test_gfnff_timing COMMAND test_gfnff_timing)
set_tests_properties(test_gfnff_timing PROPERTIES LABELS "gfnff;profile" TIMEOUT 60)
```

Test druckt nur Ergebnisse — kein Pass/Fail-Check (außer "läuft ohne
Crash"). Output wird manuell gegen die Mai-2026-Hypothese verglichen.

### Teil B — `perf record` Wrapper (optional)

Falls Microbench die 75-µs-Hypothese bestätigt: Skript für gezielten Hotspot-Locate.

`scripts/perf_bond_hotspot.sh`:

```bash
#!/usr/bin/env bash
# Wraps `perf record` around test_gfnff_timing to find the Bond-loop hotspot.
# Requires: linux-tools-common (perf), root or unprivileged perf_event_paranoid<=1.
set -e

BUILD=${1:-release}
cd "$BUILD"
perf record -g --call-graph dwarf -F 999 -o bond.perf \
    ./bin/test_gfnff_timing
perf report -i bond.perf --stdio | head -50
echo "Full report saved to bond.perf — open with 'perf report -i bond.perf'"
```

Nicht-CI, nur Entwickler-Tool.

### Teil C — Threading-Audit für CN/EEQ-Serial-Block

Erweiterung des WP-P1-Timing-Hooks: zusätzliches Feld `prep_cn` wird
**pro Thread-Count** gemessen. Sweep-Skript aus WP-P2 mit `--threads 1,4,8`-Variation:

```bash
python3 scripts/wp_profile_sweep.py --systems polymer \
        --paths cpu --threads 1,4,8 --max-time 100 \
        --output /tmp/sweep_threads
```

Ergebnis-Tabelle zeigt direkt, ob `prep_cn` linear skaliert (ideal),
sublinear (OpenMP-Overhead dominiert) oder gar nicht (Bottleneck nicht
parallelisiert).

Erweiterung des Sweep-Skripts (~30 LoC), Konsumiert von WP-P2.

### Teil D — Cross-Reference: Microbench-Ergebnis vs MD-Sweep

Neue Sektion im `docs/GFNFF_PROFILE_RESULTS_2026-05.md` (gehört eigentlich
zu WP-P2, aber Inhalt kommt aus WP-P3):

```markdown
## Bond-Hotspot-Klärung

Microbench-Ergebnis (test_gfnff_timing, repeats=10000):
- CH4 (4 bonds, isolated)   : X.X µs/bond
- triose (70 bonds)         : Y.Y µs/bond
- monosaccharide (30 bonds) : Z.Z µs/bond

Vergleich Mai-2026-Hypothese: 75 µs/bond (mixture, 4800 bonds).

Interpretation:
- Falls X.X << 75 µs: 75-µs-Wert war Cache-Miss-dominiert bei großem N.
  → Bond-Loop ist OK, Hotspot kommt aus Memory-Pressure beim großen System.
- Falls X.X ≈ 75 µs: Hotspot ist inhärent im Bond-Kernel.
  → Optimierung in `CalculateGFNFFBondContribution()` lohnt.
```

### Teil E — Falls Hypothese bestätigt: Optimierungs-Patch (optional)

Wenn der Microbench die 75-µs-Hypothese bestätigt, sind die folgenden
Patches die direkten Kandidaten:

1. **Pre-fetch CN-Werte** für alle Bonds in einen kleinen Stack-Vector
   einmal pro execute()-Call (vermeidet `m_cn(bond.i)`-Indirection pro Bond).
2. **Batch-Geometrie-Read**: statt `m_geometry.row(i)` per Bond, einmal
   alle benötigten Coords in ein lokales `std::vector<Vector3d>` kopieren.
3. **Bond-Liste sortieren** nach Atom-Index für bessere Cache-Lokalität.

Patch ist Teil dieses WP, falls Microbench-Daten ihn rechtfertigen. Sonst
verschoben in eigenes WP.

---

## Kritische Files

| Datei | Änderung | LoC |
|-------|----------|-----|
| `test_cases/test_gfnff_timing.cpp` (neu) | Microbench-Driver mit 4 Molekül-Größen | ~80 |
| `test_cases/CMakeLists.txt` | `add_executable` + `add_test` | ~6 |
| `scripts/perf_bond_hotspot.sh` (neu, optional) | `perf record`-Wrapper | ~15 |
| `scripts/wp_profile_sweep.py` (Erweiterung) | `--threads N,M,K`-Option | ~30 |
| `docs/GFNFF_PROFILE_RESULTS_2026-05.md` | Sektion "Bond-Hotspot-Klärung" + Threading-Skalierung | ~50 |
| `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` (falls Optimierung greift) | Pre-fetch CN / batch Coords | ~30 |

Build-Verify: `cd release && make -j4 && ctest -R test_gfnff_timing`.

---

## Verifikation

### a) Microbench läuft
```bash
cd /home/conrad/src/curcuma/release
MKL_ROOT=/opt/intel/oneapi/mkl/latest make -j4 test_gfnff_timing
ctest -R test_gfnff_timing --output-on-failure
```
Akzeptanz: Test läuft ohne Crash, druckt 4 Bond-Timings.

### b) Threading-Skalierung
```bash
python3 scripts/wp_profile_sweep.py --systems polymer --paths cpu \
        --threads 1,4,8 --max-time 100 --output /tmp/sweep_threads
grep -E "prep_cn|prep_eeq_solve" /tmp/sweep_threads/wp_profile_summary.md
```
Akzeptanz: Tabelle zeigt Werte für T=1, T=4, T=8. Speedup-Faktor `prep_cn`
T=1→4 sollte ≥ 2.0 sein, sonst Hinweis auf Threading-Bug.

### c) Cross-Reference im Report
```bash
grep -A 10 "Bond-Hotspot-Klärung" docs/GFNFF_PROFILE_RESULTS_2026-05.md
```
Akzeptanz: Sektion enthält Microbench-Werte und Interpretation.

### d) Regression-Check
```bash
ctest --output-on-failure
```
Akzeptanz: neuer Test `test_gfnff_timing` grün, andere unverändert.

---

## Erwartete Ergebnisse / Entscheidungs-Baum

| Microbench-Ergebnis | Interpretation | Folge-Aktion |
|---------------------|----------------|--------------|
| ≤ 5 µs/bond | Bond-Loop selbst ist OK. Mai-2026-Profil war von Memory-Pressure bei N=6200 dominiert. | Kein Bond-Patch nötig. Hotspot bei mixture liegt eher in Cache-/Memory-Bandbreite. |
| 5-30 µs/bond | Bond-Loop hat moderaten Overhead durch Indirection. | Pre-fetch CN + batch Geo lohnt; eigenes WP für Patch+Mess. |
| > 30 µs/bond | Bond-Loop ist klarer Hotspot, Optimierung dringend. | Patch in diesem WP integriert. |

| Threading-Skalierung (T=1→4) | Interpretation | Folge-Aktion |
|------------------------------|----------------|--------------|
| ≥ 3× | OpenMP wirkt. Keine Aktion. | — |
| 1.5×-3× | Skalierung suboptimal, evtl. OpenMP-Overhead. | Profilieren, kein eigenes WP. |
| < 1.5× | Serial-Bottleneck. P4b-Hypothese bestätigt. | Threading-Audit-WP starten. |

---

## Out of Scope

- **SoA-Refactoring** der CPU-Geometrie → eigenes WP nach WP-P3
- **SIMD-Vektorisierung der CN-erf()-Schleife** (P4e) → eigenes WP
- **D4-Gaussian-Weights-Vektorisierung** (P4f) → eigenes WP
- **GPU-Pfad-Bond-Term** → bereits in WP-P2 als Sweep-Eintrag, hier nicht extra

---

## Risiken

| Risiko | Mitigation |
|--------|-----------|
| Microbench-Setup zu artifiziell (nur Bond-Loop, kein Memory-Pressure) | Bewusst — Ziel ist Isolation. Mixture-Echtwert aus WP-P2-Sweep liefert das Cross-Reference |
| `perf record` braucht Root oder Setup von perf_event_paranoid | Skript ist optional, kein CI-Block |
| Threading-Audit erkennt Sub-Block-Skalierung nicht | WP-P1 hat `prep_cn` als atomares Feld; Sub-Block-Profile bleibt Out-of-Scope |
| Optimierung in Bond-Loop bricht Gradient-Test | bestehender `test_gfnff_gradients` läuft als Regression vor Merge |
| Microbench-Linkage gegen `forcefieldthread` braucht halb-private API | falls nötig, Helper-Methode `runBondsOnly()` als Public-Method in ForceFieldThread hinzufügen (~5 LoC) |
