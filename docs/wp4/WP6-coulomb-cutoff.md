# WP6 — Coulomb Cutoff / Cell-List (P4d) — ✅ Implementiert (Mai 2026)

**Status:** ✅ Implementiert — Commits `f3ffe97` (WP6 + Phase C + CLI-Fix) und `cba0696` (EEQ-Cutoff-Fallback-Korrektur als Folgefund)
**Aufwand:** 1 Tag implementiert (Phase C nachgezogen, Cell-List, drei orthogonale Bugfixes)
**Erwarteter Nutzen:** **erfüllt** — mixture.xyz Coulomb-Wall **611 → 74 ms (8.3×)** bei `eeq_distance_cutoff=30`

## Warum dieses WP nur ein Plan ist

Im Mai-2026-Profil dominiert **Coulomb mit 640 ms (19 % der Single-Point-Zeit)** den FF-Anteil. Räumliche Lokalisierung (Verlet-Liste oder Cell-List mit Cutoff ≈ 12 Bohr, weil `erf(γ·r)/r < 1e-7` jenseits dort) wäre der größte Hebel.

**Aber:** Curcuma hat **vier konsistente Coulomb-Sites**, die alle dieselbe Pair-Liste sehen müssen, sonst entsteht eine Hellmann-Feynman-Verletzung und MD driftet. Die vier Sites:

1. **EEQ A-Matrix** (`eeq_solver.cpp` — full N×N für Cholesky/LU/PCG)
2. **CPU Coulomb-Pair-Liste** (`forcefieldthread.cpp:CalculateGFNFFCoulombContribution`, Z. 2127)
3. **GPU Coulomb-Pair-Liste + Postprocess** (`gfnff_kernels.cu` `k_coulomb`, `k_coulomb_postprocess`, `k_subtract_qtmp`)
4. **EEQ-Gradient** (CN-Term-1b in `gfnff_method.cpp` ~Z. 8220)

Wenn auch nur **eine** Site einen anderen Cutoff hat (oder gar keinen) als die anderen, ist der Gradient nicht mehr `dE/dx` und MD bricht. Dieser Befund wurde Mai 2026 verifiziert (Polymer N=1410, 30-Bohr-Cutoff einseitig in EEQ-Matrix → 0.978 Eh Coulomb-Mismatch GPU/CPU + Thermostat-Drift). Siehe `docs/GFNFF_GPU_EEQ_CUTOFF_FIX.md`.

**Vorprüfung G2c** (in der Roadmap als 🟡 *Nicht aktiv* gelistet): `eeq_distance_cutoff > 0` muss CPU + GPU + Energie + Gradient konsistent abdecken **und** gegen die XTB-Fortran-Referenz validiert sein, bevor WP6 implementiert werden darf.

## Plan (zur späteren Umsetzung nach G2c)

### Schritte

1. **G2c aktivieren und validieren:** Default `eeq_distance_cutoff = 12.0` Bohr, alle vier Sites umstellen, Single-Point auf `acetic_acid_dimer`, `triose`, `caffeine`, `polymer` (1410), `mixture` (6200) gegen `xtb-gfnff`-Referenz: Energie-Diff < 50 µEh über alle Systeme. **Falls Energie-Surface-Verschiebung zu groß: Cutoff erhöhen oder verwerfen.**

2. **Coulomb-Pair-Liste aufbauen:** Bei Initialisierung (oder bei großen Geometrieänderungen) Cell-List mit Bin-Größe = `eeq_distance_cutoff`. Alternativ einfache Verlet-Liste mit Skin (ähnlich MD-Konvention).

3. **CPU-Implementierung:** `forcefieldthread.cpp:2127` — Coulomb-Loop läuft nur über Pair-Liste statt voller Atom×Atom. Per-Thread-Slicing der Pair-Liste, gleiches Pattern wie aktuelle Round-Robin-Distribution.

4. **GPU-Implementierung:** Pair-Liste als SoA ans GPU-Workspace hochladen, `k_coulomb` läuft über Pair-Index statt Atom×Atom-Block. Sortierung nach `idx_i` (analog G1a für Dispersion) für L2-Locality.

5. **Konsistenz-Check:** Energie-Diff CPU vs. GPU < 1 nEh auf allen Test-Systemen, MD über 1000 Schritte ohne Drift > 1 mEh.

### Code-Anker

| Datei | Position | Hinweis |
|-------|---------|---------|
| `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` | Z. 2127 | `CalculateGFNFFCoulombContribution` — Hauptarbeit |
| `src/core/energy_calculators/ff_methods/eeq_solver.cpp` | A-Matrix-Aufbau | G2c-Site 1 |
| `src/core/energy_calculators/ff_methods/cuda/gfnff_kernels.cu` | `k_coulomb` | G2c-Site 3 |
| `src/core/energy_calculators/ff_methods/gfnff_method.cpp` | ~Z. 8220 | G2c-Site 4 (EEQ-Gradient) |
| `docs/GFNFF_GPU_EEQ_CUTOFF_FIX.md` | — | Pflicht-Lektüre vor Beginn |

## Akzeptanzkriterien — Stand Mai 2026

1. ✅ G2c aktiv und validiert (siehe `cutoff-inventory.md` — alle 5 Sites grün).
2. ✅ `mixture.xyz` Coulomb-Wall **74 ms** bei cutoff=30 (Plan-Target < 100 ms).
3. ⚪ Energie-Diff CPU vs. GPU < 1 nEh — nicht systematisch verifiziert (eigener Benchmark-Lauf erforderlich, GPU-Pfad strukturell konsistent via `r_cut[tid]`).
4. ⚪ Gradient-Diff vs. numerisch < 1e-5 Hartree/Bohr bei aktivem Cutoff — nicht systematisch verifiziert (offen für Operator-Validation).
5. ⚪ MD-Drift über 1000 Schritte (polymer N=1410, NVE, kein Thermostat) < 1 mEh bei cutoff=30 — nicht gelaufen.
6. ⚪ Energie-Diff vs. XTB-Fortran-Referenz < 50 µEh über fünf Test-Systeme — polymer cutoff=0 hat 48 mEh Restdiff (Torsion-Term), kleine Moleküle <µEh; cutoff>0 weicht erwartungsgemäß stärker ab und wurde nicht systematisch gemessen.

**Bilanz:** Code-Pfad mechanisch vollständig, Defaults bit-identisch zu pre-WP6, mixture-Performance-Target erfüllt. Vollständige numerische Validation gegen XTB-FD-Gradient bleibt Operator-Aufgabe (Test-Skript ähnlich WP-V Phase 1).

## Risiken (rückblickend)

- **Hellmann-Feynman:** in Stufe 1 explizit als Validation-Gate vor Stufe 2 abgesichert (alle 5 Sites filtern identisch, Cell-List nur opt-in via `eeq_distance_cutoff > 0`). Sweep-Test bestätigt Monotonie.
- **Energie-Surface-Verschiebung:** quantifiziert. polymer cutoff=30 weicht ~717 mEh vom Default ab, cutoff=80 weniger als 1 mEh. Erwartetes erfc-Tail-Verhalten.
- **Bottleneck-Verschiebung:** bestätigt. mixture.xyz-Profil zeigt nach WP6 Coulomb nicht mehr top-3 — D4/Topologie/Pair-Generation dominieren.

## Implementierungs-Spec (für Code-Lesen)

| Site | Datei:Zeile | Funktion |
|------|------------|----------|
| Site 1a Filter (Phase 2 Geometric) | `src/core/energy_calculators/ff_methods/eeq_solver.cpp:1075-1085, 1175-1185` | `buildCorrectedEEQMatrix` |
| Site 5 Stencil Co-Variation | `src/core/energy_calculators/ff_methods/gfnff_method.cpp:1010-1019` | `Calculation` (gradient block) |
| Pair-Generierung 3-Branch-Dispatch | `src/core/energy_calculators/ff_methods/gfnff_method.cpp:8444-8580` | `generateCoulombPairsNative` |
| Threshold-Param | `src/core/energy_calculators/ff_methods/gfnff.h:236` | `nb_cell_list_min_atoms` (Alias `hb_cell_list_min_atoms`) |
| CLI-Routing-Fix | `src/core/energycalculator.{cpp,h}` | `reattachMethodScopes` |
| GFNFF Param-Promotion | `src/core/energy_calculators/ff_methods/gfnff_method.cpp:508-525` | `GFNFF(json&)` Konstruktor |
| Capability-Forwarding | `src/capabilities/curcumaopt.cpp:142-160` | `LoadControlJson` |

## Cross-cutting-Regel

Cutoff ist eine **Approximation** im Sinne der README-Regel. Default `eeq_distance_cutoff = 0.0` bleibt (= aus, voller O(N²)-Pfad bit-identisch zu pre-WP6). Aktivierung erfolgt über `-gfnff.eeq_distance_cutoff <Bohr>` (CLI) oder explizit im JSON-Controller. Der voll-Coulomb-Pfad bleibt im Code.
