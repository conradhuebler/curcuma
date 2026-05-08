# WP6 — Coulomb Cutoff / Cell-List (P4d) — Plan-only

**Status:** ⛔ Blockiert durch G2c — nur dokumentieren, nicht implementieren
**Aufwand:** 1 Tag (wenn G2c freigegeben)
**Erwarteter Nutzen:** Coulomb 5–20× bei N > 2000

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

## Akzeptanzkriterien (für spätere Umsetzung)

1. G2c aktiv und validiert (siehe Roadmap-Status).
2. ⚙️ `mixture.xyz` Coulomb-Wall von 640 ms → < 100 ms.
3. ⚙️ Energie-Diff CPU vs. GPU < 1 nEh.
4. ⚙️ Gradient-Diff vs. numerisch < 1e-5 Hartree/Bohr.
5. ⚙️ MD-Drift über 1000 Schritte (Polymer N=1410, NVE, kein Thermostat) < 1 mEh.
6. ⚙️ Energie-Diff vs. XTB-Fortran-Referenz < 50 µEh über fünf Test-Systeme.

## Risiken

- **Hoch.** Hellmann-Feynman-Verletzung bei einseitigem Cutoff (verifiziert Mai 2026). Der gesamte Cutoff-Pfad ist nicht-Fortran und braucht eigene Referenz-Validierung gegen XTB.
- **Energie-Surface-Verschiebung:** Selbst bei korrekter HF-Konsistenz weicht das Surface vom Fortran-Original ab. Ab welchem Cutoff der Effekt vernachlässigbar ist, muss gemessen werden.
- **Bottleneck-Verschiebung:** EEQ-Cholesky ist O(N³/6). Bei N > 5000 dominiert Cholesky, nicht Coulomb. WP6 löst nur einen Teil — siehe Hinweis in der Roadmap (G2c, "nicht der Bottleneck").

## Vorbedingung

- WP1 abgeschlossen (Threading klar)
- G2c aktiviert und validiert (separater Task, nicht Teil von WP6)
- WP3, WP4, WP5 vorzugsweise abgeschlossen, damit der Coulomb-Anteil neu gemessen werden kann (Anteilsverschiebung möglich)

## Cross-cutting-Regel

Cutoff ist eine **Approximation** im Sinne der README-Regel. Default bleibt `coulomb_distance_cutoff_bohr = 0.0` (= aus, voller O(N²)-Pfad bit-identisch zu pre-WP6). Aktivierung nur über expliziten ConfigManager-Schalter durch den Operator. Der voll-Coulomb-Pfad bleibt im Code, nicht entfernen.
