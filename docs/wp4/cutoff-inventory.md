# Cutoff-Parameter Inventar — GFN-FF + EEQ (Mai 2026)

Bestandsaufnahme aller Distanz-/Energie-/Konvergenz-Cutoff-Parameter im GFN-FF + EEQ Pfad. Stand: vor WP-C (G2c-Aktivierung). Datenquelle: Phase-1-Audit der WP-C-Plan-Phase.

## Zweck dieses Dokuments

Vor jeder Cutoff-bezogenen Optimierung (insb. WP6 Coulomb-Cell-List) muss klar sein, **welche Cutoffs es gibt und wie sie sich verhalten**. Heutige Praxis ist heterogen: Mischung aus ConfigManager-PARAMs, Struct-Defaults, Hardcoded-Konstanten — manche CPU/GPU-konsistent, manche nicht.

Wenn künftig ein neuer Cutoff dazukommt, sollte er sich an die unten beschriebene Hierarchie halten.

## Tabelle aller Cutoff-Parameter

| # | Parameter | Datei:Zeile | Term/Phase | Default (Bohr) | Quelle | CPU? | GPU? | Squared? | Status |
|---|-----------|-------------|------------|---------------:|--------|------|------|----------|--------|
| 1 | `cn_cutoff_bohr` | `gfnff.h:219` | CN-Nachbar-Liste | 6.0 | PARAM | ✅ Read | ⚠️ N/A (GPU CN nutzt anderen Weg) | linear | ✅ funktional |
| 2 | `cn_accuracy` | `gfnff.h:220` | CN-Schwelle (Fallback) | 1.0 | PARAM | ✅ Read | ⚠️ N/A | berechnet | ✅ funktional |
| 3 | **CN-Derivate-Threshold (hardcoded)** | `gfnff_method.cpp:3436, 6762, 6889` | CN-Derivate | **1600.0 Bohr² (= 40 Bohr)** | Hardcode (Konstante) | ✅ Used | ⚠️ N/A | **squared** | ❌ **Inkohärenz**: ignoriert `cn_cutoff_bohr` |
| 4 | `eeq_distance_cutoff` (gfnff) | `gfnff.h:229` | EEQ A-Matrix Sparsität | 30.0 | PARAM | ✅ Read (via forward) | ✅ Read (cutoff_sq) | linear | ❌ **Doppel-Def** |
| 5 | `eeq_distance_cutoff` (eeq_solver) | `eeq_solver.h:965` | EEQ A-Matrix Sparsität | **0.0** (matcht Fortran) | PARAM | ✅ Effektiv überschrieben durch #4 | ✅ ankommen via #4 | linear | ❌ **Doppel-Def** |
| 6 | `dispersion_cutoff` | `gfnff_method.cpp:9137,9153` | D3/D4 Pair-Generation Filter | 95.0 | `m_parameters.value()` | ✅ Used | ⚠️ N/A | linear | ✅ funktional |
| 7 | `disp.r_cut` (D3) | `gfnff_method.cpp:8994` | D3 Pair-Energie | √1500 ≈ 38.73 (Fortran-Match) | Hardcode | ✅ Used | ✅ Read aus `r_cut[tid]` SoA | linear | ✅ konsistent |
| 8 | `disp_cutoff_bohr` (D4) | `d4param_generator.cpp:1489` | D4 Pair-Generation | 60.0 | constexpr Hardcode | ✅ Used | ⚠️ N/A | linear | ⚠️ **D3≠D4** |
| 9 | `disp_cutoff_sq` (D4) | `d4param_generator.cpp:1490` | D4 Pair-Generation Filter | 3600.0 (= 60²) | constexpr Hardcode | ✅ Used | ⚠️ N/A | **squared** | ⚠️ **D3≠D4** |
| 10 | `coul.r_cut` (Struct-Default) | `gfnff_parameters.h:143` | Coulomb-Pair-Energie | **50.0** | Struct-Member-Default | ⚠️ überschrieben | ⚠️ Falls nicht überschrieben | linear | ❌ **Mismatch** |
| 11 | `coul.r_cut` (CPU-Assigned) | `gfnff_method.cpp:8467` | Coulomb-Pair-Generation | **100.0** (überschreibt #10) | Hardcode | ✅ Used | Wird mit Pair hochgeladen | linear | ❌ Mismatch zu #10 |
| 12 | `hb.r_cut` | `gfnff_parameters.h:160` / set in gfnff_method.cpp | HB-Pair-Energie | 50.0 | Struct + zugewiesen | ✅ Used | ✅ Read | linear | ✅ konsistent |
| 13 | `xb.r_cut` | `gfnff_parameters.h:186` / set in gfnff_method.cpp:7587 | XB-Pair-Energie | 20.0 (zugewiesen) | Struct + zugewiesen | ✅ Used | ✅ Read | linear | ✅ konsistent |
| 14 | `rep.r_cut` (bonded) | `gfnff_parameters.h:125` | Repulsion (bonded) | 20.0 (zugewiesen) | Struct + zugewiesen | ✅ Used | ✅ Read | linear | ✅ konsistent |
| 15 | `rep.r_cut` (non-bonded) | `gfnff_parameters.h:125` | Repulsion (non-bonded) | 20.0 (zugewiesen) | Struct + zugewiesen | ✅ Used | ✅ Read | linear | ✅ konsistent |
| 16 | `eeq_pcg_nfrag_threshold` | `eeq_solver.h:956` | EEQ-Solver-Schwelle (nfrag) | 4 | PARAM | ✅ | ⚠️ N/A | int | ✅ funktional |
| 17 | `pcg_large_threshold` | `eeq_solver.h:954` | EEQ-Solver-Schwelle (N) | 500 | PARAM | ✅ | ⚠️ N/A | int | ✅ funktional |
| 18 | `convergence_threshold` (EEQ) | `eeq_solver.h:928` | EEQ PCG-Konvergenz | 1e-6 | PARAM | ✅ | ✅ | rel. | ✅ funktional |
| 19 | `hb_cell_list_min_atoms` | `gfnff.h:233` | HB/XB Cell-List Aktivierung | 800 | PARAM | ✅ | ⚠️ N/A | int (atom-count) | ✅ funktional |

**Legende:**
- ✅ funktional / konsistent — Parameter wird wie erwartet angewendet
- ⚠️ Diskrepanz, aber semantisch akzeptabel oder dokumentationswürdig
- ❌ Inkohärenz, sollte gefixt werden

## Top-Inkohärenzen (zu fixen)

### #1 — `eeq_distance_cutoff` Doppel-Definition (kritisch)

Zwei `PARAM(eeq_distance_cutoff, ...)`-Definitionen mit **verschiedenen Defaults**:
- `gfnff.h:229` Default `30.0`
- `eeq_solver.h:965` Default `0.0`

Forwarding-Logik in `gfnff_method.cpp:533` (`forwardEEQSolverParams`) überschreibt den eeq-solver-Wert mit dem gfnff-Wert. Effektiv gilt: **30.0**, was **nicht** die Fortran-Referenz ist (XTB GFN-FF nutzt 0 = kein Cutoff für Phase-2-EEQ).

**Folge:** Curcuma-Default-EEQ ist nicht bit-identisch zu XTB. Heute durch andere Drift-Quellen kaschiert (siehe WP-V), aber sobald jene aufgeklärt sind, bleibt dieser systematische Offset.

### #2 — CN-Derivate-Cutoff effektiv tot

`gfnff_method.cpp:896` liest `cn_cutoff_bohr` aus dem Parameter (Default 6.0 Bohr). Aber dann an drei Stellen (Z. 3436, 6762, 6889) der CN-**Derivate**-Berechnung wird **hardcoded `1600.0 Bohr² = 40 Bohr`** als Schwellwert genutzt. Der Parameter ist tot — das Read in Z. 896 ist Dead Code.

**Folge:** Setzen von `-cn_cutoff_bohr 0` oder `-cn_cutoff_bohr 10` ändert nicht das Verhalten der CN-Derivate.

**Semantik:** CN selbst nutzt 6.0 Bohr (Neighbor-List), CN-Derivate brauchen größeren Cutoff (40 Bohr) wegen Gradient-Reichweite. Das ist sachlich korrekt, aber das Naming ("`cn_cutoff_bohr`") suggeriert, dass derselbe Parameter beide deckt. Der Hardcoded-Wert sollte als `constexpr cn_deriv_cutoff_sq` mit Doc-Kommentar explizit gemacht werden.

### #3 — Coulomb r_cut CPU/GPU-Mismatch (kritisch)

- Struct-Default `coul.r_cut = 50.0` (`gfnff_parameters.h:143`)
- CPU-Assignment in `gfnff_method.cpp:8467`: `c.r_cut = 100.0`

Wenn die Pair-Liste **immer** durch CPU gesetzt wird und dann hochgeladen, ist der GPU-Wert auch 100.0 — konsistent. Wenn aber irgendein Pfad den Struct-Default nutzt (z. B. wenn die Pair-Liste auf GPU regeneriert wird), gibt es 50 vs 100 Mismatch.

**Folge:** Bei Geometrien mit weitesten Pairs zwischen 50 und 100 Bohr (selten, aber bei großen Systemen wie polymer.xyz möglich) verschiedene Energien CPU vs GPU.

### #4 — D3 vs D4 Dispersion-Cutoff (mittlere Priorität)

D3 nutzt `√1500 ≈ 38.73 Bohr` (Fortran-Match), D4 nutzt `60.0 Bohr` (in `d4param_generator.cpp:1489`). Verschiedene Pair-Mengen werden generiert.

**Semantik:** D4 hat reichweitigere Polarisierbarkeitskoppelungen als D3, daher ist der größere Cutoff potenziell sachlich gerechtfertigt — aber das ist **nicht dokumentiert**. Falls beabsichtigt, als Code-Kommentar erklären; falls Versehen, den D3-Wert übernehmen.

### #5 — CN wird 3× berechnet (Performance, kein Korrektheits-Issue)

In `gfnff_method.cpp` wird CN an drei Stellen berechnet: Z. 896 (initial), Z. 3436 (für Winkel-Loop), Z. 6762 / 6889 (für CN-Derivate). Kein gemeinsames Caching. Performance-Kosten: **3× O(N²)** statt 1× — bei N=6200 sind das ~300 ms unnötig.

**Out of Scope für WP-C** — eigenes WP wert.

## Cutoff-Hierarchie für G2c-Aktivierung (Empfehlung)

Wenn `eeq_distance_cutoff > 0` aktiviert wird, sollten alle Pair-bezogenen Cutoffs konsistent sein:

```
eeq_distance_cutoff (Master)
    ├─ EEQ A-Matrix-Aufbau (Site 1)            ← skip Coulomb-Off-Diag wenn r > cutoff
    ├─ Coulomb-Pair-Liste / coul.r_cut (Site 2)  ← min(default, eeq_distance_cutoff)
    │   └─ CPU-Coulomb-Energie (Site 3)        ← liest coul.r_cut automatisch
    │   └─ GPU-Coulomb-Kernel (Site 4)         ← liest r_cut[tid] automatisch
    └─ EEQ-Gradient Term 1b (Site 5)           ← skip wenn r > cutoff
```

**Default 0** (kein Cutoff aktiv) → alle Sites verhalten sich wie heute (volle O(N²)).
**Cutoff > 0** → opt-in, alle Sites müssen identisch filtern, sonst HF-Verletzung.

Die anderen Cutoffs bleiben **unabhängig** von `eeq_distance_cutoff`:
- `cn_cutoff_bohr` (CN-Berechnung) — physikalisch von Coulomb entkoppelt
- `hb.r_cut`, `xb.r_cut`, `rep.r_cut` — Term-spezifisch
- `disp.r_cut` (D3/D4) — Dispersions-spezifisch

## Stand zu G2c-Voraussetzungen (4 Pflicht-Sites)

| Site | Komponente | CPU heute | GPU heute | Bedarf |
|------|------------|-----------|-----------|--------|
| 1 | EEQ A-Matrix-Aufbau | volle N×N (Cutoff inaktiv trotz Param 30.0?) | `cutoff_sq` Parameter wird übergeben | ⚠️ CPU effektiv 30.0 wegen Doppel-Def, sollte aber 0 sein |
| 2 | Coulomb-Pair-Liste | `coul.r_cut = 100.0` hardcoded | hochgeladen | ⚠️ Cutoff > 0 hat keinen Effekt, weil 100.0 immer übersteuert |
| 3 | CPU-Coulomb-Energie | filtert `coul.r_cut` | — | ✅ folgt Site 2 |
| 4 | GPU-`k_coulomb` | — | filtert `r_cut[tid]` | ✅ folgt Site 2 |
| 5 | EEQ-Gradient Term 1b | kein expliziter geometrischer Cutoff | kein expliziter Cutoff | ❌ muss bei Cutoff > 0 hinzugefügt werden |

**Bilanz:** G2c ist heute **strukturell halb-implementiert**. Cutoff-Mechanik existiert in Sites 2/3/4, aber:
- Site 1 hat zwei verschiedene Defaults (Inkohärenz #1)
- Site 2 wird durch hardcoded 100.0 übersteuert (Inkohärenz #3)
- Site 5 fehlt vollständig

## Beziehung zu anderen WPs

- **WP-V (Gradient-Validierung):** orthogonal. WP-V testet die Gradient-Implementierungen gegen XTB; WP-C testet Cutoff-Konsistenz. Default `eeq_distance_cutoff = 0` ist Fortran-Match, also pre-WP-C-Diff zur Fortran-Referenz bleibt durch WP-C unverändert.
- **WP6 (Coulomb-Cell-List):** vorbedingung-blockiert auf WP-C. Sobald G2c sauber, Cell-List nutzt `eeq_distance_cutoff` als Bin-Größe.
- **CN-Performance-Cleanup** (Audit-Punkt #5): unabhängiges WP, eliminiert die 3× CN-Berechnung.

## Änderungs-Historie

| Datum | Änderung | Quelle |
|-------|----------|--------|
| Mai 2026 | Initial Audit | Phase-1-Erkundung WP-C |
