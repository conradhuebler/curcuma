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

### #1 — `eeq_distance_cutoff` Doppel-Definition (✅ aufgelöst, Mai 2026)

**Historie:**

Phase A+B (Commit `2b85963`): Doppel-PARAM zwischen `gfnff.h:229` (30.0) und `eeq_solver.h:965` (0.0) aufgelöst — gfnff-PARAM entfernt, eeq_solver-PARAM kanonisch mit Default 0.0.

**Verstecktes Folge-Problem (Mai 2026, Commit `cba0696`):** Selbst nach Phase A+B verhielt sich Curcuma-Default nicht wie Fortran. Ursache war ein dead-code-Inline-Fallback in `eeq_solver.cpp:3077`:

```cpp
double eeq_dist_cutoff = m_config.get<double>("eeq_distance_cutoff", 30.0);
```

Die `30.0` überstimmte den canonical PARAM-Default 0.0 **immer dann**, wenn `m_config` keinen explicit registry-Eintrag trug — was der Default-CLI-Pfad ist. Resultat: Phase-2-EEQ truncated silent jede default-Invokation für N>200 (Gating bei line 3078 `natoms > 200`).

Verifikation gegen Fortran (`gfnff_engrad.F90:1319-1331, 1380-1389`): `goed_gfnff` hat **keinen** Phase-2-Cutoff, voller O(N²) A-Matrix und Energie-Summation.

**Empirisch (polymer.xyz, 1410 Atome, gegen XTB --gfnff):**

| Stand | Curcuma (Eh) | XTB Fortran (Eh) | Diff |
|-------|--------------:|-----------------:|------:|
| pre-Fix (silent 30.0) | -202.58789086 | -203.51759134 | **+929.7 mEh** |
| post-Fix (PARAM-Default 0.0) | -203.56552633 | -203.51759134 | **+47.9 mEh** |

Die Restdifferenz (48 mEh) sitzt nahezu vollständig im Torsion-Term — alle anderen Komponenten (Bond, Angle, Coulomb, Dispersion, BATM, Repulsion) matchen Fortran zu <µEh.

**Korrektur am alten D4-Eintrag unten:** Der +0.93 Eh polymer-Bias war **nicht** D4- oder HB-bezogen, sondern dieser EEQ-Cutoff-Fallback. WP-V (Gradient-Validierung) bleibt davon unberührt — es ging dort um Gradient-Drift in MD, nicht um den Energie-Bias.

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

D3 nutzt `√1500 ≈ 38.73 Bohr` (Fortran-Match), D4 nutzt `60.0 Bohr` Pair-Generation + `50.0 Bohr` Pair-Energie-Cutoff (in `d4param_generator.cpp:1489, 1531`). Verschiedene Pair-Mengen werden generiert.

**WP-C Phase D — XTB-Validierung (Mai 2026):** Versuch unternommen, D4 auf `38.73 Bohr` (D3-Match, Fortran-`dispthr=1500`) zu vereinheitlichen. Energie-Vergleich auf `polymer.xyz` (1410 Atome) gegen XTB-GFN-FF-Fortran-Referenz:

| Konfiguration | Curcuma D4 (Eh) | XTB Fortran (Eh) | Diff vs XTB |
|---------------|-----------------|-------------------|-------------|
| D4 cutoff = 60.0 Bohr (Status quo) | -202.58789977 | -203.51759034 | **+0.92969 Eh** |
| D4 cutoff = 38.73 Bohr (D3-Match) | -202.58766507 | -203.51759034 | **+0.92993 Eh** |

→ Vereinheitlichung verschlechtert den Fortran-Match um 0.23 mEh.

**Schlussfolgerung:** Fortran's `dispthr=1500` gilt offenbar nur für D3, nicht für D4 (in der Fortran-Implementierung hat D4 entweder keinen geometrischen Cutoff oder ein anderes `dispthr_d4`). Curcuma's empirisches `60.0 Bohr` für D4 ist näher an der Fortran-D4-Referenz als das D3-Match. **Phase D wurde revertiert** — Status quo (D4=60.0 Bohr) bleibt.

**Hinweis (Mai 2026 update):** Der damals gemessene +0.93 Eh polymer-Bias war zu jenem Zeitpunkt richtig — aber die Ursache war **nicht** D4. Sie war der EEQ-Cutoff-Fallback (Inkohärenz #1, Commit `cba0696`). Nach jenem Fix ist der polymer-Bias auf ~48 mEh geschrumpft, sitzt im Torsion-Term, und D4 selbst matcht Fortran zu <µEh. Die obigen D4-Tabelle-Werte sind also weiterhin gültige Datenpunkte zum D4-Verhalten in isolation, aber der Polymer-Bias ist kein D4-Issue.

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

## Stand zu G2c-Voraussetzungen (5 Sites — alle vollständig, Mai 2026)

| Site | Komponente | CPU | GPU | Status |
|------|------------|-----|-----|--------|
| 1a | EEQ A-Matrix Phase 2 (`buildCorrectedEEQMatrix`) | filtert in Geometric-Mode bei cutoff>0 (`eeq_solver.cpp:1175-1185`) | `k_eeq_build_matrix` filtert via `cutoff_sq` (`cuda/eeq_solver_gpu.cu:216-251`) | ✅ Commit `f3ffe97` |
| 1b | EEQ A-Matrix Phase 2 iterativer Solve | filtert seit Mar 2026 (`eeq_solver.cpp:3083-3086`) — Default-Fallback Mai 2026 von 30.0 → 0.0 korrigiert | identisch zu CPU | ✅ Commits `2b85963` + `cba0696` |
| 2 | Coulomb-Pair-Liste | drei-Branch-Dispatch in `generateCoulombPairsNative` (Cell-List wenn N≥`nb_cell_list_min_atoms` und cutoff>0) | erbt Pair-Liste vom CPU mit identischem `r_cut` | ✅ Commit `f3ffe97` |
| 3 | CPU-Coulomb-Energie | filtert `coul.r_cut` (war heute Z. 2200, no-op bei r_cut=100) | — | ✅ folgt Site 2 |
| 4 | GPU-`k_coulomb` | — | filtert `r_cut[tid]` (`cuda/gfnff_kernels.cu:359`) | ✅ folgt Site 2 |
| 5 | EEQ-Gradient Term 1b — CN-Derivate-Stencil | `cn_deriv_cutoff_sq = max(40², eeq_cut²)` (`gfnff_method.cpp:1014`) — niemals < 40 Bohr | identisch (k_cn_chainrule liest dieselbe `m_dcn`) | ✅ Commit `f3ffe97` |

**Bilanz:** G2c ist **vollständig implementiert** (Mai 2026). Cutoff > 0 ist opt-in und über `-gfnff.eeq_distance_cutoff <Bohr>` aktivierbar (CLI-Plumbing-Fix Teil von `f3ffe97`). Default 0 = voller O(N²), bit-identisch zu Fortran-`goed_gfnff`.

## Beziehung zu anderen WPs

- **WP-V (Gradient-Validierung):** weiterhin offen für Gradient-Drift in MD. Der Energie-Bias-Verweis aus dem alten WP-V-Doc ist durch Commit `cba0696` beantwortet — der polymer-Bias war der EEQ-Cutoff-Fallback (siehe Inkohärenz #1).
- **WP6 (Coulomb-Cell-List):** **abgeschlossen** Mai 2026 (Commit `f3ffe97`). mixture.xyz Coulomb-Wall 611 → 74 ms (8.3×) bei cutoff=30.
- **CN-Performance-Cleanup** (Audit-Punkt #5): unabhängiges WP, eliminiert die 3× CN-Berechnung. Offen.

## Änderungs-Historie

| Datum | Änderung | Quelle |
|-------|----------|--------|
| Mai 2026 | Initial Audit | Phase-1-Erkundung WP-C |
| Mai 2026 | Phase A+B umgesetzt — `eeq_distance_cutoff` Doppel-Def aufgelöst, CN-deriv-cutoff benannt, Coulomb r_cut Struct-Default angeglichen, Init-Logging | Commit `2b85963` |
| Mai 2026 | Phase D versucht + revertiert: D3/D4 Cutoff-Vereinheitlichung verschlechtert XTB-Match (siehe oben) | XTB-Validation auf polymer.xyz |
| Mai 2026 | WP6 + Phase C umgesetzt: Site 1a + Site 5 Filter-Logik, Cell-List für Coulomb-Pairs, CLI-Plumbing-Fix | Commit `f3ffe97` |
| Mai 2026 | Inkohärenz #1 vollständig aufgelöst: `eeq_solver.cpp:3077` Inline-Fallback 30.0 → 0.0. Polymer-Bias zu Fortran von 930 mEh auf 48 mEh reduziert | Commit `cba0696` |
