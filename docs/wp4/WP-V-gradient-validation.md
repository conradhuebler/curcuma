# WP-V — Gradient-Validierung CPU vs GPU vs Fortran-Referenz

**Status:** 🆕 Vorgeschlagen
**Aufwand:** 1–3 Tage Diagnose + Term-Lokalisation, dann je nach Befund 1–4h pro Bug-Fix
**Erwarteter Nutzen:** Reproducibility — Voraussetzung für die Behauptung "100 % Referenz-Genauigkeit" laut CLAUDE.md, MD-Long-Run-Stabilität, sauberer Optimizer-Pfad. Kein direkter Performance-Effekt.
**Voraussetzung:** keine

## Symptom

CPU- und GPU-Gradient unterscheiden sich untereinander **und** beide weichen von der Fortran-Referenz (XTB GFN-FF) ab. Beobachtbar in MD durch nicht-physikalischen Energieaustausch mit dem Thermostat. Energie ist auf < 1 μEh konsistent zwischen allen drei — der Drift sitzt im Gradient.

Konkrete vorhandene Datenpunkte (aus `docs/GPU_GFNNF_DISCREPANCIES.md`, Mar–Apr 2026):

| Größe | Native CPU | Native GPU | XTB Fortran |
|-------|-----------|------------|-------------|
| Energie (complex 231-atoms) | -37.24064797 | -37.24064797 | -37.24064869 (+0.72 μEh) |
| MD heat-bath exchange (polymer) | weiter weg | näher an XTB | Referenz (-7.91275 Eh) |
| Acetic-acid-dimer Repulsion-Gradient | — | — | bis 3.7e-3 Eh/Bohr Drift CPU↔GPU |

## Was nicht der Bug ist

Aufgrund der WP4-Reihe bekannte ausgeschlossene Quellen:

- **Nicht WP4** (CN-Derivate-Pair-List): Energie & Gradient bit-identisch nach WP4 auf kleinen Molekülen, MD-Drift war pre-WP4 schon da
- **Nicht WP-G** (RowMajor-Layout): Gradient bit-identisch zu pre-WP-G post-Fix von dangling-ref (Commit `b753f79`)
- **Nicht WP5** (D4 GW Refactor): mathematisch identisch
- **Nicht WP-G's GPU-Crash**: Commit `b753f79` hat den dangling-ref-Bug gefixt — GPU läuft jetzt stabil, aber das war ein orthogonales Problem

Die hier dokumentierte Diskrepanz war **pre-existing** vor allen WP4-Reihen-Änderungen. Erste Notiz in `GPU_GFNNF_DISCREPANCIES.md` ist Mar 27 2026.

## Existierende Diagnose-Infrastruktur

Bereits implementiert (Commit `ad10ae8`, April 2026): umfangreiche `verbosity 3` Debug-Outputs, die ohne weitere Code-Änderungen aktiviert werden können:

```bash
# Vollständige GPU/CPU Gradient-Diagnose
./curcuma -sp acetic_acid_dimer.xyz -method gfnff -gpu cuda -verbosity 3 \
  > diag_gpu.log 2>&1
./curcuma -sp acetic_acid_dimer.xyz -method gfnff -threads 1 -verbosity 3 \
  > diag_cpu.log 2>&1
```

Output-Sektionen pro Lauf:
1. **PRE-CN GRADIENT** — Gradient vor CN chain-rule (Per-Atom)
2. **POST-CN GRADIENT** — nach CN chain-rule
3. **PER-COMPONENT GRADIENT** — Zerlegung nach Term (bond, angle, repulsion, dispersion, coulomb, hb)
4. **dEdcn / qtmp / combined** — CN-Empfindlichkeit pro Atom
5. **GPU PER-KERNEL DECOMPOSITION** — pro Pipeline-Stage (rep, bonds, angles, sB_rest, 3body, coul_ph2, CN_chain)
6. **CN COMPARISON** (GPU vs CPU)
7. **HB/XB PAIR COMPARISON** (GPU pair-counts und erste 20 Paare)
8. **ENERGY TERMS** (GPU + CPU, alle 14 Komponenten)

Plus existierende Test-Cases:
- `test_cases/test_hb_gradient_isolation.cpp` — isolierter HB-Gradient-Test (FD)
- `test_cases/test_gfnff_gpu.cpp` — CPU vs GPU Vergleich
- `test_cases/test_gfnff_numgrad.cpp` — numerischer Gradient-Test (nur CPU heute)

## April-2026-Befund: GPU näher an Fortran als CPU

`GPU_GFNNF_DISCREPANCIES.md` Update vom 29.04.2026 enthält wichtige Revision:

> The GPU implementation produces MD trajectories that agree better with XTB GFN-FF than the CPU implementation, as measured by heat-bath exchange values. **The CPU gradient may contain errors** not present in the GPU path — possibly in HB gradient distribution, CN chain-rule accumulation order, or atomicAdd vs sequential summation.

Heißt: das traditionelle Annahmen-Modell ("GPU ist die unsichere Variante") ist möglicherweise umgekehrt. **CPU könnte die fehlerhafte Variante sein.** Die früheren CPU-vs-GPU-Vergleiche haben CPU als Referenz benutzt — falsche Annahme.

## Strategie

Der existierende `GPU_GFNNF_DISCREPANCIES.md`-April-Befund gibt die richtige Marschrichtung vor:

> **Priority**: Compare GPU and CPU gradients separately against XTB numerical gradients (not against each other).

WP-V macht genau das, in vier Phasen:

### Phase 1 — Numerischer XTB-Referenz-Gradient erzeugen (1 Tag)

XTB liefert nur den Total-Energy-Gradient direkt (`xtb mol.xyz --grad`). Für **Per-Term**-Vergleich brauchen wir:

1. **Total-FD-Gradient via XTB:** XTB als black-box mit ±h Geometrie-Pertubationen → numerischer Gradient. h ~5e-4 Bohr für gute FD-Genauigkeit.
2. **Per-Term FD via Energie-Decomposition:** XTB mit `-D xtb.config` und nur einem Term aktiviert → isolierter Per-Term-FD-Gradient.

Test-Set:
- `acetic_acid_dimer.xyz` (16 Atome, hat HB)
- `triose.xyz` (22 Atome, single fragment)
- `caffeine.xyz` (24 Atome, ring + HB)
- `complex.xyz` (231 Atome — der April-Datenpunkt)

Output: Tabelle mit `(atom, dim) -> (XTB FD, CPU anal, GPU anal)` pro Term. Diff CPU-XTB und GPU-XTB beide aufgetragen.

### Phase 2 — Term-Bisect (1 Tag)

Aus Phase 1: identifiziere den Term, der den größten CPU-XTB-Diff (oder GPU-XTB-Diff) hat. Kandidaten aus April-Update:

- **HB-Gradient-Distribution:** abhgfnff_eg.f90 vs C++ Implementation
- **CN-Chain-Rule:** Reduktion-Reihenfolge in `applyAdd` (post-WP4) — bit-identisch zu pre-WP4? Verify.
- **atomicAdd vs sequential:** GPU benutzt atomicAdd, CPU sequenziell. Reihenfolge der Summen-Akkumulation pro Atom unterschiedlich.
- **Repulsion (bonded + non-bonded):** CLAUDE.md sagt "⚠️ Partial" — bekannt instabil.

### Phase 3 — Single-Term-Test (variabel)

Pro identifiziertem Term: einen isolierten Test (alle anderen Terme deaktiviert via Parameter-Flags) → bit-präziser Diff zur Fortran-Referenz. Reproduziert die Diskrepanz auf einem Minimal-System (z. B. nur HB zwischen 2 Wasser-Molekülen).

### Phase 4 — Fix + Verify

Abhängig vom Term:
- Reduktion-Reihenfolge angleichen (deterministisch sortiert nach Atom-Index, nicht nach Pair-Index)
- atomicAdd-Pattern auf CPU emulieren (wenn das wirklich die Quelle ist — Per-Atom-Gradient-Buffer mit Reduktion am Ende)
- Bug in Term-Implementierung beheben (z. B. HB-Gradient-Distribution-Formel)

Nach jedem Fix: Phase-1-Test wiederholen, Diff sollte um Faktor 10+ schrumpfen.

## Critical Files (vermutete Quellen)

| Datei | Zeile / Funktion | Verdächtige Operation |
|-------|------------------|------------------------|
| `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` | `CalculateGFNFFHydrogenBondContribution` | HB-Gradient-Distribution Reduktions-Reihenfolge |
| `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` | `CalculateGFNFFNonbondedRepulsionContribution` | Bekannt "⚠️ Partial" laut CLAUDE.md |
| `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` | `CalculateGFNFFBondedRepulsionContribution` | Bekannt "⚠️ Partial" |
| `src/core/energy_calculators/ff_methods/forcefield.cpp:2688-2758` | CN chain-rule via `CNDerivStore::applyAdd` | Reduktion-Reihenfolge — Pair-Order vs SpMatrix-Order. Verify bit-equivalence pre-WP4 |
| `src/core/energy_calculators/ff_methods/cuda/gfnff_kernels.cu` | `k_cn_chainrule`, `k_hbonds`, `k_xbonds` | atomicAdd vs sequential — verschiedener Akkumulations-Stack |
| `src/core/energy_calculators/qm_methods/gfnff_gpu_method.cpp:506` | `addExternalGradient` | CPU/GPU-Boundary für externe Gradient-Beiträge |
| `external/gfnff/src/gfnff_engrad.F90` | bond-force, HB, repulsion-Sections | Fortran-Referenz für direkten Code-Diff |

## Verification

```bash
cd release && make -j4

# Phase-1-Run: XTB FD-Referenz
xtb test_cases/molecules/larger/acetic_acid_dimer.xyz --gfnff --grad
# Manuelles FD via Pertubation pro Komponente — Skript dafür schreiben:
#   for atom in 0..15; for dim in 0..2; for h in {-, +}: 
#     perturb geom[atom][dim] by ±h; run xtb --gfnff; energy[h]
#     dE/dx = (E_plus - E_minus) / (2h)

# Phase-1-Vergleich
./curcuma -sp acetic_acid_dimer.xyz -method gfnff -gpu cuda -verbosity 3 > diag_gpu.log
./curcuma -sp acetic_acid_dimer.xyz -method gfnff -threads 1 -verbosity 3 > diag_cpu.log
# Diff parsen mit existierender Per-Atom-Output-Struktur

# Phase-3-Single-Term-Test (z. B. nur HB)
./curcuma -sp acetic_acid_dimer.xyz -method gfnff \
  -gfnff.dispersion_enabled=false -gfnff.repulsion_enabled=false \
  -gfnff.coulomb_enabled=false -threads 1 -verbosity 3
# Vergleich pro Atom CPU vs FD vs GPU
```

## Akzeptanzkriterien

1. ⚙️ Phase-1-Tabelle für alle 4 Test-Moleküle vollständig (CPU-XTB-Diff + GPU-XTB-Diff pro Atom).
2. ⚙️ Hot-Term identifiziert (mind. 80 % des Total-Drifts erklärt durch einen einzigen Term).
3. ⚙️ Per-Term-Test reproduziert die Diskrepanz auf einem Minimal-System.
4. ⚙️ Fix umgesetzt, Diff CPU-XTB **und** GPU-XTB für den fixed Term um Faktor ≥10 reduziert.
5. ⚙️ MD heat-bath exchange auf polymer-system: CPU-XTB-Diff < 0.05 Eh (vs. heute > 0.2 Eh).
6. ⚙️ `test_gfnff_numgrad` und `test_gfnff_gpu` weiterhin grün.
7. Operator setzt ✅ TESTED.

## Risks

- **Hoch.** Numerischer FD via XTB ist Skripting-Aufwand. Isolation einzelner Terme via Parameter-Flags funktioniert nicht für alle Terme (manche sind hart-gekoppelt, z. B. CN-abhängige Bonds).
- **Term-spezifisch:** Falls die Quelle eine subtile Reduktions-Reihenfolge in HB ist, ist der Fix mathematisch trivial aber muss bit-präzise gegen Fortran getestet werden.
- **GPU vs CPU Asymmetrie:** Wenn GPU näher an Fortran ist als CPU (April-Befund), ist das **kein** Bug der GPU-Implementation — sondern legt einen CPU-Bug frei. Fix kann dort sein, **wo wir bisher nicht gesucht haben**.
- **Multi-Thread-Race-Confounder:** Triose T=4 zeigt 1.47-mEh-Drift zwischen Läufen (siehe WP-R). Vor WP-V muss WP-R aufgeklärt sein, sonst überlagern sich die beiden Effekte.

## Beziehung zu anderen WPs

- **Vorbedingung:** WP-R (Multi-Thread-Race) sollte zuerst gefixt sein, sonst sind alle CPU-Gradient-Daten verrauscht.
- **WP4 / WP-G / WP5:** alle ausgeschlossen als Quelle (siehe oben). Verify durch ein vorab-Test: pre-WP4-Commit auschecken, einmal MD heat-bath measuren — sollte denselben Drift zeigen.
- **Bestehende Doku:** `docs/GPU_GFNNF_DISCREPANCIES.md` (225 Zeilen, Mar–Apr 2026) und `docs/GFNFF_GRADIENTS.md` (45 Zeilen). WP-V baut darauf auf, ersetzt sie nicht.

## Out of Scope

- Performance-Optimierung der Test-Infrastruktur (z. B. paralleles XTB-FD).
- Solvation-Gradient (ALPB) — separates Issue, ALPB ist als "not validated" markiert in CLAUDE.md.
- Periodic-Boundary-Conditions — nicht implementiert.
- Mixed-Precision-Wiederaktivierung auf GPU — separates WP wenn nötig.

## Vorbedingung

- WP-R abgeschlossen (oder zumindest verifiziert dass T=4-Race nicht der Hot-Confounder hier ist)
- XTB GFN-FF binary verfügbar (`xtb` mit `--gfnff`-Flag)
- Test-Moleküle vorhanden (alle in `test_cases/molecules/larger/`)
