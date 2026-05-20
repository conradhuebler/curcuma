# WP-FF-SoA-Extension — Pass-A/B/C-Pattern auf alle GFN-FF-Inner-Loops

**Status:** 🆕 Vorgeschlagen (Mai 2026, Erweiterung von WP-FF-SoA-pair-loops.md)
**Aufwand:** ~2-3 Tage
**Erwarteter Nutzen:** Polymer N=1410, FAST_EXP=ON: **0.3-1.0 s / 100 MD-Steps** je nach Term.
**Voraussetzung:** WP-D Stage B etablierte `fast_exp_neg_sq_block`. WP-FF-DistMatrix-Sharing optional.

## Hypothese

WP-D Stage B (Mai 2026) hat das Pass-A/B/C-SoA-Pattern für den `dcn`-Loop etabliert. WP-FF-SoA-pair-loops.md hat die Erweiterung auf Repulsion (bonded + nonbonded) implementiert — Ergebnis: Repulsion-Listen zu klein, SoA-Overhead > exp-Gewinn (siehe Befund im WP).

**Aber**: einige Inner-Loops in `forcefieldthread.cpp` haben tausende von Iterationen pro Step und nicht-trivialen `exp`/`erf`-Anteil. Diese sind echte Kandidaten:

| Term | Anzahl Paare/Tripel (polymer) | exp/erf pro Iter | exp-Gesamt /Step | SoA-Hebel |
|------|-------------------------------|------------------|------------------|-----------|
| **CoulombContribution** | ~10⁶ (alle nicht-1,2/1,3-Paare) | 1 erf | ~10⁶ erf | **0.5-1.5 ms** |
| **HydrogenBondContribution** | ~10² A-H-B Tripel | 1-2 exp | ~200-500 | klein |
| **HalogenBondContribution** | ~10⁰ Tripel | 1-2 exp | ~20 | irrelevant |
| **BatmContribution** | ~10³ Tripel | 3+ exp | ~3-10k | 0.2-0.5 ms |
| **AngleContribution (gfnffdampa)** | ~5000 Angles × 2 | 1 exp | 10k | 0.1-0.3 ms |

Der größte Hebel ist **Coulomb**: pro Step ~10⁶ `erf(γ_ij · r_ij)/r_ij`-Aufrufe. Selbst SIMD-erf-Batched (oder Padé-Approximation) würde hier dramatisch helfen.

## Aufgabe

### 1. Erweiterung von `fast_exp.h` um erf

```cpp
// src/core/energy_calculators/ff_methods/fast_exp.h
namespace curcuma::gfnff {
    // existing: void fast_exp_neg_sq_block(const double* x, double* out, size_t n) noexcept;
    void fast_erf_block(const double* x, double* out, std::size_t n) noexcept;
    // Padé-Approximation der erf: AVX2 8x-Throughput vs scalar
}
```

Implementierung via Padé- oder Abramowitz/Stegun-Approximation (Genauigkeit ~1e-7 reicht für EEQ).

### 2. CoulombContribution SoA-Pattern

`forcefieldthread.cpp:CalculateGFNFFCoulombContribution` (~L2299) restrukturieren:

```cpp
// Pass A: collect surviving pairs (within cutoff)
int K = 0;
buf_x.resize(N_max);                 // γ_ij · r_ij
buf_extra.resize(N_max, 4);          // q_i, q_j, gamma_ij, r_ij
for (const auto& coul : m_gfnff_coulombs) {
    if (filter(coul)) continue;
    buf_x[K] = coul.gamma * r_ij;
    buf_extra(K, 0) = coul.q_i;
    buf_extra(K, 1) = coul.q_j;
    buf_extra(K, 2) = coul.gamma;
    buf_extra(K, 3) = r_ij;
    ++K;
}

// Pass B: SIMD erf
curcuma::gfnff::fast_erf_block(buf_x.data(), buf_erf.data(), K);

// Pass C: scatter
for (int k = 0; k < K; ++k) {
    double erfval = buf_erf[k];
    double r      = buf_extra(k, 3);
    double q_q    = buf_extra(k, 0) * buf_extra(k, 1);
    double A      = erfval / r;
    local_E += q_q * A;
    // gradient scatter
}
```

### 3. BatmContribution

Triplet-Loop mit Tang-Toennies (3+ exp pro Triplet). Analog Pass-A/B/C mit `fast_exp_neg_sq_block`. Wins ~0.3 ms.

### 4. AngleContribution `gfnffdampa`

Pro Angle 2 `exp`-Aufrufe. SoA-Pattern kann hier in der äußeren Angle-Loop kollektieren.

### 5. HBond Re-Refactor

WP-FF-SoA-pair-loops.md Status: Repulsion-SoA neutral. HBond hat verschachtelten Inner-Loop (A-H-B-Tripel + per-Tripel Damping). Wert prüfen, vielleicht zu klein.

### 6. Conditional Compilation

`USE_GFNFF_FAST_EXP=ON` aktiviert die SoA-Pfade. Default `OFF` für Bit-Identität.

### Code-Anker

| Datei | Funktion | Zeilen |
|-------|----------|--------|
| `src/core/energy_calculators/ff_methods/fast_exp.h` | `fast_erf_block`, `fast_exp_neg_block` (neu) | — |
| `src/core/energy_calculators/ff_methods/fast_exp.cpp` | Implementierung | — |
| `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` | `CalculateGFNFFCoulombContribution` | L2299+ |
| ditto | `CalculateGFNFFBatmContribution` | L3795+ |
| ditto | `CalculateGFNFFAngleContribution` (gfnffdampa-Aufrufe) | L1065+ |
| ditto | `CalculateGFNFFHydrogenBondContribution` | L2520+ |

## Akzeptanzkriterien

Pro Term separat, **nicht en bloc** (sonst Lokalisierung von Drift unmöglich):

1. ⚙️ Default `USE_GFNFF_FAST_EXP=OFF`: bit-identisch zu pre-WP.
2. ⚙️ FAST_EXP=ON Single-Point caffeine/triose/polymer: Energie-Drift ≤ 10 µEh pro Term (erf-Approximation hat ~1e-7 relative Genauigkeit).
3. ⚙️ `test_gfnff_numgrad` grün.
4. ⚙️ Per-Term-Performance (Polymer 100-Step, single-thread, FAST_EXP=ON):
   - Coulomb: ≥ 0.5 s Ersparnis vs OFF
   - BATM: ≥ 0.2 s
   - Angle: messbar (Noise nicht überschritten)
5. ⚙️ MD-Long-Run 1000 Steps: Energie-Drift gegen FAST_EXP=OFF < 0.1 mEh.

## Risiken

1. **erf-Genauigkeit für Coulomb** — γ·r liegt typisch 0-50. Padé-Approximation in diesem Bereich genau machen, ggf. range-reduction.
2. **Kumulative Float-Drift** — alle 5 Terme mit Approximationen kombiniert: pro Single-Point könnte Drift ~50-100 µEh erreichen. Pro MD-Step ähnlich kleiner Drift. Long-Run noch beobachten.
3. **Educational-First-Constraint** — Helper-Code wird komplexer. Klare Doku in `fast_exp.h` ist Pflicht (was wird approximiert, welche Domäne, welche Genauigkeit).
4. **Maintenance** — wenn der Coulomb-Term sich ändert (z.B. neue Korrekturen), muss SoA-Pfad mit-gepflegt werden.

## Verifikation

```bash
cd release
cmake -DUSE_GFNFF_FAST_EXP=ON -DUSE_MARCH_NATIVE=ON ..
make -j4

# Per-Term-Korrektheit
ctest -R test_gfnff_numgrad --output-on-failure

# Per-Term-Diff
./release/curcuma -sp test_cases/molecules/larger/caffeine.xyz -method gfnff
# Notiere Energien, vergleiche mit FAST_EXP=OFF Build

# Performance
/bin/rm -f test_cases/molecules/larger/polymer.topo.json
time ./release/curcuma -md test_cases/molecules/larger/polymer.xyz -method gfnff \
   -print 1 -maxtime 1e2 -threads 1
# Erwartet: ≤ 9.0 s mit FAST_EXP=ON
```

## Beziehung zu anderen WPs

- **WP-FF-SoA-pair-loops** (Mai 2026, abgeschlossen): Vorläufer. Repulsion-SoA war neutral. Dieser WP konzentriert sich auf größere Listen (Coulomb) und exp-intensive Tripel (BATM).
- **WP-FF-DistMatrix-Sharing**: Pass-A profitiert davon, wenn `srab` zentral ist — Pass-A wird zur reinen Index-Loop ohne `sqrt`.
- **WP-EEQ-Matrix-Cache**: Coulomb-Matrix-Aufbau in EEQ verwendet auch `erf`. Falls dort SIMD-erf greift, sparen wir auch dort 1-2 ms/Step.

## Frühere Iteration (Dokumentation)

WP-FF-SoA-pair-loops.md beschreibt den ersten Versuch (Repulsion). Befund dort: SoA neutral wegen kleiner Listen. Dieses WP erweitert auf **größere** Listen (Coulomb 10⁶ Paare) und akzeptiert das frühere Lernergebnis.
