# WP4 — CPU-GFN-FF Performance-Arbeitspakete

Aufteilung der WP4-Maßnahmen aus [GFNFF-PERFORMANCE-ROADMAP.md](../GFNFF-PERFORMANCE-ROADMAP.md) (Mai 2026, Commit `fc9954f`) in 6 unabhängig umsetzbare Arbeitspakete.

## Ausgangsprofil (Mai 2026)

`mixture.xyz`, N=6200, nfrag=1400, single-point + gradient, 4 Threads:

| Phase | Zeit | Anteil |
|-------|------|--------|
| Total energy call | 3341 ms | 100% |
| CN + EEQ (als "serial CPU" gelabelt) | 1766 ms | 53% |
| Force-Field-Energy (multi-threaded) | 1574 ms | 47% |

GPU-Vergleich: 1.5 s. Ziel: CPU auf 0.6–0.8 s (Faktor vs. GPU 2.2× → ~1.3×).

## Reihenfolge

```
WP1 (Threading-Audit) ──► WP2 (EEQ batched parallel)
                          WP3 (Bond hotspot)
                          WP4 (CN SIMD)
                          WP5 (D4 Gaussian SIMD)
                          WP6 (Coulomb cutoff — Plan, blockiert auf G2c)
```

WP1 zuerst, da es den Befund "serial CPU" entweder als Label-Bug entlarvt (dann ist WP2 trivial) oder eine echte serielle Stelle aufdeckt (dann ist WP2 die Umsetzung).
WP3–WP5 sind unabhängig und können in beliebiger Reihenfolge oder parallel angegangen werden.
WP6 nur dokumentieren, Implementierung erst nach G2c-Validierung gegen XTB.

## Verifikation pro WP

```bash
cd release && make -j4
ctest -R "gfnff" --output-on-failure
./test_cases/test_gfnff_numgrad
cd ../test_cases/cli/simplemd/10_gfnff_polymer_md
bash compare_baseline.sh 4 gfnff "after_WP<N>_<short_id>"
```

Zusätzlich ist für WP2/WP3/WP6 das Multi-Fragment-Profil auf `mixture.xyz` (N=6200) Pflicht — der WP4-Plan stützt sich genau auf diesen Datenpunkt.

## Akzeptanz pro WP

Jedes WP ist erfolgreich abgeschlossen, wenn:
1. ⚙️ alle CPU-`gfnff`-CTests grün
2. ⚙️ `test_gfnff_numgrad` < 1e-5 Hartree/Bohr
3. ⚙️ Single-Point auf `acetic_acid_dimer.xyz` und `triose.xyz` reproduziert pre-WP-Energie auf < 1 µEh
4. ⚙️ Performance-Messung in `Status`-Tabelle der Roadmap eingetragen
5. Operator setzt `✅ TESTED` (nicht der AI)

## Cross-cutting-Regeln (verbindlich für alle WPs)

### Keine Approximation ohne Referenz-Fallback

Falls ein WP eine Approximation einführt (z. B. Polynom-Approximation für `erf`, reduzierte Reihenentwicklung für `exp`, Cutoff in einer bisher exakten Loop), gilt:

1. **Die exakte Referenz-Implementierung bleibt im Code.** Nicht löschen, nicht hinter `#ifdef DEBUG`.
2. **Runtime-Schalter** (ConfigManager-Parameter, snake_case, Default = exakt): z. B. `cn_use_approx_erf=false`, `coulomb_distance_cutoff_bohr=0.0` (0 ⇒ aus).
3. **Default-Verhalten ist immer Referenzgenauigkeit.** Approximation ist Opt-in.
4. **Akzeptanz-Test nutzt beide Pfade:** CTest mit Default, und CTest mit aktivierter Approximation. Beide grün.
5. **Genauigkeitsverlust dokumentieren:** absolute und relative Abweichung gegen exakten Pfad auf den fünf Standard-Test-Systemen (`acetic_acid_dimer`, `triose`, `caffeine`, `polymer`, `mixture`).

Damit ist die Approximation eine reine Performance-Option für große Systeme, ohne dass Standardläufe die Referenzgenauigkeit verlieren. Vorbild: `cn_cutoff_bohr` aus P2b — Default 6.0 ist Approximation für Performance, `cn_cutoff_bohr=0` mit `cn_accuracy=0` schaltet auf vollen O(N²)-Referenzpfad zurück.

### Keine Single-Mechanismus-Verdrängung

Wenn ein WP einen Algorithmus austauscht (z. B. Datenstruktur, Threading-Stil), den exakten alten Pfad **wenigstens für eine Release-Iteration** als Compile-Flag oder Setter-Schalter behalten — damit Operator A/B-Messung gegen den vorherigen Code machen kann.
