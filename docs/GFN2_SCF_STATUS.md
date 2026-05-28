# GFN2 SCF-Konvergenz Status

Separates Thema vom D4/Energie-Komponenten-Alignment — die SCF-Iteration selbst.
Master: [`GFN2_NATIVE_ROADMAP.md`](GFN2_NATIVE_ROADMAP.md).

## Bekannte SCF-Probleme

### complex (231 Atome, C₈₃H₁₂₂N₁₂O₁₄) — Divergenz

```
iter 0:  -328.10 Eh   (sehr nah an tblite −329.53!)
iter 1:  -326.60 Eh
iter 2:  -307.31 Eh   (+20 Eh)
iter 3:   -94.35 Eh   (+213 Eh)
iter 4:   -99.67 Eh
iter 5: +2419.24 Eh   (Vorzeichen-Flip nach DIIS-Aktivierung)
…
iter 149: +2764.9 Eh, max|dq|~2–6 (oszilliert)
```

**Diagnose:** Klassisches Charge-Sloshing. Die SCF startet nahe der richtigen
Lösung (bare-H0 Gauß ist eine gute Initial-Guess für gut-elektronische Systeme),
aber sobald DIIS bei iter 5 anspringt (`diis_start = 5` in
`xtb_native.cpp:150`), extrapoliert es die Dichte aus dem Konvergenz-Bassin
heraus. Die Damping (`m_scf_damping = 0.4` für Iter 1..4) reicht für dieses
große, polare System nicht.

**Vergangenes Workaround (HCN, Mai 2026):** `diis_start = 5` plus density damping
während Iter 1..4 hat die Charge-Transfer-State-Falle für kleine Moleküle (HCN,
Nitrile) gelöst. Das gleiche Mechanismus reicht hier nicht für 231 Atome.

### kleine Moleküle — Konvergenz OK

- H₂: 1 Iter
- H₂O, NH₃, CH₄: 5-10 Iter
- triose (66 Atome): 16 Iter, 365 ms

## Curcuma-SCF-Aufbau

`src/core/energy_calculators/qm_methods/xtb_native.cpp:127-200`:

```cpp
// 6. SCF-Loop mit DIIS-Beschleunigung
const double damp = m_scf_damping;     // density mixing (default 0.4)
const int diis_start = 5;              // warmup vor DIIS
DIISAccelerator diis(6);               // halte letzte 6 Fock-Matrizen
Matrix P_old;

for (iter = 0; iter < max_iter; ++iter) {
    m_pot.reset();
    addCoulombShellPotential(m_pot);
    addThirdOrderPotential(m_pot);
    if (GFN2) {
        addMultipolePotential(m_pot);
        addDispersionPotential(m_pot);
    }
    Matrix F = buildFock(m_H0, m_S, m_pot);

    if (iter >= diis_start) {
        diis.push(F, m_wfn.P, m_S);
        if (diis.size() >= 2) F = diis.extrapolate();
    }

    solveEigen(F, m_S);

    if (iter > 0 && iter < diis_start) {
        m_wfn.P = damp * m_wfn.P + (1.0 - damp) * P_old;
    }
    P_old = m_wfn.P;

    updatePopulations(m_S);
    // … Energie + Konvergenzprüfung
}
```

## Mögliche Fixes (priorisiert)

### Option 1: Längere/stärkere Damping-Warmup
Erhöhe `diis_start` von 5 auf z.B. 15 für große Systeme; ggf. lasse `damp = 0.3`.
Risiko: andere Tests (HCN-Familie, kleine Moleküle) brauchen kürzeres Warmup,
also sollte das von Größe/Polarität abhängen, nicht statisch.

### Option 2: Level-Shifting im H0
Standard-Trick: addiere einen positiven Shift zu den virtuellen Orbital-Energien
im ersten Iter-Block. tblite tut das? Code-Audit nötig.

### Option 3: EEQ-Guess statt bare-H0
Statt mit `q_sh = 0` zu starten, EEQ-Ladungen als Anfangsguess benutzen → Fock
mit korrekter Coulomb-Verschiebung sofort starten. Hat das schon mal eine
Charge-Transfer-Falle vermieden?

### Option 4: ADIIS / EDIIS statt Pulay-DIIS
Robuster gegen Charge-Sloshing, etwas langsamer pro Iter. Tblite nutzt Pulay
(soweit Code-Inspektion zeigt).

### Option 5: Tblite-SCF-Verhalten exakt nachbauen
Wenn tblite für complex konvergiert, schaue genau, wie es die ersten Iter
behandelt (Mixing, DIIS-Start, Konvergenz-Kriterium). Idealerweise Iter-für-Iter
gleich.

## Diagnostische Schritte

1. **Tblite auf complex laufen lassen** (via `dump_tblite_multipole` — das
   liefert tblite's konvergierte Dichte). Wenn tblite konvergiert, hat curcuma
   einen lösbaren Bug.
2. **Iter-für-Iter Diff** curcuma vs tblite für die ersten ~10 Iter (Energie,
   max|dq|). Lokalisiert, wo curcuma anfängt zu divergieren.
3. **Pre-DIIS Verhalten:** Iter 0..4 — passt curcuma's Fock zu tblite's?
   (Damping-Faktor, P-Update-Reihenfolge).

## Aufwand-Schätzung

Phase 1 (tblite auf complex testen, falls konvergiert): ~30 Min.
Phase 2 (Iter-für-Iter Diff): erfordert ein Diff-Tool ähnlich zu diag_d4.
   ~2-3 h.
Phase 3 (Fix implementieren): je nach Lokalisierung 1-4 h.

## Verwandte Dokumente

- [`GFN2_NATIVE_ROADMAP.md`](GFN2_NATIVE_ROADMAP.md) — Master
- [`GFN2_COMPONENT_AUDIT_PLAN.md`](GFN2_COMPONENT_AUDIT_PLAN.md) — Energie-Audit
