# AP 5b — Multipol-Integral-Pulay-Gradient (`cgto_multipole_grad`)

**Status:** Offen — unmittelbar nach AP5 Schritt 1 (direkter Interaktionsgradient)
**Erstellt:** 2026-04-26
**Vorgänger:** [`NATIVE_XTB_AP5_PLAN.md`](NATIVE_XTB_AP5_PLAN.md) (Validierung), [`NATIVE_XTB_AP4_PLAN.md`](NATIVE_XTB_AP4_PLAN.md)
**Vorbedingung:** AP5 Schritt 1 abgeschlossen (direkter Interaktionsgradient in `xtb_gradient.cpp` Sektion 5)

---

## Ziel

Der GFN2-Gradient ist nach AP5 Schritt 1 **formal inkonsistent**: Die Fock-Matrix enthält die Multipol-Potenziale `v_dp` und `v_qp` als explizite Beiträge (Sektion 2b in `calculateGradient()`), aber die Geometrieabhängigkeit dieser Integrale (`∂dp_int/∂R`, `∂qp_int/∂R`) fehlt im Gradienten.

AP5b ergänzt diesen fehlenden **Pulay-Term**:

```
dG_α += -P_μν · Σ_k [ ∂D_k(μν)/∂R_A_α · v_dp(k,A) + ∂D_k(μν)/∂R_B_α · v_dp(k,B) ]
        -P_μν · Σ_k [ ∂Q_k(μν)/∂R_A_α · v_qp(k,A) + ∂Q_k(μν)/∂R_B_α · v_qp(k,B) ]
```

TBLite-Referenz: `tblite/xtb/h0.f90` Zeilen 444–447, Integralfunktion: `tblite/integral/multipole.f90:multipole_grad_cgto` (Zeilen 522–657).

---

## Geänderte Dateien

| Datei | Änderung |
|-------|----------|
| `xtb_multipole_ints.hpp` | Neue Funktion `cgto_multipole_grad()` |
| `xtb_gradient.cpp` | Sektion 2b: Multipol-Pulay-Term mit `cgto_multipole_grad` |

---

## Aufgabenliste

### 5b.1 — `cgto_multipole_grad()` in `xtb_multipole_ints.hpp`

Port von `tblite/integral/multipole.f90:multipole_grad_cgto` (Zeilen 522–657).

**Signatur:**
```cpp
// Berechnet Dipol- und Quadrupol-Integrale UND deren Gradienten
// w.r.t. Zentrum A (ddA[3][3], dqA[3][6]) und Zentrum B (ddB[3][3], dqB[3][6]).
// ddA[α][k] = ∂D_k(μν)/∂R_A_α  (3 Raumrichtungen × 3 Dipol-Komponenten)
// dqA[α][k] = ∂Q_k(μν)/∂R_A_α  (3 Raumrichtungen × 6 Quadrupol-Komponenten)
void cgto_multipole_grad(
    const CGTO::Shell& shA, const CGTO::Shell& shB,
    double Ax, double Ay, double Az,
    double Bx, double By, double Bz,
    int typeA, int typeB,
    double& Sx, double D[3], double Q[6],    // Integrale (wie cgto_multipole_ints)
    double dS[3],                              // ∂S/∂R_A (Überlapp-Gradient)
    double ddA[3][3], double ddB[3][3],       // Dipol-Gradienten
    double dqA[3][6], double dqB[3][6]);      // Quadrupol-Gradienten
```

**Implementierungsstrategie**: Der TBLite-Port ist rein analytisch (Obara-Saika-ähnlich) — keine FD-Approximation. Die bestehende `cgto_multipole_ints()` kann als Ausgangspunkt dienen; die Gradiententerme entstehen durch explizite Ableitungen der Gaussschen Integrale nach dem Zentrum.

TBLite-Referenz-Struktur:
- `multipole_grad_cgto` ruft `overlap_grad_cgto` (für ∂S/∂R) und `multipole_cgto` (für D, Q) auf
- Die Gradiententerme folgen dem Reziprozitätsprinzip: `dD/dR_A = -dD/dR_B + dD/dR_C` (Schwerpunkt-Verschiebung)

### 5b.2 — Integration in Sektion 2b von `calculateGradient()`

In der bestehenden AO-Paar-Schleife (iao/jao), nach dem Überlapp-Gradienten:

```cpp
// In Sektion 2b, innerhalb der (iao, jao)-Schleife:
if (m_method == MethodType::GFN2 && m_mp_initialized) {
    double dS_mp[3], ddA[3][3], ddB[3][3], dqA[3][6], dqB[3][6];
    double Sx_mp, D_mp[3], Q_mp[6];
    cgto_multipole_grad(sh_a, sh_b, ..., typeA, typeB,
                        Sx_mp, D_mp, Q_mp, dS_mp, ddA, ddB, dqA, dqB);

    double mp_gx = 0.0, mp_gy = 0.0, mp_gz = 0.0;
    for (int k = 0; k < 3; ++k) {
        mp_gx -= Pmn * (ddA[0][k]*m_pot.v_dp(k,iat) + ddB[0][k]*m_pot.v_dp(k,jat));
        mp_gy -= Pmn * (ddA[1][k]*m_pot.v_dp(k,iat) + ddB[1][k]*m_pot.v_dp(k,jat));
        mp_gz -= Pmn * (ddA[2][k]*m_pot.v_dp(k,iat) + ddB[2][k]*m_pot.v_dp(k,jat));
    }
    for (int k = 0; k < 6; ++k) {
        mp_gx -= Pmn * (dqA[0][k]*m_pot.v_qp(k,iat) + dqB[0][k]*m_pot.v_qp(k,jat));
        mp_gy -= Pmn * (dqA[1][k]*m_pot.v_qp(k,iat) + dqB[1][k]*m_pot.v_qp(k,jat));
        mp_gz -= Pmn * (dqA[2][k]*m_pot.v_qp(k,iat) + dqB[2][k]*m_pot.v_qp(k,jat));
    }
    G_sval[0] += mp_gx; G_sval[1] += mp_gy; G_sval[2] += mp_gz;
}
```

**Hinweis:** `m_pot.v_dp` (3×nat) und `m_pot.v_qp` (6×nat) sind nach SCF-Konvergenz in `m_pot` gespeichert (in `xtb_native.cpp:Calculation()` vor `calculateGradient()` befüllt).

### 5b.3 — Validierung

**FD-Test** (`test_cases/test_xtb_gradient.cpp`):
```
Erwartung: max|G_anal - G_num| < 5e-4 Eh/Å  (GFN2, H₂O, CH₄, NH₃)
Vorher (nach AP5-S1 only): Fehler sichtbar durch fehlenden Pulay-Term
Nachher (AP5b): deutliche Verbesserung erwartet
```

**TBLite-Vergleich** (optional, USE_TBLITE):
```
max|G_native - G_tblite| < 1e-3 Eh/Bohr  (GFN2, H₂O, CH₄, NH₃)
```

---

## Akzeptanzkriterien

- [ ] `cgto_multipole_grad()` kompiliert und liefert konsistente Ergebnisse mit `cgto_multipole_ints()`
- [ ] FD-Test `xtb_gradient_H2O` (GFN2) passiert: max|ΔG| < 5e-4 Eh/Å
- [ ] FD-Test `xtb_gradient_CH4` (GFN2) passiert
- [ ] FD-Test `xtb_gradient_NH3` (GFN2) passiert
- [ ] `-opt H2O.xyz -method gfn2` konvergiert weiterhin (Regression)

## Fortschritt

| Datum | Aufgabe | Status | Notizen |
|-------|---------|--------|---------|
| 2026-04-26 | 5b.1 cgto_multipole_grad | Offen | TBLite Ref: multipole.f90:522–657 |
| 2026-04-26 | 5b.2 Integration Sektion 2b | Offen | — |
| 2026-04-26 | 5b.3 Validierung | Offen | test_xtb_gradient.cpp vorhanden |

## Schwierigkeiten / Blocker

- `cgto_multipole_grad` ist die komplexeste Integralfunktion im Projekt — der TBLite-Port erfordert sorgfältige Behandlung der Zentrumsderivative für kontrahierte GTOs
- d-Schalen weiterhin ausgelassen (`if (typeA < 0) continue`) — getrennt in AP7

---

## Referenzen

- TBLite: `external/tblite/src/tblite/integral/multipole.f90:522–657` (`multipole_grad_cgto`)
- TBLite: `external/tblite/src/tblite/xtb/h0.f90:444–447` (Verwendung im Gradient)
- Curcuma: `src/core/energy_calculators/qm_methods/xtb_multipole_ints.hpp` (cgto_multipole_ints als Basis)
- Curcuma: `src/core/energy_calculators/qm_methods/xtb_gradient.cpp` Sektion 2b (Integrationspunkt)
