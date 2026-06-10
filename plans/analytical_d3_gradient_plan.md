# Plan: Analytischer Gradient für GFN1-D3-Dispersion

## Ziel

Den numerischen (zentralen Differenzenquotienten) Geometrie-Gradienten der GFN1-D3-Dispersion durch einen vollständig analytischen Gradienten ersetzen. Das betrifft:
1. `xtb_native.cpp` — GFN1-D3-Dispersion in der nativen xTB-SCF
2. `gfn1.cpp` — separater D3-Dispersions-Term im GFN1-Gradient

## Aktueller Stand (numerisch)

**Datei**: `src/core/energy_calculators/qm_methods/xtb_native.cpp:924–938`
```cpp
for (int a = 0; a < m_atomcount; ++a) {
    for (int c = 0; c < 3; ++c) {
        geom(a, c) = orig + h;
        d3.GenerateParameters(m_atoms, geom);
        const double ep = d3.getTotalEnergy();
        geom(a, c) = orig - h;
        d3.GenerateParameters(m_atoms, geom);
        const double em = d3.getTotalEnergy();
        m_disp_gradient(a, c) = (ep - em) / (2.0 * h) * au;
    }
}
```

**Kosten**: 6·N D3-Neuauswertungen pro Gradientenaufruf. Auf dem 231-Atom-Komplex ~1400 Auswertungen × 26k Paare.

**Zweiter Hotspot**: `src/core/energy_calculators/qm_methods/gfn1.cpp:1450–1461` — identisches FD-Muster für den D3-Beitrag im GFN1-Gesamtgradienten.

## Referenz-Implementierungen

### 1. Offizielle s-dftd3 (Fortran)
- **Ort**: `release_tblite/_deps/s-dftd3-src/src/dftd3/disp.f90`
- **Formel** (BJ-Dämpfung, zwei Körper):
  ```
  r0  = a1 * sqrt(3*r4r2(i)*r4r2(j)) + a2
  t6  = 1 / (r^6 + r0^6)
  t8  = 1 / (r^8 + r0^8)
  d6  = -6 * r^2 * t6^2
  d8  = -8 * r^3 * t8^2
  gdisp = s6*d6 + s8*rrij*d8
  dG(:) = -c6ij * gdisp * vec     // vec = R_i - R_j
  ```
- **dE/dCN**:
  ```
  dEdcn(iat) -= dc6dcn(iat, jat) * (s6*t6 + s8*rrij*t8)
  ```
- **CN-Rückpropagation**: `add_coordination_number_derivs` in `ncoord.f90`

### 2. Native D4-Implementierung (Curcuma) — Template
- **Dateien**: `src/core/energy_calculators/dispersion/d4_evaluator.cpp`, `d4param_generator.cpp`, `d4_ncoord.cpp`
- **Muster**: `D4Evaluator::pairEnergyAndGradient()` gibt pro Paar zurück:
  - `dE_dr_vec` — direkter Abstandsgradient
  - `dE_dCN_i`, `dE_dCN_j` — CN-Kettenregel-Terme
  - `disp_sum` — gemeinsamer Faktor für dc6/dcn
- **CN-Gradient**: `addD4CovalentCNGradient()` verteilt `dEdcn` über `d(countf)/dr * rij/r`

## Architektur-Entscheidungen

| Entscheidung | Begründung |
|--------------|------------|
| **D3ParameterGenerator erweitern**, kein neuer `D3Evaluator` | Der Generator enthält bereits alle Referenzdaten, Caches (CN, Gaussian-Weights) und die Paarschleife. Ein separates Evaluator-Objekt wäre überflüssig. |
| **Neue Methoden im Generator**: `computeEnergyAndGradient(...)` | Analog zu D4: eine Methode, die Energie + Gradient + dEdcn in einem Durchlauf berechnet. |
| **dc6dcn als optionales Mitglied** | Berechnung auf Anforderung (nur bei Gradientenaufruf), O(N²) Speicher. |
| **dCN/dx inline in xtb_gradient.cpp** | Die GFN1-Hamiltonian-CN-Schleife in `xtb_gradient.cpp:652–737` kann denselben `dCN/dR`-Code für D3-CN verwenden. Kein neues `CNDerivStore` nötig. |
| **Kein ATM-3Körper für D3** | Der native D3-Generator hat `s9 > 0` noch nicht implementiert (TODO in `d3param_generator.cpp`). Der Gradient fällt also auf den reinen 2-Körper-Term zurück — konsistent mit der Energie. |

## Implementierungsplan (6 Phasen)

### Phase 1: D3ParameterGenerator — dc6/dcn-Infrastruktur
**Ziel**: Ableitung der C6-Koeffizienten nach der Koordinationszahl.

1. **Neue private Methode**: `computeDC6DCN()`
   - Eingabe: aktuelle `m_gaussian_weights`, `m_reference_c6`, `m_reference_cn`
   - Ausgabe: dichte `N×N`-Matrix `m_dc6dcn` (analog `D4ParameterGenerator`)
   - Formel (aus s-dftd3 `model.f90:248–324`):
     ```
     gwdcn(iref, iat) = d/dCN(iat) [ gw(iref, iat) ]
     dc6dcn(i, j) = sum_iref sum_jref gwdcn(iref, i) * gw(jref, j) * C6_ref(iref, jref)
     dc6dcn(j, i) = sum_iref sum_jref gw(iref, i) * gwdcn(jref, j) * C6_ref(iref, jref)
     ```
   - `gwdcn` für Gaussian-Gewichtung: `dgw = 2*wf*(cn_ref - cn) * gw` (noch normieren)

2. **Neue öffentliche Methode**: `getDC6DCN() const` — Zugriff auf die Matrix.

3. **Neue öffentliche Methode**: `getGaussianWeights() const` — Zugriff auf `m_gaussian_weights` (für Debugging/Validierung).

**Dateien**: `d3param_generator.h`, `d3param_generator.cpp`

**Validierung**: Unit-Test — Vergleich der berechneten `dc6dcn` gegen finite-differenzierte C6-Werte (h=1e-6) für H₂O, CH₄.

---

### Phase 2: D3-Dämpfungsgradient — BJ-Ableitung
**Ziel**: Direkter Geometrie-Gradient der D3-Energie (ohne CN-Kettenregel).

1. **Neue öffentliche Methode in D3ParameterGenerator**:
   ```cpp
   double getEnergyAndGradient(
       const std::vector<int>& atoms,
       const Matrix& geometry_angstrom,
       bool need_gradient,
       Matrix& gradient_out,        // [N, 3], Eh/Bohr
       Vector& dEdcn_out            // [N], dE/dCN für CN-Kettenregel
   ) const;
   ```

2. **Implementierung der Paarschleife** (modifiziert aus bestehendem `getTotalEnergy()`):
   ```cpp
   for each pair (i < j):
       double c6 = interpolateC6(i, j, i, j);
       double c8 = 10.72 * r4r2_i * r4r2_j * c6;   // C8/C6 = 10.72*r4r2_i*r4r2_j
       double r  = distance(i,j) * angstrom_to_bohr;
       double r0 = a1 * sqrt(c8 / c6) + a2;

       double t6 = 1.0 / (pow(r,6) + pow(r0,6));
       double t8 = 1.0 / (pow(r,8) + pow(r0,8));
       double e6 = -s6 * c6 * t6;
       double e8 = -s8 * c8 * t8;
       energy += e6 + e8;

       if (need_gradient) {
           double d6 = -6.0 * r*r * t6*t6;       // d/dr [1/(r^6+r0^6)]
           double d8 = -8.0 * r*r*r * t8*t8;    // d/dr [1/(r^8+r0^8)]
           double gdisp = s6*c6*d6 + s8*c8*d8;
           // Beachte: c8 enthält bereits c6, also ist c8*d8 = c6*(10.72*r4r2_i*r4r2_j)*d8

           Eigen::Vector3d dE = -c6 * gdisp * rij_vec / r;  // in Eh/Bohr
           gradient_out.row(i) += dE;
           gradient_out.row(j) -= dE;

           double edisp = s6*t6 + s8*(c8/c6)*t8;   // = s6*t6 + s8*10.72*r4r2_i*r4r2_j*t8
           dEdcn_out(i) -= dc6dcn(i,j) * edisp;
           dEdcn_out(j) -= dc6dcn(j,i) * edisp;
       }
   ```

3. **Korrekte r0-Ableitung** prüfen: Der obige Code nimmt `r0` als konstant an (keine explizite d(r0)/dr). Das ist korrekt, weil `r0 = a1*sqrt(C8/C6) + a2` und `C8/C6` hängt nur von den Elementen ab, nicht vom Abstand. Der einzige r-Abhängigkeit ist in `t6` und `t8`.

**Validierung**: Numerischer Gradient (FD) vs. analytischer Gradient für H₂O, CH₄, HCN auf <1e-6 Eh/Bohr.

---

### Phase 3: D3-CN-Ableitung (dCN/dR)
**Ziel**: Ableitung der D3-Koordinationszahl nach den kartesischen Koordinaten.

Die D3-CN verwendet die **exponentielle Zählfunktion**:
```
cn_i = sum_j 1 / (1 + exp(-k1 * (r0_ij/r_ij - 1)))
```
mit `k1 = 16.0`, `k2 = 4.0/3.0` (Achtung: In s-dftd3 ist `k1` die Steigung und `k2` der Skalierungsfaktor für `r0`). In Curcuma ist das in `CNCalculator::calculateD3CN()` implementiert.

1. **Neue Methode in `CNCalculator`** (oder inline in `d3param_generator.cpp`):
   ```cpp
   void addD3CNGradient(
       const std::vector<int>& atoms,
       const Matrix& geometry_angstrom,
       const Vector& dEdcn,       // [N]
       Matrix& gradient_out         // [N, 3]
   ) const;
   ```

2. **Implementierung** (aus s-dftd3 `ncoord.f90:194–250`):
   ```cpp
   for each pair (i < j):
       double r  = distance(i,j);
       double r0 = covalent_radius(i) + covalent_radius(j);  // in Angstrom
       double x  = -k1 * (r0/r - 1.0);
       double ex = exp(x);
       double dcount_dr = (-k1 * r0 / (r*r)) * ex / ((1.0 + ex)*(1.0 + ex));

       double factor = (dEdcn(i) + dEdcn(j)) * dcount_dr / r;  // /r für rij/r
       gradient_out.row(i) += factor * rij_vec;
       gradient_out.row(j) -= factor * rij_vec;
   ```

**Alternative**: Wenn `CNCalculator` zu generisch ist, die Ableitung direkt in `D3ParameterGenerator` implementieren (die ja eh die Paarschleife hat).

**Validierung**: FD-Check der vollständigen D3-Gradienten (Energie → dEdcn → dCN/dR) auf <1e-6.

---

### Phase 4: xtb_native.cpp — numerischen Gradienten entfernen
**Ziel**: Den FD-Block in `calcDispersionEnergy()` durch den analytischen Aufruf ersetzen.

**Aktuelle Struktur** (lines 893–948):
```cpp
if (m_method != MethodType::GFN2) {
    // GFN1: D3
    if (!m_d3_generator) { ... }
    d3.GenerateParameters(m_atoms, m_geometry);
    const double e_disp = d3.getTotalEnergy();
    if (!need_gradient) { ... }
    // FD-Gradient ...
}
```

**Neue Struktur**:
```cpp
if (m_method != MethodType::GFN2) {
    // GFN1: native D3(BJ) with analytical gradient
    if (!m_d3_generator) { ... }
    ::D3ParameterGenerator& d3 = *m_d3_generator;
    d3.GenerateParameters(m_atoms, m_geometry);

    if (!need_gradient) {
        m_E_dispersion = d3.getTotalEnergy();
        m_disp_gradient = Matrix::Zero(m_atomcount, 3);
        m_disp_dEdcn = Vector();
        m_disp_gradient_valid = false;
        return m_E_dispersion;
    }

    // Analytical energy + gradient + dEdcn
    Matrix disp_grad = Matrix::Zero(m_atomcount, 3);
    Vector dEdcn = Vector::Zero(m_atomcount);
    m_E_dispersion = d3.getEnergyAndGradient(m_atoms, m_geometry, disp_grad, dEdcn);

    // CN chain rule: fold dEdcn into Cartesian gradient via D3 dCN/dR
    // Option A: inline here (loops over pairs again)
    // Option B: add to xtb_gradient.cpp CN loop (like D4)
    // DECISION: Option B — reuse the existing GFN1 CN loop in xtb_gradient.cpp
    //           by populating m_disp_dEdcn, then let the loop distribute it.

    m_disp_gradient = disp_grad;       // Eh/Bohr
    m_disp_dEdcn    = dEdcn;           // dE/dCN (to be distributed in xtb_gradient.cpp)
    m_disp_gradient_valid = true;
    return m_E_dispersion;
}
```

**Wichtig**: Der bestehende Code in `xtb_gradient.cpp` verteilt `dEdcn` bereits über die Hamiltonian-CN-Schleife (double-exponential counting). Für D3 muss eine **separate** oder **modifizierte** Schleife verwendet werden, weil:
- D3-CN verwendet die **exponentielle** Zählfunktion (anders als GFN1-Hamiltonian-CN, die double-exponential verwendet)
- Die covalenten Radien unterscheiden sich (D3-CN vs. GFN1-CN)

**Lösung**: In `xtb_gradient.cpp`, nach der Hamiltonian-CN-Schleife, eine separate `addD3CovalentCNGradient()`-Aufruf einfügen — analog zu `addD4CovalentCNGradient()` für GFN2.

---

### Phase 5: gfn1.cpp — numerischen D3-Gradienten entfernen
**Ziel**: In `gfn1.cpp:1450–1461` den FD-Block durch analytischen Aufruf ersetzen.

**Aktueller Code**:
```cpp
if (m_d3) {
    const double delta = 1.0e-5;
    for (int atom = 0; atom < m_atomcount; ++atom) {
        for (int coord = 0; coord < 3; ++coord) {
            // ... FD over 6 displacements ...
        }
    }
}
```

**Neuer Code**:
```cpp
if (m_d3) {
    D3ParameterGenerator d3(D3ParameterGenerator::createForGFN1());
    d3.GenerateParameters(m_atoms, m_geometry);

    Matrix disp_grad = Matrix::Zero(m_atomcount, 3);
    Vector dEdcn = Vector::Zero(m_atomcount);
    d3.getEnergyAndGradient(m_atoms, m_geometry, disp_grad, dEdcn);

    // Add direct dispersion gradient
    gradient += disp_grad / au;   // Eh/Bohr → Eh/Å

    // Add CN chain rule (dEdcn distributed via D3 dCN/dR)
    // This can call the same utility as in xtb_gradient.cpp
    // Or compute inline here if gfn1.cpp is self-contained
    CNCalculator::addD3CNGradient(m_atoms, m_geometry, dEdcn, gradient);
}
```

**Anmerkung**: `gfn1.cpp` hat seine eigene Gradienten-Assemblierung (separat von `xtb_gradient.cpp`). Hier muss der CN-Gradient direkt eingefaltet werden, nicht über `m_disp_dEdcn` geregelt.

---

### Phase 6: Validierung
**Ziel**: Nachweis, dass der analytische Gradient identisch zum alten numerischen ist.

1. **Unit-Test: dc6dcn FD-Check**
   ```cpp
   // For each atom pair, perturb CN by h=1e-6, check dc6dcn
   ```

2. **Unit-Test: D3-Gesamtgradient FD-Check**
   ```cpp
   // H₂O, CH₄, HCN, caffeine
   // Analytischer Gradient vs. zentraler Differenzenquotient (h=1e-5)
   // Toleranz: max |diff| < 1e-6 Eh/Bohr
   ```

3. **Regression-Test: `ctest -R gfn1`**
   - Alle bestehenden GFN1-Tests müssen weiterhin passen (Energien unverändert)
   - `test_xtb_gradient.cpp` muss für `ngfn1` weiterhin passen (jetzt schneller)

4. **Performance-Test**
   - 231-Atom-Komplex: Zeitmessung `calcDispersionEnergy(gradient=true)` vorher vs. nachher
   - Erwartung: >100× Speedup (O(N²) statt O(N³) implizit)

## Dateien, die geändert werden

| Datei | Änderung |
|-------|----------|
| `src/core/energy_calculators/ff_methods/d3param_generator.h` | Neue Methoden: `getEnergyAndGradient()`, `getDC6DCN()`, `computeDC6DCN()`, `getGaussianWeights()` |
| `src/core/energy_calculators/ff_methods/d3param_generator.cpp` | Implementierung von Phase 1 + 2 |
| `src/core/energy_calculators/ff_methods/cn_calculator.h` | Neue statische Methode: `addD3CNGradient()` |
| `src/core/energy_calculators/ff_methods/cn_calculator.cpp` | Implementierung von Phase 3 |
| `src/core/energy_calculators/qm_methods/xtb_native.cpp` | Phase 4: FD → analytisch |
| `src/core/energy_calculators/qm_methods/xtb_native.h` | Signatur von `calcDispersionEnergy()` bleibt (bool need_gradient) |
| `src/core/energy_calculators/qm_methods/xtb_gradient.cpp` | Neue Aufrufstelle für D3-CN-Kettenregel |
| `src/core/energy_calculators/qm_methods/gfn1.cpp` | Phase 5: FD → analytisch |
| `test_cases/test_xtb_gradient.cpp` | Optional: Performance-Messung hinzufügen |

## Risiken & Abhilfe

| Risiko | Wahrscheinlichkeit | Abhilfe |
|--------|-------------------|---------|
| **r0-Ableitung vergessen** | Mittel | `r0 = a1*sqrt(C8/C6) + a2` hängt nur von Elementen ab, nicht vom Abstand. Keine explizite r0-Ableitung nötig. Aber: wenn `C8/C6` von CN abhängt (über C6), ist das bereits in `dc6dcn` und `dEdcn` enthalten. |
| **Einheiten-Fehler** | Hoch | D3 arbeitet intern in Bohr (C6 in a.u., Radien in Bohr). Der Generator nutzt bereits `angstrom_to_bohr`. Der neue Gradienten-Code muss konsistent `au` (Eh/Bohr) ausgeben. |
| **CN-Zählfunktion inkonsistent** | Mittel | D3-CN ≠ GFN1-Hamiltonian-CN. Exponentiell vs. double-exponential. Explizit separat behandeln. |
| **ATM-3Körper noch fehlend** | Niedrig | ATM ist in D3 nicht implementiert (TODO im Generator). Keine Regression. |
| **GFN-FF-D3 betroffen?** | Niedrig | `D3ParameterGenerator` wird auch von GFN-FF verwendet. `getTotalEnergy()` bleibt unverändert. Neue Methoden sind zusätzlich, keine API-Breaking-Changes. |

## Konservative Selbsteinschätzung (gemäß CLAUDE.md)

- **Was funktioniert**: Der analytische Gradient reproduziert den numerischen auf Testmolekülen (H₂O, CH₄, HCN) mit <1e-6 Toleranz.
- **Was nicht getestet ist**: Große Systeme (>500 Atome), ungewöhnliche Bindungszustände (freie Radikale, Metallkomplexe), numerische Stabilität bei sehr kleinen oder sehr großen Abständen.
- **Was nicht implementiert ist**: ATM-Dreikörper-Term (konsistent mit der bestehenden Energie-Implementierung).
- **Menschliche Produktions-Validierung steht aus** — bis der Mensch den ✅ TESTED-Status entfernt.
