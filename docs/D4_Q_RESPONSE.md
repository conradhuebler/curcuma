# D4 Dispersion Charge Response — ∂E_D4/∂q · ∂q/∂x

**Status:** Phase 1 + 2 abgeschlossen und FD-validiert; Phase 3a abgeschlossen
(Mulliken-Quelle schaltbar); Phase 3b (analytischer CPSCF/Z-Vektor) offen.
**Vorgängerarbeit:** Commit `61f5846` (GFN2 D4-Dispersion + analytischer
Gradient via `D4Evaluator`). Dort wurde `zetac6` als statischer Prefaktor
behandelt — der Ladungs-Kettenregel-Term fehlte.

> ⚠️ AI-generiert, maschinell (FD) getestet. Menschliche Produktionstests stehen aus.

## 1. Motivation

Die D4-Paarenergie ist

```
E_pair = -ζc6 · C6 · ( s6·t6 + s8·r4r2·t8 ),   ζc6 = ζ(Z_i,q_i)·ζ(Z_j,q_j)
```

`ζc6` hängt über die Atomladungen `q` von der Geometrie ab. Der vollständige
Gradient enthält daher neben dem radialen und dem CN-Term auch

```
∂E_D4/∂x  ⊇  Σ_A (∂E_D4/∂q_A) · (∂q_A/∂x)
              └── D4Evaluator (Phase 1) ──┘  └── Ladungsmodell (Phase 2/3) ──┘
```

Vorher wurde `ζc6` eingefroren (statischer Prefaktor) → Sub-mEh/Bohr-Residual.
Für die finale TBLite-Übereinstimmung und saubere Geometrieoptimierung wird der
Term gebraucht.

## 2. Architektur

Der Term zerlegt sich in einen **methodenunabhängigen** Faktor `∂E_D4/∂q`
(im `D4Evaluator`) und einen **ladungsquellen-abhängigen** Faktor `∂q/∂x`.
Schalter: `d4_charge_source = "eeq" (Default) | "mulliken"`.

| Quelle | Ladungen für zetac6 | ∂q/∂x | Status |
|---|---|---|---|
| `eeq` | kanonische Einzelschritt-dftd4-EEQ | analytisch (Adjoint/Z-Vektor) | ✅ Phase 2 |
| `mulliken` | GFN2-SCF-Mulliken-Ladungen | CPSCF/Z-Vektor auf GFN2-Fock | 🔜 Phase 3b (Energie+∂E/∂q: ✅ 3a) |

GFN-FF nutzt weiterhin Zweiphasen-Topologie-Ladungen und den statischen
Prefaktor (deren ∂q/∂x ist nicht glatt differenzierbar).

## 3. Phase 1 — ∂E_D4/∂q (abgeschlossen)

- `GFNFFParameters::zetaChargeScaleDerivative(Z, q)` (`gfnff_par.h`): Closed-form
  `dζ/dq = -3·ζ·g·c·zeff/qmod²` (in geklammerten Bereichen 0).
- `D4Evaluator::computeEnergyAndGradient(..., dEdq_out, with_dEdq)`
  (`d4_evaluator.{h,cpp}`): akkumuliert pro Atom
  `∂E/∂q_i = -C6·disp_sum·(dζ_i)·ζ_j` (und symmetrisch).
- Nur `ζc6` ist ladungsabhängig (CN-only-C6-Gewichtung nicht) → exakt additiv.
- **Validierung** (`test_d4_dedq`): FD von `[E(q+δ)-E(q-δ)]/2δ` vs analytisch,
  MaxErr **1e-14 Eh/e** auf H₂O/CH₄/C₆H₆.

## 4. Phase 2 — ∂q_EEQ/∂x analytisch (abgeschlossen)

Neues `d4_charge_model.{h,cpp}`: kanonisches Einzelschritt-dftd4-EEQ (ein
glattes Linearsystem statt des GFN-FF-Zweiphasen-Solvers, dessen Floyd-Warshall-
Distanzen + diskrete Korrekturen nicht differenzierbar sind):

```
A_ii = γ_i + sqrt(2/π)/sqrt(α_i)
A_ij = erf(γ_ij·r_ij)/r_ij,   γ_ij = 1/sqrt(α_i+α_j)
b_i  = -χ_i + κ_i·sqrt(CN_i)
[ A  1 ] [ q ]   [ b ]
[ 1  0 ] [ λ ] = [ Q ]
```

Parameter: angewChem2020-Satz (`chi/gam/alpha/cnf_eeq`, α quadriert); CN: GFN-FF
log-komprimierte erf-CN. **Adjoint/Z-Vektor**: löse `M·z = [dEdq; 0]` einmal,
dann kontrahiere `z` gegen die Closed-form-Geometrieableitungen von `A` (erf-
Kernel) und `b` (über die CN-Ableitung).

- GFN2-Default (`eeq`) nutzt diese Ladungen für zetac6 (D4-Energie für
  H₂O/CH₄/C₆H₆ praktisch unverändert vs Zweiphasen).
- q-Response in `m_disp_gradient` gefaltet (`xtb_native.cpp::calcDispersionEnergy`).
- **Validierung** (`test_d4_dedq`): direkte FD der Einzelschritt-EEQ-Ladungen vs
  analytisch, MaxErr **~1e-11 Eh/Bohr**; voller GFN2-FD-Gradient **< 5e-5 Eh/Å**
  (`test_xtb_gradient`).

## 5. Phase 3a — Mulliken-Quelle (abgeschlossen)

- `XTB::setD4ChargeSource("mulliken")`; CLI-Flag `-d4_charge_source mulliken`
  (als PARAM im `xtb`-Scope registriert → `EnergyCalculator::reattachMethodScopes`
  → `GFN2Method`).
- In `calcDispersionEnergy`: bei `mulliken` werden `m_wfn.q_at` (SCF-Mulliken)
  über `setTopologyCharges` als zeta-Ladungen eingespeist; `dEdq` (Phase 1) wird
  damit w.r.t. der Mulliken-Ladungen gebildet.
- **Verifiziert**: aktiviert korrekt (H₂O: q(O)=-0.564, D4 -0.000163 vs eeq
  -0.000164 Eh).
- **Einschränkung**: der **Gradient-q-Response fehlt** (Phase 3b) — mulliken
  nutzt bis dahin den statischen Prefaktor. Für Gradienten `eeq` bevorzugen.

## 6. Roadmap — Phase 3b: ∂q_Mulliken/∂x via CPSCF/Z-Vektor

Ziel: analytischer Ladungs-Response der **GFN2-SCF-Mulliken-Ladungen**.

### 6.1 Mathematik (closed-shell, nicht-orthogonal)

Zielgröße `L = Σ_A w_A q_A`, `w_A = ∂E_D4/∂q_A` (= `m_disp_dEdq`).
Mit `q_A = n0_A − Σ_{μ∈A}(P S)_{μμ}` und `Λw = diag(w_{atom(μ)})`:

```
G = ∂L/∂x = 2·Σ_ij S̃ˣ_ji Q̃_ji  +  4·Σ_ia z_ai B̃ˣ_ai  −  Tr(P Sˣ Λw)
```

- `Q_sym = sym(S·Λw)`, `Q̃ = Cᵀ Q_sym C` (MO-Basis), `S̃ˣ = Cᵀ Sˣ C`.
- Z-Vektor: löse den Orbital-Hessian **einmal** `A·z = Q̃_ov` (occ-virt),
  dann `Σ_ia U^x_ai Q̃_ai = −Σ_ia z_ai B̃ˣ_ai`.
- `A·u`-Anwendung = `(ε_a−ε_i)u_ai + (Cᵀ δF[u] C)_ai`, wobei `δF` aus der
  SCC-Kernel-Antwort auf die von `u` erzeugte Dichteänderung kommt.
- `B̃ˣ` (CPHF-Störung) = occ-virt-Block von `H0ˣ + Vˣ − ε_i Sˣ` (+ occ-occ-
  Relaxations-Kernelterm aus `S̃ˣ`).

### 6.2 Implementierungsschritte

1. **Orbital-Hessian-Anwendung** (`A·u`): neue Routine, die aus einer occ-virt-
   Rotation `u` die Dichteänderung `δP`, daraus `δq_sh`/`δq_at`/`δ`-Multipole und
   via bestehender Kernel (`m_gamma`, 3.Ordnung, Multipol) das `δV`/`δF` bildet.
   → Refaktorierung der Potential-Routinen (`addCoulombShellPotential`,
   `addThirdOrderPotential`, `addMultipolePotential`) auf **beliebige**
   Ladungs-/Multipol-Variationen statt `m_wfn`-Member.
2. **Z-Vektor-Solver**: CG im occ-virt-Raum (Hessian SPD nahe Konvergenz),
   Diagonal-Preconditioner `(ε_a−ε_i)⁻¹`.
3. **Perturbations-/Pulay-Terme** (`B̃ˣ`, `S̃ˣ`, `Tr(P Sˣ Λw)`): Kontraktion der
   **bestehenden analytischen** Integralableitungen (`dS/dR`, `dH0/dR` in
   `xtb_gradient.cpp`) mit der **relaxierten Z-Dichte** statt mit `P`/`W`.
   → `calculateGradient` so refaktorieren, dass die Integralableitungs-
   Kontraktionen `(Dichte, energiegewichtete Dichte, Ladungen)` als Parameter
   nehmen; dann zweimal aufrufen (Energiegradient + Response).
4. **Verkabelung**: in `calcDispersionEnergy` bei `mulliken` `m_disp_dEdq` an die
   neue `computeMullikenChargeResponse(dEdq, m_disp_gradient)` übergeben.

### 6.3 Validierung (Akzeptanz)

- Unit: `A·u` vs FD der statischen Ladungsantwort auf eine Potentialstörung.
- `test_xtb_gradient --d4_charge_source mulliken`: FD **< 5e-5 Eh/Å** auf
  H₂O/CH₄/NH₃/C₆H₆.
- Cross-Check: `mulliken` vs `eeq` Energiedifferenz < 1 mEh; `mulliken` näher an
  TBLite-GFN2 (TBLite nutzt intern Mulliken).

### 6.4 Risiken

- Faktoren (Spin ×2, Symmetrie ×2), Metrik (nicht-orthogonal `S̃`),
  Vorzeichen — FD-Validierung jedes Teilschritts zwingend.
- Multipol-Kernel im Orbital-Hessian ist der heikelste Teil; ggf. zuerst
  isotrop (γ + 3.Ordnung) auf einem Molekül mit vernachlässigbaren Multipolen
  validieren, dann Multipol ergänzen.
- Umfang ~größer als Phase 1+2 zusammen — eigener AP.

## 7. Dateien

| Datei | Rolle |
|---|---|
| `dispersion/d4_evaluator.{h,cpp}` | `∂E/∂q` (Phase 1) |
| `dispersion/d4param_generator.{h,cpp}` | `getZeta(Derivative)`, `getZetaCharges`, Single-Shot-Hook |
| `dispersion/d4_charge_model.{h,cpp}` | Einzelschritt-EEQ + analytisches `∂q/∂x` (Phase 2) |
| `ff_methods/gfnff_par.h` | `zetaChargeScaleDerivative` |
| `qm_methods/xtb_native.{h,cpp}` | `d4_charge_source`, Mulliken-Quelle (3a), Response-Faltung |
| `qm_methods/xtbinterface.h` | PARAM `d4_charge_source` (CLI-Routing) |
| `qm_methods/gfn2_method.cpp` | liest `d4_charge_source` → `setD4ChargeSource` |
| `qm_methods/xtb_response.{h,cpp}` | **(Phase 3b, neu)** CPSCF/Z-Vektor |
| `test_cases/test_d4_dedq.cpp` | FD-Validierung Phase 1+2 |

## 8. Referenzen

- E. Caldeweyher et al., *J. Chem. Phys.* **150**, 154122 (2019) — D4 + EEQ.
- C. Bannwarth, S. Ehlert, S. Grimme, *JCTC* **15**, 1652 (2019) — GFN2-xTB.
- N. C. Handy, H. F. Schaefer, *J. Chem. Phys.* **81**, 5031 (1984) — Z-Vektor.
