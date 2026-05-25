# WP Phase 3b-4 (Rest) — Multipol-Interaktions-Linearisierung + expliziter Overlap-Term

**Status:** Offen (3b-1/3b-2/3b-3 erledigt, 3b-4 teilweise: Multipol-Integral-Pulay drin)
**Erstellt:** 2026-05-25
**Branch:** `feature/sqm-implementation`
**Build:** immer in `release/` (`cd release && make -j4`)

---

## Ziel

GFN2-D4 mit `d4_charge_source="mulliken"`: voller FD-Gradient **< 5e-5 Eh/Å** auf
H₂O/CH₄/NH₃/C₆H₆. Aktuell fehlt für stark-multipolare Moleküle (H₂O) die
**Multipol-Interaktions-Linearisierung** + der **explizite Overlap-Term**.

---

## Was bereits funktioniert (NICHT erneut tun)

| Phase | Inhalt | Validierung |
|-------|--------|-------------|
| 3b-1 | Parametrisierte Potential-Overloads (`addCoulombShellPotential(pot,q_sh)`, `addThirdOrderPotential(pot,q_sh,q_at)`, `addMultipolePotential(pot,q_at,dp,qp)`, `thirdOrderKernelDiag`, `buildFockFromPotential`, `mullikenFromDensity`) | ngfn2-Energien bit-identisch |
| 3b-2 | `applyOrbitalHessian` (Orbital-Hessian A·z, **inkl. Multipol-Kernel**) + `solveZVector` (PCG, relatives Kriterium) in `xtb_response.cpp` | `test_xtb_cpscf` A=1.9e-12, B: H₂O 1.6e-16 / CH₄ 3.3e-12 |
| 3b-3 | `computeMullikenChargeResponse`: isotroper Response (Pulay D_z/W_z, Coulomb δq_z linearisiert, CN-Kettenregel, Z-Vektor) | H₂ Response **exakt 0**; HCN ratio 1.105 |
| 3b-4a | Multipol-Integral-Pulay (AP5b-Analog, D_z statt P) im Pulay-Loop | HCN ratio → **1.044** |

**Schlüsselerkenntnis (gilt weiter):** Der Response ist die **Linearisierung des
Energie-Gradienten**: P→D_z (Pulay sval), W→W_z, und jedes Ladungs-/Momentprodukt
wird linearisiert (Coulomb: `q_is·q_js → δq_is·q_js + q_is·δq_js`).

**Der CPSCF-Kernel ist vollständig** (Multipole sind drin → Test B maschinengenau).
`D_z` ist also voll korrekt. Es fehlt **nur** die Gradient-Assemblierung der
Multipol-Interaktionsterme.

---

## Offene Arbeit

### (a) Multipol-Interaktions-Linearisierung — der H₂O-Fix

**Referenz:** `xtb_gradient.cpp` Section 5, Z. **523–649**
(`get_multipole_gradient_0d`-Port: SD/DD/SQ-Paarterme + mrad/CN-Kettenregel).

**Rezept:** Section 5 berechnet `dg = ∂E_mp/∂R` bei festen Momenten. Der Response
ist `∂/∂R` der **in den Momenten linearisierten** E_mp, d.h. Section 5 mit
`(q_at, dp_at, qp_at) → (q_at+δq, dp_at+δdp, qp_at+δqp)`, linearer Anteil in δ.
Konkret: in **jedem** bilinearen Term genau **einen** Momentfaktor durch sein δ
ersetzen und über alle Ersetzungen summieren. Beispiel SD:

```
Energie:  dpiqj = (vec·dp_iat)·q_jat
Response: δ(dpiqj) = (vec·δdp_iat)·q_jat + (vec·dp_iat)·δq_jat
```

analog DD (dp·dp → δdp·dp + dp·δdp) und SQ (qp·q → δqp·q + qp·δq).

**δ-Momente liegen bereits vor** in `computeMullikenChargeResponse`:
```cpp
Vector dq_sh, dq_at; Eigen::MatrixXd ddp, dqp;
mullikenFromDensity(D_z, dq_sh, dq_at, ddp, dqp);  // δq_at=dq_at, δdp=ddp, δqp=dqp
```
(`dq_at` ist nat-lang; `ddp` 3×nat; `dqp` 6×nat — gleiche Layout-Konvention wie
`m_wfn.dp_at`/`m_wfn.qp_at`.)

**mrad/CN-Kettenregel:** Section 5 akkumuliert `dEdr_mp` → `dEdcn` (Z.~647). Für den
Response analog mit den linearisierten Energiebeiträgen; in das **vorhandene**
`dEdcn`-Akkumulator in `computeMullikenChargeResponse` einspeisen (wird in Abschnitt 8
verteilt). Achtung: mrad ist rein geometrieabhängig (nicht momentabhängig), daher
kommt die δ-Variation nur über die Momentfaktoren, nicht über mrad selbst.

**Erwartung:** liefert einen **negativen** Beitrag, der H₂Os isotropen Überschuss
(`|Δanal|≈6.6e-5`) auf den FD-Wert (`|Δnum|≈2.5e-5`) zieht. HCN sollte nahe 1.0
bleiben.

### (b) Expliziter Overlap-Term −Tr(Λw·P·Sˣ)

In `computeMullikenChargeResponse` als `sval += −EXPL_FAC·lam_pair·Pmn` bereits
verdrahtet, aber `EXPL_FAC=0` (Hinzufügen allein verschlechterte die FD). Gehört zur
Sˣ-Familie wie (a) → **gemeinsam** mit (a) FD-kalibrieren. Meine Ableitung ergab
`−½(Λw_μ+Λw_ν)·P_μν` (also EXPL_FAC=0.5 mit `lam_pair=dEdq(iat)+dEdq(jat)`), aber
empirisch muss das gegen FD geprüft werden, sobald (a) drin ist.

---

## Knöpfe (oben in `computeMullikenChargeResponse`, `xtb_response.cpp`)

```cpp
constexpr double RHS_SIGN = -1.0;   // FIX (durch HCN-Vorzeichen/Skala validiert)
constexpr double EXPL_FAC = 0.0;    // (b): kalibrieren, vermutlich 0.5 oder ±1
```
`RHS_SIGN=-1` nicht ändern. Nach (a) ggf. ein globaler Skalenknopf für den
Multipol-Interaktionsterm, falls Faktor abweicht (analog zu D_z=2·(...)-Konvention).

---

## Validierung / Test

```bash
cd release && make -j4 test_xtb_cpscf
./test_cases/test_xtb_cpscf
```

**Isolations-Methode (entscheidend):** Es gibt **keinen** neutralen multipolfreien
Test mit nonzero-Response (H₂ ist homonuklear → q≡0). Daher wird der Response gegen
den **eeq-Pfad differenziert** (`testC_responseIsolation`): der vorbestehende
GFN2-Gradient-Baseline-Fehler hebt sich in `(mulliken−eeq)` weg.
```
d_anal = mull_anal − eeq_anal,  d_num = mull_FD − eeq_FD
err = max|d_anal − d_num|,  ratio = |Δanal|/|Δnum|
```
**Ziel:** ratio → 1.0 und err < 5e-5 für H₂O **und** HCN. Aktuell (vor (a)/(b)):
H₂O ratio 2.65 / err 4.2e-5; HCN ratio 1.044 / err ~5e-5. H₂ Response = exakt 0
(gegated, muss so bleiben).

Nach (a)+(b): H₂O/HCN auf err<5e-5 bringen, dann CH₄/NH₃/C₆H₆ ergänzen und auf
gegated umstellen. Regression sicherstellen:
```bash
ctest -R "xtb_gradient|d4_dedq|energy_methods|xtb_cpscf"   # muss 6/6 bleiben
```
eeq-Pfad + Energien sind unberührt (Response läuft nur im mulliken-Zweig von
`calcDispersionEnergy`, `xtb_native.cpp` ~Z.611).

---

## Kritische Dateien

| Datei | Rolle |
|-------|-------|
| `src/core/energy_calculators/qm_methods/xtb_response.cpp` | `computeMullikenChargeResponse` (hier (a)+(b) ergänzen), `applyOrbitalHessian`, `solveZVector` |
| `src/core/energy_calculators/qm_methods/xtb_gradient.cpp` Z.523–649 | Section 5 — Vorlage für die Linearisierung |
| `src/core/energy_calculators/qm_methods/xtb_native.cpp` ~Z.611 | mulliken-Zweig (Verkabelung, fertig) |
| `test_cases/test_xtb_cpscf.cpp` | Gate-Test (Test C: Isolation + H₂-Symmetrie) |

Verfügbar in `computeMullikenChargeResponse`: `D_z`, `W_z`, `dq_sh/dq_at/ddp/dqp`,
`m_wfn.{q_at,dp_at,qp_at,q_sh}`, `m_mp_amat_{sd,dd,sq}`, `m_mp_mrad`, `m_pot.v_*`,
`xyz` (Bohr), `gfn2_params::*` (mp_dmp3/dmp5, mp_kexp, p_vcn/p_rad/mp_shift/mp_rmax).

---

## Risiken / Hinweise

- **Heikelster Teil** (Plan-Wortlaut): Faktoren (Spin ×2, Symmetrie), Vorzeichen,
  mrad/CN-Kette. Jeden Teilschritt FD-isoliert prüfen (erst SD, dann DD, dann SQ).
- **Nebenbefund (separates Issue):** vorbestehender GFN2-Gradient-Baseline-Fehler
  ~4e-5 Eh/Å (h-unabhängig, in eeq UND mulliken gleich; bei `test_xtb_gradient`-Tol
  5e-4 unsichtbar). Die Isolations-Methode umgeht ihn, aber er sollte separat
  untersucht werden (evtl. fehlender Term im Energie-Gradienten selbst).
- Einheiten: alles Eh/Bohr (Overlap-Grads in Bohr); `m_disp_gradient` wird in
  `xtb_gradient.cpp` Z.669 zu `m_gradient` addiert, dann `/= au` (→ Eh/Å).

---

## Referenzen

- Handy & Schaefer, J. Chem. Phys. 81 (1984) 5031 — Z-Vektor.
- tblite `coulomb/multipole.f90::get_multipole_gradient_0d`.
- [docs/D4_Q_RESPONSE.md](D4_Q_RESPONSE.md), [AIChangelog.md](../AIChangelog.md) (Phase 3b-Einträge).
