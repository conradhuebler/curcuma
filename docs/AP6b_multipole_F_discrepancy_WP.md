# AP6b Working Plan — Native xTB F-Matrix Discrepancy gegenüber tblite

**Status:** GELÖST + **FIX IMPLEMENTIERT** (2026-05-27). Root Cause = fehlendes selbstkonsistentes D4-Potential + native D4 war eine GFN-FF-Näherung. Fix: exakter dftd4-per-Referenz-D4-Port + SCF-Kopplung. **H₂ vs tblite: E −3.3e-9, HOMO +1.9e-6** (war 1.5e-4); H₂O/CH₄/NH₃ sub-0.1 mEh. GFN-FF unverändert (7e-8). Regression 6/6 grün.

### FIX-IMPLEMENTIERUNG (2026-05-27)
- **`dispersion/d4_charge_scaling.h`** (neu): geteilte exakte dftd4-Zeta (`d4_zeta/d4_dzeta/d4_d2zeta`, Port von model.f90:725). GFN-FFs `zetaChargeScale` = formel-identischer gc=1-Spezialfall (unangetastet).
- **`D4ParameterGenerator::weightedC6Gfn2`** (neu): exakte per-Referenz-C6 + dC6/dq + dC6/dCN (+ d²C6/dq²). wf=6 + ngw-Multi-Gauss (aus `m_refcn` via set_refgw), gi=eta·gc (gc=2), qref=`m_refq[ref]`+zeff. **Kein Tabellen-Port nötig: `m_refq` IST bereits der dftd4-GFN2-Satz, `zeta_zeff`=zeff, `zeta_c`=eta.**
- **`D4Params.per_reference_charge`**: GFN2=true (exakt), GFN-FF=false (CN-only-single-zeta + CUDA-Layout intakt).
- **`D4Evaluator::computeEnergyAndGradient`** GFN2-Pfad nutzt weightedC6Gfn2; `pairEnergyAndGradient` (GFN-FF) unangetastet.
- **`XTB::addDispersionPotential`** (re-added): dE_D4/dq → `pot.v_at` in der SCF-Schleife. GFN2-D4 jetzt Mulliken-selbstkonsistent (Energie+Potential+Gradient).
- **CPSCF-Kernel (Phase E) nicht nötig**: per-Referenz-dEdq macht die mulliken-CPSCF-Antwort konsistent → `test_xtb_cpscf` grün ohne expliziten ∂²E_D4/∂q²-Term.
- **Rest** H₂O/CH₄/NH₃ ~3-9e-5: NICHT D4 (H₂ beweist D4 exakt) — separater kleinerer Effekt (gated Multipol-Response/CN), außerhalb D4-Scope. Optionale Politur: D4-CN an dftd4 angleichen.

---

## FIX-VERSUCH (2026-05-27): SCF-Kopplung implementiert, dann sauber zurückgerollt

Der strukturelle Fix (`addDispersionPotential`: ∂E_D4/∂q an Mulliken-Ladungen → `pot.v_at` in jeder SCF-Iteration, GFN2-Default „mulliken") wurde implementiert und **funktioniert grundsätzlich** (Energie mEh→μEh, H₂-HOMO-Fehler 1.5e-4→7.4e-5). Bei der Validierung zeigten sich aber zwei Probleme, die zeigen, dass die volle tblite-Genauigkeit eine **komplette dftd4-Modell-Portierung** erfordert — zu groß/risikoreich für einen einzelnen Fix. Daher zurückgerollt; **Baum ist grün** (6/6 Regression). Diagnose-Tooling bleibt.

### Warum nur ~50% des Fehlers geschlossen wurde — die native D4 ist eine GFN-FF-Näherung
- Unsere native GFN2-D4 (`d4param_generator.cpp` + `d4_evaluator.cpp`) nutzt **ein einziges per-Atom-zeta × CN-gewichtetes C6** (GFN-FF-Stil), tblite/dftd4 nutzt **per-Referenzzustand ladungsabhängige Gewichte**.
- **Dominanter Faktor 2:** unsere Zeta-Steilheit nutzt `c = eta` (= `zeta_c`), dftd4 nutzt `c = eta·gc` mit **gc=2.0** (`model.f90:454`, `gi = eta·gc`). Das halbiert unser dEdq. (`zeta_c[H]=0.47259` = dftd4 `chemical_hardness[H]` exakt; der Kommentar „×2" in `gfnff_par.h:893` ist falsch.)
- **Rest-~3%:** unsere CN-gewichtete C6 selbst ist ~3% von tblite entfernt (unsere D4-Energie ist ~3% daneben) — breitere Modell-Näherung (CN, Gauss-Gewichte).
- **Tabellen, die exakt sind:** `zeta_zeff` = dftd4 `effective_nuclear_charge` ✓; `zeta_c` = dftd4 `chemical_hardness` ✓.
- **Tabellen, die NICHT passen:** `m_refq` (Hirshfeld, `d4_refq=2`) ≠ tblites GFN2-Referenzladungen (`set_refq_gfn2`, `d4.f90:80` nutzt `ref=d4_ref%gfn2`).

### Zweites Problem: CPSCF-Gradienten-Antwort
Die SCF-Kopplung verändert den konvergierten Zustand, aber der Z-Vektor-Response-Kernel (`computeMullikenChargeResponse`) enthält den D4-Selbstkonsistenz-Term (∂²E_D4/∂q²) nicht → `test_xtb_cpscf` Test C regrediert (HCN 7.75e-6→1.17e-4; die „mulliken==eeq"-Prämisse ist durch die Selbstkonsistenz ohnehin ungültig). Voller Gradient (`test_xtb_gradient`, tol 5e-4) bleibt grün.

### Vollständige Spezifikation für den exakten Fix (Folge-AP)
1. **`addDispersionPotential`** (wieder einführen): ∂E_D4/∂q_A an `m_wfn.q_at` → `pot.v_at` in der SCF-Schleife (`xtb_native.cpp` nach `addMultipolePotential`). dEdq nur mit `with_gradient=true` aus `D4Evaluator::computeEnergyAndGradient` verfügbar.
2. **GFN2-D4-Zeta exakt**: per-Referenz-Gewichtung `W_i^ref = gw_i^ref(CN)·zeta(ga=3, gi=zeta_c[Z]·2.0, refq_i^ref+zeff, q_i+zeff)`, `C6 = ΣΣ W_i·W_j·c6ref`; zeff=`zeta_zeff[Z]`, refq = **GFN2-Referenzladungen** (dftd4 `set_refq_gfn2`, in `reference.inc`; unser `d4_refq` ist Hirshfeld — GFN2-Satz portieren). Energie + dC6/dq + dC6/dCN konsistent. GFN-FF-Pfad (gc=1, single-zeta) **nicht** anfassen → GFN2-spezifische Funktion.
3. **CPSCF-Kernel** um ∂²E_D4/∂q²-Term erweitern; `test_xtb_cpscf` Test C neu konzipieren.
4. **Validierung**: H₂/H₂O/CH₄/NH₃ Energie+Orbitale vs tblite (Dumps in `release_tblite/dumps/`, `diff_multipole_potential.py`); Regression grün.

### Diagnose-Tooling (bleibt, im Baum)
tblite-Patch exponiert Integrale + Potentiale/Momente (`*_tblite`-JSON-Keys); `CURCUMA_TBLITE_ACC`-env im Dump-Tool; `diff_multipole_ints.py` / `diff_multipole_potential.py`.

---

## ERGEBNIS Option (b) — Potential-Audit (2026-05-27): ROOT CAUSE

Der ~1.5e-4 Eh F-Shift ist **kein Multipol-Bug**. Er entsteht, weil das native GFN2-SCF das **selbstkonsistente D4-Dispersions-Potential** nicht in den Fock einbezieht.

**Beweiskette (tblite-Patch erweitert um `pot%vat/vdp/vqp` + `wfn%dpat/qpat`, alle bit-genau exponiert):**
1. `dpat`/`qpat` (atomare Momente aus tblite-P): bit-identisch (≤1e-15).
2. CN, mrad, `amat_sd/dd/sq`: tblite-Laufzeitwerte **bit-identisch** mit unseren (per Diagnose-Dump verifiziert).
3. tblites Multipol-`get_potential` auf der konvergierten Dichte liefert vat = −0.0315754 — **identisch mit unserem `vat_extra` −0.0315765** (Δ~1e-6). Der Multipol-Potential-Pfad ist also korrekt.
4. Das nach dem SCF erfasste `pot%vat` = −0.0314267. Die Differenz **+1.487e-4** entsteht in `iterator.f90` **zwischen `coulomb%get_potential` und `dispersion%get_potential`**: GFN2 nutzt **selbstkonsistentes D4**, dessen Ladungs-Potential (∂E_D4/∂q) via `add_pot_to_h1` in die Fock-Matrix und damit in tblites gespeicherte Orbitalenergien eingeht.
5. `test_xtb_scf_snapshot` (H₂): Δε_max = 1.504e-4 ≈ der D4-Potential-Beitrag (+1.487e-4). Match.

**Unser nativer Code** (`xtb_native.cpp:262`, `calcDispersionEnergy()`) addiert D4 nur als **Energie (post-SCF)** + Gradient; der Kommentar (`xtb_native.cpp:622-626`) hält den SCF-Response explizit für „deferred". Genau dieser Term fehlt.

**Fix (offen):** In jeder SCF-Iteration den D4-Atompotential-Term ∂E_D4/∂q_A auf `pot.v_at` addieren (analog tblites `dispersion%get_potential`). Die ∂E_D4/∂q-Maschinerie existiert bereits (`m_disp_dEdq`, AP7-Gradientenarbeit) — sie muss in den Fock-Build (`xtb_scf.cpp::buildFock` / Potential-Assembly) statt nur in den Gradienten einfließen. Erwartung danach: H₂ Δε_max < 1e-6, GFN2-Energie-Restfehler weiter reduziert.

**Tooling (wiederverwendbar):** tblite-Patch exponiert jetzt Integrale (`dipole/quadrupole_integral`) UND Potentiale/Momente (`atomic_dipole/quadrupole`, `potential_vat/vdp/vqp`); `dump_tblite_multipole` schreibt `*_tblite`-Keys (+ `CURCUMA_TBLITE_ACC` env-Toggle). Audits: `diff_multipole_ints.py`, `diff_multipole_potential.py`.

---

## ERGEBNIS Option (a) — Multipol-Integral-Audit (2026-05-26): Hypothese A widerlegt
**Vorgängerarbeit:** AP6 (DIIS-Integration, commit a6922de). Nach AP6 erreicht GFN2 3/7 vs TBLite. Verbleibende Fehler ~1–23 mEh werden zum großen Teil von fehlender D4-Dispersion verursacht (AP6c geplant), aber ein systematischer Anteil kommt aus einer **Orbital-Energie-Diskrepanz** die selbst für H₂ (2 AO) auftritt.

---

## ERGEBNIS Option (a) — Multipol-Integral-Audit (2026-05-26)

**Hypothese A (dp_int/qp_int weichen ab) ist widerlegt.** Die nativen GFN2-Multipol-AO-Integrale stimmen **bit-identisch** mit tblite überein.

**Methode:** tblite (Tag `9ca469a`) minimal gepatcht, um die internen `ints%dipole`/`ints%quadrupole` (post-`get_hamiltonian`/`shift_operator`, traceless, atom-zentriert — exakt die Arrays, die in den Fock-Build eingehen) über die C-API zu exponieren. `dump_tblite_multipole` schreibt sie als JSON-Keys `dipole_integral_tblite`/`quadrupole_integral_tblite` neben unsere Rekonstruktion (`gf2.dp_int`/`gf2.qp_int`, Produktionspfad via `MI::cgto_multipole`). Diff: `test_cases/sqm_reference/diff_multipole_ints.py`.

**Befund (max|Δ|, getrennt Same-Atom- vs Off-Atom-Blöcke):**

| Molekül | nao | dp same/off | qp same/off |
|---|---|---|---|
| H₂ | 2 | 2.2e-16 / 1.1e-16 | 2.2e-16 / 0 |
| H₂O | 6 | 2.2e-16 / 2.4e-16 | 4.4e-16 / 5.0e-16 |
| CH₄ | 8 | 2.2e-16 / 1.9e-16 | 4.4e-16 / 4.2e-16 |
| NH₃ | 7 | 2.2e-16 / 3.1e-16 | 4.4e-16 / 1.1e-15 |

Alle ≤ Maschinenpräzision, inkl. der zuvor verdächtigten Same-Atom-s-p-Blöcke. Der "globaler Ursprung + uniformer Spalten-Atom-Shift"-Ansatz (`xtb_multipole.cpp:setupMultipole`) ist exakt äquivalent zu tblites "Produktzentrum + `shift_operator`".

**Wo der 1.5e-4-Shift tatsächlich herkommt (neu eingegrenzt):** Für H₂ sind die Diagonal-Integrale `dp_int(0,0)=qp_int(0,0)=0` (traceless-Quadrupol einer s-Funktion ist null). Der Multipol koppelt also **nicht** über die Diagonal-Integrale, sondern über **`vat_extra`** (= `amat_sd^T·dpat + amat_sq^T·qpat`, für H₂ = −0.0316 ≠ 0) zurück auf das **isotrope Atompotential** `vat`, das diagonal in F eingeht. Der Fehler liegt damit **stromabwärts der Integrale**: in `dpat`/`qpat`, `amat_sd/dd/sq` oder der Kontraktion zu `vat_extra`/`vdp`/`vqp`.

**Nächster Schritt → Option (b):** denselben tblite-Patch erweitern, um `pot%vat`/`vdp`/`vqp` und `wfn%dpat`/`qpat` zu exponieren (Abschnitt 5, B1–B7). Direkter Vergleich entscheidet: stimmen `dpat`/`qpat` → Bug in `amat`/Kontraktion; weichen sie ab → Bug in Mulliken-Moment-Bildung/Referenz-Momenten.

---

## 1. Symptom

Im Test `sqm_scf_H2_gfn2` (sowie LiH/H₂O/CH₄/NH₃/C₆H₆ für GFN2):

```
ε[0] mine=-0.47079047  ref=-0.47094090  Δ=1.504e-04
ε[1] mine=+0.18460054  ref=+0.18445365  Δ=1.469e-04
```

`mine` = unsere F-Matrix diagonalisiert mit tblite-S; `ref` = tblite-stored `orbital_energies`. Beide Eigenwerte verschoben um **+1.5e-4** (gleichmäßig nach oben).

Bei einer 2×2-symmetrischen `F` mit `S = [[1,s],[s,1]]`:
- `ε_- = (a+b)/(1+s)`,  `ε_+ = (a-b)/(1-s)`
- Beide Eigenwerte um den gleichen Betrag verschoben ⇒ **uniformer Diagonal-Shift in F**: `ΔF(0,0) = ΔF(1,1) ≈ -1.5e-4`.

Da für H₂ `q_sh = 0` (Symmetrie) ist `v_sh = γ·q_sh = 0` und der Diagonal-Shift kann nur aus dem GFN2-Multipol-Pfad kommen → `vat_extra` oder dem `qp_int·vqp`-Beitrag.

Back-Rechnung bestätigt:
```
F_tblite(0,0) (aus ε rekonstruiert) = -0.36054949
F_ours(0,0)   (aus Test-Output)     = -0.36039965
ΔF(0,0)                              = -1.498e-04
```

## 2. Was bereits verifiziert ist (alles bit-identisch mit tblite)

| Komponente | Verifikation | Referenz |
|---|---|---|
| `mp_dmp3=3.0, mp_dmp5=4.0, mp_shift=1.2, mp_kexp=4.0, mp_rmax=5.0` | Direkter Vergleich mit `gfn2.f90:471` | `gfn2_params.hpp:947-955` |
| `p_rad`, `p_vcn`, `p_dkernel`, `p_qkernel` | Element-für-Element Diff = 0 | `gfn2_params.hpp` vs `gfn2.f90:600+` |
| CN-Formel `gfn_count` (ka=10, kb=20, r_shift=2) | Identische Formel | `xtb_params_extra.hpp:cn_gfn` vs `ncoord/gfn.f90:126` |
| Kovalenzradien `covalent_rad_2009` × 4/3 | Element-für-Element identisch | `xtb_params_extra.hpp` vs `data/covrad.f90` |
| `amat_sd/dd/sq` Konstruktion (vec, g3, g5, fdmp3, fdmp5) | Formel-identisch | `xtb_multipole.cpp:180-219` vs `coulomb/multipole.f90:429-476` |
| `vat_extra` Formel (`Σ amat_sd^T·dpat + amat_sq^T·qpat`) | Identisch | `xtb_multipole.cpp:282-294` vs `coulomb/multipole.f90:259, 264` |
| `mrad` Formel (`rad + (rmax-rad)/(1+exp(-kexp·(cn-vcn-shift)))`) | Identisch | `xtb_multipole.cpp:163-172` vs `coulomb/multipole.f90:330-361` |
| STO-3G Tabellen `pAlpha3`, `pCoeff3` | Bit-identisch (verified via Python) | `STO_CGTO.hpp:54-63` vs `basis/slater.f90:112-147` |
| CGTO-Normalisierung `(2α/π)^0.75·(4α)^(l/2)/sqrt((2l-1)!!)` | Identisch (dfactorial-Indexierung verifiziert) | `STO_CGTO.hpp:170-175` vs `basis/slater.f90:503-506` |
| Overlap `S(0,1)` für H₂ | `0.6631299012502960` vs `0.6631299012502958` (Δ=2.2e-16) | Python-Nachrechnung |
| Self-Overlap `S(0,0)` für H₂ | `1.0000000000657978` (beide gleich, kleine STO-3G-Fit-Restungenauigkeit) | Python-Nachrechnung |
| Dipol-Integral `dp_int[x](1,0)` für H₂ | `0.46428611` (Test-Output) entspricht Python-Nachrechnung exakt | siehe Reproduktion unten |

## 3. Restliche Hypothese (Quelle des 1.5e-4 Shifts)

Drei Restkandidaten, geordnet nach Wahrscheinlichkeit:

### Hypothese A — `qp_int`/`dp_int` Bra/Ket-Asymmetrie bzw. Same-Atom-Block (am wahrscheinlichsten)

tblite konstruiert `dpint(jao, iao)` und `dpint(iao, jao)` **mit separatem Origin-Shift** über `shift_operator` (`xtb/h0.f90:521-564`):
- `dpint(:, jao, iao)` ← `dtmpi(:)` (Operator centered on Ket = atom iao)
- `dpint(:, iao, jao)` ← `dtmpj(:)` aus `shift_operator(vec, s, di, qi, dj, qj)` mit `dj(k) = di(k) + vec(k)*s` (Operator centered on Bra = atom jao)

Zusätzlich erfolgt im `shift_operator` eine traceless-Re-Transformation des Quadrupol-Anteils (`xtb/h0.f90:558-563`):
```
qj(1) = qi(1) + 1.5·qj_raw(1) - tr_shifted
qj(2) = qi(2) + 1.5·qj_raw(2)
...
```

**Unsere Implementierung** (`xtb_multipole.cpp:setupMultipole`, Zeilen 117-152) berechnet stattdessen `dp_global` einmal für alle (μ,ν) und shiftet **uniform** mit der Spalten-Atom-Position:
```cpp
m_dp_int[k](mu, nu) = dp_global[k](mu, nu) - R_{atom of nu} · S(mu, nu);
```
Anschließend wird die Spur-freie Transformation auf das geshiftete Quadrupol angewandt.

**Verdacht:** Für die Same-Atom-Blöcke (iat == jat) und/oder für ν<μ-Einträge ergibt unser Ansatz numerisch unterschiedliche Werte als tblite's separate-Bra/Ket-Konstruktion mit nachträglicher traceless-Re-Transformation. Speziell der `tr`-Term in `shift_operator:554` benutzt die **geshiftete** Spur, nicht die ursprüngliche — das ist nicht trivial gleichwertig.

### Hypothese B — Integral-Cutoff
tblite verwendet eine Adjacency-List mit Cutoff (`xtb/h0.f90:215`). Wir berechnen alle (μ,ν)-Paare ohne Cutoff. Für H₂ (2 Atome) sollte das **nicht** relevant sein.

### Hypothese C — Subtile tblite-Konvention die wir übersehen
z.B. Diagonal-Behandlung in `get_hamiltonian:294-...` (die zweite OMP-Schleife — wir haben diesen Pfad noch nicht detailliert auditiert).

## 4. Working Plan — Option (a): Multipol-Integral-Audit

### Ziel
Verifizieren ob `dp_int[k](μ,ν)` und `qp_int[k](μ,ν)` zwischen unserer Implementierung und tblite bit-identisch übereinstimmen — element-für-element, inkl. Same-Atom-Blöcke und Asymmetrien.

### Schritte

**A1: Diagnostic-Dump in unserer setupMultipole hinzufügen**
- Datei: `src/core/energy_calculators/qm_methods/xtb_multipole.cpp`
- Nach Zeile 152 (Ende von Schritt 2): bei verbosity ≥ 3 eine Datei `dp_int_ours.txt` schreiben mit allen Werten `m_dp_int[k](mu, nu)` und `m_qp_int[k](mu, nu)` für alle (k, mu, nu).

**A2: Reference-Dump aus tblite extrahieren**
Tblite stores `dpint` und `qpint` **nicht** über die public C-API. Workaround:
- Build mit `release_tblite/` (existiert bereits laut git status)
- `dump_tblite_multipole` läuft schon — modifizieren so dass es **direkt** tblite's `ints%dipole` und `ints%quadrupole` extrahiert (statt sie über unsere `cgto_multipole` zu rekonstruieren)
- Aktuell extrahiert `dump_tblite_multipole.cpp:280-290` nur S, H, P. Nötig:
  - In tblite: `tblite_set_calculator_save_integrals(ctx, calc, 1)` ist schon gesetzt → `ints%dipole` und `ints%quadrupole` werden bewahrt
  - Aber **nicht** in der API exposed. Patches in tblite nötig (siehe Option (b)) ODER:
  - Alternative: Datei `external/tblite/src/tblite/xtb/singlepoint.f90:351-354` zeigt dass `ints%dipole`/`ints%quadrupole` derzeit NICHT in `results` kopiert werden. Patch hinzufügen.

**A3: Element-für-Element-Vergleich**
- Python-Script in `scripts/diff_multipole_ints.py`: liest `dp_int_ours.txt` und `dp_int_tblite.txt`, druckt max-Diff per (k, atom-pair).
- Erwartetes Resultat: Wenn Hypothese A stimmt, sind same-atom Blöcke (insbesondere quadrupol) signifikant unterschiedlich.

**A4: Falls Bug bestätigt → Fix in `setupMultipole`**
- Statt uniformen Shift: separate Konstruktion für Bra-/Ket-Centered Integrale entsprechend `shift_operator` aus `xtb/h0.f90:521-564`.
- Oder: traceless-Transformation NACH separater Bra/Ket-Behandlung anwenden, nicht vorher.

### Kritische Dateien
| Datei | Zeilen | Zweck |
|---|---|---|
| `src/core/energy_calculators/qm_methods/xtb_multipole.cpp` | 67-152 | unser setupMultipole, ggf. fix |
| `src/core/energy_calculators/qm_methods/xtb_multipole_ints.hpp` | 94-170 | `primitive_multipole`, `cgto_multipole` |
| `external/tblite/src/tblite/xtb/h0.f90` | 200-335 | tblite `get_hamiltonian` (referenz) |
| `external/tblite/src/tblite/xtb/h0.f90` | 521-564 | `shift_operator` (referenz) |
| `external/tblite/src/tblite/integral/multipole.f90` | 354-425 | `multipole_cgto` mit traceless-Transform (referenz) |
| `test_cases/sqm_reference/test_xtb_scf_snapshot.cpp` | 270-322 | Existierender Test der dp_int rekonstruiert |

### Reproduktion (current state)
```bash
cd /home/conrad/src/claude_curcuma/sqm/curcuma/release
make -j4
ctest -R "sqm_scf_H2_gfn2" --output-on-failure -V 2>&1 | grep -E "DIAG|ε\["
# Expected: Δε_max=1.504e-04 Eh
```

Python-Quick-Check (S, dp_int Nachrechnung):
```python
import json, numpy as np
from math import pi
with open('release/test_cases/sqm_reference/dumps/H2_gfn2.json') as f:
    d = json.load(f)
alpha_raw = [2.227660584e+0, 4.057711562e-1, 1.098175104e-1]
coeff_raw = [1.543289673e-1, 5.353281423e-1, 4.446345422e-1]
z2 = 1.23**2
alpha = [a*z2 for a in alpha_raw]
coeff = [c*(2*a/pi)**0.75 for c, a in zip(coeff_raw, alpha)]
R = np.array(d['atoms'][1]['xyz_bohr']) - np.array(d['atoms'][0]['xyz_bohr'])
S01 = sum(ci*cj*(pi/(ai+aj))**1.5*np.exp(-ai*aj/(ai+aj)*R@R)
          for ci,ai in zip(coeff,alpha) for cj,aj in zip(coeff,alpha))
print(f"S(0,1) ours = {S01}  tblite = {d['overlap'][0][1]}  diff = {S01 - d['overlap'][0][1]:.2e}")
```

## 5. Working Plan — Option (b): tblite C-API erweitern um vat/vdp/vqp/dpint/qpint zu exportieren

### Ziel
Direkten element-für-element-Vergleich zwischen tblite-internen Potentialen (`pot%vat`, `pot%vdp`, `pot%vqp`) und Integralen (`ints%dipole`, `ints%quadrupole`) mit unseren Werten ermöglichen.

### Schritte

**B1: Tblite Fortran-Side erweitern**
- Datei: `external/tblite/src/tblite/api/result.f90`
- Neue Routinen:
  ```fortran
  subroutine get_result_dipole_integral_api(verror, vres, dpint) bind(C, ...)
  subroutine get_result_quadrupole_integral_api(verror, vres, qpint) bind(C, ...)
  subroutine get_result_multipole_potential_api(verror, vres, vat, vdp, vqp) bind(C, ...)
  ```
- Daten kommen aus `res%results%dipole`, `res%results%quadrupole`, plus neuer Felder `res%results%pot_vat`, `pot_vdp`, `pot_vqp`.

**B2: tblite-singlepoint patchen um pot-Werte zu speichern**
- Datei: `external/tblite/src/tblite/xtb/singlepoint.f90`
- Nach SCF-Konvergenz (~Zeile 320-340): `pot%vat`, `pot%vdp`, `pot%vqp` in `results` kopieren (analog zu Z.352-353 für overlap/hamiltonian).

**B3: tblite-results-type erweitern**
- Datei: `external/tblite/src/tblite/results.f90`
- Felder hinzufügen: `dipole(:,:,:)`, `quadrupole(:,:,:)`, `pot_vat(:)`, `pot_vdp(:,:)`, `pot_vqp(:,:)`.

**B4: C-Header**
- Datei: `external/tblite/include/tblite/result.h`
- Deklarationen für die drei neuen Getter-Funktionen.

**B5: dump_tblite_multipole erweitern**
- Datei: `test_cases/sqm_reference/dump_tblite_multipole.cpp`
- Aufrufe der neuen API-Funktionen, JSON-Output erweitern um `dipole_integral`, `quadrupole_integral`, `pot_vat`, `pot_vdp`, `pot_vqp`.

**B6: Existierende JSON-Dumps neu erzeugen**
```bash
cd release
./test_cases/sqm_reference/dump_tblite_multipole gfn2 ../test_cases/water.xyz \
   test_cases/sqm_reference/dumps/H2O_gfn2.json
# usw. für alle 7 Testmoleküle
```

**B7: Diagnostic-Test schreiben**
- Datei: `test_cases/sqm_reference/test_vat_extra_comparison.cpp`
- Lädt JSON, vergleicht unsere computed vat vs tblite's stored vat element-für-element.
- Liefert direkte Antwort: stimmt vat? Wenn ja, ist das Problem dp_int. Wenn nein, ist das Problem amat_*.

### Tradeoff (a) vs (b)
- **(a)**: Schneller, weniger invasiv. Aber indirekt — leitet Bug aus dp_int-Differenzen ab, nicht direkt aus vat-Differenz.
- **(b)**: Aufwendiger (tblite-Patches, evtl. Upstream-PR), aber gibt direkte Antwort und ist wiederverwendbar für weitere Diagnosen (z.B. cn-Werte, selfenergy, etc.).

**Empfehlung:** **Erst (a)**. Wenn Hypothese A bestätigt (qp_int same-atom Diff signifikant) reicht der Fix in `xtb_multipole.cpp`. Erst falls (a) keine Diff zeigt und Mysterium bleibt, lohnt (b).

## 6. Erwartetes Ergebnis nach Fix

Falls Hypothese A korrekt und der Bug behoben:
- H₂ `Δε_max` < 1e-7 Eh (Maschinenpräzision)
- H₂O `Δε_max` < 1e-7 Eh
- CH₃OH/CH₃OCH₃ Energie-Fehler reduziert auf <0.5 mEh (bleibt: D4-Dispersion)
- C₆H₆/Koffein bleiben groß (dominiert von D4)

GFN2 sollte nach AP6b + AP6c (D4) 7/7 PASS erreichen.

## 7. Verworfene Hypothesen (nicht weiter verfolgen)

- **`updown_to_magnet`-Skalierung**: ist no-op für nspin=1 (verified an `wavefunction/spin.f90:104-114`).
- **`vat_extra`-Formel-Bug**: Formel ist identisch zu tblite `get_potential` (`coulomb/multipole.f90:243-268`).
- **Damping-Konstanten unterschiedlich**: alle 5 Konstanten identisch verifiziert.
- **STO-3G Koeffizienten**: bit-identisch verifiziert.
- **Overlap-Matrix unterschiedlich**: 2.2e-16 Diff für H₂(0,1) — Maschinenpräzision.

## 8. Kontextverknüpfungen

- AP6 commit: `a6922de` (DIIS in SCF eingebunden)
- Vorgängerstatus: `cc99a1d` (AP5b+fix: GFN2 multipole integral Pulay gradient + gradient unit fix)
- Memory-Eintrag: `project_diis_integration` (AP6-Resultate), `project_gfn2_multipole_debug` (älterer vat_extra-Verdacht — durch dieses WP teilweise korrigiert)
- Reference: tblite source unter `external/tblite/src/tblite/`
- Test infrastructure: `test_cases/sqm_reference/` mit `test_xtb_scf_snapshot`, `dump_tblite`, `dump_tblite_multipole`
