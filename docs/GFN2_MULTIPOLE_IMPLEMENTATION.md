# GFN2-xTB Multipole-Implementierung (TBLite-Referenz)

## Überblick

Die GFN2-xTB-Multipolwechselwirkung erweitert das isotrope Coulomb-Modell um
dipolare (l=1) und quadrupolarische (l=2) Beiträge. Sie wird in TBLite durch
drei Module realisiert:

1. **`integral/multipole.f90`** — Berechnung der Multipolintegrale über CGTOs
2. **`coulomb/multipole.f90`** — Dämpfungsradien und Interaktionsmatrizen (Amat)
3. **`scf/potential.f90`** — Potentialaufbau + Fock-Beitrag

---

## Schritt 1: Multipolintegrale (dp_int, qp_int)

### 1.1 Rohintegrale mit globalem Ursprung

Für jedes AO-Paar (μ,ν) werden cartesische Dipol- und Quadrupolintegrale mit
Ursprung bei (0,0,0) berechnet:

```
D_k(μ,ν) = ⟨μ| r_k |ν⟩              für k = x,y,z
Q_kl(μ,ν) = ⟨μ| r_k r_l |ν⟩          für kl = xx,xy,yy,xz,yz,zz
```

**Interface**: `MI::cgto_multipole(shell_a, shell_b, xa, ya, za, xb, yb, zb, type_a, type_b, S, D, Q)`
- Liefert `S` (Overlap), `D[3]` (Dipol), `Q[6]` (Quadrupol, roh cartesisch)
- Nutzt die STO→CGTO-Expansion aus `STO_CGTO.hpp`
- Definiert in `xtb_multipole_ints.hpp` (fertig implementiert)

### 1.2 Ursprungsverschiebung zum Atom des zweiten Index

TBLite-Konvention: Der Multipoloperator ist auf das Atom des *letzten* Index
zentriert. Die Verschiebung von Ursprung 0 → Atomposition R erfolgt via:

```
dp_int[k](μ,ν) = D_k(μ,ν) - R_k · S(μ,ν)

qp_raw[0](μ,ν) = Q_xx(μ,ν) - 2·R_x·D_x(μ,ν) + R_x²·S(μ,ν)  [xx]
qp_raw[1](μ,ν) = Q_xy(μ,ν) - R_x·D_y(μ,ν) - R_y·D_x(μ,ν) + R_x·R_y·S(μ,ν)  [xy]
...analog für yy, xz, yz, zz
```

### 1.3 Traceless-Transformation (Quadrupol)

TBLite verwendet das *spurfreie* Cartesian-Quadrupol-Paket (6 Komponenten):

```
qp_int[k] = 1.5 · qp_raw[k] - 0.5 · tr(qp_raw) · δ_{diag, k}

mit tr = 0.5 · (qp_raw[xx] + qp_raw[yy] + qp_raw[zz])
und δ_{diag, k} = 1 für k=xx,yy,zz, 0 für k=xy,xz,yz

Faktor 0.5 in tr: TBLite verwendet den normalisierten Tracedefekt:
  tr = (Q_xx + Q_yy + Q_zz) / 3
Im Code: tr = 0.5 * (qxx + qyy + qzz), dann
  qp_int[xx] = 1.5 * qxx - tr = 1.5*qxx - 0.5*(qxx+qyy+qzz)
             = (3*qxx - qxx - qyy - qzz)/2 = (2*qxx - qyy - qzz)/2
             = qxx - (qxx+qyy+qzz)/3  ✓ (spurfrei)
```

Nach diesen Schritten liegen vor:
- `dp_int[3](nao, nao)` — Dipolintegrale mit Ursprung am Spaltenatom
- `qp_int[6](nao, nao)` — tracelose Quadrupolintegrale mit Ursprung am Spaltenatom

---

## Schritt 2: Koordinationszahlen (CN)

GFN2 verwendet die `cn_gfn()`-Funktion (doppelt-exponentielle Form) aus
`xtb_params_extra.hpp`:

```
cn_i = Σ_{j≠i} 1 / (1 + exp( -k·(R_ij - r_cov(i) - r_cov(j)) ))
```

Bereits implementiert in `XTB::computeCoordinationNumbers()` → `xtb_h0.cpp:52`.

---

## Schritt 3: Dämpfungsradien (mrad)

Pro Atom wird ein CN-abhängiger Dämpfungsradius berechnet:

```
mrad[i] = rad_i + (rmax - rad_i) / (1 + exp(-mp_kexp · (cn_i - vcn_i - mp_shift)))

mit Parametern aus gfn2_params:
  p_rad[Z], p_vcn[Z]       — elementabhängig
  mp_shift, mp_kexp         — global (mp_shift ≈ -2.0, mp_kexp ≈ 2.0)
  mp_rmax                   — globaler Maximalradius
```

Zusätzlich werden zwei Kernel-Parameter pro Atom geladen:
- `vec_dkernel_iat[i] = p_dkernel[Z]` — Dipol-Selbstwechselwirkungs-Kernel
- `vec_qkernel_iat[i] = p_qkernel[Z]` — Quadrupol-Selbstwechselwirkungs-Kernel

---

## Schritt 4: Interaktionsmatrizen

### 4.1 Charge-Dipole (amat_sd)

```
r_ij = |R_i - R_j|
g1 = 1/r_ij,  g3 = g1³,  g5 = g1⁵
rr = 0.5 · (mrad_i + mrad_j) · g1
fdmp3 = 1 / (1 + 6 · rr^mp_dmp3)    (mp_dmp3 ≈ 3)
fdmp5 = 1 / (1 + 6 · rr^mp_dmp5)    (mp_dmp5 ≈ 5)

amat_sd[k](jat, iat) = R_k · g3 · fdmp3    für k = x,y,z
```

Speicherkonvention: `amat_sd[k](row=j, col=i)` = Einfluss von Atom i auf
Atom j in Richtung k.

### 4.2 Dipole-Dipole (amat_dd)

```
dd_iso   = g3 · fdmp5
dd_anis  = 3 · g5 · fdmp5

amat_dd[a][b](jat, iat) = δ_{ab} · dd_iso - v_a · v_b · dd_anis
```

Wobei `v_a = (R_i - R_j)_a`, also die a-te Komponente des Abstandsvektors.

### 4.3 Charge-Quadrupole (amat_sq)

```
amat_sq[0](jat, iat) =  vx² · g5 · fdmp5    (xx)
amat_sq[1](jat, iat) = 2·vx·vy · g5 · fdmp5  (xy)
amat_sq[2](jat, iat) =  vy² · g5 · fdmp5    (yy)
amat_sq[3](jat, iat) = 2·vx·vz · g5 · fdmp5  (xz)
amat_sq[4](jat, iat) = 2·vy·vz · g5 · fdmp5  (yz)
amat_sq[5](jat, iat) =  vz² · g5 · fdmp5    (zz)
```

Die Off-Diagonal-Komponenten (xy, xz, yz) haben Faktor 2 (tblite-Konvention
für das gepackte Quadrupol).

---

## Schritt 5: Potentialaufbau (pro SCF-Iteration)

### 5.1 Dipol-Potential vdp(3, nat)

```
vdp(k, iat) = Σ_jat  amat_sd[k](iat, jat) · q_at(jat)
            + Σ_a,jat  amat_dd[k][a](iat, jat) · dpat(a, jat)
            + 2 · dkernel[iat] · dpat(k, iat)     [on-site XC]
```

### 5.2 Quadrupol-Potential vqp(6, nat)

```
vqp(k, iat) = Σ_jat  amat_sq[k](iat, jat) · q_at(jat)
            + 2 · qkernel[iat] · qpat(k, iat) · mpscale_q[k]
```

Mit `mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0}` (Kompensation der
Off-Diagonal-Faktoren in amat_sq).

### 5.3 Atomares Zusatzpotential vat_extra(nat)

```
vat_extra(iat) = Σ_k,jat  amat_sd[k](jat, iat) · dpat(k, jat)
               + Σ_k,jat  amat_sq[k](jat, iat) · qpat(k, jat)
```

Dies ist die transponierte SD/SQ-Wechselwirkung: Der Dipol/Quadrupol auf
Atom jat erzeugt ein skalares Potential auf Atom iat.

### 5.4 AO-Potential

```
v_ao(μ) = v_sh[shell(μ)] + vat_extra[atom(μ)]
```

Wobei `v_sh` bereits das isotrope Coulomb + Third-Order enthält.

---

## Schritt 6: Fock-Beitrag

Der Multipolbeitrag zum Fock-Operator folgt dem tblite-Muster
`add_vmp_to_h1`:

```
F(μ,ν) += -0.5 · [ Σ_k dp_int[k](μ,ν) · vdp(k, atom(ν))
                  + Σ_k dp_int[k](ν,μ) · vdp(k, atom(μ))
                  + Σ_k qp_int[k](μ,ν) · vqp(k, atom(ν))
                  + Σ_k qp_int[k](ν,μ) · vqp(k, atom(μ)) ]
```

Wichtig: Jeder Term ist symmetrisch bzgl. (μ,ν), da `dp_int[k](μ,ν) ≡ dp_int[k](ν,μ)`
und `qp_int[k](μ,ν) ≡ qp_int[k](ν,μ)` (reelle Basis).

---

## Schritt 7: Multipolenergie

```
E_mp = Σ_iat Σ_jat [
    Σ_k  dpat(k, iat) · amat_sd[k](iat, jat) · q_at(jat)
  + Σ_a,b  0.5 · dpat(a, iat) · amat_dd[a][b](iat, jat) · dpat(b, jat)
  + Σ_k  qpat(k, iat) · amat_sq[k](iat, jat) · q_at(jat)
]
+ Σ_iat [
    Σ_k  dkernel[iat] · dpat(k, iat)²
  + Σ_k  qkernel[iat] · qpat(k, iat)² · mpscale_q[k]
]
```

---

## Schritt 8: Mulliken-Multipole (pro SCF-Iteration)

Nach der Dichtematrix P aus dem Eigenproblem werden die atomaren Multipole
via Mulliken-Analyse aktualisiert:

```
dpat(k, iat) = - Σ_{μ∈iat, ν} P(ν,μ) · dp_int[k](ν,μ)
qpat(k, iat) = - Σ_{μ∈iat, ν} P(ν,μ) · qp_int[k](ν,μ)
```

Das Minuszeichen: In TBLite-Konvention sind die Multipole als
*Valenz*-Abweichungen von der Referenz definiert (wie q_sh).

---

## Datenabhängigkeiten (pro SCF-Iteration)

```
v_sh(q_sh)          ← γ · q_sh                [isotrop, bereits impl.]
v_sh_third(q_sh)    ← third-order              [bereits impl.]
vdp(q_at, dpat)     ← amat_sd · q_at + amat_dd · dpat + 2·dkernel·dpat
vqp(q_at, qpat)     ← amat_sq · q_at + 2·qkernel·qpat·mpscale
vat_extra(dpat,qpat)← amat_sd^T · dpat + amat_sq^T · qpat

F(H0, S, v_sh, vdp, vqp, dp_int, qp_int)
→ diag(F, S) → P
→ Mulliken: q_sh, dpat, qpat
→ dampen, repeat
```

## Konstanten (aus gfn2_params.hpp)

| Parameter | Wert | Bedeutung |
|-----------|------|-----------|
| `gexp` | 2.0 | Exponent im γ-Kernel |
| `mp_shift` | -2.0 | CN-Shift für Dämpfung |
| `mp_kexp` | 2.0 | CN-Exponent für Dämpfung |
| `mp_dmp3` | 3 | Dämpfungsexponent für SD |
| `mp_dmp5` | 5 | Dämpfungsexponent für DD/SQ |
| `mp_rmax` | 20.0 | Maximaler Dämpfungsradius |

## Referenzcode

Die vollständige Referenzimplementierung befindet sich in
`test_cases/sqm_reference/test_xtb_scf_snapshot.cpp`, Zeilen 249-810.
Die Integral-Bibliothek in `src/core/energy_calculators/qm_methods/xtb_multipole_ints.hpp`
ist fertig und getestet.
