# AP 4 — Konkreter Implementierungsplan: Analytische Gradienten für `curcuma::xtb::XTB`

**Status:** Bereit zur Umsetzung (nach AP 3 oder parallel zu AP 3)
**Erstellt:** 2026-04-25
**Vorgängerdokument:** [`NATIVE_XTB_ROADMAP.md`](NATIVE_XTB_ROADMAP.md), Arbeitspaket 4
**Vorbedingung:** AP 1 abgeschlossen; AP 2 für End-to-end-Tests durch Wrapper, technisch jedoch unabhängig
**Implementation:** AI-generiert, nicht-getestet
**Komplexität:** Höchstes Risiko aller APs; klar in Sub-APs zerlegt

---

## Ziel

`curcuma::xtb::XTB` liefert nach AP 4 analytische Gradienten ∂E/∂R für GFN1-xTB und GFN2-xTB. Die Gradienten setzen sich additiv aus fünf Beiträgen zusammen, die in **fünf Sub-APs** unabhängig implementiert und getestet werden:

| Sub-AP | Beitrag | Datei | GFN1 | GFN2 | Komplexität |
|---|---|---|---|---|---|
| 4a | Repulsion | `xtb_h0.cpp` | ✅ | ✅ | niedrig |
| 4b | Isotrope Coulomb (ES2) | `xtb_coulomb.cpp` | ✅ | ✅ | niedrig |
| 4c | Third-order (ES3) | `xtb_thirdorder.cpp` | ✅ | ✅ | niedrig |
| 4d | H0 (Hellmann-Feynman + CN-Shift + shpoly) | `xtb_h0.cpp` | ✅ | ✅ | hoch |
| 4e | Multipol (AES2) | `xtb_multipole.cpp` | ❌ | ✅ | sehr hoch |

Die Sub-APs sind **bottom-up** geordnet: 4a/4b/4c können parallel und in beliebiger Reihenfolge erfolgen. 4d ist der größte Brocken (Hellmann-Feynman braucht die SCF-konvergente Dichte und Energy-Weighted Density `W`). 4e ist GFN2-spezifisch und kann zuletzt erfolgen.

**End-to-end-Garantie nach AP 4:**
- `./curcuma -opt H2O.xyz -method gfn2` konvergiert
- ∇E_analyt vs. ∇E_numerisch differiert maximal 1e-4 Eh/Bohr (mit FD-Schritt h = 5e-4 Bohr)

---

## Vorbedingungen

```bash
cd release && make -j4 curcuma 2>&1 | tail -3
./curcuma -sp ../test_cases/sqm_reference/molecules/H2O.xyz -method ngfn2 2>&1 | grep -i energy
```

Erwartet: build grün; ngfn2 liefert endliche Energie. Numerischer Wert wird nicht geprüft (das ist AP 5).

---

## Architektur-Entscheidungen

### A. Energy-Weighted Density `W`

Hellmann-Feynman benötigt zusätzlich zur Dichtematrix `P` die **energy-weighted density**

```
W_μν = 2 · Σ_i^occ ε_i · C_μi · C_νi
```

**Entscheidung:** `W` als Member `m_W` in `Wavefunction` cachen, in `xtb_scf.cpp::solveEigen()` berechnen — analog zu `P`. Vorteil: in jedem SCF-Schritt aktualisiert; bei Gradient-Aufruf direkt verfügbar.

```cpp
// xtb_native.h, struct Wavefunction
Matrix W;   // energy-weighted density: W_μν = 2 Σ_i^occ ε_i C_μi C_νi

// xtb_scf.cpp::solveEigen(), nach P-Berechnung
m_wfn.W = Matrix::Zero(nao, nao);
for (int i = 0; i < nocc; ++i) {
    const double eps = m_wfn.eps(i);
    for (int mu = 0; mu < nao; ++mu)
        for (int nu = 0; nu < nao; ++nu)
            m_wfn.W(mu, nu) += 2.0 * eps * m_wfn.C(mu, i) * m_wfn.C(nu, i);
}
```

### B. Gradient-Akkumulator

Eine zentrale Methode `XTB::calculateGradient()` baut den Gradienten additiv aus den Sub-AP-Methoden auf.

```cpp
// xtb_native.h
private:
    Matrix calculateGradient();   // assembled from sub-AP kernels

    // Sub-AP API (jeder schreibt in den übergebenen Akkumulator):
    void addRepulsionGradient(Matrix& grad) const;            // 4a, xtb_h0.cpp
    void addCoulombShellGradient(Matrix& grad) const;         // 4b, xtb_coulomb.cpp
    void addThirdOrderGradient(Matrix& grad) const;           // 4c, xtb_thirdorder.cpp
    void addH0Gradient(Matrix& grad) const;                   // 4d, xtb_h0.cpp
    void addMultipoleGradient(Matrix& grad) const;            // 4e, xtb_multipole.cpp (GFN2 only)
```

**`calculateGradient()` selbst:**

```cpp
// xtb_native.cpp
Matrix XTB::calculateGradient()
{
    if (!m_scf_converged) {
        CurcumaLogger::warn("XTB::calculateGradient called before SCF convergence");
    }
    Matrix grad = Matrix::Zero(m_atomcount, 3);
    addRepulsionGradient(grad);
    addCoulombShellGradient(grad);
    addThirdOrderGradient(grad);
    addH0Gradient(grad);
    if (m_method == MethodType::GFN2) {
        addMultipoleGradient(grad);
    }
    return grad;
}
```

### C. Schalter für `Calculation(gradient=true)`

Aktuell ignoriert `Calculation()` den `gradient`-Parameter. Nach AP 4:

```cpp
double XTB::Calculation(bool gradient)
{
    // ... bestehende SCF-Logik ...
    if (gradient && m_scf_converged) {
        m_gradient = calculateGradient();   // m_gradient von QMInterface
    }
    return m_E_total;
}
```

`m_gradient` ist Erbe aus `QMInterface` (Geometry-Type) — bleibt zwischen Aufrufen erhalten, wird hier überschrieben.

### D. Numerische Gradienten als Validierungs-Tool

Eine Helper-Methode für Sub-AP-Tests:

```cpp
// xtb_native.h, public:
Matrix calculateNumericalGradient(double h = 5e-4 /* in Angstrom */);
```

Implementation: zwei-seitige finite Differenz, geometriegestört, `Calculation(false)` aufrufen, alle anderen Caches korrekt invalidieren. Nutzt das `UpdateMolecule()`-Override aus AP 1.

---

## Sub-AP 4a — Repulsions-Gradient

**Datei:** `src/core/energy_calculators/qm_methods/xtb_h0.cpp`
**Referenz:** `tblite/src/tblite/classical/repulsion.f90`, alte `gfn2.cpp:855-873`
**Komplexität:** niedrig (paarweise, analytisch geschlossene Form)

### Mathematik

```
E_rep = Σ_{A<B} (Z_A·Z_B / R) · exp[ -(α_A·α_B)^kexp · R^rexp ]
                                    └────────────── E_pair ──────┘
∂E_pair/∂R = E_pair · [ -1/R - (α_A·α_B)^kexp · rexp · R^(rexp-1) ]
∂E_pair/∂R_A = (∂E_pair/∂R) · (R_A - R_B) / R
```

Beachte das Vorzeichen: ∂R/∂R_A = (R_A - R_B)/R = +n̂; ∂R/∂R_B = -n̂.

### Code-Skelett

```cpp
void XTB::addRepulsionGradient(Matrix& grad) const
{
    const int nat = m_atomcount;
    std::vector<double> xyz_bohr(3 * nat);
    for (int i = 0; i < nat; ++i) {
        xyz_bohr[3*i+0] = m_geometry(i,0) * AA_TO_AU;
        xyz_bohr[3*i+1] = m_geometry(i,1) * AA_TO_AU;
        xyz_bohr[3*i+2] = m_geometry(i,2) * AA_TO_AU;
    }
    const double kexp = (m_method == MethodType::GFN1) ? gfn1_params::rep_kexp : gfn2_params::rep_kexp;
    const double rexp = (m_method == MethodType::GFN1) ? gfn1_params::rep_rexp : gfn2_params::rep_rexp;

    auto get_alf  = [this](int z) {
        return (m_method == MethodType::GFN1)
            ? gfn1_params::rep_alpha[z-1] : gfn2_params::rep_alpha[z-1];
    };
    auto get_zeff = [this](int z) {
        return (m_method == MethodType::GFN1)
            ? gfn1_params::rep_zeff[z-1] : gfn2_params::rep_zeff[z-1];
    };

    for (int i = 0; i < nat; ++i) {
        const int zi = m_atoms[i];
        const double alf_i = get_alf(zi), zeff_i = get_zeff(zi);
        for (int j = 0; j < i; ++j) {
            const int zj = m_atoms[j];
            const double dx = xyz_bohr[3*i] - xyz_bohr[3*j];
            const double dy = xyz_bohr[3*i+1] - xyz_bohr[3*j+1];
            const double dz = xyz_bohr[3*i+2] - xyz_bohr[3*j+2];
            const double R  = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (R < 1e-12) continue;

            const double a_prod = alf_i * get_alf(zj);
            const double zz     = zeff_i * get_zeff(zj);
            const double a_pow  = std::pow(a_prod, kexp);
            const double R_pow_rexp     = std::pow(R, rexp);
            const double R_pow_rexp_m1  = std::pow(R, rexp - 1.0);

            const double E_pair = zz / R * std::exp(-a_pow * R_pow_rexp);
            const double dE_dR  = E_pair * (-1.0/R - a_pow * rexp * R_pow_rexp_m1);

            // Gradient akkumulieren — Einheit: Eh/Bohr (xyz_bohr ist Bohr).
            const double gx = dE_dR * dx / R;
            const double gy = dE_dR * dy / R;
            const double gz = dE_dR * dz / R;
            grad(i,0) += gx; grad(i,1) += gy; grad(i,2) += gz;
            grad(j,0) -= gx; grad(j,1) -= gy; grad(j,2) -= gz;
        }
    }
}
```

⚠️ **Einheitenkonvention:** Im gesamten AP 4 ist der Akkumulator **`Eh/Bohr`**. Die Konvertierung nach Eh/Å findet gegebenenfalls erst im Wrapper statt (Curcuma-Geometrie ist Å, intern Bohr für QM).

### Test 4a (eigenständig)

`test_xtb_grad_repulsion.cpp` analog zu `test_xtb_overlap.cpp`:

```cpp
auto mol = readXYZ("H2O.xyz");
curcuma::xtb::XTB x(curcuma::xtb::MethodType::GFN2);
x.QMInterface::InitialiseMolecule(mol);

// Energie + analytischer Repulsionsgradient (über calcRepulsionEnergy + addRepulsionGradient)
Matrix grad_an = Matrix::Zero(natoms, 3);
x.addRepulsionGradient(grad_an);

// Numerischer Gradient: nur Repulsionsbeitrag isoliert
Matrix grad_num = numericalRepulsionGradient(x, h=5e-4);

ASSERT_LT((grad_an - grad_num).cwiseAbs().maxCoeff(), 1e-5);
```

`addRepulsionGradient` muss dafür **public** sein oder ein Friend-Test. Empfehlung: für Test-Builds eine `friend class XTBGradientTest;` deklarieren. Alternativ: Test über `Calculation(true)` mit „nur Repulsion aktiv" — komplizierter, daher Friend-Pattern bevorzugt.

---

## Sub-AP 4b — Isotrope Coulomb-Gradient

**Datei:** `src/core/energy_calculators/qm_methods/xtb_coulomb.cpp`
**Referenz:** `tblite/src/tblite/coulomb/charge.f90:get_gradient` (Klopman-Ohno-Form)
**Komplexität:** niedrig (Klopman-Ohno-Ableitung in geschlossener Form)

### Mathematik

```
E_iso = 0.5 · Σ_{s,s'} q_sh(s) · γ(s,s') · q_sh(s')

γ_KO(R; ηs, ηs') = 1 / [ R^k + (η_avg)^(-k) ]^(1/k)        (k=2 für GFN2; gewichteter η)
∂γ/∂R          = -R^(k-1) · γ^(k+1)
∂E_iso/∂R_A = Σ_{s∈A, s'∉A} q_sh(s) · q_sh(s') · ∂γ/∂R · n̂(s,s')
```

### Code-Skelett

`xtb_coulomb.hpp` muss eine Funktion `gamma_derivative_dR` exportieren — entweder in der Header-Library oder neu in `xtb_coulomb.cpp` als file-statisch. Vorzugsweise letzteres (kein Header-Pollution):

```cpp
// xtb_coulomb.cpp
namespace {
    // Klopman-Ohno gamma + derivative w.r.t. R (Bohr).
    // Returns {gamma, dgamma/dR}.
    std::pair<double, double> gamma_and_dR(coulomb::Method m,
                                           int z_a, int z_b,
                                           int ang_a, int ang_b,
                                           double R)
    {
        // delegiert an xtb_coulomb.hpp Helpers, mit zusätzlich der Ableitung
        // ... (Implementation aus tblite/coulomb/charge.f90 portieren)
    }
}

void XTB::addCoulombShellGradient(Matrix& grad) const
{
    const auto meth = (m_method == MethodType::GFN1) ? coulomb::Method::GFN1
                                                     : coulomb::Method::GFN2;
    const int nsh = m_basis.nsh;
    const int nat = m_atomcount;

    std::vector<double> xyz_bohr(3 * nat);
    for (int i = 0; i < nat; ++i) {
        xyz_bohr[3*i+0] = m_geometry(i,0) * AA_TO_AU;
        xyz_bohr[3*i+1] = m_geometry(i,1) * AA_TO_AU;
        xyz_bohr[3*i+2] = m_geometry(i,2) * AA_TO_AU;
    }

    for (int s = 0; s < nsh; ++s) {
        const int A = m_basis.sh2at[s];
        const int z_a = m_basis.z[A], la = m_basis.ang_sh[s];
        for (int sp = s + 1; sp < nsh; ++sp) {
            const int B = m_basis.sh2at[sp];
            if (A == B) continue;        // intra-atomare Anteile haben dR = 0
            const int z_b = m_basis.z[B], lb = m_basis.ang_sh[sp];

            const double dx = xyz_bohr[3*A] - xyz_bohr[3*B];
            const double dy = xyz_bohr[3*A+1] - xyz_bohr[3*B+1];
            const double dz = xyz_bohr[3*A+2] - xyz_bohr[3*B+2];
            const double R  = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (R < 1e-12) continue;

            auto [g, dg_dR] = gamma_and_dR(meth, z_a, z_b, la, lb, R);
            (void)g;  // already used in energy

            const double q_a = m_wfn.q_sh(s), q_b = m_wfn.q_sh(sp);
            const double dE_dR = q_a * q_b * dg_dR;   // s≠s' und symmetrisch ⇒ Faktor 1, weil sp > s

            const double gx = dE_dR * dx / R;
            const double gy = dE_dR * dy / R;
            const double gz = dE_dR * dz / R;
            grad(A,0) += gx; grad(A,1) += gy; grad(A,2) += gz;
            grad(B,0) -= gx; grad(B,1) -= gy; grad(B,2) -= gz;
        }
    }
}
```

### Test 4b

```cpp
Matrix grad_an = Matrix::Zero(natoms, 3);
x.addCoulombShellGradient(grad_an);

Matrix grad_num = numericalGradientOf([&]() { return x.energyCoulombShell(); }, mol, h=5e-4);
ASSERT_LT((grad_an - grad_num).cwiseAbs().maxCoeff(), 1e-5);
```

⚠️ Numerischer Pfad muss bei jeder Verschiebung das gesamte SCF-Setup (`buildGammaMatrix`, neuer `q_sh`) neu laufen lassen — `q_sh` ist geometrieabhängig. Saubere Variante: die FD-Comparison läuft auf der **vollen** Energie nach AP 4d/4e, nicht auf `E_iso` isoliert. Daher **isolierte Tests sind nur als Smoke gedacht**, die echte Validierung erfolgt im AP-4-Schlusstest.

---

## Sub-AP 4c — Third-order-Gradient

**Datei:** `src/core/energy_calculators/qm_methods/xtb_thirdorder.cpp`
**Referenz:** `tblite/src/tblite/coulomb/thirdorder.f90:get_gradient`
**Komplexität:** niedrig (on-site, R-unabhängig)

### Mathematik

Third-order-Energie ist **on-site** — sie hängt nur von Ladungen ab, nicht von Distanzen:

```
GFN1: E3 = (1/3) Σ_A q_A^3 · Γ_A(z_A)
GFN2: E3 = (1/3) Σ_s q_sh(s)^3 · Γ_s(z_at(s), l_sh(s))
```

⇒ `∂E3/∂R = 0` direkt. Der einzige Beitrag zum Gradienten kommt aus der **q-Abhängigkeit** der Ladungen, die wiederum geometrieabhängig sind. Das wird über das **Hellmann-Feynman-Theorem** und die SCF-Konvergenz automatisch in der H0-Ableitung (4d) erfasst — explizite Behandlung in 4c **entfällt**.

### Konsequenz

`addThirdOrderGradient(grad)` ist eine **No-Op**:

```cpp
void XTB::addThirdOrderGradient(Matrix& grad) const
{
    // Third-order is on-site — no explicit geometric gradient.
    // The geometry-induced q-changes are captured via Hellmann-Feynman
    // in addH0Gradient() (Sub-AP 4d) through the SCF-converged P matrix.
    (void)grad;
}
```

Die Methode existiert dennoch, damit `calculateGradient()` symmetrisch bleibt und ein zukünftiger Beitrag (z. B. distanzabhängiger Third-order in einer Erweiterung) hier eingefügt werden kann.

---

## Sub-AP 4d — H0-Gradient (Hellmann-Feynman + CN-Shift + shpoly)

**Datei:** `src/core/energy_calculators/qm_methods/xtb_h0.cpp`
**Referenz:** `tblite/src/tblite/xtb/h0.f90:get_hamiltonian_gradient`, alte `gfn2.cpp:875-930`
**Komplexität:** **hoch** — drei verkettete Beiträge

### Mathematik

Die elektronische Energie ist `E_el = Tr(P · H0) + ES2 + ES3 + AES2`. Der H0-Beitrag zum Gradienten via Hellmann-Feynman:

```
∂(Tr(P·H0))/∂R_A = Tr( P · ∂H0/∂R_A ) - Tr( W · ∂S/∂R_A )
```

`W` ist die energy-weighted density (Architektur-Entscheidung A). Der Term `H0_μν = (ε_μ + ε_ν)/2 · π_μν · S_μν` zerfällt in:

1. **Overlap-Ableitung `∂S_μν/∂R`** — analytisch über CGTO-Derivate (`STO::calculateOverlapDerivative` existiert bereits in der alten Implementation; muss ggf. für die neue Basis-Struktur portiert werden).
2. **shpoly-Distanz-Polynom `∂π_ij/∂R`** — analytische Form aus `getHamiltonianH0()`:
   ```
   π_ij = (1 + sh_i · rr) · (1 + sh_j · rr),   rr = sqrt(R / (rad_i + rad_j))
   ∂π_ij/∂R = sh_i · (∂rr/∂R) · (1 + sh_j · rr) + (1 + sh_i · rr) · sh_j · (∂rr/∂R)
   ∂rr/∂R = 1 / (2 · sqrt(R · (rad_i + rad_j)))
   ```
3. **CN-Self-Energy-Shift `∂ε_sh/∂R`** — `ε_sh = ε⁰_sh - kcn · CN(at(sh))`, also
   ```
   ∂ε_sh/∂R_X = -kcn · ∂CN(at(sh))/∂R_X
   ```
   `∂CN/∂R` ist eine Sigmoid-Ableitung (siehe alte `gfn2.cpp:914-924` für die exakte Form bei `cn_gfn`).

### Strukturierung

`addH0Gradient` zerlegt sich in drei innere Helfer:

```cpp
void XTB::addH0Gradient(Matrix& grad) const
{
    addH0OverlapShpolyGradient(grad);   // (1) + (2): pro Shell-Paar
    addH0CnShiftGradient(grad);          // (3): pro Atom-Paar
}
```

**Helfer (1)+(2) — Schale-Schale-Schleife:**

```cpp
// Iteriere über Shell-Paare A < B, pro AO-Paar (mu, nu):
for (int sh_a < sh_b) {
    for (mu in shell sh_a) for (nu in shell sh_b) {
        double dS = STO_CGTO::overlap_derivative(...);
        double pi_ij_factor, dpi_dR;
        // ... aus getHamiltonianH0 portierte Form
        double avg_eps = 0.5*(se(sh_a) + se(sh_b));
        double k_factor = hs_factor;   // gleich wie in getHamiltonianH0
        double dH0_dR = k_factor * (dpi_dR * S(mu,nu) + pi_ij * dS) * avg_eps;
        // Hellmann-Feynman:
        double dE_pair = 2.0 * P(mu,nu) * dH0_dR  -  2.0 * W(mu,nu) * (k_factor * pi_ij * dS);
        grad(A,k) += dE_pair * n_AB[k];
        grad(B,k) -= dE_pair * n_AB[k];
    }
}
```

**Helfer (3) — CN-Shift-Schleife:**

```cpp
// ∂E/∂R_X = Σ_sh tr_sh(P_sh) · (-kcn[sh]) · ∂CN(at(sh))/∂R_X
// Σ_X' für jedes Paar (at(sh), X')
for (int sh = 0; sh < m_basis.nsh; ++sh) {
    int A = m_basis.sh2at[sh];
    double sum_P_sh = ... ;  // sum_{mu in sh} P(mu, mu)
    for (int X = 0; X < m_atomcount; ++X) {
        if (X == A) continue;
        double dCN = ... ;   // analytische Form aus tblite/data/coordnum.f90
        double dE_dR = -m_h0.kcn[sh] * sum_P_sh * dCN;
        // Symmetrische Verteilung A ↔ X
    }
}
```

⚠️ **Wichtig:** `cn_exp` (GFN1) und `cn_gfn` (GFN2) haben **unterschiedliche** Ableitungen. Beide existieren in `parameters/xtb_params_extra.hpp` als Energie-Funktionen — die Ableitungen müssen ergänzt werden:

```cpp
// xtb_params_extra.hpp (oder neue Datei xtb_cn_gradient.hpp):
inline std::vector<std::array<double,3>> cn_exp_gradient(...) { ... }
inline std::vector<std::array<double,3>> cn_gfn_gradient(...) { ... }
```

Form: Rückgabe ist die Jacobi-Matrix `∂CN_A/∂R_B[k]`. Implementation 1:1 aus tblite (`tblite/src/tblite/data/coordnum.f90:get_cn_gradient`).

### Test 4d

End-to-end FD-Test ist hier zwingend, weil 4d das größte Risiko trägt:

```bash
./curcuma -sp H2O.xyz -method ngfn2 -gradient-test 1 2>&1 | grep "max_abs_error"
```

Wo `-gradient-test 1` ein neuer CLI-Flag ist (oder direkter Test im sqm_reference-Setup). Erwartet: `max_abs_error < 1e-4 Eh/Bohr` für H₂O.

---

## Sub-AP 4e — Multipol-Gradient (GFN2 only)

**Datei:** `src/core/energy_calculators/qm_methods/xtb_multipole.cpp`
**Referenz:** `tblite/src/tblite/coulomb/multipole.f90:get_multipole_gradient`, alte `gfn2.cpp:932-972` (vereinfacht)
**Komplexität:** **sehr hoch** — vier Beiträge mit Damping-Funktionen

### Anteile

```
E_AES2 = E_dq + E_qq_diag + E_dd + E_cq
```

| Term | Form |
|---|---|
| **dipol-on-site** | `Σ_A f_A(CN_A) · |μ_A|^2` — ∂/∂R nur über `f_A(CN)` ⇒ CN-Gradient |
| **charge-dipole** | `Σ_{A≠B} q_A · μ_B · ∇γ_AB · f_dmp3(R)` |
| **dipole-dipole** | `0.5 · Σ_{A≠B} μ_A · T_AB · μ_B · f_dmp5(R)` |
| **charge-quadrupole** | `Σ_{A≠B} q_A · Θ_B : ∇∇γ_AB · f_dmp5(R)` |

`T_AB = (3·n̂⊗n̂ - I) / R^3` (Dipol-Tensor); `f_dmp_n = 1 / (1 + 6·(R/(rad_A + rad_B))^(2n))`.

### Strukturierung

```cpp
void XTB::addMultipoleGradient(Matrix& grad) const
{
    if (!m_mp_initialized) return;
    addMultipoleOnSiteCnGradient(grad);   // f_A(CN) · |μ_A|^2 — CN-only
    addMultipolePairGradient(grad);        // Charge-Dipole + Dipole-Dipole + Charge-Quadrupole
}
```

**Empfehlung:** Wegen Komplexität diesen Sub-AP mit ausführlichen Verifizierungsschritten implementieren — pro Anteil einen separaten FD-Test. Zwischendurch Energie-only validieren (TBLite-Dump-Vergleich) und erst dann Gradient anfassen.

### Realismus-Check

Die alte `gfn2.cpp:932-972`-Implementation enthält explizit den Kommentar „Simplified" (nur radiale Ableitung). Eine **vollständige** Multipol-Ableitung mit allen Tensor-Anteilen ist ein eigenständiges Projekt. Für AP 4 reicht zunächst die radiale Variante mit dokumentierter Genauigkeitsgrenze; **TBLite-genaue** Multipol-Gradienten sind eine separate Verbesserung (AP 6 oder später).

⇒ **Akzeptanzschwelle für 4e: max_abs_err < 1e-3 Eh/Bohr** (eine Größenordnung schwächer als bei 4a-d), explizit dokumentiert.

---

## Akzeptanzkriterien

### Sub-AP-Niveau

- [ ] **4a:** `addRepulsionGradient` vs. FD < 1e-5 Eh/Bohr (H₂O, CH₄, NH₃)
- [ ] **4b:** Voller `Calculation(grad)`-Pfad ohne 4d/4e — FD-Vergleich der „nur Coulomb"-Energieanteile < 1e-4 Eh/Bohr
- [ ] **4c:** No-Op verifiziert (siehe Mathematik-Erklärung)
- [ ] **4d:** End-to-end-Gradient (4a + 4b + 4d, GFN1) vs. FD < 1e-4 Eh/Bohr für H₂O, CH₄, NH₃, LiH, H₂
- [ ] **4e:** End-to-end-Gradient GFN2 (alle Beiträge) vs. FD < 1e-3 Eh/Bohr für gleiche Moleküle

### Globales Niveau

- [ ] `./curcuma -opt H2O.xyz -method gfn1` konvergiert (max ‖grad‖ < 1e-4 Eh/Bohr)
- [ ] `./curcuma -opt H2O.xyz -method gfn2` konvergiert
- [ ] CTests von AP 3 (DISABLED) bleiben DISABLED — Reaktivierung ist AP 5
- [ ] Build kompiliert ohne neue Warnings

### Nicht-Ziele

- ❌ Bit-genaue Übereinstimmung mit TBLite-Gradienten (AP 5)
- ❌ d-Funktionen (`ang_sh ≥ 2`) — `ao_to_type` returnt -1 für d in xtb_h0.cpp; daher d-Anteile im Multipol-/H0-Gradient ignoriert wie schon in der Energie. Dokumentieren als Limitation.
- ❌ DIIS, thermal smearing, open-shell

---

## Risiken und Stolperfallen

| Risiko | Eintrittswahrscheinlichkeit | Gegenmaßnahme |
|---|---|---|
| `STO::calculateOverlapDerivative` fehlt für die neue Basis-Struktur (`CGTO::Shell`) | hoch | Im Test 4d früh prüfen; ggf. neu in `STO_CGTO.hpp` implementieren — orthogonalisierte CGTOs erfordern auch im Gradient die Gram-Schmidt-Korrekturen |
| Vorzeichen- oder Faktor-2-Fehler im Hellmann-Feynman-Term | hoch | FD-Validierung pro Sub-AP; bei H₂O zuerst, danach Multi-Atom |
| `m_W` wird in `solveEigen` nicht kohärent zur SCF-Konvergenz aktualisiert (z. B. nach Damping) | mittel | `m_wfn.W` erst NACH dem Damping-Block in `Calculation()` final berechnen, nicht in `solveEigen` |
| CN-Gradienten in `parameters/xtb_params_extra.hpp` sind doppelte Implementation neben `tblite` | mittel | Direkt aus tblite-Fortran portieren; in einem Kommentar Quelldatei + Zeilenbereich vermerken |
| Multipol-Gradient liefert nur „simplified" Genauigkeit | mittel | Explizit als Limitation dokumentiert; AP 6 für vollständige Tensor-Ableitung |
| FD-Test bei kleinem h instabil (Fließkomma-Rauschen), bei großem h falsch | mittel | h = 5e-4 Bohr ≈ 2.6e-4 Å als Default; pro Test einmal h-Sweep verifizieren |
| Performance: Doppel-Schleife über AO-Paare im H0-Gradient | niedrig | Skaliert wie Energie-Aufbau (O(nao²)); vorerst kein Issue. AP 6 ggf. mit Symmetrie halbieren |

---

## Test-Plan

### Pro Sub-AP

Jeder Sub-AP bekommt einen **Standalone-Smoke-Test** im Stil der existierenden `test_xtb_*.cpp` unter `test_cases/sqm_reference/`:

- `test_xtb_grad_repulsion.cpp` (4a)
- `test_xtb_grad_coulomb.cpp` (4b)
- `test_xtb_grad_h0.cpp` (4d)
- `test_xtb_grad_multipole.cpp` (4e, GFN2 only)

Jeder Test:

1. Lädt ein Testmolekül aus `test_cases/sqm_reference/molecules/`
2. Initialisiert `XTB(MethodType::…)`
3. Ruft `Calculation(true)` auf
4. Vergleicht `Gradient()` mit `calculateNumericalGradient()`
5. Asserts `max |Δ| < threshold`

### Globaler End-to-End-Test (nach AP 4e)

Neuer CTest in `test_cases/cli/sqm/CMakeLists.txt`:

```cmake
add_cli_test(sqm 06_gfn2_gradient_h2o)
add_cli_test(sqm 07_gfn2_optimization_h2o)
add_cli_test(sqm 08_gfn1_gradient_h2o)
```

Skripte vergleichen `curcuma -gradient` mit FD-Gradient, `curcuma -opt` mit Konvergenz-Check.

### Performance-Smoke

```bash
time ./curcuma -opt benzene.xyz -method gfn2
```

Erwartet: < 30s für C₆H₆ auf 1 Thread (Vergleich mit alter monolithischer GFN2: war 0/7 → kein Referenzwert; gegen TBLite ist 5x langsamer akzeptabel).

---

## Fortschritt

| Datum | Sub-AP | Status | Notizen |
|-------|--------|--------|---------|
| 2026-04-25 | 4a | Plan | — |
| | 4b | | |
| | 4c | | |
| | 4d | | |
| | 4e | | |

## Schwierigkeiten / Blocker

- *Noch keine dokumentiert*

---

## Referenzen

- **Vorgänger:** [`NATIVE_XTB_AP3_PLAN.md`](NATIVE_XTB_AP3_PLAN.md)
- **Nachfolger:** [`NATIVE_XTB_AP5_PLAN.md`](NATIVE_XTB_AP5_PLAN.md)
- **TBLite-Quellen (Gradient-Pendants):**
  - `external/tblite/src/tblite/classical/repulsion.f90` — Repulsion-Gradient (4a)
  - `external/tblite/src/tblite/coulomb/charge.f90` — `get_gradient` (4b)
  - `external/tblite/src/tblite/xtb/h0.f90` — `get_hamiltonian_gradient` (4d)
  - `external/tblite/src/tblite/coulomb/multipole.f90` — `get_multipole_gradient` (4e)
  - `external/tblite/src/tblite/data/coordnum.f90` — `get_cn_gradient`
- **Alte Implementation als Plausibilitäts-Referenz:**
  - `src/core/energy_calculators/qm_methods/gfn2.cpp:825-981` (`calculateGradient`)
  - `gfn2.cpp:996+` (`calculateNumericalGradient`)
- **Neue Implementation (zu erweitern):**
  - `xtb_native.cpp` — Akkumulator-Methode
  - `xtb_h0.cpp` — 4a + 4d
  - `xtb_coulomb.cpp` — 4b
  - `xtb_thirdorder.cpp` — 4c (No-Op)
  - `xtb_multipole.cpp` — 4e

---

## Änderungshistorie

| Datum | Autor | Änderung |
|-------|-------|----------|
| 2026-04-25 | Claude | Erstdokument |
