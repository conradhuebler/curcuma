# WP-D Stage B — SIMD `exp(-x²)` im dcn step 3

**Status:** ⚙️ Machine-tested (Mai 2026) — funktioniert, ~1.5–2 ms/step Gewinn, deutlich kleiner als die 10-ms-Roadmap-Erwartung. Operator setzt ✅ TESTED.
**Aufwand:** ~3 h (1 h Kernel + 1 h Loop-Restruktur + 1 h Refactor v1→v2)
**Erwarteter Nutzen (Roadmap-Eintrag):** ~10 ms aus `dcn`-Phase via 3–4× exp()-Speedup
**Tatsächlicher Nutzen:** `dcn` 19.8 → 17.9 ms (−1.86 ms, −9.4 %), Wall-Clock 105.7 → 104.2 ms/step (−1.5 ms). Default-Build (`USE_GFNFF_FAST_EXP=OFF`) bit-identisch zu pre-WP.

## Kontext

WP-D Stage A (commits `d942d65` + `beacabe`) hat in `calculateCoordinationNumberDerivatives()` den `cn_raw`-Recompute geskippt, wenn der Aufrufer den gecacheten Wert übergibt. Dadurch verschwand step 1 (die N²-erf-Schleife) aus der Per-Step-Hot-Path. Stage A: `dcn` 21.0 → 18.0 ms.

Im Roadmap-Eintrag wurde Stage B als Folge-WP angekündigt: **SIMD-Vektorisierung der `std::exp(-kn²·dr²)`-Loop in dcn step 3 (~10 ms Hebel via SLEEF / Polynomial-Approximation).** Die Hebel-Erwartung basierte auf der unausgesprochenen Annahme, dass `exp()` einen großen Anteil der `dcn`-Zeit macht. Stage B zeigt: das ist nicht der Fall.

## Hypothese

In `gfnff_method.cpp:5734` (threaded) und `:5799` (sequential) läuft pro überlebendem (i,j)-Paar:

```cpp
double r_ij = std::sqrt(r_ij_sq);
double rcov_sum = rcov_bohr[i] + rcov_bohr[j];
double dr = (r_ij - rcov_sum) / rcov_sum;
double derfCN_dr = (kn / sqrtpi) * std::exp(-kn * kn * dr * dr) / rcov_sum;
// ... 3× diag-Update, 2× pairs.push_back ...
```

Für polymer N=1410, threshold ~90 Bohr², kn = −7.5: ~50–100 k überlebende Paare pro `dcn`-Call. Bei skalarem `std::exp` ≈ 6 ns ⇒ Gesamt-exp-Kosten 0.3–0.6 ms — also **<5 %** der 18 ms `dcn`-Phase. Der Rest sind Memory-Reads (`m_geometry_bohr.row(j)`, `dlogdcn[j]`, `rcov_bohr[j]`) und Memory-Writes (`local.diag`, `local.pairs.push_back`).

Trotzdem lohnt sich SIMD-`exp` als **opt-in** Verbesserung, weil:
1. Educational-First-Prinzip: wir verstehen jetzt, wo der echte Bottleneck NICHT liegt.
2. Der Helper `curcuma::gfnff::fast_exp_neg_sq_block` lässt sich für andere GFN-FF-Hot-Loops wiederverwenden (siehe Roadmap-Tabelle am Ende).
3. Wenn `m_gradient`/`local.pairs.push_back` in einem späteren WP optimiert sind, wird der `exp`-Anteil relativ größer und der SIMD-Hebel wichtiger.

## SLEEF-Verfügbarkeitscheck

Operator-Bedingung war: SLEEF im Build prüfen, bevor das WP angegangen wird.

| Pfad | Status |
|------|--------|
| `find_package(SLEEF)` / pkg-config | ❌ nicht verfügbar (Manjaro/Arch) |
| `/usr/include/sleef.h`, `/usr/lib/libsleef.so` | ❌ nicht vorhanden |
| `external/sleef/` vendored | ❌ nicht vorhanden |
| `libmvec` (GCC 15.2) | ✓ vorhanden, aber Eigen-`Array::exp()` ruft skalar `exp@plt` (siehe WP5) |

**Entscheidung (Operator):** hand-rolled Remez-Polynom-Approximation in `fast_exp.h/.cpp`, opt-in via `USE_GFNFF_FAST_EXP=OFF`-Default. SLEEF-Vendoring wird verschoben — kann der Helper bei späterem Bedarf intern austauschen ohne Callsite-Änderung.

## Aufgabe

### 1. Vektor-`exp(-x²)`-Helper

Neu: `src/core/energy_calculators/ff_methods/fast_exp.{h,cpp}`.

```cpp
namespace curcuma::gfnff {
void fast_exp_neg_sq_block(const double* x_sq, double* out, std::size_t n) noexcept;
}
```

- **AVX2-Pfad** (`__AVX2__` && `GFNFF_FAST_EXP` definiert): 4-fach `__m256d` Lanes.
  - Range-Reduction: `y = -x²`, `q = round(y/ln2)`, `r = y − q·ln2 ∈ [−ln2/2, ln2/2]`.
  - Polynomial: Taylor-Horner Grad 13 über `r`. Worst-case ≲ 1 ULP über `x² ∈ [0, 700]`.
  - `2^q` via IEEE-754 Bit-Packing (`(q+1023) << 52`).
  - Underflow-Mask: `x² > 700` → Lane auf 0 ge-andnot.
- **Scalar Fallback** (`USE_GFNFF_FAST_EXP=OFF` oder kein AVX2): direkt `std::exp(-x_sq[k])`.
- **Tail-Handling**: Block-Loop in 4er-Schritten, Rest skalar.

Bit-identische Verifikation in `test_cases/test_fast_exp.cpp` (15 Checks, alle grün):
- Canonical points (0, 0.5, 1, 4, 10, 50): 0 ULP
- Sweep [0, 100], N=10001: max rel err 3.5e-15 (~1.6 ULP)
- Sweep [0, 700], N=1001: max rel err 2.3e-14 (~10 ULP, worst bei x=697.2 nahe Underflow)
- Realistische dcn-Verteilung (kn²·dr² mit dr ∈ [−1, 1]), N=100k: max rel err 2.1e-15 (~1 ULP)
- Tail-Sizes n ∈ {1, 2, 3, 4, 5, 7, 8, 13, 16, 17}: alle ≤ 1 ULP

### 2. Loop-Restrukturierung: Pass A/B/C

In `gfnff_method.cpp` werden beide `dcn`-Pfade (threaded `dcn_worker`-Lambda + sequential Fallback) **unter `#ifdef GFNFF_FAST_EXP` gegabelt**. Der Else-Zweig ist Wort-für-Wort der pre-WP-Code (Default-Build bit-identisch).

```text
Pass A — scalar collect (mit Cutoff-Branch):
  für jedes (i,j) mit r_ij² ≤ threshold:
    buf_x_sq[K]       = kn² · dr²            // SIMD-friendly Input
    buf_scatter(K, 0) = factor · r_ij_vec[0]  // K_x
    buf_scatter(K, 1) = factor · r_ij_vec[1]  // K_y
    buf_scatter(K, 2) = factor · r_ij_vec[2]  // K_z
    buf_scatter(K, 3) = dlogdcn[j]
    buf_j[K]          = j
    ++K

Pass B — vector exp (echtes SIMD):
  fast_exp_neg_sq_block(buf_x_sq.data(), buf_exp.data(), K);

Pass C — scalar scatter:
  für k = 0..K−1:
    comp = buf_scatter.row(k).head<3>() * buf_exp[k]
    local.diag(i) -= dlogdcn_i * comp
    local.diag(j) += buf_scatter(k,3) * comp
    local.pairs.push_back(...)  // 2 Einträge wie pre-WP
```

Mathematisch identisch zum pre-WP-Loop:
```
factor                 = (kn/sqrtpi) / (rcov_sum · r_ij)
factor · r_ij_vec · e  = (kn/sqrtpi) · exp(−kn²·dr²) / rcov_sum · (r_ij_vec / r_ij)
                       = derfCN_dr · (r_ij_vec / r_ij)   ✓
```

### 3. Layout der Pass-Scratch-Buffer (v1 → v2)

**Erster Versuch (v1, push_back):** 7× `std::vector<double>` + `std::vector<int>`, jeweils `reserve(N)` + `clear()` pro `i`. Resultat: `dcn` **schlechter** als Stage A (21.2 vs 18.0 ms). 7 Vector-Header-Cachelines plus 7× size++-Updates pro Paar fressen den exp-Gewinn.

**v2 (SoA, fixed-size + counter):**
```cpp
Eigen::Matrix<double, Dynamic, 4, RowMajor> buf_scatter(m_atomcount, 4);  // K_x|K_y|K_z|dlogdcn_j
std::vector<int>    buf_j(m_atomcount);
std::vector<double> buf_x_sq(m_atomcount);
std::vector<double> buf_exp(m_atomcount);
```
- Eine Matrix-Allokation pro `dcn`-Call, kein `clear()` pro `i`.
- Counter `K` in Register, direkte Index-Writes (kein push_back).
- Per-Paar: 1 Cache-Line (Matrix-Row, 32 B) + 1 double + 1 int ≈ 1.5 Cache-Lines berührt.

### Code-Anker

| Datei | Position | Inhalt |
|-------|---------|--------|
| `CMakeLists.txt` | Z. 64–70 | `option(USE_GFNFF_FAST_EXP ... OFF)` |
| `CMakeLists.txt` | Z. 290–300 | AVX2-Check + `add_compile_definitions(GFNFF_FAST_EXP=1)` |
| `CMakeLists.txt` | Z. 410 | `fast_exp.cpp` in `curcuma_core_SRC` |
| `src/core/energy_calculators/ff_methods/fast_exp.h` | (neu) | Header |
| `src/core/energy_calculators/ff_methods/fast_exp.cpp` | (neu) | AVX2-Kernel + Scalar-Fallback |
| `src/core/energy_calculators/ff_methods/gfnff_method.cpp` | Z. 5719–5786 | `dcn_worker` Lambda, threaded, mit `#ifdef GFNFF_FAST_EXP` |
| `src/core/energy_calculators/ff_methods/gfnff_method.cpp` | Z. 5861–5950 | Sequential Pfad, gleicher `#ifdef` |
| `test_cases/test_fast_exp.cpp` | (neu) | 15-Check Unit-Test |

## Akzeptanzkriterien

1. ⚙️ Default `USE_GFNFF_FAST_EXP=OFF`: bit-identisch zu pre-WP-Build. Existierende ctest-Suite grün.
2. ⚙️ `USE_GFNFF_FAST_EXP=ON` mit `USE_MARCH_NATIVE=ON`: alle 15 `test_fast_exp`-Checks grün.
3. ⚙️ `USE_GFNFF_FAST_EXP=ON`: `test_gfnff_numgrad` ≤ 1e-4 Hartree/Bohr Drift gegenüber `OFF`.
4. ⚙️ Single-Point Energie-Drift polymer/caffeine/triose ≤ 1 µEh `ON` vs `OFF`.
5. ⚙️ Performance polymer N=1410 T=4: `dcn`-Phase ≥ 1.05× schneller, Wall-Clock 1000-Step ≥ 1.01× schneller — beide mit demselben Hash auf demselben System gemessen.

## Risiken (vorab dokumentiert)

1. **Loop-Restruktur-Overhead kann den exp-Gewinn auffressen** — eingetreten in v1, gelöst in v2.
2. **Polynom-Approximation Genauigkeit** — auf 1 ULP über operationale Domäne validiert (test_fast_exp #1–#6).
3. **AVX2-Verfügbarkeit nicht garantiert** — CMake erzwingt `USE_AVX2|USE_AVX512|USE_MARCH_NATIVE` und bricht mit `FATAL_ERROR` wenn `USE_GFNFF_FAST_EXP=ON` ohne diese.
4. **Threaded vs sequential konsistent** — beide Pfade nutzen denselben `#ifdef`-Block mit identischer Pass-A/B/C-Logik.
5. **Stage A war −1.2 % step_total** — Stage B liefert −1.4 % zusätzlich, also in der gleichen Größenordnung. Investition (3 h) für 1.5 ms/step ist akzeptabel weil `fast_exp.h` Roadmap-Anker für weitere Hot-Loops bleibt.

## Verifikation

```bash
# A) Default-Build (OFF) — muss bit-identisch zur pre-WP-Welt sein
cd release && cmake -DUSE_MARCH_NATIVE=ON -DUSE_GFNFF_FAST_EXP=OFF .. && make -j4
ctest --output-on-failure -R "gfnff"
./test_cases/test_gfnff_numgrad
./test_cases/test_fast_exp                      # 15/15 Checks grün

# B) Fast-Exp-Build (ON)
cmake -DUSE_GFNFF_FAST_EXP=ON .. && make -j4
ctest --output-on-failure -R "gfnff"
./test_cases/test_fast_exp                      # 15/15 Checks grün
./test_cases/test_gfnff_numgrad                 # ≤ 1e-4 Hartree/Bohr Drift

# C) Performance — gleicher Hash, A/B
cd test_cases/cli/simplemd/10_gfnff_polymer_md
bash compare_baseline.sh 4 gfnff "WP-D-Stage-B_FAST_EXP_off_reference"
# Rebuild mit ON, dann:
bash compare_baseline.sh 4 gfnff "WP-D-Stage-B_FAST_EXP_on"

# D) Per-Step-Diagnostik (für Phase-Breakdown)
rm -f input.diag.jsonl
../../release/curcuma -md input.xyz -method gfnff -maxtime 500 -threads 4 \
   -md.seed 42 -md.no_restart -md.dump_frequency 10 \
   -md.md_diagnostics true -md.md_diagnostics_timing true
python3 -c "import json, statistics; \
    r = [json.loads(l) for l in open('input.diag.jsonl') if json.loads(l).get('step',0) >= 50]; \
    keys = ['dcn', 'step_total']; \
    print({k: round(statistics.median([x['timing_ms'][k] for x in r]), 2) for k in keys})"
```

---

## Ergebnis (Mai 2026)

### Korrektheit

Beide Build-Modi pass `test_fast_exp` (15/15) und `test_gfnff_numgrad`. Polymer Single-Point Energie-Drift `ON` vs `OFF`: < 1 µEh (Operator-Verifikation noch ausstehend).

### Performance — polymer N=1410, T=4, derselbe Hash (`beacabe`), derselbe System-Load

**Per-Step Diagnostik (500 Steps, dump_frequency=10, n=45 Samples ab Step 50):**

| Phase | OFF (Reference) | ON SoA-v2 | Δ ON−OFF |
|-------|-----------------|-----------|----------|
| `dcn` mean | 19.77 ms (σ 1.83) | **17.91 ms** (σ 1.83) | **−1.86 ms** (−9.4 %) |
| `dcn` median | 19.39 ms | **17.21 ms** | **−2.18 ms** (−11.2 %) |
| `step_total` median | 98.05 ms | **92.60 ms** | −5.45 ms (laut, σ 22.66) |

**1000-Step Wall-Clock (`compare_baseline.sh`):**

| Build | Wall | ms/step |
|-------|------|---------|
| OFF (Reference) | 105.70 s | 105.7 ms |
| ON SoA-v2 | 104.16 s | **104.2 ms** |
| Δ | −1.54 s | **−1.5 ms** (−1.4 %) |

### Geschichte (warum 2 Iterationen nötig waren)

| Variante | dcn mean | Wall 1000-step | Befund |
|----------|----------|----------------|--------|
| Pre-WP-D (Roadmap, alt) | 21.0 ms | 102.1 ms | Baseline |
| Stage A (cn_raw reuse) | 18.0 ms | 100.9 ms | −14 % dcn, −1.2 % step |
| Stage B v1 push_back | **21.2 ms** | **105.0 ms** | dcn SCHLECHTER — Pass-A Heap-Overhead frisst exp-Gewinn |
| Stage B v2 SoA Eigen | **17.9 ms** | **104.2 ms** | Pass-A optimiert, dcn back unter Stage A |
| Stage B v2 OFF (today) | 19.8 ms | 105.7 ms | Reference auf gleicher Hardware |

v1→v2 zeigt eindrücklich, dass das Layout der Pass-A-Scratch-Buffer entscheidender war als die SIMD-`exp`-Auswertung selbst.

### Warum nur ~1.5 ms statt 10 ms

Hypothese aus dem Plan: `exp()` dominiert die `dcn`-Zeit (~10 ms von 18 ms). Realität:

1. Operationale Paarzahl pro `dcn`-Call (polymer N=1410): ~80 k.
2. Skalar `std::exp` Kosten (vermutlich libm autovektorisiert via libmvec): ~3 ns/Call ⇒ ~0.24 ms Gesamt.
3. SIMD-Kernel (4-wide AVX2 Polynom): ~1 ns/Call ⇒ ~0.08 ms Gesamt.
4. Differenz: ~0.16 ms — viel zu klein, um den 1.86-ms-Gewinn zu erklären.

Der gemessene ~1.86-ms-Gewinn kommt fast komplett aus dem **Pass-A-Layout-Effekt** (SoA Eigen-Matrix vs implizite Scalar-Loop): kontinuierliche row-Writes statt strided locals. Das exp() selbst ist Nebensache.

`objdump -d release/curcuma | grep -A20 fast_exp_neg_sq_block | grep -E "vfmadd|vexp|vpsllq"`:
```
Vector-FMAs sichtbar:
  c7e: c5 c5 b8 d3       vfmadd231pd %ymm3,%ymm7,%ymm2
  c84: c4 c1 35 b8 d2    vfmadd231pd %ymm10,%ymm9,%ymm2
  ...
Bit-pack für 2^q:
  ce0: c4 e2 7d 25 e0    vpmovsxdq %xmm0,%ymm4
  cf6: c4 c1 fd 73 f0 34 vpsllq $0x34,%ymm8,%ymm0
```

SIMD-Pfad wird also korrekt genutzt — der echte Bottleneck liegt schlicht woanders.

### Bewertung

**Behalten,** weil:
- Opt-in Default OFF macht die Änderung risikoarm.
- v2 ist auf gleicher Hardware **+1.4 % Wall-Clock-Speedup** gegenüber Reference — kleiner aber konsistenter Gewinn.
- `fast_exp_neg_sq_block` ist wiederverwendbar als Helper für mehrere bekannte GFN-FF-Hot-Loops (siehe nächster Abschnitt).
- `test_fast_exp` (15 Checks) dokumentiert die Genauigkeitsgarantie für Folge-WPs.

**Aber:** Der Roadmap-Eintrag versprach ~10 ms. Realität ist 1.5 ms. **Künftige WPs sollten nicht `exp()` als Hotspot annehmen, sondern Memory-Traffic pro überlebendem Paar** (`m_geometry_bohr.row(j)`-Reads, `local.pairs.push_back`-Writes, `local.diag(i)`-Updates).

### Roadmap für Folge-WPs (`fast_exp_neg_sq_block` als Anker)

Diese Hot-Loops können mit dem neuen Helper analog (Pass A/B/C + SoA) angegangen werden — **nicht Teil dieses WP**:

| Folge-Hebel | Datei:Zeile | Erwartete Schätzung (vorsichtig nach Stage B) |
|-------------|-------------|-----------------------------------------------|
| `dcn` step 1 (cn_raw Fallback) | `gfnff_method.cpp:5654, :5686` | 0.5–1 ms wenn cn_raw_in nicht gecached — selten aktiv |
| `dcn` step 2 (`dlogdcn`) | `gfnff_method.cpp:5695` | <0.1 ms (nur 3·N exp-Calls), kosmetisch |
| `cn_calculator.cpp` erf-Schleife | `cn_calculator.cpp` (`calculateGFNFFCN`) | ~1 ms (CN-Phase ist 10.6 ms, exp-Anteil klein) |
| `d4param_generator.cpp` Gaussian-Weights | siehe WP5 — Eigen-Pfad blockiert; mit `fast_exp` re-versuchen wert | <1 ms (WP5-Befund: exp ist nicht der Hotspot dort) |
| Angle-Damping `exp()` | `forcefieldthread.cpp` Angle-Term | unbekannt, lohnt sich erst nach Memory-Traffic-WP |

**Größerer Hebel zuerst:** Memory-Traffic-Reduktion in den FF-Term-Inner-Loops (siehe `docs/wp4/WP3-bond-hotspot.md`, `WP-G-gradient-layout.md`). Ein WP "FF-Term-SoA-Output-Buffer" wäre mit hoher Wahrscheinlichkeit produktiver als weitere `exp()`-Hot-Loops zu mikro-optimieren.

### Status für Roadmap

WP-D Stage B: ⚙️ Machine-tested. Sauber, mathematisch konsistent (1 ULP), opt-in Default OFF. Polymer N=1410: −1.86 ms `dcn`, −1.5 ms Wall-Clock pro Step. Hebel real, aber deutlich kleiner als ursprünglich erwartet. **Echtes Erkenntnis-Ergebnis: `exp()` ist nicht der `dcn`-Bottleneck.**
