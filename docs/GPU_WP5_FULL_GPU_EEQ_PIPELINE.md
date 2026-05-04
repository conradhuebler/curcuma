# WP5: Vollständig GPU-residenter EEQ-Pfad

**Kategorie**: Großer Aufwand  
**Aufwand**: ~7–12 Tage  
**Wirkung**: Transformativ — eliminiert alle CPU-Sync-Punkte aus dem heißen Pfad  
**Abhängigkeiten**: WP2 (k_build_eeq_rhs), WP3 (Pair-List CN), WP4 (Phase-2-Graph)  
**Status**: 🤖 Geplant — nur nach WP2+WP3 validiert beginnen

---

## Ziel

Nach WP2+WP3 sieht der heiße Pfad noch so aus:

```
GPU (main stream): k_cn_compute_pairs → k_build_eeq_rhs
GPU (eeq stream):  k_eeq_build_matrix → potrf (lazy) → potrs
CPU (parallel):    setCNDerivatives() → recordD4CNValues()
--- Stream-Sync ---
GPU (main):        setEEQCharges() H2D → launchChargeDependentAndFinish()
```

Das Ziel von WP5 ist, auch die CPU-Arbeiten im mittleren Block (`setCNDerivatives`, `recordD4CNValues`) in GPU-Kernel umzuformen oder zu eliminieren, und dann alles in einem einzigen zusammenhängenden GPU-Graph zu fassen.

**Endzustand:**
```
cudaGraphLaunch(full_step_graph, stream)  // Gesamter Schritt: CN → EEQ → Coulomb → Gradient
cudaStreamSynchronize(stream)             // Einmal pro Schritt
// Kein CPU-Code im heißen Pfad
```

---

## Analyse: Was hält noch CPU-Involvement aufrecht?

### 1. `m_gfnff->prepareCNAndEEQ(gradient=true, skip_eeq=true)`

**Was es tut**: Berechnet CNF-Derivate (`m_cnf`-Vektor) und CN-Ableitungen für den Coulomb-Gradienten-Term (Term 1b). Diese werden via `setCNDerivatives(cn, cnf, {})` auf die GPU hochgeladen.

**GPU-Lösung**: `cnf[i]` ist topologie-konstant (bereits nach WP2 als `d_cnf` auf GPU). Der `dCN/dr`-Term für den Coulomb-Gradienten wird bereits vom `k_cn_chainrule`-Kernel berechnet. Der einzige verbleibende CPU-Schritt ist die Einreichung von `cnf` in `d_cnf_for_coulomb` — aber `d_cnf` ist nach WP2 bereits auf GPU. Ein einfacher D2D-Copy oder Referenz-Umleitung reicht.

**Aufwand**: Klein (1–2 Stunden nach WP2).

### 2. `m_gfnff->recordD4CNValues(cn_std)`

**Was es tut**: Speichert CN-Werte für den nächsten Schritt (Skip-Check für dc6dcn Neuberechnung). Benötigt CN als `std::vector<double>` auf der CPU.

**GPU-Lösung**: Statt CN herunterzuladen und zu vergleichen, kann ein GPU-Kernel `k_check_dc6dcn_skip` die Differenz zur referenz-CN direkt auf Device berechnen und ein `int`-Flag setzen. Das Flag (1 Bool) kann dann mit `cudaMemcpyAsync` asynchron heruntergeladen werden.

**Aufwand**: Mittel (3–4 Stunden). Neuer GPU-Kernel + Device-Buffer für Referenz-CN.

### 3. CPU Schur-Komplement (EEQ, nfrag=1)

**Was es tut** (gfnff_gpu_method.cpp:512–520):
```cpp
double S = 0.0, Cz1 = 0.0;
for (int i = 0; i < N; ++i) { S += m_eeq_Z2[i]; Cz1 += m_eeq_z1[i]; }
double lambda = (Cz1 - rhs_constraints[0]) / S;
for (int i = 0; i < N; ++i) charges[i] = m_eeq_z1[i] - m_eeq_Z2[i] * lambda;
```

Das ist zwei O(N) Reduktionen + ein O(N) Schreiben. Entspricht exakt `k_eeq_reduce_sums` + `k_eeq_schur_nfrag1`, die bereits existieren in `eeq_solver_gpu.cu` (der `solveAndComputeCharges`-Pfad, der als "4s slower" markiert ist).

**Warum war `solveAndComputeCharges` langsamer?** Das Problem war nicht die Kernel selbst, sondern dass der Pfad auch den D2H-Download von Ladungen + H2D-Upload via `setEEQCharges()` inkludierte. Nach WP2 liegen EEQ-Ladungen direkt als `d_charges`-Device-Buffer vor — der `solveAndComputeCharges`-Pfad ist auf D2D-Ebene extrem schnell. Das "4s slower" bezieht sich auf den alten Pfad mit unnötigem D2H+H2D.

**GPU-Lösung**: Bestehende Kernel `k_eeq_reduce_sums` + `k_eeq_schur_nfrag1` direkt nach `potrs` auf dem EEQ-Stream starten. Ladungen bleiben als `d_charges` auf GPU. `setEEQCharges()` entfällt.

**Aufwand**: Mittel (2–3 Stunden), nach WP2.

### 4. `finalizeCNForCPU()` (Sync-Punkt)

**Nach WP2**: Nur noch für CNF-Derivate im Gradient-Pfad nötig (Punkt 1 oben). Nach vollständiger GPU-Lösung von Punkt 1 kann `finalizeCNForCPU()` vollständig entfernt werden.

---

## Implementierungsplan (sequenziell aufbauend auf WP2+WP3)

### Phase WP5-A: GPU-seitiges Schur-Komplement (2–3h)

Bestehende Kernel `k_eeq_reduce_sums` und `k_eeq_schur_nfrag1` (bereits in `eeq_solver_gpu.cu`) direkt im EEQ-Solve-Pfad starten.

```cpp
// In EEQSolverGPU::solveWithDeviceRHS() — nach potrs:
{
    // k_eeq_reduce_sums: berechnet S = sum(Z2[0..N-1]), Cz1 = sum(z1[0..N-1])
    k_eeq_reduce_sums<<<1, 256, 0, eeq_stream>>>(
        impl.d_z1.ptr, impl.d_Z2.ptr, N,
        impl.d_S.ptr, impl.d_Cz1.ptr  // Device-Skalare
    );

    // k_eeq_schur_nfrag1: charges[i] = z1[i] - Z2[i] * lambda
    // lambda = (Cz1 - rhs_constraints[0]) / S
    // Berechnung von lambda entweder via kleinen Kernel oder cudaMemcpy + CPU-Division
    // (Division von 2 Skalaren — akzeptabel als kleiner CPU-Schritt oder als Scalar-Kernel)
    double S_host, Cz1_host;
    cudaMemcpyAsync(&S_host, impl.d_S.ptr, sizeof(double), cudaMemcpyDeviceToHost, eeq_stream);
    cudaMemcpyAsync(&Cz1_host, impl.d_Cz1.ptr, sizeof(double), cudaMemcpyDeviceToHost, eeq_stream);
    cudaStreamSynchronize(eeq_stream);  // Nur 2 doubles — fast
    double lambda = (Cz1_host - rhs_constraints_host[0]) / S_host;

    k_eeq_schur_nfrag1<<<gridFor(N), 256, 0, eeq_stream>>>(
        impl.d_z1.ptr, impl.d_Z2.ptr, lambda, impl.d_charges.ptr, N
    );
    // d_charges ist jetzt der EEQ-Ladungsvektor auf Device
}
```

**Neuer Getter**: `getDeviceChargesPtr()` in `EEQSolverGPU` gibt `impl.d_charges.ptr` zurück.

### Phase WP5-B: CNF-Derivate vollständig auf GPU (3–4h)

`setCNDerivatives(cn, cnf, dcn)` in `gfnff_gpu_method.cpp` ruft `m_gpu_workspace->setCNDerivatives()` auf, was einen H2D-Upload des CNF-Vektors auslöst. Nach WP2 liegt `d_cnf` bereits auf GPU. Anpassung:

```cpp
// setCNDerivatives() ersetzen durch:
m_gpu_workspace->useCNFFromTopologyBuffer();  // d_cnf bereits da, kein Upload nötig
```

Das spart den CNF-H2D-Upload (~11 KB für N=1410) pro Schritt.

### Phase WP5-C: recordD4CNValues auf GPU (3–4h)

```cpp
// Neuer Device-Buffer: d_cn_ref_for_d4 (N doubles, einmalig alloziert)
// Nach jedem Schritt: D2D copy: d_cn_final → d_cn_ref_for_d4
// Neuer Kernel k_check_dc6dcn_skip vergleicht beide:
__global__ void k_check_dc6dcn_skip(
    int N,
    const double* d_cn_cur, const double* d_cn_ref,
    double threshold_sq, int* d_skip_flag);
// D2H download des 1-int-Flags asynchron → beeinflusst nächsten Schritt, nicht aktuellen
```

### Phase WP5-D: finalizeCNForCPU() entfernen (1h)

Nach WP5-A+B+C: Der einzige Verbraucher von CPU-seitigem CN ist `prepareCNAndEEQ(skip_eeq=false)` als Fallback. Im normalen Pfad (nfrag=1, GPU-EEQ) kann `finalizeCNForCPU()` entfernt werden.

```cpp
// gfnff_gpu_method.cpp — WP5-D:
// VORHER:
m_gpu_workspace->finalizeCNForCPU(m_gpu_cn_final);
m_gfnff->prepareCNAndEEQ(gradient, true, &m_gpu_cn_final, true);

// NACHHER (kein CPU-Sync):
// (CNF via Topology-Buffer, CN-Derivate via GPU-Kernel)
if (gradient) {
    m_gpu_workspace->launchCNDerivativeKernels(stream);  // GPU-seitig
}
```

### Phase WP5-E: Vollständiger Schritt-Graph (3–5h)

Wenn WP4 (Phase-2-Graph) und WP5-A bis WP5-D erledigt sind, kann der gesamte Schritt in einem Graph erfasst werden:

```
Graph "full_step":
  Knoten 1: k_cn_compute_pairs (CN)
  Knoten 2: k_dlogdcn (log CN)
  Knoten 3: k_build_eeq_rhs (RHS, depends on 2)
  Knoten 4: k_gaussian_weights (depends on 2)
  Knoten 5: k_dc6dcn_per_pair (depends on 4)
  Knoten 6: Bonded kernels (depends on 1)
  Knoten 7: k_eeq_build_matrix (depends on 2, lazy: skip if !force_refactor)
  Knoten 8: potrf (depends on 7, lazy)
  Knoten 9: potrs (depends on 8 or 7 wenn cached)
  Knoten 10: k_eeq_reduce_sums + k_eeq_schur_nfrag1 (depends on 9)
  Knoten 11: k_coulomb + k_coulomb_postprocess (depends on 10)
  Knoten 12: Gradient-Kernel (depends on 6, 11)
  Knoten 13: cudaMemcpyAsync energy+gradient (depends on 12)
```

**Problem mit Lazy-Refaktorisierung im Graph**: Der `force_refactor`-Flag entscheidet, ob Knoten 7+8 ausgeführt werden. CUDA-Graphs unterstützen bedingte Knoten seit CUDA 12.4 (`cudaGraphConditionalNode`). Für ältere CUDA-Versionen: zwei Graphs (mit/ohne potrf) und bedingte Wahl.

---

## CUDA-Versions-Anforderungen

| Feature | CUDA-Version |
|---------|-------------|
| `cudaStreamBeginCapture` | ≥ 10.0 |
| Multi-Stream-Capture | ≥ 10.1 |
| Graph-Param-Updates | ≥ 11.1 |
| Bedingte Graph-Knoten | ≥ 12.4 |
| Thread Block Clusters | ≥ 11.8 (SM_90+), ≥ 12.0 (SM_120) |

RTX 5080 → CUDA 12.8+, alle Features verfügbar.

---

## Verifikation

Stufenweise nach jeder Phase validieren:

```bash
# Nach jeder Phase:
ctest --output-on-failure
./curcuma -sp polymer.xyz -method gfnff -gpu cuda  # Energie-Regression

# Vollständige Pipeline (WP5-E):
nvprof ./curcuma -md polymer.xyz -method gfnff -gpu cuda -steps 100
# Erwartung: GPU Utilization >80% (vs. aktuell ~20-30%)
# nsys profile für detaillierte Timeline-Analyse:
nsys profile -o profile.nsys-rep ./curcuma -md polymer.xyz -method gfnff -gpu cuda -steps 10
nsys stats profile.nsys-rep
```

---

## Rückwärtskompatibilität

- Fallback-Pfad (CPU EEQ, nfrag>1, Cholesky failure) **bleibt erhalten**
- `m_skip_phase2`-Flag funktioniert weiterhin
- Alle CTests müssen nach jeder Phase grün bleiben
