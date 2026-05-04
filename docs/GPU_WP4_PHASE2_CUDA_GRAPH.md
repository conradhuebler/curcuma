# WP4: Phase-2-CUDA-Graph-Capture

**Kategorie**: Mittlerer Aufwand  
**Aufwand**: ~4–8 Stunden  
**Wirkung**: Mittel — eliminiert ~35 µs Kernel-Launch-Overhead + Driver-Sync in Phase 2 pro Schritt  
**Abhängigkeiten**: WP2 empfohlen (fixer d_rhs_atoms Device-Pointer vereinfacht Capture)  
**Status**: 🤖 Geplant

---

## Kontext

`m_graph_phase1` (Phase 1 — charge-independent) ist bereits implementiert und funktionsfähig. Phase 2 (`launchChargeDependentAndFinish()`) startet 7 Kernel auf 3 Streams, synchronisiert via Stream-Events, und lädt dann Energy + Gradient herunter. Diese 7 Kernel-Launch-API-Calls kosten ~35 µs, die Stream-Events-Koordination weitere ~10 µs — bei kleinen Systemen (N≤2000) ist das ein signifikanter Anteil.

**CUDA-Graphs** ersetzen den wiederkehrenden `kernelA<<<>>>(); kernelB<<<>>>()` Aufruf durch einen einzigen `cudaGraphLaunch()`, den der CUDA-Driver direkt von einer vorab optimierten DAG ausführt.

**Kernbedingung**: Das Graph-Capture funktioniert für Phase 2 weil alle Eingabe-Pointer (Koordinaten, Pair-Indizes, EEQ-Ladungen) in fixen Device-Adressen liegen, die sich zwischen Schritten nicht ändern — nur die *Werte* ändern sich, nicht die Pointer.

---

## Struct-Erweiterungen in ff_workspace_gpu.h

```cpp
// In FFWorkspaceGPUImpl (analog zu bestehenden graph_phase1 Feldern):

// Phase 2 CUDA Graph
bool            m_graph_phase2_valid    = false;
cudaGraph_t     m_graph_phase2          = nullptr;
cudaGraphExec_t m_graph_exec_phase2     = nullptr;
// Invalide wenn: EEQ-Ladungen Pointer sich ändert (unwahrscheinlich),
//               oder wenn Topology rebuilt (SoA-Pointer ungültig).
```

---

## Capture-Logik in launchChargeDependentAndFinish()

```cpp
double FFWorkspaceGPU::launchChargeDependentAndFinish(bool gradient)
{
    auto& impl = *m_impl;
    cudaStream_t stream = impl.stream;

    // WP4: Graph-Replay wenn gültig
    if (impl.m_graph_phase2_valid) {
        cudaError_t err = cudaGraphLaunch(impl.m_graph_exec_phase2, stream);
        if (err == cudaSuccess) {
            // Async-Download (bereits im Graph erfasst): Warten auf Stream
            cudaStreamSynchronize(stream);
            return impl.h_energy[0];
        }
        // Graph-Launch fehlgeschlagen → Graph invalidieren, fallthrough zu normaler Ausführung
        cudaGraphExecDestroy(impl.m_graph_exec_phase2); impl.m_graph_exec_phase2 = nullptr;
        cudaGraphDestroy(impl.m_graph_phase2); impl.m_graph_phase2 = nullptr;
        impl.m_graph_phase2_valid = false;
    }

    // Erste Ausführung oder nach Invalidierung: Capture starten
    bool do_capture = !impl.m_graph_phase2_valid && !impl.m_force_single_stream_phase2;

    if (do_capture) {
        cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
        // WICHTIG: Alle Capture-Streams (sA, sB, sC) müssen ebenfalls captured werden.
        // Dafür müssen sie über cudaEventRecord + cudaStreamWaitEvent mit stream verbunden sein
        // BEVOR BeginCapture aufgerufen wird, ODER man verwendet cudaStreamCaptureModeRelaxed
        // und fügt alle drei Streams manuell hinzu.
    }

    // === Normaler Launch-Code (unverändert) ===
    launchPhase2Kernels(gradient, stream);  // Bestehende Logik auslagern

    if (gradient) {
        k_zero_double<<<...>>>(impl.grad.ptr, 3 * impl.natoms, stream);
        launchGradientKernels(gradient, stream);
    }

    // Async Download (MUSS im Capture drin sein)
    cudaMemcpyAsync(impl.h_energy, impl.d_energy, sizeof(double),
                    cudaMemcpyDeviceToHost, stream);
    if (gradient) {
        cudaMemcpyAsync(impl.h_grad, impl.d_grad, 3 * impl.natoms * sizeof(double),
                        cudaMemcpyDeviceToHost, stream);
    }

    if (do_capture) {
        cudaError_t end_err = cudaStreamEndCapture(stream, &impl.m_graph_phase2);
        if (end_err == cudaSuccess && impl.m_graph_phase2) {
            cudaError_t inst_err = cudaGraphInstantiate(
                &impl.m_graph_exec_phase2, impl.m_graph_phase2, NULL, NULL, 0);
            if (inst_err == cudaSuccess) {
                impl.m_graph_phase2_valid = true;
            } else {
                cudaGraphDestroy(impl.m_graph_phase2);
                impl.m_graph_phase2 = nullptr;
            }
        }
    }

    cudaStreamSynchronize(stream);
    return impl.h_energy[0];
}
```

---

## Multi-Stream-Capture-Problem

Phase 2 nutzt 3 Streams (main, sA für Dispersion, sB für Coulomb). CUDA-Graph-Capture mit mehreren Streams erfordert, dass alle Streams im "capture mode" sind, wenn `cudaStreamBeginCapture` auf einem gestartet wird.

**Lösung — Fork-Join-Pattern**:

```cpp
// Vor Capture: alle Sub-Streams via Events mit main stream verbinden
cudaEvent_t fork_event, join_event_a, join_event_b;
cudaEventCreate(&fork_event);
cudaEventCreate(&join_event_a);
cudaEventCreate(&join_event_b);

// Beim Capture:
cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);

// Fork: sA und sB warten auf stream
cudaEventRecord(fork_event, stream);
cudaStreamWaitEvent(impl.sA, fork_event);
cudaStreamWaitEvent(impl.sB, fork_event);

// Kernels auf allen Streams starten...
k_dispersion<<<..., impl.sA>>>();
k_coulomb<<<..., impl.sB>>>();
k_bonds<<<..., stream>>>();

// Join: stream wartet auf sA und sB
cudaEventRecord(join_event_a, impl.sA);
cudaEventRecord(join_event_b, impl.sB);
cudaStreamWaitEvent(stream, join_event_a);
cudaStreamWaitEvent(stream, join_event_b);

// Download nach Join
cudaMemcpyAsync(impl.h_energy, ..., stream);

cudaStreamEndCapture(stream, &impl.m_graph_phase2);
```

**Alternative (einfacher)**: Für Phase 2 auf 1 Stream umstellen. Die 3 Streams in Phase 2 sparen bei N=1410 kaum Zeit (Kernel sind schnell), der Fork-Join-Overhead übersteigt den Nutzen. Einfacherer Graph mit 1 Stream.

---

## Graph-Invalidierung

```cpp
void FFWorkspaceGPU::invalidateGraph()
{
    // Bestehende Phase-1-Invalidierung (bereits vorhanden):
    if (m_impl->m_graph_exec_phase1) {
        cudaGraphExecDestroy(m_impl->m_graph_exec_phase1);
        m_impl->m_graph_exec_phase1 = nullptr;
    }
    if (m_impl->m_graph_phase1) {
        cudaGraphDestroy(m_impl->m_graph_phase1);
        m_impl->m_graph_phase1 = nullptr;
    }
    m_impl->m_graph_phase1_valid = false;

    // WP4: Phase-2-Graph ebenfalls invalidieren
    if (m_impl->m_graph_exec_phase2) {
        cudaGraphExecDestroy(m_impl->m_graph_exec_phase2);
        m_impl->m_graph_exec_phase2 = nullptr;
    }
    if (m_impl->m_graph_phase2) {
        cudaGraphDestroy(m_impl->m_graph_phase2);
        m_impl->m_graph_phase2 = nullptr;
    }
    m_impl->m_graph_phase2_valid = false;
}
```

`invalidateGraph()` wird bereits bei Topologie-Rebuild aufgerufen → Phase-2-Graph wird automatisch neu aufgebaut.

---

## Einschränkungen und Edge Cases

**gradient=false vs. gradient=true**: Zwei getrennte Graphs nötig (unterschiedliche Kernel-Sequenzen). Der aktuelle Phase-1-Code hat das identische Problem und löst es mit `need_snapshots` Flag. Analoges Pattern für Phase 2:

```cpp
// Statt einem m_graph_phase2:
cudaGraphExec_t m_graph_exec_phase2_energy  = nullptr;  // gradient=false
cudaGraphExec_t m_graph_exec_phase2_gradient = nullptr; // gradient=true
bool            m_graph_phase2_energy_valid  = false;
bool            m_graph_phase2_gradient_valid = false;
```

**EEQ-Ladungs-Upload**: `setEEQCharges()` muss **vor** `cudaGraphLaunch()` aufgerufen werden, da der Graph die Daten liest aber nicht kopiert. Die Daten liegen in `impl.d_charges.ptr` — dasselbe Device-Array, das beim Capture gelesen wird. Solange der Graph nicht den Upload selbst enthält (kein `cudaMemcpy` für Ladungen im Graph), ist die Sequenz:
```
setEEQCharges() → H2D cudaMemcpyAsync auf stream
cudaGraphLaunch(m_graph_exec_phase2, stream)  → liest die gerade geschriebenen Ladungen
```
Das funktioniert, weil cudaGraphLaunch auf demselben Stream nach der Memcpy sequenziert ist.

---

## Verifikation

```bash
# Build:
cd release && make -j4

# Erste Ausführung: Graph wird gecaptured (Schritt 1)
# Zweite Ausführung: Graph wird replayed → im Log sollte "phase2 graph replay" erscheinen
./curcuma -md ../test_cases/cli/simplemd/10_gfnff_polymer_md/input.xyz \
  -method gfnff -gpu cuda -steps 10 -verbosity 2

# Energie-Regression:
./curcuma -sp ../test_cases/molecules/larger/polymer.xyz -method gfnff -gpu cuda
# Vergleich mit vorherigem Wert

# Timing:
./curcuma -md ../test_cases/cli/simplemd/10_gfnff_polymer_md/input.xyz \
  -method gfnff -gpu cuda -steps 100 -verbosity 1
# Erwartung: Phase-2-Zeit reduziert von ~5 ms auf ~0.05 ms (Kernel-Launch-Overhead entfernt)
# Compute-Zeit der Kernel selbst bleibt gleich
```
