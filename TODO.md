# Curcuma Development TODO List

**Stand**: 2026-10-01. Die frühere Fassung (1380 Zeilen, mit Messprotokollen und erledigten Untersuchungen) steht
unverändert in [docs/TODO_ARCHIVE_2026-10.md](docs/TODO_ARCHIVE_2026-10.md); die Detailmessungen zu den Einträgen unten
sind dort unter der jeweiligen Überschrift zu finden.

Regeln für diese Datei: ein Eintrag hat höchstens drei Zeilen (Was, warum es offen ist, wo die Details stehen). Erledigtes
gehört in `AIChangelog.md` bzw. das Archiv, Fehlerberichte in [docs/KNOWN_ISSUES_ARCHIVE.md](docs/KNOWN_ISSUES_ARCHIVE.md).
"Geprüft" heißt: am angegebenen Datum gegen den Quelltext gelesen, nicht neu gemessen. Dateien und Zeilen gelten für diesen
Stand.

## 1. Korrektheit (offen)

- **Flache SCF-Namen, die auch das Optimiermodul kennt** (2026-10-26, Ursache gefunden): `-diis_start`, `-max_iterations` gehören auch zu `opt`, das `-sp` und `-opt`
  benutzen; das Befehlsmodul gewinnt, also erreichen sie den SCF nicht (keine Warnung). Wirksam: `-scf_max_iterations`, `-scf.diis_start`, `-xtb.diis_start`.
  Offen: soll `-sp` die Optimiererparameter überhaupt besitzen? [docs/PARAMETER_SYSTEM.md](docs/PARAMETER_SYSTEM.md), "Shared module scf".
- **SCF-Modul `scf` bei tblite/xtb/Ulysses nicht gelaufen** (2026-10-26): Abbildung auf `max_iterations`/`damping` nur per `-verbosity 3` gesehen, die
  TBLite-Änderung (GFN2-Anhebung auf 150/0.6 nur ohne Nutzerwert) nur syntaxgeprüft; kein Provider in `release/`. `-hessian` endet mit 0, wenn ein
  Verschiebungs-SCF scheitert (pm3, `-scf_max_iterations 3`), nicht untersucht. [docs/SCF_MODES.md](docs/SCF_MODES.md)
- **Offene Schalen über die CLI: Aliase und ctest fehlen** (2026-10-11): `-multi`, `-multiplicity`, `-uhf`, `-xtb.spin` werden still ignoriert; Zeile 2 "Ladung Multiplizität" der ORCA-xyz wird nicht gelesen;
  ungerade Elektronenzahl ohne `-spin` besetzt fraktional (tblite), xtb nimmt ein Dublett (Entscheidung offen); es gibt keinen ctest für `-spin` (Known Issue #9). Flag mit Leerzeichen und Paritätswiderspruch
  sind seit 2026-10-11 Fehler bzw. Warnung (`CLI2Json`, `buildReferenceOccupations`). Vault: Labor `curcuma GFN2-SCF-Vorkonditionierung`, Schritt 13.
- **GFN-FF-Energie hängt von der Atomnummerierung ab** (911 von 2557 Strukturen, bis 237 kcal/mol bei geladenen Mehrfragment-Systemen;
  in pprcht und curcuma). Referenzeigenschaft, kein Portierungsfehler. [docs/REV_GFNFF_TODO.md](docs/REV_GFNFF_TODO.md) #12, Known Issue #34.
- **GPU-Schrittzeit GFN-FF, Reste und nicht ROCm-spezifische Hebel** (Stand 13. Okt. 2026, Messwerte in [docs/ROCM_PERFORMANCE_ANALYSIS_2026-10-13.md](docs/ROCM_PERFORMANCE_ANALYSIS_2026-10-13.md)): (a) MD-Integrator und Thermostat auf der GPU statt `SimpleMD` auf der CPU (7320 Atome: 4-7 ms Luecke und 4,3 MB Transfer pro Schritt, etwa 10 %); (b) viele kleine Moleküle gebündelt auf der GPU (Konformerensuche, Durchsatz); (c) Coulomb (10 ms bei 7320 Atomen) und Dispersion (7,5 ms) weiter ausreizen, auch 20 % lohnen sich; (e) `-gpu auto` samt Schwelle entfernt (2026-10-26); Break-even Wasser auf ROCm zwischen 300 und 600 Atomen, andere Systemtypen offen; CUDA gemessen 17. Okt. 2026 ([GPU_TUNING.md](docs/GPU_TUNING.md)): MD-Schritt ab etwa 66 Atomen schneller, aber 0,8 s Start je Prozess, `-gpu auto` hat auf CUDA keine Schwelle (Entscheidung offen); (f) projiziertes PCG macht bei 600 Atomen ca. 11 Host-Roundtrips je Schritt, Abhilfe: Konvergenzflag auf dem Device; (g) Logzeile sagt bei kleinem Device-EEQ noch "WP2 dense Cholesky + CPU Schur" (kosmetisch); (d) danach GFN1/GFN2 auf ROCm (Host-Eigenloeser unterhalb eines Break-evens, Cache-L, Pseudo-Diagonalisierung), erst nachdem GFN-FF ausgereizt ist.
- **Periodisches GFN-FF fehlt** (zurueckgestellt, Entscheidung Okt. 2026): die CPU nutzt das Gitter nirgends, die ROCm-Kernel haben `applyMIC` und `setUnitCell`, aber `unitcell 9 9 9` aendert die Energie von `water8-cluster` weder auf CPU noch auf ROCm (8 Stellen). Plan und Stand: [docs/GFNFF-CPU.md](docs/GFNFF-CPU.md), [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md) ("Periodic boundary conditions: not implemented").
- **Frag-Ladungs-Ensemble auf der GPU**: geladene Mehrfragment-Systeme laufen seit 14. Okt. 2026 auf der CPU-Methode (GPU-Rueckfall, gezaehlt); das Ensemble selbst (`GFNFF::fragEnsembleBlend`) laeuft nur dort. Auf der GPU abbilden waere die Alternative, [docs/FRAG_CHARGE_MODEL.md](docs/FRAG_CHARGE_MODEL.md).
- **ROCm**: GFN1/GFN2 ignorieren `-scf_mixed_precision false` (G2-13); Energiesummen nehmen wave32 an (F-17). Dokumentiert in
  [docs/MULTI_GPU_GAPS.md](docs/MULTI_GPU_GAPS.md). Leistung (EEQ auf der CPU, rocSOLVER-Eigenloeser):
  [docs/ROCM_PERFORMANCE_ANALYSIS_2026-10-13.md](docs/ROCM_PERFORMANCE_ANALYSIS_2026-10-13.md).
- **Exitcode bei fehlgeschlagenem Setup**: `-md` gibt seit `a2b7c5f0` 1 zurück; prüfen, ob dasselbe für `-opt` und `-sp` gilt (nicht geprüft).
- **UFF- und QMDFF-Gradient falsch** (Eintrag vom 2026-09-27: UFF-Winkel ohne den Faktor -sin(theta); Fix vorbereitet, Betreiber: "erstmal nicht").
  Nicht erneut geprüft; Details im Archiv unter "UFF- und QMDFF-Gradient falsch".
- **`-seed` wirkt nicht auf die Anfangsgeschwindigkeiten der MD.** Geprüft 2026-10-01: `SimpleMD::InitVelocities()` (`simplemd.cpp:1141`)
  zieht aus `static std::default_random_engine generator` (Z. 1143) mit festem Standardseed. Ein Fix ändert jede MD-Startgeschwindigkeit.
- **`shouldUpdateHBXB()`: RMSD-Formel ohne Faktor `sqrt(natoms)`** (`gfnff_method.cpp:3058`, geprüft 2026-10-01). Bei 7320 Atomen
  ~13,6 Å mittlere Verschiebung bis zum Auslösen. Aus der Fortran-Referenz übernommen (`gfnff_ini2.f90:717`); vor einer Änderung gegen
  MOR41/GMTKN55 prüfen und das Abweichen von der Referenz entscheiden.
- **D4-ATM-C6 auf den Tripeln bleibt auf dem Setup-Wert eingefroren.** Ohne Wirkung, solange `dispersion_atm` aus ist (Voreinstellung).
  Ein Paar, das erst später in 60 Bohr eintritt, wurde nicht durch einen Test erzeugt. [docs/GFNFF_PAIR_LIST_REFRESH.md](docs/GFNFF_PAIR_LIST_REFRESH.md).
- **Optimierer**: (a) der eigene L-BFGS bleibt bei Koffein/GFN-FF ab Schritt 33 stehen (die Stillstandserkennung `-opt.stall_steps` beendet
  den Lauf inzwischen, Ursache nur per Codelesung); (b) `mixture2`: kein Abstieg an einer N-H···O=C-Brücke, nicht eingegrenzt;
  (c) der `auto`-Optimierer kann bei hochsymmetrischen Systemen mit nur totalsymmetrischer Mode abbrechen (Td-CH4, gfn1; `-opt.optimizer lbfgs` geht).
- **GFN-FF-NVE-Drift auf `complex`** (231 Atome): 6,1e-3 Eh über 1 ps bei dt 0,5 fs, unabhängig vom EEQ-Löser. Nicht untersucht. Auch `-gpu cuda` polymer_2x NVE 400 fs: Etot -0,0548 Eh in allen Läufen, schon vor den Optimierungen vom 2026-10-11 (n=2 je Binary).
- **SIGSEGV am Ursprung nicht gefunden** (2026-09-24; ASan, compute-sanitizer, gdb ohne Treffer). Eintrag im Archiv.
- **`DC13/c20bowl`**: Coulomb -6,29 kcal/mol gegen die Referenz, Ladungen ~6x zu groß, nicht eingegrenzt (Known Issue #14).
- **Reproduzierbarkeit**: GPU läuft nicht bitgleich (`atomicAdd`-Reihenfolge); CPU nicht bitgleich über verschiedene Threadzahlen
  (`FFWorkspace` partitioniert nach Threadzahl). Betreiber 2026-09-27: später. Erster Schritt: Repulsion als Sammel-Kernel hinter
  einem Schalter, messen auf A4500 und H200. Archiv: "CPU/GPU-Trajektorien-Divergenz" und "Bitgleichheit". Stand 17. Okt. 2026: GFN1/GFN2-SCF auf der GPU bitgleich; GFN-FF-CUDA-Ladungs- und Kraftpfad ebenfalls (Fixkomma, 9 von 9 Trajektorien identisch); Gradienten von GFN1/GFN2 ebenfalls; Energien ebenfalls; offen sind ROCm und Split gegen Ein-GPU-Lauf (1e-12, FP32 + anderer Kernel), Optionen in [docs/GPU_REPRODUCIBILITY.md](docs/GPU_REPRODUCIBILITY.md).
- **Verbosity** ist ein geteilter Static und lässt sich über `CxxThreadPool`-Worker nicht sauber scopen (Known Issue #3).
- **ConfSearch**: Wide-Hill-MTD-Aufblähung innerhalb eines Laufs; `threads=1` im Legacy-Modus mit Fortschrittsbalken (SIGSEGV nach der
  Optimierung, Eintrag von 2025, nicht erneut geprüft). [docs/CONFSEARCH_ROADMAP.md](docs/CONFSEARCH_ROADMAP.md).
- **Windows**: `-DUSE_PORTABLE_MATH=ON` nicht auf einer echten Windows-Maschine bestätigt ([docs/PORTABLE_ERF.md](docs/PORTABLE_ERF.md)).
- **Ladungsplatzierung** (Known Issue #31) nicht getestet für periodische Systeme, GPU-Laufzeit, viele Fragmente in Kontakt und lange MD mit
  `frag_charge_s_max` > 1.
- **Fensteraufloesung ueber steilen Barrieren**: NEB-MD integriert dG/ds mit Trapez/Spline ueber `nimages` Fenster; Retro-DA GFN2: 12/23/34 Fenster geben dG_TST 91.5/71.8/67.0
  (n=1). Noetig: Warnung bei grossem Sprung von dG/ds zwischen Nachbarfenstern, adaptive Verdichtung am Maximum, `nimages` mit `path_file` nicht ignorieren. [docs/NEB_MD.md](docs/NEB_MD.md).
- **Restraint-Modus-Tube** (`umbrella_kz`) laeuft ueber `applyExternalForces` (einen Schritt verzoegert, nicht in `trueGradient()`), ueberschneidet sich mit `z_wall`.
- **`Molecule::CalculateDihedral` falsch** (`src/core/molecule.cpp` ~551, benutzt von `-dihedral` in `src/main.cpp` ~523): das Skalarprodukt der unnormierten Flaechennormalen wird vor der Division durch die Betraege auf [-1,1] begrenzt
  (Agent-Nachrechnung: 116.4 Grad fuer eine +60-Grad-Geometrie), Ausgabe 0..2pi statt IUPAC. Korrekt ist `GFNFF_Geometry::calculateDihedralAngle`. Nicht repariert, nicht gegengeprueft ausser der Rechnung des Agenten.
- **MD-Startgeschwindigkeiten und Andersen**: ein prozessweit geteilter statischer Generator, Seeds mehrerer SimpleMD-Instanzen nicht unabhängig reproduzierbar
  (NEB-MD-Fenster). Vault: Projekte-Notiz `curcuma NEB-MD Pfad-Sampling`.

## 2. Leistung

- **CUDA/ROCm-Parität: ROCm-Seite** (2026-10-24, [docs/GPU_PARITY_CUDA_ROCM.md](docs/GPU_PARITY_CUDA_ROCM.md), Abschnitt A): auf ROCm fehlen residente SCF-Schleife, Geräte-Potential, Multipol on-the-fly, Solvation, Besetzungen/Broyden/Pseudodiagonalisierung, D3/D4 nach dem SCF, Gradientenschwanz; GFN-FF: Geräte-Schur, Mehr-GPU-Split, Geräte-Abstoßungsliste, tabelliertes dc6dcn, Fixed-Point.
- **CUDA/ROCm-Parität: CUDA-Seite** (Stand 2026-10-25, [docs/GPU_PARITY_CUDA_ROCM.md](docs/GPU_PARITY_CUDA_ROCM.md), Abschnitt B): erledigt 1, 10, 12 (Faktor-Probe, Stream-Reihenfolge, Fixkomma-Flag); gemessen und nicht portiert 2, 3, 4, 5, 6, 7, 8; offen nur 9 ((entschieden 2026-10-26: `-gpu auto` entfernt)) und 5 bei Einzelkarten ueber 5000 Atomen.
- **Erschoepfte SCF meldet 151 statt 150 Iterationen** (Zaehlung der Host-Schleife, nicht untersucht). Der D4-Gewichtsfallback auf CUDA ist seit 2026-10-25 erledigt (`k_d4_build_refw`, AIChangelog).
- **ROCm: rocSOLVER liefert bei mehreren Prozessen auf einer GPU sporadisch falsche Ergebnisse** (2026-10-17, [docs/ROCM_GFN_POTENTIAL_2026-10-14.md](docs/ROCM_GFN_POTENTIAL_2026-10-14.md), Reproduzierer `scripts/rocsolver_concurrency_check.cpp`): `potrf` 6 und `sygvd` 8 von je 12000 Aufrufen bei 4 Prozessen, 0 bei einem. Offen: Meldung an AMD, anderer Treiber/ROCm-Stand, Wirkung im `sygvd`-Standardpfad, Messung des Caches mit Probe gegen den Standardpfad, CUDA-Gegenprobe.
- **GFN1/GFN2 auf ROCm** (Stand 2026-10-15, [docs/ROCM_GFN_POTENTIAL_2026-10-14.md](docs/ROCM_GFN_POTENTIAL_2026-10-14.md)): erledigt sind Dichte-GEMM über ncol, gecachter Cholesky-Faktor, syrk (Opt-in) und explizites L^-1 (Opt-in seit 2026-10-17, vorher Default; peo201-chain 14.18 -> 12.85 s, n=3) mit einmaligem Fallback-Hinweis; Multipol-Potential/-Energie auf dem Gerät als Opt-in `-xtb.rocm_device_multipole` (2026-10-17, peo201-chain 14.92 -> 14.38 s, n=3).
  Offen: Host-Threads auf ROCm Default 1 (CUDA 8); Entscheidung, ob `-xtb.rocm_device_multipole` Default wird (Opt/MD, Solvens, Fallback nicht gelaufen); restliche Host-Lücke pro Iteration; D4-ATM 1.2 s; späterer FP64-Wechsel (-19 %, n=1); Pseudo-Diagonalisierung wie CUDA.
  CUDA-`w_from_goo` (`xtb_gpu_context.cu`) nimmt noch geschlossene Schale an; CUDA/Vulkan mit der Open-Shell-Gradientenkorrektur nicht gelaufen.
- **Host-Rechnung im GPU-Lauf, Rest** (Stand 17. Okt. 2026, zweite Runde): auf der Karte sind D4/ATM, D3, Multipol-Matrizen, Abstoßung, Gradient-Schwanz, Gradient-Potentiale, D4-CN, Hamilton-CN; S/H0/L/gamma werden nur auf Anforderung auf den Host kopiert (Tabelle in docs/GPU_TUNING.md). Offen: ROCm/Vulkan unberührt; Verbosity-3-Stelle (GFN1, 9,5 s, nicht gefunden); PR18 (MOR41) GFN2-Gradient springt selten um 4e-7 (GPU-Nichtreproduzierbarkeit, auch Host-gegen-Host). Mit `-threads 1` gleich schnell wie `-threads 16`: nichts davon skaliert.
- **GFN2 auf der GPU, Eigenlöser = 73 % des Laufs** (polymer_2x, 1x A4500). Offen:
  (a) erledigt 15. Okt. 2026: `-scf_pseudo_diag` ist Standard (Betreiberentscheidung); polymer_2x 4x A4500 113,6 -> 90,0 s;
  (b) `-scf_pseudo_diag_fp64` auf einer GPU mit vollem FP64 (H200) messen;
  (c) CPU-Pfad `-scf_pseudo_diag_cpu` erledigt 15. Okt. 2026 (opt-in, EXPERIMENTAL): polymer GFN2 -9 %, GFN1 -19 %, [docs/GFN2_GPU_COST_PLAN.md](docs/GFN2_GPU_COST_PLAN.md); offen: Default-Entscheidung, MOR41/GMTKN55-Lauf, FP32-GEMMs für die FP32-Phase;
  (d) Pseudo-Schritt auf mehrere GPUs verteilen: nach Messung verworfen (Schritt 2,3 s von 3,4 s je Iteration, geschätzter Gewinn 5 s von 90),
  (d) Teilspektrum: der residente Loop übergibt `n_eig = 0`, Grund unbekannt, schließt sich derzeit mit der verteilten FP64-Lösung aus;
  (e) Schwellen `gpu_eigensolver_min_nao` (4000) und `gpu_density_min_nao` nachmessen;
  (f) `purify`/`lobpcg` gibt es nur auf der CPU, auf der GPU ist es der einzige Hebel, der die Skalierung ändert;
  (g) Setup ist profiliert (~13 s, nsys, GFN2_GPU_COST_PLAN.md 0.3); Integrale/Setup nur bei Bedarf weiter.
  **`scf_fp32_stall_patience` nicht ändern**, bis die Fixpunktfrage (Schritt 1 weicht um 2,7e-6 Eh ab) geklärt ist.
  Zahlen und Begründungen: Archiv, Abschnitt "GPU-SCF GFN1/GFN2"; [docs/GFN2_GPU_COST_PLAN.md](docs/GFN2_GPU_COST_PLAN.md).
- **GPU-Gradient GFN1/GFN2**: Punkte 1 bis 6 umgesetzt; offen: der H0-Paar-Kernel (454 ms) ist einzelkartig, und `densityPatternDistributed`
  kopiert die C-Spaltenscheiben bei jedem Aufruf neu. [docs/SQM_PERFORMANCE.md](docs/SQM_PERFORMANCE.md).
- **GFN-FF-MD, EEQ-Löser** (polymer_2x): offen sind ein matrixfreier Matvec nur für Systeme > 20000 Atome sowie localinv mit offener Schale, `eeq_distance_cutoff` > 0 und MD über 2 ps
- **`localinv` bei 32 Threads (SMT) langsamer als bei 16** (2026-10-13, polymer_2x 80 auf 176 ms, urea400 34 auf 106 ms, n=3; M-Aufbau ist es nicht); Ursache offen. Messung: Vault-Labor `curcuma EEQ-Solver Benchmark`.
  (messen). Erledigtes (lokale Inverse, symmetrisches Matvec, Fernfeld-Kurzschluss CPU/CUDA) steht in `AIChangelog.md` und [docs/EEQ_SOLVER_CONTROL.md](docs/EEQ_SOLVER_CONTROL.md);
  Journal: Vault-Laborjournal "curcuma EEQ-Solver Benchmark".
- **ROCm rechnet Host-`D4DC6DCNPairs` ohne Leser**: mit `gpu_disp_pairs_on_device true` (CUDA-Standard) entfällt die Host-Liste (polymer_2x, Spitzen-RSS 3,68 GB statt 3,91 GB mit `false`, je 1 Lauf, 17. Okt. 2026); ROCm hat den Schalter aus, HIP liest `m_disp_idx_*_host` und füllt `m_h_dc6dcn_*` noch (nicht kompilierbar hier).
- **GPU-GFN-FF-MD, Schrittzerlegung polymer_2x** (Stand 21,7 ms je Schritt, vorher 38,4; erledigt: `AIChangelog.md`, [docs/GFNFF-PERFORMANCE-ROADMAP.md](docs/GFNFF-PERFORMANCE-ROADMAP.md)):
  offen: CG-Iterationen unter Last (K5); Host-Neubau überlappt nicht mit Gerätearbeit (Stand 16. Okt. 2026, 296 Schritte, 20,0 ms/Schritt: HB/XB 2,2 ms/Schritt (31 Refreshs à 25 ms: Erkennung 14,5, Upload 5,1, Bond-HB-Neuaufbau 3; die HB-Liste ändert sich bei jedem Refresh, 49,6k -> 50,4k Paare, Überspringen bei unveränderter Liste greift nicht; Erkennung skaliert schlecht mit Threads, 4/8/16/32: 24,5/21,5/12,8/10,5 ms; `hb_update_force_every` 10 -> 20 spart 1,1 ms/Schritt, ändert aber die Trajektorie, |dE_total| bis 6e-4 über 300 Schritte, n=1; Entscheidung offen), CN-Liste K1 1,0, Repulsion-Geräteneubau 0,9, D4 0,7; Repulsion lief vorher auf dem Host, 3,3 ms/Schritt). ROCm-Rückstand: siehe unten. Journal: Labor/curcuma EEQ-Solver Benchmark.
- **CUDA GFN-FF MD, Schrittpausen urea400-water1000** (Stand 2026-10-14, [docs/GPU_STEP_PAUSES_2026-10-14.md](docs/GPU_STEP_PAUSES_2026-10-14.md)): offen: ROCm-Spiegel der D4-Setup-Reihenfolge (`gfnff_rocm.hip`); D4-Sortierung und HB-Erkennung sind auf ROCm gespiegelt, beide Schalter dort Standard;
  geladene Mehrfragmentsysteme laufen auf der CPU (15-fach langsamer je MD-Schritt bei 2910 Atomen), Ladungsmodell: Vault "curcuma Ladungsträger-Heuristik", Abschnitt GPU-Aspekte. Journal: Labor/curcuma GPU-Backends.
- **ROCm-GFN-FF, Rückstand** (Hardware nötig, groß; Stand 2026-10-17, Codelesung): Geräte-ppcg, localinv, K1 (CN-Paarliste mit Verlet-Skin), Fernfeld und gestaffelter Repulsions-Upload sind vorhanden, K3 ist ein No-op; fehlen: K2 nur teilweise, Geräte-Schur für Mehrfragmentsysteme (dichter Cholesky-Pfad rechnet die Schur auf dem Host), ALPB-Solvatation und PBC auf dem Gerät (Solvens: CPU-Fallback seit 2026-10-17). [docs/EEQ_SOLVER_CONTROL.md](docs/EEQ_SOLVER_CONTROL.md), Host-/Geräteaufteilung: [docs/ROCM_HOST_PARTS_2026-10-17.md](docs/ROCM_HOST_PARTS_2026-10-17.md).
- **GPU-Fernfeld und Toleranz angleichen**: (Fernfeld in `k_eeq_build_cols` steht seit Okt. 2026 drin, mit `rsqrt`); CUDA-Schwelle 6,0 auf 6,5 wie CPU (nicht bitgleich, mit Messung);
  GPU-ppcg auf relatives `eeq_ppcg_tol` statt absolutem `pcg_tolerance` (ändert Ergebnisse).
- **CPU-MD-Schritt polymer_2x** (Stand 212 ms bei 16 Threads, vorher 294; erledigt: `AIChangelog.md`, [docs/EEQ_SOLVER_CONTROL.md](docs/EEQ_SOLVER_CONTROL.md)): offen Coulomb vektorisiert (~6 ms),
  `d4_gw` (E) (~10-15 ms); 1 Thread: `cn` +4,8 ms ungeklärt. Journal: Labor/curcuma EEQ-Solver Benchmark, Abschnitt "CPU-MD-Schritt polymer_2x zerlegt".
- **EEQ-Löser im CPU-MD** (Stand 164 ms je Schritt bei 16 Threads, vorher 240; erledigt: `AIChangelog.md`, [docs/EEQ_SOLVER_CONTROL.md](docs/EEQ_SOLVER_CONTROL.md)): offen beide Energieformen des Plausibilitätstests
  in einem Durchlauf (gemm statt 2 dgemv, 13 ms, Rundung nicht gesichert); Matrix wird jeden Schritt neu gebaut (48 ms). Nicht exakt: eps-Wiederverwendung der Matrix, Extrapolation, lockerere Toleranz.
- **GFN-FF CPU `d4_gw`** (Stand 34 ms je Schritt bei 16 Threads, vorher 103; erledigt: `AIChangelog.md`, [docs/GFNFF-PERFORMANCE-ROADMAP.md](docs/GFNFF-PERFORMANCE-ROADMAP.md)):
  offen (D) C6 und dc6 in einem Durchlauf, (E) C6 erst in `calcD4Dispersion`. Nicht exakt: Paare bei kleiner CN-Änderung auslassen (#32).
- **Genauigkeitsprofile `-accuracy loose|normal|high` (geplant, nicht begonnen, Betreiber 2026-10-10: jetzt nicht umsetzen)**: der Parameter `accuracy` (`gfnff.h`, Modul `gfnff`) wird nur
  durchgereicht, nichts liest ihn. Plan: eine zentrale Auflösung wie bei der Solver-Wahl, `normal` = heutiges Referenzverhalten bitgleich, Einzelschalter überschreiben; zuerst ein Scan der Schalter
  (`eeq_ppcg_tol`, `cn_`/`hb_accuracy`, `eeq_distance_cutoff`, `coulomb_r_cut`, `dispersion_cutoff_bohr`) mit Fehler in E, F, NVE-Drift und Zeit, dann Fehlerbudget je Stufe; Abschneidegrenzen nur mit glatter Auslaufzone und konsistenter EEQ-Matrix. Offen: Stufennamen, Budget, `gfnff` oder `revgfnff`.
- **GFN-FF Static-Mode-Initiative** (Mai 2026): die PARAMs `static_charges`/`static_cn` (WP-S1) sind im Code (`gfnff.h:323`), und der
  Changelog nennt WP-S1 bis S3 als umgesetzt. WP-S4 (Validierungssuite, 32-Konfigurationen-Matrix) nicht geprüft.
  [docs/GFNFF_STATIC_WP4_VALIDATION_SUITE.md](docs/GFNFF_STATIC_WP4_VALIDATION_SUITE.md).
- **Speicher für Systeme > 1000 Atome** (Molecule-Datenstruktur, Distanzmatrix-Cache): geplant, nicht begonnen.
- **SCF-Vorkonditionierer `-scf_precond frag`** (Okt 2026, [docs/SCF_PRECOND.md](docs/SCF_PRECOND.md)), nicht umgesetzt: Ionen und geladene
  Fragmente (nur Skalarschritt), offene Schalen, GFN1 auf der GPU; GPU siehe unten.
- **`-scf_precond` `auto`-Modus** (2026-10-12): Aufbau von M lohnt sich nur bei wiederholter SCF (Census: Einzelpunkt eines Moleküls 1,25x bis 22x langsamer; Break-even etwa 2 SCFs unter 10 Atomen, 8 bei 10-30, 62 bei 30-100). Regel (wiederholte SCF in MD/Opt/Hessian, Fragmente <= 30 Atome) und
  Volkszählung der Schalter-Standardfrage: [docs/SCF_PRECOND.md](docs/SCF_PRECOND.md), Abschnitt Census.
- **`-scf_precond frag`, ungemessen**: Systeme > 600 Atome (CPU) bzw. > 1499 (GPU); M-Erneuerung in langer MD (> 500 fs, `-scf_precond_rebuild_rms`);
  Hessian mit der eingebauten Funktion (jedes Objekt ausgelenkter Geometrien baut sein eigenes M; gemessen nur der Prototyp bis 197 Atome); Fragmente über
  `-scf_precond_max_fragment_atoms` (Standard 100; Gewinn nur bis ~30 Atomen gemessen); Elemente außer H, C, N, O; Solvatation.
- **`-scf_precond_segment_atoms` (RAS), offen** (2026-10-12): gemessen PEO (58, 212, 1410 Atome), Polystyrol 325 und ein Peptid 309 Atome (n=1, 12->5, 15->6/7 Iterationen);
  gemessen auch MD 200 fs und Hessian (M von den Hessian-Workern geteilt); offen: GPU-Hessian baut M weiterhin je Worker, lange MD (> 200 fs), ladungsgetrennte Systeme, Proteine mit Resten als Einheit. [docs/SCF_PRECOND.md](docs/SCF_PRECOND.md)
- **GFN2-MD: Etot springt in den ersten 3 Schritten um 4 bis 6 mEh** (2026-10-12, Peptid/Polystyrol 309/325 Atome, NVE, 1 fs, mit und ohne `scf_precond`); Ursache offen, 200 fs belegen keine Langzeitstabilitaet (Zeitschritt, `-scf_threshold`, laengerer Lauf pruefen).
- **`-scf_precond_jacobian response`** (2026-10-12): GFN2 und GFN1 (CPU), geschlossenschalig, gemessen n=1 an PEO/Polystyrol/Peptid; offen: Standard auf `response` umstellen?, occ-occ/virt-virt-Paare auslassen
  (etwa halbe GEMM-Kosten), analytische Potentialableitung statt Differenzen, weitere Nutzer von J (Hessian, Analyse). [docs/SCF_PRECOND.md](docs/SCF_PRECOND.md)
- **`-scf_precond frag` auf der GPU, offen** (2026-10-11): umgesetzt nur für die vollresidente CUDA-Schleife ([docs/SCF_PRECOND.md](docs/SCF_PRECOND.md), Abschnitt GPU). Mit Warnung ignoriert: ROCm, Vulkan, hostgetriebene CUDA-Schleife
  (Mischer auf dem Host, der CPU-Hook wäre dort direkt nutzbar), Mehr-GPU-SCF-Aufteilung (Eigenlöser/Dichte; ungetestet, nur ein Gerät hier).
  Ohne Wächter divergiert der unrelaxierte 599-Atom-Cluster mit Mischpräzision (151 Iterationen, harter Fehler, n=3).

## 3. Tests und Validierung

- **xtb-Referenzenergien der Teststruktur-Bibliothek** (2026-10-11, `test_cases/reference_data/xtb_library_energies.json`): 66 Strukturen GFN-FF, 64 GFN2 gerechnet (peo201-chain 1410 Atome GFN2: +1,6e-3 Eh, mit `d4_atm_cutoff 40.0` +2,6e-5). GFN-FF: 57 <= 1e-6 Eh, 8 > 1e-4 (hydronium-water2 -0,169 Eh Standard `ensemble` gegen xtb Referenzregel;
  jl22-bgal/macrocycle-bgal 173/231 Atome +0,25 bis +0,80 kcal/mol, nicht erklärt: `storsion_reference_loop_bug` ändert nichts, Term-Vergleich und pprcht/gfnff als dritte Meinung stehen aus). GFN2: max 5,2e-5 Eh, mit `-xtb.d4_atm_cutoff 40.0` auf macrocycle-bgal 5e-5 -> -2e-7 Eh. Offen: ctest gegen die Datei (Toleranzen), 4 Strukturen > 3000 Atome ohne xtb-GFN-FF (xtb wird bei 6200 und 7320 Atomen am 48-GB-Deckel getötet; curcuma 5 s), na2cl2 GFN2 xtb unkonvergiert.

- **NEB-MD Helicen-Neulauf** mit dem Stand nach `e81f1a23` (Trust-Radius, Abbruch bei toten Images): einziger vollständiger Lauf (2026-10-03, n=1) hatte tote Images in Stufe 1
  und Reaktions-dG 2.3-2.5 statt 0; zusätzlich mehrere unabhängige String-Pfade für Ethan (Pfadvarianz ungemessen). Vault: `Projekte/curcuma NEB-MD Pfad-Sampling`.
- `cli_simplemd_08/09` (Essigsäuredimer, CSVR, dt 1 fs) gegen das gemergte Binary mit der MD-Uhr-Korrektur (`ef462fcf`) erneut laufen lassen (Known Issue #32).
- `test_cg_potentials` ist seit 2026-10-01 als ctest `cg_potentials` registriert (besteht); der CG-Beispielaufruf mit VTF-Eingabe endet mit "Failed to initialize ForceField
  engine" (`docs/archive/CG.md`, geprüft 2026-10-01).
- **Zwei ctests scheitern reproduzierbar** (am 2026-10-01 auf `master` und auf `fix/remove-dead-code` gleich, `confscan_dtemplate` 3 von 3 Läufen auf `master`): `xtb_cpscf` (Gate (C), Antwort Mulliken gegen EEQ: H2 8.2e-5 gegen 1e-6, H2O und HCN je über 1.5e-4), `test_orca_interface` ("O2 triplet mult=3" in `methodToOrcaKeyword`), `confscan_dtemplate` (bis 2026-10-09 erwartet 17/2/1/305, tatsächlich 14/5/1/237 unabhängig von `CurcumaThreads` 1/4/24; die Zähler sind jetzt die aktuellen des Binaries und gleich denen von `confscan_free`/`confscan_subspace`; ob `dtemplate` tatsächlich wie `subspace` rechnet oder die alten Werte vor 2026-08 anders zustande kamen, ist nicht geklärt). Die Einordnung "flaky" für `confscan_dtemplate` in älteren Notizen trifft nicht zu.
- **Voller ctest 2026-10-11** (Paket 4, CPU-only ohne `-march=native` und CUDA-Build, je 9 rot, `confscan_dtemplate` grün): `test_orca_interface` scheitert auch bei "O2 quartet mult=5" (nicht nur Triplett); `cli_errors_01` bis `_06` und `parameter_io_tests` codieren `${PROJECT_ROOT}/release/curcuma` fest (grün mit umgebogenem Pfad, ohne `release/` rot);
  Warnung `-Wdelete-non-virtual-dtor` auf `delete md` in `src/main.cpp` (`SimpleMD`, virtuelle Funktionen, Destruktor nicht virtuell: `simplemd.h:813`); nicht gebaut: MinGW, ROCm (Compiler fehlen).
- **Strukturbibliothek** (`test_cases/structures/`, Phasen 0 bis 4 erledigt): 69 Strukturen, alle Tests und Skripte lesen sie; 65 % der Herkunftsangaben sind `unknown`
  (Altbestand, nur abgeleitet). Offen: Phase 4 (getrackte `.out`-Dateien als Beleg,
  115 untracked Laufausgaben in `test_cases/` entfernen, Berichte in `test_cases/cli/` ins Archiv), case.json je Test. Das Optimierungslevel der `unknown`-Strukturen
  kann nur der Betreiber nachtragen (`python3 scripts/structlib.py report --list-unknown`).
- Wissenschaftliche Validierung der CLI-Tests ausbauen (RMSD-Toleranzen, Energiekonvergenz); Muster für absichtlich fehlschlagende
  Tests (`03_invalid_method` in `curcumaopt`, `rmsd`, `confscan`); Performance-Benchmarks für Regressionserkennung.
- **Testfall-Audit (2026-10-09)**: Drei CLI-Tests bleiben deaktiviert: `cli_curcumaopt_02` (ruft trotz Namens `-method uff` auf, doppelt `01`), `cli_sqm_10`
  (erwartet Log-Marker, die `-opt` nicht schreibt; Konvergenz an `input.opt.xyz` prüfen), `cli_sqm_11` (braucht TBLite). `cli_sqm_04` bis `07` sind wieder aktiv (Toleranz 1e-5,
  `OH` als Radikal ohne `-spin` mit 2e-3, Abweichung 8,2e-4 Eh zur geschlossenschaligen TBLite-Referenz). `simplemd` (16 von 17 Tests) und `curcumaopt` (7 von 7) prüfen nur Exit-Code,
  Dateien oder Zeichenketten, keine Zahlen (Energieerhaltung, Temperatur, Drift).
- **CG-Testbeispiel und `mc_cg_chain`** (2026-10-09): `test_cases/cg/` (`cg_params.json`, `simple_beads.vtf`; die `.vtf` ist als `cg-spheres2` in der Strukturbibliothek, die Parameterdatei nicht) ist untracked und soll als ctest eingebunden werden. `test_cases/mc_cg_chain/` (README, Casino- und Kraftfeld-Konfiguration, 100-Perlen-Kette als VTF, Laufskript; ebenfalls untracked) muss aufgearbeitet werden: Struktur in die Bibliothek, Lauf als Test.
- **Trisaccharid als optimierte Struktur** (2026-10-09): `trisaccharide-c18h32o16` (xtb 6.6.0, Methode unvermerkt) soll durch `.gfn2-opt` ersetzt werden (und `.gfnff-opt` für GFN-FF). Blockiert: TBLite-Referenzen `sqm_reference/reference_data/triose_gfn{1,2}.ref.json` und die GFN-FF-Terme `reference_data/triose.ref.json` gelten für die alte Geometrie und müssen dort neu erzeugt werden (`USE_TBLITE` im `release/`-Build aus); danach `triose=` in `test_cases/CMakeLists.txt` (Zeilen 18, 27) umstellen. GFN-FF-Minimum liegt 1,0 Å (RMSD) von der xtb-Geometrie, Ursache offen.
- ConfScan: Accept/Reject-Meldungen bei Standard-Verbosity nicht sichtbar (Eintrag von 2025, nicht erneut geprüft).

- **Fehlschlagende ctest-Einträge nach dem Merge von origin/master** (14. Okt. 2026, 11 von 353, Ursachen und Belege in [docs/BROKEN_TESTS_2026-10-14.md](docs/BROKEN_TESTS_2026-10-14.md)): confscan 01-07 und der Thread-Override des Loggers (`set_verbosity` ist neben `set_thread_verbosity` wirkungslos, ca. 600 Aufrufe); `cli_analysis_02` (Info-Zeile, Test ohne `-verbosity 2`); `cli_simplemd_14` (keine REACT-Zeile, Ursache offen); `cli_sqm_12` (behoben 2026-10-26: Test umgeschrieben); `xtb_cpscf` (Toleranzen).

## 4. Betreiber-Prüfung offen (🤖/⚙️, nur Sie vergeben ✅)

- CODATA-2018-Einheiten (`191cebe3`; Rest 0,0036 kcal/mol gegen pprcht auf polymer_2x nicht verfolgt), Known Issue #35.
- `eeq_refactor_eps_bohr` Voreinstellung 0 (Faktor-Cache aus) und `nonbonded_skin_bohr` 2,0; Ladungsplatzierung `ensemble` (#31);
  Stillstandserkennung `-opt.stall_steps` 20; `-scf_guess fragments`; `-scf_precond frag` (Okt 2026, opt-in).

## 5. Funktionen und Refactoring (niedrige Priorität)

- **CG Phase 6**: winkelabhängige Ellipsoid-Energie (`calculateEffectiveDistance`, `calculateCGPairEnergy`, Rotationszüge in Casino). Phase 1 bis 5 sind umgesetzt.
- SimpleMD: Physik des Wandpotentials prüfen; RMSD-Strategy-Pattern Phase 3; erweiterte ConfSearch-Algorithmen; bessere Trajektorienanalyse;
  ConfSearch mit GPU und mehreren Threads (nur noch relevant, solange der GPU-Geräte-Pool inaktiv ist, Korrektur 2026-09-28).
- **Trajektorienanalyse zusammenführen** (optional): Phasen 1 bis 3 sind laut `docs/ANALYSIS_CONSOLIDATION_PLAN.md` fertig; offen sind das Umstellen von `trajectoryanalysis.cpp` und `rmsdtraj.cpp` auf die gemeinsame Ausgabe/Statistik und das Aufräumen der Geometriebefehle. Stand des Plans: Dezember 2025, nicht geprüft.
- **Molecule-Refactoring** (Phasen 3 bis 6; die zehn `setXYZComment_N`-Funktionen der Phase 2 gibt es im Code nicht mehr, `src/core/xyz_comment_parser.*` existiert; Stand der übrigen Phasen nicht geprüft): granulare Caches, O(1)-Fragmentzugriff, `ElementType`-Enum,
  SOA/AOS. Plan: `src/core/REFACTORING_ROADMAP.md`, Formate: `src/core/XYZ_COMMENT_FORMATS.md`.
- **Native QM** (Stand November 2025, nicht geprüft): GFN2-Parametererweiterung, PM3-Elementumfang (F u. a.), Validierungsmoleküle, Doku
  "Wann GFN2, GFN1 oder PM3". Der Elementumfang von PM3 im Code ist nicht verifiziert.
- **Build-System** (Stand November 2025, nicht geprüft): bedingte Kompilierung, fünf Build-Varianten, nur der Standardbuild war damals fehlerfrei.
  Seit 2026-09-30 baut `make` in `release/` wieder mit Exit 0; die übrigen Varianten vor einer Wiederaufnahme neu messen.

## 6. Entscheidungen und Fragen

Code-gebundene Entscheidungen, die bei Ihnen liegen. Fachliche Fragen und Projektstatus stehen im Vault und sind nur verlinkt.

- Branch `feature/nebmd-path-sampling` nach master mergen? (enthält überwiegend Master-Arbeit; `/pre-merge-check` fällig.) NEB-MD-Status: Vault `Projekte/curcuma NEB-MD Pfad-Sampling`.
- `-scf_pseudo_diag` als Standard einschalten? (Messungen sprechen dafür, Abschnitt 2.)
- `-scf_precond_guard`: Vorgabe `g2:4:best` kostet auf dem unrelaxierten 599-Atom-Cluster bei 8 Threads 4 Iterationen gegen `none` (Kriterium war <= 3);
  über 5 Threadzahlen liegt frag (55-62) unter none (58-138). Varianten streuen 58-109, eine konvergiert nicht. Vorgabe behalten? [docs/SCF_PRECOND.md](docs/SCF_PRECOND.md)
- `shouldUpdateHBXB()` korrigieren und damit von der Fortran-Referenz abweichen, oder als Referenzeigenschaft belassen?
- UFF-/QMDFF-Gradient reparieren (Stand: "erstmal nicht")? `-seed` für die MD-Startgeschwindigkeiten reparieren (ändert jede MD)?
- GPU-Bitgleichheit angehen (Sammel-Kernel statt `atomicAdd`) oder bei Ensemblevergleichen bleiben?
- Vault, fachliche Fragen: `Offene Fragen/` ("Valenzanteil im reaktiven GFN-FF", "Veraltete C6-Dispersion und Paarlisten in inaktiven Topologie-Ecken",
  "Lücke zwischen GFN2-xTB und Doppelhybrid-Methoden schließen", "RRHO-Näherung durch MD ... ersetzen", "Ladungsplatzierung bei geladenen Mehrfragment-Systemen in GFN-FF - was ist physikalisch besser", Known Issue #31).
- Vault, Projektstatus: `Projekte/curcuma rev-gfnff`, `curcuma GPU-Backends`, `curcuma MD-Stabilität großer Systeme`, `curcuma ConfSearch Weiterentwicklung`,
  `curcuma GFN-FF Transitionsmetalle (MOR41)`, `curcuma Reaktives GFN-FF (react-Topologie)`, `curcuma Konformerenfilter (reorder_paper)`.

## 7. Löschkandidaten

Am 2026-10-01 auf dem Branch `fix/remove-dead-code` entfernt: curcumaopt.cpp, `native_lbfgs_optimizer.*`, `optimisation/modern_optimizer_simple.*`, munkress_2.h,
vier ungenutzte `optimiser/*.h`, core/accuracy_profile.cpp (inkludierte einen nicht vorhandenen Header), `RMSDDriver::MolAlignLib()` und `Prepare*Template()`,
`D4Evaluator::launchGpuKernel()`, helpers/cli_test.cpp und gfnff_test.cpp, zwei untracked `.backup`-Dateien. Volle Testsuite danach: dieselben drei bekannten Ausfälle wie vorher.

Offen, Aufrufer vor dem Löschen erneut suchen:
- `src/capabilities/curcumaopt.h`: hält den PARAM-Block des Moduls `opt` (28 Parameter); die Klasse `CurcumaOpt` ist ohne Implementierung. Die Parameter müssen erst in ein anderes Header umziehen.
- `src/capabilities/rmsd/rmsd_costmatrix.*`, `rmsd_assignment.*` (`CostMatrixCalculator`, `MunkresAssignmentSolver`): kompiliert und von `rmsd_strategies.h` eingebunden, die Klassen selbst werden nirgends aufgerufen. Entscheidung nötig, ob die Strategien sie künftig nutzen sollen.
- dc6dcn-Skip-Pruefung (`FFWorkspaceGPU::shouldSkipDc6dcn()` in `cuda/ff_workspace_gpu.h`, hat keinen Aufrufer; ROCm fuehrt sie mit `CURCUMA_ROCM_LEAN_STEP=0` noch aus).
- `main.cpp` führt `modern_optimizer` noch in der Liste der Scope-Module; das Modul existiert nicht mehr.
- `src/helpers/`: `imagewrite.cpp`, `storage_bench.cpp`, `polymer_topo.cpp`, `gfnff_term_validator.cpp` ohne CMake-Target (eigenständige `main()`-Programme).
- PARAM `gfnff.accuracy` (loose|normal|medium|high) bleibt: kein Leser bis jetzt, soll mit den Genauigkeitsprofilen (Eintrag oben) eine Bedeutung bekommen. Nicht entfernen; der Hilfetext sagt bis dahin "no effect".

## 8. Aus der alten Liste entfernt (überholt oder erledigt)

| Eintrag | Grund |
|---|---|
| "cgfnff Parameter Generation Bug", "Missing Real GFN-FF Parameters" | GFN-FF ist vollständig, die Dateien und das Problem existieren nicht mehr |
| "`gpu_strict` fehlt" | umgesetzt, `src/core/gpu_fallback.h`, [docs/GPU_TUNING.md](docs/GPU_TUNING.md) Abschnitt 1 |
| CG Phase 1 bis 5 | umgesetzt |
| EEQ-Warmstart über die q-Loop-Durchgänge | umgesetzt und nachgemessen 2026-09-18 |
| `make` in `release/` scheitert an CUDA-Unittests | behoben 2026-09-30 |
| Dispersionspaarliste, Repulsionspaarliste | behoben, Known Issue #33 |
| Optimierer: Liniensuchfehlschlag als Konvergenz, Abbruch ohne Struktur, SCF-Unsinn als Ergebnis | behoben 2026-09-27 |
| Hessian: SCF-Schwelle erreicht die Worker nicht | behoben 2026-09-29 |
| Cholesky-Faktor-Cache in der MD nicht exakt | Cache standardmäßig aus (`eeq_refactor_eps_bohr` 0), Betreiber-Prüfung offen (Abschnitt 4) |
- **Generierte Struktur-Tests ausbauen** (2026-10-25, [docs/STRUCTURE_TESTS.md](docs/STRUCTURE_TESTS.md)): (a) Ladung/Spin fuer 57 von 80 Strukturen im Manifest bestaetigen (Vorschlagsliste `assumed_neutral_closed_shell` in `struct_test_policy.json`); (b) Groessenkonsistenz bei Clustern, Strukturen ueber 12 Atome, geladene/offene Systeme; (c) alte Ausgabetext-Tests der `cli_*`-Gruppe einordnen und ersetzen.
- **GFN-FF Inversionsterm: Knick bei C3v-Ammoniak** (Known Issue #42): Ursache und Referenzverhalten (pprcht, xtb) offen; Test `structinv_ammonia_ideal_c3v_gfnff` traegt eine Ausnahme.

