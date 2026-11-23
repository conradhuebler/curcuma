# Curcuma Development TODO List

**Stand**: 2026-10-01. Die frühere Fassung (1380 Zeilen, mit Messprotokollen und erledigten Untersuchungen) steht
unverändert in [docs/TODO_ARCHIVE_2026-10.md](docs/TODO_ARCHIVE_2026-10.md); die Detailmessungen zu den Einträgen unten
sind dort unter der jeweiligen Überschrift zu finden.

Regeln für diese Datei: ein Eintrag hat höchstens drei Zeilen (Was, warum es offen ist, wo die Details stehen). Erledigtes
gehört in `AIChangelog.md` bzw. das Archiv, Fehlerberichte in [docs/KNOWN_ISSUES_ARCHIVE.md](docs/KNOWN_ISSUES_ARCHIVE.md).
"Geprüft" heißt: am angegebenen Datum gegen den Quelltext gelesen, nicht neu gemessen. Dateien und Zeilen gelten für diesen
Stand.

## 1. Korrektheit (offen)
- **GFN-FF GPU vs CPU mit isolierten Na/Cl-Ionen in Zelle** (2026-11-23, gefunden in S3): D4-Paarterm 3,5e-9 Eh und seine CN-Kettenregel (Virial 1,1e-7 Eh, Gradient 3,5e-8 Eh/A) weichen ab, gleich im A1-Basisbuild, nicht Ewald; Ursache nicht gesucht. Test `gfnff_gpu_pbc_ewald` (a) Na2Cl2-Box hat dort gelockerte CPU-Grenzen.

- **GPU und periodische Simulation mit Druck (Betreiberziel, Nov 22, 2026)**: GFN-FF und GFN1/GFN2 sollen auf jedem Backend vollständige
  NPT-Simulationen können. Fehlend auf CUDA: Druckterme der scharfen Radien (A3); A1 (Virial/Barostat mit `-gpu cuda`) und A2/S3 (Geräte-Ewald und pcg-Operator) erledigt 2026-11-23. ROCm: Teil B.
  Zuerst alle CUDA-Lücken auf dem CUDA-Rechner, dann die ROCm-Lücken auf einem Rechner mit ROCm-GPU: [docs/GPU_PERIODIC_GAPS.md](docs/GPU_PERIODIC_GAPS.md).
- **`-large_system_mode` mit Zelle liefert still eine Vakuumenergie** (`native_xtb_method.cpp`, Rückkehr vor `setPeriodicCell`): harter Fehler nötig,
  Doku widerspricht sich (`WP-PERIODIC-XTB.md` Z. 45/167/1376). Geprüft 2026-11-22. Eintrag A5 in [docs/GPU_PERIODIC_GAPS.md](docs/GPU_PERIODIC_GAPS.md).

- **GFN-FF unter Last nicht bitgleich** (2026-11-14): zwei identische `-sp`-Läufe (Defektstruktur der Ladungsgruppen, 273 Atome, `-threads 1`) unterscheiden sich unter paralleler Last um bis zu 6.7e-16 Eh/A (Gradient bis 0.74); im Leerlauf 6 von 6 identisch, unter Last 4 von 8 verschieden, auch auf `devel/methods` pur. Die Ursache (Summationsreihenfolge trotz einem Thread) ist nicht gesucht.
  Folge: Testskripte mit Bitvergleich (`check_cn_zero_guard.py`, `check_adaptive_step.py`, `check_lbfgs_atom_step_cap.py`, `check_d4_fused_c6.py`) können unter `ctest -j` flackern; `check_charge_groups.py` vergleicht jetzt mit 1e-12.
- **FEHLER: GMTKN55 `G21EA/EA_16` (SiH3-) mit gfn1 liefert keine Energie** (2026-10-30, Ursache bekannt, vorerst gewollt): Abbruch "implausible charge distribution" (max |q| 5.55 > 5.00, auch aus H0); xtb 6.7 konvergiert auf dieselbe SCC-Energie (-7.94141906 Eh) mit Si = -5.55 e, derselbe q^3-Runaway.
  Der Fehler bleibt, die Struktur fehlt damit in jedem gfn1-Vergleich. Keine gültige Energie für dieses System.
- **PLUMED-Bias nur auf `feature/bmt-plumed-nebmd` korrigiert** (`96c29c50`, 2026-10-30): `upstream/master` hat Vorzeichenfehler und kB statt kB*T noch; biased PLUMED-Läufe davor sind ungültig. Offen: nach master bringen, betroffene Läufe neu rechnen. [docs/PLUMED_HELP.md](docs/PLUMED_HELP.md).
- **`-md -method uff` schreibt `<basename>.param.json` trotz BMT ins Arbeitsverzeichnis** (2026-10-30, gemessen, alter und neuer Stand gleich): verstößt gegen die BMT-Regel in `CLAUDE.md`.
- **MD-Instabilitätsprüfung schlägt bei wenigen Freiheitsgraden fälschlich an** (2026-11-10, gemessen n=10 Seeds): `m_unstable = T > 10000 * m_T` (`simplemd.cpp` Z. 4252/4724) vergleicht mit dem Vorschritt; Ar2 NVE (1 DOF) bricht mit Seed 1, 4, 9 an einem Umkehrpunkt (Ekin ~ 0) ab. `plumed_bias` nutzt deshalb Seed 2.
- **PLUMED rechnet Schritt 0 doppelt** (2026-10-30): Init und erster Integrationsschritt rufen beide `performPlumedCalc()`, COLVAR hat zwei Zeilen bei t = 0; Folgen für HILLS nicht geprüft. `energyUnits = 2625.5` steht als Literal statt über `CurcumaUnit`.



- **`-gfnff.charge_groups` auf der GPU** (2026-10-23): GPU-Methode faellt auf die CPU zurueck (Aenderung in `gfnff_gpu_method_impl.h` nicht uebersetzt, nicht gelaufen); mehrere Fragmente nur, wenn jede Gruppe eine Vereinigung ganzer Fragmente ist (2026-10-24) oder mit `"rest": "per_fragment"` (2026-10-25).

- **`-gfnff.charge_groups auto`: offene Grenzen** (2026-10-23): nur Kationen (Anionen, Carbenium fehlen); nur ein Fragment (Ionenpaare, S30L-CI 24/27/28 nicht erfasst); Gruppe einmal je Instanz, kein Neu-Gruppieren bei Bindungsbruch (unstetig, nicht gelaufen); GPU nicht gelaufen. Gruppenausdehnung (Radius 2, 14 Atome) offen (Schritt 3, Potentialmass). [docs/POLYMERBUILD_CHARGE.md](docs/POLYMERBUILD_CHARGE.md)

- **`-gfnff.charge_groups auto` auf hypervalenten P/S-Geruesten** (2026-10-23, n=3): PArel h4p2o71/h4p2o72/s4o41 werden erkannt (Bedingung (c) blockt P/S nicht), Energie +53 bis +178 kcal/mol; freie Gruppensumme bis -0,90 e wird auf +1 gezwungen. Ob P/S-Partner blocken sollen, ist zu entscheiden.

- **Topologie-Cache `<basename>.topo.json` traegt die Molekuelladung nicht** (2026-10-23, n=1): eine bei Ladung 2 geschriebene Datei wird bei Ladung 0 wiederverwendet (gleiche Energie -14,11678988 Eh, (AB)2). Cache-Schluessel um Ladung/Spin ergaenzen; nicht untersucht, ob die Ladung beim Laden doch gelesen wird.

- **`cli_simplemd_14_gfnff_react_h_recombination` schlaegt fehl** (2026-10-23, n=1 Binary mit und ohne Ladungsgruppen-Aenderung): "no bond formation events", MD laeuft 0,2 s ohne Fehler. Ursache nicht untersucht.

- **`polymerbuild`: Optimierungsschritte laufen mit Ladung 0, ohne Ladungsgruppen** (2026-10-23): die Formalladung steht nur im Ergebnis. `-gfnff.charge_groups` existiert, wird beim Bauen aber nicht benutzt (Gruppen der wachsenden Kette je Schritt noetig). [docs/POLYMERBUILD_CHARGE.md](docs/POLYMERBUILD_CHARGE.md)

- **`-opt` (LBFGS++): Abriss der C=O-Bindung bei 5747 Atomen** (diagnostiziert 2026-10-29, [docs/LBFGS_ATOM_STEP_CAP.md](docs/LBFGS_ATOM_STEP_CAP.md)): Schalter `-lbfgs_max_atom_step` (Default auto: 0,2 A ab 2000 Atomen) behebt es in L1. Offen: Schwelle 2000 Atome belegen (nur eine Struktur gemessen), Cap auf MOR41/GMTKN55 pruefen, Bindungscheck am Ende von `-opt`, Nebenbefunde (`alpha=`-Ausgabe, Hilfe `lbfgs_line_search`, Doppel-Default `lbfgs_max_line_search`).
- **Entscheidung: Default von `-md.seed`** (2026-10-31): seit dem Seed-Fix ist der Default -1 = Uhrzeit, eine MD ist also nicht mehr bitgleich wiederholbar (vorher war sie es, weil der Seed wirkungslos war); `cli_gfnff_07_eeq_localinv` wurde dadurch flaky und nutzt jetzt `-md.seed 1`. Alternative: Default 0 (aus Atomzahl und T0, mit der Startgeometrie im Seed reproduzierbar). Betrifft jede MD und alle Vergleichslaeufe.
- **Entscheidung: doppelte MD-Tabellenzeilen bei dt < 1 fs** (2026-11-09, F-MD Punkt 2): ganzzahliges `print_frequency` testet `int(Schritt*dt) % f`, bei dt 0,5 fs drucken Schritt 100 und 101 die 50-fs-Zeile; nicht behoben, weil `gfnff_d4_fused_c6` eine Zeile je Schritt zaehlt. Fix (erste Zeile je Rasterzeit) braucht angepassten Test. [docs/WP-PERIODIC-NONBONDED.md](docs/WP-PERIODIC-NONBONDED.md) F-MD.
- **ROCm: Ladungsgruppen und `cn_zero_guard` bleiben auf der CPU** (Luecke, 2026-10-31): CUDA hat beides auf dem Device (Commit `fix/gpu-charge-groups`, gemessen 15,5 ms/Schritt bei 5747 Atomen), ROCm nicht: kein Hardwaretest, Trait `has_charge_groups_cn_guard = false` haelt den CPU-Rueckfall, der HIP-Kern ist nie kompiliert (kein hipcc) und nie gelaufen. Zu tun auf ROCm-Hardware: GPU gegen CPU wie bei CUDA (Gruppensummen, Gradient Cl-, 200 MD-Schritte), dann Trait auf true. Vulkan: nicht betrachtet. Ausserdem: `k_subtract_qtmp` (CUDA) hat keinen Aufrufer, Loeschkandidat.
- **Rechenzeit je MD-Schritt ~ N^1,55** (diagnostiziert 2026-10-29, [docs/MD_STEPTIME_SCALING_2026-10.md](docs/MD_STEPTIME_SCALING_2026-10.md)): All-Pairs-Anteil (EEQ-Matrix/PPCG, Coulomb, D4-Liste), N^2 = 74 % bei 5750 Atomen. Offen, Hebel je 5-15 %: Constraint-Block nur fuer Schur/LU, pAp-Pruefung statt serieller Definitheitspruefung, D4-C6-Auffrischung fusionieren; Linearisierung braucht Cutoff/Ewald (aendert Energien, Entscheidung).

- **Literaturverweise der Trajektorien- und Ionen-Methoden pruefen** (Oktober 2026): die Quellen fuer Blockfehler, H-Brueckenkriterium, g(r), Shrake-Rupley, Sondenradius, Bondi-Radien, Fibonacci-Punkte und `-ionize`
  stehen in `docs/TRAJPROPS.md`, `src/core/citation_database.cpp` und `src/capabilities/trajprops.h`; sie stammen aus dem Gedaechtnis der KI und sind NICHT gegen die Arbeiten geprueft. Ebenso offen: welchen Radiensatz Min 2015 fuer die WASA verwendet.

- **Flat-Flag ohne Modulpräfix erreicht den nativen xTB nicht** (bestätigt 2026-10-11, n=1): `-diis_start 1` landet laut `-export_run` in `/opt` und `/sp`, nicht in `/xtb`; der Name ist in
  `curcumaopt.h`, `lbfgs.h` und `native_xtb_method.h` registriert, bei `-verbosity 2` kommt keine Warnung. `-xtb.diis_start 1` wirkt (seit `1bb4781b`). Ursache im Routing nicht untersucht;
  gleiches Muster für andere mehrfach registrierte Namen möglich. Vault: Labor `curcuma GFN2-SCF-Vorkonditionierung`, Schritt 19.
- **Flache SCF-Namen, die auch das Optimiermodul kennt** (2026-10-26, Ursache gefunden): `-diis_start`, `-max_iterations` gehören auch zu `opt`, das `-sp` und `-opt`
  benutzen; das Befehlsmodul gewinnt, also erreichen sie den SCF nicht (keine Warnung). Wirksam: `-scf_max_iterations`, `-scf.diis_start`, `-xtb.diis_start`.
  Offen: soll `-sp` die Optimiererparameter überhaupt besitzen? [docs/PARAMETER_SYSTEM.md](docs/PARAMETER_SYSTEM.md), "Shared module scf".
- **SCF-Modul `scf` bei tblite/xtb/Ulysses nicht gelaufen** (2026-10-26): Abbildung auf `max_iterations`/`damping` nur per `-verbosity 3` gesehen, die
  TBLite-Änderung (GFN2-Anhebung auf 150/0.6 nur ohne Nutzerwert) nur syntaxgeprüft; kein Provider in `release/`. `-hessian` endet mit 0, wenn ein
  Verschiebungs-SCF scheitert (pm3, `-scf_max_iterations 3`), nicht untersucht. [docs/SCF_MODES.md](docs/SCF_MODES.md)
- **Offene Schalen über die CLI: Aliase und ctest fehlen** (2026-10-11): `-multi`, `-multiplicity`, `-uhf`, `-xtb.spin` werden ignoriert (seit 2026-11-09 mit der Warnung "is unknown"); Zeile 2 "Ladung Multiplizität" der ORCA-xyz wird nicht gelesen;
  ungerade Elektronenzahl ohne `-spin` besetzt fraktional (tblite), xtb nimmt ein Dublett (Entscheidung offen); es gibt keinen ctest für `-spin` (Known Issue #9). Flag mit Leerzeichen und Paritätswiderspruch
  sind seit 2026-10-11 Fehler bzw. Warnung (`CLI2Json`, `buildReferenceOccupations`). Vault: Labor `curcuma GFN2-SCF-Vorkonditionierung`, Schritt 13.
- **GFN-FF-Energie hängt von der Atomnummerierung ab** (911 von 2557 Strukturen, bis 237 kcal/mol bei geladenen Mehrfragment-Systemen;
  in pprcht und curcuma). Referenzeigenschaft, kein Portierungsfehler. [docs/REV_GFNFF_TODO.md](docs/REV_GFNFF_TODO.md) #12, Known Issue #34.
- **GPU-Schrittzeit GFN-FF, Reste und nicht ROCm-spezifische Hebel** (Stand 13. Okt. 2026, Messwerte in [docs/ROCM_PERFORMANCE_ANALYSIS_2026-10-13.md](docs/ROCM_PERFORMANCE_ANALYSIS_2026-10-13.md)): (a) MD-Integrator und Thermostat auf der GPU statt `SimpleMD` auf der CPU (7320 Atome: 4-7 ms Luecke und 4,3 MB Transfer pro Schritt, etwa 10 %); (b) viele kleine Moleküle gebündelt auf der GPU (Konformerensuche, Durchsatz); (c) Coulomb (10 ms bei 7320 Atomen) und Dispersion (7,5 ms) weiter ausreizen, auch 20 % lohnen sich; (e) `-gpu auto` samt Schwelle entfernt (2026-10-26); Break-even Wasser auf ROCm zwischen 300 und 600 Atomen, andere Systemtypen offen; CUDA gemessen 17. Okt. 2026 ([GPU_TUNING.md](docs/GPU_TUNING.md)): MD-Schritt ab etwa 66 Atomen schneller, aber 0,8 s Start je Prozess, `-gpu auto` hat auf CUDA keine Schwelle (Entscheidung offen); (f) projiziertes PCG macht bei 600 Atomen ca. 11 Host-Roundtrips je Schritt, Abhilfe: Konvergenzflag auf dem Device; (g) Logzeile sagt bei kleinem Device-EEQ noch "WP2 dense Cholesky + CPU Schur" (kosmetisch); (d) danach GFN1/GFN2 auf ROCm (Host-Eigenloeser unterhalb eines Break-evens, Cache-L, Pseudo-Diagonalisierung), erst nachdem GFN-FF ausgereizt ist.
- **Periodisches natives GFN1/GFN2 (WP-PERIODIC-XTB, Stand P5 + T + P4, 2026-11-15):** `-method gfn1|gfn2` mit orthorhombischer Zelle: periodische Gamma-Punkt-Energie, Gradient (FD-Gate 8e-9 Eh/A) und Virial (Zellskalierungs-FD 1e-9 Eh), SimpleMD NVT/NVE/NPT (`-compute_pressure`, `-barostat`; nur Rauchtest, keine Dichte-Kampagne); gfn1 3-7e-11 Eh zu tblite, gfn2 bis 1e-3 Eh (8 A) von tblite entfernt, weil tblites `get_amat_sdq_rec_3d` alpha-abhaengig ist (Faktor 2 Ladung-Dipol, +1 statt -1/3 Ladung-Quadrupol; nicht gegen ein korrigiertes tblite geprueft, Entscheidung des Betreibers offen); CUDA (G1, reduziert): `-gpu cuda` mit Zelle nutzt die auf dem Host gebauten periodischen Matrizen und den geraeteresidenten SCF (Eigenloeser), Potential/Dispersion/Gradient/Virial bleiben auf dem Host (Paritaets-ctest `pbc_xtb_gpu_parity`); G3/G4a/G4b/G5a (Multipolmatrizen, ATM, D4 dE/dq, Multipolgradient+Virial) auf dem Geraet, 432 H2O MD-Schritt 6,9 -> 3,7 s; G2 (Ewald-gamma), G5b-d (W, CN-Kette, Repulsion) auf dem Geraet, 432 H2O 3,13 s; G1b (S/H0-, Multipol-Integral-Bildsummen, Cholesky) auf dem Geraet, 432 H2O 2,40 s je MD-Schritt; G1c Pseudo-Diagonalisierung auch in der hostgefuehrten periodischen Schleife, 2,20 s; E1/A1/M1/M2 (Komplement-Dichte, ATM je Tripelklasse, Multipol-k-Summen von Gradient/Virial ueber Strukturfaktoren) 1,74 s; C1 (D4-CN-Faltung) 1,65 s; H1 (H0/Pulay-Gradient) 1,53 s; H2/H3 (Ewald/KO-Gradient, CN) 1,39 s; H4/H5 (periodisches D3, D4-Paare) 1,35 s; D1 (D4 dE/dq) 1,34 s; XB mit Zelle unerreichbar (Br/I haben d-Schalen); offen: hostgefuehrte periodische SCF-Iteration (26 % des Schritts, braucht die residente Schleife mit Zelle), ROCm/Vulkan nicht geplant (CPU-Rueckfall), Strukturfaktor-Faktorisierung der Multipol-k-Summe der Energiematrizen (G3, O(nat^2 N_k)), Schub-Dehnungen im Virial-Gate (nicht orthorhombisch); Stufe T (2026-11-15): periodische D4-Paare/ATM/dE/dq und Multipolgradient blockweise parallel und thread-deterministisch, 216 H2O gfn2 10.7 -> 2.5 s (8 Threads); nicht getestet: XB mit Zelle (keine Br/I-Struktur), geladene Zellen, offene Schale, GPU-Rueckfall auf echter GPU, Energieerhaltung in der MD; d-Schalen mit Zelle abgelehnt. [docs/WP-PERIODIC-XTB.md](docs/WP-PERIODIC-XTB.md) Befund: dftd4 `get_atm_dispersion_derivs` ohne `triple_scale` im Gradienten (Selbstbild-Tripel doppelt), hier korrigiert, gegen dftd4 selbst nicht geprueft. [docs/WP-PERIODIC-XTB.md](docs/WP-PERIODIC-XTB.md)

- **Periodisches GFN-FF fehlt** (zurueckgestellt, Entscheidung Okt. 2026; ueberholt durch `feature/bulk`, siehe naechster Eintrag): die CPU nutzte das Gitter nirgends, die ROCm-Kernel haben `applyMIC` und `setUnitCell`, aber `unitcell 9 9 9` aendert die Energie von `water8-cluster` weder auf CPU noch auf ROCm (8 Stellen). Plan und Stand: [docs/GFNFF-CPU.md](docs/GFNFF-CPU.md), [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md) ("Periodic boundary conditions: not implemented").

- **Periodisches GFN-FF, Stand `feature/bulk` (2026-11-01):** Zelle, CN, D4-Paarterm und -Schwanz, Abstossung, Bindungsterme (A4a), HB/XB (A4b), Topologie-Erzeugung, reaktive Topologie und Ladungs-Ensemble-Fenster (A4c; zerrissene Molekuele werden fuer den Aufbau zusammengesetzt, ein an sein eigenes Bild gebundenes Molekuel wird abgelehnt) sind CPU-periodisch; Coulomb-Paarterm per Ewald hinter `-gfnff.coulomb_ewald` (B2, Standard aus, weil die EEQ-Ladungen noch nicht periodisch sind) ([docs/WP-PERIODIC-NONBONDED.md](docs/WP-PERIODIC-NONBONDED.md)); offen: EEQ-Matrix und -Gradient mit Ewald (B3; dann `coulomb_ewald` Standard true bei Zelle), `topology_mode auto` mit Zelle, Pass-2-Bindung ueber die Zelle zwischen getrennt platzierten Fragmenten (nur Warnung), GPU-Portierung (S0-S2 auf CUDA hinter `gpu_periodic`, s. u.; sonst CPU-Fallback bei Zelle), Barostat/NPT (N2; der atomare Virial samt Ewald- und D4-Schwanz-Anteil ist seit N1 da, s. naechster Eintrag).

- **GPU-residente MD (Vorbereitung, nicht begonnen)** (2026-11-05): Integrator/Thermostat/Wrap auf dem Geraet statt in `SimpleMD`. Nach M1 (Fragmentsuche behoben) kostet der Host-Integrator nur noch 1-3 ms/Schritt und die Transfers < 1 ms; der Energieaufruf (Host-EEQ im Hybrid) dominiert. Lohnt daher nur als Teil von S3/S4 (Geraete-EEQ/Ewald). Plan und offene Fragen: Vault `Projekte/curcuma GPU-residente MD (Vorbereitung)`.

- **Barostat / NPT (N2, 2026-11-05) und Virial (N1):** `-barostat berendsen|crescale` mit `GFNFF::updateCell` (CPU, isotrop) ist da (ctests `pbc_npt_*`, [docs/NPT.md](docs/NPT.md)); offen: (a) C-rescale-Gleichungen gegen den Originaltext pruefen (aus Auftrag/Gedaechtnis, eine Websuche ohne Treffer; Idealgas-Test stimmt mit N_eff = dof/3 + 1), (b) Virial der RATTLE-Zwangskraefte (mit `-rattle` abgelehnt; NPT-Produktion laeuft daher ohne Constraints), (c) Virial/Zellupdate im CUDA-Hybrid seit A1 (2026-11-23), keiner fuer `frag_charge_model ensemble`, ALPB, ATM, (d) molekularer Virial (COM-Korrektur) nicht berechnet, (e) Plugin-ABI 4: CUDA-Plugin muss neu gebaut werden (hier nicht gebaut), (f) Zelle in der Trajektorien-Kommentarzeile nur mit Barostat (ohne Barostat bleibt die Trajektorie md5-gleich zu N1; Entscheidung des Betreibers, ob jede periodische MD sie schreiben soll), (g) Paarlisten werden bei jeder Zellaenderung neu gebaut (kein Verlet-Argument fuer ein geaendertes Gitter; Mehraufwand gemessen in N2), (h) anisotrope/semi-isotrope Kopplung nicht implementiert. [docs/WP-PERIODIC-NONBONDED.md](docs/WP-PERIODIC-NONBONDED.md) N1/N2.
- **Druck-Fernkorrektur (V-LRC, 2026-11-09):** `pressure_lrc` (Standard true) korrigiert nur den D4-Schwanz (g(r) = 1, Grossabstandsform der Daempfung; bei 15 Bohr Rest -3,5 bar an Wasser); fuer die anderen scharfen Radien fehlt der Kreuzungsterm im Druck: gemessen (A3, 2026-11-23) HB -31 bis -38 bar, Abstossung +1,3, CN -0,45 bar; D4-Kreuzung = V-LRC auf 3 %; HB-Radius x2 konvergiert bei -51 bar (Modellaenderung); Entscheidung (a) belassen, (b) statistischer HB-Kreuzungsterm, (c) groesserer Radius als Schalter offen ([docs/WP-PERIODIC-NONBONDED.md](docs/WP-PERIODIC-NONBONDED.md) A3-P). Laeufe vor V-LRC (W1-Matrix, P1, P7, P8, P2, P10) mit altem Druck, NPT-Dichten 0,04-0,32 % zu niedrig ([docs/NPT.md](docs/NPT.md)); CUDA-Plugin auf ABI 6 neu bauen.
- **Periodisches GFN-FF auf CUDA, Stand S2** (2026-11-04, Opt-in `-gfnff.gpu_periodic true`): S1 lokale Terme, Paarbauer, HB-Erkennung, D4-Schwanz mit MIC auf dem Geraet (ctest `gfnff_gpu_pbc_terms`); S2 Hybrid mit `coulomb_ewald` (Standard): periodische EEQ und Ewald-Coulomb auf dem Host, `hasPeriodicEnergy()` true (ctest `gfnff_gpu_pbc_hybrid`), Gewinn gegen CPU nur 1,35-1,46x (MD-Schritt 114 -> 78 ms bei 2910, 324 -> 241 ms bei 5748 Atomen), weil die Host-EEQ (42 / 107 ms je Aufruf) den Schritt traegt (der SimpleMD-Aufwand ausserhalb des Energieaufrufs, 31 / 120 ms, war die Fragmentsuche des Zellwraps, seit M1 1,2 / 2,3 ms; Hybrid jetzt 54 / 135 ms je Schritt); S3 (2026-11-23): Ewald-Summe und pcg-Operator auf dem Geraet, NPT-Schritt 1,8-1,9x schneller als A1 (WP S3); offen: Verlet-Trigger der CN-Paarliste (`cnPairListStale`) nur ueber den MIC-Kernel, nicht einzeln getestet, ROCm-Spiegel (S6, `has_periodic_local_terms` false), BATM/ATM ohne MIC (wie CPU), Mehr-GPU-Split mit Zelle gesperrt. [docs/WP-PERIODIC-NONBONDED.md](docs/WP-PERIODIC-NONBONDED.md) S1/S2.
- **`hb_scale`, `repulsion_scale`, `dispersion_scale`, `coulomb_scale` auf der GPU** (2026-11-08/09, S-HB, T-SCALE): Rueckfall auf die CPU bei `s != 1` in `createNativeGfnff` ist nicht laufzeitgeprueft (kein Plugin mit ABI 7 gebaut); nach dem naechsten CUDA-Build einmal `-gpu cuda -gfnff.hb_scale 1.5 -gpu_strict true` (muss mit Exit 3 stoppen), ebenso je ein Term-Schalter. Die GPU-Kerne kennen die Faktoren nicht; `evaluateEwaldCoulomb` (GPU-Hybrid) wendet `coulomb_scale` nicht an.
- **Befunde aus S1** (2026-11-04): (a) `-gfnff.hbond false` ist CPU/GPU-asymmetrisch: die GPU behaelt die HB-CN-Modulation der X-H-Bindungen (urea197 6,5e-3 Eh/A); (b) `-gfnff.coulomb false` schaltet auf der CPU seit S1 auch Selbstterm und TERM 1b ab (vorher nur den Paarterm; GPU-Paritaet) - Entscheidung des Betreibers, ob das so bleibt (kein ctest/Skript benutzt den Schalter; die D1-Vergleichslaeufe `nve_coulomboff` in docs/WP-PERIODIC-NONBONDED.md wurden mit dem alten Verhalten gemacht).
- **GPU-Plugin-ABI, offene Reste** (2026-11-03, S0): Version gilt fuer `libcurcuma_{cuda,rocm,vulkan}.so`, nicht fuer `libcurcuma_cuda_mgpu.so`; ROCm-/Vulkan-Einstiegsdatei mit Versionssymbol ungebaut (kein SDK hier); bei jeder Aenderung an `ComputationalMethod`/Plugin-Schnittstellen `CURCUMA_GPU_PLUGIN_ABI_VERSION` erhoehen (keine automatische Pruefung, z. B. CI-Vergleich der vtable). [docs/GPU_PLUGIN_STARTUP.md](docs/GPU_PLUGIN_STARTUP.md)

- **RATTLE, Befunde M1** (2026-11-05, [docs/WP-PERIODIC-NONBONDED.md](docs/WP-PERIODIC-NONBONDED.md) M1): (a) die 1-3-Liste in `InitConstrainedBonds()` findet nur Ketten i > j > k (`for (int k = 0; k < j; ++k)`), H-O-H mit O zuerst bekommt 0 Winkel-Constraints, `-rattle_13 true` wirkt auf die Wasserboxen nicht (Fix aendert jeden `rattle_13`-Lauf); (b) Bindungen, die in der Eingabe ueber eine Zellflaeche reichen, bleiben unconstrained (nur Warnung; Abhilfe: Eingabe auspacken oder Liste/Schleifen mit MIC); (c) Standard `rattle_tol_13` 1e-3 A^2 kostet starres Wasser 100x Energieerhaltung gegen 1e-7 (n = 1); (d) `DistanceMatrix()` beim Aufbau: 2 dichte NxN-Matrizen (2 x 264 MB bei 5748 Atomen).
- **Fragmentsuche je MD-Schritt** (2026-11-05, M1): Rest 1,2 / 2,3 ms (2910/5748 Atome) fuer die Zellgitter-Suche in jedem Schritt; Wiederverwendung mit Haut-Kriterium spart sie, macht den Wrap aber verlaufsabhaengig (Entscheidung Betreiber).

- **`wall_potential pbc` ohne periodische Energie** (2026-11-03, A6): Notloesung (Freiraum 1,8 A, Reflexion) bleibt fuer alles ausser CPU-GFN-FF mit Zelle und `coulomb_ewald`; offen: Box zur Zelle machen (Option 1), Entscheidung Betreiber. Option 2 (+L/-L) erledigt in A6b. [docs/WP-PERIODIC-NONBONDED.md](docs/WP-PERIODIC-NONBONDED.md) A6.

- **MD-Restart: Datei schlaegt CLI** (2026-11-02, Befund D1): `-restart_file` uebernimmt `thermostat`, `dT`, `MaxTime`, `T` aus der Datei und ueberschreibt die CLI-Flags; ein Protokoll "CSVR, dann NVE aus dem Restart" braucht deshalb eine editierte Datei. Entscheidung offen, ob die CLI gewinnen soll (docs/WP-PERIODIC-NONBONDED.md D1). Seit F-MD (2026-11-09) gewinnt `-seed` (nur dieser Schluessel).

- **Ewald-Coulomb/EEQ (B2/B3/P1) offen** (2026-11-03): nach P1 2910 Atome 70 ms je `Calculation(true)` (ohne Zelle 46), groesster Posten Kernel-Aufbau je Geometrie (erfc, 56 ms bei 5748 Atomen) und SimpleMD mit Zelle (~37 ms je Schritt mehr als `Calculation`, nicht lokalisiert); GPU mit Zelle -> CPU-Fallback ohne `gpu_periodic`, mit Opt-in der S2-Hybrid (EEQ/Ewald weiter auf dem Host); Virial/NPT fehlt; `GAUSSIAN_CORRECTION_XMAX` 6 erzwingt r_g < L/2 (Cl/Br/I 17-18 Bohr, Ne 25); dxi liest bei 2-fach koordiniertem C einen Rohgeometrie-Winkel. [docs/WP-PERIODIC-NONBONDED.md](docs/WP-PERIODIC-NONBONDED.md) P1.
- **Solvername wirkt auf den Topologieaufbau** (2026-11-03, Befund P1, auch vor P1): `-solve_method pcg` statt `auto` aendert die periodische Energie der 2910-Atom-Box um 3,75e-6 Eh, cholesky vs pcg bei 900 Atomen 5,4e-7 Eh (Winkel-, Torsions-, HB-Parameter); vermutlich ein nichtperiodischer EEQ-Loesungsschritt im Topologieaufbau, der den Namen liest. Nicht untersucht.

- **XB-Gradient falsch (GFN-FF, CPU und CUDA, unabhaengig von der Zelle)** (2026-11-01): `calcHalogenBonds` Daempfungsableitungen: Kurzreichweite mit falschem Vorzeichen, Langreichweite r/c statt 1/r (`k_*` in `cuda/gfnff_kernels.cu:1952-1958` gleich); FD-Fehler 1,2e-4 von 5,4e-4 Eh/A (F-I...Pyridin + Essigsaeuredimer); beide Korrekturen zusammen: 4e-12 (ausprobiert, nicht uebernommen). [docs/WP-PERIODIC-NONBONDED.md](docs/WP-PERIODIC-NONBONDED.md) A4b.

- **Periodisch, kleinere Punkte** (2026-11-01): D3 (gfn1 nativ, `uff-d3`) ohne Verbraucher, da dort nirgends MIC; natives gfn1/gfn2 liest `m_has_pbc` nirgends (rechnet mit Zelle still nicht-periodisch); D4-Schwanz ohne R0-Daempfung (Fehler 2,1 % bei r_c 11,5 Bohr, 3e-5 bei 35); `calculate()` ausserhalb `calculationSingle()` behaelt den letzten Schwanz (react `dE_jump`, `NumGradFixedCharges`); NVE-Erhaltung mit Schwanz ungeprueft.

- **GFN-FF `cn_accuracy`: Abweichung von der Referenz** (2026-11-01): Standard 1 gibt 10 Bohr, die Referenz 0,1 = 12,25 Bohr (`gfnff_interface.f90:256`); Topologie-CN fest 40 Bohr statt `cnthr`; Hilfetext "10 matches the reference" und Kommentar bei `mchar` (accuracy 1) stimmen nicht. Keine messbare Energiewirkung (erf-Terme jenseits 10 Bohr exakt 0).
- **Frag-Ladungs-Ensemble auf der GPU**: geladene Mehrfragment-Systeme laufen seit 14. Okt. 2026 auf der CPU-Methode (GPU-Rueckfall, gezaehlt); das Ensemble selbst (`GFNFF::fragEnsembleBlend`) laeuft nur dort. Auf der GPU abbilden waere die Alternative, [docs/FRAG_CHARGE_MODEL.md](docs/FRAG_CHARGE_MODEL.md).
- **ROCm**: GFN1/GFN2 ignorieren `-scf_mixed_precision false` (G2-13); Energiesummen nehmen wave32 an (F-17). Dokumentiert in
  [docs/MULTI_GPU_GAPS.md](docs/MULTI_GPU_GAPS.md). Leistung (EEQ auf der CPU, rocSOLVER-Eigenloeser):
  [docs/ROCM_PERFORMANCE_ANALYSIS_2026-10-13.md](docs/ROCM_PERFORMANCE_ANALYSIS_2026-10-13.md).
- **Exitcode bei fehlgeschlagenem Setup**: `-md` gibt seit `a2b7c5f0` 1 zurück; prüfen, ob dasselbe für `-opt` und `-sp` gilt (nicht geprüft).
- **UFF- und QMDFF-Gradient falsch** (Eintrag vom 2026-09-27: UFF-Winkel ohne den Faktor -sin(theta); Fix vorbereitet, Betreiber: "erstmal nicht").
  Nicht erneut geprüft; Details im Archiv unter "UFF- und QMDFF-Gradient falsch".
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
- **Segfault: `-md_method gfnff -opt_method uff` in ConfSearch (Sep 11, 2026, not re-checked after the Oct 2026 master merge)**. Crash in
  `~EnergyCalculator` while `CxxThreadPool::~CxxThreadPool()` tears down the MD threads, after the
  run has produced correct results. Reproducible 3/3. Measured scope: ONLY this pair --
  gfnff/gfn2, gfnff/qmdff, gfnff/eht, uff/uff, uff/gfnff and gfn2/uff all run clean, so no
  production path is affected. Appeared with the Sep 2026 master merge (`611cdd6a`): pure master
  does not crash on the same command, and the confsearch branch did not before the merge, so it
  is the combination of this branch's ConfSearch with master's ForceField/FFWorkspace rework.
  Not triggered by any ConfSearch switch (`-topology_lock`, `-confgen_phase`, `-opt_feedback_bias`,
  `-relax_pes`, `-snapshot_dedup_rmsd` all still crash). `MALLOC_CHECK_=3` reports no heap
  corruption, so it is a stale/invalid object rather than an overflow -- next step is an
  AddressSanitizer build. `cli_confsearch_04` was moved to qmdff, so the suite is green and this
  issue is NOT covered by a test any more.
- **LBFGSpp prüft die Krümmung nicht** (geprüft 2026-10-18, n=1): `BFGSMat::add_correction` überspringt das Update bei `s'y <= 0` nicht; auf einer rauen Delta-ML-Fläche (alte Defaults, 120 Zentren) bricht der Optimierer `auto` mit Exit 1 ab. Beobachtet an `cli_deltaml_01` mit dem alten ACSF-Gitter, mit den neuen Defaults (6 Schalen, eta 4) konvergiert er. Details: `AIChangelog.md` (Delta-ML), Vault `Projekte/curcuma Delta-ML`.
- **`-opt.optimizer lbfgs` macht bei großem Startgradienten einen Nullschritt** (n=160 Läufe, 2026-10-19): ab etwa 0,3 Eh/Bohr Exit 1, 0 Iterationen, kein Fehlertext, auch für GFN2 und GFN-FF allein (400-K-MD-Frames); der Standardoptimierer läuft durch. Ursache nicht untersucht.
- **Delta-ML: Modellenergie am eigenen Optimum zu niedrig** (n=20 je Molekül, 2026-10-19): bei `lambda 1e-6` und 4000 Zentren liegt sie um bis zu 34 kcal/mol unter GFN2 (19/19 auf dem Trisaccharid), `lambda 1e-3` beseitigt das und verschiebt die Energie um 7 bis 16 kcal/mol nach oben; `uncertainty_abort` hilft nicht. Details: `docs/DELTA_ML_RETEST_2026-10.md`.
- **Delta-ML aktives Lernen: Abbruchkriterium** (2026-10-27, n=3+3 Läufe peo8): ersetzt durch Pool über `al_window` 3 Zyklen, zentrierte Fehler, 10-%-Quantil, Ziele 2,5/3,5 kcal/mol aus Tuning-Seeds 4-6; `both` = Validierung `al_patience`-mal + Indikator im Stoppzyklus. Löst in 2 von 3 Läufen aus (Zyklus 18, 19), Zeitpunkt vom Indikator bestimmt, JSD wie N0 (0,221).
  Offen: andere Moleküle, ob `validation` allein (stoppt bei 7-9) gleich gute Modelle gibt, Wandzeit +20-38 % ungeklärt. Details: `docs/DELTA_ML_ACTIVE_LEARNING_2026-10.md` (Addendum). Warntext `not converged ... was not met in N consecutive cycles` (`activelearn.cpp:1311`) sagt nicht, ob Validierung oder Indikator fehlte.
- **`charge_groups ionpairs`** (experimentell, Paritätsregel, 2026-10-24): setzt Formalladungen auch in Protonentransfer-TS und Komplexen mit geteiltem
  Proton (WTMAD-2 +0,47; WCPT18/PX13 bis +143 kcal/mol; Kontrollen 6/10 +39 bis +196). Offen: Wahrnehmung Guanidinium/Amidinium, Phosphonium,
  Sulfonium, protonierte Carbonyle, Dianionen; Gruppen je Instanz eingefroren (MD/react nicht geprüft). [docs/FRAG_CHARGE_MODEL.md](docs/FRAG_CHARGE_MODEL.md)
- **Kontinuierliches Ladungsmodell für Ionenpaare** (rev-gfnff, SQE/ACKS2-artig): kein geprüftes GFN-FF-Maß trennt Ammoniumhalogenid-Paare von
  Protonentransfer-TS. Falsifikator-Paket im Vault; dessen `species/*.xyz` passen bei 85/124 Spezies nicht zu den ORCA-Energien (vor Nutzung ersetzen).
- **Fragmentreihenfolge im Mehrfragment-EEQ**: IL16/144 (neutrales Ionenpaar) ändert beim Vertauschen der Fragmentreihenfolge die Energie um
  -2,4e-4 Eh und Ladungen um 7,8e-3 e, ohne Schalter, auf allen Lösern; mit `ionpairs` -4,1e-4 Eh; 230/187/008: 0 (n=1 je, 2026-10-24, nicht eingegrenzt).
- **Coulomb-CN-Kettenregel bei CN = 0, Rest** (2026-10-25): CPU-Schalter `-gfnff.cn_zero_guard` (opt-in, Known Issue #43). Offen: CUDA/ROCm-Kernel
  `k_subtract_qtmp`/`k_coulomb_postprocess` (bisher CPU-Rückfall, GPU-Zweig nicht übersetzt), Bestätigung in pprcht `gfnff_engrad.F90` (nur Vorrunde),
  Default-Entscheidung (Betreiber); Hessian (FD der Gradienten) erbt den Schalter, nicht gelaufen.
- **`frag_charge_model ensemble` bei Polymer^n+ + Cl-** (2026-10-25, n=1 je): Einheit nur +1, daher Q = +1 aus (AB)2 (+2) und Cl- als Polymer 0 / Cl +1
  (Paritätsregel), (AB)3+Cl-: Polymer +1 / Cl +1. NVE 1 ps, 50 K: Etot-Sprünge bis -0,093 Eh (4 > 1e-3), der erste bei H(N)...Cl 1,65 A;
  `reference` und Dateigruppen 1e-5 bzw. 2e-7 Eh/ps. Ursache des Sprungs in der MD nicht eingegrenzt. [docs/FRAG_CHARGE_MODEL.md](docs/FRAG_CHARGE_MODEL.md)
- **`charge_groups` `"rest": "per_fragment"`** (2026-10-25): Gruppen je Instanz eingefroren, geänderte Fragmentierung in der MD wird nicht erkannt
  (nur Codelesen; in keinem Lauf aufgetreten, Kontaktpaar-MD kam nicht unter 1,9 A, 2026-10-25); Fragmente des Skripts aus Kovalenzradien, nicht aus der GFN-FF-Topologie; periodisch, GPU, equilibrierte MD offen.
  Validierung nur PDMAEMA-Einheit/Dimer in Wasserclustern (n=16/8): Kraefte nicht signifikant besser, Relativenergien besser. [docs/POLYMERBUILD_CHARGE.md](docs/POLYMERBUILD_CHARGE.md)
- **100-ps-Probelauf (AB)5 im Wassertropfen** (2026-10-25, n=1 je T, 303/338 K): ausgewertet, offen: Kette wandert in beiden Läufen zur Tropfenoberfläche
  (COM-Versatz 1,1-1,4 -> 3,2-3,5 A, 35 A^2 unbedeckt in H2), 338 K Epot-Drift -1,25e-3 Eh/ps, Rg nicht stationär (tau 5-15 ps, 3-15 unabhängige Werte);
  (`-md.seed` wirkte in diesen Läufen nicht, seit 2026-10-30 behoben.) Läufe unter `~/curcuma_runs/probe100/`. [docs/PDMAEMA_PROBE_100PS_2026-10.md](docs/PDMAEMA_PROBE_100PS_2026-10.md)
- **`-trajprops hb_criterion energy`** baut GFN-FF ohne `charge_groups`/`cn_zero_guard` (`countHydrogenBondsByEnergy`), also nicht das
  Ladungsmodell des Laufs. [docs/PDMAEMA_PROBE_100PS_2026-10.md](docs/PDMAEMA_PROBE_100PS_2026-10.md)
- **H...Cl-Bindungserkennung als Energiestufe** (2026-10-25, starrer Scan n=1): frische GFN-FF-Topologie gibt +34 kcal/mol bei 1,7->1,6 A (ohne Gruppen),
  +50 bei 1,8->1,7 A mit Gruppen (Pass-2-Radien mit q(Cl) = -1), dazu -21 bei 1,2 A (N-Cl); ORCA/GFN2 glatt. In MD nicht erreicht; `auto` baut dort die
  Topologie ohne Bindungsterme neu, Gruppen bleiben ungeprueft (Codelesen). Entscheidung offen. [docs/POLYMERBUILD_CHARGE.md](docs/POLYMERBUILD_CHARGE.md)

## 2. Leistung
- **GFN-FF-CUDA-NPT: erzwungene Listen-Neuaufbauten nach jedem Barostatschritt** (2026-11-23, WP-PERIODIC-NONBONDED S3): bei 3456 H2O kostet ein Schritt mit Zellskalierung 396 statt 96 ms (HB/XB 101 ms, weitere Host-Listen ~200 ms), im Mittel die Hälfte des 278-ms-Schritts; groesster Hebel nach S3. Danach: reziproke Ewald-Summe FP64-gebunden (6,8 ms je pcg-Operator), groesseres `ewald_r_cut` verschiebt Arbeit in den billigen Realraum (aendert Energien ~1e-7 Eh).
- **Periodische CPU-Pfade nach dem Threading-Durchgang** (2026-11-22, WP-PERIODIC-XTB): seriell bleiben Ewald-Punktladungsgradient/-virial und KO-Gradient (~100 ms bei 432 H2O), schlecht skaliert der periodische D3-Paarterm (gfn1 216 H2O 276 -> 179 ms bei 1 -> 8 Threads: Replay, `computeDC6DCN`); der Eigenlöser trägt 38 % des Schritts bei t = 8. Details: [docs/WP-PERIODIC-XTB.md](docs/WP-PERIODIC-XTB.md) "CPU threading pass".
- **GFN-FF-CPU: periodische EEQ ist die Hälfte des Aufrufs** (2026-11-22, `CURCUMA_GFNFF_PHASE`): 432 H2O t = 8 8,1 von 20,2 ms (S8 4,6), 3456 H2O 161 von 315 ms (S8 5,5); der matrixfreie PCG ist nicht weiter zerlegt; Workspace-Vorlauf (HB-CN, Ewald-Aufbau) ~2 / 16-22 ms seriell. Erst zerlegen, dann entscheiden.
- **Eigenlöser-Kappe ohne Zelle offen** (2026-11-22): mit Zelle nutzt der Eigenlöser höchstens die physischen Kerne (-threads 32: 1,4x schneller); ohne Zelle gilt weiter `-threads` (gleicher SMT-Effekt im Basisbericht, 648 Atome), Entscheidung des Betreibers. `docs/SQM_THREADING.md` nennt noch die alte Kappe min(t, 8) (seit Sep 2026 entfernt).

- **Delta-ML-Labels auf der GPU lohnen nicht** (2026-10-27, n=2 Wiederholungen, peo8 + 90-Atom-Wirt, [docs/DELTA_ML_LABEL_GPU_2026-10.md](docs/DELTA_ML_LABEL_GPU_2026-10.md)): beste GPU 0,06x / 0,25x der besten CPU, kein `label_gpu`.
  Offen (Entscheidung): `-activelearn -gpu cuda` schickt die Labels heute schon auf die GPU (`LabelThread` übernimmt den globalen `gpu`-Schlüssel), 4-16x langsamer; Labels bei `-gpu` auf die CPU festlegen?

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
  offen: CG-Iterationen unter Last (K5); Host-Neubau überlappt nicht mit Gerätearbeit (Stand 16. Okt. 2026, 296 Schritte, 20,0 ms/Schritt: HB/XB 2,2 ms/Schritt (31 Refreshs à 25 ms: Erkennung 14,5, Upload 5,1, Bond-HB-Neuaufbau 3; die HB-Liste ändert sich bei jedem Refresh, 49,6k -> 50,4k Paare, Überspringen bei unveränderter Liste greift nicht; Erkennung skaliert schlecht mit Threads, 4/8/16/32: 24,5/21,5/12,8/10,5 ms; `hb_update_force_every` 10 -> 20 spart 1,1 ms/Schritt, ändert aber die Trajektorie, |dE_total| bis 6e-4 über 300 Schritte, n=1; Entscheidung offen; NVE-Messung D1 (2026-11-03): 420 H2O in L = 44 Bohr, 1260 Atome, dt 0,5 fs: Verlust -6,4e-4 Eh/ps mit 10, -2,9e-5 mit 1, n=1, docs/WP-PERIODIC-NONBONDED.md D1), CN-Liste K1 1,0, Repulsion-Geräteneubau 0,9, D4 0,7; Repulsion lief vorher auf dem Host, 3,3 ms/Schritt). ROCm-Rückstand: siehe unten. Journal: Labor/curcuma EEQ-Solver Benchmark.
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

- **GFN2 CPU-MD gegen xtb, Reste** (2026-11-17, [docs/SCF_THRESHOLD_2026-11.md](docs/SCF_THRESHOLD_2026-11.md), `complex`, 1 Thread, n=1-2): Ursache der Lücke sind SCF-Iterationen (9,9 gegen 8), nicht der Eigenlöser. Offen: Dipol-/Quadrupolintegrale `xtb_multipole.cpp:183-219` in einem Durchlauf mit S und H0 (Gewinn 20-30 ms, geschätzt); Precision-Stufen `scf_precision` loose 5e-5 / normal 2e-5 / tight 1e-5 / verytight 1e-6 (Vorschlag).
- **Idee, nach der GFN2-xtb-Klärung: Basiswechsel statt voller Diagonalisierung** (2026-11-17, ungemessen): Fock-Matrix in die Orbitale des Vorschritts transformieren (fast diagonal), Eigenpaare per Jacobi/Störungsrechnung, kleine Rotationen abschneiden; verwandt `-large_system_mode dc|sparse` (dünn besetzte Dichtematrix, lohnt erst bei vielen tausend AOs) und die Dichteantwort in `docs/SCF_PRECOND.md`. Zuerst prüfen, was `partial diagonalisation` und der Warmstart schon tun.
- **Idee: SCF-Vorkonditionierer `-scf_precond frag` auf `complex` messen** (2026-11-17): additive Schwarz-Blöcke je Fragment (204 + 27 Atome), `docs/SCF_PRECOND.md` hat keinen `complex`-Eintrag; Ziel Iterationen je MD-Schritt (aktuell 9,9).

- **Fehlschlagende Tests im Build `USE_MKL=ON, USE_BLAS=OFF` (CUDA, RTX 5080, 2026-11-18, 2 von 709 in `ctest -j8`)**: `cli_confsearch_05_seed_cold_bias` (stabil 2/2), `cli_deltaml_04_label_consistency` (Hilfsskript compare_labels.py fehlt, nicht getrackt). Ursachen nicht untersucht.

- **D4-Gradient der Pfade `-d4_charge_source eeq` und `cpscf` ist etwa 1e-4 Eh/A neben FD** (2026-11-18, n = 3 kleine Moleküle, nur CPU; Diagnose eines Agenten, Codestellen von mir nicht nachgelesen): H2 8,2e-5 obwohl die Antwort dort 0 ist; der Standardpfad `mulliken` liegt bei 1e-8. Vermutete Ursache: das D4-Potential steckt im SCF und in W, fehlt aber im Gradienten-Charge-Pulay dieser Pfade (`xtb_native.cpp:2550-2551`, `xtb_gradient.cpp:385`). Entscheidung: als genähert dokumentieren oder entfernen; Details `docs/BROKEN_TESTS_2026-10-14.md`.

- **Solvatation, eigenes Arbeitspaket** (2026-11-21, Audit gegen xtb 6.7.1 und tblite, [docs/SOLVATION.md](docs/SOLVATION.md) "Not implemented natively"): CPCM (tblite), ddCOSMO (`--cosmo`), CPCM-X (`--cpcmx`), Ionenstärke, Kernel- und Gitterwahl fehlen. Offen: GFN-FF `gbsa` rechnet ALPB (xtb 6.7.1: Still-GBSA), GFN2-Gradient bei `na2cl2` gegen xtb ungeklärt (GFN1: xtb weicht von FD ab), Radius für Z > 94 still 2 A.

## 3. Tests und Validierung

- **`-xtb.gpu_eigen_strategy auto` ungemessen ausser RTX 5080** (2026-11-21, WP-PERIODIC-XTB A1r): Rechenzentrums-GPU (H200/A100: erwartet fp64_direct) und
  FP64-schwache Karten ohne cc >= 12 (A4500, GTX 1660: auto = fp64_direct, ob die cuBLAS-Emulation dort lohnt, offen) je einmal `fp64_direct` gegen
  `fp32_refine` messen und die Regel in `XtbGpuContext::selectEigenStrategy` nachziehen.
- **PLUMED, ungetestet** (Stand 2026-10-30, getestet nur Ar2-RESTRAINT und Ethan-Torsion mit GFN-FF): QM-Methoden, größere Systeme, Restart, NPT/`setVirial`, PLUMED zusammen mit GPU. Gemeldete extreme PLUMED-Langsamkeit nicht nachgemessen. Well-tempered METAD lag 10 % über der freien MD (2.66 gegen 2.43 kcal/mol, Hügel nicht konvergiert).
- **GPU nach dem NEB-MD-Merge** (2026-10-30, RTX 5080, Stand vor `96c29c50`): CUDA bitgleich auf Koffein und Benzol, Vulkan nur ohne `-gpu_strict` geprüft (CPU-Fallback nicht ausgeschlossen), ROCm nicht testbar (kein SDK). NEB-MD auf CUDA nur 1 ps ohne Zahlenvergleich zur CPU.
- **Compiler-Warnungen aus dem Master** (2026-10-30, kein Neubau zum Nachzählen): `measurement.h` (`-Wsuggest-override`), `xtb_native.cpp` (`t_gp0`), `src/helpers/main.cpp` (`columns`), dazu 65 `-Winline`; die Vorgabe "warnungsfreier Build" aus `CLAUDE.md` ist nicht erfüllt.
- **`gfnff_d4_fused_c6` hängt an den Build-Flags** (2026-11-01, gemessen): der Test verlangt
  Bitgleichheit von fused- und stored-D4-Pfad; mit `-DUSE_AVX512/AVX2/BLAS=ON` weichen die
  Gradienten ab, aus demselben Quellbaum ohne die drei sind sie gleich (Details Archiv #45).

- **xtb-Referenzenergien der Teststruktur-Bibliothek** (2026-10-11, `test_cases/reference_data/xtb_library_energies.json`): 66 Strukturen GFN-FF, 64 GFN2 gerechnet (peo201-chain 1410 Atome GFN2: +1,6e-3 Eh, mit `d4_atm_cutoff 40.0` +2,6e-5). GFN-FF: 57 <= 1e-6 Eh, 8 > 1e-4 (hydronium-water2 -0,169 Eh Standard `ensemble` gegen xtb Referenzregel;
  jl22-bgal/macrocycle-bgal 173/231 Atome +0,25 bis +0,80 kcal/mol, nicht erklärt: `storsion_reference_loop_bug` ändert nichts, Term-Vergleich und pprcht/gfnff als dritte Meinung stehen aus). GFN2: max 5,2e-5 Eh, mit `-xtb.d4_atm_cutoff 40.0` auf macrocycle-bgal 5e-5 -> -2e-7 Eh. Offen: ctest gegen die Datei (Toleranzen), 4 Strukturen > 3000 Atome ohne xtb-GFN-FF (xtb wird bei 6200 und 7320 Atomen am 48-GB-Deckel getötet; curcuma 5 s), na2cl2 GFN2 xtb unkonvergiert.

- **NEB-MD Helicen-Neulauf** mit dem Stand nach `e81f1a23` (Trust-Radius, Abbruch bei toten Images): einziger vollständiger Lauf (2026-10-03, n=1) hatte tote Images in Stufe 1
  und Reaktions-dG 2.3-2.5 statt 0; zusätzlich mehrere unabhängige String-Pfade für Ethan (Pfadvarianz ungemessen). Vault: `Projekte/curcuma NEB-MD Pfad-Sampling`.
- `cli_simplemd_08/09` (Essigsäuredimer, CSVR, dt 1 fs) gegen das gemergte Binary mit der MD-Uhr-Korrektur (`ef462fcf`) erneut laufen lassen (Known Issue #32).
- `test_cg_potentials` ist seit 2026-10-01 als ctest `cg_potentials` registriert (besteht); der CG-Beispielaufruf mit VTF-Eingabe endet mit "Failed to initialize ForceField
  engine" (`docs/archive/CG.md`, geprüft 2026-10-01).
- **Zwei ctests scheitern reproduzierbar** (am 2026-10-01 auf `master` und auf `fix/remove-dead-code` gleich, `confscan_dtemplate` 3 von 3 Läufen auf `master`): `xtb_cpscf` (Gate (C), Antwort Mulliken gegen EEQ: H2 8.2e-5 gegen 1e-6, H2O und HCN je über 1.5e-4), `test_orca_interface` ("O2 triplet mult=3" in `methodToOrcaKeyword`), `confscan_dtemplate` (bis 2026-10-09 erwartet 17/2/1/305, tatsächlich 14/5/1/237 unabhängig von `CurcumaThreads` 1/4/24; die Zähler sind jetzt die aktuellen des Binaries und gleich denen von `confscan_free`/`confscan_subspace`; ob `dtemplate` tatsächlich wie `subspace` rechnet oder die alten Werte vor 2026-08 anders zustande kamen, ist nicht geklärt). Die Einordnung "flaky" für `confscan_dtemplate` in älteren Notizen trifft nicht zu.
- **Voller ctest 2026-10-11** (Paket 4, CPU-only ohne `-march=native` und CUDA-Build, je 9 rot, `confscan_dtemplate` grün): `test_orca_interface` scheitert auch bei "O2 quartet mult=5" (nicht nur Triplett); `cli_errors_01` bis `_06` und `parameter_io_tests` codieren `${PROJECT_ROOT}/release/curcuma` fest (grün mit umgebogenem Pfad, ohne `release/` rot);
  Nicht gebaut: MinGW, ROCm (Compiler fehlen).
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

- PLUMED-Vorzeichen und kB*T (`96c29c50`); CSVR chi2(Nf-1) (`521fa5e8`, Ethan 6 x 25 ps: 314.1 -> 297.0 K, Verhältnis 0.9455 gegen 18/19); NEB-MD-Pipeline, nur Ethan reproduziert (4 Läufe, AIChangelog 2026-10-30); `.md.csv`-Ausgabe und PLUMED-Routing.
- CODATA-2018-Einheiten (`191cebe3`; Rest 0,0036 kcal/mol gegen pprcht auf polymer_2x nicht verfolgt), Known Issue #35.
- `eeq_refactor_eps_bohr` Voreinstellung 0 (Faktor-Cache aus) und `nonbonded_skin_bohr` 2,0; Ladungsplatzierung `ensemble` (#31);
  Stillstandserkennung `-opt.stall_steps` 20; `-scf_guess fragments`; `-scf_precond frag` (Okt 2026, opt-in).

## 5. Funktionen und Refactoring (niedrige Priorität)

- **Open-Shell: Spin pro Fragment, GPU und Konstanten-Fragen offen** (2026-11-19, aktualisiert 2026-11-21): `-xtb.spin_polarization true` (tblite-Spinpolarisation, CPU, Broyden, `docs/SPIN_POLARISATION.md`) legt den O2-Triplett unter den Singulett; ohne die Option bleibt `-spin N` ein Orbitalsatz mit zwei Fermi-Besetzungen (Triplett über Singulett, O2 +1,65 kcal/mol). `-spin` gilt global (kein Spin pro Fragment).
  Offen: ROCm/Vulkan (heute Fehler; CUDA-Pfad seit 2026-11-22, host-getriebene Geräteschleife, mechanisch nach `rocm/xtb_hip_context.hip` portierbar, siehe docs/SPIN_POLARISATION.md Abschnitt GPU), voll geräteresidente Schleife mit Spinpolarisation (CUDA), `-large_system_mode`, DIIS/Level-Shift, `-scf_precond` mit Spinpolarisation (heute Fehler); Zitat und LGPL-Lizenzfrage der Konstanten (`parameters/spin_constants.hpp`, Operator); tblite-Pin in `CMakeLists.txt` hat alte Konstanten; mehrere SCF-Lösungen, bei 300 K oft > 500 Broyden-Zyklen.
- **PLUMED-Hilfeseiten überarbeiten** (2026-11-10, Branch `feature/plumed-help-page`, Stash-Stand vom 01.06., Basis `25d4d8d7`): `-help-plumed` und der SimpleMD-Hinweis darauf doppeln `-help metadynamics` (`printMetadynamicsHelp()` in `src/main.cpp`) und `docs/PLUMED_HELP.md`; zusammenführen, auf den aktuellen Stand bringen.
  Nur die Hilfetexte übernehmen: die Einheitenänderungen des Branches (Bohr-Positionen, kB*T, `setBox`) sind älter als `07e7544d`/`96c29c50` und würden sie teilweise zurückdrehen.

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

- **Code-Schulden aus dem Debt-Bericht** (2026-11-14, [docs/debt_2026-11/CODE.md](docs/debt_2026-11/CODE.md), nicht angefasst): zwei GFN-FF-EEQ-Pfade, einer unregistriert (H3); funktionslokale `static`-Zähler in Rechencode (H4); Riesenfunktionen und Dateigrößen (M9, M10); CUDA/ROCm-Copy-Paste und doppelte CPU-Zahlenroutinen (M11, M12);
  zwei abstrakte Basisklassen (M13); 809 `std::cout`/`cerr` am Logger vorbei (M14); 107 `getenv`-Schalter (M15); `exit()` in Bibliothekscode (M16); statische JSON-Defaults in weiteren Headern (M17).
- **Build und Doku aus dem Debt-Bericht** ([docs/debt_2026-11/BUILD_DOCS.md](docs/debt_2026-11/BUILD_DOCS.md), nicht angefasst): `curcuma_cap` kompiliert 38 Dateien ein zweites Mal (M1); ungenutzte Compile-Definitionen, `-g` in jedem Build-Typ, PCH nur teilweise wirksam, ungepinnte Abhängigkeiten (M3 bis M6); Doku-Drift; `-dipole`-Ausgabe wiederherstellen (Punkt 14).
- **Nicht verifiziert nach dem Debt-Branch** (2026-11-14): Builds mit `USE_TBLITE`, `USE_D3`, `USE_GFNFF`, `USE_XTB`, `USE_ULYSSES`, CUDA, ROCm, Vulkan; gfn1/gfn2-Regression (MOR41 nur `uff` und `gfnff`); `ctest` mit `-DCURCUMA_UNICODE_OUTPUT=ON`; Casino hat keinen Test.

- **NEB-MD nach PR #16** (2026-11-14): Nicht geprüft: PLUMED-Korrektur (Vorzeichen, `kB*T`) im Build, GPU, ROCm; NEB-MD nur auf Ethan mit GFN-FF validiert; `nebmd.cpp` hat Funktionen bis 465 Zeilen und `CoutSuppress` (tauscht `std::cout.rdbuf`).

## 6. Entscheidungen und Fragen

Code-gebundene Entscheidungen, die bei Ihnen liegen. Fachliche Fragen und Projektstatus stehen im Vault und sind nur verlinkt.

- Branch `feature/nebmd-path-sampling` nach master mergen? (enthält überwiegend Master-Arbeit; `/pre-merge-check` fällig.) NEB-MD-Status: Vault `Projekte/curcuma NEB-MD Pfad-Sampling`.
- `-scf_pseudo_diag` als Standard einschalten? (Messungen sprechen dafür, Abschnitt 2.)
- `-scf_precond_guard`: Vorgabe `g2:4:best` kostet auf dem unrelaxierten 599-Atom-Cluster bei 8 Threads 4 Iterationen gegen `none` (Kriterium war <= 3);
  über 5 Threadzahlen liegt frag (55-62) unter none (58-138). Varianten streuen 58-109, eine konvergiert nicht. Vorgabe behalten? [docs/SCF_PRECOND.md](docs/SCF_PRECOND.md)
- `eeq_ewald_dense_max_atoms` (1000) senken? Matrixfreies pcg ist bei 900 Atomen 7,2 statt 84 ms (300 Atome: 1,5 statt 16 ms), liefert aber nur `eeq_ewald_tol`-Genauigkeit statt exakt; die kleinen pbc-Tests laufen heute dicht (WP-Doku P1).
- `shouldUpdateHBXB()` korrigieren und damit von der Fortran-Referenz abweichen, oder als Referenzeigenschaft belassen?
- UFF-/QMDFF-Gradient reparieren (Stand: "erstmal nicht")?
- GPU-Bitgleichheit angehen (Sammel-Kernel statt `atomicAdd`) oder bei Ensemblevergleichen bleiben?
- Defaults je Modul für `verbosity`, `threads`, `method`, `charge`, `spin` (Casino-Verbosity 2, Opt 1, Trajektorienanalyse 2): ein Wert aus `global` für alle (ändert Ausgaben) oder die Modul-Vorgaben behalten? ([PARAMS.md](docs/debt_2026-11/PARAMS.md) M3)
- ConfScan: `-method` ist die Alignment-Methode (`confscan.h`, `rmsd.h`), die Energiemethode heißt `energy_method` und ist kein PARAM; `-method gfn2` wird dort als ungültiger Alignment-Name gemeldet (H2). `energy_method` als PARAM einführen und `-method` dafür auswerten?
- `lbfgs_max_line_search` (gemessen 2026-11-14, `scripts/lbfgs_linesearch_grid.py`): 2 stammt aus dem Legacy-`CurcumaOpt` (`LBFGS_ls_iter`), `LBFGSppOptimizer` bekam es am 2026-04-01 (`35afc6ed`), das Registry-PARAM am 2026-04-14 den Default 20 (`3ac3584e`); `-opt` läuft seitdem mit 20. 27 Kombinationen (gfnff 14, gfn1 6, gfn2 7; 8 bis 231 Atome, davon 5 gestörte Starts), Werte 2/3/5/20:
  2 macht in 19 von 27 keinen einzigen Schritt und konvergiert nie, wo 20 scheitert; 3 scheitert in 10; 5 und 20 sind gleich (0 Unterschiede im Konvergieren, größte Energiedifferenz 1e-4 Eh). Nicht verwertbar: `uff` (Gradient, siehe oben), `pm3` und `gfn2` auf dem 117-Atom-System (leere Tabelle bei jedem Wert, Ursache nicht untersucht), 1410-Atom-Kette (300 Schritte, alle gleich).
  Auf `devel/methods` bekommen ConfSearch, ConfGen, MD-Voroptimierung, Docking, PolymerBuild und RMSDTraj den Wert 2 (aus dem Code gelesen, nicht beobachtet); hier gilt überall 20, Test `optimizer_default_layers`. Bei 20 bleiben (Empfehlung) und `devel/methods` informieren?
- Vault, fachliche Fragen: `Offene Fragen/` ("Valenzanteil im reaktiven GFN-FF", "Veraltete C6-Dispersion und Paarlisten in inaktiven Topologie-Ecken",
  "Lücke zwischen GFN2-xTB und Doppelhybrid-Methoden schließen", "RRHO-Näherung durch MD ... ersetzen", "Ladungsplatzierung bei geladenen Mehrfragment-Systemen in GFN-FF - was ist physikalisch besser", Known Issue #31).
- Vault, Projektstatus: `Projekte/curcuma rev-gfnff`, `curcuma GPU-Backends`, `curcuma MD-Stabilität großer Systeme`, `curcuma ConfSearch Weiterentwicklung`,
  `curcuma GFN-FF Transitionsmetalle (MOR41)`, `curcuma Reaktives GFN-FF (react-Topologie)`, `curcuma Konformerenfilter (reorder_paper)`.

- **Parameter-Routing, Rest** (2026-11-14, [docs/debt_2026-11/PARAMS.md](docs/debt_2026-11/PARAMS.md)): erledigt sind H1 bis H7 ohne H2, M1/M7/M8/M9, L1 bis L8 (Warnungen statt Fehler für unbekannte Schlüssel, Werte außerhalb `enum`/`min`/`max` und Typfehler); Modul `global`, Golden-Test `cli_parser_golden` (24 Fälle).
  Offen: Globals in nur einer Kopie (rund 158 Leser, je Capability migrieren), M2 (Top-Level-Skalare im `ConfigManager`), `Json2KeyWord` (123 Aufrufe) und `MergeJson` in Capabilities, `Casino::LoadRestartInformation` ohne Aufrufer, Mehrfach-XYZ-`-opt` schreibt keine `.trj.xyz`.
- **Optimierer-Defaults** (2026-11-14): die Erkennung "Nutzerwert" für `-convergence_preset` vergleicht noch mit dem Default (`optimizer_factory.cpp`), ein explizit gleicher Wert wird verworfen; `OptimizerInterfaceJson` und Registry-Modul `opt` sind zwei Literale, die nur per Kommentar gleich gehalten werden.

## 7. Löschkandidaten

Am 2026-10-01 auf dem Branch `fix/remove-dead-code` entfernt: curcumaopt.cpp, `native_lbfgs_optimizer.*`, `optimisation/modern_optimizer_simple.*`, munkress_2.h,
vier ungenutzte `optimiser/*.h`, core/accuracy_profile.cpp (inkludierte einen nicht vorhandenen Header), `RMSDDriver::MolAlignLib()` und `Prepare*Template()`,
`D4Evaluator::launchGpuKernel()`, helpers/cli_test.cpp und gfnff_test.cpp, zwei untracked `.backup`-Dateien. Volle Testsuite danach: dieselben drei bekannten Ausfälle wie vorher.

Offen, Aufrufer vor dem Löschen erneut suchen:
- `src/capabilities/curcumaopt.h`: hält den PARAM-Block des Moduls `opt` (28 Parameter); die Klasse `CurcumaOpt` ist ohne Implementierung. Die Parameter müssen erst in ein anderes Header umziehen.
- `src/capabilities/rmsd/rmsd_costmatrix.*`, `rmsd_assignment.*` (`CostMatrixCalculator`, `MunkresAssignmentSolver`): kompiliert und von `rmsd_strategies.h` eingebunden, die Klassen selbst werden nirgends aufgerufen. Entscheidung nötig, ob die Strategien sie künftig nutzen sollen.
- dc6dcn-Skip-Pruefung (`FFWorkspaceGPU::shouldSkipDc6dcn()` in `cuda/ff_workspace_gpu.h`, hat keinen Aufrufer; ROCm fuehrt sie mit `CURCUMA_ROCM_LEAN_STEP=0` noch aus).
- Member-Funktionen ohne Aufrufer: Liste mit 215 sicheren Kandidaten ("kann weg"), 47 nur im eigenen Header benutzten und Bibliotheks-/GPU-Headern in `docs/debt_2026-11/DEAD_MEMBERS.md`, erzeugt von `scripts/find_dead_members.py`. Vor dem Löschen Build mit `USE_TBLITE`, `USE_D3`, `USE_GFNFF` und GPU prüfen.
- `test_cases/cli/*/run_test.sh` (81 getrackte Dateien) sind erzeugte Kopien mit altem `PROJECT_ROOT` (`claude_curcuma`); ctest benutzt die im Build erzeugten. Prüfen, ob sie aus git entfernt werden können (`.gitignore`).
- `xtb_cpscf` schlägt fehl (Gruppe C: H2 response 8.2e-5 gegen 1e-6, H2O/HCN isolation 1.7e-4/2.9e-4 gegen 1.5e-4); identisch im Stand vor dem Debt-Branch.
- PARAM `gfnff.accuracy` (loose|normal|medium|high) bleibt: kein Leser bis jetzt, soll mit den Genauigkeitsprofilen (Eintrag oben) eine Bedeutung bekommen. Nicht entfernen; der Hilfetext sagt bis dahin "no effect".

- Restliche Tabellen-Duplikate: `d4_ncoord.cpp` (`rad`, `en`, verbatim cpp-d4, Pauling-EN gleich `element_tables.h` fuer Z <= 86) und die Device-Kopien in `cuda/`/`rocm/` (laden die Host-Tabellen, nicht geprueft); Konsolidierung nur mit GPU-Test.

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

