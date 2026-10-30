# WP4 — DFT-Quadratur-Grid (Treutler-Ahlrichs / Euler-Maclaurin + Lebedev + Becke)

- **Status:** 🤖 AI-generated, ⚙️ machine-tested (`ctest -L qm_grid` 19/19, Okt 2026); menschlicher Produktionstest ausstehend
- **Abhängigkeit:** WP3 (HF-Dichte für den N_el-Test)
- **Code:** `src/core/energy_calculators/qm_methods/qm_grid.{h,cpp}`, Tabellen `qm_grid_tables.hpp`
  (erzeugt von `scripts/qm_grid_export_tables.py` aus PySCF 2.14), Engine-Anschluss `QMEngine::buildGrid`

## Was gebaut ist

- Pro Atom radiales Gitter x Lebedev-Laikov-Winkelgitter, Becke-Partition, **kein Pruning**; Bohr, a.u.
- Radial: `treutler` (Default, Treutler-Ahlrichs M4, Chebyshev 2. Art, elementspezifisches xi = PySCF `radi.treutler`)
  oder `euler_maclaurin` (Murray-Handy-Laming, R = Bragg-Radius, xcDFT `EulMac`; in PySCF nicht vorhanden).
- Winkel: alle 32 Lebedev-Gitter 6 ... 5810 Punkte; Orbit-Generatoren (code, a, b, v) aus PySCF `LebedevGrid.SphGenOh`
  mitgeschnitten, `sphGenOh` in C++ erzeugt dieselbe Punktreihenfolge.
- Partition: `becke` (Becke 1988, f dreimal) oder `stratmann` (SSF-Polynom, a = 0.64, ohne SSF-Screening);
  Größenanpassung `treutler` (sqrt Bragg-Radius, PySCF-Default), `becke` (Bragg-Radius), `none`.
  Achtung: a_ij = 1/4 (r_j/r_i - r_i/r_j) ist PySCFs Formel, nicht Beckes Anhang-Formel a = u/(u^2-1).
- Stufen `-qm.grid`: `sg0` 23x170, `sg1` 50x194 (Default), `sg2` 75x302, `sg3` 99x590, frei `<n_rad>x<n_ang>`.
  Die SG-n-Namen übernehmen nur die Punktzahlen von xcDFT `read_grid.f90`; es sind **nicht** die gepruneten SG-n-Gitter.
- Weitere Parameter: `-qm.grid_radial`, `-qm.grid_partition`, `-qm.grid_radii_adjust` (Modul `qm`).

## API für WP5/WP6 (Doxygen in `qm_grid.h`)

- `QMEngine::buildGrid(qmgrid::MolecularGrid&, why)`: Gitter der aktuellen Geometrie mit den `-qm.grid*`-Einstellungen
  (`-qm.threads` für die Partition); kein Cache in der Engine, das entscheidet WP5.
- `qmgrid::MolecularGrid`: `points` (npts x 3), `weights`, `quadrature_weights` (vor der Partition, für WP8),
  `atom`, `atom_offset`; Index = `atom_offset[A] + i*n_angular + j`.
- `qmgrid::AOEvaluator(engine.gtoBasis(), engine.sphericalTransform())`: `evaluate(points, begin, count, deriv, AOBlock&)`
  liefert phi und mit `deriv=1` dphi/dx,y,z im aktiven Basis (= Basis von `QMEngine::density()`), const und
  threadsicher, blockweise aufrufen.
- `qmgrid::densityOnBlock(ao, P, rho, &grad)`: rho und grad rho; `integrateDensity(grid, ao, P, threads)`:
  Summe in Blockreihenfolge, Ergebnis unabhängig von der Threadzahl (benzene sg1/sg3 bei 1/3/8 Threads bitgleich).

## Validierung (Okt 2026, def2-SVP, HF-Dichte mit scf_threshold 1e-10)

Elektronenzahl, `test_qm_grid nelec <xyz> --threads 4` (Treutler/Becke/Treutler):

| Molekül | N | sg0 | sg1 | sg2 | sg3 |
|---|---|---|---|---|---|
| H2O | 10 | 4.7e-05 | -3.0e-06 | 3.1e-07 | -5.6e-09 |
| CH4 | 10 | 1.4e-04 | -1.4e-05 | -3.5e-07 | 4.4e-08 |
| Benzol | 42 | 2.1e-04 | 9.3e-05 | 2.6e-05 | 1.8e-06 |
| NH3 | 10 | 1.8e-05 | -2.9e-05 | -1.2e-06 | 1.8e-08 |
| HF | 10 | 1.3e-05 | -1.0e-06 | 3.5e-08 | -5.5e-09 |
| LiH | 4 | 2.8e-04 | 3.2e-06 | -1.2e-07 | -1.1e-08 |
| Ne | 10 | 4.1e-05 | -3.0e-10 | 4.2e-13 | -1.1e-14 |

Gegen PySCF 2.14 (`gen_atomic_grids(prune=None)` + `get_partition`, gleiche Einstellungen; volle Punktliste,
einmalig gemessen): Koordinaten <= 3.6e-15 Bohr, |dw| <= 8.2e-13, Gewichtssummen und Testintegrale <= 3.5e-13
relativ, ∫rho identisch bis 2.2e-13 (10 Fälle). Rund 1e-4 der Punkte haben Partitionsfaktoren an der
Rundungsgrenze (1 - g mit g -> 1): dort ist ein Programm 0, das andere <= 4.4e-15; relative Gewichtsfehler dieser
kleinen Gewichte sind deshalb nicht aussagekräftig, der Test prüft |dw|/w_quad <= 1e-12.
Lebedev: alle reellen Y_lm mit l <= Grad auf <= 9.7e-14 exakt, bei l = Grad+1 nicht (Fehler >= 1.9); Orthonormalität
der Y_lm bis l = Grad/2 auf <= 2.3e-14 (Gitter <= 590 Punkte). AO-Gradienten gegen 5-Punkt-FD <= 7.5e-10 relativ.
rho, grad rho an Stichprobenpunkten gegen PySCF <= 2.6e-9 relativ (Konvergenz der beiden SCFs).

## Nicht getestet / nicht implementiert

- Elemente Z > 10 (Basis nur H-Ne); offene Schalen; geladene Moleküle; große Moleküle (Partition O(npts natoms^2),
  kein SSF-Screening, kein Pruning, keine Gitter-Screening-Listen für AOs).
- `euler_maclaurin` nur gegen die eigene Formel in PySCF-Umgebung geprüft (keine unabhängige Referenz).
- Gitterkonvergenz von E_xc (WP5), GPU. Gitterantwort im Gradienten: erledigt in WP8 (`qmgrid::gridWeightGradient`).
- Ein ungültiger `-qm.grid`-Name fällt erst beim Bau des Gitters auf; HF baut keins.

## Quellen

A. D. Becke, JCP 88, 2547 (1988); O. Treutler, R. Ahlrichs, JCP 102, 346 (1995); C. W. Murray, N. C. Handy,
G. J. Laming, Mol. Phys. 78, 997 (1993); R. E. Stratmann, G. E. Scuseria, M. J. Frisch, CPL 257, 213 (1996);
V. I. Lebedev, D. N. Laikov, Dokl. Math. 59, 477 (1999); PySCF `dft/gen_grid.py`, `dft/radi.py`, `dft/LebedevGrid.py`.
xcDFT-Quellen (`dft_grid.f`, `read_grid.f90`) lagen in dieser Sitzung nicht vor (`/home/conrad/src/claude_curcuma/xcDFT`
fehlt); die EM-Formel stammt aus der alten Fassung dieser Datei.
