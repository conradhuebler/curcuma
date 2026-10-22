# Streaming trajectory properties (`-trajprops`)

Properties of long trajectories (tens of thousands of atoms, thousands of frames) with one frame in
memory at a time. `Measurement` keeps all frames and `UnifiedAnalysis` all selected frames; both are
meant for small systems. Code: `src/capabilities/trajprops.{h,cpp}`, test: `test_cases/test_trajprops.cpp`.

```
curcuma -trajprops traj.xyz -selection F1 -last_fraction 0.3333     # R_g of fragment 1 over the last third
curcuma -trajprops traj.xyz -selection "1:780" -stride 10 -mass_weighted false
```

## Quantity `gyration`

- R_g^2 = sum m_i |r_i - R|^2 / sum m_i about the centre of mass R, of the selected atoms (all atoms when
  `selection` is empty). `mass_weighted false` uses unit weights.
- `selection` is resolved once, on the first frame used; the atom count must not change along the file.
- Frames: `first_frame`..`last_frame` (-1 = end), every `stride`-th; `last_fraction` in (0, 1] takes the last
  fraction of the file and overrides the range (0.3333 = the last 10 ns of 30 ns).
- Output: `<basename>_gyration.csv` (frame, R_g in Angstrom) and `Results()`: `mean`, `stddev` (over frames),
  `min`, `max`, `median`, `block_error`, `blocks`, `frames_used`, `frames_in_file`, `first_frame_used`,
  `last_frame_used`, `atoms`, `series_file`.
- `block_error` is the standard error of the mean from `blocks` consecutive blocks. Frames of an MD run are
  correlated, so `stddev / sqrt(frames_used)` understates the uncertainty of the mean; the block error does
  not assume independent frames but needs blocks longer than the correlation time.
- Coordinates are used as stored. A chain wrapped around a periodic boundary gives a wrong R_g, and an XYZ
  file carries no cell to unwrap with.

## Quantity `hbonds`

Geometric hydrogen bonds D-H...A. A hydrogen belongs to the nearest donor atom (element in `hb_elements`,
default N,O) within `hb_bond_length` (1.3 A). The bond counts when D...A <= `hb_distance` (3.5 A) and the
angle D-H...A >= `hb_angle` (150 deg). All three limits and the elements are parameters.

- `selection` (group A) and `selection_b` (group B): bonds between the groups, counted per direction
  (`a_donates`, `b_donates`); `selection_b` empty counts bonds inside A; both empty, inside the whole system.
  Each (hydrogen, acceptor) pair counts once. A leading `!` selects the complement: `-selection F1
  -selection_b "!F1"` is chain against everything else.
- Output: `<basename>_hbonds.csv` (frame, total, a_donates, b_donates) and `Results()` with the statistics of
  the total per frame plus mean, standard deviation and block error of each direction, and the criteria used.
- `hb_criterion energy` counts GFN-FF hydrogen-bond triples (A-H...B) whose energy is below `hb_energy_limit`
  (kJ/mol, negative = attractive; 0 counts every attractive triple) instead of applying distance and angle.
  The triple energies come from `GFNFF::hydrogenBondEnergies()` ([HBOND_ENERGIES.md](HBOND_ENERGIES.md)); the
  groups, directions and statistics are as above. `charge` is the total charge of the system (GFN-FF needs it for
  the EEQ charges), `threads` the threads of the GFN-FF evaluation.
- Extra output in energy mode: column `energy_sum_kj_mol` (sum of the energies of all triples in the groups) in the
  per-frame CSV, `<basename>_hbonds_energy_histogram.csv` (triples per frame in `hb_hist_bins` bins over
  `[-hb_hist_max, 0]` kJ/mol) and in `Results()` the statistics of the sum plus the counts per frame of triples
  above zero and below the histogram range. The histogram is the tool for choosing the limit: compare the
  resulting count with the geometric count on the same trajectory.
- Cost: every frame is a full GFN-FF set-up and evaluation (topology, EEQ charges), far more than the geometric
  count. One frame of a lattice of waters took 0.4 s and 0.3 GB at 1536 atoms and 2.4 s and 2.1 GB peak at 5184
  atoms (4 threads); 28,000 atoms were not tried, and the memory use grows steeply with size. Use `stride`.
- Neighbour search by a cell list; no periodic images, so no bond across a box boundary is found.
- A geminal contact inside a molecule (D...A within a few bonds) is counted when it meets the criteria; the
  angle limit removes most of them.

## Quantity `rdf`

Radial distribution function g(r) between two atom groups and the running coordination number n(r).

- Groups: `selection` (A) and `selection_b` (B), each narrowed by `rdf_elements_a` / `rdf_elements_b` (element
  symbols). With neither `selection_b` nor `rdf_elements_b` the groups are identical (A with itself, the density
  uses N-1). Two different groups must be disjoint; a partial overlap is refused.
- An XYZ file carries no cell, so give an orthorhombic one: `-cell 40,40,40` (Angstrom). It switches on the minimum
  image convention and supplies the volume for the density of B. Without it, `-rdf_volume V` (A^3) supplies the
  density and no distance is wrapped. `rdf_r_max` must not exceed half the shortest cell edge. Triclinic cells are
  not supported.
- g(r) = pairs(r) / (frames * N_A * rho_B * shell(r)), rho_B = N_B / V (N_B - 1 for identical groups), shell(r) the
  volume between the bin edges; n(r) is the cumulative number of B neighbours per A atom up to the bin's upper
  edge. The error column is the standard error of g from `blocks` blocks of consecutive frames.
- Output: `<basename>_rdf.csv` (r, g, g_error, coordination) and `Results()` with the group sizes, frames, cell
  and `coordination_at_r_max`.
- Cost: N_A * N_B distances per frame (no cell list). Measured on a generated 28,000-atom file with random
  coordinates (unphysical, timing only): 780 x 27,220 atoms, cell 60 A, r_max 10 A, 10 frames took 7.7 s in total,
  about 0.7 s per frame, at 33 MB. Water against water would be some thirty times more pairs; use `stride`.

## Quantity `sasa`

Solvent-accessible surface area of the selected atoms by the Shrake-Rupley method.

- Each atom i carries `sasa_points` test points on a sphere of radius r_i + `sasa_probe` (Fibonacci lattice); a
  point counts when it lies outside the (r_j + probe) sphere of every other atom j; the atom's area is
  4 pi (r_i + probe)^2 times the accessible fraction. Probe 1.4 A is a water molecule.
- Radii: Bondi (1964) for H 1.20, C 1.70, N 1.55, O 1.52, F 1.47, P 1.80, S 1.80, Cl 1.75, Br 1.85, I 1.98 A.
  Any other element is an error until `-sasa_radii "Na:2.27,..."` gives it. Which radius set the reference
  study used is not known here (see the TODO below).
- `sasa_occluders group` (default): only the selected atoms cover each other, i.e. the area of the solute as a
  probe sees it from outside, with the solvent removed. `all`: every atom of the system covers.
- Output: `<basename>_sasa.csv` (frame, area in A^2) and `Results()` with mean, standard deviation over frames,
  block error, the probe, the points and the radius set.
- Cost: for each selected atom, one pass over the covering atoms to find its neighbours plus `sasa_points` tests
  against them. Measured on the same generated file (780 selected atoms, 10 frames): 2.6 s with the group as
  occluders and 3.1 s with all 28,000 atoms as occluders, about 0.3 s per frame.

## References and TODO

TODO (open): the references below were written down by the AI from memory and have not been checked against the
papers; authors, volume, pages and DOIs must be verified before they are quoted. The same entries sit in
`src/core/citation_database.cpp` (printed at the end of a run for the methods that were used) and in the header of
`src/capabilities/trajprops.h`. Tracked in `TODO.md`.

| Method | Reference (as supplied, unchecked) |
|---|---|
| Block error | Flyvbjerg, H.; Petersen, H. G., J. Chem. Phys. 1989, 91, 461 (10.1063/1.457480) |
| Geometric hydrogen bond (D...A <= 3.5 A, angle >= 150 deg is the usual choice, cf. the water criterion) | Luzar, A.; Chandler, D., Nature 1996, 379, 55 (10.1038/379055a0); Kumar, R.; Schmidt, J. R.; Skinner, J. L., J. Chem. Phys. 2007, 126, 204107 (10.1063/1.2742385) |
| Hydrogen bond by GFN-FF triple energy | Spicher, S.; Grimme, S., Angew. Chem. Int. Ed. 2020 (GFN-FF; entry `gfnff` of the database) |
| g(r), n(r) | Allen, M. P.; Tildesley, D. J., Computer Simulation of Liquids, 2nd ed., Oxford University Press, 2017 |
| Shrake-Rupley | Shrake, A.; Rupley, J. A., J. Mol. Biol. 1973, 79, 351 (10.1016/0022-2836(73)90011-9) |
| Probe sphere / accessible surface | Lee, B.; Richards, F. M., J. Mol. Biol. 1971, 55, 379 (10.1016/0022-2836(71)90324-X) |
| van der Waals radii | Bondi, A., J. Phys. Chem. 1964, 68, 441 (10.1021/j100785a001) |
| Test points on a sphere (Fibonacci lattice) | Gonzalez, A., Math. Geosci. 2010, 42, 49 (10.1007/s11004-009-9257-x) |
| Counter-ions by replacing solvent (`-ionize`) | GROMACS genion: Abraham, M. J. et al., SoftwareX 2015, 1-2, 19 (10.1016/j.softx.2015.06.001) |
| The study this work reproduces | Min, S. et al., Soft Matter 2015 (10.1039/c4sm02242d) |

## Status

🤖 AI-generated, ⚙️ machine-tested.

Checked: `test_trajprops` (72 checks against hand-worked geometry: equal and unequal masses, a square,
subsets, frame ranges, statistics of 1..10, a ten-frame file with a selection, a selection matching
nothing; water dimers linear, bent and stretched against the distance and angle limits, group directions,
element restriction, the complement selection, a three-frame hbonds run; the energy criterion on a bonded and a stretched water dimer, group directions, limits, a three-frame run with histogram; g(r) and n(r) of a simple cubic lattice against an independent count of lattice vectors, two disjoint groups, the error between blocks, refusals; the surface of one atom, of two overlapping atoms against the analytic cap formula, covering by group and by all, missing radii). A generated 28,000-atom, 100-frame XYZ (71 MB) runs in 3.8 s at 33 MB peak memory; the R_g of a
rigid translation stays constant; hbonds over 30 frames of that file take 7 s at 33 MB (random
coordinates, the count itself means nothing).

Measured: the Shrake-Rupley area of two overlapping O atoms (distances 1.5, 3.0, 4.5, 5.5 A, 960 and 2000 points) deviates from the analytic value by at most 0.22 % (worst case of these eight); many-sphere geometries such as a polymer were not compared with another program. The lattice g(r) averaged over 4 to 5.5 A is 1.012 against the expected N/(N-1) = 1.005 (finite lattice, tolerance 5 %).

Not checked: g(r) or the surface area on a real solvated polymer trajectory, the choice of radii against the reference study, comparison with another program's RDF or SASA, whether the energy and the geometric criterion agree on real water (no calibration done), the energy criterion on a solvated polymer, a non-zero `charge` run, hydrogen-bond counts on real water or on a polymer in water, comparison with another program's count, real MD trajectories, VTF input, files with a changing atom count, timing and memory beyond
that size, comparison with another program's R_g. Human production testing pending.

Planned in the same capability (not implemented): water relaxation / residence times, a cell list for large
RDF groups, triclinic cells. Logging the energy-based hydrogen-bond count during MD exists as `-hb_log_frequency`
([HBOND_ENERGIES.md](HBOND_ENERGIES.md)). Work package item D in qurcuma `docs/WP-min2015-reproduction.md`.
