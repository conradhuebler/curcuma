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
- Neighbour search by a cell list; no periodic images, so no bond across a box boundary is found.
- A geminal contact inside a molecule (D...A within a few bonds) is counted when it meets the criteria; the
  angle limit removes most of them.

## Status

🤖 AI-generated, ⚙️ machine-tested.

Checked: `test_trajprops` (36 checks against hand-worked geometry: equal and unequal masses, a square,
subsets, frame ranges, statistics of 1..10, a ten-frame file with a selection, a selection matching
nothing; water dimers linear, bent and stretched against the distance and angle limits, group directions,
element restriction, the complement selection, a three-frame hbonds run). A generated 28,000-atom, 100-frame XYZ (71 MB) runs in 3.8 s at 33 MB peak memory; the R_g of a
rigid translation stays constant; hbonds over 30 frames of that file take 7 s at 33 MB (random
coordinates, the count itself means nothing).

Not checked: hydrogen-bond counts on real water or on a polymer in water, comparison with another program's count, real MD trajectories, VTF input, files with a changing atom count, timing and memory beyond
that size, comparison with another program's R_g. Human production testing pending.

Planned in the same capability (not implemented): `rdf`, `sasa` (water-accessible surface,
Shrake-Rupley); an `energy` criterion for `hbonds` from the GFN-FF hydrogen-bond terms (needs the per-triple
energies of `FFWorkspace::calcHydrogenBonds`, which are only summed today, CPU only) and logging of the
count during an MD run. Work package item D in qurcuma `docs/WP-min2015-reproduction.md`.
