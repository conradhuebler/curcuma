# Reactive-MD Example Inputs (GFN2, CPU)

> Status: AI-generated (Nov 2026), machine-tested. Human production testing pending. No reference calculation (DFT, tblite, xtb) exists for any run.
> **Forced scenarios (Haber-Bosch) are not free dynamics.** A bias potential pulls the N atoms apart; the result shows a path into the NH3 sink under that bias, no barrier, no rate, no thermodynamics.

Eleven inputs for interactive reaction MD (eight listed first, three Doebereiner variants in their own section) (GFN2, 1 fs, `csvr`, CPU) that a viewer such as qurcuma can show. Each is a JSON for `-import_config` (group `simplemd`, `_command: md`, [CLI_ROUND_TRIP.md](CLI_ROUND_TRIP.md)) in `test_cases/reaction_md/`, with a start structure `reactions-*` in `test_cases/structures/` (provenance `constructed`, atom order is load-bearing for the atom indices in the JSON).

```bash
cd <any empty directory>
R=<repo>; curcuma -md $R/test_cases/structures/clusters/reactions-wgs-gas.xyz -import_config $R/test_cases/reaction_md/reactions-wgs-gas.json
python3 $R/scripts/reaction_md_check.py wgs <output dir> --dt 1 --dump 10     # bond analysis (scenarios: h2-metal soup wgs haber fe4 doeb)
```

`_input` in each JSON is the repository-relative structure path, so `curcuma -import_config <json>` also works with the repository root as working directory. Add `-seed N`, `-threads N`, `-max_time` on the command line (CLI wins). Output goes to the BMT directory `<name>.md.<timestamp>/` (trajectory `<name>.trj.xyz`, status `<name>.md.csv`); with `-no_bmt` the same files land in `<name>.snapshots/` of the working directory (both checked).

## Scenarios

Steps/s: wall time of the verification run (incl. ~1 s start) over steps, 4 threads, machine shared with other jobs (load 5-22 on 24 threads); 1 thread gives the same rate for the small systems. Seed 1 is the stored seed; seeds 2 and 3 were run with `-threads 1`. "Event" times are from the 1.25 x covalent-radii criterion (`scripts/reaction_md_check.py`), frames every 5 or 10 fs.

| Scenario | Structure (atoms) | Input | Kind | Expectation | Steps/s | Finding of the verification | Limits |
|---|---|---|---|---|---|---|---|
| H2 on Pt4 | `reactions-h2-pt4` (6) | `reactions-h2-pt4.json`, 300 K, 1000 steps | free | H2 dissociates, Pt-H ~1.5 A, adsorption heat heats the cluster | 155 | 3/3 seeds: H-H > 2 A after 25, 25, 30 fs; Pt-H 1.5 A; Pt-Pt max 3.0 (seed 1), 4.3 and 3.8 A (seeds 2, 3): the cluster is distorted, partly loosened | Pt4 idealised, not optimised; heating not quantified (CSV prints every 100 steps) |
| H2 on Pt7 | `reactions-h2-pt7` (9) | `reactions-h2-pt7.json` | free | as above | 47 | seed 1: H-H > 2 A after 35 fs, final 2.15 A, Pt-H 1.54/1.57 A, Pt-Pt max 2.97 A | one seed only |
| Primordial soup | `reactions-soup45` (45: 3 CH4, 3 NH3, 4 H2O, 3 H2) | `reactions-soup45.json`, CSVR 4000 K, wall R 6.5 logfermi, `electronic_temperature 1500`, 6000 steps | free | bond breaking, CH3/H2N/HO, new C-C/C-O bonds from ~5 ps | 109 (4 thr), 93 (1 thr, 10000 steps) | seed 1 (6 ps): 30/32 start bonds broken, 9 new pairs, first new C-N at 5.8 ps, no C-C, no C2H2. Seed 2 (10 ps): 32/32 broken, C-C at 3.9 ps (C2H4, C2H6 at 4.0 ps), 69 new pairs. Seed 3 (10 ps): no C-C, C-O at 6.5 ps (CH2O, CHO, CO) | CPU only (on GPU the SCF aborted in 4/4 runs, earlier campaign); atoms leave the soft wall (max 7.4, 9.6, 12.2 A from origin); no C2H2 seen here; 4000 K and 0.56 g/cm3 are unphysical |
| WGS gas | `reactions-wgs-gas` (15: 3 CO, 3 H2O) | `reactions-wgs-gas.json`, 4000 K, wall R 3.5 logfermi `wall_temp 4000`, 8000 steps | free | CO2, H2, intermediates CH2O2/CHO2 | 398 | 3/3 seeds CO2 (0.65, 4.0, 4.9 ps), H2 in 2/3 (1.6, 5.7 ps); CH2O2/CHO2 in 3/3; CO2 not permanent (CO2 in 5 %, 40 %, 39 % of frames) | the earlier 12-seed statistics (12/12 CO2, 10/12 H2) were not repeated; hot H2 above 0.775 A is missed by the criterion; HOCO and formate not distinguished |
| WGS gas, compressed | same structure | `reactions-wgs-gas-compress.json`: `wall_radius 6.0`, `wall_radius_schedule "3.5:steps:3000"` | free, wall radius ramp | as above, pressure column rises | 398 | seed 1: radius 6.0 -> 3.5 A (CSV), CO2 first at 3.3 ps (59 % of frames), no H2, C2HO2 formed; instantaneous wall pressure 0 GPa at start, 1.4-4.1 GPa after the compression (noisy) | one seed; needs the binary with `wall_radius_schedule` ([WALL_RADIUS_SCHEDULE.md](WALL_RADIUS_SCHEDULE.md)); pressure is a wall force measure |
| Haber-Bosch gas | `reactions-haber-bosch-gas` (8: N2 + 3 H2) | `reactions-haber-bosch-gas.json`, 1000 K, wall R 4.5 harmonic, `distance_harmonic` atoms 1-2 (N) k 0.3 Eh/A2 r0 3.2 A, 3000 steps | **forced** | 2 NH3 | 3300 | 3/3 seeds: first N-H after 60-90 fs, first free NH3 after 620, 440, 680 fs, 2 NH3 in 30/30 final frames, N-N max 3.8 A | the pull is active the whole run; k 0.2 gave 0/10 in the earlier campaign; no free dynamics, no barrier |
| Haber-Bosch gas, control | same structure | `reactions-haber-bosch-gas-control.json` (same, no `external_potentials`) | free | N2 intact | 3300 | 3/3 seeds: N-N 1.06-1.15 A, no N-H bond, 6 free H at the end | N2 splitting needs far more than 3 ps at 1000 K |
| Haber-Bosch Fe4 | `reactions-haber-bosch-fe4` (12: Fe4, N2, 3 H2) | `reactions-haber-bosch-fe4.json`, 1000 K, `electronic_temperature 2000`, `scf_damping 0.15`, tether `centroid_harmonic` k 5 per Fe (atoms 1-4), pull atoms 5-6 k 0.3 r0 3.2 | **forced** | N-N stretched, NH/NH2 on the metal, no NH3 | 113 | 3/3 seeds: N-N max 3.85, 3.71, 3.78 A (final 3.3, 3.3, 3.1), N-H on one N (NH or NH2), no free NH3; Fe-Fe max 2.7-2.8 A; Fe displacement from the tether target max 0.15 A (seed 1) | Fe4 spin 0 (GFN2 has no spin polarisation, Fe is not magnetic here); the tether is a model for a surface; Fe13 excluded (SCF aborts) |
| Doebereiner, Pt4 | `reactions-doebereiner-pt4` (10: Pt4, 2 H2, O2) | `reactions-doebereiner-pt4.json`, 2000 K CSVR, spin 0, `electronic_temperature 2000`, wall R 6.5 -> 4.5 logfermi (`wall_radius_schedule "4.5:steps:10000"`), 20000 steps, 1 thread | free, wall radius ramp | H2 dissociates on Pt, O2 activated, OH and H2O on Pt | 86 | see "Doebereiner lighter" below | see below |
| Doebereiner, control | `reactions-doebereiner-control` (6: 2 H2, O2) | `reactions-doebereiner-control.json`, same settings, no Pt, spin 0 | free | gas phase, less reaction | 2700 | see below | closed-shell singlet O2 is reactive without a catalyst |
| Doebereiner, control (triplet) | same structure | `reactions-doebereiner-control-triplet.json`, `spin 2` | free | no stable O-H | 2700 | see below | spin 2 occupation, not a spin-polarised GFN2 |

Verification: every input run from an empty directory with `release/curcuma` built from this working tree (`ci-feature-delta-ml-170-g267dc325` plus the uncommitted wall-radius change); 24 runs, none with an SCF error or abort; CPU FP32 "false fixed point" warnings occurred and were handled by the program (continued in FP64). `-electronic_temperature` from the JSON is effective (Fe4 step-0 energy -21.0963 Eh at 2000 K against -21.0300 Eh at 300 K). Atom indices in `external_potentials` are 1-based (selection grammar): N-N pull stretches N-N in gas and Fe4, the control does not; the Fe tethers hold. Test: `ctest -R cli_simplemd_18_reaction_md_inputs` (4 s: forced vs control NH3 and H2 on Pt4).

## Doebereiner lighter (2 H2 + O2 on Pt4, with and without catalyst)

> Status: AI-generated (Nov 2026), machine-tested. Human production testing pending. No reference calculation. The statistics are small (8 seeds) and the catalytic effect on H2O formation is **not** established.

Start: Pt4 tetrahedron (as `reactions-h2-pt4`, not optimised) with 2 H2 (3.0 A) and O2 (3.2 A) above one facet, molecules 3.5 A apart, 10 atoms. Control: the same 6 gas atoms without Pt in the same wall. Everything is free dynamics (no restraint). Analysis: `python3 scripts/reaction_md_check.py doeb <run directory> --dump 10` (1.25 x covalent-radii bonds, H-H "intact" below 1.0 A, O-O "broken" above 1.65 A, "free H2O" = a fragment of formula H2O). Settings: GFN2, `thermostat csvr`, 1 fs, `electronic_temperature 2000`, `scf_damping 0.15`, `scf_max_iterations 500`, `geometry_abort_factor 0`, logfermi wall, 20000 steps (20 ps), 1 thread, CPU. No SCF abort in any of the runs (CPU FP32 "false fixed point" warnings occurred and were handled by the program: 0 in the Pt4 replay, 2 and 5 in the control replays).

Result at 2000 K with the wall compressed from R 6.5 to 4.5 A (stored inputs; seeds 1-8, 20 ps each, from the bond criterion above):

| | Pt4 (spin 0) | control spin 0 | control spin 2 |
|---|---|---|---|
| both H2 split within 1 ps | 7/8 | 3/8 | 1/8 |
| O-O above 1.65 A at least once | 8/8 (130 fs - 4.3 ps) | 7/8 | 0/8 |
| any O-H bond | 8/8 | 8/8 | 8/8 (O-H in at most 6 % of the frames) |
| O-H in at least half of the frames | 8/8 | 5/8 | 0/8 |
| O bound to two H at least once | 3/8 | 6/8 | 0/8 |
| free H2O formed | 3/8 (6.5, 8.7, 10.9 ps) | 1/8 (2.3 ps) | 0/8 |

Findings, without interpretation beyond the numbers:

- **H2 dissociation on Pt is robust**: both H-H bonds above 1.0 A within 1 ps in 22 of 23 Pt4 runs of 20 ps (all temperatures, compression and spin variants); first Pt-H contact 30-100 fs in the 5 ps runs. In the gas the H2 molecules mostly stay H2 at 1000-1500 K.
- **O2 is activated on Pt** (Pt-O contact after 80-120 fs, O-O stretched above 1.65 A in 8/8 at 2000 K compressed) and is not in the triplet gas (0/8). In the singlet gas control O-O breaks in 7/8, in H2O2-like species.
- **The singlet control is not inert**: closed-shell O2 (GFN2 has no spin polarisation, `spin 0`) reacts with H2 without catalyst (H2O2-like species, one free H2O in 8). This is a model artefact (real O2 is a triplet; not checked against a reference how reactive H2 + O2 would be under these conditions); with `spin 2` the gas shows only transient O-H contacts. `spin 2` fixes the occupation, it is not a spin-polarised calculation, and was not tried for Pt4 with compression (Pt4 `spin 2` at 2000 K without compression, 3 seeds: O-H in 0.5-1.0 of the frames, no H2O).
- **H2O formation on Pt is slow and not clearly catalytic here**: 3/8 Pt4 runs form a free H2O in 20 ps (6.5-10.9 ps), the singlet control 1/8, the triplet control 0/8; with 8 seeds the 3/8 against 1/8 is not significant. At 1000 and 1500 K Pt4 (3 seeds each, 20 ps, no compression): H2 split and O2 sometimes stretched, H2O 0/6.
- **Without compression** (R 6.5 fixed, 2000 K, 3 seeds, 20 ps): one free H2O (19.9 ps), O-H persistent in 1/3. Compression raised the frequency of OH on Pt (8/8) and of H2O (3/8) in these runs; not tested whether the wall radius or the density is the cause.
- **The Pt4 cluster does not stay a cluster**: Pt-Pt (start neighbours) up to 5.2-6.0 A in the compressed 2000 K runs, 4.8-5.1 A even at 1000 K; adsorption heats it and GFN2 over-binds Pt-H and Pt-O (not checked against a reference). Temperatures in the CSV range about 1200-3400 K at a target of 2000 K (instantaneous values, printed every 500 steps) in a 10-atom system.
- Pt7 (planar 7 + 6 atoms, 2 seeds at 1500 and 2000 K, 8000 steps): H2 split in 3/4 within 1 ps, O-O stretched in 3/4, one O-H bond (briefly, one run), no H2O; Pt-Pt up to 6-7.7 A (cluster flies apart); about 29 steps/s with 8 jobs plus foreign load running, below Pt4 (86); not included as an input. These four runs were started while the driver script was being edited (exit status 127 recorded after a complete run); not repeated.
- Speed (CPU, 1 thread, the machine carrying a foreign 18-thread job): Pt4 86 steps/s (20000 steps in 234 s, stored JSON replay, bitwise equal energy to the exploratory run), control 2700 steps/s; GPU `-gpu cuda -gpu_strict true` on Pt4: 16 steps/s (500 steps, start included, free GPU, host load about 25), slower than the CPU as expected below about 50 atoms. 4 threads were not faster than 1 under the foreign load (Pt4 87 against 107 steps/s, 300 steps).

Caveats: 2000 K and 20 ps are far from a Doebereiner lighter (room temperature, hydrogen on platinum sponge in air; the heat of H2O formation heating the metal to glowing is not modelled, the thermostat removes it); the wall compression raises the density of 10 atoms to a sphere of R 4.5 A; energies shifted by the 2000 K Fermi smearing; no reference calculation; the metal is a distorted Pt4, not a surface; seeds 1-8 are all from one start geometry (gas placement not varied).

Not tested: Pt cluster sizes other than 4 and 7, other metals, a gas placement away from the facet, temperatures between 1500 and 2000 K with compression for more than 3 seeds, GPU other than the single Pt4 throughput run, spin 2 for Pt4 with compression, and a restraint-assisted variant.

## What these are not

- 4000 K, the 0.56 g/cm3 soup and the 3.5 A WGS sphere (6 molecules in 180 A3, about 33 per nm3) are far from real conditions; reactions occur because of temperature and density.
- Pressure is only the wall radius: there is no barostat, `wall_pressure_GPa` is a force per wall area, not an equilibrium pressure. Periodic code is out of scope.
- GFN2 without spin polarisation; radicals and free atoms (H, HO, CH3) run closed shell with Fermi smearing (`electronic_temperature` 1000-2000 K is a deliberate deviation from GFN2 that shifts energies by mEh, measured in the earlier campaign); metal clusters use spin 0.
- Fe and Pt clusters are over-bound by GFN2 (earlier campaign, not checked against a reference).
- Statistics are 1-3 seeds; "expectation" columns come from the earlier campaign (12 seeds for WGS, 10 for Haber-Bosch) and were not reproduced to that extent.

## qurcuma

Read in `~/src/qurcuma` (no code changed, not run): `qurcuma <file.xyz> -md` loads the structure and starts MD with the settings of the Simulation dock. No loader for a curcuma `-import_config` JSON was found; the Recipe menu stores qurcuma's own `SimulationConfig` protocol (`src/recipe.h`), and `*.qlesson.json` bundles structures and configurations. So each scenario is started by loading the structure and setting method, temperature, thermostat, wall by hand (free scenarios: H2/Pt, soup, WGS). `electronic_temperature`, `scf_*` and `geometry_abort_factor` are not hand fields; they can be entered in the "All parameters" tab (`mdExtraParams`, passed as `simplemd` parameters). `external_potentials` is not in the hand groups (`docs/WP-ux-restructure.md` line 88), `setExternalPotentials` is not called anywhere in qurcuma, and its bundled curcuma (`external/curcuma`, Oct 14) predates `wall_radius_schedule`. Whether the generated "All parameters" tab accepts the `external_potentials` JSON was not tried. The forced Haber-Bosch scenarios therefore need code work or an untested workaround in qurcuma; TODO.md has the entry.

## Not implemented

- Time window for `external_potentials` (`start_step`/`stop_step`), unrestricted open-shell per fragment, qurcuma wiring: [TODO.md](../TODO.md).
- Not tested: Pd/Ni, Fe5-Fe7, Cu/Pt WGS (no reaction in 36 earlier runs), GPU (slower than CPU for these sizes, soup aborts), other thermostats, dt other than 1 fs, thread counts above 4.
