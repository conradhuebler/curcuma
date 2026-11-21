# Reactive-MD Example Inputs (GFN2, CPU)

> Status: AI-generated (Nov 2026), machine-tested. Human production testing pending. No reference calculation (DFT, tblite, xtb) exists for any run.
> **Forced scenarios (Haber-Bosch) are not free dynamics.** A bias potential pulls the N atoms apart; the result shows a path into the NH3 sink under that bias, no barrier, no rate, no thermodynamics.

Fifteen inputs for interactive reaction MD (eight listed first, seven Doebereiner variants in their own section) (GFN2, 1 fs, `csvr`, CPU) that a viewer such as qurcuma can show. Each is a JSON for `-import_config` (group `simplemd`, `_command: md`, [CLI_ROUND_TRIP.md](CLI_ROUND_TRIP.md)) in `test_cases/reaction_md/`, with a start structure `reactions-*` in `test_cases/structures/` (provenance `constructed`, atom order is load-bearing for the atom indices in the JSON).

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
| Doebereiner, cold Pt7 | `reactions-doebereiner-pt7` (13: Pt7, 2 H2, O2) | `reactions-doebereiner-pt7.json`, Pt region 500 K, gas region 2000 K, `spin 2`, wall R 8 -> 5.5 (6000 steps), 12000 steps | free, wall radius ramp | H2 split and O2 stretched on Pt, no H2O gain | 27 | see "cold Pt7" below | not a catalytic effect on H2O |
| Doebereiner, Pt7 control | `reactions-doebereiner-pt7-control` (6) | `reactions-doebereiner-pt7-control.json`, same gas, spin 2, same wall | free | gas phase | 2900 | see below | same wall radius as Pt7 |
| Doebereiner, cold Pt7, spin-polarised | `reactions-doebereiner-pt7` (13) | `reactions-doebereiner-pt7-spinpol.json`, as cold Pt7 but `spin 8` and `xtb.spin_polarization true` | free, wall radius ramp | H2 split, O2 stretched, no free H2O | 16 (loaded machine) | see "spin polarisation" below | CPU only; spin 2 aborts on the SCF |
| Doebereiner, Pt7 control, spin-polarised | `reactions-doebereiner-pt7-control` (6) | `reactions-doebereiner-pt7-control-spinpol.json`, spin 2 (O2 triplet) with `xtb.spin_polarization true` | free | gas phase, O2 triplet unreactive | 1400-2500 | see below | same wall radius as Pt7 |

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

## Doebereiner lighter, cold Pt7 (16 seeds, Pt against gas control)

> Status: AI-generated (Nov 2026), machine-tested. Human production testing pending. No reference calculation. Result: Pt7 splits H2 and stretches O2 significantly more often than the gas, **H2O formation is not increased** (1/16 in every variant).

Changes against the Pt4 section: (1) Pt7 pentagonal bipyramid (Pt-Pt 2.75 A) instead of Pt4; (2) the Pt atoms are held cold with a thermal region: global `temperature 500`, `temp_regions` gas atoms (8:13) at 2000 K ([TEMPERATURE_RAMP.md](TEMPERATURE_RAMP.md)), so the metal does not melt/fly apart (the earlier "cluster flies apart" came from starting all atoms at 2000 K; with this protocol Pt stays one connected cluster, Pt-Pt 3.4 A criterion, in all 64 Pt7 runs with Pt at 500 K: Pt Rg 2.12 A at start, at most 2.52 A; at Pt 1000 K up to 2.82 A, one run (F, spin 0) with two Pt components in some frames); (3) 16 start geometries: 2 H2 and O2 at random orientation 3.0-3.6 A from the cluster within a 55 degree cap (generator seed 1000+i), seed i = geometry i; (4) the control is the same gas atoms without Pt (centroid at the origin) with the same wall schedule R 8 -> 5.5 A in 6000 steps (same wall radius, so same wall; the free volume differs by the Pt7 volume); (5) 12000 steps (12 ps), spin 0 and spin 2 for Pt and control. Other settings as the Pt4 section (GFN2, CSVR, 1 fs, `electronic_temperature 2000`, `scf_damping 0.15`, `geometry_abort_factor 0`, logfermi, 1 thread). Free dynamics, no restraint. Analysis: same bond criteria as above (1.25 x covalent radii, free H2O = fragment H2O not bonded to Pt, O with two H = H2O unit also on Pt), batch evaluation, Wilson 95 % intervals and two-sided Fisher exact test (scripts were in the scratchpad, copy in the lab journal tarball, not in the repository).

| 16 seeds each | Pt7 spin 0 | control spin 0 | Pt7 spin 2 | control spin 2 |
|---|---|---|---|---|
| both H2 split within 1 ps | 12/16 | 2/16 | 14/16 | 3/16 |
| O-O above 1.65 A at least once | 15/16 | 11/16 | 13/16 | 1/16 |
| any O-H bond | 7/16 | 16/16 | 10/16 | 10/16 |
| O-H in at least half of the frames | 6/16 | 8/16 | 4/16 | 1/16 |
| O with two H (H2O unit) | 1/16 | 9/16 | 1/16 | 1/16 |
| free H2O | 1/16 (11.6 ps) | 1/16 (7.8 ps) | 1/16 (10.6 ps) | 1/16 (2.9 ps) |

Fisher (Pt against control, same spin): H2 split p = 0.001 (spin 0) and p < 0.001 (spin 2); O-O stretched p = 0.17 (spin 0, singlet gas also breaks O-O) and p < 0.001 (spin 2); H2O unit p = 0.006 (spin 0, Pt has *fewer* than the gas, which forms H2O2-like species) and p = 1.0 (spin 2); free H2O p = 1.0 in all comparisons. Pt at 1000 K (gas 2000 K, 16 seeds, same geometries): spin 2: H2 split 9/16 (against control p = 0.066), O-O stretched 14/16 (p < 0.001), H2O unit 1/16, free H2O 1/16 (4.3 ps); spin 0: H2 split 12/16, O-O 14/16, H2O unit 2/16, free H2O 1/16. Pt 500 K against 1000 K (spin 2): no difference (all p >= 0.11).

Findings: Pt7 catalyses the first steps (H2 dissociation, O2 stretching with the triplet occupation) at high significance; the product H2O does not appear more often than in the gas (1 of 16 in each of six variants, confidence interval 0.01-0.28, so an effect of about 3x or more cannot be excluded or shown). Oxygen and hydrogen stay adsorbed on the cold metal in 12 ps; with spin 0 the Pt suppresses the O-H bond formation that the singlet gas shows. A catalytic effect on H2O formation is **not established**; whether longer times, hotter metal with a stable cluster, or a larger cluster change this is not tested.

Speed (CPU, 1 thread, 8 jobs in parallel, machine carrying a foreign job at load 24-32 on 24 threads): Pt7 + 6 atoms 26.6-27.6 steps/s (below the target of 30; 12000 steps in about 450 s), gas control 2600-2900 steps/s.

Inputs: `reactions-doebereiner-pt7.json` (geometry 1, spin 2, Pt cold) and `reactions-doebereiner-pt7-control.json` (same gas, spin 2, wall identical); structures `reactions-doebereiner-pt7`, `reactions-doebereiner-pt7-control`. Geometry 1 is not representative: in this one the cold Pt7 splits both H2 (640 fs) and stretches O2 (9.1 ps) without any O-H bond, while its spin-2 gas control forms H2O at 2.9 ps. Test: `cli_simplemd_18_reaction_md_inputs` F checks only the H2 split on Pt7 in 1200 steps.

Not tested: other geometries as inputs, clusters larger than Pt7 or surfaces, Pt above 1000 K, times longer than 12 ps, wall radius against density (the wall was identical in Pt and control but the free volume is not), 1000-1500 K gas, spin-polarised O2 (spin 2 is only an occupation), comparison with a reference calculation or a different method. Caveats: GFN2 over-binding of Pt-H and Pt-O is not checked; cold Pt and hot gas is not an equilibrium situation; counts of one event in 16 runs are not suitable to rank variants.

## What `-spin` does in native GFN2, and Pt7 at spin 0, 2, 4

> Status: AI-generated (Nov 2026), machine-tested. Human production testing pending. Result: no spin treatment available in the native code makes the catalyst effect on H2O formation visible (free H2O 1/16 in every Pt7 variant, 0-1/16 in every control).

**What `-spin N` is** (`xtb_native.cpp` Z. 3151-3160, `xtb_scf.cpp` Z. 658-684, `xtb_native.h` Z. 1369-1372): N is the **number of unpaired electrons** (multiplicity - 1, tblite `nuhf`): `-spin 2` is a triplet, `-spin 4` a quintet. It sets n_alpha = (n_el - N)/2 + N and n_beta = (n_el - N)/2. The SCF keeps **one** set of orbitals and one Fock matrix and fills them with two Fermi fillings (n_alpha, n_beta), summed to one density (tblite `nspin == 1`). It is not a UHF with separate alpha and beta Fock matrices, and the code has no spin-polarisation term (no spin constants W; searched `src/core/energy_calculators/qm_methods/`, none found, as in GFN2 of xtb/tblite). `-spin 0` is the closed-shell Fermi-smeared singlet. Consequence: a triplet is only a constrained occupation of the same orbitals, so it lies **above** the singlet; the physical O2 ground state (triplet, about 22 kcal/mol below the singlet) cannot be obtained. Native GFN2 only; TBLite/xtb are not built here (`USE_TBLITE=OFF`), so a spin-polarised reference was not available.

Spin energy ladder, `-sp`, GFN2, `electronic_temperature` 2000 (O2 also 300), Eh as printed:

| System | spin 0 | spin 2 | spin 4 | spin 6 | spin 8 |
|---|---|---|---|---|---|
| O2 1.21 A, 300 K | -7.90675236 | -7.90411828 | -7.14534257 | | |
| O2 1.21 A, 2000 K | -7.92167879 | -7.90411829 | -7.14680292 | | |
| Pt7 (pentagonal bipyramid) | -32.02907763 | -32.02660670 | -32.01891122 | -32.00486430 | -31.97890440 |
| Pt7 + 2 H2 + O2, geometry 1 | -41.94941114 | -41.94733139 | -41.94092682 | -41.92960890 | -41.91194893 |
| geometry 2 | -41.96167194 | -41.95958330 | -41.95315518 | -41.94180988 | -41.92413772 |
| geometry 3 | -41.96656953 | -41.96446404 | -41.95799319 | -41.94661180 | -41.92899727 |

Singlet is lowest in every case: triplet minus singlet is +1.65 kcal/mol (O2, 300 K), +1.55 (Pt7) and +1.31 (full system, geometry 1); the energy rises monotonically with spin. At 2000 K the singlet O2 energy is lowered by the smearing (11.0 kcal/mol below the triplet), while spin 2 gives the same energy at 300 and 2000 K (fixed occupation). Consequently there is no "energetically best" physical spin state; spin 0 is the GFN2 minimum and spin 2 is the constrained triplet that blocks the artificial fractional pi* occupation of the singlet O2.

Runs (16 geometries each, same protocol, seeds and geometries as the cold-Pt7 section; Pt region 500 K, gas region 2000 K, 12000 steps), new here: **Pt7 at spin 4** (Pt7 triplet plus O2 triplet in the one-orbital-set picture) and the **gas control at spin 4**. Spin 4 for the control is not physical (O2 quintet is 0.76 Eh above the triplet); it is listed only for a like-for-like test.

| 16 seeds each | Pt7 s0 | gas s0 | Pt7 s2 | gas s2 | Pt7 s4 | gas s4 |
|---|---|---|---|---|---|---|
| both H2 split within 1 ps | 12 | 2 | 14 | 3 | 10 | 16 |
| O-O above 1.65 A | 15 | 11 | 13 | 1 | 15 | 2 |
| any O-H bond | 7 | 16 | 10 | 10 | 8 | 8 |
| O with two H | 1 | 9 | 1 | 1 | 1 | 0 |
| free H2O | 1 | 1 | 1 | 1 | 1 (0.6 ps) | 0 |

Fisher (two-sided): free H2O Pt7 against gas p = 1.0 at spin 0, 2 and 4 (1/16 against 1/16, 1/16, 0/16). Pt7 s4 against gas s2: O-O stretched p < 0.001, H2 split p = 0.029, H2O unit and free H2O p = 1.0. Pt7 s4 against Pt7 s0 and s2: all p >= 0.22. Multiple comparisons uncorrected. The only O-H-related difference is the one already reported (singlet gas forms H2O units 9/16, Pt7 1/16, p = 0.006). Spin 4 in the gas splits H2 16/16 because the forced occupation excites the H2 molecules, not a catalytic effect.

Speed: Pt7 + 6 atoms 26.7-29.3 steps/s (12000 steps in 410-449 s, 16 Pt jobs and 6 control jobs in parallel on 24 threads; machine load average 15 at the end, 22 own threads, so a foreign load of about 14 threads at times); target 30 missed again. Control 3.3 s per run (3600 steps/s). No SCF error or abort in 32 runs (settings electronic_temperature 2000, scf_damping 0.15, scf_max_iterations 500, geometry_abort_factor 0; no countermeasure needed, no "false fixed point" warning in the checked logs).

Not tested: a UHF-type or spin-polarised GFN2 (implemented later the same day, see the spin-polarisation section below), spin state per fragment (Pt7 singlet with O2 triplet is not expressible: `-spin` is global and with one orbital set the unpaired electrons go to the highest orbitals, wherever they are), other spins for the gas control, Pt7 spin 6/8, reference calculation. Caveat: because singlet GFN2 is the lowest state, the observed difference between spin 0 and spin 2 gas (O-O broken 11/16 against 1/16) is an effect of the occupation constraint, not of real O2 spin chemistry.

## Doebereiner lighter, cold Pt7 with spin polarisation (16 seeds)

> Status: AI-generated (Nov 2026), machine-tested. Human production testing pending. No reference calculation of the dynamics. Result: **no catalyst effect on H2O formation (0/16 free H2O in Pt7 and in the gas control)**; Pt7 still splits H2 and stretches O2 significantly more often than the gas.

Same protocol, geometries and seeds as the cold-Pt7 section (Pt region 500 K, gas region 2000 K, wall R 8 -> 5.5 A in 6000 steps, 12000 steps, `electronic_temperature 2000`, `scf_damping 0.15`, `scf_max_iterations 500`, `geometry_abort_factor 0`, 1 thread), but with `xtb.spin_polarization true` ([SPIN_POLARISATION.md](SPIN_POLARISATION.md): tblite spin-polarisation model, CPU and Broyden only). Spin: the `-sp` ladder uhf 0/2/4/6/8/10 of all 16 start geometries (GFN2, 2000 K, spin polarisation) has its minimum at uhf 8 in every geometry (uhf 6 0.6-4.5, uhf 4 0.06-5.6 kcal/mol higher; in geometries 2, 5, 6, 14 uhf 4 is within 0.11 kcal/mol; uhf 2 did not converge in geometries 4 and 9), so Pt7 runs use the fixed spin 8. The gas has no spin-8 state; it uses spin 2 (O2 triplet) and, for classification, spin 0 (singlet, unchanged by the option).

| 16 seeds each | Pt7 s8 | gas s2 | gas s0 (singlet) | Pt7 s2, `scf_allow_unconverged` |
|---|---|---|---|---|
| both H2 split within 1 ps | 10 [0.39-0.82] | 2 [0.03-0.36] | 2 | 15 [0.72-0.99] |
| O-O above 1.65 A | 13 [0.57-0.93] | 1 [0.01-0.28] | 11 | 16 [0.81-1.00] |
| any O-H bond | 9 [0.33-0.77] | 6 [0.18-0.61] | 16 | 9 [0.33-0.77] |
| O-H in at least half of the frames | 5 | 1 | 8 | 2 |
| O with two H | 0 | 0 | 9 | 0 |
| free H2O | 0 [0.00-0.19] | 0 [0.00-0.19] | 1 (7.8 ps) | 0 [0.00-0.19] |
| runs ended by SCF abort | 0 | 0 | 0 | 0 (15/16 with the default, see below) |

Fisher (two-sided, uncorrected for multiple comparisons): Pt7 s8 against gas s2: free H2O and H2O unit p = 1.0, any O-H p = 0.48, O-H persistent p = 0.17, O-O stretched p < 0.001, H2 split p = 0.009. Pt7 s2 against gas s2 (same total spin): free H2O and H2O unit p = 1.0, any O-H p = 0.48, O-O p < 0.001, H2 split p < 0.001. Pt7 s8 against Pt7 s2: all p >= 0.083. Gas s0 reproduces the earlier result exactly (1/16 free H2O at 7.76 ps, 9/16 H2O units, 16/16 O-H, 11/16 O-O, 2/16 H2 split), as expected since the option does not act at spin 0. Previous results without spin polarisation (Pt7 s0/s2/s4 free H2O 1/16 each; gas s0/s2/s4 1/16, 1/16, 0/16) are not different within the 16-seed intervals.

Findings: (1) with a real O2 triplet (15.8 kcal/mol below the singlet in this model) the gas is **less** reactive than the singlet gas: O-O stretched 1/16 against 11/16, any O-H 6/16 against 16/16, H2O units 0/16 against 9/16 (p <= 0.001); the earlier occupation-only spin 2 showed the same (1/16, 10/16, 1/16), so the spin polarisation changes little there. (2) Pt7 (s8 and s2) catalyses H2 dissociation and O2 stretching as before. (3) No free H2O and no H2O unit in 32 Pt7 runs; O-H bonds on Pt7 are not more frequent than in the triplet gas (p = 0.48). A catalyst effect on H2O formation within 12 ps is not shown (upper Wilson bound 0.19 per variant). Pt7 s8: Pt stays one connected cluster in 16/16 runs (3.4 A criterion, Pt Rg at most 2.45 A).

SCF: Pt7 s8 and all gas runs: no abort in 48 runs. **Pt7 s2 with spin polarisation aborted in 16/16 runs** (SCF not converged after 501 iterations, "Simulation got unstable", within the first 250 steps (aborted trajectories hold 2-24 frames of 10 steps)). Settings tried on geometry 1 (300 steps) did not remove it: `scf_damping` 0.05 and 0.3, `scf_max_iterations` 2000, `scf_guess h0`, `scf_mixed_precision false`, `scf_extrapolation aspc`, h0 with FP64; only `scf_allow_unconverged true` continues (warning at the SCF level is not visible in MD output, so the number of unconverged steps was not counted; the abort rate of about 2 % per step in the default runs is a lower guide). The table column Pt7 s2 is that run: forces from an unconverged SCF are not consistent energy derivatives, and the Pt cluster split into two components for part of the time in 6/16 runs (one run for most of the time, Pt Rg up to 4.2 A), so this column is not a clean comparison; the Pt7 s8 column is. Spin 2 on Pt7 is an excited magnetisation (0.7-5.6 kcal/mol above uhf 8) with several competing SCF solutions.

Speed: Pt7 + 6 atoms with spin polarisation 15.4-19.0 steps/s (12000 steps in 630-780 s) with 16 own jobs and a load average of 22-37 on 24 threads (foreign load about 6-20 threads); the 26 steps/s of the 1000-step test were at load 25 with 1 job; gas control 1400-2500 steps/s. The cost of the option itself was not isolated in this run.

Inputs: `reactions-doebereiner-pt7-spinpol.json` (geometry 1, spin 8; start energy -41.961434 Eh) and `reactions-doebereiner-pt7-control-spinpol.json` (spin 2); structures as the cold-Pt7 inputs. The `xtb` group in the JSON selects the option; without it the run is the old one. Test: not added to `cli_simplemd_18_reaction_md_inputs` (a 1200-step Pt7 run would take about 60 s); the two JSONs were run for 100 steps each.

Not tested: another spin per geometry (spin 8 is the minimum only of the start geometry, not along the trajectory), spin per fragment, times beyond 12 ps, other Pt cluster sizes, a reference calculation (tblite was not run on the dynamics), GFN2 over-binding of Pt-H and Pt-O, a closed-shell Pt7 with an O2 triplet (not expressible, `-spin` is global). Caveat: 16 seeds from 16 start geometries; zero events in 16 do not exclude an effect of a factor of 3.

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
