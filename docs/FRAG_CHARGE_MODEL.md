# GFN-FF fragment-charge placement (`-gfnff.frag_charge_model`)

🤖 AI-generated, ⚙️ machine-tested — human production testing pending. Known Issue #31 in `CLAUDE.md`.
Code: `src/core/energy_calculators/ff_methods/gfnff_frag_charge.cpp`; PARAMs in `gfnff.h`.
Ported (Sep 27, 2026) from the `reactff2-llm` branch, where the design and validation record is
`test_cases/revgfnff/_log/FRAG_CHARGE_STATUS.md`.

## The problem

GFN-FF constrains every perceived fragment's EEQ charge sum to an integer. For a **charged**
molecule that falls into several fragments, the reference (pprcht / xtb) puts the whole net charge
on fragment 0, i.e. on the fragment that contains **atom 1 of the input file**. Two defects follow:

1. **Index dependence.** GMTKN55 `WATER27/H3OpH2O2` lists a water oxygen first, so the +1 lands on
   a water and the hydronium is neutral. Re-ordering the file changes the energy by 105.7 kcal/mol.
2. **Discreteness.** Where the pass-1 bond threshold splits a charged pair, the charge jumps from
   delocalised to fully localised: Cl2- steps by 98.2 kcal/mol over 2e-4 A.

## The model (`ensemble`, default)

Every integer placement of the net charge over the fragments is a complete GFN-FF evaluation with
all charge-dependent parameters consistent with that placement (a `GFNFF` variant instance with a
fixed fragment override). Placements are chosen by chemistry, never by index:

1. **Electron-count parity**: no fragment stripped of all its electrons, then the fewest odd-electron
   fragments. A closed-shell ion next to closed-shell neutrals always wins (H3O+ in water).
2. Chemically **different** carriers left tied are weighted by their free Phase-1 EEQ charge,
   softmax width `frag_charge_sigma` (0.05 e). Absolute GFN-FF energies of differently charged
   fragments are not comparable, so energies are never used here.
3. Chemically **identical** carriers (the two ends of a symmetric X...X-) are weighted by their
   energies, softmax temperature `frag_charge_tau` (1 kcal/mol).

Optional **continuous window** (`frag_charge_s_max` > 1, opt-in): fragment pairs closer than
`s_max` times their pass-1 bond threshold are blended (smootherstep in r/r_thr, 2^k topology
corners, at most `frag_charge_max_edges`) with a merged corner, so the energy stays continuous
where the fragment count changes. The analytic gradient includes the blend weights.

`reference` reproduces the pprcht/xtb rule bit for bit. Engaged only for charge != 0 with two or
more fragments, on the CPU path, and not in `topology_mode react`.

## Measured on master (Sep 27, 2026)

All numbers from this port's binary vs the unmodified master binary (`c4cde7f2`), every structure
run in a fresh directory with `-gfnff.cache_topology false`.

| check | result |
|---|---|
| `reference` mode vs old binary, GMTKN55 / MOR41 / S30L-CI | 2462/2462, 95/95, 90/90 energies **and** gradients identical to 12 digits |
| default (`ensemble`, `s_max` 1.0) vs old binary, MOR41 / S30L-CI (neutral) | 95/95, 90/90 identical |
| default vs old binary, GMTKN55 | 2398 identical; 64 move (all charged): 19 by > 1 kcal/mol, 3 by 0.004-0.13, 42 by < 1e-4 |
| same effect on `reactff2-llm` (327 charged structures) | same 64 structures; per-structure agreement <= 6.2e-6 kcal/mol |
| WATER27 reaction MAD vs published reference | **58.58 -> 21.37** kcal/mol |
| BH76 / PArel | 56.44 -> 48.50 / 22.96 -> 22.55 |
| BH76RC (cost) | 43.67 -> **49.07** |
| GMTKN55 WTMAD-2 | 94.09 -> 92.60 (window `s_max` 1.1: 91.50) |

The 19 large movers are all 11 WATER27 ion clusters, 7 BH76 structures (6 SN2 complexes whose
charge sat on the CH3X leaving group, plus `hoch3fts`, where F- and OH- tie on parity and are
mixed 0.83 / 0.17 by rule 2) and `PArel/h2s2o72`.

**Why the window is off by default**: at `s_max` 1.1 it moves `SIE4x4/he2+_1.0` by -233 kcal/mol
and worsens SIE4x4 (293.4 -> 295.8), DIPCS10 and G21EA, because its merged corner inherits plain
GFN-FF's over-delocalised one-fragment charge state. None of these move at the default `s_max` 1.0.
BH76RC does get worse at the default: the old wrong placement had been partly cancelling GFN-FF's
unphysical ion energetics in those reaction energies.

## Candidate selection and memory (Oct 10, 2026)

**Defect (fixed).** The placements were found by enumerating all C(|q|+G-1, |q|) distributions of the
net charge over the G fragment groups, each a vector of length G, before the parity filter and the cap
`frag_charge_max_placements` acted. Measured with an instrumented build on a 2822-atom PEO/water/NaCl
system with 497 fragments (cgroup limit 12 GB): `-charge 1` 497 distributions, 3.7 s, 4.1 GB;
`-charge 2` 123 753, 3.9 s, 4.6 GB; `-charge 3` 2.06e7, killed after 5.1 s; `-charge 4` 2.57e9 (155 GB
seen before the OOM kill). The same input in `reference` mode runs normally.

**Selection now** (`GFNFF::fragSelectPlacements`): rank (bare nuclei x (G+1) + radicals) and score S
are sums of per-group terms, so a dynamic programme over (group, units used) finds the minimum rank
exactly, keeps only states on minimum-rank paths, and per state drops a prefix only if at least
`frag_charge_max_placements` others beat it for every possible suffix (higher-or-equal S and earlier in
the old enumeration order, or higher S by more than the rounding bound). The result is the old list,
same members and same order including the tie-break by enumeration order, at O(G |q|^2) states.
`-gfnff.frag_charge_enumerate true` runs the old enumeration (`fragEnumeratePlacements`, refused above
2e6 distributions) for validation.

**Variant memory.** Every placement is a full GFN-FF instance: VmRSS growth per variant 64-74 bytes x
N^2 (n = 12 variants: 0.48-0.51 GiB at N = 2822 on top of a 0.72 GiB master, 2.45-2.64 GiB at N = 6200,
urea400-water1000 `-charge 1`, 21 GB peak). `frag_charge_variant_memory_mb` (default 8192) bounds the
variants held at once, with an estimate of 72 bytes x N^2 per variant: above it fewer placements per corner (lowest S dropped, the
same choice as a smaller `frag_charge_max_placements`), then the lightest corners, and variants unused
in the current step are released; one warning. Below about 4400 atoms nothing changes at the default
`s_max` 1.0. With the window (`s_max` 1.1) the 497-fragment system now runs in 5.2 GB (8 corners x 1
placement); before it exceeded 16 GB.

| 2822 atoms, 497 fragments | before: time, peak RSS | after |
|---|---|---|
| `-charge 1` / `2` | 3.7 s, 4.1 GB / 3.9 s, 4.6 GB | 3.5 s, 4.1 GB / 3.8 s, 4.1 GB |
| `-charge 3` / `4` | killed at 12 GB after 5-6 s | 3.7 s, 4.1 GB / 3.9 s, 4.1 GB |

**Checked** (machine-tested): `gfnff_frag_charge_selection` (20 000 random cases G <= 12, |q| <= 4,
ties and near-ties, 0 differences from the enumeration); vs the binary before the change (`a2b7c5f0`),
energy and gradient dumps byte-identical for all 327 charged GMTKN55 structures in `ensemble`,
`reference` and `ensemble` with `s_max` 1.1 (981 runs), and for 50 runs on library-derived charged
clusters and bulk systems (up to 2910 atoms, charge -3..+3); `cli_gfnff_08_frag_charge_many_fragments`.
**Not checked**: the budget path with the window in MD (variants rebuilt when evicted), periodic systems.

## Caveats

- **The hardcoded fallback is the real default.** `GFNFF::GFNFF(const json&)` reads each setting with
  `m_parameters.value(key, FALLBACK)`, and the CLI does not merge full registry defaults into
  `controller["gfnff"]`. Changing the PARAM default in `gfnff.h` alone has no runtime effect.
- **GPU**: the GPU wrapper does not call `GFNFF::Calculation()`, so it cannot apply the ensemble. Since Oct 14, 2026 a charged
  multi-fragment system with `-gpu` and `frag_charge_model` not `reference` runs on the CPU method: `GFNFFGpuMethodImpl` sets
  `requestsCpuFallback()` after its setup, `EnergyCalculator::setMolecule()` reports a GPU fallback (counted, `-gpu_strict true`
  exits with code 3) and replaces the method. `-gfnff.frag_charge_model reference` keeps the GPU with the same rule on both sides.
  Measured on ROCm (n = 1 per structure): GMTKN55 gfnff, ROCm against CPU 0 of 2462 differ (25 before); `cl2-anion` and
  `hydronium-water2`, 100 MD steps, last frame identical. Tests `cli_gfnff_gpu_05_rocm_charged_fragments`. CUDA uses the same
  wrapper code and was not run.
- Cost: each variant is a full GFN-FF evaluation; tiny for the parity-decided cases, 2^k corners
  times placements with the window on.
- A reused calculator (optimisation, MD) on master still has the stale CN / D4-C6 defect of
  energy-only calls (fixed on `reactff2-llm`, not part of this port); the ensemble inherits it.
- Tested: GMTKN55, MOR41, S30L-CI, the Cl2- / H3O+(H2O)2 cases of `cli_gfnff_05_frag_charge_ensemble`,
  `gfnff_frag_charge_history`, and up to 503 fragments in `cli_gfnff_08_frag_charge_many_fragments`.
  **Not tested**: periodic systems, charged systems with more than a few fragments in a window, long MD
  with the window on.

## Ion pairs: `-gfnff.charge_groups ionpairs` (Oct 2026, experimental)

Status: AI-generated, machine-tested (section 5 of `test_cases/scripts/check_charge_groups.py`), no human review, human
production testing pending. Off by default. Code: `GFNFF::buildIonPairChargeGroups`, `perceiveAromaticCationRings`, called at
the end of `calculateTopologyInfo` (`gfnff_method.cpp`); PARAMs `charge_groups`, `ionpair_charge`, `ionpair_aromatic` (`gfnff.h`).

`ensemble` and `reference` place only a net charge. In a neutral system of several fragments every fragment is held at 0, so
the cation and the anion of an ion pair are both neutral while the separated ions carry +1/-1. The parity rule sets the formal charges:
1. Q = 0, closed shell (spin 0, even electron count), at least two fragments in the final topology of the reference q-loop.
2. Cation fragments hold onium centres (`perceiveOniumCentres` valence_only: N with 4 / O with 3 partners in `nb_full`, no metal
   or over-coordinated partner); with `ionpair_aromatic true` also isolated 5/6 rings of C and N with a Hueckel count of 7.
3. Anions by parity: fragments without a centre whose neutral electron count is odd; an even fragment (host, solvent) stays 0.
   Their number must equal the number of centres n_plus, otherwise nothing happens (warning, "ambiguous").
4. One EEQ group per fragment (+`ionpair_charge` per centre, -`ionpair_charge` per anion, 0 otherwise); the topology is then
   rebuilt with these groups from pass 1, through `applyUserChargeGroups` as for a group file. Once per instance; one extra
   topology build at initialisation when the rule fires, no extra cost per step; the groups are constant, so the analytic gradient
   needs no extra term. Verbosity >= 1 prints the roles and that the rule is not meant for shared-proton reactions.

Measured (n = 1 run per structure; GMTKN55/S30L-CI references; held-out: r2SCAN-3c on xtb GFN2 geometries):
- IL16, 12 onium pairs: MAD 65.58 -> 21.06 kcal/mol (11 better, 1 worse); all 16 with `ionpair_aromatic`: 71.62 -> 19.11.
- S30L-CI (n = 30): MAD 10.445 -> 9.79 (24/27/28 affected; host even, so complex and monomer get the same roles).
- Held-out, 24 ion pairs: MAD 66.99 -> 12.08 (halides n = 6: 115.43 -> 26.45, others n = 18: 50.85 -> 7.29); 4 complexes ion
  pair + neutral fragment: 15.79 -> 7.32 (5.73 aromatic). The first record (67.27 -> 13.46, 16.3 -> 7.8) mixed in 15 stale geometries.

Limits (read before use):
- Not for reactions with a shared proton: proton-transfer transition states and H3O+/OH- or NH4+/NH2- type complexes get formal
  charges. WTMAD-2 92.595 -> 93.063 (+0.47; 92.955 with `ionpair_aromatic`); WCPT18 (n = 5) MAD 22.9 -> 102.9, PX13 (n = 2)
  22.9 -> 102.4, WATER27 H3OpH2O6OHm +40 -> +112; held-out shared-proton controls: 6 of 10 grow by +39 to +196 (MAD 32.25 -> 99.60).
  No tested quantity (contact distance, free EEQ charge, GFN-FF energy) separates these from ammonium halide pairs.
- Not perceived: guanidinium/amidinium, phosphonium, sulfonium, protonated carbonyls, fused or substituted aromatic cations,
  dianions (an anion fragment gets one -q), open shells. A fragment with k onium centres gets +k (not tested for k > 1).
- Groups frozen per instance: not for MD with bond breaking or `-gfnff.react` (not tested). GPU falls back to the CPU method
  (`hasUserChargeGroups`, read in the code, not run). Fragment order: IL16/144 changes by -2.4e-4 Eh between the two fragment
  orders without the switch and by -4.1e-4 Eh with it (230, 187, 008: 0; n = 1 run each).
- Falsifier set (24 pairs, 10 shared-proton controls, 4 three-fragment complexes): vault `_artefacts/2026-10-24_curcuma_ionpair-falsifier-package.tar.gz`;
  its `species/*.xyz` are NOT the geometries of its r2SCAN-3c energies for 85 of 124 species (rewritten after ORCA started): to be replaced.
