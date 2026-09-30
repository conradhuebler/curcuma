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

## Caveats

- **The hardcoded fallback is the real default.** `GFNFF::GFNFF(const json&)` reads each setting with
  `m_parameters.value(key, FALLBACK)`, and the CLI does not merge full registry defaults into
  `controller["gfnff"]`. Changing the PARAM default in `gfnff.h` alone has no runtime effect.
- **GPU**: the GPU wrapper does not call `GFNFF::Calculation()`, so `-gpu` runs keep the reference
  rule (code inspection; not runtime-tested).
- Cost: each variant is a full GFN-FF evaluation; tiny for the parity-decided cases, 2^k corners
  times placements with the window on.
- A reused calculator (optimisation, MD) on master still has the stale CN / D4-C6 defect of
  energy-only calls (fixed on `reactff2-llm`, not part of this port); the ensemble inherits it.
- Tested: GMTKN55, MOR41, S30L-CI, the Cl2- / H3O+(H2O)2 cases of `cli_gfnff_05_frag_charge_ensemble`
  and `gfnff_frag_charge_history`. **Not tested**: periodic systems, charged systems with more than
  a few fragments in a window, long MD with the window on.
