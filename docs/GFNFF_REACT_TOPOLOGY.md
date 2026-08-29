# GFN-FF React Topology Mode (`topology_mode = react`)

🤖 AI-generated, ⚙️ machine-tested. **Not human production tested.** This mode is
uncharted territory: a reactive extension of GFN-FF exists in no reference
implementation, so there is no external result to validate against. Read the
Approximations section before using it for anything beyond qualitative dynamics.

## What it does

The native GFN-FF freezes its bonded terms (bonds, angles, torsions, inversions)
at initialisation. Bond *breaking* is smooth on that frozen surface (the bond term
is a finite Gaussian well, `E = fc * exp(-alpha (r-r0)^2)` with `fc < 0`), but bond
*formation* is impossible: a new pair has no bond term, keeps the stiff non-bonded
repulsion parameters, and no angle/torsion/inversion terms ever appear for it.

`react` is a third topology mode next to `auto` (default, adaptive caching) and
`constant` (frozen). It makes the bond list dynamic:

1. An O(N^2) hysteresis scan re-detects bonds during MD, every
   `react_check_every` energy calls or when any atom moved more than
   `react_check_disp_bohr` since the last scan.
2. A pair becomes a bond when `r < react_bond_form_factor * thr` and an existing
   bond survives until `r > react_bond_break_factor * thr`, with
   `thr = (rcov_i + rcov_j) * fat_i * fat_j` (the same fat-scaled covalent
   threshold the initial detection uses at factor 1.3).
3. When the bond set changes, **everything bonded is regenerated** through the
   normal parameter-generation path: bonds, angles, torsions (incl. extra and
   sTorsions), inversions, the bonded/non-bonded repulsion partition, Coulomb
   parameters, BATM/ATM lists and HB/XB detection. Fragments and the EEQ
   fragment-charge constraints follow the new topology in the same step.
4. The regenerated interaction lists are swapped into the live engines
   (`FFWorkspace::rebuildInteractionLists()`, legacy `ForceField` sync, and on
   CUDA/ROCm a full device-workspace reconstruction).

The react bond set owns the internal forced-bond list, so every topology consumer
(adjacency, rings, Hueckel, repulsion partition, EEQ fragments) sees exactly the
same bonds; only the hysteresis scan can change them. A rebuild reproduces the
exact state a fresh initialisation with that bond list would produce
(machine-tested to bit-identity, see Validation).

## Hysteresis defaults: optimistic formation, conservative retention

The defaults are `react_bond_form_factor = 1.6` and `react_bond_break_factor = 2.6`.
The asymmetry is deliberate:

- **Formation is optimistic** (large radius) because a colliding pair must reach
  the formation radius against the non-bonded repulsion wall. Measured for H + H:
  at the earlier 1.25 factor the formation radius (0.83 A) carries ~180 kJ/mol of
  non-bonded repulsion and no formation event occurred in 3 ps at 1000 K; at 1.6
  (1.06 A, ~105 kJ/mol) free H atoms recombine to H2 at 3000 K.
  The upper limit is set by hydrogen bonding: the H...O formation radius at 1.6 is
  ~1.57 A, which leaves normal hydrogen bonds (>= 1.7 A) untouched. Raising the
  factor much further would turn strong H-bonds into covalent bonds.
- **Retention is conservative** (break late) because the Gaussian well decays
  slowly; removing a bond earlier deletes energy the pair still has. Measured for
  H2 at 2500-3000 K: breaking at factor 1.45 removed ~90% of an intact well
  (dE_jump +482 kJ/mol) and fired spuriously from ordinary vibrations; at 2.6 a
  break event costs +21 to +34 kJ/mol and only fires after the pair has genuinely
  climbed the well.

## Approximations and known limitations

- **Energy is not conserved across rebuild events.** Forming a bond drops the
  potential energy by roughly the well depth at the formation radius (measured
  H + H -> H2: about -450 kJ/mol, the recombination energy appearing in one step);
  breaking removes the residual well (+21 to +34 kJ/mol at the default factor).
  Every rebuild logs its measured jump as `dE_jump` at verbosity >= 1.
  **The mode is NVT-only** — a thermostat must absorb the jumps; NVE runs will
  drift at every event.
- Between events the gradient is exact for the current topology (same code path
  as a normal GFN-FF run), up to the pre-existing residual below.
- **Pre-existing analytic-vs-FD gradient residual (recorded, not caused by react
  mode):** for H2 the analytic gradient deviates from central finite differences
  by up to 5.3e-2 Eh/A at equilibrium and 1.0e-1 Eh/A at 1.0 A stretch. Term
  isolation shows the bond term alone carries 2.4e-2 (the dynamic-r0 CN chain)
  and most of the rest sits in the repulsion term (whose gradient is documented
  as partial). The residual is bit-identical between a react-formed surface and a
  fresh forced-bond initialisation, so react mode adds nothing to it. Tracked by
  `test_gfnff_react_fd`, which prints both numbers.
- GFN-FF is **not parameterized for transition states**. Barrier heights, TS
  geometries and reaction energetics from this mode are qualitative at best.
  Measured for the ammonia-synthesis target: N2 + 3 H2 at 3500 K in a 3.5 A wall
  forms N-H bonds (NH2 observed) and splits N2 within 20 ps — but in the
  event-churn regime described below, so this demonstrates the machinery, not
  Haber-Bosch energetics.
- 3-/4-body terms of a fresh bond switch on through the existing distance
  damping `1/(1+r~^4)`, i.e. semi-smoothly; the repulsion re-partition
  (non-bonded alpha -> bonded alpha) jumps at the event.
- **Formation valence cap (empirical):** a new bond is only formed while both
  partners are below a hard coordination cap — 2 for H/He (allows the linear
  exchange intermediate), 6 for everything else (matching the angle generator,
  which skips centres with more than 6 neighbours). Closest candidate pairs claim
  the remaining valence first. Without the cap a hot, confined system over-bonds
  into an unphysical cluster whose energy eventually turns NaN (observed for
  N2 + 3 H2 at 3500 K in a 3.5 A wall). The cap is empirical; the principled
  replacement is an over-coordination energy (see Open refinements).
- **Event churn acts as an energy pump at extreme conditions.** A formation
  releases the well depth (-450 kJ/mol for H+H) and a later break costs only the
  decayed tail (+30 kJ/mol), so each form/break cycle injects net potential
  energy that the thermostat must drain. In a hot, dense system with frequent
  events (measured: ~300 events / 20 ps for N2 + 3 H2 at 3500 K in a 3.5 A wall)
  this drives artificial churn. At such conditions results are machinery
  demonstrations, not thermodynamics.
- Angles at a centre with more than 6 neighbours are skipped by the generator
  (pre-existing rule); transiently hypervalent atoms during an exchange are
  therefore under-described.
- A rebuild is O(N^3) (full topology incl. Floyd-Warshall). Fine for small
  systems (low ms per event); a warning is printed above 2000 atoms.
- Not tested: ALPB solvation + react, charged/radical species beyond H atoms,
  large systems, long trajectories (> 50 ps).

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `topology_mode` | `auto` | `auto` (alias `default`), `constant`, `react` |
| `react_bond_form_factor` | 1.6 | formation radius factor (optimistic) |
| `react_bond_break_factor` | 2.6 | retention radius factor (conservative) |
| `react_check_every` | 5 | scan every N energy calls (0 = displacement only) |
| `react_check_disp_bohr` | 0.25 | scan when any atom moved this far (0 = off) |

CLI: `curcuma -md in.xyz -method gfnff -gfnff.topology_mode react ...`

Guards: `react` + `static_charges`/`static_cn`/`static_all` (incl. `gfnff-fast`)
aborts initialisation — frozen CN/charges cannot follow a changing topology.
`Mol`-provided bonds (e.g. polymerbuild) seed the initial react bond set and are
owned by the scan afterwards. SimpleMD's `topo_check` / `epot_abort` must stay
off (their defaults) for reactive runs. Unknown mode strings warn and fall back
to `auto`.

## GPU (CUDA / ROCm)

On a rebuild the device workspace is destroyed and reconstructed from the
refreshed parameter set (`initGPUWorkspace()` re-run), and the EEQ solver,
fragment topology and CN pair list are reset. Cost per event: device
reallocation (~ms) plus one leaked ~100 KB parameter set (a pre-existing CUDA
heap workaround in the workspace constructor), bounded by the number of events.

## Validation (machine-tested)

- `gfnff_react_fd_gradient` (ctest): a formation rebuild reproduces the energy
  and analytic gradient of a fresh initialisation with the same forced-bond
  topology to bit-identity; records the FD residual on both surfaces.
- `cli_simplemd_13_gfnff_react_noevent`: NH3 at 300 K, react vs auto with the
  same seed — zero bond events and a bit-identical trajectory (react mode does
  not perturb non-reactive dynamics; also covers the NH3 inversion terms).
- `cli_simplemd_14_gfnff_react_h_recombination`: 4 free H atoms at 3000 K in a
  2.5 A spherical wall recombine; every rebuild logs `dE_jump`; run stays stable.
- Manual: 2x H2 at 2500 K, 3 ps — no spurious events at the default factors;
  17/17 `gfnff_val_*` Fortran-parity references and the full `cli_simplemd_*`
  suite unchanged.

## Open refinements

- Per-pair energy-based break criterion (`alpha (r-r0)^2 > threshold`) instead of
  a global geometric factor — removes each bond exactly where its own well has
  decayed.
- Replace the hard formation valence cap by an **over-coordination energy**: a
  penalty that rises with the coordination number beyond the element valence
  (ReaxFF-style). Its offset/shape could be parameterized against QM
  calculations of deliberately hyper-coordinated species — turning the empirical
  cap into a learned, quantum-mechanically anchored term.
- **Quality criterion per supported reaction:** benchmark react-mode reaction
  trajectories against NEB paths at GFN2-xTB or a QM level. Current target:
  ammonia synthesis (N2 + 3 H2); candidates for later: primordial-soup
  (Miller-Urey-type) chemistry.
- Push the force field's live bond list into the qurcuma viewer frames so drawn
  and simulated topology agree exactly during reactive runs (the viewer's own
  hysteresis redraws bonds independently at 1.25/1.45).
- ROCm/CUDA rebuild long-run stress test; ASAN pass over repeated
  `generateGFNFFParameterSet()` calls on a live object.
