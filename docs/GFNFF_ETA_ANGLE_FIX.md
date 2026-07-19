# GFN-FF Angle `feta` η-Coordination Fix (Jul 2026)

Status: 🤖 AI-generated, ⚙️ machine-tested. Human production testing pending.

## Problem

The angle-bending metal correction `feta` (reduces the angle force constant for
metal-η-coordinated π-ligands) was firing for the wrong atoms. The old code proxied
η-coordination by the outer atoms' hybridization (`hyb==1 || hyb==2`, i.e. sp/sp²).
A σ-bonded carbonyl carbon is sp AND in a π-system but is **not** η-coordinated, so it
wrongly got `feta=0.3`, softening Rh–C(O) / Rh–Cl angles.

Reference (`external/gfnff/src/gfnff_ini.f90:1392-1394`): `feta=0.3` per outer atom fires
only when that outer atom is genuinely η-coordinated —
`imetal(center)==2 .and. itag(outer)==-1 .and. piadr(outer)>0`.

Symptom: ED39 = RhCl(CO)₂ angle energy `+0.00063 Eh` vs xtb 6.6.1 `--gfnff` `+0.00410 Eh`
(too low by ~3.5 mEh), the entire remaining ED39 residual.

## Fix

### 1. η-coordination detector (ported from Fortran)
- **`GFNFF::computeEtaCoordination(neighbor_lists)`** —
  `src/core/energy_calculators/ff_methods/gfnff_method.cpp:7392` (faithful port of
  `external/gfnff/src/gfnff_ini2.f90:170-198`). Sets `itag[i]=-1` iff atom i (ati≤10, in
  practice C) satisfies: Cp `(nbf≥4 & nbm==3)` OR alkene/alkyne `(nbf==3 & nbm==2)`, with
  the metal refinement: `nm==0`→not η; `nm==1`→η only if a non-metal neighbor of i is also
  bonded to that same metal (`ncm≥1`); `nm≥2`→η. `nbf` = full bonded adjacency
  (`topo.neighbor_lists`, includes metals); `nbm` = neighbors with `metal_type==0`.
- New field **`GFNFFTopology::itag`** (`gfnff.h:316`); declaration `gfnff.h:1762`.
- Populated once in **`calculateTopologyInfo()`** right after `neighbor_lists` is built
  (`gfnff_method.cpp:8972`). Also serialized/restored in the `.topo.json` cache
  (`topo_json["itag"]`, ~`gfnff_method.cpp:2591` / `2751`).

### 2. `feta` now gates on `itag==-1`
- **`getGFNFFAngleParameters()`** — `gfnff_method.cpp:5359` region. The old
  `is_pi = hyb==1||hyb==2` test is replaced by `is_eta(atom) = (itag[atom]==-1)`, keeping
  `feta*=0.3` per η outer atom (→0.3 one, 0.09 both), `feta=1.0` otherwise, still gated on
  `is_tm_center`.

**`piadr>0` nuance (honest note):** Fortran also checks `piadr(outer)>0`, but that is
redundant with `itag==-1` there (an η atom is by construction an alkene/alkyne/Cp carbon in
a π-system). Curcuma's `pi_fragments` is derived from the FULL neighbor list (including the
metal), so metal-side-on alkene carbons come out `pi_fragments==0`, whereas Fortran's
`piadr` — computed from the metal-REMOVED list for η atoms (`gfnff_ini2.f90:199`) — is >0.
Adding a `pi_fragments>0` guard was therefore tried and **reverted**: it wrongly dropped
every genuine η carbon back to `feta=1.0` (verified on ED32). Gating on `itag==-1` alone is
the faithful reproduction of the reference behavior. (`itag` is populated but **not** wired
into the metal-bond `btyp==6` path — that changes bond strength and needs its own
validation; left as M-X for now, comment at `gfnff_method.cpp:4884`.)

## Validation (built in `release/`, run against xtb 6.6.1)

### ED39 = RhCl(CO)₂ (primary target)
| quantity | before | after (native) | xtb 6.6.1 `--gfnff` |
|----------|--------|-----------------|---------------------|
| Rh-centered angle `feta` | 0.3 / 0.09 | **1.0** | (σ ligands, no η) |
| angle energy (Eh) | +0.00063 | **+0.0040977481** | +0.004097748029 |
| total energy (Eh) | −1.05686 | **−1.05338734** | −1.053386871618 |

Angle now matches xtb to <1e-9 Eh; total to ~5e-7 Eh. (before values from the task
diagnosis.)

### ED32 = COD–Rh complex (genuine η, must NOT regress)
`computeEtaCoordination` tags the 4 COD alkene carbons (15,16,17,21), each bonded to Rh.
`feta` for Rh-centered angles with an η carbon as outer atom = **0.3** (one η) or **0.09**
(both η) — same as the old hyb-based path for these genuinely η carbons, i.e. no regression
of the η mechanism. (ED32's overall total is off ~0.28 Eh, dominated by the known
unimplemented metal bond-term / the missing xtb `phi<60°` η-angle skip — separate from this
fix.)

### Neutral suites (no metals → `feta` never fires → must be unchanged)
`cd release && ctest -R gfnff`: 75/76 passed. gfnff (21) + validation (37) + gfnff_solvation
(28) all green; H2O/caffeine/C6H6 golden-value tests bit-identical. The single failure
`gfnff_gpu_vs_cpu_h2` is pre-existing and unrelated: it aborts with
`File not found: ../../external/gfnff/test/h2.xyz` (missing test-input path), and H2 has no
metal.

## Part 2: Metal-center angle equilibrium (phi0) fix

The correct `feta=1.0` (part 1) **exposed** a masked bug: octahedral metal carbonyls
(Cr(CO)₆ = PR01, Mn(CO)₅H = PR05) had a spurious +0.02..0.03 Eh angle energy that the old
(wrong) `feta=0.09` softening had been hiding. Root cause, in two places:

1. **Metal-center hybridization** was the geometric sp3 default (`hyb=3`) instead of the
   Fortran TM value. `external/gfnff/src/gfnff_ini2.f90:319-333` assigns a TM's hyb from its
   H-excluded coordination `nni`: `nni≤2`→1 (group≤−6→2), `nni==3`→2, `nni==4`→3,
   `nni==5 & group==−3`→3, **else (octahedral `nni≥6`) stays 0**. The Fortran analyzer
   confirms Cr(CO)₆ → sp-hybrid 0. Ported into `determineHybridization()`
   (`gfnff_method.cpp`, TM block after the `neighbor_count>=4` default), gated on
   `metal_type[z]==2`.
2. **Metal-center angle branch missing.** `gfnff_ini.f90:1601-1610`: when the angle center is
   a metal (`imetal>0`), r0 is set purely by the center's hyb, **overriding all
   element-specific rules** (Fortran's final word before `vangl`): `hyb0`→r0=90 & `f2=1.35`
   ("important difference to other bends, big effect"), `hyb1`→180, `hyb2`→120, `hyb3`→109.5,
   plus the linear GEODEP override `phi>linthr(160°)`→180. Ported into
   `getGFNFFAngleParameters()` just before the `equilibrium_angle` conversion, gated on
   `imetal_ctr>0` (with the Fortran Sn/Pb/Bi small-CN de-metalization, `gfnff_ini.f90:274`).

With hyb=0 for octahedral centers, the branch gives r0=90 → the 90° cis L-M-L angles sit at
equilibrium (zero energy) and the 180° trans angles flatten to r0=180 (zero energy).

### Validation (native vs xtb 6.6.1 `--gfnff`)
| structure | angle before | angle after (native) | angle xtb | Δtotal now |
|-----------|-------------:|---------------------:|----------:|-----------:|
| PR01 Cr(CO)₆ | +0.029 Eh | **+1e-10** | +8e-11 | bond-only |
| PR05 Mn(CO)₅H | +0.019 Eh | **+0.005129** | +0.005310 | bond-only |
| ED02 Fe(CO)₄ | — | (exact) | — | 4e-7 Eh |
| ED39 RhCl(CO)₂ | (part 1) | +0.0040977 | +0.0040977 | 5e-7 Eh |

The PR01/PR05 angle regression is fully resolved. Their remaining residual (Cr-C 20 mEh,
Mn-C 13 mEh) is **entirely the metal bond term** (separate, pre-existing — the metal
bond-stretch branch reduces the ~1.95× overbinding but a per-carbonyl residual remains; see
[GFNFF_METAL_BOND_ANALYSIS.md](GFNFF_METAL_BOND_ANALYSIS.md)), not the angle.

Neutral regression: `ctest -R gfnff` validation group 37/37 pass (both changes are TM-gated,
so metal-free systems are byte-identical).

## Files changed
- `src/core/energy_calculators/ff_methods/gfnff.h` — `itag` field + `computeEtaCoordination` decl
- `src/core/energy_calculators/ff_methods/gfnff_method.cpp` — detector impl, population in
  `calculateTopologyInfo`, JSON cache round-trip, `feta` gate rewrite, metal-bond comment,
  TM-center hybridization rule (`determineHybridization`), metal-center angle branch
  (`getGFNFFAngleParameters`)
