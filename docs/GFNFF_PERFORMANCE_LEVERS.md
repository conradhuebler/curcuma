# GFN-FF + EEQ Performance Levers (June 2026)

**Question:** Where does curcuma's native GFN-FF single-point / MD time actually go on
large condensed-phase systems (thousands of atoms, many fragments), and which concrete
changes give the biggest speedups?

**Bottom line:** The earlier "~90 s on mixture2" figure was a `-v 3` artifact — verbose
logging (per-pi-system prints, per-pair C6 debug, progress bars) dominated. On the
**current** build a single point is **~6 s for 3007 atoms and ~19-21 s for 6200 atoms**
(`-v 1`, 4 threads). The cost is **not** the EEQ solve and **not** D4 dispersion as folklore
held — it is **host-side parameter generation**, dominated by the **hydrogen-bond candidate
list** (~33 % of wall) and the **energy/term distribution + evaluation** (~40 %). The EEQ
solve is already well-threaded; D4 already has a 60-Bohr cutoff on the CPU path.
**The single biggest win is a tighter, distance-pruned HB + non-bonded pair generation.**
GPU (CUDA, GTX 1660) is **not** a single-point win here — it is *slower* (22.8 s vs 19.4 s)
because the heavy host-side generation (incl. an uncapped 13.7 M-pair JSON D4 list) still
runs on the CPU and the device pipeline only pays off when amortized over many MD steps.

Measurement scheme: wall = `date` around `curcuma -sp ... -method gfnff -no_bmt`; per-phase
deltas = line-buffered stdout timestamped (`scratchpad/timing/stamp.py`); ablations toggle
PARAMs. All on this machine (8 logical cores, OpenBLAS-OpenMP, GTX 1660).

---

## Measured wall times (single point, `-v 1`)

| System | atoms | frags | T=1 | T=4 | T=8 | GPU cuda |
|---|---|---|---|---|---|---|
| mixed_3007 | 3007 | 654 | 9.13 s | 6.04 s | 5.77 s | - |
| mixture2   | 6200 | 1400 | - | 21.3 s | 19.4 s | 22.8 s |

- **Threading scales poorly: ~1.5x at 4 threads, saturates by 4.** The serial host-side
  generation (HB structs, JSON, single-threaded sections) is the ceiling, not the threaded
  energy kernels.
- **O(N^2) dominated:** 3007 -> 6200 atoms (2.06x) costs 6.0 -> 21.3 s (3.5x), close to 2.06^2.

## Per-phase breakdown — mixture2 (6200 atoms, 4 threads)

Largest wall deltas between consecutive log markers:

| Phase (log marker) | delta | what it is |
|---|---|---|
| "Force field topology summary" | 4.35 s | term distribution to threads + energy/gradient eval (4.5M Coulomb + 4.5M nonbonded-rep + 13.7M/native-disp pairs) |
| "Detected 3357808 hydrogen bonds" | 4.27 s | HB candidate list: **3.36M** `GFNFFHydrogenBond` structs, each with 3 `std::vector` members |
| "Molecule set: 6200 atoms" | 3.41 s | final energy accumulation + output |
| "Found 3600 HB hydrogens" | 2.23 s | HB donor/acceptor pre-pass (per-H acceptor scan) |
| "Two-phase EEQ charge calculation" | 1.46 s | EEQ Phase-2 (distance matrix 6200^2 + Coulomb matrix 6200^2 + SchurCholesky) |
| "D4: Loaded reference CN data" | 1.32 s | one-time file load of dftd3param CN data (cacheable across runs) |
| torsion generation | 0.99 s | |

**HB generation alone (pre-pass + detection) = ~6.5 s = ~33 % of wall.** Energy
distribution + eval + output = ~7.8 s = ~40 %. EEQ = ~1.5 s = ~8 %.

## Ablations (mixed_3007, 4 threads)

| Config | wall | E (Eh) | note |
|---|---|---|---|
| baseline | 5.99 s | -422.55322598 | |
| `-dispersion_cutoff_bohr 15` | 5.67 s | -422.43929237 | only -0.3 s; **changes E by 0.11 Eh** (help's "<1 muEh" is for sparse systems, not dense liquids) |
| `-skip_phase2 true` | 6.09 s | -422.55322598 | no gain (EEQ Phase-2 not the bottleneck) |
| `-hb_thr1/2 100/200` (mixture2) | 15.9 s (vs 19.4) | -856.4 (vs -857.18) | -3.5 s but **changes E by 0.79 Eh** — too aggressive as-is |

The ablations confirm: D4 and EEQ are *not* where the time is. HB list size is, but the
existing knob trades accuracy too coarsely.

---

## Ranked levers

| # | Lever | file:line | est. impact | effort |
|---|---|---|---|---|
| 1 | **Prune HB candidate generation by acidity/basicity AND a per-H acceptor distance, and stop materializing 3.36M heavy structs.** The `ab_pairs` cell-list pass keeps every O/N within 15.8 Bohr; for dense liquids that is ~N. Each survivor spawns a `GFNFFHydrogenBond` with 3 `std::vector` members. Pre-reserve, use flat `int` arrays for neighbor lists, and drop pairs whose `strength` is below a real (not 1e-6) threshold. | `gfnff_method.cpp:8008-8059` (ab_pairs build), `8062-8105` (struct create), `8287` (count) | ~3-5 s on mixture2 (~25 %); the single biggest SP win | M |
| 2 | **Give the GPU path the capped native D4 list, not the uncapped JSON one.** GPU calls `GenerateParameters` -> 13.7M JSON pair objects (12 string inserts each) on the host; the CPU path uses `GenerateDispersionPairsNative` with the 60-Bohr cutoff + structs. Route the GPU through the native generator (or upload geometry + cutoff and build pairs on-device). | host JSON gen `d4param_generator.cpp:255,423-564`; consumed by GPU via `gfnff_method.cpp:9764`; native (capped) path `d4param_generator.cpp:1789-1897` used at `:9575` | makes GPU SP competitive; removes ~2-4 s host JSON churn; prereq for GPU SP wins | M |
| 3 | **Cache `getChargeWeightedC6` per element-pair x CN-bin instead of per atom-pair.** For ~4 elements (H/C/N/O) the 7x7 contraction (4.3M pairs x up to 49 FMAs) is massively redundant — C6 depends only on `(Zi,Zj,CN_i,CN_j,q)`. Apply the existing `buildAtomRefW`/`contractC6Gfn2` AP2 split (already used by `D4Evaluator`) to the GFN-FF native loop too: build per-atom 7-vector once, contract in the pair loop. | `d4param_generator.cpp:1851` (per-pair call inside O(N^2) loop) vs hoistable `:1305+` `buildAtomRefW`/`contractC6Gfn2` | ~0.5-1 s; folds into lever 2 | M |
| 4 | **Make `dispersion_cutoff_bohr` actually accuracy-safe for dense systems** (or auto-tune). 15 Bohr drifts 0.1 Eh on the mixture; the safe value is system-dependent. Add a smooth taper near the cutoff or a density-aware default so the cutoff can be turned on without an energy jump. Then default it ON for N>~2000. | cutoff applied `d4param_generator.cpp:1832-1843` (native), absent in JSON `:433`; PARAM `gfnff.h:252` | ~0.5-1 s once safe | M |
| 5 | **Distance-prune the non-bonded repulsion list.** 4.5M nonbonded repulsions for 3007 atoms = uncapped O(N^2) (`rep_r_cut=50` is huge). A real cutoff (repulsion decays exp) culls most pairs. | `gfnff_method.cpp` nonbonded-rep generation (search "non-bonded repulsions"); PARAM `rep_r_cut` | ~0.5-1 s | S-M |
| 6 | **Cache the D4 reference-CN file load across the process / ship it as a static table.** "Loaded reference CN data from Fortran dftd3param.f90" costs ~1.3 s of file parsing every run; it is geometry-independent. | `d4param_generator.cpp` reference-CN init (`d4_reference_cn_fortran.cpp` include + load) | ~1.3 s one-time per SP; large relative win for small MD warmups | S |
| 7 | **Reduce the EEQ double matrix build: distance matrix AND Coulomb matrix are each built N^2 with progress bars.** They can share one pass; the Coulomb matrix can be built directly from the distance matrix without a second O(N^2) traversal. (Solve itself is fine.) | `eeq_solver.cpp` Phase-2 distance/Coulomb matrix build (search "Phase 2: distance matrix" / "Coulomb matrix") | ~0.3-0.5 s | S |
| 8 | **Thread the serial generation sections (HB struct create already parallel; torsion/inversion/topology-summary are not).** torsions 0.99 s + topology summary setup are largely single-threaded; this is why T=4 only gives 1.5x. | `gfnff_method.cpp` torsion/inversion generation, term distribution | improves the threading ceiling (toward 3-4x) | M-L |

S = hours, M = day, L = multi-day.

---

## CPU vs GPU note

- **GTX 1660 GPU GFN-FF is compiled in** (`USE_CUDA_GFNFF=ON`, verified at runtime:
  "gfnff-gpu: GPU workspace ready ... GPU EEQ enabled") and **numerically correct**
  (mixture2 SP energy -857.18483720 vs CPU -857.18483727, diff 7e-8 Eh).
- **For a single point the GPU is slower** (22.8 s vs 19.4 s) because the host-side
  generation (HB list, topology, the **uncapped 13.7M-pair JSON D4 list** at
  `gfnff_method.cpp:9764`) is unchanged and dominates; the device pipeline's one-shot
  upload/setup is pure overhead with nothing to amortize it.
- **GPU is an MD/opt lever, not an SP lever.** Per CLAUDE.md the ROCm GFN-FF MD path is
  ~2.3x faster than CPU on a 1410-atom polymer after the gather-kernel rewrite. The CUDA
  GTX 1660 path will benefit similarly only once (a) lever 2 removes the host JSON churn and
  (b) the per-step host work (HB rebuild) is throttled. For a *single point* on this
  hardware, stay on CPU with `-threads 4`.

## Biggest single win

> **IMPLEMENTED (Jun 2026):** Lever 1 done in `gfnff_method.cpp` (cell-list HB candidate
> generation: nhb2 via per-atom bonded-hydrogen lists `hyd_on[]`, nhb1 via
> `cell_list.forEachNeighbor(i, hbthr2)` instead of the full per-pair hydrogen scan).
> EXACT — identical HB set, energies bit-identical (mixed_3007/water_1002/urea_1000/mixture2
> all match to the last digit). Measured ~1.5x on mixture2 SP (clean A/B same machine/load:
> ~35.7 s -> ~23.2 s @4thr). Gated to the large-system cell-list path; small systems
> (validation suite) keep the unchanged path. Lever 2 (GPU native D4 list) remains open.

**Lever 1 (HB candidate generation).** It is ~33 % of single-point wall on a dense
H-bonding liquid, it is the same cost every MD step (HB list is rebuilt on geometry
change), and it is pure host-side bookkeeping with no accuracy trade-off if done right
(prune by real acidity/basicity/distance + flat neighbor arrays + reserve). Fixing it
helps CPU SP, CPU MD, and unblocks the GPU MD path simultaneously.

---

## Verification commands

```bash
# Per-phase timing (line-stamped):
OMP_NUM_THREADS=4 stdbuf -oL -eL ./release/curcuma -sp mixture2.xyz -method gfnff \
  -threads 4 -v 2 -no_bmt 2>&1 | sed -ur 's/\x1b\[[0-9;]*m//g' | python3 stamp.py

# Thread scaling:
for T in 1 4 8; do time OMP_NUM_THREADS=$T ./release/curcuma -sp \
  test_cases/eeq_mixture_fractions/mixed_3007.xyz -method gfnff -threads $T -v 1 -no_bmt; done

# GPU (correct but not an SP win on GTX 1660):
./release/curcuma -sp mixture2.xyz -method gfnff -gpu cuda -threads 4 -v 1 -no_bmt
```

> AI-generated performance audit, June 2026. Measurements reproducible on the machine
> above; absolute times are hardware-specific. No code changed — analysis only.
