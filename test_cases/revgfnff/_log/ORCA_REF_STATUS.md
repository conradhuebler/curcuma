# ORCA reference campaign -- session 2026-09-13 (branch reactff2-llm)

Standalone status for the three deliverables of this session. Driver `scripts/revgfnff_ref.py`,
charge report `scripts/revgfnff_hirshfeld.py`, detailed charge table `HIRSHFELD_CHARGES.md`.
Level unchanged: r2SCAN-3c (ORCA 6.1), `TightSCF`, EnGrad single points. Budget 16 cores
(`--jobs 4 --nprocs 4`, or `--jobs 2 --nprocs 4` twice when two campaigns ran side by side).
Prior campaign state: `WP2_STATUS.md`. Nothing under `src/` touched, no commits made.

Cost measured before choosing the item-1 scope: one EnGrad point **with** Hirshfeld costs
**9-10 s** at 4 cores (H2/HF/NH2OH; the class-A campaign's own average is 23 s/point over 1100
points, 7.08 h ORCA wall). Decision: **the full set** -- all 63 class-H series, 1260 points --
not far-points-only, because the per-bond charge change needs both ends of the stretch and the
marginal cost of the full set is ~1.5 h elapsed (measured: 48 min).

## Item 1 -- Hirshfeld charges at the class-A points (class H) -- COMPLETE

| | |
|---|---|
| series | 63 (32 class-A curves: 31 RKS + 32 UKS), 1086/1260 points ok |
| series fully converged | 53/63; all 10 incomplete are UKS broken-symmetry series (SCF, see below) |
| ORCA return codes | 0 for all 63 (failures are in-SCF, not crashes) |
| ORCA wall | 3.13 h (48 min elapsed at 4x4) |
| charges extracted | every converged point has a Hirshfeld list; 0 finite energies without charges |

Approach: class H re-runs the **same** points as class A (same cached reference geometries, same
`CURVE_GRID`) with `%output Print[P_Hirshfeld] 1`, into its own tree `ref/H/<system>/`, so the
existing class-A results stay valid and untouched. Verified before launch and after: for
`h2_H-H_rks` and `hf_H-F_uks` the class-H coordinates and energies are **bit-identical** to class
A's (max difference 0.0 in both), so point k of `ref/H/<sys>` is point k of `ref/A/<sys>`, same
atom order. Per-atom and per-fragment (side-of-the-stretched-bond) charges at r_eq and at 3.5
r_eq: `HIRSHFELD_CHARGES.md`.

### Headline: charge change on the two bonded atoms from r_eq to 3.5 r_eq (RKS series)

Mean over the series of each element pair, in e. The polar heavy-heavy bonds the bond-term work
is about to touch are marked `*`.

| pair | n | mean dq(atom i) | mean dq(atom j) | max abs dq | largest |
|---|---|---|---|---|---|
| C-O * | 3 | +0.034 | -0.064 | 0.151 | co_CTO_rks (C +0.087 -> +0.238, O -0.087 -> -0.238) |
| C-N * | 3 | -0.024 | +0.056 | 0.087 | hcn_CTN_rks |
| N-O * | 1 | +0.044 | -0.041 | 0.044 | nh2oh_N-O_rks |
| O-O * | 1 | +0.012 | +0.012 | 0.012 | h2o2_O-O_rks |
| C=O | 1 | -0.065 | +0.016 | 0.065 | h2co_CDO_rks (C +0.142 -> +0.078, O -0.226 -> -0.210) |
| C=N | 1 | -0.054 | +0.078 | 0.078 | ch2nh_CDN_rks |
| N-N | 3 | +0.009 | +0.009 | 0.015 | n2h2_NDN_rks |
| C-C | 3 | +0.015 | +0.015 | 0.077 | c2h2_CTC_rks |
| C-H | 2 | -0.013 | -0.013 | 0.084 | hcn_HC-H_rks |
| H-O | 2 | +0.034 | -0.047 | 0.087 | ch3oh_HO-H_rks |
| H-N | 1 | +0.079 | -0.072 | 0.079 | nh3_N-H_rks |
| H-F/H-Cl | 1+1 | +0.070/+0.051 | -0.070/-0.051 | 0.070 | hf_H-F_rks |
| C-F | 1 | +0.070 | -0.187 | 0.187 | ch3f_C-F_rks |
| C-Cl | 1 | +0.074 | -0.131 | 0.131 | ch3cl_C-Cl_rks |
| Cl-N | 1 | -0.003 | -0.208 | 0.208 | ncl3_N-Cl_rks |
| F-N | 1 | +0.024 | -0.206 | 0.206 | nf3_N-F_rks |
| F-O | 1 | -0.012 | -0.107 | 0.107 | of2_O-F_rks |
| Cl-F | 1 | +0.039 | -0.038 | 0.039 | clf_F-Cl_rks |
| Cl-O | 1 | +0.040 | -0.033 | 0.040 | hocl_O-Cl_rks |
| H-H, F-F, Cl-Cl, N#N | 1+1+1+1 | <= +0.001 | <= +0.001 | 0.001 | homonuclear |

What this measures is the **reference** charge change from r_eq to 3.5 r_eq, not the EEQ
comparison: on the four polar heavy-heavy bonds that matter here (C-O, C-N, N-O, O-O) the
closed-shell reference moves by **<= 0.09 e** over the whole stretch, with CO as the one clear
exception (0.15 e, and its C#O charges double: +0.087/-0.087 -> +0.238/-0.238). Labeled
estimate of what that is worth: for one charge pair at 4.2 A, a 0.08 e change on each atom moves
the pair Coulomb term by ~2.8 kcal/mol (332 q1q2/r, r in A) -- an order of magnitude below the
-20..-35 kcal/mol drift quoted for these bonds, so the reference charges alone do not supply it.
Whether the drift is then the EEQ side (chi(CN)) or something else in the bond/repulsion terms is
the falsifier's question and needs curcuma's EEQ at the **kept** topology (at 3.5 r_eq the GFN-FF
topology has already dropped the bond; `revgfnff_curves.py`'s `-batch_reuse_topology true` is the
right mode, and `-verbosity 4` prints per-atom `EEQ_PHASE1_CHARGES`). Not computed here.

The UKS broken-symmetry series is a different electronic state (diradical near r_eq, physical
homolysis limit at the far point) and is kept separate in `HIRSHFELD_CHARGES.md`; e.g. H-F goes
-0.266/+0.266 -> -0.003/+0.003, which is the radical limit, not a charge-model statement.

### Failures (item 1)
10 of the 32 UKS series did not converge (all `q=0`): ch3cl(0/20), co(0/20), of2(0/20),
hocl(1/20), h2co(1/20), h2o2(1/20), clf(1/20), ncl3(2/20), cl2(4/20), o2(16/20). ORCA error:
`ORCA_LEANSCF ... the SCF has not converged` -- a real convergence outcome, not a driver bug
(`rc=0`, ORCA aborts the job chain so the later points are never computed). Not retried with
`--slowconv`/`--uks-inside-out`: all 31 RKS series (the closed-shell charge reference) are
complete, and the UKS-BS solution near r_eq is a broken-symmetry artefact anyway. (The retry
described under item 3 covers the class-A of2/clf series, not these class-H ones.) Recorded, not
hidden.

## Item 2 -- rigid contact scans (class S) -- COMPLETE

New class in the driver (`CONTACTS`, `contact_dimer`), because an intermolecular separation scan
is not a bond stretch and must not go through `stretched()`. Fragment A's contact atom sits at
the origin with A's contact direction along +z, B's contact atom at (0,0,d) with B's direction
along -z; both fragments are rigid, so the scan is exact. Scan variable is the **heavy-atom**
separation (anchoring on the hydrogen instead would put the two heavy atoms 0.96 A apart at
d = 2.0).

| system | contact | pts | n_ok | rc | wall_s | d_min (A) | E_int at d_min (kcal/mol) |
|---|---|---|---|---|---|---|---|
| water_dimer_OO | O...H-O | 20 | 20/22* | 0 | 203 | 3.00 | -3.84 |
| hf_dimer_FF | F...H-F | 20 | 20/22* | 0 | 198 | 2.90 | -2.74 |
| nh3_h2o_NO | O...H-N | 20 | 20/22* | 0 | 198 | 3.25 | -1.92 |
| ch4_h2o_CO | O...H-C | 20 | 20/22* | 0 | 198 | 3.80 | -0.51 |

\* the 20 scan points all converged; the two extra points per job (the isolated monomers) were
**rejected by ORCA**: a `$new_job` chain cannot change the atom count -- GUESS fails with
"Input geometry does not match current geometry" when the monomer (3 or 5 atoms) follows the
dimer (6 or 8 atoms), since the MOs are carried over. The driver now runs the 20 scan points
only (`--classes S` re-plans as 20 pts and reports all four "done"); interaction energies are
referenced to the d = 6.0 A point, where the interaction is < 0.05 kcal/mol for these four.
The two dirs written before that fix still carry the two null monomer entries (harmless:
consumers skip null points).

The minima sit where they should for these dimers and the strength ordering is right: water dimer
-3.84 kcal/mol at O...O 3.00 A, HF dimer -2.74 at F...F 2.90, NH3...H2O -1.92 at N...O 3.25,
CH4...H2O -0.51 (weak C-H...O). These are **rigid, unrelaxed** contact curves with a linear
hydrogen bond, so they are not comparable to a relaxed-dimer literature binding energy -- the
falsifier compares them point by point against gfnff's curve at the same geometries, not against
experiment. Grid spacing near the minima is 0.1 A, so d_min and E_min carry that resolution.
Consumer note: a
class-S system is a dissociation-like curve in reverse -- the topology-defining frame is the
**largest** separation (the separated geometry), i.e. class S belongs with class C in
`_topology_index` (`revgfnff_fit.py`), which currently returns the smallest scan value for any
class other than C/D. One-line change; not made here to avoid touching a file another agent may
be editing.

## Item 3 -- the four missing class-A pairs (N-F, N-Cl, O-F, F-Cl) -- COMPLETE

Added to the driver's species registry and `CURVES`: NF3, NCl3, OF2, ClF, so they run the same
rigid-stretch grid with RKS + UKS like every other class-A curve. Reference geometries built by
the standard route (hand geometry -> curcuma GFN2 Opt -> r2SCAN-3c Opt): N-F 1.3879 A,
N-Cl 1.8029 A, O-F 1.4091 A, F-Cl 1.6556 A (symmetric in NF3/NCl3/OF2).

| pair | system | r_eq (A) | D_e (kcal/mol) | Morse D | Morse a (1/A) | Morse r_e | x50 | x90 |
|---|---|---|---|---|---|---|---|---|
| N-F | nf3_N-F_rks | 1.3879 | 92.85 | 79.60 | 1.95 | 1.3980 | 1.503 | 3.698 |
| N-Cl | ncl3_N-Cl_rks | 1.8029 | 57.93 | 52.20 | 1.70 | 1.8129 | 1.476 | 3.557 |
| O-F | of2_O-F_rks | 1.4091 | 92.87 | 90.41 | 1.90 | 1.4191 | 1.277 | 3.183 |
| F-Cl | clf_F-Cl_rks | 1.6556 | 107.44 | 104.26 | 1.75 | 1.6656 | 1.272 | 3.100 |
| N-F | nf3_N-F_uks | 1.3879 | 56.17 | 56.12 | 2.35 | 1.3779 | | |
| N-Cl | ncl3_N-Cl_uks | 1.8029 | 32.79 | 35.05 | 2.05 | 1.7929 | | |

(UKS D_e is quoted at the far point, 3.5 r_eq. The RKS rows above give the closed-shell curve's
own depth; x50/x90 are quoted for RKS only because they describe the bond term's tail.)

Conventions match `revgfnff_curves.py`: `r_eq` is the grid point with the lowest energy, `D_e` is
E(largest r) - E(min) of that series, the Morse is a least-squares fit over r <= 2.2 r_eq, and
x50/x90 are the dimensionless distances where the curve reaches 50 % / 90 % of D_e above the
minimum (1.228 / 2.970 for an ideal Morse -- the tail shape). **RKS D_e is not a physical bond
energy** (the closed-shell curve goes to the ionic limit, not the radical limit): the UKS series
is the one comparable to a bond dissociation energy, and its two values (N-F 56.2, N-Cl 32.8
kcal/mol) are of the right order for nitrogen-halogen single bonds -- quoted as a plausibility
check, no literature value was looked up for this report. The same code was run in the same pass
on three pairs that are already in the table, to show the new rows were produced identically
(RKS: C-H D_e 160.3 / a 1.65, C-F 147.7 / 1.75, O-H 170.8 / 2.00); no cross-comparison of those
values with the existing table was attempted here.

Failures: the UKS series of **of2 (1/20) and clf (0/20)** first failed outright -- BS-UKS SCF
failure, `rc=0`. Retried with the driver's documented strategy
(`run --classes A --only of2 clf --slowconv --uks-inside-out`, 561 s): **both recovered to
16/20**. The 4 points that still fail are the farthest ones in both chain directions
(of2 missing r >= 3.52 A, clf r >= 4.14 A), i.e. the dissociated broken-symmetry state is what
does not converge for these two; no UKS D_e is quoted for them, the RKS rows are complete.
nf3 and ncl3 converged 20/20 in both series.

## Per-series status

`pts` = points built into the series, `n_ok` = points with a converged energy and gradient; an INCOMPLETE series means ORCA stopped the job chain at the point where the SCF failed, so the remaining points were never computed. All 63 class-H series: `rc = 0`, charge 0, mult 1 (o2 mult 3, its O=O curve only).

### class H (Hirshfeld, item 1)

| system | q | mult | pts | n_ok | rc | wall_s |
|---|---|---|---|---|---|---|
| c2h2_CTC_rks | 0 | 1 | 20 | 20 | 0 | 192 |
| c2h2_CTC_uks | 0 | 1 | 20 | 20 | 0 | 228 |
| c2h4_CDC_rks | 0 | 1 | 20 | 20 | 0 | 193 |
| c2h4_CDC_uks | 0 | 1 | 20 | 20 | 0 | 219 |
| c2h6_C-C_rks | 0 | 1 | 20 | 20 | 0 | 188 |
| c2h6_C-C_uks | 0 | 1 | 20 | 20 | 0 | 215 |
| ch2nh_CDN_rks | 0 | 1 | 20 | 20 | 0 | 189 |
| ch2nh_CDN_uks | 0 | 1 | 20 | 20 | 0 | 209 |
| ch3cl_C-Cl_rks | 0 | 1 | 20 | 20 | 0 | 196 |
| ch3cl_C-Cl_uks | 0 | 1 | 20 | 0 **INCOMPLETE** | 0 | 20 |
| ch3f_C-F_rks | 0 | 1 | 20 | 20 | 0 | 188 |
| ch3f_C-F_uks | 0 | 1 | 20 | 20 | 0 | 217 |
| ch3nh2_C-N_rks | 0 | 1 | 20 | 20 | 0 | 194 |
| ch3nh2_C-N_uks | 0 | 1 | 20 | 20 | 0 | 210 |
| ch3oh_C-O_rks | 0 | 1 | 20 | 20 | 0 | 192 |
| ch3oh_C-O_uks | 0 | 1 | 20 | 20 | 0 | 200 |
| ch3oh_HO-H_rks | 0 | 1 | 20 | 20 | 0 | 195 |
| ch3oh_HO-H_uks | 0 | 1 | 20 | 20 | 0 | 207 |
| ch4_C-H_rks | 0 | 1 | 20 | 20 | 0 | 183 |
| ch4_C-H_uks | 0 | 1 | 20 | 20 | 0 | 189 |
| cl2_Cl-Cl_rks | 0 | 1 | 20 | 20 | 0 | 190 |
| cl2_Cl-Cl_uks | 0 | 1 | 20 | 4 **INCOMPLETE** | 0 | 68 |
| clf_F-Cl_rks | 0 | 1 | 20 | 20 | 0 | 195 |
| clf_F-Cl_uks | 0 | 1 | 20 | 1 **INCOMPLETE** | 0 | 34 |
| co_CTO_rks | 0 | 1 | 20 | 20 | 0 | 192 |
| co_CTO_uks | 0 | 1 | 20 | 0 **INCOMPLETE** | 0 | 18 |
| f2_F-F_rks | 0 | 1 | 20 | 20 | 0 | 192 |
| f2_F-F_uks | 0 | 1 | 20 | 20 | 0 | 224 |
| h2_H-H_rks | 0 | 1 | 20 | 20 | 0 | 200 |
| h2_H-H_uks | 0 | 1 | 20 | 20 | 0 | 197 |
| h2co_CDO_rks | 0 | 1 | 20 | 20 | 0 | 189 |
| h2co_CDO_uks | 0 | 1 | 20 | 1 **INCOMPLETE** | 0 | 30 |
| h2o2_O-O_rks | 0 | 1 | 20 | 20 | 0 | 185 |
| h2o2_O-O_uks | 0 | 1 | 20 | 1 **INCOMPLETE** | 0 | 27 |
| h2o_O-H_rks | 0 | 1 | 20 | 20 | 0 | 191 |
| h2o_O-H_uks | 0 | 1 | 20 | 20 | 0 | 202 |
| hcl_H-Cl_rks | 0 | 1 | 20 | 20 | 0 | 189 |
| hcl_H-Cl_uks | 0 | 1 | 20 | 20 | 0 | 196 |
| hcn_CTN_rks | 0 | 1 | 20 | 20 | 0 | 214 |
| hcn_CTN_uks | 0 | 1 | 20 | 20 | 0 | 234 |
| hcn_HC-H_rks | 0 | 1 | 20 | 20 | 0 | 209 |
| hcn_HC-H_uks | 0 | 1 | 20 | 20 | 0 | 198 |
| hf_H-F_rks | 0 | 1 | 20 | 20 | 0 | 201 |
| hf_H-F_uks | 0 | 1 | 20 | 20 | 0 | 203 |
| hocl_O-Cl_rks | 0 | 1 | 20 | 20 | 0 | 216 |
| hocl_O-Cl_uks | 0 | 1 | 20 | 1 **INCOMPLETE** | 0 | 57 |
| n2_NTN_rks | 0 | 1 | 20 | 20 | 0 | 191 |
| n2_NTN_uks | 0 | 1 | 20 | 20 | 0 | 198 |
| n2h2_NDN_rks | 0 | 1 | 20 | 20 | 0 | 188 |
| n2h2_NDN_uks | 0 | 1 | 20 | 20 | 0 | 203 |
| n2h4_N-N_rks | 0 | 1 | 20 | 20 | 0 | 192 |
| n2h4_N-N_uks | 0 | 1 | 20 | 20 | 0 | 210 |
| ncl3_N-Cl_rks | 0 | 1 | 20 | 20 | 0 | 220 |
| ncl3_N-Cl_uks | 0 | 1 | 20 | 2 **INCOMPLETE** | 0 | 69 |
| nf3_N-F_rks | 0 | 1 | 20 | 20 | 0 | 205 |
| nf3_N-F_uks | 0 | 1 | 20 | 20 | 0 | 260 |
| nh2oh_N-O_rks | 0 | 1 | 20 | 20 | 0 | 220 |
| nh2oh_N-O_uks | 0 | 1 | 20 | 20 | 0 | 249 |
| nh3_N-H_rks | 0 | 1 | 20 | 20 | 0 | 190 |
| nh3_N-H_uks | 0 | 1 | 20 | 20 | 0 | 191 |
| o2_ODO_uks | 0 | 3 | 20 | 16 **INCOMPLETE** | 0 | 174 |
| of2_O-F_rks | 0 | 1 | 20 | 20 | 0 | 204 |
| of2_O-F_uks | 0 | 1 | 20 | 0 **INCOMPLETE** | 0 | 33 |

### class S (contact scans, item 2; charge 0, mult 1)

| system | q | mult | pts | n_ok | rc | wall_s |
|---|---|---|---|---|---|---|
| ch4_h2o_CO | 0 | 1 | 22 | 20 **INCOMPLETE** | 0 | 198 |
| hf_dimer_FF | 0 | 1 | 22 | 20 **INCOMPLETE** | 0 | 198 |
| nh3_h2o_NO | 0 | 1 | 22 | 20 **INCOMPLETE** | 0 | 198 |
| water_dimer_OO | 0 | 1 | 22 | 20 **INCOMPLETE** | 0 | 203 |

### class A, the four new pairs (item 3; charge 0, mult 1)

| system | q | mult | pts | n_ok | rc | wall_s |
|---|---|---|---|---|---|---|
| clf_F-Cl_rks | 0 | 1 | 20 | 20 | 0 | 199 |
| clf_F-Cl_uks | 0 | 1 | 20 | 16 **INCOMPLETE** | 0 | 247 |
| ncl3_N-Cl_rks | 0 | 1 | 20 | 20 | 0 | 213 |
| ncl3_N-Cl_uks | 0 | 1 | 20 | 20 | 0 | 262 |
| nf3_N-F_rks | 0 | 1 | 20 | 20 | 0 | 206 |
| nf3_N-F_uks | 0 | 1 | 20 | 20 | 0 | 240 |
| of2_O-F_rks | 0 | 1 | 20 | 20 | 0 | 187 |
| of2_O-F_uks | 0 | 1 | 20 | 16 **INCOMPLETE** | 0 | 314 |

## Wall time

ORCA wall this session: class H 3.13 h + class S 796 s + new pairs 1354 s + retry 561 s + probes 54 s + 4 reference-geometry Opts (not individually logged) = **~4 h ORCA wall**; real elapsed **~1.3 h** (19:45-21:01). The 16-core budget was held throughout -- 4 jobs x 4 cores, or two campaigns x 2 jobs x 4 cores.

## Not verified

- The EEQ side of the item-1 falsifier (reference far-point charges vs curcuma EEQ) was not computed: it needs the kept-topology mode and belongs to the curcuma side.
- The UKS broken-symmetry far points of of2/clf never converge at either chain direction; no UKS D_e for those two pairs.
- Hirshfeld charges are the reference target, not a fit input; nothing here was compared against CM5.
