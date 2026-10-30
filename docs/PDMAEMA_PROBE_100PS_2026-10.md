# PDMAEMA (AB)5 in a water droplet: evaluation of the 100 ps probe run (Oct 2026)

Status: AI-generated evaluation of AI-generated runs, no human review. n = 1 trajectory per temperature; both temperatures
start from the same structure and the same velocity pattern and see the same CSVR noise sequence (`-md.seed` has no effect,
see the setup in [POLYMERBUILD_CHARGE.md](POLYMERBUILD_CHARGE.md#100-ps-probe-run-setup-oct-2026)). Statements are about
feasibility and stability of this setup only; nothing here is about LCST, collapse or experiment.

System: (AB)5 (+5, 267 atoms) + 933 H2O + 5 Cl-, 3071 atoms, spherical wall 21.7 A, GFN-FF with the group file
(units +1, Cl- -1, waters 0, `rest per_fragment`), `cn_zero_guard`, `topology_mode constant`, dt 0.5 fs, CSVR 10 fs,
303 K and 338 K; production 100 ps (t = 30-130 ps). Data: `~/curcuma_runs/probe100/` (outside the repository).
Criteria and definitions were written down before the first result (plan of the evaluation, kept with the scripts in
`~/curcuma_runs/probe100/analysis/`). Frames: 1000 (0.1 ps, t = 30.1-130.0 ps; the file holds the t = 30.0 ps
frame twice, both skipped), md.log rows 2000 (0.05 ps). Halves H1 = 30-80 ps, H2 = 80-130 ps.

## Stability (production)

| | 303 K | 338 K |
|---|---|---|
| T mean +- SD (H1 / H2), K | 303.00 +- 4.39 (302.99 / 303.01) | 338.07 +- 4.93 (338.11 / 338.04) |
| Epot drift, Eh/ps (SE from 20 blocks of 5 ps) | +1.3e-4 +- 3.8e-4 | **-1.25e-3 +- 3.1e-4** (H1 -357.924, H2 -357.994 Eh) |
| Etot+Wall drift, Eh/ps (CSVR: not conserved) | +1.6e-4 +- 4.0e-4 | -1.21e-3 +- 3.4e-4 |
| frames with an atom beyond 21.7 A / max atoms | 995/1000 / 25 | 1000/1000 / 42 |
| frames with a water O beyond 23.7 A (distinct waters) | 0 (0) | 1 (1) |
| droplet R_rms of water O about the COM, A (sqrt(5/3) R_rms) | 15.28 +- 0.06 (19.72), H1 = H2 | 15.42 +- 0.06 (19.91), H1 = H2 |
| O density 0-8 / 8-12 / 12-16 / 16-20 A, A^-3 (build 0.0333) | 0.019 / 0.031 / 0.033 / 0.023 | 0.018 / 0.030 / 0.033 / 0.023 |
| polymer COM from droplet COM, A (10-ps blocks first -> last) | 2.20 (1.37 -> 3.49) | 1.90 (1.07 -> 3.20) |
| smallest polymer-atom distance to the wall, A: mean (min) | 5.44 (2.68); H1 6.49, H2 4.39 | 5.22 (1.82); H1 5.99, H2 4.45 |
| polymer SASA with water and Cl- as occluders, A^2 (H1 / H2) | 21 (6.8 / 35.3) | 19 (2.9 / 35.2) |
| single points (20 frames): fragments 939, groups accepted, "ignored" | 20/20, 20/20, 0 | 20/20, 20/20, 0 |
| group sums: max abs(sum - Q) / water max abs(sum) | 5e-10 / 1.0e-10 | 6e-10 / 1.0e-10 |

- The droplet keeps its size (R_rms H1 = H2 to 0.02 A); the 0-8 A shell holds the chain, so its O density is low by
  construction; the 16-20 A shell is the surface.
- The chain moves off-centre in both runs (COM offset grows 2.3-3x over 100 ps) and in H2 about 35 A^2 of its surface are
  not covered by water (H1: 3-7 A^2). The data show a chain approaching the droplet surface; whether this is surface
  adsorption, a wall effect or a slow drift of one trajectory cannot be told from n = 1.
- 338 K: Epot falls by 0.07 Eh between the halves at constant T (4 sigma of the block-based slope error); not stationary.
- Single-point energy of a frame (fresh topology) minus md.log Epot of the same step: +3.1 / +2.9 mEh mean, -1.9 to +8.6 mEh
  (n = 40); same sign as the 1.2 mEh optimiser/single-point gap of the minimisation. Not root-caused.

## Polymer (mass-weighted; whole production, block error = largest SE with >= 5 blocks)

| quantity | 303 K | 338 K | 338 - 303 |
|---|---|---|---|
| Rg heavy atoms (110), A | 7.683 +- 0.028 (H1 7.617, H2 7.750) | 7.739 +- 0.038 (H1 7.805, H2 7.673) | +0.06 +- 0.05; H1 +0.19, H2 -0.08 |
| Rg all 267 atoms, A | 7.763 +- 0.027 | 7.814 +- 0.039 | +0.05 +- 0.05 |
| Rend (terminal alpha-C 1 and 240, 18 bonds), A | 19.30 +- 0.18 (H1 18.98, H2 19.62) | 19.36 +- 0.25 (H1 19.88, H2 18.84) | +0.06 +- 0.31 |
| kappa^2 (rel. shape anisotropy, heavy atoms) | 0.395 +- 0.013 | 0.413 +- 0.010 | +0.02 +- 0.02 |
| asphericity b, A^2 | 37.1 +- 0.9 | 38.4 +- 0.8 | +1.4 +- 1.2 |
| Rend / Rg_heavy | 2.51 +- 0.02 | 2.50 +- 0.02 | -0.01 +- 0.03 |
| tau_int Rg_heavy (ACF), ps; N_eff blocks / ACF | 10.1; 14.5 / 4.9 | 7.4; 12.5 / 6.8 | |

- Stationarity (H2 - H1 against the combined block error of the halves): Rg_heavy +4.7 sigma (303) and -2.1 sigma (338);
  Rend +2.5 / -3.2 sigma; kappa^2 +3.5 / -0.7 sigma. Rg drifts in both runs, in opposite directions; 10-ps block means of
  Rg_heavy span 7.55-7.81 A (303) and 7.58-7.92 A (338).
- tau_int of the shape quantities is 5-15 ps, the effective sample of 100 ps is 3-15 independent values (block and ACF
  estimates differ by a factor 2-3; both are named). Block lengths of 12.5 ps (8 blocks) set the error.
- Temperature comparison: every difference is below 3 sigma and changes sign between H1 and H2: **nicht belastbar**.

## Hydration and hydrogen bonds

Geometric criterion of `-trajprops`: H to the nearest N/O within 1.3 A, D...A <= 3.5 A, D-H...A >= 150 deg (Cl as
acceptor with the same limits). Per frame, mean +- block error (303 / 338 K):

| bonds | 303 K | 338 K |
|---|---|---|
| N+-H -> water (5 donors) | 2.63 +- 0.09 | 2.57 +- 0.09 |
| N+-H -> Cl- | 0.01 +- 0.01 | 0.01 +- 0.01 |
| water -> carbonyl O, protonated / neutral units (5 O each) | 10.9 +- 0.4 / 12.4 +- 0.3 | 10.3 +- 0.2 / 11.1 +- 0.3 |
| water -> ester O, protonated / neutral units | 0.9 +- 0.2 / 1.5 +- 0.3 | 0.4 +- 0.1 / 1.3 +- 0.1 |
| water -> amine N (5 neutral units) | 4.77 +- 0.17 | 4.75 +- 0.23 |
| water -> Cl- (5 ions) | 35.9 +- 0.2 | 34.1 +- 0.2 |
| polymer intra | 0.02 | 0.03 |

g(r) over 1000 frames, bin 0.05 A; normalised to the wall sphere (V = 4/3 pi 21.7^3), so heights are relative, positions
and n(r) are not affected. First maximum / first minimum (n within the minimum):

| pair | 303 K | 338 K |
|---|---|---|
| Cl-...O(w) | 3.28 / 3.98 A, n 9.3 | 3.28 / 3.98 A, n 9.2 |
| Cl-...H(w) | 2.33 / 3.08 A, n 9.1 | 2.33 / 3.08 A, n 9.0 |
| N+...O(w) (local max, see below) | 3.08 / 3.23 A, n 1.0; second peak 4.23, min 5.63, n 20.1 | 3.08 / 3.33 A, n 1.3; 4.18, min 5.03, n 13.6 |
| N+...Cl- | no peak (31 pairs within 4 A in 1000 frames) | no peak (36 pairs) |
| Cl-...H(N+) | 2.33 / 2.43 A, n 0.002, from 55 pairs within 4 A | 2.68 / 2.88 A, n 0.004, from 138 pairs |

- N+...O: the pre-set rule (largest g below 4 A) landed on the 4 A edge of the second peak; the table gives the first local
  maximum instead (deviation from the plan, made after seeing g). The 4.2 A peak is the methyl shell; its minimum is shallow,
  so n differs by 6.5 between the runs from where the minimum falls.
- SASA of the chain alone (Shrake-Rupley, Bondi radii H 1.20, C 1.70, N 1.55, O 1.52 A, probe 1.4 A, 960 points, occluders =
  chain): 1954 +- 10 A^2 (303, H1 1959, H2 1949) and 1924 +- 18 A^2 (338, H1 1963, H2 1885). The radius set of the WASA in
  Min et al. 2015 is not known here; no number is compared with that paper.
- Energy criterion of `-trajprops` (GFN-FF triple energies) not usable for this run: it builds GFN-FF with the total charge
  only (`trajprops.cpp`, `countHydrogenBondsByEnergy`: no charge groups, no `cn_zero_guard`), i.e. another charge model
  (Cl- neutral), and with `hb_energy_limit 0` it counts 9445 polymer-rest triples per frame against 31 geometric ones
  (2 frames, 1.1 s and 0.9 GB per frame). A limit would need a calibration that was not done.

## Cl- at the N-H

| | 303 K | 338 K |
|---|---|---|
| smallest H(N+)...Cl- over the production, A | 2.167 | 2.254 |
| frames with any contact < 2.4 / < 2.0 / < 1.8 A (of 1000) | 8 / 0 / 0 | 5 / 0 / 0 |
| ions with contacts, events (consecutive frames < 2.4 A) | 1 (atom 3066 at N-H 53), 6 events, 1-3 frames | 1 (atom 3069 at N-H 265), 4 events, 1-2 frames |
| mean H...Cl in contact frames, A | 2.26 (n = 8) | 2.32 (n = 5) |
| 1 % / 5 % / 50 % percentile of the frame-wise smallest H...Cl, A | 2.48 / 3.92 / 5.37 | 2.61 / 3.37 / 5.17 |
| heating and equilibration (context): smallest H...Cl, A | 3.62 / 3.28 | 3.66 / 4.09 |

- No frame came below 2.0 A, so the energy step of the bond detection at 1.7-1.8 A was not reached and the frozen groups were
  not tested by a change of fragmentation. Single points on the 25 frames around the closest contacts: 939 fragments,
  groups accepted, 0 "ignored"; frame-to-frame energy jumps J = abs(E_k - (E_k-1 + E_k+1)/2) 12-169 mEh against md.log Epot
  over the production (0.1 ps): median 53 / 61 mEh, 99th percentile 212 / 227 mEh, i.e. no outlier at the contacts.
- Contact distances 2.26 / 2.32 A sit where the microsolvation MD with groups put the Cl- (mean 2.31 A, GFN2 1.91 A). Two
  ions, 13 frames: consistent with the earlier finding, not a test of it. Contacts shorter than 0.1 ps are not resolved.

## Tool check (`-trajprops` on 1000 frames x 3071 atoms, 118 MB file, one thread)

| quantity | wall time | peak RSS | agreement with an independent numpy evaluation |
|---|---|---|---|
| gyration (110 heavy atoms / 267 atoms) | 0.9 s | 41 MB | max abs diff 1.2e-5 / 1.6e-5 A over 1000 frames (CSV prints 6 digits; masses not compared) |
| hbonds polymer vs rest (N, O) | 2.2 s | 41 MB | per direction identical in 1000/1000 frames, both T |
| sasa (occluders group / all) | 4.6 / 5.5 s | 40 MB | not checked against another program |

The tool's printed block error uses 5 fixed blocks (Rg 303 K: 0.035 A) against 0.028 A from the scan rule here; both are
block errors of 8 or fewer blocks. Its literature references remain unchecked (TODO.md).

## Feasibility

- Stages per temperature (two runs side by side, 8 threads each, `nice`): heating 10 ps 1869 / 1865 s (93 ms/step),
  equilibration 20 ps 3745 / 3755 s (94 ms/step), production 100 ps 18832 / 18857 s (94 ms/step); chain 14:35 -> 21:23,
  6.8 h per temperature. Minimisation 1273 steps not converged (gradient norm 3.0e-3), shared.
- Peak RSS grows with the stage length: 2.0 / 2.6 / 6.8 GB (303 K) and 2.1 / 2.6 / 7.5 GB (338 K) for 10 / 20 / 100 ps;
  cause found and fixed on Oct 26, 2026: `SimpleMD::CSVR()` appended one double per atom and step to the never-read `m_atom_temp` (removed; [MD_LARGE_SYSTEMS.md](MD_LARGE_SYSTEMS.md)). The runs described here used the binary before the fix. A follow-up with the larger droplet: [PDMAEMA_PROBE_LARGE_DROPLET_2026-10.md](PDMAEMA_PROBE_LARGE_DROPLET_2026-10.md).
- Disk: 159 / 162 MB per temperature, of which 118 MB is the production trajectory (about 1.2 GB per ns at 0.1 ps dumps).
- Health checks K1-K7 passed in all 7 stage checks; single point 4.9 s and 0.9 GB at one thread.
- Known: `-md.seed` without effect (same start velocities and CSVR noise at both T); Etot is not conserved under CSVR, so the
  jump check is a hint only; optimiser/single-point gap 1.2 mEh.

## Limits

n = 1 per temperature with correlated starts; 100 ps against tau_int 5-15 ps of the chain shape (3-15 independent values);
10 units, a droplet with a harmonic wall, no PBC/NPT; formal charges +1/-1/0; GFN-FF water and Cl- hydration not validated
(Cl- n = 9 in the first shell is a model number, not compared); g(r) normalised to the wall volume; SASA radius set chosen
here. Human production testing pending.
