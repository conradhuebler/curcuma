# GFN-FF: dependence of the energy on the atom numbering

> 🤖 AI-generated, ⚙️ machine-tested (Nov 24, 2026). Human production testing pending. Updates
> [REV_GFNFF_TODO.md](REV_GFNFF_TODO.md) #12 and [KNOWN_ISSUES_ARCHIVE.md](KNOWN_ISSUES_ARCHIVE.md) #34, #42.

## Result in short

Two mechanisms carry the whole spread seen on caffeine, benzoic acid, the acetic acid dimer and
`complex` (31 orders each). Both are in the Fortran sources, curcuma ports them faithfully
(same order -> same term energy as pprcht), so neither is a port error:

| Term | Mechanism | pprcht | xtb 6.7.1 | curcuma default | Class | Switch (default off) |
|---|---|---|---|---|---|---|
| inversion | `ssort` of pprcht stops early, neighbour roles depend on input order | spreads | no spread | = pprcht | R (pprcht translation bug; xtb is correct) | `-gfnff.inversion_sort_xtb true` |
| Coulomb | nitrogen dgam tests `piadr(i)` (compact list, an index cutoff) | spreads | spreads | = pprcht | R (both references) | `-gfnff.dgam_pi_membership true` |
| bond, repulsion, ... | summation order (rounding) | <= 6 nEh | | <= 12 nEh | T | none |

With both switches the spread is gone on these structures (<= 12 nEh, gradient <= 3e-9 Eh/A).
Known Issue #42 (NH3 kink) is a third consequence of the same design: the inversion term assigns
the neighbour roles by distance rank, and the rank changes at a distance tie (class R, also with xtb's sort).

## Method

- `scripts/gfnff_order_scan.py`: N seeded row permutations (coordinate lines textually identical),
  one `-sp -method gfnff -gfnff.cache_topology false -threads 4 -verbosity 3 -dump_gradient` per
  order, terms from the 17-digit `CPU ENERGY TERMS` block, gradient mapped back to the original
  order. `--program pprcht|xtb` runs the references; `--compare A B` gives |curcuma - reference|
  per order and term.
- References: pprcht/gfnff `0491df2` (2025-01-27, port source; built from `git archive` under
  `$HOME/tmp_build`, `-Dbuild_exe=ON`, gfortran, `OMP_NUM_THREADS=1`), instrumented copy prints
  each inversion term and the dihedral/inversion split (`PPRCHT_INVDBG`); xtb 6.7.1 `edcfbbe`
  (`/opt/bin/xtb --gfnff`). Text comparison also against the xtb source `2b5cd48` (2026-05-17).
- Diagnostic builds of pprcht (not used as reference): `fixsort` (xtb's ssort) and `fixall`
  (fixsort + `piadr2` in the dgam rule).
- curcuma: `af140bb4` (reference build) and this branch; n = 31 orders (identity + 30) per
  structure; `test_cases/structures/medium/caffeine.xyz` gives the same numbers as
  `test_cases/molecules/larger/caffeine.xyz`.

## Reproduction (span = max - min over 31 orders, nEh)

| Structure | curcuma inversion | curcuma Coulomb | pprcht inv. / Coul. | pprcht fixall | xtb 6.7.1 total | curcuma + both switches, total |
|---|---:|---:|---:|---:|---:|---:|
| caffeine | 0 | 123197 | 0 / 123197 | 0 / 0 | 123198 | 0.8 |
| C6H5COOH | 1021 | 0 | 1021 / 0 | 0 / 0 | 0.1 | 1.3 |
| acetic_acid_dimer | 914580 | 0.2 | 914580 / 0 | 0 / 0 | 0.1 | 0.6 |
| complex (231 atoms) | 911069 | 607348 | 911069 / 607348 | 0 / 0.3 | 607344 | 12.0 |
| C6H6 | 5.6 | 0 | 5.6 / 0 | 0 / 0 | 0 | 0.0 |
| triose | 0 | 0.4 | 0 / 0.1 | 0 / 0.1 | 0 | 0.7 |

- Handover findings confirmed: Coulomb spread on caffeine and complex (identical in xtb), inversion
  spread on the acetic acid dimer, benzoic acid and complex. Not confirmed: "the inversion spread
  exists only in curcuma" - pprcht has it to the last digit; only xtb (6.2.2 and 6.7.1) does not.
- Angle, dihedral, sTors, BATM, dispersion, HB, XB: <= 0.5 nEh in every structure. Bond up to
  12.4 nEh (complex), repulsion up to 1.2 nEh: summation order, also in pprcht (6.3 nEh).
- Gradient, max |g(order) - g(identity)|: default 1.1e-2 Eh/A (acetic acid dimer), 5.2e-3
  (complex), 9.4e-5 (caffeine); with both switches <= 2.9e-9 Eh/A.
- Same order, curcuma (default) vs pprcht (default): inversion <= 1.6e-9 Eh, Coulomb <= 6.7e-8 Eh
  in all 6 x 31 orders, so the spreads are ported, not created. Bond and repulsion differ by
  1e-7 .. 2.6e-6 Eh at EVERY order (unit constants, Known Issue #35); the 1e-8 per-term gate
  against pprcht is therefore not met for those terms, independent of the order.

## Source comparison (Fortran <-> C++)

| Routine | pprcht `0491df2` | xtb `2b5cd48` | curcuma | Verdict |
|---|---|---|---|---|
| inversion list: centres | `gfnff_ini.f90:1841-1845` nb=3 and (pi or N) | `gfnff_ini.f90:2041-2045` same | `gfnff_inversions.cpp:359-382` (post-Hueckel piadr) | equal |
| inversion: neighbour order | `:1851-1860` sort by distance, `call ssort(3,...)` | `:2079-2087` same | `gfnff_inversions.cpp:386-420` | equal call |
| `ssort` | `gfnff_helpers.f90:278` "if (edum(j) .gt. pp) exit", `:282` "if (k .eq. i) exit" | `gfnff_ini2.F90:1620` "go to 120", `:1624` "go to 140" (= next iteration) | `break` / `break` (= pprcht); `continue` with `inversion_sort_xtb` | pprcht differs from xtb |
| inversion parameters | `:1865-1895` (fc from repz, pbo, qa, ncarbo/nf) | same | `gfnff_inversions.cpp:435-495` | equal (nf: C++ lists F, Cl, Br, I; pprcht `group == 7` also At) |
| inversion energy | `gfnff_engrad.F90:1094-1126` `omega(i,j,k,l)`, damping over i-j, j-k, j-l | same | `ff_workspace_gfnff.cpp:686` | equal (no symmetrisation in any of them) |
| dgam, nitrogen | `gfnff_ini.f90:650` "if (piadr(i) .ne. 0) ff = -0.14", `:651` "amide(...,piadr,i)" | `gfnff_ini.f90:685-686` same | `eeq_solver.cpp:5153` (`i < npiall`), `detectAmideNitrogens` | equal (index cutoff kept) |
| `piadr` at that point | `:332` "piadr(k) = i" (k-th pi candidate -> atom), so `piadr(i) /= 0` iff i <= npi | same | `npiall` emulation | equal |
| amide-H acidity | `gfnff_ini.f90:793-798` | `:870-878` | `gfnff_method.cpp:1997` | curcuma deviates by default (#11, switch exists) |
| GEODEP sp2 -> sp3 (C, CN 3, no pi) | `:1027-1031` `omega` over unsorted nb(1..3), > 40 deg | `:1166` same | `gfnff_method.cpp:12049` | equal; order-dependent only near 40 deg, not seen here |
| rings | `getring36` with `ssort` (`gfnff_helpers.f90:129`) | correct ssort | full enumeration (`findSmallestRings`) | different code, not seen to matter here |
| fragment charge | net charge on fragment of atom 1 | same | `frag_charge_model ensemble` default | #13/#31, not part of this study |

## Cause per mechanism

- **Inversion (`ssort`)**: the out-of-plane angle of centre i is `omega(i,j,k,l)` = angle of i->l
  against the plane (i, j, k), damped over i-j, j-k and j-l. The roles must come from the
  distance sort ("such that the same inversion angle def. always results", `gfnff_ini.f90:1850`).
  pprcht's `ssort` leaves (r1, r3, r2) and (r2, r3, r1) unsorted, so the role assignment, and with
  it the energy of every non-planar centre, follows the input order. xtb's original selection sort
  does not. Same order: curcuma == pprcht to 1.6e-9 Eh; pprcht with xtb's sort: spread 0.
- **Coulomb (dgam)**: before the Hueckel step `piadr` is the compact list of pi candidates, so the
  nitrogen test `piadr(i) /= 0` asks "is the atom index <= number of pi candidates". Renumbering
  moves nitrogens across that cutoff (caffeine: N2, N4, N8, N11 change gameeq by 0.9e-3 .. 1.8e-3).
  Same in xtb. pprcht with `piadr2` there: spread 0 (caffeine, complex).
- **#42**: NH3 with N shifted along x, fresh topology per point (as the FD gate does): for dx < 0
  the roles are (j,k,l) = (H4,H3,H2), for dx >= 0 (H2,H3,H4) (xtb sort: (H2,H4,H3)); the
  inversion slope jumps from 8.4e-5 to 1.8e-6 Eh/A at dx = 0, the #42 numbers. The analytic
  gradient belongs to one fixed role set, so analytic and numeric gradient cannot agree across a
  rank change. Not fixable without changing the term (symmetrising over the three choices of l).

## Effect of the switches in the ORIGINAL order

- Six structures above: only `complex` moves (inversion -0.249, Coulomb -0.059 kcal/mol); the
  inversion part equals the curcuma(default) - xtb difference of the total (+0.249 kcal/mol).
- MOR41 (95): 36 structures change, MAD 0.0022, max 0.061 kcal/mol (ED26). GMTKN55 (2462): 598 change, MAD 0.053, max 9.73 kcal/mol (MB16-43/12; ICONF/N3P3H12_1 7.84). `scripts/refset_regression.py` against `af140bb4`; with the defaults 0 of 2557 change (bitwise).

## Canonical-order switch (`-gfnff.canonical_order`): assessed, not implemented

Feasible as a renumbering at `GFNFF::InitialiseMolecule` with back-permutation of gradient,
charges, bond orders, HB/XB lists and the topology cache key, but every consumer that reads
`GFNFF` state by atom index (MD state, `-gfnff.charge_groups` files, react mode, GPU wrapper,
per-triple HB energies, ConfSearch reuse) would need the map too. It would fix the energy to an
arbitrary canonical order, not to a more correct one, while the two switches above remove the
measured spread at their source. Left as a TODO work package.

## Tests

- ctest `gfnff_order_invariance` (`test_cases/scripts/check_gfnff_order_invariance.py`): C6H5COOH, acetic
  acid dimer, caffeine, identity + 8 seeded permutations, `-threads 1`. Default: every term except the
  reference-caused ones within 1e-8 Eh; inversion (C6H5COOH 5.9e-7, dimer 4.7e-4 Eh) and Coulomb
  (caffeine 6.6e-5 Eh) must SPREAD (> 1e-7 Eh, documents the reference behaviour). Both switches: every
  term within 1e-8 Eh (inversion <= 5e-12, Coulomb <= 1.2e-10). 1e-8 and not 1e-9, because the bond term
  spreads by up to 1.7 nEh from summation order on these structures.
- `ctest -R gfnff_param_fallbacks` passes (PARAM default == code fallback for the two switches).
- Gates (this branch vs `af140bb4`): `bitcheck_merge.sh` (6 structures x gfnff/gfn1/gfn2 x 1/4 threads + 2 MD) 38 of 38 lines identical; full `ctest -j1` 498/500, the 2 failures pre-existing (`cli_confsearch_05_seed_cold_bias`, `cli_deltaml_04_label_consistency`); MOR41 + GMTKN55 defaults bitwise unchanged.

## Found on the way

- **Undefined behaviour in the topology build, fixed**: `GFNFF::classifyBondType` (metal-metal case,
  `imetal_of`) called `getCachedTopology()` from inside `calculateTopologyInfoOnce()`, where the
  cache is still empty. With both switches on, AL2X6/al2me4, al2me5 and ALK8/li4_me4 crashed
  (SIGSEGV) depending on the input path length (heap layout); default runs did not crash in 36 of
  36 tries. Now the lists of the topology being built are passed in; default energies bitwise
  unchanged on MOR41 (95) and GMTKN55 (2462); 18 of 18 switch runs complete.
- `GFNFF::calculateDgam` (`gfnff_method.cpp`) has no caller (dead code, not removed here).

## Audit MOR41 + GMTKN55 (2557 structures, identity + 3 seeded permutations each)

Largest term spread per structure (kcal/mol), counted by the term that carries it:

| Class | Term | default n (> 1e-6) | n > 0.01 | max | both switches n (> 1e-6) | n > 0.01 | max |
|---|---|---:|---:|---:|---:|---:|---:|
| neutral | Coulomb | 515 | 482 | 1.19 (IL16/150) | 0 | 0 | - |
| neutral | inversion | 330 | 105 | 9.87 (MB16-43/12) | 0 | 0 | - |
| charged | Coulomb | 92 | 87 | 435 (DIPCS10/ch2o_2+) | 35 | 32 | 435 |
| neutral | bond | 19 | 2 | 19.6 (MB16-43/23) | 141 | 2 | 19.6 |
| charged | bond | 8 | 5 | 64.3 (WATER27/OHmH2O) | 36 | 5 | 64.3 |
| charged | inversion | 19 | 4 | 0.016 | 0 | 0 | - |

Structures with any term spread > 0.01 kcal/mol: 685 -> 39. The rise of the "bond" rows with the
switches is the rounding floor (1e-6 kcal/mol = 1.6e-9 Eh) that becomes the largest term once the
others are gone. What stays is the charged multi-fragment class (#13/#31, net charge placement) and
the MB16-43 bond cases (not traced). (The Sep 2026 audit counted 911 with another permutation set and
threshold; the numbers are not directly comparable.)

## Not done / open

- Charged Coulomb (35 structures, up to 435 kcal/mol) and bond spreads (MB16-43/23, WATER27/OHmH2O): not traced here.
- No r2SCAN-3c verdict on whether the switched values are better: they only make the energy
  independent of the numbering; the parameters were fitted with the reference behaviour.
- GPU: the inversion list and dgam are host-side inputs, so CUDA/ROCm take the switches without
  kernel changes; not run on a GPU.
