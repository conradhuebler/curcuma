# H2 dissociation on frozen Pt_n clusters (static, GFN2), Nov 21, 2026

Status: 🤖 AI-generated study, ⚙️ machine-tested only (native vs tblite, see below). Not ✅ TESTED: human production testing is pending until the operator removes this line.
Statement of the model GFN2 (native, `-xtb.spin_polarization`, [SPIN_POLARISATION.md](SPIN_POLARISATION.md)); there is no DFT reference in this build, so this is not a comparison with DFT.

## Method

- **Clusters** n = 1, 2, 3, 4, 5, 6, 7, 10, 13: `-opt` from 1 to 3 idealised motifs per n, jitter 0.05 A, GFN2, 300 K, `-spin` 0..2n (max 12), with and without spin polarisation; no global search. Frozen frame = lowest-energy geometry with spin polarisation (Pt4 planar rhombus/square; the Td tetrahedron of `reactions-doebereiner-pt4` collapses without symmetry and is run as a separate frame `pt4tet`, 17.6 kcal/mol above; Pt13: `pt13` lowest with polarisation, `pt13ico` lowest without). Pt is frozen at that one geometry for all methods.
- **Sites**: top, bridge, hollow found from the geometry (up to 3 symmetry-distinct per type); H2 parallel to the surface, axis along the first neighbour / Pt-Pt / triangle edge.
- **Surface**: rigid grid, h (centre height above the site) 3.5..1.2 A x d(H-H) 0.74..2.6 A, `-sp` at every spin 0..2n (max 12), both modes, about 70 000 points, 73 not converged (marked, not used), 707 needed retry settings (`scf_damping`, `-scf_guess h0`).
- **Barrier / product** (rigid grid is only scouting): (i) relaxed d-scan (Pt frozen, H pair: centre + 2 angles relaxed at fixed d, forward and backward, lower envelope), (ii) relaxed approach scan (h fixed, d free), (iii) free relaxation of both H from d = 2.6 (product). Barrier = highest point of the best path relative to separated Pt_n + H2, minimum of the upper bounds (ii) and max((ii) entrance, (i)). Values <= 0 mean barrierless along the paths found; they are not saddle points. No NEB, no Hessian check.
- **Methods**: A = spin 0, no spin polarisation. B = spin polarisation, lowest of spin 0, 2, 4, 6, 8 per point (spin may change along the path, spin conservation is not enforced). C = spin polarisation at the fixed ground-state spin of the bare cluster (spin conserved, H2 singlet). All energies relative to the bare frozen cluster in the method's own lowest spin + relaxed H2 (0.7765 A).
- Scripts (Python, driver around `-sp`/`-dump_gradient`) and raw data are in the vault tarball, not in `scripts/`.

## Result (300 K, lowest-barrier site per cluster, kcal/mol, see the plot in the vault `_artefacts/`)

Barrier: <= 0 (no barrier along the paths found) at the best site for every n and method (-1 to -5; the negative value is the long-range start at h = 3.5 A). Other sites: up to +4 to +10 (hollow sites, Pt5, Pt13), A and B alike. Reaction energy to two separate H: Pt1 -68, Pt2 -35, Pt3 -39, Pt4 -54, Pt5 -51 to -57, Pt6 -42 to -46, Pt7 -38 (B) / -49 (A), Pt10 -29, Pt13 -32; weaker for the larger, more saturated clusters.
Spin polarisation: barrier changes B-A are mean 0 to 0.5 kcal/mol per cluster, except Pt2 (+2.9, up to +5.7) and single sites (up to -3.6 / +1.0). It changes the reaction energy: +3 to +10 kcal/mol less exothermic for Pt5, Pt6, Pt7 (the bare high-spin cluster is stabilised by 12 to 19 kcal/mol against spin 0, and the product is low-spin), C (spin conserved) a further +4 to +10. In B the product spin differs from the bare-cluster spin in 6/7 (Pt5), 7/7 (Pt7) sites.

## Checks and caveats

- Native against tblite 0.8.0 (Pt4, Pt5, Pt7, 97 geometries: frame, d = 0.85, d at the barrier, d = 2.6, product; spin 0..8; 300 K): 90 of 97 differ by <= 1e-8 Eh, max 4.9e-3 Eh; all 7 outliers are Pt7 spin 2 or 4 (another SCF solution, as in SPIN_POLARISATION.md; native is the higher one in 5 of 7, the lower in 2). Not tested: other sizes, GFN-FF.
- 1000 K against 300 K at the 300 K geometries (d-scan only, Pt2, Pt5, Pt7): A barrier +4.0 (Pt2), +0.3 (Pt5), +2.3 kcal/mol (Pt7); B within 0.5 to 0.7 (Pt2) and 0.4 (Pt5, Pt7); the product shifts up to 2.6 kcal/mol.
- SCF has several solutions (Pt7 spin 2, multi-minima); lowest-spin selection takes the lowest found, not the proven lowest. Relaxed paths are local, hysteresis between forward and backward scans is only reduced by the lower envelope. Pt-Pt geometry is one frozen frame; the Pt clusters are not global minima. One approach path (Pt6 hollow2) left the site at h = 3.5 and was excluded.
- Not tested: dependence on the H2 orientation beyond the one parallel orientation, real saddle points, zero-point energy, spin-orbit coupling, solvent. Pt-H binding energies of GFN2 are unchecked against experiment or DFT; any literature comparison would be from memory and unchecked.
