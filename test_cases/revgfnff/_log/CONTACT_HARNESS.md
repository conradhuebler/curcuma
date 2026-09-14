# Class-S contact scans: revgfnff vs gfnff vs r2SCAN-3c (interaction curves)

AI-generated (`scripts/revgfnff_contact.py`), machine-evaluated. kcal/mol, Angstrom.

- binary: `/home/conrad/src/curcuma_branches/curcuma/build_rev/curcuma` md5 `3d76c8c9d37f6a0f9d54785d45ddc903` mtime 2026-09-14 14:27:49
- topology frame: the LARGEST-separation frame (frame 0 of the multi-frame XYZ, `-batch_reuse_topology true`); `-gfnff.cache_topology false`, fresh temp dir per run, CPU/`-gpu none`
- interaction energy E_int(d) = E(d) - E(6.00 A) (the largest separation; the two isolated-monomer ORCA jobs could not be computed, see ORCA_REF_STATUS.md)
- curve modes (both always computed): `kept` (topology of the largest separation kept over the whole scan; the verdict curve) and `fresh` (re-perceived at every point; diagnostic)
- acceptance: max |E_int(revgfnff) - E_int(gfnff)| <= 1.0 kcal/mol on `kept`

## Reference minima reproduced from the raw energies.json

| system | sampled d / E_int | parabolic d / E_int | reported d / E_int | verdict |
|---|---|---|---|---|
| ch4_h2o_CO | 3.80 / -0.513 | 3.825 / -0.514 | 3.80 / -0.51 | **AGREE** |
| hf_dimer_FF | 2.90 / -2.741 | 2.940 / -2.752 | 2.90 / -2.74 | **AGREE** |
| nh3_h2o_NO | 3.25 / -1.918 | 3.322 / -1.951 | 3.25 / -1.92 | **AGREE** |
| water_dimer_OO | 3.00 / -3.837 | 3.013 / -3.838 | 3.00 / -3.84 | **AGREE** |

## Per-system comparison

| system | topo frame idx (d) | n | pts > tol | max \|rev-gfnff\| (at d) | max \|rev-ref\| (at d) | max \|gfnff-ref\| (at d) | verdict |
|---|---|---:|---:|---|---|---|---|
| ch4_h2o_CO | 19 (6.00) | 20 | 2 | 32.838 (2.30) | 16.782 (2.30) | 16.057 (2.30) | **FAIL** |
| hf_dimer_FF | 19 (6.00) | 20 | 0 | 0.945 (2.30) | 3.745 (2.60) | 3.745 (2.60) | **PASS** |
| nh3_h2o_NO | 19 (6.00) | 20 | 1 | 15.593 (2.30) | 10.134 (2.30) | 5.458 (2.30) | **FAIL** |
| water_dimer_OO | 19 (6.00) | 20 | 1 | 5.439 (2.30) | 3.584 (2.30) | 2.271 (2.80) | **FAIL** |

**Tolerance verdict: FAIL** - worst system 32.838 kcal/mol against the 1.0 kcal/mol bound; 1/4 pass.
Secondary (informational, never instead of the verdict): the same maximum restricted to d >= 2.50 A is **0.131 kcal/mol** for all four systems - below that the rigid contact is inside the reference's own hard wall (E_int >= +1 kcal/mol there), where stage-1's over-coordination term is the whole difference (term attribution below).
`fresh` mode for reference: worst max |rev - gfnff| = 124.787 kcal/mol (perception flips included, not part of the verdict).

Absolute-energy asymptote E(d=6.00 A), rev - gfnff, in kcal/mol: ch4_h2o_CO -0.24070, hf_dimer_FF -0.02448, nh3_h2o_NO -0.15022, water_dimer_OO -0.18650 - constant and small; it is carried entirely by the over-coordination term and E_int cancels it by construction.

## Minima positions (parabolic over the sampled bracket)

| system | ref d / E_int | gfnff d / E_int | revgfnff d / E_int | rev - ref at its own min |
|---|---|---|---|
| ch4_h2o_CO | 3.825 / -0.514 | 3.623 / -0.301 | 3.623 / -0.301 | +0.213 |
| hf_dimer_FF | 2.940 / -2.752 | 2.795 / -6.055 | 2.795 / -6.055 | -3.303 |
| nh3_h2o_NO | 3.322 / -1.951 | 3.078 / -4.161 | 3.078 / -4.160 | -2.208 |
| water_dimer_OO | 3.013 / -3.838 | 2.935 / -5.935 | 2.935 / -5.930 | -2.092 |

## Term attribution at the worst point (rev - gfnff, kcal/mol, |d| > 0.001)

- **ch4_h2o_CO** (d = 2.30): OverCoord +65.109, RepulsionNonbonded -32.367, Bond -0.144 (sum +32.597)
- **hf_dimer_FF** (d = 2.30): OverCoord +2.235, RepulsionNonbonded -1.296, Bond -0.019 (sum +0.920)
- **nh3_h2o_NO** (d = 2.30): OverCoord +29.371, RepulsionNonbonded -13.834, Bond -0.095 (sum +15.443)
- **water_dimer_OO** (d = 2.30): OverCoord +10.731, RepulsionNonbonded -5.335, Bond -0.143 (sum +5.252)

## Sanity checks

| system | repulsive at smallest d | E_int at 2nd-largest d (ref/gfnff/rev) | zero agrees kept vs fresh |
|---|---|---|---|
| ch4_h2o_CO | yes (+48.4/+32.4/+65.2) | d=5.50: -0.042/-0.022/-0.022 | gfnff yes / revgfnff yes |
| hf_dimer_FF | yes (+5.2/+1.9/+2.9) | d=5.50: -0.128/-0.161/-0.161 | gfnff yes / revgfnff yes |
| nh3_h2o_NO | yes (+25.6/+20.2/+35.8) | d=5.50: -0.109/-0.095/-0.095 | gfnff yes / revgfnff yes |
| water_dimer_OO | yes (+10.7/+8.9/+14.3) | d=5.50: -0.158/-0.151/-0.151 | gfnff yes / revgfnff yes |

## Per-point curves (kept mode, kcal/mol, interaction energy vs d = 6.00 A)

| curve | 2.30 | 2.40 | 2.50 | 2.60 | 2.70 | 2.80 | 2.90 | 3.00 | 3.10 | 3.25 | 3.40 | 3.60 | 3.80 | 4.00 | 4.25 | 4.50 | 4.80 | 5.10 | 5.50 | 6.00 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ch4_h2o_CO:ref | +48.44 | +34.23 | +24.04 | +16.71 | +11.45 | +7.68 | +4.99 | +3.09 | +1.78 | +0.56 | -0.09 | -0.44 | -0.51 | -0.47 | -0.37 | -0.27 | -0.17 | -0.10 | -0.04 | +0.00 |
| ch4_h2o_CO:gfnff | +32.38 | +18.79 | +11.34 | +9.75 | +5.95 | +3.43 | +1.86 | +0.90 | +0.34 | -0.09 | -0.26 | -0.30 | -0.27 | -0.23 | -0.17 | -0.13 | -0.08 | -0.05 | -0.02 | +0.00 |
| ch4_h2o_CO:revgfnff | +65.22 | +28.10 | +11.47 | +9.74 | +5.95 | +3.43 | +1.86 | +0.90 | +0.34 | -0.09 | -0.26 | -0.30 | -0.27 | -0.23 | -0.17 | -0.13 | -0.08 | -0.05 | -0.02 | +0.00 |
| ch4_h2o_CO:**rev-gfnff** | +32.84 | +9.32 | +0.13 | -0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | -0.00 | -0.00 | +0.00 | -0.00 | +0.00 | +0.00 | -0.00 | +0.00 |
| hf_dimer_FF:ref | +5.18 | +1.77 | -0.31 | -1.55 | -2.25 | -2.61 | -2.74 | -2.73 | -2.62 | -2.38 | -2.10 | -1.73 | -1.40 | -1.13 | -0.85 | -0.63 | -0.43 | -0.28 | -0.13 | +0.00 |
| hf_dimer_FF:gfnff | +1.93 | -1.66 | -3.95 | -5.30 | -5.92 | -6.05 | -5.89 | -5.54 | -5.10 | -4.41 | -3.75 | -2.99 | -2.38 | -1.88 | -1.38 | -0.99 | -0.64 | -0.38 | -0.16 | +0.00 |
| hf_dimer_FF:revgfnff | +2.88 | -1.66 | -3.95 | -5.30 | -5.92 | -6.05 | -5.89 | -5.54 | -5.10 | -4.41 | -3.75 | -2.99 | -2.38 | -1.88 | -1.38 | -0.99 | -0.64 | -0.38 | -0.16 | +0.00 |
| hf_dimer_FF:**rev-gfnff** | +0.94 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 |
| nh3_h2o_NO:ref | +25.63 | +16.35 | +10.03 | +5.72 | +2.81 | +0.86 | -0.40 | -1.19 | -1.64 | -1.92 | -1.91 | -1.71 | -1.42 | -1.15 | -0.85 | -0.62 | -0.40 | -0.25 | -0.11 | +0.00 |
| nh3_h2o_NO:gfnff | +20.17 | +12.15 | +7.26 | +3.37 | +1.13 | -0.19 | -3.70 | -4.09 | -4.16 | -3.81 | -3.20 | -2.59 | -2.04 | -1.58 | -1.12 | -0.77 | -0.46 | -0.25 | -0.10 | +0.00 |
| nh3_h2o_NO:revgfnff | +35.77 | +13.11 | +7.25 | +3.37 | +1.13 | -0.19 | -3.70 | -4.09 | -4.15 | -3.81 | -3.20 | -2.59 | -2.04 | -1.58 | -1.12 | -0.77 | -0.46 | -0.25 | -0.10 | +0.00 |
| nh3_h2o_NO:**rev-gfnff** | +15.59 | +0.96 | -0.01 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | -0.00 | +0.00 | -0.00 | -0.00 | +0.00 | -0.00 | +0.00 |
| water_dimer_OO:ref | +10.74 | +4.90 | +1.17 | -1.15 | -2.56 | -3.35 | -3.73 | -3.84 | -3.77 | -3.49 | -3.10 | -2.56 | -2.05 | -1.62 | -1.20 | -0.87 | -0.57 | -0.36 | -0.16 | +0.00 |
| water_dimer_OO:gfnff | +8.89 | +3.35 | -0.51 | -3.18 | -4.79 | -5.62 | -5.91 | -5.86 | -5.59 | -4.84 | -4.12 | -3.29 | -2.60 | -2.03 | -1.46 | -1.03 | -0.64 | -0.37 | -0.15 | +0.00 |
| water_dimer_OO:revgfnff | +14.33 | +3.38 | -0.50 | -3.18 | -4.79 | -5.61 | -5.91 | -5.86 | -5.58 | -4.84 | -4.12 | -3.29 | -2.60 | -2.03 | -1.46 | -1.03 | -0.64 | -0.37 | -0.15 | +0.00 |
| water_dimer_OO:**rev-gfnff** | +5.44 | +0.03 | +0.01 | +0.01 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 | -0.00 | +0.00 | -0.00 | +0.00 | +0.00 | -0.00 | +0.00 |

Full precision of every point, both `kept` and `fresh`, in the `--json` output (`test_cases/revgfnff/ref/_results/contact_curves.json`).
