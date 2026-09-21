# MOR41: curcuma CPU vs GPU vs gxtb 6.7.1 — gfnff / gfn1 / gfn2

> 🤖 AI-generated evaluation, ⚙️ machine-tested only. No human production testing.
> Reference = `gxtb` (xtb **6.7.1**, `--gfnff` / `--gfn1` / `--gfn2`). Device = NVIDIA
> RTX 5080 (CUDA, `-gpu cuda`). Set = all 95 MOR41 structures, charge 0, singlet.
> Runner: `scripts/mor41_cpu_gpu_gxtb_eval.py` (energy + full analytic gradient via the
> new `curcuma -sp ... -dump_gradient <file>`). Raw per-structure table:
> `test_cases/MOR41-testset/_run_cpu_gpu_gxtb/mor41_cpu_gpu_gxtb.md`.

## Summary matrix

| axis | gfnff | gfn1 | gfn2 |
|---|---|---|---|
| **Energy CPU=GPU** (MAD / max, kcal) | main-group ok; **metals 0.91 / 8.03** | 1.9e-9 / 4.5e-8 | 8.5e-9 / 1.3e-7 |
| **Energy vs gxtb** (MAD / max, kcal) | 11.6 / 178 (method divergence) | <0.1 (94/95; PR40 = gxtb fail) | **0.0003 / 0.002** |
| **Gradient CPU=GPU** (MAD, Eh/Bohr) | main-group 7e-5; **metals 3.4e-3** | 1.7e-5 | 1.5e-6 |
| **Gradient vs gxtb** (abs MAD, Eh/Bohr) | main-group 0.2%; metals = method div. | main-group perfect; metals ~1% | **≤0.05%** (F5 fixed: variational D4 q-response) |

`gxtb` uses the **attached** selectors `--gfn1` / `--gfn2`; the spaced form `--gfn 1`
is silently ignored by this build and defaults to GFN2 (verified via the printed
`Hamiltonian GFNx-xTB` line). The old `scripts/mor41_validation.py` used the spaced
form against xtb 6.6.1 — not an issue there, but noted for future runs.

## Findings

### F1 — GPU gfnff diverged from CPU for transition-metal complexes — ✅ FIXED (Jul 2026)
**Before the fix:** CPU and GPU gfnff agreed for main-group systems (energy MAD 0.015
kcal) but **not** for metals: energy MAD 0.909, max **8.03** kcal, 45/70 metal
structures >0.1 kcal; gradient max-component MAD 3.4e-3 Eh/Bohr. Every large outlier was
a TM complex: ED29 (Pt) +8.03, PR08 (Ir) −6.18, PR22 (Ir) −5.27, ED07 (W) +4.67, …

**Root cause:** the GPU CN used a hand-copied covalent-radii table (`s_rcov_d3_bohr[]`
in `cuda/ff_workspace_gpu.cu` and `rocm/gfnff_rocm.hip`) that held the **Grimme D3**
radii, whereas GFN-FF (and the CPU, via `CNCalculator::COVALENT_RADII` ==
`GFNFFParameters::covalent_rad_d3`) uses its own **modified** covalent radii. The two are
identical for Z≤36 but differ for 47 heavier elements — e.g. Pt 1.38 vs 1.12 Å, Ir 1.38
vs 1.11, W 1.41 vs 1.23. The larger GPU radii inflated the metal CN (Pt: GPU 4.34 vs CPU
4.09), which feeds the dynamic bond r0 = (r0_base + cnfak·CN)·ff (and the CN-dependent
angle term), so the metal bond/angle energy diverged. Main-group matched because the
radii there are identical. The outlier magnitude tracked the per-element radius error
(Ir/Pt ~0.26 Å → biggest; Ru ~0.07 → small).

**Fix:** upload the GPU CN/damping radii directly from the single source of truth
`GFNFFParameters::covalent_rad_d3` (already Bohr) instead of the drifted local copy, in
both the CUDA and ROCm backends. The hand-copied array (and the stale "matches
s_rcov_d3_bohr" comment in the CPU ATM) is removed. **After:** gfnff CPU-vs-GPU energy
MAD **3.3e-6**, max **1.5e-5** kcal (0/95 above 0.01 kcal), gradient MAD ~1e-6 Eh/Bohr —
i.e. GPU == CPU to FP-reduction noise for all 95 structures, metals included. Main-group
byte-unchanged; 19/19 `gfnff_gpu` + 35/35 golden-value ctests pass.

### F2 — gfn1 / gfn2 CPU = GPU to machine precision
Energy MAD 1.9e-9 (gfn1) / 8.5e-9 (gfn2) kcal, max ≤1.3e-7 kcal; gradients agree to
~1e-5 / 1e-6 Eh/Bohr. The native xTB CUDA path is internally consistent for both
methods across all 95 structures, metals included.

### F3 — Energy fidelity vs gxtb 6.7.1
- **gfn2: essentially exact** — MAD **0.0003**, max **0.002** kcal; 95/95 within 0.1
  kcal. Native GFN2 reproduces the xtb 6.7.1 reference across the whole metal set.
- **gfn1: 94/95 within 0.1 kcal.** The one outlier, **PR40 (Ti)**, is a **gxtb GFN1
  SCF failure**, not a curcuma bug: gxtb returns E = −93.61 Eh with gradient norm 17.5
  Eh/Bohr (unphysical) where curcuma gives −24.92 Eh (consistent with its own GFN2
  −25.07 Eh and a 0.09 Eh/Bohr gradient). Reference-side divergence.
- **gfnff: MAD 11.6, max 178 kcal — the known method divergence.** curcuma tracks
  the pprcht/gfnff port; pprcht and xtb 6.7.1 themselves disagree on TMs (MAD 11.6,
  max 178 kcal — same PR10/PR04/ED10 outliers). Main-group gfnff matches gxtb (75/95
  within 0.1 kcal). This is not a curcuma defect; see Known Issue #6.

### F4 — `Gradient()` unit convention differs between backends (open, needs operator validation)
`EnergyCalculator::Gradient()` returns **Eh/Å** for gfn1/gfn2 (a deliberate
`m_gradient /= au` in `xtb_native.cpp:1416`) but **Eh/Bohr** for gfnff
(`gfnff_method.cpp:1874`, "No conversion needed"). Each is self-consistent with its own
energy (finite difference: gfnff ratio to its own dE/dBohr = 1.000; gfn1 ratio to its
own dE/dÅ = 1.000), so the two backends carry **two contradictory design decisions**.

The consumers do **not** reconcile them: `SimpleMD` sets
`m_eigen_gradient = interface.Gradient()` directly (`simplemd.cpp:3869/3886`) with the
geometry in **Å** and no per-method conversion, and the optimizer likewise takes
`Gradient()` as-is. For an Å-space integrator the dimensionally-consistent force is
Eh/Å, so **gfn2 (Eh/Å) is consistent with the MD/opt while gfnff (Eh/Bohr) is 1.88973×
too weak there** — i.e. gfnff MD forces (and the effective `gradient_threshold`, whose
help text says "Eh/Bohr") are off by that factor. It was not caught because (a) LBFGS
line search is scale-tolerant, so gfnff optimizations still reach the correct minimum
(only the convergence-norn meaning shifts), and (b) gfnff MD still produces plausible
trajectories. **This is a real but subtle inconsistency.** Fixing it means picking one
convention and (crucially) re-validating MD **energy conservation** for both methods —
not done here, because a blind 1.88973× force-scale change to the validated MD/opt paths
(gfnff MD is used inside ConfSearch) is exactly the kind of plausible-but-wrong change to
avoid without that validation. The evaluation script converts gfn1/gfn2 to Eh/Bohr
before comparing gradients to gxtb.

### F5 — gfn2 analytic gradient ~1.1–9% off vs gxtb — ✅ FIXED (Jul 2026): self-consistent-D4 q-response now variational
**Root cause: the GFN2 self-consistent-D4 charge response was handled by a separate,
approximate mechanism instead of variationally through the gradient charge-Pulay.**
GFN2 couples D4 into the SCF — the C6 are scaled by the Mulliken charges, so `E_D4 = E_D4(q)`
and `dE_D4/dq` is an extra atom-potential in the Fock. curcuma added it to `v_at` **during
the SCF** (`addDispersionPotential`, so the density is correctly D4-polarized) but **omitted
it from the gradient's `m_pot` rebuild**, then bolted the `dE_D4/dq·dq/dR` term on via a
separate single-shot **EEQ** response (`d4_charge_source="eeq"`, ∂q/∂x from the EEQ model, not
the Mulliken charges). That inconsistency left a deterministic, SCF-tolerance-independent
residual collinear with the bonds (cos +0.995 with the charge-Pulay/W terms).

tblite instead folds `dE_D4/dq` into `pot%vat` (`dispersion%get_potential`), so the standard
charge-Pulay `−P·(v_a+v_b)·dS` **plus** `W = Σ f·ε·C·Cᵀ` (whose ε already carry the in-SCF D4
potential) reproduce the full `dE_D4/dq·dq/dR` for free — exactly like the ES2/3rd-order/AES
monopole responses. **The fix mirrors this:** `addDispersionPotential` is now called in the
gradient `m_pot` rebuild and the separate eeq/CPSCF response is skipped
(`m_d4_variational_qresp`, default ON for the new default `d4_charge_source="mulliken"`).

Validation (curcuma analytic vs gxtb 6.7.1 `--gfn2 --grad`, |Δ|/|g| over all forces):

| set | eeq (old default) | mulliken/variational (new default) |
|---|---|---|
| ED39 (Rh), ED02 (Fe) | 1.15%, 1.68% | **0.020%, 0.003%** |
| 27 MOR41 ≤14 atoms (organic + TM) | 0.3–8.8% | **0.000–0.047%** |

Also FD-self-consistent vs curcuma's own tight-SCF energy derivative to ~1e-5 (the FD floor)
on ED39/ED03/ED02/PR03/ED01; GFN1 unchanged (no D4); **energy bit-identical** in both modes
(the fix is gradient-only). `d4_charge_source` now selects the gradient q-response:
`"mulliken"` (default, exact variational), `"eeq"` (legacy single-shot, ~1% on TM), `"cpscf"`
(explicit Z-vector solve — exact energy weighting but the response itself carries the same ~1%
because it re-derives what the charge-Pulay already gives; kept for validation only).

**Diagnostic path (two harness bugs corrected en route):**
1. A **stale CN cache** in the frozen-density FD harness (`computeCoordinationNumbers`
   memoises on `m_cn_cache`, invalidated only by `UpdateMolecule`; the FD set `m_geometry`
   directly) silently froze CN in every FD, making the earlier "band-term 0.04% error"
   an artifact. Fixed by force-invalidating the cache per FD step; the band term is then
   **exact (ratio 1.00000)**, and its overlap+shpoly vs H0/CN sub-parts each 1.00000.
2. The frozen-density audit **freezes the charges**, so it structurally cannot see the
   charge response — the earlier "residual isolated to the −Tr(W·dS) response" conclusion
   was wrong. A **live-q frozen-P electronic FD** proved every explicit electronic term
   (band + Coulomb + AES-direct + AP5b + charge-Pulay) exact, and a decomposition of the
   analytic gradient correlated the gxtb residual +0.995 with the charge-Pulay/W direction,
   pointing at the missing D4 charge-Pulay.

The per-term audit (`XTB::auditGfn2GradientTerms`, env `CURCUMA_MP_GRAD_AUDIT`, GFN2 only,
zero default overhead) remains as a regression tool; on ED39 every term is now 1.00000
(multipole, Coulomb, band, band:ovlp+shp, band:H0-CN) with D4 at 0.993 (the excluded direct
q-response of the frozen-density FD). `CURCUMA_D4_NONVARIATIONAL` forces the legacy path for
A/B comparison. **GPU note:** `addDispersionPotential` runs in the shared `m_pot` rebuild
before both `calculateGradient` and `calculateGradientGpu`, so the device charge-Pulay reads
the D4-inclusive `v_at` too — the fix propagates to CUDA/ROCm by construction (unvalidated
here; no active GPU build).

## What was tested / not tested
- **Tested:** single-point energy + analytic gradient, all 95 MOR41 structures, three
  methods, CPU + CUDA GPU + gxtb 6.7.1; gradient units cross-checked by finite
  difference and against the gxtb TM gradient file.
- **Not tested:** ROCm / Vulkan GPU backends (build has only CUDA on); optimization or
  MD trajectories (single points only); charged/open-shell states; systems outside the
  MOR41 element/size range; the root cause of F1 (GPU gfnff metal terms) and F5 (gfn2
  TM gradient term) — both are reported, not yet diagnosed in the source.
