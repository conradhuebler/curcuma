# Technical Debt — GFN-FF, Native GFN-xTB, QM Interfaces, EnergyCalculator

Audit of the force-field, native xTB, QM-method-interface, and method-dispatch layers.
Findings are from a code audit (file:line cited); the highest-impact correctness bugs were
spot-verified by re-reading the cited lines and are marked **verified**.

This document identifies and records debt only — no fixes are applied here. Severity:

- **high** — correctness bug, crash, UB, or silent wrong-science the user cannot see
- **medium** — maintainability / divergence / latent bug under non-default conditions
- **low** — naming, docs, dead code, style

Status labels follow CLAUDE.md. Everything here is 🤖 AI-generated audit output, not
human-reviewed. Before acting on any *correctness* finding, confirm against the reference
implementation (xtb / tblite) yourself — a plausible-looking bug can still be masked by an
upstream guard not shown in the snippet.

---

## 1. GFN-FF (`src/core/energy_calculators/ff_methods/`)

### God files / structure

| ID | Location | Sev | Finding | Verified |
|----|----------|-----|---------|---------|
| F-E1 | `ff_methods/gfnff_method.cpp` (whole, 11 738 lines) | high | Single `GFNFF` class holds topology, CN, EEQ dispatch, all param generation, `Calculation()` (~530 lines, 1686–2218), solvation, numeric gradients, diagnostics, JSON (de)serialization, 14+ accessors. Largest maintenance hazard in scope. Split by concern (`gfnff_topology/params/energy/diag`). | — |
| F-E3 | `gfnff_method.cpp:10277–10408` | medium | ~20 per-term accessors each fork `if (m_forcefield) …; if (m_workspace) …;`. Dual-backend fork smeared across API. Introduce one `GFNFFEngine` interface. | — |
| F-E2 | `gfnff_method.cpp:1849–1869` | medium | NaN-trap block calls `scanComp(...)` for 10 terms twice (workspace + forcefield). Lists already diverge in style. Iterate a `{name,getter}` table. | — |
| F-E10 | `gfnff_torsions.cpp` (2174 ln), `gfnff_inversions.cpp` (891 ln) | medium | Monolithic single-purpose TUs; torsion file embeds ring gating, pi-sp3 override, hyperconjugation, NCI in one function. Split. | — |
| F-E11 | `ff_methods/gfnff.h` (2634 ln) | medium | God header declaring the entire class + nested structs + GPU seams. Extract `gfnff_types.h`, forward-declare. | — |
| F-W1 | two files named `gfnff_method.cpp` | medium | `ff_methods/gfnff_method.cpp` (11 738 ln `GFNFF`) vs `qm_methods/gfnff_method.cpp` (201 ln adapter). Ambiguous grep. Rename adapter to `gfnff_computational_method.*`; delete stale TODO in its header. | — |

### `#ifdef` soup / magic numbers

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| F-E4 | `gfnff_method.cpp:1141,1163,1203,1300,5929,5957,6108,6121,6237,6421` | medium | `#ifdef GFNFF_CN_DCN_FUSION` / `GFNFF_FAST_EXP` interleave conditional compilation with the science in the CPU hot path (CMake options, `CMakeLists.txt:129,131`). Default-build reader cannot see control flow. Use a runtime strategy flag / template policy. |
| F-E5 | `gfnff_method.cpp:4307–4312,4553,4589,4666–4690,4752–4755,4926–4928,5735` | medium | Dozens of Fortran constants inlined as locals with `// Fortran: gen%…` comments, no central named definition. Collect into `struct GFNFFGenParams`. |
| F-E6 | `gfnff_method.cpp:11156–11234` & `11235–11309` | medium | `NumGrad` and `NumGradFixedCharges` near-identical ~80-line FD loops. One `numericGradient(dx, freeze_charges, recompute_cn)`. |

### Dead / stale code

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| ~~F-Q6~~ | ~~`qm_methods/gfnff_torsions_NEW.cpp`~~ | ~~high~~ | **RESOLVED (deleted 2026-06-26):** obsolete Dec-2025 draft, not in CMake, not `#include`d; lacked the June-2026 `!in_ring` ring-torsion fix and was an incomplete "Simplified/placeholder" port. Removed; canonical `ff_methods/gfnff_torsions.cpp` (2174 ln) is the strict superset. | — |
| F-Q1 | `gfnff_method.cpp:6907` (~250 ln) | medium | Legacy single-phase `calculateEEQCharges` behind unregistered JSON flag `use_two_phase_eeq` (default true at :8881). Divergent untested path. Delete or register as PARAM. |
| F-Q2 | `gfnff_method.cpp:11059` + `gfnff.h:2258` | low | Dead `calculateFinalCharges(TopologyInfo&,int,double)` wrapper, no callers, "DEPRECATED" comment next to it. Remove. |
| F-Q3 | `gfnff_advanced.cpp:28` | low | `GFNFFAdvanced::calculateEEQCharges` returns hardcoded H/C/N/O/F charges. Whole namespace self-referential, no callers. Drop from build (`CMakeLists.txt:504`). |
| F-E7 | `gfnff_method.cpp:2884` | low | Success log still says "Topology calculation complete (stub - always returns true)" though `calculateTopology()` is real. |
| F-E8 | `gfnff_method.cpp:3986,4063,4069,7181,7206,9363–9370,10015,10147,10224` | low | Mixed stale TODOs (some done, some open). The dxi/dgam "not implemented" at 9363 contradicts EEQ solver docs. Audit each. |
| F-W3 | `qm_methods/gfnff_method.cpp:93` | low | `getDipole()` returns `{0,0,0}` though `GFNFF` computes a dipole (`gfnff.h:436`). Forward. |
| F-W2 | `qm_methods/gfnff_method.h:84,58–60` | low | Polymorphic wrapper leaks raw `GFNFF*` + EEQ-internal PCG counter as test hooks. Move to a diagnostics interface. |

### Numerical / silent fallback

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| F-Q4 | `gfnff_method.cpp:8996,9086`; `eeq_solver.cpp:922` | high | EEQ solver silent uniform-charge fallback (`q_i = total/N`) + second `keep previous m_charges` layer; neither sets `m_has_error`. A wrong charge set propagates into Coulomb E/grad with no abort. Return a status; refuse to return E/grad when charges invalid. |
| F-Q7 | `qm_methods/gfnffinterface.cpp:304–307` (init `Vector::Zero` :97) | high | External `GFNFFInterface::Charges()` always returns zeros — the C singlepoint API has no charge-output arg, so `m_charges` is never populated. D4/analysis consumers silently get zero charges. Return empty / "unavailable" flag. | — |
| F-Q8 | `qm_methods/gfnffinterface.cpp:243–246` | medium | On `iostat != 0` logs error but returns `0.0` (success-looking) energy, no `m_has_error`. |
| F-G2 | `gfnff_method.cpp` repulsion grad | medium | `ff_methods/CLAUDE.md` table marks bonded repulsion gradient "⚠️ Partial" but `Calculation()` returns a gradient including repulsion with no runtime guard. `-opt` users get a knowingly-incomplete gradient silently. Finish or one-time-warn. |
| F-E9 | `gfnff_method.cpp:4761,4776,5025` | low | `static bool/int` debug counters shared across instances/threads. Race under molecule-level CxxThreadPool parallelism. Remove / promote to members. |
| F-Q5 | `gfnff_method.cpp:9136–9141` | medium | `std::ofstream("gfnff_diag_charges.json")` to CWD at verbosity 3 — BMT/outputPath violation (CLAUDE.md mandate). Route through `BMTUtils::outputPath()`. |

---

## 2. Native GFN-xTB (`src/core/energy_calculators/qm_methods/xtb_*`)

### SCF

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| X-S1 | `xtb_scf.cpp:543–547` | high | `q = Z - n_at` computed then immediately overwritten by `q = n0_at - n_at` (correct GFN convention). Dead line + stale comment = reorder hazard. **Delete line 543.** |
| X-S2 | `xtb_scf.cpp:664–672` vs `xtb_native.cpp:188–195` | high | Two divergence checks with 100× different thresholds (`de < threshold` vs `de < thresh*100`). The SCF loop uses the ×100 one; the stricter appears dead. One wiring mistake = 100× convergence change. Delete or unify with factor documented. |
| X-S3 | `xtb_scf.cpp:433–474` vs `594–628` | medium | Fermi-Dirac bisection + occupation assembly duplicated (`:592` comment "Keep in sync"). `solveEigen` should call `occupationsFromEps`. |
| X-S4 | `xtb_native.cpp:877,959,1082,1165` | medium | FP32 mixed-precision SCF never guarantees an FP64 landing; if `m_scf_fp32_threshold` too low it runs FP32 to `max_iter` with no polish/warn. Force a final FP64 step or warn. |
| X-S5 | `xtb_native.cpp:894–901` vs `1165–1172` | medium | Two separate convergence sites in the SCF loop (resident vs CPU/GPU); risk drifting apart. One `checkScfConvergence(...)`. |
| X-S6 | `xtb_native.cpp:1255–1263` then `1313–1322` | medium | Post-SCF potential rebuilt twice (for `m_F` then for gradient), each incl. expensive `addDispersionPotential` (D4 dE/dq). GFN2 = two D4 evals post-SCF. Build once, pass to both. |
| X-S7 | `xtb_native.cpp:1044–1075` vs `1178–1182` | low | Broyden `case` empty; actual mix 100 lines later. Obscures the science. Add pointer comment. |
| X-S8 | `xtb_scf.cpp:231,257` | low | `static bool warned_T/warned_conv` function-static mutable (latent race). `thread_local` or members. |
| X-S9 | `xtb_native.cpp:1416–1427` vs `375–454` | low | `evaluateComponentsAtFixedDensity` duplicates `Calculation()` pre-SCF setup. Extract `prepareGeometryDependentState()`. |

### Integrals (H0 / CN / gamma / multipole)

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| X-I1 | `xtb_h0.cpp:42–50` (`ao_to_type` returns -1 for ang==2), `308–313` (`if (t_a < 0) continue`) | high | d-shell AO pairs silently skipped → zeroed in S/H0. GFN2 param tables include d shells for some elements (S, Cl, metals); the H/C/N/O test set does not exercise it. CLAUDE.md doesn't document the limit. **Confirm whether the basis emits d shells**; if yes implement, if no assert + document. | — |
| X-I2 | `xtb_h0.cpp:410–416` | medium | GFN1 halogen-bond `calcHalogenBondEnergy()` returns 0.0 ("to be implemented if accuracy requires it") yet `m_E_halogen_bond` is added to the total (`xtb_native.cpp:1278`) as silent zero. Implement or document as not-implemented. |
| X-I3 | `xtb_native.cpp:375`, `xtb_multipole.cpp:167`, `xtb_response.cpp:353` | medium | Coordination numbers recomputed 3× per `Calculation`. Cache as a member after first compute; fix sequencing. |
| X-I4 | `xtb_h0.cpp:280–301` & `xtb_response.cpp:410–432` | medium | GFN1/GFN2 H0 hscale logic duplicated verbatim between H0 build and response. Factor `hscale(...)`. |
| X-I5 | `xtb_h0.cpp:27–50`, `xtb_multipole.cpp:33–54`, `xtb_response.cpp:33–48` | low | `as_cgto_shell`/`ao_to_type` copied 3×; valence-flag logic copied 3×. Move to shared `xtb_ao_utils.hpp`. |
| X-I6 | `xtb_h0.cpp:241–333`, `xtb_coulomb.cpp:21–40` | medium | No spatial cutoff in H0/overlap or gamma build (dense O(N²)/O(nsh²)). `large_system_mode` addresses SCF scaling, not integral build. >1000 atoms → H0 bottleneck. Optional cutoff; document dense limit. |
| X-I7 | `xtb_native.cpp:1872–1877,2120–2121,2129–2132` (+ ATM 1995/2184) | medium | GFN2 D4 BJ params (s6=1.0,s8=2.7,a1=0.52,a2=5.0,alp=16.0,s9=5.0) hardcoded 3×. One `gfn2D4Params()` struct. |
| X-I8 | `xtb_native.h:1098` (`m_X`), `1083` (`m_h0`) vs `1089` (`m_H0`) | low | `m_X` is the Cholesky L now (was Löwdin S⁻¹ᐟ²) — mis-named. `m_h0`/`m_H0` case-only distinction. Rename `m_cholS`, `m_h0_params`. |
| X-I9 | `xtb_native.h:936,942` | low | Header file-location comments wrong (`buildH0Data` actually in `xtb_native.cpp:1578`; `computeCoordinationNumbers` in `xtb_h0.cpp:55`). |
| X-I10 | `xtb_multipole.cpp:298 & 348` (mpscale_q), `207–208,178–181` | low | `mpscale_q` duplicated; magic damping constants undocumented. One `constexpr` array with citation. |

### Third-order / response

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| X-R1 | `xtb_response.cpp:287–676` | high | `computeMullikenChargeResponse` ~390-line god function: property gradient, Z-vector solve, relaxed density, H0-diagonal CN coupling, off-site Pulay+shpoly, multipole Pulay (gated off), Coulomb ES2 response, multipole-interaction (gated off), CN chain-rule. Decompose. |
| X-R2 | `xtb_response.cpp:312` (`constexpr bool MP_RESPONSE_ENABLED = false`), `:470,:575` | medium | Full GFN2 multipole charge-response compiled but disabled ("correct for low-multipole, mis-contracted for polar"). D4 gradient target met without it, but the full reference response is silently omitted. **Document as deferred limitation** in CLAUDE.md, not silent. |
| X-R3 | `xtb_response.cpp:295` vs `xtb_scf.cpp:470`/`xtb_native.cpp:827` | low | `std::round` vs `std::floor` for nocc. Latent for odd electron counts. One convention. |
| X-R4 | `xtb_response.cpp:233–262` | low | Z-vector CG solver: no convergence warning, hardcoded limits (max_iter=200, tol=1e-10). Silent unconverged return. Warn; make limits settable. |

### Threading

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| X-T1 | `xtb_dc.cpp:83`, `xtb_scf.cpp:231,257` | medium | `static bool` one-shot warning flags in hot functions shared across instances/threads. Safe only because calls are currently serial. `atomic` or instance members. |
| X-T2 | `xtb_fragment_scf.cpp:193` | low | Fragment loop serial (only intra-fragment MKL threads). Many-fragment solvent boxes leave parallelism unused. Parallelise or document. |
| X-T3 | `xtb_gradient.cpp:225–448` | low | Threaded gradient reduction not bit-identical across thread counts (FP reassociation). CLAUDE.md "bit-identical t1/t8" claim is for the serial path. Reconcile in SQM_THREADING.md. |
| X-T4 | `xtb_native.cpp:70–105` | low | Threading gating contract clear in code but spread over 3 files; SQM_THREADING.md describes perf not the contract. Add a block comment at `effectiveIntraThreads`. |

### Large-system modes (fragments / dc / sparse)

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| X-L3 | `xtb_sparse.cpp:195` (`m_S.inverse()`), `264` (`Heff = Sinv*F`) | high | Sparse mode computes a dense full S⁻¹ (O(N³)) + dense GEMM every iteration. In-code comment admits "S⁻¹-free O(N) route is deferred." SQM_LARGE_SYSTEMS.md markets sparse as the O(N) path — **it is not**. Prominently document as a measurement harness, or implement purification. | **verified** |
| X-L4 | `xtb_dc.cpp:345–348`, `xtb_sparse.cpp:398–400` (consumed `native_xtb_method.cpp:307–320`) | high | DC/sparse return **zero gradient** with only a warning, no `m_has_error`. `-opt`/`-md` with `-large_system_mode=dc|sparse` silently does nothing (stuck opt, no signal). Set `m_has_error` and refuse. |
| X-L1 | `xtb_dc.cpp:171–184`, `xtb_sparse.cpp:177–192` vs `xtb_native.cpp:374–184` | high | DC + sparse near-verbatim copy the global setup (CN→self-energies→H0/S→gamma→multipole→CN-assign incl. cache resets). If setup order/invalidation changes in `Calculation`, both copies silently diverge. Extract `buildGlobalSetup()`. |
| X-L2 | `xtb_dc.cpp:230–274`, `xtb_sparse.cpp:198–227` vs `xtb_native.cpp:269–302` | high | `packSCC`/`unpackSCC` duplicated 3× (SCC-vector pack/unpack incl. GFN2 multipoles). Layout change → silent desync of extrapolation history / DC / sparse mixing. Call existing `packSccState`/`unpackSccState` members. |
| X-L5 | `xtb_fragment_scf.cpp:173–178` | medium | Every fragment assigned `charge=0` regardless of `m_mol.m_charge`. Charged clusters get wrong total energy, no error flag. Refuse or partition charge. |
| X-L6 | `xtb_dc.cpp:146–163` | medium | DC core overlap (last-writer-wins) and uncovered atoms dropped from density assembly — both warn but don't abort; SCF proceeds on a wrong density. Hard error. |
| X-L7 | `xtb_fragment_scf.cpp:278` | medium | `au_per_bohr = 0.52917721092` mis-named (it is Å/Bohr); hardcoded factor violates CurcumaUnit mandate. Rename, use `CurcumaUnit`. |
| X-L8 | `xtb_dc.cpp:280`, `xtb_sparse.cpp:239` | low | `keep_diis` silently ignored in DC/sparse (fresh `BroydenMixer` per call). Document or honour with a member mixer. |
| X-L9 | `xtb_dc.cpp:245`, `xtb_sparse.cpp:229` | low | Hardcoded `3.166808e-6` K→Hartree. Use `CurcumaUnit::boltzmann_hartree_per_K`. |

### SCF extrapolation (ASPC / Gauss / XL-BOMD)

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| X-X1 | `xtb_native.cpp:502` | medium | XL-BOMD `K` clamped `min(max(order,3),7)` with no warning; predictor still uses the requested order → predictor/integrator at different orders. Warn or reject at `setMolecule`. |
| X-X2 | `xtb_native.cpp:1196–1222` | medium | On non-convergence the XL-BOMD auxiliary is not advanced → next guess is stale; persistent divergence never recovers. Re-seed from best-available density. |
| X-X3 | `xtb_native.cpp:135–136` vs `xtb_native.h:1245` | low | History/aux cleared on new molecule but `m_xlbomd_warned` not reset. |
| X-X4 | `xtb_native.cpp:547–572` | low | Extrapolation prediction has no physical-range guard (only `allFinite()`). Verify charge-sum/per-element bounds; else fall back to 1-step warm-start. |
| X-X5 | `xtb_native.cpp:1198` | low | XL-BOMD silently inert when `packSccState` returns short vector (GFN2 bootstrap). Warn once. |
| X-X6 | `xtb_native.cpp:332–342` | low | `xlbomdCoefficients` dissipation invariant `Σ c_k = 0` not checked. `assert`. |
| X-X7 | `xtb_native.cpp:1716–1718` | low | Gauss `maxlen` capped at 10 with no rationale; for `order=8` the Vandermonde fit is near-degenerate. Document why 10. |

### Gradients

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| X-G2 | `xtb_gradient.cpp:777–933` vs `551–677` & `716–765` | high | `calculateGradientGpu` duplicates ~160 lines of section 5 + section 4 ("identical to section 5/4"). Any multipole-damping / traceless / CN fix must be applied twice. Extract `addMultipoleHostGradient()` / `addCNChainRuleGradient()`. |
| X-G1 | `xtb_gradient.cpp:38–58` vs `xtb_h0.cpp:27–50` | medium | `as_cgto_shell_g`/`ao_to_type_g` duplicate H0 versions with `_g` suffix (comment admits it). Shared `xtb_cgto_utils.h`. |
| X-G3 | `xtb_gradient.cpp:94,99–102` | medium | Gradient assumes closed-shell integer occupation (`W = Cocc·diag(2·eps)·Coccᵀ`). For `electronic_temperature > 0` (Fermi smearing) density is fractional but W uses integer columns → Pulay wrong. Document 0-K only, or build W with smeared `f_i`. |
| X-G4 | `xtb_gradient.cpp:57,316` | medium | d-type AOs contribute zero gradient (`ao_to_type_g` returns -1 for ang≥2). Latent (s/p only today). `assert(ang<=1)` or implement. |
| X-G5 | `xtb_gradient.cpp:16–18` | low | Stale TODO "AP5: integral Pulay term" — implemented at 342–414 (AP5b). Update header. |
| X-G6 | `xtb_gradient.cpp:679–701` | low | D4 gradient comment says q-response "not yet implemented" — CLAUDE.md says X-AP4 made it analytic (2026-06-18). Contradicts current state. |
| X-G7 | `xtb_gradient.cpp:69` | medium | Section numbering 1,2a,2b,3,5,3b,4 ("5 must run before 4") obscures the dEdcn dependency. Renumber monotonically. |
| X-G8 | `xtb_gradient.cpp:389` | low | Unicode identifier `Δqraw`. Rename `dqraw`. |

### Wrapper / method abstraction

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| X-M1 | `xtb_method.h:115–157,163–176,209–242` | high | ~15 methods declared in `xtb_method.h` with **no definition** in `xtb_method.cpp` (which only implements the interface + `saveToFile` stub + `getSupportedMethods`). Any caller → link error. Largest dead-code surface in scope. Implement or remove declarations. |
| X-M2 | `xtb_method.cpp:97–113` vs `native_xtb_method.cpp:383–389` | medium | `getEnergyDecomposition` schema diverges (all-zero FF-style placeholder vs real QM decomposition). Backend switch → different schema, no signal. Agree on one QM schema. |
| X-M3 | four "XTB" types | medium | `XTBMethod` (external) vs `NativeXtbMethod` (native wrapper) vs `curcuma::xtb::XTB` (engine) vs `XTBInterface` (external). Conceptual collision. Rename engine → `NativeXtbEngine`. |
| X-M4 | `native_xtb_method.h:82` `solver()` | low | Raw native-engine pointer escape hatch (GPU eigensolver seam). Expose `setExternalEigensolver(...)` forwarding instead. |
| X-M5 | `native_xtb_method.h:47` | low | `getDipole()` stub returns zero though charges exist. Implement or drop override. |
| X-M6 | `xtb_method.cpp:88` | low | `isThreadSafe()==true` unsubstantiated (libxtb may have module globals). Verify/cite or return false. |
| X-M7 | `native_xtb_method.cpp:197` | low | Explicit `QMInterface::` scope bypasses virtual dispatch. Comment why. |
| X-M8 | `xtb_fragment_scf.cpp:35–42,54–58`, `native_xtb_method.cpp:113–134` | medium | Scoped-config fallback pattern hand-written 5×. Expose `getScoped(key, fallback)`. |

---

## 3. QM method interfaces (`src/core/energy_calculators/qm_methods/`)

### Confirmed correctness bugs (act-first)

| ID | Location | Sev | Finding | Verified |
|----|----------|-----|---------|---------|
| Q-1 | `xtb_method.cpp:11–19` | high | `XTBMethod` ctor never calls `m_xtb->setMethod(method_name)`. `XTBInterface` ctor (`xtbinterface.cpp:36–58`) reads accuracy/maxiter/temp/spin from config but **not** `method`, and never calls `setMethod`; `m_method_switch` defaults 0 = **GFN0**. So `xtb-gfn2`/`xtb-gfn1` via the XTB-binary backend (builds with `USE_XTB=ON`, `USE_TBLITE=OFF`) silently run GFN0-xTB. TBLite wrapper correctly calls `setMethod` (`tblite_method.cpp:37`). **One-line fix:** add `m_xtb->setMethod(method_name);` in the `#ifdef USE_XTB` block. | **verified** |
| Q-2 | `orcainterface.cpp:287` | high | ORCA multiplicity `2.0*mol.m_spin + 1.0 + 0.5`. `m_spin = multiplicity - 1 = 2S` (`abstract_interface.h:138` `m_spin = multi - 1`), so this yields `4S+1`: doublet (multi=2) → 3 (triplet), triplet → 5 (quintet). All open-shell ORCA inputs get the wrong multiplicity. **Fix:** `multiplicity = static_cast<int>(mol.m_spin) + 1;` (and drop the `+0.5`). | **verified** |
| Q-3 | `tbliteinterface.cpp:291–292` | high | `CurcumaLogger::error("GB solvation model not functional"); exit(1);` — hard-exits the **entire process** from inside a library when GB solvation is requested; line 293 is dead code. Any caller requesting GB via TBLite kills the host. Throw `std::runtime_error` instead. | **verified** |
| Q-4 | `dftd3interface.cpp:163–167` | high | `InitialiseMolecule(vector<int>)` fills `m_coord[3*i+0] = (3*i+0)/au` — divides the **index** by `au`, not a coordinate. Geometry is `0/au, 1/au, 2/au, …`. Any caller of this overload computes dispersion at garbage geometry. Take coords as a parameter or read from `m_geometry` after base init. | **verified** |
| Q-5 | `dftd3interface.cpp:72–73` | high | Destructor `delete m_coord; delete m_attyp;` — both allocated with `new[]` (lines 161–162). `delete` on `new[]` is **undefined behavior**. Use `delete[]`. | **verified** |
| Q-6 | `dftd4interface.cpp:106` | high | `UpdateAtom(i, x*factor, y*factor, z, element)` — z is **not** scaled by `factor` while x/y are. For `factor != 1` the geometry is sheared. Default `factor=1` masks it. Apply `* factor` to all three (or none) and document. | **verified** |
| Q-7 | `dftd3interface.cpp:86–96` | high | `if (m_damping.compare("bj"))` — `std::string::compare` returns nonzero on **mismatch**, so the rational-damping branch fires for every string *except* "bj"; the custom-param path loads the wrong damping model (the functional-load path at 100–110 correctly uses `== 0`). Use `if (m_damping == "bj")`. | **verified** |
| Q-8 | `tbliteinterface.cpp:166–186` | high | Re-init with a different `natoms` overwrites `m_coord[3*i+...]` without realloc → heap overflow if new natoms exceeds the old allocation. `delete[]`/`new[]` when `natoms != m_atomcount`. | — |
| Q-9 | `xtbinterface.cpp:135–147` | high | `xtb_newMolecule` error-check block is **commented out**; failure silently returns `m_initialised=true`. Uncomment and check. | — |
| Q-10 | `tblite_method.cpp:22–23,27,36,38` + `ulysses_method.cpp:34–35` | medium | Constructors call `CurcumaLogger::info/param` unconditionally — violates Level 0 silent mode (breaks opt/MD the moment a TBLite/Ulysses method is constructed). Gate behind `if (CurcumaLogger::get_verbosity() >= 2)`. | — |
| Q-11 | `ulyssesinterface.cpp:160–165` | medium | HOMO index = `m_atomcount / 2` — uses **atom count** as electron count. Reported HOMO/LUMO/gap at verbosity ≥ 2 is garbage for anything but neutral all-H. (Energy unaffected.) | — |
| ~~Q-31~~ | ~~`alpb_solvation.cpp:24` + `pch_external.h`~~ | ~~high (build break)~~ | **RESOLVED (2026-06-26):** removed the `dftd_*` (s-dftd3 / cpp-d4) includes from `pch_external.h` — `dftd_econv.h` defined global `kcaltoau`/`aatoau` that collided with `using namespace ALPBParameters;` in `alpb_solvation.cpp` (gcc-16-exposed). The D3/D4 interface TUs include those headers directly (USE_D3/USE_D4-gated), so only precompilation was lost. `curcuma_core` builds clean. Production `alpb_solvation.cpp` untouched. | v |

### Leaks / UB / resource handling

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| Q-12 | `xtbinterface.h:128–129` & dtor `xtbinterface.cpp:64–65` | medium | `double* m_coord; int* m_attyp;` not initialised to nullptr; dtor `delete[]`s them unconditionally — UB if `InitialiseMolecule` never called. Initialise to `nullptr`. |
| Q-13 | `xtbinterface.cpp:68–106` | medium | `InitialiseMolecule(Mol&)` / `(Mol*)` `new[]` `m_coord`/`m_attyp` without freeing the previous buffers when re-initialising (early-return on `m_initialised` then unconditionally `new`). Leak. `delete[]` before realloc. |
| Q-14 | `tbliteinterface.h:104` / `tbliteinterface.cpp:56–57,114–125` | medium | `char* m_solvent` allocated with `new[]` in ctor, never `delete[]`d in dtor. Leak one string per instance. |
| Q-15 | `tbliteinterface.cpp:337,424,564` | medium | `int count = 0;` never incremented; `if (count == 1) { tblite_delete_container(&m_tb_cont); }` unreachable. Container-delete path broken. |
| Q-16 | `xtbinterface.cpp:265` | medium | On XTB failure returns `return 4;` as a double energy — a caller doing `if (energy < 0)` treats 4 Eh as valid. Return 0/NaN + set `m_has_error`. |
| Q-17 | `orcainterface.cpp:163–169` | medium | `executeOrcaProcess` builds `"orca " + path + " > " + out` → `std::system()`. Unescaped paths = shell injection; hardcodes `"orca"`. Marked `@deprecated` but still public/used. Route through `runProcessAndCapture`. |
| Q-18 | `orcainterface.cpp:105–128` | medium | All ORCA files written to CWD (`m_input_path = basename + ".inp"`) — BMT/outputPath violation. No `-no_bmt` handling. Accept an output dir from the wrapper. |
| Q-19 | `orca_method.cpp:149–160` | medium | `setThreadCount` re-creates `m_orca` but `orca_basename` unchanged → all pool threads write the same `orca_calc.inp/.out`. No basename-unique helper. |
| Q-20 | `orcainterface.cpp:750–769` | medium | `calculate()` read-modify-write of input file with no locking — shared-basename instances corrupt each other. |
| Q-21 | `orcainterface.cpp:780–801` | low | natoms from re-parsing `*xyz` block via `line.find("*")` — any `*` in a comment ends the block early → wrong-sized gradient/charge. |
| Q-22 | `orcainterface.cpp:153–161` | medium | German `std::cerr << "Fehler beim Erstellen der Eingabedatei!"` — non-English, bypasses CurcumaLogger, violates the no-raw-std::cout/cerr rule. |
| Q-23 | `orcainterface.cpp:499–528` | low | `parseEnergy/Gradient/Charges/Dipole` each `catch(...)` and silently fall through — JSON schema mismatches invisible. |
| Q-24 | `dftd4interface.cpp:130–132` | medium | `clear()` is a no-op; does not free `m_mol` (only dtor calls `m_mol.FreeMemory()`). `clear()` then re-init leaks. |
| Q-25 | `dftd4interface.cpp:56–59` / `dftd3interface.h:56–62` | low | Default ctors skip `CreateParameter()`/`d4par()` → `m_param`/`m_par` uninitialised. `Calculation` on a default-constructed instance = UB. Delegate to the ConfigManager ctor. |
| Q-26 | `dftd4interface.cpp:120–128` | medium | `Calculation(bool gradient)` always computes gradient regardless of flag; if `m_gradient` uninitialised → UB. No error reporting. |
| Q-27 | `dftd3interface.cpp:192–204` | medium | Gradient flag ignored (always computed); if `m_gradient` never sized (raw-init path doesn't call base) → buffer overflow. Size in every init path, gate on flag. |
| Q-28 | `dftd3interface.cpp:184–190` vs `177–182` | medium | `UpdateGeometry(double*)` memcpy "assumed Bohr" while `UpdateAtom` divides by `au` (Å→Bohr). Opposite unit assumptions across update paths. Pick one. |
| Q-29 | `dispersion_method.cpp:30–37` | medium | `setMolecule` returns false (D3/D4 not compiled) but wrapper sets `m_initialized=true` regardless. Inconsistent. |

### Interface / abstraction / duplication

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| Q-30 | `abstract_interface.h:147` | medium | `double m_charge, m_spin = 0, m_muli = 1;` — `m_charge` left **uninitialised**. Any read before `InitialiseMolecule` = garbage. `double m_charge = 0.0, …`. |
| Q-31 | `abstract_interface.h:30–66` | medium | Four `InitialiseMolecule` + three `UpdateMolecule` overloads each copy the same 5-field bootstrap. Factor `assignFromMol(...)`. |
| Q-32 | `abstract_interface.h:127–132` | medium | Default `Charges/Dipole/BondOrders/OrbitalEnergies/OrbitalOccupations` return empty `Vector{}` — caller cannot distinguish "unsupported" from "not yet computed". Document or return optional/sentinel. |
| Q-33 | `abstract_interface.h:122` | low | `virtual bool Error()` exists; no subclass `Calculation()` consults it. Dead contract. |
| Q-34 | `computational_method.h:241` | medium | `getEnergyDecomposition()` pure-virtual → every QM wrapper returns a hardcoded zero JSON `{"Bond":0.0,…,"BATM":0.0}` (duplicated across EHT/XTB/TBLite/Ulysses/Dispersion). Misleading: callers may believe components were computed. Make it a default `json::object()`, override only in FF. |
| Q-35 | five wrappers `*_method.cpp` | medium | `getEnergyDecomposition` boilerplate + `getCharges/getBondOrders/getDipole` `#ifdef … return inner->X(); #else return Zero; #endif` pattern repeated ~5× (~150 lines). A `QMInterfaceMethod` base or templated adapter. |
| Q-36 | wrappers `hasError/clearError/getErrorMessage` | medium | No-op error plumbing: wrappers manage `m_has_error`/`m_error_message` but none propagate from the underlying interface (e.g. `XTBInterface::Calculation` returning 4 never sets `m_has_error`). Decorative. |
| Q-37 | wrappers `setParameters` | low | All do `m_parameters = params;` without re-configuring the wrapped interface. PARAM-registry path bypassed after construction → stale config. |
| Q-38 | `xtb_method.h:25` | medium | Stale include path `src/core/qm_methods/xtbinterface.h` (pre-restructure). Works only via include-path fallback. Use `xtbinterface.h`. |
| Q-39 | `tbliteinterface.cpp:615–621` vs `xtbinterface` | low | `BondOrders` returns flat NxN (TBLite) vs length-N (XTB) — two conventions for the same property. Standardise. |
| Q-40 | `tbliteinterface.h:121–122` | low | `TBLiteMethod` enum **and** `int m_method_switch` both encode the method, kept in lockstep manually in `setMethod`. Two sources of truth. |
| Q-41 | `tbliteinterface.cpp:417–488` | medium | SCF-error recovery tears down + rebuilds calculator by duplicating ~30 lines of setup (360–405). Factor `configureCalculator()`. |
| Q-42 | `tbliteinterface.cpp:175–180` | low | Raw-pointer `InitialiseMolecule` overload says "coord already in Bohr — no conversion" while every other overload divides by `au`. Intentional but fragile/undocumented. |
| Q-43 | `ulyssesinterface.h:55` | medium | Only overrides one `UpdateMolecule`; base `InitialiseMolecule(Mol&)` never calls `m_ulysses->setMolecule(...)` → base-overload callers get an unconfigured Ulysses. Override the `Mol&` overload or document the entry point. |
| Q-44 | `ulyssesinterface.h:72` | low | 26-element `m_solvents` literal in header — any solvent addition recompiles every consumer. Move to a params file. |
| Q-45 | `ulyssesinterface.cpp:144` | low | TODO "Capture and filter Ulysses output based on verbosity" — Ulysses raw stdout unfiltered at any verbosity. |
| Q-46 | `ulysses_method.cpp:98–110` | low | 11-branch if/else method→citation on every `calculateEnergy()`. Static `unordered_map`. |
| Q-47 | `ulyssesinterface.cpp:108` | low | `"C1"` point group hardcoded in `setMolecule`. Document or expose. |
| Q-48 | `eht.cpp:239` | low | `eV2Eh` name ambiguous (it is Eh/eV = 27.211, used as `/ eV2Eh`). Name `eV_per_Eh` or use `CurcumaUnit`. |
| Q-49 | `eht_method.cpp:127` | medium | `fmt::print("Warning: EHT does not provide analytical gradients")` bypasses CurcumaLogger, prints at Level 0 (violates silent mode). Use `CurcumaLogger::warn()`. |
| Q-50 | `eht.cpp:219–234` | low | Eigen/overlap failure returns `0.0` with `m_has_error` never set. Distinguish "energy zero" from "failed". |
| Q-51 | `eht.cpp:270–272` | low | `MakeBasis` stores Å, `MakeOverlap` divides by 0.529177 → Bohr. Two-step unit conversion hidden across functions. |
| Q-52 | `eht_method.cpp:160–175` | medium | `getCharges` returns `Zero(natoms)`, `getBondOrders` empty, `getDipole` zeros — different semantics than the underlying `EHT` (which returns empty). Pick one. |
| Q-53 | `eht.h:23` vs `eht.h:65` | low | Doc claims "H, C, N, O, F" but `getSupportedElements()` returns `{1,6,7,8,9,16,17}` (adds S, Cl). Doc drift. |
| Q-54 | `eht_method.cpp:367–381` | low | `initializeEHT()` defined, never called; `updateEHTParameters()` reads JSON not the PARAM registry (EHT has no PARAM definition). |
| Q-55 | all wrappers | low | Citation registration duplicated per wrapper as if/else chains. Per-method static citation table. |
| Q-56 | all wrappers | medium | `#ifdef USE_*` `#else` branches return zeros/false — a build without the dep silently returns zero energies rather than failing at MethodFactory. Mixed responsibility. |
| Q-57 | all interfaces | low | Every interface stores `mutable ConfigManager m_config;` mutated in const getters — smell; mutability works around const-correctness gaps in `ConfigManager::get`. |
| Q-58 | interfaces | medium | Verbosity plumbing inconsistent: XTB sets in `Calculation`, TBLite in ctor+`Calculation` (redundant), Ulysses passes a binary `>=3` bool (loses levels 1/2), EHT uses CurcumaLogger directly, ORCA passes the int. No single pattern. |

---

## 4. EnergyCalculator / method dispatch (`src/core/`)

### Confirmed / high-impact

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| D-1 | `energycalculator.cpp:42–56` | high | `reattachMethodScopes` hard-codes a `const char*` list of method scopes and re-merges them into `m_controller` after `ConfigManager("energycalculator", …)` stripped them, then **re-creates** the method. Every JSON-constructed `EnergyCalculator` builds the method **twice** when a method scope is present — wasted work + correctness trap if a method has ctor side effects. Fix ConfigManager to preserve scopes, or pass raw controller JSON straight to `MethodFactory::create`. | **verified** (mechanism: reattach list present) |
| D-2 | `energycalculator.cpp:153–158,237–244,547–624` | high | Silent failure cascade: method-creation failure returns soft `m_error`; constructor completes with `m_method==nullptr`; then `CalculateEnergy`/`Charges`/`Dipole`/`BondOrders`/`Energies`/`OrbitalOccuptations` all return **zeros** with no exception. Optimisers silently see `E=0`. Make the constructor throw (the factory already throws `MethodCreationException` — let it propagate). | — |
| D-3 | `energycalculator.cpp:421–430` | high | `CalculateEnergy` global-verbosity save/restore is **not RAII / not exception-safe**: the `try/catch` only catches `std::exception`; an early-return at `m_method->hasError()` (437→442) or a non-`std::exception` leaves `original_verbosity` unrestored. `EnergyCalculator` is **not** a `CurcumaMethod` subclass, so the CLAUDE.md RAII scoping does not apply here. Wrap in a guard. | — |
| D-4 | `energycalculator_enums.h:1–130` (whole file) | high | Dead duplicate `MethodRegistry` with substring-matching `detectMethodType` (`name.find(method_name)!=npos`, line 103) — `detectMethodType("gfn1")` matches `xtb-gfn1`, `GFN1L`, etc. No callers in `energycalculator.cpp`/`method_factory.cpp`. Divergent (lists gfn1/gfn2 under TBLite, contradicting AP3 native-canonical). Open TODOs at 78/126. **Delete.** | — |
| D-5 | `method_factory.cpp:460–517` vs `293–323` | high | Two divergent GFN-FF creation paths: `create()` routes `"gfnff"` to the GPU-aware native path (460–517), `"xtb-gfnff"` (566) to `createGFNFF` (External>XTB>Native, no GPU). Comment at 293–296 acknowledges the split. `xtb-gfnff -gpu cuda` silently ignores GPU. Unify. | — |
| D-6 | `curcumamethod.h:35–148` + `computational_method.h:43–360` | high | Two parallel abstract bases with overlapping verbosity/error/JSON responsibilities; **three** verbosity-management strategies coexist (CurcumaMethod RAII, EnergyCalculator manual save/restore, `setVerbosity` override flag). Root of CLAUDE.md Known Issue #3. Unify in `CurcumaLogger` with a scoped guard. | — |

### Dispatch / config

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| D-7 | `method_factory.h:155–160` vs `method_factory.cpp:66–95` | medium | Header method-list comments diverge from .cpp (e.g. `m_tblite_methods` comment says `{"ipea1","gfn1","gfn2"}`, actual `{"ipea1"}`; `m_xtb_methods` comment has `gfnff` typo vs `xtb-gfnff`). Five of six lists (`m_ff_methods`, `m_tblite_methods`, `m_xtb_methods`, `m_d3_methods`, `m_d4_methods`) are **never read** — only `m_ulysses_methods` is consumed. Delete unused + stale comments. |
| D-8 | `method_factory.cpp:418–599` | medium | ~180-line flat if/else dispatch; ADR comment (390) cites a "GCC 15 brace-init with std::function issue" as the reason — weak justification. New methods require editing this function. Registry/map. |
| D-9 | `method_factory.cpp:528–543` | medium | `xtb-gfn1/2` resolve TBLite>XTB; `xtb-gfnff` resolves External>XTB. Two different fallback orders for the same `xtb-*` family. Pick one + document. |
| D-10 | `method_factory.cpp:530–534,546–563` | medium | `if (hasTBLite())` runtime guard wrapping an `#ifdef USE_TBLITE` block — the runtime check is a constant either way, so one guard is redundant. Keep only `#ifdef`. |
| D-11 | `method_factory.cpp:493–516` | medium | GPU `#ifdef` matrix inconsistent: CUDA gates on `USE_CUDA` alone; ROCm on `USE_ROCM_GFNFF`; Vulkan is a hard-coded CPU fallback + warning regardless of `USE_VULKAN`. Three conventions in one function. |
| D-12 | `method_factory.cpp:180–228` vs `474–513` | low | `resolveNativeXtbGpuMode` duplicated inline for GFN-FF. Reuse. |
| D-13 | `method_factory.cpp:685–691` | low | `getMethodInfo` checks phantom `cgfnff` (never produced by `create()`, only in dead `energycalculator_enums.h:32`). Dead branch. |
| D-14 | `method_factory.cpp:101–158` | low | `checkCompilationFlag` string-dispatch wrapper re-lists the six `has*()` flags; unused outside the file. Delete; callers use `has*()` directly. |
| D-15 | `energycalculator.cpp:47-56` + `method_factory.cpp:418` | high | Config JSON round-tripped through `ConfigManager` then re-merged (lossy + patched by `reattachMethodScopes`). `MethodFactory::create` should accept `const ConfigManager&` directly (wrappers already build one internally). |
| D-16 | `energycalculator.cpp:143,206–209` | low | Stringly-typed config keys (`"multi"`,`"threads"`,`"gpu"`,…) scattered, no central key set. Typos = silent no-op. `constexpr` or Parameter Registry. |
| D-17 | `energycalculator.cpp:488–541` | medium | `NumGrad` uses `CurcumaLogger::get_verbosity()` directly (ignores `setVerbosity` override) and re-runs the global-verbosity mutation 6N times. Use `getEffectiveVerbosity()`. |
| D-18 | `energycalculator.cpp:571–588` | low | `BondOrders` packs the flat `getBondOrders()` vector into `result[0]` (one nested row). Old API mismatch. Document or split. |
| D-19 | `energycalculator.cpp:497` | low | Stale Chinese comment `// Step size in Angstrom (更适合 Bohr units)` — contradictory. |
| D-20 | `energycalculator.h:369,380` | medium | `Interface()` raw `ComputationalMethod*` (`@deprecated` but used) and `getCN()` duplicate of `CN()`. Delete, migrate callers. |
| D-21 | `energycalculator.h:446–447` | medium | Both `m_controller` and `m_parameter` stored; `setParameter` mirrors into the method but `m_controller` never updated — two sources of truth. |
| D-22 | `energycalculator.h:463–464` | low | `m_gpu_fallback`/`m_gpu_fallback_warned` uninitialised (others use `= false`). |
| D-23 | `energycalculator.h:184` | low | Typo `OrbitalOccuptations` (missing `a`). Rename + alias. |

### `computational_method.h` / `curcumamethod.*`

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| D-24 | `computational_method.h:43–360` | medium | God interface: ~35 virtuals mixing QM-only and FF-only accessors; pure-virtual `getEnergyDecomposition` forced on every method (even EHT). Split into `QMMethod`/`ForceFieldMethod` intermediates or a capability mixin. |
| D-25 | `computational_method.h:142–146` | low | `NumGrad` default returns `Matrix::Zero(0,3)` — "not implemented" indistinguishable from "no gradient". `EnergyCalculator::NumGrad` never calls it (reimplements FD). Wire or delete. |
| D-26 | `computational_method.h:121–123` | low | `copyGradientTo` workaround for CUDA heap corruption is now architecture. Fix the CUDA heap issue at source; remove the workaround. |
| D-27 | `computational_method.h:156` vs `energy_calculators/CLAUDE.md` | low | Doc shows `supportsGradients()` as the interface signature; real header has `hasGradient()`. Doc also claims `gfn2: TBLite→Ulysses→XTB` (contradicts AP3 native-canonical). Fix the doc. |
| D-28 | `computational_method.h:354–359` | low | `protected` `m_initialized/m_has_error/…` unused by the contract — every subclass defines its own. Enforce or drop. |
| D-29 | `curcumamethod.h:53,64` | low | `// TODO make pure virtual` left on `Initialise()`/`start()`. Decide or remove. |
| D-30 | `curcumamethod.h:107–110` | medium | `setVerbosity` writes global logger directly and updates `m_verbosity` but **not** `m_saved_global_verbosity` → a later sub-object `setVerbosity` is discarded by the dtor restore. Capture in `setVerbosity` too. |
| D-31 | `curcumamethod.cpp:43–93,96–139` | medium | Two near-identical constructors (only verbosity source differs). One delegating ctor. |
| D-32 | `curcumamethod.cpp:43–93` | medium | Legacy `silent` bool ctor + `m_silent`/`m_verbose` legacy bools kept "for backwards compatibility" but only `UpdateController` (cpp:250) reads `m_silent`. Drop. |
| D-33 | `curcumamethod.cpp:183–188` | medium | Two `catch (nlohmann::json::type_error&)` with **empty bodies** in `TriggerWriteRestart` — failed restart silently lost. `CurcumaLogger::warn`. |
| D-34 | `curcumamethod.cpp:224–237` | low | `LoadControl` unconditionally `std::cout << control`; `throw 404;` throws an integer. Use a proper exception type; remove the print. |
| D-35 | `curcumamethod.cpp:263–279` | low | `CheckStop` `#ifdef C17/_WIN32` nesting reads as a bug (dangling `std::ifstream` under `#ifdef C17`). Flatten. |
| D-36 | `curcumamethod.cpp:90,136` | medium | `CurcumaLogger::initCitationRegistry()` called in every ctor (per-instance, not once-per-process). Repeated work under `threads>1`. Call once in `main`. |

### Citation registry / threading

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| D-37 | `citation_registry.cpp:21–63,160` | medium | `clear()` wipes the global `m_seen`/`m_cited_keys`/`m_subrefs` under the lock but **does not** invalidate the `thread_local tls_seen` fast-path cache. CLAUDE.md "RESOLVED" claim is overstated — the *crash* is fixed, but the `clear()` semantic hole remains. Add a generation counter bumped in `clear()`. |
| D-38 | `citation_registry.cpp:46–48` | low | Unknown-key path re-warns on every call from every thread (not cached). A misconfigured `cite("typo")` per step spams. Cache negatives or rate-limit. |
| D-39 | `citation_registry.cpp:65–120` | low | `printTree` self-referential `std::function` capturing `[&]` — cyclically-referenced `std::function` (UB-adjacent). Plain recursive function. |
| D-40 | `citation_registry.cpp:16–19` | medium | Process-wide statics shared across all threads — two concurrent `EnergyCalculator` instances share one registry; a `clear()` from one loses the other's citations. Document or make per-instance. |

### Naming / docs / dead

| ID | Location | Sev | Finding |
|----|----------|-----|---------|
| D-41 | `method_factory.h:175–180` | low | `CURCUMA_HAS_TBLITE()` etc. macros unused. Delete. |
| D-42 | `method_factory.cpp:60` | low | `using namespace std;` in a header-included .cpp. Remove. |
| D-43 | `method_factory.cpp:56,58` | low | `#include <iostream>` unused after fmt migration. |
| D-44 | multiple | low | "Claude Generated: Big-Bang refactoring" tag in `energycalculator.h:18`, `computational_method.h:41`, `method_factory.h:41`, `method_factory.cpp:391`. Drop "Big-Bang" prefix, keep "Claude Generated". |
| D-45 | `method_factory.cpp:393–397` | low | ADR comment still says "Priority for GFN2: TBLite > Ulysses > XTB > Native" — AP3 made native canonical. Update the ADR. |

---

## 5. External third-party libraries — maintenance status

Operator directive (2026-06-26): the vendored/linked third-party computational
libraries are **legacy, kept only for tests and verification**, not for production
science. The native Curcuma implementations (GFN1/GFN2 xTB, GFN-FF) are the
canonical paths going forward.

| Library | Interface | Status | Keep? | Note |
|---------|-----------|--------|-------|------|
| **TBLite** | `tbliteinterface.cpp` (`USE_TBLITE`) | actively developed upstream, new methods | **yes** | worth maintaining — only external lib with a forward path; canonical cross-check for the native GFN1/GFN2/ALPB/GBSA validation |
| DFT-D3 (s-dftd3) | `dftd3interface.cpp` (`USE_D3`) | legacy | legacy/tests only | dispersion reference check; native D4 path exists |
| DFT-D4 (cpp-d4) | `dftd4interface.cpp` (`USE_D4`) | minimally supported | legacy/tests only | **not** used by native GFN2 — native GFN2 has its own D4 (`curcuma::dispersion::D4Evaluator` + `D4ChargeModel`, `src/core/energy_calculators/dispersion/`, built unconditionally, no cpp-d4 link); the cpp-d4 source is only transcribed into native data tables / comments |
| XTB binary | `xtbinterface.cpp` (`USE_XTB`) | legacy / frozen | legacy/tests only | upstream xtb is Fortran-frozen; superseded by native xTB |
| Ulysses | `ulyssesinterface.cpp` | legacy | legacy/tests only | NDDO reference (PM3/PM6/AM1/MNDO); native PM3/AM1/MNDO exist |
| GFN-FF (XTB) | `gfnffinterface.cpp` (`USE_GFNFF`) | legacy | legacy/tests only | native `gfnff` (`ff_methods/gfnff_method.cpp`) is canonical |

Implication for the debt items above: correctness bugs in the legacy interfaces
(Q-1 XTB setMethod, Q-3 tblite exit, Q-4..Q-7 DFT-D3/D4, Q-10 ctor logging,
Q-11 Ulysses HOMO) are worth fixing only insofar as they keep the **verification
harness** working — they are not user-facing production paths. TBLite is the
exception: treat its interface as maintained.

---

## Master priority list (cross-cutting)

Correctness / crash / UB first, then maintainability. "v" = spot-verified by re-reading cited lines.

| # | ID | Area | Finding | Verified |
|---|-----|------|---------|----------|
| 1 | Q-1 | XTB iface | `xtb-gfn2/1` via XTB-binary backend silently runs **GFN0** (setMethod never called) | v |
| 2 | Q-2 | ORCA | Open-shell inputs get wrong multiplicity (`2*m_spin+1` should be `m_spin+1`) | v |
| 3 | Q-4 | DFT-D3 | `InitialiseMolecule(vector<int>)` fills coords with `index/au` — garbage geometry | v |
| 4 | Q-5 | DFT-D3 | `delete` on `new[]` arrays in destructor — **undefined behavior** | v |
| 5 | Q-7 | DFT-D3 | `m_damping.compare("bj")` inverted — custom-param path loads wrong damping model | v |
| 6 | Q-6 | DFT-D4 | z-coordinate not scaled by `factor` while x/y are — sheared geometry for `factor!=1` | v |
| 7 | Q-3 | TBLite | `exit(1)` inside a library when GB solvation requested — kills host process | v |
| 8 | X-L4 | native xTB | DC/sparse return **zero gradient** with no error flag — `-opt`/`-md` silently stuck | — |
| 9 | X-L3 | native xTB | Sparse mode is O(N³) dense (`m_S.inverse()`) though marketed as O(N) | v |
| 10 | F-Q4 | GFN-FF | EEQ silent uniform-charge fallback, no `m_has_error` — wrong charges propagate | — |
| ~~11~~ | ~~F-Q6~~ | GFN-FF | ~~`gfnff_torsions_NEW.cpp` dead duplicate file~~ — **deleted 2026-06-26** | — |
| 12 | Q-8 | TBLite | Re-init with different natoms overflows `m_coord`/`m_attyp` (no realloc) | — |
| 13 | Q-9 | XTB | `xtb_newMolecule` error check commented out — failure silently `m_initialised=true` | — |
| 14 | X-I1 | native xTB | d-shell AOs silently zeroed in S/H0 — undocumented element-range limit | — |
| 15 | D-2 | dispatch | Silent failure cascade: bad method → `E=0`/zero properties, no exception | — |
| 16 | F-Q7 | GFN-FF | External `GFNFFInterface::Charges()` always returns zeros | — |
| 17 | D-1 | dispatch | `reattachMethodScopes` builds every method **twice** (ConfigManager scope-stripping workaround) | v |
| 18 | D-3 | dispatch | `CalculateEnergy` verbosity save/restore not RAII/exception-safe | — |
| 19 | D-4 | dispatch | `energycalculator_enums.h` whole file dead + divergent substring-matching registry | — |
| 20 | D-5 | dispatch | Two divergent GFN-FF creation paths; `xtb-gfnff -gpu` silently ignored | — |
| 21 | X-M1 | native xTB | ~15 methods declared in `xtb_method.h` with no definition → link errors | — |
| 22 | F-E1 | GFN-FF | `gfnff_method.cpp` 11 738-line god file, `Calculation()` 530 lines | — |
| 23 | X-G2 | native xTB | `calculateGradientGpu` duplicates ~160 lines of gradient sections | — |
| 24 | X-L1/L2 | native xTB | DC+sparse duplicate global setup + `packSCC`/`unpackSCC` (3 copies) | — |
| 25 | X-R1 | native xTB | `computeMullikenChargeResponse` ~390-line god function | — |
| 26 | D-6 | dispatch | Two parallel abstract bases + three verbosity strategies (Known Issue #3 root) | — |
| 27 | Q-10 | wrappers | TBLite/Ulysses ctors log unconditionally — break Level 0 silent mode | — |
| 28 | X-S2 | native xTB | Dead divergence check with 100× different threshold — wiring hazard | — |
| 29 | Q-11 | Ulysses | HOMO index = atomcount/2 — reported orbital analysis is garbage | — |
| 30 | D-37 | citation | `clear()` doesn't invalidate `thread_local` cache — "RESOLVED" overstated | — |

---

## Notes on the audit method

- Findings generated by three parallel code-reading passes (GFN-FF + native xTB; QM interfaces; dispatch), then consolidated.
- Every finding cites a `file:line`; the 11 marked **verified** were re-read directly before being recorded as fact, per the project's "show proof, don't assert" rule.
- Native xTB and GFN-FF GPU-kernel internals were deliberately **not** audited (out of scope); only CPU-path / `#ifdef`-leak / non-GPU debt there.
- This is documentation only — no fixes applied. Each row's "suggested fix" is a one-line pointer, not a reviewed patch.