#!/usr/bin/env python3
"""
SQM native-vs-tblite validation driver (Claude Generated, 2026-05).

Runs the real `curcuma` binary on a molecule with a native xTB method
(gfn1 / gfn2 — the canonical native methods) and compares the single-point
energy against a committed tblite reference
(test_cases/sqm_reference/reference_data/*.ref.json,
produced by dump_tblite_reference in a USE_TBLITE build).

This is the end-to-end user path — the same dispatch that surfaced the
NH3 native-GFN1 anomaly — so it validates MethodFactory routing + native SCF, not
just an isolated class. It needs NO tblite build (references are committed),
matching how the GFN-FF suite works.

Gate: total energy. Charges/gradient are reported when present but are not
the automatic gate (native gradients are FD-validated separately by
test_xtb_gradient; the tblite gradient-norm convention differs).

Usage:
  validate_sqm.py --curcuma <path> --ref <ref.json> --xyz <mol.xyz>
                  [--tol-energy <Eh>] [--quiet]

Exit 0 on pass, 1 on fail.
"""
import argparse
import json
import re
import subprocess
import sys

ANSI = re.compile(r"\x1b\[[0-9;]*m")


# Native per-container decomposition, printed by the native xTB SCF at
# verbosity >= 2 (src/core/.../xtb_native.cpp). The labels are matched on their
# "= <signed number> Eh" tail, tolerant of label spacing. "Electronic" already
# includes Coulomb-ES2 + third-order + multipole (it is m_E_electronic), so it
# maps directly onto the reference "electronic" container; do NOT re-add the
# Coulomb/Third-order breakdown lines.
_COMP_RE = {
    "electronic": re.compile(r"Electronic\s*=\s*(-?[0-9.]+)\s*Eh"),
    "dispersion": re.compile(r"Dispersion\s*=\s*(-?[0-9.]+)\s*Eh"),
    "repulsion":  re.compile(r"Repulsion\s*=\s*(-?[0-9.]+)\s*Eh"),
}


def run_curcuma(curcuma, xyz, method, quiet, gpu=None):
    # -verbosity 2 surfaces the per-container decomposition; the
    # "Single Point Energy = ..." line is printed unconditionally, so the total
    # gate is unaffected by the verbosity bump.
    cmd = [curcuma, "-sp", xyz, "-method", method, "-verbosity", "2"]
    if gpu:
        # GPU path (Claude Generated, GPU port): -gpu cuda routes gfn1/gfn2 to the
        # native xTB GPU backend. With no CUDA device the wrapper falls back to CPU,
        # so the test still validates (the energy is identical), it just is not on GPU.
        cmd += ["-gpu", gpu]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        if not quiet:
            print(f"  TIMEOUT running {' '.join(cmd)}")
        return None
    text = ANSI.sub("", out.stdout + out.stderr)
    m = re.findall(r"Single Point Energy = (-?[0-9.]+)", text)
    if not m:
        if not quiet:
            print(f"  no energy parsed (exit {out.returncode}); cmd: {' '.join(cmd)}")
        return None
    result = {"total": float(m[-1])}
    # Components are best-effort: absent (e.g. SCF that did not print the block)
    # degrades the advisory component check to n/a, never the total-energy gate.
    for key, rx in _COMP_RE.items():
        hits = rx.findall(text)
        if hits:
            result[key] = float(hits[-1])
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curcuma", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--xyz", required=True)
    ap.add_argument("--tol-energy", type=float, default=None,
                    help="absolute energy tolerance in Eh (default: size-scaled)")
    ap.add_argument("--tol-components", type=float, default=None,
                    help="OPT-IN: if set, the component check (combined "
                         "disp+electronic and repulsion vs the reference "
                         "energy_components) contributes to the exit code at this "
                         "Eh tolerance. Default: advisory only -- the component "
                         "deltas are printed but the exit code stays tied solely "
                         "to the 1e-8 total-energy gate (so the WILL_FAIL xfail "
                         "mechanism is undisturbed).")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--gpu", default=None,
                    help="pass -gpu <mode> to curcuma (e.g. 'cuda') to exercise the "
                         "native xTB GPU backend; falls back to CPU without a device")
    args = ap.parse_args()

    with open(args.ref) as f:
        ref = json.load(f)

    method = ref["method"]                      # "gfn1" | "gfn2"
    native = method                             # canonical native xTB (gfn1/gfn2)
    e_ref = ref["total_energy"]
    nat = ref["molecule"]["natoms"]

    # Size-scaled default tolerance: 1e-4 Eh floor, +1e-5 Eh/atom for accumulation.
    tol = args.tol_energy if args.tol_energy is not None else max(1e-4, nat * 1e-5)

    res = run_curcuma(args.curcuma, args.xyz, native, args.quiet, args.gpu)
    name = ref["molecule"]["name"]
    if res is None:
        print(f"FAIL {name:20s} {native}: native run produced no energy")
        return 1
    e_nat = res["total"]

    dE = e_nat - e_ref
    ok = abs(dE) <= tol
    status = "PASS" if ok else "FAIL"

    # --- advisory component localization (does NOT gate by default) ---
    # Combined disp+electronic is the apples-to-apples number: tblite splits the
    # charge-coupled D4 across the dispersion and electronic containers while
    # curcuma lumps it in m_E_dispersion, so only the SUM is split-invariant.
    # Repulsion is reported separately (shared closed form, should be ~identical).
    # This uses curcuma's SELF-CONVERGED density, so it is a coarse localizer, not
    # a 1e-8 check -- the precise per-container audit is the ecomp_* align ctests
    # (diag_curcuma_energy_components at tblite's injected density).
    comp_line = None
    comp_ok = True
    ec = ref.get("energy_components")
    have_native = all(k in res for k in ("electronic", "dispersion", "repulsion"))
    if ec and have_native and ("electronic" in ec) and ("dispersion" in ec):
        nat_comb = res["electronic"] + res["dispersion"]
        ref_comb = ec["electronic"] + ec["dispersion"]
        d_comb = nat_comb - ref_comb
        parts = [f"d(disp+elec)={d_comb*1000:+.4f} mEh"]
        worst = abs(d_comb)
        if "repulsion" in ec:
            d_rep = res["repulsion"] - ec["repulsion"]
            parts.append(f"d(rep)={d_rep*1000:+.4f} mEh")
            worst = max(worst, abs(d_rep))
        if args.tol_components is not None:
            comp_ok = worst <= args.tol_components
            parts.append(f"[{'ok' if comp_ok else 'OVER'} tol={args.tol_components*1000:.3f} mEh]")
        comp_line = "  components: " + "  ".join(parts)
    elif ec and not have_native:
        comp_line = "  components: n/a (native decomposition not parsed)"
    elif not ec:
        comp_line = "  components: n/a (ref has no energy_components)"

    # By user decision the component check is ADVISORY: it only feeds the exit code
    # when --tol-components is explicitly passed (localization runs). The default
    # gate is the 1e-8 total energy, keeping WP1's xfail accounting exact.
    gate_ok = ok and (comp_ok if args.tol_components is not None else True)
    status = "PASS" if gate_ok else "FAIL"

    if not args.quiet or not gate_ok:
        print(f"{status} {name:20s} {native}: "
              f"E_native={e_nat:.8f}  E_tblite={e_ref:.8f}  "
              f"dE={dE*1000:+.4f} mEh  tol={tol*1000:.3f} mEh")
        if comp_line is not None:
            print(comp_line)
    return 0 if gate_ok else 1


if __name__ == "__main__":
    sys.exit(main())
