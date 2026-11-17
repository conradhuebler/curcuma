#!/usr/bin/env python3
"""Claude Generated (Nov 2026): EEQ solver `sytrf` (-eeq_solver.solve_method sytrf), the linear solve of the Fortran
reference (LAPACK dsytrf + dsytrs, uplo 'U', on the augmented (N+nfrag) matrix). AI-generated, machine-tested.

Mode cpu (ctest eeq_sytrf_cpu): sytrf against cholesky on the CPU, fresh directory per run, one thread.
  1. PATH      - both phases name the solver at verbosity 1: "EEQ solver (topology charges): sytrf [CPU]" and
                 "EEQ solver (final charges): sytrf [CPU]"; -gfnff.solve_method sytrf is the same request.
  2. SP        - caffeine (nfrag 1), H3O+(H2O)2 with -charge 1 (nfrag 3) under frag_charge_model ensemble (default)
                 and reference, Na2Cl2 (four isolated ions): |E(sytrf) - E(cholesky)| <= 1e-10 Eh (-dump_gradient,
                 15 decimals) and Phase-1 / Phase-2 charges (-eeq_solver.dump_charges) within 1e-10 e.
  3. MD        - caffeine, short MD (5 fs): per-step charges (-md_diagnostics) of sytrf within 1e-10 e of cholesky.
Mode rocm (ctest eeq_sytrf_rocm, USE_ROCM builds; exit 77 = skip without the plugin or a HIP device): short MD with
  -gpu rocm -gpu_strict true and sytrf against the CPU sytrf MD: the device path is named in the verbosity-1 line, no
  counted GPU fallback (exit code 0 under -gpu_strict), per-step charges within 1e-9 e (caffeine, trisaccharide-c18h32o16,
  macrocycle-bgal).
Usage: check_eeq_sytrf.py CURCUMA_BINARY SOURCE_DIR [cpu|rocm]
"""
import json, os, re, shutil, subprocess, sys, tempfile

curcuma, src = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
mode = sys.argv[3] if len(sys.argv) > 3 else "cpu"
LIB = os.path.join(src, "test_cases/structures")
XYZ = {
    "caffeine": os.path.join(LIB, "medium/caffeine.xyz"),
    "trisaccharide-c18h32o16": os.path.join(LIB, "large/trisaccharide-c18h32o16.xyz"),
    "macrocycle-bgal": os.path.join(LIB, "clusters/macrocycle-bgal.xyz"),
    "hydronium-water2": os.path.join(LIB, "clusters/hydronium-water2.xyz"),
    "na2cl2-isolated-ions": os.path.join(LIB, "clusters/na2cl2-isolated-ions.xyz"),
}
S = "-eeq_solver.solve_method"
ANSI = re.compile(r"\x1b\[[0-9;]*m")
failed = 0


def check(ok, what):
    global failed
    print(("  PASS  " if ok else "  FAIL  ") + what)
    failed += 0 if ok else 1


def run(xyz, extra, md=False):
    """Fresh directory (no topology cache). Returns rc, log, energy (sp), charge dump (sp), per-step charges (md)."""
    d = tempfile.mkdtemp(prefix="curcuma_sytrf_")
    try:
        shutil.copy(xyz, os.path.join(d, "m.xyz"))
        cmd = [curcuma, "-md" if md else "-sp", "m.xyz", "-method", "gfnff", "-no_bmt", "-threads", "1", "-verbosity", "1",
               "-gfnff.cache_topology", "false"]
        if md:
            cmd += ["-max_time", "5", "-md_diagnostics", "true", "-dump_frequency", "1", "-md.seed", "7"]
        else:
            cmd += ["-dump_gradient", "grad.txt", "-eeq_solver.dump_charges", "true"]
        p = subprocess.run(cmd + extra, cwd=d, capture_output=True, text=True, timeout=900)
        log = ANSI.sub("", p.stdout + p.stderr)
        e, dump, q = None, None, []
        g = os.path.join(d, "grad.txt")
        if os.path.exists(g):
            line = next((l for l in open(g) if l.startswith("# energy")), None)
            e = float(line.split()[2]) if line else None
        dumps = [f for f in os.listdir(d) if f.startswith("charges_dump_N") and f.endswith(".json")]
        if dumps:
            dump = json.load(open(os.path.join(d, dumps[0])))
        fn = os.path.join(d, "m.diag.jsonl")
        if os.path.exists(fn):
            q = [json.loads(l)["charges"] for l in open(fn) if l.strip()]
        return p.returncode, log, e, dump, q
    finally:
        shutil.rmtree(d, ignore_errors=True)


def max_dq(a, b):
    """Largest |a_i - b_i| over two equally long lists (1.0 if they do not match in length)."""
    if len(a) != len(b) or not a:
        return 1.0
    return max(abs(x - y) for x, y in zip(a, b))


def dump_charges(dump, key):
    return [at[key] for at in dump["atoms"]] if dump else []


def md_dq(qa, qb):
    n = min(len(qa), len(qb))
    if n < 2 or len(qa) != len(qb):
        return n, 1.0
    return n, max(max_dq(qa[i], qb[i]) for i in range(n))


if mode == "cpu":
    # 1. path lines, alias
    rc, log, _, _, _ = run(XYZ["caffeine"], [S, "sytrf"])
    check(rc == 0 and "EEQ solver (topology charges): sytrf [CPU]" in log and "EEQ solver (final charges): sytrf [CPU]" in log,
          "both phases run sytrf on the CPU")
    rc, log, _, _, _ = run(XYZ["caffeine"], ["-gfnff.solve_method", "sytrf"])
    check(rc == 0 and "EEQ solver (final charges): sytrf [CPU]" in log, "-gfnff.solve_method sytrf is the same request")

    # 2. single points against cholesky
    cases = [("caffeine", []),
             ("hydronium-water2", ["-charge", "1"]),
             ("hydronium-water2", ["-charge", "1", "-gfnff.frag_charge_model", "reference"]),
             ("na2cl2-isolated-ions", [])]
    for name, extra in cases:
        tag = name + (" " + " ".join(extra) if extra else "")
        rc_s, log_s, e_s, d_s, _ = run(XYZ[name], extra + [S, "sytrf"])
        rc_c, log_c, e_c, d_c, _ = run(XYZ[name], extra + [S, "cholesky"])
        ok = rc_s == 0 and rc_c == 0 and e_s is not None and e_c is not None and d_s and d_c
        check(ok, f"{tag}: both runs finished with energy and charge dump")
        if not ok:
            continue
        nfrag = d_s.get("nfrag")
        check("falling back" not in log_s, f"{tag}: no sytrf fallback (nfrag {nfrag})")
        check(abs(e_s - e_c) <= 1e-10, f"{tag}: |dE| = {abs(e_s - e_c):.2e} Eh <= 1e-10")
        dq1 = max_dq(dump_charges(d_s, "q_phase1"), dump_charges(d_c, "q_phase1"))
        dq2 = max_dq(dump_charges(d_s, "q_phase2"), dump_charges(d_c, "q_phase2"))
        check(dq1 <= 1e-10 and dq2 <= 1e-10, f"{tag}: max |dq| phase 1 {dq1:.2e} e, phase 2 {dq2:.2e} e <= 1e-10")

    # 3. MD, per-step charges
    rc_s, _, _, _, q_s = run(XYZ["caffeine"], [S, "sytrf"], md=True)
    rc_c, _, _, _, q_c = run(XYZ["caffeine"], [S, "cholesky"], md=True)
    n, dq = md_dq(q_s, q_c)
    check(rc_s == 0 and rc_c == 0 and n >= 2 and dq <= 1e-10, f"caffeine MD: {n} steps, max |dq| = {dq:.2e} e <= 1e-10")

elif mode == "rocm":
    plugin = os.path.join(os.path.dirname(os.path.abspath(curcuma)), "libcurcuma_rocm.so")
    if not os.path.exists(plugin):
        print("SKIP  no libcurcuma_rocm.so next to the binary")
        sys.exit(77)
    rc, log, _, _, _ = run(XYZ["caffeine"], ["-gpu", "rocm", S, "sytrf"])
    if re.search(r"no HIP device|hipErrorNoDevice|no usable HIP|plugin did not provide", log):
        print("SKIP  no usable ROCm device")
        sys.exit(77)
    for name in ("caffeine", "trisaccharide-c18h32o16", "macrocycle-bgal"):
        rc_g, log_g, _, _, q_g = run(XYZ[name], ["-gpu", "rocm", "-gpu_strict", "true", S, "sytrf"], md=True)
        rc_c, _, _, _, q_c = run(XYZ[name], [S, "sytrf"], md=True)
        check(rc_g == 0 and "sytrf [GPU ROCm, augmented sytrf (upper), host dsytrs" in log_g,
              f"{name}: ROCm MD runs the device sytrf path without a GPU fallback (rc {rc_g})")
        n, dq = md_dq(q_g, q_c)
        check(rc_c == 0 and n >= 2 and dq <= 1e-9, f"{name}: ROCm sytrf = CPU sytrf, {n} steps, max |dq| = {dq:.2e} e <= 1e-9")
else:
    sys.exit(f"unknown mode '{mode}' (cpu or rocm)")

print(f"{failed} check(s) failed" if failed else "all checks passed")
sys.exit(1 if failed else 0)
