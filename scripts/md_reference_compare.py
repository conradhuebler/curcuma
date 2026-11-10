#!/usr/bin/env python3
"""Compare MD trajectories of a reference curcuma binary and a new one (Claude Generated, Oct 2026).

Runs `curcuma -md INPUT OPTIONS` with the reference binary twice (the control: the reference must agree with
itself), with the new binary without extra options, and with the new binary plus EXTRA (an option that is supposed
to leave the dynamics alone). The trajectories are compared byte for byte after the version string at the end of
each comment line is removed, since the two binaries carry different version strings.

    scripts/md_reference_compare.py --ref REF/curcuma --new NEW/curcuma system.xyz \
        --options="-max_time 100 -thermostat none -seed 1 -threads 1 -no_restart true -dump_frequency 5" \
        --extra="-hb_log_frequency 10"

Use one thread and a fixed seed: with several threads the reduction order may differ between runs. The reference
binary is built from the commit before the change (git archive of that commit, same CMake options as the new one).
Exit status 0 when the control holds and both new runs equal the reference, 1 otherwise.
"""
import argparse
import glob
import hashlib
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile


def trajectory_digest(binary, label, source, options, work):
    directory = os.path.join(work, label)
    os.makedirs(directory)
    shutil.copy(source, directory)
    run = subprocess.run([binary, "-md", os.path.basename(source)] + options, cwd=directory,
                         capture_output=True, text=True)
    stem = os.path.splitext(os.path.basename(source))[0]
    found = glob.glob(os.path.join(directory, "*", stem + ".trj.xyz")) + glob.glob(os.path.join(directory, stem + ".trj.xyz"))
    if run.returncode != 0 or not found:
        return None, 0
    text = "".join(re.sub(r"Curcuma .*$", "", line) for line in open(found[0]))
    return hashlib.md5(text.encode()).hexdigest(), len(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input")
    parser.add_argument("--ref", required=True)
    parser.add_argument("--new", required=True)
    parser.add_argument("--options", default="", help="options of the run, as one string")
    parser.add_argument("--extra", default="", help="options only the second new run gets")
    arguments = parser.parse_args()

    options = shlex.split(arguments.options)
    extra = shlex.split(arguments.extra)
    source = os.path.abspath(arguments.input)
    with tempfile.TemporaryDirectory() as work:
        runs = [("ref1", arguments.ref, []), ("ref2", arguments.ref, []),
                ("new", arguments.new, []), ("new_extra", arguments.new, extra)]
        digest = {}
        for label, binary, more in runs:
            digest[label] = trajectory_digest(os.path.abspath(binary), label, source, options + more, work)
            print(f"{label:10s} md5 {digest[label][0]}  ({digest[label][1]} characters compared)")
    control = digest["ref1"][0] is not None and digest["ref1"] == digest["ref2"]
    same_new = digest["new"] == digest["ref1"]
    same_extra = digest["new_extra"] == digest["ref1"]
    print(f"control (reference agrees with itself): {control}")
    print(f"new binary, no extra options equals reference: {same_new}")
    print(f"new binary, with the extra options equals reference: {same_extra}")
    return 0 if control and same_new and same_extra else 1


if __name__ == "__main__":
    sys.exit(main())
