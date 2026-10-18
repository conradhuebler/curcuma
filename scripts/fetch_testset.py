#!/usr/bin/env python3
# Claude Generated (Sep 2026)
"""Fetch Grimme-group benchmark test sets on demand.

The raw structures/reference data for MOR41, GMTKN55 and S30L are large and
have their own citation terms, so curcuma does not commit them to git. This
script downloads (or clones) each set into the exact directory layout the
existing validation harnesses already expect
(scripts/mor41_validation.py, scripts/s30l_*.py), so a validation run only
needs one extra command first:

    python scripts/fetch_testset.py mor41
    python scripts/mor41_validation.py

Registry entries are plain dicts (name -> {kind, url, dest, citation, note,
check}). "kind" is one of:
  tar    - download + extract a tar(.gz) archive into `dest`
  git    - shallow-clone a git repository into `dest`
  manual - no automatable source; print instructions and exit

Fetched data is gitignored (see .gitignore, "Grimme-group benchmark test-set
structures"); only the small hand-curated driver files (reactions.dat,
reference_s30l, ...) are tracked.

Usage:
    python scripts/fetch_testset.py list
    python scripts/fetch_testset.py fetch mor41
    python scripts/fetch_testset.py fetch mor41 gmtkn55 --force
"""
import argparse
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TEST_CASES = REPO / "test_cases"

REGISTRY = {
    "mor41": {
        "kind": "tar",
        "url": "https://www.chemie.uni-bonn.de/grimme/de/software/mor41/geometries-tar.gz",
        "dest": TEST_CASES / "MOR41-testset",
        "citation": "S. Dohm, A. Hansen, M. Steinmetz, S. Grimme, M. P. Checinski, "
                    "J. Chem. Theory Comput. 2018, 14, 2596. DOI: 10.1021/acs.jctc.7b01183",
        "note": "41 organometallic reactions, 95 structures. Extracts directly as "
                "<dest>/<name>/mol.xyz (reactions.dat is hand-curated and already "
                "tracked in git; the archive does not contain it, so it is untouched).",
        "check": lambda dest: len(list(dest.glob("*/mol.xyz"))) >= 90,
    },
    "gmtkn55": {
        "kind": "git",
        "url": "https://github.com/grimme-lab/GMTKN55.git",
        "dest": TEST_CASES / "GMTKN55-testset",
        "citation": "L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, "
                    "Phys. Chem. Chem. Phys. 2017, 19, 32184. DOI: 10.1039/C7CP04913G",
        "note": "55 main-group thermochemistry/kinetics/non-covalent-interaction "
                "subsets, cloned as published (geometries + reference data + the "
                "upstream eval.py). No curcuma-side validation script consumes this "
                "yet; use it directly with eval.py once curcuma output is converted.",
        "check": lambda dest: (dest / "eval.py").exists(),
    },
    "s30l": {
        "kind": "manual",
        "url": "https://pubs.acs.org/doi/10.1021/acs.jctc.5b00296",
        "dest": TEST_CASES / "s30l_test_set",
        "citation": "M. Sure, S. Grimme, J. Chem. Theory Comput. 2015, 11, 3785. "
                    "DOI: 10.1021/acs.jctc.5b00296",
        "note": "30 host-guest complexes. Not hosted on the Grimme group software "
                "page (https://www.chemie.uni-bonn.de/grimme/de/software) - only "
                "distributed as the ACS Supporting Information, which refuses "
                "automated fetches (HTTP 403). Download the SI by hand from the URL "
                "above and place it as:\n"
                "    test_cases/s30l_test_set/<1..30>/{A,B,AB}/coord   "
                "(Turbomole $coord block, Bohr)\n"
                "    test_cases/s30l_test_set/<1..30>/{A,B,AB}/.CHRG   "
                "(charged systems only, e.g. '+1')\n"
                "    test_cases/s30l_test_set/reference_s30l           "
                "(30 whitespace-separated association energies, kcal/mol)",
        "check": lambda dest: (dest / "reference_s30l").exists(),
    },
}

USER_AGENT = "curcuma-fetch-testset/1 (+https://github.com/, research use)"


def log(msg):
    print(msg, flush=True)


def write_provenance(dest, name, entry):
    dest.mkdir(parents=True, exist_ok=True)
    text = (
        f"Test set: {name}\n"
        f"Source:   {entry['url']}\n"
        f"Fetched:  {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
        f"by scripts/fetch_testset.py\n"
        f"Citation: {entry['citation']}\n"
    )
    (dest / "PROVENANCE.txt").write_text(text)


def fetch_tar(name, entry, force):
    dest = entry["dest"]
    if not force and dest.exists() and entry["check"](dest):
        log(f"[{name}] already present at {dest} (use --force to re-fetch)")
        return True
    log(f"[{name}] downloading {entry['url']}")
    req = urllib.request.Request(entry["url"], headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=120) as resp, \
                tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=False) as tmp:
            shutil.copyfileobj(resp, tmp)
            tmp_path = Path(tmp.name)
    except (urllib.error.URLError, urllib.error.HTTPError, OSError) as exc:
        log(f"[{name}] ERROR: download failed: {exc}")
        return False
    try:
        if not tarfile.is_tarfile(tmp_path):
            log(f"[{name}] ERROR: downloaded file is not a tar archive "
                f"(site layout may have changed - check {entry['url']} by hand)")
            return False
        dest.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tmp_path) as tf:
            tf.extractall(dest)
    finally:
        tmp_path.unlink(missing_ok=True)
    if not entry["check"](dest):
        log(f"[{name}] ERROR: extracted but the expected layout was not found "
            f"under {dest} (site layout may have changed)")
        return False
    write_provenance(dest, name, entry)
    log(f"[{name}] OK: {dest}")
    return True


def fetch_git(name, entry, force):
    dest = entry["dest"]
    if dest.exists():
        if not force and entry["check"](dest):
            log(f"[{name}] already present at {dest} (use --force to re-fetch)")
            return True
        if force:
            shutil.rmtree(dest)
        else:
            log(f"[{name}] ERROR: {dest} exists but does not look like a valid "
                f"checkout; remove it or pass --force")
            return False
    log(f"[{name}] cloning {entry['url']}")
    try:
        subprocess.run(["git", "clone", "--depth", "1", entry["url"], str(dest)],
                        check=True, capture_output=True, text=True, timeout=600)
    except subprocess.CalledProcessError as exc:
        log(f"[{name}] ERROR: git clone failed: {exc.stderr.strip()}")
        return False
    except subprocess.TimeoutExpired:
        log(f"[{name}] ERROR: git clone timed out")
        return False
    if not entry["check"](dest):
        log(f"[{name}] ERROR: cloned but the expected layout was not found "
            f"under {dest}")
        return False
    write_provenance(dest, name, entry)
    log(f"[{name}] OK: {dest}")
    return True


def fetch_manual(name, entry, force):
    dest = entry["dest"]
    if not force and dest.exists() and entry["check"](dest):
        log(f"[{name}] already present at {dest}")
        return True
    log(f"[{name}] no automated source available:\n{entry['note']}\n"
        f"Citation: {entry['citation']}")
    return False


FETCHERS = {"tar": fetch_tar, "git": fetch_git, "manual": fetch_manual}


def do_fetch(name, force):
    entry = REGISTRY[name]
    return FETCHERS[entry["kind"]](name, entry, force)


def do_list():
    for name, entry in REGISTRY.items():
        present = entry["dest"].exists() and entry["check"](entry["dest"])
        status = "present" if present else "missing"
        log(f"{name:10s} [{entry['kind']:6s}] [{status:7s}] {entry['dest'].relative_to(REPO)}")
        log(f"           {entry['citation']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list", help="show registered test sets and their local status")
    fp = sub.add_parser("fetch", help="fetch one or more test sets")
    fp.add_argument("names", nargs="+", choices=list(REGISTRY) + ["all"])
    fp.add_argument("--force", action="store_true",
                     help="re-fetch even if already present")
    args = ap.parse_args()

    if args.cmd == "list":
        do_list()
        return

    names = list(REGISTRY) if "all" in args.names else args.names
    ok = True
    for name in names:
        ok = do_fetch(name, args.force) and ok
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
