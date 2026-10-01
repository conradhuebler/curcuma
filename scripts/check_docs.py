#!/usr/bin/env python3
"""Documentation and repository hygiene check for curcuma.

Claude Generated (Oct 2026). Checks the rules of the root CLAUDE.md section
"Where Things Go" so that documentation does not fragment again:

  ERRORS (exit 1)
    - a CLAUDE.md over its line budget (root and sub-directory budgets below)
    - a tracked markdown file in the repository root that is not on the whitelist
    - a relative markdown link that does not resolve to a tracked file
    - a tracked backup file (*.backup, *.orig, *.bak, *.rej)
    - a line over MAX_LINE characters in a CLAUDE.md
  WARNINGS (exit 0, exit 1 with --strict)
    - a path in backticks in a CLAUDE.md / README / TODO that matches no tracked file
    - a dated status heading in a CLAUDE.md, e.g. "Foo (Sep 2026)" or "Completed ..."
    - an untracked markdown file in the repository root, an untracked backup file

Scope: tracked files only (git ls-files), so build trees and local scratch files are ignored.
Pure standard library. Usage:  python3 scripts/check_docs.py [--strict] [--all-docs]
  --all-docs  also check relative links in every docs/*.md (many old ones are known to be stale)
Enable as pre-commit hook:  git config core.hooksPath scripts/git-hooks
"""
import os
import re
import subprocess
import sys

ROOT = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip()
os.chdir(ROOT)

ROOT_MD_WHITELIST = {"README.md", "CLAUDE.md", "TODO.md", "AIChangelog.md"}
BUDGET_ROOT = 500          # lines, root CLAUDE.md
BUDGET_SUB = 120           # lines, every other CLAUDE.md
MAX_LINE = 600             # characters, CLAUDE.md lines
BACKUP_RE = re.compile(r"\.(backup|orig|bak|rej)$")
DATED_HEADING_RE = re.compile(
    r"^#{1,6} .*\((?:[A-Z][a-z]{2,8}\.? )?(?:\d{1,2}, )?20\d\d(?:-\d\d(?:-\d\d)?)?[^)]*\)\s*$"
    r"|^#{1,6} .*\b(Completed|COMPLETE|Recent Major Achievements|Status Log)\b")
# paths that are legitimately not tracked: build products, caches, external checkouts, run outputs
UNTRACKED_OK = re.compile(
    r"^(external/|release|build|_deps/|generated/|_run/|\.cache/|(Wissen|Projekte|Labor)/|test_cases/.*/(stdout|stderr)\.log|"
    r"<|\.|~|/|curcuma_|stop$)|\.topo\.json|\.param\.json|\.hessian|\.opt\.|\.trj\.|BMT|^metadata\.json$")
PATH_RE = re.compile(r"`([A-Za-z0-9_./+-]+\.(?:cpp|h|hpp|cu|cuh|hip|py|sh|md|json|txt|cmake|in|yml))`")

errors, warnings = [], []


def git(*args):
    return subprocess.check_output(["git", *args], text=True).splitlines()


tracked = git("ls-files")
tracked_set = set(tracked)
basenames = {}
for p in tracked:
    basenames.setdefault(os.path.basename(p), []).append(p)


def exists_tracked(path, rel_to):
    cand = os.path.normpath(os.path.join(os.path.dirname(rel_to), path))
    return cand in tracked_set or path in tracked_set or any(
        t.endswith("/" + path) for t in tracked_set) or (os.sep not in path and path in basenames)


claude_files = [p for p in tracked if os.path.basename(p) == "CLAUDE.md"]
prose_files = claude_files + [p for p in ("README.md", "TODO.md") if p in tracked_set]

# 1. budgets and long lines
for f in claude_files:
    with open(f, encoding="utf-8") as fh:
        lines = fh.read().split("\n")
    budget = BUDGET_ROOT if f == "CLAUDE.md" else BUDGET_SUB
    if len(lines) > budget:
        errors.append(f"{f}: {len(lines)} lines, budget {budget}")
    for i, l in enumerate(lines, 1):
        if len(l) > MAX_LINE:
            errors.append(f"{f}:{i}: line of {len(l)} characters (max {MAX_LINE})")
        if DATED_HEADING_RE.match(l):
            warnings.append(f"{f}:{i}: dated or status heading ({l.strip()[:70]})")

# 2. root markdown whitelist
for p in tracked:
    if "/" not in p and p.endswith(".md") and p not in ROOT_MD_WHITELIST:
        errors.append(f"{p}: markdown file in the repository root (allowed: {', '.join(sorted(ROOT_MD_WHITELIST))})")
for p in git("ls-files", "--others", "--exclude-standard"):
    if "/" not in p and p.endswith(".md"):
        warnings.append(f"{p}: untracked markdown file in the repository root")
    if BACKUP_RE.search(p):
        warnings.append(f"{p}: untracked backup file")

# 3. tracked backup files
for p in tracked:
    if BACKUP_RE.search(p):
        errors.append(f"{p}: tracked backup file (delete it, git keeps the history)")

# 4. relative links
LINK_RE = re.compile(r"\]\(([^)#\s]+)(?:#[^)]*)?\)")
link_scope = [p for p in tracked if p.endswith(".md") and (
    "/" not in p or os.path.basename(p) == "CLAUDE.md" or "--all-docs" in sys.argv
    or p.startswith(("docs/archive/", "docs/changelog/")) and False)]
if "--all-docs" not in sys.argv:
    link_scope += [p for p in tracked if p.endswith("_NOTES_2026-10.md") or p.startswith("docs/changelog/")
                   or p in ("docs/KNOWN_ISSUES_ARCHIVE.md", "docs/CAPABILITY_NOTES_2026.md", "docs/TODO_ARCHIVE_2026-10.md",
                            "docs/USAGE_TOOLS.md", "docs/README_METHOD_NOTES_2026.md", "docs/LOGGING_SYSTEM.md")]
for f in sorted(set(link_scope)):
    in_fence = False
    with open(f, encoding="utf-8") as fh:
        for n, line in enumerate(fh, 1):
            if line.lstrip().startswith("```"):
                in_fence = not in_fence
            if in_fence:
                continue
            for target in LINK_RE.findall(line):
                if target.startswith(("http://", "https://", "mailto:")):
                    continue
                cand = os.path.normpath(os.path.join(os.path.dirname(f), target))
                if cand not in tracked_set and not any(t.startswith(cand + "/") for t in tracked_set):
                    errors.append(f"{f}:{n}: link target not tracked: {target}")

# 5. backticked paths
for f in prose_files:
    in_fence = False
    with open(f, encoding="utf-8") as fh:
        for n, line in enumerate(fh, 1):
            if line.lstrip().startswith("```"):
                in_fence = not in_fence
            if in_fence:
                continue
            for p in PATH_RE.findall(line):
                p = re.sub(r"/\.\w+$", "", p)      # notation `x.h/.cpp` -> `x.h`
                if UNTRACKED_OK.search(p) or "*" in p or "<" in p:
                    continue
                if not exists_tracked(p, f):
                    warnings.append(f"{f}:{n}: `{p}` matches no tracked file")

for m in errors:
    print("ERROR  ", m)
for m in warnings:
    print("WARNING", m)
print(f"check_docs: {len(errors)} error(s), {len(warnings)} warning(s)")
sys.exit(1 if errors or ("--strict" in sys.argv and warnings) else 0)
