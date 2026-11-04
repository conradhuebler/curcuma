# Technical-debt analysis, Nov 4, 2026

Three read-only AI analyses (parameter routing, code level, build/tests/docs). Nothing was built or run
except where a report marks a claim [run]; all other statements are grep or code reading and are
unverified until a `release/` build confirms them. Line numbers refer to commit c3f570f.

- [PARAMS.md](PARAMS.md) - parameter routing, proposed unified design, migration steps
- [CODE.md](CODE.md) - dead code, duplicates, Unicode output, output paths, safe deletions
- [BUILD_DOCS.md](BUILD_DOCS.md) - CI, CMake, tests, doc drift, repo hygiene

Done on branch `claude/brave-lovelace-t22pd0`: strict number parsing and `-dock`/`-traj` routing in `CLI2Json`,
five unreferenced headers removed, CI smoke tests with pipefail. Everything else in the reports is open.
