# Technical-debt analysis, Nov 4, 2026

Three read-only AI analyses (parameter routing, code level, build/tests/docs). Nothing was built or run
except where a report marks a claim [run]; all other statements are grep or code reading and are
unverified until a `release/` build confirms them. Line numbers refer to commit c3f570f.

- [PARAMS.md](PARAMS.md) - parameter routing, proposed unified design, migration steps
- [CODE.md](CODE.md) - dead code, duplicates, Unicode output, output paths, safe deletions
- [BUILD_DOCS.md](BUILD_DOCS.md) - CI, CMake, tests, doc drift, repo hygiene

Done on branch `claude/brave-lovelace-t22pd0`: strict number parsing and `-dock`/`-traj` routing in `CLI2Json`,
five unreferenced headers removed, CI smoke tests with pipefail. Everything else in the reports is open.

Removed Nov 6, 2026 (S5 in CODE.md): 12 unreferenced GFNFF members and their declarations (-1164 lines in gfnff.h, gfnff_method.cpp, gfnff_torsions.cpp, gfnff_inversions.cpp). `NumGradFixedCharges` stays (used by test_gfnff_numgrad). Checked: build of curcuma, test_eeq_two_phase, test_gfnff_numgrad; `-method gfnff` energies of MOR41 (95) and GMTKN55 (2462) identical to the master binary of c3f570fb (max difference 0). Not checked: gfn1/gfn2, GPU builds, `ctest`.
