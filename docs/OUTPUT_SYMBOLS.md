# Output symbols

Output strings spell symbols through macros from `src/core/output_symbols.h`:

```cpp
CurcumaLogger::info_fmt("angle {:.2f}" SYM_DEG, angle);   // " deg" or the degree sign
```

- Default build: the Unicode characters (sigma, plus-minus, superscripts, Greek letters, degree sign, arrows).
- `cmake -DCURCUMA_UNICODE_OUTPUT=OFF`: the ASCII spelling (`sigma`, `+/-`, `^2`, `Delta`, `pi`, ` deg`, `to`), for terminals or log parsers that cannot handle UTF-8.
- The choice is made at compile time. There is no environment variable or command line switch: the texts are format strings fixed in the code.
- Letters with diacritics (u-umlaut, e-acute, A with ring above, ...) are written as UTF-8 in the source and appear in both builds.
- PARAM help texts use the same macros between the string literals, for example `"Spring constant in Eh/Å" SYM_SQ`. The parameter generator (`scripts/param_parser/main.cpp`) passes `SYM_*` names through into `generated/parameter_registry.h`, which includes `src/core/output_symbols.h`; any other token in the help argument is rejected.
- Affects terminal output, CSV and gnuplot files written by `-analysis` and the trajectory tools, and the `"unit"` field of JSON output (`" deg"` versus the degree sign).

Both variants were built and compared on `-rmsd`, `-analysis`, `-sp`, `-hessian` and `-md`: the ASCII build prints only A with ring and u-umlaut outside ASCII, the Unicode build additionally the degree sign, Delta, pi and superscripts. The The full `ctest` suite passes in both builds with the same three failures (`xtb_cpscf` and two tests that also fail on plain `devel/methods`); no test depends on the symbol spelling. Windows consoles with a legacy code page may need `-DCURCUMA_UNICODE_OUTPUT=OFF`; that was not tested.
