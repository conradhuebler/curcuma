# Output symbols

Output strings spell symbols through macros from `src/core/output_symbols.h`:

```cpp
CurcumaLogger::info_fmt("angle {:.2f}" SYM_DEG, angle);   // " deg" or the degree sign
```

- Default build: ASCII spelling (`sigma`, `+/-`, `^2`, `Delta`, `pi`, ` deg`, `to`).
- `cmake -DCURCUMA_UNICODE_OUTPUT=ON`: the Unicode characters (sigma, plus-minus, superscripts, Greek letters, degree sign, arrows).
- The choice is made at compile time. There is no environment variable or command line switch: the texts are format strings fixed in the code.
- Letters with diacritics (u-umlaut, e-acute, A with ring above, ...) are written as UTF-8 in the source and appear in both builds.
- PARAM help texts are plain ASCII symbols, because the parameter generator reads them as plain string literals.
- Affects terminal output, CSV and gnuplot files written by `-analysis` and the trajectory tools, and the `"unit"` field of JSON output (`" deg"` versus the degree sign).

Both variants were built and compared on `-rmsd`, `-analysis`, `-sp`, `-hessian` and `-md`: the ASCII build prints only A with ring and u-umlaut outside ASCII, the Unicode build additionally the degree sign, Delta, pi and superscripts. The `ctest` suite was run with the ASCII build only.
