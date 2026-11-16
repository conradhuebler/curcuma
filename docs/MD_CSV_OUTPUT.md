# MD CSV Output — Per-Step Energies & Temperatures

**Status**: 🤖 AI-generated (September 2026) — ⚙️ machine-tested, human production
testing pending.

`simplemd` writes the same per-step status data it prints to the console
(energies, temperatures, wall/virial terms) into a CSV file
`<basename>.md.csv` in the BMT output directory. The file is UTF-8 encoded
(with BOM) and uses a configurable delimiter, so it can be opened directly
in common spreadsheet programs (Excel, LibreOffice Calc).

## Activation

Enabled by default. Controlled by two PARAMs:

| PARAM | Type | Default | Description |
|-------|------|---------|-------------|
| `write_csv` | Bool | `true` | Write the CSV status file. |
| `csv_delimiter` | String | `;` | Field delimiter (`;`, `,`, or `\t` for tab). |

```bash
# Default run (semicolon-delimited CSV in the BMT dir)
./curcuma -md mol.xyz -method uff
# → mol.md.<timestamp>/mol.md.csv

# Comma-delimited instead
./curcuma -md mol.xyz -method uff -csv_delimiter ","
```

## Format

One row per console status line (frequency follows `print_frequency`).
The first line is the header; the file starts with a UTF-8 BOM
(`EF BB BF`) so Excel detects the encoding reliably.

Columns (semicolon-delimited example):

```
step;time_ps;epot;epot_avg;ekin;ekin_avg;etot;etot_avg;temperature;temperature_avg;wall_potential;wall_potential_avg;virial_correction;virial_correction_avg;remaining;dt
0;0.000000;-0.123456;-0.123456;0.010000;0.010000;-0.113456;-0.113456;298.150000;298.150000;0.000000;0.000000;0.000000;0.000000;0.000000;1.000000
```

| Column | Unit | Meaning |
|--------|------|---------|
| `step` | – | Integration step index |
| `time_ps` | ps | Simulation time |
| `epot` / `epot_avg` | Eh | Potential energy / running average |
| `ekin` / `ekin_avg` | Eh | Kinetic energy / running average |
| `etot` / `etot_avg` | Eh | Total energy / running average |
| `temperature` / `temperature_avg` | K | Instantaneous / averaged temperature |
| `wall_potential` / `wall_potential_avg` | Eh | Wall potential / average |
| `virial_correction` / `virial_correction_avg` | Eh | Virial correction / average |
| `remaining` | s or min | Estimated remaining wall time |
| `dt` | ps | Integration time step |

Optional appendix columns appear when the corresponding feature is active:
`dipole_debye` (with `-dipole`) and `n_unique` (with `-unique`).

## Notes

- The console output is unchanged; the CSV is written in addition.
- The file is written via `outputPath()`, so it lands in the BMT directory
  (`<basename>.md.<timestamp>/`) and in the current working directory if BMT
  is disabled.
- Numbers use a dot as decimal separator (C++ default). With the default
  semicolon delimiter this imports correctly into German-locale Excel /
  LibreOffice; with a comma delimiter, German-locale spreadsheets would
  misinterpret the dot as a thousands separator.
