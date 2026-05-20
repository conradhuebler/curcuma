# CLI Round-Trip: Export and Replay Curcuma Runs

This page describes how to capture a full Curcuma run as a single JSON file and replay it later — with optional CLI overrides on the replay.

The mechanism is built on three pieces that work together:

1. **Flat-flag auto-routing** — any parameter registered via the `PARAM` macros is reachable by its flat name on the CLI (`-cn_cutoff_bohr 5.5`). The Parameter Registry knows which module owns it and routes accordingly.
2. **`-export_run <file.json>`** — writes the resolved controller after all parsing and merging, plus the meta fields needed to replay the run.
3. **`-import_config <file.json>`** — performs a recursive deep merge of the JSON onto the (possibly empty) controller built from the rest of the CLI. CLI values always win.

## Quickstart

Capture a run:

```bash
curcuma -sp water.xyz -method gfnff -cn_cutoff_bohr 5.5 -export_run my_run.json
```

The resulting `my_run.json` contains:

```json
{
  "_command": "sp",
  "_input": "water.xyz",
  "method": "gfnff",
  "gfnff": {
    "accuracy": "normal",
    "cn_cutoff_bohr": 5.5,
    "eeq_max_iterations": -1,
    ...
  },
  "opt": { ... },
  "verbosity": 0,
  ...
}
```

Replay it:

```bash
curcuma -import_config my_run.json
```

Replay with one override:

```bash
curcuma -import_config my_run.json -cn_cutoff_bohr 7.0
```

CLI flags take precedence over JSON content at every nesting depth.

## Meta fields

| Field | Purpose |
|-------|---------|
| `_command` | The Curcuma command (without the leading `-`): `sp`, `opt`, `md`, `rmsd`, … |
| `_input` | The positional input file path (e.g. molecule geometry) |

Both are read when you invoke `curcuma -import_config file.json` without a CLI command. If you do pass a CLI command (`curcuma -opt mol.xyz -import_config ...`), the CLI wins and a warning is emitted noting that the JSON `_command` was ignored.

## Auto-routing in practice

The registry-driven routing means you can write:

```bash
curcuma -opt mol.xyz -method gfnff -cn_cutoff_bohr 5.5 -static_charges true
```

`-cn_cutoff_bohr` is a registered `gfnff` parameter, so it lands in `controller["gfnff"]["cn_cutoff_bohr"]` automatically. `-static_charges` is not (currently) registered, so it stays in the active command's module — write `-gfnff.static_charges true` to target gfnff explicitly.

Ambiguous names (the same `PARAM` name registered in two different modules) emit a warning and require dotted disambiguation.

## What `-export_run` writes

- The resolved controller after CLI parsing + any `-import_config` merge.
- `_command` from the active CLI keyword.
- `_input` from positional `argv[2]` (when not a flag).
- Registry defaults for every touched module, so the exported file documents every parameter and its effective value. Explicit values always win over defaults.
- Meta-flags (`import_config`, `export_run`) are stripped from the output to keep replay clean.

## What `-import_config` does

- Loads the JSON file.
- Recursively deep-merges it onto the controller. At every depth: CLI-provided values are preserved; missing keys are filled from the JSON.
- If the merge brings in `_command` or `_input` and the CLI didn't already supply them, those drive dispatch on this invocation.

## Compatibility

- Dotted notation (`-gfnff.static_charges true`) still works.
- Pre-existing CLI invocations are unchanged: unregistered flags stay in the command module exactly as before.
- The shallow 1-level merge previously used by `-import_config` is replaced by a recursive deep merge; nested overrides at any depth now work correctly.
