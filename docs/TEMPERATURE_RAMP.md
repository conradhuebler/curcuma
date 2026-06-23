# Temperature Ramps, Live Control, and Thermal Regions (SimpleMD)

> Status: AI-implemented, machine-tested (Jun 2026). Human production testing pending.
> The regional thermostat uses a `3*N` per-region degree-of-freedom count (inter-region
> constraints / COM removal are not subtracted) — validate before quantitative use.

`SimpleMD` can change the thermostat target temperature **during** a run, drive it through a
**multi-stage schedule** (ramp), and thermostat **atom subsets** (regions) to their own targets.
All three are configurable from the curcuma CLI/JSON and from the qurcuma GUI.

## 1. Live temperature

`SimpleMD::setTargetTemperature(double T)` sets the global thermostat setpoint `m_T0` between
steps (the thermostats read `m_T0` every step). It marks the run as manually overridden, so an
active **global** ramp stops touching `m_T0` for the rest of the run (manual control wins). The
qurcuma simulation dock drives this from the vertical temperature slider while the MD runs.

## 2. Global temperature ramp

Two parameters (PARAM category "Temperature Ramp"):

| Parameter       | Type   | Meaning |
|-----------------|--------|---------|
| `temp_ramp`     | bool   | Enable the multi-stage schedule. |
| `temp_schedule` | string | The schedule (grammar below). |

Grammar: `target:mode:value [; target:mode:value ...]`

- `target` — segment target temperature in Kelvin.
- `mode` — `steps` or `reach`:
  - `steps`: linearly ramp the setpoint from the segment's start value to `target` over `value`
    integration steps, then advance.
  - `reach`: hold the setpoint at `target`, advance once `|<T> - target| < value` Kelvin
    (running-mean temperature, after a short warm-up).
- `value` — step count (`steps`) or tolerance in K (`reach`).

When the schedule finishes, the last setpoint is held. A malformed schedule disables the ramp
(constant-T run) with a warning.

```bash
# Heat 300 -> 600 over 5000 steps, hold 2000 steps, then cool to 300 once equilibrated:
curcuma -md water.xyz -method gfnff -temperature 300 \
        -temp_ramp true -temp_schedule "600:steps:5000;600:steps:2000;300:reach:10"
```

The parameters round-trip through `-export_run` / `-import_config` like any other PARAM.

## 3. Thermal regions

A region thermostats an atom subset to its own (optionally ramped) target temperature. Atoms in
no region follow the global setpoint. Regions are a JSON array `simplemd.temp_regions` (best given
via a JSON input file / the GUI, since arrays are not flat CLI flags):

```json
{
  "_command": "md", "_input": "system.xyz",
  "simplemd": {
    "method": "gfnff", "temperature": 300, "thermostat": "csvr", "max_time": 5000,
    "temp_regions": [
      { "atoms": "1:50",   "temperature": 600, "temp_schedule": "800:steps:5000" },
      { "atoms": "51:100", "temperature": 100 }
    ]
  }
}
```

```bash
curcuma -import_config system.json
```

- `atoms` — selection string (`Molecule::FragString2Indicies` grammar): comma list `1,3,5`,
  ranges `1:10`, fragments `F2`, or `-1` for all (1-based indices).
- `temperature` — region start setpoint (K).
- `temp_schedule` — optional per-region ramp (same grammar as the global ramp).

Overlapping selections are resolved first-region-wins. The realized per-region temperature
converges to the DOF-weighted mix of the region targets.

### Thermostat support

| Thermostat   | Regional support |
|--------------|------------------|
| Berendsen    | yes |
| CSVR         | yes |
| Anderson     | yes (naturally per-atom) |
| Nosé-Hoover  | no — the global chain is applied to all atoms (region targets ignored, warned) |
| none         | no thermostatting |

With **no** `temp_regions` defined the thermostat path is unchanged (byte-identical to the legacy
single-thermostat behaviour), so existing runs are unaffected.

## Implementation

`external/curcuma/src/capabilities/simplemd.{h,cpp}`:
`setTargetTemperature`/`targetTemperature`, `ParseSchedule`/`StepRamp`/`UpdateTemperatureRamp`
(called at the top of `step()`), `ParseThermalRegions`/`ResolveThermalRegions`,
`ApplyThermostat`/`ApplyThermostatRegion`, `RegionTemperature`.
