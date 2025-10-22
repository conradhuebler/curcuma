# Known Bugs in CLI Tests

## curcumaopt: JSON null-Fehler in Optimierung

**Entdeckt**: 2025-10-19
**Status**: NICHT BEHOBEN - Dokumentiert für zukünftige Behebung

### Symptome
```
[ERROR] Optimization failed with lbfgspp: LBFGSpp optimization failed:
[json.exception.type_error.306] cannot use value() with null
```

### Betroffene Tests
- ❌ `cli_curcumaopt_02_gfn2_single_point` - FEHLGESCHLAGEN
- ❌ `cli_curcumaopt_03_invalid_method` - FEHLGESCHLAGEN
- ❌ `cli_curcumaopt_04_lbfgs_params` - FEHLGESCHLAGEN
- ❌ `cli_curcumaopt_05_alias_singlepoint` - FEHLGESCHLAGEN
- ❌ `cli_curcumaopt_06_opth_only` - FEHLGESCHLAGEN
- ✅ `cli_curcumaopt_01_default_uff` - BESTANDEN (Warnung: prüft nur Datei-Existenz)

### Betroffene Module
- **curcumaopt** (Hauptmodul)
- **simplemd** (nutzt möglicherweise denselben Optimierer):
  - ❌ `cli_simplemd_01_nve_short` - Trajektorie nicht erstellt
  - ❌ `cli_simplemd_02_nvt_berendsen` - Trajektorie nicht erstellt
  - ❌ `cli_simplemd_03_invalid_thermostat` - Fehlgeschlagen
  - ❌ `cli_simplemd_04_rattle` - Trajektorie nicht erstellt
  - ❌ `cli_simplemd_05_wall` - Trajektorie nicht erstellt
  - ❌ `cli_simplemd_07_restart` - Restart-Datei nicht erstellt

### Root Cause (Vermutung)
- **ConfigManager Migration**: Möglicherweise fehlen Parameter-Defaults in ParameterRegistry
- **JSON null value()**: Ein Parameter wird mit `.value()` abgefragt, ist aber `null`
- **Wahrscheinliche Quelle**: `src/capabilities/optimiser/` oder `curcumaopt.cpp`

### Reproduktion
```bash
cd /home/conrad/src/curcuma/test_cases/cli/curcumaopt/01_default_uff_opt
/home/conrad/src/curcuma/release/curcuma -opt input.xyz -method uff -verbosity 2
```

### Workaround
Keine. Test `cli_curcumaopt_01` besteht nur, weil er Datei-Existenz prüft, nicht Erfolg.

### Nächste Schritte (für zukünftige Behebung)
1. Prüfe `src/capabilities/curcumaopt.cpp` auf `.value()` Aufrufe
2. Vergleiche mit `src/capabilities/rmsd.cpp` (funktioniert korrekt)
3. Prüfe ParameterRegistry für fehlende `opt` Parameter
4. Debug mit GDB oder erweiterten Logging

---

## Funktionierende Tests (Stand 2025-10-19)

### RMSD Tests
- ✅ `cli_rmsd_01_default` - BESTANDEN
- ✅ `cli_rmsd_02_no_reorder` - BESTANDEN
- ✅ `cli_rmsd_04_template_elements` - BESTANDEN
- ✅ `cli_rmsd_05_fragment` - BESTANDEN
- ✅ `cli_rmsd_06_alias_reorder` - BESTANDEN

### ConfScan Tests
- ✅ `cli_confscan_01_default` - BESTANDEN
- ✅ `cli_confscan_02_get_rmsd` - BESTANDEN
- ✅ `cli_confscan_04_slx_logic` - BESTANDEN
- ✅ `cli_confscan_05_hybrid_elements` - BESTANDEN
- ✅ `cli_confscan_06_heavy_only` - BESTANDEN
- ✅ `cli_confscan_07_restart` - BESTANDEN

### SimpleMD Tests
- ✅ `cli_simplemd_06_alias_temperature` - BESTANDEN

**Gesamt: 13/26 Tests bestanden (50%)**
