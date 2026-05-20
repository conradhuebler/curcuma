# WP-Build-Ulysses-Cleanup — Ulysses-Reste aus CMakeLists entfernen

**Status:** 🆕 Vorgeschlagen (Mai 2026)
**Aufwand:** ~1 Stunde
**Erwarteter Nutzen:** Build-System-Hygiene; eliminiert verwirrenden `ulysses_lib`-Target, der mit `USE_ULYSSES=OFF` trotzdem existiert und nur EHT enthält.
**Hebel:** keine Laufzeit-Verbesserung; ausschließlich Code-Klarheit + ein wiederkehrender Build-Fehler.

## Hypothese

`CMakeLists.txt:540-574` definiert `ulysses_lib` **immer** als Library-Target, auch wenn `USE_ULYSSES=OFF`. In diesem Fall enthält die Library nur `src/core/energy_calculators/qm_methods/eht.cpp` — also Curcuma-eigenes EHT, **keinen Ulysses-Code**.

Konsequenzen:
- Build-Output zeigt `Linking CXX static library libulysses_lib.a`, obwohl kein Ulysses dabei ist (verwirrend für neue Entwickler).
- Bei Branch-Wechseln entstehen 0-Byte-Object-Files in `CMakeFiles/ulysses_lib.dir/` (heute beobachtet: `eht.cpp.o = 0 bytes`), was zu Linker-Fehlern `archive has no index; run ranlib` führt.
- `target_link_libraries(curcuma_core km ulysses_lib)` zwingt `curcuma_core` an einen falsch benannten Target zu linken.

## Aufgabe

### 1. EHT in `curcuma_core` integrieren

`src/core/energy_calculators/qm_methods/eht.cpp` ist eigenständiger Curcuma-Code (Extended Hückel Theory). Es gehört zu `curcuma_core`, nicht in einen "ulysses_lib"-Wrapper.

```cmake
# CMakeLists.txt — Änderung in der curcuma_core-Sourceliste

set(curcuma_core_SRC
    ...
    src/core/energy_calculators/qm_methods/eht.cpp   # <-- direkt in core
    ...
)
```

### 2. `ulysses_lib`-Target an `USE_ULYSSES` koppeln

```cmake
if(USE_ULYSSES)
    set(ulysess_SRC
        src/core/energy_calculators/qm_methods/interface/ulysses.cpp
    )
    add_library(ulysses_lib ${ulysess_SRC})
    target_compile_definitions(ulysses_lib PRIVATE USE_ULYSSES)

    if(USE_BLAS) target_link_libraries(ulysses_lib blas) endif()
    if(USE_MKL)  target_link_libraries(ulysses_lib -m64 -I${MKLROOT}/include ...) endif()

    add_executable(ulysses_helper src/helpers/ulysses_helper.cpp)
    target_link_libraries(ulysses_helper PUBLIC ulysses_lib curcuma_core)

    target_link_libraries(curcuma_core ulysses_lib)
endif()
```

Damit existiert `ulysses_lib` nur, wenn auch tatsächlich Ulysses gebaut wird.

### 3. Alte Zeile entfernen

```cmake
# REMOVE
target_link_libraries(curcuma_core km ulysses_lib)
# REPLACE WITH
target_link_libraries(curcuma_core km)
if(USE_ULYSSES)
    target_link_libraries(curcuma_core ulysses_lib)
endif()
```

### 4. CMake-Cache-Migration

Existing Build-Verzeichnisse haben `CMakeFiles/ulysses_lib.dir/eht.cpp.o`. Nach dem Patch sollte `make clean && make` oder gar `rm -rf release/CMakeFiles release/CMakeCache.txt && cmake ..` ausgeführt werden, um stale Targets zu eliminieren.

### Code-Anker

| Datei | Position | Aufgabe |
|-------|----------|---------|
| `CMakeLists.txt` | L540-574 | `ulysses_lib`-Block in `if(USE_ULYSSES)` einschließen |
| `CMakeLists.txt` | curcuma_core SRC-Liste | `eht.cpp` aufnehmen |
| `CMakeLists.txt` | `target_link_libraries(curcuma_core km ulysses_lib)` | konditional machen |

## Akzeptanzkriterien

1. ⚙️ `cmake -DUSE_ULYSSES=OFF .. && make -j4` baut ohne `ulysses_lib`-Target und ohne Verwirrung.
2. ⚙️ `cmake -DUSE_ULYSSES=ON .. && make -j4` baut Ulysses normal mit `ulysses_lib` + `ulysses_helper`.
3. ⚙️ EHT-funktional unverändert (`curcuma -sp molecule.xyz -method eht` läuft).
4. ⚙️ Keine 0-Byte-Object-Files mehr nach Branch-Wechseln.
5. ⚙️ `grep ulysses_lib release/CMakeFiles/curcuma_core.dir/` zeigt mit USE_ULYSSES=OFF: nichts.

## Risiken

1. **EHT-Dependencies** — wenn `eht.cpp` versteckte Abhängigkeiten auf alte `ulysses_lib`-Header hat, muss das identifiziert + entkoppelt werden. Ein Test-Build wird das zeigen.
2. **PCH (Precompiled Headers)** — falls `eht.cpp` Teil eines PCH-Setups wird, müssen die `target_precompile_headers(curcuma REUSE_FROM curcuma_core)`-Direktiven konsistent bleiben.
3. **Migration-Stolperstein** — bestehende Build-Verzeichnisse brauchen `rm -rf CMakeFiles/ulysses_lib.dir` vor `make`, sonst stale objects. Im PR-Beschreibung erwähnen.

## Verifikation

```bash
# Clean build mit Ulysses OFF
rm -rf release && mkdir release && cd release
cmake -DUSE_ULYSSES=OFF -DUSE_GFNFF_FAST_EXP=OFF ..
make -j4 2>&1 | grep -i "ulysses" | head -3   # erwartet: keine Treffer
./curcuma --version
./curcuma -sp ../test_cases/molecules/larger/caffeine.xyz -method eht
# erwartet: läuft normal

# Optional: Ulysses ON
cd .. && rm -rf release && mkdir release && cd release
cmake -DUSE_ULYSSES=ON ..
make -j4 ulysses_lib ulysses_helper
# erwartet: beide Targets existieren

# Regression
ctest -R cli_ --output-on-failure | tail -5
```

## Beziehung zu anderen WPs

- **Eigenständig**, keine Abhängigkeit. Kann jederzeit gemacht werden.
- Reduziert Reibungsfläche für andere WPs, die `forcefieldthread.cpp` oder `eeq_solver.cpp` ändern — saubere Build-Targets erleichtern Inkremental-Builds.
