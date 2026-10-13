# Streuungs- und Strukturanalyse - Benutzerhandbuch

**Für Nicht-Experten**: Einfache Anleitung zur Analyse von Molekülstrukturen und CG-Systemen

---

## Was macht diese Analyse?

Die Scattering-Analyse beantwortet drei Hauptfragen über Ihre Moleküle:

1. **Wie groß und welche Form?** → Formfaktor P(q) & Gyration Radius
2. **Wie sind sie angeordnet?** → Strukturfaktor S(q) & g(r)
3. **Welche Form haben sie?** → Shape-Parameter (flach, kugelförmig, stabförmig)

---

## 📊 Die drei Analysetypen

### 1. **Scattering (Streuung)** - Größe und Form messen

**Was es tut**: Simuliert ein Streuexperiment (wie SAXS/SANS im Labor)

**Wichtige Ergebnisse**:
- **P(q)** - Formfaktor: Wie das einzelne Molekül aussieht
- **S(q)** - Strukturfaktor: Wie die Moleküle zueinander stehen
- **Guinier Rg** - Gyration Radius: Ungefähre Größe des Moleküls

**Interpretation**:
```
P(q) bei kleinem q (→0): Höhe zeigt Molekülgröße
P(q) Abfall: Je schneller, desto kompakter das Molekül

S(q) = 1: Keine Ordnung (ideales Gas)
S(q) > 1: Moleküle klumpen zusammen
S(q) < 1: Moleküle meiden sich
```

---

### 2. **RDF (g(r))** - Radiale Verteilungsfunktion

**Was es tut**: Zählt, wie viele Nachbarn in welchem Abstand sind

**Wichtige Ergebnisse**:
- **Erster Peak**: Typischer Abstand der nächsten Nachbarn
- **Peak-Höhe**: Wie stark diese Ordnung ist
- **Koordinationszahl**: Wie viele Nachbarn im Schnitt

**Interpretation**:
```
g(r) = 0: Kein Atom in diesem Abstand
g(r) = 1: Durchschnittliche Dichte (keine Ordnung)
g(r) > 1: Mehr Atome als Durchschnitt (Schale/Struktur)

Erster Peak oft bei: Bindungslänge (1-2 Å)
Zweiter Peak: Nächste Nachbarschale
```

---

### 3. **Shape-Parameter** - Molekülform beschreiben

**Was es tut**: Berechnet, ob das Molekül eher flach, stabförmig oder kugelförmig ist

**Wichtige Ergebnisse**:
- **Asphericity** (b):
  - b > 0: "Prolate" (Zigarre/Stab)
  - b < 0: "Oblate" (Scheibe/flach)
  - b ≈ 0: Kugelförmig

- **Acylindricity** (c): Wie stark asymmetrisch
- **Anisotropy** (κ²): 0 = Kugel, 1 = Stab/Scheibe

**Interpretation**:
```
Benzene: asphericity = -2.0 → Flach (Scheibe) ✓
DNA: asphericity > 0 → Stabförmig
Protein (kugelförmig): asphericity ≈ 0
```

---

## 🚀 Praktische Beispiele

### Beispiel 1: Schnelle Übersicht (CG-System)

```bash
curcuma -analysis mein_system.vtf \
    -scattering_enable true \
    -scattering_q_max 1.0 \
    -scattering_q_steps 50
```

**Was Sie bekommen**:
- P(q) und S(q) Werte von q=0.01 bis q=1.0 Å⁻¹
- Guinier Rg (Größe des Moleküls)
- Automatische Erkennung: CG-System → Kugel-Formfaktoren

**Typische Ausgabe verstehen**:
```json
"scattering": {
  "guinier_rg": 5.67,        ← Molekül ist ~5.7 Å groß
  "system_type": "cg",        ← Coarse-Grained erkannt
  "P_q": [1.0, 0.95, ...]    ← Formfaktor-Werte
}
```

---

### Beispiel 2: Vollständige Strukturanalyse

```bash
curcuma -analysis polymer.vtf \
    -scattering_enable true \
    -scattering_q_max 2.0 \
    -rdf_enable true \
    -rdf_r_max 20.0 \
    -rdf_coordination_shells true \
    -shape_asphericity true \
    -shape_anisotropy true
```

**Was Sie bekommen**:
- Scattering: P(q), S(q), Guinier Rg
- RDF: g(r), Peaks, Koordinationszahlen
- Shape: Ist es flach/stabförmig/kugelförmig?

**JSON-Output analysieren**:
```json
{
  "scattering": {
    "guinier_rg": 8.5,           ← Größe
    "P_q": [...],
    "S_q": [...]
  },
  "rdf": {
    "first_peak_position": 3.8,  ← Nächster Nachbar bei 3.8 Å
    "first_peak_height": 2.45,   ← Starke Ordnung
    "coordination_number": [...] ← Wie viele Nachbarn
  },
  "shape_descriptors": {
    "asphericity": 0.123,        ← Leicht stabförmig
    "shape_type": "prolate",     ← Zigarre
    "relative_anisotropy": 0.23  ← Mäßig asymmetrisch
  }
}
```

---

### Beispiel 3: Trajektorien-Analyse (viele Frames)

```bash
curcuma -analysis md_simulation.xyz \
    -scattering_enable true \
    -rdf_enable true \
    -statistics cumulative
```

**Was passiert**:
1. Alle Frames werden automatisch durchlaufen
2. P(q), S(q), g(r) für jeden Frame berechnet
3. Mittelwert und Standardabweichung über alle Frames

**Typischer Output**:
```
Processed 100 structures with 3 metrics selected
```

Dann bekommen Sie:
- Mittelwerte über alle Frames
- Standardabweichungen (wie stark schwankt es?)
- Sie sehen, ob die Struktur stabil ist

---

## 📐 Parameter-Übersicht

### Scattering-Parameter

| Parameter | Default | Beschreibung | Typische Werte |
|-----------|---------|--------------|----------------|
| `scattering_enable` | false | Scattering aktivieren | true/false |
| `scattering_q_min` | 0.01 | Kleinster q-Wert (Å⁻¹) | 0.01 - 0.1 |
| `scattering_q_max` | 2.0 | Größter q-Wert (Å⁻¹) | 0.5 - 10.0 |
| `scattering_q_steps` | 100 | Anzahl q-Punkte | 20 - 200 |
| `scattering_form_factor` | auto | Formfaktor-Typ | auto/cromer_mann/cg_sphere |
| `scattering_cg_radius` | 3.0 | CG-Bead Radius (Å) | 2.0 - 5.0 |

**Empfehlung für verschiedene Systeme**:
```
Kleine Moleküle (< 50 Atome):
  -scattering_q_max 5.0 -scattering_q_steps 100

CG-Polymere:
  -scattering_q_max 1.0 -scattering_q_steps 50

Proteine:
  -scattering_q_max 2.0 -scattering_q_steps 100
```

---

### Q-Spacing Modes (NEU 2026) - Claude Generated

**Logarithmische q-Werte (Default)**: Bessere Auflösung im Guinier-Bereich

| Parameter | Default | Beschreibung | Optionen |
|-----------|---------|--------------|----------|
| `scattering_q_spacing` | log | Q-Wert-Verteilung | log / linear |

**Warum logarithmisch?**
- **10x bessere Auflösung** bei kleinen q-Werten (Guinier-Region)
- Entspricht **SAXS/SANS-Detektorgeometrie** in Experimenten
- **Ideal für log-log-Plots** (Standard in der Literatur)

**Formel**:
```
Logarithmisch: q[i] = q_min × (q_max/q_min)^(i/(N-1))
Linear:        q[i] = q_min + i × (q_max - q_min)/(N-1)
```

**Beispiel**:
```bash
# Logarithmisch (default, empfohlen für Guinier-Analyse)
curcuma -analysis input.xyz -scattering_enable

# Explizit linear (für direkte Abstandsraum-Interpretation)
curcuma -analysis input.xyz -scattering_enable \
    -scattering_q_spacing linear
```

**Vergleich**:
```
20 Punkte von 0.01 bis 2.0 Å⁻¹:

Logarithmisch:  0.010, 0.013, 0.017, 0.023, ... (dicht bei kleinen q)
Linear:         0.010, 0.115, 0.220, 0.325, ... (gleichmäßig verteilt)

→ Logarithmisch hat 8 Punkte < 0.1 Å⁻¹ (Guinier wichtig!)
→ Linear hat nur 1 Punkt < 0.1 Å⁻¹
```

---

### Automatische Gnuplot-Visualisierung (NEU 2026) - Claude Generated

**Curcuma erzeugt automatisch ein Gnuplot-Script** für jede Scattering-Analyse!

**Generierte Dateien**:
- `basename.scattering_statistics.csv` - Statistische Daten (P_avg, P_std, S_avg, S_std)
- `basename.scattering.gnu` - Gnuplot-Script
- `basename.scattering_plot.png` - Plot (nach `gnuplot` Ausführung)

**Verwendung**:
```bash
# 1. Analyse durchführen
curcuma -analysis input.xyz -scattering_enable -o results.json

# 2. Plot erzeugen
gnuplot input.scattering.gnu

# 3. Öffnen
xdg-open input.scattering_plot.png  # Linux
open input.scattering_plot.png      # macOS
```

**Plot-Layout** (2×2 Grid):
1. **P(q) log-log**: Form Factor (Guinier-Region sichtbar)
2. **S(q) mit Referenz**: Structure Factor vs. ideales Gas (S=1)
3. **P(q) + S(q) kombiniert**: Beide Faktoren im Vergleich
4. **Guinier-Analyse**: ln(P) vs. q² (Rg-Extraktion)

**Anpassung**:
```bash
# Script bearbeiten für eigene Farben/Größe
nano input.scattering.gnu

# Beispiel: Terminal auf SVG ändern
sed -i 's/pngcairo/svg/g' input.scattering.gnu
gnuplot input.scattering.gnu  # → SVG statt PNG
```

**Wissenschaftlicher Nutzen**:
- **Log-log-Plot**: Guinier-Plateau und Power-Law-Regime sofort sichtbar
- **Guinier-Plot**: Linearer Fit → Gyration Radius Rg
- **S(q) vs. 1**: Zeigt Abweichung vom idealen Gas (Aggregation/Repulsion)

---

### RDF-Parameter

| Parameter | Default | Beschreibung | Typische Werte |
|-----------|---------|--------------|----------------|
| `rdf_enable` | false | RDF aktivieren | true/false |
| `rdf_r_max` | 15.0 | Maximaler Abstand (Å) | 10.0 - 30.0 |
| `rdf_bin_width` | 0.05 | Histogram-Breite (Å) | 0.01 - 0.2 |
| `rdf_coordination_shells` | false | Koordinationszahlen | true/false |

**Empfehlung**:
```
Atomare Systeme:
  -rdf_r_max 15.0 -rdf_bin_width 0.05

CG-Systeme:
  -rdf_r_max 30.0 -rdf_bin_width 0.1
  -rdf_coordination_shells true
```

---

### Shape-Parameter

| Parameter | Default | Beschreibung | Bedeutung |
|-----------|---------|--------------|-----------|
| `shape_asphericity` | false | Asphericity berechnen | Flach vs. stabförmig |
| `shape_acylindricity` | false | Acylindricity berechnen | Symmetrie |
| `shape_anisotropy` | false | Anisotropie berechnen | Kugel vs. asymmetrisch |

**Empfehlung**: Alle aktivieren für vollständige Form-Analyse
```bash
-shape_asphericity true \
-shape_acylindricity true \
-shape_anisotropy true
```

---

## 💡 Häufige Fragen

### "Ich bekomme nur Zahlen - was bedeuten sie?"

**Lösung**: Nutzen Sie JSON-Output für Struktur:
```bash
curcuma -analysis mein_system.vtf \
    -scattering_enable true \
    -output_format json > results.json
```

Dann können Sie die Datei mit einem JSON-Viewer öffnen oder mit Python/Julia weiterverarbeiten.

---

### "Mein System ist sehr groß (1000+ Atome) - dauert ewig!"

**Lösung**: Reduzieren Sie q_steps und r_max:
```bash
-scattering_q_steps 30 \
-rdf_r_max 10.0
```

**Warum**:
- P(q) Berechnung: O(N² × q_steps)
- g(r) Berechnung: O(N²)

Für N=1000: ~1 Million Paar-Distanzen!

---

### "Was ist der Unterschied zwischen P(q) und S(q)?"

**P(q)** (Form factor):
- Beschreibt **einzelnes** Molekül
- "Wie sieht ein Teilchen aus?"
- Hängt nur von der Form ab

**S(q)** (Structure factor):
- Beschreibt **viele** Moleküle zueinander
- "Wie sind die Teilchen angeordnet?"
- Zeigt Ordnung/Flüssigkeitsstruktur

**Zusammenhang**: `I(q) = P(q) × S(q)` (Intensität im Experiment)

---

### "Mein CG-System hat verschiedene Bead-Größen!"

**Problem**: Aktuell wird nur ein einheitlicher Radius verwendet.

**Workaround**: Nutzen Sie den durchschnittlichen Radius:
```bash
-scattering_cg_radius 3.5  # Mittelwert Ihrer Beads
```

**Zukunft**: Element-spezifische Radien geplant (siehe TODO)

---

### "Kann ich nur g(r) ohne Scattering berechnen?"

**Ja!** Jede Analyse ist optional:
```bash
# Nur RDF
curcuma -analysis system.vtf -rdf_enable true

# Nur Shape
curcuma -analysis system.vtf -shape_asphericity true

# Nur Scattering
curcuma -analysis system.vtf -scattering_enable true
```

---

## 📈 Ergebnisse visualisieren

### Mit Python (Matplotlib)

```python
import json
import matplotlib.pyplot as plt

# JSON laden
with open('results.json') as f:
    data = json.load(f)

frame = data['timesteps'][0]

# P(q) plotten
if 'scattering' in frame:
    s = frame['scattering']
    plt.figure(figsize=(10, 4))

    plt.subplot(1, 2, 1)
    plt.loglog(s['q_values'], s['P_q'])
    plt.xlabel('q (Å⁻¹)')
    plt.ylabel('P(q)')
    plt.title('Form Factor')

    plt.subplot(1, 2, 2)
    plt.plot(s['q_values'], s['S_q'])
    plt.xlabel('q (Å⁻¹)')
    plt.ylabel('S(q)')
    plt.title('Structure Factor')

    plt.tight_layout()
    plt.savefig('scattering.png')

# g(r) plotten
if 'rdf' in frame:
    rdf = frame['rdf']
    plt.figure()
    plt.plot(rdf['r_values'], rdf['g_r'])
    plt.axhline(1.0, color='gray', linestyle='--', label='Ideal gas')
    plt.xlabel('r (Å)')
    plt.ylabel('g(r)')
    plt.title('Radial Distribution Function')
    plt.legend()
    plt.savefig('rdf.png')
```

---

### Mit gnuplot

```bash
# Extrahieren Sie die Daten zu CSV
curcuma -analysis system.vtf \
    -scattering_enable true \
    -output_format csv > scattering.csv

# Plotten
gnuplot <<EOF
set terminal png size 800,600
set output 'scattering.png'
set logscale y
set xlabel 'q (Å^{-1})'
set ylabel 'P(q)'
plot 'scattering.csv' using 1:2 with lines title 'Form Factor'
EOF
```

---

## 🎯 Schnellstart-Checkliste

1. **CG-System analysieren**:
   ```bash
   curcuma -analysis system.vtf \
       -scattering_enable true \
       -rdf_enable true \
       -output_format json > results.json
   ```

2. **Ergebnisse anschauen**:
   - `guinier_rg`: Wie groß ist mein Molekül?
   - `first_peak_position`: Abstand nächster Nachbarn
   - `shape_type`: Welche Form?

3. **Visualisieren**: Python/gnuplot (siehe oben)

4. **Bei Fragen**: Siehe "Häufige Fragen" oder CLAUDE.md

---

## 📚 Weiterführende Informationen

### Wissenschaftlicher Hintergrund

- **Form Factor**: International Tables for Crystallography, Vol. C
- **Guinier-Analyse**: A. Guinier, "La diffraction des rayons X aux très petits angles" (1939)
- **RDF**: Hansen & McDonald, "Theory of Simple Liquids" (2006)

### Curcuma-spezifische Dokumentation

- **Parameter-Details**: `PARAMETER_MIGRATION_GUIDE.md`
- **CG-Systeme**: `docs/archive/CG.md` (Plan von 2025, Stand dort vermerkt)
- **Entwicklung**: `CLAUDE.md`

---

**Erstellt**: 2026-01-06
**Autor**: Claude (unter Anleitung von Conrad Hübler)
**Version**: 1.0 (Erste Release mit Scattering-Analyse)
