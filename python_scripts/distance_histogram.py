#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import argparse

def boltzmann_weight(energy, energy_ref, temperature):
    """
    Berechnet das Boltzmann-Gewicht basierend auf der Energie.
    
    Parameters:
    -----------
    energy : float
        Energie in Hartree (Eh)
    energy_ref : float
        Referenzenergie (niedrigste Energie im Datensatz)
    temperature : float
        Temperatur in Kelvin
    
    Returns:
    --------
    float : Boltzmann-Gewicht
    """
    # Konstanten
    k_B = 3.166811429e-6  # Boltzmann-Konstante in Eh/K
    
    # Berechnung des Boltzmann-Faktors
    delta_E = energy - energy_ref  # Energiedifferenz zur Referenzenergie
    if temperature == 0:
        return 1.0 if delta_E == 0 else 0.0
    else:
        return np.exp(-delta_E / (k_B * temperature))

def main():
    parser = argparse.ArgumentParser(description='Erstellt ein Histogramm von Abständen mit verschiedenen Wichtungsoptionen.')
    parser.add_argument('file', type=str, help='Pfad zur Eingabedatei')
    parser.add_argument('--mode', '-m', type=str, choices=['boltzmann', 'uniform', 'stable'], default='boltzmann',
                        help='Wichtungsmodus: "boltzmann" (Standard), "uniform" (keine Wichtung), "stable" (nur stabilste Struktur)')
    parser.add_argument('--temperature', '-T', type=float, default=300.0, 
                        help='Temperatur in Kelvin für Boltzmann-Wichtung (Standard: 300 K)')
    parser.add_argument('--bins', type=int, default=50, 
                        help='Anzahl der Bins im Histogramm (Standard: 50)')
    parser.add_argument('--output', '-o', type=str, default='histogramm.png',
                        help='Ausgabedatei für das Histogramm (Standard: histogramm.png)')
    parser.add_argument('--show', action='store_true', 
                        help='Histogramm nach dem Speichern anzeigen')
    args = parser.parse_args()
    
    # Daten aus der Datei lesen
    data = []
    try:
        with open(args.file, 'r') as f:
            for line in f:
                values = [float(val) for val in line.strip().split()]
                if len(values) > 1:  # Sicherstellen, dass die Zeile genügend Werte enthält
                    data.append(values)
    except Exception as e:
        print(f"Fehler beim Lesen der Datei: {e}")
        return
    
    if not data:
        print("Keine gültigen Daten gefunden.")
        return
    
    # Konvertiere in NumPy-Array für effizientere Operationen
    data_array = np.array(data)
    
    # Extrahiere Energien (erste Spalte)
    energies = data_array[:, 0]
    
    # Bestimme die niedrigste Energie als Referenz
    energy_ref = np.min(energies)
    
    # Bestimme die Gewichtung basierend auf der gewählten Methode
    mode = args.mode.lower()
    
    if mode == 'boltzmann':
        # Berechne Boltzmann-Gewichte für jede Konfiguration
        weights = np.array([boltzmann_weight(e, energy_ref, args.temperature) for e in energies])
        
        # Normalisiere die Gewichte
        if np.sum(weights) > 0:
            weights = weights / np.sum(weights)
        else:
            print("Warnung: Summe der Gewichte ist null. Überprüfen Sie die Temperatur und Energiedifferenzen.")
            return
        
        title = f'Histogramm of the distances with boltzman weights (T = {args.temperature} K)'
        
    elif mode == 'uniform':
        # Gleiche Gewichtung für alle Strukturen
        weights = np.ones(len(energies)) / len(energies)
        title = 'Histogramm of the distances with uniform weights'
        
    elif mode == 'stable':
        # Finde die stabilste Struktur (mit der niedrigsten Energie)
        min_energy_index = np.argmin(energies)
        
        # Setze das Gewicht für die stabilste Struktur auf 1, alle anderen auf 0
        weights = np.zeros(len(energies))
        weights[min_energy_index] = 1.0
        
        title = 'Histogramm of the distances in the stable structure'
    print(f"Erstelle Histogramm im Modus: {mode} mit {len(energies)} Strukturen.")
    print(weights)
    # Extrahiere alle Abstände
    if mode == 'stable':
        # Nur die Abstände der stabilsten Struktur
        distances = data_array[np.argmin(energies), 1:]
        expanded_weights = None  # Keine Gewichte nötig
    else:
        # Alle Abstände
        distances = data_array[:, 1:].flatten()
        # Erstelle ein erweitertes Gewicht-Array für alle Abstände
        expanded_weights = np.repeat(weights, data_array.shape[1] - 1)
    print(distances)
    print(expanded_weights)
    # Erstelle das Histogramm
    plt.figure(figsize=(10, 6))
    
    if mode == 'stable':
        # Einfaches ungewichtetes Histogramm für die stabilste Struktur
        hist, bins, _ = plt.hist(distances, bins=args.bins, 
                                alpha=0.7, color='blue', edgecolor='black')
        
        # Statistik für die stabilste Struktur
        mean = np.mean(distances)
        std = np.std(distances)
        min_energy = energies[np.argmin(energies)]
        
        stats_text = f'average: {mean:.3f} Å\n'
        stats_text += f'standard deviation: {std:.3f} Å\n'
        stats_text += f'reference energy: {min_energy:.6f} Eh'
    else:
        # Gewichtetes Histogramm
        hist, bins, _ = plt.hist(distances, bins=args.bins, weights=expanded_weights, 
                                alpha=0.7, color='blue', edgecolor='black')
        
        # Berechne gewichtete Statistiken
        weighted_mean = np.average(distances, weights=expanded_weights)
        weighted_std = np.sqrt(np.average((distances - weighted_mean)**2, weights=expanded_weights))
        
        stats_text = f'weighted average: {weighted_mean:.3f} Å\n'
        stats_text += f'weighted deviation: {weighted_std:.3f} Å\n'
        stats_text += f'number of structures: {len(energies)}'
    
    plt.title(title)
    plt.xlabel('distance (Å)')
    plt.ylabel('count')
    plt.grid(alpha=0.3)
    
    # Statistik-Text hinzufügen
    plt.annotate(stats_text, xy=(0.97, 0.97), xycoords='axes fraction', 
                 horizontalalignment='right', verticalalignment='top',
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
    
    # Speichere das Histogramm
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Histogramm wurde als '{args.output}' gespeichert.")
    
    # Histogramm anzeigen, wenn gewünscht
    if args.show:
        plt.show()
    else:
        plt.close()

if __name__ == "__main__":
    main()

