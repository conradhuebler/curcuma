# dMatrix Migration Guide - Enhanced Topological Data Analysis

## Overview

The legacy `-dMatrix` command has been **completely integrated** into the modern unified analysis system. All functionality is preserved and enhanced while providing a cleaner, more configurable interface.

## Migration Path

### Quick Migration

**Old dMatrix usage:**
```bash
curcuma -dMatrix molecule.xyz
```

**New equivalent:**
```bash
curcuma -analysis molecule.xyz -topological.save_persistence_image true
```

### Complete Feature Mapping

| Legacy dMatrix Feature | New Analysis Equivalent |
|------------------------|-------------------------|
| Distance matrix export | `-topological.save_distance_matrix true` |
| Persistence pairs | `-topological.save_persistence_pairs true` |
| Persistence diagrams | `-topological.save_persistence_diagram true` |
| Persistence images | `-topological.save_persistence_image true` |
| Exclude hydrogen | `-topological.exclude_hydrogen true` |
| Exclude bonds | `-topological.exclude_bonds true` |
| Include elements | `-topological.print_elements true` |
| Include energy | `-topological.print_energy true` |
| Image formats | `-topological.image_format png|jpg|bmp|tga` |
| Colormaps | `-topological.colormap grayscale|jet|hot|viridis|coolwarm` |
| Image resolution | `-topological.resolution 800x800` |
| Post-processing | `-topological.post_processing none|adaptive|ring_focused` |

## Enhanced Features

The new system provides **additional capabilities** not available in the legacy dMatrix:

### Advanced Configuration
- **Multiple colormaps**: viridis, coolwarm support
- **Post-processing options**: adaptive, ring_focused algorithms
- **Flexible resolution**: Any WIDTHxHEIGHT format
- **Enhanced temperature control**: Fine-tuned image enhancement
- **Structure preservation**: Configurable balance between enhancement and original data

### Integration Benefits
- **Unified interface**: Works with all molecular formats (XYZ, VTF, MOL2, SDF, PDB)
- **Trajectory support**: Batch processing with statistical analysis
- **JSON configuration**: Consistent parameter handling
- **Comprehensive logging**: Detailed progress reporting and error handling
- **Modern architecture**: Thread-safe, memory-efficient implementation

## Complete Examples

### Basic Migration Examples

**Legacy dMatrix:**
```bash
curcuma -dMatrix structure.xyz
```

**New equivalent with same output:**
```bash
curcuma -analysis structure.xyz \
        -topological.save_distance_matrix true \
        -topological.save_persistence_pairs true \
        -topological.save_persistence_diagram true \
        -topological.save_persistence_image true
```

### Enhanced Workflows

**High-quality publication images:**
```bash
curcuma -analysis protein.pdb \
        -topological.save_persistence_image true \
        -topological.exclude_hydrogen true \
        -topological.colormap viridis \
        -topological.resolution 2048x2048 \
        -topological.image_format png \
        -topological.post_processing adaptive
```

**Trajectory analysis (new capability):**
```bash
curcuma -analysis dynamics.trj.xyz \
        -trajectory true \
        -topological.save_persistence_image true \
        -topological.temperature 3.0 \
        -output_file tda_results.json
```

**Machine learning feature generation:**
```bash
curcuma -analysis dataset.xyz \
        -topological.save_persistence_pairs true \
        -topological.exclude_hydrogen true \
        -output_format json \
        -output_file ml_features.json
```

## Output Files

The new system maintains the same file naming convention with enhanced indexing:

### Single Structure
- `prefix_0.dMat` - Distance matrix
- `prefix_0.pairs` - Persistence pairs
- `prefix_0.PD` - Persistence diagram (text)
- `prefix_0.PD.png` - Persistence diagram (image)
- `prefix_0.PI` - Persistence image (text)
- `prefix_0.PI.png` - Persistence image (visualization)

### Trajectory Analysis
- `prefix_N.dMat` - Distance matrix for frame N
- `prefix_N.pairs` - Persistence pairs for frame N
- `prefix_N.PD.png` - Persistence diagram for frame N
- `prefix_N.PI.png` - Persistence image for frame N
- Plus trajectory-wide statistical analysis in JSON output

## Technical Implementation

### Architecture
- **TDAEngine Class**: Encapsulates all dMatrix functionality
- **Modern Integration**: Seamless integration with unified analysis system
- **Research Quality**: Maintains publication-grade analysis capabilities
- **Performance**: Optimized memory usage and computational efficiency

### Backward Compatibility
- **File Formats**: All original output formats preserved
- **Algorithm Accuracy**: Identical computational results to legacy dMatrix
- **Citation Requirements**: Same research citation as original implementation

## Research Applications

The enhanced TDA system is suitable for:

- **Molecular Shape Analysis**: Quantitative comparison of molecular geometries
- **Protein Folding Studies**: Conformational analysis and folding pathway characterization
- **Materials Science**: Crystal structure analysis and phase identification
- **Drug Discovery**: Molecular recognition and binding site analysis
- **Machine Learning**: Feature generation for chemical machine learning models

## Citation

When using the enhanced TDA features, please cite:

> Townsend, J., Micucci, C.P., Hymel, J.H. et al. Representation of molecular structures with persistent homology for machine learning applications in chemistry. *Nat Commun* **11**, 3230 (2020). https://doi.org/10.1038/s41467-020-17035-5

## Migration Support

### Help System
```bash
# General analysis help
curcuma -analysis -help

# Comprehensive TDA documentation
curcuma -analysis -help-tda
```

### Configuration Files
The new system supports JSON configuration files for complex workflows:

```json
{
  "analysis": {
    "properties": "topology",
    "topological": {
      "save_persistence_image": true,
      "exclude_hydrogen": true,
      "colormap": "viridis",
      "resolution": "1024x1024",
      "image_format": "png",
      "post_processing": "adaptive",
      "temperature": 2.5
    }
  }
}
```

### Validation
To verify your migration is successful, compare output files from the legacy and new systems:
- Distance matrices should be numerically identical
- Persistence data should match within numerical precision
- Image outputs should be visually equivalent (with optional enhancements)

## Conclusion

The enhanced TDA system provides a **complete replacement** for the legacy `-dMatrix` functionality with significant improvements in usability, performance, and features. All research workflows can be migrated with minimal changes while gaining access to powerful new capabilities for molecular analysis.