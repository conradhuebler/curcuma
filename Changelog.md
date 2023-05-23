# The (incomplete) curcuma Changelog

### pre Alpha

- molalign can be used for reordering
- add forked LBFGSpp for single steps in geometry optimisation
- add tblite and forked xtb for better control of xTB calculation
- extended centroids for docking
- compare/analyse NCIPLOT's RDG vs rho plots
- add parallel batch optimisation
- add conformation statistics
- parallel docking due to CxxThreadPool
- add restriction for hydrogen bond patterns in conformational scan
- add template based reorder method for rmsd calculation
- make confscan (silently) restartable
- docking with post-optimisation and filtering (needs XTB GFN 2)
- add reordering for non-conformer/non-isomer structures
- add step-wise-rmsd for two trajectories
- find more or less unique conformers in trajectories
- scan xyz files (trajectories) for hbonds and prints each single distance for each step
- include optional XTB
- prepare NEB input geometries (reordered and aligned)
- add -led option to prepare fragment assigned xyz files (for ORCA LED calculation)
- Add conformation filter (energy, rmsd, rank) - with automatic/forced reordering of structures
- Reorder atom indicies in molecule, keeping methyl connectivitiy
- Initial docking
- Calculate RMSD for two molecules
