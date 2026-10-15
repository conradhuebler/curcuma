# Usage notes for individual tools

Per-tool usage text moved unchanged from the README on 2026-10-01 (RMSD, docking, ConfScan, trajectories,
optimisation, MD, ConfSearch, BMT output). The JSON blocks are parameter dumps from 2025 and some defaults
no longer match (for example the default `method` is `gfnff` now, not `uff`). The authoritative list for
the installed binary is `curcuma -list-modules` and `curcuma -export-config <module>`; descriptions come
from `curcuma -help [category]`.

## General
curcuma catches the Ctrl-C signals from the console if used on Linux platform. It will then create an empty file called "stop". Some methods, like **confscan**, regularly check for that file and finalise the current task. If Ctrl-C is signaled and a "stop" file already exists, curcuma will stop immediately.

## RMSD Calculator 
```sh
curcuma -rmsd file1.xyz file2.xyz
```
Computes RMSD. If the two structures are ordered differently, curcuma will automatically reorder the atom list. To force reordering use
```sh
curcuma -rmsd file1.xyz file2.xyz -reorder
```

Two basic approaches are currently implemented, one incremental (testing many, but not all possible orders, parallelised) method without any requirements and another Kuhn-Munkres like approach. 
The Kuhn-Munkres like approach needs a prior alignment (only cude) of the structures, which may be obtained using a smaller substructure (template). One way is to use a molecule in a supramolecular structure . Use
```sh
-method template -fragment 1
```
to use the second structure, eg the one bound non-covalently by the first structure, as template. With this approach a much faster reordering is obtained. Omitting -fragment, curcuma tries using the smallest fragment.

Alternatively, generate templates using the incremental approach. For this, all atoms of one element (nitrogen is the default) are used. After the templates are generated, the structures oriented according the templates are reordered using the Kuhn-Munkres like approach. Use
```sh
-method hybrid -element 8
```
for taking oxygen as template.
 
Add
```sh
-heavy
```
to perform calculation only on non-proton atoms.

A new, and all previous published (and not yet published) methods was propsed Feb 2023 by Vásquez-Pérez and coworkers. It outperforms all methods so far (and sometimes the methods implemented in curcuma). 
It can be obtained at [Github](https://github.com/qcuaeh/molalignlib) and included in curcuma for RMSD calculation conformational filtering.
```sh
-reorder -method molalign -molalignbin /anypath/molalign
```

If the method was used, please cite the authors!
[J. Chem. Inf. Model. 2023, 63, 4, 1157–1165](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.2c01187)
The method can be applied during confscan, however, problems with the random numbers occur during the test runs for larger problems, making the less ideal than the natively in curcuma implemented methods.
```sh
-rmsdmethod molalign -molalignbin /anypath/molalign
```

```json
{ "reorder", false },
{ "check", false },
{ "heavy", false },
{ "fragment", -1 },
{ "fragment_reference", -1 },
{ "fragment_target", -1 },
{ "init", -1 },
{ "pt", 0 },
{ "silent", false },
{ "storage", 1.0 },
{ "method", "incr" },
{ "noreorder", false },
{ "threads", 1 },
{ "Element", 7 },
{ "DynamicCenter", false },
{ "order", "" },
{ "check", false },
{ "topo", 0 },
{ "write", 0 },
{ "moi", false },
{ "update-rotation", false },
{ "damping", 0.8 },
{ "split", false },
{ "nomunkres", false },
{ "dmix", -1 },
{ "molalignbin", "molalign" }
```


## Docking tool
Some docking can be performed (WIP).

Use
```sh
curcuma -dock -host A.xyz -guest B.xyz
```
to perform docking of B as guest and A as host molecule, use
```sh
curcuma -dock -complex AB.xyz
```
to perform docking on a complex or use
```sh
curcuma -dock -complex AB.xyz -guest C.xyz
```
to replace substrat structures in complex AB with C.


Use
```sh
curcuma -dock -host A.xyz -guest B.xyz -Step_X X  -Step_Y Y -Step_Z Z
```
with X, Y and Z being the steps of rotation. With X = 10, 10 rotations with 360/X ° will be performed.

Use
```sh
curcuma -dock -host A.xyz -guest B.xyz -Pos_X X  -Pos_X Y -Pos_X Z
```
with {X, Y, Z} being the initial anchor position for the substrat.

After docking a PseudoFF optimisation of the docking position will be performed, where the Lennard-Jones-Potential between both structures is calculated. XTB is then used to preoptimise the unique docking structures and the results are filtered using ConfScan and the template based reordering approach.

```json
{ "Pos_X", 0.0 },
{ "Pos_Y", 0.0 },
{ "Pos_Z", 0.0 },
{ "AutoPos", true },
{ "Filter", true },
{ "PostOpt", true },
{ "Step_X", 10 },
{ "Step_Y", 10 },
{ "Step_z", 10 },
{ "Host", "none" },
{ "Guest", "none" },
{ "Complex", "none" },
{ "scaling", 1.5 },
{ "NoOpt", false },
{ "CentroidMaxDistance", 1e5 },
{ "CentroidTolDis", 1e-1 },
{ "RotationTolDis", 1e-1 },
{ "Threads", 1 },
{ "DockingThreads", 1 },
{ "Charge", 0 },
{ "Cycles", 1 },
{ "RMSDMethod", "incr" },
{ "RMSDThreads", 1 },
{ "RMSDElement", 7 }
```



## Conformation Filter
Curcuma has some conformation filter based on energy, rmsd, rotation constants and rank limitation. As some structures may be identic, yet can not be aligned due to atom ordering, depending on the difference of the energy and rotational constants and the rmsd, automatic reordering in rmsd calculation will be performed.

Use
```sh
curcuma -confscan conformation.xyz
```
to simple filter the conformation. The results will be stored in an additional xyz file.

Use
```sh
curcuma -confscan conformation.xyz -MaxHTopoDiff 0
```
to ensure, that even the rmsd is smaller than the threshold, the second molecule is only rejected if there is no difference in the hydrogen bond pattern. Set to ***-MaxHTopoDiff 2*** if two two changes are allowed. A detailed description will follow some time.
```sh
curcuma -confscan conformation.xyz -reorder
```
to force reordering for every rmsd calculation.

Add
```sh
-heavy
```
to perform rmsd calculation and reordering only on non-proton atoms. Reordering is much faster then!

Adding
```sh
-RMSDmethod template
```
the template based reordering is used, **hybrid** is also working.

With
```sh
-Useorders X
```
up to X reorder results will be reused from the last reorder calculation. Template based approaches result in only one reorder rule, however X is set to 0 in automatically, but can be changed with that argument.

By adding the argument
```sh
-accepted accepted.xyz
```
a file with already accepted structures can be passed to curcuma. Molecules in that file (accepted.xyz) will be rejected if they appear in the conformation.xyz, thus several files with conformation can be joined.

Confscan will write a restart file, finalise and quit, if a file called "stop" is found in the working directory. Such file will be generated if Ctrl-C is hit or if it is created using for example the **touch stop** command.
Within a restart file, the last energy difference and the atom indicies from reordering are stored. A restart file will automatically be read upon the start of curcuma. The content of the restart file will be used to speed up the 2nd step of the conformation filtering procedure.

Confscan supports the molalign tool. However, as too often reordering with molalign is not working, it can efficiently be used if the RMSD is only slightly above the threshold. 
```sh
-domolalign 1.1
```
Sets the threshold to 1.1*RMSDthreshold. If the molecule was accepted as to different, but the RMSD is blow 1.1*RMSDthreshold molalign will check too.

Confscan write a statistic file, where for each rejected molecule the reference alongside the energy difference and the RMSD is printed out. Furthermore, the reordered indices are given, if available. Molalign does not return the reordered indices, hence they are empty or marked **0,0** if the reordered was finally performed using molalign in a standard run.

```json
{ "noname", true },
{ "restart", true },
{ "heavy", false },
{ "rmsd", -1 },
{ "rank", -1 },
{ "writeXYZ", false },
{ "forceReorder", false },
{ "check", false },
{ "energy", 1.0 },
{ "maxenergy", -1.0 },
{ "preventreorder", false },
{ "scaleLoose", 1.5 },
{ "scaleTight", 0.1 },
{ "scaleLooseEnergy", 1.2 },
{ "scaleTightEnergy", 0.1 },
{ "scaleLooseRotational", 1.2 },
{ "scaleTightRotational", 0.1 },
{ "scaleLooseRipser", 1.2 },
{ "scaleTightRipser", 0.1 },
{ "skip", 0 },
{ "allxyz", false },
{ "update", false },
{ "MaxParam", -1 },
{ "UseOrders", -1 },
{ "RMSDMethod", "hybrid" },
{ "MaxHTopoDiff", -1 },
{ "threads", 1 },
{ "RMSDElement", 7 },
{ "accepted", "" },
{ "method", "" },
{ "lastdE", -1 },
{ "fewerFile", false },
{ "dothird", true },
{ "skipfirst", false },
{ "ignoreRotation", false },
{ "ignoreBarCode", false },
{ "skipless", false },
{ "looseThresh", 7 },
{ "tightThresh", 3 },
{ "update-rotation", false },
{ "damping", 0.8 },
{ "split", false },
{ "writefiles", false },
{ "nomunkres", false },
{ "molalignbin", "molalign" },
{ "ripser_xmax", 4 },
{ "ripser_xmin", 0 },
{ "ripser_ymax", 4 },
{ "ripser_ymin", 0 },
{ "ripser_bins", 10 },
{ "ripser_scaling", 0.1 },
{ "ripser_stdx", 10 },
{ "ripser_stdy", 10 },
{ "ripser_ratio", 1 },
{ "ripser_dimension", 2 },
{ "domolalign", -1 }
```

```cpp
/* rotational = 1
 * ripser     = 2
 * energy     = 4 */
int looseThresh = 1 * (diff_rot < m_diff_rot_threshold_loose) + 2 * (diff < m_diff_ripser_threshold_loose) + 4 * (std::abs(mol1->Energy() - mol2->Energy()) * 2625.5 < m_diff_energy_threshold_loose);
if ((looseThresh & m_looseThresh) == m_looseThresh) 
{

}
```
## Find unique structures in trajectories
xyz and trj are handled equally.
```sh
curcuma -rmsdtraj XXX.trj -writeUnique -rmsd 1.5
```


```json
{ "writeUnique", false },
{ "writeAligned", false },
{ "rmsd", 1.5 },
{ "fragment", -1 },
{ "reference", "none" },
{ "second", "none" },
{ "heavy", false },
{ "pcafile", false },
{ "allxyz", false },
{ "RefFirst", false },
{ "noreorder", true },
{ "opt", false },
{ "filter", false },
{ "writeRMSD", true },
{ "offset", 0 }
```

## Geometry optimisation (batch mode possible)
Geometry optimisation can be performed with curcuma using 
```sh
curcuma -opt XXX.xyz
```
A file called XXX.opt.xyz with the optimised structures will be written. The individual steps are stored in XXX.trj.xyz. The number of threads can be controlled with
```sh
-threads X
```

```json
{ "writeXYZ", true },
{ "printOutput", true },
{ "dE", 0.1 },
{ "dRMSD", 0.01 },
{ "method", "uff" },
{ "MaxIter", 5000 },
{ "LBFGS_eps", 1e-5 },
{ "StoreIntermediate", 2000 },
{ "SingleStep", 20 },
{ "ConvCount", 11 },
{ "GradNorm", 0.001 },
{ "Threads", 1 },
{ "Charge", 0 },
{ "Spin", 0 },
{ "SinglePoint", false },
{ "optH", false },
{ "serial", false }
```


```cpp
/*
 * Energy = 1
 * RMSD = 2
 * LBFGS Conv = 4
 * Gradient Norm = 8
 * */
converged = 1 * (abs(fun.m_energy - final_energy) * 2625.5 < dE)
    + 2 * (driver->RMSD() < dRMSD)
    + 4 * (solver.isConverged())
    + 8 * (solver.final_grad_norm() < GradNorm);
perform_optimisation = (converged != ConvCount) && (fun.isError() == 0);
```

## Reorder and Align trajectories
To reorder trajectory files with dissordered atomic indicies, for example after merging several minimum energy path files from NEB calculation, use
```sh
curcuma -rmsdtraj XXX.xyz -writeAligned
```
Reordering will be done with respect to the previouse structure in the trajectory. If the first structure should be used, add ***-reffirst*** as additional argument. The new trajectory is called XXX_aligned.xyz.

Using ***-rmsdtraj*** argument, a file **XXX_rmsd.dat** will be written, where the rmsd is stored.

## Distance and angle calculation

```sh
curcuma -distance XXX.trj atom1 atom2
```

```sh
curcuma -angle XXX.trj atom1 atom2 atom3
```

The index starts with 1. Using grep and sed via ***|grep '::' |sed 's/:://g'*** omitts unused output.

## Compare two RDG vs sign(λ<sub>2</sub>)ρ plots 
Using 
```sh
curcuma -nci file1.dat file2.dat
```
one can ''remove'' RDG vs sign(λ<sub>2</sub>)ρ points which occur in both plots (file1.dat and file2.dat). The similarity of two points is set to true, if the distance is below a threshold distance, which is defined by the averaged distance of two adjacent points.

## Molecular Dynamics and Metadynamics
Curcuma has now a Molecular Dynamics modul, which can be used with:
```sh
curcuma -md input.xyz
```

### Possible options
```json
{ "writeXYZ", true },
{ "printOutput", true },
{ "MaxTime", 5000 },
{ "T", 298.15 },
{ "dt", 1 }, // single step in fs
{ "rm_COM", 100 }, // remove translation and rotation every x fs
{ "charge", 0 },
{ "Spin", 0 },
{ "rmrottrans", 0 },
{ "nocenter", false },
{ "dump", 50 },
{ "print", 1000 },
{ "unique", false },
{ "rmsd", 1.5 },
{ "opt", false },
{ "hmass", 1 },
{ "velo", 1 },
{ "rescue", false },
{ "coupling", 10 },
{ "MaxTopoDiff", 15 },
{ "impuls", 0 },
{ "method", "uff" },
{ "impuls_scaling", 0.75 },
{ "writeinit", false },
{ "initfile", "none" },
{ "norestart", false },
{ "writerestart", 1000 },
{ "rattle", false },
{ "rattle_tolerance", 1e-6 },
{ "rattle_maxiter", 10 },
{ "thermostat", "csvr" },
{ "respa", 1 },
{ "dipole", false },
{ "seed", 1 },
{ "cleanenergy", false },
{ "wall", "none" }, // can be spheric or rect
{ "wall_type", "logfermi" }, // can be logfermi or harmonic
{ "wall_spheric_radius", 0 },
{ "wall_xl", 0 },
{ "wall_yl", 0 },
{ "wall_zl", 0 },
{ "wall_x_min", 0 },
{ "wall_x_max", 0 },
{ "wall_y_min", 0 },
{ "wall_y_max", 0 },
{ "wall_z_min", 0 },
{ "wall_z_max", 0 },
{ "wall_temp", 298.15 },
{ "wall_beta", 6 },
{ "mtd", false },
{ "plumed", "plumed.dat" }
```

For example, using 
```sh
curcuma -md input.xyz -method gfnff
``` 
the GFN-FF approach will be used.

```sh
curcuma -md input.xyz -method gfnff -T 500  -berendson 200 -dt 1 -hmass 1 -thermostat_steps 400 -velo 4 -maxtime 2e4 -dt 0.5 -impuls 500 -impuls_scaling 0.75
``` 
will perform some kind of conformational search using GFN-FF. Results are stored in **input.unique.xyz**! Repeating it will result in other conformations and the previous results stored in **input.unique.xyz** will be overwritten. Bonds may break from time to time ...

Rattle can be used to constrain (currently) all bonds, allowing larger time steps for integration. Up to 8 fs might be possible.
```sh
curcuma -md input.xyz -rattle -dt 4
``` 

The MD implementation integrates well into curcuma, hence calculation can be stopped with Ctrl-C (or a "stop" file) and will be resumed (velocities and geometries are stored) if a restart file is found.

The thermostat target temperature can follow a multi-stage **ramp** and individual atom subsets can be thermostatted as separate **regions**:
```sh
curcuma -md input.xyz -method gfnff -temperature 300 -temp_ramp true -temp_schedule "600:steps:5000;300:reach:10"
```
See [docs/TEMPERATURE_RAMP.md](TEMPERATURE_RAMP.md) for the schedule grammar (`steps`/`reach`), the `temp_regions` JSON array, live temperature control, and per-thermostat support.

On a large solvated system the default 1 fs step can be too long for a momentarily compressed X-H
bond, and a single such step heats the whole trajectory. Instead of lowering `-dt` or raising the
hydrogen mass for the whole run, the integrator can redo just that step with a subdivided one:
```sh
curcuma -md input.xyz -method gfnff -dt 1.0 -adaptive_step true
```
It is off by default and adds no constraint and no mass modification. Two things are watched: the
total energy of the step, and - because that one loses its contrast as the system grows, while a
violating step stays on a handful of atoms - the kinetic energy of the **hottest atom relative to
the per-atom mean**. On a 7320-atom solvated polymer over 300 fs the second channel is what works:
+67.59 Eh and 2175 K become **+0.53 Eh and 246 K**, and 92 of the 104 rejected steps were ones the
energy criterion accepted. See [docs/MD_LARGE_SYSTEMS.md](MD_LARGE_SYSTEMS.md) for the
mechanism (one water molecule collapsing), the calibration tables and what it does not do.

With
```sh
curcuma -md input.xyz -mtd
``` 
a metadynamics simulation can be performed using plumed. It is a ***plumed.dat*** expected, or can be set with
```sh
curcuma -md input.xyz -mtd -plumed plumed.dat
```

See [docs/PLUMED_HELP.md](PLUMED_HELP.md) for the full PLUMED integration guide (unit conversions, output files, available CVs, thermal equilibration gate, internal RMSD-MTD).

## Conformational Search (dual-method)

The MD-driven conformational search (`-confsearch`) can explore with a cheap method and refine/rank the discovered conformers with a more accurate one:

```sh
curcuma -confsearch input.xyz -md_method gfnff -opt_method gfn2
```

`-md_method` runs the MD exploration and the pre-optimization; `-opt_method` runs the per-cycle accurate re-optimization and the final ranking. Both fall back to `-method` when unset, so `curcuma -confsearch input.xyz -method gfnff` keeps the single-method behaviour. See [docs/CONFSEARCH_DUAL_METHOD.md](CONFSEARCH_DUAL_METHOD.md).

The search is **restartable** with `-restart`: a self-contained checkpoint (bias pool, cumulative conformers, seeds, energies, schedule position) is written after every MD phase and every temperature cycle, into the BMT dir and copied back to the start directory. Re-running the same command with `-restart` resumes from it (kill the process to interrupt; the checkpoint persists). See [docs/CONFSEARCH_RESTART.md](CONFSEARCH_RESTART.md).

## Output Directory System (BMT)

By default, all curcuma commands create a **Basename.Method.Timestamp** directory for their output files. For example:

```sh
curcuma -md water.xyz -method gfnff
# Output goes to: water.md.20260609_143052/
```

The BMT directory contains all trajectory files, restart data, and a `metadata.json` file with calculation details (basename, method, timestamp, input file). This keeps the working directory clean and makes it easy to compare runs.

To copy specific files back to the working directory after the calculation finishes, use the `-bak` flag:

```sh
curcuma -opt water.xyz -method gfnff -bak water.opt.xyz
# water.opt.xyz is copied from the BMT directory to the working directory
```

Multiple files can be specified: `-bak water.opt.xyz -bak water.trj.xyz`.

To disable BMT and write output to the working directory (legacy behavior):

```sh
curcuma -md input.xyz -method uff -no_bmt
```
