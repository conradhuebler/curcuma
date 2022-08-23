[![CodeFactor](https://www.codefactor.io/repository/github/conradhuebler/curcuma/badge)](https://www.codefactor.io/repository/github/conradhuebler/curcuma) [![Build](https://github.com/conradhuebler/curcuma/workflows/AutomaticBuild/badge.svg)](https://github.com/conradhuebler/curcuma/actions)  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4302722.svg)](https://doi.org/10.5281/zenodo.4302722)


# Curcuma

A simple Open Source molecular modelling tool.

## Download and requirements
git clones automatically some submodules.
- [LBFGSpp](https://github.com/conradhuebler/LBFGSpp) a fork of [yixuan/LBFGSpp](https://github.com/yixuan/LBFGSpp/) provides LBFGS optimiser, the fork allows performing single step optimisation without resetting any calculated optimsation history
- [XTB](https://github.com/conradhuebler/xtb) a fork of the [official xtb](https://github.com/grimme-lab/xtb) program for eXtended TightBinding - It is only used for GFN-FFF calculation. The adaptions in the fork suppress the output in the GFN-FF calculations. 
- [tblite](https://github.com/tblite/tblite) eXtended TightBinding via tblite for GFN1 and GFN2 calculation
- [CxxThreadPool](https://github.com/conradhuebler/CxxThreadPool) - C++ Thread Pool for parallel calculation
- [eigen](https://gitlab.com/libeigen/eigen) provides eigen C++ library for linear algebra. Eigen is not downloaded automatically, but will be fetched and updated if the build scripts in the **scripts** subdirectory are used.
- [fmt](https://github.com/fmtlib/fmt) formatted console output

Additionally, [nlohmann/json](https://github.com/nlohmann/json) is obtained via cmake.

### Using xTB in curcuma is now again WIP
xtb and tblite are automatically obtained during git clone, however it is not included in curcuma with the default compiler settings. Add '-DCOMPILE_XTB=true ' to the cmake command line to enable it. 
Using xtb calculation in curcuma can be controlled for now as follows:
- Add **-gfn 1** to run GFN1 calculation, **-gfn 66** to run GFN-FF calculation.
- GFN1 and GFN2 calculation are thread-safe now, so parallel optimisation are possible now. However, the variable **OMP_NUM_THREADS** should be set to 1. Add **-thread 12** to run optimisation after docking with 12 threads. GFN-FF calculation will fail in parallel mode.


**Please cite xtb if used within curcuma! The most recent information can be found [here](https://github.com/grimme-lab/xtb#citations)!**

## Compiling
To compile Curcuma you will need [CMake](https://cmake.org/download/) 3.15 or newer and a C++17-capable compiler, both gcc and icc (quite recent version) work.

To obtain the most recent version
```sh
git clone --recursive https://github.com/conradhuebler/curcuma
```

Compile it as follows on Unix Platform:
```sh
cd curcuma 
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

# Usage

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
{ "check", false }
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
{ "scaleLoose", 2 },
{ "scaleTight", 0.5 },
{ "skip", 0 },
{ "allxyz", false },
{ "update", false },
{ "MaxParam", -1 },
{ "UseOrders", -1 },
{ "RMSDMethod", "incr" },
{ "MaxHTopoDiff", -1 },
{ "GFN", -1 },
{ "RMSDThreads", 1 },
{ "RMSDElement", 7 },
{ "accepted", "" }
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
{ "noreorder", true }
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
{ "GFN", 2 },
{ "MaxIter", 2000 },
{ "LBFGS_eps", 1e-5 },
{ "StoreIntermediate", 600 },
{ "SingleStep", 20 },
{ "ConvCount", 4 },
{ "GradNorm", 0.001 },
{ "Threads", 1 },
{ "Charge", 0 },
{ "Spin", 0 },
{ "SinglePoint", false }
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

## Molecular Dynamics
Curcuma has now a Molecular Dynamics modul, which can be used with:
```sh
curcuma -md input.xyz
```

Naturally, there are many more better fitting tools to perform molecular dynamics. Unlike curcuma, they do respect correct units etc. However, MD calculations work and can be controlled with the following arguments:

```json
{ "GFN", 2 },
{ "MaxSteps", 5000 },
{ "T", 298.15 },
{ "dt", 1 },
{ "charge", 0 },
{ "Spin", 0 },
{ "centered", false }
```

For example, using 
```sh
curcuma -md input.xyz -gfn 66
``` 
the GFN-FF approach will be used.
For now, Velocity Rescaling will be applied to keep the Temperature at the desired level ...


The MD implementation integrates well into curcuma, hence calculation can be stopped with Ctrl-C (or a "stop" file) and will be resumed (velocities and geometries are stored) if a restart file is found.
# Citation
Please cite the software package if you obtain results:
[conradhuebler/curcuma: Curcuma Zenodo Citation](https://doi.org/10.5281/zenodo.4302722)

Have a lot of fun!
