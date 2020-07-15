# Curcuma

A simple Open Source molecular modelling tool.

## Download and requirements
git clones automatically eigen.
- [eigen](https://gitlab.com/libeigen/eigen) provides eigen C++ library for linear algebra
- [LBFGSpp](https://github.com/yixuan/LBFGSpp/) prvodies LBFGS optimiser
- [XTB](https://github.com/grimme-lab/xtb) eXtended TightBinding - Some methods use this, however it is disabled by default. Add '-DUSE_XTB=true ' to the cmake command line to enable it. GCC 8 or later has to be used.

### Using XTB in curcuma
XTB is automatically obtained during git clone, however it is not included in curcuma with default compile settings. There are two ways of including XTB:
- Compiling and linking automatically, set **COMPILE_XTB** to true ( with -DCOMPILE_XTB=true or via ccmake). XTB will be compiled with the c++ and fortran compiler in the path (using the blas and lapack libraries in the path as well). The xtb program is then quite slow.
- Linking curcuma to the official library (get if from the xtb github page). Set **LINK_XTB** to true and define the path to **XTB_DIR** where the libxtb.so has been placed. As the xtb source has been obtained already, no additional header file is needed. However, curcuma has to be compiled with an intel compiler.

- Only xtb GFN2 can be used as method right now.
- As xtb is not thread-safe (yet), xtb itself can be run in parallel, but only one xtb process can be executed. This affects the the crude optimisation after docking.

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

## RMSD Calculator
```sh
curcuma -rmsd file1.xyz file2.xyz
```
Computes RMSD. If the two structures are ordered differently, curcuma will automatically reorder the atom list. To force reordering use
```sh
curcuma -rmsd file1.xyz file2.xyz -reorder
```

Two different approaches are currently implemented, one incremental method without any requirements. The template based approaches is useful for supramolecular systems, where the substrate (fragment 1, counting begins with 0) structures are ordered equally. Use
```sh
-method template -fragment 1
```
to use the second structure, eg the one bound non-covalently by the first structure, as template. With this approach a much faster reordering is obtained. Omitting -fragment, curcuma tries used the smallest fragment.

Add
```sh
-heavy
```
to perform calculation only on non-proton atoms.

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

## Conformation Filter
Curcuma has some conformation filter based on energy, rmsd, rotation constants and rank limitation. As some structures may be identic, yet can not be aligned due to atom ordering, depending on the difference of the energy and rotational constants and the rmsd, automatic reordering in rmsd calculation will be performed.

Use
```sh
curcuma -confscan conformation.xyz
```
to simple filter the conformation. The results will be printed out.

Use
```sh
curcuma -confscan conformation.xyz -writeXYZ
```
to store the structures as xyz files

Use
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
the template based reordering is used.

With
```sh
-Useorders X
```
up to X reorder results will be reused from the last reorder calculation. Template based approaches result in only one reorder rule, however X is set to 0 in automatically, but can be changed with that argument.
## NEB Input structure preparation
Some very experimental stuff - reorder (like in -rmsd a.xyz b.xyz -reorder) and rmsd alignment and fragment moving to prepare some NEB structures.

```sh
curcuma -nebprep start.xyz end.xyz
```

## Find unique structures in trajectories
xyz and trj are handled equally.
```sh
curcuma -rmsdtraj XXX.trj -write -rmsd 1.5
```


Have a lot of fun!
