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

Add
```sh
-heavy
```
to perform calculation only on non-proton atoms.

## Docking tool
Some docking can be performed (WIP). Rotation will be added soon. (And many more.)

```sh
curcuma -dock A.xyz B.xyz XXX YYY ZZZ
```
To move molcule B to the position.

Use
```sh
curcuma -dock A.xyz B.xyz
```
to automatically move B to the center of mass of A, which is only interesting for special macrocyclic structures.

After docking a PseudoFF optimisation of the docking position will be performed, where the Lennard-Jones-Potential between both structures is calculated.

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
