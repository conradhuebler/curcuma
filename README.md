# Curcuma

A simple Open Source molecular modelling tool.

## Download and requirements
git clones automatically eigen.
- [eigen](https://gitlab.com/libeigen/eigen) provides eigen C++ library for linear algebra

## Compiling
To compile SupraFit you will need [CMake](https://cmake.org/download/) 3 or newer and a C++14-capable compiler.

To obtain the most recent version
```sh
git clone --recursive https://github.com/conradhuebler/curcuma
```

Compile it as follows on Unix Platform:
```sh
cd suprafit
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```
or use the script in the subdirectory. It should automatically update the submodules.
```sh
sh scripts/build_unix.sh 
```
On Windows Systems use for example
```sh
cd suprafit
mkdir build
cd build
```
For Visual Studio use
```sh
cmake -G "Visual Studio 15 2017 Win64" -DCMAKE_BUILD_TYPE=Release  ..
```

or for MinGW (openMP is enabled, libgomp-1.dll is expected) use

```sh
cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release ..
```

```sh
cmake --build . --config Release
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

Have a lot of fun!
