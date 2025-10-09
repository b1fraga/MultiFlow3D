# Build Instructions

## Prerequisite

To build this code and its dependencies you'll be expected to have fpm installed
and mpif90 built.

## Building MultiFlow3D

First you need to git clone the repository and its submodule using
(replace branch_name with the name of the branch your interesting in building)

```
git clone --recurse-submodules git@github.com:b1fraga/MultiFlow3D.git
cd ./MultiFlow3D
git checkout branch_name
```

To build MultiFlow3D you need json-fortran for input. This is built by executing
the following once the repo has been executed

```
cd ./json-fortran
export FPM_FC="mpif90"
fpm build
cd ..
```

Now you are in a position to build MultiFlow3D by executing

```
make -j $(nproc --all)
```
