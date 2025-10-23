# Build Instructions

## Prerequisite

To build this code and its dependencies you'll be expected to already have mpif90
and cmake already installed.

## Building MultiFlow3D

First you need to git clone the repository and its submodule using
(replace branch_name with the name of the branch your interesting in building)

```
git clone --recurse-submodules git@github.com:b1fraga/MultiFlow3D.git
cd ./MultiFlow3D
git checkout branch_name
```

To build MultiFlow3D you need json-fortran for input, and test-drive for testing/
These can be built by executing the following once the repo has been executed

```
cmake -S ./json-fortran/ -B ./build_json_fortran/ -DCMAKE_Fortran_COMPILER=mpif90
cmake --build build_json_fortran
cmake -S ./test-drive/ -B ./build_test_drive/ -DCMAKE_Fortran_COMPILER=mpif90
cmake --build build_test_drive
```

Now you are in a position to build MultiFlow3D by executing

```
make -j $(nproc --all)
```
