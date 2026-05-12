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

To build MultiFlow3D you need json-fortran for input, hdf5 for output
and test-drive for testing. These can be built by executing the following 
once the repo has been executed

```
cmake -S ./json-fortran/ -B ./build_json_fortran/ \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_INSTALL_PREFIX=$PWD/install

cmake --build build_json_fortran --parallel $(nproc --all)
cmake --install build_json_fortran


cmake -S ./test-drive/ -B ./build_test_drive/ \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_INSTALL_PREFIX=$PWD/install

cmake --build build_test_drive --parallel $(nproc --all)
cmake --install build_test_drive


cmake -S ./hdf5/ -B ./build_hdf5/ \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DHDF5_BUILD_FORTRAN=ON \
  -DCMAKE_INSTALL_PREFIX=$PWD/install

cmake --build build_hdf5 --parallel $(nproc --all)
cmake --install build_hdf5
```

Now you are in a position to build MultiFlow3D by executing

```
cmake -S . -B ./build/ -DCMAKE_Fortran_COMPILER=mpif90  -DCMAKE_PREFIX_PATH="$PWD/install/" -Dtest-drive_DIR="$PWD/install/lib/cmake/test-drive/" -DM3D_v2_USE_HDF5=ON -DM3D_v2_USE_JSON=ON -DM3D_v2_BUILD_TESTS=ON -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build --parallel $(nproc --all)
```

## Documentation

To build the documentation you must first install the tool ford. This can be installed by executing

```
python -m pip install ford
```

You then build the documentation by executing

```
ford  -d . README.md
```

Finally you open up the documentation by executing

```
open ./doc/index.html
```

which will open up the homepage of the documentation in your default browser.

