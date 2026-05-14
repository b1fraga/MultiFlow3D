F90=mpif90

OPTIONS = -cpp -O2 -g -fopenmp -DUSE_JSON=1 -DUSE_HDF5=1

ifneq (,$(findstring GNU,$(shell $(F90) --version)))
OPTIONS += -std=f2018
endif

LOPTIONS = -O2 -fopenmp

JSON_FORTRAN_INCLUDE_PATH := build_json_fortran/include/
JSON_FORTRAN_LIBRARY_PATH := build_json_fortran/lib/

TEST_DRIVE_INCLUDE_PATH := build_test_drive/include/
TEST_DRIVE_LIBRARY_PATH := build_test_drive/

HDF5_INCLUDE_PATH := build_hdf5/mod/shared/
HDF5_LIBRARY_PATH := build_hdf5/bin


objects = \
src/module_mpi.o \
src/hdf5_io.o \
src/json_io.o \
src/module_multidata.o \
src/module_vars.o \
src/module_vars_pt.o \
src/io.o \
src/module_LSM.o \
src/multiflow3d_sem.o \
src/imb.o \
src/shapes.o \
app/fdstag.o \
src/initial.o \
src/init_particle.o \
src/localparameters.o \
src/alloc_dom.o \
src/post.o \
src/flosol.o \
src/checkdt.o \
src/bounds.o \
src/bounds_keps.o \
src/sipsol.o \
src/convection.o \
src/diffusion.o \
src/newsolv_mg.o \
src/mgsolver.o \
src/wall_function.o \
src/log_law.o \
src/alloc_pt.o \
src/MPI_pt.o \
src/delta_func.o \
src/collision.o \
src/LPT.o \
src/timesig.o \
src/weno.o \
src/energy.o \
src/press.o \
src/roughness_function.o \
src/rungek.o \
src/averaging.o \
src/eddyvis_smag.o \
src/eddyvis_wale.o \
src/eddyvis_1eqn.o \
src/eddyvis_keps.o \
src/exchange_bc.o \
src/exchangep.o \
src/exchangepp.o \
src/exchangesca.o \
src/exchange.o \
src/exchangeu.o \
src/exchangev.o \
src/exchangew.o \
src/exchange_phi.o \
src/bounds_lsm.o \
src/lsm.o \
src/sediment.o


test_objects = \
src/hdf5_io.o \
src/json_io.o \
src/io.o \
tests/test_hdf5_io.o \
tests/test_json_io.o \
tests/test_io.o \
tests/main.o


all: test


test: M3D_v2.exe tests/tests.exe
	@cd tests && ./tests.exe && cd ..


M3D_v2.exe: $(objects)
	$(F90) $(objects) $(LOPTIONS) \
	-I./build_json_fortran/include/ \
	-L./build_json_fortran/lib/ \
	-Wl,-rpath,$(CURDIR)/build_json_fortran/lib/ \
	-Wl,-rpath,$(CURDIR)/build_hdf5/bin \
	-I./build_hdf5/mod/shared/ \
	-L./build_hdf5/bin \
	-ljsonfortran -lhdf5 -lhdf5_fortran \
	-o M3D_v2.exe


src/%.o: src/%.f90
	$(F90) $(OPTIONS) -c $< -o $@ \
	-I./build_json_fortran/include/ \
	-I./build_hdf5/mod/shared/


app/%.o: app/%.f90
	$(F90) $(OPTIONS) -c $< -o $@ \
	-I./src \
	-I./build_json_fortran/include/ \
	-I./build_hdf5/mod/shared/


tests/%.o: tests/%.f90
	$(F90) $(OPTIONS) -c $< -o $@ \
	-I./src \
	-I./build_json_fortran/include/ \
	-I./build_test_drive/include/ \
	-I./build_hdf5/mod/shared/


tests/tests.exe: $(test_objects)
	$(F90) $(test_objects) $(LOPTIONS) \
	-I./build_json_fortran/include/ \
	-L./build_json_fortran/lib/ \
	-I./build_test_drive/include/ \
	-L./build_test_drive/ \
	-I./build_hdf5/mod/shared/ \
	-L./build_hdf5/bin \
	-Wl,-rpath,$(CURDIR)/build_json_fortran/lib/ \
	-Wl,-rpath,$(CURDIR)/build_test_drive/ \
	-Wl,-rpath,$(CURDIR)/build_hdf5/bin \
	-ljsonfortran -ltest-drive -lhdf5 -lhdf5_fortran \
	-o tests/tests.exe


clean:
	rm -rf src/*.o src/*.mod
	rm -rf app/*.o app/*.mod
	rm -rf tests/*.o tests/tests.exe
	rm -f M3D_v2.exe
