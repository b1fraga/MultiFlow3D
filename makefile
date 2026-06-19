F90=mpif90

OPTIONS = -cpp -O2 -g -fopenmp -DUSE_JSON=1 -DUSE_HDF5=1 --coverage

ifneq (,$(findstring GNU,$(shell $(F90) --version)))
OPTIONS += -std=f2018
endif

LOPTIONS = -O2 -fopenmp --coverage

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
src/module_multidata.o \
src/module_mpi.o \
src/module_vars.o \
src/MPI_pt.o \
src/hdf5_io.o \
src/json_io.o \
src/io.o \
tests/test_hdf5_io.o \
tests/test_json_io.o \
tests/test_io.o \
tests/test_MPI_pt.o \
tests/main.o


all: test


test: M3D_v2.exe tests/tests.exe
	@cd tests && mpirun -np 4 ./tests.exe && cd ..


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

src/LPT.o : src/LPT.f90 src/module_vars_pt.o src/module_vars.o src/module_mpi.o src/module_multidata.o 
src/MPI_pt.o : src/MPI_pt.f90 src/module_mpi.o
src/alloc_dom.o : src/alloc_dom.f90 src/module_mpi.o src/module_vars.o src/module_multidata.o 
src/alloc_pt.o : src/alloc_pt.f90 src/module_vars_pt.o src/module_multidata.o src/module_mpi.o src/module_vars.o 
src/averaging.o : src/averaging.f90 src/module_LSM.o src/module_vars.o src/module_multidata.o 
src/bounds.o : src/bounds.f90 src/module_LSM.o src/imb.o src/module_mpi.o src/module_multidata.o src/module_vars.o 
src/bounds_keps.o : src/bounds_keps.f90 src/module_LSM.o src/module_multidata.o src/module_vars.o 
src/bounds_lsm.o : src/bounds_lsm.f90 src/module_multidata.o src/module_LSM.o src/module_vars.o 
src/checkdt.o : src/checkdt.f90 src/module_LSM.o src/module_multidata.o src/module_mpi.o src/module_vars.o 
src/collision.o : src/collision.f90 src/module_vars_pt.o src/module_vars.o src/module_mpi.o src/module_multidata.o 
src/convection.o : src/convection.f90 src/module_multidata.o src/module_vars.o 
src/delta_func.o : src/delta_func.f90 
src/diffusion.o : src/diffusion.f90 src/module_multidata.o src/module_mpi.o src/module_vars.o 
src/eddyvis_1eqn.o : src/eddyvis_1eqn.f90 src/module_multidata.o src/module_vars.o 
src/eddyvis_keps.o : src/eddyvis_keps.f90 src/module_multidata.o src/module_vars.o 
src/eddyvis_smag.o : src/eddyvis_smag.f90 src/module_multidata.o src/module_vars.o 
src/eddyvis_wale.o : src/eddyvis_wale.f90 src/module_multidata.o src/module_vars.o 
src/energy.o : src/energy.f90 src/module_mpi.o src/module_multidata.o src/module_vars.o 
src/exchange.o : src/exchange.f90 src/module_mpi.o src/module_multidata.o src/module_vars.o 
src/exchange_bc.o : src/exchange_bc.f90 src/module_vars.o src/module_mpi.o src/module_multidata.o 
src/exchange_bcphi.o : src/exchange_bcphi.f90 src/module_vars.o src/module_mpi.o src/module_multidata.o 
src/exchange_phi.o : src/exchange_phi.f90 src/module_vars.o src/module_mpi.o src/module_multidata.o 
src/exchangep.o : src/exchangep.f90 src/module_vars.o src/module_mpi.o src/module_multidata.o 
src/exchangepp.o : src/exchangepp.f90 src/module_vars.o src/module_mpi.o src/module_multidata.o 
src/exchangesca.o : src/exchangesca.f90 src/module_vars.o src/module_mpi.o src/module_multidata.o 
src/exchangeu.o : src/exchangeu.f90 src/module_vars.o src/module_mpi.o src/module_multidata.o 
src/exchangev.o : src/exchangev.f90 src/module_vars.o src/module_mpi.o src/module_multidata.o 
src/exchangew.o : src/exchangew.f90 src/module_vars.o src/module_mpi.o src/module_multidata.o 
src/flosol.o : src/flosol.f90 src/module_vars_pt.o src/module_multidata.o src/module_mpi.o src/module_vars.o src/MPI_pt.o
src/hdf5_io.o : src/hdf5_io.f90 
src/imb.o : src/imb.f90 src/module_mpi.o src/module_multidata.o src/module_vars.o 
src/init_particle.o : src/init_particle.f90 src/hdf5_io.o src/module_vars_pt.o src/module_vars.o src/module_mpi.o src/module_multidata.o 
src/initial.o : src/initial.f90 src/multiflow3d_sem.o src/module_LSM.o src/module_mpi.o src/module_multidata.o src/module_vars.o 
src/io.o : src/io.f90 src/json_io.o src/module_multidata.o src/module_mpi.o
src/json_io.o : src/json_io.f90 
src/localparameters.o : src/localparameters.f90 src/module_multidata.o src/module_mpi.o src/module_vars.o 
src/log_law.o : src/log_law.f90 src/module_multidata.o src/module_vars.o 
src/lsm.o : src/lsm.f90 src/module_mpi.o src/module_multidata.o src/module_LSM.o src/module_vars.o 
src/mgsolver.o : src/mgsolver.f90 src/module_LSM.o src/module_multidata.o src/module_vars.o 
src/module_LSM.o : src/module_LSM.f90 
src/module_mpi.o : src/module_mpi.f90 
src/module_multidata.o : src/module_multidata.f90 
src/module_vars.o : src/module_vars.f90 
src/module_vars_pt.o : src/module_vars_pt.f90 
src/multiflow3d_sem.o : src/multiflow3d_sem.f90 src/module_mpi.o src/module_multidata.o src/module_vars.o 
src/newsolv_mg.o : src/newsolv_mg.f90 src/module_multidata.o src/module_mpi.o src/module_vars.o 
src/post.o : src/post.f90 src/module_LSM.o src/module_vars.o src/module_multidata.o 
src/press.o : src/press.f90 src/module_LSM.o src/module_multidata.o src/module_mpi.o src/module_vars.o 
src/roughness_function.o : src/roughness_function.f90 src/module_multidata.o src/module_mpi.o src/module_vars.o 
src/rungek.o : src/rungek.f90 src/module_LSM.o src/module_multidata.o src/module_mpi.o src/module_vars.o 
src/sediment.o : src/sediment.f90 src/module_multidata.o src/module_mpi.o src/module_vars.o 
src/shapes.o : src/shapes.f90 src/module_mpi.o src/imb.o src/module_multidata.o src/module_vars.o 
src/sipsol.o : src/sipsol.f90 src/module_mpi.o src/module_multidata.o src/module_vars.o 
src/timesig.o : src/timesig.f90 src/module_mpi.o src/module_vars.o src/module_multidata.o 
src/wall_function.o : src/wall_function.f90 src/module_multidata.o src/module_vars.o 
src/weno.o : src/weno.f90 src/module_multidata.o src/module_vars.o 
app/fdstag.o : app/fdstag.f90 src/module_multidata.o src/io.o src/module_vars.o src/module_mpi.o 
tests/main.o : tests/main.f90 tests/test_hdf5_io.o tests/test_io.o tests/test_json_io.o src/module_mpi.o tests/test_MPI_pt.o 
tests/test_hdf5_io.o : tests/test_hdf5_io.f90 src/hdf5_io.o 
tests/test_io.o : tests/test_io.f90 src/io.o src/module_mpi.o
tests/test_json_io.o : tests/test_json_io.f90 src/json_io.o 
tests/test_MPI_pt.o : tests/test_MPI_pt.f90 src/MPI_pt.o
