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
module_mpi.o\
hdf5_io.o \
json_io.o \
module_multidata.o\
module_vars.o\
module_vars_pt.o\
multiflow3d_sem.o\
module_LSM.o\
imb.o\
shapes.o\
fdstag.o\
initial.o\
init_particle.o\
localparameters.o\
alloc_dom.o\
post.o\
flosol.o\
checkdt.o\
bounds.o\
bounds_keps.o\
sipsol.o\
convection.o\
diffusion.o\
newsolv_mg.o\
mgsolver.o\
wall_function.o\
log_law.o\
alloc_pt.o\
MPI_pt.o\
delta_func.o\
collision.o\
LPT.o\
timesig.o\
weno.o\
energy.o\
press.o\
roughness_function.o\
rungek.o\
averaging.o\
eddyvis_smag.o\
eddyvis_wale.o\
eddyvis_1eqn.o\
eddyvis_keps.o\
exchange_bc.o\
exchangep.o\
exchangepp.o\
exchangesca.o\
exchange.o\
exchangeu.o\
exchangev.o\
exchangew.o\
exchange_phi.o\
bounds_lsm.o\
lsm.o\
sediment.o\
io.o

test_objects = \
hdf5_io.o \
json_io.o \
io.o \
tests/test_hdf5_io.o \
tests/test_json_io.o \
tests/test_io.o \
tests/main.o

all: test

test: M3D_v2.exe tests/tests.exe
	@cd ./tests/ && ./tests.exe && cd ..

.SUFFIXES: .f90

M3D_v2.exe: $(objects) 
	$(F90) $(objects) $(LOPTIONS) \
	-I./$(JSON_FORTRAN_INCLUDE_PATH) \
	-L./$(JSON_FORTRAN_LIBRARY_PATH) \
	-Wl,-rpath,$(CURDIR)/$(JSON_FORTRAN_LIBRARY_PATH) \
	-Wl,-rpath,$(CURDIR)/$(HDF5_LIBRARY_PATH) \
	-I./$(HDF5_INCLUDE_PATH) \
	-L./$(HDF5_LIBRARY_PATH) \
	-ljsonfortran -lhdf5 -lhdf5_fortran -o M3D_v2.exe \

tests/%.o: tests/%.f90
	$(F90) $(OPTIONS) -c $< -o $@ \
	-I./$(JSON_FORTRAN_INCLUDE_PATH) \
	-I./$(TEST_DRIVE_INCLUDE_PATH) \
	-I./$(HDF5_INCLUDE_PATH)

%.o: %.f90
	$(F90) $(OPTIONS) -c $< -o $@ \
	-I./$(JSON_FORTRAN_INCLUDE_PATH) \
	-I./$(HDF5_INCLUDE_PATH)

tests/tests.exe: $(test_objects)
	$(F90) $(test_objects) $(LOPTIONS) \
	-I./$(JSON_FORTRAN_INCLUDE_PATH) -L./$(JSON_FORTRAN_LIBRARY_PATH) \
	-I./$(TEST_DRIVE_INCLUDE_PATH) -L./$(TEST_DRIVE_LIBRARY_PATH) \
	-I./$(HDF5_INCLUDE_PATH) -L./$(HDF5_LIBRARY_PATH) \
	-Wl,-rpath,$(CURDIR)/$(JSON_FORTRAN_LIBRARY_PATH) \
	-Wl,-rpath,$(CURDIR)/$(TEST_DRIVE_LIBRARY_PATH) \
	-Wl,-rpath,$(CURDIR)/$(HDF5_LIBRARY_PATH) \
	-ljsonfortran -ltest-drive -lhdf5 -lhdf5_fortran \
	-o tests/tests.exe

clean:
	rm -rf *.o *.mod tests/*.o tests/tests.exe
