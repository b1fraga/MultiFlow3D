#############################################################
F90=mpif90
OPTIONS    =  -c -cpp -O2 -g -fopenmp
ifneq (,$(findstring GNU,$(shell $(F90) --version)))
   OPTIONS += -std=f2018
endif
LOPTIONS   = -O2 -fopenmp
JSON_FORTRAN_INCLUDE_PATH := build_json_fortran/include/
JSON_FORTRAN_LIBRARY_PATH := build_json_fortran/lib/
TEST_DRIVE_INCLUDE_PATH := build_test_drive/include/
TEST_DRIVE_LIBRARY_PATH := build_test_drive/
HDF5_INCLUDE_PATH := build_hdf5/mod/shared/
HDF5_LIBRARY_PATH := build_hdf5/bin
##############################################################

objects = \
module_vars.o\
module_multidata.o\
module_mpi.o\
module_vars_pt.o\
module_SEM.o\
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
SEM.o\
sediment.o \
json_io.o \
io.o \
hdf5_io.o

test_objects = \
tests/test_json_io.o \
tests/test_io.o \
tests/test_hdf5_io.o \
json_io.o \
io.o \
hdf5_io.o \
tests/main.o

all: test

test: tests.exe
	@cd tests && \
	./tests.exe && \
	cd ..

.SUFFIXES: .f90

.f90.o:
	$(F90) $(OPTIONS) -I./$(JSON_FORTRAN_INCLUDE_PATH) -L./$(JSON_FORTRAN_LIBRARY_PATH) \
			   -I./$(HDF5_INCLUDE_PATH) -L./$(HDF5_LIBRARY_PATH) -DUSE_HDF5=0 -DUSE_JSON=1 -ljsonfortran -lhdf5 -lhdf5_fortran -o $@ $<

M3D_v2.exe: $(objects) 
	$(F90) $(objects) $(LOPTIONS) -I./$(JSON_FORTRAN_INCLUDE_PATH) -L./$(JSON_FORTRAN_LIBRARY_PATH) \
	-Wl,-rpath,$(CURDIR)/$(JSON_FORTRAN_LIBRARY_PATH) -Wl,-rpath,$(CURDIR)/$(HDF5_LIBRARY_PATH) \
				      -I./$(HDF5_INCLUDE_PATH) -L./$(HDF5_LIBRARY_PATH) -ljsonfortran -lhdf5 -lhdf5_fortran -o M3D_v2.exe \

tests/%.o: tests/%.f90
	$(F90) $(LOPTIONS) $(OPTIONS) -I./$(JSON_FORTRAN_INCLUDE_PATH) -L./$(JSON_FORTRAN_LIBRARY_PATH) -ljsonfortran \
					-I./$(TEST_DRIVE_INCLUDE_PATH) -L./$(TEST_DRIVE_LIBRARY_PATH) \
                                        -I./$(HDF5_INCLUDE_PATH) -L./$(HDF5_LIBRARY_PATH) -DUSE_HDF5=0 -DUSE_JSON=1 -ljsonfortran -ltest-drive -lhdf5 -lhdf5_fortran -c -o $@ $<

tests.exe: $(test_objects) M3D_v2.exe
	$(F90) $(test_objects) $(LOPTIONS) -I./$(JSON_FORTRAN_INCLUDE_PATH) -L./$(JSON_FORTRAN_LIBRARY_PATH) \
			-Wl,-rpath,$(CURDIR)/$(JSON_FORTRAN_LIBRARY_PATH) -Wl,-rpath,$(CURDIR)/$(TEST_DRIVE_LIBRARY_PATH) \
			-Wl,-rpath,$(CURDIR)/$(HDF5_LIBRARY_PATH) -I./$(TEST_DRIVE_INCLUDE_PATH) -L./$(TEST_DRIVE_LIBRARY_PATH) \
            -I./$(HDF5_INCLUDE_PATH) -L./$(HDF5_LIBRARY_PATH) -ljsonfortran -ltest-drive -lhdf5 -lhdf5_fortran -o tests/tests.exe

clean:
	rm -rfv *.o *.mod M3D_v2.exe


module_mpi.o : module_mpi.f90 
module_multidata.o : module_multidata.f90 
module_vars.o : module_vars.f90 
module_vars_pt.o : module_vars_pt.f90 
module_LSM.o:module_LSM.f90
module_SEM.o:module_SEM.f90
lsm.o : lsm.f90 module_mpi.o module_multidata.o module_vars.o module_LSM.o 
alloc_dom.o : alloc_dom.f90 module_mpi.o module_vars.o module_multidata.o 
alloc_pt.o : alloc_pt.f90 module_vars_pt.o module_multidata.o module_mpi.o module_vars.o 
averaging.o : averaging.f90 module_vars.o module_multidata.o 
bounds.o : bounds.f90 imb.o module_mpi.o module_multidata.o module_vars.o 
bounds_keps.o : bounds_keps.f90 module_multidata.o module_vars.o 
bounds_lsm.o : bounds_lsm.f90 module_multidata.o module_vars.o module_LSM.o 
checkdt.o : checkdt.f90 module_multidata.o module_mpi.o module_vars.o module_LSM.o  
convection.o : convection.f90 module_multidata.o module_vars.o 
delta_func.o : delta_func.f90 
diffusion.o : diffusion.f90 module_multidata.o module_mpi.o module_vars.o 
eddyvis_1eqn.o : eddyvis_1eqn.f90 module_multidata.o module_vars.o 
eddyvis_smag.o : eddyvis_smag.f90 module_multidata.o module_vars.o 
eddyvis_wale.o : eddyvis_wale.f90 module_multidata.o module_vars.o 
eddyvis_keps.o : eddyvis_keps.f90 module_multidata.o module_vars.o 
energy.o : energy.f90 module_multidata.o module_mpi.o module_vars.o 
exchange_bc.o : exchange_bc.f90 module_vars.o module_mpi.o module_multidata.o 
exchange_bcphi.o : exchange_bcphi.f90 module_vars.o module_mpi.o module_multidata.o 
exchange.o : exchange.f90 exchange_phi.o module_mpi.o module_multidata.o module_vars.o 
exchangep.o : exchangep.f90 module_vars.o module_mpi.o module_multidata.o 
exchange_phi.o : exchange_phi.f90 module_vars.o module_mpi.o module_multidata.o 
exchangepp.o : exchangepp.f90 module_vars.o module_mpi.o module_multidata.o 
exchangesca.o : exchangesca.f90 module_vars.o module_mpi.o module_multidata.o 
exchangeu.o : exchangeu.f90 module_vars.o module_mpi.o module_multidata.o 
exchangev.o : exchangev.f90 module_vars.o module_mpi.o module_multidata.o 
exchangew.o : exchangew.f90 module_vars.o module_mpi.o module_multidata.o 
fdstag.o : fdstag.f90 module_vars.o module_mpi.o io.o
flosol.o : flosol.f90 module_vars_pt.o module_multidata.o module_mpi.o module_vars.o 
imb.o : imb.f90 module_mpi.o module_multidata.o module_vars.o module_LSM.o
initial.o : initial.f90 module_mpi.o module_multidata.o module_vars.o 
init_particle.o : init_particle.f90 module_vars_pt.o module_vars.o module_mpi.o module_multidata.o 
localparameters.o : localparameters.f90 module_multidata.o module_mpi.o module_vars.o 
collision.o: LPT.f90 module_vars_pt.o module_vars.o module_mpi.o module_multidata.o 
LPT.o : LPT.f90 module_vars_pt.o module_vars.o module_mpi.o module_multidata.o 
sediment.o: sediment.f90 module_multidata.o module_vars.o module_mpi.o
mgsolver.o : mgsolver.f90 module_multidata.o module_vars.o module_LSM.o 
MPI_pt.o : MPI_pt.f90 module_vars_pt.o module_multidata.o module_mpi.o module_vars.o 
newsolv_mg.o : newsolv_mg.f90 module_multidata.o module_mpi.o module_vars.o 
post.o : post.f90 module_vars.o module_multidata.o 
press.o : press.f90 module_multidata.o module_mpi.o module_vars.o 
roughness_function.o : roughness_function.f90 module_multidata.o module_mpi.o module_vars.o 
rungek.o : rungek.f90 module_multidata.o module_mpi.o module_vars.o 
shapes.o : shapes.f90 module_mpi.o imb.o module_multidata.o module_vars.o 
sipsol.o : sipsol.f90 module_mpi.o module_multidata.o module_vars.o 
timesig.o : timesig.f90 module_mpi.o module_vars.o module_multidata.o 
wall_function.o : wall_function.f90 module_multidata.o module_vars.o 
log_law.o : log_law.f90 module_multidata.o module_vars.o 
weno.o : weno.f90 module_multidata.o module_vars.o 
SEM.o : SEM.f90 module_multidata.o module_vars.o module_SEM.o module_mpi.o
json_io.o : json_io.f90
io.o : io.f90 json_io.o
hdf5_io.o : hdf5_io.f90
tests/test_json_io.o : tests/test_json_io.f90 json_io.o
tests/test_io.o : tests/test_io.f90 io.o
tests/test_hdf5_io.o: tests/test_hdf5_io.f90 hdf5_io.o
tests/main.o : tests/main.f90 tests/test_json_io.o tests/test_io.o json_io.o io.o
