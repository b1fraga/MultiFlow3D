#############################################################
F90=mpif90
OPTIONS    =  -c -fdefault-real-8 -fdefault-double-8  -O2 -fbacktrace -fallow-argument-mismatch -g -fopenmp
LOPTIONS   = -O2 -fopenmp
INCLUDE_PATH := $(shell realpath $$(dirname $$(find . -path "*json_module.mod*")))
LIBRARY_PATH := $(shell realpath $$(dirname $$(find . -path "*json-fortran/libjson-fortran.a.log*")))
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
sediment.o

.SUFFIXES: .f90

.f90.o:
	$(F90) $(OPTIONS) -I$(INCLUDE_PATH) -L$(LIBRARY_PATH) -ljson-fortran -o $@ $<

M3D_v2.exe: $(objects) 
	$(F90) $(objects) $(LOPTIONS) -I$(INCLUDE_PATH) -L$(LIBRARY_PATH) -ljson-fortran -o M3D_v2.exe \

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
fdstag.o : fdstag.f90 module_vars.o module_mpi.o 
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
