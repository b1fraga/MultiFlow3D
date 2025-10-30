!##########################################################################
        module mpi
!##########################################################################
        implicit none
        include 'mpif.h'

        !  SAVE
        integer           :: nprocs     ! nr. of processors
        integer           :: myrank      ! nr of this processor
        integer           :: ierr
        integer           :: MPI_FLT  
        integer           :: status(MPI_STATUS_SIZE)
        real ( kind =8 )  :: wtime
        real ( kind =8 )  :: wtime2

        contains

!#########################################################################
        subroutine init_parallelisation
!#########################################################################
        implicit none

!  initialize MPI.
        call MPI_INIT ( ierr )
!  Get the number of processes.
        call MPI_COMM_SIZE ( MPI_COMM_WORLD, nprocs, ierr )
!  Get the individual process iD.
        call MPI_COMM_RANK ( MPI_COMM_WORLD, myrank, ierr )

        MPI_FLT   = MPI_DOUBLE_PRECISION

        if ( myrank == 0 ) then
           wtime = MPI_WTIME ( )
           write ( *, '(a)' ) '========================================'
           write ( *, '(a)' ) '...............MultiFlow3D..............'
           write ( *, '(a)' ) '========================================'
           write ( *, '(a)' )'~~~~~~~~~~~~~~~~~~~!~~~~~~~~~~~~~~~~~~~~'
           write ( *, '(a)' )'~~~~~~~~~~~~~~~~~~/.\~~~~~~~~~~~~~~~~~~~'                                                                                
           write ( *, '(a)' )'~~~~~~~~~~~~~~~~~/...\~~~~~~~~~~~~~~~~~~'
           write ( *, '(a)' )'~~~~~~~~~~~~~~~~{.....}~~~~~~~~~~~~~~~~~'
           write ( *, '(a)' )'~~~~~~~~~~~~~~~~~`~.~´~~~~~~~~~~~~~~~~~~'
           write ( *, '(a)' )'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
           write ( *, '(a)' ) 'Created by Bruño Fraga'
           write ( *, '(a)' ) 'Contributors: Boyang Chen'
           write ( *, '(a)' ) 'Contributors: Aleksandra Monka'
           write ( *, '(a)' ) 'Contributors: Yu Zhang'
           write ( *, '(a)' ) 'Contributors: Zhen Liu'
           write ( *, '(a)' ) '\\\\\\\University of Birmingham///////'
           write ( *, '(a)' ) ' '                
           write ( *, '(a)' ) 'Based on the 3DFM code by Mehtap Cahveci'
           write ( *, '(a)' ) ' '                      
           write ( *, '(a)' ) '3D Finite Diff Navier-Stokes solver'
           write ( *, '(a)' ) 'FORTRAN90/MPI/OpenMP'
           write ( *, '(a)' ) ' '
           write ( *, '(a,i8)' ) 'The number of processes is ', nprocs
           write ( *, '(a)' ) ' '
        end if

        end subroutine init_parallelisation
!##########################################################################
        subroutine end_parallelisation
!##########################################################################
        implicit none

        if ( myrank == 0 ) then
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'MultiFlow3D - Master process:'
           write (*,'(a)')'Normal end of execution, good luck.'
           wtime = MPI_WTIME ( ) - wtime
           write ( *, '(a)' ) ' '
           write (*,'(a,g14.6,a)') 'Elapsed wall clock time = ',wtime,
     &'seconds.'
        end if

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        call MPI_FINALIZE ( ierr )

        end subroutine end_parallelisation
!##########################################################################

        end module mpi
!##########################################################################
