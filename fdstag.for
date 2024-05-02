!##########################################################################
        program fdstag
!##########################################################################
        use mpi
        use vars
        implicit none

        call init_parallelisation

        call read_mdmap

        call read_control

        call read_infodom

        call alloc_dom

        call localparameters

        call initial


	     IF (L_LSM)  CALL initial_LSM_3D_channel

        call initflowfield


        IF (LROUGH)  THEN								!Richard2015
          IF (.not.LRESTART) THEN
          call init_rough
          ELSE
          call rough_restart
         END IF
        END IF

        if (LIMB) call imb_initial

        if (LIMB) call PartLocMPI							!Pablo2015

        call iniflux

        if(.not.LRESTART) then
           if (time_averaging) then
              call update_mean
	      if (noise.gt.0.0) call add_noise(noise)
           end if
        end if

	   if (LPT) then						!Brunho2013
			if (myrank.eq.0)  open(unit=202, file='particle_log') 
         call init_particle	
	   endif

        if ((solver.eq.2).and.(.not.L_LSM)) call coeff
      
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           write (numfile,*) '============START ITERATIONS========='
           write (6,*) '============START ITERATIONS========='
        end if

        call flosol

        call end_parallelisation

        end program
!##########################################################################
