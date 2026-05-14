!##########################################################################
      program fdstag
!##########################################################################
          use multiflow3d_mpi
          use vars
          use, intrinsic :: iso_fortran_env, only: dp => real64
#if USE_JSON == 1
          use io
          use vars
          use multidata
          use multiflow3d_mpi
#endif
          implicit none
          integer :: ib

          call init_parallelisation

          call read_mdmap
#if USE_JSON == 1
          call read_control_file("inputs/control.cin",Keyword,L_n,&
       g_dx,g_dy,g_dz,ubulk,rrey,Pr,Sc_t,beta,&
       gx,gy,gz,dens,conv_sch,diff_sch,differencing,&
       solver,ngrid_input,mg_itrsch,&
       maxcy,irestr,iproln,&
       dt,safety_factor,eps,fric,&
       L_dt,lrestart,reinitmean,LTRANSIENT,&
       sweeps,itime_end,n_out,tsteps_pt,niter,&
       nswp(1),nswp(2),nswp(3),nswp(4),bc_w,&
       bc_e,bc_s,bc_n,&
       bc_b,bc_t,&
       save_inflow,time_averaging,SGS,&
       ITMAX_PI,UPROF_SEM,ITMAX_SEM,&
       TI_SEM,t_start_averaging1,t_start_averaging2,&
       noise,Th,Tc,Tinit,SGS_model,LMR,pl_ex,&
       LIMB,LENERGY,LROUGH,LPT,l_LSM,L_LSMbase,LSCALAR,LAS,LNonNewt,&
       Tbc_w,Tbc_e,Tbc_s,Tbc_n,&
       Tbc_n,Tbc_t,n_unstpt,&
       id_unst,i_unst,j_unst,k_unst)
       Re = 1.0_dp/rrey
       if (bc_w==5) pressureforce=.TRUE.
       if (.not.LPT) np=0

       if (trim(L_n)=="n") fric=fric**0.33_dp


       do ib=1,nbp
          dom(ib)%bc_west=bc_w
          dom(ib)%bc_east=bc_e
          dom(ib)%bc_south=bc_s
          dom(ib)%bc_north=bc_n
          dom(ib)%bc_bottom=bc_b
          dom(ib)%bc_top=bc_t
          dom(ib)%Tbc_west=Tbc_w
          dom(ib)%Tbc_east=Tbc_e
          dom(ib)%Tbc_south=Tbc_s
          dom(ib)%Tbc_north=Tbc_n
          dom(ib)%Tbc_bottom=Tbc_b
          dom(ib)%Tbc_top=Tbc_t
          dom(ib)%ngrid=ngrid_input
          if (dom(ib)%bc_west==7) read_inflow=.true.
       end do

       if(differencing==3 .and. pl_ex/=2) then
          pl_ex=2
          print*,"error: you select WENO but do not assign", &
               "  correct number of ghost planes,now it is corrected to 2"
       end if

       if(SGS .and. sgs_model==3 .and. pl_ex/=2) then
          pl_ex=2
          print*,"error: you select 1-EQN model but do not assign", &
               "  correct number of ghost planes,now it is corrected to 2"
       end if

       if (L_LSM .and. solver==1) then
          if (myrank==0) then
             print*,"Error: SIP solver not presently compatible with LSM"
          end if
          stop
       end if

       if (L_LSM .and. differencing/=3) then
          if (myrank==0) then
             print*,"Error: WENO differencing must be used with LSM"
          end if
          stop
       end if
#else
          call read_control
#endif
          call read_infodom

          call alloc_dom

          call localparameters

          call initial

!        IF (L_LSM)  CALL initial_LSM_3D_channel

          call initflowfield

          IF (LROUGH)  THEN                             !Richard2015
          IF (.not.LRESTART) THEN
          call init_rough
          ELSE
          call rough_restart
          END IF
          END IF

          if (LIMB) call imb_initial
          if (LIMB) call PartLocMPI                         !Pablo2015

          call iniflux

          if(.not.LRESTART) then
          if (time_averaging) then
          call update_mean
          if (noise>0.0_dp) call add_noise(noise)
          end if
          end if

          if (LPT) then                     !Brunho2013
          if (myrank==0)  open(unit=202, file="particle_log")
          call init_particle
          end if

          if ((solver==2).and.(.not.L_LSM)) call coeff

          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
          if(myrank==0) then
          write (numfile,*) "============START ITERATIONS========="
          write (6,*) "============START ITERATIONS========="
          end if

          call flosol

          call end_parallelisation

      end program fdstag
!##########################################################################
