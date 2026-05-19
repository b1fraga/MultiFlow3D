module io
#if USE_JSON == 1
  use multidata, only: dom_id, dom_indid, imbinblk, dom_ad, dom, nbp, nbpmax, num_domains
  use vars, only: numfile
  use multiflow3d_mpi, only: ierr, mpi_comm_world, mpi_integer, mpi_max, myrank, nprocs
  use json_io, only: json_read
#else
  use multidata, only: id_unst, dom, i_unst, dom_id, dom_indid, imbinblk, j_unst, dom_ad, &
                       k_unst, nbp, nbpmax, num_domains
  use vars, only: numfile, bc_b, bc_e, bc_n, bc_s, bc_t, bc_w, beta, conv_sch, dens, diff_sch, &
                  differencing, dt, fric, g_dx, g_dy, g_dz, gx, gy, gz, iproln, irestr, itime_end, &
                  itmax_pi, itmax_sem, keyword, l_dt, l_lsm, l_lsmbase, l_n, las, lenergy, limb, &
                  lmr, lnonnewt, lpt, lrestart, lrough, lscalar, ltransient, maxcy, mg_itrsch, &
                  n_unstpt, ngrid_input, noise, np, pl_ex, pr, pressureforce, re, read_inflow, &
                  reinitmean, rrey, safety_factor, save_inflow, sc_t, sgs, sgs_model, solver, &
                  t_start_averaging1, t_start_averaging2, tbc_b, tbc_e, tbc_n, tbc_s, tbc_t, &
                  tc, th, ti_sem, time_averaging, tinit, tsteps_pt, ubulk, uprof_sem, nswp, &
                  niter, n_out, sweeps, tbc_w, eps
  use multiflow3d_mpi, only: ierr, mpi_comm_world, mpi_integer, mpi_max, myrank, nprocs, myrank
#endif
  use, intrinsic :: iso_fortran_env, only: dp => real64

  implicit none
  private
#if USE_JSON == 1
  public :: read_mdmap, read_control_file
#else
  public :: read_mdmap, read_control
#endif

contains

#if USE_JSON == 1
  subroutine read_control_file(input_file,Keyword,type_of_friction,&
       dx,dy,dz,Ubulk,kinematic_visc,Pr,turb_Schmidt,beta,&
       gx,gy,gz,dens,convection_scheme,diffusion_scheme,differencing,&
       solver,multigrid_step,multigrid_iteration_scheme,&
       multigrid_maximum_iteration_per_time_step,restriction_iter,prolongation_iter,&
       dt,safety_factor,eps,Friction_coefficient,&
       variable_dt,restart,reinitmean,LTRANSIENT,&
       sweeps,itime_end,n_out,results_output,niter,&
       nswp_1,nswp_2,nswp_3,nswp_4,West_Boundary_Condition,&
       East_Boundary_Condition,South_Boundary_Condition,North_Boundary_Condition,&
       Bottom_Boundary_Condition,Top_Boundary_Condition,&
       save_inflow_data,time_averaging,SGS_model,&
       number_of_inlets,velocity_profile,Number_inlet_profiles,&
       Turbulence_intensity,t_start_averaging1,t_start_averaging2,&
       noise,Th,Tc,Tinit,SGS_model_value,LMR,pl_ex,&
       LIMB,LENERGY,LROUGH,LPT,LSM,L_LSMbase,LSCALAR,LActiveScalar,LNonNewt,&
       West_Energy_BC,East_Energy_BC,South_Energy_BC,North_Energy_BC,&
       Bottom_Energy_BC,Top_Energy_BC,num_of_time_series_points,&
       time_series_point_1,time_series_point_2,time_series_point_3,time_series_point_4)
    !! Reads control.json assigning the variables associated output to variables of the subroutine
    character(len=*), intent(in) :: input_file
    character(len=80),  intent(out) :: Keyword,type_of_friction
    real(dp), intent(out) :: dx,dy,dz,Ubulk,kinematic_visc,Pr,turb_Schmidt,beta
    real(dp), intent(out) :: gx,gy,gz, dens
    integer, intent(out) :: convection_scheme,diffusion_scheme,differencing
    integer, intent(out) :: solver,multigrid_step,multigrid_iteration_scheme
    integer, intent(out) :: multigrid_maximum_iteration_per_time_step,restriction_iter
    integer, intent(out) :: prolongation_iter
    real(dp), intent(out) :: dt,safety_factor,eps,Friction_coefficient
    logical, intent(out) :: variable_dt,restart,reinitmean,LTRANSIENT
    integer, intent(out) :: sweeps,itime_end,n_out,results_output,niter
    integer, intent(out) :: nswp_1,nswp_2,nswp_3,nswp_4
    integer, intent(out) :: West_Boundary_Condition,East_Boundary_Condition
    integer, intent(out) :: South_Boundary_Condition,North_Boundary_Condition
    integer, intent(out) :: Bottom_Boundary_Condition,Top_Boundary_Condition
    logical, intent(out) :: save_inflow_data,time_averaging,SGS_model
    integer, intent(out) :: number_of_inlets, velocity_profile,Number_inlet_profiles
    real(dp), intent(out) :: Turbulence_intensity,t_start_averaging1,t_start_averaging2
    real(dp), intent(out) :: noise,Th,Tc, Tinit
    integer, intent(out) :: SGS_model_value,LMR,pl_ex
    logical, intent(out) :: LIMB,LENERGY,LROUGH,LPT,LSM,L_LSMbase,LSCALAR,LActiveScalar,LNonNewt
    integer, intent(out) ::  West_Energy_BC,East_Energy_BC,South_Energy_BC,North_Energy_BC
    integer, intent(out) :: Bottom_Energy_BC,Top_Energy_BC
    integer, intent(out) :: num_of_time_series_points
    integer,allocatable, intent(out) :: time_series_point_1(:),time_series_point_2(:)
    integer,allocatable, intent(out) :: time_series_point_3(:)
    integer,allocatable, intent(out) :: time_series_point_4(:)

    ! Read numerical paramters
    call json_read(input_file,"Keyword",Keyword)
    call json_read(input_file,"Ubulk",Ubulk)
    call json_read(input_file,"dx",dx)
    call json_read(input_file,"dy",dy)
    call json_read(input_file,"dz",dz)
    call json_read(input_file,"dens",dens)
    call json_read(input_file,"kinematic visc",kinematic_visc)
    call json_read(input_file,"Pr",Pr)
    call json_read(input_file,"turb Schmidt",turb_Schmidt)
    call json_read(input_file,"beta",beta)
    call json_read(input_file,"gx",gx)
    call json_read(input_file,"gy",gy)
    call json_read(input_file,"gz",gz)
    call json_read(input_file,"convection_scheme",convection_scheme)
    call json_read(input_file,"diffusion_scheme",diffusion_scheme)
    call json_read(input_file,"differencing",differencing)
    call json_read(input_file,"solver",solver)
    call json_read(input_file,"multigrid_step",multigrid_step)
    call json_read(input_file,"multigrid iteration scheme",multigrid_iteration_scheme)
    call json_read(input_file,"multigrid maximum iteration per time step",&
         multigrid_maximum_iteration_per_time_step)
    call json_read(input_file,"restriction iter",restriction_iter)
    call json_read(input_file,"prolongation iter",prolongation_iter)
    call json_read(input_file,"dt",dt)
    call json_read(input_file,"variable_dt",variable_dt)
    call json_read(input_file,"sweeps",sweeps)
    call json_read(input_file,"safety_factor",safety_factor)
    call json_read(input_file,"itime_end",itime_end)
    call json_read(input_file,"restart",restart)
    call json_read(input_file,"reinitmean",reinitmean)
    call json_read(input_file,"n_out",n_out)
    call json_read(input_file,"LTRANSIENT",LTRANSIENT)
    call json_read(input_file,"results output",results_output)
    call json_read(input_file,"niter",niter)
    call json_read(input_file,"eps",eps)
    call json_read(input_file,"nswp_1",nswp_1)
    call json_read(input_file,"nswp_2",nswp_2)
    call json_read(input_file,"nswp_3",nswp_3)
    call json_read(input_file,"nswp_4",nswp_4)
    ! Flow boundary conditions
    call json_read(input_file,"West Boundary Condition",West_Boundary_Condition)
    call json_read(input_file,"East Boundary Condition",East_Boundary_Condition)
    call json_read(input_file,"South Boundary Condition",South_Boundary_Condition)
    call json_read(input_file,"North Boundary Condition",North_Boundary_Condition)
    call json_read(input_file,"Bottom Boundary Condition",Bottom_Boundary_Condition)
    call json_read(input_file,"Top Boundary Condition",Top_Boundary_Condition)
    call json_read(input_file,"type of friction",type_of_friction)
    call json_read(input_file,"Friction coefficient",Friction_coefficient)
    call json_read(input_file,"save inflow data",save_inflow_data)
    call json_read(input_file,"number of inlets",number_of_inlets)
    ! Synthetic Eddy Method
    call json_read(input_file,"velocity profile",velocity_profile)
    call json_read(input_file,"Turbulence intensity",Turbulence_intensity)
    call json_read(input_file,"Number inlet profiles",Number_inlet_profiles)
    ! Modelling Options
    call json_read(input_file,"time_averaging",time_averaging)
    call json_read(input_file,"t_start_averaging1",t_start_averaging1)
    call json_read(input_file,"t_start_averaging2",t_start_averaging2)
    call json_read(input_file,"noise",noise)
    call json_read(input_file,"SGS-model",SGS_model)
    call json_read(input_file,"SGS-model_value",SGS_model_value)
    call json_read(input_file,"LMR",LMR)
    call json_read(input_file,"LIMB",LIMB)
    call json_read(input_file,"LENERGY",LENERGY)
    call json_read(input_file,"LROUGH",LROUGH)
    call json_read(input_file,"LPT",LPT)
    call json_read(input_file,"LSM",LSM)
    call json_read(input_file,"L_LSMbase",L_LSMbase)
    call json_read(input_file,"LSCALAR",LSCALAR)
    call json_read(input_file,"LActiveScalar",LActiveScalar)
    call json_read(input_file,"LNonNewt",LNonNewt)
    call json_read(input_file,"pl_ex",pl_ex)
    call json_read(input_file,"Th",Th)
    call json_read(input_file,"Tc",Tc)
    call json_read(input_file,"Tinit",Tinit)
    ! Energy boundary conditions
    call json_read(input_file,"West_Energy_BC",West_Energy_BC)
    call json_read(input_file,"East_Energy_BC",East_Energy_BC)
    call json_read(input_file,"South_Energy_BC",South_Energy_BC)
    call json_read(input_file,"North_Energy_BC",North_Energy_BC)
    call json_read(input_file,"Bottom_Energy_BC",Bottom_Energy_BC)
    call json_read(input_file,"Top_Energy_BC",Top_Energy_BC)
    ! Time series
    call json_read(input_file,"num of time series points",&
         num_of_time_series_points)
    allocate(time_series_point_1(num_of_time_series_points),&
         time_series_point_2(num_of_time_series_points), &
         time_series_point_3(num_of_time_series_points),&
         time_series_point_4(num_of_time_series_points))
    call json_read(input_file,"time_series_point_1",time_series_point_1)
    call json_read(input_file,"time_series_point_2",time_series_point_2)
    call json_read(input_file,"time_series_point_3",time_series_point_3)
    call json_read(input_file,"time_series_point_4",time_series_point_4)
  end subroutine read_control_file
#else
  subroutine read_control
!##########################################################################
          implicit none
          integer :: mgi,mgj,mgk,pow2,ib,i,j

          open (unit=12, file="input/control.cin")
!------DOMAIN SIZE AND DISCRETIZATION -------------------------------------
          read (12,*)
          read (12,*) keyword,ubulk
          read (12,*) g_dx,g_dy,g_dz
          read (12,*) dens,rrey,Pr,Sc_t,beta
          Re=1.d0/rrey
          read (12,*) gx,gy,gz
          read (12,*) conv_sch
          read (12,*) diff_sch
          read (12,*) differencing
          read (12,*) solver
          read (12,*) ngrid_input,mg_itrsch
          read (12,*) maxcy,irestr,iproln
          read (12,*) dt,L_dt,sweeps,safety_factor
          read (12,*) itime_end,LRESTART,reinitmean,n_out
          read (12,*) LTRANSIENT,tsteps_pt
          read (12,*) niter,eps,nswp(1),nswp(2),nswp(3),nswp(4)
          read (12,*)
          read (12,*) bc_w
          read (12,*) bc_e
          read (12,*) bc_s
          read (12,*) bc_n
          read (12,*) bc_b
          read (12,*) bc_t
          if (bc_w==5) pressureforce=.TRUE.
          read (12,*) L_n,fric
          read (12,*) save_inflow,ITMAX_PI
          read (12,*)
          read (12,*) UPROF_SEM         !Pablo 15/12/2015
          read (12,*) TI_SEM
          read (12,*) ITMAX_SEM
          read (12,*)
          read (12,*) time_averaging,t_start_averaging1, &
    t_start_averaging2,noise
          read (12,*) SGS,sgs_model
          read (12,*) LMR
          read (12,*) LIMB,LENERGY,LROUGH
          read (12,*) LPT,L_LSM,L_LSMbase
          read (12,*) LSCALAR,LAS,LNonNewt
          if (LNonNewt) then
          if (.not.LAS) then
          LAS=.TRUE.
          print*,"Density variable tracer activated"
          end if
          if (.not.LENERGY) then
          LENERGY=.TRUE.
          print*,"Energy equation activated"
          end if
          end if
          read (12,*) pl_ex
          read (12,*) Th,Tc,Tinit
          read (12,*)
          read (12,*) Tbc_w
          read (12,*) Tbc_e
          read (12,*) Tbc_s
          read (12,*) Tbc_n
          read (12,*) Tbc_b
          read (12,*) Tbc_t
          read (12,*)
          read (12,*) n_unstpt
          allocate(id_unst(n_unstpt),i_unst(n_unstpt) &
    ,j_unst(n_unstpt),k_unst(n_unstpt))
          do i=1,n_unstpt
              read (12,*)id_unst(i),i_unst(i),j_unst(i),k_unst(i)
          end do

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

      end subroutine read_control

#endif

 subroutine read_mdmap
          implicit none
          integer :: i,j,k,ib
          integer :: nptemp,myranktemp,nbtemp,cpu_no
          integer,allocatable,dimension(:) :: domtemp,buf_domindid
          character(len=80) :: dummyline

          if (myrank==0) then
          numfile=1001
          open (unit=numfile, file="output.dat")
          end if

          open (unit=12, file="input/mdmap.cin")

          read (12,*) nptemp  !number of processors

          if (nprocs /= nptemp) then
          print*, "=====ERROR====="
          print*, "number of cpus do not match map file"
          stop
          end if

          if (num_domains > 9999) then
          print*, "=====ERROR====="
          print*, "number of domains are exceeding", &
    " the limit in exchange subroutine"
          stop
          end if

          allocate (dom_id(num_domains),domtemp(num_domains))
          allocate (dom_indid(0:num_domains-1),dom_ad(0:num_domains-1))
          allocate (buf_domindid(0:num_domains-1))
          allocate (imbinblk(num_domains))  !Pablo

          dom_ad=-1
          read (12,*) dummyline

          nbpmax=0
          do i=1,nprocs
              read(12,*) myranktemp,nbtemp,domtemp(1:nbtemp)
              nbpmax=max(nbpmax,nbtemp)

              do ib=1,nbtemp
                  dom_ad(domtemp(ib))=myranktemp
              end do

              if(myranktemp==myrank) then
              nbp=nbtemp  !number of domains for this processor
              dom_id(1:nbp)=domtemp(1:nbp)

              do ib=0,num_domains-1
                  dom_indid(ib)=-1
              end do
              do ib=1,nbp
                  dom_indid(dom_id(ib))=ib
              end do

              end if
          end do

          read (12,*) dummyline
          close(12)


          buf_domindid = dom_indid

          call MPI_ALLREDUCE(buf_domindid,dom_indid,num_domains, &
    MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

          if(myrank==0) then
          do i=0,num_domains-1
              if(dom_ad(i)==-1) then
              print*,"unidentified domain in mdmap no:",i
              print*,"ERROR! check mdmap.cin"
              stop
              end if
              print*, "dom:",i," cpu:",dom_ad(i)," ib:",dom_indid(i)
              write(numfile,*) "dom:",i," cpu:",dom_ad(i)," ib:",dom_indid(i)
          end do
          end if

          allocate (dom(nbp))

          call MPI_BARRIER (MPI_COMM_WORLD,ierr)

          return
  end subroutine read_mdmap

end module io
