module io
#if USE_JSON == 1
  use json_io, only: json_read
#endif
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use multidata, only: multidom
  use multiflow3d_mpi, only: ierr, mpi_comm_world, mpi_integer, mpi_max, myrank, nprocs, myrank

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
  subroutine read_control(eps,n_out,niter,sweeps,tbc_w,nswp,th,tc,ti_sem, &
                          time_averaging,tinit,tsteps_pt,ubulk,UPROF_SEM, &
                          t_start_averaging1,t_start_averaging2,Tbc_e,Tbc_s, &
                          Tbc_n,Tbc_b,Tbc_t,reinitmean,rrey,safety_factor, &
                          save_inflow,Sc_t,SGS,sgs_model,solver,n_unstpt, &
                          ngrid_input,noise,np,pl_ex,pr,pressureforce,re, &
                          read_inflow,LMR,LNonNewt,LPT,LRESTART,LROUGH, &
                          LSCALAR,LTRANSIENT,maxcy,mg_itrsch,ITMAX_PI, &
                          ITMAX_SEM,keyword,L_dt,L_LSM,L_LSMbase,L_n, &
                          LAS,LENERGY,LIMB,differencing,dt,fric,g_dx,g_dy,&
                          g_dz,gx,gy,gz,iproln,irestr,itime_end,bc_w, &
                          bc_e,bc_s,bc_n,bc_b,bc_t,beta,conv_sch,dens,diff_sch, &
                          id_unst,i_unst,j_unst,k_unst,nbp,dom)
    integer :: ib,i
    real(dp),intent(out) :: eps,th,tc,ti_sem,tinit,ubulk,t_start_averaging1, &
                            t_start_averaging2, rrey, safety_factor,Sc_t,noise, &
                            pr,re,dt,fric,g_dx,g_dy,g_dz,gx,gy,gz, beta,dens
    integer, intent(out) :: n_out,niter,sweeps,nswp(4),tsteps_pt, &
                            UPROF_SEM, ngrid_input,np,pl_ex,itime_end, &
                            bc_w,bc_e,bc_s,bc_n,bc_b,bc_t,conv_sch,diff_sch
    integer, intent(in) :: nbp
    integer, intent(out) :: Tbc_w,Tbc_e,Tbc_s,Tbc_n,Tbc_b,Tbc_t,sgs_model, &
                            solver,n_unstpt, LMR, maxcy,mg_itrsch,ITMAX_PI,ITMAX_SEM, &
                            differencing,iproln,irestr
    logical, intent(out) :: time_averaging, reinitmean, save_inflow,SGS, &
                            pressureforce,read_inflow,LNonNewt,LPT,LRESTART,LROUGH,LSCALAR, &
                            LTRANSIENT,L_dt,L_LSM,L_LSMbase,LAS,LENERGY,LIMB
    character(len=80), intent(out) :: keyword,L_n
    integer,allocatable,dimension(:),intent(out) :: id_unst
    integer,allocatable,dimension(:),intent(out) :: i_unst,j_unst,k_unst
    type (multidom), pointer, dimension(:), intent(inout) :: dom
    integer :: file_unit

    open (newunit=file_unit, file="input/control.cin")
    !------DOMAIN SIZE AND DISCRETIZATION -------------------------------------
    read (file_unit,*)
    read (file_unit,*) keyword,ubulk
    read (file_unit,*) g_dx,g_dy,g_dz
    read (file_unit,*) dens,rrey,Pr,Sc_t,beta
    Re=1.0d0/rrey
    read (file_unit,*) gx,gy,gz
    read (file_unit,*) conv_sch
    read (file_unit,*) diff_sch
    read (file_unit,*) differencing
    read (file_unit,*) solver
    read (file_unit,*) ngrid_input,mg_itrsch
    read (file_unit,*) maxcy,irestr,iproln
    read (file_unit,*) dt,L_dt,sweeps,safety_factor
    read (file_unit,*) itime_end,LRESTART,reinitmean,n_out
    read (file_unit,*) LTRANSIENT,tsteps_pt
    read (file_unit,*) niter,eps,nswp(1),nswp(2),nswp(3),nswp(4)
    read (file_unit,*)
    read (file_unit,*) bc_w
    read (file_unit,*) bc_e
    read (file_unit,*) bc_s
    read (file_unit,*) bc_n
    read (file_unit,*) bc_b
    read (file_unit,*) bc_t
    if (bc_w==5) pressureforce=.true.
    read (file_unit,*) L_n,fric
    read (file_unit,*) save_inflow,ITMAX_PI
    read (file_unit,*)
    read (file_unit,*) UPROF_SEM         !Pablo 15/12/2015
    read (file_unit,*) TI_SEM
    read (file_unit,*) ITMAX_SEM
    read (file_unit,*)
    read (file_unit,*) time_averaging,t_start_averaging1, &
                t_start_averaging2,noise
    read (file_unit,*) SGS,sgs_model
    read (file_unit,*) LMR
    read (file_unit,*) LIMB,LENERGY,LROUGH
    read (file_unit,*) LPT,L_LSM,L_LSMbase
    read (file_unit,*) LSCALAR,LAS,LNonNewt
    if(LNonNewt) then
      if(.not.LAS) then
        LAS=.true.
        print*,"Density variable tracer activated"
      end if
      if(.not.LENERGY) then
        LENERGY=.true.
        print*,"Energy equation activated"
      end if
    end if
    read (file_unit,*) pl_ex
    read (file_unit,*) Th,Tc,Tinit
    read (file_unit,*)
    read (file_unit,*) Tbc_w
    read (file_unit,*) Tbc_e
    read (file_unit,*) Tbc_s
    read (file_unit,*) Tbc_n
    read (file_unit,*) Tbc_b
    read (file_unit,*) Tbc_t
    read (file_unit,*)
    read (file_unit,*) n_unstpt
    allocate(id_unst(n_unstpt),i_unst(n_unstpt) &
             ,j_unst(n_unstpt),k_unst(n_unstpt))
    do i=1,n_unstpt
       read (file_unit,*)id_unst(i),i_unst(i),j_unst(i),k_unst(i)
    end do

    close(file_unit)

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

    if(L_LSM .and. solver==1) then
      if(myrank==0) then
        print*,"Error: SIP solver not presently compatible with LSM"
      end if
      stop
    end if

    if(L_LSM .and. differencing/=3) then
      if(myrank==0) then
        print*,"Error: WENO differencing must be used with LSM"
      end if
      stop
    end if

  end subroutine read_control
#endif
  subroutine read_mdmap(numfile,dom_id, dom_indid,dom_ad,imbinblk,dom,nbp,num_domains)
    integer :: i,ib,file_unit
    integer :: nptemp,myranktemp,nbtemp
    integer,allocatable,dimension(:) :: domtemp,buf_domindid
    integer, intent(out) :: numfile,nbp,num_domains
    integer,allocatable,dimension(:),intent(out) :: dom_id, dom_indid, dom_ad
    integer,allocatable,dimension(:),intent(out) :: imbinblk
    type (multidom), pointer, dimension(:),intent(out) :: dom

    if(myrank==0) then
      open(newunit=numfile, file="output.dat",status="replace", action="write")
    end if

    open(newunit=file_unit, file="input/mdmap.cin", status="old", action="read")

    read(file_unit,*) num_domains  !number of domains
    read(file_unit,*) nptemp  !number of processors

    if(nprocs /= nptemp) then
      error stop "number of cpus do not match map file"
    end if

    if(num_domains > 9999) then
      error stop "number of domains are exceeding the limit in exchange subroutine"
    end if

    allocate(dom_id(num_domains),domtemp(num_domains))
    allocate(dom_indid(0:num_domains-1),dom_ad(0:num_domains-1))
    allocate(buf_domindid(0:num_domains-1))
    allocate(imbinblk(num_domains))  !Pablo

    dom_ad=-1
    read (file_unit,*)

    do i=1,nprocs
       read(file_unit,*) myranktemp,nbtemp,domtemp(1:nbtemp)

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

    close(file_unit)


    buf_domindid = dom_indid

    call MPI_ALLREDUCE(buf_domindid,dom_indid,num_domains, &
                       MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

    if(myrank==0) then
      do i=0,num_domains-1
         if (dom_ad(i) == -1) then
            write(*,*) "unidentified domain in mdmap no:", i
            error stop "check mdmap.cin"
         end if
         print*, "dom:",i," cpu:",dom_ad(i)," ib:",dom_indid(i)
         write(numfile,*) "dom:",i," cpu:",dom_ad(i)," ib:",dom_indid(i)
      end do
    end if

    allocate (dom(nbp))

    call MPI_BARRIER (MPI_COMM_WORLD,ierr)

  end subroutine read_mdmap

end module io
