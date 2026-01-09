module io
  use json_io, only: json_read
  use json_module
  use, intrinsic :: iso_fortran_env, only: dp => real64

  implicit none
  private
  public :: read_control_file


contains

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
       noise,Th,Tc,SGS_model_value,LMR,normal_ghost_velocity_interpolation,pl_ex,&
       LIMB,LENERGY,LROUGH,LPT,LSM,L_LSMbase,LSCALAR,LActiveScalar,LNonNewt,&
       West_Energy_BC,East_Energy_BC,South_Energy_BC,North_Energy_BC,&
       Bottom_Energy_BC,Top_Energy_BC,num_of_time_series_points,&
       time_series_point_1,time_series_point_2,time_series_point_3,time_series_point_4)
    !! Reads control.json assigning the variables associated output to variables of the subroutine
    character(len=*), intent(in) :: input_file
    character(kind=json_CK,len=:),allocatable, intent(out) :: Keyword,type_of_friction
    real(dp), intent(out) :: dx,dy,dz,Ubulk,kinematic_visc,Pr,turb_Schmidt,beta
    real(dp), intent(out) :: gx,gy,gz
    integer, intent(out) :: dens,convection_scheme,diffusion_scheme,differencing
    integer, intent(out) :: solver,multigrid_step,multigrid_iteration_scheme
    integer, intent(out) :: multigrid_maximum_iteration_per_time_step,restriction_iter,prolongation_iter
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
    real(dp), intent(out) :: noise,Th,Tc
    integer, intent(out) :: SGS_model_value,LMR,normal_ghost_velocity_interpolation,pl_ex
    logical, intent(out) :: LIMB,LENERGY,LROUGH,LPT,LSM,L_LSMbase,LSCALAR,LActiveScalar,LNonNewt
    integer, intent(out) ::  West_Energy_BC,East_Energy_BC,South_Energy_BC,North_Energy_BC
    integer, intent(out) :: Bottom_Energy_BC,Top_Energy_BC
    integer, intent(out) :: num_of_time_series_points
    integer, intent(out) :: time_series_point_1,time_series_point_2,time_series_point_3
    integer, intent(out) :: time_series_point_4

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
    call json_read(input_file,"normal ghost velocity interpolation",&
         normal_ghost_velocity_interpolation)
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
    call json_read(input_file,"time_series_point_1",time_series_point_1)
    call json_read(input_file,"time_series_point_2",time_series_point_2)
    call json_read(input_file,"time_series_point_3",time_series_point_3)
    call json_read(input_file,"time_series_point_4",time_series_point_4)
  end subroutine read_control_file

end module io
