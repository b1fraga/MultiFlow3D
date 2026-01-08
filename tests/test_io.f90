module test_io
  use, intrinsic :: iso_fortran_env, only:  int8,dp => real64
  use testdrive, only : error_type, unittest_type, new_unittest, check
  use json_module
  implicit none
  private

  public :: collect_io

contains

  !> Collect all exported unit tests
  subroutine collect_io(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
         new_unittest("test_read_control_file", test_read_control_file)&
         ]
  end subroutine collect_io


  subroutine test_read_control_file(error)
    use io, only : read_control_file
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    character(kind=json_CK,len=:),allocatable :: Keyword,type_of_friction
    real(dp) :: dx,dy,dz,Ubulk,kinematic_visc,Pr,turb_Schmidt,beta
    real(dp) :: gx,gy,gz
    integer :: dens,convection_scheme,diffusion_scheme,differencing
    integer :: solver,multigrid_step,multigrid_iteration_scheme
    integer :: multigrid_maximum_iteration_per_time_step,restriction_iter,prolongation_iter
    real(dp) :: dt,safety_factor,eps,Friction_coefficient
    logical :: variable_dt,restart,reinitmean,LTRANSIENT
    integer :: sweeps,itime_end,n_out,results_output,niter
    integer :: nswp_1,nswp_2,nswp_3,nswp_4
    integer :: West_Boundary_Condition,East_Boundary_Condition
    integer :: South_Boundary_Condition,North_Boundary_Condition
    integer :: Bottom_Boundary_Condition,Top_Boundary_Condition
    logical :: save_inflow_data,time_averaging,SGS_model
    integer :: number_of_inlets, velocity_profile,Number_inlet_profiles
    real(dp) :: Turbulence_intensity,t_start_averaging1,t_start_averaging2
    real(dp) :: noise,Th,Tc
    integer :: SGS_model_value,LMR,normal_ghost_velocity_interpolation,pl_ex
    logical :: LIMB,LENERGY,LROUGH,LPT,LSM,L_LSMbase,LSCALAR,LActiveScalar,LNonNewt
    integer ::  West_Energy_BC,East_Energy_BC,South_Energy_BC,North_Energy_BC
    integer :: Bottom_Energy_BC,Top_Energy_BC
    integer :: num_of_time_series_points
    integer :: time_series_point_1,time_series_point_2,time_series_point_3
    integer :: time_series_point_4

    character (len=6) :: expected_character_keyword
    character (len=1) :: expected_character_type_of_friction
    integer :: expected_integer
    real(dp) :: expected_real
    logical :: expected_logical

    call read_control_file("test_io.json",Keyword,type_of_friction,&
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

    expected_character_keyword = "column"
    call check(error, Keyword, expected_character_keyword)
    if (allocated(error)) return

    expected_real = 0.1_dp
    call check(error, UBulk, expected_real)
    if (allocated(error)) return

    expected_real = 0.004_dp
    call check(error, dx, expected_real)
    if (allocated(error)) return

    expected_real = 0.03125_dp
    call check(error, dy, expected_real)
    if (allocated(error)) return

    expected_real = 0.03125_dp
    call check(error, dz, expected_real)
    if (allocated(error)) return

    expected_integer=1000
    call check(error, dens, expected_integer)
    if (allocated(error)) return

    expected_real=0.000001_dp
    call check(error, kinematic_visc, expected_real)
    if (allocated(error)) return

    expected_real=7.0_dp
    call check(error, Pr, expected_real)
    if (allocated(error)) return

    expected_real=0.6_dp
    call check(error, turb_Schmidt, expected_real)
    if (allocated(error)) return

    expected_real=0.000207_dp
    call check(error, beta, expected_real)
    if (allocated(error)) return

    expected_real=0.0_dp
    call check(error, gx, expected_real)
    if (allocated(error)) return

    expected_real=0.0_dp
    call check(error, gy, expected_real)
    if (allocated(error)) return

    expected_real=-9.81_dp
    call check(error, gz, expected_real)
    if (allocated(error)) return

    expected_integer = 3
    call check(error, convection_scheme, expected_integer)
    if (allocated(error)) return

    expected_integer = 3
    call check(error, diffusion_scheme, expected_integer)
    if (allocated(error)) return

    expected_integer = 1
    call check(error, differencing, expected_integer)
    if (allocated(error)) return

    expected_integer = 2
    call check(error, solver, expected_integer)
    if (allocated(error)) return

    expected_integer = 1
    call check(error, multigrid_step, expected_integer)
    if (allocated(error)) return

    expected_integer = 2
    call check(error, multigrid_iteration_scheme, expected_integer)
    if (allocated(error)) return

    expected_integer = 30
    call check(error, multigrid_maximum_iteration_per_time_step, expected_integer)
    if (allocated(error)) return

    expected_integer = 2
    call check(error, restriction_iter, expected_integer)
    if (allocated(error)) return

    expected_integer = 1
    call check(error, prolongation_iter, expected_integer)
    if (allocated(error)) return

    expected_real=0.005_dp
    call check(error, dt, expected_real)
    if (allocated(error)) return

    expected_logical=.false.
    call check(error, variable_dt, expected_logical)
    if (allocated(error)) return

    expected_integer=25
    call check(error, sweeps, expected_integer)
    if (allocated(error)) return

    expected_real=0.2_dp
    call check(error, safety_factor, expected_real)
    if (allocated(error)) return

    expected_integer=6000
    call check(error, itime_end, expected_integer)
    if (allocated(error)) return

    expected_logical=.false.
    call check(error, restart, expected_logical)
    if (allocated(error)) return

    expected_logical=.false.
    call check(error, reinitmean, expected_logical)
    if (allocated(error)) return

    expected_integer=6000
    call check(error, n_out, expected_integer)
    if (allocated(error)) return

    expected_logical=.true.
    call check(error, LTRANSIENT, expected_logical)
    if (allocated(error)) return

    expected_integer=100
    call check(error, results_output, expected_integer)
    if (allocated(error)) return

    expected_integer=20
    call check(error, niter, expected_integer)
    if (allocated(error)) return

    expected_real=0.00001_dp
    call check(error, eps, expected_real)
    if (allocated(error)) return

    expected_integer=5
    call check(error, nswp_1, expected_integer)
    if (allocated(error)) return

    expected_integer=5
    call check(error, nswp_2, expected_integer)
    if (allocated(error)) return

    expected_integer=5
    call check(error, nswp_3, expected_integer)
    if (allocated(error)) return

    expected_integer=20
    call check(error, nswp_4, expected_integer)
    if (allocated(error)) return

    expected_integer=4
    call check(error, West_Boundary_Condition, expected_integer)
    if (allocated(error)) return

    expected_integer=4
    call check(error, East_Boundary_Condition, expected_integer)
    if (allocated(error)) return

    expected_integer=4
    call check(error, South_Boundary_Condition, expected_integer)
    if (allocated(error)) return

    expected_integer=4
    call check(error, North_Boundary_Condition, expected_integer)
    if (allocated(error)) return

    expected_integer=4
    call check(error, Bottom_Boundary_Condition, expected_integer)
    if (allocated(error)) return

    expected_integer=4
    call check(error, Top_Boundary_Condition, expected_integer)
    if (allocated(error)) return

    expected_character_type_of_friction = "n"
    call check(error, type_of_friction, expected_character_type_of_friction)
    if (allocated(error)) return

    expected_real = 0.03_dp
    call check(error, Friction_coefficient, expected_real)
    if (allocated(error)) return

    expected_logical = .false.
    call check(error, save_inflow_data, expected_logical)
    if (allocated(error)) return

    expected_integer = 5000
    call check(error, number_of_inlets, expected_integer)
    if (allocated(error)) return

    expected_integer = 12
    call check(error, velocity_profile, expected_integer)
    if (allocated(error)) return

    expected_real = 0.1_dp
    call check(error, Turbulence_intensity, expected_real)
    if (allocated(error)) return

    expected_integer = 1000
    call check(error, Number_inlet_profiles, expected_integer)
    if (allocated(error)) return

    expected_logical = .false.
    call check(error, time_averaging, expected_logical)
    if (allocated(error)) return

    expected_real = 0.04_dp
    call check(error, t_start_averaging1, expected_real)
    if (allocated(error)) return

    expected_real = 0.05_dp
    call check(error, t_start_averaging2, expected_real)
    if (allocated(error)) return

    expected_real = 0.0_dp
    call check(error, noise, expected_real)
    if (allocated(error)) return

    expected_logical = .true.
    call check(error, SGS_model, expected_logical)
    if (allocated(error)) return

    expected_integer = 1
    call check(error, SGS_model_value, expected_integer)
    if (allocated(error)) return

    expected_integer = 2
    call check(error, LMR, expected_integer)
    if (allocated(error)) return

    expected_integer = 2
    call check(error, normal_ghost_velocity_interpolation, expected_integer)
    if (allocated(error)) return

    expected_logical = .false.
    call check(error, LIMB, expected_logical)
    if (allocated(error)) return

    expected_logical = .false.
    call check(error, LENERGY, expected_logical)
    if (allocated(error)) return

    expected_logical = .false.
    call check(error, LROUGH, expected_logical)
    if (allocated(error)) return

    expected_logical = .false.
    call check(error, LPT, expected_logical)
    if (allocated(error)) return

    expected_logical = .false.
    call check(error, LSM, expected_logical)
    if (allocated(error)) return

    expected_logical = .false.
    call check(error, L_LSMbase, expected_logical)
    if (allocated(error)) return

    expected_logical = .true.
    call check(error, LSCALAR, expected_logical)
    if (allocated(error)) return

    expected_logical = .true.
    call check(error, LActiveScalar, expected_logical)
    if (allocated(error)) return

    expected_logical = .true.
    call check(error, LNonNewt, expected_logical)
    if (allocated(error)) return

    expected_integer = 2
    call check(error, pl_ex, expected_integer)
    if (allocated(error)) return

    expected_real = 301.11_dp
    call check(error, Th, expected_real)
    if (allocated(error)) return

    expected_real = 298.57_dp
    call check(error, Tc, expected_real)
    if (allocated(error)) return

    expected_integer = 7
    call check(error, West_Energy_BC, expected_integer)
    if (allocated(error)) return

    expected_integer = 7
    call check(error, East_Energy_BC, expected_integer)
    if (allocated(error)) return

    expected_integer = 7
    call check(error, South_Energy_BC, expected_integer)
    if (allocated(error)) return

    expected_integer = 7
    call check(error, North_Energy_BC, expected_integer)
    if (allocated(error)) return

    expected_integer = 7
    call check(error, Bottom_Energy_BC, expected_integer)
    if (allocated(error)) return

    expected_integer = 7
    call check(error, Top_Energy_BC, expected_integer)
    if (allocated(error)) return

    expected_integer = 0
    call check(error, num_of_time_series_points, expected_integer)
    if (allocated(error)) return

    expected_integer = 0
    call check(error, time_series_point_1, expected_integer)
    if (allocated(error)) return

    expected_integer = 12
    call check(error, time_series_point_2, expected_integer)
    if (allocated(error)) return

    expected_integer = 22
    call check(error, time_series_point_3, expected_integer)
    if (allocated(error)) return

    expected_integer = 22
    call check(error, time_series_point_4, expected_integer)
    if (allocated(error)) return

  end subroutine test_read_control_file



end module test_io
