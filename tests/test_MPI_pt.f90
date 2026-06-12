module test_MPI_pt
  use iso_fortran_env, only : dp => real64
  use mpi_f08

  use testdrive, only : new_unittest, unittest_type, error_type, check

  use multiflow3d_MPI_pt
  use multiflow3d_mpi, only : Myrank, nprocs, ierr

  implicit none
  private
  public :: collect_MPI_pt

contains

!==========================================================
subroutine collect_MPI_pt(tests)

  type(unittest_type), allocatable, intent(out) :: tests(:)

  tests = [ &
       new_unittest("MPI_pt distributes particles to 4 ranks", test_MPI_pt_four_ranks) &
       ]

end subroutine collect_MPI_pt

!==========================================================
subroutine test_MPI_pt_four_ranks(error)

  type(error_type), allocatable, intent(out) :: error

  logical :: Lcol
  integer :: np, np_loc, npg_loc
  integer, allocatable :: ptsinproc(:), ptsinproc_g(:), id(:)

  real(dp), allocatable :: xp_pt(:), yp_pt(:), zp_pt(:)
  real(dp), allocatable :: uop_pt(:), vop_pt(:), wop_pt(:)
  real(dp), allocatable :: dp_pt(:), rho_pt(:)

  real(dp), allocatable :: xp_loc(:), yp_loc(:), zp_loc(:)
  real(dp), allocatable :: uop_loc(:), vop_loc(:), wop_loc(:)
  real(dp), allocatable :: dp_loc(:), rhop_loc(:)

  real(dp), allocatable :: xpg_loc(:), ypg_loc(:), zpg_loc(:)
  real(dp), allocatable :: uopg_loc(:), vopg_loc(:), wopg_loc(:)
  real(dp), allocatable :: dpg_loc(:), rhopg_loc(:)

  integer :: dom_ad(4)

  allocate(ptsinproc(4), ptsinproc_g(4))

  np = 4
  np_loc = 0
  npg_loc = 0
  Lcol = .false.

  if (Myrank == 0) then
     allocate(xp_pt(np), yp_pt(np), zp_pt(np))
     allocate(uop_pt(np), vop_pt(np), wop_pt(np))
     allocate(dp_pt(np), rho_pt(np))

     xp_pt = [0.5_dp, 1.5_dp, 2.5_dp, 3.5_dp]
     yp_pt = 0.5_dp
     zp_pt = 0.5_dp
     uop_pt = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp]
     dp_pt  = 1.0_dp
     rho_pt = 1000.0_dp
  else
     allocate(xp_pt(0), yp_pt(0), zp_pt(0))
     allocate(uop_pt(0), vop_pt(0), wop_pt(0))
     allocate(dp_pt(0), rho_pt(0))
  end if

  dom_ad = [0,1,2,3]

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  call MPI_pt(Lcol, np_loc, npg_loc, &
       ptsinproc, ptsinproc_g, &
       xp_pt, yp_pt, zp_pt, &
       uop_pt, vop_pt, wop_pt, &
       dp_pt, rho_pt, &
       xp_loc, yp_loc, zp_loc, &
       uop_loc, vop_loc, wop_loc, &
       dp_loc, rhop_loc, &
       id, &
       xpg_loc, ypg_loc, zpg_loc, &
       uopg_loc, vopg_loc, wopg_loc, &
       dpg_loc, rhopg_loc, &
       0.1_dp, 0.1_dp, 0.1_dp, &
       np, &
       4.0_dp, 0.0_dp, &
       1.0_dp, 0.0_dp, &
       1.0_dp, 0.0_dp, &
       4, 1, 1, &
       dom_ad)

  call check(error, np_loc , ptsinproc(myrank+1))
  if (allocated(error)) return

end subroutine test_MPI_pt_four_ranks



end module test_MPI_pt
