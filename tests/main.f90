program tester
  use, intrinsic :: iso_fortran_env, only : error_unit
  use multiflow3d_mpi, only: init_parallelisation, end_parallelisation, myrank
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type
#if USE_JSON == 1
  use test_json_io, only : collect_json
#endif
  use test_io, only : collect_io
#if USE_HDF5 == 1
  use test_hdf5_io, only: collect_hdf5
#endif
  use test_MPI_pt, only: collect_MPI_pt
  implicit none
  integer :: stat, is, ierr
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'
#if USE_HDF5 == 1
  type(testsuite_type) :: hdf5_tests
  integer :: unit
#endif

  call init_parallelisation()

  stat = 0

#if USE_JSON == 1
  testsuites = [ &
       new_testsuite("json", collect_json),&
       new_testsuite("io", collect_io),&
       new_testsuite("MPI_pt",collect_MPI_pt)&
    ]
#else
  testsuites = [ &
       new_testsuite("io", collect_io),&
       new_testsuite("MPI_pt",collect_MPI_pt)&
    ]  
#endif
  
  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do


#if USE_HDF5 == 1
    if(myrank == 0) then
       hdf5_tests = new_testsuite("hdf5", collect_hdf5)
       call run_testsuite(hdf5_tests%collect, error_unit, stat)
       open(file="test_output.h5", newunit=unit)
       close(unit, status="delete")
    endif
#endif
    
  if (stat > 0) then
    write(error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
    error stop
  end if

  call end_parallelisation()

end program tester
