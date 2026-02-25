program tester
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type
#if USE_JSON == 1
  use test_json_io, only : collect_json
  use test_io, only : collect_io
#endif
#if USE_HDF5 == 1
  use test_hdf5_io, only: collect_hdf5
#endif
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'
#if USE_HDF5 == 1
  type(testsuite_type) :: hdf5_tests
  integer :: unit
#endif
  
  stat = 0

#if USE_JSON == 1
  testsuites = [ &
       new_testsuite("json", collect_json),&
       new_testsuite("io", collect_io)&
    ]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do
#endif

#if USE_HDF5 == 1
    hdf5_tests = new_testsuite("hdf5", collect_hdf5)
    call run_testsuite(hdf5_tests%collect, error_unit, stat)
    open(file="test_output.h5", newunit=unit)
    close(unit, status="delete")
#endif
    
  if (stat > 0) then
    write(error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
    error stop
  end if

end program tester
