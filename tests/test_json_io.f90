#if USE_JSON == 1
module test_json_io
  use, intrinsic :: iso_fortran_env, only:  int8,real64
  use testdrive, only : error_type, unittest_type, new_unittest, check
  use json_module
  implicit none
  private

  public :: collect_json

contains

  !> Collect all exported unit tests
  subroutine collect_json(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
         new_unittest("json_read_integer", test_json_read_integer),&
         new_unittest("json_read_real", test_json_read_real),&
         new_unittest("json_read_logical", test_json_read_logical),&
         new_unittest("json_read_character", test_json_read_character)&
         ]
  end subroutine collect_json


  subroutine test_json_read_integer(error)
    use json_io, only : json_read
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    integer :: expected
    integer :: output

    call json_read("test_file.json","value1",output)
    expected = 42
    call check(error, output, expected)
    if (allocated(error)) return
  end subroutine test_json_read_integer

  subroutine test_json_read_real(error)
    use json_io, only : json_read
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    real(real64) :: expected
    real(real64) :: output

    call json_read("test_file.json","value2",output)
    expected = 3.14_real64
    call check(error, output, expected)
    if (allocated(error)) return
  end subroutine test_json_read_real

  subroutine test_json_read_logical(error)
    use json_io, only : json_read
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    logical :: expected
    logical :: output

    call json_read("test_file.json","value3",output)
    expected = .true.
    call check(error, output, expected)
    if (allocated(error)) return
    expected = .false.
    call json_read("test_file.json","value4",output)
    if (allocated(error)) return
  end subroutine test_json_read_logical

  subroutine test_json_read_character(error)
    use json_io, only : json_read
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    character(len=5) :: expected
    character(len=80)  :: output

    call json_read("test_file.json","value5",output)
    expected = "Hello"
    call check(error, output, expected)
    if (allocated(error)) return

  end subroutine test_json_read_character

end module test_json_io
#endif
