#if USE_JSON == 1
module json_io
  use json_module
  use, intrinsic :: iso_fortran_env, only: int8, dp => real64

  implicit none
  private
  public :: json_read

  interface json_read
     module procedure json_read_integer
     module procedure json_read_character
     module procedure json_read_integer_array
     module procedure json_read_real
     module procedure json_read_logical
  end interface

contains

  subroutine json_read_integer(input_file,value_to_return,output)
    !! Reads an integer from a json file, which has the key
    !! given by value_to_return
    type(json_file) :: json
    character(len=*), intent(in) :: input_file, value_to_return
    logical :: found
    integer,intent(out) :: output

    call json%initialize()

    call json%load(filename = input_file)

    call json%get(value_to_return, output, found)
    if ( .not. found ) error stop "Value not found in json file"

    call json%destroy()
    if (json%failed()) error stop "Failed to close json file"

  end subroutine json_read_integer


  subroutine json_read_real(input_file,value_to_return,output)
    !! Reads a real from a json file, which has the key
    !! given by value_to_return
    type(json_file) :: json
    character(len=*), intent(in) :: input_file, value_to_return
    logical :: found
    real(dp) ,intent(out) :: output

    call json%initialize()

    call json%load(filename = input_file)

    call json%get(value_to_return, output, found)
    if ( .not. found ) error stop "Value not found in json file"

    call json%destroy()
    if (json%failed()) error stop "Failed to close json file"

  end subroutine json_read_real

  subroutine json_read_integer_array(input_file,value_to_return,output)
    !! Reads a real from a json file, which has the key
    !! given by value_to_return
    type(json_file) :: json
    character(len=*), intent(in) :: input_file, value_to_return
    logical :: found
    integer , allocatable, intent(out) :: output(:)

    call json%initialize()

    call json%load(filename = input_file)

    call json%get(value_to_return, output, found)
    if ( .not. found ) error stop "Value not found in json file"

    call json%destroy()
    if (json%failed()) error stop "Failed to close json file"

  end subroutine json_read_integer_array

  subroutine json_read_logical(input_file,value_to_return,output)
    !! Reads a logical from a json file, which has the key
    !! given by value_to_return
    type(json_file) :: json
    character(len=*), intent(in) :: input_file, value_to_return
    logical :: found
    logical,intent(out) :: output

    call json%initialize()

    call json%load(filename = input_file)

    call json%get(value_to_return, output, found)
    if ( .not. found ) error stop "Value not found in json file"

    call json%destroy()
    if (json%failed()) error stop "Failed to close json file"

  end subroutine json_read_logical

  subroutine json_read_character(input_file,value_to_return,output)
    !! Reads a character string from a json file, which has the key
    !! given by value_to_return
    type(json_file) :: json
    character(len=*), intent(in) :: input_file, value_to_return
    logical :: found
    character(kind=json_CK,len=:), allocatable :: dummy_output
    character(len=80),intent(inout) :: output

    call json%initialize()

    call json%load(filename = input_file)
    call json%get(value_to_return, dummy_output, found)
    if ( .not. found ) error stop "Value not found in json file"
    output=trim(dummy_output)

    call json%destroy()
    if (json%failed()) error stop "Failed to close json file"

  end subroutine json_read_character

end module json_io
#endif
