!> Options utilities.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune, 
!! Toni Collis, David Scott, Peter Spelt, Mike Jackson.
!! @version $Revision: 204 $
!! @copyright (c) 2013-2015, Prashant Valluri, Lenn>on O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt, The University of
!! Edinburgh, all rights reserved.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module options_utils

  implicit none

  !> Polymorphic subroutine to get options.
  interface get_typed_option
     module procedure get_logical_option
     module procedure get_integer_option
     module procedure get_double_option
     module procedure get_string_option
  end interface get_typed_option

  !> Polymorphic subroutine to put options.
  interface put_typed_option
     module procedure put_logical_option
     module procedure put_integer_option
     module procedure put_double_option
     module procedure put_string_option
  end interface put_typed_option

  !> Polymorphic function for invalid option messages.
  interface invalid_option_message
     module procedure invalid_string_option_message
     module procedure invalid_integer_option_message
     module procedure invalid_double_option_message
     module procedure invalid_logical_option_message
  end interface

  !> Option not found error.
  integer :: option_not_found_err
  parameter(option_not_found_err=100)
  !> Option format error in file - option not of form <tt>NAME VALUE</tt>.
  integer :: option_format_err
  parameter(option_format_err=101)
  !> Option invalid value error.
  integer :: option_invalid_err
  parameter(option_invalid_err=102)
  !> Option file not found error.
  integer :: option_file_not_found_err
  parameter(option_file_not_found_err=103)

contains

  !> Is an option with a given name available?
  logical function has_option(options,num_options,name)
    character(len=*),dimension(:,:),intent(in) :: options !< Name-value options.
    integer,         intent(in)  :: num_options !< Number of options.
    character(len=*),intent(in)  :: name        !< Option name.

    integer :: i

    has_option=.false.
    do i=1,num_options
       if (options(i,1)==name) then
          has_option=.true.
          exit 
       end if
    end do
  end function has_option


  !> Get value of an option with a given name.
  subroutine get_option(options,num_options,name,value,ierr)
    character(len=*),dimension(:,:),intent(in) :: options !< Name-value options.
    integer,         intent(in)  :: num_options !< Number of options.
    character(len=*),intent(in)  :: name        !< Option name.
    character(len=*),intent(out) :: value       !< Option value.
    integer,         intent(out) :: ierr        !< 0 (OK) or <tt>option_not_found_err</tt>.

    integer :: i

    ierr = option_not_found_err
    do i=1,num_options
       if (options(i,1)==name) then
          value=options(i,2)
          ierr = 0
          exit 
       end if
    end do
  end subroutine get_option


  !> Get a logical option of a given name.
  subroutine get_logical_option(options,num_options,name,value,ierr)
    character(len=*),dimension(:,:),intent(in) :: options !< Name-value options.
    integer,         intent(in)  :: num_options !< Number of options.
    character(len=*),intent(in)  :: name        !< Option name.
    logical,         intent(out) :: value       !< Option value.
    integer,         intent(out) :: ierr        !< 0 (OK) or <tt>option_not_found_err</tt>.

    character(len=60) :: value_str

    call get_option(options,num_options,name,value_str,ierr)
    if (ierr == 0) read(value_str,*) value
  end subroutine get_logical_option


  !> Get an integer option of a given name.
  subroutine get_integer_option(options,num_options,name,value,ierr)
    character(len=*),dimension(:,:),intent(in) :: options !< Name-value options.
    integer,         intent(in)  :: num_options !< Number of options.
    character(len=*),intent(in)  :: name        !< Option name.
    integer,         intent(out) :: value       !< Option value.
    integer,         intent(out) :: ierr        !< 0 (OK) or <tt>option_not_found_err</tt>.

    character(len=60) :: value_str

    call get_option(options,num_options,name,value_str,ierr)
    if (ierr == 0) read(value_str,*) value
  end subroutine get_integer_option


  !> Get a double option of a given name.
  subroutine get_double_option(options,num_options,name,value,ierr)
    character(len=*),dimension(:,:),intent(in) :: options !< Name-value options.
    integer,         intent(in)  :: num_options !< Number of options.
    character(len=*),intent(in)  :: name        !< Option name.
    double precision,intent(out) :: value       !< Option value.
    integer,         intent(out) :: ierr        !< 0 (OK) or <tt>option_not_found_err</tt>.

    character(len=60) :: value_str

    call get_option(options,num_options,name,value_str,ierr)
    if (ierr == 0) read(value_str,*) value
  end subroutine get_double_option


  !> Get an string option of a given name.
  subroutine get_string_option(options,num_options,name,value,ierr)
    character(len=*),dimension(:,:),intent(in) :: options !< Name-value options.
    integer,         intent(in)  :: num_options !< Number of options.
    character(len=*),intent(in)  :: name        !< Option name.
    character(len=*),intent(out) :: value       !< Option value.
    integer,         intent(out) :: ierr        !< 0 (OK) or <tt>option_not_found_err</tt>.

    call get_option(options,num_options,name,value,ierr)
  end subroutine get_string_option


  !> Put an option name-value pair into an array of options.
  !! If an option with the name already exists then the value is
  !! updated to be the new value.
  !! If, not then it is added. It is assumed that the options array
  !! has enough space to record the new option. <tt>num_options</tt>
  !! is incremented by 1. 
  subroutine put_option(options,num_options,name,value)
    character(len=*),dimension(:,:),allocatable,intent(inout) :: options !< Name-value options.
    integer,         intent(inout) :: num_options !< Number of options.
    character(len=*),intent(in)    :: name        !< Option name. 
    character(len=*),intent(in)    :: value       !< Option value. 

    integer :: i

    do i=1,num_options
       if (options(i,1)==name) then
          options(i,2) = value
          return
       end if
    end do
    num_options=num_options+1
    options(num_options,1) = name
    options(num_options,2) = value
  end subroutine put_option


  !> Put a name-value into an array of options.
  !! The value is converted into a string using <tt>L1</tt> format.
  subroutine put_logical_option(options,num_options,name,value)
    character(len=*),dimension(:,:),allocatable,intent(inout) :: options !< Name-value options.
    integer,         intent(inout) :: num_options !< Number of options.
    character(len=*),intent(in)    :: name        !< Option name. 
    logical,         intent(in)    :: value       !< Option value. 

    character(len=256) :: value_str
    write(value_str,'(L1)') value
    call put_option(options,num_options,name,value_str)
  end subroutine put_logical_option


  !> Put a name-value into an array of options.
  !! The value is converted into a string using <tt>I10</tt> format.
  subroutine put_integer_option(options,num_options,name,value)
    character(len=*),dimension(:,:),allocatable,intent(inout) :: options !< Name-value options.
    integer,         intent(inout) :: num_options !< Number of options.
    character(len=*),intent(in)    :: name        !< Option name. 
    integer,         intent(in)    :: value       !< Option value. 

    character(len=256) :: value_str
    write(value_str,'(I10)') value
    call put_option(options,num_options,name,value_str)
  end subroutine put_integer_option


  !> Put a name-value into an array of options.
  !! The value is converted into a string using <tt>F20.10</tt> format.
  subroutine put_double_option(options,num_options,name,value)
    character(len=*),dimension(:,:),allocatable,intent(inout) :: options !< Name-value options.
    integer,         intent(inout) :: num_options !< Number of options.
    character(len=*),intent(in)    :: name        !< Option name. 
    double precision,intent(in)    :: value       !< Option value. 

    character(len=256) :: value_str
    write(value_str,'(F20.10)') value
    call put_option(options,num_options,name,value_str)
  end subroutine put_double_option


  !> Put a name-value into an array of options.
  subroutine put_string_option(options,num_options,name,value)
    character(len=*),dimension(:,:),allocatable,intent(inout) :: options !< Name-value options.
    integer,         intent(inout) :: num_options !< Number of options.
    character(len=*),intent(in)    :: name        !< Option name. 
    character(len=*),intent(in)    :: value       !< Option value. 

    call put_option(options,num_options,name,value)
  end subroutine put_string_option


  !> Print options to current output stream.
  subroutine print_options(options,num_options)
    character(len=*),dimension(:,:),intent(in) :: options !< Name-value options.
    integer,intent(in) :: num_options !< Number of options.

    integer :: i

    do i=1,num_options
       write(*,'(A,A1,A)') trim(options(i,1)),' ',trim(options(i,2))
    end do
  end subroutine print_options


  !> Load options from a file.
  !! Each line of the file is assumed to be one of:
  !! - A comment, prefixed by <tt>#</tt>
  !! - A blank line
  !! - An option value pair of form <tt>NAME VALUE</tt>.
  !! It is assumed that lines are <=80 characters, names are <=20 characters
  !! and values are <=60 characters.
  !! If there is an error then <tt>ierr</tt> will be set to 0 if OK, 
  !! <tt>option_file_not_found_err</tt> if there is no such file, a 
  !! platform-specific value or to <tt>option_format_err</tt> if an
  !! line that is not a comment but not an option value pair is found.
  subroutine load_options(options_file,options,num_options,ierr)
    character(len=*),intent(in) :: options_file !< File name.
    character(len=*),dimension(:,:),allocatable,intent(out) :: options !< Name-value options.
    integer,intent(out) :: num_options !< Number of options.
    integer,intent(out) :: ierr        !< Error code.

    integer :: file_id,status,split
    parameter(file_id=135)
    character(1) :: separator
    parameter(separator=' ')
    character(80) :: line
    character(20) :: name
    character(60) :: value
    logical       :: file_exists

    inquire(file=options_file, exist=file_exists)
    if (.not. file_exists) then
       ierr = option_file_not_found_err
       return
    end if
    ! Count number of lines in file.
    num_options=0
    open(unit=file_id,file=options_file,action='READ',status='OLD')
    do
       read(file_id,*,iostat=status) line
       num_options=num_options+1
       if (status>0) then
          ierr = status
          close(file_id)
          return
       end if
       if (status<0) then
          exit
       end if
    end do
    close(file_id)

    ! Populate options.
    allocate(options(num_options,2))
    num_options=0
    open(unit=file_id,file=options_file,action='READ',status='OLD')
    do
       read(file_id,'(A)',iostat=status) line
       line=adjustl(line)
       line=trim(line)
       if (status > 0) then
          ierr=status
          close(file_id)
          return
       end if
       if (status<0) then
          exit
       end if
       if (line(1:1)=='#') then
          ! Ignore comment.
          cycle
       end if
       if (line(1:1)==' ') then
          ! Ignore blank linex.
          cycle
       end if
       if (len(line)==0) then
          ! Ignore blank line.
          cycle
       end if
       split=index(line,separator) 
       if (split<len(trim(line))) then
          name=trim(line(:split-1))
          value=trim(line(split+1:))
          call put_option(options,num_options,name,value)
       else
          ierr=option_format_err
          close(file_id)
          return
       end if
    end do
    ierr=0
    close(file_id)
  end subroutine load_options


  !> Save options to a file.
  !! Header lines can be provided. If so these are added first with the
  !! prefix <tt># </tt>.
  subroutine save_options(filename,options,num_options,hdr,num_hdr_lines)
    character(len=*),intent(in) :: filename !< File name.
    character(len=*),dimension(:,:),intent(in)  :: options !< Name-value options.
    integer,         intent(in)  :: num_options !< Number of options.
    character(len=*),dimension(:),intent(in)  :: hdr !< Header lines.
    integer,         intent(in)  :: num_hdr_lines !< Number of options.

    integer :: file_id
    parameter(file_id=135)
    integer :: i

    open(unit=file_id,file=filename,action='WRITE',status='REPLACE')
    do i=1,num_hdr_lines
      write(file_id,'(A2,A)') '# ',trim(hdr(i))
    end do
    close(file_id)
    call append_options(filename,options,num_options)
  end subroutine save_options


  !> Append options to a file.
  !! It is assumed that the file exists.
  subroutine append_options(filename,options,num_options)
    character(len=*),intent(in) :: filename !< File name.
    character(len=*),dimension(:,:),intent(in)  :: options !< Name-value options.
    integer,         intent(in)  :: num_options !< Number of options.

    integer :: file_id
    parameter(file_id=135)
    integer :: i

    open(unit=file_id,file=filename,action='WRITE',status='OLD',position='append')
    do i=1,num_options
      write(file_id,'(A,A1,A)') trim(options(i,1)),' ',trim(options(i,2))
    end do
    close(file_id)
  end subroutine append_options


  !> Gets options and values from the command-line. 
  !! Command-line arguments are assumed to be of form:
  !! - <tt>-n value</tt>
  !! - <tt>-n</tt>, in which case the value is an empty string.
  subroutine get_command_line_options(options,num_options)
    character(len=*),dimension(:,:),allocatable,intent(out) :: options !< Name-value options.
    integer,intent(inout) :: num_options !< Number of options.

    integer       :: num_arguments,i
    character(60) :: argument,name,value

    num_options=0
    num_arguments=command_argument_count()
    allocate(options(num_arguments,2))
    name=''
    do i=1,num_arguments
       value=''
       call get_command_argument(i,argument)
       argument=trim(argument)
       if (argument(1:1)=='-') then
          if (name/='') then
             call put_option(options,num_options,name,value)
          else
             name=argument
          end if
       else
          if (name/='') then
             call put_option(options,num_options,name,argument)
             name=''
          end if
       end if
    end do
    if (name/='') then
       value=''
       call put_option(options,num_options,name,value)
    end if
  end subroutine get_command_line_options


  !> Helper function to create an error message for a missing option.
  !! For example <tt>Error: Missing value: option</tt>.
  character(len=256) function missing_option_message(option)
    character(len=*),intent(in) :: option
    write(missing_option_message,'(A,A)') 'Error: Missing value: ',option
  end function missing_option_message


  !> Helper function to create an error message for an invalid option.
  !! For example <tt>Error: Invalid value: option value expected</tt>.
  character(len=256) function invalid_string_option_message(option,value,expected)
    character(len=*),intent(in) :: option
    character(len=*),intent(in) :: value
    character(len=*),intent(in) :: expected
    write(invalid_string_option_message,'(A,A,A1,A,A1,A)') 'Error: Invalid value: ',option,' ',value,' ',expected
  end function invalid_string_option_message


  !> Helper function to create an error message for an invalid option.
  !! For example <tt>Error: Invalid value: option I6 expected</tt>.
  character(len=256) function invalid_integer_option_message(option,value,expected)
    character(len=*),intent(in) :: option
    integer,         intent(in) :: value
    character(len=*),intent(in) :: expected
    write(invalid_integer_option_message,'(A,A,A1,I6,A1,A)') 'Error: Invalid value: ',option,' ',value,' ',expected
  end function invalid_integer_option_message


  !> Helper function to create an error message for an invalid option.
  !! For example <tt>Error: Invalid value: option F20.10 expected</tt>.
  character(len=256) function invalid_double_option_message(option,value,expected)
    character(len=*),intent(in) :: option
    double precision,intent(in) :: value
    character(len=*),intent(in) :: expected
    write(invalid_double_option_message,'(A,A,A1,F20.10,A1,A)') 'Error: Invalid value: ',option,' ',value,' ',expected
  end function invalid_double_option_message


  !> Helper function to create an error message for an invalid option.
  !! For example <tt>Error: Invalid value: option L1 expected</tt>.
  character(len=256) function invalid_logical_option_message(option,value,expected)
    character(len=*),intent(in) :: option
    logical,         intent(in) :: value
    character(len=*),intent(in) :: expected
    write(invalid_logical_option_message,'(A,A,A1,L1,A1,A)') 'Error: Invalid value: ',option,' ',value,' ',expected
  end function invalid_logical_option_message

end module options_utils
