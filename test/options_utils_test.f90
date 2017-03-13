!> grid_utils module tests.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune, David
!! Scott, Peter Spelt, Mike Jackson.
!! @version $Revision: 252 $
!! @copyright (c) 2013-2014, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, David Scott, Peter Spelt, The University of
!! Edinburgh, all rights reserved.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module options_utils_test
  use fruit
  use tpls_fruit_utils
  use options_utils

  implicit none

  character(len=256),dimension(:,:),allocatable :: options
  integer                                       :: num_options
  character(len=60)                             :: filename
  parameter(filename="test_options_utils.tmp")

contains

  !> No-op.
  subroutine setup
    num_options=0
  end subroutine setup


  subroutine teardown
    integer :: status
    logical :: exist

    if (allocated(options)) then
      deallocate(options)
    end if
    inquire(file=filename,exist=exist)
    if (exist) then
      open(unit=1234,iostat=status,file=filename,status='old')
      if (status==0) then
        close(1234, status='delete')
      end if
    end if
  end subroutine teardown


  subroutine test_has_option()
    num_options=3
    allocate(options(num_options,2))
    options(1,1)='option1'
    options(1,2)='value1'
    options(2,1)='option2'
    options(2,2)='value2'
    options(3,1)='option3'
    options(3,2)='value3'
    call assert_true(has_option(options,num_options,'option1'),'Has option')
    call assert_true(has_option(options,num_options,'option2'),'Has option')
    call assert_true(has_option(options,num_options,'option3'),'Has option')
    call assert_false(has_option(options,num_options,'option4'),'Does not have option')
  end subroutine test_has_option


  subroutine test_get_option()
    integer             :: ierr
    integer             :: int_value
    double precision    :: double_value, expected_double_value, fudge
    character(len=256 ) :: str_value
    logical             :: logical_value

    num_options=4
    allocate(options(num_options,2))
    options(1,1)='string'
    options(1,2)='blah'
    options(2,1)='int'
    options(2,2)='2'
    options(3,1)='double'
    options(3,2)='3.45'
    options(4,1)='logical'
    options(4,2)='true'

    ! String options
    call get_option(options,num_options,'string',str_value,ierr)
    call assert_equals(0,ierr,'get_option error code')
    call assert_equals('blah',str_value,'String')

    call get_string_option(options,num_options,'string',str_value,ierr)
    call assert_equals(0,ierr,'get_string_option error code')
    call assert_equals('blah',str_value,'String value')

    call get_typed_option(options,num_options,'string',str_value,ierr)
    call assert_equals(0,ierr,'get_string_option error code')
    call assert_equals('blah',str_value,'String value')

    ! Integer options
    call get_option(options,num_options,'int',str_value,ierr)
    call assert_equals('2',str_value,'Integer')

    call get_integer_option(options,num_options,'int',int_value,ierr)
    call assert_equals(0,ierr,'get_integer_option error code')
    call assert_equals(2,int_value,'Integer value')

    call get_typed_option(options,num_options,'int',int_value,ierr)
    call assert_equals(0,ierr,'get_integer_option error code')
    call assert_equals(2,int_value,'Integer value')

    ! Double options
    call get_option(options,num_options,'double',str_value,ierr)
    call assert_equals('3.45',str_value,'Double')

    call get_double_option(options,num_options,'double',double_value,ierr)
    call assert_equals(0,ierr,'get_double_option error code')
    ! Passing double as a constant gives a compilation error.
    ! Error: There is no specific subroutine for the generic 'assert_equals'.
    expected_double_value=3.45
    fudge=0.01
    call assert_equals(expected_double_value,double_value,fudge,'Double value')

    call get_typed_option(options,num_options,'double',double_value,ierr)
    call assert_equals(0,ierr,'get_double_option error code')
    call assert_equals(expected_double_value,double_value,fudge,'Double value')

    ! Logical options
    call get_option(options,num_options,'logical',str_value,ierr)
    call assert_equals('true',str_value,'Logical')

    call get_logical_option(options,num_options,'logical',logical_value,ierr)
    call assert_equals(0,ierr,'get_logical_option error code')
    call assert_true(logical_value,'Logical value')

    call get_typed_option(options,num_options,'logical',logical_value,ierr)
    call assert_equals(0,ierr,'get_logical_option error code')
    call assert_true(logical_value,'Logical value')
  end subroutine test_get_option


  subroutine test_get_no_such_option()
    integer             :: ierr
    integer             :: int_value
    double precision    :: double_value
    character(len=256 ) :: str_value
    logical             :: logical_value

    num_options=1
    allocate(options(num_options,2))
    options(1,1)='option'
    options(1,2)='value'

    call get_option(options,num_options,'blah',str_value,ierr)
    call assert_equals(option_not_found_err,ierr,'get_option not found')

    ! String options
    call get_string_option(options,num_options,'blah',str_value,ierr)
    call assert_equals(option_not_found_err,ierr,'get_string_option not found')

    call get_typed_option(options,num_options,'blah',str_value,ierr)
    call assert_equals(option_not_found_err,ierr,'get_string_option not found')

    ! Integer options
    call get_integer_option(options,num_options,'blah',int_value,ierr)
    call assert_equals(option_not_found_err,ierr,'get_integer_option not found')

    call get_typed_option(options,num_options,'blah',int_value,ierr)
    call assert_equals(option_not_found_err,ierr,'get_integer_option not found')

    ! Double options
    call get_double_option(options,num_options,'blah',double_value,ierr)
    call assert_equals(option_not_found_err,ierr,'get_double_option not found')

    call get_typed_option(options,num_options,'blah',double_value,ierr)
    call assert_equals(option_not_found_err,ierr,'get_double_option not found')

    ! Logical options
    call get_logical_option(options,num_options,'blah',logical_value,ierr)
    call assert_equals(option_not_found_err,ierr,'get_logical_option not found')

    call get_typed_option(options,num_options,'blah',logical_value,ierr)
    call assert_equals(option_not_found_err,ierr,'get_logical_option not found')
  end subroutine test_get_no_such_option


  !> Add 5 options using type-specific put subroutines.
  subroutine add_test_options()
    double precision  :: double_value

    call put_string_option(options,num_options,'string','blah')
    call put_integer_option(options,num_options,'int',2)
    ! Passing double as a constant gives a compilation error.
    ! Error: Type mismatch in argument 'value' at (1); passed REAL(4) to REAL(8). 
    double_value=3.45
    call put_double_option(options,num_options,'double',double_value)
    call put_logical_option(options,num_options,'logical',.true.)
    call put_option(options,num_options,'option','value')
  end subroutine add_test_options


  !> Add 4 options using the put_typed_option interface.
  subroutine add_typed_test_options()
    double precision  :: double_value

    call put_typed_option(options,num_options,'string_type','blah_type')
    call put_typed_option(options,num_options,'int_type',6)
    double_value=7.89
    call put_typed_option(options,num_options,'double_type',double_value)
    call put_typed_option(options,num_options,'logical_type',.false.)
  end subroutine add_typed_test_options


  !> Check options added using add_test_options are present.
  subroutine check_test_options()
    integer             :: ierr
    integer             :: int_value
    double precision    :: double_value, expected_double_value, fudge
    character(len=256 ) :: str_value
    logical             :: logical_value

    ! String options
    call get_typed_option(options,num_options,'string',str_value,ierr)
    call assert_equals(0,ierr,'get_typed_option error code')
    call assert_equals('blah',str_value,'String value')

    ! Integer options
    call get_typed_option(options,num_options,'int',int_value,ierr)
    call assert_equals(0,ierr,'get_typed_option error code')
    call assert_equals(2,int_value,'Integer value')

    ! Double options
    expected_double_value=3.45
    call get_typed_option(options,num_options,'double',double_value,ierr)
    call assert_equals(0,ierr,'get_typed_option error code')
    fudge=0.01
    call assert_equals(expected_double_value,double_value,fudge,'Double value')

    ! Logical options
    call get_typed_option(options,num_options,'logical',logical_value,ierr)
    call assert_equals(0,ierr,'get_typed_option error code')
    call assert_true(logical_value,'Logical value')
  end subroutine check_test_options

  !> Check options added using add_test_options are present.
  subroutine check_typed_test_options()
    integer             :: ierr
    integer             :: int_value
    double precision    :: double_value, expected_double_value, fudge
    character(len=256 ) :: str_value
    logical             :: logical_value

    ! String options
    call get_typed_option(options,num_options,'string_type',str_value,ierr)
    call assert_equals(0,ierr,'get_typed_option error code')
    call assert_equals('blah_type',str_value,'String value')

    ! Integer options
    call get_typed_option(options,num_options,'int_type',int_value,ierr)
    call assert_equals(0,ierr,'get_typed_option error code')
    call assert_equals(6,int_value,'Integer value')

    ! Double options
    fudge=0.01
    expected_double_value=7.89
    call get_typed_option(options,num_options,'double_type',double_value,ierr)
    call assert_equals(0,ierr,'get_typed_option error code')
    call assert_equals(expected_double_value,double_value,fudge,'Double value')

    ! Logical options
    call get_typed_option(options,num_options,'logical_type',logical_value,ierr)
    call assert_equals(0,ierr,'get_typed_option error code')
    call assert_false(logical_value,'Logical value')
  end subroutine check_typed_test_options


  subroutine test_put_option()
    num_options=0
    allocate(options(9,2))
    call add_test_options()
    call add_typed_test_options()

    call assert_equals(9, num_options,'Number of options')
    call check_test_options()
    call check_typed_test_options()
  end subroutine test_put_option


  !> Test that putting multiple options with the same name results in
  !! existing options with that name being overwritten.
  subroutine test_multiple_put_option()
    integer             :: ierr
    character(len=256 ) :: value

    allocate(options(3,2))
    call put_option(options,num_options,'key','blah')
    call get_option(options,num_options,'key',value,ierr)
    call assert_equals(1,num_options,'Number of options')
    call assert_equals(0,ierr,'get_option error code')
    call assert_equals('blah',value,'Value')

    call put_option(options,num_options,'key_two','two_blah')
    call get_option(options,num_options,'key_two',value,ierr)
    call assert_equals(2,num_options,'Number of options')
    call assert_equals(0,ierr,'get_option error code')
    call assert_equals('two_blah',value,'Value')

    ! Overwrite value
    call put_option(options,num_options,'key','tpls')
    call get_option(options,num_options,'key',value,ierr)
    call assert_equals(2,num_options,'Number of options')
    call assert_equals(0,ierr,'get_option error code')
    call assert_equals('tpls',value,'Value')

    ! Overwrite value
    call put_option(options,num_options,'key','dim')
    call get_option(options,num_options,'key',value,ierr)
    call assert_equals(2,num_options,'Number of options')
    call assert_equals(0,ierr,'get_option error code')
    call assert_equals('dim',value,'Value')

    ! Overwrite value
    call put_option(options,num_options,'key_two','69')
    call get_option(options,num_options,'key_two',value,ierr)
    call assert_equals(2,num_options,'Number of options')
    call assert_equals(0,ierr,'get_option error code')
    call assert_equals('69',value,'Value')

    call put_option(options,num_options,'key_3','This is another value')
    call get_option(options,num_options,'key_3',value,ierr)
    call assert_equals(3,num_options,'Number of options')
    call assert_equals(0,ierr,'get_option error code')
    call assert_equals('This is another value',value,'Value')
  end subroutine test_multiple_put_option


  subroutine test_save_and_load_options()
    character(len=256),dimension(3) :: header
    integer                         :: ierr

    header(1) = "options_utils_test header"
    header(2) = "Automatically created by options_utils_test"
    header(3) = "Options now follow..."

    num_options=0
    allocate(options(9,2))
    call add_test_options()
    call add_typed_test_options()
    call save_options(filename,options,num_options,header,3)
    deallocate(options)
    num_options=0

    call load_options(filename,options,num_options,ierr)
    call assert_equals(0,ierr,'Load options')
    call assert_equals(9, num_options,'Number of options')
    call check_test_options()
    call check_typed_test_options()
  end subroutine test_save_and_load_options


  subroutine test_save_append_load_options()
    character(len=256),dimension(3) :: header
    integer                         :: ierr

    header(1) = "options_utils_test header"
    header(2) = "Automatically created by options_utils_test"
    header(3) = "Options now follow..."

    num_options=0
    allocate(options(5,2))
    call add_test_options()
    call save_options(filename,options,num_options,header,3)
    deallocate(options)
  
    num_options=0
    allocate(options(4,2))
    call add_typed_test_options()
    call append_options(filename,options,num_options)
    deallocate(options)

    num_options=0
    call load_options(filename,options,num_options,ierr)
    call assert_equals(0,ierr,'Load options')
    call assert_equals(9, num_options,'Number of options')
    call check_test_options()
    call check_typed_test_options()
  end subroutine test_save_append_load_options


  subroutine test_load_options_not_found()
    integer :: ierr

    call load_options('blah.tmp',options,num_options,ierr)
    call assert_equals(option_file_not_found_err,ierr,'Option file not found')
  end subroutine test_load_options_not_found


  subroutine test_load_options_format_error()
    integer :: file_id
    parameter(file_id=135)
    integer :: ierr

    open(unit=file_id,file=filename,action='WRITE',status='REPLACE')
    write(file_id,'(A)') '# Comment'
    write(file_id,'(A)') '# Comment'
    write(file_id,'(A)') '# Comment'
    write(file_id,'(A)') 'key value'
    write(file_id,'(A)') 'integer 123'
    write(file_id,'(A)') 'k v'
    write(file_id,'(A)') 'invalidvalue'
    close(file_id)

    call load_options(filename,options,num_options,ierr)
    call assert_equals(option_format_err,ierr,'Option file format error')
  end subroutine test_load_options_format_error


  subroutine test_load_duplicates()
    integer            :: file_id
    parameter(file_id=135)
    integer            :: ierr
    character(len=256) :: value

    open(unit=file_id,file=filename,action='WRITE',status='REPLACE')
    write(file_id,'(A)') 'key 123'
    write(file_id,'(A)') 'key 456'
    write(file_id,'(A)') 'key 789'
    write(file_id,'(A)') 'key abc'
    write(file_id,'(A)') 'key check'
    close(file_id)

    call load_options(filename,options,num_options,ierr)
    call assert_equals(0,ierr,'Load options')
    call assert_equals(1,num_options,'Number of options')
    call get_option(options,num_options,'key',value,ierr)
    call assert_equals(0,ierr,'get_option error code')
    call assert_equals('check',value,'String value')
  end subroutine test_load_duplicates


  subroutine options_utils_basket()
    character(len=*) :: suite_name 
    parameter(suite_name='options_utils_test')

    call run_fruit_test_case(test_has_option,'test_has_option',&
      setup,teardown,suite_name)
    call run_fruit_test_case(test_get_option,'test_get_option',&
      setup,teardown,suite_name)
    call run_fruit_test_case(test_get_no_such_option,'test_get_no_such_option',&
      setup,teardown,suite_name)
    call run_fruit_test_case(test_put_option,'test_put_option',&
      setup,teardown,suite_name)
    call run_fruit_test_case(test_multiple_put_option,&
      'test_multiple_put_option',setup,teardown,suite_name)
    call run_fruit_test_case(test_save_and_load_options,&
      'test_save_and_load_options',setup,teardown,suite_name)
    call run_fruit_test_case(test_save_append_load_options,&
      'test_save_append_load_options',setup,teardown,suite_name)
    call run_fruit_test_case(test_load_options_not_found,&
      'test_load_options_not_found',setup,teardown,suite_name)
    call run_fruit_test_case(test_load_options_format_error,&
      'test_load_options_format_error',setup,teardown,suite_name)
    call run_fruit_test_case(test_load_duplicates,&
      'test_load_duplicates',setup,teardown,suite_name)
  end subroutine options_utils_basket

end module options_utils_test
