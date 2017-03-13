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

module grid_utils_test
  use fruit
  use tpls_fruit_utils
  use grid_utils

  implicit none

  integer,dimension(:,:),allocatable :: data

contains

  !> No-op.
  subroutine setup
  end subroutine setup


  subroutine teardown
    if (allocated(data)) then
      deallocate(data)
    end if
  end subroutine teardown


  subroutine test_get_factors()
    integer :: nf
    integer,dimension(1,2) :: expected1,expected4
    integer,dimension(2,2) :: expected2
    integer,dimension(4,2) :: expected3
    DATA expected1 / 1,1 /               ! 1x1
    DATA expected2 / 8,4,1,2/            ! 8x1,4x2
    DATA expected3 / 64,32,16,8,1,2,4,8/ ! 64x1,32x2,16x4,8x8
    DATA expected4 / 17,1 /              ! 17x1

    call get_factors(1,nf,data)
    call assert_equals(1,1,'Number of factors')
    call assert_equals(expected1,data,2,1,'Factors')
    deallocate(data)
    call get_factors(8,nf,data)
    call assert_equals(2,nf,'Number of factors')
    call assert_equals(expected2,data,2,2,'Factors')
    deallocate(data)
    call get_factors(64,nf,data)
    call assert_equals(4,nf,'Number of factors')
    call assert_equals(expected3,data,4,2,'Factors')
    deallocate(data)
    call get_factors(17,nf,data)
    call assert_equals(1,nf,'Number of factors')
    call assert_equals(expected4,data,2,1,'Factors')
    deallocate(data)
  end subroutine test_get_factors


  subroutine test_get_divisors()
    integer :: nd
    integer,dimension(1,2) :: pairs1,pairs4,pairs5,expected5
    integer,dimension(2,2) :: expected1,pairs2
    integer,dimension(4,2) :: expected2,pairs3
    integer,dimension(5,2) :: expected3
    DATA pairs1    / 1,1 /                  ! 1x1
    DATA expected1 / 1,1,1,1 /              ! 1x1,1x1
    DATA pairs2    / 8,4,1,2/               ! 8x1,4x2
    DATA expected2 / 8,1,4,2,1,8,2,4/       ! 8x1,1x8,4x2,2x4
    DATA pairs3    / 64,32,16,8,1,2,4,8/    ! 64x1,32x2,16x4,8x8
    DATA expected3 / 64,32,16,8,8,1,2,4,8,8/ ! 64x1,32x2,16x4,8x8,8x8
    DATA pairs4    / 17,1 /                  ! 17x1
    DATA pairs5    / 256,144 /               ! 256x144
    DATA expected5 / 256,144 /               ! 256x144

    call get_divisors(256,152,1,pairs1,nd,data)
    call assert_equals(2,nd,'Number of divisors')
    call assert_equals(expected1,data,2,2,'Divisors')
    deallocate(data)
    call get_divisors(256,152,2,pairs2,nd,data)
    call assert_equals(4,nd,'Number of divisors')
    call assert_equals(expected2,data,4,2,'Divisors')
    deallocate(data)
    call get_divisors(256,152,4,pairs3,nd,data)
    call assert_equals(5,nd,'Number of divisors')
    call assert_equals(expected3,data,4,2,'Divisors')
    deallocate(data)
    call get_divisors(256,144,1,pairs4,nd,data)
    call assert_equals(0,nd,'Number of divisors')
    deallocate(data)
    call get_divisors(256,144,1,pairs5,nd,data)
    call assert_equals(1,nd,'Number of divisors')
    call assert_equals(expected5,data,1,2,'Divisors')
    deallocate(data)
  end subroutine test_get_divisors


  subroutine test_get_closest_pair()
    integer :: x,y
    integer,dimension(4,2) :: pairs
    DATA pairs / 32,16,8,4,1,2,4,8/   ! 32x1,16x2,8x4,4x8

    call get_closest_pair(256,152,4,pairs,x,y)
    call assert_equals(x,8,'X grid dimension')
    call assert_equals(y,4,'Y grid dimension')
    call get_closest_pair(256,151,4,pairs,x,y)
    call assert_equals(x,8,'X grid dimension')
    call assert_equals(y,4,'Y grid dimension')
    call get_closest_pair(256,2,4,pairs,x,y)
    call assert_equals(x,32,'X grid dimension')
    call assert_equals(y,1,'Y grid dimension')
    call get_closest_pair(256,33,4,pairs,x,y)
    call assert_equals(x,16,'X grid dimension')
    call assert_equals(y,2,'Y grid dimension')
    call get_closest_pair(256,65,4,pairs,x,y)
    call assert_equals(x,8,'X grid dimension')
    call assert_equals(y,4,'Y grid dimension')
    call get_closest_pair(65,256,4,pairs,x,y)
    call assert_equals(x,4,'X grid dimension')
    call assert_equals(y,8,'Y grid dimension')   
  end subroutine test_get_closest_pair


  subroutine test_get_process_grids()
    integer :: ng
    integer,dimension(4,2) :: expected1
    integer,dimension(1,2) :: expected2
    DATA expected1 / 32,16,8,4,1,2,4,8/   ! 32x1,16x2,8x4,4x8
    DATA expected2 / 32,1 /               ! 32x1

    call get_process_grids(32,256,152,ng,data)
    call assert_equals(4,ng,'Number of grids')
    call assert_equals(expected1,data,4,2,'Grids')
    deallocate(data)
    call get_process_grids(32,256,151,ng,data)
    call assert_equals(1,ng,'Number of grids')
    call assert_equals(expected2,data,2,1,'Grids')
    deallocate(data)
    call get_process_grids(32,257,152,ng,data)
    call assert_equals(0,ng,'Number of grids')
    deallocate(data)
  end subroutine test_get_process_grids


  subroutine test_get_process_grid()
    integer :: x,y

    call get_process_grid(32,256,152,x,y)
    call assert_equals(x,8,'X grid dimension')
    call assert_equals(y,4,'Y grid dimension')
    call get_process_grid(32,256,151,x,y)
    call assert_equals(x,32,'X grid dimension')
    call assert_equals(y,1,'Y grid dimension')
    call get_process_grid(32,257,152,x,y)
    call assert_equals(x,0,'X grid dimension')
    call assert_equals(y,0,'Y grid dimension')
  end subroutine test_get_process_grid


  subroutine test_get_multiples()
    integer :: i
    integer,dimension(1,2)  :: multiples1, expected1
    integer,dimension(30,2) :: multiples30, expected30

    expected1(1,1)=3
    expected1(1,2)=4
    call get_multiples(3,4,1,multiples1)
    call assert_equals(expected1,multiples1,1,2,'Multiples')
    do i=1,30
      expected30(i,1)=i*56
      expected30(i,2)=i*78
    end do
    call get_multiples(56,78,30,multiples30)
    call assert_equals(expected30,multiples30,30,2,'Multiples')
  end subroutine test_get_multiples


  subroutine grid_utils_basket()
    character(len=*) :: suite_name 
    parameter(suite_name='grid_utils_test')

    call run_fruit_test_case(test_get_factors,'test_get_factors',&
      setup,teardown,suite_name)
    call run_fruit_test_case(test_get_divisors,'test_get_divisors',&
      setup,teardown,suite_name)
    call run_fruit_test_case(test_get_closest_pair,'test_get_closest_pair',&
       setup,teardown,suite_name)
    call run_fruit_test_case(test_get_process_grids,'test_get_process_grids',&
      setup,teardown,suite_name)
    call run_fruit_test_case(test_get_process_grid,'test_get_process_grid',&
      setup,teardown,suite_name)
    call run_fruit_test_case(test_get_multiples,'test_get_multiples',&
      setup,teardown,suite_name)
  end subroutine grid_utils_basket

end module grid_utils_test
