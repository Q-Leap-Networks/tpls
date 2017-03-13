!> twophase_io module tests.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune, David
!! Scott, Peter Spelt, Mike Jackson.
!! @version $Revision: 201 $
!! @copyright (c) 2013-2014, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, David Scott, Peter Spelt, The University of
!! Edinburgh, all rights reserved.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module tpls_io_test

  use option_names
  use tpls_io

  use fruit
  use tpls_fruit_utils

  implicit none

  double precision,dimension(:,:,:),allocatable :: out_data,in_data
  character(len=60) :: filename,description
  parameter(filename="test_tpls_io.tmp")
  parameter(description="tpls_io test file")

contains

  !> No-op.  
  subroutine setup
  end subroutine setup


  !> Deallocate arrays and delete temporary files.
  subroutine teardown
    integer :: status
    logical :: exist
    if (allocated(out_data)) then
      deallocate(out_data)
    end if
    if (allocated(in_data)) then
      deallocate(in_data)
    end if
    inquire(file=filename,exist=exist)
    if (exist) then
      open(unit=1234,iostat=status,file=filename,status='old')
      if (status==0) then
        close(1234, status='delete')
      end if
    end if
  end subroutine teardown


  !> Populate out_data with sample values.
  subroutine populate_out_data(lx,ly,mx,my,nx,ny)
    integer,intent(in) :: lx,ly,mx,my,nx,ny
    integer :: i,j,k

    allocate(out_data(lx:ly,mx:my,nx:ny))
    do k=nx,ny
      do j=mx,my
        do i=lx,ly
          out_data(i,j,k)=(i+1)*(j+2)*(k+3)
        end do
      end do
    end do
  end subroutine populate_out_data


  !> Save data and check file exists.
  subroutine test_save_data()
    integer :: lx,ly,mx,my,nx,ny
    logical :: exist
    parameter(lx=1,ly=9,mx=0,my=10,nx=1,ny=11)

    call populate_out_data(lx,ly,mx,my,nx,ny)

    call save_data(filename,description,out_data,lx,ly,mx,my,nx,ny)

    inquire(file=filename,exist=exist)
    call assert_true(exist,'File exists')
  end subroutine test_save_data


  !> Save data, load data from saved file and check loaded data equals
  !! saved data.
  subroutine test_load_data()
    integer :: lx,ly,mx,my,nx,ny
    integer :: load_lx,load_ly,load_mx,load_my,load_nx,load_ny
    integer :: ierr,i,j,k
    parameter(lx=1,ly=9,mx=0,my=10,nx=1,ny=11)

    call populate_out_data(lx,ly,mx,my,nx,ny)
    call save_data(filename,description,out_data,lx,ly,mx,my,nx,ny)

    call load_data(filename,in_data,load_lx,load_ly,load_mx,load_my,load_nx,load_ny,ierr)

    call assert_equals(0,ierr,'load_data error code')
    call assert_equals(lx,load_lx,'L dimension lower bound')
    call assert_equals(ly,load_ly,'L dimension upper bound')
    call assert_equals(mx,load_mx,'M dimension lower bound')
    call assert_equals(my,load_my,'M dimension upper bound')
    call assert_equals(nx,load_nx,'N dimension lower bound')
    call assert_equals(ny,load_ny,'N dimension upper bound')
    do k=nx,ny
      do j=mx,my
        do i=lx,ly
          call assert_equals(out_data(i,j,k),in_data(i,j,k),'Data')
        end do
      end do
    end do
  end subroutine test_load_data


  !> Save data, load data from saved file and check loaded data equals
  !! saved data.
  subroutine test_load_data_into_array()
    integer :: lx,ly,mx,my,nx,ny
    integer :: ierr,i,j,k
    parameter(lx=1,ly=9,mx=0,my=10,nx=1,ny=11)

    call populate_out_data(lx,ly,mx,my,nx,ny)
    call save_data(filename,description,out_data,lx,ly,mx,my,nx,ny)
    allocate(in_data(lx:ly,mx:my,nx:ny))

    call load_data_into_array(filename,in_data,lx,ly,mx,my,nx,ny,ierr)

    call assert_equals(0,ierr,'load_data_into_array error code')
    do k=nx,ny
      do j=mx,my
        do i=lx,ly
          call assert_equals(out_data(i,j,k),in_data(i,j,k),'Data')
        end do
      end do
    end do
  end subroutine test_load_data_into_array


  !> Save data, load data from saved file but with different expected
  !! bounds and check error codes are set.
  subroutine test_load_data_into_array_mismatch()
    integer :: lx,ly,mx,my,nx,ny
    integer :: ierr
    parameter(lx=1,ly=9,mx=0,my=10,nx=1,ny=11)

    call populate_out_data(lx,ly,mx,my,nx,ny)
    call save_data(filename,description,out_data,lx,ly,mx,my,nx,ny) 

    call load_data_into_array(filename,in_data,lx,ly+1,mx,my,nx,ny,ierr)
    call assert_equals(1,ierr,'load_data_into_array error code')

    call load_data_into_array(filename,in_data,lx,ly,mx+1,my,nx,ny,ierr)
    call assert_equals(2,ierr,'load_data_into_array error code')

    call load_data_into_array(filename,in_data,lx,ly,mx,my+1,nx,ny,ierr)
    call assert_equals(2,ierr,'load_data_into_array error code')

    call load_data_into_array(filename,in_data,lx,ly,mx,my,nx+1,ny,ierr)
    call assert_equals(3,ierr,'load_data_into_array error code')

    call load_data_into_array(filename,in_data,lx,ly,mx,my,nx,ny+1,ierr)
    call assert_equals(3,ierr,'load_data_into_array error code')
  end subroutine test_load_data_into_array_mismatch


  !> Save data and check file exists.
  subroutine test_channel_output_phi()
    integer :: l,m,n
    double precision :: dx,dy,dz
    logical :: exist
    parameter(l=9,m=5,n=5,dx=0.1,dy=0.2,dz=0.3)

    call populate_out_data(0,l,0,m,0,n)

    call channel_output_phi(out_data,l,m,n,dx,dy,dz,filename)

    inquire(file=filename,exist=exist)
    call assert_true(exist,'File exists')
  end subroutine test_channel_output_phi


  !> Save data, load data from saved file and check loaded data equals
  !! saved data.
  subroutine test_channel_input_phi()
    integer :: l,m,n
    double precision :: dx,dy,dz
    integer :: ierr,i,j,k
    parameter(l=9,m=5,n=5,dx=0.1,dy=0.2,dz=0.3)

    call populate_out_data(0,l,0,m,0,n)
    call channel_output_phi(out_data,l,m,n,dx,dy,dz,filename)
    allocate(in_data(0:l,0:m,0:n))

    call channel_input_phi(filename,in_data,l,m,n,ierr)

    call assert_equals(0,ierr,'load_data_into_array error code')
    ! Offsets in indices are due to way channel_output_phi and
    ! channel_input_phi save and load the data.
    do k=0,n-2
      do j=0,m-2
        do i=1,l-1
          call assert_equals(out_data(i+1,j,k+1),in_data(i+1,j,k+1),'Data')
        end do
      end do
    end do
  end subroutine test_channel_input_phi


  !> Save data and check file exists.
  subroutine test_channel_output_uvw()
    integer :: l,m,n
    double precision :: dx,dy,dz
    logical :: exist
    parameter(l=9,m=5,n=5,dx=0.1,dy=0.2,dz=0.3)

    call populate_out_data(0,l,0,m,0,n)

    ! Use out_data for u3,v3,w3,p
    call channel_output_uvw(out_data,out_data,out_data,out_data,l,m,n,dx,dy,dz,filename)

    inquire(file=filename,exist=exist)
    call assert_true(exist,'File exists')
  end subroutine test_channel_output_uvw


  !> Save data and check file exists.
  subroutine test_channel_output2d()
    integer :: l,m,n
    double precision :: dx,dy,dz
    logical :: exist
    parameter(l=9,m=5,n=5,dx=0.1,dy=0.2,dz=0.3)

    call populate_out_data(0,l,0,m,0,n)

    ! Use out_data for u3,v3,w3,p,phi,visc
    call channel_output2d(out_data,out_data,out_data,&
      out_data,out_data,out_data,l,m,n,dx,dy,dz,filename)

    inquire(file=filename,exist=exist)
    call assert_true(exist,'File exists')
  end subroutine test_channel_output2d


  !> Traverse a 3D domain and check row calculated is as expected for 
  !! each point.
  subroutine test_row_to_grid()
    integer :: l,m,n
    integer :: i,j,k,rl,rm,rn,row
    parameter(l=9,m=5,n=4)

    row = 1
    do k=0,n
      do j=0,m
        do i=0,l
          call row_to_grid(0,l,0,m,0,n,row,rl,rm,rn)
          call assert_equals(rl,i,'Row to L dimension')
          call assert_equals(rm,j,'Row to M dimension')
          call assert_equals(rn,k,'Row to N dimension')
          row=row+1
        end do
     end do
   end do
  end subroutine test_row_to_grid


  subroutine tpls_io_basket()
    character(len=*) :: suite_name 
    parameter(suite_name='tpls_io_test')

    call run_fruit_test_case(test_save_data,'test_save_data',&
      setup,teardown,suite_name)
    call run_fruit_test_case(test_load_data,'test_load_data',&
      setup,teardown,suite_name)
    call run_fruit_test_case(test_load_data_into_array,&
      'test_load_data_into_array',setup,teardown,suite_name)
    call run_fruit_test_case(test_load_data_into_array_mismatch,&
      'test_load_data_into_array_mismatch',setup,teardown,suite_name)
    call run_fruit_test_case(test_channel_output_phi,&
      'test_channel_output_phi',setup,teardown,suite_name)
    call run_fruit_test_case(test_channel_input_phi,&
      'test_channel_input_phi',setup,teardown,suite_name)
    call run_fruit_test_case(test_channel_output_uvw,&
      'test_channel_output_uvw',setup,teardown,suite_name)
    call run_fruit_test_case(test_channel_output2d,&
      'test_channel_output2d',setup,teardown,suite_name)
    call run_fruit_test_case(test_row_to_grid,'test_row_to_grid',&
      setup,teardown,suite_name)
  end subroutine tpls_io_basket

end module tpls_io_test
