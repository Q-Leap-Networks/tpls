!> Two-phase level set input and output subroutines.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune, 
!! Toni Collis, David Scott, Peter Spelt, Mike Jackson.
!! @version $Revision: 328 $
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, David Scott, Toni Collis, Peter Spelt,
!! The University of Edinburgh, all rights reserved.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module tpls_io

  use option_names
  use options_utils
  use netcdf
  use mpi
  use petsc

  implicit none

contains

  !> Given a directory and a local file name, return the file path.
  subroutine get_filepath(directory, filename, filepath)
    character(len=*),intent(in)  :: directory !< Directory.
    character(len=*),intent(in)  :: filename  !< Local file name.
    character(len=*),intent(out) :: filepath  !< directory/filename.

    filepath=trim(directory)//'/'//trim(filename)
  end subroutine get_filepath


  !> Save 3D data to a file.
  !! The first two lines of the file are meta-data of form:
  !! \code
  !! DATA: description
  !! DIMENSIONS: Lmin Lmax Mmin Mmax Nmin Nmax
  !! lx ly mx my nx ny
  !! \endcode
  !! This is followed by (ly-lx+1)*(my-mx+1)*(ny-nx+1) rows of data
  !! representing each value of data from data(lx,mx,nx) to data(ly,
  !! my,ny) where the rth row contains (rl,rm,rn) where: 
  !! * rl = ((row - 1) % (ly - lx + 1)) + lx
  !! * rm = (((row - 1) / (ly - lx + 1)) % (my - mx + 1)) + mx
  !! * rn = (((row - 1) / ((ly - lx + 1) * (my - mx + 1)))) + nx
  subroutine save_data(filename,description,data,lx,ly,mx,my,nx,ny)
    implicit none
    character(len=*), intent(in) :: filename !< File name.
    character(len=*), intent(in) :: description !< Contents description.
    integer,          intent(in) :: lx !< Minimum L dimension bound.
    integer,          intent(in) :: ly !< Maximum L dimension bound.
    integer,          intent(in) :: mx !< Minimum M dimension bound.
    integer,          intent(in) :: my !< Maximum M dimension bound.
    integer,          intent(in) :: nx !< Minimum N dimension bound.
    integer,          intent(in) :: ny !< Maximum N dimension bound.
    double precision, intent(in) :: data(lx:ly,mx:my,nx:ny) !< Data to save.
    integer :: i,j,k

    open(unit=9,file=filename,status='unknown')
    write(9,*) 'DATA: ',description
    write(9,*) 'DIMENSIONS: Lmin Lmax Mmin Mmax Nmin Nmax'
    write(9,*) lx,ly,mx,my,nx,ny
    do k=nx,ny
       do j=mx,my
          do i=lx,ly
             write(9,*) data(i,j,k)
          end do
       end do
    end do
    close(9)
    return
  end subroutine save_data


  !> Load 3D data from a file.
  !! This is the complement of save_data and expects the input data to
  !! be formatted in the way output by that sub-routine.
  !! The array dimensions are read from the file and the array
  !! is allocated.
  !! An error code is set if the number of rows in the file does not
  !! match the dimensions recorded within the file.
  subroutine load_data(filename,data,lx,ly,mx,my,nx,ny,ierr)
    implicit none
    character(len=*), intent(in) :: filename !< File name.
    double precision, allocatable, dimension(:,:,:), intent(out) :: data(:,:,:) !< Data loaded.
    integer, intent(out) :: lx !< Minimum L dimension bound.
    integer, intent(out) :: ly !< Maximum L dimension bound.
    integer, intent(out) :: mx !< Minimum M dimension bound.
    integer, intent(out) :: my !< Maximum M dimension bound.
    integer, intent(out) :: nx !< Minimum N dimension bound.
    integer, intent(out) :: ny !< Maximum N dimension bound.
    integer, intent(out) :: ierr !< 0 if no errors,1 otherwise.
    integer :: i,j,k,rows

    open(unit=9,file=filename,action='READ',status='OLD')
    read(9,*) ! Skip meta-data line.
    read(9,*) ! Skip meta-data line.
    read(9,*) lx,ly,mx,my,nx,ny
    allocate(data(lx:ly,mx:my,nx:ny))
    rows = 0
    do k=nx,ny
       do j=mx,my
          do i=lx,ly
             read(9,*) data(i,j,k)
             rows = rows + 1
          end do
       end do
    end do
    close(9)
    ierr = 0
    if (rows /= (ly-lx+1)*(my-mx+1)*(ny-nx+1)) then
       ierr = 1
    end if
    return
  end subroutine load_data


  !> Load 3D data from a file into an array.
  !! This is the complement of save_data and expects the input data to
  !! be formatted in the way output by that sub-routine.
  !! The array is assumed to be allocated.
  !! An error code is set:
  !! - 1 if the l bounds of the file does not match the given l bounds.
  !! - 2 if the m bounds of the file does not match the given m bounds.
  !! - 3 if the n bounds of the file does not match the given n bounds.
  subroutine load_data_into_array(filename,data,lx,ly,mx,my,nx,ny,ierr)
    implicit none
    character(len=*), intent(in) :: filename !< File name.
    double precision, allocatable, dimension(:,:,:), intent(inout) :: data(:,:,:) !< Loaded data.
    integer, intent(in)  :: lx !< Minimum L dimension bound.
    integer, intent(in)  :: ly !< Maximum L dimension bound.
    integer, intent(in)  :: mx !< Minimum M dimension bound.
    integer, intent(in)  :: my !< Maximum M dimension bound.
    integer, intent(in)  :: nx !< Minimum N dimension bound.
    integer, intent(in)  :: ny !< Maximum N dimension bound.
    integer, intent(out) :: ierr !< 0 if no errors,1 otherwise.
    integer :: flx,fly,fmx,fmy,fnx,fny,i,j,k

    open(unit=9,file=filename,action='READ',status='OLD')
    read(9,*) ! Skip meta-data line.
    read(9,*) ! Skip meta-data line.
    read(9,*) flx,fly,fmx,fmy,fnx,fny
    if ((flx /= lx).or.((fly /= ly))) then
       ierr = 1
       close(9)
       return
    end if
    if ((fmx /= mx).or.((fmy /= my))) then
       ierr = 2
       close(9)
       return
    end if
    if ((fnx /= nx).or.((fny /= ny))) then
       ierr = 3
       close(9)
       return
    end if
    do k=nx,ny
       do j=mx,my
          do i=lx,ly
             read(9,*) data(i,j,k)
          end do
       end do
    end do
    close(9)
    ierr = 0
    return
  end subroutine load_data_into_array


  !> Save info needed for output file processing
  !! The first line of the file are meta-data of form:
  !! \code
  !! DIMENSIONS: Lmin Lmax Mmin Mmax Nmin Nmax
  !! lx ly mx my nx ny
  !! DIMENSIONS: Dx,Dy,Dz
  !! dx,dy,dz
  !! \endcode
  subroutine save_backup_info(filename,lx,ly,mx,my,nx,ny,dx,dy,dz)
    implicit none
    character(len=*), intent(in) :: filename !< File name.
    integer,          intent(in) :: lx !< Minimum L dimension bound.
    integer,          intent(in) :: ly !< Maximum L dimension bound.
    integer,          intent(in) :: mx !< Minimum M dimension bound.
    integer,          intent(in) :: my !< Maximum M dimension bound.
    integer,          intent(in) :: nx !< Minimum N dimension bound.
    integer,          intent(in) :: ny !< Maximum N dimension bound.
    integer,          intent(in) :: dx !< X dimension element size.
    integer,          intent(in) :: dy !< Y dimension element size.
    integer,          intent(in) :: dz !< Z dimension element size.

    open(unit=9,file=filename,status='unknown')
    write(9,*) 'DIMENSIONS: Lmin Lmax Mmin Mmax Nmin Nmax'
    write(9,*) lx,ly,mx,my,nx,ny
    write(9,*) 'DIMENSIONS: Dx,Dy,Dz'
    write(9,*) dx,dy,dz
    close(9)
    return
  end subroutine save_backup_info


  !> Save level-set function data to a file.
  !! The first two lines of the file are meta-data of form:
  !! \code
  !! VARIABLES="X ","Y ","Z ","Phi"
  !! ZONE T="Floor",I=l-1 J=m-1 K=n-1 F=POINT
  !! \endcode
  !! This is followed by (l-1)*(m-1)*(n-1) rows of form: 
  !! \code
  !! X Y Z phi(L+1,M,N+1)
  !! \endcode
  !! where 
  !! * L:1..l-1
  !! * M:0..m-2
  !! * N:0..n-2 
  subroutine channel_output_phi(phi,l,m,n,dx,dy,dz,filename)
    implicit none
    character(len=*), intent(in) :: filename !< File name.
    integer,          intent(in) :: l !< L dimension bound.
    integer,          intent(in) :: m !< M dimension bound.
    integer,          intent(in) :: n !< N dimension bound.
    double precision, intent(in) :: phi(0:l,0:m,0:n) !< Level-set function.
    double precision, intent(in) :: dx !< X dimension element size.
    double precision, intent(in) :: dy !< Y dimension element size.
    double precision, intent(in) :: dz !< Z dimension element size.
    integer          :: i, j, k ! Loop variables.
    double precision :: x, y, z ! Temporary variables for coordinate vector.
    !write(*,*) "Opening phi output file",filename
    !write(*,*) "************************************"
    open(unit=11,file=filename,status='unknown')
    write(11,*)'VARIABLES="X ","Y ","Z ","Phi"'
    write(11,*)'ZONE T="Floor", I=',l-1,' J=',m-1,' K=',n-1,' F=POINT'
    do k=0,n-2
       z=dz/2.d0+dz*k              
       do j=0,m-2
          y=dy/2.d0+dy*j
          do i=1,l-1
             x=dx/2.d0+dx*i
             write(11,*) x,y,z,phi(i+1,j,k+1)
          end do
       end do
    end do
    close(11)
    return
  end subroutine channel_output_phi


  !> Save level-set function data to a netcdf-hdf5 file.
  subroutine output_phi_hdf5(dataname,sx,ex,sy,ey,sz_p,ez_p,&
       phi_local,nprocs_x,nprocs_y,nprocs_z,rank,iteration)
    implicit none
    character(len=*), intent(in) :: dataname !< File name.
    integer,          intent(in) :: sx !< x dimension lower bound for local array
    integer,          intent(in) :: ex !< x dimension upper bound for local array
    integer,          intent(in) :: sy !< y dimension lower bound for local array
    integer,          intent(in) :: ey !< y dimension upper bound for local array
    integer,          intent(in) :: sz_p !< z dimension lower bound for local array
    integer,          intent(in) :: ez_p !< z dimension upper bound for local array
    double precision, intent(in) :: phi_local(1:(ex-sx+1),1:(ey-sy+1),1:(ez_p-sz_p+1)) !< Level-set function local array
    integer,          intent(in) :: rank !< process ID
    integer,          intent(in) :: nprocs_x !< number of procs in x dim
    integer,          intent(in) :: nprocs_y !< number of procs in y dim
    integer,          intent(in) :: nprocs_z !< number of procs in z dim
    integer,          intent(in) :: iteration !< current iteration
    !Netcdf variables
    integer cmode, ncid, varid, dimid(3)
    integer psizes(3), gsizes(3), start(3), count(3)

    !dataname = "phi"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_p,ez_p,&
       phi_local,nprocs_x,nprocs_y,nprocs_z,rank,iteration)

    return
  end subroutine output_phi_hdf5


  !> Output 3D array with hdf5
  subroutine output_3D_hdf5(dataname,sx,ex,sy,ey,sz,ez,&
       data,nprocs_x,nprocs_y,nprocs_z,rank,iteration)
    implicit none
    character(len=*), intent(in) :: dataname !< File name.
    integer,          intent(in) :: sx !< x dimension lower bound for local array
    integer,          intent(in) :: ex !< x dimension upper bound for local array
    integer,          intent(in) :: sy !< y dimension lower bound for local array
    integer,          intent(in) :: ey !< y dimension upper bound for local array
    integer,          intent(in) :: sz !< z dimension lower bound for local array
    integer,          intent(in) :: ez !< z dimension upper bound for local array
    double precision, intent(in) :: data(1:(ex-sx+1),1:(ey-sy+1),1:(ez-sz+1)) !< Local array to be printed out to global data file
    integer,          intent(in) :: rank !< process ID
    integer,          intent(in) :: nprocs_x !< number of procs in x dim
    integer,          intent(in) :: nprocs_y !< number of procs in y dim
    integer,          intent(in) :: nprocs_z !< number of procs in z dim
    integer,          intent(in) :: iteration !< iteration number
    !Netcdf variables
    integer cmode, ncid, varid, dimid(3)
    integer psizes(3), gsizes(3), start(3), count(3)
    character(len=60) :: filename

    write(filename, '(A,A1,I4.4,A3)') trim(dataname), '_', iteration, '.nc'
    ! indicate to use PnetCDF to carry out I/O
    cmode = IOR(NF90_NETCDF4, NF90_MPIIO)
    ! Create output file
    call check(nf90_create(filename, cmode, ncid, comm=PETSC_COMM_WORLD,&
         info=MPI_INFO_NULL), 'creating file: ')
    gsizes(1) = (ex-sx+1)*nprocs_x
    gsizes(2) = (ey-sy+1)*nprocs_y
    gsizes(3) = (ez-sz+1)
    psizes(1) = nprocs_x
    psizes(2) = nprocs_y
    psizes(3) = nprocs_z

    ! Printing out a 3D array of phi array (phi_local(x,y,z))
    ! define dimensions x and y
    call check(nf90_def_dim(ncid, "x", gsizes(1), dimid(1)), 'In nf_def_dim X: ')
    call check(nf90_def_dim(ncid, "y", gsizes(2), dimid(2)), 'In nf_def_dim Y: ')
    call check(nf90_def_dim(ncid, "z", gsizes(3), dimid(3)), 'In nf_def_dim Z: ')
    ! define a 3D variable of type double
    call check(nf90_def_var(ncid, dataname, NF90_DOUBLE, dimid, varid), 'In nf_def_var: ')
    ! exit define mode
    call check(nf90_enddef(ncid), 'In nf_enddef: ')
    !Now in NETCDF Data Mode
    ! set to use MPI/PnetCDF collective I/O
    call check(nf90_var_par_access(ncid, varid, NF90_COLLECTIVE),&
         'In nf_var_par_access: ')
    !Set up start and end points for each process' data
    start(1) = sx
    start(2) = sy
    start(3) = sz+1
    count(1) = (ex-sx)+1
    count(2) = (ey-sy)+1
    count(3) = ez-sz+1
    call check(nf90_put_var(ncid, varid, data, start=start, count=count),&
         'In nf_put_vara_int: ')
    ! close the file
    call check(nf90_close(ncid), 'In nf_close: ')
    return
  end subroutine output_3D_hdf5


  !> Save velocity and pressure data to a file.
  ! Calls output_3D_hdf5 for each array
  ! Creates a new netcdf file for each data array
  subroutine output_uvw_hdf5(sx,ex,sy,ey,sz_uv,ez_uv,sz_w,ez_w,sz_p,ez_p,&
       output_u,v2,output_w,output_p,nprocs_x,nprocs_y,nprocs_z,rank, iteration)
    implicit none
    integer,          intent(in) :: sx !< x dimension lower bound for local array
    integer,          intent(in) :: ex !< x dimension upper bound for local array
    integer,          intent(in) :: sy !< y dimension lower bound for local array
    integer,          intent(in) :: ey !< y dimension upper bound for local array
    integer,          intent(in) :: sz_uv !< z dimension lower bound for local array
    integer,          intent(in) :: ez_uv !< z dimension upper bound for local array
    integer,          intent(in) :: sz_w !< z dimension lower bound for local array
    integer,          intent(in) :: ez_w !< z dimension upper bound for local array
    integer,          intent(in) :: sz_p !< z dimension lower bound for local array
    integer,          intent(in) :: ez_p !< z dimension upper bound for local array
    double precision, intent(in) :: output_u(1:(ex-sx+1),1:(ey-sy+1),1:(ez_uv-sz_uv+1)) !< u data
    double precision, intent(in) :: v2(1:(ex-sx+1),1:(ey-sy+1),1:(ez_uv-sz_uv+1)) !< v data
    double precision, intent(in) :: output_w(1:(ex-sx+1),1:(ey-sy+1),1:(ez_w-sz_w+1)) !< w data
    double precision, intent(in) :: output_p(1:(ex-sx+1),1:(ey-sy+1),1:(ez_p-sz_p+1)) !< pressure
    integer,          intent(in) :: rank !< process ID
    integer,          intent(in) :: nprocs_x !< number of procs in x dim
    integer,          intent(in) :: nprocs_y !< number of procs in y dim
    integer,          intent(in) :: nprocs_z !< number of procs in z dim
    integer,          intent(in) :: iteration !< iteration number
    character(len=60) :: dataname
   !Netcdf variables
    integer cmode, ncid, varid, dimid(3)
    integer psizes(3), gsizes(3), start(3), count(3)

    dataname = "u_global"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_uv,ez_uv,&
       output_u,nprocs_x,nprocs_y,nprocs_z,rank,iteration)

    dataname = "v_global"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_uv,ez_uv,&
       v2,nprocs_x,nprocs_y,nprocs_z,rank,iteration)

    dataname = "w_global"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_w,ez_w,&
       output_w,nprocs_x,nprocs_y,nprocs_z,rank,iteration)

    dataname = "pres_global"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_p,ez_p,&
       output_p,nprocs_x,nprocs_y,nprocs_z,rank,iteration)

    return
  end subroutine output_uvw_hdf5


  !> Save a back-up file.
  subroutine output_backup_hdf5(backup_no,dx,dy,dz,sx,ex,sy,ey,sz_p,ez_p,sz_uv,&
             ez_uv,sz_w,ez_w,u2,v2,w2,conv0_u,conv1_u,conv0_v,&
             conv1_v,conv0_w,conv1_w,pres_old,pres,phi2,conv0_phi,&
             conv1_phi,nprocs_x,nprocs_y,nprocs_z,rank)
    implicit none
    integer,          intent(in) :: backup_no
    double precision, intent(in) :: dx
    double precision, intent(in) :: dy
    double precision, intent(in) :: dz
    integer,          intent(in) :: sx !< x dimension lower bound for local array
    integer,          intent(in) :: ex !< x dimension upper bound for local array
    integer,          intent(in) :: sy !< y dimension lower bound for local array
    integer,          intent(in) :: ey !< y dimension upper bound for local array
    integer,          intent(in) :: sz_p !< z dimension lower bound for local array
    integer,          intent(in) :: ez_p !< z dimension upper bound for local array
    integer,          intent(in) :: sz_uv !< z dimension lower bound for local array
    integer,          intent(in) :: ez_uv !< z dimension upper bound for local array
    integer,          intent(in) :: sz_w !< z dimension lower bound for local array
    integer,          intent(in) :: ez_w !< z dimension lower bound for local array
    double precision, intent(in) :: u2(1:(ex-sx+1),1:(ey-sy+1),1:(ez_uv-sz_uv+1)) !<
    double precision, intent(in) :: conv0_u(1:(ex-sx+1),1:(ey-sy+1),1:(ez_uv-sz_uv+1)) !<
    double precision, intent(in) :: conv1_u(1:(ex-sx+1),1:(ey-sy+1),1:(ez_uv-sz_uv+1)) !<
    double precision, intent(in) :: v2(1:(ex-sx+1),1:(ey-sy+1),1:(ez_uv-sz_uv+1)) !<
    double precision, intent(in) :: conv0_v(1:(ex-sx+1),1:(ey-sy+1),1:(ez_uv-sz_uv+1)) !<
    double precision, intent(in) :: conv1_v(1:(ex-sx+1),1:(ey-sy+1),1:(ez_uv-sz_uv+1)) !<
    double precision, intent(in) :: w2(1:(ex-sx+1),1:(ey-sy+1),1:(ez_w-sz_w+1)) !<
    double precision, intent(in) :: conv0_w(1:(ex-sx+1),1:(ey-sy+1),1:(ez_w-sz_w+1)) !<
    double precision, intent(in) :: conv1_w(1:(ex-sx+1),1:(ey-sy+1),1:(ez_w-sz_w+1)) !<
    double precision, intent(in) :: pres(1:(ex-sx+1),1:(ey-sy+1),1:(ez_p-sz_p+1)) !<
    double precision, intent(in) :: pres_old(1:(ex-sx+1),1:(ey-sy+1),1:(ez_p-sz_p+1)) !<
    double precision, intent(in) :: phi2(1:(ex-sx+1),1:(ey-sy+1),1:(ez_p-sz_p+1)) !<
    double precision, intent(in) :: conv0_phi(1:(ex-sx+1),1:(ey-sy+1),1:(ez_p-sz_p+1)) !<
    double precision, intent(in) :: conv1_phi(1:(ex-sx+1),1:(ey-sy+1),1:(ez_p-sz_p+1)) !<
    integer,          intent(in) :: rank !< process ID
    integer,          intent(in) :: nprocs_x !< number of procs in x dim
    integer,          intent(in) :: nprocs_y !< number of procs in y dim
    integer,          intent(in) :: nprocs_z !< number of procs in z dim
    !Netcdf variables
    integer cmode, ncid, varid, dimid(3)
    integer psizes(3), gsizes(3), start(3), count(3)
    character(len=60) :: dataname

    dataname = "u2_backup"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_uv,ez_uv,&
       u2,nprocs_x,nprocs_y,nprocs_z,rank,backup_no)

    dataname = "conv0_u_backup"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_uv,ez_uv,&
       conv0_u,nprocs_x,nprocs_y,nprocs_z,rank,backup_no)

    dataname = "conv1_u_backup"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_uv,ez_uv,&
       conv1_u,nprocs_x,nprocs_y,nprocs_z,rank,backup_no)

    dataname = "v2_backup"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_uv,ez_uv,&
       v2,nprocs_x,nprocs_y,nprocs_z,rank,backup_no)

    dataname = "conv0_v_backup"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_uv,ez_uv,&
       conv0_v,nprocs_x,nprocs_y,nprocs_z,rank,backup_no)

    dataname = "conv1_v_backup"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_uv,ez_uv,&
       conv1_v,nprocs_x,nprocs_y,nprocs_z,rank,backup_no)

    dataname = "w2_backup"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_w,ez_w,&
       w2,nprocs_x,nprocs_y,nprocs_z,rank,backup_no)

    dataname = "conv0_w_backup"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_w,ez_w,&
       conv0_w,nprocs_x,nprocs_y,nprocs_z,rank,backup_no)

    dataname = "conv1_w_backup"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_w,ez_w,&
       conv1_w,nprocs_x,nprocs_y,nprocs_z,rank,backup_no)

    dataname = "phi2_backup"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_p,ez_p,&
       phi2,nprocs_x,nprocs_y,nprocs_z,rank,backup_no)

    dataname = "conv0_phi_backup"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_p,ez_p,&
       conv0_phi,nprocs_x,nprocs_y,nprocs_z,rank,backup_no)

    dataname = "conv1_phi_backup"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_p,ez_p,&
       conv1_phi,nprocs_x,nprocs_y,nprocs_z,rank,backup_no)

    dataname = "pres_backup"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_p,ez_p,&
       pres,nprocs_x,nprocs_y,nprocs_z,rank,backup_no)

    dataname = "pres_old_backup"
    call output_3D_hdf5(dataname,sx,ex,sy,ey,sz_p,ez_p,&
       pres_old,nprocs_x,nprocs_y,nprocs_z,rank,backup_no)

    if(rank==0)then
       !Backup dx,dy,dz information as well
       !TO DO
    end if

    return
  end subroutine output_backup_hdf5


  !> Load level-set function data from a file.
  !! This is the complement of channel_output_phi and expects the
  !! input data to be formatted in the way output by that
  !! sub-routine.
  subroutine channel_input_phi(filename,phi,l,m,n,ierr)
    implicit none
    character(len=*), intent(in)    :: filename !< File name.
    integer,          intent(in)    :: l !< L dimension bound.
    integer,          intent(in)    :: m !< M dimension bound.
    integer,          intent(in)    :: n !< N dimension bound.
    double precision, intent(inout) :: phi(0:l,0:m,0:n) !< Level-set function.
    integer,          intent(out)   :: ierr !< 0 if no errrors.

    integer          :: rows    ! Row count.
    integer          :: i, j, k ! Loop variables.
    double precision :: x, y, z ! Temporary variables for coordinate vector.

    open(unit=9,file=filename,status='unknown')
    read(9,*) ! Skip 1st meta-data line.
    read(9,*) ! Skip 2nd meta-data line.
    rows = 0
    do k=0,n-2
       do j=0,m-2
          do i=1,l-1
             read(9,*) x,y,z,phi(i+1,j,k+1)
             rows = rows + 1
          end do
       end do
    end do
    close(9)
    ierr = 0
    if (rows /= ((l-1)*(m-1)*(n-1))) then
       ierr = 1
    end if
    return
  end subroutine channel_input_phi


  !> Save velocity and pressure data to a file.
  !! The first two lines of the file are meta-data of form:
  !! \code
  !! VARIABLES="X ","Y ","Z ","U","V","W","P"
  !! ZONE T="Floor" I=l-1 J=m-1 K=n-1 F=POINT
  !! \endcode
  !! This is followed by (l-1)*(m-1)*(n-1) rows of form:
  !! \code
  !! X Y Z U V W p(L,M,N+1)
  !! \endcode
  !! where:
  !! * L:1..l-1
  !! * M:0..m-2
  !! * N:0..n-2 
  subroutine channel_output_uvw(u3,v3,w3,p,l,m,n,dx,dy,dz,filename)
    implicit none
    integer,          intent(in) :: l !< L dimension bound.
    integer,          intent(in) :: m !< M dimension bound.
    integer,          intent(in) :: n !< N dimension bound.
    double precision, intent(in) :: u3(0:l,0:m,0:n-2) !< Velocity.
    double precision, intent(in) :: v3(0:l,0:m,0:n-2) !< Velocity.
    double precision, intent(in) :: w3(0:l,0:m,0:n-1) !< Velocity.
    double precision, intent(in) :: p(0:l,0:m,0:n)    !< Pressure.
    double precision, intent(in) :: dx !< X dimension element size.
    double precision, intent(in) :: dy !< Y dimension element size.
    double precision, intent(in) :: dz !< Z dimension element size.
    character(len=*), intent(in) :: filename !< File name.

    integer          :: i, j, k ! Loop variables.
    double precision :: x, y, z ! Temporary variables for coordinate vector.
    double precision :: u, v, w ! Temprorary variables for velocity vector.

    open(unit=9,file=filename,status='unknown')
    write(9,*) 'VARIABLES="X ","Y ","Z ","U","V","W","P"'
    write(9,*) 'ZONE T="Floor", I=',l-1,' J=',m-1,' K=',n-1,' F=POINT'
    do k=0,n-2
       z=dz/2.d0+dz*k              
       do j=0,m-2
          y=dy/2.d0+dy*j
          do i=1,l-1
             x=dx/2.d0+dx*i
             u=(u3(i,j,k)+u3(i+1,j,k))/2.d0
             v=(v3(i,j,k)+v3(i,j+1,k))/2.d0
             w=(w3(i,j,k)+w3(i,j,k+1))/2.d0
             write(9,*) x,y,z,u,v,w,p(i,j,k+1)
          end do
       end do
    end do
    close(9)
    return
  end subroutine channel_output_uvw


  !> Save velocity, pressure, level-set and viscosity data to a file.
  !! The first two lines of the file are meta-data of form:
  !! \code
  !! VARIABLES="X ","Z ","U","V","W","P","fi","Viscosity"
  !! ZONE T="Floor" I=l-1 K=n-1 F=POINT
  !! \endcode
  !! This is followed by (l-1)*(n-1) rows of form:
  !! \code
  !! X Y Z U V W p(L,M,N+1),phi(L+1,M,N+1),visc(L,M,N+1)
  !! \endcode
  !! where:
  !! * L:1..l-1
  !! * M=(m-2)/2
  !! * N:0..n-2 
  subroutine channel_output2d(u3,v3,w3,p,phi,visc,l,m,n,dx,dy,dz,filename)
    implicit none
    integer,          intent(in) :: l !< L dimension bound.
    integer,          intent(in) :: m !< M dimension bound.
    integer,          intent(in) :: n !< N dimension bound.
    double precision, intent(in) :: u3(0:l,0:m,0:n-2) !< Velocity.
    double precision, intent(in) :: v3(0:l,0:m,0:n-2) !< Velocity.
    double precision, intent(in) :: w3(0:l,0:m,0:n-1) !< Velocity.
    double precision, intent(in) :: p(0:l,0:m,0:n)    !< Pressure.
    double precision, intent(in) :: phi(0:l,0:m,0:n)  !< Level-set function.
    double precision, intent(in) :: visc(0:l,0:m,0:n) !< Viscosity
    double precision, intent(in) :: dx !< X dimension element size.
    double precision, intent(in) :: dy !< Y dimension element size.
    double precision, intent(in) :: dz !< Z dimension element size.
    character(len=*), intent(in) :: filename !< File name.

    integer          :: i, j, k ! Loop variables.
    double precision :: x, y, z ! Temporary variables for coordinate vector.
    double precision :: u, v, w ! Temprorary variables for velocity vector.

    open(unit=9,file=filename,status='unknown')
    write(9,*)'VARIABLES="X ","Z ","U","V","W","P","fi","Viscosity"'
    write(9,*)'ZONE T="Floor", I=',l-1,' K=',n-1,' F=POINT'

    j=(m-1)/2

    do k=0,n-2
       z=dz/2.d0+dz*k              
       y=dy/2.d0+dy*j
       do i=1,l-1
          x=dx/2.d0+dx*i
          u=(u3(i,j,k)+u3(i+1,j,k))/2.d0
          v=(v3(i,j,k)+v3(i,j+1,k))/2.d0
          w=(w3(i,j,k)+w3(i,j,k+1))/2.d0
          write(9,*) x,z,u,v,w,p(i,j,k+1),phi(i+1,j,k+1),visc(i,j,k+1)
       end do
    end do
    close(9)
    return
  end subroutine channel_output2d


  !> Save a back-up file.
  subroutine backup_channel(u2,v2,w2, &
       conv0_u,conv1_u,conv0_v,conv1_v,conv0_w,conv1_w, &
       pres_old,pres,phi,conv0_phi,conv1_phi,l,m,n,filename)
    implicit none
    integer,          intent(in) :: l !< L dimension bound.
    integer,          intent(in) :: m !< M dimension bound.
    integer,          intent(in) :: n !< N dimension bound.
    double precision, intent(in) :: u2(0:l,0:m,0:n-2) !< Velocity.
    double precision, intent(in) :: v2(0:l,0:m,0:n-2) !< Velocity.
    double precision, intent(in) :: w2(0:l,0:m,0:n-1) !< Velocity.
    double precision, intent(in) :: conv0_u(0:l,0:m,0:n-2)
    double precision, intent(in) :: conv1_u(0:l,0:m,0:n-2)
    double precision, intent(in) :: conv0_v(0:l,0:m,0:n-2)
    double precision, intent(in) :: conv1_v(0:l,0:m,0:n-2)
    double precision, intent(in) :: conv0_w(0:l,0:m,0:n-1)
    double precision, intent(in) :: conv1_w(0:l,0:m,0:n-1)
    double precision, intent(in) :: pres_old(0:l,0:m,0:n) !< Old pressure.
    double precision, intent(in) :: pres(0:l,0:m,0:n) !< Pressure.
    double precision, intent(in) :: phi(0:l,0:m,0:n) !< PHI.
    double precision, intent(in) :: conv0_phi(0:l,0:m,0:n)
    double precision, intent(in) :: conv1_phi(0:l,0:m,0:n)
    character(len=*), intent(in) :: filename !< File name.

    integer          :: i, j, k ! Loop variables.

    open(unit=9,file=filename,status='unknown')

    do k=0,n-2
       do j=0,m
          do i=0,l
             write(9,*) u2(i,j,k),conv0_u(i,j,k),conv1_u(i,j,k)
          end do
       end do
    end do

    do k=0,n-2
       do j=0,m
          do i=0,l
             write(9,*) v2(i,j,k),conv0_v(i,j,k),conv1_v(i,j,k)
          end do
       end do
    end do

    do k=0,n-1
       do j=0,m
          do i=0,l
             write(9,*) w2(i,j,k),conv0_w(i,j,k),conv1_w(i,j,k)
          end do
       end do
    end do

    do k=0,n
       do j=0,m
          do i=0,l
             write(9,*) phi(i,j,k),conv0_phi(i,j,k),conv1_phi(i,j,k)
          end do
       end do
    end do

    do k=0,n
       do j=0,m
          do i=0,l
             write(9,*) pres(i,j,k),pres_old(i,j,k)
          end do
       end do
    end do

    close(9)
    return
  end subroutine backup_channel


  !> Load a back-up file.
  !! This is the complement of backup_channel.
  subroutine dataload_channel(u2,v2,w2, &
       conv0_u,conv1_u,conv0_v,conv1_v,conv0_w,conv1_w, &
       pres_old,pres,phi,conv0_phi,conv1_phi,l,m,n,filename)
    implicit none

    integer,          intent(in) :: l !< L dimension bound.
    integer,          intent(in) :: m !< M dimension bound.
    integer,          intent(in) :: n !< N dimension bound.
    double precision, intent(inout) :: u2(0:l,0:m,0:n-2) !< Velocity.
    double precision, intent(inout) :: v2(0:l,0:m,0:n-2) !< Velocity.
    double precision, intent(inout) :: w2(0:l,0:m,0:n-1) !< Velocity.
    double precision, intent(inout) :: conv0_u(0:l,0:m,0:n-2)
    double precision, intent(inout) :: conv1_u(0:l,0:m,0:n-2)
    double precision, intent(inout) :: conv0_v(0:l,0:m,0:n-2)
    double precision, intent(inout) :: conv1_v(0:l,0:m,0:n-2)
    double precision, intent(inout) :: conv0_w(0:l,0:m,0:n-1)
    double precision, intent(inout) :: conv1_w(0:l,0:m,0:n-1)
    double precision, intent(inout) :: pres_old(0:l,0:m,0:n) !< Old pressure.
    double precision, intent(inout) :: pres(0:l,0:m,0:n) !< Pressure.
    double precision, intent(inout) :: phi(0:l,0:m,0:n) !< PHI.
    double precision, intent(inout) :: conv0_phi(0:l,0:m,0:n)
    double precision, intent(inout) :: conv1_phi(0:l,0:m,0:n)
    character(len=*), intent(in)    :: filename !< File name.

    integer          :: i, j, k ! Loop variables.

    open(9,file=filename,status='old')

    do k=0,n-2
       do j=0,m
          do i=0,l
             read(9,*) u2(i,j,k),conv0_u(i,j,k),conv1_u(i,j,k)
          end do
       end do
    end do

    do k=0,n-2
       do j=0,m
          do i=0,l
             read(9,*) v2(i,j,k),conv0_v(i,j,k),conv1_v(i,j,k)
          end do
       end do
    end do

    do k=0,n-1
       do j=0,m
          do i=0,l
             read(9,*) w2(i,j,k),conv0_w(i,j,k),conv1_w(i,j,k)
          end do
       end do
    end do

    do k=0,n
       do j=0,m
          do i=0,l
             read(9,*) phi(i,j,k),conv0_phi(i,j,k),conv1_phi(i,j,k)
          end do
       end do
    end do

    do k=0,n
       do j=0,m
          do i=0,l
             read(9,*) pres(i,j,k),pres_old(i,j,k)
          end do
       end do
    end do

    close(9)
    return
  end subroutine dataload_channel


  !> Map from a 1D row number to a 3D grid index.
  !! Given grid dimensions (lx:ly,mx:my,nx:ny) that has been converted
  !! into a 1D array of form (1:(ly-lx+1)*(my-mx+1)*(ny-nx+1)) 
  !! calculate the (l,m,n) index corresponding to the given row:
  !! * rl = ((row - 1) % (ly - lx + 1)) + lx
  !! * rm = (((row - 1) / (ly - lx + 1)) % (my - mx + 1)) + mx
  !! * rn = (((row - 1) / ((ly - lx + 1) * (my - mx + 1)))) + nx
  subroutine row_to_grid(lx,ly,mx,my,nx,ny,row,l,m,n)
    integer,intent(in)  :: lx  !< Minimum L dimension bound.
    integer,intent(in)  :: ly  !< Maximum L dimension bound.
    integer,intent(in)  :: mx  !< Minimum M dimension bound.
    integer,intent(in)  :: my  !< Maximum M dimension bound.
    integer,intent(in)  :: nx  !< Minimum N dimension bound.
    integer,intent(in)  :: ny  !< Maximum N dimension bound.
    integer,intent(in)  :: row !< Row.
    integer,intent(out) :: l   !< L value.
    integer,intent(out) :: m   !< M value.
    integer,intent(out) :: n   !< N value.

    l = mod(row-1,(ly-lx+1)) + lx
    m = mod((row-1)/(ly-lx+1),(my-mx+1)) + mx
    n = (row-1)/((ly-lx+1) * (my-mx+1)) + nx
  end subroutine row_to_grid

end module tpls_io

subroutine check(err, message)
  use netcdf 
  use mpi
  implicit none
  
  integer err
  character(len=*) message
  
  ! It is a good idea to check returned value for possible error
  if (err .NE. NF90_NOERR) then
     write(*,*) "Aborting: ",trim(message), trim(nf90_strerror(err))
     call MPI_Abort(MPI_COMM_WORLD, -1, err)
  end if
end subroutine check
